# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(parallelly)
library(crew)
library(quarto)

# those packages are not strictly related to targets or quarto, I need them to
# avoid to create a custom function to do a simple filter with dplyr, for example
library(dplyr)
library(phyloseq)

# Set target options:
tar_option_set(
  # Set packages that your targets need to run. This is required for non-base packages.
  packages = c(
    "here",
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "forcats",
    "phyloseq",
    "microViz",
    "iNEXT",
    "purrr",
    "tibble"
  ), # Packages that your targets need for their tasks.
  # format = "qs", # Optionally set the default storage format. qs is fast.
  #
  # Pipelines that take a long time to run may benefit from
  # optional distributed computing. To use this capability
  # in tar_make(), supply a {crew} controller
  # as discussed at https://books.ropensci.org/targets/crew.html.
  # Choose a controller that suits your needs. For example, the following
  # sets a controller that scales up to a maximum of two workers
  # which run as local R processes. Each worker launches when there is work
  # to do and exits if 60 seconds pass with no tasks to run.
  #
  controller = crew::crew_controller_local(
    workers = parallelly::availableCores() - 1,
    seconds_idle = 60,
    tasks_max = 50
  ),
  # Set other options as needed.
  seed = 42
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# tar_source("other_functions.R") # Source other scripts as needed.

# Replace the target list below with your own:
list(
  tar_target(
    name = metadata_path,
    command = here::here("data", "metadata_bacteria_fix.tsv"),
    format = "file"
  ),
  tar_target(
    name = metadata,
    command = load_metadata(metadata_path)
  ),
  # select with dplyr where date is "nov 23"
  tar_target(
    name = nov_2023_metadata_tmp,
    command = metadata %>%
      dplyr::filter(date == "nov 23")
  ),
  # attempt to order metadata and phyloseq object using label column
  # and a custom levels order
  tar_target(
    name = custom_order_labels,
    command = c("T0", "NA", "NP", "CS", "LM")
  ),
  tar_target(
    name = nov_2023_metadata,
    command = sort_by_factor_column(
      nov_2023_metadata_tmp,
      "label",
      custom_order_labels
    )
  ),
  # open the phyloseq object
  tar_target(
    name = phyloseq_path,
    command = here::here("results-bacteria", "phyloseq", "dada2_phyloseq.rds"),
    format = "file"
  ),
  tar_target(
    name = phyloseq_obj,
    command = load_phyloseq(metadata, phyloseq_path)
  ),
  # remove chloroplast from data, same as november_2023-bacteria, so that
  # richness/Chao1 estimates are not inflated by non-bacterial ASVs
  tar_target(
    name = phyloseq_obj_pruned,
    command = prune_phyloseq(
      phyloseq_obj,
      rank = "Phylum",
      items = ("Cyanobacteria")
    )
  ),
  # subsetting samples from phyloseq object
  tar_target(
    name = nov_2023_phyloseq_obj_tmp,
    command = phyloseq::subset_samples(
      phyloseq_obj_pruned,
      date == "nov 23"
    )
  ),
  # now order the samples in the phyloseq object using the sorted_metadata
  tar_target(
    name = nov_2023_phyloseq_obj,
    command = sort_phyloseq(
      nov_2023_phyloseq_obj_tmp,
      nov_2023_metadata$sampleID
    )
  ),
  # the reference depth for the completeness check: same depth used for
  # rarefaction/diversity stats in the november_2023-bacteria project
  tar_target(
    name = nov_2023_rarefaction_depth,
    command = min(phyloseq::sample_sums(nov_2023_phyloseq_obj))
  ),
  # collect otu table
  tar_target(
    name = nov_2023_otu_table_matrix,
    command = get_otu_table(nov_2023_phyloseq_obj)
  ),
  # one abundance vector per sample, as expected by iNEXT
  tar_target(
    name = nov_2023_abundance_list,
    command = split(
      as.data.frame(nov_2023_otu_table_matrix),
      rownames(nov_2023_otu_table_matrix)
    ) %>%
      purrr::map(~as.numeric(.))
  ),
  # extrapolate far enough that the reference depth is always covered, even
  # for the shallowest sample
  tar_target(
    name = nov_2023_iNEXT_endpoint,
    command = 2 * max(phyloseq::sample_sums(nov_2023_phyloseq_obj))
  ),
  # run iNEXT (q = 0, species richness) for each sample
  tar_target(
    name = nov_2023_iNEXT_results,
    command = iNEXT(
      nov_2023_abundance_list,  # ← Passa direttamente la lista
      q = 0,
      datatype = "abundance",
      endpoint = nov_2023_iNEXT_endpoint,
      nboot = 100
    ),
    pattern = map(nov_2023_abundance_list),
    iteration = "list"
  ),
  # compare observed/extrapolated richness at nov_2023_rarefaction_depth with
  # the asymptotic Chao1 estimate, per sample
  tar_target(
    name = nov_2023_sampling_completeness,
    command = calculate_sampling_completeness(
      nov_2023_iNEXT_results,
      names(nov_2023_abundance_list),
      nov_2023_rarefaction_depth
    ) %>%
      dplyr::left_join(
        nov_2023_metadata %>%
          dplyr::select(sampleID, label),
        by = c("Assemblage" = "sampleID")
      ) %>%
      dplyr::rename(Group = label) %>%
      dplyr::relocate(Group, .after = Assemblage) %>%
      dplyr::arrange(Coverage_ratio)
  ),
  # group-level summary, used for the report's summary paragraph
  tar_target(
    name = nov_2023_completeness_group_summary,
    command = nov_2023_sampling_completeness %>%
      dplyr::group_by(Group) %>%
      dplyr::summarise(
        mean_coverage = mean(Coverage_ratio, na.rm = TRUE),
        min_coverage = min(Coverage_ratio, na.rm = TRUE),
        n_below_90 = sum(Flag_below_90pct, na.rm = TRUE),
        n_samples = dplyr::n(),
        .groups = "drop"
      )
  ),
  # render sampling completeness quarto document
  tar_quarto(
    name = november_2023_sampling_completeness,
    path = "analysis/09-november_2023-sampling_completeness.qmd",
    quiet = TRUE
  )
)
