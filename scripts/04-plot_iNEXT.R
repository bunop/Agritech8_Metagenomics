# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(parallel)
library(crew)
library(quarto)

# Set target options:
tar_option_set(
  packages = c(
    "here",
    "readr",
    "dplyr",
    "stringr",
    "phyloseq",
    "iNEXT"
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
    workers = detectCores() - 1,
    seconds_idle = 60,
    tasks_max = 50
  ),
  #
  # Alternatively, if you want workers to run on a high-performance computing
  # cluster, select a controller from the {crew.cluster} package.
  # For the cloud, see plugin packages like {crew.aws.batch}.
  # The following example is a controller for Sun Grid Engine (SGE).
  #
  #   controller = crew.cluster::crew_controller_sge(
  #     # Number of workers that the pipeline can scale up to:
  #     workers = 10,
  #     # It is recommended to set an idle time so workers can shut themselves
  #     # down if they are not running tasks.
  #     seconds_idle = 120,
  #     # Many clusters install R as an environment module, and you can load it
  #     # with the script_lines argument. To select a specific verison of R,
  #     # you may need to include a version string, e.g. "module load R/4.3.2".
  #     # Check with your system administrator if you are unsure.
  #     script_lines = "module load R"
  #   )
  #
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
  tar_target(
    name = phyloseq_path,
    command = here::here("results-bacteria", "phyloseq", "dada2_phyloseq.rds"),
    format = "file"
  ),
  tar_target(
    name = phyloseq_object,
    command = load_phyloseq(metadata, phyloseq_path)
  ),
  tar_target(
    name = otu_table_matrix,
    command = get_otu_table(phyloseq_object)
  ),
  tar_target(
    name = coverage_stats,
    command = calculate_coverage_stats(otu_table_matrix)
  ),
  tar_target(
    name = abundance_list,
    command = split(as.data.frame(otu_table_matrix), rownames(otu_table_matrix))
  ),
  tar_target(
    name = parallel_iNEXT_results,
    command = list(
      iNEXT(as.numeric(abundance_list[[1]]), q = 0, datatype = "abundance", nboot = 100)
    ),
    pattern = map(abundance_list)
  ),
  tar_target(
    name = combined_iNEXT_results,
    command = combine_iNEXT_results(parallel_iNEXT_results, abundance_list)
  ),
  tar_target(
    name = rarecurve_df,
    command = calculate_rarecurve(otu_table_matrix)
  ),
  tar_quarto(
    plot_iNEXT,
    "analysis/04-plot_iNEXT.qmd",
    quiet = TRUE
  )
)
