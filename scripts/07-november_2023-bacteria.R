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
  packages = c(
    "here",
    "readr",
    "dplyr",
    "tidyr",
    "stringr",
    "phyloseq",
    "ape",
    "microViz",
    "vegan",
    "picante",
    "ape",
    "ampvis2",
    "forcats",
    "FSA",
    "pairwiseAdonis",
    "purrr"
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
  # extract the sample names from the metadata in the same order
  # of custom labels
  tar_target(
    name = custom_order_sample_names,
    command = nov_2023_metadata %>%
      dplyr::distinct(sample_name) %>%
      dplyr::pull(sample_name) %>%
      as.character()
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
  # Add tree to phyloseq object
  tar_target(
    name = tree_path,
    command = here::here(
      "results-bacteria",
      "qiime2",
      "phylogenetic_tree",
      "tree.nwk"
    ),
    format = "file"
  ),
  tar_target(
    name = phyloseq_obj_with_tree,
    command = add_tree_to_phyloseq(phyloseq_obj, tree_path)
  ),
  # remove chloroplast from data
  tar_target(
    name = phyloseq_obj_pruned,
    command = prune_phyloseq(
      phyloseq_obj_with_tree,
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
  # determine rarefaction depth
  tar_target(
    nov_2023_rarefaction_depth,
    min(sample_sums(nov_2023_phyloseq_obj))
  ),
  # Generate a sequence of iterations
  tar_target(
    name = thousand_iterations,
    command = seq_len(1000)
  ),
  # Perform rarefaction across multiple iterations
  tar_target(
    name = nov_2023_phyloseq_obj_rarefied,
    command = rarefy_phyloseq_object(
      nov_2023_phyloseq_obj,
      nov_2023_rarefaction_depth,
    ),
    pattern = map(thousand_iterations),
    iteration = "list"
  ),
  # collect sample data to plot sequencing depth
  tar_target(
    name = nov_2023_samples_data,
    command = sort_by_factor_column(
      get_samples_data(nov_2023_phyloseq_obj),
      "label",
      custom_order_labels
    )
  ),
  # calculate bray distances to perform a permanova on technical replicates
  tar_target(
    name = nov_2023_bray_distance_matrices,
    command = calculate_distance_matrix(nov_2023_phyloseq_obj_rarefied, method = "bray"),
    pattern = map(nov_2023_phyloseq_obj_rarefied),
    iteration = "list"
  ),
  tar_target(
    name = nov_2023_bray_distance_matrix,
    command = as.dist(
      Reduce("+", nov_2023_bray_distance_matrices) / length(nov_2023_bray_distance_matrices)
    )
  ),
  # calculate distances between technical replicates
  tar_target(
    name = nov_2023_tech_rep_distances,
    command = calculate_technical_replicate_distances(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      group_column = "sample_group"
    )
  ),
  # a permanova to test if technical replicates are different
  tar_target(
    name = nov_2023_tech_rep_permanova,
    command = calculate_permanova(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      columns = c("sample_group", "technical_rep")
    )
  ),
  # merge technical replicates
  # first update metadata for merged samples
  tar_target(
    name = nov_2023_metadata_merged,
    command = {
      tmp <- merge_metadata(
        nov_2023_metadata,
        column = "sample_name"
      )
      sort_by_factor_column(
        tmp,
        "label",
        custom_order_labels
      )
    }
  ),
  # then merge phyloseq object
  tar_target(
    name = nov_2023_phyloseq_obj_merged,
    command = merge_technical_replicates(
      nov_2023_phyloseq_obj,
      nov_2023_metadata_merged,
      column = "sample_name"
    )
  ),
  # determine rarefaction depth for merged data
  tar_target(
    nov_2023_merged_rarefaction_depth,
    min(sample_sums(nov_2023_phyloseq_obj_merged))
  ),
  # Perform rarefaction on merged data
  tar_target(
    name = nov_2023_phyloseq_obj_merged_rarefied,
    command = rarefy_phyloseq_object(
      nov_2023_phyloseq_obj_merged,
      nov_2023_merged_rarefaction_depth,
    ),
    pattern = map(thousand_iterations),
    iteration = "list"
  ),
  # collect sample data to plot sequencing depth for merged data
  tar_target(
    name = nov_2023_samples_data_merged,
    command = sort_by_factor_column(
      get_samples_data(nov_2023_phyloseq_obj_merged),
      "label",
      custom_order_labels
    )
  ),
  # here are barplots
  tar_target(
    name = nov_2023_melted_phylum,
    command = agglomerate_by_taxa(
      nov_2023_phyloseq_obj_merged,
      taxrank = "Phylum",
      sample_order = nov_2023_metadata_merged$sample_name
    )
  ),
  tar_target(
    name = nov_2023_melted_class,
    command = agglomerate_by_taxa(
      nov_2023_phyloseq_obj_merged,
      taxrank = "Class",
      sample_order = nov_2023_metadata_merged$sample_name
    )
  ),
  # and here the heatmaps
  tar_target(
    name = nov_2023_ampvis2_object,
    command = phyloseq_to_ampvis2(nov_2023_phyloseq_obj_merged)
  ),
  tar_target(
    name = nov_2023_heatmap_phylum,
    command = custom_heatmap(
      nov_2023_ampvis2_object,
      group_by = "label",
      order_x_by = custom_order_labels,
      showRemainingTaxa = TRUE
    )
  ),
  tar_target(
    name = nov_2023_heatmap_class,
    command = custom_heatmap(
      nov_2023_ampvis2_object,
      group_by = "label",
      order_x_by = custom_order_labels,
      showRemainingTaxa = TRUE,
      tax_add = "Class",
      tax_show = 20
    )
  ),
  # collect otu table
  tar_target(
    name = nov_2023_otu_table_matrix,
    command = get_otu_table(nov_2023_phyloseq_obj)
  ),
  # deal with rarefaction curves
  tar_target(
    name = nov_2023_rarecurve_df,
    command = calculate_rarecurve(nov_2023_otu_table_matrix)
  ),
  # alpha diversity steps
  # calculate alpha diversity measures on each rarefied phyloseq object
  tar_target(
    name = nov_2023_samples_rarefaction,
    command = rarefy_alpha(
      nov_2023_phyloseq_obj_rarefied,
      measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "FaithPD")),
    pattern = map(nov_2023_phyloseq_obj_rarefied)
  ),
  # now transform rarefaction in a summary table
  tar_target(
    name = nov_2023_rarefaction_results,
    command = summarize_rarefactions(nov_2023_samples_rarefaction, metadata)
  ),
  # pivot data and group by sample name
  tar_target(
    name = nov_2023_rarefaction_by_sample_name,
    command = reshape_rarefaction_data(
      nov_2023_rarefaction_results,
      by_column = "sample_name"
    )
  ),
  # make plots for alpha diversity
  tar_target(
    name = nov_2023_alpha_diversity_by_sample_name,
    command = plot_alpha_diversity(
      data = nov_2023_rarefaction_by_sample_name,
      x = "sample_name",
      y = "mean",
      fill = "sample_name",
      facet = "Metric",
      title = "Alpha Diversity Metrics Across Groups",
      xlab = "Sample Name",
      ylab = "Alpha Diversity Measure",
      custom_order = custom_order_sample_names
    )
  ),
  # do the Kruskall-Wallis test
  tar_target(
    name = nov_2023_kruscal_shannon_sample_name,
    command = calculate_kruskal_wallis(
      nov_2023_rarefaction_results,
      "Shannon_mean",
      "sample_name"
    )
  ),
  # do the post-hoc tests
  tar_target(
    name = nov_2023_dunn_shannon_sample_name,
    command = calculate_dunn_test(
      nov_2023_rarefaction_results,
      "Shannon_mean",
      "sample_name"
    )
  ),
  # pivot data and group by labels
  tar_target(
    name = nov_2023_rarefaction_by_labels,
    command = reshape_rarefaction_data(
      nov_2023_rarefaction_results,
      by_column = "label"
    )
  ),
  # make plots
  tar_target(
    name = nov_2023_alpha_diversity_by_labels,
    command = plot_alpha_diversity(
      data = nov_2023_rarefaction_by_labels,
      x = "label",
      y = "mean",
      fill = "label",
      facet = "Metric",
      title = "Alpha Diversity Metrics Across Groups",
      xlab = "Sample Name",
      ylab = "Alpha Diversity Measure",
      custom_order = custom_order_labels
    )
  ),
  # do the Kruskall-Wallis test
  tar_target(
    name = nov_2023_kruscal_shannon_labels,
    command = calculate_kruskal_wallis(
      nov_2023_rarefaction_results,
      "Shannon_mean",
      "label"
    )
  ),
  # do the post-hoc tests
  tar_target(
    name = nov_2023_dunn_shannon_labels,
    command = calculate_dunn_test(
      nov_2023_rarefaction_results,
      "Shannon_mean",
      "label"
    )
  ),
  # calculate other distance metrics
  tar_target(
    name = nov_2023_jaccard_distance_matrices,
    command = calculate_distance_matrix(nov_2023_phyloseq_obj_rarefied, method = "jaccard"),
    pattern = map(nov_2023_phyloseq_obj_rarefied),
    iteration = "list"
  ),
  tar_target(
    name = nov_2023_jaccard_distance_matrix,
    command = as.dist(
      Reduce("+", nov_2023_jaccard_distance_matrices) / length(nov_2023_jaccard_distance_matrices)
    )
  ),
  tar_target(
    name = nov_2023_wunifrac_distance_matrices,
    command = phyloseq::distance(
      nov_2023_phyloseq_obj_rarefied,
      method = "wunifrac"
    ),
    pattern = map(nov_2023_phyloseq_obj_rarefied),
    iteration = "list"
  ),
  tar_target(
    name = nov_2023_wunifrac_distance_matrix,
    command = as.dist(
      Reduce("+", nov_2023_wunifrac_distance_matrices) / length(nov_2023_wunifrac_distance_matrices)
    )
  ),
  tar_target(
    name = nov_2023_unifrac_distance_matrices,
    command = phyloseq::distance(
      nov_2023_phyloseq_obj_rarefied,
      method = "unifrac"
    ),
    pattern = map(nov_2023_phyloseq_obj_rarefied),
    iteration = "list"
  ),
  tar_target(
    name = nov_2023_unifrac_distance_matrix,
    command = as.dist(
      Reduce("+", nov_2023_unifrac_distance_matrices) / length(nov_2023_unifrac_distance_matrices)
    )
  ),
  # ordinations (beta diversity)
  # PCoA
  tar_target(
    name = nov_2023_bray_pcoa_object,
    command = calculate_pcoa(nov_2023_bray_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_jaccard_pcoa_object,
    command = calculate_pcoa(nov_2023_jaccard_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_wunifrac_pcoa_object,
    command = calculate_pcoa(nov_2023_wunifrac_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_unifrac_pcoa_object,
    command = calculate_pcoa(nov_2023_unifrac_distance_matrix, nov_2023_metadata)
  ),
  # NMDS
  tar_target(
    name = nov_2023_bray_nmds_object,
    command = calculate_nmds(nov_2023_bray_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_jaccard_nmds_object,
    command = calculate_nmds(nov_2023_jaccard_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_wunifrac_nmds_object,
    command = calculate_nmds(nov_2023_wunifrac_distance_matrix, nov_2023_metadata)
  ),
  tar_target(
    name = nov_2023_unifrac_nmds_object,
    command = calculate_nmds(nov_2023_unifrac_distance_matrix, nov_2023_metadata)
  ),
  # calculate distances with and between groups
  # bray
  tar_target(
    name = nov_2023_bray_distance_by_sample_name,
    command = get_distances(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      column = "sample_name"
    )
  ),
  tar_target(
    name = nov_2023_bray_distance_by_sample_name_plot,
    command = plot_distances(nov_2023_bray_distance_by_sample_name, column = "sample_name")
  ),
  tar_target(
    name = nov_2023_bray_distance_by_labels,
    command = get_distances(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      column = "label"
    )
  ),
  tar_target(
    name = nov_2023_bray_distance_by_labels_plot,
    command = plot_distances(nov_2023_bray_distance_by_labels, column = "label")
  ),
  # unifrac
  tar_target(
    name = nov_2023_unifrac_distance_by_sample_name,
    command = get_distances(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      column = "sample_name"
    )
  ),
  tar_target(
    name = nov_2023_unifrac_distance_by_sample_name_plot,
    command = plot_distances(nov_2023_unifrac_distance_by_sample_name, column = "sample_name")
  ),
  tar_target(
    name = nov_2023_unifrac_distance_by_labels,
    command = get_distances(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      column = "label"
    )
  ),
  tar_target(
    name = nov_2023_unifrac_distance_by_labels_plot,
    command = plot_distances(nov_2023_unifrac_distance_by_labels, column = "label")
  ),
  # calculate beta dispersion
  tar_target(
    name = nov_2023_bray_sample_name_beta_dispersion,
    command = calculate_beta_dispersion(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      column = "sample_name",
      bias.adjust = TRUE,
      levels = custom_order_sample_names
    )
  ),
  tar_target(
    name = nov_2023_bray_labels_beta_dispersion,
    command = calculate_beta_dispersion(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      column = "label",
      bias.adjust = TRUE,
      levels = custom_order_labels
    )
  ),
  tar_target(
    name = nov_2023_unifrac_sample_name_beta_dispersion,
    command = calculate_beta_dispersion(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      column = "sample_name",
      bias.adjust = TRUE,
      levels = custom_order_sample_names
    )
  ),
  tar_target(
    name = nov_2023_unifrac_labels_beta_dispersion,
    command = calculate_beta_dispersion(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      column = "label",
      bias.adjust = TRUE,
      levels = custom_order_labels
    )
  ),
  # permanova on distance matrix
  ## Bray
  tar_target(
    name = nov_2023_bray_sample_name_permanova,
    command = calculate_permanova(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = nov_2023_pairwise_bray_sample_name_permanova,
    command = calculate_pairwise_permanova(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = nov_2023_bray_labels_permanova,
    command = calculate_permanova(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      columns = c("label")
    )
  ),
  tar_target(
    name = nov_2023_pairwise_bray_labels_permanova,
    command = calculate_pairwise_permanova(
      nov_2023_bray_distance_matrix,
      nov_2023_metadata,
      columns = c("label")
    )
  ),
  ## unifrac
  tar_target(
    name = nov_2023_unifrac_sample_name_permanova,
    command = calculate_permanova(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = nov_2023_pairwise_unifrac_sample_name_permanova,
    command = calculate_pairwise_permanova(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = nov_2023_unifrac_labels_permanova,
    command = calculate_permanova(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      columns = c("label")
    )
  ),
  tar_target(
    name = nov_2023_pairwise_unifrac_labels_permanova,
    command = calculate_pairwise_permanova(
      nov_2023_unifrac_distance_matrix,
      nov_2023_metadata,
      columns = c("label")
    )
  ),
  # render november 2023 quarto document
  tar_quarto(
    name = november_2023,
    path = "analysis/07-november_2023-bacteria.qmd",
    quiet = TRUE
  )
)
