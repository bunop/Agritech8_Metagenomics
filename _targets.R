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
    "tidyr",
    "stringr",
    "phyloseq",
    "vegan",
    "ape",
    "ampvis2",
    "FSA",
    "pairwiseAdonis"
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
  controller = crew::crew_controller_local(workers = detectCores() - 1, seconds_idle = 60),
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
    name = metadata,
    command = load_metadata("data", "metadata_bacteria_fix.tsv")
  ),
  tar_target(
    name = phyloseq_object,
    command = load_phyloseq(metadata, "results-bacteria/phyloseq/dada2_phyloseq.rds")
  ),
  tar_target(
    name = otu_table_matrix,
    command = get_otu_table(phyloseq_object)
  ),
  tar_target(
    name = samples_data,
    command = get_samples_data(phyloseq_object)
  ),
  tar_target(
    name = coverage_stats,
    command = calculate_coverage_stats(otu_table_matrix)
  ),
  tar_target(
    name = bray_distance_matrix,
    command = calculate_distance_matrix(
      otu_table_matrix,
      min_sequencing_depth = rarefaction_depth,
      dmethod = "bray"
    )
  ),
  tar_target(
    name = tech_rep_permanova,
    command = calculate_permanova(
      bray_distance_matrix,
      metadata,
      columns = c("sample_group", "technical_rep")
    )
  ),
  # here are barplots
  tar_target(
    name = melted_phylum,
    command = agglomerate_by_taxa(
      phyloseq_object,
      taxrank = "Phylum",
      sample_order = metadata$sampleID
    )
  ),
  tar_target(
    name = melted_class,
    command = agglomerate_by_taxa(
      phyloseq_object,
      taxrank = "Class",
      sample_order = metadata$sampleID
    )
  ),
  tar_target(
    name = ampvis2_object,
    command = phyloseq_to_ampvis2(phyloseq_object)
  ),
  tar_target(
    name = heatmap_phylum,
    command = custom_heatmap(
      ampvis2_object,
      group_by = "date_condition",
      showRemainingTaxa = TRUE
    )
  ),
  tar_target(
    name = heatmap_class,
    command = custom_heatmap(
      ampvis2_object,
      group_by = "date_condition",
      showRemainingTaxa = TRUE,
      tax_add = "Class",
      tax_show = 20
    )
  ),
  # deal with rarefaction curves
  tar_target(
    rarefaction_depth,
    min(sample_sums(phyloseq_object))
  ),
  tar_target(
    name = rarecurve_df,
    command = calculate_rarecurve(otu_table_matrix)
  ),
  tar_target(
    name = observed_features_by_sample_name,
    command = load_qiime_rarefaction(
      rarefaction_csv = "results-bacteria/qiime2/alpha-rarefaction/observed_features.csv",
      metadata = metadata,
      alpha_metric = "observed_features",
      column_name = "sample_name"
    )
  ),
  tar_target(
    name = shannon_by_sample_name,
    command = load_qiime_rarefaction(
      rarefaction_csv = "results-bacteria/qiime2/alpha-rarefaction/shannon.csv",
      metadata = metadata,
      alpha_metric = "shannon",
      column_name = "sample_name"
    )
  ),
  tar_target(
    name = faith_by_sample_name,
    command = load_qiime_rarefaction(
      rarefaction_csv = "results-bacteria/qiime2/alpha-rarefaction/faith_pd.csv",
      metadata = metadata,
      alpha_metric = "faith_pd",
      column_name = "sample_name"
    )
  ),
  # alpha diversity steps
  # Generate a sequence of iterations
  tar_target(
    name = thousand_iterations,
    command = seq_len(1000)
  ),
  # Perform rarefaction across multiple iterations
  tar_target(
    name = samples_rarefaction,
    command = rarefy_alpha(
      phyloseq_object,
      rarefaction_depth,
      measures = c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher")),
    pattern = map(thousand_iterations)
  ),
  # now transform rarefaction in a summary table
  tar_target(
    name = rarefaction_results,
    command = summarize_rarefactions(samples_rarefaction, metadata)
  ),
  # pivot data and group by sample name
  tar_target(
    name = rarefaction_by_sample_name,
    command = reshape_rarefaction_data(rarefaction_results, by_column = "sample_name")
  ),
  # make plots
  tar_target(
    name = alpha_diversity_by_sample_name,
    command = plot_alpha_diversity(
      data = rarefaction_by_sample_name,
      x = "sample_name",
      y = "mean",
      fill = "sample_name",
      facet = "Metric",
      title = "Alpha Diversity Metrics Across Groups",
      xlab = "Sample Name",
      ylab = "Alpha Diversity Measure"
    )
  ),
  # do the Kruskall-Wallis test
  tar_target(
    name = kruscal_shannon_sample_name,
    command = calculate_kruskal_wallis(rarefaction_results, "Shannon_mean", "sample_name")
  ),
  # do the post-hoc tests
  tar_target(
    name = dunn_shannon_sample_name,
    command = calculate_dunn_test(rarefaction_results, "Shannon_mean", "sample_name")
  ),
  # pivot data and group by date_condition
  tar_target(
    name = rarefaction_by_date_condition,
    command = reshape_rarefaction_data(rarefaction_results, by_column = "date_condition")
  ),
  # make plots
  tar_target(
    name = alpha_diversity_by_date_condition,
    command = plot_alpha_diversity(
      data = rarefaction_by_date_condition,
      x = "date_condition",
      y = "mean",
      fill = "date_condition",
      facet = "Metric",
      title = "Alpha Diversity Metrics Across Groups",
      xlab = "Sample Name",
      ylab = "Alpha Diversity Measure"
    )
  ),
  # do the Kruskall-Wallis test
  tar_target(
    name = kruscal_shannon_date_condition,
    command = calculate_kruskal_wallis(rarefaction_results, "Shannon_mean", "date_condition")
  ),
  # do the post-hoc tests
  tar_target(
    name = dunn_shannon_date_condition,
    command = calculate_dunn_test(rarefaction_results, "Shannon_mean", "date_condition")
  ),
  # ordinations
  tar_target(
    name = pcoa_object,
    command = calculate_pcoa(bray_distance_matrix, metadata)
  ),
  tar_target(
    name = nmds_object,
    command = calculate_nmds(bray_distance_matrix, metadata)
  ),
  # calculate distances with and between groups
  tar_target(
    name = bray_distance_by_sample_name,
    command = get_distances(
      bray_distance_matrix,
      metadata,
      column = "sample_name"
    )
  ),
  tar_target(
    name = bray_distance_by_sample_name_plot,
    command = plot_distances(bray_distance_by_sample_name, column = "sample_name")
  ),
  tar_target(
    name = bray_distance_by_date_condition,
    command = get_distances(
      bray_distance_matrix,
      metadata,
      column = "date_condition"
    )
  ),
  tar_target(
    name = bray_distance_by_date_condition_plot,
    command = plot_distances(bray_distance_by_date_condition, column = "date_condition")
  ),
  # calculate beta dispersion
  tar_target(
    name = sample_name_beta_dispersion,
    command = calculate_beta_dispersion(
      bray_distance_matrix,
      metadata,
      column = "sample_name"
    )
  ),
  tar_target(
    name = date_condition_beta_dispersion,
    command = calculate_beta_dispersion(
      bray_distance_matrix,
      metadata,
      column = "date_condition"
    )
  ),
  # permanova on distance matrix
  tar_target(
    name = sample_name_permanova,
    command = calculate_permanova(
      bray_distance_matrix,
      metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = pairwise_sample_name_permanova,
    command = calculate_pairwise_permanova(
      bray_distance_matrix,
      metadata,
      columns = c("sample_name")
    )
  ),
  tar_target(
    name = date_condition_permanova,
    command = calculate_permanova(
      bray_distance_matrix,
      metadata,
      columns = c("date_condition")
    )
  ),
  tar_target(
    name = pairwise_date_condition_permanova,
    command = calculate_pairwise_permanova(
      bray_distance_matrix,
      metadata,
      columns = c("date_condition")
    )
  ),
  # render technical replicates quarto document
  tar_quarto(
    name = technical_replicates,
    path = "analysis/02-technical_replicates.qmd",
    quiet = TRUE
  ),
  # time to join technical replicates
  tar_target(
    name = merged_metadata,
    command = merge_metadata(metadata)
  ),
  tar_target(
    name = merged_phyloseq_object,
    command = merge_technical_replicates(phyloseq_object, merged_metadata)
  )
)
