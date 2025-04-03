
# do the rarefaction curve for each samples
calculate_rarecurve <- function(otu_table_matrix) {
  rarecurve_df <- vegan::rarecurve(otu_table_matrix, step = 50, cex = 0.5, tidy = TRUE)
  return(rarecurve_df)
}

# Function to rarefy and calculate alpha diversity metrics
# This function is a wrapper around the rarefy_even_depth and estimate_richness functions
# from the phyloseq package. It make a single sub sampling like Schloss means
rarefy_alpha <- function(physeq_obj, depth, measures = NULL, rngseed = FALSE) {
  # Rarefy the phyloseq object
  physeq_rarefied <- phyloseq::rarefy_even_depth(
    physeq_obj,
    sample.size = depth,
    # sample without replacement
    replace = FALSE,
    # trim OTUs with zero reads
    trimOTUs = TRUE,
    rngseed = rngseed,
    verbose = FALSE
  )

  # Calculate alpha diversity metrics
  alpha_div <- phyloseq::estimate_richness(physeq_rarefied, measures = measures)
  alpha_div <- tibble::rownames_to_column(
    alpha_div,
    var = "sampleID"
  )
  return(alpha_div)
}

summarize_rarefactions <- function(samples_rarefaction, metadata, by_column = "sampleID") {
  summary_results <- samples_rarefaction %>%
    group_by(!!sym(by_column)) %>%
    summarise(across(everything(), list(mean = mean, sd = sd)))

  # Merge with summary results
  samples_rarefaction <- left_join(summary_results, metadata, by = by_column)
  return(samples_rarefaction)
}

reshape_rarefaction_data <- function(data, by_column) {
  # Reshape data to long format
  data_long <- data %>%
    # Select relevant columns: sampleID, by_column, and those ending with '_mean' or '_sd'
    dplyr::select(sampleID, !!sym(by_column), ends_with("_mean"), ends_with("_sd")) %>%
    tidyr::pivot_longer(
      cols = ends_with("_mean") | ends_with("_sd"),
      names_to = c("Metric", "Measure"),
      names_pattern = "(.*)_(mean|sd)$",
      values_to = "Value"
    ) %>%
    tidyr::pivot_wider(
      names_from = Measure,
      values_from = Value
    )

  return(data_long)
}

# Calculate the Kruskall-Wallis test
calculate_kruskal_wallis <- function(rarefaction_results, alpha_metric, column_name) {
  # define a new formula
  kw_formula <- as.formula(paste(alpha_metric, "~", column_name))

  # Perform the Kruskal-Wallis test
  kruskal_results <- stats::kruskal.test(
    formula = kw_formula,
    data = rarefaction_results
  )

  return(kruskal_results)
}

# Dunn's Kruskal-Wallis Multiple Comparisons
calculate_dunn_test <- function(rarefaction_results, alpha_metric, column_name, method = "bh") {
  # define a new formula
  dunn_formula <- as.formula(paste(alpha_metric, "~", column_name))

  # Perform the Dunn's test
  dunn_results <- FSA::dunnTest(
    x = dunn_formula,
    data = rarefaction_results,
    method = method
  )

  return(dunn_results)
}
