
# make rarefaction as discussed by Schloss
calculate_distance_matrix <- function(
    otu_table_matrix, min_sequencing_depth, iterations = 1000, dmethod = "bray") {
  distance_matrix <- vegan::avgdist(
    otu_table_matrix,
    sample = min_sequencing_depth,
    iterations = 1000,
    dmethod = dmethod
  )

  return(distance_matrix)
}

# calculate good's coverage as described by Schloss
calculate_coverage_stats <- function(otu_table_matrix) {
  # create a shared table to estimate the sequencing depth
  shared <- otu_table_matrix %>%
    tidyr::as_tibble(rownames = "Group") %>%
    tidyr::pivot_longer(-Group) %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(total = sum(value)) %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(total = sum(value)) %>%
    dplyr::filter(total != 0) %>%
    dplyr::ungroup() %>%
    dplyr::select(-total)

  # calculate good's coverage
  coverage_stats <- shared %>%
    dplyr::group_by(Group) %>%
    dplyr::summarize(
      n_seqs = sum(value),
      n_sings = sum(value == 1),
      goods = 100*(1 - n_sings / n_seqs))

  return(coverage_stats)
}

## calculate anova on a received function
calculate_permanova <- function(distance_matrix, metadata, columns, method = "bray", by = "term") {
  permanova_formula <- as.formula(paste0("distance_matrix ~ ", paste(columns, collapse = " + ")))
  adonis_result <- vegan::adonis2(
    permanova_formula,
    data = metadata,
    method = method,
    by = by
  )
  return(adonis_result)
}

## Pairwise permanova
calculate_pairwise_permanova <- function(distance_matrix, metadata, columns, method = "bray", by = "term") {
  permanova_formula <- as.formula(paste0("distance_matrix ~ ", paste(columns, collapse = " + ")))

  # do the pairwise permanova
  pairwise_results <- pairwise.adonis2(
    permanova_formula,
    data = metadata,
    method = method,
    by = by
  )

  # Initialize an empty list to store summaries
  summary_list <- list()

  # Loop through each pairwise comparison result
  for (comparison in names(pairwise_results)) {
    if (comparison != "parent_call") {
      # Extract the result for the current comparison
      result <- pairwise_results[[comparison]]

      # Add a column for the comparison name
      result$Comparison <- comparison

      # Append to the summary list
      summary_list[[comparison]] <- result
    }
  }

  # Combine all summaries into a single data frame
  summary_df <- do.call(rbind, summary_list)

  # Reset row names
  rownames(summary_df) <- NULL

  # Adjust p-values
  summary_df$P.adj <- p.adjust(summary_df$`Pr(>F)`, method = "BH")

  return(summary_df)
}

## calculate PCoA
calculate_pcoa <- function(distance_matrix, metadata) {
  # Calculate the PCoA
  pcoa_results <- ape::pcoa(distance_matrix, correction = "cailliez")

  # Create a data frame for plotting
  pcoa_df <- data.frame(pcoa_results$vectors[, 1:2])
  colnames(pcoa_df) <- c("PCoA1", "PCoA2")
  pcoa_df$sampleID <- rownames(pcoa_df)
  pcoa_df <- dplyr::inner_join(pcoa_df, metadata, by = "sampleID")

  # Extract the percent explained variance
  percent_explained <- pcoa_results$values$Rel_corr_eig * 100

  return(
    list(
      pcoa_results=pcoa_results,
      pcoa_df=pcoa_df,
      percent_explained=percent_explained
    )
  )
}

# do the NMDS
calculate_nmds <- function(distance_matrix, metadata, distance = "bray") {
  nmds_result <- vegan::metaMDS(distance_matrix, distance)

  nmds_tb <- vegan::scores(nmds_result) %>%
    # tibbles has not rownames, so we need to add them as a new column
    as_tibble(rownames = "sampleID") %>%
    inner_join(metadata, by = "sampleID")

  return(
    list(
      nmds_result=nmds_result,
      nmds_tb=nmds_tb
    )
  )
}

# do the rarefaction curve for each samples
calculate_rarecurve <- function(otu_table_matrix) {
  rarecurve_df <- vegan::rarecurve(otu_table_matrix, step=50, cex=0.5, tidy = TRUE)
  return(rarecurve_df)
}

#' Agglomerate by taxa and sort by sample order
#'
#' @description
#' Agglomerate the phyloseq object at a given taxonomic rank and return a
#' long-format data frame with the relative abundances of the agglomerated taxa
#' suitable for plotting with ggplot2.
#'
#' @param phyloseq_object a phyloseq object
#' @param taxrank the taxonomic rank at which to agglomerate the phyloseq object
#' ("Family", "Genus", "Species", ...)
#' @param sample_order a vector of sample IDs to specify the order in which
#' samples should be plotted
#'
#' @return a long-format data frame with the relative abundances of the
#' agglomerated taxa
#' @export
#'
#' @examples
#' melted_data <- agglomerate_by_taxa(
#'   sorted_phyloseq_object,
#'   taxrank = "Phylum",
#'   sample_order = desired_order)
agglomerate_by_taxa <- function(phyloseq_object, taxrank, sample_order = NULL) {
  # reorder phyloseq object if a specific sample order is provided
  if (!is.null(sample_order)) {
    phyloseq_object <- microViz::ps_reorder(phyloseq_object, sample_order)
  }

  # Agglomerate at the chosen taxonomic level
  physeq_aggl <- phyloseq::tax_glom(phyloseq_object, taxrank)

  # Transform to relative abundances
  physeq_aggl_rel <- phyloseq::transform_sample_counts(physeq_aggl, function(x) x / sum(x))

  # Melt the phyloseq object into a long-format data frame
  melted_data <- phyloseq::psmelt(physeq_aggl_rel)

  # If a specific sample order is provided, set the factor levels accordingly
  if (!is.null(sample_order)) {
    melted_data$sampleID <- factor(melted_data$sampleID, levels = sample_order)
  }

  return(melted_data)
}

## convert a phyloseq to ampvis2
phyloseq_to_ampvis2 <- function(phyloseq_object) {
  return(ampvis2::amp_load(phyloseq_object))
}

# Function to rarefy and calculate alpha diversity metrics
# This function is a wrapper around the rarefy_even_depth and estimate_richness functions
# from the phyloseq package. It make a single sub sampling like Schloss means
rarefy_alpha <- function(physeq_obj, depth, measures = NULL, rngseed = FALSE) {
  # Rarefy the phyloseq object
  physeq_rarefied <- phyloseq::rarefy_even_depth(
    physeq_obj,
    sample.size = depth,
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
