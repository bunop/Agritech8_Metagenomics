
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
    as_tibble(rownames = "Group") %>%
    pivot_longer(-Group) %>%
    group_by(Group) %>%
    mutate(total = sum(value)) %>%
    group_by(name) %>%
    mutate(total = sum(value)) %>%
    filter(total != 0) %>%
    ungroup() %>%
    select(-total)

  # calculate good's coverage
  coverage_stats <- shared %>%
    group_by(Group) %>%
    summarize(
      n_seqs = sum(value),
      n_sings = sum(value == 1),
      goods = 100*(1 - n_sings / n_seqs))

  return(coverage_stats)
}

## calculate anova on a received function
calculate_anova <- function(phyloseq_object, distance_matrix) {
  adonis_result <- vegan::adonis2(
    distance_matrix ~ sample_group + technical_rep,
    data = methods::as(phyloseq::sample_data(phyloseq_object), "data.frame"),
    by = "margin"
  )
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


## merge metadata relying on technical replicates
merge_metadata <- function(metadata) {
  # merging those columns make no sense
  tmp <- metadata %>%
    dplyr::select(-sample_number, -replica)

  # Identify numeric and non-numeric columns
  numeric_cols <- names(tmp)[sapply(tmp, is.numeric)]
  non_numeric_cols <- names(tmp)[!sapply(tmp, is.numeric)]

  # Remove 'sample_group' from non-numeric columns since will be
  # removed when grouping by 'sample_group'
  non_numeric_cols <- setdiff(non_numeric_cols, "sample_group")

  merged_metadata <- tmp %>%
    dplyr::group_by(sample_group) %>%
    dplyr::summarize(
      dplyr::across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(all_of(non_numeric_cols), ~ first(.x)),
      .groups = "drop"
    ) %>%
    dplyr::select(sample_name, everything()) %>%
    as.data.frame()

  # Set rownames
  rownames(merged_metadata) <- merged_metadata$sample_name

  return(merged_metadata)
}

# merge samples by replicate group, summing ASV abundances
# and divide the ASV abundances by the number of replicates to get the mean
merge_technical_replicates <- function(phyloseq_object, merged_metadata) {
  # Merge samples by replicate group, summing ASV abundances
  merged_phyloseq_object <- merge_samples(phyloseq_object, group = "sample_name")

  # replace the metadata
  phyloseq::sample_data(merged_phyloseq_object) <- phyloseq::sample_data(merged_metadata)

  # Count the number of replicates in each group
  replicate_counts <- table(sample_data(phyloseq_object)$sample_name)

  # Divide the ASV abundances by the number of replicates to get the mean
  otu_mat <- t(phyloseq::otu_table(merged_phyloseq_object))
  otu_mat <- sweep(otu_mat, 2, replicate_counts, FUN = "/")
  otu_mat_floored <- floor(otu_mat)

  # Convert the matrix back to an otu_table object
  phyloseq::otu_table(merged_phyloseq_object) <- otu_table(
    otu_mat_floored,
    taxa_are_rows = phyloseq::taxa_are_rows(
      phyloseq::otu_table(phyloseq_object)
    )
  )

  # prune empty taxa
  merged_phyloseq_object <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(merged_phyloseq_object) > 0, merged_phyloseq_object)

  return(merged_phyloseq_object)
}


