
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

calculate_rarecurve <- function(otu_table_matrix) {
  rarecurve_df <- rarecurve(otu_table_matrix, step=50, cex=0.5, tidy = TRUE)
  return(rarecurve_df)
}
