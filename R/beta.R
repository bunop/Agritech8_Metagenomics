
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

# work on beta dispersion
calculate_beta_dispersion <- function(distance_matrix, metadata, column) {
  # calculate beta dispersion
  beta_disp <- vegan::betadisper(distance_matrix, metadata[[column]])

  # Perform ANOVA on dispersion
  anova_result <- anova(beta_disp)

  # organize results for plotting
  beta_disp_df <- data.frame(
    Sample = names(beta_disp$distances),
    Distance = beta_disp$distances,
    Group = metadata[[column]]
  )

  return(list(
    beta_disp = beta_disp,
    anova_result = anova_result,
    beta_disp_df = beta_disp_df)
  )
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
      pcoa_results = pcoa_results,
      pcoa_df = pcoa_df,
      percent_explained = percent_explained
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
      nmds_result = nmds_result,
      nmds_tb = nmds_tb
    )
  )
}

get_distances <- function(distance_matrix, metadata, column) {
  # Convert distance matrix to a data frame
  distance_df <- as.data.frame(as.table(as.matrix(distance_matrix)))
  colnames(distance_df) <- c("Sample1", "Sample2", "Distance")

  # Merge metadata to get group information for each sample
  distance_df <- merge(distance_df, metadata, by.x = "Sample1", by.y = "sampleID")
  distance_df <- merge(distance_df, metadata, by.x = "Sample2", by.y = "sampleID", suffixes = c(".1", ".2"))

  # filter the column I need
  column.1 <- sym(paste0(column, ".1"))
  column.2 <- sym(paste0(column, ".2"))

  distance_df <- distance_df %>%
    dplyr::select(Sample1, Sample2, !!column.1, !!column.2, Distance)

  return(distance_df)
}
