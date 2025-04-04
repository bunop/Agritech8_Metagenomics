
# make rarefaction on a single rarefied phyloseq object
calculate_distance_matrix <- function(phyloseq_obj, method = "bray") {
  otu_table_matrix <- get_otu_table(phyloseq_obj)
  distance_matrix <- as.matrix(vegan::vegdist(otu_table_matrix, method = method))
  return(distance_matrix)
}

#' Calculate Beta Dispersion and Perform ANOVA
#'
#' This function calculates beta dispersion for a given distance matrix and conducts ANOVA
#' on the resulting dispersion values based on groups defined in the metadata. It can also
#' reorder the levels of the grouping factor in the metadata if specified.
#'
#' @param distance_matrix A dissimilarity matrix (object of class 'dist') that contains the pairwise distances between samples.
#' @param metadata A data frame that includes sample information, with one column defining the grouping of samples for dispersion analysis.
#' @param column A string indicating the name of the column in `metadata` that represents the groups to compare for beta dispersion.
#' @param levels A character vector (optional) that defines the desired order of factor levels for the specified column.
#' If provided, the levels will be reordered using `forcats::fct_relevel`.
#'
#' @return A list containing:
#'   - `beta_disp`: An object of class 'betadisper' containing the beta dispersion results.
#'   - `anova_result`: The results of the ANOVA test performed on the beta dispersion values.
#'   - `beta_disp_df`: A data frame with the beta dispersion values for each sample, including the group information.
#'
#' @examples
#' # Create a distance matrix and metadata
#' distance_matrix <- vegdist(data_matrix)
#' metadata <- data.frame(SampleID = rownames(data_matrix), Group = c("A", "A", "B", "B"))
#' # Calculate beta dispersion
#' result <- calculate_beta_dispersion(distance_matrix, metadata, column = "Group", levels = c("A", "B"))
#'
#' @importFrom vegan betadisper
#' @importFrom dplyr mutate
#' @importFrom rlang sym
#' @importFrom forcats fct_relevel
#' @importFrom stats anova
calculate_beta_dispersion <- function(distance_matrix, metadata, column, levels=NULL) {
  # calculate beta dispersion
  beta_disp <- vegan::betadisper(distance_matrix, metadata[[column]])

  # Perform ANOVA on dispersion
  anova_result <- stats::anova(beta_disp)

  # organize results for plotting
  beta_disp_df <- data.frame(
    Sample = names(beta_disp$distances),
    Distance = beta_disp$distances,
    Group = beta_disp$group
  )

  if (!is.null(levels)) {
    beta_disp_df <- beta_disp_df %>%
      dplyr::mutate(Group := forcats::fct_relevel(Group, !!!levels))
  }

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
  pairwise_results <- pairwiseAdonis::pairwise.adonis2(
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
calculate_pcoa <- function(distance_matrix, metadata, dimensions = 2) {
  # Calculate the PCoA
  pcoa_results <- ape::pcoa(distance_matrix, correction = "cailliez")

  # Create a data frame for plotting
  pcoa_df <- data.frame(pcoa_results$vectors[, 1:dimensions])
  colnames(pcoa_df) <- paste("PCoA", 1:dimensions, sep = "")
  pcoa_df$sampleID <- rownames(pcoa_df)
  pcoa_df <- dplyr::inner_join(pcoa_df, metadata, by = "sampleID")

  # Determine which relative eigenvalues to use
  if ("Rel_corr_eig" %in% names(pcoa_results$values)) {
    percent_explained <- pcoa_results$values$Rel_corr_eig * 100
  } else {
    percent_explained <- pcoa_results$values$Relative_eig * 100
  }

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
    dplyr::as_tibble(rownames = "sampleID") %>%
    dplyr::inner_join(metadata, by = "sampleID")

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
  column.1 <- rlang::sym(paste0(column, ".1"))
  column.2 <- rlang::sym(paste0(column, ".2"))

  distance_df <- distance_df %>%
    dplyr::select(Sample1, Sample2, !!column.1, !!column.2, Distance)

  return(distance_df)
}
