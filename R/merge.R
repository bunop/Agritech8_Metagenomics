#' Calculate Spearman correlations for technical replicates
#'
#' @param phyloseq_obj A phyloseq object
#' @param metadata A data frame with sample metadata
#' @param group_column Column name identifying biological replicates
#' @return A data frame with correlation statistics per group
calculate_technical_replicate_correlations <- function(phyloseq_obj, metadata, group_column = "sample_group") {

  otu_mat <- as(phyloseq::otu_table(phyloseq_obj), "matrix")
  if (!phyloseq::taxa_are_rows(phyloseq_obj)) {
    otu_mat <- t(otu_mat)
  }

  # calculate correlations per group
  unique_groups <- unique(metadata[[group_column]])

  purrr::map_dfr(unique_groups, function(grp) {
    sample_ids <- metadata |>
      dplyr::filter(.data[[group_column]] == grp) |>
      dplyr::pull(sampleID)

    if (length(sample_ids) < 2) {
      return(tibble::tibble(
        sample_group = grp,
        n_replicates = length(sample_ids),
        mean_spearman = NA_real_,
        min_spearman = NA_real_,
        max_spearman = NA_real_
      ))
    }

    # Correlation matrix
    cor_mat <- cor(otu_mat[, sample_ids, drop = FALSE], method = "spearman")

    # Extract upper triangle values (pairwise correlations)
    cor_values <- cor_mat[upper.tri(cor_mat)]

    tibble::tibble(
      sample_group = grp,
      n_replicates = length(sample_ids),
      mean_spearman = mean(cor_values),
      min_spearman = min(cor_values),
      max_spearman = max(cor_values),
      sd_spearman = sd(cor_values)
    )
  })
}

#' Calculate within-group distances for technical replicates
#'
#' @param distance_matrix A distance matrix object
#' @param metadata A data frame with sample metadata
#' @param group_column Column identifying biological replicates
#' @return A data frame with distance statistics per group
calculate_technical_replicate_distances <- function(distance_matrix, metadata, group_column = "sample_group") {

  dist_mat <- as.matrix(distance_matrix)

  unique_groups <- unique(metadata[[group_column]])

  purrr::map_dfr(unique_groups, function(grp) {
    sample_ids <- metadata |>
      dplyr::filter(.data[[group_column]] == grp) |>
      dplyr::pull(sampleID)

    if (length(sample_ids) < 2) {
      return(tibble::tibble(
        sample_group = grp,
        n_replicates = length(sample_ids),
        mean_distance = NA_real_,
        median_distance = NA_real_,
        min_distance = NA_real_,
        max_distance = NA_real_
      ))
    }

    # Subset distance matrix for the group
    group_dist <- dist_mat[sample_ids, sample_ids]

    # Extract upper triangle values (pairwise distances)
    dist_values <- group_dist[upper.tri(group_dist)]

    tibble::tibble(
      sample_group = grp,
      n_replicates = length(sample_ids),
      mean_distance = mean(dist_values),
      median_distance = median(dist_values),
      min_distance = min(dist_values),
      max_distance = max(dist_values),
      sd_distance = sd(dist_values)
    )
  })
}

## merge metadata relying on technical replicates
merge_metadata <- function(metadata, column = "sample_group") {
  # merging those columns make no sense
  tmp <- metadata %>%
    dplyr::select(-sample_number, -replica)

  # Identify numeric and non-numeric columns
  numeric_cols <- names(tmp)[sapply(tmp, is.numeric)]
  non_numeric_cols <- names(tmp)[!sapply(tmp, is.numeric)]

  # Remove _column_ from non-numeric columns since will be
  # removed when grouping by 'sample_group'
  non_numeric_cols <- setdiff(non_numeric_cols, column)

  merged_metadata <- tmp %>%
    dplyr::group_by(!!sym(column)) %>%
    dplyr::summarize(
      dplyr::across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)),
      dplyr::across(all_of(non_numeric_cols), ~ first(.x)),
      .groups = "drop"
    ) %>%
    dplyr::select(sample_name, everything()) %>%
    as.data.frame()

  # Set rownames
  rownames(merged_metadata) <- merged_metadata[[column]]

  # override sampleID column
  merged_metadata$sampleID <- merged_metadata[[column]]

  return(merged_metadata)
}

# merge samples by replicate group, summing ASV abundances
# and divide the ASV abundances by the number of replicates to get the mean
merge_technical_replicates <- function(phyloseq_object, merge_metadata, column = "sample_names") {
  # Merge samples by replicate group, calculating mean ASV abundances
  merged_phyloseq_object <- phyloseq::merge_samples(phyloseq_object, group = column, fun = mean)

  # Update sample metadata ensuring rownames match
  # phyloseq uses rownames to match samples to metadata
  merge_metadata_aligned <- merge_metadata[sample_names(merged_phyloseq_object), ]
  rownames(merge_metadata_aligned) <- sample_names(merged_phyloseq_object)

  # Add sampleID column with same values as the grouping column
  merge_metadata_aligned$sampleID <- merge_metadata_aligned[[column]]

  sample_data(merged_phyloseq_object) <- phyloseq::sample_data(merge_metadata_aligned)

  # Prune empty taxa
  merged_phyloseq_object <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(merged_phyloseq_object) > 0,
    merged_phyloseq_object
  )

  return(merged_phyloseq_object)
}
