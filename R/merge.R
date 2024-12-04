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
