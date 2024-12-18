
# load metadata and fix columns
load_metadata <- function(...) {
  # take all the arguments at once
  args <- list(...)
  metadata <- readr::read_delim(here::here(args))

  # get a new name for sample
  metadata <- metadata %>%
    # Extract the last character of sampleID to create technical_rep
    mutate(technical_rep = stringr::str_sub(sampleID, -1)) %>%
    # Convert technical_rep to a factor
    mutate(technical_rep = as.factor(technical_rep)) %>%
    # Remove the last character from sampleID to create sample_name
    mutate(sample_name = stringr::str_sub(sampleID, 1, -2)) %>%
    # Ensure sample_group is a factor
    mutate(sample_group = as.factor(sample_group)) %>%
    # Group by sample_group
    group_by(sample_group) %>%
    # Assign the first sample_name within each group
    mutate(sample_name = first(sample_name)) %>%
    # Ungroup the DataFrame
    ungroup() %>%
    # Reorder columns to place sampleID and sample_name at the beginning
    select(sampleID, sample_name, technical_rep, everything()) %>%
    # Convert to a standard data frame
    as.data.frame()

  # set rownames
  rownames(metadata) <- metadata$sampleID

  # Create a new grouping variable using only the last 20 characters of 'condition' and 'date'
  metadata$date_condition <- as.factor(with(
    metadata,
    paste0(
      gsub(" ", "", substr(condition, nchar(condition) - 12, nchar(condition))),
      "_",
      gsub(" ", "", date)
    )
  ))

  return(metadata)
}

# load and update phyloseq object
load_phyloseq <- function(metadata, ...) {
  args <- list(...)

  # load the phyloseq object
  phyloseq_object <- readRDS(here::here(args))

  # replace the metadata
  phyloseq::sample_data(phyloseq_object) <- metadata

  return(phyloseq_object)
}

# get the otu table from a phyloseq object
get_otu_table <- function(phyloseq_object) {
  otu_table_matrix <- as(phyloseq::otu_table(phyloseq_object), "matrix")

  # check if the taxa are rows
  if (taxa_are_rows(phyloseq_object)) {
    otu_table_matrix <- t(otu_table_matrix)
  }

  return(otu_table_matrix)
}

# get information on samples
get_samples_data <- function(phyloseq_object) {
  # get the sequencing depth
  sequencing_depth <- phyloseq::sample_sums(phyloseq_object)

  # get data from samples
  samples_data <- data.frame(phyloseq::sample_data(phyloseq_object))

  # add sequencing depth to samples data
  samples_data$sequencing_depth <- sequencing_depth

  return(samples_data)
}

# load the qiime2 rarefaction tables
load_qiime_rarefaction <- function(rarefaction_csv, metadata, alpha_metric, column_name) {
  # Load the rarefaction CSV file
  data <- readr::read_csv(here::here(rarefaction_csv))

  # Merge the metadata with the rarefaction data (inner join)
  # get rid of the common columns from the first dataframe
  common_columns <- intersect(names(data), names(metadata))
  data <- data  %>%
    select(-all_of(common_columns)) %>%
    dplyr::inner_join(metadata, by = c("sample-id" = "sampleID"))

  # Reshape the data to long format
  long_data <- data %>%
    pivot_longer(
      cols = starts_with("depth-"),
      names_to = c("depth", "iteration"),
      names_sep = "_iter-",
      values_to = alpha_metric
    ) %>%
    mutate(
      depth = as.numeric(gsub("depth-", "", depth))
    )

  # Convert alpha_metric and column_name to symbols for tidy evaluation
  alpha_metric_sym <- sym(alpha_metric)
  column_name_sym <- sym(column_name)

  # Group and summarize data with 95% CI
  summary_data <- long_data %>%
    # unquote the symbols
    group_by(depth, !!column_name_sym) %>%
    summarise(
      mean_value = mean(!!alpha_metric_sym, na.rm = TRUE),
      se_value = sd(!!alpha_metric_sym, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      ci_lower = mean_value - 1.96 * se_value,
      ci_upper = mean_value + 1.96 * se_value
    ) # 95% CI

  return(summary_data)
}

# select samples based on a column in metadata table
select_samples <- function(metadata, values, column_name="date_condition") {
  samples_to_keep <- metadata %>%
    dplyr::filter(!!sym(column_name) %in% values) %>%
    dplyr::select(sampleID)

  return(samples_to_keep)
}

# subsetting phyloseq object
subset_phyloseq <- function(phyloseq_object, samples_to_keep) {
  phyloseq_subset <- phyloseq::prune_samples(
    sample_names(phyloseq_object) %in% t(samples_to_keep), phyloseq_object)

  # Identify taxa with non-zero total abundance
  non_zero_taxa <- phyloseq::taxa_sums(phyloseq_subset) > 0

  # Prune taxa with zero abundance
  phyloseq_pruned <- phyloseq::prune_taxa(non_zero_taxa, phyloseq_subset)

  return(phyloseq_pruned)
}

# get metadata from a phyloseq object
get_metadata <- function(phyloseq_object, order_by="sample_number") {
  metadata <- phyloseq::sample_data(phyloseq_object) %>%
    as_tibble() %>%
    dplyr::arrange(!!sym(order_by))
  return(metadata)
}
