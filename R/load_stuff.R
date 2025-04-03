
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

# This function sorts a metadata dataframe by a specified factor column, reordering the
# levels of that factor column according to user-defined levels.
sort_by_factor_column <- function(metadata, order_by, levels) {
  # First, convert the order_by string into a symbol.
  # A symbol is a type of object in R that represents a column name in a dataframe
  # when working with dplyr in a non-standard evaluation context.
  order_by_sym <- sym(order_by)

  # Next, reorder the metadata dataframe by mutating the order_by column to reorder
  # its factor levels, then arrange the dataframe by this newly ordered factor column
  # and sampleID.
  metadata <- metadata %>%
    # Use mutate to change the order_by column into a factor with levels ordered as
    # specified by the 'levels' vector.
    # The 'forcats::fct_relevel' function is used to reorder the levels of the factor.
    # The '!!' operator is used to force evaluation of the symbol created earlier,
    # so that it is interpreted as a column name.
    mutate(!!order_by_sym := forcats::fct_relevel(!!order_by_sym, !!!levels)) %>%

    # Use arrange to sort the dataframe by the order_by column (now a factor with
    # the right order) and then by sampleID.
    arrange(!!order_by_sym, sampleID) %>%

    # Further, ensure that the 'sampleID' column is also a factor with levels
    # corresponding to the order in which they appear in the dataframe after sorting.
    # This can be useful for keeping the order consistent when using the sampleID in
    # plots or analyses later on.
    mutate(sampleID = factor(sampleID, levels = unique(sampleID)))

  # Finally, return the sorted dataframe.
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

# add tree to phyloseq object
add_tree_to_phyloseq <- function(phyloseq_object, ...) {
  # collect all the arguments in a list
  args <- list(...)

  # load the tree
  newick_file_path <- here::here(args)
  phylogenetic_tree <- ape::read.tree(newick_file_path)

  # add the tree to the phyloseq object
  phyloseq::phy_tree(phyloseq_object) <- phylogenetic_tree

  return(phyloseq_object)
}

# remove unwanted taxa from phyloseq object
prune_phyloseq <- function(phyloseq_obj, rank, items) {
  # Check if the rank provided is valid
  if (!rank %in% colnames(phyloseq::tax_table(phyloseq_obj))) {
    stop("The specified rank does not exist in the taxonomy table.")
  }

  # Ensure items is a character vector (it can be a single item or multiple items)
  if (is.vector(items) && !is.character(items)) {
    stop("Items should be a character vector.")
  }

  # Create a logical vector to identify taxa to keep (not matching any of the specified items)
  taxa_to_keep <- as.vector(!(phyloseq::tax_table(phyloseq_obj)[, rank] %in% items))

  # Use prune_taxa to remove the specified taxa
  pruned_phyloseq_obj <- phyloseq::prune_taxa(taxa_to_keep, phyloseq_obj)

  return(pruned_phyloseq_obj)
}

# order phyloseq object relying on sample ids
sort_phyloseq <- function(phyloseq_obj, sample_ids) {
  # ensure sample_ids is a character vector
  sample_ids <- as.character(sample_ids)

  # now order the samples with microViz
  sorted_phyloseq_obj <- microViz::ps_reorder(
    phyloseq_obj, sample_ids)

  return(sorted_phyloseq_obj)
}

# get the otu table from a phyloseq object
get_otu_table <- function(phyloseq_object) {
  otu_table_matrix <- as(phyloseq::otu_table(phyloseq_object), "matrix")

  # check if the taxa are rows
  if (phyloseq::taxa_are_rows(phyloseq_object)) {
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
    tidyr::pivot_longer(
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
