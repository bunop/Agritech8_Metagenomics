
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
