
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

# average good's coverage per sample across rarefaction iterations
summarize_coverage_stats <- function(coverage_stats) {
  coverage_stats %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(dplyr::across(c(n_seqs, n_sings, goods), mean), .groups = "drop")
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

#' Rarefy a phyloseq object
#'
#' Rarefy a phyloseq object to a specified depth. This function uses the
#' `rarefy_even_depth` function from the `phyloseq` package to perform
#' rarefaction. The function allows for the option to replace
#' samples with fewer reads than the specified depth with
#' samples that have more reads than the specified depth.
#' This is useful for normalizing the sequencing depth across samples.
#'
#' @param physeq_obj a phyloseq object
#' @param depth the desired sequencing depth to rarefy to (the minimum sample depth
#' across all samples in the phyloseq object as a default)
#' @param replace a logical value indicating whether to replace samples with fewer reads
#' @param trimOTUS a logical value indicating whether to trim OTUs with zero reads
#' @param rngseed a random seed for reproducibility. Set FALSE to use the current
#' @param verbose a logical value indicating whether to print verbose output
#'
#' @returns a rarefied phyloseq object
#' @export
#'
#' @examples
#' rarefied_phyloseq_object <- rarefy_phyloseq_object(
#'  physeq_obj = phyloseq_object,
#'  depth = 1000,
#'  replace = FALSE,
#'  trimOTUS = TRUE,
#'  rngseed = 123,
#'  verbose = TRUE)
rarefy_phyloseq_object <- function(physeq_obj, depth = NULL, replace = FALSE, trimOTUS = TRUE, rngseed = FALSE, verbose = FALSE) {
  # calculate the minimum sample depth if not provided
  if (is.null(depth)) {
    depth <- min(phyloseq::sample_sums(physeq_obj))
  }

  # Rarefy the phyloseq object
  physeq_rarefied <- phyloseq::rarefy_even_depth(
    physeq_obj,
    sample.size = depth,
    replace = FALSE,  # Sample without replacement
    trimOTUs = TRUE,  # Trim OTUs with zero reads
    rngseed = rngseed,
    verbose = FALSE
  )
  return(physeq_rarefied)
}

## convert a phyloseq to ampvis2
phyloseq_to_ampvis2 <- function(phyloseq_object) {
  return(ampvis2::amp_load(phyloseq_object))
}
