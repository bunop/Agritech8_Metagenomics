
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

## convert a phyloseq to ampvis2
phyloseq_to_ampvis2 <- function(phyloseq_object) {
  return(ampvis2::amp_load(phyloseq_object))
}
