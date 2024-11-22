
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
