
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
