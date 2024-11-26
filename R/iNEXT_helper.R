# combine iNEXT results lists
combine_iNEXT_results <- function(iNEXT_list, abundance_list) {
  # get the sample names
  assemblages <- names(abundance_list)

  # create an empty iNEXT object
  combined_iNEXT <- list(
    DataInfo = data.frame(),
    iNextEst = list(size_based = data.frame(), coverage_based = data.frame()),
    AsyEst = data.frame()
  )

  # setting the class
  class(combined_iNEXT) <- "iNEXT"

  # cicle on iNEXT_list items
  for (index in 1:length(iNEXT_list)) {
    # collect the current iNEXT item
    iNEXT_item <- iNEXT_list[[index]]
    assemblage <- assemblages[index]

    # skip empty items
    if (is.null(iNEXT_item)) {
      next
    }

    # get list elements for each item
    DataInfo <- iNEXT_item$DataInfo
    iNextEst <- iNEXT_item$iNextEst
    AsyEst <- iNEXT_item$AsyEst

    # fix the assemblage names
    DataInfo$Assemblage <- assemblage
    iNextEst$size_based$Assemblage <- assemblage
    iNextEst$coverage_based$Assemblage <- assemblage

    # deal with AsyEst
    n_rows <- nrow(AsyEst)
    tmp <- data.frame(
      Assemblage <- rep(assemblage, n_rows),
      Diversity <- c("Species richness", "Shannon diversity", "Simpson diversity"),
      Observed <- AsyEst$Observed,
      Estimator <- AsyEst$Estimator,
      `s.e.` <- AsyEst$`Est_s.e.`,
      LCL <- AsyEst$`95% Lower`,
      UCL <- AsyEst$`95% Upper`
    )

    # combine the changed elements into one object
    combined_iNEXT$DataInfo <- rbind(combined_iNEXT$DataInfo, DataInfo)
    combined_iNEXT$iNextEst$size_based <- rbind(combined_iNEXT$iNextEst$size_based, iNextEst$size_based)
    combined_iNEXT$iNextEst$coverage_based <- rbind(
      combined_iNEXT$iNextEst$coverage_based, iNextEst$coverage_based)
    combined_iNEXT$AsyEst <- rbind(combined_iNEXT$AsyEst, tmp)
  }

  # Return the combined iNEXT object
  return(combined_iNEXT)
}
