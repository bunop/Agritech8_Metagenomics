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

#' Compare observed/extrapolated richness at a reference depth with the
#' asymptotic Chao1 estimate, for a list of per-sample iNEXT results
#'
#' @description
#' For each sample, extracts the richness (`qD`) at the closest available
#' sample size to `depth_ref` from `iNextEst$size_based`, and the asymptotic
#' Chao1 estimate (species richness, `q = 0`) from `AsyEst`. The ratio of the
#' two gives a "sampling completeness" coverage ratio, flagged when below 90%.
#'
#' @param iNEXT_list A list of raw `iNEXT()` outputs, one per sample (as
#'   produced by `pattern = map()` over a per-sample abundance list).
#' @param assemblages Character vector of sample names, same length and order
#'   as `iNEXT_list`.
#' @param depth_ref Numeric reference sequencing depth (e.g. the rarefaction
#'   depth used elsewhere in the pipeline).
#'
#' @return A tibble with one row per `Assemblage`, columns `Assemblage`,
#'   `SampleSize_used`, `Observed_or_estimated_richness`,
#'   `Asymptotic_Chao1_estimate`, `Coverage_ratio`, `Flag_below_90pct`,
#'   sorted by `Coverage_ratio`.
#' @export
calculate_sampling_completeness <- function(iNEXT_list, abundance_list, assemblages, depth_ref) {
  results <- purrr::pmap_dfr(
    list(iNEXT_list, abundance_list, assemblages),
    function(iNEXT_item, abundance_vec, assemblage) {
      # use iNEXT::estimateD() for exact interpolation at the reference depth
      # instead of rounding to the nearest knot
      richness_at_depth <- iNEXT::estimateD(
        abundance_vec,
        q = 0,
        datatype = "abundance",
        base = "size",
        level = depth_ref,
        nboot = 0  # skip bootstrap for speed (already have it from iNEXT)
      )

      # asymptotic Chao1 estimate (species richness, q = 0)
      asy_richness <- iNEXT_item$AsyEst[iNEXT_item$AsyEst$Diversity == "Species richness", ]

      # extract F1 and F2 (singleton and doubleton counts) from DataInfo
      # for diagnosing low coverage samples
      data_info <- iNEXT_item$DataInfo
      f1 <- data_info$f1
      f2 <- data_info$f2

      tibble::tibble(
        Assemblage = assemblage,
        SampleSize_used = richness_at_depth$m,
        Observed_or_estimated_richness = richness_at_depth$qD,
        Asymptotic_Chao1_estimate = asy_richness$Estimator,
        F1_singletons = f1,
        F2_doubletons = f2,
        F1_F2_ratio = f1 / f2  # ratio for diagnostic: high = many rare species
      )
    }
  )

  results <- results %>%
    dplyr::mutate(
      Coverage_ratio = round(Observed_or_estimated_richness / Asymptotic_Chao1_estimate, 3),
      Flag_below_90pct = Coverage_ratio < 0.90
    ) %>%
    dplyr::arrange(Coverage_ratio)

  return(results)
}
