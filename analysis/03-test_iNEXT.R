
library(here)
library(targets)
library(iNEXT)
library(ggplot2)
library(future)
library(future.apply)
library(dplyr)

# determine the number of cores available
plan(multisession, workers = availableCores())

# load data
store_path <- here::here("_targets")
targets::tar_load(otu_table_matrix, store = store_path)

# combine iNEXT results lists
combine_iNEXT_results <- function(iNEXT_list) {
  # get the sample names
  assemblages <- names(iNEXT_list)

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
    # Estrazione dell'oggetto iNEXT corrente
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

# test data
## iNext
# abundance_list <- list(
#   Campione1 = c(10, 20, 30, 40, 50),
#   Campione2 = c(5, 15, 25, 35, 45),
#   Campione3 = c(8, 18, 28, 38, 48)
# )

# transform my otu table into abundance list
abundance_list <- split(as.data.frame(otu_table_matrix), rownames(otu_table_matrix))

parallel_iNEXT_results <- future_lapply(abundance_list, function(x) {
  tryCatch({
    iNEXT(as.numeric(x), q = 0, datatype = "abundance", nboot = 100)
  }, error = function(e) {
    message("Error in sample: ", e)
    return(NULL)
  })
}, future.seed = TRUE)

# merge the results
combined_iNEXT_results <- combine_iNEXT_results(parallel_iNEXT_results)

# Save the combined iNEXT object as RDS
saveRDS(combined_iNEXT_results, file = here("combined_iNEXT_results.rds"))

# Create the plot and assign it to a variable
rarefaction_plot <- ggiNEXT(combined_iNEXT_results, type = 1, se = TRUE) +
  labs(title = "Curve di Rarefazione per Tutti i Campioni",
       x = "Numero di Individui",
       y = "Ricchezza delle Specie") +
  theme_bw()

# Save the plot to a file
ggsave(
  filename = "rarefaction_plot.png",   # File name with extension
  plot = rarefaction_plot,             # Plot object
  width = 10,                          # Width in inches
  height = 6,                          # Height in inches
  units = "in",                        # Units for width and height
  dpi = 300                            # Resolution in dots per inch
)
