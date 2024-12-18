
#' Create Krona input files and return the command to run Krona
#'
#' This function creates Krona input files from a phyloseq object and returns
#' the command to run Krona. This function was adapted from the
#' \href{https://github.com/cpauvert/psadd}{cpauvert/psadd package}
#'
#' @param physeq \code{\link{phyloseq-class}} with a \code{\link{taxonomyTable-class}}
#' @param output A \code{\link{character}} stating the output filename for Krona Chart and the directory in which Krona files will be created.
#' @param variable A \code{\link{character}} indicating which variable from sample data to use.
#' @param trim Should spaces and brackets be converted to underscore automatically.
#'
#' @import phyloseq
#' @importFrom utils browseURL write.table
#' @export
#'
#' @references \href{https://github.com/marbl/Krona/wiki/KronaTools}{KronaTools}.
#' @examples
#' \dontrun{
#' require(phyloseq)
#' data(GlobalPatterns)
#' get_krona_cmd(GlobalPatterns, "GP-krona", "SampleType")# issues with brackets
#' get_krona_cmd(GlobalPatterns, "GP-krona", "SampleType", trim=T)
#' }
get_krona_cmd <- function(physeq, output, variable, trim=F){
  if (is.null(phyloseq::tax_table(physeq))) {
    stop("No taxonomy table available.")
  }

  columns <- colnames(sample_data(physeq))
  if (!variable %in% columns) {
    stop(glue::glue(
      "There's no '{variable}' in the sample data. ",
      "Available variables are: {glue::glue_collapse(columns, sep = ', ')}"
    ))
  }

  if (trim == FALSE) {
    spec.char <- grepl(" |\\(|\\)", as(phyloseq::sample_data(physeq),"data.frame")[,variable])
    if (sum(spec.char > 0)) {
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }

  # Melt the OTU table and merge associated metadata
  melted_physeq <- phyloseq::psmelt(physeq)

  # Define the columns of interest
  columns_of_interest <- c("Abundance", {{ variable }}, rank_names(physeq)[1:7])

  # Check if 'Species_exact' exists in the melted data frame
  if ("Species_exact" %in% colnames(melted_physeq)) {
    columns_of_interest <- c(columns_of_interest, "Species_exact")
  }

  df <- melted_physeq %>%
    # Fetch only column of interest
    dplyr::select(dplyr::all_of(columns_of_interest)) %>%
    # Ensure there are no spaces left in the variable of interest
    dplyr::mutate({{ variable }} := gsub(" |\\(|\\)", "", .data[[ variable ]])) %>%
    # Convert the field of interest to a factor
    dplyr::mutate({{ variable }} := as.factor(.data[[ variable ]]))

  # Create a directory for krona files
  dir.create(output)

  # For each level of the Description variable
  # Abundance and taxonomic assignations for each OTU are fetched
  # and written to a file that would be processed by Krona.
  for (lvl in levels(df[[variable]])) {
    subset <- df %>%
      dplyr::filter(.data[[variable]] == lvl & Abundance != 0) %>%
      dplyr::select(-{{ variable }})

    write.table(
      subset,
      file = paste0(output, "/" , lvl, ".txt"),
      sep = "\t", row.names = F, col.names = F, na = "", quote = F)
  }

  # Arguments for Krona command
  # taxonomic file and their associated labels.
  krona_args <- paste(output, "/", levels(df[[variable]]), ".txt", sep = "", collapse = " ")

  # Add html suffix to output
  output <- paste(output, ".html", sep = "")

  # return krona command
  return(paste(
    "ktImportText",
    krona_args,
    "-o", output,
    sep = " "))
}
