################################################################################
#' Make Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Need the installation of kronatools on the computer ([installation instruction](https://github.com/marbl/Krona/wiki/Installing)).
#'
#' @inheritParams clean_pq
#' @param file (required) the location of the html file to save
#' @param nb_seq (logical) If true, Krona set the distribution
#'   of sequences in the taxonomy. If False, Krona set the distribution of ASVs
#'   in the taxonomy.
#' @param ranks Number of the taxonomic ranks to plot
#'   (num of the column in `tax_table` slot of your `physeq` object).
#'   Default setting plot all the ranks (argument 'All').
#' @param add_unassigned_rank (int) Add unassigned for rank
#'   inferior to 'add_unassigned_rank' when necessary.
#' @param name A name for intermediary files, Useful to name
#'   your krona result files before merging using [merge_krona()]
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows" && MiscMetabar::is_krona_installed()
#' data("GlobalPatterns", package = "phyloseq")
#' GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
#' \dontrun{
#' krona(GA, "Number.of.sequences.html")
#' krona(GA, "Number.of.ASVs.html", nb_seq = FALSE)
#' merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"))
#' }
#' @return A html file
#' @export
#' @seealso \code{\link{merge_krona}}
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [Krona](https://github.com/marbl/Krona) if
#'   you use this function.
krona <-
  function(physeq,
           file = "krona.html",
           nb_seq = TRUE,
           ranks = "All",
           add_unassigned_rank = 0,
           name = NULL) {
    if (ranks[1] == "All") {
      ranks <- seq_along(physeq@tax_table[1, ])
    }

    df <- data.frame(unclass(physeq@tax_table[, ranks]))
    df$ASVs <- rownames(physeq@tax_table)

    if (is.null(name)) {
      if (nb_seq) {
        name <- "Number.of.sequences_temp"
      } else {
        name <- "Number.of.ASVs_temp"
      }
    }

    if (nb_seq) {
      df$nb_seq <- taxa_sums(physeq)
    } else {
      df$nb_otu <- rep(1, length(taxa_sums(physeq)))
    }

    df <- df[c(ncol(df), 2:ncol(df) - 1)]
    res <-
      lapply(split(df, seq_along(physeq@tax_table[, 1])), function(x) {
        as.vector(as.matrix(x))[!is.na(unlist(x))]
      })

    res <-
      lapply(res, function(x) {
        if (length(x) < add_unassigned_rank) {
          x <-
            c(x, "unassigned")[c(seq_len(x) - 1, length(x) + 1, length(x))]
        } else {
          x
        }
      })

    interm_txt <- paste(tempdir(), "/", name, ".html", sep = "")

    lapply(res,
      cat,
      "\n",
      file = interm_txt,
      append = TRUE,
      sep = "\t"
    )

    cmd <- paste("ktImportText ", interm_txt, " -o ", file, sep = "")
    # system(command = cmd)
    # system(command = paste("rm", interm_txt))
    return(cmd)
  }

###############################################################################
#' Merge Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Need the installation of kronatools on the computer
#' ([installation instruction](https://github.com/marbl/Krona/wiki/Installing)).
#'
#' Function merge_krona allows merging multiple html files in one interactive
#' krona file
#'
#' Note that you need to use the name args in `krona()` function before `merge_krona()`
#' in order to give good name to each krona pie in the output.
#' @param files (required) path to html files to merged
#' @param output path to the output file
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows" && MiscMetabar::is_krona_installed()
#' \dontrun{
#' data("GlobalPatterns", package = "phyloseq")
#' GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
#' krona(GA, "Number.of.sequences.html", name = "Nb_seq_GP_acidobacteria")
#' krona(GA, "Number.of.ASVs.html", nb_seq = FALSE, name = "Nb_asv_GP_acidobacteria")
#' merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"), "mergeKrona.html")
#' unlink(c("Number.of.sequences.html", "Number.of.ASVs.html", "mergeKrona.html"))
#' }
#' @return A html file
#' @seealso \code{\link{krona}}
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [Krona](https://github.com/marbl/Krona) if
#'   you use this function.
merge_krona <- function(files = NULL, output = "mergeKrona.html") {
  cmd <-
    paste("ktImportKrona ",
      paste(files, collapse = " "),
      " -o ",
      output,
      sep = ""
    )

  system(command = cmd)
}
################################################################################


#' Interactive Taxonomy plot with Krona from a phyloseq object
#'
#' Construct and run a Krona Chart to compare taxonomic assignations between
#' different conditions.
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
#' plot_krona(GlobalPatterns,"GP-krona", "SampleType")# issues with brackets
#' plot_krona(GlobalPatterns,"GP-krona", "SampleType",trim=T)
#' }
plot_krona<-function(physeq,output,variable, trim=F){
  # Check if KronaTools are installed.
  # if( system(command = "which ktImportText",
  #            intern = FALSE,
  #            ignore.stdout = TRUE)) {
  #   stop("KronaTools are not installed. Please see https://github.com/marbl/Krona/wiki/KronaTools.")
  # }
  if( is.null(tax_table(physeq)) ){
    stop("No taxonomy table available.")
  }
  if( ! variable %in% colnames(sample_data(physeq))){
    stop(paste(variable, "is not a variable in the sample data."))
  }
  if (trim == FALSE) {
    spec.char<- grepl(" |\\(|\\)", as(sample_data(physeq),"data.frame")[,variable] )
    if(sum(spec.char > 0 )){
      message("The following lines contains spaces or brackets.")
      print(paste(which(spec.char)))
      stop("Use trim=TRUE to convert them automatically or convert manually before re-run")
    }
  }
  # Melt the OTU table and merge associated metadata
  df<-psmelt(physeq)
  # Fetch only Abundance, Description and taxonomic rank names columns
  df<-df[ ,c("Abundance", variable, rank_names(physeq)) ]
  # Make sure there are no spaces left
  df[,2]<-gsub(" |\\(|\\)","",df[,2])
  # Convert the field of interest as factor.
  df[,2]<-as.factor(df[,2])

  # Create a directory for krona files
  dir.create(output)

  # For each level of the Description variable
  # Abundance and taxonomic assignations for each OTU are fetched
  # and written to a file that would be processed by Krona.
  for( lvl in levels(df[,2])){
    write.table(
      df[which(df[, 2] == lvl & df[,1] != 0), -2],
      file = paste0(output,"/",lvl, "taxonomy.txt"),
      sep = "\t",row.names = F,col.names = F,na = "",quote = F)
  }
  # Arguments for Krona command
  # taxonomic file and their associated labels.
  krona_args<-paste(output,"/",levels(df[,2]),
                    "taxonomy.txt,",
                    levels(df[,2]),
                    sep = "", collapse = " ")
  # Add html suffix to output
  output<-paste(output,".html",sep = "")
  # Execute Krona command
  # system(paste("ktImportText",
  #              krona_args,
  #              "-o", output,
  #              sep = " "))
  # Run the browser to visualise the output.
  # browseURL(output)
  return(paste(
    "ktImportText",
    krona_args,
    "-o", output,
    sep = " "))
}
