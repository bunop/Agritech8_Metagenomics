#!/usr/bin/env Rscript

library(here)
library(optparse)
library(phyloseq)
library(tidyverse)

source(here::here("R/load_stuff.R"))
source(here::here("R/krona.R"))

# Define command-line options
option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    type = "character",
    help = "phyloseq RDS input file [required]",
    metavar = "input"
  ),
  optparse::make_option(
    c("-m", "--metadata"),
    type = "character",
    help = "metadata input file [required]",
    metavar = "input"
  ),
  optparse::make_option(
    c("-o", "--output"),
    type = "character",
    help = "output file [required]",
    metavar = "output"
  ),
  optparse::make_option(
    c("-v", "--variable"),
    type = "character",
    help = "variable to use for the krona plot [default: %default]",
    metavar = "variable",
    default = "sample_name"
  )
)

# Parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opts <- optparse::parse_args(opt_parser)

# Check for required arguments
if (is.null(opts$input) || is.null(opts$output) || is.null(opts$metadata)) {
  print_help(opt_parser)
  stop("--input, --metadata and --output parameters are required.", call. = FALSE)
}

# Load input file
metadata <- load_metadata(opts$metadata)
phyloseq_object <- load_phyloseq(metadata, opts$input)

# Generate the Krona command
cmd <- get_krona_cmd(
  phyloseq_object,
  output = opts$output,
  variable = opts$variable
)

print(cmd)
