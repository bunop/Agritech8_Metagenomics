source("renv/activate.R")

# Load necessary libraries only if they are installed
safe_load <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    library(package, character.only = TRUE)
  }
}

safe_load("here")
safe_load("gitcreds")
safe_load("targets")
safe_load("tarchetypes")
