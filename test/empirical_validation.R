#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(root_dir, "ConsumerValidation.R"))

stopifnot(
  consumer_taxon_status(kingdom = "Animalia") == "animal",
  consumer_taxon_status(kingdom = "Metazoa") == "animal",
  consumer_taxon_status(phylum = "Arthropoda") == "animal",
  consumer_taxon_status(kingdom = "Plantae", phylum = "Streptophyta") == "non_animal",
  consumer_taxon_status(phylum = "Chlorophyta") == "non_animal",
  consumer_taxon_status() == "unresolved"
)

master <- read.csv(
  file.path(root_dir, "Data", "thermal_interactions_master.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
degrees <- read.csv(
  file.path(root_dir, "Outputs", "empirical_degrees", "consumer_degrees_combined.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)

stopifnot(
  "consumer_taxon_status" %in% names(master),
  all(master$consumer_taxon_status == "animal"),
  all(degrees$consumer_taxon_status == "animal"),
  !any(master$consumer == "Baccharis pilularis"),
  any(master$resource == "Baccharis pilularis")
)

cat("Empirical consumer validation passed.\n")
