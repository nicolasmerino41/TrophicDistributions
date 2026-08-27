script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else {
  file.path(getwd(), "test")
}
repo_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

non_animal_pattern <- paste0(
  "(^|[^a-z])(",
  paste(c("plantae", "viridiplantae", "streptophyta", "tracheophyta",
          "chlorophyta", "rhodophyta", "ochrophyta", "bacillariophyta",
          "phaeophyceae", "fungi", "ascomycota", "basidiomycota",
          "bacteria", "archaea"), collapse = "|"),
  ")([^a-z]|$)"
)

master <- read.csv(file.path(repo_root, "Data", "thermal_species_master.csv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
interactions <- read.csv(file.path(repo_root, "Data", "thermal_interactions_master.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
degree_links <- read.csv(file.path(repo_root, "Outputs", "empirical_degrees",
                                   "globi_resource_links.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(
  "is_animal" %in% names(master),
  all(master$is_animal[match(interactions$consumer, master$accepted_species)]),
  all(master$is_animal[match(interactions$resource, master$accepted_species)]),
  !any(grepl(non_animal_pattern, tolower(degree_links$source_taxon_path))),
  !any(grepl(non_animal_pattern, tolower(degree_links$target_taxon_path)))
)

cat("Animal-only empirical validation passed.\n")
