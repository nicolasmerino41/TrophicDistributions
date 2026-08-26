consumer_taxon_status <- function(kingdom = "", phylum = "", class = "", group = "", taxon_path = "") {
  fields <- list(kingdom, phylum, class, group, taxon_path)
  n <- max(1L, lengths(fields))
  fields <- lapply(fields, function(x) {
    x <- rep_len(as.character(x), n)
    x[is.na(x)] <- ""
    tolower(trimws(x))
  })
  text <- do.call(paste, c(fields, sep = " | "))

  animal <- grepl("(^|[^a-z])(animalia|metazoa)([^a-z]|$)", text) |
    grepl("(^|[^a-z])(chordata|arthropoda|mollusca|annelida|cnidaria|echinodermata|platyhelminthes|nematoda|rotifera|porifera|bryozoa|brachiopoda|tardigrada|onychophora)([^a-z]|$)", text)
  non_animal <- grepl("(^|[^a-z])(plantae|viridiplantae|streptophyta|chlorophyta|rhodophyta|ochrophyta|bacillariophyta|fungi|ascomycota|basidiomycota|bacteria|archaea)([^a-z]|$)", text)

  ifelse(animal & !non_animal, "animal",
         ifelse(non_animal, "non_animal", "unresolved"))
}

require_verified_animal_consumers <- function(status, consumer, context) {
  bad <- unique(as.character(consumer[is.na(status) | status != "animal"]))
  bad <- bad[!is.na(bad) & nzchar(bad)]
  if (length(bad)) {
    stop(sprintf(
      "%s contains consumers that are not verified animals: %s",
      context, paste(head(sort(bad), 10L), collapse = ", ")
    ))
  }
  invisible(TRUE)
}
