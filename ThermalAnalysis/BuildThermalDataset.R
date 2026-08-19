#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir <- file.path(root_dir, "Data")
output_dir <- file.path(root_dir, "Outputs", "thermal_metrics")
cache_dir <- file.path(output_dir, "globi_cache")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

TRAITS <- c("ctmax", "lt50", "ctmin", "ltmax")
SOURCE_PRIORITY <- c(ThermoFresh = 1L, ComteOlden = 2L, GlobTherm = 3L)
GLOBI_PAGE_SIZE <- 5000L
GLOBI_WORKERS <- max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L))
REFRESH_GLOBI <- identical(tolower(Sys.getenv("TD_REFRESH_GLOBI", "false")), "true")

canonical_binomial <- function(x) {
  x <- as.character(x)
  clean <- gsub("[_(),\\[\\]]", " ", x)
  clean <- gsub("\\s+", " ", trimws(clean))
  pieces <- strsplit(clean, " ", fixed = TRUE)
  vapply(pieces, function(part) {
    if (length(part) < 2L) return(NA_character_)
    genus <- gsub("[^A-Za-z-]", "", part[[1L]])
    species <- gsub("[^A-Za-z-]", "", part[[2L]])
    if (!nzchar(genus) || !nzchar(species) || tolower(species) %in% c("sp", "spp")) {
      return(NA_character_)
    }
    paste0(toupper(substr(genus, 1L, 1L)),
           tolower(substr(genus, 2L, nchar(genus))), " ", tolower(species))
  }, character(1L), USE.NAMES = FALSE)
}

numeric_value <- function(x) {
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
}

normalize_metric <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y <- gsub("[ _-]", "", y)
  out <- rep(NA_character_, length(y))
  out[grepl("ctmax|criticalthermalmaximum", y)] <- "ctmax"
  out[grepl("ctmin|criticalthermalminimum", y)] <- "ctmin"
  out[grepl("lt50", y)] <- "lt50"
  out[grepl("ltmax", y) | (grepl("upperlethal", y) & grepl("temp", y))] <- "ltmax"
  out
}

first_nonempty <- function(x) {
  x <- unique(trimws(as.character(x[!is.na(x)])))
  x <- x[nzchar(x)]
  if (length(x)) x[[1L]] else NA_character_
}

collapse_values <- function(x) {
  x <- sort(unique(trimws(as.character(x[!is.na(x)]))))
  x <- x[nzchar(x)]
  if (length(x)) paste(x, collapse = ";") else ""
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

habitat_from_thermofresh <- function(dt) {
  freshwater <- truthy(dt$freshwater) | grepl("freshwater", tolower(dt$habitat), fixed = TRUE)
  marine <- truthy(dt$marine) | grepl("marine", tolower(dt$habitat), fixed = TRUE)
  brackish <- truthy(dt$brackish) | grepl("brackish", tolower(dt$habitat), fixed = TRUE)
  terrestrial <- grepl("terrestrial", tolower(dt$habitat), fixed = TRUE)
  nrealm <- freshwater + marine + terrestrial
  fifelse(nrealm > 1L | brackish, "mixed",
          fifelse(freshwater, "freshwater",
                  fifelse(marine, "marine",
                          fifelse(terrestrial, "terrestrial", "unknown"))))
}

read_thermofresh <- function(path) {
  raw <- fread(path, encoding = "Latin-1", na.strings = c("", "NA"),
               fill = TRUE, showProgress = FALSE)
  raw[, species_raw := as.character(species)]
  raw[, species_alias := canonical_binomial(species_raw)]
  raw[, metric_norm := normalize_metric(metric)]
  raw[, trait_value := numeric_value(tol)]
  raw[, habitat_class := habitat_from_thermofresh(raw)]
  out <- raw[!is.na(species_alias) & metric_norm %chin% TRAITS & is.finite(trait_value), .(
    species_alias, species_raw, metric = metric_norm, trait_value,
    trait_source = "ThermoFresh", source_record_id = as.character(test_id),
    kingdom = as.character(kingdom), phylum = as.character(phylum),
    class = as.character(class), order = as.character(order), family = as.character(family),
    habitat_class, group = as.character(group)
  )]
  unique(out)
}

read_globtherm <- function(path) {
  # The archived GlobTherm CSV contains one irregular quoted row that causes
  # fread to stop early. Base R's parser reads all 2,133 records correctly.
  raw <- as.data.table(read.csv(
    path, encoding = "latin1", stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA"), fill = TRUE
  ))
  raw[, species_raw := trimws(paste(Genus, Species))]
  raw[, species_alias := canonical_binomial(species_raw)]
  pairs <- list(
    c("Tmax", "max_metric", "upper_primary"),
    c("Tmax_2", "max_metric_2", "upper_secondary"),
    c("tmin", "min_metric", "lower_primary"),
    c("tmin_2", "min_metric_2", "lower_secondary")
  )
  parts <- lapply(pairs, function(p) {
    dt <- raw[, .(
      species_alias, species_raw,
      metric = normalize_metric(get(p[[2L]])),
      trait_value = numeric_value(get(p[[1L]])),
      trait_source = "GlobTherm",
      source_record_id = paste0("GlobTherm_row_", .I, "_", p[[3L]]),
      kingdom = NA_character_, phylum = as.character(Phylum),
      class = as.character(Class), order = as.character(Order), family = as.character(Family),
      habitat_class = "unknown", group = NA_character_
    )]
    dt[!is.na(species_alias) & metric %chin% TRAITS & is.finite(trait_value)]
  })
  unique(rbindlist(parts, fill = TRUE))
}

read_comte_olden <- function(path) {
  raw <- fread(path, sep = ";", dec = ",", encoding = "Latin-1",
               na.strings = c("", "NA"), showProgress = FALSE)
  raw[, species_raw := as.character(Species)]
  raw[, species_alias := canonical_binomial(species_raw)]
  raw[, trait_value := numeric_value(CTmax)]
  raw[, habitat_class := fifelse(tolower(trimws(`Realm affinity`)) == "freshwater", "freshwater",
                                  fifelse(tolower(trimws(`Realm affinity`)) == "marine", "marine", "unknown"))]
  out <- raw[!is.na(species_alias) & is.finite(trait_value), .(
    species_alias, species_raw, metric = "ctmax", trait_value,
    trait_source = "ComteOlden", source_record_id = paste0("ComteOlden_row_", .I),
    kingdom = "Animalia", phylum = "Chordata", class = "Actinopterygii",
    order = NA_character_, family = as.character(Family), habitat_class,
    group = "fish"
  )]
  unique(out)
}

cat("Reading and harmonizing thermal traits.\n")
trait_records <- rbindlist(list(
  read_thermofresh(file.path(data_dir, "thermtol_comb_final.csv")),
  read_globtherm(file.path(data_dir, "GlobalTherm_upload_02_11_17.csv")),
  read_comte_olden(file.path(data_dir, "Comte_Olden_Data_Imputed.csv"))
), fill = TRUE)

trait_records[, source_priority := unname(SOURCE_PRIORITY[trait_source])]
setorder(trait_records, species_alias, metric, source_priority, trait_source)

source_traits <- trait_records[, .(
  source_trait_mean = mean(trait_value),
  source_trait_sd = if (.N > 1L) sd(trait_value) else NA_real_,
  source_trait_n = .N
), by = .(species_alias, metric, trait_source, source_priority)]

selected_traits <- source_traits[order(source_priority)][, .SD[1L], by = .(species_alias, metric)]
setnames(selected_traits,
         c("source_trait_mean", "source_trait_sd", "source_trait_n", "trait_source"),
         c("trait_value", "trait_sd", "trait_n", "selected_trait_source"))

local_species <- trait_records[, .(
  species_raw_names = collapse_values(species_raw),
  kingdom_local = first_nonempty(kingdom), phylum_local = first_nonempty(phylum),
  class_local = first_nonempty(class), order_local = first_nonempty(order),
  family_local = first_nonempty(family),
  habitat_sources = collapse_values(habitat_class[habitat_class != "unknown"]),
  group_local = first_nonempty(group),
  trait_sources = collapse_values(trait_source)
), by = species_alias]

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

read_csv_or_empty <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) return(data.table())
  tryCatch(fread(path, na.strings = c("", "NA"), showProgress = FALSE),
           error = function(e) data.table())
}

fetch_globi_page <- function(query_name, skip, limit, destination) {
  fields <- c(
    "source_taxon_external_id", "source_taxon_name", "source_taxon_path",
    "interaction_type", "target_taxon_external_id", "target_taxon_name",
    "target_taxon_path", "study_title", "study_citation", "study_external_id",
    "study_source_citation"
  )
  params <- c(
    paste0("sourceTaxon=", utils::URLencode(query_name, reserved = TRUE)),
    "interactionType=eats", "includeObservations=false",
    paste0("limit=", limit), paste0("skip=", skip), paste0("field=", fields)
  )
  url <- paste0("https://api.globalbioticinteractions.org/interaction.csv?",
                paste(params, collapse = "&"))
  old_timeout <- getOption("timeout")
  options(timeout = max(180L, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  utils::download.file(url, destination, mode = "wb", quiet = TRUE, method = "libcurl",
                       headers = c("User-Agent" = "TrophicDistributions thermal-data audit"))
  read_csv_or_empty(destination)
}

download_one_consumer <- function(query_name, cache_dir, page_size = 5000L,
                                  max_retries = 5L, refresh = FALSE) {
  cache_file <- file.path(cache_dir, paste0(safe_name(query_name), ".csv"))
  if (!refresh && file.exists(cache_file)) {
    cached <- read_csv_or_empty(cache_file)
    return(data.table(query_name, status = "cached", rows = nrow(cached), message = ""))
  }
  pages <- list()
  skip <- 0L
  status <- "ok"
  note <- ""
  repeat {
    page <- NULL
    last_error <- NULL
    for (attempt in seq_len(max_retries)) {
      page_file <- tempfile(fileext = ".csv")
      page <- tryCatch(fetch_globi_page(query_name, skip, page_size, page_file),
                       error = function(e) { last_error <<- conditionMessage(e); NULL })
      unlink(page_file)
      if (!is.null(page)) break
      Sys.sleep(min(2^(attempt - 1L), 16L))
    }
    if (is.null(page)) {
      status <- "failed"
      note <- if (is.null(last_error)) "Unknown download error" else last_error
      break
    }
    if (!nrow(page)) break
    pages[[length(pages) + 1L]] <- page
    if (nrow(page) < page_size) break
    skip <- skip + page_size
  }
  if (status == "ok") {
    result <- if (length(pages)) rbindlist(pages, fill = TRUE) else data.table()
    if (!nrow(result)) result <- data.table(query_name = character())
    result[, query_name := query_name]
    setcolorder(result, c("query_name", setdiff(names(result), "query_name")))
    fwrite(result, cache_file, na = "")
  }
  data.table(query_name, status, rows = if (status == "ok") nrow(result) else 0L,
             message = note)
}

query_names <- sort(unique(trait_records$species_alias))
pending <- query_names[REFRESH_GLOBI | !file.exists(file.path(cache_dir, paste0(safe_name(query_names), ".csv")))]
if (length(pending)) {
  cat(sprintf("Querying GloBI for %d thermal-trait species using %d workers.\n",
              length(pending), GLOBI_WORKERS))
  if (.Platform$OS.type == "windows" && GLOBI_WORKERS > 1L) {
    cluster <- parallel::makeCluster(min(GLOBI_WORKERS, length(pending)))
    parallel::clusterEvalQ(cluster, library(data.table))
    parallel::clusterExport(cluster, c(
      "download_one_consumer", "fetch_globi_page", "read_csv_or_empty", "safe_name",
      "cache_dir", "GLOBI_PAGE_SIZE", "REFRESH_GLOBI"
    ), envir = environment())
    logs <- parallel::parLapply(cluster, pending, function(sp) {
      download_one_consumer(sp, cache_dir, GLOBI_PAGE_SIZE, refresh = REFRESH_GLOBI)
    })
    parallel::stopCluster(cluster)
  } else {
    logs <- lapply(pending, download_one_consumer, cache_dir = cache_dir,
                   page_size = GLOBI_PAGE_SIZE, refresh = REFRESH_GLOBI)
  }
  query_log_new <- rbindlist(logs, fill = TRUE)
  if (any(query_log_new$status == "failed")) {
    fwrite(query_log_new, file.path(output_dir, "globi_query_log.csv"), na = "")
    stop("Some GloBI queries failed. Re-run the script; completed queries are cached.")
  }
}

cat("Combining GloBI records and constructing the taxonomic alias map.\n")
globi_parts <- lapply(query_names, function(sp) {
  path <- file.path(cache_dir, paste0(safe_name(sp), ".csv"))
  dt <- read_csv_or_empty(path)
  if (!nrow(dt)) return(NULL)
  dt[, query_name := sp]
  dt
})
globi_raw <- rbindlist(Filter(Negate(is.null), globi_parts), fill = TRUE)

for (nm in c("query_name", "source_taxon_external_id", "source_taxon_name",
             "source_taxon_path", "interaction_type", "target_taxon_external_id",
             "target_taxon_name", "target_taxon_path", "study_title", "study_citation",
             "study_external_id", "study_source_citation")) {
  if (!nm %in% names(globi_raw)) globi_raw[, (nm) := ""]
}

last_path_name <- function(path, fallback) {
  path <- as.character(path)
  fallback <- as.character(fallback)
  vapply(seq_along(path), function(i) {
    if (!is.na(path[[i]]) && nzchar(trimws(path[[i]]))) {
      part <- strsplit(path[[i]], "\\s*\\|\\s*")[[1L]]
      return(tail(part, 1L))
    }
    fallback[[i]]
  }, character(1L))
}

if (nrow(globi_raw)) {
  globi_raw[, query_alias := canonical_binomial(query_name)]
  globi_raw[, source_species_globi := canonical_binomial(last_path_name(source_taxon_path, source_taxon_name))]
  globi_raw[, resource_species_globi := canonical_binomial(last_path_name(target_taxon_path, target_taxon_name))]
}

alias_map <- data.table(species_alias = query_names, accepted_species = query_names)
if (nrow(globi_raw)) {
  accepted <- globi_raw[!is.na(query_alias) & !is.na(source_species_globi),
                        .(accepted_species = names(sort(table(source_species_globi), decreasing = TRUE))[1L]),
                        by = .(species_alias = query_alias)]
  alias_map[accepted, on = "species_alias", accepted_species := i.accepted_species]
}

# Collapse aliases that GloBI maps to the same accepted binomial.
trait_records <- alias_map[trait_records, on = "species_alias"]
trait_records[is.na(accepted_species), accepted_species := species_alias]
source_traits <- trait_records[, .(
  source_trait_mean = mean(trait_value),
  source_trait_sd = if (.N > 1L) sd(trait_value) else NA_real_,
  source_trait_n = .N,
  source_priority = min(source_priority)
), by = .(accepted_species, metric, trait_source)]
source_traits[, source_priority := unname(SOURCE_PRIORITY[trait_source])]
selected_traits <- source_traits[order(source_priority)][, .SD[1L], by = .(accepted_species, metric)]
setnames(selected_traits,
         c("source_trait_mean", "source_trait_sd", "source_trait_n", "trait_source"),
         c("trait_value", "trait_sd", "trait_n", "selected_trait_source"))
source_traits_wide <- dcast(
  source_traits, accepted_species ~ trait_source + metric,
  value.var = c("source_trait_mean", "source_trait_sd", "source_trait_n")
)

local_species <- alias_map[local_species, on = "species_alias"]
local_species[is.na(accepted_species), accepted_species := species_alias]
species_local <- local_species[, .(
  aliases = collapse_values(species_alias), species_raw_names = collapse_values(species_raw_names),
  kingdom_local = first_nonempty(kingdom_local), phylum_local = first_nonempty(phylum_local),
  class_local = first_nonempty(class_local), order_local = first_nonempty(order_local),
  family_local = first_nonempty(family_local), habitat_sources = collapse_values(habitat_sources),
  group_local = first_nonempty(group_local), trait_sources = collapse_values(trait_sources)
), by = accepted_species]

alias_lookup <- unique(rbindlist(list(
  alias_map[, .(alias = species_alias, accepted_species)],
  data.table(alias = unique(alias_map$accepted_species), accepted_species = unique(alias_map$accepted_species))
)))
setkey(alias_lookup, alias)

map_species <- function(x) {
  x <- canonical_binomial(x)
  mapped <- alias_lookup[data.table(alias = x), on = "alias", accepted_species]
  mapped[is.na(mapped)] <- x[is.na(mapped)]
  mapped
}

if (nrow(globi_raw)) {
  globi_raw[, consumer := map_species(query_alias)]
  globi_raw[, resource := map_species(resource_species_globi)]
  globi_links_all <- globi_raw[!is.na(consumer) & !is.na(resource), .(
    interaction_types = collapse_values(interaction_type),
    study_titles = collapse_values(study_title),
    study_citations = collapse_values(c(study_citation, study_source_citation)),
    study_ids = collapse_values(study_external_id),
    reported_consumer_names = collapse_values(source_taxon_name),
    reported_resource_names = collapse_values(target_taxon_name)
  ), by = .(consumer, resource)]
  globi_links_all[, GloBI := TRUE]
} else {
  globi_links_all <- data.table(consumer = character(), resource = character(),
                                interaction_types = character(), study_titles = character(),
                                study_citations = character(), study_ids = character(),
                                reported_consumer_names = character(), reported_resource_names = character(),
                                GloBI = logical())
}

globi_taxonomy <- if (nrow(globi_raw)) {
  source_tax <- globi_raw[!is.na(consumer), .(
    accepted_species = consumer,
    globi_taxon_paths = collapse_values(source_taxon_path),
    globi_taxon_ids = collapse_values(source_taxon_external_id)
  ), by = consumer][, consumer := NULL]
  resource_tax <- globi_raw[!is.na(resource), .(
    accepted_species = resource,
    globi_taxon_paths = collapse_values(target_taxon_path),
    globi_taxon_ids = collapse_values(target_taxon_external_id)
  ), by = resource][, resource := NULL]
  rbindlist(list(source_tax, resource_tax), fill = TRUE)[, .(
    globi_taxon_paths = collapse_values(globi_taxon_paths),
    globi_taxon_ids = collapse_values(globi_taxon_ids)
  ), by = accepted_species]
} else {
  data.table(accepted_species = character(), globi_taxon_paths = character(),
             globi_taxon_ids = character())
}

cat("Processing TetraEU interactions.\n")
tetra_raw <- fread(file.path(data_dir, "TetraEU_pairwise_interactions.csv"), sep = ";",
                   encoding = "Latin-1", na.strings = c("", "NA"), showProgress = FALSE)
tetra_raw[, consumer_raw := canonical_binomial(sourceTaxonName)]
tetra_raw[, resource_raw := canonical_binomial(targetTaxonName)]
tetra_raw[, consumer := map_species(consumer_raw)]
tetra_raw[, resource := map_species(resource_raw)]
tetra_links_all <- tetra_raw[!is.na(consumer) & !is.na(resource), .(
  interaction_types = collapse_values(interactionTypeName),
  study_titles = "",
  study_citations = collapse_values(c(reference1Id, reference2Id, reference3Id, reference4Id)),
  study_ids = collapse_values(c(reference1Id, reference2Id, reference3Id, reference4Id)),
  reported_consumer_names = collapse_values(sourceTaxonName),
  reported_resource_names = collapse_values(targetTaxonName)
), by = .(consumer, resource)]
tetra_links_all[, TetraEU := TRUE]
tetra_taxonomy <- rbindlist(list(
  tetra_raw[!is.na(consumer), .(accepted_species = consumer,
                                tetraeu_taxon_ids = collapse_values(sourceTaxonId)), by = consumer][, consumer := NULL],
  tetra_raw[!is.na(resource), .(accepted_species = resource,
                                tetraeu_taxon_ids = collapse_values(targetTaxonId)), by = resource][, resource := NULL]
), fill = TRUE)[, .(tetraeu_taxon_ids = collapse_values(tetraeu_taxon_ids)), by = accepted_species]

recorded <- merge(globi_links_all, tetra_links_all, by = c("consumer", "resource"), all = TRUE,
                  suffixes = c("_globi", "_tetra"))
for (flag in c("GloBI", "TetraEU")) recorded[is.na(get(flag)), (flag) := FALSE]
recorded[, interaction_sources := fifelse(GloBI & TetraEU, "GloBI;TetraEU",
                                           fifelse(GloBI, "GloBI", "TetraEU"))]
for (stem in c("interaction_types", "study_titles", "study_citations", "study_ids",
               "reported_consumer_names", "reported_resource_names")) {
  left <- paste0(stem, "_globi")
  right <- paste0(stem, "_tetra")
  recorded[, (stem) := mapply(function(a, b) collapse_values(c(a, b)), get(left), get(right))]
}
recorded[, self_link := consumer == resource]
drop_cols <- grep("_(globi|tetra)$", names(recorded), value = TRUE)
recorded[, (drop_cols) := NULL]
setorder(recorded, consumer, resource)

trait_wide <- dcast(selected_traits, accepted_species ~ metric,
                    value.var = c("trait_value", "trait_sd", "trait_n", "selected_trait_source"))
setnames(trait_wide, "accepted_species", "species")

consumer_traits <- copy(trait_wide)
setnames(consumer_traits, names(consumer_traits),
         c("consumer", paste0("consumer_", names(consumer_traits)[-1L])))
resource_traits <- copy(trait_wide)
setnames(resource_traits, names(resource_traits),
         c("resource", paste0("resource_", names(resource_traits)[-1L])))

master <- merge(recorded, consumer_traits, by = "consumer", all.x = TRUE)
master <- merge(master, resource_traits, by = "resource", all.x = TRUE)
for (metric in TRAITS) {
  cv <- paste0("consumer_trait_value_", metric)
  rv <- paste0("resource_trait_value_", metric)
  master[, (paste0("eligible_", metric)) := is.finite(get(cv)) & is.finite(get(rv)) & !self_link]
}
master[, any_metric_eligible := Reduce(`|`, .SD), .SDcols = paste0("eligible_", TRAITS)]
analysis_master <- master[any_metric_eligible == TRUE]
setorder(analysis_master, consumer, resource)

coverage_rows <- list()
for (metric in TRAITS) {
  value_col <- paste0("resource_trait_value_", metric)
  consumer_value_col <- paste0("consumer_trait_value_", metric)
  eligible_consumers <- unique(master[is.finite(get(consumer_value_col)), consumer])
  analysis_consumers <- unique(analysis_master[get(paste0("eligible_", metric)) == TRUE, consumer])
  coverage_rows[[metric]] <- master[consumer %chin% eligible_consumers & !self_link, .(
    recorded_prey_species = uniqueN(resource),
    prey_with_trait = uniqueN(resource[is.finite(get(value_col))]),
    prey_trait_coverage = uniqueN(resource[is.finite(get(value_col))]) / uniqueN(resource)
  ), by = consumer][, `:=`(
    metric = metric,
    included_in_analysis = consumer %chin% analysis_consumers
  )]
}
coverage <- rbindlist(coverage_rows, fill = TRUE)
setcolorder(coverage, c("consumer", "metric", "recorded_prey_species", "prey_with_trait",
                       "prey_trait_coverage", "included_in_analysis"))

interaction_species <- data.table(accepted_species = sort(unique(c(recorded$consumer, recorded$resource))))
species_master <- merge(interaction_species, species_local, by = "accepted_species", all = TRUE)
species_master <- merge(species_master, trait_wide,
                        by.x = "accepted_species", by.y = "species", all = TRUE)
species_master <- merge(species_master, source_traits_wide, by = "accepted_species", all = TRUE)
species_master <- merge(species_master, globi_taxonomy, by = "accepted_species", all = TRUE)
species_master <- merge(species_master, tetra_taxonomy, by = "accepted_species", all = TRUE)
species_master[, habitat_class := fifelse(grepl(";", habitat_sources), "mixed",
                                          fifelse(is.na(habitat_sources) | !nzchar(habitat_sources),
                                                  "unknown", habitat_sources))]
species_master[, recorded_as_consumer := accepted_species %chin% recorded$consumer]
species_master[, recorded_as_resource := accepted_species %chin% recorded$resource]
species_master[, available_metrics := vapply(seq_len(.N), function(i) {
  present <- TRAITS[vapply(TRAITS, function(metric) {
    col <- paste0("trait_value_", metric)
    col %in% names(species_master) && is.finite(species_master[[col]][[i]])
  }, logical(1L))]
  paste(present, collapse = ";")
}, character(1L))]
species_master[, n_available_metrics := fifelse(nzchar(available_metrics),
                                                 lengths(strsplit(available_metrics, ";", fixed = TRUE)), 0L)]
species_master[, recorded_consumer_degree := vapply(accepted_species, function(sp) {
  uniqueN(recorded[consumer == sp & !self_link, resource])
}, integer(1L))]
species_master[, recorded_resource_consumers := vapply(accepted_species, function(sp) {
  uniqueN(recorded[resource == sp & !self_link, consumer])
}, integer(1L))]
setorder(species_master, accepted_species)

query_log <- rbindlist(lapply(query_names, function(sp) {
  cache_file <- file.path(cache_dir, paste0(safe_name(sp), ".csv"))
  data.table(query_name = sp, status = if (file.exists(cache_file)) "complete" else "missing",
             rows = nrow(read_csv_or_empty(cache_file)))
}))
query_log[, queried_at_utc := format(Sys.time(), tz = "UTC", usetz = TRUE)]

cat("Validating canonical datasets.\n")
stopifnot(
  !anyDuplicated(species_master$accepted_species),
  !anyDuplicated(recorded[, .(consumer, resource)]),
  !anyDuplicated(analysis_master[, .(consumer, resource)]),
  !any(analysis_master$self_link),
  all(is.finite(coverage$prey_trait_coverage)),
  all(coverage$prey_trait_coverage >= 0 & coverage$prey_trait_coverage <= 1),
  all(analysis_master$consumer %chin% species_master$accepted_species),
  all(analysis_master$resource %chin% species_master$accepted_species)
)
for (metric in TRAITS) {
  eligible <- analysis_master[[paste0("eligible_", metric)]]
  expected <- is.finite(analysis_master[[paste0("consumer_trait_value_", metric)]]) &
    is.finite(analysis_master[[paste0("resource_trait_value_", metric)]]) &
    !analysis_master$self_link
  stopifnot(identical(eligible, expected))
}

cat("Writing canonical datasets.\n")
fwrite(trait_records, file.path(data_dir, "thermal_traits_long.csv"), na = "")
fwrite(species_master, file.path(data_dir, "thermal_species_master.csv"), na = "")
fwrite(recorded, file.path(data_dir, "thermal_recorded_interactions.csv"), na = "")
fwrite(analysis_master, file.path(data_dir, "thermal_interactions_master.csv"), na = "")
fwrite(coverage, file.path(data_dir, "thermal_consumer_coverage.csv"), na = "")
fwrite(query_log, file.path(output_dir, "globi_query_log.csv"), na = "")

empirical_degree_dir <- file.path(root_dir, "EmpiricalDegrees", "Data")
if (dir.exists(dirname(empirical_degree_dir))) {
  dir.create(empirical_degree_dir, recursive = TRUE, showWarnings = FALSE)
  thermal_consumers <- unique(analysis_master[, .(consumer)])
  thermal_consumers[, metrics := vapply(consumer, function(sp) {
    paste(TRAITS[vapply(TRAITS, function(metric) {
      any(analysis_master$consumer == sp & analysis_master[[paste0("eligible_", metric)]])
    }, logical(1L))], collapse = ";")
  }, character(1L))]
  thermal_consumers[, n_metrics := fifelse(nzchar(metrics),
                                            lengths(strsplit(metrics, ";", fixed = TRUE)), 0L)]
  setorder(thermal_consumers, consumer)
  fwrite(thermal_consumers, file.path(empirical_degree_dir, "thermal_consumers.csv"), na = "")
}

manifest <- data.table(
  item = c(
    "Build time UTC", "GloBI endpoint", "GloBI interaction umbrella",
    "GloBI queried species", "TetraEU source", "Trait sources",
    "Trait selection priority", "Species in master", "Species with usable traits",
    "Recorded consumer-resource pairs", "Analysis-ready pairs"
  ),
  value = c(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    "https://api.globalbioticinteractions.org/interaction.csv", "eats",
    length(query_names), "Data/TetraEU_pairwise_interactions.csv",
    "ThermoFresh;GlobTherm;ComteOlden",
    "ThermoFresh > ComteOlden > GlobTherm",
    nrow(species_master), sum(species_master$n_available_metrics > 0),
    nrow(recorded), nrow(analysis_master)
  )
)
fwrite(manifest, file.path(data_dir, "thermal_data_manifest.csv"), na = "")

cat(sprintf("Done: %d species, %d recorded pairs, %d analysis-ready pairs.\n",
            nrow(species_master), nrow(recorded), nrow(analysis_master)))
