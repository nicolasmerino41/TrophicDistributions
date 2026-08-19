script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
source(file.path(script_dir, "Functions.R"))

# Set to TRUE only when intentionally replacing the cached GloBI snapshot.
REFRESH_GLOBI <- FALSE
PAGE_SIZE <- 10000L
N_WORKERS <- 4L

input_file <- file.path(script_dir, "Data", "thermal_consumers.csv")
tetra_file <- file.path(root_dir, "Data", "TetraEU_pairwise_interactions.csv")
output_dir <- file.path(root_dir, "Outputs", "empirical_degrees")
cache_dir <- file.path(output_dir, "globi_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

thermal <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
consumers <- sort(unique(thermal$consumer))
cat("Thermal consumers:", length(consumers), "\n")

pending <- consumers[REFRESH_GLOBI | !file.exists(file.path(cache_dir, paste0(safe_name(consumers), ".csv")))]
if (length(pending)) {
  cat("Downloading GloBI resources for", length(pending), "consumers.\n")
  if (.Platform$OS.type == "windows" && N_WORKERS > 1L) {
    cl <- parallel::makeCluster(min(N_WORKERS, length(pending)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c("download_one_consumer", "fetch_globi_page", "read_csv_or_empty", "write_csv",
        "safe_name", "cache_dir", "PAGE_SIZE", "REFRESH_GLOBI"),
      envir = environment()
    )
    logs <- parallel::parLapply(cl, pending, function(sp) {
      download_one_consumer(sp, cache_dir, PAGE_SIZE, refresh = REFRESH_GLOBI)
    })
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
  } else {
    logs <- lapply(seq_along(pending), function(i) {
      cat(sprintf("[%d/%d] %s\n", i, length(pending), pending[[i]]))
      download_one_consumer(pending[[i]], cache_dir, PAGE_SIZE, refresh = REFRESH_GLOBI)
    })
  }
  download_log <- do.call(rbind, logs)
} else {
  download_log <- data.frame(query_name = character(), status = character(), rows = integer(),
                             pages = integer(), message = character(), stringsAsFactors = FALSE)
}

if (nrow(download_log) && any(download_log$status == "failed")) {
  download_log$queried_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  write_csv(download_log, file.path(output_dir, "globi_query_log.csv"))
  stop("Some GloBI queries failed. Re-run the script; completed consumers are cached.")
}

query_log <- do.call(rbind, lapply(consumers, function(sp) {
  cache_file <- file.path(cache_dir, paste0(safe_name(sp), ".csv"))
  cached <- read_csv_or_empty(cache_file)
  data.frame(query_name = sp, status = if (file.exists(cache_file)) "complete" else "missing",
             rows = nrow(cached), pages = NA_integer_, message = "", stringsAsFactors = FALSE)
}))
query_log$queried_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
write_csv(query_log, file.path(output_dir, "globi_query_log.csv"))

cat("Combining and deduplicating GloBI resource records.\n")
globi <- summarize_globi(consumers, cache_dir)
write_csv(globi$links, file.path(output_dir, "globi_resource_links.csv"))
write_csv(globi$degree, file.path(output_dir, "globi_consumer_degrees.csv"))

cat("Processing the shared TetraEU metaweb.\n")
tetra <- process_tetraeu(tetra_file)
write_csv(tetra$links, file.path(output_dir, "tetraeu_resource_links.csv"))
write_csv(tetra$degree, file.path(output_dir, "tetraeu_consumer_degrees.csv"))

globi_degree <- merge(thermal, globi$degree, by = "consumer", all.x = TRUE, sort = TRUE)
globi_degree$dataset <- "GloBI consumers"
tetra_degree <- tetra$degree
tetra_degree$dataset <- "TetraEU consumers"
tetra_degree$metrics <- ""
tetra_degree$n_metrics <- NA_integer_

common <- c("dataset", "consumer", "degree_species", "degree_all_taxa", "self_links", "metrics", "n_metrics")
combined <- rbind(globi_degree[, common], tetra_degree[, common])
names(combined)[names(combined) == "degree_species"] <- "degree"
write_csv(combined, file.path(output_dir, "consumer_degrees_combined.csv"))

manifest <- data.frame(
  item = c("GloBI endpoint", "GloBI interaction umbrella", "GloBI extraction time (UTC)",
           "GloBI consumers retained in thermal analysis",
           "GloBI consumers with positive species-level degree",
           "TetraEU consumers", "TetraEU consumers with positive species-level degree",
           "TetraEU source file", "Primary degree definition"),
  value = c("https://api.globalbioticinteractions.org/interaction.csv", "eats",
            format(Sys.time(), tz = "UTC", usetz = TRUE), length(consumers),
            sum(is.finite(globi_degree$degree_species) & globi_degree$degree_species > 0),
            nrow(tetra_degree),
            sum(is.finite(tetra_degree$degree_species) & tetra_degree$degree_species > 0),
            "Data/TetraEU_pairwise_interactions.csv",
            "Unique species-level non-self resources per consumer"),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(output_dir, "extraction_manifest.csv"))

cat("Done. Outputs:", output_dir, "\n")
