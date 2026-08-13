canonical_binomial <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  clean <- gsub("[_(),\\[\\]]", " ", x)
  clean <- gsub("\\s+", " ", trimws(clean))
  valid <- !is.na(clean) & grepl("^[A-Za-z-]+ [A-Za-z-]+", clean)
  genus <- sub("^([A-Za-z-]+).*$", "\\1", clean[valid])
  species <- sub("^[A-Za-z-]+ ([A-Za-z-]+).*$", "\\1", clean[valid])
  valid_species <- !tolower(species) %in% c("sp", "spp")
  idx <- which(valid)[valid_species]
  genus <- genus[valid_species]
  species <- species[valid_species]
  out[idx] <- paste0(toupper(substr(genus, 1L, 1L)),
                     tolower(substr(genus, 2L, nchar(genus))), " ", tolower(species))
  out
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

species_from_globi <- function(name, path, external_id) {
  name <- as.character(name)
  path <- as.character(path)
  external_id <- as.character(external_id)
  out <- rep(NA_character_, length(name))
  for (i in seq_along(name)) {
    if (is.na(name[[i]]) || !nzchar(trimws(name[[i]]))) next
    path_last <- ""
    if (!is.na(path[[i]]) && nzchar(trimws(path[[i]]))) {
      path_parts <- strsplit(path[[i]], "\\s*\\|\\s*")[[1]]
      path_last <- tail(path_parts, 1L)
    }
    candidate <- canonical_binomial(if (nzchar(path_last)) path_last else name[[i]])
    if (is.na(candidate)) next
    name_candidate <- canonical_binomial(name[[i]])
    path_ok <- nzchar(path_last) && !is.na(name_candidate) && identical(candidate, name_candidate)
    if (path_ok) {
      out[[i]] <- candidate
    }
  }
  out
}

read_csv_or_empty <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0L) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
           error = function(e) data.frame())
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "")
}

fetch_globi_page <- function(query_name, skip, limit, destination) {
  fields <- c("source_taxon_name", "interaction_type", "target_taxon_external_id",
              "target_taxon_name", "target_taxon_path")
  params <- c(
    paste0("sourceTaxon=", utils::URLencode(query_name, reserved = TRUE)),
    "interactionType=eats",
    "includeObservations=false",
    paste0("limit=", limit),
    paste0("skip=", skip),
    paste0("field=", fields)
  )
  url <- paste0("https://api.globalbioticinteractions.org/interaction.csv?",
                paste(params, collapse = "&"))
  utils::download.file(url, destination, mode = "wb", quiet = TRUE, method = "libcurl")
  read_csv_or_empty(destination)
}

download_one_consumer <- function(query_name, cache_dir, page_size = 10000L,
                                  max_retries = 4L, refresh = FALSE) {
  cache_file <- file.path(cache_dir, paste0(safe_name(query_name), ".csv"))
  if (!refresh && file.exists(cache_file)) {
    cached <- read_csv_or_empty(cache_file)
    return(data.frame(query_name = query_name, status = "cached", rows = nrow(cached),
                      pages = NA_integer_, message = "", stringsAsFactors = FALSE))
  }

  pages <- list()
  skip <- 0L
  page_number <- 0L
  status <- "ok"
  note <- ""
  repeat {
    page_number <- page_number + 1L
    page_file <- tempfile(fileext = ".csv")
    page <- NULL
    last_error <- NULL
    for (attempt in seq_len(max_retries)) {
      page <- tryCatch(fetch_globi_page(query_name, skip, page_size, page_file),
                       error = function(e) { last_error <<- conditionMessage(e); NULL })
      if (!is.null(page)) break
      Sys.sleep(min(2^(attempt - 1L), 8L))
    }
    unlink(page_file)
    if (is.null(page)) {
      status <- "failed"
      note <- if (is.null(last_error)) "Unknown download error" else last_error
      break
    }
    if (nrow(page) == 0L) break
    pages[[length(pages) + 1L]] <- page
    if (nrow(page) < page_size) break
    skip <- skip + page_size
  }

  if (status == "ok") {
    result <- if (length(pages)) do.call(rbind, pages) else data.frame(
      source_taxon_name = character(), interaction_type = character(),
      target_taxon_external_id = character(), target_taxon_name = character(),
      target_taxon_path = character(), stringsAsFactors = FALSE)
    result$query_name <- rep(query_name, nrow(result))
    result <- result[, c("query_name", setdiff(names(result), "query_name")), drop = FALSE]
    write_csv(result, cache_file)
  }

  data.frame(query_name = query_name, status = status,
             rows = if (status == "ok") nrow(result) else 0L,
             pages = page_number, message = note, stringsAsFactors = FALSE)
}

summarize_globi <- function(consumers, cache_dir) {
  raw <- lapply(consumers, function(x) {
    path <- file.path(cache_dir, paste0(safe_name(x), ".csv"))
    dat <- read_csv_or_empty(path)
    if (!nrow(dat)) return(NULL)
    dat$query_name <- x
    dat
  })
  raw <- Filter(Negate(is.null), raw)
  links <- if (length(raw)) do.call(rbind, raw) else data.frame()

  required <- c("query_name", "source_taxon_name", "interaction_type",
                "target_taxon_external_id", "target_taxon_name", "target_taxon_path")
  if (!nrow(links)) {
    links <- as.data.frame(setNames(replicate(length(required), character(), simplify = FALSE), required),
                           stringsAsFactors = FALSE)
  }
  for (nm in setdiff(required, names(links))) links[[nm]] <- ""
  links <- links[, required, drop = FALSE]
  links$consumer <- canonical_binomial(links$query_name)
  links$reported_consumer <- canonical_binomial(links$source_taxon_name)
  links$resource_name <- trimws(as.character(links$target_taxon_name))
  links$resource_species <- species_from_globi(
    links$target_taxon_name, links$target_taxon_path, links$target_taxon_external_id)
  links$species_level <- !is.na(links$resource_species)
  links$self_link <- links$species_level & links$resource_species == links$consumer
  links$resource_key_all <- ifelse(
    links$species_level,
    tolower(links$resource_species),
    tolower(gsub("\\s+", " ", trimws(links$resource_name)))
  )
  links <- links[nzchar(links$resource_key_all), , drop = FALSE]
  keep <- !duplicated(links[, c("consumer", "resource_key_all")])
  links <- links[keep, , drop = FALSE]
  links <- links[order(links$consumer, links$resource_key_all), , drop = FALSE]

  degree <- data.frame(consumer = consumers, stringsAsFactors = FALSE)
  degree$degree_all_taxa <- vapply(degree$consumer, function(sp) {
    length(unique(links$resource_key_all[links$consumer == sp & !links$self_link]))
  }, integer(1))
  degree$degree_species <- vapply(degree$consumer, function(sp) {
    length(unique(links$resource_species[links$consumer == sp & links$species_level & !links$self_link]))
  }, integer(1))
  degree$self_links <- vapply(degree$consumer, function(sp) {
    sum(links$consumer == sp & links$self_link)
  }, integer(1))
  degree$reported_source_names <- vapply(degree$consumer, function(sp) {
    length(unique(links$source_taxon_name[links$consumer == sp]))
  }, integer(1))
  list(links = links, degree = degree)
}

process_tetraeu <- function(path) {
  raw <- read.delim(path, sep = ";", quote = "\"", stringsAsFactors = FALSE,
                    check.names = FALSE, na.strings = c("", "NA"))
  consumers <- sort(unique(na.omit(canonical_binomial(raw$sourceTaxonName))))
  links <- data.frame(
    consumer = canonical_binomial(raw$sourceTaxonName),
    resource_species = canonical_binomial(raw$targetTaxonName),
    interaction_type = as.character(raw$interactionTypeName),
    stringsAsFactors = FALSE
  )
  links <- links[!is.na(links$consumer) & !is.na(links$resource_species), , drop = FALSE]
  links$self_link <- links$consumer == links$resource_species
  links <- unique(links)
  links <- links[order(links$consumer, links$resource_species), , drop = FALSE]
  degree <- data.frame(consumer = consumers, stringsAsFactors = FALSE)
  degree$degree_species <- vapply(degree$consumer, function(sp) {
    length(unique(links$resource_species[links$consumer == sp & !links$self_link]))
  }, integer(1))
  degree$degree_all_taxa <- degree$degree_species
  degree$self_links <- vapply(degree$consumer, function(sp) {
    sum(links$consumer == sp & links$self_link)
  }, integer(1))
  list(links = links, degree = degree)
}
