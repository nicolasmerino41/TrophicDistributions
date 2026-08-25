library(tidyverse)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
} else {
  normalizePath("SimulationsCode", mustWork = TRUE)
}

repository_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_file <- file.path(
  repository_dir, "Outputs", "degreeResults", "degree_replicate_summary.tsv"
)
output_dir <- file.path(repository_dir, "Outputs", "scenarioContrasts")
cell_file <- file.path(output_dir, "scenario_contrast_cells.tsv")
summary_file <- file.path(output_dir, "scenario_contrast_summary.tsv")
figure_file <- file.path(output_dir, "scenario_contrasts.png")

bootstrap_iterations <- 2000L
bootstrap_seed <- 20260825L

if (!file.exists(input_file)) {
  stop("Missing replicate summary: ", input_file,
       ". Run SimulationsCode/MainScript.jl first.")
}

data <- read.delim(input_file, na.strings = c("NA", "NaN"), check.names = FALSE) %>%
  filter(is.finite(mean_mismatch), n_eligible > 0) %>%
  mutate(
    environment = as.character(environment),
    regime = as.character(regime),
    target_r = as.numeric(target_r),
    degree = as.numeric(degree)
  )

environments <- c("random", "autocorr")
regimes <- c(
  "Broad + HighVar", "Broad + LowVar",
  "Narrow + HighVar", "Narrow + LowVar"
)
correlations <- sort(unique(data$target_r))
degrees <- sort(unique(data$degree))

stopifnot(
  all(environments %in% unique(data$environment)),
  all(regimes %in% unique(data$regime))
)

observed <- array(
  NA_real_,
  dim = c(length(environments), length(regimes), length(correlations), length(degrees)),
  dimnames = list(
    environment = environments,
    regime = regimes,
    target_r = as.character(correlations),
    degree = as.character(degrees)
  )
)

bootstrapped <- array(
  NA_real_,
  dim = c(
    bootstrap_iterations, length(environments), length(regimes),
    length(correlations), length(degrees)
  ),
  dimnames = list(
    bootstrap = NULL,
    environment = environments,
    regime = regimes,
    target_r = as.character(correlations),
    degree = as.character(degrees)
  )
)

set.seed(bootstrap_seed)

for (environment_value in environments) {
  for (regime_value in regimes) {
    for (correlation_value in correlations) {
      stratum <- data %>%
        filter(
          environment == environment_value,
          regime == regime_value,
          target_r == correlation_value
        )

      community_ids <- sort(unique(stratum$community_id))
      community_matrix <- matrix(
        NA_real_, nrow = length(community_ids), ncol = length(degrees),
        dimnames = list(as.character(community_ids), as.character(degrees))
      )
      row_index <- match(stratum$community_id, community_ids)
      column_index <- match(stratum$degree, degrees)
      community_matrix[cbind(row_index, column_index)] <- stratum$mean_mismatch

      cell_means <- colMeans(community_matrix, na.rm = TRUE)
      cell_means[!is.finite(cell_means)] <- NA_real_
      observed[environment_value, regime_value, as.character(correlation_value), ] <-
        cell_means

      n_communities <- nrow(community_matrix)
      bootstrap_weights <- t(rmultinom(
        bootstrap_iterations,
        size = n_communities,
        prob = rep(1 / n_communities, n_communities)
      ))
      valid <- is.finite(community_matrix)
      values <- community_matrix
      values[!valid] <- 0
      numerator <- bootstrap_weights %*% values
      denominator <- bootstrap_weights %*% (valid * 1)
      bootstrap_means <- numerator / denominator
      bootstrap_means[denominator == 0] <- NA_real_

      bootstrapped[, environment_value, regime_value,
                   as.character(correlation_value), ] <- bootstrap_means
    }
  }
}

cell_results <- list()
summary_results <- list()

add_contrast <- function(group, label, observed_difference,
                         bootstrap_difference, metadata) {
  cell_values <- as.vector(observed_difference)
  stopifnot(nrow(metadata) == length(cell_values))
  cell_table <- metadata %>%
    mutate(group = group, contrast = label, difference = cell_values) %>%
    select(group, contrast, everything()) %>%
    filter(is.finite(difference))

  bootstrap_values <- apply(
    bootstrap_difference, 1,
    function(values) mean(values, na.rm = TRUE)
  )
  bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]

  cell_results[[length(cell_results) + 1L]] <<- cell_table
  summary_results[[length(summary_results) + 1L]] <<- tibble(
    group = group,
    contrast = label,
    estimate = mean(cell_table$difference),
    lower_95 = unname(quantile(bootstrap_values, 0.025)),
    upper_95 = unname(quantile(bootstrap_values, 0.975)),
    proportion_positive = mean(cell_table$difference > 0),
    n_matched_cells = nrow(cell_table),
    bootstrap_iterations = length(bootstrap_values)
  )
}

niche_metadata <- expand.grid(
  environment = environments,
  target_r = correlations,
  degree = degrees,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% mutate(regime = NA_character_)

add_contrast(
  "Niche breadth",
  "Narrow - Broad (high variance)",
  observed[, "Narrow + HighVar", , ] - observed[, "Broad + HighVar", , ],
  bootstrapped[, , "Narrow + HighVar", , , drop = FALSE] -
    bootstrapped[, , "Broad + HighVar", , , drop = FALSE],
  niche_metadata
)

add_contrast(
  "Niche breadth",
  "Narrow - Broad (low variance)",
  observed[, "Narrow + LowVar", , ] - observed[, "Broad + LowVar", , ],
  bootstrapped[, , "Narrow + LowVar", , , drop = FALSE] -
    bootstrapped[, , "Broad + LowVar", , , drop = FALSE],
  niche_metadata
)

add_contrast(
  "Niche breadth",
  "LowVar - HighVar (broad niches)",
  observed[, "Broad + LowVar", , ] - observed[, "Broad + HighVar", , ],
  bootstrapped[, , "Broad + LowVar", , , drop = FALSE] -
    bootstrapped[, , "Broad + HighVar", , , drop = FALSE],
  niche_metadata
)

add_contrast(
  "Niche breadth",
  "LowVar - HighVar (narrow niches)",
  observed[, "Narrow + LowVar", , ] - observed[, "Narrow + HighVar", , ],
  bootstrapped[, , "Narrow + LowVar", , , drop = FALSE] -
    bootstrapped[, , "Narrow + HighVar", , , drop = FALSE],
  niche_metadata
)

for (regime_value in regimes) {
  spatial_metadata <- expand.grid(
    target_r = correlations,
    degree = degrees,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>% mutate(environment = NA_character_, regime = regime_value)

  add_contrast(
    "Spatial structure",
    paste0("Random - Autocorrelated (", regime_value, ")"),
    observed["random", regime_value, , ] - observed["autocorr", regime_value, , ],
    bootstrapped[, "random", regime_value, , , drop = FALSE] -
      bootstrapped[, "autocorr", regime_value, , , drop = FALSE],
    spatial_metadata
  )
}

overall_spatial_metadata <- expand.grid(
  regime = regimes,
  target_r = correlations,
  degree = degrees,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>% mutate(environment = NA_character_)

add_contrast(
  "Spatial structure",
  "Random - Autocorrelated (overall)",
  observed["random", , , ] - observed["autocorr", , , ],
  bootstrapped[, "random", , , , drop = FALSE] -
    bootstrapped[, "autocorr", , , , drop = FALSE],
  overall_spatial_metadata
)

cell_table <- bind_rows(cell_results)
summary_table <- bind_rows(summary_results)

contrast_order <- summary_table$contrast
cell_table <- cell_table %>%
  mutate(contrast = factor(contrast, levels = rev(contrast_order)))
summary_table <- summary_table %>%
  mutate(contrast = factor(contrast, levels = rev(contrast_order)))

group_colours <- c("Niche breadth" = "#440154", "Spatial structure" = "#21918C")

contrast_plot <- ggplot() +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.5) +
  geom_point(
    data = cell_table,
    aes(x = difference, y = contrast, colour = group),
    position = position_jitter(width = 0, height = 0.12, seed = 1),
    alpha = 0.12, size = 0.65
  ) +
  geom_errorbar(
    data = summary_table,
    aes(xmin = lower_95, xmax = upper_95, y = contrast),
    orientation = "y", width = 0.15, colour = "black", linewidth = 0.9
  ) +
  geom_point(
    data = summary_table,
    aes(x = estimate, y = contrast, fill = group),
    shape = 21, colour = "white", stroke = 0.7, size = 2.8
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = group_colours, guide = "none") +
  scale_fill_manual(values = group_colours, guide = "none") +
  labs(
    x = "Difference in mean consumer Jaccard mismatch\n(first treatment - second treatment)",
    y = NULL
  ) +
  theme_classic(base_family = "Arial", base_size = 11) +
  theme(
    plot.margin = margin(12, 20, 12, 12),
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 0, size = 11),
    panel.spacing.y = unit(1.1, "lines"),
    axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
    axis.text.y = element_text(size = 9.5, colour = "black"),
    axis.text.x = element_text(size = 9.5, colour = "black")
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(
  cell_table %>% mutate(contrast = as.character(contrast)),
  cell_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
write.table(
  summary_table %>% mutate(contrast = as.character(contrast)),
  summary_file, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
ggsave(figure_file, contrast_plot, width = 10.5, height = 7.2, dpi = 600)

cat("Saved scenario contrasts to:", output_dir, "\n")
