library(tidyverse)
library(viridis)

OUTROOT <- "Outputs"
DEGREE_FILE <- file.path(OUTROOT, "degreeResults", "degree_summary.tsv")
FIG_DIR <- file.path(OUTROOT, "heatmaps")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

DEGREE_LEVELS <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55)
ENV_LEVELS <- c("random", "autocorr")
NET_LEVELS <- c("random", "modular", "cascade")

ENV_LABELS <- c(random = "Random environment", autocorr = "Autocorrelated environment")
NET_LABELS <- c(random = "Random", modular = "Modular", cascade = "Cascade")

METRIC_LABELS <- c(
  mean_mismatch = "Mean consumer Jaccard mismatch",
  mean_q90_mismatch = "Mean within-replicate q90 mismatch",
  mean_frac_gt = "Mean fraction above mismatch threshold",
  mean_degree_correlation = "Mean degree-specific niche correlation",
  mean_max_degree_correlation_error = "Mean maximum degree-correlation error"
)

if (!file.exists(DEGREE_FILE)) {
  stop("Missing degree summary: ", DEGREE_FILE, ". Run SimulationsCode/MainScript.jl first.")
}

degree_df <- read.delim(DEGREE_FILE, na.strings = "NA", check.names = FALSE) %>%
  mutate(
    environment = factor(environment, levels = ENV_LEVELS),
    network = factor(network, levels = NET_LEVELS),
    degree = factor(degree, levels = DEGREE_LEVELS, ordered = TRUE),
    EnvironmentLabel = factor(
      ENV_LABELS[as.character(environment)],
      levels = ENV_LABELS[ENV_LEVELS]
    ),
    NetworkLabel = factor(
      NET_LABELS[as.character(network)],
      levels = NET_LABELS[NET_LEVELS]
    )
  )

theme_heat <- function() {
  theme_classic(base_family = "Arial", base_size = 12) +
    theme(
      plot.margin = margin(12, 18, 12, 18),
      strip.text = element_text(size = 11),
      strip.background = element_blank(),
      panel.spacing = unit(0.9, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 9, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10, color = "black"),
      legend.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
}

plot_degree_heatmap <- function(metric_name) {
  ggplot(degree_df, aes(x = degree, y = target_r, fill = .data[[metric_name]])) +
    geom_tile() +
    facet_grid(EnvironmentLabel + NetworkLabel ~ regime) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "viridis",
      na.value = "grey90",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "top",
        barwidth = unit(8, "cm"),
        barheight = unit(0.5, "cm")
      )
    ) +
    labs(
      x = "Focal consumer degree",
      y = "Community-level niche correlation",
      fill = METRIC_LABELS[[metric_name]]
    ) +
    theme_heat()
}

available_metrics <- names(METRIC_LABELS)[names(METRIC_LABELS) %in% names(degree_df)]
for (metric_name in available_metrics) {
  plot <- plot_degree_heatmap(metric_name)
  ggsave(
    filename = file.path(FIG_DIR, paste0("degree_heatmap_", metric_name, ".png")),
    plot = plot,
    width = 13,
    height = 14,
    dpi = 600
  )
}

if ("mean_mismatch" %in% names(degree_df)) {
  difference_df <- degree_df %>%
    select(environment, NetworkLabel, regime, target_r, degree, mean_mismatch) %>%
    pivot_wider(names_from = environment, values_from = mean_mismatch) %>%
    mutate(value = random - autocorr)

  limit <- max(abs(difference_df$value), na.rm = TRUE)
  difference_plot <- ggplot(
    difference_df,
    aes(x = degree, y = target_r, fill = value)
  ) +
    geom_tile() +
    facet_grid(NetworkLabel ~ regime) +
    scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-limit, limit), oob = scales::squish,
      guide = guide_colorbar(direction = "horizontal", title.position = "top")
    ) +
    labs(
      x = "Focal consumer degree",
      y = "Community-level niche correlation",
      fill = "Mean mismatch\n(random - autocorrelated)"
    ) +
    theme_heat()

  ggsave(
    filename = file.path(FIG_DIR, "degree_difference_mean_mismatch.png"),
    plot = difference_plot,
    width = 12,
    height = 9,
    dpi = 600
  )
}

cat("Saved degree-based heatmaps to:", FIG_DIR, "\n")
