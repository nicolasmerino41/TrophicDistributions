library(tidyverse)
library(viridis)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
} else {
  normalizePath("SimulationsCode", mustWork = TRUE)
}

repository_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
degree_file <- file.path(repository_dir, "Outputs", "degreeResults", "degree_summary.tsv")
figure_dir <- file.path(repository_dir, "Outputs", "heatmaps")
figure_file <- file.path(figure_dir, "degree_heatmap_mean_mismatch.png")

if (!file.exists(degree_file)) {
  stop("Missing degree summary: ", degree_file, ". Run SimulationsCode/MainScript.jl first.")
}

degree_levels <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55)
environment_labels <- c(random = "Random environment", autocorr = "Autocorrelated environment")

degree_df <- read.delim(degree_file, na.strings = "NA", check.names = FALSE) %>%
  mutate(
    degree = factor(degree, levels = degree_levels, ordered = TRUE),
    environment = factor(environment, levels = names(environment_labels)),
    environment_label = factor(
      environment_labels[as.character(environment)],
      levels = environment_labels
    )
  )

heatmap <- ggplot(degree_df, aes(x = degree, y = target_r, fill = mean_mismatch)) +
  geom_tile() +
  facet_grid(environment_label ~ regime) +
  scale_x_discrete(drop = FALSE, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    na.value = "grey90",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(8, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(
    x = "Consumer in-degree",
    y = "Consumer-resource niche correlation",
    fill = "Mean consumer Jaccard mismatch"
  ) +
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

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(figure_file, heatmap, width = 13, height = 7.5, dpi = 600)
cat("Saved heatmap to:", figure_file, "\n")
