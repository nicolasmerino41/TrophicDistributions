library(tidyverse)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
} else {
  normalizePath("sensitivity", mustWork = TRUE)
}

output_dir <- file.path(script_dir, "Outputs")
scenario_dir <- file.path(output_dir, "scenarios")
scenario_files <- list.files(scenario_dir, pattern = "\\.tsv$", full.names = TRUE)
if (!length(scenario_files)) {
  stop("No sensitivity scenarios found. Run sensitivity/MainSensitivity.jl first.")
}

surfaces <- map_dfr(scenario_files, ~read.delim(
  .x, na.strings = "NA", check.names = FALSE,
  colClasses = c(sensitivity_label = "character")
))
baseline <- surfaces %>%
  filter(scenario_id == "baseline") %>%
  select(environment, regime_id, target_r, degree, baseline_mismatch = mean_mismatch)
if (!nrow(baseline)) stop("The baseline scenario is missing.")

effects <- surfaces %>%
  filter(scenario_id != "baseline") %>%
  inner_join(baseline, by = c("environment", "regime_id", "target_r", "degree")) %>%
  mutate(change = mean_mismatch - baseline_mismatch)

summary <- effects %>%
  group_by(sensitivity_variable, sensitivity_value, sensitivity_label) %>%
  summarise(
    n_framework_cells = sum(is.finite(change)),
    median_change = median(change, na.rm = TRUE),
    q05_change = quantile(change, 0.05, na.rm = TRUE),
    q25_change = quantile(change, 0.25, na.rm = TRUE),
    q75_change = quantile(change, 0.75, na.rm = TRUE),
    q95_change = quantile(change, 0.95, na.rm = TRUE),
    median_absolute_change = median(abs(change), na.rm = TRUE),
    maximum_absolute_change = max(abs(change), na.rm = TRUE),
    surface_correlation = cor(mean_mismatch, baseline_mismatch, use = "complete.obs"),
    .groups = "drop"
  )

baseline_values <- tribble(
  ~sensitivity_variable, ~sensitivity_value, ~sensitivity_label,
  "Grid dimension", 60, "60 x 60",
  "Minimum patch", 50 / 3600, "1.39%",
  "Suitability threshold", 0.25, "0.25",
  "Species richness", 250, "250"
) %>%
  mutate(
    n_framework_cells = nrow(baseline),
    median_change = 0, q05_change = 0, q25_change = 0,
    q75_change = 0, q95_change = 0,
    median_absolute_change = 0, maximum_absolute_change = 0,
    surface_correlation = 1
  )

summary <- bind_rows(summary, baseline_values) %>%
  mutate(
    sensitivity_variable = factor(
      sensitivity_variable,
      levels = c("Grid dimension", "Minimum patch", "Suitability threshold", "Species richness")
    ),
    baseline = median_absolute_change == 0 & maximum_absolute_change == 0
  ) %>%
  arrange(sensitivity_variable, sensitivity_value)

plot_data <- summary %>%
  group_by(sensitivity_variable) %>%
  arrange(sensitivity_value, .by_group = TRUE) %>%
  mutate(value_display = factor(sensitivity_label, levels = sensitivity_label)) %>%
  ungroup()

sensitivity_plot <- ggplot(plot_data, aes(x = value_display, y = median_change)) +
  geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.45, linetype = 2) +
  geom_linerange(aes(ymin = q25_change, ymax = q75_change), linewidth = 2.2, colour = "#2C7FB8") +
  geom_point(aes(fill = baseline), shape = 21, size = 3.1, stroke = 0.6, colour = "black") +
  facet_wrap(~sensitivity_variable, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = c(`TRUE` = "white", `FALSE` = "#2C7FB8"), guide = "none") +
  labs(
    x = "Parameter value",
    y = "Change in mean consumer mismatch relative to baseline"
  ) +
  theme_classic(base_family = "Arial", base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.title = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(12, 18, 12, 18)
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(summary, file.path(output_dir, "sensitivity_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
ggsave(file.path(output_dir, "sensitivity_secondary_variables.png"),
       sensitivity_plot, width = 13, height = 5.3, dpi = 600)
cat("Saved sensitivity summary and combined PNG to:", output_dir, "\n")
