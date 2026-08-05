script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
source(file.path(script_dir, "PlotFunctions.R"))

output_dir <- Sys.getenv("SDM_OUTPUT_DIR", file.path(script_dir, "Outputs"))
summary <- read_sdm_table(file.path(output_dir, "final_summary.tsv"))
surface_dir <- file.path(output_dir, "all_surfaces")
dir.create(surface_dir, recursive = TRUE, showWarnings = FALSE)
stages <- c("oracle", "estimated_100", "estimated_75", "estimated_50", "estimated_25")
limit <- max(abs(summary$mean_delta_auc), na.rm = TRUE)
limit <- max(limit, 1e-6)
palette <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100)

treatments <- unique(summary[c("environment", "network", "regime")])
for (i in seq_len(nrow(treatments))) {
  treatment <- treatments[i, ]
  keep <- summary$environment == treatment$environment &
    summary$network == treatment$network & summary$regime == treatment$regime
  data <- summary[keep, ]
  safe_name <- gsub("[^A-Za-z0-9]+", "_",
                    paste(treatment$environment, treatment$network, treatment$regime, sep = "_"))
  output <- file.path(surface_dir, paste0(safe_name, ".png"))
  png(output, width = 3000, height = 720, res = 180)
  par(mfrow = c(1, 5), mar = c(4.8, 4.2, 3.3, 1.0), mgp = c(2.4, 0.7, 0), tcl = -0.25)
  for (stage in stages) {
    draw_heatmap(data[data$information_level == stage, ], "mean_delta_auc",
                 stage_label(stage), c(-limit, limit), palette,
                 show_y = stage == stages[1])
  }
  mtext(paste(treatment$environment, treatment$network, treatment$regime, sep = " | "),
        outer = TRUE, line = -1.2, font = 2)
  dev.off()
}
cat("Saved", nrow(treatments), "supplementary surface figures to", surface_dir, "\n")
