find_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg)) {
    candidate <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE))
    if (file.exists(file.path(candidate, "PlotFunctions.R"))) return(candidate)
  }
  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(frame_files)) {
    candidate <- dirname(normalizePath(tail(frame_files, 1)[[1]], mustWork = TRUE))
    if (file.exists(file.path(candidate, "PlotFunctions.R"))) return(candidate)
  }
  candidates <- c(getwd(), file.path(getwd(), "sdm"))
  matches <- candidates[file.exists(file.path(candidates, "PlotFunctions.R"))]
  if (length(matches)) return(normalizePath(matches[1], mustWork = TRUE))
  stop("Cannot locate the sdm folder containing PlotFunctions.R")
}

script_dir <- find_script_dir()
output_dir <- file.path(script_dir, "Outputs")
source(file.path(script_dir, "PlotFunctions.R"))

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
