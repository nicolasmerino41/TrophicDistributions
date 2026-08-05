script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
source(file.path(script_dir, "PlotFunctions.R"))

# Reference treatment shown in the main figure. Change these three values if a
# different architecture/environment/regime should anchor Figure 5.
reference_environment <- "autocorr"
reference_network <- "random"
reference_regime <- "Broad + LowVar"
output_dir <- Sys.getenv("SDM_OUTPUT_DIR", file.path(script_dir, "Outputs"))

sdm_summary <- read_sdm_table(file.path(output_dir, "final_summary.tsv"))
relationship <- read_sdm_table(file.path(output_dir, "relationship_summary.tsv"))
truth_path <- file.path(script_dir, "..", "Outputs", "degreeResults", "degree_summary.tsv")
truth <- read_sdm_table(truth_path)

keep_sdm <- sdm_summary$environment == reference_environment &
  sdm_summary$network == reference_network &
  sdm_summary$regime == reference_regime
sdm_reference <- sdm_summary[keep_sdm, ]
keep_truth <- truth$environment == reference_environment &
  truth$network == reference_network & truth$regime == reference_regime
truth_reference <- truth[keep_truth, ]
if (!nrow(sdm_reference) || !nrow(truth_reference)) {
  stop("The selected Figure 5 reference treatment is absent from the result tables")
}

stages <- c("oracle", "estimated_100", "estimated_75", "estimated_50")
auc_limit <- max(abs(sdm_reference$mean_delta_auc[sdm_reference$information_level %in% stages]),
                 na.rm = TRUE)
auc_limit <- max(auc_limit, 1e-6)
auc_palette <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(100)
mismatch_palette <- hcl.colors(100, "Viridis")

output_png <- file.path(output_dir, "figure5_sdm_application.png")
output_pdf <- file.path(output_dir, "figure5_sdm_application.pdf")

draw_figure <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old))
  layout(matrix(1:6, nrow = 2, byrow = TRUE))
  par(mar = c(4.5, 4.2, 3.2, 1.2), mgp = c(2.3, 0.7, 0), tcl = -0.25)
  draw_heatmap(truth_reference, "mean_mismatch", "Perfect-information result",
               c(0, 1), mismatch_palette)
  add_panel_label("A")
  for (i in seq_along(stages)) {
    stage_data <- sdm_reference[sdm_reference$information_level == stages[i], ]
    draw_heatmap(stage_data, "mean_delta_auc", stage_label(stages[i]),
                 c(-auc_limit, auc_limit), auc_palette,
                 show_y = i == 1)
    add_panel_label(LETTERS[i + 1])
  }
  draw_relationship(relationship)
  add_panel_label("F")
}

png(output_png, width = 3000, height = 1900, res = 220)
draw_figure()
dev.off()
pdf(output_pdf, width = 14, height = 9)
draw_figure()
dev.off()
cat("Saved Figure 5 to:\n", output_png, "\n", output_pdf, "\n", sep = "")
