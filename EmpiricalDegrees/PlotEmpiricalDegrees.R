script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path <- sub("^--file=", "", script_arg[[1]])
} else {
  source_paths <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)
  }, character(1))
  source_paths <- source_paths[nzchar(source_paths)]
  candidates <- c(
    if (length(source_paths)) tail(source_paths, 1) else character(),
    file.path(getwd(), "EmpiricalDegrees", "PlotEmpiricalDegrees.R"),
    file.path(getwd(), "PlotEmpiricalDegrees.R")
  )
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not locate PlotEmpiricalDegrees.R. Run the complete script with source() or Rscript.")
  }
  script_path <- candidates[[1]]
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_file <- file.path(root_dir, "Outputs", "empirical_degrees", "consumer_degrees_combined.csv")
output_file <- file.path(root_dir, "Outputs", "empirical_degrees", "Figure4.png")

dat <- read.csv(input_file, stringsAsFactors = FALSE)
dat <- dat[is.finite(dat$degree) & dat$degree > 0, , drop = FALSE]
dat$dataset[dat$dataset == "Thermal consumers (GloBI)"] <- "GloBI consumers"
groups <- c("GloBI consumers", "TetraEU consumers")
panel_a_labels <- c("GloBI\nconsumers", "TetraEU\nconsumers")
# Two well-separated colours from the Viridis palette.
cols <- c("#440154", "#21918C")
degree_classes <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55)

if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
png(output_file, width = 2800, height = 1250, res = 240, family = "Arial")
layout(matrix(c(1, 2), nrow = 1), widths = c(1.12, 1.00))
par(mar = c(5.0, 5.2, 1.6, 1.6), mgp = c(2.6, 0.65, 0), las = 1,
    cex.axis = 0.90, cex.lab = 1.05)

# Panel A: positive recorded degrees on a conventional logarithmic axis.
ticks <- 2^(0:10)
xlim <- c(min(dat$degree) / 1.12, max(dat$degree) * 1.12)
plot(NA, xlim = xlim, ylim = c(0.62, 2.38), log = "x", axes = FALSE,
     xlab = "Consumer in-degree (log scale)", ylab = "")
axis(1, at = ticks[ticks <= max(dat$degree)],
     labels = ticks[ticks <= max(dat$degree)], gap.axis = 0.15)
axis(
  2,
  at = 1:2,
  labels = panel_a_labels,
  tick = FALSE,
  line = -0.4
)
box(bty = "l")
set.seed(1)
for (i in seq_along(groups)) {
  x <- dat$degree[dat$dataset == groups[[i]]]
  y <- jitter(rep(i, length(x)), amount = 0.19)
  points(x, y, pch = 16,
         col = grDevices::adjustcolor(cols[[i]], alpha.f = 0.34), cex = 0.62)
  qs <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  segments(qs[[1]], i, qs[[3]], i, lwd = 7, col = cols[[i]])
  points(qs[[2]], i, pch = 21, bg = "white", col = cols[[i]], lwd = 2, cex = 1.08)
  text(max(dat$degree), i + 0.24, paste0("n = ", length(x)), adj = 1, cex = 0.88)
}
mtext("A", side = 3, adj = 0, line = 0.1, font = 2, cex = 1.35)

# Panel B: degree classes are categorical positions and therefore evenly spaced.
par(mar = c(5.0, 5.3, 1.6, 1.2), mgp = c(2.6, 0.65, 0), las = 1,
    cex.axis = 0.90, cex.lab = 1.05)
x_positions <- seq_along(degree_classes)
plot(NA, xlim = c(0.55, length(degree_classes) + 0.45), ylim = c(0, 1), axes = FALSE,
     xlab = "In-degree", ylab = "Proportion of consumers ≥ n in-degree")
axis(1, at = x_positions, labels = degree_classes, las = 2, cex.axis = 0.80)
axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"))
box(bty = "l")
for (i in seq_along(groups)) {
  x <- dat$degree[dat$dataset == groups[[i]]]
  prop <- vapply(degree_classes, function(k) mean(x >= k), numeric(1))
  lines(x_positions, prop, type = "o", pch = 16,
        lwd = 2.3, cex = 0.72, col = cols[[i]])
}
legend(
  "bottomleft",
  legend = groups,
  col = cols,
  lwd = 2.2,
  pch = 16,
  bty = "n",
  cex = 0.88
)
mtext("B", side = 3, adj = 0, line = 0.1, font = 2, cex = 1.35)
dev.off()

cat("Saved:", output_file, "\n")
