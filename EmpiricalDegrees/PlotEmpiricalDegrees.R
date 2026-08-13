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
output_file <- file.path(root_dir, "Outputs", "empirical_degrees", "empirical_degree_distributions.png")

dat <- read.csv(input_file, stringsAsFactors = FALSE)
dat <- dat[is.finite(dat$degree) & dat$degree >= 0, , drop = FALSE]
groups <- c("Thermal consumers (GloBI)", "TetraEU consumers")
cols <- c("#D55E00", "#0072B2")
degree_classes <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55)

if (.Platform$OS.type == "windows") grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
png(output_file, width = 2600, height = 1200, res = 240, family = "Arial")
layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 1))
par(mar = c(5.2, 8.2, 1.3, 1.2), mgp = c(3.0, 0.8, 0), las = 1,
    cex.axis = 1.05, cex.lab = 1.2)

# Panel A: all consumers, with degree shown on a log2(x + 1) scale.
ticks <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
xlim <- range(log2(dat$degree + 1), finite = TRUE)
plot(NA, xlim = xlim, ylim = c(0.55, 2.45), axes = FALSE,
     xlab = "Recorded consumer degree", ylab = "")
axis(1, at = log2(ticks + 1), labels = ticks)
axis(2, at = 1:2, labels = c("Thermal consumers\n(GloBI)", "TetraEU consumers"),
     tick = FALSE, line = -0.4)
box(bty = "l")
set.seed(1)
for (i in seq_along(groups)) {
  x <- dat$degree[dat$dataset == groups[[i]]]
  y <- jitter(rep(i, length(x)), amount = 0.22)
  points(log2(x + 1), y, pch = 16, col = grDevices::adjustcolor(cols[[i]], alpha.f = 0.28), cex = 0.65)
  qs <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  segments(log2(qs[[1]] + 1), i, log2(qs[[3]] + 1), i, lwd = 8, col = cols[[i]])
  points(log2(qs[[2]] + 1), i, pch = 21, bg = "white", col = cols[[i]], lwd = 2, cex = 1.15)
  text(xlim[[2]], i + 0.28, paste0("n = ", length(x)), adj = 1, cex = 0.9)
}
mtext("A", side = 3, adj = 0, line = 0.1, font = 2, cex = 1.35)

# Panel B: no specialist cutoff is imposed; the curve reports every possible cutoff.
plot(NA, xlim = range(degree_classes), ylim = c(0, 1), log = "x", axes = FALSE,
     xlab = "Degree threshold", ylab = "Proportion of consumers above threshold")
axis(1, at = degree_classes, labels = degree_classes, las = 2, cex.axis = 0.78)
axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"))
box(bty = "l")
for (i in seq_along(groups)) {
  x <- dat$degree[dat$dataset == groups[[i]]]
  prop <- vapply(degree_classes, function(k) mean(x > k), numeric(1))
  lines(degree_classes, prop, type = "o", pch = 16, lwd = 2.2, cex = 0.75, col = cols[[i]])
}
legend("bottomleft", legend = groups, col = cols, lwd = 2.2, pch = 16,
       bty = "n", cex = 0.9)
mtext("B", side = 3, adj = 0, line = 0.1, font = 2, cex = 1.35)
dev.off()

cat("Saved:", output_file, "\n")
