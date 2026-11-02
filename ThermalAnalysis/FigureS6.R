# CTmax alignment across terciles of recorded GloBI consumer degree.
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script <- if (length(args)) sub("^--file=", "", args[1]) else {
  files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  if (length(files)) tail(files, 1)[[1]] else "ThermalAnalysis/FigureS6.R"
}
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
out <- file.path(root, "Outputs", "thermal_metrics", "FigureS6.png")
thermal <- read.csv(file.path(root, "Outputs/thermal_metrics/ctmax/predator_level.csv"))
degrees <- read.csv(file.path(root, "Outputs/empirical_degrees/consumer_degrees_combined.csv"))
degrees <- subset(degrees, dataset == "GloBI consumers", select = c(consumer, degree))
stopifnot(!anyDuplicated(degrees$consumer), !anyDuplicated(thermal$pred))
d <- merge(thermal, degrees, by.x = "pred", by.y = "consumer", all.x = TRUE)
stopifnot(all(is.finite(d$degree)))
# Exclude inconsistent zero-degree records (currently Phoxinus lumaireul).
if (any(d$degree <= 0)) message("Excluded nonpositive degree: ", paste(d$pred[d$degree <= 0], collapse = ", "))
d <- subset(d, degree > 0)
cuts <- quantile(d$degree, c(1/3, 2/3))
d$group <- cut(d$degree, c(-Inf, cuts, Inf), labels = c("Low", "Intermediate", "High"))
set.seed(20261102)
results <- do.call(rbind, lapply(levels(d$group), function(g) {
  x <- subset(d, group == g)
  boot <- replicate(5000, {
    i <- sample.int(nrow(x), replace = TRUE)
    cor(x$pred_trait[i], x$mean_prey_trait[i])
  })
  ci <- quantile(boot, c(.025, .975), na.rm = TRUE)
  data.frame(group = g, N = nrow(x), minimum = min(x$degree), maximum = max(x$degree),
             r = cor(x$pred_trait, x$mean_prey_trait), lower = ci[1], upper = ci[2])
}))
print(results, row.names = FALSE)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
png(out, width = 1600, height = 1000, res = 180)
par(mar = c(4.5, 8, 3, 2))
y <- 3:1
plot(results$r, y, xlim = c(0, 1), ylim = c(.6, 3.4), yaxt = "n", ylab = "",
     xlab = "Pearson correlation (consumer vs mean resource CTmax)", pch = 19, cex = 1.2)
axis(2, at = y, labels = paste0(results$group, "\ndegree ", results$minimum, "-",
     results$maximum, "; n=", results$N), las = 1, tick = FALSE, cex.axis = .8)
arrows(results$lower, y, results$upper, y, angle = 90, code = 3, length = .05, lwd = 1.5)
title("Thermal alignment across consumer in-degree")
dev.off()
cat("Saved", out, "\n")
