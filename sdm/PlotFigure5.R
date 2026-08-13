find_script_dir <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)))
  }
  frame_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
  if (length(frame_files)) {
    return(dirname(normalizePath(tail(frame_files, 1)[[1]], mustWork = TRUE)))
  }
  candidates <- c(getwd(), file.path(getwd(), "sdm"))
  matches <- candidates[file.exists(file.path(candidates, "Parameters.jl"))]
  if (length(matches)) return(normalizePath(matches[1], mustWork = TRUE))
  stop("Cannot locate the sdm directory")
}

read_table <- function(path) {
  if (!file.exists(path)) stop("Missing required result file: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

weighted_cell_mean <- function(data, value, weight) {
  groups <- interaction(data$target_r, data$degree, drop = TRUE)
  pieces <- split(data, groups)
  rows <- lapply(pieces, function(x) {
    valid <- is.finite(x[[value]]) & is.finite(x[[weight]]) & x[[weight]] > 0
    data.frame(target_r = x$target_r[1], degree = x$degree[1],
               value = weighted.mean(x[[value]][valid], x[[weight]][valid]))
  })
  do.call(rbind, rows)
}

pool_sdm_cells <- function(data) {
  groups <- interaction(data$target_r, data$degree, drop = TRUE)
  pieces <- split(data, groups)
  rows <- lapply(pieces, function(x) {
    valid <- is.finite(x$mean_delta_auc) & is.finite(x$se_delta_auc) & x$n > 0
    x <- x[valid, ]
    total_n <- sum(x$n)
    grand_mean <- weighted.mean(x$mean_delta_auc, x$n)
    within_sd <- x$se_delta_auc * sqrt(x$n)
    sum_squares <- sum((x$n - 1) * within_sd^2 + x$n * (x$mean_delta_auc - grand_mean)^2)
    pooled_sd <- if (total_n > 1) sqrt(sum_squares / (total_n - 1)) else NA_real_
    data.frame(target_r = x$target_r[1], degree = x$degree[1], n = total_n,
               mean_delta_auc = grand_mean, se_delta_auc = pooled_sd / sqrt(total_n))
  })
  do.call(rbind, rows)
}

make_matrix <- function(data, value, degrees, correlations) {
  result <- matrix(NA_real_, nrow = length(correlations), ncol = length(degrees))
  for (i in seq_len(nrow(data))) {
    result[match(data$target_r[i], correlations), match(data$degree[i], degrees)] <- data[[value]][i]
  }
  result
}

script_dir <- find_script_dir()
repository_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
output_dir <- file.path(script_dir, "Outputs")

truth <- read_table(file.path(repository_dir, "Outputs", "degreeResults", "degree_summary.tsv"))
sdm <- read_table(file.path(output_dir, "oracle_summary.tsv"))

truth_average <- weighted_cell_mean(truth, "mean_mismatch", "n_communities_with_data")
sdm_average <- pool_sdm_cells(sdm)
degrees <- sort(unique(truth_average$degree))
correlations <- sort(unique(truth_average$target_r))

figure_degrees <- c(2, 6)
requested_correlations <- c(min(correlations), max(correlations))
figure_correlations <- vapply(
  requested_correlations,
  function(value) correlations[which.min(abs(correlations - value))],
  numeric(1)
)

selected <- expand.grid(
  degree = figure_degrees,
  target_r = figure_correlations
)
selected <- selected[c(1, 3, 2, 4), ]
selected$requested_correlation <- rep(requested_correlations, each = 2)[c(1, 3, 2, 4)]
selected$label <- c("Degree 2 / correlation 0", "Degree 2 / correlation 0.95",
                    "Degree 6 / correlation 0", "Degree 6 / correlation 0.95")
selected$colour <- c("#D55E00", "#E69F00", "#0072B2", "#009E73")
selected <- merge(selected, sdm_average, by = c("degree", "target_r"), sort = FALSE)
selected$order <- match(paste(selected$degree, selected$target_r),
                        paste(c(2, 2, 6, 6),
                              c(figure_correlations[1], figure_correlations[2],
                                figure_correlations[1], figure_correlations[2])))
selected <- selected[order(selected$order), ]
selected$requested_correlation <- rep(requested_correlations, each = 2)[c(1, 3, 2, 4)]
selected$label <- c("Degree 2 / correlation 0", "Degree 2 / correlation 0.95",
                    "Degree 6 / correlation 0", "Degree 6 / correlation 0.95")
selected$colour <- c("#D55E00", "#E69F00", "#0072B2", "#009E73")

truth_matrix <- make_matrix(truth_average, "value", degrees, correlations)
palette <- hcl.colors(100, "Viridis")
breaks <- seq(0, 1, length.out = length(palette) + 1)

draw_figure <- function() {
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.05, 0.13, 1.0))

  par(mar = c(5.0, 5.0, 3.5, 1.0), mgp = c(2.7, 0.8, 0), tcl = -0.25)
  image(seq_along(degrees), seq_along(correlations), t(truth_matrix),
        col = palette, breaks = breaks, axes = FALSE, useRaster = TRUE,
        xlab = "Focal consumer degree", ylab = bquote(paste("Consumer-resource niche correlation   ", italic(r))),
        cex.lab = 1.25
        # main = "Interaction relevance space"
        )
  axis(1, at = seq_along(degrees), labels = degrees, las = 2, cex.axis = 1.0)
  y_keep <- unique(round(seq(1, length(correlations), length.out = 5)))
  axis(2, at = y_keep, labels = sprintf("%.2f", correlations[y_keep]), las = 1)
  box()
  for (i in seq_len(nrow(selected))) {
    points(match(selected$degree[i], degrees), match(selected$target_r[i], correlations),
           pch = 21, bg = selected$colour[i], col = selected$colour[i], lwd = 2.8, cex = 2.0)
  }
  mtext("A", side = 3, adj = -0.10, line = 1.2, font = 2, cex = 1.35)

  par(mar = c(5.0, 0.7, 3.5, 2.8))
  image(1, seq(0, 1, length.out = 100), matrix(seq(0, 1, length.out = 100), nrow = 1),
        col = palette, axes = FALSE, xlab = "", ylab = "")
  axis(4, at = seq(0, 1, by = 0.25), labels = sprintf("%.2f", seq(0, 1, by = 0.25)), las = 1)
  mtext("Mismatch", side = 3, line = 0.5, cex = 0.78)
  box()

  par(mar = c(7.0, 5.0, 3.5, 1.0), mgp = c(2.7, 0.8, 0), tcl = -0.25)
  lower <- selected$mean_delta_auc - 1.96 * selected$se_delta_auc
  upper <- selected$mean_delta_auc + 1.96 * selected$se_delta_auc
  ylim <- range(c(0, lower, upper), finite = TRUE)
  padding <- max(diff(ylim) * 0.10, 0.002)
  plot(seq_len(nrow(selected)), selected$mean_delta_auc,
       xlim = c(0.55, 4.45), ylim = c(ylim[1] - padding, ylim[2] + padding),
       xaxt = "n", xlab = "", ylab = "Increase in AUC after adding resource distributions",
       cex.lab = 1.25,
      #  main = "sdm application: benefit of adding resource distributions",
       pch = 21, bg = selected$colour, col = "white", lwd = 2,
       cex = 2.1)
  abline(h = 0, lty = 2, col = "grey55")
  arrows(seq_len(nrow(selected)), lower, seq_len(nrow(selected)), upper,
         angle = 90, code = 3, length = 0.06, lwd = 1.8,
         col = selected$colour)
  points(seq_len(nrow(selected)), selected$mean_delta_auc,
         pch = 21, bg = selected$colour, col = "white", lwd = 2, cex = 2.1)
  degree_labels <- paste("Degree", selected$degree)
  niche_labels <- as.expression(lapply(seq_len(nrow(selected)), function(i) {
    correlation_text <- format(selected$requested_correlation[i], nsmall = 1)
    bquote(paste("niche ", italic(r), " = ", .(correlation_text)))
  }))
  axis(1, at = seq_len(nrow(selected)), labels = degree_labels,
       tick = FALSE, line = 0.1, cex.axis = 1.0)
  axis(1, at = seq_len(nrow(selected)), labels = niche_labels,
       tick = FALSE, line = 1.4, cex.axis = 1.0)
  box()
  mtext("B", side = 3, adj = -0.10, line = 1.2, font = 2, cex = 1.35)
  # mtext("Points are means; bars are 95% intervals",
  #       side = 3, line = 0.2, cex = 0.78)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(output_dir, "figure5_sdm_application.png"), width = 2500, height = 1250, res = 220)
draw_figure()
dev.off()
write.table(selected, file.path(output_dir, "figure5_selected_scenarios.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved Figure 5 and its four-scenario summary to", output_dir, "\n")
