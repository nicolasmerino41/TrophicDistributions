read_sdm_table <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

degree_order <- c(1, 2, 3, 4, 5, 6, 8, 10, 13, 17, 22, 28, 35, 44, 55)

make_matrix <- function(data, value, degrees = degree_order) {
  correlations <- sort(unique(data$target_r))
  result <- matrix(NA_real_, nrow = length(correlations), ncol = length(degrees),
                   dimnames = list(format(correlations, trim = TRUE), degrees))
  for (i in seq_len(nrow(data))) {
    row_index <- match(data$target_r[i], correlations)
    column_index <- match(data$degree[i], degrees)
    if (!is.na(row_index) && !is.na(column_index)) {
      result[row_index, column_index] <- data[[value]][i]
    }
  }
  result
}

draw_heatmap <- function(data, value, title, limits, palette,
                         show_x = TRUE, show_y = TRUE) {
  values <- make_matrix(data, value)
  breaks <- seq(limits[1], limits[2], length.out = length(palette) + 1)
  image(seq_len(ncol(values)), seq_len(nrow(values)), t(values),
        col = palette, breaks = breaks, axes = FALSE,
        xlab = if (show_x) "Focal consumer degree" else "",
        ylab = if (show_y) "Community niche correlation" else "",
        main = title, useRaster = TRUE)
  if (show_x) axis(1, at = seq_len(ncol(values)), labels = colnames(values),
                   las = 2, cex.axis = 0.68)
  if (show_y) {
    labels <- as.numeric(rownames(values))
    keep <- unique(round(seq(1, length(labels), length.out = 5)))
    axis(2, at = keep, labels = sprintf("%.2f", labels[keep]), las = 1)
  }
  box()
}

stage_label <- function(stage) {
  labels <- c(
    oracle = "True prey distributions",
    estimated_100 = "Estimated prey; 100% links",
    estimated_75 = "Estimated prey; 75% links",
    estimated_50 = "Estimated prey; 50% links",
    estimated_25 = "Estimated prey; 25% links"
  )
  unname(labels[stage])
}

draw_relationship <- function(data) {
  stages <- c("oracle", "estimated_100", "estimated_75", "estimated_50", "estimated_25")
  colours <- c("#1b1b1b", "#0072B2", "#009E73", "#E69F00", "#D55E00")
  plot(NA, xlim = c(0, 1), ylim = range(c(0, data$mean_delta_auc), na.rm = TRUE),
       xlab = "True A-AB mismatch", ylab = expression(Delta*"AUC"),
       main = "SDM gain follows interaction relevance")
  abline(h = 0, col = "grey70", lty = 2)
  for (i in seq_along(stages)) {
    subset <- data[data$information_level == stages[i] & data$n >= 2, ]
    subset <- subset[order(subset$mean_true_mismatch), ]
    if (!nrow(subset)) next
    lower <- subset$mean_delta_auc - 1.96 * subset$se_delta_auc
    upper <- subset$mean_delta_auc + 1.96 * subset$se_delta_auc
    valid <- is.finite(lower) & is.finite(upper)
    if (sum(valid) >= 2) {
      polygon(c(subset$mean_true_mismatch[valid], rev(subset$mean_true_mismatch[valid])),
              c(lower[valid], rev(upper[valid])),
              col = adjustcolor(colours[i], alpha.f = 0.13), border = NA)
    }
    lines(subset$mean_true_mismatch, subset$mean_delta_auc,
          col = colours[i], lwd = 2)
  }
  legend("topleft", legend = stage_label(stages), col = colours,
         lwd = 2, bty = "n", cex = 0.75)
}

add_panel_label <- function(label) {
  mtext(label, side = 3, adj = -0.10, line = 1.1, font = 2, cex = 1.15)
}
