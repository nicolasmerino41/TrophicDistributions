library(tidyverse)
library(viridis)

OUTROOT <- "Outputs"
MEAN_DIR <- file.path(OUTROOT, "meanMetrics")
TAIL_DIR <- file.path(OUTROOT, "tailMetrics")
FIG_DIR  <- file.path(OUTROOT, "heatmaps")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

CONNECTANCE_RANGE <- c(0.005, 0.15) # ATTENTION!!!! YOU MUST CHANGE THIS IF YOU CHANGED THE NUMBER OF CONNECTANCE VALUES
CORR_RANGE        <- c(0.0, 1.0)   # ATTENTION!!!! YOU MUST CHANGE THIS IF YOU CHANGED THE NUMBER OF CORRELATION VALUES

ENV_LEVELS <- c("random", "autocorr")
NET_LEVELS <- c("random", "modular", "heavytail", "cascade")
NET_LABELS <- c(
  random = "Random",
  modular = "Modular",
  heavytail = "Heavytail",
  cascade = "Cascade"
)

REG_LABELS <- c(
  "Narrow + LowVar",
  "Narrow + HighVar",
  "Broad + LowVar",
  "Broad + HighVar"
)

METRIC_LABELS <- c(
  dSrel = "Relative richness loss\n(1 − S[AB] / S[A])",
  mean_jaccard_mismatch = "Mean Jaccard mismatch",
  frac_affected = "Fraction of consumers affected",
  realized_overlap = "Realized prey-support overlap",
  achieved_r = "Achieved niche correlation",
  Creal = "Realized connectance (L / S²)",
  mismatch_q90 = "90th percentile of Jaccard mismatch",
  mismatch_frac_gt = "Fraction with strong mismatch\n(Jaccard mismatch > threshold)"
)

DIFF_LABEL <- "Mean Jaccard mismatch\n(random − autocorr)"

theme_heat <- function() {
  theme_classic(base_family = "Arial", base_size = 13) +
    theme(
      plot.margin = margin(12, 18, 12, 18),
      strip.text = element_text(size = 13),
      strip.background = element_blank(),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.6, "lines"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.title.x = element_text(face = "bold", size = 15, margin = margin(t = 14)),
      axis.title.y = element_text(face = "bold", size = 15, margin = margin(r = 14)),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.title = element_text(face = "bold", size = 13, margin = margin(b = 10)),
      legend.text = element_text(size = 11),
      legend.position = "bottom"
    )
}

fmt_connectance <- function(x) {
  x_round <- round(x / 0.005) * 0.005
  sapply(x_round, function(val) {
    third_decimal <- round(val * 1000) %% 10
    if (third_decimal == 0) sprintf("%.2f", val) else sprintf("%.3f", val)
  })
}

read_matrix_long <- function(file, env, net, reg, metric) {
  M <- as.matrix(read.table(file, sep = "\t"))
  nR <- nrow(M)
  nC <- ncol(M)
  Cvals <- seq(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length.out = nC)
  Rvals <- seq(CORR_RANGE[1], CORR_RANGE[2], length.out = nR)

  df <- expand.grid(r_index = 1:nR, c_index = 1:nC)
  df$value <- as.vector(M)
  df$Connectance <- Cvals[df$c_index]
  df$NicheCorr   <- Rvals[df$r_index]
  df$Environment <- env
  df$Network     <- net
  df$Regime      <- REG_LABELS[reg]
  df$Metric      <- metric
  df
}

load_dir <- function(dirpath) {
  files <- list.files(dirpath, pattern = "^mat_.*\\.tsv$", full.names = TRUE)
  if (length(files) == 0) return(tibble())

  parse_one <- function(f) {
    base <- basename(f)
    m <- str_match(base, "^mat_([^_]+)_([^_]+)_reg([0-9]+)_([^\\.]+)\\.tsv$")
    if (any(is.na(m))) return(NULL)
    env <- m[2]
    net <- m[3]
    reg <- as.integer(m[4])
    metric <- m[5]
    read_matrix_long(f, env, net, reg, metric)
  }

  bind_rows(lapply(files, parse_one))
}

heat_df <- bind_rows(
  load_dir(MEAN_DIR),
  load_dir(TAIL_DIR)
)

if (nrow(heat_df) == 0) {
  stop("No TSV matrices found in Outputs/meanMetrics or Outputs/tailMetrics.")
}

heat_df <- heat_df %>%
  mutate(
    Environment = factor(Environment, levels = ENV_LEVELS),
    Network = factor(Network, levels = NET_LEVELS),
    Regime = factor(Regime, levels = REG_LABELS),
    NetworkLabel = factor(NET_LABELS[as.character(Network)], levels = NET_LABELS[NET_LEVELS])
  )

heat_df <- heat_df %>%
  arrange(Environment, Network, Regime, Metric, NicheCorr, Connectance) %>%
  group_by(Environment, Network, Regime, Metric, NicheCorr) %>%
  mutate(
    value = {
        ok <- is.finite(value)
        if (sum(ok) >= 2) {
            approx(x = Connectance[ok], y = value[ok], xout = Connectance, rule = 2)$y
        } else if (sum(ok) == 1) {
            rep(value[ok][1], length(value))
        } else {
            value
        }
}
  ) %>%
  ungroup()

plot_heatmap_metric <- function(metric_name, env_name) {
  df <- heat_df %>% filter(Metric == metric_name, Environment == env_name)

  ggplot(df, aes(x = Connectance, y = NicheCorr, fill = value)) +
    geom_tile() +
    facet_grid(NetworkLabel ~ Regime) +
    scale_x_continuous(
      breaks = seq(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length.out = 5),
      labels = fmt_connectance,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "viridis",
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "left",
        label.position = "bottom",
        barwidth  = unit(8, "cm"),
        barheight = unit(0.5, "cm"),
        title.theme = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 0.9),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(
      x = "Connectance",
      y = "Niche correlation",
      fill = METRIC_LABELS[[metric_name]]
    ) +
    theme_heat()
}

plot_diff_mean_jaccard <- function() {
  df_diff <- heat_df %>%
    filter(Metric == "mean_jaccard_mismatch", Environment %in% c("random", "autocorr")) %>%
    select(Environment, Network, NetworkLabel, Regime, Connectance, NicheCorr, value) %>%
    pivot_wider(names_from = Environment, values_from = value) %>%
    mutate(value = random - autocorr) %>%
    select(Network, NetworkLabel, Regime, Connectance, NicheCorr, value)

  lim <- max(abs(df_diff$value), na.rm = TRUE)

  ggplot(df_diff, aes(x = Connectance, y = NicheCorr, fill = value)) +
    geom_tile() +
    facet_grid(NetworkLabel ~ Regime) +
    scale_x_continuous(
      breaks = seq(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length.out = 5),
      labels = fmt_connectance,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = c(-lim, lim),
      oob = scales::squish,
      guide = guide_colorbar(
        direction = "horizontal",
        title.position = "left",
        label.position = "bottom",
        barwidth  = unit(8, "cm"),
        barheight = unit(0.5, "cm"),
        title.theme = element_text(face = "bold", size = 13, hjust = 0.5, vjust = 0.9),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    labs(
      x = "Connectance",
      y = "Niche correlation",
      fill = DIFF_LABEL
    ) +
    theme_heat()
}

metrics_mean <- c("dSrel", "mean_jaccard_mismatch", "frac_affected", "realized_overlap", "achieved_r", "Creal")
metrics_tail <- c("mismatch_q90", "mismatch_frac_gt")
metrics_all  <- unique(c(metrics_mean, metrics_tail))
metrics_all  <- metrics_all[metrics_all %in% unique(heat_df$Metric)]

for (env in ENV_LEVELS) {
  for (metric in metrics_all) {
    p <- plot_heatmap_metric(metric, env)
    ggsave(
      filename = file.path(FIG_DIR, paste0("heatmap_", env, "_", metric, ".png")),
      plot = p,
      width = 12,
      height = 10,
      dpi = 600
    )
  }
}

if ("mean_jaccard_mismatch" %in% unique(heat_df$Metric)) {
  p_diff <- plot_diff_mean_jaccard()
  ggsave(
    filename = file.path(FIG_DIR, "Difference_mean_jaccard_mismatch.png"),
    plot = p_diff,
    width = 12,
    height = 10,
    dpi = 600
  )
}

cat("Saved heatmaps to:", FIG_DIR, "\n")
