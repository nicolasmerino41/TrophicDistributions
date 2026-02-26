library(tidyverse)

OUTROOT <- file.path("Outputs", "thermal_metrics")
OUTDIR  <- file.path(OUTROOT, "CombinedFigures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("ctmin", "ctmax", "lt50", "ltmax", "ltmin")
TRAITS <- TRAITS[TRAITS %in% list.dirs(OUTROOT, recursive = FALSE, full.names = FALSE)]

theme_nature <- function() {
  theme_classic(base_family = "Arial", base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    )
}

load_trait <- function(trait) {
  folder <- file.path(OUTROOT, trait)
  edge_file <- file.path(folder, "edge_table.csv")
  node_file <- file.path(folder, "predator_level.csv")
  null_file <- file.path(folder, "null_absdiff_values.csv")
  sum_file  <- file.path(folder, "summary.csv")

  if (!file.exists(edge_file)) return(NULL)

  edges <- read_csv(edge_file, show_col_types = FALSE) %>% mutate(trait = trait)
  nodes <- if (file.exists(node_file)) read_csv(node_file, show_col_types = FALSE) %>% mutate(trait = trait) else NULL
  nulls <- if (file.exists(null_file)) read_csv(null_file, show_col_types = FALSE) %>% mutate(trait = trait) else NULL
  summ  <- if (file.exists(sum_file))  read_csv(sum_file,  show_col_types = FALSE) %>% mutate(trait = trait) else NULL

  list(edges = edges, nodes = nodes, nulls = nulls, summ = summ)
}

all_edges <- list()
all_nodes <- list()
all_nulls <- list()
all_summ  <- list()

for (tr in TRAITS) {
  dat <- load_trait(tr)
  if (!is.null(dat)) {
    all_edges[[tr]] <- dat$edges
    if (!is.null(dat$nodes)) all_nodes[[tr]] <- dat$nodes
    if (!is.null(dat$nulls)) all_nulls[[tr]] <- dat$nulls
    if (!is.null(dat$summ))  all_summ[[tr]]  <- dat$summ
  }
}

edges_df <- bind_rows(all_edges)
nodes_df <- bind_rows(all_nodes)
null_df  <- bind_rows(all_nulls)
summ_df  <- bind_rows(all_summ)

edges_df$trait <- factor(edges_df$trait, levels = TRAITS)
nodes_df$trait <- factor(nodes_df$trait, levels = TRAITS)
null_df$trait  <- factor(null_df$trait,  levels = TRAITS)
summ_df$trait  <- factor(summ_df$trait,  levels = TRAITS)

compute_stats <- function(df, xvar, yvar) {
  df %>%
    filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]])) %>%
    group_by(trait) %>%
    summarise(
      n = n(),
      r = cor(.data[[xvar]], .data[[yvar]]),
      label = paste0("N = ", n, "\nr = ", round(r, 2)),
      .groups = "drop"
    )
}

edge_stats <- compute_stats(edges_df, "pred_trait", "prey_trait")
node_stats <- compute_stats(nodes_df, "pred_trait", "mean_prey_trait")

p_edge <- ggplot(edges_df, aes(x = pred_trait, y = prey_trait)) +
  geom_point(size = 2, alpha = 0.9, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
  geom_text(
    data = edge_stats, aes(label = label),
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE, size = 5
  ) +
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  labs(
    title = "Pairwise Predator vs Prey Thermal Limits",
    x = "Predator thermal limit (°C)",
    y = "Prey thermal limit (°C)"
  ) +
  theme_nature()

ggsave(file.path(OUTDIR, "EDGE_2x2_combined.png"), p_edge, width = 10, height = 8, dpi = 600)

p_node <- ggplot(nodes_df, aes(x = pred_trait, y = mean_prey_trait)) +
  geom_point(size = 2.5, alpha = 0.95, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") +
  geom_text(
    data = node_stats, aes(label = label),
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE, size = 5
  ) +
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  labs(
    title = "Predator vs Mean Prey Thermal Limits",
    x = "Predator thermal limit (°C)",
    y = "Mean prey thermal limit (°C)"
  ) +
  theme_nature()

ggsave(file.path(OUTDIR, "NODE_2x2_combined.png"), p_node, width = 10, height = 8, dpi = 600)

if (nrow(null_df) > 0 && nrow(summ_df) > 0) {
  summ_null <- summ_df %>% transmute(trait, observed = obs_mean_absdiff)

  p_null <- ggplot(null_df, aes(x = null)) +
    geom_histogram(bins = 30, fill = "grey70", color = "black") +
    geom_vline(data = summ_null, aes(xintercept = observed), color = "red", linewidth = 1.2) +
    facet_wrap(~ trait, ncol = 2, scales = "free") +
    labs(
      title = "Null model: mean |Predator − Prey| thermal difference",
      x = "Mean |Δ thermal limit|",
      y = "Count"
    ) +
    theme_nature()

  ggsave(file.path(OUTDIR, "EDGE_null_mean_absdiff_2x2.png"), p_null, width = 10, height = 8, dpi = 600)
}

cat("Saved figures to:", OUTDIR, "\n")
