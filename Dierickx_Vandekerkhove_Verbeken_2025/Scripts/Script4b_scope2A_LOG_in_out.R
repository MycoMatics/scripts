# =============================================================================
# Script: Script4b_scope2A_LOG_in_out.R
# Project: Deadwood fungal community analyses
# Purpose: Paired INNER–OUTER comparison within LOG cores using NMDS,
#          segment-length summaries, paired overlap/subsetness metrics,
#          alpha diversity summaries, multivariate dispersion, and paired
#          PERMANOVA blocking by natman × position
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - pick_col
# - ensure_numeric_matrix
# - fix_dist_labels
# - axis_contrib
# - subset_dist
# - shannon1
# - make_venn_plot
# - fmt_mean_sd
# - run_seglen_tests
# - ensure_dir
# - assert_objects
#
# Outputs:
# - plots/INOUT_NMDS/*
# - tables/INOUT_NMDS/*
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(glue)
  library(janitor)
  library(cowplot)
  library(readr)
  library(RColorBrewer)
  library(betapart)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c(
  "pick_col", "ensure_numeric_matrix", "fix_dist_labels", "axis_contrib",
  "subset_dist", "shannon1", "make_venn_plot", "fmt_mean_sd",
  "run_seglen_tests", "ensure_dir"
))

set.seed(SEED_GLOBAL)

# ---- Paths -------------------------------------------------------------------
data_dir <- "."
OUTDIR_PLOT <- "plots/INOUT_NMDS"
OUTDIR_TAB  <- "tables/INOUT_NMDS"

ensure_dir(data_dir)
setwd(data_dir)
ensure_dir(OUTDIR_PLOT)
ensure_dir(OUTDIR_TAB)

# ---- Parameters --------------------------------------------------------------
NMD_K_AIT <- 2
NMD_K_JAC <- 3
TRYMAX    <- 9999
N_PERM    <- 9999
DROP_SAMPLES <- c("P1_56", "P2_84", "P3_42", "P3_41", "P3_37")

shape_map_depth <- c("INNER" = 24, "OUTER" = 25)
shape_map_pos   <- c("BASE" = 15, "MIDDLE" = 16, "UPPER" = 18)

if (!exists("ds_colors")) {
  ds_colors <- setNames(RColorBrewer::brewer.pal(6, "YlOrRd"), as.character(0:5))
}

# ---- Inputs ------------------------------------------------------------------
meta <- META1 %>%
  as_tibble() %>%
  filter(!sample %in% DROP_SAMPLES) %>%
  clean_names()

Xall <- otu_matrix_filt
stopifnot(all(meta$sample %in% rownames(Xall)))

col_dw    <- pick_col(meta, c("dw_type2"))
col_pos   <- pick_col(meta, c("position", "position_2", "position2"))
col_depth <- pick_col(meta, c("depth_2"))
col_ds    <- pick_col(meta, c("ds_at_drill"))
col_nat   <- pick_col(meta, c("natman"))

meta0 <- meta %>%
  mutate(
    dw_type  = .data[[col_dw]],
    position = .data[[col_pos]],
    depth_2  = .data[[col_depth]],
    ds_raw   = .data[[col_ds]],
    natman   = .data[[col_nat]]
  ) %>%
  filter(
    dw_type == "LOG",
    depth_2 %in% c("INNER", "OUTER"),
    !is.na(natman),
    !is.na(position)
  ) %>%
  transmute(
    sample,
    natman,
    position,
    depth_2,
    ds_raw
  ) %>%
  mutate(
    ds_fac = droplevels(factor(ds_raw)),
    ds_num = suppressWarnings(as.numeric(as.character(ds_raw)))
  )

# ---- Enforce natman × position pairs -----------------------------------------
pair_tbl <- meta0 %>%
  mutate(pair_key = glue("{natman}|{position}")) %>%
  count(pair_key, depth_2, name = "n_depth") %>%
  tidyr::pivot_wider(names_from = depth_2, values_from = n_depth, values_fill = 0)

keep_pairs <- pair_tbl %>%
  filter(INNER == 1, OUTER == 1) %>%
  pull(pair_key)

meta1 <- meta0 %>%
  mutate(pair_key = glue("{natman}|{position}")) %>%
  filter(pair_key %in% keep_pairs) %>%
  mutate(agg_key = glue("{pair_key}|{depth_2}")) %>%
  arrange(pair_key, depth_2)

# ---- Build X + design --------------------------------------------------------
X <- Xall[meta1$sample, , drop = FALSE]
X <- ensure_numeric_matrix(X)
rownames(X) <- meta1$agg_key
X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]

design <- meta1 %>%
  transmute(
    agg_key  = as.character(agg_key),
    pair_key = as.character(pair_key),
    natman   = natman,
    position = position,
    depth    = factor(depth_2, levels = c("INNER", "OUTER")),
    ds_fac   = ds_fac,
    ds_num   = ds_num
  ) %>%
  arrange(match(agg_key, rownames(X)))

stopifnot(identical(rownames(X), design$agg_key))

cat("\nPaired LOG cores retained:\n")
cat(glue(" - N pairs: {length(unique(design$pair_key))}\n"))
cat(glue(" - N samples (rows): {nrow(design)} (should be 2 × N pairs)\n"))
cat(glue(" - N taxa (cols): {ncol(X)}\n"))

# ---- Distances ----------------------------------------------------------------
d_ait <- vegan::vegdist(X, method = "robust.aitchison")
d_ait <- fix_dist_labels(d_ait, samps = rownames(X))

d_jac <- vegan::vegdist((X > 0) * 1, method = "jaccard", binary = TRUE)
d_jac <- fix_dist_labels(d_jac, samps = rownames(X))

# ==============================================================================
# 1) NMDS + segment lengths: Aitchison
# ==============================================================================
set.seed(SEED_GLOBAL)
ord_ait <- vegan::metaMDS(
  d_ait,
  k = NMD_K_AIT,
  trymax = TRYMAX,
  autotransform = FALSE,
  trace = FALSE
)

cat(glue("\n[Aitchison] NMDS stress: {ord_ait$stress}\n"))
ax_ait <- axis_contrib(ord_ait, d_ait)
print(ax_ait)

sites_ait <- as.data.frame(vegan::scores(ord_ait, display = "sites")) %>%
  tibble::rownames_to_column("agg_key") %>%
  dplyr::rename_with(~ gsub("^MDS", "NMDS", .x), dplyr::starts_with("MDS")) %>%
  dplyr::left_join(design, by = "agg_key")

seg_ait <- sites_ait %>%
  dplyr::select(pair_key, depth, ds_fac, dplyr::starts_with("NMDS")) %>%
  dplyr::group_by(pair_key) %>%
  dplyr::filter(dplyr::n_distinct(depth) == 2) %>%
  dplyr::summarise(
    NMDS1_inner = NMDS1[depth == "INNER"],
    NMDS2_inner = NMDS2[depth == "INNER"],
    NMDS3_inner = if ("NMDS3" %in% names(dplyr::cur_data())) NMDS3[depth == "INNER"] else NA_real_,
    NMDS1_outer = NMDS1[depth == "OUTER"],
    NMDS2_outer = NMDS2[depth == "OUTER"],
    NMDS3_outer = if ("NMDS3" %in% names(dplyr::cur_data())) NMDS3[depth == "OUTER"] else NA_real_,
    ds_fac = ds_fac[1],
    seg_len = if (NMD_K_AIT >= 3) {
      sqrt((NMDS1_outer - NMDS1_inner)^2 +
             (NMDS2_outer - NMDS2_inner)^2 +
             (NMDS3_outer - NMDS3_inner)^2)
    } else {
      sqrt((NMDS1_outer - NMDS1_inner)^2 +
             (NMDS2_outer - NMDS2_inner)^2)
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(ds_fac = droplevels(factor(ds_fac)))

cat(glue("\n[Aitchison] Segment-length summary by decay:\n"))
print(
  seg_ait %>%
    dplyr::group_by(ds_fac) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(seg_len),
      IQR = IQR(seg_len),
      .groups = "drop"
    )
)

p_seg_ait <- ggplot2::ggplot(seg_ait, ggplot2::aes(x = ds_fac, y = seg_len, fill = ds_fac)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.2) +
  ggplot2::scale_fill_manual(values = ds_colors, guide = "none") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(p_seg_ait)
ggplot2::ggsave(
  filename = file.path(OUTDIR_PLOT, "rAIT_seglengths_quarter_A4.png"),
  plot = p_seg_ait,
  width = 68.5,
  height = 33.5,
  units = "mm",
  dpi = 900
)

xlab_ait <- paste0("NMDS1 (", round(100 * ax_ait$frac[ax_ait$axis == "NMDS1"], 1), "%)")
ylab_ait <- paste0("NMDS2 (", round(100 * ax_ait$frac[ax_ait$axis == "NMDS2"], 1), "%)")

sites_ait_inner <- sites_ait %>% dplyr::filter(depth != "OUTER")
sites_ait_outer <- sites_ait %>% dplyr::filter(depth == "OUTER")

p_nmds_ait <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = seg_ait,
    ggplot2::aes(x = NMDS1_inner, y = NMDS2_inner, xend = NMDS1_outer, yend = NMDS2_outer),
    colour = "grey70",
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    data = sites_ait_inner,
    ggplot2::aes(x = NMDS1, y = NMDS2, fill = ds_fac, shape = depth),
    colour = "black",
    stroke = 0.2,
    size = 2.0,
    alpha = 0.9
  ) +
  ggplot2::geom_point(
    data = sites_ait_outer,
    ggplot2::aes(x = NMDS1, y = NMDS2, fill = ds_fac, shape = depth),
    colour = "black",
    stroke = 0.75,
    size = 1.6,
    alpha = 0.9
  ) +
  ggplot2::scale_fill_manual(values = ds_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = shape_map_depth, guide = "none") +
  ggplot2::labs(x = xlab_ait, y = ylab_ait) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "none")

print(p_nmds_ait)
ggplot2::ggsave(
  filename = file.path(OUTDIR_PLOT, "Aitchison_NMDS_quarter_A4.png"),
  plot = p_nmds_ait,
  width = 68.5,
  height = 62.5,
  units = "mm",
  dpi = 900
)

# ==============================================================================
# 2) NMDS + segment lengths: Jaccard
# ==============================================================================
set.seed(SEED_GLOBAL)
ord_jac <- vegan::metaMDS(
  d_jac,
  k = NMD_K_JAC,
  trymax = TRYMAX,
  autotransform = FALSE,
  trace = FALSE
)

cat(glue("\n[Jaccard] NMDS stress: {ord_jac$stress}\n"))
ax_jac <- axis_contrib(ord_jac, d_jac)
print(ax_jac)

sites_jac <- as.data.frame(vegan::scores(ord_jac, display = "sites")) %>%
  tibble::rownames_to_column("agg_key") %>%
  dplyr::rename_with(~ gsub("^MDS", "NMDS", .x), dplyr::starts_with("MDS")) %>%
  dplyr::left_join(design, by = "agg_key")

seg_jac <- sites_jac %>%
  dplyr::select(pair_key, depth, ds_fac, dplyr::starts_with("NMDS")) %>%
  dplyr::group_by(pair_key) %>%
  dplyr::filter(dplyr::n_distinct(depth) == 2) %>%
  dplyr::summarise(
    NMDS1_inner = NMDS1[depth == "INNER"],
    NMDS2_inner = NMDS2[depth == "INNER"],
    NMDS3_inner = if ("NMDS3" %in% names(dplyr::cur_data())) NMDS3[depth == "INNER"] else NA_real_,
    NMDS1_outer = NMDS1[depth == "OUTER"],
    NMDS2_outer = NMDS2[depth == "OUTER"],
    NMDS3_outer = if ("NMDS3" %in% names(dplyr::cur_data())) NMDS3[depth == "OUTER"] else NA_real_,
    ds_fac = ds_fac[1],
    seg_len = if (NMD_K_JAC >= 3) {
      sqrt((NMDS1_outer - NMDS1_inner)^2 +
             (NMDS2_outer - NMDS2_inner)^2 +
             (NMDS3_outer - NMDS3_inner)^2)
    } else {
      sqrt((NMDS1_outer - NMDS1_inner)^2 +
             (NMDS2_outer - NMDS2_inner)^2)
    },
    .groups = "drop"
  ) %>%
  dplyr::mutate(ds_fac = droplevels(factor(ds_fac)))

cat(glue("\n[Jaccard] Segment-length summary by decay:\n"))
print(
  seg_jac %>%
    dplyr::group_by(ds_fac) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(seg_len),
      IQR = IQR(seg_len),
      .groups = "drop"
    )
)

p_seg_jac <- ggplot2::ggplot(seg_jac, ggplot2::aes(x = ds_fac, y = seg_len, fill = ds_fac)) +
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.2) +
  ggplot2::scale_fill_manual(values = ds_colors, guide = "none") +
  ggplot2::labs(x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

print(p_seg_jac)
ggplot2::ggsave(
  filename = file.path(OUTDIR_PLOT, "Jaccard_seglengths_quarter_A4.png"),
  plot = p_seg_jac,
  width = 68.5,
  height = 33.5,
  units = "mm",
  dpi = 900
)

xlab_jac <- paste0("NMDS1 (", round(100 * ax_jac$frac[ax_jac$axis == "NMDS1"], 1), "%)")
ylab_jac <- paste0("NMDS2 (", round(100 * ax_jac$frac[ax_jac$axis == "NMDS2"], 1), "%)")

sites_jac_inner <- sites_jac %>% dplyr::filter(depth != "OUTER")
sites_jac_outer <- sites_jac %>% dplyr::filter(depth == "OUTER")

p_nmds_jac <- ggplot2::ggplot() +
  ggplot2::geom_segment(
    data = seg_jac,
    ggplot2::aes(x = NMDS1_inner, y = NMDS2_inner, xend = NMDS1_outer, yend = NMDS2_outer),
    colour = "grey70",
    alpha = 0.7
  ) +
  ggplot2::geom_point(
    data = sites_jac_inner,
    ggplot2::aes(x = NMDS1, y = NMDS2, fill = ds_fac, shape = depth),
    colour = "black",
    stroke = 0.2,
    size = 2.0,
    alpha = 0.9
  ) +
  ggplot2::geom_point(
    data = sites_jac_outer,
    ggplot2::aes(x = NMDS1, y = NMDS2, fill = ds_fac, shape = depth),
    colour = "black",
    stroke = 0.75,
    size = 1.6,
    alpha = 0.9
  ) +
  ggplot2::scale_fill_manual(values = ds_colors, guide = "none") +
  ggplot2::scale_shape_manual(values = shape_map_depth, guide = "none") +
  ggplot2::labs(x = xlab_jac, y = ylab_jac) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "none")

print(p_nmds_jac)
ggplot2::ggsave(
  filename = file.path(OUTDIR_PLOT, "Jaccard_NMDS_quarter_A4.png"),
  plot = p_nmds_jac,
  width = 68.5,
  height = 62.5,
  units = "mm",
  dpi = 900
)

# ---- Segment-length tests ----------------------------------------------------
res_seg_ait <- run_seglen_tests(seg_ait, label = "Aitchison", decay_col = "ds_fac", seglen_col = "seg_len")
res_seg_jac <- run_seglen_tests(seg_jac, label = "Jaccard", decay_col = "ds_fac", seglen_col = "seg_len")

# ==============================================================================
# 3) Pair-level reconstruction: overlap, subsetness, and decay-gradient tests
# ==============================================================================
res_pairs <- design %>%
  mutate(
    pair_key   = as.character(pair_key),
    agg_key    = as.character(agg_key),
    ds_fac     = droplevels(factor(ds_fac)),
    ds_numeric = suppressWarnings(as.numeric(as.character(ds_fac)))
  ) %>%
  group_by(pair_key) %>%
  filter(n_distinct(depth) == 2) %>%
  summarise(
    natman   = first(natman),
    position = first(position),
    ds_fac   = first(ds_fac),
    ds_numeric = first(ds_numeric),
    nI = sum(X[agg_key[depth == "INNER"], ] > 0),
    nO = sum(X[agg_key[depth == "OUTER"], ] > 0),
    n_shared = sum((X[agg_key[depth == "INNER"], ] > 0) & (X[agg_key[depth == "OUTER"], ] > 0)),
    n_union  = sum((X[agg_key[depth == "INNER"], ] > 0) | (X[agg_key[depth == "OUTER"], ] > 0)),
    frac_inner_in_outer = ifelse(nI > 0, n_shared / nI, NA_real_),
    frac_outer_in_inner = ifelse(nO > 0, n_shared / nO, NA_real_),
    subset_score = pmin(frac_inner_in_outer, frac_outer_in_inner, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ds_factor = droplevels(as.factor(floor(ds_numeric))),
    ds_band = case_when(
      is.finite(ds_numeric) & ds_numeric <= 1 ~ "Early (1)",
      is.finite(ds_numeric) & ds_numeric <= 3 ~ "Intermediate (2-3)",
      is.finite(ds_numeric) & ds_numeric >= 4 ~ "Late (4–5)",
      TRUE ~ "Unknown"
    ),
    ds_band = factor(ds_band, levels = c("Early (1)", "Intermediate (2-3)", "Late (4–5)", "Unknown"))
  )

cat("\n[Pairs] N pairs:", nrow(res_pairs), "\n")

# ---- 3.1 Pairwise INNER–OUTER distances vs decay -----------------------------
labs_d <- attr(d_ait, "Labels")
if (is.null(labs_d)) labs_d <- rownames(as.matrix(d_ait))

M_ait <- as.matrix(d_ait)
M_jac <- as.matrix(d_jac)

pair_dist <- res_pairs %>%
  filter(is.finite(ds_numeric)) %>%
  mutate(
    inner_lab = paste0(pair_key, "|INNER"),
    outer_lab = paste0(pair_key, "|OUTER"),
    idx_inner = match(inner_lab, labs_d),
    idx_outer = match(outer_lab, labs_d)
  ) %>%
  mutate(
    dist_ait = purrr::map2_dbl(idx_inner, idx_outer, ~ if (any(is.na(c(.x, .y)))) NA_real_ else M_ait[.x, .y]),
    dist_jac = purrr::map2_dbl(idx_inner, idx_outer, ~ if (any(is.na(c(.x, .y)))) NA_real_ else M_jac[.x, .y])
  )

pair_dist_summary <- pair_dist %>%
  group_by(ds_numeric) %>%
  summarise(
    n_pairs = n(),
    median_ait = median(dist_ait, na.rm = TRUE),
    IQR_ait = IQR(dist_ait, na.rm = TRUE),
    median_jac = median(dist_jac, na.rm = TRUE),
    IQR_jac = IQR(dist_jac, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n3.1) INNER–OUTER distances by decay (numeric):\n")
print(pair_dist_summary)
readr::write_csv2(pair_dist_summary, file.path(OUTDIR_TAB, "pair_dist_by_decay_numeric.csv"))

cor_ait_decay <- suppressWarnings(cor.test(pair_dist$ds_numeric, pair_dist$dist_ait, method = "spearman", alternative = "greater"))
cor_jac_decay <- suppressWarnings(cor.test(pair_dist$ds_numeric, pair_dist$dist_jac, method = "spearman", alternative = "greater"))

cat("\nSpearman: Aitchison pair distance vs decay (greater):\n")
print(cor_ait_decay)
cat("\nSpearman: Jaccard pair distance vs decay (greater):\n")
print(cor_jac_decay)

kw_ait <- kruskal.test(dist_ait ~ ds_factor, data = pair_dist)
kw_jac <- kruskal.test(dist_jac ~ ds_factor, data = pair_dist)

cat("\nKruskal–Wallis: Aitchison pair distance ~ decay factor:\n")
print(kw_ait)
cat("\nKruskal–Wallis: Jaccard pair distance ~ decay factor:\n")
print(kw_jac)

bf_df <- pair_dist %>%
  group_by(ds_factor) %>%
  mutate(
    dev_ait = abs(dist_ait - median(dist_ait, na.rm = TRUE)),
    dev_jac = abs(dist_jac - median(dist_jac, na.rm = TRUE))
  ) %>%
  ungroup()

bf_ait <- anova(lm(dev_ait ~ ds_factor, data = bf_df))
bf_jac <- anova(lm(dev_jac ~ ds_factor, data = bf_df))

cat("\nBrown–Forsythe-style (median deviations): Aitchison ~ decay factor:\n")
print(bf_ait)
cat("\nBrown–Forsythe-style (median deviations): Jaccard ~ decay factor:\n")
print(bf_jac)

print(pairwise.t.test(bf_df$dev_jac, bf_df$ds_factor, p.adjust.method = "BH"))

# ---- 3.2 Multivariate dispersion within INNER and OUTER ----------------------
for (dep in c("INNER", "OUTER")) {
  keep <- design$depth == dep & is.finite(design$ds_num)
  if (sum(keep) < 6 || dplyr::n_distinct(design$ds_fac[keep]) < 2) next
  
  labs_dep   <- design$agg_key[keep]
  d_dep_ait  <- subset_dist(d_ait, labs_dep)
  d_dep_jac  <- subset_dist(d_jac, labs_dep)
  groups_dep <- droplevels(factor(design$ds_fac[keep]))
  
  disp_ait <- betadisper(d_dep_ait, groups_dep)
  cat("\n3.2) Aitchison betadisper for depth =", dep, "by decay stage:\n")
  print(anova(disp_ait))
  cat("\nPermutation test (Aitchison):\n")
  print(permutest(disp_ait, permutations = N_PERM))
  
  disp_jac <- betadisper(d_dep_jac, groups_dep)
  cat("\n3.2) Jaccard betadisper for depth =", dep, "by decay stage:\n")
  print(anova(disp_jac))
  cat("\nPermutation test (Jaccard):\n")
  print(permutest(disp_jac, permutations = N_PERM))
}

# ---- 3.3 Alpha diversity vs decay and depth ----------------------------------
alpha_df <- tibble(
  agg_key  = rownames(X),
  richness = rowSums(X > 0),
  shannon  = vegan::diversity(X, index = "shannon"),
  simpson  = vegan::diversity(X, index = "simpson")
) %>%
  left_join(design %>% select(agg_key, depth, ds_num, ds_fac), by = "agg_key") %>%
  filter(is.finite(ds_num)) %>%
  rename(ds_numeric = ds_num)

for (dep in c("INNER", "OUTER")) {
  df_dep <- alpha_df %>% filter(depth == dep)
  if (nrow(df_dep) < 6) next
  cat("\n3.3) Spearman alpha vs decay for depth =", dep, "\n")
  cat("Richness vs decay:\n"); print(cor.test(df_dep$ds_numeric, df_dep$richness, method = "spearman"))
  cat("Shannon vs decay:\n"); print(cor.test(df_dep$ds_numeric, df_dep$shannon, method = "spearman"))
  cat("Simpson vs decay:\n"); print(cor.test(df_dep$ds_numeric, df_dep$simpson, method = "spearman"))
}

alpha_df$depth <- droplevels(factor(alpha_df$depth))

lm_rich <- lm(richness ~ ds_numeric * depth, data = alpha_df)
lm_shan <- lm(shannon  ~ ds_numeric * depth, data = alpha_df)
lm_simp <- lm(simpson  ~ ds_numeric * depth, data = alpha_df)

cat("\n3.3) Linear models with decay × depth interaction\n")
cat("\nRichness ~ ds_numeric * depth\n"); print(anova(lm_rich))
cat("\nShannon ~ ds_numeric * depth\n"); print(anova(lm_shan))
cat("\nSimpson ~ ds_numeric * depth\n"); print(anova(lm_simp))

readr::write_csv2(alpha_df, file.path(OUTDIR_TAB, "alpha_df_paired_LOG.csv"))

# ---- 3.4 Symmetric subsetness vs decay ---------------------------------------
df_sub_sym <- res_pairs %>%
  filter(is.finite(ds_numeric), is.finite(subset_score))

sub_sym_summary <- df_sub_sym %>%
  group_by(ds_fac) %>%
  summarise(
    n_pairs = n(),
    median_subset = median(subset_score, na.rm = TRUE),
    IQR_subset = IQR(subset_score, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n3.4) Symmetric subsetness by decay band:\n")
print(sub_sym_summary)

readr::write_csv2(sub_sym_summary, file.path(OUTDIR_TAB, "subset_symmetric_by_decay_band.csv"))

cat("\nSpearman: symmetric subsetness vs decay (numeric):\n")
print(cor.test(df_sub_sym$ds_numeric, df_sub_sym$subset_score, method = "spearman"))

cat("\nKruskal–Wallis: symmetric subsetness ~ decay factor:\n")
print(kruskal.test(subset_score ~ ds_factor, data = df_sub_sym))

# ==============================================================================
# 4) Overlap, subsetness scatter, pooled venns, and per-pair venns
# ==============================================================================
subset_tol <- 0.66
venn_dir <- file.path(OUTDIR_PLOT, "venn")
ensure_dir(venn_dir)

out_tab      <- file.path(OUTDIR_TAB,  "inout_pair_overlap.tsv")
out_plot_sum <- file.path(OUTDIR_PLOT, "inout_subsetness_scatter.png")
out_venn_all <- file.path(venn_dir,    "venn_all_pairs_pooled.png")
out_venn_ms  <- file.path(venn_dir,    "venn_mean_sd.png")

pair_has <- design %>%
  distinct(pair_key) %>%
  pull(pair_key) %>%
  sort()

res_overlap <- purrr::map_dfr(pair_has, function(k) {
  ki <- paste0(k, "|INNER")
  ko <- paste0(k, "|OUTER")
  
  xi <- as.numeric(X[ki, , drop = TRUE])
  xo <- as.numeric(X[ko, , drop = TRUE])
  
  Ai <- xi > 0
  Ao <- xo > 0
  nI <- sum(Ai)
  nO <- sum(Ao)
  n_shared <- sum(Ai & Ao)
  n_union  <- sum(Ai | Ao)
  
  jaccard  <- if (n_union > 0) n_shared / n_union else NA_real_
  sorensen <- if ((nI + nO) > 0) (2 * n_shared) / (nI + nO) else NA_real_
  
  frac_IinO <- if (nI > 0) n_shared / nI else NA_real_
  frac_OinI <- if (nO > 0) n_shared / nO else NA_real_
  
  beta_SIM <- beta_SNE <- beta_SOR <- NA_real_
  beta_BRAY_BAL <- beta_BRAY_GRA <- beta_BRAY <- NA_real_
  
  Mpa <- rbind(INNER = as.integer(Ai), OUTER = as.integer(Ao))
  Mab <- rbind(INNER = xi, OUTER = xo)
  
  bp_pa <- betapart::beta.pair(Mpa, index.family = "sorensen")
  beta_SIM <- as.numeric(bp_pa$beta.sim)
  beta_SNE <- as.numeric(bp_pa$beta.sne)
  beta_SOR <- as.numeric(bp_pa$beta.sor)
  
  bp_ab <- betapart::beta.pair.abund(Mab, index.family = "bray")
  beta_BRAY_BAL <- as.numeric(bp_ab$beta.bray.bal)
  beta_BRAY_GRA <- as.numeric(bp_ab$beta.bray.gra)
  beta_BRAY     <- as.numeric(bp_ab$beta.bray)
  
  nm <- strsplit(k, "\\|")[[1]][1]
  ps <- strsplit(k, "\\|")[[1]][2]
  
  tibble(
    natman = nm,
    position = ps,
    pair_key = k,
    richness_inner = nI,
    richness_outer = nO,
    shannon_inner  = shannon1(xi),
    shannon_outer  = shannon1(xo),
    shared = n_shared,
    union = n_union,
    jaccard = jaccard,
    sorensen = sorensen,
    frac_inner_in_outer = frac_IinO,
    frac_outer_in_inner = frac_OinI,
    is_subset_I_in_O = isTRUE(all(Ai <= Ao)),
    is_subset_O_in_I = isTRUE(all(Ao <= Ai)),
    almost_I_in_O = !is.na(frac_IinO) && frac_IinO >= subset_tol,
    almost_O_in_I = !is.na(frac_OinI) && frac_OinI >= subset_tol,
    inner_only = nI - n_shared,
    outer_only = nO - n_shared,
    richness_diff_OmI = nO - nI,
    beta_SIM = beta_SIM,
    beta_SNE = beta_SNE,
    beta_SOR = beta_SOR,
    beta_BRAY_BAL = beta_BRAY_BAL,
    beta_BRAY_GRA = beta_BRAY_GRA,
    beta_BRAY = beta_BRAY
  )
}) %>%
  arrange(natman, position)

readr::write_tsv(res_overlap, out_tab)

rows_inner <- paste0(pair_has, "|INNER")
rows_outer <- paste0(pair_has, "|OUTER")

Ain_any  <- colSums(X[rows_inner, , drop = FALSE] > 0) > 0
Aout_any <- colSums(X[rows_outer, , drop = FALSE] > 0) > 0

left_only_all  <- sum(Ain_any  & !Aout_any)
overlap_all    <- sum(Ain_any  &  Aout_any)
right_only_all <- sum(Aout_any & !Ain_any)

p_all <- make_venn_plot(
  left_only  = left_only_all,
  overlap    = overlap_all,
  right_only = right_only_all,
  labels     = c("\n\n\n\n\n\nINNER\n(pooled)", "\n\n\n\n\n\nOUTER\n(pooled)"),
  subtitle   = glue("All pairs pooled (N pairs = {length(pair_has)})")
)

print(p_all)
ggsave(out_venn_all, p_all, width = 5.2, height = 4.4, dpi = 300)

p_ms <- make_venn_plot(
  left_only  = fmt_mean_sd(res_overlap$inner_only),
  overlap    = fmt_mean_sd(res_overlap$shared),
  right_only = fmt_mean_sd(res_overlap$outer_only),
  labels     = c("INNER", "OUTER"),
  subtitle   = glue("Per pair mean ± SD (N pairs = {nrow(res_overlap)})")
)

print(p_ms)
ggsave(out_venn_ms, p_ms, width = 5.2, height = 4.4, dpi = 300)

pair_cov <- design %>%
  group_by(pair_key) %>%
  summarise(
    natman = first(natman),
    position = first(position),
    ds_numeric = suppressWarnings(mean(ds_num, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    ds_band = case_when(
      is.finite(ds_numeric) & ds_numeric <= 1 ~ "Early (1)",
      is.finite(ds_numeric) & ds_numeric <= 3 ~ "Intermediate (2-3)",
      is.finite(ds_numeric) & ds_numeric >= 4 ~ "Late (4–5)",
      TRUE ~ "Unknown"
    ),
    ds_band = factor(ds_band, levels = c("Early (1)", "Intermediate (2-3)", "Late (4–5)", "Unknown")),
    ds_factor = droplevels(as.factor(floor(ds_numeric)))
  )

res_pairs2 <- res_overlap %>%
  mutate(
    shared_prop   = ifelse(union > 0, shared / union, NA_real_),
    subset_score  = pmin(frac_inner_in_outer, frac_outer_in_inner, na.rm = TRUE),
    richness_mean = (richness_inner + richness_outer) / 2
  ) %>%
  left_join(pair_cov %>% select(pair_key, ds_numeric, ds_band, ds_factor), by = "pair_key")

isoS <- c(0.2, 0.4, 0.6, 0.8)
iso_labels <- tibble(
  L = isoS,
  x = isoS / (2 - isoS),
  y = 1,
  lab = as.character(isoS)
) %>%
  filter(is.finite(x), x >= 0, x <= 1)

p_subset <- ggplot(
  res_pairs2,
  aes(
    x = frac_inner_in_outer,
    y = frac_outer_in_inner,
    shape = position,
    size = richness_outer,
    col = ds_factor
  )
) +
  geom_abline(aes(linetype = "1:1 line"), slope = 1, intercept = 0) +
  lapply(isoS, function(S) {
    stat_function(
      data = data.frame(x = 1),
      fun  = function(x) (S * x) / (2 * x - S),
      xlim = c(S / 2 + 1e-3, 1),
      inherit.aes = FALSE,
      aes(linetype = "Iso Sørensen")
    )
  }) +
  geom_text(
    data = iso_labels,
    inherit.aes = FALSE,
    aes(x, y, label = lab),
    vjust = -0.6,
    hjust = -0.1,
    size = 3
  ) +
  geom_point(alpha = 0.85) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = ds_colors, name = "Decay stage") +
  scale_shape_manual(values = shape_map_pos, name = "Position") +
  scale_linetype_manual(
    name = "Reference",
    values = c("1:1 line" = 2, "Iso Sørensen" = 3)
  ) +
  labs(
    x = "Fraction INNER in OUTER",
    y = "Fraction OUTER in INNER",
    size = "Richness (OUTER)",
    title = "INNER vs OUTER nestedness"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right") +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    shape = guide_legend(override.aes = list(size = 3)),
    size = guide_legend(),
    linetype = guide_legend(override.aes = list(color = "black"))
  )

print(p_subset)
ggsave(out_plot_sum, p_subset, width = 7, height = 6, dpi = 600)

# ---- Directional tests and correlations --------------------------------------
wilc_rich <- wilcox.test(
  res_pairs2$richness_outer,
  res_pairs2$richness_inner,
  paired = TRUE,
  alternative = "greater"
)

df_long <- res_pairs2 %>%
  mutate(
    ds_band2 = case_when(
      is.finite(ds_numeric) & ds_numeric <= 2 ~ "Early (1–2)",
      is.finite(ds_numeric) & ds_numeric == 3 ~ "Intermediate (3)",
      is.finite(ds_numeric) & ds_numeric >= 4 ~ "Late (4–5)",
      TRUE ~ "Unknown"
    ),
    ds_band2 = factor(ds_band2, levels = c("Early (1–2)", "Intermediate (3)", "Late (4–5)", "Unknown"))
  ) %>%
  select(pair_key, ds_band2, frac_inner_in_outer, frac_outer_in_inner) %>%
  pivot_longer(starts_with("frac_"), names_to = "direction", values_to = "frac") %>%
  mutate(direction = dplyr::recode(
    direction,
    frac_inner_in_outer = "INNER in OUTER",
    frac_outer_in_inner = "OUTER in INNER"
  ))

df_early <- df_long %>%
  filter(ds_band2 == "Early (1–2)") %>%
  pivot_wider(names_from = direction, values_from = frac)

wilc_dir <- wilcox.test(
  df_early$`INNER in OUTER`,
  df_early$`OUTER in INNER`,
  paired = TRUE,
  alternative = "greater"
)

cor_out1 <- cor.test(
  res_pairs2$richness_outer,
  res_pairs2$frac_inner_in_outer,
  method = "spearman",
  alternative = "greater"
)

cor_out2 <- cor.test(
  res_pairs2$richness_inner,
  res_pairs2$frac_inner_in_outer,
  method = "spearman",
  alternative = "less"
)

cat("\n-- INNER vs OUTER (same natman×position) --\n")
cat(glue("Pairs: {nrow(res_overlap)}\n"))
cat(glue("Median Jaccard:   {round(median(res_overlap$jaccard, na.rm = TRUE), 3)}\n"))
cat(glue("Median Sørensen:  {round(median(res_overlap$sorensen, na.rm = TRUE), 3)}\n"))
cat(glue("Exact INNER ⊆ OUTER: {sum(res_overlap$is_subset_I_in_O, na.rm = TRUE)} / {nrow(res_overlap)}\n"))
cat(glue("Exact OUTER ⊆ INNER: {sum(res_overlap$is_subset_O_in_I, na.rm = TRUE)} / {nrow(res_overlap)}\n"))
cat(glue("≈Subset (INNER in OUTER, tol {subset_tol}): {sum(res_overlap$almost_I_in_O, na.rm = TRUE)}\n"))
cat(glue("≈Subset (OUTER in INNER, tol {subset_tol}): {sum(res_overlap$almost_O_in_I, na.rm = TRUE)}\n"))
cat(glue("Total SHs unique to INNER: {left_only_all}\n"))
cat(glue("Total SHs unique to OUTER: {right_only_all}\n"))
cat(glue("Total SHs shared:          {overlap_all}\n"))
cat(glue("Mean unique SHs per core (INNER): {fmt_mean_sd(res_overlap$inner_only)}\n"))
cat(glue("Mean unique SHs per core (OUTER): {fmt_mean_sd(res_overlap$outer_only)}\n"))

TN <- median(res_overlap$beta_SIM, na.rm = TRUE) / median(res_overlap$beta_SNE, na.rm = TRUE)
cat(glue("Median β_SIM: {round(median(res_overlap$beta_SIM, na.rm = TRUE), 3)}\n"))
cat(glue("Median β_SNE: {round(median(res_overlap$beta_SNE, na.rm = TRUE), 3)}\n"))
cat(glue("Median β_SOR: {round(median(res_overlap$beta_SOR, na.rm = TRUE), 3)}\n"))
cat(glue("Median T:N:   {round(TN, 3)}\n"))

cat("\nPaired Wilcoxon OUTER > INNER richness:\n")
print(wilc_rich)
cat("\nDirectional Wilcoxon in early decay (INNER in OUTER > OUTER in INNER):\n")
print(wilc_dir)
cat("\nSpearman OUTER richness vs fraction INNER in OUTER (greater):\n")
print(cor_out1)
cat("\nSpearman INNER richness vs fraction INNER in OUTER (less):\n")
print(cor_out2)
cat("\n")

# ==============================================================================
# 5) PERMANOVA with pair blocking
# ==============================================================================
out_perma <- file.path(OUTDIR_TAB, "permanova_INOUT_block_by_pair.tsv")

design_pv <- design %>%
  arrange(match(agg_key, rownames(X)))

stopifnot(identical(design_pv$agg_key, rownames(X)))

form <- as.formula("d ~ depth")

set.seed(SEED_GLOBAL)
adon_ait <- adonis2(
  update(form, d_ait ~ .),
  data = design_pv,
  permutations = N_PERM,
  by = "terms",
  strata = design_pv$pair_key
)

adon_jac <- adonis2(
  update(form, d_jac ~ .),
  data = design_pv,
  permutations = N_PERM,
  by = "terms",
  strata = design_pv$pair_key
)

cat("\n-- PERMANOVA (paired by natman×position; by = 'terms') --\n")
cat(glue("N pairs: {length(pair_has)}; N rows: {nrow(design_pv)}; perms: {N_PERM}\n"))
cat("\nAitchison PERMANOVA:\n")
print(adon_ait)
cat("\nJaccard PERMANOVA:\n")
print(adon_jac)

perma_tab <- bind_rows(
  as.data.frame(adon_ait) %>% rownames_to_column("term") %>% mutate(metric = "robust.aitchison"),
  as.data.frame(adon_jac) %>% rownames_to_column("term") %>% mutate(metric = "jaccard")
)

readr::write_tsv(perma_tab, out_perma)
message("Wrote: ", out_perma)

