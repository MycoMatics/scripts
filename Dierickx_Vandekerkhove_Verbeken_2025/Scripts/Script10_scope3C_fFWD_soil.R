# =============================================================================
# Script: Script10_scope3C_fFWD_soil.R
# Project: Deadwood fungal community analyses
# Scope: fFWD_soil
# Purpose: Paired comparison of fallen fine woody debris and linked soil cores
#          using blocked PERMANOVA, per-pair null models, AUC/ECDF summaries,
#          and betapart turnover-nestedness decomposition
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - ensure_dir
# - assert_objects
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(readr)
  library(patchwork)
})

set.seed(SEED_GLOBAL)
setwd(".")

assert_objects(c("META1", "otu_matrix_filt"))

scope <- "fFWD_soil"
outdir <- file.path("plots", "FWD", "FALLEN_SOIL")
tabdir <- file.path("tables", "FWD", "FALLEN_SOIL")
ensure_dir(outdir)
ensure_dir(tabdir)

# ---- knobs ------------------------------------------------------------------
B_NULL        <- 999L
N_PERM_ANOSIM <- 9999L
N_PERM_MAIN   <- 999L
N_PERM_NULL   <- 999L

# =============================================================================
# 1) Build fFWD–SOIL pair table
# =============================================================================
pair_key <- {
  fwd_tbl <- META1 %>%
    filter(
      dw_type2 == "fFWD",
      soil_core == "yes"
    ) %>%
    transmute(
      fwd_sample = sample,
      natman
    ) %>%
    distinct()
  
  soil_tbl <- META1 %>%
    filter(
      dw_type2 == "SOIL",
      !is.na(soil_core),
      !soil_core %in% c("yes", "no")
    ) %>%
    transmute(
      soil_sample = sample,
      fwd_sample  = soil_core
    ) %>%
    distinct()
  
  fwd_tbl %>%
    left_join(soil_tbl, by = "fwd_sample") %>%
    filter(!is.na(soil_sample)) %>%
    mutate(pair_id = fct_inorder(fwd_sample)) %>%
    arrange(natman, fwd_sample)
}

cat("\nUsable fFWD–SOIL pairs:", nrow(pair_key), "\n")
stopifnot(nrow(pair_key) >= 3)

# =============================================================================
# 2) Align OTU matrices
# =============================================================================
X_fwd  <- otu_matrix_filt[pair_key$fwd_sample,  , drop = FALSE]
X_soil <- otu_matrix_filt[pair_key$soil_sample, , drop = FALSE]

all_sh <- union(colnames(X_fwd), colnames(X_soil))

align_counts <- function(M, all_cols) {
  out <- matrix(
    0L,
    nrow = nrow(M),
    ncol = length(all_cols),
    dimnames = list(rownames(M), all_cols)
  )
  out[, colnames(M)] <- as.matrix(M)
  out
}

X_fwd  <- align_counts(X_fwd, all_sh)
X_soil <- align_counts(X_soil, all_sh)

O_count <- rbind(X_fwd, X_soil)

stopifnot(all(rownames(O_count) == c(pair_key$fwd_sample, pair_key$soil_sample)))

# =============================================================================
# 3) Metadata
# =============================================================================
meta <- tibble(
  sample    = rownames(O_count),
  type      = factor(rep(c("fFWD", "SOIL"), each = nrow(pair_key))),
  pair_id   = factor(rep(pair_key$pair_id, times = 2L)),
  natman    = factor(rep(pair_key$natman, times = 2L)),
  reads     = rowSums(O_count),
  log_reads = log1p(reads)
)

# =============================================================================
# 4) Distance matrices
# =============================================================================
as_pa <- function(M) (M > 0) * 1L

D_ait <- vegdist(O_count, method = "robust.aitchison")
D_ait <- fix_dist_labels(D_ait, samps= rownames(O_count))
D_jac <- vegdist(as_pa(O_count), method = "jaccard", binary = TRUE)

# =============================================================================
# 5) Utilities
# =============================================================================
auc_from_group <- function(D, grp) {
  M <- as.matrix(D)
  n <- nrow(M)
  stopifnot(ncol(M) == n)
  
  if (length(grp) != n) {
    stop("length(grp) must equal number of rows in D")
  }
  
  if (!is.null(rownames(M)) && !is.null(names(grp))) {
    grp <- grp[rownames(M)]
  }
  
  upper_idx <- which(upper.tri(M), arr.ind = TRUE)
  same <- grp[upper_idx[, 1]] == grp[upper_idx[, 2]]
  dU <- M[cbind(upper_idx[, 1], upper_idx[, 2])]
  
  d_between <- dU[!same]
  d_within  <- dU[same]
  
  if (length(d_between) == 0L || length(d_within) == 0L) {
    return(c(AUC = NA_real_, r_rb = NA_real_))
  }
  
  r  <- rank(c(d_between, d_within))
  n1 <- length(d_between)
  n2 <- length(d_within)
  U  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  A  <- U / (n1 * n2)
  
  c(AUC = A, r_rb = 2 * A - 1)
}

pairwise_groups <- function(D, grp, upper_idx) {
  M <- as.matrix(D)
  dU <- M[cbind(upper_idx[, 1], upper_idx[, 2])]
  same <- grp[upper_idx[, 1]] == grp[upper_idx[, 2]]
  list(within = dU[same], between = dU[!same])
}

row_q <- function(mat, prob) {
  apply(mat, 1L, stats::quantile, probs = prob, na.rm = TRUE)
}

random_derangement <- function(n) {
  p <- sample.int(n, n, replace = FALSE)
  fp <- which(p == seq_len(n))
  if (length(fp) == 1L) {
    i <- fp
    j <- if (i < n) i + 1L else 1L
    tmp <- p[i]
    p[i] <- p[j]
    p[j] <- tmp
  } else if (length(fp) > 1L) {
    p[fp] <- p[fp[c(2:length(fp), 1L)]]
  }
  p
}

p_left <- function(null, obs) {
  (rowSums(null <= obs) + 1) / (ncol(null) + 1)
}

# =============================================================================
# 6) Blocked PERMANOVA for type effect
# =============================================================================
ctrl_pair <- how(nperm = N_PERM_MAIN)
setBlocks(ctrl_pair) <- meta$pair_id

ado_ait_blk <- adonis2(
  D_ait ~ log_reads + type,
  data = meta,
  permutations = ctrl_pair,
  by = "terms"
)

ado_jac_blk <- adonis2(
  D_jac ~ log_reads + type,
  data = meta,
  permutations = ctrl_pair,
  by = "terms"
)

cat("\n[Blocked PERMANOVA] type effect with pair as block\n")
print(ado_ait_blk)
print(ado_jac_blk)

# =============================================================================
# 7) PERMDISP
# =============================================================================
disp_type_Ait <- anova(betadisper(D_ait, meta$type))
disp_type_Jac <- anova(betadisper(D_jac, meta$type))

print(list(
  Aitchison_type = disp_type_Ait,
  Jaccard_type   = disp_type_Jac
))

sink(file.path(outdir, paste0("permdisp_", scope, ".txt")))
cat("Aitchison_type\n"); print(disp_type_Ait)
cat("Jaccard_type\n");   print(disp_type_Jac)
sink()

# =============================================================================
# 8) Per-pair soil imprint null tests
# =============================================================================
PA_fwd  <- as_pa(X_fwd)
PA_soil <- as_pa(X_soil)

intersect_mat <- PA_fwd %*% t(PA_soil)
n_f <- rowSums(PA_fwd)
n_s <- rowSums(PA_soil)
union_mat <- outer(n_f, n_s, "+") - intersect_mat

jac_mat_sim <- intersect_mat / pmax(union_mat, 1L)
jac_mat <- 1 - jac_mat_sim
obs_jac <- diag(jac_mat)

D_ait_mat <- as.matrix(D_ait)
idx_f <- match(pair_key$fwd_sample, rownames(O_count))
idx_s <- match(pair_key$soil_sample, rownames(O_count))
stopifnot(!anyNA(idx_f), !anyNA(idx_s))

ait_cross <- D_ait_mat[idx_f, idx_s, drop = FALSE]
obs_ait   <- diag(ait_cross)

n_pairs   <- nrow(pair_key)
nperm     <- 4999L
grp_site  <- factor(pair_key$natman)

perm_idx_blocked <- replicate(
  nperm,
  {
    unlist(
      lapply(split(seq_len(n_pairs), grp_site), function(i) sample(i, length(i))),
      use.names = FALSE
    )
  }
)

perm_idx_global <- replicate(
  nperm,
  random_derangement(n_pairs),
  simplify = FALSE
)

jac_null_blocked <- matrix(NA_real_, nrow = n_pairs, ncol = nperm)
jac_null_global  <- matrix(NA_real_, nrow = n_pairs, ncol = nperm)
ait_null_blocked <- matrix(NA_real_, nrow = n_pairs, ncol = nperm)
ait_null_global  <- matrix(NA_real_, nrow = n_pairs, ncol = nperm)

for (p in seq_len(nperm)) {
  perm_b <- perm_idx_blocked[, p]
  perm_g <- perm_idx_global[[p]]
  
  jac_null_blocked[, p] <- jac_mat[cbind(seq_len(n_pairs), perm_b)]
  jac_null_global[,  p] <- jac_mat[cbind(seq_len(n_pairs), perm_g)]
  
  ait_null_blocked[, p] <- ait_cross[cbind(seq_len(n_pairs), perm_b)]
  ait_null_global[,  p] <- ait_cross[cbind(seq_len(n_pairs), perm_g)]
}

res_pairs <- pair_key %>%
  mutate(
    obs_jaccard   = obs_jac,
    obs_ait       = obs_ait,
    p_jac_blocked = p_left(jac_null_blocked, obs_jac),
    p_jac_global  = p_left(jac_null_global,  obs_jac),
    p_ait_blocked = p_left(ait_null_blocked, obs_ait),
    p_ait_global  = p_left(ait_null_global,  obs_ait)
  ) %>%
  mutate(
    q_jac_blocked = p.adjust(p_jac_blocked, method = "BH"),
    q_jac_global  = p.adjust(p_jac_global,  method = "BH"),
    q_ait_blocked = p.adjust(p_ait_blocked, method = "BH"),
    q_ait_global  = p.adjust(p_ait_global,  method = "BH")
  )

cat("\nPer pair soil imprint results\n")
print(res_pairs)

# ---- pairwise null summaries -----------------------------------------------
jac_block_mean <- rowMeans(jac_null_blocked, na.rm = TRUE)
jac_block_lo   <- row_q(jac_null_blocked, 0.025)
jac_block_hi   <- row_q(jac_null_blocked, 0.975)

jac_glob_mean  <- rowMeans(jac_null_global, na.rm = TRUE)
jac_glob_lo    <- row_q(jac_null_global, 0.025)
jac_glob_hi    <- row_q(jac_null_global, 0.975)

ait_block_mean <- rowMeans(ait_null_blocked, na.rm = TRUE)
ait_block_lo   <- row_q(ait_null_blocked, 0.025)
ait_block_hi   <- row_q(ait_null_blocked, 0.975)

ait_glob_mean  <- rowMeans(ait_null_global, na.rm = TRUE)
ait_glob_lo    <- row_q(ait_null_global, 0.025)
ait_glob_hi    <- row_q(ait_null_global, 0.975)

annot_pairs <- pair_key %>%
  left_join(
    META1 %>%
      filter(dw_type2 == "fFWD") %>%
      transmute(
        fwd_sample   = sample,
        natman,
        decay_stage,
        size,
        decay_x_size = paste(decay_stage, size, sep = " · ")
      ) %>%
      distinct(fwd_sample, natman, decay_x_size),
    by = c("fwd_sample", "natman")
  )

tab_jaccard <- annot_pairs %>%
  transmute(
    fwd_sample          = fwd_sample,
    tree_ID             = natman,
    decay_x_size        = decay_x_size,
    metric              = "Jaccard",
    obs                 = obs_jac,
    null_mean_blocked   = jac_block_mean,
    ci_blocked          = sprintf("[%.3f, %.3f]", jac_block_lo, jac_block_hi),
    null_mean_unblocked = jac_glob_mean,
    ci_unblocked        = sprintf("[%.3f, %.3f]", jac_glob_lo, jac_glob_hi),
    p_blocked           = res_pairs$p_jac_blocked,
    p_unblocked         = res_pairs$p_jac_global
  )

tab_ait <- annot_pairs %>%
  transmute(
    fwd_sample          = fwd_sample,
    tree_ID             = natman,
    decay_x_size        = decay_x_size,
    metric              = "rAit",
    obs                 = obs_ait,
    null_mean_blocked   = ait_block_mean,
    ci_blocked          = sprintf("[%.1f, %.1f]", ait_block_lo, ait_block_hi),
    null_mean_unblocked = ait_glob_mean,
    ci_unblocked        = sprintf("[%.1f, %.1f]", ait_glob_lo, ait_glob_hi),
    p_blocked           = res_pairs$p_ait_blocked,
    p_unblocked         = res_pairs$p_ait_global
  )

pair_null_table <- bind_rows(tab_jaccard, tab_ait) %>%
  arrange(tree_ID, fwd_sample, factor(metric, levels = c("Jaccard", "rAit")))

cat("\nPer pair soil imprint summary table\n")
print(pair_null_table, n = Inf)

write_csv2(res_pairs,        file.path(tabdir, "fFWD_SOIL_soil_imprint_pairs_minimal.csv"))
write_csv2(pair_null_table,  file.path(tabdir, "fFWD_SOIL_soil_imprint_pair_summary.csv"))

# =============================================================================
# 9) Pairwise null plot data
# =============================================================================
base_df <- annot_pairs %>%
  mutate(
    obs_jac   = obs_jac,
    obs_ait   = obs_ait,
    p_jac_blk = res_pairs$p_jac_blocked,
    p_jac_glb = res_pairs$p_jac_global,
    p_ait_blk = res_pairs$p_ait_blocked,
    p_ait_glb = res_pairs$p_ait_global
  )

df_jac_long <- bind_rows(
  base_df %>%
    transmute(
      fwd_sample,
      natman,
      axis_id        = paste(natman, fwd_sample, sep = "::"),
      axis_lab       = decay_x_size      ,
      axis_lab_clean = gsub("_", " ", axis_lab),
      obs        = obs_jac,
      p_jac_blk,
      p_jac_glb,
      above_blk   = obs_jac >= jac_block_mean,
      above_glb   = obs_jac >= jac_glob_mean,
      sig_blk_dir = obs_jac > jac_block_hi,
      sig_glb_dir = obs_jac > jac_glob_hi,
      type       = "blocked null",
      null_mean  = jac_block_mean,
      lo95       = jac_block_lo,
      hi95       = jac_block_hi
    ),
  base_df %>%
    transmute(
      fwd_sample,
      natman,
      axis_id        = paste(natman, fwd_sample, sep = "::"),
      axis_lab       = decay_x_size      ,
      axis_lab_clean = gsub("_", " ", axis_lab),
      obs        = obs_jac,
      p_jac_blk,
      p_jac_glb,
      above_blk   = obs_jac >= jac_glob_mean,
      above_glb   = obs_jac >= jac_glob_mean,
      sig_blk_dir = obs_jac > jac_glob_hi,
      sig_glb_dir = obs_jac > jac_glob_hi,
      type       = "global null",
      null_mean  = jac_glob_mean,
      lo95       = jac_glob_lo,
      hi95       = jac_glob_hi
    )
) %>%
  arrange(natman, fwd_sample, type)

df_ait_long <- bind_rows(
  base_df %>%
    transmute(
      fwd_sample,
      natman,
      axis_id        = paste(natman, fwd_sample, sep = "::"),
      axis_lab       = decay_x_size      ,
      axis_lab_clean = gsub("_", " ", axis_lab),
      obs        = obs_ait,
      p_ait_blk,
      p_ait_glb,
      above_blk   = obs_ait <= ait_block_mean,
      above_glb   = obs_ait <= ait_block_mean,
      sig_blk_dir = obs_ait < ait_block_lo,
      sig_glb_dir = obs_ait < ait_block_lo,
      type       = "blocked null",
      null_mean  = ait_block_mean,
      lo95       = ait_block_lo,
      hi95       = ait_block_hi
    ),
  base_df %>%
    transmute(
      fwd_sample,
      natman,
      axis_id        = paste(natman, fwd_sample, sep = "::"),
      axis_lab       = decay_x_size      ,
      axis_lab_clean = gsub("_", " ", axis_lab),
      obs        = obs_ait,
      p_ait_blk,
      p_ait_glb,
      above_blk   = obs_ait <= ait_glob_mean,
      above_glb   = obs_ait <= ait_glob_mean,
      sig_blk_dir = obs_ait < ait_glob_lo,
      sig_glb_dir = obs_ait < ait_glob_lo,
      type       = "global null",
      null_mean  = ait_glob_mean,
      lo95       = ait_glob_lo,
      hi95       = ait_glob_hi
    )
) %>%
  arrange(natman, fwd_sample, type)

axis_info <- df_jac_long %>%
  distinct(axis_id, axis_lab_clean, natman) %>%
  arrange(natman, axis_id)

axis_levels <- axis_info$axis_id
jac_labels  <- setNames(axis_info$axis_lab_clean, axis_info$axis_id)
ait_labels  <- jac_labels

dodge_w  <- 0.5
pad_add  <- 0.4
pd       <- position_dodge2(width = dodge_w, preserve = "single")
col_blocked <- "#1b9e77"
col_global  <- "#d95f02"

gg_jac <- ggplot() +
  geom_linerange(
    data = df_jac_long,
    aes(x = factor(axis_id, levels = axis_levels), ymin = lo95, ymax = hi95, color = type, group = type),
    position = pd, linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = df_jac_long,
    aes(x = factor(axis_id, levels = axis_levels), y = null_mean, color = type, group = type),
    position = pd, size = 1.5, alpha = 0.9
  ) +
  geom_point(
    data = df_jac_long %>% distinct(axis_id, natman, obs, axis_lab_clean),
    aes(x = factor(axis_id, levels = axis_levels), y = obs),
    shape = 23, stroke = 1.0, size = 2,
    fill = "black", color = "black", alpha = 0.7
  ) +
  facet_grid(natman ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = jac_labels, expand = expansion(add = pad_add)) +
  scale_color_manual(
    name = NULL,
    values = c("blocked null" = col_blocked, "global null" = col_global),
    labels = c("blocked null" = "Tree blocked null", "global null" = "Unblocked null")
  ) +
  coord_flip() +
  labs(
    title = "FWD–SOIL Jaccard with tree blocked vs unblocked nulls\nleft = more similar",
    x = "Decay stage · Size",
    y = "Jaccard(FWD, matched SOIL)"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.spacing.y = grid::unit(8, "pt"), legend.position = "bottom")

gg_ait <- ggplot() +
  geom_linerange(
    data = df_ait_long,
    aes(x = factor(axis_id, levels = axis_levels), ymin = lo95, ymax = hi95, color = type, group = type),
    position = pd, linewidth = 0.8, alpha = 0.6
  ) +
  geom_point(
    data = df_ait_long,
    aes(x = factor(axis_id, levels = axis_levels), y = null_mean, color = type, group = type),
    position = pd, size = 1.5, alpha = 0.9
  ) +
  geom_point(
    data = df_ait_long %>% distinct(axis_id, natman, obs, axis_lab_clean),
    aes(x = factor(axis_id, levels = axis_levels), y = obs),
    shape = 23, stroke = 1.0, size = 2,
    fill = "black", color = "black", alpha = 0.7
  ) +
  facet_grid(natman ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = ait_labels, expand = expansion(add = pad_add)) +
  scale_color_manual(
    name = NULL,
    values = c("blocked null" = col_blocked, "global null" = col_global),
    labels = c("blocked null" = "Tree blocked null", "global null" = "Unblocked null")
  ) +
  coord_flip() +
  labs(
    title = "FWD–SOIL robust Aitchison\nleft = more similar",
    x = "Decay stage · Size",
    y = "Robust Aitchison distance (FWD, matched SOIL)"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.spacing.y = grid::unit(8, "pt"), legend.position = "bottom")

print(gg_jac)
print(gg_ait)

ggsave(file.path(outdir, "fwd_soil_pairs_obs_vs_twoNulls_Jaccard_dodged.png"), gg_jac, width = 7.2, height = 6.6, dpi = 300)
ggsave(file.path(outdir, "fwd_soil_pairs_obs_vs_twoNulls_rAit_dodged.png"), gg_ait, width = 7.2, height = 6.6, dpi = 300)

gg_jac2 <- gg_jac +
  labs(title = "FWD–SOIL Jaccard") +
  theme(plot.title = element_text(size = 10), legend.position = "none")

gg_ait2 <- gg_ait +
  labs(title = "") +
  theme(plot.title = element_text(size = 10), legend.position = "bottom")

gg_panel <- gg_jac2 / gg_ait2 + patchwork::plot_annotation(tag_levels = "A")

print(gg_panel)

ggsave(file.path(outdir, "fwd_soil_pairs_obs_vs_twoNulls_PANEL.png"), gg_panel, width = 17, height = 25, units = "cm", dpi = 900)

# =============================================================================
# 10) Global type effect evidence
# =============================================================================
scope_name <- "FWD_SOIL"

R2_obs_ait <- as.data.frame(ado_ait_blk)["type", "R2"]
R2_obs_jac <- as.data.frame(ado_jac_blk)["type", "R2"]

grp_type <- droplevels(meta$type)

AUC_obs_ait <- auc_from_group(D_ait, grp_type)
AUC_obs_jac <- auc_from_group(D_jac, grp_type)

anos_ait <- vegan::anosim(
  D_ait,
  grouping = grp_type,
  permutations = N_PERM_ANOSIM,
  strata = meta$pair_id
)
anos_jac <- vegan::anosim(
  D_jac,
  grouping = grp_type,
  permutations = N_PERM_ANOSIM,
  strata = meta$pair_id
)

cat("\n[ANOSIM] type effect (blocked by pair_id)\n")
print(anos_ait)
print(anos_jac)

R2_null_ait  <- numeric(B_NULL)
R2_null_jac  <- numeric(B_NULL)
AUC_null_ait <- numeric(B_NULL)
AUC_null_jac <- numeric(B_NULL)

for (b in seq_len(B_NULL)) {
  perm_type <- grp_type
  flips <- rbinom(n_pairs, size = 1L, prob = 0.5)
  
  for (i in seq_len(n_pairs)) {
    if (flips[i] == 1L) {
      i1 <- i
      i2 <- n_pairs + i
      tmp <- perm_type[i1]
      perm_type[i1] <- perm_type[i2]
      perm_type[i2] <- tmp
    }
  }
  
  meta_tmp <- meta
  meta_tmp$type <- perm_type
  
  ado_tmp_ait <- adonis2(
    D_ait ~ type,
    data = meta_tmp,
    permutations = ctrl_pair,
    by = "margin"
  )
  ado_tmp_jac <- adonis2(
    D_jac ~ type,
    data = meta_tmp,
    permutations = ctrl_pair,
    by = "margin"
  )
  
  R2_null_ait[b]  <- as.data.frame(ado_tmp_ait)["type", "R2"]
  R2_null_jac[b]  <- as.data.frame(ado_tmp_jac)["type", "R2"]
  AUC_null_ait[b] <- auc_from_group(D_ait, meta_tmp$type)["AUC"]
  AUC_null_jac[b] <- auc_from_group(D_jac, meta_tmp$type)["AUC"]
  
  if (b %% 50L == 0L) {
    cat("[type null] finished", b, "of", B_NULL, "\n")
  }
}

ciR2_ait <- stats::quantile(R2_null_ait, c(0.025, 0.975), na.rm = TRUE)
ciR2_jac <- stats::quantile(R2_null_jac, c(0.025, 0.975), na.rm = TRUE)

pR2_ait <- mean(R2_null_ait >= R2_obs_ait, na.rm = TRUE)
pR2_jac <- mean(R2_null_jac >= R2_obs_jac, na.rm = TRUE)

sesR2_ait <- (R2_obs_ait - mean(R2_null_ait, na.rm = TRUE)) / stats::sd(R2_null_ait, na.rm = TRUE)
sesR2_jac <- (R2_obs_jac - mean(R2_null_jac, na.rm = TRUE)) / stats::sd(R2_null_jac, na.rm = TRUE)

R2_corr_ait <- R2_obs_ait - mean(R2_null_ait, na.rm = TRUE)
R2_corr_jac <- R2_obs_jac - mean(R2_null_jac, na.rm = TRUE)

ciAUC_ait <- stats::quantile(AUC_null_ait, c(0.025, 0.975), na.rm = TRUE)
ciAUC_jac <- stats::quantile(AUC_null_jac, c(0.025, 0.975), na.rm = TRUE)

pAUC_ait <- mean(AUC_null_ait >= AUC_obs_ait["AUC"], na.rm = TRUE)
pAUC_jac <- mean(AUC_null_jac >= AUC_obs_jac["AUC"], na.rm = TRUE)

sesAUC_ait <- (AUC_obs_ait["AUC"] - mean(AUC_null_ait, na.rm = TRUE)) / stats::sd(AUC_null_ait, na.rm = TRUE)
sesAUC_jac <- (AUC_obs_jac["AUC"] - mean(AUC_null_jac, na.rm = TRUE)) / stats::sd(AUC_null_jac, na.rm = TRUE)

summary_type <- tibble(
  Scope        = scope_name,
  Metric       = c("robust.aitchison", "jaccard"),
  R2_obs       = c(R2_obs_ait, R2_obs_jac),
  R2_null_mean = c(mean(R2_null_ait, na.rm = TRUE), mean(R2_null_jac, na.rm = TRUE)),
  R2_null_CI_lo = c(ciR2_ait[1], ciR2_jac[1]),
  R2_null_CI_hi = c(ciR2_ait[2], ciR2_jac[2]),
  p_R2         = c(pR2_ait, pR2_jac),
  R2_SES       = c(sesR2_ait, sesR2_jac),
  R2_corrected = c(R2_corr_ait, R2_corr_jac),
  AUC_obs      = c(unname(AUC_obs_ait["AUC"]), unname(AUC_obs_jac["AUC"])),
  r_rb         = c(unname(AUC_obs_ait["r_rb"]), unname(AUC_obs_jac["r_rb"])),
  AUC_null_mean = c(mean(AUC_null_ait, na.rm = TRUE), mean(AUC_null_jac, na.rm = TRUE)),
  AUC_null_CI_lo = c(ciAUC_ait[1], ciAUC_jac[1]),
  AUC_null_CI_hi = c(ciAUC_ait[2], ciAUC_jac[2]),
  p_AUC        = c(pAUC_ait, pAUC_jac),
  AUC_SES      = c(sesAUC_ait, sesAUC_jac),
  ANOSIM_R     = c(anos_ait$statistic, anos_jac$statistic),
  ANOSIM_p     = c(anos_ait$signif, anos_jac$signif),
  n_trees      = dplyr::n_distinct(pair_key$natman),
  n_samples    = nrow(meta)
)

cat("\nSummary for table filling (fFWD vs SOIL):\n")
print(summary_type)

write_csv2(summary_type, file.path(tabdir, "fFWD_SOIL_type_evidence_summary.csv"))

# ---- null distribution plots -----------------------------------------------
base_null_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    plot.title   = element_text(hjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_text(margin = margin(r = 5))
  )

df_R2_ait <- tibble::tibble(R2 = R2_null_ait)
p_R2_ait <- ggplot(df_R2_ait, aes(x = R2)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = R2_obs_ait, linetype = "dashed") +
  labs(
    title = "Null R2 for type (robust Aitchison)",
    x = "R2 under within pair label shuffles",
    y = "Density"
  ) +
  base_null_theme

print(p_R2_ait)
ggsave(file.path(outdir, "FWD_SOIL_R2_NULL_robust_aitchison.png"), p_R2_ait, width = 6, height = 4, dpi = 300)

df_R2_jac <- tibble::tibble(R2 = R2_null_jac)
p_R2_jac <- ggplot(df_R2_jac, aes(x = R2)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = R2_obs_jac, linetype = "dashed") +
  labs(
    title = "Null R2 for type (Jaccard)",
    x = "R2 under within pair label shuffles",
    y = "Density"
  ) +
  base_null_theme

print(p_R2_jac)
ggsave(file.path(outdir, "FWD_SOIL_R2_NULL_jaccard.png"), p_R2_jac, width = 6, height = 4, dpi = 300)

df_AUC_ait <- tibble::tibble(AUC = AUC_null_ait)
p_AUC_ait <- ggplot(df_AUC_ait, aes(x = AUC)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = AUC_obs_ait["AUC"], linetype = "dashed") +
  labs(
    title = "Null AUC for type (robust Aitchison)",
    x = "AUC under within pair label shuffles",
    y = "Density"
  ) +
  base_null_theme

print(p_AUC_ait)
ggsave(file.path(outdir, "FWD_SOIL_AUC_NULL_robust_aitchison.png"), p_AUC_ait, width = 6, height = 4, dpi = 300)

df_AUC_jac <- tibble::tibble(AUC = AUC_null_jac)
p_AUC_jac <- ggplot(df_AUC_jac, aes(x = AUC)) +
  geom_density(fill = "grey80", color = "black") +
  geom_vline(xintercept = AUC_obs_jac["AUC"], linetype = "dashed") +
  labs(
    title = "Null AUC for type (Jaccard)",
    x = "AUC under within pair label shuffles",
    y = "Density"
  ) +
  base_null_theme

print(p_AUC_jac)
ggsave(file.path(outdir, "FWD_SOIL_AUC_NULL_jaccard.png"), p_AUC_jac, width = 6, height = 4, dpi = 300)

# =============================================================================
# 11) Betapart turnover vs nestedness
# =============================================================================
PA_all <- (O_count > 0) * 1L

bp <- betapart::beta.pair(PA_all, index.family = "sorensen")
DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

grp_state <- droplevels(meta$type)

upper_idx <- which(upper.tri(DM_sim), arr.ind = TRUE)
state_i <- grp_state[upper_idx[, 1]]
state_j <- grp_state[upper_idx[, 2]]
pair_class <- ifelse(state_i == state_j, "within_substrate", "between_substrate")

betapart_pairs <- tibble(
  pair_class = factor(pair_class, levels = c("within_substrate", "between_substrate")),
  beta_sim = DM_sim[upper_idx],
  beta_sne = DM_sne[upper_idx],
  beta_sor = DM_sor[upper_idx]
)

betapart_summary <- betapart_pairs %>%
  group_by(pair_class) %>%
  summarise(
    n_pairs                   = n(),
    beta_sim_median           = median(beta_sim, na.rm = TRUE),
    beta_sne_median           = median(beta_sne, na.rm = TRUE),
    beta_sor_median           = median(beta_sor, na.rm = TRUE),
    ratio_turnover_nestedness = beta_sim_median / beta_sne_median,
    .groups = "drop"
  )

cat("\n[BETAPART] Turnover vs nestedness by substrate pairing (fFWD vs SOIL)\n")
print(betapart_summary)

write_csv2(betapart_summary, file.path(tabdir, "fFWD_SOIL_betapart_within_between_substrate.csv"))

# =============================================================================
# 12) ECDFs of pairwise distances by substrate type
# =============================================================================
upper_idx <- which(upper.tri(as.matrix(D_ait)), arr.ind = TRUE)

pw_A <- pairwise_groups(D_ait, grp_type, upper_idx)
df_ecdf_A <- tibble(
  dist = c(pw_A$within, pw_A$between),
  type = factor(
    rep(c("within", "between"), c(length(pw_A$within), length(pw_A$between))),
    levels = c("within", "between")
  )
)

p_ecdf_A <- ggplot2::ggplot(df_ecdf_A, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = paste0(scope, " — ECDF of substrate type (robust Aitchison)"),
    subtitle = sprintf("AUC = %.3f", unname(AUC_obs_ait["AUC"])),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

ggplot2::ggsave(file.path(outdir, paste0("type_ECDF_robustAitchison_", scope, ".png")), p_ecdf_A, width = 6.2, height = 4.2, dpi = 300)
print(p_ecdf_A)

pw_J <- pairwise_groups(D_jac, grp_type, upper_idx)
df_ecdf_J <- tibble(
  dist = c(pw_J$within, pw_J$between),
  type = factor(
    rep(c("within", "between"), c(length(pw_J$within), length(pw_J$between))),
    levels = c("within", "between")
  )
)

p_ecdf_J <- ggplot2::ggplot(df_ecdf_J, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = paste0(scope, " — ECDF of substrate type (Jaccard)"),
    subtitle = sprintf("AUC = %.3f", unname(AUC_obs_jac["AUC"])),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

ggplot2::ggsave(file.path(outdir, paste0("type_ECDF_Jaccard_", scope, ".png")), p_ecdf_J, width = 6.2, height = 4.2, dpi = 300)
print(p_ecdf_J)
