# =============================================================================
# Script: Script8_scope3A_LOG_aFWD.R
# Project: Deadwood fungal community analyses
# Scope: LOG_aFWD
# Purpose: Analyse outer LOG cores and attached fine woody debris (aFWD)
#          using PERMANOVA, blocked PERMANOVA, natman null models for R2,
#          ANOSIM, leave-one-tree-out stability, co-location tests, ECDFs,
#          and betapart decomposition
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - fix_dist_labels
# - get_R2_term
# - loo_permanova
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(readr)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c("fix_dist_labels", "get_R2_term"))

set.seed(SEED_GLOBAL)
setwd(".")

scope  <- "LOG_aFWD"
outdir <- file.path("plots", "LOG_aFWD")
ensure_dir(outdir)

N_PERM_MAIN    <- 999L
B_NULL_AUC     <- 5000L
N_PERM_NULL_R2 <- 199L
B_NULL_R2      <- 999L
N_PERM_ANOSIM  <- 9999L
N_PERM_LOO     <- 499L

cat("\n[LOG vs aFWD ] Scope:", scope, "\n")

# -----------------------------------------------------------------------------
# 1) Data filtering
# -----------------------------------------------------------------------------
meta <- META1 %>%
  filter(dw_type2 %in% c("LOG", "aFWD")) %>%
  filter(position_2 %in% c("BASE", "MIDDLE", "UPPER", "AFWD1", "AFWD2")) %>%
  filter(is.na(depth_2) | depth_2 == "OUTER") %>%
  filter(sample %in% rownames(otu_matrix_filt)) %>%
  mutate(
    natman      = droplevels(factor(natman)),
    dw_type2    = droplevels(factor(dw_type2)),
    ds_at_drill = droplevels(factor(ds_at_drill)),
    umi         = droplevels(factor(umi))
  ) %>%
  group_by(natman) %>%
  filter(sum(dw_type2 == "aFWD") >= 2) %>%
  dplyr::rename(dw_type = dw_type2) %>%
  ungroup() %>%
  arrange(natman, dw_type, sample)

stopifnot(nrow(meta) > 0)
stopifnot(all(meta$sample %in% rownames(otu_matrix_filt)))

otu <- otu_matrix_filt[meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu), meta$sample))
otu <- otu[, colSums(otu) > 0, drop = FALSE]

if (!"log_reads" %in% names(meta)) {
  meta <- meta %>% mutate(log_reads = log1p(rowSums(otu)))
}

cat("Samples in analysis:", nrow(otu), "  SHs:", ncol(otu), "\n")
cat("trees in analysis:", length(unique(meta$natman)), "\n")
cat("dw_type counts:\n")
print(table(meta$dw_type))

# -----------------------------------------------------------------------------
# 2) Distance matrices
# -----------------------------------------------------------------------------
cat("\n[Distances] robust Aitchison and Jaccard\n")
D_robAit <- vegdist(otu, method = "robust.aitchison")
D_robAit <- fix_dist_labels(D_robAit, samps = rownames(otu))

otu_pa <- (otu > 0) * 1L
D_jacc <- vegdist(otu_pa, method = "jaccard", binary = TRUE)
D_jacc <- fix_dist_labels(D_jacc, samps = rownames(otu))

# -----------------------------------------------------------------------------
# 3) PERMANOVA unblocked
# -----------------------------------------------------------------------------
cat("\n[PERMANOVA] Unblocked (by = 'terms')\n")
perma_unblocked <- list(
  robAit = adonis2(
    D_robAit ~ log_reads + umi + natman + dw_type,
    data = meta,
    permutations = N_PERM_MAIN,
    by = "terms"
  ),
  jaccard = adonis2(
    D_jacc ~ log_reads + umi + natman + dw_type,
    data = meta,
    permutations = N_PERM_MAIN,
    by = "terms"
  )
)

print(perma_unblocked)
capture.output(perma_unblocked, file = file.path(outdir, paste0("permanova_unblocked_", scope, ".txt")))

# -----------------------------------------------------------------------------
# 4) Blocked by natman
# -----------------------------------------------------------------------------
cat("\n[PERMANOVA] Blocked by natman (restricted permutations within tree)\n")
ctrl_nat <- permute::how(nperm = N_PERM_MAIN)
permute::setBlocks(ctrl_nat) <- meta$natman

perma_blocked_natman <- list(
  robAit = adonis2(
    D_robAit ~ log_reads + umi + dw_type,
    data = meta,
    permutations = ctrl_nat,
    by = "margin"
  ),
  jaccard = adonis2(
    D_jacc ~ log_reads + umi + dw_type,
    data = meta,
    permutations = ctrl_nat,
    by = "margin"
  )
)

print(perma_blocked_natman)
capture.output(perma_blocked_natman, file = file.path(outdir, paste0("permanova_blocked_natman_", scope, ".txt")))

# -----------------------------------------------------------------------------
# 5) NATMAN R2 null model
# -----------------------------------------------------------------------------
cat("\n[NATMAN R2 nulls] label shuffle under unblocked models\n")

R2_obs_Ait_natman <- get_R2_term(perma_unblocked$robAit, "natman")
R2_obs_Jac_natman <- get_R2_term(perma_unblocked$jaccard, "natman")

cat("Observed natman R2  Aitchison:", R2_obs_Ait_natman, "  Jaccard:", R2_obs_Jac_natman, "\n")

R2_null_Ait <- numeric(B_NULL_R2)
R2_null_Jac <- numeric(B_NULL_R2)
nat <- droplevels(meta$natman)

for (b in seq_len(B_NULL_R2)) {
  nat_perm <- factor(sample(nat, replace = FALSE), levels = levels(nat))
  meta_perm <- meta
  meta_perm$natman <- nat_perm
  
  fit_Ait <- adonis2(
    D_robAit ~ log_reads + umi + natman + dw_type,
    data = meta_perm,
    permutations = N_PERM_NULL_R2,
    by = "terms"
  )
  
  fit_Jac <- adonis2(
    D_jacc ~ log_reads + umi + natman + dw_type,
    data = meta_perm,
    permutations = N_PERM_NULL_R2,
    by = "terms"
  )
  
  R2_null_Ait[b] <- get_R2_term(fit_Ait, "natman")
  R2_null_Jac[b] <- get_R2_term(fit_Jac, "natman")
  
  if (b %% 50 == 0) {
    cat("[natman R2 nulls] finished", b, "of", B_NULL_R2, "\n")
  }
}

saveRDS(R2_null_Ait, file.path(outdir, "R2_null_Ait_natman_LOGaFWD.rds"))
saveRDS(R2_null_Jac, file.path(outdir, "R2_null_Jac_natman_LOGaFWD.rds"))

ci_Ait   <- quantile(R2_null_Ait, c(0.025, 0.975), na.rm = TRUE)
ci_Jac   <- quantile(R2_null_Jac, c(0.025, 0.975), na.rm = TRUE)
mean_Ait <- mean(R2_null_Ait, na.rm = TRUE)
mean_Jac <- mean(R2_null_Jac, na.rm = TRUE)

ses_Ait <- (R2_obs_Ait_natman - mean_Ait) / sd(R2_null_Ait, na.rm = TRUE)
ses_Jac <- (R2_obs_Jac_natman - mean_Jac) / sd(R2_null_Jac, na.rm = TRUE)

p_emp_Ait <- mean(R2_null_Ait >= R2_obs_Ait_natman, na.rm = TRUE)
p_emp_Jac <- mean(R2_null_Jac >= R2_obs_Jac_natman, na.rm = TRUE)

perc_Ait <- mean(R2_null_Ait <= R2_obs_Ait_natman, na.rm = TRUE)
perc_Jac <- mean(R2_null_Jac <= R2_obs_Jac_natman, na.rm = TRUE)

natman_R2_null_summary <- tibble(
  scope         = scope,
  metric        = c("robust.aitchison", "jaccard"),
  R2_obs_natman = c(R2_obs_Ait_natman, R2_obs_Jac_natman),
  R2_null_mean  = c(mean_Ait, mean_Jac),
  R2_null_CI_lo = c(ci_Ait[1], ci_Jac[1]),
  R2_null_CI_hi = c(ci_Ait[2], ci_Jac[2]),
  delta_R2      = c(R2_obs_Ait_natman - mean_Ait, R2_obs_Jac_natman - mean_Jac),
  R2_SES        = c(ses_Ait, ses_Jac),
  p_empirical   = c(p_emp_Ait, p_emp_Jac),
  R2_percentile = c(perc_Ait, perc_Jac)
)

print(natman_R2_null_summary)
write_csv2(natman_R2_null_summary, file.path(outdir, paste0("natman_R2_null_summary_", scope, ".csv")))

pR2_Ait <- ggplot2::ggplot(
  tibble::tibble(R2_null = R2_null_Ait),
  ggplot2::aes(x = R2_null)
) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = R2_obs_Ait_natman, linetype = 2) +
  ggplot2::labs(
    title = paste0("natman R2 null vs observed (Aitchison) - ", scope),
    subtitle = sprintf(
      "R2_obs = %.3f  null_mean = %.3f  95%%CI = [%.3f, %.3f]  SES = %.2f",
      R2_obs_Ait_natman, mean_Ait, ci_Ait[1], ci_Ait[2], ses_Ait
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(outdir, paste0("natman_R2_null_Aitchison_", scope, ".png")), pR2_Ait, width = 6, height = 4, dpi = 300)
print(pR2_Ait)

pR2_Jac <- ggplot2::ggplot(
  tibble::tibble(R2_null = R2_null_Jac),
  ggplot2::aes(x = R2_null)
) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = R2_obs_Jac_natman, linetype = 2) +
  ggplot2::labs(
    title = paste0("natman R2 null vs observed (Jaccard) - ", scope),
    subtitle = sprintf(
      "R2_obs = %.3f  null_mean = %.3f  95%%CI = [%.3f, %.3f]  SES = %.2f",
      R2_obs_Jac_natman, mean_Jac, ci_Jac[1], ci_Jac[2], ses_Jac
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(file.path(outdir, paste0("natman_R2_null_Jaccard_", scope, ".png")), pR2_Jac, width = 6, height = 4, dpi = 300)
print(pR2_Jac)

# -----------------------------------------------------------------------------
# 6) PERMANOVA blocked duplicate block
# -----------------------------------------------------------------------------
cat("\n[PERMANOVA] Blocked by natman (within tree, by = 'margin')\n")
perm_block <- how(nperm = N_PERM_MAIN)
setBlocks(perm_block) <- meta$natman

perma_blocked <- list(
  robAit = adonis2(
    D_robAit ~ log_reads + umi + dw_type,
    data = meta,
    permutations = perm_block,
    by = "margin"
  ),
  jaccard = adonis2(
    D_jacc ~ log_reads + umi + dw_type,
    data = meta,
    permutations = perm_block,
    by = "margin"
  )
)

print(perma_blocked)
capture.output(perma_blocked, file = file.path(outdir, paste0("permanova_blocked_", scope, ".txt")))

# -----------------------------------------------------------------------------
# 7) Dispersion checks
# -----------------------------------------------------------------------------
cat("\n[Dispersion] betadisper by natman and dw_type\n")
disp_nat_ait  <- anova(betadisper(D_robAit, meta$natman))
disp_nat_jac  <- anova(betadisper(D_jacc,   meta$natman))
disp_type_ait <- anova(betadisper(D_robAit, meta$dw_type))
disp_type_jac <- anova(betadisper(D_jacc,   meta$dw_type))

print(list(
  Aitchison_natman = disp_nat_ait,
  Jaccard_natman   = disp_nat_jac,
  Aitchison_dwtype = disp_type_ait,
  Jaccard_dwtype   = disp_type_jac
))

sink(file.path(outdir, paste0("permdisp_", scope, ".txt")))
cat("Aitchison natman\n");  print(disp_nat_ait)
cat("Jaccard natman\n");    print(disp_nat_jac)
cat("Aitchison dw_type\n"); print(disp_type_ait)
cat("Jaccard dw_type\n");   print(disp_type_jac)
sink()

# -----------------------------------------------------------------------------
# 8) ANOSIM for natman
# -----------------------------------------------------------------------------
cat("\n[ANOSIM] natman grouping\n")
anos_nat_ait <- vegan::anosim(D_robAit, grouping = meta$natman, permutations = N_PERM_ANOSIM)
anos_nat_jac <- vegan::anosim(D_jacc,   grouping = meta$natman, permutations = N_PERM_ANOSIM)

print(anos_nat_ait)
print(anos_nat_jac)

anosim_tbl <- tibble::tibble(
  scope     = scope,
  metric    = c("robust.aitchison", "jaccard"),
  R_stat    = c(anos_nat_ait$statistic, anos_nat_jac$statistic),
  p_value   = c(anos_nat_ait$signif, anos_nat_jac$signif),
  n_perm    = N_PERM_ANOSIM,
  n_groups  = c(dplyr::n_distinct(meta$natman), dplyr::n_distinct(meta$natman)),
  n_samples = nrow(meta)
)

print(anosim_tbl)
readr::write_csv2(anosim_tbl, file.path(outdir, paste0("anosim_natman_", scope, ".csv")))

# -----------------------------------------------------------------------------
# 9) LOO PERMANOVA
# -----------------------------------------------------------------------------
cat("\n[LOO PERMANOVA] within tree effects, leave one natman out\n")

trees <- meta$natman %>%
  as.character() %>%
  unique() %>%
  setdiff(NA)

run_blocked_perma_subset <- function(D_full, meta_full, keep_idx, nperm) {
  M  <- as.matrix(D_full)
  Mk <- M[keep_idx, keep_idx, drop = FALSE]
  Dk <- stats::as.dist(Mk)
  meta_k <- meta_full[keep_idx, , drop = FALSE]
  
  permk <- permute::how(
    blocks = droplevels(meta_k$natman),
    nperm  = nperm
  )
  
  vegan::adonis2(
    Dk ~ log_reads + umi + dw_type,
    data = meta_k,
    permutations = permk,
    by = "margin"
  )
}

loo_list <- purrr::map_dfr(trees, function(tr) {
  keep <- meta$natman != tr & !is.na(meta$natman)
  if (sum(keep) < 6) return(NULL)
  
  m_ait <- run_blocked_perma_subset(D_robAit, meta, keep, N_PERM_LOO)
  dfA   <- as.data.frame(m_ait)
  dfA$term <- rownames(dfA)
  
  tabA <- dfA %>%
    tibble::as_tibble(rownames = "row") %>%
    dplyr::filter(term %in% c("log_reads", "umi", "dw_type")) %>%
    dplyr::transmute(
      tree_left_out = tr,
      metric        = "robust.aitchison",
      term,
      R2            = R2,
      p             = `Pr(>F)`
    )
  
  m_jac <- run_blocked_perma_subset(D_jacc, meta, keep, N_PERM_LOO)
  dfJ   <- as.data.frame(m_jac)
  dfJ$term <- rownames(dfJ)
  
  tabJ <- dfJ %>%
    tibble::as_tibble(rownames = "row") %>%
    dplyr::filter(term %in% c("log_reads", "umi", "dw_type")) %>%
    dplyr::transmute(
      tree_left_out = tr,
      metric        = "jaccard",
      term,
      R2            = R2,
      p             = `Pr(>F)`
    )
  
  dplyr::bind_rows(tabA, tabJ)
})

if (nrow(loo_list) > 0) {
  loo_summary <- loo_list %>%
    dplyr::group_by(metric, term) %>%
    dplyr::summarise(
      median_R2 = stats::median(R2, na.rm = TRUE),
      p05       = stats::quantile(p, 0.05, na.rm = TRUE),
      p95       = stats::quantile(p, 0.95, na.rm = TRUE),
      n_trees   = dplyr::n_distinct(tree_left_out),
      .groups   = "drop"
    )
  
  cat("\n[LOO PERMANOVA summary]\n")
  print(loo_summary)
  readr::write_csv2(loo_summary, file.path(outdir, paste0("loo_within_tree_permanova_", scope, ".csv")))
}

# -----------------------------------------------------------------------------
# 10) LOG–aFWD co-location test
# -----------------------------------------------------------------------------
cat("\n[LOG–aFWD co-location] cross type pairs only\n")

mk_log_afwd_pairs <- function(D, meta) {
  M <- as.matrix(D)
  stopifnot(identical(rownames(M), meta$sample))
  
  idx <- which(upper.tri(M), arr.ind = TRUE)
  i   <- idx[, 1]
  j   <- idx[, 2]
  
  dt_i <- as.character(meta$dw_type[i])
  dt_j <- as.character(meta$dw_type[j])
  tr_i <- as.character(meta$natman[i])
  tr_j <- as.character(meta$natman[j])
  
  is_cross <- (dt_i != dt_j) & dt_i %in% c("LOG", "aFWD") & dt_j %in% c("LOG", "aFWD")
  
  tibble(
    sample_i = meta$sample[i][is_cross],
    sample_j = meta$sample[j][is_cross],
    dw_type_i = dt_i[is_cross],
    dw_type_j = dt_j[is_cross],
    natman_i  = tr_i[is_cross],
    natman_j  = tr_j[is_cross],
    same_tree = factor(
      if_else(tr_i[is_cross] == tr_j[is_cross], "same_tree", "different_tree"),
      levels = c("same_tree", "different_tree")
    ),
    dist = M[cbind(i, j)][is_cross]
  )
}

auc_two_groups <- function(x, g) {
  g <- droplevels(g)
  stopifnot(nlevels(g) == 2L)
  
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  
  if (length(x1) == 0L || length(x2) == 0L) {
    return(c(AUC = NA_real_, r_rb = NA_real_))
  }
  
  r  <- rank(c(x1, x2))
  n1 <- length(x1)
  n2 <- length(x2)
  U  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  A  <- U / (n1 * n2)
  c(AUC = A, r_rb = 2 * A - 1)
}

anosim_two_groups <- function(x, g) {
  g <- droplevels(g)
  if (nlevels(g) != 2L) {
    return(NA_real_)
  }
  
  r <- rank(x, ties.method = "average")
  within_idx  <- g == levels(g)[1]
  between_idx <- g == levels(g)[2]
  
  if (!any(within_idx) || !any(between_idx)) {
    return(NA_real_)
  }
  
  r_within  <- mean(r[within_idx])
  r_between <- mean(r[between_idx])
  N <- length(r)
  
  (r_between - r_within) / (N / 2)
}

pairs_A <- mk_log_afwd_pairs(D_robAit, meta)
pairs_J <- mk_log_afwd_pairs(D_jacc,   meta)

cat("Number of cross type pairs by tree grouping (Aitchison):\n")
print(pairs_A %>% count(same_tree))

summ_pairs <- function(pairs_df, metric_label, D, meta, B_null) {
  sum_tab <- pairs_df %>%
    dplyr::group_by(same_tree) %>%
    dplyr::summarise(
      n_pairs     = dplyr::n(),
      median_dist = stats::median(dist),
      IQR_dist    = IQR(dist),
      .groups     = "drop"
    )
  
  eff <- auc_two_groups(pairs_df$dist, pairs_df$same_tree)
  A_obs_raw <- eff["AUC"]
  r_raw     <- eff["r_rb"]
  
  A_flipped <- 1 - A_obs_raw
  r_flipped <- 2 * A_flipped - 1
  
  R_obs <- anosim_two_groups(pairs_df$dist, pairs_df$same_tree)
  
  nat <- droplevels(meta$natman)
  M   <- as.matrix(D)
  stopifnot(identical(rownames(M), meta$sample))
  
  A_null <- numeric(B_null)
  R_null <- numeric(B_null)
  
  idx_i <- match(pairs_df$sample_i, meta$sample)
  idx_j <- match(pairs_df$sample_j, meta$sample)
  
  for (b in seq_len(B_null)) {
    nat_perm <- factor(sample(nat, replace = FALSE), levels = levels(nat))
    g_perm <- factor(
      if_else(nat_perm[idx_i] == nat_perm[idx_j], "same_tree", "different_tree"),
      levels = c("same_tree", "different_tree")
    )
    
    eff_b     <- auc_two_groups(pairs_df$dist, g_perm)
    A_null[b] <- 1 - eff_b["AUC"]
    R_null[b] <- anosim_two_groups(pairs_df$dist, g_perm)
  }
  
  ci_A_null <- stats::quantile(A_null, c(0.025, 0.975), na.rm = TRUE)
  p_A_emp   <- mean(A_null >= A_flipped, na.rm = TRUE)
  ses_A     <- (A_flipped - mean(A_null, na.rm = TRUE)) / stats::sd(A_null, na.rm = TRUE)
  
  auc_df <- tibble::tibble(AUC_null = A_null)
  p_auc <- ggplot2::ggplot(auc_df, ggplot2::aes(x = AUC_null)) +
    ggplot2::geom_density(fill = "grey80") +
    ggplot2::geom_vline(xintercept = A_flipped, linetype = 2) +
    ggplot2::labs(
      title = paste0("AUC null vs observed - ", metric_label, " - ", scope),
      subtitle = sprintf(
        "AUC_obs_flipped = %.3f  null_mean = %.3f  95%%CI = [%.3f, %.3f]  SES = %.2f",
        A_flipped, mean(A_null, na.rm = TRUE), ci_A_null[1], ci_A_null[2], ses_A
      ),
      x = "AUC_flipped (same_tree closer)",
      y = "Density"
    ) +
    ggplot2::theme_bw()
  
  print(p_auc)
  ggplot2::ggsave(
    file.path(outdir, paste0("AUC_null_", metric_label, "_", scope, ".png")),
    p_auc, width = 6, height = 4, dpi = 300
  )
  
  ci_R_null <- stats::quantile(R_null, c(0.025, 0.975), na.rm = TRUE)
  p_R_emp   <- mean(R_null >= R_obs, na.rm = TRUE)
  ses_R     <- (R_obs - mean(R_null, na.rm = TRUE)) / stats::sd(R_null, na.rm = TRUE)
  
  tibble::tibble(
    metric              = metric_label,
    n_same              = sum_tab$n_pairs[sum_tab$same_tree == "same_tree"],
    n_diff              = sum_tab$n_pairs[sum_tab$same_tree == "different_tree"],
    median_same         = sum_tab$median_dist[sum_tab$same_tree == "same_tree"],
    median_diff         = sum_tab$median_dist[sum_tab$same_tree == "different_tree"],
    IQR_same            = sum_tab$IQR_dist[sum_tab$same_tree == "same_tree"],
    IQR_diff            = sum_tab$IQR_dist[sum_tab$same_tree == "different_tree"],
    AUC_flipped         = unname(A_flipped),
    r_rank_biserial     = unname(r_flipped),
    AUC_null_mean       = mean(A_null, na.rm = TRUE),
    AUC_null_lo         = ci_A_null[1],
    AUC_null_hi         = ci_A_null[2],
    AUC_p_emp           = p_A_emp,
    AUC_SES             = ses_A,
    ANOSIM_R_obs        = unname(R_obs),
    ANOSIM_R_null_mean  = mean(R_null, na.rm = TRUE),
    ANOSIM_R_null_lo    = ci_R_null[1],
    ANOSIM_R_null_hi    = ci_R_null[2],
    ANOSIM_R_p_emp      = p_R_emp,
    ANOSIM_R_SES        = ses_R
  )
}

co_loc_tbl <- bind_rows(
  summ_pairs(pairs_A, "robust.aitchison", D_robAit, meta, B_NULL_AUC),
  summ_pairs(pairs_J, "jaccard", D_jacc, meta, B_NULL_AUC)
)

cat("\n[LOG–aFWD co-location summary]\n")
print(co_loc_tbl)

write_csv2(co_loc_tbl, file.path(outdir, paste0("LOG_aFWD_colocation_summary_", scope, ".csv")))

# -----------------------------------------------------------------------------
# 11) ECDF plots
# -----------------------------------------------------------------------------
ecdf_dir <- file.path(outdir, "ecdf")
ensure_dir(ecdf_dir)

ecdf_plot_A <- ggplot(pairs_A, aes(x = dist, colour = same_tree)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  theme_bw(base_size = 11) +
  labs(
    title  = paste0("LOG–aFWD co-location (robust Aitchison) ", scope),
    x      = "Pairwise distance",
    y      = "ECDF",
    colour = "Pair type"
  )

print(ecdf_plot_A)

ecdf_plot_J <- ggplot(pairs_J, aes(x = dist, colour = same_tree)) +
  stat_ecdf(geom = "step", linewidth = 1) +
  theme_bw(base_size = 11) +
  labs(
    title  = paste0("LOG–aFWD co-location (Jaccard) ", scope),
    x      = "Pairwise distance",
    y      = "ECDF",
    colour = "Pair type"
  )

print(ecdf_plot_J)

ggsave(file.path(ecdf_dir, paste0("ECDF_LOG_aFWD_Aitchison_", scope, ".png")), ecdf_plot_A, width = 6.2, height = 4.2, dpi = 300)
ggsave(file.path(ecdf_dir, paste0("ECDF_LOG_aFWD_Jaccard_", scope, ".png")), ecdf_plot_J, width = 6.2, height = 4.2, dpi = 300)

# -----------------------------------------------------------------------------
# 12) BETAPART LOG vs aFWD
# -----------------------------------------------------------------------------
cat("\n[BETAPART] beta.sim, beta.sne, beta.sor for LOG vs aFWD (within vs between tree)\n")
stopifnot(identical(rownames(otu_pa), meta$sample))

bp <- betapart::beta.pair(otu_pa, index.family = "sorensen")
DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

mk_log_afwd_betapart_pairs <- function(D_sim, D_sne, D_sor, meta) {
  stopifnot(
    identical(rownames(D_sim), meta$sample),
    identical(rownames(D_sne), meta$sample),
    identical(rownames(D_sor), meta$sample)
  )
  
  idx <- which(upper.tri(D_sim), arr.ind = TRUE)
  i   <- idx[, 1]
  j   <- idx[, 2]
  
  state_i <- as.character(meta$dw_type[i])
  state_j <- as.character(meta$dw_type[j])
  nat_i   <- as.character(meta$natman[i])
  nat_j   <- as.character(meta$natman[j])
  
  is_cross_state <- (state_i != state_j) &
    state_i %in% c("LOG", "aFWD") &
    state_j %in% c("LOG", "aFWD")
  
  if (!any(is_cross_state)) {
    return(tibble::tibble())
  }
  
  tibble::tibble(
    sample_i   = meta$sample[i][is_cross_state],
    sample_j   = meta$sample[j][is_cross_state],
    state_i    = state_i[is_cross_state],
    state_j    = state_j[is_cross_state],
    natman_i   = nat_i[is_cross_state],
    natman_j   = nat_j[is_cross_state],
    pair_class = factor(
      ifelse(nat_i[is_cross_state] == nat_j[is_cross_state], "within_tree", "between_tree"),
      levels = c("within_tree", "between_tree")
    ),
    beta_sim   = D_sim[cbind(i, j)][is_cross_state],
    beta_sne   = D_sne[cbind(i, j)][is_cross_state],
    beta_sor   = D_sor[cbind(i, j)][is_cross_state]
  )
}

pairs_all <- mk_log_afwd_betapart_pairs(DM_sim, DM_sne, DM_sor, meta)

if (nrow(pairs_all) == 0L) {
  cat("\n[BETAPART] No LOG vs aFWD cross state pairs, skipping\n")
} else {
  cat("\n[BETAPART LOG vs aFWD within vs between tree summary]\n")
  print(pairs_all %>% dplyr::count(pair_class, name = "n_pairs"))
  
  beta_state_within_between <- pairs_all %>%
    dplyr::group_by(pair_class) %>%
    dplyr::summarise(
      contrast   = "LOG_vs_aFWD_cross_state",
      n_pairs    = dplyr::n(),
      median_sim = stats::median(beta_sim, na.rm = TRUE),
      median_sne = stats::median(beta_sne, na.rm = TRUE),
      median_sor = stats::median(beta_sor, na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    dplyr::mutate(
      scope = scope,
      T_to_N_median = median_sim / median_sne
    ) %>%
    dplyr::select(
      scope, contrast, pair_class, n_pairs,
      median_sim, median_sne, median_sor, T_to_N_median
    )
  
  print(beta_state_within_between)
  
  readr::write_csv2(
    beta_state_within_between,
    file.path(outdir, paste0("betapart_LOG_vs_aFWD_within_between_", scope, ".csv"))
  )
}
