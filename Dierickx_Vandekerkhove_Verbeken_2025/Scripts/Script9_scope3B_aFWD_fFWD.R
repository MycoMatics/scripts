# =============================================================================
# Script: Script9_scope3B_aFWD_fFWD.R
# Project: Deadwood fungal community analyses
# Scope: aFWD_vs_fFWD
# Purpose: Compare attached and fallen fine woody debris at global and paired
#          tree-linked scales using PERMANOVA, null models, AUC/ECDF, ANOSIM,
#          PERMDISP, LOO stability, and betapart decomposition
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - d_ait_or_hel
# - fix_dist_labels
# - mk_pair_tbl
# - pairs_by
# - median_diff_boot
# - auc_once
# - upper_idx_global
# - auc_from_group
# - get_R2_term
# - state_R2_null_blocked
# - shuffle_state_within_tree
# - ensure_dir
# - assert_objects
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(betapart)
  library(readr)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c(
  "d_ait_or_hel", "fix_dist_labels", "mk_pair_tbl", "pairs_by",
  "median_diff_boot", "auc_once", "upper_idx_global", "auc_from_group",
  "get_R2_term", "state_R2_null_blocked", "shuffle_state_within_tree",
  "ensure_dir"
))

set.seed(SEED_GLOBAL)
setwd(".")

cfg <- list(
  nperm  = 999L,
  B_R2   = 999L,
  B_AUC  = 5000L,
  B_boot = 5000L
)

outdir <- "plots/FWD/aFWD_vs_fFWD"
ensure_dir(outdir)
scope <- "aFWD_vs_fFWD"

# =============================================================================
# 1) Data preparation for FWD
# =============================================================================
meta_fwd <- META1 |>
  filter(dw_type2 %in% c("aFWD", "fFWD")) |>
  mutate(
    fwd_state = case_when(
      dw_type2 == "aFWD" ~ "ATTACHED",
      dw_type2 == "fFWD" ~ "FALLEN"
    ) |> factor(levels = c("ATTACHED", "FALLEN")),
    natman      = droplevels(as.factor(natman)),
    decay_stage = droplevels(as.factor(decay_stage)),
    diameter_at_drill = suppressWarnings(as.numeric(diameter_at_drill))
  )
otu_fwd <- otu_matrix_filt[meta_fwd$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fwd), meta_fwd$sample))
otu_fwd <- otu_fwd[, colSums(otu_fwd) > 0, drop = FALSE]

meta_fwd <- meta_fwd |>
  mutate(
    reads     = rowSums(otu_fwd),
    log_reads = log1p(reads)
  ) |>
  filter(reads > 0) |>
  droplevels()

otu_fwd <- otu_fwd[meta_fwd$sample, , drop = FALSE]
otu_fwd <- otu_fwd[, colSums(otu_fwd, na.rm = TRUE) > 0, drop = FALSE]
stopifnot(identical(rownames(otu_fwd), meta_fwd$sample))

cat("GLOBAL FWD data\n")
print(meta_fwd |> count(fwd_state))
cat(
  "Samples:", nrow(meta_fwd),
  " OTUs:", ncol(otu_fwd),
  " trees:", length(unique(meta_fwd$natman)), "\n"
)

# =============================================================================
# 2) GLOBAL analysis
# =============================================================================
D_Ait_global <- d_ait_or_hel(otu_fwd)
D_Ait_global <- fix_dist_labels(D_Ait_global, samps = rownames(otu_fwd))
D_Jac_global <- vegdist((otu_fwd > 0) * 1, method = "jaccard", binary = TRUE)
D_Jac_global <- fix_dist_labels(D_Jac_global, samps = rownames(otu_fwd))

perm_unblk <- how(nperm = cfg$nperm)

ado_Ait_global <- adonis2(
  D_Ait_global ~ log_reads + umi + natman + fwd_state,
  data = meta_fwd,
  by = "terms",
  permutations = perm_unblk
)
ado_Jac_global <- adonis2(
  D_Jac_global ~ log_reads + umi + natman + fwd_state,
  data = meta_fwd,
  by = "terms",
  permutations = perm_unblk
)

cat("\n[GLOBAL] PERMANOVA Aitchison\n")
print(ado_Ait_global)
cat("\n[GLOBAL] PERMANOVA Jaccard\n")
print(ado_Jac_global)

R2_Ait_state_global <- as.data.frame(ado_Ait_global)["fwd_state", "R2"]
R2_Jac_state_global <- as.data.frame(ado_Jac_global)["fwd_state", "R2"]

pairs_global_A <- mk_pair_tbl(D_Ait_global, meta_fwd) |> pairs_by("state_global")
pairs_global_J <- mk_pair_tbl(D_Jac_global, meta_fwd) |> pairs_by("state_global")

eff_global_A <- median_diff_boot(pairs_global_A)
eff_global_J <- median_diff_boot(pairs_global_J)
auc_global_A <- auc_once(pairs_global_A)
auc_global_J <- auc_once(pairs_global_J)

summary_global <- tibble(
  scope      = "GLOBAL_unblocked",
  metric     = c("robust.aitchison", "jaccard"),
  n_samples  = nrow(meta_fwd),
  R2_state   = c(R2_Ait_state_global, R2_Jac_state_global),
  AUC_state  = c(auc_global_A$AUC, auc_global_J$AUC),
  r_rb_state = c(auc_global_A$r_rb, auc_global_J$r_rb),
  diff_median_state = c(eff_global_A$diff_median, eff_global_J$diff_median),
  diff_CI_lo_state  = c(eff_global_A$diff_CI_lo, eff_global_J$diff_CI_lo),
  diff_CI_hi_state  = c(eff_global_A$diff_CI_hi, eff_global_J$diff_CI_hi),
  n_pairs_within    = c(eff_global_A$n_within, eff_global_J$n_within),
  n_pairs_between   = c(eff_global_A$n_between, eff_global_J$n_between)
)

print(summary_global)
write_csv2(summary_global, file.path(outdir, "GLOBAL_aFWD_vs_fFWD_summary.csv"))

# =============================================================================
# 3) GLOBAL natman and state nulls + ANOSIM + PERMDISP
# =============================================================================
cat("\n[GLOBAL] natman and state nulls (R2, AUC) + ANOSIM\n")

rhs_unblocked <- c("log_reads", "umi", "natman", "fwd_state")
perm_unblk_full <- how(nperm = cfg$nperm)

perma_unblocked <- list(
  robAit = adonis2(
    D_Ait_global ~ .,
    data = meta_fwd[, rhs_unblocked, drop = FALSE],
    by = "terms",
    permutations = perm_unblk_full
  ),
  jaccard = adonis2(
    D_Jac_global ~ .,
    data = meta_fwd[, rhs_unblocked, drop = FALSE],
    by = "terms",
    permutations = perm_unblk_full
  )
)

cat("\n[GLOBAL] PERMANOVA (unblocked, natman + state)\n")
print(perma_unblocked)

# ---- PERMDISP GLOBAL --------------------------------------------------------
cat("\n[Dispersion GLOBAL] betadisper by natman and fwd_state\n")

disp_nat_Ait   <- anova(betadisper(D_Ait_global, meta_fwd$natman))
disp_nat_Jac   <- anova(betadisper(D_Jac_global, meta_fwd$natman))
disp_state_Ait <- anova(betadisper(D_Ait_global, meta_fwd$fwd_state))
disp_state_Jac <- anova(betadisper(D_Jac_global, meta_fwd$fwd_state))

print(list(
  Aitchison_natman = disp_nat_Ait,
  Jaccard_natman   = disp_nat_Jac,
  Aitchison_state  = disp_state_Ait,
  Jaccard_state    = disp_state_Jac
))

sink(file.path(outdir, paste0("permdisp_GLOBAL_", scope, ".txt")))
cat("Aitchison natman\n");    print(disp_nat_Ait)
cat("Jaccard natman\n");      print(disp_nat_Jac)
cat("Aitchison fwd_state\n"); print(disp_state_Ait)
cat("Jaccard fwd_state\n");   print(disp_state_Jac)
sink()

nat      <- droplevels(meta_fwd$natman)
state    <- droplevels(meta_fwd$fwd_state)
n_natman <- nlevels(nat)
n_state  <- nlevels(state)
n_samp   <- nrow(meta_fwd)

upper_Ait <- upper_idx_global(D_Ait_global)
upper_Jac <- upper_idx_global(D_Jac_global)

R2_nat_obs_Ait   <- get_R2_term(perma_unblocked$robAit, "natman")
R2_nat_obs_Jac   <- get_R2_term(perma_unblocked$jaccard, "natman")
R2_state_obs_Ait <- get_R2_term(perma_unblocked$robAit, "fwd_state")
R2_state_obs_Jac <- get_R2_term(perma_unblocked$jaccard, "fwd_state")

AUC_nat_obs_A   <- auc_from_group(D_Ait_global, nat, upper_Ait)
AUC_nat_obs_J   <- auc_from_group(D_Jac_global, nat, upper_Jac)
AUC_state_obs_A <- auc_from_group(D_Ait_global, state, upper_Ait)
AUC_state_obs_J <- auc_from_group(D_Jac_global, state, upper_Jac)

N_PERM_ANOSIM <- 9999L

anos_nat_Ait <- vegan::anosim(D_Ait_global, grouping = nat, permutations = N_PERM_ANOSIM)
anos_nat_Jac <- vegan::anosim(D_Jac_global, grouping = nat, permutations = N_PERM_ANOSIM)

anos_state_Ait <- vegan::anosim(D_Ait_global, grouping = state, permutations = N_PERM_ANOSIM)
anos_state_Jac <- vegan::anosim(D_Jac_global, grouping = state, permutations = N_PERM_ANOSIM)

cat("\n[GLOBAL ANOSIM] natman\n")
print(anos_nat_Ait)
print(anos_nat_Jac)

cat("\n[GLOBAL ANOSIM] state (fwd_state)\n")
print(anos_state_Ait)
print(anos_state_Jac)

# ---- natman nulls -----------------------------------------------------------
B_R2_nat       <- cfg$B_R2
B_AUC_nat      <- cfg$B_AUC
N_PERM_NULL_R2 <- 199L

R2_nat_null_Ait <- numeric(B_R2_nat)
R2_nat_null_Jac <- numeric(B_R2_nat)
AUC_nat_null_A  <- numeric(B_AUC_nat)
AUC_nat_null_J  <- numeric(B_AUC_nat)

df_terms <- meta_fwd[, rhs_unblocked, drop = FALSE]

message("[Nulls] natman R2 nulls …")
for (b in seq_len(B_R2_nat)) {
  nat_null <- factor(sample(nat, replace = FALSE), levels = levels(nat))
  df_tmp <- df_terms
  df_tmp$natman <- nat_null
  
  fit_A <- adonis2(
    D_Ait_global ~ .,
    data = df_tmp,
    by = "terms",
    permutations = N_PERM_NULL_R2
  )
  fit_J <- adonis2(
    D_Jac_global ~ .,
    data = df_tmp,
    by = "terms",
    permutations = N_PERM_NULL_R2
  )
  
  R2_nat_null_Ait[b] <- get_R2_term(fit_A, "natman")
  R2_nat_null_Jac[b] <- get_R2_term(fit_J, "natman")
  
  if (b %% 50L == 0L) {
    message("[natman R2 nulls] finished ", b, " of ", B_R2_nat)
  }
}

message("[Nulls] natman AUC nulls …")
for (b in seq_len(B_AUC_nat)) {
  nat_null <- factor(sample(nat, replace = FALSE), levels = levels(nat))
  AUC_nat_null_A[b] <- auc_from_group(D_Ait_global, nat_null, upper_Ait)["AUC"]
  AUC_nat_null_J[b] <- auc_from_group(D_Jac_global, nat_null, upper_Jac)["AUC"]
  
  if (b %% 500L == 0L) {
    message("[natman AUC nulls] finished ", b, " of ", B_AUC_nat)
  }
}

natman_R2_null <- cbind(
  robust_aitchison = R2_nat_null_Ait,
  jaccard          = R2_nat_null_Jac
)
natman_AUC_null <- cbind(
  robust_aitchison = AUC_nat_null_A,
  jaccard          = AUC_nat_null_J
)

saveRDS(natman_R2_null,  file = file.path(outdir, "NATMAN_NULL_R2.rds"))
saveRDS(natman_AUC_null, file = file.path(outdir, "NATMAN_NULL_AUC.rds"))

natman_evidence <- tibble(
  scope          = scope,
  metric         = c("robust.aitchison", "jaccard"),
  R2_obs         = c(R2_nat_obs_Ait, R2_nat_obs_Jac),
  R2_null_mean   = c(mean(R2_nat_null_Ait, na.rm = TRUE), mean(R2_nat_null_Jac, na.rm = TRUE)),
  R2_null_CI_lo  = c(quantile(R2_nat_null_Ait, 0.025, na.rm = TRUE), quantile(R2_nat_null_Jac, 0.025, na.rm = TRUE)),
  R2_null_CI_hi  = c(quantile(R2_nat_null_Ait, 0.975, na.rm = TRUE), quantile(R2_nat_null_Jac, 0.975, na.rm = TRUE)),
  p_R2_emp       = c(mean(R2_nat_null_Ait >= R2_nat_obs_Ait, na.rm = TRUE), mean(R2_nat_null_Jac >= R2_nat_obs_Jac, na.rm = TRUE)),
  R2_SES         = c(
    (R2_nat_obs_Ait - mean(R2_nat_null_Ait, na.rm = TRUE)) / sd(R2_nat_null_Ait, na.rm = TRUE),
    (R2_nat_obs_Jac - mean(R2_nat_null_Jac, na.rm = TRUE)) / sd(R2_nat_null_Jac, na.rm = TRUE)
  ),
  R2_percentile  = c(mean(R2_nat_null_Ait <= R2_nat_obs_Ait, na.rm = TRUE), mean(R2_nat_null_Jac <= R2_nat_obs_Jac, na.rm = TRUE)),
  AUC_obs        = c(AUC_nat_obs_A["AUC"], AUC_nat_obs_J["AUC"]),
  rank_biserial  = c(AUC_nat_obs_A["r_rb"], AUC_nat_obs_J["r_rb"]),
  AUC_null_mean  = c(mean(AUC_nat_null_A, na.rm = TRUE), mean(AUC_nat_null_J, na.rm = TRUE)),
  AUC_null_CI_lo = c(quantile(AUC_nat_null_A, 0.025, na.rm = TRUE), quantile(AUC_nat_null_J, 0.025, na.rm = TRUE)),
  AUC_null_CI_hi = c(quantile(AUC_nat_null_A, 0.975, na.rm = TRUE), quantile(AUC_nat_null_J, 0.975, na.rm = TRUE)),
  AUC_p_emp      = c(mean(AUC_nat_null_A >= AUC_nat_obs_A["AUC"], na.rm = TRUE), mean(AUC_nat_null_J >= AUC_nat_obs_J["AUC"], na.rm = TRUE)),
  AUC_SES        = c(
    (AUC_nat_obs_A["AUC"] - mean(AUC_nat_null_A, na.rm = TRUE)) / sd(AUC_nat_null_A, na.rm = TRUE),
    (AUC_nat_obs_J["AUC"] - mean(AUC_nat_null_J, na.rm = TRUE)) / sd(AUC_nat_null_J, na.rm = TRUE)
  ),
  ANOSIM_R        = c(anos_nat_Ait$statistic, anos_nat_Jac$statistic),
  ANOSIM_p        = c(anos_nat_Ait$signif, anos_nat_Jac$signif),
  n_natman_levels = n_natman,
  n_samples       = n_samp
)

cat("\n[GLOBAL] natman evidence summary\n")
print(natman_evidence)

write_csv2(natman_evidence, file.path(outdir, "GLOBAL_natman_evidence_R2_AUC_ANOSIM.csv"))

# ---- natman null plots ------------------------------------------------------
out_nat <- file.path(outdir, "natman")
ensure_dir(out_nat)

ciR2_nat_Ait <- quantile(R2_nat_null_Ait, c(0.025, 0.975), na.rm = TRUE)
ciR2_nat_Jac <- quantile(R2_nat_null_Jac, c(0.025, 0.975), na.rm = TRUE)
pR2_nat_Ait  <- mean(R2_nat_null_Ait >= R2_nat_obs_Ait, na.rm = TRUE)
pR2_nat_Jac  <- mean(R2_nat_null_Jac >= R2_nat_obs_Jac, na.rm = TRUE)

sesR2_nat_Ait <- (R2_nat_obs_Ait - mean(R2_nat_null_Ait, na.rm = TRUE)) / sd(R2_nat_null_Ait, na.rm = TRUE)
sesR2_nat_Jac <- (R2_nat_obs_Jac - mean(R2_nat_null_Jac, na.rm = TRUE)) / sd(R2_nat_null_Jac, na.rm = TRUE)

pR2_nat_Ait_plot <- ggplot(tibble(R2 = R2_nat_null_Ait), aes(x = R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_nat_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciR2_nat_Ait, linetype = 3) +
  labs(
    title    = paste0(scope, " natman R2 null vs obs (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_nat_obs_Ait, mean(R2_nat_null_Ait, na.rm = TRUE), ciR2_nat_Ait[1], ciR2_nat_Ait[2], pR2_nat_Ait, sesR2_nat_Ait
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

print(pR2_nat_Ait_plot)
ggsave(file.path(out_nat, paste0("natman_R2_null_vs_obs_Aitchison_", scope, ".png")), pR2_nat_Ait_plot, width = 6, height = 4, dpi = 300)

pR2_nat_Jac_plot <- ggplot(tibble(R2 = R2_nat_null_Jac), aes(x = R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_nat_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciR2_nat_Jac, linetype = 3) +
  labs(
    title    = paste0(scope, " natman R2 null vs obs (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_nat_obs_Jac, mean(R2_nat_null_Jac, na.rm = TRUE), ciR2_nat_Jac[1], ciR2_nat_Jac[2], pR2_nat_Jac, sesR2_nat_Jac
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

print(pR2_nat_Jac_plot)
ggsave(file.path(out_nat, paste0("natman_R2_null_vs_obs_Jaccard_", scope, ".png")), pR2_nat_Jac_plot, width = 6, height = 4, dpi = 300)

AUC_nat_obs_Ait <- as.numeric(AUC_nat_obs_A[["AUC"]])
AUC_nat_obs_Jac <- as.numeric(AUC_nat_obs_J[["AUC"]])

ciAUC_nat_Ait <- quantile(AUC_nat_null_A, c(0.025, 0.975), na.rm = TRUE)
ciAUC_nat_Jac <- quantile(AUC_nat_null_J, c(0.025, 0.975), na.rm = TRUE)
pAUC_nat_Ait  <- mean(AUC_nat_null_A >= AUC_nat_obs_Ait, na.rm = TRUE)
pAUC_nat_Jac  <- mean(AUC_nat_null_J >= AUC_nat_obs_Jac, na.rm = TRUE)

sesAUC_nat_Ait <- (AUC_nat_obs_Ait - mean(AUC_nat_null_A, na.rm = TRUE)) / sd(AUC_nat_null_A, na.rm = TRUE)
sesAUC_nat_Jac <- (AUC_nat_obs_Jac - mean(AUC_nat_null_J, na.rm = TRUE)) / sd(AUC_nat_null_J, na.rm = TRUE)

pAUC_nat_Ait_plot <- ggplot(tibble(AUC = AUC_nat_null_A), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_nat_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_nat_Ait, linetype = 3) +
  labs(
    title    = paste0(scope, " natman AUC null vs observed (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_nat_obs_Ait, mean(AUC_nat_null_A, na.rm = TRUE), pAUC_nat_Ait, sesAUC_nat_Ait
    ),
    x = "AUC",
    y = "Density"
  ) +
  theme_bw()

print(pAUC_nat_Ait_plot)
ggsave(file.path(out_nat, paste0("natman_AUC_null_vs_obs_Aitchison_", scope, ".png")), pAUC_nat_Ait_plot, width = 6, height = 4, dpi = 300)

pAUC_nat_Jac_plot <- ggplot(tibble(AUC = AUC_nat_null_J), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_nat_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_nat_Jac, linetype = 3) +
  labs(
    title    = paste0(scope, " natman AUC null vs observed (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_nat_obs_Jac, mean(AUC_nat_null_J, na.rm = TRUE), pAUC_nat_Jac, sesAUC_nat_Jac
    ),
    x = "AUC",
    y = "Density"
  ) +
  theme_bw()

print(pAUC_nat_Jac_plot)
ggsave(file.path(out_nat, paste0("natman_AUC_null_vs_obs_Jaccard_", scope, ".png")), pAUC_nat_Jac_plot, width = 6, height = 4, dpi = 300)

# ---- ECDF natman ------------------------------------------------------------
pairwise_from_D <- function(D, grp, upper_idx) {
  M <- as.matrix(D)
  dU <- M[cbind(upper_idx[, 1], upper_idx[, 2])]
  same <- grp[upper_idx[, 1]] == grp[upper_idx[, 2]]
  tibble::tibble(
    dist = dU,
    type = factor(ifelse(same, "within", "between"), levels = c("within", "between"))
  )
}

df_nat_Ait <- pairwise_from_D(D_Ait_global, nat, upper_Ait)
p_nat_ECDF_Ait <- ggplot2::ggplot(df_nat_Ait, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " natman ECDF of pairwise distances (Aitchison)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_nat_obs_A[["AUC"]], AUC_nat_obs_A[["r_rb"]],
      sum(df_nat_Ait$type == "within"), sum(df_nat_Ait$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_nat_ECDF_Ait)
ggplot2::ggsave(file.path(out_nat, paste0("natman_ECDF_Aitchison_", scope, ".png")), p_nat_ECDF_Ait, width = 6.2, height = 4.2, dpi = 300)

df_nat_Jac <- pairwise_from_D(D_Jac_global, nat, upper_Jac)
p_nat_ECDF_Jac <- ggplot2::ggplot(df_nat_Jac, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " natman ECDF of pairwise distances (Jaccard)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_nat_obs_J[["AUC"]], AUC_nat_obs_J[["r_rb"]],
      sum(df_nat_Jac$type == "within"), sum(df_nat_Jac$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_nat_ECDF_Jac)
ggplot2::ggsave(file.path(out_nat, paste0("natman_ECDF_Jaccard_", scope, ".png")), p_nat_ECDF_Jac, width = 6.2, height = 4.2, dpi = 300)

# ---- state nulls ------------------------------------------------------------
B_R2_state <- cfg$B_R2
B_AUC_state <- cfg$B_AUC

R2_state_null_Ait <- numeric(B_R2_state)
R2_state_null_Jac <- numeric(B_R2_state)
AUC_state_null_A  <- numeric(B_AUC_state)
AUC_state_null_J  <- numeric(B_AUC_state)

message("[Nulls] state R2 nulls …")
for (b in seq_len(B_R2_state)) {
  state_null <- factor(sample(state, replace = FALSE), levels = levels(state))
  df_tmp <- df_terms
  df_tmp$fwd_state <- state_null
  
  fit_A <- adonis2(
    D_Ait_global ~ .,
    data = df_tmp,
    by = "terms",
    permutations = N_PERM_NULL_R2
  )
  fit_J <- adonis2(
    D_Jac_global ~ .,
    data = df_tmp,
    by = "terms",
    permutations = N_PERM_NULL_R2
  )
  
  R2_state_null_Ait[b] <- get_R2_term(fit_A, "fwd_state")
  R2_state_null_Jac[b] <- get_R2_term(fit_J, "fwd_state")
  
  if (b %% 50L == 0L) {
    message("[state R2 nulls] finished ", b, " of ", B_R2_state)
  }
}

message("[Nulls] state AUC nulls …")
for (b in seq_len(B_AUC_state)) {
  state_null <- factor(sample(state, replace = FALSE), levels = levels(state))
  AUC_state_null_A[b] <- auc_from_group(D_Ait_global, state_null, upper_Ait)["AUC"]
  AUC_state_null_J[b] <- auc_from_group(D_Jac_global, state_null, upper_Jac)["AUC"]
  
  if (b %% 500L == 0L) {
    message("[state AUC nulls] finished ", b, " of ", B_AUC_state)
  }
}

state_R2_null <- cbind(
  robust_aitchison = R2_state_null_Ait,
  jaccard          = R2_state_null_Jac
)
state_AUC_null <- cbind(
  robust_aitchison = AUC_state_null_A,
  jaccard          = AUC_state_null_J
)

saveRDS(state_R2_null,  file = file.path(outdir, "STATE_NULL_R2.rds"))
saveRDS(state_AUC_null, file = file.path(outdir, "STATE_NULL_AUC.rds"))

state_evidence <- tibble(
  scope          = scope,
  metric         = c("robust.aitchison", "jaccard"),
  R2_obs         = c(R2_state_obs_Ait, R2_state_obs_Jac),
  R2_null_mean   = c(mean(R2_state_null_Ait, na.rm = TRUE), mean(R2_state_null_Jac, na.rm = TRUE)),
  R2_null_CI_lo  = c(quantile(R2_state_null_Ait, 0.025, na.rm = TRUE), quantile(R2_state_null_Jac, 0.025, na.rm = TRUE)),
  R2_null_CI_hi  = c(quantile(R2_state_null_Ait, 0.975, na.rm = TRUE), quantile(R2_state_null_Jac, 0.975, na.rm = TRUE)),
  p_R2_emp       = c(mean(R2_state_null_Ait >= R2_state_obs_Ait, na.rm = TRUE), mean(R2_state_null_Jac >= R2_state_obs_Jac, na.rm = TRUE)),
  R2_SES         = c(
    (R2_state_obs_Ait - mean(R2_state_null_Ait, na.rm = TRUE)) / sd(R2_state_null_Ait, na.rm = TRUE),
    (R2_state_obs_Jac - mean(R2_state_null_Jac, na.rm = TRUE)) / sd(R2_state_null_Jac, na.rm = TRUE)
  ),
  R2_percentile  = c(mean(R2_state_null_Ait <= R2_state_obs_Ait, na.rm = TRUE), mean(R2_state_null_Jac <= R2_state_obs_Jac, na.rm = TRUE)),
  AUC_obs        = c(AUC_state_obs_A["AUC"], AUC_state_obs_J["AUC"]),
  rank_biserial  = c(AUC_state_obs_A["r_rb"], AUC_state_obs_J["r_rb"]),
  AUC_null_mean  = c(mean(AUC_state_null_A, na.rm = TRUE), mean(AUC_state_null_J, na.rm = TRUE)),
  AUC_null_CI_lo = c(quantile(AUC_state_null_A, 0.025, na.rm = TRUE), quantile(AUC_state_null_J, 0.025, na.rm = TRUE)),
  AUC_null_CI_hi = c(quantile(AUC_state_null_A, 0.975, na.rm = TRUE), quantile(AUC_state_null_J, 0.975, na.rm = TRUE)),
  AUC_p_emp      = c(mean(AUC_state_null_A >= AUC_state_obs_A["AUC"], na.rm = TRUE), mean(AUC_state_null_J >= AUC_state_obs_J["AUC"], na.rm = TRUE)),
  AUC_SES        = c(
    (AUC_state_obs_A["AUC"] - mean(AUC_state_null_A, na.rm = TRUE)) / sd(AUC_state_null_A, na.rm = TRUE),
    (AUC_state_obs_J["AUC"] - mean(AUC_state_null_J, na.rm = TRUE)) / sd(AUC_state_null_J, na.rm = TRUE)
  ),
  ANOSIM_R       = c(anos_state_Ait$statistic, anos_state_Jac$statistic),
  ANOSIM_p       = c(anos_state_Ait$signif, anos_state_Jac$signif),
  n_state_levels = n_state,
  n_samples      = n_samp
)

cat("\n[GLOBAL] state evidence summary\n")
print(state_evidence)
write_csv2(state_evidence, file.path(outdir, "GLOBAL_state_evidence_R2_AUC_ANOSIM.csv"))

# ---- state null plots -------------------------------------------------------
out_state <- file.path(outdir, "state")
ensure_dir(out_state)

ciR2_state_Ait <- quantile(R2_state_null_Ait, c(0.025, 0.975), na.rm = TRUE)
ciR2_state_Jac <- quantile(R2_state_null_Jac, c(0.025, 0.975), na.rm = TRUE)
pR2_state_Ait  <- mean(R2_state_null_Ait >= R2_state_obs_Ait, na.rm = TRUE)
pR2_state_Jac  <- mean(R2_state_null_Jac >= R2_state_obs_Jac, na.rm = TRUE)

sesR2_state_Ait <- (R2_state_obs_Ait - mean(R2_state_null_Ait, na.rm = TRUE)) / sd(R2_state_null_Ait, na.rm = TRUE)
sesR2_state_Jac <- (R2_state_obs_Jac - mean(R2_state_null_Jac, na.rm = TRUE)) / sd(R2_state_null_Jac, na.rm = TRUE)

pR2_state_Ait_plot <- ggplot(tibble(R2 = R2_state_null_Ait), aes(x = R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_state_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciR2_state_Ait, linetype = 3) +
  labs(
    title    = paste0(scope, " state R2 null vs obs (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_state_obs_Ait, mean(R2_state_null_Ait, na.rm = TRUE), ciR2_state_Ait[1], ciR2_state_Ait[2], pR2_state_Ait, sesR2_state_Ait
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

print(pR2_state_Ait_plot)
ggsave(file.path(out_state, paste0("state_R2_null_vs_obs_Aitchison_", scope, ".png")), pR2_state_Ait_plot, width = 6, height = 4, dpi = 300)

pR2_state_Jac_plot <- ggplot(tibble(R2 = R2_state_null_Jac), aes(x = R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_state_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciR2_state_Jac, linetype = 3) +
  labs(
    title    = paste0(scope, " state R2 null vs obs (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_state_obs_Jac, mean(R2_state_null_Jac, na.rm = TRUE), ciR2_state_Jac[1], ciR2_state_Jac[2], pR2_state_Jac, sesR2_state_Jac
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

print(pR2_state_Jac_plot)
ggsave(file.path(out_state, paste0("state_R2_null_vs_obs_Jaccard_", scope, ".png")), pR2_state_Jac_plot, width = 6, height = 4, dpi = 300)

AUC_state_obs_Ait <- as.numeric(AUC_state_obs_A[["AUC"]])
AUC_state_obs_Jac <- as.numeric(AUC_state_obs_J[["AUC"]])

ciAUC_state_Ait <- quantile(AUC_state_null_A, c(0.025, 0.975), na.rm = TRUE)
ciAUC_state_Jac <- quantile(AUC_state_null_J, c(0.025, 0.975), na.rm = TRUE)
pAUC_state_Ait  <- mean(AUC_state_null_A >= AUC_state_obs_Ait, na.rm = TRUE)
pAUC_state_Jac  <- mean(AUC_state_null_J >= AUC_state_obs_Jac, na.rm = TRUE)

sesAUC_state_Ait <- (AUC_state_obs_Ait - mean(AUC_state_null_A, na.rm = TRUE)) / sd(AUC_state_null_A, na.rm = TRUE)
sesAUC_state_Jac <- (AUC_state_obs_Jac - mean(AUC_state_null_J, na.rm = TRUE)) / sd(AUC_state_null_J, na.rm = TRUE)

pAUC_state_Ait_plot <- ggplot(tibble(AUC = AUC_state_null_A), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_state_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_state_Ait, linetype = 3) +
  labs(
    title    = paste0(scope, " state AUC null vs observed (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_state_obs_Ait, mean(AUC_state_null_A, na.rm = TRUE), pAUC_state_Ait, sesAUC_state_Ait
    ),
    x = "AUC",
    y = "Density"
  ) +
  theme_bw()

print(pAUC_state_Ait_plot)
ggsave(file.path(out_state, paste0("state_AUC_null_vs_obs_Aitchison_", scope, ".png")), pAUC_state_Ait_plot, width = 6, height = 4, dpi = 300)

pAUC_state_Jac_plot <- ggplot(tibble(AUC = AUC_state_null_J), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_state_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_state_Jac, linetype = 3) +
  labs(
    title    = paste0(scope, " state AUC null vs observed (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_state_obs_Jac, mean(AUC_state_null_J, na.rm = TRUE), pAUC_state_Jac, sesAUC_state_Jac
    ),
    x = "AUC",
    y = "Density"
  ) +
  theme_bw()

print(pAUC_state_Jac_plot)
ggsave(file.path(out_state, paste0("state_AUC_null_vs_obs_Jaccard_", scope, ".png")), pAUC_state_Jac_plot, width = 6, height = 4, dpi = 300)

# ---- ECDF state -------------------------------------------------------------
df_state_Ait <- pairwise_from_D(D_Ait_global, state, upper_Ait)
p_state_ECDF_Ait <- ggplot2::ggplot(df_state_Ait, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " state ECDF of pairwise distances (Aitchison)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_state_obs_A[["AUC"]], AUC_state_obs_A[["r_rb"]],
      sum(df_state_Ait$type == "within"), sum(df_state_Ait$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_state_ECDF_Ait)
ggplot2::ggsave(file.path(out_state, paste0("state_ECDF_Aitchison_", scope, ".png")), p_state_ECDF_Ait, width = 6.2, height = 4.2, dpi = 300)

df_state_Jac <- pairwise_from_D(D_Jac_global, state, upper_Jac)
p_state_ECDF_Jac <- ggplot2::ggplot(df_state_Jac, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " state ECDF of pairwise distances (Jaccard)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_state_obs_J[["AUC"]], AUC_state_obs_J[["r_rb"]],
      sum(df_state_Jac$type == "within"), sum(df_state_Jac$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_state_ECDF_Jac)
ggplot2::ggsave(file.path(out_state, paste0("state_ECDF_Jaccard_", scope, ".png")), p_state_ECDF_Jac, width = 6.2, height = 4.2, dpi = 300)

# =============================================================================
# 4) PAIRED subset with natman blocking
# =============================================================================
tree_counts <- meta_fwd |>
  count(natman, fwd_state) |>
  tidyr::pivot_wider(
    names_from = fwd_state,
    values_from = n,
    values_fill = 0
  )

trees_keep <- tree_counts |>
  filter(ATTACHED > 0, FALLEN > 0) |>
  pull(natman)

cat("\nTrees with both aFWD and fFWD:\n")
print(trees_keep)

meta_pair <- meta_fwd |>
  filter(natman %in% trees_keep) |>
  droplevels()

stopifnot(all(meta_pair$sample %in% rownames(otu_fwd)))

otu_pair <- otu_fwd[meta_pair$sample, , drop = FALSE]
otu_pair <- otu_pair[, colSums(otu_pair, na.rm = TRUE) > 0, drop = FALSE]
stopifnot(identical(rownames(otu_pair), meta_pair$sample))

cat("\nPAIRED data\n")
print(meta_pair |> count(natman, fwd_state))
cat("Samples:", nrow(meta_pair), " OTUs:", ncol(otu_pair), "\n")

D_Ait_pair <- d_ait_or_hel(otu_pair)
D_Ait_pair <- fix_dist_labels(D_Ait_pair, samps = rownames(otu_pair))
D_Jac_pair <- vegdist((otu_pair > 0) * 1, method = "jaccard", binary = TRUE)
D_Jac_pair <- fix_dist_labels(D_Jac_pair, samps = rownames(otu_pair))

perm_blk <- how(nperm = cfg$nperm)
setBlocks(perm_blk) <- meta_pair$natman

ado_Ait_pair <- adonis2(
  D_Ait_pair ~ log_reads + umi + fwd_state,
  data = meta_pair,
  by = "margin",
  permutations = perm_blk
)
ado_Jac_pair <- adonis2(
  D_Jac_pair ~ log_reads + umi + fwd_state,
  data = meta_pair,
  by = "margin",
  permutations = perm_blk
)

cat("\n[PAIRED] blocked PERMANOVA Aitchison\n")
print(ado_Ait_pair)
cat("\n[PAIRED] blocked PERMANOVA Jaccard\n")
print(ado_Jac_pair)

# ---- paired LOO -------------------------------------------------------------
cat("\n[PAIRED LOO PERMANOVA] leave one tree out for fwd_state effect\n")

N_PERM_LOO <- 499L
trees_pair <- levels(meta_pair$natman)

run_blocked_perma_pair <- function(D_full, meta_full, drop_natman, nperm) {
  keep <- meta_full$natman != drop_natman & !is.na(meta_full$natman)
  if (sum(keep) < 6L) return(NULL)
  
  M  <- as.matrix(D_full)
  Mk <- M[keep, keep, drop = FALSE]
  Dk <- stats::as.dist(Mk)
  
  meta_k <- meta_full[keep, , drop = FALSE]
  
  perm_k <- permute::how(nperm = nperm)
  permute::setBlocks(perm_k) <- meta_k$natman
  
  vegan::adonis2(
    Dk ~ log_reads + umi + fwd_state,
    data = meta_k,
    permutations = perm_k,
    by = "margin"
  )
}

loo_pair_list <- purrr::map_dfr(trees_pair, function(tr) {
  cat("  LOO tree", tr, "\n")
  
  fitA <- run_blocked_perma_pair(D_Ait_pair, meta_pair, tr, N_PERM_LOO)
  if (is.null(fitA)) {
    tabA <- tibble()
  } else {
    dfA <- as.data.frame(fitA)
    dfA$term <- rownames(dfA)
    tabA <- dfA %>%
      as_tibble(rownames = "row") %>%
      filter(term %in% c("log_reads", "umi", "fwd_state")) %>%
      transmute(
        natman_left_out = tr,
        metric          = "robust.aitchison",
        term,
        R2              = R2,
        p               = `Pr(>F)`
      )
  }
  
  fitJ <- run_blocked_perma_pair(D_Jac_pair, meta_pair, tr, N_PERM_LOO)
  if (is.null(fitJ)) {
    tabJ <- tibble()
  } else {
    dfJ <- as.data.frame(fitJ)
    dfJ$term <- rownames(dfJ)
    tabJ <- dfJ %>%
      as_tibble(rownames = "row") %>%
      filter(term %in% c("log_reads", "umi", "fwd_state")) %>%
      transmute(
        natman_left_out = tr,
        metric          = "jaccard",
        term,
        R2              = R2,
        p               = `Pr(>F)`
      )
  }
  
  bind_rows(tabA, tabJ)
})

if (nrow(loo_pair_list) > 0L) {
  loo_pair_summary <- loo_pair_list %>%
    group_by(metric, term) %>%
    summarise(
      median_R2 = median(R2, na.rm = TRUE),
      R2_min    = min(R2, na.rm = TRUE),
      R2_max    = max(R2, na.rm = TRUE),
      p_median  = median(p, na.rm = TRUE),
      p_p05     = quantile(p, 0.05, na.rm = TRUE),
      p_p95     = quantile(p, 0.95, na.rm = TRUE),
      n_trees   = n_distinct(natman_left_out),
      .groups   = "drop"
    )
  
  cat("\n[PAIRED LOO PERMANOVA] summary\n")
  print(loo_pair_summary)
  
  write_csv2(loo_pair_list,    file.path(outdir, "PAIRED_LOO_permanova_raw.csv"))
  write_csv2(loo_pair_summary, file.path(outdir, "PAIRED_LOO_permanova_summary.csv"))
} else {
  cat("\n[PAIRED LOO PERMANOVA] skipped because no valid subsets remained\n")
}

# ---- PERMDISP paired --------------------------------------------------------
cat("\n[Dispersion PAIRED] betadisper by fwd_state and natman\n")

disp_pair_state_Ait <- anova(betadisper(D_Ait_pair, meta_pair$fwd_state))
disp_pair_state_Jac <- anova(betadisper(D_Jac_pair, meta_pair$fwd_state))
disp_pair_nat_Ait   <- anova(betadisper(D_Ait_pair, meta_pair$natman))
disp_pair_nat_Jac   <- anova(betadisper(D_Jac_pair, meta_pair$natman))

print(list(
  Aitchison_state_paired  = disp_pair_state_Ait,
  Jaccard_state_paired    = disp_pair_state_Jac,
  Aitchison_natman_paired = disp_pair_nat_Ait,
  Jaccard_natman_paired   = disp_pair_nat_Jac
))

sink(file.path(outdir, paste0("permdisp_PAIRED_", scope, ".txt")))
cat("Aitchison natman (paired)\n");    print(disp_pair_nat_Ait)
cat("Jaccard natman (paired)\n");      print(disp_pair_nat_Jac)
cat("Aitchison fwd_state (paired)\n"); print(disp_pair_state_Ait)
cat("Jaccard fwd_state (paired)\n");   print(disp_pair_state_Jac)
sink()

R2_Ait_state_pair <- as.data.frame(ado_Ait_pair)["fwd_state", "R2"]
R2_Jac_state_pair <- as.data.frame(ado_Jac_pair)["fwd_state", "R2"]

cat("\n[PAIRED] computing state R2 nulls by label shuffles within natman\n")
R2_null_Ait_state <- state_R2_null_blocked(D_Ait_pair, meta_pair, add_log_reads = TRUE)
R2_null_Jac_state <- state_R2_null_blocked(D_Jac_pair, meta_pair, add_log_reads = TRUE)

summ_R2_Ait <- summ_R2_null(R2_Ait_state_pair, R2_null_Ait_state)
summ_R2_Jac <- summ_R2_null(R2_Jac_state_pair, R2_null_Jac_state)

cat("\n[PAIRED] R2 null summary Aitchison\n")
print(summ_R2_Ait)
cat("\n[PAIRED] R2 null summary Jaccard\n")
print(summ_R2_Jac)

# ---- paired state AUC nulls -------------------------------------------------
B_AUC_state_pair <- cfg$B_AUC
AUC_state_null_Ait_pair <- numeric(B_AUC_state_pair)
AUC_state_null_Jac_pair <- numeric(B_AUC_state_pair)

cat("\n[PAIRED] computing state AUC nulls by label shuffles within natman\n")

for (b in seq_len(B_AUC_state_pair)) {
  meta_null <- meta_pair %>%
    mutate(fwd_state = shuffle_state_within_tree(fwd_state, natman))
  
  pairs_A_null <- mk_pair_tbl(D_Ait_pair, meta_null) %>% pairs_by("state_within_tree")
  pairs_J_null <- mk_pair_tbl(D_Jac_pair, meta_null) %>% pairs_by("state_within_tree")
  
  AUC_state_null_Ait_pair[b] <- auc_once(pairs_A_null)$AUC
  AUC_state_null_Jac_pair[b] <- auc_once(pairs_J_null)$AUC
  
  if (b %% 500L == 0L) {
    message("[PAIRED state AUC nulls] finished ", b, " of ", B_AUC_state_pair)
  }
}

pairs_pair_A_state <- mk_pair_tbl(D_Ait_pair, meta_pair) |> pairs_by("state_within_tree")
pairs_pair_J_state <- mk_pair_tbl(D_Jac_pair, meta_pair) |> pairs_by("state_within_tree")

eff_pair_A_state <- median_diff_boot(pairs_pair_A_state)
eff_pair_J_state <- median_diff_boot(pairs_pair_J_state)
auc_pair_A_state <- auc_once(pairs_pair_A_state)
auc_pair_J_state <- auc_once(pairs_pair_J_state)

pairs_pair_A_nat <- mk_pair_tbl(D_Ait_pair, meta_pair) |> pairs_by("natman")
pairs_pair_J_nat <- mk_pair_tbl(D_Jac_pair, meta_pair) |> pairs_by("natman")

eff_pair_A_nat <- median_diff_boot(pairs_pair_A_nat)
eff_pair_J_nat <- median_diff_boot(pairs_pair_J_nat)
auc_pair_A_nat <- auc_once(pairs_pair_A_nat)
auc_pair_J_nat <- auc_once(pairs_pair_J_nat)

AUC_obs_state_Ait <- as.numeric(auc_pair_A_state$AUC)
AUC_obs_state_Jac <- as.numeric(auc_pair_J_state$AUC)

AUC_state_null_mean_Ait_pair <- mean(AUC_state_null_Ait_pair, na.rm = TRUE)
AUC_state_null_mean_Jac_pair <- mean(AUC_state_null_Jac_pair, na.rm = TRUE)

AUC_state_null_ci_Ait <- quantile(AUC_state_null_Ait_pair, c(0.025, 0.975), na.rm = TRUE)
AUC_state_null_ci_Jac <- quantile(AUC_state_null_Jac_pair, c(0.025, 0.975), na.rm = TRUE)

AUC_state_p_emp_Ait <- mean(AUC_state_null_Ait_pair >= AUC_obs_state_Ait, na.rm = TRUE)
AUC_state_p_emp_Jac <- mean(AUC_state_null_Jac_pair >= AUC_obs_state_Jac, na.rm = TRUE)

AUC_state_SES_Ait <- (AUC_obs_state_Ait - AUC_state_null_mean_Ait_pair) / sd(AUC_state_null_Ait_pair, na.rm = TRUE)
AUC_state_SES_Jac <- (AUC_obs_state_Jac - AUC_state_null_mean_Jac_pair) / sd(AUC_state_null_Jac_pair, na.rm = TRUE)

summary_pair_state <- tibble(
  scope    = "PAIRED_within_tree_state",
  metric   = c("robust.aitchison", "jaccard"),
  n_trees  = nlevels(meta_pair$natman),
  n_samples = nrow(meta_pair),
  R2_obs_state = c(R2_Ait_state_pair, R2_Jac_state_pair),
  R2_null_state_mean  = c(summ_R2_Ait$R2_null_mean, summ_R2_Jac$R2_null_mean),
  R2_null_state_CI_lo = c(summ_R2_Ait$R2_null_CI_lo, summ_R2_Jac$R2_null_CI_lo),
  R2_null_state_CI_hi = c(summ_R2_Ait$R2_null_CI_hi, summ_R2_Jac$R2_null_CI_hi),
  R2_state_p_emp      = c(summ_R2_Ait$p_emp, summ_R2_Jac$p_emp),
  R2_state_SES        = c(summ_R2_Ait$SES, summ_R2_Jac$SES),
  AUC_state           = c(AUC_obs_state_Ait, AUC_obs_state_Jac),
  r_rb_state          = c(auc_pair_A_state$r_rb, auc_pair_J_state$r_rb),
  AUC_state_null_mean = c(AUC_state_null_mean_Ait_pair, AUC_state_null_mean_Jac_pair),
  AUC_state_null_CI_lo = c(AUC_state_null_ci_Ait[1], AUC_state_null_ci_Jac[1]),
  AUC_state_null_CI_hi = c(AUC_state_null_ci_Ait[2], AUC_state_null_ci_Jac[2]),
  AUC_state_p_emp      = c(AUC_state_p_emp_Ait, AUC_state_p_emp_Jac),
  AUC_state_SES        = c(AUC_state_SES_Ait, AUC_state_SES_Jac),
  diff_state_median = c(eff_pair_A_state$diff_median, eff_pair_J_state$diff_median),
  diff_state_CI_lo  = c(eff_pair_A_state$diff_CI_lo, eff_pair_J_state$diff_CI_lo),
  diff_state_CI_hi  = c(eff_pair_A_state$diff_CI_hi, eff_pair_J_state$diff_CI_hi),
  n_state_within    = c(eff_pair_A_state$n_within, eff_pair_J_state$n_within),
  n_state_between   = c(eff_pair_A_state$n_between, eff_pair_J_state$n_between)
)

summary_pair_nat <- tibble(
  scope    = "PAIRED_natman",
  metric   = c("robust.aitchison", "jaccard"),
  n_trees  = nlevels(meta_pair$natman),
  n_samples = nrow(meta_pair),
  AUC_nat   = c(auc_pair_A_nat$AUC, auc_pair_J_nat$AUC),
  r_rb_nat  = c(auc_pair_A_nat$r_rb, auc_pair_J_nat$r_rb),
  diff_nat_median = c(eff_pair_A_nat$diff_median, eff_pair_J_nat$diff_median),
  diff_nat_CI_lo  = c(eff_pair_A_nat$diff_CI_lo, eff_pair_J_nat$diff_CI_lo),
  diff_nat_CI_hi  = c(eff_pair_A_nat$diff_CI_hi, eff_pair_J_nat$diff_CI_hi),
  n_nat_within    = c(eff_pair_A_nat$n_within, eff_pair_J_nat$n_within),
  n_nat_between   = c(eff_pair_A_nat$n_between, eff_pair_J_nat$n_between)
)

cat("\n[PAIRED] summary state within trees (with AUC nulls)\n")
print(summary_pair_state)
cat("\n[PAIRED] summary natman within vs between\n")
print(summary_pair_nat)

write_csv2(summary_pair_state, file.path(outdir, "PAIRED_state_FWD_summary.csv"))
write_csv2(summary_pair_nat,   file.path(outdir, "PAIRED_natman_FWD_summary.csv"))

# ---- paired state null plots -----------------------------------------------
out_state_pair <- file.path(outdir, "PAIRED_state")
ensure_dir(out_state_pair)

ciR2_state_Ait_pair <- stats::quantile(R2_null_Ait_state, c(0.025, 0.975), na.rm = TRUE)
ciR2_state_Jac_pair <- stats::quantile(R2_null_Jac_state, c(0.025, 0.975), na.rm = TRUE)

pR2_state_Ait_pair <- mean(R2_null_Ait_state >= R2_Ait_state_pair, na.rm = TRUE)
pR2_state_Jac_pair <- mean(R2_null_Jac_state >= R2_Jac_state_pair, na.rm = TRUE)

sesR2_state_Ait_pair <- (R2_Ait_state_pair - mean(R2_null_Ait_state, na.rm = TRUE)) / stats::sd(R2_null_Ait_state, na.rm = TRUE)
sesR2_state_Jac_pair <- (R2_Jac_state_pair - mean(R2_null_Jac_state, na.rm = TRUE)) / stats::sd(R2_null_Jac_state, na.rm = TRUE)

pR2_state_Ait_pair_plot <- ggplot2::ggplot(tibble::tibble(R2 = R2_null_Ait_state), ggplot2::aes(x = R2)) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = R2_Ait_state_pair, linetype = 2, colour = "blue") +
  ggplot2::geom_vline(xintercept = ciR2_state_Ait_pair, linetype = 3) +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state R2 null vs obs (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_Ait_state_pair, mean(R2_null_Ait_state, na.rm = TRUE),
      ciR2_state_Ait_pair[1], ciR2_state_Ait_pair[2],
      pR2_state_Ait_pair, sesR2_state_Ait_pair
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  ggplot2::theme_bw()

print(pR2_state_Ait_pair_plot)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_R2_null_vs_obs_Aitchison_", scope, ".png")), pR2_state_Ait_pair_plot, width = 6, height = 4, dpi = 300)

pR2_state_Jac_pair_plot <- ggplot2::ggplot(tibble::tibble(R2 = R2_null_Jac_state), ggplot2::aes(x = R2)) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = R2_Jac_state_pair, linetype = 2, colour = "blue") +
  ggplot2::geom_vline(xintercept = ciR2_state_Jac_pair, linetype = 3) +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state R2 null vs obs (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
      R2_Jac_state_pair, mean(R2_null_Jac_state, na.rm = TRUE),
      ciR2_state_Jac_pair[1], ciR2_state_Jac_pair[2],
      pR2_state_Jac_pair, sesR2_state_Jac_pair
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  ggplot2::theme_bw()

print(pR2_state_Jac_pair_plot)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_R2_null_vs_obs_Jaccard_", scope, ".png")), pR2_state_Jac_pair_plot, width = 6, height = 4, dpi = 300)

ciAUC_state_Ait_pair <- stats::quantile(AUC_state_null_Ait_pair, c(0.025, 0.975), na.rm = TRUE)
ciAUC_state_Jac_pair <- stats::quantile(AUC_state_null_Jac_pair, c(0.025, 0.975), na.rm = TRUE)

pAUC_state_Ait_pair <- mean(AUC_state_null_Ait_pair >= AUC_obs_state_Ait, na.rm = TRUE)
pAUC_state_Jac_pair <- mean(AUC_state_null_Jac_pair >= AUC_obs_state_Jac, na.rm = TRUE)

sesAUC_state_Ait_pair <- (AUC_obs_state_Ait - mean(AUC_state_null_Ait_pair, na.rm = TRUE)) / stats::sd(AUC_state_null_Ait_pair, na.rm = TRUE)
sesAUC_state_Jac_pair <- (AUC_obs_state_Jac - mean(AUC_state_null_Jac_pair, na.rm = TRUE)) / stats::sd(AUC_state_null_Jac_pair, na.rm = TRUE)

pAUC_state_Ait_pair_plot <- ggplot2::ggplot(tibble::tibble(AUC = AUC_state_null_Ait_pair), ggplot2::aes(x = AUC)) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = AUC_obs_state_Ait, linetype = 2, colour = "blue") +
  ggplot2::geom_vline(xintercept = ciAUC_state_Ait_pair, linetype = 3) +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state AUC null vs obs (Aitchison)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_obs_state_Ait, mean(AUC_state_null_Ait_pair, na.rm = TRUE),
      pAUC_state_Ait_pair, sesAUC_state_Ait_pair
    ),
    x = "AUC",
    y = "Density"
  ) +
  ggplot2::theme_bw()

print(pAUC_state_Ait_pair_plot)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_AUC_null_vs_obs_Aitchison_", scope, ".png")), pAUC_state_Ait_pair_plot, width = 6, height = 4, dpi = 300)

pAUC_state_Jac_pair_plot <- ggplot2::ggplot(tibble::tibble(AUC = AUC_state_null_Jac_pair), ggplot2::aes(x = AUC)) +
  ggplot2::geom_density(fill = "grey80") +
  ggplot2::geom_vline(xintercept = AUC_obs_state_Jac, linetype = 2, colour = "blue") +
  ggplot2::geom_vline(xintercept = ciAUC_state_Jac_pair, linetype = 3) +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state AUC null vs obs (Jaccard)"),
    subtitle = sprintf(
      "obs=%.3f  null mean=%.3f  p(right-tailed)=%.3f  SES=%.2f",
      AUC_obs_state_Jac, mean(AUC_state_null_Jac_pair, na.rm = TRUE),
      pAUC_state_Jac_pair, sesAUC_state_Jac_pair
    ),
    x = "AUC",
    y = "Density"
  ) +
  ggplot2::theme_bw()

print(pAUC_state_Jac_pair_plot)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_AUC_null_vs_obs_Jaccard_", scope, ".png")), pAUC_state_Jac_pair_plot, width = 6, height = 4, dpi = 300)

# ---- paired state ECDF ------------------------------------------------------
p_state_pair_ECDF_Ait <- ggplot2::ggplot(pairs_pair_A_state, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state ECDF of pairwise distances (Aitchison)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_obs_state_Ait, auc_pair_A_state[["r_rb"]],
      sum(pairs_pair_A_state$type == "within"),
      sum(pairs_pair_A_state$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_state_pair_ECDF_Ait)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_ECDF_Aitchison_", scope, ".png")), p_state_pair_ECDF_Ait, width = 6.2, height = 4.2, dpi = 300)

p_state_pair_ECDF_Jac <- ggplot2::ggplot(pairs_pair_J_state, ggplot2::aes(x = dist, colour = type)) +
  ggplot2::stat_ecdf(geom = "step", linewidth = 1) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title    = paste0(scope, " PAIRED state ECDF of pairwise distances (Jaccard)"),
    subtitle = sprintf(
      "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
      AUC_obs_state_Jac, auc_pair_J_state[["r_rb"]],
      sum(pairs_pair_J_state$type == "within"),
      sum(pairs_pair_J_state$type == "between")
    ),
    x = "Pairwise distance",
    y = "Cumulative probability",
    colour = "Pair type"
  )

print(p_state_pair_ECDF_Jac)
ggplot2::ggsave(file.path(out_state_pair, paste0("PAIRED_state_ECDF_Jaccard_", scope, ".png")), p_state_pair_ECDF_Jac, width = 6.2, height = 4.2, dpi = 300)

# ---- paired ANOSIM state ----------------------------------------------------
N_PERM_ANOSIM_PAIR <- 9999L

cat("\n[PAIRED ANOSIM] state (fwd_state)\n")

anos_state_Ait_pair <- vegan::anosim(D_Ait_pair, grouping = meta_pair$fwd_state, permutations = N_PERM_ANOSIM_PAIR)
anos_state_Jac_pair <- vegan::anosim(D_Jac_pair, grouping = meta_pair$fwd_state, permutations = N_PERM_ANOSIM_PAIR)

print(anos_state_Ait_pair)
print(anos_state_Jac_pair)

# =============================================================================
# 5) Betapart paired subset
# =============================================================================
PA_pair <- (otu_pair > 0) * 1
beta_pair <- betapart::beta.pair(PA_pair, index.family = "sorensen")

summ_betapart <- function(D, meta, label) {
  pairs <- mk_pair_tbl(D, meta)
  by_nat   <- pairs_by(pairs, "natman")
  by_state <- pairs_by(pairs, "state_within_tree")
  
  tibble(
    component = label,
    within_nat_median   = median(by_nat$dist[by_nat$type == "within"]),
    between_nat_median  = median(by_nat$dist[by_nat$type == "between"]),
    within_state_median = median(by_state$dist[by_state$type == "within"]),
    cross_state_median  = median(by_state$dist[by_state$type == "cross"])
  )
}

beta_summary <- bind_rows(
  summ_betapart(beta_pair$beta.sim, meta_pair, "beta.SIM"),
  summ_betapart(beta_pair$beta.sne, meta_pair, "beta.SNE"),
  summ_betapart(beta_pair$beta.sor, meta_pair, "beta.SOR")
)

beta_summary <- beta_summary %>%
  mutate(
    SIM_over_SNE_nat = if_else(
      component == "beta.SIM",
      within_nat_median / beta_summary$within_nat_median[component == "beta.SNE"],
      NA_real_
    ),
    SIM_over_SNE_state = if_else(
      component == "beta.SIM",
      cross_state_median / beta_summary$cross_state_median[component == "beta.SNE"],
      NA_real_
    )
  )

cat("\n[PAIRED] betapart summary (turnover vs nestedness)\n")
print(beta_summary)

write_csv2(beta_summary, file.path(outdir, "PAIRED_betapart_turnover_nestedness_summary.csv"))

# =============================================================================
# 6) Final console summary
# =============================================================================
cat("\n✓ aFWD vs fFWD compact analysis finished. Outputs in:\n")
print(normalizePath(outdir))

cat("\n[GLOBAL] natman evidence summary\n")
print(natman_evidence)

cat("\n[GLOBAL] state evidence summary\n")
print(state_evidence)

cat("\n[PAIRED] summary state within trees\n")
print(summary_pair_state)

cat("\n[PAIRED] summary natman within vs between\n")
print(summary_pair_nat)
