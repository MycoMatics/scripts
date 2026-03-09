# =============================================================================
# Script: Script7_scope2D_fFWD.R
# Project: Deadwood fungal community analyses
# Scope: fFWD / FALLEN
# Purpose: Analyse fallen fine woody debris communities using PERMANOVA,
#          blocked PERMANOVA, LOO stability, PERMDISP, natman null models
#          for R2 and AUC, ECDF summaries, ANOSIM, betapart decomposition,
#          and alpha richness models
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - as_pa
# - fix_dist_labels
# - get_R2_table
# - loo_permanova
# - mk_within_between_tbl
# - auroc_from_ranks
# - plot_ecdf
# - boot_diff_median
# - overdisp_phi
# - check_glmm_lme4
# - plot_resid_glmm
# - ensure_dir
# - assert_objects
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(betapart)
  library(indicspecies)
  library(ggrepel)
  library(ape)
  library(MASS)
  library(glmmTMB)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c(
  "as_pa", "fix_dist_labels", "get_R2_table", "loo_permanova",
  "mk_within_between_tbl", "auroc_from_ranks", "plot_ecdf",
  "boot_diff_median", "overdisp_phi", "check_glmm_lme4",
  "plot_resid_glmm", "ensure_dir"
))

set.seed(42)
setwd(".")

# -------------------------------
# 0.1) Config
# -------------------------------
cfg <- list(
  outdir      = file.path("plots", "FWD", "FALLEN"),
  tabledir    = file.path("tables", "FWD", "FALLEN"),
  seed        = 42,
  B           = 199,
  B_R2        = 999,
  B_loo       = 499,
  dist_abund  = "robust.aitchison",
  dist_pa     = "jaccard",
  min_group_n = 3
)

ensure_dir(cfg$outdir)
ensure_dir(cfg$tabledir)

# -------------------------------
# 1) Subset to FALLEN FWD
# -------------------------------
stopifnot(exists("META1"), exists("otu_matrix_filt"))

meta <- META1 %>%
  filter(dw_type2 =="fFWD") %>%
  filter(!is.na(dw_type2))

otu <- otu_matrix_filt[match(meta$sample, rownames(otu_matrix_filt)), , drop = FALSE]

meta <- meta %>%
  mutate(
    libsize   = rowSums(otu),
    log_reads = log1p(libsize)
  )

otu <- otu[, colSums(otu) > 0, drop = FALSE]
stopifnot(nrow(otu) == nrow(meta))

keep_nat <- meta %>%
  count(natman) %>%
  filter(n >= cfg$min_group_n) %>%
  pull(natman)

meta <- meta %>% filter(natman %in% keep_nat)
otu  <- otu[match(meta$sample, rownames(otu)), , drop = FALSE]

meta <- meta %>%
  mutate(
    natman      = as.factor(natman),
    decay_stage = as.factor(decay_stage),
    log_reads   = as.numeric(log_reads),
    libsize     = rowSums(otu)
  )

# -------------------------------
# 2) Alpha diversity q0 visualisation
# -------------------------------
meta$q0 <- rowSums(otu > 0)

p_alpha <- ggplot(meta, aes(natman, q0, colour = natman)) +
  geom_boxplot() +
  labs(
    title = "FALLEN FWD — Alpha diversity (q0)",
    x = "Tree legacy (natman)",
    y = "Observed richness (q0)"
  ) +
  theme(legend.position = "none")

print(p_alpha)
ggsave(file.path(cfg$outdir, "alpha_q0_by_natman.png"), p_alpha, width = 8, height = 4.5, dpi = 300)

# -------------------------------
# 3) Distances and PERMANOVA
# -------------------------------
set.seed(cfg$seed)

D_robAit <- vegdist(otu, method = cfg$dist_abund)
D_jac    <- vegdist(as_pa(otu), method = cfg$dist_pa, binary = TRUE)

D_robAit <- fix_dist_labels(D_robAit, samps = rownames(otu))
D_jac    <- fix_dist_labels(D_jac, samps = rownames(otu))

perm_free <- how(nperm = cfg$B)

perma_unblocked_ra <- adonis2(
  D_robAit ~ log_reads + umi + natman + decay_stage + diameter_at_drill_z,
  data = meta,
  permutations = perm_free,
  by = "terms"
)
print(perma_unblocked_ra)

perma_unblocked_jac <- adonis2(
  D_jac ~ log_reads + umi + natman + decay_stage + diameter_at_drill_z,
  data = meta,
  permutations = perm_free,
  by = "terms"
)
print(perma_unblocked_jac)

perm_blocked <- how(nperm = cfg$B)
setBlocks(perm_blocked) <- meta$natman

perma_blocked_ra <- adonis2(
  D_robAit ~ log_reads + umi + decay_stage + diameter_at_drill_z,
  data = meta,
  permutations = perm_blocked,
  by = "margin"
)
print(perma_blocked_ra)

perma_blocked_jac <- adonis2(
  D_jac ~ log_reads + umi + decay_stage + diameter_at_drill_z,
  data = meta,
  permutations = perm_blocked,
  by = "margin"
)
print(perma_blocked_jac)

R2_unblocked_ra  <- get_R2_table(perma_unblocked_ra, D_robAit)
R2_unblocked_jac <- get_R2_table(perma_unblocked_jac, D_jac)
R2_blocked_ra    <- get_R2_table(perma_blocked_ra, D_robAit)
R2_blocked_jac   <- get_R2_table(perma_blocked_jac, D_jac)

readr::write_csv2(R2_unblocked_ra,  file.path(cfg$tabledir, "permanova_unblocked_robustAitchison.csv"))
readr::write_csv2(R2_unblocked_jac, file.path(cfg$tabledir, "permanova_unblocked_jaccard.csv"))
readr::write_csv2(R2_blocked_ra,    file.path(cfg$tabledir, "permanova_blocked_robustAitchison.csv"))
readr::write_csv2(R2_blocked_jac,   file.path(cfg$tabledir, "permanova_blocked_jaccard.csv"))

# -------------------------------
# 4) LOO PERMANOVA
# -------------------------------
loo_permanova(
  X = otu,
  meta = meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + decay_stage + diameter_at_drill_z",
  dist_method = "robust.aitchison",
  keep_terms = c("decay_stage", "diameter_at_drill_z", "log_reads", "umi"),
  blocked_by = "natman",
  nperm_outer = cfg$B_loo,
  by_TYPE = "margin",
  out_csv = file.path(cfg$tabledir, "loo_fallen_decay_robustAitchison.csv")
)

loo_permanova(
  X = otu,
  meta = meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + decay_stage + diameter_at_drill_z",
  dist_method = "jaccard",
  keep_terms = c("decay_stage", "diameter_at_drill_z", "log_reads", "umi"),
  blocked_by = "natman",
  nperm_outer = cfg$B_loo,
  by_TYPE = "margin",
  out_csv = file.path(cfg$tabledir, "loo_fallen_decay_jaccard.csv")
)

# -------------------------------
# 5) PERMDISP
# -------------------------------
message("[PERMDISP] Testing group dispersions …")

bd_ra <- betadisper(D_robAit, group = meta$natman)
bd_j  <- betadisper(D_jac,    group = meta$natman)

pt_ra <- permutest(bd_ra, permutations = cfg$B)
pt_j  <- permutest(bd_j,  permutations = cfg$B)

print(anova(bd_ra)); print(pt_ra)
print(anova(bd_j));  print(pt_j)

bd_rads <- betadisper(D_robAit, group = meta$decay_stage)
bd_jds  <- betadisper(D_jac,    group = meta$decay_stage)

pt_rads <- permutest(bd_rads, permutations = cfg$B)
pt_jds  <- permutest(bd_jds,  permutations = cfg$B)

print(anova(bd_rads)); print(pt_rads)
print(anova(bd_jds));  print(pt_jds)

capture.output({
  cat("\n== PERMDISP natman (Robust Aitchison) ==\n"); print(anova(bd_ra)); print(pt_ra)
  cat("\n== PERMDISP natman (Jaccard) ==\n");        print(anova(bd_j));  print(pt_j)
  cat("\n== PERMDISP decay_stage (Robust Aitchison) ==\n"); print(anova(bd_rads)); print(pt_rads)
  cat("\n== PERMDISP decay_stage (Jaccard) ==\n");          print(anova(bd_jds));  print(pt_jds)
}, file = file.path(cfg$tabledir, "PERMDISP_summary.txt"))

# -------------------------------
# 6) NATMAN NULL — R² null
# -------------------------------
set.seed(cfg$seed)

perm_nat_null_once <- function(groups) sample(groups, length(groups), replace = FALSE)
N_ADONIS_IN <- 199

R2_nat_null <- replicate(cfg$B_R2, {
  meta_tmp <- meta %>% mutate(nat_null = perm_nat_null_once(natman))
  
  ad_a <- adonis2(
    D_robAit ~ log_reads + umi + nat_null + decay_stage + diameter_at_drill_z,
    data = meta_tmp,
    by = "terms",
    permutations = N_ADONIS_IN
  )
  
  ad_j <- adonis2(
    D_jac ~ log_reads + umi + nat_null + decay_stage + diameter_at_drill_z,
    data = meta_tmp,
    by = "terms",
    permutations = N_ADONIS_IN
  )
  
  c(
    robust_aitchison = as.numeric(as.data.frame(ad_a)["nat_null", "R2"]),
    jaccard          = as.numeric(as.data.frame(ad_j)["nat_null", "R2"])
  )
}, simplify = "matrix")

saveRDS(R2_nat_null, file.path(cfg$outdir, "NATMAN_NULL_R2.rds"))

R2_nat_null <- as_tibble(t(R2_nat_null))
write_csv(R2_nat_null, file.path(cfg$tabledir, "NATMAN_NULL_R2.csv"))

R2_obs_ra  <- as.data.frame(perma_unblocked_ra)["natman", "R2"]
R2_obs_jac <- as.data.frame(perma_unblocked_jac)["natman", "R2"]

ci_ra   <- quantile(R2_nat_null$robust_aitchison, c(0.025, 0.975), na.rm = TRUE)
ci_jac  <- quantile(R2_nat_null$jaccard, c(0.025, 0.975), na.rm = TRUE)
pemp_ra <- mean(R2_nat_null$robust_aitchison >= R2_obs_ra, na.rm = TRUE)
pemp_jac <- mean(R2_nat_null$jaccard >= R2_obs_jac, na.rm = TRUE)

ses_ra <- (R2_obs_ra - mean(R2_nat_null$robust_aitchison, na.rm = TRUE)) /
  sd(R2_nat_null$robust_aitchison, na.rm = TRUE)
ses_jac <- (R2_obs_jac - mean(R2_nat_null$jaccard, na.rm = TRUE)) /
  sd(R2_nat_null$jaccard, na.rm = TRUE)

p_null_ra <- ggplot(R2_nat_null, aes(x = robust_aitchison)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_ra, colour = "blue", linetype = 2) +
  geom_vline(xintercept = ci_ra, colour = "gray30", linetype = 3) +
  theme_bw() +
  labs(
    title = "FALLEN FWD — natman R² null vs observed (robust Aitchison)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      R2_obs_ra, mean(R2_nat_null$robust_aitchison, na.rm = TRUE),
      ci_ra[1], ci_ra[2], pemp_ra, ses_ra
    ),
    x = expression(R^2), y = "Density"
  )

print(p_null_ra)
ggsave(file.path(cfg$outdir, "natman_R2_null_vs_obs_FALLEN_robustAitchison.png"), p_null_ra, width = 6, height = 4, dpi = 300)

p_null_jac <- ggplot(R2_nat_null, aes(x = jaccard)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_jac, colour = "blue", linetype = 2) +
  geom_vline(xintercept = ci_jac, colour = "grey30", linetype = 3) +
  theme_bw() +
  labs(
    title = "FALLEN FWD — natman R² null vs observed (Jaccard)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      R2_obs_jac, mean(R2_nat_null$jaccard, na.rm = TRUE),
      ci_jac[1], ci_jac[2], pemp_jac, ses_jac
    ),
    x = expression(R^2), y = "Density"
  )

print(p_null_jac)
ggsave(file.path(cfg$outdir, "natman_R2_null_vs_obs_FALLEN_jaccard.png"), p_null_jac, width = 6, height = 4, dpi = 300)

# -------------------------------
# 7) ECDF / AUC
# -------------------------------
mk_within_between_tbl <- function(DM, group){
  if (inherits(DM, "dist")) DM <- as.matrix(DM)
  stopifnot(is.matrix(DM), nrow(DM) == length(group))
  
  out <- tibble()
  for (i in seq_len(nrow(DM) - 1)){
    for (j in (i+1):nrow(DM)){
      out <- bind_rows(out, tibble(
        i = i, j = j,
        dist = DM[i,j],
        type = ifelse(group[i] == group[j], "within", "between")
      ))
    }
  }
  out
}

pairs_ra <- mk_within_between_tbl(D_robAit, meta$natman)
pairs_j  <- mk_within_between_tbl(D_jac,    meta$natman)

AUC_ra <- auroc_from_ranks(
  pairs_ra$dist[pairs_ra$type == "within"],
  pairs_ra$dist[pairs_ra$type == "between"]
)
AUC_j <- auroc_from_ranks(
  pairs_j$dist[pairs_j$type == "within"],
  pairs_j$dist[pairs_j$type == "between"]
)

rb_ra <- 2 * AUC_ra - 1
rb_j  <- 2 * AUC_j - 1

p_ecdf_ra <- plot_ecdf(pairs_ra, "Robust Aitchison – within vs between (natman)")
print(p_ecdf_ra)
ggsave(file.path(cfg$outdir, "ECDF_robust_aitchison_within_between.png"), p_ecdf_ra, width = 6, height = 4, dpi = 300)

p_ecdf_j <- plot_ecdf(pairs_j, "Jaccard – within vs between (natman)")
print(p_ecdf_j)
ggsave(file.path(cfg$outdir, "ECDF_jaccard_within_between.png"), p_ecdf_j, width = 6, height = 4, dpi = 300)

# -------------------------------
# 8) ANOSIM + AUC nulls
# -------------------------------
an_ra <- anosim(D_robAit, meta$natman, permutations = cfg$B)
print(an_ra)
an_j  <- anosim(D_jac, meta$natman, permutations = cfg$B)
print(an_j)

capture.output({
  cat("\n== ANOSIM (Robust Aitchison) ==\n"); print(an_ra)
  cat("\n== ANOSIM (Jaccard) ==\n");         print(an_j)
}, file = file.path(cfg$tabledir, "ANOSIM_natman.txt"))

R2_null_mean_ra  <- mean(R2_nat_null$robust_aitchison, na.rm = TRUE)
R2_null_mean_jac <- mean(R2_nat_null$jaccard, na.rm = TRUE)
R2_percentile_ra  <- mean(R2_nat_null$robust_aitchison <= R2_obs_ra, na.rm = TRUE)
R2_percentile_jac <- mean(R2_nat_null$jaccard <= R2_obs_jac, na.rm = TRUE)

set.seed(cfg$seed)
B_AUC <- cfg$B

perm_once_groups <- function(groups) sample(groups, length(groups), replace = FALSE)

auc_null_from_perm <- function(D, groups, B) {
  replicate(B, {
    g_null <- perm_once_groups(groups)
    pr <- mk_within_between_tbl(D, g_null)
    auroc_from_ranks(
      pr$dist[pr$type == "within"],
      pr$dist[pr$type == "between"]
    )
  })
}

AUC_null_ra  <- auc_null_from_perm(D_robAit, meta$natman, B_AUC)
AUC_null_jac <- auc_null_from_perm(D_jac, meta$natman, B_AUC)

AUC_null_mean_ra  <- mean(AUC_null_ra, na.rm = TRUE)
AUC_null_mean_jac <- mean(AUC_null_jac, na.rm = TRUE)

AUC_null_ci_ra  <- quantile(AUC_null_ra, c(0.025, 0.975), na.rm = TRUE)
AUC_null_ci_jac <- quantile(AUC_null_jac, c(0.025, 0.975), na.rm = TRUE)

pAUC_ra  <- mean(AUC_null_ra >= AUC_ra, na.rm = TRUE)
pAUC_jac <- mean(AUC_null_jac >= AUC_j, na.rm = TRUE)

AUC_SES_ra  <- (AUC_ra - AUC_null_mean_ra) / sd(AUC_null_ra, na.rm = TRUE)
AUC_SES_jac <- (AUC_j - AUC_null_mean_jac) / sd(AUC_null_jac, na.rm = TRUE)

p_AUC_ra <- ggplot(tibble(AUC = AUC_null_ra), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_ra, colour = "blue", linetype = 2) +
  geom_vline(xintercept = AUC_null_ci_ra, colour = "grey30", linetype = 3) +
  theme_bw() +
  labs(
    title = "FALLEN FWD natman AUC null vs observed (robust Aitchison)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      AUC_ra, AUC_null_mean_ra, AUC_null_ci_ra[1], AUC_null_ci_ra[2], pAUC_ra, AUC_SES_ra
    ),
    x = "AUC (within vs between natman)",
    y = "Density"
  )

print(p_AUC_ra)
ggsave(file.path(cfg$outdir, "natman_AUC_null_vs_obs_FALLEN_robustAitchison.png"), p_AUC_ra, width = 6, height = 4, dpi = 300)

p_AUC_jac <- ggplot(tibble(AUC = AUC_null_jac), aes(x = AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_j, colour = "blue", linetype = 2) +
  geom_vline(xintercept = AUC_null_ci_jac, colour = "grey30", linetype = 3) +
  theme_bw() +
  labs(
    title = "FALLEN FWD natman AUC null vs observed (Jaccard)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      AUC_j, AUC_null_mean_jac, AUC_null_ci_jac[1], AUC_null_ci_jac[2], pAUC_jac, AUC_SES_jac
    ),
    x = "AUC (within vs between natman)",
    y = "Density"
  )

print(p_AUC_jac)
ggsave(file.path(cfg$outdir, "natman_AUC_null_vs_obs_FALLEN_jaccard.png"), p_AUC_jac, width = 6, height = 4, dpi = 300)

n_trees   <- dplyr::n_distinct(meta$natman)
n_samples <- nrow(meta)

fmt_ci <- function(lo, hi) sprintf("%.3f–%.3f", lo, hi)

row_ra <- tibble::tibble(
  Scope = "fFWD",
  Metric = "robust aitchison",
  `R² obs` = round(R2_obs_ra, 3),
  `R² null mean` = round(R2_null_mean_ra, 3),
  `R² null 95% CI` = fmt_ci(ci_ra[1], ci_ra[2]),
  `p (R²)` = ifelse(pemp_ra < 0.001, "<0.001", sprintf("%.3f", pemp_ra)),
  `R² SES` = round(ses_ra, 3),
  `R² perc.` = round(R2_percentile_ra, 3),
  `AUC obs` = round(AUC_ra, 3),
  `AUC null mean` = round(AUC_null_mean_ra, 3),
  `AUC null 95% CI` = fmt_ci(AUC_null_ci_ra[1], AUC_null_ci_ra[2]),
  `p (AUC)` = ifelse(pAUC_ra < 0.001, "<0.001", sprintf("%.3f", pAUC_ra)),
  `AUC SES` = round(AUC_SES_ra, 3),
  `ANOSIM R` = round(as.numeric(an_ra$statistic), 3),
  `ANOSIM p` = ifelse(an_ra$signif < 0.001, "<0.001", sprintf("%.3f", an_ra$signif)),
  `n trees` = n_trees,
  `n samples` = n_samples
)

row_jac <- tibble::tibble(
  Scope = "fFWD",
  Metric = "jaccard",
  `R² obs` = round(R2_obs_jac, 3),
  `R² null mean` = round(R2_null_mean_jac, 3),
  `R² null 95% CI` = fmt_ci(ci_jac[1], ci_jac[2]),
  `p (R²)` = ifelse(pemp_jac < 0.001, "<0.001", sprintf("%.3f", pemp_jac)),
  `R² SES` = round(ses_jac, 3),
  `R² perc.` = round(R2_percentile_jac, 3),
  `AUC obs` = round(AUC_j, 3),
  `AUC null mean` = round(AUC_null_mean_jac, 3),
  `AUC null 95% CI` = fmt_ci(AUC_null_ci_jac[1], AUC_null_ci_jac[2]),
  `p (AUC)` = ifelse(pAUC_jac < 0.001, "<0.001", sprintf("%.3f", pAUC_jac)),
  `AUC SES` = round(AUC_SES_jac, 3),
  `ANOSIM R` = round(as.numeric(an_j$statistic), 3),
  `ANOSIM p` = ifelse(an_j$signif < 0.001, "<0.001", sprintf("%.3f", an_j$signif)),
  `n trees` = n_trees,
  `n samples` = n_samples
)

harmonized_ffwd <- dplyr::bind_rows(row_ra, row_jac)
print(harmonized_ffwd)
readr::write_csv2(harmonized_ffwd, file.path(cfg$tabledir, "natman_evidence_summary_fFWD.csv"))

# -------------------------------
# 9) BETAPART
# -------------------------------
PA <- as_pa(otu)
bp <- betapart::beta.pair(PA, index.family = "sorensen")

DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

pairs_sim <- mk_within_between_tbl(DM_sim, meta$natman)
pairs_sne <- mk_within_between_tbl(DM_sne, meta$natman)
pairs_sor <- mk_within_between_tbl(DM_sor, meta$natman)

sum_sim <- pairs_sim %>%
  group_by(type) %>%
  summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop") %>%
  mutate(metric = "beta.sim")

sum_sne <- pairs_sne %>%
  group_by(type) %>%
  summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop") %>%
  mutate(metric = "beta.sne")

sum_sor <- pairs_sor %>%
  group_by(type) %>%
  summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop") %>%
  mutate(metric = "beta.sor")

betapart_summary <- bind_rows(sum_sim, sum_sne, sum_sor) %>%
  relocate(metric, .before = type)

readr::write_csv2(betapart_summary, file.path(cfg$tabledir, "betapart_summary_basic.csv"))

mean_sim <- pairs_sim %>% group_by(type) %>% summarise(mean = mean(dist), .groups = "drop")
mean_sne <- pairs_sne %>% group_by(type) %>% summarise(mean = mean(dist), .groups = "drop")

ratio_between <- (mean_sim %>% filter(type == "between") %>% pull(mean)) /
  (mean_sne %>% filter(type == "between") %>% pull(mean))
ratio_within <- (mean_sim %>% filter(type == "within") %>% pull(mean)) /
  (mean_sne %>% filter(type == "within") %>% pull(mean))

bd_sim <- boot_diff_median(pairs_sim, nboot = 1999, seed = cfg$seed)
bd_sne <- boot_diff_median(pairs_sne, nboot = 1999, seed = cfg$seed)
bd_sor <- boot_diff_median(pairs_sor, nboot = 1999, seed = cfg$seed)

ci_sim <- quantile(bd_sim, c(0.025, 0.975), na.rm = TRUE)
ci_sne <- quantile(bd_sne, c(0.025, 0.975), na.rm = TRUE)
ci_sor <- quantile(bd_sor, c(0.025, 0.975), na.rm = TRUE)

est_sim <- with(sum_sim, median[type == "between"] - median[type == "within"])
est_sne <- with(sum_sne, median[type == "between"] - median[type == "within"])
est_sor <- with(sum_sor, median[type == "between"] - median[type == "within"])

scope_tag <- "fFWD"

mk_row <- function(metric_label, sum_tbl, est, ci, ratio_b, ratio_w) {
  tibble::tibble(
    Scope = scope_tag,
    metric = metric_label,
    est_diff_median = est,
    CI_low = ci[[1]],
    CI_high = ci[[2]],
    n_pairs_between = sum_tbl$n[sum_tbl$type == "between"],
    n_pairs_within  = sum_tbl$n[sum_tbl$type == "within"],
    mean_between    = if ("mean" %in% names(sum_tbl)) sum_tbl$mean[sum_tbl$type == "between"] else NA_real_,
    mean_within     = if ("mean" %in% names(sum_tbl)) sum_tbl$mean[sum_tbl$type == "within"] else NA_real_,
    median_between  = sum_tbl$median[sum_tbl$type == "between"],
    median_within   = sum_tbl$median[sum_tbl$type == "within"],
    IQR_between     = sum_tbl$IQR[sum_tbl$type == "between"],
    IQR_within      = sum_tbl$IQR[sum_tbl$type == "within"],
    turnover_to_nestedness_between = ratio_b,
    turnover_to_nestedness_within  = ratio_w
  )
}

sum_sim_out <- pairs_sim %>%
  group_by(type) %>%
  summarise(n = n(), mean = mean(dist), median = median(dist), IQR = IQR(dist), .groups = "drop")

sum_sne_out <- pairs_sne %>%
  group_by(type) %>%
  summarise(n = n(), mean = mean(dist), median = median(dist), IQR = IQR(dist), .groups = "drop")

sum_sor_out <- pairs_sor %>%
  group_by(type) %>%
  summarise(n = n(), mean = mean(dist), median = median(dist), IQR = IQR(dist), .groups = "drop")

betapart_enriched <- dplyr::bind_rows(
  mk_row("βSIM (turnover)",   sum_sim_out, est_sim, ci_sim, ratio_between, ratio_within),
  mk_row("βSNE (nestedness)", sum_sne_out, est_sne, ci_sne, ratio_between, ratio_within),
  mk_row("βSOR (total)",      sum_sor_out, est_sor, ci_sor, ratio_between, ratio_within)
)

print(betapart_enriched)
readr::write_csv2(betapart_enriched, file.path(cfg$tabledir, "betapart_within_between_enriched.csv"))

# -------------------------------
# 10) Alpha Diversity Modeling (q0)
# -------------------------------
scope <- "fFWD"
cat("\n[ ALPHA RICHNESS GLMM - fFWD ]\n")

min_abund <- 1
rich_ffwd <- tibble(
  sample   = rownames(otu),
  richness = rowSums(as.matrix(otu) >= min_abund, na.rm = TRUE)
)

df_fFWD <- meta %>%
  filter(dw_type2 == "fFWD") %>%
  left_join(rich_ffwd, by = "sample") %>%
  mutate(
    diameter_z  = diameter_at_drill_z,
    natman      = droplevels(factor(natman)),
    decay_stage = droplevels(factor(decay_stage)),
    umi         = droplevels(factor(umi))
  ) %>%
  filter(!is.na(richness), richness > 0)

rhs_terms_glmm_fFWD <- c("log_reads", "umi", "(1|natman)", "decay_stage", "diameter_z")
form_fFWD <- as.formula(paste("richness ~", paste(rhs_terms_glmm_fFWD, collapse = " + ")))
cat("GLMM form (fFWD):\n")
print(form_fFWD)

df_fFWD$decay_stage <- factor(df_fFWD$decay_stage)

if (nlevels(df_fFWD$decay_stage) >= 2) {
  ctrl_fFWD <- lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
  
  m_pois_fFWD <- suppressWarnings(
    lme4::glmer(
      form_fFWD,
      data = df_fFWD,
      family = poisson(link = "log"),
      control = ctrl_fFWD
    )
  )
  
  phi_p_fFWD  <- overdisp_phi(m_pois_fFWD)
  use_nb_fFWD <- is.finite(phi_p_fFWD$phi) && (phi_p_fFWD$phi > 1.2)
  
  if (use_nb_fFWD) {
    m_nb_fFWD <- suppressWarnings(
      lme4::glmer.nb(form_fFWD, data = df_fFWD, control = ctrl_fFWD)
    )
    best_fFWD <- if (AIC(m_pois_fFWD) + 2 < AIC(m_nb_fFWD)) m_pois_fFWD else m_nb_fFWD
    best_type_fFWD <- if (identical(best_fFWD, m_pois_fFWD)) "Poisson (glmer)" else "NB (glmer.nb)"
  } else {
    best_fFWD <- m_pois_fFWD
    best_type_fFWD <- "Poisson (glmer)"
  }
  
  cat(sprintf(
    "Selected fFWD alpha model: %s  Overdispersion Poisson phi = %.2f  AIC(best) = %.1f\n",
    best_type_fFWD, phi_p_fFWD$phi, AIC(best_fFWD)
  ))
  print(summary(best_fFWD))
  
  check_glmm_lme4(best_fFWD, name = "fFWD alpha")
  plot_resid_glmm(best_fFWD, file.path(cfg$outdir, "alpha_richness_fFWD_residuals.png"))
  
  an_tab_fFWD <- car::Anova(best_fFWD, type = 2)
  r2_fFWD     <- performance::r2(best_fFWD)
  
  irr_fFWD <- broom.mixed::tidy(
    best_fFWD,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      term = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "log_reads" ~ "log(reads+1)",
        term == "diameter_z" ~ "Diameter (z)",
        grepl("^decay_stage", term) ~ paste0(
          "Decay stage ",
          sub("^decay_stage", "", term),
          " vs ", levels(df_fFWD$decay_stage)[1]
        ),
        TRUE ~ term
      )
    ) %>%
    dplyr::rename(
      IRR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p = p.value
    ) %>%
    dplyr::arrange(term)
  
  emm_ds_fFWD <- emmeans::emmeans(best_fFWD, ~ decay_stage, type = "response")
  ds_df_fFWD  <- summary(emm_ds_fFWD) %>% as.data.frame()
  
  p_emm_ds_fFWD <- ggplot2::ggplot(
    ds_df_fFWD,
    ggplot2::aes(x = decay_stage, y = response)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL),
      width = 0.2
    ) +
    ggplot2::labs(
      x = "Decay stage",
      y = "Estimated richness",
      title = paste0("Richness ~ decay (", scope, ")")
    ) +
    ggplot2::theme_classic(base_size = 11)
  
  ggplot2::ggsave(
    file.path(cfg$outdir, paste0("alpha_richness_emm_decay_", scope, ".png")),
    p_emm_ds_fFWD, width = 6.2, height = 3.2, dpi = 300
  )
  print(p_emm_ds_fFWD)
  
  keep_eff_fFWD <- c("decay_stage", "diameter_z", "log_reads")
  tab_anova_fFWD <- tibble::tibble(
    Effect = rownames(an_tab_fFWD),
    Chisq  = an_tab_fFWD$`Chisq`,
    df     = an_tab_fFWD$Df,
    p      = an_tab_fFWD$`Pr(>Chisq)`
  ) %>%
    dplyr::filter(Effect %in% keep_eff_fFWD) %>%
    dplyr::mutate(p_BH = p.adjust(p, "BH")) %>%
    dplyr::arrange(match(Effect, keep_eff_fFWD))
  
  tab_r2_fFWD <- tibble::tibble(
    Metric = c("R2_marginal", "R2_conditional"),
    Value  = c(r2_fFWD$R2_marginal, r2_fFWD$R2_conditional)
  )
  
  cat("\n=== Alpha model (fFWD; richness q0) ===\n")
  cat(sprintf("Chosen family: %s\n", family(best_fFWD)$family))
  print(an_tab_fFWD)
  print(r2_fFWD)
  
  tab_q0_fix_fFWD <- broom.mixed::tidy(
    best_fFWD,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      term = dplyr::recode(
        term,
        `(Intercept)` = "Intercept",
        `log_reads`   = "log(reads+1)",
        `diameter_z`  = "Diameter (z)"
      ),
      `95% CI` = sprintf("%.3f to %.3f", conf.low, conf.high)
    ) %>%
    dplyr::transmute(
      Scope    = scope,
      Term     = term,
      IRR      = round(estimate, 3),
      `95% CI`,
      z        = round(statistic, 3),
      `p-value` = signif(p.value, 3)
    )
  
  print(tab_q0_fix_fFWD)
  readr::write_csv2(tab_q0_fix_fFWD, file.path(cfg$outdir, paste0("alpha_richness_IRR_", scope, ".csv")))
  
  get_re_sd <- function(mod, grp = "natman") {
    vc <- suppressWarnings(VarCorr(mod))
    if (!("cond" %in% names(vc))) return(NA_real_)
    vc_cond <- vc$cond
    if (!grp %in% names(vc_cond)) return(NA_real_)
    as.numeric(attr(vc_cond[[grp]], "stddev"))[1]
  }
  
  fam_str <- family(best_fFWD)$family
  AIC_val <- AIC(best_fFWD)
  re_sd   <- get_re_sd(best_fFWD, "natman")
  n_samp  <- nrow(df_fFWD)
  n_tree  <- dplyr::n_distinct(df_fFWD$natman)
  
  fit_block_fFWD <- tibble::tibble(
    `Model formula` = rlang::expr_text(formula(best_fFWD)),
    `Family` = fam_str,
    `AIC` = AIC_val,
    `Random effect SD (natman)` = re_sd,
    `R2 marginal` = unname(r2_fFWD$R2_marginal),
    `R2 conditional` = unname(r2_fFWD$R2_conditional),
    `N samples` = n_samp,
    `N trees` = n_tree
  )
  
  print(fit_block_fFWD)
  readr::write_csv2(fit_block_fFWD, file.path(cfg$outdir, paste0("alpha_richness_model_fit_", scope, ".csv")))
  
  emmeans::emmeans(best_fFWD, ~ decay_stage, type = "response")
  emmeans::emmeans(best_fFWD, ~ diameter_z, type = "response")
  emmeans::emmeans(best_fFWD, ~ umi, type = "response")
}

# -------------------------------
# 11) Outputs
# -------------------------------
get_term_val <- function(tbl, term, col) {
  v <- tbl[[col]][tbl$term == term]
  if (length(v) == 0) NA_real_ else v[1]
}

key_tbl <- tibble(
  metric = c("RobustAitchison", "Jaccard"),
  PERMANOVA_R2_unblocked = c(
    get_term_val(R2_unblocked_ra, "natman", "R2"),
    get_term_val(R2_unblocked_jac, "natman", "R2")
  ),
  PERMANOVA_R2adj_unblocked = c(
    get_term_val(R2_unblocked_ra, "natman", "R2adj"),
    get_term_val(R2_unblocked_jac, "natman", "R2adj")
  ),
  ANOSIM_R = c(an_ra$statistic, an_j$statistic),
  ECDF_AUC = c(AUC_ra, AUC_j),
  rank_biserial = c(rb_ra, rb_j)
)

write_csv(key_tbl, file.path(cfg$tabledir, "_KEY_RESULTS_FALLEN_FWD.csv"))

sink(file.path(cfg$tabledir, "_REPORT_SNIPPET.txt"))
cat("FALLEN FWD — key outcomes (natman as grouping):\n")
print(key_tbl)
cat("\nPERMANOVA (blocked by natman) — tables written to ", cfg$tabledir, "\n", sep = "")
sink()