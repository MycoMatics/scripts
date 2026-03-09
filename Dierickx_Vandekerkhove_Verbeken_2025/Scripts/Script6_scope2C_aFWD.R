# =============================================================================
# Script: Script6_scope2C_aFWD.R
# Project: Deadwood fungal community analyses
# Scope: aFWD / ATTACHED
# Purpose: Analyse attached fine woody debris communities using PERMANOVA,
#          blocked PERMANOVA, LOO stability, PERMDISP, natman null models
#          for R2 and AUC, ECDF summaries, ANOSIM, and alpha richness models
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - loo_permanova 
# - mk_within_between_tbl
# - boot_diff_median
# - wilcoxon_auc_once
# - overdisp_phi
# - check_glmm_lme4
# - plot_resid_glmm
# - ensure_dir
# - assert_objects
# - alpha_from_otu
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(forcats)
  library(vegan)
  library(permute)
  library(glmmTMB)
  library(emmeans)
  library(betapart)
  library(broom.mixed)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(performance)
  library(DHARMa)
  library(car)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c(
  "mk_within_between_tbl", "boot_diff_median", "wilcoxon_auc_once",
  "overdisp_phi", "check_glmm_lme4", "plot_resid_glmm",
  "ensure_dir", "alpha_from_otu"
))

setwd(".")

outdir <- "plots/FWD/ATTACHED"
ensure_dir(outdir)

# ---- Speed vs precision knobs -----------------------------------------------
N_PERM_OBS   <- 999
N_NULL_R2    <- 999
N_ADONIS_IN  <- 199
N_NULL_AUC   <- 5000
SEED_MAIN    <- 42

size_colors <- sc_colors
fwd_colors  <- c(ATTACHED = dw_colors[["aFWD"]])

# =============================================================================
# 1) DATA FILTERING
# =============================================================================
FWD_meta <- META1 %>%
  filter(dw_type2 == "aFWD") %>%
  filter(ds_at_drill != "0") %>% 
  filter(ds_at_drill != "5") %>% 
    mutate(
    decay_stage = factor(decay_stage, levels = c("EARLY", "AVERAGE", "LATE")),
    ds_at_drill = factor(ds_at_drill, levels = c( "1", "2", "3", "4")),
    size        = droplevels(factor(size))) %>%
  mutate(diameter_at_drill_z = as.numeric(diameter_at_drill_z)) %>% 
  filter(!is.na(diameter_at_drill))

otu_fwd <- otu_matrix_filt[FWD_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fwd), FWD_meta$sample))
otu_fwd <- otu_fwd[, colSums(otu_fwd) > 0, drop = FALSE]

FWD_meta <- FWD_meta %>%
  mutate(
    reads = rowSums(otu_fwd),
    log_reads = log1p(reads)
  ) %>%
  filter(reads > 0) %>%
  droplevels()

otu_fwd <- otu_fwd[FWD_meta$sample, , drop = FALSE]

cat(
  "ATTACHED FWD samples:", nrow(FWD_meta),
  " | OTUs:", ncol(otu_fwd),
  " | trees: ", length(unique(FWD_meta$natman)), "\n"
)

# =============================================================================
# 2) DISTANCES
# =============================================================================
D_att <- vegdist(otu_fwd, method = "robust.aitchison")
D_att <- fix_dist_labels(D_att, samps = rownames(otu_fwd))

otu_pa <- (otu_fwd > 0) * 1L
D_att_jac <- vegdist(otu_pa, method = "jaccard", binary = TRUE)

# =============================================================================
# 3) OBSERVED PERMANOVA
# =============================================================================
set.seed(SEED_MAIN)
perm_unblocked <- how(nperm = N_PERM_OBS)

ado_att_Ait <- adonis2(
  D_att ~ log_reads + umi + natman + ds_at_drill + diameter_at_drill_z,
  data = FWD_meta,
  by = "terms",
  permutations = perm_unblocked
)

ado_att_Jac <- adonis2(
  D_att_jac ~ log_reads + umi + natman + ds_at_drill + diameter_at_drill_z,
  data = FWD_meta,
  by = "terms",
  permutations = perm_unblocked
)

R2_obs_Ait <- as.data.frame(ado_att_Ait)["natman", "R2"]
R2_obs_Jac <- as.data.frame(ado_att_Jac)["natman", "R2"]

print(ado_att_Ait)
print("natman R² value Robust Aitchison")
print(R2_obs_Ait)

print(ado_att_Jac)
print("natman R² value Jaccard")
print(R2_obs_Jac)

# ---- Blocked permutations by natman -----------------------------------------
set.seed(SEED_MAIN)
nperm <- N_PERM_OBS
perm_blocked <- how(nperm = N_PERM_OBS)
setBlocks(perm_blocked) <- FWD_meta$natman

ado_att_Ait_blk <- adonis2(
  D_att ~ log_reads + umi + ds_at_drill + diameter_at_drill_z,
  data = FWD_meta,
  by = "margin",
  permutations = perm_blocked
)

ado_att_Jac_blk <- adonis2(
  D_att_jac ~ log_reads + umi + ds_at_drill + diameter_at_drill_z,
  data = FWD_meta,
  by = "margin",
  permutations = perm_blocked
)

print(ado_att_Ait_blk)
print(ado_att_Jac_blk)

# ---- LOO PERMANOVA -----------------------------------------------------------
loo_permanova(
  X = otu_fwd,
  meta = FWD_meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + ds_at_drill + diameter_at_drill_z",
  dist_method = "robust.aitchison",
  keep_terms = c("ds_at_drill", "log_reads", "umi", "diameter_at_drill_z"),
  blocked_by = "natman",
  by_TYPE = "margin",
  nperm_outer = 499,
  out_csv = file.path(outdir, "loo_attached_decay_Aitchison.csv")
)

loo_permanova(
  X = otu_fwd,
  meta = FWD_meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + ds_at_drill + diameter_at_drill_z",
  dist_method = "jaccard",
  keep_terms = c("ds_at_drill", "log_reads", "umi", "diameter_at_drill_z"),
  blocked_by = "natman",
  by_TYPE = "margin",
  nperm_outer = 499,
  out_csv = file.path(outdir, "loo_attached_decay_Jaccard.csv")
)

# =============================================================================
# 4) PERMDISP
# =============================================================================
set.seed(SEED_MAIN)

# Aitchison
bd_as_Ait <- betadisper(D_att, FWD_meta$natman)
anova(bd_as_Ait); permutest(bd_as_Ait, permutations = nperm); TukeyHSD(bd_as_Ait)
boxplot(bd_as_Ait, xlab = "natman", main = "Aitchison: distance to group centroid")

bd_as_Ait <- betadisper(D_att, FWD_meta$umi)
anova(bd_as_Ait); permutest(bd_as_Ait, permutations = nperm); TukeyHSD(bd_as_Ait)
boxplot(bd_as_Ait, xlab = "umi", main = "Aitchison: distance to group centroid")

bd_as_Ait <- betadisper(D_att, FWD_meta$ds_at_drill)
anova(bd_as_Ait); permutest(bd_as_Ait, permutations = nperm); TukeyHSD(bd_as_Ait)
boxplot(bd_as_Ait, xlab = "ds_at_drill", main = "Aitchison: distance to group centroid")

# Jaccard
bd_as_Jacc <- betadisper(D_att_jac, FWD_meta$natman)
anova(bd_as_Jacc); permutest(bd_as_Jacc, permutations = nperm); TukeyHSD(bd_as_Jacc)
boxplot(bd_as_Jacc, xlab = "natman", main = "Jaccard: distance to group centroid")

bd_as_Jacc <- betadisper(D_att_jac, FWD_meta$umi)
anova(bd_as_Jacc); permutest(bd_as_Jacc, permutations = nperm); TukeyHSD(bd_as_Jacc)
boxplot(bd_as_Jacc, xlab = "umi", main = "Jaccard: distance to group centroid")

bd_as_Jacc <- betadisper(D_att_jac, FWD_meta$ds_at_drill)
anova(bd_as_Jacc); permutest(bd_as_Jacc, permutations = nperm); TukeyHSD(bd_as_Jacc)
boxplot(bd_as_Jacc, xlab = "ds_at_drill", main = "Jaccard: distance to group centroid")

bd_nat_Ait <- betadisper(D_att, FWD_meta$natman, type = "centroid")
anova(bd_nat_Ait)

bd_nat_Jac <- betadisper(D_att_jac, FWD_meta$natman, type = "centroid")
anova(bd_nat_Jac)

pd_Ait <- permutest(bd_nat_Ait, permutations = N_PERM_OBS)
pd_Jac <- permutest(bd_nat_Jac, permutations = N_PERM_OBS)
print(pd_Ait)
print(pd_Jac)

tuk_Ait <- TukeyHSD(bd_nat_Ait)
tuk_Jac <- TukeyHSD(bd_nat_Jac)
print(tuk_Ait)
print(tuk_Jac)

capture.output(pd_Ait, file = file.path(outdir, "PERMDISP_natman_Aitchison.txt"))
capture.output(pd_Jac, file = file.path(outdir, "PERMDISP_natman_Jaccard.txt"))
capture.output(tuk_Ait, file = file.path(outdir, "PERMDISP_Tukey_natman_Aitchison.txt"))
capture.output(tuk_Jac, file = file.path(outdir, "PERMDISP_Tukey_natman_Jaccard.txt"))

# =============================================================================
# 5) R² NULLS — DESIGN-PRESERVING GROUP SIZES
# =============================================================================
set.seed(SEED_MAIN)
gs <- table(FWD_meta$natman)

shuffle_grp <- function() {
  factor(rep(seq_along(gs), times = gs)[sample(nrow(FWD_meta))])
}

R2_null_Ait <- replicate(N_NULL_R2, {
  grp <- shuffle_grp()
  meta_tmp <- FWD_meta %>% mutate(grp = grp)
  as.data.frame(
    adonis2(
      D_att ~ log_reads + umi + grp + ds_at_drill + diameter_at_drill_z,
      data = meta_tmp,
      by = "terms",
      permutations = N_ADONIS_IN
    )
  )["grp", "R2"]
})

saveRDS(R2_null_Ait, file.path(outdir, "aFWD_R2_null_Ait.rds"))

R2_null_Jac <- replicate(N_NULL_R2, {
  grp <- shuffle_grp()
  meta_tmp <- FWD_meta %>% mutate(grp = grp)
  as.data.frame(
    adonis2(
      D_att_jac ~ log_reads + umi + grp + ds_at_drill + diameter_at_drill_z,
      data = meta_tmp,
      by = "terms",
      permutations = N_ADONIS_IN
    )
  )["grp", "R2"]
})

saveRDS(R2_null_Jac, file.path(outdir, "aFWD_R2_null_Jac.rds"))

ci_Ait   <- quantile(R2_null_Ait, c(.025, .975))
ci_Jac   <- quantile(R2_null_Jac, c(.025, .975))
pemp_Ait <- mean(R2_null_Ait >= R2_obs_Ait)
pemp_Jac <- mean(R2_null_Jac >= R2_obs_Jac)
ses_Ait  <- (R2_obs_Ait - mean(R2_null_Ait)) / sd(R2_null_Ait)
ses_Jac  <- (R2_obs_Jac - mean(R2_null_Jac)) / sd(R2_null_Jac)
perc_Ait <- mean(R2_null_Ait <= R2_obs_Ait)
perc_Jac <- mean(R2_null_Jac <= R2_obs_Jac)

R2_nat_AFWD_RA <- ggplot(tibble(R2 = R2_null_Ait), aes(R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ci_Ait, linetype = 3) +
  labs(
    title = "ATTACHED – natman R² null vs observed (robust Aitchison)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      R2_obs_Ait, mean(R2_null_Ait), ci_Ait[1], ci_Ait[2], pemp_Ait, ses_Ait
    ),
    x = expression(R^2), y = "Density"
  ) +
  theme_bw()

print(R2_nat_AFWD_RA)
ggsave(file.path(outdir, "natman_R2_null_vs_observed_ATTACHED_Aitchison.png"), R2_nat_AFWD_RA, width = 6, height = 4, dpi = 300)

R2_nat_AFWD_JAC <- ggplot(tibble(R2 = R2_null_Jac), aes(R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ci_Jac, linetype = 3) +
  labs(
    title = "ATTACHED – natman R² null vs observed (Jaccard)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, 95%%CI=[%.3f,%.3f], p=%.3f, SES=%.2f",
      R2_obs_Jac, mean(R2_null_Jac), ci_Jac[1], ci_Jac[2], pemp_Jac, ses_Jac
    ),
    x = expression(R^2), y = "Density"
  ) +
  theme_bw()

print(R2_nat_AFWD_JAC)
ggsave(file.path(outdir, "natman_R2_null_vs_observed_ATTACHED_Jaccard.png"), R2_nat_AFWD_JAC, width = 6, height = 4, dpi = 300)

# =============================================================================
# 6) WITHIN vs BETWEEN DISTANCES
# =============================================================================
DmA <- as.matrix(D_att)
U <- which(upper.tri(DmA), arr.ind = TRUE)
DmJ <- as.matrix(D_att_jac)

pairs_A <- mk_within_between_tbl(DmA, FWD_meta$natman)
pairs_J <- mk_within_between_tbl(DmJ, FWD_meta$natman)

sum_A <- pairs_A %>% group_by(type) %>% summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop")
sum_J <- pairs_J %>% group_by(type) %>% summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop")

bd_A <- boot_diff_median(pairs_A, seed = SEED_MAIN)
ci_diff_A <- quantile(bd_A, c(.025, .975))
bd_J <- boot_diff_median(pairs_J, seed = SEED_MAIN)
ci_diff_J <- quantile(bd_J, c(.025, .975))

est_A <- with(sum_A, median[type == "between"] - median[type == "within"])
est_J <- with(sum_J, median[type == "between"] - median[type == "within"])

annot_ecdf <- function(dat, title_txt, xlab_txt, est, ci) {
  med_w <- median(dat$dist[dat$type == "within"])
  med_b <- median(dat$dist[dat$type == "between"])
  subtitle <- sprintf(
    "median(between)=%.3f  median(within)=%.3f  diff=%.3f  95%%CI=[%.3f, %.3f]",
    med_b, med_w, est, ci[1], ci[2]
  )
  ggplot(dat, aes(dist, color = type)) +
    stat_ecdf(linewidth = 0.7) +
    theme_bw() +
    labs(title = title_txt, subtitle = subtitle, x = xlab_txt, y = "ECDF")
}

p_ecdf_A <- annot_ecdf(
  pairs_A,
  "ATTACHED (Aitchison/Eucl.): ECDF of pairwise distances",
  "Aitchison/Euclidean distance",
  est_A, ci_diff_A
)
print(p_ecdf_A)

p_ecdf_J <- annot_ecdf(
  pairs_J,
  "ATTACHED (Jaccard): ECDF of pairwise distances",
  "Jaccard distance",
  est_J, ci_diff_J
)
print(p_ecdf_J)

ggsave(file.path(outdir, "ECDF_within_between_ATTACHED_Jaccard.png"), p_ecdf_J, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "ECDF_within_between_ATTACHED_Ait.png"), p_ecdf_A, width = 6, height = 4, dpi = 300)

# =============================================================================
# 7) AUC + DESIGN-PRESERVING AUC NULLS
# =============================================================================
upper_idx <- U
AUC_obs_Ait <- wilcoxon_auc_once(DmA, FWD_meta$natman, upper_idx)
AUC_obs_Jac <- wilcoxon_auc_once(DmJ, FWD_meta$natman, upper_idx)

set.seed(SEED_MAIN)
AUC_null_Ait <- replicate(N_NULL_AUC, wilcoxon_auc_once(DmA, shuffle_grp(), upper_idx))
AUC_null_Jac <- replicate(N_NULL_AUC, wilcoxon_auc_once(DmJ, shuffle_grp(), upper_idx))

ciAUC_Ait <- quantile(AUC_null_Ait, c(.025, .975))
pAUC_Ait <- mean(AUC_null_Ait >= AUC_obs_Ait)
ciAUC_Jac <- quantile(AUC_null_Jac, c(.025, .975))
pAUC_Jac <- mean(AUC_null_Jac >= AUC_obs_Jac)

ses_AUC_Ait <- (AUC_obs_Ait - mean(AUC_null_Ait)) / sd(AUC_null_Ait)
ses_AUC_Jac <- (AUC_obs_Jac - mean(AUC_null_Jac)) / sd(AUC_null_Jac)

rb_Ait <- 2 * AUC_obs_Ait - 1
rb_Ait_CI <- 2 * ciAUC_Ait - 1
rb_Jac <- 2 * AUC_obs_Jac - 1
rb_Jac_CI <- 2 * ciAUC_Jac - 1

AUC_nat_AFWD_RA <- ggplot(tibble(AUC = AUC_null_Ait), aes(AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_Ait, linetype = 3) +
  labs(
    title = "AUC null vs observed (ATTACHED, Aitchison)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, p=%.3f, SES=%.2f",
      AUC_obs_Ait, mean(AUC_null_Ait), pAUC_Ait, ses_AUC_Ait
    )
  ) +
  theme_bw()

print(AUC_nat_AFWD_RA)
ggsave(file.path(outdir, "AUC_null_vs_observed_ATTACHED_Aitchison.png"), AUC_nat_AFWD_RA, width = 6, height = 4, dpi = 300)

AUC_nat_AFWD_JAC <- ggplot(tibble(AUC = AUC_null_Jac), aes(AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_Jac, linetype = 3) +
  labs(
    title = "AUC null vs observed (ATTACHED, Jaccard)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, p=%.3f, SES=%.2f",
      AUC_obs_Jac, mean(AUC_null_Jac), pAUC_Jac, ses_AUC_Jac
    )
  ) +
  theme_bw()

print(AUC_nat_AFWD_JAC)
ggsave(file.path(outdir, "AUC_null_vs_observed_ATTACHED_Jaccard.png"), AUC_nat_AFWD_JAC, width = 6, height = 4, dpi = 300)

# =============================================================================
# 8) ANOSIM
# =============================================================================
set.seed(SEED_MAIN)
anos_Ait <- anosim(D_att, grouping = FWD_meta$natman, permutations = 9999)
anos_Jac <- anosim(D_att_jac, grouping = FWD_meta$natman, permutations = 9999)
print(anos_Ait)
print(anos_Jac)

# =============================================================================
# 9) ALPHA DIVERSITY MODELS (q0 only)
# =============================================================================
scope <- "aFWD"
cat("\n[ ALPHA RICHNESS GLMM - aFWD ]\n")

rich_tbl <- tibble(
  sample   = rownames(otu_fwd),
  richness = rowSums(as.matrix(otu_fwd) > 0, na.rm = TRUE)
)

df_aFWD <- FWD_meta %>%
  left_join(rich_tbl, by = "sample") %>%
  mutate(
    diameter_z  = diameter_at_drill_z,
    ds_at_drill = droplevels(ds_at_drill),
    log_reads   = if_else(is.finite(log_reads), log_reads, log1p(reads))
  ) %>%
  filter(!is.na(richness), reads > 0)

rhs_terms_glmm_aFWD <- c("log_reads", "umi", "(1|natman)", "ds_at_drill", "diameter_z")
form_aFWD <- as.formula(paste("richness ~", paste(rhs_terms_glmm_aFWD, collapse = " + ")))
cat("GLMM form (aFWD):\n")
print(form_aFWD)

if (nlevels(df_aFWD$ds_at_drill) >= 2) {
  ctrl_aFWD <- lme4::glmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
  
  m_pois_aFWD <- suppressWarnings(
    lme4::glmer(
      form_aFWD,
      data = df_aFWD,
      family = poisson(link = "log"),
      control = ctrl_aFWD
    )
  )
  
  phi_p_aFWD  <- overdisp_phi(m_pois_aFWD)
  use_nb_aFWD <- is.finite(phi_p_aFWD$phi) && (phi_p_aFWD$phi > 1.2)
  
  if (use_nb_aFWD) {
    m_nb_aFWD <- suppressWarnings(
      lme4::glmer.nb(form_aFWD, data = df_aFWD, control = ctrl_aFWD)
    )
    best_aFWD <- if (AIC(m_pois_aFWD) + 2 < AIC(m_nb_aFWD)) m_pois_aFWD else m_nb_aFWD
    best_type_aFWD <- if (identical(best_aFWD, m_pois_aFWD)) "Poisson (glmer)" else "NB (glmer.nb)"
  } else {
    best_aFWD <- m_pois_aFWD
    best_type_aFWD <- "Poisson (glmer)"
  }
  
  cat(sprintf(
    "Selected aFWD alpha model: %s  Overdispersion Poisson phi = %.2f  AIC(best) = %.1f\n",
    best_type_aFWD, phi_p_aFWD$phi, AIC(best_aFWD)
  ))
  print(summary(best_aFWD))
  
  check_glmm_lme4(best_aFWD, name = "aFWD alpha")
  plot_resid_glmm(best_aFWD, file.path(outdir, "alpha_richness_aFWD_residuals.png"))
  
  an_tab_aFWD <- car::Anova(best_aFWD, type = 2)
  r2_aFWD <- performance::r2(best_aFWD)
  
  irr_aFWD <- broom.mixed::tidy(
    best_aFWD,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      term = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "log_reads" ~ "log(reads+1)",
        term == "umi1" ~ "Dataset: UMI vs raw",
        term == "diameter_z" ~ "Diameter (z)",
        grepl("^ds_at_drill", term) ~ paste0(
          "Decay ", sub("^ds_at_drill", "", term), " vs ",
          levels(df_aFWD$ds_at_drill)[1]
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
  
  emm_ds_aFWD <- emmeans::emmeans(best_aFWD, ~ ds_at_drill, type = "response")
  ds_df_aFWD <- summary(emm_ds_aFWD) %>% as.data.frame()
  ds_df_aFWD$ds_at_drill <- factor(ds_df_aFWD$ds_at_drill, ordered = TRUE)
  
  p_emm_ds_aFWD <- ggplot2::ggplot(
    ds_df_aFWD,
    ggplot2::aes(x = ds_at_drill, y = response)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = asymp.LCL, ymax = asymp.UCL),
      width = 0.2
    ) +
    ggplot2::labs(
      x = "Decay at drill (ordered)",
      y = "Estimated richness",
      title = paste0("Richness ~ decay (", scope, ")")
    ) +
    ggplot2::theme_classic(base_size = 11)
  
  ggplot2::ggsave(
    file.path(outdir, paste0("alpha_richness_emm_DS_", scope, ".png")),
    p_emm_ds_aFWD, width = 6.2, height = 3.2, dpi = 300
  )
  print(p_emm_ds_aFWD)
  
  keep_eff_aFWD <- c("ds_at_drill", "diameter_z", "log_reads")
  tab_anova_aFWD <- tibble::tibble(
    Effect  = rownames(an_tab_aFWD),
    Chisq   = an_tab_aFWD$`Chisq`,
    df      = an_tab_aFWD$Df,
    p       = an_tab_aFWD$`Pr(>Chisq)`
  ) %>%
    dplyr::filter(Effect %in% keep_eff_aFWD) %>%
    dplyr::mutate(p_BH = p.adjust(p, "BH")) %>%
    dplyr::arrange(match(Effect, keep_eff_aFWD))
  
  tab_r2_aFWD <- tibble::tibble(
    Metric = c("R2_marginal", "R2_conditional"),
    Value  = c(r2_aFWD$R2_marginal, r2_aFWD$R2_conditional)
  )
  
  cat("\n=== Alpha model (aFWD; richness q0) ===\n")
  cat(sprintf("Chosen family: %s\n", family(best_aFWD)$family))
  print(an_tab_aFWD)
  print(r2_aFWD)
  
  tab_q0_fix_aFWD <- broom.mixed::tidy(
    best_aFWD,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    dplyr::mutate(
      term = dplyr::recode(
        term,
        `(Intercept)` = "Intercept",
        `log_reads`   = "log(reads+1)",
        `umi1`        = "Dataset: UMI vs raw",
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
  
  print(tab_q0_fix_aFWD)
  readr::write_csv2(tab_q0_fix_aFWD, file.path(outdir, paste0("alpha_richness_IRR_", scope, ".csv")))
  
  get_re_sd <- function(mod, grp = "natman") {
    vc <- suppressWarnings(VarCorr(mod))
    if (!("cond" %in% names(vc))) return(NA_real_)
    vc_cond <- vc$cond
    if (!grp %in% names(vc_cond)) return(NA_real_)
    as.numeric(attr(vc_cond[[grp]], "stddev"))[1]
  }
  
  fam_str <- family(best_aFWD)$family
  AIC_val <- AIC(best_aFWD)
  re_sd   <- get_re_sd(best_aFWD, "natman")
  n_samp  <- nrow(df_aFWD)
  n_tree  <- dplyr::n_distinct(df_aFWD$natman)
  
  fit_block_aFWD <- tibble::tibble(
    `Model formula` = rlang::expr_text(formula(best_aFWD)),
    `Family` = fam_str,
    `AIC` = AIC_val,
    `Random effect SD (natman)` = re_sd,
    `R2 marginal` = unname(r2_aFWD$R2_marginal),
    `R2 conditional` = unname(r2_aFWD$R2_conditional),
    `N samples` = n_samp,
    `N trees` = n_tree
  )
  
  print(fit_block_aFWD)
  readr::write_csv2(fit_block_aFWD, file.path(outdir, paste0("alpha_richness_model_fit_", scope, ".csv")))
  
  emmeans::emmeans(best_aFWD, ~ ds_at_drill, type = "response")
  emmeans::emmeans(best_aFWD, ~ umi, type = "response")
  emmeans::emmeans(best_aFWD, ~ diameter_z, type = "response")
}

# =============================================================================
# 10) SUMMARY CSV & COMPACT RDS BUNDLE
# =============================================================================
summary_both <- tibble(
  state = "ATTACHED",
  metric = c("robust.aitchison", "jaccard"),
  R2_obs_natman = c(R2_obs_Ait, R2_obs_Jac),
  R2_null_mean  = c(mean(R2_null_Ait), mean(R2_null_Jac)),
  R2_null_CI_lo = c(ci_Ait[1], ci_Jac[1]),
  R2_null_CI_hi = c(ci_Ait[2], ci_Jac[2]),
  p_empirical   = c(pemp_Ait, pemp_Jac),
  R2_design_extreme = c(
    R2_obs_Ait > ci_Ait[2] | R2_obs_Ait < ci_Ait[1],
    R2_obs_Jac > ci_Jac[2] | R2_obs_Jac < ci_Jac[1]
  ),
  R2_percentile = c(perc_Ait, perc_Jac),
  AUC_obs       = c(AUC_obs_Ait, AUC_obs_Jac),
  AUC_null_mean = c(mean(AUC_null_Ait), mean(AUC_null_Jac)),
  AUC_null_CI_lo = c(ciAUC_Ait[1], ciAUC_Jac[1]),
  AUC_null_CI_hi = c(ciAUC_Ait[2], ciAUC_Jac[2]),
  AUC_p_emp      = c(pAUC_Ait, pAUC_Jac),
  AUC_SES        = c(ses_AUC_Ait, ses_AUC_Jac),
  ANOSIM_R       = c(anos_Ait$statistic, anos_Jac$statistic),
  ANOSIM_p       = c(anos_Ait$signif, anos_Jac$signif),
  dist_within_median  = c(sum_A$median[sum_A$type == "within"], sum_J$median[sum_J$type == "within"]),
  dist_between_median = c(sum_A$median[sum_A$type == "between"], sum_J$median[sum_J$type == "between"]),
  dist_diff      = c(est_A, est_J),
  diff_CI_lo     = c(ci_diff_A[1], ci_diff_J[1]),
  diff_CI_hi     = c(ci_diff_A[2], ci_diff_J[2])
)

readr::write_csv(summary_both, file.path(outdir, "natman_evidence_ATTACHED_full.csv"))
print(summary_both, n = Inf)

att_results <- list(
  meta = FWD_meta,
  otu = otu_fwd,
  distances = list(Ait = D_att, Jac = D_att_jac),
  permanova = list(Ait = ado_att_Ait, Jac = ado_att_Jac),
  nulls = list(
    R2_Ait = R2_null_Ait,
    R2_Jac = R2_null_Jac,
    AUC_Ait = AUC_null_Ait,
    AUC_Jac = AUC_null_Jac
  ),
  anosim = list(Ait = anos_Ait, Jac = anos_Jac),
  alpha = list(
    q0 = list(
      model = if (exists("best_aFWD")) best_aFWD else NULL,
      anova = if (exists("an_tab_aFWD")) an_tab_aFWD else NULL,
      r2 = if (exists("r2_aFWD")) r2_aFWD else NULL,
      emm_ds = if (exists("ds_df_aFWD")) ds_df_aFWD else NULL,
      fixed = if (exists("irr_aFWD")) irr_aFWD else NULL
    )
  ),
  within_between = list(Ait = sum_A, Jac = sum_J)
)

saveRDS(att_results, file.path(outdir, "ATTACHED_FWD_results.rds"))
