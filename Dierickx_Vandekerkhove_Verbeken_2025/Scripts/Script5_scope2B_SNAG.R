# =============================================================================
# Script: Script5_scope2B_SNAG.R
# Project: Deadwood fungal community analyses
# Scope: SNAG
# Purpose: Analyse paired north-south snag samples using PERMANOVA, LOO
#          stability, PERMDISP, within-vs-between distance summaries,
#          null models for AUC and R2, betapart decomposition, and alpha
#          richness modelling
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - alpha_from_otu
# - fix_dist_labels
# - mk_within_between_tbl
# - boot_diff_median
# - boot_diff_ci
# - cliffs_delta
# - wilcoxon_auc_once
# - overdisp_phi
# - ensure_dir
# - assert_objects
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
  library(performance)
  library(broom.mixed)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(betapart)
  library(indicspecies)
})

assert_objects(c("META1", "otu_matrix_filt"))
assert_objects(c(
  "alpha_from_otu", "fix_dist_labels", "mk_within_between_tbl",
  "boot_diff_median", "boot_diff_ci", "cliffs_delta",
  "wilcoxon_auc_once", "overdisp_phi", "ensure_dir"
))

SEED_GLOBAL <- 42L
set.seed(SEED_GLOBAL)
setwd(".")

outdir <- "plots/SNAG"
ensure_dir(outdir)

# ---- Output files -------------------------------------------------------------
f_natman_r2_csv   <- file.path(outdir, "permanova_natman_R2_across_distances.csv")
f_summary_csv     <- file.path(outdir, "natman_evidence_SNAG_full.csv")
f_betapart_pairs  <- file.path(outdir, "betapart_pairs_within_between.csv")
f_betapart_sum    <- file.path(outdir, "betapart_within_between_summary_combined.csv")
f_rds             <- file.path(outdir, "SNAG_results.rds")

# ---- Speed vs precision knobs ------------------------------------------------
NPERM_MAIN    <- 999
NPERM_BLOCKED <- NPERM_MAIN
NPERM_LOO     <- 499
PERM_TERM     <- 199
B_BOOT        <- 999
B_NULL_R2     <- B_BOOT
ANOSIM_PERM   <- 9999
TRYMAX_NMDS   <- 999

# ==============================================================================
# 1) Helper functions
# ==============================================================================
loo_permanova <- function(X, meta, group_id = "natman",
                          formula_rhs = "log_reads + umi + ds_at_drill + aspect",
                          dist_method = "robust.aitchison",
                          nperm_outer = NPERM_LOO,
                          blocked_by = NULL,
                          by_TYPE = "terms",
                          keep_terms = c("aspect", "ds_at_drill"),
                          out_csv = NULL) {
  stopifnot(nrow(X) == nrow(meta), identical(rownames(X), meta$sample))
  
  groups <- droplevels(factor(meta[[group_id]]))
  trees  <- levels(groups)
  
  print(length(trees))
  print("trees")
  
  out <- purrr::map_dfr(trees, function(t) {
    keep <- groups != t & !is.na(groups)
    if (sum(keep) < 5) return(NULL)
    
    Dk <- vegdist(X[keep, , drop = FALSE], method = dist_method)
    data_k <- droplevels(meta[keep, , drop = FALSE])
    
    permk <- how(nperm = nperm_outer)
    if (!is.null(blocked_by)) setBlocks(permk) <- droplevels(data_k[[blocked_by]])
    
    fml <- reformulate(
      termlabels = all.vars(stats::terms(stats::as.formula(paste("~", formula_rhs)))),
      response = NULL
    )
    
    a2 <- adonis2(
      Dk ~ .,
      data = data_k[, all.vars(stats::terms(fml)), drop = FALSE],
      permutations = permk,
      by = by_TYPE
    )
    
    as_tibble(as.data.frame(a2), rownames = "term") |>
      filter(term %in% keep_terms) |>
      transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  
  if (nrow(out)) {
    summ <- out |>
      group_by(term) |>
      summarise(
        median_R2 = median(R2, na.rm = TRUE),
        p05 = quantile(p, 0.05, na.rm = TRUE),
        p95 = quantile(p, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    print(summ)
    
    if (!is.null(out_csv)) {
      readr::write_csv2(out, out_csv)
      readr::write_csv2(summ, sub("\\.csv$", "_summary.csv", out_csv))
    }
  }
  
  invisible(out)
}

plot_nmds <- function(nmds_obj, meta, title = "NMDS", color = "decay_stage", shape = "aspect") {
  df <- scores(nmds_obj, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample")
  
  lab <- paste0("Stress = ", round(nmds_obj$stress, 3))
  
  ggplot(df, aes(NMDS1, NMDS2, color = .data[[color]], shape = .data[[shape]])) +
    geom_point(size = 2.7, alpha = 0.85) +
    theme_minimal(base_size = 13) +
    annotate(
      "text",
      x = min(df$NMDS1, na.rm = TRUE),
      y = max(df$NMDS2, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      label = lab,
      size = 4.2,
      fontface = "italic"
    ) +
    labs(title = title, color = color, shape = shape)
}

safe_R2 <- function(mod) {
  out <- try(vegan::RsquareAdj(mod), silent = TRUE)
  if (inherits(out, "try-error")) {
    if (inherits(mod, "cca")) {
      fun <- getFromNamespace("RsquareAdj.cca", "vegan")
      out <- fun(mod)
    } else {
      warning("safe_R2: object is not a 'cca'/'capscale'; returning NA")
      return(c(r2 = NA_real_, r2_adj = NA_real_))
    }
  }
  c(r2_raw = out$r.squared, r2_adj = out$adj.r.squared)
}

make_R2_table <- function(mod_all, mod_tree, mod_micro, method) {
  R2_all   <- safe_R2(mod_all)
  R2_tree  <- safe_R2(mod_tree)
  R2_micro <- safe_R2(mod_micro)
  
  tibble(
    method       = method,
    R2_all_raw   = R2_all["r2_raw"],
    R2_all_adj   = R2_all["r2_adj"],
    R2_tree_adj  = R2_tree["r2_adj"],
    R2_micro_adj = R2_micro["r2_adj"],
    R2_shared    = R2_all["r2_adj"] - R2_tree["r2_adj"] - R2_micro["r2_adj"],
    R2_unexpl    = 1 - R2_all["r2_adj"]
  )
}

# ==============================================================================
# 2) Metadata: SNAG paired aspects
# ==============================================================================
SNAG_meta <- META1 %>%
  filter(dw_type2 == "SNAG") %>%
  mutate(
    aspect = case_when(
      position_2 %in% c("NORTH", "N") ~ "NORTH",
      position_2 %in% c("SOUTH", "S") ~ "SOUTH",
      TRUE ~ NA_character_
    ),
    decay_stage = case_when(
      sample %in% c("P1_25", "P1_26", "P2_68", "P2_72") ~ "AVERAGE",
      sample %in% c("P2_04", "P2_10") ~ "EARLY",
      sample %in% c("MI_34", "MI_35") ~ "LATE",
      TRUE ~ decay_stage
    ),
    aspect      = droplevels(factor(aspect, levels = c("NORTH", "SOUTH"))),
    decay_stage = droplevels(factor(decay_stage)),
    ds_at_drill = droplevels(factor(as.character(ds_at_drill), levels = c("0", "1", "2", "3", "4", "5"))),
    natman      = factor(gsub('"', "", as.character(natman)), ordered = FALSE)
  ) %>%
  filter(sample != "P3_33", natman != "ZF303") %>%
  droplevels()

# ==============================================================================
# 3) Align OTU table and read depth
# ==============================================================================
otu_snag <- otu_matrix_filt[SNAG_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_snag), SNAG_meta$sample))
otu_snag <- otu_snag[, colSums(otu_snag) > 0, drop = FALSE]

SNAG_meta <- SNAG_meta %>%
  mutate(
    reads = rowSums(otu_snag),
    log_reads = log1p(reads)
  ) %>%
  filter(!is.na(aspect), reads > 0) %>%
  droplevels()

otu_snag <- otu_snag[SNAG_meta$sample, , drop = FALSE]
cat("SNAG samples:", nrow(SNAG_meta), " | OTUs:", ncol(otu_snag), "\n")

# ==============================================================================
# 4) Distances
# ==============================================================================
D_robAit <- vegdist(otu_snag, method = "robust.aitchison")
D_robAit <- fix_dist_labels(D_robAit, samps = rownames(otu_snag))

otu_pa <- (otu_snag > 0) * 1L
D_jacc <- vegdist(otu_pa, method = "jaccard", binary = TRUE)
D_jacc <- fix_dist_labels(D_jacc, samps = rownames(otu_pa))

# ==============================================================================
# 5) Descriptive alpha-diversity boxplots
# ==============================================================================
alpha_snag <- tibble(
  sample   = rownames(otu_snag),
  richness = vegan::specnumber(otu_snag),
  shannon  = vegan::diversity(otu_snag, index = "shannon"),
  simpson  = vegan::diversity(otu_snag, index = "simpson")
) %>%
  left_join(SNAG_meta, by = "sample") %>%
  pivot_longer(c(richness, shannon, simpson), names_to = "metric", values_to = "value")

if (nlevels(SNAG_meta$aspect) > 1) {
  p_alpha_aspect <- alpha_snag %>%
    ggplot(aes(aspect, value, fill = aspect)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = aspect_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha by SNAG aspect (DBH)", x = "Aspect", y = "Value", fill = "Aspect")
  
  ggsave(file.path(outdir, "alpha_by_aspect.png"), p_alpha_aspect, width = 9, height = 4.8, dpi = 300)
  print(p_alpha_aspect)
}

if (nlevels(SNAG_meta$decay_stage) > 1) {
  p_alpha_ds <- alpha_snag %>%
    ggplot(aes(decay_stage, value, fill = decay_stage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = DS_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha by SNAG decay_stage", x = "decay_stage", y = "Value", fill = "decay_stage")
  
  ggsave(file.path(outdir, "alpha_by_decay_stage.png"), p_alpha_ds, width = 9, height = 4.8, dpi = 300)
  print(p_alpha_ds)
}

if (nlevels(SNAG_meta$ds_at_drill) > 1) {
  p_alpha_dsd <- alpha_snag %>%
    ggplot(aes(ds_at_drill, value, fill = ds_at_drill)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = ds_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha by SNAG ds_at_drill", x = "ds_at_drill", y = "Value", fill = "ds_at_drill")
  
  ggsave(file.path(outdir, "alpha_by_ds_at_drill.png"), p_alpha_dsd, width = 9, height = 4.8, dpi = 300)
  print(p_alpha_dsd)
}

# ==============================================================================
# 6) PERMANOVA
# ==============================================================================
set.seed(SEED_GLOBAL)
nperm <- NPERM_MAIN
perm_unblocked <- how(nperm = nperm)

perma_unblocked <- list(
  robAit = adonis2(
    D_robAit ~ log_reads + umi + natman + ds_at_drill + aspect,
    data = SNAG_meta,
    by = "terms",
    permutations = perm_unblocked
  ),
  Jacc = adonis2(
    D_jacc ~ log_reads + umi + natman + ds_at_drill + aspect,
    data = SNAG_meta,
    by = "terms",
    permutations = perm_unblocked
  )
)

cat("\nUnblocked PERMANOVA:\n")
purrr::iwalk(perma_unblocked, ~{
  cat("\n--", .y, "--\n")
  print(.x)
})

natman_R2 <- purrr::map_dfr(
  perma_unblocked,
  ~{
    df <- as.data.frame(.x)
    df$term <- rownames(df)
    as_tibble(df)
  },
  .id = "distance"
) %>%
  filter(term == "natman") %>%
  dplyr::select(distance, R2, `Pr(>F)`)

print(natman_R2)
readr::write_csv2(natman_R2, f_natman_r2_csv)

cat("\nPaired PERMANOVA (Aitchison) — aspect within snag:\n")
print(
  adonis2(
    D_robAit ~ log_reads + umi + ds_at_drill + aspect,
    data = SNAG_meta,
    by = "margin",
    permutations = how(nperm = nperm),
    strata = SNAG_meta$natman
  )
)

cat("\nPaired PERMANOVA (Jaccard) — aspect within snag:\n")
print(
  adonis2(
    D_jacc ~ log_reads + umi + ds_at_drill + aspect,
    data = SNAG_meta,
    by = "margin",
    permutations = how(nperm = nperm),
    strata = SNAG_meta$natman
  )
)

# ---- LOO PERMANOVA -----------------------------------------------------------
loo_blocked_Ait <- loo_permanova(
  X = otu_snag,
  meta = SNAG_meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + aspect + ds_at_drill",
  dist_method = "robust.aitchison",
  nperm_outer = NPERM_LOO,
  by_TYPE = "margin",
  blocked_by = "natman",
  keep_terms = c("log_reads", "umi", "aspect", "ds_at_drill"),
  out_csv = file.path(outdir, "loo_permanova_blocked_SNAG_Ait.csv")
)

loo_blocked_Jac <- loo_permanova(
  X = otu_snag,
  meta = SNAG_meta,
  group_id = "natman",
  formula_rhs = "log_reads + umi + aspect + ds_at_drill",
  dist_method = "jaccard",
  nperm_outer = NPERM_LOO,
  by_TYPE = "margin",
  blocked_by = "natman",
  keep_terms = c("log_reads", "umi", "aspect", "ds_at_drill"),
  out_csv = file.path(outdir, "loo_permanova_blocked_SNAG_Jac.csv")
)

# ==============================================================================
# 7) PERMDISP
# ==============================================================================
bd_as_Ait <- betadisper(D_robAit, SNAG_meta$aspect)
anova(bd_as_Ait); permutest(bd_as_Ait, permutations = nperm); TukeyHSD(bd_as_Ait)
boxplot(bd_as_Ait, xlab = "aspect", main = "Aitchison: distance to group centroid")

bd_ds_Ait <- betadisper(D_robAit, SNAG_meta$ds_at_drill)
anova(bd_ds_Ait); permutest(bd_ds_Ait, permutations = nperm); TukeyHSD(bd_ds_Ait)
boxplot(bd_ds_Ait, xlab = "ds_at_drill", main = "Aitchison: distance to group centroid")

bd_as_Jac <- betadisper(D_jacc, SNAG_meta$aspect)
anova(bd_as_Jac); permutest(bd_as_Jac, permutations = nperm); TukeyHSD(bd_as_Jac)
boxplot(bd_as_Jac, xlab = "aspect", main = "Jaccard: distance to group centroid")

bd_ds_Jac <- betadisper(D_jacc, SNAG_meta$ds_at_drill)
anova(bd_ds_Jac); permutest(bd_ds_Jac, permutations = nperm); TukeyHSD(bd_ds_Jac)
boxplot(bd_ds_Jac, xlab = "ds_at_drill", main = "Jaccard: distance to group centroid")

# ==============================================================================
# 8) Design-preserving natman nulls + ECDF/AUC/ANOSIM
# ==============================================================================
paired_ids <- names(which(table(SNAG_meta$natman) == 2))
SNAG_pair  <- SNAG_meta %>% filter(natman %in% paired_ids)

DmA <- as.matrix(vegdist(otu_snag[SNAG_pair$sample, , drop = FALSE], method = "robust.aitchison"))
U <- which(upper.tri(DmA), arr.ind = TRUE)
DmJ <- as.matrix(vegdist((otu_snag[SNAG_pair$sample, ] > 0) * 1, method = "jaccard", binary = TRUE))

pairs_A <- mk_within_between_tbl(DmA, SNAG_pair$natman)
pairs_J <- mk_within_between_tbl(DmJ, SNAG_pair$natman)

sum_A <- pairs_A %>% group_by(type) %>% summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop")
sum_J <- pairs_J %>% group_by(type) %>% summarise(n = n(), median = median(dist), IQR = IQR(dist), .groups = "drop")

bd_A <- boot_diff_median(pairs_A)
ci_diff_A <- quantile(bd_A, c(0.025, 0.975))
bd_J <- boot_diff_median(pairs_J)
ci_diff_J <- quantile(bd_J, c(0.025, 0.975))

est_A <- with(sum_A, median[type == "between"] - median[type == "within"])
est_J <- with(sum_J, median[type == "between"] - median[type == "within"])

label_A <- sprintf("Δmedian\n%.3f\n95%%CI\n[%.3f, %.3f]\n", est_A, ci_diff_A[1], ci_diff_A[2])
label_J <- sprintf("Δmedian\n%.3f\n95%%CI\n[%.3f, %.3f]\n", est_J, ci_diff_J[1], ci_diff_J[2])

med_w_A <- unname(sum_A$median[sum_A$type == "within"])
med_b_A <- unname(sum_A$median[sum_A$type == "between"])
med_w_J <- unname(sum_J$median[sum_J$type == "within"])
med_b_J <- unname(sum_J$median[sum_J$type == "between"])

xpos_A <- quantile(pairs_A$dist, 0.65)
ypos_A <- 0.2
xpos_J <- quantile(pairs_J$dist, 0.65)
ypos_J <- 0.2

p_ecdf_A <- ggplot(pairs_A, aes(dist, color = type)) +
  stat_ecdf() +
  theme_bw() +
  labs(
    title = "SNAG (Aitchison): ECDF of pairwise distances",
    subtitle = "Dashed lines: medians (within/between)",
    x = "Aitchison distance",
    y = "ECDF"
  )

print(p_ecdf_A)

p_ecdf_J <- ggplot(pairs_J, aes(dist, color = type)) +
  stat_ecdf() +
  theme_bw() +
  labs(
    title = "SNAG (Jaccard): ECDF of pairwise distances",
    subtitle = "Dashed lines: medians (within/between)",
    x = "Jaccard distance",
    y = "ECDF"
  )

print(p_ecdf_J)

ggsave(file.path(outdir, "ECDF_within_between_SNAG_Aitchison.png"), p_ecdf_A, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "ECDF_within_between_SNAG_Jaccard.png"), p_ecdf_J, width = 6, height = 4, dpi = 300)

n_samp <- nrow(SNAG_pair)
stopifnot(n_samp %% 2 == 0)
cat("number of paired SNAG samples:\n", n_samp)

shuffle_pairs <- function(n) factor(sample(rep(seq_len(n / 2), each = 2)))

AUC_obs_Ait <- wilcoxon_auc_once(DmA, SNAG_pair$natman, U)
AUC_obs_Jac <- wilcoxon_auc_once(DmJ, SNAG_pair$natman, U)

set.seed(SEED_GLOBAL)
AUC_null_Ait <- replicate(B_BOOT, wilcoxon_auc_once(DmA, shuffle_pairs(n_samp), U))
AUC_null_Jac <- replicate(B_BOOT, wilcoxon_auc_once(DmJ, shuffle_pairs(n_samp), U))

ciAUC_Ait <- quantile(AUC_null_Ait, c(0.025, 0.975))
pAUC_Ait <- mean(AUC_null_Ait >= AUC_obs_Ait)
ciAUC_Jac <- quantile(AUC_null_Jac, c(0.025, 0.975))
pAUC_Jac <- mean(AUC_null_Jac >= AUC_obs_Jac)

ses_AUC_Ait <- (AUC_obs_Ait - mean(AUC_null_Ait)) / sd(AUC_null_Ait)
ses_AUC_Jac <- (AUC_obs_Jac - mean(AUC_null_Jac)) / sd(AUC_null_Jac)

rb_Ait <- 2 * AUC_obs_Ait - 1
rb_Ait_CI <- 2 * ciAUC_Ait - 1
rb_Jac <- 2 * AUC_obs_Jac - 1
rb_Jac_CI <- 2 * ciAUC_Jac - 1

p_auc_ait <- ggplot(tibble(AUC = AUC_null_Ait), aes(AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_Ait, linetype = 3) +
  labs(
    title = "SNAG AUC null vs observed (Aitchison)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, p(right-tailed)=%.3f, SES=%.2f",
      AUC_obs_Ait, mean(AUC_null_Ait), pAUC_Ait, ses_AUC_Ait
    )
  ) +
  theme_bw()

print(p_auc_ait)

p_auc_jac <- ggplot(tibble(AUC = AUC_null_Jac), aes(AUC)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = AUC_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ciAUC_Jac, linetype = 3) +
  labs(
    title = "SNAG AUC null vs observed (Jaccard)",
    subtitle = sprintf(
      "obs=%.3f, null=%.3f, p(right-tailed)=%.3f, SES=%.2f",
      AUC_obs_Jac, mean(AUC_null_Jac), pAUC_Jac, ses_AUC_Jac
    )
  ) +
  theme_bw()

print(p_auc_jac)

ggsave(file.path(outdir, "AUC_null_vs_observed_SNAG_Aitchison.png"), p_auc_ait, width = 6, height = 4, dpi = 300)
ggsave(file.path(outdir, "AUC_null_vs_observed_SNAG_Jaccard.png"), p_auc_jac, width = 6, height = 4, dpi = 300)

# ---- R² nulls ----------------------------------------------------------------
D_pair_Ait <- as.dist(DmA)
attr(D_pair_Ait, "Labels") <- rownames(DmA)
D_pair_Jac <- as.dist(DmJ)
attr(D_pair_Jac, "Labels") <- rownames(DmJ)

meta_pair <- SNAG_pair

set.seed(SEED_GLOBAL)
ado_obs_Ait <- adonis2(
  D_pair_Ait ~ log_reads + umi + natman + ds_at_drill + aspect,
  data = meta_pair,
  by = "terms",
  permutations = NPERM_MAIN
)
R2_obs_Ait <- as.data.frame(ado_obs_Ait)["natman", "R2"]

R2_null_Ait <- replicate(B_NULL_R2, {
  grp <- shuffle_pairs(nrow(meta_pair))
  meta_tmp <- mutate(meta_pair, grp = grp)
  as.data.frame(
    adonis2(
      D_pair_Ait ~ log_reads + umi + grp + ds_at_drill + aspect,
      data = meta_tmp,
      by = "terms",
      permutations = PERM_TERM
    )
  )["grp", "R2"]
})

saveRDS(R2_null_Ait, file.path(outdir, "snag_R2_null_Ait.RDS"))

ci_Ait_R2   <- quantile(R2_null_Ait, c(0.025, 0.975))
pemp_Ait_R2 <- mean(R2_null_Ait >= R2_obs_Ait)
ses_Ait_R2  <- (R2_obs_Ait - mean(R2_null_Ait)) / sd(R2_null_Ait)
perc_Ait_R2 <- mean(R2_null_Ait <= R2_obs_Ait)

pR2_Ait <- ggplot(tibble(R2 = R2_null_Ait), aes(R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_Ait, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ci_Ait_R2, linetype = 3) +
  labs(
    title = "SNAG natman R² null vs observed (Aitchison)",
    subtitle = sprintf(
      "obs=%.3f; null=%.3f; 95%%CI=[%.3f,%.3f]; p=%.3f; SES=%.2f; perc=%.3f",
      R2_obs_Ait, mean(R2_null_Ait), ci_Ait_R2[1], ci_Ait_R2[2], pemp_Ait_R2, ses_Ait_R2, perc_Ait_R2
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

ggsave(file.path(outdir, "SNAG_natman_R2_null_vs_observed_Aitchison.png"), pR2_Ait, width = 6, height = 4, dpi = 300)
print(pR2_Ait)

ado_obs_Jac <- adonis2(
  D_pair_Jac ~ umi + natman + ds_at_drill + aspect,
  data = meta_pair,
  by = "terms",
  permutations = NPERM_MAIN
)
R2_obs_Jac <- as.data.frame(ado_obs_Jac)["natman", "R2"]

R2_null_Jac <- replicate(B_NULL_R2, {
  grp <- shuffle_pairs(nrow(meta_pair))
  meta_tmp <- mutate(meta_pair, grp = grp)
  as.data.frame(
    adonis2(
      D_pair_Jac ~ umi + grp + ds_at_drill + aspect,
      data = meta_tmp,
      by = "terms",
      permutations = PERM_TERM
    )
  )["grp", "R2"]
})

saveRDS(R2_null_Jac, file.path(outdir, "snag_R2_null_Jac.RDS"))

ci_Jac_R2   <- quantile(R2_null_Jac, c(0.025, 0.975))
pemp_Jac_R2 <- mean(R2_null_Jac >= R2_obs_Jac)
ses_Jac_R2  <- (R2_obs_Jac - mean(R2_null_Jac)) / sd(R2_null_Jac)
perc_Jac_R2 <- mean(R2_null_Jac <= R2_obs_Jac)

pR2_Jac <- ggplot(tibble(R2 = R2_null_Jac), aes(R2)) +
  geom_density(fill = "grey80") +
  geom_vline(xintercept = R2_obs_Jac, linetype = 2, colour = "blue") +
  geom_vline(xintercept = ci_Jac_R2, linetype = 3) +
  labs(
    title = "SNAG natman R² null vs observed (Jaccard)",
    subtitle = sprintf(
      "obs=%.3f; null=%.3f; 95%%CI=[%.3f,%.3f]; p=%.3f; SES=%.2f; perc=%.3f",
      R2_obs_Jac, mean(R2_null_Jac), ci_Jac_R2[1], ci_Jac_R2[2], pemp_Jac_R2, ses_Jac_R2, perc_Jac_R2
    ),
    x = expression(R^2),
    y = "Density"
  ) +
  theme_bw()

ggsave(file.path(outdir, "SNAG_natman_R2_null_vs_observed_Jaccard.png"), pR2_Jac, width = 6, height = 4, dpi = 300)
print(pR2_Jac)

natman_R2_nulls <- tibble(
  metric          = c("robust.aitchison", "jaccard"),
  R2_obs_natman   = c(R2_obs_Ait, R2_obs_Jac),
  R2_null_mean    = c(mean(R2_null_Ait), mean(R2_null_Jac)),
  R2_null_CI_lo   = c(ci_Ait_R2[1], ci_Jac_R2[1]),
  R2_null_CI_hi   = c(ci_Ait_R2[2], ci_Jac_R2[2]),
  p_R2_emp        = c(pemp_Ait_R2, pemp_Jac_R2),
  R2_SES          = c(ses_Ait_R2, ses_Jac_R2),
  R2_percentile   = c(perc_Ait_R2, perc_Jac_R2)
) %>%
  mutate(p_tail = "right", alternative = "within < between")

print(natman_R2_nulls)

# ---- ANOSIM ------------------------------------------------------------------
set.seed(SEED_GLOBAL)
anos_Ait <- anosim(D_pair_Ait, grouping = meta_pair$natman, permutations = ANOSIM_PERM)
anos_Jac <- anosim(D_pair_Jac, grouping = meta_pair$natman, permutations = ANOSIM_PERM)

summary_both <- tibble(
  state = "SNAG",
  metric = c("robust.aitchison", "jaccard"),
  dist_within_median  = c(sum_A$median[sum_A$type == "within"], sum_J$median[sum_J$type == "within"]),
  dist_between_median = c(sum_A$median[sum_A$type == "between"], sum_J$median[sum_J$type == "between"]),
  dist_diff           = c(est_A, est_J),
  diff_CI_lo          = c(ci_diff_A[1], ci_diff_J[1]),
  diff_CI_hi          = c(ci_diff_A[2], ci_diff_J[2]),
  AUC_obs             = c(AUC_obs_Ait, AUC_obs_Jac),
  AUC_null_mean       = c(mean(AUC_null_Ait), mean(AUC_null_Jac)),
  AUC_null_CI_lo      = c(ciAUC_Ait[1], ciAUC_Jac[1]),
  AUC_null_CI_hi      = c(ciAUC_Ait[2], ciAUC_Jac[2]),
  AUC_p_emp           = c(pAUC_Ait, pAUC_Jac),
  AUC_SES             = c(ses_AUC_Ait, ses_AUC_Jac),
  ANOSIM_R            = c(anos_Ait$statistic, anos_Jac$statistic),
  ANOSIM_p            = c(anos_Ait$signif, anos_Jac$signif)
) %>%
  left_join(natman_R2_nulls, by = "metric")

summary_both %>%
  pivot_longer(
    cols = -c(state, metric),
    names_to = "statistic",
    values_to = "value",
    values_transform = list(value = as.character)
  ) %>%
  arrange(metric, statistic) %>%
  print(n = Inf)

readr::write_csv2(summary_both, f_summary_csv)

# ==============================================================================
# 10) betapart
# ==============================================================================
PA_snag <- (otu_snag > 0) * 1
bp <- betapart::beta.pair(PA_snag, index.family = "sorensen")

DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

labs <- rownames(PA_snag)
nm   <- SNAG_meta$natman[match(labs, SNAG_meta$sample)]
asp  <- SNAG_meta$aspect[match(labs, SNAG_meta$sample)]

pairs_df <- which(upper.tri(DM_sim), arr.ind = TRUE) %>%
  as_tibble() %>%
  transmute(
    i = labs[row],
    j = labs[col],
    nat_i = nm[row],
    nat_j = nm[col],
    asp_i = asp[row],
    asp_j = asp[col],
    type = if_else(nat_i == nat_j, "within_snag", "between_snag"),
    beta_SIM = DM_sim[cbind(row, col)],
    beta_SNE = DM_sne[cbind(row, col)],
    beta_SOR = DM_sor[cbind(row, col)]
  ) %>%
  filter(!(type == "within_snag" & asp_i == asp_j))

write.csv(pairs_df, f_betapart_pairs, row.names = FALSE)

summary_long <- pairs_df %>%
  pivot_longer(c(beta_SIM, beta_SNE, beta_SOR), names_to = "metric", values_to = "value") %>%
  group_by(metric, type) %>%
  summarise(n_pairs = n(), mean = mean(value), median = median(value), IQR = IQR(value), .groups = "drop") %>%
  mutate(
    metric_label = dplyr::case_when(
      metric == "beta_SIM" ~ "βSIM (turnover)",
      metric == "beta_SNE" ~ "βSNE (nestedness)",
      metric == "beta_SOR" ~ "βSOR (total)",
      TRUE ~ metric
    )
  )

tn_ratio <- summary_long %>%
  dplyr::select(metric, type, median) %>%
  pivot_wider(names_from = metric, values_from = median) %>%
  mutate(turnover_to_nestedness = `beta_SIM` / pmax(`beta_SNE`, .Machine$double.eps)) %>%
  dplyr::select(type, turnover_to_nestedness)

tn_ratio_wide <- tn_ratio %>%
  pivot_wider(names_from = type, values_from = turnover_to_nestedness, names_glue = "turnover_to_nestedness_{type}")

effect_tab <- function(metric_col) {
  within_vals  <- pairs_df %>% pull({{ metric_col }}) %>% `[`(pairs_df$type == "within_snag")
  between_vals <- pairs_df %>% pull({{ metric_col }}) %>% `[`(pairs_df$type == "between_snag")
  
  bind_cols(
    boot_diff_ci(between_vals, within_vals),
    tibble(cliffs_delta = cliffs_delta(between_vals, within_vals))
  )
}

eff_SIM <- effect_tab(beta_SIM) %>% mutate(metric = "βSIM (turnover)")
eff_SNE <- effect_tab(beta_SNE) %>% mutate(metric = "βSNE (nestedness)")
eff_SOR <- effect_tab(beta_SOR) %>% mutate(metric = "βSOR (total)")
effects_all <- bind_rows(eff_SIM, eff_SNE, eff_SOR)

summary_wide <- summary_long %>%
  dplyr::select(metric_label, type, n_pairs, mean, median, IQR) %>%
  pivot_wider(names_from = type, values_from = c(n_pairs, mean, median, IQR), names_glue = "{.value}_{type}") %>%
  mutate(dummy = 1) %>%
  left_join(tn_ratio_wide %>% mutate(dummy = 1), by = "dummy") %>%
  dplyr::select(-dummy) %>%
  rename(metric = metric_label) %>%
  left_join(effects_all, by = "metric") %>%
  relocate(metric, est_diff_median, CI_low, CI_high, cliffs_delta)

write_csv2(summary_wide, f_betapart_sum)
print(summary_wide)

# ==============================================================================
# 11) Alpha diversity models — q0
# ==============================================================================
cat("\n[ ALPHA RICHNESS GLMM - SNAG ]\n")

df_snag <- alpha_from_otu(otu_snag) %>%
  dplyr::left_join(
    SNAG_meta %>% dplyr::select(sample, natman, aspect, ds_at_drill, umi),
    by = "sample"
  ) %>%
  dplyr::mutate(
    natman      = droplevels(natman),
    aspect      = droplevels(aspect),
    ds_at_drill = droplevels(ds_at_drill),
    umi         = umi
  ) %>%
  dplyr::filter(!is.na(richness), reads > 0)

rhs_terms_glmm_snag <- c("log_reads", "umi", "ds_at_drill", "aspect")
form_snag <- as.formula(paste("richness ~", paste(rhs_terms_glmm_snag, collapse = " + ")))
cat("GLMM form (SNAG):\n")
print(form_snag)

if (nlevels(df_snag$aspect) >= 2 && nlevels(df_snag$ds_at_drill) >= 2) {
  m_pois_snag <- suppressWarnings(
    stats::glm(form_snag, data = df_snag, family = poisson(link = "log"))
  )
  
  phi_p_snag  <- overdisp_phi(m_pois_snag)
  use_nb_snag <- is.finite(phi_p_snag$phi) && (phi_p_snag$phi > 1.2)
  
  if (use_nb_snag) {
    m_nb_snag <- suppressWarnings(MASS::glm.nb(form_snag, data = df_snag))
    best_snag <- if (AIC(m_nb_snag) + 2 < AIC(m_pois_snag)) m_nb_snag else m_pois_snag
    best_type_snag <- if (identical(best_snag, m_nb_snag)) "NB (glm.nb)" else "Poisson (glm)"
  } else {
    best_snag <- m_pois_snag
    best_type_snag <- "Poisson (glm)"
  }
  
  cat(sprintf(
    "Selected SNAG alpha model: %s  Overdispersion Poisson phi = %.2f  AIC(best)=%.1f\n",
    best_type_snag, phi_p_snag$phi, AIC(best_snag)
  ))
  print(summary(best_snag))
  
  an_tab_snag <- car::Anova(best_snag, type = 2)
  print(an_tab_snag)
  r2_snag <- performance::r2(best_snag)
  print(r2_snag)
  
  cat("adjusted pseudo R2 that is better to interpret \n 1 - (best_snag$deviance / best_snag$null.deviance")
  print(1 - (best_snag$deviance / best_snag$null.deviance))
  
  irr_snag <- broom::tidy(best_snag, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(
      term = dplyr::case_when(
        term == "(Intercept)" ~ "Intercept",
        term == "log_reads" ~ "log(reads+1)",
        grepl("^umi", term) ~ paste0("Plate: ", sub("^umi", "", term), " vs ", levels(df_snag$umi)[1]),
        grepl("^aspect", term) ~ paste0("Aspect: ", sub("^aspect", "", term), " vs ", levels(df_snag$aspect)[1]),
        grepl("^ds_at_drill", term) ~ paste0("Decay ", sub("^ds_at_drill", "", term), " vs ", levels(df_snag$ds_at_drill)[1]),
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
  
  emm_aspect <- emmeans::emmeans(best_snag, ~ aspect, type = "response")
  df_emm_as  <- as.data.frame(summary(emm_aspect))
  
  p_emm_as <- ggplot2::ggplot(
    df_emm_as,
    ggplot2::aes(x = aspect, y = response, ymin = asymp.LCL, ymax = asymp.UCL)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(width = 0.2) +
    ggplot2::labs(x = "Aspect", y = "Estimated richness", title = "Richness ~ aspect (SNAG)") +
    ggplot2::theme_classic(base_size = 11)
  
  ggplot2::ggsave(file.path(outdir, "alpha_richness_emm_aspect_SNAG.png"), p_emm_as, width = 6, height = 3.2, dpi = 300)
  print(p_emm_as)
  
  emm_ds_snag <- emmeans::emmeans(best_snag, ~ ds_at_drill, type = "response")
  df_emm_ds   <- as.data.frame(summary(emm_ds_snag))
  
  p_emm_ds_snag <- ggplot2::ggplot(
    df_emm_ds,
    ggplot2::aes(x = ds_at_drill, y = response, ymin = asymp.LCL, ymax = asymp.UCL)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(width = 0.2) +
    ggplot2::labs(x = "Decay at drill (ordered 1→5)", y = "Estimated richness", title = "Richness ~ decay (SNAG)") +
    ggplot2::theme_classic(base_size = 11)
  
  ggplot2::ggsave(file.path(outdir, "alpha_richness_emm_DS_SNAG.png"), p_emm_ds_snag, width = 6.2, height = 3.2, dpi = 300)
  print(p_emm_ds_snag)
  
  keep_eff_snag <- c("aspect", "ds_at_drill", "log_reads", "umi")
  tbl_main_snag <- tibble::tibble(
    Effect  = rownames(an_tab_snag),
    Chisq   = an_tab_snag$`Chisq`,
    df      = an_tab_snag$Df,
    p_value = an_tab_snag$`Pr(>Chisq)`
  ) %>%
    dplyr::filter(Effect %in% keep_eff_snag) %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    dplyr::arrange(match(Effect, keep_eff_snag))
  
  r2_m <- if ("R2_marginal" %in% names(r2_snag)) r2_snag$R2_marginal else NA_real_
  r2_c <- if ("R2_conditional" %in% names(r2_snag)) r2_snag$R2_conditional else NA_real_
  tbl_r2_snag <- tibble::tibble(
    Metric = c("R2_marginal", "R2_conditional"),
    Value  = c(r2_m, r2_c)
  )
  
  cat("\n=== SNAG tables ===\n")
  print(list(anova = tbl_main_snag, r2 = tbl_r2_snag, irr = irr_snag))
  
  capture.output(
    list(model = summary(best_snag), anova = an_tab_snag, r2 = r2_snag),
    file = file.path(outdir, "alpha_glm_SNAG.txt")
  )
  
  emmeans::emmeans(best_snag, ~ aspect, type = "response")
  emmeans::emmeans(best_snag, ~ ds_at_drill, type = "response")
  emmeans::emmeans(best_snag, ~ umi, type = "response")
}

# ==============================================================================
# 12) RDS bundle
# ==============================================================================
snag_results <- list(
  meta = SNAG_meta,
  otu  = otu_snag,
  distances = list(Ait = D_robAit, Jac = D_jacc),
  permanova = perma_unblocked,
  paired = list(
    Ait = "adonis2(D_robAit ~ log_reads + umi + aspect, strata=natman)",
    Jac = "adonis2(D_jacc ~ umi + aspect, strata=natman)"
  ),
  betapart = list(summary = summary_wide, effects = effects_all),
  alpha = list(
    q0 = list(
      model      = if (exists("best_snag")) best_snag else NULL,
      model_type = if (exists("best_type_snag")) best_type_snag else NULL,
      anova      = if (exists("an_tab_snag")) an_tab_snag else NULL,
      r2         = if (exists("r2_snag")) r2_snag else NULL,
      emm_aspect = if (exists("emm_aspect")) emm_aspect else NULL,
      emm_ds     = if (exists("emm_ds_snag")) emm_ds_snag else NULL
    )
  ),
  within_between = list(Ait = sum_A, Jac = sum_J),
  natman_design_nulls = list(
    AUC = list(
      obs = c(Ait = AUC_obs_Ait, Jac = AUC_obs_Jac),
      null = list(Ait = AUC_null_Ait, Jac = AUC_null_Jac),
      ci = list(Ait = ciAUC_Ait, Jac = ciAUC_Jac),
      rb = list(Ait = rb_Ait, Jac = rb_Jac)
    ),
    R2 = natman_R2_nulls
  )
)

saveRDS(snag_results, f_rds)
