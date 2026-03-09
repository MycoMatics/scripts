# =============================================================================
# Script: Script4_scope2A_LOG.R
# Project: Deadwood fungal community analyses
# Scope: LOG
# Purpose: Within-substrate analysis for logs, including PERMANOVA, null models,
#          ECDF/AUC/ANOSIM summaries, leave-one-tree-out stability, dispersion,
#          alpha-diversity GLMMs, and betapart decomposition
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
# - tax
# - threshold
#
# Required helper functions from Script0_utils.R:
# - alpha_from_otu
# - overdisp_phi
# - check_glmm_lme4
# - plot_resid_glmm
# - get_natman_R2
# - auc_from_group
# - pairwise_groups
# - mk_within_between_tbl
# - boot_diff_median
# - build_intra_meta_ds_at_drill
# =============================================================================
source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(janitor)
  library(permute)
  library(emmeans)
  library(car)
  library(broom.mixed)
  library(performance)
  library(indicspecies)
  library(eulerr)
  library(forcats)
  library(grid)
  library(patchwork)
  library(betapart)
  library(DescTools)
})

assert_objects(c("META1", "otu_matrix_filt", "tax", "threshold"))
assert_objects("build_intra_meta_ds_at_drill")

set.seed(SEED_GLOBAL)

# ---- Paths and scope ---------------------------------------------------------
data_dir <- "."
scope <- "LOG"

ensure_dir(data_dir)
setwd(data_dir)

base_outdir <- file.path("plots", "INTRA_diam")
outdir <- file.path(base_outdir, scope)
out_nat <- file.path(outdir, "TREE_IDENTITY")

ensure_dir(base_outdir)
ensure_dir(outdir)
ensure_dir(out_nat)

cat("You are using SH threshold ", threshold, "\n", sep = "")
cat(" INTRA-LOG RUNNER — scope:", scope, "\n")
cat("================================================================================\n")

# ---- Analysis parameters -----------------------------------------------------
N_PERM_MAIN   <- 999
N_PERM_ANOSIM <- 9999
N_PERM_NULL   <- 199
B_NULL        <- 999
MIN_SAMPLES   <- 3
N_PERM_LOO    <- 499

# ---- 1) Build metadata -------------------------------------------------------
meta <- build_intra_meta_ds_at_drill(scope)
print(meta)
print(table(meta$natman))

  # ---- 2) Align OTU table ----------------------------------------------------
  otu <- otu_matrix_filt[meta$sample, , drop = FALSE]
  stopifnot(identical(rownames(otu), meta$sample))
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  
  cat("Number of SH in analysis: ", ncol(otu), "\n")
  print(meta, n = Inf)
  cat("NATMAN UNIQUE:\n")
  print(unique(meta$natman))
  
  # ---- 3) Add read depth and repair factors ---------------------------------
  alpha_df <- alpha_from_otu(otu)
  
  meta_ext <- meta %>%
    left_join(alpha_df %>% select(sample, reads, log_reads), by = "sample") %>%
    mutate(
      reads = ifelse(is.na(reads), 0, reads),
      log_reads = log(reads),
      umi = droplevels(factor(umi)),
      natman = droplevels(factor(natman)),
      Position = droplevels(factor(Position)),
      size = droplevels(factor(size)),
      diameter_at_drill_z = diameter_at_drill_z,
      ds_at_drill = droplevels(factor(ds_at_drill)),
      depth_2 = if ("depth_2" %in% names(.)) {
        factor(depth_2, levels = c("INNER", "OUTER"))
      } else {
        factor(NA)
      }
    ) %>%
    filter(!is.na(diameter_at_drill_z))
  
  # ---- 4) Collinearity checks ------------------------------------------------
  cat("\n checking VIF \n")
  
  mm <- lm(
    log_reads ~ umi + Position + ds_at_drill + diameter_at_drill_z,
    data = meta_ext
  )
  
  print(summary(mm)$aliased)
  print(car::vif(mm))
  
  dat <- meta_ext %>%
    dplyr::select(
      natman,
      Position,
      umi,
      ds_at_drill,
      diameter_at_drill_z,
      log_reads
    ) %>%
    dplyr::mutate(
      diameter_at_drill_z = as.numeric(diameter_at_drill_z),
      log_reads = as.numeric(log_reads)
    )
  
  cont_vars <- c("diameter_at_drill_z", "log_reads")
  fact_vars <- c("umi", "natman", "Position", "ds_at_drill")
  
  cor_cc <- dat %>%
    dplyr::select(all_of(cont_vars)) %>%
    stats::cor(use = "pairwise.complete.obs", method = "pearson")
  
  print("Continuous–continuous Pearson correlations")
  print(cor_cc)
  
  cramer_results <- list()
  for (i in 1:(length(fact_vars) - 1)) {
    for (j in (i + 1):length(fact_vars)) {
      tbl <- table(dat[[fact_vars[i]]], dat[[fact_vars[j]]])
      cv <- DescTools::CramerV(tbl, bias.correct = TRUE)
      cramer_results[[paste(fact_vars[i], fact_vars[j], sep = "_vs_")]] <- cv
    }
  }
  
  print("Factor–factor associations (Cramér's V)")
  print(cramer_results)
  print("tree identity (natman) structures some covars but that is to be expected")
  
  model_vars <- c(
    "log_reads", "umi", "natman", "ds_at_drill",
    "depth_2", "Position", "diameter_at_drill_z"
  )
  print(model_vars)
  
  otu <- otu[meta_ext$sample, , drop = FALSE]
  rhs_terms <- model_vars
  
  mlevels <- levels(meta_ext$Position)
  micro_cols <- microhab_base_cols[names(microhab_base_cols) %in% mlevels]
  micro_cols <- micro_cols[match(mlevels, names(microhab_base_cols))]
  depth2_cols <- depth_colors
  
  # ---- Distances -------------------------------------------------------------
  D_robAit <- vegdist(otu, method = "robust.aitchison")
  D_robAit <- fix_dist_labels(D_robAit, samps=rownames(otu))
  otu_pa <- (otu > 0) * 1L
  D_jaccard <- vegdist(otu_pa, method = "jaccard", binary = TRUE)
  
  # ---- PERMANOVA -------------------------------------------------------------
  cat("\n[ PERMANOVAs ]\n")
  cat("\n RHS terms are: \n", rhs_terms)
  
  perm_unblocked <- permute::how(nperm = N_PERM_MAIN)
  
  perma_unblocked <- list(
    robAit = adonis2(
      D_robAit ~ .,
      data = meta_ext[, rhs_terms, drop = FALSE],
      by = "terms",
      permutations = perm_unblocked
    ),
    jaccard = adonis2(
      D_jaccard ~ .,
      data = meta_ext[, rhs_terms, drop = FALSE],
      by = "terms",
      permutations = perm_unblocked
    )
  )
  
  capture.output(
    perma_unblocked,
    file = file.path(outdir, paste0("permanova_unblocked_", scope, ".txt"))
  )
  cat("\npermanova UN_blocked:\n")
  print(perma_unblocked)
  
  natman_R2 <- purrr::map_dfr(
    perma_unblocked,
    ~{
      df <- as.data.frame(.x)
      df$term <- rownames(df)
      tibble::as_tibble(df)
    },
    .id = "distance"
  ) %>%
    filter(term == "natman") %>%
    select(distance, R2, `Pr(>F)`)
  
  cat("\nUnblocked PERMANOVA — natman R2 by distance:\n")
  print(natman_R2)
  
  perm_blocked <- permute::how(nperm = N_PERM_MAIN)
  permute::setBlocks(perm_blocked) <- meta_ext$natman
  
  rhs_block <- setdiff(rhs_terms, "natman")
  
  perma_blocked <- list(
    robAit = adonis2(
      D_robAit ~ .,
      data = meta_ext[, rhs_block, drop = FALSE],
      by = "margin",
      permutations = perm_blocked
    ),
    Jacc = adonis2(
      D_jaccard ~ .,
      data = meta_ext[, rhs_block, drop = FALSE],
      by = "margin",
      permutations = perm_blocked
    )
  )
  
  capture.output(
    perma_blocked,
    file = file.path(outdir, paste0("permanova_blocked_", scope, ".txt"))
  )
  cat("\npermanova BLOCKED (within tree):\n")
  print(perma_blocked)
  
  # ---- Tree identity effect sizes and nulls ----------------------------------
  cat("\n[ TREE EFFECT — effect sizes + label shuffle nulls ]\n")
  
  upper_idx <- which(upper.tri(as.matrix(D_robAit)), arr.ind = TRUE)
  nat <- droplevels(meta_ext$natman)
  
  R2_obs_Ait <- get_natman_R2(perma_unblocked$robAit)
  R2_obs_Jac <- get_natman_R2(perma_unblocked$jaccard)
  
  AUC_obs_A <- auc_from_group(D_robAit, nat, upper_idx = upper_idx)
  AUC_obs_J <- auc_from_group(D_jaccard, nat, upper_idx = upper_idx)
  
  anos_A <- vegan::anosim(D_robAit, grouping = nat, permutations = N_PERM_ANOSIM)
  anos_J <- vegan::anosim(D_jaccard, grouping = nat, permutations = N_PERM_ANOSIM)
  print(anos_A)
  print(anos_J)
  
  df_terms <- meta_ext[, rhs_terms, drop = FALSE]
  R2_null_Ait <- numeric(B_NULL)
  R2_null_Jac <- numeric(B_NULL)
  AUC_null_A <- numeric(B_NULL)
  AUC_null_J <- numeric(B_NULL)
  
  for (b in seq_len(B_NULL)) {
    nat_null <- factor(sample(nat, replace = FALSE), levels = levels(nat))
    
    df_tmp <- df_terms
    if ("natman" %in% names(df_tmp)) df_tmp$natman <- nat_null
    
    R2_null_Ait[b] <- get_natman_R2(
      adonis2(D_robAit ~ ., data = df_tmp, by = "terms", permutations = N_PERM_NULL)
    )
    R2_null_Jac[b] <- get_natman_R2(
      adonis2(D_jaccard ~ ., data = df_tmp, by = "terms", permutations = N_PERM_NULL)
    )
    
    AUC_null_A[b] <- auc_from_group(D_robAit, nat_null, upper_idx = upper_idx)["AUC"]
    AUC_null_J[b] <- auc_from_group(D_jaccard, nat_null, upper_idx = upper_idx)["AUC"]
    
    if (b %% 25 == 0) cat("[nulls] finished", b, "of", B_NULL, "\n")
  }
  
  ciR2_Ait <- quantile(R2_null_Ait, c(0.025, 0.975), na.rm = TRUE)
  ciR2_Jac <- quantile(R2_null_Jac, c(0.025, 0.975), na.rm = TRUE)
  
  pR2_Ait <- mean(R2_null_Ait >= R2_obs_Ait, na.rm = TRUE)
  pR2_Jac <- mean(R2_null_Jac >= R2_obs_Jac, na.rm = TRUE)
  
  sesR2_Ait <- (R2_obs_Ait - mean(R2_null_Ait, na.rm = TRUE)) / sd(R2_null_Ait, na.rm = TRUE)
  sesR2_Jac <- (R2_obs_Jac - mean(R2_null_Jac, na.rm = TRUE)) / sd(R2_null_Jac, na.rm = TRUE)
  
  percR2_Ait <- mean(R2_null_Ait <= R2_obs_Ait, na.rm = TRUE)
  percR2_Jac <- mean(R2_null_Jac <= R2_obs_Jac, na.rm = TRUE)
  
  ciAUC_A <- quantile(AUC_null_A, c(0.025, 0.975), na.rm = TRUE)
  ciAUC_J <- quantile(AUC_null_J, c(0.025, 0.975), na.rm = TRUE)
  
  pAUC_A <- mean(AUC_null_A >= AUC_obs_A["AUC"], na.rm = TRUE)
  pAUC_J <- mean(AUC_null_J >= AUC_obs_J["AUC"], na.rm = TRUE)
  
  sesAUC_A <- (AUC_obs_A["AUC"] - mean(AUC_null_A, na.rm = TRUE)) / sd(AUC_null_A, na.rm = TRUE)
  sesAUC_J <- (AUC_obs_J["AUC"] - mean(AUC_null_J, na.rm = TRUE)) / sd(AUC_null_J, na.rm = TRUE)
  
  pR2_A <- ggplot(tibble(R2 = R2_null_Ait), aes(R2)) +
    geom_density(fill = "grey80") +
    geom_vline(xintercept = R2_obs_Ait, linetype = 2, col = "blue") +
    geom_vline(xintercept = ciR2_Ait, linetype = 3) +
    labs(
      title = paste0(scope, " natman — R² null vs obs (Aitchison)"),
      subtitle = sprintf(
        "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
        R2_obs_Ait, mean(R2_null_Ait, na.rm = TRUE), ciR2_Ait[1], ciR2_Ait[2], pR2_Ait, sesR2_Ait
      ),
      x = expression(R^2),
      y = "Density"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_nat, paste0("natman_R2_null_vs_obs_Aitchison_", scope, ".png")),
    pR2_A, width = 6, height = 4, dpi = 300
  )
  print(pR2_A)
  
  pR2_J <- ggplot(tibble(R2 = R2_null_Jac), aes(R2)) +
    geom_density(fill = "grey80") +
    geom_vline(xintercept = R2_obs_Jac, linetype = 2, col = "blue") +
    geom_vline(xintercept = ciR2_Jac, linetype = 3) +
    labs(
      title = paste0(scope, " natman — R² null vs obs (Jaccard)"),
      subtitle = sprintf(
        "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
        R2_obs_Jac, mean(R2_null_Jac, na.rm = TRUE), ciR2_Jac[1], ciR2_Jac[2], pR2_Jac, sesR2_Jac
      ),
      x = expression(R^2),
      y = "Density"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_nat, paste0("natman_R2_null_vs_obs_Jaccard_", scope, ".png")),
    pR2_J, width = 6, height = 4, dpi = 300
  )
  print(pR2_J)
  
  # ---- ECDF plots ------------------------------------------------------------
  pw_A <- pairwise_groups(D_robAit, nat, upper_idx)
  df_ecdf_A <- tibble(
    dist = c(pw_A$within, pw_A$between),
    type = factor(
      rep(c("within", "between"), c(length(pw_A$within), length(pw_A$between))),
      levels = c("within", "between")
    )
  )
  
  p_ecdf_A <- ggplot(df_ecdf_A, aes(x = dist, colour = type)) +
    stat_ecdf(geom = "step", linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0(scope, " — ECDF of pairwise distances by natman (Aitchison)"),
      subtitle = sprintf(
        "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
        AUC_obs_A["AUC"], AUC_obs_A["r_rb"], length(pw_A$within), length(pw_A$between)
      ),
      x = "Pairwise distance",
      y = "Cumulative probability",
      colour = "Pair type"
    )
  
  ggsave(
    file.path(out_nat, paste0("natman_ECDF_Aitchison_", scope, ".png")),
    p_ecdf_A, width = 6.2, height = 4.2, dpi = 300
  )
  print(p_ecdf_A)
  
  pw_J <- pairwise_groups(D_jaccard, nat, upper_idx)
  df_ecdf_J <- tibble(
    dist = c(pw_J$within, pw_J$between),
    type = factor(
      rep(c("within", "between"), c(length(pw_J$within), length(pw_J$between))),
      levels = c("within", "between")
    )
  )
  
  p_ecdf_J <- ggplot(df_ecdf_J, aes(x = dist, colour = type)) +
    stat_ecdf(geom = "step", linewidth = 1) +
    theme_bw() +
    labs(
      title = paste0(scope, " — ECDF of pairwise distances by natman (Jaccard)"),
      subtitle = sprintf(
        "AUC = %.3f  r_rb = %.2f  n_within = %d  n_between = %d",
        AUC_obs_J["AUC"], AUC_obs_J["r_rb"], length(pw_J$within), length(pw_J$between)
      ),
      x = "Pairwise distance",
      y = "Cumulative probability",
      colour = "Pair type"
    )
  
  ggsave(
    file.path(out_nat, paste0("natman_ECDF_Jaccard_", scope, ".png")),
    p_ecdf_J, width = 6.2, height = 4.2, dpi = 300
  )
  print(p_ecdf_J)
  
  saveRDS(R2_null_Ait, file.path(outdir, paste0("R2_null_Ait_", scope, ".rds")))
  saveRDS(R2_null_Jac, file.path(outdir, paste0("R2_null_Jac_", scope, ".rds")))
  
  pAUC_A_plot <- ggplot(tibble(AUC = AUC_null_A), aes(AUC)) +
    geom_density(fill = "grey80") +
    geom_vline(xintercept = unname(AUC_obs_A["AUC"]), linetype = 2, colour = "blue") +
    geom_vline(xintercept = ciAUC_A, linetype = 3) +
    labs(
      title = paste0(scope, " natman AUC null vs obs (Aitchison)"),
      subtitle = sprintf(
        "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
        AUC_obs_A["AUC"], mean(AUC_null_A, na.rm = TRUE),
        ciAUC_A[1], ciAUC_A[2], pAUC_A, sesAUC_A
      ),
      x = "AUC (within vs between)",
      y = "Density"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_nat, paste0("natman_AUC_null_vs_obs_Aitchison_", scope, ".png")),
    pAUC_A_plot, width = 6, height = 4, dpi = 300
  )
  print(pAUC_A_plot)
  
  pAUC_J_plot <- ggplot(tibble(AUC = AUC_null_J), aes(AUC)) +
    geom_density(fill = "grey80") +
    geom_vline(xintercept = unname(AUC_obs_J["AUC"]), linetype = 2, colour = "blue") +
    geom_vline(xintercept = ciAUC_J, linetype = 3) +
    labs(
      title = paste0(scope, " natman AUC null vs obs (Jaccard)"),
      subtitle = sprintf(
        "obs=%.3f  null mean=%.3f  95%%CI=[%.3f,%.3f]  p=%.3f  SES=%.2f",
        AUC_obs_J["AUC"], mean(AUC_null_J, na.rm = TRUE),
        ciAUC_J[1], ciAUC_J[2], pAUC_J, sesAUC_J
      ),
      x = "AUC (within vs between)",
      y = "Density"
    ) +
    theme_bw()
  
  ggsave(
    file.path(out_nat, paste0("natman_AUC_null_vs_obs_Jaccard_", scope, ".png")),
    pAUC_J_plot, width = 6, height = 4, dpi = 300
  )
  print(pAUC_J_plot)
  
  cat("RDS R² objects in:\n")
  print(outdir)
  
  summary_nat <- tibble(
    scope = scope,
    metric = c("robust.aitchison", "jaccard"),
    R2_obs = c(R2_obs_Ait, R2_obs_Jac),
    R2_null_mean = c(mean(R2_null_Ait, na.rm = TRUE), mean(R2_null_Jac, na.rm = TRUE)),
    R2_null_CI_lo = c(ciR2_Ait[1], ciR2_Jac[1]),
    R2_null_CI_hi = c(ciR2_Ait[2], ciR2_Jac[2]),
    p_R2_emp = c(pR2_Ait, pR2_Jac),
    R2_SES = c(sesR2_Ait, sesR2_Jac),
    R2_percentile = c(percR2_Ait, percR2_Jac),
    AUC_obs = c(unname(AUC_obs_A["AUC"]), unname(AUC_obs_J["AUC"])),
    rank_biserial = c(unname(AUC_obs_A["r_rb"]), unname(AUC_obs_J["r_rb"])),
    AUC_null_mean = c(mean(AUC_null_A, na.rm = TRUE), mean(AUC_null_J, na.rm = TRUE)),
    AUC_null_CI_lo = c(ciAUC_A[1], ciAUC_J[1]),
    AUC_null_CI_hi = c(ciAUC_A[2], ciAUC_J[2]),
    AUC_p_emp = c(pAUC_A, pAUC_J),
    AUC_SES = c(sesAUC_A, sesAUC_J),
    ANOSIM_R = c(anos_A$statistic, anos_J$statistic),
    ANOSIM_p = c(anos_A$signif, anos_J$signif),
    n_natman_levels = length(levels(nat)),
    n_samples = nrow(meta_ext)
  )
  
  print(summary_nat)
  readr::write_csv2(summary_nat, file.path(out_nat, paste0("natman_evidence_", scope, ".csv")))
  
  # ---- LOO stability ---------------------------------------------------------
  cat("\n[ LOO STABILITY - within tree PERMANOVA ]\n")
  
  trees <- meta_ext$natman %>%
    as.character() %>%
    unique() %>%
    setdiff(NA)
  
  loo <- purrr::map_dfr(trees, function(t) {
    keep <- meta_ext$natman != t & !is.na(meta_ext$natman)
    if (sum(keep) < 6) return(NULL)
    
    Dk <- vegan::vegdist(otu[keep, , drop = FALSE], method = "robust.aitchison")
    
    permk <- permute::how(
      blocks = droplevels(factor(meta_ext$natman[keep])),
      nperm = N_PERM_LOO
    )
    
    rhs_loo <- setdiff(rhs_terms, "natman")
    if (length(rhs_loo) < 1) return(NULL)
    
    m <- vegan::adonis2(
      Dk ~ .,
      data = meta_ext[keep, rhs_loo, drop = FALSE],
      by = "margin",
      permutations = permk
    )
    
    mm <- as.data.frame(m)
    mm$term <- rownames(mm)
    
    tibble::as_tibble(mm) %>%
      filter(term %in% intersect(c("ds_at_drill", "depth_2", "Position", "diameter_at_drill_z"), rhs_loo)) %>%
      transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  
  if (nrow(loo) > 0) {
    loo_sum <- loo %>%
      group_by(term) %>%
      summarise(
        median_R2 = median(R2, na.rm = TRUE),
        p05 = quantile(p, 0.05, na.rm = TRUE),
        p95 = quantile(p, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(loo_sum)
    readr::write_csv2(
      loo_sum,
      file.path(outdir, paste0("loo_within_tree_summary_", scope, ".csv"))
    )
  }
  
  cat("\n[ LOO STABILITY - within tree PERMANOVA (Jaccard) ]\n")
  
  loo_jac <- purrr::map_dfr(trees, function(t) {
    keep <- meta_ext$natman != t & !is.na(meta_ext$natman)
    if (sum(keep) < 6) return(NULL)
    
    X_pa <- (otu[keep, , drop = FALSE] > 0) + 0
    Dk <- vegan::vegdist(X_pa, method = "jaccard", binary = TRUE)
    
    permk <- permute::how(
      blocks = droplevels(factor(meta_ext$natman[keep])),
      nperm = N_PERM_LOO
    )
    
    rhs_loo <- setdiff(rhs_terms, "natman")
    if (length(rhs_loo) < 1) return(NULL)
    
    m <- vegan::adonis2(
      Dk ~ .,
      data = meta_ext[keep, rhs_loo, drop = FALSE],
      by = "margin",
      permutations = permk
    )
    
    mm <- as.data.frame(m)
    mm$term <- rownames(mm)
    
    tibble::as_tibble(mm) %>%
      filter(term %in% intersect(c("ds_at_drill", "depth_2", "Position", "diameter_at_drill_z"), rhs_loo)) %>%
      transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  
  if (nrow(loo_jac) > 0) {
    loo_sum_jac <- loo_jac %>%
      group_by(term) %>%
      summarise(
        median_R2 = median(R2, na.rm = TRUE),
        p05 = quantile(p, 0.05, na.rm = TRUE),
        p95 = quantile(p, 0.95, na.rm = TRUE),
        .groups = "drop"
      )
    
    print(loo_sum_jac)
    readr::write_csv2(
      loo_sum_jac,
      file.path(outdir, paste0("loo_within_tree_summary_Jaccard_", scope, ".csv"))
    )
  }
  
  # ---- Dispersion checks -----------------------------------------------------
  cat("\n[ DISPERSION CHECKS ]\n")
  
  disp_nat <- anova(betadisper(D_robAit, meta_ext$natman))
  disp_mh <- anova(betadisper(D_robAit, meta_ext$Position))
  disp_ds <- anova(betadisper(D_robAit, meta_ext$ds_at_drill))
  disp_size <- anova(betadisper(D_robAit, meta_ext$diameter_at_drill_z))
  disp_depth <- anova(betadisper(D_robAit, droplevels(meta_ext$depth_2)))
  
  sink(file.path(outdir, paste0("permdisp_", scope, ".txt")))
  if (!is.null(disp_nat)) print(disp_nat); print("natman")
  if (!is.null(disp_mh)) print(disp_mh); print("position")
  if (!is.null(disp_ds)) print(disp_ds); print("DS at drill")
  if (!is.null(disp_size)) print(disp_size); print("diam at drill z")
  if (!is.null(disp_depth)) print(disp_depth); print("depth 2")
  sink()
  
  if (!is.null(disp_nat)) { cat("\nDispersion natman:\n"); print(disp_nat) }
  if (!is.null(disp_mh)) { cat("\nDispersion Position:\n"); print(disp_mh) }
  if (!is.null(disp_size)) { cat("\nDispersion size:\n"); print(disp_size) }
  if (!is.null(disp_depth)) { cat("\nDispersion depth_2:\n"); print(disp_depth) }
  if (!is.null(disp_ds)) { cat("\nDispersion ds_at_drill:\n"); print(disp_ds) }
  
  cat("\n[ DISPERSION CHECKS - JACCARD ]\n")
  
  dispJ_nat <- anova(betadisper(D_jaccard, meta_ext$natman))
  dispJ_mh <- anova(betadisper(D_jaccard, meta_ext$Position))
  dispJ_size <- anova(betadisper(D_jaccard, meta_ext$diameter_at_drill_z))
  dispJ_depth <- anova(betadisper(D_jaccard, droplevels(meta_ext$depth_2)))
  dispJ_ds <- anova(betadisper(D_jaccard, droplevels(meta_ext$ds_at_drill)))
  
  sink(file.path(outdir, paste0("permdisp_jaccard_", scope, ".txt")))
  if (!is.null(dispJ_nat)) print(dispJ_nat); print("natman")
  if (!is.null(dispJ_mh)) print(dispJ_mh); print("Position")
  if (!is.null(dispJ_size)) print(dispJ_size); print("diam et drill z")
  if (!is.null(dispJ_depth)) print(dispJ_depth); print("depth 2")
  if (!is.null(dispJ_ds)) print(dispJ_ds); print("ds at drill")
  sink()
  
  if (!is.null(dispJ_nat)) { cat("\nDispersion natman (Jaccard):\n"); print(dispJ_nat) }
  if (!is.null(dispJ_mh)) { cat("\nDispersion Position (Jaccard):\n"); print(dispJ_mh) }
  if (!is.null(dispJ_size)) { cat("\nDispersion diameter_at_drill_z (Jaccard):\n"); print(dispJ_size) }
  if (!is.null(dispJ_depth)) { cat("\nDispersion depth_2 (Jaccard):\n"); print(dispJ_depth) }
  if (!is.null(dispJ_ds)) { cat("\nDispersion ds_at_drill (Jaccard):\n"); print(dispJ_ds) }
  
  print("diameter and ds are overdispersed")
  
  bd_nat <- betadisper(D_robAit, meta_ext$natman)
  tree_disp <- data.frame(
    natman = names(bd_nat$group.distances),
    disp = as.numeric(bd_nat$group.distances)
  )
  
  attr(D_robAit, "Labels") <- meta_ext$sample
  bd_nat <- betadisper(D_robAit, meta_ext$natman)
  
  sample_disp <- data.frame(
    sample = names(bd_nat$distances),
    natman = bd_nat$group,
    dist_to_centroid = as.numeric(bd_nat$distances),
    row.names = NULL
  )
  
  readr::write_csv2(tree_disp, file.path(outdir, paste0("betadisper_tree_centroids_", scope, ".csv")))
  readr::write_csv2(sample_disp, file.path(outdir, paste0("betadisper_sample_distances_", scope, ".csv")))
  
  # ---- Alpha richness GLMM ---------------------------------------------------
  cat("\n[ ALPHA RICHNESS GLMM ]\n")
  
  df <- alpha_df %>%
    left_join(
      meta_ext %>% select(sample, natman, Position, depth_2, ds_at_drill, diameter_at_drill_z, umi),
      by = "sample"
    ) %>%
    mutate(
      natman = droplevels(natman),
      diameter_z = diameter_at_drill_z,
      depth_2 = factor(depth_2, levels = c("INNER", "OUTER")),
      log_reads = ifelse(is.finite(log_reads), log_reads, log1p(reads))
    ) %>%
    filter(!is.na(richness), reads > 0)
  
  rhs_terms_glmm <- c("log_reads", "umi", "(1|natman)", "ds_at_drill", "depth_2", "Position", "diameter_z")
  form <- as.formula(paste("richness ~", paste(rhs_terms_glmm, collapse = " + ")))
  cat("GLMM form: \n")
  print(form)
  
  if (nlevels(df$Position) >= 2 && nlevels(df$ds_at_drill) >= 2) {
    ctrl <- lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    
    m_pois <- suppressWarnings(
      lme4::glmer(form, data = df, family = poisson(link = "log"), control = ctrl)
    )
    
    phi_p <- overdisp_phi(m_pois)
    use_nb <- is.finite(phi_p$phi) && (phi_p$phi > 1.2)
    
    if (use_nb) {
      m_nb <- suppressWarnings(lme4::glmer.nb(form, data = df, control = ctrl))
      best <- if (AIC(m_pois) + 2 < AIC(m_nb)) m_pois else m_nb
      best_type <- if (identical(best, m_pois)) "Poisson (glmer)" else "NB (glmer.nb)"
    } else {
      best <- m_pois
      best_type <- "Poisson (glmer)"
    }
    
    cat(sprintf(
      "Selected alpha model: %s  Overdispersion Poisson phi = %.2f  AIC(best)=%.1f\n",
      best_type, phi_p$phi, AIC(best)
    ))
    
    print(summary(best))
    check_glmm_lme4(best, name = "LOG alpha")
    plot_resid_glmm(best, file.path(outdir, "alpha_richness_LOG_residuals.png"))
    
    an_tab <- car::Anova(best, type = 2)
    r2_tab <- performance::r2(best)
    
    irr_tab <- broom.mixed::tidy(
      best,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = TRUE
    ) %>%
      mutate(
        term = case_when(
          term == "(Intercept)" ~ "Intercept",
          term == "diameter_z" ~ "Diameter (z)",
          term == "log_reads" ~ "log(reads+1)",
          term == "umi1" ~ "Dataset: UMI vs raw",
          grepl("^Position", term) ~ paste0("Position: ", sub("^Position", "", term), " vs ", levels(df$Position)[1]),
          grepl("^depth_2", term) ~ "Depth: OUTER vs INNER",
          grepl("^ds_at_drill", term) ~ paste0("Decay ", sub("^ds_at_drill", "", term), " vs ", levels(df$ds_at_drill)[1]),
          TRUE ~ term
        )
      ) %>%
      rename(
        IRR = estimate,
        CI_low = conf.low,
        CI_high = conf.high,
        p = p.value
      ) %>%
      arrange(term)
    
    emm_pos <- emmeans::emmeans(best, ~Position, type = "response")
    p_emm_pos <- summary(emm_pos) %>%
      as.data.frame() %>%
      ggplot(aes(x = Position, y = response, ymin = asymp.LCL, ymax = asymp.UCL)) +
      geom_point(size = 3) +
      geom_errorbar(width = 0.2) +
      labs(
        x = "Microhabitat",
        y = "Estimated richness",
        title = paste0("Richness ~ microhabitat (INTRA) — ", scope)
      ) +
      theme_classic(base_size = 11)
    
    ggsave(
      file.path(outdir, paste0("alpha_richness_emm_Position_", scope, ".png")),
      p_emm_pos, width = 6, height = 3.2, dpi = 300
    )
    print(p_emm_pos)
    
    emm_ds <- emmeans::emmeans(best, ~ds_at_drill, type = "response")
    ds_df <- summary(emm_ds) %>% as.data.frame()
    
    p_emm_ds <- ggplot(ds_df, aes(x = ds_at_drill, y = response)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
      labs(
        x = "Decay at drill (ordered 1→5)",
        y = "Estimated richness",
        title = paste0("Richness ~ decay (INTRA) — ", scope)
      ) +
      theme_classic(base_size = 11)
    
    ggsave(
      file.path(outdir, paste0("alpha_richness_emm_DS_", scope, ".png")),
      p_emm_ds, width = 6.2, height = 3.2, dpi = 300
    )
    print(p_emm_ds)
    
    emm_depth <- emmeans::emmeans(best, ~depth_2, type = "response")
    depth_df <- summary(emm_depth) %>% as.data.frame()
    
    p_emm_depth <- ggplot(depth_df, aes(x = depth_2, y = response)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
      labs(
        x = "Depth (INNER/OUTER)",
        y = "Estimated richness",
        title = paste0("Richness ~ depth (INTRA) — ", scope)
      ) +
      theme_classic(base_size = 11)
    
    ggsave(
      file.path(outdir, paste0("alpha_richness_emm_depth2_", scope, ".png")),
      p_emm_depth, width = 4.6, height = 3.2, dpi = 300
    )
    print(p_emm_depth)
    
    keep_eff <- c("Position", "depth_2", "ds_at_drill", "diameter_z", "log_reads")
    
    tbl_main <- tibble(
      Effect = rownames(an_tab),
      Chisq = an_tab$`Chisq`,
      df = an_tab$Df,
      p_value = an_tab$`Pr(>Chisq)`
    ) %>%
      filter(Effect %in% keep_eff) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(match(Effect, keep_eff))
    
    tbl_r2 <- tibble(
      Metric = c("R2_marginal", "R2_conditional"),
      Value = c(r2_tab$R2_marginal, r2_tab$R2_conditional)
    )
    
    cat("\n=== Tables ===\n")
    print(list(anova = tbl_main, r2 = tbl_r2, irr = irr_tab))
    
    capture.output(
      list(model = summary(best), anova = an_tab, r2 = r2_tab, irr = irr_tab),
      file = file.path(outdir, paste0("alpha_glmm_", scope, ".txt"))
    )
    
    readr::write_csv2(tbl_main, file.path(outdir, paste0("alpha_glmm_anova_", scope, ".csv")))
    readr::write_csv2(tbl_r2, file.path(outdir, paste0("alpha_glmm_r2_", scope, ".csv")))
    readr::write_csv2(irr_tab, file.path(outdir, paste0("alpha_glmm_irr_", scope, ".csv")))
    
    print(emmeans::emmeans(best, ~depth_2, type = "response"))
    print(emmeans::emmeans(best, ~ds_at_drill, type = "response"))
    print(emmeans::emmeans(best, ~Position, type = "response"))
    print(emmeans::emmeans(best, ~umi, type = "response"))
    print(emmeans::emmeans(best, ~diameter_z, type = "response"))
  }
  
  # ---- Betapart decomposition ------------------------------------------------
  cat("\n[ BETAPART — turnover and nestedness ]\n")
  
  PA <- (otu > 0) * 1L
  bp <- betapart::beta.pair(PA, index.family = "sorensen")
  
  DM_sim <- as.matrix(bp$beta.sim)
  DM_sne <- as.matrix(bp$beta.sne)
  DM_sor <- as.matrix(bp$beta.sor)
  
  pairs_sim <- mk_within_between_tbl(DM_sim, meta_ext$natman)
  pairs_sne <- mk_within_between_tbl(DM_sne, meta_ext$natman)
  pairs_sor <- mk_within_between_tbl(DM_sor, meta_ext$natman)
  
  sum_sim <- pairs_sim %>% group_by(type) %>% summarise(n = n(), median = stats::median(dist), IQR = IQR(dist), .groups = "drop")
  sum_sne <- pairs_sne %>% group_by(type) %>% summarise(n = n(), median = stats::median(dist), IQR = IQR(dist), .groups = "drop")
  sum_sor <- pairs_sor %>% group_by(type) %>% summarise(n = n(), median = stats::median(dist), IQR = IQR(dist), .groups = "drop")
  
  bd_sim <- boot_diff_median(pairs_sim, seed = SEED_GLOBAL)
  ci_sim <- quantile(bd_sim, c(0.025, 0.975), na.rm = TRUE)
  
  bd_sne <- boot_diff_median(pairs_sne, seed = SEED_GLOBAL)
  ci_sne <- quantile(bd_sne, c(0.025, 0.975), na.rm = TRUE)
  
  bd_sor <- boot_diff_median(pairs_sor, seed = SEED_GLOBAL)
  ci_sor <- quantile(bd_sor, c(0.025, 0.975), na.rm = TRUE)
  
  est_sim <- with(sum_sim, median[type == "between"] - median[type == "within"])
  est_sne <- with(sum_sne, median[type == "between"] - median[type == "within"])
  est_sor <- with(sum_sor, median[type == "between"] - median[type == "within"])
  
  beta_table <- bind_rows(
    tibble(component = "beta.SIM", sum_sim, est = est_sim, ci_lo = ci_sim[1], ci_hi = ci_sim[2]),
    tibble(component = "beta.SNE", sum_sne, est = est_sne, ci_lo = ci_sne[1], ci_hi = ci_sne[2]),
    tibble(component = "beta.SOR", sum_sor, est = est_sor, ci_lo = ci_sor[1], ci_hi = ci_sor[2])
  )
  
  readr::write_csv2(
    beta_table,
    file.path(outdir, paste0("betapart_within_between_summary_", scope, ".csv"))
  )
  print(beta_table)

