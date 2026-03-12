# =============================================================================
# SCOPE 1: CROSS deadwood beta diversity
# =============================================================================
source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(betapart)
  library(adespatial)
  library(permute)
  library(ggplot2)
  library(readr)
  library(scales)
  library(emmeans)
  library(car)
  library(broom.mixed)
  library(performance)
  library(lme4)
})

set.seed(SEED_GLOBAL)
setwd(".")

# -----------------------------------------------------------------------------
# 0) Guards and directories
# -----------------------------------------------------------------------------
stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"))

scope <- "CROSS"

plot_dir      <- file.path("plots", "ALL")
table_dir     <- file.path("tables", "ALL")
beta_plot_dir <- file.path(plot_dir, "beta")
beta_tab_dir  <- file.path(table_dir, "beta")
null_plot_dir <- file.path(plot_dir, "nulls")
null_tab_dir  <- file.path(table_dir, "nulls")

ensure_dir(plot_dir)
ensure_dir(table_dir)
ensure_dir(beta_plot_dir)
ensure_dir(beta_tab_dir)
ensure_dir(null_plot_dir)
ensure_dir(null_tab_dir)

N_PERM_MAIN   <- get0("N_PERM_MAIN", ifnotfound = 999L)
N_PERM_LOO    <- 499L
N_PERM_NULL   <- 199L
N_PERM_ANOSIM <- 9999L
B_NULL        <- 999L

dw_levels <- c("SNAG", "aFWD", "fFWD", "LOG")
dw_cols   <- dw_colors

# -----------------------------------------------------------------------------
# 1) Metadata and OTU alignment
# -----------------------------------------------------------------------------
meta_filt <- META1 %>%
  filter(!sample %in% exclude_samples) %>%
  mutate(
    sample   = as.character(sample),
    dw_type2 = factor(dw_type2, levels = dw_levels),
    natman   = droplevels(factor(natman)),
    umi      = droplevels(factor(umi))
  ) %>%
  filter(dw_type2 %in% dw_levels) %>%
  droplevels()

X <- otu_matrix_filt[meta_filt$sample, , drop = FALSE]
stopifnot(identical(rownames(X), meta_filt$sample))

X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]

meta_filt <- meta_filt %>%
  mutate(
    reads     = rowSums(X),
    log_reads = log1p(reads)) %>%
  filter(reads > 0) %>%
  droplevels()

X <- X[meta_filt$sample, , drop = FALSE]
X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]

stopifnot(identical(rownames(X), meta_filt$sample))

pa <- (X > 0) * 1L

cat("\n[CROSS] sample counts by dw_type2\n")
print(table(meta_filt$dw_type2))
cat("Samples:", nrow(meta_filt), " | SHs:", ncol(X), "\n")

# optional taxon-labelled copy for LCBD/SCBD outputs only
stopifnot(all(c("sh_code", "lowest_taxon") %in% names(tax)))
lut <- setNames(tax$lowest_taxon, tax$sh_code)
colnames_X_tax <- paste0(
  colnames(X), "|",
  ifelse(is.na(lut[colnames(X)]), "", lut[colnames(X)])
)

# -----------------------------------------------------------------------------
# 2) Alpha diversity (q0 richness)
# -----------------------------------------------------------------------------
cat("\n[ALPHA RICHNESS GLMM - CROSS]\n")

alpha_df <- alpha_from_otu(X) %>%
  left_join(
    meta_filt %>%
      select(sample, natman, dw_type2, umi),
    by = "sample") %>%
  mutate(
    log_reads = log1p(reads),
    natman    = droplevels(factor(natman)),
    dw_type2  = droplevels(factor(dw_type2, levels = dw_levels)),
    umi       = droplevels(factor(umi))) %>%
  filter(
    !is.na(richness),
    !is.na(dw_type2),
    reads > 0) %>%
  droplevels()

print(table(alpha_df$dw_type2))

rhs_terms_glmm_cross <- c("log_reads", "umi", "(1|natman)", "dw_type2")
form_cross <- as.formula(paste("richness ~", paste(rhs_terms_glmm_cross, collapse = " + ")))
print(form_cross)

ctrl_pois <- lme4::glmerControl(
  optimizer = "bobyqa",
  optCtrl   = list(maxfun = 1e5)
)

m_pois_cross <- suppressWarnings(
  lme4::glmer(
    form_cross,
    data    = alpha_df,
    family  = poisson(link = "log"),
    control = ctrl_pois
  )
)

phi_p_cross <- overdisp_phi(m_pois_cross)

best_cross <- m_pois_cross
best_type_cross <- "Poisson (glmer)"

if (is.finite(phi_p_cross$phi) && phi_p_cross$phi > 1.2) {
  
  ctrl_nb <- lme4::glmerControl(
    optimizer = "nloptwrap",
    optCtrl   = list(maxeval = 2e5)
  )
  
  m_nb_cross <- tryCatch(
    suppressWarnings(
      lme4::glmer.nb(
        form_cross,
        data    = alpha_df,
        control = ctrl_nb
      )
    ),
    error = function(e) NULL
  )
  
  ok_nb <- !is.null(m_nb_cross) &&
    isTRUE(length(warnings()) >= 0) &&
    isFALSE(lme4::isSingular(m_nb_cross, tol = 1e-4))
  
  if (ok_nb) {
    if (AIC(m_nb_cross) + 2 < AIC(m_pois_cross)) {
      best_cross <- m_nb_cross
      best_type_cross <- "NB (glmer.nb)"
    }
  } else {
    m_nb_tmb_cross <- tryCatch(
      suppressWarnings(
        glmmTMB::glmmTMB(
          formula = form_cross,
          data    = alpha_df,
          family  = glmmTMB::nbinom2(link = "log")
        )
      ),
      error = function(e) NULL
    )
    
    if (!is.null(m_nb_tmb_cross)) {
      if (AIC(m_nb_tmb_cross) + 2 < AIC(m_pois_cross)) {
        best_cross <- m_nb_tmb_cross
        best_type_cross <- "NB (glmmTMB nbinom2)"
      }
    }
  }
}

cat(sprintf(
  "Selected CROSS alpha model: %s | Poisson phi = %.2f | AIC(best) = %.1f\n",
  best_type_cross, phi_p_cross$phi, AIC(best_cross)
))
print(summary(best_cross))
check_glmm_lme4(best_cross, name = "CROSS alpha")
plot_resid_glmm(best_cross, file.path(plot_dir, "alpha_richness_CROSS_residuals.png"))

an_tab_cross <- car::Anova(best_cross, type = 2)
r2_cross     <- performance::r2(best_cross)

print(an_tab_cross)
print(r2_cross)

tab_q0_fix_cross <- broom.mixed::tidy(
  best_cross,
  effects      = "fixed",
  conf.int     = TRUE,
  exponentiate = TRUE
) %>%
  mutate(
    term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "log_reads"   ~ "log(reads+1)",
      grepl("^umi", term)   ~ paste0("Dataset: ", sub("^umi", "", term), " vs ", levels(alpha_df$umi)[1]),
      grepl("^dw_type2", term) ~ paste0("Substrate: ", sub("^dw_type2", "", term), " vs ", levels(alpha_df$dw_type2)[1]),
      TRUE ~ term
    ),
    `95% CI` = sprintf("%.3f to %.3f", conf.low, conf.high)
  ) %>%
  transmute(
    Scope     = scope,
    Term      = term,
    IRR       = round(estimate, 3),
    `95% CI`,
    z         = round(statistic, 3),
    `p-value` = signif(p.value, 3)
  )

print(tab_q0_fix_cross)
write_csv2(tab_q0_fix_cross, file.path(table_dir, "alpha_richness_IRR_CROSS.csv"))

tab_anova_cross <- tibble(
  Effect = rownames(an_tab_cross),
  Chisq  = an_tab_cross[["Chisq"]],
  df     = an_tab_cross[["Df"]],
  p      = an_tab_cross[["Pr(>Chisq)"]]
)
write_csv2(tab_anova_cross, file.path(table_dir, "alpha_richness_Anova_CROSS.csv"))

tab_r2_cross <- tibble(
  Metric = c("R2_marginal", "R2_conditional"),
  Value  = c(r2_cross$R2_marginal, r2_cross$R2_conditional)
)
write_csv2(tab_r2_cross, file.path(table_dir, "alpha_richness_R2_CROSS.csv"))

emm_cross <- emmeans::emmeans(best_cross, ~ dw_type2, type = "response")
emm_df_cross <- summary(emm_cross, infer = TRUE) %>% as.data.frame()
write_csv2(emm_df_cross, file.path(table_dir, "alpha_richness_emmeans_dw_type2_CROSS.csv"))

pairs_cross <- pairs(emm_cross, adjust = "tukey") %>% as.data.frame()
print(pairs_cross)
write_csv2(pairs_cross, file.path(table_dir, "alpha_richness_pairwise_dw_type2_CROSS.csv"))

emm_df_cross <- emm_df_cross %>%
  mutate(dw_type2 = factor(dw_type2, levels = c("LOG", "aFWD", "fFWD", "SNAG")))

pal_emm <- dw_cols[levels(emm_df_cross$dw_type2)]

p_emm <- ggplot(
  emm_df_cross,
  aes(x = dw_type2, y = response, ymin = asymp.LCL, ymax = asymp.UCL)
) +
  geom_col(aes(fill = dw_type2)) +
  geom_errorbar(width = 0.2) +
  scale_fill_manual(values = pal_emm, guide = "none") +
  labs(
    title = "Estimated SH richness by deadwood type",
    x     = NULL,
    y     = "Predicted richness (95% CI)"
  ) +
  theme_bw(base_size = 12)

print(p_emm)
ggsave(file.path(plot_dir, "alpha_richness_emmeans_dw_type2_CROSS.png"),
       p_emm, width = 4.2, height = 4.2, dpi = 600)

# -----------------------------------------------------------------------------
# 3) Distances
# -----------------------------------------------------------------------------
cat("\n[BETA DIVERSITY - CROSS]\n")

D_robAit <- vegdist(X, method = "robust.aitchison")
D_robAit <- fix_dist_labels(D_robAit, samps = rownames(X))

X_hel  <- decostand(X, method = "hellinger")
D_jacc <- vegdist(pa, method = "jaccard", binary = TRUE)

# -----------------------------------------------------------------------------
# 4) Betapart
# -----------------------------------------------------------------------------
beta_all_pa <- betapart::beta.multi(pa, index.family = "sorensen")
cat("\nBetapart multisite partition (CROSS)\n")
print(beta_all_pa)

beta_multi_summary <- tibble(
  scope  = scope,
  metric = c("beta.SOR", "beta.SIM", "beta.SNE"),
  value  = c(beta_all_pa$beta.SOR, beta_all_pa$beta.SIM, beta_all_pa$beta.SNE)
)
write_csv2(beta_multi_summary, file.path(beta_tab_dir, "betapart_multi_CROSS.csv"))

beta_pair_pa <- betapart::beta.pair(pa, index.family = "sorensen")
labs <- rownames(pa)
dw   <- meta_filt$dw_type2[match(labs, meta_filt$sample)]

DM_sim <- as.matrix(beta_pair_pa$beta.sim)
DM_sne <- as.matrix(beta_pair_pa$beta.sne)
DM_sor <- as.matrix(beta_pair_pa$beta.sor)

upper <- which(upper.tri(DM_sim), arr.ind = TRUE)

pairs_df <- tibble(
  i         = labs[upper[, 1]],
  j         = labs[upper[, 2]],
  dw_i      = dw[upper[, 1]],
  dw_j      = dw[upper[, 2]],
  pair_class = factor(
    if_else(dw_i == dw_j, "within_substrate", "between_substrate"),
    levels = c("within_substrate", "between_substrate")
  ),
  beta_sim  = DM_sim[upper],
  beta_sne  = DM_sne[upper],
  beta_sor  = DM_sor[upper]
)

betapart_summary <- pairs_df %>%
  group_by(pair_class) %>%
  summarise(
    n_pairs         = n(),
    beta_sim_median = median(beta_sim, na.rm = TRUE),
    beta_sne_median = median(beta_sne, na.rm = TRUE),
    beta_sor_median = median(beta_sor, na.rm = TRUE),
    ratio_turnover_nestedness = beta_sim_median / pmax(beta_sne_median, .Machine$double.eps),
    .groups = "drop"
  )

print(betapart_summary)
write_csv2(betapart_summary, file.path(beta_tab_dir, "betapart_within_between_dwtype_CROSS.csv"))

# -----------------------------------------------------------------------------
# 5) PERMANOVA
# -----------------------------------------------------------------------------
perm_unblocked <- how(nperm = N_PERM_MAIN)
rhs_terms <- c("log_reads", "umi", "natman", "dw_type2")

perma_unblocked <- list(
  robAit = adonis2(
    D_robAit ~ .,
    data         = meta_filt[, rhs_terms, drop = FALSE],
    by           = "terms",
    permutations = perm_unblocked
  ),
  jaccard = adonis2(
    D_jacc ~ .,
    data         = meta_filt[, rhs_terms, drop = FALSE],
    by           = "terms",
    permutations = perm_unblocked
  )
)

cat("\nUNBLOCKED PERMANOVA (CROSS)\n")
print(perma_unblocked$robAit)
print(perma_unblocked$jaccard)

natman_R2 <- purrr::map_dfr(perma_unblocked, ~{
  df <- as.data.frame(.x)
  df$term <- rownames(df)
  tibble::as_tibble(df)
}, .id = "distance") %>%
  filter(term == "natman") %>%
  select(distance, R2, p = `Pr(>F)`) %>%
  arrange(desc(R2))

print(natman_R2)
write_csv2(natman_R2, file.path(beta_tab_dir, "permanova_natman_R2_across_distances_CROSS.csv"))

perm_blocked <- how(nperm = N_PERM_MAIN)
setBlocks(perm_blocked) <- meta_filt$natman
rhs_block <- setdiff(rhs_terms, "natman")

perma_blocked <- list(
  robAit = adonis2(
    D_robAit ~ .,
    data         = meta_filt[, rhs_block, drop = FALSE],
    by           = "margin",
    permutations = perm_blocked
  ),
  jaccard = adonis2(
    D_jacc ~ .,
    data         = meta_filt[, rhs_block, drop = FALSE],
    by           = "margin",
    permutations = perm_blocked
  )
)

cat("\nBLOCKED PERMANOVA (within tree)\n")
print(perma_blocked$robAit)
print(perma_blocked$jaccard)

perma_blocked_tab <- purrr::map_dfr(perma_blocked, ~{
  df <- as.data.frame(.x)
  df$term <- rownames(df)
  tibble::as_tibble(df)
}, .id = "distance")
write_csv2(perma_blocked_tab, file.path(beta_tab_dir, "permanova_blocked_summary_CROSS.csv"))

# -----------------------------------------------------------------------------
# 6) Leave-one-tree-out robustness
# -----------------------------------------------------------------------------
trees <- unique(as.character(meta_filt$natman))
trees <- trees[!is.na(trees)]

loo_robAit <- purrr::map_dfr(trees, function(t) {
  keep <- meta_filt$natman != t
  if (sum(keep) < 5) return(NULL)
  
  Xk <- X[keep, , drop = FALSE]
  Dk <- vegdist(Xk, method = "robust.aitchison")
  
  permk <- how(nperm = N_PERM_LOO)
  setBlocks(permk) <- droplevels(meta_filt$natman[keep])
  
  mk <- adonis2(
    Dk ~ log_reads + umi + dw_type2,
    data         = meta_filt[keep, ],
    by           = "margin",
    permutations = permk
  )
  
  as_tibble(as.data.frame(mk), rownames = "term") %>%
    filter(term %in% c("dw_type2", "umi", "log_reads")) %>%
    mutate(left_out = t) %>%
    select(left_out, term, R2, p = `Pr(>F)`)
})

if (nrow(loo_robAit) > 0) {
  loo_robAit_sum <- loo_robAit %>%
    group_by(term) %>%
    summarise(
      med_R2 = median(R2, na.rm = TRUE),
      p05    = quantile(p, 0.05, na.rm = TRUE),
      p95    = quantile(p, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
  print(loo_robAit_sum)
  write_csv2(loo_robAit, file.path(beta_tab_dir, "loo_blocked_permanova_dw_type2_robAit_CROSS.csv"))
}

loo_jacc <- purrr::map_dfr(trees, function(t) {
  keep <- meta_filt$natman != t
  if (sum(keep) < 5) return(NULL)
  
  Xk <- X[keep, , drop = FALSE]
  Dk <- vegdist((Xk > 0) * 1L, method = "jaccard", binary = TRUE)
  
  permk <- how(nperm = N_PERM_LOO)
  setBlocks(permk) <- droplevels(meta_filt$natman[keep])
  
  mk <- adonis2(
    Dk ~ log_reads + umi + dw_type2,
    data         = meta_filt[keep, ],
    by           = "margin",
    permutations = permk
  )
  
  as_tibble(as.data.frame(mk), rownames = "term") %>%
    filter(term == "dw_type2") %>%
    mutate(left_out = t) %>%
    select(left_out, term, R2, p = `Pr(>F)`)
})

if (nrow(loo_jacc) > 0) {
  loo_jacc_sum <- loo_jacc %>%
    summarise(
      med_R2 = median(R2, na.rm = TRUE),
      p05    = quantile(p, 0.05, na.rm = TRUE),
      p95    = quantile(p, 0.95, na.rm = TRUE)
    )
  print(loo_jacc_sum)
  write_csv2(loo_jacc, file.path(beta_tab_dir, "loo_blocked_permanova_dw_type2_jaccard_CROSS.csv"))
}

# -----------------------------------------------------------------------------
# 7) PERMDISP
# -----------------------------------------------------------------------------
bd_ait <- betadisper(D_robAit, meta_filt$dw_type2)
cat("\nPERMDISP (Aitchison)\n")
print(permutest(bd_ait, permutations = N_PERM_MAIN))

bd_jac <- betadisper(D_jacc, meta_filt$dw_type2)
cat("\nPERMDISP (Jaccard)\n")
print(permutest(bd_jac, permutations = N_PERM_MAIN))

# -----------------------------------------------------------------------------
# 8) LCBD / SCBD
# -----------------------------------------------------------------------------
bd_lcbd <- adespatial::beta.div(X_hel, method = "euclidean")
LCBD    <- bd_lcbd$LCBD
pLCBD   <- bd_lcbd$p.LCBD
SCBD    <- bd_lcbd$SCBD

scbd_df <- tibble(
  feature = colnames_X_tax,
  SCBD    = as.numeric(SCBD)
) %>%
  arrange(desc(SCBD))

lcbd_df <- tibble(
  sample  = rownames(X_hel),
  LCBD    = as.numeric(LCBD),
  p_LCBD  = as.numeric(pLCBD),
  dw_type = meta_filt$dw_type2,
  natman  = meta_filt$natman,
  reads   = meta_filt$reads
) %>%
  mutate(
    p_LCBD   = ifelse(is.na(p_LCBD), 1, p_LCBD),
    q_BH     = p.adjust(p_LCBD, method = "BH"),
    sig_0.05 = q_BH < 0.05,
    sig_0.10 = q_BH < 0.10
  ) %>%
  arrange(p_LCBD)

print(head(lcbd_df, 10))

scbd_topN <- scbd_df %>% slice_head(n = 20)
p_scbd <- ggplot(scbd_topN, aes(x = reorder(feature, SCBD), y = SCBD)) +
  geom_col() +
  coord_flip() +
  labs(
    x = "SH|lowest_taxon",
    y = "SCBD",
    title = "Top SCBD drivers (Hellinger)"
  ) +
  theme_bw(base_size = 11)

print(p_scbd)
ggsave(file.path(beta_plot_dir, "scbd_top20_CROSS.png"), p_scbd, width = 6, height = 8, dpi = 600)

p_lcbd <- ggplot(lcbd_df, aes(dw_type, LCBD, fill = dw_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.45) +
  scale_fill_manual(values = dw_cols[levels(meta_filt$dw_type2)]) +
  labs(
    title = "LCBD by deadwood type (Hellinger)",
    x = NULL,
    y = "LCBD"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(p_lcbd)
ggsave(file.path(beta_plot_dir, "lcbd_by_dw_type2_CROSS.png"), p_lcbd, width = 6.6, height = 4.2, dpi = 300)

write_csv2(lcbd_df, file.path(beta_tab_dir, "LCBD_by_sample_CROSS.csv"))
write_csv2(scbd_df, file.path(beta_tab_dir, "SCBD_by_feature_CROSS.csv"))

# -----------------------------------------------------------------------------
# 9) Null models for natman and dw_type2
# -----------------------------------------------------------------------------
cat("\n[NULL MODELS - CROSS]\n")

get_R2 <- function(fit, term) {
  df <- as.data.frame(fit)
  if (!(term %in% rownames(df))) return(NA_real_)
  df[term, "R2"]
}

upper_idx_from <- function(D) which(upper.tri(as.matrix(D)), arr.ind = TRUE)

pairwise_groups_local <- function(D, grp, U = NULL) {
  M <- as.matrix(D)
  if (is.null(U)) U <- which(upper.tri(M), arr.ind = TRUE)
  dU   <- M[cbind(U[, 1], U[, 2])]
  same <- grp[U[, 1]] == grp[U[, 2]]
  list(within = dU[same], between = dU[!same])
}

auc_from_group_local <- function(D, grp, U = NULL) {
  pw <- pairwise_groups_local(D, grp, U)
  if (length(pw$within) < 2 || length(pw$between) < 2) {
    return(c(AUC = NA_real_, r_rb = NA_real_))
  }
  r  <- rank(c(pw$between, pw$within))
  n1 <- length(pw$between)
  n2 <- length(pw$within)
  Ustat <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  A <- Ustat / (n1 * n2)
  c(AUC = A, r_rb = 2 * A - 1)
}

boot_diff_median_local <- function(x_within, x_between, nboot = 1999, seed = 42) {
  set.seed(seed)
  nW <- length(x_within)
  nB <- length(x_between)
  if (nW == 0 || nB == 0) return(rep(NA_real_, nboot))
  replicate(nboot, {
    med_within  <- median(sample(x_within,  nW, replace = TRUE))
    med_between <- median(sample(x_between, nB, replace = TRUE))
    med_between - med_within
  })
}

mk_ecdf_plot <- function(D, grp, title, xlab) {
  U  <- upper_idx_from(D)
  pw <- pairwise_groups_local(D, grp, U)
  
  df <- tibble(
    dist = c(pw$within, pw$between),
    type = factor(rep(c("within", "between"),
                      c(length(pw$within), length(pw$between))),
                  levels = c("within", "between"))
  )
  
  sumtab <- df %>%
    group_by(type) %>%
    summarise(n = n(), med = median(dist), .groups = "drop")
  
  bd <- boot_diff_median_local(
    df$dist[df$type == "within"],
    df$dist[df$type == "between"],
    nboot = 1999
  )
  ci  <- quantile(bd, c(0.025, 0.975), na.rm = TRUE)
  est <- unname(sumtab$med[sumtab$type == "between"] -
                  sumtab$med[sumtab$type == "within"])
  
  xpos <- quantile(df$dist, 0.70, na.rm = TRUE)
  ypos <- 0.2
  lab  <- sprintf("Δ median = %.3f\n95%% CI [%.3f, %.3f]", est, ci[1], ci[2])
  
  ggplot(df, aes(x = dist, colour = type)) +
    stat_ecdf(geom = "step", linewidth = 1) +
    geom_vline(xintercept = sumtab$med[sumtab$type == "within"], linetype = 2, alpha = 0.6) +
    geom_vline(xintercept = sumtab$med[sumtab$type == "between"], linetype = 2, alpha = 0.6) +
    annotate("text", x = xpos, y = ypos, hjust = 0, label = lab) +
    theme_bw() +
    labs(title = title, x = xlab, y = "ECDF", colour = NULL)
}

mk_density_plot <- function(null_vals, obs, title, xlab) {
  ci  <- quantile(null_vals, c(0.025, 0.975), na.rm = TRUE)
  ses <- (obs - mean(null_vals, na.rm = TRUE)) / sd(null_vals, na.rm = TRUE)
  p   <- mean(null_vals >= obs, na.rm = TRUE)
  
  ggplot(tibble(x = null_vals), aes(x)) +
    geom_density(fill = "grey85") +
    geom_vline(xintercept = obs, linetype = 2) +
    geom_vline(xintercept = ci, linetype = 3) +
    theme_bw() +
    labs(
      title = title,
      subtitle = sprintf(
        "obs = %.3f  null mean = %.3f  95%% CI [%.3f, %.3f]  p = %.3f  SES = %.2f",
        obs, mean(null_vals, na.rm = TRUE), ci[1], ci[2], p, ses
      ),
      x = xlab,
      y = "Density"
    )
}

prepare_inputs <- function(D, meta, vars) {
  m <- meta
  m$umi      <- factor(m$umi)
  m$natman   <- factor(m$natman)
  m$dw_type2 <- factor(m$dw_type2)
  m$log_reads <- ifelse(is.finite(m$log_reads), m$log_reads, log1p(m$reads))
  ok   <- complete.cases(m[, vars, drop = FALSE])
  m_ok <- droplevels(m[ok, , drop = FALSE])
  
  labs <- attr(D, "Labels")
  if (!is.null(labs) && all(m_ok$sample %in% labs)) {
    ord <- match(m_ok$sample, labs)
    M   <- as.matrix(D)[ord, ord, drop = FALSE]
    D_ok <- as.dist(M)
  } else {
    D_ok <- D
  }
  
  list(D = D_ok, meta = m_ok)
}

run_metric <- function(D_raw, meta_raw, metric_label, xlab_short) {
  vars <- c("log_reads", "umi", "natman", "dw_type2")
  pp   <- prepare_inputs(D_raw, meta_raw, vars)
  D    <- pp$D
  meta <- pp$meta
  
  use_reads  <- identical(metric_label, "robust.aitchison")
  base_terms <- c(if (use_reads) "log_reads", "umi")
  
  form_obs <- as.formula(
    sprintf("D ~ %s + natman + dw_type2", paste(base_terms, collapse = " + "))
  )
  form_nat <- as.formula(
    sprintf("D ~ %s + nat_null + dw_type2", paste(base_terms, collapse = " + "))
  )
  form_typ <- as.formula(
    sprintf("D ~ %s + natman + type_null", paste(base_terms, collapse = " + "))
  )
  
  perma_obs_nat <- adonis2(form_obs, data = meta, by = "terms",  permutations = N_PERM_NULL)
  perma_obs_typ <- adonis2(form_obs, data = meta, by = "margin", permutations = N_PERM_NULL)
  
  R2_nat_obs <- get_R2(perma_obs_nat, "natman")
  R2_typ_obs <- get_R2(perma_obs_typ, "dw_type2")
  
  null_nat <- numeric(B_NULL)
  null_typ <- numeric(B_NULL)
  
  for (b in seq_len(B_NULL)) {
    df_loop <- meta %>%
      mutate(
        nat_null  = factor(sample(natman, replace = FALSE), levels = levels(natman)),
        type_null = factor(sample(dw_type2, replace = FALSE), levels = levels(dw_type2))
      )
    
    m_nat <- tryCatch(
      adonis2(form_nat, data = df_loop, by = "terms", permutations = N_PERM_NULL),
      error = function(e) e
    )
    null_nat[b] <- if (inherits(m_nat, "error")) NA_real_ else get_R2(m_nat, "nat_null")
    
    m_typ <- tryCatch(
      adonis2(form_typ, data = df_loop, by = "margin", permutations = N_PERM_NULL),
      error = function(e) e
    )
    null_typ[b] <- if (inherits(m_typ, "error")) NA_real_ else get_R2(m_typ, "type_null")
  }
  
  U           <- upper_idx_from(D)
  auc_nat_obs <- auc_from_group_local(D, meta$natman, U)
  auc_typ_obs <- auc_from_group_local(D, meta$dw_type2, U)
  
  AUC_null_nat <- replicate(B_NULL, auc_from_group_local(D, sample(meta$natman), U)["AUC"])
  AUC_null_typ <- replicate(B_NULL, auc_from_group_local(D, sample(meta$dw_type2), U)["AUC"])
  
  anos_nat <- anosim(D, meta$natman,  permutations = N_PERM_ANOSIM)
  anos_typ <- anosim(D, meta$dw_type2, permutations = N_PERM_ANOSIM)
  
  p_ecdf_nat <- mk_ecdf_plot(D, meta$natman,
                             title = paste0("ECDF within vs between natman - ", metric_label),
                             xlab = paste(xlab_short, "distance"))
  p_ecdf_typ <- mk_ecdf_plot(D, meta$dw_type2,
                             title = paste0("ECDF within vs between dw_type2 - ", metric_label),
                             xlab = paste(xlab_short, "distance"))
  
  p_r2_nat <- mk_density_plot(null_nat, R2_nat_obs,
                              title = paste0("R2 null vs observed natman - ", metric_label),
                              xlab = "R2")
  p_r2_typ <- mk_density_plot(null_typ, R2_typ_obs,
                              title = paste0("R2 null vs observed dw_type2 - ", metric_label),
                              xlab = "R2")
  p_auc_nat <- mk_density_plot(AUC_null_nat, unname(auc_nat_obs["AUC"]),
                               title = paste0("AUC null vs observed natman - ", metric_label),
                               xlab = "AUC")
  p_auc_typ <- mk_density_plot(AUC_null_typ, unname(auc_typ_obs["AUC"]),
                               title = paste0("AUC null vs observed dw_type2 - ", metric_label),
                               xlab = "AUC")
  
  print(p_ecdf_nat); print(p_ecdf_typ)
  print(p_r2_nat);   print(p_r2_typ)
  print(p_auc_nat);  print(p_auc_typ)
  
  ggsave(file.path(null_plot_dir, paste0("ECDF_natman_", metric_label, "_CROSS.png")),
         p_ecdf_nat, width = 7.2, height = 4.2, dpi = 300)
  ggsave(file.path(null_plot_dir, paste0("ECDF_dwtype_", metric_label, "_CROSS.png")),
         p_ecdf_typ, width = 7.2, height = 4.2, dpi = 300)
  ggsave(file.path(null_plot_dir, paste0("R2_null_natman_", metric_label, "_CROSS.png")),
         p_r2_nat, width = 6.2, height = 4.2, dpi = 300)
  ggsave(file.path(null_plot_dir, paste0("R2_null_dwtype_", metric_label, "_CROSS.png")),
         p_r2_typ, width = 6.2, height = 4.2, dpi = 300)
  ggsave(file.path(null_plot_dir, paste0("AUC_null_natman_", metric_label, "_CROSS.png")),
         p_auc_nat, width = 6.2, height = 4.2, dpi = 300)
  ggsave(file.path(null_plot_dir, paste0("AUC_null_dwtype_", metric_label, "_CROSS.png")),
         p_auc_typ, width = 6.2, height = 4.2, dpi = 300)
  
  make_row <- function(term, R2_obs, null_R2, AUC_obs, null_AUC, anos) {
    ciR2  <- quantile(null_R2, c(0.025, 0.975), na.rm = TRUE)
    ciAUC <- quantile(null_AUC, c(0.025, 0.975), na.rm = TRUE)
    
    tibble(
      Scope = scope,
      Metric = metric_label,
      Term = term,
      `R2 obs` = R2_obs,
      `R2 null mean` = mean(null_R2, na.rm = TRUE),
      `R2 null 95% CI` = sprintf("[%.3f, %.3f]", ciR2[1], ciR2[2]),
      `p (R2)` = mean(null_R2 >= R2_obs, na.rm = TRUE),
      `R2 SES` = (R2_obs - mean(null_R2, na.rm = TRUE)) / sd(null_R2, na.rm = TRUE),
      `R2 percentile` = mean(null_R2 <= R2_obs, na.rm = TRUE),
      `R2 corrected (obs minus mean null)` = R2_obs - mean(null_R2, na.rm = TRUE),
      `AUC obs` = unname(AUC_obs["AUC"]),
      r_rb = unname(AUC_obs["r_rb"]),
      `AUC null mean` = mean(null_AUC, na.rm = TRUE),
      `AUC null 95% CI` = sprintf("[%.3f, %.3f]", ciAUC[1], ciAUC[2]),
      `p (AUC)` = mean(null_AUC >= AUC_obs["AUC"], na.rm = TRUE),
      `AUC SES` = (AUC_obs["AUC"] - mean(null_AUC, na.rm = TRUE)) / sd(null_AUC, na.rm = TRUE),
      `ANOSIM R` = anos$statistic,
      `ANOSIM p` = anos$signif,
      `n trees` = nlevels(meta$natman),
      `n samples` = nrow(meta)
    )
  }
  
  bind_rows(
    make_row("natman",  R2_nat_obs, null_nat, auc_nat_obs, AUC_null_nat, anos_nat),
    make_row("dw_type2", R2_typ_obs, null_typ, auc_typ_obs, AUC_null_typ, anos_typ)
  )
}

tbl_ait <- run_metric(D_robAit, meta_filt, "robust.aitchison", "Aitchison")
tbl_jac <- run_metric(D_jacc,   meta_filt, "jaccard",          "Jaccard")

final_table <- bind_rows(tbl_ait, tbl_jac)
print(final_table, n = Inf)

write_csv2(final_table, file.path(null_tab_dir, "nullmodel_summary_natman_dwtype_CROSS.csv"))