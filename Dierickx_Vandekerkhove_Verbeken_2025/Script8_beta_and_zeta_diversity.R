# =============================================================================
# Script 8 — Cross-scope β-diversity & zeta diversity
# - Preps meta (maps FWD to aFWD/fFWD), aligns OTUs (+taxa labels)
# - Betapart: turnover vs nestedness (incidence)
# - Dispersion: distances to centroid by dw_type
# - PERMANOVA (unblocked + natman-blocked) across distances
#   + console readout of natman R² across distances
# - LCBD/SCBD on Hellinger; figures + compact tables
# - Zeta diversity: within vs between trees
#   * per-tree normalisation, retention, exponential/power fits + AIC
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(betapart); library(adespatial)
  library(permute);  library(ggplot2); library(readr)
})
suppressPackageStartupMessages({
  library(glmmTMB); library(emmeans); library(performance)
})
set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# ---- Guards & dirs -----------------------------------------------------------
stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"))
if (!exists("exclude_samples")) exclude_samples <- character(0)
dir.create("plots/ALL", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/ZETA_unscaled", recursive = TRUE, showWarnings = FALSE)

# ---- Palettes (optional fallbacks) ------------------------------------------
dw_cols <- dw_colors

# ---- Prepare META with attached/fallen FWD distinction ----------------------
meta_filt <- META1 %>%
  filter(!sample %in% exclude_samples) %>%
  mutate(
    dw_type = case_when(
      dw_type == "FWD" & position %in% c("ATTACHED","ATTACHTED") ~ "aFWD",
      dw_type == "FWD" & position == "FALLEN"                    ~ "fFWD",
      TRUE ~ dw_type
    ),
    natman  = droplevels(factor(natman)),
    umi=as.factor(umi)) %>%
  filter(dw_type %in% c("SNAG","aFWD","fFWD","LOG")) %>%
  droplevels()

# ---- Align OTU table + append lowest_taxon to column names ------------------
X <- otu_matrix_filt[meta_filt$sample, , drop = FALSE]
stopifnot(identical(rownames(X), meta_filt$sample))
stopifnot(all(c("sh_code","lowest_taxon") %in% names(tax)))
X <- X[, colSums(X) > 0, drop = FALSE]

# Add read depths
meta_filt <- meta_filt %>%
  mutate(reads = rowSums(X), log_reads = log1p(reads)) %>%
  filter(reads > 0) %>% droplevels()
X <- X[meta_filt$sample, , drop = FALSE]  # keep in sync

# ---- Incidence matrix for betapart ------------------------------------------
pa <- (X > 0) * 1
grp <- (meta_filt$dw_type)
## =============================================================================
## 0b) Alpha diversity (richness) by dw_type2 + models (NB; with/without natman)
## =============================================================================


# Ensure dw_type2 exists (explicit name) and consistent factor order/colors
dw_levels <- c("SNAG","aFWD","fFWD","LOG")
meta_filt <- meta_filt %>%
  mutate(dw_type2 = factor(dw_type, levels = dw_levels)) %>%
  droplevels()

# Compute alpha metrics
S <- rowSums(pa)                           # richness (incidence)
H <- vegan::diversity(X, index = "shannon")
J <- H / log(pmax(S, 2))                   # Pielou evenness (guard S<=1)

alpha_df <- meta_filt %>%
  transmute(sample, natman, dw_type2, umi,
            reads, log_reads,
            richness = S[ sample ],
            shannon  = H[ sample ],
            pielou   = J[ sample ]) %>%
  filter(!is.na(richness) & !is.na(dw_type2)) %>%
  droplevels()

# ---- Quick exploratory plot (boxplot) ---------------------------------------
p_rich <- ggplot(alpha_df, aes(dw_type2, richness, fill = dw_type2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width= 0.8) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.25, size = 1) +
  scale_fill_manual(values = dw_cols[levels(alpha_df$dw_type2)]) +
  labs(title = "Alpha richness by dead-wood type",
       x = NULL, y = "Observed richness (SH count)") +
  theme_minimal(12) + theme(legend.position = "none")
p_rich
ggsave("plots/ALL/alpha_richness_by_dw_type2.png", p_rich, width = 4, height = 4, dpi = 600)

# ---- Modeling richness -------------------------------------------------------
# Use negative binomial to handle overdispersion; include read depth as covariate.
# (Random intercept for natman provides the 'blocked' within-tree structure.)
m_nb_blocked   <- glmmTMB(richness ~ dw_type2 + scale(log_reads) + (1|natman),
                          family = nbinom2, data = alpha_df)
m_nb_unblocked <- glmmTMB(richness ~ dw_type2 + scale(log_reads),
                          family = nbinom2, data = alpha_df)

cat("\nAlpha richness model (NB, blocked by natman):\n")
print(summary(m_nb_blocked))
cat("\nOverdispersion check (blocked):\n")
print(performance::check_overdispersion(m_nb_blocked))
cat("\nR2 (blocked):\n")
print(performance::r2(m_nb_blocked))

cat("\nAlpha richness model (NB, unblocked):\n")
print(summary(m_nb_unblocked))
cat("\nOverdispersion check (unblocked):\n")
print(performance::check_overdispersion(m_nb_unblocked))
m_nb_unblocked_null <- glmmTMB(richness ~ scale(log_reads),
                               family = nbinom2, data = alpha_df)

r2_mcfadden <- 1 - as.numeric(logLik(m_nb_unblocked) / logLik(m_nb_unblocked_null))

n <- nobs(m_nb_unblocked)
LL0 <- as.numeric(logLik(m_nb_unblocked_null))
LL1 <- as.numeric(logLik(m_nb_unblocked))

# Cox & Snell and Nagelkerke (using log-likelihoods)
r2_coxsnell  <- 1 - exp((2 / n) * (LL0 - LL1))
r2_nagelkerke <- r2_coxsnell / (1 - exp((2 / n) * LL0))

cat("\nPseudo-R2 (unblocked):\n")
print(tibble::tibble(R2_McFadden = r2_mcfadden,
                     R2_CoxSnell = r2_coxsnell,
                     R2_Nagelkerke = r2_nagelkerke))

# ---- Estimated marginal means (response scale) + pairwise tests --------------
emm_blk <- emmeans(m_nb_blocked, ~ dw_type2, type = "response")
emm_df  <- emm_blk %>% summary(infer = TRUE) %>% as.data.frame()
pairs_blk <- pairs(emm_blk, adjust = "tukey") %>% as.data.frame()

# Optional: write to disk
write_csv(emm_df,  "plots/ALL/alpha_richness_emmeans_dw_type2_blocked.csv")
write_csv(pairs_blk,"plots/ALL/alpha_richness_pairwise_dw_type2_blocked.csv")

# Ensure the factor order you want on the x-axis
dw_levels <- c("LOG","aFWD","fFWD","SNAG")
emm_df <- emm_df %>%
  dplyr::mutate(dw_type2 = factor(dw_type2, levels = dw_levels))

# Build a palette just for the levels present in emm_df
pal <- dw_colors[levels(emm_df$dw_type2)]

# Plot predicted means + 95% CI (response scale), colored by dw_type2
p_emm <- ggplot(emm_df,
                aes(x = dw_type2, y = response,
                    ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_col(aes(fill = dw_type2)) +
  geom_errorbar(width = 0.2) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "Estimated richness by dead-wood type (blocked NB model)",
       x = NULL, y = "Predicted richness (mean ± 95% CI)") +
  theme_bw(12)

p_emm
ggsave("plots/ALL/alpha_richness_emmeans_dw_type2_blocked.png",
       p_emm, width = 4, height = 4, dpi = 600)


# =============================================================================
# 1) Turnover vs nestedness (Sorensen family; incidence)
# =============================================================================
beta_all_multi <- betapart::beta.multi.abund(X)
print(beta_all_multi)

beta_all <- betapart::beta.multi(pa, index.family = "sorensen")
# turnover = beta.SIM; nestedness = beta.SNE
cat("\nSorensen partitions (all samples):\n"); print(beta_all)
# for completeness pairwise
beta_pair <- betapart::beta.pair(pa, index.family = "sorensen")

# Save compact summary
write_lines(c(
  paste0("beta.SIM (turnover): ", signif(beta_all$beta.SIM, 4)),
  paste0("beta.SNE (nestedness): ", signif(beta_all$beta.SNE, 4)),
  paste0("beta.SOR (total): ", signif(beta_all$beta.SOR, 4))
), "plots/ALL/betapart_summary.txt")


# Fungal communities are highly dissimilar across dead-wood samples, with average pairwise β-diversity values approaching unity in both abundance- and incidence-based analyses.
# Abundance partitioning of the Bray–Curtis dissimilarity (β<sub>BRAY</sub> = 0.994) showed that nearly all variation (99.1%) was attributable to balanced changes in relative abundances among taxa,
# whereas abundance gradients contributed negligibly (0.4%).
# Similarly, presence–absence partitioning of Sørensen dissimilarity (β<sub>SOR</sub> = 0.992) revealed that almost all variation (98.3%) arose from species turnover,
# with only a minor nestedness component (0.9%).
# Together, these results indicate that fungal assemblages are not structured as nested subsets of richer communities, but are instead shaped by strong taxon replacement across dead-wood types and microhabitats, consistent with deterministic filtering by substrate properties and decay stage.

# =============================================================================
# 2) PERMANOVA + Dispersion (by dw_type)
# =============================================================================
# Distances
D_robAit      <- vegdist(X, method = "robust.aitchison")
X_hel         <- decostand(X, method = "hellinger")
D_euclid_hell <- dist(X_hel, method = "euclidean")
X_rel         <- decostand(X, method = "total")
D_bray        <- vegdist(X_rel, method = "bray")
D_brayS       <- vegdist(sqrt(X_rel), method = "bray")
D_jacc        <- vegdist(pa, method = "jaccard", binary = TRUE)

# Unblocked & natman-blocked models
nperm <- 999
perm_unblocked <- how(nperm = nperm)
perm_blocked   <- how(nperm = nperm); setBlocks(perm_blocked) <- meta_filt$natman

perma_unblocked <- list(
  robAit  = adonis2(D_robAit      ~ log_reads + umi  + natman + dw_type, data = meta_filt, by = "margin", permutations = perm_unblocked),
  hellEuc = adonis2(D_euclid_hell ~ log_reads + umi  + natman + dw_type, data = meta_filt, by = "margin", permutations = perm_unblocked),
  bray    = adonis2(D_bray        ~ log_reads +umi  +  natman + dw_type, data = meta_filt, by = "margin", permutations = perm_unblocked),
  brayS   = adonis2(D_brayS       ~ log_reads +umi  +  natman + dw_type, data = meta_filt, by = "margin", permutations = perm_unblocked),
  jaccard = adonis2(D_jacc        ~ umi  + natman + dw_type,             data = meta_filt, by = "margin", permutations = perm_unblocked)
)

cat("\nUNBLOCKED PERMANOVA (key readouts):\n")
print(perma_unblocked$robAit); 
print(perma_unblocked$bray);
print(perma_unblocked$jaccard)

# natman R² across distances (console + CSV)
natman_R2 <- purrr::map_dfr(perma_unblocked, ~{
  as_tibble(as.data.frame(.x), rownames = "term")
}, .id = "distance") %>%
  filter(term == "natman") %>%
  transmute(distance, R2, p = `Pr(>F)`) %>%
  arrange(desc(R2))
cat("\nUnblocked PERMANOVA — natman R² by distance:\n"); print(natman_R2)
write_csv(natman_R2, "plots/ALL/permanova_natman_R2_across_distances.csv")

# Blocked (within-tree) effect of dw_type
perma_blocked <- adonis2(D_robAit ~ log_reads + umi+ dw_type , data = meta_filt, by = "margin", permutations = perm_blocked)
cat("\nBLOCKED PERMANOVA (robust Aitchison) — within-tree dw_type:\n"); print(perma_blocked)

# Leave-one-tree-out stability for blocked model
trees <- meta_filt$natman %>% as.character() %>% unique() %>% setdiff(NA)
loo <- purrr::map_dfr(trees, function(t) {
  keep <- meta_filt$natman != t
  if (sum(keep) < 5) return(NULL)
  Dk <- vegdist(X[keep, , drop = FALSE], method = "robust.aitchison")
  permk <- how(nperm = 499); setBlocks(permk) <- droplevels(meta_filt$natman[keep])
  m <- adonis2(Dk ~ log_reads + umi  +  dw_type, data = meta_filt[keep, ], by = "margin", permutations = permk)
  as_tibble(as.data.frame(m), rownames = "term") %>%
    filter(term == "dw_type") %>%
    mutate(left_out = t) %>%
    select(left_out, term, R2, p = `Pr(>F)`)
})
if (nrow(loo)) {
  cat("\nLOO stability — blocked PERMANOVA term dw_type:\n")
  print(loo %>% summarise(med_R2 = median(R2, na.rm=TRUE),
                          p05 = quantile(p, 0.05, na.rm=TRUE),
                          p95 = quantile(p, 0.95, na.rm=TRUE)))
  write_csv(loo, "plots/ALL/loo_blocked_permanova_dw_type.csv")
}

# Dispersion by group (distance to centroid)
bd_ait <- betadisper(D_robAit, grp)
cat("\nPERMDISP (Aitchison) by dw_type:\n"); print(permutest(bd_ait, permutations = perm_unblocked))
# quick plot
p_bd <- ggplot(tibble(dist = as.numeric(bd_ait$distances),
                      dw_type = grp),
               aes(dw_type, dist, fill = dw_type)) +
  geom_violin(trim=FALSE, alpha=0.55) +
  geom_boxplot(width=0.18, outlier.shape = NA) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.5, size = 1.4) +
  scale_fill_manual(values = dw_cols) +
  labs(title = "Distances to centroid (PERMDISP, Aitchison)", y = "Distance", x = NULL) +
  theme_minimal(12) + theme(legend.position = "none")
p_bd
ggsave("plots/ALL/permdisp_dw_type_robAit.png", p_bd, width = 6.6, height = 4.2, dpi = 300)

# Also check for Bray and Jaccard
bd_bray  <- betadisper(D_bray, group = meta_filt$dw_type)
perm_b   <- permutest(bd_bray, permutations = 999)
print(perm_b)
bd_jac   <- betadisper(D_jacc, group = meta_filt$dw_type)
perm_j   <- permutest(bd_jac, permutations = 999)
print(perm_j)
plot(bd_bray);
plot(bd_jac)
# =============================================================================
# 3) LCBD/SCBD (Hellinger)
# =============================================================================
# Local Contributions to Beta Diversity (LCBD) and Species Contributions to Beta Diversity (SCBD) metrics
bd_lcbd <- adespatial::beta.div(X_hel, method = "euclidean")
LCBD <- bd_lcbd$LCBD; pLCBD <- bd_lcbd$p.LCBD; SCBD <- bd_lcbd$SCBD
scbd_df <- tibble(feature = colnames(X_hel), SCBD = as.numeric(SCBD)) %>% arrange(desc(SCBD))
lcbd_df <- tibble(sample = rownames(X_hel),
                  LCBD = as.numeric(LCBD),
                  p_LCBD = as.numeric(pLCBD),
                  dw_type = grp,
                  natman = meta_filt$natman,
                  reads = meta_filt$reads) %>% 
  dplyr::mutate(
    p_LCBD = ifelse(is.na(p_LCBD), 1, p_LCBD),     # guard against NA p's
    q_BH   = p.adjust(p_LCBD, method = "BH"),
    sig_0.05 = q_BH < 0.05,
    sig_0.10 = q_BH < 0.10) %>% 
  arrange(p_LCBD)
lcbd_df
cat("----- NO samples are significantly contributing more to beta-div!! ---")

# Plots + tables
scbd_topN <- scbd_df %>% slice_head(n = 40)
p_scbd <- ggplot(scbd_topN, aes(x = reorder(feature, SCBD), y = SCBD)) +
  geom_col() + coord_flip() +
  labs(x = "SH|lowest_taxon", y = "SCBD", title = "Top SCBD drivers (Hellinger)") +
  theme_bw(11)
p_scbd
ggsave("plots/ALL/scbd_top40.png", p_scbd, width = 9, height = 7.5, dpi = 300)

p_lcbd <- ggplot(lcbd_df, aes(dw_type, LCBD, fill = dw_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.45) +
  scale_fill_manual(values = dw_cols) +
  labs(title = "LCBD by deadwood type (Hellinger)", x = NULL, y = "LCBD") +
  theme_minimal(12) + theme(legend.position = "none")
p_lcbd
ggsave("plots/ALL/lcbd_by_dw_type.png", p_lcbd, width = 6.6, height = 4.2, dpi = 300)

write_csv(lcbd_df, "plots/ALL/LCBD_by_sample.csv")
write_csv(scbd_df, "plots/ALL/SCBD_by_feature.csv")

# =============================================================================
# 4) ZETA DIVERSITY — normalised + retention + model fits
# =============================================================================
suppressPackageStartupMessages({ library(zetadiv) })

# ---- Run requested sets ------------------------------------------------------
run_zeta_norm(meta_filt, X, c("LOG"),                         "LOG")
run_zeta_norm(meta_filt, X, c("LOG","aFWD"),                  "LOG_aFWD")
run_zeta_norm(meta_filt, X, c("LOG","aFWD","SNAG"),           "LOG_aFWD_SNAG")
run_zeta_norm(meta_filt, X, c("LOG","aFWD","fFWD","SNAG"),    "tree-level")
run_zeta_norm(meta_filt, X, c("SNAG"),                        "SNAG")
run_zeta_norm(meta_filt, X, c("aFWD","fFWD"),                 "FWD")
run_zeta_norm(meta_filt, X, c("aFWD"),                        "aFWD")
run_zeta_norm(meta_filt, X, c("fFWD"),                        "fFWD")

run_zeta_abs(meta_filt, otu_matrix_filt, c("LOG"),                      "LOG")
run_zeta_abs(meta_filt, otu_matrix_filt, c("LOG","aFWD"),               "LOG_aFWD")
run_zeta_abs(meta_filt, otu_matrix_filt, c("LOG","aFWD","SNAG"),        "LOG_aFWD_SNAG")
run_zeta_abs(meta_filt, otu_matrix_filt, c("LOG","aFWD","fFWD","SNAG"), "tree-level")
run_zeta_abs(meta_filt, otu_matrix_filt, c("SNAG"),                     "SNAG")
run_zeta_abs(meta_filt, otu_matrix_filt, c("aFWD","fFWD"),              "FWD")
run_zeta_abs(meta_filt, otu_matrix_filt, c("aFWD"),                     "aFWD")
run_zeta_abs(meta_filt, otu_matrix_filt, c("fFWD"),                     "fFWD")


cat("\n=== Script 8 — done. Outputs in plots/ALL and plots/ZETA_unscaled ===\n")


## =============================================================================
## 5) PCA with species loadings in Aitchison geometry (CLR after robust zero-repl)
##     + depth diagnostics, partial PCA (Condition on log_reads), balanced biplot
## =============================================================================
suppressPackageStartupMessages({
  library(vegan)
  library(tidyverse)
  library(ggrepel)
  library(scales)
  library(zCompositions)
  library(compositions)
})
if (requireNamespace("zCompositions", quietly = TRUE)) {
  X_clr <- compositions::clr(zCompositions::cmultRepl(as.matrix(X), method = "CZM", label = 0))
  pca_ait <- vegan::rda(X_clr)
  eig_pca <- vegan::eigenvals(pca_ait); var_expl_pca <- eig_pca / sum(eig_pca)
  pc1 <- scales::percent(var_expl_pca[1]); pc2 <- scales::percent(var_expl_pca[2])
  
  sites_pca <- as_tibble(vegan::scores(pca_ait, display = "sites", scaling = 2), rownames = "sample") %>%
    left_join(meta_filt %>% select(sample, dw_type, natman), by = "sample")
  
  p_pca_ait <- ggplot(sites_pca, aes(PC1, PC2, color = dw_type)) +
    geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.4) +
    geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.4) +
    geom_point(alpha = 0.85, size = 2) +
    scale_color_manual(values = dw_cols) +
    labs(title = "PCA — CLR (Aitchison geometry)",
         x = paste0("PC1 (", pc1, ")"), y = paste0("PC2 (", pc2, ")")) +
    theme_minimal(12)
  p_pca_ait
  #  ggsave(file.path(outdir, "PCA_CLR_Aitchison.png"), p_pca_ait, width = 6.6, height = 4.6, dpi = 300)
} else {
  message("Optional CLR-PCA skipped: install.packages('zCompositions','compositions') to enable.")
}

# ------------------------ Parameters you can tweak -----------------------------
outdir <- "plots/ALL"
tabdir <- "tables/ALL"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir,  recursive = TRUE, showWarnings = FALSE)

scaling_choice <- 1        # vegan biplot scaling: 1 (distance) or 2 (correlation)
topN_total     <- 20L      # total species to annotate in balanced biplot (≈ 4*K below)
K_per_sign     <- 7      # per-axis-per-sign (PC1+, PC1-, PC2+, PC2-); 12*4 = 48 (deduped to ~<=40)
do_flip_axes   <- TRUE     # flip PC1/PC2 for display consistency (purely cosmetic)
prevalence_min <- 0L       # optional prevalence filter BEFORE CLR (e.g., 2 → keep taxa present in ≥2 samples)

# ----------------------------- Helpers ----------------------------------------
# Robust CLR with NO deletions; optional prevalence filter upstream
CLR_from_X <- function(X, prevalence_min = 0L) {
  X <- as.matrix(X)
  if (prevalence_min > 0L) {
    keep <- colSums(X > 0) >= prevalence_min
    X <- X[, keep, drop = FALSE]
  }
  Xz <- zCompositions::cmultRepl(
    X, method = "CZM", label = 0,
    z.warning = TRUE, z.delete = FALSE  # <-- do NOT drop rows/cols
  )
  compositions::clr(Xz)
}

# Balanced species selection by sign on PC1/PC2
balanced_species <- function(sp, K = 10L) {
  sp_pos_pc1 <- sp %>% arrange(desc(PC1)) %>% slice_head(n = K)
  sp_neg_pc1 <- sp %>% arrange(PC1)        %>% slice_head(n = K)
  sp_pos_pc2 <- sp %>% arrange(desc(PC2)) %>% slice_head(n = K)
  sp_neg_pc2 <- sp %>% arrange(PC2)        %>% slice_head(n = K)
  bind_rows(sp_pos_pc1, sp_neg_pc1, sp_pos_pc2, sp_neg_pc2) %>%
    distinct(feature, .keep_all = TRUE) %>%
    mutate(vec_len = sqrt(PC1^2 + PC2^2)) %>%
    arrange(desc(vec_len))
}

# ------------------------------ Data prep -------------------------------------
stopifnot(exists("X"), exists("meta_filt"), exists("dw_cols"))
meta_filt <- meta_filt %>% mutate(log_reads = log1p(reads))

# CLR transform (Aitchison Euclidean space)
X_clr <- CLR_from_X(X, prevalence_min = prevalence_min)

# Strict alignment of rows with metadata
common <- intersect(rownames(X_clr), meta_filt$sample)
X_clr  <- X_clr[common, , drop = FALSE]
meta2  <- meta_filt %>%
  filter(sample %in% common) %>%
  arrange(match(sample, rownames(X_clr)))
stopifnot(identical(rownames(X_clr), meta2$sample))

# ----------------------- Plain CLR-PCA (diagnostics) --------------------------
pca_ait <- vegan::rda(X_clr)  # unconstrained RDA == PCA
eig_pca   <- vegan::eigenvals(pca_ait)
var_expl  <- eig_pca / sum(eig_pca)
pc1_lab   <- percent(var_expl[1]); pc2_lab <- percent(var_expl[2])

sites_scr   <- as_tibble(vegan::scores(pca_ait, display = "sites",   scaling = scaling_choice), rownames = "sample") %>%
  left_join(dplyr::select(meta2, sample, dw_type, natman, reads, log_reads), by = "sample")
species_scr <- as_tibble(vegan::scores(pca_ait, display = "species", scaling = scaling_choice), rownames = "feature")

# Depth correlations (raw)
c_pc1 <- suppressWarnings(cor(sites_scr$PC1, sites_scr$log_reads))
c_pc2 <- suppressWarnings(cor(sites_scr$PC2, sites_scr$log_reads))
message(sprintf("CLR-PCA (raw): cor(PC1, log_reads)=%.3f; cor(PC2, log_reads)=%.3f", c_pc1, c_pc2))

# Optional axis flips for display (no effect on statistics)
if (do_flip_axes) {
  # Orient PC1 so that mean(LOG) > mean(others); edit group as desired
  grp_right <- "LOG"
  flip_pc1 <- sign(mean(sites_scr$PC1[sites_scr$dw_type == grp_right], na.rm = TRUE) -
                     mean(sites_scr$PC1[sites_scr$dw_type != grp_right], na.rm = TRUE))
  flip_pc1[is.na(flip_pc1) | flip_pc1 == 0] <- 1
  sites_scr$PC1   <- flip_pc1 * sites_scr$PC1
  species_scr$PC1 <- flip_pc1 * species_scr$PC1
  # Keep PC2 as-is; set flip_pc2 <- -1 if you prefer “upward” direction for a group.
}

# Species strength
species_scr <- species_scr %>%
  mutate(loading_len = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(loading_len))

# ----------------------------- Plots (raw) ------------------------------------
# Sites-only
p_pca_sites <- ggplot(sites_scr, aes(PC1, PC2, color = dw_type)) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = dw_cols, guide = guide_legend(override.aes = list(size = 3))) +
  labs(
    title = sprintf("PCA — CLR (Aitchison): Sites [PC1=%s, PC2=%s]", pc1_lab, pc2_lab),
    x = "PC1", y = "PC2", color = "Dead-wood type"
  ) +
  theme_minimal(12) + theme(legend.position = "right")
p_pca_sites
ggsave(file.path(outdir, "PCA_CLR_Aitchison_sites.png"), p_pca_sites, width = 6.6, height = 4.6, dpi = 300)

# Balanced biplot (top ± on each axis)
species_bal <- balanced_species(species_scr, K = K_per_sign)
# Arrow scale relative to site spread
rng <- range(c(sites_scr$PC1, sites_scr$PC2), na.rm = TRUE)
arrow_mult <- 0.9 * diff(rng) / max(1e-9, max(sqrt(species_bal$PC1^2 + species_bal$PC2^2)))

p_pca_biplot <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_point(data = sites_scr, aes(PC1, PC2, color = dw_type), alpha = 0.9, size = 2) +
  scale_color_manual(values = dw_cols, guide = guide_legend(override.aes = list(size = 3))) +
  geom_segment(
    data = species_bal,
    aes(x = 0, y = 0, xend = PC1 * arrow_mult, yend = PC2 * arrow_mult),
    arrow = arrow(length = unit(0.015, "npc")), linewidth = 0.3, alpha = 0.85
  ) +
  ggrepel::geom_text_repel(
    data = species_bal,
    aes(x = PC1 * arrow_mult, y = PC2 * arrow_mult, label = feature),
    size = 2.6, segment.size = 0.2, max.overlaps = 200, box.padding = 0.3
  ) +
  labs(
    title = sprintf("PCA — CLR (Aitchison): balanced biplot (scaling=%d)", scaling_choice),
    x = "PC1", y = "PC2", color = "Dead-wood type"
  ) +
  theme_minimal(12) + theme(legend.position = "right")
p_pca_biplot
ggsave(file.path(outdir, "PCA_CLR_Aitchison_biplot_balanced.png"), p_pca_biplot, width = 8.0, height = 6.0, dpi = 300)

# Scree plot
scree_df <- tibble(axis = seq_along(var_expl), var_expl = as.numeric(var_expl))
p_scree <- ggplot(scree_df[1:10, ], aes(axis, var_expl)) +
  geom_col() + scale_y_continuous(labels = percent_format()) +
  labs(title = "CLR-PCA: Scree (first 10 PCs)", x = "PC", y = "Variance explained") +
  theme_minimal(12)
ggsave(file.path(outdir, "PCA_CLR_Aitchison_scree.png"), p_scree, width = 5.8, height = 3.6, dpi = 300)


# ---------------- Partial CLR-PCA (remove library-size effect) ----------------
pca_part <- vegan::rda(X_clr ~ 1 + Condition(log_reads), data = meta2)
ep  <- vegan::eigenvals(pca_part); vep <- ep / sum(ep)
pc1p <- percent(vep[1]); pc2p <- percent(vep[2])

Sp <- as_tibble(vegan::scores(pca_part, display = "sites",   scaling = scaling_choice), rownames = "sample") %>%
  left_join(dplyr::select(meta2, sample, dw_type, natman, reads, log_reads), by = "sample")
Lp <- as_tibble(vegan::scores(pca_part, display = "species", scaling = scaling_choice), rownames = "feature")

# Optional axis flips applied identically for comparability
if (do_flip_axes) {
  flip_pc1 <- sign(mean(Sp$PC1[Sp$dw_type == grp_right], na.rm = TRUE) -
                     mean(Sp$PC1[Sp$dw_type != grp_right], na.rm = TRUE))
  flip_pc1[is.na(flip_pc1) | flip_pc1 == 0] <- 1
  Sp$PC1 <- flip_pc1 * Sp$PC1; Lp$PC1 <- flip_pc1 * Lp$PC1
}

# Depth correlations after partialling out (should be ~0)
cp1 <- suppressWarnings(cor(Sp$PC1, Sp$log_reads))
cp2 <- suppressWarnings(cor(Sp$PC2, Sp$log_reads))
message(sprintf("CLR-PCA (partial): cor(PC1, log_reads)=%.3f; cor(PC2, log_reads)=%.3f", cp1, cp2))

# Balanced species for partial biplot
Lp <- Lp %>% mutate(loading_len = sqrt(PC1^2 + PC2^2))
Lp_bal <- balanced_species(Lp, K = K_per_sign)

rngp <- range(c(Sp$PC1, Sp$PC2), na.rm = TRUE)
arrow_mult_p <- 0.9 * diff(rngp) / max(1e-9, max(sqrt(Lp_bal$PC1^2 + Lp_bal$PC2^2)))

p_sites_part <- ggplot(Sp, aes(PC1, PC2, color = dw_type)) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = dw_cols) +
  labs(
    title = sprintf("CLR-PCA (partial on log_reads): Sites [PC1=%s, PC2=%s]", pc1p, pc2p),
    x = "PC1", y = "PC2", color = "Dead-wood type"
  ) +
  theme_minimal(12) + theme(legend.position = "right")
p_sites_part
ggsave(file.path(outdir, "PCA_CLR_sites_partial_logreads.png"), p_sites_part, width = 6.6, height = 4.6, dpi = 300)

p_biplot_part <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_vline(xintercept = 0, linewidth = 0.3, alpha = 0.35) +
  geom_point(data = Sp, aes(PC1, PC2, color = dw_type), alpha = 0.9, size = 2) +
  scale_color_manual(values = dw_cols) +
  geom_segment(
    data = Lp_bal,
    aes(x = 0, y = 0, xend = PC1 * arrow_mult_p, yend = PC2 * arrow_mult_p),
    arrow = arrow(length = unit(0.015, "npc")), linewidth = 0.3, alpha = 0.85
  ) +
  ggrepel::geom_text_repel(
    data = Lp_bal,
    aes(x = PC1 * arrow_mult_p, y = PC2 * arrow_mult_p, label = feature),
    size = 2.6, segment.size = 0.2, max.overlaps = 200, box.padding = 0.3
  ) +
  labs(
    title = sprintf("CLR-PCA (partial on log_reads): balanced biplot (scaling=%d)", scaling_choice),
    x = "PC1", y = "PC2", color = "Dead-wood type"
  ) +
  theme_minimal(12) + theme(legend.position = "right")
p_biplot_part
ggsave(file.path(outdir, "PCA_CLR_biplot_partial_balanced.png"), p_biplot_part, width = 8.0, height = 6.0, dpi = 300)

