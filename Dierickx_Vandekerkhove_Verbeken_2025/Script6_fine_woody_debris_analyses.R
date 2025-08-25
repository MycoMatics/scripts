# =============================================================================
# Script 6 — FWD analyses (ATTACHED vs FALLEN), tied by natman
# Scopes: all FWD (attached + fallen), with fallen-only subset + matched-natman subset
#
# Adds compared to prior draft:
#   - Read-depth control everywhere (covariate or Condition(log_reads))
#   - PERMANOVA (unblocked & natman-blocked) across multiple distances
#   - Console readout of natman R² across distances
#   - Use size (factor) instead of diameter/size_class
#   - Leave-one-tree-out (LOO) stability for natman-blocked PERMANOVA
#   - Alpha GLMM model selection (Poisson vs NB), prints chosen model
#   - Alpha EMM plots for all factors (fwd_state, decay_stage, size)
#   - IRR table + compact publication tables
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(forcats)
  library(vegan);     library(permute)
  library(glmmTMB);   library(emmeans);  library(performance); library(broom.mixed)
  library(ggplot2)
})
set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")
outdir <- "plots/FWD"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Palettes (fallbacks if not defined upstream) ----
size_colors<- sc_colors
fwd_colors <- c(ATTACHED=dw_colors[["aFWD"]], FALLEN=dw_colors[["fFWD"]])

# -----------------------------------------------------------------------------
# 0) Choose OTU source (filtered preferred)
# -----------------------------------------------------------------------------
otu_source <- otu_matrix_filt
# Keep only valid FWD states (attached/fallen); stop if others creep in
fwd <- dplyr::filter(META1, dw_type == "FWD")
bad_levels <- setdiff(unique(fwd$position), c("ATTACHED","FALLEN"))
stopifnot(length(bad_levels) == 0)  # fail fast if CROWN/ENDO/etc. appear

# -----------------------------------------------------------------------------
# 1) Build FWD metadata + align OTUs
# -----------------------------------------------------------------------------
FWD_meta <- META1 %>%
  filter(dw_type == "FWD") %>%                                   # FWD only
  mutate(
    # collapse to EARLY/LATE for FWD (you noted this specifically)
    decay_stage = factor(decay_stage, levels = c("EARLY", "AVERAGE", "LATE")),
    ds_at_drill = factor(ds_at_drill, levels = c("0","1","2","3","4","5")),
    size        = droplevels(factor(size)),
    # robust ATTACHED vs FALLEN
    fwd_state = as.character(position))

# Align OTU and add reads
otu_fwd <- otu_source[FWD_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fwd), FWD_meta$sample))
otu_fwd <- otu_fwd[, colSums(otu_fwd) > 0, drop = FALSE]

FWD_meta <- FWD_meta %>%
  mutate(reads = rowSums(otu_fwd), log_reads = log1p(reads)) %>%
  filter(!is.na(fwd_state), reads > 0) %>%
  mutate(fwd_state = droplevels(factor(fwd_state)),
         natman    = droplevels(factor(natman)),
         size      = droplevels(size))
otu_fwd <- otu_fwd[FWD_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fwd), FWD_meta$sample))

cat("FWD samples:", nrow(FWD_meta), " | OTUs:", ncol(otu_fwd), "\n")
cat("Levels — natman:", nlevels(FWD_meta$natman),
    "| fwd_state:", levels(FWD_meta$fwd_state),
    "| decay:", levels(FWD_meta$decay_stage),
    "| size:", paste(levels(FWD_meta$size), collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 2) Distances (computed once, after filtering)
# -----------------------------------------------------------------------------
D_robAit      <- vegdist(otu_fwd, method = "robust.aitchison")
comm_hell     <- decostand(otu_fwd, method = "hellinger")
D_euclid_hell <- dist(comm_hell, method = "euclidean")
comm_rel      <- decostand(otu_fwd, method = "total")
D_bray_rel    <- vegdist(comm_rel, method = "bray")
D_bray_relS   <- vegdist(sqrt(comm_rel), method = "bray")  # dominants down-weighted
D_jacc        <- vegdist((otu_fwd > 0) * 1, method = "jaccard", binary = TRUE)

# -----------------------------------------------------------------------------
# 3) Alpha diversity — descriptive (by fwd_state, decay_stage, size)
# -----------------------------------------------------------------------------
alpha_fwd <- tibble(
  sample   = rownames(otu_fwd),
  richness = vegan::specnumber(otu_fwd),
  shannon  = vegan::diversity(otu_fwd, index = "shannon"),
  simpson  = vegan::diversity(otu_fwd, index = "simpson")
) %>%
  left_join(FWD_meta, by = "sample") %>%
  pivot_longer(c(richness, shannon, simpson), names_to = "metric", values_to = "value")

# state
p_alpha_state <- alpha_fwd %>%
  filter(umi == "0") %>% 
  ggplot(aes(fwd_state, value, fill = fwd_state)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
  scale_fill_manual(values = fwd_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Alpha diversity by FWD state", x = "FWD state", y = "Value", fill = "FWD state")
p_alpha_state
ggsave(file.path(outdir, "alpha_by_fwd_state.png"), p_alpha_state, width=9, height=4.8, dpi=300)

# decay
p_alpha_decay <- alpha_fwd %>%
  filter(umi == "0") %>% 
  ggplot(aes(decay_stage, value, fill = decay_stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
  scale_fill_manual(values = DS_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "Alpha diversity by decay stage", x = "Decay stage", y = "Value", fill = "Decay stage")
p_alpha_decay
ggsave(file.path(outdir, "alpha_by_decay_stage.png"), p_alpha_decay, width=9, height=4.8, dpi=300)

# size
p_alpha_size <- alpha_fwd %>%
  filter(umi == "0") %>% 
    ggplot(aes(size, value, fill = size)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = size_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha diversity by size", x = "Size", y = "Value", fill = "Size")
p_alpha_size
ggsave(file.path(outdir, "alpha_by_size.png"), p_alpha_size, width=9, height=4.8, dpi=300)



# -----------------------------------------------------------------------------
# 4) PERMANOVA — unblocked + natman-blocked; natman R² readout across distances
# -----------------------------------------------------------------------------
nperm <- 999
perm_unblocked <- how(nperm = nperm)
perm_blocked   <- how(nperm = nperm); setBlocks(perm_blocked) <- FWD_meta$natman

# Use size (factor), control for log_reads
form_unblocked <- list(
  robAit  = adonis2(D_robAit ~ log_reads + umi + natman + fwd_state + decay_stage + size,
                    data = FWD_meta, by = "margin", permutations = perm_unblocked),
  hellEuc = adonis2(D_euclid_hell ~ log_reads + umi + natman + fwd_state + decay_stage + size,
                    data = FWD_meta, by = "margin", permutations = perm_unblocked),
  relBray = adonis2(D_bray_rel ~ log_reads + umi + natman + fwd_state + decay_stage + size,
                    data = FWD_meta, by = "margin", permutations = perm_unblocked),
  relBrayS= adonis2(D_bray_relS ~ log_reads + umi + natman + fwd_state + decay_stage + size,
                    data = FWD_meta, by = "margin", permutations = perm_unblocked),
  jaccard = adonis2(D_jacc ~ umi + natman + fwd_state + decay_stage + size,  # PA: omit log_reads
                    data = FWD_meta, by = "margin", permutations = perm_unblocked)
)

# Console readout: natman R² across distances
natman_R2 <- purrr::map_dfr(form_unblocked, ~{
  df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)
}, .id = "distance") %>%
  filter(term == "natman") %>%
  transmute(distance, R2, p = `Pr(>F)`) %>%
  arrange(desc(R2))
cat("\nUnblocked PERMANOVA — natman R² by distance (controls: log_reads where applicable):\n")
print(natman_R2)

# Blocked within-tree (unique microhabitat effects; size included)
perma_blocked_fwd <- adonis2(
  D_robAit ~ log_reads +umi +  fwd_state + decay_stage + size,
  data = FWD_meta, by = "margin", permutations = perm_blocked
)
cat("\nFWD blocked PERMANOVA (robust Aitchison):\n"); print(perma_blocked_fwd)

# Presence–absence robustness (blocked), keep same terms except log_reads
perma_pa <- adonis2(D_jacc ~ umi + fwd_state + decay_stage + size,
                    data = FWD_meta, by = "margin", permutations = perm_blocked)
cat("\nPresence–absence PERMANOVA (Jaccard, blocked by natman):\n"); print(perma_pa)

# -----------------------------------------------------------------------------
# 5) Partial dbRDA + variance partitioning (Condition(log_reads))
# -----------------------------------------------------------------------------
mod_all    <- capscale(otu_fwd ~ natman + fwd_state + decay_stage + size + Condition(log_reads + umi),
                       data = FWD_meta, distance = "robust.aitchison")
mod_tree   <- capscale(otu_fwd ~ natman + Condition(log_reads+umi),
                       data = FWD_meta, distance = "robust.aitchison")
mod_micro  <- capscale(otu_fwd ~ fwd_state + decay_stage + size + Condition(log_reads)+ Condition(umi),
                       data = FWD_meta, distance = "robust.aitchison")
mod_tree_c <- capscale(otu_fwd ~ natman + Condition(fwd_state + decay_stage + size + log_reads + umi),
                       data = FWD_meta, distance = "robust.aitchison")
mod_mic_c  <- capscale(otu_fwd ~ fwd_state + decay_stage + size + Condition(natman + log_reads + umi),
                       data = FWD_meta, distance = "robust.aitchison")

cat("\nVariance partitioning (R² unadjusted; total = mod_all$r.squared):\n")
R2_all   <- RsquareAdj(mod_all)$r.squared
R2_tree  <- RsquareAdj(mod_tree)$r.squared
R2_micro <- RsquareAdj(mod_micro)$r.squared
R2_treeU <- RsquareAdj(mod_tree_c)$r.squared
R2_micU  <- RsquareAdj(mod_mic_c)$r.squared
print(c(all=R2_all, tree=R2_tree, micro=R2_micro, tree_unique=R2_treeU, micro_unique=R2_micU))

cat("\nPartial dbRDA tests (blocked by natman):\n")
print(anova.cca(mod_all, permutations = perm_blocked, by = "terms"))
# Euler / fallback bars
make_euler_plot(
  R2_all   = R2_all,
  R2_tree  = R2_tree,
  R2_micro = R2_micro,
  R2_tree_u= R2_treeU,
  R2_mic_u = R2_micU,
  title_main = "Variance partitioning (FWD, robust Aitchison dbRDA)",
  title_sub  = "Tree vs Microhab complex (fwd_state + decay + size | log_reads)",
  outfile    = file.path(outdir, "euler_variance_partitioning_FWD.png")
)

# -----------------------------------------------------------------------------
# 6) Dispersion (PERMDISP) – report & visualize
# -----------------------------------------------------------------------------
bd_state <- betadisper(D_robAit, FWD_meta$fwd_state)
bd_ds    <- betadisper(D_robAit, FWD_meta$decay_stage)
bd_size  <- betadisper(D_robAit, FWD_meta$size)

cat("\nPERMDISP fwd_state (overall & pairwise):\n")
print(permutest(bd_state, permutations = how(nperm = nperm)))
print(permutest(bd_state, pairwise = TRUE, permutations = how(nperm = nperm)))

cat("\nPERMDISP decay_stage (overall & pairwise):\n")
print(permutest(bd_ds, permutations = how(nperm = nperm)))
print(permutest(bd_ds, pairwise = TRUE, permutations = how(nperm = nperm)))

cat("\nPERMDISP size (overall & pairwise):\n")
print(permutest(bd_size, permutations = how(nperm = nperm)))
print(permutest(bd_size, pairwise = TRUE, permutations = how(nperm = nperm)))

# Distances to centroid plots
dist_df_state <- tibble(sample = rownames(otu_fwd),
                        fwd_state = FWD_meta$fwd_state,
                        distance = as.numeric(bd_state$distances))
p_dist_state <- ggplot(dist_df_state, aes(fwd_state, distance, fill=fwd_state)) +
  geom_violin(trim=FALSE, alpha=0.5) + geom_boxplot(width=0.18, outlier.shape=NA) +
  geom_jitter(aes(color=fwd_state), width=0.08, height=0, size=1.6, alpha=0.6) +
  scale_fill_manual(values=fwd_colors) + scale_color_manual(values=fwd_colors) +
  labs(title="PERMDISP: distance to centroid (fwd_state)", y="Distance", x=NULL) +
  theme_minimal(base_size=12) + theme(legend.position="none")
p_dist_state
ggsave(file.path(outdir, "permdisp_fwd_state_distances.png"), p_dist_state, width=6.2, height=4.2, dpi=300)

dist_df_ds <- tibble(sample = rownames(otu_fwd),
                     decay_stage = FWD_meta$decay_stage,
                     distance = as.numeric(bd_ds$distances))
p_dist_ds <- ggplot(dist_df_ds, aes(decay_stage, distance, fill=decay_stage)) +
  geom_violin(trim=FALSE, alpha=0.5) + geom_boxplot(width=0.18, outlier.shape=NA) +
  geom_jitter(aes(color=decay_stage), width=0.08, height=0, size=1.6, alpha=0.6) +
  scale_fill_manual(values=DS_colors) + scale_color_manual(values=DS_colors) +
  labs(title="PERMDISP: distance to centroid (decay_stage)", y="Distance", x=NULL) +
  theme_minimal(base_size=12) + theme(legend.position="none")
p_dist_ds
ggsave(file.path(outdir, "permdisp_decay_stage_distances.png"), p_dist_ds, width=6.2, height=4.2, dpi=300)

# -----------------------------------------------------------------------------
# 7) Leave-one-tree-out (LOO) stability for within-tree effects
# -----------------------------------------------------------------------------
trees <- FWD_meta$natman %>% as.character() %>% unique() %>% setdiff(NA)
loo <- purrr::map_dfr(trees, function(t) {
  keep <- !is.na(FWD_meta$natman) & FWD_meta$natman != t
  if (sum(keep) < 5) return(NULL)
  Dk <- vegdist(otu_fwd[keep, , drop = FALSE], method = "robust.aitchison")
  permk <- how(nperm = 499); setBlocks(permk) <- droplevels(factor(FWD_meta$natman[keep]))
  m <- adonis2(Dk ~ log_reads + umi + fwd_state + decay_stage + size,
               data = FWD_meta[keep, ], by = "margin", permutations = permk)
  mm <- as.data.frame(m); mm$term <- rownames(mm)
  as_tibble(mm) %>%
    filter(term %in% c("fwd_state","decay_stage","size")) %>%
    transmute(left_out = t, term, R2, p = `Pr(>F)`)
})
if (nrow(loo) > 0) {
  cat("\nLOO stability — natman-blocked PERMANOVA terms:\n")
  print(
    loo %>% group_by(term) %>%
      summarise(median_R2 = median(R2, na.rm = TRUE),
                p05 = quantile(p, 0.05, na.rm = TRUE),
                p95 = quantile(p, 0.95, na.rm = TRUE),
                .groups = "drop")
  )
}
# TO OUR SURPRISE no effect of fallen vs attached? ( because only 3 natman with fFWD)

# -----------------------------------------------------------------------------
# 8) Indicator analyses (scaled to common library size) + BH FDR
# -----------------------------------------------------------------------------
otu_scaled <- rescale_to_target(otu_fwd)
# annotate SH_code|lowest_taxon
colnames(otu_scaled) <- tibble(sh_code = colnames(otu_scaled)) %>%
  left_join(select(tax, sh_code, lowest_taxon), by="sh_code") %>%
  transmute(name = if_else(!is.na(lowest_taxon), paste0(sh_code,"|",lowest_taxon), sh_code)) %>%
  pull(name)
otu_scaled <- otu_scaled[, colSums(otu_scaled) > 0, drop = FALSE]

group_fwd    <- droplevels(FWD_meta$fwd_state)
group_natman <- droplevels(FWD_meta$natman)
group_ds     <- droplevels(FWD_meta$decay_stage)

indic_rg <- indicspecies::multipatt(otu_scaled, group_fwd, func="r.g",
                                    control = how(nperm=999), duleg=TRUE)
tab_rg <- as_tibble(indic_rg$sign, rownames="feature") %>%
  mutate(p_adj = p.adjust(p.value, "BH")) %>% arrange(p_adj)
cat("\nIndicator (r.g) – top (BH<=0.10):\n")
print(tab_rg %>% filter(p_adj <= 0.10) %>% slice(1:10))

# Presence–absence (depth-robust companion)
otu_pa <- (otu_fwd > 0) * 1
colnames(otu_pa) <- tibble(sh_code = colnames(otu_pa)) %>%
  left_join(select(tax, sh_code, lowest_taxon), by="sh_code") %>%
  transmute(name = if_else(!is.na(lowest_taxon), paste0(sh_code,"|",lowest_taxon), sh_code)) %>%
  pull(name)
otu_pa <- otu_pa[, colSums(otu_pa) > 0, drop = FALSE]
indic_pa <- indicspecies::multipatt(otu_pa, group_fwd, func="IndVal.g",
                                    control = how(nperm=999), duleg=TRUE)
tab_pa <- as_tibble(indic_pa$sign, rownames="feature") %>%
  mutate(p_adj = p.adjust(p.value, "BH")) %>% arrange(p_adj)
cat("\nIndicator (IndVal.g, PA) – top (BH<=0.10):\n")
print(tab_pa %>% filter(p_adj <= 0.10) %>% slice(1:10))

# -----------------------------------------------------------------------------
# 9) Alpha richness GLMM — model selection + EMMs + IRR tables (improved phi + DHARMa)
# -----------------------------------------------------------------------------
df_alpha <- alpha_from_otu(otu_fwd) %>%
  left_join(FWD_meta %>% select(sample, natman, fwd_state, decay_stage, size, umi), by = "sample") %>%
  mutate(across(c(natman, fwd_state, decay_stage, size), droplevels)) %>%
  filter(is.finite(richness), reads > 0)

# Fit Poisson, compute φ (+CI,+p) and decide if NB is needed; then AIC-tiebreak
fit_pois <- glmmTMB(richness ~ fwd_state + decay_stage + size + log_reads + umi+ (1|natman),
                    data = df_alpha, family = poisson)

phi_info <- overdisp_phi(fit_pois)
use_nb   <- is.finite(phi_info$phi) && (phi_info$phi > 1.2)   # threshold you already used

if (use_nb) {
  fit_nb <- glmmTMB(richness ~ fwd_state + decay_stage + size + log_reads + umi+ (1|natman),
                    data = df_alpha, family = nbinom2,
                    control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))
  # Prefer the model with lower AIC; if close, NB is usually safer when φ>1
  best <- if (AIC(fit_nb) <= AIC(fit_pois) + 2) fit_nb else fit_pois
} else {
  best <- fit_pois
}

best_type <- if (identical(best$family$family, poisson()$family)) "Poisson (glmmTMB)" else "NB (glmmTMB)"

cat(sprintf(
  "\nAlpha richness GLMM — selected model: %s\nPoisson overdispersion: phi=%.3f (95%% CI %.3f–%.3f), df=%d, one-sided p_over=%.4f\n",
  best_type, phi_info$phi, phi_info$ci[1], phi_info$ci[2], phi_info$df, phi_info$p_over
))

# --- Global tests and R2 (unchanged) ---
lrt_tab <- drop1(best, test = "Chisq")
r2_tab  <- performance::r2(best)
cat("\nAlpha richness GLMM — Type-II LRT and R²:\n"); print(lrt_tab); print(r2_tab)

# --- IRR table (unchanged) ---
coefs <- broom.mixed::tidy(best, effects = "fixed", conf.int = TRUE, exponentiate = TRUE)
irr_tab <- coefs %>%
  transmute(term,
            IRR = estimate, CI_low = conf.low, CI_high = conf.high, p_value = p.value) %>%
  arrange(term)
cat("\nIncidence Rate Ratios (IRR) on richness:\n"); print(irr_tab)

# -----------------------------------------------------------------------------
# DHARMa diagnostics (headless): dispersion + zero-inflation + residual plots
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(DHARMa))
set.seed(1)
sim <- DHARMa::simulateResiduals(best, n = 2000)  # more sims = more stable tests
cat("\nDHARMa tests:\n")
print(DHARMa::testDispersion(sim))      # complements Pearson-φ
print(DHARMa::testZeroInflation(sim))   # useful for richness
#  standard DHARMa diagnostic plots 
plot(sim)  # QQ, residual vs pred, etc.

# -----------------------------------------------------------------------------
# EMMs (response scale) for each factor (unchanged)
# -----------------------------------------------------------------------------
emm_state <- emmeans(best, ~ fwd_state,  type = "response")
emm_ds    <- emmeans(best, ~ decay_stage, type = "response")
if (nlevels(df_alpha$size) > 1) emm_size <- emmeans(best, ~ size, type = "response") else emm_size <- NULL

plot_emm <- function(emm, xlab, fname) {
  df <- as.data.frame(emm)
  if ("asymp.LCL" %in% names(df)) df <- dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  x_var <- names(df)[1]
  p <- ggplot(df, aes(x = !!sym(x_var), y = response, ymin = lower.CL, ymax = upper.CL)) +
    geom_point(size = 3) +
    geom_errorbar(width = 0.2) +
    labs(x = xlab, y = "Estimated richness",
         title = paste("Richness ~", xlab, "(EMM,", best_type, ")")) +
    theme_classic(base_size = 11)
  ggsave(file.path(outdir, fname), p, width = 5, height = 3.4, dpi = 300)
  print(p)
}
plot_emm(emm_state, "FWD state",  "alpha_emm_fwd_state.png")
plot_emm(emm_ds,    "Decay stage","alpha_emm_decay_stage.png")
if (!is.null(emm_size)) plot_emm(emm_size, "Size", "alpha_emm_size.png")

# Compact publication tables (unchanged)
publist <- list(
  model = summary(best),
  lrt   = lrt_tab,
  r2    = r2_tab,
  irr   = irr_tab,
  phi   = phi_info
)
publist
saveRDS(publist, file.path(outdir, "alpha_glmm_publication_tables.rds"))
cat("\nSaved compact publication tables to:", file.path(outdir, "alpha_glmm_publication_tables.rds"), "\n")

# Pairwise contrasts (unchanged)
emm_state_ds <- emmeans(best, ~ fwd_state + decay_stage, type = "response")
contr_pw <- contrast(emm_state_ds, method = "pairwise", type = "response", adjust = "BH")
cat("\nPairwise contrasts (EMM, response scale, BH):\n")
print(as.data.frame(summary(contr_pw)))
# for txt: “We fit Poisson and NB GLMMs (glmmTMB), selected the final model by Poisson overdispersion (Pearson φ with 95% CI and χ² test) and AIC, and verified residual dispersion and zero-inflation with DHARMa.”

# -----------------------------------------------------------------------------
# 10) Matched natman subset (trees with both ATTACHED & FALLEN)
# -----------------------------------------------------------------------------
tab_match <- with(FWD_meta, table(natman, fwd_state))
matched_ids <- rownames(tab_match)[rowSums(tab_match > 0) == 2]
keep <- FWD_meta$natman %in% matched_ids
if (sum(keep) >= 4) {
  D_rob_mat <- vegdist(otu_fwd[keep, , drop = FALSE], method = "robust.aitchison")
  perm_mat  <- how(nperm = 999); setBlocks(perm_mat) <- FWD_meta$natman[keep]
  cat("\nMatched-natman (both states present) — blocked PERMANOVA:\n")
  print(
    adonis2(D_rob_mat ~ log_reads + umi+ decay_stage + size,
            data = FWD_meta[keep, ], by = "margin", permutations = perm_mat)
  )
} else {
  cat("\nMatched-natman subset too small; skipping.\n")
}

# -----------------------------------------------------------------------------
# 11) FALLEN-only subset
# -----------------------------------------------------------------------------
sel <- FWD_meta$fwd_state == "FALLEN"
if (sum(sel) >= 4) {
  otu_fall  <- otu_fwd[sel, , drop = FALSE]
  meta_fall <- droplevels(FWD_meta[sel, ])
  cat("\nFALLEN-only: n =", nrow(meta_fall), "\n")
  D_rob_fall <- vegdist(otu_fall, method = "robust.aitchison")
  D_jac_fall <- vegdist((otu_fall > 0) * 1, method = "jaccard", binary = TRUE)
  
  # PERMANOVA (unblocked): natman as a predictor (no blocking)
  cat("\nFALLEN-only PERMANOVA (robust Aitchison):\n")
  print(adonis2(D_rob_fall ~ log_reads + decay_stage + size + rep + natman,
                data = meta_fall, by = "margin", permutations = how(nperm = 999)))
  cat("\nFALLEN-only PERMANOVA (Jaccard PA):\n")
  print(adonis2(D_jac_fall ~ decay_stage + size + rep + natman,
                data = meta_fall, by = "margin", permutations = how(nperm = 999)))
  
  # Dispersion checks
  cat("\nFALLEN-only PERMDISP checks:\n")
  print(anova(betadisper(D_rob_fall, meta_fall$natman)))
  print(anova(betadisper(D_rob_fall, meta_fall$decay_stage)))
  print(anova(betadisper(D_rob_fall, meta_fall$size)))
  if (!all(is.na(meta_fall$rep))) print(anova(betadisper(D_rob_fall, meta_fall$rep)))
  
  # Partial dbRDA (condition on log_reads and rep)
  mod_all_fall <- capscale(otu_fall ~ natman + size + decay_stage + Condition(log_reads + rep),
                           data = meta_fall, distance = "robust.aitchison")
  mod_nat_fall <- capscale(otu_fall ~ natman + Condition(log_reads + rep),
                           data = meta_fall, distance = "robust.aitchison")
  mod_mic_fall <- capscale(otu_fall ~ size + decay_stage + Condition(log_reads + rep),
                           data = meta_fall, distance = "robust.aitchison")
  cat("\nFALLEN-only variance fractions (unadjusted):\n")
  print(c(
    all = RsquareAdj(mod_all_fall)$r.squared,
    natman_unique   = RsquareAdj(mod_all_fall)$r.squared - RsquareAdj(mod_mic_fall)$r.squared,
    microhab_unique = RsquareAdj(mod_all_fall)$r.squared - RsquareAdj(mod_nat_fall)$r.squared
  ))
  cat("\nFALLEN-only mod_all term tests:\n")
  print(anova.cca(mod_all_fall, permutations = how(nperm = 999), by = "margin"))
  
  # LOOCV 1-NN on Jaccard
  D <- as.matrix(D_jac_fall)
  labs <- as.character(meta_fall$natman)
  pred <- sapply(seq_len(nrow(D)), function(i) { others <- setdiff(seq_len(nrow(D)), i)
  nn <- others[which.min(D[i, others])]; labs[nn] })
  acc <- mean(pred == labs)
  cat("\nFALLEN-only LOOCV 1-NN accuracy (Jaccard):", sprintf("%.1f%%", 100*acc), "\n")
  base_acc <- 1 / nlevels(droplevels(meta_fall$natman))
  cat(sprintf(" (chance baseline ≈ %.1f%%)", 100*base_acc), "\n")
  
}
# more variance partitioning
R_all <- RsquareAdj(mod_all_fall)$r.squared
R_nat <- RsquareAdj(mod_nat_fall)$r.squared
R_mic <- RsquareAdj(mod_mic_fall)$r.squared
# Unique fractions
nat_unique   <- R_all - R_mic
mic_unique   <- R_all - R_nat
# Shared fraction = what's in R_all but counted in both predictors
shared       <- R_all - nat_unique - mic_unique
# Residual (unexplained)
residual     <- 1 - R_all

fractions <- c(
  nat_unique   = nat_unique,
  mic_unique   = mic_unique,
  shared       = shared,
  residual     = residual,
  total_expl   = R_all
)
fractions[fractions < 0] <- 0   # truncate small negatives to 0
fractions <- fractions / sum(fractions)  # renormalize if you want them to add to 1
round(fractions, 3)

vennplot <- varpart(otu_fall, ~natman, ~size + decay_stage, data = meta_fall)
plot(vennplot, bg = c("skyblue", "orange"))

## ============================================================
## FALLEN FWD analyses + matched SOIL "signature" checks
## ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(betapart)
  library(indicspecies)
  library(permute)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

# I/O dirs ---------------------------------------------------------------
plots_dir  <- "plots/FALLEN"
tables_dir <- "tables/FALLEN"
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Helpers ---------------------------------------------------------------
tss <- function(M, target_sum = NULL) {
  # If target_sum is NULL: return compositional proportions (sums to 1).
  # If target_sum is numeric: scale each sample to that sum.
  rs <- pmax(rowSums(M), 1)
  P  <- sweep(M, 1, rs, "/")
  if (is.null(target_sum)) P else P * target_sum
}

# Blocked permutation within grouping factor levels (size-1 levels fixed)
permute_within <- function(g) {
  idx <- seq_along(g)
  unlist(lapply(split(idx, as.factor(g)),
                function(v) if (length(v) > 1) sample(v) else v),
         use.names = FALSE)
}

# =======================================================================
# 0) FALLEN FWD subset & alignment
# =======================================================================
meta_fall <- FWD_meta %>%
  filter(fwd_state == "FALLEN") %>%
  droplevels()

stopifnot("sample" %in% names(meta_fall))
stopifnot(all(meta_fall$sample %in% rownames(otu_fwd)))

otu_fall <- otu_fwd[meta_fall$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fall), meta_fall$sample))

# Optional: drop zero-columns within this subset
otu_fall <- otu_fall[, colSums(otu_fall) > 0, drop = FALSE]

# Presence–absence for betapart
otu_fall_pa <- (otu_fall > 0) * 1

# Decay stage as factor (respect your ordinal if present)
meta_fall$decay_stage <- droplevels(factor(meta_fall$decay_stage))

# =======================================================================
# 1) TURNOVER vs NESTEDNESS (betapart; Sørensen family, PA-based)
# =======================================================================
bp <- betapart::beta.pair(otu_fall_pa, index.family = "sorensen")
# bp$beta.sim (turnover), bp$beta.sne (nestedness), bp$beta.sor (total)

# Overall (all samples)
beta_all <- betapart::beta.multi(otu_fall_pa, index.family = "sorensen")
cat("\nOverall Sørensen partitions (FALLEN):\n"); print(beta_all)

# Within/between decay stage means
DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

nS <- nrow(DM_sim)
g  <- meta_fall$decay_stage

pairs_df <- expand.grid(i = seq_len(nS), j = seq_len(nS)) %>%
  filter(i < j) %>%
  mutate(stage_i = g[i], stage_j = g[j],
         type    = if_else(stage_i == stage_j, "within", "between"),
         sim     = DM_sim[cbind(i, j)],
         sne     = DM_sne[cbind(i, j)],
         sor     = DM_sor[cbind(i, j)])

within_tab <- pairs_df %>%
  filter(type == "within") %>%
  group_by(stage_i) %>%
  summarise(mean_SIM = mean(sim, na.rm = TRUE),
            mean_SNE = mean(sne, na.rm = TRUE),
            mean_SOR = mean(sor, na.rm = TRUE),
            n_pairs  = dplyr::n(), .groups = "drop") %>%
  rename(decay_stage = stage_i)

between_tab <- pairs_df %>%
  filter(type == "between") %>%
  mutate(pair = paste(pmin(as.character(stage_i), as.character(stage_j)),
                      pmax(as.character(stage_i), as.character(stage_j)), sep = "–")) %>%
  group_by(pair) %>%
  summarise(mean_SIM = mean(sim, na.rm = TRUE),
            mean_SNE = mean(sne, na.rm = TRUE),
            mean_SOR = mean(sor, na.rm = TRUE),
            n_pairs  = dplyr::n(), .groups = "drop")

cat("\nWithin-stage means (FALLEN):\n");  print(within_tab)
cat("\nBetween-stage means (FALLEN):\n"); print(between_tab)

turnover_ratio_overall <- as.numeric(beta_all$beta.SIM / beta_all$beta.SNE)
cat(sprintf("\nTurnover-to-nestedness ratio (overall): %.2f\n", turnover_ratio_overall))

# Save tables
readr::write_tsv(within_tab,  file.path(tables_dir, "betapart_within_decay.tsv"))
readr::write_tsv(between_tab, file.path(tables_dir, "betapart_between_decay.tsv"))

# =======================================================================
# 2) INDICATOR SPECIES across decay stages (IndVal.g; abundance-aware)
# =======================================================================
# Use TSS-normalized counts for summaries; for multipatt you can use raw or TSS.
# Here we keep raw counts for IndVal.g (scale-invariant) and report TSS means.
otu_fall_tss <- tss(otu_fall, target_sum = 3000)

set.seed(42)
res_ind <- indicspecies::multipatt(
  otu_fall, meta_fall$decay_stage,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)

sign_tbl <- as.data.frame(res_ind$sign)
if (!nrow(sign_tbl)) {
  cat("\nNo significant decay-stage indicators at alpha = 0.05.\n")
} else {
  sign_tbl$sh_code <- rownames(sign_tbl)
  gcols <- grep("^s\\.", names(sign_tbl), value = TRUE)
  is_single    <- rowSums(sign_tbl[, gcols, drop = FALSE] > 0) == 1
  sign_single  <- sign_tbl[is_single, , drop = FALSE]
  idx_mat      <- sign_single[, gcols, drop = FALSE] > 0
  which_col    <- apply(idx_mat, 1, function(z) which(z)[1])
  stage_lab    <- sub("^s\\.", "", gcols[which_col])
  
  ind_df <- tibble(
    sh_code     = sign_single$sh_code,
    decay_stage = factor(stage_lab, levels = levels(meta_fall$decay_stage)),
    stat        = sign_single$stat,
    p_value     = sign_single$p.value
  ) %>%
    arrange(p_value, desc(stat))
  
  # Prevalence (PA) and mean RA (TSS) by stage
  PA <- (otu_fall > 0) * 1
  prev_stage <- map_df(levels(meta_fall$decay_stage), function(st) {
    rows <- which(meta_fall$decay_stage == st)
    tibble(sh_code = colnames(PA),
           decay_stage = st,
           prevalence = colMeans(PA[rows, , drop = FALSE]))
  })
  
  mra_stage <- map_df(levels(meta_fall$decay_stage), function(st) {
    rows <- which(meta_fall$decay_stage == st)
    tibble(sh_code = colnames(otu_fall_tss),
           decay_stage = st,
           mean_rel_abund = colMeans(otu_fall_tss[rows, , drop = FALSE]))
  })
  
  # Attach taxonomy (if available)
  tax_keep <- intersect(
    c("sh_code","kingdom","phylum","class","order","family","genus","species","lowest_taxon","sh_link"),
    names(tax)
  )
  tax_slim <- if (length(tax_keep)) distinct(tax[, tax_keep, drop = FALSE]) else tibble(sh_code = character())
  
  ind_annot <- ind_df %>%
    left_join(prev_stage, by = c("sh_code","decay_stage")) %>%
    left_join(mra_stage, by = c("sh_code","decay_stage")) %>%
    left_join(tax_slim,   by = "sh_code") %>%
    arrange(decay_stage, p_value, desc(stat)) %>%
    filter(p_value < 0.05) %>%
    select(decay_stage, p_value, stat, lowest_taxon, sh_code, prevalence, mean_rel_abund)
  
  cat("\nTop decay-stage indicators (single-stage; alpha=0.05):\n")
  print(head(ind_annot, 20), n = 20)
  
  readr::write_tsv(ind_annot, file.path(tables_dir, "indicators_decay_single_stage.tsv"))
}

# =======================================================================
# 3) Partial dbRDA: decay + size | (log_reads + natman); robust Aitchison
# =======================================================================
covars <- c("log_reads", "natman")
missing_cov <- setdiff(covars, names(meta_fall))
if (length(missing_cov)) stop("Missing covariates in meta_fall: ", paste(missing_cov, collapse = ", "))

mod_decay <- capscale(otu_fall ~ decay_stage + size + Condition(log_reads + natman),
                      data = meta_fall, distance = "robust.aitchison")

print(anova.cca(mod_decay, permutations = permute::how(nperm = 999), by = "terms"))
R2s <- RsquareAdj(mod_decay)

eig_con <- as.numeric(vegan::eigenvals(mod_decay, constrained = TRUE))
tot_var <- mod_decay$tot.chi
cap1_pct_total <- if (length(eig_con) >= 1 && is.finite(eig_con[1])) 100 * eig_con[1] / tot_var else NA_real_
cap2_pct_total <- if (length(eig_con) >= 2 && is.finite(eig_con[2])) 100 * eig_con[2] / tot_var else NA_real_
constrained_var <- sum(eig_con)
cap1_pct_constr <- if (constrained_var > 0) 100 * eig_con[1] / constrained_var else NA_real_
cap2_pct_constr <- if (constrained_var > 0 && length(eig_con) >= 2) 100 * eig_con[2] / constrained_var else NA_real_

cat(sprintf("\nAxis contributions:\n  CAP1 = %.2f%% of TOTAL (%.2f%% of constrained)\n  CAP2 = %.2f%% of TOTAL (%.2f%% of constrained)\n  Constrained total (unadj. R2) = %.2f%%\n  Adjusted R2 = %.2f%%\n",
            cap1_pct_total, cap1_pct_constr,
            cap2_pct_total, cap2_pct_constr,
            100 * constrained_var / tot_var,
            100 * R2s$adj.r.squared))

sc_sites <- scores(mod_decay, display = "sites", choices = 1:2)
plot_df <- as_tibble(sc_sites, rownames = "sample") %>%
  left_join(meta_fall %>% select(sample, decay_stage, size), by = "sample")

decay_cols <- c(EARLY = "#B58900", LATE  = "#6A51A3")
x_lab <- sprintf("CAP1 (%.1f%% of total variance)", cap1_pct_total)
y_lab <- sprintf("CAP2 (%.1f%% of total variance)", cap2_pct_total)

gg1 <- ggplot(plot_df, aes(CAP1, CAP2, color = decay_stage, shape = size)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(size = 2.6, alpha = 0.95) +
  stat_ellipse(aes(group = decay_stage), type = "norm", level = 0.68,
               linewidth = 0.35, alpha = 0.35) +
  scale_color_manual(values = decay_cols, drop = FALSE, name = "Decay stage") +
  coord_equal() +
  labs(title = "FALLEN FWD — partial dbRDA",
       subtitle = "Constrained by decay stage + size; conditioned on log_reads + natman",
       x = x_lab, y = y_lab, shape = "Size class") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
print(gg1)
ggsave(file.path(plots_dir, "fallen_dbRDA_decay_size.png"), gg1, width = 6.5, height = 5.2, dpi = 300)

cat("\nINTERPRETATION NOTES:\n",
    " - If beta.SIM >> beta.SNE overall/among stages, decay signal reflects turnover (replacement).\n",
    " - Low within-stage SOR suggests relative within-stage homogeneity; high values suggest heterogeneity.\n",
    " - Indicators highlight guild shifts across decay (wood-decayers → litter/soil-associated taxa, etc.).\n",
    " - Partial dbRDA significance = decay/size gradient after controlling site & depth proxies.\n", sep = "")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(betapart)
  library(indicspecies)
  library(permute)
  library(ggplot2)
  library(ggrepel)
})

set.seed(42)

# ------------------------------------------------------------------------------
# I/O dirs + helpers
# ------------------------------------------------------------------------------
plots_dir  <- "plots/FALLEN"
tables_dir <- "tables/FALLEN"
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Total-sum scaling helper (proportions if target_sum=NULL; else scaled counts)
tss <- function(M, target_sum = NULL) {
  rs <- pmax(rowSums(M), 1)
  P  <- sweep(M, 1, rs, "/")
  if (is.null(target_sum)) P else P * target_sum
}

# Blocked permutation within grouping factor levels (size-1 levels fixed)
permute_within <- function(g) {
  idx <- seq_along(g)
  unlist(lapply(split(idx, as.factor(g)),
                function(v) if (length(v) > 1) sample(v) else v),
         use.names = FALSE)
}

# =======================================================================
# 0) FALLEN FWD subset & alignment
# =======================================================================
meta_fall <- FWD_meta %>%
  filter(fwd_state == "FALLEN") %>%
  droplevels()

stopifnot("sample" %in% names(meta_fall))
stopifnot(all(meta_fall$sample %in% rownames(otu_fwd)))

otu_fall <- otu_fwd[meta_fall$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_fall), meta_fall$sample))

# Optional: drop zero-columns within this subset
otu_fall <- otu_fall[, colSums(otu_fall) > 0, drop = FALSE]

# Presence–absence for betapart
otu_fall_pa <- (otu_fall > 0) * 1

# Decay stage as factor (respect your ordinal if present)
meta_fall$decay_stage <- droplevels(factor(meta_fall$decay_stage))

# =======================================================================
# 1) TURNOVER vs NESTEDNESS (betapart; Sørensen family, PA-based)
# =======================================================================
bp <- betapart::beta.pair(otu_fall_pa, index.family = "sorensen")
# bp$beta.sim (turnover), bp$beta.sne (nestedness), bp$beta.sor (total)

# Overall (all samples)
beta_all <- betapart::beta.multi(otu_fall_pa, index.family = "sorensen")
cat("\nOverall Sørensen partitions (FALLEN):\n"); print(beta_all)

# Within/between decay stage means
DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

nS <- nrow(DM_sim)
g  <- meta_fall$decay_stage

pairs_df <- expand.grid(i = seq_len(nS), j = seq_len(nS)) %>%
  filter(i < j) %>%
  mutate(stage_i = g[i], stage_j = g[j],
         type    = if_else(stage_i == stage_j, "within", "between"),
         sim     = DM_sim[cbind(i, j)],
         sne     = DM_sne[cbind(i, j)],
         sor     = DM_sor[cbind(i, j)])

within_tab <- pairs_df %>%
  filter(type == "within") %>%
  group_by(stage_i) %>%
  summarise(mean_SIM = mean(sim, na.rm = TRUE),
            mean_SNE = mean(sne, na.rm = TRUE),
            mean_SOR = mean(sor, na.rm = TRUE),
            n_pairs  = dplyr::n(), .groups = "drop") %>%
  rename(decay_stage = stage_i)

between_tab <- pairs_df %>%
  filter(type == "between") %>%
  mutate(pair = paste(pmin(as.character(stage_i), as.character(stage_j)),
                      pmax(as.character(stage_i), as.character(stage_j)), sep = "–")) %>%
  group_by(pair) %>%
  summarise(mean_SIM = mean(sim, na.rm = TRUE),
            mean_SNE = mean(sne, na.rm = TRUE),
            mean_SOR = mean(sor, na.rm = TRUE),
            n_pairs  = dplyr::n(), .groups = "drop")

cat("\nWithin-stage means (FALLEN):\n");  print(within_tab)
cat("\nBetween-stage means (FALLEN):\n"); print(between_tab)

turnover_ratio_overall <- as.numeric(beta_all$beta.SIM / beta_all$beta.SNE)
cat(sprintf("\nTurnover-to-nestedness ratio (overall): %.2f\n", turnover_ratio_overall))

# Save tables
readr::write_tsv(within_tab,  file.path(tables_dir, "betapart_within_decay.tsv"))
readr::write_tsv(between_tab, file.path(tables_dir, "betapart_between_decay.tsv"))

# =======================================================================
# 2) INDICATOR SPECIES across decay stages (IndVal.g; abundance-aware)
# =======================================================================
# Use TSS-normalized counts for summaries; for multipatt you can use raw or TSS.
# Here we keep raw counts for IndVal.g (scale-invariant) and report TSS means.
otu_fall_tss <- tss(otu_fall, target_sum = 3000)

set.seed(42)
res_ind <- indicspecies::multipatt(
  otu_fall, meta_fall$decay_stage,
  func = "IndVal.g",
  control = permute::how(nperm = 999)
)

sign_tbl <- as.data.frame(res_ind$sign)
if (!nrow(sign_tbl)) {
  cat("\nNo significant decay-stage indicators at alpha = 0.05.\n")
} else {
  sign_tbl$sh_code <- rownames(sign_tbl)
  gcols <- grep("^s\\.", names(sign_tbl), value = TRUE)
  is_single    <- rowSums(sign_tbl[, gcols, drop = FALSE] > 0) == 1
  sign_single  <- sign_tbl[is_single, , drop = FALSE]
  idx_mat      <- sign_single[, gcols, drop = FALSE] > 0
  which_col    <- apply(idx_mat, 1, function(z) which(z)[1])
  stage_lab    <- sub("^s\\.", "", gcols[which_col])
  
  ind_df <- tibble(
    sh_code     = sign_single$sh_code,
    decay_stage = factor(stage_lab, levels = levels(meta_fall$decay_stage)),
    stat        = sign_single$stat,
    p_value     = sign_single$p.value
  ) %>%
    arrange(p_value, desc(stat))
  
  # Prevalence (PA) and mean RA (TSS) by stage
  PA <- (otu_fall > 0) * 1
  prev_stage <- map_df(levels(meta_fall$decay_stage), function(st) {
    rows <- which(meta_fall$decay_stage == st)
    tibble(sh_code = colnames(PA),
           decay_stage = st,
           prevalence = colMeans(PA[rows, , drop = FALSE]))
  })
  
  mra_stage <- map_df(levels(meta_fall$decay_stage), function(st) {
    rows <- which(meta_fall$decay_stage == st)
    tibble(sh_code = colnames(otu_fall_tss),
           decay_stage = st,
           mean_rel_abund = colMeans(otu_fall_tss[rows, , drop = FALSE]))
  })
  
  # Attach taxonomy (if available)
  tax_keep <- intersect(
    c("sh_code","kingdom","phylum","class","order","family","genus","species","lowest_taxon","sh_link"),
    names(tax)
  )
  tax_slim <- if (length(tax_keep)) distinct(tax[, tax_keep, drop = FALSE]) else tibble(sh_code = character())
  
  ind_annot <- ind_df %>%
    left_join(prev_stage, by = c("sh_code","decay_stage")) %>%
    left_join(mra_stage, by = c("sh_code","decay_stage")) %>%
    left_join(tax_slim,   by = "sh_code") %>%
    arrange(decay_stage, p_value, desc(stat)) %>%
    filter(p_value < 0.05) %>%
    select(decay_stage, p_value, stat, lowest_taxon, sh_code, prevalence, mean_rel_abund)
  
  cat("\nTop decay-stage indicators (single-stage; alpha=0.05):\n")
  print(head(ind_annot, 20), n = 20)
  
  readr::write_tsv(ind_annot, file.path(tables_dir, "indicators_decay_single_stage.tsv"))
}

# =======================================================================
# 3) Partial dbRDA: decay + size | (log_reads + natman); robust Aitchison
# =======================================================================
covars <- c("log_reads", "natman")
missing_cov <- setdiff(covars, names(meta_fall))
if (length(missing_cov)) stop("Missing covariates in meta_fall: ", paste(missing_cov, collapse = ", "))

mod_decay <- capscale(otu_fall ~ decay_stage + size + Condition(log_reads + natman),
                      data = meta_fall, distance = "robust.aitchison")

print(anova.cca(mod_decay, permutations = permute::how(nperm = 999), by = "terms"))
R2s <- RsquareAdj(mod_decay)

eig_con <- as.numeric(vegan::eigenvals(mod_decay, constrained = TRUE))
tot_var <- mod_decay$tot.chi
cap1_pct_total <- if (length(eig_con) >= 1 && is.finite(eig_con[1])) 100 * eig_con[1] / tot_var else NA_real_
cap2_pct_total <- if (length(eig_con) >= 2 && is.finite(eig_con[2])) 100 * eig_con[2] / tot_var else NA_real_
constrained_var <- sum(eig_con)
cap1_pct_constr <- if (constrained_var > 0) 100 * eig_con[1] / constrained_var else NA_real_
cap2_pct_constr <- if (constrained_var > 0 && length(eig_con) >= 2) 100 * eig_con[2] / constrained_var else NA_real_

cat(sprintf("\nAxis contributions:\n  CAP1 = %.2f%% of TOTAL (%.2f%% of constrained)\n  CAP2 = %.2f%% of TOTAL (%.2f%% of constrained)\n  Constrained total (unadj. R2) = %.2f%%\n  Adjusted R2 = %.2f%%\n",
            cap1_pct_total, cap1_pct_constr,
            cap2_pct_total, cap2_pct_constr,
            100 * constrained_var / tot_var,
            100 * R2s$adj.r.squared))

sc_sites <- scores(mod_decay, display = "sites", choices = 1:2)
plot_df <- as_tibble(sc_sites, rownames = "sample") %>%
  left_join(meta_fall %>% select(sample, decay_stage, size), by = "sample")

decay_cols <- c(EARLY = "#B58900", LATE  = "#6A51A3")
x_lab <- sprintf("CAP1 (%.1f%% of total variance)", cap1_pct_total)
y_lab <- sprintf("CAP2 (%.1f%% of total variance)", cap2_pct_total)

gg1 <- ggplot(plot_df, aes(CAP1, CAP2, color = decay_stage, shape = size)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey70") +
  geom_point(size = 2.6, alpha = 0.95) +
  stat_ellipse(aes(group = decay_stage), type = "norm", level = 0.68,
               linewidth = 0.35, alpha = 0.35) +
  scale_color_manual(values = decay_cols, drop = FALSE, name = "Decay stage") +
  coord_equal() +
  labs(title = "FALLEN FWD — partial dbRDA",
       subtitle = "Constrained by decay stage + size; conditioned on log_reads + natman",
       x = x_lab, y = y_lab, shape = "Size class") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
print(gg1)
ggsave(file.path(plots_dir, "fallen_dbRDA_decay_size.png"), gg1, width = 6.5, height = 5.2, dpi = 300)

cat("\nINTERPRETATION NOTES:\n",
    " - If beta.SIM >> beta.SNE overall/among stages, decay signal reflects turnover (replacement).\n",
    " - Low within-stage SOR suggests relative within-stage homogeneity; high values suggest heterogeneity.\n",
    " - Indicators highlight guild shifts across decay (wood-decayers → litter/soil-associated taxa, etc.).\n",
    " - Partial dbRDA significance = decay/size gradient after controlling site & depth proxies.\n", sep = "")


