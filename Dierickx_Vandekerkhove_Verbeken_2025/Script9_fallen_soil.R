# =============================================================================
# Script 9 — FALLEN FWD + matched SOIL (pair-based tests & visuals)
#
# Scopes:
#   - FALLEN FWD subset only: community patterns across decay_stage & size
#   - Matched FWD–SOIL pairs linked via META1$soil_core (FRUITBODY excluded)
#   - Beta diversity partitioning (betapart, Sørensen): turnover vs nestedness;
#     within- and between-stage summaries.
#   - Indicator taxa by decay stage (IndVal.g) + prevalence (PA) and mean RA (TSS);
#   - Partial dbRDA (robust Aitchison): decay + size | (log_reads + natman),
#     with axis contribution readout.
#   - Pairing diagnostics for FWD↔SOIL (1:1 mapping via soil_core).
#   - Soil “signature” tests on PA Jaccard (blocked perms):
#       * Global mean paired Jaccard
#       * Per-pair Monte Carlo p-values + BH correction
#       * Hypergeometric overlap enrichment per pair
#       * Top-1 nearest-soil accuracy under blocked null
#   - Obs vs null per-pair plots (faceted by natman) with ordered y-axis labels:
#       * Jaccard similarity (PA)      — RIGHT-tail (observed > null)
#       * Bray–Curtis distance (TSS)   — LEFT-tail  (observed < null)
#       * Robust Aitchison distance    — LEFT-tail  (observed < null)
#     Labels show “decay · size”, ordered EARLY→LATE then FINE→VERY_FINE;
#     asterisks mark metric-specific significant pairs (BH q < 0.05).
#   - Taxon-level concordance across pairs with site-blocked null, effect sizes
#     (delta, z), BH-adjusted q, plus a simple McNemar co-occurrence check.
# =============================================================================
source("scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(betapart)
  library(indicspecies)
  library(permute)
  library(ggplot2)
  library(ggrepel)
  # forcats used via :: calls
})

set.seed(42)

# =============================================================================
# Setup: IO dirs + helpers
# =============================================================================
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

# =============================================================================
# 0) Choose OTU source (filtered preferred) + align FWD_meta / otu_fwd
# =============================================================================
otu_source <- otu_matrix_filt

# Keep only valid FWD states (attached/fallen); fail fast otherwise
fwd <- dplyr::filter(META1, dw_type == "FWD")
bad_levels <- setdiff(unique(fwd$position), c("ATTACHED","FALLEN"))
stopifnot(length(bad_levels) == 0)

FWD_meta <- META1 %>%
  filter(dw_type == "FWD") %>%
  mutate(
    decay_stage = factor(decay_stage, levels = c("EARLY","AVERAGE","LATE")),
    ds_at_drill = factor(ds_at_drill, levels = c("0","1","2","3","4","5")),
    size        = droplevels(factor(size)),
    fwd_state   = as.character(position)
  )

# Align OTU and add depth proxies
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
    "| fwd_state:", paste(levels(FWD_meta$fwd_state), collapse = ", "),
    "| decay:",    paste(levels(FWD_meta$decay_stage), collapse = ", "),
    "| size:",     paste(levels(FWD_meta$size),        collapse = ", "), "\n")

# Optional overview distances
D_robAit      <- vegdist(otu_fwd, method = "robust.aitchison")
comm_hell     <- decostand(otu_fwd, method = "hellinger")
D_euclid_hell <- dist(comm_hell, method = "euclidean")
comm_rel      <- decostand(otu_fwd, method = "total")
D_bray_rel    <- vegdist(comm_rel, method = "bray")
D_bray_relS   <- vegdist(sqrt(comm_rel), method = "bray")
D_jacc        <- vegdist((otu_fwd > 0) * 1, method = "jaccard", binary = TRUE)

# =============================================================================
# 1) FALLEN-only subset: betapart, indicators, partial dbRDA
# =============================================================================
sel <- FWD_meta$fwd_state == "FALLEN"
stopifnot(sum(sel) >= 4)
otu_fall  <- otu_fwd[sel, , drop = FALSE]
meta_fall <- droplevels(FWD_meta[sel, ])
cat("\nFALLEN-only: n =", nrow(meta_fall), "\n")

# ---- betapart (turnover vs nestedness; Sørensen on PA) -----------------
otu_fall_pa <- (otu_fall > 0) * 1
bp       <- betapart::beta.pair(otu_fall_pa, index.family = "sorensen")
beta_all <- betapart::beta.multi(otu_fall_pa, index.family = "sorensen")
cat("\nOverall Sørensen partitions (FALLEN):\n"); print(beta_all)

DM_sim <- as.matrix(bp$beta.sim)
DM_sne <- as.matrix(bp$beta.sne)
DM_sor <- as.matrix(bp$beta.sor)

nS <- nrow(DM_sim); g <- meta_fall$decay_stage
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

readr::write_tsv(within_tab,  file.path(tables_dir, "betapart_within_decay.tsv"))
readr::write_tsv(between_tab, file.path(tables_dir, "betapart_between_decay.tsv"))

# ---- Indicator species across decay (IndVal.g) -------------------------
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
  is_single   <- rowSums(sign_tbl[, gcols, drop = FALSE] > 0) == 1
  sign_single <- sign_tbl[is_single, , drop = FALSE]
  idx_mat     <- sign_single[, gcols, drop = FALSE] > 0
  which_col   <- apply(idx_mat, 1, function(z) which(z)[1])
  stage_lab   <- sub("^s\\.", "", gcols[which_col])
  
  ind_df <- tibble(
    sh_code     = sign_single$sh_code,
    decay_stage = factor(stage_lab, levels = levels(meta_fall$decay_stage)),
    stat        = sign_single$stat,
    p_value     = sign_single$p.value
  ) %>% arrange(p_value, desc(stat))
  
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

# ---- Partial dbRDA (robust Aitchison): decay + size | (log_reads + natman)
covars <- c("log_reads", "natman")
missing_cov <- setdiff(covars, names(meta_fall))
if (length(missing_cov)) stop("Missing covariates in meta_fall: ", paste(missing_cov, collapse = ", "))

# Joint model (both microhabitat terms together)
mod_decay <- capscale(
  otu_fall ~ decay_stage + size + Condition(log_reads + natman),
  data = meta_fall, distance = "robust.aitchison"
)
print(anova.cca(mod_decay, permutations = permute::how(nperm = 999), by = "terms"))
R2s <- RsquareAdj(mod_decay)

eig_con <- as.numeric(vegan::eigenvals(mod_decay, model = "constrained"))
tot_var <- mod_decay$tot.chi
cap1_pct_total  <- if (length(eig_con) >= 1 && is.finite(eig_con[1])) 100 * eig_con[1] / tot_var else NA_real_
cap2_pct_total  <- if (length(eig_con) >= 2 && is.finite(eig_con[2])) 100 * eig_con[2] / tot_var else NA_real_
constrained_var <- sum(eig_con)
cap1_pct_constr <- if (constrained_var > 0) 100 * eig_con[1] / constrained_var else NA_real_
cap2_pct_constr <- if (constrained_var > 0 && length(eig_con) >= 2) 100 * eig_con[2] / constrained_var else NA_real_

cat(sprintf("\nAxis contributions:\n  CAP1 = %.2f%% of TOTAL (%.2f%% of constrained)\n  CAP2 = %.2f%% of TOTAL (%.2f%% of constrained)\n  Constrained total (unadj. R2) = %.2f%%\n  Adjusted R2 = %.2f%%\n",
            cap1_pct_total, cap1_pct_constr,
            cap2_pct_total, cap2_pct_constr,
            100 * constrained_var / tot_var,
            100 * R2s$adj.r.squared))

# SEPARATE TESTS: marginal and unique effects 
perm1k <- permute::how(nperm = 999)
cat("\n===============================\n Partial dbRDA — SEPARATE TESTS\n")

## (A) Marginal (single-term) effects: term | (log_reads + natman)
mod_decay_marg <- capscale(
  otu_fall ~ decay_stage + Condition(log_reads + natman),
  data = meta_fall, distance = "robust.aitchison"
)
mod_size_marg <- capscale(
  otu_fall ~ size + Condition(log_reads + natman),
  data = meta_fall, distance = "robust.aitchison"
)

cat("\nMarginal effects (single-term | log_reads + natman):\n")
print(anova.cca(mod_decay_marg, permutations = perm1k))
print(anova.cca(mod_size_marg,  permutations = perm1k))
cat(sprintf("  decay  — R2 = %.3f (adj %.3f)\n",
            RsquareAdj(mod_decay_marg)$r.squared,
            RsquareAdj(mod_decay_marg)$adj.r.squared))
cat(sprintf("  size   — R2 = %.3f (adj %.3f)\n",
            RsquareAdj(mod_size_marg)$r.squared,
            RsquareAdj(mod_size_marg)$adj.r.squared))

## (B) Unique (partial) effects: term | (other microhabitat + log_reads + natman)
mod_decay_unique <- capscale(
  otu_fall ~ decay_stage + Condition(size + log_reads + natman),
  data = meta_fall, distance = "robust.aitchison"
)
mod_size_unique <- capscale(
  otu_fall ~ size + Condition(decay_stage + log_reads + natman),
  data = meta_fall, distance = "robust.aitchison"
)

cat("\nUnique effects (term | other microhabitat + log_reads + natman):\n")
print(anova.cca(mod_decay_unique, permutations = perm1k))
print(anova.cca(mod_size_unique,  permutations = perm1k))
cat(sprintf("  decay  — R2 = %.3f (adj %.3f)\n",
            RsquareAdj(mod_decay_unique)$r.squared,
            RsquareAdj(mod_decay_unique)$adj.r.squared))
cat(sprintf("  size   — R2 = %.3f (adj %.3f)\n",
            RsquareAdj(mod_size_unique)$r.squared,
            RsquareAdj(mod_size_unique)$adj.r.squared))

# ------------------------------------------------------------------------

sc_sites <- scores(mod_decay, display = "sites", choices = 1:2)
plot_df <- as_tibble(sc_sites, rownames = "sample") %>%
  left_join(meta_fall %>% select(sample, decay_stage, size), by = "sample")

decay_cols <- c(EARLY = "#B58900", LATE = "#6A51A3")
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


# =============================================================================
# 2) FALLEN FWD ↔ SOIL pairing: diagnostics, distances, and signature tests
# =============================================================================
# Safety: exclude FRUITBODY everywhere for pairing logic
META1_nofb <- META1 %>% filter(is.na(dw_type) | dw_type != "FRUITBODY") %>% droplevels()

# FWD that *claim* a soil core
fall_yes <- meta_fall %>%
  filter(soil_core == "yes") %>%
  select(fwd_sample = sample, decay_stage, size, natman, reads_filt) %>%
  distinct()

# SOIL rows where soil_core stores the FWD sample name
soil_map <- META1_nofb %>%
  filter(dw_type == "SOIL", !is.na(soil_core), !soil_core %in% c("yes", "no")) %>%
  transmute(soil_sample = sample, fwd_sample = soil_core) %>%
  distinct()

pair_key <- fall_yes %>% left_join(soil_map, by = "fwd_sample")

dup_fwd      <- soil_map %>% count(fwd_sample) %>% filter(n > 1)
missing_soil <- pair_key %>% filter(is.na(soil_sample)) %>% pull(fwd_sample) %>% unique()
orphan_soil  <- soil_map %>% anti_join(fall_yes, by = "fwd_sample")

cat("\n--- FWD↔SOIL pairing diagnostics ---\n")
cat("FWD FALLEN marked 'yes' (n): ", nrow(fall_yes),  "\n", sep = "")
cat("SOIL rows referencing FWD (n):", nrow(soil_map), "\n")
cat("Complete matches (n):          ", sum(!is.na(pair_key$soil_sample)), "\n", sep = "")
if (nrow(dup_fwd))      { cat("\nWARNING duplicates:\n"); print(dup_fwd) }
if (length(missing_soil)){ cat("\nMissing SOIL for these FWDs:\n"); print(missing_soil) }
if (nrow(orphan_soil))  { cat("\nSOIL orphans (map to non-FALLEN yes):\n"); print(orphan_soil) }

pair_key_clean <- pair_key %>%
  filter(!is.na(soil_sample)) %>%
  anti_join(dup_fwd %>% select(fwd_sample), by = "fwd_sample") %>%
  distinct(fwd_sample, soil_sample, .keep_all = TRUE) %>%
  mutate(pair_id = fwd_sample)

cat("\nUsable 1:1 pairs: ", nrow(pair_key_clean), "\n", sep = "")
# =============================================================================
# 2a) Sanity check: do soils cluster by site (natman)?
# =============================================================================
restrict_to_pair_soils <- TRUE  # set FALSE to use all soils available

soil_meta <- META1_nofb %>%
  filter(dw_type == "SOIL") %>%
  { if (restrict_to_pair_soils) filter(., sample %in% pair_key_clean$soil_sample) else . } %>%
  droplevels()

stopifnot(nrow(soil_meta) >= 4)
X_soil <- otu_matrix_filt[soil_meta$sample, , drop = FALSE]
X_soil <- X_soil[, colSums(X_soil) > 0, drop = FALSE]

soil_meta$reads     <- rowSums(X_soil)
soil_meta$log_reads <- log1p(soil_meta$reads)
soil_meta$natman    <- droplevels(factor(soil_meta$natman))

# Distances
X_soil_pa <- (X_soil > 0) * 1
D_soil_j  <- vegdist(X_soil_pa, method = "jaccard", binary = TRUE)
D_soil_ra <- vegdist(X_soil,    method = "robust.aitchison")
X_soil_t  <- tss(X_soil)
D_soil_b  <- vegdist(X_soil_t,  method = "bray")

cat("\n--- SOIL-by-site sanity checks ---\n")
cat("Soil samples (used here):", nrow(soil_meta),
    "| Sites (natman):", nlevels(soil_meta$natman), "\n")

# (1) Partial dbRDA (robust Aitchison): natman | log_reads
perm1k <- permute::how(nperm = 999)
mod_soil <- capscale(X_soil ~ natman + Condition(log_reads),
                     data = soil_meta, distance = "robust.aitchison")
an_soil <- anova.cca(mod_soil, permutations = perm1k, by = "terms")
print(an_soil)

R2s_soil <- RsquareAdj(mod_soil)
cat(sprintf("dbRDA (robust Aitchison): unadj R2=%.3f, adj R2=%.3f\n",
            R2s_soil$r.squared, R2s_soil$adj.r.squared))

# (2) ANOSIM on PA-Jaccard (complements dbRDA; rank-based)
anosi <- anosim(D_soil_j, grouping = soil_meta$natman, permutations = 999)
cat(sprintf("ANOSIM (Jaccard PA): R=%.3f, p=%.4f\n", anosi$statistic, anosi$signif))

# (3) Within vs between site mean similarity (PA-Jaccard), label-permutation test
nat <- soil_meta$natman
D    <- as.matrix(D_soil_j)
S    <- 1 - D  # similarity

idx <- which(row(D) < col(D))
same_site <- nat[row(D)[idx]] == nat[col(D)[idx]]
within_vals  <- S[idx][same_site]
between_vals <- S[idx][!same_site]

obs_diff <- mean(within_vals) - mean(between_vals)

set.seed(42)
nperm <- 4999L
perm_diff <- numeric(nperm)
for (p in seq_len(nperm)) {
  gperm <- sample(nat)  # permute site labels
  same_p <- gperm[row(D)[idx]] == gperm[col(D)[idx]]
  perm_diff[p] <- mean(S[idx][same_p]) - mean(S[idx][!same_p])
}
p_within <- (sum(perm_diff >= obs_diff) + 1) / (nperm + 1)  # RIGHT tail: within > between

cat(sprintf("Within–between Jaccard (PA) mean diff = %.3f; p_perm = %.4f\n",
            obs_diff, p_within))

# (4) PERMDISP (beta dispersion) — is dispersion comparable across sites?
bd <- betadisper(D_soil_j, group = soil_meta$natman, type = "centroid")
bd_aov  <- anova(bd)
bd_perm <- permutest(bd, permutations = 999)
cat("PERMDISP (Jaccard PA) — ANOVA p:", signif(bd_aov$`Pr(>F)`[1], 3),
    "| perm p:", signif(bd_perm$tab[1, "Pr(>F)"], 3), "\n")

# Optional quick summary table
soil_within_between <- tibble(
  metric      = "Jaccard_PA",
  mean_within = mean(within_vals), mean_between = mean(between_vals),
  diff        = obs_diff, p_perm = p_within
)
print(soil_within_between)

# There is a detectable site (natman) signal in soils, especially in presence/absence space,
# but with unequal dispersion across sites and small sample size.
# Consequently,proceeding with site-blocked permutations for the FWD–SOIL pairing tests is appropriate and conservative;
# any significant pairwise “soil signature” then reflects similarity beyond what is expected among soils from the same site.

# ---- Per-pair distances (Bray on TSS; Jaccard on PA) -------------------
X_fwd  <- otu_matrix_filt[pair_key_clean$fwd_sample,  , drop = FALSE]
X_soil <- otu_matrix_filt[pair_key_clean$soil_sample, , drop = FALSE]
Xf_tss <- tss(X_fwd) ; Xs_tss <- tss(X_soil)

pair_bray <- map_dbl(seq_len(nrow(Xf_tss)), function(i)
  as.numeric(vegdist(rbind(Xf_tss[i,], Xs_tss[i,]), method = "bray"))[1])
pair_jac  <- map_dbl(seq_len(nrow(Xf_tss)), function(i) {
  A <- (X_fwd[i,] > 0) * 1; B <- (X_soil[i,] > 0) * 1
  as.numeric(vegdist(rbind(A,B), method = "jaccard", binary = TRUE))[1]
})
pair_dists <- pair_key_clean %>%
  mutate(reads_fwd  = rowSums(X_fwd),
         reads_soil = rowSums(X_soil),
         bray_tss   = pair_bray,          # distance
         jaccard_pa = pair_jac) %>%       # distance
  arrange(desc(bray_tss))
cat("\nHead of per-pair distances:\n"); print(head(pair_dists, 10), n = 10)

# ---- Paired dbRDA (visualization) --------------------------------------
O <- rbind(Xf_tss, Xs_tss)
type    <- factor(rep(c("FWD","SOIL"), each = nrow(Xf_tss)))
pair_id <- factor(rep(pair_key_clean$pair_id, times = 2))
dist_bray <- vegdist(O, method = "bray")

perm <- permute::how(nperm = 999); setBlocks(perm) <- pair_id
adon <- adonis2(dist_bray ~ type, permutations = perm)
cat("\nadonis2 (Bray, TSS; perms blocked by pair):\n"); print(adon)

mod_cap <- capscale(O ~ type + Condition(pair_id), distance = "bray")
an_terms <- anova.cca(mod_cap, permutations = perm, by = "terms")
cat("\nPaired dbRDA terms test:\n"); print(an_terms)

sc <- scores(mod_cap, display = "sites", choices = 1:2)
plt_df <- as_tibble(sc, rownames = "rowname") %>%
  mutate(type = type, pair_id = pair_id) %>%
  select(pair_id, type, CAP1 = CAP1, MDS1 = MDS1)
seg_df <- plt_df %>% pivot_wider(names_from = type, values_from = c(CAP1, MDS1))

gg_pairs <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey75") +
  geom_segment(data = seg_df,
               aes(x = CAP1_FWD, y = MDS1_FWD, xend = CAP1_SOIL, yend = MDS1_SOIL),
               inherit.aes = FALSE, linewidth = 0.3, alpha = 0.6) +
  geom_point(data = plt_df, aes(CAP1, MDS1, shape = type), size = 2.6, alpha = 0.95) +
  coord_equal() +
  labs(title = "FWD (fallen) vs matched SOIL — paired dbRDA",
       subtitle = "Model: type + Condition(pair_id); distance = Bray (TSS-normalized)",
       x = "CAP1", y = "MDS1", shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
print(gg_pairs)
ggsave(file.path(plots_dir, "fwd_soil_pairs_capscale.png"),
       gg_pairs, width = 6.6, height = 5.2, dpi = 300)

# Even when each FWD is compared to its own soil core from the same micro-site,
# communities are systematically different between substrates (wood vs soil).
# The effect is not driven by a few extreme pairs:
# blocking by pair and recovering the same signal with both adonis2 and dbRDA implies a consistent,
# directional compositional shift from FWD to SOIL.

# =============================================================================
# 3) Soil “signature” tests (global, per-pair, nearest-soil) — Jaccard PA
# Even though FWD and SOIL differ strongly overall, we asked whether each fallen-wood sample retains a detectable imprint of its own local soil—i.e.,
# whether FWDᵢ shares more taxa with its matched soil than expected given the soil community at that site.
# =============================================================================
# Cross-Jaccard matrix (FWD rows × SOIL cols; PA)
A_fwd  <- (otu_matrix_filt[pair_key_clean$fwd_sample,  , drop = FALSE] > 0) * 1
B_soil <- (otu_matrix_filt[pair_key_clean$soil_sample, , drop = FALSE] > 0) * 1
keep_taxa <- colSums(rbind(A_fwd, B_soil)) > 0
A_fwd  <- A_fwd[,  keep_taxa, drop = FALSE]
B_soil <- B_soil[, keep_taxa, drop = FALSE]

intersect_mat <- A_fwd %*% t(B_soil)                 # |A∩B|
nA <- rowSums(A_fwd); nB <- rowSums(B_soil)
union_mat <- outer(nA, nB, "+") - intersect_mat      # |A∪B|
jaccard   <- intersect_mat / pmax(union_mat, 1)

obs_pair_j <- diag(jaccard)                          # similarity in [0,1]

# Global pairing test (blocked by site)
nperm <- 4999L
grp_site <- as.factor(pair_key_clean$natman)

perm_stats    <- numeric(nperm)
per_fwd_null  <- matrix(NA_real_, nrow = nrow(A_fwd), ncol = nperm)
for (p in seq_len(nperm)) {
  perm_idx <- permute_within(grp_site)
  perm_j   <- jaccard[cbind(seq_len(nrow(jaccard)), perm_idx)]
  perm_stats[p]     <- mean(perm_j, na.rm = TRUE)
  per_fwd_null[, p] <- perm_j
}
obs_global <- mean(obs_pair_j, na.rm = TRUE)
obs_global
p_global   <- (sum(perm_stats >= obs_global) + 1) / (nperm + 1)  # RIGHT-tail
p_global
# Per-pair Monte Carlo p (RIGHT-tail: matched soil is more similar than null)
# ~compare each FWD Jaccard-observed to its site-blocked null column
p_mc <- (rowSums(per_fwd_null >= obs_pair_j) + 1) / (nperm + 1)
p_mc_BH <- p.adjust(p_mc, method = "BH")
p_mc_BH
# Hypergeometric enrichment for overlap counts
# test whether the overlap count is enriched beyond random co-draws from the pool
M      <- ncol(A_fwd)
k_vec  <- diag(intersect_mat)
p_hyp  <- 1 - phyper(q = pmax(k_vec - 1, 0), m = nB, n = M - nB, k = nA)
p_hyp_adj <- p.adjust(p_hyp, method = "BH")
p_hyp_adj
# Nearest-soil accuracy (rank 1?)
rank_matrix <- apply(jaccard, 1, function(x) rank(-x, ties.method = "min")) %>% t()
is_top1  <- diag(rank_matrix) == 1
top1_obs <- sum(is_top1)
top1_perm <- integer(nperm)
for (p in seq_len(nperm)) {
  perm_idx    <- permute_within(grp_site)
  ranks_i     <- rank_matrix[cbind(seq_len(nrow(rank_matrix)), perm_idx)]
  top1_perm[p] <- sum(ranks_i == 1)
}
p_top1 <- (sum(top1_perm >= top1_obs) + 1) / (nperm + 1)

pair_results <- tibble(
  fwd_sample        = pair_key_clean$fwd_sample,
  soil_sample       = pair_key_clean$soil_sample,
  natman            = pair_key_clean$natman,
  jaccard_obs       = obs_pair_j,
  jaccard_null_mean = rowMeans(per_fwd_null),
  jaccard_null_sd   = apply(per_fwd_null, 1, sd),
  p_perm_pair       = p_mc,
  p_perm_pair_BH    = p_mc_BH,
  k_overlap         = k_vec,
  n_fwd             = nA,
  n_soil            = nB,
  M_pool            = M,
  p_hyper           = p_hyp,
  p_hyper_BH        = p_hyp_adj
) %>%
  mutate(
    z_effect = (jaccard_obs - jaccard_null_mean) / pmax(jaccard_null_sd, 1e-6)
  ) %>%
  arrange(p_perm_pair_BH, desc(z_effect))

cat("\n=== Soil signature: global pairing test (Jaccard, PA) ===\n")
cat(sprintf("Observed mean paired Jaccard = %.3f\n", obs_global))
cat(sprintf("Null model mean (±SD)  = %.3f ± %.3f\n", mean(perm_stats), sd(perm_stats)))
cat(sprintf("Permutation p-value (site-blocked) = %.4f\n", p_global))
cat(sprintf("\nTop-1 nearest-soil accuracy = %d of %d pairs (p_perm = %.4f)\n",
            top1_obs, nrow(pair_results), p_top1))
cat("\nPer-pair examples (top 10 by smallest BH-adjusted p_perm_pair):\n")
print(head(pair_results %>% arrange(p_perm_pair_BH), 11), n = 11)

readr::write_tsv(pair_results, file.path(tables_dir, "fwd_soil_signature.tsv"))

# Null vs observed plot (global stat)
df_plot <- tibble(stat = perm_stats, what = "null") |>
  add_row(stat = obs_global, what = "observed")
g_sig <- ggplot(df_plot, aes(stat, fill = what)) +
  geom_histogram(data = df_plot |> filter(what == "null"), bins = 40, alpha = 0.6) +
  geom_vline(xintercept = obs_global, linewidth = 0.5) +
  labs(title = "FALLEN FWD vs SOIL — global pairing test",
       subtitle = "Statistic = mean Jaccard(FWDᵢ, SOILᵢ); permutations blocked by site (natman)",
       x = "Mean paired Jaccard", y = "Count", fill = NULL) +
  theme_minimal(base_size = 12)
print(g_sig)
ggsave(file.path(plots_dir, "fwd_soil_signature_null_vs_obs.png"),
       g_sig, width = 6.3, height = 4.6, dpi = 300)

# Despite low absolute similarities (typical for cross-substrate comparisons), 
# there is a consistent, above-chance soil imprint at the presence/absence level: 
# many FWDs share more taxa with their own local soil than with other same-site soils.
# The agreement between the global mean, per-pair Monte-Carlo p-values, overlap enrichment,
# and nearest-soil accuracy supports a subtle but reproducible local soil signature—superimposed on a large, 
# systematic substrate difference documented by the paired Bray analyses.

# =============================================================================
# 4) Obs vs null plots by pair: Jaccard (PA), Bray distance, robust Aitchison
#     (metric-specific one-sided tests + asterisk labeling)
# =============================================================================

alpha_sig <- 0.05

# Base ordering & labels scaffold (used for all three plots)
meta_map <- meta_fall %>% select(fwd_sample = sample, decay_stage, size)
order_df <- pair_key_clean %>%
  select(natman, fwd_sample) %>%
  left_join(meta_map, by = "fwd_sample") %>%
  mutate(
    decay_stage = forcats::fct_relevel(forcats::fct_na_value_to_level(as.factor(decay_stage), "NA"),
                                       "EARLY","LATE"),
    size        = forcats::fct_relevel(forcats::fct_na_value_to_level(as.factor(size), "NA"),
                                       "FINE","VERY_FINE")
  ) %>%
  group_by(natman) %>%
  arrange(decay_stage, size, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(axis_id = paste(natman, fwd_sample, sep = "::"),
         base_lab = paste(as.character(decay_stage), "·", as.character(size)))
axis_levels <- order_df$axis_id

# ----------------------------- JACCARD (similarity, RIGHT-tail) ---------
df_jac <- tibble(
  fwd_sample = pair_key_clean$fwd_sample,
  natman     = pair_key_clean$natman,
  obs        = obs_pair_j,
  null_mean  = rowMeans(per_fwd_null),
  null_sd    = apply(per_fwd_null, 1, sd)
) %>%
  left_join(order_df, by = c("fwd_sample","natman")) %>%
  mutate(
    sig  = pair_results$p_perm_pair_BH[match(fwd_sample, pair_results$fwd_sample)] < alpha_sig,
    axis_lab = paste0(base_lab, if_else(sig, " *", ""))
  )
# ----- Jaccard percentile CIs (RIGHT-tail) -------------------------------------
lo_q_j <- apply(per_fwd_null, 1, quantile, probs = 0.025, na.rm = TRUE)
hi_q_j <- apply(per_fwd_null, 1, quantile, probs = 0.975, na.rm = TRUE)
df_jac <- df_jac %>%
  mutate(lo95 = lo_q_j, hi95 = hi_q_j)  # already in [0,1]
df_jac
gg_jac <- ggplot(df_jac, aes(x = factor(axis_id, levels = axis_levels))) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6, alpha = 0.4) +
  geom_point(aes(y = null_mean, color = "null mean"), size = 2, alpha = 0.5) +
  geom_point(aes(y = obs,       color = "observed"),  size = 2.6, alpha = 0.75) +
  scale_color_manual(name = NULL, values = c("null mean" = "grey30", "observed" = "black")) +
  facet_grid(natman ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = setNames(df_jac$axis_lab, df_jac$axis_id)) +
  coord_flip() +
  labs(title = "FWD–SOIL pairing: observed Jaccard vs site-blocked null",
       x = "Decay stage · Size", y = "Jaccard(FWD, matched SOIL)") +
  theme_minimal(base_size = 12) +
  theme(strip.text.y = element_text(face = "bold"), panel.spacing.y = unit(6, "pt"))
print(gg_jac)
ggsave(file.path(plots_dir, "fwd_soil_pairs_obs_vs_null_Jaccard.png"),
       gg_jac, width = 7.0, height = 6.2, dpi = 300)

# ----------------------------- BRAY (distance, LEFT-tail) ---------------
# Build similarity null then convert to distances for LEFT-tail testing
Af <- tss(otu_matrix_filt[pair_key_clean$fwd_sample,  , drop = FALSE])  # FWD
Bs <- tss(otu_matrix_filt[pair_key_clean$soil_sample, , drop = FALSE])  # SOIL
bray_cross <- as.matrix(vegdist(rbind(Af, Bs), method = "bray"))
nP <- nrow(Af)
bray_cross <- bray_cross[seq_len(nP), nP + seq_len(nP)]
bray_sim   <- 1 - bray_cross
obs_pair_bray <- diag(bray_sim)

nperm <- 4999L
perm_stats_bray   <- numeric(nperm)
per_fwd_null_bray <- matrix(NA_real_, nrow = nP, ncol = nperm)
for (p in seq_len(nperm)) {
  perm_idx <- permute_within(grp_site)  # block by site
  perm_s   <- bray_sim[cbind(seq_len(nrow(bray_sim)), perm_idx)]
  perm_stats_bray[p]     <- mean(perm_s, na.rm = TRUE)
  per_fwd_null_bray[, p] <- perm_s
}

# Convert to distance and compute LEFT-tail p's (observed < null)
obs_bray_dist       <- 1 - obs_pair_bray
per_fwd_null_bray_d <- 1 - per_fwd_null_bray

valid_perms_bray <- rowSums(is.finite(per_fwd_null_bray_d))
comp_bray        <- sweep(per_fwd_null_bray_d, 1, obs_bray_dist, FUN = "<=")
p_mc_brayD       <- (rowSums(comp_bray, na.rm = TRUE) + 1) / (valid_perms_bray + 1)
q_mc_brayD       <- p.adjust(p_mc_brayD, method = "BH")

perm_stats_brayD <- colMeans(per_fwd_null_bray_d, na.rm = TRUE)
obs_global_brayD <- mean(obs_bray_dist, na.rm = TRUE)
n_perm_ok_bray   <- sum(is.finite(perm_stats_brayD))
p_global_brayD   <- (sum(perm_stats_brayD <= obs_global_brayD, na.rm = TRUE) + 1) / (n_perm_ok_bray + 1)
cat(sprintf("\nBray distance — global mean = %.3f; p_perm (left tail) = %.4f\n",
            obs_global_brayD, p_global_brayD))

df_bray <- tibble(
  fwd_sample = pair_key_clean$fwd_sample,
  natman     = pair_key_clean$natman,
  obs        = obs_bray_dist,
  null_mean  = rowMeans(per_fwd_null_bray_d, na.rm = TRUE),
  null_sd    = apply(per_fwd_null_bray_d, 1, sd)
) %>%
  left_join(order_df, by = c("fwd_sample","natman")) %>%
  mutate(
    sig  = q_mc_brayD < alpha_sig,
    axis_lab = paste0(base_lab, if_else(sig, " *", ""))
  )
lo_q_b <- apply(per_fwd_null_bray_d, 1, quantile, probs = 0.025, na.rm = TRUE)
hi_q_b <- apply(per_fwd_null_bray_d, 1, quantile, probs = 0.975, na.rm = TRUE)
df_bray <- df_bray %>%
  mutate(lo95 = pmax(lo_q_b, 0), hi95 = pmin(hi_q_b, 1))

df_bray
gg_pair_bray <- ggplot(df_bray, aes(x = factor(axis_id, levels = axis_levels))) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6, alpha = 0.4) +
  geom_point(aes(y = null_mean, color = "null mean"), size = 2, alpha = 0.7) +
  geom_point(aes(y = obs,       color = "observed"),  size = 2.6) +
  scale_color_manual(name = NULL, values = c("null mean" = "grey30", "observed" = "black")) +
  facet_grid(natman ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = setNames(df_bray$axis_lab, df_bray$axis_id)) +
  coord_flip() +
  labs(title = "FWD–SOIL pairing: observed Bray–Curtis distance vs site-blocked null",
       x = "Decay stage · Size", y = "Bray–Curtis distance (TSS-normalized)") +
  theme_minimal(base_size = 12) +
  theme(strip.text.y = element_text(face = "bold"),
        panel.spacing.y = unit(6, "pt"))
print(gg_pair_bray)
ggsave(file.path(plots_dir, "fwd_soil_pairs_obs_vs_null_brayDIST.png"),
       gg_pair_bray, width = 7.0, height = 6.2, dpi = 300)

# ----------------------------- ROBUST AITCHISON (distance, LEFT-tail) ---
ait_cross <- as.matrix(vegdist(rbind(Xa, Xb), method = "robust.aitchison"))
Xa <- otu_matrix_filt[pair_key_clean$fwd_sample,  , drop = FALSE]
Xb <- otu_matrix_filt[pair_key_clean$soil_sample, , drop = FALSE]
ait_cross <- as.matrix(vegdist(rbind(Xa, Xb), method = "robust.aitchison"))

# Cross-distance matrix: FWD rows x SOIL cols
ait_full   <- as.matrix(vegdist(rbind(Xa, Xb), method = "robust.aitchison"))
nP         <- nrow(Xa)
ait_cross  <- ait_full[seq_len(nP), nP + seq_len(nP), drop = FALSE]
obs_ait    <- diag(ait_cross)

# BLOCKING: make sure the column (soil) site vector matches column order
soil_site_cols <- META1_nofb %>%
  filter(dw_type == "SOIL", sample %in% pair_key_clean$soil_sample) %>%
  slice(match(pair_key_clean$soil_sample, sample)) %>%           # reorder to columns
  pull(natman) %>% factor() %>% droplevels()

# Safety check (should be TRUE; if FALSE, use soil_site_cols for blocking)
cat("Aitchison blocking check (rows vs soil cols same site order): ",
    isTRUE(all.equal(as.character(pair_key_clean$natman), as.character(soil_site_cols))), "\n")

# Permute SOIL columns within *soil_site_cols*
set.seed(42)
nperm <- 4999L
per_fwd_null_ait <- matrix(NA_real_, nrow = nP, ncol = nperm)
for (p in seq_len(nperm)) {
  perm_idx <- permute_within(soil_site_cols)  # <-- column-wise grouping
  per_fwd_null_ait[, p] <- ait_cross[cbind(seq_len(nP), perm_idx)]
}

# Per-pair LEFT-tail p (obs smaller than null), + BH
valid_perms_ait <- rowSums(is.finite(per_fwd_null_ait))
comp_ait        <- sweep(per_fwd_null_ait, 1, obs_ait, FUN = "<=")
p_mc_ait        <- (rowSums(comp_ait, na.rm = TRUE) + 1) / (valid_perms_ait + 1)
q_mc_ait        <- p.adjust(p_mc_ait, method = "BH")

# Global LEFT-tail p (mean)
perm_stats_ait <- colMeans(per_fwd_null_ait, na.rm = TRUE)
obs_global_ait <- mean(obs_ait, na.rm = TRUE)
n_perm_ok_ait  <- sum(is.finite(perm_stats_ait))
p_global_ait   <- (sum(perm_stats_ait <= obs_global_ait, na.rm = TRUE) + 1) / (n_perm_ok_ait + 1)

# Percentile CIs (more faithful than mean±1.96*sd)
lo_q_ait <- apply(per_fwd_null_ait, 1, quantile, probs = 0.025, na.rm = TRUE)
hi_q_ait <- apply(per_fwd_null_ait, 1, quantile, probs = 0.975, na.rm = TRUE)

# Build DF + diagnostics
df_ait <- tibble(
  fwd_sample = pair_key_clean$fwd_sample,
  natman     = pair_key_clean$natman,
  obs        = obs_ait,
  null_mean  = rowMeans(per_fwd_null_ait, na.rm = TRUE),
  null_sd    = apply(per_fwd_null_ait, 1, sd),
  lo95       = lo_q_ait,
  hi95       = hi_q_ait,
  perm_n     = valid_perms_ait,
  p_perm     = p_mc_ait,
  q_perm     = q_mc_ait,
  delta      = obs_ait - rowMeans(per_fwd_null_ait, na.rm = TRUE)   # should be <0 if sig
) %>%
  left_join(order_df, by = c("fwd_sample","natman")) %>%
  mutate(sig = q_perm < alpha_sig,
         axis_lab = paste0(base_lab, if_else(sig, " *", "")))

cat(sprintf("Robust Aitchison — global mean = %.3f; null mean = %.3f ± %.3f; p_perm (left-tail) = %.4f\n",
            obs_global_ait, mean(perm_stats_ait, na.rm = TRUE), sd(perm_stats_ait, na.rm = TRUE), p_global_ait))
# Quick sanity print: direction vs p
print(df_ait %>% select(fwd_sample, natman, obs, null_mean, delta, p_perm, q_perm, perm_n) %>% arrange(q_perm))

# Plot (uses percentile CIs)
gg_pair_ait <- ggplot(df_ait, aes(x = factor(axis_id, levels = axis_levels))) +
  geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 0.6, alpha = 0.4) +
  geom_point(aes(y = null_mean, color = "null mean"), size = 2, alpha = 0.7) +
  geom_point(aes(y = obs,       color = "observed"),  size = 2.6) +
  scale_color_manual(name = NULL, values = c("null mean" = "grey30", "observed" = "black")) +
  facet_grid(natman ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = setNames(df_ait$axis_lab, df_ait$axis_id)) +
  coord_flip() +
  labs(title = "FWD–SOIL pairing: Robust Aitchison distance vs site-blocked null",
       x = "Decay stage · Size", y = "Robust Aitchison distance") +
  theme_minimal(base_size = 12) +
  theme(strip.text.y = element_text(face = "bold"),
        panel.spacing.y = unit(6, "pt"))
print(gg_pair_ait)
ggsave(file.path(plots_dir, "fwd_soil_pairs_obs_vs_null_aitDIST.png"),
       gg_pair_ait, width = 7.0, height = 6.2, dpi = 300)

# =============================================================================
# 5) Taxon-level concordance: is a taxon in FWD more often when in its matched SOIL?
# =============================================================================
# Concordance test asks: “Is it present together in FWD & its matched SOIL more than site-blocked chance?”
A <- (otu_matrix_filt[pair_key_clean$fwd_sample,  , drop = FALSE] > 0) * 1
B <- (otu_matrix_filt[pair_key_clean$soil_sample, , drop = FALSE] > 0) * 1
keep <- colSums(A) + colSums(B) >= 1
A <- A[, keep, drop = FALSE]; B <- B[, keep, drop = FALSE]

# Observed matched concurrence per taxon
obs_conc <- colSums((A == 1) & (B == 1))   # FWD=1 & SOIL=1
obs_disc <- colSums((A == 1) & (B == 0))   # FWD only

# Site-blocked null by permuting SOIL rows
nperm <- 1999L
grp_tax <- as.factor(pair_key_clean$natman)
null_conc <- matrix(0L, nrow = ncol(A), ncol = nperm)
for (p in seq_len(nperm)) {
  perm_idx <- permute_within(grp_tax)
  Bp <- B[perm_idx, , drop = FALSE]
  null_conc[, p] <- colSums((A == 1) & (Bp == 1))
}

# Per-taxon p (RIGHT-tail: observed ≥ null) + BH
hits   <- sweep(null_conc, 1, obs_conc, FUN = ">=")  # taxa x nperm logical
p_taxa <- (rowSums(hits) + 1) / (nperm + 1)
q_taxa <- p.adjust(p_taxa, method = "BH")

# Effect sizes vs null
null_mean  <- rowMeans(null_conc)
null_sd    <- sqrt(pmax(rowMeans(null_conc^2) - null_mean^2, 0))
z_taxa     <- (obs_conc - null_mean) / pmax(null_sd, 1e-9)
delta_taxa <- obs_conc - null_mean

# Attach any available taxonomy
tax_keep <- intersect(c("sh_code","lowest_taxon","sh_link"), names(tax))
tax_slim <- if (length(tax_keep)) distinct(tax[, tax_keep, drop = FALSE]) else tibble(sh_code = colnames(A))

tax_sig <- tibble(
  sh_code     = colnames(A),
  conc_pairs  = obs_conc,
  fwd_only    = obs_disc,
  null_mean   = null_mean,
  null_sd     = null_sd,
  delta_pairs = delta_taxa,
  z_pairs     = z_taxa,
  p_perm      = p_taxa,
  q_perm      = q_taxa
) %>%
  left_join(tax_slim, by = "sh_code") %>%
  arrange(q_perm, desc(z_pairs))

print(head(tax_sig, 20), n = 20)
readr::write_tsv(tax_sig, file.path(tables_dir, "fwd_soil_taxon_concordance.tsv"))

# ----------------------- Simple paired co-occurrence (McNemar) ----------
# McNemar asks: “Across pairs, is it preferentially FWD-only or SOIL-only?”
both      <- colSums((A == 1) & (B == 1))
fwd_only  <- colSums((A == 1) & (B == 0))
soil_only <- colSums((A == 0) & (B == 1))
neither   <- colSums((A == 0) & (B == 0))

mcnemar_safe <- function(a, b, c, d) {
  if ((b + c) == 0) return(1)  # no discordant pairs → no asymmetry evidence
  mcnemar.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), correct = TRUE)$p.value
}

p_mcn <- mapply(mcnemar_safe, both, fwd_only, soil_only, neither)
q_mcn <- p.adjust(p_mcn, method = "BH")
logOR_mcn <- log((fwd_only + 0.5) / (soil_only + 0.5))  # >0 ⇒ more FWD-only

tax_mcnemar <- tibble(
  sh_code    = colnames(A),
  both, fwd_only, soil_only, neither,
  logOR_mcn, p_mcnemar = p_mcn, q_mcnemar = q_mcn
) %>% arrange(q_mcnemar, desc(logOR_mcn)) %>% 
  left_join(tax, by = "sh_code") %>% 
  select(sh_code  ,both , fwd_only, soil_only, neither, logOR_mcn, p_mcnemar,q_mcnemar ,lowest_taxon) %>% 
  filter(p_mcnemar<0.01, q_mcnemar<0.1)

cat("\nMcNemar co-occurrence (top by BH q):\n")
print(tax_mcnemar, n = 30)
# Significant positive logOR → the taxon tends to be detected in FWD without the matched SOIL (FWD-leaning occurrence).
# Significant negative logOR → detected in SOIL without matched FWD (SOIL-leaning)


# readr::write_tsv(tax_mcnemar, file.path(tables_dir, "taxa_mcnemar_cooccurrence.tsv"))
tax_sig2 <- tax_sig %>%
  select(sh_code, conc_pairs, fwd_only, delta_pairs, z_pairs, p_perm, q_perm, lowest_taxon) %>%
  left_join(tax_mcnemar, by = "sh_code")
tax_sig2
cat("\n=== FALLEN FWD workflow complete. Plots:", plots_dir, " — Tables:", tables_dir, " ===\n")
