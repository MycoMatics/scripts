#### SNAG ANALYSIS – aspect (NORTH vs SOUTH) ##################################
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(forcats)
  library(vegan);     library(permute)
  library(glmmTMB);   library(emmeans);  library(performance); library(broom.mixed)
  library(ggplot2);   library(scales)   # grey_pal()
})
set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")
stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"))
outdir <- "plots/SNAG"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# 0) Choose OTU source (filtered preferred) -----------------------------------
otu_source <- otu_matrix_filt

# 1) Build SNAG metadata ------------------------------------------------------
#    Aspect = NORTH/SOUTH at DBH (paired cores)
SNAG_meta <- META1 %>%
  filter(dw_type == "SNAG") %>%
  mutate(
    # Map aspect from position (N/S); fix known mismatches
    aspect = case_when(
      position %in% c("NORTH","N") ~ "NORTH",
      position %in% c("SOUTH","S") ~ "SOUTH",
      TRUE ~ NA_character_
    ),
    # Keep decay_stage tidy (contextual)
    decay_stage = case_when(
      sample %in% c("P1_25","P1_26","P2_68","P2_72") ~ "AVERAGE",
      sample %in% c("P2_04","P2_10")                 ~ "EARLY",
      sample %in% c("MI_34","MI_35")                 ~ "LATE",
      TRUE ~ decay_stage
    ),
    aspect       = droplevels(factor(aspect, levels = c("NORTH","SOUTH"))),
    decay_stage  = droplevels(factor(decay_stage)),
    ds_at_drill  = droplevels(factor(as.character(ds_at_drill),
                                     levels = c("0","1","2","3","4","5"))),
    natman       = factor(gsub('"', "", as.character(natman)), ordered = FALSE)
  )

# 2) Align OTU table to SNAG samples + add read depth -------------------------
otu_snag <- otu_source[SNAG_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_snag), SNAG_meta$sample))
otu_snag <- otu_snag[, colSums(otu_snag) > 0, drop = FALSE]

SNAG_meta <- SNAG_meta %>%
  mutate(reads = rowSums(otu_snag), log_reads = log1p(reads)) %>%
  filter(!is.na(aspect), reads > 0) %>%
  droplevels()

otu_snag <- otu_snag[SNAG_meta$sample, , drop = FALSE]  # keep in sync
cat("SNAG n:", nrow(SNAG_meta), " | OTUs:", ncol(otu_snag), "\n")

# 3) Colors -------------------------------------------------------------------

# 4) Distances (once) ---------------------------------------------------------
D_robAit      <- vegdist(otu_snag, method = "robust.aitchison")
comm_hell     <- decostand(otu_snag, method = "hellinger")
D_euclid_hell <- dist(comm_hell, method = "euclidean")
comm_rel      <- decostand(otu_snag, method = "total")
D_bray_rel    <- vegdist(comm_rel, method = "bray")
D_bray_relS   <- vegdist(sqrt(comm_rel), method = "bray")
D_jacc        <- vegdist((otu_snag > 0) * 1, method = "jaccard", binary = TRUE)

# 5) Alpha diversity — descriptive figures -----------------------------------
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
    labs(title = "Alpha diversity by SNAG aspect (DBH)", x = "Aspect", y = "Value", fill = "Aspect")
  ggsave(file.path(outdir, "alpha_by_aspect.png"), p_alpha_aspect, width=9, height=4.8, dpi=300)
}
p_alpha_aspect
if (exists("DS_colors") && !is.null(DS_colors) && nlevels(SNAG_meta$decay_stage) > 1) {
  p_alpha_ds <- alpha_snag %>%
    ggplot(aes(decay_stage, value, fill = decay_stage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = DS_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha diversity by SNAG decay_stage", x = "decay_stage", y = "Value", fill = "decay_stage")
  ggsave(file.path(outdir, "alpha_by_decay_stage.png"), p_alpha_ds, width=9, height=4.8, dpi=300)
}
p_alpha_ds
if (!is.null(DS_colors) && nlevels(SNAG_meta$ds_at_drill) > 1) {
  p_alpha_dsd <- alpha_snag %>%
    ggplot(aes(ds_at_drill, value, fill = ds_at_drill)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = ds_colors) +
    facet_wrap(~ metric, scales = "free_y") +
    theme_minimal(base_size = 12) +
    labs(title = "Alpha diversity by SNAG ds_at_drill", x = "ds_at_drill", y = "Value", fill = "ds_at_drill")
  ggsave(file.path(outdir, "alpha_by_ds_at_drill.png"), p_alpha_dsd, width=9, height=4.8, dpi=300)
}
p_alpha_dsd
# ============================
# 6) PERMANOVA (UNBLOCKED + paired)
# ============================
nperm <- 999
perm_unblocked <- how(nperm = nperm)

SNAG_meta_perma <- droplevels(SNAG_meta)
otu_snag_perma  <- otu_snag[SNAG_meta_perma$sample, , drop = FALSE]

D_robAit_p      <- vegdist(otu_snag_perma, method = "robust.aitchison")
comm_hell_p     <- decostand(otu_snag_perma, method = "hellinger")
D_euclid_hell_p <- dist(comm_hell_p, method = "euclidean")
comm_rel_p      <- decostand(otu_snag_perma, method = "total")
D_bray_rel_p    <- vegdist(comm_rel_p, method = "bray")
D_bray_relS_p   <- vegdist(sqrt(comm_rel_p), method = "bray")
D_jacc_p        <- vegdist((otu_snag_perma > 0) * 1, method = "jaccard", binary = TRUE)

perma_unblocked <- list(
  robAit  = adonis2(D_robAit_p      ~ log_reads + aspect + ds_at_drill + natman,
                    data = SNAG_meta_perma, by = "margin", permutations = perm_unblocked),
  hellEuc = adonis2(D_euclid_hell_p ~ log_reads + aspect + ds_at_drill + natman,
                    data = SNAG_meta_perma, by = "margin", permutations = perm_unblocked),
  relBray = adonis2(D_bray_rel_p    ~ log_reads + aspect + ds_at_drill + natman,
                    data = SNAG_meta_perma, by = "margin", permutations = perm_unblocked),
  relBrayS= adonis2(D_bray_relS_p   ~ log_reads + aspect + ds_at_drill + natman,
                    data = SNAG_meta_perma, by = "margin", permutations = perm_unblocked)
)
cat("\nUnblocked PERMANOVA — Aitchison/Euclid/Bray/Bray^0.5:\n")
print(perma_unblocked$robAit); 
print(perma_unblocked$hellEuc)
print(perma_unblocked$relBray);
print(perma_unblocked$relBrayS)

# Compact readouts for key terms
aspect_R2 <- purrr::map_dfr(perma_unblocked, ~{df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)}, .id = "distance") %>%
  filter(term == "aspect") %>% select(distance, R2, `Pr(>F)`)
ds_R2 <- purrr::map_dfr(perma_unblocked, ~{df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)}, .id = "distance") %>%
  filter(term == "ds_at_drill") %>% select(distance, R2, `Pr(>F)`)
natman_R2 <- purrr::map_dfr(perma_unblocked, ~{df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)}, .id = "distance") %>%
  filter(term == "natman") %>% select(distance, R2, `Pr(>F)`)
cat("\nUnblocked PERMANOVA — aspect:\n"); print(aspect_R2)
cat("\nUnblocked PERMANOVA — ds_at_drill:\n"); print(ds_R2)
cat("\nUnblocked PERMANOVA — natman R2 across distances:\n"); print(natman_R2)
readr::write_csv(natman_R2, file.path(outdir, "permanova_natman_R2_across_distances.csv"))

# Presence–absence robustness (Jaccard) — typically omit log_reads
cat("\nPresence–absence PERMANOVA (Jaccard; unblocked):\n")
perma_pa <- adonis2(D_jacc_p ~ aspect + ds_at_drill + natman,
                    data = SNAG_meta_perma, by = "margin", permutations = perm_unblocked)
print(perma_pa)

# Paired (within-snag) aspect test via strata ---------------------------------
cat("\nPaired PERMANOVA (Aitchison) — aspect within snag (expect NS):\n")
print(adonis2(D_robAit_p ~ log_reads + aspect, data = SNAG_meta_perma,
              by = "margin", permutations = how(nperm = nperm),
              strata = SNAG_meta_perma$natman))
cat("\nPaired PERMANOVA (Jaccard) — aspect within snag (expect NS):\n")
print(adonis2(D_jacc_p ~ aspect, data = SNAG_meta_perma,
              by = "margin", permutations = how(nperm = nperm),
              strata = SNAG_meta_perma$natman))

# 7) Constrained ordination (dbRDA) + Euler var-part --------------------------
# Microhab complex = aspect + ds_at_drill (condition on log_reads)
mod_all   <- capscale(otu_snag_perma ~ natman + aspect + ds_at_drill + Condition(log_reads),
                      data = SNAG_meta_perma, distance = "robust.aitchison")
mod_tree  <- capscale(otu_snag_perma ~ natman + Condition(log_reads),
                      data = SNAG_meta_perma, distance = "robust.aitchison")
mod_micro <- capscale(otu_snag_perma ~ aspect + ds_at_drill + Condition(log_reads),
                      data = SNAG_meta_perma, distance = "robust.aitchison")
mod_tree_c  <- capscale(otu_snag_perma ~ natman + Condition(aspect + ds_at_drill + log_reads),
                        data = SNAG_meta_perma, distance = "robust.aitchison")
mod_micro_c <- capscale(otu_snag_perma ~ aspect + ds_at_drill + Condition(natman + log_reads),
                        data = SNAG_meta_perma, distance = "robust.aitchison")

R2_all   <- RsquareAdj(mod_all)$r.squared
R2_tree  <- RsquareAdj(mod_tree)$r.squared
R2_micro <- RsquareAdj(mod_micro)$r.squared
R2_treeU <- RsquareAdj(mod_tree_c)$r.squared
R2_micU  <- RsquareAdj(mod_micro_c)$r.squared
cat("\nVariance fractions (unadjusted):\n")
print(c(all=R2_all, tree=R2_tree, micro=R2_micro, tree_unique=R2_treeU, micro_unique=R2_micU))

# Term-wise tests (blocked by natman where meaningful)
cat("\nPartial dbRDA tests (Aitchison):\n")
print(anova.cca(mod_all, permutations = how(nperm = nperm), by = "terms"))

# Euler partition (Snag vs Microhab complex)
make_euler_plot(
  R2_all   = R2_all,
  R2_tree  = R2_tree,
  R2_micro = R2_micro,
  R2_tree_u= R2_treeU,
  R2_mic_u = R2_micU,
  title_main = "Variance partitioning (SNAG, robust Aitchison dbRDA)",
  title_sub  = "Snag vs Microhab complex (aspect + ds | log_reads)",
  outfile    = file.path(outdir, "euler_variance_partitioning_SNAG.png")
)

# Bray versions (align with composition tests)
mod_ds_bray    <- capscale(otu_snag_perma ~ ds_at_drill + Condition(log_reads), data = SNAG_meta_perma, distance = "bray")
mod_ds_nm_bray <- capscale(otu_snag_perma ~ ds_at_drill + Condition(log_reads) + Condition(natman), data = SNAG_meta_perma, distance = "bray")
cat("\n--- dbRDA summaries (Bray) ---\n")
cat("Model: ds | log_reads (Bray)\n")
print(anova.cca(mod_ds_bray, by = "terms", permutations = how(nperm = nperm))); print(RsquareAdj(mod_ds_bray))
cat("\nModel: ds | log_reads + natman (Bray; within-Snag)\n")
print(anova.cca(mod_ds_nm_bray, by = "terms", permutations = how(nperm = nperm))); print(RsquareAdj(mod_ds_nm_bray))

summary(mod_ds_nm_bray)
alias(mod_ds_nm_bray)
class(mod_ds_nm_bray)
suppressWarnings(print(alias.cca(mod_ds_nm_bray)))
# This indicates that ds_at_drill is a factor with levels that are perfectly predicted by combinations of  natman  variable and other levels of ds_at_drill.

# 8) Dispersion (PERMDISP) checks --------------------------------------------
bd_ds <- betadisper(D_robAit_p, SNAG_meta_perma$ds_at_drill)
bd_as <- betadisper(D_robAit_p, SNAG_meta_perma$aspect)
bd_nm <- betadisper(D_robAit_p, SNAG_meta_perma$natman)
cat("\nPERMDISP (Aitchison) by ds_at_drill:\n"); print(permutest(bd_ds, permutations = how(nperm = nperm)))
cat("\nPERMDISP (Aitchison) by aspect:\n");     print(permutest(bd_as, permutations = how(nperm = nperm)))
cat("\nPERMDISP (Aitchison) by natman:\n");     print(permutest(bd_nm, permutations = how(nperm = nperm)))

# 9) Within-snag Δds diagnostics (power) --------------------------------------
ds_by_tree <- SNAG_meta_perma %>%
  select(natman, sample, aspect, ds_at_drill, reads, log_reads) %>%
  mutate(ds_num = as.numeric(as.character(ds_at_drill))) %>%
  group_by(natman) %>%
  summarise(
    n_samples = n(),
    ds_N      = ds_num[aspect == "NORTH"][1],
    ds_S      = ds_num[aspect == "SOUTH"][1],
    delta_ds  = ds_N - ds_S,
    abs_delta = abs(delta_ds),
    .groups   = "drop"
  )
cat("\nΔds per Snag (NORTH - SOUTH) & counts by |Δds|:\n"); print(ds_by_tree %>% select(natman, ds_N, ds_S, delta_ds, abs_delta))
print(table(ds_by_tree$abs_delta, useNA = "ifany"))

# 10) snag-level ds effects (among snags) -------------------------------------
otu_tree <- rowsum(otu_snag_perma, group = SNAG_meta_perma$natman)

ds_num_per_core <- as.numeric(as.character(SNAG_meta_perma$ds_at_drill))
meta_tree_num <- SNAG_meta_perma %>%
  mutate(ds_num = ds_num_per_core) %>%
  group_by(natman) %>%
  summarise(
    ds_num_tree = mean(ds_num, na.rm = TRUE),
    reads       = sum(reads),
    log_reads   = log1p(reads),
    .groups     = "drop"
  )

meta_tree_cat <- SNAG_meta_perma %>%
  group_by(natman) %>%
  summarise(
    ds_levels = list(unique(ds_at_drill[!is.na(ds_at_drill)])),
    reads     = sum(reads),
    log_reads = log1p(reads),
    .groups   = "drop"
  ) %>%
  mutate(
    ds_tree_cat = ifelse(lengths(ds_levels) == 1,
                         purrr::map_chr(ds_levels, ~ as.character(.x[1])),
                         "mixed"),
    ds_tree_cat = factor(ds_tree_cat, levels = c("0","1","2","3","4","5","mixed"))
  ) %>% select(-ds_levels)

D_tree_robAit    <- vegdist(otu_tree, method = "robust.aitchison")
comm_tree_rel    <- decostand(otu_tree, method = "total")
D_tree_bray      <- vegdist(comm_tree_rel, method = "bray")
D_tree_braySqrt  <- sqrt(D_tree_bray)
D_tree_jacc      <- vegdist((otu_tree > 0) * 1, method = "jaccard", binary = TRUE)

cat("\nSnag-level PERMANOVA (numeric ds mean across aspects):\n")
perma_tree_num <- list(
  robAit   = adonis2(D_tree_robAit   ~ log_reads + ds_num_tree, data = meta_tree_num, by = "margin", permutations = how(nperm = nperm)),
  bray     = adonis2(D_tree_bray     ~ log_reads + ds_num_tree, data = meta_tree_num, by = "margin", permutations = how(nperm = nperm)),
  braySqrt = adonis2(D_tree_braySqrt ~ log_reads + ds_num_tree, data = meta_tree_num, by = "margin", permutations = how(nperm = nperm)),
  jaccard  = adonis2(D_tree_jacc     ~ ds_num_tree,              data = meta_tree_num, by = "margin", permutations = how(nperm = nperm))
)
print(perma_tree_num$robAit); 
print(perma_tree_num$bray)
print(perma_tree_num$braySqrt);
print(perma_tree_num$jaccard)


# 11) Indicator analyses (scaled; BH FDR) -------------------------------------
otu_scaled <- rescale_to_target(otu_snag)
colnames(otu_scaled) <- tibble(sh_code = colnames(otu_scaled)) %>%
  left_join(select(tax, sh_code, lowest_taxon), by="sh_code") %>%
  transmute(name = if_else(!is.na(lowest_taxon), paste0(sh_code,"|",lowest_taxon), sh_code)) %>%
  pull(name)
otu_scaled <- otu_scaled[, colSums(otu_scaled) > 0, drop = FALSE]
otu_pa <- (otu_snag > 0) * 1
colnames(otu_pa) <- colnames(otu_scaled)
otu_pa <- otu_pa[, colSums(otu_pa) > 0, drop = FALSE]

run_indic <- function(X, grp, func, lab) {
  if (!is.factor(grp) || nlevels(grp) < 2L) return(NULL)
  res <- indicspecies::multipatt(X, droplevels(grp), func = func, control = how(nperm=999), duleg = TRUE)
  as_tibble(res$sign, rownames="feature") %>%
    mutate(p_adj = p.adjust(p.value, "BH"), group = lab, method = func) %>%
    arrange(p_adj)
}

indic_tabs <- list(
  run_indic(otu_scaled, SNAG_meta$aspect,      "r.g",     "aspect"),
  run_indic(otu_scaled, SNAG_meta$ds_at_drill, "r.g",     "ds_at_drill"),
  run_indic(otu_scaled, SNAG_meta$natman,      "r.g",     "natman"),
  run_indic(otu_pa,     SNAG_meta$aspect,      "IndVal.g","aspect_PA"),
  run_indic(otu_pa,     SNAG_meta$ds_at_drill, "IndVal.g","ds_at_drill_PA"),
  run_indic(otu_pa,     SNAG_meta$natman,      "IndVal.g","natman_PA")
) %>% purrr::compact() %>% bind_rows()

cat("\nIndicator results (top 10 per test; BH<=0.20 shown):\n")
indic_tabs %>%
  group_by(group, method) %>%
  filter(p.value <= 0.002) %>%
  slice_min(order_by = stat, n = 10, with_ties = FALSE) %>%
  ungroup() %>% print(n = Inf)

# 12) Alpha richness GLMM — model selection + EMMs + IRR ----------------------
df_snag <- alpha_from_otu(otu_snag) %>%
  left_join(SNAG_meta %>% select(sample, natman, aspect, ds_at_drill), by = "sample") %>%
  mutate(across(c(natman, aspect, ds_at_drill), droplevels)) %>%
  filter(is.finite(richness), reads > 0)

if (nlevels(df_snag$aspect) > 1 || nlevels(df_snag$ds_at_drill) > 1) {
  fit_pois <- glmmTMB(richness ~ aspect + ds_at_drill + log_reads + (1|natman),
                      data = df_snag, family = poisson)
  rp  <- residuals(fit_pois, type = "pearson"); rdf <- df.residual(fit_pois)
  phi_p <- sqrt(sum(rp^2) / rdf)
  use_nb <- is.finite(phi_p) && (phi_p > 1.2)
  if (use_nb) {
    fit_nb <- glmmTMB(richness ~ aspect + ds_at_drill + log_reads + (1|natman),
                      data = df_snag, family = nbinom2,
                      control = glmmTMBControl(optCtrl = list(iter.max=1e4, eval.max=1e4)))
    best <- if (AIC(fit_pois) + 2 < AIC(fit_nb)) fit_pois else fit_nb
  } else best <- fit_pois
  best_type <- if (identical(best$family$family, poisson()$family)) "Poisson (glmmTMB)" else "NB (glmmTMB)"
  cat(sprintf("\nAlpha richness GLMM — selected model: %s | Poisson phi=%.2f\n", best_type, phi_p))
  # ---- Aspect × Decay interaction probe (nested LRT vs. additive "best") ----
  pick_family <- function(best){
    if ("glmmTMB" %in% class(best)) {
      f <- tolower(best$modelInfo$family$family)
    } else {
      f <- tolower(best$family$family)
    }
    if (f %in% c("poisson"))    return(poisson())
    if (f %in% c("nbinom2","negative binomial 2","nb2")) return(nbinom2())
    stop("Unsupported family in `best` (", f, ").")
  }
  

  fam_fun <- pick_family(best)
  
  # Fit interaction model with same structure + family as `best`
  fit_int <- try(
    glmmTMB(
      richness ~ aspect * ds_at_drill + log_reads + (1|natman),
      data = df_snag,
      family = fam_fun,
      control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))
    ),
    silent = TRUE
  )
  
  if (inherits(fit_int, "try-error")) {
    cat("\n[Interaction probe] glmmTMB failed to fit the aspect×decay model; keeping additive model.\n")
    lrt_int <- NULL
  } else {
    # Likelihood-ratio test (nested)
    lrt_int <- try(anova(best, fit_int), silent = TRUE)
    if (inherits(lrt_int, "try-error")) {
      cat("\n[Interaction probe] anova() failed; keeping additive model.\n")
      lrt_int <- NULL
    } else {
      cat("\n[Interaction probe] Aspect × Decay LRT (additive vs interaction):\n")
      print(lrt_int)
      
      # Save LRT table
      readr::write_csv(
        tibble::as_tibble(lrt_int, rownames = "model"),
        file.path(outdir, "alpha_glmm_aspectXdecay_LRT.csv")
      )
      
      # If interaction is significant, add a simple EMM interaction plot (facets = ds)
      p_int <- NA_real_
      if (nrow(lrt_int) >= 2 && "Pr(>Chisq)" %in% colnames(lrt_int)) {
        p_int <- suppressWarnings(as.numeric(lrt_int$`Pr(>Chisq)`[2]))
      }
      if (is.finite(p_int) && p_int < 0.05) {
        cat(sprintf("[Interaction probe] aspect×decay is significant (LRT p = %.4g). Writing EMM plot.\n", p_int))
        emm_int <- emmeans::emmeans(fit_int, ~ aspect | ds_at_drill, type = "response")
        df_emm  <- as.data.frame(emm_int)
        # Standardize CI column names if asymptotic naming is used
        if ("asymp.LCL" %in% names(df_emm)) df_emm <- dplyr::rename(df_emm, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
        
        p <- ggplot2::ggplot(
          df_emm,
          ggplot2::aes(x = aspect, y = response, ymin = lower.CL, ymax = upper.CL)
        ) +
          ggplot2::geom_point(size = 2.8) +
          ggplot2::geom_errorbar(width = 0.18) +
          ggplot2::facet_wrap(~ ds_at_drill, scales = "free_y") +
          ggplot2::labs(
            x = "Aspect",
            y = "Estimated richness",
            title = paste0("Richness ~ aspect × decay (EMM, ", best_type, ")"),
            subtitle = sprintf("LRT p = %.4g", p_int)
          ) +
          ggplot2::theme_classic(base_size = 11)
        
        ggplot2::ggsave(file.path(outdir, "alpha_emm_interaction_aspect_by_decay.png"),
                        p, width = 6.5, height = 3.8, dpi = 300)
        print(p)
      } else if (is.finite(p_int)) {
        cat(sprintf("[Interaction probe] No evidence for aspect×decay (LRT p = %.4g). Keeping additive model.\n", p_int))
      } else {
        cat("[Interaction probe] Could not compute LRT p-value. Keeping additive model.\n")
      }
    }
  }
  

  lrt_tab <- drop1(best, test="Chisq"); r2_tab <- performance::r2(best)
  cat("\nAlpha richness GLMM — Type-II LRT and R²:\n"); print(lrt_tab); print(r2_tab)
  
  # IRR table
  irr_tab <- broom.mixed::tidy(best, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    transmute(term, IRR = estimate, CI_low = conf.low, CI_high = conf.high, p_value = p.value) %>%
    arrange(term)
  cat("\nIncidence Rate Ratios (IRR):\n"); print(irr_tab)
  
  # EMMs for aspect and ds
  emm_as  <- if (nlevels(df_snag$aspect) > 1) emmeans(best, ~ aspect, type = "response") else NULL
  emm_ds  <- if (nlevels(df_snag$ds_at_drill) > 1) emmeans(best, ~ ds_at_drill, type = "response") else NULL
  
  plot_emm <- function(emm, xlab, fname) {
    if (is.null(emm)) return(invisible(NULL))
    df <- as.data.frame(emm)
    if ("asymp.LCL" %in% names(df)) df <- rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
    x_var <- names(df)[1]
    p <- ggplot(df, aes(x = !!sym(x_var), y = response, ymin = lower.CL, ymax = upper.CL)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      labs(x = xlab, y = "Estimated richness", title = paste("Richness ~", xlab, "(EMM,", best_type, ")")) +
      theme_classic(base_size = 11)
    ggsave(file.path(outdir, fname), p, width = 5, height = 3.4, dpi = 300); print(p)
  }
  plot_emm(emm_as, "Aspect",        "alpha_emm_aspect.png")
  plot_emm(emm_ds, "ds_at_drill",   "alpha_emm_ds.png")
  
  # Compact publication tables
  publist <- list(model = summary(best), lrt = lrt_tab, r2 = r2_tab, irr = irr_tab)
  saveRDS(publist, file.path(outdir, "alpha_glmm_publication_tables.rds"))
  readr::write_csv(irr_tab, file.path(outdir, "alpha_glmm_IRR.csv"))
  readr::write_csv(as_tibble(lrt_tab, rownames = "term"), file.path(outdir, "alpha_glmm_LRT.csv"))
  cat("\nSaved publication tables to:", outdir, "\n")
}

# 13) Leave-one-snag-out (LOO) stability (within-snag PERMANOVA terms) --------
trees <- SNAG_meta_perma$natman %>% as.character() %>% unique() %>% setdiff(NA)
loo <- purrr::map_dfr(trees, function(t) {
  keep <- !is.na(SNAG_meta_perma$natman) & SNAG_meta_perma$natman != t
  if (sum(keep) < 5) return(NULL)
  Dk <- vegdist(otu_snag_perma[keep, , drop = FALSE], method = "robust.aitchison")
  permk <- how(nperm = 499); setBlocks(permk) <- droplevels(factor(SNAG_meta_perma$natman[keep]))
  m <- adonis2(Dk ~ log_reads + aspect + ds_at_drill,
               data = SNAG_meta_perma[keep, ], by = "margin", permutations = permk)
  as_tibble(as.data.frame(m), rownames = "term") %>%
    filter(term %in% c("aspect","ds_at_drill")) %>%
    transmute(left_out = t, term, R2, p = `Pr(>F)`)
})
if (nrow(loo) > 0) {
  cat("\nLOO stability — within-snag PERMANOVA terms (Aitchison):\n")
  print(loo %>% group_by(term) %>%
          summarise(median_R2 = median(R2, na.rm = TRUE),
                    p05 = quantile(p, 0.05, na.rm = TRUE),
                    p95 = quantile(p, 0.95, na.rm = TRUE),
                    .groups = "drop"))
  readr::write_csv(loo, file.path(outdir, "loo_permanova_terms.csv"))
}

cat("\n=== Script 7 (SNAG) — done. Outputs in:", outdir, "===\n")
