# =============================================================================
# Script 3 — Alpha + Beta diversity (LOG samples)
#   - Alpha: richness, Shannon, Simpson; boxplots by key factors + GAM vs diameter
#   - Beta: Bray & Jaccard distances; PCoA (wcmdscale add=TRUE) & NMDS;
#           PERMANOVA (adonis2) with optional strata=natman; beta-dispersion tests
#   - Exports: plots -> plots/LOG; tables -> tables/LOG
#   - Requires objects from Script 1: META1, tax, otu_matrix_filt, color palettes
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(rstatix)
  library(tidyverse)
  library(vegan)
  library(iNEXT)
  library(janitor)
  library(indicspecies)
  library(eulerr)
  library(scales)
})

set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# Output dirs
dir.create("plots/LOG",   recursive = TRUE, showWarnings = FALSE)
dir.create("tables/LOG",  recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 0) Build LOG subset + factors
# -----------------------------------------------------------------------------
# Expecting META1, otu_matrix_filt, and color palettes from Script 1
stopifnot(exists("META1"), exists("otu_matrix_filt"))

log_meta <- META1 %>%
  filter(dw_type == "LOG", !is.na(dw_type), ds_at_drill != "0") %>%
  mutate(
    ds_at_drill = factor(ds_at_drill, levels = c("0","1","2","3","4","5")),
    decay_stage = factor(decay_stage, levels = c("LIVING","EARLY","AVERAGE","LATE")))
otu_log <- otu_matrix_filt[log_meta$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_log), log_meta$sample))
cat("for dw_type LOG, looking at samples: ", nrow(otu_log))
# -----------------------------------------------------------------------------
# 1) Alpha diversity
# -----------------------------------------------------------------------------
alpha_log <- tibble(
  sample   = rownames(otu_log),
  richness = vegan::specnumber(otu_log),
  shannon  = vegan::diversity(otu_log, index = "shannon"),
  simpson  = vegan::diversity(otu_log, index = "simpson")
) %>%
  left_join(log_meta, by = "sample")

alpha_long <- alpha_log %>%
  pivot_longer(cols = c(richness, shannon, simpson),
               names_to = "metric", values_to = "value")

# --- Boxplots by decay_stage
p_box_decay <- alpha_long %>%
  filter(!is.na(decay_stage)) %>%
  ggplot(aes(x = decay_stage, y = value, fill = decay_stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = DS_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by Decay Stage (LOG samples)",
       x = "Decay stage", y = "Value", fill = "Decay stage")
p_box_decay
ggsave("plots/LOG/boxplot_alpha_by_decay_stage.png", p_box_decay, width = 10, height = 5, dpi = 300)

# --- Boxplots by ds_at_drill
p_box_ds_drill <- alpha_long %>%
  filter(!is.na(ds_at_drill)) %>%
  ggplot(aes(x = ds_at_drill, y = value, fill = ds_at_drill)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = ds_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by Decay at Drill (LOG samples)",
       x = "Decay at drill", y = "Value", fill = "Decay at drill")
p_box_ds_drill
ggsave("plots/LOG/boxplot_alpha_by_ds_at_drill.png", p_box_ds_drill, width = 10, height = 5, dpi = 300)

# --- Boxplots by size
p_box_size <- alpha_long %>%
  filter(!is.na(size)) %>%
  ggplot(aes(x = size, y = value, fill = size)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = sc_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by Size Class (LOG samples)",
       x = "Size class", y = "Value", fill = "Size class")
p_box_size
ggsave("plots/LOG/boxplot_alpha_by_size.png", p_box_size, width = 10, height = 5, dpi = 300)

# --- GAM vs diameter_at_drill (bands mark size-class breaks)
alpha_long_cont <- alpha_long %>% filter(!is.na(diameter_at_drill))
band_df <- tibble(
  xmin  = c(-Inf, 9,19,39, 69, 99),
  xmax  = c(9,19,39, 69, 99, Inf),
  label = c("very fine","fine","small","medium","large","very large")
)

p_alpha_diam_gam <- ggplot(alpha_long_cont, aes(x = diameter_at_drill, y = value)) +
  geom_rect(data = band_df,
            aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
            inherit.aes = FALSE, fill = "grey80", alpha = 0.08) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = TRUE) +
  geom_vline(xintercept = c(39, 69, 99), linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(title = "Alpha diversity vs. diameter at drill (LOG samples)",
       x = "Diameter at drill (mm)", y = "Alpha diversity") +
  theme_minimal(base_size = 13)
p_alpha_diam_gam
ggsave("plots/LOG/alpha_vs_diameter_gam.png", p_alpha_diam_gam, width = 11, height = 4, dpi = 300)

p_alpha_diam_byclass <- ggplot(alpha_long_cont,
                               aes(x = diameter_at_drill, y = value, color = size)) +
  geom_point(alpha = 0.45, size = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = FALSE) +
  geom_vline(xintercept = c(39, 69, 99), linetype = "dashed", linewidth = 0.4, color = "grey40") +
  scale_color_manual(values = sc_colors) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(title = "Alpha diversity vs. diameter at drill (colored by size class)",
       x = "Diameter at drill (mm)", y = "Alpha diversity", color = "Size class") +
  theme_minimal(base_size = 13)
p_alpha_diam_byclass
ggsave("plots/LOG/alpha_vs_diameter_by_sizeclass.png", p_alpha_diam_byclass, width = 11, height = 4, dpi = 300)

# --- Boxplots by depth
p_box_depth <- alpha_long %>%
  mutate(depth_2 = factor(depth_2, levels = c("OUTER", "INNER"))) %>% 
  filter(!is.na(depth_2), depth_2 != "NO_DIST") %>%
  ggplot(aes(x = (depth_2), y = value, fill = depth_2)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = depth_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by Depth (LOG samples)",
       x = "Depth", y = "Value", fill = "Depth")
p_box_depth
ggsave("plots/LOG/boxplot_alpha_by_depth.png", p_box_depth, width = 10, height = 5, dpi = 300)

# --- Boxplots by position
p_box_pos <- alpha_long %>%
  filter(!is.na(position), !position %in% c("ENDO")) %>%
  ggplot(aes(x = position, y = value, fill = position)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  scale_fill_manual(values = position_colors) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by Position (LOG samples)",
       x = "Position", y = "Value", fill = "Position")
p_box_pos
ggsave("plots/LOG/boxplot_alpha_by_position.png", p_box_pos, width = 10, height = 5, dpi = 300)

# --- Boxplots by natman (tree ID), restricting to INTRA rq
dat_tree <- alpha_long %>% filter(!is.na(natman))
if ("rq" %in% names(dat_tree)) dat_tree <- dat_tree %>% filter(rq == "INTRA")
p_box_tree <- dat_tree %>%
  ggplot(aes(x = natman, y = value, fill = natman)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Alpha Diversity by INTRA-tree (LOG samples)",
       x = "tree", y = "Value", fill = "tree ID")
p_box_tree
ggsave("plots/LOG/boxplot_alpha_by_tree.png", p_box_tree, width = 10, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 2) Top-N taxa per group (relative abundance within samples)
# -----------------------------------------------------------------------------
otu_long <- otu_log %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  left_join(select(log_meta, sample, decay_stage, ds_at_drill, position, depth_2, size), by = "sample") %>%
  filter(!is.na(decay_stage), !is.na(ds_at_drill), !is.na(position), !is.na(depth_2), !is.na(size)) %>%
  filter(!position %in% c("AFWD1","AFWD2","CROWN"), depth_2 != "NO_DIST") %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = abundance / sum(abundance)) %>%
  ungroup() %>%
  left_join(select(tax, sh_code, species, genus), by = "sh_code")

top_N      <- 5
group_vars <- c("decay_stage","ds_at_drill","position","depth_2","size")
top_species_list <- list()

for (var in group_vars) {
  message("Processing ", var, " ...")
  top_tab <- otu_long %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(across(all_of(c(var, "sh_code", "species", "genus")))) %>%
    summarise(
      mean_rel_abund_sample = mean(rel_abund_sample, na.rm = TRUE),
      n_samples             = n(),
      .groups = "drop"
    ) %>%
    group_by(across(all_of(var))) %>%
    arrange(desc(mean_rel_abund_sample), .by_group = TRUE) %>%
    slice_head(n = top_N) %>%
    arrange(across(all_of(var)), desc(mean_rel_abund_sample))
  
  top_species_list[[var]] <- top_tab
#  readr::write_tsv(top_tab, paste0("tables/LOG/top_", var, "_N", top_N,"_threshold", threshold, ".tsv"))
  print(top_tab, n = top_N * length(unique(otu_long[[var]])))
}

# -----------------------------------------------------------------------------
# 3) Beta diversity — distances, ordinations, PERMANOVA, dispersion
# -----------------------------------------------------------------------------
# Distances
comm_rel  <- decostand(otu_log, method = "total")   # per-sample proportions
bray_dist <- vegdist(comm_rel, method = "bray")
jac_dist  <- vegdist((otu_log > 0) * 1, method = "jaccard", binary = TRUE)
ra_dist  <- vegdist(otu_log, method = "robust.aitchison")
attr(ra_dist, "Labels") <- rownames(otu_log)
# ---- PCoA via wcmdscale (add=TRUE to handle non-Euclidean)
pcoa_bray <- wcmdscale(bray_dist, k = 2, eig = TRUE, add = TRUE)
pcoa_jac  <- wcmdscale(jac_dist,  k = 2, eig = TRUE, add = TRUE)
str(pcoa_bray)
scores_bray <- as_tibble(pcoa_bray$points, .name_repair = "unique") %>%
  mutate(sample = rownames(pcoa_bray$points)) %>%
  rename(PCo1 = Dim1, PCo2 = Dim2) %>%
  left_join(log_meta, by = "sample")

scores_jac <- as_tibble(pcoa_jac$points, .name_repair = "unique") %>%
  mutate(sample = rownames(pcoa_jac$points)) %>%
  rename(PCo1 = Dim1, PCo2 = Dim2) %>%
  left_join(log_meta, by = "sample")
# Bray–Curtis percentages
pc1_bray <- pcoa_bray$eig[1] / sum(pcoa_bray$eig[pcoa_bray$eig > 0]) * 100
pc2_bray <- pcoa_bray$eig[2] / sum(pcoa_bray$eig[pcoa_bray$eig > 0]) * 100

# Jaccard percentages
pc1_jac <- pcoa_jac$eig[1] / sum(pcoa_jac$eig[pcoa_jac$eig > 0]) * 100
pc2_jac <- pcoa_jac$eig[2]  / sum(pcoa_jac$eig[pcoa_jac$eig > 0]) * 100
# readr::write_tsv(scores_bray, "tables/LOG/pcoa_bray_scores.tsv")
# readr::write_tsv(scores_jac,  "tables/LOG/pcoa_jaccard_scores.tsv")

# PCoA plots colored by decay stage
p_pcoa_bray <- ggplot(scores_bray, aes(PCo1, PCo2, color = decay_stage, shape = depth_2)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = DS_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "PCoA (Bray–Curtis) — LOG samples",
       color = "Decay", shape = "Depth", 
       x = sprintf("PCo1 (%.1f%%)", pc1_bray),
       y = sprintf("PCo2 (%.1f%%)", pc2_bray)
  )
p_pcoa_bray
ggsave("plots/LOG/pcoa_bray_decay_depth.png", p_pcoa_bray, width = 6.8, height = 5.4, dpi = 300)

p_pcoa_jac <- ggplot(scores_jac, aes(PCo1, PCo2, color = decay_stage, shape = depth_2)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = DS_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "PCoA (Jaccard) — LOG samples",
       x = sprintf("PCo1 (%.1f%%)", pc1_jac),
       y = sprintf("PCo2 (%.1f%%)", pc2_jac) ,
       color = "Decay", shape = "Depth")
p_pcoa_jac
ggsave("plots/LOG/pcoa_jaccard_decay_depth.png", p_pcoa_jac, width = 6.8, height = 5.4, dpi = 300)

# ---- NMDS (Bray)
set.seed(42)
nmds_bray <- metaMDS(comm_rel, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
mds_scores_bray <- as_tibble(scores(nmds_bray, display = "sites")) %>%
  mutate(sample = rownames(otu_log)) %>%
  left_join(log_meta, by = "sample")

p_nmds_bray <- ggplot(mds_scores_bray, aes(NMDS1, NMDS2, color = decay_stage, shape = depth_2)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = DS_colors) +
  theme_minimal(base_size = 13) +
  labs(title = paste0("NMDS of LOG samples (Bray) — stress = ", round(nmds_bray$stress, 3)),
       x = "NMDS1", y = "NMDS2", color = "Decay", shape = "Depth")
p_nmds_bray
ggsave("plots/LOG/nmds_bray_decay_depth.png", p_nmds_bray, width = 6.8, height = 5.4, dpi = 300)
# ---- NMDS (robust.aitchison)
nmds_ra <- metaMDS(otu_log, distance = "robust.aitchison", k = 2, trymax = 100, autotransform = FALSE)
mds_scores_ra <- as_tibble(scores(nmds_ra, display = "sites")) %>%
  mutate(sample = rownames(otu_log)) %>%
  left_join(log_meta, by = "sample")

p_nmds_ra <- ggplot(mds_scores_ra, aes(NMDS1, NMDS2, color = decay_stage, shape = depth_2)) +
  geom_point(alpha = 0.9, size = 2) +
  scale_color_manual(values = DS_colors) +
  theme_minimal(base_size = 13) +
  labs(title = paste0("NMDS of LOG samples (r.Aitchison) — stress = ", round(nmds_ra$stress, 3)),
       x = "NMDS1", y = "NMDS2", color = "Decay", shape = "Depth")
p_nmds_ra
ggsave("plots/LOG/nmds_ra_decay_depth.png", p_nmds_ra, width = 6.8, height = 5.4, dpi = 300)

# ---- PERMANOVA across key factors (Bray & Jaccard)
group_vars <- c("decay_stage","ds_at_drill","position","depth_2","size")

perm_list <- list()

  for (dist_name in c("bray","jaccard","robust.aitchison")) {
    if (dist_name == "bray") {
      D <- bray_dist
    } else if (dist_name == "jaccard") {
      D <- jac_dist
    } else if (dist_name == "robust.aitchison") {
      D <- ra_dist
    }
    
    for (var in group_vars) {
      dat <- log_meta %>%
        dplyr::select(sample, all_of(var), natman) %>%
        dplyr::filter(!is.na(.data[[var]]))
      
      # Use intersection of sample names
      common_samples <- intersect(labels(D), dat$sample)
      if (length(common_samples) < 2) next
      
      # Subset both D and dat consistently
      D_sub <- as.dist(as.matrix(D)[common_samples, common_samples])
      dat   <- dat %>% dplyr::filter(sample %in% common_samples)
      dat   <- dat[match(common_samples, dat$sample), ]  # align order
      
      strata_vec <- if ("natman" %in% names(dat)) dat$natman else NULL
      if (!is.null(strata_vec) && (dplyr::n_distinct(na.omit(strata_vec)) < 2)) strata_vec <- NULL
      
      # Guard: skip if grouping var has <2 levels
      if (dplyr::n_distinct(na.omit(dat[[var]])) < 2) next
      
      set.seed(42)
      ares <- adonis2(D_sub ~ dat[[var]], permutations = 999, strata = strata_vec)
      
      res_tbl <- as.data.frame(ares) %>%
        tibble::rownames_to_column("term") %>%
        dplyr::filter(term != "Residual") %>%
        dplyr::rename(
          df = Df,
          sumsq = SumOfSqs,
          r.squared = R2,
          statistic = F,
          p.value = `Pr(>F)`
        ) %>%
        dplyr::mutate(
          distance = dist_name,
          variable = var
        ) %>%
        dplyr::select(distance, variable, df, sumsq, r.squared, statistic, p.value)
      
      perm_list[[paste(dist_name, var, sep = "_")]] <- res_tbl
      print(res_tbl)
    }
  }

permanova_results <- bind_rows(perm_list)
permanova_results
#readr::write_tsv(permanova_results, "tables/LOG/general_permanova_results.tsv")

# ---- Beta-dispersion (homogeneity of dispersion) for each factor (Bray)
disp_list <- list()
for (var in group_vars) {
  dat <- log_meta %>% select(sample, all_of(var)) %>% filter(!is.na(.data[[var]]))
  sub_samples <- dat$sample
  if (length(sub_samples) < 3 || n_distinct(dat[[var]]) < 2) next
  
  D_sub <- as.dist(as.matrix(bray_dist)[sub_samples, sub_samples])
  grp   <- droplevels(factor(dat[[var]]))
  
  bd <- betadisper(D_sub, group = grp, type = "centroid")
  set.seed(42)
  pt <- permutest(bd, permutations = 999)
  
  disp_df <- tibble(
    variable = var,
    F        = unname(pt$tab[1, "F"]),
    p_value  = unname(pt$tab[1, "Pr(>F)"]),
    df       = paste0(unname(pt$tab[1, "Df"]), ",", unname(pt$tab[2, "Df"]))
  )
  disp_list[[var]] <- disp_df
  
  # Simple plot of distances-to-centroid
  disp_plot_df <- tibble(group = grp, dist = bd$distances)
  p_disp <- ggplot(disp_plot_df, aes(x = group, y = dist, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_blank()) +
    labs(title = paste0("Dispersion (Bray) — ", var),
         x = var, y = "Distance to centroid", fill = var)
  print(p_disp)
#  ggsave(paste0("plots/LOG/dispersion_bray_", var, ".png"), p_disp, width = 6.8, height = 5.4, dpi = 300)
}
dispersion_results_bray <- bind_rows(disp_list)
readr::write_tsv(dispersion_results_bray, "tables/LOG/general_dispersion_results_bray.tsv")
dispersion_results_bray
# position and size are overdispersed!



