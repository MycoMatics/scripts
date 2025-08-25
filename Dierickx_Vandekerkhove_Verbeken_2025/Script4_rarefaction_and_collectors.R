# =============================================================================
# Script 4 — Rarefaction & Collector's Curves (incl. iNEXT)
#   - Sample rarefaction per group
#   - Collector's (specaccum) curves by deadwood types
#   - iNEXT (incidence-frequency) for DW types, FWD split (aFWD/fFWD),
#     WOODY vs SOIL, and within-LOG contrasts (depth, position, ds_at_drill,
#     decay_stage, size)
#   - Uses objects from Script 1: META1, otu_matrix_filt, palettes, threshold
#   - Headless plotting only via ggsave() or pdf(NULL); dev.off()
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(iNEXT)
  library(janitor)
})

set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# -----------------------------------------------------------------------------
# 0) Guards, output dirs, palettes from Script 1
# -----------------------------------------------------------------------------
stopifnot(exists("META1"), exists("otu_matrix_filt"))
stopifnot(exists("threshold"))
cat("You are using SH threshold ", threshold, "%")
dir.create("plots", recursive = TRUE, showWarnings = FALSE)
dir.create("plots/LOG", recursive = TRUE, showWarnings = FALSE)

# Keep a working copy synced to META1 (should already be aligned by Script 1)
otu_work <- otu_matrix_filt[rownames(otu_matrix_filt) %in% META1$sample, , drop = FALSE]
META1    <- META1 %>% filter(sample %in% rownames(otu_work))
# -----------------------------------------------------------------------------
# 1) DW types (LOG/FWD/SNAG/SOIL) — rarefaction + specaccum + iNEXT
# -----------------------------------------------------------------------------
dw_keep <- c("SNAG","LOG","FWD","SOIL")  # adjust if needed
META_rare <- META1 %>% filter(dw_type %in% dw_keep) %>% 
  filter(sample %in% rownames(otu_work))

otu_work <- otu_work[rownames(otu_work) %in% META_rare$sample, , drop = FALSE]
otu_work <- otu_work[rowSums(otu_work) > 0, , drop = FALSE]
stopifnot(identical(rownames(otu_work), META_rare$sample))
# Sample rarefaction (all included samples)
# Remove samples with zero reads before rarefaction

sample_names <- rownames(otu_work)
pdf(NULL)  # prevent interactive plotting on headless HPC
rarecurve_list <- vegan::rarecurve(otu_work, step = 1500, label = FALSE, return = "list")
dev.off()

rare_df <- purrr::map2_dfr(
  rarecurve_list, sample_names,
  ~tibble(sample = .y,
          reads   = attr(.x, "Subsample"),
          richness= as.integer(.x))
) %>%
  left_join(META1 %>% select(sample, dw_type), by = "sample")

p_rare <- ggplot(rare_df, aes(x = reads, y = richness, group = sample, color = dw_type)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = dw_colors, na.value = "grey70") +
  theme_minimal(base_size = 13) +
  labs(title = "Rarefaction curves by deadwood type",
       x = "Reads (subsampled)", y = "Observed SH richness", color = "Deadwood type")
p_rare
ggsave(paste0("plots/rarefaction_by_dw_type_filtered_", threshold, ".png"),
       p_rare, width = 8, height = 5, dpi = 300)

# Show lowest-read samples (diagnostic)
rarefied_max <- rare_df %>%
  group_by(sample, dw_type) %>%
  summarise(max_reads = max(reads, na.rm = TRUE),
            max_richness = max(richness, na.rm = TRUE), .groups = "drop")
lowest_20_samples <- rarefied_max %>% arrange(max_reads) %>% slice_head(n = 20) %>% pull(sample)
rare_df_lowest <- rare_df %>% filter(sample %in% lowest_20_samples)
rare_df_lowest$sample <- factor(rare_df_lowest$sample, levels = lowest_20_samples)
# show samples with low read counts, the lowest are plateauing
rarefied_max %>% 
  filter(sample %in% rare_df_lowest$sample) %>% 
  arrange(max_reads)

p_rare_lowest <- ggplot(rare_df_lowest, aes(x = reads, y = richness, color = dw_type, group = sample)) +
  geom_line(alpha = 0.9, linewidth = 1) +
  scale_color_manual(values = dw_colors, na.value = "grey70") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom") +
  labs(title = "Rarefaction: 20 samples with lowest read counts",
       x = "Reads (subsampled)", y = "Observed SH richness", color = "Deadwood type")
p_rare_lowest
ggsave(paste0("plots/rarefaction_lowest_nsamples_", threshold, ".png"),
       p_rare_lowest, width = 8, height = 5, dpi = 300)

# Ensure sync before group analyses
otu_work <- otu_work[rownames(otu_work) %in% META1$sample, , drop = FALSE]
META1    <- META1 %>% filter(sample %in% rownames(otu_work))

# Collector's curves (specaccum) per dw_type
collector_curves <- lapply(dw_keep, function(dw) {
  samples <- META1 %>% filter(dw_type == dw) %>% pull(sample)
  if (length(samples) < 2) return(NULL)   # need >=2 for random perm accumulation
  otus <- otu_work[rownames(otu_work) %in% samples, , drop = FALSE]
  otus <- otus[, colSums(otus) > 0, drop = FALSE]
  acc <- specaccum(otus, method = "random", permutations = 100)
  tibble(samples = acc$sites, richness = acc$richness, sd = acc$sd, dw_type = dw)
})
collector_df <- bind_rows(collector_curves)

p_collectors <- ggplot(collector_df, aes(x = samples, y = richness, color = dw_type, fill = dw_type)) +
  geom_line(linewidth = 1.1, alpha = 0.85) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.2, color = NA) +
  scale_color_manual(values = dw_colors) +
  scale_fill_manual(values = dw_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "Collector’s curves (specaccum) by deadwood type",
       x = "Number of samples", y = "Cumulative SH richness",
       color = "Deadwood type", fill = "Deadwood type")
p_collectors
ggsave(paste0("plots/collectors_specaccum_by_dw_type_", threshold, ".png"),
       p_collectors, width = 8, height = 5, dpi = 300)

# iNEXT incidence-frequency input per dw_type
incidence_list <- lapply(dw_keep, function(dw) {
  samples <- META1 %>% filter(dw_type == dw) %>% pull(sample)
  if (length(samples) < 2) return(NULL)
  otus <- otu_work[rownames(otu_work) %in% samples, , drop = FALSE]
  otus <- otus[, colSums(otus) > 0, drop = FALSE]
  freq <- colSums(otus > 0)
  c(length(samples), freq)
})
names(incidence_list) <- dw_keep
incidence_list <- incidence_list[!sapply(incidence_list, is.null)]

endpoints <- sapply(names(incidence_list), function(dw) {
  n_samples <- incidence_list[[dw]][1]
  min(ceiling(n_samples * 2), 200)
})
iNEXT_out <- lapply(names(incidence_list), function(dw) {
  out <- iNEXT(incidence_list[[dw]], q = c(0,1,2), datatype = "incidence_freq", endpoint = endpoints[dw])
  as_tibble(out$iNextEst$size_based) %>% mutate(dw_type = dw)
})
iNEXT_df <- bind_rows(iNEXT_out)
observed_n <- sapply(names(incidence_list), function(dw) incidence_list[[dw]][1])
iNEXT_df <- iNEXT_df %>% mutate(interp_extrap = ifelse(t <= observed_n[dw_type], "interpolated", "extrapolated"))

order_labs <- c("0"="Richness (q = 0)\n(Number of SHs)",
                "1"="Exp Shannon (q = 1)\n(Number of typical SHs)",
                "2"="Inverse Simpson (q = 2)\n(Number of dominant SHs)")
interp_points <- iNEXT_df %>%
  group_by(dw_type, Order.q) %>%
  filter(t == max(t[interp_extrap == "interpolated"])) %>%
  ungroup()

p_inext_dw <- ggplot(iNEXT_df, aes(x = t, y = qD, color = dw_type, fill = dw_type, linetype = interp_extrap)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
  geom_point(data = interp_points, aes(x = t, y = qD, color = dw_type),
             inherit.aes = FALSE, size = 3, shape = 21, fill = "white", stroke = 1.2) +
  scale_color_manual(values = dw_colors) +
  scale_fill_manual(values = dw_colors) +
  scale_linetype_manual(values = c(interpolated="solid", extrapolated="dotted")) +
  facet_wrap(~ Order.q, labeller = as_labeller(order_labs), ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "iNEXT (incidence) by deadwood type",
       x = "Number of samples (observed & extrapolated)", y = "Estimated diversity (qD)",
       color = "Deadwood type", fill = "Deadwood type", linetype = "Type")
p_inext_dw
ggsave(paste0("plots/collectors_iNEXT_by_dw_type_", threshold, ".png"),
       p_inext_dw, width = 8, height = 8, dpi = 600)

# -----------------------------------------------------------------------------
# 2) Split FWD into aFWD/fFWD and repeat (rare + specaccum + iNEXT)
# -----------------------------------------------------------------------------
META1_afwd <- META1 %>%
  mutate(dw_type2 = dplyr::case_when(
    dw_type == "FWD" & (position %in% c("ATTACHED","AFWD","AFWD1","AFWD2","ATTACHTED")) ~ "aFWD",
    dw_type == "FWD" & (position %in% c("FALLEN","fFWD")) ~ "fFWD",
    TRUE ~ as.character(dw_type)
  ))

dw_keep2 <- c("SNAG","LOG","aFWD","fFWD","SOIL")

# Rarefaction
pdf(NULL)
rarecurve_list2 <- vegan::rarecurve(otu_work, step = 1500, label = FALSE, return = "list")
dev.off()
rare_df2 <- purrr::map2_dfr(
  rarecurve_list2, rownames(otu_work),
  ~tibble(sample = .y, reads = attr(.x,"Subsample"), richness = as.integer(.x))
) %>%
  left_join(META1_afwd %>% select(sample, dw_type2), by = "sample")

p_rare2 <- ggplot(rare_df2, aes(x = reads, y = richness, group = sample, color = dw_type2)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = dw_colors, na.value = "grey70") +
  theme_minimal(base_size = 13) +
  labs(title = "Rarefaction curves by deadwood type (FWD split)",
       x = "Reads (subsampled)", y = "Observed SH richness", color = "Deadwood type")
p_rare2
ggsave(paste0("plots/rarefaction_by_dw_type2_filtered_", threshold, ".png"),
       p_rare2, width = 8, height = 5, dpi = 300)

# Collector's (specaccum) for dw_type2
collector_curves2 <- lapply(dw_keep2, function(dw) {
  samples <- META1_afwd %>% filter(dw_type2 == dw) %>% pull(sample)
  if (length(samples) < 2) return(NULL)
  otus <- otu_work[rownames(otu_work) %in% samples, , drop = FALSE]
  otus <- otus[, colSums(otus) > 0, drop = FALSE]
  acc <- specaccum(otus, method = "random", permutations = 100)
  tibble(samples = acc$sites, richness = acc$richness, sd = acc$sd, dw_type2 = dw)
})
collector_df2 <- bind_rows(collector_curves2)

p_collectors2 <- ggplot(collector_df2, aes(x = samples, y = richness, color = dw_type2, fill = dw_type2)) +
  geom_line(linewidth = 1.1, alpha = 0.85) +
  geom_ribbon(aes(ymin = richness - sd, ymax = richness + sd), alpha = 0.15, color = NA) +
  scale_color_manual(values = dw_colors) +
  scale_fill_manual(values = dw_colors) +
  theme_minimal(base_size = 13) +
  labs(title = "Collector’s curves (specaccum) by deadwood type (FWD split)",
       x = "Number of samples", y = "Cumulative SH richness",
       color = "Deadwood type", fill = "Deadwood type")
p_collectors2
ggsave(paste0("plots/collectors_specaccum_by_dw_type2_", threshold, ".png"),
       p_collectors2, width = 9, height = 5, dpi = 300)

# iNEXT for dw_type2
incidence_list2 <- lapply(dw_keep2, function(dw) {
  samples <- META1_afwd %>% filter(dw_type2 == dw) %>% pull(sample)
  if (length(samples) < 2) return(NULL)
  otus <- otu_work[rownames(otu_work) %in% samples, , drop = FALSE]
  otus <- otus[, colSums(otus) > 0, drop = FALSE]
  freq <- colSums(otus > 0)
  c(length(samples), freq)
})
names(incidence_list2) <- dw_keep2
incidence_list2 <- incidence_list2[!sapply(incidence_list2, is.null)]

endpoints2 <- sapply(incidence_list2, function(x) min(ceiling(x[1] * 2), 170))
iNEXT_out2 <- lapply(names(incidence_list2), function(dw) {
  out <- iNEXT(incidence_list2[[dw]], q = c(0,1,2), datatype = "incidence_freq", endpoint = endpoints2[dw])
  as_tibble(out$iNextEst$size_based) %>% mutate(dw_type2 = dw)
})
iNEXT_df2 <- bind_rows(iNEXT_out2)
observed_n2 <- sapply(incidence_list2, function(x) x[1])
iNEXT_df2 <- iNEXT_df2 %>% mutate(interp_extrap = ifelse(t <= observed_n2[dw_type2], "interpolated", "extrapolated"))

interp_points2 <- iNEXT_df2 %>%
  group_by(dw_type2, Order.q) %>%
  filter(t == max(t[interp_extrap == "interpolated"])) %>%
  ungroup()

p_inext_dw2 <- ggplot(iNEXT_df2, aes(x = t, y = qD, color = dw_type2, fill = dw_type2, linetype = interp_extrap)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
  geom_point(data = interp_points2, aes(x = t, y = qD, color = dw_type2),
             inherit.aes = FALSE, size = 3, shape = 21, fill = "white", stroke = 1.2) +
  scale_color_manual(values = dw_colors) +
  scale_fill_manual(values = dw_colors) +
  scale_linetype_manual(values = c(interpolated="solid", extrapolated="dotted")) +
  facet_wrap(~ Order.q, labeller = as_labeller(order_labs), ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "iNEXT (incidence) by deadwood type (FWD split)",
       x = "Number of samples (observed & extrapolated)", y = "Estimated diversity (qD)",
       color = "Deadwood type", fill = "Deadwood type", linetype = "Type")
p_inext_dw2
ggsave(paste0("plots/collectors_iNEXT_by_dw_type2_", threshold, ".png"),
       p_inext_dw2, width = 6, height = 12, dpi = 600)

# -----------------------------------------------------------------------------
# 3) WOODY vs SOIL (iNEXT)
# -----------------------------------------------------------------------------
META1_ws <- META1_afwd %>% mutate(woody = ifelse(dw_type2 == "SOIL", "SOIL", "WOODY"))
incidence_combined <- list(
  WOODY = {
    s <- META1_ws %>% filter(woody == "WOODY") %>% pull(sample)
    ot <- otu_work[rownames(otu_work) %in% s, , drop = FALSE]; ot <- ot[, colSums(ot) > 0, drop = FALSE]
    c(length(s), colSums(ot > 0))
  },
  SOIL = {
    s <- META1_ws %>% filter(woody == "SOIL") %>% pull(sample)
    ot <- otu_work[rownames(otu_work) %in% s, , drop = FALSE]; ot <- ot[, colSums(ot) > 0, drop = FALSE]
    c(length(s), colSums(ot > 0))
  }
)
endpoints_combined <- sapply(incidence_combined, function(x) min(ceiling(x[1]*2), 300))
iNEXT_out_combined <- lapply(names(incidence_combined), function(g) {
  out <- iNEXT(incidence_combined[[g]], q = c(0,1,2), datatype = "incidence_freq", endpoint = endpoints_combined[g])
  as_tibble(out$iNextEst$size_based) %>% mutate(group = g)
})
iNEXT_df_combined <- bind_rows(iNEXT_out_combined)
obs_n_combined <- sapply(incidence_combined, function(x) x[1])
iNEXT_df_combined <- iNEXT_df_combined %>%
  mutate(interp_extrap = ifelse(t <= obs_n_combined[group], "interpolated", "extrapolated"))

interp_points_combined <- iNEXT_df_combined %>%
  group_by(group, Order.q) %>%
  filter(t == max(t[interp_extrap == "interpolated"])) %>%
  ungroup()

p_combined <- ggplot(iNEXT_df_combined, aes(x = t, y = qD, color = group, fill = group, linetype = interp_extrap)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
  geom_point(data = interp_points_combined, aes(x = t, y = qD, color = group),
             inherit.aes = FALSE, size = 3, shape = 21, fill = "white", stroke = 1.2) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_manual(values = c(interpolated="solid", extrapolated="dotted")) +
  facet_wrap(~ Order.q, labeller = as_labeller(order_labs), ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(title = "iNEXT: WOODY vs SOIL",
       x = "Number of samples (observed & extrapolated)", y = "Estimated diversity (qD)",
       color = "Category", fill = "Category", linetype = "Type")
p_combined
ggsave(paste0("plots/collectors_iNEXT_woody_vs_soil_", threshold, ".png"),
       p_combined, width = 10, height = 12, dpi = 300)

# -----------------------------------------------------------------------------
# 4) LOG-only contrasts: depth, position, ds_at_drill, decay_stage, size
#     (i) Rarefaction per sample per factor
#     (ii) iNEXT per factor level (with optional rarefaction to 250k reads)
# -----------------------------------------------------------------------------
# Common LOG meta w
log_meta_all <- META1 %>%
  filter(dw_type == "LOG")
str(log_meta_all)
str(otu_work)
# Depth
make_rarefaction_plot("depth_2", depth_colors, "Sample rarefaction (LOG by depth)", "depth",log_meta_all )

# Position (exclude non-longitudinal categories)
log_meta_pos <- log_meta_all %>%
  filter(is.na(position) | !position %in% c("AFWD","AFWD1","AFWD2","CROWN"))
make_rarefaction_plot("position", position_colors, "Sample rarefaction (LOG by position)", "position",log_meta_pos)

# ds_at_drill (0–5)
log_meta_ds <- META1 %>%
  filter(dw_type == "LOG", !is.na(ds_at_drill)) %>%
  mutate(ds_at_drill2 = as.character(ds_at_drill))
make_rarefaction_plot("ds_at_drill2", ds_colors, "Sample rarefaction (LOG by ds_at_drill)", "ds_at_drill",log_meta_ds)

# decay_stage (LIVING/EARLY/AVERAGE/LATE)
log_meta_dec <- META1 %>%
  filter(dw_type == "LOG", !is.na(decay_stage)) %>%
  mutate(decay_stage = as.character(decay_stage))
make_rarefaction_plot("decay_stage", DS_colors, "Sample rarefaction (LOG by decay_stage)", "decay_stage",log_meta_dec)

# ---- iNEXT builders for LOG contrasts (rarefy to 250k before incidence when specified)
# Depth
inc_depth <- build_inext_for_factor(
  meta_df = log_meta_all %>% filter(depth_2 %in% c("INNER","OUTER","MIXED"), depth_2 != "NO_DIST"),
  lvl_col = "depth_2", rarefy_to = 250000)
plot_inext_factor(inc_depth, depth_colors, "Collector’s curves (LOG by depth)", "depth", "Depth")

# Position
inc_pos <- build_inext_for_factor(
  meta_df = log_meta_all %>% filter(position_2 %in% c("BASE","MIDDLE","UPPER")),
  lvl_col = "position_2", rarefy_to = 250000)
plot_inext_factor(inc_pos, position_colors, "Collector’s curves (LOG by position)", "position", "Position")

# ds_at_drill (0–5; palette ds_colors)
inc_dsdrill <- build_inext_for_factor(
  meta_df = log_meta_ds %>% filter(sample %in% rownames(otu_work)),
  lvl_col = "ds_at_drill2", rarefy_to = 250000)
plot_inext_factor(inc_dsdrill, ds_colors, "Collector’s curves (LOG by ds_at_drill)", "ds_at_drill", "ds_at_drill")

# decay_stage (LIVING/EARLY/AVERAGE/LATE; palette DS_colors)
inc_dstg <- build_inext_for_factor(
  meta_df = log_meta_dec %>% filter(sample %in% rownames(otu_work)),
  lvl_col = "decay_stage", rarefy_to = 250000)
plot_inext_factor(inc_dstg, DS_colors, "Collector’s curves (LOG by decay stage)", "decay_stage", "Decay stage")

# size class (quartile-ish cuts; use sc_colors from Script 1)
meta_log_size <- META1 %>%
  filter(dw_type == "LOG", is.na(position) | position != "ENDO")
inc_sc <- build_inext_for_factor(
  meta_df = meta_log_size %>% filter(!is.na(size)),
  lvl_col = "size", rarefy_to = 250000)
plot_inext_factor(inc_sc, sc_colors, "Collector’s curves (LOG by size class)", "size", "Size class")
