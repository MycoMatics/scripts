# =============================================================================
# Script: Script2_venns_and_species_lists.R
# Project: Deadwood fungal community analyses
# Purpose: Build species lists, overlap summaries, and Venn-style plots for
#          deadwood fractions and woody-versus-soil comparisons
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
# Outputs:
# - tables/ALL/*.tsv
# - plots/ALL/*.png
# - plots/OVERLAP/*.png
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(eulerr)
  library(ggplot2)
  library(VennDiagram)
  library(ggVennDiagram)
  library(sf)
})

assert_objects(c("META1", "otu_matrix_filt", "tax", "threshold"))

# ---- Directories -------------------------------------------------------------
ensure_dir("tables/ALL")
ensure_dir("tables/OVERLAP")
ensure_dir("plots/ALL")
ensure_dir("plots/OVERLAP")

# ---- Reproducibility parameters ----------------------------------------------
SEED_VENN_DW   <- 123
N_ITER_DW      <- 10000L
N_DRAW_DW      <- 20L
SEED_WOODYSOIL <- 42
N_ITER_WS      <- 500L
N_DRAW_WS      <- 6L

# ---- Shared objects ----------------------------------------------------------
common_samples <- intersect(META1$sample, rownames(otu_matrix_filt))

meta <- META1 %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, rownames(otu_matrix_filt)))

otu <- otu_matrix_filt[meta$sample, , drop = FALSE]
stopifnot(identical(meta$sample, rownames(otu)))

meta_dw <- meta %>%
  mutate(dw_simple = factor(dw_type2)) %>%
  filter(!is.na(dw_simple)) %>%
  droplevels()

otu_dw <- otu[meta_dw$sample, , drop = FALSE]
stopifnot(identical(meta_dw$sample, rownames(otu_dw)))

libsize <- rowSums(otu_dw)
libsize_safe <- pmax(libsize, 1)
otu_rel <- sweep(otu_dw, 1, libsize_safe, "/") * 100

# ---- Top 10 taxa per deadwood type -------------------------------------------
types <- levels(meta_dw$dw_simple)

for (tt in types) {
  idx <- which(meta_dw$dw_simple == tt)
  if (length(idx) == 0L) next
  
  otu_rel_sub <- otu_rel[idx, , drop = FALSE]
  n_samp <- length(idx)
  
  tot_pct  <- colSums(otu_rel_sub, na.rm = TRUE)
  mean_pct <- tot_pct / n_samp
  
  ord <- order(mean_pct, decreasing = TRUE)
  ord <- ord[mean_pct[ord] > 0]
  if (length(ord) == 0L) next
  
  top_idx <- head(ord, 10L)
  top_sh  <- colnames(otu_rel)[top_idx]
  
  tax_info <- tibble(sh_code = top_sh) %>%
    left_join(select(tax, sh_code, species, genus), by = "sh_code") %>%
    mutate(taxon = coalesce(species, genus, sh_code))
  
  top_tbl <- tibble(
    rank          = seq_along(top_idx),
    sh_code       = tax_info$sh_code,
    taxon         = tax_info$taxon,
    n_samples     = n_samp,
    mean_percent  = round(mean_pct[top_idx], 3),
    total_percent = round(tot_pct[top_idx], 3)
  )
  
  cat("\n=========================================\n")
  cat("Top taxa for deadwood type:", tt, "\n")
  cat("=========================================\n")
  print(top_tbl, n = nrow(top_tbl))
}

# ---- Long table for all deadwood types ---------------------------------------
meta_all <- meta_dw
otu_all <- otu_matrix_filt[meta_all$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_all), meta_all$sample))

otu_long_all <- otu_all %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "abundance") %>%
  filter(abundance > 0) %>%
  left_join(select(meta_all, sample, dw_type2), by = "sample") %>%
  filter(!is.na(dw_type2)) %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = abundance / sum(abundance)) %>%
  ungroup() %>%
  left_join(select(tax, sh_code, species, genus), by = "sh_code")

all_species_by_dw <- summarise_alltaxa_by(otu_long_all, dw_type2)
combined_path_dw <- file.path("tables/ALL", paste0("ALL_species_by_dwtype_threshold", threshold, ".tsv"))
readr::write_tsv(all_species_by_dw, combined_path_dw)
message("Wrote: ", combined_path_dw)

message("Top 10 by dw_type2:")
print(
  all_species_by_dw %>%
    group_by(dw_type2) %>%
    slice_head(n = 10),
  n = 50
)

# ---- Bootstrap Venn for deadwood fractions -----------------------------------
dw_keep <- c("aFWD", "fFWD", "LOG", "SNAG")
dw_levels <- c("SNAG", "LOG", "fFWD", "aFWD")
region_levels <- make_region_levels(dw_levels)
print(region_levels)
length(region_levels)

otu_long2 <- otu_matrix_filt %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "n_seq") %>%
  filter(n_seq > 0) %>%
  left_join(select(META1, sample, dw_type2, position_2), by = "sample") %>%
  left_join(select(tax, sh_code, species, genus), by = "sh_code") %>%
  filter(dw_type2 %in% dw_keep) %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = n_seq / sum(n_seq)) %>%
  ungroup() %>%
  filter(!is.na(dw_type2))

species_by_sample <- otu_long2 %>%
  group_by(dw_type2, sample) %>%
  summarise(sh_codes = list(unique(sh_code)), .groups = "drop")

species_list_by_dwtype <- species_by_sample %>%
  group_by(dw_type2) %>%
  summarise(
    sh_codes = list(sort(unique(unlist(sh_codes)))),
    .groups  = "drop"
  ) %>%
  deframe()

cat("Observed unique SHs per deadwood type (all samples)\n")
print(sapply(species_list_by_dwtype[dw_levels], length))

set.seed(SEED_VENN_DW)
n_iter <- N_ITER_DW
n_draw_target <- N_DRAW_DW

boot_regions <- vector("list", n_iter)
boot_rich <- vector("list", n_iter)

for (b in seq_len(n_iter)) {
  sets_b <- setNames(vector("list", length(dw_levels)), dw_levels)
  rich_rows <- vector("list", length(dw_levels))
  
  for (i in seq_along(dw_levels)) {
    dw <- dw_levels[i]
    df_dw <- dplyr::filter(species_by_sample, dw_type2 == dw)
    
    if (nrow(df_dw) == 0L) {
      sets_b[[dw]] <- character(0)
      rich_rows[[i]] <- tibble(iter = b, dw_type2 = dw, richness = NA_integer_)
    } else {
      n_draw <- min(n_draw_target, nrow(df_dw))
      idx <- sample.int(nrow(df_dw), n_draw)
      sh_set <- unique(unlist(df_dw$sh_codes[idx]))
      sets_b[[dw]] <- sh_set
      rich_rows[[i]] <- tibble(iter = b, dw_type2 = dw, richness = length(sh_set))
    }
  }
  
  boot_rich[[b]] <- bind_rows(rich_rows)
  
  reg_counts <- region_counts_from_sets(sets_b, region_levels)
  boot_regions[[b]] <- tibble(
    iter   = b,
    region = region_levels,
    count  = as.integer(reg_counts)
  )
}

venn_boot_richness <- bind_rows(boot_rich)
venn_boot_regions  <- bind_rows(boot_regions)

venn_rich_summary <- venn_boot_richness %>%
  group_by(dw_type2) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness   = sd(richness, na.rm = TRUE),
    .groups       = "drop"
  ) %>%
  arrange(factor(dw_type2, levels = dw_levels))

cat("\nRarefied richness at ", n_draw_target, " samples per substrate across ", n_iter, " iterations\n", sep = "")
print(venn_rich_summary)

venn_region_summary <- venn_boot_regions %>%
  group_by(region) %>%
  summarise(
    mean_count = mean(count),
    sd_count   = sd(count),
    .groups    = "drop"
  ) %>%
  arrange(match(region, region_levels))

cat("\nBootstrap mean and sd counts per Venn region (", n_draw_target, " samples per substrate)\n", sep = "")
print(venn_region_summary)

dw_order_plot <- c("LOG", "aFWD", "fFWD", "SNAG")
sets_for_geometry <- species_list_by_dwtype[dw_order_plot]
sets_for_geometry <- sets_for_geometry[!vapply(sets_for_geometry, is.null, logical(1))]

venn_obj <- ggVennDiagram::Venn(sets_for_geometry)
vd <- ggVennDiagram::process_data(venn_obj)

region_label_df <- vd$regionLabel
set_data <- vd$setData
set_items <- set_data$item
names(set_items) <- set_data$name

region_combo <- character(nrow(region_label_df))
for (i in seq_len(nrow(region_label_df))) {
  items_i <- region_label_df$item[[i]]
  if (length(items_i) == 0L) {
    region_combo[i] <- NA_character_
  } else {
    sp <- items_i[1]
    in_set <- vapply(set_items, function(x) sp %in% x, logical(1))
    region_combo[i] <- make_combo(names(set_items)[in_set])
  }
}
region_label_df$combo <- region_combo

venn_region_summary2 <- venn_region_summary %>%
  mutate(region_id = sapply(strsplit(region, "&"), make_combo))

region_label_annot <- region_label_df %>%
  left_join(venn_region_summary2, by = c("combo" = "region_id")) %>%
  mutate(label = ifelse(is.na(mean_count), "", sprintf("%0.1f\n± %0.1f", mean_count, sd_count)))

cat("Region combos and labels used for plotting:\n")
print(region_label_annot[, c("combo", "label")])

p_venn_boot <- ggVennDiagram::ggVennDiagram(
  sets_for_geometry,
  label = "none"
) +
  scale_fill_gradient(
    low   = "white",
    high  = "white",
    guide = "none"
  ) +
  scale_color_manual(
    values = setNames(rep("black", length(sets_for_geometry)), names(sets_for_geometry)),
    guide = "none"
  ) +
  geom_text(
    data = subset(region_label_annot, label != ""),
    aes(x = X, y = Y, label = label),
    size = 3
  ) +
  ggtitle(
    paste0(
      "Bootstrap Venn, ", n_draw_target, " samples per substrate\n",
      "region labels are mean ± sd over ", n_iter, " iterations"
    )
  )

plotname_boot <- file.path("plots/ALL", paste0("venn_dw_bootstrap_mean_sd_n20_", threshold, ".png"))
ggsave(plotname_boot, p_venn_boot, width = 6, height = 6, dpi = 900)
message("Wrote: ", plotname_boot)

# ---- Total species overlap across deadwood fractions -------------------------
species_list_by_dwtype <- otu_long2 %>%
  filter(!is.na(dw_type2)) %>%
  group_by(dw_type2) %>%
  summarise(sh_codes = list(unique(sh_code)), .groups = "drop") %>%
  deframe()

venn_dw_types <- intersect(names(dw_colors), names(species_list_by_dwtype))
venn_colors <- dw_colors[venn_dw_types]
species_list_for_venn <- species_list_by_dwtype[venn_dw_types]

cat("Number of unique SHs per deadwood type:\n")
print(sapply(species_list_for_venn, length))

venn_plot <- venn.diagram(
  x = species_list_for_venn,
  category.names = venn_dw_types,
  filename = NULL,
  fill = venn_colors,
  col = venn_colors,
  alpha = 0.5,
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontfamily = "sans"
)

plotname <- file.path("plots/ALL", paste0("venn_dw_type_all_species_", threshold, ".png"))
ggsave(plotname, venn_plot, width = 4, height = 4, dpi = 600)
message("Wrote: ", plotname)

# ---- WOODY vs SOIL species sets ----------------------------------------------
dw_keep <- c("aFWD", "fFWD", "LOG", "SNAG", "SOIL")

otu_long_ws <- otu_matrix_filt %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "n_seq") %>%
  filter(n_seq > 0) %>%
  left_join(select(META1, sample, dw_type2, position_2), by = "sample") %>%
  filter(dw_type2 %in% dw_keep) %>%
  mutate(
    substrate2 = case_when(
      dw_type2 == "SOIL" ~ "SOIL",
      dw_type2 %in% c("aFWD", "fFWD", "LOG", "SNAG") ~ "WOODY",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(substrate2))

species_list_ws <- otu_long_ws %>%
  group_by(substrate2) %>%
  summarise(sh_codes = list(unique(sh_code)), .groups = "drop") %>%
  deframe()

cat("Number of unique SHs per group:\n")
print(sapply(species_list_ws, length))

woody_set <- species_list_ws[["WOODY"]]
soil_set  <- species_list_ws[["SOIL"]]
n_inter <- length(intersect(woody_set, soil_set))
n_union <- length(union(woody_set, soil_set))
jac_sim <- if (n_union > 0) n_inter / n_union else NA_real_

cat("\nShared SH between WOODY and SOIL:", n_inter, "\n")
cat("Jaccard similarity:", round(jac_sim, 3), "\n\n")

venn_categories <- names(species_list_ws)
venn_colors <- c(SOIL = "brown", WOODY = "#4E79A7")

venn_plot <- venn.diagram(
  x = species_list_ws,
  category.names = venn_categories,
  filename = NULL,
  fill = venn_colors,
  col = "black",
  alpha = 0.5,
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontfamily = "sans"
)

plotname <- file.path("plots/OVERLAP", paste0("venn_woody_vs_soil_threshold", threshold, ".png"))
ggplot2::ggsave(filename = plotname, plot = venn_plot, width = 4, height = 4, dpi = 600)
message("Wrote: ", plotname)

# ---- Rarefied WOODY vs SOIL comparison --------------------------------------
species_by_sample_ws <- otu_long_ws %>%
  group_by(substrate2, sample) %>%
  summarise(sh_codes = list(unique(sh_code)), .groups = "drop")

species_list_ws <- species_by_sample_ws %>%
  group_by(substrate2) %>%
  summarise(sh_codes = list(sort(unique(unlist(sh_codes)))), .groups = "drop") %>%
  deframe()

cat("Number of unique SHs per group (all samples):\n")
print(sapply(species_list_ws, length))

set.seed(SEED_WOODYSOIL)
n_iter <- N_ITER_WS
n_draw_target <- N_DRAW_WS

rows_woody <- dplyr::filter(species_by_sample_ws, substrate2 == "WOODY")
rows_soil  <- dplyr::filter(species_by_sample_ws, substrate2 == "SOIL")

cat("Samples per group for rarefaction:\n")
cat(" WOODY:", nrow(rows_woody), "samples\n")
cat(" SOIL :", nrow(rows_soil), "samples\n\n")

boot_list <- vector("list", n_iter)

for (b in seq_len(n_iter)) {
  if (nrow(rows_woody) == 0L || nrow(rows_soil) == 0L) next
  
  n_draw_w <- min(n_draw_target, nrow(rows_woody))
  n_draw_s <- min(n_draw_target, nrow(rows_soil))
  
  idx_w <- sample.int(nrow(rows_woody), n_draw_w)
  idx_s <- sample.int(nrow(rows_soil), n_draw_s)
  
  set_w <- sort(unique(unlist(rows_woody$sh_codes[idx_w])))
  set_s <- sort(unique(unlist(rows_soil$sh_codes[idx_s])))
  
  boot_list[[b]] <- tibble(
    iter = b,
    group = c("WOODY", "SOIL", "INTERSECTION", "UNION"),
    richness = c(
      length(set_w),
      length(set_s),
      length(intersect(set_w, set_s)),
      length(union(set_w, set_s))
    )
  )
}

boot_df <- bind_rows(boot_list)

rarefied_summary_ws <- boot_df %>%
  group_by(group) %>%
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    sd_richness   = sd(richness, na.rm = TRUE),
    .groups       = "drop"
  )

cat("Rarefied richness at ", n_draw_target, " samples per group across ", n_iter, " iterations:\n", sep = "")
print(rarefied_summary_ws)
