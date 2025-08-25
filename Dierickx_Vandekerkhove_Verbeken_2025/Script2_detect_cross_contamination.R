# =============================================================================
# Script 2 — Detect cross-contamination / index-bleed + build exclusion lists
#   Dependencies: run Script 1 first (META, tax, color palettes, and optionally
#   exclude_taxa_list). This script:
#     (A) derives/locks positive-control signatures
#     (B) computes contamination metrics vs POSITIVE samples
#     (C) summarizes Malassezia SHs
#     (D) (one-off) combines signature lists into exclude_taxa_list.rds
#
# be run once to create 'exclude_taxa_list.rds' used by all analyses.
# =============================================================================
source("scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# -----------------------------------------------------------------------------
# 0) Minimal guards + inputs
# -----------------------------------------------------------------------------
# META/META1 and tax should come from Script 1; load if missing (kept linear & minimal)
unique(META$dw_type)
head(tax)
# Choose the threshold you want to examine for cross-contam metrics
threshold <- "1_5"
otu_table <- read_tsv(paste0("otu_table_threshold_", threshold, ".tsv"))

# Known-problematic samples (same as Script 1). 
# exclude_samples <- c("P3_50","P1_24","P2_93","P2_86","P3_63","P1_66","P3_47")
# META1 <- META1 %>% filter(!sample %in% exclude_samples)

# OTU matrix for the chosen threshold (post sample-exclusion)
otu_matrix <- otu_table %>%
  filter(sample %in% META1$sample) %>%
  column_to_rownames("sample") %>%
  as.data.frame()

otu_matrix <- otu_matrix[, colSums(otu_matrix) > 0, drop = FALSE]

cat("Samples remaining after manual exclusion:", nrow(otu_matrix), "\n")
cat("SHs/OTUs remaining after dropping zeros:", ncol(otu_matrix), "\n")

# -----------------------------------------------------------------------------
# 1) Derive POSITIVE-control signature SHs (from POSITIVE samples)
# -----------------------------------------------------------------------------
pos_samples <- META1 %>%
  filter(dw_type == "POSITIVE") %>%
  pull(sample)

otu_pos <- otu_matrix[rownames(otu_matrix) %in% pos_samples, , drop = FALSE]

# Long format + taxonomy for inspection
otu_pos_long <- otu_pos %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "n_reads") %>%
  filter(n_reads > 0)

otu_pos_tax <- otu_pos_long %>%
  left_join(tax, by = "sh_code") %>%
  arrange(desc(n_reads)) %>%
  select(sample, sh_code, n_reads, lowest_taxon) %>%
  filter(n_reads > 20)  # threshold for visibility only

print(otu_pos_tax, n = 20)
sh_codes <- unique(otu_pos_tax$sh_code)
n_pos_filtered <- dplyr::n_distinct(otu_pos_tax$sample)

pos_sig_summary <- otu_pos_tax %>%
  distinct(sample, sh_code, .keep_all = TRUE) %>%
  group_by(sh_code) %>%
  summarise(
    n_pos_samples  = n_distinct(sample),
    prevalence_pos = n_pos_samples / n_pos_filtered,
    total_reads    = sum(n_reads),
    mean_reads     = mean(n_reads),
    median_reads   = median(n_reads),
    max_reads      = max(n_reads),
    lowest_taxon   = dplyr::first(lowest_taxon),
    .groups = "drop"
  ) %>%
  arrange(desc(n_pos_samples), desc(total_reads))

print(pos_sig_summary, n = Inf)

# If you want to lock to a specific 5-marker set for a given threshold:

pos5_0_5 <- c(  "SH1478158.10FU",   "SH1780941.10FU",  "SH1478154.10FU",  "SH1959848.10FU",   "SH1780943.10FU")
pos5_1_0 <- c(  "SH1075157.10FU",   "SH1291792.10FU", "SH1415760.10FU", "SH1212188.10FU",  "SH1212189.10FU")
pos5_1_5 <- c(  "SH0751670.10FU",  "SH0927443.10FU",   "SH0861376.10FU",    "SH0861374.10FU",    "SH1027369.10FU")
pos5_2_0 <- c(  "SH0476433.10FU",    "SH0625324.10FU",  "SH0567732.10FU",  "SH0567733.10FU",  "SH0709822.10FU")
pos5_2_5 <- c(  "SH0234034.10FU",    "SH0364181.10FU",  "SH0312698.10FU",  "SH0312695.10FU",  "SH0438243.10FU")
pos5_3_0 <- c(  "SH0015350.10FU",   "SH0131891.10FU", "SH0084987.10FU", "SH0084989.10FU", "SH0198348.10FU" )
mock_taxa_list <- list(  "0_5" = pos5_0_5,  "1_0" = pos5_1_0,  "1_5" = pos5_1_5,  "2_0" = pos5_2_0,  "2_5" = pos5_2_5,  "3_0" = pos5_3_0)

# -----------------------------------------------------------------------------
# 2) Cross-contamination metrics vs POSITIVE controls
# -----------------------------------------------------------------------------
detect_cross_contam <- function(otu_df, meta, positive_label = "POSITIVE",
                                signature_sh = NULL,
                                frac_reads_cut = 0.20,
                                jaccard_cut = 0.20,
                                prop_markers_cut = 0.40) {
  stopifnot(!is.null(rownames(otu_df)))
  stopifnot(all(meta$sample %in% rownames(otu_df)))
  otu <- as.matrix(otu_df)
  storage.mode(otu) <- "double"
  
  # POSITIVE samples
  pos_samples <- meta %>% filter(dw_type == positive_label) %>% pull(sample)
  pos_samples <- intersect(pos_samples, rownames(otu))
  if (length(pos_samples) == 0) stop("No POSITIVE samples found in META$dw_type.")
  
  pos_otu <- otu[pos_samples, , drop = FALSE]
  
  # Signature definition
  if (is.null(signature_sh)) {
    pos_signature <- colnames(pos_otu)[colSums(pos_otu > 0) > 0]  # union of positives
  } else {
    pos_signature <- intersect(signature_sh, colnames(otu))
    if (length(setdiff(signature_sh, pos_signature)) > 0) {
      warning("Some signature SHs not in OTU table: ",
              paste(setdiff(signature_sh, pos_signature), collapse = ", "))
    }
  }
  p_sh <- length(pos_signature)
  if (p_sh == 0) stop("Positive signature is empty after intersection with OTU table.")
  
  pos_pool <- colSums(pos_otu[, pos_signature, drop = FALSE])
  
  presence <- (otu > 0) * 1
  s_sh_vec <- rowSums(presence)
  a_sig    <- rowSums(presence[, pos_signature, drop = FALSE])     # markers detected
  shared_reads_vec <- rowSums(otu[, pos_signature, drop = FALSE])
  total_reads_vec  <- rowSums(otu)
  
  b <- s_sh_vec - a_sig
  c <- p_sh - a_sig
  jaccard_to_pos_union <- ifelse((a_sig + b + c) > 0, a_sig / (a_sig + b + c), NA_real_)
  prop_markers_detected <- ifelse(p_sh > 0, a_sig / p_sh, NA_real_)
  
  bray_to_pos_pool <- vapply(seq_len(nrow(otu)), function(i) {
    x <- otu[i, pos_signature, drop = TRUE]
    denom <- sum(x + pos_pool)
    if (denom == 0) return(NA_real_)
    sum(abs(x - pos_pool)) / denom
  }, numeric(1))
  names(bray_to_pos_pool) <- rownames(otu)
  
  pos_mat <- otu[, pos_signature, drop = FALSE]
  if (ncol(pos_mat) == 0) {
    top_pos_sh    <- rep(NA_character_, nrow(otu))
    top_pos_reads <- rep(NA_real_,      nrow(otu))
  } else {
    j_idx <- max.col(pos_mat, ties.method = "first")
    top_pos_reads <- pos_mat[cbind(seq_len(nrow(pos_mat)), j_idx)]
    top_pos_sh    <- colnames(pos_mat)[j_idx]
    zero_row <- rowSums(pos_mat) == 0
    top_pos_sh[zero_row]    <- NA_character_
    top_pos_reads[zero_row] <- 0
  }
  
  contam_metrics <- tibble(
    sample = rownames(otu),
    dw_type = meta$dw_type[match(rownames(otu), meta$sample)],
    n_markers_detected   = as.integer(a_sig),
    shared_reads         = as.numeric(shared_reads_vec),
    total_reads          = as.numeric(total_reads_vec),
    frac_reads_shared    = ifelse(total_reads_vec > 0, shared_reads_vec / total_reads_vec, NA_real_),
    prop_markers_detected = prop_markers_detected,
    jaccard_to_pos_union  = jaccard_to_pos_union,
    bray_to_pos_pool      = bray_to_pos_pool,
    top_pos_sh            = top_pos_sh,
    top_pos_reads         = top_pos_reads
  ) %>%
    arrange(desc(frac_reads_shared), desc(n_markers_detected))
  
  suspects <- contam_metrics %>%
    filter(dw_type != positive_label) %>%
    filter(frac_reads_shared >= frac_reads_cut |
             jaccard_to_pos_union >= jaccard_cut |
             prop_markers_detected >= prop_markers_cut) %>%
    arrange(desc(frac_reads_shared), desc(n_markers_detected))
  
  list(
    signature_SHs  = pos_signature,
    contam_metrics = contam_metrics,
    suspect        = suspects
  )
}

# Example run: lock to your 5 markers for 3_0 and use lenient cutoffs
res_cc <- detect_cross_contam(
  otu_matrix, META1,
  positive_label   = "POSITIVE",
  signature_sh     = mock_taxa_list[["1_5"]], # adjust for threshold here
  frac_reads_cut   = 0.05,   # ≥5% of reads from markers
  jaccard_cut      = 0.05,   # ≥5% overlap in presence
  prop_markers_cut = 0.40    # ≥2 of 5 markers detected
)

glimpse(res_cc$contam_metrics)
head(res_cc$suspect, 30)
length(res_cc$signature_SHs)  # should be 5

# Enrich suspects with taxonomy and a human-readable ID if available
suspect_with_tax <- res_cc$suspect %>%
  left_join(tax, by = c("top_pos_sh" = "sh_code")) %>%
  left_join(META1 %>% select(sample, sample_id), by = "sample") %>%
  select(sample, sample_id, dw_type, n_markers_detected, frac_reads_shared,
         prop_markers_detected, jaccard_to_pos_union, bray_to_pos_pool,
         top_pos_sh, top_pos_reads, lowest_taxon)

print(suspect_with_tax, n = Inf)
# In practice, most flags with tiny shared fractions (e.g., 1e-5) likely reflect
# index bleed/barcode hopping rather than material contamination - all are filtered out with minimum per sample reads of an SH of 10.

# -----------------------------------------------------------------------------
# 3) Summarize Malassezia across the chosen table (fuzzy genus/species match)
# -----------------------------------------------------------------------------
malassezia_sh_all <- tax %>%
  mutate(genus = trimws(genus), lowest_taxon = trimws(lowest_taxon)) %>%
  filter(
    (!is.na(genus)        & str_detect(genus,        regex("^malassezia", ignore_case = TRUE))) |
      (!is.na(lowest_taxon) & str_detect(lowest_taxon, regex("^malassezia", ignore_case = TRUE)))
  ) %>%
  distinct(sh_code) %>%
  pull(sh_code)

malassezia_sh_present <- intersect(malassezia_sh_all, colnames(otu_matrix))

malassezia_summary <- {
  shs <- malassezia_sh_present
  if (length(shs) == 0) {
    tibble()
  } else {
    mat <- otu_matrix[, shs, drop = FALSE]
    tibble(
      sh_code     = shs,
      prevalence  = colSums(mat > 0),
      total_reads = colSums(mat)
    ) %>%
      left_join(tax %>% select(sh_code, species, lowest_taxon, genus), by = "sh_code") %>%
      arrange(desc(total_reads), desc(prevalence))
  }
}
print(malassezia_summary, n = 50)

# Threshold-specific Malassezia reference lists 
malassezia_3_0 <- c("SH0179873.10FU","SH0145538.10FU","SH0071936.10FU","SH0145167.10FU","SH0163328.10FU")
malassezia_2_5 <- c("SH0298309.10FU","SH0378777.10FU","SH0379175.10FU","SH0399021.10FU","SH0417661.10FU")
malassezia_2_0 <- c("SH0145538.10FU","SH0551549.10FU","SH0641696.10FU","SH0642121.10FU","SH0664850.10FU","SH0686389.10FU","SH0686393.10FU")
malassezia_1_5 <- c("SH0379175.10FU","SH0842694.10FU","SH0946624.10FU","SH0947116.10FU","SH0973982.10FU","SH0999786.10FU","SH0999796.10FU")
malassezia_1_0 <- c("SH0642121.10FU","SH1189575.10FU","SH1315515.10FU","SH1316101.10FU","SH1349449.10FU","SH1381990.10FU","SH1382014.10FU")
malassezia_0_5 <- c("SH1640451.10FU","SH1640457.10FU","SH1640471.10FU","SH1815619.10FU","SH1816407.10FU",
                    "SH1864147.10FU","SH1912099.10FU","SH1912103.10FU","SH1912121.10FU","SH1912125.10FU",
                    "SH1912138.10FU","SH1912152.10FU")


malassezia_list <- list(
  "0_5" = malassezia_0_5, "1_0" = malassezia_1_0, "1_5" = malassezia_1_5,
  "2_0" = malassezia_2_0, "2_5" = malassezia_2_5, "3_0" = malassezia_3_0
)

# -----------------------------------------------------------------------------
# 4) “Benin taxa” lists per threshold  (contaminated sample, look for index bleeds)
#     NOTE: The P3_50 inspection must use an *unfiltered* table. 
# -----------------------------------------------------------------------------

otu_table_raw <- read_tsv(paste0("otu_table_threshold_", threshold, ".tsv")) %>%
  column_to_rownames("sample")

if ("P3_50" %in% rownames(otu_table_raw)) {
  abund_vec <- otu_table_raw["P3_50", , drop = TRUE]
  p3_50_abund <- tibble(
    sh_code = names(abund_vec),
    n_reads = as.numeric(abund_vec)
  ) %>%
    filter(n_reads > 0) %>%
    left_join(tax, by = "sh_code") %>%
    mutate(rel_reads = 100 * (n_reads / sum(n_reads))) %>%
    arrange(desc(n_reads)) %>%
    filter(class != "Malasseziomycetes", class != "Atractiellomycetes")
  print(p3_50_abund, n = 12)
}

# Reference lists with hopping taxa from Benin (manually constructed)
benin_taxa_0_5 <- c("SH1781110.10FU", "SH1612970.10FU", "SH1724961.10FU", "SH1556688.10FU", "SH1759069.10FU", "SH1492877.10FU", "SH1556710.10FU")
benin_taxa_1_0 <- c("SH1291895.10FU", "SH1169545.10FU", "SH1251801.10FU", "SH1127738.10FU", "SH1275826.10FU", "SH1084675.10FU","SH1127780.10FU")
benin_taxa_1_5 <- c("SH0927530.10FU", "SH0826188.10FU", "SH0894392.10FU", "SH0792195.10FU", "SH0913955.10FU", "SH0759259.10FU", "SH0883406.10FU", "SH0840439.10FU")
benin_taxa_2_0 <- c("SH0625401.10FU", "SH0537565.10FU", "SH0596667.10FU", "SH0508979.10FU", "SH0613590.10FU", "SH0482682.10FU", "SH0587089.10FU", "SH0549557.10FU")
benin_taxa_2_5 <- c("SH0364257.10FU", "SH0286004.10FU", "SH0338627.10FU", "SH0261272.10FU", "SH0353660.10FU", "SH0239487.10FU", "SH0330104.10FU", "SH0296490.10FU")
benin_taxa_3_0 <- c("SH0131965.10FU", "SH0060926.10FU", "SH0108745.10FU", "SH0038955.10FU", "SH0122430.10FU", "SH0020238.10FU", "SH0100935.10FU", "SH0070250.10FU")

benin_taxa_list <- list(
  "0_5" = benin_taxa_0_5, "1_0" = benin_taxa_1_0, "1_5" = benin_taxa_1_5,
  "2_0" = benin_taxa_2_0, "2_5" = benin_taxa_2_5, "3_0" = benin_taxa_3_0
)

# -----------------------------------------------------------------------------
# 5) Build & save the combined exclude_taxa_list  to use in import script
# -----------------------------------------------------------------------------
thresholds <- c("3_0","2_5","2_0","1_5","1_0","0_5")
exclude_taxa_list <- lapply(thresholds, function(th) {
  unique(c(benin_taxa_list[[th]], malassezia_list[[th]], mock_taxa_list[[th]]))
})
names(exclude_taxa_list) <- thresholds

# Save so Script 1 can apply at import time
saveRDS(exclude_taxa_list, file = "exclude_taxa_list.rds")

# -----------------------------------------------------------------------------
# 6) Example application to another threshold (pre-removal stats)
# -----------------------------------------------------------------------------
threshold <- "0_5"
otu_table <- read_tsv(paste0("otu_table_threshold_", threshold, ".tsv"))

otu_matrix <- otu_table %>%
  filter(sample %in% META1$sample) %>%
  column_to_rownames("sample") %>%
  as.data.frame()

present_exclude <- intersect(exclude_taxa_list[[threshold]], colnames(otu_matrix))

present_exclude_tax <- tibble(sh_code = present_exclude) %>%
  left_join(tax, by = "sh_code")

if (length(present_exclude) > 0) {
  totals <- colSums(otu_matrix[, present_exclude, drop = FALSE])
  prev   <- colSums(otu_matrix[, present_exclude, drop = FALSE] > 0)
  present_exclude_tax <- present_exclude_tax %>%
    mutate(
      total_reads = totals[sh_code],
      prevalence  = prev[sh_code]
    ) %>%
    arrange(desc(total_reads), desc(prevalence))
}

# Remove excluded taxa + drop zero columns
otu_matrix <- otu_matrix[, setdiff(colnames(otu_matrix), present_exclude), drop = FALSE]
otu_matrix <- otu_matrix[, colSums(otu_matrix) > 0, drop = FALSE]

print(present_exclude_tax, n = Inf)

# -----------------------------------------------------------------------------
# 7) EXTRA # Check which samples are dominated by Fagus genus
# -----------------------------------------------------------------------------
#### EXTRA ####
# Make sure sample names are in a column
otu_long <- otu_matrix %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "reads")

# Join taxonomy info by sh_code
otu_long <- otu_long %>%
  left_join(
    tax %>% select(sh_code, genus),
    by = "sh_code"
  )

# Aggregate reads by genus per sample
otu_genus <- otu_long %>%
  group_by(sample, genus) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE), .groups = "drop")

# Convert to relative abundance per sample
otu_genus <- otu_genus %>%
  group_by(sample) %>%
  mutate(rel_abund = total_reads / sum(total_reads, na.rm = TRUE)) %>%
  ungroup()

# Find dominant genus per sample
dominant <- otu_genus %>%
  group_by(sample) %>%
  slice_max(rel_abund, n = 1, with_ties = FALSE) %>%
  ungroup()
dominant %>%
  count(genus, sort = TRUE) %>%
  mutate(
    percentage = round(n / sum(n) * 100, 1),
    cumulative_pct = round(cumsum(percentage), 1)) %>% 
  print(n=15)
# Filter where Fagus is dominant
fagus_dominant_samples <- dominant %>%
  filter(genus == "Fagus")
fagus_dominant_samples %>% 
  arrange(desc(rel_abund)) %>% 
  left_join(
    META1 %>% select(ds_at_drill, decay_stage, sample),
    by = "sample"
  )
# add fagus-dominated samples to exclude list c("MI_33", "MI_41", "MI_42", "P1_19", "MI_45" )
