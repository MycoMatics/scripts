# =============================================================================
# Script 1 — Import & initial filtering (META, OTUs, colors)
#   - Loads metadata & OTU table for a chosen threshold
#   - Standardizes sample IDs and factor levels
#   - Excludes known-bad samples and threshold-specific contaminant SHs
#   - Applies per-sample abundance threshold and global prevalence/abundance filters
#   - (Optional) taxonomy-based OTU exclusion
#   - Produces: META1, tax, otu_matrix, otu_matrix_thr, otu_matrix_filt
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# --------------------------
# Inputs 
# --------------------------
#   5-marker set of positive  control taxa for a given threshold:
#pos5_0_5 <- c(  "SH1478158.10FU",   "SH1780941.10FU",  "SH1478154.10FU",  "SH1959848.10FU",   "SH1780943.10FU")
#pos5_1_0 <- c(  "SH1075157.10FU",   "SH1291792.10FU", "SH1415760.10FU", "SH1212188.10FU",  "SH1212189.10FU")
#pos5_1_5 <- c(  "SH0751670.10FU",  "SH0927443.10FU",   "SH0861376.10FU",    "SH0861374.10FU",    "SH1027369.10FU")
#pos5_2_0 <- c(  "SH0476433.10FU",    "SH0625324.10FU",  "SH0567732.10FU",  "SH0567733.10FU",  "SH0709822.10FU")
#pos5_2_5 <- c(  "SH0234034.10FU",    "SH0364181.10FU",  "SH0312698.10FU",  "SH0312695.10FU",  "SH0438243.10FU")
#pos5_3_0 <- c(  "SH0015350.10FU",   "SH0131891.10FU", "SH0084987.10FU", "SH0084989.10FU", "SH0198348.10FU" )
#mock_taxa_list <- list(  "0_5" = pos5_0_5,  "1_0" = pos5_1_0,  "1_5" = pos5_1_5,  "2_0" = pos5_2_0,  "2_5" = pos5_2_5,  "3_0" = pos5_3_0)
#    Threshold-specific Malassezia reference lists 
#malassezia_3_0 <- c("SH0179873.10FU","SH0145538.10FU","SH0071936.10FU","SH0145167.10FU","SH0163328.10FU")
#malassezia_2_5 <- c("SH0298309.10FU","SH0378777.10FU","SH0379175.10FU","SH0399021.10FU","SH0417661.10FU")
#malassezia_2_0 <- c("SH0145538.10FU","SH0551549.10FU","SH0641696.10FU","SH0642121.10FU","SH0664850.10FU","SH0686389.10FU","SH0686393.10FU")
#malassezia_1_5 <- c("SH0379175.10FU","SH0842694.10FU","SH0946624.10FU","SH0947116.10FU","SH0973982.10FU","SH0999786.10FU","SH0999796.10FU")
#malassezia_1_0 <- c("SH0642121.10FU","SH1189575.10FU","SH1315515.10FU","SH1316101.10FU","SH1349449.10FU","SH1381990.10FU","SH1382014.10FU")
#malassezia_0_5 <- c("SH1640451.10FU","SH1640457.10FU","SH1640471.10FU","SH1815619.10FU","SH1816407.10FU","SH1864147.10FU","SH1912099.10FU","SH1912103.10FU","SH1912121.10FU","SH1912125.10FU","SH1912138.10FU","SH1912152.10FU")
#malassezia_list <- list("0_5" = malassezia_0_5, "1_0" = malassezia_1_0, "1_5" = malassezia_1_5, "2_0" = malassezia_2_0, "2_5" = malassezia_2_5, "3_0" = malassezia_3_0)
#     Reference lists with hopping taxa from Benin (manually constructed)
#benin_taxa_0_5 <- c("SH1781110.10FU", "SH1612970.10FU", "SH1724961.10FU", "SH1556688.10FU", "SH1759069.10FU", "SH1492877.10FU", "SH1556710.10FU")
#benin_taxa_1_0 <- c("SH1291895.10FU", "SH1169545.10FU", "SH1251801.10FU", "SH1127738.10FU", "SH1275826.10FU", "SH1084675.10FU","SH1127780.10FU")
#benin_taxa_1_5 <- c("SH0927530.10FU", "SH0826188.10FU", "SH0894392.10FU", "SH0792195.10FU", "SH0913955.10FU", "SH0759259.10FU", "SH0883406.10FU", "SH0840439.10FU")
#benin_taxa_2_0 <- c("SH0625401.10FU", "SH0537565.10FU", "SH0596667.10FU", "SH0508979.10FU", "SH0613590.10FU", "SH0482682.10FU", "SH0587089.10FU", "SH0549557.10FU")
#benin_taxa_2_5 <- c("SH0364257.10FU", "SH0286004.10FU", "SH0338627.10FU", "SH0261272.10FU", "SH0353660.10FU", "SH0239487.10FU", "SH0330104.10FU", "SH0296490.10FU")
#benin_taxa_3_0 <- c("SH0131965.10FU", "SH0060926.10FU", "SH0108745.10FU", "SH0038955.10FU", "SH0122430.10FU", "SH0020238.10FU", "SH0100935.10FU", "SH0070250.10FU")
#benin_taxa_list <- list("0_5" = benin_taxa_0_5, "1_0" = benin_taxa_1_0, "1_5" = benin_taxa_1_5,  "2_0" = benin_taxa_2_0, "2_5" = benin_taxa_2_5, "3_0" = benin_taxa_3_0)
#thresholds <- c("3_0","2_5","2_0","1_5","1_0","0_5")
#exclude_taxa_list <- lapply(thresholds, function(th) {unique(c(benin_taxa_list[[th]], malassezia_list[[th]], mock_taxa_list[[th]]))})
#names(exclude_taxa_list) <- thresholds

exclude_taxa_list <- readRDS("exclude_taxa_list.rds")  # threshold -> character vector of SHs
str(exclude_taxa_list)

META <- read_delim(
  "META1.csv", delim = ";", escape_double = FALSE,
  col_types = cols(
    date = col_skip(),
    distance_base = col_integer(),
    height_from_log = col_double(),
    diameter_at_drill = col_double(),
    circum_at_drill   = col_double(),
    DS_at_drill       = col_double()
  ),
  trim_ws = TRUE
) %>%
  clean_names() %>%
  mutate(
    natman = clean_nat(natman),
    plate   = toupper(plate),
    barcode = sprintf("%02d", barcode),
    sample  = paste0(plate, "_", barcode),
    umi = as.factor(umi),
    # Standardize 'position' values + fix typos, then set consistent levels
    position = if_else(sample == "MI_33", "ENDO", position),
    position = if_else(sample == "P2_86", "SOUTH", position),
    position = if_else(sample == "P3_33", "SOUTH", position),
    position = toupper(position),
    position = dplyr::recode(position,
                             "MID" = "MIDDLE", "TOP" = "UPPER", "ATTACHTED" = "ATTACHED",
                             .default = position
    ),
    position = factor(
      position,
      levels = c("BASE","MIDDLE","UPPER", "NORTH","SOUTH","ATTACHED","ENDO","FALLEN")
    ),
    size = if_else(
      is.na(size),
      cut(
        diameter_at_drill,
        breaks = c(-Inf, 9, 19, 39, 69, 99, Inf),
        labels = c("VERY_FINE","FINE", "SMALL","MEDIUM","LARGE","VERY_LARGE"),
        right  = TRUE
      ),
      as.factor(size)  # Keep existing values as factor
    ),
    size = factor(size, levels = c("VERY_FINE","FINE", "SMALL","MEDIUM","LARGE","VERY_LARGE")),
    # Depth & decay stage factors
    depth_2 = factor(depth_2, levels = c("INNER", "MIXED","OUTER")),
    decay_stage = factor(decay_stage, levels = c("LIVING","EARLY","AVERAGE","LATE"), ordered = TRUE)) %>% 
    select(!c("helper","depth","main","sample_type"))

# --------------------------
# Color palettes (kept here so all downstream scripts can reuse them)
# --------------------------
dw_colors <- c( LOG = "#8B5C2D", FWD = "#E2B600", SNAG = "#59A14F", SOIL = "#4E79A7", aFWD = "#FFE066", fFWD = "#B89400", FRUITBODY = "#8E44AD", POSITIVE = "#C39BD3")
ds_colors <- c("0"="#66A61E","1"="#1B9E77","2"="#7570B3","3"="#E7298A","4"="#D95F02","5"="#A6761D")
DS_colors <- c(LIVING="#66A61E", EARLY="#1B9E77", AVERAGE="#E7298A", LATE="#A6761D")
group_colors <- c(WOODY="#A0522D", SOIL="#2E86AB")
depth_colors <- c(INNER="#0072B2", OUTER="#E69F00", MIXED="#CC79A7")
position_colors <- c(BASE="#009E73", MIDDLE="#D55E00", UPPER="#F0E442", ENDO="lightgreen")
sc_colors <- c("VERY_FINE"= "#5D4037","FINE"= "#BCAAA4", "SMALL"="#56B4E9", "MEDIUM"="#E15759", "LARGE"="#76B7B2", "VERY_LARGE"="#B07AA1")
microhab_base_cols <- c( BASE = position_colors[["BASE"]], MIDDLE = position_colors[["MIDDLE"]], UPPER = position_colors[["UPPER"]],ATTACHED = dw_colors[["aFWD"]], FALLEN = dw_colors[["fFWD"]], SNAG  = dw_colors[["SNAG"]])
aspect_colors <- c(NORTH = "navy", SOUTH = "orange")

# --------------------------
# Taxonomy table + threshold selection
# --------------------------
threshold <- "1_5"  # set the active per-threshold OTU table key
tax <- read_csv("/data/gent/458/vsc45818/DB/UNITE_SH_taxonomy_table.csv") %>%
  clean_names()  # expects columns like sh_code, kingdom, genus, ...

otu_table <- read_tsv(paste0("otu_table_threshold_", threshold, ".tsv"))

# --------------------------
# Build META1 (add initial richness per sample) and drop specific samples
# --------------------------
META1 <- otu_table %>%
  mutate(richness_initial = rowSums(select(., -sample) > 0)) %>%
  select(sample, richness_initial) %>%
  left_join(META, by = "sample")

# Known-problematic samples:
#  - 'P3_50': contamination (Benin species in a Belgian plate)
#  - Duplicates: keep the best (highest reads/richness)
#  - Low read count: P3_63, P1_66, P3_47
#  - almost all reads are Fagus sylvatica: MI_41  MI_33  MI_42  MI_45  P1_19  (all 99.5% or more Fagus sylvatica)
exclude_samples <- c("P3_50","P1_24","P2_93","P2_86","P3_63","P1_66","P3_47","MI_33", "MI_41", "MI_42", "P1_19", "MI_45")
META1 <- META1 %>% filter(!sample %in% exclude_samples) %>% filter(dw_type != "FRUITBODY",dw_type != "POSITIVE",dw_type != "NEGATIVE")
# --------------------------
# OTU matrix (rows = samples, cols = SHs) + remove threshold-specific contaminants
# --------------------------
otu_matrix <- otu_table %>%
  filter(sample %in% META1$sample) %>%
  column_to_rownames("sample")

#            # Total Sum Scaling (TSS) to median depth
#            # Compute library sizes per sample
#            lib_sizes <- rowSums(otu_matrix)
#            # Use the median library size as target depth
#            target_depth <- median(lib_sizes)
#            # Scale counts per sample to target depth
#            otu_matrix_tss <- sweep(otu_matrix, 1, lib_sizes, "/") * target_depth
#            # Optional: round back to integers if needed
#            otu_matrix_tss <- round(otu_matrix_tss)
#            # Check summary
#            summary(rowSums(otu_matrix_tss))
#            otu_matrix <- otu_matrix_tss
            
present_exclude <- intersect(exclude_taxa_list[[threshold]], colnames(otu_matrix))
otu_matrix <- otu_matrix[, setdiff(colnames(otu_matrix), present_exclude), drop = FALSE]
otu_matrix <- otu_matrix[, colSums(otu_matrix) > 0, drop = FALSE]

cat("Samples remaining after manual exclusion:", nrow(otu_matrix), "\n")
cat("SHs/OTUs remaining after dropping zeros:", ncol(otu_matrix), "\n")

# --------------------------
# Per-sample abundance thresholding (set counts < min_reads to 0)
# --------------------------
otu_matrix_thr <- threshold_by_sample(otu_matrix, min_reads = 10, drop_empty = FALSE)

# --------------------------
# Global OTU filters on the thresholded table
#   - Prevalence >= 1 samples
#   - Total abundance >= 10 reads
# --------------------------
prevalence_thr <- colSums(otu_matrix_thr > 0)
abundance_thr  <- colSums(otu_matrix_thr)
keep_otus <- (prevalence_thr >= 1) & (abundance_thr >= 10)

otu_matrix_filt <- otu_matrix_thr[, keep_otus, drop = FALSE]
otu_matrix_filt <- otu_matrix_filt[, colSums(otu_matrix_filt) > 0, drop = FALSE]

cat("After per-sample threshold:\n",
    "  Samples:", nrow(otu_matrix), "\n",
    "  OTUs (nonzero in ≥1 sample):", ncol(otu_matrix), "\n",
    "After global filtering on thresholded table:\n",
    "  OTUs kept:", ncol(otu_matrix_filt), "\n")

# --------------------------
# Optional taxonomy-based removal (disabled by default)
#   - Set 'bad_genera' to a non-empty vector to enable (example shown)
# --------------------------
bad_kingdom <- "Viridiplantae"
# Example to enable genus-based filtering:
# bad_genera  <- c("Kretzschmaria","Eutypa","Armillaria")
bad_genera  <- character(0)  # default: do not filter by genus

drop_sh <- tax %>%
  filter(
    sh_code %in% colnames(otu_matrix_filt),
    kingdom == bad_kingdom | genus %in% bad_genera
  ) %>%
  pull(sh_code)

if (length(drop_sh)) {
  cat("Dropping", length(drop_sh), "OTUs due to taxonomy filter\n")
  otu_matrix_filt <- otu_matrix_filt[, setdiff(colnames(otu_matrix_filt), drop_sh), drop = FALSE]
  otu_matrix_filt <- otu_matrix_filt[, colSums(otu_matrix_filt) > 0, drop = FALSE]
}
cat("Remaining OTUs after taxonomy filter:", ncol(otu_matrix_filt), "\n")

# --------------------------
# Update META1 with post-threshold richness/read depth
# --------------------------
richness_filt <- rowSums(otu_matrix_thr > 0)
reads_filt    <- rowSums(otu_matrix_thr)

META1 <- META1 %>%
  filter(sample %in% rownames(otu_matrix_thr)) %>%
  mutate(
    richness_filt = richness_filt[sample],
    reads_filt    = reads_filt[sample]
  )

# Quick sanity check:
head(META1$richness_filt)
str((otu_matrix_filt))
# append lowest_taxon to column names 
#lut <- setNames(tax$lowest_taxon, tax$sh_code)
#colnames(otu_matrix_filt) <- paste0(colnames(otu_matrix_filt), "|", ifelse(is.na(lut[colnames(otu_matrix_filt)]), "", lut[colnames(otu_matrix_filt)]))

# --------------------------
# (Optional) Save clean objects for downstream scripts
# --------------------------
saveRDS(list(META1=META1, tax=tax,
             otu_matrix=otu_matrix,
             otu_matrix_thr=otu_matrix_thr,
             otu_matrix_filt=otu_matrix_filt),
        file = paste0("import_objects_thr_", threshold, ".rds"))



