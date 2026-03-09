# =============================================================================
# Script: Script1_import_and_filter.R
# Project: Deadwood fungal community analyses
# Purpose: Import metadata and OTU tables, apply initial sample/abundance
#          filtering, and save core analysis objects for downstream scripts
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Inputs:
# - R_data/Final_tables.xlsx
# - R_data/taxonomy.csv
# - R_data/exclude_taxa_list.rds
#
# Outputs:
# - R_data/import_objects_thr_1_5.rds
# - optional session info text file
# =============================================================================

source("Scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})

set.seed(SEED_GLOBAL)

# ---- Paths and parameters ----------------------------------------------------
data_dir   <- "data"
threshold  <- "1_5"
thr_str    <- gsub("\\.", "_", as.character(threshold))
out_rds    <- file.path(data_dir, paste0("import_objects_thr_", threshold, ".rds"))
out_sess   <- file.path(data_dir, paste0("sessionInfo_import_thr_", threshold, ".txt"))

# ---- Check required input files ---------------------------------------------
required_files <- c(
  file.path(data_dir, "Final_tables.xlsx"),
  file.path(data_dir, "taxonomy.csv"),
  file.path(data_dir, "exclude_taxa_list.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "))
}

# ---- Inputs -----------------------------------------------------------------
exclude_taxa_list <- readRDS(file.path(data_dir, "exclude_taxa_list.rds"))
str(exclude_taxa_list)

otu_table <- readxl::read_excel(file.path(data_dir, "Final_tables.xlsx"), sheet = "SH_table_1_5") %>%
  rename(SHcode = 1) %>%
  mutate(SHcode = str_remove(SHcode, "\\|.*"))

META <- readxl::read_excel(file.path(data_dir, "Final_tables.xlsx"), sheet = "S1_Sample_Meta", skip = 1)
tax  <- readr::read_csv2(file.path(data_dir, "taxonomy.csv")) %>%
  clean_names()

str(otu_table)

# Known problematic samples: duplicates, low read count, Fagus sylvatica dominated (>99%)
exclude_samples <- c(
  "P3_50","P1_24","P2_93","P2_86","P3_63","P1_66","P3_47","MI_33",
  "MI_41","MI_42","P1_19","MI_45","P1_56","P2_84","P3_42","P3_41",
  "P3_37","P3_82"
)

META1 <- META %>%
  filter(!sample %in% exclude_samples) %>% 
  mutate(diameter_at_drill_z = as.numeric(diameter_at_drill_z))

# ---- OTU matrix (rows = samples, cols = SHs) --------------------------------
otu_matrix <- otu_table %>%
  as.data.frame() %>%
  column_to_rownames("SHcode") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  filter(sample %in% META1$sample) %>%
  column_to_rownames("sample")

otu_matrix <- otu_matrix[, colSums(otu_matrix) > 0, drop = FALSE]

cat("Samples remaining after manual exclusion:", nrow(otu_matrix), "\n")
cat("SHs/OTUs remaining after dropping zeros:", ncol(otu_matrix), "\n")

# ---- Per-sample abundance thresholding --------------------------------------
otu_matrix_thr <- threshold_by_sample(otu_matrix, min_reads = 10, drop_empty = FALSE)
abundance_thr  <- colSums(otu_matrix_thr)
keep_otus      <- abundance_thr >= 10
otu_matrix_filt <- otu_matrix_thr[, keep_otus, drop = FALSE]
otu_matrix_filt <- otu_matrix_filt[, colSums(otu_matrix_filt) > 0, drop = FALSE]

cat("Remaining OTUs after taxonomy filter:", ncol(otu_matrix_filt), "\n")

richness_filt <- rowSums(otu_matrix_thr > 0)
reads_filt    <- rowSums(otu_matrix_thr)

META1 <- META1 %>%
  filter(sample %in% rownames(otu_matrix_thr)) %>%
  mutate(
    richness_filt = richness_filt[sample],
    reads_filt    = reads_filt[sample],
    log_reads     = log(reads_filt)
  )

# ---- Save core objects -------------------------------------------------------
saveRDS(
  list(
    META1 = META1,
    tax = tax,
    otu_matrix = otu_matrix,
    otu_matrix_thr = otu_matrix_thr,
    otu_matrix_filt = otu_matrix_filt
  ),
  file = out_rds
)

write_session_info(out_sess)

# ----------------- QUICK RELOAD ----------------------------------------------
obj_list <- readRDS(out_rds)
META1 <- obj_list$META1
tax <- obj_list$tax
otu_matrix <- obj_list$otu_matrix
otu_matrix_thr <- obj_list$otu_matrix_thr
otu_matrix_filt <- obj_list$otu_matrix_filt

