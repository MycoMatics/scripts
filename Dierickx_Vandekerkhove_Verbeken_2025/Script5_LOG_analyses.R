# =============================================================================
# Script 5 — INTRA (within-deadwood) analyses by scope tied by natman
# Scopes:
#   - "LOG"
#   - "LOG_aFWD"
#   - "LOG_aFWD_fFWD"
#   - "LOG_aFWD_fFWD_SNAG"
#
# For each scope:
#   - Build metadata & OTU subset (rq == "INTRA", common natman)
#   - Alpha diversity (by microhab group)
#   - PERMANOVA (blocked by natman) across distances
#   - Partial dbRDA + Euler (Tree vs. Microhab complex)
#   - Dispersion checks (PERMDISP)
#   - Indicator species (indicspecies::multipatt)
#   - Alpha richness GLMM (mixed-effects; random intercept natman)
#
# Depends on Script 1 objects: META1, otu_matrix_filt, tax, palettes, threshold
# make sure functions are preloaded from Utils.R
# =============================================================================
source("scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(janitor)
  library(permute)
  library(lme4)
  library(emmeans)
  library(car)
  library(broom.mixed)
  library(performance)
  library(indicspecies)
  library(eulerr)
  library(forcats)
  library(grid)
})

set.seed(42)
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/meta_SH/R_data")

# -----------------------------------------------------------------------------
# 0) Guards, dirs
# -----------------------------------------------------------------------------
stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"))
stopifnot(exists("threshold"))
stopifnot(exists("build_intra_meta"))
stopifnot(exists("run_scope"))

cat("You are using SH threshold ", threshold, "%")

dir.create("plots/INTRA", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1) Run all scopes
# -----------------------------------------------------------------------------
scopes <- c("LOG", "LOG_aFWD", "LOG_aFWD_fFWD", "LOG_aFWD_fFWD_SNAG")
scopes <- "LOG"
# load build_intra_meta_ds_at_drill for "LOG", "LOG_aFWD",
# load build_intra_meta for "LOG_aFWD_fFWD_SNAG", "LOG_aFWD_fFWD_SNAG"
invisible(lapply(scopes, run_scope))
# Debug: can't really use LOG_aFWD_fFWD (_SNAG) because only 3 trees makes everything overdispersed  ...
invisible(lapply(scopes, run_scope_ds))
