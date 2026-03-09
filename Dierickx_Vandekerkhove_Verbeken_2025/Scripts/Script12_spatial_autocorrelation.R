# =============================================================================
# Script: Script12_spatial_mantel_tree_level.R
# Project: Deadwood fungal community analyses
# Purpose: Tree-level spatial distance decay using Mantel tests
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
#
# Required helper functions from Script0_utils.R:
# - ensure_dir
# - assert_objects
# - d_ait_or_hel
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(sf)
  library(vegan)
})

setwd(".")
set.seed(SEED_GLOBAL)

assert_objects(c("META1", "otu_matrix_filt"))
NPERM  <- 9999L
CRS_LL <- 4326
CRS_M  <- 31370
DW_TYPE_KEEP <- "LOG"   # set to NULL for all samples

ensure_dir("tables")
ensure_dir(file.path("tables", "SPATIAL_LOG"))

# =============================================================================
# 1) Read location table
# =============================================================================
loc <- META1 %>%
  select(sample, dw_type2, XCOORD, YCOORD) %>% 
  rename(xcoord=XCOORD,
         ycoord=YCOORD)
  

# expected cleaned names from XCOORD / YCOORD -> xcoord / ycoord
stopifnot(all(c("sample", "xcoord", "ycoord") %in% names(loc)))

# =============================================================================
# 2) Join metadata and align OTU matrix
# =============================================================================
meta_loc <- META1 %>%
  transmute(
    sample  = as.character(sample),
    natman  = as.character(natman),
    dw_type = as.character(dw_type2)
  )

dat <- loc %>%
  mutate(sample = as.character(sample)) %>%
  left_join(meta_loc, by = "sample") %>%
  filter(!is.na(natman))

if (!is.null(DW_TYPE_KEEP)) {
  dat <- dat %>% filter(dw_type == DW_TYPE_KEEP)
  cat("\nApplied dw_type filter:", DW_TYPE_KEEP, "\n")
}

keep <- intersect(dat$sample, rownames(otu_matrix_filt))
dat  <- dat %>%
  filter(sample %in% keep) %>%
  arrange(sample)

X <- otu_matrix_filt[dat$sample, , drop = FALSE]
X <- as.matrix(X)
X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]

stopifnot(identical(rownames(X), dat$sample))

nat <- dat$natman
tree_ids <- sort(unique(nat))

cat("\nSamples:", nrow(dat), "\n")
cat("Trees:", length(tree_ids), "\n")
cat("SHs retained:", ncol(X), "\n")

# =============================================================================
# 3) Coordinates and natman centroids
# =============================================================================
pts_ll <- sf::st_as_sf(
  dat,
  coords = c("xcoord", "ycoord"),
  crs = CRS_LL,
  remove = FALSE
)

pts_m <- sf::st_transform(pts_ll, CRS_M)
coords_m <- sf::st_coordinates(pts_m)

tree_xy_m <- tibble(natman = nat) %>%
  bind_cols(
    as_tibble(coords_m, .name_repair = "minimal") %>%
      setNames(c("x", "y"))
  ) %>%
  group_by(natman) %>%
  summarise(
    x = mean(x, na.rm = TRUE),
    y = mean(y, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(natman)

stopifnot(identical(tree_xy_m$natman, tree_ids))

D_geo_tree <- dist(tree_xy_m[, c("x", "y")], method = "euclidean")
attr(D_geo_tree, "Labels") <- tree_xy_m$natman

# =============================================================================
# 4) Build sample weights within natman
# =============================================================================
w_tbl <- dat %>%
  transmute(sample = sample, natman = natman) %>%
  group_by(natman) %>%
  mutate(n_in_tree = n()) %>%
  ungroup()

w <- 1 / w_tbl$n_in_tree
w <- w / ave(w, nat, FUN = sum)

stopifnot(all(is.finite(w)))
stopifnot(all(abs(tapply(w, nat, sum) - 1) < 1e-8))

# =============================================================================
# 5) Collapse to tree level for robust Aitchison
# =============================================================================
row_totals <- rowSums(X, na.rm = TRUE)
stopifnot(all(row_totals > 0))

P <- sweep(X, 1, row_totals, "/")
P_w <- P * w

P_tree_sum <- rowsum(P_w, group = nat, reorder = TRUE)
P_tree <- sweep(P_tree_sum, 1, rowSums(P_tree_sum), "/")

stopifnot(identical(rownames(P_tree), tree_ids))

D_ait_tree <- d_ait_or_hel(P_tree)
D_ait_tree <- fix_dist_labels(D_ait_tree, samps = tree_ids)

# =============================================================================
# 6) Collapse to tree level for Jaccard
# =============================================================================
X_pa <- (X > 0) * 1L
pa_sum <- rowsum(X_pa, group = nat, reorder = TRUE)
pa_tree <- (pa_sum > 0) * 1L

stopifnot(identical(rownames(pa_tree), tree_ids))

D_jac_tree <- vegan::vegdist(pa_tree, method = "jaccard", binary = TRUE)
D_jac_tree <- fix_dist_labels(D_jac_tree, samps = tree_ids)

# =============================================================================
# 7) Mantel tests
# =============================================================================
mt_A <- vegan::mantel(
  D_ait_tree,
  D_geo_tree,
  method = "spearman",
  permutations = NPERM
)
mt_J <- vegan::mantel(
  D_jac_tree,
  D_geo_tree,
  method = "spearman",
  permutations = NPERM
)

out <- tibble(
  level    = "tree",
  metric   = c("robust.aitchison", "jaccard"),
  mantel_r = c(as.numeric(mt_A$statistic), as.numeric(mt_J$statistic)),
  mantel_p = c(as.numeric(mt_A$signif),    as.numeric(mt_J$signif)),
  n_trees  = nrow(tree_xy_m),
  n_samples = nrow(dat),
  dw_type_filter = ifelse(is.null(DW_TYPE_KEEP), "ALL", DW_TYPE_KEEP)
)

cat("\nTree-level Mantel distance decay:\n")
print(out)

readr::write_csv2(
  out,
  file.path("tables", "SPATIAL_LOG", "spatial_mantel_tree_level_weighted_simple.csv")
)

