#!/usr/bin/env Rscript
# Four-dataset A1-filtered community sensitivity analysis for equal-ITS SINTAX.
# CD2 weights confidently assigned D1 consensuses by original `ubs` family size.
# Preserves paired barcode tests, 999 permutations and metric-specific NMDS k.

cat("\n")
cat("=====================================================================\n")
cat("SINTAX CD1-CD4 A1-FILTERED ALL-IN-ONE ANALYSIS\n")
cat("=====================================================================\n\n")

# PART 1. APPLY PREDEFINED A1 SAMPLE EXCLUSIONS

# Apply the same predefined problematic-sample exclusions used in the A1
# A1 source list from:
# Original A1 comment:
# barcode convention (41MI, 24P1), intersects with the current paired table,

ROOT <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(ROOT)) stop("Set PROJECT_ROOT to the UMIvsRAW project root.")
ROOT <- normalizePath(ROOT, mustWork = TRUE)
TABLE_DIR <- file.path(
  ROOT, "REVIEW", "SINTAX_equal_ITS_v1", "analysis", "tables"
)

INPUT <- file.path(
  TABLE_DIR,
  "SINTAX_CD1-CD4_long_SH_conf080_4datasets.csv"
)

OUTPUT <- file.path(
  TABLE_DIR,
  "SINTAX_CD1-CD4_long_SH_conf080_A1filtered_4datasets.csv"
)

AUDIT <- file.path(
  TABLE_DIR,
  "SINTAX_CD1-CD4_A1_exclusion_audit_4datasets.csv"
)

SUMMARY <- file.path(
  TABLE_DIR,
  "SINTAX_CD1-CD4_A1_exclusion_summary_4datasets.txt"
)

if (!file.exists(INPUT)) {
  stop("Input not found: ", INPUT)
}

# Exact predefined A1 exclusion list.
A1_exclude <- c(
  "P3_50", "P1_24", "P2_93", "P2_86", "P3_63", "P1_66", "P3_47",
  "MI_33", "MI_41", "MI_42", "P1_19", "MI_45", "P1_56", "P2_84",
  "P3_42", "P3_41", "P3_37", "P3_82"
)

# Convert A1 sample labels to current barcode labels:
#   MI_41 -> 41MI
a1_to_current <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  vapply(
    parts,
    function(z) {
      if (length(z) != 2L) {
        stop("Unexpected A1 sample ID: ", paste(z, collapse = "_"))
      }
      paste0(z[2], z[1])
    },
    character(1)
  )
}

exclude_map <- data.frame(
  A1_sample = A1_exclude,
  current_barcode = a1_to_current(A1_exclude),
  stringsAsFactors = FALSE
)

d <- utils::read.csv2(
  INPUT,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c("dataset", "Barcode", "SH_code", "count")
missing <- setdiff(required, names(d))
if (length(missing)) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

d$dataset <- as.character(d$dataset)
d$Barcode <- as.character(d$Barcode)

datasets <- c("CD1", "CD2", "CD3", "CD4")
if (!all(datasets %in% unique(d$dataset))) {
  stop("Expected CD1-CD4 in input.")
}

# Validate the current paired design before filtering.
sample_pairs <- unique(d[, c("dataset", "Barcode")])
sets_before <- lapply(
  datasets,
  function(ds) sort(unique(sample_pairs$Barcode[sample_pairs$dataset == ds]))
)
names(sets_before) <- datasets

if (!all(vapply(sets_before[-1], identical, logical(1), sets_before[[1]]))) {
  stop("Input is not balanced across CD1-CD4 before applying A1 exclusions.")
}

paired_before <- sets_before[[1]]

exclude_map$present_in_paired_input <- exclude_map$current_barcode %in% paired_before
exclude_map$action <- ifelse(
  exclude_map$present_in_paired_input,
  "EXCLUDE",
  "NOT_PRESENT_IN_CURRENT_PAIRED_SET"
)

matched_exclusions <- exclude_map$current_barcode[
  exclude_map$present_in_paired_input
]

# Filter the same biological samples from every dataset.
d_filt <- d[!d$Barcode %in% matched_exclusions, , drop = FALSE]

sample_pairs_after <- unique(d_filt[, c("dataset", "Barcode")])
sets_after <- lapply(
  datasets,
  function(ds) sort(unique(sample_pairs_after$Barcode[sample_pairs_after$dataset == ds]))
)
names(sets_after) <- datasets

if (!all(vapply(sets_after[-1], identical, logical(1), sets_after[[1]]))) {
  stop("Filtering broke the paired design.")
}

paired_after <- sets_after[[1]]

expected_after <- length(paired_before) - length(unique(matched_exclusions))
if (length(paired_after) != expected_after) {
  stop(
    "Unexpected retained barcode count: expected ",
    expected_after, ", observed ", length(paired_after)
  )
}

# Explicitly audit the preserved A1 exclusion set and historical barcode cohort.
expected_excluded <- sort(c("19P1", "24P1", "33MI", "41MI", "42MI", "45MI", "56P1"))
if (length(paired_before) != 97L ||
    !identical(sort(unique(matched_exclusions)), expected_excluded) ||
    length(paired_after) != 90L) {
  stop("Unexpected A1 cohort or exclusions; before=", length(paired_before),
       "; excluded=", paste(sort(unique(matched_exclusions)), collapse=","),
       "; after=", length(paired_after))
}

cat("Validated A1-filtered cohort: 90 barcodes after the predefined exclusions.\n")
cat("Equal-ITS four-dataset design: CD1, CD2, CD3, CD4.\n")
cat("CD2 abundances derive from original D1 ubs metadata.\n")

# Confirm that every matched exclusion is absent from every dataset.
still_present <- intersect(unique(d_filt$Barcode), matched_exclusions)
if (length(still_present)) {
  stop(
    "Excluded barcode(s) still present: ",
    paste(still_present, collapse = ", ")
  )
}

utils::write.csv2(
  d_filt,
  OUTPUT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

utils::write.csv2(
  exclude_map,
  AUDIT,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

summary_lines <- c(
  "SINTAX CD1-CD4 A1 sample-exclusion audit",
  "========================================",
  paste0("input_paired_barcodes\t", length(paired_before)),
  paste0("A1_predefined_exclusion_IDs\t", length(A1_exclude)),
  paste0("A1_exclusions_present_in_current_paired_set\t", length(unique(matched_exclusions))),
  paste0(
    "matched_current_barcodes\t",
    if (length(matched_exclusions)) paste(sort(unique(matched_exclusions)), collapse = ",") else "none"
  ),
  paste0("retained_paired_barcodes\t", length(paired_after)),
  paste0("input_rows\t", nrow(d)),
  paste0("output_rows\t", nrow(d_filt)),
  "paired_design_after_filtering\tPASS",
  paste0("output\t", OUTPUT),
  paste0("audit\t", AUDIT)
)

writeLines(summary_lines, SUMMARY)

cat("===== A1 EXCLUSION MAPPING =====\n")
print(exclude_map, row.names = FALSE)

cat("\n===== MATCHED EXCLUSIONS IN CURRENT PAIRED DATA =====\n")
if (length(matched_exclusions)) {
  print(sort(unique(matched_exclusions)))
} else {
  cat("none\n")
}

cat("\n===== FILTERING COMPLETE =====\n")
cat("Paired barcodes before: ", length(paired_before), "\n", sep = "")
cat("Matched A1 exclusions:  ", length(unique(matched_exclusions)), "\n", sep = "")
cat("Paired barcodes after:  ", length(paired_after), "\n", sep = "")
cat("Filtered table:\n  ", OUTPUT, "\n", sep = "")
cat("Audit:\n  ", AUDIT, "\n", sep = "")
cat("Summary:\n  ", SUMMARY, "\n", sep = "")

# PART 2. RUN COMPLETE A1-FILTERED COMMUNITY SENSITIVITY ANALYSIS

# SINTAX CD1-CD4 community sensitivity analysis, A1-filtered
# Base R + vegan only. Analysis tables and manuscript-style plots are generated.
#   90 A1-filtered barcodes paired across CD1-CD4 (97 before exclusions).
#   Global and pairwise PERMANOVA are restricted within Barcode.

options(stringsAsFactors = FALSE)
options(error = function() {
  traceback(3)
  quit(save = "no", status = 1, runLast = FALSE)
})

BASE <- file.path(ROOT, "REVIEW", "SINTAX_equal_ITS_v1", "analysis", "community_sensitivity")
TABLE_DIR <- file.path(ROOT, "REVIEW", "SINTAX_equal_ITS_v1", "analysis", "tables")
OUT_DIR <- file.path(BASE, "CD1-CD4_analysis_conf080_A1filtered_4datasets")

INPUT <- file.path(TABLE_DIR, "SINTAX_CD1-CD4_long_SH_conf080_A1filtered_4datasets.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(INPUT)) {
  stop("Input file not found: ", INPUT)
}

# R PACKAGE SETUP
required_packages <- c("vegan", "ggplot2", "tidyverse", "cowplot")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(tidyverse)
  library(cowplot)
})
cat("R version: ", R.version.string, "\n", sep = "")
for (pkg in required_packages) {
  cat(pkg, " version: ", as.character(utils::packageVersion(pkg)), "\n", sep = "")
}

set.seed(42)

N_PERM <- 999L

# Metric-specific NMDS dimensionality:
#   Bray-Curtis      k = 3
NMDS_K_BY_METRIC <- c(
  "Bray-Curtis" = 3L,
  "Hellinger" = 3L,
  "Robust-Aitchison" = 2L
)

NMDS_TRYMAX <- 1000L
SH_RE <- "^SH[0-9]{7}\\.10FU$"

# Format NMDS stress for figures.
# scientific notation so genuinely near-zero stress is not displayed as 0.000.
format_stress_value <- function(x) {
  if (length(x) != 1L || !is.finite(x)) {
    return("NA")
  }
  
  if (abs(x) < 0.001) {
    return(formatC(x, format = "e", digits = 2))
  }
  
  sprintf("%.3f", x)
}

write_csv2 <- function(x, path) {
  utils::write.csv2(
    x,
    file = path,
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

normalize_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("^barcode", "", x, ignore.case = TRUE)
  gsub("_", "", x, fixed = TRUE)
}

make_pair_name <- function(a, b) paste(a, b, sep = "_vs_")

p_from_table <- function(tab, row = 1L) {
  pcol <- grep("^Pr\\(", names(tab), value = TRUE)
  if (!length(pcol)) return(NA_real_)
  as.numeric(tab[row, pcol[[1]]])
}

set_dist_labels <- function(d, labels, metric) {
  if (!inherits(d, "dist")) {
    stop(metric, " did not return a dist object.")
  }
  
  if (length(labels) != attr(d, "Size")) {
    stop(
      metric, " distance size/label mismatch: Size=", attr(d, "Size"),
      ", labels=", length(labels)
    )
  }
  
  # vegan 2.7.x can return robust Aitchison dist objects without preserving
  # the input row names. Downstream character subsetting then fails.
  attr(d, "Labels") <- as.character(labels)
  
  observed <- attr(d, "Labels")
  if (is.null(observed) || !identical(observed, as.character(labels))) {
    stop(metric, " distance labels could not be assigned correctly.")
  }
  
  d
}

make_distance <- function(mat, metric) {
  sample_labels <- rownames(mat)
  
  if (is.null(sample_labels) || anyNA(sample_labels) ||
      any(!nzchar(sample_labels)) || anyDuplicated(sample_labels)) {
    stop("Community matrix must have unique, non-empty row names.")
  }
  
  if (metric == "Bray-Curtis") {
    rs <- rowSums(mat)
    if (any(rs <= 0)) stop("Zero-sum sample encountered for Bray-Curtis.")
    rel <- sweep(mat, 1, rs, "/")
    d <- vegan::vegdist(rel, method = "bray")
    return(set_dist_labels(d, sample_labels, metric))
  }
  
  if (metric == "Hellinger") {
    rs <- rowSums(mat)
    if (any(rs <= 0)) stop("Zero-sum sample encountered for Hellinger.")
    hell <- sqrt(sweep(mat, 1, rs, "/"))
    d <- stats::dist(hell, method = "euclidean")
    return(set_dist_labels(d, sample_labels, metric))
  }
  
  if (metric == "Robust-Aitchison") {
    ans <- try(
      vegan::vegdist(mat, method = "robust.aitchison"),
      silent = TRUE
    )
    if (inherits(ans, "try-error")) {
      stop(
        "vegan::vegdist(method='robust.aitchison') failed.\n",
        "vegan version: ", as.character(utils::packageVersion("vegan")), "\n",
        "Error: ", as.character(ans)
      )
    }
    return(set_dist_labels(ans, sample_labels, metric))
  }
  
  stop("Unknown metric: ", metric)
}

subset_dist <- function(d, ids) {
  labs <- attr(d, "Labels")
  
  if (is.null(labs)) {
    stop(
      "Distance object has no Labels attribute; cannot subset safely by sample ID."
    )
  }
  
  missing_ids <- setdiff(ids, labs)
  if (length(missing_ids)) {
    stop(
      "Distance object is missing ", length(missing_ids),
      " requested sample ID(s). First missing IDs: ",
      paste(head(missing_ids, 10), collapse = ", "),
      "\nAvailable labels begin: ",
      paste(head(labs, 10), collapse = ", ")
    )
  }
  
  m <- as.matrix(d)
  sub <- m[ids, ids, drop = FALSE]
  out <- stats::as.dist(sub)
  
  # as.dist() should preserve these, but enforce them explicitly so the
  # downstream NMDS and permutation code always has stable sample labels.
  attr(out, "Labels") <- ids
  out
}

validate_dist_for_nmds <- function(d, k, label = "distance") {
  if (!inherits(d, "dist")) {
    stop(label, " is not a dist object.")
  }
  
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop(label, " has invalid NMDS dimensionality k=", k)
  }
  
  vals <- as.numeric(d)
  
  if (!length(vals)) {
    stop(label, " contains no pairwise distances.")
  }
  if (any(!is.finite(vals))) {
    stop(label, " contains NA/NaN/Inf distances.")
  }
  if (any(vals < 0)) {
    stop(label, " contains negative distances.")
  }
  
  n_obj <- attr(d, "Size")
  if (is.null(n_obj) || n_obj <= k + 1L) {
    stop(
      label, " has too few objects for k=", k,
      ": n=", n_obj
    )
  }
  
  invisible(TRUE)
}

run_nmds <- function(d, k, label = "NMDS") {
  validate_dist_for_nmds(d, k = k, label = label)
  
  # metaMDS is supplied an already-computed dissimilarity object.
  # original community matrix and are not needed here.
  fit <- tryCatch(
    vegan::metaMDS(
      d,
      k = k,
      trymax = NMDS_TRYMAX,
      autotransform = FALSE,
      noshare = FALSE,
      wascores = FALSE,
      expand = FALSE,
      trace = FALSE
    ),
    error = function(e) {
      stop(
        label, " failed in metaMDS with k=", k, ": ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  
  scr <- vegan::scores(
    fit,
    display = "sites",
    choices = seq_len(k)
  )
  
  if (is.null(scr) || nrow(scr) != attr(d, "Size")) {
    stop(
      label, " returned an unexpected site-score matrix: expected ",
      attr(d, "Size"), " rows; observed ",
      if (is.null(scr)) 0L else nrow(scr),
      "."
    )
  }
  
  # Store explicitly for later reporting and Procrustes extraction.
  fit$k_used <- as.integer(k)
  
  fit
}

extract_adonis_row <- function(tab, term, metric, scope, test_name) {
  rn <- rownames(tab)
  i <- match(term, rn)
  
  # adonis2() returns a row named "Model" when by = NULL. For a one-predictor
  # explicitly use by = "terms" below, but retain this fallback for robustness.
  if (is.na(i) && "Model" %in% rn) {
    i <- match("Model", rn)
  }
  
  if (is.na(i)) {
    stop(
      "Term '", term, "' not found in adonis2 output. Rows were: ",
      paste(rn, collapse = ", ")
    )
  }
  
  data.frame(
    scope = scope,
    metric = metric,
    test = test_name,
    term = term,
    Df = if ("Df" %in% names(tab)) as.numeric(tab[i, "Df"]) else NA_real_,
    SumOfSqs = if ("SumOfSqs" %in% names(tab)) as.numeric(tab[i, "SumOfSqs"]) else NA_real_,
    R2 = if ("R2" %in% names(tab)) as.numeric(tab[i, "R2"]) else NA_real_,
    F = if ("F" %in% names(tab)) as.numeric(tab[i, "F"]) else NA_real_,
    p = p_from_table(tab, i),
    stringsAsFactors = FALSE
  )
}

cat("SINTAX CD1-CD4 community sensitivity analysis\n")
cat("============================================\n")
cat("Input: ", INPUT, "\n", sep = "")
cat("vegan version: ", as.character(utils::packageVersion("vegan")), "\n", sep = "")
cat("Permutations: ", N_PERM, "\n\n", sep = "")

# 1. Read and validate the balanced long table.
long <- utils::read.csv2(
  INPUT,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c("dataset", "Barcode", "SH_code", "count")
missing <- setdiff(required, names(long))
if (length(missing)) {
  stop("Missing required column(s): ", paste(missing, collapse = ", "))
}

long$dataset <- as.character(long$dataset)
long$Barcode <- normalize_barcode(long$Barcode)
long$SH_code <- as.character(long$SH_code)
long$count <- safe_num(long$count)

if (anyNA(long$count) || any(long$count <= 0) ||
    any(abs(long$count - round(long$count)) > 1e-8)) {
  stop("count must contain positive integer-like values only.")
}

if (any(!grepl(SH_RE, long$SH_code))) {
  bad <- unique(long$SH_code[!grepl(SH_RE, long$SH_code)])
  stop("Invalid SH code(s): ", paste(head(bad, 10), collapse = ", "))
}

datasets <- c("CD1", "CD2", "CD3", "CD4")
if (!all(datasets %in% unique(long$dataset))) {
  stop(
    "Expected CD1-CD4. Found: ",
    paste(sort(unique(long$dataset)), collapse = ", ")
  )
}

long <- long[long$dataset %in% datasets, , drop = FALSE]

# Collapse any accidental duplicate cells defensively.
agg <- stats::aggregate(
  count ~ dataset + Barcode + SH_code,
  data = long,
  FUN = sum
)

sample_grid <- unique(agg[, c("dataset", "Barcode")])
n_by_dataset <- table(sample_grid$dataset)

if (!all(n_by_dataset[datasets] == n_by_dataset[datasets][1])) {
  print(n_by_dataset)
  stop("The input is not balanced across datasets.")
}

barcode_sets <- lapply(datasets, function(ds) {
  sort(unique(sample_grid$Barcode[sample_grid$dataset == ds]))
})
names(barcode_sets) <- datasets

if (!all(vapply(barcode_sets[-1], identical, logical(1), barcode_sets[[1]]))) {
  stop("The four datasets do not contain identical barcode sets.")
}

barcodes <- barcode_sets[[1]]
N_BARCODE <- length(barcodes)
if (N_BARCODE != 90L) stop("Expected 90 A1-matched barcodes, got ", N_BARCODE)

if (N_BARCODE < 2L) {
  stop("Too few paired barcodes after A1 filtering: ", N_BARCODE)
}

cat("Validated paired design: 4 datasets x ", N_BARCODE,
    " barcodes = ", 4L * N_BARCODE, " community profiles\n", sep = "")
cat("Unique SHs across paired table: ", length(unique(agg$SH_code)), "\n\n", sep = "")

# 2. Build sample x SH matrix.
agg$obs_id <- paste(agg$dataset, agg$Barcode, sep = "|")

comm <- xtabs(
  count ~ obs_id + SH_code,
  data = agg,
  drop.unused.levels = TRUE
)

comm <- as.matrix(comm)
storage.mode(comm) <- "numeric"

meta <- unique(agg[, c("obs_id", "dataset", "Barcode")])
rownames(meta) <- meta$obs_id
meta <- meta[rownames(comm), , drop = FALSE]

meta$dataset <- factor(meta$dataset, levels = datasets)
meta$Barcode <- factor(meta$Barcode, levels = barcodes)

if (any(rowSums(comm) <= 0)) stop("At least one sample has zero total abundance.")
if (anyDuplicated(meta$obs_id)) stop("Duplicate observation IDs in metadata.")
if (!identical(rownames(comm), meta$obs_id)) stop("Community/meta row mismatch.")

# 3. SH overlap and abundance of newly detected CD4 taxa.
presence <- lapply(datasets, function(ds) {
  unique(agg$SH_code[agg$dataset == ds])
})
names(presence) <- datasets

overlap_rows <- do.call(
  rbind,
  lapply(datasets, function(ds) {
    data.frame(
      dataset = ds,
      n_SH = length(presence[[ds]]),
      stringsAsFactors = FALSE
    )
  })
)

focus_sets <- c("CD1", "CD3", "CD4")
union_13 <- union(presence$CD1, presence$CD3)
cd4_unique_vs_13 <- setdiff(presence$CD4, union_13)
shared_134 <- Reduce(intersect, presence[focus_sets])

cd4_rows <- agg$dataset == "CD4"
cd4_unique_reads <- sum(
  agg$count[cd4_rows & agg$SH_code %in% cd4_unique_vs_13]
)
cd4_total_reads <- sum(agg$count[cd4_rows])

overlap_summary <- data.frame(
  statistic = c(
    "CD1_SH",
    "CD3_SH",
    "CD4_SH",
    "shared_CD1_CD3_CD4_SH",
    "CD4_unique_vs_union_CD1_CD3_SH",
    "CD4_unique_vs_union_CD1_CD3_reads",
    "CD4_total_reads_A1filtered90",
    "CD4_unique_vs_union_CD1_CD3_read_fraction"
  ),
  value = c(
    length(presence$CD1),
    length(presence$CD3),
    length(presence$CD4),
    length(shared_134),
    length(cd4_unique_vs_13),
    cd4_unique_reads,
    cd4_total_reads,
    cd4_unique_reads / cd4_total_reads
  ),
  stringsAsFactors = FALSE
)

write_csv2(
  overlap_rows,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_SH_counts_conf080.csv")
)
write_csv2(
  overlap_summary,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_SH_overlap_summary_conf080.csv")
)
write_csv2(
  data.frame(SH_code = sort(cd4_unique_vs_13), stringsAsFactors = FALSE),
  file.path(OUT_DIR, "SINTAX_CD4_unique_vs_CD1_CD3_SH_conf080.csv")
)

# 4. Alpha richness: CD1, CD3, CD4 only, matching manuscript comparison.
alpha_ds <- c("CD1", "CD3", "CD4")

rich <- stats::aggregate(
  SH_code ~ dataset + Barcode,
  data = agg[agg$dataset %in% alpha_ds, , drop = FALSE],
  FUN = function(x) length(unique(x))
)
names(rich)[names(rich) == "SH_code"] <- "richness"
rich$dataset <- factor(rich$dataset, levels = alpha_ds)

rich_summary <- do.call(
  rbind,
  lapply(alpha_ds, function(ds) {
    x <- rich$richness[rich$dataset == ds]
    data.frame(
      dataset = ds,
      n = length(x),
      mean = mean(x),
      sd = stats::sd(x),
      median = stats::median(x),
      IQR = stats::IQR(x),
      min = min(x),
      max = max(x),
      stringsAsFactors = FALSE
    )
  })
)

rich_wide <- reshape(
  rich,
  idvar = "Barcode",
  timevar = "dataset",
  direction = "wide"
)

rich_cols <- paste0("richness.", alpha_ds)
if (!all(rich_cols %in% names(rich_wide))) {
  stop("Could not construct complete paired richness table.")
}
if (anyNA(rich_wide[, rich_cols])) {
  stop("Missing richness values in paired alpha-diversity table.")
}

fried_alpha <- stats::friedman.test(
  as.matrix(rich_wide[, rich_cols, drop = FALSE])
)

alpha_global <- data.frame(
  test = "Friedman",
  datasets = paste(alpha_ds, collapse = ","),
  n_barcodes = N_BARCODE,
  statistic = unname(fried_alpha$statistic),
  df = unname(fried_alpha$parameter),
  p = fried_alpha$p.value,
  kendalls_W = unname(fried_alpha$statistic) /
    (N_BARCODE * (length(alpha_ds) - 1L)),
  stringsAsFactors = FALSE
)

alpha_pairs <- combn(alpha_ds, 2, simplify = FALSE)
alpha_pairwise <- do.call(
  rbind,
  lapply(alpha_pairs, function(z) {
    a <- rich_wide[[paste0("richness.", z[1])]]
    b <- rich_wide[[paste0("richness.", z[2])]]
    wt <- suppressWarnings(stats::wilcox.test(
      a, b,
      paired = TRUE,
      exact = FALSE,
      correct = FALSE
    ))
    data.frame(
      dataset1 = z[1],
      dataset2 = z[2],
      n = length(a),
      median_difference_dataset1_minus_dataset2 = stats::median(a - b),
      V = unname(wt$statistic),
      p = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
)
alpha_pairwise$p_Bonferroni <- stats::p.adjust(alpha_pairwise$p, method = "bonferroni")

write_csv2(rich, file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_richness_by_barcode_conf080.csv"))
write_csv2(rich_summary, file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_richness_summary_conf080.csv"))
write_csv2(alpha_global, file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_Friedman_conf080.csv"))
write_csv2(alpha_pairwise, file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_pairwise_Wilcoxon_conf080.csv"))

cat("Alpha richness means (SH/barcode):\n")
print(rich_summary[, c("dataset", "mean", "sd", "median")], row.names = FALSE)
cat("\n")

# 5. Distance matrices.
metrics <- c("Bray-Curtis", "Hellinger", "Robust-Aitchison")

scopes <- list(
  CD1_CD3 = c("CD1", "CD2", "CD3"),
  CD1_CD4 = c("CD1", "CD2", "CD3", "CD4")
)

merged_dist <- setNames(vector("list", length(metrics)), metrics)
for (metric in metrics) {
  cat("Computing merged ", metric, " distance...\n", sep = "")
  merged_dist[[metric]] <- make_distance(comm, metric)
}

cat("\nDistance-label validation:\n")
for (metric in metrics) {
  labs <- attr(merged_dist[[metric]], "Labels")
  cat(
    "  ", metric,
    ": Size=", attr(merged_dist[[metric]], "Size"),
    "; labels=", length(labs),
    "; identical_to_comm=", identical(labs, rownames(comm)),
    "\n",
    sep = ""
  )
}
cat("\n")

dataset_dist <- setNames(vector("list", length(metrics)), metrics)
for (metric in metrics) {
  dataset_dist[[metric]] <- setNames(vector("list", length(datasets)), datasets)
  for (ds in datasets) {
    ids <- paste(ds, barcodes, sep = "|")
    submat <- comm[ids, , drop = FALSE]
    rownames(submat) <- barcodes
    dataset_dist[[metric]][[ds]] <- make_distance(submat, metric)
  }
}
cat("\n")

cat("NMDS dimensionality by metric:\n")
print(
  data.frame(
    metric = names(NMDS_K_BY_METRIC),
    k = as.integer(NMDS_K_BY_METRIC),
    row.names = NULL,
    check.names = FALSE
  ),
  row.names = FALSE
)
cat("\n")

# 6. NMDS stress + Mantel + Procrustes.
#    Dataset-specific NMDS uses k=3 for Bray/Hellinger and k=2 for robust Aitchison.
nmds <- setNames(vector("list", length(metrics)), metrics)
nmds_dataset_stress <- list()

for (mi in seq_along(metrics)) {
  metric <- metrics[[mi]]
  nmds[[metric]] <- setNames(vector("list", length(datasets)), datasets)
  
  k_metric <- unname(NMDS_K_BY_METRIC[[metric]])
  
  for (di in seq_along(datasets)) {
    ds <- datasets[[di]]
    
    if (is.null(dataset_dist[[metric]]) ||
        is.null(dataset_dist[[metric]][[ds]])) {
      stop(
        "Missing dataset-specific distance object for metric=",
        metric, ", dataset=", ds
      )
    }
    
    cat(
      "NMDS dataset fit: ",
      metric, " / ", ds, "\n",
      sep = ""
    )
    
    set.seed(42000 + mi * 100 + di)
    fit <- run_nmds(
      dataset_dist[[metric]][[ds]],
      k = k_metric,
      label = paste("NMDS", metric, ds, sep = " / ")
    )
    nmds[[metric]][[ds]] <- fit
    
    nmds_dataset_stress[[length(nmds_dataset_stress) + 1L]] <- data.frame(
      scope = "single_dataset",
      metric = metric,
      dataset = ds,
      k = k_metric,
      stress = fit$stress,
      converged = fit$converged > 0,
      tries = if (!is.null(fit$tries)) fit$tries else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}
nmds_dataset_stress <- do.call(rbind, nmds_dataset_stress)

scope_nmds <- setNames(vector("list", length(scopes)), names(scopes))
nmds_scope_stress <- list()

for (si in seq_along(scopes)) {
  scope_name <- names(scopes)[si]
  ds_scope <- scopes[[si]]
  ids_scope <- rownames(meta)[meta$dataset %in% ds_scope]
  
  scope_nmds[[scope_name]] <- setNames(
    vector("list", length(metrics)),
    metrics
  )
  
  for (mi in seq_along(metrics)) {
    metric <- metrics[[mi]]
    k_metric <- unname(NMDS_K_BY_METRIC[[metric]])
    d_scope <- subset_dist(merged_dist[[metric]], ids_scope)
    
    cat(
      "NMDS scope fit: ",
      scope_name, " / ", metric, "\n",
      sep = ""
    )
    
    set.seed(42500 + si * 100 + mi)
    fit <- run_nmds(
      d_scope,
      k = k_metric,
      label = paste("NMDS", scope_name, metric, sep = " / ")
    )
    scope_nmds[[scope_name]][[metric]] <- fit
    
    nmds_scope_stress[[length(nmds_scope_stress) + 1L]] <- data.frame(
      scope = scope_name,
      metric = metric,
      datasets = paste(ds_scope, collapse = ","),
      k = k_metric,
      stress = fit$stress,
      converged = fit$converged > 0,
      tries = if (!is.null(fit$tries)) fit$tries else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}
nmds_scope_stress <- do.call(rbind, nmds_scope_stress)

dataset_pairs <- combn(datasets, 2, simplify = FALSE)
concordance <- list()

for (mi in seq_along(metrics)) {
  metric <- metrics[[mi]]
  k_metric <- unname(NMDS_K_BY_METRIC[[metric]])
  
  for (pi in seq_along(dataset_pairs)) {
    z <- dataset_pairs[[pi]]
    d1 <- dataset_dist[[metric]][[z[1]]]
    d2 <- dataset_dist[[metric]][[z[2]]]
    
    set.seed(43000 + mi * 100 + pi)
    mt <- vegan::mantel(
      d1, d2,
      method = "spearman",
      permutations = N_PERM
    )
    
    x_scores <- as.matrix(
      vegan::scores(
        nmds[[metric]][[z[1]]],
        display = "sites",
        choices = seq_len(k_metric)
      )
    )
    y_scores <- as.matrix(
      vegan::scores(
        nmds[[metric]][[z[2]]],
        display = "sites",
        choices = seq_len(k_metric)
      )
    )
    
    # Dataset-specific NMDS distance objects were labelled by barcode, so
    # enforce the same barcode order before Procrustes.
    x_scores <- x_scores[barcodes, , drop = FALSE]
    y_scores <- y_scores[barcodes, , drop = FALSE]
    
    set.seed(44000 + mi * 100 + pi)
    pr <- vegan::protest(
      x_scores,
      y_scores,
      permutations = N_PERM
    )
    
    concordance[[length(concordance) + 1L]] <- data.frame(
      metric = metric,
      k = k_metric,
      dataset1 = z[1],
      dataset2 = z[2],
      Mantel_spearman_r = unname(mt$statistic),
      Mantel_p = mt$signif,
      Procrustes_correlation = unname(pr$t0),
      Procrustes_p = pr$signif,
      stringsAsFactors = FALSE
    )
  }
}

concordance <- do.call(rbind, concordance)
concordance$Mantel_p_BH <- ave(
  concordance$Mantel_p,
  concordance$metric,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)
concordance$Procrustes_p_BH <- ave(
  concordance$Procrustes_p,
  concordance$metric,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)

write_csv2(
  nmds_scope_stress,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_NMDS_stress_conf080.csv")
)
write_csv2(
  nmds_dataset_stress,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_NMDS_dataset_stress_for_Procrustes_conf080.csv")
)
write_csv2(
  concordance,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_Mantel_Procrustes_conf080.csv")
)
saveRDS(
  list(
    scope_nmds = scope_nmds,
    dataset_nmds = nmds,
    k_by_metric = NMDS_K_BY_METRIC
  ),
  file.path(OUT_DIR, "SINTAX_CD1-CD4_NMDS_fitted_objects.rds")
)
# 7. Global and pairwise stratified PERMANOVA.
global_perm <- list()
pairwise_perm <- list()

for (si in seq_along(scopes)) {
  scope_name <- names(scopes)[si]
  ds_scope <- scopes[[si]]
  
  ids_scope <- rownames(meta)[meta$dataset %in% ds_scope]
  meta_scope <- droplevels(meta[ids_scope, , drop = FALSE])
  meta_scope$dataset <- factor(meta_scope$dataset, levels = ds_scope)
  
  for (mi in seq_along(metrics)) {
    metric <- metrics[[mi]]
    d_scope <- subset_dist(merged_dist[[metric]], ids_scope)
    
    set.seed(45000 + si * 1000 + mi)
    ad <- vegan::adonis2(
      d_scope ~ dataset,
      data = meta_scope,
      permutations = N_PERM,
      strata = meta_scope$Barcode,
      by = "terms"
    )
    
    global_perm[[length(global_perm) + 1L]] <-
      extract_adonis_row(
        ad,
        term = "dataset",
        metric = metric,
        scope = scope_name,
        test_name = "stratified_PERMANOVA"
      )
  }
}

# Pairwise comparisons across all four datasets.
for (mi in seq_along(metrics)) {
  metric <- metrics[[mi]]
  
  for (pi in seq_along(dataset_pairs)) {
    z <- dataset_pairs[[pi]]
    keep <- meta$dataset %in% z
    ids <- rownames(meta)[keep]
    mp <- droplevels(meta[ids, , drop = FALSE])
    mp$dataset <- factor(mp$dataset, levels = z)
    dp <- subset_dist(merged_dist[[metric]], ids)
    
    set.seed(46000 + mi * 100 + pi)
    ad <- vegan::adonis2(
      dp ~ dataset,
      data = mp,
      permutations = N_PERM,
      strata = mp$Barcode,
      by = "terms"
    )
    
    rr <- extract_adonis_row(
      ad,
      term = "dataset",
      metric = metric,
      scope = "CD1_CD4",
      test_name = "pairwise_stratified_PERMANOVA"
    )
    rr$dataset1 <- z[1]
    rr$dataset2 <- z[2]
    pairwise_perm[[length(pairwise_perm) + 1L]] <- rr
  }
}

global_perm <- do.call(rbind, global_perm)
pairwise_perm <- do.call(rbind, pairwise_perm)

pairwise_perm$p_BH <- ave(
  pairwise_perm$p,
  pairwise_perm$metric,
  FUN = function(x) stats::p.adjust(x, method = "BH")
)

write_csv2(global_perm, file.path(OUT_DIR, "SINTAX_CD1-CD4_PERMANOVA_global_conf080.csv"))
write_csv2(pairwise_perm, file.path(OUT_DIR, "SINTAX_CD1-CD4_PERMANOVA_pairwise_conf080.csv"))

# 8. Dispersion: betadisper + permutest + paired Friedman.
dispersion_results <- list()
dispersion_distances <- list()

for (si in seq_along(scopes)) {
  scope_name <- names(scopes)[si]
  ds_scope <- scopes[[si]]
  
  ids_scope <- rownames(meta)[meta$dataset %in% ds_scope]
  ms <- droplevels(meta[ids_scope, , drop = FALSE])
  ms$dataset <- factor(ms$dataset, levels = ds_scope)
  
  for (mi in seq_along(metrics)) {
    metric <- metrics[[mi]]
    ds <- subset_dist(merged_dist[[metric]], ids_scope)
    
    bd <- vegan::betadisper(ds, ms$dataset)
    
    av <- stats::anova(bd)
    
    set.seed(47000 + si * 1000 + mi)
    pv <- vegan::permutest(bd, permutations = N_PERM)
    
    dd <- data.frame(
      scope = scope_name,
      metric = metric,
      obs_id = ids_scope,
      dataset = as.character(ms$dataset),
      Barcode = as.character(ms$Barcode),
      distance_to_centroid = bd$distances,
      stringsAsFactors = FALSE
    )
    
    # Ensure complete repeated-measures layout before Friedman.
    dwide <- reshape(
      dd[, c("Barcode", "dataset", "distance_to_centroid")],
      idvar = "Barcode",
      timevar = "dataset",
      direction = "wide"
    )
    
    dcols <- paste0("distance_to_centroid.", ds_scope)
    if (!all(dcols %in% names(dwide)) || anyNA(dwide[, dcols, drop = FALSE])) {
      stop("Incomplete dispersion repeated-measures layout for ", scope_name, " / ", metric)
    }
    
    fr <- stats::friedman.test(
      as.matrix(dwide[, dcols, drop = FALSE])
    )
    
    anova_p_col <- grep("^Pr\\(", names(av), value = TRUE)
    perm_tab <- pv$tab
    perm_p_col <- grep("^Pr\\(", names(perm_tab), value = TRUE)
    
    dispersion_results[[length(dispersion_results) + 1L]] <- data.frame(
      scope = scope_name,
      metric = metric,
      betadisper_F = as.numeric(av[1, "F value"]),
      betadisper_ANOVA_p = if (length(anova_p_col)) as.numeric(av[1, anova_p_col[1]]) else NA_real_,
      betadisper_permutest_F = if ("F" %in% names(perm_tab)) as.numeric(perm_tab[1, "F"]) else NA_real_,
      betadisper_permutest_p = if (length(perm_p_col)) as.numeric(perm_tab[1, perm_p_col[1]]) else NA_real_,
      Friedman_chisq = unname(fr$statistic),
      Friedman_df = unname(fr$parameter),
      Friedman_p = fr$p.value,
      Kendall_W = unname(fr$statistic) /
        (nrow(dwide) * (length(ds_scope) - 1L)),
      stringsAsFactors = FALSE
    )
    
    dispersion_distances[[length(dispersion_distances) + 1L]] <- dd
  }
}

dispersion_results <- do.call(rbind, dispersion_results)
dispersion_distances <- do.call(rbind, dispersion_distances)

write_csv2(
  dispersion_results,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_dispersion_tests_conf080.csv")
)
write_csv2(
  dispersion_distances,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_dispersion_distances_conf080.csv")
)

# 9. Barcode vs dataset variance partition using distance-based model.
#    Barcode is entered first, dataset second. Thus dataset R2 is the
#    methodological component remaining after accounting for barcode.
variance_partition <- list()

for (si in seq_along(scopes)) {
  scope_name <- names(scopes)[si]
  ds_scope <- scopes[[si]]
  
  ids_scope <- rownames(meta)[meta$dataset %in% ds_scope]
  ms <- droplevels(meta[ids_scope, , drop = FALSE])
  ms$dataset <- factor(ms$dataset, levels = ds_scope)
  
  for (mi in seq_along(metrics)) {
    metric <- metrics[[mi]]
    ds <- subset_dist(merged_dist[[metric]], ids_scope)
    
    # permutations = 0 because only R2 decomposition is needed here.
    vp <- vegan::adonis2(
      ds ~ Barcode + dataset,
      data = ms,
      permutations = 0,
      by = "terms"
    )
    
    rn <- rownames(vp)
    ib <- match("Barcode", rn)
    id <- match("dataset", rn)
    ir <- match("Residual", rn)
    
    if (anyNA(c(ib, id, ir))) {
      stop("Could not extract Barcode/dataset/Residual from variance model.")
    }
    
    r2_barcode <- as.numeric(vp[ib, "R2"])
    r2_dataset <- as.numeric(vp[id, "R2"])
    r2_resid <- as.numeric(vp[ir, "R2"])
    
    conditional_dataset <- r2_dataset / (r2_dataset + r2_resid)
    
    variance_partition[[length(variance_partition) + 1L]] <- data.frame(
      scope = scope_name,
      metric = metric,
      Barcode_variance_percent = 100 * r2_barcode,
      Dataset_variance_percent = 100 * r2_dataset,
      Residual_variance_percent = 100 * r2_resid,
      Dataset_percent_of_post_barcode_residual = 100 * conditional_dataset,
      sum_percent = 100 * (r2_barcode + r2_dataset + r2_resid),
      stringsAsFactors = FALSE
    )
  }
}

variance_partition <- do.call(rbind, variance_partition)

write_csv2(
  variance_partition,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_variance_partition_conf080.csv")
)

# 10. Compact master statistics table.
master <- rbind(
  data.frame(
    section = "PERMANOVA_global",
    global_perm,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  data.frame(
    section = "Dispersion",
    scope = dispersion_results$scope,
    metric = dispersion_results$metric,
    test = "Friedman_dispersion",
    term = "dataset",
    Df = dispersion_results$Friedman_df,
    SumOfSqs = NA_real_,
    R2 = NA_real_,
    F = dispersion_results$Friedman_chisq,
    p = dispersion_results$Friedman_p,
    stringsAsFactors = FALSE
  )
)

write_csv2(
  master,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_master_statistics_conf080.csv")
)

# 11. Plots.
#     ggplot2 only. Faceting replaces par()/mfrow graphics so plotting is

dataset_colours <- c(
  CD1 = "#1B9E77",
  CD2 = "#D95F02",
  CD3 = "#7570B3",
  CD4 = "#E7298A"
)

dataset_shapes <- c(
  CD1 = 16,
  CD2 = 17,
  CD3 = 15,
  CD4 = 18
)

theme_nmds <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey95"),
    strip.text = ggplot2::element_text(face = "bold")
  )

make_scope_nmds_df <- function(scope_name) {
  ds_scope <- scopes[[scope_name]]
  pieces <- vector("list", length(metrics))
  
  for (mi in seq_along(metrics)) {
    metric <- metrics[[mi]]
    fit <- scope_nmds[[scope_name]][[metric]]
    
    scr <- as.data.frame(
      vegan::scores(
        fit,
        display = "sites",
        choices = 1:2
      )
    )
    
    if (!all(c("NMDS1", "NMDS2") %in% names(scr))) {
      names(scr)[1:2] <- c("NMDS1", "NMDS2")
    }
    
    scr$obs_id <- rownames(scr)
    mp <- meta[scr$obs_id, c("dataset", "Barcode"), drop = FALSE]
    
    if (anyNA(mp$dataset) || anyNA(mp$Barcode)) {
      stop("Failed to map NMDS site scores back to metadata for ", scope_name, " / ", metric)
    }
    
    scr$dataset <- as.character(mp$dataset)
    scr$Barcode <- as.character(mp$Barcode)
    scr$metric <- metric
    scr$scope <- scope_name
    scr$stress <- fit$stress
    scr$facet_label <- paste0(
      metric,
      "\nstress = ",
      format_stress_value(fit$stress)
    )
    scr <- scr[scr$dataset %in% ds_scope, , drop = FALSE]
    pieces[[mi]] <- scr
  }
  
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out$dataset <- factor(out$dataset, levels = datasets)
  out$metric <- factor(out$metric, levels = metrics)
  out$scope <- factor(out$scope, levels = names(scopes))
  out
}

nmds_plot_data <- do.call(
  rbind,
  lapply(names(scopes), make_scope_nmds_df)
)
rownames(nmds_plot_data) <- NULL

write_csv2(
  nmds_plot_data,
  file.path(OUT_DIR, "SINTAX_CD1-CD4_NMDS_plot_coordinates_conf080.csv")
)

make_nmds_plot <- function(df, ds_scope, title_text = NULL) {
  df$dataset <- factor(as.character(df$dataset), levels = ds_scope)
  
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = NMDS1, y = NMDS2, group = Barcode)
  ) +
    ggplot2::geom_path(
      colour = "grey65",
      linewidth = 0.25,
      alpha = 0.45
    ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = dataset, shape = dataset),
      size = 1.55,
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~ facet_label,
      scales = "free",
      nrow = 1
    ) +
    ggplot2::scale_colour_manual(
      values = dataset_colours[ds_scope],
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = dataset_shapes[ds_scope],
      drop = FALSE
    ) +
    ggplot2::labs(
      x = "NMDS1",
      y = "NMDS2",
      colour = "Dataset",
      shape = "Dataset",
      title = title_text
    ) +
    theme_nmds
  
  p
}

# Separate 1 x 3 NMDS figures for each scope.
for (scope_name in names(scopes)) {
  ds_scope <- scopes[[scope_name]]
  df_scope <- nmds_plot_data[nmds_plot_data$scope == scope_name, , drop = FALSE]
  
  p_scope <- make_nmds_plot(
    df_scope,
    ds_scope = ds_scope,
    title_text = if (scope_name == "CD1_CD3") "A. CD1-CD3" else "B. CD1-CD4"
  )
  
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0("SINTAX_", scope_name, "_NMDS_conf080.pdf")),
    plot = p_scope,
    width = 12,
    height = 4.4,
    units = "in"
  )
  
  ggplot2::ggsave(
    filename = file.path(OUT_DIR, paste0("SINTAX_", scope_name, "_NMDS_conf080.png")),
    plot = p_scope,
    width = 12,
    height = 4.4,
    units = "in",
    dpi = 300
  )
}

# Combined manuscript-style 2 x 3 NMDS figure.
# IMPORTANT: use facet_wrap(), not facet_grid(), because facet_grid with
nmds_plot_data$scope_label <- ifelse(
  nmds_plot_data$scope == "CD1_CD3",
  "A. CD1-CD3",
  "B. CD1-CD4"
)

nmds_plot_data$panel_label <- paste0(
  nmds_plot_data$scope_label,
  "\n",
  as.character(nmds_plot_data$metric),
  "\nstress = ",
  vapply(
    nmds_plot_data$stress,
    format_stress_value,
    character(1)
  )
)

panel_levels <- unlist(
  lapply(
    c("CD1_CD3", "CD1_CD4"),
    function(sc) {
      vapply(
        metrics,
        function(mt) {
          ii <- which(
            nmds_plot_data$scope == sc &
              as.character(nmds_plot_data$metric) == mt
          )[1]
          
          nmds_plot_data$panel_label[ii]
        },
        character(1)
      )
    }
  )
)

nmds_plot_data$panel_label <- factor(
  nmds_plot_data$panel_label,
  levels = panel_levels
)

p_nmds_combined <- ggplot2::ggplot(
  nmds_plot_data,
  ggplot2::aes(
    x = NMDS1,
    y = NMDS2,
    group = Barcode
  )
) +
  ggplot2::geom_path(
    colour = "grey65",
    linewidth = 0.22,
    alpha = 0.40
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      colour = dataset,
      shape = dataset
    ),
    size = 1.25,
    alpha = 0.82
  ) +
  ggplot2::facet_wrap(
    ~ panel_label,
    ncol = 3,
    scales = "free"
  ) +
  ggplot2::scale_colour_manual(
    values = dataset_colours,
    drop = FALSE
  ) +
  ggplot2::scale_shape_manual(
    values = dataset_shapes,
    drop = FALSE
  ) +
  ggplot2::labs(
    x = "NMDS1",
    y = "NMDS2",
    colour = "Dataset",
    shape = "Dataset"
  ) +
  theme_nmds

ggplot2::ggsave(
  filename = file.path(
    OUT_DIR,
    "SINTAX_CD1-CD4_NMDS_6panel_conf080.pdf"
  ),
  plot = p_nmds_combined,
  width = 12,
  height = 8.2,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(
    OUT_DIR,
    "SINTAX_CD1-CD4_NMDS_6panel_conf080.png"
  ),
  plot = p_nmds_combined,
  width = 12,
  height = 8.2,
  units = "in",
  dpi = 300
)

# Paired alpha-richness plot for CD1, CD3, CD4.
rich_plot <- rich
rich_plot$dataset <- factor(as.character(rich_plot$dataset), levels = alpha_ds)

p_alpha <- ggplot2::ggplot(
  rich_plot,
  ggplot2::aes(x = dataset, y = richness, group = Barcode)
) +
  ggplot2::geom_line(
    colour = "grey70",
    linewidth = 0.25,
    alpha = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(colour = dataset),
    size = 1.2,
    alpha = 0.60,
    position = ggplot2::position_jitter(width = 0.05, height = 0)
  ) +
  ggplot2::geom_boxplot(
    ggplot2::aes(colour = dataset),
    width = 0.32,
    outlier.shape = NA,
    fill = NA,
    linewidth = 0.55
  ) +
  ggplot2::scale_colour_manual(values = dataset_colours[alpha_ds], drop = FALSE) +
  ggplot2::labs(
    x = "Community dataset",
    y = "SH richness",
    colour = "Dataset"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "none"
  )

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_richness_conf080.pdf"),
  plot = p_alpha,
  width = 6.5,
  height = 5.2,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "SINTAX_CD1-CD4_alpha_richness_conf080.png"),
  plot = p_alpha,
  width = 6.5,
  height = 5.2,
  units = "in",
  dpi = 300
)

# Three-set SH overlap plot using ggplot2 only.
A <- presence$CD1
B <- presence$CD3
C <- presence$CD4

venn_counts <- c(
  CD1_only = length(setdiff(A, union(B, C))),
  CD3_only = length(setdiff(B, union(A, C))),
  CD4_only = length(setdiff(C, union(A, B))),
  CD1_CD3_only = length(setdiff(intersect(A, B), C)),
  CD1_CD4_only = length(setdiff(intersect(A, C), B)),
  CD3_CD4_only = length(setdiff(intersect(B, C), A)),
  CD1_CD3_CD4 = length(Reduce(intersect, list(A, B, C)))
)

circle_df <- function(cx, cy, r, group, n = 360L) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta),
    group = group,
    stringsAsFactors = FALSE
  )
}

venn_circles <- rbind(
  circle_df(-0.65, 0.00, 1.15, "CD1"),
  circle_df( 0.65, 0.00, 1.15, "CD3"),
  circle_df( 0.00, 0.72, 1.15, "CD4")
)

venn_labels <- data.frame(
  x = c(-1.10, 1.10, 0.00, 0.00, -0.50, 0.50, 0.00),
  y = c(-0.05,-0.05, 1.30,-0.55,  0.55, 0.55, 0.18),
  label = as.character(venn_counts),
  stringsAsFactors = FALSE
)

dataset_labels <- data.frame(
  x = c(-1.25, 1.25, 0.00),
  y = c(-0.72,-0.72, 1.82),
  label = c("CD1", "CD3", "CD4"),
  dataset = c("CD1", "CD3", "CD4"),
  stringsAsFactors = FALSE
)

p_venn <- ggplot2::ggplot() +
  ggplot2::geom_polygon(
    data = venn_circles,
    ggplot2::aes(x = x, y = y, group = group, fill = group, colour = group),
    alpha = 0.10,
    linewidth = 0.7
  ) +
  ggplot2::geom_text(
    data = venn_labels,
    ggplot2::aes(x = x, y = y, label = label),
    size = 3.5
  ) +
  ggplot2::geom_text(
    data = dataset_labels,
    ggplot2::aes(x = x, y = y, label = label, colour = dataset),
    fontface = "bold",
    size = 4
  ) +
  ggplot2::scale_fill_manual(
    values = dataset_colours[c("CD1", "CD3", "CD4")],
    guide = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = dataset_colours[c("CD1", "CD3", "CD4")],
    guide = "none"
  ) +
  ggplot2::coord_fixed(xlim = c(-2, 2), ylim = c(-1.5, 2.0)) +
  ggplot2::labs(title = "Species Hypothesis overlap at SINTAX confidence >= 0.80") +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "SINTAX_CD1_CD3_CD4_SH_overlap_conf080.pdf"),
  plot = p_venn,
  width = 6.2,
  height = 5.8,
  units = "in"
)

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "SINTAX_CD1_CD3_CD4_SH_overlap_conf080.png"),
  plot = p_venn,
  width = 6.2,
  height = 5.8,
  units = "in",
  dpi = 300
)

write_csv2(
  data.frame(
    region = names(venn_counts),
    n_SH = as.numeric(venn_counts),
    stringsAsFactors = FALSE
  ),
  file.path(OUT_DIR, "SINTAX_CD1_CD3_CD4_SH_overlap_regions_conf080.csv")
)

# 11. Terminal summary.
cat("\n===== ANALYSIS COMPLETE =====\n")
cat("Output directory:\n  ", OUT_DIR, "\n\n", sep = "")

cat("SH overlap summary:\n")
print(overlap_summary, row.names = FALSE)

cat("\nAlpha global test:\n")
print(alpha_global, row.names = FALSE)

cat("\nAlpha pairwise tests:\n")
print(alpha_pairwise, row.names = FALSE)

cat("\nGlobal stratified PERMANOVA:\n")
print(
  global_perm[, c("scope", "metric", "F", "R2", "p")],
  row.names = FALSE
)

cat("\nPairwise stratified PERMANOVA:\n")
print(
  pairwise_perm[, c("metric", "dataset1", "dataset2", "F", "R2", "p", "p_BH")],
  row.names = FALSE
)

cat("\nDispersion tests:\n")
print(dispersion_results, row.names = FALSE)

cat("\nMantel / Procrustes:\n")
print(
  concordance[, c(
    "metric", "dataset1", "dataset2",
    "Mantel_spearman_r", "Mantel_p_BH",
    "Procrustes_correlation", "Procrustes_p_BH"
  )],
  row.names = FALSE
)

cat("\nVariance partition:\n")
print(variance_partition, row.names = FALSE)

cat("\nNMDS scope-level stress:\n")
print(nmds_scope_stress, row.names = FALSE)

cat("\nFiles written:\n")
for (f in sort(list.files(OUT_DIR, full.names = TRUE))) {
  cat("  ", f, "\n", sep = "")
}

# PART 3. CREATE A1-MATCHED NMDS FIGURES

# A1-matched NMDS plotting for the SINTAX CD1-CD4 analysis
# Styling is aligned to the A1 GitHub code and the existing UMI manuscript
#   - grey70 within-barcode segments
#   - A1 palette values

stopifnot(
  exists("scope_nmds"),
  exists("meta"),
  exists("scopes"),
  exists("metrics"),
  exists("OUT_DIR")
)

# 1. Exact dataset aesthetics
# Colours occur in the A1 GitHub palettes:

dataset_colours <- c(
  CD1 = "#7570B3",
  CD2 = "#56B4E9",
  CD3 = "#1B9E77",
  CD4 = "#D95F02"
)

dataset_shapes <- c(
  CD1 = 16,  # filled circle
  CD2 = 15,  # filled square
  CD3 = 18,  # filled diamond
  CD4 = 17   # filled triangle
)

# A1 ordination styling.
POINT_SIZE <- 1.0
POINT_ALPHA <- 0.66
SEGMENT_COLOUR <- "grey70"
SEGMENT_ALPHA <- 0.70
SEGMENT_LINEWIDTH <- 0.30
BASE_SIZE <- 8
STRESS_SIZE <- 3

cat("A1-matched dataset aesthetics:\n")
print(
  tibble::tibble(
    dataset = names(dataset_colours),
    colour = unname(dataset_colours),
    shape = unname(dataset_shapes),
    point_size = POINT_SIZE,
    point_alpha = POINT_ALPHA
  )
)
cat("\n")

# 2. Format metadata to the A1 plotting convention

meta_A1 <- meta %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::transmute(
    sample = sample,
    Barcode = as.character(Barcode),
    dataset = factor(
      as.character(dataset),
      levels = c("CD1", "CD2", "CD3", "CD4")
    )
  )

# 3. Build a plotting data frame from one scope-level NMDS object

format_nmds_A1 <- function(nmds_obj, scope_name, metric) {
  
  df <- vegan::scores(
    nmds_obj,
    display = "sites",
    choices = 1:2
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::left_join(meta_A1, by = "sample")
  
  if (!all(c("NMDS1", "NMDS2") %in% names(df))) {
    score_cols <- setdiff(
      names(df),
      c("sample", "Barcode", "dataset")
    )
    if (length(score_cols) < 2L) {
      stop("Could not identify two NMDS score columns.")
    }
    names(df)[match(score_cols[1:2], names(df))] <- c("NMDS1", "NMDS2")
  }
  
  ds_scope <- scopes[[scope_name]]
  
  df <- df %>%
    dplyr::filter(dataset %in% ds_scope) %>%
    dplyr::mutate(
      dataset = factor(as.character(dataset), levels = ds_scope),
      scope = scope_name,
      metric = metric
    ) %>%
    dplyr::arrange(Barcode, dataset)
  
  if (anyNA(df$Barcode) || anyNA(df$dataset)) {
    stop("Missing Barcode/dataset metadata after NMDS join.")
  }
  
  df
}

# 4. A1-style single NMDS panel

plot_nmds_A1 <- function(nmds_obj, scope_name, metric, show_legend = FALSE) {
  
  df <- format_nmds_A1(nmds_obj, scope_name, metric)
  ds_scope <- scopes[[scope_name]]
  
  stress_label <- paste0(
    "Stress = ",
    format_stress_value(nmds_obj$stress)
  )
  
  # Same-barcode paths reproduce the manuscript visual comparison among
  # community datasets while retaining the A1 grey-segment convention.
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = df,
      ggplot2::aes(
        x = NMDS1,
        y = NMDS2,
        group = Barcode
      ),
      colour = SEGMENT_COLOUR,
      alpha = SEGMENT_ALPHA,
      linewidth = SEGMENT_LINEWIDTH
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(
        x = NMDS1,
        y = NMDS2,
        colour = dataset,
        shape = dataset
      ),
      size = POINT_SIZE,
      alpha = POINT_ALPHA
    ) +
    ggplot2::scale_colour_manual(
      values = dataset_colours[ds_scope],
      breaks = ds_scope,
      drop = FALSE,
      name = "Dataset"
    ) +
    ggplot2::scale_shape_manual(
      values = dataset_shapes[ds_scope],
      breaks = ds_scope,
      drop = FALSE,
      name = "Dataset"
    ) +
    ggplot2::annotate(
      "text",
      x = min(df$NMDS1, na.rm = TRUE),
      y = max(df$NMDS2, na.rm = TRUE),
      label = stress_label,
      hjust = 0,
      vjust = 1,
      size = STRESS_SIZE,
      fontface = "italic"
    ) +
    ggplot2::labs(
      x = "NMDS1",
      y = "NMDS2",
      title = metric
    ) +
    ggplot2::theme_minimal(base_size = BASE_SIZE) +
    ggplot2::theme(
      legend.position = if (show_legend) "bottom" else "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      ),
      legend.title = ggplot2::element_blank()
    )
  
  p
}

# 5. Create all six panels

plot_dir <- file.path(OUT_DIR, "A1_matched_NMDS")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

panel_list <- list()

for (metric in metrics) {
  for (scope_name in names(scopes)) {
    
    key <- paste(metric, scope_name, sep = "__")
    
    panel_list[[key]] <- plot_nmds_A1(
      nmds_obj = scope_nmds[[scope_name]][[metric]],
      scope_name = scope_name,
      metric = metric,
      show_legend = FALSE
    )
    
    # Print each panel to the active graphics device.
    print(panel_list[[key]])
    
    file_stub <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      paste(scope_name, metric, sep = "__")
    )
    
    # These dimensions match the A1 quarter-A4 ordination panels.
    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(file_stub, "_A1matched.png")
      ),
      plot = panel_list[[key]],
      width = 68.5,
      height = 62.5,
      units = "mm",
      dpi = 900,
      bg = "white"
    )
    
    ggplot2::ggsave(
      filename = file.path(
        plot_dir,
        paste0(file_stub, "_A1matched.pdf")
      ),
      plot = panel_list[[key]],
      width = 68.5,
      height = 62.5,
      units = "mm",
      bg = "white"
    )
  }
}

# 6. Shared legend

legend_source <- plot_nmds_A1(
  nmds_obj = scope_nmds[["CD1_CD4"]][["Bray-Curtis"]],
  scope_name = "CD1_CD4",
  metric = "Bray-Curtis",
  show_legend = TRUE
) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1,
      override.aes = list(
        size = POINT_SIZE,
        alpha = 1
      )
    ),
    shape = ggplot2::guide_legend(
      nrow = 1,
      override.aes = list(
        size = POINT_SIZE,
        alpha = 1
      )
    )
  )

shared_legend <- cowplot::get_legend(legend_source)

# 7. Combined manuscript layout: 3 rows x 2 columns
# Existing manuscript figure:

ordered_panels <- list(
  panel_list[["Bray-Curtis__CD1_CD3"]] +
    ggplot2::labs(subtitle = "A. CD1-CD3"),
  panel_list[["Bray-Curtis__CD1_CD4"]] +
    ggplot2::labs(subtitle = "B. CD1-CD4"),
  
  panel_list[["Hellinger__CD1_CD3"]],
  panel_list[["Hellinger__CD1_CD4"]],
  
  panel_list[["Robust-Aitchison__CD1_CD3"]],
  panel_list[["Robust-Aitchison__CD1_CD4"]]
)

combined_body <- cowplot::plot_grid(
  plotlist = ordered_panels,
  ncol = 2,
  align = "hv",
  axis = "tblr"
)

combined <- cowplot::plot_grid(
  combined_body,
  shared_legend,
  ncol = 1,
  rel_heights = c(1, 0.055)
)

print(combined)

ggplot2::ggsave(
  filename = file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_3x2.png"
  ),
  plot = combined,
  width = 145,
  height = 195,
  units = "mm",
  dpi = 900,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_3x2.pdf"
  ),
  plot = combined,
  width = 145,
  height = 195,
  units = "mm",
  bg = "white"
)

# 8. Save exact plot data + style provenance

plot_data_all <- purrr::map_dfr(
  metrics,
  function(metric) {
    purrr::map_dfr(
      names(scopes),
      function(scope_name) {
        format_nmds_A1(
          scope_nmds[[scope_name]][[metric]],
          scope_name,
          metric
        )
      }
    )
  }
)

readr::write_csv2(
  plot_data_all,
  file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_coordinates.csv"
  )
)

style_table <- tibble::tibble(
  dataset = names(dataset_colours),
  colour = unname(dataset_colours),
  shape = unname(dataset_shapes),
  point_size = POINT_SIZE,
  point_alpha = POINT_ALPHA,
  segment_colour = SEGMENT_COLOUR,
  segment_alpha = SEGMENT_ALPHA,
  segment_linewidth = SEGMENT_LINEWIDTH,
  theme = paste0("theme_minimal(base_size=", BASE_SIZE, ")")
)

readr::write_csv2(
  style_table,
  file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_style.csv"
  )
)

saveRDS(
  list(
    plots = panel_list,
    combined = combined,
    plot_data = plot_data_all,
    style = style_table
  ),
  file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_plot_objects.rds"
  )
)

cat("\nA1 NMDS stress labels used in figures:\n")
stress_label_table <- do.call(
  rbind,
  lapply(names(scopes), function(scope_name) {
    do.call(
      rbind,
      lapply(metrics, function(metric) {
        x <- scope_nmds[[scope_name]][[metric]]$stress
        data.frame(
          scope = scope_name,
          metric = metric,
          k = unname(NMDS_K_BY_METRIC[[metric]]),
          stress_exact = x,
          stress_plot_label = format_stress_value(x),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
print(stress_label_table, row.names = FALSE)

readr::write_csv2(
  stress_label_table,
  file.path(
    plot_dir,
    "SINTAX_CD1-CD4_NMDS_A1matched_stress_labels.csv"
  )
)

cat("\n===== A1-MATCHED NMDS FIGURES COMPLETE =====\n")
cat("Output directory:\n  ", plot_dir, "\n\n", sep = "")
cat("Dataset mapping:\n")
print(style_table)
cat("\nCombined figure:\n  ",
    file.path(plot_dir, "SINTAX_CD1-CD4_NMDS_A1matched_3x2.png"),
    "\n", sep = "")

cat("\n")
cat("=====================================================================\n")
cat("ALL-IN-ONE WORKFLOW COMPLETED SUCCESSFULLY\n")
cat("=====================================================================\n")
cat("Filtered community table:\n  ", OUTPUT, "\n", sep = "")
cat("Analysis directory:\n  ", OUT_DIR, "\n", sep = "")
cat("A1-matched NMDS directory:\n  ", plot_dir, "\n", sep = "")
