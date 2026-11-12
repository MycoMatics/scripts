#!/usr/bin/env Rscript
# Compare original SH matching and current SINTAX on the same read pairs.
# Restrict to the A1-filtered 90-barcode cohort; write eligibility audit.


required_packages <- c("data.table")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
})

root <- Sys.getenv("PROJECT_ROOT")
if (!nzchar(root)) stop("Set PROJECT_ROOT to the UMIvsRAW project root.")
root <- normalizePath(root, mustWork = TRUE)

pair_dir <- file.path(
  root,
  "REVIEW",
  "SINTAX_equal_ITS_v1",
  "paired_by_barcode"
)

original_file <- file.path(
  root,
  "REVIEW",
  "annotation_015_rebuilt",
  "annotation_015_complete.tsv"
)

out_dir <- file.path(
  root,
  "REVIEW",
  "SINTAX_equal_ITS_v1",
  "comparison"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. Completed SINTAX barcodes

pair_files <- list.files(
  pair_dir,
  pattern = "\\.paired\\.tsv\\.gz$",
  full.names = TRUE
)

if (!length(pair_files)) {
  stop("No completed paired SINTAX files found.")
}

processed_barcodes <- sub(
  "\\.paired\\.tsv\\.gz$",
  "",
  basename(pair_files)
)

# Use the same A1-filtered 90 biological barcodes as the community analysis.
cohort_file <- file.path(root, 'REVIEW', 'SINTAX_equal_ITS_v1',
                        'analysis', 'appendixC_current', 'cohorts', 'A1filtered90.txt')
if (!file.exists(cohort_file)) stop('Run scripts/06_classification_success.py first: missing ', cohort_file)
cohort <- sort(readLines(cohort_file, warn = FALSE))
if (length(cohort) != 90L || anyDuplicated(cohort)) {
  stop('Expected 90 distinct A1-filtered barcodes')
}
if (!setequal(cohort, processed_barcodes[processed_barcodes %in% cohort])) {
  stop('Some A1-filtered barcodes have no current SINTAX paired output')
}
pair_files <- pair_files[processed_barcodes %in% cohort]
processed_barcodes <- processed_barcodes[processed_barcodes %in% cohort]
out_dir <- file.path(root, 'REVIEW', 'SINTAX_equal_ITS_v1',
                     'analysis', 'appendixC_current', 'contingency_A1filtered90')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(
  "Completed barcodes = ",
  length(processed_barcodes),
  "\n",
  sep = ""
)

# 2. Read minimal SINTAX outcome table

cat("Reading SINTAX paired outcomes...\n")

sintax <- rbindlist(
  lapply(
    pair_files,
    function(f) {

      fread(
        cmd = paste(
          "gzip -cd",
          shQuote(f)
        ),
        select = c(
          "Barcode",
          "fastq_ID",
          "annotation_080"
        ),
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE
)

cat(
  "SINTAX D3 rows = ",
  format(
    nrow(sintax),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

# 3. Validate unique read IDs

n_dup <- anyDuplicated(
  sintax[
    ,
    .(
      Barcode,
      fastq_ID
    )
  ]
)

if (n_dup > 0) {

  stop(
    paste0(
      "Duplicated Barcode + fastq_ID in SINTAX table. ",
      "First duplicated position: ",
      n_dup
    )
  )
}

# 4. Harmonize SINTAX outcome labels

sintax[
  ,
  SINTAX_outcome :=
    fcase(
      annotation_080 == "matched",
      "Matched",

      annotation_080 == "unmatched",
      "Unmatched",

      annotation_080 == "mismatched",
      "Mismatched",

      annotation_080 ==
        "excluded_D1_low_confidence",
      "D1 low confidence",

      default = NA_character_
    )
]

if (anyNA(sintax$SINTAX_outcome)) {

  print(
    sintax[
      is.na(SINTAX_outcome),
      .N,
      by = annotation_080
    ]
  )

  stop(
    "Unexpected SINTAX annotation value."
  )
}

# 5. Read original 1.5% SH annotation table
# annotation_015_complete.tsv contains only reads whose
# paired D1 consensus was eligible for the original

cat("Reading original SH annotation table...\n")

original <- fread(
  original_file,
  select = c(
    "Barcode",
    "fastq_ID",
    "annot_015"
  ),
  showProgress = TRUE
)

original <- original[
  Barcode %chin%
    processed_barcodes
]

cat(
  "Original SH-comparable reads in completed barcodes = ",
  format(
    nrow(original),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

# 6. Validate original IDs

n_dup_original <- anyDuplicated(
  original[
    ,
    .(
      Barcode,
      fastq_ID
    )
  ]
)

if (n_dup_original > 0) {

  stop(
    paste0(
      "Duplicated Barcode + fastq_ID in original table. ",
      "First duplicated position: ",
      n_dup_original
    )
  )
}

# 7. Harmonize original SH outcome labels

original[
  ,
  SH_outcome :=
    fcase(
      annot_015 == "matched",
      "Matched",

      annot_015 == "unmatched",
      "Unmatched",

      annot_015 == "mismatched",
      "Mismatched",

      default = NA_character_
    )
]

if (anyNA(original$SH_outcome)) {

  print(
    original[
      is.na(SH_outcome),
      .N,
      by = annot_015
    ]
  )

  stop(
    "Unexpected original annot_015 value."
  )
}

# 8. Exact read-level join

cat("Joining original SH and SINTAX outcomes by Barcode + fastq_ID...\n")

setkey(
  sintax,
  Barcode,
  fastq_ID
)

cmp <- sintax[
  original[
    ,
    .(
      Barcode,
      fastq_ID,
      SH_outcome
    )
  ],
  on = .(
    Barcode,
    fastq_ID
  )
]

if (nrow(cmp) != nrow(original)) {

  stop(
    paste0(
      "Unexpected join row count. Original = ",
      nrow(original),
      "; joined = ",
      nrow(cmp)
    )
  )
}

n_missing_sintax <- sum(
  is.na(
    cmp$SINTAX_outcome
  )
)

if (n_missing_sintax > 0) {

  stop(
    paste0(
      "Missing SINTAX outcome for ",
      n_missing_sintax,
      " original-comparison reads."
    )
  )
}

# 9. Extended contingency
# Original SH comparison is the row denominator.
#   D1 low confidence
# smaller the SINTAX-eligible denominator is.

extended_long <- cmp[
  ,
  .N,
  by = .(
    SH_outcome,
    SINTAX_outcome
  )
]

extended_long[
  ,
  pct_within_SH_outcome :=
    100 *
    N /
    sum(N),
  by = SH_outcome
]

extended_long[
  ,
  pct_of_all_original_comparable :=
    100 *
    N /
    nrow(cmp)
]

setorder(
  extended_long,
  SH_outcome,
  SINTAX_outcome
)

extended_counts <- dcast(
  extended_long,
  SH_outcome ~ SINTAX_outcome,
  value.var = "N",
  fill = 0
)

fwrite(
  extended_long,
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_extended_long.tsv"
  ),
  sep = "\t"
)

fwrite(
  extended_counts,
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_extended_counts.tsv"
  ),
  sep = "\t"
)

# 10. Shared eligible denominator
#   must have SH confidence >=0.8.

shared <- cmp[
  SINTAX_outcome !=
    "D1 low confidence"
]

shared_long <- shared[
  ,
  .N,
  by = .(
    SH_outcome,
    SINTAX_outcome
  )
]

shared_long[
  ,
  pct_within_SH_outcome :=
    100 *
    N /
    sum(N),
  by = SH_outcome
]

shared_long[
  ,
  pct_of_shared_denominator :=
    100 *
    N /
    nrow(shared)
]

setorder(
  shared_long,
  SH_outcome,
  SINTAX_outcome
)

shared_counts <- dcast(
  shared_long,
  SH_outcome ~ SINTAX_outcome,
  value.var = "N",
  fill = 0
)

fwrite(
  shared_long,
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_shared_long.tsv"
  ),
  sep = "\t"
)

fwrite(
  shared_counts,
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_shared_counts.tsv"
  ),
  sep = "\t"
)

# 11. Agreement statistics

n_sintax_excluded <- cmp[
  SINTAX_outcome ==
    "D1 low confidence",
  .N
]

same_outcome <- shared[
  SH_outcome ==
    SINTAX_outcome,
  .N
]

different_outcome <- (
  nrow(shared) -
    same_outcome
)

# 12. Category-specific transitions

transition_summary <- copy(
  shared_long
)

transition_summary[
  ,
  SH_total :=
    sum(N),
  by = SH_outcome
]

transition_summary[
  ,
  pct_of_SH_category :=
    100 *
    N /
    SH_total
]

fwrite(
  transition_summary,
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_transitions.tsv"
  ),
  sep = "\t"
)

# 13. Console report

cat("\n")
cat("============================================================\n")
cat("ORIGINAL SH vs SINTAX 0.8\n")
cat("============================================================\n")

cat(
  "Completed barcodes                         = ",
  length(processed_barcodes),
  "\n",
  sep = ""
)

cat(
  "Original SH-comparable D3 reads            = ",
  format(
    nrow(cmp),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "SINTAX D1 below confidence 0.8             = ",
  format(
    n_sintax_excluded,
    big.mark = ","
  ),
  " (",
  sprintf(
    "%.2f",
    100 *
      n_sintax_excluded /
      nrow(cmp)
  ),
  "% of original SH-comparable reads)\n",
  sep = ""
)

cat(
  "Shared eligible denominator                = ",
  format(
    nrow(shared),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "Same outcome category on shared denominator = ",
  format(
    same_outcome,
    big.mark = ","
  ),
  " (",
  sprintf(
    "%.2f",
    100 *
      same_outcome /
      nrow(shared)
  ),
  "%)\n",
  sep = ""
)

cat(
  "Different outcome category                 = ",
  format(
    different_outcome,
    big.mark = ","
  ),
  " (",
  sprintf(
    "%.2f",
    100 *
      different_outcome /
      nrow(shared)
  ),
  "%)\n",
  sep = ""
)

cat("\n")
cat("Extended contingency table:\n")
cat("Includes SINTAX D1-low-confidence reads.\n\n")

print(
  extended_counts
)

cat("\n")
cat("Shared-denominator 3 x 3 contingency table:\n\n")

print(
  shared_counts
)

cat("\n")
cat("Outputs:\n")

cat(
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_extended_counts.tsv"
  ),
  "\n"
)

cat(
  file.path(
    out_dir,
    "contingency_SH_vs_SINTAX080_shared_counts.tsv"
  ),
  "\n"
)
