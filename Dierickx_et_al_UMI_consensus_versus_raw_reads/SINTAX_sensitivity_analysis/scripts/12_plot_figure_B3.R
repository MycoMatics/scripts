#!/usr/bin/env Rscript
# Figure B3: Q score versus raw-consensus divergence across annotation outcomes.
# Verify all 8,432,242 read pairs against the 90-barcode contingency first.


required_packages <- c("data.table", "ggplot2", "MASS", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(MASS)
  library(scales)
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

identity_file <- file.path(
  root,
  "REVIEW",
  "raw_vs_umi_identity_harmonized",
  "raw_vs_umi_consensus_identity.tsv"
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
  "analysis",
  "appendixC_current",
  "q_divergence_figures_A1filtered90"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

SEED <- 20260909

SAMPLE_MAX <- 100000L

Q_BIN_WIDTH <- 0.25

DIV_BIN_WIDTH <- 0.10

Q_DISPLAY <- c(
  10,
  40
)

DIV_DISPLAY <- c(
  0,
  6
)

HDR_PROBS <- c(
  0.50,
  0.95
)

KDE_N_Q <- 260L

KDE_N_DIV <- 220L

annotation_levels <- c(
  "Matched",
  "Unmatched",
  "Mismatched"
)

# Cohort and cross-script validation for revised Figure B3

A1_COHORT_FILE <- file.path(
  root, "REVIEW", "SINTAX_equal_ITS_v1", "analysis",
  "appendixC_current", "cohorts", "A1filtered90.txt"
)

CONTINGENCY_FILE <- file.path(
  root, "REVIEW", "SINTAX_equal_ITS_v1", "analysis",
  "appendixC_current", "contingency_A1filtered90",
  "contingency_SH_vs_SINTAX080_shared_long.tsv"
)

EXPECTED_SHARED_N <- 8432242L

if (!file.exists(A1_COHORT_FILE)) {
  stop("Missing the validated A1 cohort file: ", A1_COHORT_FILE)
}
if (!file.exists(CONTINGENCY_FILE)) {
  stop("Missing the 90-barcode SH vs SINTAX contingency: ", CONTINGENCY_FILE)
}

A1_barcodes <- trimws(readLines(A1_COHORT_FILE, warn = FALSE))
A1_barcodes <- A1_barcodes[nzchar(A1_barcodes)]
if (length(A1_barcodes) != 90L || anyDuplicated(A1_barcodes)) {
  stop("Expected exactly 90 unique A1-filtered barcode IDs")
}

cat("Figure B3 cohort: 90 A1-filtered barcodes\n")
cat("Expected shared eligible denominator: ", EXPECTED_SHARED_N, "\n", sep = "")

# 1. Completed barcodes and SINTAX outcomes

pair_files <- list.files(
  pair_dir,
  pattern = "\\.paired\\.tsv\\.gz$",
  full.names = TRUE
)

if (!length(pair_files)) {
  stop("No completed SINTAX paired files.")
}

processed_barcodes <- sub(
  "\\.paired\\.tsv\\.gz$",
  "",
  basename(pair_files)
)

# Figure B3 must use precisely the same 90 biological barcodes as the
# Appendix C cross-algorithm contingency, not all 97 completed barcodes.
if (!all(A1_barcodes %in% processed_barcodes)) {
  stop("Missing paired SINTAX results for: ",
       paste(setdiff(A1_barcodes, processed_barcodes), collapse = ", "))
}
pair_files <- pair_files[processed_barcodes %in% A1_barcodes]
processed_barcodes <- processed_barcodes[processed_barcodes %in% A1_barcodes]
if (length(pair_files) != 90L ||
    !identical(sort(processed_barcodes), sort(A1_barcodes))) {
  stop("Figure B3 90-barcode paired-file selection failed")
}

cat(
  "Completed barcodes = ",
  length(processed_barcodes),
  "\n",
  sep = ""
)

cat("Reading SINTAX annotations...\n")

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
  stop("Unexpected SINTAX annotation.")
}

# 2. Original SH outcomes for the same completed barcodes

cat("Reading original 1.5% outcomes...\n")

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
  stop("Unexpected original SH outcome.")
}

# 3. Exact outcome join

setkey(
  sintax,
  Barcode,
  fastq_ID
)

comparison <- sintax[
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

if (anyNA(comparison$SINTAX_outcome)) {
  stop("Missing SINTAX outcomes after join.")
}

shared <- comparison[
  SINTAX_outcome !=
    "D1 low confidence"
]

# Compare against the existing Appendix C 90-barcode contingency. Both
# analyses must contain the same read-level pairs and outcome categories.
if (nrow(shared) != EXPECTED_SHARED_N) {
  stop("Figure B3 shared denominator differs from Appendix C: figure = ",
       nrow(shared), "; expected = ", EXPECTED_SHARED_N)
}

observed_cells <- shared[, .(N = .N), by = .(SH_outcome, SINTAX_outcome)]
reference_cells <- data.table::fread(CONTINGENCY_FILE,
                                     select = c("SH_outcome", "SINTAX_outcome", "N"),
                                     showProgress = FALSE)
compare_cells <- merge(observed_cells, reference_cells,
                       by = c("SH_outcome", "SINTAX_outcome"),
                       all = TRUE, suffixes = c("_figure", "_appendix"))
if (anyNA(compare_cells) ||
    any(compare_cells$N_figure != compare_cells$N_appendix)) {
  print(compare_cells)
  stop("Figure B3 outcome counts differ from Appendix C contingency")
}

cat("Figure B3 exact 90-barcode contingency check: PASS\n")
print(compare_cells[order(compare_cells$SH_outcome,
                          compare_cells$SINTAX_outcome), ])

cohort_audit <- data.table::data.table(
  measure = c("barcodes", "shared_eligible_reads",
              "SH_matched", "SH_unmatched", "SH_mismatched",
              "SINTAX_matched", "SINTAX_unmatched", "SINTAX_mismatched"),
  N = c(length(processed_barcodes), nrow(shared),
        sum(shared$SH_outcome == "Matched"),
        sum(shared$SH_outcome == "Unmatched"),
        sum(shared$SH_outcome == "Mismatched"),
        sum(shared$SINTAX_outcome == "Matched"),
        sum(shared$SINTAX_outcome == "Unmatched"),
        sum(shared$SINTAX_outcome == "Mismatched"))
)
data.table::fwrite(cohort_audit,
  file.path(out_dir, "Figure_B3_A1filtered90_denominator_audit.tsv"), sep = "\t")

# 4. Identity / fidelity table

cat("Reading harmonized Q/divergence data...\n")

identity <- fread(
  identity_file,
  select = c(
    "Barcode",
    "fastq_ID",
    "processed_avg_qscore",
    "overlap_error_pct"
  ),
  showProgress = TRUE
)

setnames(
  identity,
  c(
    "processed_avg_qscore",
    "overlap_error_pct"
  ),
  c(
    "qscore",
    "divergence"
  )
)

identity <- identity[
  is.finite(qscore) &
    is.finite(divergence)
]

setkey(
  identity,
  Barcode,
  fastq_ID
)

# 5. Common broad estimation range
# Derived once from the full harmonized identity data.

q_raw <- quantile(
  identity$qscore,
  c(0.001, 0.999),
  na.rm = TRUE
)

Q_EST <- c(
  floor(
    q_raw[1] /
      Q_BIN_WIDTH
  ) *
    Q_BIN_WIDTH,

  ceiling(
    q_raw[2] /
      Q_BIN_WIDTH
  ) *
    Q_BIN_WIDTH
)

div_raw <- quantile(
  identity$divergence,
  0.995,
  na.rm = TRUE
)

DIV_EST <- c(
  0,

  ceiling(
    div_raw /
      DIV_BIN_WIDTH
  ) *
    DIV_BIN_WIDTH
)

cat(
  "Q estimation range = ",
  paste(Q_EST, collapse = " to "),
  "\n"
)

cat(
  "Divergence estimation range = ",
  paste(DIV_EST, collapse = " to "),
  "\n"
)

# 6. Q-predicted mean base-error reference

expected_curve <- data.table(
  qscore = seq(
    Q_EST[1],
    Q_EST[2],
    length.out = 600
  )
)

expected_curve[
  ,
  predicted_error :=
    100 *
    10^(
      -qscore / 10
    )
]

# 7. Aggregate one algorithm/outcome table

set.seed(SEED)

aggregate_algorithm <- function(
  annotations,
  algorithm_name
) {

  annotations <- copy(
    annotations
  )

  annotations[
    ,
    annotation :=
      factor(
        annotation,
        levels = annotation_levels
      )
  ]

  dd <- identity[
    annotations,
    on = .(
      Barcode,
      fastq_ID
    ),
    nomatch = 0
  ]

  if (nrow(dd) != nrow(annotations)) {

    stop(
      paste(
        algorithm_name,
        "identity join mismatch:",
        nrow(annotations),
        "annotations versus",
        nrow(dd),
        "joined"
      )
    )
  }

  # Read-level descriptive statistics

  stats <- dd[
    ,
    .(
      N = .N,

      q_mean =
        mean(qscore),

      q_median =
        median(qscore),

      q_q25 =
        quantile(
          qscore,
          0.25
        ),

      q_q75 =
        quantile(
          qscore,
          0.75
        ),

      div_mean =
        mean(divergence),

      div_median =
        median(divergence),

      div_q25 =
        quantile(
          divergence,
          0.25
        ),

      div_q75 =
        quantile(
          divergence,
          0.75
        )
    ),
    by = annotation
  ]

  stats[
    ,
    percentage :=
      100 *
      N /
      sum(N)
  ]

  stats[
    ,
    algorithm :=
      algorithm_name
  ]

  # Barcode-level medians

  barcode_stats <- dd[
    ,
    .(
      N = .N,
      q_median =
        median(qscore),
      div_median =
        median(divergence)
    ),
    by = .(
      Barcode,
      annotation
    )
  ]

  matched_bc <- barcode_stats[
    annotation == "Matched" &
      N >= 50,
    .(
      Barcode,
      matched_N = N,
      matched_q =
        q_median,
      matched_div =
        div_median
    )
  ]

  unmatched_bc <- barcode_stats[
    annotation == "Unmatched" &
      N >= 50,
    .(
      Barcode,
      unmatched_N = N,
      unmatched_q =
        q_median,
      unmatched_div =
        div_median
    )
  ]

  effects <- merge(
    matched_bc,
    unmatched_bc,
    by = "Barcode"
  )

  effects[
    ,
    `:=`(
      q_delta_unmatched_minus_matched =
        unmatched_q -
        matched_q,

      div_delta_unmatched_minus_matched =
        unmatched_div -
        matched_div,

      algorithm =
        algorithm_name
    )
  ]

  # Empirical 2D bins from ALL reads

  bins <- dd[
    qscore >= Q_EST[1] &
      qscore < Q_EST[2] &
      divergence >= DIV_EST[1] &
      divergence < DIV_EST[2],
    .N,
    by = .(
      annotation,

      q_bin =
        floor(
          (
            qscore -
              Q_EST[1]
          ) /
            Q_BIN_WIDTH
        ),

      div_bin =
        floor(
          (
            divergence -
              DIV_EST[1]
          ) /
            DIV_BIN_WIDTH
        )
    )
  ]

  bins[
    ,
    q_mid :=
      Q_EST[1] +
      (
        q_bin +
          0.5
      ) *
      Q_BIN_WIDTH
  ]

  bins[
    ,
    div_mid :=
      DIV_EST[1] +
      (
        div_bin +
          0.5
      ) *
      DIV_BIN_WIDTH
  ]

  bins[
    ,
    algorithm :=
      algorithm_name
  ]

  # Balanced sample for KDE estimation

  est <- dd[
    qscore >= Q_EST[1] &
      qscore <= Q_EST[2] &
      divergence >= DIV_EST[1] &
      divergence <= DIV_EST[2],
    {

      take_n <- min(
        .N,
        SAMPLE_MAX
      )

      idx <- if (.N > take_n) {
        sample.int(
          .N,
          take_n
        )
      } else {
        seq_len(.N)
      }

      .(
        x =
          qscore[idx],

        y =
          divergence[idx]
      )

    },
    by = annotation
  ]

  est[
    ,
    algorithm :=
      algorithm_name
  ]

  rm(dd)
  gc()

  list(
    stats = stats,
    barcode_effects = effects,
    bins = bins,
    est = est
  )
}

# 8. Natural SINTAX 0.8 denominator

sintax_natural_ann <- sintax[
  SINTAX_outcome !=
    "D1 low confidence",
  .(
    Barcode,
    fastq_ID,
    annotation =
      SINTAX_outcome
  )
]

cat(
  "SINTAX natural eligible reads = ",
  format(
    nrow(sintax_natural_ann),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

sintax_natural <- aggregate_algorithm(
  sintax_natural_ann,
  "SINTAX 0.8"
)

# 9. Shared denominator for direct algorithm comparison

shared_SH_ann <- shared[
  ,
  .(
    Barcode,
    fastq_ID,
    annotation =
      SH_outcome
  )
]

shared_SINTAX_ann <- shared[
  ,
  .(
    Barcode,
    fastq_ID,
    annotation =
      SINTAX_outcome
  )
]

cat(
  "Shared eligible reads = ",
  format(
    nrow(shared),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

shared_SH <- aggregate_algorithm(
  shared_SH_ann,
  "Original SH matching"
)

shared_SINTAX <- aggregate_algorithm(
  shared_SINTAX_ann,
  "SINTAX 0.8"
)

# 10. Reflected KDE and HDR extraction

calculate_hdr <- function(
  est_data
) {

  h_common <- c(
    MASS::bandwidth.nrd(
      est_data$x
    ),

    MASS::bandwidth.nrd(
      est_data$y
    )
  )

  cat(
    "Common KDE bandwidth: ",
    h_common[1],
    ", ",
    h_common[2],
    "\n",
    sep = ""
  )

  panels <- unique(
    est_data[
      ,
      .(
        algorithm,
        annotation
      )
    ]
  )

  output <- list()

  out_id <- 1L

  for (ii in seq_len(nrow(panels))) {

    alg <- panels$algorithm[ii]

    ann <- panels$annotation[ii]

    sub <- est_data[
      algorithm == alg &
        annotation == ann
    ]

    # Boundary reflection around divergence = 0

    x_aug <- c(
      sub$x,
      sub$x
    )

    y_aug <- c(
      sub$y,
      -sub$y
    )

    kde <- MASS::kde2d(
      x = x_aug,
      y = y_aug,
      h = h_common,
      n = c(
        KDE_N_Q,
        2L *
          KDE_N_DIV -
          1L
      ),
      lims = c(
        Q_EST[1],
        Q_EST[2],
        -DIV_EST[2],
        DIV_EST[2]
      )
    )

    keep_y <- kde$y >= 0

    yy <- kde$y[
      keep_y
    ]

    zz <- kde$z[
      ,
      keep_y,
      drop = FALSE
    ] * 2

    zvals <- as.vector(
      zz
    )

    zsorted <- sort(
      zvals,
      decreasing = TRUE
    )

    mass <- cumsum(
      zsorted
    ) /
      sum(zsorted)

    for (pp in HDR_PROBS) {

      idx <- which(
        mass >= pp
      )[1]

      threshold <- zsorted[
        idx
      ]

      contours <- contourLines(
        x = kde$x,
        y = yy,
        z = zz,
        levels = threshold
      )

      for (jj in seq_along(contours)) {

        output[[out_id]] <- data.table(
          algorithm = alg,
          annotation = ann,

          hdr_level =
            ifelse(
              pp == 0.50,
              "50% HDR",
              "95% HDR"
            ),

          segment =
            paste(
              alg,
              ann,
              pp,
              jj,
              sep = "_"
            ),

          x =
            contours[[jj]]$x,

          y =
            contours[[jj]]$y
        )

        out_id <- out_id + 1L
      }
    }
  }

  rbindlist(output)
}

# 11. Generic plotting function

make_density_figure <- function(
  bins,
  est,
  stats,
  facet_mode = c(
    "single",
    "comparison"
  ),
  caption_text
) {

  facet_mode <- match.arg(
    facet_mode
  )

  bins[
    ,
    annotation :=
      factor(
        annotation,
        levels = annotation_levels
      )
  ]

  est[
    ,
    annotation :=
      factor(
        annotation,
        levels = annotation_levels
      )
  ]

  stats[
    ,
    annotation :=
      factor(
        annotation,
        levels = annotation_levels
      )
  ]

  hdr <- calculate_hdr(
    est
  )

  hdr[
    ,
    annotation :=
      factor(
        annotation,
        levels = annotation_levels
      )
  ]

  hdr[
    ,
    hdr_level :=
      factor(
        hdr_level,
        levels = c(
          "50% HDR",
          "95% HDR"
        )
      )
  ]

  labels <- stats[
    ,
    .(
      algorithm,
      annotation,

      label =
        paste0(
          ifelse(
            percentage >= 1,
            sprintf(
              "%.1f%%",
              percentage
            ),
            sprintf(
              "%.3f%%",
              percentage
            )
          ),
          "\nn = ",
          comma(N)
        )
    )
  ]

  global_max <- max(
    bins$N
  )

  p <- ggplot() +

    geom_tile(
      data = bins,

      aes(
        x = q_mid,
        y = div_mid,
        fill = N
      ),

      width =
        Q_BIN_WIDTH,

      height =
        DIV_BIN_WIDTH,

      alpha =
        0.95
    ) +

    geom_path(
      data = hdr,

      aes(
        x = x,
        y = y,
        group = segment,
        linetype = hdr_level
      ),

      colour =
        "black",

      linewidth =
        0.70
    ) +

    geom_line(
      data =
        expected_curve,

      aes(
        x =
          qscore,

        y =
          predicted_error
      ),

      colour =
        "gray40",

      linewidth =
        0.60
    ) +

    geom_label(
      data =
        labels,

      aes(
        x =
          Q_DISPLAY[1] +
          0.5,

        y =
          DIV_DISPLAY[2] -
          0.25,

        label =
          label
      ),

      hjust =
        0,

      vjust =
        1,

      size =
        2.8,

      fill =
        "white",

      label.size =
        0.20
    ) +

    scale_fill_viridis_c(
      option =
        "D",

      trans =
        "log10",

      limits =
        c(
          1,
          global_max
        ),

      oob =
        squish,

      name =
        "Reads per bin",

      labels =
        comma
    ) +

    scale_linetype_manual(
      values =
        c(
          "50% HDR" =
            "dotted",

          "95% HDR" =
            "22"
        ),

      name =
        "Highest-density region"
    ) +

    coord_cartesian(
      xlim =
        Q_DISPLAY,

      ylim =
        DIV_DISPLAY,

      expand =
        FALSE
    ) +

    labs(
      x =
        "Per-read Q-score",

      y =
        "Raw-consensus divergence (%)",

      caption =
        caption_text
    ) +

    theme_classic(
      base_size =
        11
    ) +

    theme(
      strip.text =
        element_text(
          face =
            "bold",
          size =
            10.5
        ),

      legend.position =
        "bottom",

      plot.caption =
        element_text(
          hjust =
            0,
          size =
            8
        )
    )

  if (
    facet_mode ==
      "single"
  ) {

    p <- p +
      facet_wrap(
        ~ annotation,
        nrow = 1
      )

  } else {

    p <- p +
      facet_grid(
        algorithm ~ annotation
      )

  }

  p
}

# 12. SINTAX-only 3-panel figure

caption_sintax <- paste0(
  "SINTAX outcomes were defined at SH-level confidence >=0.8 ",
  "for 90 A1-filtered barcodes. ",
  "Only raw reads whose paired UMI consensus received an SH ",
  "assignment with confidence >=0.8 are included. ",
  "Background intensity gives absolute read counts in fixed ",
  Q_BIN_WIDTH,
  " Q-score x ",
  DIV_BIN_WIDTH,
  "% divergence bins. Black contours indicate the 50% dotted ",
  "and 95% dashed highest-density regions. The solid grey curve ",
  "shows mean base-error probability predicted from the FASTQ ",
  "Q-score. Raw-consensus divergence was calculated over the ",
  "shared homologous interval."
)

p_sintax <- make_density_figure(
  bins =
    sintax_natural$bins,

  est =
    sintax_natural$est,

  stats =
    sintax_natural$stats,

  facet_mode =
    "single",

  caption_text =
    caption_sintax
)

ggsave(
  file.path(
    out_dir,
    "SINTAX080_q_divergence_3panel.png"
  ),

  p_sintax,

  width =
    250,

  height =
    130,

  units =
    "mm",

  dpi =
    600
)

ggsave(
  file.path(
    out_dir,
    "SINTAX080_q_divergence_3panel.pdf"
  ),

  p_sintax,

  width =
    250,

  height =
    130,

  units =
    "mm"
)

# 13. Direct algorithm comparison on shared denominator

comparison_bins <- rbindlist(
  list(
    shared_SH$bins,
    shared_SINTAX$bins
  ),
  use.names = TRUE
)

comparison_est <- rbindlist(
  list(
    shared_SH$est,
    shared_SINTAX$est
  ),
  use.names = TRUE
)

comparison_stats <- rbindlist(
  list(
    shared_SH$stats,
    shared_SINTAX$stats
  ),
  use.names = TRUE
)

comparison_bins[
  ,
  algorithm :=
    factor(
      algorithm,
      levels = c(
        "Original SH matching",
        "SINTAX 0.8"
      )
    )
]

comparison_est[
  ,
  algorithm :=
    factor(
      algorithm,
      levels = c(
        "Original SH matching",
        "SINTAX 0.8"
      )
    )
]

comparison_stats[
  ,
  algorithm :=
    factor(
      algorithm,
      levels = c(
        "Original SH matching",
        "SINTAX 0.8"
      )
    )
]

caption_comparison <- paste0(
  "Original SH matching and current ITS-standardized SINTAX were ",
  "compared in 90 A1-filtered barcodes using a common denominator ",
  "of ", format(nrow(shared), big.mark = ","), " read pairs. ",
  "These were reads whose paired UMI consensus was ",
  "eligible under the original 1.5% SH procedure and received a ",
  "consensus-level SINTAX SH assignment with confidence >=0.8. ",
  "Finite read Q-scores and raw-consensus divergence values were required. ",
  "Background ",
  "intensity gives absolute read counts using a common log-scaled ",
  "colour range across all six panels. Black contours indicate ",
  "the 50% dotted and 95% dashed highest-density regions. The ",
  "solid grey curve gives Q-predicted mean base-error probability."
)

p_compare <- make_density_figure(
  bins =
    comparison_bins,

  est =
    comparison_est,

  stats =
    comparison_stats,

  facet_mode =
    "comparison",

  caption_text =
    caption_comparison
)

ggsave(
  file.path(
    out_dir,
    "SH_vs_SINTAX080_q_divergence_6panel.png"
  ),

  p_compare,

  width =
    250,

  height =
    190,

  units =
    "mm",

  dpi =
    600
)

ggsave(
  file.path(
    out_dir,
    "SH_vs_SINTAX080_q_divergence_6panel.pdf"
  ),

  p_compare,

  width =
    250,

  height =
    190,

  units =
    "mm"
)

# 14. Save statistics

fwrite(
  sintax_natural$stats,

  file.path(
    out_dir,
    "SINTAX080_q_divergence_stats.tsv"
  ),

  sep =
    "\t"
)

fwrite(
  sintax_natural$barcode_effects,

  file.path(
    out_dir,
    "SINTAX080_barcode_fidelity_effects.tsv"
  ),

  sep =
    "\t"
)

fwrite(
  comparison_stats,

  file.path(
    out_dir,
    "SH_vs_SINTAX080_shared_q_divergence_stats.tsv"
  ),

  sep =
    "\t"
)

comparison_effects <- rbindlist(
  list(
    shared_SH$barcode_effects,
    shared_SINTAX$barcode_effects
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  comparison_effects,

  file.path(
    out_dir,
    "SH_vs_SINTAX080_shared_barcode_fidelity_effects.tsv"
  ),

  sep =
    "\t"
)

# 15. Barcode-level robustness report

cat("\n============================================================\n")
cat("SINTAX 0.8 BARCODE-LEVEL FIDELITY EFFECT\n")
cat("============================================================\n")

eff <- sintax_natural$barcode_effects

cat(
  "Barcodes with >=50 matched and >=50 unmatched reads = ",
  nrow(eff),
  "\n",
  sep = ""
)

if (nrow(eff)) {

  cat(
    "Unmatched lower median Q in ",
    sum(
      eff$q_delta_unmatched_minus_matched <
        0
    ),
    "/",
    nrow(eff),
    " barcodes\n",
    sep = ""
  )

  cat(
    "Unmatched greater median divergence in ",
    sum(
      eff$div_delta_unmatched_minus_matched >
        0
    ),
    "/",
    nrow(eff),
    " barcodes\n",
    sep = ""
  )

  cat(
    "Median Q delta = ",
    sprintf(
      "%.3f",
      median(
        eff$q_delta_unmatched_minus_matched
      )
    ),
    "\n",
    sep = ""
  )

  cat(
    "Median divergence delta = ",
    sprintf(
      "%.3f",
      median(
        eff$div_delta_unmatched_minus_matched
      )
    ),
    " percentage points\n",
    sep = ""
  )
}

cat("\nOutputs written to:\n")
cat(out_dir, "\n")
