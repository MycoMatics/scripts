#!/usr/bin/env Rscript
# Figure C3: SH-matching clustering thresholds versus SINTAX confidence.
# Plot proportions using method-specific input sequence denominators.

# Figure C3. Classification success across SH and SINTAX thresholds.
# Base R ONLY: no tidyverse, ggplot2, cowplot or other contributed packages.
# Figure uses the A1 palette (D1 purple, D3 green, D4 orange), white background,
#     Dataset, Confidence, Total_sequences, Classified (tab-separated)

write_csv2 <- function(x, file) utils::write.csv2(x, file, row.names = FALSE)

review_dir <- Sys.getenv("REVIEW_DIR")
if (!nzchar(review_dir)) stop("Set REVIEW_DIR to the UMIvsRAW/REVIEW directory.")
review_dir <- normalizePath(review_dir, mustWork = TRUE)
out_dir <- file.path(review_dir, "SINTAX_equal_ITS_v1", "analysis", "appendixC_current", "figure_C3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sh_file <- Sys.getenv(
  "SH_COUNTS_FILE", file.path(out_dir, "SH_classifiability_A1filtered90.csv")
)
sintax_file <- Sys.getenv(
  "SINTAX_COUNTS_FILE", file.path(review_dir, "SINTAX_equal_ITS_v1", "analysis", "appendixC_current", "Table_C2_SINTAX_A1filtered90_pooled.tsv")
)
cohort_file <- file.path(review_dir, "SINTAX_equal_ITS_v1", "analysis", "appendixC_current", "cohorts", "A1filtered90.txt")

sh_thresholds <- c(0.5, 1, 1.5, 2, 2.5, 3)
sintax_thresholds <- c(0.70, 0.80, 0.90)
datasets <- c("D1", "D3", "D4")
colours <- c(D1 = "#7570B3", D3 = "#1B9E77", D4 = "#D95F02")

if (!file.exists(sh_file)) {
  stop("Missing matched-cohort original SH-matching counts: ", sh_file,
       "\nRun scripts/08_extract_sh_matching_rates.py first.")
}
if (!file.exists(sintax_file)) {
  stop("Validated SINTAX Table C2 summary not found: ", sintax_file)
}

if (grepl("\\.tsv$", sh_file, ignore.case = TRUE)) {
  sh <- utils::read.delim(sh_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  sh <- utils::read.csv2(sh_file, check.names = FALSE, stringsAsFactors = FALSE)
}
sintax <- utils::read.delim(
  sintax_file, check.names = FALSE, stringsAsFactors = FALSE
)

required_sh <- c("Dataset", "SH_threshold_pct", "Total_sequences", "Classified")
required_si <- c("Dataset", "Confidence", "Total_sequences", "Classified")
if (!all(required_sh %in% names(sh))) {
  stop("SH input requires: ", paste(required_sh, collapse = ", "))
}
if (!all(required_si %in% names(sintax))) {
  stop("SINTAX input requires: ", paste(required_si, collapse = ", "))
}

sh$Dataset <- as.character(sh$Dataset)
sh$threshold <- as.numeric(sh$SH_threshold_pct)
sh$Total_sequences <- as.numeric(sh$Total_sequences)
sh$Classified <- as.numeric(sh$Classified)
sintax$Dataset <- as.character(sintax$Dataset)
sintax$threshold <- as.numeric(sintax$Confidence)
sintax$Total_sequences <- as.numeric(sintax$Total_sequences)
sintax$Classified <- as.numeric(sintax$Classified)

if ("Barcode" %in% names(sh)) {
  if (!file.exists(cohort_file)) {
    stop("Per-barcode SH data requires exact cohort list: ", cohort_file)
  }
  cohort <- trimws(readLines(cohort_file, warn = FALSE))
  cohort <- cohort[nzchar(cohort)]
  if (length(cohort) != 90L || anyDuplicated(cohort)) {
    stop("The current A1-filtered cohort must contain exactly 90 distinct barcodes.")
  }
  if (anyDuplicated(sh[, c("Barcode", "Dataset", "threshold")])) {
    stop("Duplicate Barcode/Dataset/threshold combination in SH input.")
  }
  for (ds in datasets) {
    for (t in sh_thresholds) {
      ix <- sh$Dataset == ds & abs(sh$threshold - t) < 1e-8
      if (!identical(sort(as.character(sh$Barcode[ix])), sort(cohort))) {
        stop("SH barcode set differs from current A1-filtered cohort: ", ds, ", threshold ", t)
      }
    }
    per_bc <- split(sh$Total_sequences[sh$Dataset == ds],
                    as.character(sh$Barcode[sh$Dataset == ds]))
    if (any(vapply(per_bc, function(x) length(unique(x)) != 1L, logical(1)))) {
      stop("SH input denominator changes across thresholds within barcodes: ", ds)
    }
  }
  sh <- stats::aggregate(
    cbind(Total_sequences, Classified) ~ Dataset + threshold,
    data = sh,
    FUN = sum
  )
  sh$Barcodes <- 90L
} else {
  if (!"Barcodes" %in% names(sh)) {
    stop("Pooled SH input needs a 'Barcodes' column (90 in every row).")
  }
  if (anyNA(sh$Barcodes) || any(as.numeric(sh$Barcodes) != 90)) {
    stop("Every pooled SH row must represent exactly 90 barcodes.")
  }
  cat("NOTE: SH values are pooled. Exact IDs cannot be audited without per-barcode input.\n")
}

validate_grid <- function(x, thresholds, method) {
  observed <- paste(x$Dataset, sprintf("%.2f", x$threshold), sep = "|")
  expected <- expand.grid(
    Dataset = datasets,
    threshold = thresholds,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  expected_keys <- paste(expected$Dataset,
                         sprintf("%.2f", expected$threshold), sep = "|")
  if (length(observed) != length(expected_keys) ||
      anyDuplicated(observed) || !setequal(observed, expected_keys)) {
    stop(method, " must contain one row per dataset × threshold (",
         length(expected_keys), " unique combinations).")
  }
}
validate_grid(sh, sh_thresholds, "SH matching")
validate_grid(sintax, sintax_thresholds, "SINTAX")

for (x in list(sh, sintax)) {
  if (anyNA(x[, c("threshold", "Total_sequences", "Classified")]) ||
      any(!is.finite(x$Total_sequences)) ||
      any(!is.finite(x$Classified)) ||
      any(x$Total_sequences <= 0) ||
      any(x$Classified < 0) ||
      any(x$Classified > x$Total_sequences) ||
      any(x$Total_sequences != floor(x$Total_sequences)) ||
      any(x$Classified != floor(x$Classified))) {
    stop("Invalid, missing or non-integer total/classified sequence counts.")
  }
  for (ds in datasets) {
    if (length(unique(x$Total_sequences[x$Dataset == ds])) != 1L) {
      stop("Sequence denominator changes across thresholds within ", ds)
    }
  }
}

sh_total <- vapply(datasets, function(ds) {
  sh$Total_sequences[which(sh$Dataset == ds)[1L]]
}, numeric(1))
si_total <- vapply(datasets, function(ds) {
  sintax$Total_sequences[which(sintax$Dataset == ds)[1L]]
}, numeric(1))
comparison <- data.frame(
  Dataset = datasets,
  SH_total = unname(sh_total),
  SINTAX_total = unname(si_total)
)
cat("\nSequence-denominator validation:\n")
print(comparison, row.names = FALSE)
if (any(sh_total != si_total)) {
  warning(
    "SH and SINTAX sequence denominators differ. The two panels report ",
    "percentages relative to method-specific processed sequence totals, NOT ",
    "necessarily identical underlying reads. See printed audit."
  )
  write_csv2(comparison, file.path(out_dir,
    "Figure_C3_SH_vs_SINTAX_denominator_audit.csv"))
}

sh$Method <- "SH matching"
sintax$Method <- "SINTAX"
sh$classified_pct <- 100 * sh$Classified / sh$Total_sequences
sintax$classified_pct <- 100 * sintax$Classified / sintax$Total_sequences

final_table <- rbind(
  sh[, c("Method", "Dataset", "threshold", "Total_sequences", "Classified", "classified_pct")],
  sintax[, c("Method", "Dataset", "threshold", "Total_sequences", "Classified", "classified_pct")]
)
final_table <- final_table[order(match(final_table$Method, c("SH matching", "SINTAX")),
                                 match(final_table$Dataset, datasets),
                                 final_table$threshold), , drop = FALSE]
rownames(final_table) <- NULL

out_stub <- file.path(out_dir, "Figure_C3_SH_vs_SINTAX_classification_thresholds")
write_csv2(final_table, paste0(out_stub, "_data.csv"))

# Device-independent base R plotting. No package installation required.
draw_figure <- function() {
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::layout(matrix(c(1L, 2L, 3L), ncol = 1L), heights = c(1, 1, 0.14))

  draw_panel <- function(x, thresholds, panel, xlabel, fmt) {
    graphics::par(mar = c(3.7, 4.5, 2.25, 0.9),
                  mgp = c(2.3, 0.65, 0), xaxs = "i", yaxs = "i",
                  family = "sans", bg = "white")
    values <- matrix(NA_real_, nrow = length(datasets),
                     ncol = length(thresholds),
                     dimnames = list(datasets, NULL))
    for (i in seq_along(datasets)) {
      for (j in seq_along(thresholds)) {
        ix <- x$Dataset == datasets[i] & abs(x$threshold - thresholds[j]) < 1e-8
        values[i, j] <- x$classified_pct[ix]
      }
    }

    # barplot(plot=FALSE) supplies the grouped-bar positions, then grid lines
    # are drawn behind all bars for a clean manuscript-style appearance.
    centers <- graphics::barplot(values, beside = TRUE, plot = FALSE,
                                 width = 0.8, space = c(0.10, 0.85))
    graphics::plot.new()
    graphics::plot.window(xlim = range(centers) + c(-0.9, 0.9),
                          ylim = c(0, 100))
    graphics::abline(h = seq(0, 100, by = 20), col = "#E7E7E7", lwd = 0.8)
    graphics::rect(xleft = as.vector(centers) - 0.40,
                   ybottom = 0,
                   xright = as.vector(centers) + 0.40,
                   ytop = as.vector(values),
                   col = rep(unname(colours[datasets]), times = ncol(values)),
                   border = NA)
    graphics::axis(2, at = seq(0, 100, 20), las = 1, tick = FALSE,
                   cex.axis = 0.78, col.axis = "#222222")
    graphics::axis(1, at = colMeans(centers), labels = fmt(thresholds),
                   tick = FALSE, cex.axis = 0.87, col.axis = "#222222")
    graphics::mtext(panel, side = 3, line = 0.55, adj = 0,
                    font = 2, cex = 0.98)
    graphics::mtext("Assigned sequences (% of method input)", side = 2, line = 2.65,
                    cex = 0.84)
    graphics::mtext(xlabel, side = 1, line = 2.15, cex = 0.84)
  }

  draw_panel(sh, sh_thresholds, "A. SH matching",
             "UNITE SH clustering threshold (%)", function(t) sprintf("%.1f", t))
  draw_panel(sintax, sintax_thresholds, "B. SINTAX",
             "Bootstrap confidence threshold", function(t) sprintf("%.2f", t))

  graphics::par(mar = c(0, 0, 0, 0), bg = "white")
  graphics::plot.new()
  graphics::legend("center", legend = datasets, fill = unname(colours[datasets]),
                   horiz = TRUE, bty = "n", cex = 0.95,
                   x.intersp = 0.6, y.intersp = 0.9)
}

grDevices::pdf(paste0(out_stub, ".pdf"), width = 145 / 25.4,
               height = 170 / 25.4, bg = "white", useDingbats = FALSE)
draw_figure()
grDevices::dev.off()

png_args <- list(filename = paste0(out_stub, ".png"),
                 width = 145, height = 170, units = "mm", res = 900,
                 bg = "white")
if (isTRUE(capabilities("cairo"))) png_args$type <- "cairo"
do.call(grDevices::png, png_args)
draw_figure()
grDevices::dev.off()

cat("\n===== FIGURE C3 INPUT SUMMARY =====\n")
print(final_table, row.names = FALSE, digits = 5)
cat("\nSaved:\n  ", out_stub, ".pdf\n  ", out_stub, ".png\n  ",
    out_stub, "_data.csv\n", sep = "")
cat("SH clustering thresholds and SINTAX confidence cutoffs are not calibrated equivalents.\n")
