# =========================================================
# HARMONISED NMDS VISUALISATION WORKFLOW
#
# Produces four PNG figures only:
#   1. 2-panel all-forest NMDS: Full dataset vs Subsetted
#   2. 10-panel site NMDS: All forests + single-forest NMDS
#   3. 10-panel species NMDS: Full dataset, species coloured by guild
#   4. 10-panel species NMDS: Subsetted dataset, species coloured by guild
#
# Additional output printed to terminal:
#   - NMDS stress summaries
#   - Forest centroid shifts
#   - Full-dataset weighted centroid shifts
#   - Duplicated abbreviated taxon labels, if present
#
# Required input objects:
#   X_use_full, meta_use_full
#   X_use_subset, meta_use_subset
#
# Optional aliases also accepted:
#   X_full, meta_full
#   X_subset, meta_subset
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(mgcv)
  library(grid)
})

set.seed(42)

# =========================================================
# 0. SETTINGS
# =========================================================

outdir_nmds <- "nmds_outputs_harmonised"
dir.create(outdir_nmds, showWarnings = FALSE, recursive = TRUE)

forest_levels <- c("Sonian", "Rajhenav", "Strødam", "Suserup")
timespan_levels <- c("previous", "current")
dataset_levels <- c("Full dataset", "Subsetted")

row_panel_levels <- c(
  "All forests",
  "Sonian",
  "Rajhenav",
  "Strødam",
  "Suserup"
)

timespan_labels <- c(
  previous = "Historical",
  current  = "Current"
)

survey_period_levels <- unname(timespan_labels)

forest_cols <- c(
  Sonian   = "#B8860B",
  Rajhenav = "#2E7D32",
  Strødam  = "#466C95",
  Suserup  = "#5FA8D3"
)

timespan_shape_vals <- c(
  previous = 21,
  current  = 24
)

decay_breaks <- seq(1.5, 4.5, by = 1)

guild_levels <- c(
  "Early ruderals",
  "Combative invader",
  "Trunk rotters",
  "Late stage specialist",
  "Cord",
  "Unknown/different"
)

guild_labeller <- c(
  "Early ruderals"        = "Early ruderals",
  "Combative invader"     = "Combative invaders",
  "Trunk rotters"         = "Trunk rotters",
  "Late stage specialist" = "Late stage specialists",
  "Cord"                  = "Cord-formers",
  "Unknown/different"     = "Unknown/different"
)

guild_cols <- c(
  "Early ruderals"        = "#D55E00",
  "Combative invader"     = "#CC79A7",
  "Trunk rotters"         = "#0072B2",
  "Late stage specialist" = "#009E73",
  "Cord"                  = "#E69F00",
  "Unknown/different"     = "grey65"
)

species_min_occ <- 5
label_top_n_per_panel <- 10
show_site_cloud <- TRUE

use_species_label_callouts <- TRUE
species_label_box_padding <- 0.10
species_label_point_padding <- 0.08

a4_width_in       <- 7.70
a4_height_2panel  <- 4.25
a4_height_10panel <- 10.80

use_contour_text <- requireNamespace("metR", quietly = TRUE)

if (!use_contour_text) {
  warning("Package `metR` is not installed. Contour labels will be skipped.")
}

lighten_col <- function(col, amount = 0.45) {
  colorRampPalette(c(col, "white"))(101)[round(amount * 100) + 1]
}

forest_time_cols <- c(
  "Sonian.previous"   = lighten_col(forest_cols[["Sonian"]], 0.45),
  "Sonian.current"    = forest_cols[["Sonian"]],
  "Rajhenav.previous" = lighten_col(forest_cols[["Rajhenav"]], 0.45),
  "Rajhenav.current"  = forest_cols[["Rajhenav"]],
  "Strødam.previous"  = lighten_col(forest_cols[["Strødam"]], 0.45),
  "Strødam.current"   = forest_cols[["Strødam"]],
  "Suserup.previous"  = lighten_col(forest_cols[["Suserup"]], 0.45),
  "Suserup.current"   = forest_cols[["Suserup"]]
)

forest_time_labels <- c(
  "Sonian.previous"   = "Sonian · Historical",
  "Sonian.current"    = "Sonian · Current",
  "Rajhenav.previous" = "Rajhenav · Historical",
  "Rajhenav.current"  = "Rajhenav · Current",
  "Strødam.previous"  = "Strødam · Historical",
  "Strødam.current"   = "Strødam · Current",
  "Suserup.previous"  = "Suserup · Historical",
  "Suserup.current"   = "Suserup · Current"
)

save_png <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(outdir_nmds, filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 900,
    bg = "transparent"
  )
}

# =========================================================
# 1. LABEL FUNCTIONS
# =========================================================

abbreviate_taxon_label <- function(x) {
  
  x <- stringr::str_squish(as.character(x))
  
  ifelse(
    stringr::str_detect(x, "^[A-Za-z][A-Za-z-]+\\s+"),
    stringr::str_replace(x, "^([A-Za-z])[A-Za-z-]+\\s+", "\\1. "),
    x
  )
}

resolve_abbreviated_label_overlaps <- function(label_df) {
  
  if (nrow(label_df) == 0) {
    return(label_df %>% mutate(Taxon_label = character(0)))
  }
  
  label_df <- label_df %>%
    mutate(
      Taxon_label_abbrev = abbreviate_taxon_label(Taxon)
    )
  
  duplicated_labels <- label_df %>%
    count(row_panel, survey_period, nmds_scope, Taxon_label_abbrev, name = "n_label_duplicates") %>%
    filter(n_label_duplicates > 1)
  
  if (nrow(duplicated_labels) > 0) {
    
    cat("\n=========================================================\n")
    cat("DUPLICATED ABBREVIATED TAXON LABELS DETECTED\n")
    cat("Full names will be used for these labels within the affected panels.\n")
    cat("=========================================================\n")
    
    print(
      label_df %>%
        semi_join(
          duplicated_labels,
          by = c("row_panel", "survey_period", "nmds_scope", "Taxon_label_abbrev")
        ) %>%
        arrange(row_panel, survey_period, Taxon_label_abbrev, Taxon) %>%
        select(row_panel, survey_period, nmds_scope, Taxon_label_abbrev, Taxon),
      n = Inf
    )
  }
  
  label_df %>%
    left_join(
      duplicated_labels %>%
        mutate(use_full_label = TRUE) %>%
        select(row_panel, survey_period, nmds_scope, Taxon_label_abbrev, use_full_label),
      by = c("row_panel", "survey_period", "nmds_scope", "Taxon_label_abbrev")
    ) %>%
    mutate(
      use_full_label = if_else(is.na(use_full_label), FALSE, use_full_label),
      Taxon_label = if_else(use_full_label, Taxon, Taxon_label_abbrev)
    )
}

# =========================================================
# 2. INPUT PREPARATION
# =========================================================

if (!exists("X_full")) {
  X_full <- X_use_full
  meta_full <- meta_use_full
}

if (!exists("X_subset")) {
  X_subset <- X_use_subset
  meta_subset <- meta_use_subset
}

standardise_nmds_inputs <- function(X, meta, dataset_name) {
  
  X <- as.data.frame(X)
  meta <- as.data.frame(meta)
  
  X[] <- lapply(X, function(z) as.integer(as.numeric(z) > 0))
  
  if (!"sample_id" %in% names(meta)) {
    meta$sample_id <- rownames(X)
  }
  
  meta <- meta %>%
    mutate(
      sample_id = as.character(sample_id),
      forest = str_replace_all(as.character(forest), "Strodam", "Strødam"),
      forest = factor(forest, levels = forest_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(
        timespan_labels[as.character(timespan)],
        levels = survey_period_levels
      ),
      log_av_ds = suppressWarnings(as.numeric(log_av_ds)),
      data_scope = factor(dataset_name, levels = dataset_levels)
    )
  
  rownames(X) <- meta$sample_id
  
  keep_samples <- rowSums(X, na.rm = TRUE) > 0
  X <- X[keep_samples, , drop = FALSE]
  meta <- meta[match(rownames(X), meta$sample_id), , drop = FALSE]
  
  keep_taxa <- colSums(X, na.rm = TRUE) > 0
  X <- X[, keep_taxa, drop = FALSE]
  
  stopifnot(
    identical(rownames(X), meta$sample_id),
    nrow(X) == nrow(meta),
    !any(names(meta) %in% names(X))
  )
  
  list(X = X, meta = meta)
}

full_prepared <- standardise_nmds_inputs(
  X = X_full,
  meta = meta_full,
  dataset_name = "Full dataset"
)

subset_prepared <- standardise_nmds_inputs(
  X = X_subset,
  meta = meta_subset,
  dataset_name = "Subsetted"
)

X_full <- full_prepared$X
meta_full <- full_prepared$meta

X_subset <- subset_prepared$X
meta_subset <- subset_prepared$meta

cat("\n=========================================================\n")
cat("INPUT CHECKS\n")
cat("=========================================================\n")

cat("\nFull dataset forest × period:\n")
print(table(meta_full$forest, meta_full$timespan, useNA = "ifany"))

cat("\nSubsetted dataset forest × period:\n")
print(table(meta_subset$forest, meta_subset$timespan, useNA = "ifany"))

cat("\nFull dataset decay-stage distribution:\n")
print(table(meta_full$log_av_ds, useNA = "ifany"))

cat("\nSubsetted dataset decay-stage distribution:\n")
print(table(meta_subset$log_av_ds, useNA = "ifany"))

# =========================================================
# 3. NMDS FUNCTIONS
# =========================================================

run_global_nmds <- function(X, meta, dataset_name, trymax = 20) {
  
  cat("\n=========================================================\n")
  cat("GLOBAL NMDS:", toupper(dataset_name), "\n")
  cat("=========================================================\n")
  cat("Samples:", nrow(X), "\n")
  cat("Taxa:", ncol(X), "\n")
  
  nmds <- metaMDS(
    X,
    distance = "jaccard",
    binary = TRUE,
    k = 2,
    trymax = trymax,
    autotransform = FALSE,
    noshare = TRUE,
    trace = TRUE
  )
  
  scores_df <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    rename(
      NMDS1 = 2,
      NMDS2 = 3
    ) %>%
    left_join(meta, by = "sample_id") %>%
    mutate(
      data_scope = factor(dataset_name, levels = dataset_levels),
      forest = factor(as.character(forest), levels = forest_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      forest_timespan = factor(
        paste(forest, timespan, sep = "."),
        levels = names(forest_time_cols)
      ),
      nmds_scope = paste(dataset_name, "All forests", sep = " | "),
      row_panel = factor("All forests", levels = row_panel_levels),
      panel_type = "Global NMDS",
      log_av_ds = suppressWarnings(as.numeric(log_av_ds))
    )
  
  centroids_df <- scores_df %>%
    group_by(data_scope, row_panel, forest, timespan, forest_timespan, nmds_scope) %>%
    summarise(
      NMDS1 = mean(NMDS1, na.rm = TRUE),
      NMDS2 = mean(NMDS2, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  shifts_df <- centroids_df %>%
    select(data_scope, row_panel, forest, nmds_scope, timespan, NMDS1, NMDS2, n) %>%
    pivot_wider(
      names_from = timespan,
      values_from = c(NMDS1, NMDS2, n)
    ) %>%
    filter(
      !is.na(NMDS1_previous),
      !is.na(NMDS1_current),
      !is.na(NMDS2_previous),
      !is.na(NMDS2_current)
    ) %>%
    mutate(
      d_NMDS1 = NMDS1_current - NMDS1_previous,
      d_NMDS2 = NMDS2_current - NMDS2_previous,
      centroid_shift_distance_2D = sqrt(d_NMDS1^2 + d_NMDS2^2),
      stress = nmds$stress
    )
  
  stress_df <- tibble(
    data_scope = factor(dataset_name, levels = dataset_levels),
    row_panel = factor("All forests", levels = row_panel_levels),
    nmds_scope = paste(dataset_name, "All forests", sep = " | "),
    panel_type = "Global NMDS",
    stress = nmds$stress,
    n_samples = nrow(X),
    n_taxa = ncol(X)
  )
  
  cat("\nStress:", nmds$stress, "\n")
  cat("\nForest centroid shifts:\n")
  print(shifts_df, n = Inf)
  
  list(
    nmds = nmds,
    X = X,
    meta = meta,
    scores = scores_df,
    centroids = centroids_df,
    shifts = shifts_df,
    stress = stress_df
  )
}

orient_nmds1_by_decay <- function(scores_df) {
  
  dat_cor <- scores_df %>%
    filter(
      !is.na(NMDS1),
      !is.na(log_av_ds)
    )
  
  if (nrow(dat_cor) < 5 || length(unique(dat_cor$log_av_ds)) < 2) {
    
    orientation <- tibble(
      nmds_scope = unique(scores_df$nmds_scope),
      n_for_orientation = nrow(dat_cor),
      n_decay_values = length(unique(dat_cor$log_av_ds)),
      spearman_cor_NMDS1_decay = NA_real_,
      NMDS1_multiplier = 1,
      orientation_action = "kept_insufficient_decay_gradient"
    )
    
    return(list(scores = scores_df, orientation = orientation))
  }
  
  cor_now <- suppressWarnings(
    cor(
      dat_cor$NMDS1,
      dat_cor$log_av_ds,
      method = "spearman",
      use = "complete.obs"
    )
  )
  
  # Desired orientation: lower decay stages toward the right.
  # Therefore, NMDS1 should correlate negatively with decay stage.
  multiplier <- ifelse(is.na(cor_now), 1, ifelse(cor_now > 0, -1, 1))
  
  scores_df <- scores_df %>%
    mutate(NMDS1 = NMDS1 * multiplier)
  
  orientation <- tibble(
    nmds_scope = unique(scores_df$nmds_scope),
    n_for_orientation = nrow(dat_cor),
    n_decay_values = length(unique(dat_cor$log_av_ds)),
    spearman_cor_NMDS1_decay = cor_now,
    NMDS1_multiplier = multiplier,
    orientation_action = case_when(
      is.na(cor_now) ~ "kept_correlation_na",
      multiplier == -1 ~ "flipped_NMDS1",
      TRUE ~ "kept_NMDS1"
    )
  )
  
  list(scores = scores_df, orientation = orientation)
}

run_forest_nmds <- function(X, meta, dataset_name, forest_now, trymax = 20) {
  
  cat("\n=========================================================\n")
  cat("FOREST-SPECIFIC NMDS\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Forest:", forest_now, "\n")
  cat("=========================================================\n")
  
  meta_f <- meta %>%
    filter(as.character(forest) == forest_now) %>%
    mutate(
      forest = factor(as.character(forest), levels = forest_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      data_scope = factor(dataset_name, levels = dataset_levels),
      log_av_ds = suppressWarnings(as.numeric(log_av_ds))
    )
  
  X_f <- X[meta_f$sample_id, , drop = FALSE]
  
  keep_samples <- rowSums(X_f, na.rm = TRUE) > 0
  X_f <- X_f[keep_samples, , drop = FALSE]
  meta_f <- meta_f[match(rownames(X_f), meta_f$sample_id), , drop = FALSE]
  
  keep_taxa <- colSums(X_f, na.rm = TRUE) > 0
  X_f <- X_f[, keep_taxa, drop = FALSE]
  
  cat("Samples:", nrow(X_f), "\n")
  cat("Taxa:", ncol(X_f), "\n")
  cat("Period table:\n")
  print(table(meta_f$timespan, useNA = "ifany"))
  
  if (nrow(X_f) < 6 || ncol(X_f) < 2) {
    warning("Skipping ", dataset_name, " / ", forest_now, ": too few samples or taxa.")
    return(NULL)
  }
  
  if (length(unique(meta_f$timespan[!is.na(meta_f$timespan)])) < 2) {
    warning("Skipping ", dataset_name, " / ", forest_now, ": only one survey period present.")
    return(NULL)
  }
  
  nmds_f <- metaMDS(
    X_f,
    distance = "jaccard",
    binary = TRUE,
    k = 2,
    trymax = trymax,
    autotransform = FALSE,
    noshare = TRUE,
    trace = TRUE
  )
  
  scores_f <- as.data.frame(scores(nmds_f, display = "sites")) %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    rename(
      NMDS1 = 2,
      NMDS2 = 3
    ) %>%
    left_join(meta_f, by = "sample_id") %>%
    mutate(
      data_scope = factor(dataset_name, levels = dataset_levels),
      forest = factor(as.character(forest), levels = forest_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(
        timespan_labels[as.character(timespan)],
        levels = survey_period_levels
      ),
      forest_timespan = factor(
        paste(forest, timespan, sep = "."),
        levels = names(forest_time_cols)
      ),
      nmds_scope = paste(dataset_name, forest_now, sep = " | "),
      row_panel = factor(forest_now, levels = row_panel_levels),
      panel_type = "Forest-specific NMDS",
      log_av_ds = suppressWarnings(as.numeric(log_av_ds))
    )
  
  orientation_result <- orient_nmds1_by_decay(scores_f)
  scores_f <- orientation_result$scores
  
  orientation_f <- orientation_result$orientation %>%
    mutate(
      data_scope = factor(dataset_name, levels = dataset_levels),
      forest = factor(forest_now, levels = forest_levels)
    )
  
  centroids_f <- scores_f %>%
    group_by(data_scope, row_panel, forest, timespan, forest_timespan, nmds_scope) %>%
    summarise(
      NMDS1 = mean(NMDS1, na.rm = TRUE),
      NMDS2 = mean(NMDS2, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  shifts_f <- centroids_f %>%
    select(data_scope, row_panel, forest, nmds_scope, timespan, NMDS1, NMDS2, n) %>%
    pivot_wider(
      names_from = timespan,
      values_from = c(NMDS1, NMDS2, n)
    ) %>%
    filter(
      !is.na(NMDS1_previous),
      !is.na(NMDS1_current),
      !is.na(NMDS2_previous),
      !is.na(NMDS2_current)
    ) %>%
    mutate(
      d_NMDS1 = NMDS1_current - NMDS1_previous,
      d_NMDS2 = NMDS2_current - NMDS2_previous,
      centroid_shift_distance_2D = sqrt(d_NMDS1^2 + d_NMDS2^2),
      stress = nmds_f$stress
    )
  
  stress_f <- tibble(
    data_scope = factor(dataset_name, levels = dataset_levels),
    forest = factor(forest_now, levels = forest_levels),
    row_panel = factor(forest_now, levels = row_panel_levels),
    nmds_scope = paste(dataset_name, forest_now, sep = " | "),
    panel_type = "Forest-specific NMDS",
    stress = nmds_f$stress,
    n_samples = nrow(X_f),
    n_taxa = ncol(X_f)
  )
  
  cat("Stress:", nmds_f$stress, "\n")
  cat("NMDS1 orientation:\n")
  print(orientation_f, n = Inf)
  cat("Centroid shift:\n")
  print(shifts_f, n = Inf)
  
  list(
    nmds = nmds_f,
    X = X_f,
    meta = meta_f,
    scores = scores_f,
    centroids = centroids_f,
    shifts = shifts_f,
    stress = stress_f,
    orientation = orientation_f
  )
}

fit_decay_surface <- function(scores_df, grouping_cols, grid_n = 120, k_max = 12) {
  
  dat_all <- scores_df %>%
    filter(
      !is.na(NMDS1),
      !is.na(NMDS2),
      !is.na(log_av_ds)
    )
  
  if (nrow(dat_all) == 0) return(tibble())
  
  split_df <- dat_all %>%
    group_by(across(all_of(grouping_cols))) %>%
    group_split()
  
  surface_list <- list()
  
  for (dat_now in split_df) {
    
    if (nrow(dat_now) < 20 || length(unique(dat_now$log_av_ds)) < 3) {
      warning(
        "Decay surface skipped for: ",
        paste(
          paste(grouping_cols, as.character(dat_now[1, grouping_cols]), sep = " = "),
          collapse = "; "
        )
      )
      next
    }
    
    k_now <- min(k_max, max(5, nrow(dat_now) - 1))
    
    gam_now <- mgcv::gam(
      log_av_ds ~ s(NMDS1, NMDS2, k = k_now),
      data = dat_now,
      method = "REML"
    )
    
    grid_now <- expand.grid(
      NMDS1 = seq(
        min(dat_now$NMDS1, na.rm = TRUE),
        max(dat_now$NMDS1, na.rm = TRUE),
        length.out = grid_n
      ),
      NMDS2 = seq(
        min(dat_now$NMDS2, na.rm = TRUE),
        max(dat_now$NMDS2, na.rm = TRUE),
        length.out = grid_n
      )
    ) %>%
      as_tibble()
    
    for (cc in grouping_cols) {
      grid_now[[cc]] <- dat_now[[cc]][1]
    }
    
    grid_now <- grid_now %>%
      mutate(
        log_av_ds_pred = as.numeric(
          predict(gam_now, newdata = ., type = "response")
        ),
        log_av_ds_pred = pmin(pmax(log_av_ds_pred, 1), 5)
      )
    
    surface_list[[length(surface_list) + 1]] <- grid_now
  }
  
  bind_rows(surface_list)
}

make_inset_label_positions <- function(scores_df, group_cols, label_df = NULL) {
  
  pos <- scores_df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      xmin = min(NMDS1, na.rm = TRUE),
      xmax = max(NMDS1, na.rm = TRUE),
      ymin = min(NMDS2, na.rm = TRUE),
      ymax = max(NMDS2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      NMDS1 = xmin + 0.025 * (xmax - xmin),
      NMDS2 = ymax - 0.035 * (ymax - ymin)
    )
  
  if (!is.null(label_df)) {
    pos <- pos %>%
      left_join(label_df, by = group_cols)
  }
  
  pos
}

# =========================================================
# 4. RUN GLOBAL AND FOREST-SPECIFIC NMDS
# =========================================================

res_full <- run_global_nmds(
  X = X_full,
  meta = meta_full,
  dataset_name = "Full dataset",
  trymax = 200
)

res_subset <- run_global_nmds(
  X = X_subset,
  meta = meta_subset,
  dataset_name = "Subsetted",
  trymax = 200
)

global_results <- list(
  "Full dataset" = res_full,
  "Subsetted" = res_subset
)

forest_nmds_results <- list()

for (dataset_now in dataset_levels) {
  
  if (dataset_now == "Full dataset") {
    X_now <- X_full
    meta_now <- meta_full
  } else {
    X_now <- X_subset
    meta_now <- meta_subset
  }
  
  for (forest_now in forest_levels) {
    
    result_name <- paste(dataset_now, forest_now, sep = "__")
    
    forest_nmds_results[[result_name]] <- run_forest_nmds(
      X = X_now,
      meta = meta_now,
      dataset_name = dataset_now,
      forest_now = forest_now,
      trymax = 200
    )
  }
}

forest_nmds_results <- forest_nmds_results[
  !vapply(forest_nmds_results, is.null, logical(1))
]

# =========================================================
# 5. COMBINE NMDS OUTPUTS
# =========================================================

global_scores_all <- bind_rows(lapply(global_results, `[[`, "scores")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor("All forests", levels = row_panel_levels)
  )

global_centroids_all <- bind_rows(lapply(global_results, `[[`, "centroids")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor("All forests", levels = row_panel_levels)
  )

global_shifts_all <- bind_rows(lapply(global_results, `[[`, "shifts")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor("All forests", levels = row_panel_levels)
  )

global_stress_all <- bind_rows(lapply(global_results, `[[`, "stress")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor("All forests", levels = row_panel_levels),
    stress_label = paste0("Stress = ", round(stress, 3))
  )

forest_scores_all <- bind_rows(lapply(forest_nmds_results, `[[`, "scores")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

forest_centroids_all <- bind_rows(lapply(forest_nmds_results, `[[`, "centroids")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

forest_shifts_all <- bind_rows(lapply(forest_nmds_results, `[[`, "shifts")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

forest_stress_all <- bind_rows(lapply(forest_nmds_results, `[[`, "stress")) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels),
    stress_label = paste0("Stress = ", round(stress, 3))
  )

forest_orientation_all <- bind_rows(lapply(forest_nmds_results, `[[`, "orientation"))

cat("\n=========================================================\n")
cat("GLOBAL NMDS STRESS SUMMARY\n")
cat("=========================================================\n")
print(global_stress_all, n = Inf)

cat("\n=========================================================\n")
cat("FOREST-SPECIFIC NMDS STRESS SUMMARY\n")
cat("=========================================================\n")
print(forest_stress_all, n = Inf)

cat("\n=========================================================\n")
cat("FOREST-SPECIFIC AXIS ORIENTATION SUMMARY\n")
cat("=========================================================\n")
print(forest_orientation_all, n = Inf)

# =========================================================
# 6. FULL-DATASET WEIGHTED CENTROID SHIFT DISTANCE
# =========================================================

cat("\n=========================================================\n")
cat("FULL-DATASET NMDS WEIGHTED CENTROID SHIFT DISTANCE\n")
cat("=========================================================\n")

full_nmds_site_scores <- res_full$scores %>%
  mutate(
    forest = factor(as.character(forest), levels = forest_levels),
    timespan = factor(as.character(timespan), levels = timespan_levels)
  )

stopifnot(
  all(c("sample_id", "forest", "timespan", "NMDS1", "NMDS2") %in% names(full_nmds_site_scores)),
  identical(rownames(X_full), full_nmds_site_scores$sample_id)
)

full_jaccard_dist <- vegdist(
  X_full,
  method = "jaccard",
  binary = TRUE
)

axis1_dist <- dist(full_nmds_site_scores$NMDS1)
axis2_dist <- dist(full_nmds_site_scores$NMDS2)

axis_fit <- tibble(
  axis = c("NMDS1", "NMDS2"),
  pearson_r = c(
    cor(as.vector(full_jaccard_dist), as.vector(axis1_dist), method = "pearson"),
    cor(as.vector(full_jaccard_dist), as.vector(axis2_dist), method = "pearson")
  )
) %>%
  mutate(
    axis_r2 = pearson_r^2,
    axis_weight = axis_r2 / sum(axis_r2, na.rm = TRUE)
  )

cat("\nPost-hoc NMDS axis weights based on correlation with original Jaccard distances:\n")
print(axis_fit, n = Inf)

w1 <- axis_fit$axis_weight[axis_fit$axis == "NMDS1"]
w2 <- axis_fit$axis_weight[axis_fit$axis == "NMDS2"]

full_centroids_for_weighted_shift <- full_nmds_site_scores %>%
  group_by(forest, timespan) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  )

full_nmds_weighted_centroid_shift <- full_centroids_for_weighted_shift %>%
  pivot_wider(
    names_from = timespan,
    values_from = c(NMDS1, NMDS2, n_samples)
  ) %>%
  filter(
    !is.na(NMDS1_previous),
    !is.na(NMDS1_current),
    !is.na(NMDS2_previous),
    !is.na(NMDS2_current)
  ) %>%
  mutate(
    d_NMDS1 = NMDS1_current - NMDS1_previous,
    d_NMDS2 = NMDS2_current - NMDS2_previous,
    unweighted_2D_shift = sqrt(d_NMDS1^2 + d_NMDS2^2),
    weighted_2D_shift = sqrt(w1 * d_NMDS1^2 + w2 * d_NMDS2^2),
    NMDS1_axis_weight = w1,
    NMDS2_axis_weight = w2,
    NMDS_stress = res_full$nmds$stress
  ) %>%
  arrange(desc(weighted_2D_shift))

cat("\nFull-dataset historical-to-current centroid shift per forest:\n")
print(full_nmds_weighted_centroid_shift, n = Inf)

cat("\nInterpretation note:\n")
cat(
  "Weighted distances use post-hoc NMDS axis weights, not true explained variance.\n",
  "This is useful as a descriptive comparison of centroid displacement, but NMDS axes remain arbitrary rotations of a rank-based ordination.\n",
  sep = ""
)

# =========================================================
# 7. DECAY SURFACES FOR SITE NMDS PANELS
# =========================================================

global_decay_surface_all <- fit_decay_surface(
  scores_df = global_scores_all,
  grouping_cols = c("data_scope", "row_panel", "nmds_scope"),
  grid_n = 140,
  k_max = 20
)

forest_decay_surface_all <- fit_decay_surface(
  scores_df = forest_scores_all,
  grouping_cols = c("data_scope", "row_panel", "forest", "nmds_scope"),
  grid_n = 120,
  k_max = 12
)

surface_10panel <- bind_rows(
  global_decay_surface_all,
  forest_decay_surface_all
) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

scores_10panel <- bind_rows(
  global_scores_all,
  forest_scores_all
) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels),
    forest_timespan = factor(
      paste(forest, timespan, sep = "."),
      levels = names(forest_time_cols)
    )
  )

centroids_10panel <- bind_rows(
  global_centroids_all,
  forest_centroids_all
) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels),
    forest_timespan = factor(
      paste(forest, timespan, sep = "."),
      levels = names(forest_time_cols)
    )
  )

shifts_10panel <- bind_rows(
  global_shifts_all,
  forest_shifts_all
) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

stress_10panel <- bind_rows(
  global_stress_all,
  forest_stress_all
) %>%
  mutate(
    data_scope = factor(as.character(data_scope), levels = dataset_levels),
    row_panel = factor(as.character(row_panel), levels = row_panel_levels)
  )

stress_pos_10panel <- make_inset_label_positions(
  scores_df = scores_10panel,
  group_cols = c("row_panel", "data_scope", "nmds_scope"),
  label_df = stress_10panel %>%
    select(row_panel, data_scope, nmds_scope, stress_label)
)

# =========================================================
# 8. THEMES AND GUIDES
# =========================================================

theme_2panel <- theme_bw(base_size = 8.8) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15, colour = "grey90"),
    axis.text = element_text(size = 7.2),
    axis.title = element_text(size = 8.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.36, "cm"),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, colour = "grey25"),
    plot.caption = element_text(hjust = 0, size = 6.8, colour = "grey25"),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

theme_nmds_10panel <- theme_bw(base_size = 7.4) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    strip.text.x = element_text(face = "bold", size = 7.4),
    strip.text.y = element_text(face = "bold", size = 7.4, angle = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.12, colour = "grey90"),
    axis.text = element_text(size = 5.8),
    axis.title = element_text(size = 7.0),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.0),
    legend.text = element_text(size = 6.3),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.35, "cm"),
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 8.3, colour = "grey25"),
    plot.caption = element_text(hjust = 0, size = 6.8, colour = "grey25"),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

legend_guides_10panel <- guides(
  shape = guide_legend(
    order = 1,
    nrow = 1,
    override.aes = list(
      shape = c(21, 24),
      fill = "grey80",
      colour = "black",
      alpha = 1,
      size = 2.6,
      stroke = 0.55
    )
  ),
  fill = guide_legend(
    order = 2,
    nrow = 2,
    override.aes = list(
      shape = rep(21, length(forest_time_cols)),
      fill = unname(forest_time_cols),
      colour = "black",
      alpha = 1,
      size = 2.4,
      stroke = 0.45
    )
  )
)

species_guides_10panel <- guides(
  fill = guide_legend(
    order = 1,
    nrow = 2,
    override.aes = list(
      shape = 21,
      colour = "grey20",
      size = 2.3,
      alpha = 1,
      stroke = 0.25
    )
  ),
  size = guide_legend(
    order = 2,
    nrow = 1,
    override.aes = list(
      shape = 21,
      fill = "grey60",
      colour = "grey20",
      alpha = 1,
      stroke = 0.25
    )
  )
)

# =========================================================
# 9. FIGURE 1: 2-PANEL GLOBAL FULL VS SUBSETTED NMDS
# =========================================================

stress_pos_2panel <- make_inset_label_positions(
  scores_df = global_scores_all,
  group_cols = c("data_scope", "nmds_scope"),
  label_df = global_stress_all %>%
    select(data_scope, nmds_scope, stress_label)
)

p_nmds_global_2panel <- ggplot(
  global_scores_all,
  aes(x = NMDS1, y = NMDS2)
) +
  geom_contour(
    data = global_decay_surface_all,
    aes(
      x = NMDS1,
      y = NMDS2,
      z = log_av_ds_pred
    ),
    inherit.aes = FALSE,
    breaks = decay_breaks,
    colour = "grey45",
    linewidth = 0.25,
    alpha = 0.75
  ) +
  geom_point(
    aes(
      fill = forest_timespan,
      shape = timespan
    ),
    colour = "grey25",
    size = 1.25,
    alpha = 0.45,
    stroke = 0.18
  ) +
  geom_segment(
    data = global_shifts_all,
    aes(
      x = NMDS1_previous,
      y = NMDS2_previous,
      xend = NMDS1_current,
      yend = NMDS2_current,
      colour = forest
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    arrow = arrow(length = unit(0.10, "cm")),
    show.legend = FALSE
  ) +
  geom_point(
    data = global_centroids_all,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = forest_timespan,
      shape = timespan
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 2.3,
    stroke = 0.75,
    show.legend = FALSE
  ) +
  geom_text(
    data = stress_pos_2panel,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = stress_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.3,
    colour = "grey15"
  ) +
  facet_grid(. ~ data_scope, scales = "fixed") +
  scale_fill_manual(
    values = forest_time_cols,
    labels = forest_time_labels,
    name = "Forest × period",
    drop = FALSE
  ) +
  scale_colour_manual(
    values = forest_cols,
    guide = "none",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = timespan_shape_vals,
    labels = timespan_labels,
    name = "Survey period",
    drop = FALSE
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 1,
      override.aes = list(
        shape = c(21, 24),
        fill = "grey80",
        colour = "black",
        alpha = 1,
        size = 2.7,
        stroke = 0.6
      )
    ),
    fill = guide_legend(
      order = 2,
      nrow = 2,
      override.aes = list(
        shape = rep(21, length(forest_time_cols)),
        fill = unname(forest_time_cols),
        colour = "black",
        alpha = 1,
        size = 2.4,
        stroke = 0.45
      )
    )
  ) +
  coord_equal() +
  labs(
    title = "Temporal shifts in deadwood-associated fungal assemblages",
    subtitle = "Global two-dimensional Jaccard NMDS for the full dataset and substrate-matched subset",
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_2panel

print(p_nmds_global_2panel)

save_png(
  plot = p_nmds_global_2panel,
  filename = "nmds_2panel_global_full_vs_subsetted.png",
  width = a4_width_in,
  height = a4_height_2panel
)

# =========================================================
# 10. FIGURE 2: 10-PANEL SITE NMDS
# =========================================================

p_nmds_10panel_sites <- ggplot(
  scores_10panel,
  aes(x = NMDS1, y = NMDS2)
) +
  geom_contour(
    data = surface_10panel,
    aes(
      x = NMDS1,
      y = NMDS2,
      z = log_av_ds_pred
    ),
    inherit.aes = FALSE,
    breaks = decay_breaks,
    colour = "grey45",
    linewidth = 0.20,
    alpha = 0.8
  ) +
  {
    if (use_contour_text) {
      metR::geom_text_contour(
        data = surface_10panel,
        aes(
          x = NMDS1,
          y = NMDS2,
          z = log_av_ds_pred,
          label = after_stat(paste0(floor(level), "|", ceiling(level)))
        ),
        inherit.aes = FALSE,
        breaks = decay_breaks,
        stroke = 0.08,
        size = 1.65,
        colour = "grey20",
        check_overlap = TRUE,
        skip = 0
      )
    }
  } +
  geom_point(
    aes(
      fill = forest_timespan,
      shape = timespan
    ),
    colour = "grey25",
    size = 0.95,
    alpha = 0.62,
    stroke = 0.16
  ) +
  geom_segment(
    data = shifts_10panel,
    aes(
      x = NMDS1_previous,
      y = NMDS2_previous,
      xend = NMDS1_current,
      yend = NMDS2_current,
      colour = forest
    ),
    inherit.aes = FALSE,
    linewidth = 0.42,
    arrow = arrow(length = unit(0.075, "cm")),
    show.legend = FALSE
  ) +
  geom_point(
    data = centroids_10panel,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = forest_timespan,
      shape = timespan
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 1.95,
    stroke = 0.65,
    show.legend = FALSE
  ) +
  geom_text(
    data = stress_pos_10panel,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = stress_label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 1.95,
    colour = "grey15"
  ) +
  facet_grid(
    row_panel ~ data_scope,
    scales = "free",
    drop = FALSE
  ) +
  scale_fill_manual(
    values = forest_time_cols,
    labels = forest_time_labels,
    name = "Forest × period",
    drop = FALSE
  ) +
  scale_colour_manual(
    values = forest_cols,
    guide = "none",
    drop = FALSE
  ) +
  scale_shape_manual(
    values = timespan_shape_vals,
    labels = timespan_labels,
    name = "Survey period",
    drop = FALSE
  ) +
  legend_guides_10panel +
  coord_cartesian(clip = "off") +
  labs(
    title = "Temporal shifts in deadwood-associated fungal assemblages",
    subtitle = "Top row shows global NMDS; lower rows show independent forest-specific NMDS solutions",
    x = "NMDS1",
    y = "NMDS2",
    caption = paste(
      "Grey isolines show fitted decay-stage boundaries at half-stage values.",
      "Large symbols are forest × period centroids.",
      "Arrows indicate centroid displacement from Historical to Current.",
      "Forest-specific axes are interpretable within panels."
    )
  ) +
  theme_nmds_10panel

print(p_nmds_10panel_sites)

save_png(
  plot = p_nmds_10panel_sites,
  filename = "nmds_10panel_sites_full_vs_subsetted.png",
  width = a4_width_in,
  height = a4_height_10panel
)

# =========================================================
# 11. TAXONOMY AND GUILD TABLE
# =========================================================

guild_df <- read_excel(
  "taxonomy_revisedJHC.xlsx",
  sheet = "eco_target"
) %>%
  mutate(
    Taxon = str_trim(as.character(Taxon)),
    GUILD = str_trim(as.character(GUILD)),
    GUILD = case_when(
      GUILD %in% c("unknown/different", "Unknown/different") ~ "Unknown/different",
      TRUE ~ GUILD
    )
  ) %>%
  filter(
    code == "L",
    TARGET == 1,
    !is.na(Taxon),
    !is.na(GUILD),
    GUILD != "-"
  ) %>%
  mutate(
    GUILD = factor(GUILD, levels = guild_levels)
  ) %>%
  filter(!is.na(GUILD)) %>%
  select(Taxon, GUILD) %>%
  distinct()

cat("\n=========================================================\n")
cat("GUILD TABLE\n")
cat("=========================================================\n")
print(table(guild_df$GUILD, useNA = "ifany"))

# =========================================================
# 12. SPECIES-SCORE PANEL FUNCTIONS
# =========================================================

extract_site_scores_for_species <- function(nmds_obj,
                                            meta,
                                            row_panel,
                                            nmds_scope,
                                            panel_type,
                                            stress_value) {
  
  as.data.frame(scores(nmds_obj, display = "sites")) %>%
    rownames_to_column("sample_id") %>%
    as_tibble() %>%
    rename(
      NMDS1 = 2,
      NMDS2 = 3
    ) %>%
    left_join(meta, by = "sample_id") %>%
    mutate(
      row_panel = factor(row_panel, levels = row_panel_levels),
      nmds_scope = nmds_scope,
      panel_type = panel_type,
      stress = stress_value,
      forest = str_replace_all(as.character(forest), "Strodam", "Strødam"),
      forest = factor(forest, levels = forest_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(
        timespan_labels[as.character(timespan)],
        levels = survey_period_levels
      ),
      log_av_ds = suppressWarnings(as.numeric(log_av_ds))
    )
}

extract_species_scores_base <- function(nmds_obj,
                                        row_panel,
                                        nmds_scope,
                                        panel_type,
                                        forest_now = NA_character_) {
  
  sp <- scores(nmds_obj, display = "species")
  
  if (is.null(sp) || nrow(sp) == 0) {
    warning("No species scores available for: ", nmds_scope)
    return(NULL)
  }
  
  as.data.frame(sp) %>%
    rownames_to_column("Taxon") %>%
    as_tibble() %>%
    rename(
      NMDS1 = 2,
      NMDS2 = 3
    ) %>%
    mutate(
      Taxon = str_trim(as.character(Taxon)),
      row_panel = factor(row_panel, levels = row_panel_levels),
      nmds_scope = nmds_scope,
      panel_type = panel_type,
      forest = factor(forest_now, levels = forest_levels)
    )
}

calculate_species_occupancy_by_period <- function(X, meta) {
  
  X %>%
    as_tibble(rownames = "sample_id") %>%
    pivot_longer(
      cols = -sample_id,
      names_to = "Taxon",
      values_to = "pa"
    ) %>%
    filter(pa > 0) %>%
    left_join(
      meta %>%
        select(sample_id, timespan),
      by = "sample_id"
    ) %>%
    mutate(
      Taxon = str_trim(as.character(Taxon)),
      timespan = factor(as.character(timespan), levels = timespan_levels)
    ) %>%
    group_by(Taxon, timespan) %>%
    summarise(
      n_occ = n(),
      .groups = "drop"
    ) %>%
    left_join(
      meta %>%
        count(timespan, name = "n_samples_period"),
      by = "timespan"
    ) %>%
    mutate(
      occupancy_prop = n_occ / n_samples_period,
      survey_period = factor(
        timespan_labels[as.character(timespan)],
        levels = survey_period_levels
      )
    )
}

fit_decay_surface_by_period <- function(site_scores_df,
                                        row_panel_now,
                                        nmds_scope_now,
                                        grid_n = 120,
                                        k_max = 12) {
  
  surface_list <- list()
  
  for (tp in timespan_levels) {
    
    dat_now <- site_scores_df %>%
      filter(
        timespan == tp,
        !is.na(NMDS1),
        !is.na(NMDS2),
        !is.na(log_av_ds)
      )
    
    if (nrow(dat_now) < 20 || length(unique(dat_now$log_av_ds)) < 3) {
      warning(
        "Decay surface skipped for ",
        nmds_scope_now,
        " / ",
        tp,
        ": insufficient decay-stage variation."
      )
      next
    }
    
    k_now <- min(k_max, max(5, nrow(dat_now) - 1))
    
    gam_now <- mgcv::gam(
      log_av_ds ~ s(NMDS1, NMDS2, k = k_now),
      data = dat_now,
      method = "REML"
    )
    
    grid_now <- expand.grid(
      NMDS1 = seq(
        min(dat_now$NMDS1, na.rm = TRUE),
        max(dat_now$NMDS1, na.rm = TRUE),
        length.out = grid_n
      ),
      NMDS2 = seq(
        min(dat_now$NMDS2, na.rm = TRUE),
        max(dat_now$NMDS2, na.rm = TRUE),
        length.out = grid_n
      )
    ) %>%
      as_tibble() %>%
      mutate(
        row_panel = factor(row_panel_now, levels = row_panel_levels),
        nmds_scope = nmds_scope_now,
        timespan = factor(tp, levels = timespan_levels),
        survey_period = factor(
          timespan_labels[tp],
          levels = survey_period_levels
        ),
        log_av_ds_pred = as.numeric(
          predict(gam_now, newdata = ., type = "response")
        ),
        log_av_ds_pred = pmin(pmax(log_av_ds_pred, 1), 5)
      )
    
    surface_list[[tp]] <- grid_now
  }
  
  bind_rows(surface_list)
}

make_species_stress_positions <- function(site_scores_df) {
  
  site_scores_df %>%
    group_by(row_panel, survey_period, nmds_scope) %>%
    summarise(
      xmin = min(NMDS1, na.rm = TRUE),
      xmax = max(NMDS1, na.rm = TRUE),
      ymin = min(NMDS2, na.rm = TRUE),
      ymax = max(NMDS2, na.rm = TRUE),
      stress = unique(stress)[1],
      .groups = "drop"
    ) %>%
    mutate(
      NMDS1 = xmin + 0.025 * (xmax - xmin),
      NMDS2 = ymax - 0.035 * (ymax - ymin),
      stress_label = paste0("Stress = ", round(stress, 3))
    )
}

make_species_n_positions <- function(site_scores_df) {
  
  site_scores_df %>%
    count(row_panel, survey_period, nmds_scope, name = "n_samples") %>%
    left_join(
      site_scores_df %>%
        group_by(row_panel, survey_period, nmds_scope) %>%
        summarise(
          xmin = min(NMDS1, na.rm = TRUE),
          xmax = max(NMDS1, na.rm = TRUE),
          ymin = min(NMDS2, na.rm = TRUE),
          ymax = max(NMDS2, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("row_panel", "survey_period", "nmds_scope")
    ) %>%
    mutate(
      NMDS1 = xmax - 0.025 * (xmax - xmin),
      NMDS2 = ymin + 0.035 * (ymax - ymin),
      n_label = paste0("n = ", n_samples)
    )
}

build_species_panel_data <- function(dataset_name) {
  
  if (dataset_name == "Full dataset") {
    X_main <- X_full
    meta_main <- meta_full
    res_global <- res_full
  } else {
    X_main <- X_subset
    meta_main <- meta_subset
    res_global <- res_subset
  }
  
  global_scope <- paste(dataset_name, "All forests", sep = " | ")
  
  site_global <- extract_site_scores_for_species(
    nmds_obj = res_global$nmds,
    meta = meta_main,
    row_panel = "All forests",
    nmds_scope = global_scope,
    panel_type = "Global NMDS",
    stress_value = res_global$nmds$stress
  )
  
  species_global_base <- extract_species_scores_base(
    nmds_obj = res_global$nmds,
    row_panel = "All forests",
    nmds_scope = global_scope,
    panel_type = "Global NMDS",
    forest_now = NA_character_
  )
  
  species_global_occ <- calculate_species_occupancy_by_period(
    X = X_main,
    meta = meta_main
  ) %>%
    mutate(
      row_panel = factor("All forests", levels = row_panel_levels),
      nmds_scope = global_scope
    )
  
  surface_global <- fit_decay_surface_by_period(
    site_scores_df = site_global,
    row_panel_now = "All forests",
    nmds_scope_now = global_scope,
    grid_n = 140,
    k_max = 20
  )
  
  forest_results <- forest_nmds_results[
    vapply(
      forest_nmds_results,
      function(z) {
        !is.null(z) &&
          as.character(z$stress$data_scope[1]) == dataset_name
      },
      logical(1)
    )
  ]
  
  site_forest_list <- list()
  species_forest_base_list <- list()
  species_forest_occ_list <- list()
  surface_forest_list <- list()
  
  for (nm in names(forest_results)) {
    
    res_now <- forest_results[[nm]]
    
    forest_now <- as.character(res_now$stress$forest[1])
    scope_now  <- as.character(res_now$stress$nmds_scope[1])
    
    meta_now <- res_now$meta %>%
      mutate(
        sample_id = as.character(sample_id),
        forest = str_replace_all(as.character(forest), "Strodam", "Strødam"),
        forest = factor(forest, levels = forest_levels),
        timespan = factor(as.character(timespan), levels = timespan_levels),
        survey_period = factor(
          timespan_labels[as.character(timespan)],
          levels = survey_period_levels
        ),
        log_av_ds = suppressWarnings(as.numeric(log_av_ds))
      )
    
    X_now <- as.data.frame(res_now$X)
    X_now[] <- lapply(X_now, function(z) as.integer(as.numeric(z) > 0))
    
    site_now <- extract_site_scores_for_species(
      nmds_obj = res_now$nmds,
      meta = meta_now,
      row_panel = forest_now,
      nmds_scope = scope_now,
      panel_type = "Forest-specific NMDS",
      stress_value = res_now$nmds$stress
    )
    
    oriented_sites <- res_now$scores %>%
      select(sample_id, NMDS1, NMDS2)
    
    site_now <- site_now %>%
      select(-NMDS1, -NMDS2) %>%
      left_join(oriented_sites, by = "sample_id") %>%
      relocate(NMDS1, NMDS2, .after = sample_id)
    
    species_base_now <- extract_species_scores_base(
      nmds_obj = res_now$nmds,
      row_panel = forest_now,
      nmds_scope = scope_now,
      panel_type = "Forest-specific NMDS",
      forest_now = forest_now
    )
    
    multiplier_now <- res_now$orientation$NMDS1_multiplier[1]
    
    species_base_now <- species_base_now %>%
      mutate(NMDS1 = NMDS1 * multiplier_now)
    
    species_occ_now <- calculate_species_occupancy_by_period(
      X = X_now,
      meta = meta_now
    ) %>%
      mutate(
        row_panel = factor(forest_now, levels = row_panel_levels),
        nmds_scope = scope_now
      )
    
    surface_now <- fit_decay_surface_by_period(
      site_scores_df = site_now,
      row_panel_now = forest_now,
      nmds_scope_now = scope_now,
      grid_n = 120,
      k_max = 12
    )
    
    site_forest_list[[nm]] <- site_now
    species_forest_base_list[[nm]] <- species_base_now
    species_forest_occ_list[[nm]] <- species_occ_now
    surface_forest_list[[nm]] <- surface_now
  }
  
  site_scores <- bind_rows(
    site_global,
    bind_rows(site_forest_list)
  ) %>%
    mutate(
      row_panel = factor(as.character(row_panel), levels = row_panel_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(as.character(survey_period), levels = survey_period_levels)
    )
  
  species_base <- bind_rows(
    species_global_base,
    bind_rows(species_forest_base_list)
  ) %>%
    mutate(
      row_panel = factor(as.character(row_panel), levels = row_panel_levels)
    )
  
  species_occ <- bind_rows(
    species_global_occ,
    bind_rows(species_forest_occ_list)
  ) %>%
    mutate(
      row_panel = factor(as.character(row_panel), levels = row_panel_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(as.character(survey_period), levels = survey_period_levels)
    )
  
  surface_species <- bind_rows(
    surface_global,
    bind_rows(surface_forest_list)
  ) %>%
    mutate(
      row_panel = factor(as.character(row_panel), levels = row_panel_levels),
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(as.character(survey_period), levels = survey_period_levels)
    )
  
  species_scores <- species_base %>%
    left_join(
      species_occ,
      by = c("Taxon", "row_panel", "nmds_scope")
    ) %>%
    left_join(
      guild_df,
      by = "Taxon"
    ) %>%
    mutate(
      timespan = factor(as.character(timespan), levels = timespan_levels),
      survey_period = factor(as.character(survey_period), levels = survey_period_levels),
      GUILD = factor(as.character(GUILD), levels = guild_levels)
    )
  
  species_scores_use <- species_scores %>%
    filter(
      !is.na(timespan),
      !is.na(GUILD),
      !is.na(n_occ),
      n_occ >= species_min_occ
    )
  
  species_labels <- species_scores_use %>%
    group_by(row_panel, survey_period, nmds_scope) %>%
    arrange(desc(n_occ), desc(occupancy_prop), .by_group = TRUE) %>%
    slice_head(n = label_top_n_per_panel) %>%
    ungroup() %>%
    resolve_abbreviated_label_overlaps()
  
  stress_pos <- make_species_stress_positions(site_scores)
  n_pos <- make_species_n_positions(site_scores)
  
  cat("\n=========================================================\n")
  cat("SPECIES PANEL SUMMARY:", toupper(dataset_name), "\n")
  cat("=========================================================\n")
  
  print(
    species_scores_use %>%
      group_by(row_panel, survey_period, GUILD) %>%
      summarise(
        n_species = n(),
        min_occ = min(n_occ, na.rm = TRUE),
        max_occ = max(n_occ, na.rm = TRUE),
        .groups = "drop"
      ),
    n = Inf
  )
  
  cat("\nUnmatched species without guild:", dataset_name, "\n")
  print(
    species_scores %>%
      filter(is.na(GUILD)) %>%
      distinct(Taxon) %>%
      arrange(Taxon),
    n = Inf
  )
  
  list(
    site_scores = site_scores,
    species_scores = species_scores,
    species_scores_use = species_scores_use,
    species_labels = species_labels,
    surface_species = surface_species,
    stress_pos = stress_pos,
    n_pos = n_pos
  )
}

plot_species_panel <- function(panel_data, dataset_name) {
  
  subtitle_text <- if (dataset_name == "Full dataset") {
    "Full dataset; columns show Historical and Current occurrences within the same row-specific ordination space"
  } else {
    "Subsetted dataset; columns show Historical and Current occurrences within the same row-specific ordination space"
  }
  
  p <- ggplot() +
    {
      if (show_site_cloud) {
        geom_point(
          data = panel_data$site_scores,
          aes(
            x = NMDS1,
            y = NMDS2
          ),
          inherit.aes = FALSE,
          colour = "grey78",
          size = 0.32,
          alpha = 0.25,
          stroke = 0
        )
      }
    } +
    geom_contour(
      data = panel_data$surface_species,
      aes(
        x = NMDS1,
        y = NMDS2,
        z = log_av_ds_pred
      ),
      inherit.aes = FALSE,
      breaks = decay_breaks,
      colour = "grey45",
      linewidth = 0.18,
      alpha = 0.65
    ) +
    {
      if (use_contour_text) {
        metR::geom_text_contour(
          data = panel_data$surface_species,
          aes(
            x = NMDS1,
            y = NMDS2,
            z = log_av_ds_pred,
            label = after_stat(paste0(floor(level), "|", ceiling(level)))
          ),
          inherit.aes = FALSE,
          breaks = decay_breaks,
          stroke = 0.08,
          size = 1.55,
          colour = "grey30",
          check_overlap = TRUE,
          skip = 0
        )
      }
    } +
    geom_point(
      data = panel_data$species_scores_use,
      aes(
        x = NMDS1,
        y = NMDS2,
        fill = GUILD,
        size = occupancy_prop
      ),
      inherit.aes = FALSE,
      shape = 21,
      colour = "grey15",
      stroke = 0.14,
      alpha = 0.88
    )
  
  if (use_species_label_callouts) {
    
    p <- p +
      ggrepel::geom_label_repel(
        data = panel_data$species_labels,
        aes(
          x = NMDS1,
          y = NMDS2,
          label = Taxon_label
        ),
        inherit.aes = FALSE,
        size = 1.35,
        colour = "black",
        fill = "white",
        alpha = 0.88,
        label.size = 0.10,
        label.padding = unit(0.06, "lines"),
        segment.size = 0.12,
        segment.alpha = 0.55,
        segment.colour = "grey25",
        box.padding = species_label_box_padding,
        point.padding = species_label_point_padding,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 42
      )
    
  } else {
    
    p <- p +
      ggrepel::geom_text_repel(
        data = panel_data$species_labels,
        aes(
          x = NMDS1,
          y = NMDS2,
          label = Taxon_label
        ),
        inherit.aes = FALSE,
        size = 1.35,
        colour = "black",
        segment.size = 0.10,
        segment.alpha = 0.42,
        box.padding = species_label_box_padding,
        point.padding = species_label_point_padding,
        min.segment.length = 0,
        max.overlaps = Inf,
        seed = 42
      )
  }
  
  p +
    geom_text(
      data = panel_data$stress_pos,
      aes(
        x = NMDS1,
        y = NMDS2,
        label = stress_label
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 1.65,
      colour = "grey15"
    ) +
    geom_text(
      data = panel_data$n_pos,
      aes(
        x = NMDS1,
        y = NMDS2,
        label = n_label
      ),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 1.60,
      colour = "grey20"
    ) +
    facet_grid(
      row_panel ~ survey_period,
      scales = "free",
      drop = FALSE
    ) +
    scale_fill_manual(
      values = guild_cols,
      labels = guild_labeller,
      drop = FALSE,
      name = "Guild"
    ) +
    scale_size_continuous(
      range = c(0.55, 2.55),
      name = "Occupancy"
    ) +
    species_guides_10panel +
    coord_cartesian(clip = "off") +
    labs(
      title = "Species-score NMDS of deadwood-associated fungal assemblages",
      subtitle = subtitle_text,
      x = "NMDS1",
      y = "NMDS2",
      caption = paste(
        "Species positions are weighted-average scores extracted from the corresponding site NMDS.",
        "The top row uses the pooled global NMDS; lower rows use independent forest-specific NMDS solutions.",
        "Species are shown only in survey periods where they occur.",
        "Point size indicates occupancy proportion within the corresponding row and survey period.",
        "Labels abbreviate genus names unless this creates duplicated labels within a panel.",
        "Grey isolines show fitted decay-stage boundaries at half-stage values."
      )
    ) +
    theme_nmds_10panel
}

# =========================================================
# 13. FIGURE 3: FULL-DATASET SPECIES 10-PANEL NMDS
# =========================================================

species_panel_full <- build_species_panel_data("Full dataset")

p_species_nmds_10panel_full <- plot_species_panel(
  panel_data = species_panel_full,
  dataset_name = "Full dataset"
)

print(p_species_nmds_10panel_full)

save_png(
  plot = p_species_nmds_10panel_full,
  filename = "nmds_10panel_species_full_dataset_by_guild.png",
  width = a4_width_in,
  height = a4_height_10panel
)

# =========================================================
# 14. FIGURE 4: SUBSETTED SPECIES 10-PANEL NMDS
# =========================================================

species_panel_subset <- build_species_panel_data("Subsetted")

p_species_nmds_10panel_subset <- plot_species_panel(
  panel_data = species_panel_subset,
  dataset_name = "Subsetted"
)

print(p_species_nmds_10panel_subset)

save_png(
  plot = p_species_nmds_10panel_subset,
  filename = "nmds_10panel_species_subsetted_by_guild.png",
  width = a4_width_in,
  height = a4_height_10panel
)

# =========================================================
# 15. FINAL FILE SUMMARY
# =========================================================

cat("\n=========================================================\n")
cat("PNG FILES WRITTEN\n")
cat("=========================================================\n")

print(
  tibble(
    file = c(
      "nmds_2panel_global_full_vs_subsetted.png",
      "nmds_10panel_sites_full_vs_subsetted.png",
      "nmds_10panel_species_full_dataset_by_guild.png",
      "nmds_10panel_species_subsetted_by_guild.png"
    ),
    path = file.path(
      outdir_nmds,
      c(
        "nmds_2panel_global_full_vs_subsetted.png",
        "nmds_10panel_sites_full_vs_subsetted.png",
        "nmds_10panel_species_full_dataset_by_guild.png",
        "nmds_10panel_species_subsetted_by_guild.png"
      )
    ),
    width_in = c(
      a4_width_in,
      a4_width_in,
      a4_width_in,
      a4_width_in
    ),
    height_in = c(
      a4_height_2panel,
      a4_height_10panel,
      a4_height_10panel,
      a4_height_10panel
    ),
    dpi = 900
  ),
  n = Inf
)
# ------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------
# =========================================================
# FIXED EXPORT OF EACH NMDS STORY PANEL SEPARATELY
#
# Fix:
#   Does NOT assume that data_scope levels are exactly
#   "Full dataset" and "Subsetted".
#
# Output:
#   Separate transparent PNGs for manual assembly in Canva.
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(grid)
})

# ---------------------------------------------------------
# 1. Output directory
# ---------------------------------------------------------

single_panel_dir <- file.path(outdir_nmds, "nmds_story_single_panels_FIXED")
dir.create(single_panel_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 2. Export settings
# ---------------------------------------------------------

single_panel_width_cm  <- 4.4
single_panel_height_cm <- 4.4
single_panel_dpi       <- 900

single_panel_width_in  <- single_panel_width_cm / 2.54
single_panel_height_in <- single_panel_height_cm / 2.54

# ---------------------------------------------------------
# Derive common plot limits from the global all-forest NMDS
# ---------------------------------------------------------

all_plot_coords <- bind_rows(
  global_scores_base %>%
    transmute(
      x = NMDS1,
      y = NMDS2
    ),
  global_centroids_base %>%
    transmute(
      x = NMDS1,
      y = NMDS2
    ),
  global_shifts_base %>%
    transmute(
      x = NMDS1_previous,
      y = NMDS2_previous
    ),
  global_shifts_base %>%
    transmute(
      x = NMDS1_current,
      y = NMDS2_current
    )
)

x_rng_story <- range(all_plot_coords$x, na.rm = TRUE)
y_rng_story <- range(all_plot_coords$y, na.rm = TRUE)

# small padding so points/arrows do not sit on the border
x_pad_story <- diff(x_rng_story) * 0.04
y_pad_story <- diff(y_rng_story) * 0.04

# in case the range is extremely small for some reason
if (!is.finite(x_pad_story) || x_pad_story == 0) x_pad_story <- 0.05
if (!is.finite(y_pad_story) || y_pad_story == 0) y_pad_story <- 0.05

x_lim_story <- c(
  x_rng_story[1] - x_pad_story,
  x_rng_story[2] + x_pad_story
)

y_lim_story <- c(
  y_rng_story[1] - y_pad_story,
  y_rng_story[2] + y_pad_story
)

cat("\nCommon x limits used for all exported panels:\n")
print(x_lim_story)

cat("\nCommon y limits used for all exported panels:\n")
print(y_lim_story)
# ---------------------------------------------------------
# 3. Detect actual dataset names
# ---------------------------------------------------------

actual_scopes <- global_scores_all %>%
  mutate(data_scope_chr = as.character(data_scope)) %>%
  distinct(data_scope_chr) %>%
  pull(data_scope_chr)

cat("\nActual data_scope values in global_scores_all:\n")
print(actual_scopes)

full_scope <- actual_scopes[str_detect(str_to_lower(actual_scopes), "full")][1]
subset_scope <- actual_scopes[str_detect(str_to_lower(actual_scopes), "subset|match")][1]

if (is.na(full_scope) || is.na(subset_scope)) {
  stop(
    "Could not automatically identify full_scope and subset_scope. ",
    "Inspect printed actual_scopes and set full_scope/subset_scope manually."
  )
}

cat("\nUsing:\n")
cat("full_scope   =", full_scope, "\n")
cat("subset_scope =", subset_scope, "\n")

# ---------------------------------------------------------
# 4. Panel map
# ---------------------------------------------------------

single_panel_map <- tibble::tribble(
  ~panel_code,          ~data_scope_raw, ~data_scope_label, ~forest_focus,  ~panel_title,
  "full_strodam",       full_scope,       "Full",           "Strødam",      "Full · Strødam",
  "full_rajhenav",      full_scope,       "Full",           "Rajhenav",     "Full · Rajhenav",
  "subset_rajhenav",    subset_scope,     "Subset",         "Rajhenav",     "Subset · Rajhenav",
  "subset_strodam",     subset_scope,     "Subset",         "Strødam",      "Subset · Strødam",
  "full_suserup",       full_scope,       "Full",           "Suserup",      "Full · Suserup",
  "full_all_forests",   full_scope,       "Full",           "All forests",  "Full · all forests",
  "subset_all_forests", subset_scope,     "Subset",         "All forests",  "Subset · all forests",
  "subset_suserup",     subset_scope,     "Subset",         "Suserup",      "Subset · Suserup",
  "full_sonian",        full_scope,       "Full",           "Sonian",       "Full · Sonian",
  "subset_sonian",      subset_scope,     "Subset",         "Sonian",       "Subset · Sonian"
)

# ---------------------------------------------------------
# 5. Clean NMDS objects
# ---------------------------------------------------------

global_scores_base <- global_scores_all %>%
  mutate(
    data_scope_chr = as.character(data_scope),
    forest_chr = as.character(forest),
    timespan_chr = as.character(timespan),
    forest_timespan = paste(forest_chr, timespan_chr, sep = ".")
  ) %>%
  select(-any_of("row_panel"))

global_centroids_base <- global_centroids_all %>%
  mutate(
    data_scope_chr = as.character(data_scope),
    forest_chr = as.character(forest),
    timespan_chr = as.character(timespan),
    forest_timespan = paste(forest_chr, timespan_chr, sep = ".")
  ) %>%
  select(-any_of("row_panel"))

global_shifts_base <- global_shifts_all %>%
  mutate(
    data_scope_chr = as.character(data_scope),
    forest_chr = as.character(forest)
  ) %>%
  select(-any_of("row_panel"))

global_surface_base <- global_decay_surface_all %>%
  mutate(
    data_scope_chr = as.character(data_scope)
  ) %>%
  select(-any_of("row_panel"))

global_stress_base <- global_stress_all %>%
  mutate(
    data_scope_chr = as.character(data_scope),
    stress_label = paste0("Stress = ", round(stress, 3))
  )

# ---------------------------------------------------------
# 6. Check colour/shape names
# ---------------------------------------------------------

cat("\nForest-time values in scores:\n")
print(sort(unique(global_scores_base$forest_timespan)))

cat("\nNames in forest_time_cols:\n")
print(names(forest_time_cols))

cat("\nTimespan values in scores:\n")
print(sort(unique(global_scores_base$timespan_chr)))

cat("\nNames in timespan_shape_vals:\n")
print(names(timespan_shape_vals))

# ---------------------------------------------------------
# 7. Export function
# ---------------------------------------------------------

make_single_story_panel <- function(
    data_scope_i,
    forest_focus_i,
    panel_title_i,
    show_stress = TRUE
) {
  
  scores_scope <- global_scores_base %>%
    filter(data_scope_chr == data_scope_i)
  
  centroids_scope <- global_centroids_base %>%
    filter(data_scope_chr == data_scope_i)
  
  shifts_scope <- global_shifts_base %>%
    filter(data_scope_chr == data_scope_i)
  
  surface_scope <- global_surface_base %>%
    filter(data_scope_chr == data_scope_i)
  
  stress_scope <- global_stress_base %>%
    filter(data_scope_chr == data_scope_i)
  
  cat("\nPanel:", panel_title_i, "\n")
  cat("  scores rows:   ", nrow(scores_scope), "\n")
  cat("  centroids rows:", nrow(centroids_scope), "\n")
  cat("  shifts rows:   ", nrow(shifts_scope), "\n")
  cat("  surface rows:  ", nrow(surface_scope), "\n")
  
  if (nrow(scores_scope) == 0) {
    stop("No NMDS scores found for data_scope_i = ", data_scope_i)
  }
  
  if (forest_focus_i == "All forests") {
    
    scores_bg <- scores_scope %>%
      slice(0)
    
    scores_fg <- scores_scope
    centroids_fg <- centroids_scope
    shifts_fg <- shifts_scope
    
  } else {
    
    scores_bg <- scores_scope
    
    scores_fg <- scores_scope %>%
      filter(forest_chr == forest_focus_i)
    
    centroids_fg <- centroids_scope %>%
      filter(forest_chr == forest_focus_i)
    
    shifts_fg <- shifts_scope %>%
      filter(forest_chr == forest_focus_i)
  }
  
  stress_label_df <- tibble(
    NMDS1 = -0.96,
    NMDS2 = 0.96,
    stress_label = ifelse(
      nrow(stress_scope) > 0,
      stress_scope$stress_label[1],
      NA_character_
    )
  )
  
  p <- ggplot() +
    geom_contour(
      data = surface_scope,
      aes(
        x = NMDS1,
        y = NMDS2,
        z = log_av_ds_pred
      ),
      inherit.aes = FALSE,
      breaks = decay_breaks,
      colour = "grey50",
      linewidth = 0.13,
      alpha = 0.55
    ) +
    geom_point(
      data = scores_bg,
      aes(
        x = NMDS1,
        y = NMDS2
      ),
      inherit.aes = FALSE,
      colour = "grey83",
      size = 0.48,
      alpha = 0.28,
      stroke = 0
    ) +
    geom_point(
      data = scores_fg,
      aes(
        x = NMDS1,
        y = NMDS2,
        fill = forest_timespan,
        shape = timespan_chr
      ),
      inherit.aes = FALSE,
      colour = "grey22",
      size = 0.90,
      alpha = 0.82,
      stroke = 0.11
    ) +
    geom_segment(
      data = shifts_fg,
      aes(
        x = NMDS1_previous,
        y = NMDS2_previous,
        xend = NMDS1_current,
        yend = NMDS2_current,
        colour = forest_chr
      ),
      inherit.aes = FALSE,
      linewidth = 0.40,
      arrow = arrow(length = unit(0.060, "cm")),
      show.legend = FALSE
    ) +
    geom_point(
      data = centroids_fg,
      aes(
        x = NMDS1,
        y = NMDS2,
        fill = forest_timespan,
        shape = timespan_chr
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 1.75,
      stroke = 0.50,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = forest_time_cols,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = forest_cols,
      guide = "none",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = timespan_shape_vals,
      drop = FALSE
    ) + coord_equal(
      xlim = x_lim_story,
      ylim = y_lim_story,
      clip = "on",
      expand = FALSE
    ) +
    labs(
      title = panel_title_i,
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 7.5) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 7.4,
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(
        colour = "grey35",
        fill = NA,
        linewidth = 0.22
      ),
      legend.position = "none",
      plot.margin = margin(1.5, 1.5, 1.5, 1.5),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA)
    )
  
  if (isTRUE(show_stress) && forest_focus_i == "All forests" && !is.na(stress_label_df$stress_label[1])) {
    p <- p +
      geom_text(
        data = stress_label_df,
        aes(
          x = NMDS1,
          y = NMDS2,
          label = stress_label
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1,
        size = 1.75,
        colour = "grey15"
      )
  }
  
  p
}

# ---------------------------------------------------------
# 8. Export all panels
# ---------------------------------------------------------

exported_single_panels <- single_panel_map %>%
  mutate(
    filename = paste0("nmds_story_panel_", panel_code, ".png"),
    path = file.path(single_panel_dir, filename)
  )

for (i in seq_len(nrow(exported_single_panels))) {
  
  panel_i <- exported_single_panels[i, ]
  
  p_i <- make_single_story_panel(
    data_scope_i = panel_i$data_scope_raw,
    forest_focus_i = panel_i$forest_focus,
    panel_title_i = panel_i$panel_title,
    show_stress = TRUE
  )
  
  ggsave(
    filename = panel_i$path,
    plot = p_i,
    width = single_panel_width_in,
    height = single_panel_height_in,
    dpi = single_panel_dpi,
    bg = "transparent"
  )
  
  print(p_i)
}

cat("\n=========================================================\n")
cat("FIXED SEPARATE NMDS STORY PANELS WRITTEN\n")
cat("=========================================================\n")

print(
  exported_single_panels %>%
    select(panel_code, data_scope_raw, forest_focus, filename, path),
  n = Inf
)