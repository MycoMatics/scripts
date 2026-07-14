# =========================================================
# FOREST-LEVEL SIMPER × OCCUPANCY-SHIFT ANALYSIS
#
# Purpose:
#   Identify fungal taxa that both:
#     1. contribute strongly to Historical–Current dissimilarity
#        within each forest according to standard SIMPER, and
#     2. show a temporal occupancy shift.
#
# Required objects:
#   X_use_full,    meta_use_full
#   X_use_subset,  meta_use_subset
#
# Main contrast:
#   Historical = 2001
#   Current    = 2021 for Sonian
#   Current    = 2022 for Rajhenav, Strødam, and Suserup
#
# Main figure taxon-selection logic:
#   For the matched subset:
#     1. rank taxa by SIMPER contribution within each forest,
#     2. rank taxa by absolute occupancy shift within each forest,
#     3. start with the intersection of the top 20 taxa from both lists,
#     4. increase the rank threshold stepwise until 20 unique taxa
#        are obtained, or until simper_top_n is reached,
#     5. retain taxa by combined rank importance,
#     6. order taxa in the plot by dominant signed occupancy change.
#
# Main outputs:
#   1. SIMPER table by dataset x forest
#   2. Occupancy-shift table by dataset x forest x taxon
#   3. Joined species-driver table
#   4. Strict threshold-based driver table
#   5. Statistically supported strict-driver table
#   6. Rank-intersection selected taxa and records
#   7. Matched-subset lollipop plot
#   8. SIMPER cumulative diagnostic plot
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(forcats)
  library(vegan)
  library(ggtext)
})

set.seed(42)

# =========================================================
# 0. SETTINGS
# =========================================================

out_dir <- "species_shift_simper_occupancy_outputs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxonomy_file <- "taxonomy_revisedJHC.xlsx"
taxonomy_sheet <- "eco_target"

forest_levels <- c("Sonian", "Rajhenav", "Strødam", "Suserup")
timespan_levels <- c("previous", "current")
dataset_levels <- c("Full data", "Matched subset")

# Strict threshold-based inclusion rule for tables
simper_cumulative_threshold <- 0.66
occupancy_shift_threshold <- 0.05
min_total_occurrence <- 5L

# SIMPER settings
# permutations = 999 adds SIMPER p-values. These are retained as
# complementary evidence, not used as the primary biological filter.
simper_permutations <- 999
simper_top_n <- 50

# Main manuscript figure controls
main_figure_dataset <- "Matched subset"
main_rank_intersection_start_n <- 20L
main_max_taxa_total <- 20L

# Plot settings
plot_dpi <- 900

forest_cols <- c(
  Sonian   = "#B8860B",
  Rajhenav = "#2E7D32",
  Strødam  = "#466C95",
  Suserup  = "#5FA8D3"
)

# =========================================================
# 1. GENERAL HELPERS
# =========================================================

standardize_forest <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_squish(x)
  
  case_when(
    x %in% c("Strødam", "Strodam") ~ "Strødam",
    TRUE ~ x
  )
}

abbreviate_taxon <- function(x, italic = FALSE) {
  x <- as.character(x)
  x <- str_replace_all(x, "_", " ")
  x <- str_squish(x)
  
  out <- map_chr(x, function(z) {
    parts <- str_split(z, "\\s+")[[1]]
    
    if (length(parts) >= 2 && str_detect(parts[1], "^[A-Za-z]+$")) {
      paste0(str_sub(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
    } else {
      z
    }
  })
  
  if (italic) {
    out <- paste0("*", out, "*")
  }
  
  out
}

align_and_prepare_dataset <- function(X, meta, dataset_id, dataset_label) {
  X <- as.data.frame(X)
  meta <- as_tibble(meta)
  
  if (!"sample_id" %in% names(meta)) {
    stop(dataset_id, ": meta must contain sample_id.")
  }
  
  if (!"forest" %in% names(meta)) {
    stop(dataset_id, ": meta must contain forest.")
  }
  
  if (!"timespan" %in% names(meta)) {
    stop(dataset_id, ": meta must contain timespan.")
  }
  
  if (!"year" %in% names(meta)) {
    meta$year <- NA_integer_
  }
  
  meta <- meta %>%
    mutate(
      sample_id = as.character(sample_id),
      forest = standardize_forest(forest),
      forest = factor(forest, levels = forest_levels),
      year = suppressWarnings(as.integer(year)),
      timespan = as.character(timespan),
      timespan = case_when(
        year == 2001L ~ "previous",
        year %in% c(2021L, 2022L) ~ "current",
        timespan %in% c("previous", "current") ~ timespan,
        TRUE ~ NA_character_
      ),
      timespan = factor(timespan, levels = timespan_levels)
    )
  
  if (is.null(rownames(X))) {
    rownames(X) <- meta$sample_id
  }
  
  if (!all(rownames(X) == meta$sample_id)) {
    idx <- match(rownames(X), meta$sample_id)
    
    if (any(is.na(idx))) {
      stop(dataset_id, ": could not align X rows to meta$sample_id.")
    }
    
    meta <- meta[idx, , drop = FALSE]
  }
  
  X[] <- lapply(X, function(z) as.integer(as.numeric(z) > 0))
  
  keep <- !is.na(meta$timespan) &
    !is.na(meta$forest) &
    meta$year != 2023L &
    rowSums(X, na.rm = TRUE) > 0
  
  X <- X[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]
  
  meta <- meta %>%
    mutate(
      dataset = dataset_id,
      dataset_scope = dataset_label
    ) %>%
    droplevels()
  
  stopifnot(nrow(X) == nrow(meta))
  stopifnot(identical(rownames(X), meta$sample_id))
  stopifnot(!any(meta$year == 2023L, na.rm = TRUE))
  
  list(X = X, meta = meta)
}

# =========================================================
# 2. OPTIONAL TAXONOMY JOIN
# =========================================================

read_taxonomy_traits <- function(taxonomy_file, taxonomy_sheet) {
  if (!file.exists(taxonomy_file)) {
    warning(
      "Taxonomy file not found: ", taxonomy_file,
      ". Trait columns will be omitted.",
      call. = FALSE
    )
    return(tibble(taxon = character()))
  }
  
  tax_raw <- readxl::read_excel(taxonomy_file, sheet = taxonomy_sheet) %>%
    as_tibble()
  
  possible_taxon_cols <- c("Taxon", "taxon", "Species", "species", "scientific_name")
  taxon_col <- intersect(possible_taxon_cols, names(tax_raw))[1]
  
  if (is.na(taxon_col)) {
    warning(
      "No taxon-name column found in taxonomy table. Expected one of: ",
      paste(possible_taxon_cols, collapse = ", "),
      call. = FALSE
    )
    return(tibble(taxon = character()))
  }
  
  tax_raw %>%
    rename(taxon = all_of(taxon_col)) %>%
    mutate(
      taxon = as.character(taxon),
      taxon = str_replace_all(taxon, "_", " "),
      taxon = str_squish(taxon)
    ) %>%
    distinct(taxon, .keep_all = TRUE)
}

taxonomy_traits <- read_taxonomy_traits(
  taxonomy_file = taxonomy_file,
  taxonomy_sheet = taxonomy_sheet
)

# =========================================================
# 3. SIMPER AND OCCUPANCY FUNCTIONS
# =========================================================

run_simper_forest <- function(X, meta, dataset_id, dataset_label,
                              top_n = 100,
                              permutations = 999) {
  
  map_dfr(forest_levels, function(fst) {
    idx <- meta$forest == fst
    
    Xf <- X[idx, , drop = FALSE]
    mf <- droplevels(meta[idx, , drop = FALSE])
    
    if (nrow(Xf) < 4 || nlevels(mf$timespan) < 2) {
      return(tibble())
    }
    
    Xf <- Xf[, colSums(Xf, na.rm = TRUE) > 0, drop = FALSE]
    
    if (ncol(Xf) == 0) {
      return(tibble())
    }
    
    sim <- vegan::simper(
      Xf,
      group = mf$timespan,
      permutations = permutations
    )
    
    sm <- summary(sim, ordered = TRUE)
    comp_name <- names(sm)[1]
    
    tb <- as.data.frame(sm[[comp_name]])
    tb$taxon <- rownames(tb)
    rownames(tb) <- NULL
    
    out <- tb %>%
      as_tibble() %>%
      mutate(
        dataset = dataset_id,
        dataset_scope = dataset_label,
        forest = fst,
        comparison = comp_name,
        taxon = as.character(taxon),
        taxon = str_replace_all(taxon, "_", " "),
        taxon = str_squish(taxon),
        simper_rank = row_number(),
        cumulative = cumsum(average) / sum(average, na.rm = TRUE)
      )
    
    if (!"p" %in% names(out)) {
      out$p <- NA_real_
    }
    
    out %>%
      select(
        dataset, dataset_scope, forest, comparison,
        simper_rank, taxon,
        average, sd, ratio, ava, avb,
        cumsum, cumulative, p,
        everything()
      ) %>%
      slice_head(n = top_n)
  })
}

occupancy_shift_forest_all <- function(X, meta, dataset_id, dataset_label) {
  X_long <- X %>%
    rownames_to_column("sample_id") %>%
    pivot_longer(
      cols = -sample_id,
      names_to = "taxon",
      values_to = "pa"
    ) %>%
    mutate(
      taxon = str_replace_all(taxon, "_", " "),
      taxon = str_squish(taxon)
    ) %>%
    left_join(
      meta %>%
        select(sample_id, forest, timespan),
      by = "sample_id"
    ) %>%
    filter(!is.na(forest), !is.na(timespan))
  
  occ_counts <- X_long %>%
    group_by(forest, taxon, timespan) %>%
    summarise(
      n_samples = n(),
      n_occ = sum(pa > 0, na.rm = TRUE),
      occupancy = n_occ / n_samples,
      .groups = "drop"
    )
  
  occ_wide <- occ_counts %>%
    select(forest, taxon, timespan, n_samples, n_occ, occupancy) %>%
    pivot_wider(
      names_from = timespan,
      values_from = c(n_samples, n_occ, occupancy),
      values_fill = list(
        n_samples = 0,
        n_occ = 0,
        occupancy = 0
      )
    )
  
  needed_cols <- c(
    "n_samples_previous", "n_samples_current",
    "n_occ_previous", "n_occ_current",
    "occupancy_previous", "occupancy_current"
  )
  
  for (nm in needed_cols) {
    if (!nm %in% names(occ_wide)) {
      occ_wide[[nm]] <- 0
    }
  }
  
  occ_wide %>%
    mutate(
      dataset = dataset_id,
      dataset_scope = dataset_label,
      previous = occupancy_previous,
      current = occupancy_current,
      delta_occ = current - previous,
      abs_delta = abs(delta_occ),
      total_occ = n_occ_previous + n_occ_current,
      total_samples = n_samples_previous + n_samples_current,
      direction = case_when(
        delta_occ > 0 ~ "increase",
        delta_occ < 0 ~ "decrease",
        TRUE ~ "no_change"
      )
    ) %>%
    select(
      dataset, dataset_scope, forest, taxon,
      n_samples_previous, n_samples_current,
      n_occ_previous, n_occ_current,
      total_occ, total_samples,
      previous, current,
      delta_occ, abs_delta, direction
    )
}

fisher_test_occupancy <- function(occ_df) {
  occ_df %>%
    rowwise() %>%
    mutate(
      fisher_p = {
        mat <- matrix(
          c(
            n_occ_previous,
            n_samples_previous - n_occ_previous,
            n_occ_current,
            n_samples_current - n_occ_current
          ),
          nrow = 2,
          byrow = TRUE
        )
        
        if (any(!is.finite(mat)) || any(mat < 0) || sum(mat) == 0) {
          NA_real_
        } else {
          fisher.test(mat)$p.value
        }
      },
      odds_ratio_current_vs_previous = {
        mat <- matrix(
          c(
            n_occ_current,
            n_samples_current - n_occ_current,
            n_occ_previous,
            n_samples_previous - n_occ_previous
          ),
          nrow = 2,
          byrow = TRUE
        )
        
        if (any(!is.finite(mat)) || any(mat < 0) || sum(mat) == 0) {
          NA_real_
        } else {
          unname(fisher.test(mat)$estimate)
        }
      }
    ) %>%
    ungroup() %>%
    group_by(dataset, forest) %>%
    mutate(
      fisher_fdr = p.adjust(fisher_p, method = "fdr")
    ) %>%
    ungroup()
}

make_driver_table <- function(simper_tbl, occ_tbl, taxonomy_traits) {
  out <- simper_tbl %>%
    left_join(
      occ_tbl,
      by = c("dataset", "dataset_scope", "forest", "taxon")
    ) %>%
    mutate(
      taxon_abbrev = abbreviate_taxon(taxon, italic = FALSE),
      taxon_label = abbreviate_taxon(taxon, italic = TRUE),
      direction = factor(direction, levels = c("decrease", "increase", "no_change")),
      dataset_scope = factor(dataset_scope, levels = dataset_levels),
      forest = factor(forest, levels = forest_levels),
      simper_p_005 = !is.na(p) & p <= 0.05,
      fisher_fdr_005 = !is.na(fisher_fdr) & fisher_fdr <= 0.05,
      stat_supported = simper_p_005 & fisher_fdr_005
    )
  
  if (nrow(taxonomy_traits) > 0) {
    out <- out %>%
      left_join(taxonomy_traits, by = "taxon")
  }
  
  out
}

# =========================================================
# 4. RANK-INTERSECTION SELECTION FOR MAIN PLOT
# =========================================================

select_driver_taxa_by_rank_intersection <- function(
    driver_all,
    dataset_to_plot = "Matched subset",
    start_n = 20L,
    target_n_plot_taxa = 20L,
    max_n = 100L,
    min_total_occurrence = 5L,
    require_directional_shift = TRUE
) {
  
  dat <- driver_all %>%
    filter(dataset_scope == dataset_to_plot) %>%
    filter(!is.na(simper_rank), !is.na(abs_delta)) %>%
    filter(total_occ >= min_total_occurrence)
  
  if (require_directional_shift) {
    dat <- dat %>%
      filter(direction %in% c("increase", "decrease"))
  }
  
  if (nrow(dat) == 0) {
    warning(
      "No taxa available for rank-intersection plot selection after filtering.",
      call. = FALSE
    )
    return(list(
      selected_taxa = tibble(),
      selected_records = tibble(),
      selection_diagnostics = tibble()
    ))
  }
  
  dat <- dat %>%
    group_by(dataset_scope, forest) %>%
    arrange(desc(abs_delta), .by_group = TRUE) %>%
    mutate(
      occupancy_rank = row_number()
    ) %>%
    ungroup()
  
  rank_intersection_at_n <- function(n_use) {
    simper_top <- dat %>%
      group_by(dataset_scope, forest) %>%
      arrange(simper_rank, .by_group = TRUE) %>%
      slice_head(n = n_use) %>%
      ungroup() %>%
      distinct(dataset, dataset_scope, forest, taxon) %>%
      mutate(in_top_simper = TRUE)
    
    occupancy_top <- dat %>%
      group_by(dataset_scope, forest) %>%
      arrange(occupancy_rank, .by_group = TRUE) %>%
      slice_head(n = n_use) %>%
      ungroup() %>%
      distinct(dataset, dataset_scope, forest, taxon) %>%
      mutate(in_top_occupancy = TRUE)
    
    inner_join(
      simper_top,
      occupancy_top,
      by = c("dataset", "dataset_scope", "forest", "taxon")
    ) %>%
      mutate(rank_threshold = n_use)
  }
  
  n_grid <- seq.int(start_n, max_n, by = 1L)
  
  diagnostic <- map_dfr(n_grid, function(n_use) {
    overlap_i <- rank_intersection_at_n(n_use)
    
    tibble(
      rank_threshold = n_use,
      n_overlap_records = nrow(overlap_i),
      n_unique_taxa = n_distinct(overlap_i$taxon)
    )
  })
  
  chosen_n <- diagnostic %>%
    filter(n_unique_taxa >= target_n_plot_taxa) %>%
    summarise(chosen_n = min(rank_threshold), .groups = "drop") %>%
    pull(chosen_n)
  
  if (length(chosen_n) == 0 || is.na(chosen_n)) {
    chosen_n <- max_n
    warning(
      "Target number of unique taxa was not reached by max_n = ",
      max_n,
      ". Using all taxa obtained at max_n.",
      call. = FALSE
    )
  }
  
  overlap_records <- rank_intersection_at_n(chosen_n)
  
  # Select taxa by combined rank importance.
  # This prevents strong decreases from being excluded simply because
  # the final plot is ordered from increases to decreases.
  selected_taxa_importance <- dat %>%
    semi_join(
      overlap_records %>% distinct(taxon),
      by = "taxon"
    ) %>%
    group_by(taxon) %>%
    summarise(
      taxon_abbrev = first(abbreviate_taxon(taxon, italic = FALSE)),
      taxon_label = first(abbreviate_taxon(taxon, italic = TRUE)),
      n_forests_total = n_distinct(forest),
      n_records_total = n(),
      
      best_simper_rank = min(simper_rank, na.rm = TRUE),
      best_occupancy_rank = min(occupancy_rank, na.rm = TRUE),
      best_combined_rank = min(simper_rank + occupancy_rank, na.rm = TRUE),
      mean_combined_rank = mean(simper_rank + occupancy_rank, na.rm = TRUE),
      
      max_abs_delta = max(abs_delta, na.rm = TRUE),
      mean_abs_delta = mean(abs_delta, na.rm = TRUE),
      
      # Dominant signed occupancy change:
      # retains the sign of the forest-level shift with the largest absolute change.
      dominant_signed_delta = delta_occ[which.max(abs_delta)][1],
      
      max_delta_occ = max(delta_occ, na.rm = TRUE),
      min_delta_occ = min(delta_occ, na.rm = TRUE),
      mean_delta_occ = mean(delta_occ, na.rm = TRUE),
      
      mean_simper_average = mean(average, na.rm = TRUE),
      n_stat_supported_records = sum(stat_supported, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(
      best_combined_rank,
      best_simper_rank,
      best_occupancy_rank,
      desc(max_abs_delta),
      desc(n_forests_total),
      taxon_abbrev
    ) %>%
    slice_head(n = target_n_plot_taxa)
  
  # Plot order only.
  # Taxa are ordered from strongest dominant increase to strongest dominant decrease.
  selected_taxa_ranked <- selected_taxa_importance %>%
    arrange(
      desc(dominant_signed_delta),
      best_combined_rank,
      desc(max_abs_delta),
      taxon_abbrev
    ) %>%
    mutate(taxon_rank = row_number())
  
  selected_records <- dat %>%
    semi_join(selected_taxa_ranked %>% distinct(taxon), by = "taxon") %>%
    select(-any_of(c("taxon_abbrev", "taxon_label"))) %>%
    left_join(
      selected_taxa_ranked %>%
        select(
          taxon, taxon_abbrev, taxon_label, taxon_rank,
          n_forests_total, n_records_total,
          best_simper_rank, best_occupancy_rank,
          best_combined_rank, dominant_signed_delta,
          max_abs_delta, mean_abs_delta
        ),
      by = "taxon"
    ) %>%
    mutate(
      selection_rank_threshold = chosen_n,
      selected_by_rank_intersection = TRUE
    )
  
  stopifnot("taxon_label" %in% names(selected_records))
  stopifnot("taxon_rank" %in% names(selected_records))
  
  list(
    selected_taxa = selected_taxa_ranked,
    selected_records = selected_records,
    selection_diagnostics = diagnostic %>%
      mutate(
        chosen = rank_threshold == chosen_n,
        target_n_plot_taxa = target_n_plot_taxa,
        start_n = start_n,
        max_n = max_n
      )
  )
}

# =========================================================
# 5. RUN ANALYSIS FOR FULL AND MATCHED SUBSET
# =========================================================

full_prepared <- align_and_prepare_dataset(
  X = X_use_full,
  meta = meta_use_full,
  dataset_id = "full",
  dataset_label = "Full data"
)

subset_prepared <- align_and_prepare_dataset(
  X = X_use_subset,
  meta = meta_use_subset,
  dataset_id = "subset",
  dataset_label = "Matched subset"
)

X_full <- full_prepared$X
meta_full <- full_prepared$meta

X_subset <- subset_prepared$X
meta_subset <- subset_prepared$meta

cat("\n====================\n")
cat("FULL DATA CHECK\n")
cat("====================\n")
cat("Samples:", nrow(X_full), "\n")
cat("Taxa:", ncol(X_full), "\n")
print(table(meta_full$forest, meta_full$timespan, useNA = "ifany"))

cat("\n====================\n")
cat("MATCHED SUBSET CHECK\n")
cat("====================\n")
cat("Samples:", nrow(X_subset), "\n")
cat("Taxa:", ncol(X_subset), "\n")
print(table(meta_subset$forest, meta_subset$timespan, useNA = "ifany"))

simper_full <- run_simper_forest(
  X = X_full,
  meta = meta_full,
  dataset_id = "full",
  dataset_label = "Full data",
  top_n = simper_top_n,
  permutations = simper_permutations
)

simper_subset <- run_simper_forest(
  X = X_subset,
  meta = meta_subset,
  dataset_id = "subset",
  dataset_label = "Matched subset",
  top_n = simper_top_n,
  permutations = simper_permutations
)

simper_all <- bind_rows(simper_full, simper_subset)

occ_full <- occupancy_shift_forest_all(
  X = X_full,
  meta = meta_full,
  dataset_id = "full",
  dataset_label = "Full data"
)

occ_subset <- occupancy_shift_forest_all(
  X = X_subset,
  meta = meta_subset,
  dataset_id = "subset",
  dataset_label = "Matched subset"
)

occ_all <- bind_rows(occ_full, occ_subset) %>%
  fisher_test_occupancy()

driver_all <- make_driver_table(
  simper_tbl = simper_all,
  occ_tbl = occ_all,
  taxonomy_traits = taxonomy_traits
)

# =========================================================
# 6. FILTER STRICT THRESHOLD-BASED DRIVER TABLE
# =========================================================

driver_strict <- driver_all %>%
  filter(
    cumulative <= simper_cumulative_threshold,
    abs_delta >= occupancy_shift_threshold,
    total_occ >= min_total_occurrence,
    direction %in% c("increase", "decrease")
  ) %>%
  arrange(dataset_scope, forest, desc(abs_delta), simper_rank)

driver_strict_stat_supported <- driver_strict %>%
  filter(stat_supported) %>%
  arrange(dataset_scope, forest, desc(abs_delta), simper_rank)

driver_presence <- driver_strict %>%
  distinct(taxon, forest, dataset_scope) %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = dataset_scope,
    values_from = present,
    values_fill = FALSE
  )

if (!"Full data" %in% names(driver_presence)) {
  driver_presence[["Full data"]] <- FALSE
}

if (!"Matched subset" %in% names(driver_presence)) {
  driver_presence[["Matched subset"]] <- FALSE
}

driver_presence <- driver_presence %>%
  mutate(
    driver_presence = case_when(
      `Full data` & `Matched subset` ~ "both",
      `Full data` & !`Matched subset` ~ "full_only",
      !`Full data` & `Matched subset` ~ "subset_only",
      TRUE ~ "not_selected"
    )
  ) %>%
  select(taxon, forest, driver_presence)

driver_strict <- driver_strict %>%
  left_join(driver_presence, by = c("taxon", "forest")) %>%
  arrange(dataset_scope, forest, delta_occ)

driver_strict_stat_supported <- driver_strict_stat_supported %>%
  left_join(driver_presence, by = c("taxon", "forest")) %>%
  arrange(dataset_scope, forest, delta_occ)

cat("\n====================\n")
cat("STRICT DRIVER TABLE STRUCTURE\n")
cat("====================\n")
str(driver_strict)

cat("\n====================\n")
cat("STRICT STATISTICALLY SUPPORTED DRIVER TABLE STRUCTURE\n")
cat("====================\n")
str(driver_strict_stat_supported)

# =========================================================
# 7. SUMMARY TABLES
# =========================================================

driver_summary_forest <- driver_strict %>%
  group_by(dataset_scope, forest, direction) %>%
  summarise(
    n_driver_taxa = n_distinct(taxon),
    mean_abs_delta = mean(abs_delta, na.rm = TRUE),
    mean_simper_average = mean(average, na.rm = TRUE),
    n_simper_p_005 = sum(simper_p_005, na.rm = TRUE),
    n_fisher_fdr_005 = sum(fisher_fdr_005, na.rm = TRUE),
    n_stat_supported = sum(stat_supported, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dataset_scope, forest, direction)

driver_summary_presence <- driver_strict %>%
  distinct(forest, taxon, driver_presence) %>%
  count(forest, driver_presence, name = "n_taxa") %>%
  arrange(forest, driver_presence)

driver_overlap_table <- driver_strict %>%
  mutate(
    dataset_scope = factor(dataset_scope, levels = dataset_levels),
    forest = factor(forest, levels = forest_levels)
  ) %>%
  group_by(dataset_scope, taxon, taxon_abbrev, taxon_label) %>%
  summarise(
    n_forests = n_distinct(forest),
    forests = paste(sort(unique(as.character(forest))), collapse = " | "),
    n_increases = sum(direction == "increase", na.rm = TRUE),
    n_decreases = sum(direction == "decrease", na.rm = TRUE),
    max_delta_occ = max(delta_occ, na.rm = TRUE),
    min_delta_occ = min(delta_occ, na.rm = TRUE),
    mean_delta_occ = mean(delta_occ, na.rm = TRUE),
    max_abs_delta = max(abs_delta, na.rm = TRUE),
    mean_abs_delta = mean(abs_delta, na.rm = TRUE),
    mean_simper_average = mean(average, na.rm = TRUE),
    min_simper_cumulative = min(cumulative, na.rm = TRUE),
    min_simper_p = min(p, na.rm = TRUE),
    min_fisher_p = min(fisher_p, na.rm = TRUE),
    min_fisher_fdr = min(fisher_fdr, na.rm = TRUE),
    n_stat_supported_records = sum(stat_supported, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    min_simper_p = ifelse(is.infinite(min_simper_p), NA_real_, min_simper_p),
    min_fisher_p = ifelse(is.infinite(min_fisher_p), NA_real_, min_fisher_p),
    min_fisher_fdr = ifelse(is.infinite(min_fisher_fdr), NA_real_, min_fisher_fdr)
  ) %>%
  arrange(dataset_scope, desc(n_forests), desc(max_abs_delta), taxon)
# =========================================================
# 8. MAIN AND SUPPLEMENTARY PLOT SELECTION BY RANK INTERSECTION
# =========================================================

prepare_rank_intersection_plot_data <- function(
    driver_all,
    dataset_to_plot,
    start_n = 20L,
    target_n_plot_taxa = 20L,
    max_n = 50L,
    min_total_occurrence = 5L,
    require_directional_shift = TRUE
) {
  
  rank_selection <- select_driver_taxa_by_rank_intersection(
    driver_all = driver_all,
    dataset_to_plot = dataset_to_plot,
    start_n = start_n,
    target_n_plot_taxa = target_n_plot_taxa,
    max_n = max_n,
    min_total_occurrence = min_total_occurrence,
    require_directional_shift = require_directional_shift
  )
  
  plot_df <- rank_selection$selected_records
  
  forest_offsets <- c(
    "Sonian"   = -0.27,
    "Rajhenav" = -0.09,
    "Strødam"  =  0.09,
    "Suserup"  =  0.27
  )
  
  if (nrow(plot_df) > 0) {
    plot_df <- plot_df %>%
      mutate(
        dataset_scope = factor(dataset_scope, levels = dataset_levels),
        forest = factor(forest, levels = forest_levels),
        y_base = taxon_rank,
        y_plot = y_base + unname(forest_offsets[as.character(forest)])
      )
  }
  
  list(
    selection = rank_selection,
    plot_df = plot_df
  )
}

matched_plot_selection <- prepare_rank_intersection_plot_data(
  driver_all = driver_all,
  dataset_to_plot = "Matched subset",
  start_n = main_rank_intersection_start_n,
  target_n_plot_taxa = main_max_taxa_total,
  max_n = simper_top_n,
  min_total_occurrence = min_total_occurrence,
  require_directional_shift = TRUE
)

full_plot_selection <- prepare_rank_intersection_plot_data(
  driver_all = driver_all,
  dataset_to_plot = "Full data",
  start_n = main_rank_intersection_start_n,
  target_n_plot_taxa = main_max_taxa_total,
  max_n = simper_top_n,
  min_total_occurrence = min_total_occurrence,
  require_directional_shift = TRUE
)

driver_main_lollipop_df <- matched_plot_selection$plot_df
driver_full_lollipop_df <- full_plot_selection$plot_df

cat("\n====================\n")
cat("RANK-INTERSECTION PLOT SELECTION: MATCHED SUBSET\n")
cat("====================\n")
print(
  matched_plot_selection$selection$selection_diagnostics %>%
    filter(chosen),
  n = Inf
)

cat("\nSelected taxa, matched subset:\n")
print(matched_plot_selection$selection$selected_taxa, n = Inf)

cat("\n====================\n")
cat("RANK-INTERSECTION PLOT SELECTION: FULL DATA\n")
cat("====================\n")
print(
  full_plot_selection$selection$selection_diagnostics %>%
    filter(chosen),
  n = Inf
)

cat("\nSelected taxa, full data:\n")
print(full_plot_selection$selection$selected_taxa, n = Inf)

cat("\nEutypa diagnostic in matched subset:\n")
print(
  driver_all %>%
    filter(
      dataset_scope == "Matched subset",
      str_detect(taxon, regex("Eutypa|spinosa", ignore_case = TRUE))
    ) %>%
    select(
      dataset_scope, forest, taxon,
      simper_rank, average, cumulative, p,
      previous, current, delta_occ, abs_delta,
      total_occ, fisher_p, fisher_fdr
    ) %>%
    arrange(forest, simper_rank),
  n = Inf
)

cat("\nEutypa diagnostic in full data:\n")
print(
  driver_all %>%
    filter(
      dataset_scope == "Full data",
      str_detect(taxon, regex("Eutypa|spinosa", ignore_case = TRUE))
    ) %>%
    select(
      dataset_scope, forest, taxon,
      simper_rank, average, cumulative, p,
      previous, current, delta_occ, abs_delta,
      total_occ, fisher_p, fisher_fdr
    ) %>%
    arrange(forest, simper_rank),
  n = Inf
)

cat("\nEutypa diagnostic in selected matched-subset plot taxa:\n")
print(
  matched_plot_selection$selection$selected_taxa %>%
    filter(str_detect(taxon, regex("Eutypa|spinosa", ignore_case = TRUE))) %>%
    select(
      taxon,
      best_simper_rank,
      best_occupancy_rank,
      best_combined_rank,
      dominant_signed_delta,
      max_abs_delta,
      taxon_rank
    ),
  n = Inf
)

cat("\nEutypa diagnostic in selected full-data plot taxa:\n")
print(
  full_plot_selection$selection$selected_taxa %>%
    filter(str_detect(taxon, regex("Eutypa|spinosa", ignore_case = TRUE))) %>%
    select(
      taxon,
      best_simper_rank,
      best_occupancy_rank,
      best_combined_rank,
      dominant_signed_delta,
      max_abs_delta,
      taxon_rank
    ),
  n = Inf
)

# =========================================================
# 9. WRITE TABLES
# =========================================================

write_csv2(
  simper_all,
  file.path(out_dir, "simper_forest_full_vs_subset.csv")
)

write_csv2(
  occ_all,
  file.path(out_dir, "occupancy_shift_forest_all_taxa_full_vs_subset.csv")
)

write_csv2(
  driver_all,
  file.path(out_dir, "species_driver_all_simper_joined_occupancy_full_vs_subset.csv")
)

write_csv2(
  driver_strict,
  file.path(out_dir, "species_driver_STRICT_full_vs_subset.csv")
)

write_csv2(
  driver_strict_stat_supported,
  file.path(out_dir, "species_driver_STRICT_statistically_supported_full_vs_subset.csv")
)

write_csv2(
  driver_summary_forest,
  file.path(out_dir, "species_driver_summary_by_forest.csv")
)

write_csv2(
  driver_summary_presence,
  file.path(out_dir, "species_driver_presence_full_subset_overlap.csv")
)

write_csv2(
  driver_overlap_table,
  file.path(out_dir, "species_driver_overlap_by_taxon_full_vs_subset.csv")
)

write_csv2(
  driver_overlap_table %>% filter(n_forests >= 2),
  file.path(out_dir, "species_driver_shared_across_forests_full_vs_subset.csv")
)

write_csv2(
  matched_plot_selection$selection$selected_taxa,
  file.path(out_dir, "plot_selection_rank_intersection_selected_taxa_matched_subset.csv")
)

write_csv2(
  matched_plot_selection$selection$selected_records,
  file.path(out_dir, "plot_selection_rank_intersection_selected_records_matched_subset.csv")
)

write_csv2(
  matched_plot_selection$selection$selection_diagnostics,
  file.path(out_dir, "plot_selection_rank_intersection_diagnostics_matched_subset.csv")
)

write_csv2(
  full_plot_selection$selection$selected_taxa,
  file.path(out_dir, "plot_selection_rank_intersection_selected_taxa_full_data.csv")
)

write_csv2(
  full_plot_selection$selection$selected_records,
  file.path(out_dir, "plot_selection_rank_intersection_selected_records_full_data.csv")
)

write_csv2(
  full_plot_selection$selection$selection_diagnostics,
  file.path(out_dir, "plot_selection_rank_intersection_diagnostics_full_data.csv")
)

# =========================================================
# 10. LOLLIPOP PLOTS: MATCHED SUBSET AND FULL DATA
# =========================================================

make_driver_lollipop_plot <- function(
    plot_df,
    plot_title,
    plot_subtitle = NULL,
    forest_cols,
    dataset_label_for_warning = "dataset"
) {
  
  if (nrow(plot_df) == 0) {
    warning(
      "No rank-intersection taxa found for ",
      dataset_label_for_warning,
      " lollipop figure.",
      call. = FALSE
    )
    return(NULL)
  }
  
  taxon_breaks <- plot_df %>%
    distinct(taxon_rank, taxon_label) %>%
    arrange(taxon_rank)
  
  taxon_bands <- taxon_breaks %>%
    mutate(
      ymin = taxon_rank - 0.5,
      ymax = taxon_rank + 0.5,
      band = taxon_rank %% 2
    ) %>%
    filter(band == 0)
  
  xlim_min <- min(plot_df$delta_occ, na.rm = TRUE) - 0.05
  xlim_max <- max(plot_df$delta_occ, na.rm = TRUE) + 0.05
  
  x_abs <- max(abs(c(xlim_min, xlim_max)), na.rm = TRUE)
  x_limits <- c(-x_abs, x_abs)
  
  ggplot(plot_df) +
    geom_rect(
      data = taxon_bands,
      aes(
        ymin = ymin,
        ymax = ymax,
        xmin = -Inf,
        xmax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey94",
      colour = NA
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.35,
      colour = "grey25"
    ) +
    geom_segment(
      aes(
        x = 0,
        xend = delta_occ,
        y = y_plot,
        yend = y_plot,
        colour = forest
      ),
      linewidth = 0.30,
      alpha = 0.60,
      lineend = "round"
    ) +
    geom_point(
      aes(
        x = delta_occ,
        y = y_plot,
        colour = forest,
        size = average
      ),
      alpha = 0.80
    ) +
    scale_colour_manual(
      values = forest_cols,
      name = "Forest"
    ) +
    scale_size_continuous(
      range = c(1.6, 4.8),
      trans = "sqrt",
      name = "SIMPER average"
    ) +
    scale_y_continuous(
      breaks = taxon_breaks$taxon_rank,
      labels = taxon_breaks$taxon_label,
      trans = "reverse",
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = c(max(taxon_breaks$taxon_rank) + 0.6, 0.4),
      clip = "off"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.18),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.y = ggtext::element_markdown(size = 7.8),
      axis.text.x = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "left",
      legend.justification = "left",
      legend.margin = margin(t = 2, r = 0, b = 2, l = 0),
      legend.box.margin = margin(t = 2, r = 0, b = 2, l = -6),
      legend.title = element_text(size = 8.5),
      legend.text = element_text(size = 8),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8.5),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 3, alpha = 1),
        order = 1,
        nrow = 1,
        byrow = TRUE
      ),
      size = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE
      )
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Change in occupancy from Historical survey",
      y = NULL
    )
}

# Matched subset figure for main manuscript
p_driver_main_matched <- make_driver_lollipop_plot(
  plot_df = driver_main_lollipop_df,
  plot_title = "Temporal taxon change in the matched subset",
  plot_subtitle = NULL,
  forest_cols = forest_cols,
  dataset_label_for_warning = "matched subset"
)

if (!is.null(p_driver_main_matched)) {
  
  print(p_driver_main_matched)
  
  main_lollipop_height_matched <- max(
    4.8,
    0.20 * n_distinct(driver_main_lollipop_df$taxon_label) + 1.8
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      "species_drivers_MAIN_matched_subset_rank_intersection_lollipop.png"
    ),
    plot = p_driver_main_matched,
    width = 7.4 / 1.1,
    height = main_lollipop_height_matched,
    dpi = plot_dpi,
    bg = "transparent"
  )
}

# Full dataset figure for supplementary material
p_driver_main_full <- make_driver_lollipop_plot(
  plot_df = driver_full_lollipop_df,
  plot_title = "Temporal taxon change in the full dataset",
  plot_subtitle = "Supplementary figure based on the full tree-level presence-absence matrix",
  forest_cols = forest_cols,
  dataset_label_for_warning = "full data"
)

if (!is.null(p_driver_main_full)) {
  
  print(p_driver_main_full)
  
  main_lollipop_height_full <- max(
    4.8,
    0.20 * n_distinct(driver_full_lollipop_df$taxon_label) + 1.8
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      "species_drivers_SUPPLEMENT_full_data_rank_intersection_lollipop.png"
    ),
    plot = p_driver_main_full,
    width = 7.4 / 1.1,
    height = main_lollipop_height_full,
    dpi = plot_dpi,
    bg = "transparent"
  )
}

