# =========================================================
# TLSR / FLSR / CLSR ANALYSIS
# Proper input sets:
#   1. X_use_full    + meta_use_full
#   2. X_use_subset  + meta_use_subset
#
# Definitions:
#   TLSR = tree-level species richness
#          mean number of taxa per deadwood tree unit
#
#   FLSR = forest-level species richness
#          pooled number of taxa per forest x timespan
#
#   CLSR = country-level species richness
#          pooled number of taxa per country x timespan
# =========================================================

setwd("/data/gent/vo/001/gvo00142/vsc45818/FRUITBODY_PAPER")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
  library(patchwork)
})

# =========================================================
# 0. SETTINGS
# =========================================================

out_dir <- "tlsr_clsr_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

country_levels <- c("Belgium", "Denmark", "Slovenia")

forest_levels <- c(
  "Sonian",
  "Strødam",
  "Suserup",
  "Rajhenav"
)

timespan_levels <- c("previous", "current")

country_cols <- c(
  "Belgium"  = "goldenrod3",
  "Denmark"  = "steelblue3",
  "Slovenia" = "darkgreen"
)

timespan_shape_vals <- c(
  previous = 21,
  current  = 24
)

country_abbrev <- c(
  Belgium  = "BE",
  Denmark  = "DK",
  Slovenia = "SI"
)

standardise_forest_name <- function(x) {
  x <- as.character(x)
  case_when(
    x == "Strodam" ~ "Strødam",
    TRUE ~ x
  )
}

# =========================================================
# 1. HELPER: ANALYSE ONE INPUT SET
# =========================================================

analyse_tlsr_flsr_clsr <- function(X, meta, data_scope) {
  
  stopifnot(is.data.frame(X) || is.matrix(X))
  stopifnot(is.data.frame(meta))
  stopifnot(nrow(X) == nrow(meta))
  stopifnot(all(c("tree_id", "country", "forest", "timespan") %in% names(meta)))
  
  X <- as.data.frame(X)
  taxa_cols <- names(X)
  
  if (!"sample_id" %in% names(meta)) {
    meta <- meta %>%
      mutate(sample_id = rownames(X))
  }
  
  if (is.null(rownames(X)) || any(rownames(X) == "")) {
    rownames(X) <- meta$sample_id
  }
  
  if (!identical(rownames(X), as.character(meta$sample_id))) {
    rownames(X) <- as.character(meta$sample_id)
  }
  
  X <- X %>%
    mutate(across(everything(), ~ as.integer(.x > 0)))
  
  meta_clean <- meta %>%
    mutate(
      data_scope = data_scope,
      country    = factor(as.character(country), levels = country_levels),
      forest     = standardise_forest_name(forest),
      forest     = factor(forest, levels = forest_levels),
      timespan   = factor(as.character(timespan), levels = timespan_levels),
      tree_id    = as.character(tree_id),
      sample_id  = as.character(sample_id)
    )
  
  if ("year" %in% names(meta_clean)) {
    meta_clean <- meta_clean %>%
      mutate(
        year = suppressWarnings(as.integer(as.character(year))),
        survey_year = factor(as.character(year))
      )
  } else {
    meta_clean <- meta_clean %>%
      mutate(
        year = NA_integer_,
        survey_year = factor(NA_character_)
      )
  }
  
  cat("\n==============================\n")
  cat("TLSR/FLSR/CLSR input:", data_scope, "\n")
  cat("==============================\n")
  cat("Samples:", nrow(X), "\n")
  cat("Taxa:", ncol(X), "\n")
  cat("Trees:", n_distinct(meta_clean$tree_id), "\n\n")
  
  cat("Country x timespan table\n")
  print(table(meta_clean$country, meta_clean$timespan, useNA = "ifany"))
  
  cat("\nForest x timespan table\n")
  print(table(meta_clean$forest, meta_clean$timespan, useNA = "ifany"))
  
  cat("\nForest x survey year table\n")
  print(table(meta_clean$forest, meta_clean$survey_year, useNA = "ifany"))
  
  # -------------------------------------------------------
  # Bind metadata and community matrix
  # -------------------------------------------------------
  
  dat_tlsr <- bind_cols(
    meta_clean %>%
      select(
        data_scope,
        sample_id,
        tree_id,
        country,
        forest,
        timespan,
        year,
        survey_year
      ),
    X
  )
  
  # -------------------------------------------------------
  # Collapse to one row per tree within forest x timespan
  # Presence-absence union is used if several sample rows occur per tree.
  # -------------------------------------------------------
  
  tree_level_pa <- dat_tlsr %>%
    filter(
      !is.na(country),
      !is.na(forest),
      !is.na(timespan)
    ) %>%
    group_by(data_scope, country, forest, timespan, tree_id) %>%
    summarise(
      across(all_of(taxa_cols), ~ as.integer(any(.x > 0, na.rm = TRUE))),
      n_sample_rows = n(),
      .groups = "drop"
    )
  
  # -------------------------------------------------------
  # TLSR: tree-level species richness
  # -------------------------------------------------------
  
  tree_richness <- tree_level_pa %>%
    mutate(
      tree_species_richness = rowSums(
        as.data.frame(across(all_of(taxa_cols))),
        na.rm = TRUE
      )
    )
  
  tlsr_by_forest_timespan <- tree_richness %>%
    group_by(data_scope, country, forest, timespan) %>%
    summarise(
      n_trees   = n(),
      TLSR_mean = mean(tree_species_richness, na.rm = TRUE),
      TLSR_sd   = sd(tree_species_richness, na.rm = TRUE),
      TLSR_se   = TLSR_sd / sqrt(n_trees),
      TLSR_min  = min(tree_species_richness, na.rm = TRUE),
      TLSR_max  = max(tree_species_richness, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(data_scope, country, forest, timespan)
  
  # -------------------------------------------------------
  # FLSR: forest-level pooled species richness
  # -------------------------------------------------------
  
  flsr_by_forest_timespan <- tree_level_pa %>%
    group_by(data_scope, country, forest, timespan) %>%
    summarise(
      FLSR = sum(colSums(as.data.frame(across(all_of(taxa_cols))), na.rm = TRUE) > 0),
      .groups = "drop"
    ) %>%
    arrange(data_scope, country, forest, timespan)
  
  tlsr_flsr_by_forest_timespan <- tlsr_by_forest_timespan %>%
    left_join(
      flsr_by_forest_timespan,
      by = c("data_scope", "country", "forest", "timespan")
    ) %>%
    arrange(data_scope, country, forest, timespan)
  
  # -------------------------------------------------------
  # CLSR: country-level pooled species richness
  # -------------------------------------------------------
  
  clsr_by_country_timespan <- tree_level_pa %>%
    group_by(data_scope, country, timespan) %>%
    summarise(
      CLSR = sum(colSums(as.data.frame(across(all_of(taxa_cols))), na.rm = TRUE) > 0),
      .groups = "drop"
    ) %>%
    arrange(data_scope, country, timespan)
  
  # -------------------------------------------------------
  # TLSR + CLSR at country level, retained for comparison
  # -------------------------------------------------------
  
  tlsr_by_country_timespan <- tree_richness %>%
    group_by(data_scope, country, timespan) %>%
    summarise(
      n_trees   = n(),
      TLSR_mean = mean(tree_species_richness, na.rm = TRUE),
      TLSR_sd   = sd(tree_species_richness, na.rm = TRUE),
      TLSR_se   = TLSR_sd / sqrt(n_trees),
      TLSR_min  = min(tree_species_richness, na.rm = TRUE),
      TLSR_max  = max(tree_species_richness, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      clsr_by_country_timespan,
      by = c("data_scope", "country", "timespan")
    ) %>%
    arrange(data_scope, country, timespan)
  
  # -------------------------------------------------------
  # Overall values by timespan only
  # -------------------------------------------------------
  
  tlsr_by_timespan <- tree_richness %>%
    group_by(data_scope, timespan) %>%
    summarise(
      n_trees   = n(),
      TLSR_mean = mean(tree_species_richness, na.rm = TRUE),
      TLSR_sd   = sd(tree_species_richness, na.rm = TRUE),
      TLSR_se   = TLSR_sd / sqrt(n_trees),
      TLSR_min  = min(tree_species_richness, na.rm = TRUE),
      TLSR_max  = max(tree_species_richness, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(data_scope, timespan)
  
  clsr_by_timespan <- tree_level_pa %>%
    group_by(data_scope, timespan) %>%
    summarise(
      CLSR = sum(colSums(as.data.frame(across(all_of(taxa_cols))), na.rm = TRUE) > 0),
      .groups = "drop"
    ) %>%
    arrange(data_scope, timespan)
  
  tlsr_clsr_by_timespan <- tlsr_by_timespan %>%
    left_join(
      clsr_by_timespan,
      by = c("data_scope", "timespan")
    ) %>%
    arrange(data_scope, timespan)
  
  list(
    data_scope = data_scope,
    tree_level_pa = tree_level_pa,
    tree_richness = tree_richness,
    tlsr_by_forest_timespan = tlsr_by_forest_timespan,
    flsr_by_forest_timespan = flsr_by_forest_timespan,
    tlsr_flsr_by_forest_timespan = tlsr_flsr_by_forest_timespan,
    clsr_by_country_timespan = clsr_by_country_timespan,
    tlsr_clsr_by_country_timespan = tlsr_by_country_timespan,
    tlsr_clsr_by_timespan = tlsr_clsr_by_timespan
  )
}

# =========================================================
# 2. RUN BOTH DATASETS
# =========================================================

tlsr_clsr_full <- analyse_tlsr_flsr_clsr(
  X = X_use_full,
  meta = meta_use_full,
  data_scope = "Full target"
)

tlsr_clsr_subset <- analyse_tlsr_flsr_clsr(
  X = X_use_subset,
  meta = meta_use_subset,
  data_scope = "Subsetted"
)

# =========================================================
# 3. COMBINE OUTPUTS
# =========================================================

tree_richness_all_sets <- bind_rows(
  tlsr_clsr_full$tree_richness,
  tlsr_clsr_subset$tree_richness
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    forest     = factor(forest, levels = forest_levels),
    timespan   = factor(timespan, levels = timespan_levels)
  )

tlsr_flsr_by_forest_timespan_all_sets <- bind_rows(
  tlsr_clsr_full$tlsr_flsr_by_forest_timespan,
  tlsr_clsr_subset$tlsr_flsr_by_forest_timespan
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    forest     = factor(forest, levels = forest_levels),
    timespan   = factor(timespan, levels = timespan_levels)
  )

tlsr_clsr_by_country_timespan_all_sets <- bind_rows(
  tlsr_clsr_full$tlsr_clsr_by_country_timespan,
  tlsr_clsr_subset$tlsr_clsr_by_country_timespan
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels)
  )

tlsr_clsr_by_timespan_all_sets <- bind_rows(
  tlsr_clsr_full$tlsr_clsr_by_timespan,
  tlsr_clsr_subset$tlsr_clsr_by_timespan
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    timespan   = factor(timespan, levels = timespan_levels)
  )

cat("\n==============================\n")
cat("TLSR and FLSR by dataset x forest x timespan\n")
cat("==============================\n")
print(tlsr_flsr_by_forest_timespan_all_sets, n = Inf)

cat("\n==============================\n")
cat("TLSR and CLSR by dataset x country x timespan\n")
cat("==============================\n")
print(tlsr_clsr_by_country_timespan_all_sets, n = Inf)

cat("\n==============================\n")
cat("TLSR and CLSR by dataset x timespan\n")
cat("==============================\n")
print(tlsr_clsr_by_timespan_all_sets, n = Inf)

write_csv2(
  tree_richness_all_sets,
  file.path(out_dir, "tree_level_richness_full_vs_subsetted.csv")
)

write_csv2(
  tlsr_flsr_by_forest_timespan_all_sets,
  file.path(out_dir, "tlsr_flsr_by_forest_timespan_full_vs_subsetted.csv")
)

write_csv2(
  tlsr_clsr_by_country_timespan_all_sets,
  file.path(out_dir, "tlsr_clsr_by_country_timespan_full_vs_subsetted.csv")
)

write_csv2(
  tlsr_clsr_by_timespan_all_sets,
  file.path(out_dir, "tlsr_clsr_by_timespan_full_vs_subsetted.csv")
)

# =========================================================
# 4. FOREST-LEVEL TABLE FOR MANUSCRIPT
# =========================================================

forest_table_long <- tlsr_flsr_by_forest_timespan_all_sets %>%
  left_join(
    tlsr_clsr_by_country_timespan_all_sets %>%
      select(data_scope, country, timespan, CLSR),
    by = c("data_scope", "country", "timespan")
  ) %>%
  mutate(
    data_scope_key = case_when(
      as.character(data_scope) == "Subsetted"   ~ "subset",
      as.character(data_scope) == "Full target" ~ "full",
      TRUE ~ as.character(data_scope)
    )
  ) %>%
  select(
    data_scope_key,
    country,
    forest,
    timespan,
    n_trees,
    TLSR_mean,
    TLSR_sd,
    TLSR_se,
    FLSR,
    CLSR
  )

forest_table_paired <- forest_table_long %>%
  pivot_wider(
    names_from = data_scope_key,
    values_from = c(
      n_trees,
      TLSR_mean,
      TLSR_sd,
      TLSR_se,
      FLSR,
      CLSR
    ),
    names_sep = "_"
  ) %>%
  mutate(
    forest_country = paste0(
      as.character(forest),
      " (",
      unname(country_abbrev[as.character(country)]),
      ")"
    ),
    
    n_trees_label = case_when(
      !is.na(n_trees_subset) & !is.na(n_trees_full) ~
        paste0(n_trees_subset, " (", n_trees_full, ")"),
      is.na(n_trees_subset) & !is.na(n_trees_full) ~
        paste0("NA (", n_trees_full, ")"),
      !is.na(n_trees_subset) & is.na(n_trees_full) ~
        paste0(n_trees_subset, " (NA)"),
      TRUE ~ NA_character_
    ),
    
    TLSR_label = case_when(
      !is.na(TLSR_mean_subset) & !is.na(TLSR_mean_full) ~
        paste0(
          sprintf("%.1f", TLSR_mean_subset),
          " ± ",
          sprintf("%.1f", TLSR_sd_subset),
          " (",
          sprintf("%.1f", TLSR_mean_full),
          ")"
        ),
      is.na(TLSR_mean_subset) & !is.na(TLSR_mean_full) ~
        paste0(
          "NA (",
          sprintf("%.1f", TLSR_mean_full),
          ")"
        ),
      !is.na(TLSR_mean_subset) & is.na(TLSR_mean_full) ~
        paste0(
          sprintf("%.1f", TLSR_mean_subset),
          " ± ",
          sprintf("%.1f", TLSR_sd_subset),
          " (NA)"
        ),
      TRUE ~ NA_character_
    ),
    
    TLSR_full_with_sd_label = case_when(
      !is.na(TLSR_mean_subset) & !is.na(TLSR_mean_full) ~
        paste0(
          sprintf("%.1f", TLSR_mean_subset),
          " ± ",
          sprintf("%.1f", TLSR_sd_subset),
          " (",
          sprintf("%.1f", TLSR_mean_full),
          " ± ",
          sprintf("%.1f", TLSR_sd_full),
          ")"
        ),
      TRUE ~ NA_character_
    ),
    
    FLSR_label = case_when(
      !is.na(FLSR_subset) & !is.na(FLSR_full) ~
        paste0(FLSR_subset, " (", FLSR_full, ")"),
      is.na(FLSR_subset) & !is.na(FLSR_full) ~
        paste0("NA (", FLSR_full, ")"),
      !is.na(FLSR_subset) & is.na(FLSR_full) ~
        paste0(FLSR_subset, " (NA)"),
      TRUE ~ NA_character_
    ),
    
    CLSR_label = case_when(
      !is.na(CLSR_subset) & !is.na(CLSR_full) ~
        paste0(CLSR_subset, " (", CLSR_full, ")"),
      is.na(CLSR_subset) & !is.na(CLSR_full) ~
        paste0("NA (", CLSR_full, ")"),
      !is.na(CLSR_subset) & is.na(CLSR_full) ~
        paste0(CLSR_subset, " (NA)"),
      TRUE ~ NA_character_
    )
  )

forest_manuscript_table <- forest_table_paired %>%
  select(
    forest,
    country,
    forest_country,
    timespan,
    n_trees = n_trees_label,
    TLSR = TLSR_label,
    TLSR_full_with_sd = TLSR_full_with_sd_label,
    FLSR = FLSR_label,
    CLSR = CLSR_label
  ) %>%
  mutate(
    forest = factor(as.character(forest), levels = forest_levels),
    timespan = factor(as.character(timespan), levels = timespan_levels)
  ) %>%
  arrange(forest, timespan)

forest_manuscript_table_wide <- forest_manuscript_table %>%
  select(
    forest,
    country,
    forest_country,
    timespan,
    n_trees,
    TLSR,
    FLSR,
    CLSR
  ) %>%
  pivot_wider(
    names_from = timespan,
    values_from = c(n_trees, TLSR, FLSR, CLSR),
    names_glue = "{.value}_{timespan}"
  ) %>%
  mutate(
    forest = factor(as.character(forest), levels = forest_levels)
  ) %>%
  arrange(forest)

cat("\n==============================\n")
cat("FOREST-LEVEL MANUSCRIPT TABLE\n")
cat("Subsetted values are shown first. Full target values are in parentheses.\n")
cat("==============================\n")
print(forest_manuscript_table_wide, n = Inf)

write_csv2(
  forest_manuscript_table,
  file.path(out_dir, "forest_level_richness_table_long_subset_with_full_parentheses.csv")
)

write_csv2(
  forest_manuscript_table_wide,
  file.path(out_dir, "forest_level_richness_table_wide_subset_with_full_parentheses.csv")
)

write_csv2(
  forest_table_paired,
  file.path(out_dir, "forest_level_richness_table_with_raw_paired_values.csv")
)