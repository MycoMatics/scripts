# =========================================================
# HARMONISED PERMANOVA ANALYSIS
# Full dataset + substrate-matched subset
#
# Core models retained:
#   1. Global unadjusted model:
#        D ~ forest + timespan + forest:timespan
#   2. Global substrate-adjusted model:
#        D ~ forest + ENV + timespan + forest:timespan
#   3. Forest-level post hoc models:
#        D_forest ~ ENV + timespan
#
# Permutations:
#   All PERMANOVA models use unrestricted permutations.
#   Forest-level models are run separately within each forest subset.
#
# Diagnostics:
#   4. Beta-dispersion tests
#   5. Forest-level beta-dispersion tests
#   6. Distance-to-centroid summaries
#   7. Pairwise beta-dispersion tests for forest x timespan
#
# Tables:
#   8. Raw PERMANOVA outputs
#   9. PERMANOVA outputs with matched dispersion diagnostics
#  10. Main manuscript table
#  11. Supplementary PERMANOVA table with matched dispersion
#  12. Supplementary standalone beta-dispersion table
#
# Required input objects:
#   X_use_full, meta_use_full
#   X_use_subset, meta_use_subset
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(vegan)
  library(permute)
  library(ggplot2)
  library(forcats)
  library(scales)
})

set.seed(42)

# =========================================================
# 0. SETTINGS
# =========================================================

nperm <- 9999

outdir <- "permanova_outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

country_levels  <- c("Belgium", "Denmark", "Slovenia")
forest_levels   <- c("Sonian", "Rajhenav", "Strødam", "Suserup")
timespan_levels <- c("previous", "current")

candidate_env_vars <- c(
  "log_av_ds",
  "log_diameter",
  "exposure",
  "mosscover",
  "snag"
)

dataset_order <- c("Full", "Substrate-matched subset")

model_level_order <- c(
  "Global",
  "Forest-level"
)

model_order <- c(
  "forest + period + forest × period",
  "forest + substrate + period + forest × period",
  "forest: Sonian",
  "forest: Rajhenav",
  "forest: Strødam",
  "forest: Suserup"
)

term_order_main <- c(
  "Survey period",
  "Forest × period"
)

term_order_supp <- c(
  "Forest",
  "Decay stage",
  "Log diameter",
  "Exposure",
  "Moss cover",
  "Snag",
  "Survey period",
  "Forest × period",
  "Residual",
  "Total"
)

# =========================================================
# 1. HELPER FUNCTIONS
# =========================================================

standardise_forest_names <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "Strodam", "Strødam")
  x
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p <= 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

format_num <- function(x, digits = 4) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(x, format = "f", digits = digits)
  )
}

clean_dataset_label <- function(x) {
  recode(
    x,
    "full" = "Full",
    "subset" = "Substrate-matched subset",
    .default = x
  )
}

clean_permanova_term <- function(x) {
  recode(
    x,
    "timespan" = "Survey period",
    "forest:timespan" = "Forest × period",
    "forest" = "Forest",
    "log_av_ds" = "Decay stage",
    "log_diameter" = "Log diameter",
    "exposure" = "Exposure",
    "mosscover" = "Moss cover",
    "snag" = "Snag",
    "Residual" = "Residual",
    "Total" = "Total",
    .default = x
  )
}

clean_dispersion_grouping <- function(x) {
  recode(
    x,
    "timespan" = "Survey period",
    "forest_timespan" = "Forest × period",
    "country_timespan" = "Country × period",
    "forest:timespan" = "Forest × period",
    .default = x
  )
}

clean_pa_inputs <- function(X, meta) {
  stopifnot(nrow(X) == nrow(meta))
  
  X <- as.data.frame(X)
  meta <- as.data.frame(meta)
  
  X[] <- lapply(X, function(z) as.integer(as.numeric(z) > 0))
  
  keep_samples <- rowSums(X, na.rm = TRUE) > 0
  X <- X[keep_samples, , drop = FALSE]
  meta <- meta[keep_samples, , drop = FALSE]
  
  keep_taxa <- colSums(X, na.rm = TRUE) > 0
  X <- X[, keep_taxa, drop = FALSE]
  
  meta <- meta %>%
    mutate(
      country = factor(country, levels = country_levels),
      forest = factor(
        standardise_forest_names(forest),
        levels = forest_levels
      ),
      timespan = factor(timespan, levels = timespan_levels)
    )
  
  if ("log_av_ds" %in% names(meta)) {
    meta$log_av_ds <- as.numeric(meta$log_av_ds)
  }
  
  if ("log_diameter" %in% names(meta)) {
    meta$log_diameter <- as.numeric(meta$log_diameter)
  }
  
  if ("exposure" %in% names(meta)) {
    meta$exposure <- factor(meta$exposure)
  }
  
  if ("mosscover" %in% names(meta)) {
    meta$mosscover <- as.numeric(meta$mosscover)
  }
  
  if ("snag" %in% names(meta)) {
    meta$snag <- factor(meta$snag)
  }
  
  list(X = X, meta = meta)
}

make_jaccard <- function(X) {
  vegan::vegdist(X, method = "jaccard", binary = TRUE)
}

make_perm_free <- function(nperm = 9999) {
  how(nperm = nperm)
}

make_formula <- function(response, rhs_terms) {
  as.formula(
    paste(response, "~", paste(rhs_terms, collapse = " + "))
  )
}

tidy_adonis <- function(x, dataset, model, permutation_scope, by) {
  as.data.frame(x) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      dataset = dataset,
      model = model,
      permutation_scope = permutation_scope,
      by = by,
      .before = 1
    )
}

run_adonis <- function(D, meta, rhs_terms, permutations, by = "terms") {
  form <- make_formula("D", rhs_terms)
  
  vegan::adonis2(
    formula = form,
    data = meta,
    permutations = permutations,
    by = by
  )
}

get_usable_env_vars <- function(meta, candidate_env_vars) {
  env_vars <- intersect(candidate_env_vars, names(meta))
  
  env_vars[
    vapply(
      env_vars,
      function(v) length(unique(na.omit(meta[[v]]))) > 1,
      logical(1)
    )
  ]
}

print_dataset_summary <- function(X, meta, dataset_name, env_vars) {
  cat("\n=========================================================\n")
  cat("DATASET:", dataset_name, "\n")
  cat("=========================================================\n")
  cat("Samples:", nrow(X), "\n")
  cat("Taxa:", ncol(X), "\n")
  cat("Environmental terms used:", paste(env_vars, collapse = ", "), "\n\n")
  
  cat("Country x timespan\n")
  print(table(meta$country, meta$timespan))
  
  cat("\nForest x timespan\n")
  print(table(meta$forest, meta$timespan))
  
  cat("\nForest x survey year, if available\n")
  if ("survey_year" %in% names(meta)) {
    print(table(meta$forest, meta$survey_year))
  } else if ("year" %in% names(meta)) {
    print(table(meta$forest, meta$year))
  } else {
    cat("No survey_year or year column found.\n")
  }
  
  cat("\nEnvironmental summaries\n")
  if (length(env_vars) > 0) {
    print(summary(meta[, env_vars, drop = FALSE]))
  } else {
    cat("No environmental variables available.\n")
  }
}

# =========================================================
# 2. PERMANOVA MODEL SET
# =========================================================

run_model_set <- function(X, meta, dataset_name, env_vars = candidate_env_vars, nperm = 9999) {
  dat <- clean_pa_inputs(X, meta)
  X <- dat$X
  meta <- dat$meta
  
  env_vars <- get_usable_env_vars(meta, env_vars)
  
  print_dataset_summary(
    X = X,
    meta = meta,
    dataset_name = dataset_name,
    env_vars = env_vars
  )
  
  D <- make_jaccard(X)
  perm_free <- make_perm_free(nperm)
  
  results <- list()
  
  # -------------------------------------------------------
  # M1. Global unadjusted model with unrestricted permutations
  # -------------------------------------------------------
  
  m1_terms <- c("forest", "timespan", "forest:timespan")
  
  m1 <- run_adonis(
    D = D,
    meta = meta,
    rhs_terms = m1_terms,
    permutations = perm_free,
    by = "terms"
  )
  
  results[["M1_global_unadjusted_free_terms"]] <- tidy_adonis(
    m1,
    dataset = dataset_name,
    model = "D ~ forest + timespan + forest:timespan",
    permutation_scope = "free",
    by = "terms"
  )
  
  cat("\nM1: global unadjusted model, unrestricted permutations, by = terms\n")
  print(m1)
  
  # -------------------------------------------------------
  # M2. Global substrate-adjusted model with unrestricted permutations
  # -------------------------------------------------------
  
  if (length(env_vars) > 0) {
    m2_terms <- c(
      "forest",
      env_vars,
      "timespan",
      "forest:timespan"
    )
    
    m2_model_text <- paste(
      "D ~",
      paste(m2_terms, collapse = " + ")
    )
    
    m2 <- run_adonis(
      D = D,
      meta = meta,
      rhs_terms = m2_terms,
      permutations = perm_free,
      by = "terms"
    )
    
    results[["M2_main_adjusted_free_terms"]] <- tidy_adonis(
      m2,
      dataset = dataset_name,
      model = m2_model_text,
      permutation_scope = "free",
      by = "terms"
    )
    
    cat("\nM2: global substrate-adjusted model, unrestricted permutations, by = terms\n")
    print(m2)
  }
  
  # -------------------------------------------------------
  # M3. Forest-level post hoc models
  # -------------------------------------------------------
  
  forest_results <- list()
  
  for (ff in levels(droplevels(meta$forest))) {
    idx <- meta$forest == ff
    
    X_f <- X[idx, , drop = FALSE]
    meta_f <- droplevels(meta[idx, , drop = FALSE])
    
    keep_taxa_f <- colSums(X_f, na.rm = TRUE) > 0
    X_f <- X_f[, keep_taxa_f, drop = FALSE]
    
    if (nrow(X_f) < 6 || length(unique(meta_f$timespan)) < 2) {
      next
    }
    
    D_f <- make_jaccard(X_f)
    env_f <- get_usable_env_vars(meta_f, env_vars)
    
    m3_terms <- c(env_f, "timespan")
    
    m3_model_text <- paste0(
      "Forest-level model for ",
      ff,
      ": D ~ ",
      paste(m3_terms, collapse = " + ")
    )
    
    m3 <- run_adonis(
      D = D_f,
      meta = meta_f,
      rhs_terms = m3_terms,
      permutations = make_perm_free(nperm),
      by = "terms"
    )
    
    forest_results[[ff]] <- tidy_adonis(
      m3,
      dataset = dataset_name,
      model = m3_model_text,
      permutation_scope = "free within forest subset",
      by = "terms"
    ) %>%
      mutate(
        forest_posthoc = ff,
        .after = dataset
      )
    
    cat("\nM3: forest-level post hoc model for", ff, "\n")
    print(m3)
  }
  
  if (length(forest_results) > 0) {
    results[["M3_forest_level_posthoc"]] <- bind_rows(forest_results)
  }
  
  bind_rows(results)
}

# =========================================================
# 3. RUN PERMANOVA MODELS
# =========================================================

permanova_full <- run_model_set(
  X = X_use_full,
  meta = meta_use_full,
  dataset_name = "full",
  env_vars = candidate_env_vars,
  nperm = nperm
)

permanova_subset <- run_model_set(
  X = X_use_subset,
  meta = meta_use_subset,
  dataset_name = "subset",
  env_vars = candidate_env_vars,
  nperm = nperm
)

permanova_all <- bind_rows(
  permanova_full,
  permanova_subset
)

cat("\n=========================================================\n")
cat("COMBINED PERMANOVA RESULTS\n")
cat("=========================================================\n")
print(permanova_all, n = Inf)

write_csv2(
  permanova_all,
  file.path(outdir, "permanova_model_set_results.csv")
)

# =========================================================
# 4. EXTRACT KEY PERMANOVA SUMMARIES
# =========================================================

key_terms <- c(
  "forest",
  "timespan",
  "forest:timespan",
  "log_av_ds",
  "log_diameter",
  "exposure",
  "mosscover",
  "snag"
)

permanova_key <- permanova_all %>%
  filter(term %in% key_terms) %>%
  arrange(dataset, model, by, term)

cat("\n=========================================================\n")
cat("KEY PERMANOVA TERMS\n")
cat("=========================================================\n")
print(permanova_key, n = Inf)

write_csv2(
  permanova_key,
  file.path(outdir, "permanova_key_terms.csv")
)

main_temporal_terms <- permanova_all %>%
  filter(
    is.na(forest_posthoc),
    term %in% c("timespan", "forest:timespan")
  ) %>%
  select(
    dataset,
    model,
    permutation_scope,
    by,
    term,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  arrange(dataset, model, term)

cat("\n=========================================================\n")
cat("MAIN TEMPORAL AND INTERACTION TERMS\n")
cat("=========================================================\n")
print(main_temporal_terms, n = Inf)

write_csv2(
  main_temporal_terms,
  file.path(outdir, "main_temporal_and_interaction_terms.csv")
)

substrate_terms <- permanova_all %>%
  filter(
    is.na(forest_posthoc),
    term %in% c(
      "log_av_ds",
      "log_diameter",
      "exposure",
      "mosscover",
      "snag"
    )
  ) %>%
  select(
    dataset,
    model,
    permutation_scope,
    by,
    term,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  arrange(dataset, desc(R2))

cat("\n=========================================================\n")
cat("SUBSTRATE TERMS\n")
cat("=========================================================\n")
print(substrate_terms, n = Inf)

write_csv2(
  substrate_terms,
  file.path(outdir, "substrate_terms.csv")
)

forest_time_effects <- permanova_all %>%
  filter(
    !is.na(forest_posthoc),
    term == "timespan"
  ) %>%
  select(
    dataset,
    forest = forest_posthoc,
    model,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`
  ) %>%
  arrange(dataset, desc(R2))

cat("\n=========================================================\n")
cat("FOREST-LEVEL TEMPORAL EFFECTS\n")
cat("=========================================================\n")
print(forest_time_effects, n = Inf)

write_csv2(
  forest_time_effects,
  file.path(outdir, "forest_level_timespan_effects.csv")
)

# =========================================================
# 5. BETA-DISPERSION DIAGNOSTICS
# =========================================================

run_betadisper_set <- function(X, meta, dataset_name, nperm = 9999) {
  dat <- clean_pa_inputs(X, meta)
  X <- dat$X
  meta <- dat$meta
  
  D <- make_jaccard(X)
  
  meta <- meta %>%
    mutate(
      forest_timespan = interaction(forest, timespan, sep = "_", drop = TRUE),
      country_timespan = interaction(country, timespan, sep = "_", drop = TRUE)
    )
  
  group_vars <- c(
    "timespan",
    "forest_timespan",
    "country_timespan"
  )
  
  out_list <- list()
  
  cat("\n=========================================================\n")
  cat("BETADISPER:", dataset_name, "\n")
  cat("=========================================================\n")
  
  for (gv in group_vars) {
    group <- factor(meta[[gv]])
    
    bd <- betadisper(D, group)
    bd_anova <- anova(bd)
    bd_perm <- permutest(bd, permutations = nperm)
    
    cat("\nGrouping:", gv, "\n")
    cat("\nANOVA\n")
    print(bd_anova)
    cat("\nPermutation test\n")
    print(bd_perm)
    
    out_list[[gv]] <- tibble(
      dataset = dataset_name,
      forest_posthoc = NA_character_,
      grouping = gv,
      Df_group = bd_anova$Df[1],
      Df_residual = bd_anova$Df[2],
      SumSq_group = bd_anova$`Sum Sq`[1],
      SumSq_residual = bd_anova$`Sum Sq`[2],
      MeanSq_group = bd_anova$`Mean Sq`[1],
      MeanSq_residual = bd_anova$`Mean Sq`[2],
      disp_F = bd_anova$`F value`[1],
      disp_anova_p = bd_anova$`Pr(>F)`[1],
      disp_perm_p = bd_perm$tab$`Pr(>F)`[1]
    )
  }
  
  bind_rows(out_list)
}

run_betadisper_by_forest <- function(X, meta, dataset_name, nperm = 9999) {
  dat <- clean_pa_inputs(X, meta)
  X <- dat$X
  meta <- dat$meta
  
  out_list <- list()
  
  cat("\n=========================================================\n")
  cat("FOREST-LEVEL BETADISPER:", dataset_name, "\n")
  cat("=========================================================\n")
  
  for (ff in levels(droplevels(meta$forest))) {
    idx <- meta$forest == ff
    
    X_f <- X[idx, , drop = FALSE]
    meta_f <- droplevels(meta[idx, , drop = FALSE])
    
    keep_taxa_f <- colSums(X_f, na.rm = TRUE) > 0
    X_f <- X_f[, keep_taxa_f, drop = FALSE]
    
    if (nrow(X_f) < 6 || length(unique(meta_f$timespan)) < 2) {
      next
    }
    
    D_f <- make_jaccard(X_f)
    group_f <- factor(meta_f$timespan)
    
    bd_f <- betadisper(D_f, group_f)
    bd_anova_f <- anova(bd_f)
    bd_perm_f <- permutest(bd_f, permutations = nperm)
    
    cat("\nForest:", ff, "\n")
    cat("\nANOVA\n")
    print(bd_anova_f)
    cat("\nPermutation test\n")
    print(bd_perm_f)
    
    out_list[[ff]] <- tibble(
      dataset = dataset_name,
      forest_posthoc = ff,
      grouping = "timespan",
      Df_group = bd_anova_f$Df[1],
      Df_residual = bd_anova_f$Df[2],
      SumSq_group = bd_anova_f$`Sum Sq`[1],
      SumSq_residual = bd_anova_f$`Sum Sq`[2],
      MeanSq_group = bd_anova_f$`Mean Sq`[1],
      MeanSq_residual = bd_anova_f$`Mean Sq`[2],
      disp_F = bd_anova_f$`F value`[1],
      disp_anova_p = bd_anova_f$`Pr(>F)`[1],
      disp_perm_p = bd_perm_f$tab$`Pr(>F)`[1]
    )
  }
  
  bind_rows(out_list)
}

betadisper_full <- run_betadisper_set(
  X = X_use_full,
  meta = meta_use_full,
  dataset_name = "full",
  nperm = nperm
)

betadisper_subset <- run_betadisper_set(
  X = X_use_subset,
  meta = meta_use_subset,
  dataset_name = "subset",
  nperm = nperm
)

betadisper_forest_full <- run_betadisper_by_forest(
  X = X_use_full,
  meta = meta_use_full,
  dataset_name = "full",
  nperm = nperm
)

betadisper_forest_subset <- run_betadisper_by_forest(
  X = X_use_subset,
  meta = meta_use_subset,
  dataset_name = "subset",
  nperm = nperm
)

betadisper_all <- bind_rows(
  betadisper_full,
  betadisper_subset
)

betadisper_forest_all <- bind_rows(
  betadisper_forest_full,
  betadisper_forest_subset
)

betadisper_all_with_forest <- bind_rows(
  betadisper_all,
  betadisper_forest_all
)

cat("\n=========================================================\n")
cat("BETA-DISPERSION SUMMARY: GLOBAL GROUPINGS\n")
cat("=========================================================\n")
print(betadisper_all, n = Inf)

cat("\n=========================================================\n")
cat("BETA-DISPERSION SUMMARY: FOREST-LEVEL TIMESPAN TESTS\n")
cat("=========================================================\n")
print(betadisper_forest_all, n = Inf)

write_csv2(
  betadisper_all,
  file.path(outdir, "betadisper_summary_global_groupings.csv")
)

write_csv2(
  betadisper_forest_all,
  file.path(outdir, "betadisper_summary_forest_level_timespan.csv")
)

write_csv2(
  betadisper_all_with_forest,
  file.path(outdir, "betadisper_summary_all.csv")
)

# =========================================================
# 5B. JOIN DISPERSION RESULTS TO PERMANOVA TERMS
# =========================================================

make_dispersion_lookup <- function(betadisper_all_with_forest) {
  betadisper_all_with_forest %>%
    mutate(
      forest_join = if_else(
        is.na(forest_posthoc),
        "Global",
        forest_posthoc
      ),
      dispersion_join_group = case_when(
        grouping == "timespan" ~ "timespan",
        grouping == "forest_timespan" ~ "forest:timespan",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(dispersion_join_group)) %>%
    select(
      dataset,
      forest_join,
      dispersion_join_group,
      disp_Df_group = Df_group,
      disp_Df_residual = Df_residual,
      disp_F,
      disp_anova_p,
      disp_perm_p
    )
}

add_dispersion_to_permanova <- function(permanova_df, dispersion_lookup) {
  permanova_df %>%
    mutate(
      forest_join = if_else(
        is.na(forest_posthoc),
        "Global",
        forest_posthoc
      ),
      dispersion_join_group = case_when(
        is.na(forest_posthoc) & term == "timespan" ~ "timespan",
        is.na(forest_posthoc) & term == "forest:timespan" ~ "forest:timespan",
        !is.na(forest_posthoc) & term == "timespan" ~ "timespan",
        TRUE ~ NA_character_
      )
    ) %>%
    left_join(
      dispersion_lookup,
      by = c(
        "dataset",
        "forest_join",
        "dispersion_join_group"
      )
    )
}

dispersion_lookup <- make_dispersion_lookup(
  betadisper_all_with_forest = betadisper_all_with_forest
)

permanova_all_with_dispersion <- add_dispersion_to_permanova(
  permanova_df = permanova_all,
  dispersion_lookup = dispersion_lookup
)

write_csv2(
  permanova_all_with_dispersion,
  file.path(outdir, "permanova_model_set_results_with_dispersion_raw.csv")
)

cat("\n=========================================================\n")
cat("PERMANOVA RESULTS WITH MATCHED DISPERSION COLUMNS\n")
cat("=========================================================\n")
print(permanova_all_with_dispersion, n = Inf)

# =========================================================
# 6. EXTRACT DISTANCES TO CENTROID
# =========================================================

extract_betadisper_distances <- function(X, meta, group_var, dataset_name) {
  dat <- clean_pa_inputs(X, meta)
  X <- dat$X
  meta <- dat$meta
  
  meta <- meta %>%
    mutate(
      forest_timespan = interaction(forest, timespan, sep = "_", drop = TRUE),
      country_timespan = interaction(country, timespan, sep = "_", drop = TRUE)
    )
  
  D <- make_jaccard(X)
  
  group <- factor(meta[[group_var]])
  bd <- betadisper(D, group)
  
  tibble(
    dataset = dataset_name,
    group_var = group_var,
    group = group,
    distance_to_centroid = bd$distances
  ) %>%
    group_by(dataset, group_var, group) %>%
    summarise(
      n = n(),
      mean_distance = mean(distance_to_centroid, na.rm = TRUE),
      sd_distance = sd(distance_to_centroid, na.rm = TRUE),
      se_distance = sd_distance / sqrt(n),
      median_distance = median(distance_to_centroid, na.rm = TRUE),
      min_distance = min(distance_to_centroid, na.rm = TRUE),
      max_distance = max(distance_to_centroid, na.rm = TRUE),
      .groups = "drop"
    )
}

dispersion_distances_all <- bind_rows(
  extract_betadisper_distances(
    X = X_use_full,
    meta = meta_use_full,
    group_var = "timespan",
    dataset_name = "full"
  ),
  extract_betadisper_distances(
    X = X_use_subset,
    meta = meta_use_subset,
    group_var = "timespan",
    dataset_name = "subset"
  ),
  extract_betadisper_distances(
    X = X_use_full,
    meta = meta_use_full,
    group_var = "forest_timespan",
    dataset_name = "full"
  ),
  extract_betadisper_distances(
    X = X_use_subset,
    meta = meta_use_subset,
    group_var = "forest_timespan",
    dataset_name = "subset"
  ),
  extract_betadisper_distances(
    X = X_use_full,
    meta = meta_use_full,
    group_var = "country_timespan",
    dataset_name = "full"
  ),
  extract_betadisper_distances(
    X = X_use_subset,
    meta = meta_use_subset,
    group_var = "country_timespan",
    dataset_name = "subset"
  )
)

cat("\n=========================================================\n")
cat("DISTANCE-TO-CENTROID SUMMARIES\n")
cat("=========================================================\n")
print(dispersion_distances_all, n = Inf)

write_csv2(
  dispersion_distances_all,
  file.path(outdir, "betadisper_group_distances.csv")
)

# =========================================================
# 7. PAIRWISE BETA-DISPERSION TESTS
# =========================================================

run_pairwise_betadisper_tukey <- function(X, meta, dataset_name, group_var) {
  dat <- clean_pa_inputs(X, meta)
  X <- dat$X
  meta <- dat$meta
  
  meta <- meta %>%
    mutate(
      forest_timespan = interaction(forest, timespan, sep = "_", drop = TRUE),
      country_timespan = interaction(country, timespan, sep = "_", drop = TRUE)
    )
  
  D <- make_jaccard(X)
  group <- factor(meta[[group_var]])
  
  bd <- betadisper(D, group)
  tuk <- TukeyHSD(bd)
  
  as.data.frame(tuk$group) %>%
    rownames_to_column("contrast") %>%
    as_tibble() %>%
    mutate(
      dataset = dataset_name,
      group_var = group_var,
      .before = 1
    )
}

pairwise_dispersion_all <- bind_rows(
  run_pairwise_betadisper_tukey(
    X = X_use_full,
    meta = meta_use_full,
    dataset_name = "full",
    group_var = "forest_timespan"
  ),
  run_pairwise_betadisper_tukey(
    X = X_use_subset,
    meta = meta_use_subset,
    dataset_name = "subset",
    group_var = "forest_timespan"
  )
)

cat("\n=========================================================\n")
cat("PAIRWISE BETA-DISPERSION TESTS: FOREST x TIMESPAN\n")
cat("=========================================================\n")
print(pairwise_dispersion_all, n = Inf)

write_csv2(
  pairwise_dispersion_all,
  file.path(outdir, "betadisper_pairwise_forest_timespan_tukey.csv")
)

# =========================================================
# 8. MAIN PERMANOVA TABLE FOR MANUSCRIPT
# Core ecological terms only, with matched dispersion tests
# =========================================================

tab_global_unadjusted <- permanova_all_with_dispersion %>%
  filter(
    is.na(forest_posthoc),
    str_detect(model, "^D ~ forest \\+ timespan \\+ forest:timespan$"),
    permutation_scope == "free",
    by == "terms",
    term %in% c("timespan", "forest:timespan")
  ) %>%
  mutate(
    model_level = "Global",
    model_type = "Unadjusted",
    model_label = "forest + period + forest × period",
    formula_used = "D ~ forest + timespan + forest:timespan",
    term_label = clean_permanova_term(term)
  )

tab_global_adjusted <- permanova_all_with_dispersion %>%
  filter(
    is.na(forest_posthoc),
    str_detect(model, "log_av_ds"),
    permutation_scope == "free",
    by == "terms",
    term %in% c("timespan", "forest:timespan")
  ) %>%
  mutate(
    model_level = "Global",
    model_type = "Substrate-adjusted",
    model_label = "forest + substrate + period + forest × period",
    formula_used = paste(
      "D ~ forest + log_av_ds + log_diameter + exposure +",
      "mosscover + snag + timespan + forest:timespan"
    ),
    term_label = clean_permanova_term(term)
  )

tab_forest_level <- permanova_all_with_dispersion %>%
  filter(
    !is.na(forest_posthoc),
    term == "timespan"
  ) %>%
  mutate(
    model_level = "Forest-level",
    model_type = "Substrate-adjusted",
    model_label = paste0("forest: ", forest_posthoc),
    formula_used = paste(
      "D_forest ~ log_av_ds + log_diameter + exposure +",
      "mosscover + snag + timespan"
    ),
    term_label = "Survey period"
  )

main_permanova_table <- bind_rows(
  tab_global_unadjusted,
  tab_global_adjusted,
  tab_forest_level
) %>%
  mutate(
    Dataset = clean_dataset_label(dataset),
    `Model level` = model_level,
    Model = model_label,
    `Formula used` = formula_used,
    Term = term_label,
    `R²` = format_num(R2, digits = 4),
    F = format_num(F, digits = 2),
    P = format_p(`Pr(>F)`),
    `Dispersion grouping` = clean_dispersion_grouping(dispersion_join_group),
    `Dispersion F` = format_num(disp_F, digits = 2),
    `Dispersion P` = format_p(disp_perm_p)
  ) %>%
  select(
    Dataset,
    `Model level`,
    Model,
    `Formula used`,
    Term,
    Df,
    `R²`,
    F,
    P,
    `Dispersion grouping`,
    `Dispersion F`,
    `Dispersion P`
  ) %>%
  arrange(
    factor(Dataset, levels = dataset_order),
    factor(`Model level`, levels = model_level_order),
    factor(Model, levels = model_order),
    factor(Term, levels = term_order_main)
  )

cat("\n=========================================================\n")
cat("MAIN PERMANOVA TABLE WITH MATCHED DISPERSION TESTS\n")
cat("=========================================================\n")
print(main_permanova_table, n = Inf)

write_csv2(
  main_permanova_table,
  file.path(outdir, "main_permanova_table_core_terms_with_dispersion.csv")
)

# =========================================================
# 9. SUPPLEMENTARY PERMANOVA TABLE WITH MATCHED DISPERSION
# Full model outputs, including substrate terms, residuals, totals
# =========================================================

supplementary_permanova_table <- permanova_all_with_dispersion %>%
  mutate(
    Dataset = clean_dataset_label(dataset),
    Forest = if_else(
      is.na(forest_posthoc),
      "Global",
      forest_posthoc
    ),
    Model = case_when(
      str_detect(model, "^D ~ forest \\+ timespan") ~
        "forest + period + forest × period",
      str_detect(model, "^D ~ forest \\+ log_av_ds") ~
        "forest + substrate + period + forest × period",
      str_detect(model, "^Forest-level model") ~
        paste0("forest: ", forest_posthoc),
      TRUE ~ model
    ),
    `Formula used` = model,
    `Permutation scope` = permutation_scope,
    `Sequential test` = by,
    Term = clean_permanova_term(term),
    `Sum of squares` = format_num(SumOfSqs, digits = 4),
    `R²` = format_num(R2, digits = 4),
    F = format_num(F, digits = 2),
    P = format_p(`Pr(>F)`),
    `Dispersion grouping` = clean_dispersion_grouping(dispersion_join_group),
    `Dispersion Df group` = disp_Df_group,
    `Dispersion Df residual` = disp_Df_residual,
    `Dispersion F` = format_num(disp_F, digits = 2),
    `Dispersion P` = format_p(disp_perm_p)
  ) %>%
  select(
    Dataset,
    Forest,
    Model,
    `Formula used`,
    `Permutation scope`,
    `Sequential test`,
    Term,
    Df,
    `Sum of squares`,
    `R²`,
    F,
    P,
    `Dispersion grouping`,
    `Dispersion Df group`,
    `Dispersion Df residual`,
    `Dispersion F`,
    `Dispersion P`
  ) %>%
  arrange(
    factor(Dataset, levels = dataset_order),
    factor(Forest, levels = c("Global", forest_levels)),
    factor(Model, levels = model_order),
    match(Term, term_order_supp)
  )

cat("\n=========================================================\n")
cat("SUPPLEMENTARY PERMANOVA TABLE WITH MATCHED DISPERSION\n")
cat("=========================================================\n")
print(supplementary_permanova_table, n = Inf)

write_csv2(
  supplementary_permanova_table,
  file.path(outdir, "supplementary_permanova_table_with_dispersion.csv")
)

# =========================================================
# 10. SUPPLEMENTARY STANDALONE BETA-DISPERSION TABLE
# =========================================================

supplementary_dispersion_table <- betadisper_all_with_forest %>%
  mutate(
    Dataset = clean_dataset_label(dataset),
    Forest = if_else(
      is.na(forest_posthoc),
      "Global",
      forest_posthoc
    ),
    Grouping = clean_dispersion_grouping(grouping),
    `Df group` = Df_group,
    `Df residual` = Df_residual,
    `Sum Sq group` = format_num(SumSq_group, digits = 4),
    `Sum Sq residual` = format_num(SumSq_residual, digits = 4),
    `Mean Sq group` = format_num(MeanSq_group, digits = 4),
    `Mean Sq residual` = format_num(MeanSq_residual, digits = 4),
    F = format_num(disp_F, digits = 2),
    `ANOVA P` = format_p(disp_anova_p),
    `Permutation P` = format_p(disp_perm_p)
  ) %>%
  select(
    Dataset,
    Forest,
    Grouping,
    `Df group`,
    `Df residual`,
    `Sum Sq group`,
    `Sum Sq residual`,
    `Mean Sq group`,
    `Mean Sq residual`,
    F,
    `ANOVA P`,
    `Permutation P`
  ) %>%
  arrange(
    factor(Dataset, levels = dataset_order),
    factor(Forest, levels = c("Global", forest_levels)),
    Grouping
  )

cat("\n=========================================================\n")
cat("SUPPLEMENTARY STANDALONE BETA-DISPERSION TABLE\n")
cat("=========================================================\n")
print(supplementary_dispersion_table, n = Inf)

write_csv2(
  supplementary_dispersion_table,
  file.path(outdir, "supplementary_betadisper_table.csv")
)

# =========================================================
# 11. PERMANOVA SUMMARY FIGURE
# Combined global + forest-level temporal effects
# =========================================================

dataset_labels <- c(
  full   = "Full dataset",
  subset = "Substrate-matched subset"
)

dataset_cols <- c(
  "Full dataset" = "grey20",
  "Substrate-matched subset" = "grey65"
)

dataset_shapes <- c(
  "Full dataset" = 21,
  "Substrate-matched subset" = 24
)

global_effects <- main_temporal_terms %>%
  mutate(
    dataset_label = recode(dataset, !!!dataset_labels),
    model_type = case_when(
      str_detect(model, "log_av_ds") ~ "Substrate-adjusted model",
      TRUE ~ "Global unadjusted model"
    ),
    term_label = recode(
      term,
      "timespan" = "Survey period",
      "forest:timespan" = "Forest × period"
    ),
    row_label = paste(model_type, term_label, sep = ": "),
    section = "Global PERMANOVA terms",
    p_label = case_when(
      `Pr(>F)` <= 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", number(`Pr(>F)`, accuracy = 0.001))
    )
  ) %>%
  select(
    dataset,
    dataset_label,
    section,
    row_label,
    term,
    model_type,
    R2,
    F,
    `Pr(>F)`,
    p_label
  )

global_row_levels <- c(
  "Global unadjusted model: Survey period",
  "Global unadjusted model: Forest × period",
  "Substrate-adjusted model: Survey period",
  "Substrate-adjusted model: Forest × period"
)

forest_order <- forest_time_effects %>%
  group_by(forest) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_R2)) %>%
  pull(forest)

forest_effects <- forest_time_effects %>%
  mutate(
    dataset_label = recode(dataset, !!!dataset_labels),
    section = "Forest-level temporal effects",
    row_label = paste0(forest, ": Survey period"),
    term = "timespan",
    model_type = "Forest-level model",
    p_label = case_when(
      `Pr(>F)` <= 0.001 ~ "p < 0.001",
      TRUE ~ paste0("p = ", number(`Pr(>F)`, accuracy = 0.001))
    )
  ) %>%
  select(
    dataset,
    dataset_label,
    section,
    row_label,
    term,
    model_type,
    R2,
    F,
    `Pr(>F)`,
    p_label
  )

forest_row_levels <- paste0(forest_order, ": Survey period")

effect_plot_df <- bind_rows(
  global_effects,
  forest_effects
) %>%
  mutate(
    dataset = factor(dataset, levels = c("full", "subset")),
    dataset_label = factor(
      dataset_label,
      levels = unname(dataset_labels[c("full", "subset")])
    ),
    section = factor(
      section,
      levels = c(
        "Global PERMANOVA terms",
        "Forest-level temporal effects"
      )
    ),
    row_label = factor(
      row_label,
      levels = rev(c(global_row_levels, forest_row_levels))
    )
  )

p_permanova_summary <- ggplot(
  effect_plot_df,
  aes(
    x = R2,
    y = row_label,
    group = row_label
  )
) +
  geom_line(
    aes(group = row_label),
    linewidth = 0.35,
    colour = "grey75"
  ) +
  geom_point(
    aes(
      shape = dataset_label,
      fill = dataset_label
    ),
    size = 3.2,
    colour = "grey20",
    stroke = 0.45
  ) +
  facet_grid(
    section ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_shape_manual(
    values = dataset_shapes,
    name = NULL
  ) +
  scale_fill_manual(
    values = dataset_cols,
    name = NULL
  ) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = expression(PERMANOVA~R^2),
    y = NULL,
    caption = paste(
      "All plotted PERMANOVA terms were significant at p < 0.001.",
      "Beta-dispersion diagnostics are reported in the associated tables."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "grey35"),
    strip.text.y = element_text(face = "bold", angle = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      colour = "grey25"
    )
  )

print(p_permanova_summary)

ggsave(
  filename = file.path(outdir, "permanova_summary_effect_sizes_combined.png"),
  plot = p_permanova_summary,
  width = 8.5,
  height = 6.8,
  dpi = 400,
  bg = "white"
)

# =========================================================
# 12. FINAL RUN SUMMARY
# =========================================================

cat("\n=========================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=========================================================\n")
cat("Output directory:", outdir, "\n\n")

cat("Main outputs written:\n")
cat("- permanova_model_set_results.csv\n")
cat("- permanova_model_set_results_with_dispersion_raw.csv\n")
cat("- main_permanova_table_core_terms_with_dispersion.csv\n")
cat("- supplementary_permanova_table_with_dispersion.csv\n")
cat("- supplementary_betadisper_table.csv\n")
cat("- betadisper_summary_all.csv\n")
cat("- betadisper_group_distances.csv\n")
cat("- betadisper_pairwise_forest_timespan_tukey.csv\n")
cat("- permanova_summary_effect_sizes_combined.png\n")