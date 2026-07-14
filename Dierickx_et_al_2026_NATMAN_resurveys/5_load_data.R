suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})
setwd("/kyukon/data/gent/vo/001/gvo00142/vsc45818/FRUITBODY_PAPER")
set.seed(42)

# =========================================================
# Shared helper functions
# =========================================================

standardize_forest <- function(x) {
  x <- str_trim(str_squish(as.character(x)))
  ifelse(x %in% c("Strødam", "Strodam"), "Strødam", x)
}

clean_chr_na <- function(x) {
  na_if(str_trim(as.character(x)), "")
}

recode_log_av_ds <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  dplyr::if_else(x_num == 6, 5, x_num)
}

join_tree_metadata <- function(dat, tree_meta_clean) {
  dat %>%
    left_join(tree_meta_clean, by = "tree_id", suffix = c("", "_meta")) %>%
    mutate(
      across(
        ends_with("_meta"),
        ~ coalesce(get(sub("_meta$", "", cur_column())), .x),
        .names = "{sub('_meta$', '', .col)}"
      )
    ) %>%
    select(-ends_with("_meta")) %>%
    mutate(
      log_av_ds = recode_log_av_ds(log_av_ds)
    )
}

drop_samples <- c("RAJ043_2001", "SON121_2001")

# =========================================================
# Read and clean full tree metadata
# =========================================================

tree_meta_clean <- read_csv2(
  "basic_data_JHC/tree_meta_clean.csv",
  show_col_types = FALSE
) %>%
  mutate(
    tree_id   = clean_chr_na(tree_id),
    forest    = standardize_forest(forest),
    country   = clean_chr_na(country),
    log_av_ds = recode_log_av_ds(log_av_ds)
  ) %>%
  distinct(tree_id, .keep_all = TRUE)

cat("\nDecay-stage check in tree_meta_clean after recoding:\n")
print(table(tree_meta_clean$log_av_ds, useNA = "ifany"))

stopifnot(!any(tree_meta_clean$log_av_ds == 6, na.rm = TRUE))

# =========================================================
# FULL DATASET
# =========================================================

dat_raw <- read_csv2(
  "basic_data_JHC/analysis_main_natman_pa_tree_year_TARGET.csv",
  show_col_types = FALSE
)

meta_cols_minimal <- c(
  "source_dataset", "country", "forest", "tree_id", "year", "sample_id"
)

taxa_cols <- setdiff(names(dat_raw), meta_cols_minimal)

taxa_keep <- taxa_cols[
  !str_detect(taxa_cols, "\\bsp\\.?$")
]

dat_raw <- dat_raw %>%
  select(all_of(meta_cols_minimal), all_of(taxa_keep)) %>%
  mutate(
    tree_id = clean_chr_na(tree_id),
    forest  = standardize_forest(forest),
    country = clean_chr_na(country),
    year    = suppressWarnings(as.integer(year))
  ) %>%
  filter(year != 2023L)

# Join all available metadata from tree_meta_clean
# log_av_ds is already recoded in tree_meta_clean and checked again after joining
dat_raw <- join_tree_metadata(dat_raw, tree_meta_clean)

# Re-identify taxa after the metadata join
taxa_cols_full <- taxa_keep[taxa_keep %in% names(dat_raw)]

meta <- dat_raw %>%
  select(-all_of(taxa_cols_full)) %>%
  mutate(
    timespan = case_when(
      year == 2001L ~ "previous",
      year %in% c(2021L, 2022L) ~ "current",
      TRUE ~ NA_character_
    ),
    timespan    = factor(timespan, levels = c("previous", "current")),
    survey_year = factor(year, levels = c(2001, 2021, 2022)),
    across(
      any_of(c("log_av_ds", "log_diameter", "mosscover")),
      ~ suppressWarnings(as.numeric(.x))
    ),
    log_av_ds = recode_log_av_ds(log_av_ds),
    across(
      any_of(c("country", "forest", "tree_id", "exposure")),
      factor
    )
  )

cat("\nDecay-stage check in full metadata before sample filtering:\n")
print(table(meta$log_av_ds, useNA = "ifany"))

stopifnot(!any(meta$log_av_ds == 6, na.rm = TRUE))

X <- dat_raw %>%
  select(all_of(taxa_cols_full)) %>%
  mutate(across(everything(), ~ as.integer(.x > 0))) %>%
  as.data.frame()

rownames(X) <- meta$sample_id

# Remove empty samples and low-richness / known bad samples
rich <- rowSums(X)

keep <- !(rownames(X) %in% drop_samples) & rich > 1

X <- X[keep, , drop = FALSE]
meta <- meta[match(rownames(X), meta$sample_id), ]

# Remove taxa absent after sample filtering
X <- X[, colSums(X) > 0, drop = FALSE]

cat("\nDecay-stage check in full active metadata after sample filtering:\n")
print(table(meta$log_av_ds, useNA = "ifany"))

stopifnot(
  identical(rownames(X), meta$sample_id),
  nrow(meta) == nrow(X),
  !any(names(meta) %in% names(X)),
  !any(meta$log_av_ds == 6, na.rm = TRUE)
)

X_use_full <- X
meta_use_full <- meta

cat("\nFull matrix dimensions:\n")
cat("Samples:", nrow(X_use_full), "\n")
cat("Active taxa:", ncol(X_use_full), "\n\n")

str(X_use_full)
str(meta_use_full)
print(table(meta_use_full$year))
print(table(meta_use_full$timespan))

# =========================================================
# SUBSET DATASET
# =========================================================

set.seed(42)

dat_raw <- read_csv2(
  "subset_multivariate_inputs/analysis_input_deterministic_single_year_2022.csv",
  show_col_types = FALSE
)

dat_raw <- dat_raw %>%
  mutate(
    tree_id = clean_chr_na(tree_id),
    forest  = standardize_forest(forest),
    country = clean_chr_na(country),
    year    = suppressWarnings(as.integer(year))
  ) %>%
  filter(year != 2023L)

# Join all available metadata from tree_meta_clean
# log_av_ds is already recoded in tree_meta_clean and checked again after joining
dat_raw <- join_tree_metadata(dat_raw, tree_meta_clean)

# ------------------------------------------------------------------
# Build community matrix
# ------------------------------------------------------------------

species_cols <- names(dat_raw)[
  vapply(dat_raw, is.numeric, logical(1)) &
    str_detect(names(dat_raw), "^[A-Z][a-z]+\\s[a-z]+")
]

species_cols <- species_cols[
  !str_detect(species_cols, "\\bsp\\.?$")
]

meta <- dat_raw %>%
  select(-all_of(species_cols)) %>%
  mutate(
    timespan = case_when(
      year == 2001L ~ "previous",
      year %in% c(2021L, 2022L) ~ "current",
      TRUE ~ NA_character_
    ),
    timespan    = factor(timespan, levels = c("previous", "current")),
    survey_year = factor(as.character(year)),
    across(
      any_of(c("log_av_ds", "log_diameter", "mosscover")),
      ~ suppressWarnings(as.numeric(.x))
    ),
    log_av_ds = recode_log_av_ds(log_av_ds),
    across(
      any_of(c("country", "forest", "tree_id", "exposure")),
      factor
    )
  )

cat("\nDecay-stage check in subset metadata before sample filtering:\n")
print(table(meta$log_av_ds, useNA = "ifany"))

stopifnot(!any(meta$log_av_ds == 6, na.rm = TRUE))

X <- dat_raw %>%
  select(all_of(species_cols)) %>%
  mutate(across(everything(), ~ as.integer(.x > 0))) %>%
  as.data.frame()

rownames(X) <- meta$sample_id

rich <- rowSums(X)

keep <- !(rownames(X) %in% drop_samples) & rich > 1

X <- X[keep, , drop = FALSE]
meta <- meta[match(rownames(X), meta$sample_id), ]

X <- X[, colSums(X) > 0, drop = FALSE]

cat("\nDecay-stage check in subset active metadata after sample filtering:\n")
print(table(meta$log_av_ds, useNA = "ifany"))

stopifnot(
  identical(rownames(X), meta$sample_id),
  nrow(meta) == nrow(X),
  !any(names(meta) %in% names(X)),
  !any(meta$log_av_ds == 6, na.rm = TRUE)
)

X_use_subset <- X
meta_use_subset <- meta

cat("\nSubset matrix dimensions:\n")
cat("Samples:", nrow(X_use_subset), "\n")
cat("Active taxa:", ncol(X_use_subset), "\n\n")

print(table(meta_use_subset$year))
print(table(meta_use_subset$timespan))

# =========================================================
# Final metadata checks
# =========================================================

cat("\nFinal decay-stage check:\n")

cat("\nFull dataset log_av_ds:\n")
print(table(meta_use_full$log_av_ds, useNA = "ifany"))

cat("\nSubset dataset log_av_ds:\n")
print(table(meta_use_subset$log_av_ds, useNA = "ifany"))

stopifnot(
  !any(meta_use_full$log_av_ds == 6, na.rm = TRUE),
  !any(meta_use_subset$log_av_ds == 6, na.rm = TRUE)
)

cat("\nCountry x log_av_ds, full dataset:\n")
print(table(meta_use_full$country, meta_use_full$log_av_ds, useNA = "ifany"))

cat("\nCountry x log_av_ds, subset dataset:\n")
print(table(meta_use_subset$country, meta_use_subset$log_av_ds, useNA = "ifany"))

# =========================================================
# Final active taxa comparison
# =========================================================

cat("\nFinal active taxa comparison:\n")
cat("Full active taxa:", ncol(X_use_full), "\n")
cat("Subset active taxa:", ncol(X_use_subset), "\n")

taxa_only_in_subset <- setdiff(names(X_use_subset), names(X_use_full))
taxa_only_in_full   <- setdiff(names(X_use_full), names(X_use_subset))

cat("\nTaxa only in subset:", length(taxa_only_in_subset), "\n")
cat("Taxa only in full:", length(taxa_only_in_full), "\n")

if (length(taxa_only_in_subset) > 0) {
  cat("\nTaxa only in subset:\n")
  print(taxa_only_in_subset)
}

if (length(taxa_only_in_full) > 0) {
  cat("\nFirst taxa only in full:\n")
  print(head(taxa_only_in_full, 50))
}

#### str ####
print(X_use_full[1:5,1:5])
str(meta_use_full)
print(X_use_subset[1:5,1:5])
str(meta_use_subset)