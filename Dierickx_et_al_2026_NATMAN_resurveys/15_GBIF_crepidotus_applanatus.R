# ------------------------------------------------------------------
# GBIF exploratory reporting trend:
# Crepidotus applanatus s.l. relative to all European Crepidotus
#
# Important:
# This tests changes in the relative share of GBIF records, not
# biological abundance or occupancy.
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------------

required_packages <- c(
  "rgbif",
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "tibble",
  "glmmTMB"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(
  lapply(required_packages, library, character.only = TRUE)
)

# ------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------

output_dir <- "gbif_crepidotus_trend"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

start_year <- 1990L
end_year   <- 2024L

# GBIF Backbone taxon keys checked 2026-07-02
crepidotus_key <- 2526680L
applanatus_key <- 5240984L
malachius_key  <- 6014641L

# Strict concept:
#   Crepidotus applanatus only
#
# Broad sensitivity concept:
#   Crepidotus applanatus + Crepidotus malachius
#
# The broad version is exploratory and checks whether any apparent
# pattern depends on the delimitation of the C. applanatus complex.

target_keys_strict <- c(applanatus_key)
target_keys_broad  <- c(applanatus_key, malachius_key)

# European countries included in the analysis.
# Russia and Turkey are excluded because their country-wide GBIF records
# cannot be separated consistently into European and Asian portions.

europe_iso2 <- c(
  "AL", "AD", "AT", "BA", "BE", "BG", "BY", "CH", "CY", "CZ",
  "DE", "DK", "EE", "ES", "FI", "FR", "GB", "GR", "HR", "HU",
  "IE", "IS", "IT", "LI", "LT", "LU", "LV", "MC", "MD", "ME",
  "MK", "MT", "NL", "NO", "PL", "PT", "RO", "RS", "SE", "SI",
  "SK", "SM", "UA", "VA", "XK"
)

retained_basis <- c(
  "HUMAN_OBSERVATION",
  "PRESERVED_SPECIMEN"
)

readr::write_csv(
  tibble::tibble(
    parameter = c(
      "start_year",
      "end_year",
      "Crepidotus_genus_key",
      "Crepidotus_applanatus_key",
      "Crepidotus_malachius_key",
      "record_basis"
    ),
    value = c(
      start_year,
      end_year,
      crepidotus_key,
      applanatus_key,
      malachius_key,
      paste(retained_basis, collapse = "; ")
    )
  ),
  file.path(output_dir, "analysis_parameters.csv")
)

# ------------------------------------------------------------------
# 2. Check GBIF credentials
# ------------------------------------------------------------------
#
# Put these in your .Renviron file, then restart R:
#
# GBIF_USER="your_GBIF_username"
# GBIF_PWD="your_GBIF_password"
# GBIF_EMAIL="your_email"
#

gbif_credentials <- c(
  Sys.getenv("GBIF_USER"),
  Sys.getenv("GBIF_PWD"),
  Sys.getenv("GBIF_EMAIL")
)

if (!all(nzchar(gbif_credentials))) {
  stop(
    "GBIF credentials are missing. Add GBIF_USER, GBIF_PWD and ",
    "GBIF_EMAIL to .Renviron, restart R, then run the script again."
  )
}

# ------------------------------------------------------------------
# 3. Request and retrieve the GBIF download
# ------------------------------------------------------------------

download_key_file <- file.path(output_dir, "gbif_download_key.txt")

if (file.exists(download_key_file)) {
  
  gbif_download_key <- readLines(download_key_file, warn = FALSE)[1]
  message("Reusing existing GBIF download key: ", gbif_download_key)
  
} else {
  
  gbif_download_key <- rgbif::occ_download(
    rgbif::pred("taxonKey", crepidotus_key),
    rgbif::pred_in("country", europe_iso2),
    rgbif::pred_gte("year", start_year),
    rgbif::pred_lte("year", end_year),
    rgbif::pred_in("basisOfRecord", retained_basis),
    rgbif::pred("occurrenceStatus", "PRESENT"),
    format = "DWCA"
  )
  
  gbif_download_key <- as.character(gbif_download_key)[1]
  
  writeLines(
    gbif_download_key,
    con = download_key_file
  )
  
  message("New GBIF download requested: ", gbif_download_key)
}

rgbif::occ_download_wait(gbif_download_key)

download_meta <- rgbif::occ_download_meta(gbif_download_key)

download_doi <- if (!is.null(download_meta$doi)) {
  as.character(download_meta$doi)
} else {
  NA_character_
}

writeLines(
  c(
    paste("GBIF download key:", gbif_download_key),
    paste("GBIF download DOI:", download_doi)
  ),
  con = file.path(output_dir, "gbif_download_citation.txt")
)

message("GBIF DOI: ", download_doi)

raw_rds_file <- file.path(output_dir, "gbif_crepidotus_raw.rds")

if (file.exists(raw_rds_file)) {
  
  gbif_raw <- readRDS(raw_rds_file)
  message("Loaded cached imported GBIF data.")
  
} else {
  
  gbif_download_file <- rgbif::occ_download_get(
    key = gbif_download_key,
    path = output_dir,
    overwrite = FALSE
  )
  
  gbif_raw <- rgbif::occ_download_import(
    gbif_download_file,
    fill = TRUE
  ) |>
    tibble::as_tibble()
  
  saveRDS(
    gbif_raw,
    raw_rds_file
  )
}

# ------------------------------------------------------------------
# 4. Check and standardise required fields
# ------------------------------------------------------------------

required_columns <- c(
  "year",
  "countryCode",
  "basisOfRecord",
  "taxonRank",
  "taxonKey",
  "datasetKey",
  "scientificName"
)

missing_columns <- setdiff(required_columns, names(gbif_raw))

if (length(missing_columns) > 0) {
  stop(
    "The GBIF download lacks required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

optional_columns <- c(
  "acceptedTaxonKey",
  "speciesKey",
  "acceptedScientificName",
  "gbifID"
)

for (column_name in optional_columns) {
  if (!column_name %in% names(gbif_raw)) {
    gbif_raw[[column_name]] <- NA
  }
}

as_integer_safely <- function(x) {
  suppressWarnings(as.integer(as.character(x)))
}

replace_missing_label <- function(x, label) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- label
  x
}

gbif_clean <- gbif_raw |>
  dplyr::mutate(
    year = as_integer_safely(year),
    countryCode = toupper(trimws(as.character(countryCode))),
    basisOfRecord = toupper(trimws(as.character(basisOfRecord))),
    taxonRank = toupper(trimws(as.character(taxonRank))),
    dataset_id = replace_missing_label(datasetKey, "missing_dataset_key"),
    
    accepted_key = dplyr::coalesce(
      as_integer_safely(acceptedTaxonKey),
      as_integer_safely(taxonKey),
      as_integer_safely(speciesKey)
    ),
    
    is_applanatus = accepted_key %in% target_keys_strict,
    is_applanatus_broad = accepted_key %in% target_keys_broad
  ) |>
  dplyr::filter(
    dplyr::between(year, start_year, end_year),
    countryCode %in% europe_iso2,
    basisOfRecord %in% retained_basis,
    taxonRank == "SPECIES",
    !is.na(accepted_key)
  )

if ("gbifID" %in% names(gbif_clean)) {
  gbif_clean <- gbif_clean |>
    dplyr::distinct(gbifID, .keep_all = TRUE)
}

if (nrow(gbif_clean) == 0) {
  stop("No usable records remained after filtering.")
}

if (sum(gbif_clean$is_applanatus) == 0) {
  stop(
    "No records matched the strict Crepidotus applanatus concept. ",
    "Inspect the taxonomic audit output."
  )
}

# ------------------------------------------------------------------
# 5. Audit the taxonomic standardisation
# ------------------------------------------------------------------

target_taxon_audit <- gbif_clean |>
  dplyr::filter(is_applanatus_broad) |>
  dplyr::count(
    accepted_key,
    acceptedScientificName,
    scientificName,
    sort = TRUE
  )

readr::write_csv(
  target_taxon_audit,
  file.path(output_dir, "target_taxon_audit.csv")
)

all_taxa_audit <- gbif_clean |>
  dplyr::count(
    accepted_key,
    acceptedScientificName,
    sort = TRUE
  )

readr::write_csv(
  all_taxa_audit,
  file.path(output_dir, "all_crepidotus_taxa_audit.csv")
)

print(target_taxon_audit)

# ------------------------------------------------------------------
# 6. Aggregate within reporting strata
# ------------------------------------------------------------------
#
# The denominator is all species-level Crepidotus records within the
# same year × country × dataset × basis-of-record stratum.
# This reduces the influence of broad temporal growth in GBIF reporting.
#

stratum_counts <- gbif_clean |>
  dplyr::group_by(
    year,
    countryCode,
    dataset_id,
    basisOfRecord
  ) |>
  dplyr::summarise(
    n_all_crepidotus = dplyr::n(),
    n_applanatus = sum(is_applanatus),
    n_applanatus_broad = sum(is_applanatus_broad),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    n_other_strict = n_all_crepidotus - n_applanatus,
    n_other_broad = n_all_crepidotus - n_applanatus_broad
  )

readr::write_csv(
  stratum_counts,
  file.path(output_dir, "crepidotus_reporting_strata.csv")
)

# ------------------------------------------------------------------
# 7. Annual descriptive summaries
# ------------------------------------------------------------------

annual_summary <- stratum_counts |>
  dplyr::select(
    year,
    basisOfRecord,
    n_all_crepidotus,
    n_applanatus,
    n_applanatus_broad
  ) |>
  tidyr::pivot_longer(
    cols = c(n_applanatus, n_applanatus_broad),
    names_to = "taxonomic_concept",
    values_to = "n_target"
  ) |>
  dplyr::group_by(
    year,
    basisOfRecord,
    taxonomic_concept
  ) |>
  dplyr::summarise(
    n_target = sum(n_target),
    n_all_crepidotus = sum(n_all_crepidotus),
    relative_reporting_share = n_target / n_all_crepidotus,
    .groups = "drop"
  ) |>
  dplyr::mutate(
    taxonomic_concept = dplyr::recode(
      taxonomic_concept,
      n_applanatus = "Strict: C. applanatus",
      n_applanatus_broad = "Broad: C. applanatus + C. malachius"
    )
  )

readr::write_csv(
  annual_summary,
  file.path(output_dir, "annual_relative_reporting_summary.csv")
)

annual_plot <- ggplot2::ggplot(
  annual_summary,
  ggplot2::aes(
    x = year,
    y = relative_reporting_share,
    linetype = basisOfRecord,
    group = basisOfRecord
  )
) +
  ggplot2::geom_line(
    colour = "black",
    linewidth = 0.65
  ) +
  ggplot2::facet_wrap(
    ~ taxonomic_concept,
    ncol = 1,
    scales = "free_y"
  ) +
  ggplot2::labs(
    x = "Year",
    y = "Relative share of species-level Crepidotus records",
    linetype = "Record basis"
  ) +
  ggplot2::theme_classic()

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "annual_relative_reporting_trend.png"
  ),
  plot = annual_plot,
  width = 180,
  height = 160,
  units = "mm",
  dpi = 600
)

print(annual_plot)

# ------------------------------------------------------------------
# 8. Beta-binomial reporting-trend models
# ------------------------------------------------------------------
#
# Response:
# target records / all other Crepidotus records
#
# The year effect estimates whether the odds that a Crepidotus record
# belongs to the target taxon changed through time.
#

make_model_data <- function(data, target_column, concept_label) {
  
  model_data <- data |>
    dplyr::transmute(
      year,
      year_c = year - mean(year),
      
      countryCode = factor(countryCode),
      dataset_id = factor(dataset_id),
      basisOfRecord = factor(basisOfRecord),
      
      target_n = .data[[target_column]],
      other_n = n_all_crepidotus - .data[[target_column]],
      
      taxonomic_concept = concept_label
    ) |>
    dplyr::filter(
      target_n >= 0,
      other_n >= 0,
      target_n + other_n > 0
    )
  
  model_data
}

strict_model_data <- make_model_data(
  stratum_counts,
  target_column = "n_applanatus",
  concept_label = "Strict: C. applanatus"
)

broad_model_data <- make_model_data(
  stratum_counts,
  target_column = "n_applanatus_broad",
  concept_label = "Broad: C. applanatus + C. malachius"
)

fit_reporting_model <- function(data, include_basis = TRUE) {
  
  if (nrow(data) < 5 ||
      sum(data$target_n) == 0 ||
      sum(data$other_n) == 0) {
    return(NULL)
  }
  
  model_formula <- if (include_basis) {
    stats::as.formula(
      "cbind(target_n, other_n) ~ year_c + basisOfRecord +
       (1 | countryCode) + (1 | dataset_id)"
    )
  } else {
    stats::as.formula(
      "cbind(target_n, other_n) ~ year_c +
       (1 | countryCode) + (1 | dataset_id)"
    )
  }
  
  tryCatch(
    glmmTMB::glmmTMB(
      formula = model_formula,
      family = glmmTMB::betabinomial(link = "logit"),
      data = data
    ),
    error = function(e) {
      message("Model failed: ", conditionMessage(e))
      NULL
    }
  )
}

models <- list(
  strict_all = fit_reporting_model(strict_model_data),
  broad_all = fit_reporting_model(broad_model_data)
)

for (record_basis in levels(strict_model_data$basisOfRecord)) {
  
  strict_subset <- strict_model_data |>
    dplyr::filter(basisOfRecord == record_basis)
  
  broad_subset <- broad_model_data |>
    dplyr::filter(basisOfRecord == record_basis)
  
  models[[paste0("strict_", record_basis)]] <- fit_reporting_model(
    strict_subset,
    include_basis = FALSE
  )
  
  models[[paste0("broad_", record_basis)]] <- fit_reporting_model(
    broad_subset,
    include_basis = FALSE
  )
}

extract_year_effect <- function(model, concept, subset_label) {
  
  if (is.null(model)) {
    return(
      tibble::tibble(
        taxonomic_concept = concept,
        subset = subset_label,
        log_odds_change_per_year = NA_real_,
        odds_ratio_per_decade = NA_real_,
        ci_lower_per_decade = NA_real_,
        ci_upper_per_decade = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  coefficient_table <- summary(model)$coefficients$cond
  
  if (!"year_c" %in% rownames(coefficient_table)) {
    stop("Could not extract the year effect from the fitted model.")
  }
  
  estimate <- coefficient_table["year_c", "Estimate"]
  se <- coefficient_table["year_c", "Std. Error"]
  p_value <- coefficient_table["year_c", "Pr(>|z|)"]
  
  tibble::tibble(
    taxonomic_concept = concept,
    subset = subset_label,
    log_odds_change_per_year = estimate,
    odds_ratio_per_decade = exp(estimate * 10),
    ci_lower_per_decade = exp((estimate - 1.96 * se) * 10),
    ci_upper_per_decade = exp((estimate + 1.96 * se) * 10),
    p_value = p_value
  )
}

model_results <- dplyr::bind_rows(
  extract_year_effect(
    models$strict_all,
    "Strict: C. applanatus",
    "All retained record bases"
  ),
  extract_year_effect(
    models$broad_all,
    "Broad: C. applanatus + C. malachius",
    "All retained record bases"
  ),
  extract_year_effect(
    models$strict_HUMAN_OBSERVATION,
    "Strict: C. applanatus",
    "Human observations only"
  ),
  extract_year_effect(
    models$strict_PRESERVED_SPECIMEN,
    "Strict: C. applanatus",
    "Preserved specimens only"
  ),
  extract_year_effect(
    models$broad_HUMAN_OBSERVATION,
    "Broad: C. applanatus + C. malachius",
    "Human observations only"
  ),
  extract_year_effect(
    models$broad_PRESERVED_SPECIMEN,
    "Broad: C. applanatus + C. malachius",
    "Preserved specimens only"
  )
)

readr::write_csv(
  model_results,
  file.path(output_dir, "beta_binomial_reporting_trend_models.csv")
)

print(model_results)

# ------------------------------------------------------------------
# 9. Dataset and country influence diagnostics
# ------------------------------------------------------------------

top_target_datasets <- gbif_clean |>
  dplyr::filter(is_applanatus_broad) |>
  dplyr::count(
    dataset_id,
    basisOfRecord,
    sort = TRUE
  ) |>
  dplyr::mutate(
    proportion_of_broad_target_records = n / sum(n)
  ) |>
  dplyr::slice_head(n = 30)

top_target_countries <- gbif_clean |>
  dplyr::filter(is_applanatus_broad) |>
  dplyr::count(
    countryCode,
    sort = TRUE
  ) |>
  dplyr::mutate(
    proportion_of_broad_target_records = n / sum(n)
  )

readr::write_csv(
  top_target_datasets,
  file.path(output_dir, "top_broad_target_datasets.csv")
)

readr::write_csv(
  top_target_countries,
  file.path(output_dir, "top_broad_target_countries.csv")
)

print(top_target_datasets)
print(top_target_countries)

# ------------------------------------------------------------------
# 10. Interpretation checklist
# ------------------------------------------------------------------
#
# A cautious claim of increased relative reporting is only defensible if:
#
# 1. Strict and broad concepts show the same direction.
# 2. The odds ratio per decade is > 1 in both concepts.
# 3. Human observations and preserved specimens broadly agree.
# 4. The pattern is not dominated by one country or one large dataset.
#
# Even then, describe the result as an increase in relative GBIF
# reporting, not as proof of a Europe-wide biological increase.
# ------------------------------------------------------------------


# ------------------------------------------------------------------
# Strict Crepidotus applanatus reporting trend:
# human observations only
# ------------------------------------------------------------------

annual_applanatus_human <- stratum_counts |>
  dplyr::filter(basisOfRecord == "HUMAN_OBSERVATION") |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    n_target = sum(n_applanatus),
    n_all_crepidotus = sum(n_all_crepidotus),
    relative_reporting_share = n_target / n_all_crepidotus,
    .groups = "drop"
  )

applanatus_human_plot <- ggplot2::ggplot(
  annual_applanatus_human,
  ggplot2::aes(
    x = year,
    y = relative_reporting_share
  )
) +
  ggplot2::geom_line(
    linewidth = 0.7,
    colour = "black"
  ) +
  ggplot2::geom_point(
    size = 1.4,
    colour = "black"
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 0.1)
  ) +
  ggplot2::labs(
    x = "Year",
    y = "Relative share Crepidotus records"
  ) +
  ggplot2::labs(
    title = "Relative GBIF reporting of Crepidotus applanatus among European Crepidotus records"
  )
print(applanatus_human_plot)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "strict_Crepidotus_applanatus_human_observation_trend.png"
  ),
  plot = applanatus_human_plot,
  width = 180,
  height = 100,
  units = "mm",
  dpi = 600
)


download_meta <- rgbif::occ_download_meta(gbif_download_key)

download_doi <- download_meta$doi
download_doi

cat(
  "GBIF.org (2026). GBIF Occurrence Download. https://doi.org/",
  download_doi,
  "\n",
  sep = ""
)
