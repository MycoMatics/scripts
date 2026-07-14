suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
  library(tibble)
})

setwd("/data/gent/vo/001/gvo00142/vsc45818/FRUITBODY_PAPER")
set.seed(42)

outdir <- "basic_data_JHC"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------
# 0. Helper functions
#-----------------------------------------------------------
invalid_species_labels <- c(
  "", "\\", "/", "/ /", "//"
)

year_levels <- c("2001", "2021", "2022", "2023", "2021-2023")

capitalize_genus <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_squish(x)
  
  parts <- stringr::str_split(x, "\\s+", simplify = TRUE)
  parts <- as.matrix(parts)
  
  first_word <- parts[, 1]
  aff_idx <- !is.na(first_word) & tolower(first_word) == "aff."
  normal_idx <- !is.na(first_word) & !aff_idx
  
  parts[normal_idx, 1] <- stringr::str_to_title(parts[normal_idx, 1])
  
  if (ncol(parts) >= 2) {
    second_word <- parts[, 2]
    aff_second_idx <- aff_idx & !is.na(second_word) & second_word != ""
    parts[aff_second_idx, 2] <- stringr::str_to_title(parts[aff_second_idx, 2])
  }
  
  out <- apply(parts, 1, function(row) {
    row <- row[!is.na(row) & row != ""]
    paste(row, collapse = " ")
  })
  
  out[out == ""] <- NA_character_
  out
}

clean_species_vector <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_squish()
}

get_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) return(NA_character_)
  hit
}

derive_tree_metadata <- function(tree_id_vec) {
  tree_id_vec <- as.character(tree_id_vec)
  tree_id_up <- toupper(stringr::str_trim(tree_id_vec))
  
  tibble(
    tree_id = tree_id_vec,
    forest = case_when(
      stringr::str_detect(tree_id_up, "^(ZF|SON)")   ~ "Sonian",
      stringr::str_detect(tree_id_up, "^(ZO)")       ~ "Sonian",
      stringr::str_detect(tree_id_up, "^(RAJ|RF)")   ~ "Rajhenav",
      stringr::str_detect(tree_id_up, "^(STR|DFST)") ~ "Strødam",
      stringr::str_detect(tree_id_up, "^(SUS|DFSU)") ~ "Suserup",
      TRUE ~ NA_character_
    ),
    country = case_when(
      stringr::str_detect(tree_id_up, "^(ZF|SON|ZO)")   ~ "Belgium",
      stringr::str_detect(tree_id_up, "^(RAJ|RF)")      ~ "Slovenia",
      stringr::str_detect(tree_id_up, "^(STR|DFST)")    ~ "Denmark",
      stringr::str_detect(tree_id_up, "^(SUS|DFSU)")    ~ "Denmark",
      TRUE ~ NA_character_
    )
  )
}

make_sample_id <- function(..., sep = "_") {
  parts <- list(...)
  parts_chr <- lapply(parts, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  out <- purrr::pmap_chr(parts_chr, function(...) {
    vals <- c(...)
    vals <- vals[vals != ""]
    paste(vals, collapse = sep)
  })
  out
}

species_cols_sorted <- function(df, meta_cols) {
  setdiff(names(df), meta_cols) %>% sort()
}

reorder_metadata_first <- function(df, meta_cols) {
  meta_keep <- meta_cols[meta_cols %in% names(df)]
  sp_keep <- setdiff(names(df), meta_keep) %>% sort()
  df %>% select(all_of(meta_keep), all_of(sp_keep))
}

coerce_binary_int <- function(df, meta_cols) {
  sp_cols <- setdiff(names(df), meta_cols)
  df %>%
    mutate(across(all_of(sp_cols), ~ as.integer(replace_na(as.numeric(.x), 0) > 0)))
}

add_missing_species <- function(df, all_species, meta_cols) {
  missing <- setdiff(all_species, names(df))
  if (length(missing) > 0) {
    for (sp in missing) df[[sp]] <- 0L
  }
  reorder_metadata_first(df, meta_cols)
}

coerce_year_factor <- function(df, year_col = "year", levels = year_levels) {
  if (!year_col %in% names(df)) return(df)
  df %>%
    mutate(
      !!year_col := factor(as.character(.data[[year_col]]), levels = levels)
    )
}

harmonize_species_pa <- function(df, tax_map, from_col, to_col) {
  stopifnot("species" %in% names(df))
  
  tax_map2 <- tax_map %>%
    transmute(
      old_name = stringr::str_trim(as.character(.data[[from_col]])),
      new_name = stringr::str_trim(as.character(.data[[to_col]]))
    ) %>%
    filter(
      !is.na(old_name), old_name != "",
      !is.na(new_name), new_name != ""
    ) %>%
    distinct(old_name, .keep_all = TRUE)
  
  df %>%
    mutate(species = stringr::str_trim(as.character(species))) %>%
    left_join(tax_map2, by = c("species" = "old_name")) %>%
    mutate(species = coalesce(new_name, species)) %>%
    select(-new_name) %>%
    group_by(species) %>%
    summarise(
      across(everything(), ~ as.integer(any(as.numeric(.x) > 0, na.rm = TRUE))),
      .groups = "drop"
    )
}

pivot_pa_matrix <- function(df, id_cols, species_col = "species") {
  id_cols <- id_cols[id_cols %in% names(df)]
  
  out <- df %>%
    distinct(across(all_of(c(id_cols, species_col)))) %>%
    mutate(pa = 1L) %>%
    pivot_wider(
      id_cols = all_of(id_cols),
      names_from = all_of(species_col),
      values_from = pa,
      values_fill = 0L
    ) %>%
    coerce_binary_int(meta_cols = id_cols) %>%
    reorder_metadata_first(meta_cols = id_cols)
  
  out
}

report_matrix <- function(df, label, meta_cols, year_col = "year") {
  sp_cols <- setdiff(names(df), meta_cols)
  cat("--------------------------------------------------\n")
  cat(label, "\n")
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
  cat("Species columns:", length(sp_cols), "\n")
  if (!is.na(year_col) && year_col %in% names(df)) {
    yrs <- unique(as.character(df[[year_col]]))
    yrs <- yrs[!is.na(yrs)]
    yrs <- yrs[yrs != ""]
    yrs <- yrs[match(year_levels, yrs, nomatch = 0)]
    cat("Years detected:", paste(yrs, collapse = ", "), "\n")
  }
}

write_matrix <- function(df, filename) {
  df <- coerce_year_factor(df)
  write_csv2(df, file.path(outdir, filename))
}

matrix_bundle_exists <- function(stem, outdir = "basic_data_JHC") {
  needed <- file.path(
    outdir,
    paste0(c(stem, paste0(stem, "_TARGET"), paste0(stem, "_WOODL")), ".csv")
  )
  all(file.exists(needed))
}

report_existing_bundle <- function(stem) {
  cat("--------------------------------------------------\n")
  cat("Existing analysis-ready matrix bundle already present.\n")
  cat("Base name:", stem, "\n")
}

#-----------------------------------------------------------
# 0A. Ecology subset helpers
#     TARGET == 1
#     WOODL  == code == "L"
#-----------------------------------------------------------
read_prepare_eco_target <- function(path = "taxonomy.xlsx", sheet = "eco_target") {
  eco <- read_excel(path, sheet = sheet) %>%
    select(Taxon, code, TARGET) %>%
    mutate(
      Taxon = clean_species_vector(Taxon),
      Taxon = capitalize_genus(Taxon),
      code = clean_species_vector(code),
      TARGET = suppressWarnings(as.integer(as.numeric(TARGET)))
    ) %>%
    filter(!is.na(Taxon), Taxon != "") %>%
    distinct(Taxon, .keep_all = TRUE)
  
  cat("eco_target rows retained after cleaning:", nrow(eco), "\n")
  cat("eco_target target taxa definitions:", sum(eco$TARGET == 1, na.rm = TRUE), "\n")
  cat("eco_target wood taxa definitions (code == 'L'):", sum(eco$code == "L", na.rm = TRUE), "\n")
  
  eco
}

derive_retained_taxa_for_matrix <- function(df, meta_cols, eco_target_tbl, matrix_label) {
  sp_cols <- setdiff(names(df), meta_cols)
  
  match_tbl <- tibble(Taxon = sp_cols) %>%
    left_join(eco_target_tbl, by = "Taxon")
  
  matched_n <- sum(!is.na(match_tbl$code) | !is.na(match_tbl$TARGET))
  unmatched_n <- sum(is.na(match_tbl$code) & is.na(match_tbl$TARGET))
  
  target_taxa <- match_tbl %>%
    filter(TARGET == 1) %>%
    pull(Taxon)
  
  wood_taxa <- match_tbl %>%
    filter(code == "L") %>%
    pull(Taxon)
  
  cat("--------------------------------------------------\n")
  cat("Ecology subset diagnostics:", matrix_label, "\n")
  cat("Species columns in matrix:", length(sp_cols), "\n")
  cat("Matched to eco_target:", matched_n, "\n")
  cat("Unmatched to eco_target:", unmatched_n, "\n")
  cat("TARGET taxa retained:", length(target_taxa), "\n")
  cat("WOODL taxa retained:", length(wood_taxa), "\n")
  
  list(
    matched_n = matched_n,
    unmatched_n = unmatched_n,
    target_taxa = target_taxa,
    wood_taxa = wood_taxa
  )
}

filter_matrix_species_subset <- function(df, meta_cols, retained_taxa) {
  meta_keep <- meta_cols[meta_cols %in% names(df)]
  sp_order <- setdiff(names(df), meta_keep)
  taxa_keep <- sp_order[sp_order %in% retained_taxa]
  
  out <- df %>%
    select(all_of(meta_keep), all_of(taxa_keep))
  
  if (length(taxa_keep) > 0) {
    out <- coerce_binary_int(out, meta_cols = meta_keep)
  }
  
  out
}

write_ecology_subset_versions <- function(df, meta_cols, stem, eco_target_tbl, label = stem) {
  df <- coerce_year_factor(df)
  
  subset_info <- derive_retained_taxa_for_matrix(
    df = df,
    meta_cols = meta_cols,
    eco_target_tbl = eco_target_tbl,
    matrix_label = label
  )
  
  df_target <- filter_matrix_species_subset(
    df = df,
    meta_cols = meta_cols,
    retained_taxa = subset_info$target_taxa
  ) %>%
    coerce_year_factor()
  
  df_wood <- filter_matrix_species_subset(
    df = df,
    meta_cols = meta_cols,
    retained_taxa = subset_info$wood_taxa
  ) %>%
    coerce_year_factor()
  
  write_matrix(df_target, paste0(stem, "_TARGET.csv"))
  write_matrix(df_wood, paste0(stem, "_WOODL.csv"))
  
  cat("Wrote:", paste0(stem, "_TARGET.csv"), "| species =", length(subset_info$target_taxa), "\n")
  cat("Wrote:", paste0(stem, "_WOODL.csv"), "| species =", length(subset_info$wood_taxa), "\n")
  
  invisible(
    list(
      TARGET = df_target,
      WOODL = df_wood,
      diagnostics = subset_info
    )
  )
}

#-----------------------------------------------------------
# 1. Read eco_target once for all later TARGET / WOODL subsets
#-----------------------------------------------------------
eco_target_tbl <- read_prepare_eco_target("taxonomy_revisedJHC.xlsx", "eco_target")

#-----------------------------------------------------------
# 2. Read main identification dataset
#-----------------------------------------------------------
dat_raw <- read_excel("master_identifications_2026.xlsx", sheet = "ID")
str(dat_raw)

tree_col <- get_first_existing_col(dat_raw, c("tree_id", "Tree_ID", "tree", "Tree"))
substrate_col <- get_first_existing_col(dat_raw, c("substrate", "Substrate", "SUBSTRATE"))
species_col <- get_first_existing_col(dat_raw, c("Species", "species", "Taxon"))

if (is.na(tree_col)) stop("No tree_id-like column found in master_identifications_2026.xlsx sheet 'ID'.")
if (is.na(substrate_col)) stop("No substrate column found in master_identifications_2026.xlsx sheet 'ID'.")
if (is.na(species_col)) stop("No species column found in master_identifications_2026.xlsx sheet 'ID'.")

year_col <- get_first_existing_col(dat_raw, c("year", "Year", "YEAR"))

if ("date" %in% names(dat_raw)) {
  if (!inherits(dat_raw$date, c("POSIXct", "POSIXt", "Date"))) {
    dat_raw$date <- suppressWarnings(as.POSIXct(dat_raw$date, tz = "UTC"))
  }
  year_vec <- suppressWarnings(as.integer(format(dat_raw$date, "%Y")))
  year_source <- "derived from date column"
} else if (!is.na(year_col)) {
  year_vec <- suppressWarnings(as.integer(dat_raw[[year_col]]))
  year_source <- paste0("explicit year column: ", year_col)
} else {
  stop("No usable year information: neither 'date' nor 'year' column available.")
}

valid_year_idx <- !is.na(year_vec)
if (!any(valid_year_idx)) {
  stop("Year detection failed: all derived year values are NA.")
}

cat("Year source:", year_source, "\n")
cat("Detected years:", paste(sort(unique(year_vec[valid_year_idx])), collapse = ", "), "\n")

tree_meta <- dat_raw %>%
  transmute(
    tree_id = as.character(.data[[tree_col]])
  ) %>%
  mutate(
    tree_id = stringr::str_trim(tree_id)
  ) %>%
  filter(!is.na(tree_id), tree_id != "") %>%
  distinct(tree_id) %>%
  left_join(
    derive_tree_metadata(.$tree_id) %>% distinct(tree_id, .keep_all = TRUE),
    by = "tree_id"
  )

#-----------------------------------------------------------
# 3. Standardize main identification dataset
#-----------------------------------------------------------
dat <- dat_raw %>%
  mutate(.year_detected = year_vec) %>%
  transmute(
    tree_id   = as.character(.data[[tree_col]]),
    year      = .year_detected,
    substrate = as.character(.data[[substrate_col]]),
    species   = capitalize_genus(as.character(.data[[species_col]]))
  ) %>%
  mutate(
    tree_id   = stringr::str_trim(tree_id),
    substrate = stringr::str_trim(substrate),
    species   = clean_species_vector(species)
  ) %>%
  filter(
    !is.na(tree_id), tree_id != "",
    !is.na(year),
    !is.na(substrate), substrate != "",
    !is.na(species), !species %in% invalid_species_labels
  ) %>%
  left_join(tree_meta, by = "tree_id") %>%
  filter(!(forest == "Sonian" & year == 2022L)) %>%
  mutate(
    sample_id = make_sample_id(tree_id, year, substrate)
  ) %>%
  select(country, forest, tree_id, year, substrate, sample_id, species)

cat("Main identification rows after cleaning:", nrow(dat), "\n")
cat("Main identification unique trees:", dplyr::n_distinct(dat$tree_id), "\n")
cat("Main identification years:", paste(sort(unique(dat$year)), collapse = ", "), "\n")

#-----------------------------------------------------------
# 4. Main matrix cascade
#-----------------------------------------------------------
meta_main_detailed <- c("country", "forest", "tree_id", "year", "substrate", "sample_id")
meta_main_year <- c("country", "forest", "tree_id", "year", "sample_id")
meta_main_tree <- c("country", "forest", "tree_id")

main_tree_year_substrate <- pivot_pa_matrix(
  df = dat,
  id_cols = meta_main_detailed
) %>%
  arrange(country, forest, tree_id, year, substrate, sample_id)

main_tree_year <- dat %>%
  distinct(country, forest, tree_id, year, species) %>%
  mutate(sample_id = make_sample_id(tree_id, year)) %>%
  select(country, forest, tree_id, year, sample_id, species) %>%
  pivot_pa_matrix(id_cols = meta_main_year) %>%
  arrange(country, forest, tree_id, year, sample_id)

main_tree <- dat %>%
  distinct(country, forest, tree_id, species) %>%
  pivot_pa_matrix(id_cols = meta_main_tree) %>%
  arrange(country, forest, tree_id)

#-----------------------------------------------------------
# 5. Save main matrix cascade
#-----------------------------------------------------------
write_matrix(main_tree_year_substrate, "main_pa_tree_year_substrate.csv")
write_matrix(main_tree_year, "main_pa_tree_year.csv")
write_matrix(main_tree, "main_pa_tree.csv")

report_matrix(main_tree_year_substrate, "Main detailed matrix: tree_id x year x substrate x species", meta_main_detailed)
report_matrix(main_tree_year, "Main intermediate matrix: tree_id x year x species", meta_main_year)
report_matrix(main_tree, "Main collapsed matrix: tree_id x species", meta_main_tree, year_col = NA_character_)

#-----------------------------------------------------------
# 6. Write TARGET and WOODL subset versions for main outputs
#-----------------------------------------------------------
main_tree_year_substrate_subsets <- write_ecology_subset_versions(
  df = main_tree_year_substrate,
  meta_cols = meta_main_detailed,
  stem = "main_pa_tree_year_substrate",
  eco_target_tbl = eco_target_tbl,
  label = "Main detailed matrix: tree_id x year x substrate"
)

main_tree_year_subsets <- write_ecology_subset_versions(
  df = main_tree_year,
  meta_cols = meta_main_year,
  stem = "main_pa_tree_year",
  eco_target_tbl = eco_target_tbl,
  label = "Main intermediate matrix: tree_id x year"
)

main_tree_subsets <- write_ecology_subset_versions(
  df = main_tree,
  meta_cols = meta_main_tree,
  stem = "main_pa_tree",
  eco_target_tbl = eco_target_tbl,
  label = "Main collapsed matrix: tree_id"
)

#-----------------------------------------------------------
# 7. Read NATMAN taxonomy mapping
#-----------------------------------------------------------
tax_natman <- read_excel("comm.xlsx", sheet = "tax_natman")

#-----------------------------------------------------------
# 8. Read NATMAN_2001 matrix (species x tree_id)
#-----------------------------------------------------------
natman_raw <- read_excel("comm.xlsx", sheet = "NATMAN_2001")
names(natman_raw)[1] <- "species"

natman_raw <- natman_raw %>%
  mutate(species = clean_species_vector(species)) %>%
  filter(
    !is.na(species),
    !species %in% invalid_species_labels
  )

#-----------------------------------------------------------
# 8A. Remove excluded NATMAN_2001 tree_ids
#-----------------------------------------------------------
drop_natman_tree_ids <- c(
  "SON172", "SON203", "SON204",
  "RAJ111", "RAJ112", "RAJ120"
)

present_drop_natman_tree_ids <- intersect(drop_natman_tree_ids, names(natman_raw))

cat("NATMAN_2001 tree_id columns requested for removal:", length(drop_natman_tree_ids), "\n")
cat("NATMAN_2001 tree_id columns found and removed:", length(present_drop_natman_tree_ids), "\n")
if (length(present_drop_natman_tree_ids) > 0) {
  print(present_drop_natman_tree_ids)
}

natman_raw <- natman_raw %>%
  select(-any_of(drop_natman_tree_ids))

#-----------------------------------------------------------
# 9. Harmonize NATMAN species names
#-----------------------------------------------------------
natman_raw_harmonized <- harmonize_species_pa(
  df = natman_raw,
  tax_map = tax_natman,
  from_col = "natman",
  to_col = "simple_NATMAN"
)

#-----------------------------------------------------------
# 10. Convert NATMAN_2001 to standardized tree_id x species
#-----------------------------------------------------------
natman_mat <- natman_raw_harmonized %>%
  column_to_rownames("species") %>%
  as.data.frame(check.names = FALSE)

natman_tree <- as.data.frame(t(natman_mat), check.names = FALSE) %>%
  rownames_to_column("tree_id") %>%
  as_tibble() %>%
  mutate(
    tree_id = as.character(tree_id),
    across(-tree_id, ~ as.integer(as.numeric(.x) > 0))
  ) %>%
  left_join(derive_tree_metadata(.$tree_id), by = "tree_id") %>%
  mutate(
    year = "2001",
    sample_id = make_sample_id(tree_id, year)
  ) %>%
  reorder_metadata_first(meta_cols = c("country", "forest", "tree_id", "year", "sample_id")) %>%
  arrange(country, forest, tree_id) %>%
  coerce_year_factor()

write_matrix(natman_tree, "NATMAN2001_pa_tree.csv")
report_matrix(natman_tree, "NATMAN_2001 standardized tree-level matrix", c("country", "forest", "tree_id", "year", "sample_id"))

cat("NATMAN species after harmonization:", ncol(natman_tree) - 5, "\n")

#-----------------------------------------------------------
# 11. Write TARGET and WOODL subset versions for NATMAN_2001
#-----------------------------------------------------------
natman_tree_subsets <- write_ecology_subset_versions(
  df = natman_tree,
  meta_cols = c("country", "forest", "tree_id", "year", "sample_id"),
  stem = "NATMAN2001_pa_tree",
  eco_target_tbl = eco_target_tbl,
  label = "NATMAN_2001 standardized tree-level matrix"
)

#-----------------------------------------------------------
# 12. Harmonize tree-level species columns across main and NATMAN
#     Collapsed workflow uses year = 2001 vs 2021-2023.
#-----------------------------------------------------------
meta_tree_std <- c("country", "forest", "tree_id", "year", "sample_id")

main_tree_std <- main_tree %>%
  mutate(
    year = "2021-2023",
    sample_id = make_sample_id(tree_id, year)
  ) %>%
  select(country, forest, tree_id, year, sample_id, everything()) %>%
  coerce_year_factor()

natman_tree_collapsed <- natman_tree %>%
  mutate(
    year = "2001",
    sample_id = make_sample_id(tree_id, year)
  ) %>%
  coerce_year_factor()

all_species_tree <- union(
  setdiff(names(main_tree_std), meta_tree_std),
  setdiff(names(natman_tree_collapsed), meta_tree_std)
) %>%
  sort()

main_tree_harmonized <- add_missing_species(main_tree_std, all_species_tree, meta_tree_std) %>%
  coerce_binary_int(meta_tree_std) %>%
  coerce_year_factor() %>%
  arrange(country, forest, tree_id)

natman_tree_harmonized <- add_missing_species(natman_tree_collapsed, all_species_tree, meta_tree_std) %>%
  coerce_binary_int(meta_tree_std) %>%
  coerce_year_factor() %>%
  arrange(country, forest, tree_id)

write_matrix(main_tree_harmonized, "harmonized_main_pa_tree.csv")
write_matrix(natman_tree_harmonized, "harmonized_NATMAN2001_pa_tree.csv")

cat("Harmonized tree-level species universe:", length(all_species_tree), "\n")
cat("Species in main tree matrix before harmonization:", ncol(main_tree) - length(meta_main_tree), "\n")
cat("Species in NATMAN_2001 before harmonization:", ncol(natman_tree_collapsed) - length(meta_tree_std), "\n")

#-----------------------------------------------------------
# 13. Write TARGET and WOODL subset versions for harmonized
#     tree-level outputs
#-----------------------------------------------------------
main_tree_harmonized_subsets <- write_ecology_subset_versions(
  df = main_tree_harmonized,
  meta_cols = meta_tree_std,
  stem = "harmonized_main_pa_tree",
  eco_target_tbl = eco_target_tbl,
  label = "Harmonized main tree-level matrix"
)

natman_tree_harmonized_subsets <- write_ecology_subset_versions(
  df = natman_tree_harmonized,
  meta_cols = meta_tree_std,
  stem = "harmonized_NATMAN2001_pa_tree",
  eco_target_tbl = eco_target_tbl,
  label = "Harmonized NATMAN_2001 tree-level matrix"
)

#-----------------------------------------------------------
# 14. Optional stacked standardized tree-level dataset
#-----------------------------------------------------------
stacked_tree_harmonized <- bind_rows(
  main_tree_harmonized %>% mutate(source_dataset = "main"),
  natman_tree_harmonized %>% mutate(source_dataset = "NATMAN_2001")
) %>%
  coerce_year_factor() %>%
  select(source_dataset, all_of(meta_tree_std), all_of(all_species_tree)) %>%
  arrange(source_dataset, country, forest, tree_id)

write_matrix(stacked_tree_harmonized, "harmonized_stacked_pa_tree_all_datasets.csv")

#-----------------------------------------------------------
# 15. Write TARGET and WOODL subset versions for stacked output
#-----------------------------------------------------------
stacked_tree_harmonized_subsets <- write_ecology_subset_versions(
  df = stacked_tree_harmonized,
  meta_cols = c("source_dataset", meta_tree_std),
  stem = "harmonized_stacked_pa_tree_all_datasets",
  eco_target_tbl = eco_target_tbl,
  label = "Harmonized stacked tree-level dataset"
)

#-----------------------------------------------------------
# 16A. Main analysis matrix = master_identification + NATMAN_2001
#      Provenance is retained via source_dataset.
#      Collapsed workflow uses year = 2001 vs 2021-2023.
#-----------------------------------------------------------
analysis_stem <- "analysis_main_natman_pa_tree"
meta_analysis <- c("source_dataset", meta_tree_std)

if (matrix_bundle_exists(analysis_stem, outdir = outdir)) {
  report_existing_bundle(analysis_stem)
} else {
  analysis_main_natman_pa_tree <- stacked_tree_harmonized %>%
    coerce_year_factor() %>%
    arrange(source_dataset, country, forest, tree_id, year, sample_id)
  
  write_matrix(analysis_main_natman_pa_tree, paste0(analysis_stem, ".csv"))
  
  report_matrix(
    analysis_main_natman_pa_tree,
    "Main analysis matrix: master_identification + NATMAN_2001 (provenance retained)",
    meta_cols = meta_analysis
  )
  
  analysis_main_natman_pa_tree_subsets <- write_ecology_subset_versions(
    df = analysis_main_natman_pa_tree,
    meta_cols = meta_analysis,
    stem = analysis_stem,
    eco_target_tbl = eco_target_tbl,
    label = "Main analysis matrix: master_identification + NATMAN_2001"
  )
  
  cat("Main analysis matrix created from stacked dataset.\n")
  cat("Base name:", analysis_stem, "\n")
  cat("Rows from main dataset:", sum(analysis_main_natman_pa_tree$source_dataset == "main"), "\n")
  cat("Rows from NATMAN_2001:", sum(analysis_main_natman_pa_tree$source_dataset == "NATMAN_2001"), "\n")
}

#-----------------------------------------------------------
# 16B. Year-resolved main analysis matrix = main_tree_year + NATMAN_2001
#      Provenance retained via source_dataset.
#      Year-resolved workflow retains 2001, 2021, 2022, 2023.
#-----------------------------------------------------------
analysis_year_stem <- "analysis_main_natman_pa_tree_year"
meta_analysis_year <- c("source_dataset", meta_main_year)

if (matrix_bundle_exists(analysis_year_stem, outdir = outdir)) {
  report_existing_bundle(analysis_year_stem)
} else {
  main_tree_year_for_analysis <- main_tree_year %>%
    mutate(
      year = as.character(year),
      sample_id = make_sample_id(tree_id, year)
    ) %>%
    coerce_year_factor()
  
  natman_tree_year_for_analysis <- natman_tree %>%
    mutate(
      year = "2001",
      sample_id = make_sample_id(tree_id, year)
    ) %>%
    coerce_year_factor()
  
  all_species_tree_year <- union(
    setdiff(names(main_tree_year_for_analysis), meta_main_year),
    setdiff(names(natman_tree_year_for_analysis), meta_main_year)
  ) %>%
    sort()
  
  main_tree_year_harmonized_for_analysis <- add_missing_species(
    main_tree_year_for_analysis,
    all_species_tree_year,
    meta_main_year
  ) %>%
    coerce_binary_int(meta_main_year) %>%
    coerce_year_factor() %>%
    arrange(country, forest, tree_id, year, sample_id)
  
  natman_tree_harmonized_for_analysis_year <- add_missing_species(
    natman_tree_year_for_analysis,
    all_species_tree_year,
    meta_main_year
  ) %>%
    coerce_binary_int(meta_main_year) %>%
    coerce_year_factor() %>%
    arrange(country, forest, tree_id, year, sample_id)
  
  analysis_main_natman_pa_tree_year <- bind_rows(
    main_tree_year_harmonized_for_analysis %>% mutate(source_dataset = "main"),
    natman_tree_harmonized_for_analysis_year %>% mutate(source_dataset = "NATMAN_2001")
  ) %>%
    coerce_year_factor() %>%
    select(source_dataset, all_of(meta_main_year), all_of(all_species_tree_year)) %>%
    arrange(source_dataset, country, forest, tree_id, year, sample_id)
  
  write_matrix(
    analysis_main_natman_pa_tree_year,
    paste0(analysis_year_stem, ".csv")
  )
  
  report_matrix(
    analysis_main_natman_pa_tree_year,
    "Year-resolved main analysis matrix: main_tree_year + NATMAN_2001 (provenance retained)",
    meta_cols = meta_analysis_year
  )
  
  analysis_main_natman_pa_tree_year_subsets <- write_ecology_subset_versions(
    df = analysis_main_natman_pa_tree_year,
    meta_cols = meta_analysis_year,
    stem = analysis_year_stem,
    eco_target_tbl = eco_target_tbl,
    label = "Year-resolved main analysis matrix: main_tree_year + NATMAN_2001"
  )
  
  cat("Year-resolved analysis matrix created from main_tree_year + NATMAN_2001.\n")
  cat("Base name:", analysis_year_stem, "\n")
  cat("Rows from main dataset:", sum(analysis_main_natman_pa_tree_year$source_dataset == "main"), "\n")
  cat("Rows from NATMAN_2001:", sum(analysis_main_natman_pa_tree_year$source_dataset == "NATMAN_2001"), "\n")
}

#-----------------------------------------------------------
# 17. Final reporting
#-----------------------------------------------------------
report_matrix(main_tree_harmonized, "Final harmonized main tree matrix", meta_tree_std)
report_matrix(natman_tree_harmonized, "Final harmonized NATMAN tree matrix", meta_tree_std)
report_matrix(stacked_tree_harmonized, "Final stacked harmonized tree-level dataset", c("source_dataset", meta_tree_std))

if (exists("analysis_main_natman_pa_tree")) {
  report_matrix(
    analysis_main_natman_pa_tree,
    "Final main analysis matrix (main + NATMAN, provenance retained)",
    meta_analysis
  )
}

if (exists("analysis_main_natman_pa_tree_year")) {
  report_matrix(
    analysis_main_natman_pa_tree_year,
    "Final year-resolved main analysis matrix (main_tree_year + NATMAN, provenance retained)",
    meta_analysis_year
  )
}

cat("--------------------------------------------------\n")
cat("Output directory:", outdir, "\n")
cat("Note: substrate-level outputs were generated only for the main dataset.\n")
cat("NATMAN_2001 was not expanded to substrate level because that information is not present in the source matrix.\n")
cat("All full files were written with metadata columns first and binary integer species columns afterward.\n")
cat("Year was coerced consistently to factor levels:", paste(year_levels, collapse = ", "), "\n")
cat("Year-resolved outputs retain survey years 2001, 2021, 2022 and 2023.\n")
cat("Collapsed tree-level outputs use 2001 vs 2021-2023, with the collapsed label propagated into sample_id.\n")
cat("Additional ecology subset outputs were also written for TARGET taxa (TARGET == 1) and WOODL taxa (code == 'L').\n")
cat("Unmatched taxa in eco_target were retained in full matrices but excluded from subset matrices.\n")
cat("If a subset matrix contained zero retained taxa, the corresponding file was still written with metadata columns only.\n")
cat("New analysis bundle basename:", analysis_stem, "\n")
cat("This analysis bundle retains provenance through source_dataset and combines the main dataset with NATMAN_2001.\n")
cat("New year-resolved analysis bundle basename:", analysis_year_stem, "\n")
cat("This year-resolved analysis bundle retains provenance through source_dataset and combines main_tree_year with NATMAN_2001.\n")