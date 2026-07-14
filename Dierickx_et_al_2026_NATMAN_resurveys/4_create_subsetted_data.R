suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})

set.seed(42)

#-----------------------------------------------------------
# 1. User options
#-----------------------------------------------------------
#comm_file_year <- "basic_data_JHC_LOGCR_only/analysis_main_natman_pa_tree_year_TARGET.csv"
#comm_file_collapsed <- "basic_data_JHC_LOGCR_only/analysis_main_natman_pa_tree_TARGET.csv"
#subset_base_dir <- "SUBSET_BY_COUNTRY_DIAMMATCH_logcronly"
#outdir <- "subset_multivariate_inputs_logcronly"
comm_file_year <- "basic_data_JHC/analysis_main_natman_pa_tree_year_TARGET.csv"
comm_file_collapsed <- "basic_data_JHC/analysis_main_natman_pa_tree_TARGET.csv"
subset_base_dir <- "SUBSET_BY_COUNTRY_DIAMMATCH"
outdir <- "subset_multivariate_inputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

min_species_per_sample <- 2L

analysis_settings <- tribble(
  ~setting_name,        ~current_sampling_mode, ~current_year_keep,
  "single_year_2022",   "single_year",          2022L,
  "single_year_2023",   "single_year",          2023L,
  "collapsed",          "collapsed",            NA_integer_,
  "year_retained",      "year_retained",        NA_integer_
)

country_lookup <- tribble(
  ~country_name, ~country_dir,
  "Belgium",     "belgium",
  "Denmark",     "denmark",
  "Slovenia",    "slovenia"
)

#-----------------------------------------------------------
# 2. Helpers
#-----------------------------------------------------------
standardize_forest <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_squish(x)
  
  case_when(
    is.na(x) ~ NA_character_,
    x %in% c("Strodam", "Strødam") ~ "Strødam",
    TRUE ~ x
  )
}

derive_period <- function(df) {
  if ("period" %in% names(df)) {
    p <- as.character(df$period)
    
  } else if ("timespan" %in% names(df)) {
    p <- as.character(df$timespan)
    
  } else if ("year" %in% names(df)) {
    y <- as.character(df$year)
    y <- str_trim(y)
    y <- str_squish(y)
    
    p <- case_when(
      y %in% c("2001", "previous", "prev", "historic", "historical") ~ "previous",
      y %in% c("2021", "2022", "2023", "2021-2023", "current", "present") ~ "current",
      TRUE ~ NA_character_
    )
    
  } else {
    stop("Could not derive period. No period, timespan, or year column found.")
  }
  
  p <- str_trim(p)
  p <- str_to_lower(p)
  p <- recode(
    p,
    previous = "previous",
    prev = "previous",
    historic = "previous",
    historical = "previous",
    current = "current",
    present = "current"
  )
  
  factor(p, levels = c("previous", "current"))
}

filter_by_min_species <- function(df, taxa_names, min_species_per_sample, sample_id_col = "sample_id") {
  X <- df %>%
    select(all_of(taxa_names)) %>%
    as.matrix()
  
  storage.mode(X) <- "numeric"
  X[is.na(X)] <- 0
  
  richness <- rowSums(X > 0)
  keep <- richness >= min_species_per_sample
  
  if (any(!keep)) {
    cat("Removing rows with <", min_species_per_sample, "species:", sum(!keep), "\n")
    if (sample_id_col %in% names(df)) {
      cat("Dropped sample IDs:\n")
      print(df[[sample_id_col]][!keep])
      cat("\n")
    }
  }
  
  df[keep, , drop = FALSE]
}

filter_unit_comm <- function(df, country_name,
                             current_sampling_mode = "single_year",
                             current_year_keep = 2022L) {
  out <- df %>%
    filter(country == country_name)
  
  if (current_sampling_mode == "single_year") {
    if (!"year" %in% names(out)) {
      stop("current_sampling_mode = 'single_year' requires a year column in the community matrix.")
    }
    
    year_chr <- as.character(out$year)
    year_chr <- str_trim(year_chr)
    
    out <- out %>%
      mutate(.year_chr = year_chr) %>%
      filter(
        period == "previous" |
          country == "Belgium" |
          (
            country %in% c("Denmark", "Slovenia") &
              (period != "current" | .year_chr == as.character(current_year_keep))
          )
      ) %>%
      select(-.year_chr)
  }
  
  out
}

prepare_subset_comm <- function(subset_df, comm_df, taxa_names, min_species_per_sample = 2L) {
  meta_sub <- subset_df %>%
    select(any_of(c(
      "tree_id", "period", "forest", "country",
      "selection_type", "decay_group", "candidate_rank",
      "son_location", "core_priority"
    ))) %>%
    distinct()
  
  join_keys <- intersect(c("tree_id", "period", "forest", "country"), names(meta_sub))
  
  out <- comm_df %>%
    inner_join(meta_sub, by = join_keys)
  
  X <- out %>%
    select(all_of(taxa_names)) %>%
    as.matrix()
  
  storage.mode(X) <- "numeric"
  X[is.na(X)] <- 0
  
  richness <- rowSums(X > 0)
  keep <- richness >= min_species_per_sample
  
  if (any(!keep)) {
    cat("Removing subset rows with <", min_species_per_sample, "species:", sum(!keep), "\n")
    if ("sample_id" %in% names(out)) {
      cat("Dropped sample IDs:\n")
      print(out$sample_id[!keep])
      cat("\n")
    }
  }
  
  out[keep, , drop = FALSE]
}

read_country_subset <- function(country_name, subset_base_dir) {
  country_dir <- file.path(
    subset_base_dir,
    str_to_lower(country_name)
  )
  
  prev_det_file <- file.path(country_dir, "previous_subset_diammatch.csv")
  curr_det_file <- file.path(country_dir, "current_subset_diammatch.csv")
  
  if (!file.exists(prev_det_file) || !file.exists(curr_det_file)) {
    stop("Missing deterministic subset files for country: ", country_name)
  }
  
  prev_det <- read_csv2(prev_det_file, show_col_types = FALSE) %>%
    mutate(
      country = as.character(country),
      forest = standardize_forest(forest),
      tree_id = as.character(tree_id),
      period = "previous"
    )
  
  curr_det <- read_csv2(curr_det_file, show_col_types = FALSE) %>%
    mutate(
      country = as.character(country),
      forest = standardize_forest(forest),
      tree_id = as.character(tree_id),
      period = "current"
    )
  
  bind_rows(prev_det, curr_det)
}

#-----------------------------------------------------------
# 3. Main builder per setting
#-----------------------------------------------------------
build_multivariate_input <- function(setting_name, current_sampling_mode, current_year_keep) {
  
  cat("\n####################################################\n")
  cat("BUILDING DATASET:", setting_name, "\n")
  cat("####################################################\n")
  
  comm_file <- if (current_sampling_mode == "collapsed") {
    comm_file_collapsed
  } else {
    comm_file_year
  }
  
  dat <- read_csv2(comm_file, show_col_types = FALSE)
  
  meta_cols_base <- c("source_dataset", "country", "forest", "tree_id", "year", "sample_id", "timespan", "period")
  meta_cols_present <- intersect(meta_cols_base, names(dat))
  
  dat <- dat %>%
    mutate(
      country = as.character(country),
      forest = standardize_forest(forest),
      tree_id = as.character(tree_id),
      year = if ("year" %in% names(.)) as.character(year) else NA_character_,
      sample_id = if ("sample_id" %in% names(.)) as.character(sample_id) else as.character(tree_id),
      period = derive_period(.)
    ) %>%
    filter(!is.na(period))
  
  taxa_cols <- setdiff(names(dat), union(meta_cols_present, "period"))
  
  cat("\nInput file:", comm_file, "\n")
  cat("Rows before global richness filter:", nrow(dat), "\n")
  cat("Taxa:", length(taxa_cols), "\n")
  
  dat <- filter_by_min_species(
    df = dat,
    taxa_names = taxa_cols,
    min_species_per_sample = min_species_per_sample,
    sample_id_col = "sample_id"
  )
  
  cat("Rows after global richness filter:", nrow(dat), "\n")
  cat("Country x forest x period after global filter:\n")
  print(dat %>% count(country, forest, period), n = Inf)
  cat("\n")
  
  combined_list <- vector("list", nrow(country_lookup))
  
  current_year_keep_out <- if (current_sampling_mode == "single_year") current_year_keep else NA_integer_
  
  for (i in seq_len(nrow(country_lookup))) {
    country_name <- country_lookup$country_name[i]
    
    cat("----------------------------------------------------\n")
    cat("Country:", country_name, "\n")
    
    dat_country <- filter_unit_comm(
      df = dat,
      country_name = country_name,
      current_sampling_mode = current_sampling_mode,
      current_year_keep = current_year_keep
    )
    
    if (nrow(dat_country) == 0) {
      cat("No community rows for this country after year filtering. Skipping.\n\n")
      next
    }
    
    det_subset <- read_country_subset(
      country_name = country_name,
      subset_base_dir = subset_base_dir
    )
    
    if (current_sampling_mode == "single_year" && country_name %in% c("Denmark", "Slovenia")) {
      current_tree_ids_present <- dat_country %>%
        filter(period == "current") %>%
        distinct(tree_id) %>%
        pull(tree_id)
      
      det_subset <- det_subset %>%
        filter(
          period == "previous" |
            (period == "current" & tree_id %in% current_tree_ids_present)
        )
    }
    
    joined_country <- prepare_subset_comm(
      subset_df = det_subset,
      comm_df = dat_country,
      taxa_names = taxa_cols,
      min_species_per_sample = min_species_per_sample
    ) %>%
      mutate(
        setting = setting_name,
        current_sampling_mode = current_sampling_mode,
        current_year_keep = current_year_keep_out,
        min_species_per_sample = min_species_per_sample,
        timespan = if ("timespan" %in% names(.)) as.character(timespan) else as.character(period)
      )
    
    cat("Rows retained for", country_name, ":", nrow(joined_country), "\n")
    cat("Forest x period counts:\n")
    print(joined_country %>% count(forest, period), n = Inf)
    
    if ("year" %in% names(joined_country)) {
      cat("Year counts:\n")
      print(table(joined_country$year, useNA = "ifany"))
    }
    cat("\n")
    
    combined_list[[i]] <- joined_country
  }
  
  combined_dat <- bind_rows(combined_list)
  
  if (nrow(combined_dat) == 0) {
    stop("Combined deterministic matched dataset is empty for setting: ", setting_name)
  }
  
  meta_first <- c(
    "setting", "current_sampling_mode", "current_year_keep", "min_species_per_sample",
    "source_dataset", "country", "forest", "tree_id", "year", "sample_id", "timespan", "period",
    "selection_type", "decay_group", "candidate_rank", "son_location", "core_priority"
  )
  
  meta_first_present <- intersect(meta_first, names(combined_dat))
  taxa_present <- setdiff(names(combined_dat), meta_first_present)
  
  combined_dat <- combined_dat %>%
    select(all_of(meta_first_present), all_of(taxa_present))
  
  cat("\n====================\n")
  cat("FINAL COMBINED DATASET\n")
  cat("====================\n")
  cat("Setting:", setting_name, "\n")
  cat("Rows:", nrow(combined_dat), "\n")
  cat("Taxa:", length(taxa_present), "\n")
  cat("Country x forest x period:\n")
  print(combined_dat %>% count(country, forest, period), n = Inf)
  cat("\n")
  
  outfile_main <- file.path(
    outdir,
    paste0("analysis_input_deterministic_", setting_name, ".csv")
  )
  
  write_csv2(combined_dat, outfile_main)
  
  summary_tab <- combined_dat %>%
    count(country, forest, period, year, name = "n_samples") %>%
    arrange(country, forest, period, year)
  
  write_csv2(
    summary_tab,
    file.path(outdir, paste0("analysis_input_deterministic_", setting_name, "_summary.csv"))
  )
  
  invisible(combined_dat)
}

#-----------------------------------------------------------
# 4. Run all settings
#-----------------------------------------------------------
all_outputs <- pmap(
  analysis_settings,
  function(setting_name, current_sampling_mode, current_year_keep) {
    build_multivariate_input(
      setting_name = setting_name,
      current_sampling_mode = current_sampling_mode,
      current_year_keep = current_year_keep
    )
  }
)

