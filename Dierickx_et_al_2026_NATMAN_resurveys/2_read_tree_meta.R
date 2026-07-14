suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
})

#-----------------------------------------------------------
# 1. Read data
#-----------------------------------------------------------
tree_meta_raw <- read_excel(
  "master_identifications_2026.xlsx",
  sheet = "tree_meta"
)

cat("Raw tree_meta:", nrow(tree_meta_raw), "rows\n")

#-----------------------------------------------------------
# 2. Standardize names
#-----------------------------------------------------------
tree_meta <- tree_meta_raw %>%
  rename_with(~ tolower(.x)) %>%
  rename(tree_id = natman_id)

#-----------------------------------------------------------
# 3. Helper functions
#-----------------------------------------------------------
clean_chr <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_squish() %>%
    na_if("") %>%
    na_if("NA") %>%
    na_if("na")
}

clean_num <- function(x) {
  x <- clean_chr(x)
  x <- gsub(",", ".", x, fixed = TRUE)
  as.numeric(x)
}

clean_pct <- function(x) {
  x <- clean_num(x)
  ifelse(is.na(x), NA_real_, pmin(pmax(x, 0), 100))
}

clean_binary01 <- function(x) {
  x <- clean_chr(x)
  case_when(
    is.na(x) ~ NA_integer_,
    x %in% c("1", "yes", "true") ~ 1L,
    x %in% c("0", "no", "false") ~ 0L,
    TRUE ~ suppressWarnings(as.integer(x))
  )
}

clean_lower <- function(x) {
  x %>%
    clean_chr() %>%
    stringr::str_to_lower()
}

#-----------------------------------------------------------
# 4. Apply cleaning
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    tree_id = dplyr::recode(tree_id, "ZFid1035" = "ZFid1035_1037")
  )
tree_meta <- tree_meta %>%
  mutate(across(everything(), ~ if (is.character(.x)) clean_chr(.x) else .x))

#-----------------------------------------------------------
# 5. Directional decay scale (0–6)
#-----------------------------------------------------------
ds_vars <- c(
  "log_av_ds", "log_min_ds", "log_max_ds",
  "crown_av_ds", "cr_min_ds", "cr_max_ds"
)

tree_meta <- tree_meta %>%
  mutate(across(any_of(ds_vars), ~ {
    x <- clean_num(.x)
    x <- pmin(pmax(x, 0), 6)
    x
  }))

#-----------------------------------------------------------
# 6. Continuous structural variables
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    log_diameter   = clean_num(log_diameter),    # meters
    log_og_diameter= clean_num(log_og_diameter),
    length_m       = clean_num(length_m),
    vol_est        = clean_num(vol_est)          # m³
  )

#-----------------------------------------------------------
# 7. Percentages (0–100)
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    mosscover   = clean_pct(mosscover),
    barkcover   = clean_pct(barkcover),
    soil_contact= clean_pct(soil_contact)
  )

#-----------------------------------------------------------
# 8. Angular variable (0–360°)
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    azimuth = clean_num(azimuth),
    azimuth = ifelse(is.na(azimuth), NA, azimuth %% 360)
  )

#-----------------------------------------------------------
# 9. Binary indicators
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    snag  = clean_binary01(snag),
    crown = clean_binary01(crown)
  )

#-----------------------------------------------------------
# 10. Categorical variables
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    type = case_when(
      clean_lower(type) %in% c("fallen") ~ "fallen",
      clean_lower(type) %in% c("broken") ~ "broken",
      TRUE ~ clean_lower(type)
    ),
    
    exposure = case_when(
      clean_lower(exposure) %in% c("sun") ~ "sun",
      clean_lower(exposure) %in% c("halfshade", "half-shade") ~ "halfshade",
      clean_lower(exposure) %in% c("shade") ~ "shade",
      TRUE ~ NA_character_
    ),
    
    son_location = case_when(
      clean_lower(son_location) %in% c("kern") ~ "core",
      clean_lower(son_location) %in% c("extern") ~ "edge",
      TRUE ~ NA_character_
    ),
    
    log_size   = factor(clean_lower(log_size),   levels = c("small", "medium", "large"), ordered = TRUE),
    crown_size = factor(clean_lower(crown_size), levels = c("small", "medium", "large"), ordered = TRUE),
    
    on_slope = case_when(
      clean_lower(on_slope) %in% c("no") ~ "no",
      clean_lower(on_slope) %in% c("partial") ~ "partial",
      clean_lower(on_slope) %in% c("yes") ~ "yes",
      TRUE ~ NA_character_
    )
  )

#-----------------------------------------------------------
# 11. Convert to factors
#-----------------------------------------------------------
tree_meta <- tree_meta %>%
  mutate(
    type        = factor(type, levels = c("fallen", "broken")),
    exposure    = factor(exposure, levels = c("shade", "halfshade", "sun"), ordered = TRUE),
    on_slope    = factor(on_slope, levels = c("no", "partial", "yes"), ordered = TRUE),
    son_location= factor(son_location),
    country     = factor(country),
    forest      = factor(forest),
    period      = factor(clean_lower(period), levels = c("previous", "current"))
  )

#-----------------------------------------------------------
# 12. Final deduplication
#-----------------------------------------------------------
tree_meta_clean <- tree_meta %>%
  distinct(tree_id, .keep_all = TRUE)

#-----------------------------------------------------------
# 13. Write output
#-----------------------------------------------------------
write_csv2(tree_meta_clean, "basic_data_JHC/tree_meta_clean.csv")

# check of the tree_ids are the same ones as in the community data
comm_check <- read_csv2("basic_data_JHC/analysis_main_natman_pa_tree_year_TARGET.csv")
length(unique(comm_check$tree_id))
length(tree_meta_clean$tree_id)
setdiff((unique(comm_check$tree_id)), (tree_meta_clean$tree_id))
str(tree_meta)

summarise_numeric <- function(x) {
  x_non_na <- x[!is.na(x)]
  n <- length(x_non_na)
  
  if (n == 0) {
    return(tibble(
      min = NA_real_,
      mean = NA_real_,
      lwr = NA_real_,
      upr = NA_real_,
      max = NA_real_
    ))
  }
  
  m <- mean(x_non_na)
  s <- sd(x_non_na)
  se <- s / sqrt(n)
  
  tibble(
    min = min(x_non_na),
    mean = m,
    lwr = m - 1.96 * se,
    upr = m + 1.96 * se,
    max = max(x_non_na)
  )
}

is_binary_numeric <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  vals <- sort(unique(x[!is.na(x)]))
  length(vals) > 0 && all(vals %in% c(0, 1))
}

summarise_binary <- function(x) {
  x_non_na <- x[!is.na(x)]
  n <- length(x_non_na)
  
  if (n == 0) return(NA_character_)
  
  n0 <- sum(x_non_na == 0)
  n1 <- sum(x_non_na == 1)
  
  p0 <- 100 * n0 / n
  p1 <- 100 * n1 / n
  
  sprintf("0 = %d (%.1f%%), 1 = %d (%.1f%%)", n0, p0, n1, p1)
}

summarise_categorical <- function(x, max_levels = 10) {
  x_non_na <- x[!is.na(x)]
  
  if (length(x_non_na) == 0) return(NA_character_)
  
  tab <- sort(table(x_non_na), decreasing = TRUE)
  
  if (length(tab) > max_levels) {
    tab_show <- tab[1:max_levels]
    extra <- length(tab) - max_levels
    out <- paste0(names(tab_show), " (", as.integer(tab_show), ")", collapse = ", ")
    out <- paste0(out, ", ... +", extra, " more")
  } else {
    out <- paste0(names(tab), " (", as.integer(tab), ")", collapse = ", ")
  }
  
  out
}

known_binary_vars <- c("snag", "crown")
summary_tbl <- tree_meta %>%
  summarise(across(everything(), ~ list(.x))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "data"
  ) %>%
  mutate(
    n_na = map_int(data, ~ sum(is.na(.x))),
    is_numeric = map_lgl(data, is.numeric),
    is_categorical = map_lgl(data, ~ is.factor(.x) || is.character(.x) || is.ordered(.x)),
    is_binary = map2_lgl(variable, data, ~ .x %in% known_binary_vars || is_binary_numeric(.y))
  ) %>%
  mutate(
    numeric_stats = map2(data, is_numeric & !is_binary, ~ {
      if (.y) {
        summarise_numeric(.x)
      } else {
        tibble(
          min = NA_real_,
          mean = NA_real_,
          lwr = NA_real_,
          upr = NA_real_,
          max = NA_real_
        )
      }
    }),
    categorical_levels = map2_chr(data, is_categorical, ~ {
      if (.y) summarise_categorical(.x) else NA_character_
    }),
    binary_levels = map2_chr(data, is_binary, ~ {
      if (.y) summarise_binary(.x) else NA_character_
    })
  ) %>%
  unnest(numeric_stats) %>%
  mutate(
    range_ci = case_when(
      is_binary ~ binary_levels,
      is_numeric ~ sprintf("(%.3f) %.3f [%.3f, %.3f] (%.3f)", min, mean, lwr, upr, max),
      is_categorical ~ categorical_levels,
      TRUE ~ NA_character_
    )
  ) %>%
  select(variable, range_ci, n_na)
print(summary_tbl, n = Inf)
