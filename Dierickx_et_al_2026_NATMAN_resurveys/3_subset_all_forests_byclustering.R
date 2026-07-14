suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

set.seed(42)

outdir_main <- "SUBSET_BY_COUNTRY_DIAMMATCH"
#outdir_main <- "SUBSET_BY_COUNTRY_DIAMMATCH_logcronly"

dir.create(outdir_main, showWarnings = FALSE, recursive = TRUE)

#-----------------------------------------------------------
# user options
#-----------------------------------------------------------
allow_difference <- 2L
n_random_reps <- 100L

#-----------------------------------------------------------
# helper: harmonize forest labels
#-----------------------------------------------------------
standardize_forest <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_squish(x)
  
  case_when(
    x %in% c("Strødam", "Strodam") ~ "Strødam",
    TRUE ~ x
  )
}

#-----------------------------------------------------------
# samples to remove globally
#-----------------------------------------------------------
drop_tree_ids <- c("SON121", "RAJ043")

#-----------------------------------------------------------
# 0. Input community matrix
#-----------------------------------------------------------
#comm_file <- "basic_data_JHC_LOGCR_only/analysis_main_natman_pa_tree_year_TARGET.csv"
comm_file <- "basic_data_JHC/analysis_main_natman_pa_tree_year_TARGET.csv"

comm_dat <- read_csv2(comm_file, show_col_types = FALSE) %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id)
  ) %>%
  filter(!tree_id %in% drop_tree_ids)

meta_cols <- c("source_dataset", "country", "forest", "tree_id", "year", "sample_id")
meta_cols_present <- intersect(meta_cols, names(comm_dat))
taxa_cols <- setdiff(names(comm_dat), meta_cols_present)

#-----------------------------------------------------------
# 1. Identify non-empty tree_id x period combinations
#-----------------------------------------------------------
nonempty_tree_period <- comm_dat %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id),
    period = case_when(
      year == 2001 ~ "previous",
      year %in% c(2021, 2022, 2023) ~ "current",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period)) %>%
  mutate(
    richness_target = rowSums(across(all_of(taxa_cols)), na.rm = TRUE)
  ) %>%
  filter(richness_target > 1) %>%
  distinct(country, forest, tree_id, period)

cat("\n====================\n")
cat("NON-EMPTY COMMUNITY FILTER\n")
cat("====================\n")
cat("Globally removed tree IDs:\n")
print(drop_tree_ids)
cat("\n")
print(nonempty_tree_period %>% count(country, forest, period), n = Inf)
cat("\n")
#-----------------------------------------------------------
# 2. Helper functions
#-----------------------------------------------------------
safe_pct <- function(num, den) {
  ifelse(den > 0, 100 * num / den, NA_real_)
}

sample_random_n <- function(df, n_take) {
  if (n_take <= 0 || nrow(df) == 0) return(df[0, , drop = FALSE])
  df %>% slice_sample(n = min(n_take, nrow(df)))
}

compute_flexible_targets <- function(n_previous, n_current, allow_difference = 2L) {
  n_previous <- as.integer(n_previous)
  n_current <- as.integer(n_current)
  allow_difference <- as.integer(allow_difference)
  
  base_n <- min(n_previous, n_current)
  
  if (n_previous > n_current) {
    target_previous <- min(n_previous, base_n + allow_difference)
    target_current  <- base_n
  } else if (n_current > n_previous) {
    target_previous <- base_n
    target_current  <- min(n_current, base_n + allow_difference)
  } else {
    target_previous <- base_n
    target_current  <- base_n
  }
  
  tibble(
    n_previous_target = target_previous,
    n_current_target = target_current
  )
}

make_cluster_plot <- function(plot_df, country_name) {
  plot_df <- plot_df %>%
    mutate(
      period = factor(period, levels = c("previous", "current")),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6"))
    )
  
  ggplot(
    plot_df,
    aes(x = decay_group, y = log_og_diameter, fill = period)
  ) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.35,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_jitter(
      aes(color = period),
      width = 0.15,
      height = 0,
      alpha = 0.65,
      size = 2
    ) +
    facet_wrap(~ forest, scales = "free_x") +
    labs(
      title = paste("Candidate pool structure by forest:", country_name),
      subtitle = "Decay group versus log original diameter before subset selection",
      x = "Decay group",
      y = "log_og_diameter"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

make_inclusion_plot <- function(plot_df, country_name) {
  plot_df <- plot_df %>%
    mutate(
      period = factor(period, levels = c("previous", "current")),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
      inclusion_status = factor(inclusion_status, levels = c("excluded", "included"))
    )
  
  ggplot(
    plot_df,
    aes(
      x = decay_group,
      y = log_og_diameter,
      color = period,
      shape = inclusion_status
    )
  ) +
    geom_jitter(
      width = 0.15,
      height = 0,
      alpha = 0.8,
      size = 2.4
    ) +
    facet_wrap(~ forest, scales = "free_x") +
    labs(
      title = paste("Final inclusion by forest:", country_name),
      subtitle = "Included versus excluded trees after diameter-aware matching",
      x = "Decay group",
      y = "log_og_diameter",
      color = "Period",
      shape = "Selection"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

make_cluster_count_plot <- function(plot_df, country_name) {
  count_df <- plot_df %>%
    count(forest, period, decay_group, name = "n")
  
  ggplot(
    count_df,
    aes(x = decay_group, y = n, fill = period)
  ) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~ forest, scales = "free_y") +
    labs(
      title = paste("Forest × decay-group counts:", country_name),
      subtitle = "Available candidate trees before subset selection",
      x = "Decay group",
      y = "Number of trees"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

#-----------------------------------------------------------
# Matching helpers
#-----------------------------------------------------------

build_greedy_pairs <- function(prev_df, curr_df, randomize_ties = FALSE) {
  if (nrow(prev_df) == 0 || nrow(curr_df) == 0) {
    return(tibble(
      prev_row = integer(),
      curr_row = integer(),
      abs_diff = numeric()
    ))
  }
  
  prev_used <- rep(FALSE, nrow(prev_df))
  curr_used <- rep(FALSE, nrow(curr_df))
  
  pair_list <- vector("list", min(nrow(prev_df), nrow(curr_df)))
  k <- 0L
  
  repeat {
    avail_prev <- which(!prev_used)
    avail_curr <- which(!curr_used)
    
    if (length(avail_prev) == 0 || length(avail_curr) == 0) break
    
    dist_grid <- expand.grid(
      prev_row = avail_prev,
      curr_row = avail_curr,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      mutate(
        abs_diff = abs(
          prev_df$log_og_diameter[prev_row] - curr_df$log_og_diameter[curr_row]
        ),
        core_score = prev_df$core_priority[prev_row] + curr_df$core_priority[curr_row],
        tie_rand = if (randomize_ties) runif(n()) else 0
      ) %>%
      arrange(abs_diff, core_score, tie_rand, prev_row, curr_row)
    
    best <- dist_grid %>% slice(1)
    
    prev_used[best$prev_row] <- TRUE
    curr_used[best$curr_row] <- TRUE
    
    k <- k + 1L
    pair_list[[k]] <- best %>%
      select(prev_row, curr_row, abs_diff)
  }
  
  bind_rows(pair_list)
}

select_matched_subset <- function(prev_df, curr_df,
                                  n_prev_target, n_curr_target,
                                  prioritize_core = FALSE) {
  
  n_prev_target <- as.integer(n_prev_target)
  n_curr_target <- as.integer(n_curr_target)
  
  prev_df <- prev_df %>%
    mutate(.row_prev = row_number())
  
  curr_df <- curr_df %>%
    mutate(.row_curr = row_number())
  
  #---------------------------------------------------------
  # 1. Build greedy diameter-based pairs
  #---------------------------------------------------------
  pair_tbl <- build_greedy_pairs(prev_df, curr_df, randomize_ties = FALSE)
  
  n_pairs_keep <- min(n_prev_target, n_curr_target)
  if (n_pairs_keep > nrow(pair_tbl)) n_pairs_keep <- nrow(pair_tbl)
  
  kept_pairs <- pair_tbl %>%
    slice_head(n = n_pairs_keep)
  
  prev_keep_rows <- kept_pairs$prev_row
  curr_keep_rows <- kept_pairs$curr_row
  
  prev_subset <- prev_df %>%
    slice(prev_keep_rows) %>%
    mutate(selection_type = "paired")
  
  curr_subset <- curr_df %>%
    slice(curr_keep_rows) %>%
    mutate(selection_type = "paired")
  
  #---------------------------------------------------------
  # 2. Remaining candidates
  #---------------------------------------------------------
  prev_remaining <- prev_df %>%
    filter(!.row_prev %in% prev_keep_rows)
  
  curr_remaining <- curr_df %>%
    filter(!.row_curr %in% curr_keep_rows)
  
  n_prev_extra <- n_prev_target - nrow(prev_subset)
  n_curr_extra <- n_curr_target - nrow(curr_subset)
  
  #---------------------------------------------------------
  # 3. Add extra trees (larger side only)
  #---------------------------------------------------------
  if (n_prev_extra > 0 && nrow(prev_remaining) > 0) {
    ref_vals <- curr_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- curr_df$log_og_diameter
    
    n_prev_remaining <- nrow(prev_remaining)
    n_prev_take <- min(n_prev_extra, n_prev_remaining)
    
    prev_extra <- prev_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L
      ) %>%
      arrange(dist_to_ref, core_rank, tree_id) %>%
      slice_head(n = n_prev_take)
    
    prev_subset <- bind_rows(
      prev_subset,
      prev_extra %>% mutate(selection_type = "extra")
    )
  }
  
  if (n_curr_extra > 0 && nrow(curr_remaining) > 0) {
    ref_vals <- prev_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- prev_df$log_og_diameter
    
    n_curr_remaining <- nrow(curr_remaining)
    n_curr_take <- min(n_curr_extra, n_curr_remaining)
    
    curr_extra <- curr_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L
      ) %>%
      arrange(dist_to_ref, core_rank, tree_id) %>%
      slice_head(n = n_curr_take)
    
    curr_subset <- bind_rows(
      curr_subset,
      curr_extra %>% mutate(selection_type = "extra")
    )
  }
  
  list(
    previous = prev_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank"))),
    current = curr_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank"))),
    pairs = kept_pairs
  )
}

select_random_matched_subset <- function(prev_df, curr_df,
                                         n_prev_target, n_curr_target,
                                         prioritize_core = FALSE) {
  
  n_prev_target <- as.integer(n_prev_target)
  n_curr_target <- as.integer(n_curr_target)
  
  prev_df <- prev_df %>%
    slice_sample(n = nrow(prev_df)) %>%
    mutate(.row_prev = row_number())
  
  curr_df <- curr_df %>%
    slice_sample(n = nrow(curr_df)) %>%
    mutate(.row_curr = row_number())
  
  pair_tbl <- build_greedy_pairs(prev_df, curr_df, randomize_ties = TRUE)
  
  n_pairs_keep <- min(n_prev_target, n_curr_target)
  if (n_pairs_keep > nrow(pair_tbl)) n_pairs_keep <- nrow(pair_tbl)
  
  kept_pairs <- pair_tbl %>%
    slice_head(n = n_pairs_keep)
  
  prev_keep_rows <- kept_pairs$prev_row
  curr_keep_rows <- kept_pairs$curr_row
  
  prev_subset <- prev_df %>%
    slice(prev_keep_rows) %>%
    mutate(selection_type = "paired")
  
  curr_subset <- curr_df %>%
    slice(curr_keep_rows) %>%
    mutate(selection_type = "paired")
  
  prev_remaining <- prev_df %>%
    filter(!.row_prev %in% prev_keep_rows)
  
  curr_remaining <- curr_df %>%
    filter(!.row_curr %in% curr_keep_rows)
  
  n_prev_extra <- n_prev_target - nrow(prev_subset)
  n_curr_extra <- n_curr_target - nrow(curr_subset)
  
  if (n_prev_extra > 0 && nrow(prev_remaining) > 0) {
    ref_vals <- curr_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- curr_df$log_og_diameter
    
    n_prev_remaining <- nrow(prev_remaining)
    n_prev_take <- min(n_prev_extra, n_prev_remaining)
    
    prev_extra <- prev_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L,
        rand_key = runif(n())
      ) %>%
      arrange(dist_to_ref, core_rank, rand_key) %>%
      slice_head(n = n_prev_take)
    
    prev_subset <- bind_rows(
      prev_subset,
      prev_extra %>% mutate(selection_type = "extra")
    )
  }
  
  if (n_curr_extra > 0 && nrow(curr_remaining) > 0) {
    ref_vals <- prev_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- prev_df$log_og_diameter
    
    n_curr_remaining <- nrow(curr_remaining)
    n_curr_take <- min(n_curr_extra, n_curr_remaining)
    
    curr_extra <- curr_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L,
        rand_key = runif(n())
      ) %>%
      arrange(dist_to_ref, core_rank, rand_key) %>%
      slice_head(n = n_curr_take)
    
    curr_subset <- bind_rows(
      curr_subset,
      curr_extra %>% mutate(selection_type = "extra")
    )
  }
  
  list(
    previous = prev_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank", "rand_key"))),
    current = curr_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank", "rand_key")))
  )
}

generate_random_matched_subsets <- function(
    candidates_previous,
    candidates_current,
    target_plan,
    strata_vars,
    n_reps = 100L,
    seed = NULL,
    prioritize_core = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  target_strata <- target_plan %>%
    filter(n_previous_target > 0 | n_current_target > 0) %>%
    select(all_of(strata_vars), n_previous_target, n_current_target)
  
  out_prev <- vector("list", n_reps)
  out_curr <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    prev_list <- vector("list", nrow(target_strata))
    curr_list <- vector("list", nrow(target_strata))
    
    for (j in seq_len(nrow(target_strata))) {
      key_row <- target_strata[j, , drop = FALSE]
      
      prev_pool <- candidates_previous
      curr_pool <- candidates_current
      
      for (v in strata_vars) {
        prev_pool <- prev_pool %>% filter(.data[[v]] == key_row[[v]][1])
        curr_pool <- curr_pool %>% filter(.data[[v]] == key_row[[v]][1])
      }
      
      sel_j <- select_random_matched_subset(
        prev_df = prev_pool,
        curr_df = curr_pool,
        n_prev_target = key_row$n_previous_target[1],
        n_curr_target = key_row$n_current_target[1],
        prioritize_core = prioritize_core
      )
      
      prev_list[[j]] <- sel_j$previous
      curr_list[[j]] <- sel_j$current
    }
    
    out_prev[[i]] <- bind_rows(prev_list) %>%
      mutate(replicate_id = i)
    
    out_curr[[i]] <- bind_rows(curr_list) %>%
      mutate(replicate_id = i)
  }
  
  list(
    previous = bind_rows(out_prev),
    current = bind_rows(out_curr)
  )
}

#-----------------------------------------------------------
# 3. Main per-country function
#-----------------------------------------------------------
run_country_subset <- function(country_name,
                               tree_meta_clean,
                               nonempty_tree_period,
                               allow_difference = 2L,
                               n_random_reps = 100L) {
  
  country_name <- as.character(country_name)
  
  cat("\n====================================================\n")
  cat("PROCESSING COUNTRY:", country_name, "\n")
  cat("====================================================\n")
  
  outdir <- file.path(
    outdir_main,
    stringr::str_replace_all(stringr::str_to_lower(country_name), "[^a-z0-9]+", "_")
  )
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  meta_country_raw <- tree_meta_clean %>%
    mutate(
      country = as.character(country),
      forest = standardize_forest(forest),
      tree_id = as.character(tree_id),
      period = as.character(period),
      son_location = if ("son_location" %in% names(.)) as.character(son_location) else NA_character_
    ) %>%
    filter(country == country_name) %>%
    select(any_of(c(
      "country", "forest", "tree_id", "period",
      "log_av_ds", "log_og_diameter", "son_location"
    )))
  
  meta_country <- meta_country_raw %>%
    inner_join(
      nonempty_tree_period %>%
        mutate(
          country = as.character(country),
          forest = standardize_forest(forest),
          tree_id = as.character(tree_id),
          period = as.character(period)
        ) %>%
        filter(country == country_name),
      by = c("country", "forest", "tree_id", "period")
    )
  
  cat("\nRows after join with non-empty community filter:\n")
  print(meta_country %>% count(forest, period), n = Inf)
  cat("\n")
  
  if (nrow(meta_country) == 0) {
    cat("No rows available after filtering. Skipping.\n")
    return(NULL)
  }
  
  if ("son_location" %in% names(meta_country) && country_name == "Belgium") {
    meta_country <- meta_country %>%
      mutate(
        core_priority = case_when(
          forest == "Sonian" & son_location == "core" ~ 1L,
          forest == "Sonian" & son_location == "edge" ~ 2L,
          TRUE ~ 3L
        )
      )
  } else {
    meta_country <- meta_country %>%
      mutate(core_priority = 1L)
  }
  
  meta_country <- meta_country %>%
    mutate(
      log_av_ds = as.integer(round(log_av_ds)),
      decay_group = case_when(
        log_av_ds == 1 ~ "1",
        log_av_ds == 2 ~ "2",
        log_av_ds == 3 ~ "3",
        log_av_ds == 4 ~ "4",
        log_av_ds %in% 5:6 ~ "5_6",
        TRUE ~ NA_character_
      ),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6"))
    )
  
  if (country_name == "Denmark") {
    strata_vars <- c("forest", "decay_group")
    cat("Matching: forest x decay_group + diameter optimization\n\n")
  } else {
    strata_vars <- c("decay_group")
    cat("Matching: decay_group + diameter optimization\n\n")
  }
  
  current_dist <- meta_country %>%
    filter(period == "current") %>%
    count(across(all_of(strata_vars)), name = "n_current")
  
  previous_dist <- meta_country %>%
    filter(period == "previous") %>%
    count(across(all_of(strata_vars)), name = "n_previous")
  
  target_plan <- full_join(current_dist, previous_dist, by = strata_vars) %>%
    mutate(
      n_current = replace_na(n_current, 0L),
      n_previous = replace_na(n_previous, 0L)
    ) %>%
    bind_cols(
      purrr::pmap_dfr(
        list(.$n_previous, .$n_current),
        ~ compute_flexible_targets(..1, ..2, allow_difference)
      )
    ) %>%
    arrange(across(all_of(strata_vars)))
  
  print(target_plan, n = Inf)
  
  candidates <- meta_country %>%
    filter(if_all(all_of(strata_vars), ~ !is.na(.x))) %>%
    arrange(period, across(all_of(strata_vars)), core_priority, tree_id) %>%
    group_by(period, across(all_of(strata_vars))) %>%
    mutate(candidate_rank = row_number()) %>%
    ungroup()
  
  candidates_previous <- candidates %>% filter(period == "previous")
  candidates_current  <- candidates %>% filter(period == "current")
  
  target_rows <- split(target_plan, seq_len(nrow(target_plan)))
  
  subset_prev_list <- vector("list", length(target_rows))
  subset_curr_list <- vector("list", length(target_rows))
  pair_diag_list   <- vector("list", length(target_rows))
  
  for (i in seq_along(target_rows)) {
    row_i <- target_rows[[i]]
    
    prev_pool <- candidates_previous
    curr_pool <- candidates_current
    
    for (v in strata_vars) {
      prev_pool <- prev_pool %>% filter(.data[[v]] == row_i[[v]][1])
      curr_pool <- curr_pool %>% filter(.data[[v]] == row_i[[v]][1])
    }
    
    sel_i <- select_matched_subset(
      prev_df = prev_pool,
      curr_df = curr_pool,
      n_prev_target = row_i$n_previous_target[1],
      n_curr_target = row_i$n_current_target[1],
      prioritize_core = (country_name == "Belgium")
    )
    
    subset_prev_list[[i]] <- sel_i$previous
    subset_curr_list[[i]] <- sel_i$current
    
    if (nrow(sel_i$pairs) > 0) {
      pair_diag_list[[i]] <- bind_cols(
        row_i[rep(1, nrow(sel_i$pairs)), strata_vars, drop = FALSE],
        sel_i$pairs %>%
          mutate(
            prev_tree_id = prev_pool$tree_id[prev_row],
            curr_tree_id = curr_pool$tree_id[curr_row],
            prev_log_og_diameter = prev_pool$log_og_diameter[prev_row],
            curr_log_og_diameter = curr_pool$log_og_diameter[curr_row]
          )
      )
    } else {
      pair_diag_list[[i]] <- tibble()
    }
  }
  
  prev_subset <- bind_rows(subset_prev_list)
  curr_subset <- bind_rows(subset_curr_list)
  pair_diagnostic <- bind_rows(pair_diag_list)
  
  comparison_table <- target_plan %>%
    left_join(
      prev_subset %>% count(across(all_of(strata_vars)), name = "n_previous_subset"),
      by = strata_vars
    ) %>%
    left_join(
      curr_subset %>% count(across(all_of(strata_vars)), name = "n_current_subset"),
      by = strata_vars
    ) %>%
    mutate(
      n_previous_subset = replace_na(n_previous_subset, 0L),
      n_current_subset = replace_na(n_current_subset, 0L)
    ) %>%
    arrange(across(all_of(strata_vars)))
  
  random_subsets <- generate_random_matched_subsets(
    candidates_previous = candidates_previous,
    candidates_current = candidates_current,
    target_plan = target_plan,
    strata_vars = strata_vars,
    n_reps = n_random_reps,
    seed = 42,
    prioritize_core = (country_name == "Belgium")
  )
  
  all_decay_original <- meta_country %>%
    count(period, decay_group) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = 0)
  
  all_decay_subset <- bind_rows(prev_subset, curr_subset) %>%
    count(period, decay_group) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = 0)
  
  summary_decay <- full_join(
    all_decay_original,
    all_decay_subset,
    by = "decay_group",
    suffix = c("", "_subset")
  ) %>%
    mutate(
      previous = replace_na(previous, 0L),
      current = replace_na(current, 0L),
      previous_subset = replace_na(previous_subset, 0L),
      current_subset = replace_na(current_subset, 0L),
      pct_previous_retained = safe_pct(previous_subset, previous),
      pct_current_retained = safe_pct(current_subset, current),
      country = country_name
    ) %>%
    select(
      country, decay_group,
      previous, current,
      previous_subset, current_subset,
      pct_previous_retained, pct_current_retained
    ) %>%
    arrange(decay_group)
  
  print(summary_decay, n = Inf)
  
  if (country_name == "Denmark") {
    forest_decay_summary <- comparison_table %>%
      group_by(forest, decay_group) %>%
      summarise(
        n_previous = sum(n_previous, na.rm = TRUE),
        n_current = sum(n_current, na.rm = TRUE),
        n_previous_target = sum(n_previous_target, na.rm = TRUE),
        n_current_target = sum(n_current_target, na.rm = TRUE),
        n_previous_subset = sum(n_previous_subset, na.rm = TRUE),
        n_current_subset = sum(n_current_subset, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(forest, decay_group)
    
    cat("\n====================\n")
    cat("DENMARK FOREST x DECAY DIAGNOSTIC\n")
    cat("====================\n")
    print(forest_decay_summary, n = Inf)
    cat("\n")
  } else {
    forest_decay_summary <- NULL
  }
  
  diameter_summary <- bind_rows(
    prev_subset %>% mutate(period = "previous_subset"),
    curr_subset %>% mutate(period = "current_subset")
  ) %>%
    group_by(period, decay_group) %>%
    summarise(
      n = n(),
      mean_log_og_diameter = mean(log_og_diameter, na.rm = TRUE),
      median_log_og_diameter = median(log_og_diameter, na.rm = TRUE),
      sd_log_og_diameter = sd(log_og_diameter, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(country = country_name) %>%
    select(
      country, period, decay_group, n,
      mean_log_og_diameter, median_log_og_diameter, sd_log_og_diameter
    )
  
  print(diameter_summary, n = Inf)
  
  included_keys <- bind_rows(
    prev_subset %>% select(country, forest, tree_id, period),
    curr_subset %>% select(country, forest, tree_id, period)
  ) %>%
    distinct() %>%
    mutate(inclusion_status = "included")
  
  plot_df_inclusion <- meta_country %>%
    left_join(
      included_keys,
      by = c("country", "forest", "tree_id", "period")
    ) %>%
    mutate(
      inclusion_status = if_else(is.na(inclusion_status), "excluded", inclusion_status)
    )
  
  p_cluster <- make_cluster_plot(meta_country, country_name)
  print(p_cluster)
  
  p_cluster_count <- make_cluster_count_plot(meta_country, country_name)
  print(p_cluster_count)
  
  p_inclusion <- make_inclusion_plot(plot_df_inclusion, country_name)
  print(p_inclusion)
  
  ggsave(
    filename = file.path(outdir, "cluster_structure_plot.png"),
    plot = p_cluster,
    width = 11,
    height = 6.5,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir, "cluster_count_plot.png"),
    plot = p_cluster_count,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(outdir, "inclusion_plot.png"),
    plot = p_inclusion,
    width = 11,
    height = 6.5,
    dpi = 300
  )
  
  write_csv2(target_plan, file.path(outdir, "target_plan_diammatch.csv"))
  write_csv2(comparison_table, file.path(outdir, "comparison_by_matching_stratum.csv"))
  write_csv2(prev_subset, file.path(outdir, "previous_subset_diammatch.csv"))
  write_csv2(curr_subset, file.path(outdir, "current_subset_diammatch.csv"))
  write_csv2(summary_decay, file.path(outdir, "country_decay_summary.csv"))
  write_csv2(diameter_summary, file.path(outdir, "diameter_summary.csv"))
  write_csv2(pair_diagnostic, file.path(outdir, "pair_diagnostic.csv"))
  write_csv2(plot_df_inclusion, file.path(outdir, "plot_inclusion_data.csv"))
  write_csv2(random_subsets$previous, file.path(outdir, "previous_random_diammatch_subsets.csv"))
  write_csv2(random_subsets$current, file.path(outdir, "current_random_diammatch_subsets.csv"))
  
  if (!is.null(forest_decay_summary)) {
    write_csv2(
      forest_decay_summary,
      file.path(outdir, "denmark_forest_decay_diagnostic.csv")
    )
  }
  
  list(
    country = country_name,
    strata_vars = strata_vars,
    target_plan = target_plan,
    comparison_table = comparison_table,
    previous_subset = prev_subset,
    current_subset = curr_subset,
    random_previous_subset = random_subsets$previous,
    random_current_subset = random_subsets$current,
    summary_decay = summary_decay,
    diameter_summary = diameter_summary,
    forest_decay_summary = forest_decay_summary,
    plot_cluster = p_cluster,
    plot_cluster_count = p_cluster_count,
    plot_inclusion = p_inclusion
  )
}

#-----------------------------------------------------------
# 4. Run all countries
#-----------------------------------------------------------
tree_meta_clean_std <- tree_meta_clean %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id),
    period = as.character(period)
  )

country_vec <- tree_meta_clean_std %>%
  filter(!is.na(country)) %>%
  distinct(country) %>%
  pull(country) %>%
  sort()

print(country_vec)

results <- map(
  country_vec,
  ~ run_country_subset(
    country_name = .x,
    tree_meta_clean = tree_meta_clean_std,
    nonempty_tree_period = nonempty_tree_period,
    allow_difference = allow_difference,
    n_random_reps = n_random_reps
  )
)

results <- results[!vapply(results, is.null, logical(1))]

#-----------------------------------------------------------
# 5. Combined country-level decay summary
#-----------------------------------------------------------
all_decay <- bind_rows(map(results, "summary_decay")) %>%
  arrange(country, decay_group)

cat("\n====================\n")
cat("COMBINED COUNTRY-LEVEL DECAY SUMMARY\n")
cat("====================\n")
print(all_decay, n = Inf)
cat("\n")

write_csv2(all_decay, file.path(outdir_main, "combined_country_decay_summary.csv"))

#-----------------------------------------------------------
# 6. Combined retained diameter summary
#-----------------------------------------------------------
all_diameter_summary <- bind_rows(map(results, "diameter_summary")) %>%
  arrange(country, period, decay_group)

cat("\n====================\n")
cat("COMBINED RETAINED DIAMETER SUMMARY\n")
cat("====================\n")
print(all_diameter_summary, n = Inf)
cat("\n")

write_csv2(
  all_diameter_summary,
  file.path(outdir_main, "combined_retained_subset_diameter_summary.csv")
)

#-----------------------------------------------------------
# 7. Optional combined Denmark forest x decay diagnostic
#-----------------------------------------------------------
all_denmark_forest_decay <- bind_rows(
  map(results, "forest_decay_summary")
)

if (nrow(all_denmark_forest_decay) > 0) {
  cat("\n====================\n")
  cat("COMBINED DENMARK FOREST x DECAY DIAGNOSTIC\n")
  cat("====================\n")
  print(all_denmark_forest_decay, n = Inf)
  cat("\n")
  
  write_csv2(
    all_denmark_forest_decay,
    file.path(outdir_main, "combined_denmark_forest_decay_diagnostic.csv")
  )
}



# all_diameter_summary is assumed to already exist
# expected columns:
# country, period, decay_group, n, mean_log_og_diameter, median_log_og_diameter, sd_log_og_diameter

#-----------------------------------------------------------
# 1. Prepare difference table
#-----------------------------------------------------------
diam_diff <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = factor(period, levels = c("previous_subset", "current_subset"))
  ) %>%
  select(country, period, decay_group, n, mean_log_og_diameter) %>%
  pivot_wider(
    names_from = period,
    values_from = c(n, mean_log_og_diameter),
    names_sep = "__"
  ) %>%
  mutate(
    diff_mean_log_og_diameter =
      mean_log_og_diameter__current_subset - mean_log_og_diameter__previous_subset
  ) %>%
  arrange(country, decay_group)

cat("\nDifference table:\n")
print(diam_diff, n = Inf)

#-----------------------------------------------------------
# 2. Figure: current - previous by decay class and country
#-----------------------------------------------------------
p_diff <- ggplot(
  diam_diff,
  aes(x = decay_group, y = diff_mean_log_og_diameter, group = country)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2) +
  geom_segment(
    aes(xend = decay_group, y = 0, yend = diff_mean_log_og_diameter),
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(aes(shape = country), size = 3) +
  facet_wrap(~ country, nrow = 1) +
  labs(
    x = "Decay group",
    y = expression(Delta * " mean log original diameter (current - previous)"),
    title = "Difference in average diameter between current and previous subsets",
    subtitle = "Shown separately for each country across decay classes"
  ) +
  theme_bw()

print(p_diff)

p_means <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = factor(period, levels = c("previous_subset", "current_subset"))
  ) %>%
  ggplot(aes(x = period, y = mean_log_og_diameter, group = decay_group)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  facet_grid(country ~ decay_group) +
  labs(
    x = NULL,
    y = "Mean log original diameter",
    title = "Average diameter in previous and current subsets",
    subtitle = "Faceted by country and decay group"
  ) +
  theme_bw()

print(p_means)

p_means_sd <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = factor(period, levels = c("previous_subset", "current_subset"))
  ) %>%
  ggplot(aes(x = period, y = mean_log_og_diameter, group = decay_group)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_errorbar(
    aes(
      ymin = mean_log_og_diameter - sd_log_og_diameter,
      ymax = mean_log_og_diameter + sd_log_og_diameter
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  facet_grid(country ~ decay_group) +
  labs(
    x = NULL,
    y = "Mean log original diameter ± SD",
    title = "Average diameter in previous and current subsets",
    subtitle = "Panels show country × decay group; bars represent ±1 SD"
  ) +
  theme_bw()

print(p_means_sd)
##################################
##################################
###### standardized graph #####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

set.seed(42)

#-----------------------------------------------------------
# user options
#-----------------------------------------------------------
allow_difference <- 2L
n_random_reps <- 100L

#-----------------------------------------------------------
# plotting options
#-----------------------------------------------------------
selection_cols <- c(
  included = "#CC6B1F",  # burnt orange
  excluded = "#C9B6E4"   # light purple
)

period_shapes <- c(
  previous = 21,         # circle
  current  = 24          # triangle
)

#-----------------------------------------------------------
# helper: harmonize forest labels
#-----------------------------------------------------------
standardize_forest <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_squish(x)
  
  case_when(
    x %in% c("Strødam", "Strodam") ~ "Strødam",
    TRUE ~ x
  )
}

#-----------------------------------------------------------
# samples to remove globally
#-----------------------------------------------------------
drop_tree_ids <- c("SON121", "RAJ043")

#-----------------------------------------------------------
# 0. Input community matrix
#-----------------------------------------------------------
#comm_file <- "basic_data_JHC_LOGCR_only/analysis_main_natman_pa_tree_year_TARGET.csv"
comm_file <- "basic_data_JHC/analysis_main_natman_pa_tree_year_TARGET.csv"

comm_dat <- read_csv2(comm_file, show_col_types = FALSE) %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id)
  ) %>%
  filter(!tree_id %in% drop_tree_ids)

meta_cols <- c("source_dataset", "country", "forest", "tree_id", "year", "sample_id")
meta_cols_present <- intersect(meta_cols, names(comm_dat))
taxa_cols <- setdiff(names(comm_dat), meta_cols_present)

#-----------------------------------------------------------
# 1. Identify non-empty tree_id x period combinations
#-----------------------------------------------------------
nonempty_tree_period <- comm_dat %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id),
    period = case_when(
      year == 2001 ~ "previous",
      year %in% c(2021, 2022, 2023) ~ "current",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period)) %>%
  mutate(
    richness_target = rowSums(across(all_of(taxa_cols)), na.rm = TRUE)
  ) %>%
  filter(richness_target > 1) %>%
  distinct(country, forest, tree_id, period)

cat("\n====================\n")
cat("NON-EMPTY COMMUNITY FILTER\n")
cat("====================\n")
cat("Globally removed tree IDs:\n")
print(drop_tree_ids)
cat("\n")
print(nonempty_tree_period %>% count(country, forest, period), n = Inf)
cat("\n")

#-----------------------------------------------------------
# 2. Helper functions
#-----------------------------------------------------------
safe_pct <- function(num, den) {
  ifelse(den > 0, 100 * num / den, NA_real_)
}

sample_random_n <- function(df, n_take) {
  if (n_take <= 0 || nrow(df) == 0) return(df[0, , drop = FALSE])
  df %>% slice_sample(n = min(n_take, nrow(df)))
}

compute_flexible_targets <- function(n_previous, n_current, allow_difference = 2L) {
  n_previous <- as.integer(n_previous)
  n_current <- as.integer(n_current)
  allow_difference <- as.integer(allow_difference)
  
  base_n <- min(n_previous, n_current)
  
  if (n_previous > n_current) {
    target_previous <- min(n_previous, base_n + allow_difference)
    target_current  <- base_n
  } else if (n_current > n_previous) {
    target_previous <- base_n
    target_current  <- min(n_current, base_n + allow_difference)
  } else {
    target_previous <- base_n
    target_current  <- base_n
  }
  
  tibble(
    n_previous_target = target_previous,
    n_current_target = target_current
  )
}

make_cluster_plot <- function(plot_df, country_name) {
  plot_df <- plot_df %>%
    mutate(
      period = factor(period, levels = c("previous", "current")),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6"))
    )
  
  ggplot(
    plot_df,
    aes(
      x = decay_group,
      y = log_og_diameter,
      fill = period,
      shape = period
    )
  ) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.25,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_jitter(
      aes(colour = period),
      width = 0.15,
      height = 0,
      alpha = 0.65,
      size = 2,
      stroke = 0.35
    ) +
    scale_shape_manual(
      values = period_shapes,
      name = "Survey period",
      labels = c(
        previous = "Previous",
        current = "Current"
      )
    ) +
    facet_wrap(~ forest, scales = "free_x") +
    labs(
      title = paste("Candidate pool structure by forest:", country_name),
      subtitle = "Decay group versus log original diameter before subset selection",
      x = "Decay group",
      y = "log_og_diameter",
      fill = "Survey period",
      colour = "Survey period"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

make_inclusion_plot <- function(plot_df, country_name) {
  plot_df <- plot_df %>%
    mutate(
      period = factor(period, levels = c("previous", "current")),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
      inclusion_status = factor(inclusion_status, levels = c("included", "excluded"))
    )
  
  ggplot(
    plot_df,
    aes(
      x = decay_group,
      y = log_og_diameter,
      fill = inclusion_status,
      colour = inclusion_status,
      shape = period
    )
  ) +
    geom_jitter(
      width = 0.15,
      height = 0,
      alpha = 0.85,
      size = 2.6,
      stroke = 0.35
    ) +
    facet_wrap(~ forest, scales = "free_x") +
    scale_fill_manual(
      values = selection_cols,
      name = "Selection",
      labels = c(
        included = "Included",
        excluded = "Excluded"
      )
    ) +
    scale_colour_manual(
      values = selection_cols,
      name = "Selection",
      labels = c(
        included = "Included",
        excluded = "Excluded"
      )
    ) +
    scale_shape_manual(
      values = period_shapes,
      name = "Survey period",
      labels = c(
        previous = "Previous",
        current = "Current"
      )
    ) +
    labs(
      title = paste("Final inclusion by forest:", country_name),
      subtitle = "Included versus excluded trees after diameter-aware matching",
      x = "Decay group",
      y = "log_og_diameter"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

make_cluster_count_plot <- function(plot_df, country_name) {
  count_df <- plot_df %>%
    mutate(
      period = factor(period, levels = c("previous", "current")),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6"))
    ) %>%
    count(forest, period, decay_group, name = "n")
  
  ggplot(
    count_df,
    aes(x = decay_group, y = n, fill = period)
  ) +
    geom_col(position = "dodge", width = 0.7) +
    facet_wrap(~ forest, scales = "free_y") +
    labs(
      title = paste("Forest × decay-group counts:", country_name),
      subtitle = "Available candidate trees before subset selection",
      x = "Decay group",
      y = "Number of trees",
      fill = "Survey period"
    ) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      legend.position = "bottom"
    )
}

#-----------------------------------------------------------
# Matching helpers
#-----------------------------------------------------------
build_greedy_pairs <- function(prev_df, curr_df, randomize_ties = FALSE) {
  if (nrow(prev_df) == 0 || nrow(curr_df) == 0) {
    return(tibble(
      prev_row = integer(),
      curr_row = integer(),
      abs_diff = numeric()
    ))
  }
  
  prev_used <- rep(FALSE, nrow(prev_df))
  curr_used <- rep(FALSE, nrow(curr_df))
  
  pair_list <- vector("list", min(nrow(prev_df), nrow(curr_df)))
  k <- 0L
  
  repeat {
    avail_prev <- which(!prev_used)
    avail_curr <- which(!curr_used)
    
    if (length(avail_prev) == 0 || length(avail_curr) == 0) break
    
    dist_grid <- expand.grid(
      prev_row = avail_prev,
      curr_row = avail_curr,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      as_tibble() %>%
      mutate(
        abs_diff = abs(
          prev_df$log_og_diameter[prev_row] - curr_df$log_og_diameter[curr_row]
        ),
        core_score = prev_df$core_priority[prev_row] + curr_df$core_priority[curr_row],
        tie_rand = if (randomize_ties) runif(n()) else 0
      ) %>%
      arrange(abs_diff, core_score, tie_rand, prev_row, curr_row)
    
    best <- dist_grid %>% slice(1)
    
    prev_used[best$prev_row] <- TRUE
    curr_used[best$curr_row] <- TRUE
    
    k <- k + 1L
    pair_list[[k]] <- best %>%
      select(prev_row, curr_row, abs_diff)
  }
  
  bind_rows(pair_list)
}

select_matched_subset <- function(prev_df, curr_df,
                                  n_prev_target, n_curr_target,
                                  prioritize_core = FALSE) {
  
  n_prev_target <- as.integer(n_prev_target)
  n_curr_target <- as.integer(n_curr_target)
  
  prev_df <- prev_df %>%
    mutate(.row_prev = row_number())
  
  curr_df <- curr_df %>%
    mutate(.row_curr = row_number())
  
  pair_tbl <- build_greedy_pairs(prev_df, curr_df, randomize_ties = FALSE)
  
  n_pairs_keep <- min(n_prev_target, n_curr_target)
  if (n_pairs_keep > nrow(pair_tbl)) n_pairs_keep <- nrow(pair_tbl)
  
  kept_pairs <- pair_tbl %>%
    slice_head(n = n_pairs_keep)
  
  prev_keep_rows <- kept_pairs$prev_row
  curr_keep_rows <- kept_pairs$curr_row
  
  prev_subset <- prev_df %>%
    slice(prev_keep_rows) %>%
    mutate(selection_type = "paired")
  
  curr_subset <- curr_df %>%
    slice(curr_keep_rows) %>%
    mutate(selection_type = "paired")
  
  prev_remaining <- prev_df %>%
    filter(!.row_prev %in% prev_keep_rows)
  
  curr_remaining <- curr_df %>%
    filter(!.row_curr %in% curr_keep_rows)
  
  n_prev_extra <- n_prev_target - nrow(prev_subset)
  n_curr_extra <- n_curr_target - nrow(curr_subset)
  
  if (n_prev_extra > 0 && nrow(prev_remaining) > 0) {
    ref_vals <- curr_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- curr_df$log_og_diameter
    
    n_prev_take <- min(n_prev_extra, nrow(prev_remaining))
    
    prev_extra <- prev_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L
      ) %>%
      arrange(dist_to_ref, core_rank, tree_id) %>%
      slice_head(n = n_prev_take)
    
    prev_subset <- bind_rows(
      prev_subset,
      prev_extra %>% mutate(selection_type = "extra")
    )
  }
  
  if (n_curr_extra > 0 && nrow(curr_remaining) > 0) {
    ref_vals <- prev_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- prev_df$log_og_diameter
    
    n_curr_take <- min(n_curr_extra, nrow(curr_remaining))
    
    curr_extra <- curr_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L
      ) %>%
      arrange(dist_to_ref, core_rank, tree_id) %>%
      slice_head(n = n_curr_take)
    
    curr_subset <- bind_rows(
      curr_subset,
      curr_extra %>% mutate(selection_type = "extra")
    )
  }
  
  list(
    previous = prev_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank"))),
    current = curr_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank"))),
    pairs = kept_pairs
  )
}

select_random_matched_subset <- function(prev_df, curr_df,
                                         n_prev_target, n_curr_target,
                                         prioritize_core = FALSE) {
  
  n_prev_target <- as.integer(n_prev_target)
  n_curr_target <- as.integer(n_curr_target)
  
  prev_df <- prev_df %>%
    slice_sample(n = nrow(prev_df)) %>%
    mutate(.row_prev = row_number())
  
  curr_df <- curr_df %>%
    slice_sample(n = nrow(curr_df)) %>%
    mutate(.row_curr = row_number())
  
  pair_tbl <- build_greedy_pairs(prev_df, curr_df, randomize_ties = TRUE)
  
  n_pairs_keep <- min(n_prev_target, n_curr_target)
  if (n_pairs_keep > nrow(pair_tbl)) n_pairs_keep <- nrow(pair_tbl)
  
  kept_pairs <- pair_tbl %>%
    slice_head(n = n_pairs_keep)
  
  prev_keep_rows <- kept_pairs$prev_row
  curr_keep_rows <- kept_pairs$curr_row
  
  prev_subset <- prev_df %>%
    slice(prev_keep_rows) %>%
    mutate(selection_type = "paired")
  
  curr_subset <- curr_df %>%
    slice(curr_keep_rows) %>%
    mutate(selection_type = "paired")
  
  prev_remaining <- prev_df %>%
    filter(!.row_prev %in% prev_keep_rows)
  
  curr_remaining <- curr_df %>%
    filter(!.row_curr %in% curr_keep_rows)
  
  n_prev_extra <- n_prev_target - nrow(prev_subset)
  n_curr_extra <- n_curr_target - nrow(curr_subset)
  
  if (n_prev_extra > 0 && nrow(prev_remaining) > 0) {
    ref_vals <- curr_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- curr_df$log_og_diameter
    
    n_prev_take <- min(n_prev_extra, nrow(prev_remaining))
    
    prev_extra <- prev_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L,
        rand_key = runif(n())
      ) %>%
      arrange(dist_to_ref, core_rank, rand_key) %>%
      slice_head(n = n_prev_take)
    
    prev_subset <- bind_rows(
      prev_subset,
      prev_extra %>% mutate(selection_type = "extra")
    )
  }
  
  if (n_curr_extra > 0 && nrow(curr_remaining) > 0) {
    ref_vals <- prev_subset$log_og_diameter
    if (length(ref_vals) == 0) ref_vals <- prev_df$log_og_diameter
    
    n_curr_take <- min(n_curr_extra, nrow(curr_remaining))
    
    curr_extra <- curr_remaining %>%
      mutate(
        dist_to_ref = map_dbl(log_og_diameter, ~ min(abs(.x - ref_vals))),
        core_rank = if (prioritize_core) core_priority else 99L,
        rand_key = runif(n())
      ) %>%
      arrange(dist_to_ref, core_rank, rand_key) %>%
      slice_head(n = n_curr_take)
    
    curr_subset <- bind_rows(
      curr_subset,
      curr_extra %>% mutate(selection_type = "extra")
    )
  }
  
  list(
    previous = prev_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank", "rand_key"))),
    current = curr_subset %>%
      select(-starts_with(".row_"), -any_of(c("dist_to_ref", "core_rank", "rand_key")))
  )
}

generate_random_matched_subsets <- function(
    candidates_previous,
    candidates_current,
    target_plan,
    strata_vars,
    n_reps = 100L,
    seed = NULL,
    prioritize_core = FALSE
) {
  if (!is.null(seed)) set.seed(seed)
  
  target_strata <- target_plan %>%
    filter(n_previous_target > 0 | n_current_target > 0) %>%
    select(all_of(strata_vars), n_previous_target, n_current_target)
  
  out_prev <- vector("list", n_reps)
  out_curr <- vector("list", n_reps)
  
  for (i in seq_len(n_reps)) {
    prev_list <- vector("list", nrow(target_strata))
    curr_list <- vector("list", nrow(target_strata))
    
    for (j in seq_len(nrow(target_strata))) {
      key_row <- target_strata[j, , drop = FALSE]
      
      prev_pool <- candidates_previous
      curr_pool <- candidates_current
      
      for (v in strata_vars) {
        prev_pool <- prev_pool %>% filter(.data[[v]] == key_row[[v]][1])
        curr_pool <- curr_pool %>% filter(.data[[v]] == key_row[[v]][1])
      }
      
      sel_j <- select_random_matched_subset(
        prev_df = prev_pool,
        curr_df = curr_pool,
        n_prev_target = key_row$n_previous_target[1],
        n_curr_target = key_row$n_current_target[1],
        prioritize_core = prioritize_core
      )
      
      prev_list[[j]] <- sel_j$previous
      curr_list[[j]] <- sel_j$current
    }
    
    out_prev[[i]] <- bind_rows(prev_list) %>%
      mutate(replicate_id = i)
    
    out_curr[[i]] <- bind_rows(curr_list) %>%
      mutate(replicate_id = i)
  }
  
  list(
    previous = bind_rows(out_prev),
    current = bind_rows(out_curr)
  )
}

#-----------------------------------------------------------
# 3. Main per-country function
#-----------------------------------------------------------
run_country_subset <- function(country_name,
                               tree_meta_clean,
                               nonempty_tree_period,
                               allow_difference = 2L,
                               n_random_reps = 100L) {
  
  country_name <- as.character(country_name)
  
  cat("\n====================================================\n")
  cat("PROCESSING COUNTRY:", country_name, "\n")
  cat("====================================================\n")
  
  meta_country_raw <- tree_meta_clean %>%
    mutate(
      country = as.character(country),
      forest = standardize_forest(forest),
      tree_id = as.character(tree_id),
      period = as.character(period),
      son_location = if ("son_location" %in% names(.)) {
        as.character(son_location)
      } else {
        NA_character_
      }
    ) %>%
    filter(country == country_name) %>%
    select(any_of(c(
      "country", "forest", "tree_id", "period",
      "log_av_ds", "log_og_diameter", "son_location"
    )))
  
  meta_country <- meta_country_raw %>%
    inner_join(
      nonempty_tree_period %>%
        mutate(
          country = as.character(country),
          forest = standardize_forest(forest),
          tree_id = as.character(tree_id),
          period = as.character(period)
        ) %>%
        filter(country == country_name),
      by = c("country", "forest", "tree_id", "period")
    )
  
  cat("\nRows after join with non-empty community filter:\n")
  print(meta_country %>% count(forest, period), n = Inf)
  cat("\n")
  
  if (nrow(meta_country) == 0) {
    cat("No rows available after filtering. Skipping.\n")
    return(NULL)
  }
  
  if ("son_location" %in% names(meta_country) && country_name == "Belgium") {
    meta_country <- meta_country %>%
      mutate(
        core_priority = case_when(
          forest == "Sonian" & son_location == "core" ~ 1L,
          forest == "Sonian" & son_location == "edge" ~ 2L,
          TRUE ~ 3L
        )
      )
  } else {
    meta_country <- meta_country %>%
      mutate(core_priority = 1L)
  }
  
  meta_country <- meta_country %>%
    mutate(
      log_av_ds = as.integer(round(log_av_ds)),
      decay_group = case_when(
        log_av_ds == 1 ~ "1",
        log_av_ds == 2 ~ "2",
        log_av_ds == 3 ~ "3",
        log_av_ds == 4 ~ "4",
        log_av_ds %in% 5:6 ~ "5_6",
        TRUE ~ NA_character_
      ),
      decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6"))
    )
  
  if (country_name == "Denmark") {
    strata_vars <- c("forest", "decay_group")
    cat("Matching: forest x decay_group + diameter optimization\n\n")
  } else {
    strata_vars <- c("decay_group")
    cat("Matching: decay_group + diameter optimization\n\n")
  }
  
  current_dist <- meta_country %>%
    filter(period == "current") %>%
    count(across(all_of(strata_vars)), name = "n_current")
  
  previous_dist <- meta_country %>%
    filter(period == "previous") %>%
    count(across(all_of(strata_vars)), name = "n_previous")
  
  target_plan <- full_join(current_dist, previous_dist, by = strata_vars) %>%
    mutate(
      n_current = replace_na(n_current, 0L),
      n_previous = replace_na(n_previous, 0L)
    ) %>%
    bind_cols(
      purrr::pmap_dfr(
        list(.$n_previous, .$n_current),
        ~ compute_flexible_targets(..1, ..2, allow_difference)
      )
    ) %>%
    arrange(across(all_of(strata_vars)))
  
  cat("\nTarget plan:\n")
  print(target_plan, n = Inf)
  
  candidates <- meta_country %>%
    filter(if_all(all_of(strata_vars), ~ !is.na(.x))) %>%
    arrange(period, across(all_of(strata_vars)), core_priority, tree_id) %>%
    group_by(period, across(all_of(strata_vars))) %>%
    mutate(candidate_rank = row_number()) %>%
    ungroup()
  
  candidates_previous <- candidates %>% filter(period == "previous")
  candidates_current  <- candidates %>% filter(period == "current")
  
  target_rows <- split(target_plan, seq_len(nrow(target_plan)))
  
  subset_prev_list <- vector("list", length(target_rows))
  subset_curr_list <- vector("list", length(target_rows))
  pair_diag_list   <- vector("list", length(target_rows))
  
  for (i in seq_along(target_rows)) {
    row_i <- target_rows[[i]]
    
    prev_pool <- candidates_previous
    curr_pool <- candidates_current
    
    for (v in strata_vars) {
      prev_pool <- prev_pool %>% filter(.data[[v]] == row_i[[v]][1])
      curr_pool <- curr_pool %>% filter(.data[[v]] == row_i[[v]][1])
    }
    
    sel_i <- select_matched_subset(
      prev_df = prev_pool,
      curr_df = curr_pool,
      n_prev_target = row_i$n_previous_target[1],
      n_curr_target = row_i$n_current_target[1],
      prioritize_core = (country_name == "Belgium")
    )
    
    subset_prev_list[[i]] <- sel_i$previous
    subset_curr_list[[i]] <- sel_i$current
    
    if (nrow(sel_i$pairs) > 0) {
      pair_diag_list[[i]] <- bind_cols(
        row_i[rep(1, nrow(sel_i$pairs)), strata_vars, drop = FALSE],
        sel_i$pairs %>%
          mutate(
            prev_tree_id = prev_pool$tree_id[prev_row],
            curr_tree_id = curr_pool$tree_id[curr_row],
            prev_log_og_diameter = prev_pool$log_og_diameter[prev_row],
            curr_log_og_diameter = curr_pool$log_og_diameter[curr_row]
          )
      )
    } else {
      pair_diag_list[[i]] <- tibble()
    }
  }
  
  prev_subset <- bind_rows(subset_prev_list)
  curr_subset <- bind_rows(subset_curr_list)
  pair_diagnostic <- bind_rows(pair_diag_list)
  
  comparison_table <- target_plan %>%
    left_join(
      prev_subset %>% count(across(all_of(strata_vars)), name = "n_previous_subset"),
      by = strata_vars
    ) %>%
    left_join(
      curr_subset %>% count(across(all_of(strata_vars)), name = "n_current_subset"),
      by = strata_vars
    ) %>%
    mutate(
      n_previous_subset = replace_na(n_previous_subset, 0L),
      n_current_subset = replace_na(n_current_subset, 0L)
    ) %>%
    arrange(across(all_of(strata_vars)))
  
  random_subsets <- generate_random_matched_subsets(
    candidates_previous = candidates_previous,
    candidates_current = candidates_current,
    target_plan = target_plan,
    strata_vars = strata_vars,
    n_reps = n_random_reps,
    seed = 42,
    prioritize_core = (country_name == "Belgium")
  )
  
  all_decay_original <- meta_country %>%
    count(period, decay_group) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = 0)
  
  all_decay_subset <- bind_rows(prev_subset, curr_subset) %>%
    count(period, decay_group) %>%
    pivot_wider(names_from = period, values_from = n, values_fill = 0)
  
  summary_decay <- full_join(
    all_decay_original,
    all_decay_subset,
    by = "decay_group",
    suffix = c("", "_subset")
  ) %>%
    mutate(
      previous = replace_na(previous, 0L),
      current = replace_na(current, 0L),
      previous_subset = replace_na(previous_subset, 0L),
      current_subset = replace_na(current_subset, 0L),
      pct_previous_retained = safe_pct(previous_subset, previous),
      pct_current_retained = safe_pct(current_subset, current),
      country = country_name
    ) %>%
    select(
      country, decay_group,
      previous, current,
      previous_subset, current_subset,
      pct_previous_retained, pct_current_retained
    ) %>%
    arrange(decay_group)
  
  cat("\nCountry-level decay summary:\n")
  print(summary_decay, n = Inf)
  
  if (country_name == "Denmark") {
    forest_decay_summary <- comparison_table %>%
      group_by(forest, decay_group) %>%
      summarise(
        n_previous = sum(n_previous, na.rm = TRUE),
        n_current = sum(n_current, na.rm = TRUE),
        n_previous_target = sum(n_previous_target, na.rm = TRUE),
        n_current_target = sum(n_current_target, na.rm = TRUE),
        n_previous_subset = sum(n_previous_subset, na.rm = TRUE),
        n_current_subset = sum(n_current_subset, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(forest, decay_group)
    
    cat("\n====================\n")
    cat("DENMARK FOREST x DECAY DIAGNOSTIC\n")
    cat("====================\n")
    print(forest_decay_summary, n = Inf)
    cat("\n")
  } else {
    forest_decay_summary <- NULL
  }
  
  diameter_summary <- bind_rows(
    prev_subset %>% mutate(period = "previous_subset"),
    curr_subset %>% mutate(period = "current_subset")
  ) %>%
    group_by(period, decay_group) %>%
    summarise(
      n = n(),
      mean_log_og_diameter = mean(log_og_diameter, na.rm = TRUE),
      median_log_og_diameter = median(log_og_diameter, na.rm = TRUE),
      sd_log_og_diameter = sd(log_og_diameter, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(country = country_name) %>%
    select(
      country, period, decay_group, n,
      mean_log_og_diameter, median_log_og_diameter, sd_log_og_diameter
    )
  
  cat("\nDiameter summary:\n")
  print(diameter_summary, n = Inf)
  
  included_keys <- bind_rows(
    prev_subset %>% select(country, forest, tree_id, period),
    curr_subset %>% select(country, forest, tree_id, period)
  ) %>%
    distinct() %>%
    mutate(inclusion_status = "included")
  
  plot_df_inclusion <- meta_country %>%
    left_join(
      included_keys,
      by = c("country", "forest", "tree_id", "period")
    ) %>%
    mutate(
      inclusion_status = if_else(is.na(inclusion_status), "excluded", inclusion_status)
    )
  
  p_cluster <- make_cluster_plot(meta_country, country_name)
  print(p_cluster)
  
  p_cluster_count <- make_cluster_count_plot(meta_country, country_name)
  print(p_cluster_count)
  
  p_inclusion <- make_inclusion_plot(plot_df_inclusion, country_name)
  print(p_inclusion)
  
  list(
    country = country_name,
    strata_vars = strata_vars,
    target_plan = target_plan,
    comparison_table = comparison_table,
    previous_subset = prev_subset,
    current_subset = curr_subset,
    random_previous_subset = random_subsets$previous,
    random_current_subset = random_subsets$current,
    summary_decay = summary_decay,
    diameter_summary = diameter_summary,
    forest_decay_summary = forest_decay_summary,
    pair_diagnostic = pair_diagnostic,
    plot_inclusion_data = plot_df_inclusion,
    plot_cluster = p_cluster,
    plot_cluster_count = p_cluster_count,
    plot_inclusion = p_inclusion
  )
}

#-----------------------------------------------------------
# 4. Run all countries
#-----------------------------------------------------------
tree_meta_clean_std <- tree_meta_clean %>%
  mutate(
    country = as.character(country),
    forest = standardize_forest(forest),
    tree_id = as.character(tree_id),
    period = as.character(period)
  )

country_vec <- tree_meta_clean_std %>%
  filter(!is.na(country)) %>%
  distinct(country) %>%
  pull(country) %>%
  sort()

cat("\nCountries to process:\n")
print(country_vec)

results <- map(
  country_vec,
  ~ run_country_subset(
    country_name = .x,
    tree_meta_clean = tree_meta_clean_std,
    nonempty_tree_period = nonempty_tree_period,
    allow_difference = allow_difference,
    n_random_reps = n_random_reps
  )
)

results <- results[!vapply(results, is.null, logical(1))]

#-----------------------------------------------------------
# 5. Combined country-level decay summary
#-----------------------------------------------------------
all_decay <- bind_rows(map(results, "summary_decay")) %>%
  arrange(country, decay_group)

cat("\n====================\n")
cat("COMBINED COUNTRY-LEVEL DECAY SUMMARY\n")
cat("====================\n")
print(all_decay, n = Inf)
cat("\n")

#-----------------------------------------------------------
# 6. Combined retained diameter summary
#-----------------------------------------------------------
all_diameter_summary <- bind_rows(map(results, "diameter_summary")) %>%
  arrange(country, period, decay_group)

cat("\n====================\n")
cat("COMBINED RETAINED DIAMETER SUMMARY\n")
cat("====================\n")
print(all_diameter_summary, n = Inf)
cat("\n")

#-----------------------------------------------------------
# 7. Optional combined Denmark forest x decay diagnostic
#-----------------------------------------------------------
all_denmark_forest_decay <- bind_rows(
  map(results, "forest_decay_summary")
)

if (nrow(all_denmark_forest_decay) > 0) {
  cat("\n====================\n")
  cat("COMBINED DENMARK FOREST x DECAY DIAGNOSTIC\n")
  cat("====================\n")
  print(all_denmark_forest_decay, n = Inf)
  cat("\n")
}

#-----------------------------------------------------------
# 8. Prepare difference table
#-----------------------------------------------------------
diam_diff <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = factor(period, levels = c("previous_subset", "current_subset"))
  ) %>%
  select(country, period, decay_group, n, mean_log_og_diameter) %>%
  pivot_wider(
    names_from = period,
    values_from = c(n, mean_log_og_diameter),
    names_sep = "__"
  ) %>%
  mutate(
    diff_mean_log_og_diameter =
      mean_log_og_diameter__current_subset - mean_log_og_diameter__previous_subset
  ) %>%
  arrange(country, decay_group)

cat("\nDifference table:\n")
print(diam_diff, n = Inf)

#-----------------------------------------------------------
# 9. Figure: current - previous by decay class and country
#-----------------------------------------------------------
p_diff <- ggplot(
  diam_diff,
  aes(x = decay_group, y = diff_mean_log_og_diameter, group = country)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = 2) +
  geom_segment(
    aes(xend = decay_group, y = 0, yend = diff_mean_log_og_diameter),
    linewidth = 0.6,
    alpha = 0.7
  ) +
  geom_point(aes(shape = country), size = 3) +
  facet_wrap(~ country, nrow = 1) +
  labs(
    x = "Decay group",
    y = expression(Delta * " mean log original diameter (current - previous)"),
    title = "Difference in average diameter between current and previous subsets",
    subtitle = "Shown separately for each country across decay classes"
  ) +
  theme_bw()

print(p_diff)

#-----------------------------------------------------------
# 10. Figure: retained mean diameter
#-----------------------------------------------------------
p_means <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = case_when(
      period == "previous_subset" ~ "previous",
      period == "current_subset" ~ "current",
      TRUE ~ as.character(period)
    ),
    period = factor(period, levels = c("previous", "current"))
  ) %>%
  ggplot(
    aes(
      x = period,
      y = mean_log_og_diameter,
      group = decay_group,
      shape = period
    )
  ) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_shape_manual(
    values = period_shapes,
    name = "Survey period",
    labels = c(
      previous = "Previous",
      current = "Current"
    )
  ) +
  facet_grid(country ~ decay_group) +
  labs(
    x = NULL,
    y = "Mean log original diameter",
    title = "Average diameter in previous and current subsets",
    subtitle = "Faceted by country and decay group"
  ) +
  theme_bw()

print(p_means)

#-----------------------------------------------------------
# 11. Figure: retained mean diameter ± SD
#-----------------------------------------------------------
p_means_sd <- all_diameter_summary %>%
  mutate(
    decay_group = factor(decay_group, levels = c("1", "2", "3", "4", "5_6")),
    period = case_when(
      period == "previous_subset" ~ "previous",
      period == "current_subset" ~ "current",
      TRUE ~ as.character(period)
    ),
    period = factor(period, levels = c("previous", "current"))
  ) %>%
  ggplot(
    aes(
      x = period,
      y = mean_log_og_diameter,
      group = decay_group,
      shape = period
    )
  ) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_errorbar(
    aes(
      ymin = mean_log_og_diameter - sd_log_og_diameter,
      ymax = mean_log_og_diameter + sd_log_og_diameter
    ),
    width = 0.12,
    linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  scale_shape_manual(
    values = period_shapes,
    name = "Survey period",
    labels = c(
      previous = "Previous",
      current = "Current"
    )
  ) +
  facet_grid(country ~ decay_group) +
  labs(
    x = NULL,
    y = "Mean log original diameter ± SD",
    title = "Average diameter in previous and current subsets",
    subtitle = "Panels show country × decay group; bars represent ±1 SD"
  ) +
  theme_bw()

print(p_means_sd)