# =========================================================
# FULL DATASET ONLY
# Global c* iNEXT workflow
#
# Outputs:
# 1. iNEXT size-based curves with interpolation + extrapolation
# 2. extrapolation endpoint = 2 x observed n trees per forest x period
# 3. one global c* across all forest x period combinations
# 4. diversity estimates at global c*
# 5. approximate period-difference tests at global c*
# 6. A1-style tables
# 7. updated p_inext_full_only plot:
#    - previous/current = lighter/darker forest colours
#    - extrapolated curve segments = dotted
#    - plot visually clipped at 200 trees
# =========================================================

setwd("/data/gent/vo/001/gvo00142/vsc45818/FRUITBODY_PAPER")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(iNEXT)
  library(colorspace)
})

# ---------------------------------------------------------
# 0. SETTINGS
# ---------------------------------------------------------

out_dir <- "tlsr_clsr_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

forest_levels <- c("Sonian", "Rajhenav", "Strødam", "Suserup")
country_levels <- c("Belgium", "Denmark", "Slovenia")
timespan_levels <- c("previous", "current")

timespan_labels <- c(
  previous = "Historical",
  current  = "Contemporary"
)

q_vals <- c(0, 1, 2)

q_labels <- c(
  "0" = "q = 0: species richness",
  "1" = "q = 1: typical diversity",
  "2" = "q = 2: dominant diversity"
)

forest_cols <- c(
  "Sonian"   = "#B8860B",
  "Rajhenav" = "#2E7D32",
  "Strødam"  = "#466C95",
  "Suserup"  = "#5FA8D3"
)

timespan_shapes <- c(
  previous = 21,
  current  = 24
)

nboot_inext <- 50
seed_inext <- 42
base_font_size <- 10

plot_x_max <- 200

# ---------------------------------------------------------
# 1. PREPARE FULL TREE-LEVEL PRESENCE-ABSENCE MATRIX
# ---------------------------------------------------------

stopifnot(exists("X_use_full"))
stopifnot(exists("meta_use_full"))

X_full <- as.data.frame(X_use_full)
meta_full <- as.data.frame(meta_use_full)

stopifnot(nrow(X_full) == nrow(meta_full))
stopifnot(all(c("tree_id", "country", "forest", "timespan") %in% names(meta_full)))

if (!"sample_id" %in% names(meta_full)) {
  meta_full <- meta_full %>%
    mutate(sample_id = rownames(X_full))
}

if (is.null(rownames(X_full)) || any(rownames(X_full) == "")) {
  rownames(X_full) <- meta_full$sample_id
}

if (!identical(rownames(X_full), as.character(meta_full$sample_id))) {
  rownames(X_full) <- as.character(meta_full$sample_id)
}

taxa_cols_full <- names(X_full)

X_full <- X_full %>%
  mutate(across(everything(), ~ as.integer(.x > 0)))

meta_full <- meta_full %>%
  mutate(
    sample_id = as.character(sample_id),
    tree_id   = as.character(tree_id),
    country   = factor(as.character(country), levels = country_levels),
    forest    = factor(as.character(forest), levels = forest_levels),
    timespan  = factor(as.character(timespan), levels = timespan_levels)
  )

dat_full <- bind_cols(
  meta_full %>%
    select(sample_id, tree_id, country, forest, timespan),
  X_full
)

full_tree_pa <- dat_full %>%
  group_by(country, forest, timespan, tree_id) %>%
  summarise(
    across(all_of(taxa_cols_full), ~ as.integer(any(.x > 0, na.rm = TRUE))),
    n_sample_rows = n(),
    .groups = "drop"
  ) %>%
  filter(
    !is.na(forest),
    !is.na(timespan)
  )

cat("\n==============================\n")
cat("Full tree-level PA matrix\n")
cat("==============================\n")
cat("Trees:", nrow(full_tree_pa), "\n")
cat("Taxa:", length(taxa_cols_full), "\n\n")
print(table(full_tree_pa$forest, full_tree_pa$timespan, useNA = "ifany"))

# ---------------------------------------------------------
# 2. MAKE INCIDENCE-FREQUENCY OBJECTS
#
# Important:
# c(n_trees, species incidence frequencies) = incidence_freq
# ---------------------------------------------------------

incidence_full_by_group <- full_tree_pa %>%
  group_by(forest, timespan) %>%
  group_split() %>%
  lapply(function(df) {
    
    this_forest <- as.character(unique(df$forest))
    this_timespan <- as.character(unique(df$timespan))
    group_name <- paste(this_forest, this_timespan, sep = "__")
    
    mat <- df %>%
      select(all_of(taxa_cols_full)) %>%
      as.data.frame()
    
    mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
    
    incidence_freq <- c(
      nrow(mat),
      colSums(mat, na.rm = TRUE)
    )
    
    tibble(
      group_name = group_name,
      forest = this_forest,
      timespan = this_timespan,
      n_trees_direct = nrow(mat),
      n_taxa_observed_direct = ncol(mat),
      incidence_freq = list(incidence_freq)
    )
  }) %>%
  bind_rows() %>%
  mutate(
    forest = factor(forest, levels = forest_levels),
    timespan = factor(timespan, levels = timespan_levels)
  ) %>%
  arrange(forest, timespan)

incidence_list_full <- incidence_full_by_group$incidence_freq
names(incidence_list_full) <- incidence_full_by_group$group_name

cat("\n==============================\n")
cat("Sanity check: incidence-frequency input objects\n")
cat("==============================\n")
print(
  incidence_full_by_group %>%
    select(group_name, forest, timespan, n_trees_direct, n_taxa_observed_direct),
  n = Inf
)

print(
  tibble(
    group_name = names(incidence_list_full),
    first_value_n_trees = sapply(incidence_list_full, function(x) x[1]),
    n_species_entries = sapply(incidence_list_full, function(x) length(x) - 1)
  ),
  n = Inf
)

# ---------------------------------------------------------
# 3. OBSERVED COVERAGE PER FOREST x TIMESPAN
# ---------------------------------------------------------

coverage_observed_full <- lapply(seq_along(incidence_list_full), function(i) {
  
  this_group <- names(incidence_list_full)[i]
  this_x <- incidence_list_full[[i]]
  
  info <- iNEXT::DataInfo(
    x = setNames(list(this_x), this_group),
    datatype = "incidence_freq"
  ) %>%
    as_tibble()
  
  tibble(
    group_name = this_group,
    n_trees = as.numeric(info$T[1]),
    observed_richness = as.numeric(info$S.obs[1]),
    observed_coverage = as.numeric(info$SC[1])
  )
}) %>%
  bind_rows() %>%
  left_join(
    incidence_full_by_group %>%
      select(group_name, forest, timespan, n_trees_direct, n_taxa_observed_direct),
    by = "group_name"
  ) %>%
  select(
    forest,
    timespan,
    group_name,
    n_trees,
    n_trees_direct,
    observed_richness,
    n_taxa_observed_direct,
    observed_coverage
  ) %>%
  arrange(forest, timespan)

cat("\n==============================\n")
cat("Observed sample coverage per forest x timespan\n")
cat("==============================\n")
print(coverage_observed_full, n = Inf)

# ---------------------------------------------------------
# 4. DEFINE ONE GLOBAL c*
#
# global c* = lowest observed sample coverage across all
# forest x timespan combinations.
# This ensures all formal diversity deltas are estimated at
# the same coverage and without extrapolation.
# ---------------------------------------------------------

global_c_star <- min(coverage_observed_full$observed_coverage, na.rm = TRUE)
global_c_star <- global_c_star-0.0001 # for safety
global_cstar_table_full <- coverage_observed_full %>%
  mutate(global_c_star = global_c_star) %>%
  select(
    forest,
    timespan,
    n_trees,
    observed_richness,
    observed_coverage,
    global_c_star
  )

cat("\n==============================\n")
cat("Global c* across all forest x timespan combinations\n")
cat("==============================\n")
print(global_c_star)
print(global_cstar_table_full, n = Inf)

write_csv2(
  coverage_observed_full,
  file.path(out_dir, "full_observed_coverage_by_forest_timespan.csv")
)

write_csv2(
  global_cstar_table_full,
  file.path(out_dir, "full_global_cstar_table.csv")
)

# ---------------------------------------------------------
# 5. RECALCULATE AND EXTRACT iNEXT SIZE-BASED CURVES
#
# Extrapolation endpoint = 2 x observed number of trees
#
# Example:
# 50 observed trees -> extrapolate to 100 trees
# 120 observed trees -> extrapolate to 240 trees
#
# This affects only the plotted iNEXT curves.
# The formal tests still use global c* and do not rely on
# extrapolation.
# ---------------------------------------------------------

set.seed(seed_inext)

inext_size_based_full_only <- lapply(seq_len(nrow(incidence_full_by_group)), function(i) {
  
  this_row <- incidence_full_by_group[i, ]
  
  this_forest <- as.character(this_row$forest)
  this_timespan <- as.character(this_row$timespan)
  this_group_name <- as.character(this_row$group_name)
  this_n_trees <- as.numeric(this_row$n_trees_direct)
  this_endpoint <- 2 * this_n_trees
  
  this_list <- list(this_row$incidence_freq[[1]])
  names(this_list) <- this_timespan
  
  out <- iNEXT::iNEXT(
    x = this_list,
    q = q_vals,
    datatype = "incidence_freq",
    endpoint = this_endpoint,
    knots = 80,
    nboot = nboot_inext,
    se = TRUE,
    conf = 0.95
  )
  
  tmp <- out$iNextEst$size_based %>%
    as_tibble()
  
  size_col <- intersect(c("m", "t"), names(tmp))[1]
  
  if (is.na(size_col)) {
    stop(
      "Could not identify the iNEXT sample-size column. Available columns are: ",
      paste(names(tmp), collapse = ", ")
    )
  }
  
  tmp %>%
    mutate(
      forest = this_forest,
      timespan = this_timespan,
      group_name = this_group_name,
      n_trees_observed = this_n_trees,
      extrapolation_endpoint = this_endpoint,
      n_trees_curve = as.numeric(.data[[size_col]])
    )
}) %>%
  bind_rows() %>%
  mutate(
    forest = factor(forest, levels = forest_levels),
    timespan = factor(timespan, levels = timespan_levels),
    q_label = factor(
      q_labels[as.character(Order.q)],
      levels = unname(q_labels[as.character(q_vals)])
    ),
    curve_segment = case_when(
      Method == "Extrapolation" ~ "Extrapolation",
      TRUE ~ "Observed/interpolated"
    )
  )

if (any(is.na(inext_size_based_full_only$n_trees_curve))) {
  stop("Some n_trees_curve values are NA. Check the iNEXT size_based output.")
}

cat("\n==============================\n")
cat("iNEXT size-based curve table\n")
cat("Extrapolation endpoint = 2 x observed number of trees\n")
cat("==============================\n")
print(
  inext_size_based_full_only %>%
    select(
      forest,
      timespan,
      Order.q,
      q_label,
      n_trees_observed,
      extrapolation_endpoint,
      n_trees_curve,
      qD,
      qD.LCL,
      qD.UCL,
      SC,
      Method
    ) %>%
    arrange(forest, timespan, Order.q, n_trees_curve),
  n = 40
)

endpoint_check_full <- inext_size_based_full_only %>%
  distinct(
    forest,
    timespan,
    n_trees_observed,
    extrapolation_endpoint
  ) %>%
  arrange(forest, timespan)

cat("\n==============================\n")
cat("Extrapolation endpoint check\n")
cat("==============================\n")
print(endpoint_check_full, n = Inf)

write_csv2(
  inext_size_based_full_only,
  file.path(out_dir, "full_inext_size_based_curves_with_2x_extrapolation.csv")
)

write_csv2(
  endpoint_check_full,
  file.path(out_dir, "full_inext_2x_extrapolation_endpoint_check.csv")
)

# ---------------------------------------------------------
# 6. OBSERVED POINTS FOR THE CURVE PLOT
# ---------------------------------------------------------

observed_points_full_only <- inext_size_based_full_only %>%
  filter(Method == "Observed") %>%
  group_by(forest, timespan, Order.q, q_label) %>%
  slice_max(order_by = n_trees_curve, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---------------------------------------------------------
# 7. UPDATED p_inext_full_only PLOT
#
# - timespan shown as lighter/darker version of forest colour
# - extrapolated curve segments dotted
# - ribbon does not inherit linetype
# - observed points retained
# - x-axis visually clipped at 200 trees
# ---------------------------------------------------------

forest_cols_timespan_inext <- c(
  setNames(
    colorspace::lighten(forest_cols, amount = 0.28),
    paste0(names(forest_cols), "__previous")
  ),
  setNames(
    colorspace::darken(forest_cols, amount = 0.12),
    paste0(names(forest_cols), "__current")
  )
)

forest_timespan_breaks <- as.vector(
  outer(forest_levels, timespan_levels, paste, sep = "__")
)

forest_timespan_labels <- setNames(
  paste(
    rep(forest_levels, each = length(timespan_levels)),
    rep(timespan_labels[timespan_levels], times = length(forest_levels)),
    sep = " - "
  ),
  forest_timespan_breaks
)

inext_size_based_full_only_plot <- inext_size_based_full_only %>%
  mutate(
    forest_timespan = paste0(as.character(forest), "__", as.character(timespan)),
    forest_timespan = factor(forest_timespan, levels = forest_timespan_breaks),
    curve_segment = factor(
      curve_segment,
      levels = c("Observed/interpolated", "Extrapolation")
    )
  )

observed_points_full_only_plot <- observed_points_full_only %>%
  mutate(
    forest_timespan = paste0(as.character(forest), "__", as.character(timespan)),
    forest_timespan = factor(forest_timespan, levels = forest_timespan_breaks)
  )

p_inext_full_only <- ggplot(
  inext_size_based_full_only_plot,
  aes(
    x = n_trees_curve,
    y = qD,
    colour = forest_timespan,
    linetype = curve_segment,
    group = interaction(forest, timespan, curve_segment)
  )
) +
  geom_ribbon(
    data = inext_size_based_full_only_plot,
    aes(
      x = n_trees_curve,
      ymin = qD.LCL,
      ymax = qD.UCL,
      fill = forest_timespan,
      group = interaction(forest, timespan)
    ),
    inherit.aes = FALSE,
    alpha = 0.10,
    colour = NA,
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  geom_line(
    linewidth = 0.50,
    na.rm = TRUE
  ) +
  geom_point(
    data = observed_points_full_only_plot,
    aes(
      x = n_trees_curve,
      y = qD,
      shape = timespan,
      colour = forest_timespan,
      fill = forest_timespan
    ),
    size = 1.9,
    stroke = 0.30,
    na.rm = TRUE,
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ q_label,
    nrow = 1,
    ncol = 3,
    scales = "fixed"
  ) +
  scale_x_continuous(
    breaks = seq(0, plot_x_max, by = 50),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  coord_cartesian(
    xlim = c(0, plot_x_max)
  ) +
  scale_colour_manual(
    values = forest_cols_timespan_inext,
    breaks = forest_timespan_breaks,
    labels = forest_timespan_labels,
    drop = FALSE
  ) +
  scale_fill_manual(
    values = forest_cols_timespan_inext,
    breaks = forest_timespan_breaks,
    labels = forest_timespan_labels,
    drop = FALSE
  ) +
  scale_linetype_manual(
    values = c(
      "Observed/interpolated" = "solid",
      "Extrapolation" = "dotted"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(
    values = timespan_shapes,
    breaks = timespan_levels,
    labels = timespan_labels,
    drop = FALSE
  ) +
  labs(
    x = "Number of deadwood trees",
    y = "Effective species diversity",
    colour = "Forest × time period",
    fill = "Forest × time period",
    linetype = "Curve segment",
    shape = "Time period"
  ) +
  guides(
    fill = "none",
    colour = guide_legend(
      override.aes = list(
        linetype = "solid",
        linewidth = 0.8,
        shape = NA
      )
    ),
    linetype = guide_legend(
      override.aes = list(
        colour = "black",
        linewidth = 0.8
      )
    )
  ) +
  theme_bw(base_size = base_font_size) +
  theme(
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

print(p_inext_full_only)

ggsave(
  file.path(out_dir, "full_inext_size_based_global_cstar_with_2x_extrapolation.png"),
  plot = p_inext_full_only,
  width = 10.5,
  height = 5.2,
  dpi = 900
)

# ---------------------------------------------------------
# 8. ESTIMATE DIVERSITY AT GLOBAL c*
# ---------------------------------------------------------

diversity_at_global_cstar_full <- lapply(forest_levels, function(this_forest) {
  
  these_groups <- incidence_full_by_group %>%
    filter(as.character(forest) == this_forest) %>%
    arrange(timespan)
  
  this_list <- these_groups$incidence_freq
  names(this_list) <- as.character(these_groups$timespan)
  
  est <- iNEXT::estimateD(
    x = this_list,
    q = q_vals,
    datatype = "incidence_freq",
    base = "coverage",
    level = global_c_star,
    conf = 0.95
  )
  
  est %>%
    as_tibble() %>%
    mutate(
      forest = this_forest,
      global_c_star = global_c_star
    )
}) %>%
  bind_rows() %>%
  rename(
    timespan = Assemblage,
    q = Order.q,
    diversity_at_global_cstar = qD,
    diversity_lcl = qD.LCL,
    diversity_ucl = qD.UCL,
    target_coverage = SC
  ) %>%
  mutate(
    forest = factor(forest, levels = forest_levels),
    timespan = factor(timespan, levels = timespan_levels),
    diversity_metric = case_when(
      q == 0 ~ "q0_species_richness",
      q == 1 ~ "q1_typical_diversity",
      q == 2 ~ "q2_dominant_diversity",
      TRUE ~ paste0("q", q)
    ),
    diversity_metric_label = case_when(
      q == 0 ~ "q = 0: species richness",
      q == 1 ~ "q = 1: typical diversity",
      q == 2 ~ "q = 2: dominant diversity",
      TRUE ~ paste0("q = ", q)
    )
  ) %>%
  select(
    forest,
    timespan,
    q,
    diversity_metric,
    diversity_metric_label,
    global_c_star,
    target_coverage,
    diversity_at_global_cstar,
    diversity_lcl,
    diversity_ucl,
    Method
  ) %>%
  arrange(forest, q, timespan)

cat("\n==============================\n")
cat("Diversity at global c* per forest x timespan\n")
cat("==============================\n")
print(diversity_at_global_cstar_full, n = Inf)

write_csv2(
  diversity_at_global_cstar_full,
  file.path(out_dir, "full_diversity_at_global_cstar_by_forest_timespan.csv")
)

# ---------------------------------------------------------
# 9. APPROXIMATE TESTS AT GLOBAL c*
#
# Delta = current - previous
# SE reconstructed from 95% CI of estimateD()
# ---------------------------------------------------------

diversity_test_at_global_cstar_full <- diversity_at_global_cstar_full %>%
  mutate(
    se_approx = (diversity_ucl - diversity_lcl) / (2 * qnorm(0.975))
  ) %>%
  select(
    forest,
    timespan,
    q,
    diversity_metric,
    diversity_metric_label,
    global_c_star,
    target_coverage,
    diversity_at_global_cstar,
    diversity_lcl,
    diversity_ucl,
    se_approx
  ) %>%
  pivot_wider(
    names_from = timespan,
    values_from = c(
      target_coverage,
      diversity_at_global_cstar,
      diversity_lcl,
      diversity_ucl,
      se_approx
    )
  ) %>%
  mutate(
    delta_current_minus_previous =
      diversity_at_global_cstar_current - diversity_at_global_cstar_previous,
    
    se_delta =
      sqrt(se_approx_current^2 + se_approx_previous^2),
    
    z_value =
      delta_current_minus_previous / se_delta,
    
    p_value_two_sided =
      2 * pnorm(abs(z_value), lower.tail = FALSE),
    
    delta_lcl_approx =
      delta_current_minus_previous - qnorm(0.975) * se_delta,
    
    delta_ucl_approx =
      delta_current_minus_previous + qnorm(0.975) * se_delta,
    
    p_label = case_when(
      p_value_two_sided < 0.001 ~ "***",
      p_value_two_sided < 0.01  ~ "**",
      p_value_two_sided < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    
    direction = case_when(
      delta_current_minus_previous > 0 ~ "increase",
      delta_current_minus_previous < 0 ~ "decrease",
      TRUE ~ "no_change"
    )
  ) %>%
  arrange(q, forest)

cat("\n==============================\n")
cat("Approximate tests at global c*\n")
cat("Delta = current - previous\n")
cat("==============================\n")
print(diversity_test_at_global_cstar_full, n = Inf)

write_csv2(
  diversity_test_at_global_cstar_full,
  file.path(out_dir, "full_test_diversity_difference_at_global_cstar_all_q_by_forest.csv")
)

# ---------------------------------------------------------
# 10. PLOT TESTED DELTAS AT GLOBAL c*
# ---------------------------------------------------------

p_diversity_test_global_cstar_all_q <- ggplot(
  diversity_test_at_global_cstar_full,
  aes(
    x = forest,
    y = delta_current_minus_previous,
    fill = forest
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dotted"
  ) +
  geom_col(
    colour = "black",
    linewidth = 0.25,
    width = 0.70
  ) +
  geom_errorbar(
    aes(
      ymin = delta_lcl_approx,
      ymax = delta_ucl_approx
    ),
    width = 0.18,
    linewidth = 0.45
  ) +
  geom_text(
    aes(
      y = ifelse(
        delta_current_minus_previous >= 0,
        delta_ucl_approx,
        delta_lcl_approx
      ),
      label = p_label
    ),
    vjust = ifelse(
      diversity_test_at_global_cstar_full$delta_current_minus_previous >= 0,
      -0.45,
      1.25
    ),
    size = 3.4
  ) +
  facet_wrap(
    ~ diversity_metric_label,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(values = forest_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = expression(Delta~"diversity at global c*"~"(current - previous)"),
    fill = "Forest"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

print(p_diversity_test_global_cstar_all_q)

ggsave(
  file.path(out_dir, "full_test_diversity_delta_at_global_cstar_all_q_by_forest.png"),
  plot = p_diversity_test_global_cstar_all_q,
  width = 10.5,
  height = 4.4,
  dpi = 400
)

# ---------------------------------------------------------
# 11. A1-STYLE TABLES AT GLOBAL c*
# ---------------------------------------------------------

format_est <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(x, format = "f", digits = digits)
  )
}

format_delta_ci <- function(delta, lcl, ucl, digits = 1) {
  ifelse(
    is.na(delta) | is.na(lcl) | is.na(ucl),
    NA_character_,
    paste0(
      formatC(delta, format = "f", digits = digits),
      " [",
      formatC(lcl, format = "f", digits = digits),
      ", ",
      formatC(ucl, format = "f", digits = digits),
      "]"
    )
  )
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

diversity_test_at_global_cstar_numeric_wide <- diversity_test_at_global_cstar_full %>%
  mutate(
    q_label = case_when(
      q == 0 ~ "q0",
      q == 1 ~ "q1",
      q == 2 ~ "q2",
      TRUE ~ paste0("q", q)
    )
  ) %>%
  select(
    forest,
    q_label,
    global_c_star,
    diversity_at_global_cstar_previous,
    diversity_at_global_cstar_current,
    delta_current_minus_previous,
    delta_lcl_approx,
    delta_ucl_approx,
    z_value,
    p_value_two_sided
  ) %>%
  pivot_wider(
    names_from = q_label,
    values_from = c(
      diversity_at_global_cstar_previous,
      diversity_at_global_cstar_current,
      delta_current_minus_previous,
      delta_lcl_approx,
      delta_ucl_approx,
      z_value,
      p_value_two_sided
    ),
    names_glue = "{q_label}_{.value}"
  ) %>%
  arrange(forest)

cat("\n==============================\n")
cat("Numeric wide diversity-test table at global c*\n")
cat("==============================\n")
print(diversity_test_at_global_cstar_numeric_wide, n = Inf)

write_csv2(
  diversity_test_at_global_cstar_numeric_wide,
  file.path(out_dir, "full_diversity_test_at_global_cstar_numeric_wide.csv")
)

diversity_test_at_global_cstar_A1_compact <- diversity_test_at_global_cstar_full %>%
  mutate(
    q_label = case_when(
      q == 0 ~ "q0",
      q == 1 ~ "q1",
      q == 2 ~ "q2",
      TRUE ~ paste0("q", q)
    ),
    
    result_fmt = paste0(
      format_delta_ci(
        delta = delta_current_minus_previous,
        lcl = delta_lcl_approx,
        ucl = delta_ucl_approx,
        digits = 1
      ),
      ", P = ",
      format_p(p_value_two_sided)
    )
  ) %>%
  select(
    forest,
    q_label,
    global_c_star,
    result_fmt
  ) %>%
  pivot_wider(
    names_from = q_label,
    values_from = result_fmt
  ) %>%
  mutate(
    global_c_star = formatC(global_c_star, format = "f", digits = 3)
  ) %>%
  rename(
    Forest = forest,
    `global c*` = global_c_star,
    `q0 richness: Δ [95% CI], P` = q0,
    `q1 typical diversity: Δ [95% CI], P` = q1,
    `q2 dominant diversity: Δ [95% CI], P` = q2
  ) %>%
  arrange(Forest)

cat("\n==============================\n")
cat("Compact A1-style table at global c*\n")
cat("==============================\n")
print(diversity_test_at_global_cstar_A1_compact, n = Inf)

write_csv2(
  diversity_test_at_global_cstar_A1_compact,
  file.path(out_dir, "full_diversity_test_at_global_cstar_A1_compact.csv")
)

# ---------------------------------------------------------
# 12. OPTIONAL MANUSCRIPT-READY METHOD NOTE
# ---------------------------------------------------------

cat("\n==============================\n")
cat("Suggested methods wording\n")
cat("==============================\n")
cat(
  "Hill diversity was estimated from tree-level incidence-frequency data using iNEXT. ",
  "For formal period comparisons, diversity was standardized to a single global sample coverage (c*) ",
  "defined as the lowest observed coverage among all forest x survey-period combinations. ",
  "Differences between contemporary and historical estimates were evaluated using approximate ",
  "two-sided z-tests, with standard errors reconstructed from the 95% confidence intervals of the ",
  "coverage-standardized estimates. ",
  "Size-based rarefaction and extrapolation curves were extrapolated to at most twice the observed ",
  "number of deadwood units within each forest x survey-period combination.\n",
  sep = ""
)




# =========================================================
# BOLT-ON BLOCK
# Number of deadwood trees corresponding to global c*
# for each forest x survey period
#
# This answers:
# At global c*, how many tree-equivalent samples are used
# in each forest x timespan combination?
# =========================================================

stopifnot(exists("incidence_full_by_group"))
stopifnot(exists("coverage_observed_full"))
stopifnot(exists("global_c_star"))
stopifnot(exists("q_vals"))
stopifnot(exists("out_dir"))

trees_at_global_cstar_full <- lapply(seq_len(nrow(incidence_full_by_group)), function(i) {
  
  this_row <- incidence_full_by_group[i, ]
  
  this_forest <- as.character(this_row$forest)
  this_timespan <- as.character(this_row$timespan)
  this_group_name <- as.character(this_row$group_name)
  
  this_list <- list(this_row$incidence_freq[[1]])
  names(this_list) <- this_group_name
  
  est <- iNEXT::estimateD(
    x = this_list,
    q = q_vals,
    datatype = "incidence_freq",
    base = "coverage",
    level = global_c_star,
    conf = 0.95
  ) %>%
    as_tibble()
  
  size_col <- intersect(
    c("m", "t", "T", "n", "sample.size", "SampleSize"),
    names(est)
  )[1]
  
  if (is.na(size_col)) {
    stop(
      "Could not identify the sample-size column in estimateD() output. Available columns are: ",
      paste(names(est), collapse = ", ")
    )
  }
  
  est %>%
    mutate(
      forest = this_forest,
      timespan = this_timespan,
      group_name = this_group_name,
      global_c_star = global_c_star,
      n_trees_at_global_cstar = as.numeric(.data[[size_col]])
    )
}) %>%
  bind_rows() %>%
  rename(
    q = Order.q,
    diversity_at_global_cstar = qD,
    diversity_lcl = qD.LCL,
    diversity_ucl = qD.UCL,
    target_coverage = SC
  ) %>%
  mutate(
    forest = factor(forest, levels = forest_levels),
    timespan = factor(timespan, levels = timespan_levels),
    diversity_metric = case_when(
      q == 0 ~ "q0_species_richness",
      q == 1 ~ "q1_typical_diversity",
      q == 2 ~ "q2_dominant_diversity",
      TRUE ~ paste0("q", q)
    )
  ) %>%
  left_join(
    coverage_observed_full %>%
      select(
        forest,
        timespan,
        n_trees_observed = n_trees,
        observed_richness,
        observed_coverage
      ),
    by = c("forest", "timespan")
  ) %>%
  select(
    forest,
    timespan,
    group_name,
    q,
    diversity_metric,
    global_c_star,
    target_coverage,
    n_trees_at_global_cstar,
    n_trees_observed,
    observed_coverage,
    observed_richness,
    diversity_at_global_cstar,
    diversity_lcl,
    diversity_ucl,
    Method
  ) %>%
  arrange(forest, timespan, q)

cat("\n==============================\n")
cat("Tree-equivalent sample size at global c*\n")
cat("==============================\n")
print(trees_at_global_cstar_full, n = Inf)

write_csv2(
  trees_at_global_cstar_full,
  file.path(out_dir, "full_tree_equivalent_sample_size_at_global_cstar_all_q.csv")
)

# ---------------------------------------------------------
# Compact version: one row per forest x period
# The tree-equivalent sample size should be identical across q
# because the target is sample coverage, not q-specific diversity.
# ---------------------------------------------------------

trees_at_global_cstar_compact <- trees_at_global_cstar_full %>%
  group_by(forest, timespan) %>%
  summarise(
    global_c_star = first(global_c_star),
    target_coverage = first(target_coverage),
    n_trees_at_global_cstar = first(n_trees_at_global_cstar),
    n_trees_at_global_cstar_rounded = round(first(n_trees_at_global_cstar), 1),
    n_trees_observed = first(n_trees_observed),
    observed_coverage = first(observed_coverage),
    observed_richness = first(observed_richness),
    Method = paste(unique(Method), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(forest, timespan)

cat("\n==============================\n")
cat("Compact tree-equivalent sample size at global c*\n")
cat("==============================\n")
print(trees_at_global_cstar_compact, n = Inf)

write_csv2(
  trees_at_global_cstar_compact,
  file.path(out_dir, "full_tree_equivalent_sample_size_at_global_cstar_compact.csv")
)

# ---------------------------------------------------------
# A1-style table
# ---------------------------------------------------------

trees_at_global_cstar_A1 <- trees_at_global_cstar_compact %>%
  mutate(
    `global c*` = formatC(global_c_star, format = "f", digits = 3),
    `target coverage` = formatC(target_coverage, format = "f", digits = 3),
    `tree-equivalent n at global c*` =
      formatC(n_trees_at_global_cstar, format = "f", digits = 1),
    `observed n` =
      formatC(n_trees_observed, format = "f", digits = 0),
    `observed coverage` =
      formatC(observed_coverage, format = "f", digits = 3),
    `observed richness` =
      formatC(observed_richness, format = "f", digits = 0)
  ) %>%
  select(
    Forest = forest,
    `Survey period` = timespan,
    `global c*`,
    `tree-equivalent n at global c*`,
    `observed n`,
    `observed coverage`,
    `observed richness`,
    Method
  )

cat("\n==============================\n")
cat("A1-style tree-equivalent table at global c*\n")
cat("==============================\n")
print(trees_at_global_cstar_A1, n = Inf)

write_csv2(
  trees_at_global_cstar_A1,
  file.path(out_dir, "full_tree_equivalent_sample_size_at_global_cstar_A1.csv")
)


# =========================================================
# SUPPLEMENTARY TABLE
# Coverage-standardized estimates, tree-equivalent sample
# sizes at C*, period differences, and test statistics
#
# Output:
# 1. Long numeric table
# 2. Compact formatted A1-style supplementary table
# =========================================================

stopifnot(exists("diversity_at_global_cstar_full"))
stopifnot(exists("diversity_test_at_global_cstar_full"))
stopifnot(exists("trees_at_global_cstar_compact"))
stopifnot(exists("global_c_star"))
stopifnot(exists("out_dir"))

# ---------------------------------------------------------
# 1. Tree-equivalent n at global C* per forest x period
# ---------------------------------------------------------

tree_n_at_cstar_for_join <- trees_at_global_cstar_compact %>%
  select(
    forest,
    timespan,
    n_trees_at_global_cstar,
    n_trees_observed,
    observed_coverage,
    observed_richness
  )

# ---------------------------------------------------------
# 2. Period-specific diversity estimates at global C*
# ---------------------------------------------------------

diversity_estimates_at_cstar_long <- diversity_at_global_cstar_full %>%
  left_join(
    tree_n_at_cstar_for_join,
    by = c("forest", "timespan")
  ) %>%
  mutate(
    se_estimate = (diversity_ucl - diversity_lcl) / (2 * qnorm(0.975)),
    q_label = case_when(
      q == 0 ~ "q0_species_richness",
      q == 1 ~ "q1_typical_diversity",
      q == 2 ~ "q2_dominant_diversity",
      TRUE ~ paste0("q", q)
    )
  ) %>%
  select(
    forest,
    timespan,
    q,
    q_label,
    global_c_star,
    target_coverage,
    n_trees_at_global_cstar,
    n_trees_observed,
    observed_coverage,
    observed_richness,
    diversity_at_global_cstar,
    diversity_lcl,
    diversity_ucl,
    se_estimate,
    Method
  ) %>%
  arrange(forest, q, timespan)

cat("\n==============================\n")
cat("Period-specific coverage-standardized estimates at global C*\n")
cat("==============================\n")
print(diversity_estimates_at_cstar_long, n = Inf)

write_csv2(
  diversity_estimates_at_cstar_long,
  file.path(out_dir, "supp_table_global_cstar_period_specific_estimates_long.csv")
)


# ---------------------------------------------------------
# 3. Full long supplementary table with tests
# One row per forest x q
# ---------------------------------------------------------

tree_n_observed_wide <- diversity_estimates_at_cstar_long %>%
  select(
    forest,
    q,
    timespan,
    n_trees_at_global_cstar,
    n_trees_observed,
    observed_coverage,
    observed_richness
  ) %>%
  pivot_wider(
    names_from = timespan,
    values_from = c(
      n_trees_at_global_cstar,
      n_trees_observed,
      observed_coverage,
      observed_richness
    )
  )

supp_global_cstar_tests_long <- diversity_test_at_global_cstar_full %>%
  mutate(
    q_label = case_when(
      q == 0 ~ "q0_species_richness",
      q == 1 ~ "q1_typical_diversity",
      q == 2 ~ "q2_dominant_diversity",
      TRUE ~ paste0("q", q)
    )
  ) %>%
  left_join(
    tree_n_observed_wide,
    by = c("forest", "q")
  ) %>%
  select(
    forest,
    q,
    q_label,
    global_c_star,
    
    n_trees_at_global_cstar_previous,
    n_trees_at_global_cstar_current,
    n_trees_observed_previous,
    n_trees_observed_current,
    observed_coverage_previous,
    observed_coverage_current,
    observed_richness_previous,
    observed_richness_current,
    
    diversity_at_global_cstar_previous,
    diversity_lcl_previous,
    diversity_ucl_previous,
    se_approx_previous,
    
    diversity_at_global_cstar_current,
    diversity_lcl_current,
    diversity_ucl_current,
    se_approx_current,
    
    delta_current_minus_previous,
    se_delta,
    delta_lcl_approx,
    delta_ucl_approx,
    z_value,
    p_value_two_sided,
    p_label,
    direction
  ) %>%
  arrange(forest, q)

cat("\n==============================\n")
cat("Supplementary global C* table with estimates and test statistics\n")
cat("==============================\n")
print(supp_global_cstar_tests_long, n = Inf)

write_csv2(
  supp_global_cstar_tests_long,
  file.path(out_dir, "supp_table_global_cstar_estimates_differences_tests_long.csv")
)

# ---------------------------------------------------------
# 4. Formatting helpers
# ---------------------------------------------------------

format_num <- function(x, digits = 1) {
  ifelse(
    is.na(x),
    NA_character_,
    formatC(x, format = "f", digits = digits)
  )
}

format_ci <- function(est, lcl, ucl, digits = 1) {
  ifelse(
    is.na(est) | is.na(lcl) | is.na(ucl),
    NA_character_,
    paste0(
      formatC(est, format = "f", digits = digits),
      " [",
      formatC(lcl, format = "f", digits = digits),
      ", ",
      formatC(ucl, format = "f", digits = digits),
      "]"
    )
  )
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ formatC(p, format = "f", digits = 3)
  )
}

# ---------------------------------------------------------
# 5. Compact A1-style supplementary table
# One row per forest x Hill number
# ---------------------------------------------------------

supp_global_cstar_A1 <- supp_global_cstar_tests_long %>%
  mutate(
    Hill_number = case_when(
      q == 0 ~ "q0 species richness",
      q == 1 ~ "q1 typical diversity",
      q == 2 ~ "q2 dominant diversity",
      TRUE ~ paste0("q", q)
    ),
    
    `C*` = formatC(global_c_star, format = "f", digits = 3),
    
    `Historical estimate [95% CI]` = format_ci(
      diversity_at_global_cstar_previous,
      diversity_lcl_previous,
      diversity_ucl_previous,
      digits = 1
    ),
    
    `Contemporary estimate [95% CI]` = format_ci(
      diversity_at_global_cstar_current,
      diversity_lcl_current,
      diversity_ucl_current,
      digits = 1
    ),
    
    `Historical tree-equivalent n at C*` =
      format_num(n_trees_at_global_cstar_previous, digits = 1),
    
    `Contemporary tree-equivalent n at C*` =
      format_num(n_trees_at_global_cstar_current, digits = 1),
    
    `Observed n historical` =
      format_num(n_trees_observed_previous, digits = 0),
    
    `Observed n contemporary` =
      format_num(n_trees_observed_current, digits = 0),
    
    `Observed coverage historical` =
      format_num(observed_coverage_previous, digits = 3),
    
    `Observed coverage contemporary` =
      format_num(observed_coverage_current, digits = 3),
    
    `Δ contemporary - historical [95% CI]` = format_ci(
      delta_current_minus_previous,
      delta_lcl_approx,
      delta_ucl_approx,
      digits = 1
    ),
    
    `SE(Δ)` = format_num(se_delta, digits = 2),
    `z` = format_num(z_value, digits = 2),
    `P` = format_p(p_value_two_sided),
    `Significance` = p_label
  ) %>%
  select(
    Forest = forest,
    `Hill number` = Hill_number,
    `C*`,
    `Historical estimate [95% CI]`,
    `Contemporary estimate [95% CI]`,
    `Historical tree-equivalent n at C*`,
    `Contemporary tree-equivalent n at C*`,
    `Observed n historical`,
    `Observed n contemporary`,
    `Observed coverage historical`,
    `Observed coverage contemporary`,
    `Δ contemporary - historical [95% CI]`,
    `SE(Δ)`,
    `z`,
    `P`,
    `Significance`
  ) %>%
  arrange(Forest, `Hill number`)

cat("\n==============================\n")
cat("A1-style supplementary table: global C* estimates and tests\n")
cat("==============================\n")
print(supp_global_cstar_A1, n = Inf)

write_csv2(
  supp_global_cstar_A1,
  file.path(out_dir, "supp_table_global_cstar_estimates_differences_tests_A1.csv")
)

# ---------------------------------------------------------
# 6. Optional wider version: one row per forest
# ---------------------------------------------------------

supp_global_cstar_A1_wide <- supp_global_cstar_A1 %>%
  mutate(
    q_short = case_when(
      `Hill number` == "q0 species richness" ~ "q0",
      `Hill number` == "q1 typical diversity" ~ "q1",
      `Hill number` == "q2 dominant diversity" ~ "q2",
      TRUE ~ `Hill number`
    )
  ) %>%
  select(
    Forest,
    q_short,
    `C*`,
    `Historical estimate [95% CI]`,
    `Contemporary estimate [95% CI]`,
    `Historical tree-equivalent n at C*`,
    `Contemporary tree-equivalent n at C*`,
    `Δ contemporary - historical [95% CI]`,
    `SE(Δ)`,
    `z`,
    `P`
  ) %>%
  pivot_wider(
    names_from = q_short,
    values_from = c(
      `Historical estimate [95% CI]`,
      `Contemporary estimate [95% CI]`,
      `Historical tree-equivalent n at C*`,
      `Contemporary tree-equivalent n at C*`,
      `Δ contemporary - historical [95% CI]`,
      `SE(Δ)`,
      `z`,
      `P`
    ),
    names_glue = "{q_short} {.value}"
  ) %>%
  arrange(Forest)

cat("\n==============================\n")
cat("Wide A1-style supplementary table: global C* estimates and tests\n")
cat("==============================\n")
print(supp_global_cstar_A1_wide, n = Inf)

write_csv2(
  supp_global_cstar_A1_wide,
  file.path(out_dir, "supp_table_global_cstar_estimates_differences_tests_A1_wide.csv")
)
