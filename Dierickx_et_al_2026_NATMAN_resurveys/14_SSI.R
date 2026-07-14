# =========================================================
# HARMONISED SSI ANALYSIS
# Full target matrix + subsetted matrix
# Plot labels use Historical / Current
# No hardcoded SSI values
# =========================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# 0. SETTINGS / KNOBS
# =========================================================

USE_TARGET_ONLY      <- TRUE
USE_LIGNICOLOUS_ONLY <- FALSE

out_dir <- "ssi_analysis_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

country_levels <- c("Belgium", "Denmark", "Slovenia")
timespan_levels <- c("previous", "current")
ssi_levels <- c("A", "B", "C")

ssi_fill_vals <- c(
  "SSI_C" = "grey75",
  "SSI_B" = "grey45",
  "SSI_A" = "black",
  "C"     = "grey75",
  "B"     = "grey45",
  "A"     = "black"
)

# =========================================================
# 1. READ AND FILTER SSI TAXONOMY
# =========================================================

ssi_df <- read_excel(
  "taxonomy_revisedJHC.xlsx",
  sheet = "eco_target"
) %>%
  mutate(
    Taxon  = str_trim(as.character(Taxon)),
#    SSI    = str_trim(as.character(SSI)),
    SSI    = str_trim(as.character(SSI_strict)),
    code   = str_trim(as.character(code)),
    TARGET = suppressWarnings(as.integer(TARGET))
  ) %>%
  filter(
    !is.na(Taxon),
    Taxon != "",
    !is.na(SSI),
    SSI != "",
    SSI != "-"
  ) %>%
  mutate(
    SSI = factor(SSI, levels = ssi_levels, ordered = TRUE)
  )

if (USE_TARGET_ONLY) {
  ssi_df <- ssi_df %>%
    filter(TARGET == 1)
}

if (USE_LIGNICOLOUS_ONLY) {
  ssi_df <- ssi_df %>%
    filter(code == "L")
}

ssi_taxa <- ssi_df %>%
  distinct(Taxon, SSI) %>%
  arrange(SSI, Taxon)

cat("\n==============================\n")
cat("SSI taxonomy after filtering\n")
cat("==============================\n")

cat("USE_TARGET_ONLY      =", USE_TARGET_ONLY, "\n")
cat("USE_LIGNICOLOUS_ONLY =", USE_LIGNICOLOUS_ONLY, "\n\n")

print(ssi_taxa %>% count(SSI, name = "n_taxa_taxonomy"), n = Inf)

write_csv2(
  ssi_taxa,
  file.path(out_dir, "ssi_taxa_used_from_taxonomy.csv")
)

# =========================================================
# 2. SMALL HELPERS
# =========================================================

make_x_id <- function(country, timespan, year_label) {
  case_when(
    country == "Belgium"  & timespan == "previous" ~ "Belgium\nHistorical",
    country == "Belgium"  & timespan == "current"  ~ "Belgium\nCurrent",
    country == "Denmark"  & timespan == "previous" ~ "Denmark\nHistorical",
    country == "Denmark"  & timespan == "current"  ~ paste0("Denmark\nCurrent ", year_label),
    country == "Slovenia" & timespan == "previous" ~ "Slovenia\nHistorical",
    country == "Slovenia" & timespan == "current"  ~ paste0("Slovenia\nCurrent ", year_label),
    TRUE ~ paste(country, year_label)
  )
}

x_id_levels <- c(
  "Belgium\nHistorical",
  "Belgium\nCurrent",
  "Denmark\nHistorical",
  "Denmark\nCurrent 2021",
  "Denmark\nCurrent 2022",
  "Slovenia\nHistorical",
  "Slovenia\nCurrent 2021",
  "Slovenia\nCurrent 2022"
)

analyse_one_ssi_set <- function(X, meta, data_scope, ssi_taxa) {
  
  stopifnot(nrow(X) == nrow(meta))
  stopifnot("sample_id" %in% names(meta))
  stopifnot("country" %in% names(meta))
  stopifnot("timespan" %in% names(meta))
  stopifnot("year" %in% names(meta))
  
  if (!identical(rownames(X), meta$sample_id)) {
    rownames(X) <- meta$sample_id
  }
  
  taxa_in_matrix <- intersect(ssi_taxa$Taxon, colnames(X))
  
  cat("\n==============================\n")
  cat(data_scope, "\n")
  cat("==============================\n")
  cat("Samples:", nrow(X), "\n")
  cat("Matrix taxa:", ncol(X), "\n")
  cat("SSI taxa present in matrix:", length(taxa_in_matrix), "\n\n")
  
  if (length(taxa_in_matrix) == 0) {
    stop("No SSI taxa from ssi_taxa were found in the columns of ", data_scope)
  }
  
  meta_small <- meta %>%
    transmute(
      sample_id  = as.character(sample_id),
      country    = as.character(country),
      forest     = as.character(forest),
      year       = as.integer(as.character(year)),
      year_label = as.character(year),
      timespan   = as.character(timespan)
    ) %>%
    mutate(
      data_scope = data_scope,
      country    = factor(country, levels = country_levels),
      timespan   = factor(timespan, levels = timespan_levels),
      year_label = factor(year_label, levels = as.character(sort(unique(year)))),
      x_id       = make_x_id(as.character(country), as.character(timespan), as.character(year_label)),
      x_id       = factor(x_id, levels = x_id_levels)
    )
  
  ssi_long <- X %>%
    select(all_of(taxa_in_matrix)) %>%
    as_tibble() %>%
    mutate(sample_id = meta_small$sample_id) %>%
    pivot_longer(
      cols = -sample_id,
      names_to = "Taxon",
      values_to = "pa"
    ) %>%
    mutate(pa = as.integer(pa > 0)) %>%
    filter(pa == 1L) %>%
    left_join(ssi_taxa, by = "Taxon") %>%
    filter(!is.na(SSI)) %>%
    left_join(meta_small, by = "sample_id") %>%
    mutate(
      SSI = factor(SSI, levels = ssi_levels, ordered = TRUE)
    )
  
  group_grid <- meta_small %>%
    distinct(data_scope, country, year, year_label, timespan, x_id) %>%
    crossing(
      SSI = factor(ssi_levels, levels = ssi_levels, ordered = TRUE)
    )
  
  richness_summary <- ssi_long %>%
    distinct(data_scope, country, year, year_label, timespan, x_id, SSI, Taxon) %>%
    count(
      data_scope, country, year, year_label, timespan, x_id, SSI,
      name = "n_species"
    ) %>%
    right_join(
      group_grid,
      by = c("data_scope", "country", "year", "year_label", "timespan", "x_id", "SSI")
    ) %>%
    mutate(
      n_species = replace_na(n_species, 0L),
      SSI_plot = factor(
        paste0("SSI_", SSI),
        levels = c("SSI_C", "SSI_B", "SSI_A")
      )
    ) %>%
    arrange(data_scope, country, year, SSI)
  
  richness_total <- richness_summary %>%
    group_by(data_scope, country, year, year_label, timespan, x_id) %>%
    summarise(
      total_ssi_species = sum(n_species),
      .groups = "drop"
    )
  
  richness_wide <- richness_summary %>%
    mutate(SSI = paste0("SSI_", SSI)) %>%
    select(data_scope, country, year, year_label, timespan, x_id, SSI, n_species) %>%
    pivot_wider(
      names_from = SSI,
      values_from = n_species,
      values_fill = 0
    ) %>%
    left_join(
      richness_total,
      by = c("data_scope", "country", "year", "year_label", "timespan", "x_id")
    ) %>%
    arrange(data_scope, country, year)
  
  sample_ssi_grid <- meta_small %>%
    select(data_scope, sample_id, country, forest, year, year_label, timespan, x_id) %>%
    crossing(
      SSI = factor(ssi_levels, levels = ssi_levels, ordered = TRUE)
    )
  
  sample_occupancy <- ssi_long %>%
    distinct(sample_id, SSI) %>%
    mutate(occupied = 1L) %>%
    right_join(
      sample_ssi_grid,
      by = c("sample_id", "SSI")
    ) %>%
    mutate(
      occupied = replace_na(occupied, 0L)
    )
  
  occupancy_summary <- sample_occupancy %>%
    group_by(data_scope, country, year, year_label, timespan, x_id, SSI) %>%
    summarise(
      n_samples_total    = n_distinct(sample_id),
      n_samples_occupied = sum(occupied == 1L, na.rm = TRUE),
      occupancy_freq     = n_samples_occupied / n_samples_total,
      .groups = "drop"
    ) %>%
    mutate(
      SSI = factor(SSI, levels = c("C", "B", "A"))
    ) %>%
    arrange(data_scope, country, year, SSI)
  
  list(
    data_scope = data_scope,
    taxa_in_matrix = taxa_in_matrix,
    ssi_long = ssi_long,
    richness_summary = richness_summary,
    richness_wide = richness_wide,
    richness_total = richness_total,
    sample_occupancy = sample_occupancy,
    occupancy_summary = occupancy_summary
  )
}

# =========================================================
# 3. RUN ANALYSIS ON BOTH INPUT SETS
# =========================================================

ssi_full <- analyse_one_ssi_set(
  X = X_use_full,
  meta = meta_use_full,
  data_scope = "Full target",
  ssi_taxa = ssi_taxa
)

ssi_subset <- analyse_one_ssi_set(
  X = X_use_subset,
  meta = meta_use_subset,
  data_scope = "Subsetted",
  ssi_taxa = ssi_taxa
)

ssi_richness_long <- bind_rows(
  ssi_full$richness_summary,
  ssi_subset$richness_summary
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    SSI_plot   = factor(SSI_plot, levels = c("SSI_C", "SSI_B", "SSI_A")),
    x_id       = factor(x_id, levels = x_id_levels)
  )

ssi_richness_wide <- bind_rows(
  ssi_full$richness_wide,
  ssi_subset$richness_wide
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    x_id       = factor(x_id, levels = x_id_levels)
  )

ssi_richness_total <- bind_rows(
  ssi_full$richness_total,
  ssi_subset$richness_total
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    x_id       = factor(x_id, levels = x_id_levels)
  )

ssi_occupancy_summary <- bind_rows(
  ssi_full$occupancy_summary,
  ssi_subset$occupancy_summary
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    SSI        = factor(SSI, levels = c("C", "B", "A")),
    SSI_plot   = factor(paste0("SSI_", SSI), levels = c("SSI_C", "SSI_B", "SSI_A")),
    x_id       = factor(x_id, levels = x_id_levels)
  )

cat("\n==============================\n")
cat("SSI richness summary, long\n")
cat("==============================\n")
print(ssi_richness_long, n = Inf)

cat("\n==============================\n")
cat("SSI richness summary, wide\n")
cat("==============================\n")
print(ssi_richness_wide, n = Inf)

cat("\n==============================\n")
cat("SSI occupancy summary\n")
cat("==============================\n")
print(ssi_occupancy_summary, n = Inf)

write_csv2(
  ssi_richness_long,
  file.path(out_dir, "ssi_richness_long_full_vs_subsetted.csv")
)

write_csv2(
  ssi_richness_wide,
  file.path(out_dir, "ssi_richness_wide_full_vs_subsetted.csv")
)

write_csv2(
  ssi_occupancy_summary,
  file.path(out_dir, "ssi_occupancy_full_vs_subsetted.csv")
)

# =========================================================
# 4. SHARED PLOT SETTINGS
# =========================================================

ssi_fill_scale <- scale_fill_manual(
  values = ssi_fill_vals,
  breaks = c("SSI_A", "SSI_B", "SSI_C"),
  labels = c("A", "B", "C"),
  name = "SSI category"
)

theme_ssi <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, lineheight = 0.9),
    legend.position = "right"
  )

theme_ssi_panel <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, lineheight = 0.9),
    legend.position = "bottom"
  )

# =========================================================
# 5. SSI ABSOLUTE RICHNESS PLOT
# =========================================================

p_ssi_richness_abs <- ggplot(
  ssi_richness_long,
  aes(x = x_id, y = n_species, fill = SSI_plot)
) +
  geom_col(
    width = 0.74,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    data = ssi_richness_total,
    aes(x = x_id, y = total_ssi_species + 0.8, label = total_ssi_species),
    inherit.aes = FALSE,
    size = 3
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  labs(
    x = NULL,
    y = "Number of SSI species"
  ) +
  theme_ssi

print(p_ssi_richness_abs)

ggsave(
  file.path(out_dir, "ssi_richness_absolute_full_vs_subsetted.png"),
  plot = p_ssi_richness_abs,
  width = 12,
  height = 7,
  dpi = 400
)

# =========================================================
# 6. SSI RELATIVE RICHNESS COMPOSITION PLOT
# =========================================================

ssi_richness_freq <- ssi_richness_long %>%
  group_by(data_scope, country, year, year_label, timespan, x_id) %>%
  mutate(
    total_ssi_species = sum(n_species, na.rm = TRUE),
    freq = if_else(total_ssi_species > 0, n_species / total_ssi_species, NA_real_),
    pct_label = if_else(
      !is.na(freq) & freq >= 0.06,
      paste0(round(100 * freq), "%"),
      NA_character_
    )
  ) %>%
  ungroup()

write_csv2(
  ssi_richness_freq,
  file.path(out_dir, "ssi_richness_relative_frequency_full_vs_subsetted.csv")
)

p_ssi_richness_freq <- ggplot(
  ssi_richness_freq,
  aes(x = x_id, y = freq, fill = SSI_plot)
) +
  geom_col(
    width = 0.74,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = pct_label),
    position = position_stack(vjust = 0.5),
    size = 2.8,
    na.rm = TRUE
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = NULL,
    y = "Relative frequency of SSI species"
  ) +
  theme_ssi

print(p_ssi_richness_freq)

ggsave(
  file.path(out_dir, "ssi_richness_relative_frequency_full_vs_subsetted.png"),
  plot = p_ssi_richness_freq,
  width = 12,
  height = 7,
  dpi = 400
)

# =========================================================
# 7. SSI ABSOLUTE OCCUPANCY PLOT
# Number of logs/samples containing at least one SSI taxon
# =========================================================

p_ssi_occupancy_abs <- ggplot(
  ssi_occupancy_summary,
  aes(x = x_id, y = n_samples_occupied, fill = SSI_plot)
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = n_samples_occupied),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 3
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  labs(
    x = NULL,
    y = "Number of occupied logs/samples"
  ) +
  theme_ssi

print(p_ssi_occupancy_abs)

ggsave(
  file.path(out_dir, "ssi_occupancy_absolute_full_vs_subsetted.png"),
  plot = p_ssi_occupancy_abs,
  width = 12,
  height = 7,
  dpi = 400
)

# =========================================================
# 8. SSI RELATIVE OCCUPANCY PLOT
# Proportion of logs/samples containing at least one SSI taxon
# =========================================================

p_ssi_occupancy_freq <- ggplot(
  ssi_occupancy_summary,
  aes(x = x_id, y = occupancy_freq, fill = SSI_plot)
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = paste0(round(100 * occupancy_freq), "%")),
    position = position_dodge(width = 0.78),
    vjust = -0.35,
    size = 3
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Proportion of occupied logs/samples"
  ) +
  theme_ssi

print(p_ssi_occupancy_freq)

ggsave(
  file.path(out_dir, "ssi_occupancy_relative_frequency_full_vs_subsetted.png"),
  plot = p_ssi_occupancy_freq,
  width = 12,
  height = 7,
  dpi = 400
)

# =========================================================
# 9. OPTIONAL SIMPLE COUNTRY X TIMESPAN OCCUPANCY TABLE
# Collapses year_label but keeps data_scope
# =========================================================

ssi_occupancy_country_timespan <- ssi_occupancy_summary %>%
  group_by(data_scope, country, timespan, SSI) %>%
  summarise(
    n_samples_total    = sum(n_samples_total),
    n_samples_occupied = sum(n_samples_occupied),
    occupancy_freq     = n_samples_occupied / n_samples_total,
    .groups = "drop"
  ) %>%
  arrange(data_scope, country, timespan, SSI)

cat("\n==============================\n")
cat("Collapsed SSI occupancy by data scope x country x timespan\n")
cat("==============================\n")
print(ssi_occupancy_country_timespan, n = Inf)

write_csv2(
  ssi_occupancy_country_timespan,
  file.path(out_dir, "ssi_occupancy_country_timespan_collapsed.csv")
)

# =========================================================
# 10. FINAL TWO-PANEL FIGURE
# A = SSI richness
# B = SSI occupancy frequency
# No numbers above bars
# =========================================================

p_ssi_richness_abs2 <- ggplot(
  ssi_richness_long,
  aes(x = x_id, y = n_species, fill = SSI_plot)
) +
  geom_col(
    width = 0.74,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  labs(
    title = "SSI richness",
    x = NULL,
    y = "Number of SSI species"
  ) +
  theme_ssi_panel

p_ssi_occupancy_freq2 <- ggplot(
  ssi_occupancy_summary,
  aes(x = x_id, y = occupancy_freq, fill = SSI_plot)
) +
  geom_col(
    position = position_dodge(width = 0.78),
    width = 0.68,
    colour = "black",
    linewidth = 0.25
  ) +
  facet_grid(
    data_scope ~ .,
    scales = "free_x",
    space = "free_x"
  ) +
  ssi_fill_scale +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "SSI occupancy",
    x = NULL,
    y = "Proportion of occupied logs/samples"
  ) +
  theme_ssi_panel

p_ssi_panel_freq_richness <- (
  p_ssi_richness_abs2 | p_ssi_occupancy_freq2
) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1)
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

print(p_ssi_panel_freq_richness)

ggsave(
  filename = file.path(out_dir, "ssi_panel_occupancy_frequency_richness_shared_legend_transparent.png"),
  plot = p_ssi_panel_freq_richness,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "transparent"
)





# =========================================================
# 3B. TAXA EXPLAINING SSI COUNT DIFFERENCES BETWEEN TIMESPANS
# Historical-only and Current-only SSI taxa per country
# =========================================================

ssi_taxa_by_country_timespan <- bind_rows(
  ssi_full$ssi_long,
  ssi_subset$ssi_long
) %>%
  distinct(
    data_scope,
    country,
    timespan,
    SSI,
    Taxon
  ) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    SSI        = factor(SSI, levels = ssi_levels, ordered = TRUE),
    present    = TRUE
  )

ssi_taxa_timespan_comparison <- ssi_taxa_by_country_timespan %>%
  tidyr::complete(
    data_scope,
    country,
    SSI,
    Taxon,
    timespan = factor(timespan_levels, levels = timespan_levels),
    fill = list(present = FALSE)
  ) %>%
  pivot_wider(
    names_from = timespan,
    values_from = present,
    values_fill = FALSE
  ) %>%
  mutate(
    previous = replace_na(previous, FALSE),
    current  = replace_na(current, FALSE),
    change_status = case_when(
      previous & !current ~ "Historical only",
      !previous & current ~ "Current only",
      previous & current  ~ "Shared",
      TRUE                ~ "Absent"
    ),
    direction = case_when(
      change_status == "Historical only" ~ "decrease",
      change_status == "Current only"    ~ "increase",
      change_status == "Shared"          ~ "shared",
      TRUE                               ~ "absent"
    )
  ) %>%
  filter(change_status %in% c("Historical only", "Current only", "Shared")) %>%
  arrange(data_scope, country, SSI, change_status, Taxon)

ssi_taxa_differ_between_timespans <- ssi_taxa_timespan_comparison %>%
  filter(change_status %in% c("Historical only", "Current only")) %>%
  arrange(data_scope, country, SSI, change_status, Taxon)

ssi_taxa_difference_summary <- ssi_taxa_differ_between_timespans %>%
  count(
    data_scope,
    country,
    SSI,
    change_status,
    name = "n_taxa"
  ) %>%
  tidyr::complete(
    data_scope,
    country,
    SSI,
    change_status = c("Historical only", "Current only"),
    fill = list(n_taxa = 0L)
  ) %>%
  arrange(data_scope, country, SSI, change_status)

cat("\n==============================\n")
cat("SSI taxa differing between Historical and Current\n")
cat("Historical only = taxa recorded in Historical but not Current\n")
cat("Current only    = taxa recorded in Current but not Historical\n")
cat("==============================\n")
print(ssi_taxa_differ_between_timespans, n = Inf)

cat("\n==============================\n")
cat("Summary of SSI taxa differing between Historical and Current\n")
cat("==============================\n")
print(ssi_taxa_difference_summary, n = Inf)

write_csv2(
  ssi_taxa_timespan_comparison,
  file.path(out_dir, "ssi_taxa_timespan_comparison_all_shared_and_different.csv")
)

write_csv2(
  ssi_taxa_differ_between_timespans,
  file.path(out_dir, "ssi_taxa_differ_between_historical_current.csv")
)

write_csv2(
  ssi_taxa_difference_summary,
  file.path(out_dir, "ssi_taxa_difference_summary_historical_current.csv")
)

# Optional: compact printed lists, easier to inspect in the terminal
ssi_taxa_difference_lists <- ssi_taxa_differ_between_timespans %>%
  group_by(data_scope, country, SSI, change_status) %>%
  summarise(
    n_taxa = n_distinct(Taxon),
    taxa = paste(sort(unique(Taxon)), collapse = " | "),
    .groups = "drop"
  ) %>%
  arrange(data_scope, country, SSI, change_status)

cat("\n==============================\n")
cat("Compact taxon lists by data scope x country x SSI x direction\n")
cat("==============================\n")
print(ssi_taxa_difference_lists, n = Inf)

write_csv2(
  ssi_taxa_difference_lists,
  file.path(out_dir, "ssi_taxa_difference_lists_historical_current.csv")
)


# =========================================================
# 3C. SPECIES-LEVEL SSI GAIN/LOSS TABLE BY COUNTRY
# Rows = taxa
# Columns = countries
# Cell format:
#   F +(delta); S +(delta)
#   F -(delta); S -(delta)
# where delta = current occurrences - historical occurrences
# =========================================================

ssi_occurrence_country_timespan <- bind_rows(
  ssi_full$ssi_long,
  ssi_subset$ssi_long
) %>%
  mutate(
    data_scope = as.character(data_scope),
    country    = as.character(country),
    timespan   = as.character(timespan),
    SSI        = as.character(SSI),
    Taxon      = as.character(Taxon)
  ) %>%
  group_by(data_scope, country, timespan, SSI, Taxon) %>%
  summarise(
    n_occurrences = n_distinct(sample_id),
    .groups = "drop"
  ) %>%
  complete(
    data_scope = c("Full target", "Subsetted"),
    country = country_levels,
    timespan = timespan_levels,
    nesting(SSI, Taxon),
    fill = list(n_occurrences = 0L)
  ) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
    country    = factor(country, levels = country_levels),
    timespan   = factor(timespan, levels = timespan_levels),
    SSI        = factor(SSI, levels = ssi_levels, ordered = TRUE)
  )

ssi_gain_loss_long <- ssi_occurrence_country_timespan %>%
  pivot_wider(
    names_from = timespan,
    values_from = n_occurrences,
    values_fill = 0
  ) %>%
  mutate(
    previous = replace_na(previous, 0L),
    current  = replace_na(current, 0L),
    delta_occurrences = current - previous,
    previous_present = previous > 0,
    current_present  = current > 0,
    change_status = case_when(
      !previous_present & current_present  ~ "gain",
      previous_present & !current_present  ~ "loss",
      previous_present & current_present   ~ "shared",
      TRUE                                 ~ "absent"
    ),
    sign = case_when(
      change_status == "gain" ~ "+",
      change_status == "loss" ~ "-",
      TRUE                    ~ ""
    ),
    cell_part = case_when(
      change_status %in% c("gain", "loss") & data_scope == "Full target" ~
        paste0("F ", sign, "(", delta_occurrences, ")"),
      change_status %in% c("gain", "loss") & data_scope == "Subsetted" ~
        paste0("S ", sign, "(", delta_occurrences, ")"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(change_status %in% c("gain", "loss")) %>%
  arrange(data_scope, country, SSI, Taxon)

ssi_gain_loss_country_cells <- ssi_gain_loss_long %>%
  select(Taxon, SSI, country, data_scope, cell_part) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted"))
  ) %>%
  arrange(Taxon, SSI, country, data_scope) %>%
  group_by(Taxon, SSI, country) %>%
  summarise(
    cell = paste(na.omit(cell_part), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    cell = na_if(cell, "")
  )

ssi_gain_loss_country_table <- ssi_gain_loss_country_cells %>%
  pivot_wider(
    names_from = country,
    values_from = cell,
    values_fill = ""
  ) %>%
  mutate(
    SSI = factor(SSI, levels = ssi_levels, ordered = TRUE)
  ) %>%
  arrange(SSI, Taxon)

cat("\n==============================\n")
cat("SSI species gain/loss table by country\n")
cat("==============================\n")
cat("Cell format: F +(delta); S +(delta)\n")
cat("F = full target dataset; S = subsetted dataset\n")
cat("+ = taxon gained in Current relative to Historical\n")
cat("- = taxon lost in Current relative to Historical\n")
cat("delta = Current occupied logs/samples minus Historical occupied logs/samples\n")
cat("==============================\n")
print(ssi_gain_loss_country_table, n = Inf)

write_csv2(
  ssi_gain_loss_long,
  file.path(out_dir, "ssi_gain_loss_long_country_timespan_full_subsetted.csv")
)

write_csv2(
  ssi_gain_loss_country_table,
  file.path(out_dir, "ssi_gain_loss_species_by_country_table.csv")
)

# =========================================================
# SSI ALL-TAXA TILE PLOT
# Full target + subsetted panel
# Gain/loss based on presence-absence
# Labels show occurrence delta (Current - Historical)
# Taxa absent everywhere are removed
# =========================================================

# ---- 1. Helper to build plot dataframe for one dataset ----

make_ssi_all_taxa_plot_df <- function(ssi_obj, dataset_label, ssi_taxa) {
  
  # Occurrence table
  ssi_occurrence_country_timespan <- ssi_obj$ssi_long %>%
    mutate(
      country  = as.character(country),
      timespan = as.character(timespan),
      SSI      = as.character(SSI),
      Taxon    = as.character(Taxon)
    ) %>%
    group_by(country, timespan, SSI, Taxon) %>%
    summarise(
      n_occurrences = n_distinct(sample_id),
      .groups = "drop"
    )
  
  # Master grid with all taxa in all country x timespan combinations
  ssi_master_grid <- tidyr::crossing(
    country  = country_levels,
    timespan = timespan_levels,
    ssi_taxa %>%
      mutate(
        SSI   = as.character(SSI),
        Taxon = as.character(Taxon)
      )
  )
  
  # Build full plotting table
  plot_df <- ssi_master_grid %>%
    left_join(
      ssi_occurrence_country_timespan,
      by = c("country", "timespan", "SSI", "Taxon")
    ) %>%
    mutate(
      n_occurrences = replace_na(n_occurrences, 0L)
    ) %>%
    pivot_wider(
      names_from = timespan,
      values_from = n_occurrences,
      values_fill = 0
    ) %>%
    mutate(
      previous = replace_na(previous, 0L),
      current  = replace_na(current, 0L),
      delta_occurrences = current - previous,
      change_status = case_when(
        previous == 0 & current  > 0 ~ "Gain",
        previous  > 0 & current == 0 ~ "Loss",
        previous  > 0 & current  > 0 ~ "Present both",
        previous == 0 & current == 0 ~ "Absent both",
        TRUE ~ NA_character_
      ),
      label = case_when(
        change_status == "Absent both" ~ NA_character_,
        delta_occurrences > 0          ~ paste0("+", delta_occurrences),
        delta_occurrences < 0          ~ as.character(delta_occurrences),
        delta_occurrences == 0         ~ "0",
        TRUE                           ~ NA_character_
      ),
      data_scope = dataset_label,
      country = factor(country, levels = country_levels),
      SSI = factor(SSI, levels = c("A", "B", "C"), ordered = TRUE),
      change_status = factor(
        change_status,
        levels = c("Loss", "Present both", "Gain", "Absent both")
      )
    )
  
  plot_df
}

# ---- 2. Build full + subsetted plot data ----

ssi_all_taxa_full_plot_df <- make_ssi_all_taxa_plot_df(
  ssi_obj = ssi_full,
  dataset_label = "Full target",
  ssi_taxa = ssi_taxa
)

ssi_all_taxa_subset_plot_df <- make_ssi_all_taxa_plot_df(
  ssi_obj = ssi_subset,
  dataset_label = "Subsetted",
  ssi_taxa = ssi_taxa
)

ssi_all_taxa_plot_df <- bind_rows(
  ssi_all_taxa_full_plot_df,
  ssi_all_taxa_subset_plot_df
) %>%
  mutate(
    data_scope = factor(data_scope, levels = c("Full target", "Subsetted"))
  )

# ---- 3. Remove taxa absent from both Historical and Current everywhere ----
# Keeps taxa if they occur at least once in either previous or current
# in at least one country and in at least one dataset.

taxa_observed_anywhere <- ssi_all_taxa_plot_df %>%
  group_by(Taxon) %>%
  summarise(
    total_occurrences = sum(previous, current, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_occurrences > 0) %>%
  pull(Taxon)

ssi_all_taxa_plot_df <- ssi_all_taxa_plot_df %>%
  filter(Taxon %in% taxa_observed_anywhere)

cat("\n==============================\n")
cat("SSI taxa retained in all-taxa tile plot\n")
cat("==============================\n")
cat("Retained taxa:", length(unique(ssi_all_taxa_plot_df$Taxon)), "\n")
cat("Removed taxa absent everywhere:",
    nrow(ssi_taxa %>% distinct(Taxon)) - length(unique(ssi_all_taxa_plot_df$Taxon)),
    "\n")

# ---- 4. Order retained taxa alphabetically within SSI ----

taxon_order <- ssi_all_taxa_plot_df %>%
  distinct(SSI, Taxon) %>%
  mutate(
    SSI   = factor(as.character(SSI), levels = c("A", "B", "C"), ordered = TRUE),
    Taxon = as.character(Taxon)
  ) %>%
  arrange(SSI, Taxon) %>%
  pull(Taxon)

ssi_all_taxa_plot_df <- ssi_all_taxa_plot_df %>%
  mutate(
    Taxon = factor(as.character(Taxon), levels = rev(taxon_order))
  )

# ---- 5. Plot ----
turnover_cols <- c(
  "Loss"         = "#8073AC",  # soft purple
  "Present both" = "grey85",
  "Gain"         = "#FDB863",  # soft orange
  "Absent both"  = "white"
)
p_ssi_all_taxa_tile <- ggplot(
  ssi_all_taxa_plot_df,
  aes(
    x = country,
    y = Taxon,
    fill = change_status
  )
) +
  geom_tile(
    colour = "grey35",
    linewidth = 0.25,
    width = 0.92,
    height = 0.92
  ) +
  geom_text(
    aes(label = label),
    size = 2.3,
    colour = "black",
    na.rm = TRUE
  ) +
  facet_grid(
    SSI ~ data_scope,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Loss"         = "#8073AC",
      "Present both" = "grey85",
      "Gain"         = "#FDB863",
      "Absent both"  = "white"
    ),
    breaks = c("Gain", "Loss", "Present both", "Absent both"),
    labels = c(
      "Gain",
      "Loss",
      "Present in both",
      "Absent in both"
    ),
    name = "Presence-absence status",
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "SSI species turnover between Historical and Current surveys") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7, face = "italic"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8)
  )

print(p_ssi_all_taxa_tile)
ggsave(
  filename = file.path(out_dir, "ssi_all_taxa_species_country_tileplot_full_and_subsetted_turnover_status_A4.png"),
  plot = p_ssi_all_taxa_tile,
  width = 7.25,
  height = 10.7,
  units = "in",
  dpi = 600
)

###########
##########
######## SSI TO ENV VARS
# =========================================================
# 11. SUBSTRATE ASSOCIATIONS OF SSI TAXA
# Which substrate characteristics are most important for SSI taxa?
# Responses:
#   1. SSI presence per deadwood unit, binomial GLM
#   2. SSI richness per deadwood unit, negative-binomial GLM where possible
#      with Poisson fallback
# Models are run for:
#   - Full target
#   - Subsetted
#   - Global models with forest and survey period adjustment
#   - Forest-level models with survey period adjustment
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(MASS)
})

ssi_substrate_out_dir <- file.path(out_dir, "ssi_substrate_models")
dir.create(ssi_substrate_out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 11.1 Helper: build sample-level SSI response table
# ---------------------------------------------------------

make_ssi_substrate_model_data <- function(X, meta, data_scope, ssi_taxa) {
  
  stopifnot(nrow(X) == nrow(meta))
  stopifnot("sample_id" %in% names(meta))
  
  if (!identical(rownames(X), meta$sample_id)) {
    rownames(X) <- meta$sample_id
  }
  
  taxa_in_matrix <- intersect(ssi_taxa$Taxon, colnames(X))
  
  if (length(taxa_in_matrix) == 0) {
    stop("No SSI taxa found in matrix columns for ", data_scope)
  }
  
  ssi_lookup <- ssi_taxa %>%
    dplyr::distinct(Taxon, SSI) %>%
    dplyr::mutate(
      Taxon = as.character(Taxon),
      SSI   = as.character(SSI)
    )
  
  X_ssi <- X %>%
    as.data.frame(check.names = FALSE) %>%
    dplyr::select(dplyr::all_of(taxa_in_matrix)) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.integer(.x > 0)))
  
  ssi_richness_total <- rowSums(X_ssi, na.rm = TRUE)
  
  ssi_richness_by_category <- X_ssi %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample_id = meta$sample_id) %>%
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "Taxon",
      values_to = "pa"
    ) %>%
    dplyr::left_join(ssi_lookup, by = "Taxon") %>%
    dplyr::filter(!is.na(SSI)) %>%
    dplyr::group_by(sample_id, SSI) %>%
    dplyr::summarise(
      n_ssi_cat = sum(pa > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      sample_id = as.character(meta$sample_id),
      SSI = ssi_levels,
      fill = list(n_ssi_cat = 0L)
    ) %>%
    dplyr::mutate(
      SSI = paste0("ssi_richness_", SSI)
    ) %>%
    tidyr::pivot_wider(
      names_from = SSI,
      values_from = n_ssi_cat,
      values_fill = 0
    )
  
  dat <- meta %>%
    dplyr::mutate(
      sample_id = as.character(sample_id),
      country   = as.character(country),
      forest    = as.character(forest),
      timespan  = as.character(timespan),
      year      = as.integer(as.character(year))
    ) %>%
    dplyr::mutate(
      forest = dplyr::recode(forest, "Strodam" = "Strødam"),
      data_scope = data_scope,
      ssi_richness_total = as.integer(ssi_richness_total),
      ssi_presence_total = as.integer(ssi_richness_total > 0)
    ) %>%
    dplyr::left_join(ssi_richness_by_category, by = "sample_id")
  
  for (nm in c("ssi_richness_A", "ssi_richness_B", "ssi_richness_C")) {
    if (!nm %in% names(dat)) dat[[nm]] <- 0L
    dat[[nm]] <- tidyr::replace_na(as.integer(dat[[nm]]), 0L)
  }
  
  dat <- dat %>%
    dplyr::mutate(
      ssi_presence_A = as.integer(ssi_richness_A > 0),
      ssi_presence_B = as.integer(ssi_richness_B > 0),
      ssi_presence_C = as.integer(ssi_richness_C > 0)
    )
  
  required_env <- c(
    "log_av_ds",
    "log_diameter",
    "exposure",
    "mosscover",
    "snag"
  )
  
  missing_env <- setdiff(required_env, names(dat))
  
  if (length(missing_env) > 0) {
    stop(
      "Missing required substrate columns in ", data_scope, ": ",
      paste(missing_env, collapse = ", ")
    )
  }
  
  dat <- dat %>%
    dplyr::mutate(
      data_scope = factor(data_scope, levels = c("Full target", "Subsetted")),
      country    = factor(country, levels = country_levels),
      forest     = factor(forest, levels = c("Sonian", "Rajhenav", "Strødam", "Suserup")),
      timespan   = factor(timespan, levels = timespan_levels),
      exposure   = factor(stringr::str_trim(as.character(exposure))),
      snag       = factor(as.integer(snag > 0), levels = c(0, 1), labels = c("log", "snag")),
      
      log_av_ds      = as.numeric(log_av_ds),
      log_diameter   = as.numeric(log_diameter),
      mosscover      = as.numeric(mosscover),
      log_diameter_l = log1p(log_diameter),
      
      log_av_ds_z    = as.numeric(scale(log_av_ds)),
      log_diameter_z = as.numeric(scale(log_diameter_l)),
      mosscover_z    = as.numeric(scale(mosscover))
    )
  
  if ("shade" %in% levels(dat$exposure)) {
    dat <- dat %>%
      dplyr::mutate(exposure = stats::relevel(exposure, ref = "shade"))
  }
  
  dat
}
ssi_substrate_full <- make_ssi_substrate_model_data(
  X = X_use_full,
  meta = meta_use_full,
  data_scope = "Full target",
  ssi_taxa = ssi_taxa
)

ssi_substrate_subset <- make_ssi_substrate_model_data(
  X = X_use_subset,
  meta = meta_use_subset,
  data_scope = "Subsetted",
  ssi_taxa = ssi_taxa
)

ssi_substrate_data <- bind_rows(
  ssi_substrate_full,
  ssi_substrate_subset
)

cat("\n==============================\n")
cat("SSI substrate model data summary\n")
cat("==============================\n")

ssi_substrate_data_summary <- ssi_substrate_data %>%
  group_by(data_scope, forest, timespan) %>%
  summarise(
    n_samples = n(),
    n_samples_with_ssi = sum(ssi_presence_total == 1L, na.rm = TRUE),
    mean_ssi_richness = mean(ssi_richness_total, na.rm = TRUE),
    max_ssi_richness = max(ssi_richness_total, na.rm = TRUE),
    .groups = "drop"
  )

print(ssi_substrate_data_summary, n = Inf)

write_csv2(
  ssi_substrate_data,
  file.path(ssi_substrate_out_dir, "ssi_substrate_model_input_sample_level.csv")
)

write_csv2(
  ssi_substrate_data_summary,
  file.path(ssi_substrate_out_dir, "ssi_substrate_model_input_summary.csv")
)

# ---------------------------------------------------------
# 11.2 Model helpers
# ---------------------------------------------------------

substrate_terms <- c(
  "log_av_ds_z",
  "log_diameter_z",
  "exposure",
  "mosscover_z",
  "snag"
)

clean_drop1_table <- function(drop_obj, response, model_type, data_scope, forest_scope) {
  
  as.data.frame(drop_obj) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    filter(term != "<none>") %>%
    rename_with(~ str_replace_all(.x, " ", "_")) %>%
    mutate(
      response = response,
      model_type = model_type,
      data_scope = data_scope,
      forest_scope = forest_scope,
      term_group = case_when(
        str_detect(term, "^log_av_ds_z$")    ~ "Decay stage",
        str_detect(term, "^log_diameter_z$") ~ "Log diameter",
        str_detect(term, "^exposure$")       ~ "Exposure",
        str_detect(term, "^mosscover_z$")    ~ "Moss cover",
        str_detect(term, "^snag$")           ~ "Snag status",
        str_detect(term, "^forest$")         ~ "Forest",
        str_detect(term, "^timespan$")       ~ "Survey period",
        TRUE                                 ~ term
      )
    ) %>%
    relocate(response, model_type, data_scope, forest_scope, term, term_group)
}

extract_glm_coefficients <- function(model, response, model_type, data_scope, forest_scope) {
  
  sm <- summary(model)
  
  as.data.frame(sm$coefficients) %>%
    rownames_to_column("coefficient") %>%
    as_tibble() %>%
    rename_with(~ str_replace_all(.x, " ", "_")) %>%
    mutate(
      response = response,
      model_type = model_type,
      data_scope = data_scope,
      forest_scope = forest_scope
    ) %>%
    relocate(response, model_type, data_scope, forest_scope, coefficient)
}

fit_binomial_ssi_model <- function(dat, response, data_scope, forest_scope, global_model = TRUE) {
  
  dat_model <- dat %>%
    filter(!is.na(.data[[response]])) %>%
    filter(
      !is.na(log_av_ds_z),
      !is.na(log_diameter_z),
      !is.na(exposure),
      !is.na(mosscover_z),
      !is.na(snag),
      !is.na(timespan)
    ) %>%
    droplevels()
  
  if (nrow(dat_model) < 20) {
    return(NULL)
  }
  
  if (length(unique(dat_model[[response]])) < 2) {
    return(NULL)
  }
  
  if (global_model) {
    model_formula <- as.formula(
      paste(
        response,
        "~ forest + timespan + log_av_ds_z + log_diameter_z + exposure + mosscover_z + snag"
      )
    )
  } else {
    model_formula <- as.formula(
      paste(
        response,
        "~ timespan + log_av_ds_z + log_diameter_z + exposure + mosscover_z + snag"
      )
    )
  }
  
  model <- glm(
    formula = model_formula,
    data = dat_model,
    family = binomial()
  )
  
  drop_tbl <- drop1(model, test = "Chisq") %>%
    clean_drop1_table(
      response = response,
      model_type = "Binomial GLM, SSI presence",
      data_scope = data_scope,
      forest_scope = forest_scope
    )
  
  coef_tbl <- extract_glm_coefficients(
    model = model,
    response = response,
    model_type = "Binomial GLM, SSI presence",
    data_scope = data_scope,
    forest_scope = forest_scope
  ) %>%
    mutate(
      odds_ratio = exp(Estimate)
    )
  
  model_summary <- tibble(
    response = response,
    model_type = "Binomial GLM, SSI presence",
    data_scope = data_scope,
    forest_scope = forest_scope,
    n_samples = nrow(dat_model),
    n_positive = sum(dat_model[[response]] == 1L, na.rm = TRUE),
    prevalence = mean(dat_model[[response]] == 1L, na.rm = TRUE),
    AIC = AIC(model),
    residual_deviance = deviance(model),
    null_deviance = model$null.deviance
  )
  
  list(
    model = model,
    drop = drop_tbl,
    coefficients = coef_tbl,
    summary = model_summary
  )
}

fit_count_ssi_model <- function(dat, response, data_scope, forest_scope, global_model = TRUE) {
  
  dat_model <- dat %>%
    dplyr::filter(!is.na(.data[[response]])) %>%
    dplyr::filter(
      !is.na(log_av_ds_z),
      !is.na(log_diameter_z),
      !is.na(exposure),
      !is.na(mosscover_z),
      !is.na(snag),
      !is.na(timespan)
    ) %>%
    droplevels()
  
  if (nrow(dat_model) < 20) {
    return(NULL)
  }
  
  if (sum(dat_model[[response]], na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  if (length(unique(dat_model[[response]])) < 2) {
    return(NULL)
  }
  
  if (global_model) {
    model_formula <- as.formula(
      paste(
        response,
        "~ forest + timespan + log_av_ds_z + log_diameter_z + exposure + mosscover_z + snag"
      )
    )
  } else {
    model_formula <- as.formula(
      paste(
        response,
        "~ timespan + log_av_ds_z + log_diameter_z + exposure + mosscover_z + snag"
      )
    )
  }
  
  model_type <- "Negative-binomial GLM, SSI richness"
  
  model <- tryCatch(
    {
      suppressWarnings(
        MASS::glm.nb(
          formula = model_formula,
          data = dat_model,
          control = glm.control(maxit = 100)
        )
      )
    },
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    model_type <- "Poisson GLM fallback, SSI richness"
    
    model <- glm(
      formula = model_formula,
      data = dat_model,
      family = poisson()
    )
  }
  
  drop_tbl <- drop1(model, test = "Chisq") %>%
    clean_drop1_table(
      response = response,
      model_type = model_type,
      data_scope = data_scope,
      forest_scope = forest_scope
    )
  
  coef_tbl <- extract_glm_coefficients(
    model = model,
    response = response,
    model_type = model_type,
    data_scope = data_scope,
    forest_scope = forest_scope
  ) %>%
    dplyr::mutate(
      incidence_rate_ratio = exp(Estimate)
    )
  
  model_summary <- tibble::tibble(
    response = response,
    model_type = model_type,
    data_scope = data_scope,
    forest_scope = forest_scope,
    n_samples = nrow(dat_model),
    total_ssi_records = sum(dat_model[[response]], na.rm = TRUE),
    mean_response = mean(dat_model[[response]], na.rm = TRUE),
    max_response = max(dat_model[[response]], na.rm = TRUE),
    AIC = AIC(model),
    residual_deviance = deviance(model),
    null_deviance = model$null.deviance
  )
  
  list(
    model = model,
    drop = drop_tbl,
    coefficients = coef_tbl,
    summary = model_summary
  )
}
run_ssi_substrate_models_one_dataset <- function(dat, data_scope) {
  
  dat_scope <- dat %>%
    filter(data_scope == !!data_scope) %>%
    droplevels()
  
  presence_responses <- c(
    "ssi_presence_total",
    "ssi_presence_A",
    "ssi_presence_B",
    "ssi_presence_C"
  )
  
  count_responses <- c(
    "ssi_richness_total",
    "ssi_richness_A",
    "ssi_richness_B",
    "ssi_richness_C"
  )
  
  model_results <- list()
  
  for (resp in presence_responses) {
    model_results[[paste(data_scope, "Global", resp, "presence", sep = "_")]] <-
      fit_binomial_ssi_model(
        dat = dat_scope,
        response = resp,
        data_scope = data_scope,
        forest_scope = "Global",
        global_model = TRUE
      )
  }
  
  for (resp in count_responses) {
    model_results[[paste(data_scope, "Global", resp, "count", sep = "_")]] <-
      fit_count_ssi_model(
        dat = dat_scope,
        response = resp,
        data_scope = data_scope,
        forest_scope = "Global",
        global_model = TRUE
      )
  }
  
  for (f in levels(droplevels(dat_scope$forest))) {
    
    dat_f <- dat_scope %>%
      filter(forest == f) %>%
      droplevels()
    
    for (resp in presence_responses) {
      model_results[[paste(data_scope, f, resp, "presence", sep = "_")]] <-
        fit_binomial_ssi_model(
          dat = dat_f,
          response = resp,
          data_scope = data_scope,
          forest_scope = f,
          global_model = FALSE
        )
    }
    
    for (resp in count_responses) {
      model_results[[paste(data_scope, f, resp, "count", sep = "_")]] <-
        fit_count_ssi_model(
          dat = dat_f,
          response = resp,
          data_scope = data_scope,
          forest_scope = f,
          global_model = FALSE
        )
    }
  }
  
  model_results <- model_results[!vapply(model_results, is.null, logical(1))]
  model_results
}

ssi_substrate_model_results <- c(
  run_ssi_substrate_models_one_dataset(ssi_substrate_data, "Full target"),
  run_ssi_substrate_models_one_dataset(ssi_substrate_data, "Subsetted")
)

# ---------------------------------------------------------
# 11.3 Combine model outputs
# ---------------------------------------------------------

ssi_substrate_drop_tests <- bind_rows(
  lapply(ssi_substrate_model_results, `[[`, "drop")
) %>%
  mutate(
    substrate_term = term_group %in% c(
      "Decay stage",
      "Log diameter",
      "Exposure",
      "Moss cover",
      "Snag status"
    )
  )

ssi_substrate_coefficients <- bind_rows(
  lapply(ssi_substrate_model_results, `[[`, "coefficients")
)

ssi_substrate_model_summaries <- bind_rows(
  lapply(ssi_substrate_model_results, `[[`, "summary")
)

cat("\n==============================\n")
cat("SSI substrate model summaries\n")
cat("==============================\n")
print(ssi_substrate_model_summaries, n = Inf)

cat("\n==============================\n")
cat("Drop-one tests for SSI substrate models\n")
cat("==============================\n")
print(ssi_substrate_drop_tests, n = Inf)

cat("\n==============================\n")
cat("SSI substrate model coefficients\n")
cat("==============================\n")
print(ssi_substrate_coefficients, n = Inf)

write_csv2(
  ssi_substrate_model_summaries,
  file.path(ssi_substrate_out_dir, "ssi_substrate_model_summaries.csv")
)

write_csv2(
  ssi_substrate_drop_tests,
  file.path(ssi_substrate_out_dir, "ssi_substrate_drop_one_tests.csv")
)

write_csv2(
  ssi_substrate_coefficients,
  file.path(ssi_substrate_out_dir, "ssi_substrate_model_coefficients.csv")
)

# ---------------------------------------------------------
# 11.4 Rank substrate importance
# ---------------------------------------------------------

ssi_substrate_importance <- ssi_substrate_drop_tests %>%
  filter(substrate_term) %>%
  mutate(
    importance_stat = case_when(
      "LRT" %in% names(.) & !is.na(LRT) ~ LRT,
      "Deviance" %in% names(.) & !is.na(Deviance) ~ Deviance,
      TRUE ~ NA_real_
    ),
    p_value = case_when(
      "Pr(>Chi)" %in% names(.) ~ .data[["Pr(>Chi)"]],
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(data_scope, forest_scope, response, model_type) %>%
  arrange(desc(importance_stat), .by_group = TRUE) %>%
  mutate(
    importance_rank = row_number()
  ) %>%
  ungroup()

cat("\n==============================\n")
cat("Ranked substrate importance for SSI models\n")
cat("==============================\n")
print(ssi_substrate_importance, n = Inf)

write_csv2(
  ssi_substrate_importance,
  file.path(ssi_substrate_out_dir, "ssi_substrate_importance_ranked.csv")
)

ssi_substrate_importance_compact <- ssi_substrate_importance %>%
  filter(response %in% c("ssi_presence_total", "ssi_richness_total")) %>%
  dplyr::select(
    data_scope,
    forest_scope,
    response,
    model_type,
    importance_rank,
    term_group,
    importance_stat,
    p_value
  ) %>%
  arrange(data_scope, response, forest_scope, importance_rank)

cat("\n==============================\n")
cat("Compact ranked substrate importance, total SSI responses\n")
cat("==============================\n")
print(ssi_substrate_importance_compact, n = Inf)

write_csv2(
  ssi_substrate_importance_compact,
  file.path(ssi_substrate_out_dir, "ssi_substrate_importance_total_responses_compact.csv")
)

# ---------------------------------------------------------
# 11.5 Plot substrate importance
# ---------------------------------------------------------

plot_importance_df <- ssi_substrate_importance %>%
  filter(
    response %in% c("ssi_presence_total", "ssi_richness_total"),
    forest_scope == "Global",
    !is.na(importance_stat)
  ) %>%
  mutate(
    response_label = case_when(
      response == "ssi_presence_total" ~ "SSI presence",
      response == "ssi_richness_total" ~ "SSI richness",
      TRUE ~ response
    ),
    term_group = factor(
      term_group,
      levels = c(
        "Decay stage",
        "Log diameter",
        "Exposure",
        "Moss cover",
        "Snag status"
      )
    )
  )

p_ssi_substrate_importance_global <- ggplot(
  plot_importance_df,
  aes(
    x = reorder(term_group, importance_stat),
    y = importance_stat
  )
) +
  geom_col(
    width = 0.72,
    colour = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  facet_grid(
    response_label ~ data_scope,
    scales = "free_x"
  ) +
  labs(
    x = NULL,
    y = "Drop-one test statistic",
    title = "Substrate importance for SSI taxa"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_ssi_substrate_importance_global)

ggsave(
  filename = file.path(ssi_substrate_out_dir, "ssi_substrate_importance_global_total_responses.png"),
  plot = p_ssi_substrate_importance_global,
  width = 8.0,
  height = 5.5,
  dpi = 600
)

# ---------------------------------------------------------
# 11.6 Plot predicted SSI responses along main continuous gradients
#      Global models, adjusted over other terms
# ---------------------------------------------------------

make_prediction_grid <- function(dat, data_scope, gradient_var) {
  
  dat_scope <- dat %>%
    filter(data_scope == !!data_scope) %>%
    droplevels()
  
  exposure_ref <- if ("shade" %in% levels(dat_scope$exposure)) {
    "shade"
  } else {
    levels(dat_scope$exposure)[1]
  }
  
  expand_grid(
    forest = levels(dat_scope$forest),
    timespan = levels(dat_scope$timespan),
    log_av_ds_z = if (gradient_var == "log_av_ds_z") {
      seq(
        min(dat_scope$log_av_ds_z, na.rm = TRUE),
        max(dat_scope$log_av_ds_z, na.rm = TRUE),
        length.out = 100
      )
    } else {
      0
    },
    log_diameter_z = if (gradient_var == "log_diameter_z") {
      seq(
        min(dat_scope$log_diameter_z, na.rm = TRUE),
        max(dat_scope$log_diameter_z, na.rm = TRUE),
        length.out = 100
      )
    } else {
      0
    },
    mosscover_z = if (gradient_var == "mosscover_z") {
      seq(
        min(dat_scope$mosscover_z, na.rm = TRUE),
        max(dat_scope$mosscover_z, na.rm = TRUE),
        length.out = 100
      )
    } else {
      0
    },
    exposure = exposure_ref,
    snag = "log"
  ) %>%
    mutate(
      data_scope = data_scope,
      forest = factor(forest, levels = levels(dat_scope$forest)),
      timespan = factor(timespan, levels = levels(dat_scope$timespan)),
      exposure = factor(exposure, levels = levels(dat_scope$exposure)),
      snag = factor(snag, levels = levels(dat_scope$snag))
    )
}

get_global_model <- function(results, data_scope, response_pattern) {
  
  ix <- names(results)[
    str_detect(names(results), fixed(data_scope)) &
      str_detect(names(results), "Global") &
      str_detect(names(results), response_pattern)
  ]
  
  if (length(ix) == 0) {
    return(NULL)
  }
  
  results[[ix[1]]]$model
}

predict_response_curve <- function(model, newdat, response_label, gradient_label, model_family) {
  
  if (is.null(model)) {
    return(NULL)
  }
  
  pred <- predict(
    model,
    newdata = newdat,
    type = "link",
    se.fit = TRUE
  )
  
  newdat %>%
    mutate(
      fit_link = pred$fit,
      se_link = pred$se.fit,
      fit = if (model_family == "binomial") {
        plogis(fit_link)
      } else {
        exp(fit_link)
      },
      lwr = if (model_family == "binomial") {
        plogis(fit_link - 1.96 * se_link)
      } else {
        exp(fit_link - 1.96 * se_link)
      },
      upr = if (model_family == "binomial") {
        plogis(fit_link + 1.96 * se_link)
      } else {
        exp(fit_link + 1.96 * se_link)
      },
      response_label = response_label,
      gradient_label = gradient_label,
      model_family = model_family
    )
}

prediction_curves <- list()

for (ds in c("Full target", "Subsetted")) {
  
  mod_presence <- get_global_model(
    results = ssi_substrate_model_results,
    data_scope = ds,
    response_pattern = "ssi_presence_total_presence"
  )
  
  mod_richness <- get_global_model(
    results = ssi_substrate_model_results,
    data_scope = ds,
    response_pattern = "ssi_richness_total_count"
  )
  
  for (grad in c("log_av_ds_z", "log_diameter_z", "mosscover_z")) {
    
    grid <- make_prediction_grid(
      dat = ssi_substrate_data,
      data_scope = ds,
      gradient_var = grad
    )
    
    gradient_label <- case_when(
      grad == "log_av_ds_z" ~ "Decay stage, standardized",
      grad == "log_diameter_z" ~ "Log diameter, standardized",
      grad == "mosscover_z" ~ "Moss cover, standardized",
      TRUE ~ grad
    )
    
    prediction_curves[[paste(ds, grad, "presence", sep = "_")]] <-
      predict_response_curve(
        model = mod_presence,
        newdat = grid,
        response_label = "SSI presence",
        gradient_label = gradient_label,
        model_family = "binomial"
      )
    
    prediction_curves[[paste(ds, grad, "richness", sep = "_")]] <-
      predict_response_curve(
        model = mod_richness,
        newdat = grid,
        response_label = "SSI richness",
        gradient_label = gradient_label,
        model_family = "count"
      )
  }
}

ssi_prediction_curves <- bind_rows(prediction_curves)

write_csv2(
  ssi_prediction_curves,
  file.path(ssi_substrate_out_dir, "ssi_substrate_prediction_curves_global_models.csv")
)

p_ssi_substrate_predictions <- ggplot(
  ssi_prediction_curves,
  aes(
    x = case_when(
      gradient_label == "Decay stage, standardized" ~ log_av_ds_z,
      gradient_label == "Log diameter, standardized" ~ log_diameter_z,
      gradient_label == "Moss cover, standardized" ~ mosscover_z,
      TRUE ~ NA_real_
    ),
    y = fit,
    group = interaction(forest, timespan),
    linetype = timespan
  )
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.12,
    linewidth = 0
  ) +
  geom_line(linewidth = 0.55) +
  facet_grid(
    response_label + data_scope ~ gradient_label,
    scales = "free_y"
  ) +
  labs(
    x = "Standardized substrate gradient",
    y = "Predicted response",
    linetype = "Survey period",
    title = "Predicted SSI responses along substrate gradients"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold", size = 8),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(p_ssi_substrate_predictions)

ggsave(
  filename = file.path(ssi_substrate_out_dir, "ssi_substrate_prediction_curves_global_models.png"),
  plot = p_ssi_substrate_predictions,
  width = 10.5,
  height = 8.0,
  dpi = 600
)

# ---------------------------------------------------------
# 11.7 Short terminal interpretation table
# ---------------------------------------------------------

ssi_top_substrate_terms <- ssi_substrate_importance %>%
  filter(
    substrate_term,
    response %in% c("ssi_presence_total", "ssi_richness_total"),
    importance_rank == 1
  ) %>%
  dplyr::select(
    data_scope,
    forest_scope,
    response,
    model_type,
    top_substrate = term_group,
    importance_stat,
    p_value
  ) %>%
  arrange(data_scope, response, forest_scope)

cat("\n==============================\n")
cat("Top-ranked substrate variable for total SSI responses\n")
cat("==============================\n")
print(ssi_top_substrate_terms, n = Inf)

write_csv2(
  ssi_top_substrate_terms,
  file.path(ssi_substrate_out_dir, "ssi_top_ranked_substrate_terms_total_responses.csv")
)
# =========================================================
# Inspect coefficient directions for main SSI substrate effects
# =========================================================

ssi_main_coefficients <- ssi_substrate_coefficients %>%
  dplyr::filter(
    response %in% c("ssi_presence_total", "ssi_richness_total"),
    coefficient %in% c(
      "log_av_ds_z",
      "log_diameter_z",
      "mosscover_z",
      "snagsnag"
    ) |
      stringr::str_detect(coefficient, "^exposure")
  ) %>%
  dplyr::mutate(
    effect_direction = dplyr::case_when(
      Estimate > 0 ~ "positive",
      Estimate < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    effect_multiplier = exp(Estimate)
  ) %>%
  dplyr::arrange(
    data_scope,
    response,
    forest_scope,
    coefficient
  )

cat("\n==============================\n")
cat("Coefficient directions for SSI substrate models\n")
cat("==============================\n")
print(ssi_main_coefficients, n = Inf)

readr::write_csv2(
  ssi_main_coefficients,
  file.path(
    ssi_substrate_out_dir,
    "ssi_substrate_main_coefficient_directions.csv"
  )
)