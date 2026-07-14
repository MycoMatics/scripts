# =========================================================
# TWO PARTIAL dbRDA SPIDER WORKFLOWS
#
# A. Residual temporal space:
#    community ~ T_Sonian + T_Rajhenav + T_Strodam + T_Suserup +
#      Condition(forest + substrate variables)
#
# B. Substrate-adjusted reference space:
#    community ~ forest + timespan + forest:timespan +
#      Condition(substrate variables)
#
# Required objects:
#   X_use_full, meta_use_full
#   X_use_subset, meta_use_subset
#
# Output:
#   PNG spider plots
#   Faceted 2-panel residual temporal-space figure
#   Faceted 2-panel substrate-adjusted reference-space figure
#   Faceted 4-panel dbRDA figure
#   CSV test tables, site scores, centroids, spider segments
#   General model summary table
#   General test summary table
#   Sonian-to-reference directional delta tables
# =========================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)
library(readr)
library(grid)

# ---------------------------------------------------------
# Settings
# ---------------------------------------------------------

# Use n_perm <- 9999 for final analyses.
# Use smaller values during testing if needed.
n_perm <- 20
n_boot_direction <- 20

out_dir <- "dbrda_combined_outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

forest_cols <- c(
  "Sonian"   = "#B8860B",
  "Rajhenav" = "#2E7D32",
  "Strødam"  = "#466C95",
  "Suserup"  = "#5FA8D3"
)

timespan_shapes <- c(
  "previous" = 21,
  "current"  = 24
)

timespan_labs <- c(
  "previous" = "Historical",
  "current"  = "Contemporary"
)

datasets <- list(
  full = list(
    X = X_use_full,
    meta = meta_use_full,
    label = "Full dataset",
    file_label = "full_dataset"
  ),
  substrate_matched = list(
    X = X_use_subset,
    meta = meta_use_subset,
    label = "Substrate-matched dataset",
    file_label = "substrate_matched_dataset"
  )
)

analysis_setups <- list(
  residual_temporal_space = list(
    label = "Residual temporal space",
    formula_type = "residual",
    x_prefix = "partial dbRDA",
    subtitle = "Residual temporal centroid shift after conditioning on forest and substrate variables",
    plot_file_suffix = "partial_dbrda_residual_centroids_spider",
    table_file_suffix = "partial_dbrda",
    use_strata = TRUE
  ),
  substrate_adjusted_reference_space = list(
    label = "Substrate-adjusted reference space",
    formula_type = "reference",
    x_prefix = "substrate-adjusted dbRDA",
    subtitle = "Substrate-adjusted centroid trajectories retaining forest structure",
    plot_file_suffix = "substrate_adjusted_reference_space_spider",
    table_file_suffix = "substrate_adjusted_dbrda",
    use_strata = FALSE
  )
)

# ---------------------------------------------------------
# Helper: prepare input matrix and metadata
# ---------------------------------------------------------

prepare_dbrda_inputs <- function(X_in, meta_in) {
  
  X_out <- as.data.frame(X_in)
  
  meta_out <- meta_in %>%
    as.data.frame() %>%
    mutate(
      sample_id = as.character(sample_id),
      forest = factor(forest, levels = c("Sonian", "Rajhenav", "Strødam", "Suserup")),
      timespan = factor(timespan, levels = c("previous", "current")),
      exposure = factor(exposure),
      snag = factor(snag),
      log_av_ds = as.numeric(log_av_ds),
      log_diameter = as.numeric(log_diameter),
      mosscover = as.numeric(mosscover)
    )
  
  rownames(X_out) <- meta_out$sample_id
  rownames(meta_out) <- meta_out$sample_id
  
  keep_complete <- complete.cases(
    meta_out[, c(
      "forest",
      "timespan",
      "log_av_ds",
      "log_diameter",
      "exposure",
      "mosscover",
      "snag"
    )]
  )
  
  X_out <- X_out[keep_complete, , drop = FALSE]
  meta_out <- meta_out[keep_complete, , drop = FALSE]
  
  X_out <- X_out[, colSums(X_out, na.rm = TRUE) > 0, drop = FALSE]
  
  keep_rich <- rowSums(X_out, na.rm = TRUE) >= 2
  
  X_out <- X_out[keep_rich, , drop = FALSE]
  meta_out <- meta_out[keep_rich, , drop = FALSE]
  
  meta_out <- meta_out %>%
    mutate(
      T_Sonian   = as.numeric(forest == "Sonian"   & timespan == "current"),
      T_Rajhenav = as.numeric(forest == "Rajhenav" & timespan == "current"),
      T_Strodam  = as.numeric(forest == "Strødam"  & timespan == "current"),
      T_Suserup  = as.numeric(forest == "Suserup"  & timespan == "current")
    )
  
  list(
    X = X_out,
    meta = meta_out
  )
}

# ---------------------------------------------------------
# Helper: safe axis percentages
# ---------------------------------------------------------

get_axis_percentages <- function(cap_mod) {
  
  eig <- eigenvals(cap_mod, model = "constrained")
  
  if (length(eig) == 0 || sum(eig, na.rm = TRUE) == 0) {
    return(
      tibble(
        CAP1_eigenvalue = NA_real_,
        CAP2_eigenvalue = NA_real_,
        CAP1_percent_of_constrained = NA_real_,
        CAP2_percent_of_constrained = NA_real_
      )
    )
  }
  
  eig_perc <- eig / sum(eig, na.rm = TRUE) * 100
  
  tibble(
    CAP1_eigenvalue = ifelse(length(eig) >= 1, eig[1], NA_real_),
    CAP2_eigenvalue = ifelse(length(eig) >= 2, eig[2], NA_real_),
    CAP1_percent_of_constrained = ifelse(length(eig_perc) >= 1, eig_perc[1], NA_real_),
    CAP2_percent_of_constrained = ifelse(length(eig_perc) >= 2, eig_perc[2], NA_real_)
  )
}

# ---------------------------------------------------------
# Themes and guides
# ---------------------------------------------------------

theme_dbrda_2panel <- theme_bw(base_size = 8.8) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15, colour = "grey90"),
    axis.text = element_text(size = 7.2),
    axis.title = element_text(size = 8.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.36, "cm"),
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, colour = "grey25"),
    plot.caption = element_text(hjust = 0, size = 6.8, colour = "grey25"),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

theme_dbrda_4panel <- theme_bw(base_size = 7.8) +
  theme(
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.25),
    strip.text.x = element_text(face = "bold", size = 7.8),
    strip.text.y = element_text(face = "bold", size = 7.8, angle = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.12, colour = "grey90"),
    axis.text = element_text(size = 6.1),
    axis.title = element_text(size = 7.2),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.5),
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.35, "cm"),
    plot.title = element_text(face = "bold", size = 10.5),
    plot.subtitle = element_text(size = 8.3, colour = "grey25"),
    plot.caption = element_text(hjust = 0, size = 6.8, colour = "grey25"),
    plot.margin = margin(5, 5, 5, 5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA)
  )

dbrda_guides <- guides(
  shape = guide_legend(
    order = 1,
    nrow = 1,
    override.aes = list(
      shape = c(21, 24),
      fill = "grey80",
      colour = "black",
      alpha = 1,
      size = 2.7,
      stroke = 0.6
    )
  ),
  fill = guide_legend(
    order = 2,
    nrow = 1,
    override.aes = list(
      shape = 21,
      colour = "black",
      alpha = 1,
      size = 2.5,
      stroke = 0.45
    )
  )
)

# ---------------------------------------------------------
# Containers
# ---------------------------------------------------------

all_results <- list()
model_summary_rows <- list()
test_summary_rows <- list()

dbrda_plot_tables <- list(
  site_scores = list(),
  centroids = list(),
  centroid_segments = list(),
  spider_segments = list(),
  overall_time_centroids = list(),
  overall_time_segment = list(),
  axis_labels = list()
)

# =========================================================
# Main loop
# =========================================================

for (analysis_name in names(analysis_setups)) {
  
  setup <- analysis_setups[[analysis_name]]
  all_results[[analysis_name]] <- list()
  
  for (dataset_name in names(datasets)) {
    
    dataset <- datasets[[dataset_name]]
    
    prepared <- prepare_dbrda_inputs(dataset$X, dataset$meta)
    X_work <- prepared$X
    meta_work <- prepared$meta
    
    print("=========================================================")
    print(paste(dataset$label, "-", setup$label))
    print("Sample counts")
    print(meta_work %>% count(forest, timespan))
    print("Community matrix dimensions after filtering")
    print(dim(X_work))
    
    # -----------------------------------------------------
    # Model
    # -----------------------------------------------------
    
    if (setup$formula_type == "residual") {
      
      cap_mod <- capscale(
        X_work ~ T_Sonian + T_Rajhenav + T_Strodam + T_Suserup +
          Condition(
            forest +
              log_av_ds +
              log_diameter +
              exposure +
              mosscover +
              snag
          ),
        data = meta_work,
        distance = "jaccard",
        binary = TRUE,
        add = "lingoes"
      )
      
    } else if (setup$formula_type == "reference") {
      
      cap_mod <- capscale(
        X_work ~ forest + timespan + forest:timespan +
          Condition(
            log_av_ds +
              log_diameter +
              exposure +
              mosscover +
              snag
          ),
        data = meta_work,
        distance = "jaccard",
        binary = TRUE,
        add = "lingoes"
      )
      
    } else {
      stop("Unknown formula_type in analysis_setups.")
    }
    
    print("capscale summary")
    print(cap_mod)
    
    # -----------------------------------------------------
    # Permutation tests
    # -----------------------------------------------------
    
    set.seed(42)
    
    if (setup$use_strata) {
      
      anova_overall <- anova.cca(
        cap_mod,
        permutations = n_perm,
        strata = meta_work$forest
      )
      
      anova_terms <- anova.cca(
        cap_mod,
        by = "terms",
        permutations = n_perm,
        strata = meta_work$forest
      )
      
      anova_margin <- anova.cca(
        cap_mod,
        by = "margin",
        permutations = n_perm,
        strata = meta_work$forest
      )
      
      anova_axes <- anova.cca(
        cap_mod,
        by = "axis",
        permutations = n_perm,
        strata = meta_work$forest
      )
      
    } else {
      
      anova_overall <- anova.cca(
        cap_mod,
        permutations = n_perm
      )
      
      anova_terms <- anova.cca(
        cap_mod,
        by = "terms",
        permutations = n_perm
      )
      
      anova_margin <- anova.cca(
        cap_mod,
        by = "margin",
        permutations = n_perm
      )
      
      anova_axes <- anova.cca(
        cap_mod,
        by = "axis",
        permutations = n_perm
      )
    }
    
    print("Overall test")
    print(anova_overall)
    
    print("Sequential term test")
    print(anova_terms)
    
    print("Marginal term test")
    print(anova_margin)
    
    print("Axis test")
    print(anova_axes)
    
    write_csv2(
      as.data.frame(anova_overall) %>% rownames_to_column("test"),
      file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_overall_test.csv"))
    )
    
    write_csv2(
      as.data.frame(anova_terms) %>% rownames_to_column("term"),
      file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_sequential_term_test.csv"))
    )
    
    write_csv2(
      as.data.frame(anova_margin) %>% rownames_to_column("term"),
      file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_marginal_term_test.csv"))
    )
    
    write_csv2(
      as.data.frame(anova_axes) %>% rownames_to_column("axis"),
      file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_axis_test.csv"))
    )
    
    # -----------------------------------------------------
    # Test summary tables
    # -----------------------------------------------------
    
    test_summary_rows[[paste(analysis_name, dataset_name, "overall", sep = "_")]] <-
      as.data.frame(anova_overall) %>%
      rownames_to_column("term") %>%
      mutate(
        analysis = setup$label,
        analysis_id = analysis_name,
        dataset = dataset$label,
        dataset_id = dataset_name,
        test_type = "overall",
        .before = 1
      )
    
    test_summary_rows[[paste(analysis_name, dataset_name, "terms", sep = "_")]] <-
      as.data.frame(anova_terms) %>%
      rownames_to_column("term") %>%
      mutate(
        analysis = setup$label,
        analysis_id = analysis_name,
        dataset = dataset$label,
        dataset_id = dataset_name,
        test_type = "sequential_terms",
        .before = 1
      )
    
    test_summary_rows[[paste(analysis_name, dataset_name, "margin", sep = "_")]] <-
      as.data.frame(anova_margin) %>%
      rownames_to_column("term") %>%
      mutate(
        analysis = setup$label,
        analysis_id = analysis_name,
        dataset = dataset$label,
        dataset_id = dataset_name,
        test_type = "marginal_terms",
        .before = 1
      )
    
    test_summary_rows[[paste(analysis_name, dataset_name, "axes", sep = "_")]] <-
      as.data.frame(anova_axes) %>%
      rownames_to_column("term") %>%
      mutate(
        analysis = setup$label,
        analysis_id = analysis_name,
        dataset = dataset$label,
        dataset_id = dataset_name,
        test_type = "axes",
        .before = 1
      )
    
    # -----------------------------------------------------
    # Site scores and centroids
    # -----------------------------------------------------
    
    site_scores <- scores(
      cap_mod,
      display = "sites",
      choices = 1:2,
      scaling = 1
    ) %>%
      as.data.frame() %>%
      rownames_to_column("sample_id") %>%
      rename(
        dbRDA1 = CAP1,
        dbRDA2 = CAP2
      ) %>%
      left_join(
        meta_work %>%
          as.data.frame() %>%
          select(-any_of("sample_id")) %>%
          rownames_to_column("sample_id"),
        by = "sample_id"
      ) %>%
      mutate(
        forest = factor(as.character(forest), levels = c("Sonian", "Rajhenav", "Strødam", "Suserup")),
        timespan = factor(as.character(timespan), levels = c("previous", "current"))
      )
    
    centroids <- site_scores %>%
      group_by(forest, timespan) %>%
      summarise(
        dbRDA1 = mean(dbRDA1, na.rm = TRUE),
        dbRDA2 = mean(dbRDA2, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        timespan_label = recode(as.character(timespan), !!!timespan_labs),
        centroid_label = recode(as.character(timespan), previous = "H", current = "C")
      )
    
    centroid_segments <- centroids %>%
      select(forest, timespan, dbRDA1, dbRDA2) %>%
      pivot_wider(
        names_from = timespan,
        values_from = c(dbRDA1, dbRDA2)
      ) %>%
      filter(
        !is.na(dbRDA1_previous),
        !is.na(dbRDA1_current),
        !is.na(dbRDA2_previous),
        !is.na(dbRDA2_current)
      )
    
    spider_segments <- site_scores %>%
      left_join(
        centroids %>%
          select(
            forest,
            timespan,
            centroid_dbRDA1 = dbRDA1,
            centroid_dbRDA2 = dbRDA2
          ),
        by = c("forest", "timespan")
      )
    
    # -----------------------------------------------------
    # Overall time centroid
    # Only interpreted as a visual aid for residual temporal space.
    # -----------------------------------------------------
    
    overall_time_centroids <- site_scores %>%
      group_by(timespan) %>%
      summarise(
        dbRDA1 = mean(dbRDA1, na.rm = TRUE),
        dbRDA2 = mean(dbRDA2, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        timespan_label = recode(as.character(timespan), !!!timespan_labs),
        centroid_label = recode(as.character(timespan), previous = "H", current = "C")
      )
    
    overall_time_segment <- overall_time_centroids %>%
      select(timespan, dbRDA1, dbRDA2) %>%
      pivot_wider(
        names_from = timespan,
        values_from = c(dbRDA1, dbRDA2)
      ) %>%
      filter(
        !is.na(dbRDA1_previous),
        !is.na(dbRDA1_current),
        !is.na(dbRDA2_previous),
        !is.na(dbRDA2_current)
      )
    
    # -----------------------------------------------------
    # Axis labels and model summary
    # -----------------------------------------------------
    
    axis_info <- get_axis_percentages(cap_mod)
    
    xlab <- paste0(
      setup$x_prefix,
      "1 (",
      round(axis_info$CAP1_percent_of_constrained[1], 1),
      "%)"
    )
    
    ylab <- paste0(
      setup$x_prefix,
      "2 (",
      round(axis_info$CAP2_percent_of_constrained[1], 1),
      "%)"
    )
    
    axis_label <- paste0(
      setup$x_prefix,
      "1: ",
      round(axis_info$CAP1_percent_of_constrained[1], 1),
      "%\n",
      setup$x_prefix,
      "2: ",
      round(axis_info$CAP2_percent_of_constrained[1], 1),
      "%"
    )
    
    anova_overall_df <- as.data.frame(anova_overall) %>%
      rownames_to_column("term")
    
    overall_model_row <- anova_overall_df %>%
      filter(term == "Model")
    
    if (nrow(overall_model_row) == 0) {
      overall_model_row <- anova_overall_df[1, , drop = FALSE]
    }
    
    total_inertia <- cap_mod$tot.chi
    conditional_inertia <- if (!is.null(cap_mod$pCCA)) cap_mod$pCCA$tot.chi else 0
    constrained_inertia <- if (!is.null(cap_mod$CCA)) cap_mod$CCA$tot.chi else 0
    unconstrained_inertia <- if (!is.null(cap_mod$CA)) cap_mod$CA$tot.chi else NA_real_
    
    conditional_rank <- if (!is.null(cap_mod$pCCA)) cap_mod$pCCA$rank else 0
    constrained_rank <- if (!is.null(cap_mod$CCA)) cap_mod$CCA$rank else 0
    unconstrained_rank <- if (!is.null(cap_mod$CA)) cap_mod$CA$rank else NA_integer_
    
    model_summary_rows[[paste(analysis_name, dataset_name, sep = "_")]] <- tibble(
      analysis = setup$label,
      analysis_id = analysis_name,
      dataset = dataset$label,
      dataset_id = dataset_name,
      n_samples = nrow(X_work),
      n_taxa = ncol(X_work),
      n_perm = n_perm,
      strata_used = setup$use_strata,
      total_inertia = total_inertia,
      conditional_inertia = conditional_inertia,
      constrained_inertia = constrained_inertia,
      unconstrained_inertia = unconstrained_inertia,
      conditional_proportion = conditional_inertia / total_inertia,
      constrained_proportion = constrained_inertia / total_inertia,
      unconstrained_proportion = unconstrained_inertia / total_inertia,
      conditional_rank = conditional_rank,
      constrained_rank = constrained_rank,
      unconstrained_rank = unconstrained_rank,
      CAP1_eigenvalue = axis_info$CAP1_eigenvalue[1],
      CAP2_eigenvalue = axis_info$CAP2_eigenvalue[1],
      CAP1_percent_of_constrained = axis_info$CAP1_percent_of_constrained[1],
      CAP2_percent_of_constrained = axis_info$CAP2_percent_of_constrained[1],
      overall_df = overall_model_row$Df[1],
      overall_sum_of_squares = overall_model_row$SumOfSqs[1],
      overall_F = overall_model_row$F[1],
      overall_p = overall_model_row$`Pr(>F)`[1]
    )
    
    # -----------------------------------------------------
    # Add facet variables and store long-format plot tables
    # -----------------------------------------------------
    
    site_scores_plot <- site_scores %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    centroids_plot <- centroids %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    centroid_segments_plot <- centroid_segments %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    spider_segments_plot <- spider_segments %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    overall_time_centroids_plot <- overall_time_centroids %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    overall_time_segment_plot <- overall_time_segment %>%
      mutate(
        analysis_id = analysis_name,
        analysis_scope = setup$label,
        dataset_id = dataset_name,
        data_scope = dataset$label
      )
    
    axis_labels_plot <- tibble(
      analysis_id = analysis_name,
      analysis_scope = setup$label,
      dataset_id = dataset_name,
      data_scope = dataset$label,
      x_axis_label = xlab,
      y_axis_label = ylab,
      axis_label = axis_label,
      axis_subtitle_label = paste0(
        dataset$label,
        ": dbRDA1 = ",
        round(axis_info$CAP1_percent_of_constrained[1], 1),
        "%, dbRDA2 = ",
        round(axis_info$CAP2_percent_of_constrained[1], 1),
        "%"
      ),
      axis_subtitle_label_short = paste0(
        dataset$label,
        " ",
        round(axis_info$CAP1_percent_of_constrained[1], 1),
        "% / ",
        round(axis_info$CAP2_percent_of_constrained[1], 1),
        "%"
      )
    )
    
    dbrda_plot_tables$site_scores[[paste(analysis_name, dataset_name, sep = "_")]] <- site_scores_plot
    dbrda_plot_tables$centroids[[paste(analysis_name, dataset_name, sep = "_")]] <- centroids_plot
    dbrda_plot_tables$centroid_segments[[paste(analysis_name, dataset_name, sep = "_")]] <- centroid_segments_plot
    dbrda_plot_tables$spider_segments[[paste(analysis_name, dataset_name, sep = "_")]] <- spider_segments_plot
    dbrda_plot_tables$overall_time_centroids[[paste(analysis_name, dataset_name, sep = "_")]] <- overall_time_centroids_plot
    dbrda_plot_tables$overall_time_segment[[paste(analysis_name, dataset_name, sep = "_")]] <- overall_time_segment_plot
    dbrda_plot_tables$axis_labels[[paste(analysis_name, dataset_name, sep = "_")]] <- axis_labels_plot
    
    # -----------------------------------------------------
    # Write plot tables
    # -----------------------------------------------------
    
    write_csv2(site_scores_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_site_scores.csv")))
    write_csv2(centroids_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_centroids.csv")))
    write_csv2(centroid_segments_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_centroid_segments.csv")))
    write_csv2(spider_segments_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_spider_segments.csv")))
    write_csv2(overall_time_centroids_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_overall_time_centroids.csv")))
    write_csv2(overall_time_segment_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_overall_time_segment.csv")))
    write_csv2(axis_labels_plot, file.path(out_dir, paste0(analysis_name, "_", dataset$file_label, "_", setup$table_file_suffix, "_axis_labels.csv")))
    
    all_results[[analysis_name]][[dataset_name]] <- list(
      model = cap_mod,
      anova_overall = anova_overall,
      anova_terms = anova_terms,
      anova_margin = anova_margin,
      anova_axes = anova_axes,
      site_scores = site_scores_plot,
      centroids = centroids_plot,
      centroid_segments = centroid_segments_plot,
      spider_segments = spider_segments_plot,
      overall_time_centroids = overall_time_centroids_plot,
      overall_time_segment = overall_time_segment_plot,
      axis_labels = axis_labels_plot
    )
  }
}

# =========================================================
# Bind long-format plotting tables
# =========================================================

dbrda_site_scores_all <- bind_rows(dbrda_plot_tables$site_scores)
dbrda_centroids_all <- bind_rows(dbrda_plot_tables$centroids)
dbrda_centroid_segments_all <- bind_rows(dbrda_plot_tables$centroid_segments)
dbrda_spider_segments_all <- bind_rows(dbrda_plot_tables$spider_segments)
dbrda_overall_time_centroids_all <- bind_rows(dbrda_plot_tables$overall_time_centroids)
dbrda_overall_time_segment_all <- bind_rows(dbrda_plot_tables$overall_time_segment)
dbrda_axis_labels_all <- bind_rows(dbrda_plot_tables$axis_labels)

write_csv2(dbrda_site_scores_all, file.path(out_dir, "all_dbrda_site_scores.csv"))
write_csv2(dbrda_centroids_all, file.path(out_dir, "all_dbrda_centroids.csv"))
write_csv2(dbrda_centroid_segments_all, file.path(out_dir, "all_dbrda_centroid_segments.csv"))
write_csv2(dbrda_spider_segments_all, file.path(out_dir, "all_dbrda_spider_segments.csv"))
write_csv2(dbrda_overall_time_centroids_all, file.path(out_dir, "all_dbrda_overall_time_centroids.csv"))
write_csv2(dbrda_overall_time_segment_all, file.path(out_dir, "all_dbrda_overall_time_segments.csv"))
write_csv2(dbrda_axis_labels_all, file.path(out_dir, "all_dbrda_axis_labels.csv"))

# =========================================================
# Helpers for subtitle axis percentages
# =========================================================

make_axis_subtitle <- function(axis_df) {
  
  axis_df %>%
    distinct(data_scope, axis_subtitle_label) %>%
    arrange(
      factor(
        data_scope,
        levels = c("Full dataset", "Substrate-matched dataset")
      )
    ) %>%
    pull(axis_subtitle_label) %>%
    paste(collapse = "; ")
}

make_axis_subtitle_4panel <- function(axis_df) {
  
  axis_df %>%
    distinct(analysis_scope, data_scope, axis_subtitle_label_short) %>%
    arrange(
      factor(
        analysis_scope,
        levels = c("Residual temporal space", "Substrate-adjusted reference space")
      ),
      factor(
        data_scope,
        levels = c("Full dataset", "Substrate-matched dataset")
      )
    ) %>%
    mutate(
      label = paste0(analysis_scope, ", ", axis_subtitle_label_short)
    ) %>%
    pull(label) %>%
    paste(collapse = "; ")
}

# =========================================================
# Faceted dbRDA plotting function
# Axis percentages are reported in subtitles, not in panels.
# =========================================================

plot_dbrda_faceted <- function(
    site_scores,
    spider_segments,
    centroids,
    centroid_segments,
    overall_time_centroids,
    overall_time_segment,
    axis_labels,
    facet_formula,
    plot_title,
    plot_subtitle,
    plot_caption = NULL,
    base_theme = theme_dbrda_2panel
) {
  
  ggplot() +
    geom_segment(
      data = spider_segments,
      aes(
        x = centroid_dbRDA1,
        y = centroid_dbRDA2,
        xend = dbRDA1,
        yend = dbRDA2,
        colour = forest
      ),
      linewidth = 0.20,
      alpha = 0.13,
      show.legend = FALSE
    ) +
    geom_point(
      data = site_scores,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        fill = forest,
        shape = timespan
      ),
      colour = "grey25",
      size = 1.20,
      alpha = 0.38,
      stroke = 0.18
    ) +
    geom_segment(
      data = centroid_segments,
      aes(
        x = dbRDA1_previous,
        y = dbRDA2_previous,
        xend = dbRDA1_current,
        yend = dbRDA2_current,
        colour = forest
      ),
      inherit.aes = FALSE,
      linewidth = 0.62,
      arrow = arrow(length = unit(0.10, "cm")),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        fill = forest,
        shape = timespan
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.35,
      stroke = 0.75,
      alpha = 0.92,
      show.legend = FALSE
    ) +
    geom_segment(
      data = overall_time_segment,
      aes(
        x = dbRDA1_previous,
        y = dbRDA2_previous,
        xend = dbRDA1_current,
        yend = dbRDA2_current
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.85,
      alpha = 0.95,
      arrow = arrow(length = unit(0.13, "cm")),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = overall_time_centroids,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        shape = timespan
      ),
      inherit.aes = FALSE,
      size = 2.55,
      fill = "grey85",
      colour = "black",
      stroke = 0.90,
      show.legend = FALSE
    ) +
    facet_grid(facet_formula, scales = "fixed") +
    scale_fill_manual(
      values = forest_cols,
      name = "Forest",
      drop = FALSE
    ) +
    scale_colour_manual(
      values = forest_cols,
      guide = "none",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = timespan_shapes,
      labels = timespan_labs,
      name = "Survey period",
      drop = FALSE
    ) +
    dbrda_guides +
    coord_equal() +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      caption = plot_caption,
      x = "dbRDA1",
      y = "dbRDA2"
    ) +
    base_theme
}

# =========================================================
# 2-panel plot: residual temporal space
# Full dataset versus substrate-matched dataset
# =========================================================

residual_id <- "residual_temporal_space"

p_dbrda_residual_2panel <- plot_dbrda_faceted(
  site_scores = dbrda_site_scores_all %>%
    filter(analysis_id == residual_id),
  spider_segments = dbrda_spider_segments_all %>%
    filter(analysis_id == residual_id),
  centroids = dbrda_centroids_all %>%
    filter(analysis_id == residual_id),
  centroid_segments = dbrda_centroid_segments_all %>%
    filter(analysis_id == residual_id),
  overall_time_centroids = dbrda_overall_time_centroids_all %>%
    filter(analysis_id == residual_id),
  overall_time_segment = dbrda_overall_time_segment_all %>%
    filter(analysis_id == residual_id),
  axis_labels = dbrda_axis_labels_all %>%
    filter(analysis_id == residual_id),
  facet_formula = . ~ data_scope,
  plot_title = "Residual temporal change after conditioning on forest and substrate variables",
  plot_subtitle = paste0(
    "Partial dbRDA spider plot for the full dataset and substrate-matched dataset\n",
    make_axis_subtitle(
      dbrda_axis_labels_all %>%
        filter(analysis_id == residual_id)
    )
  ),
  plot_caption = "Black arrows show the overall residual historical-to-contemporary shift. Coloured arrows show forest-specific residual centroid shifts.",
  base_theme = theme_dbrda_2panel
)

print(p_dbrda_residual_2panel)

ggsave(
  filename = file.path(out_dir, "faceted_dbrda_2panel_residual_full_vs_substrate_matched.png"),
  plot = p_dbrda_residual_2panel,
  width = 13,
  height = 6.5,
  dpi = 600
)

# =========================================================
# 2-panel plot: substrate-adjusted reference space
# Full dataset versus substrate-matched dataset
# =========================================================

reference_id <- "substrate_adjusted_reference_space"

p_dbrda_reference_2panel <- plot_dbrda_faceted(
  site_scores = dbrda_site_scores_all %>%
    filter(analysis_id == reference_id),
  spider_segments = dbrda_spider_segments_all %>%
    filter(analysis_id == reference_id),
  centroids = dbrda_centroids_all %>%
    filter(analysis_id == reference_id),
  centroid_segments = dbrda_centroid_segments_all %>%
    filter(analysis_id == reference_id),
  overall_time_centroids = dbrda_overall_time_centroids_all %>%
    filter(analysis_id == reference_id) %>%
    slice(0),
  overall_time_segment = dbrda_overall_time_segment_all %>%
    filter(analysis_id == reference_id) %>%
    slice(0),
  axis_labels = dbrda_axis_labels_all %>%
    filter(analysis_id == reference_id),
  facet_formula = . ~ data_scope,
  plot_title = "Substrate-adjusted reference-space dbRDA trajectories",
  plot_subtitle = paste0(
    "dbRDA spider plot retaining forest structure after conditioning on substrate variables\n",
    make_axis_subtitle(
      dbrda_axis_labels_all %>%
        filter(analysis_id == reference_id)
    )
  ),
  plot_caption = "Coloured arrows show forest-specific historical-to-contemporary centroid shifts.",
  base_theme = theme_dbrda_2panel
)

print(p_dbrda_reference_2panel)

ggsave(
  filename = file.path(out_dir, "faceted_dbrda_2panel_reference_full_vs_substrate_matched.png"),
  plot = p_dbrda_reference_2panel,
  width = 13,
  height = 6.5,
  dpi = 600
)

# =========================================================
# 4-panel plot
# Rows = dbRDA model space
# Columns = full dataset versus substrate-matched dataset
# =========================================================

p_dbrda_4panel <- plot_dbrda_faceted(
  site_scores = dbrda_site_scores_all,
  spider_segments = dbrda_spider_segments_all,
  centroids = dbrda_centroids_all,
  centroid_segments = dbrda_centroid_segments_all,
  overall_time_centroids = dbrda_overall_time_centroids_all %>%
    filter(analysis_id == residual_id),
  overall_time_segment = dbrda_overall_time_segment_all %>%
    filter(analysis_id == residual_id),
  axis_labels = dbrda_axis_labels_all,
  facet_formula = analysis_scope ~ data_scope,
  plot_title = "Partial dbRDA spider plots",
  plot_subtitle = paste0(
    "Rows show the two constrained ordination spaces; columns compare the full dataset and substrate-matched dataset\n",
    make_axis_subtitle_4panel(dbrda_axis_labels_all)
  ),
  plot_caption = "The residual temporal space conditions on forest identity and substrate variables. The reference space retains forest structure while conditioning on substrate variables. Black arrows are shown only for the residual temporal-space model.",
  base_theme = theme_dbrda_4panel
)

print(p_dbrda_4panel)

ggsave(
  filename = file.path(out_dir, "faceted_dbrda_4panel_full_vs_substrate_matched.png"),
  plot = p_dbrda_4panel,
  width = 13,
  height = 11.5,
  dpi = 600
)

# =========================================================
# Individual faceted 1-panel plots
# Useful for checking each model × dataset combination
# =========================================================

for (analysis_name in names(analysis_setups)) {
  
  for (dataset_name in names(datasets)) {
    
    analysis_label <- analysis_setups[[analysis_name]]$label
    dataset_label <- datasets[[dataset_name]]$label
    dataset_file_label <- datasets[[dataset_name]]$file_label
    setup <- analysis_setups[[analysis_name]]
    
    overall_centroids_plot <- dbrda_overall_time_centroids_all %>%
      filter(analysis_id == analysis_name, dataset_id == dataset_name)
    
    overall_segment_plot <- dbrda_overall_time_segment_all %>%
      filter(analysis_id == analysis_name, dataset_id == dataset_name)
    
    if (setup$formula_type != "residual") {
      overall_centroids_plot <- overall_centroids_plot %>% slice(0)
      overall_segment_plot <- overall_segment_plot %>% slice(0)
    }
    
    p_individual <- plot_dbrda_faceted(
      site_scores = dbrda_site_scores_all %>%
        filter(analysis_id == analysis_name, dataset_id == dataset_name),
      spider_segments = dbrda_spider_segments_all %>%
        filter(analysis_id == analysis_name, dataset_id == dataset_name),
      centroids = dbrda_centroids_all %>%
        filter(analysis_id == analysis_name, dataset_id == dataset_name),
      centroid_segments = dbrda_centroid_segments_all %>%
        filter(analysis_id == analysis_name, dataset_id == dataset_name),
      overall_time_centroids = overall_centroids_plot,
      overall_time_segment = overall_segment_plot,
      axis_labels = dbrda_axis_labels_all %>%
        filter(analysis_id == analysis_name, dataset_id == dataset_name),
      facet_formula = . ~ data_scope,
      plot_title = dataset_label,
      plot_subtitle = paste0(
        setup$subtitle,
        "\n",
        make_axis_subtitle(
          dbrda_axis_labels_all %>%
            filter(analysis_id == analysis_name, dataset_id == dataset_name)
        )
      ),
      plot_caption = NULL,
      base_theme = theme_dbrda_2panel
    )
    
    print(p_individual)
    
    ggsave(
      filename = file.path(
        out_dir,
        paste0(
          analysis_name,
          "_",
          dataset_file_label,
          "_",
          setup$plot_file_suffix,
          ".png"
        )
      ),
      plot = p_individual,
      width = 8.2,
      height = 6.2,
      dpi = 900
    )
  }
}

# =========================================================
# General summary tables
# =========================================================

general_model_summary <- bind_rows(model_summary_rows) %>%
  mutate(
    conditional_percent = 100 * conditional_proportion,
    constrained_percent = 100 * constrained_proportion,
    unconstrained_percent = 100 * unconstrained_proportion
  ) %>%
  select(
    analysis,
    analysis_id,
    dataset,
    dataset_id,
    n_samples,
    n_taxa,
    n_perm,
    strata_used,
    total_inertia,
    conditional_inertia,
    constrained_inertia,
    unconstrained_inertia,
    conditional_percent,
    constrained_percent,
    unconstrained_percent,
    conditional_rank,
    constrained_rank,
    unconstrained_rank,
    CAP1_eigenvalue,
    CAP2_eigenvalue,
    CAP1_percent_of_constrained,
    CAP2_percent_of_constrained,
    overall_df,
    overall_sum_of_squares,
    overall_F,
    overall_p
  )

general_test_summary <- bind_rows(test_summary_rows) %>%
  select(
    analysis,
    analysis_id,
    dataset,
    dataset_id,
    test_type,
    term,
    everything()
  )

print("=========================================================")
print("GENERAL dbRDA MODEL SUMMARY")
print(general_model_summary)

print("=========================================================")
print("GENERAL dbRDA TEST SUMMARY")
print(general_test_summary)

write_csv2(general_model_summary, file.path(out_dir, "general_dbrda_model_summary.csv"))
write_csv2(general_test_summary, file.path(out_dir, "general_dbrda_test_summary.csv"))

# =========================================================
# BOLT-ON: SONIAN DISTANCE-TO-REFERENCE CHANGE TABLE
#
# Uses the substrate-adjusted reference-space dbRDA scores
# already stored in all_results.
#
# Output:
#   1. Detailed Sonian-reference distance table with bootstrap
#   2. Simplified numeric manuscript table
#   3. Wide paper-style table with delta distance and percentage
#
# Interpretation:
#   Delta = contemporary distance - historical distance.
#   Negative delta means Sonian became closer to the reference forest.
# =========================================================

reference_forests <- c("Strødam", "Suserup", "Rajhenav")

calculate_sonian_reference_distance <- function(site_scores, dataset_label, n_boot = 999) {
  
  site_scores <- site_scores %>%
    mutate(
      forest = factor(as.character(forest), levels = c("Sonian", "Rajhenav", "Strødam", "Suserup")),
      timespan = factor(as.character(timespan), levels = c("previous", "current"))
    )
  
  calculate_distances_once <- function(dat) {
    
    centroids_local <- dat %>%
      group_by(forest, timespan) %>%
      summarise(
        dbRDA1 = mean(dbRDA1, na.rm = TRUE),
        dbRDA2 = mean(dbRDA2, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    
    get_centroid <- function(forest_name, period_name) {
      
      out <- centroids_local %>%
        filter(forest == forest_name, timespan == period_name) %>%
        select(dbRDA1, dbRDA2)
      
      if (nrow(out) != 1) {
        stop(paste("Missing centroid for", forest_name, period_name), call. = FALSE)
      }
      
      as.numeric(out[1, ])
    }
    
    sonian_previous <- get_centroid("Sonian", "previous")
    sonian_current <- get_centroid("Sonian", "current")
    
    bind_rows(lapply(reference_forests, function(ref) {
      
      ref_previous <- get_centroid(ref, "previous")
      ref_current <- get_centroid(ref, "current")
      
      historical_distance <- sqrt(sum((sonian_previous - ref_previous)^2))
      contemporary_distance <- sqrt(sum((sonian_current - ref_current)^2))
      distance_delta <- contemporary_distance - historical_distance
      
      tibble(
        reference_forest = ref,
        historical_distance = historical_distance,
        contemporary_distance = contemporary_distance,
        distance_delta_current_minus_historical = distance_delta,
        percent_change = 100 * distance_delta / historical_distance,
        sonian_became_closer = distance_delta < 0
      )
    }))
  }
  
  bootstrap_once <- function(dat) {
    
    dat %>%
      group_by(forest, timespan) %>%
      group_modify(~ .x[sample(seq_len(nrow(.x)), size = nrow(.x), replace = TRUE), , drop = FALSE]) %>%
      ungroup()
  }
  
  set.seed(42)
  
  observed <- calculate_distances_once(site_scores)
  
  boot_res <- bind_rows(lapply(seq_len(n_boot), function(i) {
    calculate_distances_once(bootstrap_once(site_scores)) %>%
      mutate(bootstrap_iteration = i)
  }))
  
  boot_summary <- boot_res %>%
    group_by(reference_forest) %>%
    summarise(
      boot_delta_mean = mean(distance_delta_current_minus_historical, na.rm = TRUE),
      boot_delta_ci_lower_95 = quantile(distance_delta_current_minus_historical, 0.025, na.rm = TRUE),
      boot_delta_ci_upper_95 = quantile(distance_delta_current_minus_historical, 0.975, na.rm = TRUE),
      p_boot_closer = mean(distance_delta_current_minus_historical < 0, na.rm = TRUE),
      p_boot_more_distant = mean(distance_delta_current_minus_historical > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  observed %>%
    left_join(boot_summary, by = "reference_forest") %>%
    mutate(
      dataset = dataset_label,
      interpretation = case_when(
        distance_delta_current_minus_historical < 0 ~ "Sonian closer to reference forest in contemporary survey",
        distance_delta_current_minus_historical > 0 ~ "Sonian more distant from reference forest in contemporary survey",
        TRUE ~ "No change in centroid distance"
      ),
      .before = 1
    )
}

sonian_reference_distance_full <- calculate_sonian_reference_distance(
  site_scores = all_results$substrate_adjusted_reference_space$full$site_scores,
  dataset_label = "Full dataset",
  n_boot = n_boot_direction
)

sonian_reference_distance_subset <- calculate_sonian_reference_distance(
  site_scores = all_results$substrate_adjusted_reference_space$substrate_matched$site_scores,
  dataset_label = "Substrate-matched dataset",
  n_boot = n_boot_direction
)

sonian_reference_distance_summary <- bind_rows(
  sonian_reference_distance_full,
  sonian_reference_distance_subset
) %>%
  mutate(
    dataset = factor(dataset, levels = c("Full dataset", "Substrate-matched dataset")),
    reference_forest = factor(reference_forest, levels = c("Strødam", "Suserup", "Rajhenav"))
  ) %>%
  arrange(dataset, reference_forest)

print("=========================================================")
print("SONIAN DISTANCE-TO-REFERENCE SUMMARY WITH BOOTSTRAP")
print("Delta = contemporary distance - historical distance. Negative delta means Sonian became closer.")
print(sonian_reference_distance_summary, n = Inf)

write_csv2(
  sonian_reference_distance_summary,
  file.path(out_dir, "sonian_reference_distance_summary_with_bootstrap.csv")
)

sonian_reference_simple_numeric <- sonian_reference_distance_summary %>%
  transmute(
    dataset = as.character(dataset),
    reference_forest = as.character(reference_forest),
    historical_distance = round(historical_distance, 3),
    contemporary_distance = round(contemporary_distance, 3),
    delta_distance = round(distance_delta_current_minus_historical, 3),
    percent_change = round(percent_change, 1),
    bootstrap_support_closer = round(p_boot_closer, 3)
  ) %>%
  arrange(
    factor(dataset, levels = c("Full dataset", "Substrate-matched dataset")),
    factor(reference_forest, levels = c("Strødam", "Suserup", "Rajhenav"))
  )

print("=========================================================")
print("SIMPLIFIED SONIAN DISTANCE-TO-REFERENCE TABLE: NUMERIC")
print("Delta = contemporary distance - historical distance. Negative delta means Sonian became closer.")
print(sonian_reference_simple_numeric, n = Inf)

write_csv2(
  sonian_reference_simple_numeric,
  file.path(out_dir, "sonian_reference_simple_numeric_table.csv")
)

sonian_reference_delta_percent_wide <- sonian_reference_simple_numeric %>%
  mutate(
    dataset = factor(dataset, levels = c("Full dataset", "Substrate-matched dataset")),
    reference_forest = factor(reference_forest, levels = c("Strødam", "Suserup", "Rajhenav")),
    table_value = paste0(
      sprintf("%.3f", delta_distance),
      " (",
      sprintf("%.1f", percent_change),
      "%)"
    )
  ) %>%
  select(
    dataset,
    reference_forest,
    table_value
  ) %>%
  pivot_wider(
    names_from = reference_forest,
    values_from = table_value
  ) %>%
  arrange(dataset) %>%
  mutate(dataset = as.character(dataset))

print("=========================================================")
print("PAPER-STYLE TABLE: DELTA DISTANCE AND PERCENTAGE CHANGE")
print("Values are delta distance with percentage change in parentheses.")
print("Delta = contemporary distance - historical distance.")
print("Negative values indicate reduced distance from Sonian to the reference forest.")
print(sonian_reference_delta_percent_wide, n = Inf)

write_csv2(
  sonian_reference_delta_percent_wide,
  file.path(out_dir, "sonian_reference_delta_percent_wide_table_for_manuscript.csv")
)

# =========================================================
# Notes for interpretation
# =========================================================
#
# 1. The residual temporal-space model conditions on forest identity
#    and measured substrate variables. The black arrow is therefore
#    a visual aid for the average residual historical-to-contemporary
#    displacement after removing those gradients.
#
# 2. The substrate-adjusted reference-space model conditions only on
#    measured substrate variables. It retains forest structure and is
#    therefore more appropriate for visualizing whether Sonian moves
#    toward longer-continuity reference forests in the constrained space.
#
# 3. Axis percentages differ among models and datasets because each
#    dbRDA is fitted separately. They are reported in the subtitle
#    rather than inside the ordination panels.
#
# 4. The 4-panel figure is organized like the NMDS figure:
#    rows = ordination space
#    columns = full dataset versus substrate-matched dataset
#
# 5. The Sonian-to-reference distance tables are based on the
#    substrate-adjusted reference-space dbRDA scores. Negative delta
#    values indicate that Sonian became closer to the reference forest
#    centroid in the contemporary survey.
#
# =========================================================
# =========================================================
# Publication-ready dbRDA model summary table
#
# One row per dbRDA model × dataset.
# Axis percentages are percentages of constrained variation,
# not percentages of total community variation.
# =========================================================

dbrda_publication_table <- general_model_summary %>%
  mutate(
    model_description = case_when(
      analysis_id == "residual_temporal_space" ~
        "Residual temporal space",
      analysis_id == "substrate_adjusted_reference_space" ~
        "Substrate-adjusted reference space",
      TRUE ~ analysis
    ),
    conditioned_terms = case_when(
      analysis_id == "residual_temporal_space" ~
        "Forest identity, decay stage, log diameter, exposure, moss cover, snag presence",
      analysis_id == "substrate_adjusted_reference_space" ~
        "Decay stage, log diameter, exposure, moss cover, snag presence",
      TRUE ~ NA_character_
    ),
    constrained_terms = case_when(
      analysis_id == "residual_temporal_space" ~
        "Forest-specific contemporary indicators",
      analysis_id == "substrate_adjusted_reference_space" ~
        "Forest identity, survey period, forest × survey period",
      TRUE ~ NA_character_
    ),
    total_inertia = round(total_inertia, 4),
    conditioned_percent = round(conditional_percent, 1),
    constrained_percent = round(constrained_percent, 1),
    unconstrained_percent = round(unconstrained_percent, 1),
    dbRDA1_percent_constrained = round(CAP1_percent_of_constrained, 1),
    dbRDA2_percent_constrained = round(CAP2_percent_of_constrained, 1),
    overall_F = round(overall_F, 3),
    overall_p = case_when(
      is.na(overall_p) ~ NA_character_,
      overall_p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", overall_p)
    )
  ) %>%
  transmute(
    dataset,
    model = model_description,
    conditioned_terms,
    constrained_terms,
    n_logs = n_samples,
    n_taxa,
    total_inertia,
    conditioned_percent,
    constrained_percent,
    unconstrained_percent,
    dbRDA1_percent_constrained,
    dbRDA2_percent_constrained,
    overall_F,
    overall_p
  ) %>%
  arrange(
    factor(model, levels = c(
      "Substrate-adjusted reference space",
      "Residual temporal space"
    )),
    factor(dataset, levels = c(
      "Full dataset",
      "Substrate-matched dataset"
    ))
  )

print("=========================================================")
print("PUBLICATION-READY dbRDA MODEL SUMMARY TABLE")
print("Percent columns for conditioned, constrained, and unconstrained components refer to total inertia.")
print("dbRDA1 and dbRDA2 percentages refer to constrained variation only.")
print(dbrda_publication_table, n = Inf)

write_csv2(
  dbrda_publication_table,
  file.path(out_dir, "publication_ready_dbrda_model_summary_table.csv")
)
# =========================================================
# Print each of the four dbRDA panels separately
# With wrap-style title, no legend, larger labels
# Shared axis limits within each dbRDA formula
# =========================================================

make_shared_dbrda_limits <- function(
    site_scores,
    centroids,
    centroid_segments,
    overall_time_centroids,
    overall_time_segment,
    expansion = 0.08
) {
  
  site_xy <- site_scores %>%
    select(analysis_id, dbRDA1, dbRDA2)
  
  centroid_xy <- centroids %>%
    select(analysis_id, dbRDA1, dbRDA2)
  
  centroid_segment_xy <- bind_rows(
    centroid_segments %>%
      transmute(
        analysis_id,
        dbRDA1 = dbRDA1_previous,
        dbRDA2 = dbRDA2_previous
      ),
    centroid_segments %>%
      transmute(
        analysis_id,
        dbRDA1 = dbRDA1_current,
        dbRDA2 = dbRDA2_current
      )
  )
  
  overall_centroid_xy <- overall_time_centroids %>%
    select(analysis_id, dbRDA1, dbRDA2)
  
  overall_segment_xy <- bind_rows(
    overall_time_segment %>%
      transmute(
        analysis_id,
        dbRDA1 = dbRDA1_previous,
        dbRDA2 = dbRDA2_previous
      ),
    overall_time_segment %>%
      transmute(
        analysis_id,
        dbRDA1 = dbRDA1_current,
        dbRDA2 = dbRDA2_current
      )
  )
  
  bind_rows(
    site_xy,
    centroid_xy,
    centroid_segment_xy,
    overall_centroid_xy,
    overall_segment_xy
  ) %>%
    filter(
      is.finite(dbRDA1),
      is.finite(dbRDA2)
    ) %>%
    group_by(analysis_id) %>%
    summarise(
      x_min = min(dbRDA1, na.rm = TRUE),
      x_max = max(dbRDA1, na.rm = TRUE),
      y_min = min(dbRDA2, na.rm = TRUE),
      y_max = max(dbRDA2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      x_range = if_else(x_max > x_min, x_max - x_min, 1),
      y_range = if_else(y_max > y_min, y_max - y_min, 1),
      x_min = x_min - expansion * x_range,
      x_max = x_max + expansion * x_range,
      y_min = y_min - expansion * y_range,
      y_max = y_max + expansion * y_range
    )
}

shared_dbrda_limits <- make_shared_dbrda_limits(
  site_scores = dbrda_site_scores_all,
  centroids = dbrda_centroids_all,
  centroid_segments = dbrda_centroid_segments_all,
  overall_time_centroids = dbrda_overall_time_centroids_all,
  overall_time_segment = dbrda_overall_time_segment_all,
  expansion = 0.08
)

plot_dbrda_single_clean <- function(
    site_scores,
    spider_segments,
    centroids,
    centroid_segments,
    overall_time_centroids,
    overall_time_segment,
    plot_title = NULL,
    plot_caption = NULL,
    x_limits = NULL,
    y_limits = NULL,
    base_theme = theme_dbrda_2panel
) {
  
  ggplot() +
    geom_segment(
      data = spider_segments,
      aes(
        x = centroid_dbRDA1,
        y = centroid_dbRDA2,
        xend = dbRDA1,
        yend = dbRDA2,
        colour = forest
      ),
      linewidth = 0.28,
      alpha = 0.14,
      show.legend = FALSE
    ) +
    geom_point(
      data = site_scores,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        fill = forest,
        shape = timespan
      ),
      colour = "grey25",
      size = 2.05,
      alpha = 0.42,
      stroke = 0.22,
      show.legend = FALSE
    ) +
    geom_segment(
      data = centroid_segments,
      aes(
        x = dbRDA1_previous,
        y = dbRDA2_previous,
        xend = dbRDA1_current,
        yend = dbRDA2_current,
        colour = forest
      ),
      inherit.aes = FALSE,
      linewidth = 0.90,
      arrow = arrow(length = unit(0.16, "cm")),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        fill = forest,
        shape = timespan
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 3,
      stroke = 0.95,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_segment(
      data = overall_time_segment,
      aes(
        x = dbRDA1_previous,
        y = dbRDA2_previous,
        xend = dbRDA1_current,
        yend = dbRDA2_current
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1.05,
      alpha = 0.95,
      arrow = arrow(length = unit(0.18, "cm")),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_point(
      data = overall_time_centroids,
      aes(
        x = dbRDA1,
        y = dbRDA2,
        shape = timespan
      ),
      inherit.aes = FALSE,
      size = 3.35,
      fill = "grey85",
      colour = "black",
      stroke = 1.05,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = forest_cols,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = forest_cols,
      guide = "none",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = timespan_shapes,
      labels = timespan_labs,
      drop = FALSE
    ) +
    coord_equal(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE
    ) +
    labs(
      title = plot_title,
      subtitle = NULL,
      caption = plot_caption,
      x = "dbRDA1",
      y = "dbRDA2"
    ) +
    base_theme +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        size = 17,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_blank(),
      axis.title = element_text(
        size = 15,
        face = "bold"
      ),
      axis.text = element_text(
        size = 13,
        colour = "black"
      ),
      plot.caption = element_text(
        hjust = 0,
        size = 10,
        colour = "grey25"
      ),
      strip.text = element_text(
        size = 14,
        face = "bold"
      )
    )
}

clean_single_plot_specs <- tibble::tribble(
  ~analysis_id, ~dataset_id, ~plot_title, ~file_stub,
  "residual_temporal_space", "full",
  "Residual temporal space: full dataset",
  "clean_single_residual_temporal_space_full_dataset",
  
  "residual_temporal_space", "substrate_matched",
  "Residual temporal space: substrate-matched dataset",
  "clean_single_residual_temporal_space_substrate_matched_dataset",
  
  "substrate_adjusted_reference_space", "full",
  "Substrate-adjusted reference space: full dataset",
  "clean_single_reference_space_full_dataset",
  
  "substrate_adjusted_reference_space", "substrate_matched",
  "Substrate-adjusted reference space: substrate-matched dataset",
  "clean_single_reference_space_substrate_matched_dataset"
)

for (i in seq_len(nrow(clean_single_plot_specs))) {
  
  this_analysis_id <- clean_single_plot_specs$analysis_id[i]
  this_dataset_id <- clean_single_plot_specs$dataset_id[i]
  this_plot_title <- clean_single_plot_specs$plot_title[i]
  this_file_stub <- clean_single_plot_specs$file_stub[i]
  
  this_limits <- shared_dbrda_limits %>%
    filter(analysis_id == this_analysis_id)
  
  this_x_limits <- c(this_limits$x_min[1], this_limits$x_max[1])
  this_y_limits <- c(this_limits$y_min[1], this_limits$y_max[1])
  
  clean_overall_time_centroids <- dbrda_overall_time_centroids_all %>%
    filter(
      analysis_id == this_analysis_id,
      dataset_id == this_dataset_id
    )
  
  clean_overall_time_segment <- dbrda_overall_time_segment_all %>%
    filter(
      analysis_id == this_analysis_id,
      dataset_id == this_dataset_id
    )
  
  if (this_analysis_id != "residual_temporal_space") {
    clean_overall_time_centroids <- clean_overall_time_centroids %>%
      slice(0)
    
    clean_overall_time_segment <- clean_overall_time_segment %>%
      slice(0)
  }
  
  p_clean_single <- plot_dbrda_single_clean(
    site_scores = dbrda_site_scores_all %>%
      filter(
        analysis_id == this_analysis_id,
        dataset_id == this_dataset_id
      ),
    spider_segments = dbrda_spider_segments_all %>%
      filter(
        analysis_id == this_analysis_id,
        dataset_id == this_dataset_id
      ),
    centroids = dbrda_centroids_all %>%
      filter(
        analysis_id == this_analysis_id,
        dataset_id == this_dataset_id
      ),
    centroid_segments = dbrda_centroid_segments_all %>%
      filter(
        analysis_id == this_analysis_id,
        dataset_id == this_dataset_id
      ),
    overall_time_centroids = clean_overall_time_centroids,
    overall_time_segment = clean_overall_time_segment,
    plot_title = this_plot_title,
    plot_caption = NULL,
    x_limits = this_x_limits,
    y_limits = this_y_limits,
    base_theme = theme_dbrda_2panel
  )
  
  print("=========================================================")
  print(paste("CLEAN SINGLE dbRDA PLOT:", this_file_stub))
  print("No legend; wrap-style title restored; enlarged labels; shared limits within dbRDA formula")
  print(p_clean_single)
  
  ggsave(
    filename = file.path(out_dir, paste0(this_file_stub, ".png")),
    plot = p_clean_single,
    width = 6.8 / 1.3,
    height = 6.0 / 1.3,
    dpi = 900
  )
}
