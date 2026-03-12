# =============================================================================
# CROSS + LOGS null Amplification strategy check
# =============================================================================
setwd(".")

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(lme4)
})

stopifnot(exists("META1"), exists("otu_matrix_filt"))

# -----------------------------------------------------------------------------
# 0) Directories and global settings
# -----------------------------------------------------------------------------
null_plot_dir     <- file.path("plots", "ALL", "nulls")
null_plot_rds_dir <- file.path(null_plot_dir, "rds")
null_table_dir    <- file.path("tables", "ALL", "nulls")

dir.create(null_plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(null_plot_rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(null_table_dir, recursive = TRUE, showWarnings = FALSE)

SEED_GLOBAL <- get0("SEED_GLOBAL", ifnotfound = 42L)
set.seed(SEED_GLOBAL)

B_NULL      <- get0("B_NULL", ifnotfound = 999L)
N_PERM_NULL <- get0("N_PERM_NULL", ifnotfound = 199L)
PERM_ANOSIM <- get0("PERM_ANOSIM", ifnotfound = 9999L)

# -----------------------------------------------------------------------------
# 1) Generic helper functions
# -----------------------------------------------------------------------------
boot_diff_median <- function(x_within, x_between, nboot = 1999, seed = 42) {
  set.seed(seed)
  
  x_within  <- x_within[is.finite(x_within)]
  x_between <- x_between[is.finite(x_between)]
  
  nW <- length(x_within)
  nB <- length(x_between)
  
  if (nW < 1 || nB < 1) return(rep(NA_real_, nboot))
  
  replicate(nboot, {
    med_within  <- median(sample(x_within,  size = nW, replace = TRUE))
    med_between <- median(sample(x_between, size = nB, replace = TRUE))
    med_between - med_within
  })
}

mk_ecdf_plot <- function(D, grp, title, xlab) {
  U  <- upper_idx_from(D)
  pw <- pairwise_groups(D, grp, U)
  
  df <- tibble(
    dist = c(pw$within, pw$between),
    type = factor(
      rep(c("within", "between"),
          c(length(pw$within), length(pw$between))),
      levels = c("within", "between")
    )
  )
  
  sumtab <- df %>%
    group_by(type) %>%
    summarise(n = n(), med = median(dist), .groups = "drop")
  
  bd <- boot_diff_median(
    df$dist[df$type == "within"],
    df$dist[df$type == "between"],
    nboot = 1999
  )
  ci  <- stats::quantile(bd, c(0.025, 0.975), na.rm = TRUE)
  est <- unname(sumtab$med[sumtab$type == "between"] -
                  sumtab$med[sumtab$type == "within"])
  
  xpos <- stats::quantile(df$dist, 0.70, na.rm = TRUE)
  ypos <- 0.20
  lab  <- sprintf("Δ median = %.3f\n95%% CI [%.3f, %.3f]", est, ci[1], ci[2])
  
  ggplot(df, aes(x = dist, colour = type)) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    geom_vline(
      xintercept = sumtab$med[sumtab$type == "within"],
      linetype = 2, alpha = 0.6
    ) +
    geom_vline(
      xintercept = sumtab$med[sumtab$type == "between"],
      linetype = 2, alpha = 0.6
    ) +
    annotate("text", x = xpos, y = ypos, hjust = 0, label = lab, size = 3) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_blank()
    ) +
    labs(title = title, x = xlab, y = "ECDF", colour = NULL)
}

mk_density_plot <- function(null_vals, obs, title, xlab) {
  ci  <- stats::quantile(null_vals, c(0.025, 0.975), na.rm = TRUE)
  ses <- (obs - mean(null_vals, na.rm = TRUE)) / stats::sd(null_vals, na.rm = TRUE)
  p   <- mean(null_vals >= obs, na.rm = TRUE)
  
  ggplot(tibble(x = null_vals), aes(x)) +
    geom_density(fill = "grey85") +
    geom_vline(xintercept = obs, linetype = 2) +
    geom_vline(xintercept = ci, linetype = 3) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(
      title = title,
      subtitle = sprintf(
        "obs=%.3f | null μ=%.3f | 95%%CI[%.3f, %.3f] | p=%.3f | SES=%.2f",
        obs, mean(null_vals, na.rm = TRUE), ci[1], ci[2], p, ses
      ),
      x = xlab,
      y = "Density"
    )
}

make_summary_wide_csv2 <- function(final_long_tbl,
                                   scope_label = "CROSS",
                                   keep_metric_order = c("robust.aitchison", "jaccard"),
                                   keep_dataset_order = c("UMI", "nonUMI")) {
  stopifnot(is.data.frame(final_long_tbl))
  
  req_cols <- c(
    "Scope", "Metric", "Dataset", "Term",
    "R2 obs", "R2 null mean", "R2 null 95% CI",
    "p (R2)", "R2 SES", "R2 corrected",
    "AUC obs", "AUC null mean", "AUC null 95% CI",
    "p (AUC)", "AUC SES",
    "ANOSIM R", "ANOSIM p",
    "n trees", "n samples"
  )
  
  missing <- setdiff(req_cols, names(final_long_tbl))
  if (length(missing) > 0) {
    stop("final_long_tbl missing required columns: ", paste(missing, collapse = ", "))
  }
  
  out <- final_long_tbl %>%
    mutate(
      Scope = scope_label,
      Library = Dataset,
      `design R²` = case_when(
        Term == "natman"   ~ "natman",
        Term == "dw_type2" ~ "dw_type2",
        TRUE ~ as.character(Term)
      ),
      `design AUC` = case_when(
        Term == "natman"   ~ "natman",
        Term == "dw_type2" ~ "dw_type2",
        TRUE ~ as.character(Term)
      )
    ) %>%
    transmute(
      Scope,
      Library,
      Metric,
      `design R²`,
      `R² obs`                    = `R2 obs`,
      `R² null mean`              = `R2 null mean`,
      `R² null 95% CI`            = `R2 null 95% CI`,
      `p (R²)`                    = `p (R2)`,
      `R² SES`                    = `R2 SES`,
      `R² corrected (obs-µ null)` = `R2 corrected`,
      `design AUC`,
      `AUC obs`,
      `AUC null mean`,
      `AUC null 95% CI`,
      `p (AUC)`,
      `AUC SES`,
      `ANOSIM R`,
      `ANOSIM p`,
      `n trees`,
      `n samples`
    )
  
  if (all(keep_dataset_order %in% out$Library)) {
    out <- out %>% mutate(Library = factor(Library, levels = keep_dataset_order))
  }
  if (all(keep_metric_order %in% out$Metric)) {
    out <- out %>% mutate(Metric = factor(Metric, levels = keep_metric_order))
  }
  
  out %>%
    arrange(Library, Metric, `design R²`) %>%
    mutate(
      Library = as.character(Library),
      Metric  = as.character(Metric)
    )
}

# -----------------------------------------------------------------------------
# 2) Generic null runner for one grouping term
# -----------------------------------------------------------------------------
run_metric_single_term <- function(D_raw, meta_raw, metric_label, xlab_short,
                                   dataset_tag = c("UMI", "nonUMI"),
                                   scope_label,
                                   term_var = c("natman", "dw_type2"),
                                   add_log_reads = TRUE,
                                   min_per_tree = NULL,
                                   save_prefix = NULL) {
  dataset_tag <- match.arg(dataset_tag)
  term_var    <- match.arg(term_var)
  
  stopifnot("sample" %in% names(meta_raw))
  stopifnot("umi_backup" %in% names(meta_raw))
  
  keep_flag <- if (dataset_tag == "UMI") 1 else 0
  
  meta_sub <- meta_raw %>%
    filter(umi_backup == keep_flag) %>%
    droplevels()
  
  if (!is.null(min_per_tree) && term_var == "natman") {
    meta_sub <- meta_sub %>%
      add_count(natman, name = "n_tree") %>%
      filter(n_tree >= min_per_tree) %>%
      select(-n_tree) %>%
      droplevels()
  }
  
  if (nrow(meta_sub) < 3) {
    warning(scope_label, " | ", dataset_tag, " | ", metric_label, ": too few samples after splitting.")
    return(NULL)
  }
  
  D_full <- fix_dist_labels(D_raw, samps = meta_raw$sample)
  D      <- subset_dist_to_meta(D_full, meta_sub)
  
  stopifnot(identical(attr(D, "Labels"), meta_sub$sample))
  
  if (term_var == "natman" && nlevels(droplevels(factor(meta_sub$natman))) < 2) {
    warning(scope_label, " | ", dataset_tag, " | ", metric_label, ": natman < 2 levels.")
    return(NULL)
  }
  if (term_var == "dw_type2" && nlevels(droplevels(factor(meta_sub$dw_type2))) < 2) {
    warning(scope_label, " | ", dataset_tag, " | ", metric_label, ": dw_type2 < 2 levels.")
    return(NULL)
  }
  
  if (term_var == "natman") {
    meta_sub$natman <- droplevels(factor(meta_sub$natman))
  } else {
    meta_sub$dw_type2 <- droplevels(factor(meta_sub$dw_type2))
  }
  
  use_reads  <- identical(metric_label, "robust.aitchison") && isTRUE(add_log_reads)
  base_terms <- if (use_reads) "log_reads" else character(0)
  
  focal_var <- term_var
  null_var  <- ifelse(term_var == "natman", "nat_null", "type_null")
  
  form_obs <- as.formula(sprintf(
    "D ~ %s%s%s",
    if (length(base_terms)) paste(base_terms, collapse = " + ") else "",
    if (length(base_terms)) " + " else "",
    focal_var
  ))
  
  form_null <- as.formula(sprintf(
    "D ~ %s%s%s",
    if (length(base_terms)) paste(base_terms, collapse = " + ") else "",
    if (length(base_terms)) " + " else "",
    null_var
  ))
  
  perma_obs <- vegan::adonis2(form_obs, data = meta_sub, by = "terms", permutations = N_PERM_NULL)
  print(perma_obs)
  
  R2_obs <- get_R2(perma_obs, focal_var)
  
  null_R2 <- numeric(B_NULL)
  for (b in seq_len(B_NULL)) {
    df_loop <- meta_sub
    df_loop[[null_var]] <- factor(
      sample(df_loop[[focal_var]], replace = FALSE),
      levels = levels(factor(df_loop[[focal_var]]))
    )
    
    fit_null <- tryCatch(
      vegan::adonis2(form_null, data = df_loop, by = "terms", permutations = N_PERM_NULL),
      error = function(e) e
    )
    null_R2[b] <- if (inherits(fit_null, "error")) NA_real_ else get_R2(fit_null, null_var)
    
    if (b %% 20L == 0L) {
      cat("..", scope_label, " ", dataset_tag, " ", metric_label, " ", focal_var,
          " null iteration ", b, "/", B_NULL, "\n", sep = "")
    }
  }
  
  U        <- upper_idx_from(D)
  grp_vec  <- droplevels(factor(meta_sub[[focal_var]]))
  auc_obs  <- auc_from_group(D, grp_vec, U)
  auc_null <- replicate(B_NULL, auc_from_group(D, sample(grp_vec), U)["AUC"])
  
  anos_obs <- vegan::anosim(D, grp_vec, permutations = PERM_ANOSIM)
  print(anos_obs)
  
  pretty_term <- ifelse(focal_var == "natman", "natman", "dw_type2")
  
  p_ecdf <- mk_ecdf_plot(
    D, grp_vec,
    title = sprintf("ECDF | %s | %s | %s | %s", pretty_term, scope_label, dataset_tag, metric_label),
    xlab  = paste(xlab_short, "distance")
  )
  
  p_r2 <- mk_density_plot(
    null_R2, R2_obs,
    title = sprintf("R² null | %s | %s | %s | %s", pretty_term, scope_label, dataset_tag, metric_label),
    xlab  = expression(R^2)
  )
  
  p_auc <- mk_density_plot(
    auc_null, unname(auc_obs["AUC"]),
    title = sprintf("AUC null | %s | %s | %s | %s", pretty_term, scope_label, dataset_tag, metric_label),
    xlab  = "AUC (between > within)"
  )
  
  print(p_ecdf)
  print(p_r2)
  print(p_auc)
  
  if (is.null(save_prefix)) save_prefix <- focal_var
  
  saveRDS(p_r2,  file.path(null_plot_rds_dir, sprintf("R2_%s_%s_%s_%s.rds",   save_prefix, dataset_tag, metric_label, scope_label)))
  saveRDS(p_auc, file.path(null_plot_rds_dir, sprintf("AUC_%s_%s_%s_%s.rds",  save_prefix, dataset_tag, metric_label, scope_label)))
  saveRDS(p_ecdf,file.path(null_plot_rds_dir, sprintf("ECDF_%s_%s_%s_%s.rds", save_prefix, dataset_tag, metric_label, scope_label)))
  
  ciR2  <- if (all(is.na(null_R2))) c(NA_real_, NA_real_) else quantile(null_R2, c(0.025, 0.975), na.rm = TRUE)
  ciAUC <- if (all(is.na(auc_null))) c(NA_real_, NA_real_) else quantile(auc_null, c(0.025, 0.975), na.rm = TRUE)
  
  tibble(
    Scope      = scope_label,
    Dataset    = dataset_tag,
    Metric     = metric_label,
    Term       = focal_var,
    `R2 obs`   = R2_obs,
    `R2 null mean` = mean(null_R2, na.rm = TRUE),
    `R2 null 95% CI` = if (any(!is.na(ciR2))) sprintf("[%.3f, %.3f]", ciR2[1], ciR2[2]) else NA_character_,
    `p (R2)`   = mean(null_R2 >= R2_obs, na.rm = TRUE),
    `R2 SES`   = (R2_obs - mean(null_R2, na.rm = TRUE)) / sd(null_R2, na.rm = TRUE),
    `R2 corrected` = R2_obs - mean(null_R2, na.rm = TRUE),
    `AUC obs`  = unname(auc_obs["AUC"]),
    r_rb       = unname(auc_obs["r_rb"]),
    `AUC null mean` = mean(auc_null, na.rm = TRUE),
    `AUC null 95% CI` = if (any(!is.na(ciAUC))) sprintf("[%.3f, %.3f]", ciAUC[1], ciAUC[2]) else NA_character_,
    `p (AUC)` = mean(auc_null >= auc_obs["AUC"], na.rm = TRUE),
    `AUC SES` = (auc_obs["AUC"] - mean(auc_null, na.rm = TRUE)) / sd(auc_null, na.rm = TRUE),
    `ANOSIM R` = anos_obs$statistic,
    `ANOSIM p` = anos_obs$signif,
    `n trees`  = nlevels(droplevels(factor(meta_sub$natman))),
    `n samples` = nrow(meta_sub)
  )
}

# -----------------------------------------------------------------------------
# 3) CROSS dataset split explicitly by UMI/nonUMI
# -----------------------------------------------------------------------------
# Assumes meta_filt, D_robAit, D_jacc already exist from the refactored CROSS script
stopifnot(exists("meta_filt"), exists("D_robAit"), exists("D_jacc"))

tbl_ait_umi <- run_metric_single_term(
  D_raw = D_robAit, meta_raw = meta_filt,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "UMI", scope_label = "CROSS",
  term_var = "natman", add_log_reads = TRUE, save_prefix = "natman"
)

tbl_jac_umi <- run_metric_single_term(
  D_raw = D_jacc, meta_raw = meta_filt,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "UMI", scope_label = "CROSS",
  term_var = "natman", add_log_reads = FALSE, save_prefix = "natman"
)

tbl_ait_non <- run_metric_single_term(
  D_raw = D_robAit, meta_raw = meta_filt,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "nonUMI", scope_label = "CROSS",
  term_var = "natman", add_log_reads = TRUE, save_prefix = "natman"
)

tbl_jac_non <- run_metric_single_term(
  D_raw = D_jacc, meta_raw = meta_filt,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "nonUMI", scope_label = "CROSS",
  term_var = "natman", add_log_reads = FALSE, save_prefix = "natman"
)

tbl_type_ait_umi <- run_metric_single_term(
  D_raw = D_robAit, meta_raw = meta_filt,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "UMI", scope_label = "CROSS",
  term_var = "dw_type2", add_log_reads = TRUE, save_prefix = "dwtype"
)

tbl_type_jac_umi <- run_metric_single_term(
  D_raw = D_jacc, meta_raw = meta_filt,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "UMI", scope_label = "CROSS",
  term_var = "dw_type2", add_log_reads = FALSE, save_prefix = "dwtype"
)

tbl_type_ait_non <- run_metric_single_term(
  D_raw = D_robAit, meta_raw = meta_filt,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "nonUMI", scope_label = "CROSS",
  term_var = "dw_type2", add_log_reads = TRUE, save_prefix = "dwtype"
)

tbl_type_jac_non <- run_metric_single_term(
  D_raw = D_jacc, meta_raw = meta_filt,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "nonUMI", scope_label = "CROSS",
  term_var = "dw_type2", add_log_reads = FALSE, save_prefix = "dwtype"
)

final_table_long <- bind_rows(
  tbl_ait_umi, tbl_jac_umi, tbl_ait_non, tbl_jac_non,
  tbl_type_ait_umi, tbl_type_jac_umi, tbl_type_ait_non, tbl_type_jac_non
)

final_table_csv2 <- make_summary_wide_csv2(final_table_long, scope_label = "CROSS")
print(final_table_csv2, n = Inf)

write_csv2(
  final_table_csv2,
  file.path(null_table_dir, "nullmodel_summary_natman_dwtype_CROSS_UMI_vs_nonUMI.csv")
)

# -----------------------------------------------------------------------------
# 4) Generic panel assembler
# -----------------------------------------------------------------------------
assemble_null_panels_A4 <- function(null_plot_rds_dir,
                                    out_dir = null_plot_dir,
                                    dataset_tags = c("UMI", "nonUMI"),
                                    metrics = c("robust.aitchison", "jaccard"),
                                    scope_label = "CROSS",
                                    group_key = c("natman", "dwtype"),
                                    dpi = 900) {
  group_key <- match.arg(group_key)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  metric_short <- function(x) ifelse(x == "robust.aitchison", "robAit", x)
  
  strip_numeric_annotations <- function(p, kind) {
    if (!inherits(p, "ggplot")) return(p)
    p <- p + labs(subtitle = NULL) + theme(plot.subtitle = element_blank())
    if (identical(kind, "ECDF")) {
      keep <- vapply(p$layers, function(L) !inherits(L$geom, "GeomText"), logical(1))
      p$layers <- p$layers[keep]
    }
    p
  }
  
  style_obs_vline <- function(p) {
    if (!inherits(p, "ggplot")) return(p)
    if (length(p$layers) >= 2 && inherits(p$layers[[2]]$geom, "GeomVline")) {
      p$layers[[2]]$aes_params$colour   <- "blue"
      p$layers[[2]]$aes_params$linetype <- "twodash"
      p$layers[[2]]$aes_params$linewidth <- 0.8
    }
    p
  }
  
  clean_title <- function(p, dataset, metric_label) {
    if (!inherits(p, "ggplot")) return(p)
    p + labs(title = paste0(dataset, " | ", metric_short(metric_label)))
  }
  
  drop_y_axis <- function(p) {
    if (!inherits(p, "ggplot")) return(p)
    p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  
  getp_raw <- function(kind, dataset_tag, metric_label) {
    f <- file.path(
      null_plot_rds_dir,
      sprintf("%s_%s_%s_%s_%s.rds", kind, group_key, dataset_tag, metric_label, scope_label)
    )
    if (!file.exists(f)) stop("Missing plot file: ", f)
    readRDS(f)
  }
  
  row_label <- function(txt) {
    patchwork::wrap_elements(
      grid::textGrob(txt, rot = 90, gp = grid::gpar(fontsize = 8, fontface = "bold"))
    )
  }
  
  col_specs <- tibble(
    metric  = rep(metrics, each = length(dataset_tags)),
    dataset = rep(dataset_tags, times = length(metrics))
  )
  
  row_kinds <- c("R2", "AUC", "ECDF")
  row_names <- c("R² nulls", "AUC nulls", "ECDF")
  
  all_rows <- vector("list", length(row_kinds))
  
  for (r in seq_along(row_kinds)) {
    rk <- row_kinds[r]
    
    row_plots <- purrr::pmap(col_specs, function(metric, dataset) {
      p <- getp_raw(kind = rk, dataset_tag = dataset, metric_label = metric)
      p <- strip_numeric_annotations(p, rk)
      p <- clean_title(p, dataset = dataset, metric_label = metric)
      if (rk %in% c("AUC", "R2")) p <- style_obs_vline(p)
      p + theme(plot.margin = margin(3, 3, 3, 3))
    })
    
    if (length(row_plots) == 4) {
      row_plots[2:4] <- lapply(row_plots[2:4], drop_y_axis)
    }
    
    all_rows[[r]] <- c(list(row_label(row_names[r])), row_plots)
  }
  
  grob_list <- unlist(all_rows, recursive = FALSE)
  
  title_main <- if (group_key == "natman") {
    paste0("Null + separation diagnostics: natman (", scope_label, ")")
  } else {
    paste0("Null + separation diagnostics: dw_type2 (", scope_label, ")")
  }
  
  subtitle_main <- "Columns: robAit (UMI, nonUMI) then jaccard (UMI, nonUMI). Rows: R² nulls, AUC nulls, ECDF."
  
  p_out <- patchwork::wrap_plots(
    grob_list,
    ncol = 5,
    byrow = TRUE,
    widths = c(0.08, 1, 1, 1, 1)
  ) +
    patchwork::plot_annotation(
      title = title_main,
      subtitle = subtitle_main
    ) &
    theme(plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 9))
  
  a4_w <- 9
  a4_h <- 11.69 / 2.25
  
  outfile <- sprintf("PANEL_%s_%s_A4.png", scope_label, group_key)
  ggsave(file.path(out_dir, outfile), p_out, width = a4_w, height = a4_h, dpi = dpi)
  
  print(p_out)
  invisible(p_out)
}

assemble_null_panels_A4(
  null_plot_rds_dir = null_plot_rds_dir,
  out_dir = null_plot_dir,
  scope_label = "CROSS",
  group_key = "natman"
)

assemble_null_panels_A4(
  null_plot_rds_dir = null_plot_rds_dir,
  out_dir = null_plot_dir,
  scope_label = "CROSS",
  group_key = "dwtype"
)

# -----------------------------------------------------------------------------
# 5) LOGS natman-only protocol split
# -----------------------------------------------------------------------------
scope_label <- "LOGS"
scope_dw    <- "LOG"

meta <- build_intra_meta_ds_at_drill(scope_dw)

meta_ext <- meta %>%
  filter(!is.na(diameter_at_drill_z)) %>%
  droplevels()

stopifnot("umi_backup" %in% names(meta_ext))

if (!("sample" %in% names(meta_ext))) {
  meta_ext <- meta_ext %>%
    mutate(sample = if (!is.null(rownames(meta_ext))) rownames(meta_ext) else as.character(seq_len(nrow(meta_ext))))
}

otu <- otu_matrix_filt[meta_ext$sample, , drop = FALSE]
stopifnot(identical(rownames(otu), meta_ext$sample))
otu <- otu[, colSums(otu) > 0, drop = FALSE]

cat("\nLOGS meta_ext natman counts (pre-split):\n")
print(table(meta_ext$natman))

D_robAit_log <- vegan::vegdist(otu, method = "robust.aitchison")
D_robAit_log <- fix_dist_labels(D_robAit_log, samps = meta_ext$sample)

otu_pa_log <- (otu > 0) * 1L
D_jacc_log <- vegan::vegdist(otu_pa_log, method = "jaccard", binary = TRUE)

tbl_log_ait_umi <- run_metric_single_term(
  D_raw = D_robAit_log, meta_raw = meta_ext,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "UMI", scope_label = "LOGS",
  term_var = "natman", add_log_reads = TRUE,
  min_per_tree = 3L, save_prefix = "natman"
)

tbl_log_jac_umi <- run_metric_single_term(
  D_raw = D_jacc_log, meta_raw = meta_ext,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "UMI", scope_label = "LOGS",
  term_var = "natman", add_log_reads = FALSE,
  min_per_tree = 3L, save_prefix = "natman"
)

tbl_log_ait_non <- run_metric_single_term(
  D_raw = D_robAit_log, meta_raw = meta_ext,
  metric_label = "robust.aitchison", xlab_short = "Aitchison",
  dataset_tag = "nonUMI", scope_label = "LOGS",
  term_var = "natman", add_log_reads = TRUE,
  min_per_tree = 3L, save_prefix = "natman"
)

tbl_log_jac_non <- run_metric_single_term(
  D_raw = D_jacc_log, meta_raw = meta_ext,
  metric_label = "jaccard", xlab_short = "Jaccard",
  dataset_tag = "nonUMI", scope_label = "LOGS",
  term_var = "natman", add_log_reads = FALSE,
  min_per_tree = 3L, save_prefix = "natman"
)

final_log_long <- bind_rows(
  tbl_log_ait_umi,
  tbl_log_jac_umi,
  tbl_log_ait_non,
  tbl_log_jac_non
)

final_log_csv2 <- make_summary_wide_csv2(final_log_long, scope_label = "LOGS")
print(final_log_csv2, n = Inf)

write_csv2(
  final_log_csv2,
  file.path(null_table_dir, "nullmodel_summary_natman_LOGS_UMI_vs_nonUMI.csv")
)

assemble_null_panels_A4(
  null_plot_rds_dir = null_plot_rds_dir,
  out_dir = null_plot_dir,
  scope_label = "LOGS",
  group_key = "natman"
)