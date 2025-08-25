## Harmonized BETWEEN-LOG builder
## Datasets produced:
##   1) MID-only   : otu_mid_only,  meta_mid_only   (per log, OUTER+INNER/MIXED)
##   2) LOG+aFWD   : otu_log_afwd,  meta_log_afwd   (triplet-only logs: OUTER+INNER+ATTACHED)
## Plots produced:
##   - Alpha diversity (by ds_at_drill, size) for both datasets
##   - Rarefaction (one line per sample) + single collectors curve for both datasets
## --------------------------------------------------------------------------------
source("scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tibble)
  library(tidyr); library(ggplot2); library(vegan); library(stringr)
})

# --- helpers -------------------------------------------------------------------
gm_mean <- function(x, pseudo = 1) exp(mean(log(x + pseudo))) - pseudo
combine_counts <- function(mat, rows, method = c("geo","median"), pseudo = 1) {
  method <- match.arg(method)
  X <- mat[rows, , drop = FALSE]
  if (nrow(X) == 1L) return(as.numeric(X[1, ]))
  if (method == "median") return(apply(X, 2, median))
  apply(X, 2, function(col) gm_mean(col, pseudo = pseudo))
}
mode_or_median <- function(x) {
  if (is.numeric(x)) return(median(x, na.rm = TRUE))
  ux <- na.omit(as.character(x)); if (!length(ux)) return(NA_character_)
  tab <- sort(table(ux), decreasing = TRUE)
  names(tab)[1]
}
clean_nat <- function(x) gsub('"', "", x)

# =================================================================================
# 1) Select LOG/MID samples and pick per-log representatives (OUTER/INNER/MIXED)
# =================================================================================
df0 <- META1 %>%
  filter(dw_type == "LOG",
         position %in% c("MIDDLE","MID"),
         depth_2 %in% c("OUTER","INNER","MIXED")) %>%
  select(sample, natman, position, depth_2, reads_filt, ds_at_drill, size) %>%
  mutate(natman = clean_nat(natman)) %>% 
  rename(dw_size = size)

pick_best <- function(g) {
  if (any(g$depth_2 == "MIXED")) {
    warning("Using MIXED core instead of OUTER/INNER for log ", unique(g$natman))
    g %>% filter(depth_2 == "MIXED") %>%
      arrange(desc(reads_filt)) %>%
      slice(1) %>%
      mutate(role = "MIXED", set_type = "single")
  } else {
    rows <- list()
    if (any(g$depth_2 == "OUTER")) {
      rows <- append(rows, list(g %>% filter(depth_2 == "OUTER") %>%
                                  arrange(desc(reads_filt)) %>% slice(1) %>%
                                  mutate(role = "OUTER")))
    }
    if (any(g$depth_2 == "INNER")) {
      rows <- append(rows, list(g %>% filter(depth_2 == "INNER") %>%
                                  arrange(desc(reads_filt)) %>% slice(1) %>%
                                  mutate(role = "INNER")))
    }
    if (!length(rows)) return(tibble())
    bind_rows(rows) %>%
      mutate(set_type = ifelse(n() == 2, "doublet", "singleton"))
  }
}

sel_tbl <- df0 %>%
  group_split(natman, .keep = TRUE) %>%
  map_dfr(pick_best) %>%
  arrange(natman, role) %>%
  select(natman, role, set_type, sample)

# =================================================================================
# 2) ATTACHED aFWD for the same trees + resolve duplicates by reads_filt
# =================================================================================
trees <- unique(sel_tbl$natman)

attached_raw <- META1 %>%
  mutate(natman = clean_nat(natman)) %>%
  filter(natman %in% trees, dw_type == "FWD", position == "ATTACHED") %>%
  select(natman, sample, reads_filt)

attached_tbl <- attached_raw %>%
  group_by(natman) %>%
  arrange(desc(reads_filt), .by_group = TRUE) %>%
  { 
    dup_counts <- summarise(., n = n()) %>% filter(n > 1)
    if (nrow(dup_counts)) {
      warning("Multiple ATTACHED found for: ",
              paste(dup_counts$natman, collapse=", "),
              " — choosing highest reads_filt.")
    }
    slice_head(., n = 1)
  } %>%
  ungroup() %>%
  mutate(role = "ATTACHED", set_type = "triplet") %>%
  select(natman, role, set_type, sample)

# Logs that truly have OUTER + INNER (for defining triplets later)
logs_with_pair <- sel_tbl %>%
  group_by(natman) %>%
  summarise(has_outer = any(role == "OUTER"),
            has_inner = any(role == "INNER"),
            .groups = "drop") %>%
  filter(has_outer & has_inner) %>%
  select(natman)

# =================================================================================
# 3) META/OTU for selected samples
# =================================================================================
meta_between <- META1 %>%
  mutate(natman = clean_nat(natman)) %>%
  semi_join(sel_tbl, by = "sample") %>%
  select(sample, natman, ds_at_drill, size) %>%
  left_join(sel_tbl, by = c("sample","natman")) %>%
  arrange(natman, role) %>% 
  rename(dw_size = size)

otu_between <- otu_matrix_filt[meta_between$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_between), meta_between$sample))

# =================================================================================
# 4) Aggregate to MID (combine OUTER+INNER per log; MIXED/singletons pass)
# =================================================================================
mid_map <- sel_tbl %>%
  group_by(natman) %>%
  summarise(samples_mid = list(sample), .groups = "drop")

natmans_mid <- mid_map$natman

otu_mid <- t(vapply(seq_len(nrow(mid_map)), function(i) {
  combine_counts(otu_between, rows = mid_map$samples_mid[[i]], method = "geo", pseudo = 1)
}, FUN.VALUE = numeric(ncol(otu_between))))
rownames(otu_mid) <- paste0(natmans_mid, "::MID")

meta_mid <- meta_between %>%
  semi_join(sel_tbl, by = "sample") %>%
  group_by(natman) %>%
  summarise(
    ds_at_drill = mode_or_median(ds_at_drill),
    dw_size        = mode_or_median(dw_size),
    n_cores     = n(),
    roles       = paste(sort(unique(role)), collapse = "+"),
    .groups     = "drop"
  ) %>%
  mutate(compartment = "MID",
         rowname = paste0(natman, "::MID"))

# Harmonize types (use META1 as reference)
ds_lv <- levels(factor(META1$ds_at_drill))
size_lv <- if (is.factor(META1$size)) levels(META1$size) else NULL

meta_mid <- meta_mid %>%
  mutate(
    ds_at_drill = factor(as.character(ds_at_drill), levels = ds_lv),
    dw_size = if (!is.null(size_lv)) {
      factor(as.character(dw_size), levels = size_lv)
    } else {
      suppressWarnings(as.numeric(dw_size))
    }
  )
# MID-only objects (aligned)
otu_mid_only  <- otu_mid
meta_mid_only <- meta_mid %>%
  arrange(match(rowname, rownames(otu_mid_only))) %>%
  select(-rowname)

stopifnot(identical(rownames(otu_mid_only), paste0(meta_mid$natman, "::MID")))

# =================================================================================
# 5) TRIPLET-ONLY aggregation: OUTER + INNER + ATTACHED (≈ 18 logs)
# =================================================================================
triplet_map <- sel_tbl %>%
  filter(role %in% c("OUTER","INNER")) %>%
  group_by(natman) %>%
  summarise(samples_mid = list(sample), n_mid = n(), .groups = "drop") %>%
  filter(n_mid == 2) %>%  # must have both OUTER & INNER
  inner_join(attached_tbl %>% select(natman, sample_att = sample), by = "natman") %>%
  mutate(samples_all = purrr::map2(samples_mid, sample_att, ~ c(.x, .y)))

message("Triplet logs found: ", nrow(triplet_map), " (expected ~18)")

# OTU source for all rows we will combine (only for the triplets)
sel_all_samples <- unique(c(unlist(triplet_map$samples_mid), triplet_map$sample_att))
otu_trip_src <- otu_matrix_filt[sel_all_samples, , drop = FALSE]

# Aggregate to one vector per triplet log
otu_log_afwd <- t(vapply(seq_len(nrow(triplet_map)), function(i) {
  combine_counts(otu_trip_src,
                 rows = c(triplet_map$samples_mid[[i]], triplet_map$sample_att[[i]]),
                 method = "geo", pseudo = 1)
}, FUN.VALUE = numeric(ncol(otu_trip_src))))
rownames(otu_log_afwd) <- paste0(triplet_map$natman, "::LOGaFWD")
otu_log_afwd <- otu_log_afwd[, colSums(otu_log_afwd) > 0, drop = FALSE]

# Metadata for triplet aggregation
# ds_at_drill = OUTER core's value (fallback: MID)
outer_ds_map <- sel_tbl %>%
  filter(role == "OUTER", natman %in% triplet_map$natman) %>%
  select(natman, sample_outer = sample) %>%
  left_join(META1 %>% select(sample, ds_outer = ds_at_drill), by = c("sample_outer" = "sample"))

meta_log_afwd <- meta_mid %>%
  filter(natman %in% triplet_map$natman) %>%
  select(natman, ds_mid = ds_at_drill, size_mid = dw_size) %>%
  left_join(outer_ds_map, by = "natman") %>%
  mutate(
    ds_final_chr = dplyr::coalesce(as.character(ds_outer), as.character(ds_mid)),
    ds_at_drill  = factor(ds_final_chr, levels = ds_lv),
    dw_size         = size_mid,
    n_cores      = 3L,
    roles        = "INNER+OUTER+ATTACHED",
    compartment  = "LOG_aFWD",
    rowname      = paste0(natman, "::LOGaFWD")
  ) %>%
  arrange(match(rowname, rownames(otu_log_afwd))) %>%
  select(-rowname)

stopifnot(nrow(meta_log_afwd) == nrow(otu_log_afwd))

# =================================================================================
# 6) Alpha diversity plots (MID-only and LOG+aFWD)
# =================================================================================
outdir_alpha <- "plots/BETWEEN_LOG/ALPHA"
dir.create(outdir_alpha, recursive = TRUE, showWarnings = FALSE)
has_DS_cols  <- exists("ds_colors")

compute_alpha <- function(otu_mat) {
  X <- as.matrix(otu_mat); X[X < 0] <- 0
  lib <- rowSums(X)
  S_obs <- rowSums(X > 0)
  P <- sweep(X, 1, ifelse(lib > 0, lib, 1), "/")
  H <- diversity(P, index = "shannon")
  invS <- diversity(P, index = "invsimpson")
  J <- ifelse(S_obs > 1, H / log(S_obs), NA_real_)
  tibble(
    rowname = rownames(X),
    reads   = lib,
    S_obs   = as.numeric(S_obs),
    H       = as.numeric(H),
    q1_expH = exp(H),
    invSimpson = as.numeric(invS),
    Evenness  = as.numeric(J)
  )
}

plot_alpha_by <- function(otu_mat, meta_df, covariate, dataset_label, outdir = outdir_alpha) {
  # Ensure a rowname that matches OTU
  if (!("rowname" %in% names(meta_df))) {
    meta_df <- meta_df %>%
      mutate(rowname = paste0(natman, "::", compartment))
  }
  stopifnot(identical(rownames(otu_mat), meta_df$rowname))
  
  alpha_df <- compute_alpha(otu_mat) %>%
    left_join(meta_df, by = "rowname")
  
  long <- alpha_df %>%
    select(rowname, !!covariate, reads, S_obs, H, q1_expH, invSimpson, Evenness) %>%
    pivot_longer(cols = c(reads, S_obs, H, q1_expH, invSimpson, Evenness),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
                           levels = c("reads","S_obs","H","q1_expH","invSimpson","Evenness"),
                           labels = c("Reads","Richness (S)","Shannon (H)","Hill q=1 (exp H)",
                                      "Hill q=2 (Inv. Simpson)","Evenness (Pielou)")))
  x <- long[[covariate]]
  is_num <- is.numeric(x) || is.integer(x)
  
  p <- ggplot(long, aes(x = .data[[covariate]], y = value)) +
    { if (is_num) {
      list(
        geom_point(alpha = 0.7, position = position_jitter(width = 0.05, height = 0), size = 1.8),
        suppressWarnings(geom_smooth(method = "gam", formula = y ~ s(x, k = 3), se = TRUE))
      )
    } else {
      list(
        geom_boxplot(outlier.shape = NA, width = 0.6),
        geom_jitter(width = 0.15, alpha = 0.5, size = 1.6)
      )
    }
    } +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    labs(x = covariate, y = NULL,
         title = paste0("Alpha diversity by ", covariate, " (", dataset_label, ")"),
         subtitle = if (is_num) "Points: samples; Line: GAM smooth (k=3)" else "Boxplots with jittered points") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
          axis.title.x = element_text(margin = margin(t = 6)))
  
  if (!is_num && covariate == "ds_at_drill") {
    p <- p + scale_x_discrete(drop = FALSE) +
      scale_color_manual(values = ds_colors, aesthetics = c("colour")) +
      scale_fill_manual(values = ds_colors)
  }
  
  fname <- file.path(outdir, paste0("alpha_", dataset_label, "_by_", covariate, ".png"))
  ggsave(filename = fname, plot = p, width = 10, height = 7, dpi = 300)
  message("Saved: ", fname)
  print(p)
  invisible(p)
}

# -- MID-only alpha
meta_mid_only_row <- meta_mid_only %>% mutate(rowname = paste0(natman, "::MID"),
                                              compartment = "MID")
if ("ds_at_drill" %in% names(meta_mid_only_row)) plot_alpha_by(otu_mid_only, meta_mid_only_row, "ds_at_drill", "MID_only")
if ("dw_size"       %in% names(meta_mid_only_row)) plot_alpha_by(otu_mid_only, meta_mid_only_row, "dw_size",       "MID_only")

# -- LOG+aFWD alpha (triplet-only)
meta_log_afwd_row <- meta_log_afwd %>% mutate(rowname = paste0(natman, "::LOG_aFWD"),
                                              compartment = "LOG_aFWD")
rownames(otu_log_afwd) <- paste0(gsub("::.*$", "", rownames(otu_log_afwd)), "::LOG_aFWD")
if ("ds_at_drill" %in% names(meta_log_afwd_row)) plot_alpha_by(otu_log_afwd, meta_log_afwd_row, "ds_at_drill", "LOG_aFWD")
if ("dw_size"       %in% names(meta_log_afwd_row)) plot_alpha_by(otu_log_afwd, meta_log_afwd_row, "dw_size",       "LOG_aFWD")

# =================================================================================
# 7) Rarefaction (one line per sample) + single collectors curve per dataset
# =================================================================================
outdir_rc_simple <- "plots/BETWEEN_LOG/RAREFACTION_SIMPLE"
dir.create(outdir_rc_simple, recursive = TRUE, showWarnings = FALSE)

integerize_counts <- function(X, target_sum = NULL, min_target = 2000) {
  X <- as.matrix(X); X[X < 0] <- 0
  libs <- rowSums(X)
  P <- sweep(X, 1, ifelse(libs > 0, libs, 1), "/")
  if (is.null(target_sum)) {
    target_sum <- max(min_target, round(median(libs[libs > 0], na.rm = TRUE)))
    if (!is.finite(target_sum) || target_sum <= 0) target_sum <- min_target
  }
  round(sweep(P, 1, target_sum, "*"))
}

rarecurve_df <- function(X_int, step = 500) {
  stopifnot(all(X_int == round(X_int)))
  rc_list <- vegan::rarecurve(X_int, step = step, label = FALSE, return = "list")
  names(rc_list) <- rownames(X_int)
  purrr::imap_dfr(rc_list, function(vals, rn) {
    tibble(rowname = rn,
           reads = attr(vals, "Subsample"),
           richness = as.integer(vals))
  })
}

collectors_all_df <- function(X, perms = 1000) {
  acc <- specaccum(X > 0, method = "random", permutations = perms)
  tibble(sites = acc$sites, richness = acc$richness, sd = acc$sd)
}

plot_rarecurve_lines <- function(rc_df, dataset_label, outdir = outdir_rc_simple) {
  p <- ggplot(rc_df, aes(x = reads, y = richness, group = rowname)) +
    geom_path(alpha = 0.55, linewidth = 0.45) +
    labs(title = paste0("Sample rarefaction (", dataset_label, ")"),
         x = "Reads (subsampled)", y = "Observed richness (S)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  f <- file.path(outdir, paste0("rarecurve_", dataset_label, ".png"))
  ggsave(f, p, width = 10, height = 7, dpi = 300)
  message("Saved: ", f)
  print(p); invisible(p)
}

plot_collectors_single <- function(cdf, dataset_label, outdir = outdir_rc_simple) {
  p <- ggplot(cdf, aes(x = sites, y = richness)) +
    geom_ribbon(aes(ymin = pmax(richness - sd, 0), ymax = richness + sd), alpha = 0.18) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0("Collector’s curve — all samples (", dataset_label, ")"),
         x = "Number of sampling units", y = "Cumulative richness (S)") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  f <- file.path(outdir, paste0("collectors_", dataset_label, ".png"))
  ggsave(f, p, width = 9, height = 6, dpi = 300)
  message("Saved: ", f)
  print(p); invisible(p)
}

# -- MID-only
otu_mid_only_int <- integerize_counts(otu_mid_only, target_sum = NULL)
rc_mid_only <- rarecurve_df(otu_mid_only_int, step = 1000)
plot_rarecurve_lines(rc_mid_only, dataset_label = "MID_only")
if (nrow(otu_mid_only) >= 2) {
  cdf_mid_only <- collectors_all_df(otu_mid_only, perms = 1000)
  plot_collectors_single(cdf_mid_only, dataset_label = "MID_only")
} else message("Not enough samples for collector’s curve (MID_only).")

# -- LOG+aFWD (triplet-only)
otu_log_afwd_int <- integerize_counts(otu_log_afwd, target_sum = NULL)
rc_log_afwd <- rarecurve_df(otu_log_afwd_int, step = 1000)
plot_rarecurve_lines(rc_log_afwd, dataset_label = "LOG_aFWD")
if (nrow(otu_log_afwd) >= 2) {
  cdf_log_afwd <- collectors_all_df(otu_log_afwd, perms = 1000)
  plot_collectors_single(cdf_log_afwd, dataset_label = "LOG_aFWD")
} else message("Not enough samples for collector’s curve (LOG_aFWD).")


## Beta-diversity suite for BETWEEN-LOG datasets
## Expects:
##   - otu_mid_only,  meta_mid_only      (rows "<natman>::MID")
##   - otu_log_afwd,  meta_log_afwd      (rows "<natman>::LOG_aFWD")
## Distances:
##   1) Hellinger + Euclidean
##   2) Robust Aitchison  (vegdist(method = "robust.aitchison"))
##   3) Presence/Absence + Jaccard
## Ordinations: NMDS (monoMDS) + PCoA; PCA (RDA) biplot on Hellinger; LCBD/SCBD
## Outputs: plots -> plots/BETWEEN_LOG/BETA ; stats -> results/BETWEEN_LOG/BETA
## --------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble)
  library(ggplot2); library(ggrepel)
  library(vegan); library(adespatial); library(readr)
})

# ---------------- paths ----------------
out_plot <- "plots/BETWEEN_LOG/BETA"
out_res  <- "results/BETWEEN_LOG/BETA"
dir.create(out_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(out_res,  recursive = TRUE, showWarnings = FALSE)

# ---------------- helpers ----------------
as_matrix_counts <- function(X) { X <- as.matrix(X); X[X < 0] <- 0; X }

# Transform matrices
X_hellinger <- function(X) vegan::decostand(as_matrix_counts(X), method = "hellinger")
X_presence  <- function(X) (as_matrix_counts(X) > 0) * 1
# Robust CLR via vegan (for diagnostics/plots if needed)
X_rclr      <- function(X) vegan::decostand(as_matrix_counts(X), method = "rclr")  # robust CLR

# PCoA (cmdscale)
pcoa_scores <- function(d) {
  fit <- cmdscale(d, k = 2, eig = TRUE)
  evp <- fit$eig; varex <- ifelse(evp > 0, evp / sum(evp[evp > 0]), 0)
  tibble(Axis1 = fit$points[,1], Axis2 = fit$points[,2],
         var1 = 100 * varex[1], var2 = 100 * varex[2],
         rowname = rownames(fit$points))
}

# NMDS
nmds_scores <- function(d, trymax = 100) {
  fit <- vegan::monoMDS(d, k = 2, trymax = trymax, trace = FALSE)
  tibble(Axis1 = fit$points[,1], Axis2 = fit$points[,2],
         stress = fit$stress, rowname = rownames(fit$points))
}

# PCA (RDA) biplot on Hellinger
rda_biplot_data <- function(X_hel, scaling = 2) {
  r <- vegan::rda(X_hel)  # PCA of Hellinger data
  ssit <- as_tibble(vegan::scores(r, display = "sites",   scaling = scaling), rownames = "rowname")
  sspc <- as_tibble(vegan::scores(r, display = "species", scaling = scaling), rownames = "feature")
  eig  <- r$CA$eig
  list(site = ssit %>% rename(Axis1 = PC1, Axis2 = PC2),
       species = sspc %>% rename(Axis1 = PC1, Axis2 = PC2),
       var = 100 * eig[1:2] / sum(eig))
}

# Generic scatter plot
plot_scatter <- function(df, meta, label, method_lab,
                         aes_color = "ds_at_drill", aes_shape = "dw_size",
                         xlab = "Axis 1", ylab = "Axis 2", subtitle_extra = NULL,
                         fn_prefix = "ORD") {
  stopifnot("rowname" %in% names(meta))
  df2 <- df %>% left_join(meta, by = "rowname")
  p <- ggplot(df2, aes(x = Axis1, y = Axis2)) +
    geom_point(aes(color = .data[[aes_color]], shape = .data[[aes_shape]]), size = 2, alpha = 0.9) +
    labs(title = paste0(method_lab, " — ", label),
         subtitle = subtitle_extra, x = xlab, y = ylab,
         color = aes_color, shape = aes_shape) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  fn <- file.path(out_plot, paste0(fn_prefix, "_", gsub("[^A-Za-z0-9]+","_", method_lab), "_", label, ".png"))
  ggsave(fn, p, width = 8.5, height = 6.5, dpi = 300)
  message("Saved: ", fn)
  print(p)
}

# Biplot with top-k species (by SCBD if supplied)
plot_rda_biplot <- function(X_hel, meta, label, top_k = 25, scbd = NULL) {
  dat <- rda_biplot_data(X_hel, scaling = 2)
  sites <- dat$site %>% left_join(meta, by = "rowname")
  print(head(sites))

  sp    <- dat$species
  if (!is.null(scbd)) {
    keep <- names(sort(scbd, decreasing = TRUE))[seq_len(min(top_k, length(scbd)))]
    sp <- sp %>% filter(feature %in% keep)
  } else {
    v <- apply(X_hel, 2, var); keep <- names(sort(v, decreasing = TRUE))[seq_len(min(top_k, length(v)))]
    sp <- sp %>% filter(feature %in% keep)
  }
  mult <- 1.0 * max(abs(c(sites$Axis1, sites$Axis2)), na.rm = TRUE) /
    max(abs(c(sp$Axis1, sp$Axis2)), na.rm = TRUE)
  sp <- sp %>% mutate(xend = Axis1 * mult, yend = Axis2 * mult)
  p <- ggplot() +
    geom_point(data = sites, aes(x = Axis1, y = Axis2, color = ds_at_drill, shape = dw_size), size = 2) +
    geom_segment(data = sp, aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.15, "cm")), alpha = 0.6) +
    ggrepel::geom_text_repel(data = sp, aes(x = xend, y = yend, label = feature), size = 3, max.overlaps = 60) +
    labs(title = paste0("PCA (RDA) biplot on Hellinger — ", label),
         subtitle = paste0("Axis1 ", round(dat$var[1],1), "%, Axis2 ", round(dat$var[2],1), "%"),
         x = "PC1", y = "PC2", color = "ds_at_drill", shape = "dw size") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  fn <- file.path(out_plot, paste0("BIPLOT_RDA_Hellinger_", gsub("[^A-Za-z0-9]+","_", label), ".png"))
  ggsave(fn, p, width = 9.5, height = 7.5, dpi = 300)
  message("Saved: ", fn)
}

# PERMANOVA + dispersion

run_permanova <- function(d, meta, label, dist_lab) {
  stopifnot("rowname" %in% names(meta))
  mm <- meta %>%
    mutate(
      ds_at_drill = as.factor(ds_at_drill),
      dw_size = if (is.numeric(dw_size) || is.integer(dw_size)) dw_size else as.factor(dw_size)
    )
  keep <- complete.cases(mm[, c("ds_at_drill", "dw_size")])
  if (any(!keep)) {
    dropped <- mm$rowname[!keep]
    warning("PERMANOVA (", dist_lab, " — ", label, "): dropping ",
            sum(!keep), " rows with NA in predictors: ",
            paste0(head(dropped, 10), if (length(dropped) > 10) " ..."))
  }
  mm2 <- mm[keep, , drop = FALSE]
  if (nrow(mm2) < 3 || length(unique(mm2$ds_at_drill)) < 2) {
    message("PERMANOVA (", dist_lab, " — ", label, "): skipped (need ≥3 samples and ≥2 groups).")
    return(invisible(NULL))
  }
  
  D <- as.matrix(d)
  # subset rows/cols using the *row order* of mm2 (already aligned to OTU/meta)
  D2 <- as.dist(D[match(mm2$rowname, labels(d)), match(mm2$rowname, labels(d))])
  
  mod <- vegan::adonis2(D2 ~ ds_at_drill + dw_size, data = mm2, permutations = 999, by = "margin")
  fn  <- file.path(out_res, paste0("PERMANOVA_", gsub("[^A-Za-z0-9]+","_", dist_lab), "_", label, ".csv"))
  readr::write_csv(as.data.frame(mod), fn)
  cat("\npermanova results\n")
  print(as.data.frame(mod))
  message("Saved: ", fn)
  mod
}

run_betadisper <- function(d, grouping, label, dist_lab) {
  # 'grouping' should be in the same order as the distance labels (as in your suite)
  labs <- labels(d)
  stopifnot(length(grouping) == length(labs))
  
  keep <- !is.na(grouping)
  if (any(!keep)) {
    warning("BETADISPER (", dist_lab, " — ", label, "): dropping ",
            sum(!keep), " rows with NA in grouping.")
  }
  grouping2 <- factor(grouping[keep])
  if (length(grouping2) < 3 || nlevels(grouping2) < 2) {
    message("BETADISPER (", dist_lab, " — ", label, "): skipped (need ≥3 samples and ≥2 groups).")
    return(invisible(NULL))
  }
  
  D <- as.matrix(d)
  D2 <- as.dist(D[keep, keep])
  labs2 <- labs[keep]
  
  bd <- vegan::betadisper(D2, grouping2)
  an <- anova(bd)
  
  # Save ANOVA table
  fn_tab <- file.path(out_res, paste0("BETADISPER_", gsub("[^A-Za-z0-9]+","_", dist_lab), "_", label, ".csv"))
  readr::write_csv(as.data.frame(an), fn_tab)
  message("Saved: ", fn_tab)
  cat("\n anova table from BETADISPER \n")
  print(as.data.frame(an))
  # Plot distances to centroid
  df <- tibble::tibble(
    rowname = labs2,
    distance = as.numeric(bd$distances),
    group = grouping2
  )
  p <- ggplot(df, aes(x = group, y = distance, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.12, alpha = 0.6, size = 1.6) +
    labs(title = paste0("Dispersion (", dist_lab, ") — ", label),
         x = "", y = "Distance to centroid") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "none")
  print(p)
  fn_plot <- file.path(out_plot, paste0("BETADISPER_", gsub("[^A-Za-z0-9]+","_", dist_lab), "_", label, ".png"))
  ggsave(fn_plot, p, width = 7.5, height = 5.5, dpi = 300)
  message("Saved: ", fn_plot)
  
  invisible(list(betadisper = bd, anova = an))
}


# adespatial LCBD/SCBD on Hellinger
run_beta_div_adespatial <- function(X_hel, label) {
  bd <- adespatial::beta.div(X_hel, method = "hellinger", sqrt.D = FALSE)
  lcbd <- tibble(rowname = rownames(X_hel), LCBD = bd$LCBD, p.LCBD = bd$p.LCBD)
  write_csv(lcbd, file.path(out_res, paste0("LCBD_", label, ".csv")))
  scbd <- sort(bd$SCBD, decreasing = TRUE)
  top <- head(scbd, n = min(20, length(scbd)))
  df <- tibble(feature = factor(names(top), levels = names(top)), SCBD = as.numeric(top))
  cat("df for SCBD\n")
  print(df)
  p <- ggplot(df, aes(x = feature, y = SCBD)) +
    geom_col() + coord_flip() +
    labs(title = paste0("Top species SCBD — ", label), x = NULL, y = "SCBD") +
    theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
  fn <- file.path(out_plot, paste0("SCBD_top_", gsub("[^A-Za-z0-9]+","_", label), ".png"))
  ggsave(fn, p, width = 7.5, height = 7.5, dpi = 300)
  message("Saved: ", fn)
  invisible(list(bd = bd, scbd = scbd, lcbd = lcbd))
  print(p)
}

# ------------ master runner for one dataset ------------
run_beta_suite <- function(otu, meta, label) {
  # Make sure meta has 'rowname' matching OTU rownames
  if (!("rowname" %in% names(meta))) stop("meta must contain 'rowname'.")
  stopifnot(identical(rownames(otu), meta$rowname))
  
  # ---- Distances ----
  Xhel <- X_hellinger(otu)
  d_hel <- dist(Xhel, method = "euclidean")
  
  # Robust Aitchison directly in vegan:
  d_ait <- vegan::vegdist(otu, method = "robust.aitchison")
  
  # Presence/absence + Jaccard
  d_jac <- vegan::vegdist(X_presence(otu), method = "jaccard", binary = TRUE)
  
  # Save distances
  saveRDS(list(hel = d_hel, robust_aitchison = d_ait, jaccard = d_jac),
          file = file.path(out_res, paste0("DIST_", label, ".rds")))
  
  # ---- Ordinations ----
  # Hellinger
  sc_nmds_h <- nmds_scores(d_hel)
  plot_scatter(sc_nmds_h, meta, label, "NMDS (Hellinger+Euclid)",
               subtitle_extra = paste0("Stress = ", round(unique(sc_nmds_h$stress), 3)),
               fn_prefix = "NMDS")
  sc_pcoa_h <- pcoa_scores(d_hel)
  print(plot_scatter(sc_pcoa_h, meta, label, "PCoA (Hellinger+Euclid)",
               xlab = paste0("Axis 1 (", round(unique(sc_pcoa_h$var1),1), "%)"),
               ylab = paste0("Axis 2 (", round(unique(sc_pcoa_h$var2),1), "%)"),
               fn_prefix = "PCoA"))
  
  # Robust Aitchison
  sc_nmds_a <- nmds_scores(d_ait)
  print(plot_scatter(sc_nmds_a, meta, label, "NMDS (Robust Aitchison)",
               subtitle_extra = paste0("Stress = ", round(unique(sc_nmds_a$stress), 3)),
               fn_prefix = "NMDS"))
  sc_pcoa_a <- pcoa_scores(d_ait)
  print(plot_scatter(sc_pcoa_a, meta, label, "PCoA (Robust Aitchison)",
               xlab = paste0("Axis 1 (", round(unique(sc_pcoa_a$var1),1), "%)"),
               ylab = paste0("Axis 2 (", round(unique(sc_pcoa_a$var2),1), "%)"),
               fn_prefix = "PCoA"))
  
  # Presence/Absence + Jaccard
  sc_nmds_j <- nmds_scores(d_jac)
  print(plot_scatter(sc_nmds_j, meta, label, "NMDS (Presence/Absence + Jaccard)",
               subtitle_extra = paste0("Stress = ", round(unique(sc_nmds_j$stress), 3)),
               fn_prefix = "NMDS"))
  sc_pcoa_j <- pcoa_scores(d_jac)
  print(plot_scatter(sc_pcoa_j, meta, label, "PCoA (Presence/Absence + Jaccard)",
               xlab = paste0("Axis 1 (", round(unique(sc_pcoa_j$var1),1), "%)"),
               ylab = paste0("Axis 2 (", round(unique(sc_pcoa_j$var2),1), "%)"),
               fn_prefix = "PCoA"))
  
  # ---- Species biplot (PCA on Hellinger) + LCBD/SCBD ----
  bd <- run_beta_div_adespatial(Xhel, label)
  print(plot_rda_biplot(Xhel, meta, label, top_k = 8, scbd = bd$scbd))
  
  
  # ---- Tests ----
  run_permanova(d_hel, meta, label, "Hellinger+Euclid")
  run_permanova(d_ait, meta, label, "RobustAitchison")
  run_permanova(d_jac, meta, label, "Jaccard")
  
  if ("ds_at_drill" %in% names(meta)) {
    run_betadisper(d_hel, grouping = meta$ds_at_drill, label, dist_lab = "Hellinger+Euclid")
    run_betadisper(d_ait, meta$ds_at_drill, label, "RobustAitchison")
    run_betadisper(d_jac, meta$ds_at_drill, label, "Jaccard")
    # and same for LOG_aFWD
  }
}

# ============================
# Run for your two datasets
# ============================
# Get reference levels from META1
ds_lv   <- levels(factor(META1$ds_at_drill))
size_lv <- if ("size" %in% names(META1) && is.factor(META1$size)) levels(META1$size) else NULL

# Add rowname and harmonize types
meta_mid_only_row <- meta_mid_only %>%
  mutate(rowname = paste0(natman, "::MID")) %>%
  mutate(
    ds_at_drill = factor(as.character(ds_at_drill), levels = ds_lv),
    dw_size     = if (!is.null(size_lv)) {
      factor(as.character(dw_size), levels = size_lv)
    } else {
      suppressWarnings(as.numeric(dw_size))
    }
  )

stopifnot(all(rownames(otu_mid_only) == meta_mid_only_row$rowname))
# Run the suite
run_beta_suite(otu_mid_only, meta_mid_only_row, label = "MID_only")

# LOG + aFWD (triplet-only, combined)
# Ensure rownames exactly "<natman>::LOG_aFWD"
rownames(otu_log_afwd) <- paste0(gsub("::.*$", "", rownames(otu_log_afwd)), "::LOG_aFWD")
meta_log_afwd_row <- meta_log_afwd %>% mutate(rowname = paste0(natman, "::LOG_aFWD"))
stopifnot(all(rownames(otu_log_afwd) == meta_log_afwd_row$rowname))
run_beta_suite(otu_log_afwd, meta_log_afwd_row, label = "LOG_aFWD")

