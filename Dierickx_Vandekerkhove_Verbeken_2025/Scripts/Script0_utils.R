# ===================================================
# Script: Script0_utils.R
# Project: Deadwood fungal community analyses
# Purpose: Shared utility functions, palettes, and reproducibility helpers
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# FAIR notes:
# - Central utility script reused by downstream analysis scripts
# - Defines shared colour palettes and helper functions
# - Assumes downstream scripts provide analysis-specific objects where needed
# - Recommended to source first: source("Scripts/Script0_utils.R")
#
# Core dependencies: tidyverse, janitor
# Optional downstream dependencies used by some functions:
# vegan, ggplot2, iNEXT, permute, lme4, glue, rstatix
# ===================================================
required_pkgs <- c("tidyverse", "janitor")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
})
# Reproducibility settings used across scripts
SEED_GLOBAL <- 42
set.seed(SEED_GLOBAL)
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion")
# --------------------------
# Color palettes (kept here so all downstream scripts can reuse them)
# --------------------------
dw_colors <- c( LOG = "#8B5C2D", FWD = "#E2B600", SNAG = "#59A14F", SOIL = "#4E79A7", aFWD = "#FFE066", fFWD = "#B89400", FRUITBODY = "#8E44AD", POSITIVE = "#C39BD3")
ds_colors <- c("0"="#66A61E","1"="#1B9E77","2"="#7570B3","3"="#E7298A","4"="#D95F02","5"="#A6761D")
DS_colors <- c(LIVING="#66A61E", EARLY="#1B9E77", AVERAGE="#E7298A", LATE="#A6761D")
group_colors <- c(WOODY="#A0522D", SOIL="#2E86AB")
depth_colors <- c(INNER="#0072B2", OUTER="#E69F00", MIXED="#CC79A7")
position_colors <- c(BASE="#009E73", MIDDLE="#D55E00", UPPER="#F0E442", ENDO="lightgreen")
sc_colors <- c("VERY_FINE"= "#5D4037","FINE"= "#BCAAA4", "SMALL"="#56B4E9", "MEDIUM"="#E15759", "LARGE"="#76B7B2", "VERY_LARGE"="#B07AA1")
microhab_base_cols <- c( BASE = position_colors[["BASE"]], MIDDLE = position_colors[["MIDDLE"]], UPPER = position_colors[["UPPER"]],ATTACHED = dw_colors[["aFWD"]], FALLEN = dw_colors[["fFWD"]], SNAG  = dw_colors[["SNAG"]])
aspect_colors <- c(NORTH = "navy", SOUTH = "orange")
# ---- Minimal FAIR / reproducibility helpers ----------------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

assert_objects <- function(obj_names, env = parent.frame()) {
  missing <- obj_names[!vapply(obj_names, exists, logical(1), envir = env, inherits = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required object(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

write_session_info <- function(out_file = "sessionInfo.txt") {
  readr::write_lines(capture.output(sessionInfo()), out_file)
  invisible(out_file)
}

# --------------------------
# Per-sample abundance thresholding (set counts < min_reads to 0)
# --------------------------
threshold_by_sample <- function(otu, min_reads = 10, drop_empty = TRUE) {
  stopifnot(is.matrix(otu) || is.data.frame(otu))
  x <- as.matrix(otu)
  storage.mode(x) <- "double"
  x <- x * (x >= min_reads)  # zero anything below the threshold
  if (drop_empty) {
    keep <- colSums(x) > 0
    x <- x[, keep, drop = FALSE]
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}
#---------------------------------------------------
# Helper to make sample rarefaction plots per factor
make_rarefaction_plot <- function(fac_name, pal_vals, title_lab, file_stub, input_meta) {
  assert_objects(c("otu_work", "threshold"))
  meta <- input_meta %>% filter(!is.na(.data[[fac_name]]))
  samp <- meta$sample
  otu  <- otu_work[samp, , drop = FALSE]
  # Remove zero-read samples
  sample_sums <- rowSums(otu)
  if(any(sample_sums == 0)) {
    cat("Removing", sum(sample_sums == 0), "zero-read samples\n")
    otu <- otu[sample_sums > 0, , drop = FALSE]
    valid_samples <- rownames(otu)
    meta <- meta %>% filter(sample %in% valid_samples)
  } else {
    valid_samples <- rownames(otu)
  }
  pdf(NULL)
  rc <- vegan::rarecurve(otu, step = 1000, label = FALSE, return = "list")
  dev.off()
  rc_df <- purrr::imap_dfr(
    rc,
    ~tibble(
      sample = valid_samples[.y],  # Use position index
      reads = attr(.x, "Subsample"),
      richness = as.integer(.x)
    )
  ) %>%
    left_join(meta %>% select(sample, !!sym(fac_name)), by = "sample")
  p <- ggplot(rc_df, aes(x = reads, y = richness, group = sample, color = .data[[fac_name]])) +
    geom_line(alpha = 0.7) +
    scale_color_manual(values = pal_vals, na.value = "grey70") +
    theme_minimal(base_size = 13) +
    labs(title = title_lab, x = "Reads (subsampled)", y = "Observed SH richness", color = fac_name)
  print(p)
  ensure_dir("plots/LOG")
  ggsave(paste0("plots/LOG/rarefaction_LOG_by_", file_stub, "_", threshold, ".png"), 
         p, width = 8, height = 5, dpi = 300)
}

# ---- iNEXT builders for LOG contrasts (rarefy to 250k before incidence when specified)
build_inext_for_factor <- function(meta_df, lvl_col, rarefy_to = 250000) {
  assert_objects("otu_work")
  samp <- meta_df$sample
  otu  <- otu_work[samp, , drop = FALSE]
  # Optional read-depth filter before rarefaction (only keep samples with enough reads)
  reads <- rowSums(otu)
  keep  <- reads >= rarefy_to
  if (sum(keep) < 2) return(NULL)
  otu <- otu[keep, , drop = FALSE]
  meta_df <- meta_df %>% filter(sample %in% rownames(otu))
  if (nrow(meta_df) < 2) return(NULL)
  otu_ra <- vegan::rrarefy(otu, sample = rarefy_to)
  
  lvls <- sort(unique(meta_df[[lvl_col]]))
  inc_list <- lapply(lvls, function(v) {
    s <- meta_df %>% filter(.data[[lvl_col]] == v) %>% pull(sample)
    if (length(s) < 2) return(NULL)
    ot <- otu_ra[s, , drop = FALSE]; ot <- ot[, colSums(ot) > 0, drop = FALSE]
    c(length(s), colSums(ot > 0))
  })
  names(inc_list) <- lvls
  inc_list[!sapply(inc_list, is.null)]
}
plot_inext_factor <- function(inc_list, pal_vals, title_lab, file_stub, legend_lab) {
  if (is.null(inc_list) || length(inc_list) == 0) return(invisible(NULL))
  endpoints <- sapply(inc_list, function(x) min(ceiling(x[1]*2.5), 175))
  out <- lapply(names(inc_list), function(k) {
    o <- iNEXT(inc_list[[k]], q = c(0,1,2), datatype = "incidence_freq", endpoint = endpoints[k])
    as_tibble(o$iNextEst$size_based) %>% mutate(group = k)
  })
  df <- bind_rows(out)
  obs_n <- sapply(inc_list, function(x) x[1])
  df <- df %>% mutate(interp_extrap = ifelse(t <= obs_n[group], "interpolated", "extrapolated"))
  interp_pts <- df %>% group_by(group, Order.q) %>% filter(t == max(t[interp_extrap == "interpolated"])) %>% ungroup()
  p <- ggplot(df, aes(x = t, y = qD, color = group, fill = group, linetype = interp_extrap)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
    geom_point(data = interp_pts, aes(x = t, y = qD, color = group),
               inherit.aes = FALSE, size = 3, shape = 21, fill = "white", stroke = 1.2) +
    scale_color_manual(values = pal_vals) +
    scale_fill_manual(values = pal_vals) +
    scale_linetype_manual(values = c(interpolated="solid", extrapolated="dotted")) +
    facet_wrap(~ Order.q, labeller = as_labeller(c("0"="q = 0\nRichness\n(Number of SHs)",
                                                   "1"="q = 1\nShannon diversity\n(Number of typical SHs)",
                                                   "2"="q = 2\nSimpson diversity\n(Number of dominant SHs)")),
               ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 13) +
    labs(title = title_lab, x = "Number of samples (observed & extrapolated)",
         y = "Estimated diversity (qD)", color = legend_lab, fill = legend_lab, linetype = "Type")
  ggsave(paste0("collectors_iNEXT_by_", file_stub, ".png"),
         p, width = 8, height = 12, dpi = 300)
  print(p)
}
# -----------------------------------------------------------------------------
# Per-sample alpha diversity from an OTU matrix
alpha_from_otu <- function(X) {
  reads <- rowSums(X)
  H     <- vegan::diversity(X, "shannon")
  tibble(
    sample   = rownames(X),
    reads    = reads,
    log_reads= log1p(reads),
    richness = vegan::specnumber(X),
    H        = H,
    q1       = exp(H),
    q2       = vegan::diversity(X, "invsimpson")  # Hill q2
  )}

gm_mean <- function(x, pseudo = 1) exp(mean(log(x + pseudo))) - pseudo
if (!exists("combine_counts")) {
  combine_counts <- function(X, idx, method = c("geo", "sum", "mean", "median"), pseudo = 1) {
    method <- match.arg(method)
    M <- as.matrix(X[idx, , drop = FALSE])
    if (nrow(M) == 0) {stop("`idx` selects no rows.")    }
    if (method == "sum") {return(colSums(M))}
    if (method == "mean") {return(colMeans(M))}
    if (method == "median") {return(apply(M, 2, median))}
    Mp <- M + pseudo
    exp(colMeans(log(Mp))) - pseudo
  }}
mode_or_median <- function(x) {
  if (is.numeric(x)) return(median(x, na.rm = TRUE))
  ux <- na.omit(as.character(x)); if (!length(ux)) return(NA_character_)
  tab <- sort(table(ux), decreasing = TRUE)
  names(tab)[1]
}
clean_nat <- function(x) gsub('"', "", x)
# Overdispersion for GLMM
overdisp_phi <- function(mod) {
  rp  <- residuals(mod, type = "pearson")
  X2  <- sum(rp^2, na.rm = TRUE)
  df  <- df.residual(mod)
  if (!is.finite(X2) || !is.finite(df) || df <= 0) {
    return(list(phi = NA_real_, ci = c(NA_real_, NA_real_), df = df, p_over = NA_real_))
  }
  phi <- X2 / df
  # 95% CI
  ci  <- c(qchisq(0.025, df) / df, qchisq(0.975, df) / df)
  # One-sided p
  p_over <- 1 - pchisq(X2, df)
  list(phi = phi, ci = ci, df = df, p_over = p_over)
}
boot_diff_median <- function(tbl, nboot = 10000, seed = 42) {
  set.seed(seed)
  replicate(nboot, {
    wi <- sample(tbl$dist[tbl$type=="within"],  replace=TRUE)
    bt <- sample(tbl$dist[tbl$type=="between"], replace=TRUE)
    median(bt) - median(wi)
  })
}
boot_natman_effects <- function(D, grp, B=2000, seed=42) {
  set.seed(seed)
  M <- as.matrix(D)
  upper_idx <- which(upper.tri(M), arr.ind = TRUE)
  same <- grp[upper_idx[,1]] == grp[upper_idx[,2]]
  d_within  <- M[cbind(upper_idx[,1], upper_idx[,2])][ same]
  d_between <- M[cbind(upper_idx[,1], upper_idx[,2])][!same]
  est <- median(d_between) - median(d_within)
  boot_diffs <- replicate(B, {
    wi <- sample(d_within,  replace=TRUE)
    bt <- sample(d_between, replace=TRUE)
    median(bt) - median(wi)
  })
  ci <- quantile(boot_diffs, c(0.025, 0.975))
  
  tibble(
    diff_median = est,
    diff_CI_lo  = ci[1],
    diff_CI_hi  = ci[2],
    n_within    = length(d_within),
    n_between   = length(d_between)
  )
}
has_levels <- function(x) is.factor(x) && nlevels(droplevels(x)) > 1

fix_dist_labels <- function(D, samps) {
  stopifnot(inherits(D, "dist"))
  n  <- attr(D, "Size")
  if (length(samps) != n)
    stop("Length of `samps` does not match distance size")
  
  labs <- attr(D, "Labels")
  
  # Case 1: no labels at all
  if (is.null(labs)) {
    attr(D, "Labels") <- samps
    return(D)
  }
  
  # Case 2: labels are a permutation of samps
  if (setequal(labs, samps)) {
    if (identical(labs, samps)) return(D)
    M <- as.matrix(D)
    M <- M[samps, samps, drop = FALSE]
    return(as.dist(M))
  }
  
  # Case 3: labels look like 1..n, map positionally
  if (all(grepl("^[0-9]+$", labs)) && setequal(as.integer(labs), seq_len(n))) {
    attr(D, "Labels") <- samps[as.integer(labs)]
    return(D)
  }
  
  stop("Distance labels and expected sample IDs disagree in a non recoverable way")
}
perm_free_or <- function(nperm = 999, strata = NULL) {
  if (exists("perm_free", mode = "function")) {
    return(perm_free(nperm = nperm, strata = strata))
  } else {
    ctl <- permute::how(nperm = nperm)
    if (!is.null(strata)) {
      permute::setBlocks(ctl) <- strata
    }
    return(ctl)
  }
}
d_ait_or_hel <- function(X) {
  if (inherits(try(vegdist(X, method = "robust.aitchison"), silent=TRUE), "try-error")) {
    return(dist(decostand(X, "hellinger"), method = "euclidean"))
  } else {
    return(vegdist(X, method = "robust.aitchison"))
  }
}
clr_counts_or <- function(X, pseudo = 1) {
  if (exists("clr_counts", mode = "function")) return(clr_counts(X))
  Xp <- X + pseudo
  gm <- exp(rowMeans(log(Xp)))
  sweep(log(Xp) - log(gm), 2, 0, "+")
}
mk_within_between_tbl <- function(DM, groups) {
  U <- which(upper.tri(DM), arr.ind = TRUE)
  tibble(
    type = ifelse(groups[U[,1]] == groups[U[,2]], "within","between"),
    dist = DM[cbind(U[,1], U[,2])]
  )
}

wilcoxon_auc_once <- function(Dmat, lab, upper_idx) {
  same <- lab[upper_idx[,1]] == lab[upper_idx[,2]]
  dU   <- Dmat[cbind(upper_idx[,1], upper_idx[,2])]
  r <- rank(c(dU[!same], dU[same])); n1 <- sum(!same); n2 <- sum(same)
  U  <- sum(r[seq_len(n1)]) - n1*(n1+1)/2
  U/(n1*n2)
}
ICC_from_lmer <- function(mod){
  v <- as.data.frame(VarCorr(mod))
  vb <- v$vcov[v$grp=="natman"]; vw <- v$vcov[v$grp=="Residual"]
  vb / (vb + vw)
}
plot_nmds <- function(nmds_obj, meta, title = "NMDS", col_var) {
  df <- scores(nmds_obj, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(meta, by = "sample")
  stress_label <- paste0("Stress = ", round(nmds_obj$stress, 3))
  ggplot(df, aes(x = NMDS1, y = NMDS2, color = .data[[col_var]])) +
    geom_point(size = 3, alpha = 0.8) +
    theme_minimal(base_size = 14) +
    labs(title = title) +
    annotate("text",
             x = min(df$NMDS1), y = max(df$NMDS2),
             label = stress_label, hjust = 0, vjust = 1, size = 5, fontface = "italic")
}

check_glmm_lme4 <- function(m, name = "model", tol_grad = 0.002) {
  cat("\n[Convergence check:", name, "]\n")
  
  ## Messages from optimiser
  msg <- m@optinfo$conv$lme4$messages
  if (!is.null(msg)) {
    warning("Convergence message for ", name, ": ",
            paste(msg, collapse = " | "))
    cat("Convergence messages:\n", paste(msg, collapse = "\n"), "\n")
  } else {
    cat("No optimiser convergence messages.\n")
  }
  
  ## Gradient size
  grad <- m@optinfo$derivs$gradient
  if (!is.null(grad)) {
    max_grad <- max(abs(grad))
    cat("max|grad| =", signif(max_grad, 3), "\n")
    if (is.finite(max_grad) && max_grad > tol_grad) {
      warning("Large max|grad| for ", name,
              " estimates may be unreliable.")
    }
  }
  
  ## Singularity of random effects
  sing <- lme4::isSingular(m, tol = 1e-4)
  cat("Singular fit:", sing, "\n")
  invisible(list(messages = msg, max_grad = max_grad %||% NA_real_, singular = sing))
}

plot_resid_glmm <- function(m, out_path) {
  df_res <- data.frame(
    fitted = fitted(m, type = "response"),
    resid  = residuals(m, type = "pearson")
  )
  p_res <- ggplot2::ggplot(df_res, ggplot2::aes(x = fitted, y = resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::labs(
      x = "Fitted (response)",
      y = "Pearson residuals",
      title = "Alpha richness GLMM residuals vs fitted"
    ) +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(out_path, p_res, width = 6.2, height = 4, dpi = 300)
  p_res
}
#######################################################################
## --- 0) Helpers Script 3-------------------------------------------------------
summarise_alltaxa_by <- function(long_df, group_var) {
  gv <- rlang::ensym(group_var)
  long_df %>%
    dplyr::filter(!is.na(!!gv)) %>%
    dplyr::group_by(!!gv, sh_code, species, genus) %>%
    dplyr::summarise(
      mean_rel_abund_sample = mean(rel_abund_sample, na.rm = TRUE),
      prevalence            = dplyr::n_distinct(sample),
      total_reads           = sum(abundance),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      unite_link = paste0(
        "https://unite.ut.ee/bl_forw_sh.php?sh_name=", sh_code, "#fndtn-panel1"
      )
    ) %>%
    dplyr::group_by(!!gv) %>%
    dplyr::arrange(dplyr::desc(mean_rel_abund_sample), .by_group = TRUE)
}

write_per_level <- function(tab, group_var, out_dir, stem) {
  assert_objects("threshold")
  ensure_dir(out_dir)
  lvls <- unique(dplyr::pull(tab, {{ group_var }}))
  for (lv in lvls) {
    sub <- dplyr::filter(tab, {{ group_var }} == lv)
    if (nrow(sub) == 0) next
    fn <- file.path(out_dir, paste0(stem, "_", as.character(lv), "_threshold", threshold, ".tsv"))
    readr::write_tsv(sub, fn)
  }
}
sets_from_summary <- function(df, group_var) {
  stopifnot(all(c("sh_code", group_var) %in% names(df)))
  sp <- split(df$sh_code, df[[group_var]])
  lapply(sp, function(v) unique(v[!is.na(v) & nzchar(v)]))
}
pairwise_jaccard <- function(sets) {
  labs <- names(sets)
  mat  <- matrix(NA_real_, nrow = length(sets), ncol = length(sets), dimnames = list(labs, labs))
  for (i in seq_along(sets)) for (j in seq_along(sets)) {
    a <- sets[[i]]; b <- sets[[j]]
    inter <- length(intersect(a, b)); uni <- length(union(a, b))
    mat[i, j] <- if (uni == 0) NA_real_ else inter / uni
  }
  as.data.frame(mat) %>% rownames_to_column("group")
}
parse_taxon <- function(x) {
  tibble(
    sh_code = sub("\\|.*$", "", x),
    taxon   = sub("^[^|]*\\|", "", x)
  )
}
make_region_levels <- function(dw_levels) {
  unlist(lapply(seq_along(dw_levels), function(k) {
    combn(dw_levels, k, FUN = function(v) paste(v, collapse = "&"))
  }))
}
region_counts_from_sets <- function(sets, region_levels) {
  setnames    <- names(sets)
  all_species <- unique(unlist(sets))
  if (length(all_species) == 0L) {
    out <- integer(length(region_levels))
    names(out) <- region_levels
    return(out)
  }
  mat <- sapply(setnames, function(s) all_species %in% sets[[s]])
  labels <- apply(mat, 1, function(z) {
    if (!any(z)) {
      return(NA_character_)
    }
    paste(setnames[z], collapse = "&")
  })
  fac <- factor(labels, levels = region_levels)
  counts <- tabulate(fac, nbins = length(region_levels))
  names(counts) <- region_levels
  counts
}
make_combo <- function(x) paste(sort(x), collapse = "&")
###############################################################
compare_iNEXT_df_at_coverage <- function(iNEXT_df, group_col = "dw_type2",
                                         orders = c(0,1,2), target_cov = 0.95,
                                         allow_extrapolation = TRUE) {
  df <- iNEXT_df %>%
    filter(Order.q %in% orders) %>%
    # SE from 95% CI width
    mutate(se = (qD.UCL - qD.LCL) / (2 * 1.96)) %>%
    filter(!is.na(SC), !is.na(qD))
  
  if (!allow_extrapolation && "interp_extrap" %in% names(df)) {
    df <- df %>% filter(interp_extrap == "interpolated")
  }
  
  # common coverage per order (cap at target_cov)
  Cstar_tbl <- df %>%
    group_by(Order.q, !!rlang::sym(group_col)) %>%
    summarise(maxSC = max(SC, na.rm = TRUE), .groups = "drop") %>%
    group_by(Order.q) %>%
    summarise(Cstar = min(maxSC, na.rm = TRUE), .groups = "drop")
  
  # pick nearest point to C* (with t retained)
  pts <- df %>%
    inner_join(Cstar_tbl, by = "Order.q") %>%
    group_by(Order.q, !!rlang::sym(group_col)) %>%
    slice_min(abs(SC - Cstar), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(Order.q,
              group = .data[[group_col]],
              t, SC, qD, se, Cstar)
  
  # require ≥ 2 groups per order
  valid_orders <- pts %>% count(Order.q) %>% filter(n >= 2) %>% pull(Order.q)
  pts <- pts %>% filter(Order.q %in% valid_orders)
  
  # pairwise z-tests + ratios
  contrasts <- map_dfr(split(pts, pts$Order.q), function(d) {
    cmb <- combn(seq_len(nrow(d)), 2)
    map_dfr(seq_len(ncol(cmb)), function(k) {
      i <- d[cmb[1, k], ]; j <- d[cmb[2, k], ]
      diff <- i$qD - j$qD
      se_diff <- sqrt(i$se^2 + j$se^2)
      z <- diff / se_diff
      p <- 2 * pnorm(-abs(z))
      ratio <- i$qD / j$qD
      se_log <- sqrt((i$se / i$qD)^2 + (j$se / j$qD)^2)
      tibble(
        Order.q = i$Order.q, Cstar = i$Cstar,
        group1 = i$group, group2 = j$group,
        qD1 = i$qD, qD2 = j$qD,
        diff = diff, z = z, p = p,
        ratio = ratio,
        ratio_LCL = exp(log(ratio) - 1.96 * se_log),
        ratio_UCL = exp(log(ratio) + 1.96 * se_log)
      )
    })
  }) %>%
    group_by(Order.q) %>%
    mutate(p_adj = p.adjust(p, "BH")) %>%
    ungroup()
  
  list(points = pts, contrasts = contrasts, Cstar = Cstar_tbl)
}
######################## LOG ##################################################
build_intra_meta_ds_at_drill <- function(scope,
                                         min_samples_per_tree = 3,
                                         min_position_levels = 2) {
  if (!identical(scope, "LOG")) {
    stop("This redefined builder currently supports only scope = 'LOG'.")
  }
  
  keep <- META1 %>%
    dplyr::filter(
      dw_type2 == "LOG",
      !is.na(sample),
      !is.na(natman),
      !is.na(ds_at_drill),
      ds_at_drill != "0",
      !is.na(position_2)
    ) %>%
    dplyr::mutate(
      size = droplevels(factor(size)),
      Position = dplyr::case_when(
        position_2 %in% c("BASE", "MIDDLE", "UPPER") ~ as.character(position_2),
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Position)) %>%
    dplyr::group_by(natman) %>%
    dplyr::filter(
      dplyr::n() >= min_samples_per_tree,
      dplyr::n_distinct(Position) >= min_position_levels
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      natman      = droplevels(factor(natman)),
      ds_at_drill = droplevels(factor(ds_at_drill)),
      Position    = droplevels(factor(Position, levels = c("BASE", "MIDDLE", "UPPER")))
    ) %>%
    droplevels()
  
  cat("Retained trees:", dplyr::n_distinct(keep$natman), "\n")
  cat("Retained samples:", nrow(keep), "\n")
  cat("Position levels retained:\n")
  print(table(keep$Position, useNA = "ifany"))
  cat("Samples per natman:\n")
  print(table(keep$natman))
  
  keep
}
auc_from_group <- function(D, grp, upper_idx = NULL) {
  M <- as.matrix(D)
  if (is.null(upper_idx)) upper_idx <- which(upper.tri(M), arr.ind = TRUE)
  same <- grp[upper_idx[,1]] == grp[upper_idx[,2]]
  dU <- M[cbind(upper_idx[,1], upper_idx[,2])]
  d_between <- dU[!same]; d_within <- dU[same]
  if (length(d_between) == 0 || length(d_within) == 0) return(c(AUC=NA_real_, r_rb=NA_real_))
  r <- rank(c(d_between, d_within)); n1 <- length(d_between); n2 <- length(d_within)
  U <- sum(r[seq_len(n1)]) - n1*(n1+1)/2; A <- U/(n1*n2); c(AUC = A, r_rb = 2*A - 1)
}
get_natman_R2 <- function(fit) { df <- as.data.frame(fit); if (!"natman" %in% rownames(df)) return(NA_real_); df["natman","R2"] }
pairwise_groups <- function(D, grp, upper_idx) {
  M <- as.matrix(D); dU <- M[cbind(upper_idx[,1], upper_idx[,2])]
  same <- grp[upper_idx[,1]] == grp[upper_idx[,2]]
  list(within = dU[same], between = dU[!same])
}
#################### SNAG #############################
rescale_to_target <- function(X, target = 1e4) {
  s <- rowSums(X); s[s == 0] <- 1
  sweep(X, 1, target / s, "*")
}

boot_diff_ci <- function(x_between, x_within,
                         B = get0("B_BOOT", ifnotfound = 2000),
                         seed = SEED_GLOBAL) {
  set.seed(seed)
  bd <- replicate(
    B,
    median(sample(x_between, replace = TRUE)) -
      median(sample(x_within, replace = TRUE))
  )
  tibble(
    est_diff_median = median(x_between) - median(x_within),
    CI_low  = unname(quantile(bd, 0.025)),
    CI_high = unname(quantile(bd, 0.975))
  )
}
cliffs_delta <- function(x_between, x_within) {
  nx <- length(x_between); ny <- length(x_within)
  r  <- rank(c(x_between, x_within))
  rx <- sum(r[seq_len(nx)])
  Ux <- rx - nx * (nx + 1) / 2
  (2 * Ux) / (nx * ny) - 1
}
#################################### FWD ################################
`%||%` <- function(x, y) if (!is.null(x)) x else y
safe_scale01 <- function(x) (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
as_pa <- function(M){ (M>0) * 1 }
hellinger <- function(M){ decostand(M, "hellinger") }
auroc_from_ranks <- function(within_vec, between_vec){
  n1 <- length(within_vec); n2 <- length(between_vec)
  ranks <- rank(c(within_vec, between_vec))
  U <- sum(ranks[seq_len(n1)]) - n1*(n1+1)/2
  U/(n1*n2)
}
rank_biserial <- function(x_within, x_between){
  U <- wilcox.test(x_between, x_within, exact = FALSE)$statistic[[1]]
  n1 <- length(x_between); n2 <- length(x_within)
  (2*U)/(n1*n2) - 1
}
adonis2_adjR2 <- function(adonis2_row_R2, n, p){ 1 - (1 - adonis2_row_R2) * ((n - 1)/(n - p - 1)) }
get_R2_table <- function(ad2, DM){
  n <- attr(DM, "Size") %||% nrow(as.matrix(DM))
  p <- sum(ad2$Df[!is.na(ad2$Df)])
  ad2 %>% as.data.frame() %>% rownames_to_column("term") %>%
    dplyr::filter(term != "Total") %>% mutate(R2adj = adonis2_adjR2(R2, n, p))
}
plot_ecdf <- function(df, title){
  ggplot(df, aes(x=dist, colour=type)) + stat_ecdf() +
    labs(x = "Pairwise distance", y = "ECDF", title = title, colour = "Pair type") +
    theme_bw() 
}
############## SCOPE3 ###################
get_R2_term <- function(fit, term) {
  df <- as.data.frame(fit)
  if (!term %in% rownames(df)) return(NA_real_)
  df[term, "R2"]}
## collapse pair table to simple "within vs between" for a given factor
pairs_by <- function(pairs, factor = c("state_global", "state_within_tree", "natman")) {
  factor <- match.arg(factor)
  
  if (factor == "state_global") {
    out <- pairs |>
      mutate(type = if_else(st_i == st_j, "within", "cross")) |>
      dplyr::select(type, dist)
  } else if (factor == "state_within_tree") {
    out <- pairs |>
      filter(nat_i == nat_j) |>
      mutate(type = if_else(st_i == st_j, "within", "cross")) |>
      dplyr::select(type, dist)
  } else {
    out <- pairs |>
      mutate(type = if_else(nat_i == nat_j, "within", "between")) |>
      dplyr::select(type, dist)
  }
  droplevels(out)
}

## AUC (common language effect size) and rank biserial r
auc_once <- function(df) {
  d1 <- df$dist[df$type %in% c("cross", "between")]
  d0 <- df$dist[df$type == "within"]
  r  <- rank(c(d1, d0))
  n1 <- length(d1)
  n0 <- length(d0)
  U  <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  A  <- U / (n1 * n0)
  tibble(AUC = A, r_rb = 2 * A - 1)
}
## nonparametric CI for difference in medians (group1 minus within)
median_diff_boot <- function(df, B = get0("cfg", ifnotfound = list(B_boot = 2000))$B_boot) {
  d1 <- df$dist[df$type %in% c("cross", "between")]
  d0 <- df$dist[df$type == "within"]
  est <- median(d1) - median(d0)
  boot <- replicate(
    B,
    median(sample(d1, replace = TRUE)) - median(sample(d0, replace = TRUE))
  )
  ci <- quantile(boot, c(0.025, 0.975), names = FALSE)
  tibble(
    diff_median = est,
    diff_CI_lo  = ci[1],
    diff_CI_hi  = ci[2],
    n_within    = length(d0),
    n_between   = length(d1)
  )
}

## design preserving state shuffle within natman
shuffle_state_within_tree <- function(state, natman) {
  unlist(
    tapply(
      state, natman,
      function(x) sample(x, length(x)),
      simplify = FALSE
    ),
    use.names = FALSE
  )
}

## state nulls for R2 in blocked PERMANOVA on paired subset
state_R2_null_blocked <- function(D, meta, add_log_reads,
                                  B = get0("cfg", ifnotfound = list(B_R2 = 999))$B_R2) {
  nperm_val <- get0("cfg", ifnotfound = list(nperm = 999))$nperm
  perm_blk <- permute::how(nperm = nperm_val)
  permute::setBlocks(perm_blk) <- meta$natman
  
  replicate(B, {
    meta_null <- meta |>
      mutate(fwd_null = shuffle_state_within_tree(fwd_state, natman))
    if (add_log_reads) {
      as.data.frame(
        adonis2(
          D ~ log_reads + umi + fwd_null,
          data = meta_null,
          by   = "terms",
          permutations = perm_blk
        )
      )["fwd_null", "R2"]
    } else {
      as.data.frame(
        adonis2(
          D ~ umi + fwd_null,
          data = meta_null,
          by   = "terms",
          permutations = perm_blk
        )
      )["fwd_null", "R2"]
    }
  })
}

## simple summarizer for observed vs null R2
summ_R2_null <- function(R2_obs, R2_null_vec) {
  tibble(
    R2_obs        = R2_obs,
    R2_null_mean  = mean(R2_null_vec),
    R2_null_CI_lo = quantile(R2_null_vec, 0.025),
    R2_null_CI_hi = quantile(R2_null_vec, 0.975),
    p_emp         = mean(R2_null_vec >= R2_obs),
    SES           = (R2_obs - mean(R2_null_vec)) / sd(R2_null_vec),
    perc          = mean(R2_null_vec <= R2_obs)
  )
}
## helper for AUC from grouping (within vs between) on global distances
upper_idx_global <- function(D) {
  which(upper.tri(as.matrix(D)), arr.ind = TRUE)
}
mk_pair_tbl <- function(D, meta) {
  M    <- as.matrix(D)
  labs <- rownames(M)
  nat  <- meta$natman[match(labs, meta$sample)]
  st   <- meta$fwd_state[match(labs, meta$sample)]
  U    <- which(upper.tri(M), arr.ind = TRUE)
  
  tibble(
    i    = labs[U[, 1]],
    j    = labs[U[, 2]],
    dist = M[cbind(U[, 1], U[, 2])],
    nat_i = nat[U[, 1]], nat_j = nat[U[, 2]],
    st_i  = st[U[, 1]],  st_j  = st[U[, 2]]
  )
}
######IN vs OUT ####
pick_col <- function(df, candidates) {
  nm <- names(df)
  hit <- intersect(candidates, nm)
  if (!length(hit)) stop("Missing required column(s): ", paste(candidates, collapse = " OR "))
  hit[1]}

ensure_numeric_matrix <- function(X) {
  X <- as.matrix(X)
  X <- apply(X, 2, function(z) suppressWarnings(as.numeric(z)))
  storage.mode(X) <- "numeric"
  X[is.na(X)] <- 0
  X}
axis_contrib <- function(ord, D) {
  Y <- scores(ord, display = "sites")
  labs <- rownames(Y)
  Dm <- as.matrix(D)
  Dsub <- Dm[labs, labs, drop = FALSE]
  Dvec <- as.vector(as.dist(Dsub))
  vals <- vapply(seq_len(ncol(Y)), function(j) {
    d1 <- dist(Y[, j, drop = FALSE])
    stats::cor(Dvec, as.vector(d1), use = "pairwise.complete.obs")^2
  }, numeric(1))
  frac <- vals / sum(vals, na.rm = TRUE)
  tibble(axis = colnames(Y), R2 = vals, frac = frac)}

subset_dist <- function(d, keep_labels) {
  M <- as.matrix(d)
  labs <- rownames(M)
  if (is.null(labs)) labs <- attr(d, "Labels")
  idx <- labs %in% keep_labels
  as.dist(M[idx, idx, drop = FALSE])}

shannon1 <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  v[is.na(v)] <- 0
  if (sum(v) > 0) vegan::diversity(v, index = "shannon") else NA_real_}

circle_df <- function(center = c(0, 0), r = 1.6, n = 360) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = center[1] + r * cos(t), y = center[2] + r * sin(t))}

make_venn_plot <- function(left_only, overlap, right_only,
                           labels = c("INNER", "OUTER"),
                           subtitle = NULL) {
  cL <- circle_df(c(-1, 0))
  cR <- circle_df(c( 1, 0))
  ggplot() +
    geom_path(data = cL, aes(x, y)) +
    geom_path(data = cR, aes(x, y)) +
    annotate("text", x = -2.1, y = 1.5, label = labels[1], hjust = 0, size = 4) +
    annotate("text", x =  2.1, y = 1.5, label = labels[2], hjust = 1, size = 4) +
    annotate("text", x = -1.6, y = 0,   label = left_only,  size = 5) +
    annotate("text", x =  0,   y = 0,   label = overlap,    size = 5) +
    annotate("text", x =  1.6, y = 0,   label = right_only, size = 5) +
    coord_equal(xlim = c(-3, 3), ylim = c(-2.2, 2.2), expand = FALSE) +
    labs(title = "INNER vs OUTER species overlap", subtitle = subtitle) +
    theme_void(base_size = 12)}

fmt_mean_sd <- function(x) {
  sprintf("%.1f ± %.1f", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
}
safe_str <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

run_seglen_tests <- function(df, label,
                             decay_col = "ds_fac",
                             seglen_col = "seg_len",
                             outdir = get0("OUTDIR_TAB", ifnotfound = ".")) {
  suppressPackageStartupMessages(library(rstatix))
  
  stopifnot(all(c(decay_col, seglen_col) %in% names(df)))
  
  d0 <- df %>%
    transmute(
      ds_fac_raw = .data[[decay_col]],
      seg_len    = .data[[seglen_col]]
    ) %>%
    filter(!is.na(ds_fac_raw), !is.na(seg_len), is.finite(seg_len)) %>%
    mutate(
      ds_fac_raw = droplevels(factor(ds_fac_raw)),
      ds_num     = suppressWarnings(as.numeric(as.character(ds_fac_raw))))
  
  cat("\n====================\n", label, "segment-length tests\n====================\n")
  cat("N pairs:", nrow(d0), "\n")
  cat("Decay stages (levels):", paste(levels(d0$ds_fac_raw), collapse = ", "), "\n\n")
  
  kw <- kruskal.test(seg_len ~ ds_fac_raw, data = d0)
  cat("Kruskal–Wallis (seg_len ~ decay stage)\n"); print(kw); cat("\n")
  
  dunn <- d0 %>%
    rstatix::dunn_test(seg_len ~ ds_fac_raw, p.adjust.method = "BH") %>%
    arrange(p.adj)
  cat("Dunn post hoc (BH-adjusted)\n"); print(dunn); cat("\n")
  
  if (all(!is.na(d0$ds_num))) {
    sp <- suppressWarnings(cor.test(d0$ds_num, d0$seg_len, method = "spearman", exact = FALSE))
    cat("Spearman trend (seg_len ~ numeric decay stage)\n"); print(sp); cat("\n")
  } else {
    sp <- NULL
    cat("Spearman trend skipped: ds_fac_raw not coercible to numeric reliably.\n\n")
  }
  
  bf <- d0 %>% rstatix::levene_test(seg_len ~ ds_fac_raw, center = "median")
  cat("Brown–Forsythe / Levene (median-centered): heterogeneity of spread\n"); print(bf); cat("\n")
  
  summ <- d0 %>%
    group_by(ds_fac_raw) %>%
    summarise(
      n      = n(),
      mean   = mean(seg_len),
      sd     = sd(seg_len),
      median = median(seg_len),
      iqr    = IQR(seg_len),
      .groups = "drop"
    )
  cat("Decay-stage summaries\n"); print(summ); cat("\n")
  ensure_dir(outdir)
  write_csv2(dunn, file.path(outdir, paste0(label, "_seglen_Dunn_BH.csv")))
  write_csv2(summ, file.path(outdir, paste0(label, "_seglen_stage_summaries.csv")))
  write_csv2(as.data.frame(bf), file.path(outdir, paste0(label, "_seglen_BrownForsythe.csv")))
  
  invisible(list(data = d0, kw = kw, dunn = dunn, spearman = sp, brown_forsythe = bf, summaries = summ))
}
loo_permanova <- function(X, meta, group_id = "natman",
                          formula_rhs = "log_reads + umi + aspect + ds_at_drill",
                          dist_method = "robust.aitchison",
                          nperm_outer = NPERM_LOO,
                          blocked_by = NULL, by_TYPE="margin",
                          keep_terms = c("aspect","ds_at_drill"),
                          out_csv = NULL) {
  stopifnot(nrow(X) == nrow(meta), identical(rownames(X), meta$sample))
  groups <- (factor(meta[[group_id]]))
  trees  <- levels(groups)
  
  out <- purrr::map_dfr(trees, function(t) {
    keep <- groups != t & !is.na(groups)
    if (sum(keep) < 5) return(NULL)
    
    # recompute distance on the kept rows to avoid reuse bias
    Dk <- vegdist(X[keep, , drop = FALSE], method = dist_method)
    
    data_k <- (meta[keep, , drop = FALSE])
    
    permk <- how(nperm = nperm_outer)
    if (!is.null(blocked_by)) setBlocks(permk) <- (data_k[[blocked_by]])
    
    fml <- reformulate(termlabels = all.vars(stats::terms(stats::as.formula(paste("~", formula_rhs)))),
                       response = NULL)
    # adonis2: pass D + selected columns via data=; RHS uses "."
    a2 <- adonis2(Dk ~ ., data = data_k[, all.vars(stats::terms(fml)), drop = FALSE],
                  permutations = permk, by = by_TYPE)
    
    as_tibble(as.data.frame(a2), rownames = "term") |>
      filter(term %in% keep_terms) |>
      transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  
  if (nrow(out)) {
    summ <- out |>
      group_by(term) |>
      summarise(median_R2 = median(R2, na.rm = TRUE),
                p05 = quantile(p, 0.05, na.rm = TRUE),
                p95 = quantile(p, 0.95, na.rm = TRUE),
                .groups = "drop")
    print(summ)
    if (!is.null(out_csv)) {
      readr::write_csv2(out, out_csv)
      readr::write_csv2(summ, sub("\\.csv$", "_summary.csv", out_csv))
    }
  }
  invisible(out)
}

