# ---------------------------------------------------
#  SCRIPT 0 UTILS.R 
# ---------------------------------------------------
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
clean_nat <- function(x) gsub('"', "", x)

# ---------------------------------------------------
#  SCRIPT 3
# ---------------------------------------------------
# Helper to make sample rarefaction plots for an arbitrary factor
make_rarefaction_plot <- function(fac_name, pal_vals, title_lab, file_stub, input_meta) {
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
  ggsave(paste0("plots/LOG/rarefaction_LOG_by_", file_stub, "_", threshold, ".png"), 
         p, width = 8, height = 5, dpi = 300)
}

# ---- iNEXT builders for LOG contrasts (rarefy to 250k before incidence when specified)
build_inext_for_factor <- function(meta_df, lvl_col, rarefy_to = 250000) {
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
  ggsave(paste0("plots/LOG/collectors_iNEXT_LOG_by_", file_stub, ".png"),
         p, width = 8, height = 12, dpi = 300)
  print(p)
}
# ---------------------------------------------------
#  SCRIPT 5 
# ---------------------------------------------------
# 1) Helpers  
# -----------------------------------------------------------------------------
# Per-sample alpha covariates from an OTU matrix
alpha_from_otu <- function(otu) {
  tibble(
    sample   = rownames(otu),
    richness = vegan::specnumber(otu),
    reads    = rowSums(otu),
    log_reads= log1p(reads)
  )
}

# Overdispersion for GLMM
overdisp_phi <- function(mod) {
  rp  <- residuals(mod, type = "pearson")
  X2  <- sum(rp^2, na.rm = TRUE)
  df  <- df.residual(mod)
  if (!is.finite(X2) || !is.finite(df) || df <= 0) {
    return(list(phi = NA_real_, ci = c(NA_real_, NA_real_), df = df, p_over = NA_real_))
  }
  phi <- X2 / df
  # 95% CI from chi-square quantiles on X2
  ci  <- c(qchisq(0.025, df) / df, qchisq(0.975, df) / df)
  # One-sided p-value for overdispersion (Pr[X2 >= observed] under H0: φ=1)
  p_over <- 1 - pchisq(X2, df)
  list(phi = phi, ci = ci, df = df, p_over = p_over)
}
# Euler plot helper with robust fallback to bars
make_euler_plot <- function(R2_all, R2_tree, R2_micro, R2_tree_u, R2_mic_u,
                            title_main, title_sub, outfile) {
  # sanitize
  R2_all    <- pmin(pmax(R2_all, 0), 1)
  R2_tree_u <- pmax(R2_tree_u, 0)
  R2_mic_u  <- pmax(R2_mic_u, 0)
  R2_shared <- pmax(R2_tree + R2_micro - R2_all, 0)
  
  pct_tree_u  <- 100 * R2_tree_u
  pct_micro_u <- 100 * R2_mic_u
  pct_shared  <- 100 * R2_shared
  pct_unex    <- 100 * (1 - R2_all)
  
  cat("\n PERCENTAGES IN EULER \n TOTAL: ", R2_all, "\n R2_tree_unconstrained: ",
      pct_tree_u, "% \n R2_microhab_unconstrained: ", pct_micro_u,
      "% \n R2_shared_unconstrained: ", pct_shared,
      "% \n R2_unexplained_unconstrained: ", pct_unex )
  
  fit <- try(eulerr::euler(c(
    "Tree"              = pct_tree_u + pct_shared,
    "Microhabitat"      = pct_micro_u + pct_shared,
    "Tree&Microhabitat" = pct_shared
  ), shape = "ellipse"), silent = TRUE)
  
  fallback_bar <- function() {
    dfb <- tibble(
      component = factor(c("Tree unique","Shared","Microhab unique","Unexplained"),
                         levels = c("Tree unique","Shared","Microhab unique","Unexplained")),
      percent = c(pct_tree_u, pct_shared, pct_micro_u, pct_unex)
    )
    p <- ggplot(dfb, aes(x = "Partition", y = percent, fill = component)) +
      geom_col(width = 0.6, color = "grey30") +
      coord_flip() +
      scale_fill_manual(values = c(
        "Tree unique"      = "#1b9e77",
        "Microhab unique"  = "#d95f02",
        "Shared"           = "grey60",
        "Unexplained"      = "grey85"
      )) +
      labs(x = NULL, y = "Percent of total variation",
           title = title_main, subtitle = paste0(title_sub,
                                                 sprintf("\nTree %.1f%% | Microhab %.1f%% | Shared %.1f%% | Unexplained %.1f%%",
                                                         pct_tree_u, pct_micro_u, pct_shared, pct_unex))) +
      theme_minimal(base_size = 12) +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
    ggsave(outfile, p, width = 9, height = 4, dpi = 300)
    print(p)
  }
  
  if (inherits(fit, "try-error") || is.null(fit$ellipses)) { fallback_bar(); return(invisible(NULL)) }
  
  ell <- try({
    purrr::imap_dfr(fit$ellipses, function(e, nm) {
      e <- as.list(e)
      a <- if (!is.null(e$a)) e$a else if (!is.null(e$r1)) e$r1 else if (!is.null(e$r)) e$r else NA_real_
      b <- if (!is.null(e$b)) e$b else if (!is.null(e$r2)) e$r2 else if (!is.null(e$r)) e$r else NA_real_
      tibble(h = e$h, k = e$k, a = a, b = b, phi = e$phi, set = nm, area = pi * a * b)
    })
  }, silent = TRUE)
  
  if (inherits(ell, "try-error") || any(!is.finite(ell$a)) || any(!is.finite(ell$b))) { fallback_bar(); return(invisible(NULL)) }
  
  # polygon builder
  ellipse_points <- function(h, k, a, b, phi, n = 400) {
    t <- seq(0, 2*pi, length.out = n); xp <- a * cos(t); yp <- b * sin(t)
    tibble(x = h + xp * cos(phi) - yp * sin(phi),
           y = k + xp * sin(phi) + yp * cos(phi))
  }
  poly_df <- purrr::pmap_dfr(
    list(ell$h, ell$k, ell$a, ell$b, ell$phi, ell$set),
    ~ ellipse_points(..1, ..2, ..3, ..4, ..5) %>% mutate(set = ..6)
  )
  centers <- ell %>% select(set, h, k, a, b) %>% pivot_wider(names_from = set, values_from = c(h,k,a,b))
  
  vx <- centers$h_Microhabitat - centers$h_Tree
  vy <- centers$k_Microhabitat - centers$k_Tree
  vlen <- sqrt(vx^2 + vy^2); ux <- ifelse(vlen > 0, vx/vlen, 1); uy <- ifelse(vlen > 0, vy/vlen, 0)
  
  lab_df <- tibble(
    region = c("Tree unique", "Shared", "Microhab unique"),
    x = c(centers$h_Tree - ux * (0.6 * centers$a_Tree),
          (centers$h_Tree + centers$h_Microhabitat)/2,
          centers$h_Microhabitat + ux * (0.5 * centers$a_Microhabitat)),
    y = c(centers$k_Tree - uy * (0.6 * centers$a_Tree),
          (centers$k_Tree + centers$k_Microhabitat)/2,
          centers$k_Microhabitat + uy * (0.6 * centers$a_Microhabitat)),
    value = c(pct_tree_u, pct_shared, pct_micro_u)
  ) %>%
    mutate(text = case_when(
      region == "Tree unique"     ~ sprintf("\n\nTree \nnatman|microhab\n %.1f%%", value),
      region == "Shared"          ~ sprintf("Shared\nnatman∩microhab\n%.1f%%", value),
      TRUE                        ~ sprintf("\n\nMicrohabitat\nmicrohab|natman\n %.1f%%", value)
    ))
  
  rng <- poly_df %>% summarise(xmin = min(x), xmax = max(x), ymin = min(y), ymax = max(y))
  unexp_label <- sprintf("Unexplained (1 − R²_all): %.1f%%", pct_unex)
  
  p <- ggplot(poly_df, aes(x, y)) +
    geom_polygon(aes(group = set, fill = set), alpha = 0.5, colour = "grey30", linewidth = 0.3) +
    geom_text(data = lab_df, aes(x, y, label = text), lineheight = 0.95, size = 3.6, fontface = "bold") +
    annotate("label", x = rng$xmax, y = rng$ymin, hjust = 1.5, vjust = 0,
             label = unexp_label, size = 3.3, label.r = grid::unit(0.2, "lines"), alpha = 0.9) +
    scale_fill_manual(values = c(Tree = "#1b9e77", Microhabitat = "#d95f02")) +
    coord_equal() + theme_void() +
    ggtitle(title_main, subtitle = title_sub) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "none")
  print(p)
  ggsave(outfile, p, width = 9, height = 6, dpi = 300)
}

# Indicator-friendly rescaling to common library size
rescale_to_target <- function(otu) {
  ts <- max(2000, min(rowSums(otu)))  # keep ≥2000 reads
  sweep(otu, 1, rowSums(otu), "/") * ts
}

# -----------------------------------------------------------------------------
# 2) Build scope-specific metadata (LOG / +aFWD / +fFWD / +SNAG)
# -----------------------------------------------------------------------------
build_intra_meta_ds_at_drill <- function(scope) {
  base <- META1 %>%
      filter(rq =="INTRA") %>% 
      mutate(size = droplevels(factor(size)))
  
  # what to include per scope
  add_flags <- list(
    LOG                      = list(LOG=TRUE, aFWD=FALSE, fFWD=FALSE, SNAG=FALSE),
    LOG_aFWD                 = list(LOG=TRUE, aFWD=TRUE , fFWD=FALSE, SNAG=FALSE),
    LOG_aFWD_fFWD            = list(LOG=TRUE, aFWD=TRUE , fFWD=TRUE , SNAG=FALSE),
    LOG_aFWD_fFWD_SNAG       = list(LOG=TRUE, aFWD=TRUE , fFWD=TRUE , SNAG=TRUE )
  )[[scope]]
  if (is.null(add_flags)) stop("Unknown scope: ", scope)

  keep <- base %>%
    filter(
      !is.na(natman),
      !is.na(ds_at_drill), ds_at_drill != "0"        # drop fresh wood at drill
    ) %>%
    mutate(
      # Position across substrates present in the scope
      Position = case_when(
        dw_type == "LOG" & position %in% c("BASE","MIDDLE","UPPER") ~ position,
        dw_type == "FWD" & position %in% c("AFWD1","AFWD2")         ~ "ATTACHED",
        dw_type == "FWD" & position == "ATTACHED"                   ~ "ATTACHED",
        dw_type == "FWD" & position == "FALLEN"                     ~ "FALLEN",  
        dw_type == "SNAG"                                           ~ "SNAG",
        TRUE ~ NA_character_
      )
    )
  cat("FWD samples by position:\n")
  print(keep %>% filter(dw_type == "FWD") %>% count(position, Position))
  
  cat("All dw_type values:\n")
  print(unique(keep$dw_type))
  
  # filter to the scope
  keep <- keep %>%
    filter(
      (add_flags$LOG  & dw_type == "LOG")  |
        (add_flags$aFWD & dw_type == "FWD"  & Position == "ATTACHED") |
        (add_flags$fFWD & dw_type == "FWD"  & Position == "FALLEN")   |
        (add_flags$SNAG & dw_type == "SNAG")
    ) %>%
    # drop unusable positions
    filter(!position %in% c("CROWN"),
           !(dw_type == "FWD" & as.character(position) == "CROWN")) %>%
    mutate(
      natman      = droplevels(factor(natman)),
      ds_at_drill = droplevels(factor(ds_at_drill)),
      Position   = droplevels(factor(Position, levels = c("BASE","MIDDLE","UPPER","ATTACHED","FALLEN","SNAG")))
    ) %>%
    filter(!is.na(Position)) %>%
    droplevels()
  print(unique(keep$Position))
  keep
}
build_intra_meta <- function(scope) {
  base <- META1 %>%
    filter(rq =="INTRA") %>% 
    mutate(size = droplevels(factor(size)))
  
  # what to include per scope
  add_flags <- list(
    LOG                      = list(LOG=TRUE, aFWD=FALSE, fFWD=FALSE, SNAG=FALSE),
    LOG_aFWD                 = list(LOG=TRUE, aFWD=TRUE , fFWD=FALSE, SNAG=FALSE),
    LOG_aFWD_fFWD            = list(LOG=TRUE, aFWD=TRUE , fFWD=TRUE , SNAG=FALSE),
    LOG_aFWD_fFWD_SNAG       = list(LOG=TRUE, aFWD=TRUE , fFWD=TRUE , SNAG=TRUE )
  )[[scope]]
  if (is.null(add_flags)) stop("Unknown scope: ", scope)
  
  keep <- base %>%
    filter(
      !is.na(natman),
      !is.na(decay_stage)       
    ) %>%
    mutate(
      # Position across substrates present in the scope
      Position = case_when(
        dw_type == "LOG" & position %in% c("BASE","MIDDLE","UPPER") ~ position,
        dw_type == "FWD" & position %in% c("AFWD1","AFWD2")         ~ "ATTACHED",
        dw_type == "FWD" & position == "ATTACHED"                   ~ "ATTACHED",
        dw_type == "FWD" & position == "FALLEN"                     ~ "FALLEN",  
        dw_type == "SNAG"                                           ~ "SNAG",
        TRUE ~ NA_character_
      )
    )

  # filter to the scope
  keep <- keep %>%
    filter(
      (add_flags$LOG  & dw_type == "LOG")  |
        (add_flags$aFWD & dw_type == "FWD"  & Position == "ATTACHED") |
        (add_flags$fFWD & dw_type == "FWD"  & Position == "FALLEN")   |
        (add_flags$SNAG & dw_type == "SNAG")
    ) %>%
    # drop unusable positions
    filter(!position %in% c("CROWN"),
           !(dw_type == "FWD" & as.character(position) == "CROWN")) %>%
    mutate(
      natman      = droplevels(factor(natman)),
      decay_stage = droplevels(factor(decay_stage)),
      Position   = droplevels(factor(Position, levels = c("BASE","MIDDLE","UPPER","ATTACHED","FALLEN","SNAG")))
    ) %>%
    filter(!is.na(Position)) %>%
    droplevels()
  keep
}

# -----------------------------------------------------------------------------
# 3) One scope runner
# -----------------------------------------------------------------------------
run_scope <- function(scope) {
  cat("\n ----- starting for: ", scope, "----- \n ------  Have fun!  ------- \n")
  outdir <- file.path("plots/INTRA", scope)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta <- build_intra_meta(scope)
  if (nrow(meta) < 3) { message("Skip ", scope, ": <3 samples"); return(invisible(NULL)) }
  
  # align OTU
  otu <- otu_matrix_filt[meta$sample, , drop = FALSE]
  stopifnot(identical(rownames(otu), meta$sample))
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  
  # add read-depth covariate (for controlling sequencing depth)
  alpha_df  <- alpha_from_otu(otu)
  meta_ext  <- meta %>% left_join(alpha_df %>% select(sample, reads, log_reads), by = "sample")
  
  # colors present in this scope
  mlevels    <- levels(meta_ext$Position)
  micro_cols <- microhab_base_cols[names(microhab_base_cols) %in% mlevels]
  micro_cols <- micro_cols[match(mlevels, names(micro_cols))]
  
  # ---------- Alpha diversity (plots for ALL factors) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- ALPHA plots (Position, DS, Size, Diameter smooth) ------- \n")
  
  # Position
  p_alpha_pos <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = Position, y = richness, fill = Position)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = micro_cols) +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by microhabitat — ", scope),
         x = "Microhabitat", y = "Observed SH richness", fill = "Microhabitat")
  ggsave(file.path(outdir, paste0("alpha_richness_by_microhab_", scope, "_", threshold, ".png")),
         p_alpha_pos, width = 9, height = 4.8, dpi = 300)
  print(p_alpha_pos)
  
  # decay_stage
  p_alpha_ds <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = decay_stage, y = richness, fill = decay_stage)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = DS_colors) +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by decay at drill — ", scope),
         x = "decay_stage", y = "Observed SH richness", fill = "decay_stage")
  ggsave(file.path(outdir, paste0("alpha_richness_by_ds_", scope, ".png")),
         p_alpha_ds, width = 8, height = 4.5, dpi = 300)
  print(p_alpha_ds)
  
  # size
  p_alpha_sc <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = size, y = richness, fill = size)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = sc_colors, na.value = "yellow4") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by size class — ", scope),
         x = "size", y = "Observed SH richness", fill = "size")
  ggsave(file.path(outdir, paste0("alpha_richness_by_sizeclass_", scope, ".png")),
         p_alpha_sc, width = 8, height = 4.5, dpi = 300)
  print(p_alpha_sc)
  # smooth vs diameter
  p_alpha_diam <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = diameter_at_drill, y = richness)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = TRUE) +
    geom_vline(xintercept = c(1,39,69,99), linetype = "dashed", linewidth = 0.4, color = "yellow4") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness vs diameter — ", scope),
         x = "Diameter at drill (mm)", y = "Observed SH richness")
  ggsave(file.path(outdir, paste0("alpha_richness_vs_diameter_", scope, ".png")),
         p_alpha_diam, width = 9, height = 4, dpi = 300)
  print(p_alpha_diam)
  
  # ---------- Distances ----------
  D_robAit      <- vegdist(otu, method = "robust.aitchison")
  comm_hell     <- decostand(otu, method = "hellinger")
  D_euclid_hell <- dist(comm_hell, method = "euclidean")
  comm_rel      <- decostand(otu, method = "total")
  D_bray_rel    <- vegdist(comm_rel, method = "bray")
  D_bray_relS   <- vegdist(sqrt(comm_rel), method = "bray")
  
  # ---------- PERMANOVA (unblocked + blocked by natman) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- PERMANOVAs (control read depth via log_reads; use size) ------- \n")

    perm_unblocked <- how(nperm = 999)
  perma_unblocked <- list(
    robAit  = adonis2(D_robAit ~ log_reads + natman + Position + decay_stage + size,
                      data = meta_ext, by = "margin", permutations = perm_unblocked),
    hellEuc = adonis2(D_euclid_hell ~ log_reads + natman + Position + decay_stage + size,
                      data = meta_ext, by = "margin", permutations = perm_unblocked),
    relBray = adonis2(D_bray_rel ~ log_reads + natman + Position + decay_stage + size,
                      data = meta_ext, by = "margin", permutations = perm_unblocked)
  )
  capture.output(perma_unblocked, file = file.path(outdir, paste0("permanova_unblocked_", scope, ".txt")))
  cat("\npermanova UN_blocked:\n"); print(perma_unblocked)
  
  # Readout natman R2 across distances
  natman_R2 <- purrr::map_dfr(perma_unblocked, ~{
    df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)
  }, .id = "distance") %>%
    dplyr::filter(term == "natman") %>%
    dplyr::select(distance, R2, `Pr(>F)`)
  cat("\nUnblocked PERMANOVA — natman R2 by distance:\n")
  print(natman_R2)
  
  # Blocked within tree
  perm_blocked <- how(nperm = 999); setBlocks(perm_blocked) <- meta_ext$natman
  perma_blocked <- adonis2(
    D_robAit ~ log_reads + Position + decay_stage + size,
    data = meta_ext, by = "margin", permutations = perm_blocked
  )
  capture.output(perma_blocked, file = file.path(outdir, paste0("permanova_blocked_", scope, ".txt")))
  cat("\npermanova BLOCKED (within-tree):\n"); print(perma_blocked)
  
  # ---------- LOO stability (within-tree PERMANOVA, blocked) ----------
  cat("\n ----- LOO stability (within-tree PERMANOVA) ------- \n")
  trees <- meta_ext$natman %>% as.character() %>% unique() %>% setdiff(NA)
  loo <- purrr::map_dfr(trees, function(t) {
    keep <- meta_ext$natman != t & !is.na(meta_ext$natman)
    if (sum(keep) < 6) return(NULL)
    Dk <- vegdist(otu[keep, , drop = FALSE], method = "robust.aitchison")
    permk <- how(nperm = 499); setBlocks(permk) <- droplevels(factor(meta_ext$natman[keep]))
    m <- adonis2(Dk ~ log_reads + Position + decay_stage + size,
                 data = meta_ext[keep, ], by = "margin", permutations = permk)
    mm <- as.data.frame(m); mm$term <- rownames(mm)
    as_tibble(mm) %>%
      dplyr::filter(term %in% c("Position","decay_stage","size")) %>%
      dplyr::transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  if (nrow(loo) > 0) {
    loo_sum <- loo %>% group_by(term) %>%
      summarise(median_R2 = median(R2, na.rm = TRUE),
                p05 = quantile(p, 0.05, na.rm = TRUE),
                p95 = quantile(p, 0.95, na.rm = TRUE),
                .groups = "drop")
    cat("\nLeave-one-tree-out stability (within-tree PERMANOVA):\n")
    print(loo_sum)
    #readr::write_csv(loo,     file.path(outdir, paste0("loo_within_tree_", scope, ".csv")))
    readr::write_csv(loo_sum, file.path(outdir, paste0("loo_within_tree_summary_", scope, ".csv")))
  }
  
  # ---------- Partial dbRDA + Euler (Condition log_reads) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- VARIANCE PARTITIONING (Condition log_reads) ------ \n")
  
  mod_all   <- capscale(otu ~ natman + Position + decay_stage + size + Condition(log_reads),
                        distance = "robust.aitchison", data = meta_ext)
  mod_tree  <- capscale(otu ~ natman + Condition(log_reads),
                        distance = "robust.aitchison", data = meta_ext)
  mod_micro <- capscale(otu ~ Position + decay_stage + size + Condition(log_reads),
                        distance = "robust.aitchison", data = meta_ext)
  mod_tree_c  <- capscale(otu ~ natman + Condition(Position + decay_stage + size + log_reads),
                          distance = "robust.aitchison", data = meta_ext)
  mod_micro_c <- capscale(otu ~ Position + decay_stage + size + Condition(natman + log_reads),
                          distance = "robust.aitchison", data = meta_ext)
  
  cat("\n report adj.R2 values!")
  micro_adjR2 <- RsquareAdj(mod_micro_c)$adj.r.squared
  cat("\nPartial dbRDA — microhabitat | tree + log(reads): adj.R2 =", round(micro_adjR2, 4), "\n")
  anova_micro_terms <- anova.cca(mod_micro_c, by = "terms", permutations = perm_blocked)
  cat("\nTerm-wise tests (blocked by natman & reads):\n"); print(anova_micro_terms)
  tree_adjR2 <- RsquareAdj(mod_tree_c)$adj.r.squared
  cat("\nPartial dbRDA — tree | microhabitat: adj.R2 =", round(tree_adjR2, 4), "\n")
  # Significance (cannot block by natman here — degenerate)
  anova_tree <- anova.cca(mod_tree_c, permutations = perm_unblocked)
  cat("\nOverall test for tree | microhabitat (unblocked):\n"); print(anova_tree)
  
  make_euler_plot(
    R2_all   = RsquareAdj(mod_all)$r.squared,
    R2_tree  = RsquareAdj(mod_tree)$r.squared,
    R2_micro = RsquareAdj(mod_micro)$r.squared,
    R2_tree_u= RsquareAdj(mod_tree_c)$r.squared,
    R2_mic_u = RsquareAdj(mod_micro_c)$r.squared,
    title_main = paste0("Variance partitioning (robust Aitchison dbRDA) — ", scope),
    title_sub  = "Tree vs Microhab complex (Position + ds + size | log_reads)",
    outfile    = file.path(outdir, paste0("euler_variance_partitioning_", scope, ".png"))
  )
  
  # ---------- Dispersion checks ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- DISPERSION CHECKS ------ \n")
  disp_nat  <- anova(betadisper(D_robAit, meta_ext$natman))
  disp_mh   <- anova(betadisper(D_robAit, meta_ext$Position))
  disp_size <- anova(betadisper(D_robAit, droplevels(meta_ext$size)))
  sink(file.path(outdir, paste0("permdisp_", scope, ".txt"))); print(disp_nat); print(disp_mh); print(disp_size); sink()
  cat("\nDispersion natman:\n");  print(disp_nat)
  cat("\nDispersion Position:\n");print(disp_mh)
  cat("\nDispersion size:\n");print(disp_size)
  
  # ---------- Indicator species ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- INDICATORS ------ \n")
  if (nlevels(meta_ext$Position) >= 2) {
    otu_sc <- rescale_to_target(otu)
    cn <- tibble(sh_code = colnames(otu_sc)) %>%
      left_join(tax %>% dplyr::select(sh_code, lowest_taxon), by = "sh_code") %>%
      mutate(new_name = ifelse(!is.na(lowest_taxon), paste0(sh_code, "|", lowest_taxon), sh_code)) %>%
      pull(new_name)
    colnames(otu_sc) <- cn
    otu_sc <- otu_sc[, colSums(otu_sc) > 0, drop = FALSE]
    
    g_micro <- droplevels(meta_ext$Position)
    g_ds    <- droplevels(meta_ext$decay_stage)
    g_tree  <- droplevels(meta_ext$natman)
    g_size <- droplevels(meta_ext$size)
    
#    ind_mh <- indicspecies::multipatt(otu_sc, g_micro, func = "r.g", control = how(nperm = 999), duleg = TRUE)
#    ind_ds <- indicspecies::multipatt(otu_sc, g_ds,    func = "r.g", control = how(nperm = 999), duleg = TRUE)
#    ind_tr <- indicspecies::multipatt(otu_sc, g_tree,  func = "r.g", control = how(nperm = 999), duleg = TRUE)
#    sink(file.path(outdir, paste0("indicators_corr_", scope, ".txt"))); summary(ind_mh); summary(ind_ds); summary(ind_tr);  sink()
    
    ind_mh_c <- indicspecies::multipatt(otu_sc, g_micro, control = how(nperm = 999), duleg = TRUE)
    ind_ds_c <- indicspecies::multipatt(otu_sc, g_ds,    control = how(nperm = 999), duleg = TRUE)
    ind_tr_c <- indicspecies::multipatt(otu_sc, g_tree,  control = how(nperm = 999), duleg = TRUE)
    sink(file.path(outdir, paste0("indicators_classic_", scope, ".txt"))); summary(ind_mh_c); summary(ind_ds_c); summary(ind_tr_c); sink()
    
#    cat("\nIndicator (Position, corr):\n"); print(summary(ind_mh, alpha = 0.001))
#    cat("\nIndicator (DS, corr):\n");      print(summary(ind_ds, alpha = 0.001))
#    cat("\nIndicator (tree, corr):\n"); print(summary(ind_tr, alpha = 0.001))
    
#    cat("\nIndicator (Position, classic):\n"); print(summary(ind_mh_c, alpha = 0.001))
    cat("\nIndicator (DS, classic):\n");      print(summary(ind_ds_c, alpha = 0.001))
    cat("\nIndicator (tree, classic):\n"); print(summary(ind_tr_c, alpha = 0.001))
    
  }
  
  # ---------- Alpha richness GLMM (with selection, IRR, pub tables) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- ALPHA RICHNESS  ------ \n")
  df <- alpha_df %>%
    left_join(meta_ext %>% dplyr::select(sample, natman, Position, decay_stage, diameter_at_drill), by = "sample") %>%
    mutate(
      natman     = droplevels(natman),
      ds_ord     = ordered(decay_stage, levels = sort(unique(as.character(decay_stage)))),
      diameter_z = as.numeric(scale(diameter_at_drill))
    ) %>%
    filter(!is.na(richness), reads > 0)

  if (nlevels(df$Position) >= 2) {
    ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    form <- richness ~ Position + ds_ord + diameter_z + log_reads + (1|natman)
    
    m_pois <- suppressWarnings(glmer(form, data = df, family = poisson(link = "log"), control = ctrl))
    phi_p  <- overdisp_phi(m_pois)
    use_nb <- is.finite(phi_p$phi) && (phi_p$phi > 1.2)
    if (use_nb) {
      m_nb  <- suppressWarnings(glmer.nb(form, data = df, control = ctrl))
      best  <- if (AIC(m_pois) + 2 < AIC(m_nb)) m_pois else m_nb
      best_type <- if (identical(best, m_pois)) "Poisson (glmer)" else "NB (glmer.nb)"
    } else {
      best <- m_pois; best_type <- "Poisson (glmer)"
    }
    cat(sprintf("Selected alpha model: %s | Overdispersion (Poisson) φ = %.2f | AIC(best)=%.1f\n",
                best_type, phi_p$phi, AIC(best)))
    print(summary(best))
    
    an_tab <- car::Anova(best, type = 2)
    r2_tab <- performance::r2(best)
    
    irr_tab <- broom.mixed::tidy(best, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(term = dplyr::recode(
        term,
        `(Intercept)`      = "Intercept",
        `PositionMIDDLE`   = "Position: MIDDLE vs BASE",
        `PositionUPPER`    = "Position: UPPER vs BASE",
        `PositionATTACHED` = "Position: ATTACHED vs BASE",
        `PositionFALLEN`   = "Position: FALLEN vs BASE",
        `PositionSNAG`     = "Position: SNAG vs BASE",
        `ds_ord.L`         = "Decay (linear)",
        `ds_ord.Q`         = "Decay (quadratic)",
        `ds_ord.C`         = "Decay (cubic)",
        `ds_ord^4`         = "Decay (quartic)",
        `ds_ord^5`         = "Decay (quintic)",
        `diameter_z`       = "Diameter (z)",
        `log_reads`        = "log(reads+1)"
      )) %>%
      dplyr::rename(IRR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value) %>%
      arrange(term)
    
    # EMMeans: Position + ds
    emm_pos <- emmeans(best, ~ Position, type = "response")
    emm_ds  <- emmeans(best, ~ ds_ord,  type = "response")
    
    p_emm_pos <- summary(emm_pos) %>% as.data.frame() %>%
      ggplot(aes(x = Position, y = response, ymin = asymp.LCL, ymax = asymp.UCL)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      labs(x = "Microhabitat", y = "Estimated richness",
           title = paste0("Richness ~ microhabitat (INTRA) — ", scope)) +
      theme_classic(base_size = 11)
    ggsave(file.path(outdir, paste0("alpha_richness_emm_Position_", scope, ".png")),
           p_emm_pos, width = 6, height = 3.2, dpi = 300)
    print(p_emm_pos)
    
    ds_df <- summary(emm_ds) %>% 
      as.data.frame() %>% 
      mutate(ds_ord = factor(ds_ord, levels = c("EARLY", "AVERAGE", "LATE")),
        ds_num = as.numeric(ds_ord))
    
    p_emm_ds <- ggplot(ds_df, aes(x = ds_num, y = response)) +
      geom_line() + geom_point(size = 3) +
      geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), alpha = 0.15) +
      scale_x_continuous(breaks = ds_df$ds_num, labels = as.character(ds_df$ds_ord)) +
      labs(x = "decay_stage (ordered)", y = "Estimated richness",
           title = paste0("Richness ~ decay (INTRA) — ", scope)) +
      theme_classic(base_size = 11)
    ggsave(file.path(outdir, paste0("alpha_richness_emm_DS_", scope, ".png")),
           p_emm_ds, width = 6.2, height = 3.2, dpi = 300)
    print(p_emm_ds)
    
    # Compact publication tables
    tbl_main <- tibble::tibble(
      Effect  = rownames(an_tab),
      Chisq   = an_tab$`Chisq`,
      df      = an_tab$Df,
      p_value = an_tab$`Pr(>Chisq)`
    ) %>%
      dplyr::filter(Effect %in% c("Position","ds_ord","diameter_z","log_reads")) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(match(Effect, c("Position","ds_ord","diameter_z","log_reads")))
    
    tbl_r2 <- tibble::tibble(
      Metric = c("R2_marginal", "R2_conditional"),
      Value  = c(r2_tab$R2_marginal, r2_tab$R2_conditional)
    )
    
    cat("\n=== Publication tables (objects) ===\n")
    print(list(anova = tbl_main, r2 = tbl_r2, irr = irr_tab))
    #readr::write_csv(tbl_main, file.path(outdir, paste0("alpha_pubtable_anova_", scope, ".csv")))
    #readr::write_csv(tbl_r2,   file.path(outdir, paste0("alpha_pubtable_r2_", scope, ".csv")))
    #readr::write_csv(irr_tab,  file.path(outdir, paste0("alpha_pubtable_irr_", scope, ".csv")))
    
    # Save raw model outputs
    capture.output(list(model = summary(best), anova = an_tab, r2 = r2_tab),
                   file = file.path(outdir, paste0("alpha_glmm_", scope, ".txt")))
  }
  
  cat("\n --------------  DONE FOR: ", scope, "  -------------- \n")
  invisible(NULL)
}
run_scope_ds <- function(scope) {
  cat("\n ----- starting for: ", scope, "----- \n ------  Have fun!  ------- \n")
  outdir <- file.path("plots/INTRA", scope)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  meta <- build_intra_meta_ds_at_drill(scope)
  if (nrow(meta) < 3) { message("Skip ", scope, ": <3 samples"); return(invisible(NULL)) }
  
  # align OTU
  otu <- otu_matrix_filt[meta$sample, , drop = FALSE]
  stopifnot(identical(rownames(otu), meta$sample))
  otu <- otu[, colSums(otu) > 0, drop = FALSE]
  
  # add read-depth covariate (for controlling sequencing depth)
  alpha_df  <- alpha_from_otu(otu)
  meta_ext  <- meta %>% 
    left_join(alpha_df %>% select(sample, reads, log_reads), by = "sample") %>%
    # >>> NEW (depth_2): coerce and order
    mutate(depth_2 = if ("depth_2" %in% names(.)) factor(depth_2, levels = c("INNER","OUTER")) else factor(NA))
  
  # colors present in this scope
  mlevels    <- levels(meta_ext$Position)
  micro_cols <- microhab_base_cols[names(microhab_base_cols) %in% mlevels]
  micro_cols <- micro_cols[match(mlevels, names(microhab_base_cols))]
  # >>> NEW (depth_2): simple palette
  depth2_cols <- c(INNER = "#7f8c8d", OUTER = "#2c3e50")
  
  # ---------- Alpha diversity (plots for ALL factors) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- ALPHA plots (Position, DS, Size, Diameter smooth) ------- \n")
  
  # Position
  p_alpha_pos <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = Position, y = richness, fill = Position)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = micro_cols) +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by microhabitat — ", scope),
         x = "Microhabitat", y = "Observed SH richness", fill = "Microhabitat")
  ggsave(file.path(outdir, paste0("alpha_richness_by_microhab_", scope, "_", threshold, ".png")),
         p_alpha_pos, width = 9, height = 4.8, dpi = 300)
  print(p_alpha_pos)
  
  # >>> NEW (depth_2) — alpha plot
  if ("depth_2" %in% names(meta_ext) && nlevels(droplevels(meta_ext$depth_2)) >= 2) {
    p_alpha_depth <- alpha_df %>%
      left_join(meta_ext, by = "sample") %>%
      ggplot(aes(x = depth_2, y = richness, fill = depth_2)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.65) +
      geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
      scale_fill_manual(values = depth2_cols) +
      theme_minimal(base_size = 12) +
      labs(title = paste0("Alpha richness by depth (INNER vs OUTER) — ", scope),
           x = "Depth", y = "Observed SH richness", fill = "Depth")
    ggsave(file.path(outdir, paste0("alpha_richness_by_depth2_", scope, ".png")),
           p_alpha_depth, width = 6.2, height = 4.2, dpi = 300)
    print(p_alpha_depth)
  }
  
  # ds_at_drill
  p_alpha_ds <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = ds_at_drill, y = richness, fill = ds_at_drill)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = ds_colors) +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by decay at drill — ", scope),
         x = "ds_at_drill", y = "Observed SH richness", fill = "ds_at_drill")
  ggsave(file.path(outdir, paste0("alpha_richness_by_ds_", scope, ".png")),
         p_alpha_ds, width = 8, height = 4.5, dpi = 300)
  print(p_alpha_ds)
  
  # size
  p_alpha_sc <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = size, y = richness, fill = size)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.35) +
    scale_fill_manual(values = sc_colors, na.value = "yellow4") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness by size class — ", scope),
         x = "size", y = "Observed SH richness", fill = "size")
  ggsave(file.path(outdir, paste0("alpha_richness_by_sizeclass_", scope, ".png")),
         p_alpha_sc, width = 8, height = 4.5, dpi = 300)
  print(p_alpha_sc)
  
  # smooth vs diameter
  p_alpha_diam <- alpha_df %>%
    left_join(meta_ext, by = "sample") %>%
    ggplot(aes(x = diameter_at_drill, y = richness)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 5), se = TRUE) +
    geom_vline(xintercept = c(1,39,69,99), linetype = "dashed", linewidth = 0.4, color = "yellow4") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Alpha richness vs diameter — ", scope),
         x = "Diameter at drill (mm)", y = "Observed SH richness")
  ggsave(file.path(outdir, paste0("alpha_richness_vs_diameter_", scope, ".png")),
         p_alpha_diam, width = 9, height = 4, dpi = 300)
  print(p_alpha_diam)
  
  # ---------- Distances ----------
  D_robAit      <- vegdist(otu, method = "robust.aitchison")
  comm_hell     <- decostand(otu, method = "hellinger")
  D_euclid_hell <- dist(comm_hell, method = "euclidean")
  comm_rel      <- decostand(otu, method = "total")
  D_bray_rel    <- vegdist(comm_rel, method = "bray")
  D_bray_relS   <- vegdist(sqrt(comm_rel), method = "bray")
  # Jaccard
  otu_pa        <- (otu > 0) * 1L
  D_jaccard     <- vegdist(otu_pa, method = "jaccard", binary = TRUE)
  
  # ---------- PERMANOVA (unblocked + blocked by natman) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- PERMANOVAs (control read depth via log_reads; use size) ------- \n")
  
  perm_unblocked <- how(nperm = 999)
  perma_unblocked <- list(
    robAit  = adonis2(D_robAit ~ log_reads + natman + Position + depth_2 + ds_at_drill + size,  # >>> NEW (depth_2)
                      data = meta_ext, by = "margin", permutations = perm_unblocked),
    hellEuc = adonis2(D_euclid_hell ~ log_reads + natman + Position + depth_2 + ds_at_drill + size,  # >>> NEW
                      data = meta_ext, by = "margin", permutations = perm_unblocked),
    relBray = adonis2(D_bray_rel ~ log_reads + natman + Position + depth_2 + ds_at_drill + size,     # >>> NEW
                      data = meta_ext, by = "margin", permutations = perm_unblocked),
    jaccard = adonis2(D_jaccard ~ log_reads + natman + Position + depth_2 + ds_at_drill + size,      # >>> NEW
                      data = meta_ext, by = "margin", permutations = perm_unblocked)
  )
  capture.output(perma_unblocked, file = file.path(outdir, paste0("permanova_unblocked_", scope, ".txt")))
  cat("\npermanova UN_blocked:\n"); print(perma_unblocked)
  
  # Readout natman R2 across distances
  natman_R2 <- purrr::map_dfr(perma_unblocked, ~{
    df <- as.data.frame(.x); df$term <- rownames(df); as_tibble(df)
  }, .id = "distance") %>%
    dplyr::filter(term == "natman") %>%
    dplyr::select(distance, R2, `Pr(>F)`)
  cat("\nUnblocked PERMANOVA — natman R2 by distance:\n")
  print(natman_R2)
  
  # Blocked within tree
  perm_blocked <- how(nperm = 999); setBlocks(perm_blocked) <- meta_ext$natman
  perma_blocked <- adonis2(
    D_robAit ~ log_reads + Position + depth_2 + ds_at_drill + size,     # >>> NEW (depth_2)
    data = meta_ext, by = "margin", permutations = perm_blocked
  )
  capture.output(perma_blocked, file = file.path(outdir, paste0("permanova_blocked_", scope, ".txt")))
  cat("\npermanova BLOCKED (within-tree, robust Aitchison):\n"); print(perma_blocked)
  
  # Blocked Jaccard
  perma_blocked_jacc <- adonis2(
    D_jaccard ~ log_reads + Position + depth_2 + ds_at_drill + size,    # >>> NEW (depth_2)
    data = meta_ext, by = "margin", permutations = perm_blocked
  )
  capture.output(perma_blocked_jacc, file = file.path(outdir, paste0("permanova_blocked_jaccard_", scope, ".txt")))
  cat("\npermanova BLOCKED (within-tree, Jaccard P/A):\n"); print(perma_blocked_jacc)
  
  # ---------- LOO stability (within-tree PERMANOVA, blocked) ----------
  cat("\n ----- LOO stability (within-tree PERMANOVA) ------- \n")
  trees <- meta_ext$natman %>% as.character() %>% unique() %>% setdiff(NA)
  loo <- purrr::map_dfr(trees, function(t) {
    keep <- meta_ext$natman != t & !is.na(meta_ext$natman)
    if (sum(keep) < 6) return(NULL)
    Dk <- vegdist(otu[keep, , drop = FALSE], method = "robust.aitchison")
    permk <- how(nperm = 499); setBlocks(permk) <- droplevels(factor(meta_ext$natman[keep]))
    m <- adonis2(Dk ~ log_reads + Position + depth_2 + ds_at_drill + size,   # >>> NEW (depth_2)
                 data = meta_ext[keep, ], by = "margin", permutations = permk)
    mm <- as.data.frame(m); mm$term <- rownames(mm)
    as_tibble(mm) %>%
      dplyr::filter(term %in% c("Position","depth_2","ds_at_drill","size")) %>%   # >>> NEW (depth_2)
      dplyr::transmute(left_out = t, term, R2, p = `Pr(>F)`)
  })
  if (nrow(loo) > 0) {
    loo_sum <- loo %>% group_by(term) %>%
      summarise(median_R2 = median(R2, na.rm = TRUE),
                p05 = quantile(p, 0.05, na.rm = TRUE),
                p95 = quantile(p, 0.95, na.rm = TRUE),
                .groups = "drop")
    cat("\nLeave-one-tree-out stability (within-tree PERMANOVA):\n")
    print(loo_sum)
    readr::write_csv(loo_sum, file.path(outdir, paste0("loo_within_tree_summary_", scope, ".csv")))
  }
  
  # ---------- Partial dbRDA + Euler (Condition log_reads) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- VARIANCE PARTITIONING (Condition log_reads) ------ \n")
  
  mod_all   <- capscale(otu ~ natman + Position + depth_2 + ds_at_drill + size + Condition(log_reads),    # >>> NEW
                        distance = "robust.aitchison", data = meta_ext)
  mod_tree  <- capscale(otu ~ natman + Condition(log_reads),
                        distance = "robust.aitchison", data = meta_ext)
  mod_micro <- capscale(otu ~ Position + depth_2 + ds_at_drill + size + Condition(log_reads),             # >>> NEW
                        distance = "robust.aitchison", data = meta_ext)
  mod_tree_c  <- capscale(otu ~ natman + Condition(Position + depth_2 + ds_at_drill + size + log_reads),  # >>> NEW
                          distance = "robust.aitchison", data = meta_ext)
  mod_micro_c <- capscale(otu ~ Position + depth_2 + ds_at_drill + size + Condition(natman + log_reads),  # >>> NEW
                          distance = "robust.aitchison", data = meta_ext)
  
  cat("\n report adj.R2 values!")
  micro_adjR2 <- RsquareAdj(mod_micro_c)$adj.r.squared
  cat("\nPartial dbRDA — microhabitat | tree + log(reads): adj.R2 =", round(micro_adjR2, 4), "\n")
  anova_micro_terms <- anova.cca(mod_micro_c, by = "terms", permutations = perm_blocked)
  cat("\nTerm-wise tests (blocked by natman & reads):\n"); print(anova_micro_terms)
  tree_adjR2 <- RsquareAdj(mod_tree_c)$adj.r.squared
  cat("\nPartial dbRDA — tree | microhabitat: adj.R2 =", round(tree_adjR2, 4), "\n")
  anova_tree <- anova.cca(mod_tree_c, permutations = perm_unblocked)
  cat("\nOverall test for tree | microhabitat (unblocked):\n"); print(anova_tree)
  
  make_euler_plot(
    R2_all   = RsquareAdj(mod_all)$r.squared,
    R2_tree  = RsquareAdj(mod_tree)$r.squared,
    R2_micro = RsquareAdj(mod_micro)$r.squared,
    R2_tree_u= RsquareAdj(mod_tree_c)$r.squared,
    R2_mic_u = RsquareAdj(mod_micro_c)$r.squared,
    title_main = paste0("Variance partitioning (robust Aitchison dbRDA) — ", scope),
    title_sub  = "Tree vs Microhab complex (Position + depth + ds + size | log_reads)",   # >>> NEW
    outfile    = file.path(outdir, paste0("euler_variance_partitioning_", scope, ".png"))
  )
  
  # ---------- Dispersion checks ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- DISPERSION CHECKS ------ \n")
  disp_nat   <- anova(betadisper(D_robAit, meta_ext$natman))
  disp_mh    <- anova(betadisper(D_robAit, meta_ext$Position))
  disp_size  <- anova(betadisper(D_robAit, droplevels(meta_ext$size)))
  disp_depth <- if (nlevels(droplevels(meta_ext$depth_2)) >= 2) anova(betadisper(D_robAit, droplevels(meta_ext$depth_2))) else NULL  # >>> NEW
  sink(file.path(outdir, paste0("permdisp_", scope, ".txt"))); print(disp_nat); print(disp_mh); print(disp_size); if(!is.null(disp_depth)) print(disp_depth); sink()
  cat("\nDispersion natman:\n");  print(disp_nat)
  cat("\nDispersion Position:\n");print(disp_mh)
  cat("\nDispersion size:\n");   print(disp_size)
  if(!is.null(disp_depth)) { cat("\nDispersion depth_2:\n"); print(disp_depth) }  # >>> NEW
  
  # ---------- Indicator species ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- INDICATORS ------ \n")
  if (nlevels(meta_ext$Position) >= 2) {
    otu_sc <- rescale_to_target(otu)
    cn <- tibble(sh_code = colnames(otu_sc)) %>%
      left_join(tax %>% dplyr::select(sh_code, lowest_taxon), by = "sh_code") %>%
      mutate(new_name = ifelse(!is.na(lowest_taxon), paste0(sh_code, "|", lowest_taxon), sh_code)) %>%
      pull(new_name)
    colnames(otu_sc) <- cn
    otu_sc <- otu_sc[, colSums(otu_sc) > 0, drop = FALSE]
    
    g_micro <- droplevels(meta_ext$Position)
    g_ds    <- droplevels(meta_ext$ds_at_drill)
    g_tree  <- droplevels(meta_ext$natman)
    g_size  <- droplevels(meta_ext$size)
    
    ind_mh_c <- indicspecies::multipatt(otu_sc, g_micro, control = how(nperm = 999), duleg = TRUE)
    ind_ds_c <- indicspecies::multipatt(otu_sc, g_ds,    control = how(nperm = 999), duleg = TRUE)
    ind_tr_c <- indicspecies::multipatt(otu_sc, g_tree,  control = how(nperm = 999), duleg = TRUE)
    sink(file.path(outdir, paste0("indicators_classic_", scope, ".txt"))); summary(ind_mh_c); summary(ind_ds_c); summary(ind_tr_c); sink()
    
    cat("\nIndicator (DS, classic):\n");      print(summary(ind_ds_c, alpha = 0.001))
    cat("\nIndicator (tree, classic):\n");    print(summary(ind_tr_c, alpha = 0.001))
  }
  
  # ---------- Alpha richness GLMM (ordered ds_at_drill: 1->5) ----------
  cat("\n --------------  -------------- \n")
  cat("\n ----- ALPHA RICHNESS  ------ \n")
  
  # Build analysis frame
  df <- alpha_df %>%
    left_join(
      meta_ext %>% dplyr::select(sample, natman, Position, depth_2, ds_at_drill, diameter_at_drill),  # >>> NEW (depth_2)
      by = "sample"
    ) %>%
    mutate(
      natman     = droplevels(natman),
      diameter_z = as.numeric(scale(diameter_at_drill)),
      depth_2    = factor(depth_2, levels = c("INNER","OUTER"))  # >>> NEW
    ) %>%
    filter(!is.na(richness), reads > 0)
  
  # Make ds_at_drill an ordered factor 1..5 (drops unused levels later)
  ds_num <- suppressWarnings(as.integer(as.character(df$ds_at_drill)))
  if (all(is.na(ds_num))) {
    warning("ds_at_drill not numeric-like; using its factor order as-is (ordered factor).")
    df$ds_fac <- factor(df$ds_at_drill, ordered = TRUE)
  } else {
    df$ds_fac <- factor(ds_num, levels = 1:5, ordered = TRUE)
  }
  df$ds_fac <- droplevels(df$ds_fac)
  
  # Build model formula dynamically to avoid 1-level factor issues
  has_depth2 <- "depth_2" %in% names(df) && nlevels(droplevels(df$depth_2)) >= 2
  rhs_terms  <- c("Position", if (has_depth2) "depth_2", "ds_fac", "diameter_z", "log_reads", "(1|natman)")
  form <- as.formula(paste("richness ~", paste(rhs_terms, collapse = " + ")))
  
  if (nlevels(df$Position) >= 2 && nlevels(df$ds_fac) >= 2) {
    ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    m_pois <- suppressWarnings(glmer(form, data = df, family = poisson(link = "log"), control = ctrl))
    phi_p  <- overdisp_phi(m_pois)
    use_nb <- is.finite(phi_p$phi) && (phi_p$phi > 1.2)
    
    if (use_nb) {
      m_nb  <- suppressWarnings(glmer.nb(form, data = df, control = ctrl))
      best  <- if (AIC(m_pois) + 2 < AIC(m_nb)) m_pois else m_nb
      best_type <- if (identical(best, m_pois)) "Poisson (glmer)" else "NB (glmer.nb)"
    } else {
      best <- m_pois; best_type <- "Poisson (glmer)"
    }
    cat(sprintf("Selected alpha model: %s | Overdispersion (Poisson) φ = %.2f | AIC(best)=%.1f\n",
                best_type, phi_p$phi, AIC(best)))
    print(summary(best))
    
    # Type-II tests + R2
    an_tab <- car::Anova(best, type = 2)
    r2_tab <- performance::r2(best)
    
    # IRR table (with depth_2 labels if present)
    irr_tab <- broom.mixed::tidy(best, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
      mutate(term = dplyr::case_when(
        term == "(Intercept)"         ~ "Intercept",
        term == "diameter_z"          ~ "Diameter (z)",
        term == "log_reads"           ~ "log(reads+1)",
        grepl("^Position", term)      ~ sub("^Position", "Position: ", term) %>% sub("$", paste0(" vs ", levels(df$Position)[1])),
        grepl("^depth_2", term)       ~ sub("^depth_2", "Depth: ", term) %>% sub("OUTER", "OUTER vs INNER"),  # >>> NEW
        grepl("^ds_fac\\.L$", term)   ~ "Decay (linear)",
        grepl("^ds_fac\\.Q$", term)   ~ "Decay (quadratic)",
        grepl("^ds_fac\\.C$", term)   ~ "Decay (cubic)",
        grepl("^ds_fac\\^4$", term)   ~ "Decay (quartic)",
        TRUE                          ~ term
      )) %>%
      dplyr::rename(IRR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value) %>%
      arrange(term)
    
    # EMMeans: Position + ds_fac (+ depth_2 if available)
    emm_pos <- emmeans(best, ~ Position, type = "response")
    p_emm_pos <- summary(emm_pos) %>% as.data.frame() %>%
      ggplot(aes(x = Position, y = response, ymin = asymp.LCL, ymax = asymp.UCL)) +
      geom_point(size = 3) + geom_errorbar(width = 0.2) +
      labs(x = "Microhabitat", y = "Estimated richness",
           title = paste0("Richness ~ microhabitat (INTRA) — ", scope)) +
      theme_classic(base_size = 11)
    ggsave(file.path(outdir, paste0("alpha_richness_emm_Position_", scope, ".png")),
           p_emm_pos, width = 6, height = 3.2, dpi = 300)
    print(p_emm_pos)
    
    emm_ds  <- emmeans(best, ~ ds_fac,  type = "response")
    ds_df <- summary(emm_ds) %>% as.data.frame()
    p_emm_ds <- ggplot(ds_df, aes(x = ds_fac, y = response)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
      labs(x = "Decay at drill (ordered 1→5)", y = "Estimated richness",
           title = paste0("Richness ~ decay (INTRA) — ", scope)) +
      theme_classic(base_size = 11)
    ggsave(file.path(outdir, paste0("alpha_richness_emm_DS_", scope, ".png")),
           p_emm_ds, width = 6.2, height = 3.2, dpi = 300)
    print(p_emm_ds)
    
    if (has_depth2) {  # >>> NEW (depth_2)
      emm_depth <- emmeans(best, ~ depth_2, type = "response")
      depth_df <- summary(emm_depth) %>% as.data.frame()
      p_emm_depth <- ggplot(depth_df, aes(x = depth_2, y = response)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
        labs(x = "Depth (INNER/OUTER)", y = "Estimated richness",
             title = paste0("Richness ~ depth (INTRA) — ", scope)) +
        theme_classic(base_size = 11)
      ggsave(file.path(outdir, paste0("alpha_richness_emm_depth2_", scope, ".png")),
             p_emm_depth, width = 4.6, height = 3.2, dpi = 300)
      print(p_emm_depth)
    }
    
    # Compact publication tables
    keep_eff <- c("Position","depth_2","ds_fac","diameter_z","log_reads")  # >>> NEW includes depth_2
    tbl_main <- tibble::tibble(
      Effect  = rownames(an_tab),
      Chisq   = an_tab$`Chisq`,
      df      = an_tab$Df,
      p_value = an_tab$`Pr(>Chisq)`
    ) %>%
      dplyr::filter(Effect %in% keep_eff) %>%
      mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      arrange(match(Effect, keep_eff))
    
    tbl_r2 <- tibble::tibble(
      Metric = c("R2_marginal", "R2_conditional"),
      Value  = c(r2_tab$R2_marginal, r2_tab$R2_conditional)
    )
    
    cat("\n=== Publication tables (objects) ===\n")
    print(list(anova = tbl_main, r2 = tbl_r2, irr = irr_tab))
    
    capture.output(list(model = summary(best), anova = an_tab, r2 = r2_tab),
                   file = file.path(outdir, paste0("alpha_glmm_", scope, ".txt")))
  }
  
  cat("\n --------------  DONE FOR: ", scope, "  -------------- \n")
  invisible(NULL)
}

# ---------------------------------------------------
#  SCRIPT 7
# ---------------------------------------------------

has_levels <- function(x) is.factor(x) && nlevels(droplevels(x)) > 1

# ---------------------------------------------------
#  SCRIPT 8
# ---------------------------------------------------

fit_decay_models <- function(curve_df) {
  # curve_df: columns source, order, zeta (0..1 after normalisation)
  curve_df <- curve_df %>% mutate(zeta = pmax(zeta, 1e-12)) %>% arrange(source, order)
  sources <- unique(curve_df$source)
  
  fit_one <- function(df_src) {
    src <- df_src$source[1]
    if (nrow(df_src) < 2 || length(unique(df_src$order)) < 2) {
      return(tibble(source = src, model = c("exponential","power-law"),
                    a = NA_real_, b = NA_real_, AIC = NA_real_, half_life_orders = NA_real_))
    }
    # Exponential: log(zeta) ~ (order-1)
    m_exp <- lm(log(zeta) ~ I(order - 1), data = df_src)
    exp_row <- tibble(source = src, model = "exponential",
                      a = exp(coef(m_exp)[1]), b = -coef(m_exp)[2],
                      AIC = AIC(m_exp),
                      half_life_orders = log(2) / pmax(-coef(m_exp)[2], .Machine$double.eps))
    # Power-law: log(zeta) ~ log(order)
    m_pwr <- lm(log(zeta) ~ log(order), data = df_src)
    pwr_row <- tibble(source = src, model = "power-law",
                      a = exp(coef(m_pwr)[1]), b = -coef(m_pwr)[2],
                      AIC = AIC(m_pwr), half_life_orders = NA_real_)
    bind_rows(exp_row, pwr_row)
  }
  
  fits <- curve_df %>% group_by(source) %>% group_split() %>% purrr::map_dfr(fit_one)
  ords <- sort(unique(curve_df$order))
  pred <- expand.grid(order = ords, source = sources, model = c("exponential","power-law")) %>%
    as_tibble() %>% left_join(fits, by = c("source","model")) %>%
    mutate(pred = case_when(
      model == "exponential" ~ a * exp(-b * (order - 1)),
      model == "power-law"   ~ a * (order ^ (-b)),
      TRUE ~ NA_real_
    ))
  list(fits = fits, pred = pred)
}

plot_retention <- function(curve_df, title, outdir) {
  ret_df <- curve_df %>%
    arrange(source, order) %>% group_by(source) %>%
    mutate(retention = zeta / lag(zeta)) %>% filter(order >= 2)
  p_ret <- ggplot(ret_df, aes(order, retention, color = source)) +
    geom_point() + geom_line() +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_continuous(breaks = unique(ret_df$order)) +
    labs(x = "Zeta order (k)", y = expression(zeta[k] / zeta[k-1]),
         title = paste0(title, " — Zeta retention"), color = "Comparison") +
    theme_bw(11)
  fn <- file.path(outdir, paste0("zeta_retention_", gsub("[ +]", "_", title), ".png"))
  ggsave(fn, p_ret, width = 8, height = 6, dpi = 300)
  print(p_ret)
  invisible(list(df = ret_df, file = fn))
}

plot_with_fits <- function(curve_df, fit_obj, title, outdir) {
  pred <- fit_obj$pred
  p <- ggplot(curve_df, aes(order, zeta)) +
    geom_point(aes(color = source)) + geom_line(aes(color = source)) +
    geom_line(data = pred, aes(order, pred, linetype = model), linewidth = 0.7) +
    facet_wrap(~source) +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_continuous(breaks = sort(unique(curve_df$order))) +
    labs(x = "Zeta order (k)", y = expression(zeta[k] / zeta[1]),
         title = paste0(title, " — Normalised ζ decline with fits"),
         color = "Curve", linetype = "Fit") +
    theme_bw(11)
  fn_plot  <- file.path(outdir, paste0("zeta_decline_with_fits_", gsub("[ +]", "_", title), ".png"))
  fn_stats <- file.path(outdir, paste0("zeta_fit_stats_", gsub("[ +]", "_", title), ".tsv"))
  ggsave(fn_plot, p, width = 9, height = 6, dpi = 300)
  print(p)
  readr::write_tsv(fit_obj$fits, fn_stats)
  print(fit_obj$fits); cat("Fit table:", fn_stats, "\n")
  invisible(list(file_plot = fn_plot, file_stats = fn_stats))
}
# Build one presence vector per tree using position priority,Use one presence vector per tree, built as the union of depths within the first available position in c("MIDDLE","BASE","UPPER"). If none exist for a tree, it falls back to all samples of that tree.
build_tree_presence_representative <- function(pa, meta_sub,
                                               pos_priority = c("MIDDLE","BASE","UPPER")) {
  stopifnot(identical(rownames(pa), meta_sub$sample))
  # split row indices by tree
  trees <- split(seq_len(nrow(meta_sub)), meta_sub$natman)
  pres  <- lapply(trees, function(idx) {
    m <- meta_sub[idx, , drop = FALSE]
    o <- pa[idx, , drop = FALSE]
    got <- NULL
    for (pp in pos_priority) {
      sel <- which(m$position_2 == pp)
      if (length(sel)) { got <- o[sel, , drop = FALSE]; break }
    }
    if (is.null(got)) got <- o
    as.integer(colSums(got > 0) > 0)
  })
  out <- do.call(rbind, pres)
  rownames(out) <- names(pres)
  as.data.frame(out, stringsAsFactors = FALSE)
}

# --- Core runner (per set of dw_types) ---------------------------------------
run_zeta_norm <- function(meta_filt, otu_matrix, dw_types, title, outdir = "plots/ZETA_unscaled", k_cap = 6) {
  meta_sub <- meta_filt %>%
    filter(dw_type %in% dw_types, sample %in% rownames(otu_matrix)) %>% droplevels()
  if (nrow(meta_sub) < 2) { message("[", title, "] too few samples."); return(invisible(NULL)) }
  
  Xs <- otu_matrix[meta_sub$sample, , drop = FALSE]
  pa <- (Xs > 0) * 1; pa <- pa[, colSums(pa) > 0, drop = FALSE]
  if (ncol(pa) == 0L) { message("[", title, "] no taxa."); return(invisible(NULL)) }
  
  k_max <- min(k_cap, nrow(pa))
  
  # WITHIN — compute per-tree zeta, normalise per tree by ζ1, then average at each order
  ids <- meta_sub %>% count(natman, name = "n") %>% filter(n >= 2) %>% pull(natman)
  within_norm <- purrr::map_dfr(ids, function(id) {
    rows <- meta_sub$sample[meta_sub$natman == id]
    sub  <- pa[rows, , drop = FALSE]
    kk   <- 1:min(k_max, nrow(sub))
    z    <- zetadiv::Zeta.decline.ex(as.data.frame(sub), orders = kk, rescale = FALSE)
    tibble(source = paste0("within:", id),
           order = z$zeta.order,
           zeta  = z$zeta.val / max(z$zeta.val[z$zeta.order == 1], .Machine$double.eps))
  })
  within_avg <- within_norm %>%
    group_by(order) %>%
    summarise(zeta = mean(zeta, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = "within")
  
  # BETWEEN — collapse to presence at tree level; normalise by ζ1 of between curve
  pa_between <- build_tree_presence_representative(pa, meta_sub,
                                                   pos_priority = c("MIDDLE","BASE","UPPER"))
  
  between_df <- tibble()
  if (nrow(pa_between) >= 2) {
    kkB <- 1:min(k_max, nrow(pa_between))
    zb  <- zetadiv::Zeta.decline.ex(pa_between, orders = kkB, rescale = FALSE)
    z1  <- max(zb$zeta.val[zb$zeta.order == 1], .Machine$double.eps)
    between_df <- tibble(order = zb$zeta.order, zeta = zb$zeta.val / z1, source = "between")
  }
  
  zeta_both <- bind_rows(within_avg, between_df)
  if (!nrow(zeta_both)) { message("[", title, "] no zeta curve."); return(invisible(NULL)) }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Base plot
  p_base <- ggplot(zeta_both, aes(order, zeta, color = source)) +
    geom_point() + geom_line() +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_continuous(breaks = unique(zeta_both$order)) +
    labs(x = "Zeta order (k)", y = expression(zeta[k] / zeta[1]),
         title = paste0(paste(dw_types, collapse = "+"),
                        " — Zeta decline (within vs between trees, normalised)"),
         color = "Comparison") +
    theme_bw(11)
  fn_base <- file.path(outdir, paste0("zeta_within_between_norm_", gsub("[ +]", "_", title), ".png"))
  ggsave(fn_base, p_base, width = 8, height = 6, dpi = 300)
  print(p_base)
  p_base
  # Retention + model fits
  ret_out <- plot_retention(zeta_both, title = paste(title, "(", paste(dw_types, collapse="+"), ")"), outdir = outdir)
  fit_obj <- fit_decay_models(zeta_both)
  fit_out <- plot_with_fits(zeta_both, fit_obj, title = paste(title, "(", paste(dw_types, collapse="+"), ")"), outdir = outdir)
  
  # Save curve
  write_csv(zeta_both, file.path(outdir, paste0("zeta_curves_", gsub("[ +]", "_", title), ".csv")))
  invisible(list(curve = zeta_both, base_plot = fn_base,
                 retention = ret_out$file, fit_plot = fit_out$file_plot,
                 fit_stats = fit_out$file_stats))
}
# --------------------------------------
#   absolute versions
# === ABSOLUTE versions of the plotting helpers (no normalisation in labels) ===
plot_retention_abs <- function(curve_df, title, outdir) {
  ret_df <- curve_df %>%
    arrange(source, order) %>%
    group_by(source) %>%
    mutate(retention = zeta / dplyr::lag(zeta)) %>%
    filter(order >= 2)
  
  p_ret <- ggplot(ret_df, aes(order, retention, color = source)) +
    geom_point() + geom_line() +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_continuous(breaks = unique(ret_df$order)) +
    labs(x = "Zeta order (k)",
         y = expression("Retention  " * zeta[k] / zeta[k-1]),
         title = paste0(title, " — Zeta retention"),
         color = "Comparison") +
    theme_bw(11)
  
  fn_ret <- file.path(outdir, paste0("zeta_retention_ABS_", gsub("[ +]", "_", title), ".png"))
  ggsave(fn_ret, p_ret, width = 8, height = 6, dpi = 300)
  invisible(list(df = ret_df, plot = p_ret, file = fn_ret))
  print(fn_ret)
}

plot_with_fits_abs <- function(curve_df, fit_obj, title, outdir) {
  pred <- fit_obj$pred
  p <- ggplot(curve_df, aes(order, zeta)) +
    geom_point(aes(color = source)) +
    geom_line(aes(color = source)) +
    geom_line(data = pred, aes(order, pred, linetype = model), linewidth = 0.7) +
    facet_wrap(~source) +
    coord_cartesian(ylim = c(0, max(curve_df$zeta, na.rm=TRUE) * 1.05)) +
    scale_x_continuous(breaks = sort(unique(curve_df$order))) +
    labs(x = "Zeta order (k)",
         y = expression(zeta[k]),
         title = paste0(title, " — ζ decline with model fits (absolute)"),
         color = "Curve", linetype = "Fit") +
    theme_bw(11)
  
  fn <- file.path(outdir, paste0("zeta_decline_with_fits_ABS_", gsub("[ +]", "_", title), ".png"))
  ggsave(fn, p, width = 9, height = 6, dpi = 300)
  print(p)
  stats_fn <- file.path(outdir, paste0("zeta_fit_stats_ABS_", gsub("[ +]", "_", title), ".tsv"))
  readr::write_tsv(fit_obj$fits, stats_fn)
  invisible(list(plot = p, file_plot = fn, file_stats = stats_fn))
}
# === ABSOLUTE curve version (no per-curve normalisation by ζ1) ===
run_zeta_abs <- function(meta_filt, otu_matrix, dw_types, title,
                         outdir = "plots/ZETA_unscaled", k_cap = 6) {
  
  meta_sub <- meta_filt %>%
    dplyr::filter(dw_type %in% dw_types, sample %in% rownames(otu_matrix)) %>%
    droplevels()
  if (nrow(meta_sub) < 2) { message("[", title, "] too few samples."); return(invisible(NULL)) }
  
  Xs <- otu_matrix[meta_sub$sample, , drop = FALSE]
  pa <- (Xs > 0) * 1
  pa <- pa[, colSums(pa) > 0, drop = FALSE]
  if (ncol(pa) == 0L) { message("[", title, "] no taxa."); return(invisible(NULL)) }
  
  k_max <- min(k_cap, nrow(pa))
  
  # WITHIN — per-tree ζ curves (absolute), then average ζ_k across trees
  ids <- meta_sub %>% count(natman, name = "n") %>% filter(n >= 2) %>% pull(natman)
  within_abs <- purrr::map_dfr(ids, function(id) {
    rows <- meta_sub$sample[meta_sub$natman == id]
    sub  <- pa[rows, , drop = FALSE]
    kk   <- 1:min(k_max, nrow(sub))
    z    <- zetadiv::Zeta.decline.ex(as.data.frame(sub), orders = kk, rescale = FALSE)
    tibble(source = paste0("within:", id),
           order = z$zeta.order,
           zeta  = z$zeta.val)  # <- absolute counts
  })
  
  within_avg <- within_abs %>%
    group_by(order) %>%
    summarise(zeta = mean(zeta, na.rm = TRUE), .groups = "drop") %>%
    mutate(source = "within")
  
  # BETWEEN — collapse to tree level; absolute ζ
  pa_between <- build_tree_presence_representative(pa, meta_sub,
                                                   pos_priority = c("MIDDLE","BASE","UPPER"))
  
  
  between_df <- tibble()
  if (nrow(pa_between) >= 2) {
    kkB <- 1:min(k_max, nrow(pa_between))
    zb  <- zetadiv::Zeta.decline.ex(pa_between, orders = kkB, rescale = FALSE)
    between_df <- tibble(order = zb$zeta.order, zeta = zb$zeta.val, source = "between")
  }
  
  zeta_both <- bind_rows(within_avg, between_df)
  if (!nrow(zeta_both)) { message("[", title, "] no zeta curve."); return(invisible(NULL)) }
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # Base plot (absolute)
  p_base <- ggplot(zeta_both, aes(order, zeta, color = source)) +
    geom_point() + geom_line() +
    scale_x_continuous(breaks = unique(zeta_both$order)) +
    labs(x = "Zeta order (k)", y = expression(zeta[k]),
         title = paste0(paste(dw_types, collapse = "+"),
                        " — Zeta decline (within vs between trees, absolute)"),
         color = "Comparison") +
    theme_bw(11)
  
  fn_base <- file.path(outdir, paste0("zeta_within_between_ABS_", gsub("[ +]", "_", title), ".png"))
  ggsave(fn_base, p_base, width = 8, height = 6, dpi = 300)
  print(p_base)
  
  # Retention + model fits (absolute)
  ret_out <- plot_retention_abs(zeta_both, title = paste(title, "(", paste(dw_types, collapse="+"), ")"), outdir = outdir)
  fit_obj <- fit_decay_models(zeta_both)  # works for absolute too (log-ζ safeguards inside)
  fit_out <- plot_with_fits_abs(zeta_both, fit_obj, title = paste(title, "(", paste(dw_types, collapse="+"), ")"), outdir = outdir)
  
  # Save curve data
  readr::write_csv(zeta_both, file.path(outdir, paste0("zeta_curves_ABS_", gsub("[ +]", "_", title), ".csv")))
  
  invisible(list(curve = zeta_both,
                 base_plot = fn_base,
                 retention = ret_out$file,
                 fit_plot = fit_out$file_plot,
                 fit_stats = fit_out$file_stats))
}
