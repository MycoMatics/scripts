# =============================================================================
# Script: Script11_figures_and_supplementary_panels.R
# Project: Deadwood fungal community analyses
# Purpose: Supplementary figures, ordination panel, and null-distribution figure
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
# - tax
# - threshold
# - exclude_samples
#
# Required helper functions from Script0_utils.R:
# - ensure_dir
# - assert_objects
# - as_pa
# - combine_counts
# =============================================================================

source("Scripts/Script0_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(iNEXT)
  library(janitor)
  library(vegan)
  library(ggrepel)
  library(scales)
  library(cowplot)
  library(ggridges)
  library(compositions)
})

set.seed(SEED_GLOBAL)
setwd(".")

assert_objects(c("META1", "otu_matrix_filt", "tax", "threshold", "exclude_samples"))

ensure_dir("plots")
ensure_dir(file.path("plots", "ALL"))
ensure_dir(file.path("plots", "NATMAN"))
ensure_dir(file.path("plots", "ORDINATIONS"))

# =============================================================================
# SUPPLEMENTARY FIGURE 2
# =============================================================================
beta_txt <- "
contrasts\tscope\tSIM\tSNE\tSOR\tn
between\t1: CROSS substrates\t0,786\t0,071\t0,889\t17726
within\t1: CROSS substrates\t0,691\t0,118\t0,849\t14152
between\t2A: LOG\t0,692\t0,122\t0,855\t7702
within\t2A: LOG\t0,597\t0,123\t0,745\t299
between\t2Asub: LOG: INNER–OUTER\t0\t0\t0\t0
within\t2Asub: LOG: INNER–OUTER\t0,462\t0,159\t0,697\t66
between\t2B: SNAG\t0,700\t0,114\t0,875\t220
within\t2B: SNAG\t0,464\t0,084\t0,663\t11
between\t2C: aFWD\t0,714\t0,090\t0,824\t546
within\t2C: aFWD\t0,571\t0,142\t0,628\t15
between\t2D: fFWD\t0,667\t0,067\t0,755\t408
within\t2D: fFWD\t0,608\t0,070\t0,686\t187
between\t3A: LOG_aFWD\t0,744\t0,089\t0,857\t898
within\t3A: LOG_aFWD\t0,664\t0,089\t0,784\t74
between\t3B: FWD paired\t0,679\t0,074\t0,766\t70
within\t3B: FWD paired\t0,606\t0,071\t0,685\t190
between\t3C: fFWD–SOIL\t0,681\t0,178\t0,894\t121
within\t3C: fFWD–SOIL\t0,461\t0,063\t0,553\t110
" %>%
  stringr::str_trim()

parse_num_comma_decimal <- function(x) {
  parse_number(x, locale = locale(decimal_mark = ","))
}

beta <- read_delim(
  I(beta_txt),
  delim = "\t",
  col_types = cols(.default = col_character())
) %>%
  mutate(
    contrasts = factor(contrasts, levels = c("within", "between")),
    SIM = parse_num_comma_decimal(SIM),
    SNE = parse_num_comma_decimal(SNE),
    SOR = parse_num_comma_decimal(SOR),
    n   = as.integer(n)
  )

beta_long <- beta %>%
  dplyr::select(contrasts, scope, SIM, SNE) %>%
  tidyr::pivot_longer(c(SIM, SNE), names_to = "component", values_to = "beta") %>%
  dplyr::mutate(
    component = dplyr::case_when(
      component == "SIM" ~ "Turnover (βSIM)",
      component == "SNE" ~ "Nestedness (βSNE)",
      TRUE ~ component
    ),
    component = factor(
      component,
      levels = c("Turnover (βSIM)", "Nestedness (βSNE)")
    )
  )
scope_order <- beta %>% distinct(scope) %>% pull(scope)

beta_long2 <- beta_long %>%
  mutate(
    scope = factor(scope, levels = scope_order),
    x = paste(as.character(scope), contrasts, sep = " – "),
    x = factor(
      x,
      levels = as.vector(rbind(
        paste(scope_order, "within", sep = " – "),
        paste(scope_order, "between", sep = " – ")
      ))
    )
  )

bars <- beta_long2 %>%
  select(x, scope, contrasts, component, beta)

bar_anno <- beta %>%
  mutate(
    scope = factor(scope, levels = scope_order),
    x = paste(as.character(scope), contrasts, sep = " – "),
    x = factor(
      x,
      levels = as.vector(rbind(
        paste(scope_order, "within",  sep = " – "),
        paste(scope_order, "between", sep = " – ")
      ))
    )
  ) %>%
  transmute(
    x,
    scope,
    contrasts,
    total = SIM + SNE,
    n = n,
    TN = if_else(SNE > 0, SIM / SNE, NA_real_)
  ) %>%
  mutate(
    label_n = if_else(total > 0, paste0("n=\n", format(n, big.mark = ",")), ""),
    label_TN = if_else(
      total > 0 & !is.na(TN),
      paste0("T:N=\n", formatC(TN, format = "f", digits = 1)),
      ""
    )
  )

x_labels <- gsub(" – within| – between", "", levels(beta_long2$x))

p_sup2 <- ggplot() +
  geom_col(
    data = bars,
    aes(x = x, y = beta, fill = component),
    width = 0.85,
    colour = "black",
    linewidth = 0.2
  ) +
  geom_col(
    data = bar_anno,
    aes(x = x, y = total, alpha = contrasts),
    width = 0.85,
    fill = "grey30",
    colour = NA
  ) +
  geom_text(
    data = bar_anno,
    aes(x = x, y = total * 0.6, label = label_TN),
    size = 3
  ) +
  geom_text(
    data = bar_anno,
    aes(x = x, y = total, label = label_n),
    vjust = -0.6,
    size = 3.2,
    fontface = "bold"
  ) +
  scale_alpha_manual(
    values = c(within = 0.55, between = 0.1),
    name = "Contrast"
  ) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = NULL,
    y = "Median β-diversity (Sørensen partition)",
    fill = "Beta component"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank()
  )

print(p_sup2)
ggsave(
  filename = file.path("plots", "Supplementary_Figure_2_beta_decomposition.png"),
  plot = p_sup2,
  dpi = 900,
  width = 8.5,
  height = 5.5
)

# =============================================================================
# SUPPLEMENTARY FIGURE 3
# =============================================================================
meta <- META1 %>% mutate(sample = as.character(sample))
samps <- intersect(rownames(otu_matrix_filt), meta$sample)
stopifnot(length(samps) > 0)

meta <- meta %>% filter(sample %in% samps)
otu  <- as.matrix(otu_matrix_filt[samps, , drop = FALSE])

collapse_by_group <- function(M, groups, method = c("geo", "sum", "mean", "median"), pseudo = 1) {
  method <- match.arg(method)
  groups <- as.character(groups)
  keep <- !is.na(groups)
  
  M <- as.matrix(M[keep, , drop = FALSE])
  groups <- groups[keep]
  
  if (nrow(M) == 0) stop("No rows remain after removing NA groups.")
  if (length(groups) != nrow(M)) stop("Length of `groups` must match nrow(M).")
  
  idx <- split(seq_len(nrow(M)), groups)
  
  collapsed <- t(vapply(
    idx,
    function(ii) combine_counts(M, ii, method = method, pseudo = pseudo),
    numeric(ncol(M))
  ))
  
  colnames(collapsed) <- colnames(M)
  rownames(collapsed) <- names(idx)
  collapsed
}

compute_ofd_from_mat <- function(PA_mat, label) {
  stopifnot(nrow(PA_mat) >= 1, ncol(PA_mat) >= 1)
  
  occ <- colSums(PA_mat)
  tibble(k = as.integer(occ)) %>%
    count(k, name = "n_taxa") %>%
    mutate(
      scope   = label,
      n_units = nrow(PA_mat),
      prop_k  = k / n_units,
      freq    = n_taxa / sum(n_taxa)
    ) %>%
    relocate(scope, n_units, k, prop_k, n_taxa, freq)
}

PA_sample <- as_pa(otu)
ofd_per_sample <- compute_ofd_from_mat(PA_sample, label = "per_sample")

M_natman_gm <- collapse_by_group(otu, groups = meta$natman, method = "geo", pseudo = 1)
PA_natman <- as_pa(M_natman_gm)
ofd_per_natman <- compute_ofd_from_mat(PA_natman, label = "per_natman_gm")

meta_log <- meta %>% filter(dw_type2 == "LOG")
M_log <- as.matrix(otu[meta_log$sample, , drop = FALSE])
M_natman_LOG_gm <- collapse_by_group(M_log, groups = meta_log$natman, method = "geo", pseudo = 1)
PA_natman_LOG <- as_pa(M_natman_LOG_gm)
ofd_per_natman_LOG <- compute_ofd_from_mat(PA_natman_LOG, label = "per_natman_LOG_gm")

common_scale <- list(
  scale_y_continuous(
    limits = c(0, 1100),
    breaks = seq(0, 1000, by = 250),
    expand = expansion(mult = c(0, 0))
  ),
  theme_bw(8)
)

p1 <- ggplot(ofd_per_sample, aes(k, n_taxa)) +
  geom_col() +
  labs(
    title = "Occupancy frequency\ndistribution across samples",
    x = "Occupancy (k samples)",
    y = "Number of SHs"
  ) +
  common_scale

p2 <- ggplot(ofd_per_natman, aes(k, n_taxa)) +
  geom_col() +
  labs(
    title = "OFD across trees,\ngeometrically collapsed",
    x = "Occupancy (k trees)",
    y = ""
  ) +
  common_scale

p3 <- ggplot(ofd_per_natman_LOG, aes(k, n_taxa)) +
  geom_col() +
  labs(
    title = "OFD across LOGs,\ngeometrically collapsed",
    x = "Occupancy (k trees)",
    y = ""
  ) +
  common_scale

panel_ABC <- (p1 + p3 + p2) + plot_annotation(tag_levels = "A")
print(panel_ABC)

ggsave(
  file.path("plots", "ALL", "OFD_panel.png"),
  panel_ABC,
  width = 7.8,
  height = 4,
  dpi = 900
)

# =============================================================================
# SUPPLEMENTARY FIGURE 5
# =============================================================================
meta_filt <- META1 %>%
  filter(!sample %in% exclude_samples) %>%
  mutate(
    dw_type = dw_type2,
    natman  = droplevels(factor(natman)),
    umi     = as.factor(umi)
  ) %>%
  filter(dw_type %in% c("SNAG", "aFWD", "fFWD", "LOG")) %>%
  droplevels()

X <- otu_matrix_filt[meta_filt$sample, , drop = FALSE]
stopifnot(identical(rownames(X), meta_filt$sample))
X <- X[, colSums(X) > 0, drop = FALSE]
X_mat <- as.matrix(X)

meta_use <- meta_filt %>%
  mutate(
    sample   = as.character(sample),
    natman   = as.character(natman),
    dw_type2 = as.character(dw_type2)
  ) %>%
  filter(sample %in% rownames(X_mat))

X_use <- X_mat[meta_use$sample, , drop = FALSE]
X_use <- X_use[rowSums(X_use) > 0, , drop = FALSE]
meta_use <- meta_use %>% filter(sample %in% rownames(X_use))
stopifnot(identical(rownames(X_use), meta_use$sample))

keep_ids <- c("ZF1017", "ZF401", "ZF1002", "ZF305")

log_ranked <- meta_use %>%
  filter(dw_type2 == "LOG") %>%
  count(natman, name = "n_log_samples") %>%
  arrange(desc(n_log_samples), natman)

top8_ids <- log_ranked %>% slice_head(n = 7) %>% pull(natman)
log_natman <- unique(c(top8_ids, keep_ids))
log_natman <- log_natman[seq_len(min(10, length(log_natman)))]

top_all_ranked <- meta_use %>%
  count(natman, name = "n_all_samples") %>%
  arrange(desc(n_all_samples), natman)

all_natman <- unique(c(
  log_natman,
  top_all_ranked %>%
    filter(!natman %in% log_natman) %>%
    slice_head(n = max(0, 10 - length(log_natman))) %>%
    pull(natman)
))
all_natman <- all_natman[seq_len(min(10, length(all_natman)))]

natman_levels <- all_natman

natman_sets_tbl <- tibble(
  natman = natman_levels,
  in_all_panel = natman %in% all_natman,
  in_log_panel = natman %in% log_natman,
  n_all_samples = top_all_ranked$n_all_samples[match(natman_levels, top_all_ranked$natman)],
  n_log_samples = log_ranked$n_log_samples[match(natman_levels, log_ranked$natman)]
) %>%
  mutate(
    n_all_samples = replace_na(n_all_samples, 0L),
    n_log_samples = replace_na(n_log_samples, 0L)
  )

print(natman_sets_tbl, n = Inf)

write_csv2(
  natman_sets_tbl,
  file = file.path("plots", "NATMAN", paste0("natman_top10_all_vs_LOG_threshold_", threshold, ".csv"))
)

make_inext_natman_df <- function(meta_df, X_df, natman_keep, endpoint_mult = 2, endpoint_cap = 200) {
  meta_sub <- meta_df %>%
    filter(natman %in% natman_keep) %>%
    mutate(natman = factor(natman, levels = natman_levels))
  
  X_sub <- X_df[meta_sub$sample, , drop = FALSE]
  stopifnot(identical(rownames(X_sub), meta_sub$sample))
  
  inc_list <- split(meta_sub$sample, meta_sub$natman) %>%
    lapply(function(samps) {
      ot <- X_sub[samps, , drop = FALSE]
      ot <- ot[, colSums(ot) > 0, drop = FALSE]
      freq <- colSums(ot > 0)
      c(length(samps), freq)
    })
  
  endpoints <- sapply(inc_list, function(v) {
    n <- v[1]
    min(ceiling(n * endpoint_mult), endpoint_cap)
  })
  
  out_list <- lapply(names(inc_list), function(g) {
    out <- iNEXT(
      inc_list[[g]],
      q = c(0, 1, 2),
      datatype = "incidence_freq",
      endpoint = endpoints[[g]]
    )
    as_tibble(out$iNextEst$size_based) %>% mutate(natman = g)
  })
  
  df <- bind_rows(out_list)
  obs_n <- sapply(inc_list, function(v) v[1])
  
  df <- df %>%
    mutate(
      natman = factor(natman, levels = natman_levels),
      interp_extrap = ifelse(t <= obs_n[as.character(natman)], "interpolated", "extrapolated")
    )
  
  pts <- df %>%
    group_by(natman, Order.q) %>%
    filter(t == max(t[interp_extrap == "interpolated"], na.rm = TRUE)) %>%
    ungroup()
  
  list(df = df, pts = pts, obs_n = obs_n)
}

res_all <- make_inext_natman_df(meta_df = meta_use, X_df = X_use, natman_keep = all_natman)

meta_log <- meta_use %>% filter(dw_type2 == "LOG")
X_log <- X_use[meta_log$sample, , drop = FALSE]
stopifnot(identical(rownames(X_log), meta_log$sample))

res_log <- make_inext_natman_df(meta_df = meta_log, X_df = X_log, natman_keep = log_natman)

order_labs <- c(
  "0" = "Richness (q = 0)\n(Number of SHs)",
  "1" = "Exp Shannon (q = 1)\n(Number of typical SHs)",
  "2" = "Inverse Simpson (q = 2)\n(Number of dominant SHs)"
)

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 9),
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(t = 14, r = 14, b = 10, l = 14)
  )

plot_inext_panel <- function(df, pts, title_txt, ylims_tbl) {
  df2 <- df %>% mutate(Order.q = as.numeric(as.character(Order.q)))
  
  anchors <- ylims_tbl %>% transmute(Order.q, t = 0, qD = ymin)
  anchors_top <- ylims_tbl %>% transmute(Order.q, t = 0, qD = ymax)
  
  ggplot(df2, aes(x = t, y = qD)) +
    geom_blank(data = anchors, aes(x = t, y = qD)) +
    geom_blank(data = anchors_top, aes(x = t, y = qD)) +
    geom_line(
      data = df2 %>% filter(interp_extrap == "interpolated"),
      aes(color = natman, group = natman),
      linewidth = 0.8
    ) +
    geom_line(
      data = df2 %>% filter(interp_extrap == "extrapolated"),
      aes(color = natman, group = natman),
      linewidth = 0.8,
      linetype = "dotted"
    ) +
    geom_ribbon(
      aes(ymin = qD.LCL, ymax = qD.UCL, fill = natman),
      alpha = 0.06,
      color = NA,
      show.legend = FALSE
    ) +
    geom_point(
      data = pts,
      aes(color = natman),
      size = 1.8,
      shape = 21,
      fill = "white",
      stroke = 0.7,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ Order.q,
      ncol = 1,
      labeller = as_labeller(order_labs),
      scales = "fixed"
    ) +
    coord_cartesian(xlim = c(0, 21)) +
    labs(
      title = title_txt,
      x = "Number of samples\n(observed & extrapolated)",
      y = "Estimated diversity (qD)",
      color = "natman"
    ) +
    scale_y_continuous(
      limits = c(unique(ylims_tbl$ymin), unique(ylims_tbl$ymax)),
      expand = expansion(mult = c(0, 0)),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    base_theme
}

get_legend_grob <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(NULL)
  g$grobs[[idx[1]]]
}

y_limits_by_q <- tibble(Order.q = c(0, 1, 2), ymin = 0, ymax = 750)

p_all <- plot_inext_panel(res_all$df, res_all$pts, "All deadwood types\n(dw_type2)", y_limits_by_q)
p_log <- plot_inext_panel(res_log$df, res_log$pts, "LOG only\n(dw_type2 == \"LOG\")", y_limits_by_q)

legend_grob <- get_legend_grob(p_log + theme(legend.position = "bottom"))
p_all_noleg <- p_all + theme(legend.position = "none")
p_log_noleg <- p_log + theme(legend.position = "none")

p_combo <- (p_all_noleg | p_log_noleg) /
  patchwork::wrap_elements(full = legend_grob) +
  patchwork::plot_layout(heights = c(1, 0.10))

print(p_combo)

a4_width_mm  <- 210
a4_height_mm <- 297
margin_lr_mm <- 30
margin_tb_mm <- 20

plot_width_in  <- (a4_width_mm  - 2 * margin_lr_mm) / 25.4
plot_height_in <- (a4_height_mm - 2 * margin_tb_mm) / 25.4

ggsave(
  filename = file.path("plots", "NATMAN", paste0("iNEXT_natman_all_vs_LOG_top7-ffwd_A4_threshold_", threshold, ".png")),
  plot = p_combo,
  width = plot_width_in,
  height = plot_height_in,
  units = "in",
  dpi = 900,
  bg = "white"
)

# =============================================================================
# FIGURE 4
# =============================================================================
dw_colors2 <- dw_colors
dw_colors2["aFWD"] <- "#FFBF33"
dw_colors2["LOG"]  <- "#7A4F26"
dw_cols <- dw_colors2

outdir_ord <- file.path("plots", "ORDINATIONS")
ensure_dir(outdir_ord)

meta_filt <- META1 %>%
  filter(!sample %in% exclude_samples) %>%
  mutate(
    dw_type = dw_type2,
    natman = droplevels(factor(natman)),
    umi = as.factor(umi),
    decay_stage = factor(decay_stage)
  ) %>%
  filter(dw_type %in% c("SNAG", "aFWD", "fFWD", "LOG")) %>%
  filter(decay_stage != "LIVING") %>%
  droplevels()

X <- otu_matrix_filt[meta_filt$sample, , drop = FALSE]
stopifnot(identical(rownames(X), meta_filt$sample))
X <- X[, colSums(X) > 0, drop = FALSE]

scaling_choice <- 2L
do_flip_axes   <- TRUE
prevalence_min <- 0L

CLR_from_X <- function(X, prevalence_min = 0L) {
  X <- as.matrix(X)
  if (prevalence_min > 0L) {
    keep <- colSums(X > 0) >= prevalence_min
    X <- X[, keep, drop = FALSE]
  }
  Xz <- zCompositions::cmultRepl(
    X,
    method = "CZM",
    label = 0,
    z.warning = TRUE,
    z.delete = FALSE
  )
  compositions::clr(Xz)
}

meta2 <- meta_filt %>%
  mutate(
    reads = rowSums(X),
    log_reads = log1p(reads),
    umi = factor(umi, levels = sort(unique(umi)))
  ) %>%
  filter(reads > 0)

X <- X[meta2$sample, , drop = FALSE]
X_clr <- CLR_from_X(X, prevalence_min = prevalence_min)

common <- intersect(rownames(X_clr), meta2$sample)
meta2  <- meta2 %>% filter(sample %in% common) %>% arrange(match(sample, common))
X_clr  <- X_clr[meta2$sample, , drop = FALSE]
stopifnot(identical(rownames(X_clr), meta2$sample))

pca_rAit_part <- vegan::rda(
  X_clr ~ 1 + Condition(log_reads + umi),
  data = meta2
)

eig_ra <- vegan::eigenvals(pca_rAit_part)
var_ra <- eig_ra / sum(eig_ra)
pc1ra_lb <- scales::percent(var_ra[1], accuracy = 0.1)
pc2ra_lb <- scales::percent(var_ra[2], accuracy = 0.1)

scores_sites <- vegan::scores(
  pca_rAit_part,
  display = "sites",
  scaling = scaling_choice,
  choices = c(1, 2)
)

scores_species <- vegan::scores(
  pca_rAit_part,
  display = "species",
  scaling = scaling_choice,
  choices = c(1, 2)
)

sites_ra <- as_tibble(scores_sites, rownames = "sample") %>%
  left_join(
    select(meta2, sample, dw_type, natman, reads, log_reads, umi, ds_at_drill, decay_stage),
    by = "sample"
  ) %>%
  mutate(decay_fac = factor(decay_stage))

if (isTRUE(do_flip_axes)) {
  sites_ra <- sites_ra %>%
    mutate(PC1 = -PC1, PC2 = -PC2)
  scores_species[, 1] <- -scores_species[, 1]
  scores_species[, 2] <- -scores_species[, 2]
}

sp_ra <- as_tibble(scores_species, rownames = "feature")

top_natman_tbl <- sites_ra %>%
  count(natman, sort = TRUE) %>%
  slice_head(n = 3)

top_natman <- top_natman_tbl$natman
shape_map_top <- setNames(c(13, 14, 9), as.character(top_natman))

sites_ra <- sites_ra %>%
  mutate(
    natman_shape_group = if_else(
      as.character(natman) %in% as.character(top_natman),
      as.character(natman),
      "Other"
    ),
    natman_shape_group = factor(
      natman_shape_group,
      levels = c(as.character(top_natman), "Other")
    ),
    decay_stage_grp = case_when(
      decay_stage == "EARLY"   ~ "Early (DS 1)",
      decay_stage == "AVERAGE" ~ "Intermediate (DS 2–3)",
      decay_stage == "LATE"    ~ "Late (DS 4–5)",
      TRUE ~ NA_character_
    ),
    decay_stage_grp = factor(
      decay_stage_grp,
      levels = c("Early (DS 1)", "Intermediate (DS 2–3)", "Late (DS 4–5)")
    )
  )

shape_values <- c(shape_map_top, Other = 20)

target_prefixes <- c(
  "SH0848643.10FU",
  "SH0968712.10FU",
  "SH1016964.10FU"
)
pattern_prefix <- paste0("^(", paste(target_prefixes, collapse = "|"), ")")

sp_ra_sel <- sp_ra %>%
  filter(grepl(pattern_prefix, feature)) %>%
  mutate(
    SH   = sub("\\|.*", "", feature),
    taxo = sub(".*\\|", "", feature),
    taxo = gsub("_", " ", taxo),
    genus = sub(" .*", "", taxo),
    rest  = sub("^[^ ]+\\s*", "", taxo),
    genus_initial = substr(genus, 1, 1),
    epithet_str = if_else(rest == "", genus, rest),
    epithet_parsed = gsub(" ", "~", epithet_str),
    label_expr = if_else(
      rest == "",
      SH,
      paste0(SH, "~italic(", genus_initial, ".~", epithet_parsed, ")")
    )
  )

rng_12 <- range(c(sites_ra$PC1, sites_ra$PC2), na.rm = TRUE)
arrow_mult_12 <- 0.7 * diff(rng_12) /
  max(1e-9, max(sqrt(sp_ra_sel$PC1^2 + sp_ra_sel$PC2^2)))

label_offsets <- tibble::tribble(
  ~SH,               ~dx,   ~dy,
  "SH0848643.10FU",  -3.10, -0.15,
  "SH0968712.10FU",  -1.40,  0.25,
  "SH1016964.10FU",   0.10,  0.20
)

sp_ra_sel <- sp_ra_sel %>%
  mutate(
    x_arrow = PC1 * arrow_mult_12,
    y_arrow = PC2 * arrow_mult_12
  ) %>%
  left_join(label_offsets, by = "SH") %>%
  mutate(
    dx = coalesce(dx, 0),
    dy = coalesce(dy, 0),
    x_lab = x_arrow + dx,
    y_lab = y_arrow + dy
  )

p_pc12_dw <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.25, alpha = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.25, alpha = 0.3) +
  geom_point(
    data = sites_ra,
    aes(PC1, PC2, color = dw_type, shape = natman_shape_group),
    alpha = 0.7,
    size = 2
  ) +
  scale_color_manual(
    values = dw_cols,
    guide = guide_legend(override.aes = list(size = 2.0))
  ) +
  scale_shape_manual(
    values = shape_values,
    breaks = c(as.character(top_natman), "Other"),
    name = "Tree identity"
  ) +
  geom_segment(
    data = sp_ra_sel,
    aes(x = 0, y = 0, xend = x_arrow, yend = y_arrow),
    arrow = arrow(length = unit(0.015, "npc")),
    linewidth = 0.3,
    alpha = 0.4
  ) +
  labs(
    title = "PC1 vs PC2",
    x = paste0("PC1 (", pc1ra_lb, ")"),
    y = paste0("PC2 (", pc2ra_lb, ")"),
    color = "Dead wood type"
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "right")

ds_vals <- c(
  "Early (DS 1)" = DS_colors[["EARLY"]],
  "Intermediate (DS 2–3)" = DS_colors[["AVERAGE"]],
  "Late (DS 4–5)" = DS_colors[["LATE"]]
)

p_pc12_dec <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.25, alpha = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.25, alpha = 0.3) +
  geom_point(
    data = sites_ra %>% filter(!is.na(decay_stage_grp)),
    aes(PC1, PC2, color = decay_stage_grp, shape = natman_shape_group),
    alpha = 0.7,
    size = 2
  ) +
  geom_segment(
    data = sp_ra_sel,
    aes(x = 0, y = 0, xend = x_arrow, yend = y_arrow),
    arrow = arrow(length = unit(0.015, "npc")),
    linewidth = 0.3,
    alpha = 0.4
  ) +
  scale_color_manual(
    values = ds_vals,
    breaks = names(ds_vals),
    drop = FALSE,
    guide = guide_legend(override.aes = list(size = 2.0))
  ) +
  scale_shape_manual(
    values = shape_values,
    breaks = c(as.character(top_natman), "Other"),
    name = "Tree identity"
  ) +
  labs(
    title = "",
    x = paste0("PC1 (", pc1ra_lb, ")"),
    y = paste0("PC2 (", pc2ra_lb, ")"),
    color = "Decay stage"
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "right")

leg_A <- cowplot::get_legend(p_pc12_dw + theme(legend.position = "right"))
leg_B_decay_only <- cowplot::get_legend(p_pc12_dec + guides(shape = "none") + theme(legend.position = "right"))

leg_mid <- cowplot::plot_grid(
  leg_A,
  leg_B_decay_only,
  ncol = 1,
  rel_heights = c(1, 0.38),
  align = "v"
)

pA_noleg <- p_pc12_dw + theme(legend.position = "none")
pB_noleg <- p_pc12_dec + theme(legend.position = "none")

panel_pc12 <- cowplot::plot_grid(
  pA_noleg,
  leg_mid,
  pB_noleg,
  ncol = 3,
  rel_widths = c(1, 0.46, 1),
  labels = c("A", "", "B")
)

print(panel_pc12)
ggsave(
  filename = file.path(outdir_ord, "PCA_PANEL_FINAL_3taxa.png"),
  plot = panel_pc12,
  width = 21,
  height = 10.5,
  units = "cm",
  dpi = 900
)

# =============================================================================
# FIGURE 5
# =============================================================================
df <- tibble::tribble(
  ~Scope,   ~Metric,  ~R2_obs, ~R2_null_mean, ~R2_null_CI,      ~p,       ~R2_SES,
  "LOG",    "r.Ait",    0.347,      0.175,     "0.120-0.238", "<0.001",   5.612,
  "LOG",    "Jaccard",  0.250,      0.162,     "0.153-0.171", "<0.001",  18.680,
  "SNAG",   "r.Ait",    0.595,      0.501,     "0.320-0.661", "0.162",    1.047,
  "SNAG",   "Jaccard",  0.551,      0.450,     "0.407-0.489", "<0.001",   5.052,
  "aFWD",   "r.Ait",    0.658,      0.574,     "0.427-0.697", "0.120",    1.207,
  "aFWD",   "Jaccard",  0.567,      0.534,     "0.506-0.559", "0.004",    2.190,
  "fFWD",   "r.Ait",    0.068,      0.049,     "0.011-0.109", "0.218",    0.713,
  "fFWD",   "Jaccard",  0.098,      0.058,     "0.048-0.069", "<0.001",   7.445
) %>%
  mutate(
    Metric = case_when(
      Metric %in% c("r.Ait", "RA", "robust.aitchison") ~ "robust.aitchison",
      Metric %in% c("Jaccard", "JAC", "jaccard") ~ "jaccard",
      TRUE ~ Metric
    ),
    Metric = tolower(Metric)
  )

p_fig5_sesbar <- ggplot(df, aes(x = Scope, y = R2_SES, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 1.645, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_continuous(trans = "log10") +
  labs(
    y = "Standardized effect size of R² (SES, log10)",
    x = "Scope",
    fill = "Distance metric",
    title = "Observed SES across deadwood scopes and distance metrics"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "top"
  )

print(p_fig5_sesbar)

df2 <- df %>%
  mutate(
    R2_null_CI = str_replace_all(R2_null_CI, "–", "-"),
    CI_low  = as.numeric(str_split_fixed(R2_null_CI, "-", 2)[, 1]),
    CI_high = as.numeric(str_split_fixed(R2_null_CI, "-", 2)[, 2]),
    sd_null = (CI_high - CI_low) / (2 * 1.96),
    SES = R2_SES,
    SES_low  = (R2_obs - CI_high) / sd_null,
    SES_high = (R2_obs - CI_low)  / sd_null
  )

p_fig5_sesci <- ggplot(df2, aes(y = Scope, x = SES, color = Metric)) +
  geom_errorbarh(
    aes(xmin = SES_low, xmax = SES_high),
    height = 0.25,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(size = 2.8, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 1.645, linetype = "dashed") +
  labs(
    x = "R² SES (z scale)",
    y = "Scope",
    title = "SES of tree identity by scope and distance metric"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p_fig5_sesci)

log_R2_null_Ait <- readRDS(file.path("plots", "INTRA_diam", "LOG", "R2_null_Ait_LOG.rds"))
log_R2_null_Jac <- readRDS(file.path("plots", "INTRA_diam", "LOG", "R2_null_Jac_LOG.rds"))
snag_R2_null_Ait <- readRDS(file.path("plots", "SNAG", "snag_R2_null_Ait.RDS"))
snag_R2_null_Jac <- readRDS(file.path("plots", "SNAG", "snag_R2_null_Jac.RDS"))
aFWD_R2_null_Ait <- readRDS(file.path("plots", "FWD", "ATTACHED", "aFWD_R2_null_Ait.rds"))
aFWD_R2_null_Jac <- readRDS(file.path("plots", "FWD", "ATTACHED", "aFWD_R2_null_Jac.rds"))
fFWD_R2_null <- read_csv(file.path("tables", "FWD", "FALLEN", "NATMAN_NULL_R2.csv"), show_col_types = FALSE)

get_or_na <- function(x) {
  obj <- get0(x, ifnotfound = NA)
  if (is.data.frame(obj) || is_tibble(obj)) stop("get_or_na is for atomic vectors only")
  as.numeric(obj)
}

null_R2 <- list(
  LOG = list(
    ait = get_or_na("log_R2_null_Ait"),
    jac = get_or_na("log_R2_null_Jac")
  ),
  SNAG = list(
    ait = get_or_na("snag_R2_null_Ait"),
    jac = get_or_na("snag_R2_null_Jac")
  ),
  aFWD = list(
    ait = get_or_na("aFWD_R2_null_Ait"),
    jac = get_or_na("aFWD_R2_null_Jac")
  ),
  fFWD = list(
    ait = tryCatch(as.numeric(fFWD_R2_null$robust_aitchison), error = function(e) NA_real_),
    jac = tryCatch(as.numeric(fFWD_R2_null$jaccard), error = function(e) NA_real_)
  )
)

len_tab <- bind_rows(lapply(names(null_R2), function(s) {
  tibble(
    Scope = s,
    n_ait = length(null_R2[[s]]$ait[!is.na(null_R2[[s]]$ait)]),
    n_jac = length(null_R2[[s]]$jac[!is.na(null_R2[[s]]$jac)])
  )
}))
print(len_tab)

null_R2_df <- bind_rows(lapply(names(null_R2), function(s) {
  bind_rows(
    tibble(Scope = s, Metric = "robust.aitchison", R2 = as.numeric(null_R2[[s]]$ait)),
    tibble(Scope = s, Metric = "jaccard",          R2 = as.numeric(null_R2[[s]]$jac))
  )
})) %>%
  filter(!is.na(R2)) %>%
  mutate(Metric = tolower(Metric))

scope_order <- c("LOG", "SNAG", "aFWD", "fFWD")

plot_data <- null_R2_df %>%
  left_join(
    df %>% select(Scope, Metric, R2_obs, R2_null_mean),
    by = c("Scope", "Metric")
  ) %>%
  mutate(
    Scope  = factor(Scope, levels = scope_order, ordered = TRUE),
    Metric = factor(Metric, levels = c("robust.aitchison", "jaccard"))
  ) %>%
  arrange(Metric, Scope)

obs_points <- df %>%
  select(Scope, Metric, R2_obs) %>%
  mutate(
    Scope  = factor(Scope, levels = scope_order, ordered = TRUE),
    Metric = factor(Metric, levels = c("robust.aitchison", "jaccard"))
  )

ci_ticks <- plot_data %>%
  group_by(Scope, Metric) %>%
  summarise(
    q025 = quantile(R2, 0.025, na.rm = TRUE),
    q975 = quantile(R2, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

y_map <- tibble(
  Scope = factor(scope_order, levels = scope_order, ordered = TRUE),
  y = seq_along(scope_order)
)

plot_data  <- plot_data  %>% left_join(y_map, by = "Scope")
obs_points <- obs_points %>% left_join(y_map, by = "Scope")
ci_ticks   <- ci_ticks   %>% left_join(y_map, by = "Scope")

p_twopanel <- ggplot(plot_data, aes(x = R2, y = y, fill = Scope)) +
  ggridges::geom_density_ridges(
    alpha = 0.65,
    scale = 1,
    rel_min_height = 0.005,
    color = "grey15"
  ) +
  geom_point(
    data = obs_points,
    aes(x = R2_obs, y = y),
    color = "black",
    fill = "black",
    shape = 23,
    size = 2.8,
    alpha = 0.8,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = ci_ticks,
    aes(x = q025, xend = q025, y = y - 0.08, yend = y + 0.08),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  geom_segment(
    data = ci_ticks,
    aes(x = q975, xend = q975, y = y - 0.08, yend = y + 0.08),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.3
  ) +
  facet_wrap(~ Metric, ncol = 2) +
  scale_fill_manual(values = dw_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = y_map$y,
    labels = as.character(y_map$Scope),
    expand = expansion(add = c(0.3, 0.3))
  ) +
  scale_x_continuous(limits = c(-0.02, 0.75)) +
  labs(
    x = expression(R^2),
    y = "Scope",
    title = "Null distributions with 95% CI and observed R²"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_twopanel)
ggsave(
  file.path("plots", "ALL", "fig_R2_twopanel.png"),
  p_twopanel,
  width = 6.7,
  height = 4.5,
  units = "in",
  dpi = 900
)

