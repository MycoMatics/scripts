# =============================================================================
# Script: Script3_rarefaction_collectors.R
# Project: Deadwood fungal community analyses
# Purpose: Generate rarefaction curves and iNEXT collector's curves for
#          deadwood fractions and woody-versus-soil contrasts
# Author: Glen Dierickx
# Created: October 2025
# Last updated: 2026-03-09
#
# Required in memory:
# - META1
# - otu_matrix_filt
# - threshold
#
# Outputs:
# - plots/rarefaction_by_dw_type2_filtered_<threshold>.png
# - plots/collectors_iNEXT_by_dw_type2_<threshold>.png
# - plots/collectors_iNEXT_woody_vs_soil_<threshold>.png
# =============================================================================

source("Scripts/Script0_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(iNEXT)
  library(janitor)
})

assert_objects(c("META1", "otu_matrix_filt", "threshold"))

set.seed(SEED_GLOBAL)

# ---- Paths and reproducibility -----------------------------------------------
data_dir <- "R_data"
ensure_dir(data_dir)
setwd(data_dir)

ensure_dir("plots")
ensure_dir("plots/LOG")

out_rare_dw  <- file.path("plots", paste0("rarefaction_by_dw_type2_filtered_", threshold, ".png"))
out_inext_dw <- file.path("plots", paste0("collectors_iNEXT_by_dw_type2_", threshold, ".png"))
out_inext_ws <- file.path("plots", paste0("collectors_iNEXT_woody_vs_soil_", threshold, ".png"))

# ---- Parameters ---------------------------------------------------------------
RARE_STEP       <- 1500
ENDPOINT_MULT_DW <- 2
ENDPOINT_CAP_DW  <- 170
ENDPOINT_MULT_WS <- 3
ENDPOINT_CAP_WS  <- 690

dw_keep2 <- c("SNAG", "LOG", "aFWD", "fFWD", "SOIL")

order_labs <- c(
  "0" = "q = 0\nRichness\n(Number of SHs)",
  "1" = "q = 1\nShannon diversity\n(Number of typical SHs)",
  "2" = "q = 2\nSimpson diversity\n(Number of dominant SHs)"
)

cat("You are using SH threshold ", threshold, "\n", sep = "")

# ---- Align samples -----------------------------------------------------------
otu_work <- otu_matrix_filt[rownames(otu_matrix_filt) %in% META1$sample, , drop = FALSE]
META1 <- META1 %>% filter(sample %in% rownames(otu_work))
sample_names <- rownames(otu_work)

# ------------------------------------------------------------------------------
# 1) Rarefaction curves by deadwood type
# ------------------------------------------------------------------------------
META1_afwd <- META1

pdf(NULL)
rarecurve_list2 <- vegan::rarecurve(
  otu_work,
  step = RARE_STEP,
  label = FALSE,
  return = "list"
)
dev.off()

rare_df2 <- purrr::map2_dfr(
  rarecurve_list2,
  rownames(otu_work),
  ~tibble(
    sample = .y,
    reads = attr(.x, "Subsample"),
    richness = as.integer(.x)
  )
) %>%
  left_join(META1_afwd %>% select(sample, dw_type2), by = "sample")

p_rare2 <- ggplot(rare_df2, aes(x = reads, y = richness, group = sample, color = dw_type2)) +
  geom_line(alpha = 0.7) +
  scale_color_manual(values = dw_colors, na.value = "grey70") +
  theme_minimal(base_size = 13) +
  labs(
    title = "Rarefaction curves by deadwood type (FWD split)",
    x = "Reads (subsampled)",
    y = "Observed SH richness",
    color = "Deadwood type"
  )

print(p_rare2)
ggsave(out_rare_dw, p_rare2, width = 6, height = 6, dpi = 900)
message("Wrote: ", out_rare_dw)

# ------------------------------------------------------------------------------
# 2) iNEXT by deadwood type
# ------------------------------------------------------------------------------
incidence_list2 <- lapply(dw_keep2, function(dw) {
  samples <- META1_afwd %>% filter(dw_type2 == dw) %>% pull(sample)
  if (length(samples) < 2) return(NULL)
  
  otus <- otu_work[rownames(otu_work) %in% samples, , drop = FALSE]
  otus <- otus[, colSums(otus) > 0, drop = FALSE]
  freq <- colSums(otus > 0)
  
  c(length(samples), freq)
})

names(incidence_list2) <- dw_keep2
incidence_list2 <- incidence_list2[!sapply(incidence_list2, is.null)]

endpoints2 <- sapply(
  incidence_list2,
  function(x) min(ceiling(x[1] * ENDPOINT_MULT_DW), ENDPOINT_CAP_DW)
)

iNEXT_out2 <- lapply(names(incidence_list2), function(dw) {
  out <- iNEXT(
    incidence_list2[[dw]],
    q = c(0, 1, 2),
    datatype = "incidence_freq",
    endpoint = endpoints2[dw]
  )
  as_tibble(out$iNextEst$size_based) %>% mutate(dw_type2 = dw)
})

iNEXT_df2 <- bind_rows(iNEXT_out2)
observed_n2 <- sapply(incidence_list2, function(x) x[1])

iNEXT_df2 <- iNEXT_df2 %>%
  mutate(interp_extrap = ifelse(t <= observed_n2[dw_type2], "interpolated", "extrapolated"))

interp_points2 <- iNEXT_df2 %>%
  group_by(dw_type2, Order.q) %>%
  filter(t == max(t[interp_extrap == "interpolated"])) %>%
  ungroup()

p_inext_dw2 <- ggplot(
  iNEXT_df2,
  aes(x = t, y = qD, color = dw_type2, fill = dw_type2, linetype = interp_extrap)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
  geom_point(
    data = interp_points2,
    aes(x = t, y = qD, color = dw_type2),
    inherit.aes = FALSE,
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  scale_color_manual(values = dw_colors) +
  scale_fill_manual(values = dw_colors) +
  scale_linetype_manual(values = c(interpolated = "solid", extrapolated = "dotted")) +
  facet_wrap(~Order.q, labeller = as_labeller(order_labs), ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "iNEXT (incidence) by deadwood type (FWD split)",
    x = "Number of samples (observed & extrapolated)",
    y = "Estimated diversity (qD)",
    color = "Deadwood type",
    fill = "Deadwood type",
    linetype = "Type"
  )

print(p_inext_dw2)
ggsave(out_inext_dw, p_inext_dw2, width = 6, height = 12, dpi = 900)
message("Wrote: ", out_inext_dw)

res <- compare_iNEXT_df_at_coverage(
  iNEXT_df2,
  group_col = "dw_type2",
  orders = c(0, 1, 2),
  target_cov = 1,
  allow_extrapolation = FALSE
)

print(arrange(res$Cstar, Order.q))
print(arrange(res$contrasts, Order.q, p_adj), n = Inf)
print(arrange(res$points, Order.q, group))

tab_iNEXT <- res$contrasts %>%
  arrange(Order.q, p_adj) %>%
  group_by(Order.q) %>%
  mutate(
    sig = ifelse(p_adj < 0.05, "*", ""),
    ratio_ci = sprintf("%.2f [%.2f–%.2f]", ratio, ratio_LCL, ratio_UCL),
    p_fmt = format.pval(p_adj, eps = 1e-4)
  ) %>%
  transmute(
    q = Order.q,
    Cstar,
    contrast = paste(group1, "/", group2),
    ratio_ci,
    p_adj = p_fmt,
    sig
  )

print(tab_iNEXT, n = Inf)

t_at_C <- res$points %>%
  select(Order.q, dw_type2 = group, t, SC, qD, Cstar)

print(t_at_C)

# ------------------------------------------------------------------------------
# 3) WOODY vs SOIL (iNEXT)
# ------------------------------------------------------------------------------
META1_ws <- META1_afwd %>%
  mutate(woody = ifelse(dw_type2 == "SOIL", "SOIL", "WOODY"))

incidence_combined <- list(
  WOODY = {
    s <- META1_ws %>% filter(woody == "WOODY") %>% pull(sample)
    ot <- otu_work[rownames(otu_work) %in% s, , drop = FALSE]
    ot <- ot[, colSums(ot) > 0, drop = FALSE]
    c(length(s), colSums(ot > 0))
  },
  SOIL = {
    s <- META1_ws %>% filter(woody == "SOIL") %>% pull(sample)
    ot <- otu_work[rownames(otu_work) %in% s, , drop = FALSE]
    ot <- ot[, colSums(ot) > 0, drop = FALSE]
    c(length(s), colSums(ot > 0))
  }
)

endpoints_combined <- sapply(
  incidence_combined,
  function(x) min(ceiling(x[1] * ENDPOINT_MULT_WS), ENDPOINT_CAP_WS)
)

iNEXT_out_combined <- lapply(names(incidence_combined), function(g) {
  out <- iNEXT(
    incidence_combined[[g]],
    q = c(0, 1, 2),
    datatype = "incidence_freq",
    endpoint = endpoints_combined[g]
  )
  as_tibble(out$iNextEst$size_based) %>% mutate(group = g)
})

iNEXT_df_combined <- bind_rows(iNEXT_out_combined)
obs_n_combined <- sapply(incidence_combined, function(x) x[1])

iNEXT_df_combined <- iNEXT_df_combined %>%
  mutate(interp_extrap = ifelse(t <= obs_n_combined[group], "interpolated", "extrapolated"))

interp_points_combined <- iNEXT_df_combined %>%
  group_by(group, Order.q) %>%
  filter(t == max(t[interp_extrap == "interpolated"])) %>%
  ungroup()

p_combined <- ggplot(
  iNEXT_df_combined,
  aes(x = t, y = qD, color = group, fill = group, linetype = interp_extrap)
) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha = 0.18, color = NA) +
  geom_point(
    data = interp_points_combined,
    aes(x = t, y = qD, color = group),
    inherit.aes = FALSE,
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_linetype_manual(values = c(interpolated = "solid", extrapolated = "dotted")) +
  facet_wrap(~Order.q, labeller = as_labeller(order_labs), ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 13) +
  labs(
    title = "iNEXT: WOODY vs SOIL",
    x = "Number of samples (observed & extrapolated)",
    y = "Estimated diversity (qD)",
    color = "Category",
    fill = "Category",
    linetype = "Type"
  )

print(p_combined)
ggsave(out_inext_ws, p_combined, width = 5, height = 12, dpi = 900)
message("Wrote: ", out_inext_ws)
