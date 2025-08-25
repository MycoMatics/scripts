## ===========================================================
## SCRIPT 3B species lists and overlaps
## ALL TAXA SUMMARIES + EULER OVERLAPS (clean, modular block)
## - dw_type2 (LOG, aFWD, fFWD, SNAG, SOIL)
## - ds_at_drill
## - size
## Includes: UNITE links, per-level TSVs, headless-safe saving
## ===========================================================

suppressPackageStartupMessages({
  library(tidyverse); library(janitor); library(eulerr); library(ggplot2);library(VennDiagram)
})

stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"), exists("threshold"))

## --- Directories ------------------------------------------------------
dir.create("tables/ALL",    recursive = TRUE, showWarnings = FALSE)
dir.create("plots/OVERLAP", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/OVERLAP", recursive = TRUE, showWarnings = FALSE)

## --- 0) Helpers -------------------------------------------------------
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

save_euler <- function(sets, title, outfile_base) {
  if (!"grid" %in% .packages()) suppressPackageStartupMessages(library(grid))
  fit <- eulerr::euler(sets)
  p <- plot(fit, quantities = TRUE, labels = TRUE, main = title)  # S3 dispatch -> plot.euler
  out <- paste0(outfile_base, ".png")
  
  # Headless-safe device
  open_png <- function(filename, w = 8, h = 6, dpi = 300) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(filename, width = w, height = h, units = "in", res = dpi); return(invisible(TRUE))
    }
    if (capabilities("cairo")) {
      grDevices::png(filename, width = w, height = h, units = "in", res = dpi, type = "cairo"); return(invisible(TRUE))
    }
    grDevices::png(filename, width = w * dpi, height = h * dpi, res = dpi); invisible(TRUE)
  }
  
  if (inherits(p, "ggplot")) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(out, plot = p, device = ragg::agg_png, width = 8, height = 6, units = "in", dpi = 300)
    } else if (capabilities("cairo")) {
      ggplot2::ggsave(out, plot = p, device = "png",      width = 8, height = 6, units = "in", dpi = 300)
    } else {
      open_png(out); print(p); dev.off()
    }
  } else {
    open_png(out); grid::grid.newpage(); grid::grid.draw(p); dev.off()
  }
  print(p)
  invisible(TRUE)
}

write_set_tables <- function(sets, outfile_stem) {
  sizes <- tibble(group = names(sets), n_unique_sh = vapply(sets, length, integer(1)))
  jac   <- pairwise_jaccard(sets)
  readr::write_tsv(sizes, file.path("tables/OVERLAP", paste0(outfile_stem, "_set_sizes_threshold", threshold, ".tsv")))
  readr::write_tsv(jac,   file.path("tables/OVERLAP", paste0(outfile_stem, "_pairwise_jaccard_threshold", threshold, ".tsv")))
}

## --- 1) Build dw_type2 metadata & long table -------------------------
meta_all <- META1 %>%
  filter(!is.na(dw_type)) %>%
  mutate(
    dw_type2 = dplyr::case_when(
      dw_type == "FWD" & position %in% c("ATTACHED","AFWD1","AFWD2","CROWN") ~ "aFWD",
      dw_type == "FWD" & position %in% c("FALLEN")                            ~ "fFWD",
      TRUE                                                                     ~ as.character(dw_type)
    ),
    dw_type2 = factor(dw_type2, levels = c("LOG","aFWD","fFWD","SNAG","SOIL"))
  )
otu_all <- otu_matrix_filt[meta_all$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_all), meta_all$sample))

otu_long_all <- otu_all %>%
  as.data.frame() %>% rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "sh_code", values_to = "abundance") %>%
  dplyr::filter(abundance > 0) %>%
  dplyr::left_join(dplyr::select(meta_all, sample, dw_type2), by = "sample") %>%
  dplyr::filter(!is.na(dw_type2)) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(rel_abund_sample = abundance / sum(abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(tax, sh_code, species, genus), by = "sh_code")

## --- 2) ALL taxa by dead‑wood type -----------------------------------
all_species_by_dw <- summarise_alltaxa_by(otu_long_all, dw_type2)

combined_path_dw <- file.path("tables/ALL", paste0("ALL_species_by_dwtype_threshold", threshold, ".tsv"))
readr::write_tsv(all_species_by_dw, combined_path_dw)

for (dw in levels(meta_all$dw_type2)) {
  tab_dw <- dplyr::filter(all_species_by_dw, dw_type2 == dw)
  if (nrow(tab_dw) > 0) {
    outdir_dw <- file.path("tables/ALL")
    readr::write_tsv(tab_dw, file.path(outdir_dw, paste0("ALL_species_", dw, "_threshold", threshold, ".tsv")))
  }
}

## --- 3) ALL taxa by ds_at_drill & by size (across ALL samples) -------
meta_all2 <- META1 %>%
  filter(!is.na(sample)) %>%
  mutate(ds_at_drill = factor(ds_at_drill, levels = c("0","1","2","3","4","5"))) %>%
  select(sample, ds_at_drill, size)

otu_all2 <- otu_matrix_filt[meta_all2$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_all2), meta_all2$sample))

otu_long_all2 <- otu_all2 %>%
  as.data.frame() %>% rownames_to_column("sample") %>%
  tidyr::pivot_longer(-sample, names_to = "sh_code", values_to = "abundance") %>%
  dplyr::filter(abundance > 0) %>%
  dplyr::left_join(meta_all2, by = "sample") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(rel_abund_sample = abundance / sum(abundance)) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(dplyr::select(tax, sh_code, species, genus), by = "sh_code")

all_species_by_ds   <- summarise_alltaxa_by(otu_long_all2, ds_at_drill)
all_species_by_size <- summarise_alltaxa_by(otu_long_all2, size)

path_ds_combined   <- file.path("tables/ALL", paste0("ALL_species_by_ds_at_drill_threshold", threshold, ".tsv"))
path_size_combined <- file.path("tables/ALL", paste0("ALL_species_by_size_threshold", threshold, ".tsv"))
readr::write_tsv(all_species_by_ds,   path_ds_combined)
readr::write_tsv(all_species_by_size, path_size_combined)

write_per_level(all_species_by_ds,   ds_at_drill, "tables/ALL", "ALL_species_ds")
write_per_level(all_species_by_size, size,        "tables/ALL", "ALL_species_size")

## --- 4) Euler overlaps (full categories) ------------------------------
sets_dw <- sets_from_summary(all_species_by_dw,   "dw_type2")
sets_ds <- sets_from_summary(all_species_by_ds,   "ds_at_drill")
sets_sz <- sets_from_summary(all_species_by_size, "size")

# Drop empty sets
nz <- function(x) x[vapply(x, length, integer(1)) > 0]
sets_dw <- nz(sets_dw); sets_ds <- nz(sets_ds); sets_sz <- nz(sets_sz)

save_euler(sets_dw, paste0("SH overlap by dead‑wood type (dw_type2) — threshold ", threshold),
           file.path("plots/OVERLAP", paste0("euler_dw_type2_threshold", threshold)))
write_set_tables(sets_dw, "dw_type2")

save_euler(sets_ds, paste0("SH overlap by ds_at_drill — threshold ", threshold),
           file.path("plots/OVERLAP", paste0("euler_ds_at_drill_threshold", threshold)))
write_set_tables(sets_ds, "ds_at_drill")

save_euler(sets_sz, paste0("SH overlap by size — threshold ", threshold),
           file.path("plots/OVERLAP", paste0("euler_size_threshold", threshold)))
write_set_tables(sets_sz, "size")

## --- 5) Reduced-category Euler only ----------------------------------
# (a) ds_at_drill: drop "0"
sets_ds_no0 <- sets_ds[setdiff(names(sets_ds), "0")] %>% nz()
save_euler(sets_ds_no0, paste0("SH overlap by ds_at_drill (excluding 0) — threshold ", threshold),
           file.path("plots/OVERLAP", paste0("euler_ds_at_drill_no0_threshold", threshold)))
write_set_tables(sets_ds_no0, "ds_at_drill_no0")

# (b) dw_type2: merge aFWD + fFWD → FWD
getset <- function(x, nm) if (nm %in% names(x)) unique(x[[nm]]) else character(0)
sets_dw_merged <- list(
  LOG  = getset(sets_dw, "LOG"),
  FWD  = unique(c(getset(sets_dw, "aFWD"), getset(sets_dw, "fFWD"))),
  SNAG = getset(sets_dw, "SNAG"),
  SOIL = getset(sets_dw, "SOIL")
) %>% nz()
save_euler(sets_dw_merged, paste0("SH overlap by dead‑wood type (FWD merged) — threshold ", threshold),
           file.path("plots/OVERLAP", paste0("euler_dwtype_FWDmerged_threshold", threshold)))
write_set_tables(sets_dw_merged, "dw_type2_FWDmerged")

## --- 6) Console peek --------------------------------------------------
message("Top 10 by dw_type2:"); print(all_species_by_dw %>% group_by(dw_type2) %>% slice_head(n = 10), n = 50)
message("Top 10 by ds_at_drill:"); print(all_species_by_ds %>% group_by(ds_at_drill) %>% slice_head(n = 10), n = 50)
message("Top 10 by size:"); print(all_species_by_size %>% group_by(size) %>% slice_head(n = 10), n = 50)

cat("Wrote:\n - ", combined_path_dw,
    "\n - ", path_ds_combined,
    "\n - ", path_size_combined,
    "\n - Overlaps in plots/OVERLAP/ and QC tables in tables/OVERLAP/\n", sep = "")

## -----------------------------------------------------------
## Per-sample taxa > 2% rel. abundance + prevalence + UNITE
## Requires: META1, otu_matrix_filt, tax
## -----------------------------------------------------------

stopifnot(exists("otu_matrix_filt"), exists("META1"), exists("tax"))

extract_top_taxa_per_sample <- function(otu, meta = NULL, tax_tbl, rel_cutoff = 0.02) {
  # Align to meta if provided
  if (!is.null(meta)) {
    otu <- otu[meta$sample, , drop = FALSE]
    stopifnot(identical(rownames(otu), meta$sample))
  }
  
  # 1) Prevalence across ALL samples (presence/absence)
  pa <- (otu > 0) * 1L
  prev_counts <- colSums(pa)
  prev_df <- tibble(
    sh_code        = names(prev_counts),
    prevalence_n   = as.integer(prev_counts),
    prevalence_prop = if (nrow(otu) > 0) prev_counts / nrow(otu) else NA_real_
  )
  
  # 2) Within-sample relative abundance
  rs <- pmax(rowSums(otu), 1)
  rel <- sweep(otu, 1, rs, "/")
  
  # 3) Long format & filter by cutoff
  rel_long <- rel %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::pivot_longer(-sample, names_to = "sh_code", values_to = "rel_abund_sample") %>%
    dplyr::filter(rel_abund_sample >= rel_cutoff)
  
  # 4) Join prevalence + taxonomy + UNITE link
  rel_long <- rel_long %>%
    dplyr::left_join(prev_df,   by = "sh_code") %>%
    dplyr::left_join(dplyr::select(tax_tbl, sh_code, species, genus), by = "sh_code") %>%
    dplyr::mutate(
      unite_link = paste0("https://unite.ut.ee/bl_forw_sh.php?sh_name=", sh_code, "#fndtn-panel1")
    ) %>%
    dplyr::arrange(sample, dplyr::desc(rel_abund_sample))
  
  # 5) Split to a named list of tibbles (one per sample)
  out_list <- split(rel_long, rel_long$sample)
  # Ensure empty samples still return empty tibbles (if any)
  all_samps <- rownames(otu)
  missing_samples <- setdiff(all_samps, names(out_list))
  if (length(missing_samples)) {
    empty_tbl <- tibble(sample = character(), sh_code = character(),
                        rel_abund_sample = numeric(), prevalence_n = integer(),
                        prevalence_prop = numeric(), species = character(),
                        genus = character(), unite_link = character())
    for (s in missing_samples) out_list[[s]] <- empty_tbl
  }
  
  # Order list by sample order in OTU
  out_list <- out_list[all_samps]
  return(out_list)
}
## --- Prepare OTU + META for per-sample 2% cutoff ---------------------
meta_no_fruit <- META1 %>%
  filter(!is.na(sample), dw_type != "FRUITBODY")

otu_no_fruit <- otu_matrix_filt[meta_no_fruit$sample, , drop = FALSE]
stopifnot(identical(rownames(otu_no_fruit), meta_no_fruit$sample))

## --- Run the extractor (defined earlier) ---
per_sample_top2pct <- extract_top_taxa_per_sample(
  otu       = otu_no_fruit,
  meta      = meta_no_fruit,
  tax_tbl   = tax,
  rel_cutoff = 0.02
)


## Peek at a couple of samples:
names(per_sample_top2pct)[1:5]
per_sample_top2pct[[ names(per_sample_top2pct)[1] ]] %>% print(n = 20)

## ====================================================
## Master table of all taxa with ≥ 2% rel. abund in any sample
## ====================================================

all_top2pct <- bind_rows(per_sample_top2pct, .id = "sample")

# collapse to unique taxa
all_unique_top2pct <- all_top2pct %>%
  dplyr::distinct(sh_code, species, genus, unite_link) %>%
  dplyr::left_join(
    all_top2pct %>%
      dplyr::group_by(sh_code) %>%
      dplyr::summarise(
        n_samples   = dplyr::n_distinct(sample),
        prop_samples = n_samples / length(per_sample_top2pct),
        .groups = "drop"
      ),
    by = "sh_code"
  ) %>%
  dplyr::arrange(desc(n_samples))

# Preview
print(all_unique_top2pct, n = 20)

# Write out
outpath <- file.path("tables/ALL", paste0("ALL_unique_top2pct_taxa_threshold", threshold, ".tsv"))
readr::write_tsv(all_unique_top2pct, outpath)
cat("Wrote unique-top2pct taxa table to:", outpath, "\n")



# =============================================
# Indicator analysis (indicspecies::multipatt)
# - variables: dw_type2, ds_at_drill, size
# - blocks: natman
# - methods: IndVal.g and r.g
# - excludes FRUITBODY samples
# =============================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(indicspecies)
  library(permute)
  library(janitor)
})

stopifnot(exists("META1"), exists("otu_matrix_filt"), exists("tax"), exists("threshold"))

dir.create("tables/INDICATORS", recursive = TRUE, showWarnings = FALSE)

# --- 0) Build metadata, exclude FRUITBODY, derive dw_type2 -------------
meta0 <- META1 %>%
  filter(!is.na(sample), !(dw_type %in% c("FRUITBODY", "POSITIVE", "NEGATIVE", "SOIL"))) %>% 
  mutate(
    # Split FWD into attached vs fallen (harmonized with earlier code)
    dw_type2 = case_when(
      dw_type == "FWD" & position %in% c("ATTACHED","AFWD1","AFWD2","CROWN") ~ "aFWD",
      dw_type == "FWD" & position %in% c("FALLEN")                            ~ "fFWD",
      TRUE                                                                     ~ as.character(dw_type)
    ),
    dw_type2 = factor(dw_type2, levels = c("LOG","aFWD","fFWD","SNAG","SOIL")),
    ds_at_drill = factor(ds_at_drill, levels = c("0","1","2","3","4","5"))
  )

# align OTU
X_raw <- otu_matrix_filt[meta0$sample, , drop = FALSE]
stopifnot(identical(rownames(X_raw), meta0$sample))

# --- 1) core runner ----------------------------------------------------
run_multipatt_blocked <- function(otu, meta, var, func = c("IndVal.g","r.g"),
                                  target_sum = 3000, nperm = 999, alpha = 0.05) {
  func <- match.arg(func)
  
  # select complete cases for this variable + natman
  meta_var <- meta %>% filter(!is.na(.data[[var]]), !is.na(natman))
  if (nrow(meta_var) < 3) return(tibble())  # too few
  
  X <- otu[meta_var$sample, , drop = FALSE]
  stopifnot(identical(rownames(X), meta_var$sample))
  
  # normalize per sample to target_sum (keeps relative info for IndVal.g)
  rs <- pmax(rowSums(X), 1)
  Xn <- sweep(X, 1, rs, "/") * target_sum
  
  group <- droplevels(factor(meta_var[[var]]))
  if (nlevels(group) < 2L) return(tibble())
  
  # blocked permutations by natman
  ctrl <- how(nperm = nperm)
  setBlocks(ctrl) <- meta_var$natman
  
  # run
  res <- indicspecies::multipatt(Xn, group, func = func, control = ctrl, duleg = T)
  
  # ---- extract tidy results ----
  sign_df <- as.data.frame(res$sign)
  if (!nrow(sign_df)) return(tibble())
  sign_df$sh_code <- rownames(sign_df)
  
  # Build readable group labels
  gcols <- grep("^s\\.", names(sign_df), value = TRUE)
  if (length(gcols)) {
    # derive labels from s.<level> columns
    levs <- sub("^s\\.", "", gcols)
    idx  <- as.matrix(sign_df[, gcols, drop = FALSE]) > 0
    group_label <- apply(idx, 1, function(r) {
      paste(levs[r], collapse = "+")
    })
  } else if (!is.null(res$comb) && nrow(res$comb) > 0 && "index" %in% names(sign_df)) {
    # map via comb matrix if present
    glv <- levels(group)
    comb_labels <- apply(res$comb, 1, function(row) paste(glv[as.logical(row)], collapse = "+"))
    comb_map <- tibble(index = seq_len(nrow(res$comb)), group_label = comb_labels)
    group_label <- comb_map$group_label[match(sign_df$index, comb_map$index)]
  } else {
    group_label <- NA_character_
  }
  
  out <- tibble(
    var              = var,
    method           = func,
    sh_code          = sign_df$sh_code,
    stat             = sign_df$stat,
    p_value          = sign_df$p.value,
    group_label      = group_label
  ) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    filter(is.finite(p_adj)) %>%
    arrange(p_adj, desc(stat)) %>%
    left_join(select(tax, sh_code, genus, species), by = "sh_code") %>%
    mutate(unite_link = paste0("https://unite.ut.ee/bl_forw_sh.php?sh_name=", sh_code, "#fndtn-panel1"))
}

# --- 2) run for each variable & method --------------------------------
vars <- c("dw_type2","ds_at_drill","size")
methods <- c("IndVal.g","r.g")

all_res <- list()

for (v in vars) {
  for (m in methods) {
    cat("Running multipatt for", v, "with", m, "...\n")
    res_vm <- run_multipatt_blocked(
      otu = X_raw, meta = meta0, var = v, func = m,
      target_sum = 3000, nperm = 999, alpha = 0.05
    )
    if (nrow(res_vm)) {
      # write TSV
      fn <- file.path("tables/INDICATORS",
                      paste0("indic_", v, "_", m, "_threshold", threshold, ".tsv"))
      readr::write_tsv(res_vm, fn)
      cat("  -> wrote ", fn, " (", nrow(res_vm), " rows)\n", sep = "")
      all_res[[paste(v, m, sep = "_")]] <- res_vm
    } else {
      cat("  -> no results for ", v, " / ", m, "\n", sep = "")
    }
  }
}

# --- 3) combined master table -----------------------------------------
if (length(all_res)) {
  indic_master <- bind_rows(all_res, .id = "block") %>%
    select(var, method, sh_code, genus, species, group_label, stat, p_value, p_adj, unite_link) %>%
    arrange(var, method, p_adj, desc(stat))
  out_master <- file.path("tables/INDICATORS", paste0("indic_MASTER_threshold", threshold, ".tsv"))
  readr::write_tsv(indic_master, out_master)
  cat("Master table:", out_master, " (", nrow(indic_master), " rows )\n", sep = "")
} else {
  cat("No indicator results to combine.\n")
}
indic_master <- indic_master %>% filter(p_value<=0.009) %>% filter(stat>0.6)
write_csv2(indic_master, paste0("tables/INDICATORS/indic_master_", threshold, ".csv"))
print(indic_master, n=Inf)

# --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(VennDiagram)
})

# Headless-safe PNG writer (reuses ragg if present)
ggsave_safe <- function(filename, plot, w=8, h=6, dpi=300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(filename, plot = plot, device = ragg::agg_png,
                    width = w, height = h, units = "in", dpi = dpi)
  } else if (capabilities("cairo")) {
    ggplot2::ggsave(filename, plot = plot, device = "png",
                    width = w, height = h, units = "in", dpi = dpi)
  } else {
    ggplot2::ggsave(filename, plot = plot, width = w, height = h, units = "in", dpi = dpi)
  }
}

# Convert a named list of sets (character vectors of SH codes) to a binary incidence data frame
# rows = SH codes, columns = groups (0/1 presence)
incidence_from_sets <- function(sets) {
  stopifnot(length(sets) > 0, !is.null(names(sets)))
  all_items <- unique(unlist(sets, use.names = FALSE))
  if (length(all_items) == 0) return(tibble())
  mat <- sapply(sets, function(s) as.integer(all_items %in% s))
  df  <- as_tibble(mat) %>% mutate(sh_code = all_items) %>% relocate(sh_code)
  df
}

# Build a ComplexUpset UpSet plot (ggplot-based).
# If ComplexUpset is not available, falls back to a simpler bar chart of singletons and top pairwise intersections.
save_upset <- function(sets, title, outfile_base, top_n = 25) {
  inc <- incidence_from_sets(sets)
  if (!nrow(inc)) return(invisible(NULL))
  
  out_png <- paste0(outfile_base, "_upset.png")
  
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    # ComplexUpset expects logical columns for sets
    df <- inc %>% mutate(across(-sh_code, \(x) x > 0))
    set_names <- names(df)[names(df) != "sh_code"]
    
    # Keep only the top_n largest intersections to avoid overplotting
    # ComplexUpset can sort by size; we can pass "intersections" argument = NULL and control via min_size if needed.
    p <- ComplexUpset::upset(
      df,
      set_names,
      name = "SH intersection size",
      width_ratio = 0.2,
      min_size = 1,
      themes = ComplexUpset::upset_modify_themes(
        list("intersections_matrix" = theme(axis.text.x = element_text(angle = 0)))
      )
    ) +
      ggtitle(title)
    
    ggsave_safe(out_png, p, w=10, h=6, dpi=300)
    print(p)
  } else {
    message("ComplexUpset not installed; writing a simplified ggplot summary instead.")
    # Fallback: show (i) per-set uniques; (ii) top pairwise intersections
    # (a) uniques per set
    uniques <- map_int(names(sets), function(nm) {
      length(setdiff(sets[[nm]], unique(unlist(sets[names(sets) != nm], use.names = FALSE))))
    })
    df_uni <- tibble(group = names(sets), n_unique = uniques)
    
    # (b) top pairwise intersections
    pw <- expand_grid(a = names(sets), b = names(sets)) %>%
      filter(a < b) %>%
      mutate(n_inter = map2_int(a, b, ~length(intersect(sets[[.x]], sets[[.y]])))) %>%
      arrange(desc(n_inter)) %>% slice_head(n = min(top_n, n()))
    
    p1 <- ggplot(df_uni, aes(x = reorder(group, n_unique), y = n_unique)) +
      geom_col() +
      coord_flip() +
      labs(x = NULL, y = "Unique SH count", title = paste0(title, " — uniques per group"))
    
    p2 <- ggplot(pw, aes(x = reorder(paste(a, b, sep = " ∩ "), n_inter), y = n_inter)) +
      geom_col() +
      coord_flip() +
      labs(x = "Pairwise intersection", y = "SH count",
           title = paste0(title, " — top pairwise intersections"))
    
    ggsave_safe(paste0(outfile_base, "_fallback_uniques.png"), p1, w=7, h=4, dpi=300)
    ggsave_safe(paste0(outfile_base, "_fallback_pairs.png"),  p2, w=7, h=6, dpi=300)
    print(p1); print(p2)
  }
  invisible(TRUE)
}
# Jaccard heatmap with in-cell labels
save_jaccard_heatmap <- function(sets, title, outfile_base) {
  if (length(sets) == 0) return(invisible(NULL))
  jac <- pairwise_jaccard(sets) %>%
    pivot_longer(-group, names_to = "group2", values_to = "jaccard") %>%
    filter(!is.na(jaccard))
  
  # order axes by mean similarity for stable layout
  ord <- jac %>% group_by(group) %>% summarize(m = mean(jaccard, na.rm=TRUE)) %>% arrange(desc(m)) %>% pull(group)
  jac$group  <- factor(jac$group,  levels = ord)
  jac$group2 <- factor(jac$group2, levels = ord)
  
  p <- ggplot(jac, aes(group2, group, fill = jaccard)) +
    geom_tile() +
    geom_text(aes(label = scales::percent(jaccard, accuracy = 1))) +
    scale_fill_viridis_c(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "Jaccard", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  out_png <- paste0(outfile_base, "_jaccard_heatmap.png")
  ggsave_safe(out_png, p, w=7, h=6, dpi=300)
  print(p)
  invisible(TRUE)
}

# (b) dw_type2 with FWD merged
getset <- function(x, nm) if (nm %in% names(x)) unique(x[[nm]]) else character(0)
sets_dw_merged <- list(
  LOG  = getset(sets_dw, "LOG"),
  FWD  = unique(c(getset(sets_dw, "aFWD"), getset(sets_dw, "fFWD"))),
  SNAG = getset(sets_dw, "SNAG"),
  SOIL = getset(sets_dw, "SOIL")
)
sets_dw_merged <- sets_dw_merged[vapply(sets_dw_merged, length, 1L) > 0]

save_upset(
  sets_dw_merged,
  paste0("SH intersections by dead-wood type (FWD merged) — threshold ", threshold),
  file.path("plots/OVERLAP", paste0("dwtype_FWDmerged_threshold", threshold))
)
save_jaccard_heatmap(
  sets_dw_merged,
  paste0("Pairwise Jaccard (dw_type2, FWD merged) — threshold ", threshold),
  file.path("plots/OVERLAP", paste0("dwtype_FWDmerged_threshold", threshold))
)
#### 4. Exploration: Species by Deadwood Type ####
# Transform OTU table to long format with metadata
dw_keep <- c("FWD", "LOG", "SNAG", "SOIL" )
otu_long2 <- otu_matrix_filt %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "n_seq") %>%
  filter(n_seq > 0) %>%
  left_join(select(META1, sample, dw_type, position), by = "sample") %>%
  left_join(select(tax, sh_code, species, genus), by = "sh_code") %>%
  filter(dw_type %in% dw_keep) %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = n_seq / sum(n_seq)) %>%
  ungroup()

# Compute total number of reads across all samples
total_reads <- sum(otu_long2$n_seq)

# Compute total reads per DW type
total_reads_by_dwtype <- otu_long2 %>%
  group_by(dw_type) %>%
  summarise(total_reads_dw = sum(n_seq), .groups = "drop")


#### 5. Venn Diagram: Species Overlap ####
dw_keep <- c("FWD", "LOG", "SNAG" )
str(otu_matrix_filt)
otu_long2 <- otu_matrix_filt %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  pivot_longer(-sample, names_to = "sh_code", values_to = "n_seq") %>%
  filter(n_seq > 0) %>%
  left_join(select(META1, sample, dw_type, position), by = "sample") %>%
  left_join(select(tax, sh_code, species, genus), by = "sh_code") %>%
  filter(dw_type %in% dw_keep) %>%
  group_by(sample) %>%
  mutate(rel_abund_sample = n_seq / sum(n_seq)) %>%
  ungroup()
species_list_by_dwtype <- otu_long2 %>%
  filter(!is.na(dw_type)) %>%
  mutate(
    dw_type2 = dplyr::case_when(
      dw_type == "FWD" & position %in% c("ATTACHED","AFWD1","AFWD2","CROWN") ~ "aFWD",
      dw_type == "FWD" & position %in% c("FALLEN")                            ~ "fFWD",
      TRUE                                                                     ~ as.character(dw_type)
    ),
    dw_type2 = factor(dw_type2, levels = c("LOG","aFWD","fFWD","SNAG"))
  ) %>% 
  group_by(dw_type2) %>%
  summarise(sh_codes = list(unique(sh_code)), .groups = "drop") %>%
  deframe()  # Named list

venn_dw_types <- intersect(names(dw_colors), names(species_list_by_dwtype))
venn_colors <- dw_colors[venn_dw_types]
species_list_for_venn <- species_list_by_dwtype[venn_dw_types]

cat("Number of unique SHs per deadwood type:\n")
print(sapply(species_list_for_venn, length))
venn_plot <- venn.diagram(
  x = species_list_for_venn,
  category.names = venn_dw_types,
  filename = NULL,
  fill = venn_colors,
  col = venn_colors,
  alpha = 0.5,
  cex = 1,
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontfamily = "sans"
)
plotname <- paste0("plots/ALL/","venn_dw_type_all_species_",threshold, ".png")
ggsave(plotname, venn_plot, width = 4, height = 4, dpi = 600)
