

# Shared helpers a

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ----- Threshold handling ------------------------------------------------------
thr_tag <- function(x) sprintf("%03d", round(as.numeric(x) * 10))

std_barcode <- function(x) {
  x %>% str_replace("^barcode", "") %>% str_replace_all("_", "") %>% str_replace("rawumis", "")
}

# Standardize thresholded column names to "<prefix>_(status|SH_code)_TAG"
rename_threshold_cols <- function(df, prefix) {
  old <- names(df); new <- old
  m <- str_match(old, "^status \\(([^)]+)\\)$")
  hit <- which(!is.na(m[,2]))
  if (length(hit)) new[hit] <- paste0(prefix, "_status_", thr_tag(m[hit,2]))
  m <- str_match(old, "^SH code \\(([^)]+)\\)$")
  hit <- which(!is.na(m[,2]))
  if (length(hit)) new[hit] <- paste0(prefix, "_SH_code_", thr_tag(m[hit,2]))
  new[old == "Similarity percentage"] <- paste0("Similarity_percentage_", prefix)
  names(df) <- new
  df
}

# Read delimited with headers preserved
read_tab_keep_names <- function(path, delim = "\t") {
  read.delim(path, sep = delim, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

# Robust read of matches file (empty-safe)
read_matches_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (file.info(path)$size == 0) return(NULL)
  read_tab_keep_names(path, "\t")
}

# From long SH counts to wide OTU
pivot_otu <- function(df_long, sample_cols = NULL) {
  if (is.null(df_long) || !nrow(df_long)) return(tibble())
  wide <- df_long |>
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0) |>
    arrange(SH_code) |>
    as.data.frame()
  if (!is.null(sample_cols)) {
    missing <- setdiff(sample_cols, colnames(wide))
    for (m in missing) wide[[m]] <- 0
  }
  wide
}
###### SIMPLE SH TABLES ######
# Build SH count tables per sample for D1–D4 used in the manuscript.
# Outputs: derived/mapping_tables/*_{threshold_tag}.csv
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})
set_here("/scratch/gent/vo/001/gvo00142/vsc45818/UMIvsRAW/data2")
here("data")
# Shared helpers a

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# ----- Threshold handling ------------------------------------------------------
thr_tag <- function(x) sprintf("%03d", round(as.numeric(x) * 10))

std_barcode <- function(x) {
  x %>% str_replace("^barcode", "") %>% str_replace_all("_", "") %>% str_replace("rawumis", "")
}

# Standardize thresholded column names to "<prefix>_(status|SH_code)_TAG"
rename_threshold_cols <- function(df, prefix) {
  old <- names(df); new <- old
  m <- str_match(old, "^status \\(([^)]+)\\)$")
  hit <- which(!is.na(m[,2]))
  if (length(hit)) new[hit] <- paste0(prefix, "_status_", thr_tag(m[hit,2]))
  m <- str_match(old, "^SH code \\(([^)]+)\\)$")
  hit <- which(!is.na(m[,2]))
  if (length(hit)) new[hit] <- paste0(prefix, "_SH_code_", thr_tag(m[hit,2]))
  new[old == "Similarity percentage"] <- paste0("Similarity_percentage_", prefix)
  names(df) <- new
  df
}

# Read delimited with headers preserved
read_tab_keep_names <- function(path, delim = "\t") {
  read.delim(path, sep = delim, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
}

# Robust read of matches file (empty-safe)
read_matches_safely <- function(path) {
  if (!file.exists(path)) return(NULL)
  if (file.info(path)$size == 0) return(NULL)
  read_tab_keep_names(path, "\t")
}

# From long SH counts to wide OTU
pivot_otu <- function(df_long, sample_cols = NULL) {
  if (is.null(df_long) || !nrow(df_long)) return(tibble())
  wide <- df_long |>
    tidyr::pivot_wider(names_from = sample, values_from = count, values_fill = 0) |>
    arrange(SH_code) |>
    as.data.frame()
  if (!is.null(sample_cols)) {
    missing <- setdiff(sample_cols, colnames(wide))
    for (m in missing) wide[[m]] <- 0
  }
  wide
}
###### SIMPLE SH TABLES ######
# Build SH count tables per sample for D1–D4 used in the manuscript.
# Outputs: derived/mapping_tables/*_{threshold_tag}.csv
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})
set_here("/scratch/gent/vo/001/gvo00142/vsc45818/UMIvsRAW/data2")
setwd("/scratch/gent/vo/001/gvo00142/vsc45818/UMIvsRAW/data2")
here()# CONFIG
threshold     <- 1.5
threshold_tag <- paste0("_0", gsub("\\.", "", as.character(threshold)))
status_col    <- sprintf("status (%s)", threshold)
sh_code_col   <- sprintf("SH code (%s)", threshold)

dir_consensus <- here("UMI_consensus_SHm")
dir_raw_umi   <- here( "UMI_rawinput_SHm")
dir_all_raw   <- here( "all_rawinput_SHm")
out_dir       <- here("derived", "mapping_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_mapping <- function(mapping_file) {
  lines <- readLines(mapping_file, warn = FALSE)
  tibble(raw = lines) |>
    separate(raw, into = c("barcode", "source_folder"), sep = "->", remove = TRUE) |>
    mutate(across(everything(), ~str_trim(.x))) |>
    mutate(barcode = sub("\\.fasta$", "", barcode))
}

summarise_counts <- function(df, status_col, sh_code_col, sample_name, weight_col = NULL) {
  if (is.null(df)) return(NULL)
  if (!all(c(status_col, sh_code_col) %in% colnames(df))) return(NULL)
  f <- df[df[[status_col]] == "present_in", , drop = FALSE]
  if (!nrow(f)) return(NULL)
  s <- if (is.null(weight_col)) {
    f |>
      group_by(SH_code = .data[[sh_code_col]]) |>
      summarise(count = n(), .groups = "drop")
  } else {
    f |>
      group_by(SH_code = .data[[sh_code_col]]) |>
      summarise(count = sum(.data[[weight_col]], na.rm = TRUE), .groups = "drop")
  }
  s$sample <- sample_name
  s
}

write_csv2_rel <- function(x, stem) {
  out <- file.path(out_dir, paste0(stem, threshold_tag, ".csv"))
  readr::write_delim(x, out, delim = ";")
  message("Wrote: ", out)
  invisible(out)
}

# ---------------- D1: UMI-derived consensus ----------------
mapping_file <- file.path(dir_consensus, "mapping.txt")
stopifnot(file.exists(mapping_file))
mapping <- read_mapping(mapping_file)
umi_consensus_tables <- purrr::map2(
  mapping$source_folder, mapping$barcode,
  ~{
    f <- file.path(dir_consensus, .x, "matches", "matches_out_all.csv")
    d <- read_matches_safely(f)
    summarise_counts(d, status_col, sh_code_col, sample_name = .y)
  }
) |> purrr::compact()
UMI_SHmatched <- pivot_otu(bind_rows(umi_consensus_tables), sample_cols = mapping$barcode)
file_UMI_SHmatched <- write_csv2_rel(UMI_SHmatched, "UMI_SHmatched")

# ---------------- D2: UMI binsize-weighted -----------------
umi_binsize_tables <- purrr::map2(
  mapping$source_folder, mapping$barcode,
  ~{
    f <- file.path(dir_consensus, .x, "matches", "matches_out_all.csv")
    d <- read_matches_safely(f)
    if (is.null(d) || !"seq_accno" %in% names(d)) return(NULL)
    d$multiplier <- suppressWarnings(as.numeric(sub(".*ubs=", "", d$seq_accno)))
    summarise_counts(d, status_col, sh_code_col, sample_name = .y, weight_col = "multiplier")
  }
) |> purrr::compact()
UMI_binsize_SHmatched <- pivot_otu(bind_rows(umi_binsize_tables), sample_cols = mapping$barcode)
file_UMI_binsize <- write_csv2_rel(UMI_binsize_SHmatched, "UMI_binsize_SHmatched")

# ---------------- D3: Raw reads used for UMIs --------------
if (dir.exists(dir_raw_umi)) {
  sub_dirs <- list.dirs(dir_raw_umi, full.names = TRUE, recursive = FALSE)
  sub_dirs <- sub_dirs[grepl("^raw_umi_", basename(sub_dirs))]
  sample_names <- sub("_raw_umis$", "", sub("^raw_umi_", "", basename(sub_dirs)))
  process_one <- function(sub_dir) {
    prefix      <- sub("^raw_umi_", "", basename(sub_dir))
    sample_name <- sub("_raw_umis$", "", prefix)
    path        <- file.path(sub_dir, paste0(prefix, "_all_matches_out.tsv"))
    d <- read_matches_safely(path)
    summarise_counts(d, status_col, sh_code_col, sample_name)
  }
  rawumi_long <- purrr::map(sub_dirs, process_one) |> purrr::compact() |> bind_rows()
  UMI_rawinput_SHmatched <- pivot_otu(rawumi_long, sample_cols = sample_names)
  colnames(UMI_rawinput_SHmatched) <- gsub("_", "", colnames(UMI_rawinput_SHmatched))
  file_UMI_rawinput <- write_csv2_rel(UMI_rawinput_SHmatched, "UMI_rawinput_SHmatched")
}

# ---------------- D4: All raw reads  ------------------------
if (dir.exists(dir_all_raw)) {
  sub_dirs <- list.dirs(dir_all_raw, full.names = TRUE, recursive = FALSE)
  sub_dirs <- sub_dirs[grepl("^parallel_", basename(sub_dirs))]
  sample_names <- sub("^parallel_", "", basename(sub_dirs))
  process_one <- function(sub_dir) {
    prefix      <- sub("^parallel_", "", basename(sub_dir))
    path        <- file.path(sub_dir, paste0(prefix, "_all_matches_out.tsv"))
    d <- read_matches_safely(path)
    summarise_counts(d, status_col, sh_code_col, sample_name = prefix)
  }
  allraw_long <- purrr::map(sub_dirs, process_one) |> purrr::compact() |> bind_rows()
  all_rawinput_SHmatched <- pivot_otu(allraw_long, sample_cols = sample_names)
  colnames(all_rawinput_SHmatched) <- gsub("_", "", colnames(all_rawinput_SHmatched))
  file_all_rawinput <- write_csv2_rel(all_rawinput_SHmatched, "all_rawinput_SHmatched")
}

message("Done with D1–D4 table generation.")

############# ANNOTATION SETS ##############
# Build UMI→fastq map, join raw-input and UMI-consensus matches across thresholds,
# classify input–consensus pairs, and attach taxonomy.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(stringr)
})

# CONFIG
threshold_chr <- "1.5"                     # active tag for file naming
TAG <- thr_tag(threshold_chr)              # "015"
ALL_TAGS <- c("030","025","020","015","010","005")

ssumi_dir     <- here("SSUMI")
umi_ids_dir   <- here("UMI_IDs")
rawinput_dir  <- here("UMI_rawinput_SHm")
consensus_dir <- here("UMI_consensus_SHm")
taxonomy_csv  <- here("UNITE_SH_taxonomy_table.csv")

map_out_dir   <- here("derived","mapping_tables")
final_out_dir <- here("derived","outdata")
dir.create(map_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_out_dir, recursive = TRUE, showWarnings = FALSE)

exclude_bc <- c("66P1","55MI")

# -------------- Part A: UMI → fastq IDs -------------------------
barcode_dirs <- list.files(ssumi_dir, pattern = "^barcode.*_ssumi$", full.names = TRUE)
barcode_info <- tibble(
  spath   = barcode_dirs,
  folder  = basename(barcode_dirs),
  digits  = sub("^barcode([0-9]+)(P1|P2|P3|MI)_.*$", "\\1", basename(barcode_dirs)),
  code    = sub("^barcode[0-9]+(P1|P2|P3|MI)_.*$", "\\1", basename(barcode_dirs))
) %>% mutate(barcode = paste0(digits, code)) %>% filter(!barcode %in% exclude_bc)

umi_map_list <- vector("list", nrow(barcode_info))
for (i in seq_len(nrow(barcode_info))) {
  spath <- barcode_info$spath[i]; bc <- barcode_info$barcode[i]
  out_fastq <- file.path(spath, paste0("raw_UMIs_bc", bc, ".fastq"))
  bins_tgz  <- file.path(spath, "umi_binning", "read_binning", "bins.tar.gz")
  umi_ids_txt <- file.path(umi_ids_dir, paste0("barcode", bc, "_UMI_ids"))
  map_csv_path <- file.path(spath, paste0(bc, "_mapping_umi.csv"))
  if (file.exists(map_csv_path)) { umi_map_list[[i]] <- read.csv2(map_csv_path); next }
  if (!file.exists(bins_tgz) || !file.exists(umi_ids_txt)) next
  
  tmpdir <- tempfile("bins_extract_"); dir.create(tmpdir); on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  untar(bins_tgz, exdir = tmpdir)
  all_files <- list.files(tmpdir, recursive = TRUE, full.names = TRUE)
  file_index <- tibble(basename = basename(all_files), path = all_files) %>% group_by(basename) %>% slice(1) %>% ungroup()
  
  umi_files <- readLines(umi_ids_txt, warn = FALSE)
  out_con <- file(out_fastq, "wb"); on.exit(close(out_con), add = TRUE)
  
  out_list <- vector("list", length(umi_files))
  for (j in seq_along(umi_files)) {
    uf <- umi_files[j]
    hit <- file_index$path[file_index$basename == uf]
    if (!length(hit)) next
    fq <- readLines(hit[1]); writeLines(fq, out_con)
    ids <- fq[seq(1, length(fq), by = 4)]; ids <- sub("^@", "", ids)
    out_list[[j]] <- tibble(UMI_ID = sub("\\.fastq$", "", uf), fastq_ID = ids, Barcode = bc)
  }
  close(out_con)
  umi_map <- bind_rows(out_list)
  if (nrow(umi_map)) readr::write_csv2(umi_map, map_csv_path)
  umi_map_list[[i]] <- umi_map
}
combined_mapping_umi <- bind_rows(umi_map_list) %>% distinct()
readr::write_csv2(combined_mapping_umi, file.path(map_out_dir, "combined_mapping_umi.csv"))

# -------------- Part B: Raw-input matches (TD3) -----------------
raw_subdirs <- list.dirs(rawinput_dir, full.names = TRUE, recursive = FALSE)
raw_subdirs <- raw_subdirs[!str_detect(basename(raw_subdirs), "input_rawumi")]
raw_sets <- map(raw_subdirs, function(sd) {
  files <- list.files(sd, pattern = "_all_matches_out\\.tsv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  df <- read_tab_keep_names(files[1])
  df <- rename_threshold_cols(df, "fastq")
  semi <- str_locate(df$seq_accno, ";")[,1]
  left <- ifelse(is.na(semi), df$seq_accno, substr(df$seq_accno, 1, semi - 1))
  us <- gregexpr("_", left, fixed = TRUE)
  last_us <- vapply(us, function(v) if (length(v)) tail(v, 1) else 0L, integer(1))
  bc  <- ifelse(last_us > 0, substr(left, 1, last_us - 1), NA_character_)
  fid <- ifelse(last_us > 0, substr(left, last_us + 1, nchar(left)), NA_character_)
  smp <- sub(".*sample=", "", df$seq_accno)
  df$fastq_Barcode <- std_barcode(bc); df$fastq_ID <- sub("rc$", "", fid); df$sample <- smp; df$seq_accno <- NULL
  keep <- c("fastq_Barcode","fastq_ID","sample", grep("^fastq_(status|SH_code)_\\d{3}$", names(df), value = TRUE),
            "Similarity_percentage_fastq")
  df[, intersect(keep, names(df)), drop = FALSE]
})
raw_input_all <- bind_rows(compact(raw_sets))

combined_mapping_umi_to_fastq <- combined_mapping_umi %>%
  left_join(raw_input_all, by = c("fastq_ID")) %>%
  mutate(UMI_ID_bc = paste0(UMI_ID, "_bc", Barcode)) %>%
  filter(!is.na(.data[[paste0("fastq_status_", TAG)]]))
readr::write_csv2(combined_mapping_umi_to_fastq,
                  file.path(map_out_dir, paste0("combined_mapping_umi_to_fastq_", TAG, ".csv")))

# -------------- Part C: UMI-consensus matches (TD1) ---------------
map_lines <- readLines(file.path(consensus_dir, "mapping.txt"), warn = FALSE)
umap <- tibble(raw = map_lines) |>
  separate(raw, c("barcode", "source"), "->", remove = TRUE) |>
  mutate(across(everything(), ~str_trim(.x)),
         barcode = sub("\\.fasta$", "", barcode)) |>
  filter(!sub("^barcode", "", barcode) %in% exclude_bc)

umish_list <- map(seq_len(nrow(umap)), function(i) {
  src <- umap$source[i]; bc_full <- umap$barcode[i]
  path <- file.path(consensus_dir, src, "matches", "matches_out_all.csv")
  if (!file.exists(path)) return(NULL)
  d <- read_tab_keep_names(path, "\t") %>%
    mutate(UMI_ID = sapply(strsplit(seq_accno, ";"), `[`, 1),
           binsize = sapply(strsplit(seq_accno, ";"), `[`, 2)) %>%
    select(-seq_accno)
  d <- rename_threshold_cols(d, "umishm")
  keep <- c("UMI_ID", "binsize", grep("^umishm_(status|SH_code)_\\d{3}$", names(d), value = TRUE),
            "Similarity_percentage_umishm")
  d <- d[, intersect(keep, names(d)), drop = FALSE]
  d$Barcode <- std_barcode(bc_full); d$UMI_ID_bc <- paste0(d$UMI_ID, "_bc", d$Barcode)
  d
})
UMI_combined <- bind_rows(compact(umish_list)) %>%
  filter(.data[[paste0("umishm_status_", TAG)]] == "present_in")
readr::write_csv2(UMI_combined, file.path(map_out_dir, paste0("UMI_combined_", TAG, ".csv")))

# -------------- Part D: Merge + taxonomy + per-threshold labels ---
merged <- combined_mapping_umi_to_fastq %>%
  inner_join(UMI_combined %>% select(UMI_ID_bc, starts_with("umishm_")), by = "UMI_ID_bc")
taxo <- readr::read_csv(taxonomy_csv, show_col_types = FALSE) %>% select(SH_code, lowest_taxon)

for (th in ALL_TAGS) {
  f_status <- paste0("fastq_status_", th)
  f_sh     <- paste0("fastq_SH_code_", th)
  u_status <- paste0("umishm_status_", th)
  u_sh     <- paste0("umishm_SH_code_", th)
  if (!all(c(f_sh, u_sh, u_status) %in% names(merged))) next
  annot_col <- paste0("annot_", th)
  merged <- merged %>%
    mutate(
      !!annot_col := if_else(
        .data[[u_status]] == "present_in",
        case_when(
          .data[[f_sh]] == .data[[u_sh]] ~ "matched",
          str_starts(.data[[f_sh]], "SH") & str_starts(.data[[u_sh]], "SH") & .data[[f_sh]] != .data[[u_sh]] ~ "mismatched",
          str_starts(.data[[f_sh]], "SH") & !str_starts(.data[[u_sh]], "SH") ~ "conmis",
          !str_starts(.data[[f_sh]], "SH") & str_starts(.data[[u_sh]], "SH") ~ "unmatched",
          !str_starts(.data[[f_sh]], "SH") & !str_starts(.data[[u_sh]], "SH") ~ "nomatch",
          TRUE ~ "partial_mismatch"
        ),
        NA_character_
      )
    ) %>%
    left_join(taxo %>% rename(!!paste0("lowest_taxon_umishm_", th) := lowest_taxon),
              by = setNames("SH_code", u_sh)) %>%
    left_join(taxo %>% rename(!!paste0("lowest_taxon_fastq_", th) := lowest_taxon),
              by = setNames("SH_code", f_sh))
}
core_cols <- c("UMI_ID_bc","fastq_ID","Barcode",
               grep("^Similarity_percentage_(fastq|umishm)$", names(merged), value = TRUE),
               grep("^fastq_(status|SH_code)_", names(merged), value = TRUE),
               grep("^umishm_(status|SH_code)_", names(merged), value = TRUE),
               grep("^lowest_taxon_(fastq|umishm)_", names(merged), value = TRUE),
               grep("^annot_\\d{3}$", names(merged), value = TRUE))
final_out <- merged[, unique(core_cols[core_cols %in% names(merged)]), drop = FALSE]
readr::write_csv2(final_out, file.path(final_out_dir, paste0("consensus_vs_raw_mapping_", TAG, ".csv")))
message("Wrote consensus_vs_raw_mapping_", TAG, ".csv")

################ Summarize annotation outcomes, produce Q-score summaries, compute p-distance for mismatched SH pairs at the used threshold. ##############
#
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(scales)
  library(Biostrings); library(ape); library(data.table)
})
TAG <- "015"
ALL_TAGS <- c("030","025","020","015","010","005")
in_dir   <- here("derived","outdata")
map_dir  <- here("derived","mapping_tables")
fig_dir  <- here("derived","figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

merged_csv  <- file.path(in_dir, paste0("consensus_vs_raw_mapping_", TAG, ".csv"))
qscores_rds <- file.path(map_dir, "qscores_rawumis")  # optional
refs_fasta  <- here("data_udb", "sanger_refs_sh.fasta")

stopifnot(file.exists(merged_csv))
df_merged <- readr::read_csv2(merged_csv, show_col_types = FALSE) %>%
  mutate(Barcode = if ("fastq_Barcode" %in% names(.)) fastq_Barcode else Barcode)

# Build long table
fastq_sh_long <- df_merged %>%
  dplyr::select(UMI_ID_bc, fastq_ID, Barcode,
                dplyr::matches("^fastq_SH_code_\\d{3}$")) %>%
  tidyr::pivot_longer(
    cols = dplyr::matches("^fastq_SH_code_\\d{3}$"),
    names_to   = "tag",
    values_to  = "fastq_SH_code",
    names_pattern = "^fastq_SH_code_(\\d{3})$"
  )

fastq_tax_long <- df_merged %>%
  dplyr::select(UMI_ID_bc, fastq_ID, Barcode,
                dplyr::matches("^lowest_taxon_fastq_\\d{3}$")) %>%
  tidyr::pivot_longer(
    cols = dplyr::matches("^lowest_taxon_fastq_\\d{3}$"),
    names_to   = "tag",
    values_to  = "lowest_taxon_fastq",
    names_pattern = "^lowest_taxon_fastq_(\\d{3})$"
  )

long_fastq <- fastq_sh_long %>%
  dplyr::left_join(fastq_tax_long,
                   by = c("UMI_ID_bc","fastq_ID","Barcode","tag"))


# SH codes
umi_sh_long <- df_merged %>%
  dplyr::select(UMI_ID_bc, dplyr::matches("^umishm_SH_code_\\d{3}$")) %>%
  tidyr::pivot_longer(
    cols = dplyr::matches("^umishm_SH_code_\\d{3}$"),
    names_to   = "tag",
    values_to  = "umishm_SH_code",
    names_pattern = "^umishm_SH_code_(\\d{3})$"
  ) %>%
  dplyr::group_by(UMI_ID_bc, tag) %>%
  dplyr::summarise(
    umishm_SH_code = dplyr::first(umishm_SH_code[!is.na(umishm_SH_code)]),
    .groups = "drop"
  )

# taxonomy
umi_tax_long <- df_merged %>%
  dplyr::select(UMI_ID_bc, dplyr::matches("^lowest_taxon_umishm_\\d{3}$")) %>%
  tidyr::pivot_longer(
    cols = dplyr::matches("^lowest_taxon_umishm_\\d{3}$"),
    names_to   = "tag",
    values_to  = "lowest_taxon_umishm",
    names_pattern = "^lowest_taxon_umishm_(\\d{3})$"
  ) %>%
  dplyr::group_by(UMI_ID_bc, tag) %>%
  dplyr::summarise(
    lowest_taxon_umishm = dplyr::first(lowest_taxon_umishm[!is.na(lowest_taxon_umishm)]),
    .groups = "drop"
  )

long_umi <- umi_sh_long %>%
  dplyr::left_join(umi_tax_long, by = c("UMI_ID_bc","tag"))

long_annot <- df_merged %>%
  dplyr::select(UMI_ID_bc, fastq_ID, Barcode, dplyr::starts_with("annot_")) %>%
  tidyr::pivot_longer(
    dplyr::starts_with("annot_"),
    names_to   = "tag",
    values_to  = "annotation",
    values_drop_na = TRUE
  ) %>%
  dplyr::mutate(tag = sub("^annot_", "", tag))
simcols <- df_merged %>%
  dplyr::select(
    UMI_ID_bc, fastq_ID,
    dplyr::any_of(c("Similarity_percentage_fastq", "Similarity_percentage_umishm"))
  ) %>%
  dplyr::distinct()
totaldf <- long_annot %>%
  dplyr::left_join(long_fastq, by = c("UMI_ID_bc","fastq_ID","Barcode","tag")) %>%
  dplyr::left_join(long_umi,   by = c("UMI_ID_bc","tag")) %>%
  dplyr::left_join(simcols,    by = c("UMI_ID_bc","fastq_ID")) %>%
  dplyr::mutate(threshold = paste0(as.numeric(tag)/10, "%")) %>%
  dplyr::select(
    Barcode, fastq_ID,
    fastq_SH_code, fastq_lowest_taxon = lowest_taxon_fastq,
    UMI_ID_bc, umi_SH_code = umishm_SH_code, umi_lowest_taxon = lowest_taxon_umishm,
    threshold, annotation,
    dplyr::any_of(c("Similarity_percentage_fastq", "Similarity_percentage_umishm"))
  )

if (file.exists(qscores_rds)) {
  qscores <- readRDS(qscores_rds)
  totaldf <- totaldf %>% left_join(qscores, by = c("fastq_ID" = "fastq_id"))
}

readr::write_csv2(totaldf, file.path(in_dir, "all_matchdata.csv"))

# Mismatch summaries
mismatches_df <- totaldf %>% filter(annotation == "mismatched")
readr::write_csv2(mismatches_df, file.path(in_dir, "mismatches_all.csv"))
readr::write_csv2(mismatches_df %>% filter(threshold == "1.5%"), file.path(in_dir, "mismatches_015.csv"))

# Per-barcode summaries and global proportions
match_counts <- totaldf %>%
  group_by(Barcode, threshold, annotation) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(Barcode, threshold) %>%
  mutate(total_nonNA = sum(n[!is.na(annotation)]),
         percentage  = if_else(is.na(annotation), NA_real_, 100 * n / total_nonNA)) %>%
  ungroup()
readr::write_csv2(match_counts, file.path(in_dir, "match_counts_threshold_barcode.csv"))

overall_match_counts <- totaldf %>% group_by(threshold, annotation) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(threshold) %>% mutate(total = sum(n), percentage = round(100 * n / total, 2)) %>% ungroup()
readr::write_csv2(overall_match_counts, file.path(in_dir, "match_counts_threshold_alltypes.csv"))

# Q-score plots saved only if available
if ("avg_qscore" %in% names(totaldf)) {
  annotation_colors <- c(matched="#1B9E77", unmatched="#D95F02", mismatched="#4575B4")
  p <- totaldf %>%
    ggplot(aes(x = annotation, y = avg_qscore, fill = annotation)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
    facet_wrap(~ threshold, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = annotation_colors) +
    labs(x = NULL, y = "Average Q-score") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.placement = "outside")
  print(p)
  ggsave(file.path(fig_dir, "qscore_violin_by_annotation.png"), p, width = 10, height = 3.2, dpi = 300)
}

# Reference sequences and p-distance for mismatches at 1.5%
if (file.exists(refs_fasta)) {
  dna <- readDNAStringSet(refs_fasta)
  hdr <- names(dna)
  sh_code <- vapply(strsplit(sub("^>", "", hdr), "_", fixed = TRUE),
                    function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  ref_dt <- data.table(sh_code = sh_code, sequence = as.character(dna), seq_id = hdr)[!is.na(sh_code)]
  setkey(ref_dt, sh_code); saveRDS(ref_dt, file.path(map_dir, "ref_sequences.rds"))
  
  df_015 <- totaldf %>% filter(threshold == "1.5%", annotation == "mismatched")
  if (nrow(df_015)) {
    df_015_dt <- as.data.table(df_015); setkey(ref_dt, sh_code)
    df_015_dt <- merge(df_015_dt, ref_dt[, .(sh_code, fastq_ref_seq = sequence)],
                       by.x = "fastq_SH_code", by.y = "sh_code", all.x = TRUE)
    df_015_dt <- merge(df_015_dt, ref_dt[, .(sh_code, umi_ref_seq = sequence)],
                       by.x = "umi_SH_code", by.y = "sh_code", all.x = TRUE)
    p_distance_ape <- function(seq1, seq2) {
      if (is.na(seq1) || is.na(seq2)) return(NA_real_)
      aln <- pairwiseAlignment(DNAString(seq1), DNAString(seq2),
                               substitutionMatrix = NULL, gapOpening = -10, gapExtension = -1, type = "global")
      aln_set <- DNAStringSet(c(as.character(pattern(aln)), as.character(subject(aln))))
      d <- dist.dna(as.DNAbin(aln_set), model = "raw", pairwise.deletion = TRUE)
      as.numeric(d)
    }
    pvals <- vapply(seq_len(nrow(df_015_dt)), function(i)
      p_distance_ape(df_015_dt$fastq_ref_seq[i], df_015_dt$umi_ref_seq[i]), numeric(1))
    collapsed_details <- df_015_dt %>%
      mutate(pairwise_distance = pvals) %>%
      group_by(fastq_SH_code, umi_SH_code, fastq_lowest_taxon, umi_lowest_taxon) %>%
      summarise(count = n(), avg_pairwise_distance = mean(pairwise_distance, na.rm = TRUE), .groups = "drop")
    readr::write_csv2(collapsed_details, file.path(in_dir, "mismatch015_collapsed_SHs.csv"))
  }
}
message("Annotation summaries and p-distance done.")
######## Create comparable sample set at 1.5% for CD1–CD4 and run community analyses ##########
suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(here)
})

lists_dir <- here("derived","mapping_tables")
stats_dir <- here("derived","stats")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  UMI_corrected = file.path(lists_dir, "UMI_corrected_abundance_SHmatched_list.rds"),
  UMI_binsize   = file.path(lists_dir, "UMI_binsize_SHmatched_list.rds"),
  UMI_rawinput  = file.path(lists_dir, "UMI_rawinput_SHmatched_list.rds"),
  all_rawinput  = file.path(lists_dir, "all_rawinput_SHmatched_list.rds")
)
stopifnot(all(file.exists(paths)))
L <- lapply(paths, readRDS)
thr <- "1.5"
tab_15 <- lapply(L, function(x) x[[thr]])

sample_cols <- lapply(tab_15, function(df) grep("^barcode", names(df), value = TRUE))
common_samp <- Reduce(intersect, sample_cols)
filt_common <- lapply(tab_15, function(df) df[, c("SH_code", common_samp), drop = FALSE])

to_vegan <- function(df) {
  mat <- df %>% select(-SH_code) %>% as.matrix()
  mat <- t(mat); colnames(mat) <- df$SH_code; mat
}
v_mats <- lapply(filt_common, to_vegan)
saveRDS(v_mats, file.path(stats_dir, "beta_common_samples_vegan.rds"))

# Distances
dist_bc   <- lapply(v_mats, function(m) vegdist(decostand(m, "total"), method = "bray"))
dist_hel  <- lapply(v_mats, function(m) vegdist(decostand(m, "hellinger"), method = "euclidean"))
dist_ait  <- lapply(v_mats, function(m) tryCatch(vegdist(m, method = "robust.aitchison"),
                                                 error = function(e) vegdist(decostand(m, "total"), method = "bray")))

# NMDS and Procrustes
set.seed(42)
nmds_bc  <- lapply(dist_bc,  function(d) metaMDS(d, k = 2, trymax = 100))
nmds_hel <- lapply(dist_hel, function(d) metaMDS(d, k = 2, trymax = 100))
nmds_ait <- lapply(dist_ait, function(d) metaMDS(d, k = 2, trymax = 100))

pairs <- combn(names(v_mats), 2, simplify = FALSE)
pro_tab <- function(nmds_list, label) bind_rows(lapply(pairs, function(p) {
  pr <- protest(nmds_list[[p[1]]], nmds_list[[p[2]]], permutations = 999)
  tibble(Distance = label, Comparison = paste(p, collapse = " vs "),
         Protest_correlation = pr$t0, Protest_p = pr$signif)
})) %>% mutate(Protest_adj_p = p.adjust(Protest_p, method = "BH"))
pro_bc  <- pro_tab(nmds_bc,  "Bray-Curtis")
pro_hel <- pro_tab(nmds_hel, "Hellinger")
pro_ait <- pro_tab(nmds_ait, "Robust-Aitchison")
readr::write_csv2(bind_rows(pro_bc, pro_hel, pro_ait), file.path(stats_dir, "beta_procrustes_all.csv"))

# PERMANOVA, MRPP, betadisper
comb_df <- bind_rows(lapply(v_mats, as.data.frame))
comb_df[is.na(comb_df)] <- 0
group_df <- tibble(method = factor(rep(names(v_mats), times = sapply(v_mats, nrow))))

ad_bc  <- adonis2(comb_df ~ method, data = group_df, method = "bray")
ad_hel <- adonis2(decostand(comb_df, "hellinger") ~ method, data = group_df, method = "bray")
ad_ait <- tryCatch(adonis2(comb_df ~ method, data = group_df, method = "robust.aitchison"),
                   error = function(e) adonis2(comb_df ~ method, data = group_df, method = "bray"))
broom_permanova <- function(x) as_tibble(x) %>% mutate(term = rownames(x), .before = 1)
readr::write_csv2(broom_permanova(ad_bc),  file.path(stats_dir, "beta_permanova_bray.csv"))
readr::write_csv2(broom_permanova(ad_hel), file.path(stats_dir, "beta_permanova_hellinger.csv"))
readr::write_csv2(broom_permanova(ad_ait), file.path(stats_dir, "beta_permanova_aitchison.csv"))

mr_bc  <- mrpp(comb_df, group_df$method, distance = "bray")
mr_hel <- mrpp(decostand(comb_df, "hellinger"), group_df$method, distance = "bray")
mr_ait <- tryCatch(mrpp(comb_df, group_df$method, distance = "robust.aitchison"),
                   error = function(e) mrpp(comb_df, group_df$method, distance = "bray"))
mrpp_tbl <- function(obj, label) tibble(test = label, delta = obj$delta, p_value = obj$Pvalue)
readr::write_csv2(bind_rows(mrpp_tbl(mr_bc,"bray"), mrpp_tbl(mr_hel,"hellinger"), mrpp_tbl(mr_ait,"aitchison")),
                  file.path(stats_dir, "beta_mrpp.csv"))

bd_bc  <- betadisper(vegdist(comb_df, "bray"), group_df$method)
bd_hel <- betadisper(vegdist(decostand(comb_df, "hellinger"), "bray"), group_df$method)
bd_ait <- betadisper(vegdist(comb_df, "robust.aitchison"), group_df$method)
bd_out <- bind_rows(
  tibble(test="bray",      F=anova(bd_bc)[1,"F value"],  p_value=anova(bd_bc)[1,"Pr(>F)"]),
  tibble(test="hellinger", F=anova(bd_hel)[1,"F value"], p_value=anova(bd_hel)[1,"Pr(>F)"]),
  tibble(test="aitchison", F=tryCatch(anova(bd_ait)[1,"F value"], error=function(e) NA_real_),
         p_value=tryCatch(anova(bd_ait)[1,"Pr(>F)"], error=function(e) NA_real_))
)
readr::write_csv2(bd_out, file.path(stats_dir, "beta_betadisper_anova.csv"))

# Mantel
mantel_tab <- function(d1, d2, label1, label2) {
  mt <- mantel(d1, d2, method = "spearman", permutations = 999)
  tibble(Matrix1 = label1, Matrix2 = label2, Mantel_r = mt$statistic, p = mt$signif)
}
mantel_out <- bind_rows(
  mantel_tab(dist_bc$UMI_binsize, dist_bc$UMI_corrected, "BC CD2", "BC CD1"),
  mantel_tab(dist_bc$UMI_binsize, dist_bc$UMI_rawinput,  "BC CD2", "BC CD3"),
  mantel_tab(dist_bc$UMI_corrected, dist_bc$UMI_rawinput,"BC CD1", "BC CD3")
)
readr::write_csv2(mantel_out, file.path(stats_dir, "beta_mantel.csv"))

# dbRDA variance partitioning (dataset vs barcode as strata proxy)
# This estimates the proportion of variation explained by method, with barcode differences dominating.
make_dbRDA <- function(m, method_label) {
  df <- as.data.frame(m); df$barcode <- rownames(df); df$method <- sub(".*\\.", "", df$barcode)
  cap <- capscale(df[, setdiff(names(df), c("barcode","method"))] ~ method, distance = "bray")
  rsq <- RsquareAdj(cap)$r.squared
  tibble(Distance = method_label, Dataset_variance_pct = 100 * rsq)
}
# The manuscript reports near-zero dataset variance and near-100% barcode variance.
# Here we output only the dataset component.
db <- bind_rows(
  make_dbRDA(v_mats$UMI_corrected, "Bray-Curtis"),
  make_dbRDA(v_mats$UMI_binsize,   "Bray-Curtis"),
  make_dbRDA(v_mats$UMI_rawinput,  "Bray-Curtis")
)
readr::write_csv2(db, file.path(stats_dir, "beta_dbRDA_dataset_variance.csv"))

message("Community structure tests written to: ", stats_dir)
# Exclusive CD4 taxa audit and UMI binsize histograms for mock and all data.
suppressPackageStartupMessages({
  library(tidyverse); library(here)
})
THR <- "1.5"
out_dir <- here("derived","outdata"); fig_dir <- here("derived","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

lists_dir <- here("derived","mapping_tables")
paths <- c(
  UMI_binsize   = file.path(lists_dir, "UMI_binsize_SHmatched_list.rds"),
  UMI_rawinput  = file.path(lists_dir, "UMI_rawinput_SHmatched_list.rds"),
  all_rawinput  = file.path(lists_dir, "all_rawinput_SHmatched_list.rds"),
  UMI_corrected = file.path(lists_dir, "UMI_corrected_abundance_SHmatched_list.rds")
)
stopifnot(all(file.exists(paths)))
L <- lapply(paths, readRDS); tab_15 <- lapply(L, function(x) x[[THR]])

sample_cols <- lapply(tab_15, function(df) grep("^barcode", names(df), value = TRUE))
common_cols <- Reduce(intersect, sample_cols)
df_all <- tab_15$all_rawinput[,  c("SH_code", common_cols), drop = FALSE]
df_con <- tab_15$UMI_corrected[, c("SH_code", common_cols), drop = FALSE]
df_umi <- tab_15$UMI_rawinput[,  c("SH_code", common_cols), drop = FALSE]

samp_all <- setdiff(names(df_all), "SH_code")
samp_con <- setdiff(names(df_con), "SH_code")
samp_umi <- setdiff(names(df_umi), "SH_code")

all_presence <- df_all %>%
  mutate(n_present = rowSums(across(all_of(samp_all), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE),
         total_abundance = rowSums(across(all_of(samp_all), ~ suppressWarnings(as.numeric(.))), na.rm = TRUE)) %>%
  filter(n_present > 0) %>% select(SH_code, n_present, total_abundance)

present_ref <- union(
  df_con %>% mutate(n_present = rowSums(across(all_of(samp_con), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE)) %>%
    filter(n_present > 0) %>% pull(SH_code) %>% unique(),
  df_umi %>% mutate(n_present = rowSums(across(all_of(samp_umi), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE)) %>%
    filter(n_present > 0) %>% pull(SH_code) %>% unique()
)
res <- all_presence %>% filter(!(SH_code %in% present_ref)) %>% arrange(desc(total_abundance), desc(n_present))
readr::write_csv2(res, file.path(out_dir, "exclusive_taxa_all_data_not_in_consensus_or_UMI.csv"))

total_reads_per_sample <- df_all %>% summarise(across(all_of(samp_all), ~ sum(suppressWarnings(as.numeric(.)), na.rm = TRUE)))
exclusive_reads_per_sample <- df_all %>% filter(SH_code %in% res$SH_code) %>%
  summarise(across(all_of(samp_all), ~ sum(suppressWarnings(as.numeric(.)), na.rm = TRUE)))
total_vec <- as.numeric(total_reads_per_sample[1, ]); names(total_vec) <- names(total_reads_per_sample)
excl_vec  <- if (ncol(exclusive_reads_per_sample) > 0) as.numeric(exclusive_reads_per_sample[1, ]) else rep(0, length(total_vec))
pct_vec   <- 100 * excl_vec / ifelse(total_vec > 0, total_vec, NA_real_)
per_sample_df <- tibble(sample = names(total_vec), exclusive_reads = excl_vec, total_reads = total_vec, pct_exclusive = pct_vec) %>%
  arrange(desc(pct_exclusive))
readr::write_csv2(per_sample_df, file.path(out_dir, "exclusive_taxa_share_per_sample.csv"))

# UMI binsize histograms from FASTA headers
ssumi_dir <- Sys.getenv("UMIVSRAW_SSUMI_DIR", unset = here("raw","SSUMI"))
pattern_all  <- file.path(ssumi_dir, "barcode*_ssumi", "consensus_raconx3_medakax2_raconx1.fa")
pattern_mock <- file.path(ssumi_dir, "barcode56MI_ssumi", "consensus_raconx3_medakax2_raconx1.fa")

extract_ubs <- function(fp) {
  hdr <- readr::read_lines(fp, progress = FALSE); hdr <- hdr[startsWith(hdr, ">")]
  if (!length(hdr)) return(tibble(file = character(0), ubs = integer(0)))
  tibble(file = fp, ubs = suppressWarnings(as.integer(stringr::str_match(hdr, "(?i)\\bubs=([0-9]+)\\b")[,2])))
}

plot_ubs <- function(files, out_png, x_break_every = 50) {
  if (!length(files)) return(invisible(NULL))
  ubs_tbl <- purrr::map_dfr(files, extract_ubs) %>% filter(!is.na(ubs))
  if (!nrow(ubs_tbl)) return(invisible(NULL))
  p <- ggplot(ubs_tbl, aes(x = factor(ubs))) + geom_bar() + scale_y_log10() +
    scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = x_break_every)]) +
    labs(title = "Prevalence per UMI bin size", x = "ubs", y = "Sequence count (log10)") +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(filename = file.path(fig_dir, out_png), plot = p, dpi = 900, height = 5.5, width = 7.8)
  invisible(ubs_tbl)
}
plot_ubs(Sys.glob(pattern_all),  "umi_binsize_all.png",  x_break_every = 50)
plot_ubs(Sys.glob(pattern_mock), "umi_binsize_mock.png", x_break_every = 3)
message("Exclusive CD4 tables and binsize figures written.")
)
# CONFIG
threshold     <- 1.5
threshold_tag <- paste0("_0", gsub("\\.", "", as.character(threshold)))
status_col    <- sprintf("status (%s)", threshold)
sh_code_col   <- sprintf("SH code (%s)", threshold)

dir_consensus <- here("data", "UMI_consensus_SHm")
dir_raw_umi   <- here("data", "UMI_rawinput_SHm")
dir_all_raw   <- here("data", "all_rawinput_SHm")
out_dir       <- here("derived", "mapping_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_mapping <- function(mapping_file) {
  lines <- readLines(mapping_file, warn = FALSE)
  tibble(raw = lines) |>
    separate(raw, into = c("barcode", "source_folder"), sep = "->", remove = TRUE) |>
    mutate(across(everything(), ~str_trim(.x))) |>
    mutate(barcode = sub("\\.fasta$", "", barcode))
}

summarise_counts <- function(df, status_col, sh_code_col, sample_name, weight_col = NULL) {
  if (is.null(df)) return(NULL)
  if (!all(c(status_col, sh_code_col) %in% colnames(df))) return(NULL)
  f <- df[df[[status_col]] == "present_in", , drop = FALSE]
  if (!nrow(f)) return(NULL)
  s <- if (is.null(weight_col)) {
    f |>
      group_by(SH_code = .data[[sh_code_col]]) |>
      summarise(count = n(), .groups = "drop")
  } else {
    f |>
      group_by(SH_code = .data[[sh_code_col]]) |>
      summarise(count = sum(.data[[weight_col]], na.rm = TRUE), .groups = "drop")
  }
  s$sample <- sample_name
  s
}

write_csv2_rel <- function(x, stem) {
  out <- file.path(out_dir, paste0(stem, threshold_tag, ".csv"))
  readr::write_delim(x, out, delim = ";")
  message("Wrote: ", out)
  invisible(out)
}

# ---------------- D1: UMI-derived consensus ----------------
mapping_file <- file.path(dir_consensus, "mapping.txt")
stopifnot(file.exists(mapping_file))
mapping <- read_mapping(mapping_file)
umi_consensus_tables <- purrr::map2(
  mapping$source_folder, mapping$barcode,
  ~{
    f <- file.path(dir_consensus, .x, "matches", "matches_out_all.csv")
    d <- read_matches_safely(f)
    summarise_counts(d, status_col, sh_code_col, sample_name = .y)
  }
) |> purrr::compact()
UMI_SHmatched <- pivot_otu(bind_rows(umi_consensus_tables), sample_cols = mapping$barcode)
file_UMI_SHmatched <- write_csv2_rel(UMI_SHmatched, "UMI_SHmatched")

# ---------------- D2: UMI binsize-weighted -----------------
umi_binsize_tables <- purrr::map2(
  mapping$source_folder, mapping$barcode,
  ~{
    f <- file.path(dir_consensus, .x, "matches", "matches_out_all.csv")
    d <- read_matches_safely(f)
    if (is.null(d) || !"seq_accno" %in% names(d)) return(NULL)
    d$multiplier <- suppressWarnings(as.numeric(sub(".*ubs=", "", d$seq_accno)))
    summarise_counts(d, status_col, sh_code_col, sample_name = .y, weight_col = "multiplier")
  }
) |> purrr::compact()
UMI_binsize_SHmatched <- pivot_otu(bind_rows(umi_binsize_tables), sample_cols = mapping$barcode)
file_UMI_binsize <- write_csv2_rel(UMI_binsize_SHmatched, "UMI_binsize_SHmatched")

# ---------------- D3: Raw reads used for UMIs --------------
if (dir.exists(dir_raw_umi)) {
  sub_dirs <- list.dirs(dir_raw_umi, full.names = TRUE, recursive = FALSE)
  sub_dirs <- sub_dirs[grepl("^raw_umi_", basename(sub_dirs))]
  sample_names <- sub("_raw_umis$", "", sub("^raw_umi_", "", basename(sub_dirs)))
  process_one <- function(sub_dir) {
    prefix      <- sub("^raw_umi_", "", basename(sub_dir))
    sample_name <- sub("_raw_umis$", "", prefix)
    path        <- file.path(sub_dir, paste0(prefix, "_all_matches_out.tsv"))
    d <- read_matches_safely(path)
    summarise_counts(d, status_col, sh_code_col, sample_name)
  }
  rawumi_long <- purrr::map(sub_dirs, process_one) |> purrr::compact() |> bind_rows()
  UMI_rawinput_SHmatched <- pivot_otu(rawumi_long, sample_cols = sample_names)
  colnames(UMI_rawinput_SHmatched) <- gsub("_", "", colnames(UMI_rawinput_SHmatched))
  file_UMI_rawinput <- write_csv2_rel(UMI_rawinput_SHmatched, "UMI_rawinput_SHmatched")
}

# ---------------- D4: All raw reads  ------------------------
if (dir.exists(dir_all_raw)) {
  sub_dirs <- list.dirs(dir_all_raw, full.names = TRUE, recursive = FALSE)
  sub_dirs <- sub_dirs[grepl("^parallel_", basename(sub_dirs))]
  sample_names <- sub("^parallel_", "", basename(sub_dirs))
  process_one <- function(sub_dir) {
    prefix      <- sub("^parallel_", "", basename(sub_dir))
    path        <- file.path(sub_dir, paste0(prefix, "_all_matches_out.tsv"))
    d <- read_matches_safely(path)
    summarise_counts(d, status_col, sh_code_col, sample_name = prefix)
  }
  allraw_long <- purrr::map(sub_dirs, process_one) |> purrr::compact() |> bind_rows()
  all_rawinput_SHmatched <- pivot_otu(allraw_long, sample_cols = sample_names)
  colnames(all_rawinput_SHmatched) <- gsub("_", "", colnames(all_rawinput_SHmatched))
  file_all_rawinput <- write_csv2_rel(all_rawinput_SHmatched, "all_rawinput_SHmatched")
}

message("Done with D1–D4 table generation.")

############# ANNOTATION SETS ##############
# Build UMI→fastq map, join raw-input and UMI-consensus matches across thresholds,
# classify input–consensus pairs, and attach taxonomy.
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(stringr)
})
source(here("00_utils.R"))

# CONFIG
threshold_chr <- "1.5"                     # active tag for file naming
TAG <- thr_tag(threshold_chr)              # "015"
ALL_TAGS <- c("030","025","020","015","010","005")

ssumi_dir     <- here("data","SSUMI")
umi_ids_dir   <- here("data","UMI_IDs")
rawinput_dir  <- here("data","UMI_rawinput_SHm")
consensus_dir <- here("data","UMI_consensus_SHm")
taxonomy_csv  <- here("data","UNITE_SH_taxonomy_table.csv")

map_out_dir   <- here("derived","mapping_tables")
final_out_dir <- here("derived","outdata")
dir.create(map_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_out_dir, recursive = TRUE, showWarnings = FALSE)

exclude_bc <- c("66P1","55MI")

# -------------- Part A: UMI → fastq IDs -------------------------
barcode_dirs <- list.files(ssumi_dir, pattern = "^barcode.*_ssumi$", full.names = TRUE)
barcode_info <- tibble(
  spath   = barcode_dirs,
  folder  = basename(barcode_dirs),
  digits  = sub("^barcode([0-9]+)(P1|P2|P3|MI)_.*$", "\\1", basename(barcode_dirs)),
  code    = sub("^barcode[0-9]+(P1|P2|P3|MI)_.*$", "\\1", basename(barcode_dirs))
) %>% mutate(barcode = paste0(digits, code)) %>% filter(!barcode %in% exclude_bc)

umi_map_list <- vector("list", nrow(barcode_info))
for (i in seq_len(nrow(barcode_info))) {
  spath <- barcode_info$spath[i]; bc <- barcode_info$barcode[i]
  out_fastq <- file.path(spath, paste0("raw_UMIs_bc", bc, ".fastq"))
  bins_tgz  <- file.path(spath, "umi_binning", "read_binning", "bins.tar.gz")
  umi_ids_txt <- file.path(umi_ids_dir, paste0("barcode", bc, "_UMI_ids"))
  map_csv_path <- file.path(spath, paste0(bc, "_mapping_umi.csv"))
  if (file.exists(map_csv_path)) { umi_map_list[[i]] <- read.csv2(map_csv_path); next }
  if (!file.exists(bins_tgz) || !file.exists(umi_ids_txt)) next
  
  tmpdir <- tempfile("bins_extract_"); dir.create(tmpdir); on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  untar(bins_tgz, exdir = tmpdir)
  all_files <- list.files(tmpdir, recursive = TRUE, full.names = TRUE)
  file_index <- tibble(basename = basename(all_files), path = all_files) %>% group_by(basename) %>% slice(1) %>% ungroup()
  
  umi_files <- readLines(umi_ids_txt, warn = FALSE)
  out_con <- file(out_fastq, "wb"); on.exit(close(out_con), add = TRUE)
  
  out_list <- vector("list", length(umi_files))
  for (j in seq_along(umi_files)) {
    uf <- umi_files[j]
    hit <- file_index$path[file_index$basename == uf]
    if (!length(hit)) next
    fq <- readLines(hit[1]); writeLines(fq, out_con)
    ids <- fq[seq(1, length(fq), by = 4)]; ids <- sub("^@", "", ids)
    out_list[[j]] <- tibble(UMI_ID = sub("\\.fastq$", "", uf), fastq_ID = ids, Barcode = bc)
  }
  close(out_con)
  umi_map <- bind_rows(out_list)
  if (nrow(umi_map)) readr::write_csv2(umi_map, map_csv_path)
  umi_map_list[[i]] <- umi_map
}
combined_mapping_umi <- bind_rows(umi_map_list) %>% distinct()
readr::write_csv2(combined_mapping_umi, file.path(map_out_dir, "combined_mapping_umi.csv"))

# -------------- Part B: Raw-input matches (TD3) -----------------
raw_subdirs <- list.dirs(rawinput_dir, full.names = TRUE, recursive = FALSE)
raw_subdirs <- raw_subdirs[!str_detect(basename(raw_subdirs), "input_rawumi")]
raw_sets <- map(raw_subdirs, function(sd) {
  files <- list.files(sd, pattern = "_all_matches_out\\.tsv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  df <- read_tab_keep_names(files[1])
  df <- rename_threshold_cols(df, "fastq")
  semi <- str_locate(df$seq_accno, ";")[,1]
  left <- ifelse(is.na(semi), df$seq_accno, substr(df$seq_accno, 1, semi - 1))
  us <- gregexpr("_", left, fixed = TRUE)
  last_us <- vapply(us, function(v) if (length(v)) tail(v, 1) else 0L, integer(1))
  bc  <- ifelse(last_us > 0, substr(left, 1, last_us - 1), NA_character_)
  fid <- ifelse(last_us > 0, substr(left, last_us + 1, nchar(left)), NA_character_)
  smp <- sub(".*sample=", "", df$seq_accno)
  df$fastq_Barcode <- std_barcode(bc); df$fastq_ID <- sub("rc$", "", fid); df$sample <- smp; df$seq_accno <- NULL
  keep <- c("fastq_Barcode","fastq_ID","sample", grep("^fastq_(status|SH_code)_\\d{3}$", names(df), value = TRUE),
            "Similarity_percentage_fastq")
  df[, intersect(keep, names(df)), drop = FALSE]
})
raw_input_all <- bind_rows(compact(raw_sets))

combined_mapping_umi_to_fastq <- combined_mapping_umi %>%
  left_join(raw_input_all, by = c("fastq_ID")) %>%
  mutate(UMI_ID_bc = paste0(UMI_ID, "_bc", Barcode)) %>%
  filter(!is.na(.data[[paste0("fastq_status_", TAG)]]))
readr::write_csv2(combined_mapping_umi_to_fastq,
                  file.path(map_out_dir, paste0("combined_mapping_umi_to_fastq_", TAG, ".csv")))

# -------------- Part C: UMI-consensus matches (TD1) ---------------
map_lines <- readLines(file.path(consensus_dir, "mapping.txt"), warn = FALSE)
umap <- tibble(raw = map_lines) |>
  separate(raw, c("barcode", "source"), "->", remove = TRUE) |>
  mutate(across(everything(), ~str_trim(.x)),
         barcode = sub("\\.fasta$", "", barcode)) |>
  filter(!sub("^barcode", "", barcode) %in% exclude_bc)

umish_list <- map(seq_len(nrow(umap)), function(i) {
  src <- umap$source[i]; bc_full <- umap$barcode[i]
  path <- file.path(consensus_dir, src, "matches", "matches_out_all.csv")
  if (!file.exists(path)) return(NULL)
  d <- read_tab_keep_names(path, "\t") %>%
    mutate(UMI_ID = sapply(strsplit(seq_accno, ";"), `[`, 1),
           binsize = sapply(strsplit(seq_accno, ";"), `[`, 2)) %>%
    select(-seq_accno)
  d <- rename_threshold_cols(d, "umishm")
  keep <- c("UMI_ID", "binsize", grep("^umishm_(status|SH_code)_\\d{3}$", names(d), value = TRUE),
            "Similarity_percentage_umishm")
  d <- d[, intersect(keep, names(d)), drop = FALSE]
  d$Barcode <- std_barcode(bc_full); d$UMI_ID_bc <- paste0(d$UMI_ID, "_bc", d$Barcode)
  d
})
UMI_combined <- bind_rows(compact(umish_list)) %>%
  filter(.data[[paste0("umishm_status_", TAG)]] == "present_in")
readr::write_csv2(UMI_combined, file.path(map_out_dir, paste0("UMI_combined_", TAG, ".csv")))

# -------------- Part D: Merge + taxonomy + per-threshold labels ---
merged <- combined_mapping_umi_to_fastq %>%
  inner_join(UMI_combined %>% select(UMI_ID_bc, starts_with("umishm_")), by = "UMI_ID_bc")
taxo <- readr::read_csv(taxonomy_csv, show_col_types = FALSE) %>% select(SH_code, lowest_taxon)

for (th in ALL_TAGS) {
  f_status <- paste0("fastq_status_", th)
  f_sh     <- paste0("fastq_SH_code_", th)
  u_status <- paste0("umishm_status_", th)
  u_sh     <- paste0("umishm_SH_code_", th)
  if (!all(c(f_sh, u_sh, u_status) %in% names(merged))) next
  annot_col <- paste0("annot_", th)
  merged <- merged %>%
    mutate(
      !!annot_col := if_else(
        .data[[u_status]] == "present_in",
        case_when(
          .data[[f_sh]] == .data[[u_sh]] ~ "matched",
          str_starts(.data[[f_sh]], "SH") & str_starts(.data[[u_sh]], "SH") & .data[[f_sh]] != .data[[u_sh]] ~ "mismatched",
          str_starts(.data[[f_sh]], "SH") & !str_starts(.data[[u_sh]], "SH") ~ "conmis",
          !str_starts(.data[[f_sh]], "SH") & str_starts(.data[[u_sh]], "SH") ~ "unmatched",
          !str_starts(.data[[f_sh]], "SH") & !str_starts(.data[[u_sh]], "SH") ~ "nomatch",
          TRUE ~ "partial_mismatch"
        ),
        NA_character_
      )
    ) %>%
    left_join(taxo %>% rename(!!paste0("lowest_taxon_umishm_", th) := lowest_taxon),
              by = setNames("SH_code", u_sh)) %>%
    left_join(taxo %>% rename(!!paste0("lowest_taxon_fastq_", th) := lowest_taxon),
              by = setNames("SH_code", f_sh))
}
core_cols <- c("UMI_ID_bc","fastq_ID","Barcode",
               grep("^Similarity_percentage_(fastq|umishm)$", names(merged), value = TRUE),
               grep("^fastq_(status|SH_code)_", names(merged), value = TRUE),
               grep("^umishm_(status|SH_code)_", names(merged), value = TRUE),
               grep("^lowest_taxon_(fastq|umishm)_", names(merged), value = TRUE),
               grep("^annot_\\d{3}$", names(merged), value = TRUE))
final_out <- merged[, unique(core_cols[core_cols %in% names(merged)]), drop = FALSE]
readr::write_csv2(final_out, file.path(final_out_dir, paste0("consensus_vs_raw_mapping_", TAG, ".csv")))
message("Wrote consensus_vs_raw_mapping_", TAG, ".csv")

################ Summarize annotation outcomes, produce Q-score summaries, compute p-distance for mismatched SH pairs at the used threshold. ##############
#
suppressPackageStartupMessages({
  library(tidyverse); library(here); library(scales)
  library(Biostrings); library(ape); library(data.table)
})
TAG <- "015"
ALL_TAGS <- c("030","025","020","015","010","005")
in_dir   <- here("derived","outdata")
map_dir  <- here("derived","mapping_tables")
fig_dir  <- here("derived","figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

merged_csv  <- file.path(in_dir, paste0("consensus_vs_raw_mapping_", TAG, ".csv"))
qscores_rds <- file.path(map_dir, "qscores_rawumis")  # optional
refs_fasta  <- here("data_udb", "sanger_refs_sh.fasta")

stopifnot(file.exists(merged_csv))
df_merged <- readr::read_csv2(merged_csv, show_col_types = FALSE) %>%
  mutate(Barcode = if ("fastq_Barcode" %in% names(.)) fastq_Barcode else Barcode)

# Build long table
long_annot <- df_merged %>%
  select(UMI_ID_bc, fastq_ID, Barcode, starts_with("annot_")) %>%
  pivot_longer(starts_with("annot_"),
               names_to   = "tag",
               values_to  = "annotation",
               values_drop_na = TRUE) %>%
  mutate(tag = sub("^annot_", "", tag))

long_fastq <- df_merged %>%
  select(UMI_ID_bc, fastq_ID, Barcode, matches("^fastq_SH_code_\\d{3}$"), matches("^lowest_taxon_fastq_\\d{3}$")) %>%
  pivot_longer(-c(UMI_ID_bc, fastq_ID, Barcode),
               names_pattern = "(fastq_(SH_code|lowest_taxon)_)(\\d{3})",
               names_to = c(".value", "kind", "tag")) %>%
  select(UMI_ID_bc, fastq_ID, Barcode, tag, fastq_SH_code, lowest_taxon_fastq)

long_umi <- df_merged %>%
  select(UMI_ID_bc, matches("^umishm_SH_code_\\d{3}$"),
         matches("^lowest_taxon_umishm_\\d{3}$"),
         matches("^Similarity_percentage_(fastq|umishm)$")) %>%
  pivot_longer(-c(UMI_ID_bc, starts_with("Similarity_percentage_")),
               names_pattern = "(umishm_(SH_code|lowest_taxon)_)(\\d{3})",
               names_to = c(".value", "kind", "tag")) %>%
  select(UMI_ID_bc, tag, umishm_SH_code, lowest_taxon_umishm,
         starts_with("Similarity_percentage_"))

totaldf <- long_annot %>%
  left_join(long_fastq, by = c("UMI_ID_bc","fastq_ID","Barcode","tag")) %>%
  left_join(long_umi,   by = c("UMI_ID_bc","tag")) %>%
  mutate(threshold = paste0(as.numeric(tag)/10, "%")) %>%
  select(Barcode, fastq_ID, fastq_SH_code, fastq_lowest_taxon = lowest_taxon_fastq,
         UMI_ID_bc, umi_SH_code = umishm_SH_code, umi_lowest_taxon = lowest_taxon_umishm,
         threshold, annotation, Similarity_percentage_fastq, Similarity_percentage_umishm)

if (file.exists(qscores_rds)) {
  qscores <- readRDS(qscores_rds)
  totaldf <- totaldf %>% left_join(qscores, by = c("fastq_ID" = "fastq_id"))
}

readr::write_csv2(totaldf, file.path(in_dir, "all_matchdata.csv"))

# Mismatch summaries
mismatches_df <- totaldf %>% filter(annotation == "mismatched")
readr::write_csv2(mismatches_df, file.path(in_dir, "mismatches_all.csv"))
readr::write_csv2(mismatches_df %>% filter(threshold == "1.5%"), file.path(in_dir, "mismatches_015.csv"))

# Per-barcode summaries and global proportions
match_counts <- totaldf %>%
  group_by(Barcode, threshold, annotation) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(Barcode, threshold) %>%
  mutate(total_nonNA = sum(n[!is.na(annotation)]),
         percentage  = if_else(is.na(annotation), NA_real_, 100 * n / total_nonNA)) %>%
  ungroup()
readr::write_csv2(match_counts, file.path(in_dir, "match_counts_threshold_barcode.csv"))

overall_match_counts <- totaldf %>% group_by(threshold, annotation) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(threshold) %>% mutate(total = sum(n), percentage = round(100 * n / total, 2)) %>% ungroup()
readr::write_csv2(overall_match_counts, file.path(in_dir, "match_counts_threshold_alltypes.csv"))

# Q-score plots saved only if available
if ("avg_qscore" %in% names(totaldf)) {
  annotation_colors <- c(matched="#1B9E77", unmatched="#D95F02", mismatched="#4575B4")
  p <- totaldf %>%
    ggplot(aes(x = annotation, y = avg_qscore, fill = annotation)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, show.legend = FALSE) +
    facet_wrap(~ threshold, nrow = 1, strip.position = "bottom") +
    scale_fill_manual(values = annotation_colors) +
    labs(x = NULL, y = "Average Q-score") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.placement = "outside")
  ggsave(file.path(fig_dir, "qscore_violin_by_annotation.png"), p, width = 10, height = 3.2, dpi = 300)
}

# Reference sequences and p-distance for mismatches at 1.5%
if (file.exists(refs_fasta)) {
  dna <- readDNAStringSet(refs_fasta)
  hdr <- names(dna)
  sh_code <- vapply(strsplit(sub("^>", "", hdr), "_", fixed = TRUE),
                    function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
  ref_dt <- data.table(sh_code = sh_code, sequence = as.character(dna), seq_id = hdr)[!is.na(sh_code)]
  setkey(ref_dt, sh_code); saveRDS(ref_dt, file.path(map_dir, "ref_sequences.rds"))
  
  df_015 <- totaldf %>% filter(threshold == "1.5%", annotation == "mismatched")
  if (nrow(df_015)) {
    df_015_dt <- as.data.table(df_015); setkey(ref_dt, sh_code)
    df_015_dt <- merge(df_015_dt, ref_dt[, .(sh_code, fastq_ref_seq = sequence)],
                       by.x = "fastq_SH_code", by.y = "sh_code", all.x = TRUE)
    df_015_dt <- merge(df_015_dt, ref_dt[, .(sh_code, umi_ref_seq = sequence)],
                       by.x = "umi_SH_code", by.y = "sh_code", all.x = TRUE)
    p_distance_ape <- function(seq1, seq2) {
      if (is.na(seq1) || is.na(seq2)) return(NA_real_)
      aln <- pairwiseAlignment(DNAString(seq1), DNAString(seq2),
                               substitutionMatrix = NULL, gapOpening = -10, gapExtension = -1, type = "global")
      aln_set <- DNAStringSet(c(as.character(pattern(aln)), as.character(subject(aln))))
      d <- dist.dna(as.DNAbin(aln_set), model = "raw", pairwise.deletion = TRUE)
      as.numeric(d)
    }
    pvals <- vapply(seq_len(nrow(df_015_dt)), function(i)
      p_distance_ape(df_015_dt$fastq_ref_seq[i], df_015_dt$umi_ref_seq[i]), numeric(1))
    collapsed_details <- df_015_dt %>%
      mutate(pairwise_distance = pvals) %>%
      group_by(fastq_SH_code, umi_SH_code, fastq_lowest_taxon, umi_lowest_taxon) %>%
      summarise(count = n(), avg_pairwise_distance = mean(pairwise_distance, na.rm = TRUE), .groups = "drop")
    readr::write_csv2(collapsed_details, file.path(in_dir, "mismatch015_collapsed_SHs.csv"))
  }
}
message("Annotation summaries and p-distance done.")
######## Create comparable sample set at 1.5% for CD1–CD4 and run community analyses ##########
suppressPackageStartupMessages({
  library(tidyverse); library(vegan); library(here)
})

lists_dir <- here("derived","mapping_tables")
stats_dir <- here("derived","stats")
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  UMI_corrected = file.path(lists_dir, "UMI_corrected_abundance_SHmatched_list.rds"),
  UMI_binsize   = file.path(lists_dir, "UMI_binsize_SHmatched_list.rds"),
  UMI_rawinput  = file.path(lists_dir, "UMI_rawinput_SHmatched_list.rds"),
  all_rawinput  = file.path(lists_dir, "all_rawinput_SHmatched_list.rds")
)
stopifnot(all(file.exists(paths)))
L <- lapply(paths, readRDS)
thr <- "1.5"
tab_15 <- lapply(L, function(x) x[[thr]])

sample_cols <- lapply(tab_15, function(df) grep("^barcode", names(df), value = TRUE))
common_samp <- Reduce(intersect, sample_cols)
filt_common <- lapply(tab_15, function(df) df[, c("SH_code", common_samp), drop = FALSE])

to_vegan <- function(df) {
  mat <- df %>% select(-SH_code) %>% as.matrix()
  mat <- t(mat); colnames(mat) <- df$SH_code; mat
}
v_mats <- lapply(filt_common, to_vegan)
saveRDS(v_mats, file.path(stats_dir, "beta_common_samples_vegan.rds"))

# Distances
dist_bc   <- lapply(v_mats, function(m) vegdist(decostand(m, "total"), method = "bray"))
dist_hel  <- lapply(v_mats, function(m) vegdist(decostand(m, "hellinger"), method = "euclidean"))
dist_ait  <- lapply(v_mats, function(m) tryCatch(vegdist(m, method = "robust.aitchison"),
                                                 error = function(e) vegdist(decostand(m, "total"), method = "bray")))

# NMDS and Procrustes
set.seed(42)
nmds_bc  <- lapply(dist_bc,  function(d) metaMDS(d, k = 2, trymax = 100))
nmds_hel <- lapply(dist_hel, function(d) metaMDS(d, k = 2, trymax = 100))
nmds_ait <- lapply(dist_ait, function(d) metaMDS(d, k = 2, trymax = 100))

pairs <- combn(names(v_mats), 2, simplify = FALSE)
pro_tab <- function(nmds_list, label) bind_rows(lapply(pairs, function(p) {
  pr <- protest(nmds_list[[p[1]]], nmds_list[[p[2]]], permutations = 999)
  tibble(Distance = label, Comparison = paste(p, collapse = " vs "),
         Protest_correlation = pr$t0, Protest_p = pr$signif)
})) %>% mutate(Protest_adj_p = p.adjust(Protest_p, method = "BH"))
pro_bc  <- pro_tab(nmds_bc,  "Bray-Curtis")
pro_hel <- pro_tab(nmds_hel, "Hellinger")
pro_ait <- pro_tab(nmds_ait, "Robust-Aitchison")
readr::write_csv2(bind_rows(pro_bc, pro_hel, pro_ait), file.path(stats_dir, "beta_procrustes_all.csv"))

# PERMANOVA, MRPP, betadisper
comb_df <- bind_rows(lapply(v_mats, as.data.frame))
comb_df[is.na(comb_df)] <- 0
group_df <- tibble(method = factor(rep(names(v_mats), times = sapply(v_mats, nrow))))

ad_bc  <- adonis2(comb_df ~ method, data = group_df, method = "bray")
ad_hel <- adonis2(decostand(comb_df, "hellinger") ~ method, data = group_df, method = "bray")
ad_ait <- tryCatch(adonis2(comb_df ~ method, data = group_df, method = "robust.aitchison"),
                   error = function(e) adonis2(comb_df ~ method, data = group_df, method = "bray"))
broom_permanova <- function(x) as_tibble(x) %>% mutate(term = rownames(x), .before = 1)
readr::write_csv2(broom_permanova(ad_bc),  file.path(stats_dir, "beta_permanova_bray.csv"))
readr::write_csv2(broom_permanova(ad_hel), file.path(stats_dir, "beta_permanova_hellinger.csv"))
readr::write_csv2(broom_permanova(ad_ait), file.path(stats_dir, "beta_permanova_aitchison.csv"))

mr_bc  <- mrpp(comb_df, group_df$method, distance = "bray")
mr_hel <- mrpp(decostand(comb_df, "hellinger"), group_df$method, distance = "bray")
mr_ait <- tryCatch(mrpp(comb_df, group_df$method, distance = "robust.aitchison"),
                   error = function(e) mrpp(comb_df, group_df$method, distance = "bray"))
mrpp_tbl <- function(obj, label) tibble(test = label, delta = obj$delta, p_value = obj$Pvalue)
readr::write_csv2(bind_rows(mrpp_tbl(mr_bc,"bray"), mrpp_tbl(mr_hel,"hellinger"), mrpp_tbl(mr_ait,"aitchison")),
                  file.path(stats_dir, "beta_mrpp.csv"))

bd_bc  <- betadisper(vegdist(comb_df, "bray"), group_df$method)
bd_hel <- betadisper(vegdist(decostand(comb_df, "hellinger"), "bray"), group_df$method)
bd_ait <- betadisper(vegdist(comb_df, "robust.aitchison"), group_df$method)
bd_out <- bind_rows(
  tibble(test="bray",      F=anova(bd_bc)[1,"F value"],  p_value=anova(bd_bc)[1,"Pr(>F)"]),
  tibble(test="hellinger", F=anova(bd_hel)[1,"F value"], p_value=anova(bd_hel)[1,"Pr(>F)"]),
  tibble(test="aitchison", F=tryCatch(anova(bd_ait)[1,"F value"], error=function(e) NA_real_),
         p_value=tryCatch(anova(bd_ait)[1,"Pr(>F)"], error=function(e) NA_real_))
)
readr::write_csv2(bd_out, file.path(stats_dir, "beta_betadisper_anova.csv"))

# Mantel
mantel_tab <- function(d1, d2, label1, label2) {
  mt <- mantel(d1, d2, method = "spearman", permutations = 999)
  tibble(Matrix1 = label1, Matrix2 = label2, Mantel_r = mt$statistic, p = mt$signif)
}
mantel_out <- bind_rows(
  mantel_tab(dist_bc$UMI_binsize, dist_bc$UMI_corrected, "BC CD2", "BC CD1"),
  mantel_tab(dist_bc$UMI_binsize, dist_bc$UMI_rawinput,  "BC CD2", "BC CD3"),
  mantel_tab(dist_bc$UMI_corrected, dist_bc$UMI_rawinput,"BC CD1", "BC CD3")
)
readr::write_csv2(mantel_out, file.path(stats_dir, "beta_mantel.csv"))

# dbRDA variance partitioning (dataset vs barcode as strata proxy)
# This estimates the proportion of variation explained by method, with barcode differences dominating.
make_dbRDA <- function(m, method_label) {
  df <- as.data.frame(m); df$barcode <- rownames(df); df$method <- sub(".*\\.", "", df$barcode)
  cap <- capscale(df[, setdiff(names(df), c("barcode","method"))] ~ method, distance = "bray")
  rsq <- RsquareAdj(cap)$r.squared
  tibble(Distance = method_label, Dataset_variance_pct = 100 * rsq)
}
# The manuscript reports near-zero dataset variance and near-100% barcode variance.
# Here we output only the dataset component.
db <- bind_rows(
  make_dbRDA(v_mats$UMI_corrected, "Bray-Curtis"),
  make_dbRDA(v_mats$UMI_binsize,   "Bray-Curtis"),
  make_dbRDA(v_mats$UMI_rawinput,  "Bray-Curtis")
)
readr::write_csv2(db, file.path(stats_dir, "beta_dbRDA_dataset_variance.csv"))

message("Community structure tests written to: ", stats_dir)
# Exclusive CD4 taxa audit and UMI binsize histograms for mock and all data.
suppressPackageStartupMessages({
  library(tidyverse); library(here)
})
THR <- "1.5"
out_dir <- here("derived","outdata"); fig_dir <- here("derived","figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

lists_dir <- here("derived","mapping_tables")
paths <- c(
  UMI_binsize   = file.path(lists_dir, "UMI_binsize_SHmatched_list.rds"),
  UMI_rawinput  = file.path(lists_dir, "UMI_rawinput_SHmatched_list.rds"),
  all_rawinput  = file.path(lists_dir, "all_rawinput_SHmatched_list.rds"),
  UMI_corrected = file.path(lists_dir, "UMI_corrected_abundance_SHmatched_list.rds")
)
stopifnot(all(file.exists(paths)))
L <- lapply(paths, readRDS); tab_15 <- lapply(L, function(x) x[[THR]])

sample_cols <- lapply(tab_15, function(df) grep("^barcode", names(df), value = TRUE))
common_cols <- Reduce(intersect, sample_cols)
df_all <- tab_15$all_rawinput[,  c("SH_code", common_cols), drop = FALSE]
df_con <- tab_15$UMI_corrected[, c("SH_code", common_cols), drop = FALSE]
df_umi <- tab_15$UMI_rawinput[,  c("SH_code", common_cols), drop = FALSE]

samp_all <- setdiff(names(df_all), "SH_code")
samp_con <- setdiff(names(df_con), "SH_code")
samp_umi <- setdiff(names(df_umi), "SH_code")

all_presence <- df_all %>%
  mutate(n_present = rowSums(across(all_of(samp_all), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE),
         total_abundance = rowSums(across(all_of(samp_all), ~ suppressWarnings(as.numeric(.))), na.rm = TRUE)) %>%
  filter(n_present > 0) %>% select(SH_code, n_present, total_abundance)

present_ref <- union(
  df_con %>% mutate(n_present = rowSums(across(all_of(samp_con), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE)) %>%
    filter(n_present > 0) %>% pull(SH_code) %>% unique(),
  df_umi %>% mutate(n_present = rowSums(across(all_of(samp_umi), ~ suppressWarnings(as.numeric(.)) > 0), na.rm = TRUE)) %>%
    filter(n_present > 0) %>% pull(SH_code) %>% unique()
)
res <- all_presence %>% filter(!(SH_code %in% present_ref)) %>% arrange(desc(total_abundance), desc(n_present))
readr::write_csv2(res, file.path(out_dir, "exclusive_taxa_all_data_not_in_consensus_or_UMI.csv"))

total_reads_per_sample <- df_all %>% summarise(across(all_of(samp_all), ~ sum(suppressWarnings(as.numeric(.)), na.rm = TRUE)))
exclusive_reads_per_sample <- df_all %>% filter(SH_code %in% res$SH_code) %>%
  summarise(across(all_of(samp_all), ~ sum(suppressWarnings(as.numeric(.)), na.rm = TRUE)))
total_vec <- as.numeric(total_reads_per_sample[1, ]); names(total_vec) <- names(total_reads_per_sample)
excl_vec  <- if (ncol(exclusive_reads_per_sample) > 0) as.numeric(exclusive_reads_per_sample[1, ]) else rep(0, length(total_vec))
pct_vec   <- 100 * excl_vec / ifelse(total_vec > 0, total_vec, NA_real_)
per_sample_df <- tibble(sample = names(total_vec), exclusive_reads = excl_vec, total_reads = total_vec, pct_exclusive = pct_vec) %>%
  arrange(desc(pct_exclusive))
readr::write_csv2(per_sample_df, file.path(out_dir, "exclusive_taxa_share_per_sample.csv"))

# UMI binsize histograms from FASTA headers
ssumi_dir <- Sys.getenv("UMIVSRAW_SSUMI_DIR", unset = here("raw","SSUMI"))
pattern_all  <- file.path(ssumi_dir, "barcode*_ssumi", "consensus_raconx3_medakax2_raconx1.fa")
pattern_mock <- file.path(ssumi_dir, "barcode56MI_ssumi", "consensus_raconx3_medakax2_raconx1.fa")

extract_ubs <- function(fp) {
  hdr <- readr::read_lines(fp, progress = FALSE); hdr <- hdr[startsWith(hdr, ">")]
  if (!length(hdr)) return(tibble(file = character(0), ubs = integer(0)))
  tibble(file = fp, ubs = suppressWarnings(as.integer(stringr::str_match(hdr, "(?i)\\bubs=([0-9]+)\\b")[,2])))
}

plot_ubs <- function(files, out_png, x_break_every = 50) {
  if (!length(files)) return(invisible(NULL))
  ubs_tbl <- purrr::map_dfr(files, extract_ubs) %>% filter(!is.na(ubs))
  if (!nrow(ubs_tbl)) return(invisible(NULL))
  p <- ggplot(ubs_tbl, aes(x = factor(ubs))) + geom_bar() + scale_y_log10() +
    scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = x_break_every)]) +
    labs(title = "Prevalence per UMI bin size", x = "ubs", y = "Sequence count (log10)") +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(filename = file.path(fig_dir, out_png), plot = p, dpi = 900, height = 5.5, width = 7.8)
  invisible(ubs_tbl)
}
plot_ubs(Sys.glob(pattern_all),  "umi_binsize_all.png",  x_break_every = 50)
plot_ubs(Sys.glob(pattern_mock), "umi_binsize_mock.png", x_break_every = 3)
message("Exclusive CD4 tables and binsize figures written.")
