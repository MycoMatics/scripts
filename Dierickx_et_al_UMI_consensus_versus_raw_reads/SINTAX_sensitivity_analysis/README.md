# SINTAX sensitivity analysis: standardized fungal ITS metabarcoding

Analysis code and supplementary community tables accompanying **“Do we need Unique Molecular Identifiers in Nanopore metabarcoding?”**.
This analysis was carried out as a sensitivity check following a first round of peer-review.

## Contents and scope

- `scripts/00_…` to `scripts/13_…`: final ITS preparation, SINTAX, D1–D3 concordance, four-community analysis and figure scripts. Filenames are numbered in analytical order. Run individual stages only when upstream results are available.
- `run_analysis.sh`: optional full workflow; executed with `bash`, never sourced. It **does not** run the reconstructed ITS-preparation stage automatically and does not set `errexit` or other shell-wide options.
- `SINTAX_sensitivity_CD1-CD4.xlsx`: supplied supplementary abundance data, **copied byte-for-byte without modification**. Sheets `CD1`, `CD2`, `CD3` and `CD4` each contain **90 A1-filtered barcodes** and use `dataset`, `Barcode`, `SH_code` and `count` columns.
- `standardization_audit.tsv`: 291 dataset-by-barcode standardization records, including original, standardized and unresolved counts.
- `A1_exclusions.tsv`: predefined sample exclusion map.
- `expected_final_values.tsv`: compact validation targets for the final analysis.

The workbook is an **output dataset**

## Data definitions

| Name | Biological/analytical representation |
|---|---|
| D1 / CD1 | One count per confidently assigned UMI consensus |
| CD2 | D1 SH assignments weighted by the original FASTA `ubs=` UMI-family size |
| D3 / CD3 | Independently classified raw reads linked to UMI consensuses |
| D4 / CD4 | Independently classified all raw reads from the same UMI-tagged libraries |

`umi289bins` is an identifier, **not** a family size. CD2 and CD1 have identical SH presence per barcode, but different abundance weights. Neither CD3 nor CD4 is a UMI-free library preparation.

Classification summaries and D1–D3 conditional read concordance use **97 barcodes**. Ecological statistics, cross-algorithm comparisons and Figure B3 use the **90 core barcodes**. Excluded barcodes in this study are `19P1`, `24P1`, `33MI`, `41MI`, `42MI`, `45MI` and `56P1`.

## Reproducibility and provenance

The final analysis used **Cutadapt 5.0**, **VSEARCH 2.30.0** (`--sintax_cutoff 0.0 --randseed 42 --strand both`) and **UNITE v10**. Species-rank SINTAX confidence thresholds were 0.70, 0.80 and 0.90; community analysis used **≥0.80**. VSEARCH queries are standardized to comparable ITS coverage. P1 sequence bases were retained; MI reads required primer-boundary handling, including a separate D1 branch for 35MI, 40MI and 41MI.

| Source dataset | Original sequences | Standardized sequences | Unresolved |
|---|---:|---:|---:|
| D1 | 993,738 | 993,564 | 174 |
| D3 | 12,695,049 | 12,695,049 | 0 |
| D4 | 44,242,882 | 44,211,759 | 31,123 |

**Classification-success denominators use the original counts.** Unresolved standardization inputs were not classified, but remain in the denominator. Raw read denominators, D1-eligible agreement denominators and shared SH-matching–SINTAX denominators are distinct.

Community analysis uses **vegan 2.7-6**, Bray–Curtis on relative abundances, Euclidean distance after Hellinger transformation and robust Aitchison distance (`vegdist(method="robust.aitchison")`). NMDS uses *k* = 3, 3 and 2, respectively, with `trymax=1000`. PERMANOVA uses 999 permutations constrained within barcode. Variance partitioning enters barcode before dataset. Dispersion, Friedman, Wilcoxon, Mantel and Procrustes analyses follow the analytical scripts. R packages: vegan, ggplot2, tidyverse, cowplot, data.table, MASS and scales. Python additionally requires `edlib`; install `cutadapt==5.0` and VSEARCH as external executables.

## Running the scripts

This folder contains **code and small metadata**, not the ~58 million original sequences, UNITE FASTA, UMI mappings or intermediate classification outputs. The executed scripts expect those original project files under `UMIvsRAW/REVIEW/SINTAX_equal_ITS_v1/`. 

```bash
export PROJECT_ROOT=/path/to/UMIvsRAW
export REVIEW_DIR="$PROJECT_ROOT/REVIEW"

# Optional: validate reconstructed query preparation against archived production queries.
python3 scripts/00_prepare_equal_ITS_queries.py \
  --project-root "$PROJECT_ROOT" \
  --output-root "$REVIEW_DIR/SINTAX_equal_ITS_reconstructed" \
  --reference-root "$REVIEW_DIR/SINTAX_equal_ITS_v1"

# Full classification and downstream analysis (only if rerunning the study).
bash run_analysis.sh
echo "Analysis runner exit status: $?"
```

For individual steps see `run_analysis.sh`. Dependencies must be installed beforehand; scripts do **not** silently install R packages. The runner reports each step's status and exits its own process on an error without altering the calling `tmux` shell.

## Independent validation targets

- SINTAX confident SH counts at 0.80: D1 **724,964**, D3 **9,292,552**, D4 **30,605,727** (97 barcodes).
- A1-filtered workbook richness, pooled: CD1/CD2 **470**, CD3 **501**, CD4 **932** SHs; **455** SHs shared among CD1, CD3 and CD4; CD4 sum **29,794,242** reads. The 421 CD4-exclusive SHs contribute **7,754** reads.
- D1–D3 at 0.80: **9,530,734** D1-eligible D3 reads, of which **9,031,008** matched, **462,496** unmatched and **37,230** mismatched. **36,875** mismatches involve putatively chimeric `SH0074422.10FU`. Excluding those observed conflict rows is descriptive, not a classifier rerun.
- Figure B3: **8,432,242** shared eligible read–consensus pairs for **90** barcodes. This is not the same denominator as the 97-barcode SINTAX comparison.
