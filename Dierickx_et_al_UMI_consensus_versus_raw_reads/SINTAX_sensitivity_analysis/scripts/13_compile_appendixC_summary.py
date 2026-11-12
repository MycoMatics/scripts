#!/usr/bin/env python3
"""Compile current Appendix C numeric tables and provenance into one Markdown audit.

Run after scripts 06-10, including mismatch audits for 97 and 90 barcodes. Reuses actual results,
never carries over historical SINTAX values. Original inputs are read-only.
"""
import csv
import os
from pathlib import Path

review_env = os.environ.get('REVIEW_DIR')
if not review_env:
    raise SystemExit('ERROR: set REVIEW_DIR to the UMIvsRAW/REVIEW directory')
REVIEW = Path(review_env).resolve()
SINTAX = REVIEW / 'SINTAX_equal_ITS_v1'
BASE = SINTAX / 'analysis/appendixC_current'
OUTPUT = BASE / 'AppendixC_current_numbers.md'

def read_tsv(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def fmt(value):
    try:
        x = float(value)
        return f'{x:,.3f}' if '.' in str(value) else f'{int(x):,}'
    except (ValueError, TypeError):
        return str(value)

lines = [
    '# ITS-standardized Appendix C: audited numerical results',
    '',
    'Source: VSEARCH SINTAX 2.30.0, UNITE v10, fixed bootstrap thresholds.',
    'CD2 is constructed from the current D1 SH assignments and original `ubs` family sizes.',
    'Classification denominators use original (pre-standardization) sequence counts.',
    'The primary 97-barcode classification cohort and 90-barcode A1 ecological cohort are not interchangeable.',
    '',
]

for cohort in ('all97', 'A1filtered90'):
    path = BASE / f'Table_C2_SINTAX_{cohort}_pooled.tsv'
    data = read_tsv(path)
    lines += [f'## Classification success: {cohort}', '',
              '| Dataset | Confidence | Input sequences | Confident SH | Classified (%) |',
              '|---|---:|---:|---:|---:|']
    for r in data:
        lines.append(f"| {r['Dataset']} | {r['Confidence']} | {int(r['Total_sequences']):,} | "
                     f"{int(r['Classified']):,} | {float(r['Classified_pct_original']):.3f} |")
    lines.append('')

for cohort in ('97', '90'):
    path = BASE / f'chimera_{cohort}' / 'SINTAX080_confusion_audit.tsv'
    if not path.exists():
        lines += [f'## Current chimera audit: {cohort} barcodes', '', 'Not yet run.', '']
        continue
    audits = {r['Metric']: r['Value'] for r in read_tsv(path)}
    mismatch = int(audits['Mismatched'])
    match = int(audits['Matched'])
    unmatched = int(audits['Unmatched'])
    eligible = int(audits['SINTAX D1-eligible D3 reads'])
    excluded = int(audits['D1 low-confidence excluded'])
    if eligible != match + unmatched + mismatch:
        raise ValueError('Mismatch audit does not reconcile to eligible D3 reads')
    if eligible + excluded != int(audits['Total paired D3 reads']):
        raise ValueError('Mismatch audit does not reconcile to total D3 reads')
    filtered = BASE / f'chimera_{cohort}' / 'exclude_SH0074422.10FU' / 'SINTAX080_nonchimeric_confusion_audit.tsv'
    cf = {r['Metric']: r['Value'] for r in read_tsv(filtered)}
    n_conflicts_with = int(cf['Excluded chimeric-reference conflicts'])
    nonchim = int(cf['Nonchimeric mismatches remaining'])
    if mismatch != n_conflicts_with + nonchim:
        raise ValueError('Chimeric exclusion and retained conflicts do not reconcile')
    lines += [f'## Conditional SH concordance and chimera check: {cohort} barcodes', '',
              f'- All D3 reads: **{int(audits["Total paired D3 reads"]):,}**',
              f'- Eligible D3 reads (D1 confidence ≥0.80): **{eligible:,}**',
              f'- Matched: **{match:,}** ({100*match/eligible:.3f}% of D1-eligible)',
              f'- Unmatched (D3 <0.80): **{unmatched:,}** ({100*unmatched/eligible:.3f}% of D1-eligible)',
              f'- Mismatched: **{mismatch:,}** ({100*mismatch/eligible:.3f}% of D1-eligible)',
              f'- Confident on both D1 and D3: **{match + mismatch:,}**',
              f'- Mismatches among both-confident D1–D3 comparisons: **{100*mismatch/(match+mismatch):.3f}%**',
              f'- Mismatches involving previously recognized chimeric SH0074422.10FU: '
              f'**{n_conflicts_with:,}** ({100*n_conflicts_with/mismatch if mismatch else 0:.3f}% of observed mismatches)',
              f'- Remaining observed mismatches not involving this SH: **{nonchim:,}**',
              f'- Adjusted D1-eligible mismatch rate after omitting only those observed conflict rows: '
              f'**{100*nonchim/(eligible-n_conflicts_with):.3f}%**', '',
              '**Caution:** The adjusted mismatch rate is a descriptive exclusion sensitivity, '
              'not a prediction of reclassification against a corrected reference database.', '']

    directional = BASE / f'chimera_{cohort}' / 'SINTAX080_conflicting_SH_pairs_directional.tsv'
    if directional.is_file():
        rows = read_tsv(directional)[:12]
        lines += ['| D1 consensus SH | D3 raw-read SH | Conflicting reads | % of mismatches |',
                  '|---|---|---:|---:|']
        for r in rows:
            lines.append(f"| {r['D1_consensus_SH']} | {r['D3_raw_read_SH']} | "
                         f"{int(r['Mismatched_reads']):,} | {float(r['Pct_all_SINTAX_mismatches']):.2f} |")
        lines.append('')

# Existing final community results, without recalculating or changing the model.
community = SINTAX / 'analysis/community_sensitivity/CD1-CD4_analysis_conf080_A1filtered_4datasets'
if community.is_dir():
    lines += ['## Community statistics', '',
              'Four datasets × 90 A1-filtered barcodes; these files were already computed:', '']
    for filename in (
        'SINTAX_CD1-CD4_variance_partition_conf080.csv',
        'SINTAX_CD1-CD4_PERMANOVA_global_conf080.csv',
        'SINTAX_CD1-CD4_PERMANOVA_pairwise_conf080.csv',
        'SINTAX_CD1-CD4_dispersion_tests_conf080.csv',
        'SINTAX_CD1-CD4_Mantel_Procrustes_conf080.csv',
        'SINTAX_CD1-CD4_SH_overlap_summary_conf080.csv',
    ):
        path = community / filename
        lines.append(f'- `{path}`' + ('' if path.exists() else ' **MISSING**'))
    lines.append('')

cont = BASE / 'contingency_A1filtered90/contingency_SH_vs_SINTAX080_shared_long.tsv'
if cont.is_file():
    lines += ['## Method-comparison contingency', '',
              f'The read-level SH-matching versus current SINTAX contingency is available at `{cont}`.',
              'Its shared eligible-read denominator must be reported separately from all-D3 classification success.', '']
else:
    lines += ['## Method-comparison contingency', '', 'Not yet run. Run scripts/09_compare_annotation_methods.R.', '']

sh = BASE / 'figure_C3/SH_classifiability_A1filtered90.csv'
if sh.is_file():
    lines += ['## SH-matching versus SINTAX classification thresholds', '',
              f'Original SH-matching six-threshold counts for the same 90 barcodes: `{sh}`.',
              'This comparison uses method-specific sequence denominators, which must be disclosed.', '']
else:
    lines += ['## SH-matching versus SINTAX classification thresholds', '',
              'Not yet extracted. Run scripts/08_extract_sh_matching_rates.py.', '']

lines += ['## Methodological boundary', '',
          'The historical VSEARCH 2.25.0 analysis and current 2.30.0 analysis differ in execution '
          'configuration. Identical input sequences occasionally changed SINTAX confidence outcomes. '
          'Historical-to-current differences therefore cannot be attributed exclusively to ITS coverage.',
          'All results here use the current ITS-standardized classification.', '']
OUTPUT.write_text('\n'.join(lines), encoding='utf-8')
print('\n'.join(lines))
print(f'\nWrote: {OUTPUT}')
