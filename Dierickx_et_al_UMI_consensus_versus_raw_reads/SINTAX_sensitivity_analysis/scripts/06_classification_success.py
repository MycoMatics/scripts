#!/usr/bin/env python3
"""Create validated SH classification-rate tables for current SINTAX, cohorts 97 and 90.

Reads already-generated 01_build_equal_ITS_tables.py results. No SINTAX rerun.
Writes TSV so numerical values retain a locale-independent decimal separator.
"""
import csv
import os
from collections import defaultdict
from pathlib import Path

review_env = os.environ.get('REVIEW_DIR')
if not review_env:
    raise SystemExit('ERROR: set REVIEW_DIR to the UMIvsRAW/REVIEW directory')
REVIEW = Path(review_env).resolve()
ROOT = REVIEW / 'SINTAX_equal_ITS_v1'
OUT = ROOT / 'analysis/appendixC_current'
COHORTS = OUT / 'cohorts'
COHORTS.mkdir(parents=True, exist_ok=True)
DATASETS = ('D1', 'D3', 'D4')
THRESHOLDS = (('0.70', 'SH_conf70'), ('0.80', 'SH_conf80'), ('0.90', 'SH_conf90'))

def read_tsv(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write_tsv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as handle:
        wr = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        wr.writeheader()
        wr.writerows(rows)

summary_path = ROOT / 'analysis/classification_by_barcode.tsv'
audit_path = ROOT / 'standardization_audit.tsv'
filtered_path = ROOT / 'analysis/tables/SINTAX_CD1-CD4_long_SH_conf080_A1filtered_4datasets.csv'
for path in (summary_path, audit_path, filtered_path):
    if not path.is_file():
        raise FileNotFoundError(path)

summary = read_tsv(summary_path)
audit = read_tsv(audit_path)
if len(summary) != 291 or len(audit) != 291:
    raise ValueError('Expected 291 rows in both the standardized audit and classification summary')

by_key = {}
for row in summary:
    key = (row['barcode'], row['dataset'])
    if key in by_key or key[1] not in DATASETS:
        raise ValueError(f'Duplicate or unexpected classification row {key}')
    by_key[key] = row

by_audit = {(row['barcode'], row['dataset']): row for row in audit}
if set(by_key) != set(by_audit):
    raise ValueError('Classification and original standardization audits do not have identical barcode–dataset keys')

all97 = {bc for bc, ds in by_key}
if len(all97) != 97:
    raise ValueError(f'Expected 97 barcodes, got {len(all97)}')
for bc in all97:
    if {(bc, ds) for ds in DATASETS} - set(by_key):
        raise ValueError(f'Missing D1/D3/D4 classification summary for {bc}')

for key, row in by_key.items():
    ar = by_audit[key]
    for summary_field, audit_field in (('input_reads', 'input_reads'),
                                       ('standardized_reads', 'standardized_reads'),
                                       ('unresolved_reads', 'unresolved_reads')):
        if int(row[summary_field]) != int(ar[audit_field]):
            raise ValueError(f'Inconsistent {summary_field} for {key}')
    n = int(row['input_reads'])
    m = int(row['standardized_reads'])
    hits = [int(row[col]) for _, col in THRESHOLDS]
    if not (0 <= hits[2] <= hits[1] <= hits[0] <= int(row['any_SH_rank']) <= m <= n):
        raise ValueError(f'Non-monotonic SH confidence counts for {key}')

# The final, already validated four-dataset A1-filtered community table is
# the authority for the current 90-barcode sample set.
filtered = defaultdict(set)
with filtered_path.open(newline='', encoding='utf-8-sig') as handle:
    for row in csv.DictReader(handle, delimiter=';'):
        filtered[row['dataset']].add(row['Barcode'])
if set(filtered) != {'CD1', 'CD2', 'CD3', 'CD4'}:
    raise ValueError(f'Unexpected community datasets: {sorted(filtered)}')
cohort90 = filtered['CD1']
if len(cohort90) != 90 or any(values != cohort90 for values in filtered.values()):
    raise ValueError('The A1-filtered community cohort must consist of the same 90 barcodes in all datasets')
if not cohort90 <= all97:
    raise ValueError('A1-filtered barcodes not present in the full SINTAX audit')

for label, cohort in (('all97', all97), ('A1filtered90', cohort90)):
    (COHORTS / f'{label}.txt').write_text('\n'.join(sorted(cohort)) + '\n')
    pooled = []
    per_barcode = []
    for ds in DATASETS:
        records = [by_key[(bc, ds)] for bc in sorted(cohort)]
        for cutoff, column in THRESHOLDS:
            for record in records:
                original_n = int(record['input_reads'])
                std_n = int(record['standardized_reads'])
                n = int(record[column])
                per_barcode.append(dict(Dataset=ds, Barcode=record['barcode'],
                                        Confidence=cutoff, Original_sequences=original_n,
                                        Standardized_sequences=std_n, Classified=n,
                                        Pct_original=100 * n / original_n,
                                        Pct_standardized=100 * n / std_n if std_n else '',
                                        Unresolved=int(record['unresolved_reads'])))
            originals = sum(int(rec['input_reads']) for rec in records)
            standards = sum(int(rec['standardized_reads']) for rec in records)
            classified = sum(int(rec[column]) for rec in records)
            pooled.append(dict(Dataset=ds, Confidence=cutoff,
                               Barcodes=len(cohort), Total_sequences=originals,
                               Standardized_sequences=standards,
                               Unresolved_sequences=sum(int(rec['unresolved_reads']) for rec in records),
                               Any_SH_rank=sum(int(rec['any_SH_rank']) for rec in records),
                               Classified=classified,
                               Classified_pct_original=100 * classified / originals,
                               Classified_pct_standardized=100 * classified / standards if standards else ''))
    cols = ['Dataset', 'Confidence', 'Barcodes', 'Total_sequences', 'Standardized_sequences',
            'Unresolved_sequences', 'Any_SH_rank', 'Classified', 'Classified_pct_original',
            'Classified_pct_standardized']
    path = OUT / f'Table_C2_SINTAX_{label}_pooled.tsv'
    write_tsv(path, pooled, cols)
    write_tsv(OUT / f'Table_C2_SINTAX_{label}_by_barcode.tsv', per_barcode,
              ['Dataset', 'Barcode', 'Confidence', 'Original_sequences', 'Standardized_sequences',
               'Classified', 'Pct_original', 'Pct_standardized', 'Unresolved'])
    print(f'\n===== CURRENT SINTAX CLASSIFICATION: {label} =====')
    print('Dataset\tConfidence\tTotal original\tConfident\t% original\t% standardized')
    for row in pooled:
        print(f"{row['Dataset']}\t{row['Confidence']}\t{row['Total_sequences']:,}\t"
              f"{row['Classified']:,}\t{row['Classified_pct_original']:.3f}\t"
              f"{row['Classified_pct_standardized']:.3f}")
    print('Wrote', path)
print('\nCurrent classification audits checked. Original results and SH matching unchanged.')
