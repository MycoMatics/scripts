#!/usr/bin/env python3
"""Summarise all standardized SINTAX results and build the CD1/CD3/CD4 community table.

Uses original SINTAX >=0.80 species-SH convention, and original (pre-trimming)
input counts as classification denominators. Does not modify SINTAX results.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import os
import re
import sys
import time

review_env = os.environ.get('REVIEW_DIR')
if not review_env:
    raise SystemExit('ERROR: set REVIEW_DIR to the UMIvsRAW/REVIEW directory')
ROOT = Path(review_env).resolve() / 'SINTAX_equal_ITS_v1' 
AUDIT = ROOT / 'standardization_audit.tsv'
RESULTS = ROOT / 'results'
OUT = ROOT / 'analysis'
TABLES = OUT / 'tables'
TABLES.mkdir(parents=True, exist_ok=True)

# Species-rank UNITE SH parser used throughout the final SINTAX analysis.
SH_RE = re.compile(r'(?:^|,)s:(SH[0-9]{7}\.10FU)\(([0-9]+(?:\.[0-9]+)?)\)(?:,|$)')
SPECIES_RE = re.compile(r'(?:^|,)s:([^,]+)')
DATASETS = {'D1': 'CD1', 'D3': 'CD3', 'D4': 'CD4'}
CUTOFFS = (0.70, 0.80, 0.90)

with AUDIT.open(newline='') as handle:
    audit_rows = list(csv.DictReader(handle, delimiter='\t'))
assert len(audit_rows) == 291, f'Expected 291 audit rows, found {len(audit_rows)}'

summary = []
counts_080 = defaultdict(Counter)
start = time.monotonic()

for i, audit in enumerate(sorted(audit_rows, key=lambda r: (r['barcode'], r['dataset'])), 1):
    bc = audit['barcode']
    ds = audit['dataset']
    input_n = int(audit['input_reads'])
    standardized_n = int(audit['standardized_reads'])
    unresolved_n = int(audit['unresolved_reads'])
    missing_uuid = int(audit['missing_uuid'])
    assert input_n == standardized_n + unresolved_n + missing_uuid, (bc, ds)
    assert ds in DATASETS, (bc, ds)

    path = RESULTS / f'{bc}.{ds}.sintax.tsv'
    if not path.is_file():
        raise FileNotFoundError(path)

    n_rows = 0
    any_taxonomy = 0
    any_SH = 0
    confidence_counts = Counter()
    observed_ids = set() if ds == 'D1' else None

    with path.open('rt', encoding='utf-8') as handle:
        for line_num, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f'Blank row in {path}:{line_num}')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 3:
                raise ValueError(f'Bad SINTAX TSV format: {path}:{line_num}')
            qid, taxonomy = fields[0], fields[1]
            if not qid.startswith(f'{ds}|{bc}|'):
                raise ValueError(f'Unexpected query ID in {path}:{line_num}: {qid}')
            if observed_ids is not None:
                if qid in observed_ids:
                    raise ValueError(f'Duplicate D1 query ID: {qid}')
                observed_ids.add(qid)

            n_rows += 1
            if taxonomy not in ('', '*', 'unclassified', 'Unclassified'):
                any_taxonomy += 1

            match = SH_RE.search(taxonomy)
            if match:
                sh = match.group(1)
                conf = float(match.group(2))
                if not (0 <= conf <= 1):
                    raise ValueError(f'Invalid confidence {conf}: {path}:{line_num}')
                any_SH += 1
                for cutoff in CUTOFFS:
                    if conf >= cutoff:
                        confidence_counts[cutoff] += 1
                if conf >= 0.80:
                    counts_080[(ds, bc)][sh] += 1
            elif SPECIES_RE.search(taxonomy):
                raise ValueError(f'Unrecognised species-level SH in {path}:{line_num}')

    if n_rows != standardized_n:
        raise ValueError(f'{bc} {ds}: SINTAX rows {n_rows:,}, expected {standardized_n:,}')

    row = dict(barcode=bc, dataset=ds, input_reads=input_n,
               standardized_reads=standardized_n, unresolved_reads=unresolved_n,
               missing_uuid=missing_uuid, any_taxonomy=any_taxonomy,
               any_SH_rank=any_SH)
    for cutoff in CUTOFFS:
        tag = str(int(cutoff * 100))
        n = confidence_counts[cutoff]
        row[f'SH_conf{tag}'] = n
        row[f'SH_conf{tag}_pct_original'] = f'{100 * n / input_n:.6f}' if input_n else 'NA'
        row[f'SH_conf{tag}_pct_standardized'] = f'{100 * n / standardized_n:.6f}' if standardized_n else 'NA'
    summary.append(row)
    if i % 10 == 0 or i == len(audit_rows):
        print(f'[{i}/{len(audit_rows)}] completed {bc} {ds}: {n_rows:,} rows', flush=True)

fields = list(summary[0])
with (OUT / 'classification_by_barcode.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(summary)

with (OUT / 'classification_pooled.tsv').open('w', newline='') as handle:
    cols = ['dataset', 'input_reads', 'standardized_reads', 'unresolved_reads',
            'any_taxonomy', 'any_SH_rank', 'SH_conf70', 'SH_conf80', 'SH_conf90',
            'SH_conf80_pct_original', 'SH_conf80_pct_standardized']
    writer = csv.DictWriter(handle, fieldnames=cols, delimiter='\t')
    writer.writeheader()
    for ds in DATASETS:
        rows = [r for r in summary if r['dataset'] == ds]
        totals = {c: sum(int(r[c]) for r in rows) for c in cols[1:9]}
        totals['dataset'] = ds
        totals['SH_conf80_pct_original'] = f"{100 * totals['SH_conf80'] / totals['input_reads']:.6f}"
        totals['SH_conf80_pct_standardized'] = f"{100 * totals['SH_conf80'] / totals['standardized_reads']:.6f}"
        writer.writerow(totals)
        print(ds, 'conf>=0.80:', f"{totals['SH_conf80']:,}/{totals['input_reads']:,}",
              f"({totals['SH_conf80_pct_original']}% original denominator)", flush=True)

# SH presence is defined by a confident assignment to the species-rank UNITE SH.
barcode_sets = {
    ds: {bc for (xds, bc), shs in counts_080.items() if xds == ds and shs}
    for ds in DATASETS
}
cohort = set.intersection(*barcode_sets.values())
all_barcodes = sorted({r['barcode'] for r in summary})
print('Barcodes with >=1 confident SH per dataset:',
      {ds: len(s) for ds, s in barcode_sets.items()})
print('Paired nonzero barcode intersection:', len(cohort), '/', len(all_barcodes))
if len(cohort) < 3:
    raise ValueError('Insufficient paired barcodes for community analysis')

with (OUT / 'paired_barcode_audit.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['barcode', 'CD1_nonzero', 'CD3_nonzero', 'CD4_nonzero', 'included'])
    for bc in all_barcodes:
        writer.writerow([bc, *[int(bc in barcode_sets[ds]) for ds in DATASETS], int(bc in cohort)])

for label, path, selected in [
    ('balanced', TABLES / 'SINTAX_CD1-CD4_long_SH_conf080.csv', cohort),
    ('all', TABLES / 'SINTAX_CD1-CD4_long_all_barcodes_conf080.csv', set(all_barcodes)),
]:
    n_written = 0
    with path.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter=';', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['dataset', 'Barcode', 'SH_code', 'count'])
        for ds, name in DATASETS.items():
            for bc in sorted(selected):
                for sh, n in sorted(counts_080.get((ds, bc), {}).items()):
                    writer.writerow([name, bc, sh, n])
                    n_written += 1
    print(label, 'community long table:', n_written, 'nonzero cells:', path, flush=True)

elapsed = (time.monotonic() - start) / 60
print(f'Finished all 291 SINTAX result files in {elapsed:.1f} minutes')
print('Original SINTAX and SH-matching data were not modified.')
