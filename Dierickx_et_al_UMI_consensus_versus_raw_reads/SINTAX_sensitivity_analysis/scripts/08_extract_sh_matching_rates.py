#!/usr/bin/env python3
"""Recover six-threshold UNITE SH assignment rates for the current A1-filtered 90 barcodes.

Python standard library only. Reads ORIGINAL SH matching outputs, not SINTAX data.
Run after scripts/06_classification_success.py.
By default the archived SH files are under UMIvsRAW (parent of REVIEW):
  UMI_consensus_SHm/source_XX/matches/matches_out_all.csv       (D1)
  UMI_rawinput_SHm/raw_umi_barcodeXX_P1_raw_umis/*.tsv          (D3)
  all_rawinput_SHm/parallel_barcodeXXP1/*.tsv                  (D4)
Override the original SH root if necessary:
  ORIGINAL_SH_ROOT=/path/to/UMIvsRAW python3 scripts/08_extract_sh_matching_rates.py
Override only D4 if its files were moved:
  D4_SH_ROOT=/path/to/D4_original_SH_outputs python3 scripts/08_extract_sh_matching_rates.py
If D4 files cannot be identified from their paths, make a TAB-delimited manifest
with columns Barcode and path, and set D4_SH_MANIFEST=/path/to/manifest.tsv.

SH matching output `matches_out_all.csv` contains six `SH code (x)` columns.
Only codes SH########.10FU count as assignments to existing UNITE SHs.
"""
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

review_env = os.environ.get('REVIEW_DIR')
if not review_env:
    raise SystemExit('ERROR: set REVIEW_DIR to the UMIvsRAW/REVIEW directory')
ROOT = Path(review_env).resolve()
ORIGINAL_SH_ROOT = Path(os.environ.get('ORIGINAL_SH_ROOT', str(ROOT.parent))).resolve()
OUT = ROOT / 'SINTAX_equal_ITS_v1/analysis/appendixC_current/figure_C3'
OUT.mkdir(parents=True, exist_ok=True)
COHORT = ROOT / 'SINTAX_equal_ITS_v1/analysis/appendixC_current/cohorts/A1filtered90.txt'
AUDIT = ROOT / 'annotation_015_rebuilt/annotation_015_complete_barcode_summary.tsv'
if not AUDIT.exists():
    AUDIT = ROOT / 'annotation_015_complete_barcode_summary.tsv'
THRESHOLDS = ('0.5', '1.0', '1.5', '2.0', '2.5', '3.0')
SH_CODE = re.compile(r'^SH\d{7}\.10FU$', re.IGNORECASE)
BARCODE_PATTERNS = (
    re.compile(r'barcode[_-]?(\d{1,2})[_-]?(P1|MI)(?![a-z0-9])', re.I),
    re.compile(r'(?<![a-z0-9])(\d{2})(P1|MI)(?![a-z0-9])', re.I),
)

if not COHORT.is_file():
    sys.exit(f'ERROR: Missing validated cohort list: {COHORT}\nFirst run scripts/06_classification_success.py')
cohort = sorted({s.strip() for s in COHORT.read_text().splitlines() if s.strip()})
if len(cohort) != 90:
    sys.exit('ERROR: Cohort must comprise exactly the 90 current A1-filtered barcodes.')
if not AUDIT.is_file():
    sys.exit(f'ERROR: Missing original SH comparison summary: {AUDIT}')
with AUDIT.open(newline='') as f:
    rows = list(csv.DictReader(f, delimiter='\t'))
ref = {r['Barcode']: r for r in rows if r['Barcode'] in cohort}
if set(ref) != set(cohort):
    sys.exit(f'ERROR: Annotation summary and cohort disagree: {sorted(set(cohort)-set(ref))}')

sources = {'D1': {}, 'D3': {}, 'D4': {}}
for bc in cohort:
    # The archived audit paths are relative to UMIvsRAW, NOT to REVIEW.
    sources['D1'][bc] = ORIGINAL_SH_ROOT / ref[bc]['consensus_file']
    sources['D3'][bc] = ORIGINAL_SH_ROOT / ref[bc]['raw_file']
print(f'REVIEW folder:          {ROOT}', flush=True)
print(f'Original SH data root:  {ORIGINAL_SH_ROOT}', flush=True)
print(f'D1 example: {sources["D1"][cohort[0]]}', flush=True)
print(f'D3 example: {sources["D3"][cohort[0]]}', flush=True)

def barcode_from_path(path):
    s = str(path)
    for p in BARCODE_PATTERNS:
        m = p.search(s)
        if m:
            return f'{int(m.group(1)):02d}{m.group(2).upper()}'
    return None

def is_sh_output(path):
    n = path.name.lower()
    return n in ('matches_out_all.csv', 'matches_out_all.tsv') or n.endswith('_all_matches_out.tsv')

manifest = os.environ.get('D4_SH_MANIFEST')
if manifest:
    with open(manifest, newline='') as f:
        dr = csv.DictReader(f, delimiter='\t')
        if not {'Barcode', 'path'}.issubset(dr.fieldnames or []):
            sys.exit('ERROR: D4_SH_MANIFEST requires TAB-delimited Barcode and path columns.')
        for r in dr:
            bc = r['Barcode'].strip()
            if bc in cohort:
                if bc in sources['D4']:
                    sys.exit(f'ERROR: Duplicate D4 manifest barcode: {bc}')
                sources['D4'][bc] = Path(r['path'].strip())
else:
    # D4 outputs use a different layout from D1/D3. User-confirmed example:
    # UMIvsRAW/all_rawinput_SHm/parallel_barcode25P1/barcode25P1_all_matches_out.tsv
    d4root = Path(os.environ.get('D4_SH_ROOT', str(ORIGINAL_SH_ROOT / 'all_rawinput_SHm'))).resolve()
    if not d4root.is_dir():
        sys.exit(f'ERROR: D4 original SH folder missing: {d4root}\n'
                 'Set D4_SH_ROOT to the actual all_rawinput_SHm directory.')
    print(f'D4 original SH root:  {d4root}', flush=True)
    matches = defaultdict(list)
    # Use the direct known directory structure first (avoid slow recursive
    # traversal of all original data, which may contain large FASTQ files).
    for bc in cohort:
        path = d4root / f'parallel_barcode{bc}' / f'barcode{bc}_all_matches_out.tsv'
        if path.is_file():
            matches[bc].append(path)
    # For unexpected directory naming, search ONLY within the original D4
    # folder. Do not mix in any UMI-associated raw-read SH files.
    if len(matches) != len(cohort):
        print(f'Direct D4 layout found {len(matches)}/90 files. '
              'Searching remaining barcodes within D4 folder...', flush=True)
        missing_bc = set(cohort) - set(matches)
        for p in d4root.rglob('*_all_matches_out.tsv'):
            bc = barcode_from_path(p)
            if bc in missing_bc:
                matches[bc].append(p)
        for p in d4root.rglob('matches_out_all.csv'):
            bc = barcode_from_path(p)
            if bc in missing_bc:
                matches[bc].append(p)
    duplicates = {bc: ps for bc, ps in matches.items() if len(ps) > 1}
    if duplicates:
        msg = '\n'.join(f'{bc}: ' + ' | '.join(map(str, ps[:4]))
                        for bc, ps in list(duplicates.items())[:8])
        sys.exit('ERROR: Ambiguous D4 files. Set D4_SH_MANIFEST to select '
                 'one original SH file per barcode.\n' + msg)
    sources['D4'] = {bc: ps[0] for bc, ps in matches.items()}

missing = [(ds, bc, str(sources[ds].get(bc, 'NOT FOUND')))
           for ds in ('D1', 'D3', 'D4') for bc in cohort
           if bc not in sources[ds] or not sources[ds][bc].is_file()]
if missing:
    print(f'ERROR: {len(missing)} original SH output files could not be located.', file=sys.stderr)
    for ds, bc, p in missing[:14]:
        print(f'  {ds} {bc}: {p}', file=sys.stderr)
    if any(ds == 'D4' for ds, _, _ in missing):
        print('Locate original D4 SH-matching outputs using:', file=sys.stderr)
        print("  find \"$ORIGINAL_SH_ROOT\" -type f \\( -name 'matches_out_all.csv' -o -name '*_all_matches_out.tsv' \\) | head -50", file=sys.stderr)
        print('Then rerun with D4_SH_ROOT=/path/to/D4_original_SH_folder.', file=sys.stderr)
        print('If outputs are not labelled by barcode, provide D4_SH_MANIFEST.', file=sys.stderr)
    sys.exit(1)

per_barcode = []
for ds in ('D1', 'D3', 'D4'):
    for i, bc in enumerate(cohort, 1):
        p = sources[ds][bc]
        with p.open('r', encoding='utf-8-sig', errors='replace', newline='') as fh:
            header = fh.readline().rstrip('\r\n')
            delim = '\t' if '\t' in header else ','
            columns = next(csv.reader([header], delimiter=delim))
            columns = [c.strip() for c in columns]
            indices = []
            for t in THRESHOLDS:
                wanted = f'SH code ({t})'
                ix = [j for j, n in enumerate(columns) if n.lower() == wanted.lower()]
                if len(ix) != 1:
                    sys.exit(f'ERROR: {p}: expected one {wanted} column, found {len(ix)}.\nHeader: {columns[:25]}')
                indices.append(ix[0])
            n = 0
            hits = [0] * len(THRESHOLDS)
            for line in csv.reader(fh, delimiter=delim):
                if not line or all(not x.strip() for x in line):
                    continue
                if len(line) <= max(indices):
                    sys.exit(f'ERROR: Short row {n+2} in {p}')
                n += 1
                for k, ix in enumerate(indices):
                    if SH_CODE.fullmatch(line[ix].strip()):
                        hits[k] += 1
        if n == 0:
            sys.exit(f'ERROR: SH file contains no sequence rows: {p}')
        if ds == 'D3' and n != int(ref[bc]['raw_sh_reads']):
            sys.exit(f'ERROR: D3 {bc}: found {n} SH-output rows, expected '
                     f"{ref[bc]['raw_sh_reads']} from original annotation audit. {p}")
        for t, hits_t in zip(THRESHOLDS, hits):
            per_barcode.append((ds, bc, t, n, hits_t, str(p)))
        if i % 20 == 0 or i == len(cohort):
            print(f'  {ds}: processed {i}/{len(cohort)} original SH outputs', flush=True)

# Save per-barcode audit, then the pooled CSV2 used by the existing base-R plotter.
per_file = OUT / 'SH_classifiability_A1filtered90_by_barcode.tsv'
with per_file.open('w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(('Dataset', 'Barcode', 'SH_threshold_pct', 'Total_sequences', 'Classified', 'SH_output_file'))
    w.writerows(per_barcode)
summary_file = OUT / 'SH_classifiability_A1filtered90.csv'
with summary_file.open('w', newline='') as f:
    w = csv.writer(f, delimiter=';')
    w.writerow(('Dataset','SH_threshold_pct','Barcodes','Total_sequences','Classified'))
    for ds in ('D1','D3','D4'):
        for t in THRESHOLDS:
            x = [r for r in per_barcode if r[0] == ds and r[2] == t]
            if len(x) != 90:
                sys.exit(f'ERROR: Incomplete {ds} {t} records')
            w.writerow((ds, t.replace('.', ','), 90, sum(r[3] for r in x), sum(r[4] for r in x)))
print('\nOriginal SH-matching counts successfully extracted.')
print(f'Summary: {summary_file}\nAudit:   {per_file}')
print('NOTE: SH-output row counts may differ from SINTAX input counts. See README cautions.')
