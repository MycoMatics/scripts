#!/usr/bin/env python3
"""Reconstruct equal-ITS CD2 by weighting each confidently assigned D1 SH by
its original UMI family size (the FASTA-header `ubs` attribute).

The digits in `umi289bins` are a UMI FAMILY IDENTIFIER, not a weight.
Reads current equal-ITS SINTAX results and original D1 FASTAs read-only.
Outputs a separate four-dataset table without modifying prior 3-dataset data.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import os
import re
import time

project_env = os.environ.get('PROJECT_ROOT')
if not project_env:
    raise SystemExit('ERROR: set PROJECT_ROOT to the UMIvsRAW project root')
BASE = Path(project_env).resolve()
ROOT = BASE / 'REVIEW/SINTAX_equal_ITS_v1'
SOURCE = BASE / 'UMI_consensus_SHm/input_UMI_consensus'
AUDIT = ROOT / 'standardization_audit.tsv'
RESULTS = ROOT / 'results'
TABLES = ROOT / 'analysis/tables'
TABLES.mkdir(parents=True, exist_ok=True)

THREE = TABLES / 'SINTAX_CD1-CD4_long_SH_conf080.csv'
FOUR = TABLES / 'SINTAX_CD1-CD4_long_SH_conf080_4datasets.csv'
CD2 = TABLES / 'SINTAX_CD2_consensus_uncorrected_SH_conf080_equal_ITS.csv'
AUDIT_OUT = TABLES / 'SINTAX_CD2_ubs_audit_equal_ITS.tsv'
SH_RE = re.compile(r'(?:^|,)s:(SH[0-9]{7}\.10FU)\(([0-9]+(?:\.[0-9]+)?)\)(?:,|$)')
SPECIES_RE = re.compile(r'(?:^|,)s:')
UBS_RE = re.compile(r'(?:^|[;\s])ubs=(\d+)(?=$|[;\s])')
CONFIDENCE = 0.80

with AUDIT.open(newline='') as handle:
    manifest = list(csv.DictReader(handle, delimiter='\t'))
assert len(manifest) == 291, f'Expected 291 audit rows, got {len(manifest)}'
d1_rows = {r['barcode']: r for r in manifest if r['dataset'] == 'D1'}
assert len(d1_rows) == 97, f'Expected 97 D1 barcodes, got {len(d1_rows)}'

weighted = defaultdict(Counter)
validation = []
started = time.monotonic()

for number, (barcode, row) in enumerate(sorted(d1_rows.items()), start=1):
    original = SOURCE / f'barcode{barcode}.fasta'
    result = RESULTS / f'{barcode}.D1.sintax.tsv'
    for path in (original, result):
        if not path.is_file():
            raise FileNotFoundError(f'Missing required source: {path}')

    # Read the ORIGINAL unmodified consensus headers, which retain `ubs`.
    # `umi123bins` is an identifier, NOT a family size.
    sizes = {}
    missing_metadata = set()
    duplicated_ids = set()
    fasta_records = 0
    with original.open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.startswith('>'):
                continue
            fasta_records += 1
            header = line[1:].rstrip('\r\n')
            uid = header.split(';', 1)[0].split()[0]
            if not uid:
                raise ValueError(f'Empty FASTA ID at {original}:{line_no}')
            if uid in sizes or uid in missing_metadata:
                duplicated_ids.add(uid)
                continue
            matches = UBS_RE.findall(header)
            if len(matches) == 0:
                missing_metadata.add(uid)
            elif len(matches) != 1:
                raise ValueError(f'Ambiguous ubs field for {barcode} {uid}: {header[:200]}')
            else:
                family_size = int(matches[0])
                if family_size <= 0:
                    raise ValueError(f'Nonpositive ubs for {barcode} {uid}')
                sizes[uid] = family_size

    if fasta_records != int(row['input_reads']):
        raise ValueError(f'{barcode}: original FASTA has {fasta_records} records, '
                         f"audit expects {row['input_reads']}")
    if duplicated_ids:
        raise ValueError(f'{barcode}: {len(duplicated_ids)} duplicate source UMI IDs. '
                         f'Examples: {sorted(duplicated_ids)[:6]}')

    queried = 0
    confident = 0
    confidently_missing_ubs = set()
    missing_source_ids = set()
    seen = set()
    n_larger_than_one = 0
    sum_ubs_all = 0
    sum_ubs_confident = 0
    with result.open(encoding='utf-8') as handle:
        for line_no, line in enumerate(handle, 1):
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 3:
                raise ValueError(f'Malformed SINTAX row at {result}:{line_no}')
            prefix = f'D1|{barcode}|'
            if not fields[0].startswith(prefix):
                raise ValueError(f'Unexpected SINTAX ID: {fields[0]} in {result}')
            uid = fields[0][len(prefix):]
            if uid in seen:
                raise ValueError(f'Duplicate SINTAX D1 ID: {barcode} {uid}')
            seen.add(uid)
            queried += 1

            if uid not in sizes and uid not in missing_metadata:
                missing_source_ids.add(uid)
            else:
                if uid in sizes:
                    sum_ubs_all += sizes[uid]

            match = SH_RE.search(fields[1])
            if match:
                sh = match.group(1)
                confidence = float(match.group(2))
                if not (0 <= confidence <= 1):
                    raise ValueError(f'Invalid confidence for {barcode} {uid}')
                if confidence >= CONFIDENCE:
                    confident += 1
                    if uid not in sizes:
                        confidently_missing_ubs.add(uid)
                    else:
                        ubs = sizes[uid]
                        weighted[barcode][sh] += ubs
                        sum_ubs_confident += ubs
                        n_larger_than_one += ubs > 1
            elif SPECIES_RE.search(fields[1]):
                raise ValueError(f'Unrecognised SH annotation for {barcode} {uid}')

    if queried != int(row['standardized_reads']):
        raise ValueError(f'{barcode}: {queried} SINTAX records vs '
                         f"{row['standardized_reads']} standardized reads")
    if missing_source_ids or confidently_missing_ubs:
        raise ValueError(
            f'{barcode}: missing UMI header data. '
            f'Absent original IDs={len(missing_source_ids)} examples={sorted(missing_source_ids)[:5]}; '
            f'missing ubs for confident D1={len(confidently_missing_ubs)} '
            f'examples={sorted(confidently_missing_ubs)[:5]}. '
            'Do not substitute the number in umiNNN bins or assume ubs=1.'
        )
    if sum(weighted[barcode].values()) != sum_ubs_confident:
        raise ValueError(f'{barcode}: weighted SH counts do not sum to ubs total')
    if confident > sum_ubs_confident:
        raise ValueError(f'{barcode}: sum of family sizes lower than confident D1 count')

    validation.append(dict(
        Barcode=barcode, original_D1_FASTA_records=fasta_records,
        standardized_D1_SINTAX_rows=queried,
        D1_SH_conf080=confident,
        source_D1_missing_ubs_total=len(missing_metadata),
        source_D1_missing_ubs_conf080=len(confidently_missing_ubs),
        weighted_CD2_SH_conf080=sum_ubs_confident,
        confident_D1_with_ubs_gt_1=n_larger_than_one,
        n_distinct_CD2_SH=len(weighted[barcode]),
        sum_ubs_all_standardized_with_metadata=sum_ubs_all
    ))
    if number % 10 == 0 or number == len(d1_rows):
        print(f'[{number}/{len(d1_rows)}] {barcode}: '
              f'{confident:,} confident D1 -> {sum_ubs_confident:,} weighted CD2', flush=True)

# Existing 3-dataset community table must be intact, and data must be balanced.
if not THREE.is_file():
    raise FileNotFoundError(THREE)
previous = defaultdict(Counter)
source_rows = []
with THREE.open(newline='', encoding='utf-8') as handle:
    rd = csv.DictReader(handle, delimiter=';')
    if rd.fieldnames != ['dataset', 'Barcode', 'SH_code', 'count']:
        raise ValueError(f'Unexpected 3-dataset CSV2 schema: {rd.fieldnames}')
    for r in rd:
        if r['dataset'] not in ('CD1', 'CD3', 'CD4'):
            raise ValueError(f'Invalid dataset in 3-dataset table: {r}')
        n = int(r['count'])
        if n <= 0:
            raise ValueError('Expected strictly positive counts in 3-dataset table')
        previous[(r['dataset'], r['Barcode'])][r['SH_code']] += n
        source_rows.append((r['dataset'], r['Barcode'], r['SH_code'], n))

barcodes = sorted(d1_rows)
for barcode in barcodes:
    cd1 = previous['CD1', barcode]
    cd2 = weighted[barcode]
    if not cd1 or not cd2 or not previous['CD3', barcode] or not previous['CD4', barcode]:
        raise ValueError(f'{barcode}: one or more community datasets lack confident SHs')
    if set(cd1) != set(cd2):
        raise ValueError(f'{barcode}: CD1 and CD2 have differing SH sets. '
                         'This should be impossible with positive ubs weights.')
    # Strong validation: successful SH assignments are identical at CD1 and CD2
    # and only their read-abundance weights differ.
    if any(cd2[sh] < cd1[sh] for sh in cd1):
        raise ValueError(f'{barcode}: a CD2 SH has fewer reads than its CD1 counterpart')
    if sum(cd1.values()) != next(r['D1_SH_conf080'] for r in validation if r['Barcode'] == barcode):
        raise ValueError(f'{barcode}: unexpected CD1 count in existing community table')

# Write to temporary files first; do not overwrite completed 3-dataset artifacts.
for dest in (FOUR, CD2, AUDIT_OUT):
    if dest.exists():
        raise FileExistsError(f'Output exists, refusing overwrite: {dest}')

outputs = [FOUR, CD2, AUDIT_OUT]
temps = [p.with_name(p.name + f'.tmp.{os.getpid()}') for p in outputs]
try:
    with temps[0].open('w', newline='', encoding='utf-8') as combined, \
         temps[1].open('w', newline='', encoding='utf-8') as cd2_only, \
         temps[2].open('w', newline='', encoding='utf-8') as audit_out:
        w4 = csv.writer(combined, delimiter=';')
        w2 = csv.writer(cd2_only, delimiter=';')
        wa = csv.DictWriter(audit_out, fieldnames=list(validation[0]), delimiter='\t')
        w4.writerow(['dataset', 'Barcode', 'SH_code', 'count'])
        w2.writerow(['dataset', 'Barcode', 'SH_code', 'count'])
        wa.writeheader()
        wa.writerows(validation)
        for dataset, bc, sh, n in source_rows:
            w4.writerow([dataset, bc, sh, n])
        for bc in barcodes:
            for sh, n in sorted(weighted[bc].items()):
                w4.writerow(['CD2', bc, sh, n])
                w2.writerow(['CD2', bc, sh, n])
    for temp, output in zip(temps, outputs):
        os.replace(temp, output)
finally:
    for temp in temps:
        if temp.exists():
            temp.unlink()

print('\n===== EQUAL-ITS CD2 RECONSTRUCTION COMPLETE =====')
print('Barcodes:', len(barcodes), '| Datasets: CD1, CD2, CD3, CD4')
print('Confident D1 consensus sequences:', f"{sum(r['D1_SH_conf080'] for r in validation):,}")
print('CD2 weighted classified reads:', f"{sum(r['weighted_CD2_SH_conf080'] for r in validation):,}")
print('Original D1 FASTA records lacking ubs (all confidence levels):',
      sum(r['source_D1_missing_ubs_total'] for r in validation))
print('Missing ubs among confidently classified D1 sequences: 0 (validated)')
print('Four-dataset table:', FOUR)
print('CD2-only table:', CD2)
print('Audit:', AUDIT_OUT)
print('Elapsed minutes:', round((time.monotonic()-started)/60, 2))
print('Historical SINTAX and original SH-matching files unchanged.')
