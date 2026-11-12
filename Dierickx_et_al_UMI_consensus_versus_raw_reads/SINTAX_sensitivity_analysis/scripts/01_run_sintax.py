#!/usr/bin/env python3
"""Run VSEARCH SINTAX on ITS-standardized D1, D3 and D4 query FASTAs.

The standardized query FASTAs and ``standardization_audit.tsv`` are inputs to
this script. SINTAX is run with cutoff 0 so that confidence thresholds can be
applied downstream without reclassification.
"""
import argparse
import csv
import os
import subprocess
import time
from datetime import datetime
from pathlib import Path

def count_lines(path):
    with path.open('rb') as handle:
        return sum(1 for _ in handle)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-root', type=Path, required=True,
                    help='UMIvsRAW project root')
    ap.add_argument('--sintax-root', type=Path,
                    help='Default: <project-root>/REVIEW/SINTAX_equal_ITS_v1')
    ap.add_argument('--database', type=Path,
                    help='Default: <project-root>/REVIEW/general_unite_SHv10.fasta')
    ap.add_argument('--threads', type=int, default=63)
    args = ap.parse_args()

    project = args.project_root.resolve()
    root = (args.sintax_root or project/'REVIEW'/'SINTAX_equal_ITS_v1').resolve()
    database = (args.database or project/'REVIEW'/'general_unite_SHv10.fasta').resolve()
    audit = root/'standardization_audit.tsv'
    queries = root/'queries'
    results = root/'results'
    logs = root/'logs'
    results.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)

    if args.threads < 1:
        raise SystemExit('ERROR: --threads must be >=1')
    try:
        available = len(os.sched_getaffinity(0))
        if available < args.threads:
            raise SystemExit(f'ERROR: {available} CPUs available, {args.threads} requested')
    except AttributeError:
        available = None

    for path in (database, audit):
        if not path.is_file():
            raise FileNotFoundError(path)

    with audit.open(newline='') as handle:
        jobs = list(csv.DictReader(handle, delimiter='\t'))
    if len(jobs) != 291:
        raise ValueError(f'Expected 291 barcode-dataset rows, found {len(jobs)}')
    jobs.sort(key=lambda x: int(x['standardized_reads']))

    subprocess.run(['vsearch', '--version'], check=True)
    print(f'Database: {database}')
    print(f'Jobs: {len(jobs)}; threads/job: {args.threads}')

    completed = skipped = failed = 0
    started_all = time.monotonic()

    for i, job in enumerate(jobs, 1):
        bc, ds = job['barcode'], job['dataset']
        expected = int(job['standardized_reads'])
        query = queries/f'{bc}.{ds}.fasta'
        result = results/f'{bc}.{ds}.sintax.tsv'
        temp = results/f'{bc}.{ds}.sintax.tsv.tmp.{os.getpid()}'
        log = logs/f'{bc}.{ds}.vsearch.log'

        if result.exists() and count_lines(result) == expected:
            skipped += 1
            print(f'[{i}/{len(jobs)}] SKIP {bc} {ds}: {expected:,} validated rows', flush=True)
            continue
        if not query.is_file():
            failed += 1
            print(f'[{i}/{len(jobs)}] MISSING QUERY {query}', flush=True)
            continue
        temp.unlink(missing_ok=True)

        command = [
            'vsearch', '--sintax', str(query), '--db', str(database),
            '--tabbedout', str(temp), '--sintax_cutoff', '0.0',
            '--randseed', '42', '--strand', 'both', '--threads', str(args.threads)
        ]
        started = time.monotonic()
        with log.open('w') as lh:
            lh.write(f'Started: {datetime.now().isoformat()}\n')
            lh.write(f'Command: {" ".join(command)}\n\n')
            lh.flush()
            rc = subprocess.run(command, stdout=lh, stderr=subprocess.STDOUT).returncode
        if rc != 0 or not temp.exists():
            failed += 1
            temp.unlink(missing_ok=True)
            print(f'FAILED {bc} {ds}; see {log}', flush=True)
            continue
        observed = count_lines(temp)
        if observed != expected:
            failed += 1
            temp.unlink(missing_ok=True)
            print(f'FAILED {bc} {ds}: {observed:,}/{expected:,} rows', flush=True)
            continue
        os.replace(temp, result)
        completed += 1
        print(f'COMPLETE {bc} {ds}: {observed:,} rows in {(time.monotonic()-started)/60:.1f} min', flush=True)

    print(f'New: {completed}; existing: {skipped}; failed: {failed}; total: {completed+skipped}/{len(jobs)}')
    print(f'Elapsed: {(time.monotonic()-started_all)/3600:.2f} h')
    if failed:
        raise SystemExit(1)

if __name__ == '__main__':
    main()
