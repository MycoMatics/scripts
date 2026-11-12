#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from itertools import zip_longest

EXCLUDE = {"55MI", "66P1"}
SPECIAL_MI_D1 = {"35MI", "40MI", "41MI"}
ITS1F = "CTTGGTCATTTAGAGGAAGTAA"
ITS5_CORE = "AAGTCGTAACAAGGT"
ITS4 = "GCATATCAATAAGCGGAGG"
ITS4NGS_DEGENERATE = "GCATATHANTAAGSGSAGGCGA"
UUID_RE = re.compile(
    r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-"
    r"[0-9a-fA-F]{4}-[0-9a-fA-F]{12}"
)
UMI_RE = re.compile(r"TTT[ACG]{4}TT[ACG]{4}TT[ACG]{4}TT[ACG]{4}TTT")
RC_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def read_fasta(path: Path):
    header = None
    pieces = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(pieces)
                header = line[1:]
                pieces = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                pieces.append(line.strip())
    if header is not None:
        yield header, "".join(pieces)

def write_fasta(records, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with path.open("w", encoding="utf-8") as out:
        for identifier, sequence in records:
            if not identifier or not sequence:
                raise ValueError(f"Empty FASTA ID/sequence while writing {path}")
            out.write(f">{identifier}\n{sequence}\n")
            n += 1
    return n

def count_fasta(path: Path):
    return sum(1 for _ in read_fasta(path))

def parse_d1_id(header: str):
    return header.split(";", 1)[0].split()[0]

def parse_d3_id(header: str):
    """Validated historical D3 parser retained from prepare_sintax_queries.py."""
    accession = header.split(";", 1)[0].split()[0]
    match = re.search(r"(?:^|;)sample=([^;\s]+)", header)
    if match is None:
        raise ValueError(f"No sample= field in D3 header: {header}")
    sample = match.group(1)
    prefix = sample + "_"
    if not accession.startswith(prefix):
        raise ValueError(
            "D3 accession does not begin with its exact sample prefix:\n"
            f"accession={accession}\nsample={sample}\nheader={header}"
        )
    identifier = accession[len(prefix):]
    if identifier.endswith("rc"):
        identifier = identifier[:-2]
    if not identifier:
        raise ValueError(f"Empty D3 ID parsed from: {header}")
    return identifier

def parse_d4_id(header: str, barcode: str):
    """Recover a D4 read ID from the eNano/raw-all FASTA header.

    Historical output headers are UUID-based. The parser accepts the header
    forms encountered in related eNano files: direct UUIDs, sample-prefixed
    accessions, and accessions accompanied by ``sample=`` metadata. It fails
    rather than silently inventing an identifier when no UUID can be found.
    """
    token = header.split(";", 1)[0].split()[0]
    match = UUID_RE.search(token)
    if match:
        return match.group(0)
    match = UUID_RE.search(header)
    if match:
        return match.group(0)
    raise ValueError(f"No UUID-like D4 read ID found for {barcode}: {header}")

def normalize_d1_source(source: Path, barcode: str, output: Path):
    seen = set()
    records = []
    for header, seq in read_fasta(source):
        uid = parse_d1_id(header)
        if uid in seen:
            raise ValueError(f"Duplicate D1 UMI ID in {source}: {uid}")
        seen.add(uid)
        records.append((f"D1|{barcode}|{uid}", seq))
    write_fasta(records, output)
    return len(records)

def normalize_d3(source: Path, barcode: str, output: Path, valid_fastq_ids=None):
    """Stream D3 records to normalized query IDs without changing sequences."""
    output.parent.mkdir(parents=True, exist_ok=True)
    seen = set()
    n = 0
    with output.open("w", encoding="utf-8") as out:
        for header, seq in read_fasta(source):
            fid = parse_d3_id(header)
            if fid in seen:
                raise ValueError(f"Duplicate D3 read ID in {source}: {fid}")
            seen.add(fid)
            if valid_fastq_ids is not None and fid not in valid_fastq_ids:
                raise RuntimeError(f"{barcode}: D3 read ID absent from UMI mapping: {fid}")
            out.write(f">D3|{barcode}|{fid}\n{seq}\n")
            n += 1
    return n

def normalize_d4(source: Path, barcode: str, output: Path):
    """Stream D4 records to UUID-based normalized query IDs."""
    output.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with output.open("w", encoding="utf-8") as out:
        for header, seq in read_fasta(source):
            fid = parse_d4_id(header, barcode)
            out.write(f">D4|{barcode}|{fid}\n{seq}\n")
            n += 1
    return n

def concat_fastas(inputs, output: Path):
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as out:
        for path in inputs:
            with path.open("rb") as handle:
                shutil.copyfileobj(handle, out, length=1024 * 1024)
    return count_fasta(output)

def run_cutadapt(arguments, log: Path):
    command = ["cutadapt", *map(str, arguments)]
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            command,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"Cutadapt failed ({result.returncode}); see {log}")

def normalized_copy_d4(source: Path, barcode: str, output: Path):
    return normalize_d4(source, barcode, output)

def standardize_ordinary_mi_d1(source: Path, barcode: str, query: Path,
                               unresolved: Path, work: Path, logs: Path):
    primary = work / "D1_primary.fasta"
    missed1 = work / "D1_missed1.fasta"
    fallback = work / "D1_fallback.fasta"
    missed2 = work / "D1_missed2.fasta"

    run_cutadapt([
        "-g", ITS1F, "-e", "0.30", "-O", "18",
        "--untrimmed-output", missed1,
        "-o", primary, source,
    ], logs / f"{barcode}.D1.stage1.log")

    run_cutadapt([
        "-g", ITS5_CORE, "--action=retain", "-e", "0.25", "-O", "15",
        "--untrimmed-output", missed2,
        "-o", fallback, missed1,
    ], logs / f"{barcode}.D1.stage2.log")

    standardized_n = concat_fastas([primary, fallback], query)
    shutil.copyfile(missed2, unresolved)
    return standardized_n, count_fasta(unresolved)

def standardize_special_mi_d1(source: Path, barcode: str, query: Path,
                              unresolved: Path, work: Path, logs: Path):
    primary = work / "D1_primary.fasta"
    missed1 = work / "D1_missed1.fasta"
    stage2 = work / "D1_fallback.fasta"
    missed2 = work / "D1_missed2.fasta"

    run_cutadapt([
        "--revcomp", "--rename={header}", "--action=retain",
        "-g", ITS5_CORE, "-e", "0.10", "-O", "15",
        "--untrimmed-output", missed1,
        "-o", primary, source,
    ], logs / f"{barcode}.D1.stage1.log")

    run_cutadapt([
        "-a", ITS4, "-a", ITS4NGS_DEGENERATE,
        "-e", "0.15", "-O", "16",
        "--untrimmed-output", missed2,
        "-o", stage2, primary,
    ], logs / f"{barcode}.D1.stage2.log")

    shutil.copyfile(stage2, query)
    concat_fastas([missed1, missed2], unresolved)
    return count_fasta(query), count_fasta(unresolved)

def standardize_mi_d4(source: Path, barcode: str, query: Path,
                      unresolved: Path, work: Path, logs: Path):
    primary_raw = work / "D4_primary.raw.fasta"
    missed1_raw = work / "D4_missed1.raw.fasta"
    fallback_raw = work / "D4_fallback.raw.fasta"
    missed2_raw = work / "D4_missed2.raw.fasta"

    run_cutadapt([
        "-g", ITS1F, "-e", "0.30", "-O", "18",
        "--untrimmed-output", missed1_raw,
        "-o", primary_raw, source,
    ], logs / f"{barcode}.D4.stage1.log")

    run_cutadapt([
        "-g", ITS5_CORE, "--action=retain", "-e", "0.25", "-O", "15",
        "--untrimmed-output", missed2_raw,
        "-o", fallback_raw, missed1_raw,
    ], logs / f"{barcode}.D4.stage2.log")

    primary = work / "D4_primary.fasta"
    fallback = work / "D4_fallback.fasta"
    unresolved_norm = work / "D4_unresolved.fasta"
    normalize_d4(primary_raw, barcode, primary)
    normalize_d4(fallback_raw, barcode, fallback)
    normalize_d4(missed2_raw, barcode, unresolved_norm)
    standardized_n = concat_fastas([primary, fallback], query)
    shutil.copyfile(unresolved_norm, unresolved)
    return standardized_n, count_fasta(unresolved)

def load_mapping(project: Path, barcode: str):
    path = project / "SSUMI" / f"barcode{barcode}_ssumi" / f"{barcode}_mapping_umi.csv"
    umi_ids = set()
    fastq_ids = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        required = {"UMI_ID", "fastq_ID", "Barcode"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Unexpected mapping header in {path}: {reader.fieldnames}")
        for row in reader:
            fid = row["fastq_ID"]
            uid = row["UMI_ID"]
            if fid in fastq_ids:
                raise ValueError(f"Duplicate fastq_ID in {path}: {fid}")
            fastq_ids.add(fid)
            umi_ids.add(uid)
    return umi_ids, fastq_ids

def count_remaining_umi_motifs(path: Path):
    count = 0
    for _, seq in read_fasta(path):
        seq = seq.upper()
        rc = seq.translate(RC_TABLE)[::-1]
        if UMI_RE.search(seq) or UMI_RE.search(rc):
            count += 1
    return count

def compare_fasta(observed: Path, expected: Path):
    """Compare FASTAs in a streaming, order-sensitive manner."""
    for i, pair in enumerate(zip_longest(read_fasta(observed), read_fasta(expected)), 1):
        obs, exp = pair
        if obs is None:
            return False, f"observed FASTA ended before reference at record {i:,}"
        if exp is None:
            return False, f"reference FASTA ended before observed at record {i:,}"
        oh, osq = obs
        eh, esq = exp
        oid = oh.split()[0].split(";", 1)[0]
        eid = eh.split()[0].split(";", 1)[0]
        if oid != eid:
            return False, f"record {i:,} ID differs: {oid!r} vs {eid!r}"
        if osq.upper() != esq.upper():
            return False, f"record {i:,} sequence differs for {oid}"
    return True, "record order, IDs and nucleotide sequences identical"

def discover_sources(project: Path):
    d1_dir = project / "UMI_consensus_SHm" / "input_UMI_consensus"
    d3_dir = project / "UMI_rawinput_SHm" / "input_rawumi"
    d4_dir = project / "all_rawinput_SHm" / "input_rawall"

    d1 = {}
    for path in sorted(d1_dir.glob("barcode*.fasta")):
        match = re.fullmatch(r"barcode(.+)\.fasta", path.name)
        if match:
            d1[match.group(1)] = path

    d3 = {}
    for path in sorted(d3_dir.glob("barcode*_raw_umis.fasta")):
        match = re.fullmatch(r"barcode([0-9]+)_(P1|MI)_raw_umis\.fasta", path.name)
        if match:
            d3[match.group(1) + match.group(2)] = path

    d4 = {}
    for path in sorted(d4_dir.glob("barcode*.fasta")):
        match = re.fullmatch(r"barcode(.+)\.fasta", path.name)
        if match:
            d4[match.group(1)] = path

    barcodes = sorted((set(d1) & set(d3) & set(d4)) - EXCLUDE)
    if len(barcodes) != 97:
        raise RuntimeError(
            f"Expected 97 analyzable barcodes, found {len(barcodes)}. "
            f"D1={len(d1)}, D3={len(d3)}, D4={len(d4)}"
        )
    return barcodes, d1, d3, d4

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True,
                        help="UMIvsRAW project root")
    parser.add_argument("--output-root", type=Path,
                        help="Output analysis root; default REVIEW/SINTAX_equal_ITS_reconstructed")
    parser.add_argument("--reference-root", type=Path,
                        help="Archived final SINTAX_equal_ITS_v1 root for sequence-level validation")
    parser.add_argument("--force", action="store_true",
                        help="Allow replacement of files under --output-root")
    args = parser.parse_args()

    project = args.project_root.resolve()
    root = (args.output_root or project / "REVIEW" / "SINTAX_equal_ITS_reconstructed").resolve()
    reference = args.reference_root.resolve() if args.reference_root else None
    queries = root / "queries"
    unresolved_dir = root / "unresolved"
    logs = root / "logs"
    work_root = root / "work"

    if shutil.which("cutadapt") is None:
        raise SystemExit("ERROR: cutadapt is not on PATH. Historical production used Cutadapt 5.0.")

    if root.exists() and any(root.iterdir()) and not args.force:
        raise SystemExit(f"ERROR: {root} is not empty; use a new --output-root or --force")
    for directory in (queries, unresolved_dir, logs, work_root):
        directory.mkdir(parents=True, exist_ok=True)

    barcodes, d1_sources, d3_sources, d4_sources = discover_sources(project)
    audit_rows = []
    validation_failures = []

    print(f"Analyzable barcodes: {len(barcodes)}")
    print(f"Output root: {root}")
    if reference:
        print(f"Validation reference: {reference}")

    for index, bc in enumerate(barcodes, 1):
        print(f"[{index:02d}/{len(barcodes)}] {bc}", flush=True)
        map_umi, map_fastq = load_mapping(project, bc)

        with tempfile.TemporaryDirectory(prefix=f"{bc}_", dir=work_root) as tmp:
            work = Path(tmp)

            # D1 source normalization and mapping audit.
            d1_source = work / "D1_source.fasta"
            d1_input = normalize_d1_source(d1_sources[bc], bc, d1_source)
            d1_ids = [h.split("|", 2)[2] for h, _ in read_fasta(d1_source)]
            d1_not_mapping = sum(uid not in map_umi for uid in d1_ids)

            d1_query = queries / f"{bc}.D1.fasta"
            d1_unresolved = unresolved_dir / f"{bc}.D1.fasta"
            if bc.endswith("P1"):
                shutil.copyfile(d1_source, d1_query)
                d1_unresolved.write_text("", encoding="utf-8")
                d1_standardized, d1_unresolved_n = d1_input, 0
            elif bc in SPECIAL_MI_D1:
                d1_standardized, d1_unresolved_n = standardize_special_mi_d1(
                    d1_source, bc, d1_query, d1_unresolved, work, logs
                )
            elif bc.endswith("MI"):
                d1_standardized, d1_unresolved_n = standardize_ordinary_mi_d1(
                    d1_source, bc, d1_query, d1_unresolved, work, logs
                )
            else:
                raise ValueError(f"Unknown barcode class: {bc}")

            # D3: historically unchanged nucleotide sequences, normalized IDs.
            d3_query = queries / f"{bc}.D3.fasta"
            d3_unresolved = unresolved_dir / f"{bc}.D3.fasta"
            d3_input = normalize_d3(d3_sources[bc], bc, d3_query, map_fastq)
            d3_unresolved.write_text("", encoding="utf-8")

            # D4: P1 pass-through; MI trim to the common ITS1f-side interval.
            d4_query = queries / f"{bc}.D4.fasta"
            d4_unresolved = unresolved_dir / f"{bc}.D4.fasta"
            d4_input = count_fasta(d4_sources[bc])
            if bc.endswith("P1"):
                d4_standardized = normalized_copy_d4(d4_sources[bc], bc, d4_query)
                d4_unresolved.write_text("", encoding="utf-8")
                d4_unresolved_n = 0
            else:
                d4_standardized, d4_unresolved_n = standardize_mi_d4(
                    d4_sources[bc], bc, d4_query, d4_unresolved, work, logs
                )

        if d1_input != d1_standardized + d1_unresolved_n:
            raise RuntimeError(f"{bc} D1 count conservation failed")
        if d3_input != count_fasta(d3_query):
            raise RuntimeError(f"{bc} D3 count conservation failed")
        if d4_input != d4_standardized + d4_unresolved_n:
            raise RuntimeError(f"{bc} D4 count conservation failed")

        residual_umi = count_remaining_umi_motifs(d1_query)

        for ds, input_n, standardized_n, unresolved_n in [
            ("D1", d1_input, d1_standardized, d1_unresolved_n),
            ("D3", d3_input, d3_input, 0),
            ("D4", d4_input, d4_standardized, d4_unresolved_n),
        ]:
            audit_rows.append({
                "barcode": bc,
                "dataset": ds,
                "input_reads": input_n,
                "standardized_reads": standardized_n,
                "unresolved_reads": unresolved_n,
                "missing_uuid": 0,
                "remaining_umi": "",
                "recovery_percent": f"{100 * standardized_n / input_n:.4f}".rstrip("0").rstrip("."),
                "candidate_residual_umi_motif": residual_umi if ds == "D1" else 0,
                "candidate_d1_not_mapping": d1_not_mapping if ds == "D1" else 0,
            })

        if reference:
            for ds in ("D1", "D3", "D4"):
                observed = queries / f"{bc}.{ds}.fasta"
                expected = reference / "queries" / f"{bc}.{ds}.fasta"
                if not expected.is_file():
                    validation_failures.append(f"{bc} {ds}: missing reference {expected}")
                    continue
                ok, message = compare_fasta(observed, expected)
                print(f"    validate {ds}: {'PASS' if ok else 'FAIL'} - {message}")
                if not ok:
                    validation_failures.append(f"{bc} {ds}: {message}")

    audit = root / "standardization_audit.tsv"
    fields = [
        "barcode", "dataset", "input_reads", "standardized_reads",
        "unresolved_reads", "missing_uuid", "remaining_umi", "recovery_percent"
    ]
    with audit.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(audit_rows)

    diagnostics = root / "standardization_reconstruction_diagnostics.tsv"
    with diagnostics.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "barcode", "dataset", "candidate_residual_umi_motif",
                "candidate_d1_not_mapping"
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in audit_rows:
            if row["dataset"] == "D1":
                writer.writerow({k: row[k] for k in writer.fieldnames})

    if reference:
        reference_audit = reference / "standardization_audit.tsv"
        if not reference_audit.is_file():
            validation_failures.append(f"missing reference audit: {reference_audit}")
        else:
            with reference_audit.open(encoding="utf-8", newline="") as handle:
                expected_rows = {
                    (r["barcode"], r["dataset"]): r
                    for r in csv.DictReader(handle, delimiter="\t")
                }
            motif_matches = True
            mapping_matches = True
            for row in audit_rows:
                key = (row["barcode"], row["dataset"])
                expected = expected_rows.get(key)
                if expected is None:
                    validation_failures.append(f"audit row missing from reference: {key}")
                    continue
                for col in ("input_reads", "standardized_reads", "unresolved_reads", "missing_uuid"):
                    if int(row[col]) != int(expected[col]):
                        validation_failures.append(
                            f"audit {key} {col}: reconstructed={row[col]} reference={expected[col]}"
                        )
                if row["dataset"] == "D1":
                    exp_remaining = int(expected.get("remaining_umi", 0) or 0)
                    motif_matches &= row["candidate_residual_umi_motif"] == exp_remaining
                    mapping_matches &= row["candidate_d1_not_mapping"] == exp_remaining
            if motif_matches:
                print("Historical remaining_umi field matches residual UMI-motif counts for all D1 rows.")
            elif mapping_matches:
                print("Historical remaining_umi field matches D1 IDs absent from the UMI mapping for all D1 rows.")
            else:
                print("Historical remaining_umi field matches neither reconstructed diagnostic globally; it remains provenance-only.")

    print("\nReconstructed totals:")
    for ds in ("D1", "D3", "D4"):
        rows = [r for r in audit_rows if r["dataset"] == ds]
        print(
            f"  {ds}: input={sum(r['input_reads'] for r in rows):,}; "
            f"standardized={sum(r['standardized_reads'] for r in rows):,}; "
            f"unresolved={sum(r['unresolved_reads'] for r in rows):,}"
        )

    print(f"Audit: {audit}")
    if reference:
        if validation_failures:
            failure_file = root / "reference_validation_failures.txt"
            failure_file.write_text("\n".join(validation_failures) + "\n", encoding="utf-8")
            raise SystemExit(
                f"REFERENCE VALIDATION FAILED for {len(validation_failures)} dataset(s). "
                f"See {failure_file}"
            )
        print("REFERENCE VALIDATION PASSED: all 291 query FASTAs are identical by record order, ID and sequence.")

if __name__ == "__main__":
    main()
