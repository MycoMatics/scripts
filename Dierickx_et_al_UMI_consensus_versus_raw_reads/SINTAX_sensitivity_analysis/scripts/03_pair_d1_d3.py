#!/usr/bin/env python3

from pathlib import Path
import csv
import gzip
import os
import re
from collections import Counter

project_env = os.environ.get("PROJECT_ROOT")
if not project_env:
    raise SystemExit("ERROR: set PROJECT_ROOT to the UMIvsRAW project root")
ROOT = Path(project_env).resolve()
SINTAX = ROOT / "REVIEW" / "SINTAX_equal_ITS_v1"

RESULT_DIR = SINTAX / "results"

PAIR_DIR = SINTAX / "paired_by_barcode"

SUMMARY_DIR = SINTAX / "summary_by_barcode"

SSUMI = ROOT / "SSUMI"

MANIFEST = SINTAX / "standardization_audit.tsv"

PAIR_DIR.mkdir(
    parents=True,
    exist_ok=True
)

SUMMARY_DIR.mkdir(
    parents=True,
    exist_ok=True
)

CUTOFFS = [
    0.70,
    0.80,
    0.90
]

# SINTAX species/SH parser
# Example:

SH_RE = re.compile(
    r"(?:^|,)s:(SH[^,()]+)\(([0-9.]+)\)"
)

def parse_taxonomy(taxonomy):

    m = SH_RE.search(
        taxonomy
    )

    if m is None:
        return None, None

    return (
        m.group(1),
        float(m.group(2))
    )

# Count lines efficiently

def count_lines(path):

    n = 0

    with path.open("rb") as h:

        for _ in h:
            n += 1

    return n

# Parse one complete SINTAX result

def read_sintax(path, expected_bc):

    d1 = {}
    d3 = {}

    with path.open() as h:

        for line_number, line in enumerate(
            h,
            1
        ):

            line = line.rstrip(
                "\n\r"
            )

            if not line:
                continue

            fields = line.split(
                "\t"
            )

            if len(fields) < 3:

                raise RuntimeError(
                    f"{path}:{line_number}: "
                    f"expected >=3 columns"
                )

            query = fields[0]

            taxonomy = fields[1]

            strand = fields[2]

            parts = query.split(
                "|",
                2
            )

            if len(parts) != 3:

                raise RuntimeError(
                    f"{path}:{line_number}: "
                    f"unexpected query ID: {query}"
                )

            typ, bc, seq_id = parts

            if bc != expected_bc:

                raise RuntimeError(
                    f"{path}:{line_number}: "
                    f"expected barcode {expected_bc}, "
                    f"found {bc}"
                )

            sh, conf = parse_taxonomy(
                taxonomy
            )

            rec = (
                sh,
                conf,
                strand
            )

            if typ == "D1":

                if seq_id in d1:

                    raise RuntimeError(
                        f"Duplicate D1 ID: "
                        f"{bc} {seq_id}"
                    )

                d1[seq_id] = rec

            elif typ == "D3":

                if seq_id in d3:

                    raise RuntimeError(
                        f"Duplicate D3 ID: "
                        f"{bc} {seq_id}"
                    )

                d3[seq_id] = rec

            else:

                raise RuntimeError(
                    f"Unknown query type: {typ}"
                )

    return d1, d3

# UMI mapping

def read_mapping(bc):

    path = (
        SSUMI
        / f"barcode{bc}_ssumi"
        / f"{bc}_mapping_umi.csv"
    )

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    mapping = {}

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as h:

        reader = csv.DictReader(
            h,
            delimiter=";"
        )

        needed = {
            "UMI_ID",
            "fastq_ID",
            "Barcode"
        }

        if not needed.issubset(
            set(reader.fieldnames or [])
        ):

            raise RuntimeError(
                f"Unexpected mapping columns: "
                f"{reader.fieldnames}"
            )

        for row in reader:

            fid = row[
                "fastq_ID"
            ]

            uid = row[
                "UMI_ID"
            ]

            if fid in mapping:

                raise RuntimeError(
                    f"Duplicate fastq_ID in mapping: "
                    f"{bc} {fid}"
                )

            mapping[
                fid
            ] = uid

    return mapping

# Classification

def classify_pair(
    d1_sh,
    d1_conf,
    d3_sh,
    d3_conf,
    cutoff
):

    # Consensus not confidently classified:
    # outside denominator

    if (
        d1_sh is None
        or d1_conf is None
        or d1_conf < cutoff
    ):

        return "excluded_D1_low_confidence"

    # Consensus confident, raw read not confident

    if (
        d3_sh is None
        or d3_conf is None
        or d3_conf < cutoff
    ):

        return "unmatched"

    # Both confident and same SH

    if d3_sh == d1_sh:

        return "matched"

    # Both confident but different SH

    return "mismatched"

# Process one barcode

def process_barcode(
    bc,
    expected_d1,
    expected_d3
):

    result_file_d1 = RESULT_DIR / f"{bc}.D1.sintax.tsv"
    result_file_d3 = RESULT_DIR / f"{bc}.D3.sintax.tsv"

    pair_file = (
        PAIR_DIR
        / f"{bc}.paired.tsv.gz"
    )

    summary_file = (
        SUMMARY_DIR
        / f"{bc}.summary.tsv"
    )

    # Skip if already processed

    if (
        pair_file.exists()
        and pair_file.stat().st_size > 0
        and summary_file.exists()
        and summary_file.stat().st_size > 0
    ):

        return "already_processed"

    # Result must exist

    if not result_file_d1.exists() or not result_file_d3.exists():
        return "not_finished"

    # Only process completely written result files

    expected_total = (
        expected_d1
        + expected_d3
    )

    observed_total = count_lines(result_file_d1) + count_lines(result_file_d3)

    if observed_total != expected_total:

        return (
            "incomplete",
            observed_total,
            expected_total
        )

    # Parse SINTAX

    d1, accidentally_d3 = read_sintax(result_file_d1, bc)
    accidentally_d1, d3 = read_sintax(result_file_d3, bc)
    if accidentally_d1 or accidentally_d3:
        raise RuntimeError(f"{bc}: dataset mix-up between D1 and D3 result files")

    if len(d1) != expected_d1:

        raise RuntimeError(
            f"{bc}: expected "
            f"{expected_d1:,} D1, "
            f"found {len(d1):,}"
        )

    if len(d3) != expected_d3:

        raise RuntimeError(
            f"{bc}: expected "
            f"{expected_d3:,} D3, "
            f"found {len(d3):,}"
        )

    # Mapping

    mapping = read_mapping(
        bc
    )

    missing_mapping = [
        fid
        for fid in d3
        if fid not in mapping
    ]

    if missing_mapping:

        raise RuntimeError(
            f"{bc}: "
            f"{len(missing_mapping):,} "
            f"D3 IDs missing from mapping"
        )

    missing_d1 = [
        mapping[fid]
        for fid in d3
        if mapping[fid] not in d1
    ]

    # ITS standardization excludes a small number of D1 consensus sequences.
    # Keep their underlying D3 reads in the paired denominator. Missing D1
    # consensus assignments are therefore classified as low confidence.
    missing_d1_ids = set(missing_d1)
    missing_d1_read_count = sum(1 for fid in d3 if mapping[fid] in missing_d1_ids)
    if missing_d1_ids:
        print(f"{bc}: {len(missing_d1_ids):,} mapped D1 consensus IDs missing "
              "after standardization; associated D3 reads remain in denominator", flush=True)

    # Temporary output files
    # Atomic rename prevents half-written files from being

    pair_tmp = Path(
        str(pair_file)
        + ".tmp"
    )

    summary_tmp = Path(
        str(summary_file)
        + ".tmp"
    )

    counts = {
        cutoff: Counter()
        for cutoff in CUTOFFS
    }

    # Paired table

    with gzip.open(
        pair_tmp,
        "wt"
    ) as out:

        out.write(
            "\t".join(
                [
                    "Barcode",
                    "UMI_ID",
                    "fastq_ID",

                    "D1_SH",
                    "D1_confidence",
                    "D1_strand",

                    "D3_SH",
                    "D3_confidence",
                    "D3_strand",

                    "annotation_070",
                    "annotation_080",
                    "annotation_090",
                ]
            )
            + "\n"
        )

        for fid, d3_rec in d3.items():

            uid = mapping[
                fid
            ]

            d1_rec = d1.get(uid, (None, None, ""))

            d1_sh, d1_conf, d1_strand = (
                d1_rec
            )

            d3_sh, d3_conf, d3_strand = (
                d3_rec
            )

            annotations = []

            for cutoff in CUTOFFS:

                ann = classify_pair(
                    d1_sh,
                    d1_conf,
                    d3_sh,
                    d3_conf,
                    cutoff
                )

                annotations.append(
                    ann
                )

                counts[
                    cutoff
                ][ann] += 1

            row = [
                bc,
                uid,
                fid,

                d1_sh or "",

                (
                    ""
                    if d1_conf is None
                    else f"{d1_conf:.4f}"
                ),

                d1_strand,

                d3_sh or "",

                (
                    ""
                    if d3_conf is None
                    else f"{d3_conf:.4f}"
                ),

                d3_strand,

                *annotations
            ]

            out.write(
                "\t".join(
                    row
                )
                + "\n"
            )

    # Summary

    with summary_tmp.open(
        "w"
    ) as out:

        out.write(
            "\t".join(
                [
                    "Barcode",
                    "cutoff",
                    "total_D3",
                    "eligible",
                    "matched",
                    "unmatched",
                    "mismatched",
                    "excluded_D1_low_confidence",
                    "missing_D1_after_trimming",
                    "matched_pct_eligible",
                    "unmatched_pct_eligible",
                    "mismatched_pct_eligible",
                ]
            )
            + "\n"
        )

        for cutoff in CUTOFFS:

            x = counts[
                cutoff
            ]

            eligible = (
                x["matched"]
                + x["unmatched"]
                + x["mismatched"]
            )

            if eligible:

                matched_pct = (
                    100
                    * x["matched"]
                    / eligible
                )

                unmatched_pct = (
                    100
                    * x["unmatched"]
                    / eligible
                )

                mismatched_pct = (
                    100
                    * x["mismatched"]
                    / eligible
                )

            else:

                matched_pct = 0
                unmatched_pct = 0
                mismatched_pct = 0

            out.write(
                "\t".join(
                    map(
                        str,
                        [
                            bc,
                            f"{cutoff:.2f}",
                            expected_d3,
                            eligible,
                            x["matched"],
                            x["unmatched"],
                            x["mismatched"],
                            x[
                                "excluded_D1_low_confidence"
                            ],
                            missing_d1_read_count,
                            f"{matched_pct:.6f}",
                            f"{unmatched_pct:.6f}",
                            f"{mismatched_pct:.6f}",
                        ]
                    )
                )
                + "\n"
            )

    # Atomic completion

    os.replace(
        pair_tmp,
        pair_file
    )

    os.replace(
        summary_tmp,
        summary_file
    )

    return (
        "processed",
        expected_d1,
        expected_d3
    )

# Manifest reconstructed from standardization audit

manifest_data = {}
with MANIFEST.open(newline="") as h:
    reader = csv.DictReader(h, delimiter="\t")
    for row in reader:
        bc = row["barcode"]
        ds = row["dataset"]
        manifest_data.setdefault(bc, {})[ds] = int(row["standardized_reads"])

manifest = []
for bc, datasets in sorted(manifest_data.items()):
    if not {"D1", "D3"}.issubset(datasets):
        raise RuntimeError(f"Incomplete D1/D3 audit entry: {bc}")
    manifest.append((bc, datasets["D1"], datasets["D3"]))

# Process all currently complete barcodes

n_processed = 0

n_existing = 0

n_not_finished = 0

print(
    "Scanning completed SINTAX barcodes..."
)

print()

for bc, d1_n, d3_n in manifest:

    result = process_barcode(
        bc,
        d1_n,
        d3_n
    )

    if result == "already_processed":

        n_existing += 1

        continue

    if result == "not_finished":

        n_not_finished += 1

        continue

    if isinstance(
        result,
        tuple
    ):

        if result[0] == "incomplete":

            n_not_finished += 1

            print(
                f"{bc}: currently incomplete "
                f"{result[1]:,}/{result[2]:,}, skipped"
            )

            continue

        if result[0] == "processed":

            n_processed += 1

            print(
                f"{bc}: processed "
                f"D1={result[1]:,}, "
                f"D3={result[2]:,}"
            )

print()
print(
    "============================================================"
)

print(
    "INCREMENTAL PROCESSING SUMMARY"
)

print(
    "============================================================"
)

print(
    f"Newly processed       = {n_processed}"
)

print(
    f"Previously processed  = {n_existing}"
)

print(
    f"Not yet complete      = {n_not_finished}"
)

print(
    f"Total manifest        = {len(manifest)}"
)
