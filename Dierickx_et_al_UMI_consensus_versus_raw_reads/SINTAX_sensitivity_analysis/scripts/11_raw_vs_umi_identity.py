#!/usr/bin/env python3

"""Estimate divergence and Q scores for mapped D3 read–D1 consensus pairs.

Both consensus orientations are tested. edlib HW alignment excludes unaligned
terminal sequence; divergence is calculated over the homologous interval.
Read-level Q = -10 log10(mean(10 ** (-base_Q/10))), not mean(base_Q).
Outputs support the common-denominator quality/divergence Figure B3.
"""

import argparse
import csv
import math
import re
import sys

from concurrent.futures import (
    ProcessPoolExecutor,
    as_completed,
)

from pathlib import Path

import edlib

# Sequence helpers

IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYMKWSBDHVNacgtrymkwsbdhvn",
    "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn",
)

def reverse_complement(seq):
    """Reverse-complement DNA, including common IUPAC ambiguity codes."""

    return seq.translate(
        IUPAC_COMPLEMENT
    )[::-1]

# FASTA reader

def read_fasta(path):
    """
    Stream FASTA records.

    Yields
    ------
    header : str
        Header without leading ">".

    sequence : str
        Upper-case sequence.
    """

    header = None
    seq = []

    with open(
        path,
        "r",
        encoding="utf-8"
    ) as handle:

        for line in handle:

            line = line.rstrip(
                "\r\n"
            )

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:

                    yield (
                        header,
                        "".join(seq).upper()
                    )

                header = line[1:]
                seq = []

            else:

                seq.append(
                    line.strip()
                )

        if header is not None:

            yield (
                header,
                "".join(seq).upper()
            )

# TD3 processed-read identifier parser

def processed_raw_info(header):
    """
    Recover the exact original fastq_ID from a processed TD3 FASTA header.

    This uses the sample= field as the exact prefix and therefore
    preserves suffixes such as _1 and _2.

    Examples
    --------
    >barcode32_P1_raw_umis_<ID>rc;sample=barcode32_P1_raw_umis

    >barcode35MIumi_<ID>;sample=barcode35MIumi

    Returns
    -------
    fastq_id : str

    is_rc : bool
        Whether the processed TD3 accession ended in "rc".

    Notes
    -----
    The rc flag is used only to recover the correct quality interval
    from the original FASTQ.

    It is NOT assumed to describe the orientation of the sequence
    relative to the UMI consensus.
    """

    parts = header.split(";")

    accession = (
        parts[0]
        .strip()
        .split()[0]
    )

    sample = None

    for field in parts[1:]:

        field = field.strip()

        if field.startswith(
            "sample="
        ):

            sample = field[
                len("sample="):
            ].strip()

            break

    if sample is None:

        raise ValueError(
            "No sample= field in processed "
            f"FASTA header:\n{header}"
        )

    expected_prefix = (
        sample + "_"
    )

    if not accession.startswith(
        expected_prefix
    ):

        raise ValueError(
            "Processed accession does not start "
            "with its sample prefix:\n"
            f"header={header}\n"
            f"sample={sample}\n"
            f"expected_prefix={expected_prefix}"
        )

    fastq_id = accession[
        len(expected_prefix):
    ]

    is_rc = fastq_id.endswith(
        "rc"
    )

    if is_rc:

        fastq_id = fastq_id[
            :-2
        ]

    if not fastq_id:

        raise ValueError(
            f"Empty fastq_ID parsed from: {header}"
        )

    return (
        fastq_id,
        is_rc
    )

# Mapping files

def detect_delimiter(path):
    """Detect comma versus semicolon for a small CSV header."""

    with open(
        path,
        "r",
        encoding="utf-8-sig"
    ) as handle:

        header = handle.readline()

    if header.count(";") >= header.count(","):
        return ";"

    return ","

def read_umi_mapping(path):
    """
    Read fastq_ID -> UMI_ID mapping.

    Duplicate identifiers mapping to different UMI bins are considered
    fatal because the biological pairing would be ambiguous.
    """

    delimiter = detect_delimiter(
        path
    )

    mapping = {}

    with open(
        path,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter=delimiter
        )

        required = {
            "fastq_ID",
            "UMI_ID",
        }

        missing = (
            required
            - set(
                reader.fieldnames or []
            )
        )

        if missing:

            raise ValueError(
                f"{path} is missing columns: "
                f"{sorted(missing)}"
            )

        for row in reader:

            fastq_id = (
                row["fastq_ID"]
                .strip()
            )

            umi_id = (
                row["UMI_ID"]
                .strip()
            )

            if not fastq_id:
                continue

            previous = mapping.get(
                fastq_id
            )

            if (
                previous is not None
                and previous != umi_id
            ):

                raise ValueError(
                    "fastq_ID maps to multiple "
                    "UMI_ID values:\n"
                    f"{fastq_id}\n"
                    f"{previous}\n"
                    f"{umi_id}"
                )

            mapping[
                fastq_id
            ] = umi_id

    return mapping

# UMI consensus FASTA

def read_consensus_fasta(path):
    """
    Read UMI_ID -> consensus sequence.

    Consensus header example:

        >umi26771bins;ubs=5 LN:i:952 ...
    """

    result = {}

    for header, seq in read_fasta(
        path
    ):

        umi_id = (
            header
            .split(";", 1)[0]
            .split()[0]
        )

        if umi_id in result:

            raise ValueError(
                "Duplicate UMI consensus ID "
                f"in {path}: {umi_id}"
            )

        result[
            umi_id
        ] = seq

    return result

# Processed TD3 reads

def read_processed_reads(
    path,
    limit=None
):
    """
    Read processed TD3 sequences.

    Returns
    -------
    dict

        fastq_ID ->
            {
                "sequence": str,
                "td3_rc": bool
            }
    """

    result = {}

    for header, seq in read_fasta(
        path
    ):

        fastq_id, is_rc = (
            processed_raw_info(
                header
            )
        )

        if fastq_id in result:

            raise ValueError(
                "Duplicate processed fastq_ID "
                f"in {path}: {fastq_id}"
            )

        result[
            fastq_id
        ] = {
            "sequence": seq,
            "td3_rc": is_rc,
        }

        if (
            limit is not None
            and len(result) >= limit
        ):
            break

    return result

# Nanopore Q score

def nanopore_read_qscore(
    quality_string
):
    """
    Calculate the read-level Nanopore Q score from base qualities.

    This reproduces the logic used in R/read_Q_nanopore.R:

        p_i = 10^(-Q_i / 10)

        Q_read = -10 * log10(mean(p_i))
    """

    if not quality_string:
        return None

    sum_error_probability = 0.0

    for character in quality_string:

        q = (
            ord(character)
            - 33
        )

        sum_error_probability += (
            10.0 ** (
                -q / 10.0
            )
        )

    mean_error_probability = (
        sum_error_probability
        / len(quality_string)
    )

    if mean_error_probability <= 0:
        return None

    return (
        -10.0
        * math.log10(
            mean_error_probability
        )
    )

def recover_processed_qscores(
    fastq_path,
    processed_reads
):
    """
    Recover the exact FASTQ quality interval corresponding to each
    processed TD3 sequence.

    If the TD3 accession was marked rc, the original FASTQ sequence is
    reverse-complemented and the quality string is reversed before the
    exact processed sequence is located.

    Returns
    -------
    qscores : dict
        fastq_ID -> qscore or None

    missing_fastq : int
        Processed IDs never encountered in the original FASTQ.

    not_exact : int
        IDs encountered but for which the processed sequence could not
        be recovered as an exact substring.
    """

    pending = set(
        processed_reads.keys()
    )

    qscores = {}

    not_exact = 0

    with open(
        fastq_path,
        "r",
        encoding="utf-8"
    ) as handle:

        while True:

            header = handle.readline()

            if not header:
                break

            sequence = (
                handle
                .readline()
                .rstrip("\r\n")
                .upper()
            )

            plus = handle.readline()

            quality = (
                handle
                .readline()
                .rstrip("\r\n")
            )

            if (
                not sequence
                and not quality
            ):
                continue

            if not header.startswith(
                "@"
            ):

                raise ValueError(
                    "Malformed FASTQ header in "
                    f"{fastq_path}: {header!r}"
                )

            full_id = (
                header[1:]
                .strip()
            )

            first_token = (
                full_id
                .split()[0]
            )

            fastq_id = None

            if full_id in pending:

                fastq_id = full_id

            elif first_token in pending:

                fastq_id = first_token

            if fastq_id is None:
                continue

            info = processed_reads[
                fastq_id
            ]

            processed_seq = info[
                "sequence"
            ]

            if info["td3_rc"]:

                oriented_sequence = (
                    reverse_complement(
                        sequence
                    )
                )

                oriented_quality = (
                    quality[::-1]
                )

            else:

                oriented_sequence = (
                    sequence
                )

                oriented_quality = (
                    quality
                )

            start = (
                oriented_sequence
                .find(
                    processed_seq
                )
            )

            if start < 0:

                qscores[
                    fastq_id
                ] = None

                not_exact += 1

            else:

                end = (
                    start
                    + len(
                        processed_seq
                    )
                )

                quality_slice = (
                    oriented_quality[
                        start:end
                    ]
                )

                if (
                    len(quality_slice)
                    != len(processed_seq)
                ):

                    qscores[
                        fastq_id
                    ] = None

                    not_exact += 1

                else:

                    qscores[
                        fastq_id
                    ] = (
                        nanopore_read_qscore(
                            quality_slice
                        )
                    )

            pending.remove(
                fastq_id
            )

            if not pending:
                break

    missing_fastq = len(
        pending
    )

    for fastq_id in pending:

        qscores[
            fastq_id
        ] = None

    return (
        qscores,
        missing_fastq,
        not_exact,
    )

# edlib CIGAR

CIGAR_TOKEN = re.compile(
    r"(\d+)([=XIDM])"
)

def cigar_counts(cigar):
    """
    Parse edlib CIGAR.

    Expected operations are:
        =   match
        X   substitution
        I   insertion in query relative to target
        D   deletion in query relative to target

    An M operation is rejected because it does not distinguish matches
    from substitutions.
    """

    counts = {
        "=": 0,
        "X": 0,
        "I": 0,
        "D": 0,
    }

    if cigar is None:

        raise ValueError(
            "edlib returned no CIGAR"
        )

    tokens = CIGAR_TOKEN.findall(
        cigar
    )

    if not tokens:

        raise ValueError(
            f"Could not parse CIGAR: {cigar}"
        )

    reconstructed = "".join(
        f"{n}{op}"
        for n, op in tokens
    )

    if reconstructed != cigar:

        raise ValueError(
            f"Unexpected CIGAR format: {cigar}"
        )

    for number, operation in tokens:

        if operation == "M":

            raise ValueError(
                "edlib returned CIGAR operation M. "
                "Extended =/X operations are required."
            )

        counts[
            operation
        ] += int(
            number
        )

    return counts

# Semi-global shared-interval alignment

def align_shared_interval(
    raw_seq,
    oriented_consensus
):
    """
    Align the shorter sequence completely within the longer sequence.

    Terminal sequence present only in the longer sequence is not counted
    as an error.

    All insertion/deletion counts returned by this function are
    normalized to the biological interpretation:

        insertion_in_raw
        deletion_in_raw

    regardless of which sequence was the edlib query.
    """

    raw_length = len(
        raw_seq
    )

    consensus_length = len(
        oriented_consensus
    )

    # Choose the shorter sequence as query.
    # HW aligns the complete query to a substring of the target.
    # Terminal target sequence is therefore free and does not contribute
    # to edit distance.

    if raw_length <= consensus_length:

        query = raw_seq
        target = oriented_consensus

        query_type = "raw"
        target_type = "consensus"

    else:

        query = oriented_consensus
        target = raw_seq

        query_type = "consensus"
        target_type = "raw"

    result = edlib.align(
        query,
        target,
        mode="HW",
        task="path",
    )

    if result[
        "editDistance"
    ] < 0:

        raise RuntimeError(
            "edlib failed to produce an alignment"
        )

    if not result.get(
        "locations"
    ):

        raise RuntimeError(
            "edlib returned no alignment location"
        )

    counts = cigar_counts(
        result["cigar"]
    )

    matches = counts["="]
    substitutions = counts["X"]

    # Normalize I/D to raw relative to consensus.

    if query_type == "raw":

        insertions_in_raw = (
            counts["I"]
        )

        deletions_in_raw = (
            counts["D"]
        )

    else:

        # Query = consensus, target = raw.
        # I in query means sequence is present in consensus but absent
        # D in query means sequence is present in raw but absent from
        # consensus, therefore insertion in raw.

        insertions_in_raw = (
            counts["D"]
        )

        deletions_in_raw = (
            counts["I"]
        )

    edit_distance = (
        substitutions
        + insertions_in_raw
        + deletions_in_raw
    )

    if (
        edit_distance
        != result[
            "editDistance"
        ]
    ):

        raise RuntimeError(
            "CIGAR-derived edit distance "
            "does not match edlib editDistance"
        )

    alignment_columns = (
        matches
        + substitutions
        + insertions_in_raw
        + deletions_in_raw
    )

    if alignment_columns <= 0:

        raise RuntimeError(
            "Alignment contained zero columns"
        )

    identity_pct = (
        100.0
        * matches
        / alignment_columns
    )

    error_pct = (
        100.0
        * edit_distance
        / alignment_columns
    )

    # Terminal target sequence not participating in the alignment.

    start, end = (
        result[
            "locations"
        ][0]
    )

    target_left_unaligned = (
        start
    )

    target_right_unaligned = (
        len(target)
        - 1
        - end
    )

    if target_type == "consensus":

        consensus_left_unaligned = (
            target_left_unaligned
        )

        consensus_right_unaligned = (
            target_right_unaligned
        )

        raw_left_unaligned = 0
        raw_right_unaligned = 0

        consensus_aligned_span = (
            end
            - start
            + 1
        )

        raw_aligned_span = (
            raw_length
        )

    else:

        raw_left_unaligned = (
            target_left_unaligned
        )

        raw_right_unaligned = (
            target_right_unaligned
        )

        consensus_left_unaligned = 0
        consensus_right_unaligned = 0

        raw_aligned_span = (
            end
            - start
            + 1
        )

        consensus_aligned_span = (
            consensus_length
        )

    # Coverage of one sequence by the homologous interval of the other.
    # These are sequence-extent metrics, not sequence-error metrics.

    consensus_covered_by_raw_pct = (
        100.0
        * consensus_aligned_span
        / consensus_length
    )

    raw_covered_by_consensus_pct = (
        100.0
        * raw_aligned_span
        / raw_length
    )

    return {

        "alignment_query":
            query_type,

        "matches_bp":
            matches,

        "substitutions_bp":
            substitutions,

        "insertions_in_raw_bp":
            insertions_in_raw,

        "deletions_in_raw_bp":
            deletions_in_raw,

        "edit_distance_bp":
            edit_distance,

        "alignment_columns_bp":
            alignment_columns,

        "identity_pct":
            identity_pct,

        "error_pct":
            error_pct,

        "consensus_left_unaligned_bp":
            consensus_left_unaligned,

        "consensus_right_unaligned_bp":
            consensus_right_unaligned,

        "raw_left_unaligned_bp":
            raw_left_unaligned,

        "raw_right_unaligned_bp":
            raw_right_unaligned,

        "consensus_covered_by_raw_pct":
            consensus_covered_by_raw_pct,

        "raw_covered_by_consensus_pct":
            raw_covered_by_consensus_pct,
    }

# Best orientation of the PAIRED consensus

def best_paired_consensus_alignment(
    raw_seq,
    consensus_seq
):
    """
    Compare the raw read with both orientations of its paired consensus.

    This does NOT search among different UMI consensuses.

    Selection:
        1. lower edit distance
        2. lower overlap error if tied
        3. forward orientation if still tied
    """

    forward = align_shared_interval(
        raw_seq,
        consensus_seq
    )

    reverse_consensus = (
        reverse_complement(
            consensus_seq
        )
    )

    reverse = align_shared_interval(
        raw_seq,
        reverse_consensus
    )

    if (
        reverse[
            "edit_distance_bp"
        ]
        <
        forward[
            "edit_distance_bp"
        ]
    ):

        best = reverse
        orientation = "reverse"

    elif (
        forward[
            "edit_distance_bp"
        ]
        <
        reverse[
            "edit_distance_bp"
        ]
    ):

        best = forward
        orientation = "forward"

    elif (
        reverse[
            "error_pct"
        ]
        <
        forward[
            "error_pct"
        ]
    ):

        best = reverse
        orientation = "reverse"

    else:

        best = forward
        orientation = "forward"

    result = dict(
        best
    )

    result[
        "consensus_orientation"
    ] = orientation

    return result

# Paths

def barcode_paths(
    root,
    barcode
):
    """
    Construct all input paths for one barcode.
    """

    match = re.match(
        r"^(\d+)(P1|P2|P3|MI)$",
        barcode
    )

    if not match:

        raise ValueError(
            f"Unexpected barcode format: {barcode}"
        )

    digits = match.group(1)
    suffix = match.group(2)

    processed_fasta = (
        root
        / "UMI_rawinput_SHm"
        / "input_rawumi"
        / (
            f"barcode{digits}_"
            f"{suffix}_raw_umis.fasta"
        )
    )

    mapping_csv = (
        root
        / "SSUMI"
        / f"barcode{barcode}_ssumi"
        / f"{barcode}_mapping_umi.csv"
    )

    consensus_fasta = (
        root
        / "UMI_consensus_SHm"
        / "input_UMI_consensus"
        / f"barcode{barcode}.fasta"
    )

    original_fastq = (
        root
        / "SSUMI"
        / f"barcode{barcode}_ssumi"
        / f"raw_UMIs_bc{barcode}.fastq"
    )

    return {
        "processed_fasta":
            processed_fasta,

        "mapping_csv":
            mapping_csv,

        "consensus_fasta":
            consensus_fasta,

        "original_fastq":
            original_fastq,
    }

# Output columns

OUTPUT_COLUMNS = [

    "Barcode",
    "fastq_ID",
    "UMI_ID",

    "td3_rc",

    "raw_length_bp",
    "consensus_length_bp",
    "length_difference_bp",

    "processed_avg_qscore",

    "consensus_orientation",
    "alignment_query",

    "overlap_matches_bp",
    "overlap_substitutions_bp",
    "overlap_insertions_in_raw_bp",
    "overlap_deletions_in_raw_bp",

    "overlap_edit_distance_bp",
    "overlap_alignment_columns_bp",

    "overlap_identity_pct",
    "overlap_error_pct",

    "consensus_left_unaligned_bp",
    "consensus_right_unaligned_bp",
    "consensus_total_unaligned_bp",

    "raw_left_unaligned_bp",
    "raw_right_unaligned_bp",
    "raw_total_unaligned_bp",

    "consensus_covered_by_raw_pct",
    "raw_covered_by_consensus_pct",
]

# Process one barcode

def process_barcode(
    root_string,
    output_dir_string,
    barcode,
    limit_per_barcode=None
):
    """
    Process one barcode and write one TSV.
    """

    root = Path(
        root_string
    )

    output_dir = Path(
        output_dir_string
    )

    paths = barcode_paths(
        root,
        barcode
    )

    for label, path in paths.items():

        if not path.exists():

            raise FileNotFoundError(
                f"{barcode}: missing {label}: "
                f"{path}"
            )

    mapping = read_umi_mapping(
        paths[
            "mapping_csv"
        ]
    )

    consensuses = (
        read_consensus_fasta(
            paths[
                "consensus_fasta"
            ]
        )
    )

    processed = (
        read_processed_reads(
            paths[
                "processed_fasta"
            ],
            limit=limit_per_barcode,
        )
    )

    (
        qscores,
        qscore_missing_fastq,
        qscore_not_exact,
    ) = recover_processed_qscores(
        paths[
            "original_fastq"
        ],
        processed,
    )

    outfile = (
        output_dir
        / (
            "raw_vs_consensus_"
            f"{barcode}.tsv"
        )
    )

    tmpfile = outfile.with_suffix(
        ".tsv.tmp"
    )

    processed_raw = len(
        processed
    )

    aligned = 0

    missing_mapping = 0
    missing_consensus = 0

    forward_orientation_best = 0
    reverse_orientation_best = 0

    raw_as_query = 0
    consensus_as_query = 0

    with open(
        tmpfile,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for fastq_id, info in processed.items():

            umi_id = mapping.get(
                fastq_id
            )

            if umi_id is None:

                missing_mapping += 1
                continue

            consensus_seq = (
                consensuses.get(
                    umi_id
                )
            )

            if consensus_seq is None:

                missing_consensus += 1
                continue

            raw_seq = info[
                "sequence"
            ]

            alignment = (
                best_paired_consensus_alignment(
                    raw_seq,
                    consensus_seq,
                )
            )

            orientation = alignment[
                "consensus_orientation"
            ]

            if orientation == "forward":

                forward_orientation_best += 1

            else:

                reverse_orientation_best += 1

            if (
                alignment[
                    "alignment_query"
                ]
                == "raw"
            ):

                raw_as_query += 1

            else:

                consensus_as_query += 1

            qscore = qscores.get(
                fastq_id
            )

            consensus_total_unaligned = (
                alignment[
                    "consensus_left_unaligned_bp"
                ]
                +
                alignment[
                    "consensus_right_unaligned_bp"
                ]
            )

            raw_total_unaligned = (
                alignment[
                    "raw_left_unaligned_bp"
                ]
                +
                alignment[
                    "raw_right_unaligned_bp"
                ]
            )

            row = {

                "Barcode":
                    barcode,

                "fastq_ID":
                    fastq_id,

                "UMI_ID":
                    umi_id,

                "td3_rc":
                    int(
                        info["td3_rc"]
                    ),

                "raw_length_bp":
                    len(raw_seq),

                "consensus_length_bp":
                    len(consensus_seq),

                "length_difference_bp":
                    (
                        len(raw_seq)
                        - len(consensus_seq)
                    ),

                "processed_avg_qscore":
                    (
                        "NA"
                        if qscore is None
                        else f"{qscore:.6f}"
                    ),

                "consensus_orientation":
                    orientation,

                "alignment_query":
                    alignment[
                        "alignment_query"
                    ],

                "overlap_matches_bp":
                    alignment[
                        "matches_bp"
                    ],

                "overlap_substitutions_bp":
                    alignment[
                        "substitutions_bp"
                    ],

                "overlap_insertions_in_raw_bp":
                    alignment[
                        "insertions_in_raw_bp"
                    ],

                "overlap_deletions_in_raw_bp":
                    alignment[
                        "deletions_in_raw_bp"
                    ],

                "overlap_edit_distance_bp":
                    alignment[
                        "edit_distance_bp"
                    ],

                "overlap_alignment_columns_bp":
                    alignment[
                        "alignment_columns_bp"
                    ],

                "overlap_identity_pct":
                    (
                        f"{alignment['identity_pct']:.6f}"
                    ),

                "overlap_error_pct":
                    (
                        f"{alignment['error_pct']:.6f}"
                    ),

                "consensus_left_unaligned_bp":
                    alignment[
                        "consensus_left_unaligned_bp"
                    ],

                "consensus_right_unaligned_bp":
                    alignment[
                        "consensus_right_unaligned_bp"
                    ],

                "consensus_total_unaligned_bp":
                    consensus_total_unaligned,

                "raw_left_unaligned_bp":
                    alignment[
                        "raw_left_unaligned_bp"
                    ],

                "raw_right_unaligned_bp":
                    alignment[
                        "raw_right_unaligned_bp"
                    ],

                "raw_total_unaligned_bp":
                    raw_total_unaligned,

                "consensus_covered_by_raw_pct":
                    (
                        f"{alignment['consensus_covered_by_raw_pct']:.6f}"
                    ),

                "raw_covered_by_consensus_pct":
                    (
                        f"{alignment['raw_covered_by_consensus_pct']:.6f}"
                    ),
            }

            writer.writerow(
                row
            )

            aligned += 1

    tmpfile.replace(
        outfile
    )

    return {

        "barcode":
            barcode,

        "processed_raw":
            processed_raw,

        "aligned":
            aligned,

        "missing_mapping":
            missing_mapping,

        "missing_consensus":
            missing_consensus,

        "qscore_missing_fastq":
            qscore_missing_fastq,

        "qscore_not_exact":
            qscore_not_exact,

        "forward_orientation_best":
            forward_orientation_best,

        "reverse_orientation_best":
            reverse_orientation_best,

        "raw_as_query":
            raw_as_query,

        "consensus_as_query":
            consensus_as_query,

        "outfile":
            str(outfile),
    }

# Barcode discovery

def discover_barcodes(
    root
):
    """
    Discover barcodes from processed TD3 FASTA filenames.
    """

    directory = (
        root
        / "UMI_rawinput_SHm"
        / "input_rawumi"
    )

    pattern = re.compile(
        r"^barcode(\d+)_(P1|P2|P3|MI)"
        r"_raw_umis\.fasta$"
    )

    barcodes = []

    for path in directory.glob(
        "barcode*_raw_umis.fasta"
    ):

        match = pattern.match(
            path.name
        )

        if not match:
            continue

        barcode = (
            match.group(1)
            + match.group(2)
        )

        barcodes.append(
            barcode
        )

    def sort_key(barcode):

        match = re.match(
            r"^(\d+)(.+)$",
            barcode
        )

        return (
            int(match.group(1)),
            match.group(2),
        )

    return sorted(
        set(barcodes),
        key=sort_key,
    )

# Combine per-barcode output

def combine_outputs(
    summaries,
    combined_path
):
    """
    Concatenate per-barcode TSVs without loading them into memory.
    """

    first = True

    with open(
        combined_path,
        "w",
        encoding="utf-8"
    ) as output:

        for summary in summaries:

            path = Path(
                summary[
                    "outfile"
                ]
            )

            with open(
                path,
                "r",
                encoding="utf-8"
            ) as handle:

                header = handle.readline()

                if first:

                    output.write(
                        header
                    )

                    first = False

                for line in handle:

                    output.write(
                        line
                    )

# Write processing summary

SUMMARY_COLUMNS = [

    "barcode",
    "processed_raw",
    "aligned",

    "missing_mapping",
    "missing_consensus",

    "qscore_missing_fastq",
    "qscore_not_exact",

    "forward_orientation_best",
    "reverse_orientation_best",

    "raw_as_query",
    "consensus_as_query",

    "outfile",
]

def write_processing_summary(
    summaries,
    path
):

    with open(
        path,
        "w",
        encoding="utf-8",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for summary in summaries:

            writer.writerow(
                summary
            )

# CLI

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Compare processed TD3 reads with "
            "their paired UMI-derived consensus "
            "over the shared homologous interval."
        )
    )

    parser.add_argument(
        "--root",
        default="..",
        help=(
            "UMIvsRAW root directory. "
            "Default: .."
        ),
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=20,
        help=(
            "Number of barcode workers. "
            "Default: 20"
        ),
    )

    parser.add_argument(
        "--output-dir",
        default=(
            "REVIEW/"
            "raw_vs_umi_identity_harmonized"
        ),
        help=(
            "Output directory relative to ROOT "
            "unless an absolute path is supplied."
        ),
    )

    parser.add_argument(
        "--barcodes",
        nargs="*",
        default=None,
        help=(
            "Optional barcode subset, for example "
            "--barcodes 01P1 36MI 35MI"
        ),
    )

    parser.add_argument(
        "--limit-per-barcode",
        type=int,
        default=None,
        help=(
            "Optional maximum number of processed "
            "reads per barcode for testing."
        ),
    )

    return parser.parse_args()

# Main

def main():

    args = parse_args()

    root = Path(
        args.root
    ).resolve()

    output_dir_arg = Path(
        args.output_dir
    )

    if output_dir_arg.is_absolute():

        output_dir = (
            output_dir_arg
        )

    else:

        output_dir = (
            root
            / output_dir_arg
        )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if args.barcodes:

        barcodes = list(
            args.barcodes
        )

    else:

        barcodes = (
            discover_barcodes(
                root
            )
        )

    if not barcodes:

        raise RuntimeError(
            "No TD3 barcodes found."
        )

    print(
        f"ROOT: {root}"
    )

    print(
        f"Output: {output_dir}"
    )

    print(
        f"Barcodes: {len(barcodes)}"
    )

    print(
        f"Workers: {args.workers}"
    )

    if (
        args.limit_per_barcode
        is not None
    ):

        print(
            "TEST LIMIT PER BARCODE: "
            f"{args.limit_per_barcode:,}"
        )

    summaries = []
    failures = []

    with ProcessPoolExecutor(
        max_workers=args.workers
    ) as executor:

        future_to_barcode = {}

        for barcode in barcodes:

            future = executor.submit(
                process_barcode,
                str(root),
                str(output_dir),
                barcode,
                args.limit_per_barcode,
            )

            future_to_barcode[
                future
            ] = barcode

        for future in as_completed(
            future_to_barcode
        ):

            barcode = (
                future_to_barcode[
                    future
                ]
            )

            try:

                summary = (
                    future.result()
                )

                summaries.append(
                    summary
                )

                print(
                    f"[DONE] {barcode}: "
                    f"{summary['aligned']:,} aligned"
                )

            except Exception as exc:

                failures.append(
                    (
                        barcode,
                        repr(exc),
                    )
                )

                print(
                    f"[FAILED] {barcode}: {exc}",
                    file=sys.stderr,
                )

    if failures:

        print(
            "\nFAILED BARCODES:",
            file=sys.stderr,
        )

        for barcode, error in failures:

            print(
                f"{barcode}\t{error}",
                file=sys.stderr,
            )

        raise RuntimeError(
            f"{len(failures)} barcode(s) failed. "
            "Combined output was not created."
        )

    # Restore deterministic barcode order.

    order = {
        barcode: i
        for i, barcode
        in enumerate(barcodes)
    }

    summaries.sort(
        key=lambda x: order[
            x["barcode"]
        ]
    )

    combined_path = (
        output_dir
        / (
            "raw_vs_umi_consensus_"
            "identity.tsv"
        )
    )

    summary_path = (
        output_dir
        / "processing_summary.tsv"
    )

    combine_outputs(
        summaries,
        combined_path,
    )

    write_processing_summary(
        summaries,
        summary_path,
    )

    total_processed = sum(
        x["processed_raw"]
        for x in summaries
    )

    total_aligned = sum(
        x["aligned"]
        for x in summaries
    )

    total_missing_mapping = sum(
        x["missing_mapping"]
        for x in summaries
    )

    total_missing_consensus = sum(
        x["missing_consensus"]
        for x in summaries
    )

    total_q_missing = sum(
        x["qscore_missing_fastq"]
        for x in summaries
    )

    total_q_not_exact = sum(
        x["qscore_not_exact"]
        for x in summaries
    )

    print()
    print("=" * 70)

    print(
        f"Successfully processed: "
        f"{len(summaries)} / "
        f"{len(barcodes)} barcodes"
    )

    print(
        f"Processed TD3 reads: "
        f"{total_processed:,}"
    )

    print(
        f"Aligned pairs:       "
        f"{total_aligned:,}"
    )

    print(
        f"Missing mapping:      "
        f"{total_missing_mapping:,}"
    )

    print(
        f"Missing consensus:    "
        f"{total_missing_consensus:,}"
    )

    print(
        f"Q FASTQ missing:      "
        f"{total_q_missing:,}"
    )

    print(
        f"Q region not exact:   "
        f"{total_q_not_exact:,}"
    )

    print()
    print(
        f"Combined output:\n"
        f"{combined_path}"
    )

    print()
    print(
        f"Processing summary:\n"
        f"{summary_path}"
    )

if __name__ == "__main__":

    main()
