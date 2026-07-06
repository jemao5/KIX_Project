#!/usr/bin/env python3
"""
Stage 1 of the helix-building pipeline.

Reads peptide sequences and names from an xlsx file and writes them to a
plain-text file (one "name<TAB>sequence" per line) that the Schrodinger
helix-builder can consume.

Run this with your OWN conda environment (the one with pandas/openpyxl),
NOT through the Schrodinger wrapper.

Example:
    python extract_sequences.py peptides.xlsx sequences.tsv \
        --seq-col Sequence --name-col Name
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract peptide sequences/names from an xlsx into a TSV."
    )
    parser.add_argument("xlsx_path", help="Path to the input .xlsx file")
    parser.add_argument("out_tsv", help="Path to the output .tsv file")
    parser.add_argument(
        "--seq-col",
        default="Sequence",
        help="Name of the column containing peptide sequences (default: Sequence)",
    )
    parser.add_argument(
        "--name-col",
        default="Name",
        help="Name of the column containing peptide names/IDs (default: Name)",
    )
    parser.add_argument(
        "--sheet",
        default=0,
        help="Sheet name or index to read (default: first sheet)",
    )
    return parser.parse_args()


def sanitize_name(raw_name):
    """Make a name safe to use as a filename (no spaces/slashes/odd chars)."""
    name = str(raw_name).strip()
    safe = "".join(c if (c.isalnum() or c in "-_.") else "_" for c in name)
    return safe or "unnamed"


def sanitize_sequence(raw_seq):
    """Strip whitespace and uppercase; leave letters as-is for the builder."""
    return "".join(str(raw_seq).split()).upper()


def main():
    args = parse_args()

    xlsx_path = Path(args.xlsx_path)
    if not xlsx_path.exists():
        sys.exit(f"ERROR: input file not found: {xlsx_path}")

    # read_excel accepts an int index or a string sheet name
    sheet = args.sheet
    try:
        sheet = int(sheet)
    except (TypeError, ValueError):
        pass

    df = pd.read_excel(xlsx_path, sheet_name=sheet)

    # Validate the requested columns exist; show what's available if not.
    for col in (args.seq_col, args.name_col):
        if col not in df.columns:
            sys.exit(
                f"ERROR: column '{col}' not found. "
                f"Available columns: {list(df.columns)}"
            )

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    written = 0
    seen_names = {}
    with out_tsv.open("w") as fout:
        for _, row in df.iterrows():
            raw_seq = row[args.seq_col]
            raw_name = row[args.name_col]

            # Skip blank rows (NaN sequences)
            if pd.isna(raw_seq) or str(raw_seq).strip() == "":
                continue

            seq = sanitize_sequence(raw_seq)
            name = sanitize_name(raw_name)

            # Guard against duplicate names clobbering each other's PDBs
            if name in seen_names:
                seen_names[name] += 1
                name = f"{name}_{seen_names[name]}"
            else:
                seen_names[name] = 0

            fout.write(f"{name}\t{seq}\n")
            written += 1

    print(f"Wrote {written} sequences to {out_tsv}")


if __name__ == "__main__":
    main()