#!/usr/bin/env python3

from pathlib import Path
import sys
import argparse
import pandas as pd
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description= r"Takes csv of sequences only and generates a tsv file of {name}<tab>{sequence}"
    )
    parser.add_argument("xlsx_path", help="path to input xlsx file")
    parser.add_argument(
        "--out_path",
        default="./out_full_library.tsv",
        help="output tsv path"
    )
    parser.add_argument(
        "--name_pref",
        default="seq_",
        help=r"prefix to be placed before number in each name {pref}{#}"
    )
    parser.add_argument(
        "--sheet",
        default="Full Library_Filtered",
        help="Sheet name to read (default: 'Full Library_Filtered')",
    )
    parser.add_argument(
        "--seq-col",
        default="SEQUENCE",
        help="Name of the column containing peptide sequences (default: SEQUENCE)",
    )
    return parser.parse_args()

def sanitize_sequence(raw_seq):
    """Strip whitespace and uppercase; leave letters as-is for the builder."""
    return "".join(str(raw_seq).split()).upper()

def check_seq(seq):
    result = re.search(r"\A[ARNDCEQGHILKMFPSTWYV]{10}\Z", seq)
    if result:
        return True
    else:
        return False


def main():
    total = 0
    skipped_blank = 0
    repeat_filtered = 0
    failed_sequence_filter = 0
    written = 0


    sequences = set()
    args = parse_args()
    xlsx_path = Path(args.xlsx_path)
    sheet = args.sheet
    try:
        sheet = int(sheet)
    except (ValueError, TypeError):
        pass
    if not xlsx_path.exists():
        sys.exit(f"ERROR: input file not found: {xlsx_path}")
    
    df = pd.read_excel(xlsx_path, sheet_name=sheet, engine='calamine')
    out_tsv = Path(args.out_path)

    with out_tsv.open("w") as fout:
        for orig_idx, row in df.iterrows():
            total += 1
            raw_seq = row[args.seq_col]

            # Skip blank rows (NaN sequences)
            if pd.isna(raw_seq) or str(raw_seq).strip() == "":
                skipped_blank += 1
                continue

            seq = sanitize_sequence(raw_seq)
            
            if check_seq(seq):
                if seq not in sequences:
                    name = f"{args.name_pref}{orig_idx}"
                    fout.write(f"{name}\t{seq}\n")
                    written += 1
                    sequences.add(seq)
                else:
                    repeat_filtered += 1
            else:
                failed_sequence_filter += 1

            


    print(f"Wrote {written} unique, full, sequences to {out_tsv}")
    print(f"Total rows: {total}, blank: {skipped_blank}, "
      f"repeat filtered: {repeat_filtered}, failed sequence filter: {failed_sequence_filter}, written: {written}")

if __name__ == "__main__":
    main()