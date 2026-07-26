#!/usr/bin/env python3
"""
Stage 1 for the literature control set: control_list.csv -> name<TAB>sequence TSV.

The control peptides carry non-canonical residues (Bcs, 2meF, CyHex, ...) in
`sequence_display`; `boltz_sequence` is the canonical-amino-acid stand-in that
Boltz can actually fold, so that is the column the pipeline uses.

Also writes control_face_truth.tsv (name -> literature face + positive/negative
label + measured affinity), which the scoring step joins back on.
"""
import argparse
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent


def parse_args():
    p = argparse.ArgumentParser(description="control_list.csv -> pipeline TSV")
    p.add_argument("csv_path", nargs="?",
                   default=str(SCRIPT_DIR / "control_list.csv"))
    p.add_argument("--out_dir", default=str(SCRIPT_DIR / "control_run"))
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.csv_path)
    df["name"] = df["name"].str.strip()
    df["boltz_sequence"] = df["boltz_sequence"].str.strip().str.upper()

    bad = df[~df["boltz_sequence"].str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")]
    if len(bad):
        raise SystemExit(f"Non-canonical letters in boltz_sequence:\n{bad[['name','boltz_sequence']]}")

    dupe_names = df["name"][df["name"].duplicated()].tolist()
    if dupe_names:
        raise SystemExit(f"Duplicate control names: {dupe_names}")

    tsv = out_dir / "control_list.tsv"
    df[["name", "boltz_sequence"]].to_csv(tsv, sep="\t", header=False, index=False)

    truth = out_dir / "control_face_truth.tsv"
    df[["name", "face", "label", "kd_or_ki_uM", "measurement_type",
        "sequence_display", "ncaa_swaps", "source"]].to_csv(truth, sep="\t", index=False)

    # hbond_hit_num.py routes each structure to a face's hit-residue list via a
    # `face_call` column. Controls are counted against the face they are known
    # to bind, so a negative control that mislocalizes still gets scored on the
    # right pocket instead of silently returning hit_num = 0. The independent
    # geometric call from stage 9 is kept separately for agreement checking.
    lit_face = out_dir / "control_face_literature.tsv"
    df[["name", "face"]].rename(columns={"face": "face_call"}).to_csv(
        lit_face, sep="\t", index=False)

    n_uniq = df["boltz_sequence"].nunique()
    print(f"Wrote {len(df)} controls to {tsv} ({n_uniq} unique sequences)")
    print(f"Wrote literature annotations to {truth}")
    print(df.groupby(["face", "label"]).size().to_string())


if __name__ == "__main__":
    main()
