#!/usr/bin/env python3
"""
Pick the top N per face for the dual-face cofold experiment.

Ranks on `priority_score_v2_pde` -- the best-validated blend on the literature
controls (pooled positives-vs-negatives AUC 0.788, p=0.011, vs 0.598 for the
current default). It is near-identical to `_pde` for selection purposes
(rank corr +0.985/+0.983, top-40 overlap 38/40 and 37/40).

Records source_face and within-face rank so any result traces back to where the
sequence came from.
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
METRICS = ROOT / "full_library_all_metrics"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=40, help="top N per face")
    p.add_argument("--score", default="priority_score_v2_pde")
    p.add_argument("--out", default=str(ROOT / "dual_cofold" / "dual_list.tsv"))
    p.add_argument("--meta", default=str(ROOT / "dual_cofold" / "dual_selection.csv"))
    return p.parse_args()


def main():
    a = parse_args()
    picked = []
    for face in ("cmyb", "mll"):
        d = pd.read_csv(METRICS / f"{face}_candidates.csv")
        if a.score not in d.columns:
            raise SystemExit(f"{a.score} missing from {face}_candidates.csv "
                             f"-- run analyze_and_score_all_metrics.py first")
        top = d.nlargest(a.n, a.score).reset_index(drop=True)
        top["source_face"] = face
        top["face_rank"] = top.index + 1
        picked.append(top)
    sel = pd.concat(picked, ignore_index=True)

    dup = sel["name"][sel["name"].duplicated()]
    if len(dup):
        raise SystemExit(f"duplicate names across faces: {list(dup)}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    sel[["name", "Sequence"]].to_csv(out, sep="\t", header=False, index=False)
    cols = ["name", "Sequence", "source_face", "face_rank", a.score,
            "hit_num", "hit_num_v2", "helix_score", "complex_pde",
            "interface_dG", "protein_iptm", "count"]
    sel[[c for c in cols if c in sel.columns]].to_csv(a.meta, index=False)

    print(f"selected {len(sel)} ({a.n} per face) ranked on {a.score}")
    print(f"  unique sequences: {sel['Sequence'].nunique()}")
    print(f"  lengths: {sorted(sel['Sequence'].str.len().unique())}")
    print(f"  -> {out}")
    print(f"  -> {a.meta}")


if __name__ == "__main__":
    main()
