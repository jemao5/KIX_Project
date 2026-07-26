#!/usr/bin/env python3
"""
Drive residuex_swap.py across the whole literature control set.

Reads `ncaa_swaps` from control_list.csv, runs the ResidueX graft for each
peptide that needs one, and copies the three all-canonical controls
(WT_MLL, MLL1, cMyb_native) through unchanged so downstream stages see a
uniform directory.

Writes control_run/ncaa_swap_summary.csv with per-peptide status, so a failure
is visible rather than silently producing a short structure.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path("/scratch/jem9759/ZhangWork/KIX_Project")
PY = "/scratch/jem9759/envs/residuex/bin/python3"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--control-list", default=str(ROOT / "control_list.csv"))
    p.add_argument("--smiles-csv", default=str(ROOT / "control_run/ncaa_smiles.csv"))
    p.add_argument("--boltz-dir",
                   default=str(ROOT / "control_run/boltz_out/boltz_results_yamls/predictions"))
    p.add_argument("--out-dir", default=str(ROOT / "control_run/ncaa_structures"))
    p.add_argument("--only", default=None, help="Comma-separated names to run")
    return p.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.control_list)
    if args.only:
        keep = {n.strip() for n in args.only.split(",")}
        df = df[df["name"].isin(keep)]

    rows = []
    for _, r in df.iterrows():
        name = r["name"].strip()
        src = Path(args.boltz_dir) / name / f"{name}_model_0.pdb"
        swaps_raw = str(r["ncaa_swaps"]).strip()
        dest = out_dir / f"{name}_ncaa.pdb"

        if not src.exists():
            print(f"[{name}] MISSING cofold {src}")
            rows.append({"name": name, "status": "missing_cofold", "n_swaps": 0})
            continue

        if swaps_raw in ("none", "nan", ""):
            shutil.copy(src, dest)
            print(f"[{name}] no ncAA -- copied through unchanged")
            rows.append({"name": name, "status": "copied_no_swap", "n_swaps": 0})
            continue

        specs = []
        for part in swaps_raw.split(";"):
            m = re.match(r"^([A-Z])(\d+)->(.+)$", part.strip())
            if not m:
                print(f"[{name}] UNPARSED swap spec: {part!r}")
                specs = None
                break
            specs.append(f"{m.group(2)}:{m.group(3).strip()}")
        if specs is None:
            rows.append({"name": name, "status": "unparsed_swaps", "n_swaps": 0})
            continue

        cmd = [PY, str(ROOT / "residuex_swap.py"), "--name", name,
               "--complex-pdb", str(src), "--smiles-csv", args.smiles_csv,
               "--swaps", ",".join(specs), "--out-dir", str(out_dir)]
        print(f"[{name}] {len(specs)} swap(s): {','.join(specs)}")
        p = subprocess.run(cmd, capture_output=True, text=True)
        ok = dest.exists()
        if not ok:
            tail = [l for l in p.stdout.splitlines() + p.stderr.splitlines()
                    if "Error" in l or "error" in l or "Traceback" in l
                    or "ValueError" in l or "refus" in l]
            print(f"[{name}] FAILED: {tail[-1] if tail else 'no output file'}")
        rows.append({"name": name,
                     "status": "ok" if ok else "failed",
                     "n_swaps": len(specs),
                     "swaps": ",".join(specs)})

    summary = pd.DataFrame(rows)
    out_csv = Path(args.out_dir).parent / "ncaa_swap_summary.csv"
    summary.to_csv(out_csv, index=False)
    print("\n=== SUMMARY ===")
    print(summary["status"].value_counts().to_string())
    bad = summary[~summary["status"].isin(["ok", "copied_no_swap"])]
    if len(bad):
        print("\nPROBLEMS:")
        print(bad.to_string(index=False))
    print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
