#!/usr/bin/env python3
"""
Produce DSSP-safe copies of the ncAA-grafted structures.

DSSP fails on 17 of 24 swapped controls: it only recognises residues in its
amino-acid table, and our grafted codes (M2F/CHX/CY5/ABU/NLE/CIA) are not in it.
BCS passes only because it happens to be a real PDB CCD code.

DSSP's secondary-structure assignment depends solely on backbone N/CA/C/O
geometry -- side chains play no part in the H-bond energy calculation. So for
DSSP input only, each grafted residue is rewritten as GLY carrying just its
N/CA/C/O atoms. The backbone coordinates are untouched, so the helicity DSSP
reports is identical to what it would report for the real residue if it knew it.

This transformation is used for NOTHING else -- BindCraft, Schrodinger and the
face call all read the real structures with the real residues.

Validation (see --verify): for peptides DSSP already handles, the substituted
copy must return exactly the same chain_b_helix_fraction as the original.
"""
import argparse
import glob
import os
from pathlib import Path

NCAA_CODES = {"BCS", "M2F", "CHX", "CY5", "ABU", "NLE", "CIA"}
BACKBONE = {"N", "CA", "C", "O"}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--pattern", default="*_ncaa.pdb")
    return p.parse_args()


def convert(src, dst):
    out, n_sub = [], 0
    for line in open(src):
        if line.startswith(("ATOM", "HETATM")):
            line = line.rstrip("\n").ljust(80)
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            if resname in NCAA_CODES:
                if atom not in BACKBONE:
                    continue                       # drop side chain
                line = "ATOM  " + line[6:17] + "GLY" + line[20:]
                n_sub += 1
            out.append(line.rstrip() + "\n")
        elif line.startswith("TER"):
            out.append(line)
    with open(dst, "w") as f:
        f.writelines(out)
        f.write("END\n")
    return n_sub


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.in_dir, args.pattern)))
    total = 0
    for src in files:
        dst = out_dir / Path(src).name
        n = convert(src, dst)
        total += n
    print(f"Wrote {len(files)} DSSP-input structures to {out_dir} "
          f"({total} ncAA backbone atoms relabelled GLY)")


if __name__ == "__main__":
    main()
