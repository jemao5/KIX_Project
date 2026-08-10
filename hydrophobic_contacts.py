#!/usr/bin/env python3
"""
Count, per structure, how many APOLAR hit residues the peptide is touching.

`hit_num` counts only H-bonds and pi-pi stacking, so hit residues with an apolar
side chain (MLL: PHE 27, LEU 43, ILE 75, LEU 79) can never contribute -- the
hydrophobic packing that drives MLL binding is invisible to it. This script
supplies the missing term; `hbond_hit_num.py --hydrophobic` adds the two
together as `hit_num_v2`.

Measured on the 23 MLL literature controls (raw, no length correction):
    hbond+pipi only        pos 1.00 / neg 0.27   MWU p=0.037  rho=-0.519
    hbond+pipi + apolar    pos 5.00 / neg 4.09   MWU p=0.007  rho=-0.622

Counts are emitted for BOTH faces; hbond_hit_num.py selects the one matching
each peptide's face_call.

Structure of this script deliberately mirrors full_library_face_determination.py
(same glob/ProcessPoolExecutor/--workers pattern, same BioPython parsing).

Usage:
    python3 hydrophobic_contacts.py "<glob>" out.tsv --workers 16 [--name-from parent|stem]
"""
import argparse
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kix_scoring import apolar_hit_residues

CUTOFF = 4.5          # heavy-atom; see --cutoff
TARGET_CHAIN = "A"    # KIX
PEPTIDE_CHAIN = "B"


def parse_args():
    p = argparse.ArgumentParser(description="Count apolar hit-residue contacts per structure")
    p.add_argument("input_glob")
    p.add_argument("output_path")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--cutoff", type=float, default=CUTOFF)
    p.add_argument("--residue-set", default="original")
    p.add_argument("--name-from", choices=["parent", "stem"], default="parent",
                   help="'parent' matches full_library_face_determination.py (the "
                        "Full_Library_Hit_N folder); 'stem' for a flat directory.")
    return p.parse_args()


def _parser_for(path):
    return MMCIFParser(QUIET=True) if str(path).lower().endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)


def count_contacts(args_tuple):
    path, cutoff, res_set = args_tuple
    try:
        model = _parser_for(path).get_structure("s", path)[0]
        kix, pep = model[TARGET_CHAIN], model[PEPTIDE_CHAIN]
        pep_atoms = [a for a in pep.get_atoms() if a.element != "H"]
        ns = NeighborSearch(pep_atoms)
        out = {}
        for face in ("cmyb", "mll"):
            n = 0
            for rid in apolar_hit_residues(res_set, face):
                if rid not in kix:
                    continue
                for atm in kix[rid]:
                    if atm.element == "H":
                        continue
                    if ns.search(atm.coord, cutoff):
                        n += 1
                        break          # residue counted once, not per atom
            out[face] = n
        return path, out, None
    except Exception as exc:
        return path, None, str(exc)


def main():
    args = parse_args()
    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit(f"No structures matched: {args.input_glob}")
    print(f"{len(files)} structures; residue set {args.residue_set}; cutoff {args.cutoff} A")
    for face in ("cmyb", "mll"):
        print(f"  apolar {face}: {apolar_hit_residues(args.residue_set, face)}")

    jobs = [(f, args.cutoff, args.residue_set) for f in files]
    results, errors = {}, []
    if args.workers == 1:
        for j in jobs:
            p, o, e = count_contacts(j)
            (errors.append((p, e)) if e else results.__setitem__(p, o))
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(count_contacts, j): j[0] for j in jobs}
            for fut in as_completed(futs):
                p, o, e = fut.result()
                (errors.append((p, e)) if e else results.__setitem__(p, o))

    with open(args.output_path, "w") as fh:
        fh.write("name\tcmyb_apolar_contacts\tmll_apolar_contacts\n")
        for path, o in sorted(results.items()):
            name = Path(path).parent.name if args.name_from == "parent" else Path(path).stem
            fh.write(f"{name}\t{o['cmyb']}\t{o['mll']}\n")
    print(f"Wrote {len(results)} rows to {args.output_path} ({len(errors)} errors)")
    for p, e in errors[:5]:
        print(f"  ERROR {os.path.basename(p)}: {e}")


if __name__ == "__main__":
    main()
