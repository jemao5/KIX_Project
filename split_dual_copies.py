#!/usr/bin/env python3
"""
Split each dual cofold (KIX + 2 peptide copies) into two ordinary 2-chain
complexes so the STANDARD metric pipeline can run on them unchanged.

    <name>_model_0.cif  (chains A,B,C)  ->  <name>__copyB.pdb  (A + B)
                                            <name>__copyC.pdb  (A + C renamed B)

Renaming C to B matters: every downstream tool assumes KIX = chain A and
peptide = chain B (BindCraft --binder-chain B, schrodinger_calc_hbond's
"A"/"B" pair, full_library_face_determination's model["B"]). It also avoids
score_peptide.py's sanitize_structure, which deletes chains C and D by default
and would silently drop the second copy.

The two sub-complexes keep the ORIGINAL coordinates -- nothing is re-folded or
re-minimised, so each measures that copy's interface exactly as Boltz placed it.
"""
import argparse, glob, os
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO, Select


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input_glob")
    p.add_argument("out_dir")
    p.add_argument("--target-chain", default="A")
    p.add_argument("--peptide-chains", default="B,C")
    return p.parse_args()


class Keep(Select):
    def __init__(self, keep):
        self.keep = set(keep)
    def accept_chain(self, chain):
        return chain.id in self.keep
    def accept_atom(self, atom):
        return atom.element != "H"


def main():
    a = parse_args()
    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    pep = [c.strip() for c in a.peptide_chains.split(",") if c.strip()]
    files = sorted(glob.glob(a.input_glob))
    if not files:
        raise SystemExit(f"no structures matched {a.input_glob}")

    io = PDBIO()
    n = 0
    for f in files:
        name = Path(f).stem.replace("_model_0", "")
        st = MMCIFParser(QUIET=True).get_structure(name, f)
        model = st[0]
        for ch in pep:
            if ch not in model:
                print(f"  WARNING {name}: chain {ch} absent")
                continue
            # rebuild with the peptide chain relabelled B
            import copy as _copy
            m2 = _copy.deepcopy(model)
            for c in list(m2):
                if c.id not in (a.target_chain, ch):
                    m2.detach_child(c.id)
            if ch != "B":
                if "B" in [c.id for c in m2]:
                    m2["B"].id = "_tmp"        # free up 'B' first
                m2[ch].id = "B"
            st2 = st.__class__(name); st2.add(m2)
            io.set_structure(st2)
            io.save(str(out / f"{name}__copy{ch}.pdb"), Keep([a.target_chain, "B"]))
            n += 1
    print(f"wrote {n} sub-complexes to {out}")


if __name__ == "__main__":
    main()
