#!/usr/bin/env python3
"""
Classify EACH peptide chain of a cofold independently against both KIX faces.

`full_library_face_determination.py` assumes exactly one peptide (chain B). The
dual-face experiment folds KIX with TWO copies (chains B and C) and needs to know
which face each copy landed on -- so this reports a call per chain and then a
per-structure verdict:

    DUAL     one copy on c-Myb, the other on MLL   <- the result being sought
    SAME     both copies on the same face          (competition, not dual binding)
    PARTIAL  one placed, the other on 'neither'
    NEITHER  no productive placement

Contact logic, cutoffs and the face residue sets are IDENTICAL to
full_library_face_determination.py (5.0 A heavy atom, >=2 contacting residues) --
those are contact-based and agreed with the literature on all 24 controls. Run
with `--peptide-chains B` to reproduce the single-peptide behaviour exactly;
that equivalence is the regression test.
"""
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch

CMYB_FACE = {14, 18, 21, 65, 69, 72, 73, 76}
MLL_FACE = {27, 39, 43, 46, 71, 75, 79}
CUTOFF = 5.0
MIN_CONTACTS = 2


def parse_args():
    p = argparse.ArgumentParser(description="Per-chain KIX face assignment")
    p.add_argument("input_glob")
    p.add_argument("output_path")
    p.add_argument("--target-chain", default="A")
    p.add_argument("--peptide-chains", default="B,C",
                   help="Comma-separated. Use 'B' to reproduce the single-peptide script.")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--name-from", choices=["parent", "stem"], default="parent")
    return p.parse_args()


def _parser_for(path):
    return MMCIFParser(QUIET=True) if str(path).lower().endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)


def call_one(job):
    path, target, pep_chains = job
    try:
        model = _parser_for(path).get_structure("s", path)[0]
        kix = model[target]
        out = {}
        for ch in pep_chains:
            if ch not in model:
                out[ch] = (None, None, "absent")
                continue
            atoms = [a for a in model[ch].get_atoms() if a.element != "H"]
            ns = NeighborSearch(atoms)

            def count(face):
                n = 0
                for rid in face:
                    if rid not in kix:
                        continue
                    for atm in kix[rid]:
                        if atm.element == "H":
                            continue
                        if ns.search(atm.coord, CUTOFF):
                            n += 1
                            break
                return n

            c, m = count(CMYB_FACE), count(MLL_FACE)
            hc, hm = c >= MIN_CONTACTS, m >= MIN_CONTACTS
            call = "both" if hc and hm else "cmyb" if hc else "mll" if hm else "neither"
            out[ch] = (c, m, call)
        return path, out, None
    except Exception as exc:
        return path, None, str(exc)


def verdict(calls):
    """Per-structure summary from the individual chain calls."""
    faces = [c for _, _, c in calls]
    real = [f for f in faces if f in ("cmyb", "mll", "both")]
    if len(real) < len(faces):
        return "PARTIAL" if real else "NEITHER"
    if len(set(faces)) == 1 and faces[0] in ("cmyb", "mll"):
        return "SAME"
    if {"cmyb", "mll"} <= set(faces) or "both" in faces:
        return "DUAL"
    return "SAME"


def main():
    a = parse_args()
    pep = [c.strip() for c in a.peptide_chains.split(",") if c.strip()]
    files = sorted(glob.glob(a.input_glob))
    if not files:
        raise SystemExit(f"No structures matched: {a.input_glob}")
    print(f"{len(files)} structures; target={a.target_chain}; peptide chains={pep}")

    jobs = [(f, a.target_chain, pep) for f in files]
    results, errors = {}, []
    if a.workers == 1:
        for j in jobs:
            p, o, e = call_one(j)
            (errors.append((p, e)) if e else results.__setitem__(p, o))
    else:
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(call_one, j): j[0] for j in jobs}
            for fut in as_completed(futs):
                p, o, e = fut.result()
                (errors.append((p, e)) if e else results.__setitem__(p, o))

    single = len(pep) == 1
    with open(a.output_path, "w") as fh:
        if single:
            # byte-compatible with full_library_face_determination.py
            fh.write("name\tcmyb_contacts\tmll_contacts\tface_call\n")
        else:
            cols = ["name"] + [f"chain{c}_{k}" for c in pep
                               for k in ("cmyb_contacts", "mll_contacts", "face")] + ["verdict"]
            fh.write("\t".join(cols) + "\n")
        for path, o in sorted(results.items()):
            name = Path(path).parent.name if a.name_from == "parent" else Path(path).stem
            if single:
                c, m, call = o[pep[0]]
                fh.write(f"{name}\t{c}\t{m}\t{call}\n")
            else:
                row = [name]
                for ch in pep:
                    c, m, call = o[ch]
                    row += [str(c), str(m), call]
                row.append(verdict([o[ch] for ch in pep]))
                fh.write("\t".join(row) + "\n")

    print(f"Wrote {len(results)} rows to {a.output_path} ({len(errors)} errors)")
    if not single:
        from collections import Counter
        v = Counter(verdict([o[ch] for ch in pep]) for o in results.values())
        for k, n in v.most_common():
            print(f"  {k}: {n}")
    for p, e in errors[:5]:
        print(f"  ERROR {Path(p).name}: {e}")


if __name__ == "__main__":
    main()
