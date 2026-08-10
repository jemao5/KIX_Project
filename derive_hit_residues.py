#!/usr/bin/env python3
"""
Derive the per-face `hit_num` residue sets empirically from the 2AGH complex.

Why: the hand-picked MLL hit list contains 4 residues (PHE 27, LEU 43, ILE 75,
LEU 79) that are hydrophobic and physically cannot hydrogen-bond, yet `hit_num`
counts only H-bonds and pi-pi. Across 23 MLL controls just 15 of 89 observed
H-bonds landed on hit residues. This script instead asks which KIX residues
actually contact the NATIVE peptides in 2AGH, using the same Schrodinger
workflow that computes hit_num itself.

2AGH model 1 chains:  A = c-Myb peptide (291-315)
                      B = KIX           (586-672)
                      C = MLL peptide   (839-869)

Numbering: 2AGH uses CBP numbering; the pipeline uses 1-87. offset = 585, and
chain B's sequence is identical to KIX_SEQUENCE (asserted below).

All three interaction types are reported SEPARATELY (H-bond, pi-pi, heavy-atom
contact) so the definition of a "hit" can be chosen after seeing the data.

Run with Schrodinger's python:
    /share/apps/images/run-schrodinger-2025.4.bash run python3 derive_hit_residues.py
"""
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger_calc_hbond import find_all_interactions   # reuse, do not reimplement

ROOT = os.path.dirname(os.path.abspath(__file__))
PREPARED = os.path.join(ROOT, "2agh_model1_prepared.maegz")
OFFSET = 585                 # 2AGH resnum - pipeline resnum
KIX_CHAIN = "B"
FACES = {"cmyb": "A", "mll": "C"}
HEAVY_CUTOFF = 5.0           # matches full_library_face_determination.py

KIX_SEQUENCE = ("GVRKGWHEHVTQDLRSHLVHKLVQAIFPTPDPAALKDRRMENLVAYAKKVEGDMYESANSRD"
                "EYYHLLAEKIYKIQKELEEKRRSRL")
THREE_TO_ONE = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
                "HID":"H","HIE":"H","HIP":"H"}

CURRENT = {"cmyb": {14, 18, 21, 65, 69, 72, 73, 76},
           "mll":  {27, 39, 43, 46, 71, 75, 79}}


def parse_res(tag):
    """'LEU 599 B' -> ('LEU', 599, 'B')"""
    name, num, chain = tag.split()
    return name, int(num), chain


def heavy_contacts(st, kix_chain, pep_chain):
    """Per-KIX-residue count of peptide heavy atoms within HEAVY_CUTOFF."""
    kix = [a for a in st.atom
           if a.chain == kix_chain and a.element != "H"]
    pep = [a for a in st.atom
           if a.chain == pep_chain and a.element != "H"]
    counts = defaultdict(int)
    c2 = HEAVY_CUTOFF ** 2
    for a in kix:
        ax, ay, az = a.x, a.y, a.z
        n = 0
        for b in pep:
            dx = ax - b.x; dy = ay - b.y; dz = az - b.z
            if dx*dx + dy*dy + dz*dz <= c2:
                n += 1
        if n:
            counts[(a.pdbres.strip(), a.resnum)] += n
    return counts


def main():
    if not os.path.exists(PREPARED):
        sys.exit("missing %s -- run prepwizard first" % PREPARED)
    st = next(structure.StructureReader(PREPARED))

    chains = sorted({a.chain for a in st.atom})
    print("chains in prepared structure:", chains)

    rows = []
    for face, pep_chain in FACES.items():
        if pep_chain not in chains:
            sys.exit("peptide chain %s missing" % pep_chain)
        # KIX first so hb[0] is the KIX-side residue (matches hbond_hit_num.py)
        hb, pp = find_all_interactions(PREPARED, KIX_CHAIN, pep_chain)
        hb_count, pp_count = defaultdict(int), defaultdict(int)
        for k, _ in hb:
            nm, num, _ch = parse_res(k)
            hb_count[(nm, num)] += 1
        for k, _ in pp:
            nm, num, _ch = parse_res(k)
            pp_count[(nm, num)] += 1
        hv = heavy_contacts(st, KIX_CHAIN, pep_chain)

        for key in set(hb_count) | set(pp_count) | set(hv):
            nm, num2agh = key
            pipeline = num2agh - OFFSET
            # hard fail on any numbering drift
            if not (1 <= pipeline <= len(KIX_SEQUENCE)):
                sys.exit("residue %s %d maps outside 1-87 (got %d)" % (nm, num2agh, pipeline))
            expect = KIX_SEQUENCE[pipeline - 1]
            got = THREE_TO_ONE.get(nm, "?")
            if got != expect:
                sys.exit("NUMBERING MISMATCH at 2AGH %d -> pipeline %d: "
                         "structure has %s (%s), KIX_SEQUENCE has %s"
                         % (num2agh, pipeline, nm, got, expect))
            rows.append({"face": face, "resnum_2agh": num2agh, "resnum_pipeline": pipeline,
                         "resname": nm, "n_hbond": hb_count.get(key, 0),
                         "n_pipi": pp_count.get(key, 0), "n_heavy_contacts": hv.get(key, 0),
                         "in_current_list": pipeline in CURRENT[face]})
        print("%s face: %d H-bonds, %d pi-pi, %d residues with heavy contact"
              % (face, len(hb), len(pp), len(hv)))

    out = os.path.join(ROOT, "hit_residue_derivation.tsv")
    cols = ["face", "resnum_2agh", "resnum_pipeline", "resname",
            "n_hbond", "n_pipi", "n_heavy_contacts", "in_current_list"]
    rows.sort(key=lambda r: (r["face"], -r["n_heavy_contacts"], r["resnum_pipeline"]))
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("\nnumbering check passed for all %d residues (offset %d)" % (len(rows), OFFSET))
    print("wrote", out)


if __name__ == "__main__":
    main()
