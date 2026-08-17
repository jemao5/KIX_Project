#!/usr/bin/env python3
"""Atom-level KIX<->peptide interactions, straight from Schrodinger.

`schrodinger_calc_hbond.find_all_interactions` returns Schrodinger atom objects
but keeps only residue identity, so the ATOM NAMES -- what a figure needs to draw
a bond -- were discarded. This re-runs the identical call on the structures used
for figures and keeps them.

Why this exists: ChimeraX's `hbonds` disagrees badly with Schrodinger on these
interfaces. Schrodinger measures HYDROGEN-to-acceptor <= 2.8 A on prepwizard
output that HAS hydrogens; the deposited/aligned PDBs have none, so ChimeraX
falls back to heavy-atom geometry and finds ~1 inter-chain H-bond where
Schrodinger finds several. Drawing explicit pseudobonds from these numbers makes
the figures match the CSVs exactly.

Run with Schrodinger's python:
    /share/apps/images/run-schrodinger-2025.4.bash run python3 extract_atom_interactions.py
"""
import os, sys, pickle, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schrodinger import structure
from schrodinger.structutils import analyze
from schrodinger.structutils.interactions import hbond
from schrodinger.structutils.interactions.pi import find_pi_pi_interactions

P = os.path.dirname(os.path.abspath(__file__))


def ring_center_atom(st, idxs):
    """Ring atom nearest the ring centroid -- a stable anchor for drawing."""
    pts = [(st.atom[i].x, st.atom[i].y, st.atom[i].z) for i in idxs]
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    cz = sum(p[2] for p in pts) / len(pts)
    best, bd = None, 1e9
    for i, (x, y, z) in zip(idxs, pts):
        d = (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2
        if d < bd:
            best, bd = i, d
    return st.atom[best]


def heavy(a):
    """Schrodinger returns the DONOR HYDROGEN; the aligned PDBs have no
    hydrogens, so resolve it to the heavy atom it is bonded to (HH12 -> NH1,
    H1 -> N, ...). Uses the real bond graph rather than a name table."""
    if a.element != "H":
        return a
    for b in a.bonded_atoms:
        if b.element != "H":
            return b
    return a


def atoms_for(path, c1, c2):
    """Same call as schrodinger_calc_hbond, but returning atom names."""
    st = next(structure.StructureReader(path))
    s1 = st.extract(analyze.evaluate_asl(st, f'chain.name "{c1}"'))
    s2 = st.extract(analyze.evaluate_asl(st, f'chain.name "{c2}"'))
    out = []
    for d, a in hbond.get_hydrogen_bonds(s1, st2=s2, max_dist=2.8,
                                         min_donor_angle=120.0,
                                         min_acceptor_angle=90.0):
        # orient so the KIX side is first, matching hbond_hit_num.py
        d, a = heavy(d), heavy(a)
        k, p = (d, a) if d.chain == c1 else (a, d)
        out.append(dict(interaction="hbond",
                        kix_resnum=k.resnum, kix_resname=k.pdbres.strip(),
                        kix_atom=k.pdbname.strip(),
                        pep_resnum=p.resnum, pep_resname=p.pdbres.strip(),
                        pep_atom=p.pdbname.strip(),
                        donor_is_kix=(d.chain == c1)))
    for pi in find_pi_pi_interactions(s1, struct2=s2):
        a1 = ring_center_atom(pi.struct1, pi.ring1.atoms)
        a2 = ring_center_atom(pi.struct2, pi.ring2.atoms)
        k, p = (a1, a2) if a1.chain == c1 else (a2, a1)
        out.append(dict(interaction="pi_pi",
                        kix_resnum=k.resnum, kix_resname=k.pdbres.strip(),
                        kix_atom=k.pdbname.strip(),
                        pep_resnum=p.resnum, pep_resname=p.pdbres.strip(),
                        pep_atom=p.pdbname.strip(), donor_is_kix=""))
    return out


def main():
    import pandas as pd
    k = pd.read_csv(f"{P}/full_library_all_metrics/key_residues_top10.csv")

    # name -> prepared .maegz, from the pickles that already recorded them
    lib = {os.path.basename(x["file"]).replace("_model_0_prepared.maegz", ""): x["file"]
           for x in pickle.load(open(f"{P}/interactions.pkl", "rb"))}
    dual = {os.path.basename(x["file"]).replace("_prepared.maegz", ""): x["file"]
            for x in pickle.load(open(f"{P}/dual_cofold/copies_interactions.pkl", "rb"))}

    rows = []
    for st_, name, pep, chain in k[["set", "name", "peptide", "chain"]].drop_duplicates().itertuples(index=False):
        if st_ == "top10_dual":
            src = dual.get(name)          # name is <peptide>__copyB / __copyC
        else:
            src = lib.get(name)
        if not src or not os.path.exists(src):
            print("MISSING prepared:", name); continue
        # sub-complexes always store the peptide as chain B; `chain` is where it
        # sits in the full aligned structure (B or C)
        for r in atoms_for(src, "A", "B"):
            r.update(structure_name=name, target_chain=chain, set=st_, source="model")
            rows.append(r)
        print(f"  {name}: {len([r for r in rows if r['structure_name']==name])} interactions")

    # crystal: 2AGH KIX = chain B, peptides A (c-Myb) and C (MLL); KIX +585
    for face, pc in (("cmyb", "A"), ("mll", "C")):
        for r in atoms_for(f"{P}/2agh_model1_prepared.maegz", "B", pc):
            r["kix_resnum"] -= 585        # -> pipeline numbering, as the .cxc renumbers
            r.update(structure_name=f"2AGH_native_{face}", target_chain=pc,
                     set="crystal", source="crystal")
            rows.append(r)
        print(f"  2AGH {face}: done")

    cols = ["set", "source", "structure_name", "target_chain", "interaction",
            "kix_resname", "kix_resnum", "kix_atom",
            "pep_resname", "pep_resnum", "pep_atom", "donor_is_kix"]
    out = f"{P}/full_library_all_metrics/atom_level_interactions.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows: w.writerow({c: r.get(c, "") for c in cols})
    print(f"\nwrote {out}: {len(rows)} interactions")


if __name__ == "__main__":
    main()
