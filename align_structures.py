"""
Align cofolded structures onto the 2AGH reference for viewing in ChimeraX.

Post-scoring visualisation prep. Reads one or more candidate CSVs, takes the top
N, structurally aligns each Boltz model_0 onto 2agh_model1.cif, and writes the
aligned structures plus an alignment_data.csv (rmsd/score).

Defaults reproduce the original hardcoded behaviour (top 20 per face from
{cmyb,mll}_candidates.csv in file order -> aligned_structures/). Everything is
now a CLI arg so a different ranking or a different structure set can be aligned
without editing the file:

  # individual binders, current best-validated ranking
  --candidates full_library_all_metrics/cmyb_candidates.csv \
  --candidates full_library_all_metrics/mll_candidates.csv \
  --sort-by priority_score_v2_pde --n 40 --out-dir aligned_structures_v2pde

  # dual cofolds (3-chain: KIX + 2 peptide copies)
  --candidates dual_cofold/dual_shortlist.csv --name-col peptide \
  --sort-by dual_priority_score_v2_pde --n 40 \
  --boltz-root dual_cofold/boltz_out --out-dir dual_cofold/aligned

⚠️ `--sort-by` matters: the original used df.head(n), which silently depended on
whatever order the CSV happened to be saved in. Pass the score explicitly.

Run with Schrodinger's python (needs structalign2):
    /share/apps/images/run-schrodinger-2025.4.bash run python3 align_structures.py ...
"""
import argparse
import json
from pathlib import Path

import pandas as pd
from schrodinger import structure
from schrodinger.structutils import rmsd as sch_rmsd
from schrodinger.structutils import structalign2

ROOT = Path(__file__).resolve().parent
REFERENCE_PATH = str(ROOT / "2agh_model1.cif")


def parse_args():
    p = argparse.ArgumentParser(description="Align cofolds onto 2AGH for viewing")
    p.add_argument("--candidates", action="append", default=None,
                   help="candidate CSV (repeatable). Default: the two library face lists.")
    p.add_argument("--name-col", default="name")
    p.add_argument("--sort-by", default=None,
                   help="score column to rank on. Default: file order (legacy behaviour).")
    p.add_argument("--n", type=int, default=20, help="top N per candidate file")
    p.add_argument("--boltz-root", default=str(ROOT / "boltz_out_full"))
    p.add_argument("--out-dir", default=str(ROOT / "aligned_structures"))
    p.add_argument("--reference", default=REFERENCE_PATH)
    p.add_argument("--align-on-chain", nargs=2, metavar=("REF_CHAIN", "MOBILE_CHAIN"),
                   default=None,
                   help="Superpose explicitly on this chain pair (e.g. B A = 2AGH KIX "
                        "onto our chain A) instead of letting align_many choose the "
                        "correspondence. REQUIRED for multi-chain complexes: with two "
                        "identical peptide copies present, align_many fits a peptide "
                        "instead of KIX on ~1 in 8 structures, leaving KIX 8-24 A off.")
    a = p.parse_args()
    if not a.candidates:
        a.candidates = [str(ROOT / "full_library_all_metrics" / f"{f}_candidates.csv")
                        for f in ("mll", "cmyb")]
    return a


def build_cif_index(root):
    """Map each peptide name -> its Boltz model_0 .cif path, across all chunks."""
    index = {}
    for cif in Path(root).rglob("*_model_0.cif"):
        index[cif.name[:-len("_model_0.cif")]] = cif
    return index


def ca_atoms(st, chain):
    """CA atom indices of one chain, in residue order."""
    return [a.index for a in st.atom
            if a.chain == chain and a.pdbname.strip() == "CA"]


def superpose_on_chain(ref, mobile, ref_chain, mob_chain):
    """Superpose `mobile` onto `ref` using ONE chain's CA atoms, moving the whole
    structure. Returns the CA-RMSD over that chain.

    Used instead of align_many for multi-chain complexes, where the automatic
    correspondence can latch onto the wrong chain."""
    a_ref, a_mob = ca_atoms(ref, ref_chain), ca_atoms(mobile, mob_chain)
    if len(a_ref) != len(a_mob):
        raise ValueError(f"chain {ref_chain} has {len(a_ref)} CA but mobile "
                         f"chain {mob_chain} has {len(a_mob)} -- cannot pair")
    return sch_rmsd.superimpose(ref, a_ref, mobile, a_mob)


def standardize_his_names(st):
    """Rename Schrodinger's HID/HIE/HIP protonation-state names back to standard HIS
    so ChimeraX recognizes the residues and draws continuous ribbon."""
    for atom in st.atom:
        if atom.pdbres.strip() in ("HID", "HIE", "HIP"):
            atom.pdbres = "HIS "
    return st


def main():
    a = parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    picked = []
    for csv in a.candidates:
        df = pd.read_csv(csv)
        if a.sort_by:
            if a.sort_by not in df.columns:
                raise SystemExit(f"{a.sort_by} not in {csv}")
            df = df.sort_values(a.sort_by, ascending=False)
        picked.append(df.head(a.n))
    to_align = pd.concat(picked, ignore_index=True)
    names = to_align[a.name_col].tolist()
    print(f"{len(names)} structures from {len(a.candidates)} file(s), "
          f"top {a.n} each, ranked on {a.sort_by or 'file order'}")

    index = build_cif_index(a.boltz_root)
    ref = structure.StructureReader.read(a.reference)

    mobile, missing = {}, []
    for name in names:
        if name not in index:
            missing.append(name)
            continue
        mobile[name] = structure.StructureReader.read(str(index[name]))
    if missing:
        print(f"WARNING: no Boltz cif for {len(missing)}: {missing[:5]}")

    if a.align_on_chain:
        rc, mc = a.align_on_chain
        rows = []
        for n, st in mobile.items():
            rows.append({"name": n, "rmsd": superpose_on_chain(ref, st, rc, mc),
                         "score": None})
        print(f"superposed on chain {rc} (ref) <- {mc} (mobile)")
    else:
        alignments = structalign2.align_many(ref, list(mobile.values()))
        rows = [{"name": n, "rmsd": al.rmsd, "score": al.score}
                for al, n in zip(alignments, mobile.keys())]
    pd.DataFrame(rows).to_csv(out / "alignment_data.csv", index=False)

    for name, st in mobile.items():
        with structure.StructureWriter(str(out / f"{name}.pdb")) as w:
            w.append(standardize_his_names(st))

    # provenance: without this there is no way to tell which ranking produced a
    # given aligned_structures* directory
    (out / "PROVENANCE.json").write_text(json.dumps({
        "candidates": a.candidates, "sort_by": a.sort_by, "n_per_file": a.n,
        "name_col": a.name_col, "boltz_root": a.boltz_root,
        "reference": a.reference, "n_aligned": len(mobile), "missing": missing,
        "align_on_chain": a.align_on_chain,
    }, indent=2))
    print(f"wrote {len(mobile)} structures + alignment_data.csv + PROVENANCE.json to {out}")


if __name__ == "__main__":
    main()
