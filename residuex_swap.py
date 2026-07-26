#!/usr/bin/env python3
"""
Swap non-canonical amino acids into a cofolded KIX+peptide complex using ResidueX.

Why the structure of this script:

  * ResidueX MUST see the isolated peptide, never the complex. Three of its
    functions are chain-blind -- worst is `get_pep_ready_carbon_alpha`, which
    indexes the N-th 'NCC(=O)' backbone match in the WHOLE molecule, so with KIX
    present "residue 4" resolves to KIX residue 4. We therefore split chain B
    out, swap, and recombine with the untouched chain A.

  * `NCAA_sdf_generation` emits MANY conformers (23 in the shipped example) and
    ResidueX's own `min_distance` only measures clash against the peptide -- it
    never looks at KIX. Since the whole point of the ncAA is to fill the KIX
    pocket, we score every conformer against KIX ourselves and pick from that.

  * The obabel SDF->PDB output is valid PDB but reannotates the grafted residue
    as `HETATM ... UNK A`. Residue NUMBERING survives, so we only rewrite three
    fields: record type, residue name (-> the real 3-letter code), and chain.

Run under the residuex env with ResidueX on PYTHONPATH:
    PYTHONPATH=/scratch/jem9759/ZhangWork/ResidueX \
    /scratch/jem9759/envs/residuex/bin/python3 residuex_swap.py --name MLL2 ...
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fix_graft_naming import canonicalise_residue, find_params

CLASH_CUTOFF = 2.2    # heavy-atom distance below this = steric clash with KIX
CONTACT_CUTOFF = 4.5  # heavy-atom distance below this = a pocket contact


def parse_args():
    p = argparse.ArgumentParser(description="ResidueX ncAA swap on a cofolded complex")
    p.add_argument("--name", required=True, help="Control name, e.g. MLL2")
    p.add_argument("--complex-pdb", required=True, help="Cofolded KIX+peptide PDB")
    p.add_argument("--smiles-csv", required=True, help="ncaa_smiles.csv")
    p.add_argument("--swaps", required=True,
                   help="Comma-separated pos:ncaa_name, e.g. '4:Bcs' or '4:Bcs,6:2mF'")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--peptide-chain", default="B")
    p.add_argument("--target-chain", default="A")
    p.add_argument("--obabel", default="/scratch/jem9759/envs/residuex/bin/obabel")
    p.add_argument("--keep-all-conformers", action="store_true",
                   help="Keep every conformer complex, not just the winner")
    return p.parse_args()


# --- PDB helpers (plain text; avoids Biopython reformatting the file) --------
def read_atom_lines(path):
    with open(path) as f:
        return [l for l in f if l.startswith(("ATOM", "HETATM"))]


def chain_of(line):
    return line[21]


def resnum_of(line):
    return int(line[22:26])


def resname_of(line):
    return line[17:20].strip()


def is_heavy(line):
    el = line[76:78].strip()
    if el:
        return el != "H"
    return not line[12:16].strip().startswith("H")


def coords_of(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def write_pdb(lines, path):
    with open(path, "w") as f:
        f.writelines(lines)
        f.write("END\n")


def split_chains(complex_pdb, peptide_chain, target_chain, out_dir):
    """Write peptide-only and target-only PDBs. ResidueX only ever sees the peptide."""
    lines = read_atom_lines(complex_pdb)
    pep = [l for l in lines if chain_of(l) == peptide_chain]
    tgt = [l for l in lines if chain_of(l) == target_chain]
    if not pep:
        sys.exit(f"No chain {peptide_chain} atoms in {complex_pdb}")
    if not tgt:
        sys.exit(f"No chain {target_chain} atoms in {complex_pdb}")
    pep_path = out_dir / "peptide_only.pdb"
    tgt_path = out_dir / "target_only.pdb"
    write_pdb(pep, pep_path)
    write_pdb(tgt, tgt_path)
    n_res = len({resnum_of(l) for l in pep})
    print(f"  split: chain {peptide_chain} = {n_res} residues / {len(pep)} atoms; "
          f"chain {target_chain} = {len(tgt)} atoms")
    return pep_path, tgt_path


def reannotate(pdb_path, residue_id, three_letter, chain_id, out_path):
    """obabel writes the graft as `HETATM ... UNK A`. Restore ATOM/name/chain.

    Residue numbering survives the round-trip, so only cols 1-6, 18-20 and 22
    need rewriting. Everything else is passed through byte-for-byte.
    """
    out = []
    n_fixed = 0
    for line in open(pdb_path):
        if line.startswith(("ATOM", "HETATM")):
            line = line.rstrip("\n").ljust(80)
            if int(line[22:26]) == residue_id:
                line = "ATOM  " + line[6:17] + three_letter.ljust(3) + line[20:]
                n_fixed += 1
            line = line[:21] + chain_id + line[22:]
            out.append(line.rstrip() + "\n")
    if n_fixed == 0:
        sys.exit(f"reannotate: no atoms found at residue {residue_id} in {pdb_path}")
    write_pdb(out, out_path)
    return n_fixed


def score_against_target(pep_lines, tgt_lines, residue_id):
    """Clash + contact of the grafted residue, measured on the FINAL geometry.

    Returns (min_dist_KIX, n_KIX_contacts, min_dist_intra_peptide).

    We do NOT use ResidueX's own `min_distance` return value for the clash test.
    That number measures the ncAA's atoms against the target backbone, which
    SUPERIMPOSE by design before the redundant atoms are deleted, so it sits at
    0.6-2.5 A for perfectly good conformers and cannot discriminate. Measured on
    MLL6 res6 it was <2.2 A for all 14 conformers while the true KIX clearance
    ranged 0.51-3.83 A.

    Intra-peptide clash skips residues i-1/i+1: they are covalently bonded to the
    graft (peptide bond C-N is ~1.33 A) and are legitimately close.
    """
    ncaa = np.array([coords_of(l) for l in pep_lines
                     if resnum_of(l) == residue_id and is_heavy(l)])
    tgt = np.array([coords_of(l) for l in tgt_lines if is_heavy(l)])
    if len(ncaa) == 0:
        return None, None, None

    min_kix, n_contact = None, None
    if len(tgt):
        d = np.linalg.norm(ncaa[:, None, :] - tgt[None, :, :], axis=-1)
        min_kix, n_contact = float(d.min()), int((d < CONTACT_CUTOFF).sum())

    rest = np.array([coords_of(l) for l in pep_lines
                     if is_heavy(l) and abs(resnum_of(l) - residue_id) >= 2])
    min_intra = None
    if len(rest):
        di = np.linalg.norm(ncaa[:, None, :] - rest[None, :, :], axis=-1)
        min_intra = float(di.min())
    return min_kix, n_contact, min_intra


def assert_index_alignment(pep_path, residue_id):
    """ResidueX resolves the target via `matches[residue_id-1]` over the whole
    molecule's 'NCC(=O)' backbone matches (residuex.py:183-188), which assumes
    match order == residue order. An already-grafted residue is NOT perceived as
    a backbone match by RDKit, so every residue after it shifts down one index --
    silently grafting onto the wrong residue.

    We avoid this by grafting high residue numbers first (see main), but verify
    it explicitly rather than trusting the ordering to hold.
    """
    from rdkit import Chem
    m = Chem.MolFromPDBFile(str(pep_path), removeHs=False)
    if m is None:
        raise ValueError(f"RDKit could not read {pep_path}")
    matches = m.GetSubstructMatches(Chem.MolFromSmiles("NCC(=O)"))
    if residue_id - 1 >= len(matches):
        raise ValueError(f"residue {residue_id}: only {len(matches)} backbone "
                         f"matches; ResidueX would index out of range")
    info = m.GetAtomWithIdx(matches[residue_id - 1][1]).GetPDBResidueInfo()
    got = info.GetResidueNumber()
    if got != residue_id:
        raise ValueError(
            f"index misalignment: ResidueX would graft onto residue {got} "
            f"when asked for residue {residue_id} "
            f"({len(matches)} backbone matches). Refusing to graft.")
    return got


def one_swap(name, pep_path, residue_id, smiles, three_letter, work, obabel):
    """Run ResidueX for a single position; return list of (conformer, pdb_path, min_dist_pep)."""
    from ResidueX.residuex import (split_pdb_by_residue, NCAA_sdf_generation,
                                   integrate_NCAA_into_peptide)

    stage = work / f"swap_res{residue_id}"
    stage.mkdir(parents=True, exist_ok=True)

    residue_pdb = stage / "residue.pdb"
    rest_pdb = stage / "rest.pdb"
    split_pdb_by_residue(str(pep_path), residue_id, str(residue_pdb), str(rest_pdb))

    sdf_dir = stage / "conformers"
    sdf_dir.mkdir(exist_ok=True)
    NCAA_sdf_generation(smiles, str(residue_pdb), str(sdf_dir))
    sdfs = sorted(glob.glob(str(sdf_dir / "*.sdf")))
    if not sdfs:
        sys.exit(f"{name}: NCAA_sdf_generation produced no conformers for res {residue_id}")
    print(f"  res {residue_id} ({three_letter}): {len(sdfs)} conformers generated")

    out_dir = stage / "integrated"
    out_dir.mkdir(exist_ok=True)
    results = []
    for sdf in sdfs:
        tag = Path(sdf).stem
        pep_ready = stage / f"pep_ready_{tag}.pdb"
        out_sdf = out_dir / f"{tag}.sdf"
        dist_txt = out_dir / f"{tag}_dist.txt"
        try:
            min_pep = integrate_NCAA_into_peptide(
                sdf, residue_id, str(pep_path), str(pep_ready),
                str(out_sdf), str(dist_txt))
        except Exception as e:
            print(f"    conformer {tag}: integrate FAILED ({e})")
            continue
        raw_pdb = out_dir / f"{tag}_raw.pdb"
        r = subprocess.run([obabel, "-isdf", str(out_sdf), "-O", str(raw_pdb)],
                           capture_output=True, text=True)
        if not raw_pdb.exists() or raw_pdb.stat().st_size == 0:
            print(f"    conformer {tag}: obabel FAILED ({r.stderr.strip()[:80]})")
            continue
        fixed_pdb = out_dir / f"{tag}_fixed.pdb"
        reannotate(raw_pdb, residue_id, three_letter, "B", fixed_pdb)
        results.append((tag, fixed_pdb, min_pep))
    return results


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    work = out_dir / "work" / args.name
    work.mkdir(parents=True, exist_ok=True)

    smiles_df = pd.read_csv(args.smiles_csv).set_index("ncaa_name")
    swaps = []
    for item in args.swaps.split(","):
        pos, ncaa = item.split(":")
        row = smiles_df.loc[ncaa.strip()]
        swaps.append((int(pos), ncaa.strip(), row["smiles"], row["three_letter_code"],
                      row["parent_aa"]))

    print(f"[{args.name}] {len(swaps)} swap(s): "
          + ", ".join(f"{p}->{n}({c})" for p, n, _, c, _ in swaps))

    pep_path, tgt_path = split_chains(Path(args.complex_pdb), args.peptide_chain,
                                      args.target_chain, work)
    tgt_lines = read_atom_lines(tgt_path)

    # Assert the parent residue is what control_list.csv claims, before grafting.
    pep_lines = read_atom_lines(pep_path)
    three_to_one = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
                    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
                    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    for pos, ncaa, _, _, parent_aa in swaps:
        at = {resname_of(l) for l in pep_lines if resnum_of(l) == pos}
        if not at:
            sys.exit(f"{args.name}: no residue {pos} in chain {args.peptide_chain}")
        got = three_to_one.get(at.pop(), "?")
        if got != parent_aa:
            sys.exit(f"{args.name}: residue {pos} is {got}, expected {parent_aa} "
                     f"for {ncaa} -- numbering mismatch, refusing to graft")
    print(f"  parent-residue check passed for all {len(swaps)} position(s)")

    # Sequential swaps, HIGHEST residue number first. A grafted residue stops
    # being perceived as an 'NCC(=O)' backbone match, shifting the index of every
    # residue after it; descending order keeps already-grafted residues above the
    # next target so its index is unaffected. assert_index_alignment verifies.
    swaps.sort(key=lambda s: -s[0])
    if len(swaps) > 1:
        print(f"  swap order (descending): {[s[0] for s in swaps]}")
    current_pep = pep_path
    done = []
    for pos, ncaa, smiles, code, _ in swaps:
        assert_index_alignment(current_pep, pos)
        results = one_swap(args.name, current_pep, pos, smiles, code, work, args.obabel)
        if not results:
            sys.exit(f"{args.name}: no usable conformer at residue {pos}")
        rows = []
        for tag, fixed_pdb, min_pep_residuex in results:
            lines = read_atom_lines(fixed_pdb)
            min_tgt, n_contact, min_intra = score_against_target(lines, tgt_lines, pos)
            rows.append({"conformer": tag, "pdb": str(fixed_pdb),
                         "min_dist_KIX": None if min_tgt is None else round(min_tgt, 2),
                         "min_dist_intra": None if min_intra is None else round(min_intra, 2),
                         "kix_contacts": n_contact,
                         "residuex_min_dist": round(min_pep_residuex, 2),  # recorded, not used
                         "n_atoms": len(lines)})
        df = pd.DataFrame(rows)
        df["clash_intra"] = df["min_dist_intra"] < CLASH_CUTOFF
        df["clash_kix"] = df["min_dist_KIX"] < CLASH_CUTOFF
        df["viable"] = ~df["clash_intra"] & ~df["clash_kix"]
        # Among viable conformers prefer the most pocket contact. Among
        # non-viable ones prefer the LEAST bad clash -- sorting by contacts there
        # would pick the most deeply interpenetrating conformer, since clashing
        # maximises the contact count.
        df = df.sort_values(
            ["viable", "kix_contacts" if df["viable"].any() else "min_dist_KIX",
             "min_dist_KIX"], ascending=[False, False, False])
        df.to_csv(work / f"conformer_scores_res{pos}.csv", index=False)

        n_ok = int(df["viable"].sum())
        print(f"  res {pos}: {len(df)} conformers, {n_ok} viable "
              f"(no clash <{CLASH_CUTOFF}A to peptide or KIX)")
        print(df[["conformer", "min_dist_intra", "min_dist_KIX",
                  "kix_contacts", "viable"]].head(8).to_string(index=False))
        if n_ok == 0:
            print(f"  !! WARNING: no clash-free conformer at res {pos}; taking best available")
        best = df.iloc[0]
        print(f"  -> selected {best['conformer']} "
              f"(KIX min {best['min_dist_KIX']}A, {best['kix_contacts']} contacts)")

        # Canonicalise ALL grafts done so far, not just this one. Each swap's
        # SDF->obabel round-trip reverts every previously-named residue back to
        # HETATM/UNK with generic atom names, so re-applying only the current
        # residue would silently leave earlier grafts unnamed in the final file.
        done.append((pos, code))
        src = Path(best["pdb"])
        for i, (p, c) in enumerate(sorted(done, key=lambda x: -x[0])):
            dst = work / f"named_after{pos}_{i}_{c}.pdb"
            canonicalise_residue(str(src), p, find_params(c), str(dst),
                                 chain=args.peptide_chain, resname=c)
            src = dst
        current_pep = src

    # Rebuild the peptide: grafted residues from the swap output, everything
    # else from the ORIGINAL peptide.
    #
    # ResidueX only edits the target residue, and non-grafted residues keep
    # their coordinates exactly (verified: HBS_III res-11 CA is bit-identical
    # before and after). But the SDF->obabel round-trip can still mangle their
    # *annotation* -- on HBS_III it relabelled the untouched TYR11 and TRP13 as
    # HETATM/UNK and split TYR11 across two blocks, which cost 2 residues at
    # pose load. Taking untouched residues from the original sidesteps obabel's
    # naming entirely instead of trying to repair it.
    grafted = {p for p, _ in done}
    swapped_lines = read_atom_lines(current_pep)
    orig_lines = read_atom_lines(pep_path)
    by_res = {}
    for l in swapped_lines:
        if resnum_of(l) in grafted:
            by_res.setdefault(resnum_of(l), []).append(l)
    for l in orig_lines:
        if resnum_of(l) not in grafted:
            by_res.setdefault(resnum_of(l), []).append(l)

    final_pep = []
    for rnum in sorted(by_res):
        final_pep.extend(by_res[rnum])

    n_restored = len({resnum_of(l) for l in orig_lines}) - len(grafted)
    print(f"  rebuilt peptide: {len(grafted)} grafted residue(s) + "
          f"{n_restored} restored from original")
    final = tgt_lines + final_pep
    final_path = out_dir / f"{args.name}_ncaa.pdb"
    write_pdb(final, final_path)
    print(f"[{args.name}] wrote {final_path} "
          f"({len(tgt_lines)} target + {len(final_pep)} peptide atoms)")


if __name__ == "__main__":
    main()
