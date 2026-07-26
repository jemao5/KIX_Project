#!/usr/bin/env python3
"""
Rename a ResidueX-grafted residue's atoms to Rosetta's canonical names.

Why this is needed: `integrate_NCAA_into_peptide` writes an SDF, which the
example converts to PDB with obabel. That round-trip discards per-residue atom
naming -- the grafted residue comes back with every carbon called `C` (ten of
them for BCS), split across two non-contiguous blocks, as HETATM in chain A.
Duplicate names within a residue is invalid PDB: viewers can't infer bonds
(the side chain renders "floating"), and Rosetta can't match atoms to a
ResidueType by name.

Fix: build the reference bond graph from the residue's Rosetta .params file and
solve the graph isomorphism onto the grafted fragment, which assigns each atom
its canonical name. Anchored on N/C/O, which ResidueX preserves from the
original backbone, so the search space is tiny.

Hydrogens are DROPPED. The Boltz structures are heavy-atom only, so keeping
ResidueX's hydrogens on just the grafted residue would leave the chain
inconsistently protonated. Both prepwizard and Rosetta add hydrogens from the
ResidueType anyway.

Rosetta virtual atoms (V-prefixed, e.g. VCG/VCD1 in C00/C01) are never written;
Rosetta generates them itself.
"""
import argparse
import sys
from pathlib import Path

# Rosetta atom type -> element
TYPE_ELEMENT = {
    "Nbb": "N", "CAbb": "C", "CObb": "C", "OCbb": "O",
    "CH0": "C", "CH1": "C", "CH2": "C", "CH3": "C", "aroC": "C",
    "COO": "C", "OOC": "O", "OH": "O", "ONH2": "O", "Oaro": "O",
    "S": "S", "SH1": "S", "Nlys": "N", "Nhis": "N", "Ntrp": "N",
    "NH2O": "N", "Narg": "N", "Npro": "N",
}
COVALENT = {"C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}
BOND_TOL = 0.45


def parse_args():
    p = argparse.ArgumentParser(description="Canonicalise grafted-residue atom names")
    p.add_argument("--pdb", required=True)
    p.add_argument("--residue-id", type=int, required=True)
    p.add_argument("--params", required=True, help="Rosetta .params for this ncAA")
    p.add_argument("--out", required=True)
    p.add_argument("--chain", default="B")
    return p.parse_args()


def read_params(path):
    """Heavy-atom names, elements and bonds from a Rosetta .params file."""
    names, elements, bonds = [], {}, []
    for line in open(path):
        f = line.split()
        if not f:
            continue
        if f[0] == "ATOM" and len(f) >= 3:
            name, rtype = f[1], f[2]
            if name.startswith("V"):          # Rosetta virtual atom
                continue
            el = TYPE_ELEMENT.get(rtype)
            if el is None or el == "H":
                continue
            names.append(name)
            elements[name] = el
        elif f[0] == "BOND" and len(f) >= 3:
            bonds.append((f[1], f[2]))
    heavy = set(names)
    bonds = [(a, b) for a, b in bonds if a in heavy and b in heavy]
    adj = {n: set() for n in names}
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    return names, elements, adj


def read_residue_atoms(pdb, residue_id, chain):
    """Heavy atoms of the target residue, plus every other line, in file order."""
    target, others = [], []
    for line in open(pdb):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        line = line.rstrip("\n").ljust(80)
        if line[21] == chain and int(line[22:26]) == residue_id:
            el = line[76:78].strip() or line[12:16].strip()[0]
            if el.upper() == "H":
                continue                       # drop hydrogens
            target.append(line)
        else:
            others.append(line)
    return target, others


def element_of(line):
    el = line[76:78].strip()
    return el.capitalize() if el else line[12:16].strip()[0].upper()


def coords(line):
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def infer_bonds(atom_lines):
    """Distance-based bonding using covalent radii."""
    n = len(atom_lines)
    xyz = [coords(l) for l in atom_lines]
    els = [element_of(l) for l in atom_lines]
    adj = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            ri = COVALENT.get(els[i], 0.77)
            rj = COVALENT.get(els[j], 0.77)
            d2 = sum((xyz[i][k] - xyz[j][k]) ** 2 for k in range(3))
            if d2 <= (ri + rj + BOND_TOL) ** 2:
                adj[i].add(j)
                adj[j].add(i)
    return adj, els


def match_graph(ref_names, ref_el, ref_adj, frag_el, frag_adj, seed):
    """Backtracking isomorphism: reference name -> fragment index.

    `seed` pre-assigns the backbone atoms ResidueX preserved (N/C/O), which
    collapses the search almost immediately. Symmetric groups (phenyl ring,
    carboxylate O's) admit several equivalent solutions; any is chemically
    identical, so the first is taken.
    """
    order = [n for n in ref_names if n not in seed]
    # match the most-constrained atoms first
    order.sort(key=lambda n: -len(ref_adj[n]))
    assign = dict(seed)
    used = set(seed.values())

    def ok(name, idx):
        if frag_el[idx] != ref_el[name]:
            return False
        if len(frag_adj[idx]) != len(ref_adj[name]):
            return False
        for nb in ref_adj[name]:
            if nb in assign and assign[nb] not in frag_adj[idx]:
                return False
        return True

    def bt(k):
        if k == len(order):
            return True
        name = order[k]
        # prefer candidates adjacent to something already assigned
        cands = [i for i in range(len(frag_el)) if i not in used]
        anchored = [nb for nb in ref_adj[name] if nb in assign]
        if anchored:
            cands = [i for i in cands if any(assign[nb] in frag_adj[i] for nb in anchored)]
        for idx in cands:
            if ok(name, idx):
                assign[name] = idx
                used.add(idx)
                if bt(k + 1):
                    return True
                used.discard(idx)
                del assign[name]
        return False

    return assign if bt(0) else None


def canonicalise_residue(pdb, residue_id, params, out, chain="B", resname=None,
                         verbose=True):
    """Rename+reorder one grafted residue. Returns the canonical atom names.

    Raises ValueError on any mismatch rather than exiting, so callers running a
    batch can record the failure and continue.
    """
    ref_names, ref_el, ref_adj = read_params(params)
    target, others = read_residue_atoms(pdb, residue_id, chain)
    if not target:
        raise ValueError(f"No heavy atoms at residue {residue_id} chain {chain}")

    frag_adj, frag_el = infer_bonds(target)
    if len(target) != len(ref_names):
        raise ValueError(
            f"Atom-count mismatch at residue {residue_id}: PDB has {len(target)} "
            f"heavy atoms, params {Path(params).stem} expects {len(ref_names)}")

    seed = {}
    for i, line in enumerate(target):
        nm = line[12:16].strip()
        if nm in ("N", "C", "O") and nm in ref_names:
            if nm in seed:
                seed = {}
                break
            seed[nm] = i

    assign = match_graph(ref_names, ref_el, ref_adj, frag_el, frag_adj, seed)
    if assign is None and seed:
        assign = match_graph(ref_names, ref_el, ref_adj, frag_el, frag_adj, {})
    if assign is None:
        raise ValueError(f"Graph isomorphism failed mapping residue {residue_id} "
                         f"onto {Path(params).stem}")

    rn = (resname or Path(params).stem[:3]).upper()
    renamed = []
    for name in ref_names:
        line = target[assign[name]]
        nm = name if len(name) >= 4 else f" {name:<3}"
        line = "ATOM  " + line[6:12] + nm[:4] + line[16] + rn.ljust(3) + line[20:]
        line = line[:21] + chain + line[22:]
        renamed.append(line)

    out_lines, inserted = [], False
    for line in others:
        if (not inserted and line[21] == chain
                and int(line[22:26]) > residue_id):
            out_lines.extend(renamed)
            inserted = True
        out_lines.append(line)
    if not inserted:
        out_lines.extend(renamed)

    with open(out, "w") as f:
        for i, line in enumerate(out_lines, start=1):
            f.write(f"{line[:6]}{i:5d}{line[11:].rstrip()}\n")
        f.write("END\n")

    if verbose:
        print(f"    renamed residue {residue_id} -> {rn}: "
              f"{len(renamed)} heavy atoms, contiguous")
    return ref_names


# --- params lookup ---------------------------------------------------------
ROSETTA_L_NCAA = ("/scratch/jem9759/envs/BindCraft/lib/python3.10/site-packages/"
                  "pyrosetta/database/chemical/residue_type_sets/fa_standard/"
                  "residue_types/l-ncaa")
CUSTOM_PARAMS = "/scratch/jem9759/ZhangWork/KIX_Project/control_run/rosetta_params"


def find_params(code, extra_dirs=()):
    """Locate the .params whose IO_STRING is `code`. Custom dir wins."""
    for d in (CUSTOM_PARAMS, *extra_dirs, ROSETTA_L_NCAA):
        p = Path(d)
        if not p.is_dir():
            continue
        direct = p / f"{code}.params"
        if direct.exists():
            return str(direct)
        for f in sorted(p.glob("*.params")):
            for line in open(f):
                if line.startswith("IO_STRING"):
                    if line.split()[1] == code:
                        return str(f)
                    break
    raise FileNotFoundError(f"No .params found with IO_STRING {code}")


def main():
    args = parse_args()
    ref_names, ref_el, ref_adj = read_params(args.params)
    target, others = read_residue_atoms(args.pdb, args.residue_id, args.chain)
    if not target:
        sys.exit(f"No heavy atoms at residue {args.residue_id} chain {args.chain}")

    frag_adj, frag_el = infer_bonds(target)
    if len(target) != len(ref_names):
        sys.exit(f"Atom-count mismatch: PDB has {len(target)} heavy atoms, "
                 f"params {Path(args.params).stem} expects {len(ref_names)} "
                 f"({sorted(ref_names)})")

    # Anchor on backbone atoms ResidueX preserved with correct names.
    seed = {}
    for i, line in enumerate(target):
        nm = line[12:16].strip()
        if nm in ("N", "C", "O") and nm in ref_names:
            if nm in seed:
                seed = {}      # ambiguous, fall back to unseeded search
                break
            seed[nm] = i

    assign = match_graph(ref_names, ref_el, ref_adj, frag_el, frag_adj, seed)
    if assign is None and seed:
        assign = match_graph(ref_names, ref_el, ref_adj, frag_el, frag_adj, {})
    if assign is None:
        sys.exit(f"Could not map residue {args.residue_id} onto "
                 f"{Path(args.params).stem}: graph isomorphism failed")

    resname = Path(args.params).stem[:3].upper()
    idx_to_name = {v: k for k, v in assign.items()}
    renamed = []
    for name in ref_names:                       # params order -> contiguous block
        line = target[assign[name]]
        nm = name if len(name) >= 4 else f" {name:<3}"
        line = "ATOM  " + line[6:12] + nm[:4] + line[16] + resname.ljust(3) + line[20:]
        line = line[:21] + args.chain + line[22:]
        renamed.append(line)

    # Splice the renamed block in at the right sequence position.
    out, inserted = [], False
    for line in others:
        if (not inserted and line[21] == args.chain
                and int(line[22:26]) > args.residue_id):
            out.extend(renamed)
            inserted = True
        out.append(line)
    if not inserted:
        out.extend(renamed)

    with open(args.out, "w") as f:
        for i, line in enumerate(out, start=1):
            f.write(f"{line[:6]}{i:5d}{line[11:].rstrip()}\n")
        f.write("END\n")

    print(f"  renamed residue {args.residue_id} -> {resname}: "
          f"{len(renamed)} heavy atoms, contiguous")
    print(f"  names: {' '.join(ref_names)}")


if __name__ == "__main__":
    main()
