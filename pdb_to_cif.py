#!/usr/bin/env python3
"""
Convert Schrodinger-built helix PDBs to mmCIF WITH a populated sequence.

The raw helix PDBs have no SEQRES, and gemmi's setup_entities() alone does NOT
fill in the polymer sequence -- it leaves _entity_poly_seq empty and the
one-letter code as '?'. Boltz then builds an empty sequence list and throws
IndexError in parse_polymer. The fix is to explicitly assign each polymer
entity's full_sequence from its residue names before writing the CIF.

Run in an env with gemmi installed (e.g. penv):
    python convert_helices_to_cif.py pdb_sts pdb_sts_cif
"""
import sys
from pathlib import Path

import gemmi


def convert_one(pdb_path, out_path):
    st = gemmi.read_structure(str(pdb_path))
    st.setup_entities()

    model = st[0]
    for ent in st.entities:
        if ent.entity_type == gemmi.EntityType.Polymer:
            seq = []
            for sub in ent.subchains:
                poly = model.get_subchain(sub)
                for res in poly:
                    seq.append(res.name)
            if seq:
                ent.full_sequence = seq

    doc = st.make_mmcif_document()
    doc.write_file(str(out_path))


def main():
    in_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "pdb_sts")
    out_dir = Path(sys.argv[2] if len(sys.argv) > 2 else "pdb_sts_cif")
    out_dir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(in_dir.glob("*.pdb"))
    if not pdbs:
        sys.exit(f"No PDBs found in {in_dir}")

    n_ok = 0
    for pdb in pdbs:
        try:
            convert_one(pdb, out_dir / f"{pdb.stem}.cif")
            n_ok += 1
        except Exception as e:
            print(f"  FAIL {pdb.name}: {e}")
    print(f"Converted {n_ok}/{len(pdbs)} -> {out_dir}")


if __name__ == "__main__":
    main()
