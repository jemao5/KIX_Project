#!/usr/bin/env python3
"""
Generate one Boltz-2 YAML per peptide for KIX cofolding.

KIX (chain A) is constant across all peptides: real MSA via --use_msa_server,
templated to the 2AGH KIX chain. The peptide (chain B) varies: its sequence
comes from the TSV, msa is empty (de novo), and it is templated to its own
Schrodinger-built helix PDB in pdb_sts/.

Run in your pandas/penv env (only needs PyYAML):
    python generate_yamls.py sequences.tsv --out_path ./out_yaml

Then run boltz over the whole directory:
    boltz predict ./out_yaml --use_msa_server --cache /scratch/jem9759/.boltz \
        --out_dir ./out
"""
import yaml
import argparse
from pathlib import Path
import re

SCRIPT_DIR = Path(__file__).resolve().parent

# KIX construct sequence (chain A) -- REPLACE with the exact sequence you used
# in the validated single test if different.
KIX_SEQUENCE = (
    "GVRKGWHEHVTQDLRSHLVHKLVQAIFPTPDPAALKDRRMENLVAYAKKVEGDMYESANSRDEYYHLLAEKIYKIQKELEEKRRSRL"
)

# Path to the KIX template (single model from 2AGH). KIX is chain B in 2AGH,
# which Boltz renames to "B1" for PDB input -> template_id: B1, applied to
# YAML chain A.
KIX_TEMPLATE_PDB = "/scratch/jem9759/ZhangWork/KIX_Project/2agh_model1.pdb"
KIX_TEMPLATE_ID = "B1"

# Directory holding the per-peptide helix templates AS CIF (converted from the
# Schrodinger PDBs with sequence populated via convert_helices_to_cif.py).
# The raw PDBs lack sequence metadata and cause an IndexError in Boltz; the CIF
# with a populated _entity_poly_seq fixes it. Peptide template uses chain_id
# only (B) and lets Boltz auto-match the file chain by sequence.
# THESE ARE THE TEMPLATE PARAMS
HELIX_CIF_DIR = "./pdb_sts_cif"
HELIX_TEMPLATE_ID = "Bxp"

CMYB_MAX_HIT = 35
CMYB_RESIDUES = [14, 18, 21, 65, 69, 72, 73, 76]
MLL_FACE_RESIDUES = [27, 39, 43, 46, 71, 75, 79]

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Boltz-2 YAMLs for KIX + each peptide."
    )
    parser.add_argument("tsv_path", help="Path to input .tsv (name<TAB>sequence)")
    parser.add_argument(
        "--out_path",
        default=str(SCRIPT_DIR / "out_yaml"),
        help="Output directory for the generated YAMLs",
    )
    parser.add_argument(
        "--helix_dir",
        default=HELIX_CIF_DIR,
        help="Directory with per-peptide helix CIFs (default: ./pdb_sts_cif)",
    )
    parser.add_argument(
        "--kix-msa",
        default=None,
        help="Path to a precomputed KIX .a3m. If set, chain A uses it and you "
             "run boltz WITHOUT --use_msa_server (avoids 57 redundant MSA calls).",
    )
    parser.add_argument(
        "--use_constraints",
        action="store_true",
        help="Choose whether or not to use constraints. True or False"
    )
    
    return parser.parse_args()


def read_entries(tsv_path):
    """Yield (name, sequence) tuples from the TSV, skipping blanks/comments."""
    entries = []
    with open(tsv_path, "r") as fin:
        for lineno, line in enumerate(fin, start=1):
            line = line.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"  WARNING: skipping malformed line {lineno}: {line!r}")
                continue
            name, seq = parts[0].strip(), parts[1].strip()
            if not name or not seq:
                print(f"  WARNING: skipping line {lineno} with empty field")
                continue
            entries.append((name, seq))
    return entries


def build_yaml_dict(name, seq, helix_dir, constraints, kix_msa=None):
    """Construct the Boltz YAML structure for one peptide.

    Working template configuration (verified):
      - KIX (PDB): chain_id A + template_id B1  (PDB subchain renaming).
      - Peptide (CIF): chain_id B only -> Boltz auto-matches the file chain by
        sequence. Requires the CIF to carry a populated sequence (see
        convert_helices_to_cif.py).

    If kix_msa is given, chain A references it (and you drop --use_msa_server).
    """

    peptide_template = f"{helix_dir.rstrip('/')}/{name}.cif"
    
    if constraints:
        m = re.search(r"(\d+)", name)
        num = int(m.group(1)) if m else None
        if num:
            residues = CMYB_RESIDUES if num <= CMYB_MAX_HIT else MLL_FACE_RESIDUES
        else:
            residues = None

    kix_protein = {"id": "A", "sequence": KIX_SEQUENCE}
    if kix_msa:
        kix_protein["msa"] = kix_msa  # precomputed -> no server call
    # else: omit msa -> requires --use_msa_server at runtime

    data = {
        "version": 1,
        "sequences": [
            {"protein": kix_protein},
            {
                "protein": {
                    "id": "B",
                    "sequence": seq,
                    "msa": "empty",  # de novo peptide, no homologs
                }
            },
        ],
        "templates": [
            {"pdb": KIX_TEMPLATE_PDB, "chain_id": "A", "template_id": KIX_TEMPLATE_ID},
            {"cif": peptide_template, "chain_id": "B", "template_id": HELIX_TEMPLATE_ID},
        ]
    }

    if constraints and residues:
        data["constraints"] = [
            {"pocket": {
                "binder":"B",
                "contacts":[["A",i] for i in residues],
                "max_distance": 6,
                "force": True
                }
            }
        ]
    


    return data


def main():
    args = parse_args()
    out_dir = Path(args.out_path)
    out_dir.mkdir(parents=True, exist_ok=True)

    entries = read_entries(args.tsv_path)
    print(f"Found {len(entries)} peptides")

    n_written = 0
    for name, seq in entries:
        # Sanity check: the helix CIF for this peptide should exist.
        helix_path = Path(args.helix_dir) / f"{name}.cif"
        if not helix_path.exists():
            print(f"  WARNING: helix CIF missing for {name}: {helix_path} "
                  f"(YAML still written, but boltz will fail on it)")

        data = build_yaml_dict(name, seq, args.helix_dir, args.use_constraints, kix_msa=args.kix_msa)
        out_file = out_dir / f"{name}.yaml"
        with open(out_file, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)
        n_written += 1

    print(f"Wrote {n_written} YAMLs to {out_dir}")


if __name__ == "__main__":
    main()
