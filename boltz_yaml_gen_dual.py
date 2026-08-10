#!/usr/bin/env python3
"""
Generate 3-chain Boltz YAMLs: KIX + TWO copies of the same peptide.

Purpose: test whether a single designed sequence can occupy both KIX faces at
once. Chain A is KIX; chains B and C are two copies of the peptide, each
templated on the same idealized helix. No pocket constraints -- the point is to
let Boltz place the copies freely and then ask, per chain, which face it landed
on (see the per-chain face call in dual_face_determination.py).

Everything about chain A and the peptide template config is inherited from
boltz_yaml_gen.py so the two stay consistent.

Run under boltz_env (PyYAML lives only there):
    /scratch/jem9759/envs/boltz_env/bin/python3 boltz_yaml_gen_dual.py \
        dual_cofold/dual_list.tsv --out_path dual_cofold/yamls \
        --helix_dir cif_outputs/pdb_sts_full_library_cif --kix-msa kix_msa.csv
"""
import argparse
import sys
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from boltz_yaml_gen import (KIX_SEQUENCE, KIX_TEMPLATE_PDB, KIX_TEMPLATE_ID,
                            HELIX_TEMPLATE_ID, read_entries)

PEPTIDE_CHAINS = ("B", "C")


def parse_args():
    p = argparse.ArgumentParser(description="KIX + 2 peptide copies, per peptide")
    p.add_argument("tsv_path")
    p.add_argument("--out_path", default=str(SCRIPT_DIR / "dual_cofold" / "yamls"))
    p.add_argument("--helix_dir",
                   default=str(SCRIPT_DIR / "cif_outputs" / "pdb_sts_full_library_cif"))
    p.add_argument("--kix-msa", default=None)
    return p.parse_args()


def build(name, seq, helix_dir, kix_msa):
    kix = {"id": "A", "sequence": KIX_SEQUENCE}
    if kix_msa:
        kix["msa"] = kix_msa

    template = f"{str(helix_dir).rstrip('/')}/{name}.cif"
    sequences = [{"protein": kix}]
    templates = [{"pdb": KIX_TEMPLATE_PDB, "chain_id": "A",
                  "template_id": KIX_TEMPLATE_ID}]
    for ch in PEPTIDE_CHAINS:
        # Two chains carrying the same sequence. Distinct ids keep them as
        # separate copies; both point at the same helix CIF.
        sequences.append({"protein": {"id": ch, "sequence": seq, "msa": "empty"}})
        templates.append({"cif": template, "chain_id": ch,
                          "template_id": HELIX_TEMPLATE_ID})
    return {"version": 1, "sequences": sequences, "templates": templates}


def main():
    a = parse_args()
    out = Path(a.out_path)
    out.mkdir(parents=True, exist_ok=True)
    entries = read_entries(a.tsv_path)
    print(f"{len(entries)} peptides -> KIX + 2 copies each")

    missing = 0
    for name, seq in entries:
        helix = Path(a.helix_dir) / f"{name}.cif"
        if not helix.exists():
            print(f"  WARNING: helix CIF missing for {name}: {helix}")
            missing += 1
        with open(out / f"{name}.yaml", "w") as f:
            yaml.dump(build(name, seq, a.helix_dir, a.kix_msa), f,
                      default_flow_style=False, sort_keys=False)
    if missing:
        raise SystemExit(f"{missing} helix templates missing -- fix the --helix_dir path")
    print(f"Wrote {len(entries)} YAMLs to {out} (0 missing templates)")


if __name__ == "__main__":
    main()
