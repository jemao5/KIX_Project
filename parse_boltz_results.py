#!/usr/bin/env python3
"""
Parse Boltz-2 prediction outputs into a per-peptide metrics table for reranking.

For each peptide it collects:
  - iptm, ptm, complex_plddt, complex_iplddt, complex_ipde  (from confidence JSON)
  - pair_chains_iptm cross term  (KIX <-> peptide interface ipTM)
  - mean_interface_pae           (computed from the pae .npz, cross-chain block)

Boltz output layout (per peptide <name>):
  <results_root>/predictions/<name>/confidence_<name>_model_0.json
  <results_root>/predictions/<name>/pae_<name>_model_0.npz
  <results_root>/predictions/<name>/<name>_model_0.cif

Usage:
    python parse_results.py /scratch/jem9759/boltz_out --out metrics.csv

The results root is the dir that contains a 'predictions/' subfolder (Boltz makes
one 'boltz_results_<name>' per input when run per-file, or a single results dir
when run over a directory -- this script searches recursively for predictions/).
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Collate Boltz-2 metrics per peptide.")
    p.add_argument("results_root", help="Boltz out_dir (searched recursively)")
    p.add_argument("--out", default="metrics.csv", help="Output CSV path")
    p.add_argument(
        "--n_a",
        default=None,
        help="Use if chain A has a constant length",
    )
    return p.parse_args()


def find_prediction_dirs(root):
    """Find every per-peptide prediction folder under any predictions/ dir."""
    root = Path(root)
    pred_dirs = []
    for conf in root.rglob("confidence_*_model_0.json"):
        pred_dirs.append(conf.parent)
    return sorted(set(pred_dirs))


def load_confidence(conf_json):
    with open(conf_json) as f:
        d = json.load(f)
    # pair_chains_iptm is a nested dict keyed by stringified chain indices.
    # The cross term between chain 0 (KIX) and chain 1 (peptide):
    cross_iptm = None
    pci = d.get("pair_chains_iptm")
    if isinstance(pci, dict):
        try:
            cross_iptm = pci["0"]["1"]
        except (KeyError, TypeError):
            cross_iptm = None
    return {
        "iptm": d.get("iptm"),
        "ptm": d.get("ptm"),
        "complex_plddt": d.get("complex_plddt"),
        "complex_iplddt": d.get("complex_iplddt"),
        "complex_pde": d.get("complex_pde"),
        "complex_ipde": d.get("complex_ipde"),
        "confidence_score": d.get("confidence_score"),
        "cross_chain_iptm": cross_iptm,
    }


def load_npz_array(npz_path):
    """Load the single array stored in a Boltz .npz (key name can vary)."""
    with np.load(npz_path) as data:
        keys = list(data.keys())
        # Boltz typically stores under a single key; take the first/largest.
        arr = data[keys[0]]
        if len(keys) > 1:
            # pick the 2D one if multiple
            for k in keys:
                if data[k].ndim == 2:
                    arr = data[k]
                    break
        return arr


def chain_boundary_from_plddt_cif(pred_dir, name):
    """
    Determine where chain A (KIX) ends and chain B (peptide) begins, so we can
    slice the interface block of the PAE matrix.

    Strategy: parse the CIF for the per-residue chain labels (label_asym_id) and
    count how many residues belong to the first chain. Falls back to None.
    """
    cif = pred_dir / f"{name}_model_0.cif"
    if not cif.exists():
        return None
    # Count residues per chain using CA atoms in the CIF atom_site loop.
    # Lightweight parse: find the atom_site column order, then count unique
    # (chain, residue number) for the first chain encountered.
    chain_residues = {}
    order = []
    in_loop = False
    cols = []
    with open(cif) as f:
        for line in f:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".")[1])
                in_loop = True
                continue
            if in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
                parts = s.split()
                if len(parts) < len(cols):
                    continue
                rec = dict(zip(cols, parts))
                # mmCIF uses label_asym_id for chain, label_seq_id for residue
                ch = rec.get("label_asym_id") or rec.get("auth_asym_id")
                res = rec.get("label_seq_id") or rec.get("auth_seq_id")
                if ch is None:
                    continue
                key = (ch, res)
                if ch not in chain_residues:
                    chain_residues[ch] = set()
                    order.append(ch)
                chain_residues[ch].add(res)
            elif in_loop and s.startswith("#"):
                if chain_residues:
                    break
    if not order:
        return None
    first_chain = order[0]
    return len(chain_residues[first_chain])


def mean_interface_pae(pae, n_chain_a):
    """
    Mean PAE over the cross-chain block (chain A rows x chain B cols and vice
    versa). pae is an (N, N) matrix; first n_chain_a indices are chain A.
    """
    if pae is None or n_chain_a is None:
        return None
    n = pae.shape[0]
    if n_chain_a <= 0 or n_chain_a >= n:
        return None
    a = slice(0, n_chain_a)
    b = slice(n_chain_a, n)
    block_ab = pae[a, b]
    block_ba = pae[b, a]
    return float((block_ab.mean() + block_ba.mean()) / 2.0)


def main():
    args = parse_args()
    pred_dirs = find_prediction_dirs(args.results_root)
    if not pred_dirs:
        raise SystemExit(f"No prediction outputs found under {args.results_root}")
    print(f"Found {len(pred_dirs)} prediction folders")

    rows = []
    for pred_dir in pred_dirs:
        # Derive the peptide name from the confidence file.
        conf = next(pred_dir.glob("confidence_*_model_0.json"))
        name = conf.name[len("confidence_"):-len("_model_0.json")]

        metrics = {"name": name}
        metrics.update(load_confidence(conf))

        # interface PAE
        pae_file = pred_dir / f"pae_{name}_model_0.npz"
        if pae_file.exists():
            try:
                pae = load_npz_array(pae_file)
                n_a = args.n_a
                if n_a:
                    n_a = int(n_a)
                else:
                    n_a = chain_boundary_from_plddt_cif(pred_dir, name)
                metrics["mean_interface_pae"] = mean_interface_pae(pae, n_a)
                metrics["n_chain_a_residues"] = n_a
            except Exception as e:
                print(f"  WARN: PAE parse failed for {name}: {e}")
                metrics["mean_interface_pae"] = None
        else:
            metrics["mean_interface_pae"] = None

        rows.append(metrics)

    # Write CSV
    fieldnames = [
        "name", "iptm", "ptm", "cross_chain_iptm",
        "complex_plddt", "complex_iplddt",
        "complex_pde", "complex_ipde",
        "mean_interface_pae", "confidence_score", "n_chain_a_residues",
    ]
    rows.sort(key=lambda r: (r.get("iptm") is None, -(r.get("iptm") or 0)))
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Wrote {len(rows)} rows to {args.out} (sorted by iptm desc)")


if __name__ == "__main__":
    main()
