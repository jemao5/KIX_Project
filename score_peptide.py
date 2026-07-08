#!/usr/bin/env python3
import argparse
import json
import os
import sys
import tempfile

import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBIO, PDBParser, MMCIFParser, Selection

import pyrosetta as pr
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

from functions.generic_utils import check_filters, clean_pdb
import glob

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

SCORE_LABEL_MAP = {
    "Binder_Energy_Score": "binder_score",
    "Surface_Hydrophobicity": "surface_hydrophobicity",
    "ShapeComplementarity": "interface_sc",
    "PackStat": "interface_packstat",
    "dG": "interface_dG",
    "dSASA": "interface_dSASA",
    "dG/dSASA": "interface_dG_SASA_ratio",
    "Interface_SASA_%": "interface_fraction",
    "Interface_Hydrophobicity": "interface_hydrophobicity",
    "n_InterfaceResidues": "interface_nres",
    "n_InterfaceHbonds": "interface_interface_hbonds",
    "InterfaceHbondsPercentage": "interface_hbond_percentage",
    "n_InterfaceUnsatHbonds": "interface_delta_unsat_hbonds",
    "InterfaceUnsatHbondsPercentage": "interface_delta_unsat_hbonds_percentage",
}


def init_pyrosetta(dalphaball_path):
    options = [
        "-ignore_unrecognized_res",
        "-ignore_zero_occupancy",
        "-mute all",
        "-corrections::beta_nov16 true",
        "-relax:default_repeats 1",
    ]
    if dalphaball_path:
        options.append(f"-holes:dalphaball {dalphaball_path}")
    pr.init(" ".join(options))


def parse_structure(structure_path):
    lower_path = structure_path.lower()
    if lower_path.endswith(".cif") or lower_path.endswith(".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("complex", structure_path)


def parse_chain_ids(chain_value):
    if not chain_value:
        return []
    return [chain.strip() for chain in chain_value.split(",") if chain.strip()]


def merge_target_chains(structure, target_chain_ids):
    if not target_chain_ids:
        return None, False

    keep_chain_id = target_chain_ids[0]
    if len(target_chain_ids) == 1:
        return keep_chain_id, False

    for model in structure:
        try:
            keep_chain = model[keep_chain_id]
        except KeyError as exc:
            raise ValueError(f"Target chain {keep_chain_id} not found in {structure.id}") from exc

        max_resseq = max((residue.id[1] for residue in keep_chain), default=0)
        next_resseq = max_resseq

        for chain_id in target_chain_ids[1:]:
            try:
                extra_chain = model[chain_id]
            except KeyError as exc:
                raise ValueError(f"Target chain {chain_id} not found in {structure.id}") from exc

            residues = list(extra_chain)
            for residue in residues:
                extra_chain.detach_child(residue.id)
                next_resseq += 1
                residue.id = (residue.id[0], next_resseq, " ")
                keep_chain.add(residue)

            model.detach_child(extra_chain.id)

    return keep_chain_id, True


#def sanitize_structure(structure_path, target_chain_ids=None):
#    structure = parse_structure(structure_path)
#
#    for model in structure:
#        for chain in model:
#            residues = list(chain)
#            for residue in residues:
#                if residue.get_resname().strip().upper() == "ACE":
#                    chain.detach_child(residue.id)
#
#    merged_chain = None
#    merged = False
#    if target_chain_ids:
#        merged_chain, merged = merge_target_chains(structure, target_chain_ids)
#
#    tmp_handle = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
#    tmp_handle.close()
#
#    io = PDBIO()
#    io.set_structure(structure)
#    io.save(tmp_handle.name)
#
#    return tmp_handle.name, merged_chain, merged

def sanitize_structure(structure_path, target_chain_ids=None, remove_chain_ids=("C", "D")):
    structure = parse_structure(structure_path)

    remove_set = {str(c).strip().upper() for c in (remove_chain_ids or [])}

    for model in structure:
        for chain in list(model):
            chain_id = str(chain.id).strip().upper()
            if chain_id in remove_set:
                model.detach_child(chain.id)

        for chain in model:
            for residue in list(chain):
                if residue.get_resname().strip().upper() == "ACE":
                    chain.detach_child(residue.id)

    merged_chain = None
    merged = False
    if target_chain_ids:
        merged_chain, merged = merge_target_chains(structure, target_chain_ids)

    tmp_handle = tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
    tmp_handle.close()

    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp_handle.name)
    return tmp_handle.name, merged_chain, merged


def relax_structure(pdb_path, output_path):
    pose = pr.pose_from_pdb(pdb_path)
    start_pose = pose.clone()

    mmf = MoveMap()
    mmf.set_chi(True)
    mmf.set_bb(True)
    mmf.set_jump(False)

    fastrelax = FastRelax()
    scorefxn = pr.get_fa_scorefxn()
    fastrelax.set_scorefxn(scorefxn)
    fastrelax.set_movemap(mmf)
    fastrelax.max_iter(200)
    fastrelax.min_type("lbfgs_armijo_nonmonotone")
    fastrelax.constrain_relax_to_start_coords(True)
    fastrelax.apply(pose)

    align = AlignChainMover()
    align.source_chain(0)
    align.target_chain(0)
    align.pose(start_pose)
    align.apply(pose)

    for resid in range(1, pose.total_residue() + 1):
        if pose.residue(resid).is_protein():
            bfactor = start_pose.pdb_info().bfactor(resid, 1)
            for atom_id in range(1, pose.residue(resid).natoms() + 1):
                pose.pdb_info().bfactor(resid, atom_id, bfactor)

    pose.dump_pdb(output_path)
    clean_pdb(output_path)
    return output_path


def hotspot_residues(pdb_file, binder_chain, target_chain, atom_distance_cutoff=4.0):
    structure = parse_structure(pdb_file)

    try:
        binder = structure[0][binder_chain]
    except KeyError as exc:
        raise ValueError(f"Binder chain {binder_chain} not found in {pdb_file}") from exc

    try:
        target = structure[0][target_chain]
    except KeyError as exc:
        raise ValueError(f"Target chain {target_chain} not found in {pdb_file}") from exc

    binder_atoms = Selection.unfold_entities(binder, "A")
    binder_coords = np.array([atom.coord for atom in binder_atoms])

    target_atoms = Selection.unfold_entities(target, "A")
    target_coords = np.array([atom.coord for atom in target_atoms])

    binder_tree = cKDTree(binder_coords)
    target_tree = cKDTree(target_coords)

    interacting_residues = {}
    pairs = binder_tree.query_ball_tree(target_tree, atom_distance_cutoff)

    for binder_idx, close_indices in enumerate(pairs):
        binder_residue = binder_atoms[binder_idx].get_parent()
        binder_resname = binder_residue.get_resname()

        if binder_resname in THREE_TO_ONE:
            aa_single_letter = THREE_TO_ONE[binder_resname]
            for close_idx in close_indices:
                _ = target_atoms[close_idx].get_parent()
                interacting_residues[binder_residue.id[1]] = aa_single_letter

    return interacting_residues


def score_interface(pdb_file, binder_chain="B", target_chain="A"):
    pose = pr.pose_from_pdb(pdb_file)

    iam = InterfaceAnalyzerMover()
    iam.set_interface(f"{target_chain}_{binder_chain}")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_dSASA(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.set_pack_separated(True)
    iam.apply(pose)

    interface_AA = {aa: 0 for aa in AMINO_ACIDS}
    interface_residues_set = hotspot_residues(pdb_file, binder_chain, target_chain)
    interface_residues_pdb_ids = []

    for pdb_res_num, aa_type in interface_residues_set.items():
        interface_AA[aa_type] += 1
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    interface_nres = len(interface_residues_pdb_ids)
    interface_residues_pdb_ids_str = ",".join(interface_residues_pdb_ids)

    hydrophobic_aa = set("ACFILMPVWY")
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
    else:
        interface_hydrophobicity = 0

    interfacescore = iam.get_all_data()
    interface_sc = interfacescore.sc_value
    interface_interface_hbonds = interfacescore.interface_hbonds
    interface_dG = iam.get_interface_dG()
    interface_dSASA = iam.get_interface_delta_sasa()
    interface_packstat = iam.get_interface_packstat()
    interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100
    buns_filter = XmlObjects.static_get_filter(
        '<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" '
        'ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" '
        'probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />'
    )
    interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    if interface_nres != 0:
        interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100
        interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100
    else:
        interface_hbond_percentage = None
        interface_bunsch_percentage = None

    chain_design = ChainSelector(binder_chain)
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    tem.set_residue_selector(chain_design)
    binder_score = tem.calculate(pose)

    bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    bsasa.set_residue_selector(chain_design)
    binder_sasa = bsasa.calculate(pose)

    if binder_sasa > 0:
        interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    else:
        interface_binder_fraction = 0

    binder_pose = {
        pose.pdb_info().chain(pose.conformation().chain_begin(i)): p
        for i, p in zip(range(1, pose.num_chains() + 1), pose.split_by_chain())
    }[binder_chain]

    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(binder_pose)

    exp_apol_count = 0
    total_count = 0

    for i in range(1, len(surface_res) + 1):
        if surface_res[i]:
            res = binder_pose.residue(i)
            if res.is_apolar() or res.name() in ("PHE", "TRP", "TYR"):
                exp_apol_count += 1
            total_count += 1

    surface_hydrophobicity = exp_apol_count / total_count if total_count else 0

    interface_scores = {
        "binder_score": binder_score,
        "surface_hydrophobicity": surface_hydrophobicity,
        "interface_sc": interface_sc,
        "interface_packstat": interface_packstat,
        "interface_dG": interface_dG,
        "interface_dSASA": interface_dSASA,
        "interface_dG_SASA_ratio": interface_dG_SASA_ratio,
        "interface_fraction": interface_binder_fraction,
        "interface_hydrophobicity": interface_hydrophobicity,
        "interface_nres": interface_nres,
        "interface_interface_hbonds": interface_interface_hbonds,
        "interface_hbond_percentage": interface_hbond_percentage,
        "interface_delta_unsat_hbonds": interface_delta_unsat_hbonds,
        "interface_delta_unsat_hbonds_percentage": interface_bunsch_percentage,
    }

    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores, interface_AA, interface_residues_pdb_ids_str


def average_scores(score_list):
    if not score_list:
        return {}
    keys = score_list[0].keys()
    averages = {}
    for key in keys:
        values = [scores[key] for scores in score_list if scores.get(key) is not None]
        if values:
            averages[key] = round(sum(values) / len(values), 2)
        else:
            averages[key] = None
    return averages


def average_interface_aas(aa_list):
    if not aa_list:
        return {aa: 0 for aa in AMINO_ACIDS}
    totals = {aa: 0 for aa in AMINO_ACIDS}
    for aa_counts in aa_list:
        for aa in AMINO_ACIDS:
            totals[aa] += aa_counts.get(aa, 0)
    return {aa: round(totals[aa] / len(aa_list), 2) for aa in AMINO_ACIDS}


def build_filter_metrics(scores_by_model, aas_by_model):
    metrics = {}
    avg_scores = average_scores(scores_by_model)
    avg_aas = average_interface_aas(aas_by_model)

    for label, key in SCORE_LABEL_MAP.items():
        metrics[f"Average_{label}"] = avg_scores.get(key)

    metrics["Average_InterfaceAAs"] = avg_aas

    for idx, scores in enumerate(scores_by_model, 1):
        for label, key in SCORE_LABEL_MAP.items():
            metrics[f"{idx}_{label}"] = scores.get(key)
        metrics[f"{idx}_InterfaceAAs"] = aas_by_model[idx - 1]

    return metrics


def apply_filters(metrics, filters):
    design_labels = list(metrics.keys())
    mpnn_data = [metrics[label] for label in design_labels]
    result = check_filters(mpnn_data, design_labels, filters)

    if result is True:
        return True, []
    return False, result


def load_filters(filters_path):
    with open(filters_path, "r") as handle:
        return json.load(handle)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Score protein-peptide interfaces with PyRosetta and apply filter thresholds."
    )
    parser.add_argument(
        "--input",
        "--pdb",
        dest="input_path",
        required=True,
        help="Complex structure file (.cif or .pdb). ACE residues are stripped. my edit: use a glob to batch",
    )
    parser.add_argument("--binder-chain", default="B", help="Binder (peptide) chain ID. Default: B")
    parser.add_argument(
        "--target-chain",
        default="A",
        help="Target (protein) chain ID(s). Use comma-separated list to merge into the first chain.",
    )
    parser.add_argument("--filters", default=None, help="Path to filters JSON to evaluate.")
    parser.add_argument("--dalphaball-path", default="", help="Required DAlphaBall path for PyRosetta.")
    parser.add_argument("--no-relax", dest="relax", action="store_false", help="Skip PyRosetta FastRelax.")
    parser.add_argument("--output-json", default=None, help="Optional path to write JSON output.")
    parser.add_argument("--output-relaxed-pdb", default=None, help="Optional path to write relaxed pdb output.")
    parser.set_defaults(relax=True)
    return parser.parse_args()


def main():
    args = parse_args()
    structure_paths = glob.glob(args.input_path)

    if not structure_paths:
        print(f"No structures matched: {args.input_path}", file=sys.stderr)
        return 2

    if not args.dalphaball_path:
        print("Missing --dalphaball-path (required for BuriedUnsatHbonds).", file=sys.stderr)
        return 2


    target_chain_ids = parse_chain_ids(args.target_chain)
    if not target_chain_ids:
        print("Missing --target-chain.", file=sys.stderr)
        return 2
    if args.binder_chain in target_chain_ids:
        print("Binder chain overlaps with target chains; please use distinct chain IDs.", file=sys.stderr)
        return 2
    if args.output_json:
        os.makedirs(args.output_json, exist_ok=True)
    

    init_pyrosetta(args.dalphaball_path)

    for path in structure_paths:
        try:
            sanitized_path, scoring_target_chain, target_chain_merged = sanitize_structure(
                path, target_chain_ids
            )
            relaxed_path = None
            scoring_path = sanitized_path
            relax_output_path = None
            try:
                if args.relax:
                    base_name = os.path.splitext(os.path.basename(path))[0]
                    if args.output_relaxed_pdb:
                        os.makedirs(args.output_relaxed_pdb, exist_ok=True)
                        relax_output_path = os.path.join(args.output_relaxed_pdb, f"{base_name}.relaxed.pdb")
                    else:
                        relax_output_path = os.path.join(os.getcwd(), f"{base_name}.relaxed.pdb")
                    relaxed_path = relax_structure(sanitized_path, relax_output_path)
                    scoring_path = relaxed_path

                scores, interface_aas, interface_residues = score_interface(
                    scoring_path,
                    binder_chain=args.binder_chain,
                    target_chain=scoring_target_chain,
                )
            finally:
                if sanitized_path:
                    try:
                        os.unlink(sanitized_path)
                    except OSError:
                        pass
            per_model = [
                {
                    "pdb": path,
                    "scores": scores,
                    "interface_aas": interface_aas,
                    "interface_residues": interface_residues,
                }
            ]
            scores_by_model = [scores]
            aas_by_model = [interface_aas]

            average = {
                "scores": average_scores(scores_by_model),
                "interface_aas": average_interface_aas(aas_by_model),
            }

            result = {
                "inputs": {
                    "input": path,
                    "binder_chain": args.binder_chain,
                    "target_chain": args.target_chain,
                    "target_chain_scoring": scoring_target_chain,
                    "target_chain_merged": target_chain_merged,
                    "relaxed": args.relax,
                    "relaxed_output": relax_output_path,
                    "ace_stripped": True,
                },
                "per_model": per_model,
                "average": average,
            }

            if args.filters:
                filters = load_filters(args.filters)
                metrics = build_filter_metrics(scores_by_model, aas_by_model)
                passed, unmet = apply_filters(metrics, filters)
                result["filters"] = {
                    "path": args.filters,
                    "passed": passed,
                    "unmet": unmet,
                }

            output = json.dumps(result, indent=2)
            if args.output_json:
                base_name = os.path.splitext(os.path.basename(path))[0]
                out_file = os.path.join(args.output_json, f"{base_name}.json")
                with open(out_file, "w") as handle:
                    handle.write(output + "\n")
            else:
                print(output)
        except Exception as e:
            print(f"FAILED {path}: {e}", file=sys.stderr)
            continue

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
