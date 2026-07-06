#!/usr/bin/env python3

from __future__ import annotations

import argparse
import glob
import os
import shutil
import sys
import tempfile
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import DSSP, MMCIFParser, PDBIO


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_GLOB = str(
    SCRIPT_DIR / "pep_inputs/boltz_results_pep_*/predictions/pep_*/*.cif"
)
DEFAULT_OUTPUT = SCRIPT_DIR / "chain_b_helix_results.tsv"
HELIX_CODES = {"H", "G", "I"}
RESULT_COLUMNS = [
    "file",
    "passed",
    "dssp_source",
    "dssp_warning",
    "error",
    "chain_b_residue_count",
    "chain_b_assigned_count",
    "chain_b_helix_count",
    "chain_b_helix_fraction",
    "chain_b_is_helix",
    "chain_b_dssp",
    "chain_b_residue_numbers",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use Biopython DSSP to evaluate whether chain B forms a helix."
    )
    parser.add_argument(
        "--input-glob",
        default=DEFAULT_INPUT_GLOB,
        help="Glob pattern for input CIF files.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Output TSV path.",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.30,
        help="Minimum helix fraction required for chain B to pass.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Number of parallel worker processes.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit for quick validation runs.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=25,
        help="Print progress every N processed files.",
    )
    return parser.parse_args()


def find_dssp_executable() -> str:
    env_bin = Path(sys.executable).resolve().parent
    for candidate_name in ("mkdssp", "dssp"):
        candidate_path = env_bin / candidate_name
        if candidate_path.exists():
            return str(candidate_path)

    for candidate_name in ("mkdssp", "dssp"):
        executable = shutil.which(candidate_name)
        if executable:
            return executable

    raise FileNotFoundError(
        "Could not find mkdssp or dssp. Install DSSP or add it to PATH."
    )


def collect_cif_files(pattern: str, limit: int | None) -> list[str]:
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No CIF files matched pattern: {pattern}")
    if limit is not None:
        files = files[:limit]
    return files


def _get_chain_b_residues(model) -> list:
    if "B" not in model:
        raise KeyError("Chain B not found in model 0.")

    residues = [residue for residue in model["B"] if residue.id[0] == " "]
    if not residues:
        raise ValueError("Chain B does not contain standard polymer residues.")

    return residues


def _get_standard_chain_residues(model, chain_id: str) -> list:
    if chain_id not in model:
        raise KeyError(f"Chain {chain_id} not found in model 0.")

    residues = [residue for residue in model[chain_id] if residue.id[0] == " "]
    if not residues:
        raise ValueError(f"Chain {chain_id} does not contain standard polymer residues.")

    return residues


def _atom_is_hydrogen(atom) -> bool:
    element = getattr(atom, "element", "") or ""
    if element:
        return element.strip().upper() == "H"
    return atom.get_name().strip().upper().startswith("H")


def _collect_non_hydrogen_coordinates(residues: list) -> np.ndarray:
    coordinates = [
        atom.get_coord()
        for residue in residues
        for atom in residue.get_atoms()
        if not _atom_is_hydrogen(atom)
    ]
    if not coordinates:
        return np.empty((0, 3), dtype=float)
    return np.asarray(coordinates, dtype=float)


def _compute_centroid(coordinates: np.ndarray) -> np.ndarray:
    if coordinates.size == 0:
        raise ValueError("Cannot compute centroid from an empty coordinate set.")
    return coordinates.mean(axis=0)


def analyze_target_patch_chain_b_centroid_distance(
    cif_path: str,
    target_residues: list[int],
    target_chain: str = "A",
    binder_chain: str = "B",
    centroid_distance_cutoff: float | None = None,
) -> dict:
    parser = MMCIFParser(QUIET=True)
    requested_target_residues = sorted(set(target_residues))
    result = {
        "file": cif_path,
        "target_chain": target_chain,
        "binder_chain": binder_chain,
        "target_residue_count_requested": len(requested_target_residues),
        "target_residue_count_found": 0,
        "target_residue_numbers_found": [],
        "missing_target_residue_numbers": requested_target_residues,
        "target_atom_count": 0,
        "binder_atom_count": 0,
        "centroid_distance": np.nan,
        "binds_near_target": None,
        "error": None,
    }

    try:
        structure = parser.get_structure(Path(cif_path).stem, cif_path)
        model = structure[0]

        target_chain_residues = _get_standard_chain_residues(model, target_chain)
        binder_residues = _get_standard_chain_residues(model, binder_chain)

        target_residue_lookup = {residue.id[1]: residue for residue in target_chain_residues}
        matched_target_residue_numbers = [
            residue_number
            for residue_number in requested_target_residues
            if residue_number in target_residue_lookup
        ]
        missing_target_residue_numbers = [
            residue_number
            for residue_number in requested_target_residues
            if residue_number not in target_residue_lookup
        ]
        if not matched_target_residue_numbers:
            result["error"] = (
                f"None of the requested target residues were found in chain {target_chain}."
            )
            return result

        target_coords = _collect_non_hydrogen_coordinates(
            [target_residue_lookup[residue_number] for residue_number in matched_target_residue_numbers]
        )
        binder_coords = _collect_non_hydrogen_coordinates(binder_residues)
        if target_coords.size == 0:
            result["error"] = (
                f"No non-hydrogen atoms found for requested residues in chain {target_chain}."
            )
            return result
        if binder_coords.size == 0:
            result["error"] = f"No non-hydrogen atoms found in chain {binder_chain}."
            return result

        target_centroid = _compute_centroid(target_coords)
        binder_centroid = _compute_centroid(binder_coords)
        centroid_distance = float(np.linalg.norm(target_centroid - binder_centroid))

        result.update(
            {
                "target_residue_count_found": len(matched_target_residue_numbers),
                "target_residue_numbers_found": matched_target_residue_numbers,
                "missing_target_residue_numbers": missing_target_residue_numbers,
                "target_atom_count": int(target_coords.shape[0]),
                "binder_atom_count": int(binder_coords.shape[0]),
                "centroid_distance": centroid_distance,
                "binds_near_target": (
                    None
                    if centroid_distance_cutoff is None
                    else centroid_distance <= centroid_distance_cutoff
                ),
            }
        )
    except Exception as exc:
        result["error"] = str(exc)

    return result


def _call_dssp(model, structure_path: str, dssp_executable: str, file_type: str):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=".*dictionary mmcif_ma.dic.*",
            category=UserWarning,
        )
        try:
            return DSSP(model, structure_path, dssp=dssp_executable, file_type=file_type)
        except TypeError:
            return DSSP(model, structure_path, dssp=dssp_executable)


def _run_dssp_direct(model, cif_path: str, dssp_executable: str):
    return _call_dssp(model, cif_path, dssp_executable, file_type="MMCIF")


def _run_dssp_with_fallback(model, cif_path: str, dssp_executable: str):
    direct_error = None
    try:
        dssp = _run_dssp_direct(model, cif_path, dssp_executable)
        return dssp, "mmcif", direct_error
    except Exception as exc:
        direct_error = exc

    with tempfile.TemporaryDirectory() as temp_dir:
        pdb_path = os.path.join(temp_dir, f"{Path(cif_path).stem}.pdb")
        io = PDBIO()
        io.set_structure(model)
        io.save(pdb_path)
        dssp = _call_dssp(model, pdb_path, dssp_executable, file_type="PDB")

    return dssp, "pdb_fallback", direct_error


def analyze_chain_b_helix(
    cif_path: str,
    threshold: float,
    dssp_executable: str,
) -> dict:
    parser = MMCIFParser(QUIET=True)
    result = {
        "file": cif_path,
        "passed": False,
        "dssp_source": None,
        "dssp_warning": None,
        "error": None,
        "chain_b_residue_count": np.nan,
        "chain_b_assigned_count": np.nan,
        "chain_b_helix_count": np.nan,
        "chain_b_helix_fraction": np.nan,
        "chain_b_is_helix": False,
        "chain_b_dssp": None,
        "chain_b_residue_numbers": None,
    }

    try:
        structure = parser.get_structure(Path(cif_path).stem, cif_path)
        model = structure[0]
        chain_b_residues = _get_chain_b_residues(model)
        dssp, dssp_source, direct_error = _run_dssp_with_fallback(
            model, cif_path, dssp_executable
        )

        residue_numbers = []
        dssp_codes = []
        for residue in chain_b_residues:
            residue_numbers.append(residue.id[1])
            dssp_key = ("B", residue.id)
            if dssp_key in dssp:
                dssp_code = dssp[dssp_key][2].strip() or "-"
            else:
                dssp_code = "?"
            dssp_codes.append(dssp_code)

        assigned_codes = [code for code in dssp_codes if code != "?"]
        helix_count = sum(code in HELIX_CODES for code in assigned_codes)
        assigned_count = len(assigned_codes)
        helix_fraction = helix_count / assigned_count if assigned_count else np.nan

        result.update(
            {
                "passed": bool(assigned_count and helix_fraction >= threshold),
                "dssp_source": dssp_source,
                "dssp_warning": None if direct_error is None else str(direct_error),
                "chain_b_residue_count": len(chain_b_residues),
                "chain_b_assigned_count": assigned_count,
                "chain_b_helix_count": helix_count,
                "chain_b_helix_fraction": helix_fraction,
                "chain_b_is_helix": bool(assigned_count and helix_fraction >= threshold),
                "chain_b_dssp": "".join(dssp_codes),
                "chain_b_residue_numbers": ",".join(
                    str(residue_number) for residue_number in residue_numbers
                ),
            }
        )
    except Exception as exc:
        result["error"] = str(exc)

    return result


def run_analysis(
    files: list[str],
    threshold: float,
    workers: int,
    progress_every: int,
    dssp_executable: str,
) -> pd.DataFrame:
    workers = max(1, min(workers, len(files)))
    results: list[dict] = []

    if workers == 1:
        for index, cif_path in enumerate(files, start=1):
            results.append(analyze_chain_b_helix(cif_path, threshold, dssp_executable))
            if index % progress_every == 0 or index == len(files):
                print(f"Processed {index}/{len(files)} files")
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    analyze_chain_b_helix,
                    cif_path,
                    threshold,
                    dssp_executable,
                ): cif_path
                for cif_path in files
            }

            for index, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                if index % progress_every == 0 or index == len(files):
                    print(f"Processed {index}/{len(files)} files")

    results_df = pd.DataFrame(results)
    if results_df.empty:
        results_df = pd.DataFrame(columns=RESULT_COLUMNS)
    else:
        results_df = results_df[RESULT_COLUMNS].sort_values("file").reset_index(drop=True)
    return results_df


def main() -> None:
    args = parse_args()
    dssp_executable = find_dssp_executable()
    files = collect_cif_files(args.input_glob, args.limit)

    print(f"DSSP executable: {dssp_executable}")
    print(f"Input files: {len(files)}")
    print(f"Workers: {max(1, min(args.workers, len(files)))}")
    print(f"Threshold: {args.threshold:.2f}")

    results_df = run_analysis(
        files=files,
        threshold=args.threshold,
        workers=args.workers,
        progress_every=max(1, args.progress_every),
        dssp_executable=dssp_executable,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_path, sep="\t", index=False)

    print(f"Saved TSV: {output_path}")
    print(
        "Pass summary: "
        f"{int(results_df['chain_b_is_helix'].fillna(False).sum())}/"
        f"{len(results_df)}"
    )


if __name__ == "__main__":
    main()