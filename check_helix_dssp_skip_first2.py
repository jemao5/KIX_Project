#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from check_helix_dssp import (
    DEFAULT_INPUT_GLOB,
    _get_chain_b_residues,
    _run_dssp_with_fallback,
    collect_cif_files,
    find_dssp_executable,
)
from Bio.PDB import MMCIFParser


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = SCRIPT_DIR / "chain_b_helix_results_skip_first2.tsv"
HELIX_CODES = {"H", "G", "I"}
DEFAULT_SKIP_LEADING_RESIDUES = 2
RESULT_COLUMNS = [
    "file",
    "passed",
    "dssp_source",
    "dssp_warning",
    "error",
    "chain_b_residue_count",
    "chain_b_considered_residue_count",
    "chain_b_assigned_count",
    "chain_b_helix_count",
    "chain_b_helix_fraction",
    "chain_b_is_helix",
    "chain_b_dssp",
    "chain_b_considered_dssp",
    "chain_b_residue_numbers",
    "chain_b_considered_residue_numbers",
    "skip_leading_residues",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Use Biopython DSSP to evaluate whether chain B forms a helix, "
            "ignoring the first two residues by default."
        )
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
        "--skip-leading-residues",
        type=int,
        default=DEFAULT_SKIP_LEADING_RESIDUES,
        help="How many leading chain B residues to ignore in helix scoring.",
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


def analyze_chain_b_helix_skip_leading(
    cif_path: str,
    threshold: float,
    dssp_executable: str,
    skip_leading_residues: int,
) -> dict:
    parser = MMCIFParser(QUIET=True)
    result = {
        "file": cif_path,
        "passed": False,
        "dssp_source": None,
        "dssp_warning": None,
        "error": None,
        "chain_b_residue_count": np.nan,
        "chain_b_considered_residue_count": np.nan,
        "chain_b_assigned_count": np.nan,
        "chain_b_helix_count": np.nan,
        "chain_b_helix_fraction": np.nan,
        "chain_b_is_helix": False,
        "chain_b_dssp": None,
        "chain_b_considered_dssp": None,
        "chain_b_residue_numbers": None,
        "chain_b_considered_residue_numbers": None,
        "skip_leading_residues": skip_leading_residues,
    }

    try:
        structure = parser.get_structure(Path(cif_path).stem, cif_path)
        model = structure[0]
        chain_b_residues = _get_chain_b_residues(model)
        considered_residues = chain_b_residues[skip_leading_residues:]
        if not considered_residues:
            raise ValueError(
                "No chain B residues remain after skipping leading residues."
            )

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

        considered_residue_numbers = residue_numbers[skip_leading_residues:]
        considered_dssp_codes = dssp_codes[skip_leading_residues:]
        assigned_codes = [code for code in considered_dssp_codes if code != "?"]
        helix_count = sum(code in HELIX_CODES for code in assigned_codes)
        assigned_count = len(assigned_codes)
        helix_fraction = helix_count / assigned_count if assigned_count else np.nan

        result.update(
            {
                "passed": bool(assigned_count and helix_fraction >= threshold),
                "dssp_source": dssp_source,
                "dssp_warning": None if direct_error is None else str(direct_error),
                "chain_b_residue_count": len(chain_b_residues),
                "chain_b_considered_residue_count": len(considered_residues),
                "chain_b_assigned_count": assigned_count,
                "chain_b_helix_count": helix_count,
                "chain_b_helix_fraction": helix_fraction,
                "chain_b_is_helix": bool(assigned_count and helix_fraction >= threshold),
                "chain_b_dssp": "".join(dssp_codes),
                "chain_b_considered_dssp": "".join(considered_dssp_codes),
                "chain_b_residue_numbers": ",".join(
                    str(residue_number) for residue_number in residue_numbers
                ),
                "chain_b_considered_residue_numbers": ",".join(
                    str(residue_number)
                    for residue_number in considered_residue_numbers
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
    skip_leading_residues: int,
) -> pd.DataFrame:
    workers = max(1, min(workers, len(files)))
    results: list[dict] = []

    if workers == 1:
        for index, cif_path in enumerate(files, start=1):
            results.append(
                analyze_chain_b_helix_skip_leading(
                    cif_path,
                    threshold,
                    dssp_executable,
                    skip_leading_residues,
                )
            )
            if index % progress_every == 0 or index == len(files):
                print(f"Processed {index}/{len(files)} files")
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    analyze_chain_b_helix_skip_leading,
                    cif_path,
                    threshold,
                    dssp_executable,
                    skip_leading_residues,
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
    if args.skip_leading_residues < 0:
        raise ValueError("--skip-leading-residues must be >= 0")

    dssp_executable = find_dssp_executable()
    files = collect_cif_files(args.input_glob, args.limit)

    print(f"DSSP executable: {dssp_executable}")
    print(f"Input files: {len(files)}")
    print(f"Workers: {max(1, min(args.workers, len(files)))}")
    print(f"Threshold: {args.threshold:.2f}")
    print(f"Skip leading residues: {args.skip_leading_residues}")

    results_df = run_analysis(
        files=files,
        threshold=args.threshold,
        workers=args.workers,
        progress_every=max(1, args.progress_every),
        dssp_executable=dssp_executable,
        skip_leading_residues=args.skip_leading_residues,
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