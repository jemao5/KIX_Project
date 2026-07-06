#!/usr/bin/env python3
"""
Stage 2 of the helix-building pipeline.

Reads a TSV of "name<TAB>sequence" lines (produced by extract_sequences.py)
and builds an ideal alpha-helical PDB for each peptide, forcing chain id "B".

This MUST be run with Schrodinger's bundled Python, e.g.:

    /share/apps/images/run-schrodinger-2025.4.bash run python3 \
        build_helices.py sequences.tsv

Based on the original single-peptide script; the build_helix and
write_pdb_with_chain_id logic is unchanged, just wrapped in a loop.
"""
import argparse
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
# PDB_DIR = SCRIPT_DIR / "pdb_sts"
PDB_DIR = Path("/scratch/jem9759/pdb_sts_full_library") 
CHAIN_ID = "B"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build ideal alpha-helical peptides for every entry in a TSV."
    )
    parser.add_argument(
        "tsv_path",
        help="Path to TSV with one 'name<TAB>sequence' per line",
    )
    return parser.parse_args()


def build_helix(peptide_seq, title):
    from schrodinger.protein import buildpeptide
    
    helix = buildpeptide.build_peptide(
        peptide_seq,
        secondary_structure=buildpeptide.SECONDARY_STRUCTURE.AlphaHelix,
        cap=False,
    )
    helix.title = title
    return helix


def write_pdb_with_chain_id(st, out_pdb):
    PDB_DIR.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_file:
        tmp_pdb = Path(tmp_file.name)

    try:
        st.write(str(tmp_pdb), format="pdb")

        with tmp_pdb.open("r") as fin, out_pdb.open("w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                    line = line[:21] + CHAIN_ID + line[22:]
                fout.write(line)
    finally:
        if tmp_pdb.exists():
            tmp_pdb.unlink()


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


def main():
    args = parse_args()
    entries = read_entries(args.tsv_path)
    print(f"Found {len(entries)} peptides to build")

    n_ok = 0
    n_fail = 0
    for name, seq in entries:
        out_pdb = PDB_DIR / f"{name}.pdb"
        try:
            out_pdb = PDB_DIR / f"{name}.pdb"
            if out_pdb.exists():
                n_ok += 1
                continue
            helix = build_helix(seq, name)
            write_pdb_with_chain_id(helix, out_pdb)
            print(f"  [OK]   {name}: {seq} -> {out_pdb.name}")
            n_ok += 1
        except Exception as exc:  # keep going even if one peptide fails
            print(f"  [FAIL] {name}: {seq} ({exc})")
            n_fail += 1

    print(f"\nDone. {n_ok} built, {n_fail} failed. PDBs in {PDB_DIR}")


if __name__ == "__main__":
    main()