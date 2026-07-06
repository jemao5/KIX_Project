from pathlib import Path
import argparse
from Bio.PDB import MMCIFParser, PDBIO
import glob
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

def cif_to_pdb(cif_path, pdb_path=None, structure_id=None):
    """Convert an mmCIF file to a PDB file using Biopython.

    Parameters
    ----------
    cif_path : str or Path
        Path to the input mmCIF or CIF file.
    pdb_path : str or Path, optional
        Destination PDB path. If omitted, uses the input stem with a .pdb suffix.
    structure_id : str, optional
        Structure identifier passed to the parser. Defaults to the input file stem.

    Returns
    -------
    Path
        Path to the written PDB file.
    """

    try:
        cif_path = Path(cif_path).expanduser().resolve()
        if not cif_path.exists():
            raise FileNotFoundError(f"CIF file not found: {cif_path}")

        if pdb_path is None:
            pdb_path = cif_path.with_suffix(".pdb")
        else:
            pdb_path = Path(pdb_path)
            os.makedirs(pdb_path, exist_ok=True)
            pdb_path = pdb_path / f"{cif_path.stem}.pdb"


        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(structure_id or cif_path.stem, str(cif_path))

        io = PDBIO()
        io.set_structure(structure)
        io.save(str(pdb_path))
        return pdb_path

    except Exception as e:
        print(f"{cif_path} failed conversion with exception: {e} \n")
        return None


def parse_args():
    parser = argparse.ArgumentParser(description = "Takes a glob of pdb files and creates a cif file version in the same directory or in a defined directory")
    parser.add_argument("cif_glob", help="glob for input cif files")
    parser.add_argument("--pdb_out", default = None, help="option pdb out directory")
    parser.add_argument("--workers", default = 1, type=int, help="number of workers for parallel processing")

    return parser.parse_args()

def main():
    args = parse_args()
    cif_paths = glob.glob(args.cif_glob)
    workers = args.workers


    if workers == 1:
        for cpath in cif_paths:
            cif_to_pdb(cpath, args.pdb_out)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(cif_to_pdb, cpath, args.pdb_out): cpath
                for cpath in cif_paths
            }
        for future in as_completed(futures):
                future.result()


if __name__ == "__main__":
    main()


