import argparse
from pathlib import Path
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.MMCIFParser import MMCIFParser
from concurrent.futures import ProcessPoolExecutor, as_completed
import glob

CMYB_FACE = {14, 18, 21, 65, 69, 72, 73, 76}
MLL_FACE  = {27, 39, 43, 46, 71, 75, 79}
CUTOFF = 5.0
MIN_CONTACTS = 2

def parse_args():
    parser = argparse.ArgumentParser(description="Read boltz output cif files to determine KIX face binding")
    parser.add_argument("input_glob", help="input glob for cif files")
    parser.add_argument("output_path", help="output path for output csv")
    parser.add_argument("--limit", default=None, help="Limit number of processed files to a number")
    parser.add_argument("--workers", default=1, type=int, help="number of workers")
    return parser.parse_args()

def collect_cif_files(pattern: str, limit: int | None) -> list[str]:
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No CIF files matched pattern: {pattern}")
    if limit is not None:
        limit = int(limit)
        files = files[:limit]
    return files

def classify_face(cif_path):
    structure = MMCIFParser(QUIET=True).get_structure("struct", cif_path)
    model = structure[0]
    kix = model["A"]
    peptide = model["B"]

    peptide_atoms = [a for a in peptide.get_atoms() if a.element != "H"]
    ns = NeighborSearch(peptide_atoms)

    def count_face_contact(face_residues):
        count = 0
        for res_id in face_residues:
            for atm in kix[res_id]:
                if atm.element == "H":
                    continue
                if ns.search(atm.coord, CUTOFF):
                    count += 1
                    break
        return count

    cmyb_contacts = count_face_contact(CMYB_FACE)
    mll_contacts = count_face_contact(MLL_FACE)

    c = cmyb_contacts >= MIN_CONTACTS
    m = mll_contacts >= MIN_CONTACTS
    if c and m:
        call = "both"
    elif c:
        call = "cmyb"
    elif m:
        call = "mll"
    else:
        call = "neither"
    return (cmyb_contacts, mll_contacts, call)
    
    
def error_catch_wrapper_classify_face(cif_path):
    try:
        return classify_face(cif_path)
    except Exception as e:
        return (None, None, f"ERROR: {e}")

def main():
    args = parse_args()
    files = collect_cif_files(args.input_glob, args.limit)

    workers = args.workers
    output = {}
    if workers == 1:
        for cif_path in files:
            output[cif_path] = error_catch_wrapper_classify_face(cif_path)
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(error_catch_wrapper_classify_face, cif_path): cif_path for cif_path in files
            }

            for future in as_completed(futures):
                output[futures[future]] = future.result()


    with open(args.output_path, 'w') as f:

        f.write("name\tcmyb_contacts\tmll_contacts\tface_call\n")

        for cif_path, (cmyb, mll, call) in output.items():
            name = Path(cif_path).parent.name   # the Full_Library_Hit_N folder
            f.write(f"{name}\t{cmyb}\t{mll}\t{call}\n")

    n_err = sum(1 for v in output.values() if v[2] and str(v[2]).startswith("ERROR"))
    print(f"Wrote {len(output)} rows to {args.output_path} ({n_err} errors)")


if __name__ == "__main__":
    main()
