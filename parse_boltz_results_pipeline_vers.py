import pandas as pd
from pathlib import Path
import glob
import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="takes in glob for boltz result JSON files and compiles the necessary metrics into a csv")
    parser.add_argument("input_glob", help="input glob path for JSON files")
    parser.add_argument("out_path", help="output path for combined csv")

    return parser.parse_args()





def main():
    args = parse_args()

    json_paths = glob.glob(args.input_glob)
    results = []



    for jpath in json_paths:
        with open(jpath, 'r') as f:
            data = json.load(f)
        name = Path(jpath).stem
        name = name.replace("confidence_", "")
        name = name.replace("_model_0", "")
        protein_iptm = data["protein_iptm"]
        confidence_score = data["confidence_score"]
        complex_plddt = data["complex_plddt"]
        complex_ipde = data["complex_ipde"]
        complex_pde = data["complex_pde"]
        complex_iplddt = data["complex_iplddt"]
        ptm = data["ptm"]
        pep_ptm = data["chains_ptm"]["1"]
        results.append([name, confidence_score, ptm, protein_iptm, complex_plddt, complex_iplddt, complex_ipde, complex_pde, pep_ptm])
    results_df = pd.DataFrame(results, columns=["name", "confidence_score", "ptm", "protein_iptm", "complex_plddt", "complex_iplddt", "complex_ipde", "complex_pde", "pep_ptm"])

    results_df.to_csv(args.out_path, index=False)


if __name__ == "__main__":
    main()
