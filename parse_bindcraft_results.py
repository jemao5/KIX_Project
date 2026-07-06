import pandas as pd
from pathlib import Path
import glob
import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="takes in glob for bindcraft result JSON files and compiles the necessary metrics into a csv")
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
        name = name.replace("_model_0", "")
        binder_score = data["average"]["scores"]["binder_score"]
        surface_hydrophobicity = data["average"]["scores"]["surface_hydrophobicity"]
        interface_sc = data["average"]["scores"]["interface_sc"]
        interface_packstat = data["average"]["scores"]["interface_packstat"]
        interface_dG = data["average"]["scores"]["interface_dG"]
        interface_dSASA = data["average"]["scores"]["interface_dSASA"]
        interface_dG_SASA_ratio = data["average"]["scores"]["interface_dG_SASA_ratio"]
        interface_fraction = data["average"]["scores"]["interface_fraction"]
        interface_hydrophobicity = data["average"]["scores"]["interface_hydrophobicity"]
        interface_nres = data["average"]["scores"]["interface_nres"]
        interface_interface_hbonds = data["average"]["scores"]["interface_interface_hbonds"]
        interface_hbond_percentage = data["average"]["scores"]["interface_hbond_percentage"]
        interface_delta_unsat_hbonds = data["average"]["scores"]["interface_delta_unsat_hbonds"]
        interface_delta_unsat_hbonds_percentage = data["average"]["scores"]["interface_delta_unsat_hbonds_percentage"]
        results.append([name, binder_score, surface_hydrophobicity, interface_sc, interface_packstat, interface_dG, interface_dSASA, interface_dG_SASA_ratio, interface_fraction, interface_hydrophobicity, interface_nres, interface_interface_hbonds, interface_hbond_percentage, interface_delta_unsat_hbonds, interface_delta_unsat_hbonds_percentage])
    results_df = pd.DataFrame(results, columns=["name", "binder_score", "surface_hydrophobicity", "interface_sc", "interface_packstat", "interface_dG", "interface_dSASA", "interface_dG_SASA_ratio", "interface_fraction", "interface_hydrophobicity", "interface_nres", "interface_interface_hbonds", "interface_hbond_percentage", "interface_delta_unsat_hbonds", "interface_delta_unsat_hbonds_percentage"])


    results_df.to_csv(args.out_path, index=False)


if __name__ == "__main__":
    main()
