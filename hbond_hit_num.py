import argparse
import pickle
import os
import pandas as pd

# --- the two face hit_residues (paste from the resname extraction above) ---
CMYB_HIT_RESIDUES = ['LEU 14 A', 'LEU 18 A', 'LYS 21 A', 'TYR 65 A', 'ALA 69 A', 'ILE 72 A', 'TYR 73 A', 'GLN 76 A']
MLL_HIT_RESIDUES  = ['PHE 27 A', 'ARG 39 A', 'LEU 43 A', 'TYR 46 A', 'LYS 71 A', 'ILE 75 A', 'LEU 79 A']


def parse_args():
    # Defaults reproduce the original hard-coded full-library behaviour; the
    # flags let the control / reference runs point at their own directories.
    p = argparse.ArgumentParser(
        description="Count H-bond + pi-pi contacts to the assigned face's hit residues."
    )
    p.add_argument("--interactions", default="interactions.pkl",
                   help="Pickle from schrodinger_calc_hbond.py")
    p.add_argument("--face", default="final_metric_outputs/full_library_face_assignment.tsv",
                   help="TSV with name + face_call columns")
    p.add_argument("--out", default="interactions_clean.csv",
                   help="Output TSV of name + hit_num")
    return p.parse_args()


def name_from_file(path):
    """Boltz-output basename -> join key.

    Handles both naming conventions in this project:
      Full_Library_Hit_7_model_0_prepared.maegz -> Full_Library_Hit_7
      HBS-22-A_ncaa_prepared.maegz              -> HBS-22-A     (ResidueX run)

    Peeled generically rather than by exact-suffix match: the old exact-match
    list silently failed to strip `_ncaa_prepared.maegz`, so every name missed
    the face table and *every* peptide fell through to hit_num = 0.
    """
    base = os.path.basename(path)
    for ext in (".maegz", ".mae", ".pdb", ".cif"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    for suffix in ("_prepared", "_model_0", "_model_1", "_ncaa"):
        if base.endswith(suffix):
            base = base[: -len(suffix)]
    return base


def main():
    args = parse_args()

    # --- load interactions ---
    with open(args.interactions, "rb") as f:
        all_interactions = pickle.load(f)

    # --- load face assignment: name -> face ('cmyb'/'mll'/'both'/'neither') ---
    face_df = pd.read_csv(args.face, sep="\t")
    name_to_face = dict(zip(face_df["name"], face_df["face_call"]))

    final_data = []
    for d in all_interactions:
        name = name_from_file(d["file"])
        face = name_to_face.get(name)
        # route to the correct face's hit_residues
        if face == "cmyb":
            hit_residues = CMYB_HIT_RESIDUES
        elif face == "mll":
            hit_residues = MLL_HIT_RESIDUES
        else:
            # neither/both shouldn't be in the 30,776 (filtered out), but guard anyway
            final_data.append([name, 0])
            continue
        i = 0
        for hb in d["hbond"]:
            if hb[0] in hit_residues:   # hb[0] = KIX-side residue
                i += 1
        for pp in d["pi_pi"]:
            if pp[0] in hit_residues:
                i += 1
        final_data.append([name, i])

    df = pd.DataFrame(final_data, columns=["name", "hit_num"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.out}; hit_num distribution:")
    print(df["hit_num"].value_counts().sort_index())


if __name__ == "__main__":
    main()
