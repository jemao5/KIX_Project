import argparse
import os
import pickle
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kix_scoring import HIT_RESIDUE_SETS, hit_residue_strings


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
    p.add_argument("--residue-set", default="original", choices=HIT_RESIDUE_SETS,
                   help="Which per-face hit-residue list to count against. "
                        "'original' is the hand-picked ChimeraX-contact list and is "
                        "the default so existing results stay reproducible; 'crystal' "
                        "is derived from residues that actually H-bond/pi-stack with "
                        "the native peptide in 2AGH; 'union' is both. "
                        "See kix_scoring.HIT_RESIDUE_SETS.")
    p.add_argument("--hydrophobic", default=None,
                   help="TSV from hydrophobic_contacts.py. When given, an extra "
                        "`hit_num_v2` column is emitted = hit_num + the count of "
                        "APOLAR hit residues contacted. Apolar hit residues cannot "
                        "H-bond, so they contribute nothing to hit_num; this adds "
                        "the hydrophobic packing back. `hit_num` itself is unchanged.")
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

    cmyb_hits = hit_residue_strings(args.residue_set, "cmyb")
    mll_hits = hit_residue_strings(args.residue_set, "mll")
    print(f"hit-residue set: {args.residue_set}")
    print(f"  cmyb ({len(cmyb_hits)}): {cmyb_hits}")
    print(f"  mll  ({len(mll_hits)}): {mll_hits}")

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
            hit_residues = cmyb_hits
        elif face == "mll":
            hit_residues = mll_hits
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

    if args.hydrophobic:
        hyd = pd.read_csv(args.hydrophobic, sep="\t")
        faces = df["name"].map(name_to_face)
        hyd_map = hyd.set_index("name")
        cmyb = df["name"].map(hyd_map["cmyb_apolar_contacts"])
        mll = df["name"].map(hyd_map["mll_apolar_contacts"])
        # pick the count for the face this peptide was assigned to
        picked = cmyb.where(faces == "cmyb", mll.where(faces == "mll", 0))
        missing = picked.isna().sum()
        if missing:
            raise SystemExit(
                f"{missing} structures have no hydrophobic-contact row; "
                f"the two files disagree on names (check --name-from).")
        df["hit_num_v2"] = df["hit_num"] + picked.astype(int)
        print(f"hit_num_v2 = hit_num + apolar contacts  "
              f"(mean {df['hit_num'].mean():.2f} -> {df['hit_num_v2'].mean():.2f})")

    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.out}; hit_num distribution:")
    print(df["hit_num"].value_counts().sort_index())


if __name__ == "__main__":
    main()
