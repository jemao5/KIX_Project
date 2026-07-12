import pickle
import os
import pandas as pd

# Reference-run copy of hbond_hit_num.py for the two native controls
# (kix_peptide_reference_cmyb / kix_peptide_reference_mll). Same counting logic
# as the full-library script, but all I/O is scoped to boltz_reference_structures/
# so the full-library files are never touched.

REF_DIR = "boltz_reference_structures"

# --- the two face hit_residues (identical to hbond_hit_num.py) ---
CMYB_HIT_RESIDUES = ['LEU 14 A', 'LEU 18 A', 'LYS 21 A', 'TYR 65 A', 'ALA 69 A', 'ILE 72 A', 'TYR 73 A', 'GLN 76 A']
MLL_HIT_RESIDUES  = ['PHE 27 A', 'ARG 39 A', 'LEU 43 A', 'TYR 46 A', 'LYS 71 A', 'ILE 75 A', 'LEU 79 A']

# --- load interactions ---
with open(os.path.join(REF_DIR, "reference_interactions.pkl"), "rb") as f:
    all_interactions = pickle.load(f)

# --- load face assignment: name -> face ('cmyb'/'mll'/'both'/'neither') ---
face_df = pd.read_csv(os.path.join(REF_DIR, "reference_face_assignment.tsv"), sep="\t")
name_to_face = dict(zip(face_df["name"], face_df["face_call"]))

def name_from_file(path):
    base = os.path.basename(path)
    for suffix in ["_model_0_prepared.maegz", "_model_1_prepared.maegz"]:
        base = base.replace(suffix, "")
    return base

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
out_path = os.path.join(REF_DIR, "reference_interactions_clean.csv")
df.to_csv(out_path, sep="\t", index=False)
print(f"Wrote {len(df)} rows to {out_path}; hit_num distribution:")
print(df["hit_num"].value_counts().sort_index())
print(df.to_string(index=False))
