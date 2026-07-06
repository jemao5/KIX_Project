import pickle
import os
import pandas as pd

# --- the two face hit_residues (paste from the resname extraction above) ---
CMYB_HIT_RESIDUES = ['LEU 14 A', 'LEU 18 A', 'LYS 21 A', 'TYR 65 A', 'ALA 69 A', 'ILE 72 A', 'TYR 73 A', 'GLN 76 A']
MLL_HIT_RESIDUES  = ['PHE 27 A', 'ARG 39 A', 'LEU 43 A', 'TYR 46 A', 'LYS 71 A', 'ILE 75 A', 'LEU 79 A']

# --- load interactions ---
with open("interactions.pkl", "rb") as f:
    all_interactions = pickle.load(f)

# --- load face assignment: name -> face ('cmyb'/'mll'/'both'/'neither') ---
face_df = pd.read_csv("final_metric_outputs/full_library_face_assignment.tsv", sep="\t")
# adjust column names to match your TSV (name column + face_call column)
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
df.to_csv("interactions_clean.csv", sep="\t", index=False)
print(f"Wrote {len(df)} rows; hit_num distribution:")
print(df["hit_num"].value_counts().sort_index())
