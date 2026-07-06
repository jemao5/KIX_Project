import pandas as pd
import os

METRICS_DIR = "/scratch/jem9759/ZhangWork/KIX_Project/full_library_all_metrics"
TSV_OUT = "/scratch/jem9759/ZhangWork/KIX_Project/tsv_outputs/full_library_out.tsv"

def name_from_path(path):
    base = os.path.basename(str(path))
    return base.replace("_model_0.cif", "").replace(".cif", "")

# --- 1. LOAD ALL LAYERS ---
boltz    = pd.read_csv(f"{METRICS_DIR}/boltz_metrics_full_library.csv")
bindcraft= pd.read_csv(f"{METRICS_DIR}/bindcraft_full_library.csv")
hit       = pd.read_csv(f"{METRICS_DIR}/hbond_hit_counts.csv", sep="\t")
count     = pd.read_csv(f"{METRICS_DIR}/enrichment_count.csv", sep="\t")
face      = pd.read_csv(f"{METRICS_DIR}/full_library_face_assignment.tsv", sep="\t")

# DSSP: key on 'file' path -> extract name; rename helix column
dssp = pd.read_csv(f"{METRICS_DIR}/dssp_full_library.csv", sep="\t")
dssp["name"] = dssp["file"].apply(name_from_path)
dssp = dssp[["name", "chain_b_helix_fraction"]].rename(
    columns={"chain_b_helix_fraction": "helix_score"}
)

# name -> sequence (for n_K/n_M and dedup)
seqs = pd.read_csv(TSV_OUT, sep="\t", header=None, names=["name", "Sequence"])

# --- 2. MERGE ON name ---
df = seqs
for layer in [boltz, bindcraft, hit, count, face[["name", "face_call"]], dssp]:
    df = df.merge(layer, on="name", how="inner")   # inner: keep only peptides present in all layers

print(f"Merged: {len(df)} peptides (should be ~30,776 cleanly-faced with hit_num)")

# --- 3. DERIVED COLUMNS ---
df["n_K"] = df["Sequence"].map(lambda x: str(x).count("K"))
df["n_M"] = df["Sequence"].map(lambda x: str(x).count("M"))

# --- 4. CONFIDENCE FILTER (mentor cell 45) ---
df_conf = df[
    (df["confidence_score"] > 0.85) &
    (df["pep_ptm"] > 0.80) &
    (df["complex_pde"] < 0.5)
]
print(f"After confidence filter: {len(df_conf)}")

# --- 5. PHYSICAL + everything FILTER (mentor cell 55) ---
df_filter = df_conf[
    (df_conf["binder_score"] < 0.0) &
    (df_conf["interface_dG"] < -25.0) &
    (df_conf["interface_dSASA"] > 1.0) &
    (df_conf["surface_hydrophobicity"] < 1.0) &
    (df_conf["interface_sc"] > 0.5) &
    (df_conf["interface_nres"] > 4.0) &
    (df_conf["interface_interface_hbonds"] > 1) &
    (df_conf["interface_delta_unsat_hbonds"] <= 5) &
    (df_conf["n_K"] <= 3) &
    (df_conf["n_M"] <= 3) &
    (df_conf["hit_num"] >= 1) &
    (df_conf["count"] > 1) &
    (df_conf["protein_iptm"] > 0.85)
].copy()
print(f"After physical filter: {len(df_filter)}")

# --- 6. DEDUP by sequence (keep highest protein_iptm) ---
df_filter = df_filter.sort_values("protein_iptm", ascending=False).drop_duplicates(subset=["Sequence"], keep="first")
print(f"After dedup: {len(df_filter)}")

# --- 7. priority_score (mentor cell 60) ---
def add_priority_score(df, count_col="count", hit_col="hit_num",
                       unsat_col="interface_delta_unsat_hbonds",
                       protein_col="protein_iptm", helix_score_col="helix_score",
                       count_weight=0.2, hit_weight=0.2, unsat_weight=0.1,
                       protein_weight=0.2, helix_weight=0.2, score_col="priority_score"):
    weights = count_weight + hit_weight + unsat_weight + protein_weight + helix_weight
    d = df.copy()
    count_rank   = d[count_col].rank(pct=True, method="average")
    hit_rank     = d[hit_col].rank(pct=True, method="average")
    unsat_rank   = 1.0 - d[unsat_col].rank(pct=True, method="average")
    protein_rank = d[protein_col].rank(pct=True, method="average")
    helix_rank   = d[helix_score_col].rank(pct=True, method="average")
    d[score_col] = (count_weight*count_rank + hit_weight*hit_rank + unsat_weight*unsat_rank
                    + protein_weight*protein_rank + helix_weight*helix_rank) / weights
    return d.sort_values(score_col, ascending=False)

# --- 8. SPLIT BY FACE, score each within its own population ---
cmyb = df_filter[df_filter["face_call"] == "cmyb"].copy()
mll  = df_filter[df_filter["face_call"] == "mll"].copy()
print(f"c-Myb survivors: {len(cmyb)}, MLL survivors: {len(mll)}")

cmyb_scored = add_priority_score(cmyb)
mll_scored  = add_priority_score(mll)

# --- 9. OUTPUT ---
preview_cols = ["name", "Sequence", "count", "hit_num", "helix_score",
                "interface_delta_unsat_hbonds", "interface_dG", "protein_iptm", "priority_score"]
cmyb_scored[preview_cols].to_csv(f"{METRICS_DIR}/cmyb_candidates.csv", index=False)
mll_scored[preview_cols].to_csv(f"{METRICS_DIR}/mll_candidates.csv", index=False)

print("\n=== TOP c-Myb ===")
print(cmyb_scored[preview_cols].head(10).to_string())
print("\n=== TOP MLL ===")
print(mll_scored[preview_cols].head(10).to_string())




