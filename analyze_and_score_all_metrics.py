import pandas as pd
import os

from kix_scoring import (
    CONFIDENCE_FILTERS, PHYSICAL_FILTERS, HELICITY_FILTERS,
    apply_filters, add_both_priority_scores, add_score_variants, SCORE_VARIANTS,
    SCORE_VARIANTS_V2,
)

METRICS_DIR = "/scratch/jem9759/ZhangWork/KIX_Project/full_library_all_metrics"
TSV_OUT = "/scratch/jem9759/ZhangWork/KIX_Project/tsv_outputs/full_library_out.tsv"

def name_from_path(path):
    base = os.path.basename(str(path))
    return base.replace("_model_0.cif", "").replace(".cif", "")

# --- 1. LOAD ALL LAYERS ---
boltz    = pd.read_csv(f"{METRICS_DIR}/boltz_metrics_full_library.csv")
bindcraft= pd.read_csv(f"{METRICS_DIR}/bindcraft_full_library.csv")
# v2 table carries BOTH hit_num (unchanged) and hit_num_v2 (adds the apolar
# contact term). Falls back to the original file if v2 has not been built.
import os as _os
_hitfile = f"{METRICS_DIR}/hbond_hit_counts_v2.csv"
if not _os.path.exists(_hitfile):
    _hitfile = f"{METRICS_DIR}/hbond_hit_counts.csv"
hit       = pd.read_csv(_hitfile, sep="\t")
count     = pd.read_csv(f"{METRICS_DIR}/enrichment_count.csv", sep="\t")
face      = pd.read_csv(f"{METRICS_DIR}/full_library_face_assignment.tsv", sep="\t")

# DSSP: key on 'file' path -> extract name; rename helix column
dssp = pd.read_csv(f"{METRICS_DIR}/dssp_full_library.csv", sep="\t")
dssp["name"] = dssp["file"].apply(name_from_path)
# chain_b_helix_count is carried alongside the fraction: the FRACTION is
# length-dependent (its denominator is peptide length - 2, while the helical
# core stays ~constant), so the raw COUNT is the length-robust QC statistic.
# Measured on the literature controls: fraction vs length rho=-0.824, whereas
# count vs length rho=-0.115 (n.s.).
dssp = dssp[["name", "chain_b_helix_fraction", "chain_b_helix_count"]].rename(
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
# Thresholds live in kix_scoring.py so analyze_controls.py audits the controls
# against these exact clauses.
df_conf = apply_filters(df, CONFIDENCE_FILTERS)
print(f"After confidence filter: {len(df_conf)}")

# --- 5. PHYSICAL + everything FILTER (mentor cell 55) ---
df_filter = apply_filters(df_conf, PHYSICAL_FILTERS).copy()
print(f"After physical filter: {len(df_filter)}")

# DSSP helicity filter: keep chain-B helix fraction > 0.70
df_filter = apply_filters(df_filter, HELICITY_FILTERS)
print(f"After helicity filter: {len(df_filter)}")


# --- 6. DEDUP by sequence (keep highest protein_iptm) ---
df_filter = df_filter.sort_values("protein_iptm", ascending=False).drop_duplicates(subset=["Sequence"], keep="first")
print(f"After dedup: {len(df_filter)}")

# --- 7. priority_score (mentor cell 60) -> see kix_scoring.py ---

# --- 8. SPLIT BY FACE, score each within its own population ---
cmyb = df_filter[df_filter["face_call"] == "cmyb"].copy()
mll  = df_filter[df_filter["face_call"] == "mll"].copy()
print(f"c-Myb survivors: {len(cmyb)}, MLL survivors: {len(mll)}")

cmyb_scored = add_both_priority_scores(cmyb)
mll_scored  = add_both_priority_scores(mll)

# Score variants, computed IN ADDITION to the originals above (nothing removed):
#   *_no_helix -> the helix-fraction term dropped
#   *_pde      -> that term replaced by complex_pde
# See kix_scoring.SCORE_VARIANTS and ControlPlan.md for why.
cmyb_scored = add_score_variants(cmyb_scored).sort_values("priority_score", ascending=False)
mll_scored  = add_score_variants(mll_scored).sort_values("priority_score", ascending=False)

_all_variants = [v for v in list(SCORE_VARIANTS) + list(SCORE_VARIANTS_V2) if v]
variant_cols = [f"priority_score{s}" for s in _all_variants] + \
               [f"priority_score{s}_no_enrichment" for s in _all_variants]

# --- 9. OUTPUT ---
preview_cols = ["name", "Sequence", "count", "hit_num", "hit_num_v2", "helix_score",
                "chain_b_helix_count", "complex_pde",
                "interface_delta_unsat_hbonds", "interface_dG", "protein_iptm",
                "priority_score", "priority_score_no_enrichment"] + variant_cols
preview_cols = [c for c in preview_cols if c in cmyb_scored.columns]
cmyb_scored[preview_cols].to_csv(f"{METRICS_DIR}/cmyb_candidates.csv", index=False)
mll_scored[preview_cols].to_csv(f"{METRICS_DIR}/mll_candidates.csv", index=False)

print("\n=== TOP c-Myb ===")
print(cmyb_scored[preview_cols].head(10).to_string())
print("\n=== TOP MLL ===")
print(mll_scored[preview_cols].head(10).to_string())




