#!/usr/bin/env python3
"""
Score the literature control peptides (control_list.csv) on the same scale as
the full library.

The controls have no library enrichment `count`, so `priority_score` is
undefined for them. `priority_score_no_enrichment` -- the same percentile blend
with the count term dropped and the remaining weights renormalized -- is what
makes the two populations directly comparable, and is what this script ranks on.

Controls are appended to the library's filter survivors rather than being run
through the filters, so every control gets a score even when it fails a clause.
`filter_audit` records which clauses each control would have failed.

Controls are pooled with the face they are known from the literature to bind
(`face` in control_list.csv). The independent geometric call from stage 9 is
reported as `structural_face_call` so the two can be compared.

Outputs (in control_run/):
  control_scores.csv                 -- all controls, metrics + audit + score + percentile
  cmyb_candidates_with_controls.csv  -- c-Myb library survivors + c-Myb controls
  mll_candidates_with_controls.csv   -- MLL library survivors + MLL controls
"""
import argparse
import os

import pandas as pd

from kix_scoring import (
    CONFIDENCE_FILTERS, PHYSICAL_FILTERS, HELICITY_FILTERS, ALL_FILTERS,
    LIBRARY_DESIGN_FILTERS, apply_filters, filter_audit, add_both_priority_scores,
    add_swap_sensitive_score, add_score_variants, SCORE_VARIANTS,
    SCORE_VARIANTS_V2,
)

ROOT        = "/scratch/jem9759/ZhangWork/KIX_Project"
METRICS_DIR = f"{ROOT}/full_library_all_metrics"
TSV_OUT     = f"{ROOT}/tsv_outputs/full_library_out.tsv"
CTRL_DIR    = f"{ROOT}/control_run"

# The control set exists in two variants: the original canonical-stand-in run,
# and the ResidueX ncAA-grafted re-run. Same code, different input layer files.
VARIANTS = {
    "canonical": dict(boltz="control_boltz_data.csv",
                      bindcraft="control_bindcraft_data.csv",
                      hit="control_interactions_clean.csv",
                      face="control_face_assignment.tsv",
                      dssp="control_dssp_data.csv",
                      tag=""),
    "ncaa":      dict(boltz="control_boltz_data.csv",       # fold is pre-graft, unchanged
                      bindcraft="control_bindcraft_ncaa.csv",
                      hit="control_interactions_clean_ncaa_v2.csv",
                      face="control_face_assignment_ncaa.tsv",
                      dssp="control_dssp_ncaa.csv",
                      tag="_ncaa"),
}


def parse_cli():
    p = argparse.ArgumentParser(description="Score literature controls against the library")
    p.add_argument("--variant", choices=sorted(VARIANTS), default="canonical")
    return p.parse_args()


_args = parse_cli()
V = VARIANTS[_args.variant]
TAG = V["tag"]
print(f"=== variant: {_args.variant} ===")


def name_from_path(path):
    """Basename -> join key. The ncAA variant's DSSP inputs are named
    `<name>_ncaa.cif`, so that suffix has to come off too or every control
    silently misses the join and lands with helix_score = NaN."""
    base = os.path.basename(str(path))
    base = base.replace("_model_0.cif", "").replace(".cif", "")
    if base.endswith("_ncaa"):
        base = base[: -len("_ncaa")]
    return base


def load_dssp(path):
    d = pd.read_csv(path, sep="\t")
    d["name"] = d["file"].apply(name_from_path)
    return d[["name", "chain_b_helix_fraction"]].rename(
        columns={"chain_b_helix_fraction": "helix_score"})


# --- 1. LIBRARY: reproduce analyze_and_score_all_metrics.py exactly ---------
boltz     = pd.read_csv(f"{METRICS_DIR}/boltz_metrics_full_library.csv")
bindcraft = pd.read_csv(f"{METRICS_DIR}/bindcraft_full_library.csv")
# v2 table carries both hit_num and hit_num_v2; without it the library rows
# would be NaN for hit_num_v2 and poison the pooled percentile ranks.
_hf = f"{METRICS_DIR}/hbond_hit_counts_v2.csv"
if not os.path.exists(_hf):
    _hf = f"{METRICS_DIR}/hbond_hit_counts.csv"
hit       = pd.read_csv(_hf, sep="\t")
count     = pd.read_csv(f"{METRICS_DIR}/enrichment_count.csv", sep="\t")
face      = pd.read_csv(f"{METRICS_DIR}/full_library_face_assignment.tsv", sep="\t")
dssp      = load_dssp(f"{METRICS_DIR}/dssp_full_library.csv")
seqs      = pd.read_csv(TSV_OUT, sep="\t", header=None, names=["name", "Sequence"])

lib = seqs
for layer in [boltz, bindcraft, hit, count, face[["name", "face_call"]], dssp]:
    lib = lib.merge(layer, on="name", how="inner")

lib["n_K"] = lib["Sequence"].map(lambda x: str(x).count("K"))
lib["n_M"] = lib["Sequence"].map(lambda x: str(x).count("M"))

lib_surv = apply_filters(lib, CONFIDENCE_FILTERS)
lib_surv = apply_filters(lib_surv, PHYSICAL_FILTERS).copy()
lib_surv = apply_filters(lib_surv, HELICITY_FILTERS)
lib_surv = lib_surv.sort_values("protein_iptm", ascending=False) \
                   .drop_duplicates(subset=["Sequence"], keep="first")
lib_surv["is_control"] = False
print(f"Library survivors: {len(lib_surv)}")

# --- 2. CONTROLS: same layers, scoped to control_run/ -----------------------
c_boltz     = pd.read_csv(f"{CTRL_DIR}/{V['boltz']}")
c_bindcraft = pd.read_csv(f"{CTRL_DIR}/{V['bindcraft']}")
c_hit       = pd.read_csv(f"{CTRL_DIR}/{V['hit']}", sep="\t")
c_face      = pd.read_csv(f"{CTRL_DIR}/{V['face']}", sep="\t")
c_dssp      = load_dssp(f"{CTRL_DIR}/{V['dssp']}")
c_seqs      = pd.read_csv(f"{CTRL_DIR}/control_list.tsv", sep="\t", header=None,
                          names=["name", "Sequence"])
truth       = pd.read_csv(f"{CTRL_DIR}/control_face_truth.tsv", sep="\t")

ctrl = c_seqs
for layer in [c_boltz, c_bindcraft, c_hit,
              c_face[["name", "face_call"]].rename(
                  columns={"face_call": "structural_face_call"}),
              c_dssp, truth]:
    ctrl = ctrl.merge(layer, on="name", how="left")   # left: never drop a control

missing = ctrl[ctrl["protein_iptm"].isna() | ctrl["interface_dG"].isna()]
if len(missing):
    print(f"WARNING: {len(missing)} controls missing metrics: {list(missing['name'])}")

ctrl["n_K"] = ctrl["Sequence"].map(lambda x: str(x).count("K"))
ctrl["n_M"] = ctrl["Sequence"].map(lambda x: str(x).count("M"))
ctrl["count"] = float("nan")     # no library enrichment for a literature control
ctrl["is_control"] = True
ctrl["face_call"] = ctrl["face"]  # pool with the literature face
print(f"Controls loaded: {len(ctrl)}")

# --- 3. FILTER AUDIT for the controls --------------------------------------
audit = filter_audit(ctrl, ALL_FILTERS)
ctrl["n_failed_filters"] = audit["n_failed"]
ctrl["failed_filters"] = audit["failed_filters"]
# Which clauses fail once the library-design-only clauses are set aside?
binding_only = {k: v for k, v in ALL_FILTERS.items() if k not in LIBRARY_DESIGN_FILTERS}
binding_audit = filter_audit(ctrl, binding_only)
ctrl["failed_binding_filters"] = binding_audit["failed_filters"]
ctrl["passes_binding_filters"] = binding_audit["n_failed"] == 0

# --- 4. POOL + SCORE per face ----------------------------------------------
combined = pd.concat([lib_surv, ctrl], ignore_index=True, sort=False)

scored = {}
for face_name in ["cmyb", "mll"]:
    pool = combined[combined["face_call"] == face_name].copy()
    if pool.empty:
        continue
    pool = add_both_priority_scores(pool)
    pool = add_swap_sensitive_score(pool)
    pool = add_score_variants(pool)
    # where does each row sit within this face's pool on the no-enrichment score?
    pool["pct_of_face_pool"] = pool["priority_score_no_enrichment"].rank(pct=True)
    scored[face_name] = pool
    n_c = int(pool["is_control"].sum())
    print(f"{face_name}: {len(pool)} in pool ({len(pool)-n_c} library + {n_c} controls)")

preview_cols = ["name", "Sequence", "is_control", "label", "kd_or_ki_uM",
                "count", "hit_num", "hit_num_v2", "helix_score", "interface_delta_unsat_hbonds",
                "interface_dG", "protein_iptm",
                "priority_score", "priority_score_no_enrichment",
                "priority_score_swap_sensitive", "pct_of_face_pool"] + \
               [f"priority_score{v}_no_enrichment"
                for v in list(SCORE_VARIANTS) + list(SCORE_VARIANTS_V2) if v] + \
               [f"priority_score{v}"
                for v in list(SCORE_VARIANTS) + list(SCORE_VARIANTS_V2) if v]

for face_name, out in [("cmyb", f"cmyb_candidates_with_controls{TAG}.csv"),
                       ("mll", f"mll_candidates_with_controls{TAG}.csv")]:
    if face_name not in scored:
        continue
    d = scored[face_name].sort_values("priority_score_no_enrichment", ascending=False)
    d[[c for c in preview_cols if c in d.columns]].to_csv(f"{CTRL_DIR}/{out}", index=False)

# --- 5. CONTROL-ONLY TABLE --------------------------------------------------
all_scored = pd.concat(scored.values(), ignore_index=True, sort=False)
ctrl_out = all_scored[all_scored["is_control"]].copy()
ctrl_cols = ["name", "Sequence", "face_call", "structural_face_call", "label",
             "kd_or_ki_uM", "measurement_type", "hit_num", "hit_num_v2", "helix_score",
             "confidence_score", "pep_ptm", "complex_pde", "protein_iptm",
             "binder_score", "interface_dG", "interface_dSASA", "interface_sc",
             "interface_nres", "interface_interface_hbonds",
             "interface_delta_unsat_hbonds", "n_K", "n_M",
             "priority_score_no_enrichment", "priority_score_swap_sensitive",
             "priority_score_no_helix_no_enrichment", "priority_score_pde_no_enrichment",
             "pct_of_face_pool",
             "passes_binding_filters", "n_failed_filters",
             "failed_filters", "failed_binding_filters"]
ctrl_out = ctrl_out[[c for c in ctrl_cols if c in ctrl_out.columns]]
ctrl_out = ctrl_out.sort_values(["face_call", "priority_score_no_enrichment"],
                                ascending=[True, False])
ctrl_out.to_csv(f"{CTRL_DIR}/control_scores{TAG}.csv", index=False)

pd.set_option("display.width", 250)
print("\n=== CONTROLS, ranked by priority_score_no_enrichment ===")
show = ["name", "face_call", "structural_face_call", "label", "kd_or_ki_uM",
        "hit_num", "helix_score", "protein_iptm", "interface_dG",
        "priority_score_no_enrichment", "pct_of_face_pool", "passes_binding_filters"]
print(ctrl_out[[c for c in show if c in ctrl_out.columns]].to_string(index=False))
print(f"\nWrote control_scores{TAG}.csv + per-face *_with_controls{TAG}.csv to {CTRL_DIR}")
