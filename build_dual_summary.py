#!/usr/bin/env python3
"""
Collapse the per-copy dual-cofold metrics into one row per peptide.

Carries the SAME metric columns used everywhere else in the project (the set in
{cmyb,mll}_candidates.csv), rather than ad-hoc invented names.

Columns are keyed by FACE
-------------------------
Chains B and C are two copies of one sequence with identical templates -- which
face each lands on is emergent, and comes out exactly 20/20 across the shortlist.
So `copyB_*` means c-Myb on one row and MLL on the next and cannot be read down a
column. Metrics are therefore emitted three ways:

    cmyb_*    the copy on the c-Myb face   (NaN if unoccupied)
    mll_*     the copy on the MLL face     (NaN if unoccupied)
    weaker_*  the bottleneck of the two    -- what the dual score ranks on
    copyB_/copyC_  raw, by chain id, kept for traceability

`cmyb_copy` / `mll_copy` record which chain ended up where.

Dual score
----------
A peptide only binds both faces if BOTH copies bind well, so the dual score is
the **minimum** of the two copies' scores -- the weaker interface is the
bottleneck. Taking a mean would let one excellent site hide a dead one.

    dual_priority_score_v2_pde = min(copyB, copyC) of
                                 priority_score_v2_pde_no_enrichment

That is the same blend used to pick the top 40 per face, unchanged. `_dG` and
`_v2` variants of the same idea are emitted alongside for comparison.

⚠️ The component scores are percentile ranks within the 172-sub-complex pool, so
they are NOT comparable to the library-pool scores in {cmyb,mll}_candidates.csv.
"""
import os
import numpy as np
import pandas as pd

P = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(P, "dual_cofold")

# the standard metric set, same order as {cmyb,mll}_candidates.csv
METRICS = ["hit_num", "hit_num_v2", "helix_score", "complex_pde",
           "interface_delta_unsat_hbonds", "interface_dG", "protein_iptm",
           "interface_sc", "interface_nres", "interface_interface_hbonds",
           "binder_score", "cmyb_contacts", "mll_contacts", "face_call"]
SCORES = ["priority_score_no_enrichment", "priority_score_pde_no_enrichment",
          "priority_score_v2_no_enrichment", "priority_score_v2_pde_no_enrichment"]
# lower-is-better metrics: the "weaker" copy is the one with the HIGHER value
LOWER_BETTER = {"interface_dG", "complex_pde", "interface_delta_unsat_hbonds",
                "binder_score"}


def main():
    m = pd.read_csv(f"{D}/dual_metrics_per_copy.csv")
    keep = [c for c in METRICS + SCORES if c in m.columns]

    # Columns are keyed by FACE, not by chain id. Chains B and C carry identical
    # sequence and template -- which face each lands on is emergent, and comes out
    # 20/20 across the shortlist. So `copyB_interface_dG` means c-Myb on one row
    # and MLL on the next, which makes the column meaningless to read down.
    # `cmyb_*` / `mll_*` are always the same face; `copy` records which chain it was.
    wide = None
    for cp in ("copyB", "copyC"):
        s = m[m["copy"] == cp][["peptide", "group"] + keep].copy()
        s = s.rename(columns={c: f"{cp}_{c}" for c in keep})
        wide = s if wide is None else wide.merge(s.drop(columns=["group"]), on="peptide")

    b_is_cmyb = wide["copyB_face_call"] == "cmyb"
    for c in keep:
        # the copy sitting on each face (NaN where no copy landed there)
        wide[f"cmyb_{c}"] = wide[f"copyB_{c}"].where(b_is_cmyb, wide[f"copyC_{c}"])
        wide[f"mll_{c}"] = wide[f"copyC_{c}"].where(b_is_cmyb, wide[f"copyB_{c}"])
    wide["cmyb_copy"] = np.where(b_is_cmyb, "B", "C")
    wide["mll_copy"] = np.where(b_is_cmyb, "C", "B")
    # blank the face columns where that face was never occupied
    for face in ("cmyb", "mll"):
        miss = wide[f"{face}_face_call"] != face
        for c in keep:
            if c != "face_call":
                wide.loc[miss, f"{face}_{c}"] = np.nan

    # bottleneck across the two copies, per metric
    for c in keep:
        b, cc = wide[f"copyB_{c}"], wide[f"copyC_{c}"]
        if c == "face_call":
            # SORTED: a dual binder is cmyb+mll regardless of which copy is which,
            # so the unsorted pair would silently drop half the qualifying peptides.
            wide["faces"] = ["+".join(sorted((x, y))) for x, y in zip(b, cc)]
            continue
        if c in LOWER_BETTER:
            wide[f"weaker_{c}"] = wide[[f"copyB_{c}", f"copyC_{c}"]].max(axis=1)
        else:
            wide[f"weaker_{c}"] = wide[[f"copyB_{c}", f"copyC_{c}"]].min(axis=1)

    # the dual scores: a dual binder is only as good as its weaker interface
    for s in SCORES:
        wide[f"dual_{s.replace('_no_enrichment','')}"] = \
            wide[[f"copyB_{s}", f"copyC_{s}"]].min(axis=1)

    sel = pd.read_csv(f"{D}/dual_selection.csv")[
        ["name", "Sequence", "source_face", "face_rank"]].rename(columns={"name": "peptide"})
    wide = wide.merge(sel, on="peptide", how="left")

    lead = ["peptide", "Sequence", "group", "source_face", "face_rank", "faces",
            "cmyb_copy", "mll_copy"]
    dual = [c for c in wide.columns if c.startswith("dual_")]
    byface = [c for c in wide.columns
              if (c.startswith("cmyb_") or c.startswith("mll_")) and c not in lead]
    weak = [c for c in wide.columns if c.startswith("weaker_")]
    rest = [c for c in wide.columns if c not in lead + dual + byface + weak]
    wide = wide[lead + dual + byface + weak + rest].sort_values(
        "dual_priority_score_v2_pde", ascending=False)
    wide.to_csv(f"{D}/dual_summary.csv", index=False)
    print(f"dual_summary.csv: {len(wide)} peptides x {len(wide.columns)} cols")

    # shortlist: copies on DIFFERENT faces, both interfaces past the library's own
    # dG threshold, and both beating the best decoy
    dec = wide[wide.group == "decoy"]["weaker_interface_dG"].min()
    ok = wide[(wide.group == "candidate") & (wide.faces == "cmyb+mll")
              & (wide.weaker_interface_dG < -25.0) & (wide.weaker_interface_dG < dec)]
    ok.to_csv(f"{D}/dual_shortlist.csv", index=False)
    print(f"dual_shortlist.csv: {len(ok)} candidates (decoy floor {dec:.2f})")


if __name__ == "__main__":
    main()
