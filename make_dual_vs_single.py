#!/usr/bin/env python3
"""
Write dual_vs_single.csv -- single-binder and dual metrics side by side.

Every dual candidate is also a library peptide, so it already has a single-binder
score. This puts the two next to each other: how the peptide scores alone, and
what each of its two interfaces looks like in the dual complex.

Columns, in order:
  identity        peptide, Sequence, count, source_face, single_rank
  single_*        the peptide as a single-face binder (from {cmyb,mll}_candidates.csv)
  cmyb_* / mll_*  the two interfaces in the dual complex, keyed by FACE
  d_dG            dual dG on the source face minus single dG (+ = weaker as a dual)
  dG_gap          |cmyb dG - mll dG|, how balanced the two interfaces are
  dual_*          the combined dual scores

⚠️ single_score and dual_* are percentile ranks over DIFFERENT pools (face pool of
111/87 vs the 172-sub-complex pool) and are not comparable to each other. The raw
metrics are absolute and can be compared directly.

Regenerate:  /scratch/jem9759/envs/general_penv/bin/python3 make_dual_vs_single.py
"""
import os

import numpy as np
import pandas as pd

P = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(P, "dual_cofold")

SINGLE = {"interface_dG": "single_interface_dG", "hit_num_v2": "single_hit_num_v2",
          "hit_num": "single_hit_num", "helix_score": "single_helix_score",
          "complex_pde": "single_complex_pde", "protein_iptm": "single_protein_iptm",
          "interface_sc": "single_interface_sc",
          "priority_score_v2_pde": "single_priority_score_v2_pde"}
PER_FACE = ["interface_dG", "hit_num", "hit_num_v2", "helix_score",
            "protein_iptm", "complex_pde", "interface_sc", "interface_nres",
            "interface_interface_hbonds"]


def main():
    s = pd.read_csv(f"{D}/dual_shortlist.csv")
    lib = pd.concat([pd.read_csv(f"{P}/full_library_all_metrics/{f}_candidates.csv")
                     .assign(source_face=f) for f in ("cmyb", "mll")])
    lib["single_rank"] = lib.groupby("source_face")["priority_score_v2_pde"] \
                            .rank(ascending=False, method="min").astype(int)

    # {cmyb,mll}_candidates.csv carries a preview subset, so not every metric is
    # present -- take what is there rather than failing on the rest
    have = {k: v for k, v in SINGLE.items() if k in lib.columns}
    cols = ["name", "source_face", "single_rank"] + list(have)
    single = lib[cols].rename(columns={"name": "peptide", **have})
    m = s.merge(single, on=["peptide", "source_face"], how="left")

    # dual dG on the face this peptide was originally ranked for
    same = np.where(m.source_face == "cmyb", m.cmyb_interface_dG, m.mll_interface_dG)
    m["d_dG"] = same - m.single_interface_dG          # + = weaker as a dual
    m["dG_gap"] = (m.cmyb_interface_dG - m.mll_interface_dG).abs()

    out = (["peptide", "Sequence", "count", "source_face", "single_rank"]
           + [v for k, v in SINGLE.items() if k in lib.columns]
           + [f"{f}_{c}" for f in ("cmyb", "mll") for c in PER_FACE]
           + ["d_dG", "dG_gap", "cmyb_copy", "mll_copy"]
           + [c for c in m.columns if c.startswith("dual_")])
    out = [c for c in out if c in m.columns]
    m = m[out].sort_values("dual_priority_score_v2_pde", ascending=False)
    m.to_csv(f"{D}/dual_vs_single.csv", index=False)
    print(f"wrote {D}/dual_vs_single.csv  ({len(m)} rows x {len(m.columns)} cols)")


if __name__ == "__main__":
    main()
