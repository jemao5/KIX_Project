#!/usr/bin/env python3
"""
Generate DUAL_VS_SINGLE.md -- how the best dual binders compare with the best
single-face binders.

The comparison is possible because every dual candidate is itself a library
peptide with a single-binder score already. For each one we can ask: how good is
its interface alone, and what happens to that same interface when a second copy
of itself is present on the other face?

⚠️ Percentile SCORES are not comparable across the two files -- the single scores
rank within a face pool (111 or 87), the dual scores within the 172-sub-complex
pool. Only RAW metrics (interface_dG, hit_num_v2, helix_score) are on an absolute
scale, so the comparison is built on those.

Regenerate:  /scratch/jem9759/envs/boltz_env/bin/python3 make_dual_vs_single.py
"""
import os
from datetime import date

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

P = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(P, "dual_cofold")
NOISE = 2.8   # BindCraft interface_dG run-to-run SD on identical input


def load():
    s = pd.read_csv(f"{D}/dual_shortlist.csv")
    lib = pd.concat([pd.read_csv(f"{P}/full_library_all_metrics/{f}_candidates.csv")
                     .assign(face=f) for f in ("cmyb", "mll")])
    single = lib[["name", "face", "interface_dG", "hit_num_v2", "helix_score",
                  "protein_iptm", "priority_score_v2_pde"]].rename(columns={
        "name": "peptide", "interface_dG": "single_dG", "hit_num_v2": "single_hit",
        "helix_score": "single_helix", "protein_iptm": "single_iptm",
        "priority_score_v2_pde": "single_score"})
    m = s.merge(single, left_on=["peptide", "source_face"],
                right_on=["peptide", "face"], how="left")
    m["dual_same_face_dG"] = np.where(m.source_face == "cmyb",
                                      m.cmyb_interface_dG, m.mll_interface_dG)
    m["dual_other_face_dG"] = np.where(m.source_face == "cmyb",
                                       m.mll_interface_dG, m.cmyb_interface_dG)
    m["dG_change"] = m.dual_same_face_dG - m.single_dG      # + = weaker in dual
    return m.sort_values("dual_priority_score_v2_pde", ascending=False), lib


def main():
    m, lib = load()
    L = []
    A = L.append
    A("# Dual vs Single Binders")
    A("")
    A(f"**Status:** ✅ generated from `dual_shortlist.csv` + the library candidate "
      f"lists. **Regenerate:** `make_dual_vs_single.py`. **Written:** {date.today().isoformat()}.")
    A("")
    A("*Companion to `SCORES_AND_RESULTS.md` (what each column means) and "
      "`ControlPlan.md` (the validation evidence).*")
    A("")
    A("---")
    A("")
    A("## How to read this")
    A("")
    A("Every dual candidate is also a library peptide, so it has a single-binder "
      "score already. That lets us ask the question that matters: **what happens to "
      "a peptide's interface when a second copy of itself occupies the other face?**")
    A("")
    A("⚠️ **Percentile scores are NOT comparable between the two.** Single scores rank "
      "within a face pool (111 c-Myb / 87 MLL); dual scores rank within the "
      "172-sub-complex pool. Only raw metrics (`interface_dG`, `hit_num_v2`, "
      "`helix_score`) are absolute, so this comparison is built on those.")
    A("")

    A("## 1. Dual binding costs the primary interface ~3 kcal/mol")
    A("")
    w, p = wilcoxon(m.single_dG, m.dual_same_face_dG)
    A(f"Comparing each peptide's `interface_dG` **alone** against the same face **in the "
      f"dual complex** (n={len(m)}):")
    A("")
    A("| | mean `interface_dG` |")
    A("|---|---|")
    A(f"| as a single binder | {m.single_dG.mean():.2f} |")
    A(f"| same face, with a 2nd copy present | {m.dual_same_face_dG.mean():.2f} |")
    A(f"| **change** | **{m.dG_change.mean():+.2f}** (median {m.dG_change.median():+.2f}) |")
    A("")
    A(f"Paired Wilcoxon **p = {p:.4f}** — the degradation is systematic, not noise. "
      f"{(m.dG_change > NOISE).sum()} of {len(m)} degrade by more than the "
      f"±{NOISE} kcal/mol BindCraft noise floor; {(m.dG_change < -NOISE).sum()} improve.")
    A("")
    A("Physically unsurprising — the two grooves are on opposite faces but the peptide "
      "still has to share one protein — but it means **a dual binder is not simply two "
      "independent single binders**, and the best single binder is not automatically "
      "the best dual.")
    A("")

    A("## 2. The best single binders are mostly NOT the best duals")
    A("")
    A("| face | top-10 single binders that made the dual shortlist |")
    A("|---|---|")
    for f in ("cmyb", "mll"):
        top = lib[lib.face == f].nlargest(10, "priority_score_v2_pde")
        n = sum(1 for x in top.name if x in set(m.peptide))
        A(f"| {f} | {n}/10 |")
    A("")
    A("Dual capability is a distinct property. Ranking on single-face strength would "
      "have missed most of the dual shortlist.")
    A("")

    A("## 3. Best dual binders")
    A("")
    A("Ranked on `dual_priority_score_v2_pde`. `Δ` is the change to the primary "
      "interface when the second copy is added (positive = weaker).")
    A("")
    A("| # | peptide | sequence | count | c-Myb dG | MLL dG | single dG | Δ | single rank |")
    A("|---|---|---|---|---|---|---|---|---|")
    rank = {f: {n: i + 1 for i, n in enumerate(
        lib[lib.face == f].sort_values("priority_score_v2_pde", ascending=False).name)}
        for f in ("cmyb", "mll")}
    for i, (_, r) in enumerate(m.head(15).iterrows(), 1):
        sr = rank[r.source_face].get(r.peptide, "-")
        A(f"| {i} | `{r.peptide.replace('Full_Library_Hit_','H')}` | `{r.Sequence}` | "
          f"{r['count']:.0f} | {r.cmyb_interface_dG:.1f} | {r.mll_interface_dG:.1f} | "
          f"{r.single_dG:.1f} | {r.dG_change:+.1f} | {r.source_face} #{sr} |")
    A("")

    A("## 4. Best single binders, and how they do as duals")
    A("")
    for f in ("cmyb", "mll"):
        A(f"**{f} face** — top 10 by `priority_score_v2_pde`")
        A("")
        A("| # | peptide | sequence | single dG | in dual shortlist? | dual dG (same face) |")
        A("|---|---|---|---|---|---|")
        for i, (_, r) in enumerate(
                lib[lib.face == f].nlargest(10, "priority_score_v2_pde").iterrows(), 1):
            hit = m[m.peptide == r["name"]]
            if len(hit):
                d = hit.iloc[0]
                A(f"| {i} | `{r['name'].replace('Full_Library_Hit_','H')}` | "
                  f"`{r.Sequence}` | {r.interface_dG:.1f} | **yes** | "
                  f"{d.dual_same_face_dG:.1f} ({d.dG_change:+.1f}) |")
            else:
                A(f"| {i} | `{r['name'].replace('Full_Library_Hit_','H')}` | "
                  f"`{r.Sequence}` | {r.interface_dG:.1f} | no | — |")
        A("")

    A("## 5. Most balanced duals")
    A("")
    A(f"A convincing dual binder should engage both faces comparably. "
      f"|c-Myb dG − MLL dG| within the ±{NOISE} noise floor:")
    A("")
    m2 = m.copy()
    m2["gap"] = (m2.cmyb_interface_dG - m2.mll_interface_dG).abs()
    bal = m2[m2.gap <= NOISE].sort_values("gap")
    A(f"**{len(bal)} of {len(m2)}** qualify (median gap across all: {m2.gap.median():.1f}).")
    A("")
    A("| peptide | sequence | c-Myb dG | MLL dG | gap |")
    A("|---|---|---|---|---|")
    for _, r in bal.head(10).iterrows():
        A(f"| `{r.peptide.replace('Full_Library_Hit_','H')}` | `{r.Sequence}` | "
          f"{r.cmyb_interface_dG:.1f} | {r.mll_interface_dG:.1f} | {r.gap:.2f} |")
    A("")
    A("## Caveats")
    A("")
    A(f"- `interface_dG` carries ±{NOISE} kcal/mol run-to-run noise, so single "
      f"differences below that are not meaningful; only the paired trend in §1 is.")
    A("- Dual metrics measure each copy **as placed in the 3-chain complex**, not "
      "whether that copy would bind on its own.")
    A("- The dual shortlist is drawn from the top 40 per face, so it cannot contain "
      "a peptide that ranked poorly as a single binder — the comparison in §2 is "
      "within that selected set.")
    open(f"{P}/DUAL_VS_SINGLE.md", "w").write("\n".join(L) + "\n")
    print(f"wrote {P}/DUAL_VS_SINGLE.md")


if __name__ == "__main__":
    main()
