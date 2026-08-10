#!/usr/bin/env python3
"""
Regenerate SCORES_AND_RESULTS.md.

The score formulas are read live from kix_scoring.py rather than transcribed, so
the doc cannot drift from the code. Re-run after changing SCORE_VARIANTS or
adding result files:

    /scratch/jem9759/envs/general_penv/bin/python3 make_score_reference.py

Formatting deliberately mirrors ControlPlan.md (**Status:** header, ## sections,
numbered bold entries, warning callouts) so moving between the two is seamless.
"""
import glob
import os
import sys
from datetime import date

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kix_scoring import SCORE_VARIANTS, SCORE_VARIANTS_V2, _HIT_RESIDUE_NUMBERS

P = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(P, "SCORES_AND_RESULTS.md")

PAYOFF = [
    ("full_library_all_metrics/cmyb_candidates.csv",
     "**Library shortlist, c-Myb face** (111). Carries all 12 score variants."),
    ("full_library_all_metrics/mll_candidates.csv",
     "**Library shortlist, MLL face** (87). Carries all 12 score variants."),
    ("control_run/control_scores_ncaa.csv",
     "The 24 literature controls, ncAA run, with per-clause filter audit."),
    ("control_run/mll_candidates_with_controls_ncaa.csv",
     "Controls pooled **with** library survivors — the comparison that matters."),
    ("control_run/cmyb_candidates_with_controls_ncaa.csv",
     "Same for the c-Myb face (1 control only)."),
    ("dual_cofold/dual_shortlist.csv",
     "**40 dual-face candidates**, ranked on the weaker of the two interfaces."),
    ("dual_cofold/dual_summary.csv",
     "One row per peptide: weaker/stronger copy metrics, face pair."),
    ("dual_cofold/dual_metrics_per_copy.csv",
     "One row per *copy* (172) — the full metric set for each interface."),
]

STAGES = [
    ("`*boltz*.csv`", "Boltz confidences (`protein_iptm`, `complex_pde`, `pep_ptm`)"),
    ("`*bindcraft*.csv`", "BindCraft interface metrics (`interface_dG`, `_sc`, `_nres`, hbonds)"),
    ("`*hit*.csv`, `*interactions_clean*`", "Schrödinger H-bond/pi-pi → `hit_num`, `hit_num_v2`"),
    ("`*hydrophobic*`", "apolar hit-residue contacts — the extra term in `hit_num_v2`"),
    ("`*dssp*`", "helicity (`chain_b_helix_fraction`, `chain_b_helix_count`)"),
    ("`*face_assignment*`, `*face*`", "which KIX face each peptide bound"),
]


def formula(spec, with_count):
    terms = ([("count", 0.2, False)] if with_count else []) + spec
    tot = sum(w for _, w, _ in terms)
    return " + ".join(
        f"{w/tot:.2f}·{'(1−' if inv else ''}{c}{')' if inv else ''}"
        for c, w, inv in terms)


def main():
    L = []
    A = L.append
    A("# Scores & Results Reference")
    A("")
    A(f"**Status:** ✅ auto-generated from `kix_scoring.py` — formulas are the live "
      f"definitions, not a transcription. **Regenerate:** `make_score_reference.py`. "
      f"**Written:** {date.today().isoformat()}.")
    A("")
    A("*Companion to `ControlPlan.md`, which holds the findings and the evidence "
      "behind every recommendation here. This file is the map: what each score "
      "means and which file to open.*")
    A("")
    A("---")
    A("")
    A("## How a score is built")
    A("")
    A("Every score is a **percentile-rank blend**: each metric becomes its rank within "
      "the population being scored (0–1), then the ranks are weighted and summed.")
    A("")
    A("⚠️ **Scores are relative to the pool they were computed in.** The same peptide "
      "gets a different number in the library-only pool than in the library+controls "
      "pool. Never compare a score across two files without checking they share a pool.")
    A("")
    A("Column names come from two independent switches:")
    A("")
    A("| suffix | meaning |")
    A("|---|---|")
    A("| *(none)* | 4th slot = `helix_score` — the original definition |")
    A("| `_no_helix` | 4th slot dropped, remaining weights renormalised |")
    A("| `_pde` | 4th slot = `complex_pde` (lower is better) |")
    A("| `_v2` | `hit_num_v2` in place of `hit_num` |")
    A("| `_no_enrichment` | drops the library `count` term — **required to compare controls to library** |")
    A("")
    A("## The formulas")
    A("")
    for suf, spec in {**SCORE_VARIANTS, **SCORE_VARIANTS_V2}.items():
        for wc in (True, False):
            name = "priority_score" + suf + ("" if wc else "_no_enrichment")
            A(f"- **`{name}`**")
            A(f"  `{formula(spec, wc)}`")
    A("")
    A("`(1−x)` marks metrics where **lower is better** (`complex_pde`, "
      "`interface_delta_unsat_hbonds`).")
    A("")
    A("⚠️ **`interface_dG` appears in NO composite score** — it is only ever a filter "
      "(`< -25`). On the dual-cofold controls `interface_dG` alone ordered natives / "
      "measured non-binders / nonsense correctly while the composite did not, because "
      "the composite is 58% Boltz-confidence terms. See `ControlPlan.md` §8.")
    A("")
    A("## Which score should I use?")
    A("")
    A("| purpose | column |")
    A("|---|---|")
    A("| current production ranking | `priority_score` |")
    A("| comparing controls to library | `priority_score_*_no_enrichment` |")
    A("| best validated on controls | `priority_score_v2_pde_no_enrichment` (AUC 0.788, p=0.011) |")
    A("| dual-face candidates | **`interface_dG` of the weaker copy** — beats every composite |")
    A("")
    A("## Hit-residue sets")
    A("")
    A("Selected with `hbond_hit_num.py --residue-set`:")
    A("")
    for k, v in _HIT_RESIDUE_NUMBERS.items():
        A(f"- **`{k}`** — c-Myb `{v['cmyb']}`, MLL `{v['mll']}`")
    A("")
    A("`original` is the default and the best performer. `crystal` (derived from residues "
      "that actually H-bond in 2AGH) makes the metric self-consistent but **destroys "
      "discrimination** — see `ControlPlan.md` §6.")
    A("")
    A("## Which file do I open?")
    A("")
    A("| file | what it is |")
    A("|---|---|")
    for f, d in PAYOFF:
        A(f"| `{f}` | {d} |")
    A("")
    A("⚠️ **Mind the suffixes in `control_run/`.** `_ncaa` = the ResidueX-grafted re-run "
      "(**use these**); no suffix = the original canonical stand-in run, kept for "
      "comparison; `_v2` = includes `hit_num_v2`.")
    A("")
    A("Per-stage inputs, rarely opened directly:")
    A("")
    A("| pattern | stage |")
    A("|---|---|")
    for pat, desc in STAGES:
        A(f"| {pat} | {desc} |")
    A("")
    A("## Full inventory")
    A("")
    for d in ("full_library_all_metrics", "control_run", "dual_cofold"):
        fs = sorted(glob.glob(f"{P}/{d}/*.csv") + glob.glob(f"{P}/{d}/*.tsv"))
        A(f"**`{d}/`** — {len(fs)} files")
        A("")
        for f in fs:
            try:
                n = f" ({len(pd.read_csv(f, sep=None, engine='python'))} rows)"
            except Exception:
                n = ""
            A(f"- `{os.path.basename(f)}`{n}")
        A("")
    open(OUT, "w").write("\n".join(L) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
