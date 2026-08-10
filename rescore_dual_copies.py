#!/usr/bin/env python3
"""
Re-score the per-copy dual metrics WITH the enrichment count.

The dual candidates are library mRNA-display hits and have real enrichment
counts (2-36), so the count term belongs in their score. It was dropped in the
first pass because the same pool also holds decoys and literature controls,
which have no count -- but that is the situation `_no_enrichment` exists for,
and it should not have been imposed on the library peptides too.

Both families are emitted:
  priority_score*                 includes `count`  -- for the 80 library hits
  priority_score*_no_enrichment   drops it          -- the only way to compare
                                                      against decoys/controls
Rows without a count (decoys, cMyb_native, WT_MLL) get NaN in the count-bearing
columns, exactly as the literature controls do in the library comparison.
"""
import os
import pandas as pd
import sys

P = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, P)
from kix_scoring import add_score_variants

D = os.path.join(P, "dual_cofold")
m = pd.read_csv(f"{D}/dual_metrics_per_copy.csv")
m = m.drop(columns=[c for c in m.columns if c.startswith("priority_score")])

cnt = pd.read_csv(f"{D}/dual_selection.csv")[["name", "count"]].rename(
    columns={"name": "peptide"})
m = m.merge(cnt, on="peptide", how="left")
n = m["count"].notna().sum()
print(f"count merged onto {n}/{len(m)} copies "
      f"({m[m['count'].isna()]['peptide'].nunique()} peptides have none: "
      f"decoys + natives)")

m = add_score_variants(m, with_count=True)
m.to_csv(f"{D}/dual_metrics_per_copy.csv", index=False)
sc = sorted(c for c in m.columns if c.startswith("priority_score"))
print(f"wrote {len(m)} rows, {len(sc)} score columns")
for c in sc:
    print(f"  {c}: {m[c].notna().sum()} non-null")
