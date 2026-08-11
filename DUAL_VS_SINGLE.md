# Dual vs Single Binders

**Status:** ✅ generated from `dual_shortlist.csv` + the library candidate lists. **Regenerate:** `make_dual_vs_single.py`. **Written:** 2026-08-10.

*Companion to `SCORES_AND_RESULTS.md` (what each column means) and `ControlPlan.md` (the validation evidence).*

---

## How to read this

Every dual candidate is also a library peptide, so it has a single-binder score already. That lets us ask the question that matters: **what happens to a peptide's interface when a second copy of itself occupies the other face?**

⚠️ **Percentile scores are NOT comparable between the two.** Single scores rank within a face pool (111 c-Myb / 87 MLL); dual scores rank within the 172-sub-complex pool. Only raw metrics (`interface_dG`, `hit_num_v2`, `helix_score`) are absolute, so this comparison is built on those.

## 1. Dual binding costs the primary interface ~3 kcal/mol

Comparing each peptide's `interface_dG` **alone** against the same face **in the dual complex** (n=40):

| | mean `interface_dG` |
|---|---|
| as a single binder | -35.53 |
| same face, with a 2nd copy present | -32.46 |
| **change** | **+3.07** (median +1.79) |

Paired Wilcoxon **p = 0.0015** — the degradation is systematic, not noise. 18 of 40 degrade by more than the ±2.8 kcal/mol BindCraft noise floor; 5 improve.

Physically unsurprising — the two grooves are on opposite faces but the peptide still has to share one protein — but it means **a dual binder is not simply two independent single binders**, and the best single binder is not automatically the best dual.

## 2. The best single binders are mostly NOT the best duals

| face | top-10 single binders that made the dual shortlist |
|---|---|
| cmyb | 3/10 |
| mll | 4/10 |

Dual capability is a distinct property. Ranking on single-face strength would have missed most of the dual shortlist.

## 3. Best dual binders

Ranked on `dual_priority_score_v2_pde`. `Δ` is the change to the primary interface when the second copy is added (positive = weaker).

| # | peptide | sequence | count | c-Myb dG | MLL dG | single dG | Δ | single rank |
|---|---|---|---|---|---|---|---|---|
| 1 | `H2736` | `FWGNWLEWFR` | 2 | -38.7 | -35.9 | -42.4 | +3.8 | cmyb #27 |
| 2 | `H383` | `FIGHHWDWFR` | 6 | -36.9 | -35.1 | -37.5 | +0.6 | cmyb #16 |
| 3 | `H256` | `FWGNFREWLR` | 4 | -27.1 | -33.2 | -36.4 | +9.4 | cmyb #13 |
| 4 | `H1805` | `FTGYVWDWFR` | 3 | -35.4 | -35.5 | -40.6 | +5.2 | cmyb #36 |
| 5 | `H27023` | `FYGHIWDWFR` | 3 | -29.5 | -30.0 | -33.6 | +4.0 | cmyb #21 |
| 6 | `H19844` | `FWGHISEWFR` | 2 | -24.6 | -37.0 | -36.1 | +11.5 | cmyb #35 |
| 7 | `H4875` | `FWGIFWDWFR` | 3 | -30.8 | -35.3 | -37.0 | +1.7 | mll #18 |
| 8 | `H3178` | `FRGNWVDRWK` | 2 | -29.8 | -25.1 | -36.4 | +6.5 | cmyb #33 |
| 9 | `H40119` | `FWGTLWEFYK` | 2 | -47.6 | -31.3 | -43.3 | +12.0 | mll #19 |
| 10 | `H15488` | `FTGHHWEFFR` | 2 | -30.1 | -33.4 | -37.9 | +7.8 | cmyb #18 |
| 11 | `H2478` | `FWGSFWDIHK` | 3 | -29.0 | -30.7 | -26.4 | -4.3 | mll #34 |
| 12 | `H7965` | `FYGRIWDWFR` | 2 | -39.2 | -34.5 | -35.3 | -3.9 | cmyb #24 |
| 13 | `H3019` | `FWGYWFDVFR` | 2 | -37.5 | -39.6 | -41.4 | +1.8 | mll #21 |
| 14 | `H38088` | `FWGSLHDFWK` | 2 | -27.9 | -31.7 | -28.4 | +0.4 | cmyb #39 |
| 15 | `H1969` | `FNGIFWDWFR` | 3 | -27.1 | -42.0 | -33.8 | -8.2 | mll #23 |

## 4. Best single binders, and how they do as duals

**cmyb face** — top 10 by `priority_score_v2_pde`

| # | peptide | sequence | single dG | in dual shortlist? | dual dG (same face) |
|---|---|---|---|---|---|
| 1 | `H12` | `FWGIRWDWFR` | -35.1 | no | — |
| 2 | `H200` | `FWGHHWDVYK` | -30.1 | **yes** | -32.6 (-2.4) |
| 3 | `H5186` | `FTGFHHEWFR` | -32.6 | no | — |
| 4 | `H22416` | `FWGIRWDFFK` | -39.2 | **yes** | -34.7 (+4.5) |
| 5 | `H4355` | `FWGYSHEWFR` | -31.1 | no | — |
| 6 | `H6` | `FWGNWFDRWK` | -31.6 | no | — |
| 7 | `H103` | `FNGIRWDWFR` | -35.1 | no | — |
| 8 | `H26` | `FWGIRWDWFK` | -34.0 | no | — |
| 9 | `H341` | `FTGFWHDWFR` | -36.9 | **yes** | -24.7 (+12.2) |
| 10 | `H1376` | `FHGHHWERFR` | -28.6 | no | — |

**mll face** — top 10 by `priority_score_v2_pde`

| # | peptide | sequence | single dG | in dual shortlist? | dual dG (same face) |
|---|---|---|---|---|---|
| 1 | `H10` | `FWGNFRDWSK` | -37.3 | no | — |
| 2 | `H1530` | `FNGIWPDFLR` | -34.3 | **yes** | -33.0 (+1.3) |
| 3 | `H3464` | `FWGNWHDRWK` | -38.2 | no | — |
| 4 | `H273` | `FWGIWFDAFR` | -39.7 | **yes** | -28.0 (+11.7) |
| 5 | `H1198` | `FVGFLWEWFR` | -30.7 | **yes** | -34.0 (-3.3) |
| 6 | `H23585` | `FYGHSWEWFR` | -27.1 | no | — |
| 7 | `H13284` | `FNGNWPEVFK` | -36.4 | no | — |
| 8 | `H3113` | `FVGHFPEYFR` | -43.3 | no | — |
| 9 | `H143` | `FTGIFWDWHR` | -45.6 | **yes** | -28.2 (+17.4) |
| 10 | `H4852` | `FYGIWPERYK` | -31.8 | no | — |

## 5. Most balanced duals

A convincing dual binder should engage both faces comparably. |c-Myb dG − MLL dG| within the ±2.8 noise floor:

**10 of 40** qualify (median gap across all: 4.8).

| peptide | sequence | c-Myb dG | MLL dG | gap |
|---|---|---|---|---|
| `H1805` | `FTGYVWDWFR` | -35.4 | -35.5 | 0.09 |
| `H27023` | `FYGHIWDWFR` | -29.5 | -30.0 | 0.43 |
| `H273` | `FWGIWFDAFR` | -26.9 | -28.0 | 1.17 |
| `H2712` | `FWGWFRDVFR` | -29.7 | -28.3 | 1.45 |
| `H2478` | `FWGSFWDIHK` | -29.0 | -30.7 | 1.68 |
| `H3702` | `FWGNFHEIRR` | -29.4 | -27.6 | 1.72 |
| `H383` | `FIGHHWDWFR` | -36.9 | -35.1 | 1.74 |
| `H3019` | `FWGYWFDVFR` | -37.5 | -39.6 | 2.05 |
| `H143` | `FTGIFWDWHR` | -30.5 | -28.2 | 2.30 |
| `H2736` | `FWGNWLEWFR` | -38.7 | -35.9 | 2.75 |

## Caveats

- `interface_dG` carries ±2.8 kcal/mol run-to-run noise, so single differences below that are not meaningful; only the paired trend in §1 is.
- Dual metrics measure each copy **as placed in the 3-chain complex**, not whether that copy would bind on its own.
- The dual shortlist is drawn from the top 40 per face, so it cannot contain a peptide that ranked poorly as a single binder — the comparison in §2 is within that selected set.
