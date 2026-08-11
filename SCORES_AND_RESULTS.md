# Scores & Results Reference

**Status:** ✅ auto-generated from `kix_scoring.py` — formulas are the live definitions, not a transcription. **Regenerate:** `make_score_reference.py`. **Written:** 2026-08-10.

*Companion to `ControlPlan.md` (findings and evidence) and `DUAL_VS_SINGLE.md` (how the dual binders compare with the single-face binders). This file is the map: what each score means and which file to open.*

---

## How a score is built

Every score is a **percentile-rank blend**: each metric becomes its rank within the population being scored (0–1), then the ranks are weighted and summed.

⚠️ **Scores are relative to the pool they were computed in.** The same peptide gets a different number in the library-only pool than in the library+controls pool. Never compare a score across two files without checking they share a pool.

Column names come from two independent switches:

| suffix | meaning |
|---|---|
| *(none)* | 4th slot = `helix_score` — the original definition |
| `_no_helix` | 4th slot dropped, remaining weights renormalised |
| `_pde` | 4th slot = `complex_pde` (lower is better) |
| `_v2` | `hit_num_v2` in place of `hit_num` |
| `_no_enrichment` | drops the library `count` term — **required to compare controls to library** |

## The formulas

- **`priority_score`**
  `0.22·count + 0.22·hit_num + 0.11·(1−interface_delta_unsat_hbonds) + 0.22·protein_iptm + 0.22·helix_score`
- **`priority_score_no_enrichment`**
  `0.29·hit_num + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm + 0.29·helix_score`
- **`priority_score_no_helix`**
  `0.29·count + 0.29·hit_num + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm`
- **`priority_score_no_helix_no_enrichment`**
  `0.40·hit_num + 0.20·(1−interface_delta_unsat_hbonds) + 0.40·protein_iptm`
- **`priority_score_pde`**
  `0.22·count + 0.22·hit_num + 0.11·(1−interface_delta_unsat_hbonds) + 0.22·protein_iptm + 0.22·(1−complex_pde)`
- **`priority_score_pde_no_enrichment`**
  `0.29·hit_num + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm + 0.29·(1−complex_pde)`
- **`priority_score_v2`**
  `0.22·count + 0.22·hit_num_v2 + 0.11·(1−interface_delta_unsat_hbonds) + 0.22·protein_iptm + 0.22·helix_score`
- **`priority_score_v2_no_enrichment`**
  `0.29·hit_num_v2 + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm + 0.29·helix_score`
- **`priority_score_v2_no_helix`**
  `0.29·count + 0.29·hit_num_v2 + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm`
- **`priority_score_v2_no_helix_no_enrichment`**
  `0.40·hit_num_v2 + 0.20·(1−interface_delta_unsat_hbonds) + 0.40·protein_iptm`
- **`priority_score_v2_pde`**
  `0.22·count + 0.22·hit_num_v2 + 0.11·(1−interface_delta_unsat_hbonds) + 0.22·protein_iptm + 0.22·(1−complex_pde)`
- **`priority_score_v2_pde_no_enrichment`**
  `0.29·hit_num_v2 + 0.14·(1−interface_delta_unsat_hbonds) + 0.29·protein_iptm + 0.29·(1−complex_pde)`

`(1−x)` marks metrics where **lower is better** (`complex_pde`, `interface_delta_unsat_hbonds`).

⚠️ **`interface_dG` appears in NO composite score** — it is only ever a filter (`< -25`). On the dual-cofold controls `interface_dG` alone ordered natives / measured non-binders / nonsense correctly while the composite did not, because the composite is 58% Boltz-confidence terms. See `ControlPlan.md` §8.

## Which score should I use?

| purpose | column |
|---|---|
| current production ranking | `priority_score` |
| comparing controls to library | `priority_score_*_no_enrichment` |
| best validated on controls | `priority_score_v2_pde_no_enrichment` (AUC 0.788, p=0.011) |
| dual-face candidates | **`interface_dG` of the weaker copy** — beats every composite |

## Hit-residue sets

Selected with `hbond_hit_num.py --residue-set`:

- **`original`** — c-Myb `[14, 18, 21, 65, 69, 72, 73, 76]`, MLL `[27, 39, 43, 46, 71, 75, 79]`
- **`crystal`** — c-Myb `[9, 21, 61, 76, 80, 81]`, MLL `[27, 39, 83, 84, 86]`
- **`union`** — c-Myb `[9, 14, 18, 21, 61, 65, 69, 72, 73, 76, 80, 81]`, MLL `[27, 39, 43, 46, 71, 75, 79, 83, 84, 86]`

`original` is the default and the best performer. `crystal` (derived from residues that actually H-bond in 2AGH) makes the metric self-consistent but **destroys discrimination** — see `ControlPlan.md` §6.

## Which file do I open?

| file | what it is |
|---|---|
| `full_library_all_metrics/cmyb_candidates.csv` | **Library shortlist, c-Myb face** (111). Carries all 12 score variants. |
| `full_library_all_metrics/mll_candidates.csv` | **Library shortlist, MLL face** (87). Carries all 12 score variants. |
| `control_run/control_scores_ncaa.csv` | The 24 literature controls, ncAA run, with per-clause filter audit. |
| `control_run/mll_candidates_with_controls_ncaa.csv` | Controls pooled **with** library survivors — the comparison that matters. |
| `control_run/cmyb_candidates_with_controls_ncaa.csv` | Same for the c-Myb face (1 control only). |
| `dual_cofold/dual_shortlist.csv` | **40 dual-face candidates**, ranked on the weaker of the two interfaces. |
| `dual_cofold/dual_summary.csv` | One row per peptide: weaker/stronger copy metrics, face pair. |
| `dual_cofold/dual_metrics_per_copy.csv` | One row per *copy* (172) — the full metric set for each interface. |

⚠️ **Mind the suffixes in `control_run/`.** `_ncaa` = the ResidueX-grafted re-run (**use these**); no suffix = the original canonical stand-in run, kept for comparison; `_v2` = includes `hit_num_v2`.

Per-stage inputs, rarely opened directly:

| pattern | stage |
|---|---|
| `*boltz*.csv` | Boltz confidences (`protein_iptm`, `complex_pde`, `pep_ptm`) |
| `*bindcraft*.csv` | BindCraft interface metrics (`interface_dG`, `_sc`, `_nres`, hbonds) |
| `*hit*.csv`, `*interactions_clean*` | Schrödinger H-bond/pi-pi → `hit_num`, `hit_num_v2` |
| `*hydrophobic*` | apolar hit-residue contacts — the extra term in `hit_num_v2` |
| `*dssp*` | helicity (`chain_b_helix_fraction`, `chain_b_helix_count`) |
| `*face_assignment*`, `*face*` | which KIX face each peptide bound |

## Full inventory

**`full_library_all_metrics/`** — 10 files

- `bindcraft_full_library.csv` (31392 rows)
- `boltz_metrics_full_library.csv` (31392 rows)
- `cmyb_candidates.csv` (111 rows)
- `dssp_full_library.csv` (31392 rows)
- `enrichment_count.csv` (31392 rows)
- `full_library_face_assignment.tsv` (31392 rows)
- `hbond_hit_counts.csv` (30776 rows)
- `hbond_hit_counts_v2.csv` (30776 rows)
- `hydrophobic_contacts_full_library.tsv` (31392 rows)
- `mll_candidates.csv` (87 rows)

**`control_run/`** — 22 files

- `cmyb_candidates_with_controls.csv` (112 rows)
- `cmyb_candidates_with_controls_ncaa.csv` (112 rows)
- `control_bindcraft_data.csv` (24 rows)
- `control_bindcraft_ncaa.csv` (24 rows)
- `control_boltz_data.csv` (24 rows)
- `control_dssp_data.csv` (24 rows)
- `control_dssp_ncaa.csv` (24 rows)
- `control_face_assignment.tsv` (24 rows)
- `control_face_assignment_ncaa.tsv` (24 rows)
- `control_face_literature.tsv` (24 rows)
- `control_face_truth.tsv` (24 rows)
- `control_hydrophobic.tsv` (24 rows)
- `control_interactions_clean.csv` (24 rows)
- `control_interactions_clean_ncaa.csv` (24 rows)
- `control_interactions_clean_ncaa_v2.csv` (24 rows)
- `control_list.tsv` (23 rows)
- `control_scores.csv` (24 rows)
- `control_scores_ncaa.csv` (24 rows)
- `mll_candidates_with_controls.csv` (110 rows)
- `mll_candidates_with_controls_ncaa.csv` (110 rows)
- `ncaa_smiles.csv` (7 rows)
- `ncaa_swap_summary.csv` (24 rows)

**`dual_cofold/`** — 17 files

- `copies_bindcraft.csv` (172 rows)
- `copies_dssp.csv` (172 rows)
- `copies_face.tsv` (172 rows)
- `copies_face_call.tsv` (172 rows)
- `copies_hit.csv` (172 rows)
- `copies_hydrophobic.tsv` (172 rows)
- `decoy_face_assignment.tsv` (4 rows)
- `decoy_list.tsv` (3 rows)
- `decoy_new.tsv` (1 rows)
- `dual_cofold_results.csv` (86 rows)
- `dual_controls.tsv` (1 rows)
- `dual_face_assignment.tsv` (82 rows)
- `dual_list.tsv` (79 rows)
- `dual_metrics_per_copy.csv` (172 rows)
- `dual_selection.csv` (80 rows)
- `dual_shortlist.csv` (40 rows)
- `dual_summary.csv` (86 rows)

