# KIX_Project — CLAUDE.md

Computational pipeline to design and screen **helical peptide binders against the KIX
domain of CBP**. A large library (~31k sequences) is cofolded with KIX, scored by several
structural/energetic tools, filtered, deduped, and ranked into two per-face shortlists.

Target structure: `2agh` (KIX domain). Local templates: `2agh_model1.pdb` / `2agh_model1.cif`.
In 2AGH, KIX is chain **B**, which Boltz renames to subchain **B1** on PDB input — this
matters for template configuration (see gotchas).

## The two "faces" — don't conflate them

Two peptide binding faces are screened. There are **two separate mechanisms**, and they are
easy to mix up:

1. **Design/constraint prior (by hit number)** — in `boltz_yaml_gen.py`. Parses the integer
   from the peptide name and, only when `--use_constraints` is set, chooses which pocket
   residues to constrain: `CMYB_MAX_HIT = 35` → hits ≤ 35 use c-Myb residues, hits > 35 use
   MLL residues. This is an *input-time* choice based on how the library was designed.
2. **Structural face call (post-hoc, geometric)** — in `full_library_face_determination.py`.
   Counts actual heavy-atom contacts (≤ 5.0 Å, ≥ 2 contacting residues) between the *cofolded*
   peptide and each face's residues, emitting `face_call ∈ {cmyb, mll, both, neither}`. **This
   `face_call` is what the final scoring splits on**, not the hit number.

Face residues (used by both mechanisms):
- **c-Myb face**: `[14, 18, 21, 65, 69, 72, 73, 76]`
- **MLL face**:   `[27, 39, 43, 46, 71, 75, 79]`

## Environment & cluster

SLURM cluster. `.sh` files are `sbatch` job scripts; `.py` files do the work. Standard preamble:
```bash
module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/<env>
```
Conda envs (canonical location `/scratch/jem9759/envs/`):
- `boltz_env` — Boltz-2 structure prediction (needs GPU: `#SBATCH --gres=gpu:1`)
- `BindCraft` — BindCraft interface scoring
- `dssp_env` — DSSP helicity, **and** the BioPython CIF-parsing scripts (face determination, cif2pdb)
- `general_penv` — general pandas/PyYAML/analysis work
- Schrödinger scripts use Schrödinger's own Python (`$SCHRODINGER/run`), not a conda env.

The old home env `/home/jem9759/ZhangLabWork/KIX_Project/penv` was retired (too many inodes in
`$HOME`) and replaced by `/scratch/jem9759/envs/general_penv`. Runners were repointed
accordingly. `run_schrodinger_helix_gen.sh` still `cd`s into the old home project dir so its
relative inputs resolve — that dir still exists with the files it needs, so leave it unless you
delete the home copy.

SLURM account: `torch_pr_149_chemistry`.

**Storage:** `$HOME` is only ~50GB — keep large outputs (Boltz predictions, etc.) under
`/scratch/jem9759/`. Redirect temp off the tiny `/tmp`: `export TMPDIR=/scratch/jem9759/tmp`.

## Pipeline (in order) — script · runner · env

| # | Stage | `.py` | `run_*.sh` | env |
|---|-------|-------|-----------|-----|
| 1 | xlsx → TSV (`name<TAB>seq`) | `extract_sequences.py`, `full_library_to_tsv.py`; `control_list_to_tsv.py` (controls) | `run.sh` | general |
| 2 | Build idealized helix per peptide | `make_helix_schrodinger.py` | `run_schrodinger_helix_gen.sh` | Schrödinger |
| 2b | PDB↔CIF interconvert (Boltz needs CIF w/ sequence) | `pdb_to_cif.py`, `cif2pdb.py` | `run_pdb_to_cif.sh`, `run_cif2pdb.sh` | dssp_env |
| 3 | One Boltz YAML per peptide | `boltz_yaml_gen.py` | `run_boltz_yaml_gen.sh` | penv/general |
| 3b | Split YAMLs for array jobs | — | `chunk_yamls.sh` | — |
| 4 | Run Boltz-2 (GPU) | — | `run_boltz.sh`, `run_chunked_boltz.sh` | boltz_env |
| 5 | Parse Boltz confidences | `parse_boltz_results.py` (+`_pipeline_vers`) | `run_boltz_analysis.sh` | penv/general |
| 6 | BindCraft interface scoring | `parse_bindcraft_results.py` (+ `score_peptide.py`) | `run_bindcraaft_score.sh` | BindCraft |
| 7 | Schrödinger H-bond / pi-pi + hit counts | `schrodinger_calc_hbond.py` → `hbond_hit_num.py`; merge via `merge_schrodinger_calc_hbond_pckl.py` | `run_schrodinger_prepwizard.sh`, `run_schrodinger_calc_hbond.sh` | Schrödinger |
| 8 | DSSP helicity (`chain_b_helix_fraction`) | `check_helix_dssp.py`, `check_helix_dssp_skip_first2.py` | `run_dssp.sh` | dssp_env |
| 9 | Structural face assignment (`face_call`) | `full_library_face_determination.py` | `run_full_library_face_determination.sh` | dssp_env |
| 10 | Merge + filter + rank | `analyze_and_score_all_metrics.py` | — (run directly) | general |
| 11 | Align top candidates to `2agh` for viewing | `align_structures.py` | — (run directly) | Schrödinger |

Counting helpers: `get_enrichment_count.py`, `get_aggregate_counts.py` produce the `count`
(enrichment) column consumed in scoring.

## Boltz YAML config (verified, in `boltz_yaml_gen.py`)

- `KIX_SEQUENCE` constant (line ~26) must match the exact construct used.
- KIX = chain **A**, templated from `2agh_model1.pdb` with `template_id: B1` (the PDB subchain rename).
- Peptide = chain **B**, `msa: empty` (de novo), templated from its own helix CIF in `pdb_sts_cif/`
  with `chain_id: B` + `template_id: Bxp`. Boltz auto-matches the file chain by sequence; the CIF
  **must** carry a populated `_entity_poly_seq` (raw Schrödinger PDBs lack it → IndexError in Boltz).
- KIX MSA: either `--use_msa_server` at runtime, or precompute and pass `--kix-msa kix_msa.csv`
  (avoids redundant per-peptide MSA server calls).
- `--use_constraints` adds a `pocket` constraint (`max_distance: 6`, `force: True`) on the
  face residues chosen by the hit-number rule above.

## Final scoring (`analyze_and_score_all_metrics.py`) — the payoff

Merges all metric layers **inner join on `name`**, then:

**Confidence filter:** `confidence_score > 0.85`, `pep_ptm > 0.80`, `complex_pde < 0.5`.

**Physical filter:** `binder_score < 0`, `interface_dG < -25`, `interface_dSASA > 1`,
`surface_hydrophobicity < 1`, `interface_sc > 0.5`, `interface_nres > 4`,
`interface_interface_hbonds > 1`, `interface_delta_unsat_hbonds ≤ 5`, `n_K ≤ 3`, `n_M ≤ 3`,
`hit_num ≥ 1`, `count > 1`, `protein_iptm > 0.85`.

**Helicity (DSSP) filter:** applied after the physical filter — `helix_score > 0.70`
(the DSSP `chain_b_helix_fraction` from stage 8). Tune this threshold in the script directly.

⚠️ **The thresholds and the scoring functions now live in `kix_scoring.py`**, not inline in
`analyze_and_score_all_metrics.py`. That module is imported by both the library scoring and
`analyze_controls.py`, so a control is audited against literally the same predicates the
library was filtered with — edit thresholds there, in one place.

Then **dedup by `Sequence`** (keep highest `protein_iptm`), **split by `face_call`**
(`cmyb` / `mll`), and rank each face independently by two percentile-rank blends:

- `priority_score` = `0.2·count + 0.2·hit_num + 0.1·(1−unsat_hbonds) + 0.2·protein_iptm + 0.2·helix_score`
- `priority_score_no_enrichment` = the same blend with the `count` term **dropped** and the
  remaining four weights renormalized (÷ 0.7). This exists so peptides with no library
  enrichment count — the literature controls — can be ranked on the same scale as library
  members. In `add_priority_score`, a weight of `0.0` removes its term entirely rather than
  multiplying by zero, so a NaN `count` cannot poison the no-enrichment score.

Outputs `cmyb_candidates.csv` and `mll_candidates.csv` in `full_library_all_metrics/`.
As of the `helix_score > 0.70` filter this is **198 candidates total (111 c-Myb / 87 MLL)**,
not the ~1.3k from before the helicity filter was added. Preview columns: `name, Sequence,
count, hit_num, helix_score, interface_delta_unsat_hbonds, interface_dG, protein_iptm,
priority_score, priority_score_no_enrichment`.

## Structure alignment for viewing (`align_structures.py`)

Post-scoring visualization prep, run directly with Schrödinger's Python (no `.sh` runner).
Takes the top `n_candidates` (=20) rows from each of `mll_candidates.csv` / `cmyb_candidates.csv`,
indexes every Boltz `*_model_0.cif` under `boltz_out_full/`, and structurally aligns the matched
cofolded structures onto the `2agh_model1.cif` reference via `structalign2.align_many`.

Outputs to `aligned_structures/`:
- `alignment_data.csv` — one row per candidate: `name, rmsd, score`.
- `<name>.pdb` — each aligned structure, **written as PDB (not CIF)**.

**HID/HIE/HIP → HIS fix:** `standardize_his_names()` rewrites Schrödinger's protonation-state
histidine names (`HID`/`HIE`/`HIP`) back to standard `HIS ` (space-padded `pdbres`) before writing.
Without this, ChimeraX doesn't recognize the residues and breaks the ribbon. Writing PDB (rather
than CIF) plus this rename is what makes the files open cleanly in ChimeraX.

This is **Schrödinger-specific**: the source Boltz CIF has plain `HIS`, but Schrödinger's
`StructureReader`/`StructureWriter` assign protonation-state-specific histidine names when a
structure is round-tripped through its toolkit. So the rename is only needed because these
structures pass through Schrödinger's Python here — any stage that writes structures via
Schrödinger (not just this one) can reintroduce `HID`/`HIE`/`HIP`.

## Key directories

- `full_library_all_metrics/` — per-metric CSVs + final candidate lists (the payoff).
  Note the big layer CSVs (`boltz_metrics_full_library.csv`, `bindcraft_full_library.csv`,
  `dssp_full_library.csv`) and derived (`enrichment_count.csv`, `hbond_hit_counts.csv`,
  `full_library_face_assignment.tsv`).
- `final_metric_outputs/` — where several stage-9/8 runners write first
  (`run_dssp.sh`, `run_full_library_face_determination.sh`). ⚠️ `analyze_and_score` reads from
  `full_library_all_metrics/`, so outputs must be copied/moved there before final scoring.
- `boltz_out_full/` (chunked, `chunk_0000…`), `boltz_outputs/`, `bindcraft_out/`, `bindcraft_relaxed/` — raw tool outputs.
- `aligned_structures/` — stage-11 outputs: top-20-per-face candidates aligned to `2agh`
  (`<name>.pdb`, HIS-normalized for ChimeraX) plus `alignment_data.csv` (rmsd/score).
- `yaml_chunks/`, `schrodinger_calc_hbond_chunks/` — split inputs/outputs for SLURM array jobs.
- `tsv_outputs/` — `full_library_out.tsv` is the master library (~31k rows).
- `boltz_reference_structures/` — the two **native-control** binders run through the pipeline
  alongside the library (see below).
- `control_run/` — the **literature control set** (`control_list.csv`, 24 peptides) run
  through the full pipeline (see below).
- `logs/`, `*.out` — SLURM job logs (named `*_<jobid>.out`).

## Native controls (`boltz_reference_structures/`)

Two native binders are screened as controls: `kix_peptide_reference_cmyb`
(`KEKRIKELELLLMSTENELKGQQAL`) and `kix_peptide_reference_mll`
(`SDDGNILPSDIMDFVLKNTPSMQALGESPES`). `name` = the Boltz parent folder (no `Full_Library_Hit_N`
convention here). Everything for this run is self-contained under `boltz_reference_structures/`
so the full-library files are never touched:

- Boltz/BindCraft/DSSP for the controls: `boltz_results_*/`, `bindcraft_reference_out/`,
  `dssp_reference_data.csv`, `boltz_reference_data.csv`.
- Stage-7 for the controls was run **in-session, no sbatch** (2 structures): `cif2pdb.py` →
  `run-schrodinger-2025.4.bash run prepwizard` (same flags as `run_schrodinger_prepwizard.sh`) →
  `run-schrodinger-2025.4.bash run python3 schrodinger_calc_hbond.py reference_hits.dat
  reference_interactions.pkl --start 1 --end 2`. No `merge_*` step (single chunk).
- `hbond_hit_num_reference.py` (project root) is a paths-scoped copy of `hbond_hit_num.py`:
  reads `reference_interactions.pkl` + `reference_face_assignment.tsv`, writes
  `reference_interactions_clean.csv`. `reference_face_assignment.tsv` is **hand-written**
  (cmyb→cmyb, mll→mll) — it *asserts* the face from the label rather than running stage 9.
- Comparison outputs: `cmyb_candidates_with_control.csv` / `mll_candidates_with_control.csv`
  — that face's filter survivors + its native control appended (`is_control` flag), full 32-col
  metric set. Built by reproducing the `analyze_and_score` merge/filter, then appending the
  control rows (`count`/`priority_score` blank — library-population quantities, N/A for natives).

**Finding:** both natives out-score their surviving pools on binding metrics (c-Myb `helix 0.96`,
`dG −51.4`; MLL `dG −61.3`) but each fails one *library-design* filter — c-Myb on `n_K ≤ 3`
(it has 4 lysines), MLL on `helix_score > 0.70` (only 31% helical; its native motif is a short
helix + turns + polyproline, not a clean helix). Setting aside `count > 1` (N/A), those two
clauses are the sole disqualifiers.

## Literature controls (`control_run/`)

A 24-peptide **benchmark set with measured affinities** — `control_list.csv` — run through the
full pipeline so the scoring can be validated against experiment. 23 MLL-face peptides
(12 literature positives, 11 negatives; Rooklin 2017 JACS, Modell 2022 JACS, thesis Ch. 2–3)
plus `cMyb_native`. Sources and Kd/Ki are columns in `control_list.csv`.

### ⚠️ The NCAAs are the whole point, and the first run did not model them

**`boltz_sequence` is the column the first run used, not `sequence_display`** — a
canonical-amino-acid stand-in for peptides that really carry Bcs, 2meF, CyHex, CyPent, C-IAA,
Aba, nLeu plus an HBS macrocyclic staple.

Those NCAAs mostly **remodel the KIX↔binder interface** — filling pocket space and making
contacts — rather than changing the binder's own fold. The *structural* role (the HBS staple
enforcing helicity) is already stood in for by the idealized α-helical template. So substituting
canonicals throws away precisely the chemistry that distinguishes a tight binder from a weak
one, while keeping the part that was already handled.

The cost is measurable. 21 of 23 MLL controls carry ≥1 NCAA swap, and **9 of 23 collapse onto
just 4 unique Boltz inputs**, including pairs with opposite labels:

| Same Boltz input | Members | Affinity |
|---|---|---|
| `SDGMDFILKNYP` | `HBS-PA-G` (positive) / `HBS-22-G` (negative) | 28 vs 90 µM |
| `SDIMDFVLKNTP` | `MLL1` / `MLL2` (differ only by M4→Bcs) | 1000 vs 500 µM |
| `SDAMDFILKNYP` | `HBS-22-A` / `HBS-22-Aba` / `HBS-PA-A` | 115–180 µM |
| `SDIMDFILKNYP` | `MLL4` / `MLL6` | 32 vs 22 µM |

**Use ResidueX to swap the NCAAs in** — it is a separate tool, not Boltz's built-in
`modifications:`/CCD mechanism. (Boltz 2.2.1 *does* support `modifications: [{position, ccd}]`,
`cyclic: true`, and explicit `bond` constraints, but that is not the route this project uses.)

Run with the **same settings as the full library** (idealized-helix peptide template,
`msa: empty` peptide, precomputed KIX MSA, no pocket constraints) so the numbers are directly
comparable — deliberately *unlike* `boltz_reference_structures/`, which templated its peptides
on the 2AGH crystal chains.

Stage scripts (all scoped to `control_run/`, full-library files untouched):

| Stage | Command |
|---|---|
| 1 | `control_list_to_tsv.py` → `control_list.tsv`, `control_face_truth.tsv`, `control_face_literature.tsv` |
| 2 | `make_helix_schrodinger.py --out_dir control_run/pdb_helices` (Schrödinger) |
| 2b | `pdb_to_cif.py` → `cif_helices/` |
| 3 | `boltz_yaml_gen.py --helix_dir cif_helices --kix-msa kix_msa.csv` |
| 4 | `run_control_boltz.sh` (GPU, ~5 min for 24) |
| 5 | `parse_boltz_results_pipeline_vers.py` → `control_boltz_data.csv` |
| 6 | `run_control_bindcraft.sh` → `bindcraft_out/`, then `parse_bindcraft_results.py` → `control_bindcraft_data.csv` |
| 7 | `cif2pdb.py` → prepwizard (in-session loop, ~50 s each) → `schrodinger_calc_hbond.py` → `hbond_hit_num.py` → `control_interactions_clean.csv` |
| 8 | `check_helix_dssp_skip_first2.py` → `control_dssp_data.csv` |
| 9 | `full_library_face_determination.py` → `control_face_assignment.tsv` |
| 10 | `analyze_controls.py` → `control_scores.csv`, `{cmyb,mll}_candidates_with_controls.csv` |

Two face tables exist on purpose. `control_face_literature.tsv` (from the `face` column) feeds
`hbond_hit_num.py` so each control is counted against the pocket it is *known* to bind — a
negative control that mislocalizes still gets scored on the right face instead of silently
returning `hit_num = 0`. `control_face_assignment.tsv` is the independent geometric stage-9
call, reported as `structural_face_call`. **They agree on all 24.**

`analyze_controls.py` appends the controls to the library's filter survivors rather than
filtering them, so every control gets a score, and records which clauses each *would* have
failed (`failed_filters`, and `failed_binding_filters` which excludes the library-design-only
clauses `count > 1`, `n_K ≤ 3`, `n_M ≤ 3` — see `LIBRARY_DESIGN_FILTERS` in `kix_scoring.py`).
Ranking uses **`priority_score_no_enrichment`**, since `count` is undefined for a control.

**Findings (2026-07-25):**
1. **`cMyb_native` validates the c-Myb arm.** It lands at the 100th percentile of its pool —
   above all 111 library c-Myb candidates (`score 0.915`, `helix 0.91`, `iptm 0.968`,
   `dG −47.6`, `hit_num 4`) and passes every binding filter.
2. **The score does *not* discriminate within the MLL series** — but this is a property of the
   *canonical-substituted* run, not evidence the scoring is broken. Spearman vs measured
   affinity `rho = −0.20, p = 0.39` (n=21, excluding the two 9999 µM sentinels); positives vs
   negatives Mann-Whitney `p = 0.25`.
3. **The signal is below the noise floor, and the noise floor is measurable for free.** The four
   collapsed groups above are *identical Boltz inputs*, so any score difference within a group
   is pure run-to-run noise. Observed spread: **mean 0.112, max 0.169**. The positive-vs-negative
   median gap is **0.050** — 2–3× smaller. `HBS-22-G` and `HBS-PA-G` are the same input and
   landed 0.169 apart (`hit_num` 2 vs 0). Any MLL-face ranking difference under ~0.17 is noise.
   Re-run with ResidueX before drawing conclusions here.
4. **Two of the five scoring terms are near-constant on this set**: `hit_num = 0` for 14/23 MLL
   controls, and `helix_score = 0.80` for 14/23. Each peptide *does* get its own template file
   (24 distinct CIFs, one per sequence) — but all are built with the same idealized α-helix
   backbone geometry, so DSSP largely reads the template back rather than anything
   sequence-specific. Ranking is effectively driven by `protein_iptm` alone.
6. ⚠️ **The ncAA re-run (ResidueX) changed nothing, and the per-metric findings it appeared to
   produce are CONFOUNDED BY PEPTIDE LENGTH.** See `ControlPlan.md` §4 for the full retraction.
   Grafting the real ncAAs left discrimination unchanged (ρ −0.199 → −0.204, both n.s.) — that
   part stands, since it compares the same peptides to themselves. But a per-metric scan
   *appeared* to show `helix_score` correlating with weaker binding (ρ=+0.678, p=0.0007) and
   `hit_num` with tighter binding (ρ=−0.519). **Both evaporate when peptide length is
   controlled** (partial ρ = +0.027 and −0.120). The cause: length vs log(Kd) is ρ=−0.823 and
   length vs helix_score is ρ=−0.813. The set is 14 twelve-mer MLL analogs (median 148 µM) and
   6 thirteen-mer optimised HBS compounds (median 13 µM); the 13-mers bind ~11× tighter *and*
   score lower on DSSP helix **fraction** (bigger denominator, same helical core). Any
   length-tracking metric looks predictive. **Do NOT change the helicity filter or weight on
   this evidence.** Only `interface_interface_hbonds` retains residual signal (partial ρ=−0.372,
   correct sign, n.s. at n=21). **The control set is confounded by design** — tight and weak
   binders are different chemical series of different lengths — so it cannot rank per-metric
   performance. That needs affinity variation *within* a fixed length/series.
7. ⚠️ **`hit_num` measures the wrong interaction type on the MLL face.** Of 89 peptide↔KIX
   H-bonds across the 23 MLL controls, only **15 (17%)** touch the designated MLL hit residues —
   because **4 of the 7 (PHE 27, LEU 43, ILE 75, LEU 79) are hydrophobic and cannot H-bond**.
   The MLL groove binds largely by hydrophobic packing, which `hit_num` does not count, and
   pi-stacking was detected **once across 23 structures**. The observed H-bonds concentrate on
   ARG 86 (23×), ARG 83 (18×), GLU 78 (13×), LYS 82 (11×) — none in the hit list. This is why
   genuine tight binders score `hit_num = 0`. `interface_interface_hbonds` (BindCraft, whole
   interface) is the better-behaved analogue.
8. ⚠️ **The filter set looks non-discriminating only because two clauses cancel.** Binding-filter
   pass rate is positives **17%** (2/12) vs negatives **18%** (2/11) — flat. But drop the
   helicity clause and positives jump to **42%** (5/12) while negatives stay at **18%**; drop
   `hit_num` too and it is 83% vs 73%. Per-clause failures (12 pos / 11 neg):
   `hit_num >= 1` → 5/**8** (helpful direction), `helix_score > 0.70` → **5**/1 (**harmful**),
   `interface_interface_hbonds > 1` → 1/3 (helpful), `confidence_score > 0.85` → 1/0.
   So `helix_score > 0.70` rejects five real binders and one negative — independent evidence for
   dropping helicity, from the filter side rather than the ranking side. Same root cause: the
   five are the longer 13-mers plus `WT_MLL`. **Does not implicate the library** (uniform 10-mers).
   ⚠️ `filter_audit` used to silently mark a clause FAILED when its column was absent, making a
   preview subset look like nothing passed; it now raises. Always pass the full merged metrics.
9. **Score variants added (originals untouched).** `kix_scoring.SCORE_VARIANTS` defines three
   blends sharing a `hit_num`/`unsat`/`protein_iptm` core and differing only in the 4th slot:
   helix (original), nothing (`_no_helix`), or `complex_pde` (`_pde`). Each exists with and
   without the enrichment `count` (`_no_enrichment`). On the controls, ranked in the real pool:
   current AUC 0.598 (p=0.22) → no-helix 0.750 (p=0.023) → pde 0.765 (p=0.017). **Removing
   helicity is the robust gain; `complex_pde` is a smaller further gain within noise at n=23.**
   ⚠️ Evidence is MLL-face only — the c-Myb arm has one control. Note the library is uniform
   10-mers, so `helix_fraction > 0.70` ≡ `helix_count >= 6` there (both select 2,924); the
   length problem is specific to the variable-length control set.
10. **`hit_num ≥ 1` is the harshest clause on real binders** — it fails 14/24 controls, including
   `WT_MLL` (3 µM) and `MLL6` (22 µM). `interface_interface_hbonds > 1` fails 7, `helix_score >
   0.70` fails 6. Only 5/24 pass all binding filters, and they are not the tightest binders.

## Conventions & gotchas

- **`name` is the join key** everywhere. Two extraction conventions exist, watch which you need:
  - Most metric code strips `_model_0.cif` / `.cif` from the Boltz output basename.
  - `full_library_face_determination.py` instead uses the **parent folder name**
    (`Full_Library_Hit_N`) as `name`. Make sure both resolve to the same string when joining.
- Boltz prediction paths look like
  `boltz_out_full/chunk_*/boltz_results_*/predictions/Full_Library_Hit_*/Full_Library_Hit_*_model_0.cif`.
- No formal build/test — validate by running a script on a small input
  (e.g. `tsv_outputs/test200.tsv`) and checking the output CSV.
- Prefer editing the `.py` and its matching `run_*.sh` together (paths/envs are hard-coded in the runners).
- ⚠️ **The env labels in the pipeline table are partly stale.** Verified 2026-07-25:
  `gemmi` and `PyYAML` exist **only in `boltz_env`** — so `pdb_to_cif.py` and `boltz_yaml_gen.py`
  must run under `/scratch/jem9759/envs/boltz_env/bin/python3`, despite `run_pdb_to_cif.sh` and
  `run_boltz_yaml_gen.sh` both activating `general_penv`. `general_penv` has pandas/numpy but no
  Bio; `dssp_env` adds Bio; `scipy` is only in `boltz_env`.
- ⚠️ **BindCraft metrics are STOCHASTIC — this applies to every peptide ever scored, not just
  the controls.** `relax_structure` runs `FastRelax`, a Monte Carlo protocol, with
  `-relax:default_repeats 1`. Measured 2026-07-25 by scoring one unchanged library structure
  (`Full_Library_Hit_0`) three times with the *unpatched* script:

  | run | `binder_score` | `interface_dG` | `interface_sc` |
  |---|---|---|---|
  | 1 | 0.49 | −26.45 | 0.65 |
  | 2 | 2.67 | −21.32 | 0.59 |
  | 3 | 0.43 | −25.67 | 0.72 |
  | *stored library value* | *−5.06* | *−24.94* | *0.77* |

  **`interface_dG` spans 5.1 kcal/mol on identical input** (SD ≈ 2.8 from n=3 — indicative only).
  Consequences: 74% of the 198 survivors clear `interface_dG < -25` by >2 SD and are robust, but
  **33% of all 31k library peptides sit within 1 SD of that cutoff**, so their pass/fail was
  substantially luck. False positives among survivors are limited; **false negatives among the
  ~28k rejected may be substantial**. The physical filter has ~8 BindCraft-derived clauses, each
  independently noisy. This is also part of the run-to-run noise floor measured on the controls.
  Mitigation if it matters: raise `-relax:default_repeats` and/or average several replicates —
  costly across 31k. The stored value sits outside all three fresh runs, hinting the original
  library pass used a different Rosetta build; do not compare fresh numbers to stored ones.
- ⚠️ **`score_peptide.py` lives in `KIX_Project/`**, and the BindCraft repo's copy is named
  `score_peptide_copy.py` (byte-identical). `run_bindcraaft_score.sh` `cd`s into the repo and
  calls a bare `score_peptide.py`, which no longer resolves. If you invoke it by absolute path
  instead, `sys.path[0]` becomes `KIX_Project` and `import functions` breaks — so also set
  `PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft`. See `run_control_bindcraft.sh` for the
  working invocation.
