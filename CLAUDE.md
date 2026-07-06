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
| 1 | xlsx → TSV (`name<TAB>seq`) | `extract_sequences.py`, `full_library_to_tsv.py` | `run.sh` | general |
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

Then **dedup by `Sequence`** (keep highest `protein_iptm`), **split by `face_call`**
(`cmyb` / `mll`), and rank each face independently by `priority_score` — a percentile-rank blend:
`0.2·count + 0.2·hit_num + 0.1·(1−unsat_hbonds) + 0.2·protein_iptm + 0.2·helix_score`.

Outputs `cmyb_candidates.csv` and `mll_candidates.csv` in `full_library_all_metrics/`
(~1.3k candidates each). Preview columns: `name, Sequence, count, hit_num, helix_score,
interface_delta_unsat_hbonds, interface_dG, protein_iptm, priority_score`.

## Key directories

- `full_library_all_metrics/` — per-metric CSVs + final candidate lists (the payoff).
  Note the big layer CSVs (`boltz_metrics_full_library.csv`, `bindcraft_full_library.csv`,
  `dssp_full_library.csv`) and derived (`enrichment_count.csv`, `hbond_hit_counts.csv`,
  `full_library_face_assignment.tsv`).
- `final_metric_outputs/` — where several stage-9/8 runners write first
  (`run_dssp.sh`, `run_full_library_face_determination.sh`). ⚠️ `analyze_and_score` reads from
  `full_library_all_metrics/`, so outputs must be copied/moved there before final scoring.
- `boltz_out_full/` (chunked, `chunk_0000…`), `boltz_outputs/`, `bindcraft_out/`, `bindcraft_relaxed/` — raw tool outputs.
- `yaml_chunks/`, `schrodinger_calc_hbond_chunks/` — split inputs/outputs for SLURM array jobs.
- `tsv_outputs/` — `full_library_out.tsv` is the master library (~31k rows).
- `logs/`, `*.out` — SLURM job logs (named `*_<jobid>.out`).

## Conventions & gotchas

- **`name` is the join key** everywhere. Two extraction conventions exist, watch which you need:
  - Most metric code strips `_model_0.cif` / `.cif` from the Boltz output basename.
  - `full_library_face_determination.py` instead uses the **parent folder name**
    (`Full_Library_Hit_N`) as `name`. Make sure both resolve to the same string when joining.
- Boltz prediction paths look like
  `boltz_out_full/chunk_*/boltz_results_*/predictions/Full_Library_Hit_*/Full_Library_Hit_*_model_0.cif`.
- Not a git repo. No formal build/test — validate by running a script on a small input
  (e.g. `tsv_outputs/test200.tsv`) and checking the output CSV.
- Prefer editing the `.py` and its matching `run_*.sh` together (paths/envs are hard-coded in the runners).
