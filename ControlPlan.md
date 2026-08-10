# Control Set Re-run with ResidueX ncAA Substitution

**Status:** ✅ EXECUTED. All 24 controls grafted with ResidueX and scored. **Written:** 2026-07-25.

---

# SUMMARY FOR REVIEW

*Everything below this section is the working narrative in the order it happened, including one
finding that was later retracted. This summary supersedes it. Read this; use the rest as audit
trail.*

## Question

Does the KIX screen's `priority_score` actually predict binding? Tested against 24 literature
peptides with **measured** Kd/Ki (23 MLL-face: 12 positives / 11 negatives; plus `cMyb_native`).
Sources: Rooklin 2017 JACS, Modell 2022 JACS, thesis Ch. 2–3.

## Method

Idealized helix per peptide (Schrödinger) → cofold with KIX in Boltz-2 (2AGH template,
`msa: empty` peptide, precomputed KIX MSA, no pocket constraints — identical settings to the
31k library) → **graft the real non-canonical residues in with ResidueX** → BindCraft, Schrödinger
H-bond/pi-pi, DSSP, geometric face call → score.

Seven ncAAs (Bcs, 2mF, CyHex, CyPent, C-IAA, Aba, nLeu) across 21 of 24 peptides.

## Findings

**1. Modelling the ncAAs changed essentially nothing.** Spearman vs measured affinity went
−0.199 → −0.204 (both n.s.); the positive/negative gap did not improve. The hypothesis that
canonical stand-ins were erasing the discriminating chemistry is **not supported**.

**2. The current `priority_score` does not separate positives from negatives.** Ranked in the
real pool (87 MLL library survivors + 23 controls):

| score | pos median | neg median | gap | AUC | p |
|---|---|---|---|---|---|
| **current** (incl. helix) | 0.431 | 0.399 | +0.032 | **0.598** | 0.221 |
| helix term removed | 0.507 | 0.378 | +0.129 | 0.750 | **0.023** |
| helix → `complex_pde` | 0.542 | 0.330 | +0.212 | **0.765** | **0.017** |

AUC 0.598 ≈ coin flip. The gap (0.032) is far below the measured noise floor (0.076–0.180, from
peptide pairs that are structurally *identical* after grafting). **Removing the helicity term is
the robust improvement**; substituting `complex_pde` is a smaller further gain that n=23 cannot
firmly resolve (0.750 → 0.765 is roughly one pair swapping order).

**3. The filter set as a whole does not discriminate — because two clauses pull in opposite
directions.** Binding-filter pass rate is positives 17% (2/12) vs negatives 18% (2/11),
indistinguishable. But that flat result masks a real effect:

| filter set | positives | negatives |
|---|---|---|
| all binding filters (incl. helicity) | 2/12 (17%) | 2/11 (18%) |
| **minus helicity** | **5/12 (42%)** | **2/11 (18%)** |
| minus helicity *and* `hit_num` | 10/12 (83%) | 8/11 (73%) |

Removing helicity lifts positives 17% → 42% while negatives stay at 18%. Per-clause failures
(12 positives / 11 negatives):

| clause | pos failed | neg failed | direction |
|---|---|---|---|
| `hit_num >= 1` | 5 | **8** | helpful |
| `helix_score > 0.70` | **5** | 1 | **harmful** |
| `interface_interface_hbonds > 1` | 1 | 3 | helpful |
| `confidence_score > 0.85` | 1 | 0 | neutral |

`helix_score > 0.70` rejects five positives and one negative. This is independent evidence for
dropping the helicity term, arriving from the filter side rather than the ranking side — and it
is the same root cause: the five positives it rejects are the longer 13-mers plus `WT_MLL`.
**This does not implicate the library**, which is uniformly 10-mers with no length variation.

**4. `hit_num` is measuring the wrong thing on the MLL face.** Of 89 peptide↔KIX H-bonds found
across the 23 controls, **only 15 (17%) touch the designated MLL hit residues**, because
**4 of those 7 residues (PHE 27, LEU 43, ILE 75, LEU 79) are hydrophobic and cannot H-bond at
all**. The MLL groove binds largely by hydrophobic packing, which this metric does not count.
Pi-stacking was detected **once in 23 structures**. Meanwhile the real H-bonds concentrate on
ARG 86 (23), ARG 83 (18), GLU 78 (13), LYS 82 (11) — none of which are in the hit-residue list.
This explains why tight binders score `hit_num = 0`. `interface_interface_hbonds` (BindCraft's
count over the *whole* interface) is the better-behaved analogue and had the correct sign.

⚠️ **But the obvious fix makes it WORSE — see "Re-deriving the hit residues" below.** Rebuilding
the list from residues that actually H-bond in 2AGH removes the self-inconsistency and drops the
controls failing `hit_num >= 1` from 13/24 to 1/24 — yet discrimination collapses
(positives-vs-negatives p 0.037 -> 0.400). The stringency was the point: an H-bond inside a
hydrophobic groove is rare and demands a correctly seated peptide, whereas the crystal-derived
surface arginines are salt-bridged by everything. **Keep the original list.**

**5. The c-Myb arm validates.** `cMyb_native` ranks above all 111 library c-Myb candidates
(100th percentile, helix 0.91, iptm 0.968) and passes every binding filter.

## What is verified sound

- **Fold-then-graft ordering** — ncAAs remodel the interface, not the fold; the HBS staple's
  structural role is already carried by the helix template.
- **ResidueX must see the isolated peptide** — three of its functions are chain-blind; verified
  and worked around. Its own `min_distance` cannot discriminate conformers (below threshold for
  *all* of them), so conformers are selected on measured clash/contact against KIX instead.
- **FastRelax resolves graft strain**: worst contacts 2.42 → 3.62 Å, 2.49 → 4.51 Å (vdW ≈3.4 Å).
- **Rosetta ncAA types**: 6 of 7 are database built-ins; Rosetta's own `BCS` is chemically
  S-benzylcysteine, an independent confirmation of that assignment. `CIA` derived from
  `amino-ethyl-cysteine`, charge-corrected to exactly −1.000.
- **No silent residue deletion** — the failure mode was reproduced and closed
  (`nres 98 → 99` with the fix); `--strict-residues` now makes it a loud error.
- **Geometric face calls agree with the literature on all 24.**
- **Regressions pass**: library counts unchanged (30776 → 28029 → 2551 → 198; 111/87), original
  score columns bit-identical, canonical `hit_num` reproduces exactly.
- **The library is unaffected by the helix-fraction/length issue** — all 31,392 library peptides
  are 10-mers, so `fraction > 0.70` and `count >= 6` select the identical 2,924. The length
  problem is specific to the variable-length control set.

## Known limitations (please weigh these)

1. **The HBS macrocyclic staple is not modelled.** Only side chains are grafted.
2. **10 of 23 MLL controls are linear, not stapled** (`MLL1–6` series, `WT_MLL` — Ac-capped),
   yet all received an idealized α-helix template. Free 12-mers are largely disordered in
   solution, so their folds are biased toward helix by construction.
3. **Boltz confidence metrics are computed pre-graft.** `complex_pde`, `protein_iptm`, `pep_ptm`
   describe the canonical stand-in structure — so our best-performing term never sees the ncAAs.
   Also, KIX complexes (2AGH etc.) are in the PDB, so `complex_pde` may partly reflect
   recognition of known binders rather than physics.
4. **BindCraft is stochastic**: `interface_dG` spans 5.1 kcal/mol on *identical* input
   (SD ≈2.8, n=3). ~33% of the library sits within 1 SD of the `< -25` cutoff.
5. **The control set is confounded by design**: tight binders are the longer, later-optimised
   Modell-2022 series; weak ones are the earlier 12-mer MLL analogs. Length vs log(Kd)
   ρ = −0.823.
6. **Two pairs remain structurally degenerate** after grafting (`HBS-PA-G`/`HBS-22-G`,
   `HBS-22-A`/`HBS-PA-A`) — they differ only by the X vs Z cap, which is not in `ncaa_swaps`.
   These give the noise-floor estimate.
7. **n=23, one face.** The c-Myb arm has a single control, so nothing is validated there.
8. **Conformer choice is a greedy heuristic** (most KIX contacts among clash-free), then relaxed.

## Re-deriving the hit residues from 2AGH — tested, and the original list wins

`derive_hit_residues.py` runs the same Schrodinger interaction workflow on 2AGH model 1
(chains: A = c-Myb peptide, B = KIX, C = MLL peptide; numbering offset 585, verified against
KIX_SEQUENCE) and reports H-bond, pi-pi and heavy-atom contacts per KIX residue ->
`hit_residue_derivation.tsv`. Three named sets now live in `kix_scoring.HIT_RESIDUE_SETS`,
selectable via `hbond_hit_num.py --residue-set`.

**Diagnosis confirmed.** Every residue in the current lists *does* contact the native peptide —
the ChimeraX derivation was sound. But the lists are **contact**-derived while `hit_num` counts
**H-bonds/pi-pi**, so 6 of 8 (c-Myb) and 5 of 7 (MLL) listed residues make zero H-bond/pi-pi in
the crystal and can never contribute.

| set | c-Myb | MLL |
|---|---|---|
| original | 14,18,21,65,69,72,73,76 | 27,39,43,46,71,75,79 |
| crystal | 9,21,61,76,80,81 | 27,39,83,84,86 |
| union | 12 residues | 10 residues |

GLU 78 / LYS 82 excluded by decision: frequent partners in our control structures (13x, 11x) but
5.63 A / 6.99 A from the native peptide, so including them would fit the list to the predictions.

**Result — the fix works mechanically but LOSES the signal:**

| set | controls failing `hit_num>=1` | pos mean | neg mean | MWU p | rho vs affinity | library survivors |
|---|---|---|---|---|---|---|
| **original** | 13/24 | 1.00 | 0.27 | **0.037** | **-0.519** | 198 (111/87) |
| crystal | 1/24 | 2.33 | 2.09 | 0.400 | -0.031 | 209 (107/102) |
| union | 0/24 | 2.75 | 2.09 | 0.145 | -0.167 | 220 (111/109) |

Full-score AUC moves the same way: current-score 0.598 -> 0.523 -> 0.553; `_pde` 0.765 -> 0.712
-> 0.750. **`original` is best on every discrimination measure.**

**Why the "flaw" is the feature.** The crystal residues are surface arginines (ARG 39/83/84/86);
any peptide with an Asp or Glu nearby forms a salt bridge, so *everything* scores 2-3 and the
metric stops discriminating. The original list is the hydrophobic-lined groove, where an H-bond
is rare and demands the peptide be correctly seated — so most peptides score 0 and a non-zero
score is meaningful. `hit_num = 0` for `WT_MLL` is the metric correctly reporting that the
predicted pose puts no polar group in the groove, not a bug.

**Decision: keep `original` as the default** (it already is). `crystal`/`union` are retained as
selectable alternatives so this is reproducible. Caveat: the original's rho = -0.519 is largely
length-confounded (partial rho = -0.120).

## 7. `hit_num_v2` — count hydrophobic contacts too (BEST RESULT SO FAR)

`hit_num` counts only H-bonds/pi-pi, so the 4 apolar MLL hit residues (PHE 27, LEU 43, ILE 75,
LEU 79) can never contribute — the hydrophobic packing that drives MLL binding is invisible.
Rather than change the residue list (tested in §6 and rejected), keep it and add the missing
interaction type:

```
hit_num_v2 = (H-bonds + pi-pi on hit residues)              # = hit_num, unchanged
           + (# APOLAR hit residues within 4.5 A)           # new
```

Each residue contributes via exactly one mechanism, so nothing is double counted. Apolar
membership is derived from `KIX_SEQUENCE` (`kix_scoring.apolar_hit_residues`), not hardcoded.
Contacts come from `hydrophobic_contacts.py` (mirrors `full_library_face_determination.py`;
2 min for all 31,392 library structures, 0 errors).

**Only works with the `original` residue list** — `crystal` has just 1 apolar residue left, and
`union` dilutes the signal with generic surface arginines:

| list | apolar residues | MWU p | rho vs affinity |
|---|---|---|---|
| **original** | 4 (27,43,75,79) | **0.007** | **-0.622** |
| crystal | 1 (27) | 0.400 | -0.031 |
| union | 4 | 0.075 | -0.225 |

Cutoff is insensitive: 4.5 A and 5.0 A give identical results; 4.0 A is slightly worse
(p=0.016).

**Effect on the pooled score (87 MLL library survivors + 23 controls) — improves every blend:**

| score | gap | AUC | p |
|---|---|---|---|
| current (`hit_num`, helix) | +0.032 | 0.598 | 0.221 |
| `hit_num_v2`, helix | +0.073 | 0.636 | 0.141 |
| `hit_num`, no helix | +0.129 | 0.750 | 0.023 |
| `hit_num_v2`, no helix | +0.180 | 0.773 | 0.014 |
| `hit_num`, pde | +0.212 | 0.765 | 0.017 |
| **`hit_num_v2`, pde** | **+0.234** | **0.788** | **0.011** |

The two changes are independent and additive. All six now exist as columns
(`priority_score[_v2][_no_helix|_pde][_no_enrichment]`); **`hit_num` and the original score are
untouched**, and the 198-candidate shortlist is unchanged.

⚠️ **Do not put `hit_num_v2` in the FILTER.** It is never 0 in the library (30,768 of 30,776
would pass `>= 1`), so it would silently disable that clause. It is a ranking term only.

⚠️ **The apolar term saturates**: on the controls it takes 2 values (3 or 4 of 4 contacted),
21 of 23 scoring 4. It works in combination but only p=0.074 alone. A continuous form (buried
surface area per hit residue) is the obvious follow-up. It is also the most length-correlated
metric tested (rho_len 0.64) — though the library is uniform 10-mers, where that contributes
nothing to ranking.

## 8. Dual-face cofold — the binary face call FAILS, but per-copy interface_dG works

KIX + 2 copies of each peptide (chains A/B/C), no constraints; top 40 per face by
`priority_score_v2_pde`, plus `cMyb_native`/`WT_MLL` and 4 decoys.

**The binary face call is uninformative.** 64/82 came out DUAL — but so did both single-face
natives, both measured non-binders (>9999 uM) and a scrambled sequence. Give Boltz two empty
grooves and two chains and it fills both. `pair_chains_iptm` was no better: it ranks WT_MLL #1,
cMyb_native #2 and the NON-BINDER MLL6_YA #3 of 86.

**Running the full metric pipeline per copy fixes it.** `split_dual_copies.py` splits each
3-chain complex into KIX+copyB and KIX+copyC (chain C renamed B), so the standard unmodified
pipeline runs on each — and it sidesteps `sanitize_structure`'s C/D deletion. 172 sub-complexes
through BindCraft, Schrodinger hbond/pi-pi, hydrophobic contacts, DSSP and the face call.

Ranked on the **weaker** copy (a real dual binder needs BOTH sites good), of 86:

| control | truth | by `weaker_dG` | by composite |
|---|---|---|---|
| cMyb_native | native | **#1** | #18 |
| WT_MLL | native | #3 | #3 |
| MLL6_YA | **non-binder** | #50 | **#20** |
| MLL6_CstarA | **non-binder** | #55 | #47 |
| decoy_scramble | nonsense | #75 | #84 |
| decoy_polyAla | nonsense | #86 | #85 |

`interface_dG` is the **only metric all study** that orders the controls correctly. Not a length
artefact: MLL6_YA is a 12-mer (longer than the 10-mer candidates, so more contacts available)
yet ranks below 49 of them. Length-matched candidates vs decoys (both 10-mers):
`weaker_dG` p=0.011, composite p=0.010, AUC 0.988.

⚠️ **`interface_dG` IS NOT IN ANY COMPOSITE SCORE.** The blend is
`0.29*hit_num_v2 + 0.14*(1-unsat) + 0.29*protein_iptm + 0.29*(1-complex_pde)` — 58% Boltz
confidence, the family that ranks non-binders top. `interface_dG` is used only as a filter
(`< -25`). That is why the composite is *worse* than dG alone here, and it is the clearest
candidate for a scoring change: **promote interface_dG from filter to ranking term.**

**Shortlist**: 62/80 put their copies on different faces; 40 also have both interfaces
`dG < -25` and beat both decoys -> `dual_cofold/dual_shortlist.csv`. Top by weaker interface:
Hit_3019 (-37.5), Hit_2736 (-35.9), Hit_1805 (-35.4/-35.5, notably symmetric), Hit_383 (-35.1),
Hit_7965 (-34.5) — versus cMyb_native -43.4 and WT_MLL -36.1. These come from the MIDDLE of the
per-face rankings (#16-39), so dual capability is a distinct property from single-face strength.

**Caveats**: 6 control points only; BindCraft's +/-2.8 kcal/mol noise means the top ~5 should not
be finely ranked against each other; and this measures each copy as placed in the 3-chain
complex, not whether it would bind alone.

## Decisions to discuss

- Drop the helicity term from the MLL ranking? (clear improvement, AUC 0.598 → 0.750)
- Adopt `complex_pde` in its place? (further gain but within noise at this n)
- Redefine the MLL hit-residue set, or replace `hit_num` with an interface-wide H-bond count?
- Is the idealized-helix template appropriate for the linear MLL analogs?
- Average BindCraft over replicates to cut the ±2.8 kcal/mol noise?

---

# ⚠️ WORKING NARRATIVE (superseded by the summary above)

## 1. Modelling the ncAAs did NOT improve discrimination

| run | score | ρ vs log(affinity) | p | pos−neg gap | p |
|---|---|---|---|---|---|
| canonical | `no_enrichment` | −0.199 | 0.386 | 0.050 | 0.249 |
| **ncAA-grafted** | `no_enrichment` | **−0.204** | 0.375 | 0.032 | 0.221 |
| **ncAA-grafted** | `swap_sensitive` | **−0.206** | 0.371 | 0.020 | 0.186 |

Essentially unchanged, and the positive-vs-negative gap *shrank*. The noise floor from the two
still-degenerate pairs is 0.076 and **0.180** — 4–9× the signal. **The hypothesis that canonical
substitution was what erased the discriminating chemistry is NOT supported.**

## 2. But six individual metrics DO track affinity — the composite is cancelling them

Spearman vs log10(Kd/Ki), n=21 MLL controls, ncAA run. Bonferroni threshold p<0.0036:

| metric | ρ | p | sign | verdict |
|---|---|---|---|---|
| `complex_pde` | +0.685 | 0.0006 | **correct** | survives Bonferroni |
| `helix_score` | +0.678 | 0.0007 | **BACKWARDS** | survives Bonferroni |
| `pep_ptm` | +0.573 | 0.0066 | backwards | nominal |
| `hit_num` | −0.519 | 0.0160 | correct | nominal |
| `interface_interface_hbonds` | −0.469 | 0.0319 | correct | nominal |
| `interface_sc` | +0.437 | 0.0477 | backwards | nominal |
| *`priority_score_no_enrichment`* | *−0.204* | *0.375* | — | **not significant** |

## 3. ⚠️ `helix_score > 0.70` is BACKWARDS on the MLL face

```
HBS_II  0.48uM 0.727 | WT_MLL 3uM 0.529 | HBS_VI 7.4uM 0.636 | HBS-W 13uM 0.727
HBS_IV   13uM 0.636 | HBS_III 17uM 0.636 | HBS_V 21uM 0.727
---- every peptide at helix 0.800 is >=22uM, incl. ALL negatives ----
MLL6 22uM 0.800 | HBS-PA-G 28uM 0.800 | MLL4 32uM 0.800 | MLL3 58uM 0.800 | ...
```

The pipeline **filters** `helix_score > 0.70` and **rewards** helicity at weight 0.2. On the MLL
face that rejects WT_MLL, HBS_VI, HBS_IV and HBS_III — four of the tightest binders — and keeps
every weak/negative peptide. It also explains why the composite fails: the backwards
`helix_score` term (0.2) directly cancels the correct `hit_num` term (0.2).

Chemically coherent, and **MLL-specific**: the native MLL motif is a short helix + turns +
polyproline (the native MLL reference is 31% helical), so high helicity selects *against* the
native binding mode. c-Myb is the opposite — `cMyb_native` is 0.913 helical and binds at 1.9 µM,
and it still ranks above all 111 library c-Myb candidates.

## 4. ⛔ RETRACTION — sections 2 and 3 above are CONFOUNDED BY PEPTIDE LENGTH

**Do not act on the "helix_score is backwards" conclusion.** A confound check after the fact
shows peptide length explains essentially all of it:

```
length vs log10(Kd)        rho = -0.823  p < 0.0001     <-- length IS the signal
length vs helix_score      rho = -0.813  p < 0.0001     <-- helix_score is a length proxy
```

The set is 14 twelve-mers (the MLL analog series, median 148 uM) and 6 thirteen-mers (the
optimised Modell-2022 HBS series, median 13 uM). The 13-mers bind ~11x tighter *and* score lower
on DSSP helix **fraction**, because the fraction's denominator grows with length while the
helical core does not. Any metric tracking length therefore looks predictive.

Partial Spearman controlling for length:

| metric | raw rho | partial rho (given length) | verdict |
|---|---|---|---|
| `helix_score` | +0.678 | **+0.027** | **vanishes — the section-3 claim is dead** |
| `hit_num` | −0.519 | **−0.120** | vanishes |
| `complex_pde` | +0.685 | +0.304 | much weaker, n.s. |
| `pep_ptm` | +0.573 | +0.139 | vanishes |
| `interface_sc` | +0.437 | +0.075 | vanishes |
| `interface_interface_hbonds` | −0.469 | **−0.372** | best survivor, correct sign, still n.s. (n=21) |
| `protein_iptm` | −0.321 | +0.094 | vanishes |
| `interface_dG` | −0.290 | +0.001 | vanishes |

Within-length-group tests agree: `helix_score` vs affinity is rho=+0.05 (n=14, 12-mers) and
rho=0.00 (n=6, 13-mers) — **no effect at all** once series is held constant.

**Corrected conclusions:**
1. Modelling the ncAAs did not improve discrimination (section 1 stands — that comparison is
   internal to the same peptides, so length cancels).
2. There is **no evidence** that `helix_score` is backwards, and **no basis** for changing the
   helicity filter or weight. Doing so on this evidence would be a mistake.
3. `hit_num` is **not** demonstrated to be predictive either; its apparent signal was length.
4. The only term with residual signal after controlling for length is
   `interface_interface_hbonds` (rho −0.372, correct sign), and at n=21 that is not significant.
5. **This control set cannot resolve per-metric performance.** It is confounded by design: the
   tight binders are a different, longer, separately-optimised chemical series from the weak
   ones. Deconfounding needs either affinity variation *within* a fixed length/series, or many
   more peptides.

**Caveats:** n=21, MLL face only (the c-Myb arm has one control); `helix_score` takes 5 distinct
values; BindCraft contributes ±2.8 kcal/mol run-to-run noise; the HBS staple is still unmodelled.

## 5. ✅ RESOLUTION — use helix COUNT for QC, drop helix FRACTION from the score

The retraction in §4 does not mean helicity is useless — it means the *fraction* is the wrong
statistic. Helicity legitimately serves as a **QC gate** ("the HBS peptides are stapled, they
had better come out helical"), which is separate from its role as a ranking term.

| metric | vs affinity | **vs length** |
|---|---|---|
| `chain_b_helix_fraction` | ρ=+0.678, p=0.0007 | ρ=**−0.824**, p<0.0001 |
| `chain_b_helix_count` | ρ=+0.104, p=0.65 | ρ=**−0.115**, p=0.60 |

```
length 12 (n=16):  helical residues 7-8    fraction 0.70 -0.80
length 13 (n=6):   helical residues 7-8    fraction 0.636-0.727
length 19 (n=1):   helical residues 9      fraction 0.529
```

**Every peptide has 7–9 helical residues regardless of length.** The helical core is constant —
as expected from an idealized-helix template plus an HBS staple — and the fraction moves only
because its denominator does. `helix_fraction` was measuring length, not helicity.

**Recommendation:**
- **QC gate → `chain_b_helix_count` (e.g. ≥6).** Length-robust, preserves the sanity check, and
  passes all 24 controls — whereas `fraction > 0.70` rejects 6, including WT_MLL, HBS_VI,
  HBS_IV and HBS_III.
- **Ranking score → drop the `helix_fraction` term.**

Positive-vs-negative separation, MLL controls (AUC 0.5 = none; noise floor ≈0.18):

| score | full set (n=23) gap / AUC / p | 12-mers only (n=14) gap / AUC / p |
|---|---|---|
| A current (with helix) | 0.071 / 0.553 / 0.345 | 0.012 / 0.473 / 0.590 |
| B helix term removed | 0.167 / 0.731 / **0.032** | 0.070 / 0.536 / 0.432 |
| C helix → `complex_pde` | **0.298 / 0.792 / 0.010** | 0.124 / 0.627 / 0.231 |

C is the only variant whose gap clears the noise floor, but **length-matched it is
underpowered (5 pos vs 11 neg) and n.s.** — promising, not established. `complex_pde` also
retains the most residual signal after controlling for length (partial ρ=+0.304).

---

## Context

`control_list.csv` is a 24-peptide benchmark with **measured affinities** (12 MLL positives,
11 MLL negatives, 1 c-Myb native), intended to validate whether the KIX screen's scoring
actually predicts binding.

A first pass ran all 24 through the full pipeline using the `boltz_sequence` column — a
**canonical-amino-acid stand-in** for peptides that really carry Bcs, 2mF, CyHex, CyPent,
C-IAA, Aba and nLeu. That was the wrong call. Those ncAAs are what remodel the KIX↔binder
interface (filling pocket space, making contacts); the *structural* role of the HBS staple is
already stood in for by the idealized α-helical template. Substituting canonicals therefore
erased exactly the chemistry that distinguishes a tight binder from a weak one.

This plan re-runs the set correctly: **build idealized helices → cofold in Boltz-2 with the
helix + 2AGH → swap ncAAs in with ResidueX → then score.**

---

## Part 1 — What the first pass already did (session recap)

Work completed and on disk, all self-contained under `control_run/`:

| Stage | Output | Status |
|---|---|---|
| 1 Build TSV from `control_list.csv` | `control_list.tsv`, `control_face_truth.tsv`, `control_face_literature.tsv` | done |
| 2 Idealized helices (Schrödinger) | `pdb_helices/` 24 PDBs | done, **reusable** |
| 2b PDB→CIF | `cif_helices/` 24 CIFs | done, **reusable** |
| 3 Boltz YAMLs | `yamls/` 24 | done, **reusable** |
| 4 Boltz-2 cofold (GPU, 4m43s) | `boltz_out/` 24 predictions | done, **reusable — this is the ResidueX input** |
| 5 Parse Boltz confidences | `control_boltz_data.csv` | done, still valid |
| 6 BindCraft interface scoring | `control_bindcraft_data.csv` | **must redo after swap** |
| 7 Schrödinger prepwizard + H-bond/pi-pi | `control_interactions_clean.csv` | **must redo after swap** |
| 8 DSSP helicity | `control_dssp_data.csv` | redo (cheap), expected ~unchanged |
| 9 Geometric face call | `control_face_assignment.tsv` | redo (cheap), expected ~unchanged |
| 10 Scoring/comparison | `control_scores.csv`, `{cmyb,mll}_candidates_with_controls.csv` | redo |

**Boltz does not need re-running.** The cofold uses the canonical stand-in sequence by design;
ResidueX grafts onto that finished structure.

Scripts created: `control_list_to_tsv.py`, `analyze_controls.py`, `kix_scoring.py`,
`run_control_boltz.sh`, `run_control_bindcraft.sh`.
Scripts modified: `analyze_and_score_all_metrics.py`, `hbond_hit_num.py` (now takes CLI paths),
`make_helix_schrodinger.py` (added `--out_dir`), `CLAUDE.md`.

**Scoring change delivered:** `priority_score_no_enrichment` — the same percentile blend as
`priority_score` with the `count` term dropped and the remaining four weights renormalized
(÷0.7). Lets controls (which have no library enrichment count) rank on the library's scale.
Thresholds + scoring functions were pulled into `kix_scoring.py`, imported by both the library
scoring and the control comparison so they cannot drift. Refactor verified to reproduce library
numbers exactly (30776 → 28029 → 2551 → 198; 111 c-Myb / 87 MLL).

**Results from the canonical run:**
- `cMyb_native` ranks above **all 111** library c-Myb candidates (score 0.915). c-Myb arm validates.
- MLL series shows **no discrimination**: Spearman vs affinity `rho=−0.20, p=0.39`;
  positives-vs-negatives `p=0.25`.
- **Noise floor measured for free**: 9 of 23 MLL controls collapsed onto 4 *identical* Boltz
  inputs, so within-group score spread is pure noise → **mean 0.112, max 0.169**. The
  positive-vs-negative median gap is **0.050**. Signal is 2–3× below noise.

Two environment traps were hit and are now documented in `CLAUDE.md`: `gemmi`/`PyYAML`/`scipy`
live only in `boltz_env` despite the runners activating `general_penv`; and `score_peptide.py`
moved to `KIX_Project/` (the BindCraft copy is `score_peptide_copy.py`), so it needs
`PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft` when invoked by absolute path.

---

## Part 2 — Findings that shape this plan

**1. ResidueX API** (github.com/XDaiNYU/ResidueX): three functions —
`split_pdb_by_residue`, `NCAA_sdf_generation`, `integrate_NCAA_into_peptide`.
Input = PDB + numeric `residue_id` + **SMILES**. Output = modified PDB, ncAA SDF, distance file.
It grafts onto an existing structure, confirming the fold-then-swap ordering.

**2. Residue numbering verified.** In the cofolded PDB, chain B is numbered 1..N matching
`boltz_sequence` exactly (checked `HBS-PA-G`: SER1 ASP2 GLY3 MET4 ASP5 PHE6 … PRO12).
Chain A (KIX) is 1..87. So `ncaa_swaps` positions map **directly** to ResidueX `residue_id`
with no offset. The `-NH2`/`-OH`/`Ac-` termini and the X/Z cap are absent from chain B.

**3. Seven distinct ncAAs, each at a fixed position:**

| ncAA | uses | position | parent |
|---|---|---|---|
| 2mF | 17 | 6 | PHE |
| Bcs | 16 | 4 | MET (stand-in) |
| CyHex | 3 | 4 | PHE |
| C-IAA | 3 | 12 | GLU |
| CyPent | 2 | 11 | PHE |
| Aba | 1 | 3 | ALA |
| nLeu | 1 | 3 | LEU |

`WT_MLL`, `MLL1`, `cMyb_native` need no swap.

**4. ⚠️ BindCraft silently deletes ncAAs.** `score_peptide.py:49-59` calls `pyrosetta.init`
with `-ignore_unrecognized_res` **and** `-mute all`. A residue with no Rosetta params is
**deleted from the pose with no visible warning**, then `FastRelax` (`:182`) relaxes the
truncated peptide and `pose.dump_pdb` (`:196`) writes it out as the scoring structure. Every
interface metric is then computed on a peptide with a hole where the ncAA was. Nothing in
Python strips it — `clean_pdb` (`generic_utils.py:303-310`) is a line-type filter that
explicitly *keeps* HETATM. There is **no `-extra_res_fa` hook anywhere**. Also note
`-corrections::beta_nov16 true` (`:54`) — params must use beta_nov16 atom typing. A mid-chain
deletion additionally fragments chain B, and the dict comprehension at `:308-311` keyed by
chain letter silently keeps only the last fragment.

**5. Two scoring terms are frozen across the swap.** Boltz confidences (`protein_iptm`,
`pep_ptm`, `complex_pde`) are computed during folding, *before* the swap. `helix_score` reads
the backbone, which ResidueX preserves. In `priority_score_no_enrichment` those are 0.4 of the
0.7 total weight (~57%) and cannot respond to the substitution — only `hit_num` and
`unsat_hbonds` react. → **Add a swap-sensitive score variant.**

**6. ⚠️ The swap resolves 5 of 9 degeneracies, not all 9.** Identity after swap =
`boltz_sequence` + `ncaa_swaps`; unique structures go 19 → 22 of 24. Two pairs remain identical
because they differ **only by the HBS cap (X vs Z), which is not recorded in `ncaa_swaps`**:

| Still identical | Members | Affinity |
|---|---|---|
| `SDGMDFILKNYP + M4->Bcs; F6->2mF` | `HBS-PA-G` (X, **positive**) / `HBS-22-G` (Z, **negative**) | 28 vs 90 µM |
| `SDAMDFILKNYP + M4->Bcs; F6->2mF` | `HBS-22-A` (Z) / `HBS-PA-A` (X) | 115 vs 180 µM |

A positive/negative pair the pipeline still cannot separate. Treat as a **built-in noise
control** for the re-run (they give a fresh noise-floor estimate), and flag the X/Z cap as an
open modelling question.

**7. Environment: ResidueX is not installed.** No clone anywhere; no `residuex` conda env.
No env has the full stack — `openbabel` and `rmsd` are missing from *all* envs; `boltz_env`
is closest (has rdkit + biopython) but shouldn't be perturbed. `~/.condarc` has **only the
`defaults` channel**, so rdkit/openbabel/py3Dmol need explicit `-c conda-forge`. No mamba.
Quota: 5 TB / 5 M inodes, currently 1.1% space but **21% inodes** — conda envs are inode-heavy,
so build in `/scratch/jem9759/envs/`, never `$HOME` (30k inode cap).

---

## Part 3 — Step-by-step plan

### Step 1 — Install ResidueX
```
git clone https://github.com/XDaiNYU/ResidueX /scratch/jem9759/ZhangWork/ResidueX
conda create -p /scratch/jem9759/envs/residuex python=3.10
conda install -p /scratch/jem9759/envs/residuex -c conda-forge rdkit openbabel py3Dmol biopython
pip install rmsd && pip install -e /scratch/jem9759/ZhangWork/ResidueX
```
Prefer the repo's `minimal_environment.yml` if it solves; fall back to the explicit list above.
Smoke-test the three imports before proceeding.

### Step 2 — ncAA SMILES table (**blocked on you**)

**ResidueX needs the FULL amino acid SMILES, not just the side chain.** Confirmed from source
(`ResidueX/residuex.py`): it matches `Chem.MolFromSmarts('NCC(=O)')` on the ncAA and derives the
side chain by *subtracting* the backbone match; on the peptide side it keeps only
`allowed_atoms = ['N', 'H', 'CA', 'C', 'O']`. The README example is a **capped** residue —
`CN[C@@H](CC1=CN(C)C=N1)C(=O)C` — i.e. N-methyl + methyl ketone, not free NH2/COOH.

Template: `CN[C@@H](<side chain>)C(=O)C`  (`[C@@H]` = L in this atom ordering).

**RESOLVED — table written to `control_run/ncaa_smiles.csv`, all 7 RDKit-validated.**

| ncAA | full SMILES | code | pos | source |
|---|---|---|---|---|
| Bcs | `CN[C@@H](CSCc1ccccc1)C(C)=O` | BCS | 4 | **user-confirmed** — S-benzyl-Cys |
| C-IAA | `CN[C@@H](CSCC(=O)O)C(C)=O` | CIA | 12 | **user-confirmed** — S-carboxymethyl-Cys (acid form) |
| 2mF | `CN[C@@H](Cc1ccccc1C)C(C)=O` | M2F | 6 | drafted — 2-methyl-Phe |
| CyHex | `CN[C@@H](CC1CCCCC1)C(C)=O` | CHX | 4 | drafted — cyclohexyl-Ala |
| CyPent | `CN[C@@H](CC1CCCC1)C(C)=O` | CPN | 11 | drafted — cyclopentyl-Ala |
| Aba | `CN[C@@H](CC)C(C)=O` | ABA | 3 | drafted — 2-aminobutyric acid |
| nLeu | `CN[C@@H](CCCC)C(C)=O` | NLE | 3 | drafted — norleucine |

Validation run (RDKit): all 7 parse, and the `NCC(=O)` backbone SMARTS matches **exactly once**
in each — so no ambiguous graft. Formulas/MW sane.

⚠️ **CIP descriptors differ but stereochemistry is uniform.** Bcs and C-IAA report as *R*, the
other five as *S*. This is a CIP-priority artifact of the β-sulfur (same reason L-cysteine is
*R* while most L-amino acids are *S*), **not** an error — all seven share the same `[C@@H]`
spatial arrangement and are L. Do not "fix" this.

Still worth a glance before the run: that **2mF is the *ortho* isomer** (assumed, not confirmed).
`BCS`, `ABA` and `NLE` happen to match real PDB CCD codes for these residues, which is a
convenient cross-check.

### Step 3 — Build the swap table
New script `build_ncaa_swap_table.py`: parse the `ncaa_swaps` column of `control_list.csv`
(regex `^([A-Z])(\d+)->(.+)$`), join to `ncaa_smiles.csv`, and emit
`control_run/ncaa_swap_table.csv` = `name, residue_id, parent_aa, ncaa_name, smiles`.
**Assert** the parent one-letter code matches the actual residue at that position in the
cofolded PDB — this catches any numbering drift loudly rather than grafting onto the wrong site.

### Step 4 — Run ResidueX on the existing cofolds

Read of `ResidueX/residuex.py` (cloned to `/scratch/jem9759/ZhangWork/ResidueX`) surfaced four
constraints that dictate the design here. **Do not pass the KIX+peptide complex to ResidueX.**

**(a) ⚠️ ResidueX must be run on the ISOLATED PEPTIDE (chain B), never the complex.** Three
places are chain-blind and would silently corrupt KIX:
- `split_pdb_by_residue` selects on `residue.id[1] == residue_id` with no chain filter (`:26`)
  — matches residue 4 of chain A *and* chain B.
- `filter_residue_atoms` filters on `int(line[22:26])` only (`:156`) — would strip KIX residue 4's
  side chain too.
- **Worst:** `get_pep_ready_carbon_alpha` does `matches = pep_ready.GetSubstructMatches(Chem.MolFromSmiles('NCC(=O)'))`
  then `matches[residue_id-1]` (`:183-188`) — it indexes the *N-th backbone match in the whole
  molecule*. With KIX present, `matches[3]` is KIX residue 4, not peptide residue 4.

  The package's own example confirms this usage: `example/6ox2_Z/` ships `ranked_100_sp_pep.pdb`
  and `ranked_100_sp_pro.pdb` as **separate** peptide/protein files.

  → Split chain B out, swap, then recombine with the untouched chain A from the original cofold.

**(b) ⚠️ Multi-swap is not supported upstream and must be chained by us.** The tutorial and
`ResidueX_example.py` demonstrate exactly one swap (`residue_id = 8`). Our load:

| swaps needed | peptides |
|---|---|
| 1 | 6 |
| 2 | 9 |
| 3 | 5 |
| 4 | 1 (`HBS_II`) |

15 of 21 need ≥2. Each swap round-trips SDF→PDB via `obabel`, which mangles residue identity,
and the next round depends on `Chem.MolFromPDBFile` still yielding backbone matches **in residue
order**. This is the single biggest technical risk in the whole plan — validate the intermediate
after every swap (residue count, backbone-match count, and that already-swapped positions survive).

**(c) `min_distance` ignores KIX.** It is computed between ncAA side-chain atoms and the
*peptide* only (`:387-395`), i.e. intra-peptide clash. Since the entire point is the ncAA
filling the KIX pocket, do our own clash/contact check **after recombination**: generate all
conformers, integrate each, recombine with KIX, then select on (i) no heavy-atom clash with KIX
(min distance > ~2.2 Å) and (ii) best pocket contact. Keep the selection criterion in one place
and log the rejected conformers.

**(d) Output annotation needs fixing — but only three fields.** `integrate_NCAA_into_peptide`
writes an SDF; the example converts it with `obabel -isdf ... -O ....pdb`. The result **is** a
valid PDB, but the grafted residue comes back as:

```
HETATM  127  N   UNK A   8      28.922  33.712  17.828     <- was: ATOM ... GLU B 7
```

Three fields differ from the input: record type `ATOM`→`HETATM`, residue name →`UNK`
(Open Babel's "unknown"; `rename_UNL_pdb_file` at `:435` uses `UNL`), and chain `B`→`A`.

✅ **Residue numbering survives intact** (verified on the shipped example: 1,2,…,7,`UNK 8`,9,…).
So this is a text rewrite of three fields, not a reconstruction — and it means the multi-swap
chaining in (b) is far lower risk than feared, since it depends on numbering surviving each
round-trip.

Rewrite to the real 3-letter codes from `ncaa_smiles.csv` (BCS/CIA/M2F/CHX/CPN/ABA/NLE) so the
Step-5 Rosetta params match, and restore chain `B`. Downstream consumers that key on these:
Rosetta (residue type by name), DSSP (only reads recognized `ATOM` residues), and both
`schrodinger_calc_hbond.py` and `full_library_face_determination.py` (select `chain B`).

**(e) A/B test: Met vs Phe as the Bcs stand-in — TESTED, REJECTED (n=3).**
Bcs (CH2-S-CH2-phenyl) is the one ncAA whose canonical stand-in is a clear shape mismatch
(Met = short linear thioether, no ring), and it is used by 16 of 23 MLL controls. Hypothesis:
folding with Phe at position 4 gives an aromatic-shaped pocket, so the benzyl graft fits better.
Tested by refolding MLL2/MLL3/MLL4 with `F` at position 4 (`control_run/ab_test_standin/`).

| peptide | viable MET | viable PHE | best KIX dist MET | PHE | contacts MET | PHE | Δiptm |
|---|---|---|---|---|---|---|---|
| MLL2 | 1/29 | 4/26 | 2.25 | **2.41** | 61 | 65 | **+0.024** |
| MLL3 | 2/31 | 6/28 | **2.32** | 2.24 | 60 | 61 | −0.013 |
| MLL4 | 6/28 | 6/26 | 2.34 | 2.33 | 47 | 47 | −0.019 |

**Verdict: keep Met; do not switch.** On n=1 (MLL2) Phe looked clearly better on every axis, but
it did not replicate — the iptm/pde gain reversed on MLL3 and MLL4, and selected-structure
clearance and contacts are a wash. The only consistent effect is *more viable conformers to
choose from*, which does not translate into a better chosen structure. Not worth refolding 16
peptides.

**The durable finding is that ~2.2–2.4 Å clearance to KIX is INVARIANT to the stand-in.** That
strain is inherent to grafting a rigid pre-generated rotamer into a pocket that folded around a
different side chain, and no canonical stand-in fixes it. → **A post-graft relaxation step is
required**, not a better stand-in. We need Rosetta params for FastRelax anyway (Step 5/6), so
verify there that relax pulls these contacts out to a physical distance (~3.4 Å for C–C vdW).

**Implementation:** `run_residuex_swap.py` + `run_control_residuex.sh`. Per peptide: extract
chain B → apply each `(residue_id, smiles)` in sequence → rename to the proper code → recombine
with chain A → write `control_run/ncaa_structures/<name>_ncaa.pdb`. Copy the 3 no-swap peptides
(`WT_MLL`, `MLL1`, `cMyb_native`) through unchanged. **Keep the per-ncAA SDFs — Step 5 needs them.**

### Step 5 — Rosetta params: 6 of 7 ALREADY EXIST (scope collapsed)

**`molfile_to_params.py` is not needed for 6 of the 7, and is not even shipped with PyRosetta.**
Rosetta's `fa_standard` set already contains our residues in `database/chemical/
residue_type_sets/fa_standard/residue_types/l-ncaa/`, and all six are listed in the default
`residue_types.txt`. **Verified live** under `score_peptide.py`'s exact init flags (including
`-corrections::beta_nov16 true`): all load, all `polymer=True`, **no `-extra_res_fa` required**.

| our ncAA | Rosetta params file | code | heavy atoms (canonical names) |
|---|---|---|---|
| Bcs | `BCS.params` | **BCS** | N CA C O CB SG CD CE CZ1 CZ2 CT1 CT2 CI |
| 2mF | `2-methyl-phenylalanine.params` | **A48** | N CA C O CB CG CD1 CD2 CE1 CE2 CE3 CZ |
| CyHex | `beta-cyclohexyl-alanine.params` | **C00** | … CB CG CD1 CD2 CE1 CE2 CZ + *VCG VCD1* |
| CyPent | `beta-cyclopentyl-alanine.params` | **C01** | … CB CG CD1 CD2 CE1 CE2 + *VCG VCD2* |
| Aba | `ABA.params` | **ABA** | N CA C O CB CG |
| nLeu | `NLU.params` | **NLU** | N CA C O CB CG CD CE |
| **C-IAA** | — | **CIA** | ⚠️ **absent — the only custom params needed** |

Rosetta's `BCS` is chemically exactly S-benzylcysteine (`CB–SG–CD` then a 6-carbon aromatic
ring), an independent confirmation of the user's Bcs assignment.

⚠️ `V`-prefixed atoms in C00/C01 (`VCG`, `VCD1`, `VCD2`) are Rosetta **virtual** atoms used to
close the aliphatic ring. Do **not** write them into the PDB — Rosetta adds them itself.

**Remaining work here is only C-IAA** (S-carboxymethyl-cysteine; used 3× — HBS_II, HBS_V,
HBS_VI, all at position 12). The database has just three cysteine derivatives
(`amino-ethyl-cysteine`, `homocysteine`, `tert-butyl-cysteine`) and none matches. Best route is
to hand-derive from `BCS.params`, which already has the correct `CB–SG–CD` linkage plus
`LOWER_CONNECT N` / `UPPER_CONNECT C`: keep the backbone and CB/SG/CD, replace the phenyl ring
with a carboxyl. That reuses a verified polymer backbone rather than trusting
`molfile_to_params` (a ligand-oriented tool) to get polymer connectivity right.

**This also solves the Step-4(d) naming defect**: the canonical atom names above are exactly
what the grafted residue must be renamed to, so names and ResidueTypes agree by construction.

### Step 5b — Atom-name canonicalisation (`fix_graft_naming.py`) — DONE & VERIFIED

Solves Step-4(d). Builds the reference bond graph from the residue's `.params` and solves the
graph isomorphism onto the grafted fragment, assigning canonical names. Anchored on N/C/O
(which ResidueX preserves), so the search collapses immediately; symmetric groups (phenyl ring,
carboxylate oxygens) admit equivalent solutions and any is taken. Also writes the residue as one
contiguous block in sequence position, as `ATOM` in chain B.

**Hydrogens are dropped.** Boltz structures are heavy-atom only, so keeping ResidueX's hydrogens
on just the grafted residue would leave the chain inconsistently protonated; prepwizard and
Rosetta both add them from the ResidueType. V-prefixed virtual atoms are never written.

**Verified end-to-end on MLL2** (`pose_from_pdb`, *without* `-ignore_unrecognized_res`):

```
BEFORE fix: nres=98  peptide = SER ASP ILE ---  ASP PHE VAL LEU LYS ASN THR PRO
AFTER  fix: nres=99  peptide = SER ASP ILE BCS ASP PHE VAL LEU LYS ASN THR PRO
```

99 = 87 KIX + 12 peptide. The 98 is the silent deletion predicted in (4a)/(6) actually happening.
This closes it.

### Step 5c — ⚠️ PyRosetta ships NO ncAA rotamer libraries — DONE

The database params for 5 of our 7 request
`NCAA_ROTLIB_PATH ncaa_rotamer_libraries/alpha_amino_acid/<x>.rotlib`, but PyRosetta ships that
directory **empty** (18 `.rotlib` files total, none ours). FastRelax therefore died with
"Could not open rotamer library file" on **17 of 24** controls — only canonical and BCS-only
peptides survived.

Fix keyed off the database itself: **`BCS` works and has neither a rotlib path nor rotamer
bins.** So `control_run/rosetta_params/` now holds local copies of the five with both
`NCAA_ROTLIB_*` lines stripped, renamed to avoid colliding with the database entries:

| ncAA | db code | our code | source |
|---|---|---|---|
| 2mF | A48 | **M2F** | `2-methyl-phenylalanine.params` |
| CyHex | C00 | **CHX** | `beta-cyclohexyl-alanine.params` |
| CyPent | C01 | **CY5** | `beta-cyclopentyl-alanine.params` (CPN taken by the PNA set) |
| Aba | ABA | **ABU** | `ABA.params` |
| nLeu | NLU | **NLE** | `NLU.params` |

Chemistry, atom names, ICOOR tree and connect records unchanged. `BCS` is used straight from the
database; `CIA` is ours. All 6 non-BCS params are passed via `--extra-res-fa`.

⚠️ **Race to avoid:** relabelling structures while the swap batch is still writing silently
restores old codes. `HBS-22-nLeu` was rewritten after relabelling, reverted to `NLU`, and Rosetta
fell back to the database entry and died. Always confirm the batch has exited before relabelling.

### Step 5d — DSSP: reuse the canonical run, do not recompute

**The graft preserves the backbone EXACTLY**: max deviation of chain-B `N/CA/C/O` between the
pre-graft Boltz structure and the final ncAA structure is **0.0000 Å across all 24**, with no
backbone atoms missing. DSSP secondary structure depends only on backbone H-bond geometry, so
`helix_score` is *provably* identical to the canonical run — reuse `control_dssp_data.csv`.

(Attempting to recompute fails anyway: gemmi/DSSP only recognises `BCS` as a real PDB chemical
component; our invented codes drop out of the polymer and break the chain, so only 7/24 pass.)

This also converts Part-2 finding 5 from an argument into a measurement: `helix_score` and the
Boltz confidences are **exactly** frozen across the swap, so only `hit_num` and the BindCraft
terms can respond.

### Step 5e — FastRelax DOES resolve the graft strain — RESOLVED

The open question from the Met/Phe A/B test. Measured on `HBS_II` (4 grafts), min heavy-atom
distance to KIX before vs after `relax_structure`:

| graft | before | after |
|---|---|---|
| CHX4 | 2.42 Å | **3.62 Å** |
| M2F6 | 3.37 Å | **3.53 Å** |
| CY511 | 2.95 Å | **3.57 Å** |
| CIA12 | 2.49 Å | **4.51 Å** |

All reach proper van der Waals contact (C–C vdW ≈ 3.4 Å) despite
`constrain_relax_to_start_coords(True)`. The ~2.3 Å contacts were a starting-geometry artifact,
confirming that keeping Met over Phe as the Bcs stand-in cost nothing.

### Step 6 — Patch `score_peptide.py` for ncAA support
Four changes, all backward-compatible (defaults preserve current library behaviour):
1. Add `--extra-res-fa` (repeatable / dir) → append `-extra_res_fa <...>` to the options list at `:50-56`.
2. Add `--strict-residues` that **drops `-ignore_unrecognized_res`** so a missing params file crashes loudly.
3. After `pose_from_pdb` (`:167`, `:240`) assert pose residue count == input count; abort on mismatch.
4. Extend `THREE_TO_ONE` (`:24-29`) with the 7 codes mapped to their parent letter. Do **not**
   map to a letter outside `AMINO_ACIDS` (`:23`) or `interface_AA[aa_type] += 1` (`:259`) raises.
Temporarily drop `-mute all` during the first run to read Rosetta's warnings.

### Step 7 — Re-run downstream stages on the swapped structures
- **6** BindCraft — `run_control_bindcraft.sh` pointed at `ncaa_structures/`, with the new params flags.
  (Runner already carries the `PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft` fix.)
- **7** Schrödinger — cif2pdb not needed (already PDB); prepwizard → `schrodinger_calc_hbond.py`
  → `hbond_hit_num.py`. ~50 s/structure, ~20 min for 24, in-session.
- **8** DSSP and **9** face determination — cheap; re-run to confirm backbone is unperturbed.

### Step 8 — Scoring with a swap-sensitive variant
Add to `kix_scoring.py` a third score alongside `priority_score` /
`priority_score_no_enrichment`, weighted toward the terms that actually move post-swap
(`hit_num`, `interface_delta_unsat_hbonds`, `interface_dG`), excluding the frozen
`protein_iptm` / `helix_score`. Report all three in `analyze_controls.py`.

### Step 9 — Re-measure and compare
Rebuild `control_scores.csv` and the per-face `*_with_controls.csv`, then recompute:
- Spearman(score, log10 affinity) and positives-vs-negatives, canonical run vs ncAA run
- the noise floor from the 2 surviving identical pairs
- per-control before/after delta in `hit_num`, `interface_dG`, and each score
Keep the canonical outputs (copy to `control_run/canonical_run/`) so the comparison is direct.

---

## Verification

1. **Swap correctness** — for 2–3 peptides, confirm in PyMOL/ChimeraX that the ncAA sits at the
   right position with sane geometry, and that chain A (KIX) is untouched (RMSD ≈ 0 vs input).
2. **No silent deletion** — the Step-6 residue-count assertion must pass for all 24; verify the
   relaxed PDB still contains the ncAA resname.
3. **Degeneracy check** — confirm 22 unique structures, and that the 2 known-identical pairs
   still produce near-identical metrics (they are the noise control).
4. **Regression** — re-run `analyze_and_score_all_metrics.py` and confirm the library numbers
   are still 30776 → 28029 → 2551 → 198 (111 / 87). The control work must not perturb the library.
5. **Payoff** — does Spearman vs affinity improve materially over `rho=−0.20`, and does the
   positive-vs-negative gap rise above the measured noise floor?

## Risks / open items

- **`molfile_to_params.py` on 7 ncAAs is the biggest execution risk** — beta_nov16 typing,
  backbone connectivity for a polymer residue (needs `--polymer` handling), and correct
  N/C connect atoms. If a residue resists params generation, fall back to reporting
  Schrödinger + Boltz metrics for that peptide and flagging BindCraft as N/A rather than wrong.
- **X/Z HBS cap is unmodelled**, leaving 2 pairs degenerate including one positive/negative pair.
  Decide whether to model the caps or accept them as noise controls.
- **Boltz confidences are frozen pre-swap** by construction — the swap-sensitive score is the
  mitigation, not a fix.
- ResidueX docs list no supported-ncAA registry, so SMILES quality is entirely on the input table.
