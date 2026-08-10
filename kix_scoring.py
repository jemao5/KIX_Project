#!/usr/bin/env python3
"""
Shared scoring/filter definitions for the KIX screen.

Both the full-library ranking (`analyze_and_score_all_metrics.py`) and the
literature-control comparison (`analyze_controls.py`) import from here, so the
thresholds a control is audited against are literally the same objects the
library was filtered with -- they cannot drift apart.
"""
import pandas as pd

# --- filter clauses, in pipeline order -------------------------------------
# name -> predicate(DataFrame) -> boolean Series
CONFIDENCE_FILTERS = {
    "confidence_score > 0.85": lambda d: d["confidence_score"] > 0.85,
    "pep_ptm > 0.80":          lambda d: d["pep_ptm"] > 0.80,
    "complex_pde < 0.5":       lambda d: d["complex_pde"] < 0.5,
}

PHYSICAL_FILTERS = {
    "binder_score < 0":                    lambda d: d["binder_score"] < 0.0,
    "interface_dG < -25":                  lambda d: d["interface_dG"] < -25.0,
    "interface_dSASA > 1":                 lambda d: d["interface_dSASA"] > 1.0,
    "surface_hydrophobicity < 1":          lambda d: d["surface_hydrophobicity"] < 1.0,
    "interface_sc > 0.5":                  lambda d: d["interface_sc"] > 0.5,
    "interface_nres > 4":                  lambda d: d["interface_nres"] > 4.0,
    "interface_interface_hbonds > 1":      lambda d: d["interface_interface_hbonds"] > 1,
    "interface_delta_unsat_hbonds <= 5":   lambda d: d["interface_delta_unsat_hbonds"] <= 5,
    "n_K <= 3":                            lambda d: d["n_K"] <= 3,
    "n_M <= 3":                            lambda d: d["n_M"] <= 3,
    "hit_num >= 1":                        lambda d: d["hit_num"] >= 1,
    "count > 1":                           lambda d: d["count"] > 1,
    "protein_iptm > 0.85":                 lambda d: d["protein_iptm"] > 0.85,
}

HELICITY_FILTERS = {
    "helix_score > 0.70": lambda d: d["helix_score"] > 0.70,
}

ALL_FILTERS = {**CONFIDENCE_FILTERS, **PHYSICAL_FILTERS, **HELICITY_FILTERS}

# Clauses that describe the *library design*, not binding quality. They are
# meaningless or unfair for literature controls: `count` is a library
# enrichment tally the controls were never part of, and n_K/n_M were synthesis
# constraints on the designed library.
LIBRARY_DESIGN_FILTERS = ["count > 1", "n_K <= 3", "n_M <= 3"]


def apply_filters(df, filters):
    """Return the subset of df passing every clause in `filters`."""
    mask = pd.Series(True, index=df.index)
    for pred in filters.values():
        mask &= pred(df)
    return df[mask]


def filter_audit(df, filters=None):
    """One boolean column per filter clause, plus n_failed / failed_filters."""
    filters = ALL_FILTERS if filters is None else filters
    audit = pd.DataFrame(index=df.index)
    missing = []
    for label, pred in filters.items():
        try:
            audit[label] = pred(df).fillna(False).astype(bool)
        except KeyError:
            # Previously this silently recorded the clause as FAILED for every
            # row, which makes a dataframe that merely lacks the column look
            # like nothing passes. Collect and raise instead.
            missing.append(label)
    if missing:
        raise KeyError(
            f"filter_audit: {len(missing)} clause(s) reference columns absent from the "
            f"dataframe: {missing}. Pass the full merged metrics, not a preview subset.")
    audit["n_failed"] = (~audit[list(filters)]).sum(axis=1)
    audit["failed_filters"] = audit[list(filters)].apply(
        lambda r: "; ".join(c for c in filters if not r[c]), axis=1
    )
    return audit


# --- priority score --------------------------------------------------------
def add_priority_score(df, count_col="count", hit_col="hit_num",
                       unsat_col="interface_delta_unsat_hbonds",
                       protein_col="protein_iptm", helix_score_col="helix_score",
                       count_weight=0.2, hit_weight=0.2, unsat_weight=0.1,
                       protein_weight=0.2, helix_weight=0.2,
                       score_col="priority_score"):
    """Percentile-rank blend of the ranking metrics.

    A weight of 0.0 drops its term entirely rather than multiplying it by zero,
    so a zero-weighted column that is NaN (e.g. `count` for the literature
    controls, which have no library enrichment) does not poison the score.
    """
    d = df.copy()
    terms = [
        (count_weight,   d[count_col].rank(pct=True, method="average")),
        (hit_weight,     d[hit_col].rank(pct=True, method="average")),
        (unsat_weight,   1.0 - d[unsat_col].rank(pct=True, method="average")),
        (protein_weight, d[protein_col].rank(pct=True, method="average")),
        (helix_weight,   d[helix_score_col].rank(pct=True, method="average")),
    ]
    terms = [(w, r) for w, r in terms if w != 0.0]
    weights = sum(w for w, _ in terms)
    d[score_col] = sum(w * r for w, r in terms) / weights
    return d.sort_values(score_col, ascending=False)


# --- score variants -------------------------------------------------------
# Each entry is a list of (column, weight, invert). `invert=True` means LOWER
# raw values are better, so the percentile rank is flipped (1 - rank).
#
# All three share the same hit_num / unsat / protein_iptm core and differ ONLY
# in the fourth slot:
#   helix : the original -- DSSP chain_b_helix_fraction
#   nohelix: that slot removed entirely, remaining weights renormalised
#   pde   : that slot replaced by complex_pde (Boltz predicted distance error)
#
# Motivation (see ControlPlan.md): on the 23 MLL literature controls the
# helix-fraction term does not separate positives from negatives (AUC 0.553,
# p=0.345); removing it gives AUC 0.731 (p=0.032) and replacing it with
# complex_pde gives AUC 0.792 (p=0.010) -- the only variant whose
# positive/negative gap exceeds the measured noise floor.
#
# ⚠️ helix_fraction is length-dependent (rho=-0.824 with peptide length): the
# helical core is ~8 residues regardless of length, so the fraction falls as the
# denominator grows. Use chain_b_helix_count for QC instead of the fraction.
_CORE = [("hit_num", 0.2, False),
         ("interface_delta_unsat_hbonds", 0.1, True),
         ("protein_iptm", 0.2, False)]

SCORE_VARIANTS = {
    "":          _CORE + [("helix_score", 0.2, False)],   # original
    "_no_helix": _CORE,
    "_pde":      _CORE + [("complex_pde", 0.2, True)],
}

# Same three blends but with `hit_num` swapped for `hit_num_v2`, which adds a
# contact term for the apolar hit residues that cannot H-bond (see
# hydrophobic_contacts.py). On the 23 MLL controls hit_num_v2 separates
# positives from negatives at MWU p=0.007 / rho=-0.622 vs p=0.037 / rho=-0.519
# for hit_num. Variants are only built when the column is present, so nothing
# breaks on tables predating it.
_CORE_V2 = [("hit_num_v2", 0.2, False)] + _CORE[1:]
SCORE_VARIANTS_V2 = {
    "_v2":          _CORE_V2 + [("helix_score", 0.2, False)],
    "_v2_no_helix": _CORE_V2,
    "_v2_pde":      _CORE_V2 + [("complex_pde", 0.2, True)],
}


def _blend(d, spec):
    terms = [(w, (1.0 - d[c].rank(pct=True, method="average")) if inv
              else d[c].rank(pct=True, method="average"))
             for c, w, inv in spec if c in d.columns]
    total = sum(w for w, _ in terms)
    return sum(w * r for w, r in terms) / total


def add_score_variants(df, with_count=True):
    """Add every score variant. Existing columns are NOT modified or removed.

    Produces, for each variant V in SCORE_VARIANTS:
      priority_score{V}                  -- includes the enrichment `count` term
      priority_score{V}_no_enrichment    -- `count` dropped, weights renormalised

    `priority_score` and `priority_score_no_enrichment` (V="") are the original
    definitions and must remain bit-identical to the previous implementation.
    """
    d = df.copy()
    variants = dict(SCORE_VARIANTS)
    if "hit_num_v2" in d.columns:
        variants.update(SCORE_VARIANTS_V2)
    for suffix, spec in variants.items():
        d[f"priority_score{suffix}_no_enrichment"] = _blend(d, spec)
        if with_count and "count" in d.columns:
            d[f"priority_score{suffix}"] = _blend(d, [("count", 0.2, False)] + spec)
    return d


def add_swap_sensitive_score(df, score_col="priority_score_swap_sensitive"):
    """Percentile blend restricted to terms that RESPOND to an ncAA graft.

    `protein_iptm` and `pep_ptm` are computed by Boltz during folding, i.e.
    BEFORE ResidueX runs, and `helix_score` reads the backbone, which the graft
    leaves untouched. Those terms are therefore frozen across the swap and are
    excluded here -- in `priority_score_no_enrichment` they carry 0.4 of the 0.7
    total weight, so ~57% of that score cannot react to the substitution at all.

    Terms kept, and why:
      hit_num                        0.4  H-bond / pi-stacking count to the face's
                                          hit residues -- specificity-weighted
      interface_dG                   0.4  Rosetta interface energy
      interface_delta_unsat_hbonds   0.2  buried unsatisfied H-bonds (lower better)

    ⚠️ Interpret `interface_dG` with care on grafted structures: it rewards
    buried hydrophobic surface, so adding a benzyl or cyclohexyl improves it
    whether or not binding improves. Measured on this control set, grafting
    shifted it by -6.0 kcal/mol on average and the largest gains went to
    *negative* controls.
    """
    d = df.copy()
    hit = d["hit_num"].rank(pct=True, method="average")
    dG = (-d["interface_dG"]).rank(pct=True, method="average")   # more negative = better
    unsat = 1.0 - d["interface_delta_unsat_hbonds"].rank(pct=True, method="average")
    d[score_col] = (0.4 * hit + 0.4 * dG + 0.2 * unsat)
    return d.sort_values(score_col, ascending=False)


def add_both_priority_scores(df):
    """priority_score (with enrichment `count`) + priority_score_no_enrichment.

    The no-enrichment variant re-normalizes the remaining four weights
    (hit_num .2, unsat .1, protein_iptm .2, helix .2 -> /0.7). It exists so
    peptides with no library enrichment count -- the literature controls --
    can be ranked on the same scale as library members.
    """
    d = add_priority_score(df)
    d = add_priority_score(d, count_weight=0.0,
                           score_col="priority_score_no_enrichment")
    return d.sort_values("priority_score", ascending=False)


# --- hit_num residue sets -------------------------------------------------
# `hit_num` counts peptide<->KIX H-bonds and pi-pi stacking, but ONLY those
# landing on a per-face list of "hit residues". Three lists are available.
#
#   original : the hand-picked list, derived from ChimeraX VDW CONTACTS on 2AGH.
#              Correct as a description of the binding face, but mismatched to
#              the metric: contacts include hydrophobic packing, H-bonds do not,
#              so 6 of 8 (c-Myb) and 5 of 7 (MLL) of these residues make ZERO
#              H-bond/pi-pi in the crystal and can never contribute to hit_num.
#   crystal  : residues that actually H-bond or pi-stack with the NATIVE peptide
#              in 2AGH model 1, derived by derive_hit_residues.py using the same
#              Schrodinger workflow that computes hit_num.
#   union    : original | crystal.
#
# Cross-validation: of the 99 interactions our 23 MLL controls actually make,
# `original` captures 15%, `crystal` 52%, `union` 57%. The three residues our
# controls H-bond to most (ARG 86 x23, ARG 83 x18, ARG 39 x9) are all in
# `crystal`; two are absent from `original`.
#
# GLU 78 and LYS 82 are deliberately EXCLUDED despite being frequent partners in
# our control structures (13x and 11x): neither contacts the native peptide in
# 2AGH (5.63 A and 6.99 A minimum heavy-atom distance), and including them would
# mean fitting the residue list to the predictions we are trying to validate.
#
# Residue numbers are in pipeline numbering (1-87). 2AGH uses CBP numbering;
# offset is 585.
KIX_SEQUENCE = ("GVRKGWHEHVTQDLRSHLVHKLVQAIFPTPDPAALKDRRMENLVAYAKKVEGDMYESANSRD"
                "EYYHLLAEKIYKIQKELEEKRRSRL")

_ONE_TO_THREE = {"A":"ALA","R":"ARG","N":"ASN","D":"ASP","C":"CYS","Q":"GLN","E":"GLU",
                 "G":"GLY","H":"HIS","I":"ILE","L":"LEU","K":"LYS","M":"MET","F":"PHE",
                 "P":"PRO","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR"}

_HIT_RESIDUE_NUMBERS = {
    "original": {"cmyb": [14, 18, 21, 65, 69, 72, 73, 76],
                 "mll":  [27, 39, 43, 46, 71, 75, 79]},
    "crystal":  {"cmyb": [9, 21, 61, 76, 80, 81],
                 "mll":  [27, 39, 83, 84, 86]},
}
_HIT_RESIDUE_NUMBERS["union"] = {
    face: sorted(set(_HIT_RESIDUE_NUMBERS["original"][face])
                 | set(_HIT_RESIDUE_NUMBERS["crystal"][face]))
    for face in ("cmyb", "mll")
}

HIT_RESIDUE_SETS = sorted(_HIT_RESIDUE_NUMBERS)


def hit_residue_strings(set_name, face, chain="A"):
    """['LEU 14 A', ...] for one face. Residue NAMES are derived from
    KIX_SEQUENCE rather than transcribed, so the numbers are the only thing
    that can be got wrong."""
    if set_name not in _HIT_RESIDUE_NUMBERS:
        raise KeyError(f"unknown hit-residue set {set_name!r}; "
                       f"choose from {HIT_RESIDUE_SETS}")
    out = []
    for i in _HIT_RESIDUE_NUMBERS[set_name][face]:
        if not 1 <= i <= len(KIX_SEQUENCE):
            raise ValueError(f"residue {i} outside KIX 1-{len(KIX_SEQUENCE)}")
        out.append(f"{_ONE_TO_THREE[KIX_SEQUENCE[i - 1]]} {i} {chain}")
    return out


# Side chains that carry no H-bond donor or acceptor. A hit residue of one of
# these types can NEVER contribute to hit_num, which counts only H-bonds and
# pi-pi -- on the MLL face that is 4 of the 7 listed residues (PHE 27, LEU 43,
# ILE 75, LEU 79), i.e. the hydrophobic packing that actually drives MLL binding
# is invisible to the metric. `hit_num_v2` adds a contact term for exactly these
# residues; see hydrophobic_contacts.py.
APOLAR_AA = set("AVLIMFWPCG")


def apolar_hit_residues(set_name, face):
    """Hit residues whose side chain cannot H-bond, so they can only ever be
    detected as a contact. Membership is derived from KIX_SEQUENCE rather than
    hardcoded, so it stays correct if the residue list changes."""
    return [i for i in _HIT_RESIDUE_NUMBERS[set_name][face]
            if KIX_SEQUENCE[i - 1] in APOLAR_AA]
