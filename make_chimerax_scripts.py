#!/usr/bin/env python3
"""Generate the ChimeraX figure scripts.

Interactions are DRAWN FROM SCHRODINGER, not recomputed by ChimeraX.

ChimeraX's `hbonds` disagrees badly here: Schrodinger measures hydrogen-to-
acceptor <= 2.8 A on prepwizard output that has explicit hydrogens, while the
aligned PDBs have none, so ChimeraX falls back to heavy-atom geometry. On
Full_Library_Hit_1530 that is 4 H-bonds vs 1 -- and none of ChimeraX's land on a
hit residue. `hbonds` also cannot detect pi-pi at all, which is most of the MLL
interface. So every polar/aromatic line is an explicit `pbond` taken from
atom_level_interactions.csv, and the figures match the scored CSVs exactly.

Apolar packing still uses ChimeraX `contacts`: it is pure heavy-atom proximity
with no chemistry criterion to disagree about, and drawing every pair as its own
pseudobond would be unreadable.

    python3 make_chimerax_scripts.py [--base <windows path>]
"""
import argparse, os, shutil
import pandas as pd

P = os.path.dirname(os.path.abspath(__file__))
M = f"{P}/full_library_all_metrics"
OUT = f"{P}/chimerax_scripts"
DEFAULT_BASE = r"C:\Users\jonat\Documents\NYU_Research\new\chimerax_scripts"
NCH = {"cmyb": "A", "mll": "C"}          # 2AGH peptide chains
OFFSET = 585


def digits(s):
    return int("".join(c for c in s if c.isdigit()))


def scene(title, base, pdbfile, kc, pc, hits, apolar, show_pep, inter, png,
          renumber_kix=False):
    polar = [r for r in hits if r not in apolar]
    hs = ",".join(map(str, hits))
    ap = ",".join(map(str, apolar))
    lines = [f"# {title}", "close session",
             "# bgColor: change to 'black' or 'dark gray' for dark slides.",
             "# The save line at the bottom is transparent either way.",
             "set bgColor white", "lighting soft",
             "graphics silhouettes true width 1.5",
             f"open {base}\\structures\\{pdbfile}", "",
             "# Hydrogens are removed for a clean figure. Interactions below are",
             "# drawn from Schrodinger, so nothing depends on them being present.",
             "delete H", ""]
    if renumber_kix:
        lines += ["# KIX crystal numbering is +585 vs the pipeline. Renumber to 1-87 so",
                  "# these figures compare to the binder panels residue-for-residue.",
                  "# Chain B is contiguous 586-672 (87 residues), so `start 1` is exact.",
                  f"renumber /{kc} start 1", ""]
    lines += ["hide atoms", "show cartoons",
              f"color /{kc} gray(180)", f"color /{pc} cornflower blue", "",
              "# key hit residues: gold = can H-bond, tan = apolar",
              f"show /{kc}:{hs} atoms", f"style /{kc}:{hs} stick",
              f"color /{kc}:{','.join(map(str,polar))} goldenrod",
              f"color /{kc}:{ap} tan",
              f'label /{kc}:{hs} residues text "{{0.name}}{{0.number}}" height 1.1 color black', "",
              "# contacting peptide residues",
              f"show /{pc}:{show_pep} atoms", f"style /{pc}:{show_pep} stick",
              f'label /{pc}:{show_pep} residues text "{{0.name}}{{0.number}}" height 1.0 color blue', ""]

    # HIT RESIDUES ONLY. A pbond to a non-hit KIX residue would never render
    # (that residue is not `show`n) and the comparison is defined on the hit set,
    # so emitting it is dead code that only confuses anyone reading the script.
    hset = set(hits)
    inter = [r for r in inter if r.kix_resnum in hset]
    hb = [r for r in inter if r.interaction == "hbond"]
    pp = [r for r in inter if r.interaction == "pi_pi"]
    lines += ["# ---- interactions AS SCHRODINGER FOUND THEM (not recomputed) ----",
              "# red = H-bond.  These come from atom_level_interactions.csv, so the",
              "# figure and the scored CSVs cannot disagree."]
    if hb:
        for r in hb:
            lines.append(f"pbond /{kc}:{r.kix_resnum}@{r.kix_atom} /{pc}:{r.pep_resnum}@{r.pep_atom}"
                         f" color red radius 0.12 dashes 6"
                         f"   # {r.kix_resname}{r.kix_resnum} .. {r.pep_resname}{r.pep_resnum}"
                         .split("   #")[0])
    else:
        lines.append("# (none on hit residues for this structure)")
    lines += ["", "# purple = pi-pi stacking. ChimeraX has NO command for this;",
              "# without these lines the aromatic packing is invisible."]
    if pp:
        for r in pp:
            lines.append(f"pbond /{kc}:{r.kix_resnum}@{r.kix_atom} /{pc}:{r.pep_resnum}@{r.pep_atom}"
                         f" color medium purple radius 0.12 dashes 3")
    else:
        lines.append("# (none)")
    lines += ["", "# orange = van der Waals packing, on ALL hit residues.",
              "# NOT just the apolar ones: TYR46 (31 atom pairs) and TYR65 (44, the",
              "# largest contact on the c-Myb face) pack hydrophobically while making",
              "# no H-bond or stack, so restricting this to apolar residues made them",
              "# vanish from the figure entirely.",
              "# Geometric proximity only, so ChimeraX and Schrodinger agree here.",
              f"contacts /{kc}:{hs} restrict /{pc} overlapCutoff -0.8 hbondAllowance 0.0"
              f" color orange radius 0.09 dashes 2", "",
              f"view /{kc}:{hs}",
              f"# save {base}\\png\\{png} width 2400 supersample 3 transparentBackground true"]
    return "\n".join(lines) + "\n"


def main():
    ap_ = argparse.ArgumentParser()
    ap_.add_argument("--base", default=DEFAULT_BASE)
    a = ap_.parse_args()
    base = a.base.rstrip("\\")

    os.makedirs(f"{OUT}/structures", exist_ok=True)
    os.makedirs(f"{OUT}/png", exist_ok=True)
    k = pd.read_csv(f"{M}/key_residues_top10.csv")
    nat = pd.read_csv(f"{M}/key_residues_native.csv")
    ai = pd.read_csv(f"{M}/atom_level_interactions.csv")

    lines = []
    for (st, rank, pep, name), g in k.groupby(["set", "rank", "peptide", "name"]):
        face = g.face.iloc[0]
        if face not in NCH:
            continue
        pc = g.chain.iloc[0]
        src = (f"{P}/dual_cofold/aligned/{pep}.pdb" if st == "top10_dual"
               else f"{P}/aligned_structures_v2pde/{name}.pdb")
        pdbfile = f"{pep}.pdb" if st == "top10_dual" else f"{name}.pdb"
        shutil.copy(src, f"{OUT}/structures/{pdbfile}")

        hits = sorted(g.kix_resnum_model.unique())
        apolar = sorted(g[g.is_apolar_hit].kix_resnum_model.unique())
        inter = list(ai[ai.structure_name == name].itertuples())
        # every pbond endpoint must be displayed, so union them into the show list
        peps = {digits(p) for s in g.binder_partners.dropna() for p in s.split("+")}
        peps |= {r.pep_resnum for r in inter}
        stem = f"{st}_r{int(rank):02d}_{name}"
        open(f"{OUT}/{stem}.cxc", "w", newline="").write(
            scene(f"{name} -- {face} face, rank {int(rank)} ({st})", base, pdbfile,
                  "A", pc, hits, apolar, ",".join(map(str, sorted(peps))),
                  inter, f"{stem}.png"))
        lines.append((st, int(rank), f"{stem}.cxc"))

    shutil.copy(f"{P}/2agh_model1.pdb", f"{OUT}/structures/2agh_model1.pdb")
    for face, pc in NCH.items():
        g = nat[nat.face == face]
        hits = sorted(g.kix_resnum_model.unique())
        apolar = sorted(g[g.is_apolar_hit].kix_resnum_model.unique())
        inter = list(ai[ai.structure_name == f"2AGH_native_{face}"].itertuples())
        peps = {digits(p) for s in g.binder_partners.dropna() for p in s.split("+")}
        peps |= {r.pep_resnum for r in inter}
        open(f"{OUT}/crystal_{face}.cxc", "w", newline="").write(
            scene(f"2AGH native -- {face} face", base, "2agh_model1.pdb", "B", pc,
                  hits, apolar, ",".join(map(str, sorted(peps))), inter,
                  f"crystal_{face}.png", renumber_kix=True))
        lines.append(("crystal", 0, f"crystal_{face}.cxc"))

    order = {"crystal": 0, "top10_cmyb": 1, "top10_mll": 2, "top10_dual": 3}
    lines.sort(key=lambda x: (order[x[0]], x[1]))
    with open(f"{OUT}/OPEN_THESE.txt", "w", newline="\r\n") as f:
        f.write("# Paste one line at a time into the ChimeraX command line.\n"
                "# Each script starts with `close session`, so the next line replaces the last.\n")
        cur = None
        for st, rank, fn in lines:
            if st != cur:
                f.write(f"\n# ---- {st} ----\n"); cur = st
            f.write(f"open {base}\\{fn}\n")
    print(f"wrote {len(lines)} scripts to {OUT} (base: {base})")


if __name__ == "__main__":
    main()
