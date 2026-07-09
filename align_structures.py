import argparse
from pathlib import Path
from schrodinger import structure
from schrodinger.structutils import structalign2
import pandas as pd

BOLTZ_ROOT = Path("boltz_out_full")
OUTPUT_DIR = Path("/scratch/jem9759/ZhangWork/KIX_Project/aligned_structures")
n_candidates = 20
REFERENCE_PATH = "/scratch/jem9759/ZhangWork/KIX_Project/2agh_model1.cif"
MLL_STRUCTURE_LIST = "/scratch/jem9759/ZhangWork/KIX_Project/full_library_all_metrics/mll_candidates.csv"
CMYB_STRUCTURE_LIST = "/scratch/jem9759/ZhangWork/KIX_Project/full_library_all_metrics/cmyb_candidates.csv"


def build_cif_index(root):
    """Map each peptide name -> its Boltz model_0 .cif path, across all chunks."""
    index = {}
    for cif in root.rglob("*_model_0.cif"):
        name = cif.name[:-len("_model_0.cif")]   # 'Full_Library_Hit_42'
        index[name] = cif
    return index

def standardize_his_names(st):
    """Rename Schrodinger's HID/HIE/HIP protonation-state names back to standard HIS
    so ChimeraX recognizes the residues and draws continuous ribbon."""
    for atom in st.atom:
        resname = atom.pdbres.strip()
        if resname in ("HID", "HIE", "HIP"):
            atom.pdbres = "HIS "   # note: pdbres is typically space-padded to 4 chars
    return st

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

mll_df = pd.read_csv(MLL_STRUCTURE_LIST)
cmyb_df = pd.read_csv(CMYB_STRUCTURE_LIST)

mll_top_candidates = mll_df.head(n_candidates)
cmyb_top_candidates = cmyb_df.head(n_candidates)

to_be_aligned_df = pd.concat([mll_top_candidates, cmyb_top_candidates], ignore_index=True)


path_dict = build_cif_index(BOLTZ_ROOT)

ref_struct = structure.StructureReader.read(REFERENCE_PATH)

mobile_struct_dict = {}
for name in to_be_aligned_df["name"].tolist():
    if name not in path_dict:
        print(f"WARNING: no Boltz cif for {name}, skipping")
        continue
    mobile_struct_dict[name] = structure.StructureReader.read(str(path_dict[name]))

alignments = structalign2.align_many(ref_struct, list(mobile_struct_dict.values()))

names_w_alignments = zip(alignments, list(mobile_struct_dict.keys()))
alignment_data = []
for alignment, name in names_w_alignments:
    alignment_data.append({"name": name, "rmsd":alignment.rmsd, "score":alignment.score})

alignment_df = pd.DataFrame(alignment_data)
alignment_df.to_csv(str(OUTPUT_DIR / "alignment_data.csv"), index=False)



for name in mobile_struct_dict:
    with structure.StructureWriter(str(OUTPUT_DIR / f"{name}.pdb")) as writer:
        writer.append(standardize_his_names(mobile_struct_dict[name]))



