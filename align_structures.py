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

for name in mobile_struct_dict:
    structure.StructureWriter.write(mobile_struct_dict[name], str(OUTPUT_DIR / f"{name}.cif"))


