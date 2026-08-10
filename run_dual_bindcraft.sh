#!/bin/bash
#SBATCH --job-name=dual_bindcraft
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold/logs/dual_bindcraft_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --qos=cpu48

# BindCraft on the 172 per-copy sub-complexes from the dual-face cofold.
# Each is an ordinary 2-chain KIX+peptide complex (chain C was renamed B by
# split_dual_copies.py), so no --extra-res-fa and no chain-C override is needed.
# --output-relaxed-pdb is REQUIRED: without it relaxed PDBs land in the cwd,
# which is the BindCraft repo.

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"
cd /scratch/jem9759/ZhangWork/BindCraft
module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/BindCraft
P=/scratch/jem9759/ZhangWork/KIX_Project
export PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft:$PYTHONPATH

python -u "$P/score_peptide.py" \
    --input "$P/dual_cofold/copies/*.pdb" \
    --binder-chain B --target-chain A \
    --filters settings_filters/peptide_filters.json \
    --dalphaball-path functions/DAlphaBall.gcc \
    --output-json "$P/dual_cofold/bindcraft_out" \
    --output-relaxed-pdb "$P/dual_cofold/bindcraft_relaxed"
