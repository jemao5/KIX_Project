#!/bin/bash
#SBATCH --job-name=dual_decoy
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold/logs/dual_decoy_%j.out
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --account=torch_pr_149_chemistry
# Negative controls for the dual-face assay: poly-Ala, a scrambled hit, and two
# MEASURED non-binders (MLL6_YA / MLL6_CstarA, both >9999 uM). If these come out
# DUAL, the assay is docking anything into the two empty grooves and the DUAL
# call carries no information about the sequence.
module purge; module load anaconda3/2025.06
source activate /scratch/jem9759/envs/boltz_env
export TMPDIR=/scratch/jem9759/tmp
D=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold
boltz predict "$D/decoy/yamls" --cache /scratch/jem9759/.boltz --out_dir "$D/decoy/boltz_out"
