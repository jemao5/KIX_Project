#!/bin/bash
#SBATCH --job-name=dual_cofold
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold/logs/dual_cofold_%j.out
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --account=torch_pr_149_chemistry

# KIX + 2 copies of each peptide (chains A/B/C) for the dual-face test.
# 82 structures: top 40 per face by priority_score_v2_pde, plus cMyb_native and
# WT_MLL as single-face sanity controls (neither should classify DUAL).

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/boltz_env
export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

D=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold
boltz predict "$D/yamls" --cache /scratch/jem9759/.boltz --out_dir "$D/boltz_out"
