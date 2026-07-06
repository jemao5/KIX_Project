#!/bin/bash
#SBATCH --array=0-62
#SBATCH --cpus-per-task=1        # (or a few; relax is single-threaded mostly)
#SBATCH --mem=8G
#SBATCH --time=04:00:00          # ~2.3hr/chunk + margin
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --qos=cpu48
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/bindcraft_logs/bindcraft_%A_%a.out
# NO GPU (CPU-only)

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

CHUNK=$(printf "chunk_%04d" $SLURM_ARRAY_TASK_ID)
cd /scratch/jem9759/ZhangWork/BindCraft   # repo dir for functions/ + settings_filters/

# activate BindCraft env by full scratch path (avoid the activation bug)
# source $(conda info --base)/etc/profile.d/conda.sh
#
module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/BindCraft

python -u score_peptide.py \
    --input "/scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full/${CHUNK}/boltz_results_${CHUNK}/predictions/*/*_model_0.cif" \
    --binder-chain B --target-chain A \
    --filters settings_filters/peptide_filters.json \
    --dalphaball-path functions/DAlphaBall.gcc \
    --output-json /scratch/jem9759/ZhangWork/KIX_Project/bindcraft_out/${CHUNK} \
    --output-relaxed-pdb /scratch/jem9759/ZhangWork/KIX_Project/bindcraft_relaxed/${CHUNK}
