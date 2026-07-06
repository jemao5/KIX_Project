#!/bin/bash
#SBATCH --job-name=boltz_full_library
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/logs/boltz_full_%A_%a.out
#SBATCH --array=0-62%24
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/boltz_env

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

CHUNK=$(printf "chunk_%04d" $SLURM_ARRAY_TASK_ID)
CHUNK_DIR=/scratch/jem9759/ZhangWork/KIX_Project/yaml_chunks/$CHUNK
OUT_DIR=/scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full/$CHUNK

mkdir -p "$OUT_DIR"
echo "Task $SLURM_ARRAY_TASK_ID: $CHUNK_DIR -> $OUT_DIR"

boltz predict "$CHUNK_DIR" \
    --cache /scratch/jem9759/.boltz \
    --out_dir "$OUT_DIR"