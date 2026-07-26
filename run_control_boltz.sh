#!/bin/bash
#SBATCH --job-name=boltz_controls
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/control_run/logs/boltz_controls_%j.out
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry

# Stage 4 for the literature control set (control_list.csv, 24 peptides).
# Same settings as the full library (run_chunked_boltz.sh): precomputed KIX MSA
# in the YAMLs, idealized-helix peptide templates, no pocket constraints.

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/boltz_env

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

P=/scratch/jem9759/ZhangWork/KIX_Project

boltz predict "$P/control_run/yamls" \
    --cache /scratch/jem9759/.boltz \
    --out_dir "$P/control_run/boltz_out"
