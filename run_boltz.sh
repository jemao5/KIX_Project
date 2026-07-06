#!/bin/bash
#SBATCH --job-name=boltz_kix
#SBATCH --output=boltz_kix_%j.out
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --account=torch_pr_149_chemistry

# --- environment (match what worked in your interactive session) ---
module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/boltz_env

# Redirect temp off the tiny /tmp (in case anything needs scratch temp space)
export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

# --- run boltz over ALL yamls in the directory, sequentially ---
# Outputs go to scratch (home is only 50GB; 57 predictions will be large).
boltz predict ./out_yaml_constraints \
    --cache /scratch/jem9759/.boltz \
    --out_dir /scratch/jem9759/boltz_out_constraints_forced
