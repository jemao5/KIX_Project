#!/bin/bash
#SBATCH --job-name=dssp_full_library
#SBATCH --output=dssp_full_library_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06
cd /scratch/jem9759/ZhangWork/KIX_Project
source activate /scratch/jem9759/envs/dssp_env

echo "=== which python3: ==="
which python3
echo "=== sys.executable: ==="
python3 -c "import sys; print(sys.executable)"
echo "=== CONDA_PREFIX: $CONDA_PREFIX ==="


python3 ./check_helix_dssp_skip_first2.py --input-glob "/scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full/chunk_*/*/predictions/*/*.cif" --output "/scratch/jem9759/ZhangWork/KIX_Project/final_metric_outputs/dssp_full_library.csv" --workers 16
