#!/bin/bash
#SBATCH --job-name=pdb_to_cif_full
#SBATCH --output=pdb_to_cif_full_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06

cd /scratch/jem9759/ZhangWork/KIX_Project
source activate /scratch/jem9759/envs/general_penv

python3 ./pdb_to_cif.py /scratch/jem9759/pdb_sts_full_library /scratch/jem9759/pdb_sts_full_library_cif
