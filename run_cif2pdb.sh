#!/bin/bash
#SBATCH --job-name=cif2pdb
#SBATCH --output=cif2pdb_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06

source activate /scratch/jem9759/envs/dssp_env

python3 cif2pdb.py "/scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full/chunk_*/boltz_results_chunk_*/predictions/Full_Library_Hit_*/Full_Library_Hit_*_model_*.cif" --workers 16
