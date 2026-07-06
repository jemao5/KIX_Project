#!/bin/bash

#SBATCH --job-name=full_library_face_determination
#SBATCH --output=full_library_face_determination_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06

source activate /scratch/jem9759/envs/dssp_env

python3 /scratch/jem9759/ZhangWork/KIX_Project/full_library_face_determination.py "/scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full/chunk_*/*/predictions/*/*.cif" /scratch/jem9759/ZhangWork/KIX_Project/final_metric_outputs/full_library_face_assignment.tsv  --workers 16
