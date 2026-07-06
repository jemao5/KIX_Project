#!/bin/bash
#SBATCH --job-name=boltz_analysis_full
#SBATCH --output=boltz_analysis_full_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06

source activate /scratch/jem9759/envs/general_penv


python parse_boltz_results.py /scratch/jem9759/ZhangWork/KIX_Project/boltz_out_full --out /scratch/jem9759/ZhangWork/KIX_Project/metrics_full_library.csv --n_a 87
