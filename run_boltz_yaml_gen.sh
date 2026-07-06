#!/bin/bash
#SBATCH --job-name=boltz_yaml_gen_full
#SBATCH --output=boltz_yaml_gen_full_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --account=torch_pr_149_chemistry

module purge
module load anaconda3/2025.06

source activate /scratch/jem9759/envs/general_penv

python3 ./boltz_yaml_gen.py /scratch/jem9759/ZhangWork/KIX_Project/tsv_outputs/full_library_out.tsv --out_path /scratch/jem9759/ZhangWork/KIX_Project/yaml_outputs/out_yaml_full_library --helix_dir /scratch/jem9759/ZhangWork/KIX_Project/cif_outputs/pdb_sts_full_library_cif --kix-msa /scratch/jem9759/ZhangWork/KIX_Project/kix_msa.csv
