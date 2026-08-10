#!/bin/bash
#SBATCH --job-name=hydrophobic_contacts
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/logs/hydrophobic_contacts_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry

# Counts APOLAR hit-residue contacts for the full library -- the term hit_num
# is blind to (see hydrophobic_contacts.py). Mirrors
# run_full_library_face_determination.sh: same glob, same env, same worker count.

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/dssp_env
export TMPDIR=/scratch/jem9759/tmp

P=/scratch/jem9759/ZhangWork/KIX_Project
python3 $P/hydrophobic_contacts.py \
    "$P/boltz_out_full/chunk_*/*/predictions/*/*.cif" \
    "$P/full_library_all_metrics/hydrophobic_contacts_full_library.tsv" \
    --workers 16 --name-from parent
