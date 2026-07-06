#!/bin/bash
#SBATCH --job-name=build_helices_full
#SBATCH --output=build_helices_full_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry

export SCHROD_WORK=/scratch/jem9759/schrod_work
export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$SCHROD_WORK" "$TMPDIR"

cd /home/jem9759/ZhangLabWork/KIX_Project   # so relative paths resolve

/share/apps/images/run-schrodinger-2025.4.bash run python3 \
    make_helix_schrodinger.py ./full_library_out.tsv
