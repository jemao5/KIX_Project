#!/bin/bash
#SBATCH --job-name=schrodinger_hbond_calc
#SBATCH --array=0-123
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/logs/schrodinger_hbond_calc_%A_%a.out

export SCHROD_WORK=/scratch/jem9759/schrod_work
export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$SCHROD_WORK" "$TMPDIR"
mkdir -p /scratch/jem9759/ZhangWork/KIX_Project/schrodinger_calc_hbond_chunks

CHUNK_SIZE=250
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
END=$(( START + CHUNK_SIZE - 1 ))
list_file="/scratch/jem9759/ZhangWork/KIX_Project/binding_face_filtered_hits.dat"
out_pkl="/scratch/jem9759/ZhangWork/KIX_Project/schrodinger_calc_hbond_chunks/schrodinger_calc_hbond_chunk_${SLURM_ARRAY_TASK_ID}.pkl"

/share/apps/images/run-schrodinger-2025.4.bash run python3 \
    /scratch/jem9759/ZhangWork/KIX_Project/schrodinger_calc_hbond.py \
    "${list_file}" "${out_pkl}" --start ${START} --end ${END}


