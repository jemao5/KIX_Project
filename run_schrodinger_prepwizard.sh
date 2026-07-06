#!/bin/bash

#SBATCH --job-name=pwdizard_pdbbind
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-123
#SBATCH --mem=8GB
#SBATCH --time=4:00:00
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --output=logs/prepwizard_%A_%a.out

export SCHROD_WORK=/scratch/jem9759/schrod_work
export TMPDIR=/scratch/jem9759/tmp

CHUNK_SIZE=250
START=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE + 1 ))
END=$(( START + CHUNK_SIZE - 1 ))

list_file="/scratch/jem9759/ZhangWork/KIX_Project/binding_face_filtered_hits.dat"

sed -n "${START},${END}p" ${list_file} | while read pdb_path; do
    /share/apps/images/run-schrodinger-2025.4.bash run prepwizard \
    -disulfides -nobondorders -rehtreat -noepik -noprotassign -rmsd 0.5 -watdist 0.01 \
    -NOJOBID \
    ${pdb_path} ${pdb_path%.pdb}_prepared.maegz
done




