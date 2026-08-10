#!/bin/bash
#SBATCH --job-name=dual_prepwizard
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold/logs/dual_prep_%A_%a.out
#SBATCH --array=0-11
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry

# prepwizard on the 172 per-copy sub-complexes, 15 per array task.
# ⚠️ prepwizard drops a <jobname>-001 working dir into the CURRENT directory
# (SCHROD_WORK does NOT redirect it), so cd into a scratch dir first and sweep
# afterwards -- otherwise these accumulate in the project/scratch root.

export SCHROD_WORK=/scratch/jem9759/schrod_work
export TMPDIR=/scratch/jem9759/tmp
WORK=/scratch/jem9759/schrod_work/dual_prep/task_${SLURM_ARRAY_TASK_ID}
mkdir -p "$WORK" && cd "$WORK"

CHUNK=15
START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))
LIST=/scratch/jem9759/ZhangWork/KIX_Project/dual_cofold/copies_list.dat

sed -n "${START},${END}p" $LIST | while read pdb; do
  out="${pdb%.pdb}_prepared.maegz"
  [ -f "$out" ] && continue
  /share/apps/images/run-schrodinger-2025.4.bash run prepwizard \
    -disulfides -nobondorders -rehtreat -noepik -noprotassign -rmsd 0.5 -watdist 0.01 \
    -NOJOBID "$pdb" "$out" >/dev/null 2>&1
done

cd /; rm -rf "$WORK"
