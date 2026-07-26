#!/bin/bash
#SBATCH --job-name=bindcraft_controls
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/control_run/logs/bindcraft_controls_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --qos=cpu48
# NO GPU (CPU-only)

# Stage 6 for the literature control set. Same flags as run_bindcraaft_score.sh,
# just pointed at control_run/ instead of the chunked full-library output.

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

cd /scratch/jem9759/ZhangWork/BindCraft   # repo dir for functions/ + settings_filters/

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/BindCraft

P=/scratch/jem9759/ZhangWork/KIX_Project

# score_peptide.py lives in KIX_Project (the BindCraft repo copy is named
# score_peptide_copy.py); cwd stays in the repo so functions/ and
# settings_filters/ resolve relatively. Running it by absolute path puts
# KIX_Project on sys.path instead of the repo, so `import functions` needs
# PYTHONPATH pointed back at the repo.
export PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft:$PYTHONPATH

python -u "$P/score_peptide.py" \
    --input "$P/control_run/boltz_out/*/predictions/*/*_model_0.cif" \
    --binder-chain B --target-chain A \
    --filters settings_filters/peptide_filters.json \
    --dalphaball-path functions/DAlphaBall.gcc \
    --output-json "$P/control_run/bindcraft_out" \
    --output-relaxed-pdb "$P/control_run/bindcraft_relaxed"
