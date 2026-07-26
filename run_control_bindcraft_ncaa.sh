#!/bin/bash
#SBATCH --job-name=bindcraft_ncaa
#SBATCH --output=/scratch/jem9759/ZhangWork/KIX_Project/control_run/logs/bindcraft_ncaa_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --account=torch_pr_149_chemistry
#SBATCH --qos=cpu48
# NO GPU (CPU-only)

# Stage 6 for the ncAA-swapped literature controls.
#
# Differences from run_control_bindcraft.sh:
#   --input points at control_run/ncaa_structures/ (PDB, post-ResidueX) rather
#     than the raw Boltz CIFs.
#   --extra-res-fa supplies six params files. Only BCS is used straight from the
#     Rosetta database. The other five (M2F/CHX/CY5/ABU/NLE) are LOCAL renamed
#     copies of the database entries with the NCAA_ROTLIB_* lines stripped:
#     the shipped params request rotamer libraries under
#     ncaa_rotamer_libraries/alpha_amino_acid/, but PyRosetta ships that dir
#     EMPTY, so FastRelax died with "Could not open rotamer library file" on 17
#     of 24 controls. They are renamed because a duplicate NAME would collide
#     with the database copy. CIA is our derived S-carboxymethyl-cysteine.
#   --strict-residues drops -ignore_unrecognized_res so a missing ResidueType
#     is a loud crash rather than the residue being silently deleted from the
#     pose (which is what happened before the naming fix, and would have
#     produced plausible-looking but wrong interface metrics).

export TMPDIR=/scratch/jem9759/tmp
mkdir -p "$TMPDIR"

cd /scratch/jem9759/ZhangWork/BindCraft   # repo dir for functions/ + settings_filters/

module purge
module load anaconda3/2025.06
source activate /scratch/jem9759/envs/BindCraft

P=/scratch/jem9759/ZhangWork/KIX_Project
export PYTHONPATH=/scratch/jem9759/ZhangWork/BindCraft:$PYTHONPATH

python -u "$P/score_peptide.py" \
    --input "$P/control_run/ncaa_structures/*_ncaa.pdb" \
    --binder-chain B --target-chain A \
    --filters settings_filters/peptide_filters.json \
    --dalphaball-path functions/DAlphaBall.gcc \
    --extra-res-fa "$P/control_run/rosetta_params/M2F.params" \
    --extra-res-fa "$P/control_run/rosetta_params/CHX.params" \
    --extra-res-fa "$P/control_run/rosetta_params/CY5.params" \
    --extra-res-fa "$P/control_run/rosetta_params/ABU.params" \
    --extra-res-fa "$P/control_run/rosetta_params/NLE.params" \
    --extra-res-fa "$P/control_run/rosetta_params/CIA.params" \
    --strict-residues \
    --output-json "$P/control_run/bindcraft_ncaa_out" \
    --output-relaxed-pdb "$P/control_run/bindcraft_ncaa_relaxed"
