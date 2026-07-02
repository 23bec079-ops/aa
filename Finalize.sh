#!/bin/bash
# =========================================================
# finalize.sh — run master_organizer() ONCE, after every array
# task has finished, to build Final_Docking_Results.xlsx from
# the shared outputs/ directory.
#
# Submit it chained to the array job so it only starts once every
# task has completed successfully:
#
#   ARRAY_JOBID=$(sbatch --parsable run_docking.sh)
#   sbatch --dependency=afterok:$ARRAY_JOBID finalize.sh
#
# (If you already ran run_docking.sh separately, just check
# `squeue -u <username>` shows no more batch1_dock tasks, then
# `sbatch finalize.sh` directly.)
# =========================================================

#SBATCH --job-name=batch1_finalize
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --output=logs/finalize_%j.out
#SBATCH --error=logs/finalize_%j.err

set -euo pipefail
mkdir -p logs

module load DL-Conda_3.7
eval "$(conda shell.bash hook)"
conda activate docking
cd "$SLURM_SUBMIT_DIR"

export DOCK_INPUT_XLSX="input.xlsx"   # full file — main() won't run again
export RUN_ORGANIZER=1

# We only want master_organizer(), not a second full main() pass.
# Easiest: call it directly via -c rather than running code.py's
# __main__ block (which would re-run main() over the full input.xlsx —
# harmless since everything is skip-if-done, but slow and pointless).
python -c "
import code
code.master_organizer()
"

echo "Final Excel written: $(pwd)/Final_Docking_Results.xlsx"
