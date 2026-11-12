#!/bin/bash
#SBATCH --partition=mmb_cpu_zen3,mmb_cpu_sphr
#SBATCH --job-name=nanobodies_array
#SBATCH --output=array_logs/mwf_%A_%a.out
#SBATCH --error=array_logs/mwf_%A_%a.err
#SBATCH --array=0-219%30
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --time=24:00:00
#SBATCH --mail-user=<your_email>

# Submit from the dataset root (the folder holding nanobodies.db and vhh_run/):
#   mkdir -p array_logs && sbatch job_array.sh
# Override the range/concurrency at submit time, e.g. a small pilot:
#   sbatch --array=0-2%1 job_array.sh
# Total directories: 220 -> array is 0-219. Recount with: ls -d vhh_run/*/ | wc -l

# `module` is not available on the IRB compute nodes (see references/irb-cluster.md)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate mwf_env

# Store the submission directory (dataset root)
SUBMIT_DIR=$(pwd)
UTILS=$SUBMIT_DIR/utils

# Build the list of project directories (must match the --array range length)
DIRS=($(ls -d vhh_run/*/))

# Pick this task's project directory
DIR=${DIRS[$SLURM_ARRAY_TASK_ID]}
echo "Processing directory: $DIR"
DIR_NAME=$(basename "$DIR")

# Move into the project directory
cd "$DIR" || exit 1

# Skip projects already finished, so re-submitting the array is cheap
STATE=$(mwf dataset status ../../nanobodies.db -p .)
CURRENT_STATE=$(echo "$STATE" | grep "^State:" | awk '{print $2}')
echo "Current state: $CURRENT_STATE"
if [ "$CURRENT_STATE" == "done" ]; then
    echo "Skipping $DIR (already done)."
    exit 0
fi

# (Optional) regenerate this project's inputs.yaml just before running.
# Skip if you already ran `mwf dataset inputs ...` over the whole dataset beforehand.
# mwf dataset inputs ../../nanobodies.db -it "$UTILS/inputs_template.yaml" -ig "$UTILS/nanobodies_generator.py" -p "*${DIR_NAME}" -o || exit 1

# Run the workflow
# -e: skip the slowest analyses; -lb: local BLAST (see references/launching.md)
mwf run -img -fit -e energies clusters pockets -lb

# Copy the SLURM logs into the project's own logs directory
mkdir -p logs
cp "${SUBMIT_DIR}/array_logs/mwf_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out" logs/ 2>/dev/null
cp "${SUBMIT_DIR}/array_logs/mwf_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" logs/ 2>/dev/null
