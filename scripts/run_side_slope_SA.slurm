#!/bin/bash
#SBATCH --job-name=waterroute-side-slope
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=OD-226458

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Each array task gets its own working directory because routing writes
# intermediate pickle files with fixed names in the current directory.
TASK_TMP_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TASK_WORK_DIR="${TASK_TMP_ROOT}/waterroute-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASK_WORK_DIR}"
cd "${TASK_WORK_DIR}"

# Load/activate the Python environment using your cluster's setup if needed.
# For conda, for example:
#   source "${HOME}/miniforge3/etc/profile.d/conda.sh"
#   conda activate channelroute

PYTHON="${PYTHON:-python}"
exec "${PYTHON}" "${SCRIPT_DIR}/run_side_slope_sensitivity.py" \
    --task-index "${SLURM_ARRAY_TASK_ID}"
