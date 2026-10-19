
#!/bin/bash
#SBATCH --job-name=waterroute-side-slope
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-10
#SBATCH --account=OD-239352

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Each array task gets its own working directory
TASK_TMP_ROOT="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
TASK_WORK_DIR="${TASK_TMP_ROOT}/waterroute-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
mkdir -p "${TASK_WORK_DIR}"
cd "${TASK_WORK_DIR}"

# Activate your Python environment here if needed.
# Example:
# source "${HOME}/miniforge3/etc/profile.d/conda.sh"
# conda activate channelroute

PYTHON="${PYTHON:-python}"

exec "${PYTHON}" "${SCRIPT_DIR}/run_side_slope_sensitivity.py" \
    --task-index "${SLURM_ARRAY_TASK_ID}"
