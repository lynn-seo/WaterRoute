#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python}"
CLAIM_FILE="${TASKS_SA_CLAIM_FILE:-${SCRIPT_DIR}/.TASKS_SA_claimed.jsonl}"
LOCK_FILE="${CLAIM_FILE}.lock"

mkdir -p "$(dirname -- "${CLAIM_FILE}")"
exec 9>"${LOCK_FILE}"
flock -x 9

SELECTION="$(
    "${PYTHON}" - "${SCRIPT_DIR}" "${CLAIM_FILE}" <<'PY'
import json
import os
import random
import sys
from pathlib import Path

script_dir = Path(sys.argv[1])
claim_file = Path(sys.argv[2])
sys.path.insert(0, str(script_dir))

import TASKS_SA

tasks = TASKS_SA.TASKS
claimed = set()
if claim_file.exists():
    with claim_file.open(encoding="utf-8") as claims:
        for line_number, line in enumerate(claims, start=1):
            try:
                claim = json.loads(line)
                claimed.add(claim["task_key"])
            except (json.JSONDecodeError, KeyError) as error:
                raise ValueError(
                    f"Invalid task claim at {claim_file}:{line_number}"
                ) from error

available = [
    (index, task, json.dumps(task, separators=(",", ":")))
    for index, task in enumerate(tasks)
    if json.dumps(task, separators=(",", ":")) not in claimed
]
if not available:
    raise SystemExit("All tasks in TASKS_SA.TASKS have already been claimed.")

index, task, task_key = random.SystemRandom().choice(available)
claim = {"index": index, "task": task, "task_key": task_key}
with claim_file.open("a", encoding="utf-8") as claims:
    claims.write(json.dumps(claim, separators=(",", ":")) + "\n")
    claims.flush()
    os.fsync(claims.fileno())

print(index)
print(json.dumps(task))
PY
)"

mapfile -t SELECTION_LINES <<< "${SELECTION}"
TASK_INDEX="${SELECTION_LINES[0]}"
TASK_DETAILS="${SELECTION_LINES[1]}"

echo "Submitting task ${TASK_INDEX}: ${TASK_DETAILS}"
if ! JOB_ID="$(
    sbatch --parsable \
        --array="${TASK_INDEX}-${TASK_INDEX}" \
        "${SCRIPT_DIR}/submit_SA.slurm"
)"; then
    echo "sbatch failed; releasing task ${TASK_INDEX} claim." >&2
    "${PYTHON}" - "${CLAIM_FILE}" "${TASK_INDEX}" <<'PY'
import json
import sys
from pathlib import Path

claim_file = Path(sys.argv[1])
task_index = int(sys.argv[2])
lines = claim_file.read_text(encoding="utf-8").splitlines()
kept = []
removed = False
for line in lines:
    claim = json.loads(line)
    if not removed and claim["index"] == task_index:
        removed = True
    else:
        kept.append(line)
if not removed:
    raise RuntimeError(f"Could not release task claim {task_index} from {claim_file}")
claim_file.write_text("\n".join(kept) + "\n", encoding="utf-8")
PY
    exit 1
fi

echo "Submitted Slurm job ${JOB_ID} for task ${TASK_INDEX}."
