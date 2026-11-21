#!/usr/bin/env bash
# Run every stage in order. Stage 4 and 8 fold with ESMFold and are resumable,
# so the script can be interrupted and restarted without losing work.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python}"
export HF_HOME="${HF_HOME:-$PWD/.hf}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"

for stage in 01_reference_structure 02_tether_model 03_construct_library \
             04_fold_constructs 07_redesign_masks 08_fold_redesigns \
             05_figures 06_report; do
  echo "=== $stage ==="
  "$PY" "pipeline/${stage}.py"
done
echo "=== done ==="
