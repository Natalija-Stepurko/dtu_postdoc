#!/usr/bin/env bash
# Every stage, in the order it runs. Stage numbers are execution order: a
# stage that reads another stage's output has the higher number, which is why
# verification is 19 and not 11 - it consumes the ranking from 13.
#
# The folding stages (04, 08, 16) and the structure predictions (09, 17, 19)
# are the expensive ones. Folding is resumable and keyed on the sequence in
# each file, so the script can be interrupted and restarted without losing
# work or serving a stale model.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python}"
export HF_HOME="${HF_HOME:-$PWD/.hf}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$(nproc)}"

stages=(
  01_reference_structure     # the complex: paratope, termini, contacts
  02_tether_model            # what mask affinity both states allow
  02b_linkers                # design linkers: cleavable by construction
  03a_mask_design            # design masks: seeded from the contact face
  03_construct_library       # assemble: every mask on every linker
  04_fold_constructs         # fold, once, shared by all three scores
  16_linker_length_sweep     # does the mask reach, and at what length
  07_redesign_masks          # redesign the mask on the folded backbone
  08_fold_redesigns          # fold the redesigns
  09_complex_prediction      # the alignment control, kept as a diagnostic
  10_trigger_accessibility   # score: can a protease reach the linker
  12_feasibility             # score: can the construct be made
  13_rank                    # the nondominated set over the three scores
  14_ensemble                # all designs on one backbone, for the viewer
  15_interface_controls      # the same measures on independent complexes
  17_mask_binding            # does the mask bind untethered
  18_cleavage_evidence       # the linkers against the measured library
  19_verify_front            # predict the whole assembly for the front
  05_figures
  06_report
)

for stage in "${stages[@]}"; do
  echo "=== $stage ==="
  "$PY" "pipeline/${stage}.py"
done
echo "=== done ==="
