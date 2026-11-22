#!/usr/bin/env python
"""Stage 18. Check the cleavability score against the measurements, not the model.

CleaveNet proposes the substrates and CleaveNet scores them. Both heads are
trained on the same measured library, so their agreement says only that the
model is self-consistent. This stage scores the same substrates against that
library directly, by two methods that are not neural networks: a position
weight matrix built from the measured extremes, and a lookup of what the most
similar measured substrates were measured to do.

The measurements are the ones CleaveNet learned from, so this is an independent
estimator rather than independent evidence. It can show that the neural scores
are recoverable from the data by simpler means, which is worth knowing. It
cannot show that any of it transfers to a snake venom metalloproteinase.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from nbswitch import cleavage_evidence as E
from nbswitch.config import load_json, save_json

linkers = load_json("02b_linkers.json")["linkers"]
calibration = E.calibrate()
rows = E.evidence([x["substrate"] for x in linkers])

frame = E.measured_library()
for x, r in zip(linkers, rows):
    r["id"] = x["id"]
    r["kind"] = x.get("kind", "generated")
    r["cleavenet_selectivity"] = x.get("selectivity")
    r["measured_percentile"] = round(float(
        (frame.measured_selectivity < (x.get("selectivity") or -9)).mean()) * 100, 1)

generated = [r for r in rows if r["kind"] == "generated"]
controls = [r for r in rows if r["kind"] != "generated"]


def _mean(items, key):
    vals = [i[key] for i in items if i.get(key) is not None]
    return float(np.mean(vals)) if vals else float("nan")


agreement = float("nan")
pairs = [(r["cleavenet_selectivity"], r["pwm_score"]) for r in rows
         if r.get("cleavenet_selectivity") is not None]
if len(pairs) > 2:
    a, b = zip(*pairs)
    agreement = float(np.corrcoef(a, b)[0, 1])

save_json("18_cleavage_evidence.json", {
    "question": "do the generated substrates hold up against the measured library",
    "library": str(E.MEASURED_CSV), "n_measured": int(len(frame)),
    "calibration": calibration,
    "cleavenet_vs_pwm_r": round(agreement, 3),
    "results": rows})

print(f"measured library: {len(frame)} substrates")
print(f"PWM recovers measured selectivity out of sample at r = "
      f"{calibration['pwm_vs_measured_r']:+.3f} (n={calibration['n_test']})")
print(f"CleaveNet against the PWM across our linkers: r = {agreement:+.3f}\n")
print(f"{'id':<11s} {'substrate':<11s} {'CleaveNet':>9s} {'PWM':>7s} "
      f"{'neighbours':>11s}  nearest measured")
for r in rows:
    print(f"{r['id']:<11s} {r['substrate']:<11s} "
          f"{r['cleavenet_selectivity']:+9.2f} {r['pwm_score']:+7.2f} "
          f"{r['neighbour_measured_selectivity']:+11.2f}  "
          f"{r['best_neighbour']} ({r['best_neighbour_selectivity']:+.2f})")

print(f"\ngenerated   PWM {_mean(generated, 'pwm_score'):+.2f}  "
      f"neighbours {_mean(generated, 'neighbour_measured_selectivity'):+.2f}")
print(f"controls    PWM {_mean(controls, 'pwm_score'):+.2f}  "
      f"neighbours {_mean(controls, 'neighbour_measured_selectivity'):+.2f}")
if any(r["in_measured_library"] for r in rows):
    print("\nWARNING: a generated substrate is present verbatim in the measured "
          "library, so it was reproduced rather than designed")
else:
    print("\nNo generated substrate appears verbatim in the measured library.")
