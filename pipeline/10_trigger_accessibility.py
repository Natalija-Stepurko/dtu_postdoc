#!/usr/bin/env python
"""Stage 10. The half of the trigger score that does not need specificity data.

A protease must thread its substrate through the active-site cleft in an
extended conformation. A linker folded against the nanobody cannot be cut
whatever its sequence, so accessibility is a necessary condition for the
trigger and it is computable from the folded constructs alone.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import trigger
from nbswitch.config import load_json, save_json, STRUCTURES

lib = load_json("03_constructs.json")
rows = trigger.score_library(lib["constructs"], STRUCTURES)

try:
    rows += trigger.score_library(load_json("07_redesign.json")["constructs"], STRUCTURES)
except FileNotFoundError:
    pass

save_json("10_trigger.json", {
    "measures": "solvent accessibility and extension of the linker in the folded construct",
    "exposed_threshold_A2": trigger.EXPOSED_A2,
    "note": "a necessary condition for cleavage, not a sufficient one; sequence "
            "recognition and selectivity need specificity data for the proteases in question",
    "n_scored": len(rows), "results": rows})

if not rows:
    raise SystemExit("no folded constructs with linkers yet - run stage 4 first")

by_len = {}
for r in rows:
    by_len.setdefault(r["linker_residues"], []).append(r)
print(f"\n{'linker':>7} {'n':>4} {'exposed':>9} {'extension':>10} {'score':>7}")
for n in sorted(by_len):
    g = by_len[n]
    print(f"{n:>7} {len(g):>4} "
          f"{np.mean([x['fraction_exposed'] for x in g]):>9.2f} "
          f"{np.mean([x['extension_ratio'] for x in g]):>10.2f} "
          f"{np.mean([x['accessibility_score'] for x in g]):>7.2f}")
best = max(rows, key=lambda r: r["accessibility_score"])
worst = min(rows, key=lambda r: r["accessibility_score"])
print(f"\nmost accessible : {best['id']:22s} {best['accessibility_score']:.2f}")
print(f"least accessible: {worst['id']:22s} {worst['accessibility_score']:.2f}")
