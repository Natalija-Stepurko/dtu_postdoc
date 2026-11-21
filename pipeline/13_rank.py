#!/usr/bin/env python
"""Stage 13. Combine the scores and take the Pareto front."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import ranking
from nbswitch.config import load_json, save_json

feas = load_json("12_feasibility.json")["results"]
trig = load_json("10_trigger.json")["results"]
try:
    unmask = load_json("04_folding.json")["results"] + \
             load_json("08_redesign_folding.json")["results"]
except FileNotFoundError:
    unmask = []

rows = ranking.combine(feas, trig, unmask)
result = ranking.rank(rows)

save_json("13_ranking.json", {
    "objectives": result["objectives"],
    "note": "a design is kept unless another is at least as good on every objective "
            "and better on one; no weights are applied",
    "n_total": result["n_total"], "n_front": result["n_front"],
    "front": result["front"], "rows": result["rows"]})

print(f"objectives used: {', '.join(result['objectives'])}")
print(f"{result['n_front']} of {result['n_total']} designs are on the Pareto front\n")
print(f"{'construct':24s} " + " ".join(f"{o:>12s}" for o in result["objectives"]))
for r in result["front"][:12]:
    print(f"{r['id']:24s} " + " ".join(f"{r[o]:>12.3f}" for o in result["objectives"]))
dom = [r for r in result["rows"] if not r["on_front"]]
if dom:
    print(f"\n{len(dom)} designs are dominated. Example:")
    d = max(dom, key=lambda r: r.get("feasibility", 0))
    print(f"  {d['id']} scores " +
          ", ".join(f"{o} {d[o]:.3f}" for o in result["objectives"]) +
          " and is beaten on every one by something else.")
