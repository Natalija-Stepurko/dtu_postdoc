#!/usr/bin/env python
"""Stage 13. Constraints, then the nondominated set, then a stated ordering.

Three steps, kept apart because they answer different questions. Constraints
say what cannot work at all. Non-domination says which of the rest are not
beaten outright. The ordering rule says which to look at first, and it is a
choice, so it is named in the output rather than buried in a sort key.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import ranking
from nbswitch.config import load_json, save_json

#: How heavily to discount a score for its own uncertainty. Nothing, for now:
#: no calibrated uncertainty exists for any of the three, and inventing one
#: would make the ordering look better founded than it is.
UNCERTAINTY_PENALTY = 0.0

#: How many of the front to carry forward for verification in stage 19.
N_SHORTLIST = 6

feas = load_json("12_feasibility.json")["results"]
trig = load_json("10_trigger.json")["results"]
try:
    unmask = load_json("04_folding.json")["results"] + \
             load_json("08_redesign_folding.json")["results"]
except FileNotFoundError:
    unmask = []

rows = ranking.combine(feas, trig, unmask)
result = ranking.rank(rows, penalty=UNCERTAINTY_PENALTY, n_shortlist=N_SHORTLIST)

save_json("13_ranking.json", {
    "objectives": result["objectives"],
    "objective_basis": result.get("objective_basis", {}),
    "note": "hard constraints first, then the nondominated set, then maximin on "
            "the normalised objectives. No weights are applied at any step, and "
            "the ordering is a stated rule rather than a property of the front.",
    "constraints": result.get("constraints", []),
    "normalisation": result.get("normalisation", ""),
    "front_order": result.get("front_order", ""),
    "uncertainty_penalty": UNCERTAINTY_PENALTY,
    "n_total": result.get("n_total", 0),
    "n_skipped": result.get("n_skipped", 0),
    "n_excluded": result.get("n_excluded", 0),
    "excluded_by": result.get("excluded_by", {}),
    "n_feasible": result.get("n_feasible", 0),
    "n_front": result.get("n_front", 0),
    "shortlist": result.get("shortlist", []),
    "front": result["front"], "rows": result["rows"]})

objs = result["objectives"]
if not objs:
    raise SystemExit("no construct carries all three objectives; nothing to rank")

print(f"objectives: {', '.join(objs)}")
print(f"{result['n_total']} scored, {result['n_skipped']} missing an objective")
print(f"{result['n_excluded']} excluded by hard constraints"
      + (f" ({result['excluded_by']})" if result["excluded_by"] else ""))
print(f"{result['n_front']} of {result['n_feasible']} feasible designs are "
      f"nondominated\n")

print(f"{'construct':26s} {'maximin':>8s} " + " ".join(f"{o:>18s}" for o in objs))
for r in result["front"][:12]:
    print(f"{r['id']:26s} {r['robust_score']:>8.3f} "
          + " ".join(f"{r[o]:>18.3f}" for o in objs))

print(f"\nshortlist for verification, one per distinct mask: "
      f"{', '.join(result['shortlist'])}")

dom = [r for r in result["rows"] if not r.get("on_front")
       and not r.get("failed_constraints")]
if dom:
    d = max(dom, key=lambda r: r.get("feasibility", 0))
    print(f"\n{len(dom)} feasible designs are dominated. Example: {d['id']} scores "
          + ", ".join(f"{o} {d[o]:.3f}" for o in objs)
          + ", and something beats it on all three.")
