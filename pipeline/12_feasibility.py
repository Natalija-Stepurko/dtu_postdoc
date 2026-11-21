#!/usr/bin/env python
"""Stage 12. Score 01: whether the construct could be made."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import feasibility
from nbswitch.config import load_json, save_json, STRUCTURES

lib = load_json("03_constructs.json")
constructs = list(lib["constructs"])
try:
    constructs += load_json("07_redesign.json")["constructs"]
except FileNotFoundError:
    pass

plddt = {}
for path in STRUCTURES.glob("*.plddt.npy"):
    plddt[path.name.replace(".plddt.npy", "")] = np.load(path)

rows = feasibility.score_library(constructs, STRUCTURES, plddt)
if not rows:
    raise SystemExit("no folded constructs yet; run stage 4 first")

save_json("12_feasibility.json", {
    "components": ["f_fold", "f_preserved", "f_surface", "f_cysteine"],
    "note": "computed from sequence and predicted structure; no measured expression data",
    "n_scored": len(rows), "results": rows})

print(f"\n{'component':>26}  {'mean':>6}  {'min':>6}")
for key in ("f_fold", "f_preserved", "f_surface", "f_cysteine", "feasibility_score"):
    v = [r[key] for r in rows]
    print(f"{key:>26}  {np.mean(v):6.2f}  {min(v):6.2f}")
worst = sorted(rows, key=lambda r: r["feasibility_score"])[:4]
print("\nlowest scoring:")
for r in worst:
    print(f"  {r['id']:22s} {r['feasibility_score']:.2f}   "
          f"unpaired Cys {r['unpaired_cysteines']}  "
          f"exposed hydrophobic {r['exposed_hydrophobic_fraction']:.2f}  "
          f"pLDDT {r['plddt_all']:.0f}")
