#!/usr/bin/env python
"""Stage 4. Fold every construct and measure paratope occlusion.

Runs on CPU; a GPU only makes it faster. Resumable - already-folded constructs
are reused, so the stage can be interrupted and restarted.
"""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from nbswitch import folding, occlusion as O
from nbswitch.config import load_json, save_json, STRUCTURES

lib = load_json("03_constructs.json")
library = lib["constructs"]
plddt = folding.fold_library(library, STRUCTURES, device=None)

para_idx = O.paratope_indices(lib["vhh_residue_labels"], lib["paratope_labels"])
free = O.read_pdb(STRUCTURES / "VHH_only.pdb")
free_sasa = float(__import__("nbswitch.structure", fromlist=["residue_sasa"])
                  .residue_sasa(free)[para_idx].sum())

rows = []
for c in library:
    if c["id"] == "VHH_only":
        continue
    row = O.occlusion(c, STRUCTURES, para_idx, free_sasa, plddt)
    row["vhh_core_rmsd_A"] = O.vhh_core_rmsd(c, STRUCTURES)
    rows.append(row)

save_json("04_folding.json", {
    "model": folding.MODEL_NAME, "note": "single-sequence prediction, no MSA",
    "free_vhh_paratope_sasa_A2": round(free_sasa, 1),
    "free_vhh_plddt": round(float(plddt["VHH_only"].mean()), 2),
    "n_folded": len(library), "results": rows})

print(f"\nfree VHH paratope SASA {free_sasa:.1f} A2, "
      f"pLDDT {plddt['VHH_only'].mean():.1f}")
by_kind = {}
for r in rows:
    by_kind.setdefault(r["kind"], []).append(r["occlusion"])
print("\nocclusion by mask kind")
for kind, vals in sorted(by_kind.items()):
    print(f"  {kind:20s} n={len(vals):2d}  mean {np.mean(vals):+.3f}  "
          f"range {min(vals):+.3f} to {max(vals):+.3f}")
