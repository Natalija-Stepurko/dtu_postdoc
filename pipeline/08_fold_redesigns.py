#!/usr/bin/env python
"""Stage 8. Fold the redesigned constructs and score them the same way."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import folding, occlusion as O
from nbswitch.structure import residue_sasa
from nbswitch.config import load_json, save_json, STRUCTURES

lib = load_json("03_constructs.json")
redesign = load_json("07_redesign.json")
constructs = redesign["constructs"]

plddt = folding.fold_library(constructs, STRUCTURES, device=None)

para_idx = O.paratope_indices(lib["vhh_residue_labels"], lib["paratope_labels"])
free_sasa = float(residue_sasa(O.read_pdb(STRUCTURES / "VHH_only.pdb"))[para_idx].sum())

rows = []
for c in constructs:
    row = O.occlusion(c, STRUCTURES, para_idx, free_sasa, plddt)
    row["vhh_core_rmsd_A"] = O.vhh_core_rmsd(c, STRUCTURES)
    row["mpnn_score"] = c["mpnn_score"]
    row["parent"] = c["parent"]
    row["mask_seq"] = c["mask_seq"]
    rows.append(row)

save_json("08_redesign_folding.json", {
    "model": folding.MODEL_NAME, "n_folded": len(rows),
    "free_vhh_paratope_sasa_A2": round(free_sasa, 1), "results": rows})

occ = [r["occlusion"] * 100 for r in rows]
print(f"\n{len(rows)} redesigns folded")
print(f"  occlusion mean {np.mean(occ):+.1f}%  best {max(occ):+.1f}%  worst {min(occ):+.1f}%")
best = max(rows, key=lambda r: r["occlusion"])
print(f"  best: {best['id']}  mask {best['mask_seq']}  "
      f"{best['occlusion']*100:+.1f}%  pLDDT {best['plddt_vhh']:.0f}")
