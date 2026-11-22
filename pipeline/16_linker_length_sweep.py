#!/usr/bin/env python
"""Stage 16. Does the mask fail to occlude because it cannot reach?

The library fixes every linker at sixteen residues, which is not a design
choice: CleaveNet emits ten-residue substrates and the flanks are three a side.
Length and cleavability are therefore the same number, and the masking result
cannot be attributed to either.

This stage separates them. The substrate is held constant and the flanks are
lengthened, giving total linkers of 16 to 32 residues that are equally
cleavable by construction. Two linkers are swept: the most selective generated
one, whose substrate is hydrophobic, and the Gly/Ser control, which is not. If
reach is the limit, occlusion rises with length. If the limit is that the mask
has no affinity for the paratope, the designed masks gain no more than the
inert control does, at any length.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from nbswitch import constructs as C, folding, linkers as L, occlusion as O
from nbswitch.config import load_json, save_json, STRUCTURES
from nbswitch.structure import residue_sasa

lib = load_json("03_constructs.json")
all_linkers = load_json("02b_linkers.json")["linkers"]


def _pick(predicate, what):
    hit = next((x for x in all_linkers if predicate(x)), None)
    if hit is None:
        raise SystemExit(f"no {what} in 02b_linkers.json")
    return hit


generated = max((x for x in all_linkers if x.get("kind") == "generated"),
                key=lambda x: x["selectivity"])
gs_control = _pick(lambda x: x["id"] == "ctrl_gs", "Gly/Ser control linker")

series = []
for base in (generated, gs_control):
    series += L.length_variants(base["substrate"], base["id"],
                                kind=base.get("kind", "generated"),
                                selectivity=base.get("selectivity"))

masks = {name: m for name, m in lib["masks"].items()}
library = [c for c in C.build_library(lib["vhh_sequence"], masks, series)
           if c["id"] != "VHH_only"]
print(f"{len(masks)} masks x {len(series)} linkers = {len(library)} constructs "
      f"over lengths {sorted({v['total_length'] for v in series})}")

plddt = folding.fold_library(library, STRUCTURES, device=None)

para_idx = O.paratope_indices(lib["vhh_residue_labels"], lib["paratope_labels"])
free = O.read_pdb(STRUCTURES / "VHH_only.pdb")
free_sasa = float(residue_sasa(free)[para_idx].sum())

rows = []
for c in library:
    row = O.occlusion(c, STRUCTURES, para_idx, free_sasa, plddt)
    row["linker_residues"] = c["linker_residues"]
    row["linker_id"] = c["linker_id"]
    row["substrate_kind"] = ("generated" if c["linker_id"].startswith(generated["id"])
                             else "Gly/Ser control")
    rows.append(row)

save_json("16_length_sweep.json", {
    "question": "does occlusion improve when the mask is given more reach",
    "substrate_held_constant": True,
    "lengths": sorted({c["linker_residues"] for c in library}),
    "free_vhh_paratope_sasa_A2": round(free_sasa, 1),
    "n_constructs": len(rows), "results": rows})

print("\nmean occlusion (%), by linker length and mask kind")
kinds = sorted({r["kind"] for r in rows})
print(f"  {'length':>8s}  {'substrate':<16s}" + "".join(f"{k:>20s}" for k in kinds))
for n in sorted({r["linker_residues"] for r in rows}):
    for sub in ("generated", "Gly/Ser control"):
        sel = [r for r in rows if r["linker_residues"] == n
               and r["substrate_kind"] == sub]
        if not sel:
            continue
        cells = ""
        for k in kinds:
            v = [r["occlusion"] * 100 for r in sel if r["kind"] == k]
            cells += f"{np.mean(v):+19.2f}%" if v else f"{'-':>20s}"
        print(f"  {n:>6d}    {sub:<16s}" + cells)

print("\ndoes any length let the grafts beat the inert control?")
for n in sorted({r["linker_residues"] for r in rows}):
    sel = [r for r in rows if r["linker_residues"] == n]
    graft = [r["occlusion"] * 100 for r in sel if r["kind"] == "epitope graft"]
    inert = [r["occlusion"] * 100 for r in sel if r["kind"] == "inert control"]
    if graft and inert:
        gap = np.mean(graft) - np.mean(inert)
        print(f"  {n:2d} residues: graft {np.mean(graft):+.2f}%  "
              f"inert {np.mean(inert):+.2f}%  gap {gap:+.2f} "
              f"{'GRAFT WINS' if gap > 0 else 'control wins'}")
