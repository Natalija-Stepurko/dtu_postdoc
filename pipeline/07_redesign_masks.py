#!/usr/bin/env python
"""Stage 7. Rewrite the mask on each folded graft with ProteinMPNN.

The grafted mask is borrowed from the antigen, so nothing about it is tuned to
this paratope. Holding the nanobody and linker fixed and redesigning only the
mask positions turns the graft into a starting point.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import design
from nbswitch.config import load_json, save_json, STRUCTURES, RANDOM_SEED

PARENTS = ["epi1_L5", "epi1_L10", "epi1_L20", "epi2_L5", "epi2_L10", "epi1_half_L5"]
TEMPERATURES = (0.1, 0.2, 0.3)
PER_TEMPERATURE = 8
PYTHON = sys.executable

if not design.available():
    raise SystemExit(f"ProteinMPNN not found at {design.PROTEINMPNN}")

lib = load_json("03_constructs.json")
by_id = {c["id"]: c for c in lib["constructs"]}
vhh = lib["vhh_sequence"]

designs, native_scores = [], {}
for parent_id in PARENTS:
    parent = by_id.get(parent_id)
    pdb = STRUCTURES / f"{parent_id}.pdb"
    if parent is None or not pdb.exists():
        print(f"  skip {parent_id} (not folded yet)")
        continue
    for temp in TEMPERATURES:
        rows = design.design_mask(pdb, parent["mask_range"],
                                  STRUCTURES, n_sequences=PER_TEMPERATURE,
                                  temperature=temp, seed=RANDOM_SEED, python=PYTHON)
        for row in rows:
            mask_seq = row["sequence"][parent["mask_range"][0] - 1:parent["mask_range"][1]]
            if row["is_native"]:
                native_scores[parent_id] = row["score"]
                continue
            designs.append({"parent": parent_id, "temperature": temp,
                            "mask_seq": mask_seq, "score": row["score"],
                            "global_score": row["global_score"],
                            "linker_residues": parent["linker_residues"],
                            "linker_seq": parent["linker_seq"]})
    print(f"  {parent_id}: native {native_scores.get(parent_id, float('nan')):.3f}, "
          f"{len([d for d in designs if d['parent'] == parent_id])} designs")

# keep the best-scoring instance of each distinct mask sequence per parent
best = {}
for d in designs:
    key = (d["parent"], d["mask_seq"])
    if key not in best or d["score"] < best[key]["score"]:
        best[key] = d
unique = sorted(best.values(), key=lambda d: d["score"])

constructs = []
for i, d in enumerate(unique):
    mask_end = len(d["mask_seq"])
    n = d["linker_residues"]
    constructs.append({
        "id": f"mpnn{i:03d}_{d['parent']}", "mask": f"mpnn_{d['parent']}",
        "mask_seq": d["mask_seq"], "mask_kind": "ProteinMPNN redesign",
        "parent": d["parent"], "temperature": d["temperature"],
        "mpnn_score": d["score"], "native_score": native_scores.get(d["parent"]),
        "linker_residues": n, "linker_seq": d["linker_seq"],
        "sequence": d["mask_seq"] + d["linker_seq"] + vhh,
        "length": mask_end + n + len(vhh),
        "mask_range": [1, mask_end], "linker_range": [mask_end + 1, mask_end + n],
        "vhh_range": [mask_end + n + 1, mask_end + n + len(vhh)]})

save_json("07_redesign.json", {
    "method": "ProteinMPNN, nanobody and linker fixed",
    "parents": PARENTS, "temperatures": list(TEMPERATURES),
    "per_temperature": PER_TEMPERATURE, "native_scores": native_scores,
    "n_designs": len(designs), "n_unique": len(unique),
    "constructs": constructs})

print(f"\n{len(designs)} designs, {len(unique)} distinct mask sequences")
improved = [d for d in unique if d["score"] < native_scores.get(d["parent"], 9e9)]
print(f"{len(improved)} score better than their parent graft")
for d in unique[:8]:
    print(f"  {d['mask_seq']:<10s} from {d['parent']:<14s} T={d['temperature']}  "
          f"score {d['score']:.3f}  (native {native_scores.get(d['parent'], float('nan')):.3f})")
