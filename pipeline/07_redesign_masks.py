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

#: Parents are derived from the library rather than named, so renaming or
#: regenerating linkers cannot silently leave this stage redesigning nothing.
#: One construct per mask, on the most selective linker available.
TEMPERATURES = (0.1, 0.2, 0.3)
PER_TEMPERATURE = 8
PYTHON = sys.executable

if not design.available():
    raise SystemExit(f"ProteinMPNN not found at {design.PROTEINMPNN}")

lib = load_json("03_constructs.json")
by_id = {c["id"]: c for c in lib["constructs"]}
vhh = lib["vhh_sequence"]


def _parents(constructs, reference_length=None):
    """One construct per mask: the most selective linker, at the reference length.

    Length variants of a linker carry the same selectivity as the linker they
    came from, so selectivity alone leaves ties and the parent set would depend
    on dictionary order. Redesign holds the linker fixed and rewrites the mask,
    so the parent should be the reference length and the tie is broken on that.
    """
    best = {}
    for c in constructs:
        if not c.get("mask") or not c.get("linker_id"):
            continue
        score = c.get("linker_selectivity")
        score = -1e9 if score is None else score
        at_reference = (reference_length is None
                        or c.get("linker_residues") == reference_length)
        key = (score, at_reference)
        current = best.get(c["mask"])
        if current is None or key > current[0]:
            best[c["mask"]] = (key, c["id"])
    return [cid for _, cid in sorted(best.values(), key=lambda kv: kv[0], reverse=True)]


try:
    _ref_len = load_json("02b_linkers.json").get("reference_length")
except FileNotFoundError:
    _ref_len = None
PARENTS = _parents(lib["constructs"], _ref_len)
if not PARENTS:
    raise SystemExit("no constructs to redesign: the library is empty")
print(f"  redesigning {len(PARENTS)} parents: {', '.join(PARENTS)}")

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
                            "linker_seq": parent["linker_seq"],
                            # the linker is held fixed through the redesign, so
                            # its identity and its scores carry over. Dropping
                            # them left every redesign without a cleavability
                            # score, which the ranking then quietly filled in
                            # from reach alone.
                            "linker_id": parent.get("linker_id"),
                            "linker_kind": parent.get("linker_kind"),
                            "linker_substrate": parent.get("linker_substrate", ""),
                            "linker_selectivity": parent.get("linker_selectivity")})
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
        "linker_id": d.get("linker_id"), "linker_kind": d.get("linker_kind"),
        "linker_substrate": d.get("linker_substrate", ""),
        "linker_selectivity": d.get("linker_selectivity"),
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
