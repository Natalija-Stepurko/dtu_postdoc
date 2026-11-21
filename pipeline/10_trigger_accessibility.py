#!/usr/bin/env python
"""Stage 10. Cleavability: sequence recognition, and reach in the folded construct.

Two questions that a single score would conflate. Whether a protease
recognises the linker sequence is answered by CleaveNet's predictor, as a
contrast between the bite-site proteases and the circulating ones. Whether it
can physically reach the site is answered from the coordinates, against three
necessary conditions, and no sequence model can express it.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import cleavage_context as CC
from nbswitch.occlusion import read_pdb
from nbswitch.config import load_json, save_json, STRUCTURES

lib = load_json("03_constructs.json")
by_linker = {item["id"]: item for item in load_json("02b_linkers.json")["linkers"]}

constructs = list(lib["constructs"])
try:
    constructs += load_json("07_redesign.json")["constructs"]
except FileNotFoundError:
    pass

rows = []
for item in constructs:
    if not item.get("linker_range"):
        continue
    path = STRUCTURES / f"{item['id']}.pdb"
    if not path.exists():
        continue
    context = CC.cleavage_accessibility(read_pdb(path), tuple(item["linker_range"]))
    linker = by_linker.get(item.get("linker_id"), {})
    selectivity = linker.get("selectivity")
    row = {"id": item["id"], "linker_id": item.get("linker_id"),
           "linker_kind": item.get("linker_kind", ""),
           "linker_residues": item["linker_residues"],
           "selectivity": selectivity,
           "target_z": linker.get("target_z"),
           "off_target_z": linker.get("off_target_z"), **context}
    # recognition and reach are separate requirements, so they multiply
    row["cleavability_score"] = (
        None if selectivity is None
        else round(max(selectivity, 0.0) * context["context_score"], 4))
    rows.append(row)

save_json("10_trigger.json", {
    "measures": "CleaveNet selectivity for sequence recognition; extension, "
                "exposure and catalytic-domain clearance for reach in the "
                "folded construct",
    "catalytic_domain_radius_A": CC.CATALYTIC_DOMAIN_RADIUS_A,
    "n": len(rows), "results": rows})

scored = [r for r in rows if r["cleavability_score"] is not None]
print(f"  {len(rows)} constructs scored ({len(scored)} with a generated linker)")
if scored:
    best = max(scored, key=lambda r: r["cleavability_score"])
    print(f"  best: {best['id']} selectivity {best['selectivity']:+.2f} "
          f"x context {best['context_score']:.3f} = {best['cleavability_score']:.4f}")
