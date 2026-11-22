#!/usr/bin/env python
"""Stage 02b. Generate protease-cleavable linkers, before any construct exists.

Cleavability is a constraint here, not a filter later: CleaveNet is asked for
substrates that the bite-site proteases cut and the circulating ones do not, so
every linker entering the library is cleavable by construction.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import linkers as L
from nbswitch.config import save_json, RESULTS

N_GENERATE = 120

#: How many generated substrates to carry into the library. The substrates
#: are near-equivalent on selectivity, so each extra one multiplies the
#: construct count by a factor the scores cannot use. The mask axis is where
#: the design space is; this one only has to be represented.
N_KEEP = 4
PYTHON = "/scratch/.venv-dtu/bin/python"

#: Linker length is a design variable, not a constant. The reference length
#: comes from the flank setting and is the one every substrate is built at, so
#: that sequences can be compared with length held constant. A ladder is then
#: added for one generated substrate and one control, so that length varies
#: with sequence held constant and the ranking can trade the two against each
#: other.
#:
#: It matters because the masks are not all the same size. A six-residue mask
#: and a twenty-residue one attach at the same point and have to cover the same
#: surface, so they do not need the same amount of linker, and a single length
#: for the whole library is a guess for at least one of them.
LENGTH_LADDER = tuple(n for n in L.LENGTH_SERIES if n != L.REFERENCE_LENGTH)

work = RESULTS / "linkers"
out = L.cleavable_linkers(n_generate=N_GENERATE, n_keep=N_KEEP,
                          work_dir=work, python=PYTHON)

# two controls, so the generated ones have something to be better than
best = out["linkers"][0]["substrate"]
import random
rng = random.Random(42)
scrambled = "".join(rng.sample(best, len(best)))
# Controls carry a substrate of the same length as the generated ones and are
# built the same way, so every linker in the library is the same size and only
# the sequence differs.
controls = [
    {"id": "ctrl_gs", "substrate": "GGGGSGGGGS",
     "sequence": L.build_linker("GGGGSGGGGS"),
     "kind": "protease-resistant control", "selectivity": None,
     "target_z": None, "off_target_z": None},
    {"id": "ctrl_scram", "substrate": scrambled,
     "sequence": L.build_linker(scrambled),
     "kind": "scrambled control", "selectivity": None,
     "target_z": None, "off_target_z": None},
]
for item in out["linkers"]:
    item["kind"] = "generated"
out["linkers"] = out["linkers"] + controls
out["n_controls"] = len(controls)

# score the controls with the same predictor, so the comparison is like for like
frame = L.selectivity(L.predict([c["substrate"] for c in controls],
                                work_dir=work / "ctrl", python=PYTHON))
for c, row in zip(controls, frame.itertuples()):
    c["target_z"] = round(float(row.target_z), 3)
    c["off_target_z"] = round(float(row.off_target_z), 3)
    c["selectivity"] = round(float(row.selectivity), 3)

# length variants, so the library carries sequence and length as separate axes
ladder = []
_ladder_bases = [out["linkers"][0]] + controls          # best generated + both controls
for base in _ladder_bases:
    for v in L.length_variants(base["substrate"], base["id"], lengths=LENGTH_LADDER):
        if len(v["sequence"]) == len(base["sequence"]):
            continue                      # the reference length is already in
        ladder.append({**v, "kind": base.get("kind", "generated"),
                       "selectivity": base.get("selectivity"),
                       "target_z": base.get("target_z"),
                       "off_target_z": base.get("off_target_z")})
out["linkers"] = out["linkers"] + ladder
out["n_length_variants"] = len(ladder)
out["reference_length"] = L.REFERENCE_LENGTH
out["length_ladder"] = list(LENGTH_LADDER)

save_json("02b_linkers.json", out)
print(f"  generated {out['n_generated']}, kept {out['n_kept']}, "
      f"plus {out['n_controls']} controls and {out['n_length_variants']} length "
      f"variants = {len(out['linkers'])} linkers")
for item in out["linkers"]:
    print(f"    {item['id']:<12s} {len(item['sequence']):3d} aa  "
          f"{item['sequence']:<34s} selectivity {item['selectivity']}")
