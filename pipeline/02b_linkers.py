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
N_KEEP = 8
PYTHON = "/scratch/.venv-dtu/bin/python"

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

save_json("02b_linkers.json", out)
print(f"  generated {out['n_generated']}, kept {out['n_kept']}, "
      f"plus {out['n_controls']} controls = {len(out['linkers'])} linkers")
for item in out["linkers"]:
    print(f"    {item['id']:<12s} {item['sequence']:<20s} "
          f"selectivity {item['selectivity']}")
