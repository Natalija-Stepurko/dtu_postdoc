#!/usr/bin/env python
"""Stage 14. One file holding every construct aligned on the nanobody."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import ensemble
from nbswitch.config import load_json, save_json, STRUCTURES, DOCS

lib = load_json("03_constructs.json")
constructs = list(lib["constructs"])
try:
    constructs += load_json("07_redesign.json")["constructs"]
except FileNotFoundError:
    pass

out = ensemble.superimpose_on_nanobody(
    constructs, STRUCTURES, "VHH_only", DOCS / "structures" / "ensemble.pdb")
save_json("14_ensemble.json", {"n_models": out["n_models"],
                               "constructs": out["constructs"],
                               "file": "structures/ensemble.pdb"})
print(f"{out['n_models']} models written ({out['n_models']-1} masks on one nanobody)")
