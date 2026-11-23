#!/usr/bin/env python
"""Stage 19b. The two structures the final screen comes down to, in one frame.

Stage 19 predicts the shortlisted constructs with their antigen and reports
numbers. The question those numbers answer is visual: does the mask sit where
the antigen binds, and is the antigen pushed off that surface? Two pictures
answer it directly, but only if they are in the same orientation - two
separately predicted complexes come out of Boltz in arbitrary frames, and a
reader comparing them would be comparing viewpoints.

So the construct's model is superposed onto the control's by the nanobody
alpha carbons alone. The nanobody is the part both structures share and the
part neither prediction was asked to change, which makes it the right frame:
anything that moves between the two pictures moves relative to the nanobody,
which is what the comparison is about.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from biotite.structure import superimpose
from biotite.structure.info import one_letter_code
from biotite.structure.io.pdb import PDBFile

from nbswitch.config import RESULTS, STRUCTURES, load_json, save_json

#: Where the viewer expects them.
REFERENCE_OUT = STRUCTURES / "final_reference.pdb"
MASKED_OUT = STRUCTURES / "final_masked.pdb"

#: Which chain carries the nanobody in a Boltz prediction of these cases. It
#: is the first chain in both, because that is the order the inputs are given.
VHH_CHAIN = "A"

ver = load_json("19_verification.json")
rows = ver["results"]

control = next((r for r in rows if r["construct"] is None and r["status"] == "ok"), None)
picked = next((r for r in rows if r.get("verdict") == "worth testing"), None)
if control is None or picked is None:
    raise SystemExit("stage 19 has no control and no construct worth testing; "
                     "nothing to pair")


def _ca(atoms, chain):
    sel = atoms[(atoms.chain_id == chain) & (atoms.atom_name == "CA")]
    return "".join(one_letter_code(r) or "X" for r in sel.res_name), sel


def _read(rel):
    path = RESULTS / rel
    if not path.exists():
        raise SystemExit(f"missing predicted structure: {path}")
    return PDBFile.read(str(path)).get_structure(model=1)


ref = _read(control["structure"])
msk = _read(picked["structure"])

seq_ref, ca_ref = _ca(ref, VHH_CHAIN)
seq_msk, ca_msk = _ca(msk, VHH_CHAIN)
n = min(len(seq_ref), len(seq_msk))
if seq_ref[:n] != seq_msk[:n]:
    raise SystemExit("the two models do not share a nanobody sequence; the "
                     "superposition would be meaningless")

fitted, transform = superimpose(ca_ref[:n], ca_msk[:n])
rmsd = float(np.sqrt(np.mean(np.sum(
    (fitted.coord - ca_ref[:n].coord) ** 2, axis=1))))

STRUCTURES.mkdir(parents=True, exist_ok=True)
out = PDBFile()
out.set_structure(ref)
out.write(str(REFERENCE_OUT))
out = PDBFile()
out.set_structure(transform.apply(msk))
out.write(str(MASKED_OUT))


def _chains(atoms):
    return {str(c): int(len({int(i) for i in atoms[atoms.chain_id == c].res_id}))
            for c in sorted(set(atoms.chain_id))}


save_json("19b_final_pair.json", {
    "note": "the construct's prediction superposed onto the control's by the "
            "nanobody alpha carbons, so the two can be shown in one frame",
    "superposed_on": f"chain {VHH_CHAIN}, {n} alpha carbons",
    "vhh_ca_rmsd_A": round(rmsd, 2),
    "reference": {"case": control["case"], "file": REFERENCE_OUT.name,
                  "chains": _chains(ref),
                  "fraction_of_contacts_on_paratope":
                      control["fraction_of_contacts_on_paratope"]},
    "masked": {"case": picked["case"], "construct": picked["construct"],
               "file": MASKED_OUT.name, "chains": _chains(msk),
               "fraction_of_contacts_on_paratope":
                   picked["fraction_of_contacts_on_paratope"],
               "mask_fraction_on_paratope": picked["mask_fraction_on_paratope"],
               "antigen_contacts_displaced": picked["antigen_contacts_displaced"]},
})

print(f"superposed on {n} nanobody alpha carbons, RMSD {rmsd:.2f} A")
print(f"  reference {control['case']}: {_chains(ref)}")
print(f"  masked    {picked['case']}: {_chains(msk)}")
print(f"  wrote {REFERENCE_OUT.name} and {MASKED_OUT.name}")
