#!/usr/bin/env python
"""Stage 17. Does the mask bind the nanobody at all, with no tether?

Every masking measure in the pipeline so far has been made on the tethered
construct, where the linker holds the mask against the paratope whether or not
it has any affinity for it. Proximity there is not evidence of binding, so none
of those measures can separate "the mask binds and blocks the site" from "the
mask is simply nearby".

This asks the question on its own. The mask is predicted as a free peptide
against the bare nanobody: no linker, no antigen, nothing holding it in place.
A peptide that cannot find the paratope unaided will not block it when tethered
either, however long the linker.

Two controls decide whether the answer means anything. The antigen itself is a
positive control, since the crystal structure already says it binds the
paratope, and the inert poly-glycine mask is a negative one.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import complex_prediction as CP
from nbswitch.config import load_json, save_json, RESULTS, RANDOM_SEED

BOLTZ = "/scratch/.venv-dtu/bin/boltz"
WORK = Path("/scratch/boltz_jobs")
OUT = RESULTS / "binding"
WORK.mkdir(parents=True, exist_ok=True)

lib = load_json("03_constructs.json")
vhh, toxin = lib["vhh_sequence"], lib["toxin_sequence"]
paratope = [i + 1 for i, l in enumerate(lib["vhh_residue_labels"])
            if l in set(lib["paratope_labels"])]

cases = [("antigen_positive_control", toxin,
          "the antigen itself; the crystal structure already answers this")]
for name, mask in lib["masks"].items():
    cases.append((name, mask["seq"], f"{mask['kind']}: {mask['note']}"))

results = []
for name, peptide, note in cases:
    start = time.time()
    yaml = CP.write_binding_job(WORK / f"bind_{name}.yaml", vhh, peptide)
    try:
        pdb = CP.predict_with_msa(yaml, OUT, binary=BOLTZ, seed=RANDOM_SEED)
    except Exception as exc:                                    # noqa: BLE001
        print(f"  {name}: FAILED {type(exc).__name__}")
        results.append({"case": name, "peptide": peptide, "note": note,
                        "status": "failed", "error": str(exc)[:300]})
        continue
    row = {"case": name, "peptide": peptide, "length": len(peptide), "note": note,
           "status": "ok", "seconds": round(time.time() - start, 1),
           "structure": str(pdb.relative_to(RESULTS))}
    row |= CP.interface_contacts(pdb, "A", "B")
    row |= CP.toxin_paratope_overlap(pdb, "A", "B", paratope)
    # a graft cut from a disulfide-rich antigen can re-form its native bond,
    # pair up cysteines whose partners were left behind, or leave them free
    row["disulfides"] = CP.peptide_disulfides(pdb)
    results.append(row)
    print(f"  {name:22s} {row['seconds']/60:4.1f} min  "
          f"on paratope {row['fraction_of_contacts_on_paratope']:.2f}  "
          f"buried {row['buried_surface_area_A2']:6.0f} A2")

save_json("17_mask_binding.json", {
    "question": "does the mask find the paratope with no linker holding it there",
    "method": "Boltz-2, two chains, alignment for the public nanobody only",
    "n_cases": len(results), "results": results})

ok = [r for r in results if r["status"] == "ok"]
pos = next((r for r in ok if r["case"] == "antigen_positive_control"), None)
if pos:
    print(f"\npositive control (antigen): {pos['fraction_of_contacts_on_paratope']:.0%} "
          f"of contacts on the paratope, {pos['buried_surface_area_A2']:.0f} A2 buried")
    masks = [r for r in ok if r["case"] != "antigen_positive_control"]
    beat = [r for r in masks
            if r["fraction_of_contacts_on_paratope"] >= pos["fraction_of_contacts_on_paratope"]]
    print(f"masks reaching the control's on-paratope fraction: "
          f"{len(beat)} of {len(masks)}"
          + (f" ({', '.join(r['case'] for r in beat)})" if beat else ""))
    inert = next((r for r in masks if r["case"] == "polyG"), None)
    graft = [r for r in masks if "graft" in r["note"]]
    if inert and graft:
        best = max(graft, key=lambda r: r["fraction_of_contacts_on_paratope"])
        print(f"best graft {best['case']} {best['fraction_of_contacts_on_paratope']:.2f} "
              f"against inert control {inert['fraction_of_contacts_on_paratope']:.2f}")
