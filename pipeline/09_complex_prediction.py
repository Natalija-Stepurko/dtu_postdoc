#!/usr/bin/env python
"""Stage 9. Can the toxin still bind when the mask is in place?

ESMFold says whether the mask folds back onto the paratope. Only a two-chain
prediction says whether that matters. The free nanobody with toxin is the
control: the answer is already known from the crystal structure, so it measures
how much to trust the rest.
"""
import sys, time, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import complex_prediction as CP
from nbswitch.constructs import complex_cases
from nbswitch.config import load_json, save_json, RESULTS, RANDOM_SEED

BOLTZ = "/scratch/.venv-dtu/bin/boltz"
WORK = Path("/scratch/boltz_jobs")
OUT = RESULTS / "complexes"
WORK.mkdir(parents=True, exist_ok=True)

lib = load_json("03_constructs.json")
by_id = {c["id"]: c for c in lib["constructs"]}
toxin = lib["toxin_sequence"]

cases = [("control_free_vhh", "VHH_only", "free nanobody with toxin — the crystal "
          "structure already answers this, so it calibrates the rest")]
cases += [(f"masked_{cid}", cid, note)
          for cid, note in complex_cases(lib["constructs"])]
if len(cases) == 1:
    raise SystemExit("no masked construct to predict: the library has no mask+linker pair")

try:
    redesign = load_json("08_redesign_folding.json")
    best = max(redesign["results"], key=lambda r: r["occlusion"])
    rc = {c["id"]: c for c in load_json("07_redesign.json")["constructs"]}[best["id"]]
    by_id[best["id"]] = rc
    cases.append(("masked_best_redesign", best["id"],
                  f"best-occluding ProteinMPNN redesign ({rc['mask_seq']})"))
except (FileNotFoundError, KeyError, ValueError):
    print("  (stage 8 not finished; skipping the redesign case)")

paratope_positions = [i + 1 for i, l in enumerate(lib["vhh_residue_labels"])
                      if l in set(lib["paratope_labels"])]

results = []
for name, cid, note in cases:
    item = by_id.get(cid)
    if item is None:
        print(f"  skip {name}")
        continue
    offset = item["vhh_range"][0] - 1
    positions = [offset + p for p in paratope_positions]
    yaml = CP.write_job(WORK / f"{name}.yaml", {"A": item["sequence"], "B": toxin})
    start = time.time()
    try:
        pdb = CP.predict(yaml, OUT, binary=BOLTZ, seed=RANDOM_SEED)
    except Exception as exc:                                   # noqa: BLE001
        print(f"  {name}: FAILED ({type(exc).__name__})")
        results.append({"case": name, "construct": cid, "note": note,
                        "status": "failed", "error": str(exc)[:300]})
        continue
    elapsed = time.time() - start
    row = {"case": name, "construct": cid, "note": note, "status": "ok",
           "seconds": round(elapsed, 1),
           "construct_length": item["length"], "toxin_length": len(toxin),
           "structure": str(pdb.relative_to(RESULTS))}
    row |= CP.interface_contacts(pdb, "A", "B")
    row |= CP.toxin_paratope_overlap(pdb, "A", "B", positions)
    results.append(row)
    print(f"  {name:22s} {elapsed/60:5.1f} min  contacts {row['atom_contacts']:5d}  "
          f"BSA {row['buried_surface_area_A2']:7.1f}  "
          f"on-paratope {row['fraction_of_contacts_on_paratope']:.2f}")

save_json("09_complexes.json", {
    "method": "Boltz-2, single-sequence (no MSA server; designed sequences stay local)",
    "accelerator": "cpu", "n_cases": len(results), "results": results})
print(f"\n{sum(1 for r in results if r['status']=='ok')}/{len(results)} cases predicted")
