#!/usr/bin/env python
"""Stage 11. The same complexes, with alignments for the public chains.

Stage 9 ran every chain without an alignment, because an alignment over the
whole fusion would send the designed mask to a remote server. That turned out
to be why it failed. Splitting the construct into three chains joined by an
explicit bond lets the nanobody and the antigen take real alignments, since
both sequences are already public in the PDB, while the designed peptide takes
none and stays on this machine.
"""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import complex_prediction as CP
from nbswitch.config import load_json, save_json, RESULTS, RANDOM_SEED

BOLTZ = "/scratch/.venv-dtu/bin/boltz"
WORK = Path("/scratch/boltz_msa_stage")
OUT = RESULTS / "complexes_msa"
WORK.mkdir(parents=True, exist_ok=True)

lib = load_json("03_constructs.json")
vhh, toxin = lib["vhh_sequence"], lib["toxin_sequence"]
by_id = {c["id"]: c for c in lib["constructs"]}
paratope = [i + 1 for i, l in enumerate(lib["vhh_residue_labels"])
            if l in set(lib["paratope_labels"])]

cases = [("control_free_vhh_msa", None,
          "free nanobody with its toxin; the crystal structure already answers this")]
for cid in ("epi1_L5", "epi1_L15", "polyG_L5"):
    if cid in by_id:
        cases.append((f"{cid}_msa", cid, f"{by_id[cid]['mask_kind']}, "
                      f"{by_id[cid]['linker_residues']}-residue linker"))

results = []
for name, cid, note in cases:
    start = time.time()
    if cid is None:
        yaml = CP.write_job(WORK / f"{name}.yaml", {"A": vhh, "B": toxin},
                            )
        # the two-chain job wants alignments too, so rewrite without msa: empty
        yaml.write_text("version: 1\nsequences:\n"
                        f"  - protein:\n      id: A\n      sequence: {vhh}\n"
                        f"  - protein:\n      id: B\n      sequence: {toxin}\n")
        chains = ("A", "B")
        offset = 0
    else:
        item = by_id[cid]
        yaml = CP.write_tethered_job(WORK / f"{name}.yaml", vhh,
                                     item["mask_seq"] + item["linker_seq"], toxin)
        chains = ("A", "C")
        offset = 0
    try:
        pdb = CP.predict_with_msa(yaml, OUT, binary=BOLTZ, seed=RANDOM_SEED)
    except Exception as exc:                                    # noqa: BLE001
        print(f"  {name}: FAILED {type(exc).__name__}")
        results.append({"case": name, "construct": cid, "note": note,
                        "status": "failed", "error": str(exc)[:300]})
        continue
    row = {"case": name, "construct": cid, "note": note, "status": "ok",
           "seconds": round(time.time() - start, 1),
           "structure": str(pdb.relative_to(RESULTS))}
    row |= CP.interface_contacts(pdb, chains[0], chains[1])
    row |= CP.toxin_paratope_overlap(pdb, chains[0], chains[1],
                                     [offset + p for p in paratope])
    if cid is not None:
        mask = CP.toxin_paratope_overlap(pdb, "A", "B", paratope)
        row["mask_paratope_contacts"] = mask["paratope_contacts"]
        row["mask_fraction_on_paratope"] = mask["fraction_of_contacts_on_paratope"]
    results.append(row)
    print(f"  {name:24s} {row['seconds']/60:4.1f} min  "
          f"toxin on paratope {row['fraction_of_contacts_on_paratope']:.2f}"
          + (f"   mask on paratope {row['mask_fraction_on_paratope']:.2f}"
             if cid else ""))

save_json("11_complexes_msa.json", {
    "method": "Boltz-2, three chains joined by a bond; alignments for the public "
              "nanobody and antigen chains only",
    "n_cases": len(results), "results": results})
print(f"\n{sum(1 for r in results if r['status']=='ok')}/{len(results)} predicted")
