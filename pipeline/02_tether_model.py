#!/usr/bin/env python
"""Stage 2. Turn the measured geometry into bounds on mask affinity."""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from nbswitch import tether as T
from nbswitch.config import (load_json, save_json, MAX_SYSTEMIC_LEAK, MIN_TOXIN_RECOVERY,
                             MASK_CONC_M, TOXIN_CONC_M, VHH_TOXIN_KD_M,
                             PERSISTENCE_LENGTH_A, RISE_PER_RESIDUE_A)

geom = load_json("01_reference.json")["terminal_geometry"]
spans = {"N": geom["n_term_to_paratope_centroid_A"], "C": geom["c_term_to_paratope_centroid_A"]}
lengths = np.arange(2, 61)

scan, best = [], {}
for term, span in spans.items():
    rows = T.scan_linkers(span, lengths)
    for r in rows:
        r["terminus"] = term
    scan += rows
    best[term] = T.widest_window(rows)

sensitivity = [{"KT_M": kt, "T_M": t, "km_min_for_release_M": T.release_floor(toxin_kd=kt, toxin=t)}
               for kt in (1e-10, 1e-9, 1e-8) for t in (1e-8, 1e-7, 1e-6)]

save_json("02_tether.json", {
    "model": {"persistence_length_A": PERSISTENCE_LENGTH_A,
              "rise_per_residue_A": RISE_PER_RESIDUE_A, "leak_max": MAX_SYSTEMIC_LEAK,
              "recovery_min": MIN_TOXIN_RECOVERY, "mask_conc_M": MASK_CONC_M,
              "toxin_conc_M": TOXIN_CONC_M, "KT_M": VHH_TOXIN_KD_M},
    "spans_A": spans, "km_min_for_release_M": T.release_floor(),
    "scan": scan, "best_per_terminus": best, "sensitivity": sensitivity})

print(f"release floor on mask KD: {T.release_floor():.3e} M")
for term, b in best.items():
    closed = [r["linker_residues"] for r in scan
              if r["terminus"] == term and not r["window_open"]]
    note = f"closed at or below {max(closed)} residues" if closed else "open at every length"
    print(f"  {term}-terminal: widest {b['window_decades']:.2f} decades at "
          f"{b['linker_residues']} residues; {note}")
