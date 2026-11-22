#!/usr/bin/env python
"""Stage 2. Turn the measured geometry into bounds on mask affinity."""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from nbswitch import tether as T
from nbswitch.config import (load_json, save_json, MAX_SYSTEMIC_LEAK, MIN_TOXIN_RECOVERY,
                             MASK_CONC_M, TOXIN_CONC_M, VHH_TOXIN_KD_M,
                             VHH_TOXIN_KD_RANGE_M, VHH_TOXIN_KOFF_MAX_S,
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

# The nanobody's own affinity is reported as a range, not measured for this
# pair, so the window is reported across that range rather than at one point.
lo, hi = VHH_TOXIN_KD_RANGE_M
# the grid carries the value the rest of the page uses and the point where
# the window shuts, so the figure and the text cannot disagree about either
best_n = best["N"]["linker_residues"]
kd_closes = T.affinity_at_which_window_closes(spans["N"], best_n)
kd_grid = sorted(set(list(np.logspace(np.log10(lo), np.log10(hi), 9))
                     + [VHH_TOXIN_KD_M, kd_closes * 0.995]))
affinity_scan = T.window_across_affinity(spans["N"], best_n, kd_grid)

save_json("02_tether.json", {
    "model": {"persistence_length_A": PERSISTENCE_LENGTH_A,
              "rise_per_residue_A": RISE_PER_RESIDUE_A, "leak_max": MAX_SYSTEMIC_LEAK,
              "recovery_min": MIN_TOXIN_RECOVERY, "mask_conc_M": MASK_CONC_M,
              "toxin_conc_M": TOXIN_CONC_M, "KT_M": VHH_TOXIN_KD_M},
    "spans_A": spans, "km_min_for_release_M": T.release_floor(),
    "scan": scan, "best_per_terminus": best, "sensitivity": sensitivity,
    "affinity": {"reported_range_M": list(VHH_TOXIN_KD_RANGE_M),
                 "koff_max_s": VHH_TOXIN_KOFF_MAX_S,
                 "used_M": VHH_TOXIN_KD_M,
                 "linker_residues": best_n,
                 "window_closes_above_M": kd_closes,
                 "scan": affinity_scan}})

print(f"release floor on mask KD: {T.release_floor():.3e} M")
print(f"  window closes once the nanobody binds its antigen more weakly than "
      f"{kd_closes*1e9:.1f} nM")
for r in affinity_scan:
    d = r["window_decades"]
    print(f"    KT {r['vhh_toxin_kd_M']*1e9:7.2f} nM -> floor "
          f"{r['km_min_for_release_M']*1e9:9.1f} nM, "
          + (f"{d:.2f} decades" if np.isfinite(d) else "CLOSED"))
for term, b in best.items():
    closed = [r["linker_residues"] for r in scan
              if r["terminus"] == term and not r["window_open"]]
    note = f"closed at or below {max(closed)} residues" if closed else "open at every length"
    print(f"  {term}-terminal: widest {b['window_decades']:.2f} decades at "
          f"{b['linker_residues']} residues; {note}")
