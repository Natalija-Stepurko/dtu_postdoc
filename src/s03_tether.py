"""Stage 2. Tethered-mask thermodynamics from the measured 9GCN geometry.

A tethered mask is held near the paratope at an effective concentration C_eff set
by the linker. Modelling the Gly/Ser linker as a worm-like chain (Zhou 2001,
J Phys Chem B; persistence length for flexible GS linkers ~4.5 A, contour length
3.8 A/residue), the end-to-end distribution gives

    <r^2> = 2 lp Lc [1 - (lp/Lc)(1 - exp(-Lc/lp))]
    C_eff(d) = 1660.54 * (3/(2 pi <r^2>))^{3/2} * exp(-3 d^2 / (2 <r^2>))   [M]

with 1660.54 converting molecules/A^3 to mol/L. d is the distance the linker must
span, taken from the structure. Masking and release are then

    P_masked   = C_eff / (C_eff + KM)                     (before cleavage)
    P_NT       = (T/KT) / (1 + T/KT + M/KM)               (after cleavage)

Two requirements pull KM in opposite directions, which defines a design window.
"""
import numpy as np
from common import save_json, load_json, style, TOL, FIGS

LP = 4.5       # persistence length, A
RISE = 3.8     # contour length per residue, A
CONV = 1660.54 # molecules/A^3 -> M

def mean_sq_end_to_end(n_res, lp=LP, rise=RISE):
    Lc = n_res * rise
    return 2 * lp * Lc * (1 - (lp / Lc) * (1 - np.exp(-Lc / lp)))

def c_eff(n_res, d):
    r2 = mean_sq_end_to_end(n_res)
    return CONV * (3 / (2 * np.pi * r2)) ** 1.5 * np.exp(-3 * d ** 2 / (2 * r2))

geom = load_json("s02_interface.json")["terminal_geometry"]
D_N = geom["n_term_to_paratope_centroid_A"]
D_C = geom["c_term_to_paratope_centroid_A"]

# Requirements
LEAK_MAX = 0.01   # <=1% of construct unmasked in circulation
RECOVERY = 0.90   # >=90% of cleaved nanobody toxin-bound
M_FREE = 1e-6     # released mask concentration, M (= construct concentration)
KT = 1e-9         # VHH-toxin KD, M (swept separately)
T_FREE = 1e-7     # free toxin concentration, M

def km_upper(ce, leak=LEAK_MAX):
    """Masking requires P_masked >= 1-leak  ->  KM <= C_eff * leak/(1-leak)."""
    return ce * leak / (1 - leak)

def km_lower(kt=KT, t=T_FREE, m=M_FREE, rec=RECOVERY):
    """Release requires P_NT >= rec  ->  KM >= m / ((t/kt)(1/rec - 1) - 1)."""
    denom = (t / kt) * (1 / rec - 1) - 1
    return np.inf if denom <= 0 else m / denom

ns = np.arange(2, 61)
rows = []
for term, d in (("N", D_N), ("C", D_C)):
    for n in ns:
        ce = c_eff(n, d)
        hi = km_upper(ce)
        lo = km_lower()
        rows.append({"terminus": term, "linker_residues": int(n),
                     "span_A": round(d, 2),
                     "mean_end_to_end_A": round(float(np.sqrt(mean_sq_end_to_end(n))), 2),
                     "c_eff_M": float(ce),
                     "km_max_for_masking_M": float(hi),
                     "km_min_for_release_M": float(lo),
                     "window_open": bool(hi > lo),
                     "window_decades": float(np.log10(hi / lo)) if hi > lo else float("-inf")})

best = {}
for term in ("N", "C"):
    sub = [r for r in rows if r["terminus"] == term and r["window_open"]]
    best[term] = max(sub, key=lambda r: r["window_decades"]) if sub else None

# sensitivity of the lower bound to toxin affinity and concentration
sens = []
for kt in (1e-10, 1e-9, 1e-8):
    for t in (1e-8, 1e-7, 1e-6):
        sens.append({"KT_M": kt, "T_M": t, "km_min_for_release_M": float(km_lower(kt=kt, t=t))})

out = {"model": {"persistence_length_A": LP, "rise_per_residue_A": RISE,
                 "leak_max": LEAK_MAX, "recovery_min": RECOVERY,
                 "mask_conc_M": M_FREE, "toxin_conc_M": T_FREE, "KT_M": KT},
       "span_N_term_A": D_N, "span_C_term_A": D_C,
       "km_min_for_release_M": float(km_lower()),
       "scan": rows, "best_per_terminus": best, "sensitivity": sens}
save_json("s03_tether.json", out)

print(f"span to paratope centroid: N-term {D_N} A, C-term {D_C} A")
print(f"release floor on KM (KT={KT:.0e} M, T={T_FREE:.0e} M, M={M_FREE:.0e} M): "
      f"{km_lower():.3e} M")
print("\n  term  n_res  <r> (A)     C_eff (M)   KM_max(mask)   window (decades)")
for term in ("N", "C"):
    for n in (4, 8, 12, 16, 20, 30, 40, 60):
        r = next(x for x in rows if x["terminus"] == term and x["linker_residues"] == n)
        w = f"{r['window_decades']:.2f}" if r["window_open"] else "CLOSED"
        print(f"  {term}   {n:5d}  {r['mean_end_to_end_A']:7.1f}  {r['c_eff_M']:12.3e}  "
              f"{r['km_max_for_masking_M']:12.3e}   {w}")
print()
for term, b in best.items():
    if b: print(f"  widest window, {term}-terminal fusion: {b['linker_residues']} residues, "
                f"{b['window_decades']:.2f} decades of KM")
    else: print(f"  {term}-terminal fusion: no open window at any linker length")
