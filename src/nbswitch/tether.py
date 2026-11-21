"""What a tether is worth.

A mask held near the paratope by a flexible linker behaves as though it were
present at an effective concentration set by the linker's end-to-end
distribution. Modelling a Gly/Ser linker as a worm-like chain,

    <r^2> = 2 lp Lc [1 - (lp/Lc)(1 - exp(-Lc/lp))]
    C_eff(d) = 1660.54 (3 / 2 pi <r^2>)^(3/2) exp(-3 d^2 / 2 <r^2>)   [M]

where d is the distance the linker must span, taken from the structure. Masking
and release then bound the mask-nanobody dissociation constant from opposite
sides, which is what defines a design window.
"""
from __future__ import annotations
import numpy as np

from .config import (PERSISTENCE_LENGTH_A, RISE_PER_RESIDUE_A, CONC_CONVERSION,
                     MAX_SYSTEMIC_LEAK, MIN_TOXIN_RECOVERY, MASK_CONC_M,
                     TOXIN_CONC_M, VHH_TOXIN_KD_M)


def mean_square_end_to_end(n_residues, persistence_A=PERSISTENCE_LENGTH_A,
                           rise_A=RISE_PER_RESIDUE_A):
    """Worm-like-chain <r^2> in square angstroms."""
    contour = np.asarray(n_residues, dtype=float) * rise_A
    return (2 * persistence_A * contour
            * (1 - (persistence_A / contour) * (1 - np.exp(-contour / persistence_A))))


def effective_concentration(n_residues, span_A, **kw):
    """Molar effective concentration of a tethered mask at distance ``span_A``."""
    r2 = mean_square_end_to_end(n_residues, **kw)
    return CONC_CONVERSION * (3 / (2 * np.pi * r2)) ** 1.5 * np.exp(-3 * span_A ** 2 / (2 * r2))


def fraction_masked(c_eff, mask_kd):
    """Fraction of intact construct with the paratope occupied by its own mask."""
    return c_eff / (c_eff + mask_kd)


def fraction_toxin_bound(mask_kd, toxin_kd=VHH_TOXIN_KD_M, toxin=TOXIN_CONC_M,
                         mask=MASK_CONC_M):
    """Fraction of released nanobody bound to toxin, mask competing for the site."""
    return (toxin / toxin_kd) / (1 + toxin / toxin_kd + mask / mask_kd)


def masking_ceiling(c_eff, max_leak=MAX_SYSTEMIC_LEAK):
    """Loosest mask that still keeps systemic leak below ``max_leak``."""
    return c_eff * max_leak / (1 - max_leak)


def release_floor(toxin_kd=VHH_TOXIN_KD_M, toxin=TOXIN_CONC_M, mask=MASK_CONC_M,
                  min_recovery=MIN_TOXIN_RECOVERY):
    """Tightest mask that still lets ``min_recovery`` of the nanobody bind toxin."""
    denominator = (toxin / toxin_kd) * (1 / min_recovery - 1) - 1
    return float("inf") if denominator <= 0 else mask / denominator


def scan_linkers(span_A, lengths, **kw) -> list[dict]:
    """Design window across a ladder of linker lengths for one attachment point."""
    floor = release_floor(**kw)
    rows = []
    for n in lengths:
        c_eff = float(effective_concentration(n, span_A))
        ceiling = masking_ceiling(c_eff)
        open_ = ceiling > floor
        rows.append({"linker_residues": int(n), "span_A": round(float(span_A), 2),
                     "mean_end_to_end_A": round(float(np.sqrt(mean_square_end_to_end(n))), 2),
                     "c_eff_M": c_eff,
                     "km_max_for_masking_M": float(ceiling),
                     "km_min_for_release_M": float(floor),
                     "window_open": bool(open_),
                     "window_decades": float(np.log10(ceiling / floor)) if open_
                                       else float("-inf")})
    return rows


def widest_window(rows):
    """The linker length giving the most room for mask affinity, or None."""
    openings = [r for r in rows if r["window_open"]]
    return max(openings, key=lambda r: r["window_decades"]) if openings else None
