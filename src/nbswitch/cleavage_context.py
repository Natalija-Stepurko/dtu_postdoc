"""Can a protease reach and engage the linker once the construct is folded?

A sequence model scores a peptide floating free. In the construct the same
peptide is tethered at both ends, between a mask and a nanobody, and may be
packed against either. A substrate a protease recognises perfectly is not
cleaved if the enzyme cannot physically reach it, so sequence recognition and
accessibility are separate questions and are kept as separate scores.

Three requirements are checked, each taken from how a metalloproteinase engages
a substrate rather than from a fitted model:

**Extension.** The catalytic cleft binds substrate in an extended conformation
across roughly P3 to P3'. A segment folded back on itself, or wound into a
turn, cannot lie in the groove. Measured as the span of the six residues
around the scissile bond against the length they would have fully extended.

**Exposure.** Buried backbone cannot be engaged. Measured as solvent-accessible
surface over the substrate, relative to the same residues in isolation.

**Clearance.** This is the requirement a sequence model cannot express and the
reason this module exists. The catalytic domain is a globular protein roughly
forty angstroms across, not a point. It has to occupy space next to the
scissile bond. Clearance asks how much of the volume a catalytic domain would
need is already filled by the construct itself: a linker lying in a groove
between the mask and the nanobody can be extended and solvent-exposed by the
usual measures and still be unreachable by anything larger than a water
molecule.

None of this is a prediction of rate. It is a set of necessary conditions, and
a construct failing any of them will not be cleaved whatever its sequence
score says.
"""
from __future__ import annotations

import numpy as np
import biotite.structure as struc

from .structure import residues, residue_sasa

#: Rise per residue in an extended backbone.
EXTENDED_RISE_A = 3.3

#: Radius of a sphere standing in for an MMP catalytic domain. Catalytic
#: domains are about 160 residues; a globular protein of that size has a radius
#: of gyration near 15 A and an excluded radius a little larger.
CATALYTIC_DOMAIN_RADIUS_A = 17.0

#: How far from the scissile bond the probe sphere is centred: touching the
#: substrate, not overlapping it.
PROBE_OFFSET_A = 10.0

#: Directions sampled when looking for an approach path.
N_PROBE_DIRECTIONS = 64


def _fibonacci_directions(n: int) -> np.ndarray:
    """Roughly even directions on a sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


#: Residues either side of the scissile bond that the catalytic cleft holds,
#: about P3 to P3'. Extension over this window is what decides whether the
#: groove can bind; extension over the whole linker is a different quantity
#: and a long flexible linker will look coiled whatever its cleavage site does.
CLEFT_HALF_WIDTH = 3


def extension(atoms, linker_range: tuple[int, int]) -> float:
    """Span of a segment against its fully extended length, 0 to 1."""
    start, stop = linker_range
    ca = atoms[(atoms.atom_name == "CA") &
               (atoms.res_id >= start) & (atoms.res_id <= stop)]
    n = ca.array_length()
    if n < 2:
        return 0.0
    span = float(np.linalg.norm(ca.coord[-1] - ca.coord[0]))
    return min(span / ((n - 1) * EXTENDED_RISE_A), 1.0)


def exposure(atoms, linker_range: tuple[int, int]) -> float:
    """Fraction of the linker's surface that is solvent accessible."""
    start, stop = linker_range
    res = residues(atoms)
    sasa_in_construct = residue_sasa(atoms)
    mask = np.array([start <= r <= stop for r in res.res_id])
    if not mask.any():
        return 0.0
    isolated = atoms[(atoms.res_id >= start) & (atoms.res_id <= stop)]
    sasa_alone = residue_sasa(isolated)
    total_alone = float(sasa_alone.sum())
    if total_alone <= 0:
        return 0.0
    return float(min(sasa_in_construct[mask].sum() / total_alone, 1.0))


def clearance(atoms, scissile_resid: int,
              radius: float = CATALYTIC_DOMAIN_RADIUS_A,
              offset: float = PROBE_OFFSET_A,
              n_directions: int = N_PROBE_DIRECTIONS) -> dict:
    """Fraction of approach directions with room for a catalytic domain.

    A sphere the size of a catalytic domain is placed just off the scissile
    bond in many directions. A direction counts as open when no atom of the
    construct lies inside that sphere. The best direction is reported too,
    since one clear approach is enough for an enzyme.
    """
    ca = atoms[(atoms.atom_name == "CA") & (atoms.res_id == scissile_resid)]
    if ca.array_length() == 0:
        return {"open_fraction": 0.0, "min_occupancy": 1.0}

    centre = ca.coord[0]
    heavy = atoms[atoms.element != "H"]
    # the linker's own atoms do not block the enzyme that is binding them
    body = heavy[np.abs(heavy.res_id - scissile_resid) > 3]
    if body.array_length() == 0:
        return {"open_fraction": 1.0, "min_occupancy": 0.0}

    directions = _fibonacci_directions(n_directions)
    probes = centre + directions * (offset + radius)
    d = np.linalg.norm(body.coord[None, :, :] - probes[:, None, :], axis=-1)
    inside = (d < radius).sum(axis=1)
    per_direction = inside / max(body.array_length(), 1)
    return {"open_fraction": float((inside == 0).mean()),
            "min_occupancy": float(per_direction.min())}


def cleavage_accessibility(atoms, linker_range: tuple[int, int]) -> dict:
    """The three necessary conditions, and their product.

    The product is deliberate: these are requirements, not contributions. A
    construct that fails any one of them is not cleavable, so a high score
    elsewhere must not compensate.

    Extension is measured over the cleft window around the scissile bond, not
    over the whole linker. A protease holds about six residues; whether the
    other twenty are coiled is a fact about the tether and not about whether
    this bond can be bound. The whole-linker value is reported alongside,
    because it is what decides how far the mask reaches.
    """
    start, stop = linker_range
    scissile = (start + stop) // 2
    cleft = (max(start, scissile - CLEFT_HALF_WIDTH),
             min(stop, scissile + CLEFT_HALF_WIDTH))
    ext_local = extension(atoms, cleft)
    ext_whole = extension(atoms, linker_range)
    exp = exposure(atoms, cleft)
    exp_whole = exposure(atoms, linker_range)
    clear = clearance(atoms, scissile)
    score = ext_local * exp * clear["open_fraction"]
    return {
        "extension": round(ext_local, 3),
        "extension_whole_linker": round(ext_whole, 3),
        "exposure": round(exp, 3),
        "exposure_whole_linker": round(exp_whole, 3),
        "clearance_open_fraction": round(clear["open_fraction"], 3),
        "clearance_min_occupancy": round(clear["min_occupancy"], 4),
        "scissile_residue": scissile,
        "cleft_window": list(cleft),
        "context_score": round(score, 4),
    }
