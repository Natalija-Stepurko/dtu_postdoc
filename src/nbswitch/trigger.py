"""Scoring whether a linker could be cut at all.

Cleavage needs two independent things: a sequence the protease recognises, and
a linker the protease can physically reach. The first needs specificity data for
the enzymes in question. The second does not - it is a property of the folded
construct, and a linker packed against the nanobody surface cannot be cut
whatever its sequence, because proteases bind their substrate threaded through
the active-site cleft in an extended conformation.

This module scores the second. It is a necessary condition, not a sufficient
one, and the page says so.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import biotite.structure as struc

from .structure import residue_sasa, residues
from .occlusion import read_pdb

#: A residue is counted as exposed above this side-chain-inclusive SASA.
EXPOSED_A2 = 40.0
#: Rise per residue for a fully extended backbone, as used by the tether model.
RISE_PER_RESIDUE_A = 3.8


def linker_accessibility(construct: dict, structure_dir: Path) -> dict:
    """How reachable and how extended the linker is in the folded construct."""
    if not construct.get("linker_range"):
        return {}
    atoms = read_pdb(structure_dir / f"{construct['id']}.pdb")
    sasa = residue_sasa(atoms)
    start, stop = construct["linker_range"]            # 1-based inclusive
    idx = np.arange(start - 1, stop)
    if len(idx) == 0:
        return {}
    linker_sasa = sasa[idx]

    ca = atoms[atoms.atom_name == "CA"]
    coords = ca.coord[idx]
    span = float(np.linalg.norm(coords[-1] - coords[0])) if len(coords) > 1 else 0.0
    contour = (len(idx) - 1) * RISE_PER_RESIDUE_A
    extension = span / contour if contour else 0.0

    # how close the linker sits to the rest of the construct
    other = np.setdiff1d(np.arange(ca.array_length()), idx)
    if len(other):
        d = np.linalg.norm(coords[:, None, :] - ca.coord[other][None, :, :], axis=-1)
        contacts = int((d < 8.0).sum())
        nearest = float(d.min())
    else:
        contacts, nearest = 0, float("nan")

    exposed = float((linker_sasa > EXPOSED_A2).mean())
    return {
        "linker_residues": len(idx),
        "linker_mean_sasa_A2": round(float(linker_sasa.mean()), 1),
        "linker_min_sasa_A2": round(float(linker_sasa.min()), 1),
        "fraction_exposed": round(exposed, 3),
        "extension_ratio": round(extension, 3),
        "contacts_with_rest_8A": contacts,
        "nearest_approach_A": round(nearest, 2),
        "accessibility_score": round(float(exposed * min(extension / 0.7, 1.0)), 3),
    }


def score_library(constructs, structure_dir: Path, log=print) -> list[dict]:
    """Accessibility for every construct that has a linker and a structure."""
    rows = []
    for item in constructs:
        path = structure_dir / f"{item['id']}.pdb"
        if not path.exists() or not item.get("linker_range"):
            continue
        row = {"id": item["id"], "mask": item["mask"], "kind": item["mask_kind"]}
        row |= linker_accessibility(item, structure_dir)
        if row.get("linker_residues"):
            rows.append(row)
    log(f"scored {len(rows)} linkers")
    return rows
