"""Score 01: can this protein be made at all.

Nothing downstream matters if the construct will not express, folds badly, or
aggregates, so this is scored first and independently of what the mask is for.
Everything here is computed from the sequence and the predicted structure.

The components are deliberately simple and each is a known liability rather
than a learned quantity: confidence that the fold was predicted at all, whether
adding the mask disturbed the nanobody, how much hydrophobic surface is left
exposed, and whether the construct carries unpaired cysteines.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import biotite.structure as struc

from .structure import residue_sasa, residues
from .occlusion import read_pdb, vhh_core_rmsd

#: Kyte-Doolittle hydrophobic set, used for exposed-hydrophobic surface.
HYDROPHOBIC = set("AILMFWVY")
#: Above this fraction of exposed surface being hydrophobic, aggregation risk rises.
HYDROPHOBIC_REFERENCE = 0.35


def feasibility(construct: dict, structure_dir: Path, plddt: dict,
                reference_id: str = "VHH_only") -> dict:
    """Per-construct folding confidence, fold preservation and surface liabilities."""
    path = structure_dir / f"{construct['id']}.pdb"
    if not path.exists():
        return {}
    atoms = read_pdb(path)
    res = residues(atoms)
    sasa = residue_sasa(atoms)
    scores = plddt.get(construct["id"])
    if scores is None:
        return {}

    offset = construct["vhh_range"][0] - 1
    hydrophobic = np.array([c in HYDROPHOBIC for c in res.code])
    exposed_total = float(sasa.sum())
    exposed_phobic = float(sasa[hydrophobic].sum())
    phobic_fraction = exposed_phobic / exposed_total if exposed_total else 0.0

    sequence = construct["sequence"]
    cys = sequence.count("C")
    # the nanobody's own pair is expected; anything beyond it is unpaired
    unpaired_cys = max(0, cys - 2)

    try:
        rmsd = vhh_core_rmsd(construct, structure_dir, reference_id)
    except (FileNotFoundError, ValueError):
        rmsd = float("nan")

    plddt_all = float(scores.mean())
    plddt_vhh = float(scores[offset:].mean())

    # each term is 0 to 1, higher being better, and they are kept visible
    f_fold = min(plddt_all / 90.0, 1.0)
    f_preserved = float(np.clip(1.0 - (rmsd / 3.0), 0, 1)) if rmsd == rmsd else 0.0
    f_surface = float(np.clip(1.0 - (phobic_fraction / HYDROPHOBIC_REFERENCE), 0, 1))
    f_cys = 1.0 if unpaired_cys == 0 else max(0.0, 1.0 - 0.25 * unpaired_cys)

    return {
        "id": construct["id"], "mask": construct["mask"], "kind": construct["mask_kind"],
        "linker_residues": construct["linker_residues"],
        "plddt_all": round(plddt_all, 1), "plddt_vhh": round(plddt_vhh, 1),
        "vhh_core_rmsd_A": round(rmsd, 3) if rmsd == rmsd else None,
        "exposed_hydrophobic_fraction": round(phobic_fraction, 3),
        "unpaired_cysteines": unpaired_cys,
        "f_fold": round(f_fold, 3), "f_preserved": round(f_preserved, 3),
        "f_surface": round(f_surface, 3), "f_cysteine": round(f_cys, 3),
        "feasibility_score": round(float(f_fold * f_preserved * f_surface * f_cys), 3),
    }


def score_library(constructs, structure_dir: Path, plddt: dict, log=print) -> list[dict]:
    rows = [r for r in (feasibility(c, structure_dir, plddt) for c in constructs) if r]
    log(f"scored {len(rows)} constructs for feasibility")
    return rows
