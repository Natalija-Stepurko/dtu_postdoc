"""Score 01: can this protein be made at all.

Nothing downstream matters if the construct will not express, folds badly, or
aggregates, so this is scored first and independently of what the mask is for.
Everything here is computed from the sequence and the predicted structure.

The components are deliberately simple and each is a known liability rather
than a learned quantity: confidence that the fold was predicted at all, whether
adding the mask disturbed the nanobody, how much hydrophobic surface is left
exposed, and how much cysteine the mask adds beyond the nanobody's own pair.

On that last point: an earlier version counted every cysteine past two as
unpaired, which sent any cysteine-rich graft to a feasibility of zero. That is
an assumption, not a measurement, and here it is the wrong one. The graft this
pipeline takes from 9RIU is a three-finger-toxin loop whose cysteines are
disulfide-bonded in the antigen, two of them to partners outside the excised
segment. Added cysteine is a real risk to folding and to the nanobody's own
disulfide, but whether a given one is unpaired is something the predicted
structure can be asked rather than assumed, so it is measured here and the
score degrades with the count instead of collapsing.
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

#: The VHH's own conserved intradomain disulfide.
NATIVE_VHH_CYSTEINES = 2

#: Longest S-S distance counted as a disulfide in a predicted model. Real bonds
#: are near 2.05 A; the allowance is for prediction error, not for contacts.
DISULFIDE_MAX_A = 2.5

#: Penalty per unpaired cysteine, and the floor it cannot fall below. A free
#: thiol is a genuine liability for expression and for disulfide shuffling, but
#: it is one liability among several and should not by itself zero a product.
CYS_PENALTY = 0.15
CYS_FLOOR = 0.25


def unpaired_cysteines(atoms, max_distance: float = DISULFIDE_MAX_A) -> int:
    """Cysteines in the predicted model whose SG has no partner within bonding range."""
    sg = atoms[(atoms.atom_name == "SG") & (atoms.res_name == "CYS")]
    n = sg.array_length()
    if n == 0:
        return 0
    d = np.linalg.norm(sg.coord[:, None, :] - sg.coord[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return int((d.min(axis=1) > max_distance).sum())


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
    cys_total = sequence.count("C")
    added_cys = max(0, cys_total - NATIVE_VHH_CYSTEINES)
    unpaired_cys = unpaired_cysteines(atoms)

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
    # a graded penalty on cysteine the construct adds, measured against the
    # predicted pairing rather than assumed from the count
    f_cys = float(np.clip(1.0 - CYS_PENALTY * unpaired_cys, CYS_FLOOR, 1.0))

    return {
        "id": construct["id"], "mask": construct["mask"], "kind": construct["mask_kind"],
        "linker_residues": construct["linker_residues"],
        "plddt_all": round(plddt_all, 1), "plddt_vhh": round(plddt_vhh, 1),
        "vhh_core_rmsd_A": round(rmsd, 3) if rmsd == rmsd else None,
        "exposed_hydrophobic_fraction": round(phobic_fraction, 3),
        "cysteines_total": cys_total,
        "added_cysteines": added_cys,
        "unpaired_cysteines": unpaired_cys,
        "f_fold": round(f_fold, 3), "f_preserved": round(f_preserved, 3),
        "f_surface": round(f_surface, 3), "f_cysteine": round(f_cys, 3),
        "feasibility_score": round(float(f_fold * f_preserved * f_surface * f_cys), 3),
    }


def score_library(constructs, structure_dir: Path, plddt: dict, log=print) -> list[dict]:
    rows = [r for r in (feasibility(c, structure_dir, plddt) for c in constructs) if r]
    log(f"scored {len(rows)} constructs for feasibility")
    return rows
