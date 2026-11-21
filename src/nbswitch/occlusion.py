"""Measuring whether a mask covers the paratope.

The quantity reported is paratope occlusion: how much of the paratope's
solvent-accessible surface in the folded free nanobody is lost in the folded
fusion. It is a property of a predicted structure, so it is a hypothesis about
the masked state and not a masking efficiency or a dissociation constant.
"""
from __future__ import annotations
import io
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb

from .structure import residue_sasa, residues


def read_pdb(path: Path):
    atoms = pdb.PDBFile.read(str(path)).get_structure(model=1)
    return atoms[struc.filter_amino_acids(atoms)]


def paratope_indices(vhh_residue_labels, paratope_labels) -> np.ndarray:
    wanted = set(paratope_labels)
    return np.array([i for i, label in enumerate(vhh_residue_labels) if label in wanted])


def occlusion(construct, structure_dir: Path, para_idx, free_paratope_sasa, plddt):
    """Occlusion and contact geometry for one folded construct."""
    atoms = read_pdb(structure_dir / f"{construct['id']}.pdb")
    offset = construct["vhh_range"][0] - 1
    sasa = residue_sasa(atoms)[offset + para_idx]
    scores = plddt[construct["id"]]
    row = {"id": construct["id"], "mask": construct["mask"],
           "kind": construct["mask_kind"], "linker": construct["linker_residues"],
           "paratope_sasa_A2": round(float(sasa.sum()), 1),
           "occlusion": round(float(1 - sasa.sum() / free_paratope_sasa), 4),
           "plddt_all": round(float(scores.mean()), 2),
           "plddt_vhh": round(float(scores[offset:].mean()), 2)}
    if construct["mask_range"]:
        ca = atoms[atoms.atom_name == "CA"]
        mask_idx = np.arange(construct["mask_range"][0] - 1, construct["mask_range"][1])
        d = np.linalg.norm(ca.coord[mask_idx][:, None]
                           - ca.coord[offset + para_idx][None], axis=-1)
        row["min_mask_paratope_CA_A"] = round(float(d.min()), 2)
        row["n_mask_paratope_contacts_8A"] = int((d < 8.0).sum())
        row["plddt_mask"] = round(float(scores[mask_idx].mean()), 2)
    return row


def vhh_core_rmsd(construct, structure_dir: Path, reference_id="VHH_only") -> float:
    """CA RMSD of the nanobody inside the fusion against the folded free nanobody."""
    fusion = read_pdb(structure_dir / f"{construct['id']}.pdb")
    free = read_pdb(structure_dir / f"{reference_id}.pdb")
    offset = construct["vhh_range"][0] - 1
    a = fusion[fusion.atom_name == "CA"].coord[offset:]
    b = free[free.atom_name == "CA"].coord
    n = min(len(a), len(b))
    fitted, _ = struc.superimpose(b[:n], a[:n])
    return round(float(np.sqrt(((fitted - b[:n]) ** 2).sum(axis=1).mean())), 3)
