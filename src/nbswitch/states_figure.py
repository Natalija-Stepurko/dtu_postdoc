"""Render the three states from real coordinates, in one orientation.

The schematic version drew the nanobody as a rounded box. These panels use the
deposited complex and the predicted construct instead, superimposed on the
nanobody so that the paratope faces the viewer identically in all three.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

from .occlusion import read_pdb


def _write(atoms, path: Path) -> Path:
    f = pdbio.PDBFile()
    f.set_structure(atoms)
    f.write(str(path))
    return path


def aligned_panels(crystal_pdb: Path, construct_pdb: Path, vhh_range,
                   out_dir: Path) -> dict:
    """Superimpose the predicted construct's nanobody onto the crystal nanobody.

    Returns the paths of three coordinate files: the masked construct, the
    nanobody alone with the mask moved away, and the crystal complex.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    crystal = read_pdb(crystal_pdb)
    construct = read_pdb(construct_pdb)

    vhh_crystal = crystal[crystal.chain_id == "A"]
    offset = vhh_range[0] - 1
    con_res = struc.get_residue_starts(construct)
    # CA atoms of the construct's nanobody region
    con_ca = construct[construct.atom_name == "CA"]
    cry_ca = vhh_crystal[vhh_crystal.atom_name == "CA"]
    n = min(len(con_ca) - offset, len(cry_ca))
    fitted, transform = struc.superimpose(cry_ca[:n], con_ca[offset:offset + n])
    construct = transform.apply(construct)

    masked = _write(construct, out_dir / "state_a_masked.pdb")
    complexed = _write(crystal, out_dir / "state_c_bound.pdb")

    # the released state: translate the mask and linker clear of the nanobody
    released = construct.copy()
    res_ids = np.array([a.res_id for a in released])
    mask_sel = res_ids <= (vhh_range[0] - 1)
    if mask_sel.any():
        direction = np.array([1.0, 0.6, 0.0])
        direction /= np.linalg.norm(direction)
        released.coord[mask_sel] += direction * 34.0
    cut = _write(released, out_dir / "state_b_cleaved.pdb")
    return {"masked": masked, "cleaved": cut, "bound": complexed}
