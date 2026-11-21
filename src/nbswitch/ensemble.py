"""Superimpose every construct on the nanobody and write one multi-model file.

Each construct is predicted independently, so their coordinate frames are
unrelated. Superimposing the nanobody of each onto a common reference puts them
in one frame, and the masks then show where the generator actually places them
relative to the antigen-binding surface.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdbio

from .occlusion import read_pdb


def superimpose_on_nanobody(constructs, structure_dir: Path, reference_id: str,
                            out_path: Path, max_models: int = 60) -> dict:
    """Write a multi-model PDB: one nanobody, and every construct's mask and linker.

    Model 1 holds the reference nanobody. Each later model holds only the mask
    and linker of one construct, already in the reference frame.
    """
    reference = read_pdb(structure_dir / f"{reference_id}.pdb")
    ref_ca = reference[reference.atom_name == "CA"]

    kept, written = [], []
    for item in constructs:
        path = structure_dir / f"{item['id']}.pdb"
        if not path.exists() or not item.get("mask_range"):
            continue
        atoms = read_pdb(path)
        offset = item["vhh_range"][0] - 1
        ca = atoms[atoms.atom_name == "CA"]
        n = min(len(ca) - offset, len(ref_ca))
        if n < 50:
            continue
        _, transform = struc.superimpose(ref_ca[:n], ca[offset:offset + n])
        atoms = transform.apply(atoms)
        res_ids = atoms.res_id
        part = atoms[res_ids <= item["linker_range"][1]]
        if part.array_length() == 0:
            continue
        part.chain_id[:] = "M"
        written.append((item["id"], part))
        kept.append(item["id"])
        if len(written) >= max_models:
            break

    stack = [reference] + [p for _, p in written]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for n, model in enumerate(stack, 1):
            fh.write(f"MODEL     {n:4d}\n")
            tmp = pdbio.PDBFile()
            tmp.set_structure(model)
            for line in tmp.lines:
                if line.startswith(("ATOM", "HETATM")):
                    fh.write(line + "\n")
            fh.write("ENDMDL\n")
        fh.write("END\n")
    return {"path": out_path, "n_models": len(stack), "constructs": kept}
