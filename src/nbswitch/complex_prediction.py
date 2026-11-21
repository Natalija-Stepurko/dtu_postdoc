"""Predicting the nanobody-toxin complex with Boltz-2.

ESMFold answers whether the mask folds back onto the paratope. It cannot answer
the question that actually decides the design: with the mask in place, can the
toxin still bind? That needs a two-chain prediction.

Sequences are run in single-sequence mode. Boltz can fetch alignments from a
remote MSA server, but designed sequences would then leave the machine, so the
alignment is left empty and the accuracy cost is accepted and recorded.
"""
from __future__ import annotations
import subprocess
from pathlib import Path

import numpy as np
import biotite.structure as struc

CACHE = Path("/scratch/.boltz")


def write_job(path: Path, chains: dict[str, str]) -> Path:
    """One YAML per case. ``chains`` maps chain id to sequence."""
    lines = ["version: 1", "sequences:"]
    for chain_id, seq in chains.items():
        lines += [f"  - protein:", f"      id: {chain_id}",
                  f"      sequence: {seq}", f"      msa: empty"]
    path.write_text("\n".join(lines) + "\n")
    return path


def predict(job_yaml: Path, out_dir: Path, binary: str = "boltz",
            recycling: int = 3, samples: int = 1, seed: int = 42,
            accelerator: str = "cpu") -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "predict", str(job_yaml), "--out_dir", str(out_dir),
           "--cache", str(CACHE), "--accelerator", accelerator, "--devices", "1",
           "--recycling_steps", str(recycling), "--diffusion_samples", str(samples),
           "--seed", str(seed), "--output_format", "pdb", "--override"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    hits = sorted((out_dir / f"boltz_results_{job_yaml.stem}" / "predictions"
                   / job_yaml.stem).glob("*.pdb"))
    if not hits:
        raise FileNotFoundError(f"no prediction written for {job_yaml.stem}")
    return hits[0]


def interface_contacts(pdb_path: Path, chain_a: str, chain_b: str,
                       cutoff: float = 5.0) -> dict:
    """How much of the two chains actually touch in the predicted pose."""
    import biotite.structure.io.pdb as pdb
    atoms = pdb.PDBFile.read(str(pdb_path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    a = atoms[atoms.chain_id == chain_a]
    b = atoms[atoms.chain_id == chain_b]
    if a.array_length() == 0 or b.array_length() == 0:
        return {"atom_contacts": 0, "buried_surface_area_A2": 0.0}
    d = np.linalg.norm(a.coord[:, None, :] - b.coord[None, :, :], axis=-1)
    free = struc.sasa(a).sum() + struc.sasa(b).sum()
    bound = struc.sasa(a + b).sum()
    return {"atom_contacts": int((d < cutoff).sum()),
            "min_distance_A": round(float(d.min()), 2),
            "buried_surface_area_A2": round(float(np.nan_to_num(free - bound)), 1)}


def toxin_paratope_overlap(pdb_path: Path, fusion_chain: str, toxin_chain: str,
                           paratope_positions) -> dict:
    """Does the toxin land on the paratope, or somewhere else on the nanobody?"""
    import biotite.structure.io.pdb as pdb
    atoms = pdb.PDBFile.read(str(pdb_path)).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    fusion = atoms[atoms.chain_id == fusion_chain]
    toxin = atoms[atoms.chain_id == toxin_chain]
    if toxin.array_length() == 0:
        return {"paratope_contacts": 0, "fraction_of_contacts_on_paratope": 0.0}
    d = np.linalg.norm(fusion.coord[:, None, :] - toxin.coord[None, :, :], axis=-1)
    close = d.min(axis=1) < 5.0
    contacting = set(fusion.res_id[close].tolist())
    para = set(int(p) for p in paratope_positions)
    on_para = contacting & para
    return {"paratope_contacts": len(on_para),
            "total_contacting_residues": len(contacting),
            "fraction_of_contacts_on_paratope":
                round(len(on_para) / max(len(contacting), 1), 3)}
