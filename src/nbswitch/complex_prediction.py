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


def write_tethered_job(path: Path, vhh: str, mask_linker: str, toxin: str,
                       vhh_msa: str = "server", peptide_msa: str = "empty") -> Path:
    """The construct as three chains joined by an explicit covalent bond.

    Modelling the fusion as one chain forces a single alignment over the whole
    sequence, which would mean sending the designed mask to an alignment server.
    Splitting it lets the nanobody and antigen - both public, both already in the
    PDB - take real alignments, while the designed peptide takes none and never
    leaves the machine. The bond constraint keeps the tether.

    Chain A nanobody, chain B mask+linker, chain C antigen. The bond joins the
    peptide's C-terminus to the nanobody's N-terminus.
    """
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {vhh}"]
    if vhh_msa != "server":
        lines.append(f"      msa: {vhh_msa}")
    lines += ["  - protein:", "      id: B", f"      sequence: {mask_linker}",
              f"      msa: {peptide_msa}",
              "  - protein:", "      id: C", f"      sequence: {toxin}"]
    if vhh_msa != "server":
        lines.append(f"      msa: {vhh_msa}")
    lines += ["constraints:",
              "  - bond:",
              f"      atom1: [B, {len(mask_linker)}, C]",
              "      atom2: [A, 1, N]"]
    path.write_text("\n".join(lines) + "\n")
    return path


def predict_with_msa(job_yaml: Path, out_dir: Path, binary: str = "boltz",
                     recycling: int = 3, samples: int = 1, seed: int = 42) -> Path:
    """Run with the alignment server enabled for the chains that request it."""
    import subprocess
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [binary, "predict", str(job_yaml), "--out_dir", str(out_dir),
           "--cache", str(CACHE), "--accelerator", "cpu", "--devices", "1",
           "--recycling_steps", str(recycling), "--diffusion_samples", str(samples),
           "--seed", str(seed), "--output_format", "pdb", "--override",
           "--use_msa_server"]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    hits = sorted((out_dir / f"boltz_results_{job_yaml.stem}" / "predictions"
                   / job_yaml.stem).glob("*.pdb"))
    if not hits:
        raise FileNotFoundError(f"no prediction written for {job_yaml.stem}")
    return hits[0]


def write_binding_job(path: Path, vhh: str, peptide: str,
                      vhh_msa: str = "server", peptide_msa: str = "empty") -> Path:
    """A free peptide against the bare nanobody: two chains, no tether, no antigen.

    This asks the question the tethered prediction cannot. In the construct the
    mask is held at the paratope by its linker whether or not it has any
    affinity for it, so proximity there is not evidence of binding. Cut the
    tether and the peptide has to find the site on its own.

    As elsewhere, the nanobody is a public sequence and takes a real alignment,
    while the designed peptide takes none and never leaves the machine.
    """
    lines = ["version: 1", "sequences:",
             "  - protein:", "      id: A", f"      sequence: {vhh}"]
    if vhh_msa != "server":
        lines.append(f"      msa: {vhh_msa}")
    lines += ["  - protein:", "      id: B", f"      sequence: {peptide}",
              f"      msa: {peptide_msa}"]
    path.write_text("\n".join(lines) + "\n")
    return path


def peptide_disulfides(pdb_path, chain: str = "B", max_distance: float = 2.5) -> dict:
    """Which cysteine pairs are bonded in a predicted peptide.

    A graft cut out of a disulfide-rich antigen can re-form the bonds it still
    has both partners for, pair up cysteines that were bonded to residues left
    behind, or leave them free. Those are three different molecules, and only
    the first is the thing the graft was meant to be.
    """
    import numpy as np
    from .occlusion import read_pdb

    atoms = read_pdb(pdb_path)
    sel = atoms[atoms.chain_id == chain]
    sg = sel[(sel.atom_name == "SG") & (sel.res_name == "CYS")]
    n = sg.array_length()
    if n < 2:
        return {"n_cysteines": int(n), "bonded_pairs": [], "n_unpaired": int(n)}
    d = np.linalg.norm(sg.coord[:, None, :] - sg.coord[None, :, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    pairs, seen = [], set()
    for i in range(n):
        j = int(np.argmin(d[i]))
        if d[i, j] <= max_distance and (j, i) not in seen:
            seen.add((i, j))
            pairs.append([int(sg.res_id[i]), int(sg.res_id[j]),
                          round(float(d[i, j]), 2)])
    bonded = {r for p in pairs for r in p[:2]}
    return {"n_cysteines": int(n), "bonded_pairs": pairs,
            "n_unpaired": int(n - len(bonded))}
