"""Reading the reference complex and measuring its interface.

The deposited VHH is Kabat-numbered, so residues carry insertion codes: 122
residues share only 112 distinct ``res_id`` values. Every routine here keys
residues on ``(res_id, ins_code)``; keying on ``res_id`` alone silently merges
residues and truncates CDR3.
"""
from __future__ import annotations
import gzip, hashlib, re, urllib.request
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

from .config import (DATA, PDB_ID, BURIAL_THRESHOLD_A2, SASA_PROBE_A, SASA_POINTS)

RCSB = "https://files.rcsb.org/download/{pdb_id}.cif.gz"

#: Conserved VHH framework motifs bracketing each CDR. Anchoring on these makes
#: the CDR boundaries reproducible from sequence alone. CDR2 is bracketed by the
#: end of framework 2 and the start of framework 3; the framework-2 motif varies
#: between VHHs, so its anchor is found by searching back from framework 3 for
#: the last "VA" rather than by matching one family's sequence.
CDR_ANCHORS = {"CDR1": (r"SC[AVTG][AVTGS]S", r"W[YFVG]RQ"),
               "CDR2": (None, r"RF[TAVS][IVL]S"),
               "CDR3": (r"TA[VILMF][YFH][YFHILV]C", r"WG[QRKP]G")}


def fetch_structure(pdb_id: str = PDB_ID, directory: Path = DATA) -> Path:
    """Download a mmCIF from RCSB once and cache it."""
    path = directory / f"{pdb_id}.cif"
    if not path.exists():
        url = RCSB.format(pdb_id=pdb_id)
        with urllib.request.urlopen(url, timeout=60) as response:
            path.write_bytes(gzip.decompress(response.read()))
    return path


def file_digest(path: Path, length: int = 16) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:length]


def read_metadata(path: Path) -> dict:
    block = pdbx.CIFFile.read(path).block
    return {"title": block["struct"]["title"].as_item().strip(),
            "method": block["exptl"]["method"].as_item(),
            "resolution_A": float(block["refine"]["ls_d_res_high"].as_item()),
            "deposition_date":
                block["pdbx_database_status"]["recvd_initial_deposition_date"].as_item()}


def read_atoms(path: Path):
    """Amino-acid atoms of model 1, author numbering."""
    atoms = pdbx.get_structure(pdbx.CIFFile.read(path), model=1, use_author_fields=True)
    return atoms[struc.filter_amino_acids(atoms)]


@dataclass
class Residues:
    """Residues of one chain, keyed on (res_id, ins_code)."""
    res_id: np.ndarray
    ins_code: np.ndarray
    code: list[str]
    atom_index: np.ndarray = field(repr=False)

    def __len__(self) -> int:
        return len(self.code)

    @property
    def labels(self) -> list[str]:
        return [f"{c}{r}{i}" for c, r, i in zip(self.code, self.res_id, self.ins_code)]

    @property
    def sequence(self) -> str:
        return "".join(self.code)


def residues(atoms) -> Residues:
    starts = struc.get_residue_starts(atoms)
    index = np.zeros(atoms.array_length(), dtype=int)
    for n, start in enumerate(starts):
        end = starts[n + 1] if n + 1 < len(starts) else atoms.array_length()
        index[start:end] = n
    return Residues(
        res_id=atoms.res_id[starts], ins_code=atoms.ins_code[starts],
        code=[struc.info.one_letter_code(name) or "X" for name in atoms.res_name[starts]],
        atom_index=index)


def residue_sasa(atoms) -> np.ndarray:
    """Solvent-accessible surface per residue (Shrake-Rupley)."""
    per_atom = np.nan_to_num(
        struc.sasa(atoms, probe_radius=SASA_PROBE_A, point_number=SASA_POINTS))
    res = residues(atoms)
    out = np.zeros(len(res))
    np.add.at(out, res.atom_index, per_atom)
    return out


def buried_on_binding(partner_a, partner_b) -> tuple[np.ndarray, np.ndarray]:
    """Surface each residue of both partners loses when the complex forms."""
    free_a, free_b = residue_sasa(partner_a), residue_sasa(partner_b)
    complexed = residue_sasa(partner_a + partner_b)
    return free_a - complexed[:len(free_a)], free_b - complexed[len(free_a):]


def pair_chains(atoms, vhh_chains, toxin_chains, cutoff_A: float = 5.0) -> dict:
    """Match each VHH copy to the toxin copy it contacts most."""
    pairs = {}
    for vhh in vhh_chains:
        a = atoms[atoms.chain_id == vhh]
        counts = {}
        for toxin in toxin_chains:
            b = atoms[atoms.chain_id == toxin]
            d = np.linalg.norm(a.coord[:, None, :] - b.coord[None, :, :], axis=-1)
            counts[toxin] = int((d < cutoff_A).sum())
        best = max(counts, key=counts.get)
        pairs[vhh] = {"toxin_chain": best, "atom_contacts_5A": counts[best]}
    return pairs


def interface(atoms, vhh_chain: str, toxin_chain: str) -> dict:
    a, b = atoms[atoms.chain_id == vhh_chain], atoms[atoms.chain_id == toxin_chain]
    res_a, res_b = residues(a), residues(b)
    delta_a, delta_b = buried_on_binding(a, b)
    para = [i for i in range(len(res_a)) if delta_a[i] > BURIAL_THRESHOLD_A2]
    epi = [i for i in range(len(res_b)) if delta_b[i] > BURIAL_THRESHOLD_A2]
    total = float(delta_a.sum() + delta_b.sum())
    return {
        "toxin_chain": toxin_chain,
        "buried_surface_area_A2": round(total, 1),
        "interface_area_A2": round(total / 2, 1),
        "n_paratope": len(para), "n_epitope": len(epi),
        "paratope_labels": [res_a.labels[i] for i in para],
        "epitope_labels": [res_b.labels[i] for i in epi],
        "paratope_delta_sasa": {res_a.labels[i]: round(float(delta_a[i]), 2) for i in para},
        "vhh_delta_sasa_all": [round(float(x), 2) for x in delta_a],
        "vhh_residue_labels": res_a.labels,
        "toxin_residue_labels": res_b.labels,
        "toxin_delta_sasa_all": [round(float(x), 2) for x in delta_b],
    }


def find_cdrs(sequence: str, anchors: dict = CDR_ANCHORS) -> dict:
    """CDR spans as half-open 0-based intervals, from framework anchors.

    A ``None`` start anchor means the span begins at the last ``VA`` before the
    end anchor, which is how framework 2 closes in every VHH family checked.
    """
    spans = {}
    for name, (before, after) in anchors.items():
        end_match = None
        if before is None:
            end_match = re.search(after, sequence)
            if end_match is None:
                raise ValueError(f"no framework anchor {after!r} for {name}")
            prefix = sequence[:end_match.start()]
            va = list(re.finditer(r"VA", prefix))
            if not va:
                raise ValueError(f"no framework-2 anchor before {name}")
            start = va[-1].end()
        else:
            m = re.search(before, sequence)
            if m is None:
                raise ValueError(f"no framework anchor {before!r} for {name}")
            start = m.end()
        if end_match is None:
            m = re.search(after, sequence[start:])
            if m is None:
                raise ValueError(f"no framework anchor {after!r} for {name}")
            stop = m.start() + start
        else:
            stop = end_match.start()
        spans[name] = (start, stop)
    return spans


def terminal_geometry(atoms, vhh_chain: str, paratope_labels) -> dict:
    """Where the termini sit relative to the paratope - this sets linker length."""
    chain = atoms[atoms.chain_id == vhh_chain]
    ca = chain[chain.atom_name == "CA"]
    labels = residues(ca).labels
    mask = np.array([label in set(paratope_labels) for label in labels])
    patch = ca.coord[mask]
    centroid = patch.mean(axis=0)
    n_term, c_term = ca.coord[0], ca.coord[-1]
    return {
        "reference_chain": vhh_chain,
        "n_term_residue": labels[0], "c_term_residue": labels[-1],
        "n_term_to_paratope_centroid_A": round(float(np.linalg.norm(n_term - centroid)), 2),
        "c_term_to_paratope_centroid_A": round(float(np.linalg.norm(c_term - centroid)), 2),
        "n_term_to_nearest_paratope_A":
            round(float(np.min(np.linalg.norm(patch - n_term, axis=1))), 2),
        "c_term_to_nearest_paratope_A":
            round(float(np.min(np.linalg.norm(patch - c_term, axis=1))), 2),
        "n_to_c_term_A": round(float(np.linalg.norm(n_term - c_term)), 2),
        "paratope_max_extent_A": round(float(
            np.max(np.linalg.norm(patch[:, None] - patch[None], axis=-1))), 2),
        "radius_of_gyration_A": round(float(struc.gyration_radius(chain)), 2),
    }
