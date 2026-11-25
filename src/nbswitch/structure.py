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

from .config import (
    CONTACT_CUTOFF_A, MASK_CONTACT_CUTOFF_A, CDR_SCHEME, DATA, PDB_ID,
    BURIAL_THRESHOLD_A2, SASA_PROBE_A, SASA_POINTS)

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


def contact_residues(partner_a, partner_b, cutoff: float) -> tuple[np.ndarray, np.ndarray]:
    """Minimum heavy-atom distance from each residue of each partner to the other.

    Hydrogens are excluded: deposited models may carry riding hydrogens, and a
    donor-acceptor criterion is defined between heavy atoms. A residue is in
    contact when this distance is at or below the cutoff.
    """
    def _min_dist(this, other):
        this_heavy = this[this.element != "H"]
        other_heavy = other[other.element != "H"]
        res = residues(this)
        heavy_mask = this.element != "H"
        index = res.atom_index[heavy_mask]
        d = np.linalg.norm(
            this_heavy.coord[:, None, :] - other_heavy.coord[None, :, :], axis=-1).min(axis=1)
        out = np.full(len(res), np.inf)
        np.minimum.at(out, index, d)
        return out

    return _min_dist(partner_a, partner_b), _min_dist(partner_b, partner_a)


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
    """Two descriptions of the same interface, kept separate because they answer
    different questions.

    The polar contacts are residues with a heavy atom within CONTACT_CUTOFF_A of
    the partner: hydrogen-bond and salt-bridge range, so this is the set that is
    "interacting" in the chemical sense. The footprint is every residue that
    loses solvent-accessible surface on binding, which also counts van der Waals
    packing. A mask has to cover the footprint, not just the polar core.
    """
    a, b = atoms[atoms.chain_id == vhh_chain], atoms[atoms.chain_id == toxin_chain]
    res_a, res_b = residues(a), residues(b)
    delta_a, delta_b = buried_on_binding(a, b)
    dist_a, dist_b = contact_residues(a, b, CONTACT_CUTOFF_A)
    para = [i for i in range(len(res_a)) if delta_a[i] > BURIAL_THRESHOLD_A2]
    epi = [i for i in range(len(res_b)) if delta_b[i] > BURIAL_THRESHOLD_A2]
    contact_a = [i for i in range(len(res_a)) if dist_a[i] <= CONTACT_CUTOFF_A]
    contact_b = [i for i in range(len(res_b)) if dist_b[i] <= CONTACT_CUTOFF_A]
    # the wider cutoff the mask is designed from, reported on both partners so
    # the figures and the viewer can show the same set the mask was taken from
    wide_a, wide_b = contact_residues(a, b, MASK_CONTACT_CUTOFF_A)
    near_a = [i for i in range(len(res_a)) if wide_a[i] <= MASK_CONTACT_CUTOFF_A]
    near_b = [i for i in range(len(res_b)) if wide_b[i] <= MASK_CONTACT_CUTOFF_A]
    total = float(delta_a.sum() + delta_b.sum())
    return {
        "toxin_chain": toxin_chain,
        "contact_cutoff_A": CONTACT_CUTOFF_A,
        "n_paratope_contacts": len(contact_a),
        "n_epitope_contacts": len(contact_b),
        "paratope_contact_labels": [res_a.labels[i] for i in contact_a],
        "epitope_contact_labels": [res_b.labels[i] for i in contact_b],
        "mask_contact_cutoff_A": MASK_CONTACT_CUTOFF_A,
        "n_paratope_near": len(near_a), "n_epitope_near": len(near_b),
        "paratope_near_labels": [res_a.labels[i] for i in near_a],
        "epitope_near_labels": [res_b.labels[i] for i in near_b],
        "paratope_min_distance": {res_a.labels[i]: round(float(dist_a[i]), 2)
                                  for i in range(len(res_a)) if np.isfinite(dist_a[i])
                                  and dist_a[i] <= 8.0},
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


#: Framework 4 closes every VHH. Anything after it in a deposited construct is
#: cloning, not antibody: here a TEV remnant left when a C-terminal tag was cut.
FRAMEWORK4_END = re.compile(r"WG[QRKPE]G[TA][QLM]?[VLM]TVSS")


def trim_expression_tags(sequence: str) -> tuple[str, str]:
    """Cut a deposited VHH back to its own C-terminus.

    Crystallised constructs carry purification tags, and what is left after
    cleavage is deposited as though it were part of the chain. 9RIU ends
    ...VTVSSENLY, where ENLY is the start of the TEV site ENLYFQ. Scoring a
    construct that carries it measures four residues of vector.

    Returns the trimmed sequence and whatever was removed.
    """
    match = None
    for match in FRAMEWORK4_END.finditer(sequence):
        pass
    if match is None:
        return sequence, ""
    return sequence[:match.end()], sequence[match.end():]


def find_cdrs(sequence: str, anchors: dict = CDR_ANCHORS,
              scheme: str = CDR_SCHEME) -> dict:
    """CDR spans under a published numbering scheme.

    Numbering is delegated to riot_na, which implements Kabat, Chothia and
    IMGT. The schemes disagree substantially: on the worked example Kabat
    CDR2 is 17 residues and IMGT CDR2 is 7, so the scheme has to be named
    whenever a CDR count is quoted.

    The regex anchors are a fallback for when riot_na is unavailable. They
    were written before the numbering tool was installed and they do not
    reproduce any single scheme: they match Kabat at CDR2 and IMGT at CDR3.
    """
    try:
        spans = _cdr_spans_riot(sequence, scheme)
    except Exception:
        spans = _cdr_spans_anchors(sequence, anchors)
    return spans


def _cdr_spans_riot(sequence: str, scheme: str) -> dict:
    from riot_na import create_riot_aa, Scheme

    wanted = {"kabat": Scheme.KABAT, "chothia": Scheme.CHOTHIA,
              "imgt": Scheme.IMGT}[scheme.lower()]
    result = create_riot_aa().run_on_sequence("vhh", sequence, scheme=wanted)
    spans = {}
    for attr, name in (("cdr1_aa", "CDR1"), ("cdr2_aa", "CDR2"), ("cdr3_aa", "CDR3")):
        sub = getattr(result, attr, None)
        if not sub:
            raise ValueError(f"{name} not assigned under {scheme}")
        start = sequence.find(sub)
        if start < 0:
            raise ValueError(f"{name} sequence not found in the chain")
        spans[name] = (start, start + len(sub))
    return spans


def _cdr_spans_anchors(sequence: str, anchors: dict) -> dict:
    """Framework-motif fallback. See find_cdrs for why this is not preferred."""
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


#: Longest S-S separation still counted as a bond. A disulphide sits near
#: 2.05 A; the next-nearest sulphur pairs in these structures are past 4 A, so
#: the cut is wide of the bonded distances and well clear of the unbonded ones.
DISULPHIDE_MAX_A = 2.5


def disulphides(atoms, max_distance: float = DISULPHIDE_MAX_A) -> list[dict]:
    """Cysteine pairs joined by a disulphide, measured SG to SG.

    Returned sorted, each as the two chains and author residue numbers and the
    distance between the sulphurs, so the viewer draws the bonds the structure
    actually contains rather than every cysteine it happens to hold.
    """
    sg = atoms[(atoms.res_name == "CYS") & (atoms.atom_name == "SG")]
    found = []
    for i in range(sg.array_length()):
        for j in range(i + 1, sg.array_length()):
            d = float(np.linalg.norm(sg.coord[i] - sg.coord[j]))
            if d <= max_distance:
                found.append({
                    "chain_a": str(sg.chain_id[i]), "res_a": int(sg.res_id[i]),
                    "chain_b": str(sg.chain_id[j]), "res_b": int(sg.res_id[j]),
                    "distance_A": round(d, 2)})
    return sorted(found, key=lambda b: (b["chain_a"], b["res_a"]))


#: Widest SG-SG separation still drawn as a bridge in a predicted structure.
#: The deposited complex holds its five bridges at 2.01-2.04 A, but neither
#: ESMFold nor Boltz-2 constrains the bond, and the same bridges come out of a
#: prediction between 2.1 and 2.5 A. Past this the two sulphurs are too far
#: apart for a stick drawn between them to describe anything.
DISULPHIDE_PREDICTED_MAX_A = 2.6


def map_disulphides(atoms, bonds, chain_map=None, offset: int = 0,
                    max_distance: float = DISULPHIDE_PREDICTED_MAX_A) -> list[dict]:
    """Reference disulphides located in another structure of the same chains.

    Which cysteine pairs are bridged is carried over from the reference rather
    than rediscovered here. Rediscovering it from SG-SG distance alone finds
    pairs a prediction has fused together by accident: a mask carries
    cysteines out of the antigen loop it was cut from, their partners stay
    behind in the rest of the antigen, and a predicted structure is free to
    leave the orphans an angstrom apart. Those are not bridges, and drawing
    them would put a covalent bond on the page that the construct does not
    have.

    ``chain_map`` renames reference chains onto the target's, ``offset`` shifts
    residue numbering for a construct carrying a mask and linker ahead of the
    nanobody. Each mapped pair is then measured in the target, so a bridge the
    prediction left open is dropped instead of drawn across empty space.
    """
    names = chain_map or {}
    sg = atoms[(atoms.res_name == "CYS") & (atoms.atom_name == "SG")]
    at = {(str(c), int(r)): sg.coord[i]
          for i, (c, r) in enumerate(zip(sg.chain_id, sg.res_id))}
    found = []
    for b in bonds:
        key_a = (names.get(b["chain_a"], b["chain_a"]), b["res_a"] + offset)
        key_b = (names.get(b["chain_b"], b["chain_b"]), b["res_b"] + offset)
        if key_a not in at or key_b not in at:
            continue
        d = float(np.linalg.norm(at[key_a] - at[key_b]))
        if d > max_distance:
            continue
        found.append({"chain_a": key_a[0], "res_a": key_a[1],
                      "chain_b": key_b[0], "res_b": key_b[1],
                      "distance_A": round(d, 2)})
    return sorted(found, key=lambda b: (b["chain_a"], b["res_a"]))
