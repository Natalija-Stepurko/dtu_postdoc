#!/usr/bin/env python
"""Stage 15. Is the worked example's interface typical of a VHH complex?

The pipeline measures one interface and draws conclusions from it. This stage
measures the same quantities on independent VHH-antigen complexes, so the
worked example can be read against them rather than on its own.
"""
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import biotite.structure as struc

from nbswitch import structure as S
from nbswitch.config import save_json

#: Independent complexes. 9GCN and 9S9E are VHH-three-finger-toxin structures
#: from other groups; 9EA9 is a VHH bound to a much larger toxin, included so
#: the range is not all one antigen class.
COMPARATORS = ("9RIU", "9GCN", "9S9E", "9EA9")

VHH_START = re.compile(r"^[QEDKA]V[QKE]L")
VHH_END = re.compile(r"WG[QRKPE]G")


def vhh_chains(atoms):
    """Chains that look like a VHH: right length, and both framework anchors."""
    found = []
    for chain in sorted(set(atoms.chain_id)):
        residues = S.residues(atoms[atoms.chain_id == chain])
        if not 100 <= len(residues) <= 145:
            continue
        seq = "".join(residues.code)
        if VHH_START.search(seq) and VHH_END.search(seq):
            found.append(chain)
    return found


def best_partner(atoms, vhh):
    """The chain burying the most atoms against this one."""
    a = atoms[(atoms.chain_id == vhh) & (atoms.element != "H")]
    best, most = None, 0
    for chain in sorted(set(atoms.chain_id)):
        if chain == vhh:
            continue
        other = atoms[(atoms.chain_id == chain) & (atoms.element != "H")]
        if other.array_length() == 0:
            continue
        d = np.linalg.norm(a.coord[:, None, :] - other.coord[None, :, :], axis=-1)
        n = int((d < 5.0).sum())
        if n > most:
            best, most = chain, n
    return best, most


rows = []
for pdb in COMPARATORS:
    atoms = S.read_atoms(S.fetch_structure(pdb))
    chains = vhh_chains(atoms)
    if not chains:
        continue
    vhh = chains[0]
    antigen, _ = best_partner(atoms, vhh)
    if antigen is None:
        continue
    iface = S.interface(atoms, vhh, antigen)
    a = atoms[atoms.chain_id == vhh]
    residues = S.residues(a)
    sse = struc.annotate_sse(a)
    cdrs = S.find_cdrs("".join(residues.code))
    cdr3 = {residues.labels[i] for i in range(*cdrs["CDR3"])}
    delta = iface["paratope_delta_sasa"]
    total = sum(delta.values()) or 1.0
    rows.append({
        "pdb_id": pdb,
        "vhh_chain": vhh,
        "antigen_chain": antigen,
        "interface_area_A2": iface["interface_area_A2"],
        "n_paratope": iface["n_paratope"],
        "n_contacts": iface["n_paratope_contacts"],
        "cdr3_percent_bsa": round(
            100 * sum(v for k, v in delta.items() if k in cdr3) / total, 1),
        "loop_percent_bsa": round(
            100 * sum(v for k, v in delta.items()
                      if sse[residues.labels.index(k)] != "b") / total, 1),
    })
    r = rows[-1]
    print(f"  {pdb} {vhh}:{antigen}  {r['interface_area_A2']:>6.0f} A2  "
          f"{r['n_paratope']:>2d} paratope  {r['n_contacts']:>2d} within 3 A  "
          f"CDR3 {r['cdr3_percent_bsa']:>4.0f}%  loops {r['loop_percent_bsa']:>4.0f}%")

save_json("15_interface_controls.json", {"comparators": rows})
print(f"{len(rows)} complexes compared")
