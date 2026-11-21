#!/usr/bin/env python
"""Stage 1. Measure the reference nanobody-antigen complex."""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from nbswitch import structure as S
from nbswitch.config import (PDB_ID, VHH_CHAINS, TOXIN_CHAINS, REFERENCE_VHH_CHAIN,
                             save_json)

path = S.fetch_structure(PDB_ID)
atoms = S.read_atoms(path)
meta = S.read_metadata(path) | {"pdb_id": PDB_ID, "sha256_16": S.file_digest(path),
                                "source": S.RCSB.format(pdb_id=PDB_ID)}

chains = {}
for chain in sorted(set(atoms.chain_id)):
    res = S.residues(atoms[atoms.chain_id == chain])
    chains[chain] = {"n_residues": len(res), "sequence": res.sequence,
                     "first": res.labels[0], "last": res.labels[-1]}
meta["chains"] = chains

pairs = S.pair_chains(atoms, VHH_CHAINS, TOXIN_CHAINS)
interfaces = {v: S.interface(atoms, v, p["toxin_chain"]) | {"atom_contacts_5A": p["atom_contacts_5A"]}
              for v, p in pairs.items()}

ref = interfaces[REFERENCE_VHH_CHAIN]
res = S.residues(atoms[atoms.chain_id == REFERENCE_VHH_CHAIN])
spans = S.find_cdrs(res.sequence)
cdrs = {name: {"seq": res.sequence[a:b], "length": b - a,
               "residues": res.labels[a:b]} for name, (a, b) in spans.items()}

by_region = {}
para = set(ref["paratope_labels"])
for name, (a, b) in spans.items():
    by_region[name] = [l for l in res.labels[a:b] if l in para]
in_cdr = {l for v in by_region.values() for l in v}
by_region["framework"] = [l for l in ref["paratope_labels"] if l not in in_cdr]

geometry = S.terminal_geometry(atoms, REFERENCE_VHH_CHAIN, ref["paratope_labels"])
bsa = [i["buried_surface_area_A2"] for i in interfaces.values()]

save_json("01_reference.json", {
    "metadata": meta, "sequence": res.sequence, "n_residues": len(res),
    "cdrs": cdrs, "chain_pairs": interfaces, "terminal_geometry": geometry,
    "paratope_by_region": {k: {"residues": v, "n": len(v)} for k, v in by_region.items()},
    "reproducibility": {"n_copies": len(bsa), "bsa_mean_A2": round(float(np.mean(bsa)), 1),
                        "bsa_sd_A2": round(float(np.std(bsa, ddof=1)), 1)}})

print(f"{meta['title']}")
print(f"  {meta['method'].lower()}, {meta['resolution_A']} A, deposited {meta['deposition_date']}")
print(f"  {len(res)} residues, CDR3 length {cdrs['CDR3']['length']}")
print(f"  interface {ref['interface_area_A2']} A2, paratope {ref['n_paratope']}, "
      f"epitope {ref['n_epitope']}")
print(f"  BSA across {len(bsa)} copies: {np.mean(bsa):.1f} +/- {np.std(bsa, ddof=1):.1f} A2")
print(f"  paratope by region: " + ", ".join(f"{k} {len(v)}" for k, v in by_region.items()))
print(f"  N-term to paratope {geometry['n_term_to_paratope_centroid_A']} A, "
      f"C-term {geometry['c_term_to_paratope_centroid_A']} A")
