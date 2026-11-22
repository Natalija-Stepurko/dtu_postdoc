#!/usr/bin/env python
"""Stage 3a. Design the masks, from the structure alone.

Mask design and linker design answer different questions and were previously
tangled in one stage. A mask has to occupy a particular surface on a particular
nanobody, so it is designed from the complex. A linker has to be cut by a
protease and reach far enough for the mask to arrive, so it is designed from
substrate data and from the tether geometry. Neither constrains the other, and
holding them apart means a change to one does not silently alter the other.
This stage designs masks only; stage 02b designs linkers; stage 03 assembles
every combination.

The starting sequence is what the antigen itself presents at the paratope:
every antigen residue with a heavy atom within 3.5 A of the nanobody, grouped
into runs along the chain. That is wider than the 3 A used elsewhere for polar
contacts, because a mask has to reproduce a whole contacting face rather than
only its hydrogen bonds.

Taking the contact face rather than a surrounding window matters here for a
reason specific to this antigen. The contacts span C42 to C53, and C42-C53 is
itself a disulfide. Cutting there keeps the staple that holds the loop in the
conformation the nanobody binds. The earlier window, C40 to T55, took four
cysteines instead of two, and the two extra ones are bonded to C17 and C59,
outside the graft: they have no partner once excised, and the predicted models
show them unpaired.
"""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from nbswitch import constructs as C
from nbswitch import mask_design as M
from nbswitch import structure as S
from nbswitch.config import (load_json, save_json, REFERENCE_VHH_CHAIN, PDB_ID,
                             MASK_CONTACT_CUTOFF_A, MASK_MAX_GAP, MASK_MIN_RUN)

ref_data = load_json("01_reference.json")
ref = ref_data["chain_pairs"][REFERENCE_VHH_CHAIN]
toxin_chain = ref["toxin_chain"]

atoms = S.read_atoms(S.fetch_structure(PDB_ID))
segments = C.contact_segments(atoms, REFERENCE_VHH_CHAIN, toxin_chain)
if not segments:
    raise SystemExit(f"no antigen residue lies within {MASK_CONTACT_CUTOFF_A} A "
                     f"of the nanobody: nothing to graft")

res = S.residues(atoms[atoms.chain_id == toxin_chain])
by_id = {int(res.res_id[i]): res.code[i] for i in range(len(res))}
label = {int(res.res_id[i]): res.labels[i] for i in range(len(res))}


def window(start, end):
    """Antigen sequence over an inclusive residue range, clipped to the chain."""
    ids = [r for r in range(max(start, min(by_id)), min(end, max(by_id)) + 1)
           if r in by_id]
    return "".join(by_id[r] for r in ids), ids


hotspots = C.hotspot_indices(segments[0])
masks = M.design_masks(window, segments, hotspot_positions=hotspots)


def _disulfide_report(atoms, chain, lo, hi):
    """Which cysteines in a span keep their partner, and which the cut orphans."""
    b = atoms[atoms.chain_id == chain]
    sg = b[(b.atom_name == "SG") & (b.res_name == "CYS")]
    ids = [int(i) for i in sg.res_id]
    inside = set(range(lo, hi + 1))
    internal, orphaned = [], []
    for i, rid in enumerate(ids):
        if rid not in inside:
            continue
        partner = next((ids[j] for j in range(len(ids)) if j != i and
                        float(np.linalg.norm(sg.coord[i] - sg.coord[j])) < 2.5), None)
        if partner is not None and partner in inside:
            internal.append(tuple(sorted((rid, partner))))
        else:
            orphaned.append(rid)
    return sorted(set(internal)), sorted(orphaned)

primary = segments[0]
internal, orphaned = _disulfide_report(atoms, toxin_chain,
                                       primary["start"], primary["end"])

save_json("03a_masks.json", {
    "method": f"antigen residues within {MASK_CONTACT_CUTOFF_A} A of the nanobody, "
              f"runs merged across gaps of up to {MASK_MAX_GAP} residues",
    "cutoff_A": MASK_CONTACT_CUTOFF_A, "max_gap": MASK_MAX_GAP,
    "min_run": MASK_MIN_RUN,
    "toxin_chain": toxin_chain,
    "segments": segments, "masks": masks,
    "hotspot_positions": hotspots,
    "hotspot_labels": [label[segments[0]["start"] + i] for i in hotspots
                       if segments[0]["start"] + i in label],
    "methods_run": sorted({m.get("method", "graft") for m in masks.values()}),
    "methods_available": ["RFdiffusion", "RFantibody", "BindCraft",
                          "ESM-2 masked resampling", "ProGen2"],
    "primary_disulfides_internal": [list(p) for p in internal],
    "primary_cysteines_orphaned": orphaned})

print(f"contact segments at {MASK_CONTACT_CUTOFF_A} A: {len(segments)}")
for s in segments:
    print(f"  {s['start_label']}-{s['end_label']:<6s} {s['length']:2d} residues  "
          f"{s['seq']:<16s} {s['n_contact']} in contact, buries {s['buried_A2']:.0f} A2")
print(f"\nprimary segment {primary['start_label']}-{primary['end_label']}:")
print(f"  disulfides kept whole: {internal if internal else 'none'}")
print(f"  cysteines orphaned by the cut: {orphaned if orphaned else 'none'}")
print(f"\nhotspots held fixed: "
      f"{[label[segments[0]['start'] + i] for i in hotspots if segments[0]['start'] + i in label]}")
print(f"\nmask design space: {len(masks)} candidates seeded by the contact face")
for name, m in masks.items():
    span = ""
    if m.get("start"):
        kept, orph = _disulfide_report(atoms, toxin_chain, m["start"], m["end"])
        span = f"  disulfides {len(kept)}, orphaned {len(orph)}"
    print(f"  {name:18s} {m['seq']:<26s} {len(m['seq']):2d} aa  "
          f"{m.get('method', 'graft'):<14s}{m['kind']:<18s}{span}")
