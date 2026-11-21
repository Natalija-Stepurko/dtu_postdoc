#!/usr/bin/env python
"""Stage 3. Derive mask candidates from the structure and enumerate constructs."""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
from nbswitch import constructs as C
from nbswitch.structure import trim_expression_tags
from nbswitch.config import load_json, save_json, REFERENCE_VHH_CHAIN

ref_data = load_json("01_reference.json")
ref = ref_data["chain_pairs"][REFERENCE_VHH_CHAIN]
toxin_seq = ref_data["metadata"]["chains"][ref["toxin_chain"]]["sequence"]

# the deposited chain can carry what is left of a purification tag; that is
# vector, not antibody, and must not be folded or scored as part of the binder
vhh_seq, removed = trim_expression_tags(ref_data["sequence"])
if removed:
    print(f"  trimmed {removed!r} from the C-terminus "
          f"({len(ref_data['sequence'])} -> {len(vhh_seq)} residues)")

segments = C.epitope_segments(ref["toxin_residue_labels"], toxin_seq, ref["epitope_labels"])
masks = C.mask_candidates(segments)
library = C.build_library(vhh_seq, masks)

save_json("03_constructs.json", {
    "vhh_sequence": vhh_seq, "vhh_sequence_deposited": ref_data["sequence"],
    "trimmed_from_c_terminus": removed, "toxin_sequence": toxin_seq,
    "paratope_labels": ref["paratope_labels"],
    "vhh_residue_labels": ref["vhh_residue_labels"],
    "epitope_segments": segments, "masks": masks,
    "linker_lengths": list(C.DEFAULT_LINKER_LENGTHS), "linker_motif": C.LINKER_MOTIF,
    "n_constructs": len(library), "constructs": library})

print(f"epitope segments (>=3 residues): {len(segments)}")
for s in segments:
    print(f"  {s['start_label']}-{s['end_label']:<6} {s['seq']:<8} "
          f"{s['n_epitope']} contacts")
print(f"\nmasks: {len(masks)}")
for name, m in masks.items():
    print(f"  {name:12s} {m['seq']:<10s} {m['kind']}")
print(f"\n{len(library)} constructs, "
      f"{min(c['length'] for c in library)}-{max(c['length'] for c in library)} residues")
