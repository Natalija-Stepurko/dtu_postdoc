#!/usr/bin/env python
"""Stage 3. Assemble constructs from the designed masks and the designed linkers.

This stage designs nothing. Masks come from stage 03a and linkers from stage
02b, and the only work here is enumerating the combinations and recording where
each part sits in the fused sequence. Keeping assembly separate means a change
to how masks are designed cannot quietly alter the linkers, or the reverse.
"""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
from nbswitch import constructs as C
from nbswitch import linkers as L

#: Linker sequences carried through the whole length range. One
#: generated substrate and both controls, so each length is represented
#: in every sequence class.
LENGTH_LINKERS = ("cn01", "ctrl_gs", "ctrl_scram")

#: Masks carried through the whole length range. One graft, its scramble and
#: the inert control: enough to say whether a gain with length belongs to the
#: mask or to anything of that size being nearby.
LENGTH_MASKS = ("epi_core", "epi_core_scram", "polyG")
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

mask_data = load_json("03a_masks.json")
segments, masks = mask_data["segments"], mask_data["masks"]
linker_data = load_json("02b_linkers.json")
linkers = linker_data["linkers"]
reference_length = linker_data.get("reference_length") or L.REFERENCE_LENGTH
# every mask against every linker at every length. The subsetted matrix that
# stood here traded coverage for folding time, which is the wrong trade when
# the question is what the library actually contains.
pairs = [(name, link) for name in masks for link in linkers]
print(f"  full factorial: {len(masks)} masks x {len(linkers)} linkers "
      f"= {len(pairs)} pairs, about {(len(pairs) + 1) * 50 / 60:.0f} minutes "
      f"to fold on CPU")

library = C.build_library(vhh_seq, masks, linkers, pairs)

save_json("03_constructs.json", {
    "vhh_sequence": vhh_seq, "vhh_sequence_deposited": ref_data["sequence"],
    "trimmed_from_c_terminus": removed, "toxin_sequence": toxin_seq,
    "paratope_labels": ref["paratope_labels"],
    "vhh_residue_labels": ref["vhh_residue_labels"],
    "epitope_segments": segments, "masks": masks,
    "linkers": linkers, "n_linkers": len(linkers),
    "n_substrates_generated": linker_data["n_generated"],
    "target_proteases": linker_data["target_proteases"],
    "off_target_proteases": linker_data["off_target_proteases"],
    "reference_length": reference_length,
    "length_linkers": list(LENGTH_LINKERS),
    "n_pairs": len(pairs),
    "n_constructs": len(library), "constructs": library})

print(f"{len(masks)} masks, {len(linkers)} linkers, {len(pairs)} designed pairs "
      f"+ the free nanobody = {len(library)} constructs, "
      f"{min(c['length'] for c in library)}-{max(c['length'] for c in library)} residues")
print(f"  masks from stage 03a ({mask_data['cutoff_A']} A contact seed), "
      f"linkers from stage 02b ({linker_data['n_generated']} substrates generated)")
