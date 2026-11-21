"""Stage 4. Build the N-L-M construct library to be folded.

Masks are derived from the structure itself: alpha-cobratoxin's own
VHH-contacting segment is the one peptide known to engage this paratope, so a
grafted copy of it is a biomimetic mask candidate requiring no generative model.
Scrambled and truncated variants act as negative and intermediate controls.
Linker ladder spans the lengths where stage 2 says the design window is open.
"""
import numpy as np, random, json
from common import DATA, PDB_ID, save_json, load_json
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

iface = load_json("s02_interface.json")
VHH = iface["sequence"]
ref = iface["chain_pairs"]["A"]

f = pdbx.CIFFile.read(DATA / f"{PDB_ID}.cif")
arr = pdbx.get_structure(f, model=1, use_author_fields=True)
arr = arr[struc.filter_amino_acids(arr)]
ta = arr[arr.chain_id == "X"]
starts = struc.get_residue_starts(ta)
tox = [(int(ta.res_id[i]), struc.info.one_letter_code(ta.res_name[i]) or "X") for i in starts]
tox_seq = "".join(c for _, c in tox)
epi = set(ref["epitope_labels"])
is_epi = [f"{c}{r}" in epi for r, c in tox]

# contiguous epitope segments, allowing single-residue gaps
segs, cur = [], []
gap = 0
for i, e in enumerate(is_epi):
    if e:
        cur.append(i); gap = 0
    elif cur:
        gap += 1
        if gap > 1:
            segs.append((cur[0], cur[-1])); cur = []; gap = 0
if cur: segs.append((cur[0], cur[-1]))
segs = [(a, b) for a, b in segs if b - a + 1 >= 3]
seg_info = [{"start_res": tox[a][0], "end_res": tox[b][0], "seq": tox_seq[a:b+1],
             "length": b - a + 1,
             "n_epitope": sum(is_epi[a:b+1])} for a, b in segs]
seg_info.sort(key=lambda s: -s["n_epitope"])

rng = random.Random(42)   # seed fixed for reproducibility
def scramble(s):
    l = list(s); rng.shuffle(l); return "".join(l)

primary = seg_info[0]["seq"]
masks = {
    "epi1": {"seq": primary, "kind": "epitope graft",
             "note": f"cobratoxin {seg_info[0]['start_res']}-{seg_info[0]['end_res']}, "
                     f"{seg_info[0]['n_epitope']} epitope residues"},
    "epi1_scram": {"seq": scramble(primary), "kind": "scrambled control",
                   "note": "same composition, shuffled (seed 42)"},
    "epi1_half": {"seq": primary[: max(3, len(primary) // 2)], "kind": "truncated",
                  "note": "N-terminal half of the graft"},
}
if len(seg_info) > 1:
    masks["epi2"] = {"seq": seg_info[1]["seq"], "kind": "epitope graft",
                     "note": f"cobratoxin {seg_info[1]['start_res']}-{seg_info[1]['end_res']}, "
                             f"{seg_info[1]['n_epitope']} epitope residues"}
masks["polyG"] = {"seq": "GGGGGGGG", "kind": "inert control",
                  "note": "no designed interaction; masking should fail"}

def linker(n):
    return ("GGGGS" * ((n // 5) + 1))[:n]

LINKERS = [5, 10, 15, 20, 30]
constructs = []
for mk, m in masks.items():
    for n in LINKERS:
        L = linker(n)
        constructs.append({
            "id": f"{mk}_L{n}", "mask": mk, "mask_seq": m["seq"], "mask_kind": m["kind"],
            "linker_residues": n, "linker_seq": L,
            "sequence": m["seq"] + L + VHH,
            "length": len(m["seq"]) + n + len(VHH),
            "mask_range": [1, len(m["seq"])],
            "linker_range": [len(m["seq"]) + 1, len(m["seq"]) + n],
            "vhh_range": [len(m["seq"]) + n + 1, len(m["seq"]) + n + len(VHH)],
        })
constructs.append({"id": "VHH_only", "mask": None, "mask_seq": "", "mask_kind": "reference",
                   "linker_residues": 0, "linker_seq": "", "sequence": VHH,
                   "length": len(VHH), "mask_range": None, "linker_range": None,
                   "vhh_range": [1, len(VHH)]})

out = {"vhh_sequence": VHH, "vhh_length": len(VHH),
       "toxin_sequence": tox_seq,
       "epitope_segments": seg_info, "masks": masks,
       "linker_lengths": LINKERS, "linker_motif": "GGGGS repeat",
       "paratope_labels": ref["paratope_labels"],
       "n_constructs": len(constructs), "constructs": constructs}
save_json("s05_constructs.json", out)

print(f"toxin chain X, {len(tox_seq)} residues; {sum(is_epi)} epitope residues")
print("\ncontiguous epitope segments (>=3 residues):")
for s in seg_info:
    print(f"  {s['start_res']:>3}-{s['end_res']:<3} len {s['length']:2d}  "
          f"{s['n_epitope']} contacts   {s['seq']}")
print("\nmask candidates:")
for k, m in masks.items():
    print(f"  {k:12s} {len(m['seq']):2d} aa  {m['seq']:<22s} {m['kind']} — {m['note']}")
print(f"\n{len(constructs)} constructs "
      f"({len(masks)} masks x {len(LINKERS)} linkers + free VHH), "
      f"lengths {min(c['length'] for c in constructs)}-{max(c['length'] for c in constructs)} aa")
