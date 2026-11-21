"""Stage 1b. Paratope/epitope from buried surface area, CDR assignment, terminal geometry.

The deposited VHH uses Kabat insertion codes, so residues are keyed on
(res_id, ins_code); keying on res_id alone silently merges residues and
truncates CDR3. Paratope: VHH residue burying >1 A^2 of solvent-accessible
surface on complex formation (Shrake-Rupley, probe 1.4 A). CDR boundaries are
anchored on conserved VHH framework motifs, reproducible from sequence alone.
"""
import numpy as np, re
from common import DATA, PDB_ID, save_json
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

f = pdbx.CIFFile.read(DATA / f"{PDB_ID}.cif")
arr = pdbx.get_structure(f, model=1, use_author_fields=True)
arr = arr[struc.filter_amino_acids(arr)]
VHH_CHAINS, TOX_CHAINS = ["A", "B", "C", "D"], ["V", "X", "Y", "Z"]

def res_table(a):
    """Ordered unique residues as (res_id, ins_code, one_letter), plus a per-atom index."""
    starts = struc.get_residue_starts(a)
    labels = [(int(a.res_id[i]), str(a.ins_code[i]),
               struc.info.one_letter_code(a.res_name[i]) or "X") for i in starts]
    idx = np.zeros(a.array_length(), dtype=int)
    for n, s in enumerate(starts):
        e = starts[n + 1] if n + 1 < len(starts) else a.array_length()
        idx[s:e] = n
    return labels, idx

def res_sasa(a):
    v = struc.sasa(a, probe_radius=1.4, point_number=500)
    labels, idx = res_table(a)
    out = np.zeros(len(labels))
    np.add.at(out, idx, np.nan_to_num(v))
    return labels, out

def lab(x):
    return f"{x[2]}{x[0]}{x[1]}"

results = {}
for v in VHH_CHAINS:
    va = arr[arr.chain_id == v]
    best, bestn = None, -1
    for t in TOX_CHAINS:
        ta = arr[arr.chain_id == t]
        d = np.linalg.norm(va.coord[:, None, :] - ta.coord[None, :, :], axis=-1)
        n = int((d < 5.0).sum())
        if n > bestn: best, bestn = t, n
    t = best
    ta = arr[arr.chain_id == t]
    lv, fv = res_sasa(va)
    lt, ft = res_sasa(ta)
    comp = va + ta
    lc, fc = res_sasa(comp)
    bv, bt = fc[:len(lv)], fc[len(lv):]
    dv, dt = fv - bv, ft - bt
    para = [i for i in range(len(lv)) if dv[i] > 1.0]
    epi = [i for i in range(len(lt)) if dt[i] > 1.0]
    results[v] = {
        "toxin_chain": t, "atom_contacts_5A": bestn,
        "buried_surface_area_A2": round(float(dv.sum() + dt.sum()), 1),
        "interface_area_A2": round(float((dv.sum() + dt.sum()) / 2), 1),
        "n_residues_vhh": len(lv), "n_residues_toxin": len(lt),
        "paratope_labels": [lab(lv[i]) for i in para],
        "epitope_labels": [lab(lt[i]) for i in epi],
        "n_paratope": len(para), "n_epitope": len(epi),
        "paratope_delta_sasa": {lab(lv[i]): round(float(dv[i]), 2) for i in para},
        "vhh_delta_sasa_all": [round(float(x), 2) for x in dv],
        "vhh_residue_labels": [lab(x) for x in lv],
    }

ref = "A"
va = arr[arr.chain_id == ref]
lv, _ = res_table(va)
seq = "".join(x[2] for x in lv)

def span(p0, p1):
    i = re.search(p0, seq).end()
    j = re.search(p1, seq[i:]).start() + i
    return i, j  # 0-based [i, j)
cdr_raw = {"CDR1": span(r"SC[AV][AV]S", r"W[YF]RQ"),
           "CDR2": span(r"ELVA", r"RFTIS"),
           "CDR3": span(r"YYC", r"WGQG")}
cdrs = {k: {"seq": seq[i:j], "length": j - i,
            "residues": [lab(lv[n]) for n in range(i, j)]}
        for k, (i, j) in cdr_raw.items()}

para_idx = [n for n in range(len(lv)) if lab(lv[n]) in set(results[ref]["paratope_labels"])]
region = {}
for k, (i, j) in cdr_raw.items():
    region[k] = [lab(lv[n]) for n in para_idx if i <= n < j]
region["framework"] = [lab(lv[n]) for n in para_idx
                       if not any(i <= n < j for i, j in cdr_raw.values())]

ca = va[va.atom_name == "CA"]
lca, _ = res_table(ca)
ca_lab = [lab(x) for x in lca]
pmask = np.array([l in set(results[ref]["paratope_labels"]) for l in ca_lab])
pc = ca.coord[pmask].mean(axis=0)
nt, ct = ca.coord[0], ca.coord[-1]
pw = ca.coord[pmask]
geom = {
    "reference_chain": ref,
    "n_term_residue": ca_lab[0], "c_term_residue": ca_lab[-1],
    "n_term_to_paratope_centroid_A": round(float(np.linalg.norm(nt - pc)), 2),
    "c_term_to_paratope_centroid_A": round(float(np.linalg.norm(ct - pc)), 2),
    "n_term_to_nearest_paratope_A": round(float(np.min(np.linalg.norm(pw - nt, axis=1))), 2),
    "c_term_to_nearest_paratope_A": round(float(np.min(np.linalg.norm(pw - ct, axis=1))), 2),
    "n_to_c_term_A": round(float(np.linalg.norm(nt - ct)), 2),
    "paratope_max_extent_A": round(float(np.max(np.linalg.norm(pw[:, None] - pw[None], axis=-1))), 2),
    "radius_of_gyration_A": round(float(struc.gyration_radius(va)), 2),
}
save_json("s02_interface.json",
          {"sequence": seq, "n_residues": len(seq), "cdrs": cdrs,
           "chain_pairs": results, "terminal_geometry": geom,
           "paratope_by_region": {k: {"residues": r, "n": len(r)} for k, r in region.items()}})

print(f"VHH {len(seq)} residues (Kabat-numbered with insertion codes)")
for k, c in cdrs.items(): print(f"  {k}: len {c['length']:2d}  {c['seq']}")
print("\nper-copy interface")
for v, r in results.items():
    print(f"  {v}:{r['toxin_chain']}  BSA {r['buried_surface_area_A2']:7.1f} A2  "
          f"interface {r['interface_area_A2']:6.1f} A2  paratope {r['n_paratope']:2d}  epitope {r['n_epitope']:2d}")
bsas = [r["buried_surface_area_A2"] for r in results.values()]
print(f"  BSA across 4 copies: mean {np.mean(bsas):.1f}  sd {np.std(bsas, ddof=1):.1f} A2")
print("\nparatope by region: " + ", ".join(f"{k} {len(r)}" for k, r in region.items()))
print("paratope: " + " ".join(results[ref]["paratope_labels"]))
print("\nterminal geometry (A)")
for k, x in geom.items():
    if k.endswith("_A"): print(f"  {k[:-2]}: {x}")
