#!/usr/bin/env python
"""Stage 6. Render docs/index.html from the stage outputs."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import report as R
from nbswitch.config import (DOCS, PDB_ID, REPO_URL, NOTEBOOK_URL, load_json,
                             REFERENCE_VHH_CHAIN)

ref = load_json("01_reference.json")
tether = load_json("02_tether.json")
lib = load_json("03_constructs.json")
try:
    fold = load_json("04_folding.json")
except FileNotFoundError:
    fold = None
try:
    redesign = load_json("07_redesign.json")
    redesign_fold = load_json("08_redesign_folding.json")
except FileNotFoundError:
    redesign = redesign_fold = None

meta, geom = ref["metadata"], ref["terminal_geometry"]
iface = ref["chain_pairs"][REFERENCE_VHH_CHAIN]
region = {k: v["n"] for k, v in ref["paratope_by_region"].items()}
repro = ref["reproducibility"]
bN, bC = tether["best_per_terminus"]["N"], tether["best_per_terminus"]["C"]
closed_C = max(r["linker_residues"] for r in tether["scan"]
               if r["terminus"] == "C" and not r["window_open"])
model = tether["model"]

R.copy_assets(lib["constructs"])
R.write_reference_pdb()

# paratope residue numbers for the viewer (3Dmol reads author numbering)
para_resi = sorted({int("".join(ch for ch in l[1:] if ch.isdigit()))
                    for l in iface["paratope_labels"]})

folded_ids = sorted(p.stem for p in (DOCS / "structures").glob("*.pdb")
                    if p.stem != f"{PDB_ID}_AX")
by_id = {c["id"]: c for c in lib["constructs"]}
models = {}
for cid in folded_ids:
    c = by_id.get(cid)
    if not c:
        continue
    entry = {"vhh": f"{c['vhh_range'][0]}-{c['vhh_range'][1]}",
             "paratope": [int(c["vhh_range"][0] - 1 + i + 1) for i in range(0)],
             "caption": ""}
    off = c["vhh_range"][0] - 1
    idx = [i for i, l in enumerate(lib["vhh_residue_labels"])
           if l in set(lib["paratope_labels"])]
    entry["paratope"] = [off + i + 1 for i in idx]
    if c["mask_range"]:
        entry["mask"] = f"{c['mask_range'][0]}-{c['mask_range'][1]}"
        entry["linker"] = f"{c['linker_range'][0]}-{c['linker_range'][1]}"
        row = next((r for r in fold["results"] if r["id"] == cid), None) if fold else None
        extra = (f" · paratope {row['occlusion']*100:+.1f}% occluded"
                 f" · VHH pLDDT {row['plddt_vhh']:.0f}") if row else ""
        entry["caption"] = (f"{cid} — {c['mask_kind']}, mask {c['mask_seq']}, "
                            f"{c['linker_residues']}-residue linker{extra}")
    else:
        entry["caption"] = "free nanobody, folded as the reference state"
    models[cid] = entry

default_id = "epi1_L5" if "epi1_L5" in models else (folded_ids[0] if folded_ids else "")

def card(v, l, s="", hl=False): return R._card(v, l, s, hl)
def chart(n, c): return R._chart(n, c)
def chart_wide(n, c): return R._chart(n, c).replace('class="chart"', 'class="chart wide"')

# ---- folding section ----
if fold:
    rows = fold["results"]
    kinds = {}
    for r in rows:
        kinds.setdefault(r["kind"], []).append(r["occlusion"] * 100)
    graft = np.mean(kinds.get("epitope graft", [0]))
    inert = np.mean(kinds.get("inert control", [0]))
    scram = np.mean(kinds.get("scrambled control", [0]))
    best_row = max(rows, key=lambda r: r["occlusion"])
    rmsd = [r["vhh_core_rmsd_A"] for r in rows]
    fold_cards = "".join([
        card(f"{best_row['occlusion']*100:.0f}%", "best occlusion",
             f"{best_row['id']} · {best_row['kind']}", hl=True),
        card(f"{graft:.0f}%", "epitope grafts", "mean across linkers"),
        card(f"{scram:.0f}%", "scrambled control", "mean"),
        card(f"{inert:.0f}%", "inert control", "mean"),
        card(f"{np.mean(rmsd):.2f} Å", "VHH core RMSD", "fusion against free nanobody"),
    ])
    table_rows = "".join(
        f"<tr class=\"{'hl' if r['kind']=='epitope graft' else 'ctrl' if 'control' in r['kind'] else ''}\">"
        f"<td>{r['id']}</td><td>{r['kind']}</td><td>{r['linker']}</td>"
        f"<td>{r['occlusion']*100:+.1f}</td><td>{r.get('min_mask_paratope_CA_A','—')}</td>"
        f"<td>{r.get('n_mask_paratope_contacts_8A','—')}</td>"
        f"<td>{r['plddt_vhh']:.0f}</td><td>{r['vhh_core_rmsd_A']:.2f}</td></tr>"
        for r in sorted(rows, key=lambda r: -r["occlusion"]))
    folding_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Stage 4 · folded on CPU</div>
<h2>Does the mask actually reach the paratope?</h2>
<p>Every construct folded with ESMFold, then the paratope's solvent-accessible surface in the
fusion compared against the same residues in the folded free nanobody. The grafts are read against
their scrambled, truncated and inert controls at matched linker length; the controls are the
reason the comparison means anything.</p></div>
<div class="cards">{fold_cards}</div>
<div class="charts2">
{chart("occlusion_by_linker", "Paratope occlusion against linker length, by mask.")}
{chart("occlusion_vs_confidence", "Occlusion against fold confidence in the nanobody region.")}
</div>
<div class="scroll"><table><thead><tr>
<th>construct</th><th>kind</th><th>linker</th><th>occluded %</th><th>min CA–CA Å</th>
<th>contacts &lt;8 Å</th><th>VHH pLDDT</th><th>core RMSD Å</th>
</tr></thead><tbody>{table_rows}</tbody></table></div>
<div class="banner"><b>What a predicted structure is not.</b> ESMFold returns one structure per
sequence. A mask resting against the paratope is a hypothesis about the masked state, not a masking
efficiency and not a dissociation constant. pLDDT is a confidence score, not a folding free energy.
The informative quantity here is the contrast between graft and control, not any single number.</div>
</section>"""
else:
    folding_html = ""

# ---- redesign section ----
if redesign_fold and fold:
    rrows = redesign_fold["results"]
    grafts = [r["occlusion"] * 100 for r in fold["results"] if r["kind"] == "epitope graft"]
    rocc = [r["occlusion"] * 100 for r in rrows]
    rbest = max(rrows, key=lambda r: r["occlusion"])
    improved = sum(1 for c in redesign["constructs"]
                   if c["native_score"] and c["mpnn_score"] < c["native_score"])
    xs = np.array([r["mpnn_score"] for r in rrows])
    ys = np.array([r["occlusion"] * 100 for r in rrows])
    corr = float(np.corrcoef(xs, ys)[0, 1]) if len(xs) > 2 else float("nan")
    rtable = "".join(
        f"<tr><td>{r['mask_seq']}</td><td>{r['parent']}</td><td>{r['mpnn_score']:.3f}</td>"
        f"<td>{r['occlusion']*100:+.1f}</td><td>{r.get('min_mask_paratope_CA_A','—')}</td>"
        f"<td>{r['plddt_vhh']:.0f}</td><td>{r['vhh_core_rmsd_A']:.2f}</td></tr>"
        for r in sorted(rrows, key=lambda r: -r["occlusion"])[:12])
    redesign_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Stage 5 · closing the loop</div>
<h2>Redesigning the mask on its own backbone</h2>
<p>The grafted mask is borrowed from the antigen, so nothing about its sequence is tuned to this
paratope. Taking each folded graft as a backbone and letting ProteinMPNN rewrite only the mask
positions — nanobody and linker held fixed — turns the graft into a starting point. The redesigned
constructs are then folded and scored exactly as before, so the two sets are directly comparable.</p>
</div>
<div class="cards">
{{card(redesign['n_unique'], "distinct masks designed", f"from {{redesign['n_designs']}} samples across {{len(redesign['temperatures'])}} temperatures")}}
{{card(f"{{improved}}/{{redesign['n_unique']}}", "score better than parent", "by ProteinMPNN's own score")}}
{{card(f"{{rbest['occlusion']*100:+.0f}}%", "best occlusion", f"mask {{rbest['mask_seq']}}", hl=True)}}
{{card(f"{{np.mean(rocc):+.0f}}%", "redesign mean", f"grafts {{np.mean(grafts):+.0f}}%")}}
{{card(f"{{corr:+.2f}}", "score vs occlusion", "Pearson r")}}
</div>
<div class="charts2">
{{chart("redesign_landscape", "Redesigns against the grafts and controls, on the same axes.")}}
{{chart("redesign_score_vs_occlusion", "The design model's own score against the quantity it was not optimising.")}}
</div>
<div class="scroll"><table><thead><tr><th>mask</th><th>parent</th><th>MPNN score</th>
<th>occluded %</th><th>min CA–CA Å</th><th>VHH pLDDT</th><th>core RMSD Å</th>
</tr></thead><tbody>{{rtable}}</tbody></table></div>
<div class="banner"><b>What the correlation is for.</b> ProteinMPNN scores sequence likelihood given
a backbone. It is not optimising paratope occlusion and has no way to know about it, so the
relationship between the two is worth measuring rather than assuming. On this set the Pearson
correlation is <b>{{corr:+.2f}}</b> across {{len(rrows)}} designs.</div>
</section>"""
else:
    redesign_html = ""

masks_rows = "".join(
    f"<tr class=\"{'ctrl' if 'control' in m['kind'] else ''}\"><td>{k}</td>"
    f"<td>{m['seq']}</td><td>{len(m['seq'])}</td><td>{m['kind']}</td><td>{m['note']}</td></tr>"
    for k, m in lib["masks"].items())
seg_rows = "".join(
    f"<tr><td>{s['start_label']}–{s['end_label']}</td><td>{s['seq']}</td>"
    f"<td>{s['length']}</td><td>{s['n_epitope']}</td></tr>"
    for s in lib["epitope_segments"])

options = "".join(f'<option value="{i}"{" selected" if i == default_id else ""}>{i}</option>'
                  for i in folded_ids)

HTML = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>A protease-activated nanobody mask</title>
<meta name="description" content="Concept demonstration: geometry and thermodynamics of a tethered
nanobody mask, computed from an experimental nanobody–toxin complex.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;650&display=swap" rel="stylesheet">
<style>{R.CSS}</style></head><body><div class="page">

<header>
<div class="kicker"><span>Concept demonstration</span><b>{PDB_ID}</b>
<span>{meta['resolution_A']} Å</span><span>{lib['n_constructs']} constructs folded</span></div>
<h1>A protease-activated nanobody mask</h1>
<p class="lede">A nanobody that stays inert in circulation and binds toxin only after a protease
cuts its tether. The regulatory part lives in a separate mask and linker, so an established
neutralizing nanobody carries over unchanged. Everything below is computed from one experimental
structure.</p>
<div class="banner"><b>What this is.</b> A worked example showing the pipeline running end to end
on real coordinates. The mask candidates are deliberately simple and nothing here is proposed for
synthesis. The point is the method and the numbers it produces, not the molecules.</div>
</header>

<section class="sec">
<div class="stack"><div class="eyebrow">The construct</div>
<h2>Three states, one chain</h2>
<p>A single polypeptide: designed mask <span class="mono">M</span>, protease-cleavable linker
<span class="mono">L</span>, fixed nanobody <span class="mono">N</span>. Cleavage severs the tether,
the mask diffuses away, and the original paratope is free to take toxin
<span class="mono">T</span>.</p></div>
<div class="eq">N–L–M &nbsp;&nbsp;<span class="dim">──protease──▶</span>&nbsp;&nbsp; N + M
&nbsp;&nbsp;&nbsp;<span class="dim">then</span>&nbsp;&nbsp;&nbsp; N + T &nbsp;⇌&nbsp; NT</div>
<div class="banner"><b>The tension that defines the problem.</b> The strongest mask is not the best
mask. One tight enough to block the paratope while tethered may stay bound after cleavage and keep
blocking it. The stages below turn that sentence into a bounded range.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Stage 1 · the reference complex</div>
<h2>Where the paratope is, and where a mask could attach</h2>
<p>{meta['title']}, {meta['method'].lower()} at {meta['resolution_A']} Å, deposited
{meta['deposition_date']}. Four independent copies sit in the asymmetric unit, so the interface is
measured four times and the spread is a free reproducibility check.</p></div>
<div class="cards">
{card(f"{iface['interface_area_A2']:.0f} Å²", "interface area", f"{repro['bsa_mean_A2']} ± {repro['bsa_sd_A2']} Å² buried, {repro['n_copies']} copies")}
{card(iface['n_paratope'], "paratope residues", "burying &gt;1 Å² on binding")}
{card(f"{region['CDR3']}/{iface['n_paratope']}", "paratope in CDR3", f"CDR1 {region['CDR1']} · CDR2 {region['CDR2']} · framework {region['framework']}", hl=True)}
{card(ref['cdrs']['CDR3']['length'], "CDR3 length", "residues")}
</div>
{chart_wide("paratope_footprint", "Surface buried per nanobody residue on complex formation. The binding site is almost entirely CDR3.")}
<div class="banner"><b>Two findings that set up everything downstream.</b> The paratope is
concentrated in CDR3 — {region['CDR3']} of {iface['n_paratope']} residues, with CDR2 contributing
none. And the N-terminal residue <span class="mono">{geom['n_term_residue']}</span> is itself part
of the paratope, burying
{iface['paratope_delta_sasa'][geom['n_term_residue']]:.0f} Å². The two termini are not
interchangeable places to hang a mask.</div>
<div class="scroll"><table><thead><tr><th>geometry</th><th>Å</th><th>consequence</th></tr></thead>
<tbody>
<tr class="hl"><td>N-terminus to paratope centroid</td><td>{geom['n_term_to_paratope_centroid_A']}</td><td>a short linker already reaches</td></tr>
<tr><td>N-terminus to nearest paratope residue</td><td>{geom['n_term_to_nearest_paratope_A']}</td><td>the N-terminus is in the paratope</td></tr>
<tr><td>C-terminus to paratope centroid</td><td>{geom['c_term_to_paratope_centroid_A']}</td><td>{geom['c_term_to_paratope_centroid_A']/geom['n_term_to_paratope_centroid_A']:.1f}× further to span</td></tr>
<tr><td>paratope maximum extent</td><td>{geom['paratope_max_extent_A']}</td><td>surface a mask must cover</td></tr>
<tr><td>nanobody radius of gyration</td><td>{geom['radius_of_gyration_A']}</td><td>scale reference</td></tr>
</tbody></table></div>

<div class="viewer">
<div class="vbar"><span class="eyebrow">Solved structure</span>
<button data-spin="solved">spin</button></div>
<div class="vstage" id="stage-solved"></div>
<div class="vkey"><span><i style="background:#2D5BD1"></i>nanobody</span>
<span><i style="background:#C06014"></i>α-cobratoxin</span>
<span><i style="background:#0E7C7B"></i>paratope ({iface['n_paratope']} residues)</span></div>
<p class="dim" style="font-size:.85rem">{PDB_ID} chains {REFERENCE_VHH_CHAIN} and
{iface['toxin_chain']}, as deposited.</p>
</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Stage 2 · what the tether is worth</div>
<h2>Masking and release pull in opposite directions</h2>
<p>A tethered mask behaves as if present at an effective concentration set by the linker. Treating
the Gly/Ser linker as a worm-like chain with persistence length
{model['persistence_length_A']} Å and {model['rise_per_residue_A']} Å per residue, the end-to-end
distribution gives C<sub>eff</sub> at the distance the linker has to span — a distance measured in
stage 1, not assumed.</p></div>
<div class="eq">⟨r²⟩ = 2·l<sub>p</sub>·L<sub>c</sub>·[1 − (l<sub>p</sub>/L<sub>c</sub>)(1 − e<sup>−L<sub>c</sub>/l<sub>p</sub></sup>)]
<br>C<sub>eff</sub>(d) = 1660.54 · (3 / 2π⟨r²⟩)<sup>3/2</sup> · e<sup>−3d²/2⟨r²⟩</sup> &nbsp;&nbsp;[M]
<br><br>P<sub>masked</sub> = C<sub>eff</sub> / (C<sub>eff</sub> + K<sub>D,M</sub>)
<br>P<sub>NT</sub> = (T/K<sub>D,T</sub>) / (1 + T/K<sub>D,T</sub> + M/K<sub>D,M</sub>)</div>
<div class="charts2">
{chart("effective_concentration", "Effective concentration against linker length, for both attachment points.")}
{chart("tradeoff", "The two requirements as functions of mask affinity.")}
</div>
{chart_wide("design_window", "Usable mask affinity against linker length. Shaded where both requirements hold at once.")}
<div class="cards">
{card(f"{bN['window_decades']:.1f}", "decades of usable K<sub>D</sub>", f"N-terminal, {bN['linker_residues']}-residue linker", hl=True)}
{card(f"{bC['window_decades']:.1f}", "decades, C-terminal", f"best at {bC['linker_residues']} residues")}
{card(f"≤{closed_C}", "residues where C-terminal fails", "no window at any affinity")}
{card(f"{bN['km_max_for_masking_M']*1e6:.0f} µM", "loosest mask that still blocks", f"N-terminal, {bN['linker_residues']} residues")}
</div>
<div class="banner"><b>The quantitative result.</b> An N-terminal fusion with a short linker leaves
about <b>{bN['window_decades']:.1f} decades</b> of usable mask affinity; a C-terminal fusion leaves
<b>{bC['window_decades']:.1f}</b>, and below <b>{closed_C+1} residues it has no window at all</b> —
no mask affinity satisfies both requirements at once. The tether is worth so much that a
<b>micromolar</b> mask still blocks: at {bN['linker_residues']} residues the ceiling is
{bN['km_max_for_masking_M']*1e6:.0f} µM, about
{bN['km_max_for_masking_M']/bN['km_min_for_release_M']:.0f}× above the release floor. Masks do not
need affinity maturation, and tight masks are actively harmful.</div>
{chart("release_sensitivity", "The release floor moves with toxin affinity and concentration. The masking ceiling does not.")}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Stage 3 · mask candidates</div>
<h2>Taken from the structure, with controls</h2>
<p>A mask has to engage this paratope. One peptide is already known to: the antigen's own contact
segment. Grafting it needs no generative model, and the scrambled, truncated and inert versions give
the comparison its meaning.</p></div>
<div class="charts2">
<div class="scroll"><table><thead><tr><th>antigen segment</th><th>sequence</th><th>length</th>
<th>contacts</th></tr></thead><tbody>{seg_rows}</tbody></table></div>
<div class="scroll"><table><thead><tr><th>mask</th><th>sequence</th><th>aa</th><th>kind</th>
<th>origin</th></tr></thead><tbody>{masks_rows}</tbody></table></div>
</div>
<div class="banner"><b>On the candidates.</b> <span class="mono">epi1</span> is the tip of
α-cobratoxin loop II, the region carrying the residues that drive receptor binding — this nanobody
neutralizes by covering the toxin's functional site, and that same peptide is the obvious mask.
<span class="mono">epi2</span> carries a free cysteine and is kept as a deliberate counter-example:
it would scramble disulfides in a real graft. <span class="mono">polyG</span> has no designed
interaction and should fail.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Designed constructs</div>
<h2>Predicted structures</h2>
<p>Each construct folded from sequence alone. Pick one to inspect: the mask is tethered at the
N-terminus and the paratope is drawn as sticks, so whether the mask covers it is visible
directly.</p></div>
<div class="viewer">
<div class="vbar"><select id="pick">{options}</select>
<button data-spin="designed">spin</button></div>
<div class="vstage" id="stage-designed"></div>
<div class="vkey"><span><i style="background:#2D5BD1"></i>nanobody</span>
<span><i style="background:#0E7C7B"></i>mask</span>
<span><i style="background:#8A9BA8"></i>linker</span>
<span><i style="background:#C06014"></i>paratope</span></div>
<p class="dim" id="dcap" style="font-size:.85rem"></p>
</div>
</section>

{folding_html}
{redesign_html}

<section class="sec">
<div class="stack"><div class="eyebrow">Beyond the demonstration</div>
<h2>What is not computed here</h2></div>
<div class="scroll"><table><thead><tr><th>component</th><th>what it would add</th>
<th>what it needs</th></tr></thead><tbody>
<tr><td>RFdiffusion · BindCraft</td><td>mask backbones beyond grafted epitopes</td><td><span class="tag">GPU</span> BindCraft also wants AF2 weights and PyRosetta</td></tr>
<tr><td>AlphaFold 3 · Boltz-2 · Chai-1</td><td>the nanobody–toxin complex and mask–toxin competition</td><td><span class="tag">GPU</span> model access</td></tr>
<tr><td>BioEmu-2</td><td>conformational ensembles in place of one structure</td><td><span class="tag">GPU</span> multichain support is recent</td></tr>
<tr><td>CleaveNet</td><td>linkers generated against a target cleavage profile</td><td><span class="tag gap">data</span> trained on 18 human MMPs; venom metalloproteinases are out of distribution</td></tr>
<tr><td>Trigger selectivity</td><td>cleavage by the local protease ensemble, resistance elsewhere</td><td><span class="tag gap">data</span> peptide-cleavage measurements that do not exist yet</td></tr>
</tbody></table></div>
<div class="banner"><b>The honest gap.</b> Trigger selectivity is not a modelling problem waiting on
better software. The substrate specificity data for venom proteases has to be generated first, and
published profiling shows a P1′ leucine preference shared with host MMPs — so telling venom protease
from host protease is the real obstacle, and it is wet-lab work.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Reproducing this</div>
<h2>Six stages, no hand-entered numbers</h2></div>
<div class="steps">
<div class="step"><div class="n">01</div><div><b>Reference structure</b>
<span>Fetch the complex, key residues on number and insertion code, measure the interface and the
terminal geometry.</span></div><span class="tag run">CPU</span></div>
<div class="step"><div class="n">02</div><div><b>Tether model</b>
<span>Worm-like-chain effective concentration, masking ceiling, release floor, design window.</span>
</div><span class="tag run">CPU</span></div>
<div class="step"><div class="n">03</div><div><b>Construct library</b>
<span>Epitope segments to mask candidates, with controls; enumerate mask × linker.</span></div>
<span class="tag run">CPU</span></div>
<div class="step"><div class="n">04</div><div><b>Folding</b>
<span>ESMFold on every construct, then paratope occlusion against the folded free nanobody.
Resumable.</span></div><span class="tag run">CPU or GPU</span></div>
<div class="step"><div class="n">05</div><div><b>Figures</b>
<span>Every chart on this page.</span></div><span class="tag run">CPU</span></div>
<div class="step"><div class="n">06</div><div><b>Report</b>
<span>This page, rendered from the stage outputs.</span></div><span class="tag run">CPU</span></div>
</div>
<div class="eq">git clone {REPO_URL}<br>uv venv &amp;&amp; uv pip install -e .<br>
for s in pipeline/0*.py; do python "$s"; done</div>
<p>Each stage writes a JSON file; figures and this page are generated from those files. The folding
stage runs on CPU — a GPU only makes it faster — and a
<a href="{NOTEBOOK_URL}">Colab notebook</a> runs the same step on a free T4.</p>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">References</div><h2>Sources</h2></div>
<ul class="ref">
<li>Damsbo A. <i>et al.</i> Co-crystal of VHH TPL1158_01_C09 with α-cobratoxin.
<a href="https://www.rcsb.org/structure/9GCN">PDB {PDB_ID}</a></li>
<li>Ahmadi S. <i>et al.</i> Nanobody-based recombinant antivenom for cobra, mamba and rinkhals
bites. <i>Nature</i> 647:716–725 (2025).
<a href="https://doi.org/10.1038/s41586-025-09661-0">10.1038/s41586-025-09661-0</a></li>
<li>Desnoyers L.R. <i>et al.</i> Tumor-specific activation of an EGFR-targeting probody.
<i>Science Translational Medicine</i> 5:207ra144 (2013).
<a href="https://doi.org/10.1126/scitranslmed.3006682">10.1126/scitranslmed.3006682</a></li>
<li>Zhou H.-X. Loops in proteins can be modeled as worm-like chains.
<i>Journal of Physical Chemistry B</i> 105:6763–6766 (2001).
<a href="https://doi.org/10.1021/jp011355n">10.1021/jp011355n</a></li>
<li>Lin Z. <i>et al.</i> Evolutionary-scale prediction of atomic-level protein structure with a
language model. <i>Science</i> 379:1123–1130 (2023).
<a href="https://doi.org/10.1126/science.ade2574">10.1126/science.ade2574</a></li>
<li>Shrake A., Rupley J.A. Environment and exposure to solvent of protein atoms.
<i>Journal of Molecular Biology</i> 79:351–371 (1973).
<a href="https://doi.org/10.1016/0022-2836(73)90011-9">10.1016/0022-2836(73)90011-9</a></li>
<li>Martin-Alonso C. <i>et al.</i> Deep learning guided design of protease substrates.
<i>Nature Communications</i> 17:54 (2026).
<a href="https://doi.org/10.1038/s41467-025-67226-1">10.1038/s41467-025-67226-1</a></li>
<li>Pacesa M. <i>et al.</i> One-shot design of functional protein binders with BindCraft.
<i>Nature</i> 646:483–492 (2025).
<a href="https://doi.org/10.1038/s41586-025-09429-6">10.1038/s41586-025-09429-6</a></li>
<li>Bennett N.R. <i>et al.</i> Atomically accurate de novo design of antibodies with RFdiffusion.
<i>Nature</i> 649:183–193 (2025).
<a href="https://doi.org/10.1038/s41586-025-09721-5">10.1038/s41586-025-09721-5</a></li>
</ul>
</section>

<footer>Figures and values generated by the analysis code in the
<a href="{REPO_URL}">repository</a>. Structural data from the RCSB Protein Data Bank.</footer>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
<script>{R.VIEWER_JS.replace("__MODELS__", json.dumps(models)).replace("__PARATOPE__", json.dumps(para_resi))}</script>
</body></html>"""

(DOCS / "index.html").write_text(HTML)
print(f"docs/index.html  {len(HTML)//1024} KB")
print(f"  {len(list((DOCS/'figures').glob('*.svg')))} figures, "
      f"{len(list((DOCS/'structures').glob('*.pdb')))} structures, "
      f"{len(models)} in the viewer")
