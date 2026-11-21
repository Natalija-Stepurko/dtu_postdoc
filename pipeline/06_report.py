#!/usr/bin/env python
"""Stage 6. Render docs/index.html from the stage outputs."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
from nbswitch import report as R
from nbswitch.config import (DOCS, PDB_ID, REPO_URL, NOTEBOOK_URL, load_json,
                             REFERENCE_VHH_CHAIN, TARGET_ID, TARGET, TARGETS,
                             results_dir)

ref = load_json("01_reference.json")
tether = load_json("02_tether.json")
lib = load_json("03_constructs.json")
try:
    fold = load_json("04_folding.json")
except FileNotFoundError:
    fold = None
for _name, _key in (("07_redesign.json", "redesign"),
                    ("08_redesign_folding.json", "redesign_fold"),
                    ("09_complexes.json", "complexes")):
    try:
        globals()[_key] = load_json(_name)
    except FileNotFoundError:
        globals()[_key] = None

meta, geom = ref["metadata"], ref["terminal_geometry"]
iface = ref["chain_pairs"][REFERENCE_VHH_CHAIN]
region = {k: v["n"] for k, v in ref["paratope_by_region"].items()}
repro = ref["reproducibility"]
bN, bC = tether["best_per_terminus"]["N"], tether["best_per_terminus"]["C"]
closed_C = max(r["linker_residues"] for r in tether["scan"]
               if r["terminus"] == "C" and not r["window_open"])
model = tether["model"]

region_summary = ", ".join(
    f"{v['n']} in {k}" for k, v in ref["paratope_by_region"].items() if v["n"])
_n_in_paratope = geom["n_term_residue"] in set(iface["paratope_labels"])
n_term_verdict = (
    "The N-terminal residue is itself part of the paratope, which makes masking easy and "
    "measuring it impossible — see the set-aside case above."
    if _n_in_paratope else
    "Because that gap exists, a peptide only reaches the paratope if its linker is long enough "
    "and it folds back, which is what makes the controls able to fail.")

schematic_facts = dict(
    n_paratope=iface["n_paratope"], cdr3_frac=str(region["CDR3"]),
    span_n=geom["n_term_to_paratope_centroid_A"],
    span_c=geom["c_term_to_paratope_centroid_A"],
    win_n=f"{bN['window_decades']:.1f}", closed_c=closed_C,
    best_linker=bN["linker_residues"],
    km_ceiling_uM=bN["km_max_for_masking_M"] * 1e6,
    n_masks=len(lib["masks"]), n_linkers=len(lib["linker_lengths"]),
    n_folded=(fold["n_folded"] if fold else lib["n_constructs"]),
    n_redesigns=(redesign["n_unique"] if redesign else 0),
    boltz_tag=("RUN" if complexes else "RUNNING"))

# the complex the pipeline was run on first and then set aside
CAUTIONARY = "9GCN"
cautionary_html = ""
if TARGET_ID != CAUTIONARY:
    import json as _json
    _d = results_dir(CAUTIONARY)
    try:
        _ref = _json.loads((_d / "01_reference.json").read_text())
        _teth = _json.loads((_d / "02_tether.json").read_text())
        _fold = _json.loads((_d / "04_folding.json").read_text())
        try:
            _cplx = _json.loads((_d / "09_complexes.json").read_text())
        except FileNotFoundError:
            _cplx = None
        cautionary_html = R.cautionary_section(
            _ref, _teth, _fold, _cplx,
            TARGETS[CAUTIONARY]["label"], CAUTIONARY,
            TARGETS[CAUTIONARY]["reference_vhh_chain"])
    except FileNotFoundError:
        cautionary_html = ""

background_html = R.background_section(
    bN["km_min_for_release_M"] * 1e9, bN["km_max_for_masking_M"] * 1e6)

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
    import statistics as _st
    graft = kinds.get("epitope graft", [])
    inert = kinds.get("inert control", [])
    scram = kinds.get("scrambled control", [])
    separation = (np.mean(graft) - np.mean(inert)) if graft and inert else float("nan")
    discriminates = separation > 2 * (np.std(graft, ddof=1) if len(graft) > 1 else 0)
    verdict = ("the designed masks separate from the controls"
               if discriminates else
               "the designed masks do <b>not</b> separate from the controls")
    fold_cards = "".join([
        card(f"{np.mean(graft):+.0f}%", "epitope grafts", "mean occlusion"),
        card(f"{np.mean(scram):+.0f}%", "scrambled", "same residues, shuffled"),
        card(f"{np.mean(inert):+.0f}%", "poly-glycine", "no designed interaction", hl=True),
        card(f"{separation:+.1f} pp", "graft minus inert", "the quantity that matters"),
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
<div class="chart wide">{{R.occlusion_schematic()}}</div>
<div class="cards">{fold_cards}</div>
<div class="charts2">
{chart("occlusion_by_linker", "Paratope occlusion against linker length, by mask.")}
{chart("occlusion_vs_confidence", "Occlusion against fold confidence in the nanobody region.")}
</div>
<div class="scroll"><table><thead><tr>
<th>construct</th><th>kind</th><th>linker</th><th>occluded %</th><th>min CA–CA Å</th>
<th>contacts &lt;8 Å</th><th>VHH pLDDT</th><th>core RMSD Å</th>
</tr></thead><tbody>{table_rows}</tbody></table></div>
<div class="banner" style="border-left-color:#C06014"><b>The result is negative, and that is the
point of running the controls.</b> {{verdict}}: the grafts occlude {{np.mean(graft):+.1f}}% of the
paratope on average, the scrambled masks {{np.mean(scram):+.1f}}% and the inert poly-glycine mask
{{np.mean(inert):+.1f}}% — a separation of {{separation:+.1f}} percentage points, with the inert
control <b>above</b> the designed one. Any peptide tethered at this N-terminus lands near the
paratope, so this measurement is reading <b>tether geometry, not design</b>.</div>
<div class="banner"><b>What follows from that.</b> Occlusion measured this way cannot rank masks,
and the module needs a discriminating readout before it can feed a ranking step — the two-chain
competition below is the first candidate. Had the controls been left out, the grafts'
{{np.mean(graft):+.1f}}% would have looked like a working design.</div>
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

# ---- complex-prediction section ----
if complexes:
    ok = [r for r in complexes["results"] if r["status"] == "ok"]
    ctrl = next((r for r in ok if r["case"] == "control_free_vhh"), None)
    crows = "".join(
        f"<tr class=\"{'hl' if r['case']=='control_free_vhh' else ''}\"><td>{r['case']}</td>"
        f"<td>{r['note']}</td><td>{r['atom_contacts']}</td>"
        f"<td>{r['buried_surface_area_A2']:.0f}</td>"
        f"<td>{r['fraction_of_contacts_on_paratope']:.2f}</td>"
        f"<td>{r['seconds']/60:.0f}</td></tr>" for r in ok)
    complex_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Module 02b · the other half</div>
<h3>With the mask on, can the antigen still bind?</h3>
<p>Occlusion says the mask covers the site. It does not say the antigen is excluded — those are
different claims, and only a two-chain prediction separates them. The free nanobody with its toxin
is the calibration case: the crystal structure already answers it, so how well the prediction
recovers that answer sets how much the masked cases are worth.</p></div>
<div class="scroll"><table><thead><tr><th>case</th><th>what it tests</th><th>atom contacts</th>
<th>buried Å²</th><th>fraction on paratope</th><th>min</th></tr></thead>
<tbody>{{crows}}</tbody></table></div>
<div class="banner"><b>How to read the control.</b> Predicted in single-sequence mode, with no
alignment server, so designed sequences never leave the machine — accuracy was traded for that
deliberately. If the free nanobody does not reproduce a toxin pose on its known paratope, the
masked numbers describe the model rather than the molecule.</div>
</section>"""
else:
    complex_html = ""

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
<title>Designing a binder that is switched off until it arrives</title>
<meta name="description" content="A design method for conditionally active binders, and a worked
example on an experimental nanobody-toxin complex.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;650&display=swap" rel="stylesheet">
<style>{R.CSS}</style></head><body><div class="page">

<header>
<div class="kicker"><span>Concept demonstration</span><b>Method + worked example</b></div>
<h1>Designing a binder that is switched off until it arrives</h1>
<p class="lede">How do you design a protein that does nothing in circulation and starts working only
where you want it to? This page sets out the question, proposes a design method, and then runs that
method end to end on one experimental structure to see what it produces.</p>
<div class="banner"><b>Read this as a method, not a molecule.</b> The worked example exists to
exercise the pipeline on real coordinates. Its candidates are deliberately simple and none is
proposed for synthesis.</div>
</header>

<!-- ============================ PART 1 ============================ -->
<section class="sec">
<div class="stack"><div class="eyebrow">Part 1 · the question</div>
<h2>A binder that is inert until something local switches it on</h2>
<p>Most therapeutic binders work everywhere in the body at once, which is why their side effects so
often come from hitting the right target in the wrong place. A conditionally active binder is meant
to be inert while it circulates and active only where some local cue — a protease, an acidic
pocket, a metal ion — tells it to switch on.</p>
<p>The usual way to build one is to <b>tether a mask over the binding site</b> and give the tether a
cleavable linker. Clinically this architecture has a long history and a mixed record; the biology is
not the hard part. The hard part is a thermodynamic one, and it is easy to state.</p></div>
<div class="chart wide">{{R.states_schematic()}}</div>
<div class="banner"><b>The tension.</b> The tether holds the mask at a high effective concentration
near the site, so even a weak mask blocks well. Once the protease cuts, that advantage disappears —
and a mask chosen for how firmly it blocks may simply stay put and keep blocking. <b>A mask can be
too good.</b></div>
<p>So the design target is not affinity. It is a <em>window</em>: tight enough to hold state A,
loose enough to lose state B. That window may be wide, narrow, or empty, and which of those it is
depends on geometry that can be measured before any candidate exists. That observation is what the
method below is built around.</p>
<div class="cards">
{{card("2", "states to satisfy", "masked and released, at once")}}
{{card("opposed", "the two requirements", "one bounds K<sub>D</sub> above, one below")}}
{{card("geometry", "what decides the window", "measurable before any design")}}
</div>
</section>

<!-- ============================ PART 2 ============================ -->
<section class="sec">
<div class="stack"><div class="eyebrow">Part 2 · the method</div>
<h2>A pipeline built around that window</h2>
<p>The ordering is the design decision. Generative models are usually put first, and the structure
used afterwards to filter what they produce. Here it runs the other way: <b>the complex is measured
and the tether priced before anything is generated</b>, so the attachment point, the linker length
and the usable range of mask affinities are already fixed when generation starts. Generation then
searches a space the geometry has narrowed, rather than a space that has to be filtered down
afterwards.</p></div>
<div class="chart wide">{{R.architecture_schematic()}}</div>

<h3>Why each module is separate</h3>
<div class="scroll"><table><thead><tr><th>module</th><th>what it answers</th>
<th>why it is its own box</th></tr></thead><tbody>
<tr><td><b>00 Constrain</b></td><td>Where could a mask attach, and what mask affinity could
possibly work?</td><td>It runs once per target and bounds everything downstream. Skipping it means
generating candidates that no affinity could rescue.</td></tr>
<tr><td><b>01 Generate</b></td><td>What masks and linkers are worth trying?</td>
<td>Grafts, language models and diffusion models are interchangeable here. Nothing downstream needs
to know which was used.</td></tr>
<tr><td><b>02a Trigger</b></td><td>Does the cue open the switch here and nowhere else?</td>
<td>Selectivity, not cleavability, is the quantity. The cue itself is a parameter — protease, pH or
metal ion all plug into the same slot.</td></tr>
<tr><td><b>02b Binding</b></td><td>Is the site covered while tethered, and free once cut?</td>
<td>Two states, scored separately. Collapsing them to one number is what hides the failure mode.</td></tr>
<tr><td><b>03 Rank</b></td><td>Which designs survive all the constraints at once?</td>
<td>Trade-offs are kept visible. Averaging the scores would let a design that cannot release
compensate with a high masking score.</td></tr>
<tr><td><b>04 Measure</b></td><td>What is actually true?</td><td>The only module producing data the
others can learn from, and the slowest.</td></tr>
</tbody></table></div>

<div class="banner"><b>What "universal" can honestly mean.</b> A single mask sequence that works on
every binder is unlikely: paratopes differ. What can transfer is the <em>architecture</em> — a fixed
attachment geometry, a reusable mask family, an interchangeable trigger linker and one scoring
protocol. That is a cassette, and it is a realistic target.</div>

<div class="banner" style="border-left-color:#C06014"><b>The rate limit is the lab, not the
compute.</b> Every module except 04 runs in minutes to hours. Module 04 returns a handful of
measurements a month, and conventional supervised learning wants 20–30 results before it can even
validate, around 100 to optimise against, and hundreds more to train on. At three measurements a
month that is years, which is why the computational modules are posed as a way to <em>spend</em>
those measurements well, not as a substitute for them.</div>
</section>

<!-- ============================ PART 3 ============================ -->
<section class="sec">
<div class="stack"><div class="eyebrow">Part 3 · the method in action</div>
<h2>Running it on a real complex</h2>
<p>The example below applies every built module to one deposited structure: {meta['title']},
{meta['method'].lower()} at {meta['resolution_A']} Å, deposited {meta['deposition_date']}. Four
independent copies of the complex sit in the asymmetric unit, so the interface is measured four
times and the spread between them is a free reproducibility check.</p></div>
<div class="chart wide">{R.pipeline_schematic(schematic_facts)}</div>
<p>This complex was not the first one tried. The constraint module is also a screening step: it
reports, before anything is generated, whether a target can support a measurable mask at all. The
first candidate failed that test, and the record of it is kept below.</p>
{cautionary_html}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Module 00 · constrain</div>
<h3>Measuring the complex</h3>
<p>The paratope is defined as the residues that lose solvent-accessible surface when the complex
forms, and the geometry that matters for a tether is how far each terminus sits from that
surface.</p></div>
<div class="cards">
{card(f"{{iface['interface_area_A2']:.0f}} Å²", "interface area", f"{{repro['bsa_mean_A2']}} ± {{repro['bsa_sd_A2']}} Å² buried across {{repro['n_copies']}} copies")}
{card(iface['n_paratope'], "paratope residues", "burying &gt;1 Å² on binding")}
{card(f"{{region['CDR3']}}/{{iface['n_paratope']}}", "of it in CDR3", f"CDR1 {{region['CDR1']}} · CDR2 {{region['CDR2']}} · framework {{region['framework']}}", hl=True)}
{card(ref['cdrs']['CDR3']['length'], "CDR3 length", "residues")}
</div>
{chart_wide("paratope_footprint", "Surface buried per nanobody residue on complex formation.")}
<div class="banner"><b>What this fixes before anything is designed.</b> The binding site spans
{region_summary}. The N-terminus sits <b>{geom['n_term_to_nearest_paratope_A']} Å</b> from the
nearest paratope residue and the C-terminus <b>{geom['c_term_to_nearest_paratope_A']} Å</b> — so a
mask tethered at the N-terminus has a real gap to span, and the two termini are not
interchangeable. {n_term_verdict}</div>
<div class="scroll"><table><thead><tr><th>geometry</th><th>Å</th><th>consequence</th></tr></thead>
<tbody>
<tr class="hl"><td>N-terminus to paratope centroid</td><td>{geom['n_term_to_paratope_centroid_A']}</td><td>a short linker already reaches</td></tr>
<tr><td>N-terminus to nearest paratope residue</td><td>{geom['n_term_to_nearest_paratope_A']}</td><td>the N-terminus is in the paratope</td></tr>
<tr><td>C-terminus to paratope centroid</td><td>{geom['c_term_to_paratope_centroid_A']}</td><td>{geom['c_term_to_paratope_centroid_A']/geom['n_term_to_paratope_centroid_A']:.1f}× further to span</td></tr>
<tr><td>paratope maximum extent</td><td>{geom['paratope_max_extent_A']}</td><td>surface a mask must cover</td></tr>
</tbody></table></div>
<div class="viewer">
<div class="vbar"><span class="eyebrow">The deposited complex</span>
<button data-spin="solved">spin</button></div>
<div class="vstage" id="stage-solved"></div>
<div class="vkey"><span><i style="background:#2D5BD1"></i>nanobody</span>
<span><i style="background:#C06014"></i>α-cobratoxin</span>
<span><i style="background:#C06014"></i>paratope sticks ({iface['n_paratope']} residues)</span>
<span><b style="color:#12191F">N</b> / <b style="color:#6E7880">C</b> termini · primes mark the toxin</span>
<span class="dim">sticks coloured by element</span></div>
</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Module 00 · constrain</div>
<h3>Pricing the tether</h3>
<p>A tethered mask behaves as though present at an effective concentration set by the linker.
Treating the Gly/Ser linker as a worm-like chain with persistence length
{model['persistence_length_A']} Å and {model['rise_per_residue_A']} Å per residue gives
C<sub>eff</sub> at the span measured above — a distance taken from the structure, not assumed.</p>
</div>
<div class="eq">⟨r²⟩ = 2·l<sub>p</sub>·L<sub>c</sub>·[1 − (l<sub>p</sub>/L<sub>c</sub>)(1 − e<sup>−L<sub>c</sub>/l<sub>p</sub></sup>)]
<br>C<sub>eff</sub>(d) = 1660.54 · (3 / 2π⟨r²⟩)<sup>3/2</sup> · e<sup>−3d²/2⟨r²⟩</sup> &nbsp;[M]
<br><br>P<sub>masked</sub> = C<sub>eff</sub> / (C<sub>eff</sub> + K<sub>D,M</sub>)
<br>P<sub>NT</sub> = (T/K<sub>D,T</sub>) / (1 + T/K<sub>D,T</sub> + M/K<sub>D,M</sub>)</div>
<div class="charts2">
{chart("effective_concentration", "Effective concentration against linker length, for both attachment points.")}
{chart("tradeoff", "The two requirements as functions of mask affinity.")}
</div>
{chart_wide("design_window", "Usable mask affinity against linker length. Shaded where both requirements hold at once.")}
<div class="cards">
{card(f"{{bN['window_decades']:.1f}}", "decades of usable K<sub>D</sub>", f"N-terminal, {{bN['linker_residues']}}-residue linker", hl=True)}
{card(f"{{bC['window_decades']:.1f}}", "decades, C-terminal", f"best at {{bC['linker_residues']}} residues")}
{card(f"≤{{closed_C}}", "residues where C-terminal fails", "no window at any affinity")}
{card(f"{{bN['km_max_for_masking_M']*1e6:.0f}} µM", "loosest mask that still blocks", f"N-terminal, {{bN['linker_residues']}} residues")}
</div>
<div class="banner"><b>What the constraint module returns.</b> Attach at the N-terminus, use a
{bN['linker_residues']}-residue linker, and aim for a mask K<sub>D</sub> somewhere below
{bN['km_max_for_masking_M']*1e6:.0f} µM but above {bN['km_min_for_release_M']*1e9:.0f} nM. The window
is <b>{bN['window_decades']:.1f} decades</b> wide at the N-terminus and
<b>{bC['window_decades']:.1f}</b> at the C-terminus, where it closes entirely below
{closed_C+1} residues. A <b>micromolar</b> mask is sufficient — about
{bN['km_max_for_masking_M']/bN['km_min_for_release_M']:.0f}× looser than the release floor — so
affinity maturation would be the wrong move.</div>
{chart("release_sensitivity", "The release floor moves with toxin affinity and concentration; the masking ceiling does not.")}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Module 01 · generate</div>
<h3>Candidates, and the controls that make them legible</h3>
<p>One peptide is already known to engage this paratope: the antigen's own contact segment. Grafting
it needs no generative model and brings its own negative controls.</p></div>
<div class="charts2">
<div class="scroll"><table><thead><tr><th>antigen segment</th><th>sequence</th><th>length</th>
<th>contacts</th></tr></thead><tbody>{seg_rows}</tbody></table></div>
<div class="scroll"><table><thead><tr><th>mask</th><th>sequence</th><th>aa</th><th>kind</th>
<th>origin</th></tr></thead><tbody>{masks_rows}</tbody></table></div>
</div>
<div class="banner"><b>On the candidates.</b> <span class="mono">epi1</span> is the tip of
α-cobratoxin loop II, the region carrying the residues that drive receptor binding — this nanobody
neutralizes by covering the toxin's functional site, so that same peptide is the obvious mask.
<span class="mono">epi2</span> carries a free cysteine and is kept as a deliberate counter-example.
<span class="mono">polyG</span> has no designed interaction and should fail.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Module 02b · score the binding</div>
<h3>Is the site covered?</h3></div>
<div class="chart wide">{{R.occlusion_schematic()}}</div>
{folding_html}
<div class="viewer">
<div class="vbar"><select id="pick">{options}</select>
<button data-spin="designed">spin</button></div>
<div class="vstage" id="stage-designed"></div>
<div class="vkey"><span><i style="background:#2D5BD1"></i>nanobody</span>
<span><i style="background:#0E7C7B"></i>mask</span>
<span><i style="background:#8A9BA8"></i>linker</span>
<span><i style="background:#C06014"></i>paratope</span>
<span><b style="color:#12191F">N</b> / <b style="color:#6E7880">C</b> termini</span>
<span class="dim">sticks coloured by element</span></div>
<p class="dim" id="dcap" style="font-size:.85rem"></p>
</div>
</section>

{redesign_html}

{complex_html}

<section class="sec">
<div class="stack"><div class="eyebrow">What the example does not settle</div>
<h2>The gaps, named</h2></div>
<div class="scroll"><table><thead><tr><th>module</th><th>what is missing</th>
<th>what it would take</th></tr></thead><tbody>
<tr><td><b>02a Trigger</b></td><td>no selectivity score at all</td><td><span class="tag gap">data</span>
cleavage measurements for the relevant proteases against a negative panel. Published profiling shows
the obvious candidates share a P1′ preference with host enzymes, so selectivity, not cleavability,
is the obstacle</td></tr>
<tr><td><b>01 Generate</b></td><td>no de novo backbones</td><td><span class="tag">GPU</span>
RFdiffusion or BindCraft would widen the search beyond grafted and redesigned peptides</td></tr>
<tr><td><b>02b Binding</b></td><td>one structure per design, no ensemble</td>
<td><span class="tag">GPU</span> BioEmu-2 would sample conformations instead of a single pose</td></tr>
<tr><td><b>03 Rank</b></td><td>scores computed, selection not automated</td>
<td>a Pareto front with uncertainty attached, once 02a exists to rank against</td></tr>
<tr><td><b>04 Measure</b></td><td>nothing has been made</td><td>the loop is open until it is</td></tr>
</tbody></table></div>
<div class="banner"><b>The honest summary.</b> The constraint and binding modules are real and
produce numbers. The trigger module is empty, and it is empty for want of data rather than for want
of a model — which is the kind of gap that no amount of compute closes.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Provenance</div>
<h2>Where the numbers come from</h2>
<p>Every value and chart here is produced by the analysis code, run in stages against one deposited
structure. Each stage writes its results to a file and this page is generated from those files, so
nothing is typed in by hand. All of it runs on CPU.</p></div>
<div class="scroll"><table><thead><tr><th>stage</th><th>method</th><th>produces</th></tr></thead>
<tbody>
<tr><td>Measure the complex</td><td>Shrake-Rupley surface burial</td><td>paratope, CDRs, terminal geometry</td></tr>
<tr><td>Price the tether</td><td>worm-like chain</td><td>effective concentration, masking ceiling, release floor</td></tr>
<tr><td>Build the library</td><td>antigen epitope graft, with controls</td><td>mask × linker constructs</td></tr>
<tr><td>Fold and measure</td><td>ESMFold</td><td>paratope occlusion per construct</td></tr>
<tr><td>Redesign the mask</td><td>ProteinMPNN</td><td>new masks on the folded backbone</td></tr>
<tr><td>Test the competition</td><td>Boltz-2</td><td>two-chain complex with the antigen</td></tr>
</tbody></table></div>
<p class="dim" style="font-size:.88rem">Code and stage outputs:
<a href="{REPO_URL}">{REPO_URL.split("//")[1]}</a></p>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">Prior work</div>
<h2>Where this problem is familiar from</h2>
<p>The two-state, tethered-affinity problem on this page is one I characterised experimentally
during the PhD, in a protein that solves it naturally.</p></div>
{background_html}
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
<li>Dauparas J. <i>et al.</i> Robust deep learning-based protein sequence design using ProteinMPNN.
<i>Science</i> 378:49–56 (2022).
<a href="https://doi.org/10.1126/science.add2187">10.1126/science.add2187</a></li>
<li>Wohlwend J. <i>et al.</i> Boltz-1 and Boltz-2, open biomolecular interaction models.
<a href="https://github.com/jwohlwend/boltz">github.com/jwohlwend/boltz</a></li>
<li>Shrake A., Rupley J.A. Environment and exposure to solvent of protein atoms.
<i>Journal of Molecular Biology</i> 79:351–371 (1973).
<a href="https://doi.org/10.1016/0022-2836(73)90011-9">10.1016/0022-2836(73)90011-9</a></li>
<li>Martin-Alonso C. <i>et al.</i> Deep learning guided design of protease substrates.
<i>Nature Communications</i> 17:54 (2026).
<a href="https://doi.org/10.1038/s41467-025-67226-1">10.1038/s41467-025-67226-1</a></li>
</ul>
</section>

<footer>Figures and values generated by the analysis code. Structural data from the RCSB Protein
Data Bank.</footer>
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
