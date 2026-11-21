#!/usr/bin/env python
"""Stage 6. Render docs/index.html from the stage outputs."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import numpy as _np
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
for _name, _key in (("13_ranking.json", "ranking_data"),
                    ("12_feasibility.json", "feas_data"),
                    ("11_complexes_msa.json", "complexes_msa"),
                    ("10_trigger.json", "trigger"),
                    ("07_redesign.json", "redesign"),
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
<div class="stack"><div class="eyebrow">Step 04 · does it mask</div>
<h2>Does the mask actually reach the paratope?</h2>
<p>Every construct folded with ESMFold, then the paratope's solvent-accessible surface in the
fusion compared against the same residues in the folded free nanobody. The grafts are read against
their scrambled, truncated and inert controls at matched linker length; the controls are the
reason the comparison means anything.</p></div>
<div class="chart wide">{R.occlusion_schematic()}</div>
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
point of running the controls.</b> {verdict}: the grafts occlude {np.mean(graft):+.1f}% of the
paratope on average, the scrambled masks {np.mean(scram):+.1f}% and the inert poly-glycine mask
{np.mean(inert):+.1f}% — a separation of {separation:+.1f} percentage points, with the inert
control <b>above</b> the designed one. Any peptide tethered at this N-terminus lands near the
paratope, so this measurement is reading <b>tether geometry, not design</b>.</div>
<div class="banner"><b>What follows from that.</b> Occlusion measured this way cannot rank masks,
and the module needs a discriminating readout before it can feed a ranking step — the two-chain
competition below is the first candidate. Had the controls been left out, the grafts'
{np.mean(graft):+.1f}% would have looked like a working design.</div>
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
<div class="stack"><div class="eyebrow">Step 02 · generate, second method</div>
<h2>Redesigning the mask on its own backbone</h2>
<p>The grafted mask is borrowed from the antigen, so nothing about its sequence is tuned to this
paratope. Taking each folded graft as a backbone and letting ProteinMPNN rewrite only the mask
positions — nanobody and linker held fixed — turns the graft into a starting point. The redesigned
constructs are then folded and scored exactly as before, so the two sets are directly comparable.</p>
</div>
<div class="cards">
{card(redesign['n_unique'], "distinct masks designed", f"from {redesign['n_designs']} samples across {len(redesign['temperatures'])} temperatures")}
{card(f"{improved}/{redesign['n_unique']}", "score better than parent", "by ProteinMPNN's own score")}
{card(f"{rbest['occlusion']*100:+.0f}%", "best occlusion", f"mask {rbest['mask_seq']}", hl=True)}
{card(f"{np.mean(rocc):+.0f}%", "redesign mean", f"grafts {np.mean(grafts):+.0f}%")}
{card(f"{corr:+.2f}", "score vs occlusion", "Pearson r")}
</div>
<div class="charts2">
{chart("redesign_landscape", "Redesigns against the grafts and controls, on the same axes.")}
{chart("redesign_score_vs_occlusion", "The design model's own score against the quantity it was not optimising.")}
</div>
<div class="scroll"><table><thead><tr><th>mask</th><th>parent</th><th>MPNN score</th>
<th>occluded %</th><th>min CA–CA Å</th><th>VHH pLDDT</th><th>core RMSD Å</th>
</tr></thead><tbody>{rtable}</tbody></table></div>
<div class="banner"><b>What the correlation is for.</b> ProteinMPNN scores sequence likelihood given
a backbone. It is not optimising paratope occlusion and has no way to know about it, so the
relationship between the two is worth measuring rather than assuming. On this set the Pearson
correlation is <b>{corr:+.2f}</b> across {len(rrows)} designs.</div>
</section>"""
else:
    redesign_html = ""

# ---- ranking ----
if ranking_data and ranking_data.get("objectives"):
    _objs = ranking_data["objectives"]
    _front = ranking_data["front"]
    _rows = "".join(
        "<tr><td>" + r["id"] + "</td>" +
        "".join(f"<td>{r[o]:.3f}</td>" for o in _objs) + "</tr>"
        for r in _front[:10])
    _hdr = "".join(f"<th>{o}</th>" for o in _objs)
    ranking_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Step 05 · combine the three scores</div>
<h2>Keeping the three scores separate</h2>
<p>A design has to satisfy all three requirements. A construct that cannot be expressed is of no
use however well it masks, and one that cannot be cut is of no use however well it expresses. If
the three scores were averaged, a high score on one would hide a failure on another, so instead I
keep a design unless some other design is at least as good on every score and better on at least
one. No weights are applied.</p>
<p>Of {ranking_data['n_total']} scored designs, {ranking_data['n_front']} survive this test. The
rest are beaten outright by something else.</p></div>
<div class="scroll"><table class="wrap"><thead><tr><th>construct</th>{_hdr}</tr></thead>
<tbody>{_rows}</tbody></table></div>
{chart("pareto_panels", "The three scores plotted against each other. Designs on the front are "
       "marked; the rest are beaten on every axis by something else.")}
<div class="note"><b>The front shows the conflict rather than resolving it.</b> Designs scoring
highest on the trigger score sit at the shortest linker lengths, and those are the designs that
score lowest on masking. The pipeline does not choose between them, because choosing requires
knowing which failure is more expensive, and that is a question for measurement.</div>
<div class="note warn"><b>One of the three inputs is unreliable.</b> The masking axis here is the
occlusion measure, which I showed above does not separate designed masks from controls. The
combining stage works as a mechanism, and it will give sensible answers when that axis is replaced
by the tethered complex prediction. Until then the front should be read as a demonstration of the
method rather than as a list of candidates.</div>
</section>"""
else:
    ranking_html = ""

# ---- trigger accessibility ----
if trigger:
    import numpy as _np
    _by = {}
    for _r in trigger["results"]:
        _by.setdefault(_r["linker_residues"], []).append(_r)
    _win = {r["linker_residues"]: r["window_decades"]
            for r in tether["scan"] if r["terminus"] == "N"}
    _rows = "".join(
        f"<tr><td>{n} residues</td><td>{_win.get(n, float('nan')):.2f}</td>"
        f"<td>{_np.mean([x['fraction_exposed'] for x in g]):.2f}</td>"
        f"<td>{_np.mean([x['extension_ratio'] for x in g]):.2f}</td>"
        f"<td>{_np.mean([x['accessibility_score'] for x in g]):.2f}</td></tr>"
        for n, g in sorted(_by.items()))
    _short = min(_by)
    _bestwin = max(_win, key=_win.get)
    _acc_short = _np.mean([x["accessibility_score"] for x in _by[_short]])
    _acc_long = _np.mean([x["accessibility_score"] for x in _by[max(_by)]])
    trigger_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Step 03 · can it be cut</div>
<h2>Can the protease reach the linker at all?</h2>
<p>A protease does not simply recognise a sequence; it threads the substrate through its
active-site cleft, which requires the linker to be both exposed and extended. That is a
<b>necessary condition</b> for cleavage and it is independent of sequence — a linker packed
against the nanobody cannot be cut however good its motif. It is also computable from the folded
constructs, with no specificity data at all.</p></div>
<div class="scroll"><table class="wrap"><colgroup><col style="width:22%"><col style="width:22%">
<col style="width:19%"><col style="width:19%"><col style="width:18%"></colgroup>
<thead><tr><th>linker</th><th>design window (decades)</th><th>fraction exposed</th>
<th>extension ratio</th><th>accessibility</th></tr></thead><tbody>{_rows}</tbody></table></div>
<div class="note"><b>Every linker is exposed; not every linker is extended.</b> Solvent exposure
is near 1.0 throughout, so nothing is buried. What collapses is extension — the end-to-end
distance as a fraction of contour length falls from {_np.mean([x['extension_ratio'] for x in _by[_short]]):.2f}
at {_short} residues to {_np.mean([x['extension_ratio'] for x in _by[max(_by)]]):.2f} at
{max(_by)}. A long flexible linker balls up, and a balled-up linker is poorly presented to a
protease.</div>
{chart("trigger_tradeoff", "The two computable scores plotted against linker length.")}
<div class="note warn"><b>The two scores disagree, and that is the point.</b> The tether model
wants a <b>{_bestwin}-residue</b> linker, where the usable affinity window is widest at
{_win[_bestwin]:.2f} decades. Accessibility wants <b>{_short}</b>, where it scores
{_acc_short:.2f} against {_acc_long:.2f} at the longest. Neither is wrong. A single combined
number would hide the conflict and pick one arbitrarily; keeping them apart makes the trade-off
visible and leaves the choice to be made deliberately, or measured.</div>
<div class="note"><b>What is still missing.</b> This is the accessibility half of the trigger
score. The other half — does <em>this</em> protease cut <em>this</em> sequence, and leave other
sequences alone — needs substrate specificity data for the enzymes in question. Accessibility can
rule a linker out; it cannot rule one in.</div>
</section>"""
else:
    trigger_html = ""

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
    _m = {r["case"]: r for r in (complexes_msa["results"] if complexes_msa else [])
          if r["status"] == "ok"}
    _masked = [r for r in _m.values() if r.get("construct")]
    _msa_lo = min((r["fraction_of_contacts_on_paratope"] for r in _masked), default=0.0)
    _msa_hi = max((r["fraction_of_contacts_on_paratope"] for r in _masked), default=0.0)
    _mrows = "".join(
        f"<tr><td>{r['construct'] or 'free nanobody'}</td>"
        f"<td>{r['fraction_of_contacts_on_paratope']:.2f}</td>"
        f"<td>{r.get('mask_fraction_on_paratope', float('nan')):.2f}</td>"
        f"<td>{r['buried_surface_area_A2']:.0f}</td></tr>"
        if r.get("construct") else
        f"<tr class='good'><td>free nanobody (control)</td>"
        f"<td>{r['fraction_of_contacts_on_paratope']:.2f}</td><td>—</td>"
        f"<td>{r['buried_surface_area_A2']:.0f}</td></tr>"
        for r in _m.values())
    complex_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Step 04 · does it mask, second method</div>
<h2>Asking whether the antigen can still bind</h2>
<p>Occlusion measures whether the mask covers the site. It does not measure whether the antigen is
excluded, and those are different things. To separate them I predicted the whole assembly: the
nanobody, the tethered mask and the toxin together.</p>
<p>The first attempt gave nothing, because every chain was predicted without an alignment. The
control, which is the free nanobody with its toxin and whose answer the crystal structure already
gives, put only {ctrl['fraction_of_contacts_on_paratope']:.0%} of its contacts on the known
paratope. Giving the nanobody and the toxin alignments raised that to
{_m['control_free_vhh_msa']['fraction_of_contacts_on_paratope']:.0%}. The failure had been the
missing alignment rather than the model.</p>
<p>An alignment over the whole fusion would have sent the designed mask to a remote server, so I
split the construct into three chains joined by an explicit covalent bond. The nanobody and the
toxin are both already public in the PDB and take real alignments. The designed peptide takes
none and stays on this machine.</p></div>
<div class="scroll"><table class="wrap">
<colgroup><col style="width:34%"><col style="width:24%"><col style="width:22%">
<col style="width:20%"></colgroup>
<thead><tr><th>construct</th><th>toxin contacts on paratope</th>
<th>mask contacts on paratope</th><th>buried Å²</th></tr></thead>
<tbody>{_mrows}</tbody></table></div>
<div class="note warn"><b>None of these masks block the toxin.</b> In every masked construct the
toxin still reaches its own epitope, at {_msa_lo:.0%} to {_msa_hi:.0%} of contacts on the paratope
against {_m['control_free_vhh_msa']['fraction_of_contacts_on_paratope']:.0%} for the free
nanobody. The designed graft sits off the paratope altogether, and the inert poly-glycine control
has more contact with it than the graft does. The masks are tethered near the site but are not
occupying it.</div>
<div class="note"><b>This is what a working score is for.</b> Occlusion could not separate a
designed mask from a random one, so it could not have told me this. This measure reproduces a
known epitope when it is given one, and when it is given the current designs it says they do not
work. A negative answer from a calibrated measure is more useful than a positive one from an
uncalibrated measure.</div>
<div class="note"><b>Two limits worth stating.</b> Each case is a single prediction rather than a
sampled ensemble, and the model is given the mask and the toxin at the same time, so it describes
an assembly at equilibrium rather than a competition that unfolds over time. Both would need more
samples and a kinetic treatment before any design decision rests on them.</div>
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


_graft = _np.mean([r["occlusion"] * 100 for r in fold["results"]
                   if r["kind"] == "epitope graft"]) if fold else float("nan")
_inert = _np.mean([r["occlusion"] * 100 for r in fold["results"]
                   if r["kind"] == "inert control"]) if fold else float("nan")
_acc_rows = trigger["results"] if trigger else []
_acc_by = {}
for _r in _acc_rows:
    _acc_by.setdefault(_r["linker_residues"], []).append(_r["accessibility_score"])
_msa_ctrl = None
if complexes_msa:
    _msa_ctrl = next((r for r in complexes_msa["results"]
                      if r["case"].startswith("control") and r["status"] == "ok"), None)
_msa_masked = [r for r in (complexes_msa["results"] if complexes_msa else [])
               if r["status"] == "ok" and r.get("construct")]
_msa_lo = min((r["fraction_of_contacts_on_paratope"] for r in _msa_masked), default=0)
_msa_hi = max((r["fraction_of_contacts_on_paratope"] for r in _msa_masked), default=0)
_ctrl = None
if complexes:
    _ctrl = next((r for r in complexes["results"]
                  if r["case"] == "control_free_vhh" and r["status"] == "ok"), None)
ceff_mM = bN["c_eff_M"] * 1e3
from nbswitch.tether import fraction_toxin_bound as _ftb
pt_1nM = _ftb(1e-9) * 100

status_html = R.status_table([
    ("Measure the complex", "run",
     f"{meta['resolution_A']} Å structure, {iface['n_paratope']} paratope residues, "
     f"N-terminus {geom['n_term_to_nearest_paratope_A']} Å from the nearest of them."),
    ("Price the tether", "run",
     f"Gives a usable mask affinity between {bN['km_min_for_release_M']*1e9:.0f} nM and "
     f"{bN['km_max_for_masking_M']*1e6:.0f} µM, and says to attach at the N-terminus."),
    ("Build a construct library", "run",
     f"{lib['n_constructs']} constructs: grafted masks from the antigen, plus scrambled, "
     f"truncated and inert controls."),
    ("Fold them (ESMFold)", "run",
     f"{fold['n_folded'] if fold else 0} structures predicted on CPU."),
    ("Score masking by occlusion", "failed",
     f"Grafts occlude {_graft:+.1f}% of the paratope and the inert control {_inert:+.1f}%. "
     f"The control scores as well as the design, so the measure does not separate them."),
    ("Redesign masks (ProteinMPNN)", "failed",
     f"{redesign['n_unique'] if redesign else 0} masks, all scoring better than their parent by "
     f"ProteinMPNN's own measure, but that score does not predict occlusion (r = {corr:+.2f})."
     if redesign else "not run"),
    ("Score the trigger: can a protease reach the linker", "run",
     f"Works and separates the library. Accessibility falls from "
     f"{_np.mean(_acc_by[min(_acc_by)]):.2f} at {min(_acc_by)} linker residues to "
     f"{_np.mean(_acc_by[max(_acc_by)]):.2f} at {max(_acc_by)}."
     if _acc_by else "not run"),
    ("Score the trigger: will this protease cut this sequence", "proposed",
     "Needs substrate specificity data for the proteases in question. No such dataset exists "
     "for snake venom metalloproteinases."),
    ("Predict the masked complex (Boltz-2)", "run" if _msa_ctrl else "failed",
     (f"Without alignments the control put only "
      f"{_ctrl['fraction_of_contacts_on_paratope']:.0%} of its toxin contacts on the known "
      f"paratope, so the prediction was worthless. Giving the nanobody and antigen alignments "
      f"raises that to {_msa_ctrl['fraction_of_contacts_on_paratope']:.0%}."
      if _msa_ctrl and _ctrl else
      f"The control put only {_ctrl['fraction_of_contacts_on_paratope']:.0%} of its toxin "
      f"contacts on the known paratope, so it cannot be trusted." if _ctrl else "not run")),
    ("Ask whether the masks block the toxin", "run",
     f"They do not. With the nanobody and antigen aligned the toxin still reaches its own "
     f"epitope in every masked construct ({_msa_lo:.0%} to {_msa_hi:.0%} of contacts on the "
     f"paratope, against {_msa_ctrl['fraction_of_contacts_on_paratope']:.0%} for the free "
     f"nanobody), and the designed mask sits off the paratope entirely."
     if _msa_ctrl else "not run"),
    ("Rank designs across the scores", "proposed",
     "Two of the three scores work, so a Pareto front over all three cannot be built yet."),
    ("Make anything", "proposed",
     "Nothing has been synthesised. No measured affinity, expression or cleavage exists."),
])

HTML = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Designing a binder that is switched off until it arrives</title>
<meta name="description" content="A design method for conditionally active binders, and a worked
example on an experimental nanobody-toxin complex.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;650&display=swap" rel="stylesheet">
<style>{R.CSS}</style></head><body><div class="page">

<header>
<div class="kicker"><span>DTU Bioengineering</span><b>Interview case</b>
<span>Computational design</span></div>
<h1>A universal, protease-activated mask for nanobodies</h1>
<p class="lede">The case asked how to design a nanobody that stays inactive in circulation and
switches on at a snakebite site. I chose the protease as the trigger and made the regulatory part
reusable, so the same design works on other nanobodies. This page gives the design, then what
happened when I ran it.</p>
</header>

<section class="sec">
<div class="stack"><div class="eyebrow">Background</div>
<h2>Where I have met this problem before</h2>
<p>My PhD was on follistatin, a protein that carries its own tethered mask. The C-terminal tail of
FST315 folds back onto the domain that binds heparan sulphate and competes with it, so the molecule
has two states and the switch between them is set by how tightly the tail binds. I measured those
affinities and resolved the two states by SEC-SAXS. The design problem on this page is the same
problem, approached from the other direction.</p></div>
{background_html}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">The case</div>
<h2>The question, and what I chose to work on</h2>
<p>The brief described a viper bite site as having three local signals that are absent from
circulation: an acidic pH of around 6.0 to 6.5, raised protease activity from tissue damage and
from the venom itself, and local hypoxia. It asked for a computational pipeline that designs a
nanobody masked in circulation and activated at the bite site by one or more of these signals, and
asked whether a universal switch and standardised methods are possible.</p>
<p>I made two choices before designing anything.</p></div>

<div class="scroll"><table class="wrap">
<colgroup><col style="width:22%"><col style="width:46%"><col style="width:32%"></colgroup>
<thead><tr><th>signal</th><th>what the evidence shows</th><th>decision</th></tr></thead><tbody>
<tr><td><b>Protease</b></td><td>Gelatinolytic activity in exudate next to envenomed muscle peaks
about an hour after the bite and comes from the venom's own metalloproteinases. Host enzymes
arrive later. The signal is therefore early, abundant and proportional to the thing being
neutralised.</td><td><b>Chosen as the trigger.</b></td></tr>
<tr><td>pH</td><td>No study has measured pH at a bite site. The nearest proxies, inflamed and
infected tissue, sit near pH 6.8 to 7.0 rather than 6.0. A mask switching at 6.0 would rarely
open.</td><td>Not used. I would measure it before designing against it.</td></tr>
<tr><td>Hypoxia</td><td>Argued from microvascular damage, never measured in envenomed tissue. No
validated protein switch responds directly to oxygen.</td><td>Not used.</td></tr>
</tbody></table></div>

<div class="note"><b>The second choice was to make the regulatory part reusable.</b> A single mask
sequence that works on every nanobody is unlikely, because paratopes differ. What can transfer is
the architecture: a fixed attachment point, a mask family, an interchangeable cleavable linker and
one scoring protocol. The nanobody stays unchanged and only the mask is redesigned, so the same
pipeline can be pointed at a different toxin without rebuilding it.</p></div>

<div class="chart wide">{R.states_schematic()}</div>
<div class="note"><b>Why the mask has to be weak.</b> The tether holds the mask near the site at
an effective concentration of {ceff_mM:.0f} mM, so a mask with a dissociation constant of
{bN['km_max_for_masking_M']*1e6:.0f} µM still blocks almost all of the time. After cleavage that
advantage is gone and the mask has to leave. A mask at 1 nM would block beautifully and then stay
put, so only {pt_1nM:.0f}% of the nanobody would ever reach toxin. Every protein design tool
maximises affinity, which means the obvious approach produces masks that score as successes and do
not work.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">The design</div>
<h2>Five modules</h2>
<p>The pipeline generates constructs, scores them three separate ways, and combines the three
scores into one ranked list. Each module is a slot. The methods listed under it are
interchangeable, so a better generator or a better cleavage model can be dropped in without
changing anything else.</p></div>
<div class="chart wide">{R.architecture_schematic()}</div>
{R.module_tables()}
<div class="note"><b>The three scores are kept apart on purpose.</b> A construct that cannot be
expressed is of no use however well it masks, and one that cannot be cut is of no use however well
it expresses. Averaging them would let a high score on one hide a failure on another, so the
combining stage takes the Pareto front instead.</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">What I ran</div>
<h2>The example, step by step</h2>
<p>I ran the pipeline on an experimental nanobody-toxin complex from the same laboratory.
Everything below is computational and nothing has been made, so no affinity, expression level or
cleavage rate here has been measured.</p></div>
{status_html}
</section>

{cautionary_html}

<section class="sec">
<div class="stack"><div class="eyebrow">Step 01 · the reference structure</div>
<h3>Measuring the complex</h3>
<p>The paratope is defined as the residues that lose solvent-accessible surface when the complex
forms, and the geometry that matters for a tether is how far each terminus sits from that
surface.</p></div>
<div class="cards">
{card(f"{iface['interface_area_A2']:.0f} Å²", "interface area", f"{repro['bsa_mean_A2']} ± {repro['bsa_sd_A2']} Å² buried across {repro['n_copies']} copies")}
{card(iface['n_paratope'], "paratope residues", "burying &gt;1 Å² on binding")}
{card(f"{region['CDR3']}/{iface['n_paratope']}", "of it in CDR3", f"CDR1 {region['CDR1']} · CDR2 {region['CDR2']} · framework {region['framework']}", hl=True)}
{card(ref['cdrs']['CDR3']['length'], "CDR3 length", "residues")}
</div>
{chart_wide("paratope_footprint", "Surface buried per nanobody residue on complex formation.")}
<div class="banner"><b>What this fixes before anything is designed.</b> The binding site spans
{region_summary}. The N-terminus sits <b>{geom['n_term_to_nearest_paratope_A']} Å</b> from the
nearest paratope residue and the C-terminus <b>{geom['c_term_to_nearest_paratope_A']} Å</b> — so a
mask tethered at the N-terminus has a real gap to span, and the two termini are not
interchangeable. {n_term_verdict}</div>
<div class="scroll"><table class="wrap"><thead><tr><th>geometry</th><th>Å</th><th>consequence</th></tr></thead>
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
<div class="stack"><div class="eyebrow">Step 01 · the design window</div>
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
{card(f"{bN['window_decades']:.1f}", "decades of usable K<sub>D</sub>", f"N-terminal, {bN['linker_residues']}-residue linker", hl=True)}
{card(f"{bC['window_decades']:.1f}", "decades, C-terminal", f"best at {bC['linker_residues']} residues")}
{card(f"≤{closed_C}", "residues where C-terminal fails", "no window at any affinity")}
{card(f"{bN['km_max_for_masking_M']*1e6:.0f} µM", "loosest mask that still blocks", f"N-terminal, {bN['linker_residues']} residues")}
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
<div class="stack"><div class="eyebrow">Step 02 · generate constructs</div>
<h3>Candidates, and the controls that make them legible</h3>
<p>One peptide is already known to engage this paratope: the antigen's own contact segment. Grafting
it needs no generative model and brings its own negative controls.</p></div>
<div class="charts2">
<div class="scroll"><table><thead><tr><th>antigen segment</th><th>sequence</th><th>length</th>
<th>contacts</th></tr></thead><tbody>{seg_rows}</tbody></table></div>
<div class="scroll"><table class="wrap"><colgroup><col style="width:15%"><col style="width:30%">
<col style="width:8%"><col style="width:20%"><col style="width:27%"></colgroup>
<thead><tr><th>mask</th><th>sequence</th><th>aa</th><th>kind</th>
<th>origin</th></tr></thead><tbody>{masks_rows}</tbody></table></div>
</div>
<div class="banner"><b>On the candidates.</b> <span class="mono">epi1</span> is the tip of
α-cobratoxin loop II, the region carrying the residues that drive receptor binding — this nanobody
neutralizes by covering the toxin's functional site, so that same peptide is the obvious mask.
<span class="mono">epi2</span> carries a free cysteine and is kept as a deliberate counter-example.
<span class="mono">polyG</span> has no designed interaction and should fail.</div>
</section>

{redesign_html}

{trigger_html}

<section class="sec">
<div class="stack"><div class="eyebrow">Step 04 · does it mask</div>
<h3>Measuring whether the site is covered</h3></div>
<div class="chart wide">{R.occlusion_schematic()}</div>
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

{complex_html}

{ranking_html}

<section class="sec">
<div class="stack"><div class="eyebrow">What the example does not settle</div>
<h2>The gaps, named</h2></div>
<div class="scroll"><table class="wrap"><thead><tr><th>module</th><th>what is missing</th>
<th>what it would take</th></tr></thead><tbody>
<tr><td><b>02a Trigger</b></td><td>no selectivity score at all</td><td><span class="tag gap">data</span>
cleavage measurements for the relevant proteases against a host negative panel. Published substrate
profiling puts a P1&prime; leucine preference on both the venom metalloproteinases and the host MMPs,
so a linker designed to be cut by one would likely be cut by the other. Selectivity, not
cleavability, is what has to be measured and what no current model predicts</td></tr>
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
<div class="scroll"><table class="wrap"><thead><tr><th>stage</th><th>method</th><th>produces</th></tr></thead>
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

{R.scale_section()}

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


import re

# ---- slide framing -------------------------------------------------------
_titles = re.findall(r'<div class="eyebrow">(.*?)</div>', HTML)
_sections = HTML.count('<section class="sec"')
_n = 0
def _stamp(m):
    global _n
    _n += 1
    return f'<section class="sec" id="s{_n}" data-slide="{_n:02d} / {_sections:02d}"'
HTML = re.sub(r'<section class="sec"', _stamp, HTML)

_dots = "".join(f'<a href="#s{i}" data-s="{i}">{i:02d}</a>' for i in range(1, _sections + 1))
_nav = ('<nav class="navbar"><div class="inner"><div class="dots">' + _dots +
        '</div><span class="hint">\u2190 \u2192 to move</span></div></nav>'
        '<script>(function(){'
        'var secs=[].slice.call(document.querySelectorAll("section.sec"));'
        'var dots=[].slice.call(document.querySelectorAll(".navbar .dots a"));'
        'function mark(){var y=window.scrollY+120,k=0;'
        'secs.forEach(function(s,i){if(s.offsetTop<=y)k=i;});'
        'dots.forEach(function(d,i){d.setAttribute("aria-current",i===k?"true":"false");});}'
        'function go(d){var y=window.scrollY+120,k=0;'
        'secs.forEach(function(s,i){if(s.offsetTop<=y)k=i;});'
        'var t=Math.max(0,Math.min(secs.length-1,k+d));'
        'window.scrollTo({top:secs[t].offsetTop-40,behavior:"smooth"});}'
        'document.addEventListener("keydown",function(e){'
        'if(e.target.tagName==="INPUT"||e.target.tagName==="SELECT")return;'
        'if(e.key==="ArrowRight"||e.key==="PageDown"){e.preventDefault();go(1);}'
        'if(e.key==="ArrowLeft"||e.key==="PageUp"){e.preventDefault();go(-1);}});'
        'window.addEventListener("scroll",mark,{passive:true});mark();})();</script>')
HTML = HTML.replace("</div>\n\n<script src=", _nav + "</div>\n\n<script src=", 1)

import re as _re

# Build-time guard. Scripts and styles legitimately contain braces, so they are
# removed before the scan; what is left is markup, where a surviving "{name[" or
# "{R.fn()" means an f-string was escaped when it should have interpolated.
_markup = _re.sub(r"<script\b.*?</script>", "", HTML, flags=_re.S | _re.I)
_markup = _re.sub(r"<style\b.*?</style>", "", _markup, flags=_re.S | _re.I)
_markup = _re.sub(r"<svg\b.*?</svg>", "", _markup, flags=_re.S | _re.I)
_LEAK = _re.compile(r"\{\s*(?:R\.[a-z_]+\(|[A-Za-z_][A-Za-z0-9_]*\s*\[['\"])")
_hits = list(_LEAK.finditer(_markup))
if _hits:
    _spots = [_markup[m.start():m.start() + 70].replace("\n", " ") for m in _hits[:6]]
    raise SystemExit("unresolved template expressions in the output:\n  "
                     + "\n  ".join(_spots))

(DOCS / "index.html").write_text(HTML)
print(f"docs/index.html  {len(HTML)//1024} KB")
print(f"  {len(list((DOCS/'figures').glob('*.svg')))} figures, "
      f"{len(list((DOCS/'structures').glob('*.pdb')))} structures, "
      f"{len(models)} in the viewer")
