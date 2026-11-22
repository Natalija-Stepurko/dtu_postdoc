#!/usr/bin/env python
"""Stage 6. Render docs/index.html from the stage outputs."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import numpy as _np
from nbswitch import report as R
from nbswitch.report import mri
from nbswitch.config import (DOCS, PDB_ID, REPO_URL, NOTEBOOK_URL, load_json,
                             REFERENCE_VHH_CHAIN, TARGET_ID, TARGET, TARGETS,
                             RESULTS, results_dir)

ref = load_json("01_reference.json")

def _check_library_matches_structures():
    """Refuse to build a page whose numbers and structures disagree.

    The construct library and the folded coordinates are separate artefacts.
    Changing the library without refolding leaves residue ranges pointing at
    the wrong atoms, which renders as a page that looks fine and is wrong.
    """
    from nbswitch.config import STRUCTURES
    lib = load_json("03_constructs.json")
    for item in lib["constructs"]:
        path = STRUCTURES / f"{item['id']}.pdb"
        if not path.exists():
            continue
        last = 0
        for line in path.read_text().splitlines():
            if line.startswith("ATOM"):
                last = max(last, int(line[22:26]))
        if last and last != item["vhh_range"][1]:
            raise SystemExit(
                f"library and structures disagree: {item['id']} ends at residue "
                f"{item['vhh_range'][1]} in 03_constructs.json but the folded "
                f"model has {last}. Refold (stages 04 and 08) before building "
                f"the page.")


_check_library_matches_structures()
tether = load_json("02_tether.json")
lib = load_json("03_constructs.json")
try:
    fold = load_json("04_folding.json")
except FileNotFoundError:
    fold = None
for _name, _key in (("02b_linkers.json", "linkers"),
                    ("13_ranking.json", "ranking_data"),
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
    n_masks=len(lib["masks"]), n_linkers=lib["n_linkers"],
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

background_html = R.background_section()

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

# the picker offers what the library actually contains. Building it from the
# directory listing offered leftovers from earlier runs, whose ids are absent
# from MODELS, so selecting one threw and the viewer stayed empty.
viewer_ids = sorted(models)
_masked = [i for i in viewer_ids if models[i].get("mask")]
default_id = (_masked[0] if _masked else (viewer_ids[0] if viewer_ids else ""))

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
    _rows_all = trigger["results"]
    _gen = [r for r in _rows_all if r.get("linker_kind") == "generated"]
    _ctrl = [r for r in _rows_all if "control" in str(r.get("linker_kind", ""))]
    _win = {r["linker_residues"]: r["window_decades"]
            for r in tether["scan"] if r["terminus"] == "N"}
    _bestwin = max(_win, key=_win.get) if _win else 0
    _len = _gen[0]["linker_residues"] if _gen else 0

    def _mean(rows, key):
        vals = [r[key] for r in rows if r.get(key) is not None]
        return float(_np.mean(vals)) if vals else float("nan")

    # if reach tracked recognition, one of the two measures would be redundant
    _wl = [r for r in _rows_all if r.get("linker_id")]
    _r_sel_ctx = (float(_np.corrcoef([r["selectivity"] for r in _wl],
                                     [r["context_score"] for r in _wl])[0, 1])
                  if len(_wl) > 2 else float("nan"))
    _worst = min(_wl, key=lambda r: r["cleavability_score"]) if _wl else None

    _parts = "".join(
        f"<tr><td>{label}</td><td>{len(rows)}</td>"
        f"<td>{_mean(rows, 'extension'):.2f}</td>"
        f"<td>{_mean(rows, 'exposure'):.2f}</td>"
        f"<td>{_mean(rows, 'clearance_open_fraction'):.2f}</td>"
        f"<td><b>{_mean(rows, 'context_score'):.3f}</b></td></tr>"
        for label, rows in (("generated substrates", _gen),
                            ("controls", _ctrl)) if rows)

    trigger_html = f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Step 03 · can it be cut</div>
<h2>Can the protease reach the linker at all?</h2>
<p>Recognition and reach are different questions. A protease threads its substrate through an
active-site cleft, so the linker has to lie extended, stay exposed, and leave room beside the
scissile bond for a catalytic domain roughly forty &#8491; across. Each is a
<b>necessary condition</b>: a linker packed against the nanobody cannot be cut however good its
motif. All three are computable from the folded constructs alone.</p>
<p>The two axes turn out to be close to independent across the library
(<b>r&nbsp;=&nbsp;{_r_sel_ctx:+.2f}</b> between sequence selectivity and structural reach,
n&nbsp;=&nbsp;{len(_wl)}). Neither one substitutes for the other, which is the reason they are
multiplied rather than averaged into a single number.</p></div>
<div class="scroll"><table class="wrap">
<thead><tr><th>linkers</th><th>n</th><th>extension</th><th>exposure</th><th>clearance</th>
<th>product</th></tr></thead><tbody>{_parts}</tbody></table></div>
<div class="note"><b>Exposure is not the binding constraint.</b> Solvent exposure sits near
{_mean(_rows_all, 'exposure'):.2f} across the library, so nothing is buried. What separates
constructs is extension and clearance &mdash; whether the segment lies in a groove rather than
balled up, and whether an enzyme-sized object can get to it.</div>
<div class="note"><b>Both axes have to hold at once.</b> The worst construct in the library,
{_worst['id']}, is physically reachable &mdash; its context score of {_worst['context_score']:.2f}
is mid-pack &mdash; but its linker is a scramble of a good substrate, so recognition collapses to
{_worst['selectivity']:+.2f} and the product with it, to {_worst['cleavability_score']:.3f}. Same
composition as something cleavable, same accessibility, no site.</div>
<div class="note"><b>The generated linkers are one length, and it is close to the right one.</b>
Every substrate is ten residues on three-residue flanks, so {_len} in all. The tether model wants
<b>{_bestwin}</b> residues, where the usable affinity window is widest at {_win.get(_bestwin, float('nan')):.2f}
decades. That agreement was not designed in: the substrate length is CleaveNet's, and the window
comes from the worm-like chain model, which knows nothing about proteases.</div>
<div class="note warn"><b>What this still cannot tell you.</b> Generation and verification share a
model family, so a high selectivity score is weaker evidence than two independent models would
give. More importantly, CleaveNet is trained on human matrix metalloproteinases, and the enzyme
that arrives first at a bite site is a snake venom metalloproteinase, which is not in its training
set. The gelatinases are a documented stand-in, not the enzyme itself.</div>
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
                  for i in viewer_ids)


_graft = _np.mean([r["occlusion"] * 100 for r in fold["results"]
                   if r["kind"] == "epitope graft"]) if fold else float("nan")
_inert = _np.mean([r["occlusion"] * 100 for r in fold["results"]
                   if r["kind"] == "inert control"]) if fold else float("nan")
_acc_rows = trigger["results"] if trigger else []
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



# ---------------------------------------------------------------- figures
N = {}
N["resolution"] = meta["resolution_A"]
N["paratope"] = iface["n_paratope"]
N["iface_area"] = iface["interface_area_A2"]
N["copies"] = repro["n_copies"]
N["bsa_mean"], N["bsa_sd"] = repro["bsa_mean_A2"], repro["bsa_sd_A2"]
N["cdr3"] = region["CDR3"]
N["n_term"] = geom["n_term_to_nearest_paratope_A"]
N["c_term"] = geom["c_term_to_nearest_paratope_A"]
N["window"] = bN["window_decades"]
N["linker_best"] = bN["linker_residues"]
N["km_lo"] = bN["km_min_for_release_M"] * 1e9
N["km_hi"] = bN["km_max_for_masking_M"] * 1e6
N["ceff"] = bN["c_eff_M"] * 1e3
N["n_constructs"] = lib["n_constructs"]
N["n_folded"] = fold["n_folded"] if fold else 0

_k = {}
for _r in (fold["results"] if fold else []):
    _k.setdefault(_r["kind"], []).append(_r["occlusion"] * 100)
N["occ_graft"] = _np.mean(_k.get("epitope graft", [0]))
N["occ_scram"] = _np.mean(_k.get("scrambled control", [0]))
N["occ_inert"] = _np.mean(_k.get("inert control", [0]))

N["n_redesign"] = redesign["n_unique"] if redesign else 0
N["mpnn_corr"] = corr if redesign else float("nan")

_ctx = [_r["context_score"] for _r in (trigger["results"] if trigger else [])]
N["ctx_mean"] = float(_np.mean(_ctx)) if _ctx else 0.0
N["ctx_max"] = float(max(_ctx)) if _ctx else 0.0

# recognition against reach. If these tracked each other one measure would be
# redundant; the product would then be doing no more work than either alone.
_lin_rows = [_r for _r in (trigger["results"] if trigger else []) if _r.get("linker_id")]
N["n_cleav"] = len(_lin_rows)
N["sel_ctx_r"] = (float(_np.corrcoef([_r["selectivity"] for _r in _lin_rows],
                                     [_r["context_score"] for _r in _lin_rows])[0, 1])
                  if len(_lin_rows) > 2 else float("nan"))
_least = min(_lin_rows, key=lambda r: r["cleavability_score"]) if _lin_rows else None
N["least_id"] = _least["id"] if _least else ""
N["least_ctx"] = _least["context_score"] if _least else 0.0
N["least_sel"] = _least["selectivity"] if _least else 0.0
N["least_cleav"] = _least["cleavability_score"] if _least else 0.0

_mm = {r["case"]: r for r in (complexes_msa["results"] if complexes_msa else [])
       if r["status"] == "ok"}
_mk = [r for r in _mm.values() if r.get("construct")]
N["boltz_ctrl_noalign"] = (_ctrl["fraction_of_contacts_on_paratope"] if _ctrl else 0) * 100
N["boltz_ctrl"] = (_mm["control_free_vhh_msa"]["fraction_of_contacts_on_paratope"] * 100
                   if "control_free_vhh_msa" in _mm else 0)
N["boltz_lo"] = min((r["fraction_of_contacts_on_paratope"] for r in _mk), default=0) * 100
N["boltz_hi"] = max((r["fraction_of_contacts_on_paratope"] for r in _mk), default=0) * 100
N["boltz_rows"] = "".join(
    (f"<tr class='good'><td>free nanobody, control</td>"
     f"<td>{r['fraction_of_contacts_on_paratope']:.2f}</td><td>—</td></tr>"
     if not r.get("construct") else
     f"<tr><td>{r['construct']}</td>"
     f"<td>{r['fraction_of_contacts_on_paratope']:.2f}</td>"
     f"<td>{r.get('mask_fraction_on_paratope', 0):.2f}</td></tr>")
    for r in _mm.values())

try:
    N["n_ensemble"] = load_json("14_ensemble.json")["n_models"] - 1
except FileNotFoundError:
    N["n_ensemble"] = 0
try:
    _ens = load_json("14_ensemble.json")
    N["n_ensemble"] = _ens["n_models"] - 1
    ENSEMBLE_IDS = _ens["constructs"]
except FileNotFoundError:
    N["n_ensemble"] = 0
    ENSEMBLE_IDS = []

def _norm(values):
    vals = [v for v in values.values() if v is not None]
    lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
    span = (hi - lo) or 1.0
    return {k: (None if v is None else (v - lo) / span) for k, v in values.items()}

_cdr_member = {r: k for k, v in ref["cdrs"].items() for r in v["residues"]}
_para_labels = sorted(iface["paratope_labels"],
                      key=lambda l: int("".join(c for c in l[1:] if c.isdigit())))
N["n_paratope"] = len(_para_labels)
N["n_cdr_para"] = sum(1 for l in _para_labels if l in _cdr_member)
N["n_fw_para"] = N["n_paratope"] - N["n_cdr_para"]
N["fw2_para"] = ", ".join(l for l in _para_labels
                          if l not in _cdr_member
                          and 39 <= int("".join(c for c in l[1:] if c.isdigit())) <= 52)

cdr_resi = {name: [int("".join(c for c in l[1:] if c.isdigit())) for l in v["residues"]]
            for name, v in ref["cdrs"].items()}

top_paratope = [
    {"resi": int("".join(c for c in lbl[1:] if c.isdigit())), "label": lbl}
    for lbl in iface["paratope_contact_labels"]]
N["n_contacts"] = iface["n_paratope_contacts"]
N["contact_cutoff"] = iface["contact_cutoff_A"]
N["contact_labels"] = ", ".join(iface["paratope_contact_labels"])

SCORES = {}
_feas_by = {r["id"]: r["feasibility_score"] for r in (feas_data["results"] if feas_data else [])}
_occ_by = {}
for _src in (fold, redesign_fold):
    for _r in (_src["results"] if _src else []):
        _occ_by[_r["id"]] = _r["occlusion"]
_acc_by_id = {r["id"]: r["context_score"] for r in (trigger["results"] if trigger else [])}
for _key, _d in (("feasibility", _feas_by), ("masking", _occ_by), ("cleavability", _acc_by_id)):
    _n = _norm({it["id"]: _d.get(it["id"]) for it in ENSEMBLE_IDS})
    for _i, _v in _n.items():
        SCORES.setdefault(_i, {})[_key] = _v

controls = load_json("15_interface_controls.json")
_by = {r["pdb_id"]: r for r in controls["comparators"]}
N["cdr3_share"] = _by.get(PDB_ID, {}).get("cdr3_percent_bsa", 0)
N["iface_area"] = iface["interface_area_A2"]
_lat = _by.get(PDB_ID, {}).get("lattice_contact", {})
N["lattice_area"] = _lat.get("interface_area_A2", 0)
N["lattice_n"] = _lat.get("n_paratope", 0)
N["cdr3_other"] = _by.get("9GCN", {}).get("cdr3_percent_bsa", 0)
def _count_rows():
    """Read the per-stage counts back out of the stage outputs."""
    rows = []
    try:
        lk = load_json("02b_linkers.json")
        rows.append({"stage": "02", "what": "protease substrates generated, "
                     "conditioned on the target profile", "made": lk["n_generated"],
                     "kept": lk["n_kept"]})
        rows.append({"stage": "02", "what": "controls: a protease-resistant "
                     "Gly/Ser repeat and a scramble of the best substrate",
                     "made": lk.get("n_controls", 0),
                     "kept": lk.get("n_controls", 0)})
        rows.append({"stage": "02", "what": "linkers entering the library",
                     "made": len(lk["linkers"]), "kept": len(lk["linkers"])})
    except FileNotFoundError:
        pass
    try:
        cl = load_json("03_constructs.json")
        rows.append({"stage": "02", "what": "masks: epitope grafts, a truncation, "
                     "a scramble and an inert poly-glycine",
                     "made": len(cl["masks"]), "kept": len(cl["masks"])})
        rows.append({"stage": "02", "what": "constructs: every mask on every "
                     "linker, plus the free nanobody",
                     "made": cl["n_constructs"], "kept": cl["n_constructs"]})
    except FileNotFoundError:
        pass
    for name, stage, what in (
            ("04_folding.json", "03", "constructs folded with ESMFold"),
            ("07_redesign.json", "02", "mask sequences redesigned with ProteinMPNN"),
            ("08_redesign_folding.json", "03", "redesigns folded with ESMFold"),
            ("11_complexes_msa.json", "04", "complexes predicted with alignments"),
    ):
        try:
            d = load_json(name)
        except FileNotFoundError:
            # a stage that has not run has no output at all; a stage that has
            # run but is named wrong here looks identical, and silently drops
            # a row from the table. Tell them apart.
            if list(RESULTS.glob(f"{name.split('_')[0]}_*.json")):
                raise SystemExit(
                    f"{name} is missing but stage {name.split('_')[0]} has output "
                    f"({[f.name for f in RESULTS.glob(name.split('_')[0] + '_*.json')]}). "
                    f"Fix the filename in _count_rows.")
            continue
        # n_folded counts the free-nanobody reference too; results does not
        n = d.get("n_folded") or len(d.get("results", d.get("constructs", [])))
        rows.append({"stage": stage, "what": what, "made": n, "kept": None})
    try:
        rk = load_json("13_ranking.json")
        rows.append({"stage": "06", "what": "carried all three scores and were "
                     "ranked", "made": rk["n_total"],
                     "kept": rk["n_front"]})
    except (FileNotFoundError, KeyError):
        pass
    rows.sort(key=lambda r: r["stage"])
    return rows


try:
    _lk = load_json("02b_linkers.json")
    _gen = [x for x in _lk["linkers"] if x.get("kind") == "generated"]
    _sel = [x["selectivity"] for x in _gen if x.get("selectivity") is not None]
    _ctrl = next((x["selectivity"] for x in _lk["linkers"]
                  if x["id"] == "ctrl_gs" and x.get("selectivity") is not None), 0.0)
    N["n_linkers"] = len(_lk["linkers"])
    N["n_generated"] = _lk["n_generated"]
    N["sel_min"] = min(_sel) if _sel else 0.0
    N["sel_max"] = max(_sel) if _sel else 0.0
    N["sel_ctrl"] = _ctrl
    N["target_proteases"] = ", ".join(_lk["target_proteases"])
    N["off_target_proteases"] = ", ".join(_lk["off_target_proteases"])
    LINKERS = _lk["linkers"]
except FileNotFoundError:
    N.update({"n_linkers": 0, "n_generated": 0, "sel_min": 0.0, "sel_max": 0.0,
              "sel_ctrl": 0.0, "target_proteases": "", "off_target_proteases": ""})
    LINKERS = []

COUNTS = _count_rows()
N["rank_total"] = ranking_data["n_total"] if ranking_data else 0
N["rank_front"] = ranking_data["n_front"] if ranking_data else 0

_feas = {r["id"]: r for r in (feas_data["results"] if feas_data else [])}
_graft_feas = [r for r in _feas.values() if r["mask"] == "epi1"]
N["feas_graft"] = _np.mean([r["feasibility_score"] for r in _graft_feas]) if _graft_feas else 0
N["feas_cys"] = _graft_feas[0]["unpaired_cysteines"] if _graft_feas else 0
N["rmsd_mean"] = _np.mean([r["vhh_core_rmsd_A"] for r in _feas.values()
                           if r.get("vhh_core_rmsd_A")]) if _feas else 0


masks_rows_compact = "".join(
    f"<tr><td>{k}</td><td class='mono'>{m['seq']}</td><td>{m['kind']}</td></tr>"
    for k, m in lib["masks"].items())

HTML = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Computational design of protease-responsive nanobodies</title>
<meta name="description" content="A computational approach to designing nanobodies whose antigen
binding is restored by proteolysis, with a partial implementation on an experimental complex.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;650&display=swap" rel="stylesheet">
<style>{R.CSS}</style></head><body><div class="page">

<header>
<div class="kicker"><span>Natalija Stepurko</span><b>DTU · postdoctoral interview</b></div>
<h1>Computational design of protease-responsive nanobodies</h1>
<p class="lede">A proposed pipeline for designing nanobodies that are masked in circulation and
restored to antigen binding by proteolysis, together with a partial implementation on an
experimental nanobody–antigen complex.</p>
</header>

<section class="sec">
<div class="stack"><div class="eyebrow">01</div><h2>Previous research</h2>
<p>During my PhD at the University of Cambridge I investigated how structural components of
TGF-β superfamily receptors and their extracellular inhibitors contribute to biological function.
I used SAXS-based structural modelling and biophysical assays to characterise the C-terminal tail
of follistatin, showing that it interacts with the FS1 domain and alters heparan sulphate binding,
and that SEC-SAXS data are best described by two distinct solution conformations. I also
investigated the contribution of receptor transmembrane helices to assembly and signalling, for
which I developed a bacterial signalling assay for transmembrane helix heterodimerisation.</p>
<p>The follistatin work is directly relevant here, since it concerns an intramolecular interaction
that regulates access to a protein-binding surface, which is the arrangement the present design
introduces deliberately.</p></div>
{background_html}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">02</div><h2>Research question and molecular architecture</h2>
<p>The question concerned nanobodies against snake venom toxins that would remain inactive in
circulation and become active at the site of a bite, where the local environment differs from
plasma. Three environmental cues were proposed: acidic pH, raised protease activity and local
hypoxia.</p></div>

<div class="scroll"><table class="wrap">
<colgroup><col style="width:22%"><col style="width:52%"><col style="width:26%"></colgroup>
<thead><tr><th>environmental cue</th><th>consideration</th><th>decision</th></tr></thead><tbody>
<tr><td>Acidic pH</td><td>The physiological range at a bite site would require further
characterisation before a switching threshold could be chosen.</td><td>Not selected</td></tr>
<tr class="good"><td>Protease activity</td><td>Proteolysis is irreversible, which makes activation a
single committed event rather than an equilibrium, and venom metalloproteinases are themselves part
of the pathology being treated.</td><td>Selected</td></tr>
<tr><td>Hypoxia</td><td>Would require an oxygen-sensitive molecular mechanism suitable for a small
secreted protein.</td><td>Not selected</td></tr>
</tbody></table></div>

<p>I propose to develop a computational approach for designing protease-responsive nanobodies
while preserving their original antigen-binding properties. Responsiveness is introduced by
attaching a masking domain through a protease-cleavable linker. The mask should prevent antigen
binding while attached and dissociate once the linker is cleaved. To make the approach applicable
to more than one nanobody, the antigen-binding domain is left unchanged and only the regulatory
components are designed, so that the architecture and the computational methods transfer even
though individual mask sequences will not.</p>
<div class="chart wide">{R.states_schematic()}</div>
<p>The two states impose opposing requirements on the mask. Because the mask is tethered it is held
near the paratope at a high local concentration, estimated at {N['ceff']:.0f} mM for this
construct, so a weak interaction is sufficient to block antigen binding. After cleavage that
local concentration is lost and the mask must dissociate, which an interaction of high affinity
would not do. The design target is therefore a bounded range of mask affinity rather than a
maximum, which is the opposite of the objective most protein design methods optimise.</p>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">03</div><h2>Proposed computational pipeline</h2>
<p>The pipeline consists of five modules. Candidate constructs are generated, three independent
scoring functions assess structural feasibility, the interaction between mask and nanobody, and
the susceptibility of the linker to proteolysis, and the three scores are combined by multiobjective
ranking. The output is a list of candidate sequences with their associated scores. Each module can
use several computational methods, so individual methods can be compared or replaced without
altering the pipeline.</p></div>
<div class="chart wide">{R.architecture_schematic()}</div>
{R.module_summary_table()}
<div class="eq">S(x) = [ S<sub>structure</sub>(x), S<sub>binding</sub>(x), S<sub>cleavage</sub>(x) ]</div>
<p>The three scores are not combined into a single value, since a high score on one objective
would otherwise compensate for failure on another. Pareto optimisation identifies the nondominated
candidates, which is a set rather than an ordered list, so a further criterion would be required to
order candidates within a front.</p>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">04</div><h2>Worked computational example</h2>
<p>The pipeline was implemented in part and applied to {meta['title'][:66]}, a {N['resolution']} Å
X-ray structure deposited in {meta['deposition_date'][:4]}. All results below are computational.
No construct has been expressed and no affinity, expression level or cleavage rate has been
measured.</p>
<p>What was produced at each stage, and how much of it survived:</p>
{R.generation_table(COUNTS)}
<p class="dim" style="font-size:.9rem">Every figure above is read from the output of the stage
that produced it.</p></div>

<h3 class="mod"><span>Reference</span>Structural characterisation of the complex</h3>
{mri("Solvent-accessible surface was calculated for the nanobody alone and in complex, and "
     "residues losing more than 1 Å<sup>2</sup> on binding were defined as the paratope. "
     "Distances from each terminus to that surface were measured.",
     f"The paratope comprises {N['paratope']} residues burying {N['iface_area']:.0f} Å<sup>2</sup>. "
     f"The N-terminus lies {N['n_term']} Å from the nearest paratope residue and the C-terminus "
     f"{N['c_term']} Å. Across the {N['copies']} VHH–antigen pairs in the entry the buried area is "
     f"{N['bsa_mean']:.0f} ± {N['bsa_sd']:.0f} Å<sup>2</sup>.",
     "The N-terminus is closer to the antigen-binding surface and is therefore the attachment "
     "point used throughout. The agreement between the two pairs indicates the interface is "
     "consistently resolved within this crystal.")}
<div class="viewer">
<div class="vbar"><span class="eyebrow">Deposited structure · {PDB_ID}</span>
<button data-spin="solved">spin</button></div>
<div class="vstage" id="stage-solved"></div>
<div class="vkey"><span><i style="background:#337D79"></i>framework</span>{R.CDR_KEY}
<span><i style="background:#C99555"></i>short neurotoxin 1</span>
<span><i style="background:#C06014"></i>buried on binding, sticks</span>
<span><i style="background:#C06014"></i>within {N['contact_cutoff']:.0f} &#8491;, labelled</span>
<span><b style="color:#12191F">N</b> / <b style="color:#6E7880">C</b> termini</span></div>
</div>
{chart_wide("paratope_footprint",
  "Surface buried per nanobody residue on complex formation. The six residues burying the most "
  "surface are labelled in the viewer above.")}
<p>Two cutoffs describe this interface and they answer different questions. Taking a heavy-atom
distance of {N['contact_cutoff']:.0f} &#8491;, the range of hydrogen bonds and salt bridges, gives
{N['n_contacts']} residues: {N['contact_labels']}. These are labelled in the viewer and marked on the
plot. Taking every residue that loses solvent-accessible surface gives {N['n_paratope']}, because
van der Waals packing reaches further than a hydrogen bond. The smaller set is what is interacting
chemically; the larger one is the surface a mask has to cover, so the pipeline targets the larger
one. On one of the other complexes examined the {N['contact_cutoff']:.0f} &#8491; criterion returns a
single residue, which is too little to define a footprint from.</p>
<p>One interface measured on its own says little, so the same quantities were measured on
independent VHH complexes: two further nanobody&ndash;three-finger-toxin structures from other
groups, and a nanobody bound to a much larger toxin.</p>
{R.interface_controls_table(controls, PDB_ID)}
<p>One of the comparators is biparatopic: the same VHH binds through two surfaces that share a
single residue, so its row describes one of two sites. That architecture is a poor candidate for
this pipeline, since a mask would have to cover both surfaces at once, and it is listed here as a
comparator rather than a target for that reason. The worked example is not biparatopic: its second
toxin copy in the asymmetric unit buries {N['lattice_area']:.0f} &#8491;&#178; against
{N['lattice_n']} residues, which is lattice packing next to a {N['iface_area']:.0f} &#8491;&#178;
interface.</p>
<p>The worked example sits inside the range on every quantity except one. Interface area and
paratope size match the other toxin complexes closely, and binding is loop-mediated in all four.
What is atypical is the CDR3 share: {N['cdr3_share']:.0f}% of the buried surface here against
{N['cdr3_other']:.0f}% for the comparator with the same antigen class. This nanobody binds
substantially through framework&nbsp;2 and CDR2 rather than predominantly through CDR3. That is a
property of this binder, not an artefact, and it widens the surface a mask has to cover.</p>
<p>Of the {N['n_paratope']} residues that bury surface on binding, {N['n_cdr_para']} lie in the
CDR loops and {N['n_fw_para']} in framework. Six of those framework residues
({N['fw2_para']}) sit in framework 2, the face that in a conventional antibody is buried against
the light chain. A VHH has no light chain, so that face is exposed and available to the antigen.
The consequence for the design is that a mask has to cover framework 2 as well as the loops, which
is more surface than a graft taken from the CDR contacts alone would address.</p>

<h3 class="mod"><span>Module 01</span>Estimation of the usable mask affinity</h3>
{mri("The Gly/Ser linker was treated as a worm-like chain with persistence length "
     f"{model['persistence_length_A']} Å and {model['rise_per_residue_A']} Å per residue, giving "
     "the effective concentration of the tethered mask at the distance measured above. An upper "
     "bound on the mask dissociation constant follows from requiring that little of the construct "
     "is unmasked before cleavage, and a lower bound from requiring that most of it binds antigen "
     "afterwards.",
     f"For an N-terminal attachment the effective concentration is {N['ceff']:.0f} mM and the "
     f"usable range of mask dissociation constants is {N['km_lo']:.0f} nM to {N['km_hi']:.0f} µM, "
     f"widest at a {N['linker_best']}-residue linker. A C-terminal attachment gives a narrower "
     f"range.",
     "Candidate masks should be designed towards micromolar rather than nanomolar affinity. The "
     "values depend on the assumptions of the polymer model and on the assumed antigen "
     "concentration, so they constrain the search rather than specify a target.")}
<div class="two">
{chart("design_window", "Usable mask affinity against linker length.")}
{chart("tradeoff", "The two requirements as functions of mask affinity.")}
</div>
<details class="tech"><summary>Model and equations</summary><div class="techbody">
<div class="eq">⟨r²⟩ = 2·l<sub>p</sub>·L<sub>c</sub>·[1 − (l<sub>p</sub>/L<sub>c</sub>)(1 − e<sup>−L<sub>c</sub>/l<sub>p</sub></sup>)]
<br>C<sub>eff</sub>(d) = 1660.54 · (3 / 2π⟨r²⟩)<sup>3/2</sup> · e<sup>−3d²/2⟨r²⟩</sup>
<br><br>P<sub>masked</sub> = C<sub>eff</sub> / (C<sub>eff</sub> + K<sub>D,M</sub>)
<br>P<sub>bound</sub> = (A/K<sub>D,A</sub>) / (1 + A/K<sub>D,A</sub> + M/K<sub>D,M</sub>)</div>
<p>The bounds correspond to at most {model['leak_max']*100:.0f}% of the construct unmasked before
cleavage and at least {model['recovery_min']*100:.0f}% antigen-bound afterwards.</p>
{chart("release_sensitivity", "Dependence of the lower bound on antigen affinity and concentration.")}
</div></details>

<h3 class="mod"><span>Module 02</span>Generation of candidate constructs</h3>
{mri("Masks were grafted from the antigen's own contact segments, with scrambled, truncated and "
     "poly-glycine controls. Linkers were not fixed sequences but generated: CleaveNet was "
     "conditioned on a profile that asks for substrates the bite-site proteases cut and the "
     "circulating ones do not, so cleavability constrains generation rather than filtering it "
     "afterwards. ProteinMPNN then redesigned the mask positions on each folded construct with "
     "the nanobody and linker held fixed.",
     f"{N['n_constructs']} constructs in the initial library and {N['n_redesign']} distinct "
     f"redesigned masks. Every redesign scored better than its parent graft by ProteinMPNN's own "
     f"measure, but that score correlated with paratope occlusion at r = {N['mpnn_corr']:+.2f}.",
     "ProteinMPNN generated sequence diversity, but sequence compatibility with a backbone is not "
     "a substitute for functional scoring and should not be used as one.")}
<div class="two">
{chart("redesign_score_vs_occlusion", "ProteinMPNN score against paratope occlusion.")}
<div class="scroll"><table class="wrap"><thead><tr><th>mask</th><th>sequence</th><th>kind</th>
</tr></thead><tbody>{masks_rows_compact}</tbody></table></div>
</div>

<h3 class="mod"><span>Module 03</span>Structural feasibility</h3>
{mri("Each construct was folded with ESMFold as one chain: mask, linker and nanobody "
     "together, not the mask alone. Preservation of the nanobody was assessed by "
     "superimposing its core on the separately folded free nanobody, and sequence-based "
     "descriptors were calculated over the whole construct for exposed hydrophobic "
     "surface and unpaired cysteines.",
     f"{N['n_folded']} constructs were folded and {N['n_ensemble']} superimposed on the nanobody. "
     f"The mean core RMSD against the free nanobody is {N['rmsd_mean']:.2f} Å. The mask grafted "
     f"from the antigen carries {N['feas_cys']} unpaired cysteines and scores "
     f"{N['feas_graft']:.2f} on this module.",
     "The nanobody fold is retained, which is a necessary condition rather than evidence of "
     "expression or thermodynamic stability. Cysteine content is the dominant liability and is "
     "detected by no other module.")}
{chart("feasibility_components", "Feasibility components by construct family. A low score on any "
       "one component is enough to make a construct unusable.")}
{R.ensemble_viewer("ens-feas", "feasibility",
   "Each construct is a single chain and was predicted on its own, so each was superimposed on "
   "the nanobody backbone to bring them into a common frame. The mask and linker carry the "
   "feasibility colour and the nanobody is drawn once, in teal, with a short run of it kept "
   "past each junction to show the chain is continuous.",
   R.ramp_legend("low", "high"))}
<details class="tech"><summary>Conformational heterogeneity</summary><div class="techbody">
<p>A single predicted structure does not describe the conformational distribution of a tethered
peptide. Ensemble methods such as BioEmu, or molecular dynamics, would be the appropriate
extension of this module and are not implemented here.</p></div></details>

<h3 class="mod"><span>Module 04</span>Mask binding and release</h3>
{mri("Two approaches were applied. Paratope occlusion was calculated from the folded constructs "
     "by comparing the accessible surface of the paratope with that of the free nanobody. The "
     "assembly of nanobody, tethered mask and antigen was then predicted with Boltz-2, with "
     "alignments supplied for the nanobody and the antigen.",
     f"Occlusion was {N['occ_graft']:+.1f}% for the grafted masks, {N['occ_scram']:+.1f}% for "
     f"scrambled controls and {N['occ_inert']:+.1f}% for the poly-glycine control. In the complex "
     f"predictions the antigen retained {N['boltz_lo']:.0f}–{N['boltz_hi']:.0f}% of its contacts "
     f"on the paratope in the masked constructs, against {N['boltz_ctrl']:.0f}% for the free "
     f"nanobody.",
     "Occlusion does not distinguish designed masks from controls and cannot be used to rank them. "
     "The complex prediction reproduces the known interface when given alignments, and under that "
     "measure none of the present masks excludes the antigen.")}
<div class="two">
{chart("occlusion_by_linker", "Paratope occlusion against linker length, by mask.")}
<div class="scroll"><table class="wrap"><thead><tr><th>construct</th>
<th>antigen contacts on paratope</th><th>mask contacts on paratope</th></tr></thead>
<tbody>{N['boltz_rows']}</tbody></table></div>
</div>
{R.ensemble_viewer("ens-mask", "masking",
   "The same constructs, coloured by paratope occlusion. The masks do not concentrate over the "
   "antigen-binding residues, which is the result reported above.",
   R.ramp_legend("least occluding", "most"))}
<div class="viewer">
<div class="vbar"><select id="pick">{options}</select>
<button data-spin="designed">spin</button></div>
<div class="vstage" id="stage-designed"></div>
<div class="vkey"><span><i style="background:#337D79"></i>framework</span>{R.CDR_KEY}
<span><i style="background:#8C72AD"></i>mask</span>
<span><i style="background:#9AA7B1"></i>linker</span>
<span><i style="background:#C06014"></i>antigen-binding residues</span></div>
<p class="dim" id="dcap" style="font-size:.85rem"></p>
</div>
<details class="tech"><summary>Complex prediction without alignments</summary><div class="techbody">
<p>When every chain was predicted without an alignment the free nanobody control placed only
{N['boltz_ctrl_noalign']:.0f}% of its antigen contacts on the known paratope. An alignment over the
whole fusion would require submitting the designed sequence to a remote server, so the construct
was modelled as three chains joined by an explicit covalent bond, allowing alignments for the
nanobody and the antigen while the designed peptide takes none.</p>
<p>Each case is a single prediction rather than a sampled ensemble, and the mask and antigen are
supplied simultaneously, so the result describes an assembly at equilibrium rather than competition
over time.</p></div></details>

<h3 class="mod"><span>Module 05</span>Linker cleavability</h3>
{mri("Two questions, kept apart. Sequence recognition is verified with CleaveNet's predictor and "
     "scored as a contrast, since a linker cut by everything would unmask in circulation. Reach is "
     "then checked in the folded construct against three necessary conditions: the segment lies "
     "extended enough for a catalytic groove, its backbone is exposed, and a sphere the size of a "
     "catalytic domain fits beside the scissile bond. Their product is taken, because a "
     "requirement that can be compensated for is not a requirement.",
     f"{N['n_linkers']} linkers were carried forward from {N['n_generated']} generated. "
     f"Selectivity spans {N['sel_min']:+.2f} to {N['sel_max']:+.2f}; the Gly/Ser control scores "
     f"{N['sel_ctrl']:+.2f}. The two halves are close to independent across the "
     f"{N['n_cleav']} linkered constructs (r = {N['sel_ctx_r']:+.2f}), so neither substitutes "
     f"for the other and the most selective sequence is not the most reachable.",
     "Generation and verification share a model family, so their agreement is weaker evidence "
     "than two independent models would be. CleaveNet is trained on human MMPs; the protease that "
     "arrives first at a bite site is a snake venom metalloproteinase, which is not in its "
     "training set, so the gelatinases stand in for it.")}
<p>The linkers carried forward, and the two controls they are read against. Target is the mean
predicted cleavage by {N['target_proteases']}; off-target is the mean over the membrane-type and
constitutive enzymes.</p>
{R.linker_table(LINKERS)}
<div class="two">
{chart("cleavability_panel", "Recognition against reach. The controls say whether either axis "
       "discriminates.")}
{chart("cleavability_components", "The three necessary conditions and their product, by construct.")}
</div>
{R.ensemble_viewer("ens-cleave", "cleavability",
   "The same constructs, coloured by how reachable their linker is in the folded model.",
   R.ramp_legend("least reachable", "most"))}
<p>Every linker is the same length, so length is held constant and only sequence varies. Reach in
the construct averages {N['ctx_mean']:.2f} and reaches {N['ctx_max']:.2f} at best; the affinity
window of Module 01 is widest at {N['linker_best']} residues, which the generated linkers sit
close to.</p>
<p><b>Both axes have to hold at once.</b> The lowest-scoring construct in the library,
{N['least_id']}, is physically reachable &mdash; its reach score of {N['least_ctx']:.2f} is
mid-pack &mdash; but its linker is a scramble of a good substrate, so recognition falls to
{N['least_sel']:+.2f} and the product with it, to {N['least_cleav']:.3f}. Same composition as
something cleavable, same accessibility, no site. That is the case a single averaged score would
have hidden.</p>

<h3 class="mod"><span>Module 06</span>Multiobjective ranking</h3>
{mri("The three per-construct scores were combined by identifying nondominated candidates. No "
     "weights were applied.",
     f"Of {N['rank_total']} scored constructs, {N['rank_front']} are nondominated.",
     "This demonstrates the ranking procedure. It does not identify functional switches, because "
     "the masking score does not separate designs from controls and the cleavability score covers "
     "accessibility only.")}
{chart("pareto_panels", "Demonstration of the ranking algorithm. The three scores plotted pairwise.")}

<h3 class="mod"><span>Summary</span>Implementation status</h3>
{R.implementation_table([
  ("01 Affinity estimation", "Interface geometry and tether model",
   "Validation against measured mask affinities"),
  ("02 Generation", "Candidate library and ProteinMPNN redesign",
   "De novo masking architectures"),
  ("03 Structural feasibility", "Static predictions and structural descriptors",
   "Conformational modelling; experimental expression data"),
  ("04 Mask binding", "Occlusion and complex prediction",
   "Reliable affinity and release prediction"),
  ("05 Linker cleavability", "Structural accessibility",
   "Protease-specific recognition and selectivity"),
  ("06 Ranking", "Pareto analysis on preliminary scores",
   "Ranking on validated functional scores"),
])}
{cautionary_html}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">05</div><h2>Experimental work</h2>
<p>I would work with experimental PhD students to establish measurements corresponding to the
three scoring functions: protein stability and developability, mask interaction and dissociation,
and protease responsiveness. Purification behaviour and gel analysis address the first,
biolayer interferometry in the presence and absence of the protease addresses the second and
third.</p>
<p>These data would allow the predictive ability of individual computational methods to be
assessed and would identify which parts of the pipeline require development. As more measurements
accumulate, I would investigate whether sequence and structural representations of the constructs
can be used to train more accurate scoring functions and improve the selection of subsequent
candidates.</p></div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">References</div><h2>Sources</h2></div>
<ul class="ref">
<li>Ahmadi S., Burlet N.J., Benard-Valle M. <i>et al.</i> Nanobody-based recombinant antivenom for
cobra, mamba and rinkhals bites. <i>Nature</i> 647:716–725 (2025). Source of
<a href="https://www.rcsb.org/structure/9RIU">PDB 9RIU</a>, used here.</li>
<li>Damsbo A., Burlet N.J., Fernández-Quintero M.L. <i>et al.</i> Structural mechanisms behind the
neutralisation of long-chain α-neurotoxins. <i>Communications Chemistry</i> (2025). Source of
<a href="https://www.rcsb.org/structure/9GCN">PDB 9GCN</a>.</li>
<li>Lin Z. <i>et al.</i> Evolutionary-scale prediction of atomic-level protein structure with a
language model. <i>Science</i> 379:1123–1130 (2023).</li>
<li>Dauparas J. <i>et al.</i> Robust deep learning-based protein sequence design using ProteinMPNN.
<i>Science</i> 378:49–56 (2022).</li>
<li>Wohlwend J. <i>et al.</i> Boltz-1 and Boltz-2, open biomolecular interaction models.
<a href="https://github.com/jwohlwend/boltz">github.com/jwohlwend/boltz</a></li>
<li>Zhou H.-X. Loops in proteins can be modeled as worm-like chains. <i>Journal of Physical
Chemistry B</i> 105:6763–6766 (2001).</li>
<li>Shrake A., Rupley J.A. Environment and exposure to solvent of protein atoms. <i>Journal of
Molecular Biology</i> 79:351–371 (1973).</li>
<li>Stepurko N. Molecular insights into signalling mechanisms within the TGF-β superfamily.
PhD thesis, University of Cambridge (2023).
<a href="https://doi.org/10.17863/CAM.99087">10.17863/CAM.99087</a></li>
</ul>
</section>

<footer>Figures and values generated by the analysis code. Structural data from the RCSB Protein
Data Bank.</footer>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
<script>{R.VIEWER_JS.replace("__MODELS__", json.dumps(models)).replace("__PARATOPE__", json.dumps(para_resi)).replace("__SOLVED__", json.dumps(f"structures/{PDB_ID}_{REFERENCE_VHH_CHAIN}{iface['toxin_chain']}.pdb")).replace("__ENSEMBLE__", json.dumps(ENSEMBLE_IDS)).replace("__TOPPARA__", json.dumps(top_paratope)).replace("__CDRS__", json.dumps(cdr_resi)).replace("__SCORES__", json.dumps(SCORES))}</script>
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
