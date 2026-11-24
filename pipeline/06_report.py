#!/usr/bin/env python
"""Stage 6. Render docs/index.html from the stage outputs."""
import sys, json, os
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
    """Constructs whose folded model does not match the library.

    Refuse to build a page whose numbers and structures disagree.

    The construct library and the folded coordinates are separate artefacts.
    Changing the library without refolding leaves residue ranges pointing at
    the wrong atoms, which renders as a page that looks fine and is wrong.
    """
    from nbswitch.config import STRUCTURES
    lib = load_json("03_constructs.json")
    stale = []
    for item in lib["constructs"]:
        path = STRUCTURES / f"{item['id']}.pdb"
        if not path.exists():
            continue
        last = 0
        for line in path.read_text().splitlines():
            if line.startswith("ATOM"):
                last = max(last, int(line[22:26]))
        if last and last != item["vhh_range"][1]:
            stale.append(item["id"])
    return stale


#: Building while a refold is in flight mixes new design stages with scores
#: computed on the previous library. That is publishable only if the page says
#: so, so the guard is a hard stop unless NBSWITCH_ALLOW_STALE is set, and then
#: it reports instead of refusing.
STALE = _check_library_matches_structures()
ALLOW_STALE = os.environ.get("NBSWITCH_ALLOW_STALE") == "1"
if STALE and not ALLOW_STALE:
    raise SystemExit(
        f"library and structures disagree for {len(STALE)} constructs "
        f"({', '.join(STALE[:3])}...). Refold stages 04 and 08, or set "
        f"NBSWITCH_ALLOW_STALE=1 to publish with the mismatch declared on the page.")
tether = load_json("02_tether.json")
lib = load_json("03_constructs.json")
try:
    fold = load_json("04_folding.json")
except FileNotFoundError:
    fold = None
for _name, _key in (("02b_linkers.json", "linkers"),
                    ("13_ranking.json", "ranking_data"),
                    ("12_feasibility.json", "feas_data"),
                    ("19_verification.json", "verification"),
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
_c_in_paratope = geom["c_term_residue"] in set(iface["paratope_labels"])
#: whether a mask could be hung off either end. 9GCN fails this, which is why
#: it is a set-aside case, so the claim is derived rather than asserted.
termini_free = (
    "neither terminus is part of it, so either could in principle carry a mask. Which one does "
    "is decided by the distances below"
    if not (_n_in_paratope or _c_in_paratope) else
    ("the N-terminal residue is itself part of it, which is why this complex is a set-aside case"
     if _n_in_paratope else
     "the C-terminal residue is itself part of it, so only the N-terminus can carry a mask"))
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
def chart(n, c):
    if not (DOCS / "figures" / f"{n}.svg").exists():
        return R.pending(f"the figure <span class='mono'>{n}</span>")
    return R._chart(n, c)
def chart_wide(n, c):
    if not (DOCS / "figures" / f"{n}.svg").exists():
        return R.pending(f"the figure <span class='mono'>{n}</span>")
    return R._chart(n, c).replace('class="chart"', 'class="chart wide"')


# ---- redesign section ----
# names this block defines that the live page reads, defined first so that a
# stage which has not run leaves a missing value rather than a NameError
corr = float("nan")
rrows, rocc, grafts, improved = [], [], [], 0
not_improved, not_improved_parents = [], []
best_native_parent, best_native_score = "", float("nan")
rbest = None

if redesign_fold and fold:
    rrows = redesign_fold["results"]
    grafts = [r["occlusion"] * 100 for r in fold["results"] if r["kind"] == "epitope graft"]
    rocc = [r["occlusion"] * 100 for r in rrows]
    rbest = max(rrows, key=lambda r: r["occlusion"])
    improved = sum(1 for c in redesign["constructs"]
                   if c["native_score"] and c["mpnn_score"] < c["native_score"])
    # which parents ProteinMPNN could not improve on, and by how much they
    # already beat it. The page claimed every redesign improved its parent;
    # the stage output says 88 of 96, and the 8 share one parent.
    _not_improved = [c for c in redesign["constructs"]
                     if c["native_score"] and c["mpnn_score"] >= c["native_score"]]
    not_improved = _not_improved
    not_improved_parents = sorted({c["parent"] for c in _not_improved})
    best_native_parent = min(redesign["native_scores"],
                             key=redesign["native_scores"].get)
    best_native_score = redesign["native_scores"][best_native_parent]
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
<p>Of {ranking_data['n_total']} scored designs, {ranking_data.get('n_excluded', 0)} are excluded
by hard constraints before any comparison: a nanobody that is not confidently predicted, a fold
that has moved, a linker cut no more at the site than in circulation, or one a protease cannot
reach. None of those can be compensated for by another score, which is what a front would
otherwise allow. Of the {ranking_data.get('n_feasible', 0)} that remain,
{ranking_data['n_front']} survive this test. The
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
    _m = {r["case"]: r for r in (verification["results"] if verification else [])
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
if verification:
    _msa_ctrl = next((r for r in verification["results"]
                      if r["case"].startswith("control") and r["status"] == "ok"), None)
_msa_masked = [r for r in (verification["results"] if verification else [])
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
N["window_fold"] = 10 ** bN["window_decades"]
N["linker_best"] = bN["linker_residues"]
N["km_lo"] = bN["km_min_for_release_M"] * 1e9
N["km_hi"] = bN["km_max_for_masking_M"] * 1e6
N["ceff"] = bN["c_eff_M"] * 1e3
#: the tether model spans the paratope centroid, not the nearest paratope atom,
#: so both distances are reported: the nearest-atom pair says which terminus is
#: closer, the centroid pair is what C_eff is actually evaluated at.
N["span_n"] = geom["n_term_to_paratope_centroid_A"]
N["span_c"] = geom["c_term_to_paratope_centroid_A"]
# one linker length worked through in the prose, computed from the same model
# the figure uses so the worked example cannot drift from the curve
import math as _math
from nbswitch import tether as _tether
N["ex_n"] = 14
N["ex_lc"] = N["ex_n"] * model["rise_per_residue_A"]
N["ex_ree"] = _math.sqrt(_tether.mean_square_end_to_end(N["ex_n"]))
N["ex_ceff"] = _tether.effective_concentration(N["ex_n"], N["span_n"]) * 1e3

# the two Module 01 equations, typeset rather than written in running text
from nbswitch import mathfig as _mf
EQ_WLC = _mf.equation(
    r"\langle r^2 \rangle = 2 l_p L_c "
    r"\left[ 1 - \frac{l_p}{L_c}\left(1 - e^{-L_c/l_p}\right) \right]",
    "Mean square end-to-end distance equals two l p L c, times one minus "
    "l p over L c, times one minus e to the minus L c over l p.",
    "wlc")
EQ_CEFF = _mf.equation(
    r"C_{\mathrm{eff}}(d) = \frac{10^{27}}{N_A} "
    r"\left( \frac{3}{2\pi\langle r^2\rangle} \right)^{3/2} \, "
    r"e^{-3d^2/(2\langle r^2\rangle)}",
    "C eff of d equals ten to the twenty-seventh over Avogadro's number, "
    "times three over two pi mean square end-to-end distance, to the power "
    "three halves, times e to the minus three d squared over twice the mean "
    "square end-to-end distance.",
    "ceff")

N["ree_best"] = bN["mean_end_to_end_A"]
N["linker_best_c"] = bC["linker_residues"]
N["window_c"] = bC["window_decades"]
N["n_constructs"] = lib["n_constructs"]
N["n_folded"] = fold["n_folded"] if fold else 0

_k = {}
for _r in (fold["results"] if fold else []):
    _k.setdefault(_r["kind"], []).append(_r["occlusion"] * 100)
N["occ_scram"] = _np.mean(_k.get("scrambled control", [0]))
N["occ_inert"] = _np.mean(_k.get("inert control", [0]))

N["n_redesign"] = redesign["n_unique"] if redesign else 0
# the figure plots every construct carrying a mask, so its correlation is over
# all of them and not only the redesigns the earlier version drew
LIBSC = None
try:
    LIBSC = load_json("07b_library_scores.json")
except FileNotFoundError:
    pass
ALL_CORR, ALL_N, ALL_SD = float("nan"), 0, float("nan")
if LIBSC and fold:
    _occ = {r["id"]: r["occlusion"] * 100
            for src in (fold, redesign_fold) for r in (src or {}).get("results", [])}
    _pairs = [(v["mask_score"], _occ[k]) for k, v in LIBSC["scores"].items() if k in _occ]
    if len(_pairs) > 2:
        ALL_CORR = float(np.corrcoef([p[0] for p in _pairs], [p[1] for p in _pairs])[0, 1])
        ALL_N = len(_pairs)
    ALL_SD = float(np.mean([v["mask_score_sd"] for v in LIBSC["scores"].values()]))

N["mpnn_corr"] = corr if redesign else float("nan")
N["n_improved"] = improved
N["n_not_improved"] = len(not_improved)
N["not_improved_parents"] = not_improved_parents
N["best_native_parent"] = best_native_parent
N["best_native_score"] = best_native_score

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

_mm = {r["case"]: r for r in (verification["results"] if verification else [])
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

# the set the mask is designed from, shown on both partners
top_paratope = [
    {"resi": int("".join(c for c in lbl[1:] if c.isdigit())), "label": lbl}
    for lbl in iface.get("paratope_near_labels", iface["paratope_contact_labels"])]
toxin_contact = sorted(int("".join(c for c in lbl[1:] if c.isdigit()))
                       for lbl in iface.get("epitope_near_labels", []))
N["mask_cutoff"] = iface.get("mask_contact_cutoff_A", 3.5)
N["n_near"] = len(top_paratope)
N["n_epitope_near"] = len(iface.get("epitope_near_labels", []))
N["n_epitope_contacts"] = iface.get("n_epitope_contacts", 0)
N["n_epitope"] = iface.get("n_epitope", 0)
try:
    _seg = load_json("03a_masks.json")["segments"][0]
    N["mask_face"] = f"{_seg['start_label']}–{_seg['end_label']}"
except (FileNotFoundError, KeyError, IndexError):
    N["mask_face"] = "not yet designed"
N["n_contacts"] = iface["n_paratope_contacts"]
N["contact_cutoff"] = iface["contact_cutoff_A"]

# what a pH trigger would have to work with, computed rather than asserted
from nbswitch import triggers as _trig
PH = _trig.switch_window()
PH_SCAN = _trig.sensitivity()

_ref_feas = next((r for r in (feas_data["results"] if feas_data else [])
                  if r["mask"] is None), None)
N["ref_surface"] = _ref_feas["f_surface"] if _ref_feas else float("nan")
from nbswitch.config import VHH_TOXIN_KD_M as _VHH_KD
N["vhh_toxin_kd_nM"] = _VHH_KD * 1e9
N["contact_labels"] = ", ".join(iface["paratope_contact_labels"])

SCORES = {}
_feas_by = {r["id"]: r["feasibility_score"] for r in (feas_data["results"] if feas_data else [])}
_occ_by = {}
for _src in (fold, redesign_fold):
    for _r in (_src["results"] if _src else []):
        _occ_by[_r["id"]] = _r["occlusion"]
# the viewer's "cleavability" ramp held the context score, which is reach in
# the fold and only half of what the page calls cleavability. It is the score
# the ranking uses: recognition multiplied by accessibility.
_cleav_by_id = {r["id"]: r["cleavability_score"]
                for r in (trigger["results"] if trigger else [])}
for _key, _d in (("feasibility", _feas_by), ("masking", _occ_by), ("cleavability", _cleav_by_id)):
    _n = _norm({it["id"]: _d.get(it["id"]) for it in ENSEMBLE_IDS})
    for _i, _v in _n.items():
        SCORES.setdefault(_i, {})[_key] = _v

# the page must build from whatever has finished, so that its structure can
# be reviewed while the expensive stages are still running
try:
    controls = load_json("15_interface_controls.json")
except FileNotFoundError:
    controls = None
_by = {r["pdb_id"]: r for r in (controls["comparators"] if controls else [])}
N["cdr3_share"] = _by.get(PDB_ID, {}).get("cdr3_percent_bsa", 0)
N["iface_area"] = iface["interface_area_A2"]
_lat = _by.get(PDB_ID, {}).get("lattice_contact", {})
N["lattice_area"] = _lat.get("interface_area_A2", 0)
N["lattice_n"] = _lat.get("n_paratope", 0)
N["cdr3_other"] = _by.get("9GCN", {}).get("cdr3_percent_bsa", 0)
def _count_rows():
    """What each module produced, labelled by the module in the diagram.

    One result column, not an in/out pair. Most steps are not selections:
    masks are built from contact residues, redesigns from parents, and the
    units differ on each side, so an arrow between two counts said something
    untrue about every row except the two that really do discard candidates.
    Those two say so in words.
    """
    rows = []

    def row(module, what, result, feeds, selected=False):
        rows.append({"stage": module, "what": what, "result": result,
                     "feeds": feeds, "selected": selected})

    try:
        teth = load_json("02_tether.json")
        row("01", "the deposited complex measured: paratope, termini and "
                  "contact faces, then a tether model applied to the geometry",
            f"a {teth['best_per_terminus']['N']['window_decades']:.2f}-decade "
            f"window on mask affinity",
            "02a and 02b")
    except (FileNotFoundError, KeyError):
        pass

    try:
        lk = load_json("02b_linkers.json")
        row("02b", "protease substrates generated against the target profile "
                   "and scored for selectivity",
            f"{lk['n_generated']} generated, <b>{lk['n_kept']} kept</b>",
            "the linker set", selected=True)
        row("02b", "controls added, and the kept substrates rebuilt at each "
                   "modelled linker length",
            f"<b>{len(lk['linkers'])} linkers</b>", "assembly")
    except FileNotFoundError:
        pass

    try:
        mk = load_json("03a_masks.json")
        n_contact = sum(seg["n_contact"] for seg in mk["segments"])
        row("02a", f"masks built from the {n_contact} antigen residues within "
                   f"{mk['cutoff_A']} \u00c5 of the nanobody: window and "
                   f"substitution variants, and matched controls",
            f"<b>{len(mk['masks'])} masks</b>", "assembly")
    except FileNotFoundError:
        pass

    try:
        cl = load_json("03_constructs.json")
        n_mask, n_link = len(cl["masks"]), cl["n_linkers"]
        row("\u2014", "every mask paired with every linker, at every length, so no "
                      "comparison rests on a combination that was not built",
            f"{n_mask} \u00d7 {n_link} = {n_mask * n_link} pairs, plus the free "
            f"nanobody = <b>{cl['n_constructs']}</b>",
            "structure prediction")
        f = load_json("04_folding.json")
        row("\u2014", "each construct folded once with ESMFold, as a single chain",
            f"<b>{f.get('n_folded', 0)} structures</b>", "03, 04 and 05")
    except FileNotFoundError:
        pass

    try:
        rd = load_json("07_redesign.json")
        n_par = len(rd.get("parents", [])) or 1
        row("02a", f"the mask rewritten by ProteinMPNN on {n_par} of those "
                   f"folded backbones, nanobody and linker held fixed",
            f"<b>{len(rd['constructs'])} designs</b>", "structure prediction")
        rf = load_json("08_redesign_folding.json")
        row("\u2014", "redesigns folded with ESMFold",
            f"<b>{rf.get('n_folded') or len(rf['results'])} structures</b>",
            "03, 04 and 05")
    except FileNotFoundError:
        pass

    try:
        n = len(load_json("12_feasibility.json")["results"])
        row("03", "scored for fold retention, exposed surface and disulfides",
            f"<b>{n} scored</b>", "06")
    except FileNotFoundError:
        pass
    try:
        n = len(load_json("10_trigger.json")["results"])
        row("05", "scored for protease recognition and for reach in the fold",
            f"<b>{n} scored</b>", "06")
    except FileNotFoundError:
        pass
    try:
        n = len(load_json("04_folding.json")["results"])
        row("04", "paratope surface occlusion measured against the free nanobody",
            f"<b>{n} scored</b>", "06")
    except FileNotFoundError:
        pass
    try:
        mb = load_json("17_mask_binding.json")
        ok = [r for r in mb["results"] if r.get("status") == "ok"]
        row("04", "each mask predicted against the bare nanobody, untethered, "
                  "with the antigen as a positive control",
            f"<b>{len(ok)} predictions</b>", "read alongside 06")
    except FileNotFoundError:
        pass

    try:
        rk = load_json("13_ranking.json")
        row("06", "hard constraints, then the nondominated set, then ordered "
                  "by weakest objective and deduplicated by mask",
            f"{rk.get('n_total', 0)} scored, {rk.get('n_excluded', 0)} ruled out, "
            f"{rk.get('n_front', 0)} nondominated, "
            f"<b>{len(rk.get('shortlist', []))} shortlisted</b>",
            "04 verification", selected=True)
    except (FileNotFoundError, KeyError):
        pass
    try:
        cx = load_json("19_verification.json")
        ok = [r for r in cx["results"] if r.get("status") == "ok" and r.get("construct")]
        worth = cx.get("n_worth_testing")
        row("04", "the shortlist predicted as a whole assembly with its antigen",
            f"{len(ok)} predicted, <b>{worth if worth is not None else len(ok)} "
            f"worth making</b>", "the reported candidates", selected=True)
    except FileNotFoundError:
        pass
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
except FileNotFoundError:
    N.update({"n_linkers": 0, "n_generated": 0, "sel_min": 0.0, "sel_max": 0.0,
              "sel_ctrl": 0.0, "target_proteases": "", "off_target_proteases": ""})


COUNTS = _count_rows()


N["rank_total"] = ranking_data["n_total"] if ranking_data else 0
N["rank_front"] = ranking_data["n_front"] if ranking_data else 0
RANK_FRONT = ranking_data.get("front", []) if ranking_data else []
RANK_OBJECTIVES = ranking_data.get("objectives", []) if ranking_data else []
RANK_SHORTLIST = set(ranking_data.get("shortlist", [])) if ranking_data else set()
N["rank_shortlist"] = len(RANK_SHORTLIST)
# the shortlist in front order, and the rest of the front behind it, so the
# page never reorders what stage 13 decided
RANK_SHORTLIST_ROWS = [r for r in RANK_FRONT if r["id"] in RANK_SHORTLIST]
RANK_FRONT_REST = [r for r in RANK_FRONT if r["id"] not in RANK_SHORTLIST]
if RANK_FRONT and len(RANK_SHORTLIST_ROWS) != len(RANK_SHORTLIST):
    raise SystemExit("a shortlisted construct is missing from the front; "
                     "stage 13's output is inconsistent")

_feas = {r["id"]: r for r in (feas_data["results"] if feas_data else [])}
# the primary graft, whatever it is called: masks are renamed when the
# design changes, and a literal name here silently selects nothing
_primary_mask = next((k for k, m in lib["masks"].items()
                      if m["kind"] == "epitope graft"), None)
_graft_feas = [r for r in _feas.values() if r["mask"] == _primary_mask]
# no rows for the current graft means the scores predate it, which must not
# render as a score of zero
N["has_graft_feas"] = bool(_graft_feas)
N["feas_graft"] = _np.mean([r["feasibility_score"] for r in _graft_feas]) if _graft_feas else float("nan")
N["feas_cys"] = _graft_feas[0]["unpaired_cysteines"] if _graft_feas else -1
N["rmsd_mean"] = _np.mean([r["vhh_core_rmsd_A"] for r in _feas.values()
                           if r.get("vhh_core_rmsd_A")]) if _feas else 0


masks_rows_compact = "".join(
    f"<tr><td>{k}</td><td class='mono'>{m['seq']}</td><td>{m['kind']}</td></tr>"
    for k, m in lib["masks"].items())

_LK = load_json("02b_linkers.json")
# the key is off_target_proteases. Reading "off_target" found nothing, and the
# `or 15` behind it then published a number no stage had produced: 15 is every
# protease that is not a target, but only 9 are pushed down. The other 6 are
# left at zero, which is a different statement. Indexed, not .get, so a rename
# raises here instead of substituting a plausible-looking constant.
#: stage 02b records which proteases it asked for and which it asked against,
#: but not the full set CleaveNet scores, so the total comes from the module
#: that built the profile rather than from a literal.
from nbswitch.linkers import PROTEASES as _PROTEASES, conditioning_profile
_TGT, _OFF = _LK["target_proteases"], _LK["off_target_proteases"]
_NEUTRAL = [p for p in _PROTEASES if p not in _TGT and p not in _OFF]
_PROFILE = conditioning_profile().iloc[0]
LINKER_GEN = {"n_generated": _LK["n_generated"], "n_kept": _LK["n_kept"],
              "n_controls": _LK["n_controls"],
              "n_target": len(_TGT), "n_off": len(_OFF),
              "n_neutral": len(_NEUTRAL),
              "target_names": _TGT, "off_names": _OFF,
              "neutral_names": _NEUTRAL,
              # the z-scores the profile actually carries, not restated here
              "target_z": float(_PROFILE[_TGT[0]]),
              "off_z": float(_PROFILE[_OFF[0]])}
# how the library was folded, read from the folding stage rather than named
FOLD_MODEL = fold["model"] if fold else "facebook/esmfold_v1"
FREE_PLDDT = fold["free_vhh_plddt"] if fold else float("nan")
FREE_SASA = fold["free_vhh_paratope_sasa_A2"] if fold else float("nan")

# ---- the mask side: what step 2 of module 02 is built from
_MD = load_json("03a_masks.json")
MASKS, HOTSPOTS = _MD["masks"], {int(l[1:]) for l in _MD["hotspot_labels"]}
CORE_ID = "epi_core"
_CORE = MASKS[CORE_ID]
CORE_SEQ, CORE_LEN = _CORE["seq"], len(_CORE["seq"])
CORE_RANGE = f"{_CORE['seq'][0]}{_CORE['start']}–{_CORE['seq'][-1]}{_CORE['end']}"
N["n_masks"] = len(lib["masks"])
# the construct the schematic draws: the reference mask on the first CleaveNet
# linker, which is the parent every redesign in step 3 came from
REF_CONSTRUCT = next(c for c in lib["constructs"] if c["id"] == f"{CORE_ID}_cn01")
# the worked-example funnel, every count read from the stage that produced it
_RANK = load_json("13_ranking.json")
try:
    _VER = load_json("19_verification.json")
except FileNotFoundError:
    _VER = {"n_worth_testing": 0}
# the 21 linkers are six sequences at the reference length plus length
# variants of three of them, not four motifs and two controls repeated
_REF_LINKERS = [x for x in _LK["linkers"] if "_n" not in x["id"]]
_LEN_VARIANTS = [x for x in _LK["linkers"] if "_n" in x["id"]]
# the library file carries the free nanobody alongside the constructs, so the
# construct count is one more than the mask x linker product
_N_REDESIGNS = redesign["n_unique"] if redesign else 0
_N_PARENTS = len(redesign["parents"]) if redesign else 0

FUNNEL = {
    "n_masks": len(lib["masks"]), "n_linkers": lib["n_linkers"],
    "n_pairs": lib["n_pairs"], "n_motifs": _LK["n_kept"],
    "n_linker_controls": _LK["n_controls"], "n_generated": _LK["n_generated"],
    "n_ref_linkers": len(_REF_LINKERS), "n_len_variants": len(_LEN_VARIANTS),
    "ref_linker_len": len(_REF_LINKERS[0]["sequence"]) if _REF_LINKERS else 0,
    "n_redesigns": _N_REDESIGNS, "n_parents": _N_PARENTS,
    "n_per_parent": _N_REDESIGNS // _N_PARENTS if _N_PARENTS else 0,
    # every construct folded once, plus the free nanobody reference
    "n_first_round": lib["n_constructs"],
    "n_structures": lib["n_constructs"] + _N_REDESIGNS,
    "n_evaluated": _RANK["n_total"], "n_front": _RANK["n_front"],
    "n_shortlist": len(_RANK["shortlist"]), "n_worth": _VER["n_worth_testing"],
}
if FUNNEL["n_pairs"] + FUNNEL["n_redesigns"] != FUNNEL["n_evaluated"]:
    raise SystemExit(
        f"the funnel does not add up: {FUNNEL['n_pairs']} pairs + "
        f"{FUNNEL['n_redesigns']} redesigns is not {FUNNEL['n_evaluated']} "
        "scored candidates")
N["n_evaluated"] = FUNNEL["n_evaluated"]
N["n_structures"] = FUNNEL["n_structures"]

# ---- module 04: the verification table and the numbers around it
_VER = load_json("19_verification.json")
_VROWS = _VER["results"]
_CTL = next(r for r in _VROWS if r["construct"] is None)
# the matched comparison: one antigen-derived mask against the three controls
# built from it, each over the same linker set
_CONTROL_Q = {
    "epi_core": ("the antigen's main contact loop, grafted",
                 "does reproducing a contact surface help cover the paratope?"),
    "epi_core_scram": ("the same residues in a different order",
                       "does the order matter, or only the composition?"),
    "epi_core_half": ("the first half of that loop",
                      "does the whole contact region beat a fragment of it?"),
    "polyG": ("twelve glycines, matched in length",
              "how much coverage comes from a flexible peptide alone?"),
}
_occ_by_mask = {}
for _r in (fold["results"] if fold else []):
    _occ_by_mask.setdefault(_r["mask"], []).append(_r["occlusion"] * 100)
def _occ_row(mask, what, question):
    cls = " class='hl'" if mask == CORE_ID else ""
    mean = _np.mean(_occ_by_mask.get(mask, [0]))
    return (f"<tr{cls}><td class='mono'>{mask}</td><td>{what}</td>"
            f"<td>{question}</td><td><b>{mean:.1f}%</b></td></tr>")

OCC_CONTROL_ROWS = "".join(_occ_row(m, w, q) for m, (w, q) in _CONTROL_Q.items())

N["occ_trunc"] = float(_np.mean([r["occlusion"] * 100 for r in fold["results"]
                                 if r["kind"] == "truncated control"])) if fold else 0.0

def _verdict(r):
    """Why a construct passed or failed, in the terms the criteria are in."""
    if r["construct"] is None:
        return "reference"
    displaced = r["antigen_contacts_displaced"] >= _VER["blocking_margin"]
    on_site = r["mask_fraction_on_paratope"] >= _VER["mask_on_paratope_min"]
    if displaced and on_site:
        return "<b>antigen displaced, the mask is bound</b>"
    if displaced:
        return "antigen displaced, but the mask is off site"
    return "antigen keeps its binding site"

_ordered = [_CTL] + sorted((r for r in _VROWS if r["construct"]),
                           key=lambda r: r["fraction_of_contacts_on_paratope"])
VERIFY_ROWS = "".join(
    "<tr{cls}><td class='mono'>{name}</td><td>{ag:.0f}%</td><td>{mk}</td>"
    "<td>{why}</td></tr>".format(
        cls=" class='hl'" if r.get("verdict") == "worth testing" else "",
        name=r["construct"] or "free nanobody",
        ag=r["fraction_of_contacts_on_paratope"] * 100,
        mk=("&mdash;" if r["construct"] is None
            else f"{r['mask_fraction_on_paratope'] * 100:.0f}%"),
        why=_verdict(r))
    for r in _ordered)
# ---- module 05: the two requirements, and the constants behind the second
from nbswitch.cleavage_context import (CATALYTIC_DOMAIN_RADIUS_A as PROBE_R,
                                       N_PROBE_DIRECTIONS as PROBE_N,
                                       PROBE_OFFSET_A as PROBE_OFF,
                                       CLEFT_HALF_WIDTH as CLEFT_HALF,
                                       EXTENDED_RISE_A as EXT_RISE)
N["len_lo"] = min(len(l["sequence"]) for l in _LK["linkers"])
N["len_hi"] = max(len(l["sequence"]) for l in _LK["linkers"])

# the feasibility stage's own constants and the family it scores worst, so
# the module text and its method panel cannot drift from the code
from nbswitch import feasibility as _FEAS
_FZ = load_json("12_feasibility.json")
FEAS = {"plddt_ref": 90, "rmsd_ref": 3, "phobic_ref": _FEAS.HYDROPHOBIC_REFERENCE,
        "cys_penalty": _FEAS.CYS_PENALTY, "cys_floor": _FEAS.CYS_FLOOR,
        "ss_max": _FEAS.DISULFIDE_MAX_A}
N["feas_n"] = _FZ["n_scored"]

# the final screen's two structures, already in one frame by stage 19b
try:
    _FP = load_json("19b_final_pair.json")
except FileNotFoundError:
    _FP = None
FINAL_PAIR_JS = None
if _FP:
    def _roles(chains, has_mask):
        """A Boltz prediction names chains in input order: nanobody, then the
        designed peptide if there is one, then the antigen."""
        names = sorted(chains)
        return ({"vhh": names[0], "mask": names[1], "antigen": names[2]}
                if has_mask and len(names) > 2 else
                {"vhh": names[0], "mask": None, "antigen": names[-1]})
    FINAL_PAIR_JS = {
        "reference": {"file": _FP["reference"]["file"],
                      **_roles(_FP["reference"]["chains"], False)},
        "masked": {"file": _FP["masked"]["file"],
                   **_roles(_FP["masked"]["chains"], True)},
    }
FP_ID = _FP["masked"]["construct"] if _FP else ""
FP_RMSD = _FP["vhh_ca_rmsd_A"] if _FP else float("nan")
FP_MASK_ON = (_FP["masked"]["mask_fraction_on_paratope"] * 100) if _FP else 0
FP_ANTIGEN_ON = (_FP["masked"]["fraction_of_contacts_on_paratope"] * 100) if _FP else 0
FP_REF_ON = (_FP["reference"]["fraction_of_contacts_on_paratope"] * 100) if _FP else 0
N_PARENTS = len(redesign["parents"]) if redesign else 0

# the contacting antigen residues are reported as runs along the chain, and
# for this complex they are not one run. Phrase it from the segments rather
# than describing the shape they happen to have.
_SEGS = sorted(_MD["segments"], key=lambda g: -g["n_contact"])
_big = _SEGS[0]
_rest = sum(g["n_contact"] for g in _SEGS[1:])
SEG_SPLIT = (f"{R.word(_big['n_contact']).capitalize()} lie in a single run, "
             f"{_big['start_label']}–{_big['end_label']}, and "
             + (f"one sits apart at {_SEGS[1]['start_label']}"
                if _rest == 1 else
                f"{R.word(_rest)} sit apart in {R.word(len(_SEGS) - 1)} shorter runs")
             + ".")

# the extensions pick up cysteines bonded to partners outside the graft. Read
# the partners off the structure-derived report so the claim cannot drift.
_EXT = MASKS.get("epi_core_ext2", {})
_ORPH = _EXT.get("cysteines_orphaned") or []
_PART = sorted({v for v in (_EXT.get("orphan_partners") or {}).values()})
ORPHAN_TEXT = (f"{R.word(len(_ORPH))} more cysteines whose partners, "
               + " and ".join(f"C{v}" for v in _PART)
               + ", stay behind in the antigen"
               if _ORPH else "no further cysteines")

CLEAVE_TARGETS = ", ".join(f"MMP-{p[3:]}" for p in
                           sorted(_TGT, key=lambda p: int(p[3:])))
CLEAVE_TARGETS = " and ".join(CLEAVE_TARGETS.rsplit(", ", 1))
_BY_ID = {l["id"]: l for l in _LK["linkers"]}
CN_BEST = max(l["selectivity"] for l in _LK["linkers"]
              if l.get("kind") == "generated")
CTRL_GS = _BY_ID["ctrl_gs"]["selectivity"]
CTRL_SCRAM = _BY_ID["ctrl_scram"]["selectivity"]
# the claim is that the scramble holds composition fixed, so check it rather
# than describe it: a scramble that is not one would make the control vacuous.
_BEST_ID = max((l for l in _LK["linkers"] if l.get("kind") == "generated"),
               key=lambda l: l["selectivity"])["substrate"]
assert sorted(_BY_ID["ctrl_scram"]["substrate"]) == sorted(_BEST_ID), \
    "the scrambled control is not a permutation of the best generated substrate"
# the reference-length linkers are not length variants and carry no
# total_length, so taking that field alone drops the reference from the series
LENGTHS = sorted({len(x["sequence"]) for x in _LK["linkers"]})
SUB_LEN = len(_LK["linkers"][0]["substrate"])
REF_LEN = lib["reference_length"]

HTML = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Computational design of environment-responsive nanobodies</title>
<meta name="description" content="A computational approach to designing nanobodies whose antigen
binding is restored by proteolysis, with a partial implementation on an experimental complex.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;650&display=swap" rel="stylesheet">
<style>{R.CSS}</style></head><body><div class="page">

<header>
<div class="kicker"><span>Natalija Stepurko</span><b>DTU · postdoctoral interview</b></div>
<h1>Computational design of environment-responsive nanobodies</h1>
<div class="banner"><b>Environment-responsive nanobodies</b><br>
The bite site of a viper envenomation is marked by several local cues absent from systemic
circulation: acidic pH (~6.0&ndash;6.5), elevated protease activity (tissue-damaged, MMPs/venom
proteases), and localized hypoxia. Propose a computational pipeline to design a nanobody-based
&ldquo;conditional switch&rdquo; that stays masked systemically but activates at the bite site
using one or more of these cues. Consider whether there would be possibilities for making an
integrated universal switch and/or robust standardized methods. Also consider how you would
collaborate with wet lab PhD students for experimental testing and validation.</div>
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
{R.rebuild_banner(STALE, len(lib['constructs']))}
{background_html}
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">02</div><h2>Research question and molecular architecture</h2>
<p>The question concerned nanobodies against snake venom toxins that would remain inactive in
circulation and become active at the site of a bite, where the local environment differs from
plasma. Three environmental cues were proposed: acidic pH, raised protease activity and local
hypoxia.</p></div>

<h4 class="subh">Choosing an environmental trigger</h4>
<div class="scroll"><table class="wrap">
<colgroup><col style="width:13%"><col style="width:30%"><col style="width:44%">
<col style="width:13%"></colgroup>
<thead><tr><th>environmental cue</th><th>possible switching mechanism</th>
<th>advantages and limitations</th><th>decision</th></tr></thead><tbody>
<tr><td>Acidic pH</td>
<td>Design a mask whose affinity for the nanobody decreases at pH
{PH['site_ph'][0]:g}&ndash;{PH['site_ph'][1]:g}. Histidine substitutions and other ionisable
residues could be used to engineer a pH-dependent interface.</td>
<td>A fully protein-based, reversible switch. However, the pH difference from blood
({PH['blood_ph'][0]}&ndash;{PH['blood_ph'][1]}) is relatively small. For a histidine with pKa
{PH['pka']:.1f}, protonation increases from approximately {PH['blood_worst']:.0%} in blood to
{PH['site_worst']:.0%} at pH {PH['site_ph'][1]:g}; achieving a large change in binding may require
several coordinated interactions.</td>
<td><b>Alternative</b></td></tr>
<tr><td>Protease activity</td>
<td>Attach the mask through a cleavable peptide linker. Local host proteases or snake venom
metalloproteinases (SVMPs) release the mask, allowing the nanobody to bind its toxin.</td>
<td>Cleavage provides an irreversible activation mechanism using a genetically encoded linker.
Sequence specificity and structural accessibility can be designed separately. The challenge is
avoiding cleavage in circulation while ensuring activation by the relevant venom and host
proteases.</td>
<td><b>Selected</b></td></tr>
<tr><td>Hypoxia</td>
<td>Attach the mask through an oxygen-sensitive chemical trigger, such as a
<b>4-nitrobenzyl-based self-immolative linker</b>. Enzymatic reduction of the nitro group under
suitable low-oxygen conditions can initiate linker breakdown and release the mask.</td>
<td>Could provide irreversible activation without proteolysis. However, the chemistry requires
reduction by enzymes and cofactors that are often intracellular. Reliable activation around an
extracellular nanobody at a bite site is therefore uncertain and would require chemical
conjugation and experimental validation.</td>
<td><b>Alternative</b></td></tr>
</tbody></table></div>

<div class="stack"><p><b>Why proteolysis?</b> All three approaches are plausible, but proteolysis
offers a relatively direct route to a modular, extracellular switch. It also builds on the
principle of intramolecular regulation I investigated during my PhD: a tethered segment controls
access to a functional protein surface.</p></div>


<div class="chart wide">{R.states_schematic()}</div>
<p class="dim" style="font-size:.9rem">The nanobody&ndash;toxin complex this work builds on is
<a href="https://www.rcsb.org/structure/{PDB_ID}">PDB {PDB_ID}</a>: the broadly neutralising
nanobody VHH5 a-sNTx bound to short neurotoxin 1 (UniProt P01426) from the red spitting cobra
<i>Naja pallida</i>. Deposited with: Ahmadi S.,
Burlet N.J., Benard-Valle M. <i>et al.</i> Nanobody-based recombinant antivenom for cobra, mamba
and rinkhals bites. <i>Nature</i> 647:716&ndash;725 (2025).
<a href="https://doi.org/10.1038/s41586-025-09661-0">10.1038/s41586-025-09661-0</a></p>
</section>


<section class="sec">
<div class="stack"><div class="eyebrow">03</div><h2>Proposed computational pipeline</h2>
</div>
<div class="chart wide">{R.architecture_schematic()}</div>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">04</div><h2>Worked computational example</h2>
<p>The pipeline was implemented and applied to {meta['title']} ({PDB_ID}), a {N['resolution']} Å X-ray structure deposited in {meta['deposition_date'][:4]}.</p>
<p>It was chosen because it is a crystallographically characterised nanobody&ndash;toxin
complex with a single binding site, which is what the geometry and the masking measurements
need. The architecture and the methods are what transfer; a viper application would start
from a viper toxin complex and redesign the mask for it.</p>
<p>What the run generated, evaluated and selected. Each of the
{N['n_evaluated']} candidates carries all three scores, computed from its own predicted
structure; the free nanobody is folded as the reference occlusion is measured against and is
not itself a candidate, which is why {N['n_structures']} structures give {N['n_evaluated']}
candidates.</p></div>
<figure class="fig wide"><div class="figbody">{R.pipeline_funnel(FUNNEL)}</div></figure>

<h3 class="mod"><span>Module 01</span>Structural and affinity constraints</h3>
<div class="stack">
<p><b>Reference structure.</b> I used the experimental VHH–toxin complex (PDB {PDB_ID}) to
locate the paratope and the possible linker attachment sites; which antigen residues face it, and
what the masks are grafted from, is Module 02. The paratope contains {N['n_paratope']} residues,
and {termini_free}.</p></div>

<div class="two">
<div class="stackfig">
{chart("affinity_by_length",
   "<b>Predicted mask-affinity range for different linker lengths.</b> The green region shows "
   "combinations of linker length and mask affinity that allow at least "
   f"{100 - model['leak_max'] * 100:.0f}% masking before cleavage and "
   f"{model['recovery_min'] * 100:.0f}% antigen binding afterwards. The upper curve represents "
   "the weakest mask affinity that still provides sufficient masking while tethered; the lower "
   "line represents the strongest affinity that permits antigen binding after cleavage. The "
   f"widest range occurs at {N['linker_best']} residues. Triangles indicate linker lengths used "
   "in the computational designs.")}
</div>
<div class="stackfig">
{chart("usable_affinity",
   f"The same constraint at the widest point, {N['linker_best']} residues."
   f"<br>N-terminus distance to the paratope {N['span_n']:.1f} Å."
   f"<br>The fully extended {N['ex_n']}-residue linker length is {N['ex_lc']:.1f} Å, but its "
   f"characteristic end-to-end distance is {N['ex_ree']:.1f} Å because the chain bends."
   f"<br>Estimated local concentration of the mask is {N['ex_ceff']:.1f} mM.")}
</div>
</div>

<div class="stack">
<p><b>Estimating linker length and mask affinity.</b> I used established models from polymer
physics and binding equilibrium to estimate which linker lengths and mask affinities could allow
both masking and release. The worm-like-chain treatment of a flexible linker follows
<a href="https://doi.org/10.1021/jp011355n">Zhou (2001)</a>; the use of an effective concentration
to describe what a disordered linker buys follows
<a href="https://doi.org/10.1073/pnas.1904813116">Sørensen and Kjaergaard (2019)</a> and
<a href="https://doi.org/10.1016/bs.mie.2020.09.012">Kjaergaard and colleagues (2021)</a>.</p>

<p><b>1 &middot; Effective concentration of the tethered mask</b></p>

<p><b>Choosing the attachment point.</b> The experimental structure shows that the N-terminus is
closer to the centre of the paratope than the C-terminus ({N['span_n']:.1f} Å versus
{N['span_c']:.1f} Å). I therefore selected N-terminal attachment.</p>

<p><b>Step 1 &mdash; Estimating how far the linker can reach.</b> The Gly/Ser linker is flexible
and can adopt many conformations rather than remaining fully extended. I used a
<i>worm-like-chain</i> model to calculate its characteristic end-to-end distance at different
linker lengths.</p>
<div class="eq">{EQ_WLC}</div>
<ul class="where">
<li><span class="sym">L<sub>c</sub></span> is the fully extended linker length, calculated as
{model['rise_per_residue_A']:.1f} Å × number of residues.</li>
<li><span class="sym">l<sub>p</sub></span> is the persistence length, set to
{model['persistence_length_A']:.1f} Å, which describes how readily the chain bends.</li>
<li><span class="sym">√⟨r²⟩</span> is the root-mean-square end-to-end distance, describing the
linker's characteristic size.</li>
</ul>
<p>For a {N['ex_n']}-residue linker, the fully extended length is {N['ex_lc']:.1f} Å, but its
characteristic end-to-end distance is only {N['ex_ree']:.1f} Å because the chain bends.</p>

<p><b>Step 2 &mdash; Estimating the local concentration of the mask.</b> Knowing the linker's
characteristic size, I used a Gaussian approximation to estimate the probability density of its
free end at the paratope, {N['span_n']:.1f} Å from the attachment point. This was converted into
an effective local concentration.</p>
<div class="eq">{EQ_CEFF}</div>
<ul class="where">
<li><span class="sym">d</span> is the distance from the attachment point to the paratope centre
({N['span_n']:.1f} Å).</li>
<li><span class="sym">⟨r²⟩</span> comes from the worm-like-chain calculation above.</li>
<li><span class="sym">10²⁷/N<sub>A</sub></span> converts probability density from Å⁻³ into molar
concentration, using Avogadro's constant. It evaluates to 1660.54.
<ul>
<li><b>Probability density:</b> How likely is the linker's free end to occupy a small volume?
Measured in Å⁻³.</li>
<li><b>Molar concentration:</b> How many molecules would occupy that volume at an equivalent
concentration? Measured in mol/L (M).</li>
</ul></li>
<li><span class="sym">C<sub>eff</sub></span> is the estimated effective local concentration of
the tethered mask, in mol/L.</li>
</ul>
<p>For the {N['ex_n']}-residue linker, this gives C<sub>eff</sub> = {N['ex_ceff']:.1f} mM. It
describes how tethering increases the mask's local availability near the paratope compared with a
freely diffusing peptide.</p>
<p>The resulting effective concentration is then used to estimate how strongly the mask must bind
to occupy the paratope before cleavage.</p>

<p><b>2 &middot; Masking before cleavage, the upper curve.</b> The fraction of construct whose
paratope is occupied by its own tethered mask is</p>
<div class="eq">P<sub>masked</sub> = C<sub>eff</sub> / (C<sub>eff</sub> + K<sub>D,mask</sub>)
<br><br>P<sub>masked</sub> ≥ {100 - model['leak_max'] * 100:.0f}%   requires   K<sub>D,mask</sub> ≤ C<sub>eff</sub> / {(1 - model['leak_max']) / model['leak_max']:.0f}</div>
<p>Requiring at least {100 - model['leak_max'] * 100:.0f}% masking therefore caps the mask
K<sub>D</sub> at a nineteenth of C<sub>eff</sub>, which is the upper curved boundary in the figure.
A short linker cannot reach the paratope, and a long one spreads the mask through a larger volume,
so C<sub>eff</sub> — and with it the loosest mask that still works — rises to a maximum and falls
again.</p>

<p><b>3 &middot; Antigen binding after cleavage, the lower line.</b> Once the linker is cut, the
released mask and the antigen compete for the same site:</p>
<div class="eq">P<sub>antigen</sub> = (A/K<sub>D,A</sub>) / (1 + A/K<sub>D,A</sub> + M/K<sub>D,mask</sub>)</div>
<p>Requiring at least {model['recovery_min'] * 100:.0f}% antigen occupancy sets a minimum
acceptable mask K<sub>D</sub>. That competition no longer involves the linker that used to hold the
mask, so the bound is the same at every length: it is the horizontal lower boundary.</p>

<p><b>Result.</b> The green region between the two boundaries is the set of linker lengths and mask
affinities predicted to satisfy both requirements. A {N['linker_best']}-residue N-terminal linker
gives the widest range, <b>{N['km_lo']:.0f} nM to {N['km_hi']:.0f} µM</b>. That is the length
whose characteristic end-to-end distance almost matches the distance the linker has to span, which
is where C<sub>eff</sub> peaks. A relatively weak mask could therefore be sufficient, provided the
linker holds it close to the paratope.</p>
<p>These estimates use a V<sub>H</sub>H–antigen K<sub>D</sub> of {N['vhh_toxin_kd_nM']:.1f} nM,
{model['toxin_conc_M'] * 1e9:.0f} nM antigen and {model['mask_conc_M'] * 1e6:.0f} µM released mask.
They are initial design constraints, not predictions of binding kinetics or of experimental
switching behaviour.</p>

</div>


<h3 class="mod"><span>Module 02</span>Generation of candidate constructs</h3>

<h4 class="subh"><span>Step 1</span>Design the protease-cleavable linkers</h4>
<figure class="fig wide"><div class="figbody">{R.linker_schematic(LINKER_GEN, LENGTHS, SUB_LEN, REF_LEN)}</div>
</figure>
<div class="stack">
<p>The {R.word(LINKER_GEN['n_controls'])} controls are a Gly/Ser linker with no designed cleavage
motif and a scrambled version of the highest-scoring sequence.</p>
<p><b>What this does not show.</b> CleaveNet was trained on human MMPs, not on snake
venom metalloproteinases. These results demonstrate how to design a desired human-MMP cleavage
profile; they do not establish selectivity at a viper bite site.</p></div>

<h4 class="subh"><span>Step 2</span>Identify the interface and design the masks</h4>
<div class="two">
<div class="viewer">
<div class="vbar"><span class="eyebrow">Deposited structure &middot; {PDB_ID}</span>
<button data-spin="solved">spin</button></div>
<div class="vstage" id="stage-solved"></div>
<div class="vkey"><span><i style="background:#337D79"></i>framework</span>{R.CDR_KEY}
<span><i style="background:#C99555"></i>short neurotoxin 1</span>
<span><i style="background:#C06014"></i>buried on binding, sticks</span>
<span><i style="background:#C06014"></i>within {N['mask_cutoff']:g} &#8491;, sticks</span>
<span><i style="background:#B8541F"></i>antigen face within {N['mask_cutoff']:g} &#8491;</span>
<span><b style="color:#12191F">N</b> / <b style="color:#6E7880">C</b> termini</span></div>
</div>
<div class="stackfig">
{chart("paratope_footprint",
  f"<b>Nanobody&ndash;antigen interface in PDB {PDB_ID}.</b> Bars show the surface buried by each "
  f"nanobody residue on antigen binding. The {N['n_near']} marked residues lie within "
  f"{N['mask_cutoff']:g} Å of the antigen; the same cutoff identifies the "
  f"{N['n_epitope_near']} antigen residues the masks are grafted from. A closer "
  f"{N['contact_cutoff']:g} Å cutoff picks out the {N['n_contacts']} residues making polar "
  f"contacts, reported so this interface can be compared with others. The broader paratope of "
  f"{N['n_paratope']} residues, each losing more than 1 Å² of surface, is the footprint "
  f"masking is measured against in Modules 03 and 04.")}
</div>
</div>
<div class="stack">
<p>The {N['n_epitope_near']} contacting antigen residues are not one continuous peptide.
{SEG_SPLIT} Both regions are carried into the mask design below, separately and together, because a mask that
reproduces only the larger one may or may not be enough.</p></div>

<h5 class="subh2">Designing the masks</h5>
<div class="refmask">
<div class="eyebrow">Reference mask &middot; <span class="mono">{CORE_ID}</span></div>
<div class="seq">{CORE_SEQ}</div>
<div class="sub">{CORE_LEN} amino acids &middot; antigen residues {CORE_RANGE}</div>
<p>The two terminal cysteines are bonded to each other in the native antigen, so cutting here
keeps the staple that holds the loop in the shape the nanobody binds. That is also the argument
against simply taking more sequence: the extensions reach past this pair and pick up
{ORPHAN_TEXT}. Longer is not automatically better.</p>
</div>
<figure class="fig wide"><div class="figbody">{R.mask_families(MASKS, HOTSPOTS)}</div>
<figcaption>The {N['n_masks']} masks, grouped by the question each was built to answer and drawn
against the antigen residues they come from. Each group moves one variable: how much of the
contact region to take, whether the second region helps, whether the sequence can change while
the buried residues stay, and whether the disulphide matters. The last group is the controls the
others have to beat.</figcaption></figure>

<h5 class="subh2">Combining them into constructs</h5>
<div class="sharedbox">
<div class="stack">
<p>Every mask was paired with every linker: {N['n_masks']} masks
&times; {lib['n_linkers']} linkers = {lib['n_pairs']} pairs, and the free nanobody as a reference,
giving <b>{N['n_constructs']} constructs</b> to fold.</p>
<p><b>Structure prediction and paratope occlusion.</b> I used ESMFold to predict the
structures of the complete mask–linker–nanobody constructs. The free nanobody was
predicted separately as a reference (pLDDT {FREE_PLDDT:.1f}).</p>
<p>I then compared the solvent-accessible surface area of the {N['n_paratope']} paratope residues
in each construct with that of the free nanobody ({FREE_SASA:,.0f} Å&sup2;). The difference gives
the predicted percentage of paratope occlusion.</p>
<p>This measures surface coverage in a single predicted structure, not mask affinity or
functional inhibition.</p></div>
<figure class="fig wide"><div class="figbody">{R.construct_schematic(REF_CONSTRUCT)}</div>
</figure>
<div class="stack">
<p><b>Reading a construct name.</b> <span class="mono">{CORE_ID}_cn01</span> is the main
antigen-derived mask joined to CleaveNet linker 01, and
<span class="mono">mpnn00_{CORE_ID}_cn01</span> is a ProteinMPNN redesign of that construct. The
ids stay as the files are named; the figure above gives each mask a readable label.</p></div>
<details class="tech"><summary>Any construct in the library, folded</summary>
<div class="techbody">
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
</div></details>
</div>

<h4 class="subh"><span>Step 3</span>Redesign the mask sequences</h4>
<div class="stack">
<p>Each folded construct was passed to ProteinMPNN with the nanobody and the linker held fixed,
so only the mask positions could change. That produced <b>{N['n_redesign']} distinct redesigned
masks</b> across {N_PARENTS} parent constructs.</p>
<p>{N['n_improved']} of the {N['n_redesign']} redesigns scored better than the graft they came
from on ProteinMPNN's own measure. The {N['n_not_improved']} that did not are all redesigns of
<span class="mono">{N['not_improved_parents'][0]}</span>, the poly-glycine control, which already
holds the best native score of any parent at {N['best_native_score']:.2f}: a run of glycines fits
any backbone, so there is nothing there to improve on.</p>
<p>That measure is how well a sequence fits the backbone it was given, which is not the same
question as whether the mask covers the paratope, and the two barely move together:
r = {N['mpnn_corr']:+.2f} across all {N['n_redesign']} redesigns.</p></div>
{chart_wide("redesign_score_vs_occlusion",
   f"ProteinMPNN score against paratope occlusion for all {ALL_N} constructs carrying a mask, "
   "coloured by the family the mask belongs to and shaped by whether it is a library construct "
   f"or a redesign. Across the whole library the two are unrelated, r = {ALL_CORR:+.2f}: a "
   "sequence that fits the backbone better does not occlude the paratope more, so the design "
   "model's own score cannot stand in for a functional one. Each score is the mean over "
   f"{LIBSC['n_orders'] if LIBSC else 0} decoding orders, which carry a spread of "
   f"±{ALL_SD:.2f}.")}
<div class="stack">
<p><b>Horizontal axis &mdash; ProteinMPNN score.</b> How compatible ProteinMPNN considers the
sequence with the backbone it was given. Lower is better, so constructs further left have the
better ProteinMPNN scores.</p>
<p><b>Vertical axis &mdash; paratope occlusion.</b> Calculated from the predicted structures: how
much less solvent-accessible surface the paratope has in the fusion than in the free nanobody.
Higher values mean more of the paratope is predicted to be covered.</p></div>

</section>

<section class="sec">
<div class="stack"><div class="eyebrow">05</div><h2>Scoring and computational results</h2>
<p>Each construct is scored on three independent requirements. Two of the three work, in the sense
that they separate designs from matched controls. The third does not, and that is the result this
section is mainly about.</p></div>

<h3 class="mod"><span>Module 03</span>Structural feasibility</h3>
{mri("I used ESMFold to predict each complete mask–linker–nanobody construct. "
     "Structural feasibility was assessed on four criteria: prediction confidence, preservation "
     "of the V<sub>H</sub>H fold, exposed hydrophobic surface and unpaired cysteines. The four "
     "components were combined into a single heuristic score, with penalties for poor "
     "performance on any one of them.",
     f"All <b>{N['feas_n']} predicted structures</b> were scored, the "
     f"{N['n_redesign']} redesigned masks and the free V<sub>H</sub>H reference included. The "
     f"mean V<sub>H</sub>H core RMSD was <b>{N['rmsd_mean']:.2f} Å</b>, so the predicted "
     "nanobody fold changes little.",
     "The predicted V<sub>H</sub>H fold is largely preserved once the mask and linker are "
     "added. Exposed hydrophobic residues and unpaired cysteines remain potential expression "
     "and folding liabilities.")}
<details class="tech"><summary>How is structural feasibility calculated?</summary>
<div class="techbody">
<div class="quad">
<div><div class="k">1 &middot; Fold confidence</div><b>ESMFold pLDDT</b>
<span>Is the complete fusion confidently predicted?</span></div>
<div><div class="k">2 &middot; Fold preservation</div><b>V<sub>H</sub>H core RMSD</b>
<span>Does the nanobody keep its original structure?</span></div>
<div><div class="k">3 &middot; Surface liability</div><b>Exposed hydrophobic residues</b>
<span>How much of the exposed surface is hydrophobic?</span></div>
<div><div class="k">4 &middot; Cysteine liability</div><b>Unpaired cysteines</b>
<span>Are any cysteines left with no predicted disulphide partner?</span></div>
</div>
<div class="eq">S<sub>feasibility</sub> = S<sub>fold</sub> &times; S<sub>preservation</sub>
&times; S<sub>surface</sub> &times; S<sub>cysteine</sub></div>
<p>All four components run from 0 to 1, and a higher score means fewer predicted liabilities.
What each one computes:</p>
<table class="tbl">
<colgroup><col style="width:30%"><col style="width:70%"></colgroup>
<thead><tr><th>component</th><th>calculation</th></tr></thead><tbody>
<tr><td>fold confidence</td><td>mean pLDDT of the whole construct divided by
{FEAS['plddt_ref']:g}, capped at 1</td></tr>
<tr><td>V<sub>H</sub>H preservation</td><td>1 &minus; RMSD / {FEAS['rmsd_ref']:g}, clipped to
0&ndash;1</td></tr>
<tr><td>surface liability</td><td>1 &minus; f<sub>hydrophobic</sub> /
{FEAS['phobic_ref']:g}, clipped to 0&ndash;1</td></tr>
<tr><td>cysteine liability</td><td>1 &minus; {FEAS['cys_penalty']:g} n<sub>unpaired</sub>,
with a floor of {FEAS['cys_floor']:g}</td></tr>
</tbody></table>
<p>For the cysteine term the code counts cysteine sulfur atoms with no other cysteine sulfur
within {FEAS['ss_max']:g} Å in the predicted structure. That is an approximate structural check,
not experimental confirmation of disulphide bonding.</p>
</div></details>
{chart_wide("feasibility_components", "The four components by mask family, in the same colours "
       "the families carry elsewhere on this page. A low score on any one component is enough to "
       "make a construct unusable, since the four are multiplied. The dashed line is the free "
       f"nanobody, which already scores {N['ref_surface']:.2f} on exposed hydrophobic residues, "
       "so that component is read as a change against it rather than on its own.")}
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

<h3 class="mod"><span>Module 04</span>Paratope occlusion</h3>
<div class="stack">
<p>I calculated how much of the nanobody's antigen-binding surface becomes buried when the mask
and linker are attached. For each ESMFold-predicted construct, I compared the solvent-accessible
area of the {N['n_paratope']} paratope residues with that of the free nanobody.</p>
<p>The designed masks showed no greater occlusion than their controls, indicating that surface
coverage alone cannot establish specific mask binding.</p></div>
<div class="stack">
<p><b>Controls for paratope occlusion.</b> To test whether surface coverage depends on the
antigen-derived sequence, I compared the main contact-loop mask
(<span class="mono">{CORE_ID}</span>) with scrambled, truncated and poly-glycine controls, each
built from it. All four were paired with the same {lib['n_linkers']} linker variants, so each
mask below stands for {lib['n_linkers']} constructs, and the mask sequence is the only thing
that differs between them.</p></div>
<table class="tbl">
<colgroup><col style="width:21%"><col style="width:36%"><col style="width:33%"><col style="width:10%"></colgroup>
<thead><tr><th>mask</th><th>what it is</th><th>the question it asks</th>
<th>mean</th></tr></thead><tbody>{OCC_CONTROL_ROWS}</tbody></table>
<div class="stack">
<p><b>Effect of linker length on paratope occlusion.</b> Predicted coverage across linker
lengths of {N['len_lo']}&ndash;{N['len_hi']} residues, with the mask and the cleavage sequence
held constant along each line, so the only thing varying is the length.</p></div>
<div class="two">
<div class="stackfig">
{chart("occlusion_by_group", "Every construct as a point, with the mask's mean marked. "
   "The antigen-derived mask does not outperform its controls.")}
<details class="tech"><summary>Every antigen-derived mask</summary><div class="techbody">
<p>The other antigen-derived masks, each against the same
{lib['n_linkers']} linkers. Extensions, substitutions and the cysteine-free variant all sit in
the same band.</p>
{chart("occlusion_by_mask", "Mean occlusion per antigen-derived mask.")}
</div></details>
</div>
<div class="stackfig">
{chart("occlusion_by_linker", "One antigen-derived mask against its two controls, all on "
   "the same cleavage motif. Longer linkers generally raise predicted occlusion, but the "
   "scrambled and poly-glycine controls rise with them, so part of the effect is geometric: a "
   "flexible peptide covers some of the surface without a specific interaction.")}
<details class="tech"><summary>Every mask with a length series</summary><div class="techbody">
{chart("occlusion_by_linker_all", "The same comparison for every mask carrying a length "
   "series, coloured by family.")}
</div></details>
</div>
</div>
{R.ensemble_viewer("ens-mask", "masking",
   "Every folded construct, coloured by paratope occlusion. The masks do not concentrate over "
   "the antigen-binding residues, which is the result reported above.",
   R.ramp_legend("least occluding", "most"))}

<h3 class="mod"><span>Module 05</span>Linker cleavability</h3>
<div class="stack">
<p><b>A linker must contain a suitable protease cleavage sequence and be accessible in the
folded construct.</b></p></div>
{mri("I combined two measures. <b>Sequence selectivity:</b> the CleaveNet score calculated "
     f"during linker design (Module 02), measuring predicted preference for {CLEAVE_TARGETS} "
     f"over {R.word(LINKER_GEN['n_off'])} off-target MMPs. <b>Structural accessibility:</b> "
     "calculated from the ESMFold structure, based on cleavage-site extension, solvent exposure "
     "and space for a protease to approach. The cleavability score is the product of the "
     "non-negative selectivity score and the three accessibility measures.",
     "The designed sequences showed greater predicted selectivity than the controls, but "
     "accessibility varied substantially between constructs. Sequence selectivity and "
     "accessibility were moderately negatively correlated "
     f"(r = {N['sel_ctx_r']:+.2f}).",
     "A suitable cleavage sequence is not necessarily accessible to a protease. The combined "
     "score accounts for both requirements, but does not predict a cleavage rate. CleaveNet is "
     "trained on human MMPs, so cleavage by snake venom metalloproteinases remains untested.")}
<h5 class="subh2">Can a protease reach the cleavage site?</h5>
<figure class="fig wide"><div class="figbody">{R.accessibility_checks(
   PROBE_R, PROBE_N, PROBE_OFF, CLEFT_HALF, EXT_RISE)}</div></figure>
<details class="tech"><summary>How the clearance probe is placed</summary><div class="techbody">
<p>The sphere stands in for an MMP catalytic domain, which is about 160 residues and so roughly
{PROBE_R * 2:.0f} Å across: it is given a radius of {PROBE_R:g} Å. Its centre is put
{PROBE_R + PROBE_OFF:g} Å from the scissile C&#945;, so its surface stands {PROBE_OFF:g} Å clear
of the bond rather than overlapping it, and the {PROBE_N} directions are spread over the sphere
on a Fibonacci lattice. A direction counts as open when no heavy atom of the construct lies
inside the probe; the cleavage segment's own atoms are excluded, since the enzyme is binding
them rather than being blocked by them.</p></div></details>
{R.annotated_equation([
   ("S<sub>cleavability</sub>", "Combined cleavability score", "used in Module 06"), "=",
   ("max(0, S<sub>selectivity</sub>)", "Sequence selectivity", "from Module 02"), "&times;",
   ("S<sub>extension</sub> &middot; S<sub>exposure</sub> &middot; S<sub>clearance</sub>",
    "Structural accessibility", "calculated in Module 05")])}
<div class="stack">
<p>The three accessibility measures are multiplied, so a cleavage site must be sufficiently
extended, exposed and accessible to a protease-sized object.</p></div>
<div class="two top">
<div class="stackfig">
{chart("cleavability_panel", "Sequence selectivity versus structural accessibility. CleaveNet "
   "selectivity was calculated during linker design (selectivity = mean target score &minus; "
   f"mean off-target score). Higher is a stronger predicted preference for {CLEAVE_TARGETS} "
   f"over the {R.word(LINKER_GEN['n_off'])} off-target MMPs.<br>Accessibility was assessed "
   "separately in each predicted fusion structure.<br>One number per construct, measured on "
   "its own ESMFold model, coloured by mask family. The product of the three checks above: "
   "extension, exposure and clearance, each from 0 to 1.")}
</div>
<div class="stackfig">
{chart("cleavability_components", "Cleavage-site accessibility across "
   f"{N['rank_total']} constructs, ordered from lowest to highest overall score. Extension, "
   "exposure and clearance are multiplied to give overall accessibility (grey). Limited "
   "protease clearance is the main constraint.")}
</div>
</div>
{R.ensemble_viewer("ens-cleave", "cleavability",
   "Every folded construct aligned on the nanobody, with each mask and linker coloured by its "
   "predicted cleavability: sequence recognition multiplied by structural accessibility.",
   R.ramp_legend("least cleavable", "most"))}

<h3 class="mod"><span>Module 06</span>Candidate selection</h3>
{mri(f"I ranked {N['rank_total']} constructs using three criteria: structural feasibility, "
     "paratope occlusion and predicted linker cleavability. Pareto optimisation identified "
     "designs that were not outperformed across all three criteria. I then ranked these by "
     "their weakest normalised score and selected different mask sequences to retain "
     "diversity.",
     f"Of {N['rank_total']} constructs, <b>{N['rank_front']}</b> were nondominated. "
     f"<b>{R.word(N['rank_shortlist']).capitalize()}</b> were shortlisted for further "
     "structural analysis.",
     "This provides a systematic way to prioritise designs without relying on a weighted "
     "average of the three scores. The ranking remains preliminary, particularly because "
     "surface occlusion does not reliably distinguish designed masks from controls.")}
<figure class="fig wide"><div class="figbody">{R.selection_funnel(FUNNEL)}</div></figure>
{chart_wide("pareto_panels", "The three scores plotted pairwise, with the nondominated set "
   "marked.")}
<details class="tech"><summary>How candidates were selected</summary><div class="techbody">
<p>Each construct was assessed for structural feasibility, paratope occlusion and linker
cleavability. These are three separate objectives; occlusion measures surface coverage, not
mask-binding affinity.</p>
<p><b>1 &middot; Basic requirements.</b> Before comparing designs, I checked V<sub>H</sub>H
confidence, preservation of the nanobody fold, predicted protease selectivity and cleavage-site
accessibility. These checks prevent candidates with clearly unfavourable properties from being
promoted by other scores.</p>
<p><b>2 &middot; Pareto selection.</b> A design was removed if another performed at least as
well on all three objectives and better on at least one. This left {N['rank_front']}
nondominated candidates, representing different trade-offs between the three scores.</p>
<p><b>3 &middot; Shortlisting.</b> I normalised each score to the range observed across the
library and ranked the {N['rank_front']} candidates by their lowest normalised score (maximin).
This favours balanced designs over those that perform well in two criteria but poorly in the
third. I then retained one construct per distinct mask sequence, producing a shortlist of
{R.word(N['rank_shortlist'])} candidates for further structural assessment.</p>
<div class="eq">R(x) = min<sub>i</sub> Ŝ<sub>i</sub>(x)</div>
<p>Here Ŝ<sub>i</sub>(x) is each normalised objective. The lowest of the three determines the
ranking score.</p></div></details>
{R.normalised_table(RANK_SHORTLIST_ROWS, RANK_OBJECTIVES)}
<div class="stack">
<p class="dim">Scores are min&ndash;max normalised across the {N['rank_total']} candidates.
These are relative ranking values, not experimental probabilities.</p></div>
<details class="tech"><summary>The rest of the nondominated set</summary><div class="techbody">
<p>The other {N['rank_front'] - N['rank_shortlist']} candidates on the front, in the same order.
Each was dropped from the shortlist because a higher-ranked candidate already carried its
mask.</p>
{R.normalised_table(RANK_FRONT_REST, RANK_OBJECTIVES)}</div></details>
<div class="note warn"><b>A poly-glycine control survives the selection.</b>
<span class="mono">polyG_cn01</span> is on the shortlist. That is why the shortlist needs the
Boltz-2 assessment that follows, and why the occlusion score should eventually be replaced or
improved.</div>

<h3 class="mod"><span>Final screen</span>Structural assessment with Boltz-2</h3>
<div class="stack">
<p>This is not a fourth score. The three preliminary scores selected {N['rank_shortlist']} constructs out of {N['rank_total']}; this is a closer look at those {R.word(N["rank_shortlist"])}.</p></div>
{mri("I first tested whether Boltz-2 could reproduce the experimental V<sub>H</sub>H&ndash;antigen interface. With evolutionary alignments, " f"{N['boltz_ctrl']:.0f}% of antigen-contacting nanobody residues lay in the known "
     f"paratope, against {N['boltz_ctrl_noalign']:.0f}% without. I then predicted the "
     f"{R.word(N['rank_shortlist'])} shortlisted constructs with their antigen and "
     "examined whether the mask occupied the paratope while the antigen was displaced.",
     f"One construct, <span class='mono'>{FP_ID}</span>, met both structural criteria: "
     f"<b>{FP_MASK_ON:.0f}%</b> of its mask-contacting nanobody residues were on the "
     f"paratope, while antigen contacts there fell from {FP_REF_ON:.0f}% in the free "
     f"nanobody reference to <b>{FP_ANTIGEN_ON:.0f}%</b>.",
     "That construct was prioritised for experimental testing. The prediction suggests "
     "interference with antigen binding; it does not demonstrate mask affinity or "
     "release after cleavage.")}
<div class="two">
<div class="viewer">
<div class="vbar"><span class="eyebrow">Reference &middot; free nanobody with antigen</span>
<button data-spin="finalref">spin</button></div>
<div class="pairstat"><b>{FP_REF_ON:.0f}%</b><span>antigen on paratope</span></div>
<div class="vstage" id="stage-final-ref"></div>
<div class="vkey"><span><i style="background:#337D79"></i>nanobody</span>{R.CDR_KEY}
<span><i style="background:#C06014"></i>paratope, sticks</span>
<span><i style="background:#E0C49A"></i>antigen</span></div>
</div>
<div class="viewer">
<div class="vbar"><span class="eyebrow">Selected &middot; <span class="mono">{FP_ID}</span> with
antigen</span><button data-spin="finalmask">spin</button></div>
<div class="pairstat"><b>{FP_ANTIGEN_ON:.0f}%</b><span>antigen on paratope</span>
<b class="ok">{FP_MASK_ON:.0f}%</b><span>mask on paratope</span></div>
<div class="vstage" id="stage-final-mask"></div>
<div class="vkey"><span><i style="background:#337D79"></i>nanobody</span>{R.CDR_KEY}
<span><i style="background:#C06014"></i>paratope, sticks</span>
<span><i style="background:#8C72AD"></i>mask and linker</span>
<span><i style="background:#E0C49A"></i>antigen</span></div>
</div>
</div>
<div class="stack">
<p class="dim">Both predictions are shown in the same frame: the construct&rsquo;s model is
superposed on the control&rsquo;s by the nanobody alpha carbons alone, which agree to
{FP_RMSD:.2f} Å. Anything that moves between the two views moves relative to the nanobody.</p>
</div>
<table class="tbl">
<colgroup><col style="width:26%"><col style="width:17%"><col style="width:16%"><col style="width:41%"></colgroup>
<thead><tr><th>construct</th><th>antigen on paratope</th><th>mask on paratope</th>
<th>what it means</th></tr></thead><tbody>{VERIFY_ROWS}</tbody></table>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">06</div><h2>Experimental collaboration and further development</h2></div>

<h4 class="subh"><span>First order</span>Three constructs, for diagnosis</h4>
<div class="stack">
<p>The first order is three constructs, and its purpose is diagnosis. Each
measurement fails against a different part of the pipeline, so three is enough to start.</p></div>
<table class="tbl">
<colgroup><col style="width:38%"><col style="width:22%"><col style="width:40%"></colgroup>
<thead><tr><th>measurement</th><th>reads out</th><th>if it fails, the problem is</th></tr></thead>
<tbody>
<tr><td>Purification &mdash; &Auml;KTA trace, single band on the gel</td>
<td>can it be made at all</td>
<td>the feasibility score: fold, aggregation, developability</td></tr>
<tr><td>Biolayer interferometry, antigen affinity with the protease absent</td>
<td>is it masked</td><td>the masking half of the unmasking score</td></tr>
<tr><td>Biolayer interferometry, antigen affinity with the protease present, with a band shift
on the gel when cleaved</td>
<td>does it activate</td><td>the trigger score, or the release half</td></tr>
</tbody></table>

<h4 class="subh"><span>Then twenty to thirty</span>The pipeline itself gets scored</h4>
<figure class="fig wide"><div class="figbody">{R.validation_matrix()}</div></figure>
<div class="note warn"><b>The batch has to include designs the pipeline rejected.</b> If only
top-ranked candidates are made, every measurement lands in the left-hand column: true and false
positives can be counted, so precision is computable, but false negatives cannot be, so recall is
not. Spending perhaps a quarter of the batch on designs the score rejected is what turns a hit
list into a validation set.</div>

<h4 class="subh"><span>Beyond 100 constructs</span>Learning from experimental results</h4>
<div class="stack">
<p>With more than 100 experimentally characterised constructs, the measured expression, masking
and cleavage responses could be combined with the sequence and structural features already
calculated by the pipeline. This would provide a dataset for training and testing predictive
models, potentially replacing the current heuristic scores with models informed by experimental
results.</p></div>
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
Chemistry B</i> 105:6763–6766 (2001).
<a href="https://doi.org/10.1021/jp011355n">10.1021/jp011355n</a></li>
<li>Sørensen C.S., Kjaergaard M. Effective concentrations enforced by intrinsically disordered
linkers are governed by polymer physics. <i>PNAS</i> 116:23124–23131 (2019).
<a href="https://doi.org/10.1073/pnas.1904813116">10.1073/pnas.1904813116</a></li>
<li>Kjaergaard M., Glavina J., Chemes L.B. Predicting the effect of disordered linkers on effective
concentrations and avidity. <i>Methods in Enzymology</i> 647:145–171 (2021).
<a href="https://doi.org/10.1016/bs.mie.2020.09.012">10.1016/bs.mie.2020.09.012</a></li>
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
<script>{R.VIEWER_JS.replace("__MODELS__", json.dumps(models)).replace("__PARATOPE__", json.dumps(para_resi)).replace("__SOLVED__", json.dumps(f"structures/{PDB_ID}_{REFERENCE_VHH_CHAIN}{iface['toxin_chain']}.pdb")).replace("__ENSEMBLE__", json.dumps(ENSEMBLE_IDS)).replace("__TOPPARA__", json.dumps(top_paratope)).replace("__TOXCONTACT__", json.dumps(toxin_contact)).replace("__CDRS__", json.dumps(cdr_resi)).replace("__SCORES__", json.dumps(SCORES)).replace("__FINALPAIR__", json.dumps(FINAL_PAIR_JS))}</script>
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

# Each module becomes a <details>, so a section can be put away once read.
HTML = R.fold_modules(HTML)

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
