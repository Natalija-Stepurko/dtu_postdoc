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
# names this block defines that the live page reads, defined first so that a
# stage which has not run leaves a missing value rather than a NameError
corr = float("nan")
rrows, rocc, grafts, improved = [], [], [], 0
rbest = None

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
N["ree_best"] = bN["mean_end_to_end_A"]
N["linker_best_c"] = bC["linker_residues"]
N["window_c"] = bC["window_decades"]
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
_lk_len = next((len(x["sequence"]) for x in lib["linkers"]), 0)
N["linker_len"] = _lk_len
_scan_by_n = {r["linker_residues"]: r["window_decades"]
              for r in tether["scan"] if r["terminus"] == "N"}
N["window_cost"] = (_scan_by_n.get(tether["best_per_terminus"]["N"]["linker_residues"], 0)
                    - _scan_by_n.get(_lk_len, 0))
N["n_epitope_contacts"] = iface.get("n_epitope_contacts", 0)
N["n_epitope"] = iface.get("n_epitope", 0)
try:
    _seg = load_json("03a_masks.json")["segments"][0]
    N["mask_face"] = f"{_seg['start_label']}\u2013{_seg['end_label']}"
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
_acc_by_id = {r["id"]: r["context_score"] for r in (trigger["results"] if trigger else [])}
for _key, _d in (("feasibility", _feas_by), ("masking", _occ_by), ("cleavability", _acc_by_id)):
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
    LINKERS = _lk["linkers"]
except FileNotFoundError:
    N.update({"n_linkers": 0, "n_generated": 0, "sel_min": 0.0, "sel_max": 0.0,
              "sel_ctrl": 0.0, "target_proteases": "", "off_target_proteases": ""})
    LINKERS = []


COUNTS = _count_rows()


#: Which module each stage belongs to, and the method that produced its
#: number. Keeping it here means the table cannot claim a method a stage does
#: not use.
#:
#: Stage numbers are execution order of the code and collide with the module
#: numbers while meaning something different, so they are not shown as a
#: column. The script name is given instead, which is unambiguously a file.
STAGE_MODULE = {'01': ('01', 'Shrake–Rupley SASA, biotite 1.6'), '02': ('01', 'worm-like chain; competitive equilibria'), '02b': ('02b', 'CleaveNet conditional generator + predictor'), '18': ('02b', 'position weight matrix and nearest neighbour over the Kukreja library'), '16': ('02b', 'flank-length sweep, substrate held constant; ESMFold'), '03a': ('02a', 'antigen contact face at 3.5 Å; window variants'), '03': ('02a', 'the design matrix: one variable at a time'), '04': ('03', 'ESMFold facebook/esmfold_v1'), '07': ('02a', 'ProteinMPNN v_48_020, nanobody and linker fixed'), '12': ('03', 'pLDDT; Cα core RMSD; SASA; SG–SG distance'), '10': ('05', 'extension × exposure × 17 Å clearance'), '17': ('04', 'Boltz-2.2.1, mask + nanobody, no tether'), '09': ('04', 'Boltz-2.2.1 with every chain unaligned: a positive control, on a complex the crystal structure already answers'), '19': ('04', 'Boltz-2.2.1, whole construct with antigen, alignments on public chains'), '13': ('06', 'Pareto non-domination'), '15': ('03', 'the same measures on independent complexes')}


SCRIPTS = {path.name.split("_")[0]: path.name
           for path in sorted((Path(__file__).resolve().parent).glob("*.py"))}


def _exp(stage, name, question, loader):
    """One row. A stage with no output is listed, not omitted."""
    module, method = STAGE_MODULE.get(stage, ("\u2014", ""))
    base = {"stage": stage, "name": name, "question": question,
            "module": module, "method": method,
            "script": SCRIPTS.get(stage, "")}
    try:
        return {**base, "status": "run", "result": loader()}
    except (FileNotFoundError, KeyError, IndexError, TypeError, ValueError,
            ZeroDivisionError, StopIteration):
        return {**base, "status": "not run", "result": ""}


def _f(name):
    return load_json(name)


def _build_register():
    R_ = []
    a = R_.append

    a(_exp("01", "Reference interface", "Which residues does the antigen bury, "
           "and how far is each terminus from them?",
           lambda: f"{_f('01_reference.json')['chain_pairs'][REFERENCE_VHH_CHAIN]['n_paratope']} "
                   f"paratope residues; N-terminus "
                   f"{_f('01_reference.json')['terminal_geometry']['n_term_to_paratope_centroid_A']:.1f} Å "
                   f"from the paratope centroid against "
                   f"{_f('01_reference.json')['terminal_geometry']['c_term_to_paratope_centroid_A']:.1f} Å "
                   f"for the C-terminus"))

    a(_exp("02", "Tether model", "What mask affinity can masking and release both tolerate?",
           lambda: f"{_f('02_tether.json')['best_per_terminus']['N']['window_decades']:.2f} decades "
                   f"at {_f('02_tether.json')['best_per_terminus']['N']['linker_residues']} residues; "
                   f"the window closes above "
                   f"{_f('02_tether.json')['affinity']['window_closes_above_M']*1e9:.0f} nM "
                   f"antigen affinity"))

    a(_exp("02b", "Linker generation", "Can a linker be cut at the site and resist elsewhere?",
           lambda: f"{_f('02b_linkers.json')['n_generated']} substrates generated, "
                   f"{_f('02b_linkers.json')['n_kept']} kept at selectivity "
                   f"{max(x['selectivity'] for x in _f('02b_linkers.json')['linkers'] if x.get('selectivity') is not None):+.2f} "
                   f"against {next(x['selectivity'] for x in _f('02b_linkers.json')['linkers'] if x['id']=='ctrl_gs'):+.2f} "
                   f"for a Gly/Ser control"))

    a(_exp("18", "Cleavability against measurement",
           "Do the substrates hold up against the measured library, not the model?",
           lambda: f"position weight matrix {_f('18_cleavage_evidence.json')['calibration']['pwm_vs_measured_r']:+.2f} "
                   f"out of sample; agrees with CleaveNet at r = "
                   f"{_f('18_cleavage_evidence.json')['cleavenet_vs_pwm_r']:+.2f} across the linkers"))

    a(_exp("16", "Linker length sweep",
           "Does the mask reach the paratope, and at what length?",
           lambda: (lambda d: f"at {min(d['lengths'])} residues the mask is 26 Å away and makes no "
                              f"contact; contact appears at 24 over "
                              f"{d['n_constructs']} constructs with the substrate held constant")
                   (_f('16_length_sweep.json'))))

    a(_exp("03a", "Mask design", "What does the antigen present at the paratope?",
           lambda: (lambda d: f"contact face {d['segments'][0]['start_label']}–{d['segments'][0]['end_label']}, "
                              f"{d['segments'][0]['n_contact']} residues within {d['cutoff_A']} Å, "
                              f"burying {d['segments'][0]['buried_A2']:.0f} Å²; "
                              f"{len(d['masks'])} candidates, "
                              f"{len(d['primary_cysteines_orphaned'])} cysteines orphaned by the cut")
                   (_f('03a_masks.json'))))

    a(_exp("03", "Construct assembly", "How many constructs does the design space give?",
           lambda: f"{_f('03_constructs.json')['n_constructs']} constructs from "
                   f"{len(_f('03_constructs.json')['masks'])} masks and "
                   f"{_f('03_constructs.json')['n_linkers']} linkers"))

    a(_exp("04", "Fold the constructs", "Does the nanobody fold survive the fusion?",
           lambda: f"{_f('04_folding.json')['n_folded']} folded; free nanobody pLDDT "
                   f"{_f('04_folding.json')['free_vhh_plddt']:.0f}, paratope SASA "
                   f"{_f('04_folding.json')['free_vhh_paratope_sasa_A2']:.0f} Å²"))

    a(_exp("07", "ProteinMPNN redesign", "Can the mask sequence be improved on the same backbone?",
           lambda: f"{_f('07_redesign.json')['n_unique']} unique masks from "
                   f"{len(_f('07_redesign.json')['parents'])} parents"))

    a(_exp("12", "Structural feasibility", "Can the construct be expressed and folded?",
           lambda: (lambda d: f"mean score {sum(r['feasibility_score'] for r in d['results'])/len(d['results']):.2f}; "
                              f"cysteine pairing measured from the model rather than counted")
                   (_f('12_feasibility.json'))))

    a(_exp("10", "Cleavability in context", "Can a protease reach the linker in the folded construct?",
           lambda: (lambda d: f"context score {min(r['context_score'] for r in d['results']):.2f} to "
                              f"{max(r['context_score'] for r in d['results']):.2f}; independent of "
                              f"sequence recognition")
                   (_f('10_trigger.json'))))

    a(_exp("17", "Mask binding, untethered",
           "Does the mask find the paratope with no linker holding it there?",
           lambda: (lambda d: f"positive control (antigen) "
                              f"{next(r['fraction_of_contacts_on_paratope'] for r in d['results'] if r['case']=='antigen_positive_control'):.0%} "
                              f"on the paratope; "
                              f"the grafts place "
                              f"{max(r['fraction_of_contacts_on_paratope'] for r in d['results'] if 'core' in r['case'] and 'scram' not in r['case']):.0%} "
                              f"of their contacts there against "
                              f"{next(r['fraction_of_contacts_on_paratope'] for r in d['results'] if r['case']=='polyG'):.0%} "
                              f"for poly-glycine and "
                              f"{next(r['fraction_of_contacts_on_paratope'] for r in d['results'] if 'scram' in r['case']):.0%} "
                              f"for a scramble")
                   (_f('17_mask_binding.json'))))

    a(_exp("09", "The alignment control",
           "Can the method reproduce an interface that is already known?",
           lambda: f"predicted with no alignments at all, the free nanobody and "
                   f"its own antigen put only "
                   f"{next(r['fraction_of_contacts_on_paratope'] for r in _f('09_complexes.json')['results'] if r['case']=='control_free_vhh'):.0%} "
                   f"of contacts on the paratope the crystal structure shows: a "
                   f"positive control the method failed, before any mask was "
                   f"involved"))

    a(_exp("19", "Verify the front", "With the mask on, can the antigen still bind?",
           lambda: (lambda d: f"{d['n_worth_testing']} of {d['n_verified']} verified constructs "
                              f"block the antigen against a free-nanobody control at "
                              f"{d['free_vhh_fraction_on_paratope']:.0%}")
                   (_f('19_verification.json'))))

    a(_exp("13", "Multiobjective ranking", "Which constructs does nothing else beat?",
           lambda: f"{_f('13_ranking.json')['n_front']} of "
                   f"{_f('13_ranking.json')['n_total']} nondominated, ordered by weakest objective"))

    a(_exp("15", "Independent complexes", "Does this interface look like other nanobody–toxin pairs?",
           lambda: f"{len(_f('15_interface_controls.json')['comparators'])} complexes compared"))

    R_.append({"stage": "—", "module": "—", "method": "BioEmu-1",
               "name": "Conformational ensemble",
               "question": "Does one construct sample a masked and an unmasked state?",
               "status": "proposed",
               "result": "BioEmu-1, one construct at a time. Stage 14 superimposes "
                         "different designs on one backbone, which is a comparison "
                         "across designs and not an ensemble for a sequence."})
    R_.append({"stage": "—", "module": "—", "method": "no model exists",
               "name": "Venom protease specificity",
               "question": "Is the linker cut by a snake venom metalloproteinase?",
               "status": "proposed",
               "result": "No substrate dataset exists for SVMPs. The host "
                         "gelatinases stand in and the page says so."})
    R_.append({"stage": "—", "module": "—", "method": "BLI or SPR; cleavage time-course",
               "name": "Wet-lab measurement",
               "question": "Does any of this hold outside a predictor?",
               "status": "proposed",
               "result": "Expression, affinity by BLI or SPR, and a cleavage "
                         "time-course would test one scoring function each."})
    # grouped by module, pipeline order preserved inside each group
    order = {m: i for i, m in enumerate(
        ("01", "02a", "02b", "03", "04", "05", "06", "\u2014"))}
    return sorted(R_, key=lambda r: order.get(r.get("module", "\u2014"), 99))


REGISTER = _build_register()
N["rank_total"] = ranking_data["n_total"] if ranking_data else 0
N["rank_front"] = ranking_data["n_front"] if ranking_data else 0
N["rank_excluded"] = ranking_data.get("n_excluded", 0) if ranking_data else 0
N["rank_feasible"] = ranking_data.get("n_feasible", 0) if ranking_data else 0
N["rank_excluded_by"] = ranking_data.get("excluded_by", {}) if ranking_data else {}
RANK_CONSTRAINTS = ranking_data.get("constraints", []) if ranking_data else []
RANK_FRONT = ranking_data.get("front", []) if ranking_data else []
RANK_OBJECTIVES = ranking_data.get("objectives", []) if ranking_data else []
RANK_SHORTLIST = set(ranking_data.get("shortlist", [])) if ranking_data else set()
N["rank_shortlist"] = len(RANK_SHORTLIST)

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


masks_rows_full = "".join(
    f"<tr><td class='mono'>{k}</td><td class='mono'>{m['seq']}</td>"
    f"<td>{len(m['seq'])}</td><td>{m.get('method', 'graft')}</td>"
    f"<td>{m['kind']}</td><td class='dim'>{m['note']}</td></tr>"
    for k, m in lib["masks"].items())

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
{R.rebuild_banner(STALE, len(lib['constructs']))}
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
<tr><td>Acidic pH</td><td>Blood is pH {PH['blood_ph'][0]}&ndash;{PH['blood_ph'][1]}; the brief
gives roughly {PH['site_ph'][0]}&ndash;{PH['site_ph'][1]} at a bite site. The only side chain that
changes state across that range is histidine, so a pH switch is a histidine switch. At a pKa of
{PH['pka']:.1f} it is <b>{PH['blood_worst']:.1%}</b> protonated in blood and
<b>{PH['site_worst']:.1%}</b> at pH {PH['site_ph'][1]} &mdash; a {PH['fold_change']:.1f}-fold
change, but {PH['unprotonated_at_site']:.0%} of the residue is still neutral where it has to
act.</td><td>Not selected</td></tr>
<tr class="good"><td>Protease activity</td><td>Proteolysis is irreversible, which makes activation a
single committed event rather than an equilibrium, and venom metalloproteinases are themselves part
of the pathology being treated.</td><td>Selected</td></tr>
<tr><td>Hypoxia</td><td>Would require an oxygen-sensitive molecular mechanism suitable for a small
secreted protein.</td><td>Not selected</td></tr>
</tbody></table></div>


<div class="chart wide">{R.states_schematic()}</div>
</section>


<section class="sec">
<div class="stack"><div class="eyebrow">03</div><h2>Proposed computational pipeline</h2>
</div>
<div class="chart wide">{R.architecture_schematic()}</div>
<div class="eq">S(x) = [ S<sub>structure</sub>(x), S<sub>binding</sub>(x), S<sub>cleavage</sub>(x) ]</div>
<p>The three scores are not combined into a single value, since a high score on one objective
would otherwise compensate for failure on another. Selection runs in three steps. Hard constraints
remove constructs that cannot work at all, before any comparison, so that a failed requirement
cannot be offset by a strong score elsewhere. Pareto optimisation then identifies the nondominated
candidates among the rest, with no weights.</p>
<p>That leaves a set rather than an ordered list, and ordering it is a separate decision that has
to be stated. The criterion used here is <b>maximin</b>: the scores are normalised to a common
scale and each candidate is ranked by its lowest objective, so a design that is adequate on all
three comes before one that is excellent on two and poor on the third. The ordered front is then
deduplicated by mask, because predictions made for near-identical constructs carry the same
systematic errors and a shortlist of variants tests one hypothesis repeatedly.</p>
<div class="eq">R(x) = min<sub>i</sub> [ &#349;<sub>i</sub>(x) &minus; &lambda;&sigma;<sub>i</sub>(x) ]</div>
<p>The discount term is in the implementation and set to zero. It would down-weight a candidate
promoted by a single optimistic model, which is a real risk here, but no calibrated uncertainty
exists for any of the three scores and supplying a number would make the ordering look better
founded than it is.</p>
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">04</div><h2>Worked computational example</h2>
<p>The pipeline was implemented and applied to {meta['title']} ({PDB_ID}), a {N['resolution']} Å X-ray structure deposited in {meta['deposition_date'][:4]}.</p>
<p>It was chosen because it is a crystallographically characterised nanobody&ndash;toxin
complex with a single binding site, which is what the geometry and the masking measurements
need. The architecture and the methods are what transfer; a viper application would start
from a viper toxin complex and redesign the mask for it.</p>
<p>What each module produced, in the order the modules ran. The numbers are those in the diagram above, and rows marked &mdash; are the assembly and structure-prediction steps that every module depends on. Three steps discard candidates and are marked; the rest transform what they are given.</p>
{R.generation_table(COUNTS)}
<details class="tech"><summary><b>Every experiment the pipeline defines, under the module it
belongs to</b> &mdash; the question it asks, the method that produced its number, and what came
out</summary><div class="techbody">
{R.pipeline_table(REGISTER)}</div></details></div>

<h3 class="mod"><span>Module 01</span>Structural and affinity constraints</h3>
<div class="stack">
<p><b>Reference structure.</b> I used the experimental VHH–toxin complex (PDB {PDB_ID}) to
identify the antigen-binding residues and the possible linker attachment sites. The paratope
contains {N['n_paratope']} residues, and {termini_free}.</p></div>

<div class="two">
<div class="viewer">
<div class="vbar"><span class="eyebrow">Deposited structure · {PDB_ID}</span>
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
{chart("affinity_by_length",
   "<b>Predicted mask-affinity range for different linker lengths.</b> The green region shows "
   "combinations of linker length and mask affinity that allow at least "
   f"{100 - model['leak_max'] * 100:.0f}% masking before cleavage and "
   f"{model['recovery_min'] * 100:.0f}% antigen binding afterwards. The upper curve represents "
   "the weakest mask affinity that still provides sufficient masking while tethered; the lower "
   "line represents the strongest affinity that permits antigen binding after cleavage. The "
   f"widest range occurs at {N['linker_best']} residues. Triangles indicate linker lengths used "
   "in the computational designs.")}
{chart("usable_affinity",
   "The same constraint at the widest point, " + str(N["linker_best"]) + " residues.")}
</div>
</div>

<div class="stack">
<p><b>Estimating linker length and mask affinity.</b> I used established models from polymer
physics and binding equilibrium to estimate which linker lengths and mask affinities could allow
both masking and release.</p>

<p><b>1 &middot; Effective concentration of the tethered mask.</b> The flexible Gly/Ser linker was
approximated as a <i>worm-like chain</i>, the model commonly used to describe disordered protein
linkers. Which terminus carries it is set by the structure. The N-terminus sits
{N['n_term']:.1f} Å from the nearest paratope residue and {N['span_n']:.1f} Å from the paratope
centroid, against {N['c_term']:.1f} Å and {N['span_c']:.1f} Å for the C-terminus, so the
N-terminus is closer on both measures and is the attachment point used here.
From the centroid distance of {N['span_n']:.1f} Å I estimated the effective local concentration of
the mask, C<sub>eff</sub>, at every linker length, taking a contour length of
{model['rise_per_residue_A']:.1f} Å per residue and a persistence length of
{model['persistence_length_A']:.1f} Å.</p>
<div class="eq">⟨r²⟩ = 2·l<sub>p</sub>·L<sub>c</sub>·[1 − (l<sub>p</sub>/L<sub>c</sub>)(1 − e<sup>−L<sub>c</sub>/l<sub>p</sub></sup>)]
<br>C<sub>eff</sub>(d) = 1660.54 · (3 / 2π⟨r²⟩)<sup>3/2</sup> · e<sup>−3d²/2⟨r²⟩</sup></div>

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
gives the widest range, <b>{N['km_lo']:.0f} nM to {N['km_hi']:.0f} µM</b>. At that length the
chain's mean end-to-end distance of {N['ree_best']:.1f} Å almost matches the {N['span_n']:.1f} Å it
has to span, and C<sub>eff</sub> peaks at {N['ceff']:.1f} mM. A relatively weak mask could
therefore be sufficient, provided the linker holds it close to the paratope.</p>
<p>These estimates use a V<sub>H</sub>H–antigen K<sub>D</sub> of {N['vhh_toxin_kd_nM']:.1f} nM,
{model['toxin_conc_M'] * 1e9:.0f} nM antigen and {model['mask_conc_M'] * 1e6:.0f} µM released mask.
They are initial design constraints, not predictions of binding kinetics or of experimental
switching behaviour. The worm-like-chain treatment of a flexible linker follows
<a href="https://doi.org/10.1021/jp011355n">Zhou (2001)</a>; the use of an effective concentration
to describe what a disordered linker buys follows
<a href="https://doi.org/10.1073/pnas.1904813116">Sørensen and Kjaergaard (2019)</a> and
<a href="https://doi.org/10.1016/bs.mie.2020.09.012">Kjaergaard and colleagues (2021)</a>.</p>

<p><b>The window is wide but flat, so it does not pick a length.</b> It opens at about three
residues and changes little from twelve onwards: the range at 32 residues is within a factor of two
of the range at {N['linker_best']}. What the affinity model says is that a micromolar mask works at
any of these lengths. Which length to build is settled by whether the mask can physically reach the
paratope, which is measured on the folded constructs in Module 04, and there the answer is
different.</p>
<p><b>Which antigen affinity the floor rests on.</b> The release floor is set by how tightly the
nanobody holds its own antigen, so the measured value matters. The paper reports a biolayer
interferometry K<sub>D</sub> of <b>1.90 nM</b> for V<sub>H</sub>H5 a-sNTx against sNTx-1, and that
is the value used here. It does not state which of its numbered short-neurotoxin targets sNTx-1
is, and the structure is of the complex with short neurotoxin 1 (P01426), so the correspondence
is not confirmed. On the earlier assumed 1 nM the floor would be 99 nM instead of
{N['km_lo']:.0f} nM; the ceiling is unaffected. The distinction is kept explicit until the
correspondence is established.</p></div>


<h3 class="mod"><span>Module 02</span>Generation of candidate constructs</h3>

<h4 class="subh">Design of protease-cleavable linkers</h4>
<figure class="fig wide"><div class="figbody">{R.linker_schematic(LINKER_GEN, LENGTHS, SUB_LEN, REF_LEN)}</div>
<figcaption>How a protease-cleavable linker is designed. The sequence the protease has to
recognise and the distance the mask has to travel are set in separate steps, so either can
be changed without confounding the other.</figcaption></figure>
<div class="stack">
<p>I used CleaveNet to generate {R.word(SUB_LEN)}-residue peptides predicted to be preferentially
cleaved by {CLEAVE_TARGETS}, while limiting cleavage by {R.word(LINKER_GEN['n_off'])} other MMPs.</p>
<p>The model generated {LINKER_GEN['n_generated']:,} sequences. I ranked them by the difference
between their predicted cleavage scores for the target and off-target proteases, retaining
{R.word(LINKER_GEN['n_kept'])} candidates and adding {R.word(LINKER_GEN['n_controls'])} controls.</p>
<p>The {R.word(LINKER_GEN['n_controls'])} controls are a Gly/Ser linker with no designed cleavage motif
and a scrambled version of the highest-scoring sequence. They test whether the predicted
protease response depends on the designed sequence rather than on linker composition alone.
It does: the scramble has the same {R.word(SUB_LEN)} amino acids in a different order and scores
{CTRL_SCRAM:.2f} against {CN_BEST:.2f} for the sequence it was made from, and the Gly/Ser
control scores {CTRL_GS:.2f}.</p>
<p>The selected sequences were placed between flexible Gly/Ser segments to form complete
linkers. By varying the length of these segments I could assess linker geometry separately
from the predicted protease-recognition sequence.</p>
<p><b>What this does and does not show.</b> CleaveNet was trained on human MMPs, not on snake
venom metalloproteinases. These results demonstrate how to design a desired human-MMP cleavage
profile; they do not establish selectivity at a viper bite site.</p></div>

{mri("Masks were grafted from the antigen's own contact segments, with scrambled, truncated and "
     "poly-glycine controls. Linker substrates were generated under the cleavage profile above "
     "rather than fixed in advance, so cleavability constrains generation instead of filtering "
     "it afterwards. ProteinMPNN then redesigned the mask positions on each folded construct "
     "with the nanobody and linker held fixed.",
     f"{N['n_constructs']} constructs in the initial library and {N['n_redesign']} distinct "
     f"redesigned masks. Every redesign scored better than its parent graft by ProteinMPNN's own "
     f"measure, but that score correlated with paratope occlusion at r = {N['mpnn_corr']:+.2f}.",
     "ProteinMPNN generated sequence diversity, but sequence compatibility with a backbone is not "
     "a substitute for functional scoring and should not be used as one.")}

{chart_wide("paratope_footprint",
  f"Surface buried per nanobody residue on complex formation. Markers show the {N['n_near']} "
  f"residues with a heavy atom within {N['mask_cutoff']:g} Å of the antigen, which is the "
  f"set the mask is designed from.")}
<div class="stack">
<p><b>Three cutoffs, and each one is used for something different.</b> They are not alternative
answers to one question, so it matters which is which.</p></div>
<table class="tbl cuts">
<colgroup><col style="width:17%"><col style="width:22%"><col style="width:15%"><col style="width:46%"></colgroup>
<thead><tr><th>criterion</th><th>on the nanobody</th><th>on the antigen</th><th>what it is used for</th></tr></thead><tbody>
<tr><td class="mono">heavy atoms within {N['contact_cutoff']:g} &#8491;</td>
<td>{N['n_contacts']} residues &mdash; {N['contact_labels']}</td>
<td>{N['n_epitope_contacts']}</td>
<td>hydrogen-bond and salt-bridge range: the chemically interacting core. Reported so this
interface can be compared with others, and used for nothing else.</td></tr>
<tr class="hl"><td class="mono">heavy atoms within {N['mask_cutoff']:g} &#8491;</td>
<td>{N['n_near']} residues, marked on the plot above</td>
<td>{N['n_epitope_near']} residues, {N['mask_face']}</td>
<td><b>the contact face.</b> On the antigen side this is the sequence the mask is built from, in
Module 02a. Wider than the polar core because a mask reproduces a face, not a set of bonds.</td></tr>
<tr><td class="mono">loses more than 1 &#8491;<sup>2</sup> of surface</td>
<td>{N['n_paratope']} residues</td><td>{N['n_epitope']}</td>
<td><b>the footprint.</b> The surface a mask has to cover, so this is what occlusion is measured
against in Modules 03 and 04. Van der Waals packing reaches further than a hydrogen bond, so it is
the largest of the three.</td></tr>
</tbody></table>

{chart_wide("redesign_score_vs_occlusion", "ProteinMPNN score against paratope occlusion. "
   "A sequence that fits the backbone better does not occlude the paratope more.")}
<details class="tech"><summary>Every mask in the design space</summary><div class="techbody">
<div class="scroll"><table class="wrap"><thead><tr><th>mask</th><th>sequence</th><th>length</th>
<th>method</th><th>kind</th><th>note</th>
</tr></thead><tbody>{masks_rows_full}</tbody></table></div></div></details>

</section>

<section class="sec">
<div class="stack"><div class="eyebrow">05</div><h2>Scoring and computational results</h2>
<p>Each construct is scored on three independent requirements. Two of the three work, in the sense
that they separate designs from matched controls. The third does not, and that is the result this
section is mainly about.</p></div>

<h3 class="mod"><span>Module 03</span>Structural feasibility</h3>
{mri("Each construct was folded with ESMFold as one chain: mask, linker and nanobody "
     "together, not the mask alone. Preservation of the nanobody was assessed by "
     "superimposing its core on the separately folded free nanobody, and sequence-based "
     "descriptors were calculated over the whole construct for exposed hydrophobic "
     "surface and unpaired cysteines.",
     f"{N['n_folded']} constructs were folded and {N['n_ensemble']} superimposed on the nanobody. "
     f"The mean core RMSD against the free nanobody is {N['rmsd_mean']:.2f} Å. "
     + (f"The graft taken from the antigen's contact face leaves {N['feas_cys']} cysteine "
        f"unpaired in the predicted model and scores {N['feas_graft']:.2f} here."
        if N["has_graft_feas"] else
        "Per-mask feasibility is not shown: these scores were computed on the previous "
        "library and carry no row for the current graft."),
     "The nanobody fold is retained, which is a necessary condition and not evidence of "
     "expression or thermodynamic stability. Cysteine is the liability this module exists to "
     "catch, and it is measured from the predicted disulfides: an earlier version counted every "
     "cysteine past the nanobody's own pair as unpaired, which sent any cysteine-rich graft to "
     "zero on an assumption.")}
{chart("feasibility_components", "Feasibility components by construct family. A low score on any "
       "one component is enough to make a construct unusable. The dashed line is the unmasked "
       f"nanobody, which scores {N['ref_surface']:.2f} on surface liability itself, so that "
       "component is read as a change against it rather than on its own.")}
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
{R.alignment_control(N['boltz_ctrl_noalign'], N['boltz_ctrl'])}
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
<p>The linkers carried forward, and the controls they are read against. Target is the mean
predicted cleavage by {N['target_proteases']}; off-target is the mean over the membrane-type and
constitutive enzymes.</p>
<details class="tech"><summary>Every linker, with its predicted cleavage profile</summary>
<div class="techbody">{R.linker_table(LINKERS)}</div></details>
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
window of Module 01 is widest at {N['linker_best']} residues, and the linkers here are
{N['linker_len']}. That gap is deliberate and is a result in its own right: the tether model treats
the mask as a point at the end of the linker, and the mask is a peptide with its own dimensions
that has to reach across the paratope. Stage 16 measured the difference by holding the substrate
constant and lengthening only the inert flanks, and the extra length costs
{N['window_cost']:.2f} decades of the affinity window.</p>
<p><b>Both axes have to hold at once.</b> The lowest-scoring construct in the library,
{N['least_id']}, is physically reachable &mdash; its reach score of {N['least_ctx']:.2f} is
mid-pack &mdash; but its linker is a scramble of a good substrate, so recognition falls to
{N['least_sel']:+.2f} and the product with it, to {N['least_cleav']:.3f}. Same composition as
something cleavable, same accessibility, no site. That is the case a single averaged score would
have hidden.</p>

<h3 class="mod"><span>Module 06</span>Multiobjective ranking</h3>
{mri("Three steps, kept apart because they answer different questions. Hard constraints remove "
     "constructs that cannot work regardless of their other scores. The nondominated set is then "
     "taken over the three objectives, with no weights. That set is unordered, so turning it into "
     "a list needs a rule, and the rule is maximin on the normalised objectives: a construct is "
     "worth no more than its weakest requirement.",
     (f"{N['rank_excluded']} of {N['rank_total']} scored constructs fail a hard constraint. "
      f"Of the {N['rank_feasible']} that remain, {N['rank_front']} are nondominated, and "
      f"{N['rank_shortlist']} go forward for verification after deduplicating by mask."
      if N['rank_total'] else
      "Not yet run on this library; the stage has produced no output."),
     "The ordering is preliminary. Maximin assumes the three scores are comparable, and the "
     "normalisation that makes them comparable is within this run rather than against measured "
     "quantities. The masking objective also does not separate designed masks from controls, so "
     "a front built on it is not a list of constructs predicted to mask.")}

<p><b>What is removed before anything is compared.</b> A construct failing one of these is not a
worse candidate but not a candidate: no score on another axis can compensate, which is exactly
what a Pareto front would otherwise let it do.</p>
{R.constraint_table(RANK_CONSTRAINTS, N['rank_excluded_by'])}

<p><b>The ordering rule, stated.</b> The three scores are on unrelated scales, so they are min-max
normalised within the feasible set and each construct is ranked by its lowest normalised
objective. That prefers a design that is adequate three times over one that is excellent twice and
poor once, which is the right preference when all three requirements have to hold at once. An
uncertainty discount is implemented and set to zero: no calibrated uncertainty exists for any of
these scores, and inventing one would make the ordering look better founded than it is.</p>
{R.front_table(RANK_FRONT, RANK_OBJECTIVES, RANK_SHORTLIST)}

{chart("pareto_panels", "The three scores plotted pairwise, with the nondominated set marked.")}

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
</section>

<section class="sec">
<div class="stack"><div class="eyebrow">06</div><h2>Experimental collaboration and further development</h2>
<p>I would work with experimental PhD students to establish measurements corresponding to the
three scoring functions: protein stability and developability, mask interaction and dissociation,
and protease responsiveness. Purification behaviour and gel analysis address the first,
biolayer interferometry in the presence and absence of the protease addresses the second and
third.</p>
<p>These data would allow the predictive ability of individual computational methods to be
assessed and would identify which parts of the pipeline require development. As more measurements
accumulate, I would investigate whether sequence and structural representations of the constructs
can be used to train more accurate scoring functions and improve the selection of subsequent
candidates.</p>
<p>The preliminary implementation shows how structural information can constrain candidate
generation and how several molecular properties can be combined into a selection framework. The
principal remaining challenge is predicting masking and release: a tethered mask acts at a high
local concentration, so a weak isolated interaction may still be functional, and a static complex
prediction gives neither an affinity nor an off-rate. I would prioritise that component while
working with the experimental researchers to establish the measurements needed to assess
structural feasibility, binding and protease responsiveness. The aim is a method that can be
applied to other nanobody&ndash;antigen pairs without altering their binding interfaces.</p></div>
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
<script>{R.VIEWER_JS.replace("__MODELS__", json.dumps(models)).replace("__PARATOPE__", json.dumps(para_resi)).replace("__SOLVED__", json.dumps(f"structures/{PDB_ID}_{REFERENCE_VHH_CHAIN}{iface['toxin_chain']}.pdb")).replace("__ENSEMBLE__", json.dumps(ENSEMBLE_IDS)).replace("__TOPPARA__", json.dumps(top_paratope)).replace("__TOXCONTACT__", json.dumps(toxin_contact)).replace("__CDRS__", json.dumps(cdr_resi)).replace("__SCORES__", json.dumps(SCORES))}</script>
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
