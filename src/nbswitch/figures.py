"""Charts, drawn to match the house style: flat colour, mono labels, no chart junk.

Figures are emitted as SVG so they stay sharp and inline well in the report.
"""
from __future__ import annotations
import numpy as np

from .config import FIGURES, load_json

#: House palette. Sequence-side blue, structure-side orange, muted grey, teal accent.
SEQ, STR, NULL, CHROME = "#2D5BD1", "#C06014", "#6E7880", "#0E7C7B"
INK, INK2, INK3 = "#12191F", "#3A464F", "#6E7880"
RULE, SEQ_WASH, STR_WASH = "#DDE4E7", "#EAEFFB", "#FAEEE5"
#: DejaVu first so matplotlib has real metrics; the SVG keeps the whole stack as
#: a font-family string, which the browser resolves to the nicest available.
MONO = ["DejaVu Sans Mono", "SF Mono", "Menlo", "Consolas", "monospace"]


def _style():
    import logging, matplotlib
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "svg.fonttype": "none", "font.family": MONO, "font.size": 8.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": True, "axes.spines.bottom": True,
        "axes.edgecolor": RULE, "axes.linewidth": 1.0,
        "axes.labelcolor": INK3, "axes.labelsize": 8,
        "xtick.color": INK3, "ytick.color": INK3,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 1.0, "ytick.major.width": 1.0,
        "axes.grid": True, "grid.color": "#EBF0F2", "grid.linewidth": 1.0,
        "legend.frameon": False, "figure.facecolor": "none",
        "axes.facecolor": "none", "savefig.facecolor": "none",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.08,
    })
    return plt


def _save(fig, name):
    path = FIGURES / f"{name}.svg"
    fig.savefig(path, format="svg")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return path


CDR_COLOUR = {"CDR1": "#2E7DBE", "CDR2": "#B8860B", "CDR3": "#A8324A"}
FRAMEWORK_COLOUR = "#7E8A93"


def paratope_footprint(ref, cdrs, chain="A"):
    """Surface buried per residue, each bar coloured by the region it sits in."""
    plt = _style()
    data = ref["chain_pairs"][chain]
    delta = np.array(data["vhh_delta_sasa_all"])
    labels = data["vhh_residue_labels"]
    index = {l: i for i, l in enumerate(labels)}
    x = np.arange(len(delta))
    top = delta.max()

    region = {}
    for name in ("CDR1", "CDR2", "CDR3"):
        for r in cdrs[name]["residues"]:
            region[r] = name

    paratope = set(data["paratope_labels"])
    # the set the mask is designed from, so the plot and the mask agree
    contacts = set(data.get("paratope_near_labels")
                   or data.get("paratope_contact_labels", []))

    fig, ax = plt.subplots(figsize=(8.6, 2.7))
    for name in ("CDR1", "CDR2", "CDR3"):
        res = cdrs[name]["residues"]
        a, b = index[res[0]], index[res[-1]]
        ax.axvspan(a - .5, b + .5, color=CDR_COLOUR[name], alpha=0.08, zorder=0,
                   linewidth=0)
        ax.text((a + b) / 2, top * 1.11, name, ha="center", va="bottom",
                fontsize=7.5, color=CDR_COLOUR[name])

    # one bar colour per region, pale where the residue buries nothing
    for name, colour in list(CDR_COLOUR.items()) + [("framework", FRAMEWORK_COLOUR)]:
        sel = np.array([region.get(l, "framework") == name for l in labels])
        buried = sel & np.array([l in paratope for l in labels])
        if (sel & ~buried).any():
            ax.bar(x[sel & ~buried], delta[sel & ~buried], width=.8, color=colour,
                   alpha=0.22, linewidth=0)
        if buried.any():
            ax.bar(x[buried], delta[buried], width=.8, color=colour, linewidth=0,
                   label=name if name != "framework" else "framework")

    for label in contacts:
        if label in index:
            i = index[label]
            ax.plot([i], [delta[i] + top * 0.05], marker="v", markersize=4.5,
                    color="#12191F", linestyle="none", zorder=5)
    if contacts:
        cut = data.get("mask_contact_cutoff_A",
                       data.get("contact_cutoff_A", 3.5))
        ax.plot([], [], marker="v", markersize=4.5, color="#12191F", linestyle="none",
                label=f"within {cut:g} \u00c5")

    for label, value in sorted(data["paratope_delta_sasa"].items(),
                               key=lambda kv: -kv[1])[:5]:
        ax.annotate(label, (index[label], value), textcoords="offset points",
                    xytext=(0, 9), ha="center", fontsize=6.5, color=INK)

    ax.legend(fontsize=6.8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.30),
              frameon=False, handlelength=1.1, columnspacing=1.2)
    ax.set_xlim(-1, len(delta)); ax.set_ylim(0, top * 1.30)
    ax.set_xlabel("VHH residue, N to C"); ax.set_ylabel("buried on binding  (\u00c5\u00b2)")
    ticks = list(range(0, len(labels), 10))
    ax.set_xticks(ticks); ax.set_xticklabels([labels[t] for t in ticks], fontsize=6.5)
    return _save(fig, "paratope_footprint")


def effective_concentration(tether):
    plt = _style()
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    for term, colour, label in (("N", SEQ, "N-terminal"), ("C", STR, "C-terminal")):
        rows = [r for r in tether["scan"] if r["terminus"] == term]
        n = [r["linker_residues"] for r in rows]
        y = [r["c_eff_M"] * 1e3 for r in rows]
        ax.plot(n, y, lw=2, color=colour, solid_capstyle="round")
        j = n.index(40)
        ax.annotate(label, (n[j], y[j]), textcoords="offset points",
                    xytext=(4, 6 if term == "N" else -12), fontsize=7.5, color=colour)
    ax.set_yscale("log"); ax.set_xlim(0, 60)
    ax.set_xlabel("linker length  (residues)")
    ax.set_ylabel("$C_{eff}$  (mM)")
    return _save(fig, "effective_concentration")


def design_window(tether):
    plt = _style()
    floor = tether["km_min_for_release_M"] * 1e9
    lo, hi_lim = 1e-1, 1e7
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.0), sharey=True)
    for ax, term, colour, wash, title in (
            (axes[0], "N", SEQ, SEQ_WASH, f"N-terminal · spans {tether['spans_A']['N']} Å"),
            (axes[1], "C", STR, STR_WASH, f"C-terminal · spans {tether['spans_A']['C']} Å")):
        rows = [r for r in tether["scan"] if r["terminus"] == term]
        n = np.array([r["linker_residues"] for r in rows])
        ceiling = np.array([r["km_max_for_masking_M"] for r in rows]) * 1e9
        ax.fill_between(n, np.clip(np.full_like(ceiling, floor), lo, hi_lim),
                        np.clip(ceiling, lo, hi_lim), where=ceiling > floor,
                        color=wash, linewidth=0)
        ax.plot(n, np.clip(ceiling, lo * .9, hi_lim), lw=2, color=colour)
        ax.axhline(floor, lw=1.4, color=NULL, ls=(0, (5, 3)))
        ax.set_yscale("log"); ax.set_ylim(lo, hi_lim); ax.set_xlim(0, 60)
        ax.set_xlabel("linker length  (residues)")
        ax.set_title(title, fontsize=8, color=INK2, pad=7)
        closed = [r["linker_residues"] for r in rows if not r["window_open"]]
        if closed:
            ax.axvspan(0, max(closed) + .5, color="#F1F3F4", zorder=0, linewidth=0)
            ax.text((max(closed) + .5) / 2, hi_lim * .18, "no\nwindow", ha="center",
                    va="top", fontsize=7, color=INK3)
        best = max((r for r in rows if r["window_open"]), key=lambda r: r["window_decades"])
        ax.annotate(f"{best['window_decades']:.1f} decades",
                    (best["linker_residues"], min(best["km_max_for_masking_M"] * 1e9, hi_lim)),
                    textcoords="offset points", xytext=(9, -16), fontsize=7.5, color=colour)
    axes[0].set_ylabel("mask–VHH $K_D$  (nM)")
    axes[0].annotate("masking ceiling", (30, 2.2e5), fontsize=7.5, color=SEQ)
    axes[0].annotate("release floor", (30, floor * 2.2), fontsize=7.5, color=INK3)
    return _save(fig, "design_window")


def tradeoff(tether):
    plt = _style()
    model = tether["model"]
    km = np.logspace(-11, -2, 400)
    fig, ax = plt.subplots(figsize=(5.6, 3.3))
    handles = []
    for n_res, colour in ((5, SEQ), (20, CHROME)):
        row = next(r for r in tether["scan"]
                   if r["terminus"] == "N" and r["linker_residues"] == n_res)
        c_eff = row["c_eff_M"]
        line, = ax.plot(km * 1e9, c_eff / (c_eff + km), lw=2, color=colour,
                        label=f"masked, {n_res}-residue linker")
        handles.append(line)
    bound = ((model["toxin_conc_M"] / model["KT_M"])
             / (1 + model["toxin_conc_M"] / model["KT_M"] + model["mask_conc_M"] / km))
    line, = ax.plot(km * 1e9, bound, lw=2, color=STR,
                    label="toxin-bound after cleavage")
    handles.append(line)
    for level in (1 - model["leak_max"], model["recovery_min"]):
        ax.axhline(level, lw=1, color=NULL, ls=(0, (2, 3)), zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 1e7)
    ax.set_ylim(0, 1.18)
    ax.set_yticks([0, .25, .5, .75, 1.0])
    # thresholds labelled in the margin, clear of every curve
    ax.text(1.4e7, 1 - model["leak_max"], f"{(1-model['leak_max'])*100:.0f}% masked",
            fontsize=7, color=INK3, va="center")
    ax.text(1.4e7, model["recovery_min"], f"{model['recovery_min']*100:.0f}% recovered",
            fontsize=7, color=INK3, va="center")
    ax.set_xlabel("mask\u2013VHH $K_D$  (nM)")
    ax.set_ylabel("fraction of nanobody")
    ax.legend(handles=handles, fontsize=7.2, loc="lower center",
              bbox_to_anchor=(0.46, -0.02), ncol=1, handlelength=1.6)
    return _save(fig, "tradeoff")


def release_sensitivity(tether):
    plt = _style()
    fig, ax = plt.subplots(figsize=(5.4, 3.0))
    colours = {1e-10: SEQ, 1e-9: STR, 1e-8: CHROME}
    for kt in (1e-10, 1e-9, 1e-8):
        pts = [s for s in tether["sensitivity"] if s["KT_M"] == kt]
        x = [p["T_M"] * 1e9 for p in pts]
        y = [p["km_min_for_release_M"] * 1e9 for p in pts]
        ax.plot(x, y, lw=2, marker="o", ms=4.5, color=colours[kt],
                markeredgecolor="white", markeredgewidth=1)
        ax.annotate(f"$K_D$ {kt*1e9:g} nM", (x[-1], y[-1]), textcoords="offset points",
                    xytext=(-4, 8), fontsize=7.5, color=colours[kt], ha="right")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("free toxin  (nM)"); ax.set_ylabel("minimum mask $K_D$  (nM)")
    return _save(fig, "release_sensitivity")


def occlusion_by_linker(folding):
    plt = _style()
    rows = folding["results"]
    order = ["epi1", "epi2", "epi1_half", "epi1_scram", "polyG"]
    colours = {"epi1": SEQ, "epi2": CHROME, "epi1_half": "#8A9BA8",
               "epi1_scram": STR, "polyG": NULL}
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    for mask in order:
        pts = sorted([r for r in rows if r["mask"] == mask], key=lambda r: r["linker"])
        if not pts:
            continue
        x = [p["linker"] for p in pts]
        y = [p["occlusion"] * 100 for p in pts]
        ax.plot(x, y, lw=2, marker="o", ms=4.5, color=colours[mask],
                markeredgecolor="white", markeredgewidth=1)
        ax.annotate(mask, (x[-1], y[-1]), textcoords="offset points", xytext=(5, 0),
                    fontsize=7.5, color=colours[mask], va="center")
    ax.axhline(0, lw=1.1, color=NULL, ls=(0, (2, 3)))
    ax.set_xlim(2, 38)
    ax.set_xlabel("linker length  (residues)")
    ax.set_ylabel("paratope occluded  (%)")
    return _save(fig, "occlusion_by_linker")


def occlusion_vs_confidence(folding):
    plt = _style()
    rows = folding["results"]
    colours = {"epitope graft": SEQ, "scrambled control": STR,
               "truncated control": "#8A9BA8", "inert control": NULL}
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    for kind, colour in colours.items():
        pts = [r for r in rows if r["kind"] == kind]
        if not pts:
            continue
        ax.scatter([p["plddt_vhh"] for p in pts], [p["occlusion"] * 100 for p in pts],
                   s=34, color=colour, edgecolor="white", linewidth=1.1, label=kind)
    ax.axhline(0, lw=1.1, color=NULL, ls=(0, (2, 3)))
    ax.set_xlabel("VHH-region pLDDT"); ax.set_ylabel("paratope occluded  (%)")
    ax.legend(fontsize=7, loc="best")
    return _save(fig, "occlusion_vs_confidence")


def build_all():
    ref = load_json("01_reference.json")
    tether = load_json("02_tether.json")
    made = [paratope_footprint(ref, ref["cdrs"]),
            effective_concentration(tether), design_window(tether),
            tradeoff(tether), release_sensitivity(tether)]
    window = affinity_window(tether)
    if window is not None:
        made.append(window)
    fold = None
    try:
        fold = load_json("04_folding.json")
        made += [occlusion_by_linker(fold), occlusion_vs_confidence(fold)]
    except FileNotFoundError:
        pass
    try:
        trig = load_json("10_trigger.json")
        for fig in (cleavability_panel(trig), cleavability_components(trig)):
            if fig is not None:
                made.append(fig)
    except FileNotFoundError:
        pass
    try:
        made.append(feasibility_components(load_json("12_feasibility.json")))
    except FileNotFoundError:
        pass
    try:
        made.append(pareto_panels(load_json("13_ranking.json")))
    except FileNotFoundError:
        pass
    try:
        trig = load_json("10_trigger.json")
    except FileNotFoundError:
        pass
    try:
        redesign = load_json("08_redesign_folding.json")
        if fold is not None:
            made.append(redesign_landscape(redesign, fold))
        made.append(redesign_score_vs_occlusion(redesign))
    except FileNotFoundError:
        pass
    return made


def redesign_landscape(redesign_fold, base_fold):
    """Where the redesigned masks land against the grafts and the controls."""
    plt = _style()
    rows = redesign_fold["results"]
    base = base_fold["results"]
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    groups = [("inert control", NULL, "o"), ("scrambled control", STR, "o"),
              ("truncated control", "#8A9BA8", "o"), ("epitope graft", SEQ, "o")]
    for kind, colour, marker in groups:
        pts = [r for r in base if r["kind"] == kind]
        if pts:
            ax.scatter([p["plddt_vhh"] for p in pts], [p["occlusion"] * 100 for p in pts],
                       s=38, color=colour, marker=marker, edgecolor="white",
                       linewidth=1.1, label=kind, zorder=3)
    ax.scatter([r["plddt_vhh"] for r in rows], [r["occlusion"] * 100 for r in rows],
               s=26, color=CHROME, marker="D", edgecolor="white", linewidth=.9,
               label="ProteinMPNN redesign", zorder=2, alpha=.85)
    ax.axhline(0, lw=1.1, color=NULL, ls=(0, (2, 3)))
    ax.set_xlabel("VHH-region pLDDT"); ax.set_ylabel("paratope occluded  (%)")
    ax.legend(fontsize=7, loc="best")
    return _save(fig, "redesign_landscape")


def redesign_score_vs_occlusion(redesign_fold):
    """Does the design model's own score predict the thing we care about?"""
    plt = _style()
    rows = redesign_fold["results"]
    parents = sorted({r["parent"] for r in rows})
    colours = [SEQ, STR, CHROME, "#8A9BA8", "#7B5EA7", "#2E7D32"]
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    for parent, colour in zip(parents, colours):
        pts = [r for r in rows if r["parent"] == parent]
        ax.scatter([p["mpnn_score"] for p in pts], [p["occlusion"] * 100 for p in pts],
                   s=30, color=colour, edgecolor="white", linewidth=.9, label=parent)
    xs = np.array([r["mpnn_score"] for r in rows])
    ys = np.array([r["occlusion"] * 100 for r in rows])
    if len(xs) > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        ax.annotate(f"Pearson r = {r:+.2f}  (n = {len(xs)})",
                    xy=(.03, .95), xycoords="axes fraction", fontsize=7.5, color=INK2,
                    va="top")
    ax.set_xlabel("ProteinMPNN score  (lower is better)")
    ax.set_ylabel("paratope occluded  (%)")
    ax.legend(fontsize=6.5, loc="best", ncol=2)
    return _save(fig, "redesign_score_vs_occlusion")


def pareto_panels(ranking):
    """The three objectives, pairwise, with the non-dominated designs marked."""
    plt = _style()
    rows = ranking["rows"]
    objs = ranking["objectives"]
    if len(objs) < 2:
        return None
    pairs = [(objs[i], objs[j]) for i in range(len(objs)) for j in range(i + 1, len(objs))]
    fig, axes = plt.subplots(1, len(pairs), figsize=(3.1 * len(pairs), 3.0))
    axes = np.atleast_1d(axes)
    for ax, (a, b) in zip(axes, pairs):
        dom = [r for r in rows if not r.get("on_front")]
        front = [r for r in rows if r.get("on_front")]
        ax.scatter([r[a] for r in dom], [r[b] for r in dom], s=16, color="#C9D1D6",
                   linewidth=0, label="dominated", zorder=2)
        ax.scatter([r[a] for r in front], [r[b] for r in front], s=34, color=SEQ,
                   edgecolor="white", linewidth=1.1, label="on the front", zorder=3)
        ax.set_xlabel(a); ax.set_ylabel(b)
    axes[0].legend(fontsize=7, loc="best")
    fig.tight_layout()
    return _save(fig, "pareto_panels")


def feasibility_components(feas):
    """Where each construct family loses feasibility."""
    plt = _style()
    rows = feas["results"]
    #: the unmasked nanobody carries no mask, so it has no family name. It is
    #: the baseline the masked constructs are read against, not a design, and
    #: is drawn as a reference rather than as one more bar.
    REFERENCE = "free nanobody (reference)"
    fams = {}
    for r in rows:
        if r["mask"] is None:
            key = REFERENCE
        elif str(r["mask"]).startswith("mpnn"):
            key = "ProteinMPNN"
        else:
            key = str(r["mask"])
        fams.setdefault(key, []).append(r)
    keys = ("f_fold", "f_preserved", "f_surface", "f_cysteine")
    labels = ("fold\nconfidence", "nanobody\npreserved", "surface\nliability", "cysteine\npairing")
    designs = [k for k in fams if k != REFERENCE]
    order = sorted(designs, key=lambda k: -np.mean([x["feasibility_score"] for x in fams[k]]))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    x = np.arange(len(keys))
    width = 0.8 / max(len(order), 1)
    colours = [SEQ, CHROME, STR, "#8A9BA8", "#7B5EA7", "#2E7D32", NULL]
    for i, fam in enumerate(order):
        vals = [np.mean([r[k] for r in fams[fam]]) for k in keys]
        pos = x + i * width - 0.4 + width / 2
        colour = colours[i % len(colours)]
        ax.bar(pos, vals, width * 0.92, color=colour, linewidth=0, label=fam)
        # a bar of height zero draws nothing, and zero is the value that
        # disqualifies a construct, so it has to be visible as a value rather
        # than read as missing data
        for px, v in zip(pos, vals):
            if v < 0.01:
                ax.plot([px], [0.012], marker="_", markersize=3.2,
                        markeredgewidth=1.6, color=colour, zorder=4)
                ax.annotate("0", (px, 0.035), ha="center", va="bottom",
                            fontsize=5.6, color=colour)
    if REFERENCE in fams:
        ref = [np.mean([r[k] for r in fams[REFERENCE]]) for k in keys]
        for j, v in enumerate(ref):
            ax.hlines(v, j - 0.44, j + 0.44, color="#12191F", linewidth=1.4,
                      linestyle=(0, (3, 2)), zorder=5,
                      label=REFERENCE if j == 0 else None)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("component score"); ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.22))
    return _save(fig, "feasibility_components")


def cleavability_panel(trigger):
    """Recognition against reach: the two halves of the cleavability score.

    A linker has to be both recognised by the protease and physically
    reachable in the folded construct. Plotting them against each other shows
    whether the two agree, and the controls say whether either axis is
    discriminating at all.
    """
    plt = _style()
    rows = [r for r in trigger["results"] if r.get("selectivity") is not None]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(5.2, 3.6))

    groups = {}
    for r in rows:
        groups.setdefault(r.get("linker_kind", "generated"), []).append(r)
    style = {"generated": (SEQ, "o", "generated"),
             "protease-resistant control": (NULL, "s", "Gly/Ser control"),
             "scrambled control": (CHROME, "^", "scrambled control")}
    for kind, items in groups.items():
        colour, marker, label = style.get(kind, (STR, "D", kind))
        ax.scatter([r["selectivity"] for r in items],
                   [r["context_score"] for r in items],
                   s=34, c=colour, marker=marker, linewidths=0, alpha=.85,
                   label=label, zorder=3)

    ax.axvline(0, color=RULE, lw=1, zorder=1)
    ax.set_xlabel("sequence recognition  (selectivity, target − off-target)")
    ax.set_ylabel("reach in the construct")
    ax.legend(fontsize=7, loc="best", frameon=False)
    return _save(fig, "cleavability_panel")


def cleavability_components(trigger):
    """The three necessary conditions, per construct, sorted by their product."""
    plt = _style()
    rows = sorted(trigger["results"], key=lambda r: -r["context_score"])[:28]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    x = np.arange(len(rows))
    for key, colour, label in (("extension", SEQ, "extension"),
                               ("exposure", CHROME, "exposure"),
                               ("clearance_open_fraction", STR, "clearance")):
        ax.plot(x, [r[key] for r in rows], lw=1.6, color=colour, label=label,
                solid_capstyle="round")
    ax.bar(x, [r["context_score"] for r in rows], width=.7, color="#E3E9EC",
           zorder=0, label="product")
    ax.set_xticks(x)
    ax.set_xticklabels([r["id"] for r in rows], rotation=90, fontsize=5.5)
    ax.set_ylabel("score"); ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, ncol=4, loc="upper right", frameon=False)
    return _save(fig, "cleavability_components")


def affinity_window(tether):
    """How the usable mask affinity depends on the nanobody's own affinity.

    One y-axis, not two. The floor and the ceiling are both mask dissociation
    constants, so they belong on the same scale, and the window is the gap
    between them rather than a separate series needing its own axis. Drawn
    this way the quantity that matters is the height of the band, and the
    point where it closes is where the two lines meet.
    """
    scan = tether.get("affinity", {}).get("scan")
    if not scan:
        return None
    plt = _style()

    kd = np.array([r["vhh_toxin_kd_M"] for r in scan]) * 1e9
    floor = np.array([r["km_min_for_release_M"] for r in scan]) * 1e9
    ceiling = np.array([r["km_max_for_masking_M"] for r in scan]) * 1e9
    used = tether["affinity"]["used_M"] * 1e9
    closes = tether["affinity"]["window_closes_above_M"] * 1e9

    fig, ax = plt.subplots(figsize=(6.2, 3.8))

    open_ = floor < ceiling
    ax.fill_between(kd, floor, ceiling, where=open_, color=SEQ_WASH,
                    linewidth=0, zorder=1)
    ax.plot(kd, floor, color=STR, linewidth=2, zorder=3,
            label="release floor: the mask must be looser than this")
    ax.plot(kd, ceiling, color=CHROME, linewidth=2, zorder=3,
            label="masking ceiling: the mask must be tighter than this")

    # the value the rest of the page uses
    ax.axvline(used, color=INK, linewidth=1, linestyle=(0, (3, 2)), zorder=2)
    ax.annotate(f"used here, {used:g} nM", (used, ceiling[0] * 0.55),
                ha="center", va="top", fontsize=7, color=INK,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="none", alpha=0.9))

    # where the band pinches shut
    if closes < kd.max():
        ax.axvline(closes, color=STR, linewidth=1, linestyle=(0, (1, 2)), zorder=2)
        ax.annotate(f"window closes\nat {closes:.0f} nM", (closes, ceiling[0] * 0.35),
                    ha="left", va="top", fontsize=7, color=STR,
                    xytext=(5, 0), textcoords="offset points")

    # the width of the band, which is the number the table reported
    finite = [i for i, (lo, hi) in enumerate(zip(floor, ceiling)) if lo < hi]
    last_open = finite[-1] if finite else 0
    for i, (r, x, lo, hi) in enumerate(zip(scan, kd, floor, ceiling)):
        if r["window_decades"] != r["window_decades"] or lo >= hi:
            continue
        label_here = (abs(x - used) / used < 0.01 or i == 0 or i == last_open)
        if label_here:
            ax.annotate(f"{r['window_decades']:.1f} decades", (x, (lo * hi) ** 0.5),
                        ha="center", va="center", fontsize=7, color=INK2,
                        bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                                  edgecolor="none", alpha=0.88))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("V$_H$H–antigen K$_D$  (nM)")
    ax.set_ylabel("usable mask K$_D$  (nM)")
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, 1.17),
              frameon=False, handlelength=1.6, columnspacing=1.6)
    ax.set_xlim(kd.min(), kd.max())
    ax.set_ylim(floor.min() * 0.45, ceiling.max() * 3)
    return _save(fig, "affinity_window")
