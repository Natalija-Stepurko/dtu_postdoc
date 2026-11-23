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

    # drawn at the width it is displayed at. It now sits in a half-page column
    # beside the structure viewer, and at the old 8.6 inches the browser scaled
    # it to about half, taking its 6.5 pt labels down to roughly 4 pt.
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    for name in ("CDR1", "CDR2", "CDR3"):
        res = cdrs[name]["residues"]
        a, b = index[res[0]], index[res[-1]]
        ax.axvspan(a - .5, b + .5, color=CDR_COLOUR[name], alpha=0.08, zorder=0,
                   linewidth=0)
        ax.text((a + b) / 2, top * 1.13, name, ha="center", va="bottom",
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
                    xytext=(0, 9), ha="center", fontsize=7, color=INK)

    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.32),
              frameon=False, handlelength=1.1, columnspacing=1.1)
    ax.set_xlim(-1, len(delta)); ax.set_ylim(0, top * 1.34)
    ax.set_xlabel("VHH residue, N to C"); ax.set_ylabel("buried on binding  (\u00c5\u00b2)")
    ticks = list(range(0, len(labels), 10))
    ax.set_xticks(ticks); ax.set_xticklabels([labels[t] for t in ticks], fontsize=7)
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
    # only what the page shows. effective_concentration, design_window,
    # tradeoff, release_sensitivity and affinity_window all draw the same
    # constraint that affinity_by_length and usable_affinity now carry, so
    # building them left five figures published that nothing linked to. The
    # functions stay; re-add a call here if a page needs one again.
    made = [affinity_by_length(tether), usable_affinity(tether),
            paratope_footprint(ref, ref["cdrs"])]
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
        made.append(feasibility_components(load_json("12_feasibility.json"),
                                            load_json("03_constructs.json")))
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
        try:
            made.append(redesign_score_vs_occlusion(
                load_json("07b_library_scores.json"), fold, redesign,
                load_json("03_constructs.json")))
        except FileNotFoundError:
            pass
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


def redesign_score_vs_occlusion(scores, fold, redesign_fold, lib):
    """Does the design model's own score predict the thing the design is for?

    Every construct carrying a mask, not only the redesigns. The earlier
    version zipped twelve parents against six colours, so it drew 48 of 96
    points while reporting a correlation over all 96: the picture and the
    number described different sets. Both now cover the same points.

    The x value is ProteinMPNN's score for the construct's own mask sequence
    on its own folded backbone, averaged over decoding orders, so library
    constructs and redesigns sit on one scale.
    """
    from .mask_design import (MASK_FAMILY_COLOUR, MASK_FAMILY_SHORT,
                              family_by_construct)
    plt = _style()
    mask_of = {c["id"]: c.get("mask") for c in lib["constructs"]}
    occ, parent_of, mask_name = {}, {}, {}
    for src in (fold, redesign_fold):
        for r in (src or {}).get("results", []):
            occ[r["id"]] = r["occlusion"] * 100
            parent_of[r["id"]] = r.get("parent")
            mask_name[r["id"]] = r.get("mask")

    pts = {}
    for cid, rec in scores["scores"].items():
        if cid not in occ:
            continue
        fam = family_by_construct(mask_name.get(cid), parent_of.get(cid), mask_of)
        if fam is None:
            continue
        pts.setdefault(fam, {"x": [], "y": [], "redesign": []})
        pts[fam]["x"].append(rec["mask_score"])
        pts[fam]["y"].append(occ[cid])
        pts[fam]["redesign"].append(rec["source"] == "redesign")

    fig, ax = plt.subplots(figsize=(8.6, 3.8))
    fig.subplots_adjust(right=0.70)
    order = [f for f, _ in MASK_FAMILY_SHORT.items() if f in pts]
    for fam in order:
        d = pts[fam]
        colour = MASK_FAMILY_COLOUR[fam]
        x = np.array(d["x"]); y = np.array(d["y"]); red = np.array(d["redesign"])
        if (~red).any():
            ax.scatter(x[~red], y[~red], s=26, color=colour, marker="o",
                       edgecolor="white", linewidth=.7, alpha=.9,
                       label=MASK_FAMILY_SHORT[fam], zorder=3)
        if red.any():
            ax.scatter(x[red], y[red], s=30, color=colour, marker="^",
                       edgecolor="white", linewidth=.7, alpha=.9,
                       label=None if (~red).any() else MASK_FAMILY_SHORT[fam],
                       zorder=3)

    xs = np.array([v for d in pts.values() for v in d["x"]])
    ys = np.array([v for d in pts.values() for v in d["y"]])
    n_red = sum(sum(d["redesign"]) for d in pts.values())
    if xs.size > 2:
        r = float(np.corrcoef(xs, ys)[0, 1])
        ax.annotate(f"Pearson r = {r:+.2f}   n = {xs.size}"
                    f"  ({xs.size - n_red} library, {n_red} redesigns)",
                    xy=(.02, .96), xycoords="axes fraction", fontsize=7.5,
                    color=INK2, va="top")

    handles, labels = ax.get_legend_handles_labels()
    shape = [plt.Line2D([], [], marker=m, ls="none", color=INK3, ms=5.4,
                        markeredgecolor="white", markeredgewidth=.7, label=t)
             for m, t in (("o", "library construct"), ("^", "ProteinMPNN redesign"))]
    leg = ax.legend(handles + shape, labels + [h.get_label() for h in shape],
                    fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, handletextpad=.5, labelspacing=.65,
                    borderaxespad=0)
    leg.set_title("mask family", prop={"size": 7.5, "weight": "semibold"})
    leg._legend_box.align = "left"

    # the spread belongs in the caption, not over the points it describes
    ax.set_xlabel("ProteinMPNN score for the mask  (lower is better)")
    ax.set_ylabel("paratope occluded  (%)")
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


def feasibility_components(feas, lib=None):
    """Where each construct family loses feasibility.

    Grouped and coloured by the same six mask families as the score figure,
    so a family keeps one colour across the page. It grouped by individual
    mask before, which gave fourteen bars per component in colours that meant
    nothing anywhere else.
    """
    from .mask_design import (MASK_FAMILY_COLOUR, MASK_FAMILY_SHORT,
                              family_of, MASK_FAMILIES)
    plt = _style()
    rows = feas["results"]
    #: the unmasked nanobody carries no mask, so it has no family. It is the
    #: baseline the masked constructs are read against, not a design.
    REFERENCE = "free nanobody (reference)"
    mask_of = {c["id"]: c.get("mask") for c in (lib or {}).get("constructs", [])}

    def family(row):
        name = row["mask"]
        if name is None:
            return REFERENCE
        if str(name).startswith("mpnn_"):
            # a redesign's mask is named for its parent construct
            name = mask_of.get(str(name)[5:], name)
        return family_of(name)

    fams = {}
    for r in rows:
        key = family(r)
        if key is None:
            continue
        fams.setdefault(key, []).append(r)

    keys = ("f_fold", "f_preserved", "f_surface", "f_cysteine")
    # named for the quantity each one measures rather than for the liability
    # it stands in for, so a reader can match a bar to the method panel
    labels = ("ESMFold\npLDDT", "V$_H$H fold\npreservation",
              "exposed hydrophobic\nresidues", "unpaired\ncysteines")
    order = [q for q, _ in MASK_FAMILIES if q in fams]

    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    fig.subplots_adjust(right=0.70)
    x = np.arange(len(keys))
    width = 0.82 / max(len(order), 1)
    for i, fam in enumerate(order):
        vals = [float(np.mean([r[k] for r in fams[fam]])) for k in keys]
        pos = x + i * width - 0.41 + width / 2
        colour = MASK_FAMILY_COLOUR[fam]
        ax.bar(pos, vals, width * 0.9, color=colour, linewidth=0,
               label=MASK_FAMILY_SHORT[fam])
        # a bar of height zero draws nothing, and zero is the value that
        # disqualifies a construct, so it has to read as a value rather than
        # as missing data
        for px, v in zip(pos, vals):
            if v < 0.01:
                ax.plot([px], [0.012], marker="_", markersize=3.2,
                        markeredgewidth=1.6, color=colour, zorder=4)
                ax.annotate("0", (px, 0.035), ha="center", va="bottom",
                            fontsize=5.6, color=colour)

    if REFERENCE in fams:
        ref = [float(np.mean([r[k] for r in fams[REFERENCE]])) for k in keys]
        for j, v in enumerate(ref):
            ax.hlines(v, j - 0.46, j + 0.46, color="#12191F", linewidth=1.4,
                      linestyle=(0, (3, 2)), zorder=5,
                      label="free nanobody" if j == 0 else None)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("component score")
    ax.set_ylim(0, 1.05)
    leg = ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, handlelength=1.3, handletextpad=.6,
                    labelspacing=.7, borderaxespad=0)
    leg.set_title("mask family", prop={"size": 7.5, "weight": "semibold"})
    leg._legend_box.align = "left"
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


#: Washes for the three regions of the affinity axis. Each region also carries
#: its own text label, so the reading never rests on colour alone.
BAD_WASH, BAD_INK = "#F6E4E4", "#9C3B3B"
OK_WASH, OK_INK = "#E2F0E8", "#1F6B4A"


def _si_molar(value):
    """A molar quantity as a number and a unit, e.g. 2.06e-7 -> '206 nM'."""
    for scale, unit in ((1e-3, "mM"), (1e-6, "µM"), (1e-9, "nM"),
                        (1e-12, "pM")):
        if value >= scale:
            n = value / scale
            return (f"{n:.0f} {unit}" if n >= 10 or abs(n - round(n)) < 0.05
                    else f"{n:.1f} {unit}")
    return f"{value:.1e} M"


def affinity_by_length(tether, built_lengths=(12, 16, 20, 24, 28, 32)):
    """The design constraint across every linker length, in three regions.

    The release floor is set by competition with the antigen and does not
    depend on the tether, so it is one horizontal line. The masking ceiling
    is the effective concentration the tether produces, which rises as the
    linker reaches the paratope and falls again as it grows slack. Between
    them a mask does both jobs; outside them it does neither.
    """
    plt = _style()
    rows = sorted((r for r in tether["scan"] if r["terminus"] == "N"),
                  key=lambda r: r["linker_residues"])
    rows = [r for r in rows if r["linker_residues"] <= 48]
    n = np.array([r["linker_residues"] for r in rows], float)
    floor = np.array([r["km_min_for_release_M"] for r in rows], float)
    ceil = np.array([r["km_max_for_masking_M"] for r in rows], float)
    best = tether["best_per_terminus"]["N"]
    lo_y, hi_y = 1e-9, 1e-1

    fig, ax = plt.subplots(figsize=(5.4, 4.1))
    ax.grid(False)
    ax.set_yscale("log")
    ax.set_xlim(n.min(), n.max())
    ax.set_ylim(lo_y, hi_y)

    open_ = ceil > floor
    ax.fill_between(n, lo_y, floor, color=BAD_WASH, lw=0, zorder=1)
    ax.fill_between(n, np.maximum(ceil, floor), hi_y, color=BAD_WASH, lw=0, zorder=1)
    ax.fill_between(n, floor, ceil, where=open_, color=OK_WASH, lw=0, zorder=1,
                    interpolate=True)
    ax.plot(n, floor, color=BAD_INK, lw=1.6, zorder=3)
    ax.plot(n[open_], ceil[open_], color=BAD_INK, lw=1.6, zorder=3)

    ax.text(n.max() - 1.0, np.sqrt(lo_y * floor[0]), "too tight",
            ha="right", va="center", color=BAD_INK, fontsize=9, zorder=4)
    ax.text(n.max() - 1.0, 2.6e-2, "too weak", ha="right",
            va="center", color=BAD_INK, fontsize=9, zorder=4)
    mid = n.searchsorted(26)
    ax.text(n[mid], np.sqrt(floor[mid] * ceil[mid]), "usable range",
            ha="center", va="center", color=OK_INK, fontsize=12,
            fontweight="bold", zorder=4)
    ax.text(n[mid], np.sqrt(floor[mid] * ceil[mid]) / 7.0,
            "masks, then releases", ha="center",
            va="center", color=OK_INK, fontsize=8, zorder=4)

    b = best["linker_residues"]
    ax.plot([b, b], [floor[0], best["km_max_for_masking_M"]], color=INK2, lw=1.1,
            ls=(0, (4, 3)), zorder=5)
    ax.annotate(f"widest at {b}\n{_si_molar(best['km_min_for_release_M'])}"
                f" to {_si_molar(best['km_max_for_masking_M'])}",
                xy=(b, best["km_max_for_masking_M"]), xytext=(b + 3.0, 1.6e-3),
                color=INK, fontsize=8.6, zorder=6,
                arrowprops=dict(arrowstyle="-", color=INK2, lw=0.9))

    for x in built_lengths:
        ax.plot([x], [lo_y * 1.45], marker="^", ms=5, color=CHROME,
                clip_on=False, zorder=6)
    ax.text(float(np.mean(built_lengths)), lo_y * 2.6,
            "lengths built", ha="center", va="bottom", color=CHROME,
            fontsize=8, zorder=6)

    ax.set_xlabel("linker length (residues)", color=INK3, fontsize=8.5)
    ax.set_ylabel("mask dissociation constant K$_D$", color=INK3, fontsize=8.5)
    ax.set_yticks([1e-9, 1e-7, 1e-5, 1e-3, 1e-1])
    ax.set_yticklabels([_si_molar(t) for t in (1e-9, 1e-7, 1e-5, 1e-3, 1e-1)])
    ax.tick_params(labelsize=8, colors=INK3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, "affinity_by_length")


def usable_affinity(tether):
    """The design constraint as one horizontal axis, split into three regions.

    Everything the reader needs is the position of two boundaries on a
    dissociation-constant axis: below the lower one the mask would not come
    off after cleavage, above the upper one it would not mask in the first
    place, and between them it does both.
    """
    plt = _style()
    best = tether["best_per_terminus"]["N"]
    lo = best["km_min_for_release_M"]
    hi = best["km_max_for_masking_M"]
    left, right = 1e-9, 1e-1

    fig, ax = plt.subplots(figsize=(5.4, 2.1))
    ax.grid(False)
    ax.set_xscale("log")
    ax.set_xlim(left, right)
    ax.set_ylim(0, 1)

    bot, top = 0.42, 0.92
    for x0, x1, wash, edge in ((left, lo, BAD_WASH, BAD_INK),
                               (lo, hi, OK_WASH, OK_INK),
                               (hi, right, BAD_WASH, BAD_INK)):
        ax.fill_between([x0, x1], bot, top, color=wash, lw=0, zorder=1)
        ax.plot([x0, x0], [bot, top], color="white", lw=2.0, zorder=3)

    def region(x0, x1, head, sub, ink, head_size=11):
        mid = 10 ** ((np.log10(x0) + np.log10(x1)) / 2)
        ax.text(mid, 0.70, head, ha="center", va="center", color=ink,
                fontsize=head_size, fontweight="bold", zorder=4)
        ax.text(mid, 0.54, sub, ha="center", va="center", color=ink,
                fontsize=7.4, zorder=4)

    region(left, lo, "Too tight", "stays bound", BAD_INK, 9.5)
    region(lo, hi, "Usable range", "masks, then releases", OK_INK, 11)
    region(hi, right, "Too weak", "does not mask", BAD_INK, 9.5)

    for x in (lo, hi):
        ax.plot([x, x], [bot - 0.09, top], color=INK2, lw=1.1, ls=(0, (4, 3)),
                zorder=5)
        ax.text(x, bot - 0.14, _si_molar(x), ha="center", va="top", color=INK,
                fontsize=9.5, fontweight="bold", zorder=5)

    ax.set_yticks([])
    for side in ("left", "top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_position(("data", 0.06))
    ax.spines["bottom"].set_color(INK3)
    ticks = [10 ** e for e in range(-9, 0)]
    ax.set_xticks(ticks)
    ax.set_xticklabels([_si_molar(t) if t in (1e-9, 1e-5, 1e-1) else ""
                        for t in ticks])
    ax.tick_params(axis="x", length=3, pad=3, labelsize=8, colors=INK3)
    ax.set_xlabel("mask dissociation constant K$_D$ (log scale) "
                  "→ weaker binding", color=INK3, fontsize=8.5, labelpad=6)
    return _save(fig, "usable_affinity")
