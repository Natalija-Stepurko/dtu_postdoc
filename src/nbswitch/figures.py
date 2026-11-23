"""Charts, drawn to match the house style: flat colour, mono labels, no chart junk.

Figures are emitted as SVG so they stay sharp and inline well in the report.
"""
from __future__ import annotations
import numpy as np

from .config import FIGURES, load_json

#: House palette. Sequence-side blue, structure-side orange, muted grey, teal accent.
#:
#: The teal was #0E7C7B, which carries a chroma of 0.088 and so falls under the
#: 0.1 floor the palette checker sets: it reads as grey. That is a problem
#: wherever it is a series beside the grey the same figures use for a
#: reference or a total. Lifted to the same hue at a chroma that clears the
#: floor; all three categorical slots then pass on every check, worst CVD
#: separation 13.8 (deutan) and worst normal-vision separation 21.8.
SEQ, STR, NULL, CHROME = "#2D5BD1", "#C06014", "#6E7880", "#0B9488"
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


#: The matched comparison: one antigen-derived mask against the three
#: controls built from it. Every one of the four was paired with the same 21
#: linkers, so the only thing that differs along the axis is the mask.
#:
#: The axis carries the library's mask ids, not descriptions of them, so that
#: a column on the plot and a row in the table beside it are the same name.
OCCLUSION_CONTROLS = (("epi_core", "epi_core", SEQ),
                      ("epi_core_scram", "epi_core_scram", STR),
                      ("epi_core_half", "epi_core_half", "#8A9BA8"),
                      ("polyG", "polyG", NULL))


def occlusion_by_group(fold, masks=OCCLUSION_CONTROLS,
                       name="occlusion_by_group", label_fn=None):
    """Whether covering the paratope depends on the antigen-derived sequence.

    One column per mask rather than per design category. Averaging nine
    different designs into a single "antigen-derived" bar compared a group of
    nine against groups of one, which is not the comparison the controls were
    built for: each control is derived from epi_core, so epi_core is what it
    has to be read against.

    Every construct is a point, because the spread is the result - the
    designed mask sits inside it.
    """
    plt = _style()
    rows = fold["results"]
    groups = [(label, colour, [r["occlusion"] * 100 for r in rows if r["mask"] == m])
              for m, label, colour in masks]
    groups = [g for g in groups if g[2]]
    if not groups:
        raise ValueError("no construct matches any mask in the comparison; "
                         "the names have drifted from the library")

    # four columns leave room for an upright label under each one and for the
    # mean to be written beside its line; nine do not, and both then overprint
    # their neighbours. Past five the labels go on an angle and the means move
    # on top of their own lines. The figure is kept narrow deliberately: it is
    # served into one half of a two-column grid, where a wide figure is scaled
    # down and takes its type with it.
    crowded = len(groups) > 5
    fig, ax = plt.subplots(figsize=(5.6, 4.6) if crowded else (6.4, 3.2))
    rng = np.random.default_rng(42)
    for i, (label, colour, vals) in enumerate(groups):
        jitter = rng.uniform(-0.17, 0.17, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=16, color=colour,
                   alpha=0.45, linewidth=0, zorder=2)
        m = float(np.mean(vals))
        ax.hlines(m, i - 0.32, i + 0.32, color=colour, linewidth=2.2, zorder=3)
        if crowded:
            ax.annotate(f"{m:.1f}%", (i, m), fontsize=7.5, color=colour,
                        ha="center", va="bottom", textcoords="offset points",
                        xytext=(0, 3), fontweight="semibold")
        else:
            ax.annotate(f"{m:.1f}%", (i + 0.36, m), fontsize=8.5, color=colour,
                        va="center", fontweight="semibold")
    ax.axhline(0, color=INK3, linewidth=0.8, alpha=0.4, zorder=1)
    ax.set_xticks(range(len(groups)))
    if crowded:
        ax.set_xticklabels([g[0].replace("\n", " ") for g in groups],
                           fontsize=7, rotation=35, ha="right",
                           rotation_mode="anchor")
    else:
        ax.set_xticklabels([g[0].replace(" ", "\n", 1) for g in groups],
                           fontsize=8)
    ax.set_ylabel("paratope occluded  (%)")
    # the uncrowded form writes each mean to the right of its line, so the
    # last column needs room past it; the crowded form writes it on top
    ax.set_xlim(-0.6, len(groups) - 0.4 if crowded else len(groups) - 0.2)
    return _save(fig, name)


def occlusion_by_mask(fold):
    """Every antigen-derived mask, for the panel behind the main comparison."""
    from .mask_design import MASK_FAMILY_COLOUR, MASK_LABELS, family_of
    kinds = {r["mask"]: r["kind"] for r in fold["results"] if r["mask"]}
    grafts = sorted(m for m, k in kinds.items() if k == "epitope graft")
    masks = tuple((m, MASK_LABELS.get(m, m).replace(", ", "\n"),
                   MASK_FAMILY_COLOUR.get(family_of(m), NULL)) for m in grafts)
    return occlusion_by_group(fold, masks, "occlusion_by_mask")


#: The three series the length comparison is built on: one antigen-derived
#: mask against the two controls that remove its sequence and then any
#: designed interaction at all. Grey reads as grey to the palette checker,
#: which is legal here because every line carries its own label.
LENGTH_SERIES_MASKS = (("epi_core", "main contact loop", SEQ),
                       ("epi_core_scram", "scrambled control", STR),
                       ("polyG", "inert control", NULL))


def occlusion_by_linker(folding, motif="cn01", series=LENGTH_SERIES_MASKS,
                        name="occlusion_by_linker"):
    """Paratope occlusion against linker length, one cleavage motif throughout.

    The mask names here were the library's previous ones - epi1, epi2,
    epi1_half, epi1_scram - and only polyG survived the rename, so the figure
    had been drawing a single family out of twelve while its caption said
    "by mask". It also joined every construct of a family in linker-length
    order regardless of cleavage motif, which drew vertical jumps between
    constructs that share only an x value.

    Holding the motif fixed makes each line one mask across six lengths, so
    the only thing varying along a line is the length.
    """
    plt = _style()
    rows = folding["results"]
    # three lines can each carry their name at their end; eleven cannot. Past
    # four the end labels are stacked apart so far that they leave the axes,
    # and they are long enough to run off the right of a half-column figure,
    # so the names move to a legend below the plot instead.
    many = len(series) > 4
    fig, ax = plt.subplots(figsize=(5.4, 4.6) if many else (5.4, 3.2))
    drawn, ends = 0, []
    for mask, label, colour in series:
        pts = sorted((r for r in rows
                      if r["mask"] == mask and r["id"].startswith(f"{mask}_{motif}")),
                     key=lambda r: r["linker"])
        if not pts:
            continue
        drawn += 1
        x = [p["linker"] for p in pts]
        y = [p["occlusion"] * 100 for p in pts]
        ax.plot(x, y, lw=2, marker="o", ms=5, color=colour,
                markeredgecolor="white", markeredgewidth=1.1, zorder=3,
                label=label.replace("\n", " "))
        ends.append((y[-1], x[-1], label, colour))
    if not drawn:
        raise ValueError(f"no construct matches the series for motif {motif!r}; "
                         "the mask names have drifted from the library again")
    if not many:
        # two series can finish within a line's height of each other, so the
        # end labels are pushed apart rather than drawn on top of one another
        ends.sort()
        span = (max(e[0] for e in ends) - min(e[0] for e in ends)
                if len(ends) > 1 else 1)
        gap = max(span * 0.11, 0.9)
        placed = []
        for value, xend, label, colour in ends:
            at = value if not placed else max(value, placed[-1] + gap)
            placed.append(at)
            ax.annotate(label, (xend, value), xytext=(xend + 1.2, at),
                        textcoords="data", fontsize=8, color=colour,
                        va="center",
                        arrowprops=(dict(arrowstyle="-", color=colour, lw=0.8,
                                         shrinkA=2, shrinkB=0)
                                    if abs(at - value) > 0.35 else None))
    ax.axhline(0, lw=1.1, color=INK3, alpha=0.4, ls=(0, (2, 3)))
    lens = sorted({r["linker"] for r in rows if r["linker"]})
    ax.set_xticks(lens)
    ax.set_xlim(min(lens) - 2, max(lens) + (1 if many else 9))
    ax.set_xlabel("linker length  (residues)")
    ax.set_ylabel("paratope occluded  (%)")
    if many:
        # several lines share a colour because the colour is the family, so
        # the legend is folded to one entry per name. Listing eleven masks
        # against six colours would promise an identification the figure
        # cannot make.
        handles, labels = ax.get_legend_handles_labels()
        seen, h, lab = set(), [], []
        for handle, text in zip(handles, labels):
            if text in seen:
                continue
            seen.add(text)
            h.append(handle)
            lab.append(text)
        ax.legend(h, lab, loc="upper center", bbox_to_anchor=(0.5, -0.22),
                  ncol=2, frameon=False, fontsize=7.5, handlelength=1.6,
                  columnspacing=1.2, labelspacing=0.5, borderaxespad=0)
    return _save(fig, name)


def occlusion_by_linker_all(folding, motif="cn01"):
    """The same comparison across every mask that has a length series.

    Lines are named for their family, not for the mask, because the colour is
    the family: three substitution masks are one green, and a legend that gave
    them three names would say they can be told apart on the plot.
    """
    from .mask_design import MASK_FAMILY_COLOUR, MASK_FAMILY_SHORT, family_of
    rows = folding["results"]
    masks = sorted({r["mask"] for r in rows if r["mask"]
                    and any(x["mask"] == r["mask"] and x["id"].startswith(f"{r['mask']}_{motif}")
                            and x["linker"] != 24 for x in rows)})
    series = tuple((m, MASK_FAMILY_SHORT.get(family_of(m), m),
                    MASK_FAMILY_COLOUR.get(family_of(m), NULL)) for m in masks)
    return occlusion_by_linker(folding, motif, series, "occlusion_by_linker_all")


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
        made += [occlusion_by_group(fold), occlusion_by_mask(fold),
                 occlusion_by_linker(fold), occlusion_by_linker_all(fold)]
    except FileNotFoundError:
        pass
    try:
        trig = load_json("10_trigger.json")
        # the panel colours by mask family, which the trigger stage does not
        # record: the library names each construct's mask and the two folding
        # stages name each redesign's parent
        panel = cleavability_panel(trig, load_json("03_constructs.json"),
                                   load_json("04_folding.json"),
                                   load_json("08_redesign_folding.json"))
        for fig in (panel, cleavability_components(trig)):
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
        # redesign_landscape and occlusion_vs_confidence are not on the page:
        # the first was never placed, the second lost its home when the
        # folding section that held it turned out never to be inserted.
        # The functions stay; add a call back here if a page needs one.
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


#: How a linker is drawn on the cleavability panel. The shape is the linker,
#: because the colour is already spoken for by the mask family.
LINKER_MARKER = (("generated", "o", "generated motif"),
                 ("protease-resistant control", "s", "Gly/Ser control"),
                 ("scrambled control", "^", "scrambled control"))


def cleavability_panel(trigger, lib=None, fold=None, redesign_fold=None):
    """Recognition against reach: the two halves of the cleavability score.

    A linker has to be both recognised by the protease and physically
    reachable in the folded construct. Plotting them against each other shows
    whether the two agree, and the controls say whether either axis is
    discriminating at all.

    Every construct is drawn, coloured by mask family like the rest of the
    page. The colour used to be the linker kind, which is the one thing the
    marker shape already says, so the same six families were six different
    colours here and three flat ones there.

    The trigger stage records no mask, so the family is resolved through the
    construct library and the two folding stages, the same route the score
    figure takes.
    """
    from .mask_design import (MASK_FAMILY_COLOUR, MASK_FAMILY_SHORT,
                              MASK_FAMILIES, family_by_construct)
    plt = _style()
    rows = [r for r in trigger["results"] if r.get("selectivity") is not None]
    if not rows:
        return None

    mask_of = {c["id"]: c.get("mask") for c in (lib or {}).get("constructs", [])}
    fam_of = {}
    for source in (fold, redesign_fold):
        for r in (source or {}).get("results", []):
            fam_of[r["id"]] = family_by_construct(r.get("mask"),
                                                  r.get("parent"), mask_of)
    unplaced = [r["id"] for r in rows if fam_of.get(r["id"]) is None]
    if unplaced:
        raise ValueError(
            f"{len(unplaced)} constructs have no mask family on the "
            f"cleavability panel, starting with {unplaced[0]!r}; the folding "
            "results no longer cover what the trigger stage scored")

    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    fig.subplots_adjust(right=0.70)
    # selectivity is a property of the linker, and there are five linkers, so
    # 348 constructs fall in five vertical lines. Drawn on the exact value the
    # columns are solid and the colour cannot be read, which is the whole
    # point of the figure, so the points are spread by a small random offset.
    # It is smaller than half the gap between the two closest columns, so no
    # point crosses into a neighbour.
    rng = np.random.default_rng(7)
    for kind, marker, _label in LINKER_MARKER:
        items = [r for r in rows if r.get("linker_kind") == kind]
        if not items:
            continue
        # one call per marker with the colours carried per point, and the
        # order shuffled, so no family is systematically drawn over another
        order = rng.permutation(len(items))
        x = np.array([items[i]["selectivity"] for i in order])
        y = np.array([items[i]["context_score"] for i in order])
        c = [MASK_FAMILY_COLOUR[fam_of[items[i]["id"]]] for i in order]
        ax.scatter(x + rng.uniform(-0.022, 0.022, len(x)), y,
                   s=26, color=c, marker=marker, edgecolor="white",
                   linewidth=.6, alpha=.85, zorder=3)

    ax.axvline(0, color=RULE, lw=1, zorder=1)
    ax.set_xlabel("sequence recognition  (CleaveNet selectivity)")
    ax.set_ylabel("structural accessibility")

    from matplotlib.lines import Line2D
    present = [f for f, _ in MASK_FAMILIES
               if any(fam_of[r["id"]] == f for r in rows)]
    fams = ax.legend(
        handles=[Line2D([], [], marker="o", linestyle="none", markersize=5,
                        color=MASK_FAMILY_COLOUR[f],
                        label=MASK_FAMILY_SHORT[f]) for f in present],
        fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=False, handlelength=1.3, handletextpad=.6, labelspacing=.7,
        borderaxespad=0)
    fams.set_title("mask family", prop={"size": 7.5, "weight": "semibold"})
    fams._legend_box.align = "left"
    ax.add_artist(fams)

    kinds = ax.legend(
        handles=[Line2D([], [], marker=m, linestyle="none", markersize=5,
                        color=INK3, label=lab) for _, m, lab in LINKER_MARKER],
        fontsize=7, loc="lower left", bbox_to_anchor=(1.02, 0.0),
        frameon=False, handlelength=1.3, handletextpad=.6, labelspacing=.7,
        borderaxespad=0)
    kinds.set_title("linker", prop={"size": 7.5, "weight": "semibold"})
    kinds._legend_box.align = "left"
    return _save(fig, "cleavability_panel")


#: The three measurements accessibility is built from, in the order they are
#: described on the page.
CONTEXT_PARTS = (("extension", SEQ, "extension"),
                 ("exposure", CHROME, "exposure"),
                 ("clearance_open_fraction", STR, "clearance"))


def cleavability_components(trigger):
    """The three accessibility requirements, over the whole library.

    Every construct is drawn, ordered along x by the product the three make,
    so the product is the rising edge and each component can be read against
    it. It used to show the twenty-eight best by name, which told the reader
    what the best looked like and nothing about where the library sits or
    which of the three is the binding constraint.

    Three requirements multiplied means the product is bounded above by the
    smallest of them, and clearance is the smallest almost everywhere: the
    figure's job is to make that visible rather than asserted.
    """
    plt = _style()
    rows = sorted(trigger["results"], key=lambda r: r["context_score"])
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(8.6, 3.4))
    fig.subplots_adjust(right=0.78)
    x = np.arange(1, len(rows) + 1)

    ax.fill_between(x, [r["context_score"] for r in rows], color="#E3E9EC",
                    linewidth=0, zorder=0, label="product")
    for key, colour, label in CONTEXT_PARTS:
        ax.plot(x, [r[key] for r in rows], lw=1.3, color=colour, label=label,
                solid_capstyle="round", zorder=3)

    ax.set_xlim(1, len(rows))
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(f"all {len(rows)} constructs, ordered by the product")
    ax.set_ylabel("score")
    leg = ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=False, handlelength=1.3, handletextpad=.6,
                    labelspacing=.7, borderaxespad=0)
    leg._legend_box.align = "left"
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
