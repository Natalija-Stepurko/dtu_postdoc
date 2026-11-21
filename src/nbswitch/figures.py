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


def paratope_footprint(ref, cdrs, chain="A"):
    plt = _style()
    data = ref["chain_pairs"][chain]
    delta = np.array(data["vhh_delta_sasa_all"])
    labels = data["vhh_residue_labels"]
    index = {l: i for i, l in enumerate(labels)}
    x = np.arange(len(delta))
    fig, ax = plt.subplots(figsize=(8.6, 2.5))
    for name in ("CDR1", "CDR2", "CDR3"):
        res = cdrs[name]["residues"]
        a, b = index[res[0]], index[res[-1]]
        ax.axvspan(a - .5, b + .5, color=SEQ_WASH, zorder=0, linewidth=0)
        ax.text((a + b) / 2, delta.max() * 1.10, name, ha="center", va="bottom",
                fontsize=7.5, color=INK3)
    paratope = set(data["paratope_labels"])
    hit = np.array([l in paratope for l in labels])
    ax.bar(x[~hit], delta[~hit], width=.8, color="#C9D1D6", linewidth=0)
    ax.bar(x[hit], delta[hit], width=.8, color=SEQ, linewidth=0)
    for label, value in sorted(data["paratope_delta_sasa"].items(), key=lambda kv: -kv[1])[:5]:
        ax.annotate(label, (index[label], value), textcoords="offset points",
                    xytext=(0, 3), ha="center", fontsize=6.5, color=INK)
    ax.set_xlim(-1, len(delta)); ax.set_ylim(0, delta.max() * 1.28)
    ax.set_xlabel("VHH residue, N to C"); ax.set_ylabel("buried on binding  (Å²)")
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
    fig, ax = plt.subplots(figsize=(5.0, 3.0))
    for n_res, colour in ((5, SEQ), (20, CHROME)):
        row = next(r for r in tether["scan"]
                   if r["terminus"] == "N" and r["linker_residues"] == n_res)
        c_eff = row["c_eff_M"]
        y = c_eff / (c_eff + km)
        ax.plot(km * 1e9, y, lw=2, color=colour)
        ax.annotate(f"masked · {n_res}-residue linker", (km[195] * 1e9, y[195]),
                    textcoords="offset points", xytext=(4, 5), fontsize=7.5, color=colour)
    bound = ((model["toxin_conc_M"] / model["KT_M"])
             / (1 + model["toxin_conc_M"] / model["KT_M"] + model["mask_conc_M"] / km))
    ax.plot(km * 1e9, bound, lw=2, color=STR)
    ax.annotate("toxin-bound after cleavage", (km[255] * 1e9, bound[255]),
                textcoords="offset points", xytext=(-8, -15), fontsize=7.5,
                color=STR, ha="right")
    for level, text in ((1 - model["leak_max"], "99% masked"),
                        (model["recovery_min"], "90% recovered")):
        ax.axhline(level, lw=1.1, color=NULL, ls=(0, (2, 3)))
        ax.text(1.5e-2, level, f" {text}", fontsize=6.8, color=INK3, va="center")
    ax.set_xscale("log"); ax.set_ylim(0, 1.06)
    ax.set_xlabel("mask–VHH $K_D$  (nM)"); ax.set_ylabel("fraction of nanobody")
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
    fold = None
    try:
        fold = load_json("04_folding.json")
        made += [occlusion_by_linker(fold), occlusion_vs_confidence(fold)]
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
