"""Stage 3. Figures. Every value plotted is read from the stage JSON outputs."""
import numpy as np, json
from common import load_json, FIGS, TOL, style
plt = style()
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

iface = load_json("s02_interface.json")
teth = load_json("s03_tether.json")
C = TOL

# ---------- Figure 1: paratope footprint ----------
ref = iface["chain_pairs"]["A"]
labels = ref["vhh_residue_labels"]
d = np.array(ref["vhh_delta_sasa_all"])
x = np.arange(len(d))
cdrs = iface["cdrs"]
pos = {l: i for i, l in enumerate(labels)}

fig, ax = plt.subplots(figsize=(9.0, 3.0))
for name in ("CDR1", "CDR2", "CDR3"):
    r = cdrs[name]["residues"]
    a, b = pos[r[0]], pos[r[-1]]
    ax.axvspan(a - 0.5, b + 0.5, color=C["band"], zorder=0)
    ax.text((a + b) / 2, d.max() * 1.07, name, ha="center", va="bottom",
            fontsize=8, color=C["ink2"])
isp = np.array([l in set(ref["paratope_labels"]) for l in labels])
ax.bar(x[~isp], d[~isp], width=0.78, color=C["grid"], linewidth=0)
ax.bar(x[isp], d[isp], width=0.78, color=C["blue"], linewidth=0)
top = sorted(ref["paratope_delta_sasa"].items(), key=lambda kv: -kv[1])[:5]
for lbl, val in top:
    ax.annotate(lbl, (pos[lbl], val), textcoords="offset points", xytext=(0, 3),
                ha="center", fontsize=7, color=C["ink"])
ax.set_xlim(-1, len(d))
ax.set_ylim(0, d.max() * 1.22)
ax.set_xlabel("VHH residue (Kabat numbering, N → C)")
ax.set_ylabel("buried on binding (Å²)")
ax.set_title(f"Paratope footprint of VHH TPL1158_01_C09 on α-cobratoxin  ·  PDB 9GCN chain A:X")
ticks = [pos[l] for l in labels if l.endswith(("1", "0")) and len(l) <= 4][::3]
ax.set_xticks(ticks); ax.set_xticklabels([labels[t] for t in ticks], fontsize=6.5, rotation=90)
fig.savefig(FIGS / "fig1_paratope.png"); plt.close(fig)

# ---------- Figure 2: schematic of the three states ----------
fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.5))
def vhh(ax, cx, cy, highlight=True):
    ax.add_patch(FancyBboxPatch((cx - .55, cy - .42), 1.1, .84,
        boxstyle="round,pad=0.05,rounding_size=0.12", fc="#e9eef5", ec=C["ink2"], lw=1.0))
    ax.text(cx, cy - .02, "VHH", ha="center", va="center", fontsize=9, color=C["ink"])
    ax.plot([cx - .30, cx + .30], [cy + .40, cy + .40], lw=3.2,
            color=C["blue"] if highlight else C["grid"], solid_capstyle="round")
def mask(ax, cx, cy, col):
    ax.add_patch(Circle((cx, cy), .30, fc=col, ec="white", lw=1.6))
    ax.text(cx, cy, "M", ha="center", va="center", fontsize=8.5, color="white")
def toxin(ax, cx, cy):
    ax.add_patch(FancyBboxPatch((cx - .30, cy - .26), .60, .52,
        boxstyle="round,pad=0.04,rounding_size=0.18", fc=C["orange"], ec="white", lw=1.6))
    ax.text(cx, cy, "T", ha="center", va="center", fontsize=8.5, color="white")

titles = ["A — masked", "B — cleaved", "C — toxin bound"]
for ax, t in zip(axes, titles):
    ax.set_xlim(0, 3); ax.set_ylim(0, 2.3); ax.axis("off")
    ax.set_title(t, fontsize=9.5, color=C["ink"])
ax = axes[0]; vhh(ax, 1.5, .9); mask(ax, 1.5, 1.62, C["aqua"])
ax.plot([1.2, 1.2, 1.35], [1.3, 1.52, 1.62], lw=1.4, color=C["ink2"])
ax.text(0.95, 1.46, "L", fontsize=8, color=C["ink2"])
ax.text(1.5, .18, "paratope occluded", ha="center", fontsize=7.5, color=C["ink2"])
ax = axes[1]; vhh(ax, 1.5, .9); mask(ax, 2.35, 1.72, C["aqua"])
ax.add_patch(FancyArrowPatch((1.62, 1.45), (2.08, 1.66), arrowstyle="-|>",
    mutation_scale=9, lw=1.2, color=C["ink2"]))
ax.text(1.12, 1.44, "✂ protease", fontsize=7.5, color=C["orange"])
ax.text(1.5, .18, "tether removed; mask diffuses", ha="center", fontsize=7.5, color=C["ink2"])
ax = axes[2]; vhh(ax, 1.35, .9); toxin(ax, 1.35, 1.66); mask(ax, 2.62, .55, C["aqua"])
ax.text(1.5, .18, "paratope available", ha="center", fontsize=7.5, color=C["ink2"])
fig.suptitle("Design states of the N–L–M construct", fontsize=10, y=1.02)
fig.savefig(FIGS / "fig2_states.png"); plt.close(fig)

# ---------- Figure 3: effective concentration ----------
scan = teth["scan"]
fig, ax = plt.subplots(figsize=(5.4, 3.4))
for term, col, lab in (("N", C["blue"], "N-terminal fusion"), ("C", C["orange"], "C-terminal fusion")):
    s = [r for r in scan if r["terminus"] == term]
    n = [r["linker_residues"] for r in s]; y = [r["c_eff_M"] * 1e3 for r in s]
    ax.plot(n, y, lw=2.0, color=col, label=lab)
    j = n.index(30)
    ax.annotate(lab, (n[j], y[j]), textcoords="offset points",
                xytext=(6, 10 if term == "N" else -16), fontsize=8, color=col)
ax.set_yscale("log"); ax.set_xlabel("linker length (residues)")
ax.set_ylabel("effective concentration $C_{eff}$ (mM)")
ax.set_title("Tethered mask: effective concentration at the paratope")
ax.legend(loc="lower right", fontsize=8)
fig.savefig(FIGS / "fig3_ceff.png"); plt.close(fig)

# ---------- Figure 4: the design window ----------
fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
floor = teth["km_min_for_release_M"] * 1e9
YLO, YHI = 1e-1, 1e7
for ax, term, col, name in ((axes[0], "N", C["blue"], "N-terminal fusion  (span 11.8 Å)"),
                            (axes[1], "C", C["orange"], "C-terminal fusion  (span 42.0 Å)")):
    s_ = [r for r in scan if r["terminus"] == term]
    n = np.array([r["linker_residues"] for r in s_])
    hi = np.array([r["km_max_for_masking_M"] for r in s_]) * 1e9
    lo = np.full_like(hi, floor)
    ax.fill_between(n, np.clip(lo, YLO, YHI), np.clip(hi, YLO, YHI),
                    where=hi > lo, color=col, alpha=0.16, linewidth=0)
    ax.plot(n, np.clip(hi, YLO * 0.9, YHI), lw=2.0, color=col, clip_on=True)
    ax.axhline(floor, lw=1.5, color=C["ink2"], ls=(0, (5, 3)))
    ax.set_yscale("log"); ax.set_ylim(YLO, YHI)
    ax.set_xlim(0, 60); ax.set_xlabel("linker length (residues)")
    ax.set_title(name, fontsize=9.5)
    closed = [r["linker_residues"] for r in s_ if not r["window_open"]]
    if closed:
        ax.axvspan(0, max(closed) + 0.5, color="#efeeec", zorder=0)
        ax.text((max(closed) + 0.5) / 2, YHI * 0.25, "window\nclosed", ha="center",
                va="top", fontsize=7.5, color=C["ink2"])
    ok = [r for r in s_ if r["window_open"]]
    b = max(ok, key=lambda r: r["window_decades"])
    ax.annotate(f"widest {b['window_decades']:.1f} decades\nat {b['linker_residues']} residues",
                (b["linker_residues"], min(b["km_max_for_masking_M"] * 1e9, YHI)),
                textcoords="offset points", xytext=(10, -26), fontsize=7.5, color=col)
axes[0].set_ylabel("mask–VHH $K_D$ (nM)")
axes[0].annotate("masking ceiling", (42, 6e4), fontsize=8, color=C["blue"])
axes[0].annotate("release floor", (42, floor * 1.9), fontsize=8, color=C["ink2"])
axes[0].annotate("usable mask affinities", (14, 2e2), fontsize=8, color=C["ink2"])
fig.suptitle("Design window for mask affinity, from the measured 9GCN geometry", fontsize=10.5, y=1.02)
fig.savefig(FIGS / "fig4_window.png"); plt.close(fig)

# ---------- Figure 5: the two opposing requirements ----------
m = teth["model"]
km = np.logspace(-11, -2, 400)
fig, ax = plt.subplots(figsize=(5.8, 3.6))
for nres, col, ls in ((5, C["blue"], "-"), (20, C["aqua"], "-")):
    r = next(x for x in scan if x["terminus"] == "N" and x["linker_residues"] == nres)
    ce = r["c_eff_M"]
    ax.plot(km * 1e9, ce / (ce + km), lw=2.0, color=col, ls=ls)
    ax.annotate(f"masked, {nres}-residue linker", (km[180] * 1e9, ce / (ce + km[180])),
                textcoords="offset points", xytext=(4, 6), fontsize=7.5, color=col)
pnt = (m["toxin_conc_M"] / m["KT_M"]) / (1 + m["toxin_conc_M"] / m["KT_M"] + m["mask_conc_M"] / km)
ax.plot(km * 1e9, pnt, lw=2.0, color=C["orange"])
ax.annotate("toxin-bound after cleavage", (km[250] * 1e9, pnt[250]),
            textcoords="offset points", xytext=(-10, -16), fontsize=7.5, color=C["orange"])
ax.axhline(1 - m["leak_max"], lw=1.1, color=C["ink2"], ls=(0, (4, 2)))
ax.axhline(m["recovery_min"], lw=1.1, color=C["ink2"], ls=(0, (1, 2)))
ax.text(1.2e-2, 1 - m["leak_max"], " 99% masked", fontsize=7, color=C["ink2"], va="center")
ax.text(1.2e-2, m["recovery_min"], " 90% recovered", fontsize=7, color=C["ink2"], va="center")
ax.set_xscale("log"); ax.set_xlabel("mask–VHH $K_D$ (nM)")
ax.set_ylabel("fraction of nanobody")
ax.set_ylim(0, 1.05)
ax.set_title("The mask must be tight enough to block, loose enough to leave")
fig.savefig(FIGS / "fig5_tradeoff.png"); plt.close(fig)

# ---------- Figure 6: sensitivity of the release floor ----------
fig, ax = plt.subplots(figsize=(5.4, 3.3))
cols = {1e-10: C["blue"], 1e-9: C["orange"], 1e-8: C["aqua"]}
for kt in (1e-10, 1e-9, 1e-8):
    pts = [s for s in teth["sensitivity"] if s["KT_M"] == kt]
    xs = [p["T_M"] * 1e9 for p in pts]; ys = [p["km_min_for_release_M"] * 1e9 for p in pts]
    ax.plot(xs, ys, lw=2.0, marker="o", ms=5, color=cols[kt])
    ax.annotate(f"$K_{{D,T}}$ = {kt*1e9:g} nM", (xs[-1], ys[-1]),
                textcoords="offset points", xytext=(-6, 8), fontsize=7.5,
                color=cols[kt], ha="right")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("free toxin concentration (nM)")
ax.set_ylabel("minimum mask $K_D$ for 90% recovery (nM)")
ax.set_title("Release floor depends on how tightly the VHH binds toxin")
fig.savefig(FIGS / "fig6_sensitivity.png"); plt.close(fig)
print("figures written:", sorted(p.name for p in FIGS.glob("*.png")))
