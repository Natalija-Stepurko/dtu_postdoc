"""Shared constants and helpers. Every number quoted on the site comes from these scripts."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGS = RESULTS / "figures"
for d in (DATA, RESULTS, FIGS):
    d.mkdir(parents=True, exist_ok=True)

# Reference complex: VHH TPL1158_01_C09 bound to alpha-cobratoxin.
PDB_ID = "9GCN"

# Validated categorical palette, slots 1-4 (blue, orange, aqua, yellow).
# Checked with the data-viz validator: all-pairs pass for slots 1-3, adjacent
# pass for 1-4, light surface #fcfcfb. Aqua and yellow sit below 3:1 contrast,
# so every chart using them carries visible direct labels (the relief rule).
TOL = {
    "blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
    "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#d8d7d2",
    "surface": "#fcfcfb", "band": "#e8eef7",
}

def save_json(name, obj):
    p = RESULTS / name
    p.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str))
    print(f"wrote {p.relative_to(ROOT)}")
    return p

def load_json(name):
    return json.loads((RESULTS / name).read_text())

def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "axes.labelsize": 9, "axes.titlesize": 10, "legend.frameon": False,
        "figure.facecolor": "white", "savefig.bbox": "tight",
    })
    return plt
