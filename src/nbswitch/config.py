"""Paths, physical constants and the validated plotting palette."""
from __future__ import annotations
import json, os
from pathlib import Path

ROOT = Path(os.environ.get("NBSWITCH_ROOT", Path(__file__).resolve().parents[2]))
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
STRUCTURES = RESULTS / "structures"
DOCS = ROOT / "docs"

for _d in (DATA, RESULTS, FIGURES, STRUCTURES, DOCS):
    _d.mkdir(parents=True, exist_ok=True)

#: Reference complex: VHH TPL1158_01_C09 bound to alpha-cobratoxin.
PDB_ID = "9GCN"
VHH_CHAINS = ("A", "B", "C", "D")
TOXIN_CHAINS = ("V", "X", "Y", "Z")
REFERENCE_VHH_CHAIN = "A"

#: Worm-like-chain parameters for a Gly/Ser linker (Zhou 2001, J Phys Chem B).
PERSISTENCE_LENGTH_A = 4.5
RISE_PER_RESIDUE_A = 3.8
#: molecules per cubic angstrom -> mol/L
CONC_CONVERSION = 1660.54

#: Design requirements used to bound the mask affinity.
MAX_SYSTEMIC_LEAK = 0.01     # at most 1% of construct unmasked in circulation
MIN_TOXIN_RECOVERY = 0.90    # at least 90% of cleaved nanobody toxin-bound
MASK_CONC_M = 1e-6           # released mask concentration = construct concentration
TOXIN_CONC_M = 1e-7
VHH_TOXIN_KD_M = 1e-9

#: Residue is counted in the interface if it buries more than this on binding.
BURIAL_THRESHOLD_A2 = 1.0
SASA_PROBE_A = 1.4
SASA_POINTS = 500

RANDOM_SEED = 42

REPO_URL = "https://github.com/Natalija-Stepurko/dtu_postdoc"
PAGES_URL = "https://natalija-stepurko.github.io/dtu_postdoc"
RAW_URL = "https://raw.githubusercontent.com/Natalija-Stepurko/dtu_postdoc/main"
NOTEBOOK_URL = ("https://colab.research.google.com/github/Natalija-Stepurko/"
                "dtu_postdoc/blob/main/notebooks/esmfold_masked_vhh.ipynb")

#: Categorical palette, slots 1-4. Checked with the data-viz validator against
#: the light surface #fcfcfb: all-pairs pass for slots 1-3, adjacent pass for
#: 1-4. Aqua and yellow fall below 3:1 contrast, so charts using them carry
#: visible direct labels.
PALETTE = {
    "blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100",
    "ink": "#0b0b0b", "ink2": "#52514e", "grid": "#d8d7d2",
    "surface": "#fcfcfb", "band": "#e8eef7",
}


def save_json(name: str, obj) -> Path:
    path = RESULTS / name
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str))
    return path


def load_json(name: str):
    return json.loads((RESULTS / name).read_text())


def plt_style():
    """Matplotlib with the house style applied. Returns pyplot."""
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
