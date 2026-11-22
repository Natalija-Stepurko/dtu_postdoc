"""Paths, physical constants and the validated plotting palette."""
from __future__ import annotations
import json, os
from pathlib import Path

ROOT = Path(os.environ.get("NBSWITCH_ROOT", Path(__file__).resolve().parents[2]))
DATA = ROOT / "data"
DOCS = ROOT / "docs"

#: Complexes the pipeline can be pointed at. Each is an experimental VHH-toxin
#: structure from the same laboratory, which makes them comparable, and they
#: differ in the one property that decides whether masking is even testable:
#: how far the N-terminus sits from the binding site.
TARGETS = {
    "9RIU": {
        "pdb_id": "9RIU",
        "vhh_chains": ("A", "C"), "toxin_chains": ("B", "D"),
        "reference_vhh_chain": "A",
        "label": "broadly neutralizing VHH with short neurotoxin 1",
        "note": "one paratope, consistent between copies, N-terminus clear of it",
    },
    "9S9E": {
        "pdb_id": "9S9E",
        "vhh_chains": ("A", "B"), "toxin_chains": ("C", "D"),
        "reference_vhh_chain": "A",
        "label": "VHH TPL0870_0 with cytotoxin 1",
        "note": "biparatopic: the two copies share one paratope residue, so there is "
                "no single site to mask - rejected as a demonstration system",
    },
    "9GCN": {
        "pdb_id": "9GCN",
        "vhh_chains": ("A", "B", "C", "D"), "toxin_chains": ("V", "X", "Y", "Z"),
        "reference_vhh_chain": "A",
        "label": "VHH TPL1158_01_C09 with alpha-cobratoxin",
        "note": "N-terminal residue is itself part of the paratope",
    },
}

#: Which complex this run uses. Override with NBSWITCH_TARGET.
TARGET_ID = os.environ.get("NBSWITCH_TARGET", "9RIU")
if TARGET_ID not in TARGETS:
    raise SystemExit(f"unknown target {TARGET_ID}; choose from {sorted(TARGETS)}")
TARGET = TARGETS[TARGET_ID]

PDB_ID = TARGET["pdb_id"]
VHH_CHAINS = TARGET["vhh_chains"]
TOXIN_CHAINS = TARGET["toxin_chains"]
REFERENCE_VHH_CHAIN = TARGET["reference_vhh_chain"]

#: Worm-like-chain parameters for a Gly/Ser linker (Zhou 2001, J Phys Chem B).
PERSISTENCE_LENGTH_A = 4.5
RISE_PER_RESIDUE_A = 3.8
#: molecules per cubic angstrom -> mol/L
CONC_CONVERSION = 1660.54

#: Design requirements used to bound the mask affinity. These are choices, not
#: measurements, and the window moves with them: see the sensitivity note in
#: Module 01. The leak sets the ceiling and the recovery sets the floor, so
#: each bound answers to exactly one of them.
MAX_SYSTEMIC_LEAK = 0.05     # at most 5% of construct unmasked in circulation
MIN_TOXIN_RECOVERY = 0.90    # at least 90% of cleaved nanobody toxin-bound
MASK_CONC_M = 1e-6           # released mask concentration = construct concentration
TOXIN_CONC_M = 1e-7

#: Nanobody-toxin dissociation constant. VHH5 a-sNTx, the nanobody of the
#: worked example, is reported as low nanomolar against its short-neurotoxin
#: targets (Ahmadi et al. 2025, Fig. 2c and Supplementary Table 5), which is
#: where this value comes from. The entry's own complex is with P01426; the
#: paper does not label which of its sNTx targets that is, so this is the
#: reported range for the nanobody rather than the measured value for this pair.
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


#: Per-target output directories, so two complexes never overwrite each other.
RESULTS = ROOT / "results" / TARGET_ID
FIGURES = RESULTS / "figures"
STRUCTURES = RESULTS / "structures"
for _d in (DATA, DOCS, RESULTS, FIGURES, STRUCTURES):
    _d.mkdir(parents=True, exist_ok=True)


def results_dir(target_id: str) -> Path:
    """Where another target's outputs live, for cross-target comparison."""
    return ROOT / "results" / target_id


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

#: Heavy-atom distance defining a chemical contact: hydrogen bonds and salt
#: bridges fall at or below this. Van der Waals packing extends further, which
#: is why the masking footprint uses buried surface rather than this cutoff.
CONTACT_CUTOFF_A = 3.0

#: Antibody numbering scheme for CDR assignment: "kabat", "chothia" or "imgt".
#: The schemes disagree: on 9RIU, Kabat CDR2 is 17 residues and IMGT CDR2 is 7,
#: so anything quoting a CDR count has to name the scheme alongside it.
CDR_SCHEME = "imgt"
