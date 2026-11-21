"""Building the construct library.

A mask has to engage this particular paratope. One peptide is already known to
do so: the antigen's own contact segment. Grafting it needs no generative model
and brings its own negative controls, which is what makes it a usable starting
point for a demonstration.
"""
from __future__ import annotations
import random

from .config import RANDOM_SEED

LINKER_MOTIF = "GGGGS"
DEFAULT_LINKER_LENGTHS = (5, 10, 15, 20, 30)


def linker(n_residues: int, motif: str = LINKER_MOTIF) -> str:
    return (motif * (n_residues // len(motif) + 1))[:n_residues]


def contiguous_segments(flags, min_length=3, max_gap=1):
    """Index spans where ``flags`` is True, tolerating short interruptions."""
    spans, current, gap = [], [], 0
    for i, flag in enumerate(flags):
        if flag:
            current.append(i); gap = 0
        elif current:
            gap += 1
            if gap > max_gap:
                spans.append((current[0], current[-1])); current, gap = [], 0
    if current:
        spans.append((current[0], current[-1]))
    return [(a, b) for a, b in spans if b - a + 1 >= min_length]


def epitope_segments(toxin_labels, toxin_sequence, epitope_labels, **kw):
    """Contiguous stretches of antigen that touch the nanobody, longest contact first."""
    in_epitope = [label in set(epitope_labels) for label in toxin_labels]
    segments = []
    for a, b in contiguous_segments(in_epitope, **kw):
        segments.append({"start_label": toxin_labels[a], "end_label": toxin_labels[b],
                         "seq": toxin_sequence[a:b + 1], "length": b - a + 1,
                         "n_epitope": sum(in_epitope[a:b + 1])})
    segments.sort(key=lambda s: -s["n_epitope"])
    return segments


def mask_candidates(segments, seed: int = RANDOM_SEED) -> dict:
    """Grafted epitope masks plus the controls that make them interpretable."""
    rng = random.Random(seed)
    primary = segments[0]
    shuffled = list(primary["seq"]); rng.shuffle(shuffled)
    masks = {
        "epi1": {"seq": primary["seq"], "kind": "epitope graft",
                 "note": f"antigen {primary['start_label']}-{primary['end_label']}, "
                         f"{primary['n_epitope']} contact residues"},
        "epi1_scram": {"seq": "".join(shuffled), "kind": "scrambled control",
                       "note": f"same composition, shuffled (seed {seed})"},
        "epi1_half": {"seq": primary["seq"][:max(3, primary["length"] // 2)],
                      "kind": "truncated control", "note": "leading half of the graft"},
    }
    if len(segments) > 1:
        second = segments[1]
        masks["epi2"] = {"seq": second["seq"], "kind": "epitope graft",
                         "note": f"antigen {second['start_label']}-{second['end_label']}, "
                                 f"{second['n_epitope']} contact residues"}
    masks["polyG"] = {"seq": "GGGGGGGG", "kind": "inert control",
                      "note": "no designed interaction; masking should fail"}
    return masks


def build_library(vhh_sequence, masks, lengths=DEFAULT_LINKER_LENGTHS) -> list[dict]:
    """Every mask x linker combination, plus the free nanobody as reference.

    Ranges are 1-based inclusive positions in the construct sequence.
    """
    library = []
    for name, mask in masks.items():
        for n in lengths:
            link = linker(n)
            mask_end = len(mask["seq"])
            library.append({
                "id": f"{name}_L{n}", "mask": name, "mask_seq": mask["seq"],
                "mask_kind": mask["kind"], "linker_residues": n, "linker_seq": link,
                "sequence": mask["seq"] + link + vhh_sequence,
                "length": mask_end + n + len(vhh_sequence),
                "mask_range": [1, mask_end],
                "linker_range": [mask_end + 1, mask_end + n],
                "vhh_range": [mask_end + n + 1, mask_end + n + len(vhh_sequence)]})
    library.append({"id": "VHH_only", "mask": None, "mask_seq": "",
                    "mask_kind": "reference", "linker_residues": 0, "linker_seq": "",
                    "sequence": vhh_sequence, "length": len(vhh_sequence),
                    "mask_range": None, "linker_range": None,
                    "vhh_range": [1, len(vhh_sequence)]})
    return library
