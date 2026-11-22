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
    # A grafted segment taken from a disulfide-rich antigen can carry free
    # cysteines, which would scramble the nanobody's own disulfides. The
    # substituted variant keeps the shape of the graft without that liability.
    for name in [k for k in list(masks) if "C" in masks[k]["seq"]
                 and masks[k]["kind"] == "epitope graft"]:
        masks[f"{name}_noCys"] = {
            "seq": masks[name]["seq"].replace("C", "S"),
            "kind": "epitope graft",
            "note": f"{name} with cysteines substituted to serine"}

    masks["polyG"] = {"seq": "GGGGGGGG", "kind": "inert control",
                      "note": "no designed interaction; masking should fail"}
    return masks


def build_library(vhh_sequence, masks, linkers) -> list[dict]:
    """Every mask x linker combination, plus the free nanobody as reference.

    ``linkers`` are explicit sequences, not lengths: they carry a generated
    protease substrate, so a length alone no longer describes one.

    Ranges are 1-based inclusive positions in the construct sequence.
    """
    library = []
    for name, mask in masks.items():
        for link in linkers:
            seq = link["sequence"]
            n = len(seq)
            mask_end = len(mask["seq"])
            library.append({
                "id": f"{name}_{link['id']}", "mask": name,
                "mask_seq": mask["seq"], "mask_kind": mask["kind"],
                "linker_id": link["id"], "linker_residues": n, "linker_seq": seq,
                "linker_substrate": link.get("substrate", ""),
                "linker_kind": link.get("kind", "generated"),
                "linker_selectivity": link.get("selectivity"),
                "sequence": mask["seq"] + seq + vhh_sequence,
                "length": mask_end + n + len(vhh_sequence),
                "mask_range": [1, mask_end],
                "linker_range": [mask_end + 1, mask_end + n],
                "vhh_range": [mask_end + n + 1, mask_end + n + len(vhh_sequence)]})
    library.append({"id": "VHH_only", "mask": None, "mask_seq": "",
                    "mask_kind": "reference", "linker_id": None,
                    "linker_residues": 0, "linker_seq": "",
                    "linker_substrate": "", "linker_kind": "reference",
                    "linker_selectivity": None,
                    "sequence": vhh_sequence, "length": len(vhh_sequence),
                    "mask_range": None, "linker_range": None,
                    "vhh_range": [1, len(vhh_sequence)]})
    return library


def complex_cases(constructs, n_linkers: int = 2):
    """Which constructs to spend a complex prediction on.

    Derived from the library rather than named, because construct ids change
    whenever the linkers are regenerated: a hardcoded list silently selects
    nothing, or worse, reports results for molecules no longer in the library.

    The comparison needs a grafted mask and an inert one on the same linker, so
    that a difference between them is attributable to the mask, plus the same
    graft on a second linker to show whether the linker matters.
    """
    masked = [c for c in constructs if c.get("mask") and c.get("linker_id")]
    if not masked:
        return []

    def selectivity(c):
        s = c.get("linker_selectivity")
        return -1e9 if s is None else s

    linkers, seen = [], set()
    for c in sorted(masked, key=selectivity, reverse=True):
        if c["linker_id"] not in seen:
            seen.add(c["linker_id"])
            linkers.append(c["linker_id"])
    linkers = linkers[:max(1, n_linkers)]

    graft = next((c["mask"] for c in masked if c.get("mask_kind") == "epitope graft"), None)
    inert = next((c["mask"] for c in masked if c.get("mask_kind") == "inert control"), None)

    cases = []
    for mask, label in ((graft, "grafted mask"), (inert, "inert mask, should not block")):
        if mask is None:
            continue
        for i, lid in enumerate(linkers):
            if mask is inert and i > 0:
                continue        # the inert control only needs the one linker
            hit = next((c for c in masked
                        if c["mask"] == mask and c["linker_id"] == lid), None)
            if hit:
                cases.append((hit["id"], f"{label}, linker {lid}"))
    return cases
