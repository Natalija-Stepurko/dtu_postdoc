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


def contact_segments(atoms, vhh_chain: str, toxin_chain: str,
                     cutoff: float = None, max_gap: int = None,
                     min_run: int = None) -> list[dict]:
    """Stretches of antigen that actually touch the nanobody.

    A mask is meant to stand in for the antigen at the paratope, so the
    sequence to start from is the one the antigen presents there. That is read
    off the structure directly: every antigen residue with a heavy atom within
    ``cutoff`` of the nanobody, grouped into runs along the chain.

    Runs separated by a short gap are merged. The gap residues are not in
    contact, but they are part of the same loop, and a graft that skipped them
    would be a sequence the antigen does not contain and could not adopt.
    """
    from .structure import residues, contact_residues, buried_on_binding
    from .config import (MASK_CONTACT_CUTOFF_A, MASK_MAX_GAP, MASK_MIN_RUN)
    import numpy as np

    cutoff = MASK_CONTACT_CUTOFF_A if cutoff is None else cutoff
    max_gap = MASK_MAX_GAP if max_gap is None else max_gap
    min_run = MASK_MIN_RUN if min_run is None else min_run

    a = atoms[atoms.chain_id == vhh_chain]
    b = atoms[atoms.chain_id == toxin_chain]
    res_b = residues(b)
    _, dist = contact_residues(a, b, cutoff)
    _, delta = buried_on_binding(a, b)

    touching = sorted(int(res_b.res_id[i]) for i in range(len(res_b))
                      if dist[i] <= cutoff)
    if not touching:
        return []

    groups, current = [], [touching[0]]
    for rid in touching[1:]:
        if rid - current[-1] - 1 <= max_gap:
            current.append(rid)
        else:
            groups.append(current)
            current = [rid]
    groups.append(current)

    index = {int(res_b.res_id[i]): i for i in range(len(res_b))}
    segments = []
    for g in groups:
        span = list(range(g[0], g[-1] + 1))
        rows = [index[r] for r in span if r in index]
        if len(g) < min_run:
            continue
        segments.append({
            "start": g[0], "end": g[-1], "length": len(span),
            "seq": "".join(res_b.code[i] for i in rows),
            "start_label": res_b.labels[rows[0]],
            "end_label": res_b.labels[rows[-1]],
            "n_contact": len(g),
            # the older footprint-derived segments call this n_epitope, and
            # mask_candidates reads that name
            "n_epitope": len(g),
            "contact_labels": [res_b.labels[index[r]] for r in g if r in index],
            "buried_A2": round(float(sum(delta[i] for i in rows)), 1),
            "cutoff_A": cutoff})
    return sorted(segments, key=lambda s: -s["buried_A2"])


#: How far to extend the contact seed on each side, in residues. Zero is the
#: contact face itself; the wider windows carry flanking antigen that is not in
#: contact but may hold the face in shape.
MASK_EXTENSIONS = (0, 2, 4)


def mask_family(atoms, toxin_chain: str, segments, seed: int = RANDOM_SEED,
                extensions=MASK_EXTENSIONS) -> dict:
    """A space of masks seeded by the contact face, not a single candidate.

    The contact face is where to start, not what to test. How much antigen to
    carry with it is an open question: too little and the peptide has nothing
    holding its shape, too much and it brings residues that do no work, add
    length and, on this antigen, add cysteines whose partners were left behind.
    The window is therefore a variable, and each setting becomes a candidate.

    Controls are generated from the seed itself so that composition, length and
    attachment are matched and only the designed property differs.
    """
    from .structure import residues
    import numpy as np

    b = atoms[atoms.chain_id == toxin_chain]
    res = residues(b)
    by_id = {int(res.res_id[i]): res.code[i] for i in range(len(res))}
    label = {int(res.res_id[i]): res.labels[i] for i in range(len(res))}
    lo_bound, hi_bound = min(by_id), max(by_id)

    def window(start, end):
        start, end = max(start, lo_bound), min(end, hi_bound)
        ids = [r for r in range(start, end + 1) if r in by_id]
        return "".join(by_id[r] for r in ids), ids

    primary = segments[0]
    masks: dict[str, dict] = {}

    for ext in extensions:
        seq, ids = window(primary["start"] - ext, primary["end"] + ext)
        if not seq:
            continue
        name = "epi_core" if ext == 0 else f"epi_core_ext{ext}"
        cys = seq.count("C")
        masks[name] = {
            "seq": seq, "kind": "epitope graft",
            "note": (f"antigen {label[ids[0]]}-{label[ids[-1]]}, contact face"
                     + ("" if ext == 0 else f" plus {ext} residues each side")
                     + f", {cys} cysteines"),
            "extension": ext, "start": ids[0], "end": ids[-1]}

    # a second contact face, where there is one, for whether both are needed.
    # The first segment is already the core, so it is not repeated.
    for i, seg in enumerate(segments[1:3], 2):
        masks[f"epi_face{i}"] = {
            "seq": seg["seq"], "kind": "epitope graft",
            "note": f"second contact face, {seg['start_label']}-{seg['end_label']}, "
                    f"burying {seg['buried_A2']:.0f} A2",
            "extension": 0, "start": seg["start"], "end": seg["end"]}

    core = masks["epi_core"]["seq"]
    if "C" in core:
        masks["epi_core_noCys"] = {
            "seq": core.replace("C", "S"), "kind": "epitope graft",
            "note": "contact face with cysteines substituted to serine, which "
                    "removes the disulfide that holds it in shape",
            "extension": 0}

    rng = random.Random(seed)
    shuffled = list(core)
    rng.shuffle(shuffled)
    masks["epi_core_scram"] = {
        "seq": "".join(shuffled), "kind": "scrambled control",
        "note": f"the contact face shuffled (seed {seed}): same composition, "
                f"same length, no order", "extension": 0}
    masks["epi_core_half"] = {
        "seq": core[:max(3, len(core) // 2)], "kind": "truncated control",
        "note": "leading half of the contact face", "extension": 0}
    masks["polyG"] = {
        "seq": "G" * len(core), "kind": "inert control",
        "note": f"no designed interaction, matched to the face at "
                f"{len(core)} residues", "extension": 0}
    return masks
