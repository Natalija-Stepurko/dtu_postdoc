"""Generating masks from the antigen's contact face.

The mask has to occupy the paratope while it is tethered and let go once the
linker is cut. The antigen already solves the first half of that problem, so
its contacting face is where this starts: it is a sequence known to be
complementary to this surface, taken from a structure rather than proposed.

A face is a starting point and not a candidate, because three things about it
are unknown and each is a design variable:

**How much to carry.** The contact residues alone may not hold their
conformation; the surrounding antigen may hold them and may also add bulk,
cysteines and length that do no work.

**Which residues matter.** Some contacts bury most of the surface and others
graze. Holding the first fixed and varying the rest samples sequence space
without discarding what makes the face fit.

**Whether it needs stapling.** A loop cut out of a disulfide-rich antigen can
keep a native bond, form one the antigen does not have, or come apart.

Each is turned into a family of candidates here. Everything in this module is
sequence-level and runs in milliseconds: no model is consulted, so these are
hypotheses to be scored downstream and not predictions.

**Methods this module does not implement.** Inverse folding on the folded
construct is stage 07 (ProteinMPNN), which is the one generative method in the
pipeline. De novo binder design against the paratope (RFdiffusion, RFantibody,
BindCraft) would generate masks that owe nothing to the antigen, which is the
obvious next source of diversity; RFdiffusion is installed without weights and
wants a GPU. A masked language model (ESM-2) could resample non-hotspot
positions in place of the substitution rules used here.
"""
from __future__ import annotations

import random

from .config import RANDOM_SEED

#: Substitutions that preserve the character of a side chain. Used to vary
#: non-hotspot positions without changing what the face is made of.
CONSERVATIVE = {
    "A": "GS", "G": "AS", "S": "TAG", "T": "SA", "C": "SA",
    "V": "IL", "I": "LV", "L": "IVM", "M": "LI", "F": "YW",
    "Y": "FW", "W": "FY", "P": "AG", "D": "EN", "E": "DQ",
    "N": "DQS", "Q": "ENK", "K": "RQ", "R": "KQ", "H": "QNY",
}

#: How far around a contact face to look, in residues each side. Zero is the
#: face itself.
WINDOW_EXTENSIONS = (0, 2, 4)

#: A secondary face can be a single hotspot residue, which is not a peptide.
#: Pad it with flanking antigen until it is at least this long.
MIN_FACE_LENGTH = 8

#: Spacer joining two faces, in place of the antigen between them when those
#: residues make no contact.
FACE_SPACER = "GSG"


def window_variants(window, face, extensions=WINDOW_EXTENSIONS) -> dict:
    """The face itself, and the face carrying flanking antigen.

    ``window`` is a callable returning (sequence, residue ids) for an inclusive
    residue range, so this function does not need the structure.
    """
    out = {}
    for ext in extensions:
        seq, ids = window(face["start"] - ext, face["end"] + ext)
        if not seq:
            continue
        name = "epi_core" if ext == 0 else f"epi_core_ext{ext}"
        out[name] = {
            "seq": seq, "kind": "epitope graft", "method": "graft",
            "note": (f"contact face"
                     + ("" if ext == 0 else f" plus {ext} residues of antigen each side")
                     + f", {seq.count('C')} cysteines"),
            "extension": ext, "start": ids[0], "end": ids[-1]}
    return out


def secondary_face(window, face, index: int = 2,
                   min_length: int = MIN_FACE_LENGTH) -> dict:
    """A further contact face, padded until it is long enough to be a peptide."""
    pad = max(0, (min_length - (face["end"] - face["start"] + 1) + 1) // 2)
    seq, ids = window(face["start"] - pad, face["end"] + pad)
    if len(seq) < 3:
        return {}
    return {f"epi_face{index}": {
        "seq": seq, "kind": "epitope graft", "method": "graft",
        "note": f"second contact face, around a hotspot burying "
                f"{face['buried_A2']:.0f} A2",
        "extension": pad, "start": ids[0], "end": ids[-1]}}


def combined_faces(masks: dict, spacer: str = FACE_SPACER) -> dict:
    """Both faces in one mask, joined in the order the antigen has them.

    A mask of one face leaves the other uncovered, and the antigen presents
    both at once. The residues between them make no contact, so carrying them
    would add bulk for nothing; a short spacer stands in.
    """
    faces = [m for name, m in masks.items()
             if name == "epi_core" or name.startswith("epi_face")]
    if len(faces) < 2:
        return {}
    first, second = sorted(faces, key=lambda m: m.get("start", 0))[:2]
    return {"epi_both": {
        "seq": f"{first['seq']}{spacer}{second['seq']}",
        "kind": "epitope graft", "method": "graft",
        "note": "both contact faces, joined by a spacer in place of the antigen "
                "residues between them, which make no contact",
        "extension": 0}}


def hotspot_substitutions(face_seq: str, hotspot_positions, n: int = 3,
                          seed: int = RANDOM_SEED) -> dict:
    """Vary the residues that are not doing the binding, keep the ones that are.

    The positions burying the most surface are held fixed and the rest take
    conservative substitutions. This asks whether the face needs to be the
    antigen's exact sequence, or only needs its hotspots in the right places.

    ``hotspot_positions`` are 0-based indices into ``face_seq``.
    """
    rng = random.Random(seed)
    # cysteines are held as well as hotspots: substituting one removes the
    # disulfide, which is the variable staple_variants exists to test, and
    # changing two things at once makes neither attributable
    held = set(hotspot_positions) | {i for i, c in enumerate(face_seq) if c == "C"}
    free = [i for i in range(len(face_seq)) if i not in held]
    out, seen = {}, {face_seq}
    attempts = 0
    while len(out) < n and attempts < 50 * n:
        attempts += 1
        seq = list(face_seq)
        for i in free:
            options = CONSERVATIVE.get(seq[i], "")
            if options:
                seq[i] = rng.choice(options)
        candidate = "".join(seq)
        if candidate in seen:          # draws repeat; a duplicate is not a variant
            continue
        seen.add(candidate)
        k = len(out) + 1
        out[f"epi_sub{k:02d}"] = {
            "seq": candidate, "kind": "epitope graft", "method": "substitution",
            "note": f"contact face with its hotspots and cysteines held, and the "
                    f"other {len(free)} positions conservatively substituted",
            "extension": 0}
    return out


def staple_variants(face_seq: str) -> dict:
    """With and without the disulfide the face carries.

    A graft cut out of a disulfide-rich antigen may keep a native bond. Whether
    that bond helps is a question the structure can answer, so both are built:
    removing the cysteines removes the staple and nothing else.
    """
    out = {}
    if face_seq.count("C") >= 2:
        out["epi_core_noCys"] = {
            "seq": face_seq.replace("C", "S"),
            "kind": "epitope graft", "method": "staple variant",
            "note": "cysteines substituted to serine, which removes the "
                    "disulfide holding the face in shape and changes nothing else",
            "extension": 0}
    return out


def controls(face_seq: str, seed: int = RANDOM_SEED) -> dict:
    """What the designed masks have to beat.

    Each is matched to the face on everything except the property being
    tested, so a difference is attributable to that property.
    """
    rng = random.Random(seed)
    shuffled = list(face_seq)
    rng.shuffle(shuffled)
    return {
        "epi_core_scram": {
            "seq": "".join(shuffled), "kind": "scrambled control",
            "method": "control",
            "note": f"the face shuffled, seed {seed}: same composition, same "
                    f"length, no order",
            "extension": 0},
        "epi_core_half": {
            "seq": face_seq[:max(3, len(face_seq) // 2)],
            "kind": "truncated control", "method": "control",
            "note": "leading half of the face: right residues, not enough of them",
            "extension": 0},
        "polyG": {
            "seq": "G" * len(face_seq), "kind": "inert control",
            "method": "control",
            "note": f"no designed interaction, matched to the face at "
                    f"{len(face_seq)} residues",
            "extension": 0},
    }


def design_masks(window, segments, hotspot_positions=None,
                 n_substitutions: int = 3, seed: int = RANDOM_SEED) -> dict:
    """Every mask candidate, from the contact faces outward.

    ``window(start, end)`` returns (sequence, residue ids) over the antigen.
    ``segments`` are contact faces, most buried first.
    """
    if not segments:
        return {}
    masks = window_variants(window, segments[0])
    for i, face in enumerate(segments[1:3], 2):
        masks |= secondary_face(window, face, index=i)
    masks |= combined_faces(masks)

    core = masks["epi_core"]["seq"]
    masks |= staple_variants(core)
    if hotspot_positions:
        masks |= hotspot_substitutions(core, hotspot_positions,
                                       n=n_substitutions, seed=seed)
    masks |= controls(core, seed=seed)
    return masks
