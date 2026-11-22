"""The mask design space: what it generates and what it must not."""
import pytest

from nbswitch import mask_design as M
from nbswitch import constructs as C


ANTIGEN = "LECHNQQSSQPPTTKTCPGETNCYKKVWRDHRGTIIERGCGCPTVKPGIKLNCCTTDKCNN"
#: 1-based residue ids as the deposited chain numbers them
FIRST_ID = 1


def _window(start, end):
    """The antigen over an inclusive residue range, as stage 03a supplies it."""
    lo = max(start, FIRST_ID)
    hi = min(end, FIRST_ID + len(ANTIGEN) - 1)
    if hi < lo:
        return "", []
    ids = list(range(lo, hi + 1))
    return ANTIGEN[lo - FIRST_ID:hi - FIRST_ID + 1], ids


FACE = {"start": 42, "end": 53, "buried_A2": 580.7,
        "per_residue_buried_A2": [71.2, 24.7, 133.9, 47.0, 43.3, 51.3,
                                  8.0, 6.0, 5.0, 62.9, 63.4, 43.1]}
SECOND = {"start": 22, "end": 22, "buried_A2": 42.7,
          "per_residue_buried_A2": [42.7]}


def test_the_face_itself_is_the_first_candidate():
    masks = M.window_variants(_window, FACE)
    assert "epi_core" in masks
    assert masks["epi_core"]["seq"] == _window(42, 53)[0]
    assert masks["epi_core"]["extension"] == 0


def test_window_variants_carry_flanking_antigen():
    masks = M.window_variants(_window, FACE, extensions=(0, 2))
    assert len(masks["epi_core_ext2"]["seq"]) == len(masks["epi_core"]["seq"]) + 4


def test_a_single_residue_face_is_padded_into_a_peptide():
    """A hotspot can be one residue, which cannot fold into anything."""
    masks = M.secondary_face(_window, SECOND, index=2, min_length=8)
    assert len(masks["epi_face2"]["seq"]) >= 8


def test_both_faces_are_joined_in_antigen_order():
    masks = M.window_variants(_window, FACE)
    masks |= M.secondary_face(_window, SECOND, index=2)
    both = M.combined_faces(masks)["epi_both"]
    # the second face starts at residue 22, before the core at 42
    assert both["seq"].startswith(masks["epi_face2"]["seq"])
    assert both["seq"].endswith(masks["epi_core"]["seq"])
    assert M.FACE_SPACER in both["seq"]


def test_substitutions_hold_the_hotspots():
    face = _window(42, 53)[0]
    hotspots = C.hotspot_indices(FACE)
    masks = M.hotspot_substitutions(face, hotspots, n=3)
    assert masks, "no substitution variants were generated"
    for m in masks.values():
        assert len(m["seq"]) == len(face)
        for i in hotspots:
            assert m["seq"][i] == face[i], "a hotspot was substituted"


def test_substitutions_hold_cysteines_so_the_staple_is_not_a_second_variable():
    face = _window(42, 53)[0]
    masks = M.hotspot_substitutions(face, C.hotspot_indices(FACE), n=3)
    for m in masks.values():
        assert m["seq"].count("C") == face.count("C")


def test_substitution_variants_are_distinct():
    face = _window(42, 53)[0]
    masks = M.hotspot_substitutions(face, C.hotspot_indices(FACE), n=3)
    sequences = [m["seq"] for m in masks.values()]
    assert len(set(sequences)) == len(sequences)
    assert face not in sequences


def test_the_staple_variant_changes_only_the_cysteines():
    face = _window(42, 53)[0]
    variant = M.staple_variants(face)["epi_core_noCys"]["seq"]
    assert variant.count("C") == 0
    for a, b in zip(face, variant):
        assert a == b or (a == "C" and b == "S")


def test_controls_are_matched_to_the_face():
    face = _window(42, 53)[0]
    ctrl = M.controls(face)
    assert sorted(ctrl["epi_core_scram"]["seq"]) == sorted(face)
    assert len(ctrl["polyG"]["seq"]) == len(face)
    assert set(ctrl["polyG"]["seq"]) == {"G"}
    assert len(ctrl["epi_core_half"]["seq"]) < len(face)


def test_hotspots_are_the_residues_burying_most_surface():
    hotspots = C.hotspot_indices(FACE, fraction=0.6)
    per = FACE["per_residue_buried_A2"]
    assert 2 in hotspots, "the highest-burying position must be held"
    assert sum(per[i] for i in hotspots) >= 0.6 * sum(per)
    # the positions burying almost nothing are free to vary
    assert 8 not in hotspots


def test_design_masks_produces_one_sequence_per_candidate():
    masks = M.design_masks(_window, [FACE, SECOND],
                           hotspot_positions=C.hotspot_indices(FACE))
    assert len(masks) >= 10
    sequences = [m["seq"] for m in masks.values()]
    assert len(set(sequences)) == len(sequences), "duplicate masks in the space"
    kinds = {m["kind"] for m in masks.values()}
    assert {"epitope graft", "scrambled control", "inert control"} <= kinds
