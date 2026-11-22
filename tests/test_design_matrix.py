"""The design matrix: which mask and linker pairs get built, and how many.

Folding is the cost of this pipeline, so the size of the library is a design
decision. These pin the shape of that decision rather than its exact number.
"""
import pytest

from nbswitch import constructs as C


REF = 24


def _masks(n):
    return {f"m{i}": {"seq": "X" * 12, "kind": "epitope graft"} for i in range(n)}


def _linkers(n_reference, ladder_bases=(), ladder_lengths=()):
    out = [{"id": f"cn{i:02d}", "sequence": "X" * REF} for i in range(1, n_reference + 1)]
    for base in ladder_bases:
        for length in ladder_lengths:
            out.append({"id": f"{base}_n{length}", "sequence": "X" * length})
    return out


def test_without_a_ladder_it_is_a_full_cross_product():
    pairs = C.design_matrix(_masks(3), _linkers(4), REF)
    assert len(pairs) == 12


def test_length_variants_only_reach_the_chosen_masks():
    masks = _masks(5)
    linkers = _linkers(4, ladder_bases=("cn01",), ladder_lengths=(16, 32))
    pairs = C.design_matrix(masks, linkers, REF,
                            length_masks=("m0", "m1"), length_linkers=("cn01",))
    off_reference = [(m, l) for m, l in pairs if len(l["sequence"]) != REF]
    assert len(off_reference) == 2 * 2          # two masks, two extra lengths
    assert {m for m, _ in off_reference} == {"m0", "m1"}


def test_a_ladder_without_chosen_masks_adds_nothing():
    linkers = _linkers(4, ladder_bases=("cn01",), ladder_lengths=(16, 32))
    pairs = C.design_matrix(_masks(5), linkers, REF)
    assert all(len(l["sequence"]) == REF for _, l in pairs)


def test_the_matrix_is_far_smaller_than_the_cross_product():
    """The point of the matrix: vary one thing at a time, not everything."""
    masks = _masks(12)
    linkers = _linkers(6, ladder_bases=("cn01", "ctrl_gs", "ctrl_scram"),
                       ladder_lengths=(12, 16, 20, 28, 32))
    pairs = C.design_matrix(masks, linkers, REF,
                            length_masks=("m0", "m1", "m2"),
                            length_linkers=("cn01", "ctrl_gs", "ctrl_scram"))
    full = len(masks) * len(linkers)
    assert len(pairs) < full / 2
    assert len(pairs) + 1 <= C.MAX_CONSTRUCTS


def test_the_summary_counts_what_will_be_built():
    masks = _masks(12)
    linkers = _linkers(6, ladder_bases=("cn01", "ctrl_gs", "ctrl_scram"),
                       ladder_lengths=(12, 16, 20, 28, 32))
    plan = C.matrix_summary(masks, linkers, REF,
                            length_masks=("m0", "m1", "m2"),
                            length_linkers=("cn01", "ctrl_gs", "ctrl_scram"))
    pairs = C.design_matrix(masks, linkers, REF,
                            length_masks=("m0", "m1", "m2"),
                            length_linkers=("cn01", "ctrl_gs", "ctrl_scram"))
    assert plan["n_pairs"] == len(pairs)
    assert plan["n_constructs"] == len(pairs) + 1
    assert plan["n_at_reference_length"] + plan["n_length_variants"] == plan["n_pairs"]
    assert plan["estimated_fold_minutes"] > 0


def test_build_library_honours_the_matrix():
    masks = _masks(3)
    linkers = _linkers(2)
    pairs = C.design_matrix(masks, linkers, REF)[:2]
    library = C.build_library("QVQL" * 5, masks, linkers, pairs)
    assert len(library) == 3                    # two pairs plus the reference
    assert library[-1]["id"] == "VHH_only"
    for item in library[:-1]:
        a, b = item["vhh_range"]
        assert item["sequence"][a - 1:b] == "QVQL" * 5
