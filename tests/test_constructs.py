import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nbswitch import constructs as C


def test_linker_length_is_exact():
    for n in (0, 1, 5, 7, 20):
        assert len(C.linker(n)) == n


def test_segments_tolerate_a_single_gap():
    flags = [True, True, False, True, True]
    assert C.contiguous_segments(flags, min_length=3) == [(0, 4)]


def test_segments_split_on_a_wider_gap():
    flags = [True, True, True, False, False, True, True, True]
    assert C.contiguous_segments(flags, min_length=3) == [(0, 2), (5, 7)]


def test_ranges_index_the_sequence_correctly():
    """Linkers are explicit sequences now, since a generated substrate is not
    described by its length."""
    masks = {"m": {"seq": "RGKRV", "kind": "graft", "note": ""}}
    linkers = [{"id": "cn01", "sequence": "GGSHMGSNLTGMYGGS",
                "substrate": "HMGSNLTGMY", "kind": "generated",
                "selectivity": 1.5}]
    lib = C.build_library("QVQLX" * 4, masks, linkers)
    item = next(c for c in lib if c["id"] == "m_cn01")
    seq = item["sequence"]
    a, b = item["mask_range"]
    assert seq[a - 1:b] == "RGKRV"
    a, b = item["linker_range"]
    assert seq[a - 1:b] == "GGSHMGSNLTGMYGGS"
    a, b = item["vhh_range"]
    assert seq[a - 1:b] == "QVQLX" * 4
    assert item["linker_substrate"] == "HMGSNLTGMY"


def test_the_reference_carries_no_linker():
    masks = {"m": {"seq": "RGKRV", "kind": "graft", "note": ""}}
    lib = C.build_library("QVQLX" * 4, masks, [])
    assert [c["id"] for c in lib] == ["VHH_only"]
    assert lib[0]["linker_range"] is None


def test_scramble_is_a_permutation():
    segments = [{"seq": "RGKRV", "length": 5, "n_epitope": 5,
                 "start_label": "R33", "end_label": "V37"}]
    masks = C.mask_candidates(segments)
    assert sorted(masks["epi1_scram"]["seq"]) == sorted("RGKRV")
    assert masks["polyG"]["kind"] == "inert control"


def test_generate_clears_the_output_directory_first():
    """CleaveNet appends to its output file.

    A work directory that survives between runs would accumulate every draw
    ever made in it, so the generated count would depend on how many times
    the stage had been run rather than on what was asked for.
    """
    import inspect
    from nbswitch import linkers as L
    source = inspect.getsource(L.generate)
    assert "rmtree" in source, "the output directory must be cleared before generating"
