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
    masks = {"m": {"seq": "RGKRV", "kind": "graft", "note": ""}}
    lib = C.build_library("QVQLX" * 4, masks, lengths=(5,))
    item = next(c for c in lib if c["id"] == "m_L5")
    seq = item["sequence"]
    a, b = item["mask_range"]
    assert seq[a - 1:b] == "RGKRV"
    a, b = item["linker_range"]
    assert seq[a - 1:b] == C.linker(5)
    a, b = item["vhh_range"]
    assert seq[a - 1:b] == "QVQLX" * 4


def test_scramble_is_a_permutation():
    segments = [{"seq": "RGKRV", "length": 5, "n_epitope": 5,
                 "start_label": "R33", "end_label": "V37"}]
    masks = C.mask_candidates(segments)
    assert sorted(masks["epi1_scram"]["seq"]) == sorted("RGKRV")
    assert masks["polyG"]["kind"] == "inert control"
