"""The score-against-occlusion figure plots what its correlation describes.

The earlier version zipped twelve parents against six colours, so it drew 48
of 96 points and reported a correlation over all 96. Nothing failed; the
picture and the number simply described different sets, and the six parents
dropped were the tail of the alphabet, which is where most of the controls
live.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import json
import re

import pytest

from nbswitch.config import RESULTS
from nbswitch.mask_design import (MASK_FAMILIES, MASK_FAMILY_COLOUR,
                                  MASK_FAMILY_SHORT, family_by_construct,
                                  family_of)


def _load(name):
    path = RESULTS / name
    if not path.exists():
        pytest.skip(f"{name} has not been produced for this target")
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def data():
    return (_load("07b_library_scores.json"), _load("04_folding.json"),
            _load("08_redesign_folding.json"), _load("03_constructs.json"))


def test_every_family_has_a_colour_and_a_short_name():
    questions = [q for q, _ in MASK_FAMILIES]
    assert set(questions) == set(MASK_FAMILY_COLOUR)
    assert set(questions) == set(MASK_FAMILY_SHORT)
    assert len(set(MASK_FAMILY_COLOUR.values())) == len(questions), "a colour is reused"


def test_a_truncated_control_is_not_filed_under_the_graft_it_extends():
    """epi_core_half starts with epi_core_ and is a different mask."""
    assert family_of("epi_core_half") != family_of("epi_core")
    assert family_of("epi_core_half") == "Do the designed masks beat controls?"


def test_a_redesign_takes_its_family_from_the_parent_it_came_from(data):
    _, _, _, lib = data
    mask_of = {c["id"]: c.get("mask") for c in lib["constructs"]}
    assert family_by_construct("mpnn_epi_core_half_cn01", "epi_core_half_cn01",
                               mask_of) == "Do the designed masks beat controls?"
    assert family_by_construct("mpnn_epi_core_cn01", "epi_core_cn01",
                               mask_of) == family_of("epi_core")


def test_the_figure_draws_every_construct_its_correlation_covers(data):
    from nbswitch import figures as F
    scores, fold, redesign_fold, lib = data
    # the figure function returns the path it wrote, not the markup
    svg = Path(F.redesign_score_vs_occlusion(scores, fold, redesign_fold, lib)).read_text()

    occ = {r["id"] for src in (fold, redesign_fold) for r in src["results"]}
    expected = sum(1 for cid in scores["scores"] if cid in occ)
    m = re.search(r"n = (\d+)", svg)
    assert m, "the figure does not report how many points it plots"
    assert int(m.group(1)) == expected, (
        f"the figure says n = {m.group(1)} but {expected} constructs have both "
        "a score and an occlusion")

    # every marker is one point; <use> elements reference the shared path
    drawn = len(re.findall(r"<use ", svg))
    assert drawn >= expected, f"{drawn} markers drawn for {expected} points"


def test_the_legend_names_every_family_present(data):
    from nbswitch import figures as F
    scores, fold, redesign_fold, lib = data
    svg = Path(F.redesign_score_vs_occlusion(scores, fold, redesign_fold, lib)).read_text()
    mask_of = {c["id"]: c.get("mask") for c in lib["constructs"]}
    mask_name, parent_of = {}, {}
    for src in (fold, redesign_fold):
        for r in src["results"]:
            mask_name[r["id"]] = r.get("mask")
            parent_of[r["id"]] = r.get("parent")
    present = {family_by_construct(mask_name.get(c), parent_of.get(c), mask_of)
               for c in scores["scores"] if c in mask_name}
    present.discard(None)
    assert present, "no construct resolved to a family"
    text = " ".join(re.findall(r">([^<>]+)</text>", svg))
    for fam in present:
        assert MASK_FAMILY_SHORT[fam] in text, f"{fam} missing from the legend"


def test_library_and_redesign_scores_are_on_one_scale(data):
    """Both sides are scored the same way, or the axis means two things."""
    scores, _, _, _ = data
    assert "score_only" in scores["method"]
    assert scores["n_orders"] > 1, "a single decoding order is a noisy draw"
    kinds = {v["source"] for v in scores["scores"].values()}
    assert kinds == {"library", "redesign"}
    for v in scores["scores"].values():
        assert v["n_orders"] == scores["n_orders"]
