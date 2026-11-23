"""The mask-family figure: grouped by question, aligned on the antigen."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import json
import pytest

from nbswitch import report as R
from nbswitch.config import RESULTS


@pytest.fixture(scope="module")
def md():
    path = RESULTS / "03a_masks.json"
    if not path.exists():
        pytest.skip("mask design has not been run for this target")
    return json.loads(path.read_text())


def test_every_mask_appears_under_exactly_one_question():
    named = [n for _, names in R.MASK_FAMILIES for n in names]
    assert len(named) == len(set(named)), "a mask is filed under two questions"


def test_every_mask_in_the_run_is_in_the_figure(md):
    named = {n for _, names in R.MASK_FAMILIES for n in names}
    assert set(md["masks"]) == named, (
        "the figure and the run disagree about which masks exist: "
        f"{set(md['masks']) ^ named}")
    assert named <= set(R.MASK_LABELS), "a mask has no readable label"


def test_derived_masks_are_drawn_under_the_loop_they_came_from(md):
    """They carry no antigen coordinates, so without this they land at zero."""
    masks = md["masks"]
    core = masks["epi_core"]
    for name in ("epi_sub01", "epi_core_noCys", "epi_core_scram", "polyG"):
        blocks = R.mask_blocks(name, masks[name], masks)
        assert blocks == [(core["start"], masks[name]["seq"])]


def test_the_joined_mask_is_drawn_as_the_two_regions_it_joins(md):
    """Drawn as one run, its second half would sit under the wrong residues."""
    masks = md["masks"]
    blocks = R.mask_blocks("epi_both", masks["epi_both"], masks)
    assert len(blocks) == 3, "expected region, spacer, region"
    (p1, s1), (pspacer, spacer), (p2, s2) = blocks
    assert (p1, s1) == (masks["epi_face2"]["start"], masks["epi_face2"]["seq"])
    assert (p2, s2) == (masks["epi_core"]["start"], masks["epi_core"]["seq"])
    assert pspacer is None and s1 + spacer + s2 == masks["epi_both"]["seq"]


def test_the_extensions_report_the_cysteines_the_cut_orphans(md):
    """The argument against a longer graft, so it has to come from the data."""
    masks = md["masks"]
    assert masks["epi_core"]["cysteines_orphaned"] == []
    assert masks["epi_core"]["disulfides_kept"] == [[42, 53]]
    for name in ("epi_core_ext2", "epi_core_ext4"):
        m = masks[name]
        assert m["cysteines_orphaned"] == [40, 54]
        assert sorted(m["orphan_partners"].values()) == [17, 59]
        assert m["disulfides_kept"] == [[42, 53]], "the core pair must survive"


def test_the_figure_draws_a_box_for_every_residue_of_every_mask(md):
    import re
    svg = R.mask_families(md["masks"], {int(l[1:]) for l in md["hotspot_labels"]})
    drawn = len(re.findall(r'class="m-t"', svg))
    expected = sum(len(m["seq"]) for m in md["masks"].values())
    spacer = len(md["masks"]["epi_both"]["seq"]) - (
        len(md["masks"]["epi_face2"]["seq"]) + len(md["masks"]["epi_core"]["seq"]))
    assert drawn == expected - spacer, "the spacer is labelled once, not per residue"
