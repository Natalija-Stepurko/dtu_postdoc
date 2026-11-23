"""What the page publishes alongside itself.

Both of these were silent faults rather than failures: the page looked right
while carrying files nothing linked to, and every caption was announced twice
to anyone not reading it with their eyes.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from nbswitch import report as R


def test_a_caption_is_not_repeated_in_the_alt_text():
    """The figcaption sits beside the image, so alt must not restate it.

    With both filled in, a screen reader reads the sentence twice and a copied
    selection contains it twice. It also put caption markup inside an
    attribute, where a quote would have broken the tag.
    """
    html = R._chart("some_figure", '<b>A caption</b> with "quotes" in it.')
    assert 'alt=""' in html
    assert html.count("A caption") == 1


def test_copying_figures_removes_the_ones_no_longer_generated(tmp_path, monkeypatch):
    """docs/figures mirrors the results directory instead of accumulating.

    Five figures stayed published after the page stopped linking to them,
    because the copy step only ever added.
    """
    figures, docs = tmp_path / "figures", tmp_path / "docs"
    figures.mkdir()
    (docs / "figures").mkdir(parents=True)
    (docs / "structures").mkdir()
    (figures / "kept.svg").write_text("<svg/>")
    (docs / "figures" / "kept.svg").write_text("old")
    (docs / "figures" / "dropped.svg").write_text("stale")
    monkeypatch.setattr(R, "FIGURES", figures)
    monkeypatch.setattr(R, "DOCS", docs)
    monkeypatch.setattr(R, "STRUCTURES", tmp_path / "structures")

    R.copy_assets([])

    assert {p.name for p in (docs / "figures").glob("*.svg")} == {"kept.svg"}
    assert (docs / "figures" / "kept.svg").read_text() == "<svg/>"


def test_an_unbuilt_figure_stage_does_not_wipe_the_published_figures(tmp_path, monkeypatch):
    """A page build in a checkout that never ran the figure stage must not
    delete the figures already published, so the prune is conditional."""
    figures, docs = tmp_path / "figures", tmp_path / "docs"
    figures.mkdir()
    (docs / "figures").mkdir(parents=True)
    (docs / "structures").mkdir()
    (docs / "figures" / "published.svg").write_text("keep me")
    monkeypatch.setattr(R, "FIGURES", figures)
    monkeypatch.setattr(R, "DOCS", docs)
    monkeypatch.setattr(R, "STRUCTURES", tmp_path / "structures")

    R.copy_assets([])

    assert (docs / "figures" / "published.svg").exists()


def test_small_counts_spell_out():
    """The suite passed while the page build died on a missing _WORDS, so the
    helper that reads it is checked directly."""
    assert R.word(2) == "two"
    assert R.word(12) == "twelve"
    assert R.word(714) == "714", "past the table it falls back to digits"
    assert R.word(0) == "zero"
