"""Protonation arithmetic for the pH trigger."""
import pytest

from nbswitch import triggers as T


def test_half_protonated_at_the_pka():
    assert T.fraction_protonated(6.0, pka=6.0) == pytest.approx(0.5)


def test_a_unit_above_the_pka_is_about_a_tenth():
    assert T.fraction_protonated(7.0, pka=6.0) == pytest.approx(1 / 11, abs=1e-6)


def test_protonation_falls_as_ph_rises():
    values = [T.fraction_protonated(p) for p in (5.0, 6.0, 7.0, 8.0)]
    assert values == sorted(values, reverse=True)


def test_blood_and_bite_site_for_a_free_histidine():
    """The numbers the page quotes, pinned."""
    assert T.fraction_protonated(7.35) == pytest.approx(0.0428, abs=5e-4)
    assert T.fraction_protonated(6.50) == pytest.approx(0.2403, abs=5e-4)
    assert T.fraction_protonated(6.00) == pytest.approx(0.5000, abs=5e-4)


def test_the_comparison_is_conservative():
    """Worst case both ends: most protonated in blood, least at the site."""
    w = T.switch_window()
    # the reported values are rounded for display, so compare at that precision
    assert w["blood_worst"] == pytest.approx(T.fraction_protonated(7.35), abs=1e-4)
    assert w["site_worst"] == pytest.approx(T.fraction_protonated(6.50), abs=1e-4)
    assert w["fold_change"] > 1


def test_raising_the_pka_buys_site_protonation_and_costs_systemic_leak():
    """The trade-off the page rests on: it is not free to tune the pKa up."""
    scan = T.sensitivity()
    by_pka = {r["pka"]: r for r in scan}
    assert by_pka[7.0]["site_worst"] > by_pka[5.5]["site_worst"]
    assert by_pka[7.0]["blood_worst"] > by_pka[5.5]["blood_worst"]
    # and the fold change shrinks as the baseline rises
    assert by_pka[7.0]["fold_change"] < by_pka[5.5]["fold_change"]


def test_most_of_the_residue_is_still_neutral_at_the_site():
    w = T.switch_window()
    assert w["unprotonated_at_site"] == pytest.approx(1 - w["site_worst"], abs=1e-4)
    assert w["unprotonated_at_site"] > 0.5
