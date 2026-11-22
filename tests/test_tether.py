"""The tether model carries the headline result, so it is checked against
limits that can be reasoned about independently of the implementation."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pytest

from nbswitch import tether as T


def test_end_to_end_grows_with_length():
    lengths = np.array([5, 10, 20, 40])
    r2 = T.mean_square_end_to_end(lengths)
    assert np.all(np.diff(r2) > 0)


def test_long_chain_approaches_gaussian_limit():
    """For Lc >> lp the worm-like chain tends to <r^2> = 2 lp Lc."""
    n = 2000
    contour = n * T.RISE_PER_RESIDUE_A if hasattr(T, "RISE_PER_RESIDUE_A") else None
    from nbswitch.config import PERSISTENCE_LENGTH_A, RISE_PER_RESIDUE_A
    expected = 2 * PERSISTENCE_LENGTH_A * n * RISE_PER_RESIDUE_A
    assert T.mean_square_end_to_end(n) == pytest.approx(expected, rel=0.01)


def test_effective_concentration_falls_with_span():
    c = [T.effective_concentration(20, d) for d in (10, 20, 40)]
    assert c[0] > c[1] > c[2]


def test_effective_concentration_is_physically_plausible():
    """Tethered systems sit in the micromolar-to-millimolar range, not molar."""
    c = T.effective_concentration(10, 12.0)
    assert 1e-6 < c < 1.0


def test_masking_and_release_are_complementary():
    """Each bound is exactly the affinity that just satisfies its own requirement.

    Read the budgets from the config rather than pinning the numbers here: they
    are design choices and get changed, and a test that hardcodes them fails on
    the change instead of checking the relation that has to hold.
    """
    from nbswitch.config import MAX_SYSTEMIC_LEAK, MIN_TOXIN_RECOVERY

    c_eff = T.effective_concentration(5, 11.78)
    ceiling = T.masking_ceiling(c_eff)
    assert T.fraction_masked(c_eff, ceiling) == pytest.approx(1 - MAX_SYSTEMIC_LEAK,
                                                              abs=1e-6)
    floor = T.release_floor()
    assert T.fraction_toxin_bound(floor) == pytest.approx(MIN_TOXIN_RECOVERY, abs=1e-6)


def test_window_closes_for_a_long_span_and_short_linker():
    """A linker too short to reach cannot mask, whatever its affinity."""
    rows = T.scan_linkers(42.0, [2, 3, 4])
    assert not any(r["window_open"] for r in rows)


def test_window_opens_for_a_short_span():
    rows = T.scan_linkers(11.78, [5, 10, 20])
    assert all(r["window_open"] for r in rows)
    assert all(r["window_decades"] > 2 for r in rows)


def test_release_floor_rises_with_weaker_toxin_binding():
    tight = T.release_floor(toxin_kd=1e-10)
    loose = T.release_floor(toxin_kd=1e-8)
    assert loose > tight
