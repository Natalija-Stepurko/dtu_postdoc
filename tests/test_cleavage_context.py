"""The in-context checks are necessary conditions, so they must behave like one."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pytest

from nbswitch import cleavage_context as CC


def test_clearance_is_open_in_free_space():
    """A segment with nothing around it must be fully approachable."""
    import biotite.structure as struc
    n = 7
    atoms = struc.AtomArray(n)
    atoms.coord = np.column_stack([np.arange(n) * 3.3, np.zeros(n), np.zeros(n)])
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)
    out = CC.clearance(atoms, scissile_resid=4)
    assert out["open_fraction"] == 1.0


def test_extension_is_one_for_a_fully_extended_chain():
    import biotite.structure as struc
    n = 6
    atoms = struc.AtomArray(n)
    atoms.coord = np.column_stack([np.arange(n) * CC.EXTENDED_RISE_A,
                                   np.zeros(n), np.zeros(n)])
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)
    assert CC.extension(atoms, (1, n)) == pytest.approx(1.0, abs=1e-6)


def test_extension_is_near_zero_for_a_closed_loop():
    """A segment returning to its start cannot lie in a catalytic groove."""
    import biotite.structure as struc
    n = 8
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    atoms = struc.AtomArray(n)
    atoms.coord = np.column_stack([np.cos(angle) * 5, np.sin(angle) * 5,
                                   np.zeros(n)])
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)
    assert CC.extension(atoms, (1, n)) < 0.2


def test_the_score_is_a_product_so_one_failure_is_fatal():
    """A requirement that can be compensated for is not a requirement."""
    import biotite.structure as struc
    n = 8
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    atoms = struc.AtomArray(n)
    atoms.coord = np.column_stack([np.cos(angle) * 5, np.sin(angle) * 5,
                                   np.zeros(n)])
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)
    out = CC.cleavage_accessibility(atoms, (1, n))
    assert out["context_score"] < 0.2
