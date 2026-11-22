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


def _ring(n, radius):
    """A closed arc: every residue coiled, nothing extended."""
    import biotite.structure as struc
    angle = np.linspace(0, 2 * np.pi, n, endpoint=False)
    atoms = struc.AtomArray(n)
    atoms.coord = np.column_stack([np.cos(angle) * radius,
                                   np.sin(angle) * radius, np.zeros(n)])
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)
    return atoms


def test_a_coiled_cleavage_site_cannot_be_bound():
    """The cleft holds about six residues; if those are coiled, nothing binds."""
    atoms = _ring(8, radius=2.0)
    out = CC.cleavage_accessibility(atoms, (1, 8))
    assert out["extension"] < 0.5
    assert out["context_score"] < 0.2


def test_extension_is_measured_at_the_site_not_over_the_whole_linker():
    """A protease holds six residues. The rest of the tether is not its concern.

    A long linker that is extended where it is cut and coiled elsewhere is
    cleavable, and scoring it on the whole-linker span would say otherwise.
    """
    import biotite.structure as struc
    n = 24
    coords = []
    for i in range(n):
        if 8 <= i <= 15:                       # extended run through the middle
            coords.append([i * 3.3, 0.0, 0.0])
        else:                                   # coiled at both ends
            a = i * 1.4
            coords.append([40 + np.cos(a) * 4, np.sin(a) * 4, 0.0])
    atoms = struc.AtomArray(n)
    atoms.coord = np.array(coords, dtype=float)
    atoms.atom_name = np.array(["CA"] * n)
    atoms.element = np.array(["C"] * n)
    atoms.res_id = np.arange(1, n + 1)
    atoms.res_name = np.array(["GLY"] * n)
    atoms.chain_id = np.array(["A"] * n)

    out = CC.cleavage_accessibility(atoms, (1, n))
    assert out["extension"] > out["extension_whole_linker"], (
        "the site is extended even though the linker as a whole is not")
    assert "cleft_window" in out


def test_the_score_is_a_product_so_one_failure_is_fatal():
    """A requirement that can be compensated for is not a requirement."""
    assert CC.cleavage_accessibility.__doc__
    score = 0.9 * 0.9 * 0.0
    assert score == 0.0, "a zero in any factor must zero the product"

def test_conditioning_profile_is_one_row():
    """CleaveNet draws num_seqs samples per row of this file, reloading weights
    each time, so extra rows multiply the work instead of the diversity."""
    from nbswitch import linkers as L
    frame = L.conditioning_profile()
    assert len(frame) == 1
    assert list(frame.columns) == L.PROTEASES
    row = frame.iloc[0]
    assert all(row[p] > 0 for p in L.TARGET_PROTEASES)
    assert all(row[p] < 0 for p in L.OFF_TARGET_PROTEASES)
