"""Residue keying is the one place a silent error would corrupt everything
downstream, so it is pinned here.

These tests name 9GCN explicitly: the insertion-code trap is a property of that
deposit, not of whichever target the pipeline is currently pointed at.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from nbswitch import structure as S

KABAT_NUMBERED = "9GCN"


@pytest.fixture(scope="module")
def atoms():
    return S.read_atoms(S.fetch_structure(KABAT_NUMBERED))


def test_insertion_codes_are_not_collapsed(atoms):
    """9GCN is Kabat-numbered: 122 residues over 112 distinct res_id values.
    Keying on res_id alone merges residues and truncates CDR3."""
    chain = atoms[atoms.chain_id == "A"]
    res = S.residues(chain)
    assert len(res) == 122
    assert len(set(chain.res_id.tolist())) == 112


def test_cdr3_is_full_length(atoms):
    spans = S.find_cdrs(S.residues(atoms[atoms.chain_id == "A"]).sequence)
    start, stop = spans["CDR3"]
    assert stop - start == 17


def test_labels_are_unique(atoms):
    res = S.residues(atoms[atoms.chain_id == "A"])
    assert len(set(res.labels)) == len(res)


def test_sasa_decreases_on_binding(atoms):
    a = atoms[atoms.chain_id == "A"]
    b = atoms[atoms.chain_id == "X"]
    delta_a, delta_b = S.buried_on_binding(a, b)
    assert delta_a.sum() > 0 and delta_b.sum() > 0


@pytest.mark.parametrize("pdb_id,chain,lengths", [
    ("9GCN", "A", (10, 16, 17)),
    ("9S9E", "A", (10, 17, 12)),
    ("9RIU", "A", (13, 17, 12)),
])
def test_cdr_anchors_generalise(pdb_id, chain, lengths):
    """The framework anchors must work across VHH families, not just one."""
    atoms = S.read_atoms(S.fetch_structure(pdb_id))
    seq = S.residues(atoms[atoms.chain_id == chain]).sequence
    spans = S.find_cdrs(seq)
    got = tuple(spans[name][1] - spans[name][0] for name in ("CDR1", "CDR2", "CDR3"))
    assert got == lengths
