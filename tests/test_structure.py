"""Residue keying is the one place a silent error would corrupt everything
downstream, so it is pinned here."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from nbswitch import structure as S
from nbswitch.config import PDB_ID, REFERENCE_VHH_CHAIN


@pytest.fixture(scope="module")
def atoms():
    return S.read_atoms(S.fetch_structure(PDB_ID))


def test_insertion_codes_are_not_collapsed(atoms):
    """9GCN is Kabat-numbered: 122 residues over 112 distinct res_id values.
    Keying on res_id alone merges residues and truncates CDR3."""
    chain = atoms[atoms.chain_id == REFERENCE_VHH_CHAIN]
    res = S.residues(chain)
    assert len(res) == 122
    assert len(set(chain.res_id.tolist())) == 112


def test_cdr3_is_full_length(atoms):
    chain = atoms[atoms.chain_id == REFERENCE_VHH_CHAIN]
    spans = S.find_cdrs(S.residues(chain).sequence)
    start, stop = spans["CDR3"]
    assert stop - start == 17


def test_labels_are_unique(atoms):
    res = S.residues(atoms[atoms.chain_id == REFERENCE_VHH_CHAIN])
    assert len(set(res.labels)) == len(res)


def test_sasa_decreases_on_binding(atoms):
    a = atoms[atoms.chain_id == "A"]
    b = atoms[atoms.chain_id == "X"]
    da, db = S.buried_on_binding(a, b)
    assert da.sum() > 0 and db.sum() > 0
