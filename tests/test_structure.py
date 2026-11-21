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
    ("9GCN", "A", (8, 7, 17)),
    ("9S9E", "A", (8, 8, 12)),
    ("9RIU", "A", (11, 7, 12)),
])
def test_cdr_lengths_under_imgt(pdb_id, chain, lengths):
    """CDR spans come from a published scheme, so they are pinned per scheme.

    These are IMGT lengths. Changing CDR_SCHEME changes them, which is the
    point: a CDR count is meaningless without the scheme that produced it.
    """
    atoms = S.read_atoms(S.fetch_structure(pdb_id))
    seq = S.residues(atoms[atoms.chain_id == chain]).sequence
    spans = S.find_cdrs(seq, scheme="imgt")
    got = tuple(spans[name][1] - spans[name][0] for name in ("CDR1", "CDR2", "CDR3"))
    assert got == lengths


@pytest.mark.parametrize("pdb_id,chain,lengths", [
    ("9GCN", "A", (10, 16, 17)),
    ("9S9E", "A", (10, 17, 12)),
    ("9RIU", "A", (13, 17, 12)),
])
def test_anchor_fallback_still_works(pdb_id, chain, lengths):
    """The regex fallback is used when the numbering tool is unavailable.

    It does not reproduce any single published scheme, so it is pinned only
    to catch it silently breaking, never used to make a claim.
    """
    atoms = S.read_atoms(S.fetch_structure(pdb_id))
    seq = S.residues(atoms[atoms.chain_id == chain]).sequence
    spans = S._cdr_spans_anchors(seq, S.CDR_ANCHORS)
    got = tuple(spans[name][1] - spans[name][0] for name in ("CDR1", "CDR2", "CDR3"))
    assert got == lengths


def test_schemes_disagree_so_the_scheme_must_be_named():
    """Kabat and IMGT give different CDR2 lengths on the worked example."""
    atoms = S.read_atoms(S.fetch_structure("9RIU"))
    seq = S.residues(atoms[atoms.chain_id == "A"]).sequence
    kabat = S.find_cdrs(seq, scheme="kabat")
    imgt = S.find_cdrs(seq, scheme="imgt")
    k2 = kabat["CDR2"][1] - kabat["CDR2"][0]
    i2 = imgt["CDR2"][1] - imgt["CDR2"][0]
    assert k2 == 17 and i2 == 7, (k2, i2)


def test_contact_cutoff_excludes_hydrogens_and_is_tighter_than_burial():
    """The 3 A contact set is a subset of the burial footprint, and riding
    hydrogens in the deposited model must not affect it."""
    import numpy as np
    from nbswitch import structure as S
    from nbswitch.config import load_json
    ref = load_json("01_reference.json")
    pair = ref["chain_pairs"]["A"]
    contacts = set(pair["paratope_contact_labels"])
    buried = set(pair["paratope_labels"])
    assert contacts, "a 3 A criterion should still select something"
    assert contacts <= buried, "a contacting residue must also bury surface"
    assert len(contacts) < len(buried), "3 A is tighter than burial"
    from nbswitch.config import PDB_ID
    atoms = S.read_atoms(S.fetch_structure(PDB_ID))
    assert (atoms.element == "H").any(), "this model carries riding hydrogens"
    a = atoms[atoms.chain_id == "A"]
    b = atoms[atoms.chain_id == pair["toxin_chain"]]
    d_with, _ = S.contact_residues(a, b, 3.0)
    d_without, _ = S.contact_residues(a[a.element != "H"], b[b.element != "H"], 3.0)
    assert np.allclose(d_with, d_without, equal_nan=True)
