"""The fold cache must key on the sequence, not on the construct id.

Construct ids are reused across rebuilds. ``VHH_only`` keeps its name while the
sequence under it changes - trimming an expression tag shortened it from 127 to
123 residues - and an id-keyed cache then serves a model of a different
molecule. That model is the reference every occlusion value is measured
against, so the error does not announce itself: it shifts a baseline.
"""
from pathlib import Path

import numpy as np
import pytest

from nbswitch import folding


SEQ = "QVQLQESGGGLVRPGGSLRLSCAVS"


def _write_model(directory: Path, name: str, sequence: str) -> None:
    """A minimal CA-only PDB plus its pLDDT sidecar, as the cache expects."""
    from biotite.sequence import ProteinSequence

    lines = []
    for i, aa in enumerate(sequence, 1):
        three = ProteinSequence.convert_letter_1to3(aa)
        lines.append(f"ATOM  {i:5d}  CA  {three} A{i:4d}    "
                     f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 80.00           C")
    (directory / f"{name}.pdb").write_text("\n".join(lines) + "\nEND\n")
    np.save(directory / f"{name}.plddt.npy", np.full(len(sequence), 80.0))


def test_matching_sequence_is_served_from_cache(tmp_path):
    _write_model(tmp_path, "VHH_only", SEQ)
    item = {"id": "VHH_only", "sequence": SEQ, "length": len(SEQ)}

    def refuse(*args, **kwargs):
        raise AssertionError("a cached model that matches must not be refolded")

    plddt = folding.fold_library([item], tmp_path, model=refuse, tokenizer=None,
                                 device="cpu", log=lambda *a: None)
    assert plddt["VHH_only"].shape == (len(SEQ),)


def test_a_changed_sequence_under_the_same_id_is_refolded(tmp_path, monkeypatch):
    # on disk: the untrimmed construct, four residues of expression tag longer
    _write_model(tmp_path, "VHH_only", SEQ + "ENLY")
    item = {"id": "VHH_only", "sequence": SEQ, "length": len(SEQ)}

    calls = []

    def fake_fold(sequence, model, tokenizer, device):
        calls.append(sequence)
        return "END\n", np.full(len(sequence), 75.0)

    monkeypatch.setattr(folding, "fold", fake_fold)
    plddt = folding.fold_library([item], tmp_path, model=object(),
                                 tokenizer=None, device="cpu",
                                 log=lambda *a: None)
    assert calls == [SEQ], "the stale model was served instead of being refolded"
    assert plddt["VHH_only"].shape == (len(SEQ),)


def test_pdb_sequence_reads_what_is_in_the_file(tmp_path):
    _write_model(tmp_path, "x", SEQ)
    assert folding.pdb_sequence(tmp_path / "x.pdb") == SEQ
