"""Folding the constructs with ESMFold.

The construct is a single polypeptide, so a single-sequence predictor applies
directly: no MSA, no complex assembly. ESMFold runs on CPU given enough memory,
which is why this demonstration can be reproduced without a GPU, and on a GPU it
is simply faster.
"""
from __future__ import annotations
import time
from pathlib import Path

import numpy as np

MODEL_NAME = "facebook/esmfold_v1"


def load_model(device: str | None = None, chunk_size: int = 64, half: bool | None = None):
    """Return (model, tokenizer, device). ``half`` defaults to True only on CUDA."""
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if half is None:
        half = device == "cuda"

    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    model = EsmForProteinFolding.from_pretrained(MODEL_NAME, low_cpu_mem_usage=True)
    model = model.to(device).eval()
    if half:
        model.esm = model.esm.half()
    if device == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = True
    if chunk_size:
        model.trunk.set_chunk_size(chunk_size)
    return model, tokenizer, device


def fold(sequence: str, model, tokenizer, device: str = "cpu") -> tuple[str, np.ndarray]:
    """Fold one sequence. Returns (PDB text, per-residue CA pLDDT)."""
    import torch
    from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
    from transformers.models.esm.openfold_utils.protein import Protein, to_pdb

    with torch.no_grad():
        tokens = tokenizer([sequence], return_tensors="pt",
                           add_special_tokens=False)["input_ids"].to(device)
        output = model(tokens)
    positions = atom14_to_atom37(output["positions"][-1], output).cpu().numpy()
    arrays = {k: v.cpu().numpy() for k, v in output.items()}
    protein = Protein(aatype=arrays["aatype"][0], atom_positions=positions[0],
                      atom_mask=arrays["atom37_atom_exists"][0],
                      residue_index=arrays["residue_index"][0] + 1,
                      b_factors=arrays["plddt"][0],
                      chain_index=arrays.get("chain_index", [None])[0])
    plddt = arrays["plddt"][0, :, 1]
    if float(np.nanmax(plddt)) <= 1.0:      # some builds return 0-1, normalise to 0-100
        plddt = plddt * 100.0
    return to_pdb(protein), plddt


def fold_library(library, out_dir: Path, model=None, tokenizer=None, device=None,
                 skip_existing: bool = True, log=print) -> dict:
    """Fold every construct, writing one PDB each. Resumable."""
    out_dir.mkdir(parents=True, exist_ok=True)
    if model is None:
        model, tokenizer, device = load_model(device)
    plddt = {}
    start = time.time()
    for i, item in enumerate(library, 1):
        path = out_dir / f"{item['id']}.pdb"
        npy = out_dir / f"{item['id']}.plddt.npy"
        if skip_existing and path.exists() and npy.exists():
            plddt[item["id"]] = np.load(npy)
            log(f"[{i:2d}/{len(library)}] {item['id']:16s} cached")
            continue
        text, scores = fold(item["sequence"], model, tokenizer, device)
        path.write_text(text)
        np.save(npy, scores)
        plddt[item["id"]] = scores
        log(f"[{i:2d}/{len(library)}] {item['id']:16s} {item['length']:3d} aa  "
            f"pLDDT {scores.mean():5.1f}  ({time.time() - start:6.0f}s)")
    return plddt
