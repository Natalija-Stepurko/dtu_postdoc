"""Redesigning the mask on a folded backbone with ProteinMPNN.

The grafted mask fixes the design problem before it starts: its sequence comes
from the antigen, not from anything about this paratope. Taking the folded
fusion as a backbone and letting ProteinMPNN rewrite only the mask positions -
with the nanobody and linker held fixed - turns the graft into a starting point
rather than the answer.

ProteinMPNN is small enough to run on CPU in seconds per design.
"""
from __future__ import annotations
import json, subprocess, tempfile
from pathlib import Path

PROTEINMPNN = Path("/scratch/ProteinMPNN")


def available() -> bool:
    return (PROTEINMPNN / "protein_mpnn_run.py").exists()


def design_mask(pdb_path: Path, mask_range, out_dir: Path, chain: str = "A",
                n_sequences: int = 8, temperature: float = 0.2, seed: int = 42,
                python: str = "python") -> list[dict]:
    """Rewrite only ``mask_range`` (1-based inclusive) on a folded construct.

    Everything outside that span is held fixed, so the nanobody fold and the
    linker are carried over untouched.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = pdb_path.stem
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        designed = list(range(mask_range[0], mask_range[1] + 1))
        # ProteinMPNN's fixed-position dict lists the positions to KEEP.
        n_res = _n_residues(pdb_path, chain)
        fixed = [i for i in range(1, n_res + 1) if i not in designed]
        (tmp / "fixed.jsonl").write_text(json.dumps({stem: {chain: fixed}}) + "\n")
        (tmp / "chains.jsonl").write_text(json.dumps({stem: [[chain], []]}) + "\n")
        pdb_dir = tmp / "pdbs"; pdb_dir.mkdir()
        (pdb_dir / pdb_path.name).write_text(pdb_path.read_text())
        cmd = [python, str(PROTEINMPNN / "protein_mpnn_run.py"),
               "--pdb_path", str(pdb_dir / pdb_path.name),
               "--pdb_path_chains", chain,
               "--fixed_positions_jsonl", str(tmp / "fixed.jsonl"),
               "--out_folder", str(tmp / "out"),
               "--num_seq_per_target", str(n_sequences),
               "--sampling_temp", str(temperature),
               "--seed", str(seed), "--batch_size", "1"]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        fasta = next((tmp / "out" / "seqs").glob("*.fa"))
        return _parse_fasta(fasta.read_text(), stem)


def _n_residues(pdb_path: Path, chain: str) -> int:
    seen = set()
    for line in pdb_path.read_text().splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            seen.add(line[22:27])
    return len(seen)


def _parse_fasta(text: str, stem: str) -> list[dict]:
    """ProteinMPNN writes the native sequence first, then the designs."""
    out, header = [], None
    for line in text.splitlines():
        if line.startswith(">"):
            header = line[1:]
        elif line.strip() and header is not None:
            fields = dict(
                part.split("=", 1) for part in header.split(", ") if "=" in part)
            out.append({"parent": stem, "sequence": line.strip(),
                        "score": float(fields.get("score", "nan")),
                        "global_score": float(fields.get("global_score", "nan")),
                        "sample": fields.get("sample", "native"),
                        "is_native": "sample" not in fields})
            header = None
    return out
