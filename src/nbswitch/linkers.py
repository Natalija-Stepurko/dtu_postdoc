"""Protease-cleavable linkers, generated under a cleavability constraint.

The first version of this pipeline used (GGGGS)n linkers and then scored how
reachable they were. That measures the wrong thing: a Gly/Ser repeat is chosen
in protein engineering precisely because proteases do not cut it, so no score
applied afterwards could have found a cleavable one. Nothing in the library
could have worked.

Cleavability is therefore a constraint on generation, not a filter applied
after it. CleaveNet is a conditional generator: it takes a target profile of
z-scores across eighteen matrix metalloproteinases and emits substrates
matching that profile. Asking it for substrates that the bite-site proteases
cut and the circulating ones do not means every candidate entering the library
is cleavable by construction. The scoring module downstream then verifies and
ranks, rather than searching a space that contains no solutions.

The requirement is a contrast, not a single number. A linker cut by everything
is useless: the construct would unmask in circulation, which is the failure the
whole design exists to prevent. What is wanted is a linker cut at the site and
resistant elsewhere, so the criterion is the gap between the two groups.

**Boundary of validity.** CleaveNet is trained on human MMPs. The protease that
actually arrives first at a bite site is a snake venom metalloproteinase, which
is not in its training set and is not a human MMP, so no score here is a
prediction about SVMP cleavage. The host gelatinases are used as a stand-in
because they are documented to rise in envenomed tissue, on a slower clock than
the venom's own proteases. Any number this module produces is a statement about
human MMPs and has to be read that way.
"""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

#: The eighteen proteases CleaveNet models, in the order its weights expect.
PROTEASES = ["MMP1", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP15",
             "MMP16", "MMP17", "MMP19", "MMP2", "MMP20", "MMP24", "MMP25",
             "MMP3", "MMP7", "MMP8", "MMP9"]

#: Wanted at the bite site. MMP2 and MMP9 are the gelatinases whose activity
#: rises in envenomed tissue; MMP8 is neutrophil collagenase, released by the
#: infiltrate. These stand in for the venom metalloproteinases, which CleaveNet
#: cannot model. See the module docstring.
TARGET_PROTEASES = ("MMP2", "MMP9", "MMP8")

#: Not wanted. The membrane-type MMPs are cell-surface and widely distributed,
#: and MMP1, MMP3 and MMP7 are broadly present in healthy tissue, so cleavage
#: by these is what would unmask the construct before it arrives.
OFF_TARGET_PROTEASES = ("MMP1", "MMP3", "MMP7", "MMP14", "MMP15", "MMP16",
                        "MMP17", "MMP24", "MMP25")

#: Flanks between the mask, the substrate and the nanobody. Short Gly/Ser runs
#: keep the substrate accessible without contributing a cleavage site of their
#: own, which is the one job a Gly/Ser repeat is right for.
FLANK = "GGS"

CLEAVENET = Path("/scratch/cleavenet")


def conditioning_profile(target_z: float = 2.0,
                         off_target_z: float = -1.0) -> pd.DataFrame:
    """The z-score profile to generate against: cut by the first group, not the second.

    Exactly one row. CleaveNet's generator loops over every row of this file
    and draws ``--num-seqs`` samples for each, reloading its weights on every
    draw, so a file with n rows asks for n x num_seqs sequences rather than
    num_seqs of them.

    Proteases in neither group are left at zero, meaning no preference
    expressed, rather than being pushed either way on evidence we do not have.
    """
    row = {}
    for protease in PROTEASES:
        if protease in TARGET_PROTEASES:
            row[protease] = target_z
        elif protease in OFF_TARGET_PROTEASES:
            row[protease] = off_target_z
        else:
            row[protease] = 0.0
    return pd.DataFrame([row])


def generate(n: int, temperature: float = 0.9, seed_dir: Path | None = None,
             python: str = "python") -> list[str]:
    """Generate substrates conditioned on the profile. Returns unique sequences."""
    work = Path(seed_dir or tempfile.mkdtemp(prefix="cleavenet-"))
    work.mkdir(parents=True, exist_ok=True)
    profile = work / "profile.csv"
    conditioning_profile().to_csv(profile, index=False)
    out = work / "generated"
    subprocess.run(
        [python, "src/generate.py", "--num-seqs", str(n),
         "--temperature", str(temperature), "--z-scores", str(profile),
         "--output-dir", str(out)],
        cwd=CLEAVENET, check=True, capture_output=True)

    seen, sequences = set(), []
    for path in sorted(out.rglob("*.csv")):
        for line in path.read_text().splitlines():
            peptide = line.strip().split(",")[0].strip().upper()
            if peptide and peptide.isalpha() and peptide not in seen:
                seen.add(peptide)
                sequences.append(peptide)
    return sequences


def predict(sequences, work_dir: Path | None = None,
            python: str = "python") -> pd.DataFrame:
    """Predicted z-score per protease for each substrate, from the same model family."""
    work = Path(work_dir or tempfile.mkdtemp(prefix="cleavenet-pred-"))
    work.mkdir(parents=True, exist_ok=True)
    seq_csv = work / "sequences.csv"
    seq_csv.write_text("\n".join(sequences) + "\n")
    out = work / "predictions"
    subprocess.run(
        [python, "src/predict.py", "--path-to-sequence-csv", str(seq_csv),
         "--no-csv-header", "--save-dir", str(out)],
        cwd=CLEAVENET, check=True, capture_output=True)

    scores = next(out.rglob("all_scores.csv"))
    frame = pd.read_csv(scores, header=None)
    frame.columns = PROTEASES[:frame.shape[1]]
    frame.insert(0, "sequence", sequences[:len(frame)])
    return frame


def selectivity(frame: pd.DataFrame) -> pd.DataFrame:
    """Score each substrate on the contrast the design actually needs.

    ``target`` is mean predicted cleavage by the bite-site group and ``off``
    the mean over the circulating group. ``selectivity`` is their difference:
    a substrate cut by everything scores near zero however high its absolute
    cleavage, which is the correct outcome, since it would unmask in transit.
    """
    frame = frame.copy()
    targets = [p for p in TARGET_PROTEASES if p in frame.columns]
    off = [p for p in OFF_TARGET_PROTEASES if p in frame.columns]
    frame["target_z"] = frame[targets].mean(axis=1)
    frame["off_target_z"] = frame[off].mean(axis=1)
    frame["selectivity"] = frame["target_z"] - frame["off_target_z"]
    return frame.sort_values("selectivity", ascending=False).reset_index(drop=True)


def build_linker(substrate: str, flank: str = FLANK) -> str:
    """A substrate presented on short flexible flanks."""
    return f"{flank}{substrate}{flank}"


def cleavable_linkers(n_generate: int = 120, n_keep: int = 10,
                      temperature: float = 0.9, work_dir: Path | None = None,
                      python: str = "python") -> dict:
    """Generate, verify, and keep the most selective substrates.

    Generation and verification use different CleaveNet heads, so the ranking
    is not simply reading back the conditioning tag, but they share training
    data and the agreement between them is weaker evidence than two
    independent models would be.
    """
    work = Path(work_dir or tempfile.mkdtemp(prefix="linkers-"))
    sequences = generate(n_generate, temperature=temperature,
                         seed_dir=work / "gen", python=python)
    if not sequences:
        raise RuntimeError("CleaveNet produced no sequences")
    scored = selectivity(predict(sequences, work_dir=work / "pred", python=python))
    kept = scored.head(n_keep)
    return {
        "n_generated": len(sequences),
        "n_kept": len(kept),
        "target_proteases": list(TARGET_PROTEASES),
        "off_target_proteases": list(OFF_TARGET_PROTEASES),
        "flank": FLANK,
        "linkers": [
            {"id": f"cn{i:02d}",
             "substrate": row.sequence,
             "sequence": build_linker(row.sequence),
             "target_z": round(float(row.target_z), 3),
             "off_target_z": round(float(row.off_target_z), 3),
             "selectivity": round(float(row.selectivity), 3)}
            for i, row in enumerate(kept.itertuples(), 1)],
        "all_scored": [
            {"substrate": row.sequence,
             "selectivity": round(float(row.selectivity), 3)}
            for row in scored.itertuples()],
    }
