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

import shutil
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

#: The model's positional embedding is sized for this; longer input fails
#: inside TensorFlow with a shape error rather than a useful message.
SUBSTRATE_LENGTH = 10

#: Flanks between the mask, the substrate and the nanobody. Inert Gly/Ser runs,
#: which contribute no cleavage site of their own: being uncuttable is the one
#: job a Gly/Ser repeat is right for.
#:
#: Length is set by measurement, not by convention. Stage 16 holds the
#: substrate constant and varies only these flanks: with three residues a side,
#: giving a 16-residue linker, the mask sits 26 A from the paratope in the
#: predicted fold and never contacts it, so every masking score computed on
#: such a construct is measuring nothing. Contact appears at seven residues a
#: side, a 24-residue linker, and costs 0.08 decades of the affinity window
#: against the worm-like-chain optimum at 14.
#:
#: The gap between those two numbers is itself a result: the tether model puts
#: the mask at the end of the linker as a point, and the mask is a peptide with
#: its own dimensions that has to reach across the paratope rather than touch
#: its centre.
FLANK = "GGSGGSG"        # seven a side, giving REFERENCE_LENGTH

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
    # CleaveNet opens its output file in append mode, so a work directory that
    # survives between runs accumulates every draw ever made there and the
    # reported count becomes a function of how often the stage was run.
    if out.exists():
        shutil.rmtree(out)
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
    wrong = sorted({len(x) for x in sequences} - {SUBSTRATE_LENGTH})
    if wrong:
        raise ValueError(
            f"CleaveNet scores {SUBSTRATE_LENGTH}-residue substrates only; got "
            f"lengths {wrong}. Pass the substrate, not the flanked linker.")
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
    # the file carries a header and its own sequence column, and the model
    # writes its start token back out in front of each peptide
    frame = pd.read_csv(scores)
    first = frame.columns[0]
    frame = frame.rename(columns={first: "sequence"})
    frame["sequence"] = frame["sequence"].astype(str).str.lstrip("!").str.strip()
    missing = [p for p in PROTEASES if p not in frame.columns]
    if missing:
        raise ValueError(f"predictor output is missing {missing}")
    return frame[["sequence"] + PROTEASES]


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


#: Total linker lengths to model. The substrate is fixed at ten residues by
#: CleaveNet, so length is varied in the flanks, which are inert Gly/Ser and
#: contribute no cleavage site of their own. This separates two things that
#: were previously one number: whether a protease recognises the linker, and
#: whether the mask can reach the paratope.
#:
#: The range is bounded at both ends by the two models that disagree about it.
#: Below 16 the worm-like chain says the mask is held too tightly to reach past
#: its own attachment point, and at 16 the folded constructs show it 26 A from
#: the paratope. Above 32 the affinity window has lost a fifth of its width and
#: the chain is long enough to wrap the nanobody, which masks by crowding
#: rather than by binding.
LENGTH_SERIES = (12, 16, 20, 24, 28, 32)

#: The length every substrate is built at, so sequences can be compared with
#: length held constant. Measured: stage 16 finds no mask-paratope contact
#: below this.
REFERENCE_LENGTH = 24


def flank_of_length(n: int) -> str:
    """An inert Gly/Ser flank of exactly ``n`` residues."""
    return ("GGS" * (n // 3 + 1))[:n]


def length_variants(substrate: str, linker_id: str,
                    lengths=LENGTH_SERIES, **extra) -> list[dict]:
    """One linker per total length, all carrying the same substrate.

    Only the flanks change, so any difference in masking between these is
    attributable to how far the mask can reach, not to its cleavability.
    """
    out = []
    for total in lengths:
        pad = total - len(substrate)
        if pad < 0 or pad % 2:
            raise ValueError(
                f"total length {total} cannot flank a {len(substrate)}-residue "
                f"substrate symmetrically")
        flank = flank_of_length(pad // 2)
        out.append({"id": f"{linker_id}_n{total}",
                    "substrate": substrate,
                    "sequence": f"{flank}{substrate}{flank}",
                    "flank": flank, "total_length": total, **extra})
    return out
