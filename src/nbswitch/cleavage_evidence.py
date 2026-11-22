"""An empirical cleavability score, to check the one the generator believes.

CleaveNet generates substrates and CleaveNet scores them. The two heads are
trained on the same data, so their agreement is not confirmation: a motif the
model likes will be generated and then rewarded, whether or not a protease
would cut it. Everything the pipeline says about cleavability so far rests on
that single model family.

This module scores the same substrates against the measurements the model was
trained on, by two methods that are not neural networks:

**Nearest measured neighbours.** For a candidate, find the substrates in the
measured library most similar to it and report what those were measured to do.
This makes no model of sequence at all. It answers a narrow question honestly:
when peptides of this shape were put in front of these enzymes, what happened?

**A position weight matrix.** Built as log-odds from the measured extremes,
the most selective substrates against the least. A PWM assumes positions act
independently, which is wrong for a protease that reads an extended groove,
but it is wrong in a completely different way from a transformer, so agreement
between them is informative.

**What this is not.** The measurements are the same ones CleaveNet learned
from, so this is an independent *estimator*, not independent *evidence*. It can
show that the neural scores are reproducible from the data by simpler means. It
cannot show that the data generalise to a snake venom metalloproteinase, which
is not in it.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from .linkers import PROTEASES, TARGET_PROTEASES, OFF_TARGET_PROTEASES, SUBSTRATE_LENGTH

#: The measured substrate library CleaveNet was trained on: 10-mers profiled
#: against eighteen MMPs (Kukreja et al.), shipped with the model.
MEASURED_CSV = Path("/scratch/cleavenet/data/kukreja.csv")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

#: Quantiles defining the two classes the weight matrix is built from.
TOP_QUANTILE = 0.95
BOTTOM_QUANTILE = 0.50

#: Neighbours averaged when reporting what similar measured substrates did.
N_NEIGHBOURS = 25


@lru_cache(maxsize=1)
def measured_library(path: Path = MEASURED_CSV) -> pd.DataFrame:
    """The measured substrates, with the same selectivity contrast used elsewhere."""
    frame = pd.read_csv(path)
    frame = frame.dropna(subset=["sequence"]).copy()
    frame["sequence"] = frame["sequence"].astype(str).str.strip().str.upper()
    frame = frame[frame.sequence.str.len() == SUBSTRATE_LENGTH]
    missing = [p for p in PROTEASES if p not in frame.columns]
    if missing:
        raise ValueError(f"measured library is missing {missing}")
    targets = [p for p in TARGET_PROTEASES if p in frame.columns]
    off = [p for p in OFF_TARGET_PROTEASES if p in frame.columns]
    frame["measured_target_z"] = frame[targets].mean(axis=1)
    frame["measured_off_target_z"] = frame[off].mean(axis=1)
    frame["measured_selectivity"] = (frame.measured_target_z
                                     - frame.measured_off_target_z)
    return frame.reset_index(drop=True)


@lru_cache(maxsize=1)
def _blosum62() -> tuple[dict, np.ndarray]:
    from biotite.sequence.align import SubstitutionMatrix
    from biotite.sequence import ProteinSequence
    matrix = SubstitutionMatrix.std_protein_matrix()
    alphabet = matrix.get_alphabet1()
    index = {str(sym): i for i, sym in enumerate(alphabet)}
    return index, matrix.score_matrix()


def similarity(a: str, b: str) -> float:
    """Ungapped BLOSUM62 similarity between two equal-length peptides."""
    index, score = _blosum62()
    return float(sum(score[index.get(x, 0), index.get(y, 0)] for x, y in zip(a, b)))


def nearest_measured(sequence: str, frame: pd.DataFrame | None = None,
                     k: int = N_NEIGHBOURS) -> dict:
    """What the most similar measured substrates were measured to do."""
    frame = measured_library() if frame is None else frame
    index, score = _blosum62()
    seq_idx = np.array([index.get(c, 0) for c in sequence])
    lib = np.array([[index.get(c, 0) for c in s] for s in frame.sequence])
    sims = score[seq_idx[None, :], lib].sum(axis=1)
    order = np.argsort(-sims)[:k]
    near = frame.iloc[order]
    return {
        "n_neighbours": int(k),
        "mean_similarity": round(float(sims[order].mean()), 2),
        "neighbour_measured_selectivity": round(float(near.measured_selectivity.mean()), 3),
        "neighbour_selectivity_sd": round(float(near.measured_selectivity.std()), 3),
        "best_neighbour": str(near.sequence.iloc[0]),
        "best_neighbour_selectivity": round(float(near.measured_selectivity.iloc[0]), 3),
    }


def weight_matrix(frame: pd.DataFrame | None = None,
                  top_q: float = TOP_QUANTILE,
                  bottom_q: float = BOTTOM_QUANTILE) -> np.ndarray:
    """Log-odds position weight matrix: selective substrates against the rest.

    Positions are treated as independent, which a protease groove is not. The
    point is that this assumption is wrong differently from a transformer's.
    """
    frame = measured_library() if frame is None else frame
    hi = frame[frame.measured_selectivity >= frame.measured_selectivity.quantile(top_q)]
    lo = frame[frame.measured_selectivity <= frame.measured_selectivity.quantile(bottom_q)]

    def counts(subset):
        c = np.ones((SUBSTRATE_LENGTH, len(AMINO_ACIDS)))      # Laplace prior
        for seq in subset.sequence:
            for i, aa in enumerate(seq[:SUBSTRATE_LENGTH]):
                j = AMINO_ACIDS.find(aa)
                if j >= 0:
                    c[i, j] += 1
        return c / c.sum(axis=1, keepdims=True)

    return np.log2(counts(hi) / counts(lo))


def pwm_score(sequence: str, matrix: np.ndarray | None = None) -> float:
    """Log-odds that a substrate came from the selective class."""
    matrix = weight_matrix() if matrix is None else matrix
    total = 0.0
    for i, aa in enumerate(sequence[:SUBSTRATE_LENGTH]):
        j = AMINO_ACIDS.find(aa)
        if j >= 0:
            total += float(matrix[i, j])
    return total


def evidence(sequences, k: int = N_NEIGHBOURS) -> list[dict]:
    """Both empirical scores for each substrate, plus the measured percentile."""
    frame = measured_library()
    matrix = weight_matrix(frame)
    rows = []
    for seq in sequences:
        row = {"substrate": seq,
               "pwm_score": round(pwm_score(seq, matrix), 3),
               "in_measured_library": bool((frame.sequence == seq).any())}
        row |= nearest_measured(seq, frame, k)
        rows.append(row)
    return rows


def calibrate(frame: pd.DataFrame | None = None, n: int = 2000,
              seed: int = 42) -> dict:
    """How well the empirical scores recover measured selectivity, held out.

    Without this the scores are just numbers. The matrix is built on one half
    of the library and asked about the other, so the correlation reported is
    out of sample.
    """
    frame = measured_library() if frame is None else frame
    rng = np.random.default_rng(seed)
    shuffled = frame.sample(frac=1.0, random_state=seed).reset_index(drop=True)
    cut = len(shuffled) // 2
    train, test = shuffled.iloc[:cut], shuffled.iloc[cut:]
    matrix = weight_matrix(train)
    sample = test.sample(n=min(n, len(test)), random_state=seed)
    scores = np.array([pwm_score(s, matrix) for s in sample.sequence])
    measured = sample.measured_selectivity.to_numpy()
    r = float(np.corrcoef(scores, measured)[0, 1])
    return {"n_train": int(len(train)), "n_test": int(len(sample)),
            "pwm_vs_measured_r": round(r, 3)}
