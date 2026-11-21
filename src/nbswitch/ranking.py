"""The final stage: combining the three scores without collapsing them.

A design has to satisfy three separate requirements, and a high score on one
cannot make up for failing another. A construct that cannot be expressed is
useless however well it masks, and a construct that cannot be cut is useless
however well it expresses. Averaging the three would let one compensate for
another, so instead the Pareto front is taken: a design is kept unless some
other design is at least as good on every score and better on one.
"""
from __future__ import annotations

import numpy as np

#: The three objectives, in the order they are reported.
OBJECTIVES = ("feasibility", "trigger", "unmasking")


def pareto_front(points: np.ndarray) -> np.ndarray:
    """Boolean mask of the non-dominated rows. Higher is better in every column."""
    n = len(points)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        better_or_equal = np.all(points >= points[i], axis=1)
        strictly_better = np.any(points > points[i], axis=1)
        dominated_by = better_or_equal & strictly_better
        if dominated_by.any():
            keep[i] = False
    return keep


def combine(feasibility_rows, trigger_rows, unmasking_rows=None) -> list[dict]:
    """Join the per-construct scores into one table, keyed on construct id."""
    table: dict[str, dict] = {}
    for r in feasibility_rows:
        table.setdefault(r["id"], {"id": r["id"], "mask": r["mask"],
                                   "kind": r["kind"],
                                   "linker_residues": r["linker_residues"]})
        table[r["id"]]["feasibility"] = r["feasibility_score"]
    for r in trigger_rows:
        if r["id"] in table:
            # cleavability is recognition times reach; fall back to reach alone
            # for constructs whose linker carries no selectivity score
            table[r["id"]]["trigger"] = (
                r.get("cleavability_score")
                if r.get("cleavability_score") is not None
                else r.get("context_score"))
    for r in (unmasking_rows or []):
        if r["id"] in table:
            table[r["id"]]["unmasking"] = r.get("occlusion")
    return list(table.values())


def rank(rows, objectives=None) -> dict:
    """Pareto front over the objectives, using only designs that carry all of them.

    A construct with no linker has no trigger score, so it is excluded rather
    than allowed to dominate on the objectives it happens to have.
    """
    if objectives is None:
        objectives = [o for o in OBJECTIVES
                      if sum(1 for r in rows if r.get(o) is not None) > len(rows) / 2]
    complete = [r for r in rows
                if all(r.get(o) is not None for o in objectives)]
    skipped = len(rows) - len(complete)
    if not objectives or not complete:
        return {"objectives": [], "front": [], "rows": rows, "n_skipped": skipped}
    rows = complete
    points = np.array([[float(r[o]) for o in objectives] for r in rows])
    keep = pareto_front(points)
    for r, on_front in zip(rows, keep):
        r["on_front"] = bool(on_front)
    front = [r for r, k in zip(rows, keep) if k]
    front.sort(key=lambda r: -sum(r[o] for o in objectives))
    return {"objectives": list(objectives), "n_total": len(rows),
            "n_front": len(front), "n_skipped": skipped,
            "front": front, "rows": rows}
