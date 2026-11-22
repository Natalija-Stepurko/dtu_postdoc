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
OBJECTIVES = ("feasibility", "cleavability", "surface_occlusion")

#: What each objective is actually computed from, as opposed to what it is
#: meant to capture. These are not the same thing here, and the difference is
#: the main limitation of the ranking. The third objective is named for the
#: quantity it computes, surface occlusion, and not for the property it was
#: meant to stand in for. The controls show it does not separate designed
#: masks from scrambled or inert ones, so a front built on it is not a list of
#: constructs predicted to mask.
OBJECTIVE_BASIS = {
    "feasibility": "a heuristic combination of predicted fold confidence, "
                   "nanobody preservation, exposed hydrophobic surface and "
                   "measured disulfide pairing",
    "cleavability": "predicted human-MMP selectivity times structural reach of "
                    "the linker; not a cleavage rate, and not specificity for "
                    "a snake venom metalloproteinase",
    "surface_occlusion": "paratope surface lost in the folded construct. It is "
                         "named for what it measures: the controls show it does "
                         "not separate designed masks from scrambled or inert "
                         "ones, so it is not a masking score and a front built "
                         "on it is not a list of masked nanobodies",
}


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
            # recognition times reach. A construct whose linker carries no
            # selectivity score has no value for this objective: substituting
            # reach alone would compare it against others on a different
            # quantity, which is worse than leaving it out, since rank()
            # already excludes designs that are missing an objective.
            table[r["id"]]["cleavability"] = r.get("cleavability_score")
    for r in (unmasking_rows or []):
        if r["id"] in table:
            table[r["id"]]["surface_occlusion"] = r.get("occlusion")
    return list(table.values())


def rank(rows, objectives=None) -> dict:
    """Pareto front over the objectives, using only designs that carry all of them.

    A construct with no linker has no trigger score, so it is excluded rather
    than allowed to dominate on the objectives it happens to have.
    """
    # All three or none. Dropping an objective because it is missing compares
    # candidates on different requirements and calls the result a front.
    if objectives is None:
        missing = [o for o in OBJECTIVES
                   if not any(r.get(o) is not None for r in rows)]
        objectives = [] if missing else list(OBJECTIVES)
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
    # The front is a set, not a ranking. Ordering it by the sum applies equal
    # weights after saying none are used; ordering it by the weakest objective
    # compares scores on different scales. Neither is defensible, so it is
    # returned in a stable arbitrary order and labelled unordered. Choosing
    # between its members needs a stated rule: a feasibility threshold,
    # predictive uncertainty, or diversity.
    front.sort(key=lambda r: r["id"])
    return {"objectives": list(objectives), "n_total": len(rows),
            "n_front": len(front), "n_skipped": skipped,
            "front_order": "unordered; sorted by id for stability only",
            "objective_basis": {o: OBJECTIVE_BASIS[o] for o in objectives},
            "front": front, "rows": rows}
