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
        # carried so the hard constraints can be checked on the molecule
        # rather than on the heuristic that summarises it
        for key in ("plddt_vhh", "vhh_core_rmsd_A", "unpaired_cysteines"):
            if key in r:
                table[r["id"]][key] = r[key]
    for r in trigger_rows:
        if r["id"] in table:
            # recognition times reach. A construct whose linker carries no
            # selectivity score has no value for this objective: substituting
            # reach alone would compare it against others on a different
            # quantity, which is worse than leaving it out, since rank()
            # already excludes designs that are missing an objective.
            table[r["id"]]["cleavability"] = r.get("cleavability_score")
            for key in ("selectivity", "context_score"):
                if key in r:
                    table[r["id"]][key] = r[key]
    for r in (unmasking_rows or []):
        if r["id"] in table:
            table[r["id"]]["surface_occlusion"] = r.get("occlusion")
    return list(table.values())


def rank(rows, objectives=None, penalty: float = 0.0, n_shortlist: int = 6) -> dict:
    """Hard constraints, then the nondominated set, then a stated ordering.

    Three steps, each doing one thing:

    **Constraints.** Drop what cannot work regardless of its other scores. A
    Pareto front would otherwise let a construct survive on two objectives
    while failing a requirement outright.

    **Non-domination.** Keep a design unless another is at least as good on
    every objective and better on one. No weights, so no score compensates.

    **Ordering.** The front is a set, so turning it into a list needs a rule.
    The rule is maximin on the normalised scores: a construct is worth no more
    than its weakest requirement. It is reported as the rule it is.

    A design with any objective missing is excluded rather than ranked on the
    rest, which would compare candidates on different requirements.
    """
    if objectives is None:
        missing = [o for o in OBJECTIVES
                   if not any(r.get(o) is not None for r in rows)]
        objectives = [] if missing else list(OBJECTIVES)

    scored = [r for r in rows if all(r.get(o) is not None for o in objectives)]
    n_unscored = len(rows) - len(scored)
    if not objectives or not scored:
        return {"objectives": [], "front": [], "rows": rows,
                "n_skipped": n_unscored}

    feasible, excluded = apply_hard_constraints(scored)
    if not feasible:
        return {"objectives": list(objectives), "front": [], "rows": scored,
                "n_total": len(scored), "n_front": 0,
                "n_skipped": n_unscored, "n_excluded": len(excluded),
                "excluded": excluded}

    points = np.array([[float(r[o]) for o in objectives] for r in feasible])
    keep = pareto_front(points)
    for r, on_front in zip(feasible, keep):
        r["on_front"] = bool(on_front)
    for r in excluded:
        r["on_front"] = False
    front = [r for r, k in zip(feasible, keep) if k]

    scale = normalise(feasible, objectives)
    for r in feasible:
        r["robust_score"] = round(robust_score(r, scale, objectives, penalty), 4)
    ordered = sorted(front, key=lambda r: -r["robust_score"])
    shortlist = select_diverse(ordered, n=n_shortlist)

    from collections import Counter
    why = Counter(k for r in excluded for k in r["failed_constraints"])

    return {
        "objectives": list(objectives),
        "objective_basis": {o: OBJECTIVE_BASIS[o] for o in objectives},
        "n_total": len(scored), "n_skipped": n_unscored,
        "n_feasible": len(feasible), "n_excluded": len(excluded),
        "excluded_by": dict(why),
        "constraints": [{"field": f, "op": o, "limit": l, "why": w}
                        for f, o, l, w in HARD_CONSTRAINTS],
        "n_front": len(front),
        "front_order": "maximin on the normalised objectives: ordered by the "
                       "weakest of the three. The front itself is a set; this "
                       "is a stated rule applied to it, not a property of it.",
        "normalisation": "min-max within the feasible set, so the three scales "
                         "are comparable. Not a calibration against measured "
                         "quantities, which do not exist for these scores.",
        "uncertainty_penalty": penalty,
        "front": ordered,
        "shortlist": [r["id"] for r in shortlist],
        "rows": feasible + excluded,
    }


# --------------------------------------------------------------------------
# Constrained Pareto selection, then a stated ordering rule
#
# A front is a set. Turning it into a list needs a rule, and the rule has to
# be said out loud rather than hidden in a sort key. What follows is three
# steps, in order: drop what cannot work at all, take the nondominated set of
# what remains, then order that set by its weakest normalised objective.
# --------------------------------------------------------------------------

#: Requirements a construct has to meet before its scores are worth comparing.
#: Each is a property of the molecule rather than a quantile of the observed
#: distribution, so the set does not move when the library changes.
HARD_CONSTRAINTS = (
    ("plddt_vhh", "ge", 70.0,
     "the nanobody itself must be confidently predicted, or nothing measured "
     "on its paratope means anything"),
    ("vhh_core_rmsd_A", "le", 2.0,
     "the fold has to survive the fusion; past about 2 A the core has moved"),
    ("selectivity", "gt", 0.0,
     "a linker cut no more at the site than in circulation unmasks in transit, "
     "which is the failure the design exists to prevent"),
    ("context_score", "gt", 0.0,
     "a linker a protease cannot reach is not cleavable at any sequence"),
)


def apply_hard_constraints(rows, constraints=HARD_CONSTRAINTS):
    """Split candidates into those worth ranking and those ruled out.

    A construct that fails one of these is not a worse candidate, it is not a
    candidate: no score on another axis can compensate, which is exactly what
    a Pareto front would otherwise let it do.
    """
    ops = {"ge": lambda a, b: a >= b, "le": lambda a, b: a <= b,
           "gt": lambda a, b: a > b, "lt": lambda a, b: a < b}
    feasible, excluded = [], []
    for r in rows:
        failed = []
        for key, op, limit, _why in constraints:
            value = r.get(key)
            if value is None:
                continue                      # not measured is not a failure
            if not ops[op](float(value), limit):
                failed.append(key)
        (excluded if failed else feasible).append(
            {**r, "failed_constraints": failed} if failed else r)
    return feasible, excluded


def normalise(rows, objectives=OBJECTIVES) -> dict:
    """Min-max each objective over the candidates being compared.

    The three scores are on unrelated scales, so comparing their raw values
    would make the ordering an artefact of which happens to be numerically
    smallest. This puts them on 0 to 1 within the set being ranked.

    It is a within-run normalisation and not a calibration: it says where a
    construct sits among those tested, not whether any of them is good.
    """
    out = {}
    for o in objectives:
        vals = [float(r[o]) for r in rows if r.get(o) is not None]
        lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
        span = hi - lo
        out[o] = (lo, span if span > 0 else 1.0)
    return out


def robust_score(row, scale, objectives=OBJECTIVES, penalty: float = 0.0,
                 uncertainty=None) -> float:
    """The weakest normalised objective, optionally discounted for uncertainty.

    Maximin: a construct is worth no more than its worst requirement, because
    all three have to hold at once. Ordering by this prefers the balanced
    candidate over one that is excellent twice and poor once.

    ``penalty`` multiplies a per-objective uncertainty before taking the
    minimum, which stops a design being promoted by one optimistic model. No
    calibrated uncertainty exists for these scores yet, so it defaults to
    nothing and the hook is here for when one does.
    """
    uncertainty = uncertainty or {}
    worst = None
    for o in objectives:
        value = row.get(o)
        if value is None:
            return float("-inf")
        lo, span = scale[o]
        z = (float(value) - lo) / span
        z -= penalty * float(uncertainty.get(o, 0.0))
        worst = z if worst is None else min(worst, z)
    return float(worst if worst is not None else float("-inf"))


def select_diverse(rows, n: int = 6, key=lambda r: r.get("mask")):
    """Take the best of each distinct design, then fill up by order.

    Predictions made for near-identical constructs share whatever systematic
    error they carry, so a shortlist of six variants of one mask tests one
    hypothesis six times. ``rows`` must already be ordered.
    """
    seen, chosen = set(), []
    for r in rows:
        k = key(r)
        if k not in seen:
            seen.add(k)
            chosen.append(r)
        if len(chosen) >= n:
            return chosen
    for r in rows:
        if r not in chosen:
            chosen.append(r)
        if len(chosen) >= n:
            break
    return chosen
