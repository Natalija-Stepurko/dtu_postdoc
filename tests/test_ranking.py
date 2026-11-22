"""Constrained Pareto selection, then a stated ordering rule."""
import numpy as np
import pytest

from nbswitch import ranking as R


def _row(name, s1, s2, s3, **over):
    base = dict(id=name, mask=name, kind="epitope graft", linker_residues=24,
                feasibility=s1, cleavability=s2, surface_occlusion=s3,
                plddt_vhh=85.0, vhh_core_rmsd_A=0.3,
                selectivity=1.2, context_score=0.3)
    base.update(over)
    return base


#: The comparison in the review: B dominates D, and A, B, C each keep a
#: different trade-off.
WORKED = [_row("A", 0.95, 0.40, 0.85), _row("B", 0.85, 0.78, 0.80),
          _row("C", 0.80, 0.85, 0.60), _row("D", 0.75, 0.70, 0.70)]


def test_domination_is_all_objectives_at_least_as_good_and_one_better():
    out = R.rank(list(WORKED))
    front = {r["id"] for r in out["front"]}
    assert front == {"A", "B", "C"}, "D is dominated by B and must be dropped"


def test_maximin_prefers_the_balanced_candidate():
    out = R.rank(list(WORKED))
    assert out["front"][0]["id"] == "B"


def test_a_constraint_failure_is_not_a_low_score():
    """A broken fold cannot be rescued by the other two objectives."""
    rows = list(WORKED) + [_row("perfect", 0.99, 0.99, 0.99, plddt_vhh=55.0)]
    out = R.rank(rows)
    assert "perfect" not in {r["id"] for r in out["front"]}
    assert out["excluded_by"]["plddt_vhh"] == 1


def test_every_constraint_is_checked():
    for field, op, limit, _why in R.HARD_CONSTRAINTS:
        bad = limit - 1 if op in ("ge", "gt") else limit + 1
        rows = list(WORKED) + [_row("bad", 0.9, 0.9, 0.9, **{field: bad})]
        out = R.rank(rows)
        assert out["excluded_by"].get(field) == 1, f"{field} was not enforced"


def test_a_missing_measurement_is_not_a_constraint_failure():
    rows = list(WORKED)
    rows.append({k: v for k, v in _row("partial", 0.9, 0.9, 0.9).items()
                 if k != "plddt_vhh"})
    out = R.rank(rows)
    assert out["n_excluded"] == 0


def test_all_three_objectives_are_required():
    rows = [dict(r) for r in WORKED]
    del rows[0]["cleavability"]
    out = R.rank(rows)
    assert rows[0]["id"] not in {r["id"] for r in out["front"]}
    assert out["n_skipped"] == 1


def test_normalisation_puts_the_objectives_on_one_scale():
    scale = R.normalise(WORKED)
    assert set(scale) == set(R.OBJECTIVES)
    for lo, span in scale.values():
        assert span > 0


def test_the_uncertainty_penalty_lowers_a_score():
    scale = R.normalise(WORKED)
    plain = R.robust_score(WORKED[1], scale)
    docked = R.robust_score(WORKED[1], scale, penalty=1.0,
                            uncertainty={"feasibility": 0.5})
    assert docked < plain


def test_diversity_takes_one_of_each_design_first():
    rows = [_row("a1", 0.9, 0.9, 0.9), _row("a2", 0.8, 0.8, 0.8),
            _row("b1", 0.7, 0.7, 0.7)]
    for r in rows[:2]:
        r["mask"] = "same"
    rows[2]["mask"] = "other"
    chosen = R.select_diverse(rows, n=2)
    assert {r["mask"] for r in chosen} == {"same", "other"}


def test_the_ordering_rule_is_recorded_in_the_output():
    out = R.rank(list(WORKED))
    assert "maximin" in out["front_order"]
    assert "calibration" in out["normalisation"]
    assert out["constraints"], "the constraints applied must be reported"
