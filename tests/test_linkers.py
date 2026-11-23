

def _gen_facts():
    """The schematic's inputs, derived the way the page derives them."""
    from nbswitch import linkers as L
    prof = L.conditioning_profile().iloc[0]
    neutral = [p for p in L.PROTEASES
               if p not in L.TARGET_PROTEASES and p not in L.OFF_TARGET_PROTEASES]
    return {"n_generated": 714, "n_kept": 4, "n_controls": 2,
            "n_target": len(L.TARGET_PROTEASES),
            "n_off": len(L.OFF_TARGET_PROTEASES),
            "n_neutral": len(neutral),
            "target_names": list(L.TARGET_PROTEASES),
            "off_names": list(L.OFF_TARGET_PROTEASES),
            "neutral_names": neutral,
            "target_z": float(prof[L.TARGET_PROTEASES[0]]),
            "off_z": float(prof[L.OFF_TARGET_PROTEASES[0]])}


def _schematic():
    from nbswitch import report as R, linkers as L
    return R.linker_schematic(_gen_facts(), L.LENGTH_SERIES, 10, L.REFERENCE_LENGTH)


def test_the_schematic_shows_three_protease_groups_not_two():
    """It said "15 others down". 15 is everything that is not a target; the
    profile pushes 9 down and leaves 6 at zero, which is no preference rather
    than resistance. Three groups, and the diagram has to show three."""
    from nbswitch import linkers as L
    facts, svg = _gen_facts(), _schematic()
    for n, want in ((facts["n_target"], "High cleavage"),
                    (facts["n_off"], "Low cleavage"),
                    (facts["n_neutral"], "Unspecified")):
        assert f"{n} MMPs" in svg, f"the group of {n} is missing"
        assert want in svg, f"'{want}' is missing"
    n_not_target = len(L.PROTEASES) - len(L.TARGET_PROTEASES)
    assert facts["n_neutral"], "a three-group profile needs a third group"
    assert f"{n_not_target} others down" not in svg
    assert f"{n_not_target} MMPs" not in svg, "the two down-weighted groups are merged"


def test_the_schematic_shows_the_z_scores_the_profile_carries():
    """The three groups differ by the number asked for, so show that number."""
    facts, svg = _gen_facts(), _schematic()
    assert f"z = {facts['target_z']:+.0f}" in svg
    assert f"z = {facts['off_z']:+.0f}" in svg
    assert "z = 0" in svg, "the unconstrained group must read as zero, not +0"


def test_the_schematic_names_the_proteases_in_each_group():
    from nbswitch import report as R
    facts, svg = _gen_facts(), _schematic()
    for names in (facts["target_names"], facts["off_names"], facts["neutral_names"]):
        # the names wrap across lines, so check the first of each run
        first = R._mmp_names(names).split(",")[0]
        assert first in svg, f"{first} missing from the diagram"


def test_the_schematic_counts_match_the_run():
    facts, svg = _gen_facts(), _schematic()
    assert f"{facts['n_generated']:,}" in svg
    assert f"{facts['n_kept']} highest-scoring" in svg
    assert f"{facts['n_controls']} controls" in svg


def test_the_schematic_reports_the_real_flank_range():
    """The flanks shown must be what flank_of_length would build."""
    from nbswitch import linkers as L
    svg = _schematic()
    lo, hi = min(L.LENGTH_SERIES), max(L.LENGTH_SERIES)
    assert f"{(lo - 10) // 2} to {(hi - 10) // 2} residues a side" in svg
    assert f"linkers of {lo} to {hi} residues" in svg
    assert f"{L.REFERENCE_LENGTH} residues, shown" in svg


def test_the_profile_pushes_only_the_two_named_groups():
    """Everything else stays at zero, so the selectivity contrast is a
    comparison between two stated sets and not between a set and the rest."""
    from nbswitch import linkers as L
    row = L.conditioning_profile(target_z=2.0, off_target_z=-1.0).iloc[0]
    assert set(row.index) == set(L.PROTEASES)
    assert [row[p] for p in L.TARGET_PROTEASES] == [2.0] * len(L.TARGET_PROTEASES)
    assert [row[p] for p in L.OFF_TARGET_PROTEASES] == [-1.0] * len(L.OFF_TARGET_PROTEASES)
    rest = [p for p in L.PROTEASES
            if p not in L.TARGET_PROTEASES and p not in L.OFF_TARGET_PROTEASES]
    assert rest, "a three-group profile needs a third group"
    assert all(row[p] == 0.0 for p in rest)
