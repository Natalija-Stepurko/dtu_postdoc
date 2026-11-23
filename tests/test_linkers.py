

def test_the_schematic_reads_the_real_length_series():
    """The diagram must not drift from the pipeline it describes."""
    from nbswitch import report as R, linkers as L
    svg = R.linker_schematic({"n_generated": 714, "n_kept": 4, "n_off": 15},
                             L.LENGTH_SERIES, 10, L.REFERENCE_LENGTH)
    for n in L.LENGTH_SERIES:
        assert f"n = {n}" in svg, f"length {n} missing from the diagram"
    assert str(L.REFERENCE_LENGTH) in svg
    # the flanks shown must be what flank_of_length would build
    for n in L.LENGTH_SERIES:
        assert f"{(n - 10) // 2} Gly/Ser each side" in svg
