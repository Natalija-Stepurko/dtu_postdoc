"""Equations typeset as inline SVG.

The page is read on machines whose maths fonts we cannot predict, and a
browser with no maths font turns a stacked fraction into a slash. The glyph
outlines are therefore embedded in the page at build time: the equations look
the same for every reader and nothing is fetched when the page loads.

Each equation is a separate SVG, and matplotlib names its glyph definitions
after the font slot they came from (``Cmmi10-72``). Inlined side by side in one
document those names collide, so every id is prefixed per equation before the
markup is handed to the page.
"""

from __future__ import annotations

import re

#: Rendered at this point size, then scaled to the reader's text size.
_BASE_PT = 10.0

#: Size of the equation's base glyphs, in rem. Computer Modern has a small
#: x-height next to the page's sans face, so matching the nominal size would
#: leave the maths looking shrunken; it is set a little larger instead.
#: In rem rather than em so the panel's own font-size does not scale it.
_SCALE = 1.28

_DROP = (
    # a <style> block inside an inlined SVG is document-wide in HTML, so the
    # "*" selector matplotlib writes would restyle every stroke on the page
    re.compile(r"<style[^>]*>.*?</style>", re.S),
    re.compile(r"<metadata>.*?</metadata>", re.S),  # holds a build timestamp
    re.compile(r"<\?xml[^>]*\?>\s*", re.S),
    re.compile(r"<!DOCTYPE[^>]*>\s*", re.S),
    re.compile(r"<!--.*?-->", re.S),
)


def _plt():
    import logging
    import matplotlib

    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "mathtext.fontset": "cm",   # Computer Modern, the usual maths face
        "svg.fonttype": "path",     # embed outlines, do not trust the reader
        "figure.facecolor": "none",
        "savefig.facecolor": "none",
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.015,
    })
    return plt


def equation(latex: str, label: str, key: str) -> str:
    """Typeset ``latex`` and return inline SVG markup.

    ``latex`` is a mathtext expression without the surrounding ``$``.
    ``label`` is read aloud in place of the image. ``key`` namespaces the
    glyph ids so several equations can share one page.
    """
    import io

    plt = _plt()
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, f"${latex}$", fontsize=_BASE_PT)
    buf = io.StringIO()
    fig.savefig(buf, format="svg", transparent=True)
    plt.close(fig)
    return _clean(buf.getvalue(), label, key)


def _clean(svg: str, label: str, key: str) -> str:
    for pat in _DROP:
        svg = pat.sub("", svg)

    # namespace every id, and the <use> references that point at them
    ids = sorted(set(re.findall(r'id="([^"]+)"', svg)), key=len, reverse=True)
    for i in ids:
        svg = svg.replace(f'id="{i}"', f'id="{key}-{i}"')
        svg = svg.replace(f'xlink:href="#{i}"', f'xlink:href="#{key}-{i}"')

    m = re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg)
    if not m:  # pragma: no cover - matplotlib always writes one
        raise ValueError("no viewBox in the rendered equation")
    h_rem = float(m.group(2)) / _BASE_PT * _SCALE

    # the glyphs carry no fill of their own, so the svg's own fill reaches
    # them and the equation takes the colour of the text around it
    svg = re.sub(
        r"<svg ",
        f'<svg role="img" aria-label="{label}" fill="currentColor" '
        f'style="height:{h_rem:.3f}rem;width:auto" ',
        svg, count=1,
    )
    svg = re.sub(r'\swidth="[\d.]+pt"', "", svg, count=1)
    svg = re.sub(r'\sheight="[\d.]+pt"', "", svg, count=1)
    return re.sub(r"\n\s*\n", "\n", svg).strip()
