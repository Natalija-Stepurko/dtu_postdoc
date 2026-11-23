"""Folding a module must not break the markup or strand a viewer.

A module runs from its own heading to the next one or to the end of its
section, and the wrapper only holds. The risk is at the edges: a span that
is not balanced leaves the page malformed, and a 3Dmol viewer created inside
a closed <details> has no size, draws nothing, and stays blank when the
section is opened.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import re

from nbswitch import report as R


VOID = {"img", "br", "hr", "meta", "link", "input", "col", "source"}


def _balanced(fragment):
    stack = []
    for m in re.finditer(r"<(/?)([a-zA-Z][a-zA-Z0-9]*)([^>]*)>", fragment):
        close, tag, rest = m.group(1), m.group(2).lower(), m.group(3)
        if tag in VOID or rest.rstrip().endswith("/"):
            continue
        if not close:
            stack.append(tag)
        elif stack and stack[-1] == tag:
            stack.pop()
        else:
            return False
    return not stack


PAGE = ('<section class="sec"><p>intro</p>'
        '<h3 class="mod"><span>Module 01</span>First</h3>'
        '<div class="stack"><p>one</p></div>'
        '<details class="tech"><summary>inner</summary><p>x</p></details>'
        '<h3 class="mod"><span>Module 02</span>Second</h3>'
        '<figure><div>two</div></figure>'
        '</section><section class="sec"><p>after</p></section>')


def test_every_module_heading_becomes_a_summary():
    out = R.fold_modules(PAGE)
    assert out.count('<details class="modfold" open>') == 2
    assert out.count("</details>") == 3          # two modules plus the inner one
    for title in ("First", "Second"):
        assert f"<summary><h3 class=\"mod\">" in out
        assert title in out


def test_a_module_stops_at_the_next_one_and_at_the_section_end():
    out = R.fold_modules(PAGE)
    first = out[out.index("Module 01"):out.index("Module 02")]
    assert "one" in first and "two" not in first
    tail = out[out.index("Module 02"):]
    assert "two" in tail
    assert "after" not in tail[:tail.index("</section>")]


def test_the_wrapped_page_is_still_balanced():
    assert _balanced(PAGE), "the fixture itself must be balanced"
    assert _balanced(R.fold_modules(PAGE))


def test_a_page_with_no_modules_is_returned_unchanged():
    plain = '<section class="sec"><p>nothing to fold</p></section>'
    assert R.fold_modules(plain) == plain


def test_modules_open_by_default():
    """Folding is for putting a section away once read, not for hiding the
    page behind its headings."""
    assert "<details class=\"modfold\" open>" in R.fold_modules(PAGE)


def test_every_viewer_goes_through_the_registry():
    """A viewer created directly cannot be resized when its module opens, so
    it renders blank for anyone who collapsed and reopened the section."""
    js = R.VIEWER_JS
    assert js.count("$3Dmol.createViewer(") == 1, (
        "only makeViewer may call createViewer directly")
    assert "VIEWERS.push" in js
    assert js.count("makeViewer(") >= 5, "a viewer is bypassing the registry"


def test_opening_a_module_resizes_the_viewers_inside_it():
    js = R.VIEWER_JS
    assert "details.modfold" in js and "toggle" in js
    assert "wakeViewers" in js
    assert ".resize()" in js and ".render()" in js
