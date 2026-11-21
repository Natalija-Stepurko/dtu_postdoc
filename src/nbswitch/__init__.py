"""nbswitch - structure-grounded design of protease-activated nanobody masks.

A concept demonstration. The construct is a single chain, mask-linker-nanobody,
that is intended to be inert while the tether holds and active once a protease
cuts it. The package computes the geometry and thermodynamics of that tether
from an experimental nanobody-antigen complex, builds a construct library, and
folds it.
"""
__version__ = "0.2.0"
from . import config, structure, tether, constructs, occlusion  # noqa: F401
