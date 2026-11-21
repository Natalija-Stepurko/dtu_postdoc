#!/usr/bin/env python
"""Stage 5. Draw every chart."""
import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[1] / "src"))
from nbswitch import figures
for path in figures.build_all():
    print(f"  {path.name}")
