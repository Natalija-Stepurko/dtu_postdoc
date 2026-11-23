#!/usr/bin/env python
"""Stage 7b. Score every folded construct's own mask sequence with ProteinMPNN.

Stage 07 recorded a ProteinMPNN score for the 96 redesigns and for the 12
parents it redesigned, which is 108 of the 348 folded constructs. Asking
whether that score predicts masking on 108 of them answers the question on
less than a third of the library.

Two things are fixed here.

**Coverage.** Every construct with a mask is scored, library and redesigns
alike, over its own mask positions with the nanobody and linker held fixed -
the same span stage 07 designed into, so the number means the same thing.

**Noise.** The score depends on the order ProteinMPNN decodes the sequence in,
and a single draw carries a standard deviation of about 0.03 to 0.05. Stage 07
reported single draws. Each construct is scored over ``N_ORDERS`` decoding
orders here and the mean is reported with its spread, so a weak correlation
cannot be an artefact of comparing one noisy draw against another.

Scoring is cheap next to folding: one model load covers the whole library.
"""
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from nbswitch.design import PROTEINMPNN, available, _n_residues
from nbswitch.config import STRUCTURES, load_json, save_json

#: How many decoding orders to average over. 20 puts the standard error on the
#: mean near 0.008, an order of magnitude below the spread across constructs.
N_ORDERS = 20
CHAIN = "A"

if not available():
    raise SystemExit(f"ProteinMPNN not found at {PROTEINMPNN}")

lib = load_json("03_constructs.json")
targets = {}
for c in lib["constructs"]:
    if c.get("mask") and c.get("mask_range"):
        targets[c["id"]] = (c["mask_range"], "library", c["mask"])
try:
    for c in load_json("07_redesign.json")["constructs"]:
        targets[c["id"]] = (c["mask_range"], "redesign", c["parent"])
except FileNotFoundError:
    pass

missing = [cid for cid in targets if not (STRUCTURES / f"{cid}.pdb").exists()]
for cid in missing:
    targets.pop(cid)
print(f"{len(targets)} constructs to score"
      + (f", {len(missing)} skipped with no folded structure" if missing else ""))

t0 = time.time()
with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    pdbs = tmp / "pdbs"
    pdbs.mkdir()
    fixed = {}
    for cid, (rng, _kind, _parent) in targets.items():
        src = STRUCTURES / f"{cid}.pdb"
        (pdbs / src.name).write_text(src.read_text())
        lo, hi = rng
        n_res = _n_residues(src, CHAIN)
        # ProteinMPNN's fixed-position dict lists the positions to KEEP, so
        # the loss is taken over the mask alone
        fixed[cid] = {CHAIN: [i for i in range(1, n_res + 1) if not (lo <= i <= hi)]}
    (tmp / "fixed.jsonl").write_text(json.dumps(fixed) + "\n")

    subprocess.run([sys.executable,
                    str(PROTEINMPNN / "helper_scripts/parse_multiple_chains.py"),
                    "--input_path", str(pdbs),
                    "--output_path", str(tmp / "parsed.jsonl")],
                   check=True, capture_output=True, text=True)
    subprocess.run([sys.executable, str(PROTEINMPNN / "protein_mpnn_run.py"),
                    "--jsonl_path", str(tmp / "parsed.jsonl"),
                    "--fixed_positions_jsonl", str(tmp / "fixed.jsonl"),
                    "--out_folder", str(tmp / "out"),
                    "--score_only", "1",
                    "--num_seq_per_target", str(N_ORDERS),
                    "--seed", "42", "--batch_size", "1"],
                   check=True, capture_output=True, text=True)

    scored = {}
    for path in sorted((tmp / "out" / "score_only").glob("*.npz")):
        cid = path.stem[:-4] if path.stem.endswith("_pdb") else path.stem
        if cid not in targets:
            continue
        draws = np.atleast_1d(np.load(path)["score"]).ravel().astype(float)
        rng, kind, parent = targets[cid]
        scored[cid] = {
            "mask_score": float(draws.mean()),
            "mask_score_sd": float(draws.std(ddof=1)) if draws.size > 1 else 0.0,
            "n_orders": int(draws.size),
            "mask_residues": rng[1] - rng[0] + 1,
            "source": kind, "mask": parent}

elapsed = time.time() - t0
if not scored:
    raise SystemExit("ProteinMPNN produced no scores")

save_json("07b_library_scores.json", {
    "method": f"ProteinMPNN --score_only over the mask positions, mean of "
              f"{N_ORDERS} decoding orders, nanobody and linker fixed",
    "n_orders": N_ORDERS, "n_scored": len(scored),
    "note": "the same span stage 07 designs into, so library constructs and "
            "redesigns are measured on one scale. Reported as a mean with its "
            "spread because a single decoding order carries a standard "
            "deviation of roughly 0.03 to 0.05.",
    "seconds": round(elapsed, 1),
    "scores": scored})

sds = [v["mask_score_sd"] for v in scored.values()]
vals = [v["mask_score"] for v in scored.values()]
print(f"scored {len(scored)} constructs in {elapsed:.0f} s")
print(f"  mask score {min(vals):.2f} to {max(vals):.2f}")
print(f"  decoding-order sd {np.mean(sds):.3f} mean, {max(sds):.3f} worst")
by_source = {}
for v in scored.values():
    by_source.setdefault(v["source"], []).append(v["mask_score"])
for src, xs in sorted(by_source.items()):
    print(f"  {src:9s} n = {len(xs):3d}  mean {np.mean(xs):.3f}")
