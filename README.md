# nbswitch — a protease-activated nanobody mask

A concept demonstration: a nanobody that stays inert in circulation and binds its antigen
only after a protease cuts its tether. The regulatory function lives in a separately
designed mask and linker, so an established neutralizing nanobody carries over unchanged.

**Results page → https://natalija-stepurko.github.io/dtu_postdoc/**

The construct is a single polypeptide — mask `M`, cleavable linker `L`, nanobody `N`:

```
N–L–M   ──protease──▶   N + M        then        N + T ⇌ NT
```

Everything is computed from one experimental complex,
[PDB 9GCN](https://www.rcsb.org/structure/9GCN) (VHH TPL1158_01_C09 bound to
α-cobratoxin). Each stage writes a JSON file, and the figures and the page are generated
from those files, so no number is entered by hand.

## What it shows

**The two termini are not interchangeable.** The nanobody's N-terminus sits 11.8 Å from
the paratope centroid and the C-terminus 42.0 Å. Feeding those measured distances through
a worm-like-chain tether model bounds the mask affinity from both sides — tight enough to
block while tethered, loose enough to leave after cleavage. An N-terminal fusion with a
5-residue linker leaves about **3.9 decades** of usable mask affinity; a C-terminal fusion
leaves 2.2, and below 10 residues it has **no window at all**.

**Tethering is worth more than affinity.** At 5 residues the masking ceiling is ~750 µM,
roughly 7600× above the release floor. A micromolar mask still blocks, so masks do not
need affinity maturation — and tight masks are actively harmful, because they fail to
leave once the tether is cut.

**The paratope is almost entirely CDR3** (11 of 16 residues; CDR2 contributes none), and
the N-terminal residue is itself part of it.

## Running it

```bash
uv venv && uv pip install -e ".[folding]"
for s in pipeline/0*.py; do python "$s"; done
pytest
```

Stages 1–3, 5 and 6 take under a minute on a laptop. Stage 4 folds 26 constructs with
ESMFold; on 8 CPU cores that is about 20 s each, and it is resumable, so it can be
interrupted and restarted. A GPU only makes it faster.

| stage | what it does | needs |
|---|---|---|
| `01_reference_structure.py` | interface, paratope, CDRs, terminal geometry | CPU |
| `02_tether_model.py` | effective concentration, masking ceiling, release floor, design window | CPU |
| `03_construct_library.py` | mask candidates from the structure, with controls | CPU |
| `04_fold_constructs.py` | ESMFold on every construct, paratope occlusion | CPU or GPU |
| `05_figures.py` | every chart on the page | CPU |
| `06_report.py` | the page itself | CPU |

`notebooks/esmfold_masked_vhh.ipynb` runs stage 4 on a free Colab T4 instead
([open in Colab](https://colab.research.google.com/github/Natalija-Stepurko/dtu_postdoc/blob/main/notebooks/esmfold_masked_vhh.ipynb)).

## Layout

```
src/nbswitch/      config · structure · tether · constructs · folding · occlusion · figures · report
pipeline/          one script per stage, each writing JSON to results/
tests/             17 tests, including the residue-keying trap below
notebooks/         the Colab version of the folding stage
results/           stage outputs, figures and predicted structures
docs/              the generated page
```

## A trap worth knowing about

9GCN is Kabat-numbered, so residues carry insertion codes: 122 residues share only 112
distinct `res_id` values. Keying residues on `res_id` alone silently merges them and
truncates CDR3 from 17 residues to 10 — which would quietly corrupt the paratope, the
geometry and everything downstream. Every routine keys on `(res_id, ins_code)`, and
`tests/test_structure.py` pins it.

## What is not computed

Mask backbone generation (RFdiffusion, BindCraft), complex prediction (AlphaFold 3,
Boltz-2, Chai-1) and conformational ensembles (BioEmu-2) are not run here.

Trigger selectivity is a data problem before it is a modelling one: published profiling of
snake venom metalloproteinases shows a P1′ leucine preference shared with host MMPs, so
peptide-cleavage measurements on venom proteases against a host negative panel would have
to be generated before a linker could be designed to discriminate between them.

## Licence

MIT for the code. Structural data belongs to the RCSB PDB and its depositors.
