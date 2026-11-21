# Computational design of a universal protease-activated nanobody switch

A structure-grounded pipeline for designing a nanobody that stays masked in systemic
circulation and binds venom toxin only after local proteolysis. The regulatory function
lives in a separately designed mask and linker, so an established neutralizing VHH is
carried over unchanged.

**Results page: https://natalija-stepurko.github.io/dtu_postdoc/**

> **Concept demonstration.** A worked example showing the pipeline running end to end on a
> real structure. Not a therapeutic design campaign — the mask candidates are deliberately
> simple and no candidate here is proposed for synthesis.

The construct is a single chain — mask `M`, protease-cleavable linker `L`, fixed
nanobody `N`:

```
N–L–M  --protease-->  N + M      then      N + T ⇌ NT
```

## What is computed here

Everything on the results page is produced by the code in `src/`, run against
[PDB 9GCN](https://www.rcsb.org/structure/9GCN) (VHH TPL1158_01_C09 bound to
α-cobratoxin). Each stage writes a JSON file; figures and the page are generated from
those files.

| stage | script | output |
|---|---|---|
| 1a | `s01_fetch.py` | structure provenance, chains, sequences |
| 1b | `s02_interface.py` | paratope and epitope from buried surface area, CDRs, terminal geometry |
| 2 | `s03_tether.py` | effective concentration, masking ceiling, release floor, design window |
| 3 | `s04_figures.py` | figures |
| 4 | `s05_constructs.py` | mask candidates derived from the structure, construct library |
| 5 | `s06_page.py` | the results page |

## Headline results

- The paratope is concentrated in CDR3 (11 of 16 residues; CDR2 contributes none), and
  the N-terminal residue is itself part of the paratope.
- The two termini are not interchangeable: the N-terminus sits 11.8 Å from the paratope
  centroid, the C-terminus 42.0 Å.
- An N-terminal fusion with a 5-residue linker gives ~3.9 decades of usable mask
  affinity. A C-terminal fusion gives ~2.2, and below 10 residues it has **no window at
  all** — no mask affinity satisfies both masking and release.
- The tether is worth enough that a **micromolar** mask still blocks. Tight masks are
  actively harmful, because they fail to leave after cleavage.

## Running it

```bash
uv venv && uv pip install -r requirements.txt
for s in src/s0*.py; do python "$s"; done
```

Runs on CPU in under a minute; it fetches 9GCN from RCSB on first use.

## GPU step

`colab/esmfold_masked_vhh.ipynb` folds every construct with ESMFold and measures whether
the mask reaches and covers the paratope, contrasting the grafted mask against scrambled
and poly-glycine controls.
[Open in Colab](https://colab.research.google.com/github/Natalija-Stepurko/dtu_postdoc/blob/main/colab/esmfold_masked_vhh.ipynb).
The fusion is a single chain, so ESMFold needs no MSA and a ~160-residue construct folds
in seconds on a free T4.

## What is not computed

Mask backbone generation (RFdiffusion, BindCraft), complex prediction (AlphaFold 3,
Boltz-2, Chai-1) and ensembles (BioEmu-2) need a GPU and are not run here. The protease
trigger score is a data problem before it is a modelling one: published SVMP profiling
shows a P1′ leucine preference shared with host MMPs, so peptide-cleavage data on venom
proteases against a host negative panel would have to be generated first.

## Licence

MIT for the code. Structural data belongs to the RCSB PDB and its depositors.
