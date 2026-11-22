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
[PDB 9RIU](https://www.rcsb.org/structure/9RIU) (V<sub>H</sub>H5 a-sNTx bound to short
neurotoxin 1, P01426, from *Naja pallida*; Ahmadi et al., *Nature* 2025). Each stage writes a
JSON file, and the figures and the page are generated from those files, so no number is entered
by hand.

[PDB 9GCN](https://www.rcsb.org/structure/9GCN) (VHH TPL1158_01_C09 bound to α-cobratoxin) was
the first target and is kept only as a set-aside case: its own N-terminus sits in its paratope,
so masking it is trivially easy and the measure says nothing. The target is selected by
`NBSWITCH_TARGET`, and results are written per target under `results/<id>/`.

## What the pipeline does

Six modules. One estimates what is achievable, two design independently, three score,
and a ranking selects; a final step verifies what survived.

**01 — Usable mask affinity.** The interface is measured on the deposited complex, and
the distance from each terminus to the paratope goes into a worm-like-chain tether model.
Because the mask is tethered it acts at a high local concentration, so a weak mask blocks;
after cleavage that advantage is gone and the mask has to leave. The two requirements
bound mask affinity from opposite sides, which gives a usable window rather than a target
to maximise. The N-terminus is the attachment point: it is 21.2 Å from the paratope
centroid against 34.2 Å for the C-terminus.

**02a — Mask design.** The mask starts from what the antigen presents at the paratope:
every antigen residue with a heavy atom within 3.5 Å of the nanobody, grouped into contact
faces. That face seeds a space rather than being the candidate — window variants carrying
more or less flanking antigen, a variant with the disulfide removed, substitution variants
holding the highest-burying positions fixed, a mask covering both contact faces, and
matched scrambled, truncated and inert controls.

**02b — Cleavable linker design.** Cleavability is a constraint on generation, not a
filter afterwards: a Gly/Ser repeat is chosen in protein engineering precisely because
proteases do not cut it, so no score applied later could have found a cleavable one.
CleaveNet generates substrates against a target profile and they are scored as a contrast
between bite-site and circulating proteases, because a linker cut by everything would
unmask in transit. Length is a separate axis, since the substrate is fixed at ten residues
and the flanks are free.

**03, 04, 05 — Three scores.** Can the construct be made; does the mask block the antigen
and release when cut; is the linker both recognised and reachable. Each is computed for
every construct from one shared structure prediction.

**06 — Ranking.** The nondominated set over the three scores, with no weights, so no score
compensates for another. The front is then verified by predicting the whole assembly with
its antigen, which is the only step that asks directly whether the antigen can still bind.

## Running it

```bash
uv venv && uv pip install -e ".[folding]"
NBSWITCH_TARGET=9RIU PYTHON=python bash pipeline/run_all.sh
pytest
```

`run_all.sh` runs all twenty stages in order. Stage numbers are execution order, which is
why verification is 19 and not 11 — it consumes the ranking from 13.

Most stages take seconds. The cost is folding and structure prediction: about fifty
seconds a construct for ESMFold on eight CPU cores, and two to three minutes for each
Boltz-2 prediction. The whole pipeline is a few hours on CPU; a GPU only makes it faster.
Folding is resumable and keyed on the sequence in each file, so an interrupted run
restarts without losing work and cannot serve a model of a construct that has changed.

The library is a design matrix rather than a full cross product: every mask meets every
linker at a reference length, and a chosen few masks meet a chosen few linkers at the
other lengths. Stage 03 prints the size and the estimated fold time before building
anything and refuses to exceed `MAX_CONSTRUCTS`, so growing the search space is a priced
decision.

| stage | what it does | cost |
|---|---|---|
| `01_reference_structure.py` | interface, paratope, CDRs, terminal geometry | seconds |
| `02_tether_model.py` | effective concentration, masking ceiling, release floor, window | seconds |
| `02b_linkers.py` | generate cleavable substrates, score them, build the length series | minutes |
| `03a_mask_design.py` | contact faces, and the mask space seeded from them | seconds |
| `03_construct_library.py` | assemble the design matrix | seconds |
| `04_fold_constructs.py` | ESMFold on every construct | ~50 s each |
| `16_linker_length_sweep.py` | reach against linker length, substrate held constant | ~50 s each |
| `07_redesign_masks.py` | ProteinMPNN on the folded backbones, mask positions only | minutes |
| `08_fold_redesigns.py` | ESMFold on the redesigns | ~50 s each |
| `09_complex_prediction.py` | Boltz-2 without alignments — the control that shows why they are needed | ~3 min each |
| `10_trigger_accessibility.py` | extension, exposure and catalytic-domain clearance | seconds |
| `12_feasibility.py` | pLDDT, fold preservation, exposed hydrophobic surface, measured disulfides | seconds |
| `13_rank.py` | the nondominated set | seconds |
| `14_ensemble.py` | every design superimposed on one backbone, for the viewers | seconds |
| `15_interface_controls.py` | the same measures on independent nanobody–toxin complexes | seconds |
| `17_mask_binding.py` | does the mask bind untethered, with the antigen as positive control | ~3 min each |
| `18_cleavage_evidence.py` | the linkers against the measured substrate library, without a neural net | seconds |
| `19_verify_front.py` | the whole assembly with its antigen, for the ranked front | ~3 min each |
| `05_figures.py` | every figure, from stage outputs | seconds |
| `06_report.py` | `docs/index.html`, from stage outputs | seconds |

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
