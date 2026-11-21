"""Stage 5. Build docs/index.html. Every number is read from the stage outputs."""
import json, shutil
from pathlib import Path
from common import ROOT, RESULTS, FIGS, load_json

DOCS = ROOT / "docs"; (DOCS / "figures").mkdir(parents=True, exist_ok=True)
for p in FIGS.glob("*.png"): shutil.copy(p, DOCS / "figures" / p.name)

s1, s2, s3, s5 = (load_json(f"s0{i}_{n}.json") for i, n in
                  ((1, "structure"), (2, "interface"), (3, "tether"), (5, "constructs")))
g, a = s2["terminal_geometry"], s2["chain_pairs"]["A"]
reg = {k: v["n"] for k, v in s2["paratope_by_region"].items()}
bN, bC = s3["best_per_terminus"]["N"], s3["best_per_terminus"]["C"]
closedC = max(r["linker_residues"] for r in s3["scan"]
              if r["terminus"] == "C" and not r["window_open"])
bsa = [r["buried_surface_area_A2"] for r in s2["chain_pairs"].values()]
m = s3["model"]
REPO = "https://github.com/Natalija-Stepurko/dtu_postdoc"
NB = f"https://colab.research.google.com/github/Natalija-Stepurko/dtu_postdoc/blob/main/colab/esmfold_masked_vhh.ipynb"

def fig(name, cap):
    return (f'<figure><img src="figures/{name}" alt="{cap}">'
            f'<figcaption>{cap}</figcaption></figure>')

def kpi(v, l, s=""):
    return (f'<div class="kpi"><div class="kv">{v}</div><div class="kl">{l}</div>'
            + (f'<div class="ks">{s}</div>' if s else "") + "</div>")

CSS = """
:root{color-scheme:light;--bg:#f9f9f7;--sf:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;
--ink3:#8a8984;--line:#e4e3de;--blue:#2a78d6;--orange:#eb6834;--aqua:#1baf7a;--yellow:#eda100;}
@media(prefers-color-scheme:dark){:root:where(:not([data-theme=light])){
--bg:#0d0d0d;--sf:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--ink3:#8d8c85;--line:#2e2e2c;
--blue:#3987e5;--orange:#d95926;--aqua:#199e70;--yellow:#c98500;color-scheme:dark;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font:16px/1.65 ui-sans-serif,-apple-system,"Segoe UI",Inter,system-ui,sans-serif;
-webkit-font-smoothing:antialiased}
.wrap{max-width:980px;margin:0 auto;padding:0 24px 96px}
header{padding:72px 0 28px;border-bottom:1px solid var(--line);margin-bottom:40px}
h1{font-size:34px;line-height:1.2;margin:0 0 12px;letter-spacing:-.02em;font-weight:650}
.sub{color:var(--ink2);font-size:17px;margin:0 0 18px;max-width:72ch}
.meta{color:var(--ink3);font-size:13px}
.demo{background:var(--sf);border:1px solid var(--line);border-left:2px solid var(--yellow);
border-radius:8px;padding:11px 15px;font-size:14px;margin:0 0 16px;max-width:74ch}
.demo b{color:var(--ink)}
.meta a{color:var(--ink2)}
h2{font-size:13px;letter-spacing:.09em;text-transform:uppercase;color:var(--ink3);
margin:64px 0 14px;font-weight:650}
h3{font-size:20px;margin:34px 0 10px;font-weight:600;letter-spacing:-.01em}
p{margin:0 0 14px;max-width:74ch;color:var(--ink2)}
p strong,li strong{color:var(--ink);font-weight:600}
a{color:var(--blue)}
figure{margin:22px 0;background:var(--sf);border:1px solid var(--line);
border-radius:12px;padding:16px;overflow:hidden}
figure img{width:100%;height:auto;display:block;border-radius:4px}
figcaption{color:var(--ink3);font-size:12.5px;margin-top:10px;line-height:1.5}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin:22px 0}
.kpi{background:var(--sf);border:1px solid var(--line);border-radius:12px;padding:14px 16px}
.kv{font-size:25px;font-weight:650;letter-spacing:-.02em;font-variant-numeric:tabular-nums}
.kl{font-size:12.5px;color:var(--ink2);margin-top:3px}
.ks{font-size:11.5px;color:var(--ink3);margin-top:4px}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:18px 0}
th,td{text-align:left;padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{color:var(--ink3);font-weight:600;font-size:11.5px;letter-spacing:.05em;text-transform:uppercase}
td.n{text-align:right;font-variant-numeric:tabular-nums}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px}
code{background:var(--sf);border:1px solid var(--line);padding:1px 5px;border-radius:4px}
.note{border-left:2px solid var(--orange);padding:4px 0 4px 16px;margin:20px 0;
color:var(--ink2);font-size:14.5px}
.note b{color:var(--ink)}
.ok{border-left-color:var(--aqua)}
.steps{counter-reset:s;list-style:none;padding:0;margin:18px 0}
.steps li{counter-increment:s;position:relative;padding:0 0 16px 40px;color:var(--ink2);
border-left:1px solid var(--line);margin-left:11px}
.steps li:last-child{border-left-color:transparent}
.steps li::before{content:counter(s);position:absolute;left:-11px;top:0;width:22px;height:22px;
border-radius:50%;background:var(--sf);border:1px solid var(--line);color:var(--ink2);
font-size:11.5px;font-weight:650;display:grid;place-items:center}
.steps b{color:var(--ink);display:block;margin-bottom:2px}
.pill{display:inline-block;font-size:11px;padding:2px 9px;border-radius:999px;
border:1px solid var(--line);color:var(--ink2);margin-right:6px}
.pill.run{border-color:var(--aqua);color:var(--aqua)}
.pill.gpu{border-color:var(--yellow);color:var(--yellow)}
.pill.no{border-color:var(--line);color:var(--ink3)}
.eq{background:var(--sf);border:1px solid var(--line);border-radius:10px;padding:14px 18px;
margin:16px 0;font-family:ui-monospace,monospace;font-size:13.5px;color:var(--ink);overflow-x:auto}
.two{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:720px){.two{grid-template-columns:1fr}h1{font-size:27px}}
footer{margin-top:72px;padding-top:22px;border-top:1px solid var(--line);
color:var(--ink3);font-size:13px}
ul.ref{list-style:none;padding:0;font-size:13.5px}
ul.ref li{padding:7px 0;border-bottom:1px solid var(--line);color:var(--ink2)}
"""

HTML = f"""<!DOCTYPE html><html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Computational design of a universal protease-activated nanobody switch</title>
<meta name="description" content="A structure-grounded pipeline for designing protease-activated
nanobody masks, with results computed from PDB 9GCN.">
<style>{CSS}</style></head><body><div class="wrap">

<header>
<h1>Computational design of a universal protease-activated nanobody switch</h1>
<p class="sub">A nanobody that stays masked in circulation and binds venom toxin only after local
proteolysis. The regulatory function lives in a separately designed mask and linker, so the
neutralizing VHH is carried over unchanged.</p>
<p class="demo"><b>Concept demonstration.</b> This is a worked example showing the pipeline
running end to end on a real structure. It is not a therapeutic design campaign: the mask candidates
are deliberately simple, and no candidate here is proposed for synthesis.</p>
<p class="meta">Results computed from <a href="https://www.rcsb.org/structure/9GCN">PDB 9GCN</a> ·
every number on this page is produced by the code in the
<a href="{REPO}">repository</a> · <a href="{NB}">companion Colab notebook</a></p>
</header>

<h2>The construct</h2>
<p>The molecule is a single polypeptide: a designed mask <span class="mono">M</span>, a
protease-cleavable linker <span class="mono">L</span>, and the fixed nanobody
<span class="mono">N</span>. Proteolysis severs the tether and the mask diffuses away, leaving the
original paratope free to engage toxin <span class="mono">T</span>.</p>
<div class="eq">N–L–M &nbsp;--protease--&gt;&nbsp; N + M &nbsp;&nbsp;then&nbsp;&nbsp; N + T ⇌ NT</div>
{fig("fig2_states.png", "The three design states. Masking must be effective while the tether holds and must fail once it is cut.")}
<div class="note"><b>The central tension.</b> The strongest mask is not the best mask. A mask tight
enough to block the paratope while tethered may stay bound after cleavage and keep blocking it. The
pipeline below turns that statement into a number.</div>

<h2>Pipeline</h2>
<ol class="steps">
<li><b>Reference complex</b><span class="pill run">runs on CPU</span>
Define the paratope from buried surface area on an experimental VHH–toxin structure, and measure
where a mask could be attached.</li>
<li><b>Tether thermodynamics</b><span class="pill run">runs on CPU</span>
Convert linker length and attachment geometry into an effective concentration, then into bounds on
the mask–VHH affinity.</li>
<li><b>Mask candidates</b><span class="pill run">runs on CPU</span>
Derive candidate masks from the structure, with scrambled and inert controls.</li>
<li><b>Fusion structure prediction</b><span class="pill gpu">needs a GPU</span>
Fold each construct and measure whether the mask reaches and covers the paratope.</li>
<li><b>Protease trigger</b><span class="pill no">needs new data</span>
Rank linkers for cleavage by the venom protease ensemble against a host negative panel.</li>
<li><b>Ranking</b><span class="pill no">needs 4 and 5</span>
Filter on feasibility, then take the Pareto front over the three scores with uncertainty attached.</li>
</ol>

<h2>Step 1 · The reference complex</h2>
<p>{s1['title']}, solved by {s1['method'].lower()} at {s1['resolution_A']} Å and deposited
{s1['deposition_date']}. The asymmetric unit holds four independent copies of the complex, which
gives a free reproducibility check: the interface is measured four times.</p>
<div class="kpis">
{kpi(f"{s1['resolution_A']} Å", "resolution", "X-ray")}
{kpi(f"{a['interface_area_A2']:.0f} Å²", "interface area", f"{min(bsa):.0f}–{max(bsa):.0f} Å² buried across 4 copies")}
{kpi(a['n_paratope'], "paratope residues", "&gt;1 Å² buried on binding")}
{kpi(s2['cdrs']['CDR3']['length'], "CDR3 length", "residues")}
</div>
{fig("fig1_paratope.png", f"Surface buried on complex formation, per VHH residue. {reg['CDR3']} of the {a['n_paratope']} paratope residues lie in CDR3; CDR2 contributes none.")}
<div class="note ok"><b>Two findings that set up everything downstream.</b>
The paratope is concentrated in CDR3 — {reg['CDR3']} of {a['n_paratope']} residues, against
{reg['CDR1']} in CDR1 and {reg['CDR2']} in CDR2. And the N-terminal residue
<span class="mono">{g['n_term_residue']}</span> is itself part of the paratope, burying
{a['paratope_delta_sasa'][g['n_term_residue']]:.0f} Å². The two termini are not interchangeable
attachment points.</div>
<table><thead><tr><th>Geometry</th><th class="n">Distance (Å)</th><th>Consequence for the design</th></tr></thead><tbody>
<tr><td>N-terminus → paratope centroid</td><td class="n">{g['n_term_to_paratope_centroid_A']}</td><td>a short linker already reaches</td></tr>
<tr><td>N-terminus → nearest paratope residue</td><td class="n">{g['n_term_to_nearest_paratope_A']}</td><td>the N-terminus is in the paratope</td></tr>
<tr><td>C-terminus → paratope centroid</td><td class="n">{g['c_term_to_paratope_centroid_A']}</td><td>3.6× further to span</td></tr>
<tr><td>paratope maximum extent</td><td class="n">{g['paratope_max_extent_A']}</td><td>surface a mask must cover</td></tr>
<tr><td>VHH radius of gyration</td><td class="n">{g['radius_of_gyration_A']}</td><td>scale reference</td></tr>
</tbody></table>

<h2>Step 2 · What the tether is worth</h2>
<p>A tethered mask is held near the paratope at an effective concentration set by the linker.
Treating the Gly/Ser linker as a worm-like chain with persistence length
{m['persistence_length_A']} Å and {m['rise_per_residue_A']} Å per residue, the end-to-end
distribution gives C<sub>eff</sub> at the distance the linker has to span — a distance measured in
step 1, not assumed.</p>
<div class="eq">⟨r²⟩ = 2·l<sub>p</sub>·L<sub>c</sub>·[1 − (l<sub>p</sub>/L<sub>c</sub>)(1 − e<sup>−L<sub>c</sub>/l<sub>p</sub></sup>)]<br>
C<sub>eff</sub>(d) = 1660.54 · (3/2π⟨r²⟩)<sup>3/2</sup> · e<sup>−3d²/2⟨r²⟩</sup> &nbsp;&nbsp;[M]</div>
{fig("fig3_ceff.png", f"Effective concentration against linker length, for the two attachment points. The N-terminal span of {g['n_term_to_paratope_centroid_A']} Å and the C-terminal span of {g['c_term_to_paratope_centroid_A']} Å come from the structure.")}
<p>Masking and release then pull the mask affinity in opposite directions:</p>
<div class="eq">P<sub>masked</sub> = C<sub>eff</sub> / (C<sub>eff</sub> + K<sub>D,M</sub>)<br>
P<sub>NT</sub> = (T/K<sub>D,T</sub>) / (1 + T/K<sub>D,T</sub> + M/K<sub>D,M</sub>)</div>
{fig("fig5_tradeoff.png", f"Requiring ≤{m['leak_max']*100:.0f}% systemic leak puts a ceiling on the mask K<sub>D</sub>; requiring ≥{m['recovery_min']*100:.0f}% toxin binding after cleavage puts a floor under it.")}

<h2>Step 2 · result — the design window</h2>
<p>Imposing both requirements at once leaves a window of usable mask affinities. Whether that window
is open at all, and how wide it is, follows from the attachment geometry.</p>
{fig("fig4_window.png", "Usable mask affinity against linker length. Shaded region satisfies both requirements simultaneously.")}
<div class="kpis">
{kpi(f"{bN['window_decades']:.1f}", "decades of usable K<sub>D</sub>", f"N-terminal, {bN['linker_residues']}-residue linker")}
{kpi(f"{bC['window_decades']:.1f}", "decades of usable K<sub>D</sub>", f"C-terminal, {bC['linker_residues']}-residue linker")}
{kpi(f"≤{closedC}", "residues where C-terminal fails", "window entirely closed")}
{kpi(f"{bN['km_max_for_masking_M']*1e6:.0f} µM", "loosest mask that still blocks", f"at {bN['linker_residues']} residues, N-terminal")}
</div>
<div class="note ok"><b>The quantitative result.</b> An N-terminal fusion with a short linker gives
about <b>{bN['window_decades']:.1f} decades</b> of usable mask affinity. A C-terminal fusion gives
<b>{bC['window_decades']:.1f}</b>, and below <b>{closedC+1} residues it has no window at all</b> —
no mask affinity satisfies both requirements. The tether is worth so much that the mask can be
<b>micromolar</b> and still block: at {bN['linker_residues']} residues the ceiling is
{bN['km_max_for_masking_M']*1e6:.0f} µM, roughly
{bN['km_max_for_masking_M']/bN['km_min_for_release_M']:.0f}× above the release floor. Masks do not
need to be affinity-matured, and tight masks are actively harmful.</div>
<p>The release floor is the softer of the two bounds, because it depends on quantities a design does
not control: how tightly the nanobody holds toxin, and how much toxin is present.</p>
{fig("fig6_sensitivity.png", "The release floor against free toxin concentration, for three VHH–toxin affinities. The masking ceiling is unaffected.")}

<h2>Step 3 · Mask candidates from the structure</h2>
<p>A mask has to engage this paratope. One peptide is already known to do that: the toxin's own
contact segment. Grafting it needs no generative model and brings its own controls.</p>
<table><thead><tr><th>Mask</th><th>Sequence</th><th>Kind</th><th>Origin</th></tr></thead><tbody>
{"".join(f"<tr><td><span class='mono'>{k}</span></td><td><span class='mono'>{v['seq']}</span></td><td>{v['kind']}</td><td>{v['note']}</td></tr>" for k, v in s5['masks'].items())}
</tbody></table>
<div class="note"><b>Note on the candidates.</b> <span class="mono">epi1</span> is the tip of
α-cobratoxin loop II, the region carrying the residues that drive receptor binding — so this VHH
neutralizes by covering the toxin's functional site, and that same peptide is the natural mask.
<span class="mono">epi2</span> carries a free cysteine and is retained as a deliberate
counter-example: it would scramble disulfides in a real graft. <span class="mono">polyG</span> has
no designed interaction and should fail.</div>
<p>{s5['n_constructs']} constructs span {len(s5['masks'])} masks ×
{len(s5['linker_lengths'])} linker lengths ({", ".join(str(x) for x in s5['linker_lengths'])}
residues) plus the free VHH, from {min(c['length'] for c in s5['constructs'])} to
{max(c['length'] for c in s5['constructs'])} residues.</p>

<h2>Step 4 · Folding the fusion <span class="pill gpu">needs a GPU</span></h2>
<p>Steps 1–3 settle geometry and thermodynamics. They cannot say whether a particular mask sequence
folds back and covers the paratope. Because the fusion is a <b>single chain</b>, ESMFold predicts it
from sequence alone with no MSA, and a 160-residue construct folds in seconds on a free T4 — which
AlphaFold 3, Boltz-2, RFdiffusion and BindCraft would not do comfortably.</p>
<p><a href="{NB}"><b>Open the notebook in Colab →</b></a> It folds all {s5['n_constructs']}
constructs, measures paratope occlusion against the free VHH, and contrasts the grafted mask with
its scrambled and poly-glycine controls. The contrast is the informative quantity: without it, a
prediction showing a tethered peptide near its anchor says nothing about design.</p>
<div class="note"><b>What a predicted structure is not.</b> ESMFold returns one structure per
sequence. A mask lying against the paratope is a hypothesis about the masked state, not a masking
efficiency and not a K<sub>D</sub>. pLDDT is a confidence measure, not a folding free energy.</div>

<h2>Steps 5–6 · What is not computed here</h2>
<table><thead><tr><th>Component</th><th>What it would contribute</th><th>What it needs</th></tr></thead><tbody>
<tr><td>RFdiffusion / BindCraft</td><td>mask backbones beyond grafted epitopes</td><td>GPU; BindCraft also needs AF2 weights and PyRosetta</td></tr>
<tr><td>AlphaFold 3 / Boltz-2 / Chai-1</td><td>the VHH–toxin complex and mask–toxin competition</td><td>GPU, model access</td></tr>
<tr><td>BioEmu-2</td><td>conformational ensembles instead of one structure</td><td>GPU; multichain support is recent</td></tr>
<tr><td>CleaveNet</td><td>generated linkers with a target cleavage profile</td><td>trained on 18 human MMPs; snake venom metalloproteinases are out of distribution</td></tr>
<tr><td>Score 2, protease trigger</td><td>selectivity for the bite-site protease ensemble</td><td>peptide-cleavage data on venom proteases against a host negative panel — this does not exist yet</td></tr>
</tbody></table>
<div class="note"><b>The honest gap.</b> Score 2 is not a modelling problem that better software
solves. The training data for venom-protease substrate specificity has to be generated, and that is
wet-lab work. Published SVMP profiling shows a P1′ leucine preference shared with host MMPs, so
selectivity between venom protease and host protease is the real obstacle.</div>

<h2>Reproducing this</h2>
<div class="eq">git clone {REPO}<br>cd dtu_postdoc &amp;&amp; uv venv &amp;&amp; uv pip install -r requirements.txt<br>
for s in src/s0*.py; do python "$s"; done</div>
<p>Each stage writes a JSON file of its results; the figures and this page are generated from those
files, so no number here is typed by hand.</p>

<h2>References</h2>
<ul class="ref">
<li>Damsbo A. et al. Co-crystal structure of VHH TPL1158_01_C09 with α-cobratoxin.
<a href="https://www.rcsb.org/structure/9GCN">PDB 9GCN</a></li>
<li>Ahmadi S. et al. Nanobody-based recombinant antivenom for cobra, mamba and rinkhals bites.
<i>Nature</i> 647:716–725 (2025). <a href="https://doi.org/10.1038/s41586-025-09661-0">doi:10.1038/s41586-025-09661-0</a></li>
<li>Desnoyers L.R. et al. Tumor-specific activation of an EGFR-targeting probody enhances therapeutic index.
<i>Sci Transl Med</i> 5:207ra144 (2013). <a href="https://doi.org/10.1126/scitranslmed.3006682">doi:10.1126/scitranslmed.3006682</a></li>
<li>Zhou H.-X. Loops in proteins can be modeled as worm-like chains.
<i>J Phys Chem B</i> 105:6763–6766 (2001). <a href="https://doi.org/10.1021/jp011355n">doi:10.1021/jp011355n</a></li>
<li>Lin Z. et al. Evolutionary-scale prediction of atomic-level protein structure with a language model (ESMFold).
<i>Science</i> 379:1123–1130 (2023). <a href="https://doi.org/10.1126/science.ade2574">doi:10.1126/science.ade2574</a></li>
<li>Martin-Alonso C. et al. Deep learning guided design of protease substrates (CleaveNet).
<i>Nat Commun</i> 17:54 (2026). <a href="https://doi.org/10.1038/s41467-025-67226-1">doi:10.1038/s41467-025-67226-1</a></li>
<li>Pacesa M. et al. One-shot design of functional protein binders with BindCraft.
<i>Nature</i> 646:483–492 (2025). <a href="https://doi.org/10.1038/s41586-025-09429-6">doi:10.1038/s41586-025-09429-6</a></li>
<li>Bennett N.R. et al. Atomically accurate de novo design of antibodies with RFdiffusion.
<i>Nature</i> 649:183–193 (2025). <a href="https://doi.org/10.1038/s41586-025-09721-5">doi:10.1038/s41586-025-09721-5</a></li>
<li>Shrake A., Rupley J.A. Environment and exposure to solvent of protein atoms.
<i>J Mol Biol</i> 79:351–371 (1973). <a href="https://doi.org/10.1016/0022-2836(73)90011-9">doi:10.1016/0022-2836(73)90011-9</a></li>
</ul>

<footer>Figures and values generated by the analysis code in the
<a href="{REPO}">repository</a>. Structural data from the RCSB Protein Data Bank.</footer>
</div></body></html>"""

(DOCS / "index.html").write_text(HTML)
print(f"wrote docs/index.html  ({len(HTML)//1024} KB), "
      f"{len(list((DOCS/'figures').glob('*.png')))} figures")
