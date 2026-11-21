"""Rendering the results page.

Values are read from the stage outputs, so nothing on the page is typed by hand.
"""
from __future__ import annotations
import shutil
from pathlib import Path

from .config import (DOCS, FIGURES, STRUCTURES, DATA, PDB_ID, REPO_URL,
                     NOTEBOOK_URL, load_json)

CSS = """
:root{
  --ink:#12191F; --ink-2:#3A464F; --ink-3:#6E7880;
  --ground:#FFFFFF; --panel:#FFFFFF; --rule:#DDE4E7; --rule-2:#EBF0F2;
  --seq:#2D5BD1; --str:#C06014; --null:#6E7880; --chrome:#0E7C7B;
  --seq-wash:#EAEFFB; --str-wash:#FAEEE5;
  --sans:Inter,-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace;
}
*{box-sizing:border-box}
body{margin:0;background:var(--ground);color:var(--ink);font-family:var(--sans);
     font-size:16.5px;line-height:1.62;-webkit-font-smoothing:antialiased}
.page{max-width:1140px;margin:0 auto;padding:0 36px}
h1,h2,h3{text-wrap:balance;line-height:1.16;margin:0}
h1{font-size:clamp(2rem,5vw,3.1rem);font-weight:650;letter-spacing:-.025em}
h2{font-size:clamp(1.25rem,2.3vw,1.7rem);font-weight:650;letter-spacing:-.018em}
h3{font-size:1.06rem;font-weight:600}
p{margin:0;max-width:74ch}
a{color:var(--seq)}
.num{font-family:var(--mono);font-variant-numeric:tabular-nums}
.mono{font-family:var(--mono)}
.dim{color:var(--ink-3)}
.sec{padding:52px 0;border-top:1px solid var(--rule);display:flex;flex-direction:column;gap:24px}
.stack{display:flex;flex-direction:column;gap:12px}
.eyebrow{font-family:var(--mono);font-size:.67rem;letter-spacing:.16em;text-transform:uppercase;
         color:var(--chrome);font-weight:600}
header{padding:76px 0 50px;display:flex;flex-direction:column;gap:22px}
.kicker{font-family:var(--mono);font-size:.7rem;letter-spacing:.16em;text-transform:uppercase;
        color:var(--ink-3);display:flex;gap:12px;flex-wrap:wrap}
.kicker b{color:var(--seq);font-weight:600}
.lede{font-size:1.2rem;color:var(--ink-2);line-height:1.52;font-weight:400}
.banner{background:var(--panel);border:1px solid var(--rule);border-left:3px solid var(--chrome);
        border-radius:3px;padding:16px 20px;font-size:.96rem;color:var(--ink-2);max-width:none}
.banner b{color:var(--ink)}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:14px}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:18px 20px;
      display:flex;flex-direction:column;gap:5px}
.card .big{font-family:var(--mono);font-variant-numeric:tabular-nums;font-size:1.72rem;
           font-weight:500;letter-spacing:-.02em;line-height:1.05}
.card .cl{font-family:var(--mono);font-size:.63rem;letter-spacing:.1em;text-transform:uppercase;
          color:var(--ink-3)}
.card .sub{font-size:.85rem;color:var(--ink-2);line-height:1.45}
.card.hl{border-top:3px solid var(--seq)}
.chart{background:var(--panel);border:1px solid var(--rule);border-radius:3px;
       padding:18px 18px 12px;display:flex;flex-direction:column;gap:9px;overflow-x:auto}
.chart img,.chart svg{width:100%;height:auto;display:block;max-width:760px;margin:0 auto}
.chart figcaption{font-size:.85rem;color:var(--ink-3);line-height:1.5}
.charts2{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:16px}
.charts2 .chart img{max-width:none}
.chart.wide img{max-width:none}
.scroll{overflow-x:auto;background:var(--panel);border:1px solid var(--rule);border-radius:3px}
table{border-collapse:collapse;width:100%;font-family:var(--mono);
      font-variant-numeric:tabular-nums;font-size:.76rem}
th,td{padding:9px 13px;text-align:right;white-space:nowrap;border-bottom:1px solid var(--rule-2)}
th:first-child,td:first-child{text-align:left;font-family:var(--sans);font-size:.88rem;
      white-space:normal}
thead th{font-size:.62rem;letter-spacing:.08em;text-transform:uppercase;color:var(--ink-3);
         font-weight:600;border-bottom:1px solid var(--rule)}
tbody tr:last-child td{border-bottom:0}
tr.hl td{background:var(--seq-wash)}
tr.ctrl td{color:var(--ink-3)}
.tag{font-family:var(--mono);font-size:.59rem;letter-spacing:.07em;text-transform:uppercase;
     color:var(--ink-3);border:1px solid var(--rule);border-radius:2px;padding:1px 7px}
.tag.run{color:var(--chrome);border-color:#A9CFCE}
.tag.gap{color:var(--str);border-color:#E3C3A8}
.steps{display:grid;gap:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px}
.step{display:grid;grid-template-columns:2.2rem 1fr auto;gap:16px;align-items:start;
      padding:15px 20px;border-bottom:1px solid var(--rule-2)}
.step:last-child{border-bottom:0}
.step .n{font-family:var(--mono);font-size:.78rem;color:var(--ink-3);padding-top:2px}
.step b{display:block;font-size:.98rem}
.step span{font-size:.87rem;color:var(--ink-2);line-height:1.45}
.eq{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:15px 20px;
    font-family:var(--mono);font-size:.82rem;color:var(--ink);overflow-x:auto;line-height:1.9}
.viewer{background:var(--panel);border:1px solid var(--rule);border-radius:3px;padding:18px;
        display:flex;flex-direction:column;gap:12px}
.vstage{position:relative;width:100%;height:420px;border:1px solid var(--rule-2);
        border-radius:2px;background:#FFFFFF;overflow:hidden}
.vbar{display:flex;gap:7px;flex-wrap:wrap;align-items:center}
.vbar select,.vbar button{font-family:var(--mono);font-size:.72rem;padding:5px 10px;
    border:1px solid var(--rule);border-radius:2px;background:var(--panel);color:var(--ink-2);
    cursor:pointer}
.vbar button:hover,.vbar select:hover{border-color:var(--ink-3)}
.vkey{display:flex;gap:16px;flex-wrap:wrap;font-family:var(--mono);font-size:.66rem;
      color:var(--ink-3)}
.vkey i{display:inline-block;width:12px;height:3px;border-radius:2px;margin-right:6px;
        vertical-align:middle}
ul.ref{list-style:none;padding:0;margin:0;font-size:.88rem}
ul.ref li{padding:9px 0;border-bottom:1px solid var(--rule-2);color:var(--ink-2)}
ul.ref li:last-child{border-bottom:0}
footer{padding:36px 0 72px;color:var(--ink-3);font-size:.85rem}
@media(max-width:720px){.page{padding:0 20px}.vstage{height:320px}}
"""

VIEWER_JS = """
const MODELS = __MODELS__;
const PARATOPE = __PARATOPE__;
let v1 = null, v2 = null;

// Sticks are drawn with element colouring on a region-tinted carbon, so
// nitrogen, oxygen and sulfur read as themselves while the residue still
// says which part of the construct it belongs to.
const STICK_PARATOPE = {stick: {colorscheme: 'orangeCarbon', radius: 0.19}};
const STICK_MASK = {stick: {colorscheme: 'greenCarbon', radius: 0.17}};

function termini(view, model, sel, labelSuffix) {
  const cas = model.selectedAtoms(Object.assign({atom: 'CA'}, sel));
  if (!cas.length) return;
  const first = cas[0], last = cas[cas.length - 1];
  const mk = (atom, text, colour) => {
    view.addSphere({center: {x: atom.x, y: atom.y, z: atom.z}, radius: 1.1,
                    color: colour});
    view.addLabel(text, {position: {x: atom.x, y: atom.y, z: atom.z},
      backgroundColor: colour, backgroundOpacity: 0.92, fontColor: 'white',
      fontSize: 11, borderThickness: 0, inFront: true, alignment: 'center'});
  };
  mk(first, 'N' + labelSuffix, '#12191F');
  mk(last, 'C' + labelSuffix, '#6E7880');
}

function boot() {
  v1 = $3Dmol.createViewer('stage-solved', {backgroundColor: '#FFFFFF'});
  jQuery.get('structures/9GCN_AX.pdb', function (data) {
    const m = v1.addModel(data, 'pdb');
    v1.setStyle({chain: 'A'}, {cartoon: {color: '#2D5BD1'}});
    v1.setStyle({chain: 'X'}, {cartoon: {color: '#C06014'}});
    v1.addStyle({chain: 'A', resi: PARATOPE}, STICK_PARATOPE);
    termini(v1, m, {chain: 'A'}, '');
    termini(v1, m, {chain: 'X'}, '\u2032');
    v1.zoomTo(); v1.zoom(1.75); v1.render();
  });

  v2 = $3Dmol.createViewer('stage-designed', {backgroundColor: '#FFFFFF'});
  loadDesigned(document.getElementById('pick').value);
  document.getElementById('pick').addEventListener('change',
    e => loadDesigned(e.target.value));
  document.querySelectorAll('[data-spin]').forEach(b =>
    b.addEventListener('click', () => {
      const v = b.dataset.spin === 'solved' ? v1 : v2;
      v.spin(v.isSpinning && v.isSpinning() ? false : 'y');
    }));
}

function loadDesigned(id) {
  const spec = MODELS[id];
  v2.clear();
  jQuery.get('structures/' + id + '.pdb', function (data) {
    const m = v2.addModel(data, 'pdb');
    v2.setStyle({}, {cartoon: {color: '#C9D1D6'}});
    if (spec.vhh) v2.setStyle({resi: spec.vhh}, {cartoon: {color: '#2D5BD1'}});
    if (spec.linker) v2.setStyle({resi: spec.linker}, {cartoon: {color: '#8A9BA8'}});
    if (spec.mask) {
      v2.setStyle({resi: spec.mask}, {cartoon: {color: '#0E7C7B'}});
      v2.addStyle({resi: spec.mask}, STICK_MASK);
    }
    if (spec.paratope) v2.addStyle({resi: spec.paratope}, STICK_PARATOPE);
    termini(v2, m, {}, '');
    v2.zoomTo(); v2.zoom(1.75); v2.render();
    document.getElementById('dcap').textContent = spec.caption;
  });
}
window.addEventListener('load', boot);
"""


def _card(value, label, sub="", highlight=False):
    return (f'<div class="card{" hl" if highlight else ""}">'
            f'<div class="cl">{label}</div><div class="big">{value}</div>'
            + (f'<div class="sub">{sub}</div>' if sub else "") + "</div>")


def _chart(name, caption):
    return (f'<figure class="chart"><img src="figures/{name}.svg" alt="{caption}">'
            f'<figcaption>{caption}</figcaption></figure>')


def copy_assets(library):
    """Copy figures and the structures the viewers need into docs/."""
    (DOCS / "figures").mkdir(exist_ok=True)
    (DOCS / "structures").mkdir(exist_ok=True)
    for svg in FIGURES.glob("*.svg"):
        shutil.copy(svg, DOCS / "figures" / svg.name)
    for item in library:
        src = STRUCTURES / f"{item['id']}.pdb"
        if src.exists():
            shutil.copy(src, DOCS / "structures" / src.name)
    return DOCS


def write_reference_pdb(chains=("A", "X")):
    """One VHH-antigen pair from the deposited file, for the structure viewer."""
    import biotite.structure.io.pdb as pdb
    from .structure import read_atoms, fetch_structure
    atoms = read_atoms(fetch_structure(PDB_ID))
    subset = atoms[__import__("numpy").isin(atoms.chain_id, list(chains))]
    out = DOCS / "structures" / f"{PDB_ID}_{''.join(chains)}.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    f = pdb.PDBFile()
    f.set_structure(subset)
    f.write(str(out))
    return out


# ---------------------------------------------------------------- schematic
SEQ, STR, CHROME, NULL = "#2D5BD1", "#C06014", "#0E7C7B", "#6E7880"
INK, INK2, INK3, RULE = "#12191F", "#3A464F", "#6E7880", "#DDE4E7"


def _box(x, y, w, h, accent, number, title, lines, status, status_colour):
    """One module of the architecture diagram."""
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="#FFFFFF" '
           f'stroke="{RULE}" stroke-width="1"/>',
           f'<rect x="{x}" y="{y}" width="{w}" height="3" fill="{accent}"/>',
           f'<circle cx="{x+20}" cy="{y+26}" r="9.5" fill="none" stroke="{accent}" '
           f'stroke-width="1.2"/>',
           f'<text x="{x+20}" y="{y+29.5}" class="d-num" fill="{accent}" '
           f'text-anchor="middle">{number}</text>',
           f'<text x="{x+38}" y="{y+30}" class="d-title" fill="{INK}">{title}</text>']
    ty = y + 52
    for line in lines:
        bullet = line.startswith("• ")
        text = line[2:] if bullet else line
        cls = "d-item" if bullet else "d-sub"
        fill = INK2 if bullet else INK3
        if bullet:
            out.append(f'<circle cx="{x+17}" cy="{ty-3.5}" r="1.8" fill="{accent}"/>')
        out.append(f'<text x="{x + (25 if bullet else 16)}" y="{ty}" class="{cls}" '
                   f'fill="{fill}">{text}</text>')
        ty += 15 if bullet else 13
    out.append(f'<text x="{x+w-14}" y="{y+h-13}" class="d-tag" fill="{status_colour}" '
               f'text-anchor="end">{status}</text>')
    return "".join(out)


def _arrow(x1, y1, x2, y2, colour=NULL, dashed=False, label=None, lx=None, ly=None):
    dash = ' stroke-dasharray="5 4"' if dashed else ""
    out = [f'<path d="M {x1} {y1} L {x2} {y2}" stroke="{colour}" stroke-width="1.4" '
           f'fill="none" marker-end="url(#ah)"{dash}/>']
    if label:
        out.append(f'<text x="{lx}" y="{ly}" class="d-edge" fill="{INK3}" '
                   f'text-anchor="middle">{label}</text>')
    return "".join(out)


def pipeline_schematic(facts: dict) -> str:
    """The pipeline as built, in the order it runs. ``facts`` supplies the
    numbers so the diagram cannot drift from the stage outputs."""
    W, H = 1120, 650
    f = facts
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Pipeline as built">',
         '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
         f'markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.d-title{font:600 13.5px Inter,system-ui,sans-serif}'
         '.d-num{font:600 9.5px ui-monospace,monospace}'
         '.d-item{font:12px Inter,system-ui,sans-serif}'
         '.d-sub{font:10.5px ui-monospace,monospace;letter-spacing:.03em}'
         '.d-tool{font:600 10px ui-monospace,monospace;letter-spacing:.05em}'
         '.d-tag{font:600 9px ui-monospace,monospace;letter-spacing:.09em}'
         '.d-edge{font:10px ui-monospace,monospace;letter-spacing:.04em}'
         '.d-head{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '.d-band{font:600 11px Inter,system-ui,sans-serif}'
         '</style>']

    def stage(x, y, w, h, accent, num, title, tool, out_lines, tag, tag_col):
        o = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>',
             f'<rect x="{x}" y="{y}" width="{w}" height="3" fill="{accent}"/>',
             f'<text x="{x+14}" y="{y+24}" class="d-num" fill="{accent}">{num}</text>',
             f'<text x="{x+36}" y="{y+25}" class="d-title" fill="{INK}">{title}</text>',
             f'<text x="{x+14}" y="{y+43}" class="d-tool" fill="{accent}">{tool}</text>']
        ty = y + 62
        for line in out_lines:
            o.append(f'<text x="{x+14}" y="{ty}" class="d-item" fill="{INK2}">{line}</text>')
            ty += 16
        o.append(f'<text x="{x+w-12}" y="{y+h-11}" class="d-tag" fill="{tag_col}" '
                 f'text-anchor="end">{tag}</text>')
        return "".join(o)

    # band 1 - the structure constrains the design space before anything is generated
    s.append(f'<text x="0" y="13" class="d-head" fill="{CHROME}">'
             'STEP 1 — LET THE STRUCTURE SET THE LIMITS</text>')
    s.append(stage(0, 26, 345, 136, SEQ, "01", "Measure the complex",
                   "biotite · Shrake-Rupley", [
                       f"paratope: {f['n_paratope']} residues, {f['cdr3_frac']} in CDR3",
                       f"N-terminus {f['span_n']} Å from it, C-term {f['span_c']} Å",
                       "4 copies in the crystal agree",
                   ], "RUN", CHROME))
    s.append(stage(388, 26, 345, 136, SEQ, "02", "Price the tether",
                   "worm-like chain", [
                       "C_eff, masking ceiling, release floor",
                       f"N-terminal window {f['win_n']} decades wide",
                       f"C-terminal closed below {f['closed_c']+1} residues",
                   ], "RUN", CHROME))
    s.append(f'<rect x="776" y="26" width="344" height="136" rx="3" fill="#FFFFFF" '
             f'stroke="{CHROME}" stroke-width="1.4"/>')
    s.append(f'<rect x="776" y="26" width="344" height="3" fill="{CHROME}"/>')
    s.append(f'<text x="790" y="47" class="d-head" fill="{CHROME}">WHAT THIS DECIDES</text>')
    for n, line in enumerate([
            f"attach the mask at the N-terminus",
            f"use a {f['best_linker']}-residue linker",
            f"aim for a mask K_D near {f['km_ceiling_uM']:.0f} µM — not tighter"]):
        s.append(f'<text x="790" y="{70+n*19}" class="d-item" fill="{INK}">• {line}</text>')
    s.append(f'<text x="790" y="143" class="d-sub" fill="{INK3}">'
             'decided before any candidate exists</text>')
    s.append(_arrow(345, 94, 382, 94))
    s.append(_arrow(733, 94, 770, 94))

    # band 2 - generate and score
    s.append(f'<text x="0" y="203" class="d-head" fill="{CHROME}">'
             'STEP 2 — BUILD CANDIDATES AND SCORE THEM</text>')
    s.append(stage(0, 216, 262, 150, SEQ, "03", "Build the library",
                   "antigen epitope graft", [
                       f"{f['n_masks']} masks x {f['n_linkers']} linker lengths",
                       "scrambled, truncated and",
                       "inert masks as controls",
                   ], "RUN", CHROME))
    s.append(stage(286, 216, 262, 150, SEQ, "04", "Fold and measure",
                   "ESMFold · CPU", [
                       f"{f['n_folded']} constructs folded",
                       "paratope occlusion against",
                       "the folded free nanobody",
                   ], "RUN", CHROME))
    s.append(stage(572, 216, 262, 150, SEQ, "05", "Redesign the mask",
                   "ProteinMPNN · CPU", [
                       "nanobody and linker held fixed,",
                       "mask positions rewritten",
                       f"{f['n_redesigns']} distinct masks, then refolded",
                   ], "RUN", CHROME))
    s.append(stage(858, 216, 262, 150, SEQ, "06", "Test the competition",
                   "Boltz-2 · CPU", [
                       "two-chain prediction: with the",
                       "mask in place, can the toxin",
                       "still reach the paratope?",
                   ], f['boltz_tag'], CHROME if f['boltz_tag'] == "RUN" else STR))
    for x in (262, 548, 834):
        s.append(_arrow(x, 291, x + 20, 291))

    # the missing module
    s.append(f'<rect x="180" y="404" width="470" height="118" rx="3" fill="#FFFFFF" '
             f'stroke="{STR}" stroke-dasharray="5 4"/>')
    s.append(f'<rect x="180" y="404" width="470" height="3" fill="{STR}"/>')
    s.append(f'<text x="194" y="426" class="d-num" fill="{STR}">2a</text>')
    s.append(f'<text x="216" y="427" class="d-title" fill="{INK}">Score the trigger</text>')
    s.append(f'<text x="194" y="445" class="d-tool" fill="{STR}">NOT BUILT</text>')
    for n, line in enumerate([
            "Is the linker cut by the protease at the target site, and left alone elsewhere?",
            "The module is swappable — pH, metal ion or protease would each slot in here.",
            "It is missing because the cleavage data for these proteases does not exist yet."]):
        s.append(f'<text x="194" y="{465+n*17}" class="d-item" fill="{INK2}">{line}</text>')

    # decide
    s.append(stage(690, 404, 430, 118, SEQ, "07", "Rank, then hand over",
                   "Pareto over the scores", [
                       "keep the trade-offs separate, do not average them",
                       "keep structurally different designs",
                       "top candidates go to the lab",
                   ], "PARTIAL", STR))
    s.append(_arrow(989, 366, 989, 400))
    s.append(f'<path d="M 650 463 L 684 463" stroke="{STR}" stroke-width="1.4" '
             f'stroke-dasharray="5 4" fill="none" marker-end="url(#ah)"/>')

    # feedback
    s.append(f'<path d="M 905 522 L 905 548 L 70 548 L 70 372" stroke="{CHROME}" '
             f'stroke-width="1.4" fill="none" stroke-dasharray="5 4" '
             f'marker-end="url(#ah)"/>')
    s.append(f'<text x="500" y="542" class="d-edge" fill="{CHROME}" text-anchor="middle">'
             'measurements come back and retrain the scorers</text>')

    # cost band
    s.append(f'<rect x="0" y="566" width="{W}" height="76" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="18" y="588" class="d-head" fill="{CHROME}">'
             'WHAT THE FEEDBACK LOOP COSTS AT 3 MEASUREMENTS A MONTH</text>')
    for i2, (a, b) in enumerate([("20–30", "a validation set — 7 to 10 months"),
                                 ("~100", "optimising against it — nearly 3 years"),
                                 ("100–1000", "training on it directly — longer still")]):
        x = 18 + i2 * 366
        s.append(f'<text x="{x}" y="{614}" class="d-band" fill="{INK}">{a} measurements</text>')
        s.append(f'<text x="{x}" y="{632}" class="d-sub" fill="{STR}">{b}</text>')
    s.append("</svg>")
    return "".join(s)


def _unused_pipeline_schematic() -> str:
    """The proposed architecture, with what is built marked as built."""
    W, H = 1120, 600
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Pipeline architecture">',
         '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
         f'markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.d-title{font:600 14px Inter,system-ui,sans-serif}'
         '.d-num{font:600 10px ui-monospace,monospace}'
         '.d-item{font:12.5px Inter,system-ui,sans-serif}'
         '.d-sub{font:11px ui-monospace,monospace;letter-spacing:.04em}'
         '.d-tag{font:600 9.5px ui-monospace,monospace;letter-spacing:.09em}'
         '.d-edge{font:10.5px ui-monospace,monospace;letter-spacing:.05em}'
         '.d-head{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '.d-note{font:12px Inter,system-ui,sans-serif}'
         '</style>']

    s.append(f'<text x="0" y="14" class="d-head" fill="{CHROME}">'
             'WHAT THE PIPELINE DOES</text>')

    # 1 generator
    s.append(_box(0, 34, 300, 182, SEQ, "1", "Generate candidates", [
        "• unconditional — plausible nanobodies",
        "• conditional — biased toward function",
        "• evolutionary — vary a starting point",
        "",
        "HERE: epitope graft + ProteinMPNN",
        "redesign on the folded backbone",
    ], "BUILT", CHROME))

    # 2a trigger
    s.append(_box(410, 34, 300, 182, STR, "2a", "Score the trigger", [
        "• does the local protease cut the linker",
        "• and leave it alone everywhere else",
        "",
        "swappable: pH, metal ion, protease",
        "NEEDS: cleavage data for the",
        "relevant proteases",
    ], "NOT BUILT", STR))

    # 2b binding
    s.append(_box(410, 244, 300, 182, SEQ, "2b", "Score the binding", [
        "• masked state — is the paratope covered",
        "• released state — can the target bind",
        "",
        "HERE: tether thermodynamics,",
        "ESMFold occlusion, Boltz-2 complex",
    ], "BUILT", CHROME))

    # 3 ranking
    s.append(_box(820, 139, 300, 182, SEQ, "3", "Rank and choose", [
        "• combine scores without averaging them",
        "• keep structurally diverse designs",
        "",
        "HERE: scores computed per design;",
        "Pareto selection not yet a stage",
    ], "PARTIAL", STR))

    # edges
    s.append(_arrow(300, 110, 404, 110, label="candidates", lx=352, ly=100))
    s.append(_arrow(300, 140, 404, 320, label=None))
    s.append(_arrow(710, 110, 814, 190, label="trigger score", lx=762, ly=128))
    s.append(_arrow(710, 320, 814, 260, label="binding score", lx=762, ly=310))

    # lab + feedback
    s.append(f'<rect x="820" y="366" width="300" height="56" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}" stroke-dasharray="5 4"/>')
    s.append(f'<text x="836" y="390" class="d-title" fill="{INK2}">Top candidates to the lab</text>')
    s.append(f'<text x="836" y="408" class="d-sub" fill="{INK3}">express, measure, feed back</text>')
    s.append(_arrow(970, 321, 970, 362))
    s.append(f'<path d="M 970 422 L 970 450 L 150 450 L 150 222" stroke="{CHROME}" '
             f'stroke-width="1.4" fill="none" stroke-dasharray="5 4" '
             f'marker-end="url(#ah)"/>')
    s.append(f'<text x="560" y="443" class="d-edge" fill="{CHROME}" text-anchor="middle">'
             'measured results retrain the generator and the scorers</text>')

    # the data-rate band
    s.append(f'<rect x="0" y="482" width="{W}" height="104" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="18" y="506" class="d-head" fill="{CHROME}">'
             'WHAT THE FEEDBACK LOOP COSTS</text>')
    marks = [(0, "20–30 measurements", "a validation set, treated with care"),
             (1, "~100 measurements", "proper optimisation against that set"),
             (2, "100–1000 measurements", "enough to train on directly")]
    for i, (n, a, b) in enumerate(marks):
        x = 18 + i * 366
        s.append(f'<text x="{x}" y="{534}" class="d-title" fill="{INK}">{a}</text>')
        s.append(f'<text x="{x}" y="{552}" class="d-note" fill="{INK2}">{b}</text>')
        months_lo = {0: 7, 1: 33, 2: 33}[n]
        months_hi = {0: 10, 1: 33, 2: 333}[n]
        span = (f"{months_lo} months" if months_lo == months_hi
                else f"{months_lo}–{months_hi} months")
        s.append(f'<text x="{x}" y="{572}" class="d-sub" fill="{STR}">'
                 f'at 3 per month: {span}</text>')
    s.append("</svg>")
    return "".join(s)


def states_schematic() -> str:
    """What the molecule is supposed to do, and why that is hard."""
    W, H = 1120, 300
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="The three states">',
         '<defs><marker id="ah2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
         f'markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.s-t{font:600 13px Inter,system-ui,sans-serif}'
         '.s-i{font:12px Inter,system-ui,sans-serif}'
         '.s-m{font:10.5px ui-monospace,monospace;letter-spacing:.04em}'
         '.s-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '.s-lbl{font:600 11px ui-monospace,monospace}'
         '</style>',
         f'<text x="0" y="13" class="s-h" fill="{CHROME}">THE THREE STATES</text>']

    def nanobody(cx, cy, paratope_colour):
        return (f'<path d="M {cx-46} {cy-30} q -10 30 0 60 q 46 14 92 0 q 10 -30 0 -60 '
                f'q -46 -14 -92 0 z" fill="#EEF2F7" stroke="{INK3}" stroke-width="1"/>'
                f'<text x="{cx} " y="{cy+5}" class="s-lbl" fill="{INK2}" '
                f'text-anchor="middle">VHH</text>'
                f'<path d="M {cx-26} {cy-33} q 26 -9 52 0" stroke="{paratope_colour}" '
                f'stroke-width="4" fill="none" stroke-linecap="round"/>')

    def panel(x, title, note):
        return (f'<rect x="{x}" y="34" width="340" height="212" rx="3" fill="#FFFFFF" '
                f'stroke="{RULE}"/>'
                f'<text x="{x+16}" y="{58}" class="s-t" fill="{INK}">{title}</text>'
                f'<text x="{x+16}" y="{230}" class="s-i" fill="{INK2}">{note}</text>')

    # A - masked
    s.append(panel(0, "A · masked", "tethered mask sits over the paratope"))
    s.append(nanobody(170, 150, SEQ))
    s.append(f'<circle cx="170" cy="104" r="19" fill="{CHROME}" stroke="#FFFFFF" '
             f'stroke-width="2"/>')
    s.append(f'<text x="170" y="109" class="s-lbl" fill="#FFFFFF" text-anchor="middle">M</text>')
    s.append(f'<path d="M 140 128 q -12 -18 12 -26" stroke="{INK3}" stroke-width="1.6" '
             f'fill="none"/>')
    s.append(f'<text x="112" y="116" class="s-m" fill="{INK3}">L</text>')

    # B - cleaved
    s.append(panel(390, "B · cleaved", "protease cuts; the mask drifts off"))
    s.append(nanobody(560, 150, SEQ))
    s.append(f'<circle cx="648" cy="92" r="19" fill="{CHROME}" stroke="#FFFFFF" '
             f'stroke-width="2"/>')
    s.append(f'<text x="648" y="97" class="s-lbl" fill="#FFFFFF" text-anchor="middle">M</text>')
    s.append(f'<path d="M 576 118 L 622 96" stroke="{NULL}" stroke-width="1.4" '
             f'stroke-dasharray="4 3" marker-end="url(#ah2)" fill="none"/>')
    s.append(f'<text x="531" y="100" class="s-m" fill="{STR}">protease</text>')
    s.append(f'<path d="M 524 112 l 14 10 M 538 112 l -14 10" stroke="{STR}" '
             f'stroke-width="1.8"/>')

    # C - bound
    s.append(panel(780, "C · active", "the original paratope takes the toxin"))
    s.append(nanobody(950, 160, SEQ))
    s.append(f'<rect x="926" y="86" width="48" height="34" rx="6" fill="{STR}" '
             f'stroke="#FFFFFF" stroke-width="2"/>')
    s.append(f'<text x="950" y="108" class="s-lbl" fill="#FFFFFF" text-anchor="middle">T</text>')

    s.append(f'<rect x="0" y="258" width="{W}" height="36" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="16" y="281" class="s-i" fill="{INK}">'
             'A mask tight enough to hold state A may refuse to leave in state B. '
             'That single sentence is the whole design problem, and stage 02 turns it '
             'into a range of allowed mask affinities.</text>')
    s.append("</svg>")
    return "".join(s)


def occlusion_schematic(free_sasa: float | None = None) -> str:
    """How the masked state is actually measured."""
    W, H = 1120, 238
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="How occlusion is measured">',
         '<defs><marker id="ah3" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
         f'markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.o-t{font:600 12.5px Inter,system-ui,sans-serif}'
         '.o-i{font:11.5px Inter,system-ui,sans-serif}'
         '.o-m{font:10.5px ui-monospace,monospace;letter-spacing:.04em}'
         '.o-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>',
         f'<text x="0" y="13" class="o-h" fill="{CHROME}">HOW THE MASKED STATE IS MEASURED</text>']
    steps = [("fold the free nanobody", "ESMFold, sequence only",
              "measure the paratope's exposed surface"),
             ("fold the whole fusion", "mask + linker + nanobody",
              "measure the same residues again"),
             ("take the difference", "occlusion = 1 − fused / free",
              "0% means the mask changed nothing")]
    x = 0
    for i, (title, tool, note) in enumerate(steps):
        s.append(f'<rect x="{x}" y="30" width="330" height="118" rx="3" fill="#FFFFFF" '
                 f'stroke="{RULE}"/>')
        s.append(f'<rect x="{x}" y="30" width="330" height="3" fill="{SEQ}"/>')
        s.append(f'<text x="{x+16}" y="{56}" class="o-t" fill="{INK}">{title}</text>')
        s.append(f'<text x="{x+16}" y="{76}" class="o-m" fill="{SEQ}">{tool}</text>')
        s.append(f'<text x="{x+16}" y="{100}" class="o-i" fill="{INK2}">{note}</text>')
        if i < 2:
            s.append(f'<path d="M {x+330} 89 L {x+390} 89" stroke="{NULL}" '
                     f'stroke-width="1.4" marker-end="url(#ah3)" fill="none"/>')
        x += 395
    s.append(f'<rect x="0" y="162" width="{W}" height="64" rx="3" fill="#FFFFFF" '
             f'stroke="{STR}" stroke-dasharray="5 4"/>')
    s.append(f'<text x="16" y="184" class="o-t" fill="{INK}">'
             'Why the controls carry the argument</text>')
    s.append(f'<text x="16" y="203" class="o-i" fill="{INK2}">'
             'A tethered peptide sits near its own anchor whatever its sequence, so a '
             'single occlusion number means very little on its own.</text>')
    s.append(f'<text x="16" y="219" class="o-i" fill="{INK2}">'
             'Only the gap between the designed mask and its scrambled, truncated and '
             'poly-glycine twins separates a design from a geometric accident.</text>')
    s.append("</svg>")
    return "".join(s)
