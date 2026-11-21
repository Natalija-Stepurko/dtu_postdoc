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
.sec{padding:64px 0 72px;border-top:1px solid var(--rule);display:flex;
  flex-direction:column;gap:24px;scroll-margin-top:52px;position:relative}
.sec::after{content:attr(data-slide);position:absolute;top:22px;right:0;
  font-family:var(--mono);font-size:.64rem;letter-spacing:.14em;color:var(--ink-3)}
.navbar{position:fixed;left:0;right:0;bottom:0;z-index:60;background:rgba(255,255,255,.95);
  backdrop-filter:blur(8px);border-top:1px solid var(--rule)}
.navbar .inner{max-width:1140px;margin:0 auto;padding:7px 36px;display:flex;gap:6px;
  align-items:center}
.navbar .dots{display:flex;gap:5px;flex:1;flex-wrap:wrap}
.navbar .dots a{width:22px;height:22px;display:grid;place-items:center;border-radius:3px;
  font-family:var(--mono);font-size:.6rem;color:var(--ink-3);text-decoration:none;
  border:1px solid transparent}
.navbar .dots a:hover{border-color:var(--rule);color:var(--ink)}
.navbar .dots a[aria-current=true]{background:var(--ink);color:#fff}
.navbar .hint{font-family:var(--mono);font-size:.6rem;letter-spacing:.1em;color:var(--ink-3)}
body{padding-bottom:44px}
@media print{
  .navbar{display:none}
  body{padding-bottom:0;font-size:15px}
  .sec{padding:26px 0 30px;break-inside:auto}
  .sec::after{top:8px}
  figure.chart,.card,.viewer,.scroll,.kpis,.steps,details.fold,details.tech{
    break-inside:avoid;page-break-inside:avoid}
  h2,h3{break-after:avoid;page-break-after:avoid}
  .tfig{break-inside:avoid}
  .vstage{height:360px}
  a{text-decoration:none;color:inherit}
}
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
table.wrap{table-layout:fixed;width:100%}
table.wrap th,table.wrap td{white-space:normal;text-align:left;vertical-align:top;
  overflow-wrap:anywhere;
  font-family:var(--sans);font-size:.92rem;line-height:1.45;padding:12px 14px}
table.wrap thead th{font-family:var(--mono);font-size:.6rem}
table.wrap td:first-child{font-weight:600;color:var(--ink)}
.scroll:has(table.wrap){overflow-x:visible}
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
details.fold{background:var(--panel);border:1px solid var(--rule);border-radius:3px;
    margin:4px 0}
details.fold>summary{cursor:pointer;list-style:none;padding:16px 20px;display:flex;
    gap:14px;align-items:baseline;flex-wrap:wrap}
details.fold>summary::-webkit-details-marker{display:none}
details.fold>summary::before{content:"+";font-family:var(--mono);font-size:1rem;
    color:var(--str);width:14px;flex:none;line-height:1}
details.fold[open]>summary::before{content:"−"}
details.fold>summary b{font-size:1.02rem;font-weight:600;color:var(--ink)}
details.fold>summary span{font-size:.9rem;color:var(--ink-2);flex:1 1 24ch}
details.fold>summary:hover b{color:var(--str)}
.foldbody{padding:0 20px 22px;display:flex;flex-direction:column;gap:18px;
    border-top:1px solid var(--rule-2);margin-top:-2px;padding-top:20px}
.tfig{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;
    padding:14px;display:flex;flex-direction:column;gap:9px}
.tfig img{width:100%;height:auto;max-height:260px;object-fit:contain;display:block}
.tfig.tall img{max-height:340px}
.tfig.full img{max-height:none;max-width:980px;margin:0 auto}
.tfig.portrait img{max-height:720px;max-width:620px;margin:0 auto}
.tfig.wide img{max-height:190px}
.tfig figcaption{font-size:.82rem;color:var(--ink-3);line-height:1.45}
.tgrid{display:grid;gap:14px;margin:4px 0}
.tgrid.c2{grid-template-columns:repeat(auto-fit,minmax(280px,1fr))}
.tgrid.c3{grid-template-columns:repeat(auto-fit,minmax(215px,1fr))}
.tstep{display:flex;gap:12px;align-items:baseline;margin:26px 0 2px}
.tstep b{font-family:var(--mono);font-size:.68rem;letter-spacing:.14em;color:var(--chrome)}
.tstep span{font-size:1rem;color:var(--ink);font-weight:600}
.stat-tag{display:inline-block;font-family:var(--mono);font-size:.6rem;letter-spacing:.07em;
  text-transform:uppercase;padding:2px 8px;border-radius:999px;border:1px solid currentColor;
  white-space:nowrap}
h2 .stat-tag,h3 .stat-tag{vertical-align:middle;margin-left:12px;font-size:.55rem}
.modblock{margin:22px 0}
.modhead{display:flex;gap:16px;align-items:baseline;margin-bottom:10px}
.modnum{font-family:var(--mono);font-size:.72rem;color:var(--chrome);font-weight:650;
  border:1px solid var(--rule);border-radius:3px;padding:3px 9px}
.modhead b{font-size:1.08rem;display:block;color:var(--ink)}
.modhead span{font-size:.95rem;color:var(--ink-2);display:block;margin-top:2px}
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
SEQ_WASH, STR_WASH = "#EAEFFB", "#FAEEE5"


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


def architecture_schematic() -> str:
    """The computational pipeline. Three scores, then one stage that combines them."""
    W, H = 1120, 430
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Pipeline architecture">',
         '<defs><marker id="ah4" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
         f'markerHeight="6" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.a-t{font:600 13.5px Inter,system-ui,sans-serif}'
         '.a-n{font:600 9.5px ui-monospace,monospace}'
         '.a-i{font:12px Inter,system-ui,sans-serif}'
         '.a-s{font:10.5px ui-monospace,monospace;letter-spacing:.03em}'
         '.a-g{font:600 9px ui-monospace,monospace;letter-spacing:.09em}'
         '.a-e{font:10px ui-monospace,monospace;letter-spacing:.04em}'
         '.a-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>']

    def mod(x, y, w, h, accent, num, title, role, slots, tag, tagcol, dashed=False):
        d = ' stroke-dasharray="5 4"' if dashed else ""
        o = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="3" fill="#FFFFFF" '
             f'stroke="{accent if dashed else RULE}"{d}/>',
             f'<rect x="{x}" y="{y}" width="{w}" height="3" fill="{accent}"/>',
             f'<text x="{x+14}" y="{y+24}" class="a-n" fill="{accent}">{num}</text>',
             f'<text x="{x+38}" y="{y+25}" class="a-t" fill="{INK}">{title}</text>']
        ty = y + 46
        for line in role:
            o.append(f'<text x="{x+14}" y="{ty}" class="a-i" fill="{INK2}">{line}</text>')
            ty += 16
        ty += 4
        for line in slots:
            o.append(f'<text x="{x+14}" y="{ty}" class="a-s" fill="{accent}">{line}</text>')
            ty += 14
        o.append(f'<text x="{x+w-12}" y="{y+h-11}" class="a-g" fill="{tagcol}" '
                 f'text-anchor="end">{tag}</text>')
        return "".join(o)

    s.append(f'<text x="0" y="13" class="a-h" fill="{CHROME}">'
             'THE COMPUTATIONAL PIPELINE</text>')

    s.append(mod(0, 26, 538, 120, CHROME, "00", "Constrain", [
        "Measure the antigen complex, then work out what",
        "mask affinity and linker length could work at all.",
    ], ["input: an experimental or predicted complex"], "RUN", CHROME))

    s.append(mod(582, 26, 538, 120, SEQ, "01", "Generate", [
        "Propose mask sequences and linkers inside the",
        "range the previous stage allows.",
    ], ["graft \u00b7 ProteinMPNN \u00b7 RFdiffusion \u00b7 BindCraft"], "RUN", CHROME))

    # three scores, equal weight, side by side
    for n, (num, title, role, slots, tag, col) in enumerate([
        ("02", "Score: can it be made", [
            "Does it fold, does the nanobody survive the",
            "fusion, is the surface a liability?",
        ], ["ESMFold confidence \u00b7 core RMSD",
            "exposed hydrophobic \u00b7 free cysteine"], "RUN", CHROME),
        ("03", "Score: can it be cut", [
            "Is the linker reachable and extended, and does",
            "this protease recognise this sequence?",
        ], ["accessibility: run",
            "sequence recognition: needs data"], "HALF BUILT", STR),
        ("04", "Score: does it mask", [
            "Is the site blocked while tethered, and does the",
            "antigen return once the tether is cut?",
        ], ["occlusion: does not separate",
            "tethered complex: run"], "PART RUN", STR),
    ]):
        s.append(mod(n * 380, 176, 358, 134, SEQ, num, title, role, slots, tag, col))

    # Generate feeds all three scores through one bus
    s.append(f'<path d="M 851 146 L 851 158" stroke="{NULL}" stroke-width="1.4" fill="none"/>')
    s.append(f'<path d="M 179 158 L 939 158" stroke="{NULL}" stroke-width="1.4" '
             f'fill="none"/>')
    for n in range(3):
        s.append(f'<path d="M {n*380+179} 158 L {n*380+179} 170" stroke="{NULL}" '
                 f'stroke-width="1.4" fill="none" marker-end="url(#ah4)"/>')
    s.append(mod(291, 340, 538, 76, CHROME, "05", "Combine", [
        "Keep a design unless another beats it on every score. No weights.",
    ], [], "RUN", CHROME))
    for n in range(3):
        s.append(f'<path d="M {n*380+179} 310 L {n*380+179} 326 L 560 326 L 560 336" '
                 f'stroke="{NULL}" stroke-width="1.4" fill="none" '
                 f'marker-end="url(#ah4)"/>' if n != 1 else
                 f'<path d="M 560 310 L 560 336" stroke="{NULL}" stroke-width="1.4" '
                 f'fill="none" marker-end="url(#ah4)"/>')
    s.append(f'<path d="M 538 86 L 576 86" stroke="{NULL}" stroke-width="1.4" '
             f'fill="none" marker-end="url(#ah4)"/>')
    s.append("</svg>")
    return "".join(s)


def cautionary_section(ref: dict, tether: dict, fold: dict, complexes: dict | None,
                       label: str, pdb_id: str, chain: str) -> str:
    """A target the pipeline was run on and then set aside, folded away by default.

    Keeping it is the point: the same code on a different complex returns
    nothing usable, and the reason is geometry rather than code.
    """
    import numpy as np
    iface = ref["chain_pairs"][chain]
    geom = ref["terminal_geometry"]
    region = {k: v["n"] for k, v in ref["paratope_by_region"].items()}
    kinds: dict[str, list[float]] = {}
    for row in fold["results"]:
        kinds.setdefault(row["kind"], []).append(row["occlusion"] * 100)
    graft = float(np.mean(kinds.get("epitope graft", [0])))
    inert = float(np.mean(kinds.get("inert control", [0])))
    scram = float(np.mean(kinds.get("scrambled control", [0])))
    ctrl = None
    if complexes:
        ctrl = next((r for r in complexes["results"]
                     if r["case"] == "control_free_vhh" and r["status"] == "ok"), None)

    rows = "".join(
        f"<tr><td>{k}</td><td>{np.mean(v):+.1f}</td><td>{len(v)}</td></tr>"
        for k, v in sorted(kinds.items(), key=lambda kv: -float(np.mean(kv[1]))))
    ctrl_html = ""
    if ctrl:
        ctrl_html = f"""<p><b>The second control failed too.</b> Predicting the free nanobody with
its toxin should have reproduced a pose the crystal structure already shows. Only
<b>{ctrl['fraction_of_contacts_on_paratope']:.0%}</b> of the predicted contacts landed on the known
paratope, so the complex module could not adjudicate masking here either.</p>"""

    return f"""
<details class="fold">
<summary><b>Set aside: {pdb_id}, where the geometry defeats the measurement</b>
<span>The pipeline was run end to end on this complex first. Both controls failed, for a
reason worth keeping.</span></summary>
<div class="foldbody">
<p>{label}, {ref['metadata']['resolution_A']} Å. The paratope sits almost entirely in CDR3
({region['CDR3']} of {iface['n_paratope']} residues, with CDR2 contributing none) — and the
N-terminal residue <span class="mono">{geom['n_term_residue']}</span> is <b>itself part of the
paratope</b>, burying {iface['paratope_delta_sasa'][geom['n_term_residue']]:.0f} Å².</p>
<p>That makes masking trivially easy and measuring it impossible. Any peptide tethered at that
terminus starts on top of the binding site, whatever its sequence, so occlusion cannot tell a
designed mask from a random one.</p>
<div class="scroll"><table><thead><tr><th>mask kind</th><th>mean occlusion %</th><th>n</th>
</tr></thead><tbody>{rows}</tbody></table></div>
<p>The inert poly-glycine control occluded <b>{inert:+.1f}%</b> of the paratope against
<b>{graft:+.1f}%</b> for the designed grafts and {scram:+.1f}% for the scrambled ones. The control
scored <em>highest</em>. There is no separation to rank on.</p>
{ctrl_html}
<p><b>Why it is still here.</b> The constraint module had already flagged it: the N-terminus was
0.0 Å from the nearest paratope residue, against {geom['c_term_to_paratope_centroid_A']:.0f} Å for
the C-terminus. That number was visible before a single construct was folded, and it is the
criterion the replacement system was chosen on. A pipeline that measures the geometry first can
tell you that a target is unsuitable before you spend anything on it.</p>
</div>
</details>"""


def background_section(window_lo_nM: float, window_hi_nM: float) -> str:
    """Prior work, told as a sequence rather than a gallery.

    Follistatin's C-terminal tail is an endogenous tethered competitor for the
    same domain that binds heparan sulphate - the arrangement this pipeline
    designs on purpose - so the link is architectural rather than decorative.
    """
    def f(name, cap, cls=""):
        return (f'<figure class="tfig {cls}"><img src="thesis/{name}.png" alt="{cap}">'
                f'<figcaption>{cap}</figcaption></figure>')

    def step(n, title):
        return f'<div class="tstep"><b>{n}</b><span>{title}</span></div>'

    return f"""
<details class="fold">
<summary><b>Background: a tethered competitive regulator, measured</b>
<span>PhD work on follistatin, which solves this same two-state problem naturally.</span>
</summary>
<div class="foldbody">

<p>The construct on this page needs a tethered module that binds its own protein tightly enough to
block and weakly enough to let go. Follistatin already does this. Its isoforms differ only at the
C-terminus, and the longest, FST315, carries an acidic tail that folds back onto the FS1 domain —
the same surface that binds heparan sulphate. The tail is an <b>endogenous mask</b>.</p>

{step("01", "The molecule, and where the isoforms differ")}
{f("fst_domains", "Domain organisation, charge and conservation. FST288, FST303 and FST315 differ "
   "only in how much C-terminal tail they carry; the tail is strongly acidic (red, right) against "
   "a basic heparan-sulphate binding site (HSB) on FS1.", "full")}
{f("fst_topology", "Topology of the full-length molecule, from the N-terminus to the tail.", "wide")}

{step("02", "The tail binds back, and the strength can be measured")}
<div class="tgrid c2">
{f("fig4_1_3_a_isoforms", "Biolayer interferometry: tail peptides binding each isoform.")}
{f("fig4_1_8_d", "Fluorescence polarisation, tail against FS12 at three ionic strengths. "
   "K<sub>D</sub> = 145 nM at 50 mM NaCl, 8.2 µM at 150 mM, 34.6 µM at 300 mM — the interaction "
   "is largely electrostatic.")}
</div>

{step("03", "What the tail competes with")}
{f("fig4_1_4_e", "Isoform binding to heparan sulphate, the ligand the tail displaces from that "
   "same surface.", "wide")}

<div class="note good"><b>Those affinities land inside the window this pipeline computes.</b>
Masking and release together allow a mask dissociation constant between roughly
<b>{window_lo_nM:.0f} nM and {window_hi_nM:.0f} µM</b> for the nanobody here. The follistatin tail
works at 145 nM to 34.6 µM depending on ionic strength — inside that range at every salt
concentration measured. A natural tethered competitor sits where the model says a designed one
should. Different protein and different geometry, so this says the regime is physically sensible
rather than that the model is right.</div>

{step("04", "Two states, resolved experimentally rather than predicted")}
{f("fst_saxs_abc", "From crystal structure to solution ensemble. (a) FST288 in the "
   "activin A complex, fitted to the SAXS envelope. (b) A model of FST315 adding the tails. "
   "(c) The solution structure is not one conformation but two — 65% open and 35% closed, "
   "interconverting, with the tail engaged in the closed state.", "portrait")}
<div class="tgrid c2">
{f("fig4_1_14_f", "SEC-SAXS scattering with fits from each model. The two-state model gives the "
   "lowest chi-squared, 1.058, beating either single structure.", "tall")}
{f("fig4_1_23_a_isoforms", "Conformational space explored for each isoform, radius of gyration "
   "against maximum dimension. The longer the tail, the more closed states appear.", "tall")}
</div>

<div class="note"><b>Why this is the relevant experience.</b> The problem on this page is a
two-state one: hold a conformation while tethered, lose it when cut, at an affinity inside a narrow
window. That is the system characterised above — by SEC-SAXS, ensemble modelling and direct
affinity measurement — rather than predicted.</div>
</div>
</details>"""




def scale_section() -> str:
    """How large the campaign has to be, and how the lab makes it smaller.

    Counts are projections from published hit rates; runtimes are measured on
    this machine. Both are labelled as such on the page.
    """
    W, H = 1120, 352
    rows = [
        ("Candidate complexes", 6, 2, "geometry screen: is a terminus clear of the paratope?"),
        ("Mask backbones", 10000, 10000, "RFdiffusion, per surviving complex"),
        ("Sequences", 80000, 80000, "ProteinMPNN, 8 per backbone"),
        ("Fold and developability", 80000, 4000, "ESMFold; filters cut cost, not failure rate"),
        ("Two-state scoring", 4000, 400, "tethered complex: masked and released"),
        ("Pareto selection", 400, 30, "diverse, not just top-ranked"),
        ("Ordered", 30, 3, "first cycle"),
    ]
    import math
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Design funnel">',
         '<style>'
         '.n-t{font:600 12.5px Inter,system-ui,sans-serif}'
         '.n-s{font:11px Inter,system-ui,sans-serif}'
         '.n-v{font:600 12px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.n-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>',
         f'<text x="0" y="13" class="n-h" fill="{CHROME}">'
         'HOW BIG THE CAMPAIGN HAS TO BE, PER COMPLEX</text>']
    y, maxw = 32, 560
    for name, a, b, note in rows:
        wa = max(4, maxw * math.log10(a + 1) / math.log10(80001))
        wb = max(4, maxw * math.log10(b + 1) / math.log10(80001))
        s.append(f'<rect x="250" y="{y}" width="{wa:.0f}" height="20" rx="2" '
                 f'fill="{SEQ_WASH}"/>')
        s.append(f'<rect x="250" y="{y}" width="{wb:.0f}" height="20" rx="2" fill="{SEQ}"/>')
        s.append(f'<text x="242" y="{y+14}" class="n-t" fill="{INK}" '
                 f'text-anchor="end">{name}</text>')
        s.append(f'<text x="{250+max(wa,wb)+10}" y="{y+14}" class="n-v" fill="{INK}">'
                 f'{b:,}</text>')
        s.append(f'<text x="{250+max(wa,wb)+90}" y="{y+14}" class="n-s" fill="{INK3}">'
                 f'{note}</text>')
        y += 30
    s.append(f'<text x="250" y="{y+14}" class="n-s" fill="{INK3}">'
             'Pale bar: entering the stage. Solid bar: leaving it. Log scale.</text>')
    s.append(f'<rect x="0" y="{y+30}" width="{W}" height="56" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="16" y="{y+52}" class="n-t" fill="{INK}">'
             'About 660 CPU-hours per complex on eight cores, or roughly 33 GPU-hours.</text>')
    s.append(f'<text x="16" y="{y+72}" class="n-s" fill="{INK2}">'
             'Scaled from measured runtimes here: 0.4 s a sequence for ProteinMPNN, 22 s a '
             'construct for ESMFold, 148 s for a tethered three-chain complex.</text>')
    s.append("</svg>")
    funnel = "".join(s)

    cm = _confusion_matrix()
    return f"""
<section class="sec">
<div class="stack"><div class="eyebrow">Beyond the pipeline</div>
<h2>Taking it to the laboratory</h2>
<p>Everything above is computational. This section is separate from the pipeline because none of it
has been done, and none of it is part of the design. It sets out how large the campaign would have
to be and how the first measurements would be used. The counts come from published hit rates for
binder design, which for conditional binders are far lower than for ordinary binders, and the
compute is scaled from runtimes measured on this machine.</p></div>
<div class="chart wide">{funnel}</div>
<div class="note"><b>Filtering buys time, not success.</b> The RFdiffusion antibody work reports
that design success was no higher for filtered than for unfiltered designs. The in-silico stages
here are therefore sized to make the campaign affordable, not to raise the hit rate. What raises
the hit rate is the next part.</p></div>

<div class="stack"><h3>Three first, and what they are for</h3>
<p>The first order is three constructs, and its purpose is diagnosis rather than success. Each
tells you which stage of the pipeline is wrong:</p></div>
<div class="scroll"><table><thead><tr><th>measurement</th><th>reads out</th>
<th>if it fails, the problem is</th></tr></thead><tbody>
<tr><td>Purification — ÄKTA trace, single band on the gel</td><td>can it be made at all</td>
<td>the feasibility score: fold, aggregation, developability</td></tr>
<tr><td>Octet, toxin affinity with the protease absent</td><td>is it masked</td>
<td>the masking half of the unmasking score</td></tr>
<tr><td>Octet, toxin affinity with the protease present</td><td>does it activate</td>
<td>the trigger score, or the release half</td></tr>
</tbody></table></div>
<div class="note good"><b>Each failure mode points at a different module.</b> That is why three is
enough to start: a construct that will not express, one that expresses but never masks, and one
that masks but never releases are three different diagnoses, and each sends the work back to a
different part of the pipeline. Two or three cycles of this should give a pipeline that reliably
proposes three producible, toxin-binding, switching variants.</div>

<div class="stack"><h3>Then twenty to thirty, and the pipeline gets scored</h3>
<p>Once the loop reliably produces working candidates, the point of a larger batch changes. It is
no longer about finding a molecule; it is about measuring whether the scores predict anything. With
twenty to thirty measured constructs, each score can be scored in turn.</p></div>
<div class="chart wide">{cm}</div>
<div class="note warn"><b>The batch has to include designs the pipeline rejected.</b> If only
top-ranked candidates are made, every measurement lands in the left-hand column: TP and FP can be
counted, so precision is computable, but FN cannot, so recall is not. Spending perhaps a quarter of
the batch on designs the score rejected is what turns a hit list into a validation set.</div>
<div class="note"><b>Beyond that, the labels become training data.</b> At around a hundred measured
constructs the measurements stop being only a scorecard and start being a dataset: a binary label —
does it switch — against the sequence and structure representations already computed for every
candidate. That is the point at which a learned scorer can replace a hand-built one, and it is the
reason the representations are stored for designs that were never made.</div>
</section>"""


def _metrics(tp, fn, fp, tn):
    recall = tp / (tp + fn) if tp + fn else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return recall, precision, f1


#: Three outcomes from the same 24-construct batch. The numbers are illustrative;
#: the metrics beneath them are computed, so the arithmetic cannot drift.
VALIDATION_SCENARIOS = [
    ("Perfect scorer", 8, 0, 0, 16, "every switch found, nothing wasted"),
    ("Cautious scorer", 2, 6, 2, 14, "keeps few designs and misses most real switches"),
    ("Permissive scorer", 6, 2, 15, 1, "finds most real switches, but nearly everything "
                                       "it keeps is a dud"),
]


def _confusion_matrix() -> str:
    W, H = 1120, 360
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Validation outcomes">',
         '<style>'
         '.c-t{font:600 13px Inter,system-ui,sans-serif}'
         '.c-n{font:600 15px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.c-k{font:9.5px ui-monospace,monospace;letter-spacing:.06em}'
         '.c-s{font:11px Inter,system-ui,sans-serif}'
         '.c-m{font:600 11.5px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.c-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>',
         f'<text x="0" y="13" class="c-h" fill="{CHROME}">'
         'WHAT 24 MEASURED CONSTRUCTS SAY ABOUT A SCORE</text>']
    for n, (name, tp, fn, fp, tn, note) in enumerate(VALIDATION_SCENARIOS):
        x = n * 380
        recall, precision, f1 = _metrics(tp, fn, fp, tn)
        s.append(f'<text x="{x}" y="44" class="c-t" fill="{INK}">{name}</text>')
        s.append(f'<text x="{x}" y="62" class="c-s" fill="{INK3}">{note}</text>')
        s.append(f'<text x="{x+121}" y="84" class="c-k" fill="{INK3}" '
                 'text-anchor="middle">score keeps</text>')
        s.append(f'<text x="{x+211}" y="84" class="c-k" fill="{INK3}" '
                 'text-anchor="middle">score rejects</text>')
        for val, key, colour, col, row in ((tp, "TP", SEQ, 0, 0), (fn, "FN", STR, 1, 0),
                                           (fp, "FP", STR, 0, 1), (tn, "TN", SEQ, 1, 1)):
            cx, cy = x + 78, 92 + row * 64
            cx += col * 90
            good = key in ("TP", "TN")
            s.append(f'<rect x="{cx}" y="{cy}" width="86" height="60" rx="2" '
                     f'fill="{"#F2F7F6" if good else "#FCF4EC"}" stroke="{RULE}"/>')
            s.append(f'<text x="{cx+10}" y="{cy+19}" class="c-k" fill="{colour}">{key}</text>')
            s.append(f'<text x="{cx+10}" y="{cy+47}" class="c-n" fill="{INK}">{val}</text>')
        s.append(f'<text x="{x+72}" y="{126}" class="c-k" fill="{INK3}" '
                 'text-anchor="end">switches</text>')
        s.append(f'<text x="{x+72}" y="{190}" class="c-k" fill="{INK3}" '
                 'text-anchor="end">does not</text>')
        yy = 240
        for label, value in (("recall", recall), ("precision", precision), ("F1", f1)):
            bold = label == "F1"
            if bold:
                s.append(f'<line x1="{x+78}" y1="{yy-15}" x2="{x+254}" y2="{yy-15}" '
                         f'stroke="{RULE}"/>')
            s.append(f'<text x="{x+78}" y="{yy}" class="c-s" '
                     f'fill="{INK if bold else INK2}">{label}</text>')
            s.append(f'<text x="{x+254}" y="{yy}" class="c-m" '
                     f'fill="{INK if bold else INK2}" text-anchor="end">{value:.2f}</text>')
            yy += 20
    s.append(f'<rect x="0" y="300" width="{W}" height="54" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="16" y="322" class="c-t" fill="{INK}">'
             'F1 rises only when recall and precision are high together.</text>')
    s.append(f'<text x="16" y="342" class="c-s" fill="{INK2}">'
             'Both imperfect scorers look defensible on one metric and fail on F1. A score '
             'that keeps almost nothing and one that keeps almost everything are both useless, '
             'and a single headline number hides which failure you have.</text>')
    s.append("</svg>")
    return "".join(s)


#: Every claim on the page is one of these. The tag is shown beside the heading
#: so the reader never has to infer whether something was done or proposed.
STATUS = {
    "run": ("I ran this", "#0E7C7B"),
    "failed": ("I ran this and it did not work", "#C06014"),
    "proposed": ("Proposed, not done", "#6E7880"),
}


def status_tag(kind: str) -> str:
    label, colour = STATUS[kind]
    return (f'<span class="stat-tag" style="color:{colour};border-color:{colour}">'
            f'{label}</span>')


def status_table(rows) -> str:
    """One table, near the top, saying what exists and what does not.

    Each row is (what, status key, outcome in one sentence).
    """
    body = ""
    for what, kind, outcome in rows:
        label, colour = STATUS[kind]
        body += (f'<tr><td>{what}</td>'
                 f'<td><span class="stat-tag" style="color:{colour};'
                 f'border-color:{colour}">{label}</span></td>'
                 f'<td>{outcome}</td></tr>')
    return (f'<div class="scroll"><table class="wrap">'
            f'<colgroup><col style="width:27%"><col style="width:21%">'
            f'<col style="width:52%"></colgroup>'
            f'<thead><tr><th>step</th><th>status</th><th>what happened</th></tr></thead>'
            f'<tbody>{body}</tbody></table></div>')


#: Candidate methods for each module. Listing several makes the point that the
#: module is a slot, and that the pipeline does not depend on any one tool.
MODULE_METHODS = [
    ("Generate constructs",
     "Propose mask sequences and linkers inside the range the geometry allows.",
     [("Epitope graft", "take the antigen's own contact segment", "used here"),
      ("ProteinMPNN", "redesign a mask onto a fixed backbone", "used here"),
      ("RFdiffusion / RFantibody", "generate new mask backbones", "needs a GPU"),
      ("BindCraft", "one-shot binder design against the paratope", "needs a GPU"),
      ("Peptide language models", "sequence-only proposals for short masks", "not tried")]),
    ("Score: can it fold and be produced",
     "Reject constructs that will not express or will aggregate, before anything else.",
     [("ESMFold confidence", "whether a fold was predicted at all", "used here"),
      ("Core RMSD against the free nanobody", "whether the fusion disturbs it", "used here"),
      ("Exposed hydrophobic surface", "aggregation liability", "used here"),
      ("Unpaired cysteine count", "disulfide scrambling", "used here"),
      ("AlphaFold 3 / Boltz-2", "a second opinion on the fold", "GPU, or slow on CPU"),
      ("Rosetta or FoldX", "approximate stability", "not tried")]),
    ("Score: does the mask bind and release",
     "Two states. The site must be blocked while tethered and free once the tether is cut.",
     [("Worm-like-chain tether model", "the affinity range that can work", "used here"),
      ("Paratope occlusion from a predicted fold", "whether the site is covered",
       "used here, did not separate designs from controls"),
      ("Tethered three-chain complex (Boltz-2)", "whether the antigen still binds",
       "used here"),
      ("BioEmu or molecular dynamics", "an ensemble instead of one structure", "needs a GPU"),
      ("Constant-pH or free-energy methods", "quantitative affinity", "not tried")]),
    ("Score: can the linker be cleaved",
     "The protease must reach the linker and recognise it.",
     [("Linker accessibility and extension", "whether a protease can reach it", "used here"),
      ("CleaveNet", "generate linkers with a target cleavage profile",
       "trained on human MMPs"),
      ("ProsperousPlus / UniZyme", "cleavage-site prediction", "not tried"),
      ("Substrate profiling data", "the specificity the above models lack",
       "does not exist for these proteases")]),
    ("Combine into a ranking",
     "Turn three scores into one ordered list without letting them hide each other.",
     [("Pareto front", "keep a design unless another beats it on every score", "used here"),
      ("Hard filters first", "drop anything that cannot fold or cannot be cut", "used here"),
      ("Diversity selection", "keep structurally different designs, not one family",
       "not tried"),
      ("Uncertainty-aware ranking", "penalise designs the models disagree about",
       "not tried")]),
]


def module_tables() -> str:
    """One table per module, listing the methods that could fill it."""
    out = []
    for n, (title, blurb, methods) in enumerate(MODULE_METHODS):
        rows = "".join(
            f"<tr><td>{name}</td><td>{what}</td><td>{status}</td></tr>"
            for name, what, status in methods)
        out.append(
            f'<div class="modblock"><div class="modhead">'
            f'<span class="modnum">{n+1:02d}</span>'
            f'<div><b>{title}</b><span>{blurb}</span></div></div>'
            f'<div class="scroll"><table class="wrap">'
            f'<colgroup><col style="width:30%"><col style="width:43%">'
            f'<col style="width:27%"></colgroup>'
            f'<thead><tr><th>method</th><th>what it gives</th><th>status</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div></div>')
    return "".join(out)
