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
.vstage{position:relative;width:100%;height:460px;border:1px solid var(--rule-2);
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

function style(view, spec) {
  view.setStyle({}, {cartoon: {color: '#C9D1D6', opacity: 0.95}});
  if (spec.vhh) view.setStyle({resi: spec.vhh}, {cartoon: {color: '#2D5BD1'}});
  if (spec.toxin) view.setStyle({chain: spec.toxin}, {cartoon: {color: '#C06014'}});
  if (spec.mask) {
    view.setStyle({resi: spec.mask}, {cartoon: {color: '#0E7C7B'}});
    view.addStyle({resi: spec.mask}, {stick: {color: '#0E7C7B', radius: 0.17}});
  }
  if (spec.linker) view.setStyle({resi: spec.linker}, {cartoon: {color: '#8A9BA8'}});
  if (spec.paratope) view.addStyle({resi: spec.paratope},
      {stick: {color: '#C06014', radius: 0.19}});
}

function boot() {
  v1 = $3Dmol.createViewer('stage-solved', {backgroundColor: '#FFFFFF'});
  jQuery.get('structures/9GCN_AX.pdb', function (data) {
    v1.addModel(data, 'pdb');
    v1.setStyle({chain: 'A'}, {cartoon: {color: '#2D5BD1'}});
    v1.setStyle({chain: 'X'}, {cartoon: {color: '#C06014'}});
    v1.addStyle({chain: 'A', resi: PARATOPE}, {stick: {color: '#0E7C7B', radius: 0.19}});
    v1.zoomTo(); v1.zoom(1.3); v1.render();
  });

  v2 = $3Dmol.createViewer('stage-designed', {backgroundColor: '#FFFFFF'});
  loadDesigned(document.getElementById('pick').value);
  document.getElementById('pick').addEventListener('change', e => loadDesigned(e.target.value));
  document.querySelectorAll('[data-spin]').forEach(b => b.addEventListener('click', () => {
    const v = b.dataset.spin === 'solved' ? v1 : v2;
    v.spin(v.isSpinning && v.isSpinning() ? false : 'y');
  }));
}

function loadDesigned(id) {
  const m = MODELS[id];
  v2.clear();
  jQuery.get('structures/' + id + '.pdb', function (data) {
    v2.addModel(data, 'pdb');
    style(v2, {vhh: m.vhh, mask: m.mask, linker: m.linker, paratope: m.paratope});
    v2.zoomTo(); v2.zoom(1.3); v2.render();
    document.getElementById('dcap').textContent = m.caption;
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
