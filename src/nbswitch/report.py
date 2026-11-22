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
p{margin:0;max-width:none}
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
.lede{font-size:1.2rem;color:var(--ink-2);line-height:1.52;font-weight:400;max-width:none}
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
.chart img,.chart svg{width:100%;height:auto;display:block;max-width:none}
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
.tfig.full img{max-height:none;max-width:none}
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
.mri{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:0;
  border:1px solid var(--rule);border-radius:3px;overflow:hidden;margin:16px 0}
.mri>div{padding:15px 18px;border-right:1px solid var(--rule-2)}
.mri>div:last-child{border-right:0}
.mri .k{font-family:var(--mono);font-size:.6rem;letter-spacing:.12em;text-transform:uppercase;
  color:var(--ink-3);margin-bottom:6px}
.mri .v{font-size:.95rem;color:var(--ink);line-height:1.5}
.modnote{font-family:var(--mono);font-size:.64rem;letter-spacing:.1em;text-transform:uppercase;
  color:var(--chrome);margin-bottom:10px}
h3.mod{font-size:1.18rem;margin:38px 0 8px}
h3.mod span{font-family:var(--mono);font-size:.7rem;color:var(--ink-3);margin-right:12px}
.states{display:grid;grid-template-columns:repeat(3,1fr);gap:18px;margin:18px 0 10px}
.stfig{margin:0;background:var(--panel);border:1px solid var(--rule);border-radius:3px;
  padding:14px;display:flex;flex-direction:column;gap:10px}
.stfig img{width:100%;height:210px;object-fit:contain;display:block}
.stfig figcaption b{font-size:.95rem;display:block;color:var(--ink)}
.stfig figcaption span{font-size:.85rem;color:var(--ink-2);display:block;margin-top:3px;
  line-height:1.45}
.stkey{display:flex;gap:18px;flex-wrap:wrap;font-size:.82rem;color:var(--ink-3)}
.stkey i{display:inline-block;width:12px;height:12px;border-radius:2px;margin-right:6px;
  vertical-align:-1px}
@media(max-width:800px){.states{grid-template-columns:1fr}}
ul.ref{list-style:none;padding:0;margin:0;font-size:.88rem}
ul.ref li{padding:9px 0;border-bottom:1px solid var(--rule-2);color:var(--ink-2)}
ul.ref li:last-child{border-bottom:0}
footer{padding:36px 0 72px;color:var(--ink-3);font-size:.85rem}
@media(max-width:720px){.page{padding:0 20px}.vstage{height:320px}}
"""

VIEWER_JS = """
const MODELS = __MODELS__;
const PARATOPE = __PARATOPE__;
const SOLVED_FILE = __SOLVED__;
const ENSEMBLE = __ENSEMBLE__;
const CDRS = __CDRS__;
const TOP_PARATOPE = __TOPPARA__;
const SCORES = __SCORES__;
let v1 = null, v2 = null;

// Sticks are drawn with element colouring on a region-tinted carbon, so
// nitrogen, oxygen and sulfur read as themselves while the residue still
// says which part of the construct it belongs to.
const STICK_PARATOPE = {stick: {colorscheme: 'orangeCarbon', radius: 0.22}};
const STICK_MASK = {stick: {colorscheme: 'purpleCarbon', radius: 0.19}};

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
  jQuery.get(SOLVED_FILE, function (data) {
    const m = v1.addModel(data, 'pdb');
    v1.setStyle({chain: 'A'}, {cartoon: {colorfunc: function (atom) {
      return nanobodyColour(atom.resi, 0, '#337D79');
    }}});
    v1.setStyle({chain: 'B'}, {cartoon: {color: '#C99555'}});
    v1.addStyle({chain: 'A', resi: PARATOPE}, STICK_PARATOPE);
    (typeof TOP_PARATOPE !== 'undefined' ? TOP_PARATOPE : []).forEach(function (r) {
      const sel = m.selectedAtoms({chain: 'A', resi: r.resi, atom: 'CA'});
      if (!sel.length) return;
      const a = sel[0];
      v1.addLabel(r.label, {position: {x: a.x, y: a.y, z: a.z},
        backgroundColor: '#C06014', backgroundOpacity: 0.9, fontColor: 'white',
        fontSize: 10, borderThickness: 0, inFront: true});
    });
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
    // One setStyle over the whole chain. Restyling residue ranges one at a
    // time makes 3Dmol draw each range as its own cartoon, so the construct
    // renders as separate pieces even though it is a single chain.
    // ranges arrive as "a-b" strings
    const rangeEnd = function (r) {
      if (!r) return 0;
      const parts = String(r).split('-');
      return parseInt(parts[parts.length - 1], 10) || 0;
    };
    const maskEnd = rangeEnd(spec.mask);
    const linkEnd = rangeEnd(spec.linker) || maskEnd;
    const vhhStart = spec.vhh ? parseInt(String(spec.vhh).split('-')[0], 10) : 1;
    v2.setStyle({}, {cartoon: {colorfunc: function (atom) {
      if (maskEnd && atom.resi <= maskEnd) return '#8C72AD';
      if (linkEnd && atom.resi <= linkEnd) return '#9AA7B1';
      return nanobodyColour(atom.resi, vhhStart - 1, '#337D79');
    }}});
    // a thin backbone trace, so the linker reads as connected even where it is
    // coil and the cartoon renderer draws almost nothing
    v2.addStyle({atom: ['N', 'CA', 'C']}, {stick: {radius: 0.12, color: '#7E8A93'}});
    if (spec.mask) v2.addStyle({resi: spec.mask}, STICK_MASK);
    if (spec.paratope) v2.addStyle({resi: spec.paratope}, STICK_PARATOPE);
    termini(v2, m, {}, '');
    v2.zoomTo(); v2.zoom(1.15); v2.render();
    document.getElementById('dcap').textContent = spec.caption;
  });
}
// One aligned ensemble, mounted several times and coloured by a different
// score each time. Model 0 is the nanobody; every later model is one mask.
const RAMP = ['#2C5F8A', '#4E8FB0', '#8FBDC9', '#D8C9A3', '#D79A5B', '#C06014'];

function colourFor(value) {
  if (value === null || value === undefined || isNaN(value)) return '#B9C2C8';
  const i = Math.max(0, Math.min(RAMP.length - 1, Math.round(value * (RAMP.length - 1))));
  return RAMP[i];
}

const CDR_COLOUR = {CDR1: '#2E7DBE', CDR2: '#B8860B', CDR3: '#A8324A'};

// Which CDR a nanobody residue belongs to, built once from the numbering of
// the free nanobody. offset shifts construct numbering back onto it.
const CDR_OF = (function () {
  const m = {};
  if (typeof CDRS !== 'undefined') {
    Object.keys(CDRS).forEach(function (name) {
      CDRS[name].forEach(function (r) { m[r] = name; });
    });
  }
  return m;
})();

function nanobodyColour(resi, offset, base) {
  const name = CDR_OF[resi - (offset || 0)];
  return name ? CDR_COLOUR[name] : (base || '#337D79');
}

function mountEnsemble(hostId, scoreKey) {
  const host = document.getElementById(hostId);
  if (!host || typeof ENSEMBLE === 'undefined') return;
  const v = $3Dmol.createViewer(hostId, {backgroundColor: '#FFFFFF'});
  jQuery.get('structures/ensemble.pdb', function (data) {
    v.addModels(data, 'pdb');
    v.setStyle({model: 0}, {cartoon: {colorfunc: function (atom) {
      return nanobodyColour(atom.resi, 0, '#337D79');
    }}});
    v.addStyle({model: 0, resi: PARATOPE},
               {stick: {colorscheme: 'orangeCarbon', radius: 0.22}});
    ENSEMBLE.forEach(function (item, n) {
      const score = scoreKey ? (SCORES[item.id] || {})[scoreKey] : null;
      // the mask and linker carry the score colour; the short run of nanobody
      // past the junction is drawn as nanobody, so the chain reads continuous
      const end = item.linker_end;
      const c = colourFor(score);
      v.setStyle({model: n + 1}, {cartoon: {opacity: 0.85, thickness: 0.3,
        colorfunc: function (atom) {
          return atom.resi <= end ? c : nanobodyColour(atom.resi - end, 0, '#337D79');
        }}});
      v.addStyle({model: n + 1, atom: ['N', 'CA', 'C']},
                 {stick: {radius: 0.1, color: '#8A949B'}});
    });
    v.zoomTo(); v.zoom(1.55); v.render();
  });
}

window.addEventListener('load', function () {
  boot();
  document.querySelectorAll('[data-ensemble]').forEach(function (el) {
    mountEnsemble(el.id, el.dataset.ensemble || null);
  });
});
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
    wanted = set()
    for item in library:
        src = STRUCTURES / f"{item['id']}.pdb"
        if src.exists():
            shutil.copy(src, DOCS / "structures" / src.name)
            wanted.add(src.name)
    # docs/ mirrors the current library rather than accumulating. Renaming the
    # constructs left thirty models from an earlier run in place, and the
    # viewer offered them in its picker, where selecting one did nothing.
    for old in (DOCS / "structures").glob("*.pdb"):
        if old.name in wanted or old.stem == "ensemble" or old.stem[0].isdigit():
            continue
        old.unlink()
    return DOCS


def write_reference_pdb(chains=None):
    """One VHH-antigen pair from the deposited file, for the structure viewer."""
    import biotite.structure.io.pdb as pdb
    from .structure import read_atoms, fetch_structure
    from .config import REFERENCE_VHH_CHAIN, load_json
    if chains is None:
        toxin = load_json("01_reference.json")["chain_pairs"][REFERENCE_VHH_CHAIN]["toxin_chain"]
        chains = (REFERENCE_VHH_CHAIN, toxin)
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
    """The three states, drawn from coordinates rather than as blocks.

    All three panels are superimposed on the nanobody, so the antigen-binding
    surface faces the viewer identically. Panels A and B are the predicted
    construct; panel C is the deposited complex.
    """
    panels = [
        ("A", "a.png", "Intact construct",
         "The mask and linker are attached at the N-terminus of the nanobody."),
        ("B", "b.png", "After cleavage",
         "Proteolysis separates the mask, which then diffuses away."),
        ("C", "c.png", "Antigen bound",
         "The antigen engages the same surface, as in the deposited complex."),
    ]
    cells = "".join(
        f'<figure class="stfig"><img src="states/{img}" alt="{title}">'
        f'<figcaption><b>{tag} · {title}</b><span>{text}</span></figcaption></figure>'
        for tag, img, title, text in panels)
    return (f'<div class="states">{cells}</div>'
            f'<div class="stkey">'
            f'<span><i style="background:#337D79"></i>nanobody</span>'
            f'<span><i style="background:#8C72AD"></i>mask</span>'
            f'<span><i style="background:#9AA7B1"></i>linker</span>'
            f'<span><i style="background:#C06014"></i>antigen-binding residues</span>'
            f'<span><i style="background:#C99555"></i>antigen</span></div>')


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
    """Six modules: one estimation, one generation, three scores and a ranking."""
    W, H = 980, 838
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Computational pipeline">',
         '<defs><marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" '
         f'markerHeight="7" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{NULL}"/>'
         '</marker></defs>',
         '<style>'
         '.p-t{font:600 19px Inter,system-ui,sans-serif}'
         '.p-n{font:600 13px ui-monospace,monospace}'
         '.p-b{font:15px Inter,system-ui,sans-serif}'
         '.p-m{font:13.5px ui-monospace,monospace;letter-spacing:.02em}'
         '.p-io{font:600 14px ui-monospace,monospace;letter-spacing:.1em}'
         '.p-s{font:600 12px ui-monospace,monospace;letter-spacing:.09em}'
         '</style>']

    def box(x, y, w, h, accent, num, title, lines, methods, tag=None, tagcol=None):
        o = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="#FFFFFF" '
             f'stroke="{RULE}" stroke-width="1.2"/>',
             f'<rect x="{x}" y="{y}" width="{w}" height="4" fill="{accent}"/>',
             f'<text x="{x+20}" y="{y+34}" class="p-n" fill="{accent}">{num}</text>',
             f'<text x="{x+56}" y="{y+35}" class="p-t" fill="{INK}">{title}</text>']
        ty = y + 60
        for line in lines:
            o.append(f'<text x="{x+20}" y="{ty}" class="p-b" fill="{INK2}">{line}</text>')
            ty += 21
        ty += 4
        for line in methods:
            o.append(f'<text x="{x+20}" y="{ty}" class="p-m" fill="{accent}">{line}</text>')
            ty += 19
        if tag:
            o.append(f'<text x="{x+w-16}" y="{y+h-14}" class="p-s" fill="{tagcol}" '
                     f'text-anchor="end">{tag}</text>')
        return "".join(o)

    def band(y, text):
        return (f'<rect x="0" y="{y}" width="{W}" height="40" rx="4" fill="#F4F8F8" '
                f'stroke="{CHROME}" stroke-width="1"/>'
                f'<text x="{W/2}" y="{y+26}" class="p-io" fill="{CHROME}" '
                f'text-anchor="middle">{text}</text>')

    def shared(y, title, note, methods):
        """A step every downstream module depends on, drawn across the width."""
        return (f'<rect x="0" y="{y}" width="{W}" height="76" rx="4" fill="#FBF8F3" '
                f'stroke="{STR}" stroke-width="1.2" stroke-dasharray="5 3"/>'
                f'<text x="24" y="{y+30}" class="p-t" fill="{INK}">{title}</text>'
                f'<text x="24" y="{y+54}" class="p-b" fill="{INK2}">{note}</text>'
                f'<text x="{W-24}" y="{y+30}" class="p-m" fill="{STR}" '
                f'text-anchor="end">{methods}</text>'
                f'<text x="{W-24}" y="{y+54}" class="p-s" fill="{STR}" '
                f'text-anchor="end">RUN</text>')

    s.append(band(0, "INPUT · EXPERIMENTAL OR PREDICTED NANOBODY–ANTIGEN STRUCTURE"))
    s.append(f'<path d="M {W/2} 40 L {W/2} 58" stroke="{NULL}" stroke-width="1.6" '
             f'fill="none" marker-end="url(#ar)"/>')

    s.append(box(0, 64, 478, 148, CHROME, "01", "Usable mask affinity",
                 ["Interface geometry sets the attachment", "point, linker length and the range of",
                  "mask affinity that both states allow."],
                 ["interface analysis · tether model"], "RUN", CHROME))
    s.append(box(502, 64, 478, 148, SEQ, "02", "Generate constructs",
                 ["Propose mask sequences and cleavable", "linkers inside that range."],
                 ["RFdiffusion · ProteinMPNN · BindCraft",
                  "CleaveNet for linkers"], "PART RUN", STR))
    s.append(f'<path d="M 478 138 L 496 138" stroke="{NULL}" stroke-width="1.6" '
             f'fill="none" marker-end="url(#ar)"/>')
    # every construct is folded once, and all three scores read those coordinates
    s.append(f'<path d="M 741 212 L 741 228" stroke="{NULL}" stroke-width="1.6" '
             f'fill="none" marker-end="url(#ar)"/>')
    s.append(shared(232, "Structure prediction",
                    "Every construct folded once, as a single chain: mask, linker "
                    "and nanobody together.",
                    "ESMFold \u00b7 CPU"))
    s.append(f'<path d="M 490 308 L 490 324" stroke="{NULL}" stroke-width="1.6" '
             f'fill="none"/>')
    s.append(f'<path d="M 160 324 L 820 324" stroke="{NULL}" stroke-width="1.6" fill="none"/>')
    for x in (160, 490, 820):
        s.append(f'<path d="M {x} 324 L {x} 350" stroke="{NULL}" stroke-width="1.6" '
                 f'fill="none" marker-end="url(#ar)"/>')

    for n, (num, title, lines, methods, tag, col) in enumerate([
        ("03", "Structural feasibility",
         ["Is the fold kept, and can", "it be produced?"],
         ["pLDDT \u00b7 core RMSD", "surface \u00b7 cysteines"], "RUN", CHROME),
        ("04", "Mask binding and release",
         ["Is the site blocked, and", "freed after cleavage?"],
         ["Boltz-2 with alignments", "interface analysis"], "RUN", CHROME),
        ("05", "Linker cleavability",
         ["Is the linker recognised,", "and can a protease reach it?"],
         ["CleaveNet selectivity", "extension \u00b7 clearance"], "RUN", CHROME),
    ]):
        s.append(box(n * 330, 356, 320, 176, SEQ, num, title, lines, methods, tag, col))
        s.append(f'<text x="{n*330+300}" y="{550}" class="p-m" fill="{INK3}" '
                 f'text-anchor="end">S<tspan dy="4" font-size="11">{n+1}</tspan></text>')
        s.append(f'<path d="M {n*330+160} 532 L {n*330+160} 566 L 490 566 L 490 596" '
                 f'stroke="{NULL}" stroke-width="1.6" fill="none"'
                 + (' marker-end="url(#ar)"/>' if n == 1 else '/>'))

    s.append(box(166, 602, 648, 140, CHROME, "06", "Multiobjective ranking",
                 ["Hard constraints, then the nondominated set.",
                  "No weights, so no score compensates for another."],
                 ["Pareto front · uncertainty · diversity"], "RUN", CHROME))
    s.append(f'<path d="M 490 742 L 490 778" stroke="{NULL}" stroke-width="1.6" '
             f'fill="none" marker-end="url(#ar)"/>')
    s.append(band(784, "OUTPUT · RANKED CANDIDATE SEQUENCES WITH THEIR THREE SCORES"))
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
<figure class="tfig full"><img src="figures/paratope_footprint_9GCN.svg"
alt="Surface buried per residue for the set-aside complex">
<figcaption>Surface buried per nanobody residue for this complex. Almost all of it lies in CDR3,
and the N-terminal residue itself contributes, which is the geometry that defeats the
measurement.</figcaption></figure>
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
{f("fig4_1_3_a_isoforms", "Biolayer interferometry: tail peptides binding each isoform.", "full")}
{f("fig4_1_8_d", "Fluorescence polarisation, tail against FS12 at three ionic strengths. "
   "K<sub>D</sub> = 145 nM at 50 mM NaCl, 8.2 µM at 150 mM and 34.6 µM at 300 mM, indicating a "
   "largely electrostatic interaction.", "full")}

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
{f("fig4_1_14_f", "SEC-SAXS scattering with fits from each model. The two-state model gives the "
   "lowest chi-squared, 1.058, beating either single structure.", "full")}
{f("fig4_1_23_a_isoforms", "Conformational space explored for each isoform, radius of gyration "
   "against maximum dimension. The longer the tail, the more closed states appear.", "full")}

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
    "partial": ("I ran this, within a stated limit", "#8A6A1F"),
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


def mri(method: str, result: str, interpretation: str) -> str:
    """Method, result, interpretation. The order used throughout the thesis."""
    return (f'<div class="mri">'
            f'<div><div class="k">Method</div><div class="v">{method}</div></div>'
            f'<div><div class="k">Result</div><div class="v">{result}</div></div>'
            f'<div><div class="k">Interpretation</div><div class="v">{interpretation}</div></div>'
            f'</div>')


#: One compact table in place of five, as a module summary.
MODULE_SUMMARY = [
    ("01", "Estimation of the usable mask affinity",
     "Define the attachment point, linker length and the range of mask affinity "
     "that can satisfy both states",
     "Interface analysis, polymer models of the tether, competitive binding equilibria"),
    ("02", "Generation of candidate constructs",
     "Generate candidate masking architectures",
     "RFdiffusion, ProteinMPNN, BindCraft, peptide language models, evolutionary approaches"),
    ("03", "Structural feasibility",
     "Assess preservation of the VHH fold and structural plausibility",
     "ESMFold, AlphaFold, Boltz, BioEmu, Rosetta, FoldX"),
    ("04", "Mask binding and release",
     "Assess binding in the tethered and cleaved states",
     "Complex prediction, interface analysis, thermodynamic modelling"),
    ("05", "Linker cleavability",
     "Assess protease responsiveness and selectivity",
     "CleaveNet, MEROPS, substrate specificity models, accessibility analysis"),
    ("06", "Multiobjective ranking",
     "Select candidates satisfying the three objectives",
     "Hard filtering, Pareto optimisation, uncertainty and diversity-based selection"),
]


def module_summary_table() -> str:
    rows = "".join(
        f"<tr><td><b>{num}. {title}</b></td><td>{objective}</td><td>{methods}</td></tr>"
        for num, title, objective, methods in MODULE_SUMMARY)
    return ('<div class="scroll"><table class="wrap">'
            '<colgroup><col style="width:26%"><col style="width:30%">'
            '<col style="width:44%"></colgroup>'
            '<thead><tr><th>module</th><th>objective</th><th>methods</th></tr></thead>'
            f'<tbody>{rows}</tbody></table></div>')


def implementation_table(rows) -> str:
    """Implemented against remaining, one row per module."""
    body = "".join(f"<tr><td><b>{m}</b></td><td>{done}</td><td>{todo}</td></tr>"
                   for m, done, todo in rows)
    return ('<div class="scroll"><table class="wrap">'
            '<colgroup><col style="width:24%"><col style="width:38%">'
            '<col style="width:38%"></colgroup>'
            '<thead><tr><th>module</th><th>implemented</th><th>remaining work</th></tr></thead>'
            f'<tbody>{body}</tbody></table></div>')


def ensemble_viewer(element_id: str, score_key: str | None, caption: str,
                    legend: str) -> str:
    """A viewer over the aligned ensemble, coloured by one score."""
    attr = f' data-ensemble="{score_key}"' if score_key else ' data-ensemble=""'
    return (f'<div class="viewer">'
            f'<div class="vbar"><span class="eyebrow">All constructs, aligned on the nanobody'
            f'</span></div>'
            f'<div class="vstage" id="{element_id}"{attr}></div>'
            f'<div class="vkey">{legend}</div>'
            f'<p class="dim" style="font-size:.85rem">{caption}</p></div>')


CDR_KEY = ('<span><i style="background:#2E7DBE"></i>CDR1</span>'
           '<span><i style="background:#B8860B"></i>CDR2</span>'
           '<span><i style="background:#A8324A"></i>CDR3</span>')


def ramp_legend(low: str, high: str) -> str:
    stops = ["#2C5F8A", "#4E8FB0", "#8FBDC9", "#D8C9A3", "#D79A5B", "#C06014"]
    swatches = "".join(f'<i style="background:{c};width:18px;border-radius:0"></i>'
                       for c in stops)
    return (f'<span><i style="background:#337D79"></i>framework</span>{CDR_KEY}'
            f'<span><i style="background:#C06014"></i>antigen-binding residues</span>'
            f'<span style="display:flex;align-items:center;gap:0">{low}&nbsp;{swatches}'
            f'&nbsp;{high}</span>')


def interface_controls_table(controls, target: str) -> str:
    """The worked example read against independent VHH complexes."""
    rows = controls["comparators"]
    head = ("<tr><th>complex</th><th>chains</th><th>interface</th>"
            "<th>paratope</th><th>within 3 &#8491;</th><th>CDR3 share</th>"
            "<th>from loops</th></tr>")
    body = ""
    for r in rows:
        here = r["pdb_id"] == target
        mark = ' class="here"' if here else ""
        name = f"<b>{r['pdb_id']}</b>" if here else r["pdb_id"]
        if r.get("biparatopic"):
            ss = r["second_site"]
            name += (f' <span class="flag" title="second site {ss["chains"]}, '
                     f'{ss["interface_area_A2"]:.0f} &#8491;&#178;, '
                     f'{ss["shared_residues"]} residue shared">biparatopic</span>')
        body += (f'<tr{mark}><td>{name}</td>'
                 f'<td class="mono">{r["vhh_chain"]}:{r["antigen_chain"]}</td>'
                 f'<td>{r["interface_area_A2"]:.0f} &#8491;&#178;</td>'
                 f'<td>{r["n_paratope"]}</td><td>{r["n_contacts"]}</td>'
                 f'<td>{r["cdr3_percent_bsa"]:.0f}%</td>'
                 f'<td>{r["loop_percent_bsa"]:.0f}%</td></tr>')
    return f'<table class="tbl controls">{head}{body}</table>'


def generation_table(counts) -> str:
    """How many of what, at each stage. Every figure read from a stage output."""
    head = ("<tr><th>module</th><th>what is produced</th><th>how many</th>"
            "<th>kept</th></tr>")
    body = ""
    for row in counts:
        kept = row.get("kept")
        kept = "&mdash;" if kept is None else f"<b>{kept}</b>"
        body += (f'<tr><td class="mono">{row["stage"]}</td><td>{row["what"]}</td>'
                 f'<td>{row["made"]}</td><td>{kept}</td></tr>')
    return f'<table class="tbl gen">{head}{body}</table>'


def linker_table(linkers) -> str:
    """The generated substrates and their controls, with the contrast they scored."""
    head = ("<tr><th>linker</th><th>sequence</th><th>kind</th>"
            "<th>target</th><th>off&#8209;target</th><th>selectivity</th></tr>")
    body = ""
    for item in linkers:
        sel = item.get("selectivity")
        sel = "&mdash;" if sel is None else f"{sel:+.2f}"
        tz = item.get("target_z")
        oz = item.get("off_target_z")
        control = "control" in str(item.get("kind", ""))
        row_class = ' class="ctrl"' if control else ""
        body += (f'<tr{row_class}><td class="mono">{item["id"]}</td>'
                 f'<td class="mono">{item["sequence"]}</td>'
                 f'<td>{item.get("kind", "")}</td>'
                 f'<td>{"&mdash;" if tz is None else f"{tz:+.2f}"}</td>'
                 f'<td>{"&mdash;" if oz is None else f"{oz:+.2f}"}</td>'
                 f'<td><b>{sel}</b></td></tr>')
    return f'<table class="tbl linkers">{head}{body}</table>'
