"""Rendering the results page.

Values are read from the stage outputs, so nothing on the page is typed by hand.
"""
from __future__ import annotations
import hashlib
import shutil
from pathlib import Path

from .mask_design import (MASK_FAMILIES, MASK_LABELS, MASK_NOTES,
                          family_of, family_by_construct)
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
tr.warnrow td{background:#FDF4EC}
.await{border:1px dashed var(--rule);border-radius:3px;padding:13px 16px;margin:10px 0 18px;font-size:.88rem;color:var(--ink-3);background:#FAFBFB}
table.tbl{table-layout:fixed;width:100%}
table.tbl th,table.tbl td{white-space:normal;text-align:left;vertical-align:top;font-family:var(--sans);font-size:.88rem;line-height:1.45;padding:11px 13px;overflow-wrap:break-word}
table.tbl td.mono,table.tbl td:first-child{font-family:var(--mono);font-size:.78rem}
table.tbl thead th,table.tbl tr:first-child th{font-family:var(--mono);font-size:.6rem;letter-spacing:.08em;text-transform:uppercase;color:var(--ink-3);font-weight:600}
table.exp td{overflow-wrap:normal;word-break:normal}
.register{border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.modblock+.modblock{border-top:1px solid var(--rule)}
.modbar{display:flex;gap:18px;align-items:baseline;justify-content:space-between;background:#F4F8F8;padding:11px 14px}
.modname{font-weight:600;font-size:.95rem;color:var(--ink)}
.script{display:block;font-family:var(--mono);font-size:.66rem;color:var(--ink-3);font-weight:400;margin-top:3px}
.register table.exp{border:0}
.register table.exp tbody tr:last-child td{border-bottom:0}
.alts{font-family:var(--mono);font-size:.72rem;color:var(--ink-3);font-weight:400;text-align:right;flex:1 1 auto}
.rebuild{background:#FDF4EC;border:1px solid #C06014;border-left:4px solid #C06014;border-radius:3px;padding:14px 18px;margin:0 0 26px;font-size:.92rem;line-height:1.55;color:var(--ink-2)}
.pending{color:var(--ink-3);font-style:italic}
.tag{font-family:var(--mono);font-size:.6rem;letter-spacing:.06em;text-transform:uppercase;color:var(--chrome);border:1px solid var(--chrome);border-radius:2px;padding:0 5px;margin-left:6px;white-space:nowrap}
table.exp td{vertical-align:top}
td.dead{color:var(--chrome);font-weight:600}
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
.two{display:grid;grid-template-columns:1fr 1fr;gap:20px;align-items:stretch;margin:18px 0}
.two>*{margin:0;min-width:0}
.stackfig{display:flex;flex-direction:column;gap:14px;justify-content:center}
.stackfig figure{margin:0}
.two .viewer{display:flex;flex-direction:column}
.two .vstage{flex:1;height:auto;min-height:330px}
.two figure{display:flex;flex-direction:column;justify-content:center}
@media(max-width:860px){.two{grid-template-columns:1fr}.two .vstage{min-height:300px}}
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
.tfig.small img{max-height:230px;max-width:560px;margin:0 auto}
.tfig.medium img{max-height:470px;max-width:820px;margin:0 auto}
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
h4.subh{font-size:1.02rem;margin:32px 0 6px;color:var(--ink)}
h4.subh span{font-family:var(--mono);font-size:.64rem;letter-spacing:.1em;
  text-transform:uppercase;color:var(--chrome);margin-right:11px}
h5.subh2{font-size:.95rem;margin:30px 0 4px;color:var(--ink);font-weight:600}
.note{border:1px solid var(--rule);border-left:3px solid var(--chrome);border-radius:3px;
  background:var(--panel);padding:13px 18px;margin:16px 0;font-size:.92rem;color:var(--ink-2)}
.note b{color:var(--ink)}
.note.good{border-left-color:#1F6B4A}
.note.warn{border-left-color:#C06014}
.refmask{border:1px solid var(--rule);border-left:3px solid var(--chrome);border-radius:3px;
  background:var(--panel);padding:16px 20px;margin:18px 0}
.refmask .seq{font-family:var(--mono);font-size:1.32rem;letter-spacing:.07em;
  color:var(--ink);margin:6px 0 4px;word-break:break-all}
.refmask .sub{font-family:var(--mono);font-size:.76rem;color:var(--ink-3)}
.refmask p{margin:10px 0 0;font-size:.92rem;color:var(--ink-2)}
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
const TOXIN_CONTACT = __TOXCONTACT__;
const SCORES = __SCORES__;
const FINAL_PAIR = __FINALPAIR__;
let v1 = null, v2 = null, vRef = null, vMask = null;

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
    // the antigen's contacting face is shown by colouring its cartoon, not by
    // drawing side chains: what matters here is which part of the fold meets
    // the nanobody, and sticks on both partners obscure exactly that
    const toxContact = (typeof TOXIN_CONTACT !== 'undefined' ? TOXIN_CONTACT : []);
    v1.setStyle({chain: 'B'}, {cartoon: {colorfunc: function (atom) {
      return toxContact.indexOf(atom.resi) >= 0 ? '#B8541F' : '#E0C49A';
    }}});
    v1.addStyle({chain: 'A', resi: PARATOPE}, STICK_PARATOPE);
    termini(v1, m, {chain: 'A'}, '');
    termini(v1, m, {chain: 'X'}, '\u2032');
    v1.zoomTo(); v1.zoom(1.75); v1.render();
  });

  bootFinalPair();

  v2 = $3Dmol.createViewer('stage-designed', {backgroundColor: '#FFFFFF'});
  loadDesigned(document.getElementById('pick').value);
  document.getElementById('pick').addEventListener('change',
    e => loadDesigned(e.target.value));
  document.querySelectorAll('[data-spin]').forEach(b =>
    b.addEventListener('click', () => {
      const named = {solved: v1, designed: v2, finalref: vRef, finalmask: vMask};
      const v = named[b.dataset.spin] || v2;
      v.spin(v.isSpinning && v.isSpinning() ? false : 'y');
    }));
}

// The final screen: the same nanobody in two predictions, already
// superposed by stage 19b, so one camera serves both. Each viewer keeps its
// own camera object, so the zoom is applied identically rather than linked -
// linking them would let a drag on one silently reframe the other.
function bootFinalPair() {
  if (!FINAL_PAIR || !document.getElementById('stage-final-ref')) return;
  const paint = function (view, spec) {
    jQuery.get('structures/' + spec.file, function (data) {
      view.addModel(data, 'pdb');
      view.setStyle({chain: spec.vhh}, {cartoon: {colorfunc: function (atom) {
        return nanobodyColour(atom.resi, 0, '#337D79');
      }}});
      view.addStyle({chain: spec.vhh, resi: PARATOPE}, STICK_PARATOPE);
      view.setStyle({chain: spec.antigen}, {cartoon: {color: '#E0C49A',
                                                      opacity: 0.95}});
      if (spec.mask) {
        view.setStyle({chain: spec.mask}, {cartoon: {color: '#8C72AD'}});
        view.addStyle({chain: spec.mask}, STICK_MASK);
      }
      view.zoomTo({chain: spec.vhh});
      view.zoom(1.5);
      view.render();
    });
  };
  vRef = $3Dmol.createViewer('stage-final-ref', {backgroundColor: '#FFFFFF'});
  vMask = $3Dmol.createViewer('stage-final-mask', {backgroundColor: '#FFFFFF'});
  paint(vRef, FINAL_PAIR.reference);
  paint(vMask, FINAL_PAIR.masked);
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


def asset(path):
    """A served path with a short content hash on it.

    An <img src> that never changes is cached by the browser for as long as it
    likes, so a figure redrawn under the same name keeps showing the old
    picture while the page around it updates. The page looks published and is
    not, which is the same failure as a stale deploy and harder to see.
    Hashing the bytes means the URL changes exactly when the content does, and
    an unchanged figure still comes from cache.
    """
    local = DOCS / path
    if not local.exists():
        return path
    digest = hashlib.sha256(local.read_bytes()).hexdigest()[:8]
    return f"{path}?v={digest}"


def _chart(name, caption):
    # alt is empty on purpose. The caption sits next to the image and describes
    # it, so repeating it in alt made a screen reader read every caption twice
    # and put the same sentence twice into any copied selection. It also meant
    # a caption carrying markup or a quote leaked into the attribute.
    return (f'<figure class="chart"><img src="{asset(f"figures/{name}.svg")}" alt="">'
            f'<figcaption>{caption}</figcaption></figure>')


def copy_assets(library):
    """Copy figures and the structures the viewers need into docs/."""
    (DOCS / "figures").mkdir(exist_ok=True)
    (DOCS / "structures").mkdir(exist_ok=True)
    # docs/figures mirrors the results directory for the same reason docs/
    # structures mirrors the library: copying without pruning left five figures
    # published after the page stopped linking to them. Skip the prune if the
    # figure stage has not run in this checkout, so a page build cannot wipe
    # figures it did not generate.
    current = {svg.name for svg in FIGURES.glob("*.svg")}
    for svg in FIGURES.glob("*.svg"):
        shutil.copy(svg, DOCS / "figures" / svg.name)
    if current:
        for old in (DOCS / "figures").glob("*.svg"):
            if old.name not in current:
                old.unlink()
    wanted = set()
    for item in library:
        src = STRUCTURES / f"{item['id']}.pdb"
        if src.exists():
            shutil.copy(src, DOCS / "structures" / src.name)
            wanted.add(src.name)
    # docs/ mirrors the current library rather than accumulating. Renaming the
    # constructs left thirty models from an earlier run in place, and the
    # viewer offered them in its picker, where selecting one did nothing.
    # structures that belong to a stage rather than to the library: the
    # reference complex, the ensemble, and the final screen's superposed pair
    keep = {"ensemble", "final_reference", "final_masked"}
    for extra in ("final_reference.pdb", "final_masked.pdb"):
        src = STRUCTURES / extra
        if src.exists():
            shutil.copy(src, DOCS / "structures" / extra)
    for old in (DOCS / "structures").glob("*.pdb"):
        if old.name in wanted or old.stem in keep or old.stem[0].isdigit():
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
        f'<figure class="stfig"><img src="{asset(f"states/{img}")}" alt="">'
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
    """One estimation, two design stages, an assembly, three scores, a selection."""
    W, H = 980, 1196
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

    def box(x, y, w, h, accent, num, title, lines, methods, tag=None, tagcol=None,
            available=()):
        """One module. Methods that produced a number on this page are in the
        module's colour; methods the module accepts but that were not run are
        grey, so the difference is visible without a legend."""
        o = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="4" fill="#FFFFFF" '
             f'stroke="{RULE}" stroke-width="1.2"/>',
             f'<rect x="{x}" y="{y}" width="{w}" height="4" fill="{accent}"/>',
             f'<text x="{x+20}" y="{y+34}" class="p-n" fill="{accent}">{num}</text>',
             f'<text x="{x+62}" y="{y+35}" class="p-t" fill="{INK}">{title}</text>']
        ty = y + 60
        for line in lines:
            o.append(f'<text x="{x+20}" y="{ty}" class="p-b" fill="{INK2}">{line}</text>')
            ty += 21
        ty += 4
        for line in methods:
            o.append(f'<text x="{x+20}" y="{ty}" class="p-m" fill="{accent}">{line}</text>')
            ty += 19
        for line in available:
            o.append(f'<text x="{x+20}" y="{ty}" class="p-m" fill="{NULL}">{line}</text>')
            ty += 19
        if tag:
            o.append(f'<text x="{x+w-16}" y="{y+h-14}" class="p-s" fill="{tagcol}" '
                     f'text-anchor="end">{tag}</text>')
        return "".join(o)

    def band(y, text, h=40):
        return (f'<rect x="0" y="{y}" width="{W}" height="{h}" rx="4" fill="#F4F8F8" '
                f'stroke="{CHROME}" stroke-width="1"/>'
                f'<text x="{W/2}" y="{y+h/2+6}" class="p-io" fill="{CHROME}" '
                f'text-anchor="middle">{text}</text>')

    def shared(y, title, note, methods, h=76):
        """A step every downstream module depends on, drawn across the width."""
        return (f'<rect x="0" y="{y}" width="{W}" height="{h}" rx="4" fill="#FBF8F3" '
                f'stroke="{STR}" stroke-width="1.2" stroke-dasharray="5 3"/>'
                f'<text x="24" y="{y+30}" class="p-t" fill="{INK}">{title}</text>'
                f'<text x="24" y="{y+54}" class="p-b" fill="{INK2}">{note}</text>'
                f'<text x="{W-24}" y="{y+30}" class="p-m" fill="{STR}" '
                f'text-anchor="end">{methods}</text>'
                )

    def arrow(x1, y1, x2, y2):
        return (f'<path d="M {x1} {y1} L {x2} {y2}" stroke="{NULL}" stroke-width="1.6" '
                f'fill="none" marker-end="url(#ar)"/>')

    s.append(band(0, "INPUT · EXPERIMENTAL OR PREDICTED NANOBODY–ANTIGEN STRUCTURE"))
    s.append(arrow(W/2, 40, W/2, 56))

    # one estimation, which constrains both design stages below it
    s.append(box(0, 60, W, 150, CHROME, "01", "Usable mask affinity and linker length",
                 ["Interface geometry sets where to attach, how long the linker must be, and the",
                  "range of mask affinity that masking and release can both tolerate."],
                 ["interface analysis · worm-like chain · competition equilibria"],
                 available=["CALVADOS 2 · Martini 3 · Rosetta ddG · FoldX 5"]))

    # the two design stages are independent: neither constrains the other
    s.append(arrow(245, 210, 245, 234))
    s.append(arrow(735, 210, 735, 234))
    s.append(box(0, 238, 478, 196, SEQ, "02a", "Mask design",
                 ["What does the antigen present at the", "paratope, and how much of it to carry?"],
                 ["no model: contact face at 3.5 Å, window variants,",
                  "substitutions off the hotspots, matched controls",
                  "then ProteinMPNN v_48_020 on the folded backbone"],
                 available=["RFdiffusion · RFantibody · BindCraft",
                            "ESM-2 · ProGen2 · AfDesign"]))
    s.append(box(502, 238, 478, 196, SEQ, "02b", "Cleavable linker design",
                 ["Which sequences does the protease cut,",
                  "and how long must the linker be to reach?"],
                 ["CleaveNet, conditioned on a target profile",
                  "measured substrate library · flank length"],
                 available=["MEROPS · PROSPERous · DeepCleave",
                            "Procleave · iProt-Sub"]))

    # they only meet here
    s.append(f'<path d="M 239 434 L 239 458 L 490 458" stroke="{NULL}" stroke-width="1.6" fill="none"/>')
    s.append(f'<path d="M 741 434 L 741 458 L 490 458" stroke="{NULL}" stroke-width="1.6" fill="none"/>')
    s.append(arrow(490, 458, 490, 476))
    s.append(band(480, "ASSEMBLE · MASK + CLEAVABLE LINKER + NANOBODY, AS ONE CHAIN", 44))
    s.append(arrow(W/2, 524, W/2, 542))

    s.append(shared(546, "Structure prediction",
                    "Every construct folded once, as a single chain: mask, linker "
                    "and nanobody together.",
                    "ESMFold"))

    s.append(f'<path d="M 490 622 L 490 640" stroke="{NULL}" stroke-width="1.6" fill="none"/>')
    s.append(f'<path d="M 160 640 L 820 640" stroke="{NULL}" stroke-width="1.6" fill="none"/>')
    for x in (160, 490, 820):
        s.append(arrow(x, 640, x, 666))

    for n, (num, title, lines, methods, tag, col, avail) in enumerate([
        ("03", "Structural feasibility",
         ["Is the fold kept, and can", "it be produced?"],
         ["pLDDT · core RMSD", "surface · measured disulfides"], None, None,
         ["AlphaFold3 · Boltz-2 · Chai-1", "BioEmu-1 · Rosetta · FoldX 5"]),
        ("04", "Mask binding and release",
         ["Is the site occupied before", "cleavage, and free after?"],
         ["Boltz-2, mask alone and tethered", "interface analysis"], None, None,
         ["AlphaFold3 · Chai-1 · PRODIGY", "Rosetta · MM-GBSA · BioEmu-1"]),
        ("05", "Linker cleavability",
         ["Is the linker recognised,", "and can a protease reach it?"],
         ["CleaveNet · measured library", "extension · clearance"], None, None,
         ["MEROPS · PROSPERous · DeepCleave", "SitePrediction · MD of the site"]),
    ]):
        s.append(box(n * 330, 670, 320, 196, SEQ, num, title, lines, methods,
                     tag, col, available=avail))
        s.append(f'<text x="{n*330+300}" y="{884}" class="p-m" fill="{INK3}" '
                 f'text-anchor="end">S<tspan dy="4" font-size="11">{n+1}</tspan></text>')
        s.append(f'<path d="M {n*330+160} 866 L {n*330+160} 900 L 490 900 L 490 930" '
                 f'stroke="{NULL}" stroke-width="1.6" fill="none"'
                 + (' marker-end="url(#ar)"/>' if n == 1 else '/>'))

    # the alternatives listed here were NSGA-II, hypervolume and
    # uncertainty-based selection, none of which this pipeline runs, beside
    # hard thresholds and diversity selection, which it does. The methods line
    # now names what the stage actually does and the list is gone. The old
    # summary also said no score can compensate for another; that is what the
    # hard constraints do, not non-domination, which keeps a design that is
    # poor on one objective as long as nothing beats it on all three.
    s.append(box(166, 936, 648, 128, CHROME, "06", "Multiobjective selection",
                 ["Apply minimum requirements, take the nondominated",
                  "designs, and shortlist a diverse set of them."],
                 ["hard constraints · Pareto non-domination",
                  "maximin · one per distinct mask"]))
    s.append(arrow(490, 1064, 490, 1082))
    s.append(band(1086, "FINAL STRUCTURAL SCREEN · BOLTZ-2 ON THE MASKED "
                        "CONSTRUCT WITH ITS ANTIGEN", 44))
    s.append(arrow(490, 1130, 490, 1148))
    s.append(band(1152, "OUTPUT · CANDIDATES PRIORITISED FOR EXPERIMENTAL TESTING"))
    s.append("</svg>")
    return "".join(s)


def background_section() -> str:
    """Prior work, told as a sequence rather than a gallery.

    Follistatin's C-terminal tail is an endogenous tethered competitor for the
    same domain that binds heparan sulphate - the arrangement this pipeline
    designs on purpose - so the link is architectural rather than decorative.
    """
    def f(name, cap, cls=""):
        return (f'<figure class="tfig {cls}"><img src="{asset(f"thesis/{name}.png")}" alt="">'
                f'<figcaption>{cap}</figcaption></figure>')

    def step(n, title):
        return f'<div class="tstep"><b>{n}</b><span>{title}</span></div>'

    return f"""
<details class="fold">
<summary><b>Background in solution behaviour modelling of flexible multidomain signalling
inhibitor</b>
</summary>
<div class="foldbody">

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

{step("04", "Two states, resolved experimentally rather than predicted")}
{f("fst_saxs_abc", "From crystal structure to solution ensemble. (a) FST288 in the "
   "activin A complex, fitted to the SAXS envelope. (b) A model of FST315 adding the tails. "
   "(c) The solution structure is not one conformation but two — 65% open and 35% closed, "
   "interconverting, with the tail engaged in the closed state.", "portrait")}
{f("fig4_1_14_f", "SEC-SAXS scattering with fits from each model. The two-state model gives the "
   "lowest chi-squared, 1.058, beating either single structure.", "medium")}
{f("fig4_1_23_a_isoforms", "Conformational space explored for each isoform, radius of gyration "
   "against maximum dimension. The longer the tail, the more closed states appear.", "full")}

</div>
</details>"""




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


def validation_matrix() -> str:
    """Three outcomes from one batch of measured constructs.

    The counts are illustrative and the metrics under them are computed, so
    the arithmetic cannot drift from the cells above it. Its CSS classes are
    v- rather than c-: an inline <style> is not scoped to its own svg, and
    the construct schematic already owns c-.
    """
    W, H = 1120, 360
    s = [f'<svg viewBox="0 0 {W} {H}" width="100%" xmlns="http://www.w3.org/2000/svg" '
         f'role="img" aria-label="Validation outcomes">',
         '<style>'
         '.v-t{font:600 13px Inter,system-ui,sans-serif}'
         '.v-n{font:600 15px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.v-k{font:9.5px ui-monospace,monospace;letter-spacing:.06em}'
         '.v-s{font:11px Inter,system-ui,sans-serif}'
         '.v-m{font:600 11.5px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.v-h{font:600 10.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>',
         f'<text x="0" y="13" class="v-h" fill="{CHROME}">'
         'WHAT 24 MEASURED CONSTRUCTS SAY ABOUT A SCORE</text>']
    for n, (name, tp, fn, fp, tn, note) in enumerate(VALIDATION_SCENARIOS):
        x = n * 380
        recall, precision, f1 = _metrics(tp, fn, fp, tn)
        s.append(f'<text x="{x}" y="44" class="v-t" fill="{INK}">{name}</text>')
        s.append(f'<text x="{x}" y="62" class="v-s" fill="{INK3}">{note}</text>')
        s.append(f'<text x="{x+121}" y="84" class="v-k" fill="{INK3}" '
                 'text-anchor="middle">score keeps</text>')
        s.append(f'<text x="{x+211}" y="84" class="v-k" fill="{INK3}" '
                 'text-anchor="middle">score rejects</text>')
        for val, key, colour, col, row in ((tp, "TP", SEQ, 0, 0), (fn, "FN", STR, 1, 0),
                                           (fp, "FP", STR, 0, 1), (tn, "TN", SEQ, 1, 1)):
            cx, cy = x + 78, 92 + row * 64
            cx += col * 90
            good = key in ("TP", "TN")
            s.append(f'<rect x="{cx}" y="{cy}" width="86" height="60" rx="2" '
                     f'fill="{"#F2F7F6" if good else "#FCF4EC"}" stroke="{RULE}"/>')
            s.append(f'<text x="{cx+10}" y="{cy+19}" class="v-k" fill="{colour}">{key}</text>')
            s.append(f'<text x="{cx+10}" y="{cy+47}" class="v-n" fill="{INK}">{val}</text>')
        s.append(f'<text x="{x+72}" y="{126}" class="v-k" fill="{INK3}" '
                 'text-anchor="end">switches</text>')
        s.append(f'<text x="{x+72}" y="{190}" class="v-k" fill="{INK3}" '
                 'text-anchor="end">does not</text>')
        yy = 240
        for label, value in (("recall", recall), ("precision", precision), ("F1", f1)):
            bold = label == "F1"
            if bold:
                s.append(f'<line x1="{x+78}" y1="{yy-15}" x2="{x+254}" y2="{yy-15}" '
                         f'stroke="{RULE}"/>')
            s.append(f'<text x="{x+78}" y="{yy}" class="v-s" '
                     f'fill="{INK if bold else INK2}">{label}</text>')
            s.append(f'<text x="{x+254}" y="{yy}" class="v-m" '
                     f'fill="{INK if bold else INK2}" text-anchor="end">{value:.2f}</text>')
            yy += 20
    s.append(f'<rect x="0" y="300" width="{W}" height="54" rx="3" fill="#FFFFFF" '
             f'stroke="{RULE}"/>')
    s.append(f'<text x="16" y="322" class="v-t" fill="{INK}">'
             'F1 rises only when recall and precision are high together.</text>')
    s.append(f'<text x="16" y="342" class="v-s" fill="{INK2}">'
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


#: What each module could use, beyond what was run. Attached to the module and
#: shown once per group, since it does not vary by experiment.
MODULE_ALTERNATIVES = {
    "01": "CALVADOS 2; Martini 3; Rosetta ddG; FoldX 5",
    "02a": "RFdiffusion; RFantibody; BindCraft; ESM-2 650M; ProGen2; AfDesign",
    "02b": "MEROPS; PROSPERous; DeepCleave; Procleave; iProt-Sub",
    "03": "AlphaFold3; Boltz-2; Chai-1; BioEmu-1; Rosetta; FoldX 5; ThermoMPNN",
    "04": "AlphaFold3; Chai-1; Rosetta InterfaceAnalyzer; PRODIGY; MM-GBSA; BioEmu-1",
    "05": "MEROPS; PROSPERous; DeepCleave; SitePrediction; MD of the cleavage site",
    "06": "NSGA-II; hypervolume; hard thresholds; uncertainty- and diversity-based selection",
    "\u2014": "",
}

MODULE_TITLES = {
    "01": "Usable mask affinity",
    "02a": "Mask design",
    "02b": "Cleavable linker design",
    "03": "Structural feasibility",
    "04": "Mask binding and release",
    "05": "Linker cleavability",
    "06": "Multiobjective ranking",
    "\u2014": "Not yet implemented",
}


def pipeline_table(rows) -> str:
    """Every experiment, under the module it belongs to.

    One table per module with its heading above it. A single table with
    spanning header rows is the obvious construction and does not survive
    contact with a fixed column layout: the spanning cell collapses into the
    first column and the heading wraps one character per line.
    """
    label_of = {k: v[0] for k, v in EXP_STATUS.items()}
    colour_of = {k: v[1] for k, v in EXP_STATUS.items()}

    groups = []
    for r in rows:
        mod = r.get("module", "\u2014")
        if not groups or groups[-1][0] != mod:
            groups.append((mod, []))
        groups[-1][1].append(r)

    out = []
    for mod, items in groups:
        title = MODULE_TITLES.get(mod, "")
        alts = MODULE_ALTERNATIVES.get(mod, "")
        out.append(f'<div class="modblock"><div class="modbar">'
                   f'<span class="modname">{mod}. {title}</span>'
                   + (f'<span class="alts">also available: {alts}</span>' if alts else "")
                   + "</div>")
        body = ""
        for r in items:
            label = label_of.get(r["status"], "proposed")
            colour = colour_of.get(r["status"], "#6E7880")
            pending_row = r["status"] in ("not run", "proposed", "running")
            result = r.get("result") or ""
            if pending_row and not result:
                result = '<span class="pending">no data yet</span>'
            tag = (f'<span class="stat-tag" style="color:{colour};'
                   f'border-color:{colour}">{label}</span>')
            cls = ' class="ctrl"' if pending_row else ""
            script = r.get("script") or ""
            name_cell = (f'<b>{r["name"]}</b>'
                         + (f'<span class="script">{script}</span>' if script else ""))
            body += (f'<tr{cls}><td>{name_cell}</td><td>{r["question"]}</td>'
                     f'<td class="mono">{r.get("method", "")}</td>'
                     f'<td>{tag} {result}</td></tr>')
        out.append('<table class="tbl exp wrap">'
                   '<colgroup><col style="width:19%"><col style="width:24%">'
                   '<col style="width:26%"><col style="width:31%"></colgroup>'
                   f'<tbody>{body}</tbody></table></div>')
    return '<div class="register">' + "".join(out) + "</div>"

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


def pending(what: str) -> str:
    """A slot whose stage has not produced its output yet.

    The page is generated from stage outputs, and the expensive stages take
    hours. Rendering a labelled gap keeps the structure reviewable while they
    run, and makes a missing result visible as missing.
    """
    return (f'<div class="await"><b>Awaiting results.</b> {what} has not '
            f'produced an output yet. This slot fills from the stage output '
            f'when it does.</div>')


def interface_controls_table(controls, target: str) -> str:
    """The worked example read against independent VHH complexes."""
    if not controls:
        return pending("the independent-complex comparison, stage 15")
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
    """What each module produced, in the order the modules ran.

    One result column rather than an in/out pair. Most steps transform rather
    than select, and the units differ on each side of such a step, so an arrow
    between two counts would claim a funnel where there is none. Rows that do
    discard candidates say so in words and are marked.
    """
    head = ('<colgroup><col style="width:8%"><col style="width:46%">'
            '<col style="width:24%"><col style="width:22%"></colgroup>'
            "<tr><th>module</th><th>what happens</th><th>result</th>"
            "<th>feeds</th></tr>")
    body = ""
    for r in counts:
        feeds = r.get("feeds", "")
        dead = feeds == "nothing downstream"
        cell = f'<td class="dead">{feeds}</td>' if dead else f"<td>{feeds}</td>"
        tr = '<tr class="warnrow">' if dead else "<tr>"
        # a merit-based cut and a decision not to build a combination are
        # different things, and only the first discards a candidate
        sel = r.get("selected")
        label = sel if isinstance(sel, str) else ("selects" if sel else "")
        mark = f' <span class="tag">{label}</span>' if label else ""
        body += (f'{tr}<td class="mono">{r["stage"]}</td>'
                 f'<td>{r["what"]}{mark}</td>'
                 f'<td class="mono">{r["result"]}</td>{cell}</tr>')
    return f'<table class="tbl gen">{head}{body}</table>'


def linker_table(linkers) -> str:
    """The generated substrates and their controls, with the contrast they scored."""
    if not linkers:
        return pending("linker generation, stage 02b")
    head = ('<colgroup><col style="width:11%"><col style="width:34%">'
            '<col style="width:21%"><col style="width:11%">'
            '<col style="width:12%"><col style="width:11%"></colgroup>'
            "<tr><th>linker</th><th>sequence</th><th>kind</th>"
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


def assumptions_table(rows) -> str:
    """Every constant the design window rests on, with where it came from.

    The window is a model output, not a specification, and which parts of it
    are measured and which are chosen should be visible next to the number
    rather than in a methods note further down.
    """
    head = ('<colgroup><col style="width:21%"><col style="width:12%">'
            '<col style="width:30%"><col style="width:37%"></colgroup>'
            "<tr><th>quantity</th><th>value</th><th>where it comes from</th>"
            "<th>effect if it moves</th></tr>")
    body = ""
    for r in rows:
        cls = ' class="ctrl"' if r.get("kind") == "measured" else ""
        body += (f'<tr{cls}><td>{r["what"]}</td><td class="mono"><b>{r["value"]}</b></td>'
                 f'<td>{r["source"]}</td><td>{r["effect"]}</td></tr>')
    return f'<table class="tbl assume">{head}{body}</table>'


def affinity_scan_table(scan, used_M, closes_M) -> str:
    """The design window as a function of the nanobody's own antigen affinity."""
    head = ('<colgroup><col style="width:30%"><col style="width:35%">'
            '<col style="width:35%"></colgroup>'
            "<tr><th>V<sub>H</sub>H&ndash;antigen K<sub>D</sub></th>"
            "<th>release floor on the mask</th><th>usable window</th></tr>")
    body = ""
    for r in scan:
        kd = r["vhh_toxin_kd_M"] * 1e9
        here = abs(kd - used_M * 1e9) < 1e-6
        dec = r["window_decades"]
        shown = f"{dec:.2f} decades" if dec == dec and dec != float("-inf") else "closed"
        cls = ' class="hl"' if here else ""
        mark = " &larr; used here" if here else ""
        body += (f'<tr{cls}><td class="mono">{kd:.2f} nM{mark}</td>'
                 f'<td class="mono">{r["km_min_for_release_M"] * 1e9:,.0f} nM</td>'
                 f'<td class="mono">{shown}</td></tr>')
    return (f'<table class="tbl affin">{head}{body}</table>'
            f'<div class="note warn"><b>The architecture has a hard limit here.</b> '
            f'Above a V<sub>H</sub>H&ndash;antigen K<sub>D</sub> of about '
            f'<b>{closes_M * 1e9:.0f} nM</b> the window closes: the antigen can no longer '
            f'outcompete the released mask, and no mask affinity satisfies both '
            f'requirements. That is a property of the design, not of any particular mask.'
            f'</div>')


#: Status vocabulary for the experiment register. "run" means the output on
#: this page came from it; "not run" means the code exists and has not been
#: executed; "proposed" means it has not been written.
EXP_STATUS = {
    "run": ("run", "#0E7C7B"),
    "running": ("running", "#8A6A1F"),
    "not run": ("not run", "#6E7880"),
    "proposed": ("proposed", "#6E7880"),
    "failed": ("ran, did not work", "#C06014"),
}


def rebuild_banner(stale, n_total: int) -> str:
    """Say, at the top of the page, that a rebuild is in flight.

    Published mid-rebuild, the page mixes design stages that have been rerun
    with scores computed on the previous library. Every number is real; some
    of them answer a question about a construct set that no longer exists.
    """
    if not stale:
        return ""
    return (
        f'<div class="rebuild"><b>This page was built while a rebuild was '
        f'running.</b> The mask design, linker design and construct library have '
        f'been regenerated: masks are now taken from the antigen’s contact '
        f'face and linkers are 24 residues rather than 16. '
        f'{len(stale)} of {n_total} constructs are still being folded, so the '
        f'three scoring modules, the ranking and the Pareto front below were '
        f'computed on the previous library and do not yet describe these '
        f'constructs. The structural characterisation, the affinity window, the '
        f'linker generation, the length sweep and the untethered mask-binding '
        f'test are current.</div>')


def protonation_table(scan) -> str:
    """Histidine protonation in blood against a bite site, across its pKa range."""
    head = ('<colgroup><col style="width:17%"><col style="width:27%">'
            '<col style="width:27%"><col style="width:29%"></colgroup>'
            "<tr><th>histidine pK<sub>a</sub></th><th>protonated in blood</th>"
            "<th>protonated at the bite site</th><th>change</th></tr>")
    body = ""
    for r in scan:
        here = abs(r["pka"] - 6.0) < 1e-9
        cls = ' class="hl"' if here else ""
        mark = " &larr; free histidine" if here else ""
        body += (f'<tr{cls}><td class="mono">{r["pka"]:.1f}{mark}</td>'
                 f'<td class="mono">{r["blood_worst"]:.1%}</td>'
                 f'<td class="mono">{r["site_worst"]:.1%}</td>'
                 f'<td class="mono">{r["fold_change"]:.1f}&times;</td></tr>')
    return f'<table class="tbl proton">{head}{body}</table>'


def constraint_table(constraints, excluded_by) -> str:
    """The requirements applied before any comparison, and what each removed."""
    if not constraints:
        return pending("the ranking, stage 13")
    ops = {"ge": "&ge;", "le": "&le;", "gt": "&gt;", "lt": "&lt;"}
    head = ('<colgroup><col style="width:22%"><col style="width:12%">'
            '<col style="width:52%"><col style="width:14%"></colgroup>'
            "<tr><th>quantity</th><th>limit</th><th>why it is a requirement</th>"
            "<th>removed</th></tr>")
    body = ""
    for c in constraints:
        n = excluded_by.get(c["field"], 0)
        body += (f'<tr><td class="mono">{c["field"]}</td>'
                 f'<td class="mono">{ops.get(c["op"], c["op"])} {c["limit"]:g}</td>'
                 f'<td>{c["why"]}</td>'
                 f'<td class="mono">{n if n else "&mdash;"}</td></tr>')
    return f'<table class="tbl cons">{head}{body}</table>'


def front_table(front, objectives, shortlist=(), limit: int = 12) -> str:
    """The nondominated set in the stated order, with the shortlist marked."""
    if not front:
        return pending("the ranking, stage 13")
    label = {"feasibility": "feasibility", "cleavability": "cleavability",
             "surface_occlusion": "surface occlusion"}
    cols = "".join(f"<th>{label.get(o, o)}</th>" for o in objectives)
    head = ('<colgroup><col style="width:30%"><col style="width:14%">'
            + '<col style="width:18%">' * len(objectives) + "</colgroup>"
            f"<tr><th>construct</th><th>weakest, normalised</th>{cols}</tr>")
    body = ""
    for r in front[:limit]:
        picked = r["id"] in shortlist
        cls = ' class="hl"' if picked else ""
        mark = ' <span class="tag">verified</span>' if picked else ""
        cells = "".join(f'<td class="mono">{r.get(o, float("nan")):.3f}</td>'
                        for o in objectives)
        body += (f'<tr{cls}><td class="mono">{r["id"]}{mark}</td>'
                 f'<td class="mono"><b>{r.get("robust_score", float("nan")):.3f}</b></td>'
                 f"{cells}</tr>")
    more = ""
    if len(front) > limit:
        more = (f'<p class="dim" style="font-size:.85rem">'
                f'{len(front) - limit} further nondominated designs not shown.</p>')
    return f'<table class="tbl front">{head}{body}</table>{more}'


def alignment_control(without: float, with_msa: float) -> str:
    """The positive control that decided how the complexes are predicted.

    Before asking whether a mask blocks the antigen, the method has to
    reproduce the interface the crystal structure already shows. It is a
    positive control, and the first attempt failed it.
    """
    if not (without or with_msa):
        return pending("the alignment control, stages 09 and 19")
    return (
        f'<p><b>The method was checked against a known answer first.</b> Before '
        f'assessing the effect of the mask, I tested whether Boltz-2 could '
        f'reproduce the experimentally determined V<sub>H</sub>H&ndash;antigen '
        f'interface. Without evolutionary alignments, only {without:.0f}% of '
        f'predicted antigen contacts involved the known paratope. Supplying '
        f'alignments for the V<sub>H</sub>H and the antigen, which are public '
        f'sequences, raised this to {with_msa:.0f}%, while the designed mask took '
        f'no alignment and never left this machine.</p>'
        f'<p>I therefore used the revised setup for the rest of the analysis. Two '
        f'qualifications. That figure is recovery of the reference interface, not '
        f'a probability that the proteins bind, and reproducing a known interface '
        f'does not establish that predicted mask binding or release is accurate. '
        f'For the control the comparison is clean, since both runs modelled the '
        f'same two chains and differed only in the alignment; for the masked '
        f'constructs it is not, because the revised setup also splits the fusion '
        f'into separate chains joined by a bond, so the improvement there cannot '
        f'be attributed to the alignment alone.</p>')

_WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven",
          "eight", "nine", "ten", "eleven", "twelve")


def word(n):
    """Small counts spelled out, so prose reads as prose.

    The count still comes from the data; only its rendering changes, which is
    the difference between "nine other MMPs" and a literal nine in the text.
    """
    return _WORDS[n] if 0 <= n < len(_WORDS) else f"{n:,}"


def _mmp_names(names):
    """['MMP2','MMP9','MMP8'] -> 'MMP-2, -8, -9', in numeric order."""
    nums = sorted(int(n.replace("MMP", "")) for n in names)
    return "MMP-" + ", -".join(str(n) for n in nums)


def _wrap(text, width):
    out, line = [], ""
    for word in text.split(" "):
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    if line:
        out.append(line)
    return out


def mask_blocks(name, m, masks, core="epi_core", face="epi_face2"):
    """Where a mask's residues sit on the antigen, as (position, text) blocks.

    Only the grafts cut straight out of the antigen carry coordinates. The
    substitutions and controls are variants of the main loop, so they are
    drawn under it; the mask that joins both regions is drawn as the two
    regions it joins, with its spacer in the gap, because drawing it as one
    run would put its second half under the wrong residues.
    """
    if m.get("start") is not None:
        return [(m["start"], m["seq"])]
    if name == "epi_both":
        a, b = masks[face], masks[core]
        spacer = m["seq"][len(a["seq"]):len(m["seq"]) - len(b["seq"])]
        if a["seq"] + spacer + b["seq"] == m["seq"]:
            return [(a["start"], a["seq"]), (None, spacer), (b["start"], b["seq"])]
    return [(masks[core]["start"], m["seq"])]


def construct_schematic(c, free_vhh_len=None):
    """One assembled construct, end to end, drawn to scale.

    The arithmetic above says how many constructs there are; this says what
    one of them is. Drawn to scale because the proportions are the point: the
    mask and linker together are a short N-terminal addition to a nanobody
    that is most of the molecule, and the protease has to find ten residues
    inside that.

    Every label sits outside the bar. Drawn inside, the substrate band covered
    the linker's own name and range, because the band is in the middle of the
    linker and so is its label.
    """
    MASKC, FLANKC, CUTC, VHHC = "#C06014", "#9AA6AE", "#2D5BD1", "#337D79"
    W, bh = 980, 40
    m_lo, m_hi = c["mask_range"]
    l_lo, l_hi = c["linker_range"]
    v_lo, v_hi = c["vhh_range"]
    total = v_hi
    px = lambda n: n * (W - 2) / total

    o = [None,
         '<style>'
         '.c-h{font:600 13px Inter,system-ui,sans-serif}'
         '.c-b{font:12px Inter,system-ui,sans-serif}'
         '.c-m{font:600 11.5px ui-monospace,monospace}'
         '.c-n{font:10.5px ui-monospace,monospace}'
         '</style>']

    def label(x, y, text, cls="c-b", fill=INK2, anchor="start"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}" fill="{fill}" '
                 f'text-anchor="{anchor}">{text}</text>')

    label(0, 13, "N", "c-m", INK)
    label(W, 13, "C", "c-m", INK, "end")
    o.append(f'<line x1="13" y1="9" x2="{W - 13}" y2="9" stroke="{INK3}" '
             f'stroke-width="1" stroke-opacity="0.35"/>')

    by = 42
    segs = [("mask", MASKC, m_lo, m_hi), ("cleavable linker", FLANKC, l_lo, l_hi),
            ("nanobody", VHHC, v_lo, v_hi)]
    for name, colour, lo, hi in segs:
        x, w = px(lo - 1), px(hi - lo + 1)
        o.append(f'<rect x="{x}" y="{by}" width="{w}" height="{bh}" rx="3" '
                 f'fill="{colour}" fill-opacity="0.18" stroke="{colour}" '
                 f'stroke-width="1.4"/>')
        label(x + w / 2, by - 9, name, "c-h", colour, "middle")
        label(x + w / 2, by + bh + 16, f"{lo}\u2013{hi}, {hi - lo + 1} aa",
              "c-n", INK3, "middle")

    # the ten residues the protease reads, inside the linker and unlabelled:
    # its callout hangs below the bar so nothing is drawn over anything
    sub = c.get("linker_substrate") or ""
    lseq = c.get("linker_seq") or ""
    if sub and sub in lseq:
        off = lseq.index(sub)
        s_lo = l_lo + off
        x, w = px(s_lo - 1), px(len(sub))
        o.append(f'<rect x="{x}" y="{by + 1.5}" width="{w}" height="{bh - 3}" rx="2" '
                 f'fill="{CUTC}" fill-opacity="0.85"/>')
        cx = x + w / 2
        o.append(f'<path d="M {cx} {by + bh + 24} L {cx} {by + bh + 34} '
                 f'L {cx + 26} {by + bh + 34}" fill="none" stroke="{CUTC}" '
                 f'stroke-width="1.1" stroke-dasharray="3 2"/>')
        label(cx + 32, by + bh + 38,
              f"protease site, {len(sub)} residues", "c-n", CUTC)

    y = by + bh + 70
    rows = [("mask", c.get("mask", ""), c.get("mask_seq", ""), MASKC),
            ("linker", c.get("linker_id", ""), lseq, FLANKC)]
    for name, ident, seq, colour in rows:
        o.append(f'<rect x="0" y="{y - 11}" width="3" height="14" rx="1.5" fill="{colour}"/>')
        label(12, y, name, "c-b", INK)
        label(78, y, ident, "c-n", INK3)
        label(206, y, seq, "c-m", INK)
        y += 22
    o.append(f'<rect x="0" y="{y - 11}" width="3" height="14" rx="1.5" fill="{VHHC}"/>')
    label(12, y, "nanobody", "c-b", INK)
    label(206, y,
          f"{v_hi - v_lo + 1} residues, unchanged in every construct", "c-n", INK3)
    y += 26
    label(0, y, f'<tspan class="c-m">{c["id"]}</tspan>  '
                f'\u00b7  {total} residues, folded as one chain', "c-b", INK3)
    H = int(y + 10)
    o[0] = (f'<svg viewBox="0 0 {W} {H}" width="100%" '
            'xmlns="http://www.w3.org/2000/svg" role="img" '
            'aria-label="One assembled construct, drawn to scale">')
    o.append("</svg>")
    return "".join(o)

def selection_funnel(f) -> str:
    """How 348 scored designs became one worth making.

    The three pairwise score plots showed the front; they did not show the
    narrowing, which is what the section is about. The final structural screen
    is drawn apart from the ranking because it answers a different question:
    the ranking says which designs are worth examining, Boltz-2 says whether
    the mask actually blocks the antigen.
    """
    W = 980
    ACC, OK = "#2D5BD1", "#1F6B4A"
    o = [None,
         '<defs><marker id="sf" viewBox="0 0 10 10" refX="9" refY="5" '
         'markerWidth="6" markerHeight="6" orient="auto">'
         f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{INK3}"/></marker></defs>',
         '<style>'
         '.s-n{font:600 30px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.s-t{font:600 13px Inter,system-ui,sans-serif}'
         '.s-s{font:11px Inter,system-ui,sans-serif}'
         '.s-h{font:600 9.5px ui-monospace,monospace;letter-spacing:.14em}'
         '</style>']

    def label(x, y, text, cls="s-s", fill=INK2, anchor="middle"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}" fill="{fill}" '
                 f'text-anchor="{anchor}">{text}</text>')

    def step(x, w, y, num, title, sub, accent=INK, wash="#FFFFFF"):
        o.append(f'<rect x="{x}" y="{y}" width="{w}" height="104" rx="6" '
                 f'fill="{wash}" stroke="{RULE}" stroke-width="1.2"/>')
        cx = x + w / 2
        label(cx, y + 42, f"{num:,}", "s-n", accent)
        label(cx, y + 62, title, "s-t", INK)
        label(cx, y + 80, sub, "s-s", INK3)

    y = 22
    label(0, 11, "RANKING", "s-h", CHROME, "start")
    # the screen card carries the longest sub-label, so it gets the width:
    # at 236 a side the three ranking cards left it 158 px and it overran
    w, gap = 200.0, 26.0
    cells = [(f["n_evaluated"], "scored constructs",
              "all three scores", INK),
             (f["n_front"], "Pareto candidates",
              f"{f['n_evaluated'] - f['n_front']} beaten on all three", INK),
             (f["n_shortlist"], "shortlisted", "maximin, one per mask", INK)]
    for i, (num, title, sub, col) in enumerate(cells):
        x = i * (w + gap)
        step(x, w, y, num, title, sub, col)
        if i:
            o.append(f'<line x1="{x - gap + 6}" y1="{y + 52}" x2="{x - 6}" '
                     f'y2="{y + 52}" stroke="{INK3}" stroke-width="1.2" '
                     f'marker-end="url(#sf)"/>')

    # the screen is a separate question, so it is drawn across a divider
    xd = 3 * w + 2 * gap + 44
    o.append(f'<line x1="{xd - 22}" y1="{y - 14}" x2="{xd - 22}" y2="{y + 118}" '
             f'stroke="{RULE}" stroke-width="1.2" stroke-dasharray="4 4"/>')
    label(xd, 11, "FINAL STRUCTURAL SCREEN · BOLTZ-2", "s-h", CHROME, "start")
    step(xd, W - xd, y, f["n_worth"], "prioritised for testing",
         "mask on the paratope, antigen displaced", OK, "#F3F8F5")
    o.append(f'<line x1="{xd - 38}" y1="{y + 52}" x2="{xd - 6}" y2="{y + 52}" '
             f'stroke="{INK3}" stroke-width="1.2" marker-end="url(#sf)"/>')

    H = int(y + 104 + 6)
    o[0] = (f'<svg viewBox="0 0 {W} {H}" width="100%" '
            'xmlns="http://www.w3.org/2000/svg" role="img" '
            'aria-label="From scored constructs to one worth making">')
    o.append("</svg>")
    return "".join(o)


def pipeline_funnel(f) -> str:
    """What the worked example generated, evaluated and selected.

    Three rows, one per phase, each a line of counts with the phase named in
    the left gutter. The first version boxed every count and ran to 708 px of
    height for twelve numbers; this fits the same numbers in a third of that,
    which is the right size for something read once on the way past.
    """
    W = 980
    ACC, OK = "#2D5BD1", "#1F6B4A"
    GUT, X0 = 104, 124
    span = W - X0

    o = [None,
         '<defs><marker id="fn" viewBox="0 0 10 10" refX="9" refY="5" '
         'markerWidth="6" markerHeight="6" orient="auto">'
         f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{INK3}"/></marker></defs>',
         '<style>'
         '.f-p{font:600 9.5px ui-monospace,monospace;letter-spacing:.14em}'
         '.f-n{font:600 23px ui-monospace,monospace;font-variant-numeric:tabular-nums}'
         '.f-t{font:600 11.5px Inter,system-ui,sans-serif}'
         '.f-s{font:10px Inter,system-ui,sans-serif}'
         '.f-o{font:9.5px ui-monospace,monospace}'
         '</style>']

    def label(x, y, text, cls="f-s", fill=INK2, anchor="middle"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}" fill="{fill}" '
                 f'text-anchor="{anchor}">{text}</text>')

    def item(cx, y, num, title, sub="", over="", accent=INK):
        if over:
            label(cx, y - 12, over, "f-o", INK3)
        label(cx, y + 18, f"{num:,}", "f-n", accent)
        label(cx, y + 34, title, "f-t", INK)
        if sub:
            label(cx, y + 48, sub, "f-s", INK3)

    def link(x0, x1, y, glyph=""):
        mid = (x0 + x1) / 2
        if glyph:
            label(mid, y + 4, glyph, "f-t", INK3)
        else:
            o.append(f'<line x1="{x0 + 6}" y1="{y}" x2="{x1 - 6}" y2="{y}" '
                     f'stroke="{INK3}" stroke-width="1.1" marker-end="url(#fn)"/>')

    def row(y, h, tag):
        o.append(f'<rect x="0" y="{y}" width="{W}" height="{h}" rx="5" '
                 f'fill="#FBFBFA" stroke="{RULE}" stroke-width="1"/>')
        label(14, y + h / 2 + 3, tag, "f-p", CHROME, "start")

    # --- 1 generate: two inputs, their product, and what was folded --------
    r1, h1 = 0, 92
    row(r1, h1, "1 GENERATE")
    cells = [(f["n_masks"], "mask sequences", "grafts, variants, controls",
              "antigen interface", ACC),
             (f["n_linkers"], "linker variants",
              f"{f['n_motifs']} motifs + {f['n_linker_controls']} controls",
              f"CleaveNet \u00b7 {f['n_generated']:,} peptides", ACC),
             (f["n_pairs"], "mask\u2013linker pairs", "every mask, every linker", "", INK),
             (f["n_structures"], "predicted structures", "+ the free nanobody", "", INK)]
    cw = span / len(cells)
    for i, (num, title, sub, over, col) in enumerate(cells):
        cx = X0 + i * cw + cw / 2
        item(cx, r1 + 26, num, title, sub, over, col)
        if i:
            link(X0 + (i - 1) * cw + cw / 2 + 62, cx - 62, r1 + 44,
                 "\u00d7" if i == 1 else "")

    # --- 2 evaluate: one set, three scores, all of them on all of it -------
    r2, h2 = h1 + 10, 78
    row(r2, h2, "2 EVALUATE")
    item(X0 + cw / 2, r2 + 20, f["n_evaluated"], "candidate constructs",
         "the free nanobody is the reference")
    chips = ("structural feasibility", "mask binding", "linker cleavability")
    cx0 = X0 + cw + 24
    chw = (W - 16 - cx0) / len(chips)
    for i, name in enumerate(chips):
        x = cx0 + i * chw
        o.append(f'<rect x="{x}" y="{r2 + 18}" width="{chw - 12}" height="42" rx="4" '
                 f'fill="#FFFFFF" stroke="{RULE}" stroke-width="1"/>')
        label(x + (chw - 12) / 2, r2 + 36, name, "f-t", ACC)
        label(x + (chw - 12) / 2, r2 + 51,
              f"{f['n_evaluated']} of {f['n_evaluated']}", "f-o", INK3)
    link(X0 + cw / 2 + 62, cx0 - 6, r2 + 38)

    # --- 3 select: three narrowings ---------------------------------------
    r3, h3 = r2 + h2 + 10, 86
    row(r3, h3, "3 SELECT")
    steps = [(f["n_front"], "nondominated", "Pareto, no weights", INK),
             (f["n_shortlist"], "shortlisted", "one per distinct mask", INK),
             (f["n_worth"], "worth making", "verified against the antigen", OK)]
    sw = span / len(steps)
    for i, (num, title, sub, col) in enumerate(steps):
        cx = X0 + i * sw + sw / 2
        item(cx, r3 + 20, num, title, sub, "", col)
        if i:
            link(X0 + (i - 1) * sw + sw / 2 + 74, cx - 74, r3 + 38)

    H = int(r3 + h3 + 2)
    o[0] = (f'<svg viewBox="0 0 {W} {H}" width="100%" '
            'xmlns="http://www.w3.org/2000/svg" role="img" '
            'aria-label="What the worked example generated, evaluated and selected">')
    o.append("</svg>")
    return "".join(o)

def mask_families(masks, hotspots, core="epi_core", families=MASK_FAMILIES):
    """Every mask, grouped by the question it was built to answer.

    A list of twelve sequences says what was made, not why. Grouped by
    question and drawn against the loop they came from, each row says which
    variable it moves. The cysteines are coloured because they are the reason
    a longer graft is not automatically a better one: past C42-C53 the
    extensions pick up cysteines whose partners stay behind in the antigen.
    """
    # Kept and orphaned cysteines are a green/red pair, which the palette
    # validator puts at dE 4.5 under deuteranopia: hue alone cannot carry this
    # distinction, and it is the one the figure exists to make. So it is drawn
    # as shape first - a closed bracket where the partner is inside the mask,
    # an open dangling stub where it is not - with colour only reinforcing it.
    # The hotspot tick is blue rather than the teal used elsewhere, which sat
    # at dE 7.9 against the kept-cysteine green even in normal vision.
    KEPT, ORPH, SUBST, PLAIN = "#1F6B4A", "#9C3B3B", "#C06014", "#8C979E"
    HOT = "#2D5BD1"
    core_seq = masks[core]["seq"]
    core_start, core_end = masks[core]["start"], masks[core]["end"]
    blocks = {n: mask_blocks(n, masks[n], masks, core)
              for _, names in families for n in names}

    # the two contact regions are far apart on the antigen, so the ruler skips
    # the stretch between them rather than spending half its width on residues
    # no mask uses
    used = sorted({pos + i for bs in blocks.values() for pos, txt in bs
                   if pos is not None for i in range(len(txt))})
    runs = []
    for pos in used:
        if runs and pos == runs[-1][1] + 1:
            runs[-1][1] = pos
        else:
            runs.append([pos, pos])
    col, c = {}, 0
    for i, (lo, hi) in enumerate(runs):
        if i:
            c += 1                                   # one column for the break
        for pos in range(lo, hi + 1):
            col[pos] = c
            c += 1
    ncol = c
    # x0 clears the longest row label; 'Conservative substitutions III'
    # ran into the first column at 206
    cell, x0, notew = 19.0, 228.0, 250.0
    W = x0 + ncol * cell + notew

    o = [None,
         '<style>'
         '.m-h{font:600 13px Inter,system-ui,sans-serif}'
         '.m-b{font:12px Inter,system-ui,sans-serif}'
         '.m-q{font:600 11px ui-monospace,monospace;letter-spacing:.09em}'
         '.m-t{font:600 12px ui-monospace,monospace}'
         '.m-n{font:10.5px ui-monospace,monospace}'
         '</style>']

    def label(x, y, text, cls="m-b", fill=INK2, anchor="start"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}" fill="{fill}" '
                 f'text-anchor="{anchor}">{text}</text>')

    y = 15
    for lo, hi in runs:
        for pos in (lo, hi):
            label(x0 + col[pos] * cell + cell / 2, y, str(pos), "m-n", INK3, "middle")
        o.append(f'<rect x="{x0 + col[lo] * cell}" y="{y + 6}" '
                 f'width="{(hi - lo + 1) * cell - 2}" height="3" rx="1.5" '
                 f'fill="{INK3}" fill-opacity="0.35"/>')
    for i in range(len(runs) - 1):
        bx = x0 + (col[runs[i][1]] + 1) * cell + cell / 2
        label(bx, y + 10, "//", "m-n", INK3, "middle")
    label(0, y, "antigen residue", "m-n", INK3)
    y += 28

    for question, names in families:
        y += 20
        label(0, y, question.upper(), "m-q", INK3)
        y += 8
        for name in names:
            m = masks[name]
            orph = set(m.get("cysteines_orphaned") or [])
            partners = m.get("orphan_partners") or {}
            # readable name over the id, because the id is what the construct
            # picker, the shortlist and the files are all called: a reader who
            # meets "Main contact loop" here has to recognise epi_core there
            label(0, y + 13, MASK_LABELS.get(name, name), "m-h", INK)
            label(0, y + 27, name, "m-n", INK3)
            for pos, text in blocks[name]:
                for i, aa in enumerate(text):
                    p = None if pos is None else pos + i
                    if pos is None:                  # the spacer, drawn in the gap
                        bx = x0 + (col[runs[0][1]] + 1) * cell
                        o.append(f'<rect x="{bx}" y="{y + 11}" width="{cell - 2}" '
                                 f'height="9" rx="2" fill="{PLAIN}" fill-opacity="0.2"/>')
                        label(bx + (cell - 2) / 2, y + 19, text, "m-n", INK3, "middle")
                        break
                    if aa == "C":
                        fill = ORPH if p in orph else KEPT
                    elif core_start <= p <= core_end and core_seq[p - core_start] != aa:
                        fill = SUBST
                    else:
                        fill = PLAIN
                    solid = aa == "C" or fill == SUBST
                    bx = x0 + col[p] * cell
                    o.append(f'<rect x="{bx}" y="{y + 5}" width="{cell - 2}" height="21" '
                             f'rx="2" fill="{fill}" '
                             f'fill-opacity="{0.9 if solid else 0.16}"/>')
                    if p in hotspots:
                        o.append(f'<rect x="{bx}" y="{y}" width="{cell - 2}" '
                                 f'height="2.5" rx="1" fill="{HOT}"/>')
                    label(bx + (cell - 2) / 2, y + 20, aa, "m-t",
                          "#FFFFFF" if solid else INK2, "middle")
            # a closed bracket says the partner is inside the mask; an open
            # stub says it was left behind in the antigen
            cy = y + 29
            for lo_c, hi_c in (m.get("disulfides_kept") or []):
                if lo_c in col and hi_c in col:
                    xa = x0 + col[lo_c] * cell + (cell - 2) / 2
                    xb = x0 + col[hi_c] * cell + (cell - 2) / 2
                    o.append(f'<path d="M {xa} {cy - 4} L {xa} {cy} L {xb} {cy} '
                             f'L {xb} {cy - 4}" fill="none" stroke="{KEPT}" '
                             f'stroke-width="1.3"/>')
            for c in sorted(orph):
                if c not in col:
                    continue
                xa = x0 + col[c] * cell + (cell - 2) / 2
                out = -1 if c < core_start else 1
                o.append(f'<path d="M {xa} {cy - 4} L {xa} {cy} '
                         f'L {xa + out * (cell + 3)} {cy}" fill="none" '
                         f'stroke="{ORPH}" stroke-width="1.3" '
                         f'stroke-dasharray="3 2.4"/>')
            if orph:
                note = ("cysteine cut from C" + ", C".join(str(v) for v in partners.values())
                        if len(orph) == 1 else
                        f"{len(orph)} cysteines cut from C"
                        + ", C".join(str(v) for v in partners.values()))
            elif name == core:
                note = f"C{core_start}\u2013C{core_end} disulphide kept"
            else:
                note = MASK_NOTES.get(name, "")
            label(x0 + ncol * cell + 14, y + 20, note, "m-n", INK3)
            y += 38

    y += 14
    o.append(f'<path d="M 2 {y - 9} L 2 {y - 4} L 16 {y - 4} L 16 {y - 9}" '
             f'fill="none" stroke="{KEPT}" stroke-width="1.3"/>')
    label(24, y, "cysteines bonded to each other inside the mask", "m-b", INK3)
    y += 20
    o.append(f'<path d="M 2 {y - 9} L 2 {y - 4} L 16 {y - 4}" fill="none" '
             f'stroke="{ORPH}" stroke-width="1.3" stroke-dasharray="3 2.4"/>')
    label(24, y, "cysteine whose partner stays behind in the antigen", "m-b", INK3)
    y += 20
    o.append(f'<rect x="2" y="{y - 10}" width="12" height="12" rx="2" '
             f'fill="{SUBST}" fill-opacity="0.9"/>')
    label(24, y, "residue differing from the main contact loop", "m-b", INK3)
    y += 20
    o.append(f'<rect x="2" y="{y - 6}" width="12" height="2.5" rx="1" fill="{HOT}"/>')
    label(24, y, "one of the residues burying the most surface", "m-b", INK3)
    y += 6
    H = int(y - 2)
    o[0] = (f'<svg viewBox="0 0 {int(W)} {H}" width="100%" '
            'xmlns="http://www.w3.org/2000/svg" role="img" '
            'aria-label="Mask families, grouped by the question each answers">')
    o.append("</svg>")
    return "".join(o)


def linker_schematic(gen, lengths, substrate_len, reference) -> str:
    """How a protease-cleavable linker was designed, in four steps.

    Sequence specificity and linker geometry are separate decisions, and the
    diagram keeps them apart because the pipeline does. A generator picks the
    residues the protease has to recognise; the flanks around those residues
    set how far the mask can reach. Varying one does not confound the other.
    """
    W = 980
    CUT, FLANKC = "#2D5BD1", "#9AA6AE"
    GOOD, BAD, NEUT = "#1F6B4A", "#9C3B3B", "#6E7880"
    o = [None,                               # the svg tag, once H is known
         '<defs><marker id="lk" viewBox="0 0 10 10" refX="9" refY="5" '
         'markerWidth="7" markerHeight="7" orient="auto">'
         f'<path d="M 0 0 L 10 5 L 0 10 z" fill="{INK3}"/></marker>'
         '<marker id="lk2" viewBox="0 0 10 10" refX="5" refY="5" '
         'markerWidth="6" markerHeight="6" orient="auto">'
         f'<path d="M 10 0 L 0 5 L 10 10 z" fill="{INK3}"/></marker></defs>',
         '<style>'
         '.l-h{font:600 15px Inter,system-ui,sans-serif}'
         '.l-b{font:13px Inter,system-ui,sans-serif}'
         '.l-m{font:12.5px ui-monospace,monospace}'
         '.l-s{font:600 11px ui-monospace,monospace;letter-spacing:.11em}'
         '.l-n{font:600 12px ui-monospace,monospace}'
         '.l-xl{font:600 30px ui-monospace,monospace}'
         '</style>']

    def label(x, y, text, cls="l-b", fill=INK2, anchor="start"):
        o.append(f'<text x="{x}" y="{y}" class="{cls}" fill="{fill}" '
                 f'text-anchor="{anchor}">{text}</text>')

    def head(n, title, y):
        o.append(f'<rect x="0" y="{y}" width="24" height="24" rx="4" fill="{INK}"/>')
        label(12, y + 17, str(n), "l-n", "#FFFFFF", "middle")
        label(36, y + 17, title, "l-h", INK)
        return y + 24

    def panel(x, y, w, h, colour, opacity="0.05"):
        o.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="5" '
                 f'fill="{colour}" fill-opacity="{opacity}" stroke="{colour}" '
                 f'stroke-opacity="0.45" stroke-width="1.2"/>')

    # --- 1 · the profile the generator is asked for ----------------------
    y = head(1, "Define the desired cleavage profile", 0)
    cw, gap, cy, ch = (W - 32) / 3, 16, y + 14, 106
    groups = [(GOOD, gen["n_target"], "High cleavage", gen["target_names"],
               gen["target_z"]),
              (BAD, gen["n_off"], "Low cleavage", gen["off_names"],
               gen["off_z"]),
              (NEUT, gen["n_neutral"], "Unspecified", gen["neutral_names"], 0.0)]
    for i, (col, n, want, names, z) in enumerate(groups):
        x = i * (cw + gap)
        panel(x, cy, cw, ch, col)
        o.append(f'<rect x="{x}" y="{cy}" width="4" height="{ch}" rx="2" fill="{col}"/>')
        label(x + 16, cy + 27, f"{n} MMPs", "l-h", col)
        label(x + cw - 14, cy + 27,
              "z = 0" if not z else f"z = {z:+.0f}", "l-n", col, "end")
        label(x + 16, cy + 48, want, "l-b", INK2)
        ty = cy + 70
        for line in _wrap(_mmp_names(names), 32):
            label(x + 16, ty, line, "l-m", INK3)
            ty += 15
    ny = cy + ch + 22
    for line in _wrap(
            f"The {gen['n_neutral']} unspecified proteases are left at zero, which "
            "asks for no preference rather than for resistance, and they are "
            "excluded from the score in step 3.", 108):
        label(0, ny, line, "l-b", INK2)
        ny += 19

    # --- 2 · generation ---------------------------------------------------
    y = head(2, "CleaveNet \u2014 generate sequences", ny + 8)
    py, phh = y + 14, 64
    panel(0, py, W, phh, CUT, "0.06")
    n_str = f"{gen['n_generated']:,}"
    label(24, py + 44, n_str, "l-xl", CUT)
    nx = 24 + len(n_str) * 30 * 0.602 + 26
    label(nx, py + 32, f"unique {substrate_len}-residue substrates", "l-h", INK)
    label(nx, py + 51,
          "drawn in one conditional run against the profile above", "l-m", INK3)

    # --- 3 · scoring and selection ---------------------------------------
    y = head(3, "Predict cleavage and select candidates", py + phh + 26)
    ey, ehh = y + 14, 50
    panel(0, ey, 430, ehh, INK3, "0.04")
    o.append(f'<text x="24" y="{ey + 32}" class="l-m" fill="{INK}">'
             f'S  =  \u27e8z\u27e9<tspan dy="4" font-size="10">target</tspan>'
             f'<tspan dy="-4">  \u2212  \u27e8z\u27e9</tspan>'
             f'<tspan dy="4" font-size="10">off-target</tspan></text>')
    label(452, ey + 22,
          f"Retain the {gen['n_kept']} highest-scoring sequences", "l-h", INK)
    label(452, ey + 42,
          f"and add {gen['n_controls']} controls", "l-b", INK2)

    # --- 4 · assembly, and the length series cut from it -----------------
    y = head(4, "Assemble linkers, and vary their length", ey + ehh + 26)
    unit, by, bh = 16.0, y + 18, 38
    fl = (reference - substrate_len) // 2
    segs = [("Gly/Ser flank", FLANKC, fl * unit, 0.34),
            ("cleavage sequence", CUT, substrate_len * unit, 0.85),
            ("Gly/Ser flank", FLANKC, fl * unit, 0.34)]
    x = 0
    for name, col, w, op in segs:
        o.append(f'<rect x="{x}" y="{by}" width="{w}" height="{bh}" rx="3" '
                 f'fill="{col}" fill-opacity="{op}"/>')
        label(x + w / 2, by + 24, name, "l-m",
              "#FFFFFF" if col == CUT else INK2, "middle")
        x += w
    label(x + 18, by + 24, f"{reference} residues", "l-m", INK3)

    # the series itself: only the flanks change, so the substrate column is
    # the same width in every row and the comparison is reach, not sequence
    lo, hi = min(lengths), max(lengths)
    runit, rx0, ry = 9.4, 92.0, by + bh + 26
    longest = max(lengths)
    for i, n in enumerate(sorted(lengths)):
        flank = (n - substrate_len) // 2
        ry_i = ry + i * 30
        label(0, ry_i + 16, f"n = {n}", "l-n", INK if n != reference else CUT)
        fw, cw = flank * runit, substrate_len * runit
        o.append(f'<rect x="{rx0}" y="{ry_i}" width="{fw}" height="20" rx="2" '
                 f'fill="{FLANKC}" fill-opacity="0.30"/>')
        o.append(f'<rect x="{rx0 + fw}" y="{ry_i}" width="{cw}" height="20" rx="2" '
                 f'fill="{CUT}" fill-opacity="0.85"/>')
        o.append(f'<rect x="{rx0 + fw + cw}" y="{ry_i}" width="{fw}" height="20" '
                 f'rx="2" fill="{FLANKC}" fill-opacity="0.30"/>')
        label(rx0 + longest * runit + 18, ry_i + 15,
              f"{flank} Gly/Ser each side", "l-m", INK3)
    o.append(f'<line x1="{rx0}" y1="{ry - 4}" x2="{rx0}" '
             f'y2="{ry + (len(lengths) - 1) * 30 + 24}" stroke="{RULE}" '
             f'stroke-width="1"/>')

    ly = ry + (len(lengths) - 1) * 30 + 48
    for line in _wrap(
            f"Varying the flanks alone gives linkers of {lo} to {hi} residues. "
            f"The {substrate_len}-residue substrate is identical in every one, so a "
            "difference in masking across the series is reach and not cleavability.",
            108):
        label(0, ly, line, "l-b", INK2)
        ly += 19

    H = int(ly + 4)
    o[0] = (f'<svg viewBox="0 0 {W} {H}" width="100%" '
            'xmlns="http://www.w3.org/2000/svg" role="img" '
            'aria-label="Four-step design of a protease-cleavable linker">')
    o.append("</svg>")
    return "".join(o)
