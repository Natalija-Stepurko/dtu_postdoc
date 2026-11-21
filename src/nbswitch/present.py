"""The presentation layer: a talk that can be given, with the detail kept underneath.

Two levels. The surface carries large structures, one question per section and a
short answer. Everything technical - equations, model names, limitations,
numbers - sits in disclosure blocks that stay shut while someone is talking.
"""
from __future__ import annotations

#: One colour per molecular component, used identically in every figure,
#: every viewer and every legend on the page.
COMPONENT = {
    "vhh": "#337D79",      # the nanobody, fixed throughout
    "toxin": "#C99555",    # the antigen
    "mask": "#8C72AD",     # the designed module
    "linker": "#9AA7B1",   # the tether
    "struct": "#527FA4",   # measured or predicted structure, as data
    "ink": "#1B2838",
}

CSS = """
:root{
  --vhh:#337D79; --toxin:#C99555; --mask:#8C72AD; --linker:#9AA7B1; --struct:#527FA4;
  --ink:#1B2838; --ink-2:#44586B; --ink-3:#7B8B99;
  --bg:#FFFFFF; --panel:#FFFFFF; --rule:#DFE5EA; --rule-2:#EFF3F5; --wash:#F5F8F9;
  --sans:Inter,"Source Sans 3",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  --mono:ui-monospace,"SF Mono",SFMono-Regular,Menlo,Consolas,monospace;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--sans);
     font-size:19px;line-height:1.6;-webkit-font-smoothing:antialiased}
.page{max-width:1180px;margin:0 auto;padding:0 40px}
h1,h2,h3{margin:0;line-height:1.14;text-wrap:balance;letter-spacing:-.02em}
h1{font-size:clamp(2.3rem,5vw,3.1rem);font-weight:650}
h2{font-size:clamp(1.8rem,3.4vw,2.35rem);font-weight:650}
h3{font-size:1.25rem;font-weight:600;letter-spacing:-.01em}
p{margin:0;max-width:66ch;color:var(--ink-2);font-size:1.05rem}
p.big{font-size:1.2rem;color:var(--ink-2);max-width:62ch}
a{color:var(--struct)}
.mono{font-family:var(--mono)}
.dim{color:var(--ink-3)}

/* ---- chrome ---- */
.topbar{position:sticky;top:0;z-index:40;background:rgba(255,255,255,.93);
  backdrop-filter:blur(8px);border-bottom:1px solid var(--rule)}
.topbar .inner{max-width:1180px;margin:0 auto;padding:11px 40px;display:flex;
  align-items:center;gap:18px}
.brand{font-family:var(--mono);font-size:.68rem;letter-spacing:.16em;
  text-transform:uppercase;color:var(--ink-3);flex:1}
.modes{display:flex;gap:0;border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.modes button{font-family:var(--mono);font-size:.68rem;letter-spacing:.08em;
  text-transform:uppercase;padding:7px 14px;border:0;background:#fff;color:var(--ink-3);
  cursor:pointer}
.modes button[aria-pressed=true]{background:var(--ink);color:#fff}

/* ---- sections ---- */
section.slide{padding:84px 0 24px;border-top:1px solid var(--rule);scroll-margin-top:60px}
section.slide:first-of-type{border-top:0}
.snum{font-family:var(--mono);font-size:.68rem;letter-spacing:.16em;color:var(--ink-3);
  margin-bottom:16px;display:flex;gap:14px;align-items:baseline}
.snum b{color:var(--vhh);font-weight:650}
.lede{font-size:1.3rem;line-height:1.5;color:var(--ink-2);max-width:60ch;margin-top:14px}
.stack{display:flex;flex-direction:column;gap:14px}
.gap{height:26px}

/* ---- question / approach / output / limitation ---- */
.qaol{display:grid;grid-template-columns:repeat(auto-fit,minmax(215px,1fr));gap:0;
  border:1px solid var(--rule);border-radius:3px;overflow:hidden;margin-top:8px}
.qaol>div{padding:17px 19px;border-right:1px solid var(--rule-2)}
.qaol>div:last-child{border-right:0}
.qaol .k{font-family:var(--mono);font-size:.62rem;letter-spacing:.12em;
  text-transform:uppercase;color:var(--ink-3);margin-bottom:6px}
.qaol .v{font-size:.98rem;color:var(--ink);line-height:1.45}
.qaol .lim .k{color:var(--toxin)}

/* ---- molecular viewers ---- */
.mol{margin:22px 0}
.mol .frame{position:relative;width:100%;height:clamp(420px,56vh,620px);
  border:1px solid var(--rule);border-radius:3px;background:#fff;overflow:hidden}
.mol.pair{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.mol.pair .frame{height:clamp(340px,42vh,460px)}
.molbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:10px}
.molbar .ttl{font-family:var(--mono);font-size:.68rem;letter-spacing:.12em;
  text-transform:uppercase;color:var(--ink-3);flex:1}
.molbar button,.molbar select{font-family:var(--mono);font-size:.7rem;padding:5px 11px;
  border:1px solid var(--rule);border-radius:3px;background:#fff;color:var(--ink-2);
  cursor:pointer}
.molbar button:hover,.molbar select:hover{border-color:var(--ink-3);color:var(--ink)}
.molbar button[aria-pressed=true]{background:var(--ink);color:#fff;border-color:var(--ink)}
.legend{display:flex;gap:18px;flex-wrap:wrap;margin-top:11px;font-size:.84rem;
  color:var(--ink-2)}
.legend i{display:inline-block;width:13px;height:13px;border-radius:2px;
  margin-right:7px;vertical-align:-2px}
.cap{margin-top:9px;font-size:.88rem;color:var(--ink-3);max-width:76ch}
.badge{display:inline-block;font-family:var(--mono);font-size:.6rem;letter-spacing:.1em;
  text-transform:uppercase;padding:2px 8px;border-radius:999px;margin-right:9px;
  border:1px solid currentColor}
.badge.exp{color:var(--vhh)} .badge.pred{color:var(--struct)}
.badge.concept{color:var(--ink-3)}

/* ---- figures ---- */
figure.fig{margin:20px 0;border:1px solid var(--rule);border-radius:3px;padding:18px;
  background:#fff}
figure.fig img,figure.fig svg{width:100%;height:auto;display:block}
figure.fig figcaption{margin-top:11px;font-size:.86rem;color:var(--ink-3)}
.two{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:16px}

/* ---- numbers ---- */
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(168px,1fr));gap:14px;
  margin:20px 0}
.stat{border-left:2px solid var(--rule);padding:3px 0 3px 15px}
.stat.on{border-left-color:var(--vhh)}
.stat.off{border-left-color:var(--toxin)}
.stat .v{font-family:var(--mono);font-variant-numeric:tabular-nums;font-size:1.75rem;
  font-weight:500;letter-spacing:-.02em;line-height:1.1}
.stat .l{font-size:.85rem;color:var(--ink-2);margin-top:3px}
.stat .s{font-size:.78rem;color:var(--ink-3);margin-top:2px}

/* ---- technical disclosure ---- */
details.tech{border:1px solid var(--rule);border-radius:3px;margin:16px 0;background:#fff}
details.tech>summary{cursor:pointer;list-style:none;padding:13px 18px;
  font-family:var(--mono);font-size:.7rem;letter-spacing:.11em;text-transform:uppercase;
  color:var(--ink-3);display:flex;gap:11px;align-items:center}
details.tech>summary::-webkit-details-marker{display:none}
details.tech>summary::before{content:"+";font-size:.95rem;color:var(--struct);width:11px}
details.tech[open]>summary::before{content:"−"}
details.tech>summary:hover{color:var(--ink)}
.techbody{padding:0 18px 20px;border-top:1px solid var(--rule-2);padding-top:17px;
  display:flex;flex-direction:column;gap:15px}
.techbody p{font-size:.97rem}
body[data-mode="presentation"] details.tech{display:none}

/* ---- tables ---- */
.scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:3px;background:#fff}
table{border-collapse:collapse;width:100%;font-family:var(--mono);font-size:.8rem;
  font-variant-numeric:tabular-nums}
th,td{padding:9px 13px;text-align:right;white-space:nowrap;
  border-bottom:1px solid var(--rule-2)}
th:first-child,td:first-child{text-align:left;font-family:var(--sans);font-size:.92rem;
  white-space:normal}
thead th{font-size:.62rem;letter-spacing:.09em;text-transform:uppercase;
  color:var(--ink-3);font-weight:600;border-bottom:1px solid var(--rule)}
tbody tr:last-child td{border-bottom:0}
tr.good td{background:#F2F7F6} tr.bad td{background:#FBF5EE}

/* ---- equations ---- */
.eq{background:var(--wash);border:1px solid var(--rule);border-radius:3px;
  padding:15px 19px;font-family:var(--mono);font-size:.82rem;line-height:1.95;
  overflow-x:auto;color:var(--ink)}

/* ---- verdict banner ---- */
.note{border-left:2px solid var(--struct);padding:3px 0 3px 17px;margin:18px 0;
  font-size:1.02rem;color:var(--ink-2);max-width:72ch}
.note b{color:var(--ink)}
.note.warn{border-left-color:var(--toxin)}
.note.good{border-left-color:var(--vhh)}

footer{padding:40px 0 80px;color:var(--ink-3);font-size:.85rem;
  border-top:1px solid var(--rule);margin-top:60px}
ul.ref{list-style:none;padding:0;margin:0;font-size:.88rem}
ul.ref li{padding:9px 0;border-bottom:1px solid var(--rule-2);color:var(--ink-2)}
ul.ref li:last-child{border-bottom:0}
@media(max-width:860px){.page{padding:0 20px}.mol.pair{grid-template-columns:1fr}
  body{font-size:17px}}
"""
