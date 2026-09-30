/* ------------------------------------------------------------------ utils */
const $ = s => document.querySelector(s);
let CVD = false;
const el = (t, a = {}, ...kids) => { const n = document.createElement(t);
  for (const [k, v] of Object.entries(a)) { if (k === 'class') n.className = v;
    else if (k === 'html') n.innerHTML = v; else if (k.startsWith('on')) n[k] = v;
    else n.setAttribute(k, v); }
  kids.flat().forEach(c => n.append(c && c.nodeType ? c : document.createTextNode(c ?? '')));
  return n; };
const fmt = (x, d = 2) => (x === null || x === undefined || Number.isNaN(x)) ? '–'
  : (typeof x === 'number' ? (Number.isInteger(x) ? x : x.toFixed(d)) : x);
const nice = s => String(s).replace(/_/g, ' ');
/* A <select> reports the first option's value in a browser, but code that
   assumes it can crash on an empty value. Always fall back to a known-good
   member of the list. */
const pick = (sel, list) => (sel && list.includes(sel.value)) ? sel.value : list[0];
const dl = (name, text) => { const b = new Blob([text], {type: 'text/tab-separated-values'});
  const u = URL.createObjectURL(b); const a = el('a', {href: u, download: name}); a.click();
  URL.revokeObjectURL(u); };

const SP = D.species.map(s => s.name);
const SPI = Object.fromEntries(SP.map((s, i) => [s, i]));
const FAM = D.families.map(f => f.family);
const FMETA = Object.fromEntries(D.families.map(f => [f.family, f]));
const TRAITS = [...new Set(D.species.map(s => s.trait).filter(x => x))].sort();
const TCOL = ['#bc4749','#e9a13b','#2a9d8f','#1b4965','#7b2cbf','#6a994e','#9c6644'];
let traitColour = t => t == null ? '#c9ced4' : TCOL[TRAITS.indexOf(t) % TCOL.length];
const CATS = [...new Set(D.families.map(f => f.category))];
const CCOL = ['#1b4965','#2a9d8f','#8ecae6','#bc4749','#e9a13b','#7b2cbf','#6a994e','#9c6644'];
const catColour = c => CCOL[CATS.indexOf(c) % CCOL.length];

const tip = $('#tip');
function showTip(e, html) { tip.innerHTML = html; tip.style.opacity = 1;
  const pad = 14; let x = e.pageX + pad, y = e.pageY + pad;
  if (x + tip.offsetWidth > window.innerWidth + window.scrollX) x = e.pageX - tip.offsetWidth - pad;
  tip.style.left = x + 'px'; tip.style.top = y + 'px'; }
function hideTip() { tip.style.opacity = 0; }

/* matrix access: 'raw' | 'per10k' */
function matrix(kind) { return (kind === 'per10k' && D.per10k) ? D.per10k : D.counts; }
function col(kind, fam) { const m = matrix(kind); return m[fam] || SP.map(() => 0); }
function stats(v) { const n = v.length, mean = v.reduce((a, b) => a + b, 0) / n;
  const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / (n - 1 || 1));
  const s = [...v].sort((a, b) => a - b);
  const q = p => { const i = (s.length - 1) * p, lo = Math.floor(i), hi = Math.ceil(i);
    return s[lo] + (s[hi] - s[lo]) * (i - lo); };
  return {mean, sd, min: s[0], max: s[s.length - 1], q1: q(.25), med: q(.5), q3: q(.75)}; }
function zcol(kind, fam) { const v = col(kind, fam), st = stats(v);
  return st.sd === 0 ? v.map(() => 0) : v.map(x => (x - st.mean) / st.sd); }

/* families that carry no variance cannot be z-scored; say so rather than hide them */
const INVARIANT = {};
FAM.forEach(f => { const raw = D.counts[f]; if (!raw) return;
  const sum = raw.reduce((a, b) => a + b, 0);
  const allSame = raw.every(x => x === raw[0]);
  if (sum === 0) INVARIANT[f] = 'never detected in any species';
  else if (allSame) INVARIANT[f] = 'constant at ' + raw[0] + ' in every species'; });
const VARFAM = FAM.filter(f => !INVARIANT[f]);

/* diverging + sequential colour ramps, no dependencies */
function rdbu(t) { t = Math.max(0, Math.min(1, t));
  const a = [[5,48,97],[33,102,172],[67,147,195],[146,197,222],[209,229,240],
             [247,247,247],[253,219,199],[244,165,130],[214,96,77],[178,24,43],[103,0,31]];
  const i = Math.min(a.length - 2, Math.floor(t * (a.length - 1))), f = t * (a.length - 1) - i;
  const c = a[i].map((v, k) => Math.round(v + (a[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`; }
function viridis(t) { t = Math.max(0, Math.min(1, t));
  const a = [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
             [31,158,137],[53,183,121],[110,206,88],[181,222,43],[253,231,37]];
  const i = Math.min(a.length - 2, Math.floor(t * (a.length - 1))), f = t * (a.length - 1) - i;
  const c = a[i].map((v, k) => Math.round(v + (a[i + 1][k] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`; }

/* sortable table with search and TSV export */
function makeTable(host, cols, rows, opts = {}) {
  let sortI = opts.sortBy ?? 0, dir = opts.desc ? -1 : 1, filter = '';
  const search = el('input', {type: 'search', placeholder: 'filter…',
    oninput: e => { filter = e.target.value.toLowerCase(); draw(); }});
  const dlBtn = el('button', {class: 'act', onclick: () => dl((opts.name || 'table') + '.tsv',
    [cols.join('\t'), ...view().map(r => r.join('\t'))].join('\n'))}, 'download TSV');
  const count = el('span', {class: 'hint'});
  const wrap = el('div', {class: 'scroll'});
  host.append(el('div', {class: 'row'}, el('div', {class: 'ctl'},
    el('label', {}, 'search'), search), dlBtn, count), wrap);
  const view = () => { let r = rows;
    if (filter) r = r.filter(x => x.some(v => String(v).toLowerCase().includes(filter)));
    return [...r].sort((a, b) => { const x = a[sortI], y = b[sortI];
      if (typeof x === 'number' && typeof y === 'number') return (x - y) * dir;
      return String(x).localeCompare(String(y)) * dir; }); };
  function draw() { const r = view(); count.textContent = r.length + ' rows';
    const cap = opts.cap || 800;
    const t = el('table', {},
      el('thead', {}, el('tr', {}, cols.map((c, i) => el('th', {
        onclick: () => { if (sortI === i) dir = -dir; else { sortI = i; dir = 1; } draw(); }},
        c + (sortI === i ? (dir > 0 ? ' ▲' : ' ▼') : ''))))),
      el('tbody', {}, r.slice(0, cap).map(row => el('tr', {}, row.map((v, i) =>
        el('td', {class: opts.mono && opts.mono.includes(i) ? 'mono' : ''}, fmt(v)))))));
    wrap.innerHTML = ''; wrap.append(t);
    if (r.length > cap) wrap.append(el('div', {class: 'note'},
      `showing the first ${cap} of ${r.length} rows. Narrow with the filter box, or download the full TSV.`)); }
  draw();
}

/* ------------------------------------------------------------- OVERVIEW */
function tabOverview(host) {
  const nInv = Object.keys(INVARIANT).length;
  const flagged = D.species.filter(s => s.flag && s.flag !== 'ok').length;
  const totCalls = FAM.reduce((a, f) => a + (D.counts[f] || []).reduce((x, y) => x + y, 0), 0);
  host.append(el('div', {class: 'grid g3'},
    kpi(SP.length, 'species'), kpi(FAM.length, 'gene families'),
    kpi(totCalls.toLocaleString(), 'defensome genes called'),
    kpi(VARFAM.length, 'families with variance'),
    kpi(nInv, 'families invariant or absent'),
    kpi(flagged, 'species flagged by QC')));

  if (nInv) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, 'Families that cannot be z-scored'),
      el('p', {class: 'hint'}, 'A family with no variance has no z-score, so it never appears in a heatmap. Both cases below are findings, not gaps, so they are listed rather than dropped.'));
    p.append(el('div', {class: 'scroll'}, el('table', {},
      el('thead', {}, el('tr', {}, ['family', 'Pfam', 'why', 'what to do'].map(c => el('th', {}, c)))),
      el('tbody', {}, Object.entries(INVARIANT).map(([f, why]) => el('tr', {},
        el('td', {}, f), el('td', {class: 'mono'}, FMETA[f] ? FMETA[f].pfam : ''),
        el('td', {}, why),
        el('td', {}, why.startsWith('never')
          ? (D.qc_zero && (D.qc_zero.find(z => z.family === f) || {}).verdict
              ? 'verdict: ' + D.qc_zero.find(z => z.family === f).verdict.replace(/_/g, ' ').toLowerCase() + ' (see Map tab)'
              : 'run qc: it writes qc_zero_families.tsv with a verdict for every empty family')
          : 'single-copy control: the pipeline is calibrated')))))));
    host.append(p);
  }

  if (D.qc_families) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, 'Are the control families behaving?'),
      el('p', {class: 'hint'}, 'Observed median across species against the literature expectation for insect genomes. A high ratio means the map rule is catching more than the family name implies, which is a map problem rather than a genome problem.'));
    makeTable(p, ['family', 'expected', 'observed median', 'min', 'max', 'ratio', 'verdict'],
      D.qc_families.map(r => [r.family, r.expected, r.observed_median, r.observed_min,
        r.observed_max, r.median_vs_expected, r.verdict]), {sortBy: 5, desc: true, name: 'qc_families'});
    host.append(p);
  }

  const p2 = el('div', {class: 'panel'});
  p2.append(el('h2', {}, 'Species'), el('p', {class: 'hint'},
    'max_ratio compares each species against the observed median of the low-copy control families. Anything above about 3 has duplicated gene models.'));
  const cols = ['species', 'proteins', D.trait_name || 'trait', 'CORE total (raw)',
    'CORE per 10k', 'QC max ratio', 'worst family', 'flag'];
  const coreV = VARFAM.filter(f => FMETA[f] && FMETA[f].tier === 'CORE');
  const rows = SP.map((s, i) => [s, D.species[i].n_proteins, D.species[i].trait ?? '–',
    coreV.reduce((a, f) => a + D.counts[f][i], 0),
    D.per10k ? +coreV.reduce((a, f) => a + D.per10k[f][i], 0).toFixed(1) : null,
    D.species[i].max_ratio ?? null, D.species[i].worst_family ?? '–', D.species[i].flag ?? '–']);
  makeTable(p2, cols, rows, {sortBy: 5, desc: true, name: 'species_overview'});
  host.append(p2);
}
const kpi = (v, l) => el('div', {class: 'kpi'}, el('div', {class: 'v'}, String(v)),
  el('div', {class: 'l'}, l));

/* -------------------------------------------------------------- HEATMAP */
function tabHeatmap(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', scale: 'z', tiers: new Set(['CORE']),
    cats: new Set(CATS), sort: 'trait', sortFam: null};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Species by family'), el('p', {class: 'hint'},
    'Click a column header to sort species by that family. Click a cell for the family detail. Raw counts track proteome size; per-10k asks what share of the gene budget goes to defence.'));
  const ctrls = el('div', {class: 'row'});
  const mk = (label, node) => el('div', {class: 'ctl'}, el('label', {}, label), node);
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const scaleSel = el('select', {onchange: e => { S.scale = e.target.value; draw(); }},
    el('option', {value: 'z'}, 'z-score per family'),
    el('option', {value: 'abs'}, 'absolute value'),
    el('option', {value: 'log'}, 'log10(x+1)'));
  const sortSel = el('select', {onchange: e => { S.sort = e.target.value; S.sortFam = null; draw(); }},
    el('option', {value: 'trait'}, D.trait_name || 'trait'),
    el('option', {value: 'name'}, 'species name'),
    el('option', {value: 'total'}, 'total defensome'),
    el('option', {value: 'size'}, 'proteome size'));
  ctrls.append(mk('values', kindSel), mk('scale', scaleSel), mk('sort species by', sortSel));
  const tierChips = el('div', {class: 'chips'});
  ['CORE', 'BROAD'].forEach(t => tierChips.append(el('span', {
    class: 'chip' + (S.tiers.has(t) ? ' on' : ''),
    onclick: e => { S.tiers.has(t) ? S.tiers.delete(t) : S.tiers.add(t);
      e.target.classList.toggle('on'); draw(); }}, t)));
  const catChips = el('div', {class: 'chips'});
  CATS.forEach(c => catChips.append(el('span', {class: 'chip on',
    onclick: e => { S.cats.has(c) ? S.cats.delete(c) : S.cats.add(c);
      e.target.classList.toggle('on'); draw(); }}, c)));
  ctrls.append(mk('tier', tierChips), mk('category', catChips));
  p.append(ctrls);
  const cv = el('canvas'); const holder = el('div', {style: 'overflow:auto'});
  holder.append(cv); p.append(figure(cv, 'defensome_heatmap'));
  const cbar = el('div', {style: 'margin-top:10px'}); p.append(cbar);
  const leg = el('div', {class: 'legend'}); p.append(leg);
  host.append(p);

  function draw() {
    const fams = FAM.filter(f => S.tiers.has(FMETA[f].tier) && S.cats.has(FMETA[f].category)
      && !(S.scale === 'z' && INVARIANT[f]));
    let order = SP.map((s, i) => i);
    const coreV = fams.length ? fams : VARFAM;
    const tot = i => coreV.reduce((a, f) => a + col(S.kind, f)[i], 0);
    if (S.sortFam) { const v = col(S.kind, S.sortFam); order.sort((a, b) => v[b] - v[a]); }
    else if (S.sort === 'name') order.sort((a, b) => SP[a].localeCompare(SP[b]));
    else if (S.sort === 'total') order.sort((a, b) => tot(b) - tot(a));
    else if (S.sort === 'size') order.sort((a, b) => D.species[b].n_proteins - D.species[a].n_proteins);
    else order.sort((a, b) => String(D.species[a].trait).localeCompare(String(D.species[b].trait))
      || tot(b) - tot(a));

    const CW = 26, CH = 17, L = 210, T = 130, R = 30;
    cv.width = L + fams.length * CW + R; cv.height = T + order.length * CH + 30;
    const g = cv.getContext('2d'); g.clearRect(0, 0, cv.width, cv.height);
    const vals = {}; fams.forEach(f => {
      vals[f] = S.scale === 'z' ? zcol(S.kind, f)
        : S.scale === 'log' ? col(S.kind, f).map(x => Math.log10(x + 1))
        : col(S.kind, f); });
    let lo, hi;
    if (S.scale === 'z') { lo = -2.5; hi = 2.5; }
    else { const all = fams.flatMap(f => vals[f]); lo = 0; hi = Math.max(...all) || 1; }
    g.font = '11px sans-serif';
    fams.forEach((f, j) => { g.save(); g.translate(L + j * CW + CW / 2, T - 8);
      g.rotate(-Math.PI / 2.4); g.textAlign = 'left';
      g.fillStyle = INVARIANT[f] ? '#9aa1a8' : '#1c1f23'; g.fillText(f, 0, 0); g.restore(); });
    order.forEach((si, r) => {
      g.fillStyle = traitColour(D.species[si].trait);
      g.fillRect(L - 12, T + r * CH + 2, 7, CH - 4);
      g.fillStyle = '#1c1f23'; g.textAlign = 'right'; g.font = 'italic 11px sans-serif';
      g.fillText(nice(SP[si]), L - 18, T + r * CH + CH - 5);
      fams.forEach((f, j) => { const v = vals[f][si];
        g.fillStyle = S.scale === 'z' ? rdbu((v - lo) / (hi - lo)) : viridis((v - lo) / (hi - lo));
        g.fillRect(L + j * CW, T + r * CH, CW - 1, CH - 1); }); });
    cv.onmousemove = e => { const b = cv.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW), r = Math.floor((e.clientY - b.top - T) / CH);
      if (j < 0 || r < 0 || j >= fams.length || r >= order.length) return hideTip();
      const f = fams[j], si = order[r];
      showTip(e, `<b>${nice(SP[si])}</b><br>${f} (${FMETA[f].category}, ${FMETA[f].tier})<br>
        raw ${D.counts[f][si]}${D.per10k ? ' · per10k ' + D.per10k[f][si].toFixed(2) : ''}<br>
        z ${zcol(S.kind, f)[si].toFixed(2)}${INVARIANT[f] ? '<br><i>' + INVARIANT[f] + '</i>' : ''}
        ${D.species[si].trait ? '<br>' + D.trait_name + ': ' + D.species[si].trait : ''}`); };
    cv.onmouseleave = hideTip;
    cv.onclick = e => { const b = cv.getBoundingClientRect();
      const j = Math.floor((e.clientX - b.left - L) / CW), r = Math.floor((e.clientY - b.top - T) / CH);
      if (j >= 0 && j < fams.length) { if (e.clientY - b.top < T) { S.sortFam = fams[j]; draw(); }
        else if (r >= 0 && r < order.length) { go('family'); famSel.value = fams[j];
          famSel.dispatchEvent(new Event('change')); } } };
    cbar.innerHTML = '';
    cbar.append(colourbar(S.scale === 'z' ? rdbu : viridis, lo, hi,
      S.scale === 'z' ? `z-score, ${S.kind === 'per10k' ? `copies per 10k ${D.size_unit || 'proteins'}` : 'raw copies'}`
                      : (S.scale === 'log' ? 'log10(x+1)' : (S.kind === 'per10k' ? `copies per 10k ${D.size_unit || 'proteins'}` : 'raw copies')),
      280, 12));
    leg.innerHTML = '';
    TRAITS.forEach(t => leg.append(el('span', {}, el('i', {class: 'sw',
      style: 'background:' + traitColour(t)}), t)));
    const hidden = FAM.filter(f => S.scale === 'z' && INVARIANT[f]
      && S.tiers.has(FMETA[f].tier) && S.cats.has(FMETA[f].category));
    if (hidden.length) leg.append(el('span', {style: 'color:#6b7280'},
      '· hidden (no variance, cannot be z-scored): ' + hidden.join(', ')));
  }
  draw();
}

/* --------------------------------------------------------------- FAMILY */
let famSel;
function tabFamily(host) {
  const p = el('div', {class: 'panel'});
  famSel = el('select', {}, ...FAM.map(f => el('option', {value: f}, f)));
  const kindSel = el('select', {}, ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  p.append(el('h2', {}, 'Family detail'), el('p', {class: 'hint'},
    'Everything the run knows about one gene family: how it varies, how it splits by trait, and whether its domain architectures are intact.'),
    el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'family'), famSel),
      el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel)));
  const body = el('div'); p.append(body); host.append(p);
  const draw = () => { body.innerHTML = ''; renderFamily(body, pick(famSel, FAM), pick(kindSel, ['per10k','raw'])); };
  famSel.onchange = draw; kindSel.onchange = draw; draw();
}

function renderFamily(host, fam, kind) {
  const m = FMETA[fam];
  if (!m) { host.append(el('div', {class: 'note'}, 'Unknown family: ' + fam)); return; }
  const v = col(kind, fam), st = stats(v);
  host.append(el('div', {class: 'note'},
    `${m.pfam} · rule ${m.rule} · min_cov ${m.min_cov} · min_len ${m.min_len} · tier ${m.tier} · ${m.category}`
    + (m.notes ? ' — ' + m.notes : '')));
  if (INVARIANT[fam]) host.append(el('div', {class: 'note',
    style: 'border-left-color:var(--warm)'}, 'This family is ' + INVARIANT[fam] +
    ', so it has no z-score and is excluded from the z-scored panels.'));
  host.append(el('div', {class: 'grid g3'},
    kpi(fmt(st.med), 'median'), kpi(fmt(st.min) + ' – ' + fmt(st.max), 'range'),
    kpi(st.mean ? (st.sd / st.mean).toFixed(2) : '–', 'coefficient of variation')));

  const g2 = el('div', {class: 'grid g2'});
  /* distribution by trait */
  if (TRAITS.length) {
    const c = el('div', {class: 'panel'});
    c.append(el('h2', {}, 'By ' + D.trait_name), el('p', {class: 'hint'},
      'Species are not independent, so treat any separation here as a screen and not a result.'));
    const groups = TRAITS.map(t => ({t, v: SP.map((s, i) => [D.species[i].trait, v[i]])
      .filter(x => x[0] === t).map(x => x[1])}));
    c.append(figure(boxPlot(groups, 520, 300), fam + '_by_trait'));
    const tbl = el('table', {}, el('thead', {}, el('tr', {},
      [D.trait_name, 'n', 'mean', 'median', 'min', 'max'].map(h => el('th', {}, h)))),
      el('tbody', {}, groups.map(gr => { const s = gr.v.length ? stats(gr.v) : null;
        return el('tr', {}, el('td', {}, gr.t), el('td', {}, gr.v.length),
          el('td', {}, s ? fmt(s.mean) : '–'), el('td', {}, s ? fmt(s.med) : '–'),
          el('td', {}, s ? fmt(s.min) : '–'), el('td', {}, s ? fmt(s.max) : '–')); })));
    c.append(tbl); g2.append(c);
  }
  /* domain completeness */
  if (D.completeness && D.completeness[fam]) {
    const c = el('div', {class: 'panel'});
    const cc = D.completeness[fam], tot = Object.values(cc).reduce((a, b) => a + b, 0);
    c.append(el('h2', {}, 'Domain completeness'), el('p', {class: 'hint'},
      'PARTIAL means the protein was called but a required domain is truncated below the family threshold. FRAGMENT means it is far shorter than the family norm.'));
    const bar = el('div', {style: 'display:flex;height:26px;border-radius:4px;overflow:hidden;margin:8px 0'});
    const cols = {COMPLETE: '#2a9d8f', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#adb5bd'};
    Object.entries(cc).forEach(([k, n]) => { if (!n) return;
      bar.append(el('div', {style: `width:${n / tot * 100}%;background:${cols[k]}`,
        title: `${k}: ${n}`})); });
    c.append(bar, el('div', {class: 'legend'}, Object.entries(cc).map(([k, n]) =>
      el('span', {}, el('i', {class: 'sw', style: 'background:' + cols[k]}),
        `${k} ${n} (${(n / tot * 100).toFixed(1)}%)`))));
    if (D.architectures) { const ar = D.architectures.filter(r => r[0] === fam).slice(0, 12);
      if (ar.length) { c.append(el('h2', {style: 'margin-top:14px'}, 'Domain architectures'));
        c.append(el('table', {}, el('thead', {}, el('tr', {},
          ['architecture', 'proteins'].map(h => el('th', {}, h)))),
          el('tbody', {}, ar.map(r => el('tr', {}, el('td', {class: 'mono'}, r[1]),
            el('td', {}, r[2])))))); } }
    g2.append(c);
  }
  host.append(g2);

  /* per-species bar */
  const c3 = el('div', {class: 'panel'});
  c3.append(el('h2', {}, 'Per species'));
  const order = SP.map((s, i) => i).sort((a, b) => v[b] - v[a]);
  c3.append(figure(barChart(order.map(i => ({label: nice(SP[i]), value: v[i],
    colour: traitColour(D.species[i].trait),
    extra: `raw ${D.counts[fam][i]} · ${D.species[i].n_proteins.toLocaleString()} proteins`})),
    880, Math.max(220, SP.length * 17)), fam + '_per_species'));
  host.append(c3);

  /* completeness per species for this family */
  if (D.comp_sf) {
    const rows = D.comp_sf.filter(r => r[1] === fam)
      .map(r => [r[0], r[2], r[3], r[4], r[5], +(r[2] / (r[2] + r[3] + r[4] + r[5]) * 100).toFixed(1)]);
    if (rows.length) { const c4 = el('div', {class: 'panel'});
      c4.append(el('h2', {}, 'Completeness per species'));
      makeTable(c4, ['species', 'complete', 'partial', 'fragment', 'missing', '% complete'],
        rows, {sortBy: 5, name: fam + '_completeness'});
      host.append(c4); }
  }
}

/* ------------------------------------------------------------- SPECIES */
function tabSpecies(host) {
  const p = el('div', {class: 'panel'});
  const sel = el('select', {}, ...SP.map(s => el('option', {value: s}, nice(s))));
  const kindSel = el('select', {}, ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  p.append(el('h2', {}, 'Species detail'), el('p', {class: 'hint'},
    'How one genome compares with the rest, family by family. The z-score is the number of standard deviations from the cross-species mean.'),
    el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'species'), sel),
      el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel)));
  const body = el('div'); p.append(body); host.append(p);
  const draw = () => { body.innerHTML = ''; renderSpecies(body, pick(sel, SP), pick(kindSel, ['per10k','raw'])); };
  sel.onchange = draw; kindSel.onchange = draw; draw();
}

function renderSpecies(host, sp, kind) {
  const i = SPI[sp], s = D.species[i];
  if (s === undefined) { host.append(el('div', {class: 'note'}, 'Unknown species: ' + sp)); return; }
  host.append(el('div', {class: 'grid g3'},
    kpi(s.n_proteins.toLocaleString(), `${D.size_unit || 'proteins'} in proteome`),
    kpi(s.trait ?? '–', D.trait_name || 'trait'),
    kpi(fmt(s.max_ratio), 'QC max ratio' + (s.flag && s.flag !== 'ok' ? ' — ' + s.flag : ''))));
  if (s.flag && s.flag !== 'ok') host.append(el('div', {class: 'note',
    style: 'border-left-color:var(--warm)'},
    `Flagged ${s.flag}. Worst control family: ${s.worst_family}. Counts from this genome are likely inflated by duplicated gene models; use max_ratio as a covariate or exclude it.`));
  const rows = FAM.map(f => { const v = col(kind, f), z = zcol(kind, f);
    return [f, FMETA[f].category, FMETA[f].tier, D.counts[f][i],
      D.per10k ? +D.per10k[f][i].toFixed(2) : null, +z[i].toFixed(2),
      +stats(v).med.toFixed(2), INVARIANT[f] ? INVARIANT[f] : '']; });
  const c = el('div', {class: 'panel'});
  c.append(el('h2', {}, 'Family profile'), el('p', {class: 'hint'},
    'Sorted by z-score: the top rows are the families where this genome is most unusual.'));
  makeTable(c, ['family', 'category', 'tier', 'raw', 'per 10k', 'z', 'median across species', 'note'],
    rows, {sortBy: 5, desc: true, name: sp + '_profile'});
  host.append(c);
  const c2 = el('div', {class: 'panel'});
  c2.append(el('h2', {}, 'Deviation from the cross-species mean'));
  const zr = VARFAM.map(f => ({label: f, value: zcol(kind, f)[i],
    colour: zcol(kind, f)[i] >= 0 ? '#bc4749' : '#1b4965', extra: 'raw ' + D.counts[f][i]}))
    .sort((a, b) => b.value - a.value);
  c2.append(figure(barChart(zr, 880, Math.max(220, zr.length * 16), true), sp + '_deviation'));
  host.append(c2);
}

/* --------------------------------------------------------- SVG PRIMITIVES */
const svgEl = (t, a = {}) => { const n = document.createElementNS('http://www.w3.org/2000/svg', t);
  for (const [k, v] of Object.entries(a)) n.setAttribute(k, v); return n; };

function boxPlot(groups, W, H) {
  const all = groups.flatMap(g => g.v); if (!all.length) return el('div', {}, 'no data');
  const lo = Math.min(...all), hi = Math.max(...all), pad = (hi - lo) * .1 || 1;
  const y0 = lo - pad, y1 = hi + pad;
  const ML = 56, MB = 44, MT = 12, MR = 12;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = v => MT + (H - MT - MB) * (1 - (v - y0) / (y1 - y0));
  const bw = (W - ML - MR) / groups.length;
  for (let k = 0; k <= 4; k++) { const v = y0 + (y1 - y0) * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(v), y2: Y(v), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(v) + 4, 'text-anchor': 'end',
      'font-size': 10, fill: '#6b7280'}); t.textContent = v.toFixed(1); sv.append(t); }
  groups.forEach((g, j) => { const cx = ML + bw * (j + .5), col = traitColour(g.t);
    if (g.v.length) { const s = stats(g.v);
      sv.append(svgEl('line', {x1: cx, x2: cx, y1: Y(s.min), y2: Y(s.max), stroke: '#556'}));
      sv.append(svgEl('rect', {x: cx - bw * .22, y: Y(s.q3), width: bw * .44,
        height: Math.max(1, Y(s.q1) - Y(s.q3)), fill: col, 'fill-opacity': .3,
        stroke: col}));
      sv.append(svgEl('line', {x1: cx - bw * .22, x2: cx + bw * .22, y1: Y(s.med),
        y2: Y(s.med), stroke: col, 'stroke-width': 2}));
      g.v.forEach((v, q) => sv.append(svgEl('circle', {
        cx: cx + ((q * 2654435761 % 100) / 100 - .5) * bw * .30, cy: Y(v), r: 2.6,
        fill: col, 'fill-opacity': .8, stroke: '#fff', 'stroke-width': .5}))); }
    const t = svgEl('text', {x: cx, y: H - MB + 16, 'text-anchor': 'middle', 'font-size': 11});
    t.textContent = g.t; sv.append(t);
    const t2 = svgEl('text', {x: cx, y: H - MB + 30, 'text-anchor': 'middle',
      'font-size': 10, fill: '#6b7280'}); t2.textContent = 'n=' + g.v.length; sv.append(t2); });
  return sv;
}

function barChart(items, W, H, signed = false) {
  const ML = 200, MR = 60, MT = 8;
  const rows = items.length, rh = Math.max(9, (H - MT) / rows);
  const hh = MT + rows * rh;
  const mx = Math.max(...items.map(d => Math.abs(d.value))) || 1;
  const zero = signed ? ML + (W - ML - MR) / 2 : ML;
  const span = signed ? (W - ML - MR) / 2 : (W - ML - MR);
  const sv = svgEl('svg', {width: W, height: hh, viewBox: `0 0 ${W} ${hh}`});
  if (signed) sv.append(svgEl('line', {x1: zero, x2: zero, y1: MT, y2: hh, stroke: '#ccc'}));
  items.forEach((d, i) => { const y = MT + i * rh, w = Math.abs(d.value) / mx * span;
    const x = signed ? (d.value >= 0 ? zero : zero - w) : ML;
    const r = svgEl('rect', {x, y: y + 1, width: Math.max(1, w), height: rh - 2,
      fill: d.colour || '#1b4965'});
    r.addEventListener('mousemove', e => showTip(e,
      `<b>${d.label}</b><br>${fmt(d.value, 3)}${d.extra ? '<br>' + d.extra : ''}`));
    r.addEventListener('mouseleave', hideTip); sv.append(r);
    if (rh >= 9) { const t = svgEl('text', {x: ML - 8, y: y + rh / 2 + 3.5,
      'text-anchor': 'end', 'font-size': Math.min(11, rh - 1), fill: '#1c1f23'});
      t.textContent = d.label; sv.append(t); }
    const t2 = svgEl('text', {x: (signed && d.value < 0 ? x - 5 : x + w + 5),
      y: y + rh / 2 + 3.5, 'text-anchor': signed && d.value < 0 ? 'end' : 'start',
      'font-size': Math.min(10, rh - 1), fill: '#6b7280'});
    t2.textContent = fmt(d.value, 2); sv.append(t2); });
  return sv;
}

/* ------------------------------------------------------------ ORDINATION */
function jacobiEig(A, n) {   /* symmetric eigendecomposition, enough for 40x40 */
  const V = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => i === j ? 1 : 0));
  const a = A.map(r => r.slice());
  for (let sweep = 0; sweep < 60; sweep++) {
    let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] ** 2;
    if (off < 1e-12) break;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (Math.abs(a[i][j]) < 1e-14) continue;
      const th = (a[j][j] - a[i][i]) / (2 * a[i][j]);
      const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < n; k++) { const aik = a[i][k], ajk = a[j][k];
        a[i][k] = c * aik - s * ajk; a[j][k] = s * aik + c * ajk; }
      for (let k = 0; k < n; k++) { const aki = a[k][i], akj = a[k][j];
        a[k][i] = c * aki - s * akj; a[k][j] = s * aki + c * akj; }
      for (let k = 0; k < n; k++) { const vki = V[k][i], vkj = V[k][j];
        V[k][i] = c * vki - s * vkj; V[k][j] = s * vki + c * vkj; } }
  }
  const idx = Array.from({length: n}, (_, i) => i).sort((x, y) => a[y][y] - a[x][x]);
  return {values: idx.map(i => a[i][i]), vectors: idx.map(i => V.map(r => r[i]))};
}

function tabOrdination(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', fams: new Set(
    VARFAM.filter(f => FMETA[f].tier === 'CORE'))};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Ordination'), el('p', {class: 'hint'},
    'Each family is z-scored first so a 150-copy family does not outweigh a 5-copy one, then the covariance matrix is decomposed. The loadings panel says which families drive the separation. If PC1 explains only 30 to 40 percent, do not read this as a single defence gradient.'));
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const chips = el('div', {class: 'chips'});
  VARFAM.forEach(f => chips.append(el('span', {class: 'chip' + (S.fams.has(f) ? ' on' : ''),
    onclick: e => { S.fams.has(f) ? S.fams.delete(f) : S.fams.add(f);
      e.target.classList.toggle('on'); draw(); }}, f)));
  p.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(VARFAM.filter(f =>
      FMETA[f].tier === 'CORE')); refreshChips(); draw(); }}, 'CORE only'),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(VARFAM); refreshChips(); draw(); }}, 'all'),
    el('button', {class: 'act', onclick: () => { S.fams = new Set(); refreshChips(); draw(); }}, 'none')));
  p.append(el('div', {class: 'ctl'}, el('label', {}, 'families included'), chips));
  const out = el('div'); p.append(out); host.append(p);
  const refreshChips = () => [...chips.children].forEach(c =>
    c.classList.toggle('on', S.fams.has(c.textContent)));
  function draw() { out.innerHTML = '';
    const fams = VARFAM.filter(f => S.fams.has(f));
    if (fams.length < 3) { out.append(el('div', {class: 'note'}, 'Select at least three families.')); return; }
    const Z = SP.map((_, i) => fams.map(f => zcol(S.kind, f)[i]));
    const n = fams.length, N = SP.length;
    const C = Array.from({length: n}, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      let s = 0; for (let k = 0; k < N; k++) s += Z[k][i] * Z[k][j];
      C[i][j] = C[j][i] = s / (N - 1); }
    const {values, vectors} = jacobiEig(C, n);
    const tot = values.reduce((a, b) => a + Math.max(0, b), 0);
    const scores = Z.map(r => [0, 1].map(k => r.reduce((a, v, q) => a + v * vectors[k][q], 0)));
    out.append(el('div', {class: 'grid g2'},
      wrapPanel('Species scores', figure(scatter(scores.map((s, i) => ({x: s[0], y: s[1],
        label: nice(SP[i]), colour: traitColour(D.species[i].trait),
        extra: (D.species[i].trait ?? '') + ' · ' + D.species[i].n_proteins.toLocaleString() + ' proteins'})),
        620, 520, `PC1 (${(values[0] / tot * 100).toFixed(1)}%)`,
        `PC2 (${(values[1] / tot * 100).toFixed(1)}%)`, true), 'pca_scores')),
      wrapPanel('Family loadings', figure(scatter(fams.map((f, q) => ({x: vectors[0][q],
        y: vectors[1][q], label: f, colour: catColour(FMETA[f].category),
        extra: FMETA[f].category})), 620, 520, 'PC1 loading', 'PC2 loading', true),
        'pca_loadings'))));
    const varRows = values.slice(0, 8).map((v, k) => [`PC${k + 1}`,
      +(v / tot * 100).toFixed(2), +(values.slice(0, k + 1).reduce((a, b) => a + b, 0) / tot * 100).toFixed(2)]);
    const vp = el('div', {class: 'panel'});
    vp.append(el('h2', {}, 'Variance explained'));
    makeTable(vp, ['component', '% variance', 'cumulative %'], varRows, {name: 'pca_variance'});
    out.append(vp);
  }
  draw();
}
const wrapPanel = (title, node) => el('div', {class: 'panel'}, el('h2', {}, title), node);

function scatter(pts, W, H, xl, yl, labels) {
  const ML = 56, MB = 44, MT = 12, MR = 16;
  const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
  const pad = v => { const lo = Math.min(...v), hi = Math.max(...v), d = (hi - lo) * .12 || 1;
    return [lo - d, hi + d]; };
  const [x0, x1] = pad(xs), [y0, y1] = pad(ys);
  const X = v => ML + (W - ML - MR) * (v - x0) / (x1 - x0);
  const Y = v => MT + (H - MT - MB) * (1 - (v - y0) / (y1 - y0));
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  sv.append(svgEl('line', {x1: X(0), x2: X(0), y1: MT, y2: H - MB, stroke: '#e3e6ea'}));
  sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(0), y2: Y(0), stroke: '#e3e6ea'}));
  pts.forEach(p => { const c = svgEl('circle', {cx: X(p.x), cy: Y(p.y), r: 5,
    fill: p.colour, 'fill-opacity': .85, stroke: '#fff', 'stroke-width': 1});
    c.addEventListener('mousemove', e => showTip(e, `<b>${p.label}</b>` +
      (p.extra ? '<br>' + p.extra : '') + `<br>${p.x.toFixed(3)}, ${p.y.toFixed(3)}`));
    c.addEventListener('mouseleave', hideTip); sv.append(c);
    if (labels) { const t = svgEl('text', {x: X(p.x) + 7, y: Y(p.y) + 3,
      'font-size': 8.5, fill: '#444'}); t.textContent = p.label; sv.append(t); } });
  const tx = svgEl('text', {x: (ML + W - MR) / 2, y: H - 8, 'text-anchor': 'middle',
    'font-size': 11, fill: '#6b7280'}); tx.textContent = xl; sv.append(tx);
  const ty = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'});
  ty.textContent = yl; sv.append(ty);
  return sv;
}

/* ---------------------------------------------------------- CORRELATION */
function tabCorrelation(host) {
  const S = {kind: D.per10k ? 'per10k' : 'raw', tier: 'CORE'};
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Do families expand together?'), el('p', {class: 'hint'},
    'Pearson correlation across species between every pair of families. Blocks of red are families that expand and contract as a unit. This is the same information the ordination rotates, without the rotation.'));
  const kindSel = el('select', {onchange: e => { S.kind = e.target.value; draw(); }},
    ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const tierSel = el('select', {onchange: e => { S.tier = e.target.value; draw(); }},
    el('option', {value: 'CORE'}, 'CORE only'), el('option', {value: 'ALL'}, 'all families'));
  p.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel),
    el('div', {class: 'ctl'}, el('label', {}, 'tier'), tierSel)));
  const cv = el('canvas'); p.append(figure(cv, 'family_correlation'));
  const cbar = el('div', {style: 'margin-top:10px'}); p.append(cbar);
  const pairs = el('div'); p.append(pairs); host.append(p);
  function draw() {
    const fams = VARFAM.filter(f => S.tier === 'ALL' || FMETA[f].tier === 'CORE');
    const Zs = fams.map(f => zcol(S.kind, f)), n = fams.length, N = SP.length;
    const R = Array.from({length: n}, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) for (let j = i; j < n; j++) {
      let s = 0; for (let k = 0; k < N; k++) s += Zs[i][k] * Zs[j][k];
      R[i][j] = R[j][i] = s / (N - 1); }
    const ord = Array.from({length: n}, (_, i) => i)
      .sort((a, b) => R[b].reduce((x, y) => x + y, 0) - R[a].reduce((x, y) => x + y, 0));
    const CS = 20, L = 130, T = 130;
    cv.width = L + n * CS + 20; cv.height = T + n * CS + 20;
    const g = cv.getContext('2d'); g.clearRect(0, 0, cv.width, cv.height);
    g.font = '10.5px sans-serif';
    ord.forEach((fi, r) => { g.textAlign = 'right'; g.fillStyle = '#1c1f23';
      g.fillText(fams[fi], L - 6, T + r * CS + CS - 6);
      g.save(); g.translate(L + r * CS + CS - 6, T - 6); g.rotate(-Math.PI / 2);
      g.textAlign = 'left'; g.fillText(fams[fi], 0, 0); g.restore(); });
    ord.forEach((fi, r) => ord.forEach((fj, c) => {
      g.fillStyle = rdbu((R[fi][fj] + 1) / 2);
      g.fillRect(L + c * CS, T + r * CS, CS - 1, CS - 1); }));
    cv.onmousemove = e => { const b = cv.getBoundingClientRect();
      const c = Math.floor((e.clientX - b.left - L) / CS), r = Math.floor((e.clientY - b.top - T) / CS);
      if (c < 0 || r < 0 || c >= n || r >= n) return hideTip();
      showTip(e, `<b>${fams[ord[r]]}</b> vs <b>${fams[ord[c]]}</b><br>r = ${R[ord[r]][ord[c]].toFixed(3)}`); };
    cv.onmouseleave = hideTip;
    cbar.innerHTML = '';
    cbar.append(colourbar(rdbu, -1, 1, 'Pearson r', 280, 12));
    const rows = [];
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++)
      rows.push([fams[i], fams[j], +R[i][j].toFixed(3)]);
    pairs.innerHTML = ''; pairs.append(el('h2', {style: 'margin-top:16px'}, 'Pairs'));
    makeTable(pairs, ['family A', 'family B', 'r'], rows, {sortBy: 2, desc: true, name: 'family_correlations'});
  }
  draw();
}

/* -------------------------------------------------------------- DOMAINS */
function tabDomains(host) {
  if (!D.completeness) { host.append(el('div', {class: 'panel'},
    el('h2', {}, 'Domains'), el('p', {class: 'hint'},
    'No domain output found. Run: defensome.py domains --out <results>'))); return; }
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Architecture completeness'), el('p', {class: 'hint'},
    'A protein is COMPLETE when every domain its family requires is present above the family coverage threshold, using merged non-overlapping HMM segments. Families whose model matches in two pieces would look truncated if scored on the best single segment.'));
  const fams = Object.keys(D.completeness).sort((a, b) => {
    const ta = D.completeness[a], tb = D.completeness[b];
    const pa = ta.COMPLETE / Object.values(ta).reduce((x, y) => x + y, 0);
    const pb = tb.COMPLETE / Object.values(tb).reduce((x, y) => x + y, 0);
    return pa - pb; });
  const cols = {COMPLETE: '#2a9d8f', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#adb5bd'};
  const box = el('div');
  fams.forEach(f => { const cc = D.completeness[f], tot = Object.values(cc).reduce((a, b) => a + b, 0);
    const row = el('div', {style: 'display:flex;align-items:center;gap:10px;margin:3px 0'});
    row.append(el('div', {style: 'width:170px;font-size:12px'}, `${f} (${tot})`));
    const bar = el('div', {style: 'flex:1;display:flex;height:16px;border-radius:3px;overflow:hidden'});
    Object.entries(cc).forEach(([k, v]) => { if (!v) return;
      const d = el('div', {style: `width:${v / tot * 100}%;background:${cols[k]}`});
      d.addEventListener('mousemove', e => showTip(e, `<b>${f}</b><br>${k}: ${v} (${(v / tot * 100).toFixed(1)}%)`));
      d.addEventListener('mouseleave', hideTip); bar.append(d); });
    row.append(bar, el('div', {style: 'width:56px;text-align:right;font-size:12px;color:#6b7280'},
      (cc.COMPLETE / tot * 100).toFixed(1) + '%'));
    box.append(row); });
  p.append(box, el('div', {class: 'legend'}, Object.entries(cols).map(([k, v]) =>
    el('span', {}, el('i', {class: 'sw', style: 'background:' + v}), k))));
  host.append(p);

  if (D.architectures) { const p2 = el('div', {class: 'panel'});
    p2.append(el('h2', {}, 'Domain architectures'), el('p', {class: 'hint'},
      'Ordered Pfam domain strings observed per family. Unexpected architectures are usually fused or split gene models.'));
    makeTable(p2, ['family', 'architecture', 'proteins'], D.architectures,
      {sortBy: 2, desc: true, mono: [1], name: 'architectures', cap: 1200});
    host.append(p2); }

  if (D.flagged_proteins) { const p3 = el('div', {class: 'panel'});
    p3.append(el('h2', {}, 'Proteins that are not COMPLETE'), el('p', {class: 'hint'},
      'The individual gene models worth eyeballing. domain_coverages gives the merged HMM coverage per required domain.'));
    makeTable(p3, ['species', 'protein', 'family', 'status', 'length', 'length ratio',
      'min domain cov', 'domain coverages'], D.flagged_proteins,
      {sortBy: 6, mono: [1, 7], name: 'flagged_proteins', cap: 1000});
    host.append(p3); }
}

/* ---------------------------------------------------------------- CLANS */
function tabClans(host) {
  const keys = Object.keys(D.clans || {});
  if (!keys.length) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Clans'),
    el('p', {class: 'hint'}, 'No clan output found. Run: defensome.py cyp --out <results> --refs <clan-labelled FASTA>, then genetree.'))); return; }
  const p = el('div', {class: 'panel'});
  const sel = el('select', {}, ...keys.map(k => el('option', {value: k}, k)));
  const modeSel = el('select', {}, el('option', {value: 'frac'}, 'fraction of family'),
    el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`), el('option', {value: 'raw'}, 'raw counts'));
  p.append(el('h2', {}, 'Clan composition'), el('p', {class: 'hint'},
    'Fractions answer whether the balance between clans shifts; per-10k answers whether a clan expanded. They can disagree, and the disagreement is usually the interesting part.'),
    el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'family'), sel),
      el('div', {class: 'ctl'}, el('label', {}, 'scale'), modeSel)));
  const body = el('div'); p.append(body); host.append(p);
  const draw = () => { body.innerHTML = ''; renderClans(body, pick(sel, keys), pick(modeSel, ['frac','per10k','raw'])); };
  sel.onchange = draw; modeSel.onchange = draw; draw();
}

function renderClans(host, fam, mode) {
  const C = D.clans[fam];
  if (!C) { host.append(el('div', {class: 'note'}, 'No clan data for ' + fam)); return; }
  const clans = C.clans;
  const ccol = c => CCOL[clans.indexOf(c) % CCOL.length];
  const rows = SP.filter(s => C.data[s]).map(s => { const i = SPI[s], raw = C.data[s];
    const tot = raw.reduce((a, b) => a + b, 0);
    const v = mode === 'frac' ? raw.map(x => tot ? x / tot : 0)
      : mode === 'per10k' ? raw.map(x => D.species[i].n_proteins ? x / D.species[i].n_proteins * 1e4 : 0)
      : raw;
    return {sp: s, i, raw, v, tot}; });
  rows.sort((a, b) => (b.v[0] || 0) - (a.v[0] || 0));
  const W = 900, rh = 18, H = rows.length * rh + 30;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const ML = 200, MR = 30;
  const mx = Math.max(...rows.map(r => r.v.reduce((a, b) => a + b, 0))) || 1;
  rows.forEach((r, k) => { let x = ML; const y = k * rh + 6;
    r.v.forEach((v, q) => { const w = v / mx * (W - ML - MR);
      const rect = svgEl('rect', {x, y, width: Math.max(0, w), height: rh - 3, fill: ccol(clans[q])});
      rect.addEventListener('mousemove', e => showTip(e,
        `<b>${nice(r.sp)}</b><br>${clans[q]}: ${r.raw[q]} genes` +
        `<br>${(r.raw[q] / r.tot * 100).toFixed(1)}% of ${fam}` +
        `<br>${(r.raw[q] / D.species[r.i].n_proteins * 1e4).toFixed(2)} per 10k ${D.size_unit || 'proteins'}`));
      rect.addEventListener('mouseleave', hideTip); sv.append(rect); x += w; });
    const t = svgEl('text', {x: ML - 8, y: y + rh - 6, 'text-anchor': 'end', 'font-size': 11,
      'font-style': 'italic', fill: traitColour(D.species[r.i].trait)});
    t.textContent = nice(r.sp); sv.append(t); });
  host.append(figure(sv, fam + '_clans'), el('div', {class: 'legend'},
    clans.map(c => el('span', {}, el('i', {class: 'sw', style: 'background:' + ccol(c)}), c)),
    TRAITS.map(t => el('span', {style: 'color:' + traitColour(t)}, '● ' + t))));

  if (TRAITS.length) {
    const agg = {}; TRAITS.forEach(t => agg[t] = clans.map(() => 0));
    rows.forEach(r => { const t = D.species[r.i].trait; if (!t) return;
      r.raw.forEach((v, q) => agg[t][q] += v); });
    const tp = el('div', {class: 'panel'});
    tp.append(el('h2', {}, `${fam} clans by ${D.trait_name}`), el('p', {class: 'hint'},
      'Totals summed across species, then the same table as within-group fractions. If the fractions are near identical, clan balance does not track the trait, which is a real negative result.'));
    const mkT = (title, fmtf) => { const t = el('table', {},
      el('thead', {}, el('tr', {}, [D.trait_name, ...clans, 'total'].map(h => el('th', {}, h)))),
      el('tbody', {}, TRAITS.map(tr => { const tot = agg[tr].reduce((a, b) => a + b, 0);
        return el('tr', {}, el('td', {}, tr), ...agg[tr].map(v => el('td', {}, fmtf(v, tot))),
          el('td', {}, tot)); })));
      return el('div', {}, el('h2', {style: 'font-size:13px;margin:10px 0 4px'}, title), t); };
    tp.append(mkT('counts', v => v), mkT('fraction within group', (v, t) => t ? (v / t).toFixed(3) : '–'));
    host.append(tp);
  }
  const tbl = el('div', {class: 'panel'});
  tbl.append(el('h2', {}, 'Per species'));
  makeTable(tbl, ['species', D.trait_name || 'trait', ...clans, 'total'],
    rows.map(r => [r.sp, D.species[r.i].trait ?? '–', ...r.raw, r.tot]),
    {sortBy: 2, desc: true, name: fam + '_clans'});
  host.append(tbl);
}

/* ---------------------------------------------------------------- TREES */
/* Gene trees are read as evidence about the METHOD comparison: when DToL,
   genome-guided and de novo proteomes each recover a gene, their copies
   should fall together; a clade built from one method alone is either a
   transcript-only gene, a fragment family, or a contaminant, and is worth
   finding. Everything here is computed from the embedded Newick. */
function parseNewick(txt) {
  txt = txt.replace(/\[[^\]]*\]/g, '').trim().replace(/;\s*$/, '');
  let i = 0;
  function node() {
    const n = {name: '', len: 0, kids: [], support: null};
    if (txt[i] === '(') { i++;
      for (;;) { const k = node(); k.parent = n; n.kids.push(k);
        if (txt[i] === ',') i++; else break; }
      if (txt[i] === ')') i++; }
    let s = i;
    if (txt[i] === "'") { i++; while (i < txt.length && txt[i] !== "'") i++; i++; }
    else while (i < txt.length && !',():;'.includes(txt[i])) i++;
    const label = txt.slice(s, i).trim().replace(/^'|'$/g, '');
    if (txt[i] === ':') { i++; s = i;
      while (i < txt.length && /[-\d.eE+]/.test(txt[i])) i++;
      n.len = parseFloat(txt.slice(s, i)) || 0; }
    if (!n.kids.length) n.name = label;
    else if (label !== '' && !isNaN(+label)) n.support = +label;   // FastTree local support
    return n; }
  return node();
}
const walk = r => { const out = [], st = [r];
  while (st.length) { const n = st.pop(); out.push(n); st.push(...n.kids); } return out; };
const leaves = r => walk(r).filter(n => !n.kids.length);
function pruneTree(root, keep) {
  const rec = n => { if (!n.kids.length) return keep.has(n.name) ? n : null;
    const k = n.kids.map(rec).filter(Boolean);
    if (!k.length) return null;
    if (k.length === 1) { k[0].len += n.len; return k[0]; }
    n.kids = k; return n; };
  const r = rec(root); if (r) r.len = 0; return r; }

/* ---- per-tip annotation ------------------------------------------------ */
function tipInfo(name) {
  const bar = name.indexOf('|');
  const sid = bar >= 0 ? name.slice(0, bar) : name;
  const prot = bar >= 0 ? name.slice(bar + 1) : '';
  const s = (D.samples || {})[sid] || {};
  const i = SPI[sid];
  const row = i !== undefined ? D.species[i] : {};
  return {sid, prot,
    species: s.species || row.species_name || sid,
    method: s.method || row.method || '',
    pipeline: s.pipeline || row.pipeline || '',
    clan: (D.tip_clan || {})[name] || '',
    status: prot ? ((D.tip_status || {})[name] || 'COMPLETE') : '',
    trait: row.trait ?? null};
}

/* ---- palettes: every ring gets its own, so no colour means two things --- */
const OKABE = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#F0E442', '#000000'];
const RING_FIXED = {
  method: {DToL: '#1d3557', genome_guided: '#457b9d', denovo: '#e63946'},
  clan: {CYP2: '#7b2cbf', CYP3: '#f4a261', CYP4: '#2a9d8f', MITO: '#8d6e63'},
  status: {COMPLETE: '#ced4da', PARTIAL: '#e9a13b', FRAGMENT: '#bc4749', MISSING_DOMAIN: '#343a40'}};
function ringPalette(field, values) {
  const vals = [...new Set(values.filter(v => v !== '' && v != null))].sort();
  const pal = {};
  if (CVD) { vals.forEach((v, i) => pal[v] = OKABE[i % OKABE.length]); return pal; }
  const fixed = RING_FIXED[field] || {};
  let h = field === 'species' ? 12 : 200, extra = 0;
  vals.forEach(v => {
    if (fixed[v]) pal[v] = fixed[v];
    else if (field === 'trait') pal[v] = traitColour(v);
    else { pal[v] = `hsl(${(h + extra * 137.5) % 360},55%,${46 + (extra % 3) * 8}%)`; extra++; }
  });
  return pal;
}
const METHOD_ABBR = {DToL: 'DToL', genome_guided: 'GG', denovo: 'dn'};
function tipText(t, mode, isSpecies) {
  if (isSpecies) return nice(t.name);
  const i = t.info;
  if (mode === 'prot') return i.prot || t.name;
  if (mode === 'sid') return i.sid;
  if (mode === 'sp') return nice(i.species);
  return `${nice(i.species)} ${METHOD_ABBR[i.method] || i.method || ''}`.trim();
}
const RING_LABEL = {method: 'method', species: 'species', clan: 'CYP clan',
                    status: 'domain completeness', trait: () => D.trait_name || 'trait'};

/* ---- model: parse once, ladderise, index tips in drawing order --------- */
const TREE_CACHE = {};
function treeModel(key) {
  if (TREE_CACHE[key]) return TREE_CACHE[key];
  let root = parseNewick(D.trees[key]);
  const isSpecies = leaves(root).every(n => SPI[n.name] !== undefined);
  if (isSpecies) { const pr = pruneTree(root, new Set(SP)); if (pr) root = pr; }
  (function lad(n) { if (!n.kids.length) { n.nt = 1; return 1; }
    n.nt = n.kids.reduce((a, k) => a + lad(k), 0);
    n.kids.sort((a, b) => a.nt - b.nt); return n.nt; })(root);
  (function par(n, p) { n.parent = p; n.kids.forEach(k => par(k, n)); })(root, null);
  const lv = []; (function dfs(n) { if (!n.kids.length) lv.push(n); else n.kids.forEach(dfs); })(root);
  const all = []; (function po(n) { n.kids.forEach(po); all.push(n); })(root);
  (function dd(n, dep, dist) { n.dep = dep; n.dist = dist;
    n.kids.forEach(k => dd(k, dep + 1, dist + Math.max(0, k.len))); })(root, 0, 0);
  lv.forEach((n, i) => { n.idx = i; n.info = tipInfo(n.name); });
  // clade tips are a contiguous range in DFS order, so a clade is [lo, hi]
  all.forEach(n => { if (!n.kids.length) { n.lo = n.hi = n.idx; }
    else { n.lo = n.kids[0].lo; n.hi = n.kids[n.kids.length - 1].hi; } });
  let maxDep = 0, maxDist = 0;
  lv.forEach(n => { if (n.dep > maxDep) maxDep = n.dep; if (n.dist > maxDist) maxDist = n.dist; });
  const m = {key, root, lv, all, isSpecies, maxDep: maxDep || 1, maxDist: maxDist || 1,
             hasSupport: all.some(n => n.support != null),
             fields: ['method', 'species', 'clan', 'status', 'trait'].filter(f =>
               lv.some(t => t.info[f] !== '' && t.info[f] != null) &&
               !(f === 'trait' && D.trait_name === 'method'))};
  TREE_CACHE[key] = m; return m;
}

function cladeStats(m, n) {
  const tips = m.lv.slice(n.lo, n.hi + 1);
  const by = f => { const o = {}; tips.forEach(t => { const v = t.info[f] || '–'; o[v] = (o[v] || 0) + 1; }); return o; };
  return {n: tips.length, tips, method: by('method'), species: by('species'),
          clan: by('clan'), status: by('status')};
}

/* maximal clades whose every tip shares one value of `field`. For method,
   these are the clades only one annotation route produced. */
function exclusiveClades(m, field, minTips = 3) {
  const val = new Map();
  m.all.forEach(n => {
    if (!n.kids.length) { val.set(n, n.info[field] || '–'); return; }
    const vs = n.kids.map(k => val.get(k));
    val.set(n, vs.every(v => v !== null && v === vs[0]) ? vs[0] : null); });
  const out = [];
  m.all.forEach(n => {
    if (!n.kids.length) return;
    const v = val.get(n); if (v === null) return;
    const size = n.hi - n.lo + 1; if (size < minTips) return;
    if (n.parent && val.get(n.parent) === v) return;           // not maximal
    const sp = new Set(m.lv.slice(n.lo, n.hi + 1).map(t => t.info.species));
    out.push({node: n, value: v, n: size, species: sp.size, support: n.support}); });
  return out.sort((a, b) => b.n - a.n);
}

/* ---- geometry ---------------------------------------------------------- */
const TG = {S: 940, RT: 0.29, RW: 13, GAP: 2, TOP: 30, LEFT: 30, TW: 360, CW: 13};
function place(m, layout, useLen) {
  const nr = n => useLen ? n.dist / m.maxDist : (n.kids.length ? n.dep / m.maxDep : 1);
  const n = m.lv.length;
  if (layout === 'radial') {
    const span = 2 * Math.PI * 0.97, a0 = -Math.PI / 2 + Math.PI * 0.015;
    const step = span / Math.max(n, 1), RT = TG.S * TG.RT;
    m.lv.forEach(t => t.a = a0 + t.idx * step);
    m.all.forEach(x => { if (x.kids.length) x.a = (x.kids[0].a + x.kids[x.kids.length - 1].a) / 2; });
    m.all.forEach(x => { x.r = RT * nr(x); x.x = x.r * Math.cos(x.a); x.y = x.r * Math.sin(x.a); });
    return {layout, step, a0, span, RT, unit: RT / (useLen ? m.maxDist : 1)};
  }
  const H = Math.max(640, Math.min(n * 11, 2400)), dy = (H - 2 * TG.TOP) / Math.max(n - 1, 1);
  m.lv.forEach(t => t.y = TG.TOP + t.idx * dy);
  m.all.forEach(x => { if (x.kids.length) x.y = (x.kids[0].y + x.kids[x.kids.length - 1].y) / 2; });
  m.all.forEach(x => x.x = TG.LEFT + TG.TW * nr(x));
  return {layout, dy, H, unit: TG.TW / (useLen ? m.maxDist : 1)};
}
function branchPath(m, g) {
  let d = '';
  const f = v => v.toFixed(1);
  m.all.forEach(x => {
    if (!x.kids.length) return;
    if (g.layout === 'radial') {
      x.kids.forEach(k => { d += `M${f(x.r * Math.cos(k.a))},${f(x.r * Math.sin(k.a))}L${f(k.x)},${f(k.y)}`; });
      const a0 = x.kids[0].a, a1 = x.kids[x.kids.length - 1].a;
      if (x.r > 0 && a1 > a0) d += `M${f(x.r * Math.cos(a0))},${f(x.r * Math.sin(a0))}A${f(x.r)},${f(x.r)} 0 ${a1 - a0 > Math.PI ? 1 : 0},1 ${f(x.r * Math.cos(a1))},${f(x.r * Math.sin(a1))}`;
    } else {
      d += `M${f(x.x)},${f(x.kids[0].y)}L${f(x.x)},${f(x.kids[x.kids.length - 1].y)}`;
      x.kids.forEach(k => { d += `M${f(x.x)},${f(k.y)}L${f(k.x)},${f(k.y)}`; });
    }
  });
  return d;
}
function ringPaths(m, g, field, k, pal) {
  const byCol = {}, f = v => v.toFixed(1);
  m.lv.forEach(t => {
    const v = field === '__value' ? t.__v : t.info[field];
    const col = field === '__value' ? t.__c : (pal[v] || '#f1f3f5');
    let d;
    if (g.layout === 'radial') {
      const r0 = g.RT + 8 + k * (TG.RW + TG.GAP), r1 = r0 + TG.RW, h = g.step / 2 * 0.96;
      const p = (r, a) => `${f(r * Math.cos(a))},${f(r * Math.sin(a))}`;
      d = `M${p(r0, t.a - h)}A${f(r0)},${f(r0)} 0 0,1 ${p(r0, t.a + h)}L${p(r1, t.a + h)}A${f(r1)},${f(r1)} 0 0,0 ${p(r1, t.a - h)}Z`;
    } else {
      const x0 = TG.LEFT + TG.TW + 8 + k * (TG.CW + TG.GAP), h = Math.max(g.dy / 2 * 0.96, .3);
      d = `M${f(x0)},${f(t.y - h)}h${TG.CW}v${f(2 * h)}h${-TG.CW}Z`;
    }
    (byCol[col] = byCol[col] || []).push(d);
  });
  return byCol;
}

/* ---- the tab ----------------------------------------------------------- */
function tabTrees(host) {
  const keys = Object.keys(D.trees || {});
  if (!keys.length) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Trees'),
    el('p', {class: 'hint'}, 'No Newick files embedded. Run `trees`, and pass --species-tree to `dashboard` for a species tree.'))); return; }
  const gene = keys.filter(k => k !== 'species');
  const S = {key: gene.includes('CYP') ? 'CYP' : keys[0], layout: 'radial', useLen: false,
             rings: new Set(), labels: 'auto', txt: 'sm', support: false, famRing: '', query: '',
             excl: 'method', sel: null, tip: null, k: 1, tx: 0, ty: 0};

  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Trees'), el('p', {class: 'hint'},
    'Scroll to zoom, drag to pan, double-click to reset. Click a branch point to inspect that clade; click a tip for its record. ' +
    'In a gene tree, copies of one gene recovered by several methods should sit together. A clade built from a single method is a transcript-only gene, a fragment family or a contaminant, and the table below lists them.'));
  const mk = (l, n) => el('div', {class: 'ctl'}, el('label', {}, l), n);
  const selTree = el('select', {}, ...keys.map(k => el('option', {value: k}, k === 'species' ? 'species tree' : k + ' gene tree')));
  selTree.value = S.key;
  const selLay = el('select', {}, el('option', {value: 'radial'}, 'radial'), el('option', {value: 'rect'}, 'rectangular'));
  const selLen = el('select', {}, el('option', {value: 'clad'}, 'cladogram'), el('option', {value: 'phylo'}, 'phylogram'));
  const selLab = el('select', {}, el('option', {value: 'auto'}, 'labels: when readable'),
    el('option', {value: 'on'}, 'labels: always'), el('option', {value: 'off'}, 'labels: off'));
  const selTxt = el('select', {}, el('option', {value: 'sm'}, 'species + method'),
    el('option', {value: 'prot'}, 'protein ID'), el('option', {value: 'sid'}, 'sample'),
    el('option', {value: 'sp'}, 'species'));
  const search = el('input', {type: 'search', placeholder: 'find protein, sample or species'});
  const famSel = el('select', {}, el('option', {value: ''}, '(none)'), ...VARFAM.map(f => el('option', {value: f}, f)));
  const chips = el('div', {class: 'chips'});
  const supChip = el('span', {class: 'chip'}, 'support');
  p.append(el('div', {class: 'row'}, mk('tree', selTree), mk('layout', selLay), mk('branches', selLen),
    mk('labels', selLab), mk('label text', selTxt), mk('search', search),
    mk('value ring (species tree)', famSel)),
    el('div', {class: 'row'}, mk('rings', chips), mk('nodes', supChip)));
  const zbar = el('div', {class: 'row', style: 'margin:0 0 6px'});
  const zb = (t, fn) => el('button', {class: 'act', onclick: fn}, t);
  zbar.append(zb('+', () => zoomBy(1.5)), zb('−', () => zoomBy(1 / 1.5)), zb('reset', () => resetView()),
              zb('zoom to selection', () => S.sel && fitClade(S.sel)), zb('clear selection', () => { S.sel = null; S.tip = null; paint(); }));
  const stage = el('div', {style: 'display:grid;grid-template-columns:minmax(0,1fr) 300px;gap:14px;align-items:start'});
  const figHost = el('div'), side = el('div', {class: 'panel', style: 'margin:0;max-height:900px;overflow:auto'});
  stage.append(figHost, side);
  const legend = el('div', {class: 'legend', style: 'margin-top:8px'});
  const exHost = el('div', {class: 'panel'});
  p.append(zbar, stage, legend); host.append(p, exHost);

  let m, g, svg, vp, gBranch, gSel, gRing, gNode, gHit, gLab, gBar;
  const NS = 'http://www.w3.org/2000/svg';

  function build() {
    m = treeModel(S.key);
    if (!S.rings.size) (m.isSpecies ? ['trait'] : ['method', 'clan']).forEach(f => m.fields.includes(f) && S.rings.add(f));
    chips.innerHTML = '';
    m.fields.forEach(f => chips.append(el('span', {class: 'chip' + (S.rings.has(f) ? ' on' : ''),
      onclick: e => { S.rings.has(f) ? S.rings.delete(f) : S.rings.add(f); e.target.classList.toggle('on'); render(); }},
      typeof RING_LABEL[f] === 'function' ? RING_LABEL[f]() : RING_LABEL[f])));
    supChip.style.display = m.hasSupport ? '' : 'none';
    famSel.disabled = !m.isSpecies;
    S.sel = null; S.tip = null; resetView(false); render(); renderExclusive();
  }

  function render() {
    hideTip();
    g = place(m, S.layout, S.useLen);
    const rings = [...S.rings].filter(f => m.fields.includes(f));
    if (S.famRing && m.isSpecies) rings.push('__value');
    if (S.famRing && m.isSpecies) {
      const kind = D.per10k ? 'per10k' : 'raw', z = zcol(kind, S.famRing), v = col(kind, S.famRing);
      m.lv.forEach(t => { const i = SPI[t.name]; t.__v = i === undefined ? null : v[i];
        t.__c = i === undefined ? '#f1f3f5' : rdbu((z[i] + 2.5) / 5); });
    }
    // Reserve room for tip labels when they will be drawn at the starting
    // zoom, so the outermost text is never cut off at the canvas edge.
    const ringsOuter = S.layout === 'radial' ? g.RT + 8 + rings.length * (TG.RW + TG.GAP)
                                            : TG.LEFT + TG.TW + 8 + rings.length * (TG.CW + TG.GAP);
    const sp0 = S.layout === 'radial' ? (ringsOuter + 4) * g.step : g.dy;
    const labelsAt1 = S.labels === 'on' || (S.labels === 'auto' && sp0 >= 6);
    const fs0 = Math.min(11, Math.max(4, sp0 * .82));
    const maxChars = labelsAt1 ? Math.max(...m.lv.map(t => tipText(t, S.txt, m.isSpecies).length)) : 0;
    const labW = labelsAt1 ? maxChars * fs0 * .56 + 10 : 0;
    let vb;
    if (S.layout === 'radial') { const half = Math.max(TG.S / 2, ringsOuter + 4 + labW + 12); vb = [-half, -half, 2 * half, 2 * half]; }
    else vb = [0, 0, Math.max(TG.S, ringsOuter + 4 + labW + 12), g.H];
    g.vb = vb;
    const W = TG.S, H = Math.round(vb[3] * TG.S / vb[2]);
    svg = svgEl('svg', {width: W, height: H, viewBox: vb.join(' '),
      style: 'background:#fff;cursor:grab;touch-action:none;max-width:100%;height:auto'});
    vp = svgEl('g'); svg.append(vp);
    gBranch = svgEl('path', {d: branchPath(m, g), fill: 'none', stroke: '#5c636a',
      'stroke-width': m.lv.length > 800 ? .45 : .9, 'vector-effect': 'non-scaling-stroke'});
    gSel = svgEl('path', {fill: 'none', stroke: '#e63946', 'stroke-width': 2.2, 'vector-effect': 'non-scaling-stroke'});
    gRing = svgEl('g'); gNode = svgEl('g'); gHit = svgEl('g'); gLab = svgEl('g'); gBar = svgEl('g');
    vp.append(gBranch, gSel, gRing, gNode, gHit, gLab);
    svg.append(gBar);
    const pals = {};
    rings.forEach((f, k) => {
      if (f !== '__value') pals[f] = ringPalette(f, m.lv.map(t => t.info[f]));
      const byCol = ringPaths(m, g, f, k, pals[f] || {});
      Object.entries(byCol).forEach(([c, ds]) => gRing.append(svgEl('path', {d: ds.join(''), fill: c, stroke: 'none'})));
    });
    g.ringsOuter = S.layout === 'radial' ? g.RT + 8 + rings.length * (TG.RW + TG.GAP)
                                        : TG.LEFT + TG.TW + 8 + rings.length * (TG.CW + TG.GAP);
    if (S.support && m.hasSupport) m.all.forEach(x => {
      if (x.support == null || !x.kids.length) return;
      const c = x.support >= .9 ? '#1d3557' : x.support >= .7 ? '#8d99ae' : '#e63946';
      gNode.append(svgEl('circle', {cx: x.x.toFixed(1), cy: x.y.toFixed(1), r: 1.8, fill: c}));
    });
    if (S.useLen) {
      const raw = m.maxDist / 5, pow = Math.pow(10, Math.floor(Math.log10(raw)));
      const L = [1, 2, 5, 10].map(q => q * pow).find(q => q >= raw) || raw, px = L * g.unit;
      const x0 = g.vb[0] + 20, y0 = g.vb[1] + g.vb[3] - 16;
      gBar.append(svgEl('path', {d: `M${x0},${y0}h${px.toFixed(1)}M${x0},${y0 - 4}v8M${(x0 + px).toFixed(1)},${y0 - 4}v8`, stroke: '#333', 'stroke-width': 1.2, fill: 'none'}));
      const tt = svgEl('text', {x: x0 + px / 2, y: y0 - 7, 'text-anchor': 'middle', 'font-size': 10, fill: '#333'});
      tt.textContent = `${+L.toPrecision(2)} substitutions/site`; gBar.append(tt);
    }
    wire(); paint(); labels();
    figHost.innerHTML = '';
    figHost.append(figure(svg, `tree_${S.key}_${S.layout}`));
    legend.innerHTML = '';
    rings.forEach(f => {
      if (f === '__value') { legend.append(colourbar(rdbu, -2.5, 2.5, `${S.famRing}, z-score`, 220, 10)); return; }
      const pal = pals[f], lab = typeof RING_LABEL[f] === 'function' ? RING_LABEL[f]() : RING_LABEL[f];
      const cnt = {}; m.lv.forEach(t => { const v = t.info[f]; if (v) cnt[v] = (cnt[v] || 0) + 1; });
      legend.append(el('span', {style: 'font-weight:600'}, lab + ':'),
        ...Object.keys(pal).slice(0, 14).map(v => el('span', {},
          el('i', {class: 'sw', style: 'background:' + pal[v]}), `${nice(v)} (${cnt[v] || 0})`)));
      if (Object.keys(pal).length > 14) legend.append(el('span', {class: 'hint'}, `+${Object.keys(pal).length - 14} more`));
    });
    if (S.support && m.hasSupport) legend.append(el('span', {style: 'font-weight:600'}, 'support:'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#1d3557'}), '≥0.9'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#8d99ae'}), '0.7–0.9'),
      el('span', {}, el('i', {class: 'sw', style: 'background:#e63946'}), '<0.7'));
    legend.append(el('span', {class: 'hint'}, `${m.lv.length.toLocaleString()} tips`));
  }

  /* labels are rebuilt on zoom so they appear exactly when they become legible */
  function labels() {
    if (!gLab) return;
    gLab.innerHTML = '';
    const spacing = (S.layout === 'radial' ? (g.ringsOuter + 4) * g.step : g.dy) * S.k;
    if (S.labels === 'off' || (S.labels === 'auto' && spacing < 6)) return;
    const fs = Math.min(11, Math.max(4, spacing * .82)) / S.k;
    const q = S.query.toLowerCase();
    m.lv.forEach(t => {
      const txt = tipText(t, S.txt, m.isSpecies);
      const hit = q && t.name.toLowerCase().includes(q);
      let node;
      if (S.layout === 'radial') {
        const r = g.ringsOuter + 4, deg = t.a * 180 / Math.PI, flip = deg > 90 && deg < 270;
        const x = r * Math.cos(t.a), y = r * Math.sin(t.a);
        node = svgEl('text', {x: x.toFixed(1), y: y.toFixed(1), 'font-size': fs.toFixed(2),
          'text-anchor': flip ? 'end' : 'start', 'dominant-baseline': 'middle',
          transform: `rotate(${(flip ? deg + 180 : deg).toFixed(2)} ${x.toFixed(1)} ${y.toFixed(1)})`,
          fill: hit ? '#e63946' : '#212529', 'font-weight': hit ? 700 : 400});
      } else {
        node = svgEl('text', {x: (g.ringsOuter + 4).toFixed(1), y: t.y.toFixed(1), 'font-size': fs.toFixed(2),
          'dominant-baseline': 'middle', fill: hit ? '#e63946' : '#212529', 'font-weight': hit ? 700 : 400});
      }
      node.textContent = txt; gLab.append(node);
    });
  }

  /* selection highlight, search markers and the side panel */
  function paint() {
    if (!gSel) return;
    let d = '';
    if (S.sel) {
      const f = v => v.toFixed(1), sub = [];
      (function rec(x) { sub.push(x); x.kids.forEach(rec); })(S.sel);
      sub.forEach(x => { if (!x.kids.length) return;
        if (g.layout === 'radial') {
          x.kids.forEach(k => { d += `M${f(x.r * Math.cos(k.a))},${f(x.r * Math.sin(k.a))}L${f(k.x)},${f(k.y)}`; });
          const a0 = x.kids[0].a, a1 = x.kids[x.kids.length - 1].a;
          if (x.r > 0) d += `M${f(x.r * Math.cos(a0))},${f(x.r * Math.sin(a0))}A${f(x.r)},${f(x.r)} 0 ${a1 - a0 > Math.PI ? 1 : 0},1 ${f(x.r * Math.cos(a1))},${f(x.r * Math.sin(a1))}`;
        } else {
          d += `M${f(x.x)},${f(x.kids[0].y)}L${f(x.x)},${f(x.kids[x.kids.length - 1].y)}`;
          x.kids.forEach(k => { d += `M${f(x.x)},${f(k.y)}L${f(k.x)},${f(k.y)}`; }); } });
    }
    gSel.setAttribute('d', d);
    gHit.innerHTML = '';
    const q = S.query.toLowerCase(), hits = q ? m.lv.filter(t => t.name.toLowerCase().includes(q) ||
      String(t.info.species).toLowerCase().includes(q)) : [];
    hits.slice(0, 2000).forEach(t => {
      const r = S.layout === 'radial' ? g.ringsOuter + 2 : 0;
      const x = S.layout === 'radial' ? r * Math.cos(t.a) : g.ringsOuter + 2, y = S.layout === 'radial' ? r * Math.sin(t.a) : t.y;
      gHit.append(svgEl('circle', {cx: x.toFixed(1), cy: y.toFixed(1), r: (3.2 / S.k).toFixed(2), fill: '#ffb703', stroke: '#333', 'stroke-width': (.5 / S.k).toFixed(2)}));
    });
    if (S.tip) { const t = S.tip;
      const x = S.layout === 'radial' ? t.x : t.x, y = t.y;
      gHit.append(svgEl('circle', {cx: x.toFixed(1), cy: y.toFixed(1), r: (5 / S.k).toFixed(2), fill: 'none', stroke: '#e63946', 'stroke-width': (2 / S.k).toFixed(2)})); }
    renderSide(hits);
  }

  function renderSide(hits) {
    side.innerHTML = '';
    if (S.query) side.append(el('div', {class: 'note'}, `${hits.length} match${hits.length === 1 ? '' : 'es'} for “${S.query}”`),
      ...(hits.length ? [el('button', {class: 'act', onclick: () => fitTips(hits)}, 'zoom to matches')] : []));
    if (S.tip) {
      const i = S.tip.info;
      side.append(el('h2', {}, 'Tip'), el('table', {}, el('tbody', {},
        ...[['sample', i.sid], ['species', nice(i.species)], ['method', i.method], ['pipeline', i.pipeline],
            ['protein', i.prot], ['CYP clan', i.clan], ['domains', i.status],
            ['branch length', S.tip.len.toPrecision(3)]]
          .filter(r => r[1]).map(r => el('tr', {}, el('td', {style: 'color:#6b7280'}, r[0]),
            el('td', {class: r[0] === 'protein' ? 'mono' : ''}, String(r[1])))))));
    }
    if (!S.sel) {
      if (!S.tip && !S.query) side.append(el('h2', {}, 'Clade'), el('p', {class: 'hint'},
        'Click a branch point to see what that clade is made of: how many tips each method contributed, which species, and which clans.'));
      return;
    }
    const st = cladeStats(m, S.sel);
    side.append(el('h2', {}, `Clade: ${st.n} tips`),
      el('div', {class: 'hint'}, S.sel.support != null ? `support ${S.sel.support}` : ''));
    const bars = (title, obj, field) => {
      const pal = ringPalette(field, Object.keys(obj)), tot = Object.values(obj).reduce((a, b) => a + b, 0);
      const box = el('div', {style: 'margin:8px 0'}, el('div', {style: 'font-size:12px;font-weight:600'}, title));
      Object.entries(obj).sort((a, b) => b[1] - a[1]).slice(0, 12).forEach(([k, v]) => box.append(
        el('div', {style: 'display:flex;align-items:center;gap:6px;font-size:11.5px;margin:2px 0'},
          el('span', {style: 'width:92px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap'}, nice(k)),
          el('span', {style: `height:10px;width:${(v / tot * 120).toFixed(0)}px;background:${pal[k] || '#adb5bd'};border-radius:2px`}),
          el('span', {style: 'color:#6b7280'}, String(v)))));
      return box; };
    side.append(bars('by method', st.method, 'method'));
    if (Object.keys(st.clan).some(k => k !== '–')) side.append(bars('by clan', st.clan, 'clan'));
    side.append(bars('by species', st.species, 'species'));
    if (Object.keys(st.status).some(k => k !== 'COMPLETE' && k !== '–')) side.append(bars('domain completeness', st.status, 'status'));
    const methods = Object.keys(st.method).filter(k => k !== '–');
    if (methods.length === 1 && st.n >= 2) side.append(el('div', {class: 'note', style: 'border-left-color:var(--warm)'},
      `Every tip here comes from ${methods[0]}. If the other methods sampled these species, this is a ${methods[0]}-only gene, a fragment family, or a contaminant.`));
    side.append(el('div', {class: 'row'},
      el('button', {class: 'act', onclick: () => fitClade(S.sel)}, 'zoom to clade'),
      el('button', {class: 'act', onclick: () => dl(`clade_${S.key}_${st.n}tips.tsv`,
        ['sample\tspecies\tmethod\tprotein\tclan\tdomains', ...st.tips.map(t =>
          [t.info.sid, t.info.species, t.info.method, t.info.prot, t.info.clan, t.info.status].join('\t'))].join('\n'))}, 'download tips')));
  }

  function renderExclusive() {
    exHost.innerHTML = '';
    if (m.isSpecies) return;
    const fields = ['method', 'species', 'clan'].filter(f => m.fields.includes(f));
    if (!fields.length) return;
    const sel = el('select', {}, ...fields.map(f => el('option', {value: f}, f)));
    sel.value = fields.includes(S.excl) ? S.excl : fields[0];
    exHost.append(el('h2', {}, 'Single-source clades'), el('p', {class: 'hint'},
      'Maximal clades of three or more tips in which every tip shares one value. By method, these are the genes only one annotation route recovered, which is where fragments and contaminants hide. By species, they are lineage-specific expansions.'),
      el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'exclusive by'), sel)));
    const body = el('div'); exHost.append(body);
    const draw = () => { S.excl = sel.value; body.innerHTML = '';
      const ex = exclusiveClades(m, S.excl, 3);
      const tot = {}; ex.forEach(e => { tot[e.value] = (tot[e.value] || 0) + e.n; });
      body.append(el('div', {class: 'hint'}, ex.length ? `${ex.length} clades; tips in them by ${S.excl}: ` +
        Object.entries(tot).map(([k, v]) => `${nice(k)} ${v}`).join(', ') : 'none'));
      if (!ex.length) return;
      const wrap = el('div', {class: 'scroll', style: 'max-height:360px'});
      wrap.append(el('table', {}, el('thead', {}, el('tr', {},
        [S.excl, 'tips', 'species', 'support', ''].map(h => el('th', {}, h)))),
        el('tbody', {}, ex.slice(0, 400).map(e => el('tr', {},
          el('td', {}, nice(e.value)), el('td', {}, e.n), el('td', {}, e.species),
          el('td', {}, e.support == null ? '–' : e.support),
          el('td', {}, el('button', {class: 'act', onclick: () => { S.sel = e.node; S.tip = null; paint(); fitClade(e.node); }}, 'show')))))));
      body.append(wrap);
      body.append(el('div', {class: 'row'}, el('button', {class: 'act', onclick: () => dl(`single_source_${S.key}_${S.excl}.tsv`,
        [`${S.excl}\ttips\tspecies\tsupport`, ...ex.map(e => [e.value, e.n, e.species, e.support ?? ''].join('\t'))].join('\n'))}, 'download TSV')));
    };
    sel.onchange = draw; draw();
  }

  /* ---- zoom and pan -------------------------------------------------- */
  function applyView() { if (vp) vp.setAttribute('transform', `translate(${S.tx.toFixed(2)},${S.ty.toFixed(2)}) scale(${S.k.toFixed(4)})`); }
  let labTimer = null;
  function afterZoom() { applyView(); clearTimeout(labTimer); labTimer = setTimeout(() => { labels(); paint(); }, 110); }
  function resetView(redraw = true) { S.k = 1; S.tx = 0; S.ty = 0; if (redraw) { applyView(); labels(); paint(); } }
  function toModel(e) {
    const b = svg.getBoundingClientRect(), vb = svg.getAttribute('viewBox').split(' ').map(Number);
    const sx = vb[2] / (b.width || vb[2]), sy = vb[3] / (b.height || vb[3]);
    const vx = vb[0] + (e.clientX - b.left) * sx, vy = vb[1] + (e.clientY - b.top) * sy;
    return {vx, vy, x: (vx - S.tx) / S.k, y: (vy - S.ty) / S.k};
  }
  function zoomBy(f, at) {
    const c = at || {vx: g.vb[0] + g.vb[2] / 2, vy: g.vb[1] + g.vb[3] / 2};
    const k2 = Math.max(.5, Math.min(60, S.k * f)), r = k2 / S.k;
    S.tx = c.vx - (c.vx - S.tx) * r; S.ty = c.vy - (c.vy - S.ty) * r; S.k = k2; afterZoom();
  }
  function fitBox(x0, y0, x1, y1) {
    const W = g.vb[2], H = g.vb[3], pad = W * .06;
    const k = Math.max(.5, Math.min(60, Math.min((W - pad) / Math.max(x1 - x0, 1), (H - pad) / Math.max(y1 - y0, 1))));
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    const vcx = g.vb[0] + W / 2, vcy = g.vb[1] + H / 2;
    S.k = k; S.tx = vcx - cx * k; S.ty = vcy - cy * k; afterZoom();
  }
  function fitTips(tips) {
    const xs = [], ys = [];
    tips.forEach(t => { if (S.layout === 'radial') { const r = g.ringsOuter; xs.push(r * Math.cos(t.a), t.x); ys.push(r * Math.sin(t.a), t.y); }
      else { xs.push(t.x, g.ringsOuter + 120); ys.push(t.y); } });
    fitBox(Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys));
  }
  function fitClade(n) { const tips = m.lv.slice(n.lo, n.hi + 1);
    fitTips(tips.length > 3000 ? tips.filter((_, i) => i % Math.ceil(tips.length / 3000) === 0) : tips.concat([n])); }

  /* hit testing is geometric: one listener, no per-tip DOM nodes */
  function tipAt(pt) {
    if (S.layout === 'radial') {
      const r = Math.hypot(pt.x, pt.y);
      if (r < g.RT * .6 || r > g.ringsOuter + 260) return null;
      let a = Math.atan2(pt.y, pt.x); while (a < g.a0) a += 2 * Math.PI;
      const i = Math.round((a - g.a0) / g.step);
      return (i >= 0 && i < m.lv.length) ? m.lv[i] : null;
    }
    if (pt.x < TG.LEFT + TG.TW * .5) return null;
    const i = Math.round((pt.y - TG.TOP) / g.dy);
    return (i >= 0 && i < m.lv.length) ? m.lv[i] : null;
  }
  function nodeAt(pt) {
    let best = null, bd = 9 / S.k;
    m.all.forEach(x => { if (!x.kids.length) return;
      const d = Math.hypot(x.x - pt.x, x.y - pt.y); if (d < bd) { bd = d; best = x; } });
    return best;
  }
  function wire() {
    let drag = null;
    svg.addEventListener('wheel', e => { e.preventDefault(); zoomBy(e.deltaY < 0 ? 1.2 : 1 / 1.2, toModel(e)); }, {passive: false});
    svg.addEventListener('mousedown', e => { drag = {x: e.clientX, y: e.clientY, tx: S.tx, ty: S.ty, moved: false};
      svg.style.cursor = 'grabbing'; });
    window.addEventListener('mouseup', () => { if (svg) svg.style.cursor = 'grab'; setTimeout(() => drag = null, 0); });
    svg.addEventListener('mousemove', e => {
      if (drag && (e.buttons & 1)) {
        const b = svg.getBoundingClientRect(), vb = svg.getAttribute('viewBox').split(' ').map(Number);
        const s = vb[2] / (b.width || vb[2]);
        const dx = (e.clientX - drag.x) * s, dy = (e.clientY - drag.y) * s;
        if (Math.abs(dx) + Math.abs(dy) > 2) drag.moved = true;
        S.tx = drag.tx + dx; S.ty = drag.ty + dy; applyView(); hideTip(); return;
      }
      const pt = toModel(e), n = nodeAt(pt);
      if (n) { const st = cladeStats(m, n);
        showTip(e, `<b>clade, ${st.n} tips</b>` + (n.support != null ? `<br>support ${n.support}` : '') +
          '<br>' + Object.entries(st.method).map(([k, v]) => `${k} ${v}`).join(' · ') + '<br><i>click to inspect</i>'); return; }
      const t = tipAt(pt);
      if (t) { const i = t.info;
        showTip(e, `<b>${nice(i.species)}</b>${i.method ? ' · ' + i.method : ''}` +
          (i.prot ? `<br><span class="mono">${i.prot}</span>` : '') +
          (i.clan ? `<br>clan ${i.clan}` : '') + (i.status && i.status !== 'COMPLETE' ? `<br>${i.status}` : '')); }
      else hideTip();
    });
    svg.addEventListener('mouseleave', hideTip);
    svg.addEventListener('click', e => { if (drag && drag.moved) return;
      const pt = toModel(e), n = nodeAt(pt);
      if (n) { S.sel = n; S.tip = null; } else { const t = tipAt(pt); if (t) { S.tip = t; } else { S.sel = null; S.tip = null; } }
      paint(); });
    svg.addEventListener('dblclick', e => { e.preventDefault(); resetView(); });
    applyView();
  }

  selTree.onchange = () => { S.key = selTree.value; S.rings = new Set(); S.famRing = ''; famSel.value = ''; build(); };
  selLay.onchange = () => { S.layout = selLay.value; resetView(false); render(); };
  selLen.onchange = () => { S.useLen = selLen.value === 'phylo'; render(); };
  selLab.onchange = () => { S.labels = selLab.value; render(); };
  selTxt.onchange = () => { S.txt = selTxt.value; render(); };
  famSel.onchange = () => { S.famRing = famSel.value; render(); };
  supChip.onclick = e => { S.support = !S.support; e.target.classList.toggle('on'); render(); };
  let qTimer = null;
  search.oninput = e => { clearTimeout(qTimer); qTimer = setTimeout(() => { S.query = e.target.value.trim(); paint(); labels(); }, 150); };
  build();
  /* exposed for the headless tests */
  tabTrees._state = S; tabTrees._model = () => m; tabTrees._geom = () => g;
  tabTrees._tipAt = tipAt; tabTrees._nodeAt = nodeAt;
  tabTrees._setLayout = l => { S.layout = l; resetView(false); render(); };
  tabTrees._svg = () => svg; tabTrees._select = n => { S.sel = n; paint(); };
  tabTrees._zoomBy = zoomBy; tabTrees._fitClade = fitClade;
}

/* ----------------------------------------------------------- DATA TABLES */
function tabData(host) {
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Count matrices'), el('p', {class: 'hint'},
    'The matrices every other tab is built from. Sort, filter and download.'));
  const kindSel = el('select', {}, ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  const body = el('div');
  p.append(el('div', {class: 'row'}, el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel)), body);
  const draw = () => { body.innerHTML = '';
    const m = matrix(pick(kindSel, ['per10k', 'raw']));
    makeTable(body, ['species', D.trait_name || 'trait', 'proteins', ...FAM],
      SP.map((s, i) => [s, D.species[i].trait ?? '–', D.species[i].n_proteins,
        ...FAM.map(f => m[f] ? (Number.isInteger(m[f][i]) ? m[f][i] : +m[f][i].toFixed(2)) : 0)]),
      {name: 'counts_' + pick(kindSel, ['per10k', 'raw'])}); };
  kindSel.onchange = draw; draw(); host.append(p);

  if (D.gene_calls) { const p2 = el('div', {class: 'panel'});
    p2.append(el('h2', {}, 'Gene calls'), el('p', {class: 'hint'},
      `${D.gene_calls.length.toLocaleString()} called proteins. Filter by species, family or protein ID.`));
    makeTable(p2, ['species', 'family', 'protein', 'length', 'HMM coverage'],
      D.gene_calls.map(r => [SP[r[0]], FAM[r[1]], r[2], r[3], r[4]]),
      {sortBy: 0, mono: [2], name: 'gene_calls', cap: 600});
    host.append(p2); }
}

/* -------------------------------------------------------------- MAP TAB */
function tabMap(host) {
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Defensome map'), el('p', {class: 'hint'},
    'The rules that produced everything else. rule ALL means every listed domain must occur in the same protein; ANY means at least one. min_cov is the fraction of the HMM that must be aligned, summed over non-overlapping segments.'));
  makeTable(p, ['family', 'category', 'tier', 'Pfam', 'rule', 'min_cov', 'min_len',
    'median (raw)', 'notes'],
    D.families.map(f => [f.family, f.category, f.tier, f.pfam, f.rule, f.min_cov, f.min_len,
      +stats(D.counts[f.family]).med.toFixed(1), f.notes]),
    {sortBy: 0, mono: [3], name: 'defensome_map', cap: 200});
  host.append(p);
  if (D.qc_zero && D.qc_zero.length) { const p2 = el('div', {class: 'panel'});
    p2.append(el('h2', {}, 'Why each empty family is empty'), el('p', {class: 'hint'},
      'Every family with no calls at Pfam\u2019s gathering threshold, with a verdict. A zero can mean the HMM was never in the database, that a map rule filtered everything, that real copies score just below the threshold, that the input cannot contain proteins that short, or genuine absence. Those need different responses.'));
    const vc = {ACCESSION_NOT_IN_DATABASE: '#6a040f', FILTERED_BY_MAP: '#e9a13b', FOUND_BELOW_GA: '#2a9d8f',
      BELOW_GA_FAILED_VALIDATION: '#f4a261', INPUT_LACKS_SHORT_PROTEINS: '#7b2cbf',
      COMPOSITION_CANDIDATES_ONLY: '#8d99ae', NOT_DETECTED: '#6c757d', NOT_DETECTED_AT_GA_ONLY: '#6c757d'};
    const list = el('div');
    D.qc_zero.forEach(z => list.append(el('div', {style: 'display:grid;grid-template-columns:220px 1fr;gap:10px;padding:8px 0;border-bottom:1px solid var(--line)'},
      el('div', {}, el('b', {}, z.family), el('div', {class: 'mono', style: 'color:#6b7280'}, z.pfam_ids)),
      el('div', {},
        el('span', {class: 'badge', style: `background:${vc[z.verdict] || '#6c757d'};color:#fff`},
          String(z.verdict).replace(/_/g, ' ').toLowerCase().replace(/ ga\b/g, ' GA')),
        el('div', {style: 'font-size:12.5px;margin-top:4px'}, z.what_it_means || ''),
        el('div', {class: 'hint', style: 'margin:3px 0 0'},
          `GA hits ${z.raw_hits_at_GA ?? '–'} · below-GA hits ${z.raw_hits_below_GA ?? '–'} · rescued ${z.rescued_validated ?? '–'}` +
          (z['proteomes_with_<=5_proteins_under_60aa'] ? ` · proteomes lacking short proteins ${z['proteomes_with_<=5_proteins_under_60aa']}` : ''))))));
    p2.append(list);
    host.append(p2); }
  if (D.rescued && Object.keys(D.rescued).length) { const p3 = el('div', {class: 'panel'});
    p3.append(el('h2', {}, 'Found below the gathering threshold'), el('p', {class: 'hint'},
      'Genes from the rescue pass that passed the map\u2019s validation rules. They are kept out of counts.tsv and every other tab; this is the only place they are added in.'));
    const fams = Object.keys(D.rescued);
    makeTable(p3, ['sample', ...fams.flatMap(f => [f + ' at GA', f + ' rescued'])],
      SP.map((s, i) => [s, ...fams.flatMap(f => [(D.counts[f] || [])[i] ?? 0, D.rescued[f][i] ?? 0])]),
      {name: 'rescued_counts'});
    host.append(p3); }
}

/* ------------------------------------------------- FIGURE EXPORT + LEGENDS */
/* Every figure gets PNG and, where the source is SVG, vector SVG too.
   A publication figure that cannot leave the dashboard is not much use. */
function savePNG(canvas, name) {
  canvas.toBlob(b => { const u = URL.createObjectURL(b);
    const a = el('a', {href: u, download: name + '.png'}); a.click();
    setTimeout(() => URL.revokeObjectURL(u), 1000); });
}
function svgSource(svg) {
  const c = svg.cloneNode(true);
  c.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  c.setAttribute('style', 'background:#fff');
  return '<?xml version="1.0" encoding="UTF-8"?>\n' + new XMLSerializer().serializeToString(c);
}
function saveSVG(svg, name) {
  const b = new Blob([svgSource(svg)], {type: 'image/svg+xml'});
  const u = URL.createObjectURL(b);
  const a = el('a', {href: u, download: name + '.svg'}); a.click();
  setTimeout(() => URL.revokeObjectURL(u), 1000);
}
function svgToPNG(svg, name, scale = 3) {
  const w = +svg.getAttribute('width') || 800, h = +svg.getAttribute('height') || 600;
  const img = new Image();
  img.onload = () => { const c = document.createElement('canvas');
    c.width = w * scale; c.height = h * scale;
    const g = c.getContext('2d');
    g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height);
    g.setTransform(scale, 0, 0, scale, 0, 0); g.drawImage(img, 0, 0);
    savePNG(c, name); };
  img.src = 'data:image/svg+xml;base64,' + btoa(unescape(encodeURIComponent(svgSource(svg))));
}
/* wrap any figure with a download toolbar */
function figure(node, name, extra) {
  const bar = el('div', {class: 'row', style: 'margin:6px 0 0'});
  if (node.tagName === 'CANVAS') {
    bar.append(el('button', {class: 'act', onclick: () => savePNG(node, name)}, 'PNG'));
  } else {
    bar.append(el('button', {class: 'act', onclick: () => svgToPNG(node, name)}, 'PNG (3×)'),
               el('button', {class: 'act', onclick: () => saveSVG(node, name)}, 'SVG (vector)'));
  }
  if (extra) bar.append(extra);
  return el('div', {}, node, bar);
}

/* continuous colour scale legend */
function colourbar(ramp, lo, hi, label, w = 240, h = 12, ticks = 5) {
  const H = h + 34;
  const sv = svgEl('svg', {width: w, height: H, viewBox: `0 0 ${w} ${H}`});
  const N = 120;
  for (let i = 0; i < N; i++) {
    sv.append(svgEl('rect', {x: i * w / N, y: 14, width: w / N + .6, height: h,
      fill: ramp(i / (N - 1)), stroke: 'none'}));
  }
  sv.append(svgEl('rect', {x: 0, y: 14, width: w, height: h, fill: 'none', stroke: '#c9ced4'}));
  for (let k = 0; k < ticks; k++) {
    const x = w * k / (ticks - 1), v = lo + (hi - lo) * k / (ticks - 1);
    sv.append(svgEl('line', {x1: x, x2: x, y1: 14 + h, y2: 14 + h + 4, stroke: '#6b7280'}));
    const t = svgEl('text', {x: Math.min(w - 10, Math.max(10, x)), y: 14 + h + 15,
      'text-anchor': 'middle', 'font-size': 10, fill: '#6b7280'});
    t.textContent = Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(Math.abs(v) < 1 ? 2 : 1);
    sv.append(t);
  }
  const lb = svgEl('text', {x: 0, y: 9, 'font-size': 10.5, fill: '#1c1f23'});
  lb.textContent = label; sv.append(lb);
  return sv;
}

/* categorical legend with a colour-vision-safe toggle */
const CVD_SAFE = ['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#F0E442','#000000'];
traitColour = t => { if (t == null) return '#c9ced4';
  const i = TRAITS.indexOf(t);
  return CVD ? CVD_SAFE[i % CVD_SAFE.length] : TCOL[i % TCOL.length]; };

/* ---------------------------------------------------- STATISTICS + METHODS */
function tabStats(host) {
  if (!D.stats || !Object.keys(D.stats).length) {
    host.append(el('div', {class: 'panel'}, el('h2', {}, 'Statistics'),
      el('p', {class: 'hint'}, 'No group tests found. Run report with --metadata and --group-by.')));
  } else {
    for (const [grp, rows] of Object.entries(D.stats)) {
      const p = el('div', {class: 'panel'});
      p.append(el('h2', {}, `Kruskal-Wallis by ${grp}`));
      p.append(el('div', {class: 'note', style: 'border-left-color:var(--warm)'},
        'These tests treat species as independent observations. They are not. ' +
        'Diet breadth is phylogenetically conserved in Lepidoptera, so congeners ' +
        'contribute almost no independent information while inflating the degrees ' +
        'of freedom. Read this table as a screen for what to look at, never as a ' +
        'result. A defensible test needs a phylogenetic comparative model ' +
        '(PGLS, or a phylogenetic GLMM) on a time-calibrated tree, and a count of ' +
        'independent trait transitions rather than a count of tips.'));
      makeTable(p, ['family', 'H', 'p', 'p (BH-adjusted)'],
        rows.map(r => [r.family, +Number(r.H).toFixed(3), +Number(r.p).toFixed(4),
                       +Number(r.p_bh).toFixed(4)]),
        {sortBy: 2, name: 'kruskal_' + grp});
      const sig = rows.filter(r => Number(r.p_bh) < 0.05).length;
      p.append(el('div', {class: 'hint'},
        `${sig} of ${rows.length} families pass a 5% false discovery rate. ` +
        `With ${SP.length} species and strong phylogenetic structure, expect this ` +
        `to shrink under a phylogenetic model.`));
      host.append(p);
    }
  }
  for (const [grp, gm] of Object.entries(D.group_means || {})) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, `Group means by ${grp}`));
    makeTable(p, ['family', ...gm.groups],
      Object.entries(gm.data).map(([f, v]) => [f, ...v.map(x => x === null ? null : +x.toFixed(2))]),
      {sortBy: 0, name: 'group_means_' + grp});
    host.append(p);
  }
}

function tabMethods(host) {
  const p = el('div', {class: 'panel'});
  p.append(el('h2', {}, 'Provenance'), el('p', {class: 'hint'},
    'Which pipeline outputs this dashboard was built from, and which are absent. Empty rows tell you which tabs will be thin.'));
  if (D.provenance) {
    makeTable(p, ['output', 'path', 'present', 'files', 'last modified'],
      D.provenance.map(r => [r.item, r.path, r.found ? 'yes' : 'no',
        r.n_files ?? '–', r.modified ?? '–']), {sortBy: 0, name: 'provenance'});
  }
  p.append(el('div', {class: 'note'},
    `defensome ${D.tool_version || '?'} · generated ${D.generated} · source ${D.out_dir}`));
  host.append(p);

  const m = el('div', {class: 'panel'});
  m.append(el('h2', {}, 'Methods text'), el('p', {class: 'hint'},
    'A first draft of the methods paragraph, filled in from this run. Check every number before it goes near a thesis.'));
  const nCalls = FAM.reduce((a, f) => a + D.counts[f].reduce((x, y) => x + y, 0), 0);
  const inv = Object.entries(INVARIANT).map(([f, w]) => `${f} (${w})`).join('; ');
  const txt =
`Proteomes for ${SP.length} species were searched against Pfam-A with hmmsearch
using the curated per-family gathering thresholds (--cut_ga), which are used in
preference to a fixed E-value because E-values depend on database size and are
therefore not comparable across proteomes of different sizes. Domain hits were
assigned to ${FAM.length} defensome families using explicit architecture rules:
families defined by a single domain require that domain at or above a
family-specific fraction of the HMM length, and families defined by a domain
combination require every listed domain within the same protein. HMM coverage
was computed by merging non-overlapping segments, so models that match in two
pieces are not scored on their best segment alone. This yielded ${nCalls.toLocaleString()}
gene calls. Counts were normalised to copies per 10,000 predicted proteins to
remove the dependence of raw counts on proteome size and annotation
completeness; both raw and normalised matrices are reported. Low-copy families
were used as an internal calibration check${inv ? ' (' + inv + ')' : ''}.
Analyses were performed with defensome ${D.tool_version || 'x.y.z'}.`;
  m.append(el('div', {class: 'mono',
    style: 'white-space:pre-wrap;background:#f6f8fa;padding:12px;border-radius:6px'}, txt));
  m.append(el('div', {class: 'row'}, el('button', {class: 'act',
    onclick: () => dl('methods_draft.txt', txt)}, 'download methods draft')));
  host.append(m);

  const c = el('div', {class: 'panel'});
  c.append(el('h2', {}, 'Display'), el('p', {class: 'hint'},
    'The default palette is not safe for the commonest forms of colour vision deficiency. Switch to the Okabe-Ito set for figures that go into print.'));
  c.append(el('div', {class: 'row'},
    el('button', {class: 'act' + (CVD ? ' on' : ''), onclick: e => {
      CVD = !CVD; e.target.classList.toggle('on');
      Object.keys(built).forEach(k => { if (k !== 'methods') {
        const h = $('#t-' + k); h.innerHTML = ''; delete built[k]; } });
      go('methods');
    }}, 'colour-vision-safe palette')));
  const sw = el('div', {class: 'legend'});
  TRAITS.forEach(t => sw.append(el('span', {},
    el('i', {class: 'sw', style: 'background:' + traitColour(t)}), t)));
  c.append(sw);
  host.append(c);
}

/* --------------------------------------------------------- SPECIES COMPARE */
function tabCompare(host) {
  const p = el('div', {class: 'panel'});
  const a1 = el('select', {}, ...SP.map(s => el('option', {value: s}, nice(s))));
  const b1 = el('select', {}, ...SP.map(s => el('option', {value: s}, nice(s))));
  if (SP.length > 1) b1.value = SP[1];
  const kindSel = el('select', {}, ...(D.per10k ? [el('option', {value: 'per10k'}, `per 10k ${D.size_unit || 'proteins'}`)] : []),
    el('option', {value: 'raw'}, 'raw copies'));
  p.append(el('h2', {}, 'Compare two genomes'), el('p', {class: 'hint'},
    'Family-by-family difference between two species, sorted by effect. Bars point toward whichever genome has more. Check the QC flags first: a difference between an inflated and a clean annotation is not biology.'),
    el('div', {class: 'row'},
      el('div', {class: 'ctl'}, el('label', {}, 'species A'), a1),
      el('div', {class: 'ctl'}, el('label', {}, 'species B'), b1),
      el('div', {class: 'ctl'}, el('label', {}, 'values'), kindSel)));
  const body = el('div'); p.append(body); host.append(p);
  const draw = () => { body.innerHTML = '';
    const A = pick(a1, SP), B = pick(b1, SP), kind = pick(kindSel, ['per10k', 'raw']);
    const ia = SPI[A], ib = SPI[B];
    const sa = D.species[ia], sb = D.species[ib];
    body.append(el('div', {class: 'grid g2'},
      el('div', {class: 'kpi'}, el('div', {class: 'v'}, nice(A)),
        el('div', {class: 'l'}, `${sa.n_proteins.toLocaleString()} proteins · ${sa.trait ?? '–'} · QC ${sa.flag ?? '–'}`)),
      el('div', {class: 'kpi'}, el('div', {class: 'v'}, nice(B)),
        el('div', {class: 'l'}, `${sb.n_proteins.toLocaleString()} proteins · ${sb.trait ?? '–'} · QC ${sb.flag ?? '–'}`))));
    if ((sa.flag && sa.flag !== 'ok') || (sb.flag && sb.flag !== 'ok'))
      body.append(el('div', {class: 'note', style: 'border-left-color:var(--warm)'},
        'At least one of these genomes is QC-flagged. Differences below are confounded with annotation quality.'));
    const items = VARFAM.map(f => { const v = col(kind, f);
      return {label: f, value: v[ia] - v[ib],
        colour: v[ia] >= v[ib] ? '#1b4965' : '#bc4749',
        extra: `${nice(A)} ${fmt(v[ia])} vs ${nice(B)} ${fmt(v[ib])}`}; })
      .sort((x, y) => y.value - x.value);
    const sv = barChart(items, 880, Math.max(240, items.length * 16), true);
    body.append(figure(sv, `compare_${A}_vs_${B}`));
    const rows = FAM.map(f => { const v = col(kind, f);
      return [f, FMETA[f].category, +fmt(v[ia], 2), +fmt(v[ib], 2),
        +(v[ia] - v[ib]).toFixed(2),
        v[ib] ? +(v[ia] / v[ib]).toFixed(2) : null]; });
    const t = el('div', {class: 'panel'});
    t.append(el('h2', {}, 'All families'));
    makeTable(t, ['family', 'category', nice(A), nice(B), 'difference', 'ratio A/B'],
      rows, {sortBy: 4, desc: true, name: `compare_${A}_${B}`});
    body.append(t);
  };
  [a1, b1, kindSel].forEach(s => s.onchange = draw); draw();
}


/* ---------------------------------------------------- ANNOTATION METHODS */
/* The method comparison is the central result of a multi-method run, so it
   gets its own tab: one panel per question, each with the numbers beside it. */
const METHOD_ORDER = ['DToL', 'genome_guided', 'denovo'];
const methodColour = m => (CVD ? {DToL: '#0072B2', genome_guided: '#009E73', denovo: '#D55E00'}
                               : RING_FIXED.method)[m] || '#8d99ae';
const methodsIn = recs => { const s = new Set(recs.map(r => r.method));
  return [...METHOD_ORDER.filter(m => s.has(m)), ...[...s].filter(m => !METHOD_ORDER.includes(m)).sort()]; };

function groupedBars(groups, series, valueOf, W, H, ylab) {
  /* groups: x categories; series: bar within group; valueOf(g,s) -> number|null */
  const ML = 64, MB = 86, MT = 14, MR = 12;
  const vals = []; groups.forEach(g => series.forEach(s => { const v = valueOf(g, s); if (v != null) vals.push(v); }));
  const top = (Math.max(...vals, 1)) * 1.08;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = v => MT + (H - MT - MB) * (1 - v / top), gw = (W - ML - MR) / groups.length;
  for (let k = 0; k <= 4; k++) { const v = top * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(v), y2: Y(v), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(v) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#6b7280'});
    t.textContent = v >= 1000 ? (v / 1000).toFixed(v >= 10000 ? 0 : 1) + 'k' : v.toFixed(0); sv.append(t); }
  const yl = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'}); yl.textContent = ylab; sv.append(yl);
  groups.forEach((g, i) => {
    const bw = gw * .78 / series.length, x0 = ML + gw * i + gw * .11;
    series.forEach((s, j) => { const v = valueOf(g, s); if (v == null) return;
      const r = svgEl('rect', {x: x0 + j * bw, y: Y(v), width: bw - 1, height: Math.max(.5, Y(0) - Y(v)), fill: methodColour(s)});
      r.addEventListener('mousemove', e => showTip(e, `<b>${nice(g)}</b><br>${s}: ${v.toLocaleString()}`));
      r.addEventListener('mouseleave', hideTip); sv.append(r); });
    const t = svgEl('text', {x: x0 + gw * .39, y: H - MB + 14, 'font-size': 10.5, 'font-style': 'italic',
      'text-anchor': 'end', transform: `rotate(-40 ${x0 + gw * .39} ${H - MB + 14})`});
    t.textContent = nice(g); sv.append(t); });
  return sv;
}

function pairedDots(rows, species, methods, W, H, ylab) {
  /* rows: {species, method, value}; one grey line per species joins its methods */
  const ML = 64, MB = 86, MT = 14, MR = 12;
  const v = rows.map(r => r.value), top = Math.max(...v, 1) * 1.1, bot = Math.min(...v, 0);
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const Y = x => MT + (H - MT - MB) * (1 - (x - bot) / (top - bot)), gw = (W - ML - MR) / species.length;
  for (let k = 0; k <= 4; k++) { const x = bot + (top - bot) * k / 4;
    sv.append(svgEl('line', {x1: ML, x2: W - MR, y1: Y(x), y2: Y(x), stroke: '#eceff2'}));
    const t = svgEl('text', {x: ML - 8, y: Y(x) + 4, 'text-anchor': 'end', 'font-size': 10, fill: '#6b7280'});
    t.textContent = x.toFixed(0); sv.append(t); }
  const yl = svgEl('text', {x: 14, y: (MT + H - MB) / 2, 'font-size': 11, fill: '#6b7280',
    transform: `rotate(-90 14 ${(MT + H - MB) / 2})`, 'text-anchor': 'middle'}); yl.textContent = ylab; sv.append(yl);
  species.forEach((sp, i) => {
    const pts = methods.map((m, j) => { const r = rows.find(q => q.species === sp && q.method === m);
      return r ? {x: ML + gw * i + gw * (.25 + .5 * j / Math.max(methods.length - 1, 1)), y: Y(r.value), m, v: r.value} : null; })
      .filter(Boolean);
    if (pts.length > 1) sv.append(svgEl('path', {d: pts.map((p, k) => `${k ? 'L' : 'M'}${p.x},${p.y}`).join(''),
      stroke: '#c9ced4', 'stroke-width': 1.5, fill: 'none'}));
    pts.forEach(p => { const c = svgEl('circle', {cx: p.x, cy: p.y, r: 5.5, fill: methodColour(p.m), stroke: '#fff', 'stroke-width': 1.2});
      c.addEventListener('mousemove', e => showTip(e, `<b>${nice(sp)}</b><br>${p.m}: ${p.v.toFixed(0)}`));
      c.addEventListener('mouseleave', hideTip); sv.append(c); });
    const t = svgEl('text', {x: ML + gw * (i + .5), y: H - MB + 14, 'font-size': 10.5, 'font-style': 'italic',
      'text-anchor': 'end', transform: `rotate(-40 ${ML + gw * (i + .5)} ${H - MB + 14})`});
    t.textContent = nice(sp); sv.append(t); });
  return sv;
}

function recoveryHeatmap(recs, methods) {
  /* median recovery per family x method, on a log2 scale centred on parity */
  const fams = [...new Set(recs.map(r => r.family))];
  const med = {}; fams.forEach(f => methods.forEach(m => {
    const v = recs.filter(r => r.family === f && r.method === m && r.recovery != null && isFinite(r.recovery)).map(r => r.recovery);
    if (v.length) { const s = [...v].sort((a, b) => a - b); med[f + '|' + m] = s[Math.floor(s.length / 2)]; } }));
  fams.sort((a, b) => (med[a + '|' + methods[methods.length - 1]] ?? 9) - (med[b + '|' + methods[methods.length - 1]] ?? 9));
  const CW = 110, CH = 20, L = 130, T = 30;
  const W = L + methods.length * CW + 20, H = T + fams.length * CH + 10;
  const sv = svgEl('svg', {width: W, height: H, viewBox: `0 0 ${W} ${H}`});
  const col = r => rdbu(1 - (Math.max(-2, Math.min(2, Math.log2(r))) + 2) / 4);
  methods.forEach((m, j) => { const t = svgEl('text', {x: L + j * CW + CW / 2, y: T - 10, 'text-anchor': 'middle', 'font-size': 11.5, 'font-weight': 600});
    t.textContent = m; sv.append(t); });
  fams.forEach((f, i) => {
    const t = svgEl('text', {x: L - 8, y: T + i * CH + CH * .7, 'text-anchor': 'end', 'font-size': 11}); t.textContent = f; sv.append(t);
    methods.forEach((m, j) => { const r = med[f + '|' + m];
      const x = L + j * CW, y = T + i * CH;
      const rect = svgEl('rect', {x: x + 1, y: y + 1, width: CW - 2, height: CH - 2, fill: r == null ? '#f1f3f5' : col(r)});
      rect.addEventListener('mousemove', e => showTip(e, `<b>${f}</b> · ${m}<br>median recovery ${r == null ? '–' : r.toFixed(2)}`));
      rect.addEventListener('mouseleave', hideTip); sv.append(rect);
      if (r != null) { const tx = svgEl('text', {x: x + CW / 2, y: y + CH * .7, 'text-anchor': 'middle', 'font-size': 10.5,
        fill: Math.abs(Math.log2(r)) > 1.1 ? '#fff' : '#212529'}); tx.textContent = r.toFixed(2); sv.append(tx); } }); });
  return sv;
}

function tabAMethods(host) {
  const C = D.cmp || {}, long = C.long || [], rec = C.recovery || [], col = C.collapse || [];
  if (!long.length && !col.length) {
    host.append(el('div', {class: 'panel'}, el('h2', {}, 'Annotation methods'),
      el('p', {class: 'hint'}, 'No method comparison in this run. It appears when the sample sheet assigns several methods to the same species and `compare` has run.')));
    return; }
  const methods = methodsIn(long.length ? long : col);
  const species = [...new Set((long.length ? long : col).map(r => r.species))].sort();
  const ref = (rec.length && !methods.includes(rec[0].method)) ? '' : 'DToL';
  const lab = long.length && long[0] ? (D.cmp_label || 'CORE defensome genes') : '';
  host.append(el('div', {class: 'grid g3'}, kpi(species.length, 'species'), kpi(methods.length, 'annotation methods'),
    kpi((long.length ? new Set(long.map(r => r.sample_id)).size : col.length), 'proteomes')));
  host.append(el('div', {class: 'legend', style: 'margin:8px 0 14px'},
    methods.map(m => el('span', {}, el('i', {class: 'sw', style: 'background:' + methodColour(m)}), m))));

  if (col.length) { const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, '1 · Genes after isoform collapse'), el('p', {class: 'hint'},
      'Proteins become genes once isoforms are merged. Transcriptome "genes" still run above the genome annotation of the same species because assemblies split one real gene across several components; that inflates any per-gene denominator, which is why the comparison below uses counts.'));
    p.append(figure(groupedBars(species, methods, (g, s) => { const r = col.find(q => q.species === g && q.method === s); return r ? r.n_genes_out : null; },
      Math.max(960, species.length * 80), 380, 'genes'), 'collapse_by_method'));
    makeTable(p, ['sample', 'species', 'method', 'proteins in', 'genes out', 'collapse %', 'utrorf dropped'],
      col.map(r => [r.sample_id, r.species, r.method, r.n_proteins_in, r.n_genes_out, r.collapse_pct, r.dropped_utrorf ?? 0]),
      {sortBy: 1, name: 'collapse_stats'});
    host.append(p); }

  if (long.length) { const p = el('div', {class: 'panel'});
    const tot = {}; long.forEach(r => { const k = r.species + '|' + r.method; tot[k] = (tot[k] || 0) + (r.value || 0); });
    const rows = Object.entries(tot).map(([k, v]) => ({species: k.split('|')[0], method: k.split('|')[1], value: v}));
    p.append(el('h2', {}, '2 · The same genome, measured by each method'), el('p', {class: 'hint'},
      'Total CORE defensome per species under each method. A grey line joins one species, so the slope is the method effect with the genome held fixed. Flat lines mean the methods agree.'));
    p.append(figure(pairedDots(rows, species, methods, Math.max(960, species.length * 80), 400, 'CORE defensome genes'), 'paired_totals'));
    host.append(p); }

  if (rec.length) { const others = methods.filter(m => m !== ref && rec.some(r => r.method === m));
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, `3 · Which families each method recovers, relative to ${ref || 'the reference'}`),
      el('p', {class: 'hint'}, 'Median over species of (method count ÷ reference count). 1.00 is parity; red is loss, blue is excess. Families are ordered by the last column, so the worst-recovered sit at the top.'));
    const hm = recoveryHeatmap(rec, others);
    p.append(el('div', {style: 'display:flex;gap:18px;align-items:flex-start;flex-wrap:wrap'},
      figure(hm, 'recovery_heatmap'),
      el('div', {}, colourbar(t => rdbu(1 - t), -2, 2, 'log2 recovery (0 = parity)', 220, 12),
        el('div', {class: 'hint', style: 'max-width:260px;margin-top:8px'},
          'Recovery above 1 is rarely extra biology. It is usually residual isoform inflation in the transcriptome, or members the genome annotation missed; check those families in the Trees tab.'))));
    const pipes = [...new Set(rec.map(r => r.ref_pipeline).filter(p => p && p !== 'unknown'))];
    if (pipes.length > 1) {
      const cells = pipes.map(pp => [pp, ...others.map(m => { const v = rec.filter(r => r.ref_pipeline === pp && r.method === m && r.recovery != null).map(r => r.recovery).sort((a, b) => a - b);
        return v.length ? +v[Math.floor(v.length / 2)].toFixed(3) : null; })]);
      p.append(el('h2', {style: 'margin-top:16px'}, 'Split by the pipeline that produced the reference'),
        el('p', {class: 'hint'}, 'If these rows differ, "recovery against DToL" averages over two different denominators and should be reported separately.'));
      makeTable(p, ['reference pipeline', ...others], cells, {name: 'recovery_by_pipeline'}); }
    host.append(p); }

  if ((C.tests || []).length) { const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, '4 · Which differences survive testing'), el('p', {class: 'hint'},
      'Wilcoxon signed-rank tests paired within species, Benjamini-Hochberg adjusted. Pairing within species removes phylogeny from this comparison: the same genome is measured twice.'));
    const nPairs = Math.min(...C.tests.map(r => +r.n_pairs));
    const floor = 2 / Math.pow(2, nPairs);
    p.append(el('div', {class: 'note', style: floor >= 0.05 ? 'border-left-color:var(--warm)' : ''},
      `With ${nPairs} pairs, the smallest p-value a two-sided signed-rank test can return is 2/2^${nPairs} = ${floor.toPrecision(2)}.` +
      (floor >= 0.05 ? ' That is above 0.05, so no difference can reach significance here whatever the data; add species before reading this table.' : '')));
    makeTable(p, ['family', 'method A', 'method B', 'pairs', 'median A', 'median B', 'p', 'p (BH)'],
      C.tests.map(r => [r.family, r.method_a, r.method_b, r.n_pairs, +(+r.median_a).toFixed(2), +(+r.median_b).toFixed(2),
        +(+r.p).toExponential(2), +(+r.p_bh).toExponential(2)]), {sortBy: 7, name: 'paired_tests'});
    host.append(p); }

  /* small proteins: where metallothioneins and MATE went */
  const so = D.short_orfs || [], comp = D.mt_composition || {}, resc = D.rescued || {};
  if (so.length || Object.keys(comp).length || resc.Metallothionein || resc.MATE) {
    const p = el('div', {class: 'panel'});
    p.append(el('h2', {}, '5 · Small proteins, and why they go missing'), el('p', {class: 'hint'},
      'Insect metallothioneins are 40–64 residues and no Pfam model covers them. TransDecoder keeps ORFs of at least 100 residues by default, so a transcriptome proteome cannot contain them. Each column is an independent line of evidence; none of them is merged into counts.tsv.'));
    const idx = s => SPI[s];
    const rows = SP.map((s, i) => { const r = so.find(q => q.sample_id === s);
      return [s, D.species[i].method || '', (D.counts.Metallothionein || [])[i] ?? 0,
        resc.Metallothionein ? (resc.Metallothionein[i] ?? 0) : '–', comp[s] ?? 0,
        r ? r.mt_like_genes : '–', (D.counts.MATE || [])[i] ?? 0, resc.MATE ? (resc.MATE[i] ?? 0) : '–']; });
    makeTable(p, ['sample', 'method', 'MT at Pfam GA', 'MT rescued below GA', 'MT by composition',
      'MT from transcripts (short-orfs)', 'MATE at GA', 'MATE rescued'], rows, {sortBy: 1, name: 'small_proteins'});
    if (so.length) {
      const bm = {}; so.forEach(r => { (bm[r.method] = bm[r.method] || []).push(r.mt_like_genes); });
      p.append(el('div', {class: 'note'}, 'Metallothionein-like genes recovered from transcripts, median by method: ' +
        Object.entries(bm).map(([m, v]) => { const s = [...v].sort((a, b) => a - b); return `${m} ${s[Math.floor(s.length / 2)]}`; }).join(' · ') +
        '. Composition-screen candidates: confirm a few by alignment before citing counts.')); }
    host.append(p); }
}

/* ------------------------------------------------------------------ BOOT */
const TABS = [['overview', 'Overview', tabOverview], ['heatmap', 'Heatmap', tabHeatmap],
  ['family', 'Families', tabFamily], ['species', 'Species', tabSpecies],
  ['ordination', 'Ordination', tabOrdination], ['correlation', 'Correlations', tabCorrelation],
  ['domains', 'Domains', tabDomains], ['clans', 'Clans', tabClans],
  ['amethods', 'Annotation methods', tabAMethods],
  ['trees', 'Trees', tabTrees], ['compare', 'Compare', tabCompare],
  ['stats', 'Statistics', tabStats], ['data', 'Data', tabData], ['map', 'Map', tabMap],
  ['methods', 'Methods', tabMethods]];
const built = {};
function go(id) {
  TABS.forEach(([k]) => { $('#t-' + k).classList.toggle('on', k === id); });
  [...$('#nav').children].forEach(b => b.classList.toggle('on', b.dataset.k === id));
  if (!built[id]) { const host = $('#t-' + id);
    try { TABS.find(t => t[0] === id)[2](host); }
    catch (e) { host.append(el('div', {class: 'panel'}, el('h2', {}, 'Error'),
      el('pre', {class: 'mono'}, String(e && e.stack || e)))); }
    built[id] = 1; }
  location.hash = id;
}
TABS.forEach(([k, label]) => $('#nav').append(
  el('button', {'data-k': k, onclick: () => go(k)}, label)));
$('#hdr').textContent = (D.samples && Object.keys(D.samples).length
  ? `${SP.length} proteomes · ${new Set(Object.values(D.samples).map(s => s.species)).size} species · ${new Set(Object.values(D.samples).map(s => s.method)).size} methods`
  : `${SP.length} species`) + ` · ${FAM.length} families · generated ${D.generated}`;
go((location.hash || '#overview').slice(1));
