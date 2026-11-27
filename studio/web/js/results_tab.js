// Results workspace: run list + views (summary, modes 3-D, dCp, flutter, time
// responses, curves, TSD, figures, text files). All data comes from output/**
// through /api/result (FLAAP result text files) and /api/file (figures, raw text).

import { api, fileUrl } from './api.js';
import { el, $, $$, clear, toast, fmt, seriesColor } from './util.js';
import { Plot, linkPlots } from './plot.js';
import { View3D } from './view3d.js';
import { diverging, sequential, colorbar } from './colormap.js';

const R = { runs: [], run: null, view: 'summary', cleanup: [], sub: {} };
let hooks = {};

// ---------------------------------------------------------------- helpers
const fpath = (run, f) => `${run.path}/${f}`;
const getJson = (run, f) => api.result(fpath(run, f), run.mtime);
const isFlutterJson = (f) => /^flutter_.*\.txt$/i.test(f);

/** Flutter result files: from results.txt artifacts (current contract), else glob. */
function flutterFiles(run) {
  const art = (run.artifacts || []).filter((a) => a.kind === 'flutter' && /\.txt$/i.test(a.file)).map((a) => a.file);
  if (art.length) return art.filter((f) => run.files.includes(f));
  return run.files.filter(isFlutterJson);
}
const byPrefix = (run, rx) => run.files.filter((f) => rx.test(f));
const textFiles = (run) => run.files.filter((f) => /\.(txt|csv)$/i.test(f));

function onCleanup(fn) { R.cleanup.push(fn); }
function runCleanup() { for (const fn of R.cleanup.splice(0)) { try { fn(); } catch { /* ignore */ } } }

function plotBox(parent, title, cls = '') {
  const box = el('div.plotbox' + cls, {}, title ? el('div.ptitle', { text: title }) : null);
  parent.append(box);
  return box;
}
function newPlot(box, opts) { const p = new Plot(box, opts); onCleanup(() => p.destroy()); return p; }

function stat(label, value, unit = '', hl = false) {
  return el('div.stat' + (hl ? '.hl' : ''), {}, el('div.l', { text: label }), el('div', {}, el('span.v', { text: value }), unit ? el('span.u', { text: unit }) : null));
}
function ctrlRow(label, ctrl) { return el('div.form-row', {}, el('label', { text: label }), ctrl); }
function selectEl(options, value, onchange) {
  const s = el('select', {}, options.map(([v, t]) => el('option', { value: String(v), text: t, selected: String(v) === String(value) })));
  s.addEventListener('change', () => onchange(s.value));
  return s;
}
function segEl(options, value, onchange) {
  const seg = el('div.seg');
  for (const [v, t] of options) {
    const b = el('button' + (String(v) === String(value) ? '.active' : ''), { text: t });
    b.addEventListener('click', () => { for (const x of seg.children) x.classList.toggle('active', x === b); onchange(v); });
    seg.append(b);
  }
  return seg;
}
const loading = (box, what) => box.append(el('div.empty-note', { text: `loading ${what}…` }));

// ---------------------------------------------------------------- views registry
const VIEWS = [
  { id: 'summary', label: 'Summary', has: () => true, render: renderSummary },
  { id: 'modes', label: 'Modes 3D', has: (r) => r.files.includes('modes.txt'), render: renderModes },
  { id: 'cp', label: 'ΔCp', has: (r) => r.files.includes('cp.txt') && r.files.includes('aero_mesh.txt'), render: renderCp },
  { id: 'flutter', label: 'Flutter', has: (r) => flutterFiles(r).length > 0, render: renderFlutter },
  { id: 'time', label: 'Time', has: (r) => byPrefix(r, /^timehist_.*\.txt$/i).length > 0, render: renderTime },
  { id: 'curves', label: 'Curves', has: (r) => byPrefix(r, /^curves_.*\.txt$/i).length > 0, render: renderCurves },
  { id: 'tsd', label: 'TSD', has: (r) => r.files.includes('tsd_steady.txt'), render: renderTsd },
  { id: 'figures', label: 'Figures', has: (r) => r.figures.length > 0, render: renderFigures },
  { id: 'files', label: 'Files', has: (r) => textFiles(r).length > 0, render: renderFiles },
];
const KIND_VIEW = { modes: 'modes', cp: 'cp', aero_mesh: 'cp', flutter: 'flutter', timehist: 'time', curves: 'curves',
  tsd: 'tsd', text: 'files', config: 'files', csv: 'files', figure: 'figures' };

// ---------------------------------------------------------------- run list
export async function loadRuns(selectName) {
  try {
    const r = await api.outputs();
    R.runs = r.outputs.sort((a, b) => {
      const sa = a.name.startsWith('studio_'), sb = b.name.startsWith('studio_');
      if (sa !== sb) return sa ? -1 : 1;
      return sa ? b.mtime - a.mtime : a.name.localeCompare(b.name);
    });
    renderRunList();
    if (selectName) await openRun(selectName);
    else if (R.run) {
      const again = R.runs.find((x) => x.name === R.run.name);
      if (again && again.mtime !== R.run.mtime) await openRun(again.name, R.view);
    }
  } catch (e) { toast(`outputs: ${e.message}`, true); }
}

function renderRunList() {
  const box = clear($('#runList'));
  const f = $('#runFilter').value.trim().toLowerCase();
  let last = null;
  for (const r of R.runs) {
    const hay = `${r.name} ${r.title || ''} ${r.model || ''} ${r.mode || ''}`.toLowerCase();
    if (f && !hay.includes(f)) continue;
    const grp = r.name.startsWith('studio_') ? 'studio runs' : 'output/';
    if (grp !== last) { box.append(el('div.grp', { text: grp })); last = grp; }
    const ok = r.status === 'ok';
    box.append(el('div.item' + (R.run && R.run.name === r.name ? '.active' : ''), {
      title: r.title || r.name, onclick: () => openRun(r.name),
    }, el('div', {}, r.name, el('span.chip' + (ok ? '.ok' : '.error'), { text: r.status || '?' })),
    el('div.sub', { text: `${r.model || '?'} / ${r.mode || '?'}${r.created ? ' · ' + r.created.slice(5, 16) : ''}` })));
  }
  if (!box.childElementCount) box.append(el('div.empty-note', { style: { padding: '6px 8px' }, text: 'no runs in output/' }));
}

export async function openRun(name, view) {
  const run = R.runs.find((r) => r.name === name);
  if (!run) { toast(`run ${name} not found`, true); return; }
  R.run = run;
  R.sub = {};
  renderRunList();
  $('#runTitle').textContent = run.title || run.name;
  $('#runPath').textContent = run.path;
  // view buttons
  const seg = clear($('#viewSeg'));
  const avail = VIEWS.filter((v) => v.has(run));
  for (const v of avail) {
    seg.append(el('button' + (v.id === view ? '.active' : ''), { text: v.label, 'data-view': v.id, onclick: () => showView(v.id) }));
  }
  // log
  const lg = $('#runLog');
  lg.textContent = '';
  if (run.files.includes('log.txt')) api.text(fpath(run, 'log.txt')).then((t) => { lg.textContent = t; }).catch(() => {});
  const want = avail.some((v) => v.id === view) ? view : 'summary';
  await showView(want);
}

export async function showView(id, sub) {
  runCleanup();
  R.view = id;
  if (sub) Object.assign(R.sub, sub);
  for (const b of $$('#viewSeg button')) b.classList.toggle('active', b.dataset.view === id);
  const main = clear($('#resultView'));
  const ctrls = clear($('#viewCtrls'));
  const v = VIEWS.find((x) => x.id === id);
  $('#ctrlTitle').textContent = v ? v.label : 'View';
  if (hooks.onRoute && R.run) hooks.onRoute(R.run.name, id);
  if (!R.run || !v) return;
  try { await v.render(R.run, main, ctrls); } catch (e) {
    console.error(e);
    main.append(el('div.card', {}, el('h2', { text: 'Could not render view' }), el('pre.text', { text: String(e.stack || e) })));
  }
}

// ---------------------------------------------------------------- summary
function flatten(obj, prefix = '', out = []) {
  for (const [k, v] of Object.entries(obj || {})) {
    const key = prefix ? `${prefix}.${k}` : k;
    if (v && typeof v === 'object' && !Array.isArray(v)) flatten(v, key, out);
    else out.push([key, v]);
  }
  return out;
}
const showVal = (v) => (Array.isArray(v) ? (v.length > 12 ? `[${v.slice(0, 12).map((x) => fmt(x, 5)).join(', ')}, … ${v.length}]` : `[${v.map((x) => fmt(x, 5)).join(', ')}]`)
  : typeof v === 'number' ? fmt(v, 6) : String(v));

async function renderSummary(run, main, ctrls) {
  let res = null;
  try { res = await getJson(run, 'results.txt'); } catch (e) { main.append(el('div.card', { text: e.message })); return; }
  const ok = res.status === 'ok';
  const card = el('div.card', {},
    el('h2', { text: res.title || run.name }),
    el('div.dim', {}, `${res.model} / ${res.mode}`, el('span.chip' + (ok ? '.ok' : '.error'), { text: res.status }),
      `  ·  created ${res.created || '–'}  ·  ${fmt(res.elapsed_s, 3)} s`));
  if (res.error) card.append(el('pre.text', { style: { marginTop: '8px', color: '#b84d48' }, text: res.error }));
  main.append(card);
  const sm = res.summary || {};
  const tiles = el('div.stats');
  // critical points (any summary entry with a velocity: flutter_pk / flutter_k / ..., old *_critical)
  const crit = Object.entries(sm).filter(([, v]) => v && typeof v === 'object' && !Array.isArray(v) && 'velocity' in v);
  for (const [k, v] of crit) {
    const label = k.replace(/^flutter_/, '').replace(/_critical$/, '');
    tiles.append(stat(`${v.kind || 'critical'} · ${label}`, v.velocity === null ? 'none' : fmt(v.velocity, 5), v.velocity === null ? '' : 'm/s', true));
    if (v.freq_hz !== undefined && v.freq_hz !== null) tiles.append(stat(`frequency · ${label}`, fmt(v.freq_hz, 4), 'Hz'));
  }
  if (Array.isArray(sm.frequencies_hz)) sm.frequencies_hz.slice(0, 6).forEach((f, i) => tiles.append(stat(`mode ${i + 1}`, fmt(f, 4), 'Hz')));
  for (const [k, v] of Object.entries(sm)) if (typeof v === 'number' && tiles.childElementCount < 14 && k !== 'n_boxes') tiles.append(stat(k, fmt(v, 5)));
  if (tiles.childElementCount) main.append(tiles);
  const rows = flatten(sm);
  if (rows.length) {
    main.append(el('div.card', {}, el('h4', { text: 'summary' }),
      el('table.kv', {}, rows.map(([k, v]) => el('tr', {}, el('td', { text: k }), el('td.num', { text: showVal(v) }))))));
  }
  // artifacts
  const arts = res.artifacts || [];
  ctrls.append(el('div.desc', { text: `${arts.length} artifacts — click to open` }));
  const t = el('table.grid');
  t.append(el('tr', {}, el('th', { text: 'file' }), el('th', { text: 'kind' })));
  for (const a of arts) {
    const view = KIND_VIEW[a.kind];
    const okView = view && VIEWS.find((v) => v.id === view && v.has(run));
    t.append(el('tr' + (okView ? '.click' : ''), {
      title: a.title,
      onclick: okView ? () => {
        if (a.kind === 'figure') showLightbox(fpath(run, a.file), a.title);
        else showView(view, a.kind === 'flutter' ? { flutter: a.file } : a.kind === 'timehist' ? { time: a.file } : a.kind === 'curves' ? { curves: a.file } : /\.(txt|csv)$/.test(a.file) ? { file: a.file } : {});
      } : null,
    }, el('td', { text: a.file }), el('td', { text: a.kind })));
  }
  ctrls.append(t);
}

// ---------------------------------------------------------------- modes 3-D
async function renderModes(run, main, ctrls) {
  loading(main, 'modes');
  const md = await getJson(run, 'modes.txt');
  const mesh = run.files.includes('aero_mesh.txt') ? await getJson(run, 'aero_mesh.txt').catch(() => null) : null;
  clear(main);
  const nm = md.freq_hz.length;
  const st = R.sub.modes || (R.sub.modes = { mode: 0, amp: 0.12, play: true, aero: true, struct: true, undeformed: true, mirror: false, color: 'disp', persp: false });
  const hasAero = !!(mesh && md.aero_h);
  const sym = mesh ? String(mesh.symmetry || 'none') : 'none';
  const wrap = el('div.fill');
  const vbox = el('div.grow');
  wrap.append(vbox);
  main.append(wrap);
  const v3 = new View3D(vbox, {
    perspective: st.persp,
    onHover: (f) => f.data || '',
  });
  onCleanup(() => v3.destroy());
  const cb = colorbar();
  vbox.append(cb.el);
  cb.set(0, 1, '|w| / max (aero boxes)', 'seq', (x) => x.toFixed(1));

  // reference length for amplitude
  const all = [...md.nodes];
  if (hasAero) for (const c of mesh.corners) all.push(...c);
  let lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
  for (const p of all) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[k]); hi[k] = Math.max(hi[k], p[k]); }
  const L = Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]) || 1;
  const nodesOnly = !md.quads && !md.tris && !md.edges;

  const build = (phase) => {
    const faces = [], lines = [], points = [];
    const m = st.mode, a = st.amp * L * phase;
    const signs = st.mirror ? [1, -1] : [1];
    for (const sg of signs) {
      const ysg = (p) => [p[0], sg * p[1], p[2]];
      const hs = sg < 0 && /^anti/i.test(sym) ? -1 : 1;
      if (st.aero && hasAero) {
        const H = md.aero_h[m];
        mesh.corners.forEach((c, i) => {
          const h = H[i];
          const pts = c.map((p, k) => ysg([p[0], p[1], p[2] + hs * a * h[k]]));
          const mag = (Math.abs(h[0]) + Math.abs(h[1]) + Math.abs(h[2]) + Math.abs(h[3])) / 4;
          const fill = st.color === 'disp' ? sequential(mag) : '#cfe3d8';
          faces.push({ pts, fill, stroke: 'rgba(18,40,32,0.25)', lineWidth: 0.5, shade: true, alpha: sg < 0 ? 0.55 : 0.95,
            data: `box ${i + 1} (${mesh.surfaces[mesh.surface_id[i]] || ''})  w/max ${fmt(mag * Math.sign(h[0] + h[1] + h[2] + h[3]), 3)}` });
          if (st.undeformed && a !== 0) faces.push({ pts: c.map(ysg), fill: null, stroke: 'rgba(125,136,132,0.35)', lineWidth: 0.5 });
        });
      }
      if (st.struct) {
        const U = md.shapes[m];
        const P = md.nodes.map((p, i) => ysg([p[0] + a * U[i][0], p[1] + a * U[i][1], p[2] + a * U[i][2]]));
        const P0 = md.nodes.map(ysg);
        for (const q of md.quads || []) {
          faces.push({ pts: q.map((i) => P[i]), fill: '#9fbfd6', stroke: '#3f6f91', lineWidth: 0.6, shade: true, alpha: st.aero && hasAero ? 0.55 : 0.95, data: `element [${q.join(', ')}]` });
        }
        for (const t of md.tris || []) {
          faces.push({ pts: t.map((i) => P[i]), fill: '#9fbfd6', stroke: '#3f6f91', lineWidth: 0.5, shade: true, alpha: 0.9, data: `tri [${t.join(', ')}]` });
        }
        for (const [i, j] of md.edges || []) {
          lines.push({ a: P[i], b: P[j], color: '#3f6f91', width: 3, bias: 0.05 });
          if (st.undeformed && a !== 0) lines.push({ a: P0[i], b: P0[j], color: 'rgba(125,136,132,0.6)', width: 1, dash: [4, 3] });
        }
        if (nodesOnly || (md.edges && md.nodes.length < 400)) P.forEach((p) => points.push({ p, color: '#3f6f91', r: nodesOnly ? 2.5 : 2 }));
      }
    }
    return { faces, lines, points };
  };
  let first = true;
  const draw = (phase) => {
    v3.setScene(build(phase));
    if (first) { v3.setScene(build(0)); v3.fit(); v3.setScene(build(phase)); first = false; }
    v3.overlay.textContent = `${md.name || ''} · ${md.labels[st.mode]} · ${fmt(md.freq_hz[st.mode], 5)} Hz` + (st.play ? '' : ' (peak)');
    v3.render();
  };
  // animation
  let raf = 0; const t0 = performance.now();
  const tick = () => {
    if (!st.play) return;
    draw(Math.sin(2 * Math.PI * (performance.now() - t0) / 1600));
    raf = requestAnimationFrame(tick);
  };
  onCleanup(() => cancelAnimationFrame(raf));
  const restart = () => { cancelAnimationFrame(raf); if (st.play) tick(); else draw(1); };

  // controls
  const tbl = el('table.grid');
  tbl.append(el('tr', {}, el('th', { text: '#' }), el('th', { text: 'f [Hz]' }), el('th', { text: 'gen. mass' })));
  const rows = [];
  for (let i = 0; i < nm; i++) {
    const tr = el('tr.click' + (i === st.mode ? '.active' : ''), {
      title: md.labels[i],
      onclick: () => { st.mode = i; rows.forEach((r, k) => r.classList.toggle('active', k === i)); restart(); },
    }, el('td', { text: md.labels[i] }), el('td', { text: fmt(md.freq_hz[i], 5) }), el('td', { text: fmt(md.gen_mass[i], 4) }));
    rows.push(tr); tbl.append(tr);
  }
  const amp = el('input', { type: 'range', min: 0, max: 0.4, step: 0.005, value: st.amp });
  const ampLbl = el('span.mono.small', { text: `${(st.amp * 100).toFixed(1)} % L` });
  amp.addEventListener('input', () => { st.amp = Number(amp.value); ampLbl.textContent = `${(st.amp * 100).toFixed(1)} % L`; if (!st.play) draw(1); });
  const playBtn = el('button.btn.small', { text: st.play ? '❚❚ Pause' : '▶ Play' });
  playBtn.addEventListener('click', () => { st.play = !st.play; playBtn.textContent = st.play ? '❚❚ Pause' : '▶ Play'; restart(); });
  const chk = (key, label, disabled) => {
    const c = el('input', { type: 'checkbox', checked: st[key], disabled });
    c.addEventListener('change', () => { st[key] = c.checked; if (key === 'persp') v3.opts.perspective = c.checked; if (!st.play) draw(1); });
    return el('label.chk', {}, c, label);
  };
  ctrls.append(tbl,
    el('h4', { text: 'Animation' }),
    el('div.row.gap', {}, playBtn, amp, ampLbl),
    el('h4', { text: 'Show' }),
    el('div.row.gap.wrap', {}, chk('aero', 'aero boxes', !hasAero), chk('struct', 'structure'), chk('undeformed', 'undeformed'),
      chk('mirror', `mirror (${sym})`), chk('persp', 'perspective')),
    ctrlRow('aero colour', selectEl([['disp', '|w| (sequential)'], ['plain', 'plain']], st.color, (v) => { st.color = v; if (!st.play) draw(1); })),
    el('h4', { text: 'View' }),
    segEl([['iso', 'Iso'], ['top', 'Top'], ['front', 'Front'], ['side', 'Side']], '', (v) => v3.setView(v)),
    el('div.desc', { style: { marginTop: '8px' }, text: `drag: orbit · shift/right-drag: pan · wheel: zoom · double-click: fit. Structure: ${md.structure || '?'}, ${md.nodes.length} nodes${md.quads ? `, ${md.quads.length} quads` : ''}${md.edges ? `, ${md.edges.length} beam edges` : ''}${nodesOnly ? ' (no connectivity: nodes only)' : ''}.` }));
  if (!hasAero) cb.el.hidden = true;
  draw(1);
  restart();
}

// ---------------------------------------------------------------- dCp
async function renderCp(run, main, ctrls) {
  loading(main, 'ΔCp');
  const [cp, mesh] = await Promise.all([getJson(run, 'cp.txt'), getJson(run, 'aero_mesh.txt')]);
  clear(main);
  const st = R.sub.cp || (R.sub.cp = { k: Math.min(1, cp.k.length - 1), mode: 0, part: 're', mirror: false, view: 'top' });
  if (st.k >= cp.k.length) st.k = 0;
  if (st.mode >= cp.modes.length) st.mode = 0;
  const wrap = el('div.fill'); const vbox = el('div.grow'); wrap.append(vbox); main.append(wrap);
  let vals = [];
  const v3 = new View3D(vbox, { view: st.view, onHover: (f) => (f.data === undefined ? '' : `box ${f.data + 1} (${mesh.surfaces[mesh.surface_id[f.data]] || ''})  ΔCp ${partLabel()} = ${fmt(vals[f.data], 4)}`) });
  onCleanup(() => v3.destroy());
  const cb = colorbar(); vbox.append(cb.el);
  const statsBox = el('div.derived');
  const partLabel = () => ({ re: 'Re', im: 'Im', abs: '|·|' }[st.part]);
  const ctrlNames = Object.entries(mesh.controls || {}).filter(([, m]) => Array.isArray(m));
  let fitted = false;
  const draw = () => {
    const re = cp.re[st.k][st.mode], im = cp.im[st.k][st.mode];
    vals = st.part === 're' ? re : st.part === 'im' ? im : re.map((r, i) => Math.hypot(r, im[i]));
    let lim = 0; for (const x of vals) if (Number.isFinite(x)) lim = Math.max(lim, Math.abs(x));
    lim = lim || 1;
    const faces = [];
    const signs = st.mirror ? [1, -1] : [1];
    for (const sg of signs) {
      mesh.corners.forEach((c, i) => {
        const t = st.part === 'abs' ? vals[i] / lim : vals[i] / lim;
        faces.push({ pts: c.map((p) => [p[0], sg * p[1], p[2]]), fill: st.part === 'abs' ? sequential(t) : diverging(t),
          stroke: 'rgba(23,33,30,0.18)', lineWidth: 0.5, data: sg > 0 ? i : undefined, alpha: sg > 0 ? 1 : 0.6 });
      });
    }
    v3.setScene({ faces });
    if (!fitted) { v3.fit(); fitted = true; }
    v3.render();
    if (st.part === 'abs') cb.set(0, lim, `|ΔCp| · k = ${fmt(cp.k[st.k], 4)} · ${cp.modes[st.mode]}`, 'seq', (x) => fmt(x, 3));
    else cb.set(-lim, lim, `${partLabel()} ΔCp · k = ${fmt(cp.k[st.k], 4)} · ${cp.modes[st.mode]}`, 'div', (x) => fmt(x, 3));
    v3.overlay.textContent = `Mach ${fmt(cp.mach, 3)} · ${mesh.corners.length} boxes · symmetry ${mesh.symmetry}`;
    let mn = Infinity, mx = -Infinity; for (const x of vals) { mn = Math.min(mn, x); mx = Math.max(mx, x); }
    clear(statsBox).append('min', el('b', { text: fmt(mn, 4) }), 'max', el('b', { text: fmt(mx, 4) }), 'k', el('b', { text: fmt(cp.k[st.k], 4) }),
      'b_ref', el('b', { text: `${fmt(cp.bref, 4)} m` }));
  };
  ctrls.append(
    ctrlRow('reduced freq.', selectEl(cp.k.map((k, i) => [i, `k = ${fmt(k, 4)}`]), st.k, (v) => { st.k = Number(v); draw(); })),
    ctrlRow('mode / input', selectEl(cp.modes.map((m, i) => [i, m]), st.mode, (v) => { st.mode = Number(v); draw(); })),
    ctrlRow('part', segEl([['re', 'Re'], ['im', 'Im'], ['abs', '|ΔCp|']], st.part, (v) => { st.part = v; draw(); })),
    (() => { const c = el('input', { type: 'checkbox', checked: st.mirror }); c.addEventListener('change', () => { st.mirror = c.checked; draw(); }); return el('label.chk', {}, c, 'mirror about y = 0'); })(),
    el('h4', { text: 'View' }),
    segEl([['top', 'Top'], ['iso', 'Iso'], ['front', 'Front']], st.view, (v) => { st.view = v; v3.setView(v); }),
    statsBox,
    ...(ctrlNames.length ? [el('div.desc', { text: `controls: ${ctrlNames.map(([n, m]) => `${n} (${m.filter(Boolean).length} boxes)`).join(', ')}` })] : []),
    el('div.desc', { text: 'Diverging scale is symmetric about 0 (blue < 0 < red), limits = max |value| of the current selection.' }));
  draw();
}

// ---------------------------------------------------------------- flutter
async function renderFlutter(run, main, ctrls) {
  const files = flutterFiles(run);
  loading(main, 'flutter results');
  const datas = await Promise.all(files.map((f) => getJson(run, f).catch(() => null)));
  clear(main);
  const items = files.map((f, i) => ({ file: f, d: datas[i] })).filter((x) => x.d);
  if (!items.length) { main.append(el('div.card', { text: 'no readable flutter_*.txt' })); return; }
  const label = (it) => it.d.method_title || it.d.method || it.file.replace(/^flutter_|\.txt$/g, '');
  let cur = Math.max(0, items.findIndex((x) => x.file === R.sub.flutter));
  const st = R.sub.flutterOpts || (R.sub.flutterOpts = { gmin: '', gmax: '' });

  const tiles = el('div.stats');
  const plots = el('div.plots');
  main.append(tiles, plots);
  const bg = plotBox(plots, 'V–g  (damping)');
  const bf = plotBox(plots, 'V–f  (frequency)');
  const pg = newPlot(bg, { height: 270, xlabel: 'V [m/s]', ylabel: 'g', padX: false });
  const pf = newPlot(bf, { height: 270, xlabel: 'V [m/s]', ylabel: 'f [Hz]', padX: false, legend: false });
  linkPlots([pg, pf]);
  // all-methods table
  const tcard = el('div.card', {}, el('h4', { text: 'Critical points (all methods)' }));
  const tbl = el('table.grid');
  tbl.append(el('tr', {}, ['method', 'kind', 'branch', 'V [m/s]', 'f [Hz]', 'k', 'ρ [kg/m³]'].map((h) => el('th', { text: h }))));
  for (const it of items) {
    for (const p of it.d.points || []) {
      tbl.append(el('tr', {}, el('td', { text: label(it) }), el('td', { text: p.kind }), el('td', { text: (it.d.names || [])[p.branch] || `#${p.branch + 1}` }),
        el('td', { text: fmt(p.velocity, 6) }), el('td', { text: fmt(p.freq_hz, 5) }), el('td', { text: fmt(p.k, 4) }), el('td', { text: fmt(p.density, 4) })));
    }
    if (!(it.d.points || []).length) tbl.append(el('tr', {}, el('td', { text: label(it) }), el('td', { text: 'no crossing in range', colspan: 6 })));
  }
  tcard.append(tbl);
  main.append(tcard);

  const draw = () => {
    const d = items[cur].d;
    R.sub.flutter = items[cur].file;
    const names = d.names || d.V[0].map((_, i) => `branch ${i + 1}`);
    const col = (A, j) => A.map((r) => r[j]);
    const nm = names.length;
    const sg = [], sf = [];
    for (let j = 0; j < nm; j++) {
      sg.push({ name: names[j], x: col(d.V, j), y: col(d.g, j), color: seriesColor(j) });
      sf.push({ name: names[j], x: col(d.V, j), y: col(d.freq, j), color: seriesColor(j) });
    }
    const pts = d.points || [];
    const markCol = (p) => (p.kind === 'flutter' ? '#b84d48' : p.kind === 'divergence' ? '#b76b2a' : '#4a5753');
    const short = (p) => (p.kind === 'flutter' ? 'F' : p.kind === 'divergence' ? 'D' : 'H');
    pg.setSeries(sg).setMarks(pts.map((p) => ({ x: p.velocity, y: d.g_struct || 0, label: `${short(p)} ${fmt(p.velocity, 4)}`, color: markCol(p) })));
    const hl = [{ y: 0, label: 'g = 0', fit: true }];
    if (d.g_struct) hl.push({ y: d.g_struct, label: `g_struct ${fmt(d.g_struct, 3)}`, dash: '6 3', color: '#b76b2a' });
    pg.setHLines(hl);
    const lim = {};
    if (st.gmin !== '' && Number.isFinite(Number(st.gmin))) lim.ymin = Number(st.gmin);
    if (st.gmax !== '' && Number.isFinite(Number(st.gmax))) lim.ymax = Number(st.gmax);
    pg.setOptions({ limits: lim }).render();
    pf.setSeries(sf).setMarks(pts.map((p) => ({ x: p.velocity, y: p.freq_hz, label: `${short(p)} ${fmt(p.freq_hz, 4)} Hz`, color: markCol(p) }))).render();
    // keep legend toggles in sync: hiding in V-g hides in V-f
    pg.setOptions({ onToggle: (s) => { const t = pf.series.find((x) => x.name === s.name); if (t) { t.hidden = s.hidden; pf.render(); } } });
    clear(tiles);
    const fl = pts.filter((p) => p.kind === 'flutter').sort((a, b) => a.velocity - b.velocity)[0];
    const dv = pts.filter((p) => p.kind === 'divergence').sort((a, b) => a.velocity - b.velocity)[0];
    tiles.append(stat(`flutter speed · ${label(items[cur])}`, fl ? fmt(fl.velocity, 5) : 'none', fl ? 'm/s' : '', true));
    if (fl) tiles.append(stat('flutter frequency', fmt(fl.freq_hz, 4), 'Hz'), stat('flutter k', fmt(fl.k, 4)), stat('branch', names[fl.branch] || String(fl.branch + 1)));
    if (dv) tiles.append(stat('divergence', fmt(dv.velocity, 5), 'm/s'));
    tiles.append(stat('density', fmt(d.density, 4), 'kg/m³'), stat('V points', String(d.V.length)));
  };
  if (items.length > 1) {
    const mlist = el('div.list');
    items.forEach((it, i) => mlist.append(el('div.item' + (i === cur ? '.active' : ''), {
      title: it.file,
      onclick: () => { cur = i; [...mlist.children].forEach((c, k) => c.classList.toggle('active', k === i)); draw(); },
    }, label(it), el('div.sub', { text: `${it.file} · ${(it.d.points || []).filter((p) => p.kind === 'flutter').length} flutter pt(s)` }))));
    ctrls.append(el('h4', { text: 'Method' }), mlist);
  } else ctrls.append(el('div.desc', { text: `method: ${label(items[0])}` }));
  const gmin = el('input', { value: st.gmin, placeholder: 'auto' }), gmax = el('input', { value: st.gmax, placeholder: 'auto' });
  gmin.addEventListener('change', () => { st.gmin = gmin.value.trim(); draw(); });
  gmax.addEventListener('change', () => { st.gmax = gmax.value.trim(); draw(); });
  ctrls.append(el('h4', { text: 'V–g axis' }), ctrlRow('g min', gmin), ctrlRow('g max', gmax),
    el('div.desc', { text: 'Click legend entries to hide branches. Hover for a linked crosshair. F = flutter, D = divergence, H = hump recovery.' }),
    el('div.desc.mono', { text: items.map((it) => it.file).join('\n') }));
  draw();
}

// ---------------------------------------------------------------- grouped line charts (time / curves)
function groupOf(name) { return name.includes(':') ? name.split(':')[0] : name; }

function groupedCharts(main, ctrls, { x, series, xlabel, ylabel = '', logx = false, logy = false, key, height = 190 }) {
  const names = Object.keys(series);
  const st = R.sub[key] || (R.sub[key] = { hidden: {} });
  const groups = [];
  for (const n of names) { const g = groupOf(n); let G = groups.find((q) => q.name === g); if (!G) groups.push(G = { name: g, names: [] }); G.names.push(n); }
  const plotsBox = el('div.plots'); main.append(plotsBox);
  const colorOf = (n) => { const G = groups.find((q) => q.names.includes(n)); return seriesColor(G.names.indexOf(n)); };
  let plots = [];
  const draw = () => {
    for (const p of plots) p.destroy();
    clear(plotsBox); plots = [];
    for (const G of groups) {
      const vis = G.names.filter((n) => !st.hidden[n]);
      if (!vis.length) continue;
      const title = G.names.length > 1 ? `${G.name}  (${G.names.length})` : G.names[0];
      const box = plotBox(plotsBox, title);
      const p = new Plot(box, { height, xlabel, ylabel: G.names.length > 1 ? G.name : ylabel || G.names[0], logx, logy, legend: G.names.length > 1, padX: false });
      p.setSeries(vis.map((n) => ({ name: G.names.length > 1 ? n.split(':').slice(1).join(':') || n : n, x, y: series[n], color: colorOf(n) }))).render();
      plots.push(p);
    }
    linkPlots(plots);
    if (!plots.length) plotsBox.append(el('div.empty-note', { text: 'all series hidden' }));
  };
  onCleanup(() => plots.forEach((p) => p.destroy()));
  const tog = el('div.series-toggles');
  for (const G of groups) {
    if (G.names.length > 1) tog.append(el('div.small.dim', { style: { marginTop: '4px' }, text: G.name }));
    for (const n of G.names) {
      const c = el('input', { type: 'checkbox', checked: !st.hidden[n] });
      c.addEventListener('change', () => { st.hidden[n] = !c.checked; draw(); });
      tog.append(el('label', {}, c, el('span.sw', { style: { background: colorOf(n) } }), G.names.length > 1 ? n.split(':').slice(1).join(':') || n : n));
    }
  }
  const allBtn = el('button.btn.small', { text: 'all', onclick: () => { st.hidden = {}; for (const c of $$('input', tog)) c.checked = true; draw(); } });
  const noneBtn = el('button.btn.small', { text: 'none', onclick: () => { for (const n of names) st.hidden[n] = true; for (const c of $$('input', tog)) c.checked = false; draw(); } });
  ctrls.append(el('h4', {}, 'Series ', allBtn, ' ', noneBtn), tog);
  draw();
}

async function renderTime(run, main, ctrls) {
  const files = byPrefix(run, /^timehist_.*\.txt$/i);
  let cur = Math.max(0, files.indexOf(R.sub.time));
  const show = async () => {
    runCleanup(); clear(main); clear(ctrls);
    R.sub.time = files[cur];
    if (files.length > 1) ctrls.append(ctrlRow('file', selectEl(files.map((f, i) => [i, f]), cur, (v) => { cur = Number(v); show(); })));
    loading(main, files[cur]);
    const d = await getJson(run, files[cur]);
    clear(main);
    const meta = Object.entries(d.meta || {});
    main.append(el('div.stats', {}, stat('file', files[cur].replace(/^timehist_|\.txt$/g, '')), stat('samples', String(d.t.length)),
      stat('duration', fmt(d.t[d.t.length - 1] - d.t[0], 4), 's'), ...meta.map(([k, v]) => stat(k, typeof v === 'number' ? fmt(v, 5) : String(v)))));
    groupedCharts(main, ctrls, { x: d.t, series: d.series, xlabel: 't [s]', key: 'time:' + files[cur] });
  };
  await show();
}

async function renderCurves(run, main, ctrls) {
  const files = byPrefix(run, /^curves_.*\.txt$/i);
  let cur = Math.max(0, files.indexOf(R.sub.curves));
  const show = async () => {
    runCleanup(); clear(main); clear(ctrls);
    R.sub.curves = files[cur];
    if (files.length > 1) ctrls.append(ctrlRow('file', selectEl(files.map((f, i) => [i, f]), cur, (v) => { cur = Number(v); show(); })));
    loading(main, files[cur]);
    const d = await getJson(run, files[cur]);
    clear(main);
    const key = 'curves:' + files[cur];
    const o = R.sub[key + ':axes'] || (R.sub[key + ':axes'] = {
      logy: !!d.logy,
      logx: d.x.length > 2 && d.x[0] > 0 && d.x[d.x.length - 1] / d.x[0] > 300,
    });
    const lx = el('input', { type: 'checkbox', checked: o.logx }), ly = el('input', { type: 'checkbox', checked: o.logy });
    lx.addEventListener('change', () => { o.logx = lx.checked; show(); });
    ly.addEventListener('change', () => { o.logy = ly.checked; show(); });
    ctrls.append(el('div.row.gap', {}, el('label.chk', {}, lx, 'log x'), el('label.chk', {}, ly, 'log y')));
    main.append(el('div.stats', {}, stat('file', files[cur].replace(/^curves_|\.txt$/g, '')), stat('points', String(d.x.length)),
      stat('curves', String(Object.keys(d.curves).length))));
    groupedCharts(main, ctrls, { x: d.x, series: d.curves, xlabel: d.xlabel || 'x', ylabel: d.ylabel || '', logx: o.logx, logy: o.logy, key, height: 210 });
  };
  await show();
}

// ---------------------------------------------------------------- TSD
async function renderTsd(run, main, ctrls) {
  loading(main, 'TSD');
  const d = await getJson(run, 'tsd_steady.txt');
  clear(main);
  const ny = d.y.length;
  const st = R.sub.tsd || (R.sub.tsd = { norm: true, every: Math.max(1, Math.round(ny / 6)) });
  main.append(el('div.stats', {}, stat('Mach', fmt(d.mach, 4), '', true), stat('Cp*', fmt(d.cp_star, 4)),
    stat('converged', d.converged ? 'yes' : 'no'), stat('iterations', String((d.residuals || []).length)), stat('stations', String(ny))));
  const row1 = el('div.plots', { style: { gridTemplateColumns: 'minmax(0,1.3fr) minmax(0,1fr)' } });
  main.append(row1);
  const b1 = plotBox(row1, 'Chordwise Cp (upper / mean surface)');
  const b2 = el('div.plotbox', {}, el('div.ptitle', { text: 'Cp on the planform' }));
  row1.append(b2);
  const vbox = el('div', { style: { height: '300px', position: 'relative' } }); b2.append(vbox);
  const v3 = new View3D(vbox, { view: 'top', triad: false, onHover: (f) => f.data || '' });
  onCleanup(() => v3.destroy());
  const cb = colorbar(); vbox.append(cb.el);
  const p1 = newPlot(b1, { height: 300, xlabel: st.norm ? 'x/c' : 'x [m]', ylabel: 'Cp', invertY: true, padX: false });
  const draw = () => {
    const ser = [];
    for (let j = 0; j < ny; j += st.every) {
      const xs = d.x[j];
      const x0 = Math.min(...xs), x1 = Math.max(...xs);
      ser.push({ name: `y = ${fmt(d.y[j], 3)}`, x: st.norm ? xs.map((x) => (x - x0) / ((x1 - x0) || 1)) : xs, y: d.cp[j], color: seriesColor(ser.length) });
    }
    p1.setOptions({ xlabel: st.norm ? 'x/c (station extent)' : 'x [m]' }).setSeries(ser)
      .setHLines([{ y: d.cp_star, label: `Cp* ${fmt(d.cp_star, 3)}`, dash: '6 3', color: '#b84d48' }]).render();
  };
  // planform map: cells between neighbouring grid points
  let lim = 0; for (const r of d.cp) for (const v of r) lim = Math.max(lim, Math.abs(v));
  lim = lim || 1;
  const faces = [];
  for (let j = 0; j + 1 < ny; j++) {
    const n = Math.min(d.x[j].length, d.x[j + 1].length);
    for (let i = 0; i + 1 < n; i++) {
      const v = 0.25 * (d.cp[j][i] + d.cp[j][i + 1] + d.cp[j + 1][i] + d.cp[j + 1][i + 1]);
      faces.push({ pts: [[d.x[j][i], d.y[j], 0], [d.x[j][i + 1], d.y[j], 0], [d.x[j + 1][i + 1], d.y[j + 1], 0], [d.x[j + 1][i], d.y[j + 1], 0]],
        fill: diverging(-v / lim), stroke: null, data: `y ${fmt(d.y[j], 3)}  x ${fmt(d.x[j][i], 3)}  Cp ${fmt(v, 4)}` });
    }
  }
  v3.setScene({ faces }).fit(); v3.render();
  cb.set(lim, -lim, 'Cp (red = suction)', 'div', (x) => fmt(x, 3));
  if (d.residuals && d.residuals.length > 1) {
    const b3 = plotBox(main, 'Convergence');
    newPlot(b3, { height: 180, xlabel: 'iteration', ylabel: 'residual', logy: true, legend: false })
      .setSeries([{ name: 'residual', x: d.residuals.map((_, i) => i + 1), y: d.residuals }]).render();
  }
  const norm = el('input', { type: 'checkbox', checked: st.norm });
  norm.addEventListener('change', () => { st.norm = norm.checked; draw(); });
  ctrls.append(el('label.chk', {}, norm, 'normalise x by station extent'),
    ctrlRow('every n-th station', selectEl([1, 2, 3, 4, 5, 6, 8, 10].filter((k) => k <= ny).map((k) => [k, String(k)]), st.every, (v) => { st.every = Number(v); draw(); })),
    el('div.desc', { text: 'Cp axis is inverted (suction up). Dashed line: critical Cp* (sonic).' }));
  draw();
}

// ---------------------------------------------------------------- figures / files
function showLightbox(path, cap) {
  const lb = $('#lightbox');
  $('img', lb).src = fileUrl(path);
  $('.cap', lb).textContent = cap || path;
  lb.hidden = false;
}

async function renderFigures(run, main, ctrls) {
  const g = el('div.gallery');
  const titles = Object.fromEntries((run.artifacts || []).filter((a) => a.kind === 'figure').map((a) => [a.file, a.title]));
  for (const f of run.figures) {
    const cap = titles[f] || f.replace(/^figures\//, '');
    g.append(el('figure', { onclick: () => showLightbox(fpath(run, f), cap) },
      el('img', { src: fileUrl(fpath(run, f)) + `&t=${Math.round(run.mtime)}`, alt: cap, loading: 'lazy' }), el('figcaption', { text: cap })));
  }
  main.append(g);
  ctrls.append(el('div.desc', { text: `${run.figures.length} PNG figure(s) written by the run (matplotlib). Click to enlarge.` }));
}

async function renderFiles(run, main, ctrls) {
  const files = textFiles(run);
  let cur = files.indexOf(R.sub.file);
  if (cur < 0) cur = Math.max(0, files.indexOf('results.txt'));
  const pre = el('pre.text', { style: { height: '100%' } });
  main.append(pre);
  const list = el('div.list');
  const show = async (i) => {
    cur = i; R.sub.file = files[i];
    for (const [k, it] of [...list.children].entries()) it.classList.toggle('active', k === i);
    pre.textContent = 'loading…';
    try {
      const t = await api.text(fpath(run, files[i]));
      const MAX = 400000;   // large result tables (cp.txt, modes.txt) - show the head, open raw for the rest
      pre.textContent = t.length > MAX ? t.slice(0, MAX) + `\n\n% ... ${((t.length - MAX) / 1e6).toFixed(1)} MB more - open the raw file` : t;
    } catch (e) { pre.textContent = e.message; }
  };
  files.forEach((f, i) => list.append(el('div.item', { text: f, onclick: () => show(i) })));
  const raw = el('a.small', { href: fileUrl(fpath(run, files[cur] || 'results.txt')), target: '_blank', text: 'open raw file', style: { display: 'block', marginTop: '8px' } });
  list.addEventListener('click', () => { raw.href = fileUrl(fpath(run, files[cur])); });
  ctrls.append(list, raw);
  await show(cur);
}

// ---------------------------------------------------------------- init
export function initResultsTab(h = {}) {
  hooks = h;
  $('#btnReloadRuns').addEventListener('click', () => loadRuns());
  $('#runFilter').addEventListener('input', renderRunList);
  $('#lightbox').addEventListener('click', () => { $('#lightbox').hidden = true; });
  window.addEventListener('keydown', (e) => { if (e.key === 'Escape') $('#lightbox').hidden = true; });
}

export const currentRun = () => R.run;
export const currentView = () => R.view;
