// Model workspace: config library, model tree, planform editor, property panel,
// analysis form, two-way config text, validation and runs.

import { api } from './api.js';
import { el, $, $$, clear, debounce, toast, fmt, SURF_COLORS } from './util.js';
import { parseConfig, generateConfig, emptyModel, getPlate, setPlate, pfDerived, controlBoxes, num, fnum, isIndexedKey } from './config.js';
import { PlanformEditor } from './editor.js';

const S = {
  models: null, presets: {}, configs: [],
  m: emptyModel(), source: null, sel: null,
  linkPlate: false, openGroups: new Set(['common']),
  job: null, jobTimer: null, lastResult: null, lastDiag: null,
};
let editor = null;
let hooks = {};

const PLATE_KEYS = ['plate_root_le', 'plate_root_chord', 'plate_tip_le', 'plate_tip_chord', 'plate_nchord', 'plate_nspan'];
const NEW_TEMPLATE = `% FLAAP config - new model (FLAAP Studio)
model        flutter
mode         pk
title        New plate wing
structure    plate
plate_root_le     0.0, 0.0, 0.0
plate_root_chord  0.5
plate_tip_le      0.2, 1.2, 0.0
plate_tip_chord   0.3
thickness    0.003
material     aluminum
clamp        y, 0.0
nmodes       4
mach         0.3
density      1.225
velocities   10, 150, 36
surface1     wing, 0, 0, 0, 0.5, 0.2, 1.2, 0, 0.3, 8, 12
`;

// ---------------------------------------------------------------- text <-> model
const scheduleValidate = debounce(() => validate(false), 700);

function modelChanged({ text = true, tree = true, props = true, form = false, validate: v = true } = {}) {
  if (text) {
    $('#cfgText').value = generateConfig(S.m);
  }
  if (tree) renderTree();
  if (props) renderProps();
  if (form) renderForm();
  editor.render();
  if (v) scheduleValidate();
}

function loadText(text, source = null, { fit = true } = {}) {
  const { model } = parseConfig(text);
  S.m = model;
  S.source = source;
  $('#cfgText').value = text;
  $('#cfgSource').textContent = source ? `source: ${source}` : 'unsaved model';
  S.sel = null;
  editor.select(null);
  editor.setModel(S.m);
  syncRouteSelectors();
  renderTree(); renderProps(); renderForm();
  if (fit) editor.fit();
  const nm = source ? source.split('/').pop().replace(/\.txt$/, '').replace(/^config_/, '') : '';
  $('#runName').value = nm || slug(S.m.values.title) || 'model';
  renderConfigList();
  validate(false);
}

const slug = (s) => String(s || '').toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_|_$/g, '').slice(0, 40);

const onTextInput = debounce(() => {
  const { model } = parseConfig($('#cfgText').value);
  S.m = model;
  editor.setModel(S.m);
  if (S.sel && !selValid(S.sel)) { S.sel = null; editor.select(null); }
  syncRouteSelectors();
  renderTree(); renderProps(); renderForm();
  validate(false);
}, 400);

function selValid(sel) {
  if (!sel) return true;
  if (sel.kind === 'surface') return !!S.m.surfaces[sel.i];
  if (sel.kind === 'control') return !!S.m.controls[sel.i];
  if (sel.kind === 'station') return !!S.m.stations[sel.i];
  if (sel.kind === 'plate') return !!editor.plate();
  if (sel.kind === 'ghost') return !!editor.ghosts()[sel.i];
  return true;
}

// ---------------------------------------------------------------- selection
function select(sel) {
  S.sel = sel;
  editor.select(sel);
  renderTree(); renderProps();
}

// ---------------------------------------------------------------- config library
async function loadConfigList() {
  try {
    const r = await api.configs();
    S.configs = r.configs;
    renderConfigList();
  } catch (e) { toast(`configs: ${e.message}`, true); }
}

function renderConfigList() {
  const box = clear($('#cfgList'));
  const f = $('#cfgFilter').value.trim().toLowerCase();
  let lastGroup = null;
  for (const c of S.configs) {
    if (f && !c.path.toLowerCase().includes(f)) continue;
    const grp = c.group === 'studio_runs' ? 'studio runs' : c.folder || c.group;
    if (grp !== lastGroup) { box.append(el('div.grp', { text: grp })); lastGroup = grp; }
    box.append(el('div.item' + (c.path === S.source ? '.active' : ''), {
      title: c.path,
      onclick: async () => {
        try { const r = await api.config(c.path); loadText(r.text, r.path); } catch (e) { toast(e.message, true); }
      },
    }, c.name.replace(/\.txt$/, '')));
  }
  if (!box.childElementCount) box.append(el('div.empty-note', { style: { padding: '6px 8px' }, text: 'no configs' }));
}

// ---------------------------------------------------------------- tree
function renderTree() {
  const t = clear($('#tree'));
  const node = (label, sel, opts = {}) => {
    const active = S.sel && sel && S.sel.kind === sel.kind && S.sel.i === sel.i;
    return el('div.node' + (active ? '.active' : ''), { onclick: () => select(sel), title: opts.title },
      opts.color ? el('span.sw', { style: { background: opts.color, border: opts.dashed ? '1px dashed #4a5753' : '' } }) : null,
      el('span', { text: label }), opts.chip ? el('span.chip' + (opts.chipCls || ''), { text: opts.chip }) : null);
  };
  const v = S.m.values;
  t.append(el('div.cat', { text: 'Analysis' }));
  t.append(el('div.kids', {}, node(`${v.model || '?'} / ${v.mode || '(default)'}`, null, { chip: v.title ? '' : null })));
  t.append(el('div.cat', { text: 'Structure' }));
  const st = editor.structure();
  const sk = el('div.kids');
  if (st === 'plate') sk.append(node('plate (FE planform)', { kind: 'plate' }, { color: '#3f6f91' }));
  else sk.append(node(st + (v.structure ? '' : ' (default)'), null, { chip: v.structure_file ? 'file' : null }));
  t.append(sk);
  t.append(el('div.cat', { text: 'Lifting surfaces' }));
  const sf = el('div.kids');
  S.m.surfaces.forEach((s, i) => {
    sf.append(node(s.name, { kind: 'surface', i }, { color: SURF_COLORS[i % SURF_COLORS.length], chip: `${s.nchord}×${s.nspan}` }));
    const kids = el('div.kids');
    S.m.controls.forEach((c, k) => { if (c.surface === s.name) kids.append(node(c.name, { kind: 'control', i: k }, { color: '#b76b2a', chip: 'control' })); });
    if (kids.childElementCount) sf.append(kids);
  });
  const orphan = S.m.controls.map((c, k) => [c, k]).filter(([c]) => !S.m.surfaces.some((s) => s.name === c.surface));
  for (const [c, k] of orphan) sf.append(node(`${c.name} → ${c.surface}?`, { kind: 'control', i: k }, { color: '#b84d48', chip: 'no surface', chipCls: '.error' }));
  editor.ghosts().forEach((g, i) => sf.append(node(`${g.name} (preset)`, { kind: 'ghost', i }, { color: '#e9e6dc', dashed: true, chip: `${g.nchord}×${g.nspan}` })));
  const imp = editor.aeroSurfaces().find((s) => s.implicit);
  if (imp) sf.append(node('plate planform (implicit)', { kind: 'plate' }, { color: '#e9e6dc', dashed: true, chip: `${imp.nchord}×${imp.nspan}` }));
  if (!sf.childElementCount) sf.append(el('div.empty', { text: 'none - draw one or add +Surface' }));
  t.append(sf);
  t.append(el('div.cat', { text: 'Monitor stations' }));
  const ss = el('div.kids');
  S.m.stations.forEach((s, i) => ss.append(node(s.name, { kind: 'station', i }, { color: '#6a4ce2' })));
  if (!ss.childElementCount) ss.append(el('div.empty', { text: 'none' }));
  t.append(ss);
}

// ---------------------------------------------------------------- property panel
function numInput(value, onval, { step = 'any', int = false, min } = {}) {
  const inp = el('input', { type: 'number', step, value: int ? String(value) : fnum(value) });
  if (min !== undefined) inp.min = min;
  inp.addEventListener('input', () => {
    const v = Number(inp.value);
    const ok = inp.value.trim() !== '' && Number.isFinite(v) && (!int || Number.isInteger(v)) && (min === undefined || v >= min);
    inp.classList.toggle('invalid', !ok);
    if (ok) onval(v);
  });
  return inp;
}
function vec3(arr, onval) {
  return el('div.vec3', {}, arr.map((v, k) => numInput(v, (x) => { arr[k] = x; onval(); })));
}
const row = (label, ctrl, help) => el('div.form-row', {}, el('label', { text: label, title: label }), ctrl, help ? el('div.help', { text: help }) : null);

function derivedBox(s) {
  const d = pfDerived(s);
  return el('div.derived', {},
    'span', el('b', { text: `${fmt(d.span, 4)} m` }), 'area', el('b', { text: `${fmt(d.area, 4)} m²` }),
    'aspect ratio', el('b', { text: fmt(d.ar, 3) }), 'taper', el('b', { text: fmt(d.taper, 3) }),
    'sweep LE', el('b', { text: `${fmt(d.sweepLE, 3)}°` }), 'sweep c/4', el('b', { text: `${fmt(d.sweepQC, 3)}°` }),
    'dihedral', el('b', { text: `${fmt(d.dihedral, 3)}°` }), 'MAC', el('b', { text: `${fmt(d.mac, 4)} m` }),
    'boxes', el('b', { text: String(d.boxes) }));
}

const liveChange = (kind) => {
  if (kind === 'surface' && S.linkPlate && S.sel && S.sel.i === 0) setPlate(S.m, S.m.surfaces[0]);
  if (kind === 'plate' && S.linkPlate && S.m.surfaces[0]) Object.assign(S.m.surfaces[0], pick(getPlate(S.m)));
  modelChanged({ props: false, tree: kind !== 'station' && kind !== 'plate' });
  refreshDerived();
};
const pick = (p) => ({ rle: [...p.rle], rc: p.rc, tle: [...p.tle], tc: p.tc });

function refreshDerived() {
  const box = $('#props .derived');
  if (!box || !S.sel) return;
  const s = S.sel.kind === 'surface' ? S.m.surfaces[S.sel.i] : S.sel.kind === 'plate' ? getPlate(S.m) : null;
  if (s) box.replaceWith(derivedBox(s));
}

function renderProps() {
  const p = clear($('#props'));
  const sel = S.sel;
  if (!sel) return renderOverview(p);
  if (sel.kind === 'surface') {
    const s = S.m.surfaces[sel.i];
    if (!s) return renderOverview(p);
    const nameIn = el('input', { value: s.name });
    nameIn.addEventListener('change', () => {
      const nv = nameIn.value.trim().replace(/[\s,]+/g, '_');
      if (!nv) return;
      for (const c of S.m.controls) if (c.surface === s.name) c.surface = nv;
      s.name = nv; modelChanged();
    });
    p.append(el('h4', { text: `Surface ${sel.i + 1}` }),
      row('name', nameIn),
      row('root LE x,y,z', vec3(s.rle, () => liveChange('surface'))),
      row('root chord', numInput(s.rc, (v) => { s.rc = v; liveChange('surface'); }, { min: 0 })),
      row('tip LE x,y,z', vec3(s.tle, () => liveChange('surface'))),
      row('tip chord', numInput(s.tc, (v) => { s.tc = v; liveChange('surface'); }, { min: 0 })),
      row('nchord', numInput(s.nchord, (v) => { s.nchord = v; liveChange('surface'); }, { int: true, step: 1, min: 1 })),
      row('nspan', numInput(s.nspan, (v) => { s.nspan = v; liveChange('surface'); }, { int: true, step: 1, min: 1 })),
      derivedBox(s),
      el('div.row.gap.wrap', {},
        el('button.btn.small', { onclick: () => addControl(s.name) }, '+ Control on this'),
        el('button.btn.small', { onclick: () => duplicateSurface(sel.i) }, 'Duplicate'),
        el('button.btn.small', { onclick: () => { s.tle = [s.tle[0], -s.tle[1], s.tle[2]]; s.rle = [s.rle[0], -s.rle[1], s.rle[2]]; modelChanged(); } }, 'Flip y'),
        editor.structure() === 'plate' && sel.i === 0
          ? el('button.btn.small', { onclick: () => { setPlate(S.m, s); modelChanged({ form: true }); } }, 'Plate ← this') : null,
        el('button.btn.small.danger', { onclick: () => deleteSel() }, 'Delete')));
    return;
  }
  if (sel.kind === 'control') {
    const c = S.m.controls[sel.i];
    if (!c) return renderOverview(p);
    const sfSel = el('select', {}, S.m.surfaces.map((s) => el('option', { value: s.name, text: s.name, selected: s.name === c.surface })));
    if (!S.m.surfaces.some((s) => s.name === c.surface)) sfSel.prepend(el('option', { value: c.surface, text: `${c.surface} (missing)`, selected: true }));
    sfSel.addEventListener('change', () => { c.surface = sfSel.value; modelChanged(); });
    const nameIn = el('input', { value: c.name });
    nameIn.addEventListener('change', () => { c.name = nameIn.value.trim().replace(/[\s,]+/g, '_') || c.name; modelChanged(); });
    const sf = S.m.surfaces.find((s) => s.name === c.surface);
    const nb = sf ? controlBoxes(sf, c).length : 0;
    p.append(el('h4', { text: `Control ${sel.i + 1}` }),
      row('name', nameIn), row('surface', sfSel),
      row('hinge x/c', numInput(c.hinge, (v) => { c.hinge = v; liveChange('control'); }, { step: 0.01, min: 0 }), 'boxes with centre aft of the hinge'),
      row('eta0', numInput(c.eta0, (v) => { c.eta0 = v; liveChange('control'); }, { step: 0.01, min: 0 })),
      row('eta1', numInput(c.eta1, (v) => { c.eta1 = v; liveChange('control'); }, { step: 0.01, min: 0 }), 'span fractions of the surface'),
      el('div.derived', {}, 'control boxes', el('b', { text: String(nb) })),
      el('div.row.gap', {}, el('button.btn.small.danger', { onclick: () => deleteSel() }, 'Delete')));
    return;
  }
  if (sel.kind === 'station') {
    const s = S.m.stations[sel.i];
    if (!s) return renderOverview(p);
    const nameIn = el('input', { value: s.name });
    nameIn.addEventListener('change', () => { s.name = nameIn.value.trim().replace(/[\s,]+/g, '_') || s.name; modelChanged(); });
    p.append(el('h4', { text: `Station ${sel.i + 1}` }), row('name', nameIn),
      row('point x,y,z', vec3(s.p, () => liveChange('station'))),
      row('axis x,y,z', vec3(s.a, () => liveChange('station')), 'cut-plane normal (loads about this axis)'),
      el('div.row.gap', {}, el('button.btn.small.danger', { onclick: () => deleteSel() }, 'Delete')));
    return;
  }
  if (sel.kind === 'plate') {
    const pf = getPlate(S.m);
    const upd = () => { setPlate(S.m, pf); liveChange('plate'); };
    const setv = (k, v) => { if (v === '' || v === undefined) delete S.m.values[k]; else S.m.values[k] = String(v); modelChanged({ props: false }); };
    const mats = (S.models && S.models.choices.material) || ['aluminum'];
    const matSel = el('select', {}, el('option', { value: '', text: '(default)' }),
      mats.map((m) => el('option', { value: m, text: m, selected: S.m.values.material === m })));
    matSel.addEventListener('change', () => setv('material', matSel.value));
    const thick = el('input', { value: S.m.values.thickness ?? '', placeholder: '0.002' });
    thick.addEventListener('change', () => setv('thickness', thick.value.trim()));
    const clamp = el('input', { value: S.m.values.clamp ?? '', placeholder: 'y, 0.0' });
    clamp.addEventListener('change', () => setv('clamp', clamp.value.trim()));
    const nc = el('input', { type: 'number', step: 1, min: 1, value: S.m.values.plate_nchord ?? '', placeholder: '8' });
    nc.addEventListener('change', () => setv('plate_nchord', nc.value.trim()));
    const ns = el('input', { type: 'number', step: 1, min: 1, value: S.m.values.plate_nspan ?? '', placeholder: '12' });
    ns.addEventListener('change', () => setv('plate_nspan', ns.value.trim()));
    const link = el('input', { type: 'checkbox', checked: S.linkPlate });
    link.addEventListener('change', () => {
      S.linkPlate = link.checked;
      if (S.linkPlate && S.m.surfaces[0]) { Object.assign(S.m.surfaces[0], pick(getPlate(S.m))); modelChanged(); }
    });
    p.append(el('h4', { text: 'Plate structure (FE)' }),
      row('root LE x,y,z', vec3(pf.rle, upd)),
      row('root chord', numInput(pf.rc, (v) => { pf.rc = v; upd(); }, { min: 0 })),
      row('tip LE x,y,z', vec3(pf.tle, upd)),
      row('tip chord', numInput(pf.tc, (v) => { pf.tc = v; upd(); }, { min: 0 })),
      row('plate_nchord', nc), row('plate_nspan', ns),
      row('thickness', thick, 'm'), row('material', matSel), row('clamp', clamp, "'y, 0.0' or none"),
      derivedBox(pf),
      el('label.chk', { style: { margin: '4px 0' } }, link, 'surface1 follows plate planform'),
      el('div.row.gap.wrap', {},
        S.m.surfaces[0] ? el('button.btn.small', { onclick: () => { Object.assign(S.m.surfaces[0], pick(getPlate(S.m))); modelChanged(); } }, 'surface1 ← plate') : null,
        S.m.surfaces[0] ? el('button.btn.small', { onclick: () => { setPlate(S.m, S.m.surfaces[0]); modelChanged({ form: true }); } }, 'plate ← surface1') : null,
        !S.m.surfaces.length ? el('button.btn.small', { onclick: () => {
          S.m.surfaces.push({ name: 'wing', ...pick(getPlate(S.m)), nchord: Math.round(num(S.m.values.nchord, 8)), nspan: Math.round(num(S.m.values.nspan, 12)) });
          S.linkPlate = true; select({ kind: 'surface', i: 0 }); modelChanged();
        } }, 'Make aero surface1 from plate') : null),
      el('div.help.small.dim', { style: { marginTop: '6px' }, text: 'Without surface1 the plate planform is also the DLM surface (nchord / nspan from the aero keys).' }));
    return;
  }
  if (sel.kind === 'ghost') {
    const g = editor.ghosts();
    p.append(el('h4', { text: `Preset planform: ${editor.structure()}` }),
      el('div.desc', { text: 'This structure brings its own aerodynamic surface(s). They are shown read-only. Convert them to editable surfaceN keys to modify the planform (the preset structure stays).' }),
      el('table.kv', {}, g.map((s) => el('tr', {}, el('td', { text: s.name }), el('td.num', { text: `c ${fmt(s.rc, 4)} → ${fmt(s.tc, 4)}, span ${fmt(pfDerived(s).span, 4)}, ${s.nchord}×${s.nspan}` })))),
      el('div.row.gap', {}, el('button.btn.primary.small', { onclick: convertGhosts }, 'Convert to editable surfaces')));
  }
}

function renderOverview(p) {
  const aero = editor.aeroSurfaces();
  let area = 0, boxes = 0;
  for (const s of aero) { const d = pfDerived(s); area += d.area; boxes += d.boxes; }
  p.append(el('h4', { text: S.m.values.title || 'Model' }),
    el('div.derived', {}, 'structure', el('b', { text: editor.structure() || '–' }),
      'aero surfaces', el('b', { text: String(aero.length) }), 'controls', el('b', { text: String(S.m.controls.length) }),
      'stations', el('b', { text: String(S.m.stations.length) }), 'DLM boxes', el('b', { text: String(boxes) }),
      'planform area', el('b', { text: `${fmt(area, 4)} m²` })),
    el('div.desc', { html: 'Select an item in the tree or the canvas.<br><b>Draw surface</b>: four clicks (root LE, tip LE, tip TE, root TE).<br>Keys: <span class="mono">V</span> select, <span class="mono">S</span> draw, <span class="mono">M</span> station, <span class="mono">F</span> fit, <span class="mono">Del</span> delete, <span class="mono">Esc</span> cancel.' }));
}

// ---------------------------------------------------------------- add / delete
function uniqueName(base, list) {
  const names = new Set(list.map((x) => x.name));
  if (!names.has(base)) return base;
  for (let i = 2; ; i++) if (!names.has(`${base}${i}`)) return `${base}${i}`;
}
function addSurface(s) {
  const n = S.m.surfaces.length;
  if (!s) {
    const last = S.m.surfaces[n - 1];
    s = last ? { ...JSON.parse(JSON.stringify(last)), rle: [last.rle[0] + 2 * last.rc, last.rle[1], last.rle[2]], tle: [last.tle[0] + 2 * last.rc, last.tle[1], last.tle[2]] }
      : { rle: [0, 0, 0], rc: 1, tle: [0.3, 3, 0], tc: 0.6, nchord: 8, nspan: 10 };
  }
  s.name = s.name || uniqueName(n ? 'surf' : 'wing', S.m.surfaces);
  S.m.surfaces.push(s);
  select({ kind: 'surface', i: S.m.surfaces.length - 1 });
  modelChanged();
}
function duplicateSurface(i) {
  const s = JSON.parse(JSON.stringify(S.m.surfaces[i]));
  s.name = uniqueName(s.name, S.m.surfaces);
  s.rle[0] += s.rc * 1.5; s.tle[0] += s.rc * 1.5;
  addSurface(s);
}
function addControl(surfaceName) {
  const sf = surfaceName || (S.sel && S.sel.kind === 'surface' && S.m.surfaces[S.sel.i] ? S.m.surfaces[S.sel.i].name : S.m.surfaces[0] && S.m.surfaces[0].name);
  if (!sf) { toast('add a surface first (controls need an explicit surfaceN)', true); return; }
  S.m.controls.push({ name: uniqueName('flap', S.m.controls), surface: sf, hinge: 0.75, eta0: 0.5, eta1: 0.9 });
  select({ kind: 'control', i: S.m.controls.length - 1 });
  modelChanged();
}
function addStation(st) {
  const s = st || { p: [0, 0, 0], a: [0, 1, 0] };
  if (!st) {
    const sf = editor.aeroSurfaces()[0];
    if (sf) s.p = [sf.rle[0] + 0.4 * sf.rc, sf.rle[1], sf.rle[2]];
  }
  s.name = s.name || uniqueName('sta', S.m.stations);
  S.m.stations.push(s);
  select({ kind: 'station', i: S.m.stations.length - 1 });
  modelChanged();
}
function deleteSel() {
  const sel = S.sel; if (!sel) return;
  if (sel.kind === 'surface') {
    const s = S.m.surfaces[sel.i];
    S.m.surfaces.splice(sel.i, 1);
    const n0 = S.m.controls.length;
    S.m.controls = S.m.controls.filter((c) => c.surface !== s.name);
    if (n0 !== S.m.controls.length) toast(`removed ${n0 - S.m.controls.length} control(s) on ${s.name}`);
  } else if (sel.kind === 'control') S.m.controls.splice(sel.i, 1);
  else if (sel.kind === 'station') S.m.stations.splice(sel.i, 1);
  else return;
  select(null);
  modelChanged();
}
function convertGhosts() {
  const g = editor.ghosts();
  if (!g.length) return;
  for (const s of g) {
    S.m.surfaces.push({ name: s.name, rle: [...s.rle], rc: s.rc, tle: [...s.tle], tc: s.tc, nchord: s.nchord, nspan: s.nspan });
    for (const c of s.controls || []) S.m.controls.push({ name: c.name, surface: s.name, hinge: c.hinge, eta0: c.eta0, eta1: c.eta1 });
  }
  select({ kind: 'surface', i: 0 });
  modelChanged();
  toast(`${g.length} preset surface(s) converted to surfaceN keys`);
}

// ---------------------------------------------------------------- analysis form
function route() { return S.models && S.models.routes[String(S.m.values.model || '').toLowerCase()]; }

function syncRouteSelectors() {
  if (!S.models) return;
  const sm = $('#selModel'), sd = $('#selMode');
  const model = String(S.m.values.model || '').toLowerCase();
  clear(sm);
  if (!S.models.routes[model]) sm.append(el('option', { value: model, text: model ? `${model} (unknown)` : '(choose)' }));
  for (const k of S.models.order) sm.append(el('option', { value: k, text: k, selected: k === model }));
  const r = route();
  clear(sd);
  const mode = String(S.m.values.mode || '').toLowerCase();
  if (r) {
    if (!mode) sd.append(el('option', { value: '', text: '(default)', selected: true }));
    if (mode && !(mode in r.modes)) sd.append(el('option', { value: mode, text: `${mode} (unknown)`, selected: true }));
    for (const [k, d] of Object.entries(r.modes)) sd.append(el('option', { value: k, text: k, title: d, selected: k === mode }));
    $('#routeDesc').textContent = r.description + (r.modes[mode] ? ` — ${mode}: ${r.modes[mode]}` : '');
  } else $('#routeDesc').textContent = '';
  editor.defaults = r ? r.defaults : {};
}

function renderForm() {
  const box = clear($('#routeForm'));
  const r = route();
  if (!r) { box.append(el('div.empty-note', { text: 'Choose a model to see its keys.' })); return; }
  const mode = String(S.m.values.mode || '').toLowerCase();
  const req = new Set((r.required_by_mode[mode] || []).concat(r.required_by_mode['*'] || []));
  const known = new Set(['model', 'mode']);
  for (const g of r.groups) {
    const keys = g.keys.filter((k) => !known.has(k.key));
    keys.forEach((k) => known.add(k.key));
    if (!keys.length) continue;
    const nset = keys.filter((k) => k.key in S.m.values).length;
    const det = el('details.grp', { open: S.openGroups.has(g.name) || keys.some((k) => req.has(k.key)) });
    det.addEventListener('toggle', () => { if (det.open) S.openGroups.add(g.name); else S.openGroups.delete(g.name); });
    det.append(el('summary', {}, g.name, nset ? el('span.n', { text: `  ${nset} set` }) : null));
    const inner = el('div');
    for (const k of keys) inner.append(formRow(k, req.has(k.key)));
    det.append(inner);
    box.append(det);
  }
  const other = Object.keys(S.m.values).filter((k) => !known.has(k) && !isIndexedKey(k));
  if (other.length) {
    const det = el('details.grp', { open: true }, el('summary', {}, 'Other keys ', el('span.n', { text: `${other.length} not used by this model` })));
    const inner = el('div');
    for (const k of other) inner.append(formRow({ key: k, desc: 'not a key of this model (see diagnostics)', type: 'any' }, false));
    det.append(inner); box.append(det);
  }
}

function formRow(k, required) {
  const v = S.m.values;
  const cur = v[k.key];
  let ctrl;
  const commit = (val) => {
    const s = String(val ?? '').trim();
    if (s === '') delete v[k.key]; else v[k.key] = s;
    rowEl.classList.toggle('set', s !== '');
    const structural = ['structure', 'symmetry', 'nchord', 'nspan', 'clamp', ...PLATE_KEYS].includes(k.key);
    modelChanged({ props: structural, tree: structural });
    if (structural && k.key === 'structure') { editor.fit(); }
  };
  const dflt = k.default === null || k.default === undefined ? '' : Array.isArray(k.default) ? k.default.join(', ') : String(k.default);
  if (k.choices || k.type === 'bool') {
    const opts = k.type === 'bool' ? ['true', 'false'] : k.choices;
    ctrl = el('select', {}, el('option', { value: '', text: dflt ? `(default: ${dflt})` : '(unset)' }),
      opts.map((o) => el('option', { value: o, text: o, selected: String(cur ?? '').toLowerCase() === String(o).toLowerCase() })));
    if (cur !== undefined && !opts.some((o) => String(o).toLowerCase() === String(cur).toLowerCase())) {
      ctrl.append(el('option', { value: cur, text: `${cur} (?)`, selected: true }));
    }
    ctrl.addEventListener('change', () => commit(ctrl.value));
  } else {
    ctrl = el('input', { value: cur ?? '', placeholder: dflt, spellcheck: false });
    ctrl.addEventListener('change', () => commit(ctrl.value));
  }
  const rowEl = el('div.form-row' + (required ? '.req' : '') + (cur !== undefined ? '.set' : ''), { title: `${k.key} (${k.type}): ${k.desc}` },
    el('label', { text: k.key }), ctrl, el('div.help', { text: k.desc }));
  return rowEl;
}

// ---------------------------------------------------------------- validation
async function validate(show = true) {
  const text = $('#cfgText').value;
  try {
    const d = await api.validate(text, S.source);
    S.lastDiag = d;
    renderDiag(d);
    if (show) showPane('diag');
    return d;
  } catch (e) {
    if (show) toast(`validate: ${e.message}`, true);
    return null;
  }
}

function renderDiag(d) {
  const view = clear($('#diagView'));
  const { errors, warnings, notices } = d.summary;
  const cnt = $('#diagCount');
  cnt.textContent = errors ? String(errors) : warnings ? String(warnings) : '✓';
  cnt.className = 'count ' + (errors ? 'err' : warnings ? 'warn' : 'ok');
  view.append(el('div.verdict' + (d.ok ? '.ok' : '.bad'), {
    text: d.ok ? `Config is valid — ${warnings} warning(s), ${notices} notice(s). Ready to run.`
      : `${errors} error(s), ${warnings} warning(s) — fix errors before running.`,
  }));
  if (!d.diagnostics.length) return;
  const tb = el('tbody');
  for (const g of d.diagnostics) {
    const line = g.location && g.location.line;
    const tr = el('tr' + (line ? '.clickable' : ''), { onclick: line ? () => gotoLine(line) : null },
      el('td', {}, el('span.chip.' + g.severity, { text: g.severity })),
      el('td.mono', { text: g.code }),
      el('td', {}, el('div', { text: g.message }), g.hint ? el('div.hint-t', { text: `→ ${g.hint}` }) : null),
      el('td.mono', { text: g.field || '' }),
      el('td.mono', { text: line ? `line ${line}` : '' }));
    tb.append(tr);
  }
  view.append(el('table', {}, el('thead', {}, el('tr', {}, ['severity', 'code', 'message / hint', 'key', 'where'].map((h) => el('th', { text: h })))), tb));
}

function gotoLine(line) {
  showPane('config');
  const ta = $('#cfgText');
  const lines = ta.value.split('\n');
  let a = 0;
  for (let i = 0; i < line - 1 && i < lines.length; i++) a += lines[i].length + 1;
  const b = a + (lines[line - 1] || '').length;
  ta.focus(); ta.setSelectionRange(a, b);
  ta.scrollTop = Math.max(0, (line - 4) * 19.5);
}

function showPane(name) {
  for (const b of $$('#bottomSeg button')) b.classList.toggle('active', b.dataset.pane === name);
  for (const p of $$('#ws-model .pane')) p.hidden = p.dataset.pane !== name;
}

// ---------------------------------------------------------------- runs
function logLine(text) {
  const cls = /error|traceback|exception|failed/i.test(text) ? 'e' : /warn/i.test(text) ? 'w' : /^\[studio\]|done|finished|critical|flutter at/i.test(text) ? 's' : '';
  return cls ? el('span.' + cls, { text: text + '\n' }) : document.createTextNode(text + '\n');
}

async function run() {
  if (S.job && S.job.status === 'running') { toast('a run is already in progress'); return; }
  const text = $('#cfgText').value;
  const d = await validate(false);
  if (d && !d.ok) { showPane('diag'); toast('config has errors — see Diagnostics', true); return; }
  let name = $('#runName').value.trim() || slug(S.m.values.title) || 'model';
  name = name.replace(/[^A-Za-z0-9_-]+/g, '_');
  $('#runName').value = name;
  try {
    const r = await api.run(text, name, S.source);
    S.job = { id: r.id, name: r.name, next: 0, status: 'running', output_dir: r.output_dir };
    clear($('#logView'));
    showPane('log');
    setJobUi('running');
    toast(`started ${r.name} → ${r.output_dir}`);
    pollJob();
  } catch (e) { toast(`run: ${e.message}`, true); }
}

async function pollJob() {
  clearTimeout(S.jobTimer);
  const j = S.job; if (!j) return;
  try {
    const r = await api.job(j.id, j.next);
    const lv = $('#logView');
    const atBottom = lv.parentElement.scrollTop + lv.parentElement.clientHeight >= lv.parentElement.scrollHeight - 30;
    const frag = document.createDocumentFragment();
    for (const ln of r.lines) frag.append(logLine(ln));
    lv.append(frag);
    if (atBottom) lv.parentElement.scrollTop = lv.parentElement.scrollHeight;
    j.next = r.next; j.status = r.status; j.has_results = r.has_results;
    if (r.status === 'running' || r.status === 'starting') {
      $('#jobBadge').textContent = `running ${j.name} · ${r.elapsed}s`;
      S.jobTimer = setTimeout(pollJob, 600);
      return;
    }
    setJobUi(r.status, r);
  } catch (e) {
    setJobUi('failed');
    toast(`job: ${e.message}`, true);
  }
}

function setJobUi(status, r) {
  const badge = $('#jobBadge');
  const running = status === 'running';
  $('#btnRun').disabled = running;
  $('#btnCancel').hidden = !running;
  if (running) { badge.className = 'badge run'; badge.textContent = `running ${S.job.name}`; $('#btnOpenResult').hidden = true; return; }
  const ok = status === 'done';
  badge.className = 'badge ' + (ok ? 'ok' : 'err');
  badge.textContent = `${S.job ? S.job.name : ''}: ${status}${r ? ` (${r.elapsed}s)` : ''}`;
  $('#btnOpenResult').hidden = !(S.job && S.job.has_results);
  if (ok) toast(`run ${S.job.name} finished`);
  else toast(`run ${S.job ? S.job.name : ''} ${status.replace('_', ' ')}`, true);
  if (hooks.onRunFinished) hooks.onRunFinished(S.job);
}

// ---------------------------------------------------------------- init
export async function initModelTab(h = {}) {
  hooks = h;
  editor = new PlanformEditor($('#editorCanvas'), {
    onChange: (kind) => liveChange(kind),
    onSelect: (sel) => { S.sel = sel; renderTree(); renderProps(); },
    onCursor: (x, y) => { $('#cursorPos').textContent = `x ${x.toFixed(3)}  y ${y.toFixed(3)} m`; },
    onAdd: (kind, obj) => (kind === 'surface' ? addSurface(obj) : addStation(obj)),
    onHint: (t) => { $('#editorHint').textContent = t; },
    onError: (msg) => toast(msg, true),
    onConvertGhosts: convertGhosts,
  });
  editor.setTool('select');

  // toolbar
  for (const b of $$('#toolSeg button')) {
    b.addEventListener('click', () => {
      for (const x of $$('#toolSeg button')) x.classList.toggle('active', x === b);
      editor.setTool(b.dataset.tool);
    });
  }
  $('#optSnap').addEventListener('change', (e) => { editor.snap = e.target.checked; });
  $('#optSnapStep').addEventListener('change', (e) => { editor.snapStep = e.target.value; });
  $('#optBoxes').addEventListener('change', (e) => { editor.showBoxes = e.target.checked; editor.render(); });
  $('#optMirror').addEventListener('change', (e) => { editor.mirror = e.target.checked; editor.render(); });
  $('#btnFit').addEventListener('click', () => editor.fit());

  $('#btnAddSurface').addEventListener('click', () => addSurface());
  $('#btnAddControl').addEventListener('click', () => addControl());
  $('#btnAddStation').addEventListener('click', () => addStation());
  $('#btnNew').addEventListener('click', () => loadText(NEW_TEMPLATE, null));
  $('#btnReloadCfgs').addEventListener('click', loadConfigList);
  $('#cfgFilter').addEventListener('input', renderConfigList);

  for (const b of $$('#bottomSeg button')) b.addEventListener('click', () => showPane(b.dataset.pane));
  $('#cfgText').addEventListener('input', onTextInput);
  $('#btnValidate').addEventListener('click', () => validate(true));
  $('#btnRun').addEventListener('click', run);
  $('#btnCancel').addEventListener('click', async () => { if (S.job) { try { await api.cancel(S.job.id); } catch (e) { toast(e.message, true); } } });
  $('#btnOpenResult').addEventListener('click', () => { if (S.job && hooks.openResults) hooks.openResults(S.job.output_dir.split('/').pop()); });

  $('#selModel').addEventListener('change', (e) => {
    S.m.values.model = e.target.value;
    const r = route();
    if (r && S.m.values.mode && !(String(S.m.values.mode).toLowerCase() in r.modes)) S.m.values.mode = Object.keys(r.modes)[0];
    syncRouteSelectors(); modelChanged({ form: true, props: true });
  });
  $('#selMode').addEventListener('change', (e) => {
    if (e.target.value) S.m.values.mode = e.target.value; else delete S.m.values.mode;
    syncRouteSelectors(); modelChanged({ form: true });
  });

  window.addEventListener('keydown', (e) => {
    if ($('#ws-model').hidden) return;
    const tag = (document.activeElement && document.activeElement.tagName) || '';
    if (['INPUT', 'TEXTAREA', 'SELECT'].includes(tag)) return;
    const tool = (t) => { const b = $(`#toolSeg button[data-tool="${t}"]`); if (b) b.click(); };
    if (e.key === 'Escape') { if (!editor.cancelDraw()) tool('select'); }
    else if (e.key === 'Delete' || e.key === 'Backspace') { deleteSel(); e.preventDefault(); }
    else if (e.key === 'v' || e.key === 'V') tool('select');
    else if (e.key === 's' || e.key === 'S') tool('draw');
    else if (e.key === 'm' || e.key === 'M') tool('station');
    else if (e.key === 'f' || e.key === 'F') editor.fit();
  });

  const [models, presets] = await Promise.all([
    api.models().catch((e) => { toast(`models: ${e.message}`, true); return null; }),
    api.presets().catch(() => ({ presets: {} })),
  ]);
  S.models = models;
  S.presets = presets.presets || {};
  editor.setPresets(S.presets);
  await loadConfigList();
  const first = S.configs.find((c) => /config_goland_pk\.txt$/.test(c.path)) || S.configs.find((c) => c.group === 'configs');
  if (h.initialConfig) {
    try { const r = await api.config(h.initialConfig); loadText(r.text, r.path); return; } catch { /* fall through */ }
  }
  if (first) {
    try { const r = await api.config(first.path); loadText(r.text, r.path); return; } catch { /* fall through */ }
  }
  loadText(NEW_TEMPLATE, null);
}

export function modelTabShown() { editor && editor.render(); }
