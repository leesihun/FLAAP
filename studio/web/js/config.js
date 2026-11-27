// Flat `key value` config <-> editor model (mirrors flaap/launcher/parser.py grammar).
//
// Model state:
//   { header: ['% ...'], values: {key: raw}, comments: {key: '# ...'},
//     surfaces: [{name, rle:[x,y,z], rc, tle:[x,y,z], tc, nchord, nspan}],
//     controls: [{name, surface, hinge, eta0, eta1}],
//     stations: [{name, p:[x,y,z], a:[x,y,z]}] }

const RX = { surface: /^surface(\d+)$/, control: /^control(\d+)$/, station: /^station(\d+)$/ };
const FIRST = ['model', 'mode', 'title'];

export function emptyModel() {
  return { header: ['% FLAAP config - created in FLAAP Studio'], values: { model: 'flutter', mode: 'pk' },
    comments: {}, surfaces: [], controls: [], stations: [] };
}

export function num(s, dflt = NaN) {
  if (typeof s === 'number') return s;
  if (s === undefined || s === null) return dflt;
  const v = Number(String(s).trim());
  return Number.isFinite(v) && String(s).trim() !== '' ? v : dflt;
}

/** Compact number text (max 6 decimals, no trailing zeros). */
export function fnum(v) {
  if (!Number.isFinite(v)) return '0';
  const r = Math.round(v * 1e6) / 1e6;
  return Object.is(r, -0) ? '0' : String(r);
}

export function splitList(raw) {
  const s = String(raw ?? '').trim();
  if (!s) return [];
  return s.includes(',') ? s.split(',').map((x) => x.trim()) : s.split(/\s+/);
}

/** Parse flat config text. Returns {model, errors:[{line, msg}]}. */
export function parseConfig(text) {
  const m = { header: [], values: {}, comments: {}, surfaces: [], controls: [], stations: [] };
  const errors = [];
  const idx = { surface: [], control: [], station: [] };
  let seenKey = false;
  String(text).replace(/^﻿/, '').split(/\r?\n/).forEach((line, i) => {
    const ln = i + 1;
    const s = line.trim();
    if (!s) return;
    if (s.startsWith('%')) { if (!seenKey) m.header.push(s); return; }
    const hash = s.indexOf('#');
    const body = (hash >= 0 ? s.slice(0, hash) : s).trim();
    const comment = hash >= 0 ? s.slice(hash).trim() : '';
    if (!body) return;
    const sp = body.search(/\s/);
    if (sp < 0) { errors.push({ line: ln, msg: `expected "key value": ${body}` }); return; }
    seenKey = true;
    const key = body.slice(0, sp).toLowerCase();
    const raw = body.slice(sp).trim();
    for (const [kind, rx] of Object.entries(RX)) {
      const mm = key.match(rx);
      if (!mm) continue;
      const p = splitList(raw);
      const n = { surface: 11, control: 5, station: 7 }[kind];
      if (p.length !== n) { errors.push({ line: ln, msg: `${key} needs ${n} values` }); break; }
      let item;
      if (kind === 'surface') {
        item = { name: p[0], rle: [num(p[1], 0), num(p[2], 0), num(p[3], 0)], rc: num(p[4], 1),
          tle: [num(p[5], 0), num(p[6], 1), num(p[7], 0)], tc: num(p[8], 1), nchord: Math.round(num(p[9], 8)),
          nspan: Math.round(num(p[10], 10)) };
      } else if (kind === 'control') {
        item = { name: p[0], surface: p[1], hinge: num(p[2], 0.75), eta0: num(p[3], 0), eta1: num(p[4], 1) };
      } else {
        item = { name: p[0], p: [num(p[1], 0), num(p[2], 0), num(p[3], 0)], a: [num(p[4], 0), num(p[5], 1), num(p[6], 0)] };
      }
      if (comment) item.comment = comment;
      idx[kind].push([Number(mm[1]), item]);
      return;
    }
    if (Object.keys(RX).some((k) => RX[k].test(key))) { m.values[key] = raw; return; } // malformed indexed: keep raw
    if (key in m.values) errors.push({ line: ln, msg: `duplicate key ${key}` });
    m.values[key] = raw;
    if (comment) m.comments[key] = comment;
  });
  for (const k of Object.keys(idx)) {
    idx[k].sort((a, b) => a[0] - b[0]);
    m[k === 'surface' ? 'surfaces' : k === 'control' ? 'controls' : 'stations'] = idx[k].map((x) => x[1]);
  }
  return { model: m, errors };
}

/** Generate config text from the editor model. */
export function generateConfig(m) {
  const keys = [...FIRST.filter((k) => k in m.values), ...Object.keys(m.values).filter((k) => !FIRST.includes(k))];
  const rows = [];
  for (const k of keys) {
    const v = m.values[k];
    if (v === undefined || v === null || String(v).trim() === '') continue;
    rows.push([k, String(v).trim(), m.comments[k]]);
  }
  m.surfaces.forEach((s, i) => rows.push([`surface${i + 1}`,
    [s.name, ...s.rle.map(fnum), fnum(s.rc), ...s.tle.map(fnum), fnum(s.tc), s.nchord, s.nspan].join(', '), s.comment]));
  m.controls.forEach((c, i) => rows.push([`control${i + 1}`,
    [c.name, c.surface, fnum(c.hinge), fnum(c.eta0), fnum(c.eta1)].join(', '), c.comment]));
  m.stations.forEach((s, i) => rows.push([`station${i + 1}`,
    [s.name, ...s.p.map(fnum), ...s.a.map(fnum)].join(', '), s.comment]));
  const w = Math.max(13, ...rows.map((r) => r[0].length + 2));
  const lines = [...m.header];
  for (const [k, v, c] of rows) {
    let line = k.padEnd(w) + v;
    if (c) line = line.padEnd(Math.max(line.length + 3, w + 24)) + c;
    lines.push(line);
  }
  return lines.join('\n') + '\n';
}

// ----------------------------------------------------------------- plate structure
export function getPlate(m) {
  const v = m.values;
  const vec = (k, d) => { const p = splitList(v[k]); return p.length >= 2 ? [num(p[0], d[0]), num(p[1], d[1]), num(p[2], 0)] : d; };
  return { name: 'plate', rle: vec('plate_root_le', [0, 0, 0]), rc: num(v.plate_root_chord, 0.5),
    tle: vec('plate_tip_le', [0.2, 1.0, 0]), tc: num(v.plate_tip_chord, 0.3),
    nchord: Math.round(num(v.plate_nchord, 8)), nspan: Math.round(num(v.plate_nspan, 12)) };
}

export function setPlate(m, pf) {
  const v = m.values;
  v.plate_root_le = pf.rle.map(fnum).join(', ');
  v.plate_root_chord = fnum(pf.rc);
  v.plate_tip_le = pf.tle.map(fnum).join(', ');
  v.plate_tip_chord = fnum(pf.tc);
  if (pf.nchord !== undefined && v.plate_nchord !== undefined) v.plate_nchord = String(pf.nchord);
  if (pf.nspan !== undefined && v.plate_nspan !== undefined) v.plate_nspan = String(pf.nspan);
}

// ----------------------------------------------------------------- planform geometry
/** Point on a planform at span fraction eta and chord fraction xi. */
export function pfPoint(s, eta, xi) {
  const le = [0, 1, 2].map((i) => s.rle[i] + eta * (s.tle[i] - s.rle[i]));
  const ch = s.rc + eta * (s.tc - s.rc);
  return [le[0] + xi * ch, le[1], le[2]];
}

export function pfDerived(s) {
  const span = Math.hypot(s.tle[1] - s.rle[1], s.tle[2] - s.rle[2]);
  const area = 0.5 * (s.rc + s.tc) * span;
  const sweepLE = span > 0 ? Math.atan2(s.tle[0] - s.rle[0], span) * 180 / Math.PI : 0;
  const xq = (e) => pfPoint(s, e, 0.25)[0];
  const sweepQC = span > 0 ? Math.atan2(xq(1) - xq(0), span) * 180 / Math.PI : 0;
  const dihedral = span > 0 ? Math.atan2(s.tle[2] - s.rle[2], Math.abs(s.tle[1] - s.rle[1])) * 180 / Math.PI : 0;
  const taper = s.rc > 0 ? s.tc / s.rc : 0;
  const ar = area > 0 ? span * span / area : 0;
  const mac = (2 / 3) * s.rc * (1 + taper + taper * taper) / (1 + taper || 1);
  return { span, area, sweepLE, sweepQC, dihedral, taper, ar, mac, boxes: s.nchord * s.nspan };
}

/** Box indices that belong to a control (same rule as flaap.geometry.build_mesh, linear spacing). */
export function controlBoxes(s, c) {
  const out = [];
  for (let j = 0; j < s.nspan; j++) {
    const etaC = (j + 0.5) / s.nspan;
    if (etaC < c.eta0 - 1e-9 || etaC > c.eta1 + 1e-9) continue;
    for (let i = 0; i < s.nchord; i++) if ((i + 0.5) / s.nchord > c.hinge) out.push([j, i]);
  }
  return out;
}

export const isIndexedKey = (k) => Object.values(RX).some((rx) => rx.test(k));
