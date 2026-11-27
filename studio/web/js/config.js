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

// ----------------------------------------------------------------- structure orientation
// Mirrors flaap/structures/transform.py: struct_mirror (x|y|z), struct_rotate (axis, deg, ...),
// struct_translate (dx, dy, dz); applied in that order.
const AX = { x: 0, y: 1, z: 2 };
const I3 = () => [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const mul3 = (A, B) => A.map((r) => [0, 1, 2].map((j) => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const mulv = (A, v) => A.map((r) => r[0] * v[0] + r[1] * v[1] + r[2] * v[2]);

export function rotMatrix(axis, deg) {
  const a = deg * Math.PI / 180;
  let c = Math.cos(a), s = Math.sin(a);
  if (Math.abs(c) < 1e-12) c = 0;
  if (Math.abs(s) < 1e-12) s = 0;
  const i = AX[axis], j = (i + 1) % 3, k = (i + 2) % 3;
  const R = I3();
  R[j][j] = c; R[j][k] = -s; R[k][j] = s; R[k][k] = c;
  return R;
}

export const orientIdentity = (T) => !T.mirror && !T.rotations.length && !T.translate.some((t) => t);

export function orientFromValues(v) {
  const T = { mirror: null, rotations: [], translate: [0, 0, 0] };
  const mir = String(v.struct_mirror ?? '').trim().toLowerCase();
  if (mir in AX) T.mirror = mir;
  const rot = splitList(v.struct_rotate);
  for (let i = 0; i + 1 < rot.length; i += 2) {
    const ax = rot[i].toLowerCase(), d = num(rot[i + 1], NaN);
    if (ax in AX && Number.isFinite(d)) T.rotations.push([ax, d]);
  }
  const tr = splitList(v.struct_translate).map((x) => num(x, 0));
  if (tr.length === 3) T.translate = tr;
  return T;
}

export function orientToValues(T, v) {
  const set = (k, val) => { if (val === null) delete v[k]; else v[k] = val; };
  set('struct_mirror', T.mirror || null);
  set('struct_rotate', T.rotations.length ? T.rotations.map(([a, d]) => `${a}, ${fnum(d)}`).join(', ') : null);
  set('struct_translate', T.translate.some((t) => t) ? T.translate.map(fnum).join(', ') : null);
}

export function orientMatrix(T) {
  let Q = I3();
  if (T.mirror) Q[AX[T.mirror]][AX[T.mirror]] = -1;
  for (const [a, d] of T.rotations) Q = mul3(rotMatrix(a, d), Q);
  return Q;
}

export const orientPoint = (T, p) => mulv(orientMatrix(T), [p[0], p[1], p[2] || 0]).map((x, i) => x + T.translate[i]);

export function orientDescribe(T) {
  const p = [];
  if (T.mirror) p.push(`mirror ${T.mirror}`);
  for (const [a, d] of T.rotations) p.push(`rotate ${a} ${fnum(d)}°`);
  if (T.translate.some((t) => t)) p.push(`shift ${T.translate.map(fnum).join(', ')}`);
  return p.join(' · ') || 'as modelled';
}

/** Transformed planform (same fields as a surface); null when the chord direction would turn
 *  (a Surface has its chord along +x, so a preset aero surface cannot follow). */
export function orientSurface(s, T) {
  const Q = orientMatrix(T);
  if (Math.abs(Q[0][0] - 1) > 1e-9 || Math.abs(Q[1][0]) > 1e-9 || Math.abs(Q[2][0]) > 1e-9) return null;
  let r = orientPoint(T, s.rle), t = orientPoint(T, s.tle), rc = s.rc, tc = s.tc;
  let controls = (s.controls || []).map((c) => ({ ...c }));
  if (t[1] < r[1] - 1e-12) {
    [r, t] = [t, r]; [rc, tc] = [tc, rc];
    controls = controls.map((c) => ({ ...c, eta0: 1 - c.eta1, eta1: 1 - c.eta0 }));
  }
  return { ...s, rle: r, tle: t, rc, tc, controls };
}

/** Rotation sequence (<= 3 quarter turns) that maps the file frame to the aero frame, given the
 *  direction the leading->trailing edge and the root->tip span point to in the file.
 *  Returns null for an invalid pair of directions. */
export function solveOrientation(chordDir, spanDir) {
  const unit = (s) => { const v = [0, 0, 0]; v[AX[s[1]]] = s[0] === '-' ? -1 : 1; return v; };
  const a = unit(chordDir), b = unit(spanDir);
  if (a.some((x, i) => x && b[i])) return null;                     // same axis
  const c = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
  const want = [a, b, c];                                          // Q rows
  const same = (Q) => Q.every((r, i) => r.every((x, j) => Math.abs(x - want[i][j]) < 1e-9));
  if (same(I3())) return [];
  const moves = [];
  for (const ax of ['x', 'y', 'z']) for (const d of [90, 180, -90]) moves.push([ax, d]);
  let frontier = [[]];
  for (let depth = 1; depth <= 3; depth++) {
    const next = [];
    for (const seq of frontier) for (const m of moves) {
      const s2 = [...seq, m];
      let Q = I3();
      for (const [ax, d] of s2) Q = mul3(rotMatrix(ax, d), Q);
      if (same(Q)) return s2;
      next.push(s2);
    }
    frontier = next;
  }
  return null;
}
