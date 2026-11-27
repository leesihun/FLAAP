// Planform editor: top view of the lifting surfaces (x = chord, right; y = span, up).
//
// Works directly on the config model object (see config.js). Tools:
//   select  - click to select, drag corner handles / body, drag empty space to pan
//   draw    - four clicks: root LE, tip LE, tip TE, root TE  (Esc cancels)
//   station - click to drop a monitor station
// Wheel zooms around the cursor, middle/right drag pans.

import { fitCanvas, SURF_COLORS } from './util.js';
import { getPlate, setPlate, pfPoint, controlBoxes, num, splitList } from './config.js';
import { niceStep } from './plot.js';

const HANDLE = 6;
const AMBER = '#b76b2a', BLUE = '#3f6f91', INK = '#17211e', GREEN = '#126b57';

function inPoly(pts, x, y) {
  let inside = false;
  for (let a = 0, b = pts.length - 1; a < pts.length; b = a++) {
    const [xa, ya] = pts[a], [xb, yb] = pts[b];
    if ((ya > y) !== (yb > y) && x < (xb - xa) * (y - ya) / (yb - ya) + xa) inside = !inside;
  }
  return inside;
}
const corners = (s) => [[s.rle[0], s.rle[1]], [s.tle[0], s.tle[1]], [s.tle[0] + s.tc, s.tle[1]], [s.rle[0] + s.rc, s.rle[1]]];
const r6 = (v) => Math.round(v * 1e6) / 1e6;

export class PlanformEditor {
  constructor(canvas, opts = {}) {
    this.cv = canvas;
    this.wrap = canvas.parentElement;
    this.opts = opts; // onChange(kind), onSelect(sel), onCursor(x, y), onAdd(kind, obj), onHint(text)
    this.m = null;
    this.presets = {};
    this.sel = null;
    this.tool = 'select';
    this.snap = true; this.snapStep = 'auto';
    this.showBoxes = true; this.mirror = false;
    this.view = { scale: 200, ox: 80, oy: 400 };
    this.draw = null; // {pts: [[x,y],...]}
    this.hover = null;
    this.mouse = null;
    this.defaults = {}; // route defaults (structure, symmetry, nchord, ...)
    this._bind();
    new ResizeObserver(() => this.render()).observe(this.wrap);
  }

  // ------------------------------------------------------------ public API
  setModel(m) { this.m = m; this.render(); }
  setPresets(p) { this.presets = p || {}; this.render(); }
  select(sel, notify = false) {
    this.sel = sel; this.render();
    if (notify && this.opts.onSelect) this.opts.onSelect(sel);
  }
  setTool(t) {
    this.tool = t; this.draw = null;
    this.cv.style.cursor = t === 'select' ? 'default' : 'crosshair';
    this._hint(); this.render();
  }

  /** Read-only planforms that a structure preset brings when no surfaceN is given. */
  ghosts() {
    if (!this.m || this.m.surfaces.length) return [];
    const st = this.structure();
    const p = this.presets[st];
    if (!p) return [];
    const nc = Math.round(num(this.m.values.nchord, NaN)), ns = Math.round(num(this.m.values.nspan, NaN));
    return p.map((s) => ({ ...s, nchord: Number.isFinite(nc) ? nc : s.nchord, nspan: Number.isFinite(ns) ? ns : s.nspan }));
  }
  structure() {
    return String(this.m.values.structure || this.defaults.structure || '').trim().toLowerCase();
  }
  plate() {
    if (!this.m) return null;
    return this.structure() === 'plate' ? getPlate(this.m) : null;
  }
  /** Aero surfaces actually used (explicit, else preset, else plate planform). */
  aeroSurfaces() {
    if (!this.m) return [];
    if (this.m.surfaces.length) return this.m.surfaces;
    const g = this.ghosts(); if (g.length) return g;
    const p = this.plate();
    if (p) {
      return [{ ...p, nchord: Math.round(num(this.m.values.nchord, 8)), nspan: Math.round(num(this.m.values.nspan, 12)), implicit: true }];
    }
    return [];
  }

  fit() {
    const pts = [];
    if (this.m) {
      for (const s of [...this.m.surfaces, ...this.ghosts()]) pts.push(...corners(s));
      const p = this.plate(); if (p) pts.push(...corners(p));
      for (const st of this.m.stations) pts.push([st.p[0], st.p[1]]);
    }
    if (this.mirror || this._sym()) pts.push(...pts.map(([x, y]) => [x, -y]));
    if (!pts.length) pts.push([-0.2, -0.2], [2, 2]);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const [x, y] of pts) { x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
    const r = this.wrap.getBoundingClientRect();
    const w = Math.max(100, r.width), h = Math.max(100, r.height);
    const dx = Math.max(x1 - x0, 1e-3), dy = Math.max(y1 - y0, 1e-3);
    const sc = Math.min((w - 120) / dx, (h - 90) / dy);
    this.view.scale = Math.max(1e-3, sc);
    this.view.ox = w / 2 - sc * (x0 + x1) / 2;
    this.view.oy = h / 2 + sc * (y0 + y1) / 2;
    this.render();
  }

  // ------------------------------------------------------------ coordinates
  w2s(x, y) { return [this.view.ox + x * this.view.scale, this.view.oy - y * this.view.scale]; }
  s2w(px, py) { return [(px - this.view.ox) / this.view.scale, (this.view.oy - py) / this.view.scale]; }
  gridStep() { return niceStep(110 / this.view.scale, 1); }
  snapVal(v) {
    if (!this.snap) return r6(v);
    const st = this.snapStep === 'auto' ? this.gridStep() / 5 : Number(this.snapStep);
    return r6(Math.round(v / st) * st);
  }
  _sym() { return !!this.m && /^(sym|anti)/i.test(String(this.m.values.symmetry || this.defaults.symmetry || '')); }

  // ------------------------------------------------------------ rendering
  render() {
    const { w, h, ctx } = fitCanvas(this.cv);
    this.W = w; this.H = h;
    ctx.fillStyle = '#fbfaf6'; ctx.fillRect(0, 0, w, h);
    this._grid(ctx, w, h);
    if (!this.m) return;
    const mirror = this.mirror || this._sym();
    const passes = mirror ? [1, -1] : [1];
    for (const sgn of passes) {
      ctx.globalAlpha = sgn > 0 ? 1 : 0.35;
      this._plate(ctx, sgn);
      this.ghosts().forEach((s, i) => this._surface(ctx, s, i, sgn, true));
      const imp = this.aeroSurfaces().find((s) => s.implicit);
      if (imp) this._surface(ctx, imp, 0, sgn, true, true);
      this.m.surfaces.forEach((s, i) => this._surface(ctx, s, i, sgn, false));
      ctx.globalAlpha = 1;
    }
    this._clamp(ctx, w);
    this._stations(ctx);
    this._selection(ctx);
    this._drawPreview(ctx);
  }

  _grid(ctx, w, h) {
    const step = this.gridStep(), minor = step / 5;
    const [x0, y1] = this.s2w(0, 0), [x1, y0] = this.s2w(w, h);
    const line = (x, y, X, Y) => { ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(X, Y); ctx.stroke(); };
    ctx.lineWidth = 1;
    if (this.view.scale * minor > 7) {
      ctx.strokeStyle = '#f0ede5';
      for (let x = Math.ceil(x0 / minor) * minor; x <= x1; x += minor) { const [px] = this.w2s(x, 0); line(Math.round(px) + 0.5, 0, Math.round(px) + 0.5, h); }
      for (let y = Math.ceil(y0 / minor) * minor; y <= y1; y += minor) { const [, py] = this.w2s(0, y); line(0, Math.round(py) + 0.5, w, Math.round(py) + 0.5); }
    }
    ctx.strokeStyle = '#e2ded3';
    ctx.fillStyle = '#7d8884'; ctx.font = '10.5px Consolas, monospace';
    const dec = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    for (let x = Math.ceil(x0 / step) * step; x <= x1; x += step) {
      const [px] = this.w2s(x, 0); line(Math.round(px) + 0.5, 0, Math.round(px) + 0.5, h);
      ctx.fillText(x.toFixed(dec), px + 3, h - 6);
    }
    for (let y = Math.ceil(y0 / step) * step; y <= y1; y += step) {
      const [, py] = this.w2s(0, y); line(0, Math.round(py) + 0.5, w, Math.round(py) + 0.5);
      ctx.fillText(y.toFixed(dec), 4, py - 3);
    }
    // axes
    ctx.strokeStyle = '#9aa39f'; ctx.lineWidth = 1.2;
    const [ax, ay] = this.w2s(0, 0);
    line(ax, 0, ax, h); line(0, ay, w, ay);
    ctx.fillStyle = '#4a5753'; ctx.font = '600 11px Inter, system-ui, sans-serif';
    ctx.fillText('x (chord) →', w - 84, ay - 6);
    ctx.save(); ctx.translate(ax + 14, 14); ctx.fillText('y (span) ↑', 0, 0); ctx.restore();
    ctx.fillStyle = '#7d8884'; ctx.font = '10.5px Inter, system-ui, sans-serif';
    ctx.fillText('flow →', 12, 18);
  }

  _poly(ctx, pts, sgn) {
    ctx.beginPath();
    pts.forEach(([x, y], i) => { const [px, py] = this.w2s(x, sgn * y); i ? ctx.lineTo(px, py) : ctx.moveTo(px, py); });
    ctx.closePath();
  }

  _surface(ctx, s, i, sgn, ghost, implicit = false) {
    const col = ghost ? '#7d8884' : SURF_COLORS[i % SURF_COLORS.length];
    const c = corners(s);
    // controls on this surface
    const ctrls = this.m.controls.map((cc, k) => [cc, k]).filter(([cc]) => cc.surface === s.name);
    for (const ctrl of (s.controls || [])) ctrls.push([{ ...ctrl, surface: s.name }, -1]);
    ctx.save();
    this._poly(ctx, c, sgn);
    ctx.fillStyle = col; ctx.globalAlpha *= ghost ? 0.07 : 0.12; ctx.fill();
    ctx.restore();
    // control boxes
    for (const [cc, k] of ctrls) {
      const boxes = controlBoxes(s, cc);
      ctx.save();
      const selC = this.sel && this.sel.kind === 'control' && this.sel.i === k;
      ctx.fillStyle = AMBER; ctx.globalAlpha *= selC ? 0.5 : 0.3;
      for (const [j, ii] of boxes) {
        const q = [[j / s.nspan, ii / s.nchord], [(j + 1) / s.nspan, ii / s.nchord], [(j + 1) / s.nspan, (ii + 1) / s.nchord], [j / s.nspan, (ii + 1) / s.nchord]]
          .map(([e, xi]) => pfPoint(s, e, xi));
        this._poly(ctx, q, sgn); ctx.fill();
      }
      ctx.restore();
      // hinge line over eta0..eta1
      const a = pfPoint(s, cc.eta0, cc.hinge), b = pfPoint(s, cc.eta1, cc.hinge);
      const [ax, ay] = this.w2s(a[0], sgn * a[1]), [bx, by] = this.w2s(b[0], sgn * b[1]);
      ctx.strokeStyle = AMBER; ctx.lineWidth = selC ? 2.5 : 1.6; ctx.setLineDash([6, 4]);
      ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by); ctx.stroke(); ctx.setLineDash([]);
      if (sgn > 0) {
        ctx.fillStyle = AMBER; ctx.font = '600 11px Inter, system-ui, sans-serif';
        const t = pfPoint(s, 0.5 * (cc.eta0 + cc.eta1), 1);
        const [tx, ty] = this.w2s(t[0], t[1]);
        ctx.fillText(cc.name, tx + 6, ty + 4);
      }
    }
    // box grid
    if (this.showBoxes && s.nchord > 0 && s.nspan > 0 && s.nchord * s.nspan < 20000) {
      ctx.save();
      ctx.strokeStyle = col; ctx.globalAlpha *= 0.4; ctx.lineWidth = 0.7;
      ctx.beginPath();
      for (let i2 = 1; i2 < s.nchord; i2++) {
        const a = pfPoint(s, 0, i2 / s.nchord), b = pfPoint(s, 1, i2 / s.nchord);
        const [ax, ay] = this.w2s(a[0], sgn * a[1]), [bx, by] = this.w2s(b[0], sgn * b[1]);
        ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
      }
      for (let j = 1; j < s.nspan; j++) {
        const a = pfPoint(s, j / s.nspan, 0), b = pfPoint(s, j / s.nspan, 1);
        const [ax, ay] = this.w2s(a[0], sgn * a[1]), [bx, by] = this.w2s(b[0], sgn * b[1]);
        ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
      }
      ctx.stroke();
      ctx.restore();
    }
    // outline
    this._poly(ctx, c, sgn);
    const selS = !ghost && this.sel && this.sel.kind === 'surface' && this.sel.i === i;
    const hovS = this.hover && ((this.hover.kind === 'surface' && !ghost && this.hover.i === i) || (this.hover.kind === 'ghost' && ghost && this.hover.i === i));
    ctx.strokeStyle = col; ctx.lineWidth = selS ? 2.6 : hovS ? 2 : 1.5;
    if (ghost) ctx.setLineDash([5, 4]);
    ctx.stroke(); ctx.setLineDash([]);
    // quarter-chord line
    if (sgn > 0) {
      const a = pfPoint(s, 0, 0.25), b = pfPoint(s, 1, 0.25);
      const [ax, ay] = this.w2s(a[0], a[1]), [bx, by] = this.w2s(b[0], b[1]);
      ctx.strokeStyle = col; ctx.globalAlpha *= 0.5; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by); ctx.stroke(); ctx.setLineDash([]);
      ctx.globalAlpha = 1;
      const [lx, ly] = this.w2s(s.rle[0] + s.rc / 2, s.rle[1]);
      ctx.fillStyle = ghost ? '#4a5753' : col; ctx.font = '600 12px Inter, system-ui, sans-serif';
      const label = implicit ? `${s.name} (aero = plate planform)` : ghost ? `${s.name} (preset)` : s.name;
      ctx.textAlign = 'center'; ctx.fillText(label, lx, ly + 16); ctx.textAlign = 'left';
    }
  }

  _plate(ctx, sgn) {
    const p = this.plate();
    if (!p) return;
    const c = corners(p);
    if (sgn > 0) {
      // FE grid
      ctx.save();
      ctx.strokeStyle = BLUE; ctx.globalAlpha *= 0.35; ctx.lineWidth = 0.7; ctx.setLineDash([1.5, 2.5]);
      ctx.beginPath();
      for (let i = 1; i < p.nchord; i++) {
        const a = pfPoint(p, 0, i / p.nchord), b = pfPoint(p, 1, i / p.nchord);
        const [ax, ay] = this.w2s(a[0], a[1]), [bx, by] = this.w2s(b[0], b[1]); ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
      }
      for (let j = 1; j < p.nspan; j++) {
        const a = pfPoint(p, j / p.nspan, 0), b = pfPoint(p, j / p.nspan, 1);
        const [ax, ay] = this.w2s(a[0], a[1]), [bx, by] = this.w2s(b[0], b[1]); ctx.moveTo(ax, ay); ctx.lineTo(bx, by);
      }
      ctx.stroke(); ctx.restore();
    }
    this._poly(ctx, c, sgn);
    const selP = this.sel && this.sel.kind === 'plate';
    ctx.strokeStyle = BLUE; ctx.lineWidth = selP ? 2.6 : 1.6; ctx.setLineDash([8, 3, 2, 3]); ctx.stroke(); ctx.setLineDash([]);
    if (sgn > 0) {
      const [lx, ly] = this.w2s(p.tle[0] + p.tc / 2, p.tle[1]);
      ctx.fillStyle = BLUE; ctx.font = '600 11.5px Inter, system-ui, sans-serif'; ctx.textAlign = 'center';
      ctx.fillText(`plate FE ${p.nchord}×${p.nspan}`, lx, ly - 8); ctx.textAlign = 'left';
    }
  }

  _clamp(ctx, w) {
    const p = this.plate(); if (!p) return;
    const raw = String(this.m.values.clamp ?? 'y, 0.0').trim().toLowerCase();
    if (raw === 'none' || raw === 'free') return;
    const parts = splitList(raw);
    if (parts[0] !== 'y') return;
    const yv = num(parts[1], 0);
    const x0 = Math.min(p.rle[0], p.tle[0]) - 0.1 * p.rc, x1 = Math.max(p.rle[0] + p.rc, p.tle[0] + p.tc) + 0.1 * p.rc;
    const [ax, ay] = this.w2s(x0, yv), [bx] = this.w2s(x1, yv);
    ctx.strokeStyle = INK; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, ay); ctx.stroke();
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let x = ax; x < bx; x += 7) { ctx.moveTo(x, ay); ctx.lineTo(x - 6, ay + 7); }
    ctx.stroke();
    ctx.fillStyle = INK; ctx.font = '10.5px Inter, system-ui, sans-serif'; ctx.fillText('clamp', bx + 4, ay + 4);
  }

  _stations(ctx) {
    this.m.stations.forEach((st, i) => {
      const [px, py] = this.w2s(st.p[0], st.p[1]);
      const sel = this.sel && this.sel.kind === 'station' && this.sel.i === i;
      // axis arrow (projected)
      const L = 28, ax = st.a[0], ay = st.a[1], n = Math.hypot(ax, ay);
      if (n > 1e-9) {
        const ex = px + (ax / n) * L, ey = py - (ay / n) * L;
        ctx.strokeStyle = '#6a4ce2'; ctx.lineWidth = 1.5;
        ctx.beginPath(); ctx.moveTo(px, py); ctx.lineTo(ex, ey); ctx.stroke();
        const ang = Math.atan2(ey - py, ex - px);
        ctx.beginPath(); ctx.moveTo(ex, ey);
        ctx.lineTo(ex - 7 * Math.cos(ang - 0.4), ey - 7 * Math.sin(ang - 0.4));
        ctx.lineTo(ex - 7 * Math.cos(ang + 0.4), ey - 7 * Math.sin(ang + 0.4)); ctx.closePath();
        ctx.fillStyle = '#6a4ce2'; ctx.fill();
      }
      const r = sel ? 8 : 6;
      ctx.beginPath(); ctx.moveTo(px, py - r); ctx.lineTo(px + r, py); ctx.lineTo(px, py + r); ctx.lineTo(px - r, py); ctx.closePath();
      ctx.fillStyle = '#6a4ce2'; ctx.fill(); ctx.strokeStyle = '#fbfaf6'; ctx.lineWidth = 2; ctx.stroke();
      ctx.fillStyle = INK; ctx.font = '600 11px Inter, system-ui, sans-serif';
      ctx.fillText(st.name, px + 10, py - 8);
    });
  }

  /** Handles of the selected object: [{id, x, y}] in world coordinates. */
  handles() {
    const s = this.sel; if (!s || !this.m) return [];
    let pf = null;
    if (s.kind === 'surface') pf = this.m.surfaces[s.i];
    else if (s.kind === 'plate') pf = this.plate();
    if (pf) {
      const c = corners(pf);
      return [{ id: 'rle', x: c[0][0], y: c[0][1] }, { id: 'tle', x: c[1][0], y: c[1][1] },
        { id: 'tte', x: c[2][0], y: c[2][1] }, { id: 'rte', x: c[3][0], y: c[3][1] }];
    }
    if (s.kind === 'control') {
      const c = this.m.controls[s.i]; if (!c) return [];
      const sf = this.m.surfaces.find((q) => q.name === c.surface); if (!sf) return [];
      const em = 0.5 * (c.eta0 + c.eta1);
      const h = pfPoint(sf, em, c.hinge), a = pfPoint(sf, c.eta0, 0.5 * (1 + c.hinge)), b = pfPoint(sf, c.eta1, 0.5 * (1 + c.hinge));
      return [{ id: 'hinge', x: h[0], y: h[1] }, { id: 'eta0', x: a[0], y: a[1] }, { id: 'eta1', x: b[0], y: b[1] }];
    }
    return [];
  }

  _selection(ctx) {
    for (const hd of this.handles()) {
      const [px, py] = this.w2s(hd.x, hd.y);
      const hot = this.hover && this.hover.handle === hd.id;
      ctx.fillStyle = hot ? '#ccf48d' : '#fff'; ctx.strokeStyle = GREEN; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.rect(px - HANDLE, py - HANDLE, 2 * HANDLE, 2 * HANDLE); ctx.fill(); ctx.stroke();
    }
  }

  _drawPreview(ctx) {
    if (this.tool !== 'draw' || !this.draw) return;
    const pts = [...this.draw.pts];
    const mp = this.mouse ? this._drawPoint(this.mouse) : null;
    if (mp) pts.push(mp);
    ctx.strokeStyle = GREEN; ctx.lineWidth = 1.6; ctx.setLineDash([5, 4]);
    ctx.beginPath();
    pts.forEach(([x, y], i) => { const [px, py] = this.w2s(x, y); i ? ctx.lineTo(px, py) : ctx.moveTo(px, py); });
    if (pts.length === 4) ctx.closePath();
    ctx.stroke(); ctx.setLineDash([]);
    for (const [x, y] of pts) {
      const [px, py] = this.w2s(x, y);
      ctx.beginPath(); ctx.arc(px, py, 4, 0, 2 * Math.PI); ctx.fillStyle = GREEN; ctx.fill();
    }
  }

  /** Constrained draw point: 3rd point shares the tip y, 4th the root y. */
  _drawPoint([x, y]) {
    const n = this.draw ? this.draw.pts.length : 0;
    if (n === 2) return [Math.max(x, this.draw.pts[1][0] + 1e-3), this.draw.pts[1][1]];
    if (n === 3) return [Math.max(x, this.draw.pts[0][0] + 1e-3), this.draw.pts[0][1]];
    return [x, y];
  }

  _hint() {
    let t = '';
    if (this.tool === 'draw') {
      const n = this.draw ? this.draw.pts.length : 0;
      t = ['Draw surface: click ROOT leading edge', 'click TIP leading edge', 'click TIP trailing edge', 'click ROOT trailing edge'][n] + '   (Esc cancels)';
    } else if (this.tool === 'station') t = 'Click to place a monitor station (axis +y). Esc: back to select';
    else t = 'Drag handles to reshape · drag body to move · drag empty space / right-drag to pan · wheel to zoom';
    if (this.opts.onHint) this.opts.onHint(t);
  }

  // ------------------------------------------------------------ hit testing
  hit(px, py) {
    if (!this.m) return null;
    const [x, y] = this.s2w(px, py);
    for (const hd of this.handles()) {
      const [hx, hy] = this.w2s(hd.x, hd.y);
      if (Math.abs(hx - px) <= HANDLE + 2 && Math.abs(hy - py) <= HANDLE + 2) return { ...this.sel, handle: hd.id };
    }
    for (let i = this.m.stations.length - 1; i >= 0; i--) {
      const [sx, sy] = this.w2s(this.m.stations[i].p[0], this.m.stations[i].p[1]);
      if (Math.hypot(sx - px, sy - py) < 9) return { kind: 'station', i };
    }
    for (let k = this.m.controls.length - 1; k >= 0; k--) {
      const c = this.m.controls[k];
      const sf = this.m.surfaces.find((q) => q.name === c.surface); if (!sf) continue;
      const poly = [pfPoint(sf, c.eta0, c.hinge), pfPoint(sf, c.eta1, c.hinge), pfPoint(sf, c.eta1, 1), pfPoint(sf, c.eta0, 1)];
      if (inPoly(poly, x, y)) return { kind: 'control', i: k };
    }
    for (let i = this.m.surfaces.length - 1; i >= 0; i--) if (inPoly(corners(this.m.surfaces[i]), x, y)) return { kind: 'surface', i };
    const p = this.plate(); if (p && inPoly(corners(p), x, y)) return { kind: 'plate' };
    const g = this.ghosts();
    for (let i = g.length - 1; i >= 0; i--) if (inPoly(corners(g[i]), x, y)) return { kind: 'ghost', i };
    return null;
  }

  // ------------------------------------------------------------ editing
  _apply(drag, wx, wy) {
    const o = drag.orig;
    const sx = this.snapVal(wx), sy = this.snapVal(wy);
    const sel = drag.target;
    if (sel.kind === 'station') {
      const st = this.m.stations[sel.i];
      st.p[0] = this.snapVal(o.p[0] + wx - drag.x0); st.p[1] = this.snapVal(o.p[1] + wy - drag.y0);
      return 'station';
    }
    if (sel.kind === 'control') {
      const c = this.m.controls[sel.i];
      const sf = this.m.surfaces.find((q) => q.name === c.surface); if (!sf) return null;
      const dy = sf.tle[1] - sf.rle[1];
      const eta = Math.abs(dy) > 1e-12 ? (wy - sf.rle[1]) / dy : 0;
      const cl = (v) => Math.max(0, Math.min(1, r6(Math.round(v * 200) / 200)));
      if (drag.handle === 'hinge' || !drag.handle) {
        const e = Math.max(0, Math.min(1, eta));
        const le = pfPoint(sf, e, 0), ch = sf.rc + e * (sf.tc - sf.rc);
        c.hinge = Math.min(0.98, Math.max(0.02, cl((wx - le[0]) / ch)));
      } else if (drag.handle === 'eta0') c.eta0 = Math.min(cl(eta), c.eta1 - 0.005);
      else if (drag.handle === 'eta1') c.eta1 = Math.max(cl(eta), c.eta0 + 0.005);
      return 'control';
    }
    const isPlate = sel.kind === 'plate';
    const s = isPlate ? getPlate(this.m) : this.m.surfaces[sel.i];
    if (!s) return null;
    const minC = 1e-3;
    if (!drag.handle) {
      const dx = this.snapVal(o.rle[0] + wx - drag.x0) - o.rle[0];
      const dy = this.snapVal(o.rle[1] + wy - drag.y0) - o.rle[1];
      s.rle = [r6(o.rle[0] + dx), r6(o.rle[1] + dy), o.rle[2]];
      s.tle = [r6(o.tle[0] + dx), r6(o.tle[1] + dy), o.tle[2]];
    } else if (drag.handle === 'rle') {
      const te = o.rle[0] + o.rc;
      const x = Math.min(sx, te - minC);
      s.rle = [x, sy, o.rle[2]]; s.rc = r6(te - x);
    } else if (drag.handle === 'rte') {
      s.rc = r6(Math.max(minC, sx - o.rle[0])); s.rle = [o.rle[0], sy, o.rle[2]];
    } else if (drag.handle === 'tle') {
      const te = o.tle[0] + o.tc;
      const x = Math.min(sx, te - minC);
      s.tle = [x, sy, o.tle[2]]; s.tc = r6(te - x);
    } else if (drag.handle === 'tte') {
      s.tc = r6(Math.max(minC, sx - o.tle[0])); s.tle = [o.tle[0], sy, o.tle[2]];
    }
    if (isPlate) { setPlate(this.m, s); return 'plate'; }
    return 'surface';
  }

  _bind() {
    const c = this.cv;
    let drag = null;
    const pos = (e) => { const b = c.getBoundingClientRect(); return [e.clientX - b.left, e.clientY - b.top]; };
    c.addEventListener('contextmenu', (e) => e.preventDefault());
    c.addEventListener('pointerdown', (e) => {
      const [px, py] = pos(e);
      c.setPointerCapture(e.pointerId);
      if (e.button === 1 || e.button === 2) { drag = { pan: true, px, py, ox: this.view.ox, oy: this.view.oy }; return; }
      const [wx, wy] = this.s2w(px, py);
      if (this.tool === 'draw') {
        if (!this.draw) this.draw = { pts: [] };
        const p = this._drawPoint([this.snapVal(wx), this.snapVal(wy)]);
        if (this.draw.pts.length >= 2) { p[0] = this.snapVal(p[0]); }
        this.draw.pts.push(p);
        if (this.draw.pts.length === 4) {
          const [rle, tle, tte, rte] = this.draw.pts;
          const s = { name: '', rle: [rle[0], rle[1], 0], rc: r6(Math.max(1e-3, rte[0] - rle[0])),
            tle: [tle[0], tle[1], 0], tc: r6(Math.max(1e-3, tte[0] - tle[0])), nchord: 8, nspan: 10 };
          this.draw = null;
          if (Math.abs(tle[1] - rle[1]) < 1e-6) { if (this.opts.onError) this.opts.onError('root and tip must have different y'); }
          else if (this.opts.onAdd) this.opts.onAdd('surface', s);
          this.setTool('select');
        }
        this._hint(); this.render();
        return;
      }
      if (this.tool === 'station') {
        if (this.opts.onAdd) this.opts.onAdd('station', { name: '', p: [this.snapVal(wx), this.snapVal(wy), 0], a: [0, 1, 0] });
        return;
      }
      const h = this.hit(px, py);
      if (h && h.kind !== 'ghost') {
        const same = this.sel && h.kind === this.sel.kind && h.i === this.sel.i;
        if (!same) this.select({ kind: h.kind, i: h.i }, true);
        const target = { kind: h.kind, i: h.i };
        let orig;
        if (h.kind === 'station') orig = JSON.parse(JSON.stringify(this.m.stations[h.i]));
        else if (h.kind === 'control') orig = { ...this.m.controls[h.i] };
        else if (h.kind === 'plate') orig = getPlate(this.m);
        else orig = JSON.parse(JSON.stringify(this.m.surfaces[h.i]));
        drag = { target, handle: h.handle || null, orig, x0: wx, y0: wy, moved: false };
      } else {
        if (h && h.kind === 'ghost') this.select({ kind: 'ghost', i: h.i }, true);
        else if (this.sel) this.select(null, true);
        drag = { pan: true, px, py, ox: this.view.ox, oy: this.view.oy };
      }
    });
    c.addEventListener('pointermove', (e) => {
      const [px, py] = pos(e);
      const [wx, wy] = this.s2w(px, py);
      this.mouse = [wx, wy];
      if (this.opts.onCursor) this.opts.onCursor(wx, wy);
      if (drag && drag.pan) {
        this.view.ox = drag.ox + px - drag.px; this.view.oy = drag.oy + py - drag.py;
        this.render(); return;
      }
      if (drag) {
        drag.moved = true;
        const kind = this._apply(drag, wx, wy);
        this.render();
        if (kind && this.opts.onChange) this.opts.onChange(kind, true);
        return;
      }
      if (this.tool === 'draw') { this.render(); return; }
      if (this.tool === 'select') {
        const h = this.hit(px, py);
        const key = (x) => (x ? `${x.kind}:${x.i}:${x.handle || ''}` : '');
        if (key(h) !== key(this.hover)) {
          this.hover = h;
          c.style.cursor = h ? (h.handle ? 'move' : 'pointer') : 'default';
          this.render();
        }
      }
    });
    c.addEventListener('pointerup', () => {
      if (drag && !drag.pan && drag.moved && this.opts.onChange) this.opts.onChange(drag.target.kind, false);
      drag = null;
    });
    c.addEventListener('pointerleave', () => { this.mouse = null; if (this.tool === 'draw') this.render(); });
    c.addEventListener('dblclick', (e) => {
      const [px, py] = pos(e);
      const h = this.hit(px, py);
      if (h && h.kind === 'ghost' && this.opts.onConvertGhosts) this.opts.onConvertGhosts();
    });
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      const [px, py] = pos(e);
      const [wx, wy] = this.s2w(px, py);
      const f = Math.exp(-e.deltaY * 0.0015);
      this.view.scale = Math.max(1e-2, Math.min(1e6, this.view.scale * f));
      this.view.ox = px - wx * this.view.scale; this.view.oy = py + wy * this.view.scale;
      this.render();
    }, { passive: false });
  }

  cancelDraw() { if (this.draw) { this.draw = null; this._hint(); this.render(); return true; } return false; }
}
