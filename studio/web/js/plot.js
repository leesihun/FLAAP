// Reusable SVG line-plot module: axes, nice ticks (lin/log), multiple series,
// legend with toggles, markers / reference lines, hover crosshair + tooltip,
// linked cursors between plots. No external libraries.
//
//   const p = new Plot(container, {xlabel:'V [m/s]', ylabel:'g', height: 260});
//   p.setSeries([{name, x:[...], y:[...], color, dash, hidden}]);
//   p.setMarks([{x, y, label, color}]); p.setHLines([{y, label}]);
//   p.render();

import { fmt, seriesColor } from './util.js';

const NS = 'http://www.w3.org/2000/svg';
const svgEl = (tag, attrs = {}, parent) => {
  const n = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== undefined && v !== null) n.setAttribute(k, v);
  if (parent) parent.appendChild(n);
  return n;
};
const finite = (v) => typeof v === 'number' && Number.isFinite(v);

/** "Nice" step for a range divided into about n intervals. */
export function niceStep(range, n) {
  if (!(range > 0)) return 1;
  const raw = range / Math.max(1, n);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const f = raw / p;
  return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
}

export function linearTicks(lo, hi, n) {
  const st = niceStep(hi - lo, n);
  const out = [];
  for (let v = Math.ceil(lo / st - 1e-9) * st; v <= hi + st * 1e-9; v += st) out.push(Math.abs(v) < st * 1e-9 ? 0 : v);
  return { ticks: out, step: st };
}

export function logTicks(lo, hi, n) {
  const a = Math.floor(Math.log10(lo)), b = Math.ceil(Math.log10(hi));
  const every = Math.max(1, Math.ceil((b - a) / Math.max(2, n)));
  const out = [];
  for (let e = a; e <= b; e += every) { const v = Math.pow(10, e); if (v >= lo * 0.999 && v <= hi * 1.001) out.push(v); }
  const minor = [];
  if (every === 1 && b - a <= 6) {
    for (let e = a; e < b; e++) for (let m = 2; m < 10; m++) { const v = m * Math.pow(10, e); if (v > lo && v < hi) minor.push(v); }
  }
  return { ticks: out, minor };
}

function tickLabel(v, step, log) {
  if (log) {
    const e = Math.round(Math.log10(v));
    return e >= -2 && e <= 3 ? String(Number(v.toPrecision(3))) : `1e${e}`;
  }
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || (a < 1e-3)) return v.toExponential(1).replace('e+', 'e');
  const dec = Math.max(0, Math.min(8, -Math.floor(Math.log10(step) + 1e-9)));
  return v.toFixed(dec);
}

export class Plot {
  constructor(container, opts = {}) {
    this.box = container;
    this.opts = { height: 260, legend: true, logx: false, logy: false, invertY: false, xlabel: '', ylabel: '',
      limits: {}, legendMax: 14, ...opts };
    this.series = []; this.marks = []; this.hlines = []; this.vlines = [];
    this.svg = svgEl('svg', { class: 'plot' });
    this.box.appendChild(this.svg);
    this.tip = document.createElement('div');
    this.tip.className = 'plottip'; this.tip.hidden = true;
    this.box.appendChild(this.tip);
    this.hoverCbs = [];
    this.ro = new ResizeObserver(() => { if (this.box.clientWidth !== this._lastW) this.render(); });
    this.ro.observe(this.box);
  }

  setSeries(s) { this.series = s.map((x, i) => ({ color: seriesColor(i), ...x })); return this; }
  setMarks(m) { this.marks = m || []; return this; }
  setHLines(h) { this.hlines = h || []; return this; }
  setVLines(v) { this.vlines = v || []; return this; }
  setOptions(o) { Object.assign(this.opts, o); return this; }
  onHover(cb) { this.hoverCbs.push(cb); return this; }
  destroy() { this.ro.disconnect(); this.box.innerHTML = ''; }

  _ok(v, log) { return finite(v) && (!log || v > 0); }

  _range() {
    const { logx, logy, limits } = this.opts;
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (const s of this.series) {
      if (s.hidden) continue;
      const n = Math.min(s.x.length, s.y.length);
      for (let i = 0; i < n; i++) {
        const x = s.x[i], y = s.y[i];
        if (!this._ok(x, logx) || !this._ok(y, logy)) continue;
        if ((finite(limits.xmin) && x < limits.xmin) || (finite(limits.xmax) && x > limits.xmax)) continue;
        if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
      }
    }
    for (const m of this.marks) {
      if (this._ok(m.x, logx) && this._ok(m.y, logy) && m.inRange !== false) {
        x0 = Math.min(x0, m.x); x1 = Math.max(x1, m.x); y0 = Math.min(y0, m.y); y1 = Math.max(y1, m.y);
      }
    }
    for (const h of this.hlines) if (this._ok(h.y, logy) && h.fit !== false) { y0 = Math.min(y0, h.y); y1 = Math.max(y1, h.y); }
    if (!finite(x0)) { x0 = logx ? 1 : 0; x1 = logx ? 10 : 1; }
    if (!finite(y0)) { y0 = logy ? 1 : 0; y1 = logy ? 10 : 1; }
    const pad = (a, b, log) => {
      if (log) { if (a === b) return [a / 2, b * 2]; return [a, b]; }
      if (a === b) { const d = Math.abs(a) * 0.1 || 1; return [a - d, b + d]; }
      const d = (b - a) * 0.04; return [a - d, b + d];
    };
    [x0, x1] = this.opts.padX === false ? (x0 === x1 ? pad(x0, x1, logx) : [x0, x1]) : pad(x0, x1, logx);
    [y0, y1] = pad(y0, y1, logy);
    if (finite(limits.xmin)) x0 = limits.xmin; if (finite(limits.xmax)) x1 = limits.xmax;
    if (finite(limits.ymin)) y0 = limits.ymin; if (finite(limits.ymax)) y1 = limits.ymax;
    if (logx) { x0 = Math.max(x0, 1e-300); x1 = Math.max(x1, x0 * 10); }
    if (logy) { y0 = Math.max(y0, 1e-300); y1 = Math.max(y1, y0 * 1.0001); }
    if (x1 <= x0) x1 = x0 + 1; if (y1 <= y0) y1 = y0 + 1;
    return { x0, x1, y0, y1 };
  }

  render() {
    const o = this.opts;
    const W = Math.max(200, this.box.clientWidth - 2);
    this._lastW = this.box.clientWidth;
    const H = o.height;
    const svg = this.svg;
    while (svg.firstChild) svg.removeChild(svg.firstChild);
    svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
    svg.setAttribute('height', H);

    // legend layout (top, wrapping)
    const legendItems = o.legend && this.series.length > 1 ? this.series.slice(0, o.legendMax) : [];
    const lrows = []; let row = [], rw = 0;
    const maxW = W - 70;
    for (const s of legendItems) {
      const w = 22 + Math.min(180, s.name.length * 6.4) + 14;
      if (rw + w > maxW && row.length) { lrows.push(row); row = []; rw = 0; }
      row.push([s, rw, w]); rw += w;
    }
    if (row.length) lrows.push(row);
    const legendH = lrows.length * 16 + (lrows.length ? 4 : 0);

    const r = this._range();
    const fx = o.logx ? Math.log10 : (v) => v, fy = o.logy ? Math.log10 : (v) => v;
    const yt = o.logy ? logTicks(r.y0, r.y1, 5) : linearTicks(r.y0, r.y1, Math.max(3, Math.round((H - legendH) / 50)));
    const ylabels = yt.ticks.map((v) => tickLabel(v, yt.step, o.logy));
    const ml = 18 + Math.max(24, ...ylabels.map((s) => s.length * 6.3)) + (o.ylabel ? 14 : 0);
    const m = { l: ml, r: 14, t: 8 + legendH, b: 36 };
    const pw = W - m.l - m.r, ph = Math.max(30, H - m.t - m.b);
    const X0 = fx(r.x0), X1 = fx(r.x1), Y0 = fy(r.y0), Y1 = fy(r.y1);
    const sx = (v) => m.l + (fx(v) - X0) / (X1 - X0) * pw;
    const sy = o.invertY ? (v) => m.t + (fy(v) - Y0) / (Y1 - Y0) * ph : (v) => m.t + ph - (fy(v) - Y0) / (Y1 - Y0) * ph;
    this._g = { m, pw, ph, sx, sy, r, fx, fy, X0, X1, W, H };
    const xt = o.logx ? logTicks(r.x0, r.x1, Math.max(3, Math.round(pw / 70))) : linearTicks(r.x0, r.x1, Math.max(3, Math.round(pw / 85)));

    const clipId = 'c' + Math.random().toString(36).slice(2, 9);
    const defs = svgEl('defs', {}, svg);
    const cp = svgEl('clipPath', { id: clipId }, defs);
    svgEl('rect', { x: m.l, y: m.t, width: pw, height: ph }, cp);

    // grid
    const grid = svgEl('g', { class: 'grid' }, svg);
    for (const v of (yt.minor || [])) svgEl('line', { x1: m.l, x2: m.l + pw, y1: sy(v), y2: sy(v), opacity: 0.45 }, grid);
    for (const v of (xt.minor || [])) svgEl('line', { y1: m.t, y2: m.t + ph, x1: sx(v), x2: sx(v), opacity: 0.45 }, grid);
    for (const v of yt.ticks) svgEl('line', { x1: m.l, x2: m.l + pw, y1: sy(v), y2: sy(v) }, grid);
    for (const v of xt.ticks) svgEl('line', { y1: m.t, y2: m.t + ph, x1: sx(v), x2: sx(v) }, grid);
    // zero line
    if (!o.logy && r.y0 < 0 && r.y1 > 0) svgEl('line', { class: 'zero', x1: m.l, x2: m.l + pw, y1: sy(0), y2: sy(0) }, svg);

    // axes + tick labels
    const ax = svgEl('g', { class: 'axis' }, svg);
    svgEl('line', { x1: m.l, x2: m.l + pw, y1: m.t + ph, y2: m.t + ph }, ax);
    svgEl('line', { x1: m.l, x2: m.l, y1: m.t, y2: m.t + ph }, ax);
    yt.ticks.forEach((v, i) => {
      svgEl('text', { x: m.l - 6, y: sy(v) + 3.5, 'text-anchor': 'end' }, ax).textContent = ylabels[i];
    });
    let lastX = -1e9;
    for (const v of xt.ticks) {
      const x = sx(v), lab = tickLabel(v, xt.step, o.logx);
      if (x - lastX < lab.length * 6.5 + 6) continue;
      lastX = x;
      svgEl('text', { x, y: m.t + ph + 15, 'text-anchor': 'middle' }, ax).textContent = lab;
    }
    if (o.xlabel) svgEl('text', { class: 'xlabel', x: m.l + pw / 2, y: H - 4, 'text-anchor': 'middle' }, svg).textContent = o.xlabel;
    if (o.ylabel) {
      svgEl('text', { class: 'ylabel', x: 0, y: 0, 'text-anchor': 'middle',
        transform: `translate(12 ${m.t + ph / 2}) rotate(-90)` }, svg).textContent = o.ylabel;
    }

    // reference lines
    const ref = svgEl('g', { 'clip-path': `url(#${clipId})` }, svg);
    for (const h of this.hlines) {
      if (!this._ok(h.y, o.logy)) continue;
      svgEl('line', { x1: m.l, x2: m.l + pw, y1: sy(h.y), y2: sy(h.y), stroke: h.color || '#4a5753',
        'stroke-dasharray': h.dash || '5 4', 'stroke-width': 1.2 }, ref);
      if (h.label) svgEl('text', { x: m.l + pw - 4, y: sy(h.y) - 4, 'text-anchor': 'end' }, ref).textContent = h.label;
    }
    for (const v of this.vlines) {
      if (!this._ok(v.x, o.logx)) continue;
      svgEl('line', { y1: m.t, y2: m.t + ph, x1: sx(v.x), x2: sx(v.x), stroke: v.color || '#4a5753',
        'stroke-dasharray': v.dash || '5 4', 'stroke-width': 1.2 }, ref);
      if (v.label) svgEl('text', { x: sx(v.x) + 4, y: m.t + 12 }, ref).textContent = v.label;
    }

    // series
    const sg = svgEl('g', { class: 'series', 'clip-path': `url(#${clipId})` }, svg);
    for (const s of this.series) {
      if (s.hidden) continue;
      let d = '', pen = false;
      const n = Math.min(s.x.length, s.y.length);
      for (let i = 0; i < n; i++) {
        const x = s.x[i], y = s.y[i];
        if (!this._ok(x, o.logx) || !this._ok(y, o.logy)) { pen = false; continue; }
        d += (pen ? 'L' : 'M') + sx(x).toFixed(1) + ' ' + sy(y).toFixed(1);
        pen = true;
      }
      svgEl('path', { d, stroke: s.color, 'stroke-dasharray': s.dash, 'stroke-width': s.width || 2 }, sg);
      if (s.points) {
        for (let i = 0; i < n; i++) {
          if (this._ok(s.x[i], o.logx) && this._ok(s.y[i], o.logy)) svgEl('circle', { cx: sx(s.x[i]), cy: sy(s.y[i]), r: 2.5, fill: s.color }, sg);
        }
      }
    }

    // marks
    const mg = svgEl('g', { class: 'mark' }, svg);
    for (const mk of this.marks) {
      if (!this._ok(mk.x, o.logx) || !this._ok(mk.y, o.logy)) continue;
      const x = sx(mk.x), y = sy(mk.y);
      if (x < m.l - 1 || x > m.l + pw + 1 || y < m.t - 1 || y > m.t + ph + 1) continue;
      const sz = 6;
      svgEl('path', { d: `M${x} ${y - sz}L${x + sz} ${y}L${x} ${y + sz}L${x - sz} ${y}Z`, fill: mk.color || '#b84d48',
        stroke: '#fbfaf6', 'stroke-width': 2 }, mg);
      if (mk.label) {
        const right = x < m.l + pw - 120;
        svgEl('text', { x: right ? x + 9 : x - 9, y: y - 8, 'text-anchor': right ? 'start' : 'end' }, mg).textContent = mk.label;
      }
    }

    // legend
    const lg = svgEl('g', { class: 'legend' }, svg);
    lrows.forEach((rowItems, ri) => {
      for (const [s, x0, w] of rowItems) {
        const g = svgEl('g', { transform: `translate(${m.l + x0} ${8 + ri * 16})`, class: s.hidden ? 'off' : '' }, lg);
        svgEl('rect', { x: -2, y: -1, width: w - 6, height: 15, fill: 'transparent' }, g);
        svgEl('line', { x1: 0, x2: 16, y1: 7, y2: 7, stroke: s.color, 'stroke-width': 3, 'stroke-dasharray': s.dash,
          'stroke-linecap': 'round', opacity: s.hidden ? 0.3 : 1 }, g);
        const t = svgEl('text', { x: 21, y: 11 }, g);
        t.textContent = s.name.length > 28 ? s.name.slice(0, 27) + '…' : s.name;
        g.addEventListener('click', () => { s.hidden = !s.hidden; this.render(); if (this.opts.onToggle) this.opts.onToggle(s); });
        const tt = svgEl('title', {}, g); tt.textContent = s.name + ' (click to toggle)';
      }
    });

    // hover layer
    this.cross = svgEl('g', { class: 'cross-g' }, svg);
    const ov = svgEl('rect', { x: m.l, y: m.t, width: pw, height: ph, fill: 'transparent' }, svg);
    ov.addEventListener('mousemove', (ev) => {
      const b = svg.getBoundingClientRect();
      const px = (ev.clientX - b.left) * (W / b.width);
      const xd = this._invX(px);
      this._drawCursor(xd, true, ev);
      for (const cb of this.hoverCbs) cb(xd);
    });
    ov.addEventListener('mouseleave', () => { this._clearCursor(); for (const cb of this.hoverCbs) cb(null); });
  }

  _invX(px) {
    const g = this._g;
    const f = g.X0 + (px - g.m.l) / g.pw * (g.X1 - g.X0);
    return this.opts.logx ? Math.pow(10, f) : f;
  }

  _clearCursor() { if (this.cross) while (this.cross.firstChild) this.cross.removeChild(this.cross.firstChild); this.tip.hidden = true; }

  /** Show the crosshair at data x (linked cursor); null hides it. */
  showCursor(x) { if (x === null || x === undefined) this._clearCursor(); else this._drawCursor(x, false); }

  _drawCursor(xd, withTip, ev) {
    const g = this._g; if (!g) return;
    this._clearCursor();
    const o = this.opts;
    const px = g.sx(xd);
    if (!(px >= g.m.l && px <= g.m.l + g.pw)) return;
    svgEl('line', { class: 'cross', x1: px, x2: px, y1: g.m.t, y2: g.m.t + g.ph }, this.cross);
    const rows = [];
    for (const s of this.series) {
      if (s.hidden) continue;
      let best = -1, bd = Infinity;
      const n = Math.min(s.x.length, s.y.length);
      for (let i = 0; i < n; i++) {
        if (!this._ok(s.x[i], o.logx) || !this._ok(s.y[i], o.logy)) continue;
        const d = Math.abs(g.sx(s.x[i]) - px);
        if (d < bd) { bd = d; best = i; }
      }
      if (best < 0 || bd > 30) continue;
      const cy = g.sy(s.y[best]);
      if (cy >= g.m.t - 2 && cy <= g.m.t + g.ph + 2) {
        svgEl('circle', { cx: g.sx(s.x[best]), cy, r: 4, fill: s.color, stroke: '#fbfaf6', 'stroke-width': 2 }, this.cross);
      }
      rows.push([s, s.x[best], s.y[best]]);
    }
    if (!withTip || !rows.length) { this.tip.hidden = true; return; }
    const xl = o.xlabel ? o.xlabel.split(' [')[0] : 'x';
    let html = `<b>${xl} = ${fmt(xd, 5)}</b>`;
    for (const [s, x, y] of rows.slice(0, 14)) {
      const xs = Math.abs(x - xd) > Math.abs(xd) * 0.01 + 1e-12 ? ` @${fmt(x, 4)}` : '';
      html += `<br><span class="k" style="background:${s.color}"></span>${s.name.slice(0, 24)}: ${fmt(y, 5)}${xs}`;
    }
    if (rows.length > 14) html += `<br>… ${rows.length - 14} more`;
    this.tip.innerHTML = html;
    this.tip.hidden = false;
    const bb = this.box.getBoundingClientRect();
    const sb = this.svg.getBoundingClientRect();
    const scale = sb.width / g.W;
    let left = sb.left - bb.left + px * scale + 14;
    const tw = this.tip.offsetWidth;
    if (left + tw > bb.width - 4) left = sb.left - bb.left + px * scale - tw - 14;
    const top = ev ? Math.max(4, Math.min(ev.clientY - bb.top - 20, bb.height - this.tip.offsetHeight - 4)) : 10;
    this.tip.style.left = `${Math.max(2, left)}px`;
    this.tip.style.top = `${top}px`;
  }
}

/** Link the crosshairs of several plots (same x quantity). */
export function linkPlots(plots) {
  for (const p of plots) p.onHover((x) => { for (const q of plots) if (q !== p) q.showCursor(x); });
}
