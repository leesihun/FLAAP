// Minimal 3-D viewer on a 2-D canvas: orthographic / perspective camera,
// orbit (left drag), pan (right or shift drag), zoom (wheel), painter's
// algorithm for faces, lines and points, face picking, axis triad.
//
//   const v = new View3D(container, {onHover: (face) => ...});
//   v.setScene({faces:[{pts:[[x,y,z],...], fill, stroke, data}], lines:[{a,b,color,width}], points:[{p,color,r}]});
//   v.fit(); v.render();

import { fitCanvas } from './util.js';

const D2R = Math.PI / 180;
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };

export const VIEWS = { iso: [-60, 30], top: [-90, 89.9], front: [-90, 0.01], side: [180, 0.01], rear: [0, 25] };

export class View3D {
  constructor(container, opts = {}) {
    this.box = container;
    this.box.classList.add('view3d');
    this.canvas = document.createElement('canvas');
    this.overlay = document.createElement('div'); this.overlay.className = 'overlay';
    this.readout = document.createElement('div'); this.readout.className = 'readout';
    this.box.append(this.canvas, this.overlay, this.readout);
    this.opts = { perspective: false, background: null, triad: true, ...opts };
    this.scene = { faces: [], lines: [], points: [] };
    this.target = [0, 0, 0]; this.radius = 1;
    [this.yaw, this.pitch] = VIEWS[opts.view || 'iso'];
    this.zoom = 1;
    this._bind();
    this.ro = new ResizeObserver(() => this.render());
    this.ro.observe(this.box);
  }

  destroy() { this.ro.disconnect(); this.box.innerHTML = ''; this.box.classList.remove('view3d'); }

  setScene(s) { this.scene = { faces: [], lines: [], points: [], ...s }; return this; }

  /** Fit camera target/radius to the given points (or the current scene). */
  fit(pts) {
    const all = pts || [...this.scene.faces.flatMap((f) => f.pts), ...this.scene.lines.flatMap((l) => [l.a, l.b]),
      ...this.scene.points.map((p) => p.p)];
    if (!all.length) return this;
    const lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (const p of all) for (let k = 0; k < 3; k++) { if (p[k] < lo[k]) lo[k] = p[k]; if (p[k] > hi[k]) hi[k] = p[k]; }
    this.target = [0, 1, 2].map((k) => 0.5 * (lo[k] + hi[k]));
    this.radius = Math.max(1e-9, 0.5 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]));
    this.zoom = 1;
    return this;
  }

  setView(name) { [this.yaw, this.pitch] = VIEWS[name] || VIEWS.iso; this.zoom = 1; this.render(); }

  _camera(w, h) {
    const y = this.yaw * D2R, p = this.pitch * D2R;
    const d = [Math.cos(p) * Math.cos(y), Math.cos(p) * Math.sin(y), Math.sin(p)];
    const f = [-d[0], -d[1], -d[2]];
    let r = cross(f, [0, 0, 1]);
    if (Math.hypot(...r) < 1e-6) r = [1, 0, 0];
    r = norm(r);
    const u = cross(r, f);
    const scale = Math.min(w, h) * 0.46 / this.radius * this.zoom;
    const dist = this.radius * 4;
    const persp = this.opts.perspective;
    const cx = w / 2, cy = h / 2, t = this.target;
    return {
      d, r, u, scale,
      proj: (P) => {
        const v = sub(P, t);
        const xr = dot(v, r), yu = dot(v, u), zd = dot(v, d);
        const k = persp ? dist / Math.max(dist - zd, dist * 0.05) : 1;
        return [cx + xr * scale * k, cy - yu * scale * k, -zd];
      },
    };
  }

  render() {
    const { w, h, ctx } = fitCanvas(this.canvas);
    this._w = w; this._h = h;
    ctx.clearRect(0, 0, w, h);
    if (this.opts.background) { ctx.fillStyle = this.opts.background; ctx.fillRect(0, 0, w, h); }
    const cam = this._camera(w, h);
    this._cam = cam;
    const items = [];
    const light = norm([0.3, -0.4, 0.85]);
    for (const f of this.scene.faces) {
      const sp = f.pts.map(cam.proj);
      let z = 0; for (const q of sp) z += q[2];
      let fill = f.fill;
      if (f.shade && f.pts.length >= 3) {
        const n = norm(cross(sub(f.pts[2], f.pts[0]), sub(f.pts[f.pts.length - 1], f.pts[1])));
        fill = shadeColor(f.fill, 0.72 + 0.28 * Math.abs(dot(n, light)));
      }
      items.push({ t: 0, z: z / sp.length, sp, f, fill });
    }
    for (const l of this.scene.lines) {
      const a = cam.proj(l.a), b = cam.proj(l.b);
      items.push({ t: 1, z: 0.5 * (a[2] + b[2]) - (l.bias || 0) * this.radius, a, b, l });
    }
    for (const p of this.scene.points) {
      const q = cam.proj(p.p);
      items.push({ t: 2, z: q[2] - this.radius * 0.01, q, p });
    }
    items.sort((A, B) => B.z - A.z);
    this._faceItems = items.filter((i) => i.t === 0);
    ctx.lineJoin = 'round';
    for (const it of items) {
      if (it.t === 0) {
        ctx.beginPath();
        it.sp.forEach((q, i) => (i ? ctx.lineTo(q[0], q[1]) : ctx.moveTo(q[0], q[1])));
        ctx.closePath();
        if (it.fill) { ctx.globalAlpha = it.f.alpha ?? 1; ctx.fillStyle = it.fill; ctx.fill(); ctx.globalAlpha = 1; }
        if (it.f.stroke) { ctx.strokeStyle = it.f.stroke; ctx.lineWidth = it.f.lineWidth || 0.6; ctx.stroke(); }
        if (it === this._hoverItem) { ctx.strokeStyle = '#17211e'; ctx.lineWidth = 2; ctx.stroke(); }
      } else if (it.t === 1) {
        ctx.beginPath(); ctx.moveTo(it.a[0], it.a[1]); ctx.lineTo(it.b[0], it.b[1]);
        ctx.strokeStyle = it.l.color || '#17211e'; ctx.lineWidth = it.l.width || 1;
        ctx.setLineDash(it.l.dash || []); ctx.stroke(); ctx.setLineDash([]);
      } else {
        ctx.beginPath(); ctx.arc(it.q[0], it.q[1], it.p.r || 2.5, 0, 2 * Math.PI);
        ctx.fillStyle = it.p.color || '#17211e'; ctx.fill();
      }
    }
    if (this.opts.triad) this._triad(ctx, cam, h);
  }

  _triad(ctx, cam, h) {
    const o = [44, h - 44], L = 26;
    const axes = [[[1, 0, 0], 'x', '#b84d48'], [[0, 1, 0], 'y', '#126b57'], [[0, 0, 1], 'z', '#3f6f91']];
    ctx.font = '600 11px Inter, system-ui, sans-serif';
    for (const [a, n, c] of axes) {
      const ex = o[0] + dot(a, cam.r) * L, ey = o[1] - dot(a, cam.u) * L;
      ctx.beginPath(); ctx.moveTo(o[0], o[1]); ctx.lineTo(ex, ey);
      ctx.strokeStyle = c; ctx.lineWidth = 2; ctx.stroke();
      ctx.fillStyle = c; ctx.fillText(n, ex + (ex - o[0]) * 0.25 - 3, ey + (ey - o[1]) * 0.25 + 4);
    }
  }

  _pick(mx, my) {
    const items = this._faceItems || [];
    for (let i = items.length - 1; i >= 0; i--) {
      const sp = items[i].sp;
      let inside = false;
      for (let a = 0, b = sp.length - 1; a < sp.length; b = a++) {
        const [xa, ya] = sp[a], [xb, yb] = sp[b];
        if ((ya > my) !== (yb > my) && mx < (xb - xa) * (my - ya) / (yb - ya) + xa) inside = !inside;
      }
      if (inside) return items[i];
    }
    return null;
  }

  _bind() {
    const c = this.canvas;
    let drag = null;
    c.addEventListener('contextmenu', (e) => e.preventDefault());
    c.addEventListener('pointerdown', (e) => {
      c.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, pan: e.button === 2 || e.button === 1 || e.shiftKey };
      c.style.cursor = 'grabbing';
    });
    c.addEventListener('pointerup', () => { drag = null; c.style.cursor = ''; });
    c.addEventListener('pointermove', (e) => {
      if (drag) {
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        drag.x = e.clientX; drag.y = e.clientY;
        if (drag.pan && this._cam) {
          const s = this._cam.scale;
          this.target = [0, 1, 2].map((k) => this.target[k] - (dx * this._cam.r[k] - dy * this._cam.u[k]) / s);
        } else {
          this.yaw -= dx * 0.45;
          this.pitch = Math.max(-89.9, Math.min(89.9, this.pitch + dy * 0.45));
        }
        this.render();
        return;
      }
      if (this.opts.onHover) {
        const b = c.getBoundingClientRect();
        const it = this._pick(e.clientX - b.left, e.clientY - b.top);
        if (it !== this._hoverItem) {
          this._hoverItem = it;
          this.readout.textContent = it ? (this.opts.onHover(it.f) || '') : '';
          this.render();
        }
      }
    });
    c.addEventListener('pointerleave', () => {
      if (this._hoverItem) { this._hoverItem = null; this.readout.textContent = ''; this.render(); }
    });
    c.addEventListener('wheel', (e) => {
      e.preventDefault();
      this.zoom *= Math.exp(-e.deltaY * 0.0015);
      this.zoom = Math.max(0.05, Math.min(80, this.zoom));
      this.render();
    }, { passive: false });
    c.addEventListener('dblclick', () => { this.fit(); this.render(); });
  }
}

function shadeColor(col, k) {
  let r, g, b;
  if (col.startsWith('#')) {
    const n = parseInt(col.slice(1), 16);
    r = (n >> 16) & 255; g = (n >> 8) & 255; b = n & 255;
  } else {
    const m = col.match(/\d+/g); if (!m) return col;
    [r, g, b] = m.map(Number);
  }
  return `rgb(${Math.round(r * k)},${Math.round(g * k)},${Math.round(b * k)})`;
}
