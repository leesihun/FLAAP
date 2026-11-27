// Diverging colour map (cool blue - neutral - warm red) and a DOM colour bar.

const STOPS = [
  [-1.0, [33, 80, 130]],
  [-0.5, [111, 158, 196]],
  [0.0, [236, 233, 226]],
  [0.5, [219, 140, 118]],
  [1.0, [160, 44, 40]],
];
const SEQ = [ // sequential single-hue (light -> dark green)
  [0.0, [233, 243, 226]], [0.5, [104, 170, 140]], [1.0, [13, 74, 60]],
];

function interp(stops, t) {
  if (!Number.isFinite(t)) return [200, 200, 200];
  if (t <= stops[0][0]) return stops[0][1];
  for (let i = 1; i < stops.length; i++) {
    if (t <= stops[i][0]) {
      const [a, ca] = stops[i - 1], [b, cb] = stops[i];
      const f = (t - a) / (b - a);
      return [0, 1, 2].map((k) => Math.round(ca[k] + f * (cb[k] - ca[k])));
    }
  }
  return stops[stops.length - 1][1];
}

/** t in [-1, 1] -> 'rgb(...)'. */
export function diverging(t) { const c = interp(STOPS, t); return `rgb(${c[0]},${c[1]},${c[2]})`; }
/** t in [0, 1] -> 'rgb(...)'. */
export function sequential(t) { const c = interp(SEQ, t); return `rgb(${c[0]},${c[1]},${c[2]})`; }

/** Colour bar element: {el, set(lo, hi, label, kind)} */
export function colorbar() {
  const el = document.createElement('div');
  el.className = 'cbar';
  const lbl = document.createElement('div'); lbl.className = 'lbl';
  const cv = document.createElement('canvas'); cv.width = 256; cv.height = 12;
  const ticks = document.createElement('div'); ticks.className = 'ticks';
  el.append(lbl, cv, ticks);
  const set = (lo, hi, label = '', kind = 'div', fmtf = (v) => Number(v.toPrecision(3))) => {
    const ctx = cv.getContext('2d');
    for (let i = 0; i < 256; i++) {
      const t = i / 255;
      ctx.fillStyle = kind === 'div' ? diverging(2 * t - 1) : sequential(t);
      ctx.fillRect(i, 0, 1, 12);
    }
    lbl.textContent = label;
    const mid = 0.5 * (lo + hi);
    ticks.innerHTML = '';
    for (const v of [lo, mid, hi]) { const s = document.createElement('span'); s.textContent = fmtf(v); ticks.appendChild(s); }
  };
  return { el, set };
}
