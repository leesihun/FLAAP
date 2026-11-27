"""Supersonic lifting-surface aerodynamics: piston theory and the Mach box method.

Both methods return a :class:`~flaap.aero.gaf.GAFTable` in exactly the
convention of :meth:`flaap.aero.gaf.AeroModel.table`, so the result can be
passed unchanged to :func:`flaap.flutter.pk_method`, :func:`k_method`, the
Roger RFA and the time-domain / gust routines::

    W(k)   = dh/dx + i (k / b) h                (normalwash per unit U)
    dCp    = Cp_lower - Cp_upper                (positive = lift)
    Q(k)   = H^T diag(A) dCp(k)                 (generalized force / q)
    Qg(k)  : gust column for w_g = exp(i w (t - (x - x0)/U)) per unit W_g/U

Launcher entry point (``aero_method piston | machbox``)::

    supersonic_gaf_table(mesh, modes, mach, bref, kvals, method="machbox",
                         gust=True, keep_cp=True, x_gust0=None, **options)

* ``mesh``  - the DLM :class:`~flaap.geometry.AeroMesh` (planar surfaces;
  ``mesh.symmetry`` decides the image about ``y = 0``)
* ``modes`` - :class:`~flaap.aero.gaf.AeroModes` on that mesh (piston theory
  uses them box by box; the Mach box method interpolates them onto its own
  box grid with :func:`aeromodes_field`), or for ``machbox`` a callable
  ``f(points (P, 3)) -> (h (P, nm), hx (P, nm))``
* options: piston - ``order`` (1 or 3), ``van_dyke`` (bool), ``thickness_slope``,
  ``gamma``, ``slope``; machbox - ``nx`` (boxes along the overall chord), ``n_gauss``.

Piston theory
-------------
Local "piston" relation of Lighthill (1953) / Ashley & Zartarian (1956).
First order (``van_dyke=False``)::

    dCp = -(4 / M) [hx + i (k/b) h]

Van Dyke's supersonic quasi-steady form (default; the first two terms of
the small-k expansion of the exact 2-D theory, correct as M -> 1)::

    dCp = -(4 / beta) [hx + (M^2 - 2)/(M^2 - 1) i (k/b) h],  beta = sqrt(M^2 - 1)

Third order with a mean thickness slope ``tau' = d(thickness)/dx`` per box::

    dCp = -(4 / M) [1 + (gamma+1)/2 v0 + (gamma+1)/4 v0^2] [hx + i (k/b) h],
    v0 = M tau' / 2

Mach box method
---------------
Upper-surface potential of a planar thin wing in linearised supersonic flow
(Pines, Dugundji & Neuringer 1955; Zartarian & Hsu 1955)::

    phi(x, y) = -(1/pi) iint w(xi, eta) exp(-i lam (x - xi)) cos(kap R) / R dxi deta
    R = sqrt((x - xi)^2 - beta^2 (y - eta)^2),  lam = w M^2/(U beta^2),
    kap = w M/(U beta^2)

over the forward Mach cone.  The plane is divided into boxes of
``dx x dy`` with ``dx = beta dy`` (square in the scaled coordinates
``(x, beta y)``, diagonals along the Mach lines).  Each box carries a
constant normalwash; the box integrals are evaluated with the substitution
``beta (y - eta) = (x - xi) sin(theta)`` (``dY / R = d theta``), which
removes the Mach-line singularity, and Gauss quadrature.  The self box gives
``pi dx / 2`` in the steady case.  Off the wing, in the plane ``z = 0``:

* diaphragm boxes (beside the wing): ``phi = 0`` (continuous potential,
  no pressure jump) - the normalwash there is unknown and solved for,
* wake boxes (aft of the trailing edge, inside the span): ``dp = 0`` ->
  ``phi = phi_TE exp(-i w (x - x_TE) / U)``.

Because the influence of a box only reaches downstream rows, the solution is
an explicit row-by-row march.  Pressures follow from
``dCp = 4 (dphi/dx + i (k/b) phi) / U``, evaluated per box with the potential at
the box edges (``phi = 0`` at the leading edge, ``phi`` interpolated at the
trailing edge), so that the chordwise integral of ``dCp`` is exact for a
piecewise linear potential.

Validated (``tests/test_supersonic.py``, ``validation/supersonic_validation.py``)
against Ackeret 2-D theory, the rectangular- and delta-wing lift slopes of
linear supersonic theory, the Garrick-Rubinow 2-D oscillating-airfoil
integral (:func:`garrick_rubinow_2d`), the piston-theory panel-flutter
boundaries (Hedgepeth / Dowell) and the AGARD 445.6 supersonic flutter points.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..geometry import AeroMesh
from .gaf import AeroModes, GAFTable

GAMMA_AIR = 1.4

__all__ = ["piston_theory_table", "aeromodes_field", "MachBox", "machbox_gaf_table",
           "supersonic_gaf_table", "garrick_rubinow_2d"]


def _beta(mach: float) -> float:
    if mach <= 1.0:
        raise ValueError(f"supersonic aerodynamics need mach > 1 (got {mach})")
    return float(np.sqrt(mach * mach - 1.0))


def _gust_w(x, gamma, k, bref, x0):
    return -np.cos(gamma) * np.exp(-1j * (k / bref) * (x - x0))


# =====================================================================
# piston theory
# =====================================================================
def piston_theory_table(mesh: AeroMesh, modes: AeroModes, mach: float, bref: float, kvals,
                        order: int = 1, van_dyke: bool = True, thickness_slope=None,
                        gust: bool = True, keep_cp: bool = True, x_gust0: float | None = None,
                        gamma: float = GAMMA_AIR, slope: str = "central") -> GAFTable:
    """GAF table from (first- or third-order) piston theory on the DLM boxes.

    The box displacement is taken at the box centre (mean of the 1/4- and
    3/4-chord values).  The slope is evaluated at the box centre too
    (``slope='central'``: ``(h_c - h_f) / (x_c - x_f)``, exact for a quadratic
    displacement); boxes where a mode has no displacement (``alpha`` and rate
    modes, ``h = 0``) fall back to the 3/4-chord slope ``hx_c``, which is also
    used everywhere with ``slope='colloc'`` (first-order accurate in the box
    length).  ``thickness_slope`` (scalar or per box) is the streamwise
    derivative of the full airfoil thickness, used only with ``order=3``.
    """
    beta = _beta(mach)
    kvals = np.atleast_1d(np.asarray(kvals, float))
    h = 0.5 * (modes.h_c + modes.h_f)
    hx = modes.hx_c
    if slope == "central":
        dxs = (mesh.colloc[:, 0] - mesh.force_point[:, 0])[:, None]
        has_h = (modes.h_c != 0.0) | (modes.h_f != 0.0)
        hx = np.where(has_h & (dxs > 0), (modes.h_c - modes.h_f) / np.where(dxs > 0, dxs, 1.0), hx)
    elif slope != "colloc":
        raise ValueError("slope must be 'central' or 'colloc'")
    if order == 1:
        if van_dyke:
            a0 = np.full(mesh.n, 4.0 / beta)
            a1 = a0 * (mach * mach - 2.0) / (mach * mach - 1.0)
        else:
            a0 = a1 = np.full(mesh.n, 4.0 / mach)
    elif order == 3:
        ts = np.zeros(mesh.n) if thickness_slope is None else np.broadcast_to(
            np.asarray(thickness_slope, float), (mesh.n,))
        v0 = 0.5 * mach * ts
        a0 = a1 = (4.0 / mach) * (1.0 + 0.5 * (gamma + 1) * v0 + 0.25 * (gamma + 1) * v0 * v0)
    else:
        raise ValueError("piston theory order must be 1 or 3")
    x0 = float(mesh.corners[:, :, 0].min()) if x_gust0 is None else x_gust0
    HA = h * mesh.area[:, None]
    nk, nm = kvals.size, modes.nm
    Q = np.zeros((nk, nm, nm), complex)
    Qg = np.zeros((nk, nm), complex) if gust else None
    cp = np.zeros((nk, mesh.n, nm), complex) if keep_cp else None
    cpg = np.zeros((nk, mesh.n), complex) if (keep_cp and gust) else None
    for i, k in enumerate(kvals):
        c = -(a0[:, None] * hx + 1j * (k / bref) * a1[:, None] * h)
        Q[i] = HA.T @ c
        if gust:
            cg = -a0 * _gust_w(mesh.center[:, 0], mesh.gamma, k, bref, x0)
            Qg[i] = HA.T @ cg
            if keep_cp:
                cpg[i] = cg
        if keep_cp:
            cp[i] = c
    return GAFTable(kvals, Q, Qg, float(mach), float(bref), list(modes.names), cp, cpg)


# =====================================================================
# mode interpolation onto other points
# =====================================================================
def aeromodes_field(mesh: AeroMesh, modes: AeroModes) -> Callable:
    """Return ``f(points) -> (h, hx)`` interpolating ``modes`` off the DLM mesh.

    ``h`` is an infinite-plate spline through the 1/4- and 3/4-chord box
    values; ``hx`` is its x-derivative plus a spline of the slope residual at
    the 3/4-chord points, so modes whose slope is not the derivative of the
    displacement (``alpha``, rate modes) are reproduced exactly.
    """
    from ..structures.spline import ips_matrices
    xd = np.vstack([mesh.force_point[:, :2], mesh.colloc[:, :2]])
    vals = np.vstack([modes.h_f, modes.h_c])
    _, Gx_c = ips_matrices(xd, mesh.colloc[:, :2])
    resid = modes.hx_c - Gx_c @ vals

    def f(points):
        p = np.asarray(points, float)[:, :2]
        G, Gx = ips_matrices(xd, p)
        G2, _ = ips_matrices(mesh.colloc[:, :2], p)
        return G @ vals, Gx @ vals + G2 @ resid

    return f


# =====================================================================
# Mach box method
# =====================================================================
_WING, _DIAPHRAGM, _WAKE = 0, 1, 2


@dataclass
class _Trap:
    x_le0: float
    c0: float
    x_le1: float
    c1: float
    y0: float
    y1: float

    def le(self, y):
        t = (y - self.y0) / (self.y1 - self.y0)
        return self.x_le0 + t * (self.x_le1 - self.x_le0)

    def te(self, y):
        t = (y - self.y0) / (self.y1 - self.y0)
        return self.le(y) + self.c0 + t * (self.c1 - self.c0)


def _trapezoids(mesh_or_surfaces):
    surfs = mesh_or_surfaces.surfaces if isinstance(mesh_or_surfaces, AeroMesh) else list(mesh_or_surfaces)
    out = []
    for s in surfs:
        r, t = np.asarray(s.root_le, float), np.asarray(s.tip_le, float)
        if abs(t[2] - r[2]) > 1e-6 * max(s.span, 1e-9):
            raise ValueError(f"Mach box method needs planar surfaces (surface {s.name} has dihedral)")
        if t[1] >= r[1]:
            out.append(_Trap(r[0], s.root_chord, t[0], s.tip_chord, r[1], t[1]))
        else:
            out.append(_Trap(t[0], s.tip_chord, r[0], s.root_chord, t[1], r[1]))
    return out


def _gauss(n):
    t, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (t + 1.0), 0.5 * w


class MachBox:
    """Mach box grid for the planar surfaces of an aero mesh.

    Parameters
    ----------
    mesh : AeroMesh or list of Surface (planar, ``z`` constant)
    mach : free-stream Mach number (> 1)
    nx : approximate number of boxes along the overall streamwise extent;
        the box size is adjusted so that the planform span is an integer
        number of box widths ``dy = dx / beta``
    symmetry : 'symmetric' | 'antisymmetric' | 'none' (default: ``mesh.symmetry``)
    n_gauss : Gauss points per integration interval for the box integrals
    """

    def __init__(self, mesh, mach: float, nx: int = 60, symmetry: str | None = None,
                 n_gauss: int = 10):
        self.mach = float(mach)
        self.beta = _beta(mach)
        self.traps = _trapezoids(mesh)
        if symmetry is None:
            symmetry = mesh.symmetry if isinstance(mesh, AeroMesh) else "none"
        self.symmetry = symmetry
        self.sign = {"symmetric": 1.0, "antisymmetric": -1.0}.get(symmetry, 0.0)
        self.n_gauss = int(n_gauss)
        b = self.beta
        x_min = min(t.x_le0 if t.x_le0 < t.x_le1 else t.x_le1 for t in self.traps)
        x_max = max(max(t.x_le0 + t.c0, t.x_le1 + t.c1) for t in self.traps)
        y_lo = min(t.y0 for t in self.traps)
        y_hi = max(t.y1 for t in self.traps)
        dx0 = (x_max - x_min) / nx
        span = y_hi - y_lo
        nsp = max(1, int(round(span * b / dx0)))
        self.dy = span / nsp
        self.dx = b * self.dy
        ext = (x_max - x_min) / (2 * b) + 2 * self.dy        # outboard diaphragm reach
        if self.sign != 0.0:
            if y_lo < -1e-9 * span:
                raise ValueError("symmetric Mach box model needs the planform at y >= 0")
            y_start = 0.0
        else:
            y_start = y_lo - self.dy * np.ceil(ext / self.dy)
        y_end = y_hi + ext
        self.ny = int(np.ceil((y_end - y_start) / self.dy))
        self.nrow = int(np.ceil((x_max - x_min) / self.dx)) + 2
        self.x0, self.y0 = x_min, y_start
        self.xc = x_min + (np.arange(self.nrow) + 0.5) * self.dx
        self.yc = y_start + (np.arange(self.ny) + 0.5) * self.dy
        self._classify()
        self._coef_cache: dict[float, np.ndarray] = {}

    # ------------------------------------------------------------ geometry
    def _classify(self):
        X, Y = np.meshgrid(self.xc, self.yc, indexing="ij")
        kind = np.full(X.shape, _DIAPHRAGM, int)
        # per box: strip chord limits of the surface under it (for edge handling)
        self.x_le = np.full(X.shape, np.nan)
        self.x_te = np.full(X.shape, np.nan)
        for t in self.traps:
            inspan = (Y > t.y0) & (Y < t.y1)
            le, te = t.le(Y), t.te(Y)
            wing = inspan & (X >= le) & (X <= te)
            wake = inspan & (X > te) & (kind != _WING)
            kind[wake & (kind == _DIAPHRAGM)] = _WAKE
            kind[wing] = _WING
            self.x_le = np.where(wing, le, self.x_le)
            self.x_te = np.where(wing, te, self.x_te)
        self.kind = kind
        self.wing = kind == _WING
        # chordwise segment of each wing box: first/last wing box of a strip run
        # extends to the true leading/trailing edge
        fwd = np.ones_like(self.wing)
        fwd[1:] = ~self.wing[:-1]
        aft = np.ones_like(self.wing)
        aft[:-1] = ~self.wing[1:]
        self.first = self.wing & fwd
        self.last = self.wing & aft
        xf = np.where(self.first, self.x_le, X - 0.5 * self.dx)
        xa = np.where(self.last, self.x_te, X + 0.5 * self.dx)
        self.seg_fwd = np.where(self.wing, xf, np.nan)
        self.seg_aft = np.where(self.wing, xa, np.nan)
        self.seg_len = np.where(self.wing, xa - xf, 0.0)
        self.area = self.seg_len * self.dy
        self.X, self.Y = X, Y

    @property
    def n_wing(self) -> int:
        return int(self.wing.sum())

    def wing_points(self, where: str = "center") -> np.ndarray:
        """(Nw, 3) wing-box points: 'center' (normalwash) or 'segment' (force) midpoints."""
        x = self.X[self.wing] if where == "center" else 0.5 * (self.seg_fwd + self.seg_aft)[self.wing]
        return np.stack([x, self.Y[self.wing], np.zeros(x.size)], 1)

    # ------------------------------------------------------------ influence
    def coefficients(self, k: float, bref: float) -> np.ndarray:
        """Box influence coefficients ``C[m, n]`` (n >= 0), ``phi/U = C * w/U``.

        ``m`` = row offset downstream, ``n`` = |spanwise box offset|.
        """
        key = (round(float(k) / bref * self.dx, 12))
        if key in self._coef_cache:
            return self._coef_cache[key]
        b = self.beta
        M = self.mach
        h = self.dx
        lam = key * M * M / (b * b)           # (w/U) dx M^2/beta^2   (u in box units)
        kap = key * M / (b * b)
        nmax = min(self.nrow + 1, 2 * self.ny + 1)
        tq, wq = _gauss(self.n_gauss)
        C = np.zeros((self.nrow, nmax + 1), complex)
        m = np.arange(self.nrow)[:, None]
        n = np.arange(nmax + 1)[None, :]
        m, n = np.broadcast_arrays(m, n)
        m, n = m.ravel(), n.ravel()
        ok = n <= m + 1
        m, n = m[ok], n[ok]
        ulo_box = np.maximum(m - 0.5, 0.0)
        uhi_box = m + 0.5
        vlo, vhi = n - 0.5, n + 0.5
        # breakpoints in u where the clipping of the v-range changes
        bps = np.stack([ulo_box, np.clip(np.abs(vlo), ulo_box, uhi_box),
                        np.clip(np.abs(vhi), ulo_box, uhi_box), uhi_box], 1)
        bps.sort(axis=1)
        total = np.zeros(m.size, complex)
        for s in range(3):
            a, bb = bps[:, s], bps[:, s + 1]
            L = bb - a
            act = L > 1e-14
            if not act.any():
                continue
            # u = a + L t^2 (clusters points at the segment start: sqrt behaviour)
            u = a[act, None] + L[act, None] * tq[None, :] ** 2
            du = 2 * L[act, None] * tq[None, :] * wq[None, :]
            lo = np.maximum(vlo[act, None], -u)
            hi = np.minimum(vhi[act, None], u)
            valid = hi > lo
            with np.errstate(invalid="ignore", divide="ignore"):
                th_lo = np.arcsin(np.clip(np.where(valid, lo / u, 0.0), -1, 1))
                th_hi = np.arcsin(np.clip(np.where(valid, hi / u, 0.0), -1, 1))
            dth = np.where(valid, th_hi - th_lo, 0.0)
            if kap == 0.0:
                inner = dth
            else:
                th = th_lo[..., None] + dth[..., None] * tq[None, None, :]
                inner = dth * np.sum(wq * np.cos(kap * u[..., None] * np.cos(th)), axis=-1)
            total[act] += np.sum(du * np.exp(-1j * lam * u) * inner, axis=1)
        C[m, n] = -(h / (np.pi * b)) * total
        self._coef_cache[key] = C
        return C

    # ------------------------------------------------------------ solution
    def solve(self, k: float, bref: float, W: np.ndarray):
        """March the potential for wing normalwash ``W`` (Nw, nc) per unit U.

        Returns ``(phi, dcp)``: potential/U on the full grid (nrow, ny, nc) and
        the pressure jump on the wing boxes (Nw, nc).  The spanwise sums are
        convolutions (the coefficients depend on box offsets only), done with
        FFTs; the image half-wing is added for symmetric / antisymmetric grids.
        """
        from scipy.fft import fft, ifft, next_fast_len
        W = np.asarray(W, complex)
        if W.ndim == 1:
            W = W[:, None]
        nc = W.shape[1]
        ny = self.ny
        C = self.coefficients(k, bref)
        c00 = C[0, 0]
        nmax = C.shape[1] - 1
        nw = 2 * ny if self.sign != 0.0 else ny
        L = next_fast_len(nw + 2 * nmax + 1)
        Chat = fft(np.concatenate([C[:, :0:-1], C], axis=1), n=L, axis=1)     # (nrow, L)
        off = (ny if self.sign != 0.0 else 0) + nmax
        w = np.zeros((self.nrow, ny, nc), complex)
        phi = np.zeros_like(w)
        pte = np.zeros_like(w)                 # trailing-edge potential of the last wing boxes
        What = np.zeros((self.nrow, L, nc), complex)
        w[self.wing] = W
        kq = k / bref
        decay = np.exp(-1j * kq * self.dx)
        for i in range(self.nrow):
            if i > 0:
                acc = np.einsum("ml,mlc->lc", Chat[1:i + 1], What[i - 1::-1])
                known = ifft(acc, axis=0)[off:off + ny]
            else:
                known = np.zeros((ny, nc), complex)
            kd = self.kind[i]
            wing = kd == _WING
            phi[i, wing] = known[wing] + c00 * w[i, wing]
            dia = kd == _DIAPHRAGM
            w[i, dia] = -known[dia] / c00                       # phi = 0 off the wing
            wk = kd == _WAKE
            if wk.any():
                target = phi[i - 1] * decay
                new = wk & self.last[i - 1]                     # first wake box behind the TE
                if new.any():
                    target[new] = pte[i - 1, new] * np.exp(-1j * kq * (self.X[i, new] - self.x_te[i - 1, new]))[:, None]
                w[i, wk] = (target[wk] - known[wk]) / c00
                phi[i, wk] = target[wk]
            last = self.last[i]
            if last.any():
                xl = self.X[i, last]
                pw = self.wing[i - 1, last] if i > 0 else np.zeros(xl.size, bool)
                # linear extrapolation from the two last box centres; a strip of a
                # single box keeps its centre value
                slope = np.where(pw[:, None], (phi[i, last] - phi[i - 1, last]) / self.dx, 0.0)
                pte[i, last] = phi[i, last] + slope * (self.x_te[i, last] - xl)[:, None]
            wr = w[i] if self.sign == 0.0 else np.concatenate([self.sign * w[i, ::-1], w[i]], axis=0)
            What[i] = fft(wr, n=L, axis=0)
        # pressure: dCp = 4 (dphi/dx + i k/b phi), phi at the box edges
        prev = np.zeros_like(phi)
        prev[1:] = phi[:-1]
        nxt = np.zeros_like(phi)
        nxt[:-1] = phi[1:]
        pf = np.where(self.first[..., None], 0.0, 0.5 * (prev + phi))
        pa = np.where(self.last[..., None], pte, 0.5 * (phi + nxt))
        Lseg = np.where(self.wing, self.seg_len, 1.0)[..., None]
        dcp = 4.0 * ((pa - pf) / Lseg + 1j * kq * phi)
        return phi, dcp[self.wing]

    def table(self, mode_field: Callable, bref: float, kvals, names=None, gust: bool = True,
              x_gust0: float | None = None, verbose: bool = False):
        """GAF table for mode functions ``mode_field(points) -> (h, hx)``.

        Returns ``(GAFTable, dcp)`` with ``dcp`` (nk, Nw, nm[+1]) on the wing boxes.
        """
        kvals = np.atleast_1d(np.asarray(kvals, float))
        hc, hxc = mode_field(self.wing_points("center"))
        hf, _ = mode_field(self.wing_points("segment"))
        nm = hc.shape[1]
        x = self.X[self.wing]
        x0 = self.x0 if x_gust0 is None else x_gust0
        HA = hf * self.area[self.wing][:, None]
        Q = np.zeros((kvals.size, nm, nm), complex)
        Qg = np.zeros((kvals.size, nm), complex) if gust else None
        dcps = np.zeros((kvals.size, self.n_wing, nm + int(gust)), complex)
        for i, k in enumerate(kvals):
            W = hxc + 1j * (k / bref) * hc
            if gust:
                W = np.hstack([W, -np.exp(-1j * (k / bref) * (x - x0))[:, None]])
            _, dcp = self.solve(k, bref, W)
            dcps[i] = dcp
            Q[i] = HA.T @ dcp[:, :nm]
            if gust:
                Qg[i] = HA.T @ dcp[:, nm]
            if verbose:
                print(f"  Mach box k={k:.4f} done")
        names = list(names) if names else [f"mode{i + 1}" for i in range(nm)]
        return GAFTable(kvals, Q, Qg, self.mach, float(bref), names), dcps

    def sample(self, values: np.ndarray, points: np.ndarray) -> np.ndarray:
        """Nearest wing-box value at ``points`` (for plotting on the DLM boxes)."""
        idx = -np.ones(self.wing.shape, int)
        idx[self.wing] = np.arange(self.n_wing)
        i = np.clip(np.floor((points[:, 0] - self.x0) / self.dx).astype(int), 0, self.nrow - 1)
        yy = np.abs(points[:, 1]) if self.sign != 0.0 else points[:, 1]
        j = np.clip(np.floor((yy - self.y0) / self.dy).astype(int), 0, self.ny - 1)
        sel = idx[i, j]
        if np.any(sel < 0):   # point on a jagged edge: nearest wing box
            wi = np.argwhere(self.wing)
            for p in np.where(sel < 0)[0]:
                d = (wi[:, 0] - i[p]) ** 2 + (wi[:, 1] - j[p]) ** 2
                sel[p] = idx[tuple(wi[np.argmin(d)])]
        return values[..., sel, :] if values.ndim == 3 else values[sel]

    def summary(self) -> str:
        return (f"Mach box grid: M = {self.mach}, beta = {self.beta:.4f}, dx = {self.dx:.4g}, "
                f"dy = {self.dy:.4g}, {self.nrow} x {self.ny} boxes ({self.n_wing} on the wing, "
                f"{int((self.kind == _DIAPHRAGM).sum())} diaphragm, {int((self.kind == _WAKE).sum())} wake), "
                f"symmetry {self.symmetry}")


def machbox_gaf_table(mesh: AeroMesh, modes, mach: float, bref: float, kvals, nx: int = 60,
                      gust: bool = True, keep_cp: bool = True, x_gust0: float | None = None,
                      n_gauss: int = 10, verbose: bool = False, names=None) -> GAFTable:
    """Mach-box GAF table.  ``modes``: AeroModes on ``mesh`` or ``f(points) -> (h, hx)``.

    With ``keep_cp`` the Mach-box pressures are sampled (nearest box) at the
    DLM box centres so that ``table.cp`` matches ``mesh`` like a DLM table.
    """
    mb = MachBox(mesh, mach, nx=nx, n_gauss=n_gauss)
    if verbose:
        print("  " + mb.summary())
    if isinstance(modes, AeroModes):
        field = aeromodes_field(mesh, modes)
        names = names or modes.names
    else:
        field = modes
    x0 = float(mesh.corners[:, :, 0].min()) if (x_gust0 is None and isinstance(mesh, AeroMesh)) else x_gust0
    tab, dcps = mb.table(field, bref, kvals, names, gust, x0, verbose)
    if keep_cp and isinstance(mesh, AeroMesh):
        s = mb.sample(dcps, mesh.center)            # (nk, N, nm+g)
        nm = tab.Q.shape[1]
        tab.cp = s[:, :, :nm]
        tab.cpg = s[:, :, nm] if gust else None
    tab.machbox = mb
    return tab


def supersonic_gaf_table(mesh: AeroMesh, modes, mach: float, bref: float, kvals,
                         method: str = "machbox", gust: bool = True, keep_cp: bool = True,
                         x_gust0: float | None = None, verbose: bool = False, **options) -> GAFTable:
    """Supersonic GAF table, ``method`` 'piston' or 'machbox' (see module docstring)."""
    method = str(method).lower()
    if method == "piston":
        return piston_theory_table(mesh, modes, mach, bref, kvals, gust=gust, keep_cp=keep_cp,
                                   x_gust0=x_gust0, **options)
    if method in ("machbox", "mach_box"):
        return machbox_gaf_table(mesh, modes, mach, bref, kvals, gust=gust, keep_cp=keep_cp,
                                 x_gust0=x_gust0, verbose=verbose, **options)
    raise ValueError(f"unknown supersonic method {method!r} (piston | machbox)")


# =====================================================================
# exact 2-D reference (independent of the Mach box discretisation)
# =====================================================================
def garrick_rubinow_2d(mach: float, k: float, bref: float, chord: float, x_ref: float,
                       n: int = 200) -> np.ndarray:
    """Exact linear 2-D supersonic GAFs of a flat plate in plunge and pitch.

    Garrick & Rubinow (1946): upper-surface potential of a supersonic
    oscillating airfoil (leading edge at ``x = 0``)::

        phi(x) = -(1/beta) int_0^x w(xi) exp(-i lam X) J0(kap X) dxi,  X = x - xi

    with ``lam = (k/b) M^2/beta^2`` and ``kap = (k/b) M/beta^2``.  Returns the
    2x2 matrix ``Q[i, j] = int_0^c h_i dCp_j dx`` per unit span for the modes
    plunge ``h = 1`` and pitch ``h = -(x - x_ref)``, with
    ``dCp = 4 (dphi/dx + i (k/b) phi)`` integrated by parts.  Gauss-Legendre
    quadrature with ``n`` points per integral.
    """
    from scipy.special import j0
    beta = _beta(mach)
    kq = k / bref
    lam, kap = kq * mach * mach / beta ** 2, kq * mach / beta ** 2
    t, w = np.polynomial.legendre.leggauss(n)
    t, w = 0.5 * (t + 1.0), 0.5 * w

    def normalwash(xx):
        return np.stack([1j * kq * np.ones_like(xx), -1.0 - 1j * kq * (xx - x_ref)], -1)

    def phi(xv):                                       # (P,) -> (P, 2)
        xv = np.atleast_1d(xv)[:, None]
        xi = xv * t[None, :]
        X = xv - xi
        ker = (xv * w[None, :]) * np.exp(-1j * lam * X) * j0(kap * X)
        return -(1.0 / beta) * np.einsum("pq,pqj->pj", ker, normalwash(xi))

    x, wx = chord * t, chord * w
    ph = phi(x)
    ph_te = phi(np.array([chord]))[0]
    h = np.stack([np.ones_like(x), -(x - x_ref)], -1)
    h_te = np.array([1.0, -(chord - x_ref)])
    hx = np.array([0.0, -1.0])
    return 4.0 * (h_te[:, None] * ph_te[None, :]
                  - hx[:, None] * (wx @ ph)[None, :]
                  + 1j * kq * np.einsum("xi,xj,x->ij", h, ph, wx))
