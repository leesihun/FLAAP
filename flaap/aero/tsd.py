"""Transonic Small-Disturbance (TSD) aerodynamics for thin wings.

A 3-D, frequency-domain, time-linearised TSD solver in the spirit of
CAP-TSD / XTRAN3S / ATRAN3S, intended for flutter of thin swept wings
(e.g. AGARD 445.6, NACA 65A004, M up to ~0.96).

References
----------
* Batina, J.T., "Efficient algorithm for solution of the unsteady
  transonic small-disturbance equation", J. Aircraft 25(7), 1988, and the
  CAP-TSD code description, NASA TM-101565.
* Guruswamy, G.P. & Goorjian, P.M. (ATRAN3S) and Borland, C.J. &
  Rizzetta, D.P. (XTRAN3S), time-accurate 3-D TSD codes for wing flutter.
* Murman, E.M. & Cole, J.D., "Calculation of plane steady transonic
  flows", AIAA J. 9(1), 1971; Murman, E.M., "Analysis of embedded shock
  waves calculated by relaxation methods", AIAA J. 12(5), 1974 (fully
  conservative switching / shock-point operator).
* Berenger, J.-P., J. Comput. Phys. 114, 1994; Chew, W.C. & Weedon, W.H.,
  Microw. Opt. Tech. Lett. 7, 1994 (PML / complex coordinate stretching).

Non-dimensionalisation (all internal quantities)
------------------------------------------------
* lengths by the reference length ``Lref`` = root chord of the surface:
  ``x^ = x / Lref`` etc.;  time by ``Lref / U``
* perturbation potential by ``U Lref`` (so ``phi_x^ = u / U``)
* reduced frequency on the reference length: ``kbar = w Lref / U``
  (FLAAP's ``k = w bref / U`` is converted with ``kbar = k Lref / bref``)

Equations
---------
Steady (classical TSD, conservative form)::

    d/dx[(1 - M^2) phi_x - (g+1)/2 M^2 phi_x^2] + phi_yy + phi_zz = 0
    Cp = -2 phi_x

Unsteady, harmonic ``phi = phi0 + phi1 exp(i kbar t)`` linearised about
the steady mean flow ``phi0``::

    d/dx[A phi1_x] + phi1_yy + phi1_zz + M^2 kbar^2 phi1 - 2 i kbar M^2 phi1_x = 0
    A = 1 - M^2 - (g+1) M^2 phi0_x
    Cp1 = -2 (phi1_x + i kbar phi1)

With ``phi0 = 0`` this is the exact linear compressible potential
equation, so at subsonic Mach numbers and zero thickness the TSD results
must reproduce the doublet-lattice method (see ``tests/test_tsd.py``).

Symmetry in z
-------------
Only the upper half domain ``z >= 0`` is discretised:

* thickness (mean flow) problem - symmetric in z: Neumann
  ``phi_z = dz_t/dx`` on wing nodes, ``phi_z = 0`` elsewhere on ``z=0``.
  (So the nonlinear steady solution is restricted to symmetric sections
  at zero incidence; this is the mean flow used for flutter.)
* lifting (unsteady) problem - antisymmetric in z: Neumann
  ``phi1_z = w/U = dh/dx + i (w/U) h`` on the wing, ``phi1 = 0`` on
  ``z=0`` upstream and outboard of the wing, and in the wake
  ``phi1(x,y,0) = phi1_TE(y) exp(-i kbar (x - x_te))`` (half the convected
  circulation jump; this carries the Kutta condition, the pressure jump
  vanishes in the wake).  Pressure jump (upper minus lower Cp sign as in
  :mod:`flaap.aero.dlm`): ``dCp = Cp_l - Cp_u = 4 (phi1_x + i kbar phi1)``
  on the upper surface, positive for lift.

The spanwise symmetry plane is ``y = 0`` (symmetric half model).

Grid and discretisation
-----------------------
A sheared grid follows the leading/trailing edges::

    x^(u, y) = a(y) + b(y) u ,   physical LE at u = u_le, TE at u = u_te

(``a``, ``b`` linear in ``y`` on the wing and smoothly frozen outboard of
the tip).  Lines of constant ``y`` are streamwise, so upwinding in ``u`` is
upwinding in ``x``.  The equation multiplied by the Jacobian ``b`` is
written in conservation form::

    d_u(F - S G) + d_y(b G) + d_z(b phi_z) = 0,
    F = (1-M^2) p - (g+1)/2 M^2 p^2 ,  p = phi_u / b  (= phi_x)
    G = phi_y = phi_eta - (S/b) phi_u ,  S = a' + b' u

and discretised by a vertex-centred finite-volume scheme (nodes on
``z = 0``; half cells on ``z = 0`` and ``y = 0``).  The x-flux uses the
fully conservative Murman (1974) switching ``(1-mu_i) [F_{i+1/2} -
F_{i-1/2}] + mu_{i-1} [F_{i-1/2} - F_{i-3/2}]`` with ``mu = 1`` where
``A < 0`` (central at subsonic points, backward at supersonic points,
shock-point and sonic-point operators at transitions); the unsteady
``-2 i kbar M^2 phi_x`` flux uses the same switches (frozen from the mean
flow).  The cross-derivative (shear) terms are central.  The steady
problem is solved by Newton's method (switches re-evaluated after each
step, back-tracking line search).

Chordwise wing nodes are uniform, ``u = (m + 1/2)/nchord``, spanwise nodes
``y = j dy`` with the tip half-way between nodes, ``(nspan - 1/2) dy``.
With the cell-face (u = 0, 1) edge placement the discrete solution behaves
as if the leading edge were ~0.85 cells upstream and the trailing edge
~0.15 cells upstream of the faces (first-order edge error caused by the
square-root LE singularity; calibrated against thin-airfoil theory and
the DLM).  The mapping therefore places the physical LE/TE at
``u_le = -0.85/nchord`` and ``u_te = 1 - 0.15/nchord`` (``edge_offset``
argument; ``None`` for the plain face placement).  This removes the
O(1/nchord) centre-of-pressure shift (~0.7 chordwise cells forward).

Linear systems: every operator couples only neighbouring u-planes (plus
the i-2 plane of the supersonic backward operator), so it is solved
exactly by block-banded Gaussian elimination over the u-planes with dense
LAPACK kernels on the (ny*nz)^2 blocks (``linear_solver='block'``, much
faster than SuperLU for these 3-D grids; ``'splu'`` is kept as a check).

Generalised forces are integrated by parts on every strip (physical
``x``, ``h``), which avoids the leading-edge pressure singularity::

    int h dCp dx = 4 Lref [ h_TE phi_TE - int (dh/dx - i (kbar/Lref) h) phi dx ]

(``integration='direct'`` sums nodal dCp and is much less accurate.)

Far field: ``phi = 0`` upstream, ``phi_x + i kbar phi = 0`` downstream,
Robin (``d phi/dn = -phi (n.r)/r^2``, a cross-flow dipole decay) on the
outboard and upper boundaries.  For ``kbar > 0`` perfectly matched layers
(complex coordinate stretching, ``x~ = x - i n sigma(x)``, Berenger 1994;
Chew & Weedon 1994) occupy the outer ``pml_frac`` of the upstream,
outboard and upper far-field regions (the boundary rows use the same
stretched coordinates).  Without them, acoustic waves reflected from the
truncated boundaries spoil high-subsonic / high-k results (e.g. 30-40 %
errors at M = 0.9, k = 0.6 on the AGARD planform that also depend on the
domain size; with the PML: ~2 % and domain independent).  No PML is used
downstream (convected waves can have negative x phase velocity there).

Limitations
-----------
* thin planar trapezoidal wing, root on the symmetry plane ``y = 0``
* this class solves the nonlinear steady flow only for the symmetric
  (thickness) problem; lift enters linearly through the time-linearised
  problem (k = 0 gives the steady lifting increment about the thickness
  mean flow).  :class:`flaap.aero.tsd_lifting.LiftingTSDSolver` solves the
  nonlinear mean flow with incidence and camber (shock positions vs alpha)
  and the flutter problem about it
* frozen shock position in the time-linearised problem (shock motion is
  represented only through the shock-point operator); cross-derivative
  shear terms are not upwinded, so strongly swept supersonic pockets are
  handled approximately
* classical TSD coefficients (no NASA Ames modified terms, no viscous /
  boundary-layer coupling); TSD Cp* differs from the isentropic value
* the edge-offset calibration assumes the default ``dz1 = 0.6/nchord``
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence
import time
from types import SimpleNamespace

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import brentq

from ..geometry import Surface

GAMMA = 1.4

__all__ = ["TSDSolver", "SteadyResult", "tsd_gaf_table", "modes_from_structure",
           "rigid_mode_funcs", "parabolic_arc", "naca4_thickness", "naca65a004",
           "cp_critical"]


# ---------------------------------------------------------------- sections
def parabolic_arc(tc: float) -> Callable:
    """Biconvex (parabolic-arc) section, total thickness ``t(xi) = 4 tc xi(1-xi)``."""
    def t(xi):
        xi = np.asarray(xi, float)
        return np.where((xi >= 0) & (xi <= 1), 4 * tc * xi * (1 - xi), 0.0)

    def slope(xi):   # d(t/2)/dxi  (upper-surface slope)
        xi = np.asarray(xi, float)
        return np.where((xi >= 0) & (xi <= 1), 2 * tc * (1 - 2 * xi), 0.0)
    t.slope = slope
    return t


def naca4_thickness(tc: float, closed_te: bool = True) -> Callable:
    """NACA 4-digit symmetric thickness ``t(xi)`` (00xx)."""
    a4 = -0.1036 if closed_te else -0.1015
    co = np.array([0.2969, -0.1260, -0.3516, 0.2843, a4])

    def t(xi):
        xi = np.clip(np.asarray(xi, float), 0, 1)
        return 10 * tc * (co[0] * np.sqrt(xi) + co[1] * xi + co[2] * xi ** 2
                          + co[3] * xi ** 3 + co[4] * xi ** 4)

    def slope(xi):
        xi = np.clip(np.asarray(xi, float), 1e-12, 1)
        return 5 * tc * (0.5 * co[0] / np.sqrt(xi) + co[1] + 2 * co[2] * xi
                         + 3 * co[3] * xi ** 2 + 4 * co[4] * xi ** 3)
    t.slope = slope
    return t


# NACA 65-010 half-thickness ordinates (Abbott & von Doenhoff), percent chord.
# NACA 65A004 half-thickness ordinates, % chord (Abbott & von Doenhoff; the
# "A" modification gives straight sides aft of ~80 % chord, finite TE)
_N65A_X = np.array([0, 0.5, 0.75, 1.25, 2.5, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                    55, 60, 65, 70, 75, 80, 85, 90, 95, 100]) / 100
_N65A_Y = np.array([0.0, 0.311, 0.378, 0.481, 0.656, 0.877, 1.062, 1.216, 1.463, 1.649,
                    1.790, 1.894, 1.962, 1.997, 1.996, 1.952, 1.867, 1.742, 1.584, 1.400,
                    1.193, 0.966, 0.728, 0.490, 0.249, 0.009]) / 100


def naca65a004() -> Callable:
    """NACA 65A004 total thickness ``t(xi)`` from the tabulated ordinates
    (pchip in ``sqrt(xi)`` to follow the round leading edge), with
    ``.slope(xi)`` = upper-surface slope d(t/2)/dxi."""
    from scipy.interpolate import PchipInterpolator
    f = PchipInterpolator(np.sqrt(_N65A_X), 2 * _N65A_Y)
    fd = f.derivative()

    def t(xi):
        xi = np.clip(np.asarray(xi, float), 0, 1)
        return f(np.sqrt(xi))

    def slope(xi):
        xi = np.clip(np.asarray(xi, float), 1e-8, 1)
        return 0.5 * fd(np.sqrt(xi)) * 0.5 / np.sqrt(xi)
    t.slope = slope
    return t


def cp_critical(mach: float, gamma: float = GAMMA, tsd: bool = False) -> float:
    """Critical (sonic) pressure coefficient; ``tsd=True`` gives the TSD value
    ``-2 (1 - M^2) / ((g + 1) M^2)`` consistent with ``A = 0``."""
    m2 = mach * mach
    if tsd:
        return -2 * (1 - m2) / ((gamma + 1) * m2)
    return 2 / (gamma * m2) * (((2 + (gamma - 1) * m2) / (gamma + 1)) ** (gamma / (gamma - 1)) - 1)


# ---------------------------------------------------------------- mode funcs
def rigid_mode_funcs(x_ref: float, kinds: Sequence[str] = ("plunge", "pitch")):
    """Rigid-body mode functions ``f(x, y) -> (h, dh/dx)``.

    plunge ``h = 1``; pitch (nose up) ``h = -(x - x_ref)``; roll ``h = y``.
    Returns ``(funcs, names)``.
    """
    funcs = []
    for k in kinds:
        if k == "plunge":
            funcs.append(lambda x, y: (np.ones_like(np.asarray(x, float)),
                                       np.zeros_like(np.asarray(x, float))))
        elif k == "pitch":
            funcs.append(lambda x, y, xr=x_ref: (-(np.asarray(x, float) - xr),
                                                 -np.ones_like(np.asarray(x, float))))
        elif k == "roll":
            funcs.append(lambda x, y: (np.asarray(y, float) * np.ones_like(np.asarray(x, float)),
                                       np.zeros_like(np.asarray(x, float))))
        else:
            raise ValueError(k)
    return funcs, list(kinds)


class _SplineMode:
    """Mode function built from nodal ``uz`` with the infinite-plate spline."""

    def __init__(self, owner, col):
        self.owner = owner
        self.col = col

    def __call__(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
        G, Gx = self.owner._mats(x, y)
        u = self.owner.uz[:, self.col]
        return (G @ u).reshape(x.shape), (Gx @ u).reshape(x.shape)


class _SplineModeSet:
    def __init__(self, nodes, uz):
        self.nodes = np.asarray(nodes, float)[:, :2]
        self.uz = np.asarray(uz, float)
        self._cache = {}

    def _mats(self, x, y):
        from ..structures.spline import ips_matrices
        pts = np.column_stack([np.ravel(x), np.ravel(y)])
        key = (pts.shape, hash(pts.tobytes()))
        if key not in self._cache:
            if len(self._cache) > 8:
                self._cache.clear()
            self._cache[key] = ips_matrices(self.nodes, pts)
        return self._cache[key]


def modes_from_structure(nodes, shapes=None, names=None, mode_idx=None):
    """Mode functions from a structural modal model (IPS spline of ``uz``).

    ``nodes`` may be a :class:`flaap.structures.modal.ModalModel` (then
    ``shapes`` is taken from it) or an (Nn, 3) node array with ``shapes``
    (Nn, 6, nm) or (Nn, nm) (normal displacement).  Returns ``(funcs,
    names)`` where each ``f(x, y) -> (h, dh/dx)`` (value and analytic
    x-slope of the infinite-plate spline, :func:`ips_matrices`).
    """
    if shapes is None and hasattr(nodes, "shapes"):
        mm = nodes
        nodes, shapes = mm.nodes, mm.shapes
        if names is None and getattr(mm, "labels", None):
            names = list(mm.labels)
    shapes = np.asarray(shapes, float)
    uz = shapes[:, 2, :] if shapes.ndim == 3 else shapes
    if mode_idx is not None:
        uz = uz[:, list(mode_idx)]
        if names is not None:
            names = [names[i] for i in mode_idx]
    ms = _SplineModeSet(nodes, uz)
    funcs = [_SplineMode(ms, c) for c in range(uz.shape[1])]
    names = list(names) if names is not None else [f"mode{i + 1}" for i in range(uz.shape[1])]
    return funcs, names


# ---------------------------------------------------------------- grid utils
def _stretch(d0: float, length: float, n: int) -> np.ndarray:
    """``n`` geometric increments starting with ``d0`` summing to ``length``."""
    if n <= 0:
        return np.zeros(0)
    if abs(n * d0 - length) < 1e-12 * length:
        return np.full(n, d0)

    def f(r):
        return d0 * (r ** n - 1) / (r - 1) - length if abs(r - 1) > 1e-12 else n * d0 - length
    lo, hi = (1.0 + 1e-10, 10.0) if n * d0 < length else (1e-3, 1.0 - 1e-10)
    r = brentq(f, lo, hi)
    return d0 * r ** np.arange(n)


def _d1_faces(x):
    """(n-1, n) face difference operator d/dx."""
    n = x.size
    h = np.diff(x)
    r = np.arange(n - 1)
    return sp.csr_matrix((np.r_[-1 / h, 1 / h], (np.r_[r, r], np.r_[r, r + 1])), shape=(n - 1, n))


def _avg_faces(n):
    r = np.arange(n - 1)
    return sp.csr_matrix((np.full(2 * (n - 1), 0.5), (np.r_[r, r], np.r_[r, r + 1])), shape=(n - 1, n))


def _central(x, zero_first=False):
    """(n, n) node derivative: central inside, one-sided at the ends."""
    n = x.size
    rows, cols, vals = [], [], []
    for i in range(n):
        if i == 0:
            if zero_first:
                continue
            d = x[1] - x[0]
            rows += [0, 0]; cols += [0, 1]; vals += [-1 / d, 1 / d]   # noqa: E702
        elif i == n - 1:
            d = x[-1] - x[-2]
            rows += [i, i]; cols += [i - 1, i]; vals += [-1 / d, 1 / d]   # noqa: E702
        else:
            d = x[i + 1] - x[i - 1]
            rows += [i, i]; cols += [i - 1, i + 1]; vals += [-1 / d, 1 / d]   # noqa: E702
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n))


def _nd_order(shape, leaf: int = 40, sep_u: int = 2) -> np.ndarray:
    """Geometric nested-dissection ordering of a structured (u, y, z) grid.

    The box is bisected recursively along its longest index direction by
    separator planes (two planes thick in ``u`` because the supersonic
    backward operator reaches ``i-2``), separators ordered last.  This gives
    far less LU fill than COLAMD/MMD for the 3-D stencil.
    """
    idx = np.arange(int(np.prod(shape))).reshape(shape)
    out = []
    stack = [(0, shape[0], 0, shape[1], 0, shape[2], False)]
    # iterative post-order: (box, emitted?) ; separators are pushed so they
    # come out after both halves
    while stack:
        i0, i1, j0, j1, k0, k1, is_sep = stack.pop()
        if is_sep:
            out.append(idx[i0:i1, j0:j1, k0:k1].ravel())
            continue
        ext = [i1 - i0, j1 - j0, k1 - k0]
        n = ext[0] * ext[1] * ext[2]
        if n <= 0:
            continue
        d = int(np.argmax(ext))
        w = sep_u if d == 0 else 1
        if n <= leaf or ext[d] <= w + 2:
            out.append(idx[i0:i1, j0:j1, k0:k1].ravel())
            continue
        a = [i0, j0, k0][d] + (ext[d] - w) // 2
        lo = [i0, i1, j0, j1, k0, k1]
        hi = [i0, i1, j0, j1, k0, k1]
        sep = [i0, i1, j0, j1, k0, k1]
        lo[2 * d + 1] = a
        hi[2 * d] = a + w
        sep[2 * d], sep[2 * d + 1] = a, a + w
        # stack is LIFO: push separator first, then high, then low
        stack.append((*sep, True))
        stack.append((*hi, False))
        stack.append((*lo, False))
    return np.concatenate(out)


def _div(n, shift=0):
    """(n, n-1) node <- faces: row i = f_{i+1/2-shift} - f_{i-1/2-shift}."""
    rows, cols, vals = [], [], []
    for i in range(n):
        fp, fm = i - shift, i - 1 - shift
        if 0 <= fp <= n - 2:
            rows.append(i); cols.append(fp); vals.append(1.0)   # noqa: E702
        if 0 <= fm <= n - 2:
            rows.append(i); cols.append(fm); vals.append(-1.0)  # noqa: E702
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n - 1))


def _dual(x):
    """Dual-cell widths (half widths at the ends)."""
    v = np.empty_like(x)
    v[1:-1] = 0.5 * (x[2:] - x[:-2])
    v[0] = 0.5 * (x[1] - x[0])
    v[-1] = 0.5 * (x[-1] - x[-2])
    return v


def _blas_threads(P: int) -> int:
    """BLAS thread count for dense P x P block operations.  Small blocks gain
    nothing from many threads, and on a loaded machine an oversubscribed
    OpenBLAS pool spin-waits so badly that a 700 x 700 LU can take seconds
    instead of milliseconds.  ``FLAAP_BLAS_THREADS`` overrides."""
    import os
    env = os.environ.get("FLAAP_BLAS_THREADS")
    if env:
        return max(1, int(env))
    return int(min(os.cpu_count() or 1, 8, max(1, P // 400)))


def _block_solve(K, R, P: int, nb: int, max_bytes: float = 1.5e9):
    """:func:`_block_solve_impl` with the BLAS pool limited to
    :func:`_blas_threads` (when ``threadpoolctl`` is available)."""
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:                            # pragma: no cover
        return _block_solve_impl(K, R, P, nb, max_bytes)
    with threadpool_limits(limits=_blas_threads(P), user_api="blas"):
        return _block_solve_impl(K, R, P, nb, max_bytes)


def _block_solve_impl(K, R, P: int, nb: int, max_bytes: float = 1.5e9):
    """Exact solution of ``K x = R`` for an operator that is block banded over
    ``nb`` u-planes of size ``P``: lower block bandwidth <= 2 (the i-2 entries
    come from the supersonic backward operator), upper bandwidth 1.

    Block Gaussian elimination (generalised Thomas algorithm) with LU
    (partial pivoting) of the dense diagonal blocks; the right-hand sides are
    processed during the forward sweep so only ``D_c^-1 C_c`` is stored.  If
    storing all ``nb`` dense blocks would exceed ``max_bytes`` the forward
    sweep is check-pointed every ~sqrt(nb) planes and recomputed segment by
    segment during the back substitution (about twice the work, memory
    ~2 sqrt(nb) blocks).
    """
    import scipy.linalg as sla
    K = K.tocsr()
    R = np.asarray(R)
    vec = R.ndim == 1
    R = R.reshape(K.shape[0], -1)
    dtype = np.result_type(K.dtype, R.dtype)
    slabs = [K[c * P:(c + 1) * P] for c in range(nb)]
    for c in range(nb - 2):
        if slabs[c][:, (c + 2) * P:].nnz:
            raise RuntimeError("operator is not block banded (upper bandwidth > 1)")

    def step(c, prev1, prev2):
        """Eliminate plane c given (Cp, dp) of planes c-1 and c-2."""
        slab = slabs[c]
        D = slab[:, c * P:(c + 1) * P].toarray().astype(dtype)
        d = np.array(R[c * P:(c + 1) * P], dtype=dtype)
        if c >= 1:
            A = slab[:, (c - 1) * P:c * P]
            if c >= 2:
                E = slab[:, (c - 2) * P:(c - 1) * P]
                if E.nnz:
                    d -= E @ prev2[1]
                    A = A.toarray() - E @ prev2[0]
            D -= A @ prev1[0]
            d -= A @ prev1[1]
        lu = sla.lu_factor(D, overwrite_a=True, check_finite=False)
        Cp = None
        if c + 1 < nb:
            C = slab[:, (c + 1) * P:(c + 2) * P].toarray().astype(dtype)
            Cp = sla.lu_solve(lu, C, overwrite_b=True, check_finite=False)
        return Cp, sla.lu_solve(lu, d, overwrite_b=True, check_finite=False)

    x = np.empty((nb * P, R.shape[1]), dtype)
    block_bytes = P * P * np.dtype(dtype).itemsize
    none = (None, None)
    if nb * block_bytes <= max_bytes:
        st = [None] * nb
        for c in range(nb):
            st[c] = step(c, st[c - 1] if c >= 1 else none, st[c - 2] if c >= 2 else none)
        x[(nb - 1) * P:] = st[nb - 1][1]
        for c in range(nb - 2, -1, -1):
            x[c * P:(c + 1) * P] = st[c][1] - st[c][0] @ x[(c + 1) * P:(c + 2) * P]
    else:
        seg = max(2, int(np.ceil(np.sqrt(nb))))
        starts = list(range(0, nb, seg))
        ckpt = {}
        p1 = p2 = none
        for c in range(nb):                       # forward sweep, keep check-points
            if c in starts:
                ckpt[c] = (p1, p2)
            cur = step(c, p1, p2)
            p1, p2 = cur, p1
        x[(nb - 1) * P:] = p1[1]
        for c0 in reversed(starts):               # recompute and back-substitute
            c1 = min(c0 + seg, nb)
            p1, p2 = ckpt[c0]
            st = []
            for c in range(c0, c1):
                cur = step(c, p1, p2)
                st.append(cur)
                p1, p2 = cur, p1
            for c in range(c1 - 1, c0 - 1, -1):
                if c == nb - 1:
                    continue
                Cp, dp = st[c - c0]
                x[c * P:(c + 1) * P] = dp - Cp @ x[(c + 1) * P:(c + 2) * P]
    return x[:, 0] if vec else x


class _PermLU:
    """LU of ``K[p][:, p]`` presented as a solver for ``K``."""

    def __init__(self, lu, p):
        self.lu, self.p = lu, p
        self.nnz = lu.L.nnz + lu.U.nnz

    def solve(self, b):
        b = np.asarray(b)
        x = np.empty(b.shape, np.result_type(b.dtype, self.lu.L.dtype))
        x[self.p] = self.lu.solve(np.ascontiguousarray(b[self.p]))
        return x


@dataclass
class SteadyResult:
    phi: np.ndarray            # (NU, NY, NZ) mean-flow potential
    cp: np.ndarray             # (nspan, nchord) upper-surface Cp = -2 phi_x on the wing
    x: np.ndarray              # (nspan, nchord) physical x of wing nodes
    y: np.ndarray              # (nspan,)
    mach: float
    n_supersonic: int
    converged: bool
    iterations: int
    residuals: list = field(default_factory=list)
    cp_star: float = np.nan    # TSD critical Cp


class TSDSolver:
    """3-D frequency-domain TSD solver for a planar trapezoidal half wing.

    Parameters
    ----------
    surface : flaap.geometry.Surface (root on ``y = 0``, planar)
    mach : free-stream Mach number (< 1)
    thickness : callable ``t(xi)`` total thickness / chord vs chord fraction
        (optionally with a ``.slope(xi)`` attribute = d(t/2)/dxi), or None
    nchord, nspan : grid nodes on the wing chord / half span
    n_up, n_down : chordwise nodes upstream / downstream of the wing
    n_out : spanwise nodes outboard of the tip
    nz : nodes in z (including z = 0)
    x_up, x_down : upstream / downstream extent in local chords
    y_out : extent outboard of the tip, z_top : height (units of root chord;
        default ``max(3, 1.5 s)`` and ``max(3, 2 s)`` with ``s`` = semispan/c_r)
    dz1 : first z spacing (root chords); default ``0.6 / nchord``
    symmetry : only 'symmetric' (half model, symmetric about y = 0)
    linear_solver : 'block' (default, exact block-banded elimination over
        u-planes) or 'splu' (SuperLU with nested dissection; slower)
    max_bytes : memory budget of the block solver's stored dense blocks;
        above it the elimination is check-pointed (~2x time)
    edge_offset : (d_le, d_te) in chordwise cells, computational position of
        the physical LE/TE (see module notes); None = plain face placement
    pml, pml_frac, pml_sigma : perfectly matched layers for kbar > 0 in the
        outer ``pml_frac`` of the upstream / outboard / top regions, with
        total imaginary stretch ``pml_sigma`` x layer thickness
    """

    def __init__(self, surface: Surface, mach: float, thickness: Callable | None = None, *,
                 nchord: int = 24, nspan: int = 10, n_up: int = 12, n_down: int = 14,
                 n_out: int = 10, nz: int = 22, x_up: float = 4.0, x_down: float = 6.0,
                 y_out: float | None = None, z_top: float | None = None,
                 dz1: float | None = None, symmetry: str = "symmetric",
                 gamma: float = GAMMA, linear_solver: str = "block",
                 pivot_thresh: float = 0.1, edge_offset=(0.85, 0.15),
                 pml: bool = True, pml_frac: float = 0.4, pml_sigma: float = 1.0,
                 max_bytes: float = 1.5e9,
                 verbose: bool = False):
        if not 0.0 <= mach < 1.0:
            raise ValueError("TSDSolver requires 0 <= M < 1")
        if symmetry != "symmetric":
            raise NotImplementedError("only the symmetric half model (y = 0 plane) is supported")
        r = np.asarray(surface.root_le, float)
        t = np.asarray(surface.tip_le, float)
        if abs(r[1]) > 1e-9 * max(1.0, surface.span):
            raise ValueError("TSD solver expects the wing root on the symmetry plane y = 0")
        if abs(t[2] - r[2]) > 1e-6 * surface.span:
            raise ValueError("TSD solver expects a planar wing (no dihedral)")
        self.surface = surface
        self.mach = float(mach)
        self.gamma = float(gamma)
        self.thickness = thickness
        self.symmetry = symmetry
        if linear_solver not in ("block", "splu"):
            raise ValueError("linear_solver must be 'block' or 'splu'")
        self.linear_solver = linear_solver
        self.pivot_thresh = pivot_thresh
        self._perm = None
        self.max_bytes = float(max_bytes)
        self.pml, self.pml_frac, self.pml_sigma = bool(pml), float(pml_frac), float(pml_sigma)
        self._pml_ops = None
        self.verbose = verbose
        self.Lref = float(surface.root_chord)
        L = self.Lref
        self.s = (t[1] - r[1]) / L                       # nondim semispan
        self.a0 = r[0] / L
        self.tan_le = (t[0] - r[0]) / (t[1] - r[1])
        self.lam = surface.tip_chord / surface.root_chord
        self.bp0 = (self.lam - 1.0) / self.s
        self.blend = self.lam if self.bp0 >= 0 else min(self.lam, 0.5 * self.lam / abs(self.bp0))
        self.nchord, self.nspan = int(nchord), int(nspan)
        y_out = max(3.0, 1.5 * self.s) if y_out is None else y_out
        z_top = max(3.0, 2.0 * self.s) if z_top is None else z_top
        dz1 = 0.6 / nchord if dz1 is None else dz1
        # computational positions of the physical LE / TE (see module notes)
        dle, dte = (0.0, 0.0) if edge_offset is None else map(float, edge_offset)
        self.edge_offset = (dle, dte)
        self.u_le = -dle / self.nchord
        self.u_te = 1.0 - dte / self.nchord
        self.gmap = 1.0 / (self.u_te - self.u_le)
        self._build_grid(n_up, n_down, n_out, nz, x_up, x_down, y_out, z_top, dz1)
        self._build_operators()
        self._mean = None           # SteadyResult

    # ------------------------------------------------------------------ grid
    def _map(self, y):
        y = np.asarray(y, float)
        s, Lb = self.s, self.blend
        out = y > s
        e = np.exp(-np.where(out, y - s, 0.0) / Lb)
        a = np.where(out, self.a0 + self.tan_le * (s + Lb * (1 - e)), self.a0 + self.tan_le * y)
        ap = np.where(out, self.tan_le * e, self.tan_le)
        b = np.where(out, 1 + self.bp0 * (s + Lb * (1 - e)), 1 + self.bp0 * y)
        bp = np.where(out, self.bp0 * e, self.bp0)
        # x = a + b u with the physical LE at u_le and the TE at u_te
        g, u0 = self.gmap, self.u_le
        b, bp = b * g, bp * g
        return a - b * u0, ap - bp * u0, b, bp

    def _build_grid(self, n_up, n_down, n_out, nz, x_up, x_down, y_out, z_top, dz1):
        nc, ns = self.nchord, self.nspan
        du = 1.0 / nc
        uw = (np.arange(nc) + 0.5) * du
        up = -0.5 * du - np.cumsum(_stretch(du, x_up - 0.5 * du, n_up))
        dn = 1 + 0.5 * du + np.cumsum(_stretch(du, x_down - 0.5 * du, n_down))
        self.u = np.r_[up[::-1], uw, 1 + 0.5 * du, dn]
        self.iw0 = n_up
        self.iw1 = n_up + nc - 1            # last wing node (wake starts after)
        self.du = du
        dy = self.s / (ns - 0.5)
        yw = np.arange(ns) * dy
        yo = self.s + 0.5 * dy + np.r_[0.0, np.cumsum(_stretch(dy, y_out - 0.5 * dy, n_out - 1))]
        self.eta = np.r_[yw, yo]
        self.dy = dy
        self.zeta = np.r_[0.0, np.cumsum(_stretch(dz1, z_top, nz - 1))]
        self.shape = (self.u.size, self.eta.size, self.zeta.size)
        self.N = int(np.prod(self.shape))
        a, ap, b, bp = self._map(self.eta)
        self.a, self.ap, self.b, self.bp = a, ap, b, bp
        # wing node physical coordinates (nspan, nchord)
        jw = np.arange(ns)
        self.wing_x = (a[jw, None] + b[jw, None] * uw[None, :]) * self.Lref
        self.wing_y = self.eta[jw] * self.Lref
        self.wing_xi = np.broadcast_to((uw - self.u_le) * self.gmap, (ns, nc)).copy()
        self.te_x = (a[jw] + b[jw] * self.u_te) * self.Lref
        wy = np.full(ns, dy)
        wy[0] = 0.5 * dy
        self.strip_width = wy * self.Lref                    # physical strip widths
        self.wing_area = (wy[:, None] * b[jw, None] / self.gmap * du * self.Lref ** 2 * np.ones(nc))

    @property
    def wing_points(self) -> np.ndarray:
        """(nspan*nchord, 2) physical (x, y) of the wing nodes (strip-major)."""
        return np.column_stack([self.wing_x.ravel(),
                                np.repeat(self.wing_y, self.nchord)])

    def grid_summary(self) -> str:
        NU, NY, NZ = self.shape
        return (f"TSD grid {NU} x {NY} x {NZ} = {self.N} nodes; wing {self.nchord} x {self.nspan}; "
                f"x/c in [{self.u[0]:.2f}, {self.u[-1]:.2f}], y/c_r <= {self.eta[-1]:.2f}, "
                f"z/c_r <= {self.zeta[-1]:.2f}, dz1 = {self.zeta[1]:.4f}")

    # ------------------------------------------------------------- operators
    def _operators(self, u, eta, zeta):
        """Finite-volume operators on the (possibly complex-stretched) grid
        coordinates ``u, eta, zeta``; the shear mapping a(y), b(y) is always
        evaluated at the real spanwise stations."""
        NU, NY, NZ = self.shape
        Iu, Iz = sp.identity(NU, format="csr"), sp.identity(NZ, format="csr")
        Iyz = sp.identity(NY * NZ, format="csr")
        Iuy = sp.identity(NU * NY, format="csr")

        def ku(A):
            return sp.kron(A, Iyz, format="csr")

        def ky(A):
            return sp.kron(Iu, sp.kron(A, Iz), format="csr")

        def kz(A):
            return sp.kron(Iuy, A, format="csr")

        o = SimpleNamespace()
        o.Du = ku(_d1_faces(u))
        o.Au = ku(_avg_faces(NU))
        o.Cu = ku(_central(u))
        o.DivU = ku(_div(NU))
        o.DivUs = ku(_div(NU, shift=1))
        Dy = ky(_d1_faces(eta))
        Ay = ky(_avg_faces(NY))
        Cy = ky(_central(eta, zero_first=True))
        DivY = ky(_div(NY))
        Dz = kz(_d1_faces(zeta))
        DivZ = kz(_div(NZ))

        Vu, Vy, Vz = _dual(u), _dual(eta), _dual(zeta)
        ap, b, bp = self.ap, self.b, self.bp
        er = self.eta
        uf = 0.5 * (u[1:] + u[:-1])
        af, apf, bf, bpf = self._map(0.5 * (er[1:] + er[:-1]))
        ones_z = np.ones(NZ)
        # u-face fields (NU-1, NY, NZ)
        S_uf = (ap[None, :] + bp[None, :] * uf[:, None])
        S_uf[:, 0] = 0.0                          # G = phi_y = 0 on the symmetry plane
        S_uf = (S_uf[:, :, None] * ones_z).ravel()
        o.b_uf = ((np.ones(NU - 1)[:, None] * b[None, :])[:, :, None] * ones_z).ravel()
        # eta-face fields (NU, NY-1, NZ)
        S_yf = ((apf[None, :] + bpf[None, :] * u[:, None])[:, :, None] * ones_z).ravel()
        b_yf = ((np.ones(NU)[:, None] * bf[None, :])[:, :, None] * ones_z).ravel()
        # node fields
        o.b_n = (np.ones(NU)[:, None, None] * b[None, :, None] * ones_z).ravel()
        o.Vfac_u = (np.ones(NU)[:, None, None] * Vy[None, :, None] * Vz[None, None, :]).ravel()
        Vfac_y = (Vu[:, None, None] * np.ones(NY)[None, :, None] * Vz[None, None, :]).ravel()
        o.Vfac_z = (Vu[:, None, None] * Vy[None, :, None] * b[None, :, None] * ones_z).ravel()
        o.vol_b = (Vu[:, None, None] * Vy[None, :, None] * Vz[None, None, :]
                   * b[None, :, None]).ravel()
        dg = sp.diags
        Lshear = -dg(o.Vfac_u) @ o.DivU @ dg(S_uf) @ (o.Au @ Cy - dg(S_uf / o.b_uf) @ o.Du)
        Ly = dg(Vfac_y) @ DivY @ (dg(b_yf) @ Dy - dg(S_yf) @ Ay @ o.Cu)
        Lz = dg(o.Vfac_z) @ DivZ @ Dz
        o.Llin = (Lshear + Ly + Lz).tocsr()
        return o

    def _pml_coords(self):
        """Complex-stretched coordinates (perfectly matched layers) upstream,
        outboard and on top: x~ = x - i n sigma0 L xi^3 in the outer fraction
        ``pml_frac`` of each far-field region (n = outward direction)."""
        f, s0 = self.pml_frac, self.pml_sigma

        def layer(x, x_start, x_end):
            L = x_end - x_start
            xi = np.clip((x - x_start) / L, 0.0, None)
            return np.sign(L) * s0 * abs(L) * xi ** 3

        u = self.u.astype(complex)
        u0 = self.u[0]
        u = u + 1j * layer(self.u, (1 - f) * u0, u0)          # upstream: Im(u) grows upstream
        eta = self.eta.astype(complex)
        e1 = self.eta[-1]
        eta = eta - 1j * layer(self.eta, self.s + (1 - f) * (e1 - self.s), e1)
        zeta = self.zeta.astype(complex)
        z1 = self.zeta[-1]
        zeta = zeta - 1j * layer(self.zeta, (1 - f) * z1, z1)
        return u, eta, zeta

    def _unsteady_ops(self):
        if not self.pml:
            return self
        if self._pml_ops is None:
            self._pml_ops = self._operators(*self._pml_coords())
        return self._pml_ops

    def _build_operators(self):
        NU, NY, NZ = self.shape
        o = self._operators(self.u, self.eta, self.zeta)
        for k, v in vars(o).items():
            setattr(self, k, v)
        # masks
        I, J, K = np.meshgrid(np.arange(NU), np.arange(NY), np.arange(NZ), indexing="ij")
        self._I, self._J, self._K = I.ravel(), J.ravel(), K.ravel()
        interior = (I > 0) & (I < NU - 1) & (J < NY - 1) & (K < NZ - 1)
        self.interior = interior.ravel()
        wing = np.zeros((NU, NY), bool)
        wing[self.iw0:self.iw1 + 1, :self.nspan] = True
        self.wing2d = wing
        wake = np.zeros((NU, NY), bool)
        wake[self.iw1 + 1:, :self.nspan] = True
        self.wake2d = wake
        # wing nodes on z = 0, flattened index arrays (nspan, nchord)
        jj, mm = np.meshgrid(np.arange(self.nspan), np.arange(self.nchord), indexing="ij")
        self.wing_idx = ((self.iw0 + mm) * NY + jj) * NZ
        # typical FV diagonal magnitude for boundary-row scaling
        d = np.abs(self.Llin.diagonal())
        self.bc_scale = float(np.median(d[self.interior & (d > 0)]))

    def _solve(self, K, rhs):
        """Solve ``K x = rhs`` (rhs (N,) or (N, nrhs)) for a grid operator.

        ``linear_solver='block'`` (default): exact block-banded Gaussian
        elimination over the chordwise u-planes (see :func:`_block_solve`),
        dense LAPACK kernels.  ``'splu'``: SuperLU with a geometric
        nested-dissection ordering (slower: scipy's SuperLU has no
        optimised BLAS, 3-D fill is large).
        """
        if self.linear_solver == "block":
            P = self.shape[1] * self.shape[2]
            return _block_solve(K, rhs, P, self.shape[0], self.max_bytes)
        if self._perm is None:
            self._perm = _nd_order(self.shape)
        p = self._perm
        Kp = K.tocsc()[p][:, p].tocsc()
        lu = spla.splu(Kp, permc_spec="NATURAL", diag_pivot_thresh=self.pivot_thresh,
                       options=dict(SymmetricMode=True))
        return _PermLU(lu, p).solve(rhs)

    def _node(self, i, j, k):
        NU, NY, NZ = self.shape
        return (np.asarray(i) * NY + np.asarray(j)) * NZ + np.asarray(k)

    def _bc_rows(self, kbar: float, antisym: bool):
        """Boundary-condition rows (sparse) and mask of rows they own."""
        NU, NY, NZ = self.shape
        if self.pml and kbar > 0:      # far-field rows in the stretched coordinates
            u, eta, zeta = self._pml_coords()
        else:
            u, eta, zeta = self.u, self.eta, self.zeta
        sc = self.bc_scale
        I, J, K = self._I, self._J, self._K
        owned = ~self.interior
        rows, cols, vals = [], [], []

        def add(r, c, v):
            rows.append(np.atleast_1d(r)); cols.append(np.atleast_1d(c))   # noqa: E702
            vals.append(np.broadcast_to(np.asarray(v, complex), np.atleast_1d(r).shape))

        done = np.zeros(self.N, bool)
        # z = 0 plane rows for the antisymmetric problem
        if antisym:
            on0 = (K == 0) & ~done
            wake0 = on0 & self.wake2d[I, J]
            off0 = on0 & ~self.wing2d[I, J] & ~self.wake2d[I, J]
            n = np.where(off0)[0]
            add(n, n, sc)
            n = np.where(wake0)[0]
            jn = J[n]
            # recursive form of phi(x) = phi_TE exp(-i kbar (x - x_te)): couples
            # only neighbouring u-planes (keeps the operator block-tridiagonal)
            ph = np.exp(-1j * kbar * self.b[jn] * (u[I[n]] - u[I[n] - 1]))
            add(n, n, sc)
            add(n, self._node(I[n] - 1, jn, 0), -sc * ph)
            done |= off0 | wake0
            owned = owned | off0 | wake0
        # upstream: phi = 0
        n = np.where((I == 0) & ~done)[0]
        add(n, n, sc)
        done[n] = True
        # downstream: phi_x + i kbar phi = 0
        n = np.where((I == NU - 1) & ~done)[0]
        h = (u[-1] - u[-2]) * self.b[J[n]]
        m = self._node(NU - 2, J[n], K[n])
        add(n, n, sc * (1 + 0.5j * kbar * h))
        add(n, m, sc * (-1 + 0.5j * kbar * h))
        done[n] = True
        # top: d phi/dz + phi z / r^2 = 0
        n = np.where((K == NZ - 1) & ~done)[0]
        hz = zeta[-1] - zeta[-2]
        r2 = eta[J[n]] ** 2 + zeta[-1] ** 2
        add(n, n, sc * (1 + hz * zeta[-1] / r2))
        add(n, self._node(I[n], J[n], NZ - 2), -sc * np.ones(n.size))
        done[n] = True
        # outboard: d phi/dy + phi y / r^2 = 0
        n = np.where((J == NY - 1) & ~done)[0]
        hy = eta[-1] - eta[-2]
        r2 = eta[-1] ** 2 + zeta[K[n]] ** 2
        add(n, n, sc * (1 + hy * eta[-1] / r2))
        add(n, self._node(I[n], NY - 2, K[n]), -sc * np.ones(n.size))
        done[n] = True
        assert np.all(done[owned])
        B = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.N, self.N))
        fv = ~owned
        return B, fv

    def _switched_div(self, mu):
        """Murman fully conservative switched divergence ``DivU_sw``."""
        NU, NY, NZ = self.shape
        mu3 = mu.reshape(self.shape)
        mus = np.zeros_like(mu3)
        mus[1:] = mu3[:-1]
        return (sp.diags(1.0 - mu) @ self.DivU + sp.diags(mus.ravel()) @ self.DivUs).tocsr()

    # ---------------------------------------------------------------- steady
    def _thickness_slope(self, xi):
        t = self.thickness
        if t is None:
            return np.zeros_like(xi)
        if hasattr(t, "slope"):
            return np.asarray(t.slope(xi), float)
        h = 1e-6
        return 0.25 * (np.asarray(t(xi + h)) - np.asarray(t(xi - h))) / h

    def _flux(self, p):
        M2, g = self.mach ** 2, self.gamma
        return (1 - M2) * p - 0.5 * (g + 1) * M2 * p * p, (1 - M2) - (g + 1) * M2 * p

    def _node_A(self, phi):
        M2, g = self.mach ** 2, self.gamma
        px = (self.Cu @ phi) / self.b_n
        return 1 - M2 - (g + 1) * M2 * px

    def steady(self, tol: float = 1e-9, maxiter: int = 40, verbose: bool | None = None,
               freeze_after: int = 25) -> SteadyResult:
        """Nonlinear steady thickness (symmetric) flow at zero incidence."""
        verbose = self.verbose if verbose is None else verbose
        B, fv = self._bc_rows(0.0, antisym=False)
        B = B.real.tocsr()
        Pfv = sp.diags(fv.astype(float))
        rhs = np.zeros(self.N)
        slope = self._thickness_slope(self.wing_xi)
        rhs[self.wing_idx.ravel()] = self.Vfac_z[self.wing_idx.ravel()] * slope.ravel()
        phi = np.zeros(self.N)
        interior = self.interior

        def parts(phi, mu):
            p = (self.Du @ phi) / self.b_uf
            F, Af = self._flux(p)
            Dsw = self._switched_div(mu)
            R = fv * (self.Vfac_u * (Dsw @ F) + self.Llin @ phi) + B @ phi - fv * rhs
            return R, Dsw, Af

        def switches(phi):
            return ((self._node_A(phi) < 0) & interior).astype(float)

        mu = switches(phi)
        R, Dsw, Af = parts(phi, mu)
        res0 = max(np.linalg.norm(rhs), 1e-300)
        hist = [np.linalg.norm(R) / res0]
        converged = hist[-1] < tol
        it = 0
        t0 = time.perf_counter()
        while not converged and it < maxiter:
            it += 1
            Jx = sp.diags(self.Vfac_u) @ Dsw @ sp.diags(Af / self.b_uf) @ self.Du
            Jm = (Pfv @ (Jx + self.Llin) + B).tocsc()
            dphi = self._solve(Jm, -R)
            # back-tracking on the residual (switches recomputed unless frozen)
            lam, rn_old = 1.0, np.linalg.norm(R)
            for _ in range(8):
                trial = phi + lam * dphi
                mu_t = mu if it > freeze_after else switches(trial)
                Rt, Dt, At = parts(trial, mu_t)
                if np.linalg.norm(Rt) < rn_old or lam < 0.05:
                    break
                lam *= 0.5
            phi, mu, R, Dsw, Af = trial, mu_t, Rt, Dt, At
            hist.append(np.linalg.norm(R) / res0)
            if verbose:
                print(f"  TSD steady it {it:2d}: |R|/|R0| = {hist[-1]:.3e}  lam = {lam:.3f}  "
                      f"supersonic = {int(mu.sum())}  ({time.perf_counter() - t0:.1f} s)")
            converged = hist[-1] < tol or (np.linalg.norm(lam * dphi) < tol * max(1.0, np.abs(phi).max())
                                           and hist[-1] < 1e-5)
        phi3 = phi.reshape(self.shape)
        px = ((self.Cu @ phi) / self.b_n).reshape(self.shape)
        cp = -2 * px[self.iw0:self.iw1 + 1, :self.nspan, 0].T
        res = SteadyResult(phi3, cp, self.wing_x.copy(), self.wing_y.copy(), self.mach,
                           int(mu.sum()), bool(converged), it, hist,
                           cp_critical(self.mach, self.gamma, tsd=True))
        self._mean = res
        return res

    @property
    def mean_flow(self) -> SteadyResult:
        if self._mean is None:
            if self.thickness is None:
                z = np.zeros(self.shape)
                self._mean = SteadyResult(z, np.zeros((self.nspan, self.nchord)), self.wing_x.copy(),
                                          self.wing_y.copy(), self.mach, 0, True, 0, [],
                                          cp_critical(self.mach, self.gamma, tsd=True) if self.mach > 0 else np.nan)
            else:
                self.steady()
        return self._mean

    # -------------------------------------------------------------- unsteady
    def unsteady_matrix(self, kbar: float):
        """Assembled time-linearised (antisymmetric) operator at ``kbar``;
        returns (K, fv) where ``fv`` marks finite-volume rows."""
        phi0 = self.mean_flow.phi.ravel()
        M2 = self.mach ** 2
        p = (self.Du @ phi0) / self.b_uf
        _, Af = self._flux(p)
        mu = ((self._node_A(phi0) < 0) & self.interior).astype(float)
        Dsw = self._switched_div(mu)
        B, fv = self._bc_rows(kbar, antisym=True)
        o = self._unsteady_ops() if kbar > 0 else self
        X = sp.diags(o.Vfac_u) @ Dsw @ (sp.diags(Af / o.b_uf) @ o.Du
                                       - 2j * kbar * M2 * o.Au)
        Kmat = sp.diags(fv.astype(float)) @ (X + o.Llin + sp.diags(M2 * kbar ** 2 * o.vol_b)) + B
        Kmat = Kmat.tocsr()
        Kmat.eliminate_zeros()
        return Kmat, fv

    def _eval_modes(self, h_funcs, x, y):
        H, HX = [], []
        for f in h_funcs:
            h, hx = f(x, y)
            H.append(np.broadcast_to(np.asarray(h, float), np.shape(x)))
            HX.append(np.broadcast_to(np.asarray(hx, float), np.shape(x)))
        return np.stack(H, -1), np.stack(HX, -1)

    def solve_normalwash(self, kbar: float, w: np.ndarray) -> np.ndarray:
        """Solve the antisymmetric problem for wing normalwash ``w/U``
        (nspan*nchord, nrhs) complex; returns phi (N, nrhs)."""
        Kmat, fv = self.unsteady_matrix(kbar)
        w = np.asarray(w, complex).reshape(self.nspan * self.nchord, -1)
        rhs = np.zeros((self.N, w.shape[1]), complex)
        idx = self.wing_idx.ravel()
        rhs[idx] = self.Vfac_z[idx, None] * w
        rhs *= fv[:, None]
        t0 = time.perf_counter()
        phi = self._solve(Kmat, rhs)
        if self.verbose:
            print(f"  TSD unsteady solve kbar={kbar:.4f}: {time.perf_counter() - t0:.2f} s "
                  f"({w.shape[1]} rhs)")
        return phi

    def _surface_phi(self, phi):
        """Upper-surface potential on z = 0 along wing+wake rows: (NU, nspan, nrhs)."""
        NU, NY, NZ = self.shape
        return phi.reshape(NU, NY, NZ, -1)[:, :self.nspan, 0, :]

    def dcp_from_phi(self, kbar, phi):
        """Pressure jump ``dCp = 4 (phi_x + i kbar phi)`` at the wing nodes,
        (nspan, nchord, nrhs)."""
        ps = self._surface_phi(phi)
        u = self.u
        i = np.arange(self.iw0, self.iw1 + 1)
        dphidu = (ps[i + 1] - ps[i - 1]) / (u[i + 1] - u[i - 1])[:, None, None]
        b = self.b[:self.nspan][None, :, None]
        dcp = 4 * (dphidu / b + 1j * kbar * ps[i])
        return np.transpose(dcp, (1, 0, 2))

    def phi_te(self, kbar, phi):
        """Trailing-edge upper potential per strip (nspan, nrhs): the last wing
        node value convected to the TE, ``phi_last exp(-i kbar b (u_te - u_last))``."""
        ps = self._surface_phi(phi)
        d = self.u_te - self.u[self.iw1]
        return ps[self.iw1] * np.exp(-1j * kbar * self.b[:self.nspan] * d)[:, None]

    def unsteady(self, k: float, bref: float, h_func, return_phi: bool = False):
        """Pressure jump for one mode (or a list of modes).

        ``h_func(x, y) -> (h, dh/dx)`` at physical wing points.  Returns
        ``dCp`` with shape (nspan, nchord) (or (nspan, nchord, nm) for a list).
        """
        single = callable(h_func)
        funcs = [h_func] if single else list(h_func)
        kbar = float(k) * self.Lref / bref
        H, HX = self._eval_modes(funcs, self.wing_x, self.wing_y[:, None] * np.ones(self.nchord))
        w = HX + 1j * (kbar / self.Lref) * H
        phi = self.solve_normalwash(kbar, w.reshape(-1, len(funcs)))
        dcp = self.dcp_from_phi(kbar, phi)
        out = dcp[..., 0] if single else dcp
        return (out, phi) if return_phi else out

    def generalized_forces(self, k: float, bref: float, mode_funcs, gust: bool = False,
                           x_gust0: float | None = None, integration: str = "parts"):
        """GAF matrix ``Q`` (nm, nm) [+ gust column] and dCp (nspan, nchord, nm[+1]).

        ``Q_ij = int int h_i dCp_j dA`` over the half wing (physical units,
        same convention as :class:`flaap.aero.gaf.AeroModel`).
        """
        funcs = list(mode_funcs)
        nm = len(funcs)
        kbar = float(k) * self.Lref / bref
        ns, nc = self.nspan, self.nchord
        Y = self.wing_y[:, None] * np.ones(nc)
        H, HX = self._eval_modes(funcs, self.wing_x, Y)          # (ns, nc, nm)
        W = HX + 1j * (kbar / self.Lref) * H
        cols = [W.reshape(-1, nm)]
        if gust:
            x0 = float(self.wing_x.min()) if x_gust0 is None else x_gust0
            cols.append((-np.exp(-1j * (kbar / self.Lref) * (self.wing_x - x0))).reshape(-1, 1))
        phi = self.solve_normalwash(kbar, np.hstack(cols))
        dcp = self.dcp_from_phi(kbar, phi)                       # (ns, nc, nrhs)
        if integration == "direct":
            Qall = np.einsum("jmi,jmr->ir", H * self.wing_area[..., None], dcp)
        else:
            Hte, _ = self._eval_modes(funcs, self.te_x, self.wing_y)   # (ns, nm)
            ps = self._surface_phi(phi)[self.iw0:self.iw1 + 1]          # (nc, ns, nrhs)
            ps = np.transpose(ps, (1, 0, 2))                             # (ns, nc, nrhs)
            pte = self.phi_te(kbar, phi)                                 # (ns, nrhs)
            dx = (self.b[:ns] * self.du * self.Lref)[:, None]            # physical dx
            wts = np.ones(nc)
            wts[-1] = (self.u_te - (1.0 - self.du)) / self.du          # cell cut at the TE
            G = (HX - 1j * (kbar / self.Lref) * H) * (dx * wts)[..., None]  # (ns, nc, nm)
            per_strip = 4 * self.Lref * (Hte[:, :, None] * pte[:, None, :]
                             - np.einsum("jmi,jmr->jir", G, ps))         # (ns, nm, nrhs)
            Qall = np.einsum("j,jir->ir", self.strip_width, per_strip)
        Q = Qall[:, :nm]
        Qg = Qall[:, nm] if gust else None
        return Q, Qg, dcp

    def cp_lower_upper(self, kbar, phi):
        """Unsteady upper-surface Cp1 = -dCp/2 at the wing nodes."""
        return -0.5 * self.dcp_from_phi(kbar, phi)


# ---------------------------------------------------------------- GAF table
def tsd_gaf_table(surface: Surface, mach: float, kvals, bref: float, mode_funcs,
                  names: Sequence[str] | None = None, thickness: Callable | None = None,
                  solver: TSDSolver | None = None, gust: bool = False,
                  x_gust0: float | None = None, keep_cp: bool = True,
                  verbose: bool = False, **grid_kw):
    """Tabulate TSD generalized aerodynamic forces as a :class:`GAFTable`.

    The steady thickness mean flow is computed once (if ``thickness`` is
    given), then for each ``k`` one block elimination of the time-linearised
    operator serves all modes (multiple right-hand sides).  ``cp`` holds dCp on the TSD wing grid (nk, nspan*nchord, nm).
    """
    from .gaf import GAFTable
    kvals = np.atleast_1d(np.asarray(kvals, float))
    if solver is None:
        solver = TSDSolver(surface, mach, thickness, verbose=verbose, **grid_kw)
    if solver.thickness is not None and solver._mean is None:
        t0 = time.perf_counter()
        st = solver.steady(verbose=verbose)
        if verbose:
            print(f"  TSD mean flow: {st.iterations} it, converged={st.converged}, "
                  f"supersonic pts={st.n_supersonic}, {time.perf_counter() - t0:.1f} s")
    funcs = list(mode_funcs)
    nm = len(funcs)
    names = list(names) if names is not None else [f"mode{i + 1}" for i in range(nm)]
    nw = solver.nspan * solver.nchord
    Q = np.zeros((kvals.size, nm, nm), complex)
    Qg = np.zeros((kvals.size, nm), complex) if gust else None
    cp = np.zeros((kvals.size, nw, nm), complex) if keep_cp else None
    cpg = np.zeros((kvals.size, nw), complex) if (keep_cp and gust) else None
    for i, k in enumerate(kvals):
        t0 = time.perf_counter()
        q, qg, dcp = solver.generalized_forces(k, bref, funcs, gust=gust, x_gust0=x_gust0)
        Q[i] = q
        if gust:
            Qg[i] = qg
        if keep_cp:
            cp[i] = dcp[..., :nm].reshape(nw, nm)
            if gust:
                cpg[i] = dcp[..., nm].reshape(nw)
        if verbose:
            print(f"  TSD GAF k={k:.4f} done ({time.perf_counter() - t0:.1f} s)")
    tab = GAFTable(kvals, Q, Qg, mach, bref, names, cp, cpg)
    tab.tsd_solver = solver
    return tab
