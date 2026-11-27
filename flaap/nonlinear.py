"""Control-surface freeplay and other concentrated structural nonlinearities.

The module adds a single nonlinear spring to a linear FLAAP aeroelastic
model (:class:`flaap.timedomain.AeroelasticSystem`) and predicts the
resulting limit-cycle oscillations (LCO) in two ways:

1. **Time marching** of the RFA state space with the nonlinear restoring
   moment as a feedback force.
2. **Describing-function** (harmonic-balance) analysis, in which the
   nonlinear spring is replaced by its amplitude-dependent equivalent
   stiffness ``K_eq(A)`` and the linear flutter problem is solved for the
   stiffness.

Model
-----
The structural coordinates ``eta`` (``ns``) of the linear system contain
an elastic hinge rotation ``delta = d . eta``.  For the pitch-plunge-flap
section ``d = e_beta`` (see :mod:`flaap.models.typical_section`).  For a
3-D wing whose modal basis is augmented with a rigid control rotation
mode (:func:`flaap.aero.gaf.control_aero_mode`) with hinge stiffness
``K_h``, hinge inertia ``I_h`` and inertial coupling to the elastic modes,
``d`` is the unit vector of that mode (see :mod:`flaap.models.goland_flap`).
The linear system contains a linear hinge spring ``K_lin d d^T`` in ``K``
(``K_lin = d^T K d`` by default). This spring is removed and replaced by
the nonlinear restoring moment ``M_h(delta)``::

    M eta'' + C eta' + (K - K_lin d d^T) eta + d M_h(d . eta) = q Q(p) xi + ...

With the RFA (Roger) state space ``x = [eta, eta', lag states]`` of
:meth:`AeroelasticSystem.state_space` (matrices ``A_lin``, ``B``) this is::

    x' = A_free x - b_h M_h(c . x) + B u,       A_free = A_lin + K_lin b_h c^T
    b_h = [0; Mbar^-1 d; 0],  c = [d; 0; 0],    Mbar = M - (rho b^2 / 2) A_2

Nonlinearities (restoring hinge moment ``M_h``, positive = trailing edge down)
-------------------------------------------------------------------------------
* :class:`Freeplay` ``(K, delta0, offset=0)``: zero stiffness for
  ``|delta - offset| < delta0`` and ``K`` outside, continuous.
  :func:`preloaded_freeplay` gives the gap shifted by a preload angle.
* :class:`Bilinear` ``(K, delta0, K_inner)``: stiffness ``K_inner``
  inside the gap and ``K`` outside.  ``K_inner = 0`` is freeplay.
* :class:`PiecewiseLinear`: the general continuous piecewise-linear
  spring that the three above derive from.
* :class:`Cubic` ``(K, K3)``: ``M = K delta + K3 delta^3``, hardening
  for ``K3 > 0``.

Time marching (:meth:`NonlinearAeroelastic.simulate`)
-----------------------------------------------------
* ``method="exact"`` (default for piecewise-linear springs): the system
  is linear and time-invariant inside every segment of the spring
  characteristic.  The solution is advanced with the exact matrix
  exponential of the augmented affine system.  Switching instants
  (``delta`` reaching a corner) are located by Newton iteration on the
  exact solution, started from a cubic Hermite estimate.  No integration
  error is introduced apart from the root tolerance, and the time step
  only sets the output sampling.
* ``method="events"``: ``scipy.integrate.solve_ivp`` (RK45 / DOP853 /
  ...) per segment with terminal events at the corners.  This method
  accepts arbitrary inputs ``u(t)``, such as gusts.
* ``method="smooth"``: a single ``solve_ivp`` run on a smoothed
  characteristic. The ramp ``max(0, x)`` is replaced by
  ``(x + sqrt(x^2 + eps^2)) / 2``.  This is the method for smooth
  springs (cubic).

LCO post-processing (:func:`lco_from_history`) uses the last part of
the history.  It returns the amplitude (half peak-to-peak and
first-harmonic), the mean, the frequency from mean-crossings and a
classification: ``'decay'``, ``'lco'``, ``'growing'`` or
``'irregular'``.  :func:`bifurcation_sweep` performs velocity
continuation (the end state of one speed is the initial condition of
the next).

Describing function
-------------------
For a sinusoidal hinge motion ``delta = B + A sin(theta)``, the
first-harmonic (dual-input) describing function is::

    N1(A, B) = 1/(pi A) int_0^{2pi} M_h(B + A sin th) sin th d th    (equivalent stiffness K_eq)
    N0(A, B) = 1/(2 pi)  int_0^{2pi} M_h(B + A sin th) d th          (mean moment)

For a continuous piecewise-linear spring ``M_h = c0 + s0 x +
sum Ds_i ramp(x - b_i)``, both are analytic.  With ``s = (b_i - B)/A``
and ``|s| < 1``, every ramp contributes::

    N1 += Ds_i (1/pi) (pi/2 - asin s - s sqrt(1 - s^2))
    N0 += Ds_i (A/pi) (sqrt(1 - s^2) - s (pi/2 - asin s))

For freeplay with ``B = 0`` this is the classical
``K_eq/K = 1 - (2/pi) (asin(d0/A) + (d0/A) sqrt(1 - (d0/A)^2))``.

LCO condition (:func:`describing_function_lco`): at every speed ``V``,
find the hinge stiffnesses ``K*`` for which the linear system
``A_free - K* b_h c^T`` has a neutrally stable oscillatory root, i.e. the
inverse of the flutter-speed curve ``V_f(K)``.  The LCO amplitude then
solves ``N1(A) = K*``.  An LCO is stable when increasing the amplitude
stabilises the equivalent linear system,
``d sigma / dA = (d sigma / dK)(dN1/dA) < 0``.  For freeplay, ``N1``
is independent of ``delta0`` when written in terms of ``A / delta0``.
Freeplay LCO amplitudes are therefore exactly proportional to
``delta0``, and no LCO exists below the zero-stiffness flutter speed
``V_f(K = 0)``.  Above the nominal flutter speed ``V_f(K)`` the motion
is unbounded.

Conventions follow FLAAP: SI units, ``e^{i w t}``, ``k = w b / U``,
control rotations positive trailing-edge down, moment ``M_h`` restoring
(acting as ``-M_h`` on the flap).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

__all__ = [
    "HingeNonlinearity", "PiecewiseLinear", "Freeplay", "Bilinear", "Cubic",
    "preloaded_freeplay", "NonlinearAeroelastic", "SimulationResult", "LCOInfo",
    "lco_from_history", "bifurcation_sweep", "describing_function_lco", "DFResult",
    "flutter_speed_vs_stiffness", "freeplay_keq", "hinge_compliance",
]


_trapz = getattr(np, "trapezoid", None) or np.trapz


# ======================================================================
# nonlinear springs
# ======================================================================
def _ramp(x):
    return np.maximum(x, 0.0)


def _ramp_smooth(x, eps):
    return 0.5 * (x + np.sqrt(x * x + eps * eps))


def _ramp_df(A, b):
    """Describing functions ``(N1, N0)`` of ``ramp(A sin th - b)``."""
    A = np.asarray(A, float)
    b = np.asarray(b, float)
    A_, b_ = np.broadcast_arrays(A, b)
    N1 = np.zeros(A_.shape)
    N0 = np.zeros(A_.shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(A_ > 0, b_ / np.where(A_ > 0, A_, 1.0), np.sign(b_) * np.inf)
    inside = np.abs(s) < 1
    below = s <= -1                       # ramp always active
    si = np.clip(s, -1, 1)
    c = np.sqrt(1 - si * si)
    th = 0.5 * np.pi - np.arcsin(si)
    N1 = np.where(inside, (th - si * c) / np.pi, np.where(below, 1.0, 0.0))
    N0 = np.where(inside, A_ * (c - si * th) / np.pi, np.where(below, -b_, 0.0))
    return N1, N0


class HingeNonlinearity:
    """Base class of a nonlinear restoring hinge moment ``M_h(delta)``.

    ``K`` is the reference (nominal) linear stiffness; describing-function
    results are usually reported as ``K_eq / K``.
    """

    K: float = 1.0
    piecewise_linear: bool = False

    def moment(self, x):                                  # pragma: no cover - abstract
        raise NotImplementedError

    def smooth_moment(self, x, eps: float | None = None):
        return self.moment(x)

    def stiffness(self, x, h: float = 1e-7):
        """Local tangent stiffness ``dM/d delta``."""
        return (self.moment(np.asarray(x) + h) - self.moment(np.asarray(x) - h)) / (2 * h)

    @property
    def is_odd(self) -> bool:
        x = np.linspace(0.01, 1.0, 7) * self.scale
        return bool(np.allclose(self.moment(-x), -self.moment(x), rtol=1e-9, atol=1e-12 * self.K))

    @property
    def scale(self) -> float:
        """Characteristic angle (gap size) used for default tolerances."""
        return 0.01

    def describing_function(self, A, B=0.0, n: int = 256):
        """``(N1, N0)`` for ``delta = B + A sin(th)`` (numerical quadrature;
        analytic in subclasses)."""
        th = (np.arange(n) + 0.5) * 2 * np.pi / n
        A = np.atleast_1d(np.asarray(A, float))
        B = np.broadcast_to(np.asarray(B, float), A.shape)
        x = B[..., None] + A[..., None] * np.sin(th)
        Mx = self.moment(x)
        N0 = Mx.mean(-1)
        with np.errstate(divide="ignore", invalid="ignore"):
            N1 = np.where(A > 0, 2 * (Mx * np.sin(th)).mean(-1) / np.where(A > 0, A, 1), self.stiffness(B))
        return N1, N0

    def keq(self, A, B=0.0):
        """Equivalent (first-harmonic) stiffness ``N1(A, B)``."""
        return self.describing_function(A, B)[0]


class PiecewiseLinear(HingeNonlinearity):
    """Continuous piecewise-linear spring.

    ``breaks`` (sorted, length n) are the corner angles and ``slopes``
    (length n+1) are the stiffnesses of the segments.  ``point = (x_ref,
    M_ref)`` fixes the level.  Written as
    ``M(x) = c0 + s0 x + sum_i (s_i - s_{i-1}) ramp(x - b_i)``.
    """

    piecewise_linear = True

    def __init__(self, breaks, slopes, point=(0.0, 0.0), K: float | None = None):
        self.breaks = np.atleast_1d(np.asarray(breaks, float))
        self.slopes = np.atleast_1d(np.asarray(slopes, float))
        if self.slopes.size != self.breaks.size + 1:
            raise ValueError("need len(slopes) == len(breaks) + 1")
        if np.any(np.diff(self.breaks) <= 0):
            raise ValueError("breaks must be strictly increasing")
        self.dslope = np.diff(self.slopes)
        xr, mr = point
        self.c0 = mr - self.slopes[0] * xr - np.sum(self.dslope * _ramp(xr - self.breaks))
        self.K = float(np.max(np.abs(self.slopes))) if K is None else float(K)

    @property
    def scale(self) -> float:
        if self.breaks.size > 1:
            return float(np.min(np.diff(self.breaks)) / 2)
        return max(abs(float(self.breaks[0])), 1e-3) if self.breaks.size else 0.01

    def moment(self, x):
        x = np.asarray(x, float)
        return self.c0 + self.slopes[0] * x + np.sum(
            self.dslope * _ramp(x[..., None] - self.breaks), axis=-1)

    def smooth_moment(self, x, eps: float | None = None):
        eps = 0.02 * self.scale if eps is None else eps
        x = np.asarray(x, float)
        r = _ramp_smooth(x[..., None] - self.breaks, eps) - _ramp_smooth(-self.breaks, eps) + _ramp(-self.breaks)
        return self.c0 + self.slopes[0] * x + np.sum(self.dslope * r, axis=-1)

    def segments(self):
        """List of ``(lo, hi, slope, intercept)`` with ``M = slope x + intercept``."""
        edges = np.r_[-np.inf, self.breaks, np.inf]
        out = []
        for j, s in enumerate(self.slopes):
            if np.isfinite(edges[j]):
                xr = edges[j]
            elif np.isfinite(edges[j + 1]):
                xr = edges[j + 1]
            else:
                xr = 0.0
            m = float(self.moment(xr)) - s * xr
            out.append((edges[j], edges[j + 1], float(s), m))
        return out

    def segment_of(self, x, xdot=0.0) -> int:
        """Segment index of ``x``; at a corner the direction ``xdot`` decides."""
        j = int(np.searchsorted(self.breaks, x, side="right"))
        tol = 1e-12 * max(1.0, self.scale)
        for i, b in enumerate(self.breaks):
            if abs(x - b) <= tol:
                return i + 1 if xdot > 0 else i
        return j

    def describing_function(self, A, B=0.0, n=None):
        A = np.atleast_1d(np.asarray(A, float))
        B = np.broadcast_to(np.asarray(B, float), A.shape)
        N1 = np.full(A.shape, self.slopes[0])
        N0 = self.c0 + self.slopes[0] * B
        for b, ds in zip(self.breaks, self.dslope):
            n1, n0 = _ramp_df(A, b - B)
            # A = 0: tangent stiffness / static value
            n1 = np.where(A > 0, n1, (B > b).astype(float))
            n0 = np.where(A > 0, n0, _ramp(B - b))
            N1 = N1 + ds * n1
            N0 = N0 + ds * n0
        return N1, N0


class Freeplay(PiecewiseLinear):
    """Freeplay ``+-delta0`` (rad) around ``offset`` with stiffness ``K`` outside.

    ``M = 0`` for ``|x - offset| <= delta0`` and ``K (x - offset -+ delta0)``
    outside.
    """

    def __init__(self, K: float, delta0: float, offset: float = 0.0):
        self.delta0 = float(delta0)
        self.offset = float(offset)
        super().__init__([offset - delta0, offset + delta0], [K, 0.0, K], (offset, 0.0), K)

    def __repr__(self):
        return f"Freeplay(K={self.K:.6g}, delta0={self.delta0:.6g}, offset={self.offset:.6g})"


class Bilinear(PiecewiseLinear):
    """Bilinear spring: ``K_inner`` for ``|x| < delta0``, ``K`` outside."""

    def __init__(self, K: float, delta0: float, K_inner: float = 0.0):
        self.delta0 = float(delta0)
        self.K_inner = float(K_inner)
        super().__init__([-delta0, delta0], [K, K_inner, K], (0.0, 0.0), K)

    def __repr__(self):
        return f"Bilinear(K={self.K:.6g}, delta0={self.delta0:.6g}, K_inner={self.K_inner:.6g})"


def preloaded_freeplay(K: float, delta0: float, preload: float) -> Freeplay:
    """Freeplay whose gap is centred at ``preload`` (rad) instead of 0.

    This models a hinge whose spring is preloaded (Conner et al. 1997):
    at rest the flap sits on the stiff branch if ``|preload| > delta0``,
    and the freeplay is reached only by motions with amplitude above
    about ``|preload| - delta0``.  The characteristic is not odd, so the
    LCO has a mean (bias) component and the describing-function analysis
    of :func:`describing_function_lco` (which assumes ``B = 0``) does not
    apply.  Use time marching.
    """
    return Freeplay(K, delta0, offset=preload)


class Cubic(HingeNonlinearity):
    """Cubic spring ``M = K x + K3 x^3`` (hardening for ``K3 > 0``)."""

    def __init__(self, K: float, K3: float, scale: float = 0.01):
        self.K = float(K)
        self.K3 = float(K3)
        self._scale = scale

    @property
    def scale(self):
        return self._scale

    def moment(self, x):
        x = np.asarray(x, float)
        return self.K * x + self.K3 * x ** 3

    def describing_function(self, A, B=0.0, n=None):
        A = np.atleast_1d(np.asarray(A, float))
        B = np.broadcast_to(np.asarray(B, float), A.shape)
        N1 = self.K + self.K3 * (3 * B * B + 0.75 * A * A)
        N0 = self.K * B + self.K3 * (B ** 3 + 1.5 * B * A * A)
        return N1, N0

    def __repr__(self):
        return f"Cubic(K={self.K:.6g}, K3={self.K3:.6g})"


def freeplay_keq(A_over_delta0):
    """``K_eq / K`` of symmetric freeplay versus ``A / delta0`` (0 for A <= delta0)."""
    r = np.asarray(A_over_delta0, float)
    s = np.clip(1.0 / np.maximum(r, 1e-300), 0, 1)
    return np.where(r > 1, 1 - 2 / np.pi * (np.arcsin(s) + s * np.sqrt(1 - s * s)), 0.0)


# ======================================================================
# results
# ======================================================================
@dataclass
class LCOInfo:
    """Steady-state characteristics of a (window of a) time history."""

    amplitude: float          # half peak-to-peak in the window
    amplitude_h1: float       # first-harmonic amplitude at the detected frequency
    mean: float
    frequency: float          # Hz (nan if no oscillation)
    kind: str                 # 'decay' | 'lco' | 'growing' | 'irregular'
    growth: float             # amplitude ratio (last half / first half of the window)
    peaks: np.ndarray = field(default_factory=lambda: np.zeros(0))   # local maxima in the window
    period_scatter: float = 0.0

    @property
    def omega(self) -> float:
        return 2 * np.pi * self.frequency


def lco_from_history(t, y, window: float | None = None, nwindow_periods: int = 12,
                     tol_lco: float = 0.03, amp_floor: float = 0.0) -> LCOInfo:
    """LCO amplitude / frequency from the tail of a time history.

    ``window`` (s) defaults to the last 40 % of the record.  The window is
    split in two halves.  ``growth`` is the ratio of their peak-to-peak
    amplitudes, and the history is classified as:

    * ``'decay'`` if the amplitude is below ``amp_floor`` (or at round-off
      level, ``1e-8 max|y|``) or ``growth < 1 - tol_lco``;
    * ``'growing'`` if ``growth > 1 + tol_lco``;
    * ``'irregular'`` if the period scatter is above 5 % (quasi-periodic or
      chaotic, or several peaks per cycle);
    * ``'lco'`` otherwise.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    T = t[-1] - t[0]
    w = 0.4 * T if window is None else min(window, T)
    sel = t >= t[-1] - w
    tw, yw = t[sel], y[sel]
    mean = float(0.5 * (yw.max() + yw.min()))
    amp = float(0.5 * (yw.max() - yw.min()))
    half = tw.size // 2
    a1 = 0.5 * (yw[:half].max() - yw[:half].min()) if half > 2 else amp
    a2 = 0.5 * (yw[half:].max() - yw[half:].min()) if half > 2 else amp
    growth = float(a2 / a1) if a1 > 0 else (1.0 if a2 == 0 else np.inf)
    # mean crossings (upward) for the period
    z = yw - np.mean(yw)
    idx = np.nonzero((z[:-1] < 0) & (z[1:] >= 0))[0]
    if idx.size >= 3:
        tc = tw[idx] - z[idx] * (tw[idx + 1] - tw[idx]) / (z[idx + 1] - z[idx])
        P = np.diff(tc)
        per = float(np.mean(P))
        scatter = float(np.std(P) / per)
        freq = 1.0 / per
        # first harmonic over an integer number of periods
        t0, t1 = tc[0], tc[-1]
        m = (tw >= t0) & (tw <= t1)
        wv = 2 * np.pi * freq
        if m.sum() > 8:
            tt, yy = tw[m], yw[m] - np.mean(yw[m])
            cs = 2 * _trapz(yy * np.cos(wv * tt), tt) / (t1 - t0)
            sn = 2 * _trapz(yy * np.sin(wv * tt), tt) / (t1 - t0)
            a_h1 = float(np.hypot(cs, sn))
        else:
            a_h1 = amp
    else:
        freq, scatter, a_h1 = np.nan, 0.0, amp
    pk = np.nonzero((yw[1:-1] > yw[:-2]) & (yw[1:-1] >= yw[2:]))[0] + 1
    peaks = yw[pk]
    # oscillation at round-off level about a (possibly non-zero) static state
    tiny = 1e-8 * float(np.max(np.abs(y))) if y.size else 0.0
    if amp <= max(amp_floor, tiny) or growth < 1 - tol_lco or not np.isfinite(freq):
        kind = "decay"
    elif growth > 1 + tol_lco:
        kind = "growing"
    elif scatter > 0.05:
        kind = "irregular"
    else:
        kind = "lco"
    return LCOInfo(amp, a_h1, mean, freq, kind, growth, peaks, scatter)


@dataclass
class SimulationResult:
    """Time history of a nonlinear simulation."""

    t: np.ndarray
    x: np.ndarray                 # (nt, n) full state
    eta: np.ndarray               # (nt, ns) structural coordinates
    etad: np.ndarray              # (nt, ns)
    delta: np.ndarray             # (nt,) hinge rotation d . eta
    names: list[str]
    V: float
    method: str
    nswitch: int = 0
    y: np.ndarray | None = None   # (nt, nout) outputs of the linear state space (if requested)
    outputs: list[str] = field(default_factory=list)

    def dof(self, name: str) -> np.ndarray:
        if name in ("delta", "hinge"):
            return self.delta
        return self.eta[:, self.names.index(name)]

    def rate(self, name: str) -> np.ndarray:
        return self.etad[:, self.names.index(name)]

    def output(self, label: str) -> np.ndarray:
        return self.y[:, self.outputs.index(label)]

    def lco(self, name: str = "delta", **kw) -> LCOInfo:
        return lco_from_history(self.t, self.dof(name), **kw)


# ======================================================================
# model
# ======================================================================
class NonlinearAeroelastic:
    """Linear aeroelastic system with one nonlinear hinge spring.

    Parameters
    ----------
    system : AeroelasticSystem
        Linear model whose structural coordinates contain the hinge
        rotation as an elastic DOF.
    hinge : int, str or array (ns,)
        Hinge coordinate ``delta = d . eta``: an index, a DOF name or the
        vector ``d``.
    nonlinearity : HingeNonlinearity
        Restoring moment ``M_h(delta)`` that replaces the linear spring.
    K_lin : float, optional
        Linear hinge stiffness contained in ``system.K`` (default
        ``d^T K d``) which is removed.
    """

    def __init__(self, system, hinge, nonlinearity: HingeNonlinearity, K_lin: float | None = None):
        self.system = system
        ns = system.ns
        if isinstance(hinge, str):
            hinge = system.names.index(hinge)
        if np.ndim(hinge) == 0:
            d = np.zeros(ns)
            d[int(hinge)] = 1.0
        else:
            d = np.asarray(hinge, float)
        self.d = d
        self.nl = nonlinearity
        self.K_lin = float(d @ system.K @ d) if K_lin is None else float(K_lin)
        self._cache = {}

    # ------------------------------------------------------------------
    @property
    def names(self):
        return list(self.system.names[:self.system.ns])

    def linear(self, V: float, rho: float, gust: bool = False, outputs: bool = False) -> dict:
        """Matrices of the piecewise-affine model at speed ``V``.

        Returns ``ss`` (linear state space with ``K_lin``), ``A_free``
        (hinge spring removed), ``bh`` and ``c`` (``delta = c . x``) and
        ``cd`` (``delta' = cd . x``).
        """
        key = (float(V), float(rho), bool(gust), bool(outputs))
        if key in self._cache:
            return self._cache[key]
        sysm = self.system
        ss = sysm.state_space(V, rho, gust=gust, outputs=outputs)
        ns = sysm.ns
        S = slice(0, ns)
        tau = sysm.bref / V
        q = 0.5 * rho * V * V
        Mbar = sysm.M - q * tau * tau * sysm.rfa.A[2][S, S]
        n = ss.A.shape[0]
        bh = np.zeros(n)
        bh[ns:2 * ns] = np.linalg.solve(Mbar, self.d)
        c = np.zeros(n)
        c[:ns] = self.d
        cd = np.zeros(n)
        cd[ns:2 * ns] = self.d
        A_free = ss.A + self.K_lin * np.outer(bh, c)
        out = dict(ss=ss, A_free=A_free, bh=bh, c=c, cd=cd, n=n, ns=ns)
        self._cache[key] = out
        return out

    def equivalent_linear(self, V: float, rho: float, K: float) -> np.ndarray:
        """State matrix with a linear hinge spring ``K`` (e.g. ``K_eq``)."""
        L = self.linear(V, rho)
        return L["A_free"] - K * np.outer(L["bh"], L["c"])

    def initial_state(self, V: float, rho: float, gust: bool = False, **disp) -> np.ndarray:
        """State vector with prescribed displacements (``name=value``),
        rates (``etad_name=value``) or ``delta=value``; lag states zero."""
        L = self.linear(V, rho, gust)
        x = np.zeros(L["n"])
        ns = L["ns"]
        for k, v in disp.items():
            if k in ("delta", "hinge"):
                x[:ns] += v * self.d / (self.d @ self.d)
            elif k.startswith("etad_"):
                x[ns + self.names.index(k[5:])] = v
            else:
                x[self.names.index(k)] = v
        return x

    # ------------------------------------------------------------------
    def simulate(self, V: float, rho: float, t_end: float, x0=None, dt: float | None = None,
                 method: str | None = None, u=None, gust: bool | None = None,
                 integrator: str = "RK45", rtol: float = 1e-8, atol: float = 1e-11,
                 max_step: float | None = None, eps: float | None = None,
                 outputs: bool = False) -> SimulationResult:
        """Time history from the initial state ``x0``.

        ``x0``: full state vector or a dict for :meth:`initial_state`
        (default: hinge rotation ``2 x`` the nonlinearity scale).  ``dt``:
        output sampling (default: 1/100 of the shortest in-vacuo period).
        ``method``: ``'exact'``, ``'events'`` or ``'smooth'`` (see the
        module docstring).  ``u``: optional input function ``u(t) ->
        (m,)`` of the linear state-space inputs, e.g. ``[w_g, w_g']``
        when ``gust=True``.  This requires ``'events'`` or ``'smooth'``.
        ``outputs=True`` stores the outputs of the linear state space.
        Note that outputs involving ``eta''`` use the linear hinge
        spring, so use ``etadd`` from the result only when the spring
        is linear.
        """
        nl = self.nl
        if method is None:
            method = "exact" if (nl.piecewise_linear and u is None) else ("events" if nl.piecewise_linear else "smooth")
        if gust is None:
            gust = u is not None and self.system.rfa.Ag is not None
        L = self.linear(V, rho, gust, outputs)
        n, ns = L["n"], L["ns"]
        if x0 is None:
            x0 = self.initial_state(V, rho, gust, delta=2 * nl.scale)
        elif isinstance(x0, dict):
            x0 = self.initial_state(V, rho, gust, **x0)
        x0 = np.asarray(x0, float)
        if x0.size != n:
            raise ValueError(f"x0 has {x0.size} entries, state has {n}")
        if dt is None:
            w = np.sqrt(np.abs(np.linalg.eigvals(np.linalg.solve(self.system.M, self.system.K + np.eye(ns) * 1e-9))))
            dt = 2 * np.pi / max(w.max(), 1.0) / 100
        t = np.arange(0.0, t_end + 0.5 * dt, dt)
        if method == "exact":
            if u is not None:
                raise ValueError("method='exact' supports free responses only; use 'events'")
            if not nl.piecewise_linear:
                raise ValueError("method='exact' needs a piecewise-linear spring")
            X, nsw = self._march_exact(L, t, x0)
        elif method == "events":
            X, nsw = self._march_events(L, t, x0, u, integrator, rtol, atol, max_step)
        elif method == "smooth":
            X, nsw = self._march_smooth(L, t, x0, u, integrator, rtol, atol, max_step, eps)
        else:
            raise ValueError(method)
        y = None
        labels = []
        if outputs:
            ss = L["ss"]
            U = np.zeros((t.size, ss.B.shape[1])) if u is None else np.array([u(tt) for tt in t])
            y = X @ ss.C.T + U @ ss.D.T
            labels = list(ss.outputs)
        return SimulationResult(t, X, X[:, :ns], X[:, ns:2 * ns], X @ L["c"], self.names, V, method,
                                nsw, y, labels)

    # ---- exact piecewise-affine marching -----------------------------
    def _segment_mats(self, L):
        segs = self.nl.segments()
        n = L["n"]
        mats = []
        for lo, hi, s, m in segs:
            Aa = np.zeros((n + 1, n + 1))
            Aa[:n, :n] = L["A_free"] - s * np.outer(L["bh"], L["c"])
            Aa[:n, n] = -m * L["bh"]
            mats.append(Aa)
        return segs, mats

    def _march_exact(self, L, t, x0):
        segs, mats = self._segment_mats(L)
        dt = t[1] - t[0]
        n = L["n"]
        c, cd = L["c"], L["cd"]
        Phi = [sla.expm(Aa * dt) for Aa in mats]
        nl = self.nl
        X = np.zeros((t.size, n))
        X[0] = x0
        xa = np.r_[x0, 1.0]
        seg = nl.segment_of(c @ x0, cd @ x0)
        nsw = 0
        for k in range(1, t.size):
            rem = dt
            full = True
            while True:
                if full:
                    xn = Phi[seg] @ xa
                else:
                    xn = sla.expm(mats[seg] * rem) @ xa
                lo, hi = segs[seg][0], segs[seg][1]
                dn = c @ xn[:n]
                if lo <= dn <= hi:
                    xa = xn
                    break
                # switching inside the step: locate delta(tau) = corner
                bnd = hi if dn > hi else lo
                Aa = mats[seg]
                d0, v0 = c @ xa[:n], cd @ xa[:n]
                v1 = cd @ xn[:n]
                tau = _hermite_root(d0 - bnd, v0, dn - bnd, v1, rem)
                for _ in range(8):
                    xt = sla.expm(Aa * tau) @ xa
                    f = c @ xt[:n] - bnd
                    fp = cd @ xt[:n]
                    if abs(f) < 1e-13 * max(nl.scale, 1e-12):
                        break
                    step = f / fp if fp != 0 else 0.0
                    tau_new = tau - step
                    if not (0 < tau_new < rem):
                        tau_new = brentq(lambda s: c @ (sla.expm(Aa * s) @ xa)[:n] - bnd, 0.0, rem,
                                         xtol=1e-14 * rem)
                        xt = sla.expm(Aa * tau_new) @ xa
                        tau = tau_new
                        break
                    tau = tau_new
                xa = xt
                xa[:n] += (bnd - c @ xa[:n]) * c / (c @ c)     # snap exactly onto the corner
                seg = nl.segment_of(bnd, cd @ xa[:n])
                nsw += 1
                rem -= tau
                full = False
                if rem <= 1e-14 * dt:
                    break
            X[k] = xa[:n]
        return X, nsw

    # ---- solve_ivp with events ---------------------------------------
    def _march_events(self, L, t, x0, u, integrator, rtol, atol, max_step):
        nl = self.nl
        if not nl.piecewise_linear:
            return self._march_smooth(L, t, x0, u, integrator, rtol, atol, max_step, eps=0.0)
        segs = nl.segments()
        n = L["n"]
        A_free, bh, c, cd = L["A_free"], L["bh"], L["c"], L["cd"]
        B = L["ss"].B
        As = [A_free - s * np.outer(bh, c) for (_, _, s, _) in segs]
        bs = [-m * bh for (_, _, _, m) in segs]
        X = np.zeros((t.size, n))
        X[0] = x0
        tcur, x = 0.0, np.asarray(x0, float).copy()
        seg = nl.segment_of(c @ x, cd @ x)
        nsw = 0
        T = t[-1]
        ms = np.inf if max_step is None else max_step
        while tcur < T - 1e-14:
            Aj, bj = As[seg], bs[seg]
            lo, hi = segs[seg][0], segs[seg][1]
            if u is None:
                f = lambda tt, xx, Aj=Aj, bj=bj: Aj @ xx + bj
            else:
                f = lambda tt, xx, Aj=Aj, bj=bj: Aj @ xx + bj + B @ u(tt)
            evs = []
            if np.isfinite(hi):
                e_hi = lambda tt, xx, hi=hi: c @ xx - hi
                e_hi.terminal, e_hi.direction = True, 1
                evs.append(e_hi)
            if np.isfinite(lo):
                e_lo = lambda tt, xx, lo=lo: c @ xx - lo
                e_lo.terminal, e_lo.direction = True, -1
                evs.append(e_lo)
            k0 = np.searchsorted(t, tcur, side="right")
            jac = {} if integrator in ("RK45", "RK23", "DOP853") else {"jac": lambda tt, xx, Aj=Aj: Aj}
            sol = solve_ivp(f, (tcur, T), x, method=integrator, t_eval=t[k0:], events=evs or None,
                            rtol=rtol, atol=atol, max_step=ms, **jac)
            m = sol.t.size
            X[k0:k0 + m] = sol.y.T
            if sol.status == 1:      # event
                te = [e for e in sol.t_events if e.size]
                ye = [e for e in sol.y_events if e.size]
                tcur = float(te[0][0])
                x = ye[0][0].copy()
                dlt = c @ x
                bnd = hi if abs(dlt - hi) < abs(dlt - lo) else lo
                x += (bnd - dlt) * c / (c @ c)
                seg = nl.segment_of(bnd, cd @ x)
                nsw += 1
            else:
                break
        return X, nsw

    # ---- smoothed characteristic --------------------------------------
    def _march_smooth(self, L, t, x0, u, integrator, rtol, atol, max_step, eps):
        A_free, bh, c = L["A_free"], L["bh"], L["c"]
        B = L["ss"].B
        nl = self.nl
        if eps == 0.0:
            mom = nl.moment
        else:
            mom = lambda x: nl.smooth_moment(x, eps)
        if u is None:
            f = lambda tt, xx: A_free @ xx - bh * float(mom(c @ xx))
        else:
            f = lambda tt, xx: A_free @ xx - bh * float(mom(c @ xx)) + B @ u(tt)
        ms = np.inf if max_step is None else max_step
        sol = solve_ivp(f, (t[0], t[-1]), x0, method=integrator, t_eval=t, rtol=rtol, atol=atol,
                        max_step=ms)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y.T, 0

    # ------------------------------------------------------------------
    def stability(self, V: float, rho: float, K, wmin: float = 1e-3):
        """Largest real part (and its frequency) of the oscillatory roots of
        the equivalent linear system for hinge stiffness(es) ``K``."""
        L = self.linear(V, rho)
        K = np.atleast_1d(np.asarray(K, float))
        sig = np.zeros(K.size)
        om = np.zeros(K.size)
        P = np.outer(L["bh"], L["c"])
        for i, k in enumerate(K):
            lam = np.linalg.eigvals(L["A_free"] - k * P)
            lam = lam[lam.imag > wmin]
            j = int(np.argmax(lam.real))
            sig[i], om[i] = lam[j].real, lam[j].imag
        return sig, om


# ======================================================================
# helpers
# ======================================================================
def _hermite_root(f0, v0, f1, v1, h):
    """Root in (0, h) of the cubic Hermite interpolant (fallback: secant)."""
    # p(s) = h00 f0 + h10 h v0 + h01 f1 + h11 h v1,  s in [0, 1]
    a = 2 * f0 + h * v0 - 2 * f1 + h * v1
    b = -3 * f0 - 2 * h * v0 + 3 * f1 - h * v1
    cc = h * v0
    d = f0
    r = np.roots([a, b, cc, d])
    r = r[(np.abs(r.imag) < 1e-9) & (r.real > 0) & (r.real < 1)].real
    if r.size:
        return float(r.min() * h)
    return float(h * f0 / (f0 - f1)) if f0 != f1 else 0.5 * h


# ======================================================================
# velocity sweeps
# ======================================================================
def bifurcation_sweep(model: NonlinearAeroelastic, velocities, rho: float, t_end: float,
                      x0=None, dof: str = "delta", continuation: bool = True, method=None,
                      dt=None, window: float | None = None, perturb: float | None = None,
                      verbose: bool = False, keep_histories: bool = False) -> dict:
    """LCO amplitude versus speed by time marching.

    For each ``V`` the model is integrated over ``t_end``, and the last
    ``window`` (default 40 %) is analysed with :func:`lco_from_history`.
    With ``continuation`` the final state of one speed is the initial
    state of the next.  If that state has decayed below ``perturb``
    (default the nonlinearity scale), it is re-seeded with ``x0``, so
    that LCOs which need a finite disturbance are still found.

    Returns a dict of arrays ``V``, ``amp`` (half peak-to-peak of
    ``dof``), ``amp_h1``, ``freq`` (Hz), ``mean``, ``kind``,
    ``amp_dofs`` (nV, ns: amplitudes of all structural DOFs), ``peaks``
    (list of arrays of local maxima, for bifurcation diagrams) and
    optionally ``histories``.
    """
    velocities = np.asarray(velocities, float)
    ns = model.system.ns
    out = {k: np.full(velocities.size, np.nan) for k in ("amp", "amp_h1", "freq", "mean", "growth")}
    out["V"] = velocities
    out["kind"] = []
    out["peaks"] = []
    out["amp_dofs"] = np.zeros((velocities.size, ns))
    out["nswitch"] = np.zeros(velocities.size, int)
    hist = []
    xprev = None
    seed = perturb if perturb is not None else model.nl.scale
    for i, V in enumerate(velocities):
        xi = x0
        if continuation and xprev is not None:
            xi = xprev
        r = model.simulate(V, rho, t_end, x0=xi, dt=dt, method=method)
        y = r.delta if dof in ("delta", "hinge") else r.dof(dof)
        info = lco_from_history(r.t, y, window)
        out["amp"][i] = info.amplitude
        out["amp_h1"][i] = info.amplitude_h1
        out["freq"][i] = info.frequency
        out["mean"][i] = info.mean
        out["growth"][i] = info.growth
        out["kind"].append(info.kind)
        out["peaks"].append(info.peaks)
        out["nswitch"][i] = r.nswitch
        sel = r.t >= r.t[-1] - (0.4 * r.t[-1] if window is None else window)
        out["amp_dofs"][i] = 0.5 * (r.eta[sel].max(0) - r.eta[sel].min(0))
        xprev = r.x[-1].copy()
        if np.max(np.abs(r.delta[sel])) < seed or info.kind == "growing" or not np.all(np.isfinite(r.x[-1])):
            xprev = None                              # decayed or diverged: re-seed the next run
        if keep_histories:
            hist.append(r)
        if verbose:
            print(f"  V = {V:8.3f}  {info.kind:9s}  amp = {info.amplitude:.4e}  f = {info.frequency:.3f} Hz"
                  f"  growth = {info.growth:.3f}  switches = {r.nswitch}")
    if keep_histories:
        out["histories"] = hist
    return out


def flutter_speed_vs_stiffness(model: NonlinearAeroelastic, rho: float, K_values, V_range,
                               nV: int = 120) -> dict:
    """Lowest flutter speed ``V_f(K)`` of the equivalent linear system for
    each hinge stiffness ``K`` (scan on ``nV`` speeds in ``V_range`` and
    ``brentq`` refinement).  NaN if stable over the range."""
    K_values = np.atleast_1d(np.asarray(K_values, float))
    Vg = np.linspace(V_range[0], V_range[1], nV)
    Vf = np.full(K_values.size, np.nan)
    wf = np.full(K_values.size, np.nan)
    for i, K in enumerate(K_values):
        s = np.array([model.stability(V, rho, K)[0][0] for V in Vg])
        j = np.nonzero((s[:-1] < 0) & (s[1:] >= 0))[0]
        if j.size:
            a, b = Vg[j[0]], Vg[j[0] + 1]
            V0 = brentq(lambda V: model.stability(V, rho, K)[0][0], a, b, xtol=1e-6 * b)
            Vf[i] = V0
            wf[i] = model.stability(V0, rho, K)[1][0]
    return {"K": K_values, "V_f": Vf, "omega_f": wf}




# ======================================================================
# describing function
# ======================================================================
def hinge_compliance(model: NonlinearAeroelastic, V: float, rho: float, omega) -> np.ndarray:
    """Dynamic hinge compliance ``G(i w)`` of the aeroelastic system with
    the hinge spring removed: hinge rotation per unit applied hinge
    moment (TE down).

    With a hinge stiffness ``K`` the equivalent linear system has the
    characteristic equation ``1 - K G(s) = 0``.  It is therefore
    neutrally stable at ``w`` when ``G(i w)`` is real and
    ``K = 1 / G(i w)``.
    """
    L = model.linear(V, rho)
    A = L["A_free"]
    H, Q = sla.hessenberg(A, calc_q=True)
    bq = Q.T @ L["bh"]
    cq = L["c"] @ Q
    omega = np.atleast_1d(np.asarray(omega, float))
    I = np.eye(A.shape[0])
    return np.array([-(cq @ np.linalg.solve(1j * w * I - H, bq)) for w in omega])


def _hinge_neutral_points(model, V, rho, wgrid, Kmax):
    """All ``(w, K)`` with ``K = 1/G(i w)`` real, ``0 < K <= Kmax``."""
    G = hinge_compliance(model, V, rho, wgrid)
    im = G.imag
    pts = []
    for j in np.nonzero(np.sign(im[:-1]) * np.sign(im[1:]) < 0)[0]:
        f = lambda w: hinge_compliance(model, V, rho, w)[0].imag
        try:
            w0 = brentq(f, wgrid[j], wgrid[j + 1], xtol=1e-12 * wgrid[j + 1])
        except ValueError:
            continue
        g0 = hinge_compliance(model, V, rho, w0)[0]
        if g0.real <= 0:
            continue
        K = 1.0 / g0.real
        if K <= Kmax:
            pts.append((w0, K))
    return pts


@dataclass
class DFResult:
    """Describing-function LCO branches.

    Entry ``i`` is one LCO at speed ``V[i]``, with:

    * ``A[i]``: hinge amplitude (rad);
    * ``K_eq[i]``: equivalent stiffness;
    * ``omega[i]``: frequency (rad/s);
    * ``stable[i]``: stability flag;
    * ``branch[i]``: branch id. Points of one branch are continuous in
      ``V`` and ``omega``;
    * ``others_stable[i]``: whether all other roots of the equivalent
      linear system are stable. If not, a second mode is unstable at
      this ``K_eq`` and the single-harmonic prediction is doubtful.

    ``V_f0`` and ``V_fK`` are the flutter speeds at zero and reference
    stiffness (the lowest unstable speed in the scanned range).
    """

    V: np.ndarray
    A: np.ndarray
    K_eq: np.ndarray
    omega: np.ndarray
    stable: np.ndarray
    branch: np.ndarray
    K_ref: float
    V_f0: float = np.nan
    V_fK: float = np.nan
    others_stable: np.ndarray | None = None
    dof_amplitude: np.ndarray | None = None     # (n, ns) amplitudes of all structural DOFs
    dof_phase: np.ndarray | None = None         # (n, ns) phase relative to the hinge rotation [rad]
    names: list = field(default_factory=list)

    def dof(self, name: str):
        """Amplitudes of structural DOF ``name`` for every entry."""
        return self.dof_amplitude[:, self.names.index(name)]

    @property
    def freq(self):
        return self.omega / (2 * np.pi)

    def branches(self):
        """List of index arrays (sorted by ``V``), one per branch."""
        out = []
        for b in np.unique(self.branch):
            idx = np.nonzero(self.branch == b)[0]
            out.append(idx[np.argsort(self.V[idx])])
        return out

    def _interp(self, V, stable_only):
        out = []
        for idx in self.branches():
            Vs, As, Ws, St = self.V[idx], self.A[idx], self.omega[idx], self.stable[idx]
            for i in range(Vs.size - 1):
                if Vs[i] <= V <= Vs[i + 1] and Vs[i + 1] > Vs[i]:
                    if stable_only and not (St[i] and St[i + 1]):
                        continue
                    f = (V - Vs[i]) / (Vs[i + 1] - Vs[i])
                    out.append((As[i] + f * (As[i + 1] - As[i]), Ws[i] + f * (Ws[i + 1] - Ws[i])))
                    break
        return sorted(out)

    def amplitude(self, V: float, stable_only: bool = True) -> np.ndarray:
        """LCO amplitude(s) at speed ``V`` (linear interpolation along each branch)."""
        return np.array([a for a, _ in self._interp(V, stable_only)])

    def frequency(self, V: float, stable_only: bool = True) -> np.ndarray:
        """LCO frequency(ies) in rad/s at ``V``, ordered like :meth:`amplitude`."""
        return np.array([w for _, w in self._interp(V, stable_only)])


def describing_function_lco(model: NonlinearAeroelastic, rho: float, velocities,
                            omega=None, K_max: float | None = None, A_max: float | None = None,
                            nA: int = 400, link_tol: float = 0.15) -> DFResult:
    """Describing-function (first-harmonic balance) LCO amplitude versus speed.

    For each speed ``V``, the neutral points of the equivalent linear
    system are found from the hinge compliance (:func:`hinge_compliance`):
    ``Im G(i w) = 0`` with ``K* = 1 / G(i w)`` in ``[0, K_max]``.  The
    compliance is scanned on ``omega`` (default 400 points from 0.2 x the
    lowest to 3 x the highest in-vacuo frequency) and the roots are
    refined with ``brentq``.  These points are exactly the
    flutter-speed-versus-stiffness curves ``V_f(K)`` of all modes.  Each
    one gives LCO amplitude(s) ``A`` with ``N1(A) = K*`` (odd
    nonlinearity, no bias).

    Stability: the crossing root moves with ``d lam/dK``, from first-order
    eigenvalue perturbation.  The LCO is stable when
    ``Re(d lam/dK) dN1/dA < 0`` and all other roots of the equivalent
    system are stable.
    """
    nl = model.nl
    sysm = model.system
    velocities = np.asarray(velocities, float)
    Kref = nl.K
    K_max = 3.0 * Kref if K_max is None else K_max
    if omega is None:
        w = np.sqrt(np.abs(np.linalg.eigvals(np.linalg.solve(sysm.M, sysm.K))))
        w = w[w > 1e-6]
        omega = np.geomspace(0.2 * w.min(), 3 * w.max(), 400)
    omega = np.asarray(omega, float)
    A_max = 500 * nl.scale if A_max is None else A_max
    Ag = np.geomspace(1e-4 * nl.scale, A_max, nA)
    N1g = nl.describing_function(Ag)[0]
    rows = []
    shapes = []
    ns = sysm.ns
    for V in velocities:
        L = model.linear(V, rho)
        P = np.outer(L["bh"], L["c"])
        for w0, Ks in _hinge_neutral_points(model, V, rho, omega, K_max):
            AK = L["A_free"] - Ks * P
            lam, X = np.linalg.eig(AK)
            i0 = int(np.argmin(np.abs(lam - 1j * w0)))
            tol = 1e-7 * max(w0, 1.0)
            others = np.delete(lam, i0)
            others = others[np.abs(others + 1j * w0) > 1e-6 * max(w0, 1.0)]
            others_stable = bool(np.all(others.real < tol))
            Y = np.linalg.inv(X)
            dlam = -(Y[i0] @ P @ X[:, i0])          # d lam / dK
            g = N1g - Ks
            for a in np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]:
                Aa = brentq(lambda A: nl.keq(A)[0] - Ks, Ag[a], Ag[a + 1], xtol=1e-13 * Ag[a + 1])
                dA = 1e-6 * Aa
                dNdA = (nl.keq(Aa + dA)[0] - nl.keq(Aa - dA)[0]) / (2 * dA)
                stable = dlam.real * dNdA < 0
                rows.append([V, Aa, Ks, w0, float(stable), float(others_stable)])
                phi = X[:ns, i0] / (model.d @ X[:ns, i0])
                shapes.append(phi * Aa)
    rows = np.array(rows, float) if rows else np.zeros((0, 6))
    # branch linking by continuity in (V, omega)
    branch = np.full(rows.shape[0], -1)
    last = {}            # branch id -> (V, omega, K_eq)
    nb = 0
    for iv, V in enumerate(velocities):
        idx = np.nonzero(rows[:, 0] == V)[0] if rows.size else []
        used = set()
        prevV = velocities[iv - 1] if iv > 0 else None
        for i in idx:
            best, bd = None, link_tol
            for b, (Vb, wb, ab) in last.items():
                if b in used or prevV is None or Vb != prevV:
                    continue
                dd = abs(rows[i, 3] - wb) / wb + abs(rows[i, 2] - ab) / max(Kref, 1e-300)
                if dd < bd:
                    best, bd = b, dd
            if best is None:
                best = nb
                nb += 1
            used.add(best)
            branch[i] = best
            last[best] = (V, rows[i, 3], rows[i, 2])
    vr = (velocities.min(), velocities.max())
    f = flutter_speed_vs_stiffness(model, rho, [0.0, Kref], vr, nV=max(60, velocities.size))
    return DFResult(rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3], rows[:, 4].astype(bool),
                    branch, Kref, f["V_f"][0], f["V_f"][1], rows[:, 5].astype(bool),
                    np.abs(np.array(shapes)).reshape(-1, ns), np.angle(np.array(shapes)).reshape(-1, ns),
                    model.names)
