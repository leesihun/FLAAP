"""Rational function approximation of GAFs (port of ``generalizedforce.f90``
/ ``generate_A_matrices`` and the Roger least-squares in ``calculation.f90``).

Roger (1977)::

    Q(p) ~ A0 + A1 p + A2 p^2 + sum_l A_{l+2} p / (p + beta_l),   p = i k

with lag roots ``beta_l = 1.7 k_max (l / (n_lag + 1))^2`` (the Fortran
default).  Least squares is done element-wise on the real/imag parts
of all tabulated ``k`` (the Fortran ``C``/``ANG`` arrays).  Options:

* ``constrain_k0``: force exact match at ``k = 0`` (steady GAF),
* gust column with ``A2 = 0`` (no gust acceleration term),
* Karpel minimum-state approximation (:class:`MinimumStateRFA`).

Time-domain realisation (``s`` Laplace variable, ``p = s b / U``)::

    q Q(p) eta = q [A0 + A1 (b/U) s + A2 (b/U)^2 s^2] eta + q sum_l x_l
    x_l' = A_{l+2} eta' - (U/b) beta_l x_l
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .gaf import GAFTable


def default_lags(kmax: float, nlag: int) -> np.ndarray:
    return np.array([1.7 * kmax * (l / (nlag + 1)) ** 2 for l in range(1, nlag + 1)])


@dataclass
class RogerRFA:
    A: np.ndarray            # (3+nlag, nm, nm)
    beta: np.ndarray         # (nlag,)
    Ag: np.ndarray | None    # (3+nlag, nm)  gust columns (Ag[2] = 0)
    bref: float
    mach: float
    names: list[str]
    fit_error: float = 0.0
    kmax: float = np.inf        # highest tabulated reduced frequency of the fit

    @property
    def nlag(self) -> int:
        return self.beta.size

    def evaluate(self, k) -> tuple[np.ndarray, np.ndarray | None]:
        p = 1j * np.asarray(k, float)
        p = p if p.ndim else p[None]
        terms = [np.ones_like(p), p, p * p] + [p / (p + b) for b in self.beta]
        Q = np.einsum("tk,tij->kij", np.array(terms), self.A)
        Qg = None if self.Ag is None else np.einsum("tk,ti->ki", np.array(terms), self.Ag)
        return Q, Qg


def _basis(k, beta, with_acc=True):
    p = 1j * k
    cols = [np.ones_like(p), p] + ([p * p] if with_acc else []) + [p / (p + b) for b in beta]
    return np.array(cols).T   # (nk, nterm)


def _lsq(k, data, beta, with_acc, constrain_k0, weights=None, constrain_slope=False):
    """Solve for real coefficient vectors so that sum_t c_t B_t(k) ~ data(k).

    ``constrain_k0``: ``A0 = Re Q(0)`` exactly.  ``constrain_slope``
    (requires ``constrain_k0``): also match the low-frequency slope
    ``dQ/dp(0) = A1 + sum_l A_l / beta_l`` to ``Im Q(k1) / k1`` at the
    smallest positive tabulated ``k1`` - this pins the quasi-steady
    damping / rate derivatives (rigid-body dynamics, short period).
    Returns coefficients in the order [A0, A1, (A2), A_lag...].
    """
    if weights is not None and np.ndim(weights) == 2:
        # per-column weights: independent least-squares problems
        cols = [_lsq(k, data[:, [j]], beta, with_acc, constrain_k0, weights[:, j], constrain_slope)
                for j in range(data.shape[1])]
        return np.hstack(cols)
    B = _basis(k, beta, with_acc)
    w = np.ones(k.size) if weights is None else np.asarray(weights, float)
    W = w[:, None]
    if not constrain_k0:
        Ar = np.vstack([(B.real * W), (B.imag * W)])
        dr = np.vstack([(data.real * W), (data.imag * W)])
        sol, *_ = np.linalg.lstsq(Ar, dr, rcond=None)
        return sol
    i0 = int(np.argmin(np.abs(k)))
    a0 = data[i0].real
    res = data - a0[None, :]
    if constrain_slope:
        pos = np.where(k > 1e-12)[0]
        i1 = pos[np.argmin(k[pos])]
        slope = res[i1].imag / k[i1]
        p = 1j * k
        res = res - p[:, None] * slope[None, :]
        # remaining basis: (A2 p^2), lag_l: p/(p+b) - p/b
        cols = ([p * p] if with_acc else []) + [p / (p + b) - p / b for b in beta]
        Bc = np.array(cols).T
        Ar = np.vstack([Bc.real * W, Bc.imag * W])
        dr = np.vstack([res.real * W, res.imag * W])
        sol, *_ = np.linalg.lstsq(Ar, dr, rcond=None)
        lag = sol[(1 if with_acc else 0):]
        a1 = slope - np.sum(lag / beta[:, None], axis=0)
        return np.vstack([a0[None, :], a1[None, :], sol])
    Ar = np.vstack([B.real[:, 1:] * W, B.imag[:, 1:] * W])
    dr = np.vstack([res.real * W, res.imag * W])
    sol, *_ = np.linalg.lstsq(Ar, dr, rcond=None)
    return np.vstack([a0[None, :], sol])


def auto_weights(k, data, k_w: float = 0.1, floor: float = 0.1) -> np.ndarray:
    """Per-element weights ``1 / ((1 + k/k_w) (|Q_e(k)| + floor max|Q_e|))``.

    Makes the fit approximately *relative* for every element (small
    cross-coupling terms are not swamped by large ones) and emphasises the
    low reduced frequencies that govern rigid-body / short-period dynamics
    and gust response.  Flutter-relevant frequencies are still fitted well
    (the Goland flutter speed changes by < 1 %).
    """
    k = np.asarray(k, float)
    d = np.abs(np.asarray(data).reshape(k.size, -1))
    scale = d.max(axis=0)
    scale[scale == 0] = 1.0
    return 1.0 / ((1.0 + k / k_w)[:, None] * (d + floor * scale[None, :]))


def _weights(weights, k, data):
    if isinstance(weights, str):
        if weights != "auto":
            raise ValueError(f"unknown weighting {weights!r}")
        return auto_weights(k, data)
    return weights


def fit_roger(table: GAFTable, nlag: int = 4, beta: np.ndarray | None = None,
              constrain_k0: bool = True, weights="auto", constrain_slope: bool = True) -> RogerRFA:
    """Roger RFA of a GAF table.

    By default ``Q(0)`` and the low-frequency slope ``dQ/dp(0)`` are
    matched exactly (see :func:`_lsq`); this is essential for rigid-body
    modes (flight dynamics, manoeuvres, gust response with free aircraft).
    ``weights``: ``"auto"`` (default, :func:`auto_weights`), ``None``
    (uniform), an array (nk,) per frequency or (nk, nm*nm) per element.
    """
    k = table.k
    kmax = float(k.max())
    beta = default_lags(kmax, nlag) if beta is None else np.asarray(beta, float)
    nk, nm, _ = table.Q.shape
    data = table.Q.reshape(nk, nm * nm)
    c0 = constrain_k0 and np.min(np.abs(k)) < 1e-9
    cs = constrain_slope and c0 and np.sum(k > 1e-12) > 0
    sol = _lsq(k, data, beta, True, c0, _weights(weights, k, data), cs)
    A = sol.reshape(3 + beta.size, nm, nm)
    Ag = None
    if table.Qg is not None:
        solg = _lsq(k, table.Qg, beta, False, c0, _weights(weights, k, table.Qg), cs)
        Ag = np.zeros((3 + beta.size, nm))
        Ag[:2] = solg[:2]
        Ag[3:] = solg[2:]
    rfa = RogerRFA(A, beta, Ag, table.bref, table.mach, list(table.names), kmax=kmax)
    Qf, _ = rfa.evaluate(k)
    rfa.fit_error = float(np.linalg.norm(Qf - table.Q) / max(np.linalg.norm(table.Q), 1e-30))
    return rfa


def fit_rows(k, data, beta, constrain_k0: bool = True, with_acc: bool = True,
             weights="auto", constrain_slope: bool = True) -> np.ndarray:
    """Roger fit of arbitrary tabulated transfer rows with given lag roots.

    ``data`` has shape (nk, ...); returns coefficients of shape
    (3 + nlag, ...) with the acceleration term zero if ``with_acc`` is False.
    Used for gust columns and for aerodynamic load / pressure outputs
    that must share the lag states of the GAF RFA.
    """
    k = np.asarray(k, float)
    beta = np.asarray(beta, float)
    data = np.asarray(data)
    shp = data.shape[1:]
    flat = data.reshape(k.size, -1)
    c0 = constrain_k0 and np.min(np.abs(k)) < 1e-9
    cs = constrain_slope and c0 and np.sum(k > 1e-12) > 0
    sol = _lsq(k, flat, beta, with_acc, c0, _weights(weights, k, flat), cs)
    out = np.zeros((3 + beta.size, flat.shape[1]))
    if with_acc:
        out[:] = sol
    else:
        out[:2] = sol[:2]
        out[3:] = sol[2:]
    return out.reshape((3 + beta.size,) + shp)


def evaluate_rows(coef, beta, k) -> np.ndarray:
    """Evaluate Roger rows ``coef`` (3+nlag, ...) at reduced frequencies k."""
    p = 1j * np.atleast_1d(np.asarray(k, float))
    terms = np.array([np.ones_like(p), p, p * p] + [p / (p + b) for b in beta])
    return np.tensordot(terms.T, coef, axes=(1, 0))


@dataclass
class MinimumStateRFA:
    """Karpel minimum-state RFA  Q(p) = A0 + A1 p + A2 p^2 + D (pI - R)^-1 E p."""

    A: np.ndarray       # (3, nm, nm)
    D: np.ndarray       # (nm, nlag)
    E: np.ndarray       # (nlag, nm)
    R: np.ndarray       # (nlag,) negative real roots (diag)
    bref: float
    mach: float
    names: list[str]
    fit_error: float = 0.0

    def evaluate(self, k):
        p = 1j * np.atleast_1d(np.asarray(k, float))
        out = []
        for pp in p:
            lag = self.D @ np.diag(1.0 / (pp - self.R)) @ self.E * pp
            out.append(self.A[0] + self.A[1] * pp + self.A[2] * pp * pp + lag)
        return np.array(out), None


def fit_minimum_state(table: GAFTable, nlag: int = 4, iterations: int = 30,
                      R: np.ndarray | None = None) -> MinimumStateRFA:
    """Karpel (1982) minimum-state fit by alternating least squares on D, E."""
    k = table.k
    nk, nm, _ = table.Q.shape
    kmax = float(k.max())
    R = -default_lags(kmax, nlag) if R is None else np.asarray(R, float)
    p = 1j * k
    # initial E: Roger-like (identity-ish blocks)
    rng = np.random.default_rng(0)
    E = rng.standard_normal((nlag, nm)) * 0.1 + np.eye(nlag, nm)
    i0 = int(np.argmin(np.abs(k)))
    A0 = table.Q[i0].real

    def solve_D(E):
        # unknowns: A1, A2, D (A0 fixed) per row of Q
        cols = []
        for n in range(nk):
            lag = (np.diag(p[n] / (p[n] - R)) @ E)            # (nlag, nm)
            # Q_row = A1 p + A2 p^2 + D lag  -> for each output row i: x = [A1_i, A2_i, D_i]
            Phi = np.vstack([p[n] * np.eye(nm), p[n] ** 2 * np.eye(nm), lag])   # (2nm+nlag, nm)
            cols.append(Phi)
        Phi = np.array(cols)                                     # (nk, 2nm+nlag, nm)
        Mr = np.concatenate([Phi.real.transpose(0, 2, 1).reshape(-1, Phi.shape[1]),
                             Phi.imag.transpose(0, 2, 1).reshape(-1, Phi.shape[1])])
        Y = table.Q - A0[None]
        Yr = np.concatenate([Y.real.transpose(0, 2, 1).reshape(-1, nm),
                             Y.imag.transpose(0, 2, 1).reshape(-1, nm)])
        X, *_ = np.linalg.lstsq(Mr, Yr, rcond=None)              # (2nm+nlag, nm) -> rows i
        A1 = X[:nm].T
        A2 = X[nm:2 * nm].T
        D = X[2 * nm:].T
        return A1, A2, D

    def solve_E(A1, A2, D):
        # Q - A0 - A1 p - A2 p^2 = D diag(p/(p-R)) E  -> linear in E (column-wise)
        blocks, rhs = [], []
        for n in range(nk):
            G = D @ np.diag(p[n] / (p[n] - R))                  # (nm, nlag)
            Y = table.Q[n] - A0 - A1 * p[n] - A2 * p[n] ** 2
            blocks.append(G)
            rhs.append(Y)
        G = np.concatenate(blocks)
        Y = np.concatenate(rhs)
        Gr = np.concatenate([G.real, G.imag])
        Yr = np.concatenate([Y.real, Y.imag])
        E, *_ = np.linalg.lstsq(Gr, Yr, rcond=None)
        return E

    for _ in range(iterations):
        A1, A2, D = solve_D(E)
        E = solve_E(A1, A2, D)
    A1, A2, D = solve_D(E)
    ms = MinimumStateRFA(np.array([A0, A1, A2]), D, E, R, table.bref, table.mach, list(table.names))
    Qf, _ = ms.evaluate(k)
    ms.fit_error = float(np.linalg.norm(Qf - table.Q) / max(np.linalg.norm(table.Q), 1e-30))
    return ms
