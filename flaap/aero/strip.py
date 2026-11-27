"""Theodorsen strip-theory aerodynamics for beam (stick) models.

Used as an independent check of the flutter solvers (the classical Goland
wing result 137 m/s / 70.7 rad/s is a strip-theory result) and as a cheap
aerodynamic option.  For mode shapes given by the bending deflection
``w(y)`` (up) and twist ``th(y)`` (nose-up) about an elastic axis located
at ``a`` semichords aft of mid-chord::

    L/q = 2 pi b [k^2 w/b + i k th + a k^2 th]
        + 4 pi b C(k) [-i k w/b + th + (1/2 - a) i k th]
    M/q = 2 pi b^2 [a k^2 w/b - (1/2 - a) i k th + (1/8 + a^2) k^2 th]
        + 4 pi b^2 (a + 1/2) C(k) [-i k w/b + th + (1/2 - a) i k th]

and ``Q_ij = Int (L_j w_i + M_j th_i) dy`` (``q`` = dynamic pressure).
Compressibility can be approximated with a Prandtl-Glauert factor on the
circulatory lift slope.  Gusts use the Sears function.
"""
from __future__ import annotations

import numpy as np
from scipy.special import hankel2

from .gaf import GAFTable


def theodorsen(k):
    """Theodorsen function C(k) (C(0) = 1)."""
    k = np.asarray(k, float)
    ks = np.where(k > 0, k, 1.0)
    H1, H0 = hankel2(1, ks), hankel2(0, ks)
    C = H1 / (H1 + 1j * H0)
    return np.where(k > 0, C, 1.0 + 0j)


def sears(k):
    """Sears function referenced to mid-chord (gust front at mid-chord)."""
    k = np.asarray(k, float)
    ks = np.where(k > 0, k, 1.0)
    from scipy.special import jv
    C = theodorsen(ks)
    S = C * (jv(0, ks) - 1j * jv(1, ks)) + 1j * jv(1, ks)
    return np.where(k > 0, S, 1.0 + 0j)


def strip_gaf_table(y_nodes, w_modes, th_modes, chord, x_ea_frac, kvals, mach=0.0,
                    names=None, nstrip=200, bref=None, y_tip=None, lift_slope_factor=None,
                    x_le=0.0) -> GAFTable:
    """Strip-theory GAF table.

    ``w_modes``, ``th_modes``: (Nn, nm) bending deflection at the elastic
    axis and nose-up twist at the beam nodes ``y_nodes``.  ``chord`` and
    ``x_ea_frac`` (elastic axis, fraction of chord from LE) may be scalars
    or callables of ``y``.  Gust column (Sears, gust reference at LE x_le).
    """
    y_nodes = np.asarray(y_nodes, float)
    y_tip = y_nodes[-1] if y_tip is None else y_tip
    edges = np.linspace(y_nodes[0], y_tip, nstrip + 1)
    ys = 0.5 * (edges[1:] + edges[:-1])
    dy = np.diff(edges)
    ch = np.array([chord(y) if callable(chord) else chord for y in ys])
    xf = np.array([x_ea_frac(y) if callable(x_ea_frac) else x_ea_frac for y in ys])
    b = 0.5 * ch
    a = 2 * xf - 1.0
    bref = float(np.mean(b)) if bref is None else bref
    W = np.array([np.interp(ys, y_nodes, w_modes[:, j]) for j in range(w_modes.shape[1])]).T
    T = np.array([np.interp(ys, y_nodes, th_modes[:, j]) for j in range(th_modes.shape[1])]).T
    nm = W.shape[1]
    fac = (1.0 / np.sqrt(1 - mach ** 2)) if lift_slope_factor is None else lift_slope_factor
    kvals = np.atleast_1d(np.asarray(kvals, float))
    Q = np.zeros((kvals.size, nm, nm), complex)
    Qg = np.zeros((kvals.size, nm), complex)
    for n, kr in enumerate(kvals):
        kl = kr * b / bref                 # local reduced frequency
        C = theodorsen(kl) * fac
        S = sears(kl) * fac
        wj, tj = W, T
        Lnc = 2 * np.pi * b[:, None] * (kl[:, None] ** 2 * wj / b[:, None] + 1j * kl[:, None] * tj
                                        + a[:, None] * kl[:, None] ** 2 * tj)
        circ = (-1j * kl[:, None] * wj / b[:, None] + tj + (0.5 - a[:, None]) * 1j * kl[:, None] * tj)
        L = Lnc + 4 * np.pi * b[:, None] * C[:, None] * circ
        Mnc = 2 * np.pi * b[:, None] ** 2 * (a[:, None] * kl[:, None] ** 2 * wj / b[:, None]
                                             - (0.5 - a[:, None]) * 1j * kl[:, None] * tj
                                             + (0.125 + a[:, None] ** 2) * kl[:, None] ** 2 * tj)
        M = Mnc + 4 * np.pi * b[:, None] ** 2 * (a[:, None] + 0.5) * C[:, None] * circ
        Q[n] = (W * dy[:, None]).T @ L + (T * dy[:, None]).T @ M
        # gust: unit wg/U, phase referenced to the gust front at x_le (mid-chord at x_le+b)
        ph = np.exp(-1j * (kr / bref) * b)   # mid-chord lies b behind the gust reference (LE)
        Lg = 4 * np.pi * b * S * ph
        Mg = 4 * np.pi * b ** 2 * (a + 0.5) * S * ph
        Qg[n] = (W * dy[:, None]).T @ Lg + (T * dy[:, None]).T @ Mg
    names = names or [f"mode{i + 1}" for i in range(nm)]
    return GAFTable(kvals, Q, Qg, mach, bref, list(names))
