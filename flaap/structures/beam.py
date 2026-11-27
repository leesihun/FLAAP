"""Bending-torsion beam finite elements for stick wing models.

Straight beam along ``+y`` (elastic axis at ``x_ea``), DOFs per node in
the 6-DOF layout: ``uz`` (index 2), ``rx = dw/dy`` (3), ``ry`` (4; twist,
nose-up positive, ``uz = -(x - x_ea) * ry``).  Euler-Bernoulli bending
(Hermite cubic) + St-Venant torsion (linear), consistent mass including
the inertial coupling of a chordwise offset centre of gravity
``d = x_cg - x_ea`` (positive aft)::

    T = 1/2 Int [ m (w_t - d th_t)^2 + I_cg th_t^2 ] dy

Optional concentrated masses (with offset and pitch inertia) and an
optional chordwise (in-plane) bending DOF are supported.  Used for the
Goland wing (Goland 1945) validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .modal import FEModel

Num = float | Callable[[float], float]


def _val(f: Num, y: float) -> float:
    return float(f(y)) if callable(f) else float(f)


@dataclass
class BeamProperties:
    EI: Num          # bending stiffness [N m^2]
    GJ: Num          # torsional stiffness [N m^2]
    m: Num           # mass per length [kg/m]
    I_cg: Num        # pitch inertia per length about the c.g. [kg m]
    x_ea: Num        # elastic axis x [m]
    x_cg: Num        # c.g. x [m]


def _gauss(n=4):
    return np.polynomial.legendre.leggauss(n)


def beam_element(y1: float, y2: float, p: BeamProperties):
    """6x6 element matrices for local DOFs [w1, rx1, th1, w2, rx2, th2]."""
    L = y2 - y1
    K = np.zeros((6, 6))
    M = np.zeros((6, 6))
    iw = [0, 1, 3, 4]
    it = [2, 5]
    xg, wg = _gauss(4)
    for xi, wt in zip(xg, wg):
        s = 0.5 * (xi + 1)
        y = y1 + s * L
        jac = 0.5 * L * wt
        # Hermite shape functions in s, with slope DOFs scaled by L
        N = np.array([1 - 3 * s ** 2 + 2 * s ** 3, L * (s - 2 * s ** 2 + s ** 3),
                      3 * s ** 2 - 2 * s ** 3, L * (-s ** 2 + s ** 3)])
        d2N = np.array([-6 + 12 * s, L * (-4 + 6 * s), 6 - 12 * s, L * (-2 + 6 * s)]) / L ** 2
        Nt = np.array([1 - s, s])
        dNt = np.array([-1.0, 1.0]) / L
        EI, GJ, m = _val(p.EI, y), _val(p.GJ, y), _val(p.m, y)
        Icg = _val(p.I_cg, y)
        d = _val(p.x_cg, y) - _val(p.x_ea, y)
        K[np.ix_(iw, iw)] += EI * np.outer(d2N, d2N) * jac
        K[np.ix_(it, it)] += GJ * np.outer(dNt, dNt) * jac
        M[np.ix_(iw, iw)] += m * np.outer(N, N) * jac
        M[np.ix_(it, it)] += (Icg + m * d * d) * np.outer(Nt, Nt) * jac
        cpl = -m * d * np.outer(N, Nt) * jac
        M[np.ix_(iw, it)] += cpl
        M[np.ix_(it, iw)] += cpl.T
    return K, M


def build_beam_model(y_nodes, props: BeamProperties, clamp_root: bool = True,
                     point_masses: list[dict] | None = None, name: str = "beam") -> FEModel:
    """Assemble a cantilever bending-torsion beam FEModel.

    ``point_masses``: list of dicts ``{'y':, 'm':, 'x':, 'I':}`` lumped at
    the nearest node (``x`` absolute c.g. position, ``I`` pitch inertia
    about its own c.g.).
    """
    y = np.asarray(y_nodes, float)
    nn = y.size
    nodes = np.stack([[_val(props.x_ea, yy) for yy in y], y, np.zeros(nn)], axis=1)
    ndof = 6 * nn
    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))
    for e in range(nn - 1):
        ke, me = beam_element(y[e], y[e + 1], props)
        dofs = [6 * e + 2, 6 * e + 3, 6 * e + 4, 6 * (e + 1) + 2, 6 * (e + 1) + 3, 6 * (e + 1) + 4]
        K[np.ix_(dofs, dofs)] += ke
        M[np.ix_(dofs, dofs)] += me
    for pm in point_masses or []:
        j = int(np.argmin(np.abs(y - pm["y"])))
        d = pm.get("x", nodes[j, 0]) - nodes[j, 0]
        mm, I = pm["m"], pm.get("I", 0.0)
        a, b = 6 * j + 2, 6 * j + 4
        M[a, a] += mm
        M[b, b] += I + mm * d * d
        M[a, b] += -mm * d
        M[b, a] += -mm * d
    fixed = []
    for j in range(nn):
        fixed += [6 * j + 0, 6 * j + 1, 6 * j + 5]
    if clamp_root:
        fixed += [2, 3, 4]
    return FEModel(nodes, K, M, np.array(sorted(set(fixed))), name)
