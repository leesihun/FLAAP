"""Analytic / semi-analytic reference solutions for structural validation.

* :data:`LEISSA_CFFF` - frequency parameters ``lambda = w a^2 sqrt(rho h / D)``
  of a square cantilever (CFFF) Kirchhoff plate, nu = 0.3 (Leissa, NASA SP-160,
  table 4.89; first five modes 3.492, 8.525, 21.43, 27.33, 31.11).
* :func:`eb_cantilever_hz` - Euler-Bernoulli cantilever,
  ``f_n = (beta_n L)^2 / (2 pi L^2) sqrt(E I / (rho A))`` with
  ``beta_n L = 1.87510, 4.69409, 7.85476, 10.99554, ...``.
* :func:`timoshenko_cantilever_hz` - Timoshenko cantilever (shear + rotary
  inertia) from a converged 1-D finite-element solution (linear elements,
  one-point shear integration, 400 elements: relative discretisation error
  < 1e-5 for the first modes).  Shear coefficient of a rectangle
  ``kappa = 10 (1 + nu) / (12 + 11 nu)`` (Cowper).
* :func:`plate_lambda` - converts plate frequencies [rad/s] to Leissa's
  ``lambda``.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

__all__ = ["LEISSA_CFFF", "BETA_L_CANTILEVER", "eb_cantilever_hz", "timoshenko_cantilever_hz",
           "plate_lambda", "kappa_rectangle"]

LEISSA_CFFF = np.array([3.492, 8.525, 21.43, 27.33, 31.11])
BETA_L_CANTILEVER = np.array([1.8751040687, 4.6940911330, 7.8547574382, 10.9955407349,
                              14.1371683910])


def eb_cantilever_hz(L: float, E: float, I: float, rho: float, A: float, n: int = 3) -> np.ndarray:
    """Euler-Bernoulli cantilever natural frequencies [Hz] (first ``n``)."""
    bl = BETA_L_CANTILEVER[:n]
    return bl ** 2 / (2 * np.pi * L ** 2) * np.sqrt(E * I / (rho * A))


def kappa_rectangle(nu: float) -> float:
    """Cowper shear coefficient of a rectangular section."""
    return 10 * (1 + nu) / (12 + 11 * nu)


def timoshenko_cantilever_hz(L: float, E: float, nu: float, I: float, rho: float, A: float,
                             n: int = 3, kappa: float | None = None, nel: int = 400) -> np.ndarray:
    """Timoshenko cantilever bending frequencies [Hz] (first ``n``)."""
    kappa = kappa_rectangle(nu) if kappa is None else kappa
    G = E / (2 * (1 + nu))
    le = L / nel
    nd = 2 * (nel + 1)
    K = np.zeros((nd, nd))
    M = np.zeros((nd, nd))
    kb = E * I / le * np.array([[0, 0, 0, 0], [0, 1, 0, -1], [0, 0, 0, 0], [0, -1, 0, 1]], float)
    b = np.array([-1 / le, -0.5, 1 / le, -0.5])
    ks = kappa * G * A * le * np.outer(b, b)
    m2 = le / 6 * np.array([[2, 1], [1, 2]], float)
    me = np.zeros((4, 4))
    me[np.ix_([0, 2], [0, 2])] = rho * A * m2
    me[np.ix_([1, 3], [1, 3])] = rho * I * m2
    for e in range(nel):
        d = np.arange(2 * e, 2 * e + 4)
        K[np.ix_(d, d)] += kb + ks
        M[np.ix_(d, d)] += me
    free = np.arange(2, nd)
    w2 = eigh(K[np.ix_(free, free)], M[np.ix_(free, free)], eigvals_only=True,
              subset_by_index=[0, n - 1])
    return np.sqrt(w2) / (2 * np.pi)


def plate_lambda(omega, a: float, h: float, E: float, nu: float, rho: float) -> np.ndarray:
    """``lambda = w a^2 sqrt(rho h / D)``, ``D = E h^3 / (12 (1 - nu^2))``."""
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    return np.asarray(omega, float) * a ** 2 * np.sqrt(rho * h / D)
