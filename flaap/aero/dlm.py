"""Subsonic Doublet-Lattice Method (Albano & Rodden 1969).

Replaces the Fortran ``DLM.f90`` / ``AIC_calculation`` / ``A_matrix``
routines (only their object files survived).  Formulation:

* steady part ``D0`` from horseshoe vortices (:mod:`flaap.aero.vlm`)
* oscillatory increment ``D1`` from the Landahl (1967) kernel, with the
  Laschka 11-term exponential approximation of the integrals ``I1``,
  ``I2`` and the quartic spanwise approximation of Rodden, Taylor &
  McIntosh (1998, J. Aircraft 35(5))::

      D = D0 + (dx_s / 8 pi) * int_{-e}^{e} [ P1/r1^2 + P2/r1^4 ] d(eta)

      P1 = (K1 exp(-i w x0/U) - K10) T1
      P2 = (K2 exp(-i w x0/U) - K20) T2*

Normalwash convention::

    w_r / U = sum_s D[r, s] * dCp_s ,  dCp = Cp_lower - Cp_upper

with simple-harmonic motion ``exp(+i w t)``.  The spanwise integration is
done in closed form; the weights that map the five kernel samples to the
integral depend only on geometry and are cached, so evaluating ``D(k)``
for many reduced frequencies is cheap.
"""
from __future__ import annotations

import numpy as np

from ..geometry import AeroMesh
from .vlm import steady_aic

# Laschka coefficients (Rodden et al. 1998, Table / Desmarais refinement not used)
_LA = np.array([0.24186198, -2.7918027, 24.991079, -111.59196, 271.43549,
                -305.75288, -41.18363, 545.98537, -644.78155, 328.72755,
                -64.279511])
_LC = 0.372
_NC = _LC * np.arange(1, 12)

_SNODE = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
_VINV = np.linalg.inv(np.vander(_SNODE, 5, increasing=True))  # b = VINV @ P
_GL_X, _GL_W = np.polynomial.legendre.leggauss(12)
_LAG_GL = np.vander(_GL_X, 5, increasing=True) @ _VINV          # (12,5) basis at GL pts

# relation between the Landahl-kernel integral and the horseshoe D0
# (verified numerically in tests/test_dlm.py::test_kernel_sign)
KERNEL_SIGN = 1.0

_PLANAR_TOL = 1e-6   # |zbar| / e below which the sending/receiving pair is coplanar
_FAR = 3.0           # receivers farther than _FAR * e use Gauss quadrature weights


# ----------------------------------------------------------------- kernel
def _laschka(u: np.ndarray, k1: np.ndarray):
    """I1 and 3*I2 for u >= 0 (arrays broadcastable); real arithmetic inside."""
    E = _LA * np.exp(-np.multiply.outer(u, _NC))            # (...,11)
    kk = k1[..., None]
    den = 1.0 / (_NC * _NC + kk * kk)
    ir = _NC * den                                           # Re 1/s_n
    ii = -kk * den                                           # Im 1/s_n
    S1r = np.sum(E * ir, axis=-1)
    S1i = np.sum(E * ii, axis=-1)
    uu = u[..., None]
    NE = _NC * E
    S2r = np.sum(NE * (uu * ir + ir * ir - ii * ii), axis=-1)
    S2i = np.sum(NE * (uu * ii + 2 * ir * ii), axis=-1)
    sq = np.sqrt(1.0 + u * u)
    g = 1.0 - u / sq
    eu = np.exp(-1j * k1 * u)
    I1 = eu * ((g + k1 * S1i) - 1j * (k1 * S1r))
    Tr = 2 * S1r - S2r
    Ti = 2 * S1i - S2i
    I2_3 = eu * ((2 * g - u / sq ** 3 + k1 * Ti) - 1j * (k1 * Tr))
    return I1, I2_3


def landahl_kernel(x0, r1, k1, mach):
    """Landahl's K1, K2 and their steady values K10, K20.

    ``x0`` streamwise distance receiver - sender, ``r1`` lateral distance
    (> 0), ``k1 = w r1 / U``.  Arrays broadcast elementwise.
    """
    x0 = np.asarray(x0, float)
    r1 = np.asarray(r1, float)
    k1 = np.asarray(k1, float) * np.ones_like(r1)
    b2 = 1.0 - mach ** 2
    R = np.sqrt(x0 * x0 + b2 * r1 * r1)
    u1 = (mach * R - x0) / (b2 * r1)
    au = np.abs(u1)
    I1, I23 = _laschka(au, k1)
    neg = u1 < 0
    if np.any(neg):
        z = np.zeros_like(au)
        I10, I230 = _laschka(z, k1)
        I1 = np.where(neg, 2 * I10.real - I1.real + 1j * I1.imag, I1)
        I23 = np.where(neg, 2 * I230.real - I23.real + 1j * I23.imag, I23)
    sq = np.sqrt(1 + u1 * u1)
    eu = np.exp(-1j * k1 * u1)
    K1 = I1 + mach * r1 * eu / (R * sq)
    K2 = (-I23 - 1j * k1 * mach ** 2 * r1 ** 2 * eu / (R * R * sq)
          - (mach * r1 / R) * ((1 + u1 * u1) * b2 * r1 * r1 / (R * R) + 2
                               + mach * r1 * u1 / R) * eu / sq ** 3)
    K10 = 1.0 + x0 / R
    K20 = -2.0 - (x0 / R) * (2.0 + b2 * r1 * r1 / (R * R))
    return K1, K2, K10, K20


# ------------------------------------------------------------ integration
def _moments(ybar, zbar, e, planar):
    """Spanwise moments for the quartic integration (near field).

    Returns weights W1 (…,5) for int P1/r1^2 and W2 (…,5) for int P2/r1^4.
    """
    ta = -e - ybar
    tb = e - ybar
    z = np.abs(zbar)
    zs = np.where(planar, 1.0, z)
    z2 = np.where(planar, 0.0, z * z)
    with np.errstate(divide="ignore", invalid="ignore"):
        J0 = np.where(planar, 1.0 / ta - 1.0 / tb,
                      (np.arctan(tb / zs) - np.arctan(ta / zs)) / zs)
        J1 = np.where(planar, np.log(np.abs(tb / ta)),
                      0.5 * np.log((tb * tb + z2) / (ta * ta + z2)))
    J2 = (tb - ta) - z2 * J0
    J3 = 0.5 * (tb ** 2 - ta ** 2) - z2 * J1
    J4 = (tb ** 3 - ta ** 3) / 3.0 - z2 * J2
    J = [J0, J1, J2, J3, J4]
    # second family (only non-planar)
    with np.errstate(divide="ignore", invalid="ignore"):
        F = lambda t: t / (t * t + z2)                       # noqa: E731
        L0 = np.where(planar, 0.0, (F(tb) - F(ta) + J0) / (2 * np.where(planar, 1.0, z2)))
        L1 = np.where(planar, 0.0, -0.5 * (1.0 / (tb * tb + z2) - 1.0 / (ta * ta + z2)))
    L2 = J0 - z2 * L0
    L3 = J1 - z2 * L1
    L4 = J2 - z2 * L2
    L = [L0, L1, L2, L3, L4]
    # moments in eta: int eta^m f = sum_k C(m,k) ybar^(m-k) J_k
    from math import comb
    W1 = np.zeros(ybar.shape + (5,))
    W2 = np.zeros(ybar.shape + (5,))
    for m in range(5):
        Im = sum(comb(m, kk) * ybar ** (m - kk) * J[kk] for kk in range(m + 1))
        Lm = sum(comb(m, kk) * ybar ** (m - kk) * L[kk] for kk in range(m + 1))
        sc = e ** m
        W1 += (Im / sc)[..., None] * _VINV[m][None, :]
        W2 += (Lm / sc)[..., None] * _VINV[m][None, :]
    W2 = np.where(planar[..., None], 0.0, W2)
    return W1, W2


def _gauss_weights(ybar, zbar, e):
    """Far-field weights by 12-point Gauss quadrature of the Lagrange basis."""
    eta = e[..., None] * _GL_X                              # (…,12)
    r2 = (ybar[..., None] - eta) ** 2 + zbar[..., None] ** 2
    w1 = (e[..., None] * _GL_W) / r2                        # (…,12)
    w2 = w1 / r2
    return w1 @ _LAG_GL, w2 @ _LAG_GL


class _PairBlock:
    """Geometry of all receiver/sender pairs for one set of sending boxes."""

    def __init__(self, rec, nrm_gamma, snd: AeroMesh, factor: float):
        self.factor = factor
        mid = 0.5 * (snd.A + snd.B)
        d = rec[:, None, :] - mid[None, :, :]
        cg, sg = np.cos(snd.gamma)[None, :], np.sin(snd.gamma)[None, :]
        ybar = d[..., 1] * cg + d[..., 2] * sg
        zbar = -d[..., 1] * sg + d[..., 2] * cg
        e = np.broadcast_to(snd.e[None, :], ybar.shape)
        eta = e[..., None] * _SNODE                          # (Nr,Ns,5)
        self.x0 = d[..., 0][..., None] - eta * snd.tan_sweep[None, :, None]
        y0 = ybar[..., None] - eta
        self.r1 = np.sqrt(y0 * y0 + zbar[..., None] ** 2)
        dg = nrm_gamma[:, None] - snd.gamma[None, :]
        self.T1 = np.cos(dg)
        self.T2 = zbar[..., None] * (zbar[..., None] * np.cos(dg)[..., None]
                                     - y0 * np.sin(dg)[..., None])
        planar = np.abs(zbar) < _PLANAR_TOL * e
        far = np.hypot(ybar, zbar) > _FAR * e
        W1 = np.empty(ybar.shape + (5,))
        W2 = np.empty(ybar.shape + (5,))
        if np.any(far):
            W1[far], W2[far] = _gauss_weights(ybar[far], zbar[far], e[far])
        nf = ~far
        if np.any(nf):
            W1[nf], W2[nf] = _moments(ybar[nf], zbar[nf], e[nf], planar[nf])
        W2[planar] = 0.0
        self.planar = planar
        self.W1 = W1
        self.W2 = W2
        self.dx = snd.chord[None, :]

    def increment(self, w_over_u: float, mach: float) -> np.ndarray:
        if w_over_u == 0.0:
            return np.zeros(self.x0.shape[:2], complex)
        out = np.empty(self.x0.shape[:2], complex)
        nr = self.x0.shape[0]
        chunk = max(1, int(4e5 // max(1, self.x0.shape[1] * 5)))
        for i0 in range(0, nr, chunk):
            sl = slice(i0, min(nr, i0 + chunk))
            x0 = self.x0[sl]
            r1 = self.r1[sl]
            zero = r1 < 1e-12 * (np.abs(x0) + 1.0)
            r1s = np.where(zero, 1.0, r1)
            K1, K2, K10, K20 = landahl_kernel(x0, r1s, w_over_u * r1s, mach)
            ph = np.exp(-1j * w_over_u * x0)
            P1 = (K1 * ph - K10) * self.T1[sl][..., None]
            P1 = np.where(zero, np.where(x0 > 0, 2 * ph - 2, 0.0) * self.T1[sl][..., None], P1)
            P2 = (K2 * ph - K20) * self.T2[sl]
            P2 = np.where(zero, 0.0, P2)
            s = np.sum(self.W1[sl] * P1, axis=-1) + np.sum(self.W2[sl] * P2, axis=-1)
            out[sl] = s * self.dx / (8 * np.pi)
        return KERNEL_SIGN * self.factor * out


class DLM:
    """Doublet-lattice AIC generator for one mesh and Mach number.

    >>> dlm = DLM(mesh, mach=0.5)
    >>> D = dlm.aic(k=0.3, bref=0.5)    # complex (N, N)
    """

    def __init__(self, mesh: AeroMesh, mach: float,
                 receivers: np.ndarray | None = None,
                 receiver_gamma: np.ndarray | None = None):
        if not 0.0 <= mach < 1.0:
            raise ValueError("DLM requires 0 <= M < 1")
        self.mesh = mesh
        self.mach = mach
        rec = mesh.colloc if receivers is None else receivers
        gam = mesh.gamma if receiver_gamma is None else receiver_gamma
        nrm = np.stack([np.zeros_like(gam), -np.sin(gam), np.cos(gam)], axis=1)
        self.D0 = steady_aic(mesh, mach, rec, nrm)
        self.blocks = [_PairBlock(rec, gam, mesh, 1.0)]
        if mesh.symmetry in ("symmetric", "antisymmetric"):
            f = 1.0 if mesh.symmetry == "symmetric" else -1.0
            self.blocks.append(_PairBlock(rec, gam, mesh.image(), f))
        self._cache: dict[float, np.ndarray] = {}

    def aic(self, k: float, bref: float) -> np.ndarray:
        """Normalwash AIC ``D(k)`` with ``k = w * bref / U``."""
        key = round(float(k) / bref, 14)
        if key in self._cache:
            return self._cache[key]
        wu = float(k) / bref
        D = self.D0.astype(complex)
        for blk in self.blocks:
            D = D + blk.increment(wu, self.mach)
        self._cache[key] = D
        return D

    def clear_cache(self):
        self._cache.clear()
