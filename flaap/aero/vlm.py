"""Steady subsonic vortex-lattice method (horseshoe vortices).

The influence matrix ``D0`` relates the box pressure jumps (positive =
lift, ``dCp = Cp_lower - Cp_upper``) to the normalwash at the collocation
points, normalised by the free-stream speed::

    w_r / U = sum_s D0[r, s] * dCp_s

Compressibility is included through the Prandtl-Glauert stretch
``x -> x / beta`` of the geometry (the steady planar/non-planar subsonic
kernel equals the incompressible one evaluated at ``x0/beta``).
This is the steady part of the Doublet-Lattice method (``D0`` in the
Albano-Rodden decomposition ``D = D0 + D1(k)``).
"""
from __future__ import annotations

import numpy as np

from ..geometry import AeroMesh

_EPS = 1e-12


def _segment(P, P1, P2):
    """Velocity at P (…,3) induced by a unit-strength segment P1->P2."""
    r1 = P - P1
    r2 = P - P2
    r0 = P2 - P1
    cr = np.cross(r1, r2)
    cr2 = np.sum(cr * cr, axis=-1)
    n1 = np.linalg.norm(r1, axis=-1)
    n2 = np.linalg.norm(r2, axis=-1)
    scale0 = np.sum(r0 * r0, axis=-1)
    bad = (cr2 < _EPS * scale0 * scale0) | (n1 < 1e-14) | (n2 < 1e-14)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (np.sum(r0 * r1, axis=-1) / n1 - np.sum(r0 * r2, axis=-1) / n2) / (4 * np.pi * cr2)
    f = np.where(bad, 0.0, f)
    return cr * f[..., None]


def _semi_infinite(P, P0, d):
    """Unit vortex starting at P0 running to infinity along unit vector d."""
    r = P - P0
    dxr = np.cross(np.broadcast_to(d, r.shape), r)
    c2 = np.sum(dxr * dxr, axis=-1)
    nr = np.linalg.norm(r, axis=-1)
    scale = np.maximum(nr, 1e-30)
    bad = c2 < _EPS * scale * scale
    with np.errstate(divide="ignore", invalid="ignore"):
        f = (1.0 + np.sum(r * d, axis=-1) / nr) / (4 * np.pi * c2)
    f = np.where(bad, 0.0, f)
    return dxr * f[..., None]


def horseshoe_velocity(points: np.ndarray, A: np.ndarray, B: np.ndarray,
                       wake_dir=(1.0, 0.0, 0.0)) -> np.ndarray:
    """Velocity (Npts, Nhs, 3) from unit horseshoes A->B with legs to +inf."""
    P = points[:, None, :]
    A_ = A[None, :, :]
    B_ = B[None, :, :]
    d = np.asarray(wake_dir, float)
    v = _segment(P, A_, B_)
    v += _semi_infinite(P, B_, d)
    v -= _semi_infinite(P, A_, d)
    return v


def _stretch(p: np.ndarray, beta: float) -> np.ndarray:
    q = np.array(p, dtype=float, copy=True)
    q[..., 0] /= beta
    return q


def steady_aic(mesh: AeroMesh, mach: float = 0.0,
               receivers: np.ndarray | None = None,
               normals: np.ndarray | None = None) -> np.ndarray:
    """Steady normalwash AIC ``D0`` (Nr x Ns), symmetry included."""
    if mach >= 1.0:
        raise ValueError("VLM/DLM are subsonic methods (M < 1)")
    beta = np.sqrt(1.0 - mach ** 2)
    pts = mesh.colloc if receivers is None else receivers
    nrm = mesh.normal if normals is None else normals
    P = _stretch(pts, beta)

    def block(m: AeroMesh):
        v = horseshoe_velocity(P, _stretch(m.A, beta), _stretch(m.B, beta))
        # Gamma = U * dCp * dx / 2  ->  w/U = dx/2 * n.v
        return 0.5 * m.chord[None, :] * np.einsum("rk,rsk->rs", nrm, v)

    D = block(mesh)
    if mesh.symmetry in ("symmetric", "antisymmetric"):
        Di = block(mesh.image())
        D = D + Di if mesh.symmetry == "symmetric" else D - Di
    return D


def induced_velocity_matrix(mesh: AeroMesh, points: np.ndarray, mach: float = 0.0):
    """(Npts, N, 3) velocity per unit dCp (for drag / visualisation)."""
    beta = np.sqrt(1.0 - mach ** 2)
    P = _stretch(points, beta)
    v = horseshoe_velocity(P, _stretch(mesh.A, beta), _stretch(mesh.B, beta))
    v = v * 0.5 * mesh.chord[None, :, None]
    if mesh.symmetry in ("symmetric", "antisymmetric"):
        im = mesh.image()
        vi = horseshoe_velocity(P, _stretch(im.A, beta), _stretch(im.B, beta))
        vi = vi * 0.5 * im.chord[None, :, None]
        v = v + vi if mesh.symmetry == "symmetric" else v - vi
    v[..., 0] /= beta
    return v
