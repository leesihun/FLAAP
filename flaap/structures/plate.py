"""Flat-plate finite elements for wings lying in the x-y plane.

Element
-------
MITC4 Reissner-Mindlin bending quadrilateral (Bathe & Dvorkin 1985):
bilinear ``w``, ``phi_x``, ``phi_y`` with the assumed (mixed-interpolated)
covariant transverse-shear field tied at the edge mid-points, which removes
shear locking for thin plates.  Material: orthotropic plane-stress lamina
whose principal axis 1 lies in the plate plane at an angle ``angle`` (deg)
measured from +x towards +y (counter-clockwise seen from +z).  For a swept
right wing with sweep ``L`` (positive aft) the spanwise line has
``angle = 90 - L``; see :meth:`OrthotropicMaterial.along_sweep`.

Thickness is given per node and interpolated bilinearly; bending
(``D = t^3/12 Qbar``), transverse shear (``k t G``) and consistent mass
(translational ``rho t`` and rotary ``rho t^3/12``) are integrated with a
3x3 Gauss rule, so tapered / airfoil-shaped thickness is captured inside an
element.  An optional bilinear plane-stress membrane (``ux``, ``uy``) can be
included; by default the in-plane DOFs are simply constrained.

DOF / sign conventions (FLAAP FE frame)
---------------------------------------
Global node DOFs are ``ux, uy, uz, rx, ry, rz`` with right-handed rotations
about the axes.  Mindlin section rotations are ``u = z phi_x``,
``v = z phi_y``, hence::

    ry = +phi_x      (Kirchhoff limit: ry = -dw/dx, nose-up pitch positive,
                      uz = -(x - x0) ry)
    rx = -phi_y      (Kirchhoff limit: rx = +dw/dy)

so the plate rotations are directly compatible with
:func:`flaap.structures.modal.rigid_body_shapes` and the beam spline.
The drilling rotation ``rz`` is always constrained.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from ..geometry import _spacing
from .modal import FEModel

__all__ = ["OrthotropicMaterial", "plate_mesh_trapezoid", "build_plate_model",
           "mitc4_matrices", "membrane_matrices"]


# ----------------------------------------------------------------- material
@dataclass
class OrthotropicMaterial:
    """Plane-stress orthotropic material.

    ``angle`` [deg] is the direction of axis 1 in the x-y plane measured
    from +x towards +y.  ``G13``/``G23`` (transverse shear moduli) default
    to ``G12``; ``kappa`` is the shear correction factor.
    """

    E1: float
    E2: float
    G12: float
    nu12: float
    rho: float
    angle: float = 0.0
    G13: float | None = None
    G23: float | None = None
    kappa: float = 5.0 / 6.0

    @classmethod
    def isotropic(cls, E: float, nu: float, rho: float, **kw) -> "OrthotropicMaterial":
        G = E / (2 * (1 + nu))
        return cls(E, E, G, nu, rho, **kw)

    @staticmethod
    def along_sweep(sweep_deg: float) -> float:
        """Material angle (deg from +x) of a line swept back by ``sweep_deg``
        on the right (y > 0) wing."""
        return 90.0 - sweep_deg

    @property
    def nu21(self) -> float:
        return self.nu12 * self.E2 / self.E1

    def Q(self) -> np.ndarray:
        """Reduced plane-stress stiffness (3x3, engineering shear) in x-y axes."""
        d = 1.0 - self.nu12 * self.nu21
        Q11, Q22 = self.E1 / d, self.E2 / d
        Q12, Q66 = self.nu12 * self.E2 / d, self.G12
        th = np.deg2rad(self.angle)
        c, s = np.cos(th), np.sin(th)
        c2, s2 = c * c, s * s
        q11 = Q11 * c2 * c2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * s2 * s2
        q22 = Q11 * s2 * s2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * c2 * c2
        q12 = (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2)
        q16 = (Q11 - Q12 - 2 * Q66) * s * c2 * c + (Q12 - Q22 + 2 * Q66) * s2 * s * c
        q26 = (Q11 - Q12 - 2 * Q66) * s2 * s * c + (Q12 - Q22 + 2 * Q66) * s * c2 * c
        q66 = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (s2 * s2 + c2 * c2)
        return np.array([[q11, q12, q16], [q12, q22, q26], [q16, q26, q66]])

    def Gs(self) -> np.ndarray:
        """Transverse shear modulus matrix (2x2) in x-y axes (without kappa)."""
        G13 = self.G12 if self.G13 is None else self.G13
        G23 = self.G12 if self.G23 is None else self.G23
        th = np.deg2rad(self.angle)
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        return R @ np.diag([G13, G23]) @ R.T


# ----------------------------------------------------------------- shape fns
_XI = np.array([-1.0, 1.0, 1.0, -1.0])
_ETA = np.array([-1.0, -1.0, 1.0, 1.0])


def _shape(xi, eta):
    N = 0.25 * (1 + _XI * xi) * (1 + _ETA * eta)
    dxi = 0.25 * _XI * (1 + _ETA * eta)
    deta = 0.25 * _ETA * (1 + _XI * xi)
    return N, dxi, deta


def _gauss(n):
    g, w = np.polynomial.legendre.leggauss(n)
    return [(a, b, wa * wb) for a, wa in zip(g, w) for b, wb in zip(g, w)]


def _jac(X, dxi, deta):
    """X (ne,4,2) -> J (ne,2,2) with rows [x_xi, y_xi], [x_eta, y_eta]."""
    J = np.empty((X.shape[0], 2, 2))
    J[:, 0, :] = np.einsum("i,eik->ek", dxi, X)
    J[:, 1, :] = np.einsum("i,eik->ek", deta, X)
    return J


def _shear_row(X, xi, eta):
    """Covariant shear strains (gamma_xi, gamma_eta) at (xi, eta): (ne,2,12)."""
    N, dxi, deta = _shape(xi, eta)
    J = _jac(X, dxi, deta)
    ne = X.shape[0]
    B = np.zeros((ne, 2, 12))
    for a, dN in ((0, dxi), (1, deta)):
        B[:, a, 0::3] = dN[None, :]
        B[:, a, 1::3] = N[None, :] * J[:, a, 0:1]
        B[:, a, 2::3] = N[None, :] * J[:, a, 1:2]
    return B


def mitc4_matrices(X: np.ndarray, t: np.ndarray, mat: OrthotropicMaterial,
                   ngauss: int = 3):
    """MITC4 bending element matrices.

    ``X`` (ne,4,2) nodal x-y coordinates (counter-clockwise), ``t`` (ne,4)
    nodal thickness.  Element DOF order per node: ``w, phi_x, phi_y``.
    Returns ``Ke, Me`` of shape (ne,12,12).
    """
    X = np.asarray(X, float)
    t = np.asarray(t, float)
    ne = X.shape[0]
    Q = mat.Q()
    Gs = mat.kappa * mat.Gs()
    # shear tying points: gamma_xi at A(0,-1), C(0,1); gamma_eta at D(-1,0), B(1,0)
    gA = _shear_row(X, 0.0, -1.0)[:, 0]
    gC = _shear_row(X, 0.0, 1.0)[:, 0]
    gD = _shear_row(X, -1.0, 0.0)[:, 1]
    gB = _shear_row(X, 1.0, 0.0)[:, 1]
    K = np.zeros((ne, 12, 12))
    M = np.zeros((ne, 12, 12))
    for xi, eta, w in _gauss(ngauss):
        N, dxi, deta = _shape(xi, eta)
        J = _jac(X, dxi, deta)
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        if np.any(det <= 0):
            raise ValueError("non-positive Jacobian: quads must be counter-clockwise "
                             "(seen from +z) and convex")
        Ji = np.empty_like(J)
        Ji[:, 0, 0] = J[:, 1, 1] / det
        Ji[:, 1, 1] = J[:, 0, 0] / det
        Ji[:, 0, 1] = -J[:, 0, 1] / det
        Ji[:, 1, 0] = -J[:, 1, 0] / det
        dN = np.einsum("eab,bi->eai", Ji, np.stack([dxi, deta]))  # (ne,2,4) d/dx, d/dy
        tg = t @ N
        wd = w * det
        # bending
        Bb = np.zeros((ne, 3, 12))
        Bb[:, 0, 1::3] = dN[:, 0]
        Bb[:, 1, 2::3] = dN[:, 1]
        Bb[:, 2, 1::3] = dN[:, 1]
        Bb[:, 2, 2::3] = dN[:, 0]
        K += np.einsum("e,eai,ab,ebj->eij", wd * tg ** 3 / 12.0, Bb, Q, Bb)
        # assumed transverse shear
        gxi = 0.5 * (1 - eta) * gA + 0.5 * (1 + eta) * gC
        geta = 0.5 * (1 - xi) * gD + 0.5 * (1 + xi) * gB
        Bs = np.einsum("eab,ebi->eai", Ji, np.stack([gxi, geta], axis=1))
        K += np.einsum("e,eai,ab,ebj->eij", wd * tg, Bs, Gs, Bs)
        # consistent mass
        NN = np.outer(N, N)
        mt = mat.rho * tg * wd
        mr = mat.rho * tg ** 3 / 12.0 * wd
        M[:, 0::3, 0::3] += mt[:, None, None] * NN
        M[:, 1::3, 1::3] += mr[:, None, None] * NN
        M[:, 2::3, 2::3] += mr[:, None, None] * NN
    return K, M


def membrane_matrices(X: np.ndarray, t: np.ndarray, mat: OrthotropicMaterial,
                      ngauss: int = 2):
    """Bilinear plane-stress membrane quad; DOFs ``ux, uy`` per node.
    Returns ``Ke, Me`` (ne,8,8)."""
    X = np.asarray(X, float)
    t = np.asarray(t, float)
    ne = X.shape[0]
    Q = mat.Q()
    K = np.zeros((ne, 8, 8))
    M = np.zeros((ne, 8, 8))
    for xi, eta, w in _gauss(ngauss):
        N, dxi, deta = _shape(xi, eta)
        J = _jac(X, dxi, deta)
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        Ji = np.empty_like(J)
        Ji[:, 0, 0] = J[:, 1, 1] / det
        Ji[:, 1, 1] = J[:, 0, 0] / det
        Ji[:, 0, 1] = -J[:, 0, 1] / det
        Ji[:, 1, 0] = -J[:, 1, 0] / det
        dN = np.einsum("eab,bi->eai", Ji, np.stack([dxi, deta]))
        tg = t @ N
        B = np.zeros((ne, 3, 8))
        B[:, 0, 0::2] = dN[:, 0]
        B[:, 1, 1::2] = dN[:, 1]
        B[:, 2, 0::2] = dN[:, 1]
        B[:, 2, 1::2] = dN[:, 0]
        K += np.einsum("e,eai,ab,ebj->eij", w * det * tg, B, Q, B)
        mt = mat.rho * tg * w * det
        NN = np.outer(N, N)
        M[:, 0::2, 0::2] += mt[:, None, None] * NN
        M[:, 1::2, 1::2] += mt[:, None, None] * NN
    return K, M


# ----------------------------------------------------------------- meshing
def plate_mesh_trapezoid(root_le, root_chord: float, tip_le, tip_chord: float,
                         n_chord: int, n_span: int, chord_spacing="linear",
                         span_spacing="linear"):
    """Structured quad mesh of a trapezoidal planform.

    ``root_le`` / ``tip_le`` are (x, y[, z]) leading-edge points.  Node
    ``j*(n_chord+1) + i`` is at span station ``j`` (root = 0) and chord
    station ``i`` (LE = 0).  Quads are counter-clockwise seen from +z for
    a right wing (tip y > root y).  Returns ``nodes`` (Nn,3), ``quads``
    (ne,4).
    """
    r = np.zeros(3)
    tp = np.zeros(3)
    r[:len(root_le)] = root_le
    tp[:len(tip_le)] = tip_le
    eta = _spacing(n_span, span_spacing)
    xi = _spacing(n_chord, chord_spacing)
    le = r[None, :] + eta[:, None] * (tp - r)[None, :]
    ch = root_chord + eta * (tip_chord - root_chord)
    nodes = np.repeat(le[:, None, :], n_chord + 1, axis=1)
    nodes[:, :, 0] += ch[:, None] * xi[None, :]
    nodes = nodes.reshape(-1, 3)
    nc1 = n_chord + 1
    j, i = np.meshgrid(np.arange(n_span), np.arange(n_chord), indexing="ij")
    n0 = (j * nc1 + i).ravel()
    quads = np.stack([n0, n0 + 1, n0 + 1 + nc1, n0 + nc1], axis=1)
    if tp[1] < r[1]:                       # left wing: keep counter-clockwise
        quads = quads[:, ::-1]
    return nodes, quads


# ----------------------------------------------------------------- assembly
_BEND_DOF = np.array([2, 4, 3])            # element (w, phi_x, phi_y) -> (uz, ry, rx)
_BEND_SGN = np.array([1.0, 1.0, -1.0])     # phi_y = -rx


def _assemble(n, dofs, Ke):
    ne, nd = dofs.shape
    rows = np.repeat(dofs, nd, axis=1).ravel()
    cols = np.tile(dofs, (1, nd)).ravel()
    return sp.coo_matrix((Ke.ravel(), (rows, cols)), shape=(n, n)).tocsr()


def build_plate_model(nodes: np.ndarray, quads: np.ndarray, thickness_at_nodes,
                      material: OrthotropicMaterial, clamp_y: float | None = 0.0,
                      membrane: bool = False, ngauss: int = 3, clamp_tol: float = 1e-6,
                      name: str = "plate") -> FEModel:
    """Assemble an MITC4 plate model into an :class:`FEModel`.

    ``thickness_at_nodes`` scalar or (Nn,).  All 6 DOFs of nodes with
    ``|y - clamp_y| < clamp_tol * size`` are fixed (``clamp_y=None``: free).
    ``rz`` is always fixed; ``ux, uy`` are fixed unless ``membrane=True``.
    """
    nodes = np.asarray(nodes, float)
    quads = np.asarray(quads, int)
    nn = nodes.shape[0]
    t = np.broadcast_to(np.asarray(thickness_at_nodes, float), (nn,)).copy()
    X = nodes[quads][:, :, :2]
    te = t[quads]
    Kb, Mb = mitc4_matrices(X, te, material, ngauss)
    sgn = np.tile(_BEND_SGN, 4)
    Kb *= sgn[None, :, None] * sgn[None, None, :]
    Mb *= sgn[None, :, None] * sgn[None, None, :]
    dofs = (6 * quads[:, :, None] + _BEND_DOF[None, None, :]).reshape(-1, 12)
    ndof = 6 * nn
    K = _assemble(ndof, dofs, Kb)
    M = _assemble(ndof, dofs, Mb)
    fixed = [6 * np.arange(nn) + 5]
    if membrane:
        Km, Mm = membrane_matrices(X, te, material)
        dm = (6 * quads[:, :, None] + np.array([0, 1])[None, None, :]).reshape(-1, 8)
        K = K + _assemble(ndof, dm, Km)
        M = M + _assemble(ndof, dm, Mm)
    else:
        fixed += [6 * np.arange(nn), 6 * np.arange(nn) + 1]
    if clamp_y is not None:
        size = max(np.ptp(nodes[:, 0]), np.ptp(nodes[:, 1]), 1e-12)
        cl = np.where(np.abs(nodes[:, 1] - clamp_y) < clamp_tol * size)[0]
        if cl.size == 0:
            raise ValueError(f"no nodes found at clamp_y={clamp_y}")
        fixed.append((6 * cl[:, None] + np.arange(6)[None, :]).ravel())
    fixed = np.unique(np.concatenate(fixed))
    K = 0.5 * (K + K.T)
    M = 0.5 * (M + M.T)
    return FEModel(nodes, K.tocsr(), M.tocsr(), fixed, name=name,
                   elements={"quad4": quads, "thickness": t, "material": material})
