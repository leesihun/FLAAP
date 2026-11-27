"""General flat-shell finite elements in 3-D space (quadrilateral and triangle).

Every element is flat: it is given a local Cartesian frame ``(e1, e2, e3)``
(``e3`` = element normal), its nodes are projected onto the element plane,
plate-bending, membrane and drilling matrices are formed in the local frame
and rotated to the global 6-DOF layout ``ux, uy, uz, rx, ry, rz``::

    K_global = T^T K_local T,   T = blockdiag(R, R, ..., R),   R = [e1; e2; e3]

Warped quadrilaterals are projected onto their mean plane (the warping
correction is neglected; keep warp angles small).

Quadrilateral (``quad4``, 24 DOF)
---------------------------------
* **Bending:** MITC4 Reissner-Mindlin (Bathe & Dvorkin 1985), identical
  to :func:`flaap.structures.plate.mitc4_matrices` (assumed covariant
  transverse shear tied at the edge mid-points; no shear locking), but
  with a per-element constitutive matrix so that an orthotropic material
  axis can be defined on curved / arbitrarily oriented shells.
* **Membrane:** QM6 bilinear quadrilateral with Wilson-Taylor incompatible
  modes (``1 - xi^2``, ``1 - eta^2``, Jacobian frozen at the centre,
  statically condensed).  This removes the parasitic shear of the
  bilinear Q4 in in-plane bending, which matters for box/skin structures.
* **Drilling:** Hughes-Brezzi penalty
  ``gamma G t Int (theta_z - 1/2 (v,x - u,y))^2 dA`` (2x2 Gauss), with
  ``gamma = drill`` (default 1e-3).  This gives the drilling DOF a small,
  mesh-consistent stiffness and removes the zero-energy ``theta_z`` modes
  of coplanar meshes without measurably stiffening the element.

Triangle (``tri3``, 18 DOF)
---------------------------
* **Bending:** DKT, discrete Kirchhoff triangle (Batoz, Bathe & Ho 1980).
  The rotation (slope) field ``g = grad w`` is quadratic; it equals the
  nodal slopes at the vertices, at the edge mid-points its tangential
  component is the derivative of the cubic Hermite edge deflection and its
  normal component is linear.  The element is built numerically from this
  definition (3-point Gauss rule, exact).  Thin-shell (Kirchhoff) theory:
  no transverse shear deformation.
* **Membrane:** Allman (1984) triangle with vertex drilling rotations,
  obtained from the 6-node LST by ``u_mid = (u_i + u_j)/2 + (theta_j -
  theta_i)/8 * (y_j - y_i, x_i - x_j)``; plus the Hughes-Brezzi penalty
  at the centroid that removes Allman's spurious constant-rotation mode.

Mass (both elements)
--------------------
Consistent: translational ``rho t N^T N`` (bilinear / linear ``N``; DKT
uses the linear deflection field for its mass) and rotary inertia
``rho t^3 / 12`` on all three rotations (the drilling inertia is
fictitious but keeps the mass matrix frame-invariant).  Lumped: diagonal,
``rho t A / n`` per node and ``rho t^3 A / (12 n)`` per rotation.

Materials
---------
:class:`flaap.structures.plate.OrthotropicMaterial` (plane-stress lamina,
use ``OrthotropicMaterial.isotropic(E, nu, rho)`` for metals).  The
material axis 1 is the projection of the global ``ref_dir`` (default
``+x``) on the element plane, rotated by ``material.angle`` [deg] about
the element normal - for a wing lying in the x-y plane with normals ``+z``
this is the same convention as :mod:`flaap.structures.plate`.

Rotations follow the global right-hand rule; for a plate in the x-y plane
with ``e3 = +z`` the Kirchhoff limit gives ``rx = dw/dy``, ``ry = -dw/dx``,
compatible with :mod:`flaap.structures.plate` and the splines.
"""
from __future__ import annotations

import numpy as np

from .plate import OrthotropicMaterial, _gauss, _jac, _shape, _shear_row

__all__ = ["local_frames", "material_matrices", "shell_quad_matrices", "shell_tri_matrices",
           "dkt_bending", "allman_membrane", "mitc4_bending", "qm6_membrane"]


# ----------------------------------------------------------------- frames
def local_frames(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Local frames of flat elements.

    ``X`` (ne, nn, 3) with nn = 3 or 4 (counter-clockwise about the normal).
    Returns ``R`` (ne, 3, 3) whose rows are ``e1, e2, e3`` and the local
    in-plane coordinates ``xl`` (ne, nn, 2) relative to the centroid.
    ``e3`` is the normal (``(x3-x1) x (x4-x2)`` for quads), ``e1`` points
    from node 1 to node 2 (projected).
    """
    X = np.asarray(X, float)
    nn = X.shape[1]
    if nn == 4:
        n = np.cross(X[:, 2] - X[:, 0], X[:, 3] - X[:, 1])
    else:
        n = np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0])
    nl = np.linalg.norm(n, axis=1)
    if np.any(nl <= 0):
        raise ValueError("degenerate shell element (zero area)")
    e3 = n / nl[:, None]
    a = X[:, 1] - X[:, 0]
    a = a - np.sum(a * e3, 1)[:, None] * e3
    e1 = a / np.linalg.norm(a, axis=1)[:, None]
    e2 = np.cross(e3, e1)
    R = np.stack([e1, e2, e3], axis=1)
    c = X.mean(axis=1)
    xl = np.einsum("eab,enb->ena", R[:, :2, :], X - c[:, None, :])
    return R, xl


def _q_rot(mat: OrthotropicMaterial, theta: np.ndarray) -> np.ndarray:
    """Vectorized :meth:`OrthotropicMaterial.Q` for angles ``theta`` [rad]."""
    d = 1.0 - mat.nu12 * mat.nu21
    Q11, Q22 = mat.E1 / d, mat.E2 / d
    Q12, Q66 = mat.nu12 * mat.E2 / d, mat.G12
    c, s = np.cos(theta), np.sin(theta)
    c2, s2 = c * c, s * s
    q = np.empty(theta.shape + (3, 3))
    q[..., 0, 0] = Q11 * c2 * c2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * s2 * s2
    q[..., 1, 1] = Q11 * s2 * s2 + 2 * (Q12 + 2 * Q66) * s2 * c2 + Q22 * c2 * c2
    q[..., 0, 1] = q[..., 1, 0] = (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (s2 * s2 + c2 * c2)
    q[..., 0, 2] = q[..., 2, 0] = ((Q11 - Q12 - 2 * Q66) * s * c2 * c
                                   + (Q12 - Q22 + 2 * Q66) * s2 * s * c)
    q[..., 1, 2] = q[..., 2, 1] = ((Q11 - Q12 - 2 * Q66) * s2 * s * c
                                   + (Q12 - Q22 + 2 * Q66) * s * c2 * c)
    q[..., 2, 2] = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (s2 * s2 + c2 * c2)
    return q


def material_matrices(mat: OrthotropicMaterial, R: np.ndarray, ref_dir=(1.0, 0.0, 0.0)):
    """Per-element in-plane stiffness ``Q`` (ne,3,3), transverse shear
    ``kappa*Gs`` (ne,2,2) and the material angle in the local frame [rad]."""
    e1, e2, e3 = R[:, 0], R[:, 1], R[:, 2]
    a = np.broadcast_to(np.asarray(ref_dir, float), e1.shape)
    ap = a - np.sum(a * e3, 1)[:, None] * e3
    bad = np.linalg.norm(ap, axis=1) < 1e-3 * np.linalg.norm(a, axis=1)
    if np.any(bad):   # reference direction normal to the element: fall back to +y
        alt = np.array([0.0, 1.0, 0.0])
        ap[bad] = alt - np.sum(alt * e3[bad], 1)[:, None] * e3[bad]
    th = np.arctan2(np.sum(ap * e2, 1), np.sum(ap * e1, 1)) + np.deg2rad(mat.angle)
    Q = _q_rot(mat, th)
    G13 = mat.G12 if mat.G13 is None else mat.G13
    G23 = mat.G12 if mat.G23 is None else mat.G23
    c, s = np.cos(th), np.sin(th)
    Gs = np.empty(th.shape + (2, 2))
    Gs[:, 0, 0] = c * c * G13 + s * s * G23
    Gs[:, 1, 1] = s * s * G13 + c * c * G23
    Gs[:, 0, 1] = Gs[:, 1, 0] = c * s * (G13 - G23)
    return Q, mat.kappa * Gs, th


def _inv2(J):
    det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
    Ji = np.empty_like(J)
    Ji[:, 0, 0] = J[:, 1, 1] / det
    Ji[:, 1, 1] = J[:, 0, 0] / det
    Ji[:, 0, 1] = -J[:, 0, 1] / det
    Ji[:, 1, 0] = -J[:, 1, 0] / det
    return Ji, det


# ----------------------------------------------------------------- quad parts
def mitc4_bending(xl, t, Q, Gs, ngauss: int = 2):
    """MITC4 bending stiffness (ne,12,12), DOFs ``w, phi_x, phi_y`` per node
    (see :func:`flaap.structures.plate.mitc4_matrices`), per-element ``Q``,
    ``Gs`` and nodal thickness ``t`` (ne,4)."""
    ne = xl.shape[0]
    gA = _shear_row(xl, 0.0, -1.0)[:, 0]
    gC = _shear_row(xl, 0.0, 1.0)[:, 0]
    gD = _shear_row(xl, -1.0, 0.0)[:, 1]
    gB = _shear_row(xl, 1.0, 0.0)[:, 1]
    K = np.zeros((ne, 12, 12))
    for xi, eta, w in _gauss(ngauss):
        N, dxi, deta = _shape(xi, eta)
        Ji, det = _inv2(_jac(xl, dxi, deta))
        if np.any(det <= 0):
            raise ValueError("non-positive Jacobian in a shell quad (non-convex or "
                             "inconsistently ordered nodes)")
        dN = np.einsum("eab,bi->eai", Ji, np.stack([dxi, deta]))
        tg = t @ N
        wd = w * det
        Bb = np.zeros((ne, 3, 12))
        Bb[:, 0, 1::3] = dN[:, 0]
        Bb[:, 1, 2::3] = dN[:, 1]
        Bb[:, 2, 1::3] = dN[:, 1]
        Bb[:, 2, 2::3] = dN[:, 0]
        K += np.einsum("e,eai,eab,ebj->eij", wd * tg ** 3 / 12.0, Bb, Q, Bb, optimize=True)
        gxi = 0.5 * (1 - eta) * gA + 0.5 * (1 + eta) * gC
        geta = 0.5 * (1 - xi) * gD + 0.5 * (1 + xi) * gB
        Bs = np.einsum("eab,ebi->eai", Ji, np.stack([gxi, geta], axis=1))
        K += np.einsum("e,eai,eab,ebj->eij", wd * tg, Bs, Gs, Bs, optimize=True)
    return K


def qm6_membrane(xl, t, Q, G, drill: float = 1e-3):
    """QM6 membrane + Hughes-Brezzi drilling for quads.

    Returns (ne,12,12) with DOFs ``u, v, theta_z`` per node.  ``G`` (ne,)
    is the in-plane shear modulus used for the drilling penalty.
    """
    ne = xl.shape[0]
    J0 = _jac(xl, *_shape(0.0, 0.0)[1:])
    J0i, det0 = _inv2(J0)
    Kcc = np.zeros((ne, 8, 8))
    Kci = np.zeros((ne, 8, 4))
    Kii = np.zeros((ne, 4, 4))
    Kd = np.zeros((ne, 12, 12))
    for xi, eta, w in _gauss(2):
        N, dxi, deta = _shape(xi, eta)
        Ji, det = _inv2(_jac(xl, dxi, deta))
        dN = np.einsum("eab,bi->eai", Ji, np.stack([dxi, deta]))
        tg = t @ N
        B = np.zeros((ne, 3, 8))
        B[:, 0, 0::2] = dN[:, 0]
        B[:, 1, 1::2] = dN[:, 1]
        B[:, 2, 0::2] = dN[:, 1]
        B[:, 2, 1::2] = dN[:, 0]
        # incompatible modes P1 = 1 - xi^2, P2 = 1 - eta^2 (Taylor: J0, det0/det)
        dP = np.array([[-2 * xi, 0.0], [0.0, -2 * eta]])        # d/dxi, d/deta (2 modes)
        dPx = np.einsum("eab,bm->eam", J0i, dP) * (det0 / det)[:, None, None]
        Bi = np.zeros((ne, 3, 4))
        Bi[:, 0, 0::2] = dPx[:, 0]
        Bi[:, 1, 1::2] = dPx[:, 1]
        Bi[:, 2, 0::2] = dPx[:, 1]
        Bi[:, 2, 1::2] = dPx[:, 0]
        f = w * det * tg
        Kcc += np.einsum("e,eai,eab,ebj->eij", f, B, Q, B, optimize=True)
        Kci += np.einsum("e,eai,eab,ebj->eij", f, B, Q, Bi, optimize=True)
        Kii += np.einsum("e,eai,eab,ebj->eij", f, Bi, Q, Bi, optimize=True)
        # drilling penalty: theta_z - 1/2 (v,x - u,y)
        bd = np.zeros((ne, 12))
        bd[:, 0::3] = 0.5 * dN[:, 1]          # +1/2 u,y
        bd[:, 1::3] = -0.5 * dN[:, 0]         # -1/2 v,x
        bd[:, 2::3] = N[None, :]
        Kd += (drill * G * f)[:, None, None] * np.einsum("ei,ej->eij", bd, bd)
    Km = Kcc - np.einsum("eij,ejk->eik", Kci, np.linalg.solve(Kii, np.swapaxes(Kci, 1, 2)))
    K = Kd
    idx = np.array([0, 1, 3, 4, 6, 7, 9, 10])
    K[:, idx[:, None], idx[None, :]] += Km
    return K


def _quad_mass(xl, t, rho, lumped):
    """(ne,24,24) local mass, DOFs u v w tx ty tz per node."""
    ne = xl.shape[0]
    M = np.zeros((ne, 24, 24))
    for xi, eta, w in _gauss(3):
        N, dxi, deta = _shape(xi, eta)
        _, det = _inv2(_jac(xl, dxi, deta))
        tg = t @ N
        NN = np.outer(N, N)
        mt = rho * tg * w * det
        mr = rho * tg ** 3 / 12.0 * w * det
        for c in range(6):
            M[:, c::6, c::6] += (mt if c < 3 else mr)[:, None, None] * NN
    if lumped:
        d = M.sum(axis=2)
        M = np.zeros_like(M)
        M[:, np.arange(24), np.arange(24)] = d
    return M


# ----------------------------------------------------------------- triangle parts
_TRI3_GP = (np.array([[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]]),
            np.full(3, 1 / 3))          # area coordinates, weights (x area)
_EDGES = ((0, 1), (1, 2), (2, 0))       # midside nodes 3, 4, 5


def _lst_dN(L, xl):
    """Derivatives of the 6 quadratic triangle shape functions at area
    coordinates ``L`` (3,), straight-sided element; returns (ne,2,6) and area."""
    x, y = xl[:, :, 0], xl[:, :, 1]
    A2 = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1) / A2[:, None]
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1) / A2[:, None]
    # N_i = L_i (2 L_i - 1), N_{3+k} = 4 L_i L_j for edge k = (i, j)
    dNdL = np.zeros((6, 3))
    for i in range(3):
        dNdL[i, i] = 4 * L[i] - 1
    for k, (i, j) in enumerate(_EDGES):
        dNdL[3 + k, i] = 4 * L[j]
        dNdL[3 + k, j] = 4 * L[i]
    dx = np.einsum("ni,ei->en", dNdL, b)
    dy = np.einsum("ni,ei->en", dNdL, c)
    return np.stack([dx, dy], 1), 0.5 * A2


def dkt_bending(xl, t, Q):
    """DKT bending stiffness (ne,9,9); DOFs ``w, rx, ry`` per node with the
    Kirchhoff relations ``rx = w,y``, ``ry = -w,x`` (local frame)."""
    ne = xl.shape[0]
    # C: nodal dofs (9) -> slope field g = grad w at the 6 LST nodes (6 x 2 = 12)
    C = np.zeros((ne, 6, 2, 9))
    for i in range(3):
        C[:, i, 0, 3 * i + 2] = -1.0     # w,x = -ry
        C[:, i, 1, 3 * i + 1] = 1.0      # w,y = rx
    for k, (i, j) in enumerate(_EDGES):
        d = xl[:, j] - xl[:, i]
        Lk = np.linalg.norm(d, axis=1)
        s = d / Lk[:, None]
        n = np.stack([s[:, 1], -s[:, 0]], 1)
        gi, gj = C[:, i], C[:, j]                              # (ne,2,9)
        gs_i = np.einsum("ea,eaq->eq", s, gi)
        gs_j = np.einsum("ea,eaq->eq", s, gj)
        gn = 0.5 * (np.einsum("ea,eaq->eq", n, gi) + np.einsum("ea,eaq->eq", n, gj))
        gs = -0.25 * (gs_i + gs_j)
        gs[:, 3 * j] += 1.5 / Lk
        gs[:, 3 * i] -= 1.5 / Lk
        C[:, 3 + k] = s[:, :, None] * gs[:, None, :] + n[:, :, None] * gn[:, None, :]
    K = np.zeros((ne, 9, 9))
    for L, w in zip(*_TRI3_GP):
        dN, A = _lst_dN(L, xl)
        gx, gy = C[:, :, 0, :], C[:, :, 1, :]                 # (ne,6,9)
        B = np.stack([np.einsum("en,enq->eq", dN[:, 0], gx),
                      np.einsum("en,enq->eq", dN[:, 1], gy),
                      np.einsum("en,enq->eq", dN[:, 1], gx) + np.einsum("en,enq->eq", dN[:, 0], gy)],
                     1)
        K += np.einsum("e,eai,eab,ebj->eij", w * A * t ** 3 / 12.0, B, Q, B, optimize=True)
    return K


def allman_membrane(xl, t, Q, G, drill: float = 1e-3):
    """Allman membrane triangle with drilling DOFs (ne,9,9); DOFs
    ``u, v, theta_z`` per node; plus the centroid Hughes-Brezzi penalty."""
    ne = xl.shape[0]
    T = np.zeros((ne, 12, 9))           # LST dofs (u,v at 6 nodes) <- (u,v,tz) x 3
    for i in range(3):
        T[:, 2 * i, 3 * i] = 1.0
        T[:, 2 * i + 1, 3 * i + 1] = 1.0
    for k, (i, j) in enumerate(_EDGES):
        m = 3 + k
        dx = xl[:, j, 0] - xl[:, i, 0]
        dy = xl[:, j, 1] - xl[:, i, 1]
        for a in (i, j):
            T[:, 2 * m, 3 * a] = 0.5
            T[:, 2 * m + 1, 3 * a + 1] = 0.5
        T[:, 2 * m, 3 * j + 2] += dy / 8.0
        T[:, 2 * m, 3 * i + 2] -= dy / 8.0
        T[:, 2 * m + 1, 3 * j + 2] -= dx / 8.0
        T[:, 2 * m + 1, 3 * i + 2] += dx / 8.0
    K = np.zeros((ne, 9, 9))
    for L, w in zip(*_TRI3_GP):
        dN, A = _lst_dN(L, xl)
        B = np.zeros((ne, 3, 12))
        B[:, 0, 0::2] = dN[:, 0]
        B[:, 1, 1::2] = dN[:, 1]
        B[:, 2, 0::2] = dN[:, 1]
        B[:, 2, 1::2] = dN[:, 0]
        Bt = np.einsum("eai,eij->eaj", B, T)
        K += np.einsum("e,eai,eab,ebj->eij", w * A * t, Bt, Q, Bt, optimize=True)
    dN, A = _lst_dN(np.full(3, 1 / 3), xl)
    rot = np.zeros((ne, 12))
    rot[:, 0::2] = 0.5 * dN[:, 1]
    rot[:, 1::2] = -0.5 * dN[:, 0]
    bd = np.einsum("ei,eij->ej", rot, T)
    bd[:, 2::3] += 1.0 / 3.0
    K += (drill * G * A * t)[:, None, None] * np.einsum("ei,ej->eij", bd, bd)
    return K


def _tri_mass(xl, t, rho, lumped):
    ne = xl.shape[0]
    x, y = xl[:, :, 0], xl[:, :, 1]
    A = 0.5 * ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    if lumped:
        NN = np.eye(3) / 3.0
    else:
        NN = (np.ones((3, 3)) + np.eye(3)) / 12.0
    M = np.zeros((ne, 18, 18))
    mt = rho * t * A
    mr = rho * t ** 3 / 12.0 * A
    for c in range(6):
        M[:, c::6, c::6] = (mt if c < 3 else mr)[:, None, None] * NN
    return M


# ----------------------------------------------------------------- global elements
def _rotate(Kl, R):
    nn = Kl.shape[1] // 6
    T = np.zeros((R.shape[0], 6 * nn, 6 * nn))
    for a in range(2 * nn):
        T[:, 3 * a:3 * a + 3, 3 * a:3 * a + 3] = R
    return np.einsum("eki,ekl,elj->eij", T, Kl, T, optimize=True)


def _shear_modulus(Q):
    return Q[:, 2, 2]


def shell_quad_matrices(X, t, mat: OrthotropicMaterial, lumped: bool = False,
                        drill: float = 1e-3, ref_dir=(1.0, 0.0, 0.0), ngauss: int = 2):
    """Global stiffness and mass of flat shell quads.

    ``X`` (ne,4,3) node coordinates, ``t`` scalar, (ne,) or (ne,4) nodal
    thickness.  Returns ``K, M`` (ne,24,24) in the global DOF order
    ``ux uy uz rx ry rz`` per node.
    """
    X = np.asarray(X, float)
    ne = X.shape[0]
    t = np.asarray(t, float)
    t = np.broadcast_to(t[:, None] if t.ndim == 1 else t, (ne, 4)).astype(float)
    R, xl = local_frames(X)
    Q, Gs, _ = material_matrices(mat, R, ref_dir)
    Kl = np.zeros((ne, 24, 24))
    Kb = mitc4_bending(xl, t, Q, Gs, ngauss)
    # (w, phi_x, phi_y) -> (w, ty, -tx)
    bidx = np.array([[6 * a + 2, 6 * a + 4, 6 * a + 3] for a in range(4)]).ravel()
    sg = np.tile([1.0, 1.0, -1.0], 4)
    Kl[:, bidx[:, None], bidx[None, :]] += Kb * sg[None, :, None] * sg[None, None, :]
    Km = qm6_membrane(xl, t, Q, _shear_modulus(Q), drill)
    midx = np.array([[6 * a, 6 * a + 1, 6 * a + 5] for a in range(4)]).ravel()
    Kl[:, midx[:, None], midx[None, :]] += Km
    Ml = _quad_mass(xl, t, mat.rho, lumped)
    return _rotate(Kl, R), _rotate(Ml, R)


def shell_tri_matrices(X, t, mat: OrthotropicMaterial, lumped: bool = False,
                       drill: float = 1e-3, ref_dir=(1.0, 0.0, 0.0)):
    """Global stiffness and mass of flat shell triangles (DKT + Allman).

    ``X`` (ne,3,3), ``t`` scalar, (ne,) or (ne,3) (the element mean is
    used).  Returns ``K, M`` (ne,18,18), DOFs ``ux uy uz rx ry rz``.
    """
    X = np.asarray(X, float)
    ne = X.shape[0]
    t = np.asarray(t, float)
    t = np.broadcast_to(t.mean(axis=1) if t.ndim == 2 else t, (ne,)).astype(float)
    R, xl = local_frames(X)
    Q, _, _ = material_matrices(mat, R, ref_dir)
    Kl = np.zeros((ne, 18, 18))
    bidx = np.array([[6 * a + 2, 6 * a + 3, 6 * a + 4] for a in range(3)]).ravel()
    Kl[:, bidx[:, None], bidx[None, :]] += dkt_bending(xl, t, Q)
    midx = np.array([[6 * a, 6 * a + 1, 6 * a + 5] for a in range(3)]).ravel()
    Kl[:, midx[:, None], midx[None, :]] += allman_membrane(xl, t, Q, _shear_modulus(Q), drill)
    Ml = _tri_mass(xl, t, mat.rho, lumped)
    return _rotate(Kl, R), _rotate(Ml, R)
