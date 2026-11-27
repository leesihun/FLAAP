"""3-D continuum (solid) finite elements for modal analysis.

Elements (isotropic linear elasticity, consistent mass, vectorized over
elements):

``tet4``
    linear tetrahedron (constant strain).  Exact 1-point stiffness,
    exact consistent mass ``rho V / 20 (1 + delta_ij)``.  Very stiff in
    bending (locks): use only for chunky parts or as a fallback.
``tet10``
    quadratic tetrahedron (node order of VTK / meshio / Radioss:
    corners 0-3, mid-edge nodes on 01, 12, 02, 03, 13, 23).  Curved
    (iso-parametric) edges allowed.  Stiffness with the 4-point (degree 2)
    rule, consistent mass with a 64-point collapsed-Gauss (conical
    product) rule (exact for straight edges).  Recommended for meshes of
    general geometry (gmsh / STL); converges quickly in bending.
``hex8``
    trilinear hexahedron with Wilson-Taylor incompatible modes (``1 -
    xi^2`` etc., 9 internal DOFs statically condensed, Jacobian frozen at
    the centroid so the patch test is passed), 2x2x2 Gauss.  The
    incompatible modes remove shear locking in bending, so a single
    element through the thickness bends almost exactly.
    ``incompatible=False`` gives the plain (locking) trilinear brick.

Solid nodes carry translations only.  In the 6-DOF FLAAP layout the
rotation DOFs of solid-only nodes have no stiffness and no mass and are
removed before the eigen-solution (see :mod:`flaap.structures.fem`).
Nodal rotations of the mode shapes are then recovered with
:func:`rotations_from_translations` (least-squares displacement gradient,
rotation = skew part = 1/2 curl u) or set to zero; see
:func:`flaap.structures.fem.solve_modes`.

Voigt order: ``[xx, yy, zz, xy, yz, zx]`` with engineering shear strains.
"""
from __future__ import annotations

import numpy as np

__all__ = ["iso_D", "tet4_matrices", "tet10_matrices", "hex8_matrices", "fix_tet10_order",
           "TET10_EDGES", "box_mesh", "tet4_to_tet10", "tet_volumes", "rotations_from_translations",
           "SOLID_FACES"]

TET10_EDGES = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))

# boundary faces of each element type (outward normals for positive volume)
SOLID_FACES = {
    "tet4": ((0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)),
    "hex8": ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}
SOLID_FACES["tet10"] = SOLID_FACES["tet4"]


def iso_D(E: float, nu: float) -> np.ndarray:
    """Isotropic 6x6 elasticity matrix (engineering shear strains)."""
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    G = E / (2 * (1 + nu))
    D = np.zeros((6, 6))
    D[:3, :3] = lam
    D[np.arange(3), np.arange(3)] += 2 * G
    D[np.arange(3, 6), np.arange(3, 6)] = G
    return D


def _B(dN):
    """dN (ne,3,nn) cartesian derivatives -> B (ne,6,3nn)."""
    ne, _, nn = dN.shape
    B = np.zeros((ne, 6, 3 * nn))
    B[:, 0, 0::3] = dN[:, 0]
    B[:, 1, 1::3] = dN[:, 1]
    B[:, 2, 2::3] = dN[:, 2]
    B[:, 3, 0::3] = dN[:, 1]
    B[:, 3, 1::3] = dN[:, 0]
    B[:, 4, 1::3] = dN[:, 2]
    B[:, 4, 2::3] = dN[:, 1]
    B[:, 5, 0::3] = dN[:, 2]
    B[:, 5, 2::3] = dN[:, 0]
    return B


def _expand_mass(Ms):
    """scalar mass (ne,nn,nn) -> (ne,3nn,3nn)."""
    ne, nn, _ = Ms.shape
    M = np.zeros((ne, 3 * nn, 3 * nn))
    for c in range(3):
        M[:, c::3, c::3] = Ms
    return M


def tet_volumes(X):
    """Signed volumes of tetrahedra from their 4 corner nodes (ne,>=4,3)."""
    X = np.asarray(X, float)
    return np.einsum("ei,ei->e", np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]),
                     X[:, 3] - X[:, 0]) / 6.0


# ----------------------------------------------------------------- tet4
def tet4_matrices(X, E: float, nu: float, rho: float, lumped: bool = False):
    """Linear tetrahedron: ``K, M`` (ne,12,12), DOFs ux uy uz per node."""
    X = np.asarray(X, float)
    J = np.stack([X[:, 1] - X[:, 0], X[:, 2] - X[:, 0], X[:, 3] - X[:, 0]], axis=1)  # rows
    det = np.linalg.det(J)
    if np.any(det <= 0):
        raise ValueError("tet4 with non-positive volume (check node ordering)")
    Ji = np.linalg.inv(J)                         # d(xi)/dx: columns
    dNdxi = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], float)
    dN = np.einsum("eab,bi->eai", Ji, dNdxi)
    B = _B(dN)
    V = det / 6.0
    K = np.einsum("e,eai,ab,ebj->eij", V, B, iso_D(E, nu), B, optimize=True)
    Ms = (np.ones((4, 4)) + np.eye(4)) / 20.0 if not lumped else np.eye(4) / 4.0
    M = _expand_mass(rho * V[:, None, None] * Ms[None])
    return K, M


# ----------------------------------------------------------------- tet10
def _tet10_N(L):
    """Shape functions and d/d(xi,eta,zeta) at barycentric L (4,), with
    xi = L1, eta = L2, zeta = L3 and L0 = 1 - xi - eta - zeta."""
    N = np.empty(10)
    dL = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], float)  # dL_k/dxi_a
    dN = np.zeros((3, 10))
    for i in range(4):
        N[i] = L[i] * (2 * L[i] - 1)
        dN[:, i] = (4 * L[i] - 1) * dL[:, i]
    for k, (i, j) in enumerate(TET10_EDGES):
        N[4 + k] = 4 * L[i] * L[j]
        dN[:, 4 + k] = 4 * (L[i] * dL[:, j] + L[j] * dL[:, i])
    return N, dN


def _tet_rule(order: str):
    if order == "4":
        a, b = 0.5854101966249685, 0.1381966011250105
        L = np.array([[a, b, b, b], [b, a, b, b], [b, b, a, b], [b, b, b, a]])
        return L, np.full(4, 1 / 24)
    # collapsed Gauss (conical product), 4^3 points: x = u, y = v(1-u), z = w(1-u)(1-v)
    g, w = np.polynomial.legendre.leggauss(4)
    g = 0.5 * (g + 1)
    w = 0.5 * w
    pts, wts = [], []
    for u, wu in zip(g, w):
        for v, wv in zip(g, w):
            for s, ws in zip(g, w):
                x, y, z = u, v * (1 - u), s * (1 - u) * (1 - v)
                pts.append([1 - x - y - z, x, y, z])
                wts.append(wu * wv * ws * (1 - u) ** 2 * (1 - v))
    return np.array(pts), np.array(wts)


def tet10_matrices(X, E: float, nu: float, rho: float, lumped: bool = False):
    """Quadratic tetrahedron: ``K, M`` (ne,30,30), DOFs ux uy uz per node.

    Lumped mass uses the HRZ (diagonal scaling) procedure.
    """
    X = np.asarray(X, float)
    ne = X.shape[0]
    D = iso_D(E, nu)
    K = np.zeros((ne, 30, 30))
    for L, w in zip(*_tet_rule("4")):
        _, dNr = _tet10_N(L)
        J = np.einsum("ai,eik->eak", dNr, X)       # rows d x/d xi_a
        det = np.linalg.det(J)
        if np.any(det <= 0):
            raise ValueError("tet10 with non-positive Jacobian (check node order: "
                             "corners, then mid-edges 01 12 02 03 13 23)")
        dN = np.linalg.solve(J, np.broadcast_to(dNr, (ne, 3, 10)))
        B = _B(dN)
        K += np.einsum("e,eai,ab,ebj->eij", w * det, B, D, B, optimize=True)
    Ms = np.zeros((ne, 10, 10))
    for L, w in zip(*_tet_rule("mass")):
        N, dNr = _tet10_N(L)
        det = np.linalg.det(np.einsum("ai,eik->eak", dNr, X))
        Ms += (rho * w * det)[:, None, None] * np.outer(N, N)[None]
    if lumped:   # HRZ
        tot = Ms.sum(axis=(1, 2))
        d = np.einsum("eii->ei", Ms)
        d = d * (tot / d.sum(axis=1))[:, None]
        Ms = np.zeros_like(Ms)
        Ms[:, np.arange(10), np.arange(10)] = d
    return K, _expand_mass(Ms)


def fix_tet10_order(nodes: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """Reorder the mid-edge nodes of tet10 connectivity to :data:`TET10_EDGES`.

    Works for any input convention (gmsh's 8/9 swap, Abaqus, ...): each
    mid-edge slot is assigned the node nearest to that edge's midpoint.
    Also flips negatively oriented elements.
    """
    conn = np.asarray(conn, int).copy()
    X = nodes[conn]
    neg = tet_volumes(X) < 0
    if np.any(neg):     # swap corners 1 and 2 (and their mid-edges via the matching below)
        conn[neg, 1], conn[neg, 2] = conn[neg, 2].copy(), conn[neg, 1].copy()
        X = nodes[conn]
    mids = X[:, 4:10]                                       # (ne,6,3) candidates
    out = conn.copy()
    for k, (i, j) in enumerate(TET10_EDGES):
        m = 0.5 * (X[:, i] + X[:, j])
        d = np.linalg.norm(mids - m[:, None, :], axis=2)
        out[:, 4 + k] = conn[np.arange(conn.shape[0]), 4 + np.argmin(d, axis=1)]
    return out


# ----------------------------------------------------------------- hex8
_HXI = np.array([-1, 1, 1, -1, -1, 1, 1, -1], float)
_HETA = np.array([-1, -1, 1, 1, -1, -1, 1, 1], float)
_HZETA = np.array([-1, -1, -1, -1, 1, 1, 1, 1], float)


def _hex_N(a, b, c):
    N = 0.125 * (1 + _HXI * a) * (1 + _HETA * b) * (1 + _HZETA * c)
    dN = 0.125 * np.stack([_HXI * (1 + _HETA * b) * (1 + _HZETA * c),
                           _HETA * (1 + _HXI * a) * (1 + _HZETA * c),
                           _HZETA * (1 + _HXI * a) * (1 + _HETA * b)])
    return N, dN


def hex8_matrices(X, E: float, nu: float, rho: float, lumped: bool = False,
                  incompatible: bool = True):
    """Trilinear brick (optionally with incompatible modes): ``K, M`` (ne,24,24)."""
    X = np.asarray(X, float)
    ne = X.shape[0]
    D = iso_D(E, nu)
    g, wg = np.polynomial.legendre.leggauss(2)
    _, dN0 = _hex_N(0.0, 0.0, 0.0)
    J0 = np.einsum("ai,eik->eak", dN0, X)
    det0 = np.linalg.det(J0)
    Kcc = np.zeros((ne, 24, 24))
    Kci = np.zeros((ne, 24, 9))
    Kii = np.zeros((ne, 9, 9))
    Ms = np.zeros((ne, 8, 8))
    for a, wa in zip(g, wg):
        for b, wb in zip(g, wg):
            for c, wc in zip(g, wg):
                N, dNr = _hex_N(a, b, c)
                J = np.einsum("ai,eik->eak", dNr, X)
                det = np.linalg.det(J)
                if np.any(det <= 0):
                    raise ValueError("hex8 with non-positive Jacobian (check node ordering)")
                dN = np.linalg.solve(J, np.broadcast_to(dNr, (ne, 3, 8)))
                B = _B(dN)
                w = wa * wb * wc * det
                Kcc += np.einsum("e,eai,ab,ebj->eij", w, B, D, B, optimize=True)
                Ms += (rho * w)[:, None, None] * np.outer(N, N)[None]
                if incompatible:
                    dP = np.diag([-2 * a, -2 * b, -2 * c])          # d P_m / d xi_a
                    dPx = np.linalg.solve(J0, np.broadcast_to(dP, (ne, 3, 3)))
                    dPx *= (det0 / det)[:, None, None]
                    Bi = _B(dPx)
                    Kci += np.einsum("e,eai,ab,ebj->eij", w, B, D, Bi, optimize=True)
                    Kii += np.einsum("e,eai,ab,ebj->eij", w, Bi, D, Bi, optimize=True)
    K = Kcc
    if incompatible:
        K = Kcc - np.einsum("eij,ejk->eik", Kci, np.linalg.solve(Kii, np.swapaxes(Kci, 1, 2)))
    if lumped:
        d = Ms.sum(axis=2)
        Ms = np.zeros_like(Ms)
        Ms[:, np.arange(8), np.arange(8)] = d
    return K, _expand_mass(Ms)


# ----------------------------------------------------------------- meshing helper
def box_mesh(lx: float, ly: float, lz: float, nx: int, ny: int, nz: int,
             kind: str = "hex8", origin=(0.0, 0.0, 0.0)):
    """Structured mesh of the box ``[0,lx]x[0,ly]x[0,lz]`` (+ origin).

    ``kind`` = ``hex8``, ``tet4`` (6 tets per brick, conforming) or
    ``tet10`` (``tet4`` + mid-edge nodes).  Returns ``nodes``, ``conn``.
    """
    xs = np.linspace(0, lx, nx + 1) + origin[0]
    ys = np.linspace(0, ly, ny + 1) + origin[1]
    zs = np.linspace(0, lz, nz + 1) + origin[2]
    Z, Y, Xg = np.meshgrid(zs, ys, xs, indexing="ij")
    nodes = np.stack([Xg.ravel(), Y.ravel(), Z.ravel()], 1)
    nid = lambda i, j, k: (k * (ny + 1) + j) * (nx + 1) + i   # noqa: E731
    I, J, Kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    I, J, Kk = I.ravel(), J.ravel(), Kk.ravel()
    hexes = np.stack([nid(I, J, Kk), nid(I + 1, J, Kk), nid(I + 1, J + 1, Kk), nid(I, J + 1, Kk),
                      nid(I, J, Kk + 1), nid(I + 1, J, Kk + 1), nid(I + 1, J + 1, Kk + 1),
                      nid(I, J + 1, Kk + 1)], 1)
    if kind == "hex8":
        return nodes, hexes
    # 6 tetrahedra around the 0-6 diagonal (conforming for a structured grid)
    split = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))
    tets = np.concatenate([hexes[:, list(s)] for s in split], 0)
    vol = tet_volumes(nodes[tets])
    neg = vol < 0
    tets[neg, 1], tets[neg, 2] = tets[neg, 2].copy(), tets[neg, 1].copy()
    if kind == "tet4":
        return nodes, tets
    if kind != "tet10":
        raise ValueError(kind)
    return tet4_to_tet10(nodes, tets)


def tet4_to_tet10(nodes: np.ndarray, tets: np.ndarray):
    """Add straight mid-edge nodes to a tet4 mesh -> (nodes, tet10 conn)."""
    nodes = np.asarray(nodes, float)
    tets = np.asarray(tets, int)
    edges = np.concatenate([tets[:, [i, j]] for i, j in TET10_EDGES], 0)
    es = np.sort(edges, axis=1)
    uniq, inv = np.unique(es, axis=0, return_inverse=True)
    mid = 0.5 * (nodes[uniq[:, 0]] + nodes[uniq[:, 1]])
    mid_id = nodes.shape[0] + inv.ravel().reshape(len(TET10_EDGES), -1).T
    return np.vstack([nodes, mid]), np.hstack([tets, mid_id])


# ----------------------------------------------------------------- rotations
def rotations_from_translations(nodes: np.ndarray, conns, translations: np.ndarray,
                                nodes_mask=None) -> np.ndarray:
    """Nodal rotations from a translation field by local least squares.

    For each node ``i`` the displacement gradient ``G`` is fitted to the
    neighbours ``j`` (nodes sharing an element):
    ``u_j - u_i ~= G (x_j - x_i)``, and the infinitesimal rotation is the
    skew part ``theta = (G32 - G23, G13 - G31, G21 - G12) / 2 = curl(u)/2``.
    Exact for rigid-body motions and linear fields.  Neighbour sets that do
    not span 3-D (e.g. planar shells) use the pseudo-inverse (in-plane
    gradient only).

    ``conns`` list of (ne, nn) connectivity arrays, ``translations``
    (Nn, 3, nm).  Returns (Nn, 3, nm) rotations (zero where
    ``nodes_mask`` is False).
    """
    nodes = np.asarray(nodes, float)
    nn = nodes.shape[0]
    pairs = []
    for c in conns:
        c = np.asarray(c, int)
        k = c.shape[1]
        a, b = np.meshgrid(np.arange(k), np.arange(k), indexing="ij")
        m = a != b
        pairs.append(np.stack([c[:, a[m]].ravel(), c[:, b[m]].ravel()], 1))
    P = np.unique(np.concatenate(pairs, 0), axis=0)
    i, j = P[:, 0], P[:, 1]
    dx = nodes[j] - nodes[i]
    A = np.zeros((nn, 3, 3))
    np.add.at(A, i, np.einsum("pa,pb->pab", dx, dx))
    Ai = np.linalg.pinv(A, rcond=1e-8, hermitian=True)
    cvec = np.einsum("pab,pb->pa", Ai[i], dx)             # (np,3) weights
    U = np.asarray(translations, float)
    du = U[j] - U[i]                                      # (np,3,nm)
    G = np.zeros((nn, 3, 3) + U.shape[2:])
    np.add.at(G, i, np.einsum("pa...,pb->pab...", du, cvec))
    rot = 0.5 * np.stack([G[:, 2, 1] - G[:, 1, 2], G[:, 0, 2] - G[:, 2, 0],
                          G[:, 1, 0] - G[:, 0, 1]], 1)
    if nodes_mask is not None:
        rot[~np.asarray(nodes_mask, bool)] = 0.0
    return rot

