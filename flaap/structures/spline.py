"""Structure-to-aerodynamics interpolation (splines).

Port / replacement of

* ``Thin_plate_spline.f90`` (Burkardt RBF library: ``phi1`` multiquadric,
  ``phi2`` inverse multiquadric, ``phi3`` thin-plate, ``phi4`` gaussian,
  ``rbf_weight`` / ``rbf_interp_nd`` with SVD solve) -> :func:`rbf_matrix`
* the inverse-distance 4-corner interpolation in ``calculation.f90``
  (``Zcp``/``Zxcp``/``Zf``) -> :func:`legacy_corner_matrices`
* new: infinite-plate spline (Harder & Desmarais 1972, Nastran SPLINE1)
  with analytic slopes, and a beam spline (Nastran SPLINE2-like) for
  stick models (Goland wing).

All splines return linear operators acting on the full structural DOF
vector ``u`` (6 per node: ux, uy, uz, rx, ry, rz)::

    h(points)    = G  @ u      (normal displacement)
    dh/dx(points)= Gx @ u      (streamwise slope)

so that for a modal model ``h_modes = G @ shapes.reshape(-1, nm)``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------- RBF kernels
def phi1(r, r0):
    """Multiquadric sqrt(r^2 + r0^2)."""
    return np.sqrt(r * r + r0 * r0)


def phi2(r, r0):
    """Inverse multiquadric."""
    return 1.0 / np.sqrt(r * r + r0 * r0)


def phi3(r, r0):
    """Thin-plate spline r^2 log(r / r0) (0 at r = 0)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        v = r * r * np.log(r / r0)
    return np.where(r > 0, v, 0.0)


def phi4(r, r0):
    """Gaussian exp(-r^2 / (2 r0^2))."""
    return np.exp(-0.5 * r * r / (r0 * r0))


_PHI = {"multiquadric": phi1, "inverse_multiquadric": phi2, "tps": phi3, "gaussian": phi4,
        1: phi1, 2: phi2, 3: phi3, 4: phi4}


def rbf_matrix(xd: np.ndarray, xi: np.ndarray, kind="tps", r0: float | None = None,
               polynomial: bool = True, rcond: float = 1e-12) -> np.ndarray:
    """Interpolation operator ``fi = G @ fd`` (rbf_weight + rbf_interp_nd).

    ``xd`` (nd, m) data sites, ``xi`` (ni, m) evaluation sites.  With
    ``polynomial=True`` a linear polynomial is appended (exact for rigid
    body motion), which the Fortran version lacked.
    """
    xd = np.atleast_2d(xd)
    xi = np.atleast_2d(xi)
    phi = _PHI[kind]
    if r0 is None:
        ext = xd.max(axis=0) - xd.min(axis=0)
        r0 = float(np.prod(ext[ext > 0]) / xd.shape[0]) ** (1.0 / max(1, np.sum(ext > 0)))
    rdd = np.linalg.norm(xd[:, None, :] - xd[None, :, :], axis=-1)
    rid = np.linalg.norm(xi[:, None, :] - xd[None, :, :], axis=-1)
    A = phi(rdd, r0)
    B = phi(rid, r0)
    nd = xd.shape[0]
    if polynomial:
        P = np.hstack([np.ones((nd, 1)), xd])
        Pi = np.hstack([np.ones((xi.shape[0], 1)), xi])
        npoly = P.shape[1]
        big = np.zeros((nd + npoly, nd + npoly))
        big[:nd, :nd] = A
        big[:nd, nd:] = P
        big[nd:, :nd] = P.T
        inv = np.linalg.pinv(big, rcond=rcond)
        return np.hstack([B, Pi]) @ inv[:, :nd]
    return B @ np.linalg.pinv(A, rcond=rcond)


# ----------------------------------------------------------- infinite plate
def _ips_kernel(r2):
    with np.errstate(divide="ignore", invalid="ignore"):
        v = r2 * np.log(r2)
    return np.where(r2 > 0, v, 0.0)


def ips_matrices(xd: np.ndarray, xi: np.ndarray, rcond: float = 1e-13):
    """Infinite-plate spline on the (x, y) plane: returns (G, Gx)
    acting on nodal values (nd,), giving value and d/dx at ``xi``."""
    xd = np.asarray(xd, float)[:, :2]
    xi = np.asarray(xi, float)[:, :2]
    scale = max(np.ptp(xd[:, 0]), np.ptp(xd[:, 1]), 1e-12)
    xs, xis = xd / scale, xi / scale
    nd = xs.shape[0]
    d = xs[:, None, :] - xs[None, :, :]
    A = _ips_kernel(np.sum(d * d, -1))
    P = np.hstack([np.ones((nd, 1)), xs])
    big = np.zeros((nd + 3, nd + 3))
    big[:nd, :nd] = A
    big[:nd, nd:] = P
    big[nd:, :nd] = P.T
    inv = np.linalg.pinv(big, rcond=rcond)[:, :nd]
    di = xis[:, None, :] - xs[None, :, :]
    r2 = np.sum(di * di, -1)
    B = np.hstack([_ips_kernel(r2), np.ones((xi.shape[0], 1)), xis])
    with np.errstate(divide="ignore", invalid="ignore"):
        dk = np.where(r2 > 0, 2 * di[..., 0] * (np.log(r2) + 1.0), 0.0)
    Bx = np.hstack([dk, np.zeros((xi.shape[0], 1)), np.ones((xi.shape[0], 1)),
                    np.zeros((xi.shape[0], 1))]) / scale
    return B @ inv, Bx @ inv


# ------------------------------------------------------------ spline object
@dataclass
class SplineOperator:
    """Linear maps from the structural DOF vector to aero box quantities."""

    h_colloc: np.ndarray      # (N, 6Nn) normal displacement at 3/4-chord points
    hx_colloc: np.ndarray     # (N, 6Nn) streamwise slope at 3/4-chord points
    h_force: np.ndarray       # (N, 6Nn) normal displacement at 1/4-chord (force) points
    force_to_struct: np.ndarray | None = None  # (6Nn, 3N) transfer of box forces (optional)

    def modal(self, shapes: np.ndarray):
        """Return (h_c, hx_c, h_f), each (N, nm), for mode shapes (Nn,6,nm)."""
        U = shapes.reshape(-1, shapes.shape[-1])
        return self.h_colloc @ U, self.hx_colloc @ U, self.h_force @ U


def _expand(Gnode: np.ndarray, nn: int, comp: int) -> np.ndarray:
    G = np.zeros((Gnode.shape[0], 6 * nn))
    G[:, comp::6] = Gnode
    return G


def surface_spline(mesh, nodes: np.ndarray, method: str = "ips", box_mask=None,
                   node_mask=None, eps_x: float | None = None) -> SplineOperator:
    """Plate-type spline (IPS or RBF) from nodal translations to the boxes.

    The normal displacement is ``h = n_y u_y + n_z u_z``; the slope is taken
    analytically (IPS) or by central differences (RBF).
    """
    nodes = np.asarray(nodes, float)
    nn = nodes.shape[0]
    bm = np.ones(mesh.n, bool) if box_mask is None else np.asarray(box_mask)
    nm = np.ones(nn, bool) if node_mask is None else np.asarray(node_mask)
    xd = nodes[nm]
    Hc = np.zeros((mesh.n, 6 * nn))
    Hx = np.zeros((mesh.n, 6 * nn))
    Hf = np.zeros((mesh.n, 6 * nn))
    # planar coordinates: project onto each surface's (x, s) plane where s is the
    # in-plane spanwise coordinate (handles dihedral/vertical surfaces)
    gam = np.median(mesh.gamma[bm])
    sproj = lambda p: np.stack([p[:, 0], p[:, 1] * np.cos(gam) + p[:, 2] * np.sin(gam)], 1)  # noqa
    xs = sproj(xd)
    pc, pf = sproj(mesh.colloc[bm]), sproj(mesh.force_point[bm])
    if method == "ips":
        Gc, Gcx = ips_matrices(xs, pc)
        Gf, _ = ips_matrices(xs, pf)
    else:
        kind = "tps" if method == "rbf" else method
        Gc = rbf_matrix(xs, pc, kind)
        Gf = rbf_matrix(xs, pf, kind)
        h = eps_x or 1e-4 * np.ptp(xs[:, 0])
        Gcx = (rbf_matrix(xs, pc + [h, 0], kind) - rbf_matrix(xs, pc - [h, 0], kind)) / (2 * h)
    idx = np.where(nm)[0]
    ny = mesh.normal[bm, 1][:, None]
    nz = mesh.normal[bm, 2][:, None]
    rows = np.where(bm)[0]
    for comp, nvec in ((1, ny), (2, nz)):
        cols = 6 * idx + comp
        Hc[np.ix_(rows, cols)] += nvec * Gc
        Hx[np.ix_(rows, cols)] += nvec * Gcx
        Hf[np.ix_(rows, cols)] += nvec * Gf
    return SplineOperator(Hc, Hx, Hf)


def beam_spline(mesh, nodes: np.ndarray, box_mask=None, axis: int = 1,
                hermite: bool = True) -> SplineOperator:
    """Beam (stick-model) spline: boxes are rigidly attached to the beam
    line, ``h(x, y) = w(y) - (x - x_ea(y)) * ry(y)``.

    ``w`` is interpolated with cubic Hermite polynomials using the bending
    slope ``rx = dw/dy`` (``hermite=True``) or linearly, ``ry`` (torsion,
    nose-up positive) linearly.  ``nodes`` must be ordered along the span.
    """
    nodes = np.asarray(nodes, float)
    nn = nodes.shape[0]
    bm = np.ones(mesh.n, bool) if box_mask is None else np.asarray(box_mask)
    # the interpolation walks the beam from small to large s; a mirrored / rotated model may list the
    # nodes the other way round -> work on the sorted line and map the columns back
    order = np.argsort(nodes[:, axis], kind="stable")
    col = np.empty(nn, int)
    col[order] = np.arange(nn)
    nodes = nodes[order]
    s_n = nodes[:, axis]
    Hc = np.zeros((mesh.n, 6 * nn))
    Hx = np.zeros((mesh.n, 6 * nn))
    Hf = np.zeros((mesh.n, 6 * nn))

    def fill(H, pts, rows, slope=False):
        for r, p in zip(rows, pts):
            s = p[axis]
            j = int(np.clip(np.searchsorted(s_n, s) - 1, 0, nn - 2))
            L = s_n[j + 1] - s_n[j]
            t = (s - s_n[j]) / L
            t = min(max(t, 0.0), 1.0) if s <= s_n[-1] else t
            xea = nodes[j, 0] * (1 - t) + nodes[j + 1, 0] * t
            dx = p[0] - xea
            if slope:
                H[r, 6 * j + 4] += -(1 - t)
                H[r, 6 * (j + 1) + 4] += -t
                continue
            if hermite:
                h00 = 2 * t ** 3 - 3 * t ** 2 + 1
                h10 = t ** 3 - 2 * t ** 2 + t
                h01 = -2 * t ** 3 + 3 * t ** 2
                h11 = t ** 3 - t ** 2
                H[r, 6 * j + 2] += h00
                H[r, 6 * j + 3] += h10 * L
                H[r, 6 * (j + 1) + 2] += h01
                H[r, 6 * (j + 1) + 3] += h11 * L
            else:
                H[r, 6 * j + 2] += 1 - t
                H[r, 6 * (j + 1) + 2] += t
            H[r, 6 * j + 4] += -dx * (1 - t)
            H[r, 6 * (j + 1) + 4] += -dx * t

    rows = np.where(bm)[0]
    fill(Hc, mesh.colloc[bm], rows)
    fill(Hf, mesh.force_point[bm], rows)
    fill(Hx, mesh.colloc[bm], rows, slope=True)
    if not np.array_equal(order, np.arange(nn)):        # sorted-node columns -> input-node columns
        cols = (6 * order[:, None] + np.arange(6)[None, :]).ravel()
        Hc, Hx, Hf = Hc[:, np.argsort(cols)], Hx[:, np.argsort(cols)], Hf[:, np.argsort(cols)]
    return SplineOperator(Hc, Hx, Hf)


def legacy_corner_matrices(mesh, nodes: np.ndarray) -> SplineOperator:
    """Fortran ``calculation.f90`` scheme: inverse-distance weighting of the
    four lattice-corner nodes for displacement, and the finite-difference
    chordwise slope ``-((zur+zul)/2-(zlr+zll)/2)/chord``.

    Kept for one-to-one regression against the original program.  Nodes
    must coincide with box corners (legacy mode-shape grids do); corners
    without a node (clamped root row) are taken as zero displacement.
    """
    nodes = np.asarray(nodes, float)
    nn = nodes.shape[0]
    Hc = np.zeros((mesh.n, 6 * nn))
    Hx = np.zeros((mesh.n, 6 * nn))
    Hf = np.zeros((mesh.n, 6 * nn))
    tol = 1e-5 * max(1.0, np.ptp(nodes[:, :2]))

    def find(p):
        d = np.hypot(nodes[:, 0] - p[0], nodes[:, 1] - p[1])
        j = int(np.argmin(d))
        return j if d[j] < tol else None

    for b in range(mesh.n):
        c = mesh.corners[b]
        ids = [find(c[i]) for i in range(4)]
        for H, p in ((Hc, mesh.colloc[b]), (Hf, mesh.force_point[b])):
            r = np.hypot(c[:, 0] - p[0], c[:, 1] - p[1])
            w = (1.0 / r) / np.sum(1.0 / r)
            for i, j in enumerate(ids):
                if j is not None:
                    H[b, 6 * j + 2] += w[i]
        ch = mesh.chord[b]
        # LE corners 0,1 ; TE corners 2,3
        for i, j in enumerate(ids):
            if j is not None:
                Hx[b, 6 * j + 2] += (0.5 if i in (2, 3) else -0.5) / ch
    return SplineOperator(Hc, Hx, Hf)
