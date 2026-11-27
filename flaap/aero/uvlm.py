"""Vortex-lattice aerodynamics in the time domain.

Three models share one vortex-ring discretisation of the aerodynamic boxes
of an :class:`~flaap.geometry.AeroMesh` (box ``b`` <-> ring ``b``):

* :func:`steady_vlm` - classical steady VLM (horseshoes, :mod:`flaap.aero.vlm`)
  with lift, pitching moment, span load and the Trefftz-plane induced drag.
* :class:`UVLM` - the nonlinear unsteady vortex-lattice method of Katz &
  Plotkin (2001, ch. 13): vortex rings with the leading segment on the box
  1/4 chord and the collocation point at 3/4 chord, explicit wake shedding
  from the trailing edge, prescribed (``U t``) or free (force-free, rolled-up)
  wake.  Prescribed rigid kinematics (incidence step, harmonic plunge /
  pitch) and frozen gusts.  Loads follow Simpson, Palacios & Murua (2013,
  *AIAA J.* 51(7)): Kutta-Joukowski on every bound segment with the local
  velocity (free stream + body motion + gust + induced) plus the unsteady
  term ``rho dGamma/dt A n``; this gives lift, induced drag / thrust
  (leading-edge suction included), side force and pitching moment versus
  time.
* :class:`LinearUVLM` - the UVLM linearised about a flat, frozen wake
  (wake panels ``U dt`` long, Murua, Palacios & Graham 2012).  It gives
  (i) frequency-domain GAF tables compatible with :class:`~flaap.aero.gaf.GAFTable`
  (PK flutter, gust, RFA ...), and (ii) the monolithic discrete-time
  aeroelastic state space ``E x[n+1] = F x[n]`` (UVLM + Newmark structure)
  for flutter eigenvalues and time responses.

Axes: ``x`` downstream, ``y`` right, ``z`` up; the air flows along ``+x``
relative to the wing, so lift is ``F_z`` and drag ``F_x``.  Pitch is nose-up
positive, plunge ``h`` is upward positive.  Incompressible (``M = 0``).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from ..geometry import AeroMesh

_FOUR_PI = 4.0 * np.pi


# ======================================================================
# Biot-Savart kernels
# ======================================================================
def _seg_kernel(P, A, B, core):
    """(n,m,3) velocity at P (n,3) from unit segments A->B (m,3).

    ``core`` is a cut-off radius (Rankine-like regularisation
    ``d^2 -> d^2 + core^2``); points on the segment line give zero.
    """
    r1 = P[:, None, :] - A[None, :, :]
    r2 = P[:, None, :] - B[None, :, :]
    r0 = (B - A)[None, :, :]
    cr = np.cross(r1, r2)
    cr2 = np.einsum("nmk,nmk->nm", cr, cr)
    l2 = np.einsum("nmk,nmk->nm", r0, r0)
    n1 = np.maximum(np.sqrt(np.einsum("nmk,nmk->nm", r1, r1)), 1e-300)
    n2 = np.maximum(np.sqrt(np.einsum("nmk,nmk->nm", r2, r2)), 1e-300)
    num = np.einsum("nmk,nmk->nm", r0, r1) / n1 - np.einsum("nmk,nmk->nm", r0, r2) / n2
    f = num / (_FOUR_PI * (cr2 + core * core * l2) + 1e-300)
    return cr * f[..., None]


def _mirror(p):
    q = np.array(p, float, copy=True)
    q[..., 1] *= -1.0
    return q


def _with_images(A, B, G, symmetry):
    """Append the mirror segments (y -> -y, reversed direction)."""
    if symmetry not in ("symmetric", "antisymmetric"):
        return A, B, G
    s = 1.0 if symmetry == "symmetric" else -1.0
    return (np.concatenate([A, _mirror(B)]), np.concatenate([B, _mirror(A)]),
            np.concatenate([G, s * G]))


def segment_velocity(P, A, B, G, core=1e-9, symmetry="none", chunk=1_500_000):
    """Velocity (n,3) at points P from segments A->B with strengths G."""
    P = np.atleast_2d(np.asarray(P, float))
    A, B, G = _with_images(np.asarray(A, float), np.asarray(B, float), np.asarray(G, float), symmetry)
    out = np.zeros_like(P)
    if A.shape[0] == 0 or P.shape[0] == 0:
        return out
    step = max(1, chunk // max(A.shape[0], 1))
    for i in range(0, P.shape[0], step):
        out[i:i + step] = np.einsum("nmk,m->nk", _seg_kernel(P[i:i + step], A, B, core), G)
    return out


def ring_normal_matrix(P, nrm, rings, core=1e-9, symmetry="none", chunk=1_500_000):
    """Normal velocity (n, Nr) at P (n,3) along nrm (n,3) per unit ring strength.

    ``rings`` is (Nr,4,3) with corners ordered so that a positive strength
    produces lift (inboard-LE, outboard-LE, outboard-TE, inboard-TE).
    """
    P = np.atleast_2d(P)
    Nr = rings.shape[0]
    out = np.zeros((P.shape[0], Nr))
    if Nr == 0:
        return out
    sgn = 1.0 if symmetry == "symmetric" else -1.0
    step = max(1, chunk // (4 * Nr))
    for i in range(0, P.shape[0], step):
        Pi, ni = P[i:i + step], nrm[i:i + step]
        acc = np.zeros((Pi.shape[0], Nr))
        for a in range(4):
            A, B = rings[:, a], rings[:, (a + 1) % 4]
            acc += np.einsum("nmk,nk->nm", _seg_kernel(Pi, A, B, core), ni)
            if symmetry in ("symmetric", "antisymmetric"):
                acc += sgn * np.einsum("nmk,nk->nm", _seg_kernel(Pi, _mirror(B), _mirror(A), core), ni)
        out[i:i + step] = acc
    return out


# ======================================================================
# Lattice geometry
# ======================================================================
def corner_grids(mesh: AeroMesh):
    """Per surface: (panel corner grid (nc+1, ns+1, 3), box index (nc, ns)).

    Rows run chordwise (LE -> TE), columns spanwise (inboard -> outboard).
    """
    out = []
    for sid in np.unique(mesh.surface_id):
        idx = np.where(mesh.surface_id == sid)[0]
        strips = list(np.unique(mesh.strip_id[idx]))
        nc = int(mesh.chord_index[idx].max()) + 1
        ns = len(strips)
        G = np.zeros((nc + 1, ns + 1, 3))
        box = np.full((nc, ns), -1, int)
        for b in idx:
            j = strips.index(mesh.strip_id[b])
            i = int(mesh.chord_index[b])
            box[i, j] = b
            c = mesh.corners[b]
            G[i, j], G[i, j + 1], G[i + 1, j + 1], G[i + 1, j] = c[0], c[1], c[2], c[3]
        if (box < 0).any():
            raise ValueError("UVLM needs a structured (nchord x nspan) box grid per surface")
        out.append((G, box))
    return out


def ring_grid(G: np.ndarray) -> np.ndarray:
    """Vortex-ring nodes: panel grid shifted by 1/4 panel chord downstream."""
    R = np.empty_like(G)
    R[:-1] = G[:-1] + 0.25 * (G[1:] - G[:-1])
    R[-1] = G[-1] + 0.25 * (G[-1] - G[-2])
    return R


def grid_rings(R: np.ndarray) -> np.ndarray:
    """(nr*ns, 4, 3) ring corners from a node grid (nr+1, ns+1, 3), row-major."""
    c = np.stack([R[:-1, :-1], R[:-1, 1:], R[1:, 1:], R[1:, :-1]], axis=2)
    return c.reshape(-1, 4, 3)


def grid_segments(R: np.ndarray, Gam: np.ndarray):
    """Unique segments (A, B, strength) of a ring grid with strengths Gam (nr, ns)."""
    nr, ns = Gam.shape
    z_row = np.zeros((1, ns))
    s_span = np.vstack([Gam, z_row]) - np.vstack([z_row, Gam])          # (nr+1, ns)
    z_col = np.zeros((nr, 1))
    s_chord = np.hstack([z_col, Gam]) - np.hstack([Gam, z_col])         # (nr, ns+1)
    A = np.concatenate([R[:, :-1].reshape(-1, 3), R[:-1, :].reshape(-1, 3)])
    B = np.concatenate([R[:, 1:].reshape(-1, 3), R[1:, :].reshape(-1, 3)])
    return A, B, np.concatenate([s_span.ravel(), s_chord.ravel()])


def _quad_normals(c):
    n = np.cross(c[:, 2] - c[:, 0], c[:, 1] - c[:, 3])
    return n / np.linalg.norm(n, axis=1)[:, None]


def _quad_area(c):
    return 0.5 * np.linalg.norm(np.cross(c[:, 2] - c[:, 0], c[:, 1] - c[:, 3]), axis=1)


@dataclass
class Lattice:
    """Vortex-ring lattice of an AeroMesh (arrays indexed like the mesh boxes)."""

    mesh: AeroMesh
    grids: list            # per surface: dict(G, R, box)
    rings: np.ndarray      # (N,4,3)
    colloc: np.ndarray     # (N,3)
    normal: np.ndarray     # (N,3)
    area: np.ndarray       # (N,) panel area
    ring_area: np.ndarray  # (N,)
    ring_center: np.ndarray
    width: np.ndarray      # (N,) span width of the ring leading segment (y-z projection)
    prev: np.ndarray       # (N,) index of the ring ahead (-1 at the LE)
    te: np.ndarray         # (Nstrip,) TE ring index per strip (surface order)

    @classmethod
    def from_mesh(cls, mesh: AeroMesh) -> "Lattice":
        grids, N = [], mesh.n
        rings = np.zeros((N, 4, 3))
        prev = np.full(N, -1, int)
        te = []
        for G, box in corner_grids(mesh):
            R = ring_grid(G)
            rings[box.ravel()] = grid_rings(R)
            prev[box[1:].ravel()] = box[:-1].ravel()
            te.extend(box[-1].tolist())
            grids.append({"G": G, "R": R, "box": box})
        pc = mesh.corners
        colloc = mesh.colloc.copy()
        normal = _quad_normals(pc)
        lead = rings[:, 1] - rings[:, 0]
        return cls(mesh, grids, rings, colloc, normal, _quad_area(pc), _quad_area(rings),
                   rings.mean(axis=1), np.hypot(lead[:, 1], lead[:, 2]), prev, np.array(te, int))

    @property
    def symmetry(self):
        return self.mesh.symmetry

    @property
    def s_ref(self) -> float:
        return float(self.area.sum())

    def strip_y(self):
        """Mid-span y of every strip (TE ring order) and its width."""
        r = self.rings[self.te]
        return 0.5 * (r[:, 0, 1] + r[:, 1, 1]), self.width[self.te]


# ======================================================================
# Steady VLM loads
# ======================================================================
@dataclass
class SteadyVLMResult:
    alpha_deg: np.ndarray
    CL: np.ndarray
    CDi: np.ndarray
    Cm: np.ndarray
    CL_alpha: float           # per rad
    e: float                  # span efficiency CL^2 / (pi A CDi)
    aspect_ratio: float
    s_ref: float
    c_ref: float
    x_ref: float
    y_strip: np.ndarray
    cl_c: np.ndarray          # (nalpha, nstrip) local lift coefficient x chord / c_ref
    cl_local: np.ndarray      # (nalpha, nstrip) local section lift coefficient
    x_ac: float               # aerodynamic centre (x of the CL_alpha line of action)


def steady_vlm(mesh: AeroMesh, alpha_deg, mach: float = 0.0, x_ref: float | None = None,
               c_ref: float | None = None) -> SteadyVLMResult:
    """Steady VLM (horseshoes on the box 1/4 chord, Prandtl-Glauert for M > 0).

    The induced drag is evaluated in the Trefftz plane,
    ``D_i = -rho/2 sum Gamma w_T dy`` with ``w_T`` the normalwash of the
    trailing legs far downstream.  Coefficients refer to the modelled
    area (the half wing for symmetric meshes, i.e. the full-wing values).
    """
    from .vlm import steady_aic
    alpha = np.atleast_1d(np.asarray(alpha_deg, float))
    D0 = steady_aic(mesh, mach)
    nz = mesh.normal[:, 2]
    W = -np.outer(nz, np.radians(alpha))                 # normalwash for incidence alpha (w/U)
    dcp = np.linalg.solve(D0, W)                          # (N, na)
    S = mesh.total_area
    strips = [idx for _, idx in mesh.strips()]
    width = np.array([2 * mesh.e[idx[0]] for idx in strips])
    span_half = float(width.sum())
    full = mesh.symmetry in ("symmetric", "antisymmetric")
    A_ratio = (2 * span_half) ** 2 / (2 * S) if full else span_half ** 2 / S
    c_ref = S / span_half if c_ref is None else c_ref
    x_ref = float(mesh.corners[:, 0, 0].min()) if x_ref is None else x_ref
    F = dcp * (mesh.area * nz)[:, None]                   # lift / q per box
    CL = F.sum(0) / S
    Cm = -((mesh.force_point[:, 0] - x_ref)[:, None] * F).sum(0) / (S * c_ref)
    # Trefftz plane: Gamma/U per strip, piecewise-linear sheet (see trefftz_drag)
    gam_box = 0.5 * mesh.chord[:, None] * dcp
    gam = np.array([gam_box[idx].sum(0) for idx in strips])            # (ns, na)
    edges = [(mesh.A[idx[0], 1:], mesh.B[idx[0], 1:], int(mesh.surface_id[idx[0]])) for idx in strips]
    CDi = np.array([trefftz_drag(edges, gam[:, j], mesh.symmetry) for j in range(alpha.size)]) / S
    y_strip = np.array([mesh.force_point[idx, 1].mean() for idx in strips])
    chord_strip = np.array([mesh.chord[idx].sum() for idx in strips])
    lstrip = np.array([F[idx].sum(0) for idx in strips]) / width[:, None]     # cl * c
    a = np.radians(alpha)
    if alpha.size > 1:
        cla = float(np.polyfit(a, CL, 1)[0])
        cma = float(np.polyfit(a, Cm, 1)[0])
    else:
        cla = float(CL[0] / a[0]) if a[0] != 0 else float("nan")
        cma = float(Cm[0] / a[0]) if a[0] != 0 else float("nan")
    x_ac = x_ref - cma / cla * c_ref if np.isfinite(cla) and cla != 0 else float("nan")
    with np.errstate(divide="ignore", invalid="ignore"):
        e = CL ** 2 / (np.pi * A_ratio * CDi)
    e_val = float(np.nanmedian(e[np.abs(alpha) > 1e-9])) if np.any(np.abs(alpha) > 1e-9) else float("nan")
    return SteadyVLMResult(alpha, CL, CDi, Cm, cla, e_val, A_ratio, S, c_ref, x_ref, y_strip,
                           (lstrip / c_ref).T, (lstrip / chord_strip[:, None]).T, x_ac)


def trefftz_drag(edges, gam, symmetry="none", refine: int = 24) -> float:
    """Induced drag / (q) in the Trefftz plane for strip circulations ``gam`` (Gamma / U).

    ``edges`` = [(inboard (y, z), outboard (y, z), surface id)] per strip.
    The strip values are joined by a piecewise-linear ``Gamma(s)`` along each
    surface trace (``Gamma -> 0`` at free ends, zero slope on the symmetry
    plane), split into ``refine`` sub-panels per strip and integrated with
    discrete trailing vortices, ``D_i / q = -sum Gamma w_T ds``.  This removes
    the O(1/n) bias of evaluating a step-wise loading at the strip midpoints.
    """
    gam = np.asarray(gam, float)
    sym = symmetry in ("symmetric", "antisymmetric")
    vort_p, vort_g, pts, pg, pds, pn = [], [], [], [], [], []
    sids = [e[2] for e in edges]
    for sid in dict.fromkeys(sids):
        idx = [i for i, s_ in enumerate(sids) if s_ == sid]
        P0 = np.array([edges[i][0] for i in idx])
        P1 = np.array([edges[i][1] for i in idx])
        mid = 0.5 * (P0 + P1)
        nodes = np.vstack([P0[:1], mid, P1[-1:]])
        on_plane0 = sym and abs(P0[0, 0]) < 1e-6 * max(1.0, abs(P1[-1, 0]))
        g = np.r_[gam[idx[0]] if on_plane0 else 0.0, gam[idx], 0.0]
        seg = np.diff(nodes, axis=0)
        L = np.linalg.norm(seg, axis=1)
        s_nodes = np.r_[0.0, np.cumsum(L)]
        s_f = np.unique(np.concatenate([np.linspace(s_nodes[i], s_nodes[i + 1], refine // 2 + 1)
                                        for i in range(len(L))]))
        P_f = np.column_stack([np.interp(s_f, s_nodes, nodes[:, 0]), np.interp(s_f, s_nodes, nodes[:, 1])])
        G_f = np.interp(s_f, s_nodes, g)
        # constant Gamma_j on each fine sub-panel; trailing vortex (along +x) at
        # node k carries Gamma_inboard - Gamma_outboard.  The symmetry-plane
        # vortex is cancelled (symmetric) or doubled (antisymmetric) by its image.
        Gj = 0.5 * (G_f[1:] + G_f[:-1])
        vort_p.append(P_f)
        vort_g.append(np.r_[0.0, Gj] - np.r_[Gj, 0.0])
        cm = 0.5 * (P_f[1:] + P_f[:-1])
        d = np.diff(P_f, axis=0)
        ds = np.linalg.norm(d, axis=1)
        t = d / ds[:, None]
        pts.append(cm)
        pg.append(Gj)
        pds.append(ds)
        pn.append(np.column_stack([-t[:, 1], t[:, 0]]))       # (y, z) normal: (0,1) for a flat wing
    V, Gv = np.vstack(vort_p), np.concatenate(vort_g)
    if sym:
        V = np.vstack([V, V * np.array([-1.0, 1.0])])
        Gv = np.r_[Gv, -Gv if symmetry == "symmetric" else Gv]
    C, Gc, ds, nn = np.vstack(pts), np.concatenate(pg), np.concatenate(pds), np.vstack(pn)
    dy = C[:, None, 0] - V[None, :, 0]
    dz = C[:, None, 1] - V[None, :, 1]
    r2 = dy * dy + dz * dz
    vy = (-dz / (2 * np.pi * r2)) @ Gv
    vz = (dy / (2 * np.pi * r2)) @ Gv
    wn = vy * nn[:, 0] + vz * nn[:, 1]
    return float(-(Gc * wn * ds).sum())


# ======================================================================
# Kinematics and gusts
# ======================================================================
def harmonic(amplitude: float, omega: float, phase: float = 0.0, mean: float = 0.0):
    """t -> (value, rate) of ``mean + amplitude sin(omega t + phase)``."""
    return lambda t: (mean + amplitude * np.sin(omega * t + phase), amplitude * omega * np.cos(omega * t + phase))


def sharp_edged_gust(w0: float, U: float, x_front: float):
    """Frozen sharp-edged gust whose front reaches ``x_front`` at t = 0."""
    return lambda x, t: w0 * ((U * t - (np.asarray(x) - x_front)) >= 0.0)


def one_minus_cos_gust(w0: float, H: float, U: float, x_front: float):
    """Frozen 1-cos gust of gradient H (peak w0 at penetration s = H)."""
    def g(x, t):
        s = U * t - (np.asarray(x) - x_front)
        return np.where((s >= 0) & (s <= 2 * H), 0.5 * w0 * (1 - np.cos(np.pi * s / H)), 0.0)
    return g


@dataclass
class Kinematics:
    """Prescribed rigid motion: incidence + pitch(t) about x_ref, plunge(t), gust(x, t)."""

    alpha0: float = 0.0                     # rad, applied at t = 0+ (impulsive start)
    x_ref: float = 0.0
    pitch: Callable | None = None           # t -> (theta, theta_dot)  nose up
    plunge: Callable | None = None          # t -> (h, h_dot)  up
    gust: Callable | None = None            # (x, t) -> vertical gust velocity

    def state(self, t):
        th, thd = self.pitch(t) if self.pitch else (0.0, 0.0)
        h, hd = self.plunge(t) if self.plunge else (0.0, 0.0)
        return self.alpha0 + th, thd, h, hd


def _rot_y(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


# ======================================================================
# Nonlinear UVLM
# ======================================================================
@dataclass
class UVLMResult:
    t: np.ndarray
    CL: np.ndarray
    CD: np.ndarray            # F_x / (q S): induced drag (negative = thrust)
    CY: np.ndarray
    Cm: np.ndarray            # about the (moving) reference point, nose up
    L: np.ndarray
    D: np.ndarray
    My: np.ndarray
    cl_strip: np.ndarray      # (nt, nstrip) local lift / (q c_ref width)
    y_strip: np.ndarray
    gamma: np.ndarray         # final bound circulation (N,)
    wake: list                # final wake node grids per surface
    s_ref: float
    c_ref: float
    U: float
    rho: float
    dt: float
    meta: dict = field(default_factory=dict)

    @property
    def s(self):
        """Distance travelled in semi-chords (2 U t / c_ref)."""
        return 2 * self.U * self.t / self.c_ref


class UVLM:
    """Nonlinear unsteady vortex-lattice method (rigid lifting surfaces)."""

    def __init__(self, mesh: AeroMesh, core: float | None = None, wake_core: float | None = None):
        self.lat = Lattice.from_mesh(mesh)
        self.mesh = mesh
        dx = np.mean([np.mean(g["G"][-1, :, 0] - g["G"][-2, :, 0]) for g in self.lat.grids])
        self.dx_te = float(dx)
        self.c_ref = self.lat.s_ref / float(self.lat.width[self.lat.te].sum())
        self.core = 1e-6 * self.c_ref if core is None else core
        self.wake_core = 0.05 * self.dx_te if wake_core is None else wake_core

    def default_dt(self, U: float) -> float:
        return self.dx_te / U

    def simulate(self, U: float, t_end: float, kin: Kinematics | None = None, *, rho: float = 1.225,
                 dt: float | None = None, free_wake: bool = False, wake_length: float | None = None,
                 x_ref: float | None = None, verbose: bool = False) -> UVLMResult:
        """March in time; ``wake_length`` (in reference chords) truncates the wake."""
        kin = kin or Kinematics()
        lat = self.lat
        sym = lat.symmetry
        dt = self.default_dt(U) if dt is None else float(dt)
        nsteps = int(round(t_end / dt))
        xr = kin.x_ref if x_ref is None else x_ref
        pref = np.array([xr, 0.0, 0.0])
        q = 0.5 * rho * U * U
        S, cref = lat.s_ref, self.c_ref
        max_rows = None if wake_length is None else max(2, int(np.ceil(wake_length * cref / (U * dt))))
        # body-frame geometry relative to the reference point
        col0 = lat.colloc - pref
        nrm0 = lat.normal
        cen0 = lat.ring_center - pref
        Rg0 = [g["R"] - pref for g in lat.grids]
        boxes = [g["box"] for g in lat.grids]
        # AIC is invariant under rigid motion (image symmetry preserved by pitch/plunge)
        AIC = ring_normal_matrix(lat.colloc, lat.normal, lat.rings, self.core, sym)
        from scipy.linalg import lu_factor, lu_solve
        lu = lu_factor(AIC)
        N = lat.mesh.n
        gam_prev = np.zeros(N)
        # wake per surface: nodes (rows, ns+1, 3) and strengths (rows-1, ns)
        th, thd, h, hd = kin.state(0.0)
        Rm = _rot_y(th)
        off = pref + np.array([0.0, 0.0, h])
        wakes = [(Rg[-1:] @ Rm.T + off).copy() for Rg in Rg0]
        wgam = [np.zeros((0, b.shape[1])) for b in boxes]
        te_prev = [np.zeros(b.shape[1]) for b in boxes]
        nt = nsteps + 1
        out = {k: np.zeros(nt) for k in ("CL", "CD", "CY", "Cm", "L", "D", "My")}
        nstrip = lat.te.size
        cl_strip = np.zeros((nt, nstrip))
        strip_of = np.zeros(N, int)
        k0 = 0
        for b in boxes:
            for j in range(b.shape[1]):
                strip_of[b[:, j]] = k0 + j
            k0 += b.shape[1]
        widths = lat.width[lat.te]
        ex = np.array([1.0, 0.0, 0.0])
        for n in range(1, nt):
            t = n * dt
            th, thd, h, hd = kin.state(t)
            Rm = _rot_y(th)
            off = pref + np.array([0.0, 0.0, h])
            omega = np.array([0.0, thd, 0.0])
            col = col0 @ Rm.T + off
            nrm = nrm0 @ Rm.T
            cen = cen0 @ Rm.T + off
            Rg = [R @ Rm.T + off for R in Rg0]
            # shed: prepend the current TE row; its panel carries Gamma_TE(t - dt)
            for s in range(len(Rg)):
                wakes[s] = np.concatenate([Rg[s][-1:], wakes[s]])
                wgam[s] = np.concatenate([te_prev[s][None, :], wgam[s]])
                if max_rows is not None and wgam[s].shape[0] > max_rows:
                    wgam[s] = wgam[s][:max_rows]
                    wakes[s] = wakes[s][:max_rows + 1]
            WA, WB, WG = self._wake_segments(wakes, wgam)

            def vbody(p):
                return np.array([0.0, 0.0, hd]) + np.cross(omega, p - off)

            def vfree(p):
                v = np.tile(U * ex, (p.shape[0], 1)) - vbody(p)
                if kin.gust is not None:
                    v[:, 2] += kin.gust(p[:, 0], t)
                return v

            vw = segment_velocity(col, WA, WB, WG, self.core, sym)
            rhs = -np.einsum("nk,nk->n", vfree(col) + vw, nrm)
            gam = lu_solve(lu, rhs)
            # ---- loads (Simpson et al. 2013)
            BA, BB, BG = [], [], []
            for s, b in enumerate(boxes):
                A_, B_, G_ = grid_segments(Rg[s], gam[b])
                nc, ns = b.shape
                n_span = (nc + 1) * ns
                keep = np.r_[np.arange(nc * ns), n_span + np.arange(nc * (ns + 1))]   # drop TE spanwise row
                BA.append(A_[keep]); BB.append(B_[keep]); BG.append(G_[keep])
            BA, BB, BG = np.concatenate(BA), np.concatenate(BB), np.concatenate(BG)
            mid = 0.5 * (BA + BB)
            bA, bB, bG = self._bound_all(Rg, boxes, gam)
            allA, allB, allG = np.concatenate([bA, WA]), np.concatenate([bB, WB]), np.concatenate([bG, WG])
            vmid = vfree(mid) + segment_velocity(mid, allA, allB, allG, self.core, sym)
            Fseg = rho * BG[:, None] * np.cross(vmid, BB - BA)
            dgdt = (gam - gam_prev) / dt
            Fu = rho * (dgdt * lat.ring_area)[:, None] * nrm
            Ftot = Fseg.sum(0) + Fu.sum(0)
            Mtot = np.cross(mid - off, Fseg).sum(0) + np.cross(cen - off, Fu).sum(0)
            out["L"][n], out["D"][n], out["My"][n] = Ftot[2], Ftot[0], Mtot[1]
            out["CL"][n] = Ftot[2] / (q * S)
            out["CD"][n] = Ftot[0] / (q * S)
            out["CY"][n] = Ftot[1] / (q * S)
            out["Cm"][n] = Mtot[1] / (q * S * cref)
            # span load: unsteady part per ring + KJ lift of the leading segments
            lift_ring = Fu[:, 2].copy()
            lift_ring += self._kj_lift_per_ring(Rg, boxes, gam, vfree, allA, allB, allG, rho, sym)
            cl_strip[n] = np.bincount(strip_of, lift_ring, nstrip) / (q * cref * widths)
            # ---- wake convection
            for s in range(len(wakes)):
                pts = wakes[s].reshape(-1, 3)
                if free_wake:              # air velocity (no body motion): free stream + gust + induced
                    v = np.tile(U * ex, (pts.shape[0], 1))
                    if kin.gust is not None:
                        v[:, 2] += kin.gust(pts[:, 0], t)
                    v += segment_velocity(pts, allA, allB, allG, self.wake_core, sym)
                else:
                    v = np.tile(U * ex, (pts.shape[0], 1))
                wakes[s] = (pts + v * dt).reshape(wakes[s].shape)
            te_prev = [gam[b[-1]] for b in boxes]
            gam_prev = gam
            if verbose and n % max(1, nt // 10) == 0:
                print(f"  UVLM step {n}/{nsteps}  CL = {out['CL'][n]:.4f}")
        t = np.arange(nt) * dt
        y_strip, _ = lat.strip_y()
        return UVLMResult(t[1:], out["CL"][1:], out["CD"][1:], out["CY"][1:], out["Cm"][1:], out["L"][1:],
                          out["D"][1:], out["My"][1:], cl_strip[1:], y_strip, gam, wakes, S, cref, U, rho, dt,
                          {"free_wake": free_wake, "nsteps": nsteps, "x_ref": xr})

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _wake_segments(wakes, wgam):
        A, B, G = [np.zeros((0, 3))], [np.zeros((0, 3))], [np.zeros(0)]
        for W, g in zip(wakes, wgam):
            if g.shape[0] == 0:
                continue
            a, b, s = grid_segments(W, g)
            A.append(a); B.append(b); G.append(s)
        return np.concatenate(A), np.concatenate(B), np.concatenate(G)

    @staticmethod
    def _bound_all(Rg, boxes, gam):
        A, B, G = [], [], []
        for R, b in zip(Rg, boxes):
            a, bb, s = grid_segments(R, gam[b])
            A.append(a); B.append(bb); G.append(s)
        return np.concatenate(A), np.concatenate(B), np.concatenate(G)

    def _kj_lift_per_ring(self, Rg, boxes, gam, vfree, allA, allB, allG, rho, sym):
        """KJ lift of each ring's leading segment (strength Gamma_i - Gamma_ahead)."""
        lift = np.zeros(self.lat.mesh.n)
        for R, b in zip(Rg, boxes):
            nc, ns = b.shape
            A_ = R[:-1, :-1].reshape(-1, 3)
            B_ = R[:-1, 1:].reshape(-1, 3)
            g = gam[b]
            s = (g - np.vstack([np.zeros((1, ns)), g[:-1]])).ravel()
            mid = 0.5 * (A_ + B_)
            v = vfree(mid) + segment_velocity(mid, allA, allB, allG, self.core, sym)
            F = rho * s[:, None] * np.cross(v, B_ - A_)
            lift[b.ravel()] += F[:, 2]
        return lift


# ======================================================================
# Linearised UVLM (flat frozen wake)
# ======================================================================
class LinearUVLM:
    """UVLM linearised about a flat wake of ``wake_chords`` reference chords.

    The wake rows are ``dx_wake = U dt`` long; with ``dt = dx_wake / U`` the
    wake strengths are a delay line ``Gamma_w[j](n+1) = Gamma_w[j-1](n)``,
    ``Gamma_w[0](n+1) = Gamma_TE(n)``.  Pressure on ring ``i``::

        dp_i = rho [ U (Gamma_i - Gamma_ahead) w_i / A_i + dGamma_i/dt ]

    (Katz & Plotkin 2001, eq. 13.37, linearised); the first part acts at the
    box 1/4 chord (``h_f``), the second at the ring centre (``h_c``).
    """

    def __init__(self, mesh: AeroMesh, modes=None, bref: float = 1.0, wake_chords: float = 30.0,
                 dx_wake: float | None = None, core: float | None = None):
        self.mesh = mesh
        self.lat = lat = Lattice.from_mesh(mesh)
        self.modes = modes
        self.bref = float(bref)
        dx = np.mean([np.mean(g["G"][-1, :, 0] - g["G"][-2, :, 0]) for g in lat.grids])
        self.dx_wake = float(dx) if dx_wake is None else float(dx_wake)
        self.c_ref = lat.s_ref / float(lat.width[lat.te].sum())
        self.nrow = max(2, int(np.ceil(wake_chords * self.c_ref / self.dx_wake)))
        self.core = 1e-6 * self.c_ref if core is None else core
        # wake rings: per surface rows x ns, ordered (surface, row, strip)
        wr, wstrip = [], []
        k0 = 0
        for g in lat.grids:
            R_te = g["R"][-1]
            ns = R_te.shape[0] - 1
            j = np.arange(self.nrow + 1)
            W = R_te[None, :, :] + (j * self.dx_wake)[:, None, None] * np.array([1.0, 0.0, 0.0])
            wr.append(grid_rings(W))
            wstrip.append((k0 + np.tile(np.arange(ns), self.nrow), np.repeat(np.arange(self.nrow), ns)))
            k0 += ns
        self.wrings = np.concatenate(wr)
        self.w_strip = np.concatenate([w[0] for w in wstrip])
        self.w_row = np.concatenate([w[1] for w in wstrip])
        N = mesh.n
        self.A = ring_normal_matrix(lat.colloc, lat.normal, lat.rings, self.core, mesh.symmetry)
        self.Aw = ring_normal_matrix(lat.colloc, lat.normal, self.wrings, self.core, mesh.symmetry)
        Ks = np.diag(lat.width.copy())
        ahead = lat.prev >= 0
        Ks[np.where(ahead)[0], lat.prev[ahead]] = -lat.width[ahead]
        self.Ks = Ks                                     # F_KJ / (rho U) = Ks Gamma
        self.N = N
        self.Nw = self.wrings.shape[0]

    # ---------------------------------------------------------- frequency domain
    def _wake_sum(self, z):
        """N x N contribution of the harmonic wake (Gamma_w[j] = Gamma_TE z^-(j+1))."""
        fac = z ** (-(self.w_row + 1.0))
        cols = self.Aw * fac[None, :]
        out = np.zeros((self.N, self.N), complex)
        te_box = self.lat.te[self.w_strip]
        np.add.at(out.T, te_box, cols.T)
        return out

    def gaf_table(self, kvals, modes=None, gust: bool = True, keep_cp: bool = True, x_gust0=None):
        """GAF table Q(k) (same convention as :meth:`AeroModel.table`)."""
        from .gaf import GAFTable
        modes = modes or self.modes
        lat, b = self.lat, self.bref
        kvals = np.atleast_1d(np.asarray(kvals, float))
        nm = modes.nm
        x0 = float(self.mesh.corners[:, :, 0].min()) if x_gust0 is None else x_gust0
        Q = np.zeros((kvals.size, nm, nm), complex)
        Qg = np.zeros((kvals.size, nm), complex) if gust else None
        cp = np.zeros((kvals.size, self.N, nm), complex) if keep_cp else None
        cpg = np.zeros((kvals.size, self.N), complex) if (keep_cp and gust) else None
        Hf = modes.h_f
        Hc = modes.h_c
        for i, k in enumerate(kvals):
            z = np.exp(1j * k * self.dx_wake / b)
            Aeff = self.A + self._wake_sum(z)
            W = modes.hx_c + 1j * (k / b) * modes.h_c
            if gust:
                wg = -np.cos(self.mesh.gamma) * np.exp(-1j * (k / b) * (lat.colloc[:, 0] - x0))
                W = np.hstack([W, wg[:, None]])
            G = np.linalg.solve(Aeff, W)                     # Gamma / U
            Fkj = 2.0 * (self.Ks @ G)                        # F / q, at 1/4 chord
            Fu = 2.0 * (1j * k / b) * lat.ring_area[:, None] * G   # F / q, ring centre
            Qall = Hf.T @ Fkj + Hc.T @ Fu
            Q[i] = Qall[:, :nm]
            if gust:
                Qg[i] = Qall[:, nm]
            if keep_cp:
                dcp = (Fkj + Fu) / self.mesh.area[:, None]
                cp[i] = dcp[:, :nm]
                if gust:
                    cpg[i] = dcp[:, nm]
        return GAFTable(kvals, Q, Qg, 0.0, b, list(modes.names), cp, cpg)

    # ---------------------------------------------------------- discrete time
    def state_space(self, U: float, rho: float, M, K, C=None, modes=None,
                    beta: float = 0.25, gamma_nm: float = 0.5):
        """Monolithic discrete-time system ``E x[n+1] = F x[n] + Bg wg[n+1]``.

        ``x = [Gamma_b (N), Gamma_w (Nw), q, q_dot, q_ddot]``; Newmark-beta
        (average acceleration) structure; gust input ``wg`` = uniform vertical
        gust velocity at all collocation points per unit (frozen-gust delays are
        applied by :meth:`simulate`).  Returns ``(E, F, dt)``.
        """
        modes = modes or self.modes
        M = np.atleast_2d(M); K = np.atleast_2d(K)
        m = M.shape[0]
        C = np.zeros_like(M) if C is None else np.atleast_2d(C)
        N, Nw = self.N, self.Nw
        dt = self.dx_wake / U
        n = N + Nw + 3 * m
        E = np.zeros((n, n)); F = np.zeros((n, n))
        ib, iw, iq, iv, ia = 0, N, N + Nw, N + Nw + m, N + Nw + 2 * m
        Hf, Hc, Hx = modes.h_f, modes.h_c, modes.hx_c
        # non-penetration
        E[ib:iw, ib:iw] = self.A
        E[ib:iw, iw:iq] = self.Aw
        E[ib:iw, iq:iv] = -U * Hx
        E[ib:iw, iv:ia] = -Hc
        # wake delay line
        E[iw:iq, iw:iq] = np.eye(Nw)
        row0 = np.where(self.w_row == 0)[0]
        F[iw + row0, ib + self.lat.te[self.w_strip[row0]]] = 1.0
        later = np.where(self.w_row > 0)[0]
        key = {(s, r): i for i, (s, r) in enumerate(zip(self.w_strip, self.w_row))}
        src = np.array([key[(self.w_strip[i], self.w_row[i] - 1)] for i in later], int)
        F[iw + later, iw + src] = 1.0
        # structure  M a + C v + K q - Q_aero = 0
        Ar = self.lat.ring_area
        E[iq:iv, iq:iv] = K
        E[iq:iv, iv:ia] = C
        E[iq:iv, ia:] = M
        E[iq:iv, ib:iw] = -(rho * U * Hf.T @ self.Ks + rho / dt * (Hc.T * Ar[None, :]))
        F[iq:iv, ib:iw] = -(rho / dt) * (Hc.T * Ar[None, :])
        # Newmark
        I = np.eye(m)
        E[iv:ia, iq:iv] = I
        E[iv:ia, ia:] = -dt * dt * beta * I
        F[iv:ia, iq:iv] = I
        F[iv:ia, iv:ia] = dt * I
        F[iv:ia, ia:] = dt * dt * (0.5 - beta) * I
        E[ia:, iv:ia] = I
        E[ia:, ia:] = -dt * gamma_nm * I
        F[ia:, iv:ia] = I
        F[ia:, ia:] = dt * (1 - gamma_nm) * I
        self._idx = (ib, iw, iq, iv, ia, m)
        return E, F, dt

    def eigen(self, U, rho, M, K, C=None, modes=None, guess=None, nev: int = 4):
        """Structural eigenvalues ``lambda = ln(z) / dt`` of the coupled discrete-time system.

        Sparse shift-invert around ``z0 = exp(lambda_guess dt)`` for each
        structural root (``guess`` defaults to the in-vacuo frequencies
        ``i omega_j``; pass the previous speed's roots to track branches).
        Returns ``lam (m,)`` ordered like ``guess``.
        """
        import scipy.sparse as sp
        from scipy.sparse.linalg import LinearOperator, eigs, splu
        E, F, dt = self.state_space(U, rho, M, K, C, modes)
        ib, iw, iq, iv, ia, m = self._idx
        if guess is None:
            from scipy.linalg import eigh
            w2 = eigh(np.atleast_2d(K), np.atleast_2d(M), eigvals_only=True)
            guess = 1j * np.sqrt(np.maximum(w2, 1e-12))
        Es, Fs = sp.csc_matrix(E), sp.csc_matrix(F)
        n = E.shape[0]
        out = np.zeros(len(guess), complex)
        for j, lg in enumerate(guess):
            z0 = np.exp(lg * dt)
            lu = splu((Fs - z0 * Es).tocsc())
            op = LinearOperator((n, n), matvec=lambda x, lu=lu: lu.solve(Es @ x), dtype=complex)
            mu, vec = eigs(op, k=min(nev, n - 2), which="LM", tol=1e-10)
            z = z0 + 1.0 / mu
            lam = np.log(z) / dt
            struct = np.linalg.norm(vec[iq:iv], axis=0) * np.abs(lam)
            # nearest root to the prediction; wake-dominated roots are rejected
            score = np.abs(lam - lg) / max(abs(lg), 1e-9)
            score[struct < 0.1 * np.max(struct)] += 10.0
            score[lam.imag <= 0] += 10.0
            for i in np.argsort(score):                  # skip roots taken by another branch
                if not np.any(np.abs(out[:j] - lam[i]) < 1e-6 * abs(lam[i])):
                    out[j] = lam[i]
                    break
        return out

    def simulate(self, U, rho, M, K, C=None, t_end=1.0, q0=None, qd0=None, gust: Callable | None = None,
                 modes=None, x_gust0=None):
        """Time response of the coupled linear system (initial condition and/or frozen gust).

        ``gust(s)`` is the vertical gust velocity versus penetration distance
        ``s = U t - (x - x_gust0)``.  Returns dict with t, q, CL (rigid
        lift coefficient of the modelled half, positive up) and Gamma_TE sum.
        """
        from scipy.linalg import lu_factor, lu_solve
        modes = modes or self.modes
        E, F, dt = self.state_space(U, rho, M, K, C, modes)
        ib, iw, iq, iv, ia, m = self._idx
        lu = lu_factor(E)
        nt = int(round(t_end / dt)) + 1
        x = np.zeros(E.shape[0])
        if q0 is not None:
            x[iq:iv] = q0
        if qd0 is not None:
            x[iv:ia] = qd0
        x0 = float(self.mesh.corners[:, :, 0].min()) if x_gust0 is None else x_gust0
        nz = self.mesh.normal[:, 2]
        xc = self.lat.colloc[:, 0]
        q = 0.5 * rho * U * U
        S = self.lat.s_ref
        Ar = self.lat.ring_area
        hist_q = np.zeros((nt, m))
        CL = np.zeros(nt)
        hist_q[0] = x[iq:iv]
        for n in range(1, nt):
            rhs = F @ x
            if gust is not None:
                rhs[ib:iw] -= gust(U * n * dt - (xc - x0)) * nz
            x_new = lu_solve(lu, rhs)
            g_new, g_old = x_new[ib:iw], x[ib:iw]
            lift = rho * U * (self.Ks @ g_new) + rho * Ar * (g_new - g_old) / dt
            CL[n] = float((lift * nz).sum()) / (q * S)
            x = x_new
            hist_q[n] = x[iq:iv]
        return {"t": np.arange(nt) * dt, "q": hist_q, "CL": CL, "dt": dt}


def uvlm_gaf_table(mesh: AeroMesh, modes, bref: float, kvals, gust: bool = True, keep_cp: bool = True,
                   wake_chords: float = 30.0):
    """Frequency-domain GAFs of the linearised UVLM (drop-in for :meth:`AeroModel.table`)."""
    return LinearUVLM(mesh, modes, bref, wake_chords).gaf_table(kvals, gust=gust, keep_cp=keep_cp)


def flutter_eigen_sweep(lin: LinearUVLM, rho, M, K, velocities, C=None, modes=None):
    """Discrete-time eigenvalue sweep with branch tracking. Returns (V, lam (nV, m)).

    Each speed is predicted by linear extrapolation of the two previous
    roots; when a root lands far from its prediction (a likely branch swap
    between close frequencies) the speed step is bisected, down to 1/16 of
    the requested spacing, and the intermediate roots are used only for
    tracking.
    """
    V = np.asarray(velocities, float)
    tol, max_depth = 0.08, 3

    def solve(v, guess):
        return lin.eigen(v, rho, M, K, C, modes, guess=guess)

    def advance(v0, l0, slope, v1, depth=0):
        pred = None if l0 is None else l0 + (0 if slope is None else slope * (v1 - v0))
        g = solve(v1, pred)
        if pred is not None and depth < max_depth and np.any(np.abs(g - pred) > tol * np.abs(pred)):
            vm = 0.5 * (v0 + v1)
            gm = advance(v0, l0, slope, vm, depth + 1)
            return advance(vm, gm, (gm - l0) / (vm - v0), v1, depth + 1)
        return g

    lam, prev, slope, vprev = [], None, None, None
    for v in V:
        g = solve(v, None) if prev is None else advance(vprev, prev, slope, v)
        if prev is not None:
            slope = (g - prev) / (v - vprev)
        lam.append(g)
        prev, vprev = g, v
    return V, np.array(lam)


def flutter_crossing(V, lam):
    """Linear interpolation of the first zero crossing of max Re(lambda)."""
    g = lam.real
    best = None
    for j in range(g.shape[1]):
        idx = np.where((g[:-1, j] < 0) & (g[1:, j] >= 0))[0]
        if idx.size:
            i = idx[0]
            f = g[i, j] / (g[i, j] - g[i + 1, j])
            Vf = V[i] + f * (V[i + 1] - V[i])
            w = lam[i, j].imag + f * (lam[i + 1, j].imag - lam[i, j].imag)
            if best is None or Vf < best[0]:
                best = (float(Vf), float(w), j)
    return best
