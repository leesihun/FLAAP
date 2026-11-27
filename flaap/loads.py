"""Integrated loads at monitor stations.

A :class:`MonitorStation` is a cutting plane through ``point`` with
normal ``normal``; everything on the positive side (``(r - point) . normal
> 0``) is the "outboard" free body.  The station load is the resultant of
the external forces on that free body, expressed in global axes about
``point``::

    [Fx, Fy, Fz, Mx, My, Mz]

For a right wing root cut (normal ``+y``) ``Fz`` is the vertical shear,
``Mx`` the bending moment (positive = tip up) and ``My`` the torque
(positive = nose up... about +y, i.e. nose *down* is negative).

Two recovery methods are provided:

* **summation of forces** (force summation, FSM): aerodynamic box forces
  ``q dCp A n`` plus nodal inertia forces ``-M_fe Phi eta''`` (and any
  gravity / rigid-body acceleration field),
* **mode displacement** (MDM): elastic forces ``K_fe Phi eta``.

Both are linear in the modal / aerodynamic quantities and therefore are
supplied as matrices that can be attached to frequency-domain responses
or to the time-domain state-space model as output rows
(:class:`flaap.timedomain.OutputGroup`).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

COMPONENTS = ("Fx", "Fy", "Fz", "Mx", "My", "Mz")


@dataclass
class MonitorStation:
    name: str
    point: np.ndarray
    normal: np.ndarray = None
    components: tuple = COMPONENTS
    tol: float = 1e-9

    def __post_init__(self):
        self.point = np.asarray(self.point, float)
        self.normal = np.array([0.0, 1.0, 0.0]) if self.normal is None else np.asarray(self.normal, float)
        self.normal = self.normal / np.linalg.norm(self.normal)

    def outboard(self, pts: np.ndarray) -> np.ndarray:
        return (np.asarray(pts, float) - self.point) @ self.normal > self.tol

    @property
    def labels(self) -> list[str]:
        return [f"{self.name}.{c}" for c in self.components]

    def _sel(self):
        return [COMPONENTS.index(c) for c in self.components]


def _cross_rows(d: np.ndarray, f_dir: np.ndarray) -> np.ndarray:
    """(N,3) moment arms x (N,3) unit forces -> (N,3) moments."""
    return np.cross(d, f_dir)


def aero_load_matrix(mesh, station: MonitorStation) -> np.ndarray:
    """(ncomp, N) matrix T with station loads = q T dCp (dCp per box)."""
    on = station.outboard(mesh.center)
    F = (mesh.area[:, None] * mesh.normal) * on[:, None]          # per unit dCp
    Mo = _cross_rows(mesh.force_point - station.point, F)
    T = np.vstack([F.T, Mo.T])
    return T[station._sel()]


def nodal_load_matrix(nodes: np.ndarray, station: MonitorStation) -> np.ndarray:
    """(ncomp, 6 Nn) matrix summing nodal forces/moments of outboard nodes."""
    nodes = np.asarray(nodes, float)
    nn = nodes.shape[0]
    on = station.outboard(nodes)
    d = nodes - station.point
    T = np.zeros((6, 6 * nn))
    for j in np.where(on)[0]:
        b = 6 * j
        T[0:3, b:b + 3] = np.eye(3)
        # moment of a unit force e_i at arm d: d x e_i
        dx, dy, dz = d[j]
        T[3:6, b:b + 3] = np.array([[0.0, -dz, dy], [dz, 0.0, -dx], [-dy, dx, 0.0]])
        T[3:6, b + 3:b + 6] = np.eye(3)
    return T[station._sel()]


def _matmul(A, B):
    return np.asarray(A @ B) if not sp.issparse(A) else np.asarray((A @ B))


def modal_load_matrices(modal, station: MonitorStation, fe=None):
    """Inertia and stiffness load matrices for a modal model.

    Returns ``(T_inertia, T_stiffness)`` (ncomp, nm) such that the station
    load is ``T_inertia @ eta''`` (FSM inertia part, ``-T M Phi``) and
    ``T_stiffness @ eta`` (MDM, ``T K Phi``).  ``fe`` defaults to
    ``modal.fe``.
    """
    fe = fe or modal.fe
    if fe is None:
        raise ValueError("modal_load_matrices needs the FE model (modal.fe)")
    Phi = modal.dof_vectors()
    T = nodal_load_matrix(fe.nodes, station)
    MPhi = _matmul(fe.M, Phi)
    KPhi = _matmul(fe.K, Phi)
    if len(fe.fixed):
        # forces at constrained DOFs go straight into the supports and are
        # not loads on the free body (consistent-mass inertia at a clamped
        # node, reactions)
        KPhi = KPhi.copy()
        MPhi = MPhi.copy()
        KPhi[fe.fixed] = 0.0
        MPhi[fe.fixed] = 0.0
    return -T @ MPhi, T @ KPhi


def gravity_load(fe, station: MonitorStation, nz: float = 1.0, g: float = 9.80665) -> np.ndarray:
    """Station load of the outboard mass under a vertical load factor
    ``nz`` (inertia + gravity: nodal force ``-nz g m`` along z)."""
    nn = fe.nodes.shape[0]
    a = np.zeros(6 * nn)
    a[2::6] = -nz * g          # downward body force per unit mass
    T = nodal_load_matrix(fe.nodes, station)
    return T @ _matmul(fe.M, a)


def aero_load_table(table, mesh, stations):
    """Tabulated aerodynamic station loads per unit q.

    Returns ``(labels, rows (nk, nr, nm), rows_g (nk, nr) | None)`` from the
    box pressures stored in a :class:`GAFTable` (``keep_cp=True``).
    """
    if table.cp is None:
        raise ValueError("GAF table has no box pressures (build it with keep_cp=True)")
    T = np.vstack([aero_load_matrix(mesh, s) for s in stations])
    labels = sum((s.labels for s in stations), [])
    rows = np.einsum("rn,knm->krm", T, table.cp)
    rows_g = None if table.cpg is None else np.einsum("rn,kn->kr", T, table.cpg)
    return labels, rows, rows_g


def station_output_group(system, table, mesh, modal, stations, method: str = "fsm"):
    """Attach station loads to an :class:`~flaap.timedomain.AeroelasticSystem`.

    ``method``: ``'fsm'`` (aero + inertia, summation of forces),
    ``'mdm'`` (mode displacement) or ``'both'`` (two groups, suffixes
    ``.fsm`` / ``.mdm``).
    """
    from .timedomain import OutputGroup
    labels = sum((s.labels for s in stations), [])
    ns = system.ns
    groups = []
    if method in ("fsm", "both"):
        lab, rows, rows_g = aero_load_table(table, mesh, stations)
        Ti = np.vstack([modal_load_matrices(modal, s)[0] for s in stations])[:, :ns]
        lab = [l + (".fsm" if method == "both" else "") for l in labels]
        groups.append(system.add_aero_output(lab, table.k, rows, rows_g, inertia=Ti))
    if method in ("mdm", "both"):
        Tk = np.vstack([modal_load_matrices(modal, s)[1] for s in stations])[:, :ns]
        lab = [l + (".mdm" if method == "both" else "") for l in labels]
        groups.append(system.add_output(OutputGroup(lab, stiffness=Tk)))
    return groups


def correlated_loads(y: np.ndarray, labels, pair) -> dict:
    """Peak values and correlated (time-consistent) values of two load
    components from a time history ``y`` (nt, nout), e.g. bending vs torque
    for loads envelopes / potato plots."""
    i, j = labels.index(pair[0]), labels.index(pair[1])
    out = {}
    for name, idx in (("max_" + pair[0], np.argmax(y[:, i])), ("min_" + pair[0], np.argmin(y[:, i])),
                      ("max_" + pair[1], np.argmax(y[:, j])), ("min_" + pair[1], np.argmin(y[:, j]))):
        out[name] = (float(y[idx, i]), float(y[idx, j]))
    return out
