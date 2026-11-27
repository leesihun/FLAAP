"""Time-domain manoeuvre loads of a flexible aircraft.

The structural model combines rigid-body modes (plunge, pitch, roll) with
elastic modes through the *full* FE mass matrix (fixed-root elastic modes
+ rigid modes is the classical "attached-modes" model; free-free elastic
modes are handled the same way, the rigid/elastic mass coupling simply
vanishes).  Aerodynamics: the RFA state-space model of
:mod:`flaap.timedomain` with the control surfaces as prescribed inputs.

Linear perturbation about 1-g level flight:

* plunge ``z`` (up) and pitch ``theta`` (nose up) are inertial
  displacements; the angle of attack follows from the normalwash of the
  pitch mode (``alpha = theta - z'/V``) - short-period dynamics,
* the incremental load factor is ``dn_z = z'' / g`` at the reference point,
* station loads are incremental (add the 1-g trim loads from
  :mod:`flaap.trim` for totals).

Implements the CS-25.331(c) pitch manoeuvres:

``(c)(1)`` abrupt (unchecked) pitch control displacement at V_A,
``(c)(2)`` checked manoeuvre ``delta(t) = delta_1 sin(w t)`` for
``0 <= t <= t_max = 3 pi / (2 w)`` (then returned linearly to neutral),
with ``delta_1`` scaled so that the limit load factor is just reached.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .structures.modal import rigid_body_shapes
from .trim import G0


@dataclass
class FlexibleAircraft:
    """Generalized structural matrices of rigid + elastic modes."""

    M: np.ndarray
    K: np.ndarray
    C: np.ndarray
    shapes: np.ndarray           # (Nn, 6, ns)
    names: list[str]
    nodes: np.ndarray
    fe: object = None
    n_rigid: int = 0
    extra: dict = field(default_factory=dict)

    @property
    def ns(self) -> int:
        return self.M.shape[0]

    def dof_vectors(self):
        return self.shapes.reshape(-1, self.ns)


def flexible_aircraft(elastic, ref, kinds=("plunge", "pitch"), body_mass: dict | None = None,
                      damping: float | None = None) -> FlexibleAircraft:
    """Combine rigid modes about ``ref`` with the elastic modes of ``elastic``.

    ``body_mass``: lumped rigid-body mass added to the rigid modes,
    ``{'m': kg, 'x': cg x, 'z': cg z, 'Iyy': pitch inertia about its cg,
    'Ixx': roll inertia}`` (for a half model give half values).
    """
    fe = elastic.fe
    if fe is None:
        raise ValueError("elastic ModalModel must carry its FE model")
    ref = np.asarray(ref, float)
    R = rigid_body_shapes(fe.nodes, ref, kinds).reshape(-1, len(kinds))
    E = elastic.dof_vectors()
    Phi = np.hstack([R, E])
    M = np.asarray(Phi.T @ (fe.M @ Phi))
    K = np.asarray(Phi.T @ (fe.K @ Phi))
    nr = len(kinds)
    K[:nr, :] = 0.0
    K[:, :nr] = 0.0
    if body_mass:
        m = body_mass.get("m", 0.0)
        dx = body_mass.get("x", ref[0]) - ref[0]
        for i, a in enumerate(kinds):
            for j, b in enumerate(kinds):
                v = 0.0
                if a == b == "plunge":
                    v = m
                elif {a, b} == {"plunge", "pitch"}:
                    v = -m * dx
                elif a == b == "pitch":
                    v = body_mass.get("Iyy", 0.0) + m * dx * dx
                elif a == b == "roll":
                    v = body_mass.get("Ixx", 0.0)
                M[i, j] += v
    C = np.zeros_like(M)
    zeta = elastic.damping if damping is None else np.full(elastic.nmodes, damping)
    C[nr:, nr:] = np.diag(2 * zeta * elastic.omega * elastic.gen_mass)
    names = list(kinds) + (list(elastic.labels) if elastic.labels else
                           [f"mode{i + 1}" for i in range(elastic.nmodes)])
    shapes = Phi.reshape(fe.nodes.shape[0], 6, -1)
    return FlexibleAircraft(M, K, C, shapes, names, fe.nodes, fe, nr)


# ----------------------------------------------------------------------
# control inputs
# ----------------------------------------------------------------------
def checked_pitch_input(t, omega: float, amp: float = 1.0, return_to_neutral: bool = True):
    """CS-25.331(c)(2) checked control input and its two derivatives."""
    t = np.asarray(t, float)
    tmax = 1.5 * np.pi / omega
    d = np.where(t <= tmax, amp * np.sin(omega * t), 0.0)
    dd = np.where(t <= tmax, amp * omega * np.cos(omega * t), 0.0)
    if return_to_neutral:
        tr = 0.5 * np.pi / omega
        on = (t > tmax) & (t <= tmax + tr)
        d = np.where(on, -amp * (1 - (t - tmax) / tr), d)
        dd = np.where(on, amp / tr, dd)
    ddd = np.gradient(dd, t)
    return d, dd, ddd


def abrupt_input(t, amp: float = 1.0, t_ramp: float = 0.1, t_hold: float | None = None):
    """Unchecked (c)(1) control input: ramp to ``amp`` in ``t_ramp`` and hold."""
    t = np.asarray(t, float)
    d = amp * np.clip(t / t_ramp, 0, 1)
    if t_hold is not None:
        d = np.where(t > t_ramp + t_hold, amp * np.clip(1 - (t - t_ramp - t_hold) / t_ramp, 0, 1), d)
    dd = np.gradient(d, t)
    ddd = np.gradient(dd, t)
    return d, dd, ddd


def short_period(system, V: float, rho: float):
    """Short-period root ``(wn, zeta, lambda)`` of the state space.

    Among the oscillatory roots, the one whose displacement eigenvector
    is dominated by the rigid-body modes (largest rigid fraction of
    ``|eta|``; ties broken by pitch participation) - the elastic and
    aerodynamic-lag roots have little rigid content.
    """
    ss = system.state_space(V, rho, gust=False, outputs=False)
    lam, vec = np.linalg.eig(ss.A)
    ns = system.ns
    nr = getattr(system, "n_rigid", None)
    if nr is None:
        nr = sum(1 for n in system.names[:ns] if n in ("plunge", "pitch", "roll", "yaw", "lateral", "surge"))
    ip = system.names.index("pitch")
    best, score = None, -1.0
    for l, v in zip(lam, vec.T):
        if l.imag <= 1e-6:
            continue
        eta = v[:ns]
        nrm = np.linalg.norm(eta)
        if nrm == 0:
            continue
        rigid = np.linalg.norm(eta[:max(nr, ip + 1)]) / nrm
        sc = rigid + 1e-3 * abs(eta[ip]) / nrm
        if sc > score:
            best, score = l, sc
    lam0 = best if best is not None else lam[np.argmax(lam.real)]
    wn = abs(lam0)
    return wn, -lam0.real / wn, lam0


@dataclass
class ManeuverResult:
    t: np.ndarray
    y: np.ndarray
    labels: list[str]
    delta: np.ndarray
    nz: np.ndarray
    scale: float = 1.0

    def series(self, name):
        return self.y[:, self.labels.index(name)]

    def peak(self, name):
        s = self.series(name)
        return float(s.max()), float(s.min())


def simulate_control(system, V: float, rho: float, control: str, d, dd, ddd, t,
                     g: float = G0) -> ManeuverResult:
    ss = system.state_space(V, rho, gust=False)
    u = np.zeros((t.size, len(ss.inputs)))
    i = ss.inputs.index(control)
    u[:, i], u[:, i + 1], u[:, i + 2] = d, dd, ddd
    y, _ = ss.simulate(t, u)
    nz = 1.0 + y[:, ss.outputs.index("etadd:plunge")] / g
    return ManeuverResult(t, y, list(ss.outputs), d, nz)


def checked_pitch_maneuver(system, V: float, rho: float, control: str, nz_limit: float = 2.5,
                           omega: float | None = None, V_A: float | None = None,
                           t_end: float | None = None, dt: float | None = None,
                           g: float = G0) -> ManeuverResult:
    """CS-25.331(c)(2) checked manoeuvre scaled to reach ``nz_limit``.

    ``omega`` defaults to the short-period frequency, but not less than
    ``pi V / (4 V_A)`` when ``V_A`` is given.
    """
    if omega is None:
        omega = short_period(system, V, rho)[0]
        if V_A is not None:
            omega = max(omega, np.pi * V / (4 * V_A))
    t_end = t_end or (3.5 * np.pi / omega + 2.0)
    dt = dt or min(0.005, 2 * np.pi / omega / 200)
    t = np.arange(0.0, t_end + dt / 2, dt)
    d, dd, ddd = checked_pitch_input(t, omega, 1.0)
    r = simulate_control(system, V, rho, control, d, dd, ddd, t, g)
    dn = r.nz - 1.0
    peak = dn[np.argmax(np.abs(dn))]
    scale = (nz_limit - 1.0) / peak
    r.y *= scale
    r.delta = d * scale
    r.nz = 1.0 + dn * scale
    r.scale = scale
    return r
