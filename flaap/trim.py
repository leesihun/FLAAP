"""Static aeroelasticity: trim, divergence, control effectiveness / reversal.

Uses the steady (k = 0) generalized aerodynamic forces of a set of
"variables" ``xi`` that are columns of the GAF matrix: structural modes
(elastic and rigid-body, ``ns`` of them, first) followed by trim
parameters (``alpha``, ``pitch_rate``, ``roll_rate``, control surfaces,
...).  The steady equations (one per structural mode) are::

    K eta - q Q0[s, :] xi = -M a                     (a: rigid accelerations)

where the acceleration vector ``a`` contains gravity + inertia, e.g. for
a symmetric manoeuvre with plunge mode ``h = 1`` and load factor
``n_z``: ``a_plunge = n_z g`` (gravity is the uniform field ``-g`` along the
plunge mode).  For free-free elastic modes the right-hand side vanishes
on the elastic rows (inertia relief is automatic); for clamped
("fixed-root") elastic modes coupled with rigid modes through a full
mass matrix it does not.

Each variable is either *free* (solved for) or *fixed* (given);
the number of free variables must equal the number of equations
(``ns``).  This is the Nastran SOL 144 formulation in modal coordinates.

Also provided:

* :func:`divergence` - lowest positive dynamic pressure with
  ``det(K - q Q0_ee) = 0``,
* :func:`control_effectiveness` - flexible-to-rigid ratio of a control
  derivative (aileron effectiveness), and :func:`reversal` (zero
  effectiveness dynamic pressure).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

G0 = 9.80665


@dataclass
class TrimResult:
    q: float
    names: list[str]
    xi: np.ndarray                     # all variables
    cp: np.ndarray | None = None       # box dCp (steady)
    loads: dict = field(default_factory=dict)
    forces: np.ndarray | None = None   # generalized aero forces q Q0 xi (structural rows)

    def __getitem__(self, name):
        return float(self.xi[self.names.index(name)])

    def summary(self) -> str:
        lines = [f"Static aeroelastic trim, q = {self.q:.6g} Pa"]
        for n, v in zip(self.names, self.xi):
            lines.append(f"  {n:>16s} = {v: .6e}")
        for k, v in self.loads.items():
            lines.append(f"  {k:>16s} = {v: .6e}")
        return "\n".join(lines)


class StaticAeroelastic:
    """Steady aeroelastic model.

    Parameters
    ----------
    M, K : (ns, ns) generalized mass / stiffness of the structural modes
    Q0   : (ns, nv) steady GAF rows of the structural modes for all
           variables (``nv >= ns``; structural modes first)
    names : variable names (``nv``)
    cp0  : optional (N, nv) steady box pressures per unit variable
    load_rows : optional dict ``label -> (nv,)`` aero load rows per unit q
    inertia_rows : optional dict ``label -> (ns,)`` mapping rigid/elastic
           accelerations ``a`` to station loads (``-T M_fe Phi``)
    """

    def __init__(self, M, K, Q0, names, cp0=None, load_rows=None, inertia_rows=None,
                 stiffness_rows=None):
        self.M = np.atleast_2d(np.asarray(M, float))
        self.K = np.atleast_2d(np.asarray(K, float))
        self.Q0 = np.real(np.asarray(Q0))
        self.names = list(names)
        self.ns = self.M.shape[0]
        self.nv = self.Q0.shape[1]
        if len(self.names) != self.nv:
            raise ValueError("names must label every GAF column")
        self.cp0 = None if cp0 is None else np.real(np.asarray(cp0))
        self.load_rows = load_rows or {}
        self.inertia_rows = inertia_rows or {}
        self.stiffness_rows = stiffness_rows or {}

    @classmethod
    def from_table(cls, M, K, table, **kw):
        i0 = int(np.argmin(np.abs(table.k)))
        if abs(table.k[i0]) > 1e-9:
            raise ValueError("GAF table needs k = 0 for static aeroelasticity")
        ns = np.atleast_2d(M).shape[0]
        cp0 = None if table.cp is None else table.cp[i0]
        return cls(M, K, table.Q[i0][:ns], table.names, cp0=cp0, **kw)

    def index(self, name) -> int:
        return self.names.index(name)

    # ------------------------------------------------------------------
    def solve(self, q: float, fixed: dict | None = None, free=None, accel: dict | None = None,
              stations=True) -> TrimResult:
        """Solve the trim equations.

        ``fixed``: ``{name: value}`` for prescribed variables (others
        default to 0 unless listed in ``free``).  ``free``: names of the
        unknown trim variables; all structural modes that are not fixed are
        free automatically.  ``accel``: ``{mode name: acceleration}`` (e.g.
        ``{'plunge': nz * g}``).
        """
        fixed = dict(fixed or {})
        free = list(free or [])
        unknown = [n for n in self.names[:self.ns] if n not in fixed] + \
                  [n for n in free if n not in self.names[:self.ns]]
        if len(unknown) != self.ns:
            raise ValueError(f"{len(unknown)} unknowns {unknown} for {self.ns} equations")
        a = np.zeros(self.ns)
        for n, v in (accel or {}).items():
            a[self.index(n)] = v
        Kfull = np.zeros((self.ns, self.nv))
        Kfull[:, :self.ns] = self.K
        Afull = Kfull - q * self.Q0                 # residual = Afull xi + M a = 0
        iu = [self.index(n) for n in unknown]
        xf = np.zeros(self.nv)
        for n, v in fixed.items():
            xf[self.index(n)] = v
        rhs = -self.M @ a - Afull @ xf
        sol = np.linalg.solve(Afull[:, iu], rhs)
        xi = xf.copy()
        xi[iu] = sol
        res = TrimResult(q, self.names, xi)
        res.forces = q * self.Q0 @ xi
        if self.cp0 is not None:
            res.cp = self.cp0 @ xi
        if stations:
            for lab, row in self.load_rows.items():
                v = q * np.real(row) @ xi
                if lab in self.inertia_rows:
                    v += self.inertia_rows[lab] @ a
                res.loads[lab] = float(v)
            for lab, row in self.stiffness_rows.items():
                res.loads[lab + ".mdm"] = float(row @ xi[:self.ns])
        return res

    # ------------------------------------------------------------------
    def elastic_subset(self, elastic):
        idx = [self.index(n) if isinstance(n, str) else int(n) for n in elastic]
        return idx

    def divergence(self, elastic=None) -> float | None:
        """Divergence dynamic pressure of the (clamped / fixed-trim) system."""
        idx = list(range(self.ns)) if elastic is None else self.elastic_subset(elastic)
        return divergence(self.K[np.ix_(idx, idx)], self.Q0[np.ix_(idx, idx)])

    def effectiveness(self, q, control: str, row, elastic=None, rigid_vars=()):
        """Ratio of flexible to rigid derivative of output ``row`` (nv,) w.r.t.
        ``control``, with the variables in ``rigid_vars`` held at zero and
        all structural modes in ``elastic`` free."""
        idx = list(range(self.ns)) if elastic is None else self.elastic_subset(elastic)
        return control_effectiveness(self.K[np.ix_(idx, idx)], self.Q0[np.ix_(idx, idx)],
                                     self.Q0[idx, self.index(control)], np.real(row)[idx],
                                     float(np.real(row)[self.index(control)]), q)


# ----------------------------------------------------------------------
def divergence(K, Q0) -> float | None:
    """Smallest positive real ``q`` with ``det(K - q Q0) = 0``."""
    lam = np.linalg.eigvals(np.linalg.solve(K, np.real(Q0)))     # K^-1 Q0 v = (1/q) v
    lam = lam[np.abs(lam.imag) < 1e-8 * np.maximum(np.abs(lam), 1e-30)].real
    lam = lam[lam > 0]
    return float(1.0 / lam.max()) if lam.size else None


def control_effectiveness(Kee, Qee, Qec, row_e, row_c, q):
    """Flexible / rigid ratio of a control derivative.

    ``Kee, Qee``: elastic stiffness and steady GAF; ``Qec``: elastic GAF
    column of the control; ``row_e``, ``row_c``: the output's aero row
    (e.g. rolling moment) for elastic modes and the control.
    """
    q = np.atleast_1d(np.asarray(q, float))
    out = []
    for qq in q:
        eta = np.linalg.solve(Kee - qq * Qee, qq * Qec)
        out.append(1.0 + (row_e @ eta) / row_c)
    out = np.array(out)
    return out if out.size > 1 else float(out[0])


def reversal(Kee, Qee, Qec, row_e, row_c, q_max: float, n: int = 400) -> float | None:
    """Control reversal dynamic pressure (effectiveness = 0) below ``q_max``."""
    qs = np.linspace(0.0, q_max, n)
    eff = control_effectiveness(Kee, Qee, Qec, row_e, row_c, qs)
    s = np.where(np.sign(eff[1:]) != np.sign(eff[:-1]))[0]
    if not s.size:
        return None
    i = s[0]
    f = lambda qq: control_effectiveness(Kee, Qee, Qec, row_e, row_c, qq)  # noqa: E731
    return float(brentq(f, qs[i], qs[i + 1]))


def lift_coefficient(mesh, cp: np.ndarray, sref: float) -> float:
    """CL of the modelled surfaces from box pressures (z-force / q S)."""
    return float(np.sum(cp * mesh.area * mesh.normal[:, 2]) / sref)


def load_rows_from_table(table, mesh, stations) -> dict:
    """Steady aero station rows ``label -> (nv,)`` per unit q from a table."""
    from .loads import aero_load_matrix
    i0 = int(np.argmin(np.abs(table.k)))
    rows = {}
    for st in stations:
        T = aero_load_matrix(mesh, st)
        R = T @ np.real(table.cp[i0])
        for lab, r in zip(st.labels, R):
            rows[lab] = r
    return rows
