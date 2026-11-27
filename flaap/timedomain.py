"""Time-domain aeroelastic state-space models built on a Roger RFA.

Equations of motion (modal coordinates ``eta``, prescribed control
rotations ``delta``, vertical gust ``w_g`` at the gust reference point)::

    M eta'' + C eta' + K eta + M_c delta'' = q Q(p) xi + q Q_g(p) w_g / V
    xi = [eta; delta],  p = s b / V,  q = rho V^2 / 2

With the Roger approximation ``Q(p) = A0 + A1 p + A2 p^2 +
sum_l A_{l+3} p / (p + beta_l)`` the aerodynamic lag states are taken
*unscaled* (one per mode and lag)::

    y_l = p / (p + beta_l) xi      ->   y_l' = xi' - (V / b) beta_l y_l
    z_l = p / (p + beta_l) w_g / V ->   z_l' = w_g' / V - (V / b) beta_l z_l

so that any further aerodynamic output fitted with the same lag roots
(box pressures, integrated loads at monitor stations, ...) is a linear
combination of the same states.  The state vector is
``x = [eta, eta', y_1 .. y_L, z_1 .. z_L]`` and the input vector is
``u = [w_g, w_g', (delta_c, delta_c', delta_c'') for each control]``.

This replaces the Fortran ``generate_A_matrices`` / state-space part of
``calculation.f90`` and provides the basis for

* p-method flutter (eigenvalues of ``A(V)``; :func:`flaap.flutter.p_method`),
* discrete and continuous gust responses (:mod:`flaap.gust`),
* time-domain manoeuvre loads (:mod:`flaap.maneuver`).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla

from .aero.rfa import RogerRFA, fit_rows


# ----------------------------------------------------------------------
# generic LTI system
# ----------------------------------------------------------------------
@dataclass
class StateSpace:
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    states: list[str] = field(default_factory=list)
    inputs: list[str] = field(default_factory=list)
    outputs: list[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return self.A.shape[0]

    def index(self, names, kind: str = "outputs") -> list[int]:
        lst = getattr(self, kind)
        names = [names] if isinstance(names, str) else list(names)
        return [lst.index(n) for n in names]

    def eig(self):
        return np.linalg.eigvals(self.A)

    def is_stable(self, tol: float = 1e-9) -> bool:
        return bool(np.all(self.eig().real < tol))

    def freq_response(self, omega, inputs=None, outputs=None) -> np.ndarray:
        """H(i w), shape (nw, nout, nin)."""
        omega = np.atleast_1d(np.asarray(omega, float))
        iu = slice(None) if inputs is None else self.index(inputs, "inputs")
        io = slice(None) if outputs is None else self.index(outputs, "outputs")
        B = self.B[:, iu]
        C = self.C[io]
        D = self.D[io][:, iu]
        I = np.eye(self.n)
        # Hessenberg form speeds up repeated solves
        H, Qh = sla.hessenberg(self.A, calc_q=True)
        Bh = Qh.T @ B
        Ch = C @ Qh
        out = np.empty((omega.size, C.shape[0], B.shape[1]), complex)
        for i, w in enumerate(omega):
            X = np.linalg.solve(1j * w * I - H, Bh)
            out[i] = Ch @ X + D
        return out

    def discretize(self, dt: float):
        """First-order-hold discretisation: returns (Phi, G1, G2) with
        ``x[k+1] = Phi x[k] + G1 u[k] + G2 (u[k+1] - u[k]) / dt``."""
        n, m = self.n, self.B.shape[1]
        Mx = np.zeros((n + 2 * m, n + 2 * m))
        Mx[:n, :n] = self.A
        Mx[:n, n:n + m] = self.B
        Mx[n:n + m, n + m:] = np.eye(m)
        E = sla.expm(Mx * dt)
        return E[:n, :n], E[:n, n:n + m], E[:n, n + m:]

    def simulate(self, t, u, x0=None):
        """Response to input samples ``u`` (nt, nin) on a uniform grid ``t``.

        Exact for piecewise-linear inputs.  Returns ``(y, x)``.
        """
        t = np.asarray(t, float)
        u = np.asarray(u, float)
        if u.ndim == 1:
            u = u[:, None]
        dt = t[1] - t[0]
        if not np.allclose(np.diff(t), dt, rtol=1e-6, atol=1e-12):
            raise ValueError("simulate needs a uniform time grid")
        Phi, G1, G2 = self.discretize(dt)
        x = np.zeros((t.size, self.n))
        if x0 is not None:
            x[0] = x0
        du = np.diff(u, axis=0) / dt
        for k in range(t.size - 1):
            x[k + 1] = Phi @ x[k] + G1 @ u[k] + G2 @ du[k]
        y = x @ self.C.T + u @ self.D.T
        return y, x


# ----------------------------------------------------------------------
# aeroelastic system
# ----------------------------------------------------------------------
@dataclass
class OutputGroup:
    """Extra outputs ``y = q R(p) xi + q R_g(p) w_g/V + inertia eta'' + stiffness eta``.

    ``R`` (3+nlag, nr, nm_total) and ``Rg`` (3+nlag, nr) are Roger
    coefficients with the *same lag roots* as the GAF RFA.
    """
    labels: list[str]
    R: np.ndarray | None = None
    Rg: np.ndarray | None = None
    inertia: np.ndarray | None = None      # (nr, ns)
    stiffness: np.ndarray | None = None    # (nr, ns)
    constant: np.ndarray | None = None     # (nr,) e.g. 1-g trim loads (added by the user)


class AeroelasticSystem:
    """Generalised aeroelastic model ``M, C, K`` + Roger RFA -> state space.

    The RFA must cover ``ns`` structural (incl. rigid-body) modes followed
    by ``n_ctrl`` control-surface modes.  ``M_ctrl`` (ns, n_ctrl) is the
    optional inertial coupling of the control rotations.
    """

    def __init__(self, M, K, rfa: RogerRFA, C=None, n_ctrl: int = 0, M_ctrl=None,
                 names=None):
        self.M = np.atleast_2d(np.asarray(M, float))
        self.K = np.atleast_2d(np.asarray(K, float))
        ns = self.M.shape[0]
        self.C = np.zeros((ns, ns)) if C is None else np.asarray(C, float)
        self.rfa = rfa
        self.nm = rfa.A.shape[1]
        self.ns = ns
        self.nc = n_ctrl
        if self.nm != ns + n_ctrl:
            raise ValueError(f"RFA has {self.nm} modes, expected {ns} + {n_ctrl}")
        self.M_ctrl = np.zeros((ns, n_ctrl)) if M_ctrl is None else np.asarray(M_ctrl, float)
        self.names = list(names) if names is not None else list(rfa.names)
        self.groups: list[OutputGroup] = []

    # ------------------------------------------------------------------
    @property
    def beta(self) -> np.ndarray:
        return self.rfa.beta

    @property
    def bref(self) -> float:
        return self.rfa.bref

    def add_output(self, group: OutputGroup):
        self.groups.append(group)
        return group

    def add_aero_output(self, labels, k, rows, rows_g=None, constrain_k0=True,
                        inertia=None, stiffness=None):
        """Fit tabulated aerodynamic output rows (nk, nr, nm) / (nk, nr) with
        the RFA lag roots and register them as outputs."""
        R = fit_rows(k, rows, self.beta, constrain_k0, True)
        Rg = None if rows_g is None else fit_rows(k, rows_g, self.beta, constrain_k0, False)
        return self.add_output(OutputGroup(list(labels), R, Rg, inertia, stiffness))

    # ------------------------------------------------------------------
    def state_space(self, V: float, rho: float, gust: bool = True,
                    outputs: bool = True) -> StateSpace:
        ns, nc, nm = self.ns, self.nc, self.nm
        L = self.beta.size
        tau = self.bref / V
        q = 0.5 * rho * V * V
        A = self.rfa.A
        Ag = self.rfa.Ag if (gust and self.rfa.Ag is not None) else None
        gust = Ag is not None
        nz = L if gust else 0
        n = 2 * ns + L * nm + nz
        m = (2 if gust else 0) + 3 * nc
        iE, iEd = slice(0, ns), slice(ns, 2 * ns)
        iY = [slice(2 * ns + l * nm, 2 * ns + (l + 1) * nm) for l in range(L)]
        iZ = [2 * ns + L * nm + l for l in range(nz)]
        ug = 0
        uc = (2 if gust else 0)
        iD = [uc + 3 * c for c in range(nc)]          # delta, delta', delta''

        S = slice(0, ns)
        Mbar = self.M - q * tau * tau * A[2][S, S]
        Mi = np.linalg.inv(Mbar)

        # generalized forces (right-hand side for eta'') as linear maps of x and u
        Fx = np.zeros((ns, n))
        Fu = np.zeros((ns, m))
        Fx[:, iE] = -(self.K - q * A[0][S, S])
        Fx[:, iEd] = -(self.C - q * tau * A[1][S, S])
        for l in range(L):
            Fx[:, iY[l]] = q * A[3 + l][S, :]
        if gust:
            for l in range(L):
                Fx[:, iZ[l]] = q * Ag[3 + l][S]
            Fu[:, ug] = q * Ag[0][S] / V
            Fu[:, ug + 1] = q * tau * Ag[1][S] / V
        for c in range(nc):
            j = ns + c
            Fu[:, iD[c]] = q * A[0][S, j]
            Fu[:, iD[c] + 1] = q * tau * A[1][S, j]
            Fu[:, iD[c] + 2] = q * tau * tau * A[2][S, j] - self.M_ctrl[:, c]
        Accx = Mi @ Fx
        Accu = Mi @ Fu

        As = np.zeros((n, n))
        Bs = np.zeros((n, m))
        As[iE, iEd] = np.eye(ns)
        As[iEd] = Accx
        Bs[iEd] = Accu
        for l in range(L):
            rows = np.arange(iY[l].start, iY[l].stop)
            As[rows[:ns], np.arange(ns, 2 * ns)] += 1.0        # xi_s' = eta'
            for c in range(nc):
                Bs[rows[ns + c], iD[c] + 1] += 1.0              # xi_c' = delta'
            As[rows, rows] += -self.beta[l] / tau
        if gust:
            for l in range(L):
                Bs[iZ[l], ug + 1] = 1.0 / V
                As[iZ[l], iZ[l]] = -self.beta[l] / tau

        states = ([f"eta:{s}" for s in self.names[:ns]] + [f"etad:{s}" for s in self.names[:ns]]
                  + [f"lag{l + 1}:{s}" for l in range(L) for s in self.names[:nm]]
                  + [f"glag{l + 1}" for l in range(nz)])
        inputs = (["wg", "wg_dot"] if gust else []) + [
            f"{nm_}{suf}" for nm_ in self.names[ns:nm] for suf in ("", "_dot", "_ddot")]

        # ---- outputs: eta, eta', eta'' and the registered groups
        Cl, Dl, labels = [], [], []
        if outputs:
            I = np.eye(n)
            Cl += [I[iE], I[iEd], Accx]
            Dl += [np.zeros((ns, m)), np.zeros((ns, m)), Accu]
            labels += ([f"eta:{s}" for s in self.names[:ns]] + [f"etad:{s}" for s in self.names[:ns]]
                       + [f"etadd:{s}" for s in self.names[:ns]])
            # xi, xi', xi'' as maps
            X0x = np.zeros((nm, n)); X0u = np.zeros((nm, m))
            X1x = np.zeros((nm, n)); X1u = np.zeros((nm, m))
            X2x = np.zeros((nm, n)); X2u = np.zeros((nm, m))
            X0x[:ns, iE] = np.eye(ns)
            X1x[:ns, iEd] = np.eye(ns)
            X2x[:ns] = Accx
            X2u[:ns] = Accu
            for c in range(nc):
                X0u[ns + c, iD[c]] = 1.0
                X1u[ns + c, iD[c] + 1] = 1.0
                X2u[ns + c, iD[c] + 2] = 1.0
            for g in self.groups:
                nr = len(g.labels)
                Cg = np.zeros((nr, n))
                Dg = np.zeros((nr, m))
                if g.R is not None:
                    R = g.R
                    Cg += q * (R[0] @ X0x + tau * R[1] @ X1x + tau * tau * R[2] @ X2x)
                    Dg += q * (R[0] @ X0u + tau * R[1] @ X1u + tau * tau * R[2] @ X2u)
                    for l in range(L):
                        Cg[:, iY[l]] += q * R[3 + l]
                if gust and g.Rg is not None:
                    Rg = g.Rg
                    Dg[:, ug] += q * Rg[0] / V
                    Dg[:, ug + 1] += q * tau * Rg[1] / V
                    for l in range(L):
                        Cg[:, iZ[l]] += q * Rg[3 + l]
                if g.inertia is not None:
                    Cg += g.inertia @ Accx
                    Dg += g.inertia @ Accu
                if g.stiffness is not None:
                    Cg[:, iE] += g.stiffness
                Cl.append(Cg)
                Dl.append(Dg)
                labels += list(g.labels)
        C = np.vstack(Cl) if Cl else np.zeros((0, n))
        D = np.vstack(Dl) if Dl else np.zeros((0, m))
        return StateSpace(As, Bs, C, D, states, inputs, labels)

    def build_A(self, rho: float, gust: bool = False):
        """Callable ``V -> A(V)`` for :func:`flaap.flutter.p_method`."""
        return lambda V: self.state_space(V, rho, gust=gust, outputs=False).A

    def flutter(self, velocities, rho: float, nstruct: int | None = None, w_initial=None):
        from .flutter import p_method
        ns = self.ns if nstruct is None else nstruct
        if w_initial is None:
            w_initial = np.sqrt(np.maximum(np.real(np.linalg.eigvals(np.linalg.solve(self.M, self.K))), 0))
            w_initial = np.sort(w_initial)
            # rigid-body (zero frequency) roots are not tracked as oscillatory branches
            w_initial = w_initial[w_initial > 1e-6][:ns] if np.any(w_initial > 1e-6) else None
        return p_method(self.build_A(rho), velocities, rho,
                        ns if w_initial is None else len(w_initial), self.bref,
                        names=self.names[:ns], w_initial=w_initial,
                        k_valid=getattr(self.rfa, "kmax", np.inf))


def modal_system(modal, rfa: RogerRFA, n_ctrl: int = 0, damping: bool = True, **kw) -> AeroelasticSystem:
    """Convenience: build an :class:`AeroelasticSystem` from a ModalModel."""
    C = modal.C() if damping else None
    return AeroelasticSystem(modal.M, modal.K, rfa, C=C, n_ctrl=n_ctrl,
                             names=list(rfa.names), **kw)
