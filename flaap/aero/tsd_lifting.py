"""Lifting (non-symmetric) TSD mean flow: angle of attack and camber.

:class:`flaap.aero.tsd.TSDSolver` solves the nonlinear steady flow only for
the symmetric thickness problem (zero incidence); lift enters linearly.
:class:`LiftingTSDSolver` removes that restriction.  The steady mean flow
about a wing at angle of attack ``alpha`` with thickness *and* camber is
solved with the full TSD nonlinearity on both surfaces, so the upper- and
lower-surface shocks move with incidence.  The time-linearised flutter
problem is then solved about that non-symmetric mean flow.

Formulation (upper / lower fields on one half grid)
---------------------------------------------------
The domain ``z >= 0`` carries the upper potential ``U(x, y, z) = phi(x, y, z)``.
The mirrored lower potential ``L(x, y, zeta) = phi(x, y, -zeta)`` lives on the
same half grid.  The TSD equation is invariant under ``z -> -z``, so ``U``
and ``L`` each satisfy it with their own coefficient
``A = 1 - M^2 - (g+1) M^2 phi_x`` and their own Murman switches.  They are
coupled only on ``z = 0``:

* **Wing** (thin-wing boundary conditions on the mean plane):

  - upper: ``U_z = dz_u/dx = t'/2 + z_c' - alpha``
  - lower: ``L_zeta = -dz_l/dx = t'/2 - z_c' + alpha``

* **Off the wing, upstream and outboard:** ``U = L`` (continuity), and the
  upper and lower half-cell flux balances are added, so the normal flux is
  continuous.
* **Wake:** normal flux is continuous (summed half cells).  The jump
  ``U - L`` is convected from the trailing edge:
  ``(U-L)(x) = (U-L)_TE exp(-i kbar (x - x_TE))``.  This keeps the pressure
  continuous and carries the Kutta condition (``kbar = 0`` in the steady
  problem: constant circulation jump).

Pressures: ``Cp_u = -2 U_x``, ``Cp_l = -2 L_x``.  Jump:
``dCp = Cp_l - Cp_u = 2 (U - L)_x``, and in the unsteady problem
``2 [(U1 - L1)_x + i kbar (U1 - L1)]``.  The antisymmetric part
``(U1 - L1)/2`` therefore plays the role of the potential of the symmetric
solver.  The generalised forces, ``dCp`` and the trailing-edge potential of
:class:`TSDSolver` are reused unchanged.

Checks:

* With ``alpha = 0`` and no camber, ``U = L``: the steady flow and the GAFs
  equal those of :class:`TSDSolver` to round-off.
* At small incidence, ``CL / alpha`` tends to the linear (``k = 0``)
  lifting solution.
* ``alpha -> -alpha`` swaps the upper and lower surfaces.

Unknowns are ordered plane by plane, ``[U_i, L_i]``, so the block-banded
elimination of :mod:`flaap.aero.tsd` still applies, with planes of size
``2 ny nz``: about 8 times the dense work and 4 times the memory of the
symmetric solver.

Limitations (in addition to those of :mod:`flaap.aero.tsd`):

* inviscid, so shocks are usually stronger and further aft than in
  experiments, especially for thick sections or large incidence;
* thin-wing boundary conditions on ``z = 0``;
* no twist;
* the shock is frozen in the time-linearised problem.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .tsd import SteadyResult, TSDSolver, _block_solve, cp_critical, naca4_thickness

__all__ = ["LiftingTSDSolver", "LiftingSteadyResult", "naca4_camber", "naca4_section",
           "shock_positions"]


# ---------------------------------------------------------------- sections
def naca4_camber(m: float, p: float) -> Callable:
    """NACA 4-digit mean line ``z_c/c`` (max camber ``m`` at ``p``) with ``.slope``."""
    m, p = float(m), float(p)

    def zc(xi):
        xi = np.asarray(xi, float)
        if m == 0.0:
            return np.zeros_like(xi)
        return np.where(xi < p, m / p ** 2 * (2 * p * xi - xi ** 2),
                        m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * xi - xi ** 2))

    def slope(xi):
        xi = np.asarray(xi, float)
        if m == 0.0:
            return np.zeros_like(xi)
        return np.where(xi < p, 2 * m / p ** 2 * (p - xi), 2 * m / (1 - p) ** 2 * (p - xi))

    zc.slope = slope
    zc.m, zc.p = m, p
    return zc


def naca4_section(code: str):
    """``'2412'`` -> (thickness callable, camber callable or None)."""
    code = str(code).lower().removeprefix("naca")
    if len(code) != 4 or not code.isdigit():
        raise ValueError(f"not a NACA 4-digit code: {code!r}")
    m, p, tc = int(code[0]) / 100.0, int(code[1]) / 10.0, int(code[2:]) / 100.0
    if m > 0 and p == 0:
        raise ValueError(f"NACA {code}: camber without a camber position")
    return naca4_thickness(tc), (naca4_camber(m, p) if m > 0 else None)


def shock_positions(xi, cp, cp_star, min_jump: float = 0.05):
    """Shock chord fraction and Cp jump per strip.

    A shock is a recompression from ``Cp < Cp*`` to ``Cp > Cp*`` along the chord
    (the Murman shock-point operator spreads it over 1-3 nodes).  The position
    is the midpoint of the steepest node interval of the strongest such
    recompression.  ``xi, cp``: (nspan, nchord).  Returns ``(x_s, dcp)``, each
    (nspan,), with NaN / 0 where there is no shock.
    """
    xi = np.asarray(xi, float)
    cp = np.asarray(cp, float)
    ns, nc = cp.shape
    xs = np.full(ns, np.nan)
    jump = np.zeros(ns)
    for j in range(ns):
        c = cp[j]
        best = (0.0, None)
        for i in np.nonzero((c[:-1] < cp_star) & (c[1:] >= cp_star))[0]:
            lo, hi = max(i - 1, 0), min(i + 3, nc - 1)
            dj = float(c[hi] - c[lo])
            if dj > best[0]:
                k = lo + int(np.argmax(np.diff(c[lo:hi + 1])))
                best = (dj, 0.5 * (xi[j, k] + xi[j, k + 1]))
        if best[1] is not None and best[0] >= min_jump:
            jump[j], xs[j] = best
    return xs, jump


# ---------------------------------------------------------------- results
@dataclass
class LiftingSteadyResult(SteadyResult):
    phi_lower: np.ndarray = None       # (NU, NY, NZ) mirrored lower potential L
    cp_lower: np.ndarray = None        # (nspan, nchord) lower-surface Cp = -2 L_x
    xi: np.ndarray = None              # (nspan, nchord) chord fraction of the wing nodes
    alpha_deg: float = 0.0
    cl: np.ndarray = None              # (nspan,) section lift coefficient
    CL: float = np.nan                 # half-wing lift coefficient (planform area)
    shock_upper: np.ndarray = None     # (nspan,) x/c, NaN = none
    shock_lower: np.ndarray = None
    jump_upper: np.ndarray = None      # (nspan,) Cp rise across the shock
    jump_lower: np.ndarray = None
    n_supersonic_lower: int = 0

    def table(self) -> str:
        """Per-strip table: y, cl, upper / lower shock x/c."""
        lines = [f"TSD lifting mean flow: M = {self.mach}, alpha = {self.alpha_deg:g} deg, "
                 f"CL = {self.CL:.4f}, Cp* = {self.cp_star:.4f}",
                 "      y        cl     x_s/c upper  dCp_s   x_s/c lower  dCp_s"]
        for j in range(self.cl.size):
            lines.append(f"  {self.y[j]:8.4f}  {self.cl[j]:8.4f}   {self.shock_upper[j]:9.3f}  "
                         f"{self.jump_upper[j]:6.3f}   {self.shock_lower[j]:9.3f}  {self.jump_lower[j]:6.3f}")
        return "\n".join(lines)


# ---------------------------------------------------------------- solver
class LiftingTSDSolver(TSDSolver):
    """TSD solver with a lifting (non-symmetric) nonlinear mean flow.

    Same parameters as :class:`flaap.aero.tsd.TSDSolver`, plus:

    alpha_deg : angle of attack [deg] (uniform; positive nose up)
    camber : callable ``z_c(xi)`` (mean line / chord, optional ``.slope``) or None
    """

    lifting = True

    def __init__(self, surface, mach: float, thickness: Callable | None = None, *,
                 alpha_deg: float = 0.0, camber: Callable | None = None, **kw):
        super().__init__(surface, mach, thickness, **kw)
        self.alpha_deg = float(alpha_deg)
        self.alpha = np.radians(self.alpha_deg)
        self.camber = camber
        self._masks_cache = {}

    # -------------------------------------------------------------- helpers
    def _camber_slope(self, xi):
        c = self.camber
        if c is None:
            return np.zeros_like(xi)
        if hasattr(c, "slope"):
            return np.asarray(c.slope(xi), float)
        h = 1e-6
        return 0.5 * (np.asarray(c(xi + h)) - np.asarray(c(xi - h))) / h

    def surface_slopes(self):
        """(dz_u/dx, dz_l/dx) at the wing nodes, (nspan, nchord)."""
        xi = self.wing_xi
        t = self._thickness_slope(xi)
        c = self._camber_slope(xi)
        return t + c - self.alpha, -t + c - self.alpha

    @property
    def trivial(self) -> bool:
        """No thickness, camber or incidence: the mean flow is zero."""
        return self.thickness is None and self.camber is None and self.alpha == 0.0

    def _masks(self, fv):
        key = fv.tobytes()
        if key not in self._masks_cache:
            I, J, K = self._I, self._J, self._K
            z0 = K == 0
            wing0 = z0 & self.wing2d[I, J]
            wake0 = z0 & self.wake2d[I, J]
            lcont = z0 & ~wing0                      # L rows: continuity / wake jump
            m = SimpleMasks(
                mU=fv.astype(float),                 # U rows: own FV balance ...
                msum=(fv & lcont).astype(float),     # ... plus the lower half cell
                mL=(fv & ~lcont).astype(float),      # L rows: own FV balance
                bL=(~fv & ~lcont).astype(float),     # L rows: far-field condition
                lcont=lcont, wake0=wake0)
            self._masks_cache = {key: m}
        return self._masks_cache[key]

    def _continuity(self, kbar: float, masks):
        """Rows ``C`` (N x N) acting on the jump ``U - L`` on the z = 0 plane
        off the wing: ``U - L = 0``.  Wake: ``J_i - J_{i-1} exp(-i kbar dx) = 0``."""
        sc = self.bc_scale
        u = self._pml_coords()[0] if (self.pml and kbar > 0) else self.u
        I, J = self._I, self._J
        off = np.nonzero(masks.lcont & ~masks.wake0)[0]
        wk = np.nonzero(masks.wake0)[0]
        jn = J[wk]
        ph = np.exp(-1j * kbar * self.b[jn] * (u[I[wk]] - u[I[wk] - 1]))
        rows = np.r_[off, wk, wk]
        cols = np.r_[off, wk, self._node(I[wk] - 1, jn, 0)]
        vals = np.r_[np.full(off.size + wk.size, sc, complex), -sc * ph]
        C = sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N))
        return C if kbar > 0 else C.real.tocsr()

    def _solve2(self, Kmat, rhs):
        """Solve the coupled (U, L) system, unknowns [U; L]."""
        N = self.N
        if self.linear_solver == "splu":
            lu = spla.splu(Kmat.tocsc())
            return lu.solve(np.asarray(rhs))
        NU = self.shape[0]
        P = self.shape[1] * self.shape[2]
        base = np.arange(P)
        perm = np.concatenate([np.r_[i * P + base, N + i * P + base] for i in range(NU)])
        Kp = Kmat.tocsr()[perm][:, perm]
        xp = _block_solve(Kp, np.asarray(rhs)[perm], 2 * P, NU, self.max_bytes)
        x = np.empty_like(xp)
        x[perm] = xp
        return x

    def _assemble(self, KU, KL, B, C, masks):
        """Coupled operator from the single-field FV operators ``KU``, ``KL``."""
        dg = sp.diags
        JUU = dg(masks.mU) @ KU + B
        JUL = dg(masks.msum) @ KL
        JLL = dg(masks.mL) @ KL + dg(masks.bL) @ B - C
        return sp.bmat([[JUU, JUL], [C, JLL]], format="csr")

    # ---------------------------------------------------------------- steady
    def set_alpha(self, alpha_deg: float) -> None:
        """Change the incidence, keeping the operators and the last solution
        (used as the initial guess by ``steady(warm_start=True)``)."""
        self.alpha_deg = float(alpha_deg)
        self.alpha = np.radians(self.alpha_deg)
        self._mean = None

    def steady(self, tol: float = 1e-9, maxiter: int = 40, verbose: bool | None = None,
               freeze_after: int = 25, warm_start: bool = False, ptc: bool = True,
               ptc_maxiter: int = 150, ptc_cfl0: float = 100.0) -> LiftingSteadyResult:
        """Nonlinear lifting steady flow (thickness + camber + incidence).

        Newton's method with sonic-switch updates and a backtracking line
        search.  ``warm_start`` starts from the previous solution of this
        solver (incidence sweeps).  If that does not converge within
        ``maxiter`` iterations or the line search stalls (four steps with
        lambda < 0.1: strong shocks that must travel many cells, the switches
        chatter), the solve is restarted from rest with pseudo-transient
        continuation (always from rest: started from a converged sharp-shock
        state of another incidence the first pseudo-time steps blow up): a local
        pseudo-time term ``J_ii / CFL`` on the flow-equation rows, CFL grown
        by switched evolution relaxation (CFL ~ R0 / R), switches updated
        every step and frozen only once the residual is below 1e-4.  Early
        steps behave like implicit time marching (the shock moves to its
        place), late ones like Newton (quadratic convergence)."""
        verbose = self.verbose if verbose is None else verbose
        N = self.N
        B, fv = self._bc_rows(0.0, antisym=False)
        B = B.real.tocsr()
        mk = self._masks(fv)
        C = self._continuity(0.0, mk)
        su, sl = self.surface_slopes()
        widx = self.wing_idx.ravel()
        rhsU = np.zeros(N)
        rhsL = np.zeros(N)
        rhsU[widx] = self.Vfac_z[widx] * su.ravel()
        rhsL[widx] = -self.Vfac_z[widx] * sl.ravel()
        interior = self.interior

        def field(phi, mu, rhs):
            p = (self.Du @ phi) / self.b_uf
            F, Af = self._flux(p)
            Dsw = self._switched_div(mu)
            return self.Vfac_u * (Dsw @ F) + self.Llin @ phi - rhs, Dsw, Af

        def switches(phi):
            return ((self._node_A(phi) < 0) & interior).astype(float)

        def parts(X, muU, muL):
            U, L = X[:N], X[N:]
            FU, DU, AU = field(U, muU, rhsU)
            FL, DL, AL = field(L, muL, rhsL)
            RU = mk.mU * FU + mk.msum * FL + B @ U
            RL = mk.mL * FL + mk.bL * (B @ L) + C @ (U - L)
            return np.r_[RU, RL], (DU, AU, DL, AL)

        res0 = max(np.linalg.norm(np.r_[rhsU, rhsL]), 1e-300)
        t0 = time.perf_counter()

        def newton(X, maxit):
            muU, muL = switches(X[:N]), switches(X[N:])
            R, lin = parts(X, muU, muL)
            hist = [np.linalg.norm(R) / res0]
            converged = hist[-1] < tol
            it, small = 0, 0
            while not converged and it < maxit:
                if ptc and small >= 4:          # line search stalled: hand over to PTC
                    break
                it += 1
                DU, AU, DL, AL = lin
                KU = sp.diags(self.Vfac_u) @ DU @ sp.diags(AU / self.b_uf) @ self.Du + self.Llin
                KL = sp.diags(self.Vfac_u) @ DL @ sp.diags(AL / self.b_uf) @ self.Du + self.Llin
                Jm = self._assemble(KU, KL, B, C, mk)
                dX = self._solve2(Jm, -R)
                lam, rn_old = 1.0, np.linalg.norm(R)
                for _ in range(8):
                    trial = X + lam * dX
                    if it > freeze_after:
                        mU_t, mL_t = muU, muL
                    else:
                        mU_t, mL_t = switches(trial[:N]), switches(trial[N:])
                    Rt, lin_t = parts(trial, mU_t, mL_t)
                    if np.linalg.norm(Rt) < rn_old or lam < 0.05:
                        break
                    lam *= 0.5
                X, muU, muL, R, lin = trial, mU_t, mL_t, Rt, lin_t
                small = small + 1 if lam < 0.1 else 0
                hist.append(np.linalg.norm(R) / res0)
                if verbose:
                    print(f"  TSD lifting it {it:2d}: |R|/|R0| = {hist[-1]:.3e}  lam = {lam:.3f}  "
                          f"supersonic upper/lower = {int(muU.sum())}/{int(muL.sum())}  "
                          f"({time.perf_counter() - t0:.1f} s)")
                converged = hist[-1] < tol or (np.linalg.norm(lam * dX) < tol * max(1.0, np.abs(X).max())
                                               and hist[-1] < 1e-5)
            return X, muU, muL, it, hist, converged

        rowmask = np.r_[mk.mU, mk.mL] > 0

        def pseudo_transient(X, maxit):
            muU, muL = switches(X[:N]), switches(X[N:])
            R, lin = parts(X, muU, muL)
            hist = [np.linalg.norm(R) / res0]
            cfl, it, frozen, converged = ptc_cfl0, 0, False, hist[-1] < tol
            while not converged and it < maxit:
                it += 1
                DU, AU, DL, AL = lin
                KU = sp.diags(self.Vfac_u) @ DU @ sp.diags(AU / self.b_uf) @ self.Du + self.Llin
                KL = sp.diags(self.Vfac_u) @ DL @ sp.diags(AL / self.b_uf) @ self.Du + self.Llin
                Jm = self._assemble(KU, KL, B, C, mk)
                # pseudo-time term with the sign of the (Laplacian-like,
                # negative) diagonal: J + diag(J)/CFL damps the Newton step
                dJ = Jm.diagonal()
                sgn = np.sign(np.median(dJ[rowmask])) or -1.0
                floor = 1e-8 * np.abs(dJ[rowmask]).max()
                dJ = np.where(rowmask, np.where(sgn * dJ > floor, dJ, sgn * floor), 0.0)
                rn_old = np.linalg.norm(R)
                for _ in range(6):
                    dX = self._solve2(Jm + sp.diags(dJ / cfl), -R)
                    trial = X + dX
                    mU_t, mL_t = (muU, muL) if frozen else (switches(trial[:N]), switches(trial[N:]))
                    Rt, lin_t = parts(trial, mU_t, mL_t)
                    rn = np.linalg.norm(Rt)
                    if rn < 3.0 * rn_old or cfl <= 1e-3:
                        break
                    cfl *= 0.2
                X, muU, muL, R, lin = trial, mU_t, mL_t, Rt, lin_t
                hist.append(rn / res0)
                cfl = float(np.clip(cfl * np.clip(rn_old / max(rn, 1e-300), 0.3, 4.0), 1e-3, 1e12))
                frozen = frozen or hist[-1] < 1e-4
                if verbose:
                    print(f"  TSD lifting PTC it {it:3d}: |R|/|R0| = {hist[-1]:.3e}  CFL = {cfl:.2e}  "
                          f"supersonic upper/lower = {int(muU.sum())}/{int(muL.sum())}  "
                          f"({time.perf_counter() - t0:.1f} s)")
                converged = hist[-1] < tol or (np.linalg.norm(dX) < tol * max(1.0, np.abs(X).max())
                                               and hist[-1] < 1e-5)
            return X, muU, muL, it, hist, converged

        X0 = np.zeros(2 * N)
        if warm_start and getattr(self, "_X_last", None) is not None:
            X0 = self._X_last.copy()
        X, muU, muL, it, hist, converged = newton(X0, maxiter)
        if not converged and ptc:
            if verbose:
                print(f"  TSD lifting: Newton not converged ({hist[-1]:.1e}); pseudo-transient continuation")
            X, muU, muL, itp, h, converged = pseudo_transient(np.zeros(2 * N), ptc_maxiter)
            it += itp
            hist += h
        self._X_last = X.copy()
        res = self._result(X[:N], X[N:], int(muU.sum()), int(muL.sum()), bool(converged), it, hist)
        self._mean = res
        return res

    def _result(self, U, L, nsu, nsl, converged, it, hist) -> LiftingSteadyResult:
        ns = self.nspan
        sl = (slice(self.iw0, self.iw1 + 1), slice(0, ns), 0)
        pxU = ((self.Cu @ U) / self.b_n).reshape(self.shape)
        pxL = ((self.Cu @ L) / self.b_n).reshape(self.shape)
        cpu = -2 * pxU[sl].T
        cpl = -2 * pxL[sl].T
        # section lift: int dCp dx = 2 Lref (U - L)_TE  (dCp = 2 (U - L)_x)
        jump = (U - L).reshape(self.shape)[self.iw1, :ns, 0]
        chord = self.b[:ns] / self.gmap * self.Lref
        cl = 2 * self.Lref * jump / chord
        w = self.strip_width * chord
        CL = float(np.sum(w * cl) / np.sum(w))
        cps = cp_critical(self.mach, self.gamma, tsd=True) if self.mach > 0 else np.nan
        xi = self.wing_xi
        su, ju = shock_positions(xi, cpu, cps)
        sw, jl = shock_positions(xi, cpl, cps)
        return LiftingSteadyResult(
            U.reshape(self.shape), cpu, self.wing_x.copy(), self.wing_y.copy(), self.mach,
            nsu, converged, it, hist, cps,
            phi_lower=L.reshape(self.shape), cp_lower=cpl, xi=xi.copy(), alpha_deg=self.alpha_deg,
            cl=cl, CL=CL, shock_upper=su, shock_lower=sw, jump_upper=ju, jump_lower=jl,
            n_supersonic_lower=nsl)

    @property
    def mean_flow(self) -> LiftingSteadyResult:
        if self._mean is None:
            if self.trivial:
                z = np.zeros(self.N)
                self._mean = self._result(z, z, 0, 0, True, 0, [])
            else:
                self.steady()
        return self._mean

    # -------------------------------------------------------------- unsteady
    def unsteady_matrix(self, kbar: float):
        """Coupled time-linearised (U1, L1) operator about the lifting mean flow;
        returns (K, fv) with ``K`` of size 2N."""
        st = self.mean_flow
        M2 = self.mach ** 2
        B, fv = self._bc_rows(kbar, antisym=False)
        if kbar == 0:
            B = B.real.tocsr()
        mk = self._masks(fv)
        C = self._continuity(kbar, mk)
        o = self._unsteady_ops() if kbar > 0 else self

        def op(phi0):
            phi0 = phi0.ravel()
            p = (self.Du @ phi0) / self.b_uf
            _, Af = self._flux(p)
            mu = ((self._node_A(phi0) < 0) & self.interior).astype(float)
            Dsw = self._switched_div(mu)
            X = sp.diags(o.Vfac_u) @ Dsw @ (sp.diags(Af / o.b_uf) @ o.Du - 2j * kbar * M2 * o.Au)
            return X + o.Llin + sp.diags(M2 * kbar ** 2 * o.vol_b)

        K = self._assemble(op(st.phi), op(st.phi_lower), B, C, mk)
        K.eliminate_zeros()
        return K, fv

    def solve_normalwash(self, kbar: float, w: np.ndarray) -> np.ndarray:
        """Solve for the wing normalwash ``w/U`` (nspan*nchord, nrhs).

        Both surfaces move with the thin wing: ``U1_z = w``, ``L1_zeta = -w``.
        Returns the antisymmetric part ``(U1 - L1)/2``, which has the role of
        ``phi1`` in :class:`TSDSolver` (N, nrhs).
        """
        N = self.N
        Kmat, fv = self.unsteady_matrix(kbar)
        mk = self._masks(fv)
        w = np.asarray(w, complex).reshape(self.nspan * self.nchord, -1)
        idx = self.wing_idx.ravel()
        rhs = np.zeros((2 * N, w.shape[1]), complex)
        rhs[idx] = mk.mU[idx, None] * self.Vfac_z[idx, None] * w
        rhs[N + idx] = -mk.mL[idx, None] * self.Vfac_z[idx, None] * w
        t0 = time.perf_counter()
        X = self._solve2(Kmat, rhs)
        if self.verbose:
            print(f"  TSD lifting unsteady solve kbar={kbar:.4f}: {time.perf_counter() - t0:.2f} s "
                  f"({w.shape[1]} rhs)")
        return 0.5 * (X[:N] - X[N:])


@dataclass
class SimpleMasks:
    mU: np.ndarray
    msum: np.ndarray
    mL: np.ndarray
    bL: np.ndarray
    lcont: np.ndarray
    wake0: np.ndarray = field(default=None)
