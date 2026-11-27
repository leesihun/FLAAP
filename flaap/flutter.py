"""Flutter solutions and flutter / divergence identification.

Replaces the Fortran modules ``flutteranalysis`` (``flutter_k_method``,
``flutter_kmethod``, ``flutter_V_g_method``, ``vg_method_old``),
``PKflutter`` (``flutter_P-K_method``) and ``saveV_g``.

Methods
-------
``k_method``   classical V-g method (artificial structural damping g)::

        [K (1 + i g) - w^2 (M + rho b^2 Q(k) / (2 k^2))] eta = 0

``pk_method``  British / Nastran PK method (iterated on k per branch)::

        [M p^2 + (C - rho b V Q_I / (2k)) p + (K - q Q_R)] eta = 0,
        p = w (gamma + i),   g = 2 gamma

``p_method``   exact eigenvalues of the time-domain state-space model
               built on a Roger RFA (see :mod:`flaap.timedomain`).

Sign convention: ``g > 0`` is unstable.  Each result carries tracked
branches (``V``, ``freq``, ``g``, ``k`` arrays of shape (npts, nbranch))
and the automatically identified critical points (flutter onset, hump
modes, static divergence).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .aero.gaf import GAFTable
from .linalg import mac


@dataclass
class FlutterPoint:
    kind: str            # 'flutter' | 'divergence' | 'hump-recovery'
    branch: int
    velocity: float
    freq_hz: float
    k: float
    density: float
    g_slope: float = 0.0

    @property
    def omega(self) -> float:
        return 2 * np.pi * self.freq_hz

    def __str__(self):
        return (f"{self.kind:>11s}  mode {self.branch + 1:2d}  V = {self.velocity:9.3f} m/s  "
                f"f = {self.freq_hz:8.3f} Hz ({self.omega:8.3f} rad/s)  k = {self.k:7.4f}")


@dataclass
class FlutterResult:
    method: str
    V: np.ndarray            # (npts, nb)
    freq: np.ndarray         # (npts, nb) [Hz]
    g: np.ndarray            # (npts, nb)
    k: np.ndarray            # (npts, nb)
    density: float
    bref: float
    names: list[str] = field(default_factory=list)
    vectors: np.ndarray | None = None   # (npts, nm, nb) complex
    points: list[FlutterPoint] = field(default_factory=list)
    g_struct: float = 0.0
    k_valid: float = np.inf          # max reduced frequency covered by the GAF table

    @property
    def nbranch(self) -> int:
        return self.V.shape[1]

    def critical(self) -> FlutterPoint | None:
        fl = [p for p in self.points if p.kind in ("flutter", "divergence")]
        return min(fl, key=lambda p: p.velocity) if fl else None

    def summary(self) -> str:
        lines = [f"Flutter analysis ({self.method}), rho = {self.density:.5g} kg/m^3"]
        if not self.points:
            lines.append("  no flutter / divergence found in the analysed range")
        for p in sorted(self.points, key=lambda p: p.velocity):
            lines.append("  " + str(p))
        c = self.critical()
        if c is not None:
            lines.append(f"  => critical: {c.kind} at V = {c.velocity:.3f} m/s, "
                         f"f = {c.freq_hz:.3f} Hz")
        return "\n".join(lines)

    # ----------------------------------------------------------- legacy IO
    def write_legacy(self, directory: str = ".", nmodes: int | None = None):
        """Write ``Vg_save.dat``, ``Vf_save.dat`` and ``pk_results_mode*.dat``
        in the format of the Fortran program (list-directed output)."""
        import os
        nb = self.nbranch if nmodes is None else min(nmodes, self.nbranch)
        with open(os.path.join(directory, "Vg_save.dat"), "w") as f:
            f.write(" For each mode, velocity, damping\n")
            for i in range(nb):
                f.write(f" Mode number : {i + 1:12d}\n")
                for v, g in zip(self.V[:, i], self.g[:, i]):
                    f.write(f"  {v:22.15f}  {g:22.15E}\n")
        with open(os.path.join(directory, "Vf_save.dat"), "w") as f:
            f.write(" For each mode, velocity, frequency\n")
            for i in range(nb):
                f.write(f" Mode number : {i + 1:12d}\n")
                for v, fr in zip(self.V[:, i], self.freq[:, i]):
                    f.write(f"  {v:22.15f}  {fr:22.15E}\n")
        if self.method.startswith("pk") or self.method.startswith("p-"):
            for i in range(nb):
                with open(os.path.join(directory, f"pk_results_mode{i + 1}.dat"), "w") as f:
                    for v, g, fr in zip(self.V[:, i], self.g[:, i], self.freq[:, i]):
                        f.write(f"  {v:20.14f}  {g:22.15E}  {fr:20.12f}\n")   # Hz, as the Fortran program

    def plot(self, ax=None, vunit=1.0, title=None, show_points=True):
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
        for i in range(self.nbranch):
            lab = self.names[i] if i < len(self.names) else f"mode {i + 1}"
            ax[0].plot(self.V[:, i] / vunit, self.g[:, i], lw=1.3, label=lab)
            ax[1].plot(self.V[:, i] / vunit, self.freq[:, i], lw=1.3)
        ax[0].axhline(self.g_struct, color="k", lw=0.8, ls="--")
        if show_points:
            for p in self.points:
                ax[0].plot(p.velocity / vunit, self.g_struct, "ro" if p.kind == "flutter" else "ks")
                ax[1].plot(p.velocity / vunit, p.freq_hz, "ro" if p.kind == "flutter" else "ks")
        ax[0].set_ylabel("damping g")
        ax[1].set_ylabel("frequency [Hz]")
        ax[1].set_xlabel("velocity" + (" [m/s]" if vunit == 1.0 else ""))
        ax[0].legend(fontsize=7, ncol=2)
        ax[0].grid(alpha=0.3)
        ax[1].grid(alpha=0.3)
        ax[0].set_title(title or f"{self.method}  (rho={self.density:.4g})")
        return ax


# ----------------------------------------------------------------------
# identification
# ----------------------------------------------------------------------
def identify(res: FlutterResult, g_struct: float | None = None,
             div_freq_tol: float = 1e-3) -> list[FlutterPoint]:
    """Find g zero crossings (flutter), hump-mode recoveries and divergence."""
    gs = res.g_struct if g_struct is None else g_struct
    pts = []
    for b in range(res.nbranch):
        V, g, f, k = res.V[:, b], res.g[:, b], res.freq[:, b], res.k[:, b]
        order = np.argsort(V)
        V, g, f, k = V[order], g[order], f[order], k[order]
        ok = np.isfinite(V) & np.isfinite(g)
        V, g, f, k = V[ok], g[ok], f[ok], k[ok]
        for i in range(len(V) - 1):
            g0, g1 = g[i] - gs, g[i + 1] - gs
            if g0 < 0 <= g1 or g0 > 0 >= g1:
                t = g0 / (g0 - g1) if g0 != g1 else 0.0
                Vc = V[i] + t * (V[i + 1] - V[i])
                fc = f[i] + t * (f[i + 1] - f[i])
                kc = k[i] + t * (k[i + 1] - k[i])
                if kc > res.k_valid * 1.001:
                    continue      # outside the tabulated k range (extrapolated GAFs)
                up = g1 > g0
                fmax = max(np.nanmax(res.freq), 1e-12)
                if up:
                    kind = "divergence" if fc < div_freq_tol * fmax else "flutter"
                else:
                    kind = "hump-recovery"
                slope = (g1 - g0) / max(V[i + 1] - V[i], 1e-30)
                pts.append(FlutterPoint(kind, b, float(Vc), float(fc), float(kc), res.density, slope))
    res.points = pts
    return pts


# ----------------------------------------------------------------------
# K (V-g) method
# ----------------------------------------------------------------------
def k_method(M, K, table: GAFTable, density: float, kvals=None, track: bool = True,
             g_struct: float = 0.0) -> FlutterResult:
    """V-g method.  ``kvals`` defaults to a fine grid inside the table range
    (descending ``k`` == increasing speed)."""
    M = np.asarray(M, float)
    K = np.asarray(K, float)
    b = table.bref
    if kvals is None:
        kpos = table.k[table.k > 0]
        kvals = np.geomspace(kpos.max(), max(kpos.min(), 1e-3), 200)
    kvals = np.asarray(kvals, float)
    Kinv = np.linalg.inv(K)
    nm = M.shape[0]
    V = np.zeros((kvals.size, nm))
    F = np.zeros_like(V)
    G = np.zeros_like(V)
    Kk = np.zeros_like(V)
    vecs = np.zeros((kvals.size, nm, nm), complex)
    prev = None
    for n, k in enumerate(kvals):
        Q, _ = table.interp(k)
        A = Kinv @ (M + density * b * b * Q / (2 * k * k))
        lam, vec = np.linalg.eig(A)
        if track and prev is not None:
            Mc = mac(prev, vec)
            perm = _assign(Mc)
            lam, vec = lam[perm], vec[:, perm]
        elif prev is None:
            order = np.argsort(-lam.real)  # largest Re(lambda) -> lowest frequency
            lam, vec = lam[order], vec[:, order]
        prev = vec
        with np.errstate(invalid="ignore", divide="ignore"):
            w = 1.0 / np.sqrt(lam.real)
        g = lam.imag / lam.real
        V[n] = w * b / k
        F[n] = w / (2 * np.pi)
        G[n] = g
        Kk[n] = k
        vecs[n] = vec
    res = FlutterResult("k-method", V, F, G, Kk, density, b, list(table.names), vecs,
                        g_struct=g_struct, k_valid=float(table.k.max()))
    identify(res)
    return res


def _assign(Mc: np.ndarray) -> np.ndarray:
    """Column permutation maximizing MAC with previous vectors (Hungarian)."""
    from scipy.optimize import linear_sum_assignment
    r, c = linear_sum_assignment(-Mc)
    perm = np.empty_like(c)
    perm[r] = c
    return perm


# ----------------------------------------------------------------------
# PK method
# ----------------------------------------------------------------------
def _pk_state(M, C, K, Q, q, rho, V, b, k):
    QR, QI = Q.real, Q.imag
    Ceff = C - (rho * V * b / (2 * k)) * QI
    Keff = K - q * QR
    nm = M.shape[0]
    Minv = np.linalg.inv(M)
    A = np.zeros((2 * nm, 2 * nm))
    A[:nm, nm:] = np.eye(nm)
    A[nm:, :nm] = -Minv @ Keff
    A[nm:, nm:] = -Minv @ Ceff
    return A


def pk_method(M, K, table: GAFTable, density: float, velocities, C=None,
              tol: float = 1e-6, maxit: int = 60, kmin: float = 1e-4,
              g_struct: float = 0.0, nbranch: int | None = None) -> FlutterResult:
    """PK method with per-branch k iteration and continuation tracking."""
    M = np.asarray(M, float)
    K = np.asarray(K, float)
    nm = M.shape[0]
    C = np.zeros((nm, nm)) if C is None else np.asarray(C, float)
    b = table.bref
    velocities = np.asarray(velocities, float)
    nb = nm if nbranch is None else nbranch
    # initial guesses: in-vacuo roots
    lam0, vec0 = np.linalg.eig(np.linalg.solve(M, K))
    order = np.argsort(lam0.real)[:nb]
    w0 = np.sqrt(np.maximum(lam0.real[order], 0))
    v_prev = [vec0[:, i].astype(complex) for i in order]
    p_prev = [1j * w for w in w0]
    out_V = np.zeros((velocities.size, nb))
    out_F = np.zeros_like(out_V)
    out_G = np.zeros_like(out_V)
    out_K = np.zeros_like(out_V)
    out_vec = np.zeros((velocities.size, nm, nb), complex)
    hist = [[p] for p in p_prev]
    for n, V in enumerate(velocities):
        q = 0.5 * density * V * V
        taken = []
        for j in range(nb):
            # predictor (linear extrapolation in V)
            if len(hist[j]) >= 2 and n >= 2:
                dv = velocities[n - 1] - velocities[n - 2]
                pg = hist[j][-1] + (hist[j][-1] - hist[j][-2]) * (V - velocities[n - 1]) / dv
            else:
                pg = hist[j][-1]
            k = max(abs(pg.imag) * b / V, kmin)
            pj, vj = pg, None
            for it in range(maxit):
                Q, _ = table.interp(k)
                A = _pk_state(M, C, K, Q, q, density, V, b, k)
                lam, vec = np.linalg.eig(A)
                cand = np.where(lam.imag >= -1e-9 * np.abs(lam))[0]
                # roots already claimed by another branch at this speed are
                # excluded (prevents two branches collapsing onto one root)
                free = [c for c in cand if all(abs(lam[c] - t) > 2e-3 * max(abs(t), 1e-6) for t in taken)]
                cand = np.array(free) if free else cand
                # eigenvalue proximity + mode-shape correlation (MAC) with the
                # branch's previous vector: keeps close modes apart
                d = np.abs(lam[cand] - pg) / max(abs(pg), 1e-6)
                V_ = vec[:nm, cand]
                vp = v_prev[j]
                mac = np.abs(vp.conj() @ V_) ** 2 / (np.real(vp.conj() @ vp) *
                                                     np.sum(np.abs(V_) ** 2, axis=0) + 1e-300)
                i = cand[np.argmin(d - 0.2 * mac)]
                pj, vj = lam[i], vec[:nm, i]
                k_new = max(abs(pj.imag) * b / V, kmin)
                if abs(k_new - k) < tol * max(k, 1e-3):
                    k = k_new
                    break
                k = 0.5 * (k + k_new) if it > 20 else k_new
            hist[j].append(pj)
            taken.append(pj)
            if vj is not None:
                v_prev[j] = vj
            w = abs(pj.imag)
            out_V[n, j] = V
            out_F[n, j] = w / (2 * np.pi)
            out_G[n, j] = 2 * pj.real / w if w > 1e-9 * (abs(pj.real) + 1) else np.sign(pj.real) * 1.0
            out_K[n, j] = w * b / V
            out_vec[n, :, j] = vj
    res = FlutterResult("pk-method", out_V, out_F, out_G, out_K, density, b,
                        list(table.names[:nb]), out_vec, g_struct=g_struct,
                        k_valid=float(table.k.max()))
    identify(res)
    # aperiodic root crossing the origin: divergence (k -> 0, real root > 0)
    for p in res.points:
        if p.freq_hz < 1e-3 * max(np.max(res.freq), 1e-9) and p.kind == "flutter":
            p.kind = "divergence"
    return res


# ----------------------------------------------------------------------
# p method on a state-space (RFA) model
# ----------------------------------------------------------------------
def p_method(build_A, velocities, density: float, nstruct: int, bref: float,
             names=None, w_initial=None, k_valid: float = np.inf) -> FlutterResult:
    """Eigenvalues of ``build_A(V)`` (continuous state matrix) tracked over V.

    Only the ``nstruct`` structural branches are tracked (by continuity
    from the in-vacuo roots); aerodynamic lag roots are ignored unless
    they become unstable, in which case they are reported as extra points.
    """
    velocities = np.asarray(velocities, float)
    prev = None
    out = {k: np.zeros((velocities.size, nstruct)) for k in ("V", "F", "G", "K")}
    for n, V in enumerate(velocities):
        A = build_A(V)
        lam = np.linalg.eigvals(A)
        lam = lam[lam.imag >= -1e-9]
        if prev is None:
            if w_initial is not None:
                sel = [lam[np.argmin(np.abs(lam - 1j * w))] for w in w_initial]
            else:
                osc = lam[lam.imag > 1e-6]
                sel = list(osc[np.argsort(osc.imag)][:nstruct])
            prev = np.array(sel)
        else:
            new = []
            used = set()
            for p0 in prev:
                d = np.abs(lam - p0)
                for u in used:
                    d[u] = np.inf
                i = int(np.argmin(d))
                used.add(i)
                new.append(lam[i])
            prev = np.array(new)
        w = np.abs(prev.imag)
        out["V"][n] = V
        out["F"][n] = w / (2 * np.pi)
        with np.errstate(divide="ignore", invalid="ignore"):
            out["G"][n] = np.where(w > 1e-9, 2 * prev.real / w, np.sign(prev.real))
        out["K"][n] = w * bref / max(V, 1e-12)
    res = FlutterResult("p-method (RFA state space)", out["V"], out["F"], out["G"], out["K"],
                        density, bref, list(names or []))
    res.k_valid = k_valid
    identify(res)
    return res


# ----------------------------------------------------------------------
# convenience wrappers
# ----------------------------------------------------------------------
def structural_matrices(modal, nmodes=None, include_damping=True):
    mm = modal if nmodes is None else modal.select(range(nmodes))
    C = mm.C() if include_damping else np.zeros((mm.nmodes, mm.nmodes))
    return mm.M, mm.K, C


def matched_point(solve_at_density, mach: float, h_lo: float = -1000.0, h_hi: float = 20000.0,
                  tol: float = 1.0):
    """Altitude where the flutter speed equals ``M * a(h)`` (matched flutter
    point at fixed Mach).  ``solve_at_density(rho) -> flutter speed``."""
    from .atmosphere import isa

    def f(h):
        a = isa(h)
        vf = solve_at_density(a.density)
        return (vf if vf is not None else np.inf) - mach * a.speed_of_sound

    flo, fhi = f(h_lo), f(h_hi)
    if np.sign(flo) == np.sign(fhi):
        return None
    for _ in range(60):
        hm = 0.5 * (h_lo + h_hi)
        fm = f(hm)
        if np.sign(fm) == np.sign(flo):
            h_lo, flo = hm, fm
        else:
            h_hi, fhi = hm, fm
        if h_hi - h_lo < tol:
            break
    return 0.5 * (h_lo + h_hi)


def flutter_margin_check(vf_eas: float, vd_eas: float, factor: float = 1.15) -> dict:
    """CS/FAR 25.629: no flutter up to 1.15 V_D (EAS)."""
    return {"V_flutter": vf_eas, "required": factor * vd_eas,
            "margin": vf_eas / vd_eas - 1.0, "compliant": vf_eas >= factor * vd_eas}
