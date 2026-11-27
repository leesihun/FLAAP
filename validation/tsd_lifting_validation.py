"""Validation: lifting (incidence / camber) nonlinear TSD mean flow.

1. Consistency checks (exact):
   - alpha = 0, no camber: identical to the symmetric solver (flow and GAFs);
   - alpha -> -alpha swaps the surfaces;
   - small incidence: CL/alpha equals the linear (k = 0) TSD lift and is close
     to the Helmbold lifting-line estimate.
2. NACA 0012 at M = 0.8, incidence sweep on a rectangular wing of aspect
   ratio 16 (root strip ~ 2-D).  Reference at alpha = 1.25 deg: the inviscid
   Euler consensus of AGARD-AR-211 (1985), CL ~ 0.35, upper shock at
   x/c ~ 0.62, weak lower shock at x/c ~ 0.35.  Classical (conservative)
   TSD without the Ames modified terms is known to put the shock further
   aft and give more lift than Euler.  The 3-D root strip also sees the
   downwash of the finite wing (effective incidence ~ CL/(pi A) lower), so
   the comparison is made both at equal alpha and at equal section cl.
3. AGARD 445.6 at M = 0.954: flutter speed index vs mean incidence (no
   experiment at incidence; trend only, the alpha = 0 point is Yates 1987).

Run:  python validation/tsd_lifting_validation.py [--skip-agard] [--alphas-agard 0 1 2] [--plot out.png]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from flaap.aero.tsd import (TSDSolver, modes_from_structure, naca4_thickness, naca65a004,  # noqa: E402
                            rigid_mode_funcs, tsd_gaf_table)
from flaap.aero.tsd_lifting import LiftingTSDSolver  # noqa: E402
from flaap.geometry import Surface  # noqa: E402

AR211 = {"alpha": 1.25, "CL": 0.35, "xs_upper": 0.62, "xs_lower": 0.35}


def consistency(verbose=True):
    surf = Surface("w", (0, 0, 0), 1.0, (0, 4.0, 0), 1.0, nchord=10, nspan=10)
    g = dict(nchord=16, nspan=8, nz=18, n_up=10, n_down=10, n_out=8)
    th = naca4_thickness(0.12)
    s0, s1 = TSDSolver(surf, 0.8, th, **g), LiftingTSDSolver(surf, 0.8, th, **g)
    r0, r1 = s0.steady(), s1.steady()
    f, n = rigid_mode_funcs(0.25)
    qa = tsd_gaf_table(surf, 0.8, [0.0, 0.3], 0.5, f, n, thickness=th, solver=s0).Q
    qb = tsd_gaf_table(surf, 0.8, [0.0, 0.3], 0.5, f, n, thickness=th, solver=s1).Q
    rp = LiftingTSDSolver(surf, 0.8, th, alpha_deg=1.0, **g).steady()
    rm = LiftingTSDSolver(surf, 0.8, th, alpha_deg=-1.0, **g).steady()
    lin = TSDSolver(surf, 0.5, None, **g)
    a = np.radians(1.0)
    phi = lin.solve_normalwash(0.0, np.full((lin.nspan * lin.nchord, 1), -a))
    cl_lin = 4 * phi.reshape(*lin.shape, -1)[lin.iw1, :lin.nspan, 0, 0].real * lin.gmap / lin.b[:lin.nspan]
    rl = LiftingTSDSolver(surf, 0.5, None, alpha_deg=1.0, **g).steady()
    AR, beta = 8.0, np.sqrt(0.75)
    cla_h = 2 * np.pi * AR / (2 + np.sqrt(AR ** 2 * beta ** 2 + 4))
    out = {"alpha0_phi": float(np.abs(r1.phi - r0.phi).max()),
           "alpha0_gaf": float(np.abs(qb - qa).max() / np.abs(qa).max()),
           "sign_cp": float(np.abs(rp.cp - rm.cp_lower).max()),
           "sign_CL": float(abs(rp.CL + rm.CL) / rp.CL),
           "linear_cl": float(np.abs(rl.cl - cl_lin).max() / np.abs(cl_lin).max()),
           "CLa": rl.CL / a, "CLa_helmbold": cla_h}
    if verbose:
        print("1. consistency")
        print(f"   alpha = 0 vs symmetric solver: max |dphi| = {out['alpha0_phi']:.2e}, "
              f"GAF rel. diff = {out['alpha0_gaf']:.2e}")
        print(f"   alpha = +-1 deg: |Cp_u(+) - Cp_l(-)| = {out['sign_cp']:.2e}, CL antisymmetry {out['sign_CL']:.2e}")
        print(f"   M 0.5 flat, alpha 1 deg: cl vs linear TSD rel. diff {out['linear_cl']:.2e}; "
              f"CL_alpha = {out['CLa']:.3f} /rad (Helmbold {cla_h:.3f}, {100 * (out['CLa'] / cla_h - 1):+.1f} %)")
    return out


def naca0012_sweep(alphas=(0.0, 0.5, 1.0, 1.25, 1.5, 2.0), nchord=40, verbose=True):
    surf = Surface("w", (0, 0, 0), 1.0, (0, 8.0, 0), 1.0, nchord=10, nspan=10)
    sol = LiftingTSDSolver(surf, 0.8, naca4_thickness(0.12), nchord=nchord, nspan=8, verbose=False)
    rows = []
    if verbose:
        print(f"2. NACA 0012, M = 0.8, rectangular A = 16, {nchord} chord points (root strip)")
        print("   alpha    CL      cl_root  xs_upper  dCp_u   xs_lower  dCp_l   it conv [s]")
    res = []
    for i, a in enumerate(alphas):
        sol.set_alpha(a)
        t0 = time.perf_counter()
        r = sol.steady(warm_start=i > 0)
        res.append(r)
        rows.append((a, r.CL, r.cl[0], r.shock_upper[0], r.jump_upper[0], r.shock_lower[0], r.jump_lower[0]))
        if verbose:
            print(f"   {a:5.2f}  {r.CL:7.4f}  {r.cl[0]:7.4f}  {r.shock_upper[0]:7.3f}  {r.jump_upper[0]:6.3f}  "
                  f"{r.shock_lower[0]:7.3f}  {r.jump_lower[0]:6.3f}  {r.iterations:3d}  {'y' if r.converged else 'N'}  "
                  f"[{time.perf_counter() - t0:.0f}]")
    rows = np.array(rows, float)
    if verbose:
        i = int(np.argmin(np.abs(rows[:, 0] - AR211["alpha"])))
        print(f"   AGARD-AR-211 Euler, alpha {AR211['alpha']} deg: cl ~ {AR211['CL']}, "
              f"xs upper ~ {AR211['xs_upper']}, lower ~ {AR211['xs_lower']} (weak)")
        print(f"   TSD root strip, alpha {rows[i, 0]} deg: cl {rows[i, 2]:.3f}, xs upper {rows[i, 3]:.3f}, "
              f"lower {rows[i, 5]:.3f}")
        cl = rows[:, 2]
        if cl.min() <= AR211["CL"] <= cl.max():
            xs = np.interp(AR211["CL"], cl, rows[:, 3])
            print(f"   at equal section cl = {AR211['CL']}: TSD upper shock x/c = {xs:.3f} "
                  f"(alpha {np.interp(AR211['CL'], cl, rows[:, 0]):.2f} deg)")
    return rows, res


def agard_alpha(alphas=(0.0, 1.0, 2.0), mach=0.954, verbose=True):
    from flaap.flutter import identify, pk_method
    from flaap.models.agard445 import (B_REF, EXPERIMENT_FREQ_HZ, FLUTTER_DATA, OMEGA_ALPHA, agard445_fe,
                                       agard445_surface, frustum_volume)
    fe = agard445_fe(n_chord=16, n_span=16)
    modal = fe.modes(4, name="AGARD 445.6")
    modal.omega = 2 * np.pi * EXPERIMENT_FREQ_HZ.copy()
    u = np.zeros(fe.ndof)
    u[2::6] = 1.0
    m_ratio = float(u @ (fe.M @ u)) / (14.315 * frustum_volume())
    surf = agard445_surface(8, 10)
    funcs, names = modes_from_structure(modal)
    d = FLUTTER_DATA[mach]
    rho = d["rho"] * m_ratio
    K = np.array([1e-4, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.9])
    rows = []
    if verbose:
        print(f"3. AGARD 445.6, M = {mach}: flutter vs mean incidence (experiment at alpha = 0: FSI {d['FSI']})")
        print("   alpha    CL     xs_up(root)  FSI     f/f_alpha  conv [s]")
    for a in alphas:
        t0 = time.perf_counter()
        sol = LiftingTSDSolver(surf, mach, naca65a004(), alpha_deg=a)
        tab = tsd_gaf_table(surf, mach, K, B_REF, funcs, names, thickness=naca65a004(), solver=sol, keep_cp=False)
        mf = sol.mean_flow
        res = pk_method(modal.M, modal.K, tab, rho, np.linspace(50.0, 700.0, 261))
        identify(res)
        fl = [p for p in res.points if p.kind == "flutter"]
        c = min(fl, key=lambda p: p.velocity) if fl else None
        fsi = np.nan if c is None else c.velocity / (B_REF * OMEGA_ALPHA * np.sqrt(d["mu"]))
        fr = np.nan if c is None else 2 * np.pi * c.freq_hz / OMEGA_ALPHA
        rows.append((a, mf.CL, mf.shock_upper[0], fsi, fr))
        if verbose:
            print(f"   {a:5.2f}  {mf.CL:7.4f}  {mf.shock_upper[0]:9.3f}  {fsi:8.4f}  {fr:8.4f}   {'y' if mf.converged else 'N'}  "
                  f"[{time.perf_counter() - t0:.0f}]")
    return np.array(rows, float)


def make_figure(rows, res, ag, path):
    """Root Cp, shock position vs alpha and (if ``ag``) AGARD FSI vs mean alpha."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = 2 if ag is None else 3
    fig, ax = plt.subplots(1, n, figsize=(5 * n, 4))
    cmap = plt.get_cmap("viridis")
    for i, r in enumerate(res):
        c = cmap(i / max(1, len(res) - 1))
        ax[0].plot(r.xi[0], r.cp[0], color=c, label=f"{r.alpha_deg:g} deg")
        ax[0].plot(r.xi[0], r.cp_lower[0], color=c, ls="--")
    ax[0].axhline(res[0].cp_star, color="k", ls=":", lw=0.8)
    ax[0].invert_yaxis(); ax[0].set_xlabel("x/c"); ax[0].set_ylabel("root Cp (solid upper, dashed lower)")
    ax[0].set_title("NACA 0012, M 0.8, A = 16"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    ax[1].plot(rows[:, 0], rows[:, 3], "o-", label="TSD upper shock (root)")
    ax[1].plot(rows[:, 0], rows[:, 5], "s--", label="TSD lower shock (root)")
    ax[1].plot([AR211["alpha"]] * 2, [AR211["xs_upper"], AR211["xs_lower"]], "k*", ms=10,
               label="AGARD-AR-211 Euler (2-D)")
    ax[1].set_xlabel("alpha [deg]"); ax[1].set_ylabel("shock x/c"); ax[1].set_ylim(0, 1)
    ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    if ag is not None:
        ax[2].plot(ag[:, 0], ag[:, 3], "r^-", label="TSD lifting mean flow + PK")
        ax[2].plot([0.0], [0.3059], "ko", label="experiment, alpha = 0 (Yates 1987)")
        ax[2].set_xlabel("mean alpha [deg]"); ax[2].set_ylabel("flutter speed index")
        ax[2].set_title("AGARD 445.6, M 0.954"); ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=120)
    plt.close(fig)
    print("saved", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-agard", action="store_true")
    ap.add_argument("--alphas-agard", type=float, nargs="*", default=[0.0, 1.0, 2.0])
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    consistency()
    rows, res = naca0012_sweep()
    ag = None if a.skip_agard else agard_alpha(tuple(a.alphas_agard))
    if a.plot:
        make_figure(rows, res, ag, a.plot)


if __name__ == "__main__":
    main()
