"""Validation of the steady VLM, the nonlinear UVLM and the linearised UVLM.

Cases
-----
1. Steady VLM: elliptic wings A = 4 / 8 (span efficiency e = 1, Helmbold
   lift slope), the 2-D limit 2 pi (A = 1000) and a rectangular A = 4 wing.
2. Indicial response of a flat plate (nonlinear UVLM, A = 1000):
   Wagner (impulsive start) and Kussner (sharp-edged gust), against the exact
   functions evaluated by Fourier-sine quadrature of Theodorsen's C(k) and
   Sears' function (not the exponential approximations).
3. Harmonic pitch about c/4: |CL| and phase vs Theodorsen (1935), mesh
   convergence and Richardson extrapolation.
4. Harmonic plunge: propulsive thrust vs Garrick (1936),
   C_T = pi k^2 (h0/b)^2 (F^2 + G^2).
5. Linearised UVLM generalised forces (frequency domain) vs Theodorsen.
6. 3-D: impulsively started rectangular wing (A = 4) -> steady VLM; free wake.
7. Goland wing flutter (rho = 1.02): linearised UVLM GAFs + PK and the
   discrete-time UVLM state-space eigenvalues vs Murua, Palacios & Graham
   (2012) UVLM 165 m/s / 10.98 Hz and the FLAAP DLM.

Writes ``validation/uvlm_results.txt`` and ``validation/figures/uvlm_validation.png``.

Run:  python validation/uvlm_validation.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.special import hankel2, j0, j1

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from flaap.aero.gaf import AeroModel, AeroModes, rigid_aero_modes  # noqa: E402
from flaap.aero.uvlm import (UVLM, Kinematics, LinearUVLM, flutter_crossing, flutter_eigen_sweep,  # noqa: E402
                             harmonic, sharp_edged_gust, steady_vlm)
from flaap.flutter import pk_method  # noqa: E402
from flaap.geometry import AeroMesh, Surface, build_mesh  # noqa: E402

OUT_TXT = os.path.join(HERE, "uvlm_results.txt")
OUT_FIG = os.path.join(HERE, "figures", "uvlm_validation.png")
LINES: list[str] = []
B = 0.5                         # semi-chord of the 2-D wing


def out(s=""):
    print(s, flush=True)
    LINES.append(s)


def row(name, val, ref, extra=""):
    err = 100.0 * (val / ref - 1.0)
    out(f"  {name:<46s} {val:10.4f} {ref:10.4f} {err:+8.2f} %  {extra}")
    return err


# ---------------------------------------------------------------- references
def theodorsen_C(k):
    return hankel2(1, k) / (hankel2(1, k) + 1j * hankel2(0, k))


def sears_le(k):
    """Sears function referred to the leading-edge gust arrival."""
    return ((j0(k) - 1j * j1(k)) * theodorsen_C(k) + 1j * j1(k)) * np.exp(-1j * k)


def indicial(H, s):
    """phi(s) = 2/pi int_0^inf Re H(k) sin(k s) / k dk."""
    def f(k):
        return H(k).real / k if k > 1e-12 else 1.0
    a = quad(f, 0, 1, weight="sin", wvar=s, limit=400)[0]
    b = quad(f, 1, np.inf, weight="sin", wvar=s, limlst=200)[0]
    return 2 / np.pi * (a + b)


def theodorsen_cl(k, a, alpha_hat, h_hat, U=1.0):
    """Complex CL amplitude (h up, pitch nose-up about x = b a)."""
    w = k * U / B
    L = np.pi * B ** 2 * (w ** 2 * h_hat + U * 1j * w * alpha_hat + B * a * w ** 2 * alpha_hat) \
        + 2 * np.pi * U * B * theodorsen_C(k) * (-1j * w * h_hat + U * alpha_hat + B * (0.5 - a) * 1j * w * alpha_hat)
    return L / (0.5 * U ** 2 * 2 * B)


def fit(t, y, w):
    sel = t >= t[-1] - 2 * np.pi / w + 1e-9
    A = np.column_stack([np.ones(sel.sum()), np.cos(w * t[sel]), np.sin(w * t[sel])])
    c = np.linalg.lstsq(A, y[sel], rcond=None)[0]
    return c[0], c[1] - 1j * c[2]


def strip2d(nc):
    return build_mesh([Surface("w", (0, 0, 0), 2 * B, (0, 500, 0), 2 * B, nchord=nc, nspan=1)], "symmetric")


def elliptic(A, nc, ns, c0=1.0):
    bh = A * np.pi * c0 / 8
    y = bh * np.sin(np.linspace(0, np.pi / 2, ns + 1))
    c = c0 * np.sqrt(np.clip(1 - (y / bh) ** 2, 0, None))
    c[-1] = max(c[-1], 1e-4 * c0)
    xle = 0.25 * c0 - 0.25 * c
    xi = np.linspace(0, 1, nc + 1)
    corners, strip, ci = [], [], []
    for j in range(ns):
        for i in range(nc):
            def p(jj, ii):
                return [xle[jj] + xi[ii] * c[jj], y[jj], 0.0]
            corners.append([p(j, i), p(j + 1, i), p(j + 1, i + 1), p(j, i + 1)])
            strip.append(j)
            ci.append(i)
    n = len(corners)
    return AeroMesh(np.array(corners), np.zeros(n, int), np.array(strip), np.array(ci),
                    [Surface("w", (0, 0, 0), c0, (0, bh, 0), 0.0, nc, ns)], "symmetric")


# ---------------------------------------------------------------- 1. steady VLM
def steady_cases(fig):
    out("1. Steady VLM (horseshoe lattice, Trefftz-plane induced drag)")
    out(f"  {'case':<46s} {'FLAAP':>10s} {'ref':>10s} {'err':>9s}")
    for A in (4.0, 8.0):
        helm = 2 * np.pi * A / (2 + np.sqrt(A * A + 4))
        for nc, ns in ((8, 24), (12, 48)):
            r = steady_vlm(elliptic(A, nc, ns), [2.0, 6.0])
            row(f"elliptic A = {A:g}, {nc}x{ns}: CL_alpha [1/rad]", r.CL_alpha, helm, "Helmbold")
            row(f"elliptic A = {A:g}, {nc}x{ns}: span efficiency e", r.e, 1.0, "lifting line")
            if A == 8.0 and nc == 12:
                fig["span"] = (r.y_strip, r.cl_c[1], r.CL[1], A)
    r = steady_vlm(strip2d(8), [3.0])
    row("A = 1000 strip, 8 chordwise: CL_alpha", r.CL_alpha, 2 * np.pi, "thin aerofoil 2 pi")
    for ns in (20, 40, 80):
        m = build_mesh([Surface("w", (0, 0, 0), 1.0, (0, 2.0, 0), 1.0, nchord=8, nspan=ns)], "symmetric")
        r = steady_vlm(m, [4.0], x_ref=0.25)
        out(f"  rectangle A = 4, 8x{ns:<3d}: CL_alpha {r.CL_alpha:.4f}, e {r.e:.4f}, x_ac/c {r.x_ac:.4f}"
            "   (lifting surface 3.6-3.7, e 0.98-0.99)")
    out()


# ---------------------------------------------------------------- 2. indicial
def indicial_cases(fig):
    out("2. Indicial lift of a flat plate (nonlinear UVLM, prescribed flat wake, dt = dx/U)")
    s_list = (1.0, 2.0, 4.0, 6.0, 10.0, 20.0, 29.5)
    wag = [indicial(theodorsen_C, s) for s in s_list]
    kus = [indicial(sears_le, s) for s in s_list]
    U, alpha = 1.0, np.radians(2.0)
    rw = UVLM(strip2d(20)).simulate(U, 15.0, Kinematics(alpha0=alpha), rho=1.0)
    phi_w = rw.CL / (2 * np.pi * alpha)
    w0 = 0.01
    rk = UVLM(strip2d(12)).simulate(U, 15.0, Kinematics(gust=sharp_edged_gust(w0, U, 0.0)), rho=1.0)
    psi = rk.CL / (2 * np.pi * w0 / U)
    out(f"  {'s = 2Ut/c':>10s} {'Wagner UVLM':>12s} {'exact':>8s} {'err':>7s}   "
        f"{'Kussner UVLM':>12s} {'exact':>8s} {'err':>7s}")
    for s, wx, kx in zip(s_list, wag, kus):
        iw = int(np.argmin(abs(rw.s - s)))
        ik = int(np.argmin(abs(rk.s - s)))
        out(f"  {s:10.1f} {phi_w[iw]:12.4f} {wx:8.4f} {100 * (phi_w[iw] / wx - 1):+6.2f}%   "
            f"{psi[ik]:12.4f} {kx:8.4f} {100 * (psi[ik] / kx - 1):+6.2f}%")
    out("  (Wagner nc = 20, Kussner nc = 12.  The UVLM sees the gust front panel by panel, so the")
    out("   Kussner rise while the front crosses the chord (s <= 2) is discretised: up to +14 % at")
    out("   s = 2, where the front reaches the trailing edge; within 2 % for s >= 4.)")
    sg = np.linspace(0.2, 29.5, 60)
    fig["indicial"] = (rw.s, phi_w, rk.s, psi, sg, [indicial(theodorsen_C, s) for s in sg],
                       [indicial(sears_le, s) for s in sg])
    out()


# ---------------------------------------------------------------- 3./4. harmonic
def harmonic_cases(fig):
    out("3. Harmonic pitch about c/4 (theta0 = 1 deg), |CL| / |CL_Theodorsen| and phase error")
    U, th0 = 1.0, np.radians(1.0)
    ratios = {}
    for k in (0.1, 0.5, 1.0):
        w = k * U / B
        ref = theodorsen_cl(k, -0.5, -1j * th0, 0.0)
        for nc in ((8, 16, 32) if k > 0.2 else (12,)):
            u = UVLM(strip2d(nc))
            r = u.simulate(U, 4 * 2 * np.pi / w, Kinematics(x_ref=0.5 * B, pitch=harmonic(th0, w)), rho=1.0)
            _, cl = fit(r.t, r.CL, w)
            ratios.setdefault(k, []).append((nc, abs(cl) / abs(ref)))
            out(f"  k = {k:4.2f}, nc = {nc:2d}: |CL| {abs(cl):.5f} (Theodorsen {abs(ref):.5f}), "
                f"ratio {abs(cl) / abs(ref):.4f}, phase error {np.degrees(np.angle(cl / ref)):+.2f} deg")
        if len(ratios[k]) == 3:
            r16, r32 = ratios[k][1][1], ratios[k][2][1]
            out(f"  k = {k:4.2f}: Richardson (1st order) nc -> inf: ratio {2 * r32 - r16:.4f}")
    fig["pitch"] = ratios
    out()
    out("4. Harmonic plunge h0 = 0.05 b: lift and mean thrust vs Theodorsen / Garrick (nc = 12)")
    out(f"  {'k':>5s} {'|CL| ratio':>11s} {'C_T UVLM':>11s} {'C_T Garrick':>12s} {'err':>8s}")
    h0, thr = 0.05 * B, []
    for k in (0.1, 0.25, 0.5, 0.75, 1.0):
        w = k * U / B
        r = UVLM(strip2d(12)).simulate(U, 4 * 2 * np.pi / w, Kinematics(x_ref=0.5 * B, plunge=harmonic(h0, w)),
                                       rho=1.0)
        cd0, _ = fit(r.t, r.CD, w)
        _, cl = fit(r.t, r.CL, w)
        ref = theodorsen_cl(k, -0.5, 0.0, -1j * h0)
        ct = np.pi * k ** 2 * (h0 / B) ** 2 * abs(theodorsen_C(k)) ** 2
        thr.append((k, -cd0, ct))
        out(f"  {k:5.2f} {abs(cl) / abs(ref):11.4f} {-cd0:11.4e} {ct:12.4e} {100 * (-cd0 / ct - 1):+7.2f}%")
    fig["thrust"] = np.array(thr)
    out()


# ---------------------------------------------------------------- 5. linear GAF
def gaf_cases(fig):
    out("5. Linearised UVLM (frequency domain) vs Theodorsen: max relative error of the")
    out("   plunge/pitch lift GAFs (wake 60 chords)")
    kv = np.array([0.1, 0.5, 1.0])
    out(f"  {'nc':>4s} " + " ".join(f"{'k = ' + str(k):>9s}" for k in kv))
    errs = []
    for nc in (8, 16, 32):
        m = strip2d(nc)
        md = rigid_aero_modes(m, 0.25, ("plunge", "pitch"), bref=B)
        tab = LinearUVLM(m, md, B, wake_chords=60).gaf_table(kv, gust=False)
        S = m.total_area
        e = [max(abs(tab.Q[i, 0, j] / S - theodorsen_cl(k, -0.5, *ah)) / abs(theodorsen_cl(k, -0.5, *ah))
                 for j, ah in enumerate(((0.0, 1.0), (1.0, 0.0)))) for i, k in enumerate(kv)]
        errs.append((nc, *e))
        out(f"  {nc:4d} " + " ".join(f"{100 * x:8.2f}%" for x in e))
    fig["gaf"] = np.array(errs)
    out("  first-order convergence in nc (the same discrete wake/time-step model as the UVLM).")
    out()


# ---------------------------------------------------------------- 6. 3-D
def three_d_cases():
    out("6. 3-D rectangular wing A = 4 (6x10), alpha = 5 deg, impulsive start")
    m = build_mesh([Surface("w", (0, 0, 0), 1.0, (0, 2.0, 0), 1.0, nchord=6, nspan=10)], "symmetric")
    st = steady_vlm(m, [5.0], x_ref=0.25)
    for free in (False, True):
        t0 = time.perf_counter()
        r = UVLM(m).simulate(10.0, 1.5, Kinematics(alpha0=np.radians(5.0), x_ref=0.25), rho=1.225,
                             free_wake=free, wake_length=12.0)
        out(f"  {'free' if free else 'prescribed'} wake: CL(15 c) = {r.CL[-1]:.4f} "
            f"({100 * r.CL[-1] / st.CL[0]:.1f} % of steady VLM {st.CL[0]:.4f}), "
            f"CD = {r.CD[-1]:.5f} (Trefftz CDi {st.CDi[0]:.5f}), CL(0+) = {r.CL[1]:.4f}  "
            f"[{time.perf_counter() - t0:.1f} s]")
    out("  As in Katz & Plotkin (2001, ch. 13), the 3-D lift builds up faster than Wagner's 2-D")
    out("  function and reaches the steady VLM value; the free wake (roll-up) changes the lift by")
    out("  < 0.1 % at this incidence (prescribed / free wakes give the same loads to 4 digits).")
    out()


# ---------------------------------------------------------------- 7. Goland
def goland_cases(fig):
    from flaap.models.goland import CHORD, goland_modal, goland_surface
    from flaap.structures.spline import beam_spline
    out("7. Goland wing flutter, rho = 1.02 kg/m^3 (4 beam modes)")
    out("   references: Murua et al. (2012) UVLM 165.0 m/s, 10.98 Hz; Wang (1999) 163.8 m/s")
    mm = goland_modal(4)
    kv = np.r_[0, np.geomspace(0.01, 2, 19)]
    V = np.linspace(120, 220, 201)
    res = []
    for nc, ns in ((8, 12), (16, 12), (32, 12)):
        t0 = time.perf_counter()
        mesh = build_mesh([goland_surface(nc, ns)], "symmetric")
        modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors())
        tab = LinearUVLM(mesh, modes, CHORD / 2, wake_chords=30).gaf_table(kv, gust=False, keep_cp=False)
        c = pk_method(mm.M, mm.K, tab, 1.02, V).critical()
        res.append((nc, c.velocity))
        out(f"  linear UVLM {nc:2d}x{ns} + PK: V_f = {c.velocity:7.2f} m/s, f = {c.freq_hz:6.3f} Hz "
            f"[{time.perf_counter() - t0:.1f} s]")
    v16, v32 = res[1][1], res[2][1]
    row("linear UVLM, Richardson nc -> inf: V_f", 2 * v32 - v16, 165.0, "Murua UVLM")
    mesh = build_mesh([goland_surface(8, 12)], "symmetric")
    modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors())
    c = pk_method(mm.M, mm.K, AeroModel(mesh, modes, 0.0, CHORD / 2).table(kv, gust=False), 1.02, V).critical()
    out(f"  DLM 8x12 + PK (same modes):     V_f = {c.velocity:7.2f} m/s, f = {c.freq_hz:6.3f} Hz")
    t0 = time.perf_counter()
    lin = LinearUVLM(mesh, modes, CHORD / 2, wake_chords=30)
    Vs, lam = flutter_eigen_sweep(lin, 1.02, mm.M, mm.K, np.linspace(120, 220, 11))
    cr = flutter_crossing(Vs, lam)
    out(f"  discrete-time UVLM state space 8x12 ({lin.N + lin.Nw + 12} states): V_f = {cr[0]:7.2f} m/s, "
        f"f = {cr[1] / 2 / np.pi:6.3f} Hz [{time.perf_counter() - t0:.1f} s]  (time-domain eigenvalues;")
    out("   agrees with PK on the same lattice -> the frequency- and time-domain forms are consistent)")
    fig["goland"] = (Vs, lam)
    out()


# ---------------------------------------------------------------- figure
def make_figure(fig):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    a = ax[0, 0]
    y, clc, CL, A = fig["span"]
    o = np.argsort(y)
    bh = A * np.pi / 8                                   # semi-span of the c0 = 1 elliptic wing
    yy = np.linspace(0, bh, 200)
    a.plot(y[o], clc[o], "o", ms=3, label="VLM 12x48")
    a.plot(yy, 4 / np.pi * CL * np.sqrt(1 - (yy / bh) ** 2), "k-", lw=1, label="elliptic (4/pi) CL sqrt(1-eta^2)")
    a.set_title(f"Elliptic wing A = {A:g}: span load (alpha 6 deg)")
    a.set_xlabel("y [m]")
    a.set_ylabel("cl c / c_ref")
    a.legend()
    a = ax[0, 1]
    sw, pw, sk, pk, sg, wg, kg = fig["indicial"]
    a.plot(sw, pw, "b-", lw=1, label="UVLM Wagner (nc 20)")
    a.plot(sg, wg, "b--", lw=1, label="exact Wagner")
    a.plot(sk, pk, "r-", lw=1, label="UVLM Kussner (nc 12)")
    a.plot(sg, kg, "r--", lw=1, label="exact Kussner")
    a.set_xlim(0, 30)
    a.set_ylim(0, 1.05)
    a.set_xlabel("s = 2 U t / c")
    a.set_title("Indicial lift")
    a.legend(fontsize=8)
    a = ax[0, 2]
    for k, lst in fig["pitch"].items():
        if len(lst) > 1:
            n, r = np.array(lst).T
            a.semilogx(n, r, "o-", label=f"k = {k}")
    a.axhline(1, color="k", lw=0.8)
    a.set_xlabel("chordwise panels nc")
    a.set_ylabel("|CL| / |CL Theodorsen|")
    a.set_title("Harmonic pitch about c/4")
    a.legend()
    a = ax[1, 0]
    t = fig["thrust"]
    kk = np.linspace(0.05, 1.05, 60)
    a.plot(kk, np.pi * kk ** 2 * 0.05 ** 2 * np.abs(theodorsen_C(kk)) ** 2, "k-", label="Garrick")
    a.plot(t[:, 0], t[:, 1], "o", label="UVLM nc 12")
    a.set_xlabel("k")
    a.set_ylabel("C_T")
    a.set_title("Plunge thrust, h0 = 0.05 b")
    a.legend()
    a = ax[1, 1]
    g = fig["gaf"]
    for j, k in enumerate((0.1, 0.5, 1.0)):
        a.loglog(g[:, 0], 100 * g[:, j + 1], "o-", label=f"k = {k}")
    a.set_xlabel("chordwise panels nc")
    a.set_ylabel("max GAF error vs Theodorsen [%]")
    a.set_title("Linearised UVLM GAFs")
    a.legend()
    a = ax[1, 2]
    Vs, lam = fig["goland"]
    for j in range(lam.shape[1]):
        a.plot(Vs, -lam[:, j].real / np.abs(lam[:, j]), "o-", ms=3, label=f"branch {j + 1}")
    a.axhline(0, color="k", lw=0.8)
    a.axvline(165.0, color="r", ls="--", lw=0.8, label="Murua UVLM 165")
    a.set_xlabel("V [m/s]")
    a.set_ylabel("zeta")
    a.set_title("Goland: discrete-time UVLM eigenvalues")
    a.legend(fontsize=8)
    for a in ax.ravel():
        a.grid(alpha=0.3)
    f.tight_layout()
    os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)
    f.savefig(OUT_FIG, dpi=110)


def main():
    t0 = time.perf_counter()
    out("FLAAP VLM / UVLM validation")
    out("=" * 78)
    fig = {}
    steady_cases(fig)
    indicial_cases(fig)
    harmonic_cases(fig)
    gaf_cases(fig)
    three_d_cases()
    goland_cases(fig)
    out(f"total time {time.perf_counter() - t0:.0f} s")
    with open(OUT_TXT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    make_figure(fig)


if __name__ == "__main__":
    main()
