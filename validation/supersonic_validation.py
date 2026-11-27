"""Validation of the supersonic aerodynamics (piston theory and Mach box).

Cases
-----
1. Ackeret 2-D lift slope 4/beta (piston theory; Mach box on an A = 40 wing).
2. Rectangular wing, supersonic tips: CL_a = 4/beta (1 - 1/(2 beta A)).
3. Delta wing with subsonic leading edges: CL_a = 2 pi tan(eps) / E(sqrt(1 - m^2)),
   m = beta tan(eps) (Stewart 1946; Jones & Cohen).
4. 2-D oscillating airfoil (plunge/pitch) vs the Garrick-Rubinow (1946)
   integral, evaluated by independent quadrature.
5. Panel flutter with piston theory (coalescence, lambda = 2 q a^3 / (beta D)):
   2-D simply supported 343.3, 2-D clamped 636.6 (Hedgepeth 1957; Dowell 1975),
   square simply supported plate 512.6 (Dowell 1970).
6. AGARD 445.6 weakened model 3 at M = 1.072 and 1.141 vs Yates (1987) FSI
   0.3201 / 0.4031 (piston theory and Mach box, PK method, matched mass ratio).

Writes ``validation/supersonic_results.txt`` and
``validation/figures/supersonic_validation.png``.

Run:  python validation/supersonic_validation.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
from scipy.special import ellipe

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from flaap.aero.gaf import AeroModes, default_k_list, rigid_aero_modes  # noqa: E402
from flaap.aero.supersonic import (MachBox, garrick_rubinow_2d, piston_theory_table,  # noqa: E402
                                   supersonic_gaf_table)
from flaap.flutter import identify, pk_method  # noqa: E402
from flaap.geometry import Surface, build_mesh  # noqa: E402

OUT_TXT = os.path.join(HERE, "supersonic_results.txt")
OUT_FIG = os.path.join(HERE, "figures", "supersonic_validation.png")

LINES: list[str] = []


def out(s=""):
    print(s)
    LINES.append(s)


def row(name, val, ref, extra=""):
    err = 100.0 * (val / ref - 1.0)
    out(f"  {name:<52s} {val:10.4f} {ref:10.4f} {err:+8.2f} %  {extra}")
    return err


# ---------------------------------------------------------------- steady
def cl_alpha(surfaces, mach, nx, symmetry="symmetric"):
    mb = MachBox(surfaces, mach, nx=nx, symmetry=symmetry)
    t0 = time.perf_counter()
    _, d = mb.table(lambda p: (np.zeros((len(p), 1)), -np.ones((len(p), 1))), 1.0, [0.0],
                    gust=False)
    a = mb.area[mb.wing]
    return float(np.real(d[0, :, 0] @ a) / a.sum()), time.perf_counter() - t0, mb


def steady_cases(fig_data):
    M = 2.0
    b = np.sqrt(M * M - 1)
    out("1-3. Steady lift slopes at M = 2 (beta = %.4f)" % b)
    out(f"  {'case':<52s} {'FLAAP':>10s} {'theory':>10s} {'error':>10s}")
    mesh = build_mesh([Surface("w", (0, 0, 0), 1.0, (0, 1.0, 0), 1.0, nchord=8, nspan=2)])
    tab = piston_theory_table(mesh, rigid_aero_modes(mesh, 0.25, ("plunge", "alpha")), M, 0.5,
                              [0.0], gust=False)
    row("piston (Van Dyke), Ackeret 4/beta", tab.Q[0, 0, 1].real / mesh.total_area, 4 / b)
    tab = piston_theory_table(mesh, rigid_aero_modes(mesh, 0.25, ("plunge", "alpha")), M, 0.5,
                              [0.0], gust=False, van_dyke=False)
    out(f"  {'piston 1st order 4/M (hypersonic limit)':<52s} "
        f"{tab.Q[0, 0, 1].real / mesh.total_area:10.4f}   (4/M = {4 / M:.4f})")
    errs = {}
    for name, surf, ref in (
            ("Mach box rectangle A = 40", Surface("w", (0, 0, 0), 1.0, (0, 20.0, 0), 1.0),
             4 / b * (1 - 1 / (2 * b * 40))),
            ("Mach box rectangle A = 2", Surface("w", (0, 0, 0), 1.0, (0, 1.0, 0), 1.0),
             4 / b * (1 - 1 / (2 * b * 2))),
            ("Mach box delta tan(eps) = 0.4 (m = 0.69)", Surface("d", (0, 0, 0), 1.0, (1.0, 0.4, 0), 0.0),
             2 * np.pi * 0.4 / ellipe(1 - (b * 0.4) ** 2))):
        errs[name] = []
        for nx in (20, 40, 80):
            v, t, mb = cl_alpha([surf], M, nx)
            errs[name].append((nx, row(f"{name}, nx = {nx}", v, ref,
                                       f"[{mb.n_wing} wing boxes, {t:.2f} s]")))
    fig_data["steady"] = errs
    out()


# ---------------------------------------------------------------- 2-D unsteady
def machbox_2d(mach, k, b, c, xr, nx):
    mb = MachBox([Surface("w", (0, 0, 0), c, (0, 50 * c, 0), c)], mach, nx=nx, symmetry="symmetric")
    kq = k / b
    pc, pf = mb.wing_points("center"), mb.wing_points("segment")
    W = np.stack([1j * kq * np.ones(len(pc)), -1 - 1j * kq * (pc[:, 0] - xr)], -1)
    _, d = mb.solve(k, b, W)
    sel = mb.Y[mb.wing] < mb.dy                  # root strip (2-D region)
    h = np.stack([np.ones(len(pf)), -(pf[:, 0] - xr)], -1)[sel]
    return np.einsum("xi,xj,x->ij", h, d[sel], mb.seg_len[mb.wing][sel])


def piston_2d(mach, k, b, c, xr, n=200):
    mesh = build_mesh([Surface("w", (0, 0, 0), c, (0, 1.0, 0), c, nchord=n, nspan=1)])
    tab = piston_theory_table(mesh, rigid_aero_modes(mesh, xr, ("plunge", "pitch")), mach, b, [k],
                              gust=False)
    return tab.Q[0]


def unsteady_cases(fig_data):
    b, c, xr = 0.5, 1.0, 0.5
    out("4. 2-D oscillating flat plate, plunge h = 1 and pitch h = -(x - c/2), b = c/2")
    out("   Q[i,j] = int h_i dCp_j dx; error = max|Q - Q_GR| / max|Q_GR| over the 2x2 matrix")
    out(f"  {'M':>5s} {'k':>5s}  {'Q_GR pitch->lift':>22s}  {'MB nx=40':>9s} {'MB nx=80':>9s} "
        f"{'piston':>9s}")
    curves = {}
    for M in (1.2, 2.0):
        for k in (0.1, 0.5, 1.0):
            ref = garrick_rubinow_2d(M, k, b, c, xr)
            e = [np.abs(machbox_2d(M, k, b, c, xr, nx) - ref).max() / np.abs(ref).max()
                 for nx in (40, 80)]
            ep = np.abs(piston_2d(M, k, b, c, xr) - ref).max() / np.abs(ref).max()
            out(f"  {M:5.2f} {k:5.2f}  {ref[0, 1].real:+10.4f}{ref[0, 1].imag:+10.4f}j  "
                f"{100 * e[0]:8.3f}% {100 * e[1]:8.3f}% {100 * ep:8.2f}%")
        ks = np.linspace(0.02, 1.5, 25)
        curves[M] = (ks, np.array([garrick_rubinow_2d(M, k, b, c, xr)[0, 1] for k in ks]),
                     np.array([machbox_2d(M, k, b, c, xr, 40)[0, 1] for k in ks]),
                     np.array([piston_2d(M, k, b, c, xr)[0, 1] for k in ks]))
    fig_data["unsteady"] = curves
    out("   (Van Dyke piston theory is the small-k / large-M expansion: its error grows with k")
    out("    and as M -> 1, as expected.)")
    out()


# ---------------------------------------------------------------- panel flutter
_BL = [4.730040745, 7.853204624, 10.99560784, 14.13716549, 17.27875966,
       20.42035225, 23.56194490, 26.70353756]


def panel_modes(x, y, kind, n, width):
    h, hx, K = [], [], []
    sy = np.sin(np.pi * y / width) if kind == "square" else np.ones_like(x)
    for i in range(n):
        if kind in ("ss", "square"):
            a = (i + 1) * np.pi
            h.append(np.sin(a * x) * sy)
            hx.append(a * np.cos(a * x) * sy)
            K.append(a ** 4 / 2 if kind == "ss" else (a * a + (np.pi / width) ** 2) ** 2 * width / 4)
        else:
            b = _BL[i]
            s = (np.cosh(b) - np.cos(b)) / (np.sinh(b) - np.sin(b))
            h.append(np.cosh(b * x) - np.cos(b * x) - s * (np.sinh(b * x) - np.sin(b * x)))
            hx.append(b * (np.sinh(b * x) + np.sin(b * x) - s * (np.cosh(b * x) - np.cos(b * x))))
            K.append(b ** 4)
    m = {"ss": 0.5, "clamped": 1.0, "square": width / 4}[kind]
    return np.array(h).T, np.array(hx).T, np.diag(K), np.eye(n) * m


def panel_flutter(kind, n=8, nc=100, mach=2.0, loci=False):
    """lambda_cr = 2 q a^3 / (beta D) by eigenvalue coalescence (a = 1, D = 1),
    piston theory on one face (half of the two-sided lifting pressure)."""
    ns = 10 if kind == "square" else 1
    mesh = build_mesh([Surface("p", (0, 0, 0), 1.0, (0, 1.0, 0), 1.0, nchord=nc, nspan=ns)])
    hc, hxc, K, Mm = panel_modes(mesh.colloc[:, 0], mesh.colloc[:, 1], kind, n, 1.0)
    hf = panel_modes(mesh.force_point[:, 0], mesh.force_point[:, 1], kind, n, 1.0)[0]
    tab = piston_theory_table(mesh, AeroModes(hc, hxc, hf), mach, 1.0, [0.0], gust=False)
    A = 0.5 * tab.Q[0].real
    beta = np.sqrt(mach * mach - 1)

    def eig(lam):
        return np.linalg.eigvals(np.linalg.solve(Mm, K - 0.5 * lam * beta * A))

    def unstable(lam):
        ev = eig(lam)
        return np.abs(ev.imag).max() > 1e-8 * np.abs(ev).max()

    lo, hi = 10.0, 3000.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    lam_cr = 0.5 * (lo + hi)
    if not loci:
        return lam_cr
    ls = np.linspace(0, 1.3 * lam_cr, 200)
    return lam_cr, ls, np.array([np.sort_complex(eig(x))[:3] for x in ls])


def panel_cases(fig_data):
    out("5. Panel flutter, piston theory, lambda_cr = 2 q a^3 / (beta D) (eigenvalue coalescence)")
    out(f"  {'case':<52s} {'FLAAP':>10s} {'reference':>10s} {'error':>10s}")
    for kind, ref, label, nc in (("ss", 343.3, "2-D simply supported (8 modes)", 100),
                                 ("clamped", 636.6, "2-D clamped (8 modes)", 100),
                                 ("square", 512.6, "square plate, simply supported (8 modes)", 40)):
        t0 = time.perf_counter()
        lam = panel_flutter(kind, nc=nc)
        row(f"{label}, {nc} boxes", lam, ref, f"[{time.perf_counter() - t0:.2f} s]")
    lam2 = panel_flutter("ss", n=2)
    out(f"  {'2-D simply supported, 2-mode Galerkin':<52s} {lam2:10.4f}   "
        f"(analytic 2-mode: 3 (k2 - k1)/8 = {45 * np.pi ** 4 / 16:.4f})")
    out("   References: Hedgepeth (1957) / Dowell, Aeroelasticity of Plates and Shells (1975);")
    out("   the often-quoted 512 is the square simply supported plate, 343.3 the 2-D panel.")
    fig_data["panel"] = panel_flutter("ss", loci=True)
    out()


# ---------------------------------------------------------------- AGARD 445.6
def agard_cases(fig_data):
    from flaap.models.agard445 import (B_REF, EXPERIMENT_FREQ_HZ, FLUTTER_DATA, OMEGA_ALPHA,
                                       agard445_fe, agard445_surface, frustum_volume)
    from flaap.structures.spline import surface_spline
    out("6. AGARD 445.6 weakened model 3, supersonic points (PK, matched mass ratio,")
    out("   MITC4 16x16 FE mode shapes; DLM mesh 12x16 for the mode interpolation)")
    fe = agard445_fe(n_chord=16, n_span=16)
    modal = fe.modes(4)
    fe_omega = modal.omega.copy()
    u = np.zeros(fe.ndof)
    u[2::6] = 1.0
    m_model = float(u @ (fe.M @ u))
    m_test = 14.315 * frustum_volume()
    mesh = build_mesh([agard445_surface(12, 16)], "symmetric")
    sp = surface_spline(mesh, modal.nodes, method="ips")
    modes = AeroModes.from_spline(sp, modal.dof_vectors(), [f"m{i + 1}" for i in range(4)])
    kv = default_k_list(1.2, 14)
    kv[0] = 1e-4
    V = np.linspace(20.0, 1500.0, 371)
    methods = (("piston, Van Dyke", "piston", {}),
               ("piston, 1st order 4/M", "piston", {"van_dyke": False}),
               ("Mach box nx = 40", "machbox", {"nx": 40}),
               ("Mach box nx = 80", "machbox", {"nx": 80}),
               ("Mach box nx = 160", "machbox", {"nx": 160}))
    res_tab = {}
    for freq in ("measured", "FE"):
        modal.omega = 2 * np.pi * EXPERIMENT_FREQ_HZ.copy() if freq == "measured" else fe_omega.copy()
        out(f"  natural frequencies: {freq} ({np.round(modal.omega / 2 / np.pi, 2)} Hz)")
        out(f"  {'M':>6s} {'method':<24s} {'FSI':>8s} {'FSI exp':>8s} {'err':>8s} "
            f"{'wf/wa':>7s} {'exp':>7s}  time")
        for M in (1.072, 1.141):
            d = FLUTTER_DATA[M]
            rho_eq = d["rho"] * m_model / m_test
            for label, meth, opt in methods:
                t0 = time.perf_counter()
                tab = supersonic_gaf_table(mesh, modes, M, B_REF, kv, method=meth, gust=False,
                                           keep_cp=False, **opt)
                res = pk_method(modal.M, modal.K, tab, rho_eq, V)
                identify(res)
                fl = [p for p in res.points if p.kind == "flutter"]
                dt = time.perf_counter() - t0
                if not fl:
                    out(f"  {M:6.3f} {label:<24s} {'none':>8s} {d['FSI']:8.4f} "
                        f"{'':>8s} {'':>7s} {d['freq_ratio']:7.4f}  {dt:.1f} s"
                        f"  (no flutter up to {V[-1]:.0f} m/s)")
                    res_tab[(freq, M, label)] = np.nan
                    continue
                c = min(fl, key=lambda p: p.velocity)
                fsi = c.velocity / (B_REF * OMEGA_ALPHA * np.sqrt(d["mu"]))
                fr = 2 * np.pi * c.freq_hz / OMEGA_ALPHA
                res_tab[(freq, M, label)] = fsi
                out(f"  {M:6.3f} {label:<24s} {fsi:8.4f} {d['FSI']:8.4f} "
                    f"{100 * (fsi / d['FSI'] - 1):+7.1f}% {fr:7.4f} {d['freq_ratio']:7.4f}  {dt:.1f} s")
        out()
    fig_data["agard"] = (res_tab, {M: FLUTTER_DATA[M]["FSI"] for M in (1.072, 1.141)})
    out("  Notes: piston theory is a high-Mach / high-frequency approximation and is not")
    out("  applicable at M = 1.07-1.14 (Van Dyke: 4/beta stiffness with a negative damping factor")
    out("  (M^2-2)/(M^2-1): no flutter found; 1st order 4/M: none at 1.072, FSI 3.6x at 1.141).")
    out("  The Mach box (linear potential theory) is 8-14 % above the test with the measured")
    out("  frequencies; the experiment is at very low supersonic Mach where thickness and")
    out("  viscous (shock-boundary layer) effects not in linear theory are significant.")
    out()


# ---------------------------------------------------------------- figure
def make_figure(fig_data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    a = ax[0, 0]
    for name, e in fig_data["steady"].items():
        nx, err = np.array(e).T
        a.loglog(nx, np.abs(err), "o-", label=name)
    a.set_xlabel("boxes along the chord nx")
    a.set_ylabel("|CL_alpha error| [%]")
    a.set_title("Steady lift slope, M = 2")
    a.legend(fontsize=8)
    a = ax[0, 1]
    for M, (ks, gr, mb, ps) in fig_data["unsteady"].items():
        line, = a.plot(ks, gr.real, "-", label=f"Garrick-Rubinow M={M} (Re)")
        a.plot(ks, gr.imag, "--", color=line.get_color(), label="(Im)")
        a.plot(ks, mb.real, "o", ks, mb.imag, "s", color=line.get_color(), ms=4, mfc="none")
        a.plot(ks, ps.real, ":", ks, ps.imag, ":", color=line.get_color(), lw=1)
    a.set_xlabel("k = w b / U")
    a.set_ylabel("Q lift due to pitch")
    a.set_title("2-D oscillating plate (o/s Mach box nx=40, dotted piston)")
    a.legend(fontsize=7)
    a = ax[1, 0]
    lam_cr, ls, ev = fig_data["panel"]
    for j in range(ev.shape[1]):
        a.plot(ls, ev[:, j].real, "k-", lw=1)
        a.plot(ls, ev[:, j].imag, "r--", lw=1)
    a.axvline(lam_cr, color="b", ls=":")
    a.text(lam_cr, a.get_ylim()[1] * 0.9, f"  lambda_cr = {lam_cr:.1f}\n  (exact 343.3)", color="b")
    a.set_xlabel("lambda = 2 q a^3 / (beta D)")
    a.set_ylabel("Re / Im of omega^2 (nondim.)")
    a.set_title("2-D simply supported panel, piston theory")
    a = ax[1, 1]
    res_tab, exp = fig_data["agard"]
    labels = [lab for lab in dict.fromkeys(k[2] for k in res_tab)
              if "4/M" not in lab and any(np.isfinite(v) for kk, v in res_tab.items() if kk[2] == lab)]
    Ms = sorted(exp)
    w = 0.8 / (len(labels) + 1)
    xs = np.arange(len(Ms))
    a.bar(xs - 0.4 + w / 2, [exp[M] for M in Ms], w, color="k", label="experiment (Yates 1987)")
    for i, lab in enumerate(labels):
        a.bar(xs - 0.4 + w * (i + 1.5), [res_tab.get(("measured", M, lab), np.nan) for M in Ms], w,
              label=lab)
    a.set_xticks(xs)
    a.set_xticklabels([f"M = {M}" for M in Ms])
    a.set_ylabel("flutter speed index")
    a.set_title("AGARD 445.6 (measured frequencies)")
    a.legend(fontsize=8)
    for x in ax.ravel():
        x.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_FIG), exist_ok=True)
    fig.savefig(OUT_FIG, dpi=110)
    plt.close(fig)


def main():
    t0 = time.perf_counter()
    out("FLAAP supersonic aerodynamics validation (flaap.aero.supersonic)")
    out("=" * 78)
    out()
    fig_data: dict = {}
    steady_cases(fig_data)
    unsteady_cases(fig_data)
    panel_cases(fig_data)
    agard_cases(fig_data)
    make_figure(fig_data)
    out(f"total run time {time.perf_counter() - t0:.1f} s")
    with open(OUT_TXT, "w", encoding="utf-8") as f:
        f.write("\n".join(LINES) + "\n")
    print("wrote", OUT_TXT, "and", OUT_FIG)


if __name__ == "__main__":
    main()
