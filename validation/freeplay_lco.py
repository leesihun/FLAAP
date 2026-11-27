"""Validation: control-surface freeplay limit-cycle oscillations.

Case 1 is the typical section of Vasconcellos et al. (2016), a replication of
the Conner, Tang and Dowell (1997) pitch-plunge-flap experiment.  It uses the
exact Theodorsen RFA with Jones lags and the measured Conner damping, with flap
freeplay of +-delta0.  The checks are:

* the linear flutter speed with the nominal hinge stiffness (published 9.2 m/s)
  and with zero hinge stiffness (V_f0, the lower edge of the LCO band);
* the describing-function (DF) LCO amplitude and frequency against exact
  piecewise-linear time marching, at several speeds;
* bifurcation diagrams from upward and downward speed sweeps with continuation
  (subcritical fold and hysteresis);
* homogeneity: the LCO amplitude scales with delta0 (delta0 = 1 deg vs 2 deg);
* preloaded freeplay: a gap centred at 0.5 delta0 gives a biased LCO, while a
  gap centred at 1.5 delta0 (flap resting on the stiff branch) gives none;
* time histories and phase planes.

Case 2 is the Goland wing with a 50-83 % span trailing-edge flap
(:mod:`flaap.models.goland_flap`), 4 wing modes, a DLM GAF table, and 1 %
modal damping.  It compares the DF and time marching of the 3-D model.

Run:  python validation/freeplay_lco.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
FIG = os.path.join(HERE, "figures")

from flaap.models.typical_section import (VASCONCELLOS2016_REFERENCE, section_system,  # noqa: E402
                                          vasconcellos2016)
from flaap.nonlinear import (Bilinear, Cubic, Freeplay, NonlinearAeroelastic,  # noqa: E402
                             bifurcation_sweep, describing_function_lco, flutter_speed_vs_stiffness,
                             preloaded_freeplay)


def _section(verbose):
    sec = vasconcellos2016()
    sysm, _ = section_system(sec)
    d0 = np.radians(1.0)
    m = NonlinearAeroelastic(sysm, "beta", Freeplay(sec.K_beta, d0))
    out = {"delta0": d0}

    t0 = time.perf_counter()
    df = describing_function_lco(m, sec.rho, np.linspace(0.5, 11.0, 85))
    out["t_df"] = time.perf_counter() - t0
    out["df"] = df
    out["V_f_linear"] = df.V_fK
    out["V_f0"] = df.V_f0
    st = [i for b in df.branches() for i in b if df.stable[i]]
    out["V_fold"] = float(min(df.V[i] for i in st))

    # time marching vs DF
    rows = []
    t0 = time.perf_counter()
    for V in (3.5, 5.0, 7.0, 8.5):
        r = m.simulate(V, sec.rho, 80.0, x0={"delta": 3 * d0}, dt=0.005)
        info = r.lco()
        a_df, w_df = df.amplitude(V), df.frequency(V)
        rows.append(dict(V=V, kind=info.kind, A_tm=info.amplitude / d0, A1_tm=info.amplitude_h1 / d0,
                         w_tm=info.omega, A_df=float(a_df[0]) / d0, w_df=float(w_df[0]),
                         err_A=info.amplitude / a_df[0] - 1, err_w=info.omega / w_df[0] - 1,
                         alpha_deg=np.degrees(r.lco("alpha").amplitude), h_mm=1e3 * r.lco("h").amplitude))
    out["t_tm"] = time.perf_counter() - t0
    out["tm_vs_df"] = rows

    # bifurcation sweeps (continuation up and down)
    Vs = np.arange(1.0, 9.51, 0.5)
    t0 = time.perf_counter()
    out["sweep_up"] = bifurcation_sweep(m, Vs, sec.rho, 60.0, x0={"delta": 3 * d0}, dt=0.005)
    out["sweep_down"] = bifurcation_sweep(m, Vs[::-1], sec.rho, 60.0, x0={"delta": 3 * d0}, dt=0.005)
    out["t_sweep"] = time.perf_counter() - t0

    # amplitude scaling with delta0
    m2 = NonlinearAeroelastic(sysm, "beta", Freeplay(sec.K_beta, 2 * d0))
    sc = []
    for V in (4.0, 7.0):
        r1 = m.simulate(V, sec.rho, 60.0, x0={"delta": 3 * d0}, dt=0.005)
        r2 = m2.simulate(V, sec.rho, 60.0, x0={"delta": 6 * d0}, dt=0.005)
        sc.append((V, r2.lco().amplitude / r1.lco().amplitude, r2.lco("alpha").amplitude / r1.lco("alpha").amplitude))
    out["scaling"] = sc

    # below the fold: decay for any initial amplitude
    Vb = 0.8 * out["V_fold"]
    out["below_fold"] = [(Vb, a0, m.simulate(Vb, sec.rho, 60.0, x0={"delta": a0 * d0, "alpha": 0.02},
                                             dt=0.005).lco().kind) for a0 in (1.5, 4.0, 10.0)]

    # histories for the figures
    out["hist"] = {V: m.simulate(V, sec.rho, 30.0, x0={"delta": 3 * d0}, dt=0.002) for V in (2.5, 5.0, 8.5)}

    # preloaded freeplay and other nonlinearities (time marching, V = 6 m/s)
    out["preload"] = []
    for pl in (0.5, 1.5):
        pre = NonlinearAeroelastic(sysm, "beta", preloaded_freeplay(sec.K_beta, d0, pl * d0))
        rp = pre.simulate(6.0, sec.rho, 60.0, x0={"delta": 3 * d0}, dt=0.005)
        out["preload"].append((pl, rp.lco()))
    out["section"] = sec
    out["model"] = m
    return out


def _goland(verbose):
    from flaap.models.goland_flap import goland_flap
    t0 = time.perf_counter()
    g = goland_flap(zeta=0.01)
    out = {"t_build": time.perf_counter() - t0, "K_h": g["K_h"], "I_h": g["I_h"],
           "omega_flap": g["omega_flap"], "M_qd": g["M_qd"]}
    rho = g["rho"]
    d0 = np.radians(0.5)
    m = NonlinearAeroelastic(g["hinge"], "flap", Freeplay(g["K_h"], d0))
    out["delta0"] = d0
    Kr = np.array([0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0])
    out["vf_K"] = flutter_speed_vs_stiffness(m, rho, g["K_h"] * Kr, (20.0, 250.0))
    out["vf_K"]["K_ratio"] = Kr
    t0 = time.perf_counter()
    df = describing_function_lco(m, rho, np.linspace(130.0, 160.0, 61))
    out["t_df"] = time.perf_counter() - t0
    out["df"] = df
    rows = []
    t0 = time.perf_counter()
    for V in (140.0, 145.0, 150.0, 153.0):
        r = m.simulate(V, rho, 8.0, x0={"flap": 3 * d0, "mode2": 1e-3}, dt=2e-4)
        info = r.lco()
        a, w = df.amplitude(V), df.frequency(V)
        rows.append(dict(V=V, kind=info.kind, A_tm=info.amplitude / d0, w_tm=info.omega,
                         A_df=float(a[0]) / d0 if a.size else np.nan, w_df=float(w[0]) if w.size else np.nan))
        if V == 150.0:
            out["hist"] = r
    out["t_tm"] = time.perf_counter() - t0
    out["tm_vs_df"] = rows
    return out


def _figures(S, G):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    d0 = S["delta0"]
    df = S["df"]
    files = []

    # bifurcation diagram
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for b in df.branches():
        b = np.asarray(b)
        stb = df.stable[b]
        for mask, ls, lab in ((stb, "-", "DF stable"), (~stb, "--", "DF unstable")):
            if mask.any():
                ax[0].plot(df.V[b][mask], df.A[b][mask] / d0, "k" + ls, lw=1.4, label=lab)
                ax[1].plot(df.V[b][mask], df.omega[b][mask], "k" + ls, lw=1.4)
    for sw, mk, lab in ((S["sweep_up"], "o", "time marching, V up"), (S["sweep_down"], "x", "time marching, V down")):
        k = np.array(sw["kind"]) == "lco"
        ax[0].plot(sw["V"][k], sw["amp"][k] / d0, mk, ms=6, mfc="none", label=lab)
        ax[0].plot(sw["V"][~k], 0 * sw["V"][~k], mk, ms=6, mfc="none", color="gray")
        ax[1].plot(sw["V"][k], 2 * np.pi * sw["freq"][k], mk, ms=6, mfc="none")
    for a in ax:
        a.axvline(S["V_f0"], color="tab:red", ls=":", lw=1)
        a.axvline(S["V_f_linear"], color="tab:blue", ls=":", lw=1)
        a.set_xlabel("V [m/s]")
        a.grid(alpha=0.3)
    ax[0].text(S["V_f0"], 0.2, " V_f(K=0)", color="tab:red", fontsize=8)
    ax[0].text(S["V_f_linear"], 0.2, " V_f(K_beta)", color="tab:blue", fontsize=8, ha="right")
    ax[0].set_ylabel("flap LCO amplitude / delta0")
    ax[1].set_ylabel("LCO frequency [rad/s]")
    h, l = ax[0].get_legend_handles_labels()
    u = dict(zip(l, h))
    ax[0].legend(u.values(), u.keys(), fontsize=8, loc="upper left")
    ax[0].set_ylim(0, 4)
    fig.suptitle("Pitch-plunge-flap section with flap freeplay (Vasconcellos 2016 / Conner 1997)")
    fig.tight_layout()
    files.append(os.path.join(FIG, "freeplay_bifurcation.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    # time histories
    fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    for V, r in S["hist"].items():
        ax[0].plot(r.t, np.degrees(r.dof("beta")), lw=0.8, label=f"V = {V} m/s")
        ax[1].plot(r.t, np.degrees(r.dof("alpha")), lw=0.8)
        ax[2].plot(r.t, 1e3 * r.dof("h"), lw=0.8)
    ax[0].axhline(np.degrees(d0), color="k", lw=0.5, ls=":")
    ax[0].axhline(-np.degrees(d0), color="k", lw=0.5, ls=":")
    ax[0].set_ylabel("flap [deg]")
    ax[1].set_ylabel("pitch [deg]")
    ax[2].set_ylabel("plunge [mm]")
    ax[2].set_xlabel("t [s]")
    ax[0].legend(fontsize=8, ncol=3)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    files.append(os.path.join(FIG, "freeplay_time_history.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    # phase planes
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for V, r in S["hist"].items():
        k = r.t > 15.0
        ax[0].plot(np.degrees(r.dof("beta")[k]), np.degrees(r.rate("beta")[k]), lw=0.7, label=f"V = {V} m/s")
        ax[1].plot(np.degrees(r.dof("alpha")[k]), np.degrees(r.rate("alpha")[k]), lw=0.7)
    for x in (-1, 1):
        ax[0].axvline(x * np.degrees(d0), color="k", lw=0.5, ls=":")
    ax[0].set_xlabel("flap [deg]")
    ax[0].set_ylabel("flap rate [deg/s]")
    ax[1].set_xlabel("pitch [deg]")
    ax[1].set_ylabel("pitch rate [deg/s]")
    ax[0].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle("Phase planes, t > 15 s")
    fig.tight_layout()
    files.append(os.path.join(FIG, "freeplay_phase_plane.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    # hinge characteristics and equivalent stiffness
    K = 1.0
    nls = {"freeplay": Freeplay(K, 1.0), "bilinear (K_in = 0.3 K)": Bilinear(K, 1.0, 0.3),
           "preloaded freeplay (1.5 delta0)": preloaded_freeplay(K, 1.0, 1.5),
           "cubic (K3 = 0.1 K)": Cubic(K, 0.1, 1.0)}
    x = np.linspace(-5, 5, 401)
    A = np.linspace(0.01, 8, 300)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for name, nl in nls.items():
        ax[0].plot(x, nl.moment(x), label=name)
        ax[1].plot(A, nl.keq(A) / K, label=name)
    ax[0].set_xlabel("delta / delta0")
    ax[0].set_ylabel("M_h / (K delta0)")
    ax[1].set_xlabel("A / delta0")
    ax[1].set_ylabel("K_eq / K (B = 0)")
    ax[0].legend(fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    files.append(os.path.join(FIG, "freeplay_keq.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    if G is not None:
        dg, gd = G["df"], G["delta0"]
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        vk = G["vf_K"]
        ax[0].plot(vk["K_ratio"], vk["V_f"], "o-")
        ax[0].set_xlabel("K / K_h")
        ax[0].set_ylabel("V_f [m/s]")
        ax[0].set_title("Flutter speed vs hinge stiffness")
        for b in dg.branches():
            b = np.asarray(b)
            ax[1].plot(dg.V[b], dg.A[b] / gd, "k-", lw=1.4, label="DF")
        for row in G["tm_vs_df"]:
            ax[1].plot(row["V"], row["A_tm"] if row["kind"] == "lco" else 0.0, "o", mfc="none",
                       color="tab:orange", label="time marching")
        ax[1].axvline(dg.V_f0, color="tab:red", ls=":")
        ax[1].axvline(dg.V_fK, color="tab:blue", ls=":")
        ax[1].set_xlabel("V [m/s]")
        ax[1].set_ylabel("flap LCO amplitude / delta0")
        h, l = ax[1].get_legend_handles_labels()
        u = dict(zip(l, h))
        ax[1].legend(u.values(), u.keys(), fontsize=8)
        ax[1].set_title("Goland + flap, delta0 = 0.5 deg")
        r = G["hist"]
        ax[2].plot(r.t, np.degrees(r.dof("flap")), lw=0.6)
        ax[2].set_xlabel("t [s]")
        ax[2].set_ylabel("flap [deg]")
        ax[2].set_title("V = 150 m/s")
        for a in ax:
            a.grid(alpha=0.3)
        fig.tight_layout()
        files.append(os.path.join(FIG, "goland_freeplay_lco.png"))
        fig.savefig(files[-1], dpi=110)
        plt.close(fig)
    return files


def run(verbose: bool = True, goland: bool = True, figures: bool = True) -> dict:
    """Run the freeplay validation and return a dict of results and arrays."""
    t0 = time.perf_counter()
    S = _section(verbose)
    G = _goland(verbose) if goland else None
    files = _figures(S, G) if figures else []
    res = {
        "section": {k: v for k, v in S.items() if k not in ("model", "section", "hist")},
        "goland": G,
        "figures": files,
        "runtime": time.perf_counter() - t0,
        "reference": VASCONCELLOS2016_REFERENCE,
    }
    res["max_err_amplitude"] = max(abs(r["err_A"]) for r in S["tm_vs_df"])
    res["max_err_frequency"] = max(abs(r["err_w"]) for r in S["tm_vs_df"])
    if verbose:
        _print(res, S, G)
    return res


def _print(res, S, G):
    d0 = S["delta0"]
    print("=" * 78)
    print("Freeplay LCO validation: Vasconcellos (2016) / Conner (1997) section, delta0 = 1 deg")
    print("=" * 78)
    print(f"linear flutter speed (K = K_beta): {S['V_f_linear']:.3f} m/s   "
          f"(published {VASCONCELLOS2016_REFERENCE['U_f_linear']} m/s)")
    print(f"flutter speed with K = 0         : {S['V_f0']:.3f} m/s")
    print(f"DF stable LCO branch from        : {S['V_fold']:.2f} m/s (subcritical fold)")
    print("\n  V     kind   A_TM/d0  A1_TM/d0  A_DF/d0   err_A   w_TM    w_DF    err_w  pitch[deg] h[mm]")
    for r in S["tm_vs_df"]:
        print(f"{r['V']:5.1f}  {r['kind']:5s}  {r['A_tm']:7.3f}  {r['A1_tm']:7.3f}  {r['A_df']:7.3f}  "
              f"{100 * r['err_A']:+5.1f}%  {r['w_tm']:6.2f}  {r['w_df']:6.2f}  {100 * r['err_w']:+5.1f}%  "
              f"{r['alpha_deg']:7.3f}  {r['h_mm']:6.2f}")
    print("\nbifurcation sweep (time marching, continuation):")
    up, dn = S["sweep_up"], S["sweep_down"]
    print("  V     up: kind  A/d0     down: kind  A/d0")
    amp_dn = dict(zip(dn["V"], zip(dn["kind"], dn["amp"])))
    for V, k, a in zip(up["V"], up["kind"], up["amp"]):
        kd, ad = amp_dn[V]
        fa = f"{a / d0:6.3f}" if k != "growing" else "   -  "
        fd = f"{ad / d0:6.3f}" if kd != "growing" else "   -  "
        print(f"{V:5.1f}   {k:9s} {fa}   {kd:9s} {fd}")
    print("\namplitude scaling (delta0 = 2 deg vs 1 deg):")
    for V, rf, ra in S["scaling"]:
        print(f"  V = {V:4.1f} m/s: flap ratio {rf:.4f}, pitch ratio {ra:.4f}")
    print("below the fold:", ", ".join(f"V={V:.2f} A0={a}d0 -> {k}" for V, a, k in S["below_fold"]))
    for pl, p in S["preload"]:
        print(f"preloaded freeplay (gap centre {pl} delta0) at 6 m/s: {p.kind}, amplitude "
              f"{p.amplitude / d0:.3f} d0, mean {np.degrees(p.mean):.3f} deg")
    print(f"timing: DF {S['t_df']:.2f} s (85 speeds), TM {S['t_tm']:.2f} s (4 x 80 s), sweeps {S['t_sweep']:.1f} s")
    if G is not None:
        print("\n" + "-" * 78)
        print("Goland wing + 50-83 % span flap (hinge 75 % c), zeta = 1 %, delta0 = 0.5 deg")
        print(f"K_h = {G['K_h']:.1f} N m/rad, I_h = {G['I_h']:.4f} kg m^2, omega_flap = {G['omega_flap']:.2f} rad/s")
        vk = G["vf_K"]
        print("V_f(K):", ", ".join(f"{k:.2f}K_h:{v:.1f}" for k, v in zip(vk["K_ratio"], vk["V_f"])))
        print(f"DF: V_f0 = {G['df'].V_f0:.2f} m/s, V_f(K_h) = {G['df'].V_fK:.2f} m/s")
        print("  V      kind   A_TM/d0  A_DF/d0   w_TM    w_DF")
        for r in G["tm_vs_df"]:
            print(f"{r['V']:6.1f}  {r['kind']:5s}  {r['A_tm']:7.3f}  {r['A_df']:7.3f}  {r['w_tm']:6.2f}  {r['w_df']:6.2f}")
        print(f"timing: build {G['t_build']:.1f} s, DF {G['t_df']:.2f} s, TM {G['t_tm']:.1f} s (4 x 8 s)")
    print(f"\nfigures: {', '.join(os.path.basename(f) for f in res['figures'])}")
    print(f"total runtime {res['runtime']:.1f} s")


if __name__ == "__main__":
    run(verbose=True)
