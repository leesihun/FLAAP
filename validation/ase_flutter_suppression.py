"""Validation: aeroservoelastic flutter suppression and gust-load alleviation.

Case 1: active flutter suppression (AFS) on a typical section.  The plant is
the Vasconcellos (2016) / Conner (1997) pitch-plunge section, with the flap
as the control surface (exact Theodorsen RFA, Jones lags, Conner damping).
The actuator is second order (15 Hz, zeta = 0.7) and the sensors read pitch,
pitch rate, plunge and plunge rate.

* Zero gain must reproduce the open-loop p-method flutter speed of
  ``AeroelasticSystem.flutter`` (flap locked by the actuator).
* Two static output-feedback laws ``beta_cmd = K y`` are compared: law A
  uses pitch and pitch rate; law B is full output feedback, with gains from a
  small random search and Nelder-Mead optimisation of V_f, then rounded.
  Their closed-loop V-g / root-locus diagrams are plotted.
* With either law the open-loop flutter mode is stabilised.  The
  closed-loop limit is a mode whose frequency falls to zero, i.e. a
  closed-loop divergence (see the V-g plot).
* Gain/phase margins and the Nyquist criterion are evaluated at several
  speeds and checked against the closed-loop eigenvalues.  The loop is
  scaled by 0.97 and 1.03 times the gain margin to confirm the root crossing.
  Bode and Nyquist plots are produced.

Case 2: gust-load alleviation (GLA) on the Goland wing with a 50-83 % span
trailing-edge flap (:mod:`flaap.models.goland_flap`, 2 % modal damping).  A
1-cos gust (Uds = 5 m/s) acts at V = 100 m/s.  The controller is a root
bending strain-gauge loop: ``flap_cmd = k M_x,sg``, where ``M_x,sg`` is the
mode-displacement root bending moment, with a 20 Hz actuator.  Root bending
moment (summation of forces) is compared with and without control, and the
closed-loop flutter speed is checked.

Run:  python validation/ase_flutter_suppression.py
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
FIG = os.path.join(HERE, "figures")

from flaap.ase import ASEModel, Actuator, Sensor, bode_data, gain  # noqa: E402
from flaap.models.typical_section import section_system, vasconcellos2016  # noqa: E402

LAW_SIGNALS = ["pitch", "pitch_rate", "plunge", "plunge_rate"]
LAWS = {
    "A: pitch + pitch rate": np.array([[-1.0, -0.03, 0.0, 0.0]]),
    "B: full output feedback": np.array([[1.6, 0.16, -18.0, -2.3]]),
}


def _section():
    sec = vasconcellos2016()
    sysc, _ = section_system(sec, flap="control")
    act = Actuator.second_order("beta", 2 * np.pi * 15.0, 0.7)
    sens = [Sensor("pitch", [0.0, 1.0], "disp"), Sensor("pitch_rate", [0.0, 1.0], "rate"),
            Sensor("plunge", [1.0, 0.0], "disp"), Sensor("plunge_rate", [1.0, 0.0], "rate")]
    ase = ASEModel(sysc, [act], sens, sec.rho)
    out = {}
    V = np.linspace(1.0, 25.0, 241)
    t0 = time.perf_counter()
    ol = sysc.flutter(V, sec.rho)
    out["V_f_open_pmethod"] = ol.critical().velocity
    out["f_open"] = ol.critical().freq_hz
    zero = gain(np.zeros((1, 4)), LAW_SIGNALS, "beta_cmd")
    r0 = ase.flutter(V, zero)
    out["V_f_zero_gain"] = r0.critical().velocity
    out["V_f_zero_gain_bisect"] = ase.flutter_speed((1.0, 25.0), zero, n=97)
    out["vg"] = {"open": r0}
    out["laws"] = {}
    for name, Kg in LAWS.items():
        K = gain(Kg, LAW_SIGNALS, "beta_cmd")
        r = ase.flutter(V, K)
        vf = ase.flutter_speed((1.0, 25.0), K, n=97)
        out["vg"][name] = r
        out["laws"][name] = {"gains": Kg.ravel(), "V_f": vf, "increase": vf / out["V_f_open_pmethod"] - 1}
    out["t_flutter"] = time.perf_counter() - t0

    # margins of law B
    K = gain(LAWS["B: full output feedback"], LAW_SIGNALS, "beta_cmd")
    rows = []
    t0 = time.perf_counter()
    for Vm in (5.0, 9.0, 12.0, 15.0, 17.0, 18.0):
        m = ase.margins(Vm, K)
        chk = []
        for f in (m.gm_lower, m.gm_upper):
            if 0 < f < np.inf:
                lo = np.max(np.linalg.eigvals(ase.closed_loop(Vm, K, 0.97 * f, gust=False, extra_outputs=[]).A).real)
                hi = np.max(np.linalg.eigvals(ase.closed_loop(Vm, K, 1.03 * f, gust=False, extra_outputs=[]).A).real)
                chk.append(bool(lo * hi < 0))
        rows.append(dict(V=Vm, P=m.P, N=m.N_ccw, nyq=m.nyquist_stable, eig=m.closed_loop_stable,
                         gm_lo=m.gm_lower_db, gm_up=m.gm_upper_db, pm=m.pm, w_pm=m.w_pm,
                         gm_check=all(chk) if chk else None, margins=m))
    out["margins"] = rows
    out["t_margins"] = time.perf_counter() - t0
    w, L, _ = ase.loop(12.0, K)
    out["loop12"] = {"w": w, "L": L, **bode_data(L, w)}
    out["ase"] = ase
    return out


def _gla():
    from flaap.loads import MonitorStation, modal_load_matrices, station_output_group
    from flaap.models.goland import SPAN, X_EA
    from flaap.models.goland_flap import goland_flap
    from flaap.ase import beam_sensor
    t0 = time.perf_counter()
    g = goland_flap(zeta=0.02)
    sysc = g["ctrl"]
    st = [MonitorStation("root", [0.0, 0.01, 0.0], [0.0, 1.0, 0.0], ("Fz", "Mx", "My"))]
    station_output_group(sysc, g["table"], g["mesh"], g["modal"], st, "fsm")
    Tk = modal_load_matrices(g["modal"], st[0])[1][:, :sysc.ns]
    act = Actuator.second_order("flap", 2 * np.pi * 20.0, 0.7)
    sens = [Sensor("Mx_sg", Tk[1], "disp"),
            beam_sensor(g["modal"], X_EA, 0.95 * SPAN, X_EA, "acc", "acc_tip")]
    ase = ASEModel(sysc, [act], sens, g["rho"])
    out = {"t_build": time.perf_counter() - t0}
    k = -2.0e-6                                   # rad per N m
    K = gain(k, "Mx_sg", "flap_cmd")
    out["gain"] = k
    V, Uds = 100.0, 5.0
    out["V"], out["Uds"] = V, Uds
    rows = []
    hist = {}
    t0 = time.perf_counter()
    for H in (5.0, 15.0, 30.0):
        r0 = ase.gust_response(V, H, Uds, None, t_end=1.0, dt=5e-4, outputs=["root.Mx", "root.Fz"])
        r1 = ase.gust_response(V, H, Uds, K, t_end=1.0, dt=5e-4, outputs=["root.Mx", "root.Fz"])
        rows.append(dict(H=H, Mx_open=r0["peak"]["root.Mx"], Mx_cl=r1["peak"]["root.Mx"],
                         Fz_open=r0["peak"]["root.Fz"], Fz_cl=r1["peak"]["root.Fz"],
                         flap_deg=np.degrees(r1["peak"]["flap"]), acc_open=r0["peak"]["acc_tip"],
                         acc_cl=r1["peak"]["acc_tip"]))
        hist[H] = (r0, r1)
    out["t_gust"] = time.perf_counter() - t0
    out["rows"] = rows
    out["hist"] = hist
    Vs = np.linspace(40.0, 250.0, 106)
    out["V_f_open"] = sysc.flutter(Vs, g["rho"]).critical().velocity
    out["V_f_zero"] = ase.flutter(Vs, None).critical().velocity
    out["V_f_closed"] = ase.flutter_speed((40.0, 250.0), K, n=85)
    m = ase.margins(V, K)
    out["margins"] = m
    return out


def _figures(S, G):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(FIG, exist_ok=True)
    files = []
    cols = {"open": "k", "A: pitch + pitch rate": "tab:blue", "B: full output feedback": "tab:red"}

    # V-g and root locus
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for name, r in S["vg"].items():
        lab = "open loop (K = 0)" if name == "open" else name
        for j in range(r.nbranch):
            ax[0].plot(r.V[:, j], r.g[:, j], color=cols[name], lw=1.2, label=lab if j == 0 else None)
            ax[1].plot(r.V[:, j], r.freq[:, j], color=cols[name], lw=1.2)
            w = 2 * np.pi * r.freq[:, j]
            if name == "open":
                ax[2].plot(r.g[:, j] * w / 2, w, "k-", lw=0.8)
            elif name.startswith("B"):
                sc = ax[2].scatter(r.g[:, j] * w / 2, w, c=r.V[:, j], s=4, cmap="viridis")
    ax[0].axhline(0, color="gray", lw=0.6)
    ax[0].set_ylim(-0.4, 0.15)
    ax[0].set_ylabel("damping g")
    ax[1].set_ylabel("frequency [Hz]")
    for a in ax[:2]:
        a.set_xlabel("V [m/s]")
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    ax[2].axvline(0, color="gray", lw=0.6)
    ax[2].set_xlabel("Re(lambda) [1/s]")
    ax[2].set_ylabel("Im(lambda) [rad/s]")
    ax[2].set_xlim(-8, 4)
    ax[2].set_title("root locus vs V: law B (colour = V), open loop (black)", fontsize=10)
    fig.colorbar(sc, ax=ax[2], label="V [m/s]")
    ax[2].grid(alpha=0.3)
    fig.tight_layout()
    files.append(os.path.join(FIG, "ase_vg_rootlocus.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    # Bode + Nyquist at V = 12 m/s (open loop unstable, closed loop stable)
    lp = S["loop12"]
    m12 = [r for r in S["margins"] if r["V"] == 12.0][0]["margins"]
    fig = plt.figure(figsize=(13, 5))
    a1 = fig.add_subplot(2, 2, 1)
    a2 = fig.add_subplot(2, 2, 3, sharex=a1)
    a3 = fig.add_subplot(1, 2, 2)
    a1.semilogx(lp["w"], lp["mag_db"])
    a1.axhline(0, color="gray", lw=0.6)
    a1.set_ylabel("|L| [dB]")
    a2.semilogx(lp["w"], lp["phase_deg"])
    ph = lp["phase_deg"]
    for k in range(-3, 3):
        if ph.min() - 30 <= -180 + 360 * k <= ph.max() + 30:
            a2.axhline(-180 + 360 * k, color="gray", lw=0.6, ls=":")
    a2.set_ylabel("phase L [deg]")
    a2.set_xlabel("omega [rad/s]")
    a1.set_title("Loop transfer L = -K G, law B, V = 12 m/s")
    L = lp["L"]
    a3.plot(L.real, L.imag, "b-", lw=1, label="omega > 0")
    a3.plot(L.real, -L.imag, "b--", lw=0.6, label="omega < 0")
    a3.plot([-1], [0], "r+", ms=12, mew=2)
    th = np.linspace(0, 2 * np.pi, 200)
    a3.plot(np.cos(th), np.sin(th), color="gray", lw=0.5)
    a3.set_aspect("equal")
    lim = max(3.0, 1.2 / max(m12.gm_lower, 1e-3)) if m12.gm_lower > 0 else 3.0
    a3.set_xlim(-lim, lim)
    a3.set_ylim(-lim, lim)
    a3.set_title(f"Nyquist: P = {m12.P}, N_ccw = {m12.N_ccw}, PM = {m12.pm:.1f} deg, "
                 f"GM_low = {m12.gm_lower_db:.1f} dB", fontsize=9)
    a3.legend(fontsize=8)
    for a in (a1, a2, a3):
        a.grid(alpha=0.3)
    fig.tight_layout()
    files.append(os.path.join(FIG, "ase_bode_nyquist.png"))
    fig.savefig(files[-1], dpi=110)
    plt.close(fig)

    if G is not None:
        fig, ax = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
        r0, r1 = G["hist"][15.0]
        io = r0["outputs"].index("root.Mx")
        ax[0].plot(r0["t"], r0["wg"], "k")
        ax[0].set_ylabel("w_g [m/s]")
        ax[1].plot(r0["t"], r0["y"][:, io] / 1e3, "k", label="open loop")
        ax[1].plot(r1["t"], r1["y"][:, r1["outputs"].index("root.Mx")] / 1e3, "tab:red", label="GLA")
        ax[1].set_ylabel("root M_x [kN m]")
        ax[1].legend(fontsize=8)
        ax[2].plot(r1["t"], np.degrees(r1["y"][:, r1["outputs"].index("flap")]), "tab:red")
        ax[2].set_ylabel("flap [deg]")
        ax[2].set_xlabel("t [s]")
        ax[0].set_title(f"Goland + flap, 1-cos gust H = 15 m, Uds = {G['Uds']} m/s, V = {G['V']} m/s")
        for a in ax:
            a.grid(alpha=0.3)
        fig.tight_layout()
        files.append(os.path.join(FIG, "ase_gla_goland.png"))
        fig.savefig(files[-1], dpi=110)
        plt.close(fig)
    return files


def run(verbose: bool = True, gla: bool = True, figures: bool = True) -> dict:
    """Run the ASE validation and return a dict of results and arrays."""
    t0 = time.perf_counter()
    S = _section()
    G = _gla() if gla else None
    files = _figures(S, G) if figures else []
    res = {"section": {k: v for k, v in S.items() if k != "ase"}, "gla": G, "figures": files,
           "runtime": time.perf_counter() - t0}
    if verbose:
        _print(res, S, G)
    return res


def _print(res, S, G):
    print("=" * 78)
    print("ASE validation: flutter suppression on the Vasconcellos/Conner section")
    print("=" * 78)
    print(f"open-loop p-method (AeroelasticSystem.flutter): V_f = {S['V_f_open_pmethod']:.4f} m/s, "
          f"f = {S['f_open']:.3f} Hz")
    print(f"ASE closed loop with K = 0:                     V_f = {S['V_f_zero_gain']:.4f} m/s "
          f"(bisection {S['V_f_zero_gain_bisect']:.4f})")
    for name, d in S["laws"].items():
        print(f"law {name:28s} gains {np.array2string(d['gains'], precision=3)}: V_f = {d['V_f']:.3f} m/s "
              f"({100 * d['increase']:+.1f} %)")
    print("\nmargins, law B (GM in dB; lower = gain reduction that destabilises):")
    print("   V    P  N_ccw  Nyquist  eig     GM_low   GM_up     PM[deg] @w[rad/s]  GM check")
    for r in S["margins"]:
        print(f"{r['V']:5.1f}  {r['P']:2d}  {r['N']:4d}   {str(r['nyq']):7s} {str(r['eig']):7s} "
              f"{r['gm_lo']:7.2f}  {r['gm_up']:7.2f}   {r['pm']:7.2f} @ {r['w_pm']:6.2f}   {r['gm_check']}")
    print(f"timing: flutter {S['t_flutter']:.1f} s, margins {S['t_margins']:.1f} s")
    if G is not None:
        print("\n" + "-" * 78)
        print(f"GLA, Goland + flap: root-bending strain-gauge feedback k = {G['gain']:.1e} rad/(N m), "
              f"V = {G['V']} m/s, Uds = {G['Uds']} m/s")
        print("   H[m]   Mx open    Mx GLA   change   Fz open   Fz GLA   change   flap[deg]")
        for r in G["rows"]:
            print(f"{r['H']:6.1f}  {r['Mx_open']:9.0f} {r['Mx_cl']:9.0f}  {100 * (r['Mx_cl'] / r['Mx_open'] - 1):+6.1f}%"
                  f"  {r['Fz_open']:8.0f} {r['Fz_cl']:8.0f}  {100 * (r['Fz_cl'] / r['Fz_open'] - 1):+6.1f}%  "
                  f"{r['flap_deg']:6.2f}")
        print(f"flutter: open {G['V_f_open']:.2f} m/s, ASE K = 0 {G['V_f_zero']:.2f} m/s, "
              f"with GLA {G['V_f_closed']:.2f} m/s")
        m = G["margins"]
        print("GLA loop at 100 m/s: " + m.summary().replace("\n", "; "))
        print(f"timing: build {G['t_build']:.1f} s, gust runs {G['t_gust']:.1f} s")
    print(f"\nfigures: {', '.join(os.path.basename(f) for f in res['figures'])}")
    print(f"total runtime {res['runtime']:.1f} s")


if __name__ == "__main__":
    run(verbose=True)
