"""Launcher routes for the nonlinear and aeroservoelastic extensions.

* ``model freeplay``: control-surface freeplay LCO (:mod:`flaap.nonlinear`):
  describing-function branches, time marching and bifurcation sweeps.
* ``model ase``: actuator + sensors + output feedback on the aeroelastic
  state space (:mod:`flaap.ase`): flutter suppression, loop margins and
  gust-load alleviation.

Both routes run on two built-in plants, chosen with ``structure``:

* ``section``: Vasconcellos (2016) / Conner (1997) pitch-plunge-flap
  typical section (exact Theodorsen RFA);
* ``goland_flap``: Goland wing with a 50-83 % span trailing-edge flap
  (DLM GAFs, Roger RFA with gust columns).

Defaults that depend on the plant (speed ranges, gains, gap size) are
chosen here and written to the log and to ``summary["defaults"]``.
"""
from __future__ import annotations

import time

import numpy as np

from .builder import velocity_range
from .runner import Output, _num, _plt, jlist, write_curves, write_flutter, write_timehist

SECTION_SENSORS = ("pitch", "pitch_rate", "plunge", "plunge_rate")
GOLAND_SENSORS = ("root_bending", "tip_accel")

PLANT_DEFAULTS = {
    "section": {
        "freeplay": dict(freeplay_deg=1.0, velocities=[0.5, 11.0, 85], t_end=60.0, dt=0.005),
        "ase": dict(velocities=[1.0, 25.0, 241], velocity=12.0, actuator_frequency=15.0,
                    sensors=list(SECTION_SENSORS), feedback_gains=[1.6, 0.16, -18.0, -2.3]),
    },
    "goland_flap": {
        "freeplay": dict(freeplay_deg=0.5, velocities=[130.0, 160.0, 61], t_end=8.0, dt=2e-4),
        "ase": dict(velocities=[40.0, 250.0, 106], velocity=100.0, actuator_frequency=20.0,
                    sensors=["root_bending"], feedback_gains=[-2.0e-6]),
    },
}


def _defaults(out: Output, v: dict, route: str) -> dict:
    """Fill plant-dependent defaults (only for keys the user did not give)."""
    v = dict(v)
    used = {}
    for k, d in PLANT_DEFAULTS[v["structure"]][route].items():
        if k not in v:
            v[k] = d
            used[k] = d
    if used:
        out.log("plant defaults: " + ", ".join(f"{k} = {d}" for k, d in used.items()))
        out.summary["defaults"] = used
    return v


def _as_list(x):
    return list(x) if isinstance(x, (list, tuple)) else [x]


# ======================================================================
# plants
# ======================================================================
def _plant(v: dict, flap: str):
    """(system, rho, info) for structure section | goland_flap.

    ``flap='structural'``: hinge as an elastic DOF (freeplay);
    ``flap='control'``: flap as a control input (ASE)."""
    if v["structure"] == "section":
        from ..models.typical_section import section_system, vasconcellos2016
        sec = vasconcellos2016()
        sysm, _ = section_system(sec, flap=flap)
        info = dict(name="Vasconcellos/Conner typical section", hinge="beta", K_h=sec.K_beta, sec=sec)
        return sysm, float(v.get("density", sec.rho)), info
    from ..models.goland_flap import goland_flap
    zeta = float(v.get("damping") or 0.0)
    g = goland_flap(nmodes=int(v.get("nmodes", 4)), zeta=zeta)
    sysm = g["hinge"] if flap == "structural" else g["ctrl"]
    info = dict(name="Goland wing + 50-83 % span flap", hinge="flap", K_h=g["K_h"], g=g)
    return sysm, float(v.get("density", g["rho"])), info


def _nonlinearity(v: dict, K: float):
    from ..nonlinear import Bilinear, Freeplay, preloaded_freeplay
    d0 = np.radians(float(v["freeplay_deg"]))
    kind = v.get("nonlinearity", "freeplay")
    if kind == "bilinear":
        return Bilinear(K, d0, float(v.get("inner_stiffness_ratio", 0.3)) * K), d0
    if kind == "preloaded":
        return preloaded_freeplay(K, d0, np.radians(float(v.get("preload_deg", 0.0)))), d0
    return Freeplay(K, d0), d0


# ======================================================================
# freeplay
# ======================================================================
def _x0(v, info, d0, sysm):
    amp = np.radians(float(v["initial_deflection_deg"])) if "initial_deflection_deg" in v else 3.0 * d0
    x0 = {"delta": amp}
    if info["hinge"] == "flap" and "mode2" in sysm.names:
        x0["mode2"] = 1e-3                    # seed the torsion mode as in the validation case
    return x0


def _write_history(out: Output, r, name: str, d0: float, title: str):
    series = {"hinge [deg]": np.degrees(r.delta)}
    for nm in r.names:
        if nm not in ("beta", "flap"):
            series[f"{nm}"] = r.dof(nm)
    write_timehist(out, name, r.t, series, title, meta={"V": r.V, "delta0_deg": np.degrees(d0)},
                   units={"hinge [deg]": "deg"})


def route_freeplay(out: Output, v: dict, base):
    from ..nonlinear import NonlinearAeroelastic, bifurcation_sweep, describing_function_lco
    v = _defaults(out, v, "freeplay")
    sysm, rho, info = _plant(v, "structural")
    nl, d0 = _nonlinearity(v, info["K_h"])
    m = NonlinearAeroelastic(sysm, info["hinge"], nl)
    out.log(f"plant: {info['name']}; hinge '{info['hinge']}', K_h = {info['K_h']:.5g} N m/rad, rho = {rho}")
    out.log(f"nonlinearity: {nl!r}")
    out.summary.update({"plant": info["name"], "rho": rho, "nonlinearity": repr(nl),
                        "delta0_deg": float(v["freeplay_deg"]), "K_h": info["K_h"]})
    t_end, dt = float(v["t_end"]), float(v["dt"])
    x0 = _x0(v, info, d0, sysm)
    mode = v["mode"]
    plots = v.get("plots", True)

    if mode == "describing_function":
        if v.get("nonlinearity") == "preloaded":
            raise ValueError("the describing-function analysis assumes an odd nonlinearity; "
                             "use mode time or bifurcation for a preloaded freeplay")
        t0 = time.perf_counter()
        Vs = velocity_range(v)
        df = describing_function_lco(m, rho, Vs)
        out.log(f"describing function over {Vs.size} speeds: {time.perf_counter() - t0:.1f} s")
        st = np.nonzero(df.stable)[0]
        V_fold = float(min(df.V[i] for i in st)) if st.size else float("nan")
        out.summary.update({"V_f_linear": _num(df.V_fK), "V_f0": _num(df.V_f0), "V_fold": _num(V_fold)})
        out.log(f"V_f(K_h) = {df.V_fK:.4g} m/s, V_f(K=0) = {df.V_f0:.4g} m/s, stable LCO from {V_fold:.4g} m/s")
        # amplitude / frequency versus V on the sweep grid (null where no branch)
        nmax = max([df.amplitude(V, False).size for V in Vs] + [1])
        cur = {}
        for lab, only in (("stable", True), ("all", False)):
            amps = [df.amplitude(V, only) for V in Vs]
            frq = [df.frequency(V, only) for V in Vs]
            for j in range(nmax if lab == "all" else max([a.size for a in amps] + [1])):
                cur[f"A/delta0:{lab} {j + 1}"] = [a[j] / d0 if a.size > j else np.nan for a in amps]
                cur[f"f [Hz]:{lab} {j + 1}"] = [w[j] / (2 * np.pi) if w.size > j else np.nan for w in frq]
        out.data("curves_lco_df", {"x": jlist(Vs, 6), "curves": {k: jlist(c, 6) for k, c in cur.items()},
                                        "xlabel": "V [m/s]", "ylabel": "A/delta0", "logy": False,
                                        "meta": {"V_f_linear": _num(df.V_fK), "V_f0": _num(df.V_f0),
                                                 "V_fold": _num(V_fold)}},
                 "curves", "Describing-function LCO amplitude and frequency")
        # time-marching check at a few speeds inside the stable band
        tm_V = _as_list(v["tm_velocities"]) if "tm_velocities" in v else (
            list(np.linspace(V_fold, min(df.V_fK, Vs[-1]), 5)[1:-1]) if np.isfinite(V_fold) else [])
        rows = []
        for V in tm_V:
            t0 = time.perf_counter()
            r = m.simulate(float(V), rho, t_end, x0=x0, dt=dt)
            li = r.lco()
            a_df, w_df = df.amplitude(float(V)), df.frequency(float(V))
            rows.append({"V": _num(V), "A_tm_deg": _num(np.degrees(li.amplitude)),
                         "A1_tm_deg": _num(np.degrees(li.amplitude_h1)), "kind": li.kind,
                         "f_tm_hz": _num(li.frequency),
                         "A_df_deg": _num(np.degrees(a_df.max())) if a_df.size else None,
                         "f_df_hz": _num(w_df[np.argmax(a_df)] / (2 * np.pi)) if a_df.size else None,
                         "t_s": round(time.perf_counter() - t0, 2)})
            out.log(f"  V = {V:7.3f}: time marching {li.kind}, A = {rows[-1]['A_tm_deg']} deg, "
                    f"f = {rows[-1]['f_tm_hz']} Hz | DF A = {rows[-1]['A_df_deg']} deg ({rows[-1]['t_s']} s)")
            if V == tm_V[len(tm_V) // 2]:
                _write_history(out, r, "lco", d0, f"Freeplay LCO at V = {V:.3g} m/s")
        out.summary["time_marching"] = rows
        if rows:
            write_curves(out, "lco_tm_vs_df", [r_["V"] for r_ in rows],
                         {"A [deg]:time marching": [r_["A_tm_deg"] for r_ in rows],
                          "A [deg]:describing function": [r_["A_df_deg"] for r_ in rows]},
                         "LCO amplitude: time marching vs describing function", "V [m/s]", "A [deg]")
        if plots:
            _plot_df(out, df, d0, rows)

    elif mode == "time":
        V = float(v["velocity"])
        r = m.simulate(V, rho, t_end, x0=x0, dt=dt)
        li = r.lco()
        out.summary["lco"] = {"V": V, "kind": li.kind, "amplitude_deg": _num(np.degrees(li.amplitude)),
                              "amplitude_h1_deg": _num(np.degrees(li.amplitude_h1)),
                              "mean_deg": _num(np.degrees(li.mean)), "freq_hz": _num(li.frequency),
                              "growth": _num(li.growth), "switches": int(r.nswitch)}
        out.log(f"V = {V}: {li.kind}, A = {np.degrees(li.amplitude):.4g} deg, f = {li.frequency:.4g} Hz, "
                f"{r.nswitch} gap crossings")
        _write_history(out, r, "lco", d0, f"Freeplay response at V = {V:.4g} m/s")
        if plots:
            plt = _plt()
            fig, ax = plt.subplots(figsize=(5, 4.2))
            n = r.t.size
            ax.plot(np.degrees(r.delta[n // 2:]), np.degrees(r.etad[n // 2:, r.names.index(info["hinge"])]
                                                            if info["hinge"] in r.names else np.gradient(r.delta[n // 2:], r.t[n // 2:])),
                    lw=0.8)
            ax.axvspan(-np.degrees(d0), np.degrees(d0), color="0.9")
            ax.set_xlabel("hinge rotation [deg]"); ax.set_ylabel("hinge rate [deg/s]"); ax.grid(alpha=0.3)
            ax.set_title(f"Phase plane, V = {V:.4g} m/s (second half)")
            out.savefig(fig, "phase_plane.png", "Hinge phase plane")

    elif mode == "bifurcation":
        Vs = velocity_range(v)
        t0 = time.perf_counter()
        up = bifurcation_sweep(m, Vs, rho, t_end, x0=x0, dt=dt)
        down = bifurcation_sweep(m, Vs[::-1], rho, t_end, x0=x0, dt=dt)
        out.log(f"bifurcation sweeps (up + down, {Vs.size} speeds): {time.perf_counter() - t0:.1f} s")
        a_dn = np.asarray(down["amp"])[::-1]
        f_dn = np.asarray(down["freq"])[::-1]
        write_curves(out, "bifurcation", Vs,
                     {"A [deg]:sweep up": np.degrees(up["amp"]), "A [deg]:sweep down": np.degrees(a_dn),
                      "f [Hz]:sweep up": up["freq"], "f [Hz]:sweep down": f_dn},
                     "Bifurcation diagram (time marching, continuation in V)", "V [m/s]", "A [deg]")
        lco = [float(V) for V, k in zip(Vs, up["kind"]) if k == "lco"]
        out.summary["bifurcation"] = {"first_lco_up": lco[0] if lco else None,
                                      "kinds_up": list(up["kind"]), "kinds_down": list(down["kind"])[::-1]}
    else:
        raise ValueError(f"unknown freeplay mode {mode!r}")


def _plot_df(out: Output, df, d0, rows):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for idx in df.branches():
        s = df.stable[idx]
        A = df.A[idx] / d0
        for a, w, key in ((ax[0], A, "A"), (ax[1], df.freq[idx], "f")):
            a.plot(df.V[idx][s], w[s], "o", ms=2.5, color="tab:blue")
            a.plot(df.V[idx][~s], w[~s], "o", ms=2.0, mfc="none", color="tab:red")
    for r in rows:
        if r["A_tm_deg"] is not None:
            ax[0].plot(r["V"], np.radians(r["A_tm_deg"]) / d0, "k^", ms=7)
            if r["f_tm_hz"] is not None:
                ax[1].plot(r["V"], r["f_tm_hz"], "k^", ms=7)
    for a in ax:
        for Vx, c in ((df.V_f0, "tab:red"), (df.V_fK, "tab:blue")):
            if np.isfinite(Vx):
                a.axvline(Vx, color=c, ls=":", lw=1)
        a.grid(alpha=0.3); a.set_xlabel("V [m/s]")
    ax[0].set_ylabel("hinge amplitude A / delta0"); ax[1].set_ylabel("LCO frequency [Hz]")
    ax[0].plot([], [], "o", color="tab:blue", label="DF stable"); ax[0].plot([], [], "o", mfc="none", color="tab:red",
                                                                           label="DF unstable")
    ax[0].plot([], [], "k^", label="time marching"); ax[0].legend(fontsize=8)
    fig.suptitle("Freeplay LCO: describing function vs time marching")
    out.savefig(fig, "lco_bifurcation.png", "LCO amplitude / frequency vs speed")


# ======================================================================
# ASE
# ======================================================================
def _ase_model(out: Output, v: dict):
    from ..ase import ASEModel, Actuator, Sensor, gain
    sysm, rho, info = _plant(v, "control")
    sensors = [str(s) for s in _as_list(v["sensors"])]
    K = np.atleast_1d(np.asarray(v["feedback_gains"], float))
    if K.size != len(sensors):
        raise ValueError(f"feedback_gains has {K.size} values for {len(sensors)} sensors ({', '.join(sensors)})")
    wa = 2 * np.pi * float(v["actuator_frequency"])
    za = float(v.get("actuator_damping", 0.7))
    if v["structure"] == "section":
        act = Actuator.second_order("beta", wa, za)
        lib = {"pitch": Sensor("pitch", [0.0, 1.0], "disp"), "pitch_rate": Sensor("pitch_rate", [0.0, 1.0], "rate"),
               "plunge": Sensor("plunge", [1.0, 0.0], "disp"), "plunge_rate": Sensor("plunge_rate", [1.0, 0.0], "rate")}
        cmd = "beta_cmd"
    else:
        from ..ase import beam_sensor
        from ..loads import MonitorStation, modal_load_matrices, station_output_group
        from ..models.goland import SPAN, X_EA
        g = info["g"]
        st = [MonitorStation("root", [0.0, 0.01, 0.0], [0.0, 1.0, 0.0], ("Fz", "Mx", "My"))]
        station_output_group(sysm, g["table"], g["mesh"], g["modal"], st, "fsm")
        Tk = modal_load_matrices(g["modal"], st[0])[1][:, :sysm.ns]
        act = Actuator.second_order("flap", wa, za)
        lib = {"root_bending": Sensor("root_bending", Tk[1], "disp"),
               "tip_accel": beam_sensor(g["modal"], X_EA, 0.95 * SPAN, X_EA, "acc", "tip_accel")}
        cmd = "flap_cmd"
    bad = [s for s in sensors if s not in lib]
    if bad:
        raise ValueError(f"unknown sensor(s) {bad} for structure {v['structure']}; choose from {list(lib)}")
    ase = ASEModel(sysm, [act], [lib[s] for s in sensors], rho)
    Kc = gain(K[None, :], sensors, cmd)
    out.log(f"plant: {info['name']}; actuator {v['actuator_frequency']} Hz / zeta {za}; "
            f"{cmd} = K y, y = [{', '.join(sensors)}], K = {K.tolist()}")
    out.summary.update({"plant": info["name"], "rho": rho, "sensors": sensors, "feedback_gains": K.tolist(),
                        "actuator_hz": float(v["actuator_frequency"])})
    return ase, Kc, sysm, rho, info


def route_ase(out: Output, v: dict, base):
    v = _defaults(out, v, "ase")
    if v["structure"] == "section" and v["mode"] == "gla":
        raise ValueError("mode gla needs gust columns: use structure goland_flap")
    ase, Kc, sysm, rho, info = _ase_model(out, v)
    mode = v["mode"]
    Vs = velocity_range(v)
    plots = v.get("plots", True)

    if mode == "flutter_suppression":
        t0 = time.perf_counter()
        r_ol = ase.flutter(Vs, None)
        r_cl = ase.flutter(Vs, Kc)
        Vlo, Vhi = float(Vs[0]), float(Vs[-1])
        vf_ol = ase.flutter_speed((Vlo, Vhi), None, n=97)
        vf_cl = ase.flutter_speed((Vlo, Vhi), Kc, n=97)
        out.log(f"open / closed loop p-method: {time.perf_counter() - t0:.1f} s")
        write_flutter(out, r_ol, None, "open", "Open loop (actuator locked)")
        write_flutter(out, r_cl, None, "closed", "Closed loop")
        inc = vf_cl / vf_ol - 1 if np.isfinite(vf_ol) and np.isfinite(vf_cl) else None
        out.summary.update({"V_f_open": _num(vf_ol), "V_f_closed": _num(vf_cl), "increase": _num(inc)})
        out.log(f"V_f open loop = {vf_ol:.4g} m/s, closed loop = {vf_cl:.4g} m/s"
                + (f" ({100 * inc:+.1f} %)" if inc is not None else ""))
        if plots:
            _plot_rootlocus(out, r_ol, r_cl)

    elif mode == "margins":
        V = float(v["velocity"])
        w, L, _ = ase.loop(V, Kc)
        from ..ase import bode_data
        bd = bode_data(L, w)
        mg = ase.margins(V, Kc)
        out.log(mg.summary())
        out.summary["margins"] = {"V": V, "P": int(mg.P), "N_ccw": int(mg.N_ccw),
                                  "nyquist_stable": bool(mg.nyquist_stable),
                                  "closed_loop_stable": None if mg.closed_loop_stable is None
                                  else bool(mg.closed_loop_stable),
                                  "gm_upper_db": _num(mg.gm_upper_db), "gm_lower_db": _num(mg.gm_lower_db),
                                  "pm_deg": _num(mg.pm), "w_pm": _num(mg.w_pm)}
        f = np.asarray(w) / (2 * np.pi)
        write_curves(out, "bode", f, {"|L| [dB]:loop": bd["mag_db"], "phase [deg]:loop": bd["phase_deg"]},
                     f"Loop transfer L(i w) at V = {V:.4g} m/s", "f [Hz]", "", meta=out.summary["margins"])
        # margins versus speed
        rows = []
        for Vm in np.linspace(Vs[0], Vs[-1], 9)[1:]:
            mm = ase.margins(float(Vm), Kc)
            rows.append({"V": _num(Vm), "gm_upper_db": _num(mm.gm_upper_db), "gm_lower_db": _num(mm.gm_lower_db),
                         "pm_deg": _num(mm.pm), "nyquist_stable": bool(mm.nyquist_stable)})
        out.summary["margins_vs_V"] = rows
        if plots:
            plt = _plt()
            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            ax[0].plot(L.real, L.imag, lw=1.1); ax[0].plot(L.real, -L.imag, lw=0.8, ls="--", color="0.5")
            ax[0].plot([-1], [0], "r+", ms=12, mew=2)
            ax[0].set_xlim(-3, 3); ax[0].set_ylim(-3, 3); ax[0].set_aspect("equal"); ax[0].grid(alpha=0.3)
            ax[0].set_title(f"Nyquist (P = {mg.P}, N_ccw = {mg.N_ccw})")
            ax[1].semilogx(f, bd["mag_db"], lw=1.2); ax[1].set_ylabel("|L| [dB]"); ax[1].grid(alpha=0.3)
            ax2 = ax[1].twinx(); ax2.semilogx(f, bd["phase_deg"], lw=0.9, color="tab:orange")
            ax2.set_ylabel("phase [deg]"); ax[1].set_xlabel("f [Hz]"); ax[1].set_title("Bode")
            fig.suptitle(f"Loop at V = {V:.4g} m/s")
            out.savefig(fig, "nyquist_bode.png", "Nyquist / Bode of the loop")

    elif mode == "gla":
        V = float(v["velocity"])
        Uds = float(v.get("gust_amplitude", 5.0))
        Hs = velocity_range(v, "gust_gradients") if "gust_gradients" in v else np.array([5., 15., 30.])
        t_end, dt = float(v.get("t_end", 1.0)), float(v.get("dt", 5e-4))
        outs = ["root.Mx", "root.Fz"]
        rows = []
        for i, H in enumerate(Hs):
            r0 = ase.gust_response(V, float(H), Uds, None, t_end=t_end, dt=dt, outputs=outs)
            r1 = ase.gust_response(V, float(H), Uds, Kc, t_end=t_end, dt=dt, outputs=outs)
            p0, p1 = r0["peak"], r1["peak"]
            rows.append({"H": _num(H), "Mx_open": _num(p0["root.Mx"]), "Mx_closed": _num(p1["root.Mx"]),
                         "Mx_reduction": _num(1 - p1["root.Mx"] / p0["root.Mx"]),
                         "Fz_open": _num(p0["root.Fz"]), "Fz_closed": _num(p1["root.Fz"]),
                         "flap_peak_deg": _num(np.degrees(p1["flap"]))})
            out.log(f"  H = {H:6.2f} m: root Mx {p0['root.Mx']:.4g} -> {p1['root.Mx']:.4g} N m "
                    f"({100 * rows[-1]['Mx_reduction']:.1f} % less), flap peak {rows[-1]['flap_peak_deg']} deg")
            if i == len(Hs) // 2 or len(Hs) == 1:
                y0 = r0["y"][:, r0["outputs"].index("root.Mx")]
                y1 = r1["y"][:, r1["outputs"].index("root.Mx")]
                write_timehist(out, "gla", r0["t"], {"root Mx [N m]:open loop": y0, "root Mx [N m]:closed loop": y1,
                                                     "flap [deg]:closed loop":
                                                         np.degrees(r1["y"][:, r1["outputs"].index("flap")]),
                                                     "w_g [m/s]": r0["wg"]},
                               f"1-cos gust H = {H:.3g} m, Uds = {Uds} m/s, V = {V} m/s", meta={"H": H, "V": V})
        out.summary["gla"] = {"V": V, "Uds": Uds, "rows": rows}
        vf_ol = ase.flutter_speed((float(Vs[0]), float(Vs[-1])), None, n=85)
        vf_cl = ase.flutter_speed((float(Vs[0]), float(Vs[-1])), Kc, n=85)
        out.summary.update({"V_f_open": _num(vf_ol), "V_f_closed": _num(vf_cl)})
        out.log(f"flutter speed open / closed loop: {vf_ol:.4g} / {vf_cl:.4g} m/s")
    else:
        raise ValueError(f"unknown ase mode {mode!r}")


def _plot_rootlocus(out: Output, r_ol, r_cl):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for r, c, lab in ((r_ol, "0.4", "open loop"), (r_cl, "tab:red", "closed loop")):
        V = np.asarray(r.V)
        for b in range(np.asarray(r.g).shape[1]):
            Vb = V[:, b] if V.ndim == 2 else V
            ax[0].plot(Vb, np.asarray(r.g)[:, b], color=c, lw=1.1, label=lab if b == 0 else None)
            w = 2 * np.pi * np.asarray(r.freq)[:, b]
            lam = np.asarray(r.g)[:, b] * 0.5 * w
            ok = w > 0.5                           # drop aperiodic (divergence / lag) roots
            ax[1].plot(np.where(ok, lam, np.nan), np.where(ok, w, np.nan), color=c, lw=1.0)
    ax[0].axhline(0, color="k", lw=0.7); ax[0].set_ylim(-1.0, 0.4)
    ax[0].set_xlabel("V [m/s]"); ax[0].set_ylabel("damping g"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].axvline(0, color="k", lw=0.7); ax[1].set_xlabel("Re(lambda) ~ g w / 2"); ax[1].set_ylabel("Im(lambda) [rad/s]")
    ax[1].grid(alpha=0.3); ax[1].set_title("root locus with V")
    fig.suptitle("Active flutter suppression: open vs closed loop")
    out.savefig(fig, "ase_vg_rootlocus.png", "Open- vs closed-loop V-g and root locus")
