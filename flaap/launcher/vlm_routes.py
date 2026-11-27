"""``model uvlm``: steady VLM and time-domain UVLM (lift, drag, moment histories).

Outputs (flaap-txt, read by FLAAP Studio)::

    vlm_steady.txt            alpha sweep: CL, CDi, Cm, span efficiency, span loads
    curves_vlm_lift.txt       CL, Cm versus alpha
    curves_vlm_polar.txt      CDi versus CL (+ elliptic CL^2 / (pi A))
    curves_vlm_spanload.txt   local cl c / c_ref versus y
    timehist_uvlm.txt         CL, CD, Cm, L, D, My (+ motion / gust) versus t
    curves_uvlm_spanload.txt  span load at the end of the run
    uvlm_wake.txt             final wake node grids (one table per surface)
    flutter_pk.txt            mode flutter: PK on linearised-UVLM GAFs
    curves_uvlm_ss_*.txt      mode flutter: discrete-time state-space damping / frequency
"""
from __future__ import annotations

import time
from types import SimpleNamespace

import numpy as np

from .builder import build_case, density, surfaces_from_values, velocity_range
from .runner import (Output, _num, _plt, gaf_table, jlist, write_aero_mesh, write_curves, write_flutter,
                     write_modes, write_timehist)


def _geometry(out: Output, v: dict, base):
    """Mesh from ``surfaceN`` keys (no structure needed) or from the structure preset."""
    from ..geometry import build_mesh
    surfs = surfaces_from_values(v)
    if surfs is not None:
        mesh = build_mesh(surfs, v.get("symmetry", "symmetric"))
        bref = float(v["bref"]) if "bref" in v else 0.25 * float(np.mean([s.root_chord + s.tip_chord
                                                                          for s in surfs]))
    else:
        case = build_case(v, base, need_aero=True, log=out.log)
        mesh, bref = case.mesh, case.bref
    write_aero_mesh(out, SimpleNamespace(mesh=mesh, bref=bref, mach=float(v.get("mach", 0.0))))
    return mesh, bref


def _reference(mesh, v):
    """S, c_ref (= S / span), x_ref (default: 1/4 of c_ref behind the root leading edge)."""
    strips = [idx for _, idx in mesh.strips()]
    span = float(sum(2 * mesh.e[idx[0]] for idx in strips))
    S = float(mesh.total_area)
    c_ref = S / span
    root = int(np.argmin(np.abs(mesh.colloc[:, 1])))
    x_le = float(mesh.corners[mesh.strip_id == mesh.strip_id[root]][:, :, 0].min())
    x_ref = float(v["x_ref"]) if "x_ref" in v else x_le + 0.25 * c_ref
    return S, c_ref, x_ref


def route_uvlm(out: Output, v: dict, base):
    mode = v["mode"]
    if mode == "flutter":
        return _flutter(out, v, base)
    mesh, bref = _geometry(out, v, base)
    S, c_ref, x_ref = _reference(mesh, v)
    out.summary.update({"s_ref": _num(S), "c_ref": _num(c_ref), "x_ref": _num(x_ref), "n_panels": int(mesh.n)})
    if mode == "steady":
        return _steady(out, v, mesh, x_ref, c_ref)
    return _unsteady(out, v, mesh, bref, S, c_ref, x_ref)


# ----------------------------------------------------------------- steady VLM
def _steady(out, v, mesh, x_ref, c_ref):
    from ..aero.uvlm import steady_vlm
    a = v.get("alphas", [-4.0, 10.0, 8])
    alphas = np.linspace(float(a[0]), float(a[1]), int(a[2]))
    mach = float(v.get("mach", 0.0))
    t0 = time.perf_counter()
    r = steady_vlm(mesh, alphas, mach, x_ref=x_ref, c_ref=c_ref)
    out.log(f"steady VLM: {mesh.n} panels, {alphas.size} incidences, M = {mach} ({time.perf_counter() - t0:.2f} s)")
    out.log(f"  CL_alpha = {r.CL_alpha:.4f} /rad ({np.radians(r.CL_alpha):.5f} /deg), A = {r.aspect_ratio:.3f}, "
            f"e = {r.e:.4f}, x_ac = {r.x_ac:.4f} m")
    out.data("vlm_steady", {"alpha_deg": jlist(alphas), "CL": jlist(r.CL), "CDi": jlist(r.CDi), "Cm": jlist(r.Cm),
                            "CL_alpha": _num(r.CL_alpha), "span_efficiency": _num(r.e),
                            "aspect_ratio": _num(r.aspect_ratio), "s_ref": _num(r.s_ref), "c_ref": _num(r.c_ref),
                            "x_ref": _num(r.x_ref), "x_ac": _num(r.x_ac), "mach": mach,
                            "y_strip": jlist(r.y_strip), "cl_c": jlist(r.cl_c), "cl_local": jlist(r.cl_local)},
             "vlm_steady", "Steady VLM polar and span loads")
    out.summary.update({"CL_alpha_per_rad": _num(r.CL_alpha), "CL_alpha_per_deg": _num(np.radians(r.CL_alpha)),
                        "span_efficiency": _num(r.e), "aspect_ratio": _num(r.aspect_ratio), "x_ac": _num(r.x_ac),
                        "CL_max_alpha": _num(r.CL[-1]), "CDi_max_alpha": _num(r.CDi[-1])})
    if v.get("plots", True):
        write_curves(out, "vlm_lift", alphas, {"CL": r.CL, "Cm": r.Cm}, "VLM lift and moment", "alpha [deg]",
                     "coefficient")
        ell = r.CL ** 2 / (np.pi * r.aspect_ratio)
        o = np.argsort(r.CL)
        write_curves(out, "vlm_polar", r.CL[o], {"CDi VLM (Trefftz)": r.CDi[o], "CL^2/(pi A) (e = 1)": ell[o]},
                     "Induced drag polar", "CL", "CDi")
        sel = np.unique(np.linspace(0, alphas.size - 1, min(alphas.size, 5)).round().astype(int))
        oy = np.argsort(r.y_strip)
        write_curves(out, "vlm_spanload", r.y_strip[oy],
                     {f"alpha {alphas[i]:g} deg": r.cl_c[i][oy] for i in sel if abs(alphas[i]) > 1e-12},
                     "Span load cl c / c_ref", "y [m]", "cl c / c_ref")


# ----------------------------------------------------------------- time-domain UVLM
def _unsteady(out, v, mesh, bref, S, c_ref, x_ref):
    from ..aero.uvlm import UVLM, Kinematics, harmonic, one_minus_cos_gust, sharp_edged_gust, steady_vlm
    mode = v["mode"]
    U = float(v["velocity"])
    rho = density(v)
    u = UVLM(mesh)
    dt = float(v["dt"]) if "dt" in v else u.default_dt(U) / max(1, int(v.get("chord_steps", 1)))
    alpha = float(v.get("alpha", 0.0))
    kin = Kinematics(alpha0=np.radians(alpha), x_ref=x_ref)
    x_le = float(mesh.corners[:, :, 0].min())
    meta = {"mode": mode, "velocity": U, "density": rho, "alpha_deg": alpha, "x_ref": x_ref, "c_ref": c_ref,
            "s_ref": S, "dt": dt, "free_wake": bool(v.get("free_wake", False))}
    if mode == "impulsive":
        t_end = float(v.get("t_end", 20.0 * c_ref / U))
    elif mode == "harmonic":
        k = float(v["reduced_frequency"])
        w = k * U / bref
        th0 = np.radians(float(v.get("pitch_amplitude_deg", 0.0)))
        h0 = float(v.get("plunge_amplitude", 0.0))
        ph = np.radians(float(v.get("phase_deg", 90.0)))
        if th0:
            kin.pitch = harmonic(th0, w, ph)
        if h0:
            kin.plunge = harmonic(h0, w)
        T = 2 * np.pi / w
        t_end = float(v.get("t_end", float(v.get("ncycles", 4.0)) * T))
        meta.update({"reduced_frequency": k, "omega": w, "bref": bref, "pitch_amplitude_deg": np.degrees(th0),
                     "plunge_amplitude": h0, "phase_deg": np.degrees(ph)})
    else:   # gust
        w0 = float(v["gust_velocity"])
        shape = v.get("gust_shape", "1-cos")
        if shape == "sharp":
            kin.gust = sharp_edged_gust(w0, U, x_le)
            t_end = float(v.get("t_end", 30.0 * c_ref / U))
        else:
            H = float(v["gust_length"])
            kin.gust = one_minus_cos_gust(w0, H, U, x_le)
            t_end = float(v.get("t_end", (2 * H + 10.0 * c_ref) / U))
            meta["gust_length"] = H
        meta.update({"gust_shape": shape, "gust_velocity": w0})
    nsteps = int(round(t_end / dt))
    wl = v.get("wake_length")
    out.log(f"UVLM {mode}: {mesh.n} rings, V = {U} m/s, dt = {dt:.4g} s, {nsteps} steps, "
            f"{'free' if v.get('free_wake') else 'prescribed'} wake" + (f" ({wl} chords)" if wl else ""))
    if nsteps > 4000:
        out.log(f"  warning: {nsteps} time steps - consider t_end / wake_length to limit the run time")
    t0 = time.perf_counter()
    r = u.simulate(U, t_end, kin, rho=rho, dt=dt, free_wake=bool(v.get("free_wake", False)),
                   wake_length=None if wl is None else float(wl), x_ref=x_ref,
                   verbose=bool(v.get("verbose", True)))
    out.log(f"UVLM done in {time.perf_counter() - t0:.1f} s")
    series = {"CL": r.CL, "CD": r.CD, "Cm": r.Cm, "L": r.L, "D": r.D, "My": r.My}
    base_run = None
    if mode == "gust" and abs(alpha) > 0:
        # the impulsive start at alpha != 0 is removed with a gust-free run (same wake, same steps)
        out.log("  baseline run without gust (alpha != 0) for the gust increments")
        kb = Kinematics(alpha0=np.radians(alpha), x_ref=x_ref)
        base_run = u.simulate(U, t_end, kb, rho=rho, dt=dt, free_wake=bool(v.get("free_wake", False)),
                              wake_length=None if wl is None else float(wl), x_ref=x_ref)
        series = {"CL": r.CL, "dCL_gust": r.CL - base_run.CL, "CD": r.CD, "Cm": r.Cm,
                  "dL_gust": r.L - base_run.L, "dMy_gust": r.My - base_run.My, "L": r.L, "D": r.D, "My": r.My}
    units = {"CL": "-", "CD": "-", "Cm": "-", "L": "N", "D": "N", "My": "N m", "t": "s"}
    if kin.pitch is not None:
        series["pitch_deg"] = np.degrees([kin.pitch(t)[0] for t in r.t])
        units["pitch_deg"] = "deg"
    if kin.plunge is not None:
        series["plunge"] = np.array([kin.plunge(t)[0] for t in r.t])
        units["plunge"] = "m"
    if kin.gust is not None:
        series["gust_at_LE"] = np.array([float(kin.gust(np.array([x_le]), t)[0]) for t in r.t])
        units["gust_at_LE"] = "m/s"
    write_timehist(out, "uvlm", r.t, series, f"UVLM {mode}: lift, drag, moment", meta, units)
    oy = np.argsort(r.y_strip)
    out.summary.update({"velocity": U, "density": _num(rho), "dt": _num(dt), "nsteps": int(nsteps),
                        "CL_final": _num(r.CL[-1]), "CD_final": _num(r.CD[-1]), "Cm_final": _num(r.Cm[-1])})
    if mode == "impulsive" and abs(alpha) > 0:
        st = steady_vlm(mesh, [alpha], 0.0, x_ref=x_ref, c_ref=c_ref)
        out.summary.update({"CL_steady_vlm": _num(st.CL[0]), "CDi_steady_vlm": _num(st.CDi[0]),
                            "CL_final_over_steady": _num(r.CL[-1] / st.CL[0])})
        out.log(f"  CL(t_end) = {r.CL[-1]:.4f}, steady VLM CL = {st.CL[0]:.4f} "
                f"({100 * r.CL[-1] / st.CL[0]:.1f} %); CD(t_end) = {r.CD[-1]:.5f}, steady CDi = {st.CDi[0]:.5f}")
    if mode == "harmonic":
        w = meta["omega"]
        sel = r.t >= r.t[-1] - 2 * np.pi / w + 1e-12
        A = np.column_stack([np.ones(sel.sum()), np.cos(w * r.t[sel]), np.sin(w * r.t[sel])])
        res = {}
        for key in ("CL", "CD", "Cm"):
            c = np.linalg.lstsq(A, series[key][sel], rcond=None)[0]
            res[key] = (c[0], complex(c[1], -c[2]))
        out.summary.update({"CL_mean": _num(res["CL"][0]), "CL_1st_harmonic_amplitude": _num(abs(res["CL"][1])),
                            "CL_1st_harmonic_phase_deg": _num(np.degrees(np.angle(res["CL"][1]))),
                            "CD_mean_last_cycle": _num(res["CD"][0]), "CT_mean_last_cycle": _num(-res["CD"][0]),
                            "Cm_1st_harmonic_amplitude": _num(abs(res["Cm"][1]))})
        out.log(f"  last cycle: CL1 = {abs(res['CL'][1]):.4f} at {np.degrees(np.angle(res['CL'][1])):.1f} deg "
                f"(phase vs sin(w t)), mean CD = {res['CD'][0]:.3e} (thrust CT = {-res['CD'][0]:.3e})")
    if mode == "gust":
        dCL = r.CL if base_run is None else r.CL - base_run.CL
        dL = r.L if base_run is None else r.L - base_run.L
        dMy = r.My if base_run is None else r.My - base_run.My
        i = int(np.argmax(np.abs(dCL)))
        out.summary.update({"dCL_gust_peak": _num(dCL[i]), "t_peak": _num(r.t[i]), "dL_gust_peak": _num(dL[i]),
                            "dMy_gust_peak": _num(dMy[i]), "CL_at_peak": _num(r.CL[i])})
        out.log(f"  gust peak: dCL = {dCL[i]:.4f} at t = {r.t[i]:.4f} s (dL = {dL[i]:.1f} N)")
    out.data("uvlm_wake", {f"surface{i + 1}": jlist(W.reshape(-1, 3), 6) for i, W in enumerate(r.wake)} |
             {"shape": [list(W.shape) for W in r.wake], "free_wake": bool(v.get("free_wake", False))},
             "uvlm_wake", "Final wake node grids")
    if v.get("plots", True):
        write_curves(out, "uvlm_spanload", r.y_strip[oy], {"cl c / c_ref (t_end)": r.cl_strip[-1][oy]},
                     "UVLM span load at t_end", "y [m]", "cl c / c_ref")
        _plot_wake(out, mesh, r)


def _plot_wake(out, mesh, r):
    plt = _plt()
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for W in r.wake:
        for j in range(W.shape[1]):
            ax[0].plot(W[:, j, 0], W[:, j, 2], "b-", lw=0.4)
            ax[1].plot(W[:, j, 0], W[:, j, 1], "b-", lw=0.4)
    for c in mesh.corners:
        q = np.vstack([c, c[:1]])
        ax[0].plot(q[:, 0], q[:, 2], "k-", lw=0.6)
        ax[1].plot(q[:, 0], q[:, 1], "k-", lw=0.6)
    ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("z [m]"); ax[0].set_title("side view")
    ax[1].set_xlabel("x [m]"); ax[1].set_ylabel("y [m]"); ax[1].set_title("top view"); ax[1].set_aspect("equal")
    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle("UVLM wake at t_end" + (" (free wake)" if r.meta.get("free_wake") else ""))
    out.savefig(fig, "uvlm_wake.png", "UVLM wake geometry")


# ----------------------------------------------------------------- flutter
def _flutter(out, v, base):
    from ..aero.uvlm import flutter_crossing, flutter_eigen_sweep
    from ..flutter import identify, pk_method
    v = dict(v, aero_method="uvlm")
    case = build_case(v, base)
    case.values = v
    write_modes(out, case)
    write_aero_mesh(out, case)
    tab = gaf_table(out, case)
    rho = density(v)
    V = velocity_range(v, default=(60.0, 260.0, 41))
    out.summary["density"] = rho
    t0 = time.perf_counter()
    r = pk_method(case.M, case.K, tab, rho, V, C=case.C())
    identify(r)
    out.log(f"PK on linearised-UVLM GAFs: {time.perf_counter() - t0:.1f} s")
    write_flutter(out, r, case, "pk", "PK method (linearised UVLM GAFs)")
    lin = case.extra["linear_uvlm"]
    t0 = time.perf_counter()
    Vs, lam = flutter_eigen_sweep(lin, rho, case.M, case.K, V, C=case.C())
    out.log(f"discrete-time UVLM state space ({lin.N + lin.Nw + 3 * case.M.shape[0]} states): "
            f"{V.size} speeds in {time.perf_counter() - t0:.1f} s")
    zeta = -lam.real / np.abs(lam)
    fr = lam.imag / (2 * np.pi)
    names = case.names
    write_curves(out, "uvlm_ss_damping", Vs, {names[j]: zeta[:, j] for j in range(lam.shape[1])},
                 "UVLM state space: damping ratio", "V [m/s]", "zeta = -Re(lambda)/|lambda|")
    write_curves(out, "uvlm_ss_frequency", Vs, {names[j]: fr[:, j] for j in range(lam.shape[1])},
                 "UVLM state space: frequency", "V [m/s]", "f [Hz]")
    cr = flutter_crossing(Vs, lam)
    out.summary["flutter_uvlm_state_space"] = None if cr is None else {
        "velocity": _num(cr[0]), "freq_hz": _num(cr[1] / (2 * np.pi)), "mode": int(cr[2]) + 1}
    out.log("state-space flutter: " + ("none in the speed range" if cr is None else
                                       f"V = {cr[0]:.2f} m/s, f = {cr[1] / (2 * np.pi):.3f} Hz"))
