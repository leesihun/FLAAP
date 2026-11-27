"""Execute a resolved flat config and write the FLAAP output contract.

Output directory layout (read by FLAAP Studio, see ``docs/OUTPUT_FORMAT.md``).
Every result file is plain text in the config grammar (``key value`` lines and
``table ... end`` blocks, :mod:`flaap.io.results_txt`)::

    results.txt             title, model, mode, status, elapsed, summary, artifacts
    config_resolved.txt     every key incl. defaults
    log.txt                 run log
    modes.txt               modal table, nodes, connectivity, mode shapes, aero-box shapes
    aero_mesh.txt           aerodynamic boxes (corners, surfaces, control surfaces)
    cp.txt                  dCp per mode (and gust) at each reduced frequency
    flutter_<method>.txt    identified points and V-g / V-f per branch
    timehist_<name>.txt     time histories (gust, manoeuvre, freeplay, ASE)
    curves_<name>.txt       generic x-y curves (gust envelope, PSD, LCO amplitude, root locus)
    tsd_steady.txt          steady TSD Cp on the wing grid
    vlm_steady.txt          steady VLM polar and span loads (model uvlm)
    figures/*.png           matplotlib figures

``output_json true`` additionally writes the same content as ``*.json``.
"""
from __future__ import annotations

import json
import math
import time
import traceback
from pathlib import Path

import numpy as np

from ..io.results_txt import write_result
from .builder import (Case, build_case, density, k_values, tsd_alpha, tsd_camber, tsd_lifting, tsd_thickness,
                      velocity_range)
from .parser import write_config


# ----------------------------------------------------------------- number helpers
def _num(x, sig=7):
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    x = float(x)
    return float(f"{x:.{sig}g}") if math.isfinite(x) else None


def jlist(a, sig=7):
    """numpy array -> nested list with NaN/inf -> null and limited digits."""
    a = np.asarray(a)
    if a.dtype.kind in "iub":
        return a.tolist()
    if a.ndim == 0:
        return _num(a, sig)
    return [jlist(r, sig) for r in a] if a.ndim > 1 else [_num(v, sig) for v in a.tolist()]


class Output:
    def __init__(self, directory: Path, verbose=True, also_json=False):
        self.dir = Path(directory)
        self.also_json = also_json
        self.fig = self.dir / "figures"
        self.fig.mkdir(parents=True, exist_ok=True)
        self.artifacts: list[dict] = []
        self.summary: dict = {}
        self.verbose = verbose
        self._log = open(self.dir / "log.txt", "w", encoding="utf-8")

    def log(self, *msg):
        s = " ".join(str(m) for m in msg)
        self._log.write(s + "\n")
        self._log.flush()
        if self.verbose:
            print(s, flush=True)

    def data(self, name: str, data: dict, kind: str, title: str):
        """Write the result file ``<name>.txt`` (contract dictionary ``data`` of ``kind``)."""
        stem = name[:-5] if name.endswith(".json") else name[:-4] if name.endswith(".txt") else name
        write_result(self.dir / f"{stem}.txt", kind, data, title)
        self.artifacts.append({"file": f"{stem}.txt", "kind": kind, "title": title})
        if self.also_json:
            (self.dir / f"{stem}.json").write_text(json.dumps(data, separators=(",", ":"), allow_nan=False),
                                                    encoding="utf-8")

    def file(self, name: str, kind: str, title: str):
        self.artifacts.append({"file": name, "kind": kind, "title": title})

    def savefig(self, fig, name: str, title: str):
        import matplotlib.pyplot as plt
        fig.savefig(self.fig / name, dpi=110, bbox_inches="tight")
        plt.close(fig)
        self.file(f"figures/{name}", "figure", title)

    def close(self):
        self._log.close()


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# ----------------------------------------------------------------- writers
def _modal_arrays(modal):
    """(nodes, shapes (Nn,6,nm), freq_hz, gen_mass, labels) for ModalModel or FlexibleAircraft."""
    if hasattr(modal, "freq_hz"):
        return modal.nodes, modal.shapes, np.asarray(modal.freq_hz), np.asarray(modal.gen_mass), \
            list(modal.labels or [f"mode {i + 1}" for i in range(modal.nmodes)])
    from scipy.linalg import eigh
    M, K = np.asarray(modal.M), np.asarray(modal.K)
    w2 = np.sort(np.real(eigh(K, M, eigvals_only=True)))
    return modal.nodes, modal.shapes, np.sqrt(np.maximum(w2, 0)) / (2 * np.pi), np.diag(M), list(modal.names)


def aero_corner_shapes(case: Case) -> np.ndarray:
    """(nm, N, 4) normal displacement of each aero mode at the box corners
    (linear chordwise extrapolation from the force and collocation points)."""
    m, md = case.mesh, case.modes
    xf, xc = m.force_point[:, 0], m.colloc[:, 0]
    slope = (md.h_c - md.h_f) / (xc - xf)[:, None]
    cx = m.corners[:, :, 0]                                    # (N,4)
    h = md.h_f[:, None, :] + (cx - xf[:, None])[:, :, None] * slope[:, None, :]
    h = np.transpose(h, (2, 0, 1))                              # (nm, N, 4)
    # boxes sharing a corner point extrapolate independently -> average them so
    # the deformed surface is continuous (no gaps between box strips)
    pts = m.corners.reshape(-1, 3)
    tol = 1e-6 * max(np.ptp(pts, axis=0).max(), 1e-12)
    _, inv = np.unique(np.round(pts / tol).astype(np.int64), axis=0, return_inverse=True)
    inv = inv.ravel()
    flat = h.reshape(h.shape[0], -1)
    acc = np.zeros((flat.shape[0], inv.max() + 1), dtype=flat.dtype)
    np.add.at(acc.T, inv, flat.T)
    cnt = np.bincount(inv)
    return (acc / cnt)[:, inv].reshape(h.shape)


def write_modes(out: Output, case: Case):
    nodes, shapes, f, gm, labels = _modal_arrays(case.modal)
    nm = shapes.shape[2]
    u = np.transpose(shapes[:, :3, :], (2, 0, 1)).copy()         # (nm, Nn, 3)
    for i in range(nm):
        s = np.abs(u[i]).max()
        if s > 0:
            u[i] /= s
    data = {"name": getattr(case.modal, "name", "structure"), "labels": labels, "freq_hz": jlist(f),
            "gen_mass": jlist(gm), "nodes": jlist(nodes), "shapes": jlist(u, 5), "structure": case.values.get("structure")}
    for key in ("quads", "tris", "edges"):
        if key in case.struct_vis:
            data[key] = case.struct_vis[key]
    if case.mesh is not None and case.modes is not None:
        h = aero_corner_shapes(case)[: len(labels)]
        for i in range(h.shape[0]):
            s = np.abs(h[i]).max()
            if s > 0:
                h[i] /= s
        data["aero_h"] = jlist(h, 5)
    out.data("modes", data, "modes", "Mode shapes and frequencies")
    lines = ["mode  label                  freq [Hz]      omega [rad/s]   gen. mass"]
    for i, lab in enumerate(labels):
        lines.append(f"{i + 1:4d}  {lab:<20s} {f[i]:12.5f} {2 * np.pi * f[i]:14.5f} {gm[i]:12.5g}")
    out.summary["frequencies_hz"] = jlist(f, 6)
    out.log("\n".join(lines))


def write_aero_mesh(out: Output, case: Case):
    m = case.mesh
    data = {"corners": jlist(m.corners, 6), "surface_id": m.surface_id.tolist(),
            "surfaces": [s.name for s in m.surfaces], "symmetry": m.symmetry,
            "colloc": jlist(m.colloc, 6), "area": jlist(m.area, 6),
            "controls": {n: np.isfinite(x).tolist() for n, x in m.control.items()},
            "bref": case.bref, "mach": case.mach}
    out.data("aero_mesh", data, "aero_mesh", "Aerodynamic boxes")
    out.summary["n_boxes"] = int(m.n)


def write_cp(out: Output, tab, title="Unsteady pressure dCp"):
    if tab.cp is None:
        return
    names = list(tab.names)
    cp = np.transpose(tab.cp, (0, 2, 1))                         # (nk, nm, N)
    if tab.cpg is not None:
        cp = np.concatenate([cp, tab.cpg[:, None, :]], axis=1)
        names.append("gust")
    out.data("cp", {"k": jlist(tab.k), "modes": names, "re": jlist(cp.real, 5), "im": jlist(cp.imag, 5),
                         "mach": tab.mach, "bref": tab.bref}, "cp", title)


def method_key(method: str) -> str:
    """'pk-method' -> 'pk', 'p-method (RFA state space)' -> 'p', 'k-method' -> 'k'."""
    m = str(method).strip().lower()
    for key in ("pk", "k", "p"):
        if m == key or m.startswith(key + "-") or m.startswith(key + " "):
            return key
    return "".join(ch if ch.isalnum() else "_" for ch in m).strip("_") or "flutter"


def write_flutter(out: Output, r, case: Case, key: str | None = None, title: str | None = None):
    key = key or method_key(r.method)
    title = title or r.method
    pts = [{"kind": p.kind, "branch": int(p.branch), "velocity": _num(p.velocity), "freq_hz": _num(p.freq_hz),
            "k": _num(p.k), "density": _num(p.density)} for p in r.points]
    data = {"method": key, "method_title": title, "V": jlist(r.V), "freq": jlist(r.freq), "g": jlist(r.g),
            "k": jlist(r.k), "names": list(r.names), "density": r.density, "bref": r.bref,
            "g_struct": r.g_struct, "points": pts}
    out.data(f"flutter_{key}", data, "flutter", f"{title}: V-g / V-f")
    c = r.critical()
    out.summary[f"flutter_{key}"] = None if c is None else {
        "kind": c.kind, "velocity": _num(c.velocity), "freq_hz": _num(c.freq_hz), "mode": int(c.branch) + 1}
    out.log(r.summary())
    _plt()
    ax = r.plot(title=title)
    out.savefig(ax[0].figure, f"flutter_{key}.png", f"{title}: V-g / V-f")
    d = out.dir / f"flutter_{key}"
    d.mkdir(exist_ok=True)
    r.write_legacy(str(d))


def write_timehist(out: Output, name: str, t, series: dict, title: str, meta=None, units=None):
    out.data(f"timehist_{name}", {"t": jlist(t, 6), "series": {k: jlist(v, 6) for k, v in series.items()},
                                        "meta": meta or {}, "units": units or {}}, "timehist", title)
    plt = _plt()
    keys = list(series)[:6]
    fig, axs = plt.subplots(len(keys), 1, sharex=True, figsize=(7, 1.6 * len(keys) + 0.8), squeeze=False)
    for a, k in zip(axs[:, 0], keys):
        a.plot(t, series[k], lw=1.2)
        a.set_ylabel(k, fontsize=8)
        a.grid(alpha=0.3)
    axs[-1, 0].set_xlabel("t [s]")
    fig.suptitle(title)
    out.savefig(fig, f"timehist_{name}.png", title)


def write_curves(out: Output, name: str, x, curves: dict, title: str, xlabel="", ylabel="", logy=False, meta=None):
    out.data(f"curves_{name}", {"x": jlist(x, 6), "curves": {k: jlist(v, 6) for k, v in curves.items()},
                                      "xlabel": xlabel, "ylabel": ylabel, "logy": logy, "meta": meta or {}},
             "curves", title)
    plt = _plt()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for k, v in curves.items():
        ax.plot(x, v, lw=1.3, label=k)
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.grid(alpha=0.3)
    if len(curves) > 1:
        ax.legend(fontsize=8)
    ax.set_title(title)
    out.savefig(fig, f"curves_{name}.png", title)


def plot_modes(out: Output, case: Case):
    """3-D figure of the first modes on the aero surface (or structure nodes)."""
    plt = _plt()
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    nodes, shapes, f, _, labels = _modal_arrays(case.modal)
    n = min(len(labels), 6)
    ncol = 2 if n in (2, 4) else min(n, 3)
    nrow = (n + ncol - 1) // ncol
    fig = plt.figure(figsize=(4.6 * ncol, 3.6 * nrow))
    h = aero_corner_shapes(case) if case.mesh is not None else None
    span = np.ptp(nodes, axis=0).max() if h is None else np.ptp(case.mesh.corners.reshape(-1, 3), axis=0).max()
    for i in range(n):
        ax = fig.add_subplot(nrow, ncol, i + 1, projection="3d")
        if h is not None:
            hi = h[i] / max(np.abs(h[i]).max(), 1e-30) * 0.15 * span
            c = case.mesh.corners.copy()
            c[:, :, 2] += hi * np.cos(case.mesh.gamma)[:, None]
            c[:, :, 1] -= hi * np.sin(case.mesh.gamma)[:, None]
            col = plt.cm.coolwarm(0.5 + 0.5 * hi.mean(axis=1) / (0.15 * span))
            ax.add_collection3d(Poly3DCollection(c, facecolors=col, edgecolor="k", linewidths=0.15))
            pts = c.reshape(-1, 3)
        else:
            u = shapes[:, :3, i] / max(np.abs(shapes[:, :3, i]).max(), 1e-30) * 0.15 * span
            pts = nodes + u
            ax.scatter(*pts.T, s=4, c=u[:, 2], cmap="coolwarm")
        lo, hi_ = pts.min(0), pts.max(0)
        ext = np.maximum(hi_ - lo, 0.25 * (hi_ - lo).max())
        ctr = 0.5 * (lo + hi_)
        ax.set_xlim(ctr[0] - ext[0] / 2, ctr[0] + ext[0] / 2)
        ax.set_ylim(ctr[1] - ext[1] / 2, ctr[1] + ext[1] / 2)
        ax.set_zlim(ctr[2] - ext[2] / 2, ctr[2] + ext[2] / 2)
        ax.set_box_aspect(tuple(ext))
        ax.set_title(f"{labels[i]}  {f[i]:.3f} Hz", fontsize=9)
        ax.set_axis_off()
        ax.view_init(28, -125)
    fig.tight_layout()
    out.savefig(fig, "modes.png", "Mode shapes")


def plot_cp(out: Output, case: Case, tab, ik=None):
    if tab.cp is None:
        return
    plt = _plt()
    from matplotlib.collections import PolyCollection
    ik = (tab.k.size - 1) // 2 if ik is None else ik
    nm = min(tab.cp.shape[2], 4)
    fig, axs = plt.subplots(2, nm, figsize=(3.2 * nm, 5.8), squeeze=False)
    poly = case.mesh.corners[:, :, :2]
    for j in range(nm):
        for r, part in enumerate(("real", "imag")):
            v = getattr(tab.cp[ik, :, j], part)
            pc = PolyCollection(poly, array=v, cmap="RdBu_r", edgecolor="none")
            lim = np.abs(v).max() or 1.0
            pc.set_clim(-lim, lim)
            a = axs[r, j]
            a.add_collection(pc); a.autoscale(); a.set_aspect("equal")
            a.set_title(f"{tab.names[j]} {'Re' if r == 0 else 'Im'} dCp", fontsize=8)
            fig.colorbar(pc, ax=a, shrink=0.7)
    fig.suptitle(f"dCp at k = {tab.k[ik]:.3f}, M = {tab.mach}")
    fig.tight_layout()
    out.savefig(fig, "cp.png", f"dCp at k = {tab.k[ik]:.3f}")


def plot_gaf(out: Output, tab):
    plt = _plt()
    nm = min(tab.Q.shape[1], 4)
    fig, axs = plt.subplots(nm, nm, figsize=(2.4 * nm, 2.2 * nm), squeeze=False)
    for i in range(nm):
        for j in range(nm):
            a = axs[i, j]
            a.plot(tab.k, tab.Q[:, i, j].real, "b.-", lw=1, ms=3, label="Re")
            a.plot(tab.k, tab.Q[:, i, j].imag, "r.-", lw=1, ms=3, label="Im")
            a.set_title(f"Q[{i + 1},{j + 1}]", fontsize=8); a.tick_params(labelsize=6); a.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=6)
    fig.suptitle("Generalized aerodynamic forces vs k")
    fig.tight_layout()
    out.savefig(fig, "gaf.png", "GAF Q(k)")


# ----------------------------------------------------------------- GAF construction
def gaf_table(out: Output, case: Case, keep_cp=True):
    v = case.values
    kv = k_values(v)
    if v.get("aero_method", "dlm") == "tsd":
        from ..aero.tsd import modes_from_structure, tsd_gaf_table
        if not hasattr(case.modal, "freq_hz"):
            raise ValueError("aero_method tsd needs a clamped modal structure (not demo_aircraft)")
        funcs, names = modes_from_structure(case.modal)
        surf = case.mesh.surfaces[0]
        grid = {k2: int(v[k1]) for k1, k2 in (("tsd_nx", "nchord"), ("tsd_ny", "nspan")) if k1 in v}
        th = tsd_thickness(v.get("tsd_thickness"))
        solver = None
        msg = ""
        if tsd_lifting(v):
            from ..aero.tsd_lifting import LiftingTSDSolver
            a = tsd_alpha(v)
            solver = LiftingTSDSolver(surf, case.mach, th, alpha_deg=a, camber=tsd_camber(v.get("tsd_camber")),
                                      verbose=v.get("verbose", True), **grid)
            msg = f", lifting mean flow alpha = {a:g} deg, camber {v.get('tsd_camber', 'flat')}"
        out.log(f"TSD GAFs: M = {case.mach}, thickness {v.get('tsd_thickness', 'flat')}{msg}, {kv.size} k values")
        t0 = time.perf_counter()
        tab = tsd_gaf_table(surf, case.mach, kv, case.bref, funcs, case.names, thickness=th, solver=solver,
                            keep_cp=False, verbose=v.get("verbose", True), **grid)
        if solver is not None:
            mf = solver.mean_flow
            out.log(f"TSD lifting mean flow: CL = {mf.CL:.4f}, converged {mf.converged}; "
                    f"root shock x/c upper {mf.shock_upper[0]:.3f}, lower {mf.shock_lower[0]:.3f}")
            out.summary.update({"tsd_alpha_deg": mf.alpha_deg, "tsd_CL": _num(mf.CL),
                                "tsd_mean_converged": bool(mf.converged)})
        out.log(f"TSD GAF table done in {time.perf_counter() - t0:.1f} s")
        return tab
    method = v.get("aero_method", "dlm")
    if method in ("machbox", "piston"):
        from ..aero.supersonic import supersonic_gaf_table
        opts = {}
        if method == "machbox" and "machbox_nx" in v:
            opts["nx"] = int(v["machbox_nx"])
        if method == "piston":
            opts["van_dyke"] = int(v.get("piston_order", 2)) != 1
        t0 = time.perf_counter()
        tab = supersonic_gaf_table(case.mesh, case.modes, case.mach, case.bref, kv, method=method, gust=True,
                                   keep_cp=keep_cp, verbose=False, **opts)
        out.log(f"{'Mach box' if method == 'machbox' else 'piston theory'} GAF table: {kv.size} k values, "
                f"M = {case.mach} ({time.perf_counter() - t0:.1f} s)")
        return tab
    if method == "uvlm":
        from ..aero.uvlm import LinearUVLM
        if case.mach > 0.3:
            out.log(f"  note: the linearised UVLM is incompressible; mach {case.mach} is ignored")
        t0 = time.perf_counter()
        lin = LinearUVLM(case.mesh, case.modes, case.bref, float(v.get("wake_chords", 30.0)))
        tab = lin.gaf_table(kv, gust=True, keep_cp=keep_cp)
        case.extra["linear_uvlm"] = lin
        out.log(f"linearised UVLM GAF table: {kv.size} k values, {case.mesh.n} rings, {lin.nrow} wake rows "
                f"({time.perf_counter() - t0:.1f} s)")
        return tab
    from ..aero.gaf import AeroModel
    am = AeroModel(case.mesh, case.modes, case.mach, case.bref)
    t0 = time.perf_counter()
    tab = am.table(kv, gust=True, keep_cp=keep_cp)
    out.log(f"DLM GAF table: {kv.size} k values, {case.mesh.n} boxes, M = {case.mach} "
            f"({time.perf_counter() - t0:.1f} s)")
    return tab


def _stations(v):
    from ..loads import MonitorStation
    st = []
    for k in sorted(k for k in v if k.startswith("station") and k[7:].isdigit()):
        s = v[k]
        st.append(MonitorStation(str(s[0]), [float(x) for x in s[1:4]], [float(x) for x in s[4:7]],
                                 ("Fz", "Mx", "My")))
    return st


def _system(case: Case, tab, n_ctrl=0):
    from ..aero.rfa import fit_roger
    from ..timedomain import AeroelasticSystem
    rfa = fit_roger(tab, int(case.values.get("nlag", 4)))
    names = case.names + case.controls[:n_ctrl]
    return AeroelasticSystem(case.M, case.K, rfa, C=case.C(), n_ctrl=n_ctrl, names=names)


# ----------------------------------------------------------------- routes
def route_modal(out, v, base):
    case = build_case(v, base, need_aero=any(k.startswith("surface") for k in v) or
                      v.get("structure") in ("goland", "goland_ctrl", "agard445", "plate", "demo_aircraft"))
    write_modes(out, case)
    if case.mesh is not None:
        write_aero_mesh(out, case)
    if v.get("plots", True):
        plot_modes(out, case)
    out.summary["nmodes"] = len(case.names)
    return case


def route_aero(out, v, base):
    case = build_case(v, base, with_controls=True)
    write_modes(out, case)
    write_aero_mesh(out, case)
    if v["mode"] == "steady":
        v = dict(v, k_list=[0.0])
        case.values = v
    tab = gaf_table(out, case)
    tab.save(str(out.dir / "gaf_table.npz")); out.file("gaf_table.npz", "gaf", "GAF table (numpy)")
    tab.write_text(str(out.dir / "gaf_table.txt")); out.file("gaf_table.txt", "text", "GAF table (text)")
    write_cp(out, tab)
    if v.get("plots", True):
        plot_cp(out, case, tab, 0 if tab.k.size == 1 else None)
        if tab.k.size > 1:
            plot_gaf(out, tab)
        plot_modes(out, case)
    out.summary.update({"nk": int(tab.k.size), "mach": case.mach, "bref": case.bref})
    if tab.k[0] == 0.0 and case.mesh is not None and tab.cp is not None:
        out.summary["Q0_diag"] = jlist(np.real(np.diag(tab.Q[0])), 5)


def route_flutter(out, v, base):
    from ..flutter import identify, k_method, pk_method
    case = build_case(v, base)
    write_modes(out, case)
    write_aero_mesh(out, case)
    tab = gaf_table(out, case)
    tab.save(str(out.dir / "gaf_table.npz")); out.file("gaf_table.npz", "gaf", "GAF table (numpy)")
    write_cp(out, tab)
    rho = density(v)
    V = velocity_range(v)
    methods = list(v.get("methods", ["pk"])) if v["mode"] == "all" else [v["mode"]]
    if isinstance(v.get("methods"), str):
        methods = [v["methods"]]
    out.summary["density"] = rho
    for m in methods:
        t0 = time.perf_counter()
        if m == "pk":
            r = pk_method(case.M, case.K, tab, rho, V, C=case.C())
        elif m == "k":
            r = k_method(case.M, case.K, tab, rho)
        elif m == "p":
            r = _system(case, tab).flutter(V, rho)
        else:
            raise ValueError(f"unknown flutter method {m!r}")
        identify(r)
        out.log(f"{m}: {time.perf_counter() - t0:.1f} s")
        write_flutter(out, r, case)
    if v.get("plots", True):
        plot_modes(out, case)
        plot_cp(out, case, tab)


def route_gust(out, v, base):
    from .. import gust as G
    case = build_case(v, base)
    write_modes(out, case)
    write_aero_mesh(out, case)
    tab = gaf_table(out, case)
    sysm = _system(case, tab)
    st = _stations(v)
    if st and getattr(case.modal, "fe", None) is not None:
        from ..loads import station_output_group
        station_output_group(sysm, tab, case.mesh, case.modal, st, "fsm")
    V, rho = float(v["velocity"]), density(v)
    alt = float(v.get("altitude", 0.0))
    ss = sysm.state_space(V, rho)
    outputs = [o for o in ss.outputs if not o.startswith("etad:")]
    out.summary.update({"velocity": V, "density": rho})
    if v["mode"] == "discrete":
        gr = v.get("gust_gradients")
        Hs = np.linspace(*gr[:2], int(gr[2])) if gr else G.gust_gradients()
        fg = float(v.get("fg", 1.0))
        env = G.gust_envelope(sysm, V, rho, Hs, lambda H: G.design_gust_velocity(H, alt, fg), outputs=outputs)
        lab = env["labels"]
        amp = np.maximum(np.abs(env["max"]), np.abs(env["min"]))
        write_curves(out, "gust_envelope", env["H"], {l: amp[:, j] for j, l in enumerate(lab) if not l.startswith("etadd")},
                     "1-cos gust peak response vs gradient H", "H [m]", "peak |y|")
        key = next((l for l in lab if ".F" in l or ".M" in l), lab[0])
        Hc = env["critical_H"][key]
        r = G.discrete_gust(ss, V, Hc, G.design_gust_velocity(Hc, alt, fg))
        series = {"wg [m/s]": r.wg}
        series.update({l: r.y[:, r.labels.index(l)] for l in outputs[:8]})
        write_timehist(out, "gust", r.t, series, f"1-cos gust H = {Hc:.1f} m, U_ds = {r.Uds:.2f} m/s",
                       {"H": Hc, "Uds": r.Uds, "V": V})
        rows = ["label, H_critical [m], max, min"]
        for j, l in enumerate(lab):
            rows.append(f"{l}, {env['critical_H'][l]:.2f}, {env['max'][:, j].max():.6e}, {env['min'][:, j].min():.6e}")
        (out.dir / "gust_envelope.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        out.file("gust_envelope.csv", "text", "Gust envelope")
        out.summary["critical_H"] = {k: _num(x) for k, x in env["critical_H"].items()}
        out.summary["peak"] = {l: _num(amp[:, j].max()) for j, l in enumerate(lab)}
        out.log("\n".join(rows))
    else:
        L = float(v.get("turbulence_length", 762.0))
        res = G.turbulence_from_ss(ss, V, outputs, L=L)
        psd = np.asarray(res["psd"])
        write_curves(out, "turbulence_psd", res["omega"], {l: psd[:, j] for j, l in enumerate(res["labels"])},
                     "Output PSD for unit-sigma von Karman turbulence", "omega [rad/s]", "PSD", logy=True)
        out.summary["A_bar"] = {l: _num(a) for l, a in zip(res["labels"], np.atleast_1d(res["Abar"]))}
        out.summary["N0"] = {l: _num(a) for l, a in zip(res["labels"], np.atleast_1d(res["N0"]))}
        out.log("A_bar: " + json.dumps(out.summary["A_bar"]))


def _demo_parts(v, base, out):
    from ..aero.gaf import AeroModel
    from ..loads import MonitorStation
    if v.get("structure") != "demo_aircraft":
        raise ValueError("trim/static and maneuver need structure demo_aircraft (free-flying model)")
    case = build_case(v, base, with_controls=True)
    write_modes(out, case)
    write_aero_mesh(out, case)
    tab = AeroModel(case.mesh, case.modes, case.mach, case.bref).table(k_values(v))
    write_cp(out, tab)
    st = _stations(v) or [MonitorStation("root", [0, 0.05, 0], [0, 1, 0], ("Fz", "Mx", "My"))]
    return case, tab, st


def route_trim(out, v, base):
    from ..trim import G0, StaticAeroelastic, load_rows_from_table
    from ..loads import modal_load_matrices
    if v["mode"] == "divergence":
        from ..aero.gaf import AeroModel
        case = build_case(v, base)
        write_modes(out, case)
        write_aero_mesh(out, case)
        tab = AeroModel(case.mesh, case.modes, case.mach, case.bref).table([0.0])
        write_cp(out, tab, "Steady dCp per mode")
        sa = StaticAeroelastic.from_table(case.M, case.K, tab)
        qd = sa.divergence()
        rho = density(v)
        out.summary["q_divergence"] = _num(qd) if qd else None
        out.summary["V_divergence"] = _num(math.sqrt(2 * qd / rho)) if qd else None
        out.log(f"divergence dynamic pressure: {qd} Pa" + (f", V_D = {math.sqrt(2 * qd / rho):.2f} m/s at rho {rho}"
                                                          if qd else " (none)"))
        if v.get("plots", True):
            plot_modes(out, case); plot_cp(out, case, tab, 0)
        return
    case, tab, st = _demo_parts(v, base, out)
    ac = case.modal
    rows = load_rows_from_table(tab, case.mesh, st)
    inert = {}
    for s in st:
        inert.update({lab: r for lab, r in zip(s.labels, modal_load_matrices(ac, s)[0])})
    sa = StaticAeroelastic.from_table(ac.M, ac.K, tab, load_rows=rows, inertia_rows=inert)
    V, rho, nz = float(v["velocity"]), density(v), float(v.get("nz", 1.0))
    free = v.get("trim_free", [v.get("control", "elevator")])
    free = [free] if isinstance(free, str) else list(free)
    tr = sa.solve(0.5 * rho * V * V, fixed={"plunge": 0.0}, free=free, accel={"plunge": nz * G0})
    out.log(tr.summary())
    out.summary.update({"velocity": V, "density": rho, "nz": nz,
                        "alpha_deg": _num(np.degrees(tr["pitch"])),
                        **{f"{f}_deg": _num(np.degrees(tr[f])) for f in free},
                        "loads": {k: _num(x) for k, x in tr.loads.items()}})
    (out.dir / "trim.txt").write_text(tr.summary() + "\n", encoding="utf-8")
    out.file("trim.txt", "text", "Trim result")
    if tr.cp is not None:
        cp = np.real(np.asarray(tr.cp)).ravel()
        out.data("cp", {"k": [0.0], "modes": ["trim"], "re": [[jlist(cp, 5)]], "im": [[[0.0] * cp.size]],
                             "mach": case.mach, "bref": case.bref}, "cp", "Trimmed steady dCp")
        if v.get("plots", True):
            plt = _plt()
            from matplotlib.collections import PolyCollection
            fig, ax = plt.subplots(figsize=(7, 4))
            pc = PolyCollection(case.mesh.corners[:, :, :2], array=cp, cmap="viridis", edgecolor="none")
            ax.add_collection(pc); ax.autoscale(); ax.set_aspect("equal"); fig.colorbar(pc, ax=ax, label="dCp")
            ax.set_title(f"Trimmed dCp, n_z = {nz}, V = {V} m/s")
            out.savefig(fig, "trim_cp.png", "Trimmed pressure")


def route_maneuver(out, v, base):
    from ..loads import station_output_group
    from ..maneuver import abrupt_input, checked_pitch_maneuver, short_period, simulate_control
    case, tab, st = _demo_parts(v, base, out)
    sysm = _system(case, tab, n_ctrl=len(case.controls))
    station_output_group(sysm, tab, case.mesh, case.modal, st, "both")
    V, rho = float(v["velocity"]), density(v)
    ctrl = v.get("control", "elevator")
    wn, z, _ = short_period(sysm, V, rho)
    out.summary.update({"velocity": V, "density": rho, "short_period_wn": _num(wn), "short_period_zeta": _num(z)})
    if v["mode"] == "checked":
        r = checked_pitch_maneuver(sysm, V, rho, ctrl, nz_limit=float(v.get("nz_limit", 2.5)))
    else:
        t = np.arange(0, float(v.get("t_end", 4.0)), float(v.get("dt", 0.005)))
        # nose-up command = trailing edge up = negative deflection
        d, dd, ddd = abrupt_input(t, -np.radians(float(v.get("amplitude_deg", 1.0))))
        r = simulate_control(sysm, V, rho, ctrl, d, dd, ddd, t)
    series = {f"{ctrl} [deg]": np.degrees(r.delta), "n_z": r.nz}
    for l in r.labels:
        if l.endswith(".fsm") or l.startswith("eta:pitch"):
            series[l] = r.series(l)
    write_timehist(out, "maneuver", r.t, series, f"{v['mode']} pitch manoeuvre, V = {V} m/s",
                   {"scale": r.scale, "control": ctrl})
    out.summary["peak_control_deg"] = _num(np.degrees(np.abs(r.delta).max()))
    out.summary["peak_nz"] = _num(r.nz.max())
    out.summary["peak_loads"] = {l: _num(np.abs(r.series(l)).max()) for l in r.labels if l.endswith(".fsm")}
    out.log(f"short period wn = {wn:.3f} rad/s, zeta = {z:.3f}; peak n_z = {r.nz.max():.3f}")


def route_tsd(out, v, base):
    from ..aero.tsd import TSDSolver
    if v["mode"] == "flutter":
        return route_flutter(out, dict(v, mode="pk", aero_method="tsd"), base)
    if v["mode"] == "gaf":
        return route_aero(out, dict(v, mode="gaf", aero_method="tsd"), base)
    case = build_case(v, base, need_aero=True)
    write_aero_mesh(out, case)
    th = tsd_thickness(v.get("tsd_thickness", "naca65a004"))
    grid = {k2: int(v[k1]) for k1, k2 in (("tsd_nx", "nchord"), ("tsd_ny", "nspan")) if k1 in v}
    if tsd_lifting(v) or "alphas" in v:
        return _tsd_steady_lifting(out, v, case, th, grid)
    t0 = time.perf_counter()
    sol = TSDSolver(case.mesh.surfaces[0], case.mach, th, verbose=v.get("verbose", True), **grid)
    st = sol.steady()
    out.log(f"TSD steady: {st.iterations} iterations, converged {st.converged}, "
            f"{st.n_supersonic} supersonic points, Cp* = {st.cp_star:.4f} ({time.perf_counter() - t0:.1f} s)")
    xc = (st.x - st.x[:, :1]) / np.maximum(np.ptp(st.x, axis=1, keepdims=True), 1e-12)
    out.data("tsd_steady", {"x": jlist(st.x), "y": jlist(st.y), "cp": jlist(st.cp), "cp_star": _num(st.cp_star),
                                  "mach": case.mach, "converged": bool(st.converged),
                                  "residuals": jlist(st.residuals)}, "tsd_steady", "TSD steady Cp")
    out.summary.update({"converged": bool(st.converged), "iterations": int(st.iterations),
                        "n_supersonic": int(st.n_supersonic), "cp_min": _num(np.nanmin(st.cp)), "cp_star": _num(st.cp_star)})
    if v.get("plots", True):
        plt = _plt()
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for j in range(0, st.cp.shape[0], max(1, st.cp.shape[0] // 5)):
            ax[0].plot(xc[j], st.cp[j], label=f"y = {st.y[j]:.3f}")
        ax[0].axhline(st.cp_star, color="k", ls="--", lw=0.8, label="Cp*")
        ax[0].invert_yaxis(); ax[0].set_xlabel("x/c"); ax[0].set_ylabel("Cp (upper)"); ax[0].legend(fontsize=7)
        ax[0].grid(alpha=0.3)
        Y = np.repeat(st.y[:, None], st.x.shape[1], 1)
        cs = ax[1].contourf(st.x, Y, st.cp, 30, cmap="RdBu")
        fig.colorbar(cs, ax=ax[1], label="Cp"); ax[1].set_aspect("equal"); ax[1].set_title("Upper-surface Cp")
        fig.suptitle(f"TSD steady flow, M = {case.mach}")
        out.savefig(fig, "tsd_steady.png", "TSD steady Cp")
        if len(st.residuals):
            write_curves(out, "tsd_residual", np.arange(1, len(st.residuals) + 1), {"residual": st.residuals},
                         "Newton residual", "iteration", "residual", logy=True)


def _tsd_steady_data(st, mach):
    return {"x": jlist(st.x), "y": jlist(st.y), "cp": jlist(st.cp), "cp_lower": jlist(st.cp_lower),
            "cp_star": _num(st.cp_star), "mach": mach, "converged": bool(st.converged),
            "residuals": jlist(st.residuals), "alpha_deg": _num(st.alpha_deg), "CL": _num(st.CL),
            "cl": jlist(st.cl), "shock_upper": jlist(st.shock_upper), "jump_upper": jlist(st.jump_upper),
            "shock_lower": jlist(st.shock_lower), "jump_lower": jlist(st.jump_lower)}


def _tsd_steady_lifting(out, v, case, th, grid):
    """Lifting TSD mean flow (incidence / camber), optionally an incidence sweep ``alphas``."""
    from ..aero.tsd_lifting import LiftingTSDSolver
    alphas = list(np.linspace(*map(float, v["alphas"][:2]), int(v["alphas"][2]))) if "alphas" in v \
        else [tsd_alpha(v)]
    camber = tsd_camber(v.get("tsd_camber"))
    sol = LiftingTSDSolver(case.mesh.surfaces[0], case.mach, th, alpha_deg=alphas[0], camber=camber,
                           verbose=v.get("verbose", True), **grid)
    res = []
    for i, a in enumerate(alphas):
        sol.set_alpha(a)
        t0 = time.perf_counter()
        st = sol.steady(warm_start=i > 0)            # continuation from the previous incidence
        res.append(st)
        out.log(f"TSD lifting steady alpha = {a:g} deg: CL = {st.CL:.4f}, {st.iterations} iterations, "
                f"converged {st.converged}, supersonic points upper/lower {st.n_supersonic}/{st.n_supersonic_lower}, "
                f"root shock x/c upper {st.shock_upper[0]:.3f} lower {st.shock_lower[0]:.3f} "
                f"({time.perf_counter() - t0:.1f} s)")
        if len(alphas) > 1:
            out.data(f"tsd_steady_a{i + 1:02d}", _tsd_steady_data(st, case.mach), "tsd_steady",
                     f"TSD steady Cp, alpha = {a:g} deg")
    st = res[-1]
    out.data("tsd_steady", _tsd_steady_data(st, case.mach), "tsd_steady", f"TSD steady Cp, alpha = {st.alpha_deg:g} deg")
    out.log(st.table())
    A = np.array(alphas)
    CL = np.array([r.CL for r in res])
    jm = st.cl.size // 2
    sw = {"CL": CL,
          "xs_upper_root": [r.shock_upper[0] for r in res], "xs_lower_root": [r.shock_lower[0] for r in res],
          f"xs_upper_y{jm + 1}": [r.shock_upper[jm] for r in res],
          f"xs_lower_y{jm + 1}": [r.shock_lower[jm] for r in res]}
    out.summary.update({"converged": all(r.converged for r in res), "alpha_deg": jlist(A), "CL": jlist(CL),
                        "n_supersonic": int(st.n_supersonic), "n_supersonic_lower": int(st.n_supersonic_lower),
                        "cp_min": _num(np.nanmin(st.cp)), "cp_star": _num(st.cp_star),
                        "shock_upper_root": jlist(sw["xs_upper_root"]),
                        "shock_lower_root": jlist(sw["xs_lower_root"]),
                        "camber": v.get("tsd_camber", "flat"), "thickness": v.get("tsd_thickness", "naca65a004")})
    if len(alphas) > 1:
        out.data("curves_tsd_alpha_sweep", {"x": jlist(A, 6), "curves": {k: jlist(c, 6) for k, c in sw.items()},
                                            "xlabel": "alpha [deg]", "ylabel": "CL, shock x/c", "logy": False,
                                            "meta": {"mach": case.mach, "y_mid": _num(st.y[jm])}},
                 "curves", "TSD incidence sweep: CL and shock positions")
    if not v.get("plots", True):
        return
    plt = _plt()
    xi = st.xi
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    for j in range(0, st.cp.shape[0], max(1, st.cp.shape[0] // 4)):
        l, = ax[0].plot(xi[j], st.cp[j], label=f"y = {st.y[j]:.3f}")
        ax[0].plot(xi[j], st.cp_lower[j], ls="--", color=l.get_color())
    ax[0].axhline(st.cp_star, color="k", ls=":", lw=0.8, label="Cp*")
    ax[0].invert_yaxis(); ax[0].set_xlabel("x/c"); ax[0].set_ylabel("Cp (solid upper, dashed lower)")
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    Y = np.repeat(st.y[:, None], st.x.shape[1], 1)
    for a_, cp, ttl in ((ax[1], st.cp, "Upper-surface Cp"), (ax[2], st.cp_lower, "Lower-surface Cp")):
        cs = a_.contourf(st.x, Y, cp, 30, cmap="RdBu")
        fig.colorbar(cs, ax=a_, label="Cp"); a_.set_aspect("equal"); a_.set_title(ttl)
    fig.suptitle(f"TSD lifting steady flow, M = {case.mach}, alpha = {st.alpha_deg:g} deg, CL = {st.CL:.4f}")
    out.savefig(fig, "tsd_steady.png", "TSD steady Cp (upper / lower)")
    if len(st.residuals):
        write_curves(out, "tsd_residual", np.arange(1, len(st.residuals) + 1), {"residual": st.residuals},
                     "Newton residual", "iteration", "residual", logy=True)
    if len(alphas) > 1:
        fig, ax = plt.subplots(1, 3, figsize=(15, 4))
        cmap = plt.get_cmap("viridis")
        for i, r in enumerate(res):
            c = cmap(i / max(1, len(res) - 1))
            ax[0].plot(xi[0], r.cp[0], color=c, label=f"{r.alpha_deg:g} deg")
            ax[0].plot(xi[0], r.cp_lower[0], color=c, ls="--")
        ax[0].axhline(st.cp_star, color="k", ls=":", lw=0.8)
        ax[0].invert_yaxis(); ax[0].set_xlabel("x/c"); ax[0].set_ylabel("root Cp (solid upper, dashed lower)")
        ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=0.3)
        ax[1].plot(A, CL, "o-"); ax[1].set_xlabel("alpha [deg]"); ax[1].set_ylabel("CL"); ax[1].grid(alpha=0.3)
        for k, ls in (("xs_upper_root", "o-"), ("xs_lower_root", "s--"), (f"xs_upper_y{jm + 1}", "o-"),
                      (f"xs_lower_y{jm + 1}", "s--")):
            ax[2].plot(A, sw[k], ls, label=k.replace("xs_", "").replace(f"y{jm + 1}", f"y = {st.y[jm]:.2f}"))
        ax[2].set_xlabel("alpha [deg]"); ax[2].set_ylabel("shock x/c"); ax[2].set_ylim(0, 1)
        ax[2].legend(fontsize=7); ax[2].grid(alpha=0.3)
        fig.suptitle(f"TSD incidence sweep, M = {case.mach}")
        out.savefig(fig, "tsd_alpha_sweep.png", "TSD incidence sweep: Cp, CL, shock positions")


def route_uvlm(out, v, base):
    from . import vlm_routes
    return vlm_routes.route_uvlm(out, v, base)


def route_freeplay(out, v, base):
    from . import ext_routes
    return ext_routes.route_freeplay(out, v, base)


def route_ase(out, v, base):
    from . import ext_routes
    return ext_routes.route_ase(out, v, base)


ROUTE_FUNCS = {"modal": route_modal, "aero": route_aero, "flutter": route_flutter, "gust": route_gust,
               "trim": route_trim, "maneuver": route_maneuver, "tsd": route_tsd, "uvlm": route_uvlm,
               "freeplay": route_freeplay, "ase": route_ase}


def run(values: dict, base: Path, *, verbose=True, config_path: Path | None = None) -> dict:
    """Run a resolved config. Returns the results dictionary (also written to results.txt)."""
    outdir = Path(values["output_dir"])
    if not outdir.is_absolute():
        outdir = (base / outdir).resolve()
    out = Output(outdir, verbose, also_json=bool(values.get("output_json", False)))
    write_config(values, outdir / "config_resolved.txt",
                 header=f"FLAAP resolved config ({config_path or 'in-memory'}) - defaults included")
    out.file("config_resolved.txt", "config", "Resolved config")
    t0 = time.perf_counter()
    status, err = "ok", None
    out.log(f"== FLAAP {values['model']}/{values['mode']}: {values.get('title', '')}")
    try:
        ROUTE_FUNCS[values["model"]](out, values, base)
    except Exception as exc:                       # noqa: BLE001 - recorded in results.txt
        status, err = "error", f"{type(exc).__name__}: {exc}"
        out.log(traceback.format_exc())
    res = {"title": values.get("title", values["model"]), "model": values["model"], "mode": values["mode"],
           "status": status, "error": err, "elapsed_s": round(time.perf_counter() - t0, 2),
           "created": time.strftime("%Y-%m-%d %H:%M:%S"), "config": str(config_path) if config_path else None,
           "output_dir": str(outdir),
           "summary": out.summary, "artifacts": out.artifacts}
    write_result(outdir / "results.txt", "results", res, res["title"])
    if out.also_json:
        (outdir / "results.json").write_text(json.dumps(res, indent=2, allow_nan=False, default=_num),
                                             encoding="utf-8")
    out.log(f"== {status} in {res['elapsed_s']} s -> {outdir}")
    out.close()
    return res
