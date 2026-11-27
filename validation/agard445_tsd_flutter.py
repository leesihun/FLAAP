"""Validation: AGARD 445.6 transonic flutter - DLM vs TSD aerodynamics.

Same structural model, mass-ratio scaling and PK solution as
``agard445_flutter.py`` (measured natural frequencies, FE mode shapes), but
the generalized aerodynamic forces come from the transonic small-disturbance
solver :mod:`flaap.aero.tsd` with the NACA 65A004 thickness mean flow
("TSD") and, for reference, a zero-thickness TSD run ("TSD flat", which
should reproduce the DLM).

Run:  python validation/agard445_tsd_flutter.py [--machs 0.901 0.954 0.960] [--plot out.png]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from flaap.aero.gaf import AeroModel, AeroModes  # noqa: E402
from flaap.aero.tsd import modes_from_structure, naca65a004, tsd_gaf_table  # noqa: E402
from flaap.flutter import identify, pk_method  # noqa: E402
from flaap.geometry import build_mesh  # noqa: E402
from flaap.models.agard445 import (B_REF, EXPERIMENT_FREQ_HZ, FLUTTER_DATA, OMEGA_ALPHA,  # noqa: E402
                                   agard445_fe, agard445_surface, frustum_volume)
from flaap.structures.spline import surface_spline  # noqa: E402

K_LIST = np.array([1e-4, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.9])


def _flutter(modal, table, rho):
    V = np.linspace(50.0, 700.0, 261)
    res = pk_method(modal.M, modal.K, table, rho, V)
    identify(res)
    fl = [p for p in res.points if p.kind == "flutter"]
    return min(fl, key=lambda p: p.velocity) if fl else None


def run(machs=(0.901, 0.954, 0.960), flat=True, verbose=True):
    fe = agard445_fe(n_chord=16, n_span=16)
    modal = fe.modes(4, name="AGARD 445.6")
    modal.omega = 2 * np.pi * EXPERIMENT_FREQ_HZ.copy()
    u = np.zeros(fe.ndof)
    u[2::6] = 1.0
    m_ratio = float(u @ (fe.M @ u)) / (14.315 * frustum_volume())
    surf = agard445_surface(8, 10)
    mesh = build_mesh([surf], "symmetric")
    modes = AeroModes.from_spline(surface_spline(mesh, modal.nodes, method="ips"), modal.dof_vectors())
    funcs, names = modes_from_structure(modal)
    rows = []
    if verbose:
        print("   M    FSI_exp   FSI_DLM  FSI_TSDflat  FSI_TSD    err_DLM%  err_TSD%   [time]")
    for M in machs:
        d = FLUTTER_DATA[M]
        rho = d["rho"] * m_ratio
        t0 = time.perf_counter()
        out = {}
        tabs = {"DLM": AeroModel(mesh, modes, M, B_REF).table(K_LIST, gust=False, keep_cp=False)}
        if flat:
            tabs["TSD flat"] = tsd_gaf_table(surf, M, K_LIST, B_REF, funcs, names, keep_cp=False)
        tabs["TSD"] = tsd_gaf_table(surf, M, K_LIST, B_REF, funcs, names, thickness=naca65a004(),
                                    keep_cp=False)
        for key, tab in tabs.items():
            c = _flutter(modal, tab, rho)
            out[key] = (np.nan, np.nan) if c is None else (
                c.velocity / (B_REF * OMEGA_ALPHA * np.sqrt(d["mu"])), 2 * np.pi * c.freq_hz / OMEGA_ALPHA)
        rows.append((M, d["FSI"], out["DLM"][0], out.get("TSD flat", (np.nan,))[0], out["TSD"][0],
                     d["freq_ratio"], out["DLM"][1], out["TSD"][1]))
        if verbose:
            r = rows[-1]
            print(f"{M:6.3f}  {r[1]:7.4f}  {r[2]:8.4f}  {r[3]:10.4f}  {r[4]:8.4f}   "
                  f"{100 * (r[2] / r[1] - 1):7.1f}  {100 * (r[4] / r[1] - 1):7.1f}   "
                  f"[{time.perf_counter() - t0:.0f} s]")
    return np.array(rows, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--machs", type=float, nargs="*", default=[0.901, 0.954, 0.960])
    ap.add_argument("--no-flat", action="store_true")
    ap.add_argument("--plot", default=None)
    a = ap.parse_args()
    r = run(tuple(a.machs), flat=not a.no_flat)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(r[:, 0], r[:, 1], "ko-", label="experiment (Yates 1987)")
        ax.plot(r[:, 0], r[:, 2], "bs--", label="DLM + PK")
        ax.plot(r[:, 0], r[:, 4], "r^-", label="TSD (65A004 mean flow) + PK")
        ax.set_xlabel("Mach"); ax.set_ylabel("flutter speed index"); ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout(); fig.savefig(a.plot, dpi=120)
        print("saved", a.plot)


if __name__ == "__main__":
    main()
