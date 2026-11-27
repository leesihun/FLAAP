"""Validation: time-linearised TSD (flaap.aero.tsd) vs. the doublet-lattice
method (flaap.aero.gaf) for rigid plunge / pitch generalized forces.

In the linear (zero-thickness, subsonic) limit the TSD solution must
reproduce linear potential theory, so the DLM is the reference.  The last
case adds NACA 65A004 thickness at M = 0.96 (AGARD 445.6 planform), where the
two methods are *expected* to differ (transonic mean-flow effects).

Run:  python validation/tsd_vs_dlm.py [--fine]
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from flaap.geometry import Surface, build_mesh  # noqa: E402
from flaap.aero.gaf import AeroModel, rigid_aero_modes  # noqa: E402
from flaap.aero.tsd import naca65a004, rigid_mode_funcs, tsd_gaf_table  # noqa: E402
from flaap.models.agard445 import agard445_surface  # noqa: E402


def agard_planform(nchord=10, nspan=10):
    return agard445_surface(nchord=nchord, nspan=nspan)


def rectangular(c=1.0, s=2.0, nchord=10, nspan=16):
    return Surface("rect", (0, 0, 0), c, (0, s, 0), c, nchord=nchord, nspan=nspan)


def compare(title, surf_dlm, surf_tsd, mach, kvals, thickness=None, **grid_kw):
    c_r = surf_tsd.root_chord
    bref, xref = 0.5 * c_r, 0.25 * c_r
    mesh = build_mesh([surf_dlm], "symmetric")
    t0 = time.perf_counter()
    dlm = AeroModel(mesh, rigid_aero_modes(mesh, xref), mach, bref).table(kvals, gust=False)
    t_dlm = time.perf_counter() - t0
    funcs, names = rigid_mode_funcs(xref)
    t0 = time.perf_counter()
    tsd = tsd_gaf_table(surf_tsd, mach, kvals, bref, funcs, names, thickness=thickness, **grid_kw)
    t_tsd = time.perf_counter() - t0
    sol = tsd.tsd_solver
    print(f"\n=== {title}: M = {mach}, bref = {bref:.4f}, x_ref = {xref:.4f} ===")
    print(f"    {sol.grid_summary()}")
    if thickness is not None:
        mf = sol.mean_flow
        print(f"    mean flow: {mf.iterations} Newton it, converged={mf.converged}, "
              f"supersonic nodes={mf.n_supersonic}, min Cp={mf.cp.min():.3f}, Cp*={mf.cp_star:.3f}")
    print(f"    DLM {t_dlm:.1f} s, TSD {t_tsd:.1f} s ({len(kvals)} k values)")
    print(f"    {'k':>5} {'term':>16} {'DLM':>22} {'TSD':>22} {'|TSD-DLM|/|DLM|':>16}")
    scale = np.abs(dlm.Q).max(axis=(1, 2))
    for i, k in enumerate(kvals):
        for a in range(2):
            for b in range(2):
                qd, qt = dlm.Q[i, a, b], tsd.Q[i, a, b]
                if abs(qd) < 1e-9 * scale[i]:
                    continue
                err = abs(qt - qd) / abs(qd)
                print(f"    {k:5.2f} {names[a] + ',' + names[b]:>16} {qd.real:10.4f}{qd.imag:+10.4f}j "
                      f"{qt.real:10.4f}{qt.imag:+10.4f}j {err:16.3f}")
    return dlm, tsd


def main(fine=False):
    kvals = [0.0, 0.1, 0.3, 0.6]
    g = dict(nchord=32, nspan=14, nz=26) if fine else {}
    compare("Rectangular AR 4 flat plate", rectangular(), rectangular(), 0.5, kvals, **g)
    compare("AGARD 445.6 planform, flat plate", agard_planform(12, 12), agard_planform(), 0.5, kvals, **g)
    compare("AGARD 445.6 planform, flat plate", agard_planform(12, 12), agard_planform(), 0.9, kvals, **g)
    compare("AGARD 445.6, NACA 65A004 thickness (TSD nonlinear mean flow)",
            agard_planform(12, 12), agard_planform(), 0.96, kvals, thickness=naca65a004(), **g)


if __name__ == "__main__":
    main(fine="--fine" in sys.argv)
