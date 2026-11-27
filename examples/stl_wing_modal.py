"""STL wing -> tet10 modes -> IPS spline -> DLM PK flutter.

Steps
-----
1. Generate a closed, tapered, swept NACA 0012 wing surface and write it as
   ``examples/geometry/wing.stl`` (:func:`naca_wing_surface`, binary STL).
2. Mesh the STL volume with gmsh (tet4 -> tet10, target ``--size``), clamp
   the root (``y <= 0``) and compute the normal modes
   (:func:`flaap.structures.geometry_modal.geometry_modal`).
3. Spline the modes from a planform-uniform node subset to a doublet-lattice
   mesh of the same planform (infinite-plate spline).
4. Generalized aerodynamic forces (DLM) and PK flutter sweep at sea level.
5. Figures: ``validation/figures/stl_wing_modes.png`` (mode shapes) and
   ``validation/figures/stl_wing_vg.png`` (V-g / V-f).

The default material is a light balsa/foam-like isotropic solid
(E = 0.4 GPa, rho = 120 kg/m^3) so that this solid scale-model wing flutters
(bending-torsion, about 160 m/s at sea level, M = 0.3) inside the sweep; use
``--material aluminium_7075`` for a (much stiffer) solid metal wing.

Run:  python examples/stl_wing_modal.py [--size 0.03] [--nmodes 6] [--no-plot]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, ROOT)

from flaap.aero.gaf import AeroModel, AeroModes, default_k_list  # noqa: E402
from flaap.flutter import identify, pk_method  # noqa: E402
from flaap.geometry import Surface, build_mesh  # noqa: E402
from flaap.structures.geometry_modal import (geometry_modal, naca_wing_surface,  # noqa: E402
                                             planform_nodes, plot_modes, write_stl)
from flaap.structures.modal import ModalModel  # noqa: E402
from flaap.structures.spline import surface_spline  # noqa: E402

# planform [m]
SEMISPAN, ROOT_CHORD, TIP_CHORD, SWEEP_LE, TC = 1.2, 0.40, 0.25, 10.0, 0.12
LIGHT = dict(E=0.4e9, nu=0.3, rho=120.0, name="balsa-like isotropic")
STL_PATH = os.path.join(HERE, "geometry", "wing.stl")
FIG_DIR = os.path.join(ROOT, "validation", "figures")


def make_stl(path: str = STL_PATH, n_chord: int = 20, n_span: int = 10) -> str:
    V, F = naca_wing_surface(SEMISPAN, ROOT_CHORD, TIP_CHORD, SWEEP_LE, TC, n_chord, n_span)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    write_stl(path, V, F)
    return path


def aero_surface(nchord: int = 8, nspan: int = 12) -> Surface:
    tip_x = SEMISPAN * np.tan(np.deg2rad(SWEEP_LE))
    return Surface("wing", (0.0, 0.0, 0.0), ROOT_CHORD, (tip_x, SEMISPAN, 0.0), TIP_CHORD,
                   nchord, nspan)


def run(size: float = 0.03, nmodes: int = 6, material=None, mach: float = 0.3,
        rho: float = 1.225, velocities=None, plot: bool = True, verbose: bool = True,
        stl_path: str = STL_PATH, element: str = "tet10") -> dict:
    t0 = time.perf_counter()
    make_stl(stl_path)
    res = geometry_modal(stl_path, nmodes, kind="solid", element=element, size=size,
                         material=material or LIGHT, clamp=dict(below=("y", 0.0)), verbose=verbose)
    modal = res.modal
    t1 = time.perf_counter()

    # spline from a planform-uniform node subset (avoids a huge dense IPS system)
    mask = planform_nodes(modal.nodes, 24, 8, axes=(1, 0), exclude=res.clamped)
    mask[res.clamped[::max(1, res.clamped.size // 8)]] = True          # a few root nodes (u = 0)
    sub = ModalModel(modal.nodes[mask], modal.shapes[mask], modal.omega, modal.gen_mass,
                     name=modal.name, labels=modal.labels)
    mesh = build_mesh([aero_surface()], "symmetric")
    sp = surface_spline(mesh, sub.nodes, method="ips")
    modes = AeroModes.from_spline(sp, sub.shapes, [f"m{i + 1}" for i in range(nmodes)])
    bref = 0.5 * 0.5 * (ROOT_CHORD + TIP_CHORD)
    kv = default_k_list(1.5, 14)
    table = AeroModel(mesh, modes, mach, bref).table(kv, gust=False, keep_cp=False)
    V = np.linspace(10.0, 300.0, 146) if velocities is None else np.asarray(velocities, float)
    fl = pk_method(modal.M, modal.K, table, rho, V)
    identify(fl)
    t2 = time.perf_counter()

    out = dict(freq_hz=modal.freq_hz, info=res.info, n_spline_nodes=int(mask.sum()),
               n_boxes=mesh.n, flutter=fl, critical=fl.critical(),
               t_modal=t1 - t0, t_flutter=t2 - t1, result=res)
    if verbose:
        print(f"spline nodes: {mask.sum()}, DLM boxes: {mesh.n}, bref = {bref:.4f} m, M = {mach}")
        print(fl.summary())
        print(f"time: modal {t1 - t0:.2f} s, aero + flutter {t2 - t1:.2f} s")
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        os.makedirs(FIG_DIR, exist_ok=True)
        f1 = os.path.join(FIG_DIR, "stl_wing_modes.png")
        plot_modes(modal, range(min(nmodes, 6)), path=f1,
                   suptitle=f"STL wing, {res.info['nodes']} nodes {element}", view=(30, -125))
        plt.close("all")
        ax = fl.plot(title=f"STL wing PK flutter (M={mach}, rho={rho})")
        f2 = os.path.join(FIG_DIR, "stl_wing_vg.png")
        ax[0].figure.savefig(f2, dpi=120)
        plt.close("all")
        out["figures"] = [f1, f2]
        if verbose:
            print("figures:", f1, f2)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--size", type=float, default=0.03)
    p.add_argument("--nmodes", type=int, default=6)
    p.add_argument("--material", default=None)
    p.add_argument("--mach", type=float, default=0.3)
    p.add_argument("--element", default="tet10", choices=["tet4", "tet10"])
    p.add_argument("--no-plot", action="store_true")
    a = p.parse_args()
    run(a.size, a.nmodes, a.material, a.mach, plot=not a.no_plot, element=a.element)
