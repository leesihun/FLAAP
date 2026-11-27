"""Command line interface: ``python -m flaap <command> ...``

Commands
--------
``run CONFIG.txt [-o DIR]``
    preflight and run a flat ``key value`` config file - the same as
    ``python FLAAP_main.py --config CONFIG.txt`` (see ``docs/CONFIGURATION.md``).
``legacy [DIR] [--method pk|k] [--vmin --vmax --nv] [-o DIR]``
    run a Fortran-FLAAP input directory (``LEE_input.txt``, ``M_K.txt``,
    ``MODESHAPE_input.txt``, ``atmosphere_input.txt``) and write
    ``Vg_save.dat``, ``Vf_save.dat``, ``pk_results_mode*.dat`` to
    ``output/legacy_<DIR name>`` (or ``-o DIR``) - the replacement of ``FLAAP.exe`` (the original prompted for the method and
    ``Vmin / Vmax / Vinc`` interactively).
``mesh FILE [--size H --order 1|2 --recombine --dim 2|3 -o OUT.msh]``
    mesh a gmsh ``.geo`` / STEP / IGES / BREP (or read a ``.msh``) and print
    the FE element blocks and physical groups (``structure gmsh`` uses the same mesher).
``f06 FILE``
    print the modes read from an MSC Nastran SOL 103 ``.f06``.
``validate [goland|agard|tsd|all]``
    quick validation runs.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np


def _legacy(a):
    from .io.legacy import run_legacy
    V = np.arange(a.vmin, a.vmax + 0.5 * a.vinc, a.vinc)
    # default output folder: output/legacy_<deck folder>, never the deck folder itself, so the
    # original Fortran results kept next to a deck (legacy_fortran/Vg_save.dat, ...) stay untouched
    name = os.path.basename(os.path.abspath(a.directory)) or "deck"
    out = a.output or os.path.join("output", f"legacy_{name}")
    res, case, _ = run_legacy(a.directory, velocities=V, method=a.method, output=out,
                              spline=a.spline, verbose=a.verbose)
    print(f"results -> {out}")
    print(case.modal().table())
    print(res.summary())


def _f06(a):
    from .io.nastran import read_f06_modal
    m = read_f06_modal(a.file, nmodes=a.nmodes)
    print(f"{m.nodes.shape[0]} grids")
    print(m.table())


def _validate(a):
    which = a.case
    if which in ("goland", "all"):
        from .aero.gaf import AeroModel, AeroModes
        from .flutter import pk_method
        from .geometry import build_mesh
        from .models.goland import CHORD, REFERENCE, goland_modal, goland_surface
        from .structures.spline import beam_spline
        mm = goland_modal(4)
        mesh = build_mesh([goland_surface(6, 12)], "symmetric")
        modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors())
        tab = AeroModel(mesh, modes, 0.0, CHORD / 2).table(np.r_[0, np.geomspace(0.01, 2, 19)], keep_cp=False)
        r = pk_method(mm.M, mm.K, tab, 1.225, np.linspace(60, 260, 101))
        c = r.critical()
        print("Goland wing, DLM 6x12 + PK:")
        print(f"  natural frequencies {np.round(mm.omega, 1)} rad/s (exact {REFERENCE['omega_rad_s']})")
        print(f"  flutter V = {c.velocity:.1f} m/s, w = {c.omega:.1f} rad/s "
              f"(strip theory 137.2 / 70.7 at rho 1.225; 3-D lifting-surface results "
              f"{REFERENCE['dlm_flutter_speed_range']} m/s at rho {REFERENCE['dlm_density']})")
    if which in ("agard", "all"):
        import os
        import runpy
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "validation", "agard445_flutter.py")
        sys.argv = [path, "--exp-freq"]
        runpy.run_path(path, run_name="__main__")
    if which == "tsd":
        import os
        import runpy
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "validation", "agard445_tsd_flutter.py")
        sys.argv = [path]
        runpy.run_path(path, run_name="__main__")


def _mesh(a):
    from .structures.gmsh_mesher import gmsh_mesh
    gm = gmsh_mesh(a.file, size=a.size, order=a.order, recombine=a.recombine or None, dim=a.dim,
                   length_scale=a.length_scale, save=a.output)
    lo, hi = gm.nodes.min(axis=0), gm.nodes.max(axis=0)
    print(f"{a.file}: {'meshed' if gm.info['generated'] else 'read'} ({gm.dim}-D)")
    print(f"  {len(gm.nodes)} nodes, {gm.nelem} elements")
    for kind, conn, grp in gm.blocks:
        print(f"    {kind:6s} {len(conn):8d}  group {grp}")
    print(f"  bounding box [m]: {np.round(lo, 6)} .. {np.round(hi, 6)}")
    if gm.groups:
        print("  physical groups: " + ", ".join(f"{k} ({len(v)} nodes)" for k, v in sorted(gm.groups.items())))
    if a.output:
        print(f"  written: {a.output}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="flaap", description="FLAAP aeroelastic analysis (Python)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="run a key-value .txt config file")
    p.add_argument("case")
    p.add_argument("-o", "--output", default=None)
    p = sub.add_parser("legacy", help="run a Fortran-FLAAP input directory")
    p.add_argument("directory", nargs="?", default=".")
    p.add_argument("--method", choices=("pk", "k"), default="pk")
    p.add_argument("--vmin", type=float, default=5.0)
    p.add_argument("--vmax", type=float, default=300.0)
    p.add_argument("--vinc", type=float, default=1.0)
    p.add_argument("--spline", choices=("legacy", "ips"), default="legacy")
    p.add_argument("-o", "--output", default=None)
    p.add_argument("-v", "--verbose", action="store_true")
    p = sub.add_parser("f06", help="read modes from a Nastran .f06")
    p.add_argument("file")
    p.add_argument("--nmodes", type=int, default=None)
    p = sub.add_parser("validate", help="validation runs")
    p.add_argument("case", nargs="?", default="all", choices=("goland", "agard", "tsd", "all"))
    p = sub.add_parser("mesh", help="mesh a .geo / STEP (or read a .msh) with gmsh and report the FE mesh")
    p.add_argument("file")
    p.add_argument("--size", type=float, default=None, help="Mesh.MeshSizeMax [m]")
    p.add_argument("--order", type=int, choices=(1, 2), default=None)
    p.add_argument("--recombine", action="store_true", help="quad4 surface elements")
    p.add_argument("--dim", type=int, choices=(2, 3), default=None)
    p.add_argument("--length-scale", type=float, default=1.0)
    p.add_argument("-o", "--output", default=None, help="write the mesh (.msh, .vtk, .inp, ... by extension)")
    a = ap.parse_args(argv)
    if a.cmd == "mesh":
        _mesh(a)
        return
    if a.cmd == "run":
        from .launcher.main import main as launcher_main
        args = ["--config", a.case] + (["--output-dir", a.output] if a.output else [])
        sys.exit(launcher_main(args))
    elif a.cmd == "legacy":
        _legacy(a)
    elif a.cmd == "f06":
        _f06(a)
    elif a.cmd == "validate":
        _validate(a)


if __name__ == "__main__":
    main()
