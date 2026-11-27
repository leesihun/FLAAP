"""Validation: 3-D structural FE, ANSYS / OpenRadioss readers, geometry -> modes.

Tables (``run()`` returns ``{name: list of row dicts}``; every row has
``computed``, ``reference`` and ``err_pct`` where a reference exists):

``shell_leissa``     MITC4 quads / DKT-Allman triangles, square CFFF plate,
                     a/h = 100, 8x8 and 16x16 meshes, vs Leissa lambda_1..5
``solid_beam``       hex8 (incompatible modes), tet10, tet4 cantilever
                     1 x 0.06 x 0.04 m steel vs Timoshenko (converged 1-D FE)
                     and Euler-Bernoulli
``rigid_modes``      free-free models: the 6 rigid-body frequencies (max) and
                     the first elastic frequency
``openradioss``      the four example decks vs Leissa / Timoshenko
``ansys``            MAPDL listing / PRNSOL / Workbench round trips vs the
                     source model (frequency error, min MAC)
``calculix``         CalculiX .frd (long SI + .dat, short mm-t-s with a static
                     step and rigid modes) round trips vs the source model;
                     with ``flutter``: PK flutter of the .frd config vs the
                     same plate built by FLAAP (configs/flutter/ex11)
``gmsh``             gmsh mesher: .geo plate (transfinite quad4) vs the plate FE
                     reference, .msh round trip, OCC .geo bar (tet10) vs
                     Euler-Bernoulli; with ``flutter``: PK of configs/flutter/ex12
                     vs the CalculiX .frd model (ex11)
``geometry``         STL box beam (tet10), STEP box (gmsh OCC, mm -> m),
                     flat STL remeshed into DKT shells, NACA wing STL
                     (solid tet10 vs thin-skin shell) - with runtimes
``wing_flutter``     examples/stl_wing_modal.py: STL -> tet10 -> IPS -> DLM PK

Figures (validation/figures/): ``structures_shell_convergence.png``,
``structures_beam_modes.png``, ``structures_plate_modes.png``,
``stl_wing_modes.png``, ``stl_wing_vg.png``.

Run:  python validation/structures_validation.py [--no-flutter]
"""
from __future__ import annotations

import os
import sys
import tempfile
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from flaap.io import ansys as A  # noqa: E402
from flaap.io.openradioss import radioss_modal  # noqa: E402
from flaap.structures.fem import FEMesh, dofs_of_nodes, modal_analysis, select_nodes  # noqa: E402
from flaap.structures.plate import plate_mesh_trapezoid  # noqa: E402
from flaap.structures.reference import (LEISSA_CFFF, eb_cantilever_hz, plate_lambda,  # noqa: E402
                                        timoshenko_cantilever_hz)
from flaap.structures.solid import box_mesh  # noqa: E402

FIG = os.path.join(HERE, "figures")
AL = dict(E=70e9, nu=0.3, rho=2700.0)
ST = dict(E=210e9, nu=0.3, rho=7850.0)
BEAM = (1.0, 0.06, 0.04)


def _row(case, qty, comp, ref=None, **kw):
    r = dict(case=case, quantity=qty, computed=float(comp))
    if ref is not None:
        r["reference"] = float(ref)
        r["err_pct"] = 100.0 * (float(comp) - float(ref)) / float(ref)
    r.update(kw)
    return r


def _plate(n, kind):
    nodes, quads = plate_mesh_trapezoid((0, 0), 1.0, (0, 1.0), 1.0, n, n)
    nodes = np.c_[nodes, np.zeros(len(nodes))] if nodes.shape[1] == 2 else nodes
    conn = quads if kind == "quad4" else np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    return FEMesh(nodes, name=f"plate {kind} {n}x{n}").add(kind, conn, AL, 0.01)


def _beam_refs(mat, n_z=2, n_y=1):
    L, b, h = BEAM
    E, nu, rho = mat["E"], mat["nu"], mat["rho"]
    ti = [(f, "bend-z") for f in timoshenko_cantilever_hz(L, E, nu, b * h ** 3 / 12, rho, b * h, n_z)]
    ti += [(f, "bend-y") for f in timoshenko_cantilever_hz(L, E, nu, h * b ** 3 / 12, rho, b * h, n_y)]
    eb = [f for f in eb_cantilever_hz(L, E, b * h ** 3 / 12, rho, b * h, n_z)]
    eb += [f for f in eb_cantilever_hz(L, E, h * b ** 3 / 12, rho, b * h, n_y)]
    o = np.argsort([t[0] for t in ti])
    return [ti[i] for i in o], [eb[i] for i in o]


# ----------------------------------------------------------------- tables
def shell_leissa(conv=False):
    rows, curves = [], {}
    meshes = (4, 8, 12, 16, 24) if conv else (8, 16)
    for kind in ("quad4", "tri3"):
        for n in meshes:
            mesh = _plate(n, kind)
            t0 = time.perf_counter()
            mm = modal_analysis(mesh, 5, fixed=dofs_of_nodes(select_nodes(mesh.nodes, below=("y", 0.0))))
            dt = time.perf_counter() - t0
            lam = plate_lambda(mm.omega, 1.0, 0.01, AL["E"], AL["nu"], AL["rho"])
            curves[(kind, n)] = lam
            if n in (8, 16):
                for i in range(5):
                    rows.append(_row(f"{kind} {n}x{n}", f"lambda{i + 1}", lam[i], LEISSA_CFFF[i],
                                     time_s=dt))
    return rows, curves, mm


def solid_beam():
    rows = []
    ti, eb = _beam_refs(ST)
    L, b, h = BEAM
    last = None
    for kind, (nx, ny, nz) in (("hex8", (20, 2, 2)), ("tet10", (10, 2, 2)), ("tet10", (16, 2, 2)),
                               ("tet4", (40, 4, 4))):
        nodes, conn = box_mesh(L, b, h, nx, ny, nz, kind)
        mesh = FEMesh(nodes, name=f"{kind} {nx}x{ny}x{nz}").add(kind, conn, ST)
        t0 = time.perf_counter()
        mm = modal_analysis(mesh, 3, fixed=dofs_of_nodes(select_nodes(nodes, below=("x", 0.0)), "123"))
        dt = time.perf_counter() - t0
        for i in range(3):
            rows.append(_row(mesh.name, f"f{i + 1} [Hz] ({ti[i][1]})", mm.freq_hz[i], ti[i][0],
                             euler_bernoulli=eb[i],
                             err_eb_pct=100 * (mm.freq_hz[i] - eb[i]) / eb[i], time_s=dt,
                             dofs=int(3 * mesh.nn)))
        if kind == "tet10" and nx == 16:
            last = mm
    return rows, last


def rigid_modes():
    rows = []
    for kind in ("quad4", "tri3"):
        f = modal_analysis(_plate(8, kind), 7).freq_hz
        rows.append(_row(f"free plate {kind}", "max rigid f [Hz]", f[:6].max(), first_elastic=f[6]))
    for kind in ("hex8", "tet10", "tet4"):
        nodes, conn = box_mesh(0.5, 0.06, 0.04, 8, 2, 2, kind)
        f = modal_analysis(FEMesh(nodes).add(kind, conn, ST), 7).freq_hz
        rows.append(_row(f"free beam {kind}", "max rigid f [Hz]", f[:6].max(), first_elastic=f[6]))
    return rows


def openradioss():
    rows = []
    d = os.path.join(ROOT, "examples", "openradioss")
    ti, _ = _beam_refs(ST, 2, 2)
    for deck in ("plate_cantilever_0000.rad", "plate_tria_0000.rad", "beam_brick_0000.rad",
                 "beam_tetra10_0000.rad"):
        t0 = time.perf_counter()
        mm = radioss_modal(os.path.join(d, deck), nmodes=4 if deck.startswith("beam") else 5)
        dt = time.perf_counter() - t0
        if deck.startswith("plate"):
            lam = plate_lambda(mm.omega, 1.0, 0.01, AL["E"], AL["nu"], AL["rho"])
            for i in range(5):
                rows.append(_row(deck, f"lambda{i + 1}", lam[i], LEISSA_CFFF[i], time_s=dt))
        else:
            for i in range(4):
                rows.append(_row(deck, f"f{i + 1} [Hz] ({ti[i][1]})", mm.freq_hz[i], ti[i][0],
                                 time_s=dt))
    return rows


def ansys():
    d = os.path.join(ROOT, "examples", "ansys")
    r = np.load(os.path.join(d, "reference_modes.npz"))
    f_ref, sh = r["freq_hz"], r["shapes"]

    def mac(a, b):
        a, b = a.ravel(), b.ravel()
        return (a @ b) ** 2 / ((a @ a) * (b @ b))

    rows = []
    wb = os.path.join(d, "workbench")
    cases = {
        "MAPDL SET,LIST + PRNSOL U/ROT + CDB": lambda: A.read_mapdl_modal(
            os.path.join(d, "plate_modal.out"), nodes=os.path.join(d, "plate.cdb")),
        "MAPDL PRNSOL,DOF (wrapped)": lambda: A.read_mapdl_modal(
            os.path.join(d, "plate_prnsol_dof.out"), nodes=os.path.join(d, "plate.cdb")),
        "MAPDL listing in mm-t units": lambda: A.read_mapdl_modal(
            os.path.join(d, "plate_modal_mm.out"), length_scale=1e-3, mass_scale=1e3),
        "Workbench per-mode text (UZ)": lambda: A.read_ansys_modal(wb),
    }
    for name, fn in cases.items():
        t0 = time.perf_counter()
        mm = fn()
        dt = time.perf_counter() - t0
        comps = slice(2, 3) if "Workbench" in name else slice(0, 6)
        macs = [mac(mm.shapes[:, comps, m], sh[:, comps, m]) for m in range(mm.nmodes)]
        err = np.abs(mm.freq_hz - f_ref[:mm.nmodes]) / f_ref[:mm.nmodes]
        rows.append(_row(name, "max |df/f| [%]", 100 * err.max(), min_mac=float(min(macs)),
                         nmodes=mm.nmodes, time_s=dt))
    return rows


def calculix(flutter=True):
    from flaap.io.calculix import read_calculix_modal
    d = os.path.join(ROOT, "examples", "calculix")
    r = np.load(os.path.join(d, "reference_modes.npz"))
    f_ref, sh = r["freq_hz"], r["shapes"]

    def mac(a, b):
        a, b = a.ravel(), b.ravel()
        return (a @ b) ** 2 / ((a @ a) * (b @ b))

    rows = []
    cases = {
        ".frd long format + .dat (SI)": lambda: read_calculix_modal(os.path.join(d, "plate.frd")),
        ".frd short format, mm-t-s, static step + 3 rigid modes": lambda: read_calculix_modal(
            os.path.join(d, "plate_mm_short.frd"), length_scale=1e-3, mass_scale=1e3),
    }
    for name, fn in cases.items():
        t0 = time.perf_counter()
        mm = fn()
        dt = time.perf_counter() - t0
        macs = [mac(mm.shapes[:, :3, m], sh[:, :3, m]) for m in range(mm.nmodes)]
        err = np.abs(mm.freq_hz - f_ref[:mm.nmodes]) / f_ref[:mm.nmodes]
        # generalized-mass consistency: phi^T M phi / gen_mass of the source FE model
        rows.append(_row(name, "max |df/f| [%]", 100 * err.max(), min_mac=float(min(macs)),
                         nmodes=mm.nmodes, gen_mass=float(mm.gen_mass[0]), time_s=dt))
    if flutter:
        from flaap.launcher import run_config
        cfg = os.path.join(ROOT, "configs", "flutter", "ex11", "config_calculix_pk.txt")
        tmp = tempfile.mkdtemp(prefix="flaap_ccx_")
        a = run_config(cfg, output_dir=os.path.join(tmp, "ccx"), verbose=False, overrides={"plots": False})
        plate = os.path.join(tmp, "plate.txt")
        with open(plate, "w") as fh:
            fh.write("\n".join([
                "model flutter", "mode pk", "structure plate", "plate_root_le 0, 0, 0",
                "plate_root_chord 0.5", "plate_tip_le 0, 1.0, 0", "plate_tip_chord 0.5",
                "plate_nchord 6", "plate_nspan 12", "thickness 0.004", "material custom",
                "young 71.7e9", "poisson 0.33", "material_density 2810", "nmodes 4",
                "surface1 plate, 0, 0, 0, 0.5, 0, 1.0, 0, 0.5, 6, 12", "symmetry symmetric",
                "mach 0.0", "kmax 2.0", "nk 16", "density 1.225", "velocities 10, 150, 141",
                "plots false", ""]))
        b = run_config(plate, output_dir=os.path.join(tmp, "plate"), verbose=False)
        va = float(a["summary"]["flutter_pk"]["velocity"])
        vb = float(b["summary"]["flutter_pk"]["velocity"])
        rows.append(_row("PK flutter: .frd (mm-t-s) vs FLAAP plate FE", "V_f [m/s]", va, vb,
                         f_hz=float(a["summary"]["flutter_pk"]["freq_hz"]),
                         f_ref_hz=float(b["summary"]["flutter_pk"]["freq_hz"])))
    return rows


def gmsh_mesher(flutter=True):
    """gmsh mesher: .geo plate (transfinite quads) vs the plate FE reference, solid bar vs beam theory,
    .msh round trip and the PK flutter speed of the .geo model vs the CalculiX .frd model (ex11)."""
    from flaap.structures.gmsh_mesher import gmsh_modal
    d = os.path.join(ROOT, "examples", "gmsh")
    al = dict(E=71.7e9, nu=0.33, rho=2810.0)
    f_ref = np.load(os.path.join(ROOT, "examples", "calculix", "reference_modes.npz"))["freq_hz"]
    rows = []
    t0 = time.perf_counter()
    r = gmsh_modal(os.path.join(d, "plate_wing.geo"), 6, material=al, thickness=0.004, clamp_group="Root")
    err = np.abs(r.modal.freq_hz - f_ref) / f_ref
    rows.append(_row(".geo plate, transfinite 6x12 quad4 (Recombine)", "max |df/f| [%]", 100 * err.max(),
                     nodes=r.info["nodes"], clamped=r.info["clamped_nodes"], time_s=time.perf_counter() - t0))
    tmp = tempfile.mkdtemp(prefix="flaap_gmsh_")
    msh = os.path.join(tmp, "plate.msh")
    r2 = gmsh_modal(os.path.join(d, "plate_wing.geo"), 6, material=al, thickness=0.004, clamp_group="Root",
                    save=msh)
    r3 = gmsh_modal(msh, 6, material=al, thickness=0.004, clamp_group="Root")
    rows.append(_row(".msh (MSH 4.1) round trip", "max |df/f| [%]",
                     100 * float(np.max(np.abs(r3.modal.freq_hz / r2.modal.freq_hz - 1)))))
    E, rho, b, h = al["E"], al["rho"], 0.10, 0.02
    f_eb = 1.875104 ** 2 / (2 * np.pi) * np.sqrt(E * b * h ** 3 / 12 / (rho * b * h))
    t0 = time.perf_counter()
    s = gmsh_modal(os.path.join(d, "box_beam.geo"), 2, material=al, clamp_group="Root", order=2)
    rows.append(_row(".geo OCC bar 1.0x0.10x0.02 m, tet10", "f1 [Hz] (Euler-Bernoulli)", s.modal.freq_hz[0], f_eb,
                     nodes=s.info["nodes"], time_s=time.perf_counter() - t0))
    if flutter:
        from flaap.launcher import run_config
        a = run_config(os.path.join(ROOT, "configs", "flutter", "ex12", "config_gmsh_pk.txt"),
                       output_dir=os.path.join(tmp, "gmsh"), verbose=False, overrides={"plots": False})
        c = run_config(os.path.join(ROOT, "configs", "flutter", "ex11", "config_calculix_pk.txt"),
                       output_dir=os.path.join(tmp, "ccx"), verbose=False, overrides={"plots": False})
        rows.append(_row("PK flutter: gmsh .geo plate vs CalculiX .frd", "V_f [m/s]",
                         float(a["summary"]["flutter_pk"]["velocity"]),
                         float(c["summary"]["flutter_pk"]["velocity"]),
                         f_hz=float(a["summary"]["flutter_pk"]["freq_hz"])))
    return rows


def geometry(plot=True):
    import gmsh
    import trimesh
    from flaap.structures.geometry_modal import geometry_modal, naca_wing_surface, write_stl
    rows = []
    L, b, h = BEAM
    tmp = tempfile.mkdtemp(prefix="flaap_geo_")
    ti, _ = _beam_refs(ST)
    # STL box beam
    tm = trimesh.creation.box(extents=(L, b, h))
    tm.apply_translation((L / 2, 0, 0))
    p = os.path.join(tmp, "beam.stl")
    tm.export(p)
    for size in (0.03, 0.02):
        r = geometry_modal(p, 3, kind="solid", size=size, material=ST, clamp=dict(below=("x", 0.0)))
        for i in range(3):
            rows.append(_row(f"STL box tet10 h={size}", f"f{i + 1} [Hz]", r.modal.freq_hz[i], ti[i][0],
                             nodes=r.info["nodes"], time_s=r.info["t_total"]))
    # STEP box (mm) via gmsh OCC
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("box")
    gmsh.model.occ.addBox(0, -b / 2 * 1e3, -h / 2 * 1e3, L * 1e3, b * 1e3, h * 1e3)
    gmsh.model.occ.synchronize()
    ps = os.path.join(tmp, "beam.step")
    gmsh.write(ps)
    gmsh.finalize()
    r = geometry_modal(ps, 3, kind="solid", size=0.02, material=ST, clamp=dict(below=("x", 0.0)))
    for i in range(3):
        rows.append(_row("STEP box (mm) tet10 h=0.02", f"f{i + 1} [Hz]", r.modal.freq_hz[i], ti[i][0],
                         nodes=r.info["nodes"], time_s=r.info["t_total"]))
    # flat STL -> remeshed DKT shells
    pp = os.path.join(tmp, "plate.stl")
    write_stl(pp, np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float), np.array([[0, 1, 2], [0, 2, 3]]))
    r = geometry_modal(pp, 5, kind="shell", thickness=0.01, size=0.06, remesh=True, material=AL,
                       clamp=dict(below=("y", 0.0)))
    lam = plate_lambda(r.modal.omega, 1.0, 0.01, AL["E"], AL["nu"], AL["rho"])
    for i in range(5):
        rows.append(_row("STL plate -> gmsh tri3 h=0.06", f"lambda{i + 1}", lam[i], LEISSA_CFFF[i],
                         nodes=r.info["nodes"], time_s=r.info["t_total"]))
    plate_modal = r.modal
    # NACA wing: solid vs skin shell (no analytic reference; runtimes + frequencies)
    V, F = naca_wing_surface(1.2, 0.4, 0.25, 10.0, 0.12, 20, 10)
    pw = os.path.join(tmp, "wing.stl")
    write_stl(pw, V, F)
    for kind, kw in (("solid", dict(element="tet10", size=0.03)), ("solid", dict(element="tet4", size=0.03)),
                     ("shell", dict(thickness=0.003))):
        r = geometry_modal(pw, 4, kind=kind, material="aluminium_7075", clamp=dict(below=("y", 0.0)), **kw)
        lab = f"NACA wing {kind} {kw.get('element', 'tri3 t=3mm')}"
        for i in range(4):
            rows.append(_row(lab, f"f{i + 1} [Hz]", r.modal.freq_hz[i], nodes=r.info["nodes"],
                             time_s=r.info["t_total"]))
    return rows, plate_modal


def wing_flutter(plot=True):
    sys.path.insert(0, os.path.join(ROOT, "examples"))
    import stl_wing_modal as W
    o = W.run(plot=plot, verbose=False)
    c = o["critical"]
    rows = [_row("STL wing tet10", f"f{i + 1} [Hz]", f) for i, f in enumerate(o["freq_hz"])]
    rows.append(_row("STL wing PK (M=0.3, sea level)", "flutter V [m/s]",
                     c.velocity if c else np.nan, freq_hz=c.freq_hz if c else np.nan,
                     time_s=o["t_modal"] + o["t_flutter"]))
    return rows


# ----------------------------------------------------------------- driver
def _print(name, rows):
    print(f"\n== {name}")
    for r in rows:
        s = f"  {r['case']:<38s} {r['quantity']:<24s} {r['computed']:12.5g}"
        if "reference" in r:
            s += f"  ref {r['reference']:12.5g}  err {r['err_pct']:+7.3f} %"
        extra = {k: v for k, v in r.items() if k not in ("case", "quantity", "computed", "reference", "err_pct")}
        if extra:
            s += "  " + ", ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                                  for k, v in extra.items())
        print(s)


def run(verbose: bool = True, plot: bool = True, flutter: bool = True) -> dict:
    t0 = time.perf_counter()
    out = {}
    out["shell_leissa"], curves, plate_mm = shell_leissa(conv=plot)
    out["solid_beam"], beam_mm = solid_beam()
    out["rigid_modes"] = rigid_modes()
    out["openradioss"] = openradioss()
    out["ansys"] = ansys()
    out["calculix"] = calculix(flutter)
    out["gmsh"] = gmsh_mesher(flutter)
    out["geometry"], stl_plate = geometry(plot)
    if flutter:
        out["wing_flutter"] = wing_flutter(plot)
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from flaap.structures.geometry_modal import plot_modes
        os.makedirs(FIG, exist_ok=True)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for j, kind in enumerate(("quad4", "tri3")):
            ns = sorted(n for k, n in curves if k == kind)
            for i in range(5):
                err = [100 * (curves[(kind, n)][i] - LEISSA_CFFF[i]) / LEISSA_CFFF[i] for n in ns]
                ax[j].plot(ns, err, "o-", label=f"lambda{i + 1}")
            ax[j].axhline(0, color="k", lw=0.7)
            ax[j].set_title(f"{kind}: error vs Leissa (CFFF, a/h=100)")
            ax[j].set_xlabel("elements per side")
            ax[j].set_ylabel("error [%]")
            ax[j].grid(alpha=0.3)
            ax[j].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "structures_shell_convergence.png"), dpi=120)
        plt.close(fig)
        plot_modes(beam_mm, range(3), path=os.path.join(FIG, "structures_beam_modes.png"),
                   suptitle="tet10 16x2x2 steel cantilever", view=(20, -60))
        plt.close("all")
        plot_modes(stl_plate, range(5), path=os.path.join(FIG, "structures_plate_modes.png"),
                   suptitle="STL plate -> gmsh DKT/Allman shells (CFFF)", view=(30, -60))
        plt.close("all")
    out["runtime_s"] = time.perf_counter() - t0
    if verbose:
        for k, v in out.items():
            if isinstance(v, list):
                _print(k, v)
        print(f"\ntotal runtime {out['runtime_s']:.1f} s")
    return out


if __name__ == "__main__":
    run(flutter="--no-flutter" not in sys.argv)
