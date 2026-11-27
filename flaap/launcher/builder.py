"""Build structure + aerodynamic model from resolved flat-config values."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

MATERIALS = {   # E [Pa], nu, rho [kg/m^3]
    "aluminum": (71.0e9, 0.33, 2700.0),
    "steel": (210.0e9, 0.30, 7850.0),
    "titanium": (114.0e9, 0.34, 4430.0),
    "cfrp_quasi": (60.0e9, 0.30, 1600.0),
}


@dataclass
class Case:
    """Everything a route needs."""
    values: dict
    modal: object                      # ModalModel or FlexibleAircraft
    mesh: object                       # AeroMesh
    modes: object                      # AeroModes (structural [+ control] columns)
    bref: float
    mach: float
    names: list[str]
    controls: list[str] = field(default_factory=list)   # control columns appended to `modes`
    struct_vis: dict = field(default_factory=dict)       # nodes / quads / edges / tris for viewers
    extra: dict = field(default_factory=dict)

    @property
    def M(self):
        return self.modal.M

    @property
    def K(self):
        return self.modal.K

    def C(self):
        c = self.modal.C
        return c() if callable(c) else c

    @property
    def nstruct(self) -> int:
        return self.M.shape[0]


def _path(values, key, base: Path) -> Path:
    p = Path(values[key])
    return p if p.is_absolute() else (base / p)


def _material(v):
    from ..structures.plate import OrthotropicMaterial
    name = v.get("material", "aluminum")
    E, nu, rho = MATERIALS.get(name, MATERIALS["aluminum"])
    E = float(v.get("young", E))
    nu = float(v.get("poisson", nu))
    rho = float(v.get("material_density", rho))
    return OrthotropicMaterial.isotropic(E, nu, rho), (E, nu, rho)


def _clamp(v):
    c = v.get("clamp", "y, 0.0")
    if c in (None, "none", False):
        return None
    if isinstance(c, list) and len(c) == 2:
        return str(c[0]), float(c[1])
    return "y", 0.0


def surfaces_from_values(v: dict):
    """``surfaceN`` / ``controlN`` keys -> list of Surface (None if none given)."""
    from ..geometry import Surface
    idx = sorted(int(k[7:]) for k in v if k.startswith("surface") and k[7:].isdigit())
    if not idx:
        return None
    surfs = []
    for i in idx:
        s = v[f"surface{i}"]
        name = str(s[0])
        surfs.append(Surface(name, (float(s[1]), float(s[2]), float(s[3])), float(s[4]),
                             (float(s[5]), float(s[6]), float(s[7])), float(s[8]),
                             nchord=int(s[9]), nspan=int(s[10])))
    add_controls(surfs, v)
    return surfs


def add_controls(surfs, v):
    from ..geometry import ControlSurface
    for k in sorted(k for k in v if k.startswith("control") and k[7:].isdigit()):
        c = v[k]
        for s in surfs:
            if s.name == str(c[1]):
                s.control_surfaces.append(ControlSurface(str(c[0]), float(c[2]), (float(c[3]), float(c[4]))))


def _vis_from_fe(fe, nodes):
    vis = {"nodes": np.asarray(nodes).tolist()}
    if fe is not None:
        el = getattr(fe, "elements", {}) or {}
        for key, out in (("quad4", "quads"), ("quads", "quads"), ("tri3", "tris"), ("tris", "tris"),
                         ("tet4", "tets"), ("edges", "edges"), ("beam2", "edges")):
            if key in el:
                vis[out] = np.asarray(el[key]).astype(int).tolist()
    if len(vis) == 1:                                   # beam / stick: chain the nodes
        n = len(nodes)
        vis["edges"] = [[i, i + 1] for i in range(n - 1)]
    return vis


def _surface_tris(vis):
    """Boundary triangles of tets for display."""
    if "tets" in vis and "tris" not in vis:
        tets = np.asarray(vis["tets"])
        faces = np.vstack([tets[:, [0, 1, 2]], tets[:, [0, 1, 3]], tets[:, [0, 2, 3]], tets[:, [1, 2, 3]]])
        key = np.sort(faces, axis=1)
        _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        vis["tris"] = faces[cnt[inv.ravel()] == 1].tolist()
        del vis["tets"]
    return vis


def _vis_from_mesh(mm):
    """Viewer connectivity from an FEMesh attached to a ModalModel (shell / solid)."""
    mesh = (getattr(getattr(mm, "fe", None), "elements", None) or {}).get("mesh")
    if mesh is None:
        return _vis_from_fe(None, mm.nodes)
    from ..structures.geometry_modal import surface_triangulation
    return {"nodes": mm.nodes.tolist(), "tris": np.asarray(surface_triangulation(mesh)).astype(int).tolist()}


def _auto_surfaces(nodes, v):
    """Default lifting surface for file-based structures: the planform of the nodes.

    Root/tip chords come from the x-extent of the nodes near the minimum/maximum
    spanwise coordinate (y), so a swept or tapered wing gets a matching trapezoid.
    """
    from ..geometry import Surface
    X = np.asarray(nodes, float)
    y0, y1 = X[:, 1].min(), X[:, 1].max()
    if y1 - y0 <= 1e-9:
        return None
    band = 0.05 * (y1 - y0)
    r = X[X[:, 1] <= y0 + band]
    t = X[X[:, 1] >= y1 - band]
    zc = float(np.median(X[:, 2]))
    return [Surface("wing", (r[:, 0].min(), y0, zc), float(np.ptp(r[:, 0])) or 1.0,
                    (t[:, 0].min(), y1, zc), float(np.ptp(t[:, 0])) or 1.0,
                    nchord=int(v.get("nchord") or 8), nspan=int(v.get("nspan") or 12))]


def _spline_subset(modal, extra, max_nodes=600):
    """Planform-uniform node subset for IPS/TPS splining of large FE models (O(n^3))."""
    n = modal.nodes.shape[0]
    if n <= max_nodes:
        return modal.nodes, modal.dof_vectors()
    from ..structures.geometry_modal import planform_nodes
    excl = extra.get("spline_exclude")
    mask = planform_nodes(modal.nodes, 24, 10, axes=(1, 0), exclude=excl)
    if excl is not None and len(excl):
        mask[np.asarray(excl)[::max(1, len(excl) // 8)]] = True      # a few clamped nodes (u = 0)
    return modal.nodes[mask], modal.shapes[mask].reshape(-1, modal.shapes.shape[2])


# ----------------------------------------------------------------- structures
_AUTO_SURFACE_STRUCTURES = ("ansys", "openradioss", "openradioss_eig", "calculix", "gmsh", "geometry", "stl")


def build_structure(v: dict, base: Path, log=print):
    """Returns (modal, default_surfaces, default_bref, spline_kind, vis, extra).

    The structure is built in its own axes, then re-oriented by ``struct_mirror`` / ``struct_rotate`` /
    ``struct_translate`` (see :mod:`flaap.structures.transform`) into the aerodynamic frame.
    """
    from ..structures.transform import transform_from_values, transform_modal, transform_surface
    modal, surfs, bref0, spl, vis, extra = _build_structure_raw(v, base, log)
    T = transform_from_values(v)
    if T.identity:
        return modal, surfs, bref0, spl, vis, extra
    s = v.get("structure", "goland")
    if s == "demo_aircraft" or spl == "demo":
        raise ValueError("struct_mirror / struct_rotate / struct_translate do not apply to demo_aircraft "
                         "(free-flying model with its own aero)")
    if s == "legacy" and "surface1" not in v:
        raise ValueError("a legacy deck brings its own aero panels: give surface1 ... to re-orient it")
    modal = transform_modal(modal, T)
    vis = dict(vis)
    vis["nodes"] = modal.nodes.tolist()
    extra = dict(extra, transform=T)
    log(f"  structure transform: {T.describe()}")
    if surfs is not None and "surface1" not in v:
        if s in _AUTO_SURFACE_STRUCTURES:               # planform recomputed from the re-oriented nodes
            surfs = _auto_surfaces(modal.nodes, v)
        else:                                           # preset planform follows the transform
            surfs = [transform_surface(sf, T) for sf in surfs]
    return modal, surfs, bref0, spl, vis, extra


def _build_structure_raw(v: dict, base: Path, log=print):
    s = v.get("structure", "goland")
    nm = int(v.get("nmodes", 4))
    zeta = float(v.get("damping", 0.0))
    nchord, nspan = v.get("nchord"), v.get("nspan")
    extra = {}
    if s in ("goland", "goland_ctrl"):
        from ..models.goland import CHORD, goland_modal, goland_surface
        mm = goland_modal(nm, nelem=int(v.get("nelem", 20)), damping=zeta)
        surf = goland_surface(int(nchord or 6), int(nspan or 12))
        if s == "goland_ctrl":
            from ..geometry import ControlSurface
            surf.control_surfaces.append(ControlSurface("aileron", 0.8, (0.6, 1.0)))
        mm.labels = mm.labels or [f"mode {i + 1}" for i in range(nm)]
        return mm, [surf], CHORD / 2, "beam", _vis_from_fe(None, mm.nodes), extra
    if s == "agard445":
        from ..models.agard445 import B_REF, EXPERIMENT_FREQ_HZ, agard445_modal, agard445_surface
        mm = agard445_modal(nm)
        if v.get("experimental_frequencies", False):
            n = min(4, nm)
            mm.omega[:n] = 2 * np.pi * EXPERIMENT_FREQ_HZ[:n]
        mm.damping[:] = zeta
        return mm, [agard445_surface(int(nchord or 8), int(nspan or 10))], B_REF, "ips", \
            _vis_from_fe(mm.fe, mm.nodes), extra
    if s == "plate":
        from ..geometry import Surface
        from ..structures.plate import build_plate_model, plate_mesh_trapezoid
        rle = v.get("plate_root_le", [0.0, 0.0, 0.0])
        tle = v.get("plate_tip_le", [0.2, 1.0, 0.0])
        rc = float(v.get("plate_root_chord", 0.5))
        tc = float(v.get("plate_tip_chord", 0.3))
        nodes, quads = plate_mesh_trapezoid(rle, rc, tle, tc, int(v.get("plate_nchord", 8)),
                                            int(v.get("plate_nspan", 12)))
        mat, (E, nu, rho) = _material(v)
        cl = _clamp(v)
        fe = build_plate_model(nodes, quads, float(v.get("thickness", 0.002)), mat,
                               clamp_y=None if cl is None else cl[1], name="plate")
        mm = fe.modes(nm, name="plate", damping=zeta)
        mm.labels = [f"mode {i + 1}" for i in range(nm)]
        r3 = list(rle) + [0.0] * (3 - len(rle))
        t3 = list(tle) + [0.0] * (3 - len(tle))
        surf = Surface("plate", tuple(r3), rc, tuple(t3), tc, nchord=int(nchord or 8), nspan=int(nspan or 12))
        extra["material"] = {"E": E, "nu": nu, "rho": rho}
        return mm, [surf], rc / 2, "ips", _vis_from_fe(fe, nodes), extra
    if s == "demo_aircraft":
        from ..models.demo_aircraft import demo_aircraft
        ac, mesh, modes = demo_aircraft(nm, damping=zeta if zeta > 0 else 0.02)
        extra["demo"] = (mesh, modes)
        vis = _vis_from_fe(None, ac.nodes)
        return ac, None, 0.9144, "demo", vis, extra
    if s == "legacy":
        from ..io.legacy import read_case
        case = read_case(str(_path(v, "structure_file", base)))
        mm = case.modal(float(v.get("length_scale", 1.0)), v.get("nmodes"))
        mm.damping[:] = zeta
        extra["legacy_case"] = case
        return mm, None, None, "legacy", _vis_from_fe(None, mm.nodes), extra
    if s == "nastran":
        from ..io.nastran import read_f06_modal
        bulk = _path(v, "bulk_file", base) if "bulk_file" in v else None
        mm = read_f06_modal(str(_path(v, "structure_file", base)), None if bulk is None else str(bulk),
                            nmodes=nm, length_scale=float(v.get("length_scale", 1.0)),
                            mass_scale=float(v.get("mass_scale", 1.0)))
        mm.damping[:] = zeta
        return mm, None, None, "ips", {"nodes": mm.nodes.tolist()}, extra
    if s == "ansys":
        from ..io.ansys import read_ansys_modal
        kw = {}
        if "bulk_file" in v:                                 # .cdb / NLIST node file
            kw["nodes"] = str(_path(v, "bulk_file", base))
        mm = read_ansys_modal(str(_path(v, "structure_file", base)), nmodes=nm,
                              length_scale=float(v.get("length_scale", 1.0)),
                              mass_scale=float(v.get("mass_scale", 1.0)), **kw)
        mm.damping[:] = zeta
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", _vis_from_fe(None, mm.nodes), extra
    if s == "openradioss":
        from ..io.openradioss import radioss_modal
        mm = radioss_modal(str(_path(v, "structure_file", base)), nmodes=nm, damping=zeta)
        if float(v.get("length_scale", 1.0)) != 1.0:
            mm.nodes = mm.nodes * float(v["length_scale"])
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", _vis_from_mesh(mm), extra
    if s == "openradioss_eig":                               # Radioss /EIG modes (anim_to_vtk files)
        from ..io.openradioss_results import read_radioss_eigen
        opt = {k: str(_path(v, k, base)) for k in ("mode_files", "bulk_file", "frequency_file") if k in v}
        units = v.get("radioss_units")
        mm = read_radioss_eigen(str(_path(v, "structure_file", base)), vtk_files=opt.get("mode_files"),
                                deck=opt.get("bulk_file"), listing=opt.get("frequency_file"), nmodes=nm,
                                units=None if units is None else tuple(str(u) for u in units), damping=zeta)
        if float(v.get("length_scale", 1.0)) != 1.0:
            mm.nodes = mm.nodes * float(v["length_scale"])
        mm.gen_mass = mm.gen_mass * float(v.get("mass_scale", 1.0))
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", _vis_from_mesh(mm), extra
    if s == "calculix":                                      # CalculiX *FREQUENCY .frd (+ .dat)
        from ..io.calculix import frd_surface, read_calculix_modal, read_frd
        path = str(_path(v, "structure_file", base))
        frd = read_frd(path)
        ff = str(_path(v, "frequency_file", base)) if "frequency_file" in v else None
        mm = read_calculix_modal(path, nmodes=nm, frequency_file=ff, frd=frd,
                                 length_scale=float(v.get("length_scale", 1.0)),
                                 mass_scale=float(v.get("mass_scale", 1.0)), damping=zeta)
        surf = frd_surface(frd, mm.node_ids)
        vis = {"nodes": mm.nodes.tolist()}
        if len(surf["tris"]):
            vis["tris"] = surf["tris"].tolist()
        if len(surf["edges"]):
            vis["edges"] = surf["edges"].tolist()
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", vis, extra
    if s == "gmsh":                                          # gmsh .geo / .msh / CAD -> FE mesh
        from ..structures.fem import Material
        from ..structures.gmsh_mesher import gmsh_modal
        _, (E, nu, rho) = _material(v)
        cg = v.get("clamp_group")
        cg = None if cg in (None, "none") else ([str(c) for c in cg] if isinstance(cg, list) else [str(cg)])
        cl = None if cg else _clamp(v)
        gt = v.get("group_thickness")
        gt = gt if isinstance(gt, list) or gt is None else [gt]
        res = gmsh_modal(str(_path(v, "structure_file", base)), nm,
                         material=Material(E, nu, rho, v.get("material", "custom")),
                         thickness=float(v.get("thickness", 0.002)),
                         group_thickness={str(a): float(b) for a, b in zip(gt[0::2], gt[1::2])} if gt else None,
                         clamp=None if cl is None else dict(below=cl), clamp_group=cg,
                         length_scale=float(v.get("length_scale", 1.0)), damping=zeta,
                         size=v.get("mesh_size"), order=v.get("mesh_order"), recombine=v.get("recombine"),
                         dim=v.get("mesh_dim"), algorithm=v.get("mesh_algorithm"),
                         save=str(_path(v, "mesh_save", base)) if v.get("mesh_save") else None)
        mm = res.modal
        mm.labels = [f"mode {i + 1}" for i in range(nm)]
        gm = res.info.pop("gmsh_mesh")
        extra["material"] = {"E": E, "nu": nu, "rho": rho}
        extra["fe_info"] = {k: (float(x) if isinstance(x, (int, float, np.floating)) else x)
                            for k, x in res.info.items()}
        extra["spline_exclude"] = res.clamped
        log(f"  gmsh FE: {gm.info['nodes']} nodes, {gm.info['elements']} {gm.info['kinds']} elements "
            f"({'meshed' if gm.info.get('generated') else 'read'}), groups {sorted(gm.groups)}, "
            f"{res.info['clamped_nodes']} clamped, mass {res.info['mass']:.4g} kg")
        vis = {"nodes": mm.nodes.tolist(), "tris": np.asarray(res.surface).astype(int).tolist()}
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", vis, extra
    if s in ("geometry", "stl"):
        from ..structures.fem import Material
        from ..structures.geometry_modal import geometry_modal
        _, (E, nu, rho) = _material(v)
        cl = _clamp(v)
        kind = v.get("fe_type", "shell")
        res = geometry_modal(str(_path(v, "structure_file", base)), nm, kind=kind,
                             element=v.get("element", "tet10"), size=v.get("mesh_size"),
                             thickness=float(v.get("thickness", 0.002)) if kind == "shell" else None,
                             material=Material(E, nu, rho, v.get("material", "custom")),
                             clamp=None if cl is None else dict(below=cl),
                             length_scale=float(v.get("length_scale", 1.0)), damping=zeta)
        mm = res.modal
        mm.labels = [f"mode {i + 1}" for i in range(nm)]
        extra["material"] = {"E": E, "nu": nu, "rho": rho}
        extra["fe_info"] = {k: (float(x) if isinstance(x, (int, float, np.floating)) else x)
                            for k, x in res.info.items()}
        extra["spline_exclude"] = res.clamped
        log(f"  geometry FE: {res.info['nodes']} nodes, {res.info['elements']} {res.info['kinds']} elements, "
            f"mass {res.info['mass']:.4g} kg")
        vis = {"nodes": mm.nodes.tolist(), "tris": np.asarray(res.surface).astype(int).tolist()}
        return mm, _auto_surfaces(mm.nodes, v), None, "ips", vis, extra
    raise ValueError(f"unknown structure {s!r}")


def build_case(v: dict, base: Path, *, need_aero=True, with_controls=False, log=print) -> Case:
    from ..aero.gaf import AeroModes, control_aero_mode
    from ..geometry import build_mesh
    from ..structures.spline import beam_spline, surface_spline
    modal, surfs, bref0, spl, vis, extra = build_structure(v, base, log)
    names = list(getattr(modal, "names", None) or modal.labels or [f"mode {i + 1}" for i in range(modal.M.shape[0])])
    if not need_aero:
        return Case(v, modal, None, None, bref0 or 1.0, float(v.get("mach", 0.0)), names, [], vis, extra)
    mach = float(v.get("mach", 0.0))
    controls = []
    if spl == "demo":                                       # demo aircraft: modes already built
        mesh, modes = extra["demo"]
        n = modal.M.shape[0]
        controls = list(modes.names[n:])
        if not with_controls:
            modes = modes.select(range(n))
            controls = []
        return Case(v, modal, mesh, modes, float(v.get("bref", bref0)), mach, names, controls, vis, extra)
    if spl == "legacy" and "surface1" not in v:
        from ..io.legacy import legacy_aero
        am = legacy_aero(extra["legacy_case"], modal, v.get("spline", "auto") if v.get("spline") != "auto"
                         else "legacy", v.get("mach"))
        return Case(v, modal, am.mesh, am.modes, am.bref, am.mach, names, [], vis, extra)
    user = surfaces_from_values(v)
    if user is not None:
        surfs = user
    elif surfs is not None:
        add_controls(surfs, v)
    mesh = build_mesh(surfs, v.get("symmetry", "symmetric"))
    kind = v.get("spline", "auto")
    kind = spl if kind == "auto" else kind
    if kind == "beam":
        sp = beam_spline(mesh, modal.nodes)
        vecs = modal.dof_vectors()
    else:
        snodes, vecs = _spline_subset(modal, extra)
        if snodes.shape[0] < modal.nodes.shape[0]:
            log(f"  spline: {snodes.shape[0]} of {modal.nodes.shape[0]} structural nodes (planform subset)")
        sp = surface_spline(mesh, snodes, method="ips" if kind == "legacy" else kind)
    modes = AeroModes.from_spline(sp, vecs, names)
    if with_controls:
        for cname in mesh.control:
            modes = modes + control_aero_mode(mesh, cname)
            controls.append(cname)
    bref = float(v["bref"]) if "bref" in v else bref0
    if bref is None:                                    # half the mean geometric chord
        bref = 0.25 * float(np.mean([sf.root_chord + sf.tip_chord for sf in surfs]))
    return Case(v, modal, mesh, modes, bref, mach, names, controls, vis, extra)


def k_values(v: dict) -> np.ndarray:
    if "k_list" in v:
        kl = v["k_list"]
        return np.atleast_1d(np.asarray(kl if isinstance(kl, list) else [kl], float))
    kmax, nk = float(v.get("kmax", 2.0)), int(v.get("nk", 20))
    return np.r_[0.0, np.geomspace(min(0.01, kmax / 100), kmax, nk - 1)]


def density(v: dict) -> float:
    from ..atmosphere import isa
    if "density" in v:
        return float(v["density"])
    return isa(float(v.get("altitude", 0.0))).density


def velocity_range(v: dict, key="velocities", default=(10.0, 300.0, 146)) -> np.ndarray:
    r = v.get(key, list(default))
    return np.linspace(float(r[0]), float(r[1]), int(r[2]))


def tsd_thickness(name):
    from ..aero import tsd
    name = str(name or "flat").lower()
    if name in ("flat", "none", "0"):
        return None
    if name == "naca65a004":
        return tsd.naca65a004()
    if name.startswith("naca") and name[4:].isdigit():
        return tsd.naca4_thickness(int(name[-2:]) / 100.0)
    if name.startswith("biconvex"):
        return tsd.parabolic_arc(float(name[8:]) / 100.0)
    raise ValueError(f"unknown tsd_thickness {name!r}")


def tsd_camber(name):
    """``flat`` / ``none`` -> None; ``naca2412`` (or ``naca24``) -> NACA 4-digit mean line."""
    from ..aero.tsd_lifting import naca4_camber
    name = str(name or "flat").lower()
    if name in ("flat", "none", "0"):
        return None
    d = name[4:] if name.startswith("naca") else ""
    if d.isdigit() and len(d) in (2, 4):
        m, p = int(d[0]) / 100.0, int(d[1]) / 10.0
        if m == 0:
            return None
        if p == 0:
            raise ValueError(f"tsd_camber {name!r}: camber without a camber position")
        return naca4_camber(m, p)
    raise ValueError(f"unknown tsd_camber {name!r}")


def tsd_alpha(v) -> float:
    """TSD mean-flow incidence [deg]: ``tsd_alpha``, else ``alpha`` (default 0)."""
    return float(v.get("tsd_alpha", v.get("alpha", 0.0)))


def tsd_lifting(v) -> bool:
    """True when the TSD mean flow is non-symmetric (incidence or camber)."""
    return tsd_alpha(v) != 0.0 or tsd_camber(v.get("tsd_camber")) is not None
