"""Route specifications: which keys each ``model`` / ``mode`` accepts, defaults
and validation rules (AI-CAE4ALL ``MethodSpec`` analogue)."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from .diagnostics import DiagnosticReport, Severity
from .parser import ParsedConfig

# ----------------------------------------------------------------- key catalogue
# key -> (description, type) ; type in {"str","int","float","bool","list","vec","path","any"}
COMMON_KEYS: dict[str, tuple[str, str]] = {
    "model": ("analysis route (modal, aero, flutter, gust, trim, maneuver, tsd, uvlm, freeplay, ase)", "str"),
    "mode": ("operation inside the route", "str"),
    "title": ("case title (keeps case)", "str"),
    "output_dir": ("output directory (default output/<config name>)", "path"),
    "plots": ("write PNG figures", "bool"),
    "verbose": ("print progress", "bool"),
    "output_json": ("also write every result file as JSON (default: text only)", "bool"),
}
STRUCTURE_KEYS: dict[str, tuple[str, str]] = {
    "structure": ("goland | goland_ctrl | agard445 | demo_aircraft | plate | nastran | ansys | openradioss | "
                  "openradioss_eig | calculix | geometry | stl | gmsh | legacy", "str"),
    "structure_file": ("modal result file: .f06 (nastran), ANSYS modal text, OpenRadioss starter deck, "
                       "Radioss /EIG mode VTK folder / <run>A001.vtk / engine listing (openradioss_eig), "
                       "CalculiX *FREQUENCY .frd (calculix), STL/STEP (geometry), gmsh .geo/.msh/STEP (gmsh), "
                       "legacy input directory", "path"),
    "bulk_file": ("Nastran bulk file with GRID cards (if not echoed in the .f06); openradioss_eig: "
                  "starter deck *_0000.rad (SI coordinates, units, generalized masses)", "path"),
    "mode_files": ("openradioss_eig: folder or <run>A001.vtk of the anim_to_vtk mode files "
                   "(when structure_file is the listing)", "path"),
    "frequency_file": ("openradioss_eig: engine listing / 'mode frequency' table overriding the VTK TIME "
                       "frequencies (solver units); calculix: .dat with the eigenvalue table "
                       "(default: the .dat next to the .frd)", "path"),
    "radioss_units": ("openradioss_eig: mass, length, time units of the run when no starter deck gives "
                      "them (e.g. kg, mm, ms; 't' for tonne; default SI)", "list"),
    "nmodes": ("number of structural modes", "int"),
    "damping": ("modal damping ratio zeta", "float"),
    "length_scale": ("multiply coordinates (e.g. 0.0254 for inch -> m)", "float"),
    "mass_scale": ("multiply generalized masses", "float"),
    "struct_mirror": ("mirror the structure about the plane normal to x | y | z (y = left <-> right); "
                      "applied first", "str"),
    "struct_rotate": ("rotate the structure: axis, degrees (right-hand rule), e.g. 'z, 90'; pairs may be "
                      "chained ('z, 90, x, 180'); applied after the mirror", "list"),
    "struct_translate": ("shift the structure dx, dy, dz [m] (after length_scale, mirror and rotation)", "vec"),
    "experimental_frequencies": ("AGARD 445.6: replace FE frequencies by measured ones", "bool"),
    "nelem": ("Goland beam elements", "int"),
    # plate / geometry FE
    "plate_root_le": ("plate root leading edge x, y, z", "vec"),
    "plate_root_chord": ("plate root chord", "float"),
    "plate_tip_le": ("plate tip leading edge x, y, z", "vec"),
    "plate_tip_chord": ("plate tip chord", "float"),
    "plate_nchord": ("plate FE elements chordwise", "int"),
    "plate_nspan": ("plate FE elements spanwise", "int"),
    "thickness": ("plate / shell thickness [m]", "float"),
    "material": ("aluminum | steel | titanium | cfrp_quasi | custom", "str"),
    "young": ("Young's modulus [Pa]", "float"),
    "poisson": ("Poisson ratio", "float"),
    "material_density": ("material density [kg/m^3]", "float"),
    "clamp": ("clamp condition: 'y, 0.0' (axis, coordinate) or none", "any"),
    "fe_type": ("geometry FE: shell (midsurface) | solid (tetra)", "str"),
    "mesh_size": ("target FE element size [m] for geometry meshing", "float"),
    "element": ("solid element for geometry FE: tet10 | tet4", "str"),
    # gmsh mesher
    "clamp_group": ("gmsh: physical group name(s) whose nodes are clamped (replaces clamp)", "list"),
    "mesh_order": ("gmsh: element order 1 (tri3/quad4/tet4/hex8) | 2 (tet10)", "int"),
    "recombine": ("gmsh: recombine surface triangles into quad4 (Mesh.RecombineAll)", "bool"),
    "mesh_dim": ("gmsh: mesh dimension 2 (shells) | 3 (solids); default: highest entity dimension", "int"),
    "mesh_algorithm": ("gmsh: Mesh.Algorithm (2-D: 1 MeshAdapt, 5 Delaunay, 6 Frontal-Delaunay, 8 quads)", "int"),
    "group_thickness": ("gmsh: per physical group shell thickness [m]: name, t, name, t, ...", "list"),
    "mesh_save": ("gmsh: also write the FE mesh to this .msh file (relative to the config)", "path"),
}
AERO_KEYS: dict[str, tuple[str, str]] = {
    "mach": ("free-stream Mach number", "float"),
    "bref": ("reference semi-chord b [m] (k = w b / V)", "float"),
    "symmetry": ("none | symmetric | antisymmetric (about y = 0)", "str"),
    "kmax": ("largest reduced frequency in the GAF table", "float"),
    "nk": ("number of reduced frequencies", "int"),
    "k_list": ("explicit reduced frequency list", "list"),
    "spline": ("ips | tps | rbf | beam | legacy", "str"),
    "aero_method": ("dlm | tsd | uvlm (linearised UVLM, M = 0) | machbox | piston (supersonic, M > 1)", "str"),
    "machbox_nx": ("Mach box: boxes along the root chord (default 60)", "int"),
    "piston_order": ("piston theory: 1 (4/M) or Van Dyke 2nd order (default 2)", "int"),
    "wake_chords": ("UVLM: wake length in reference chords (default 30)", "float"),
    "tsd_thickness": ("TSD mean flow thickness: flat | naca65a004 | nacaMPXX (thickness XX %) | biconvexXX", "str"),
    "tsd_camber": ("TSD mean flow camber line: flat | nacaMPXX / nacaMP (NACA 4-digit mean line, e.g. naca2412)",
                   "str"),
    "tsd_alpha": ("TSD mean flow angle of attack [deg]: nonlinear lifting mean flow for the GAFs / flutter "
                  "(default 0 = symmetric thickness flow)", "float"),
    "nchord": ("default chordwise boxes for preset surfaces", "int"),
    "nspan": ("default spanwise boxes for preset surfaces", "int"),
}
FLIGHT_KEYS: dict[str, tuple[str, str]] = {
    "density": ("air density [kg/m^3] (overrides altitude)", "float"),
    "altitude": ("ISA altitude [m]", "float"),
    "velocities": ("velocity sweep start, stop, count [m/s]", "list"),
    "velocity": ("analysis velocity [m/s]", "float"),
    "nlag": ("Roger RFA lag terms", "int"),
    "methods": ("flutter methods for mode 'all': pk, k, p", "list"),
}
GUST_KEYS = {
    "gust_gradients": ("gust gradient sweep H start, stop, count [m]", "list"),
    "fg": ("CS-25 flight profile alleviation factor", "float"),
    "turbulence_length": ("von Karman scale length [m]", "float"),
}
TRIM_KEYS = {
    "nz": ("load factor for trim", "float"),
    "trim_free": ("free trim variables (e.g. elevator)", "list"),
    "control": ("control surface used for trim / manoeuvre / ASE", "str"),
    "nz_limit": ("manoeuvre limit load factor", "float"),
    "amplitude_deg": ("abrupt manoeuvre control amplitude [deg]", "float"),
    "t_end": ("simulation end time [s]", "float"),
    "dt": ("time step [s]", "float"),
}
TSD_KEYS = {
    "tsd_nx": ("TSD grid points chordwise on the wing", "int"),
    "tsd_ny": ("TSD grid points spanwise on the wing", "int"),
    "alpha": ("mean incidence [deg] of the lifting TSD mean flow (same as tsd_alpha)", "float"),
    "alphas": ("steady TSD incidence sweep start, stop, count [deg]: CL and shock positions vs alpha", "list"),
}
UVLM_KEYS = {
    "alpha": ("incidence [deg] (impulsive start / mean angle of the motion)", "float"),
    "alphas": ("steady VLM incidence sweep start, stop, count [deg]", "list"),
    "velocity": ("free-stream speed [m/s]", "float"),
    "x_ref": ("pitch axis / moment reference x [m] (default 1/4 of the mean chord)", "float"),
    "pitch_amplitude_deg": ("harmonic pitch amplitude [deg]", "float"),
    "plunge_amplitude": ("harmonic plunge amplitude [m] (positive up)", "float"),
    "reduced_frequency": ("harmonic motion reduced frequency k = w bref / V", "float"),
    "phase_deg": ("pitch phase lead over plunge [deg]", "float"),
    "ncycles": ("harmonic: number of motion cycles", "float"),
    "gust_shape": ("sharp | 1-cos", "str"),
    "gust_velocity": ("gust vertical velocity w0 [m/s]", "float"),
    "gust_length": ("1-cos gust gradient H [m] (half wavelength)", "float"),
    "free_wake": ("force-free (relaxed) wake roll-up instead of a flat prescribed wake", "bool"),
    "wake_length": ("truncate the time-marching wake after this many reference chords", "float"),
    "chord_steps": ("time steps per trailing-edge panel passage (dt = dx_TE / (V chord_steps)); default 1", "int"),
    "t_end": ("simulation end time [s]", "float"),
    "dt": ("time step [s] (default dx_TE / V)", "float"),
}
NONLINEAR_KEYS = {
    "nonlinearity": ("hinge spring: freeplay | bilinear | preloaded", "str"),
    "freeplay_deg": ("freeplay half-gap delta0 [deg] (section 1.0, goland_flap 0.5)", "float"),
    "inner_stiffness_ratio": ("bilinear: stiffness inside the gap / K_h", "float"),
    "preload_deg": ("preloaded freeplay: gap centre offset [deg]", "float"),
    "tm_velocities": ("describing_function: speeds for the time-marching check (default 3 in the LCO band)", "list"),
    "initial_deflection_deg": ("initial hinge rotation [deg] (default 3 delta0)", "float"),
    "t_end": ("simulation end time [s]", "float"),
    "dt": ("time step / output sampling [s]", "float"),
}
ASE_KEYS = {
    "actuator_frequency": ("actuator bandwidth [Hz] (2nd order)", "float"),
    "actuator_damping": ("actuator damping ratio", "float"),
    "sensors": ("section: pitch, pitch_rate, plunge, plunge_rate | goland_flap: root_bending, tip_accel", "list"),
    "feedback_gains": ("static output feedback, one gain per sensor (u_cmd = K y)", "list"),
    "gust_amplitude": ("GLA: 1-cos gust design velocity Uds [m/s]", "float"),
}
NL_PLANT_KEYS = {
    "structure": ("section | goland_flap", "str"),
    "nmodes": ("goland_flap: number of wing modes", "int"),
    "damping": ("goland_flap: modal damping ratio (the section uses Conner damping)", "float"),
    "density": ("air density [kg/m^3] (default: the plant's)", "float"),
    "velocities": ("speed sweep start, stop, count [m/s]", "list"),
    "velocity": ("analysis speed [m/s]", "float"),
}
INDEXED_PATTERNS = {
    "surface": (re.compile(r"^surface(\d+)$"),
                "name, rlex, rley, rlez, root_chord, tlex, tley, tlez, tip_chord, nchord, nspan"),
    "control": (re.compile(r"^control(\d+)$"), "name, surface_name, hinge, eta0, eta1"),
    "station": (re.compile(r"^station(\d+)$"), "name, px, py, pz, ax, ay, az"),
}


@dataclass(frozen=True)
class RouteSpec:
    name: str
    description: str
    modes: dict[str, str]
    key_groups: tuple[dict, ...]
    defaults: dict[str, Any] = field(default_factory=dict)
    required_by_mode: dict[str, tuple[str, ...]] = field(default_factory=dict)
    validators: tuple[Callable[[ParsedConfig, dict, DiagnosticReport], None], ...] = ()

    @property
    def known_keys(self) -> dict[str, tuple[str, str]]:
        out: dict[str, tuple[str, str]] = {}
        for g in (COMMON_KEYS,) + self.key_groups:
            out.update(g)
        return out


STRUCTURES = ("goland", "goland_ctrl", "agard445", "demo_aircraft", "plate", "nastran", "ansys", "openradioss",
              "openradioss_eig", "calculix", "geometry", "stl", "gmsh", "legacy")
NL_PLANTS = ("section", "goland_flap")
BASE_DEFAULTS = {"plots": True, "verbose": True, "structure": "goland", "nmodes": 4, "damping": 0.0,
                 "mach": 0.0, "symmetry": "symmetric", "kmax": 2.0, "nk": 20, "spline": "auto",
                 "aero_method": "dlm", "nlag": 4, "length_scale": 1.0, "mass_scale": 1.0,
                 "material": "aluminum", "thickness": 0.002, "plate_nchord": 8, "plate_nspan": 12,
                 "clamp": "y, 0.0", "fe_type": "shell"}
STRUCT_AERO = (STRUCTURE_KEYS, AERO_KEYS)
FLIGHT = (FLIGHT_KEYS,)


# ----------------------------------------------------------------- validators
def _positive(keys):
    def v(pc: ParsedConfig, vals: dict, rep: DiagnosticReport):
        for k in keys:
            if k in vals and isinstance(vals[k], (int, float)) and not isinstance(vals[k], bool) and vals[k] <= 0:
                rep.add("CFG-VALUE-002", Severity.ERROR, f"{k} must be positive (got {vals[k]}).",
                        field_name=k, location=pc.location(k))
    return v


def _range3(keys):
    def v(pc, vals, rep):
        for k in keys:
            if k not in vals:
                continue
            x = vals[k]
            if not (isinstance(x, list) and len(x) == 3 and all(isinstance(t, (int, float)) for t in x)):
                rep.add("CFG-VALUE-003", Severity.ERROR, f"{k} needs three numbers: start, stop, count.",
                        field_name=k, location=pc.location(k), hint=f"e.g. {k} 50, 300, 126")
            elif x[1] <= x[0] or x[2] < 2:
                rep.add("CFG-VALUE-004", Severity.ERROR, f"{k}: stop must exceed start and count >= 2.",
                        field_name=k, location=pc.location(k))
    return v


def _structure_rules(pc, vals, rep):
    s = vals.get("structure")
    needs_file = {"nastran", "ansys", "openradioss", "openradioss_eig", "calculix", "geometry", "stl",
                  "gmsh", "legacy"}
    if s in needs_file and "structure_file" not in vals:
        rep.add("CFG-REQ-002", Severity.ERROR, f"structure {s} needs structure_file.",
                field_name="structure_file", location=pc.location("structure"),
                hint="structure_file path/to/file (relative to the config)")
    if "structure_file" in vals:
        p = Path(vals["structure_file"])
        if not p.is_absolute():
            p = pc.base_dir / p
        if not p.exists():
            rep.add("CFG-PATH-001", Severity.ERROR, f"structure_file does not exist: {p}",
                    field_name="structure_file", location=pc.location("structure_file"))
    custom_surfaces = any(k.startswith("surface") and k[7:].isdigit() for k in vals)
    if s in ("nastran", "ansys", "openradioss", "openradioss_eig", "calculix", "geometry", "stl", "gmsh",
             "plate") \
            and not custom_surfaces \
            and vals.get("model") not in ("modal",):
        if s == "plate":
            rep.add("CFG-DEFAULT-002", Severity.NOTICE,
                    "No surface1 given: the aero surface follows the plate planform.", field_name="surface1")
        else:
            rep.add("CFG-REQ-003", Severity.ERROR, f"structure {s} needs at least one aero surface (surface1 ...).",
                    field_name="surface1",
                    hint="surface1 wing, 0, 0, 0, 1.0, 0.2, 5.0, 0, 0.6, 8, 12")
    if s in ("nastran", "ansys", "openradioss", "openradioss_eig", "calculix", "geometry", "stl", "gmsh") \
            and "bref" not in vals \
            and vals.get("model") not in ("modal",):
        rep.add("CFG-REQ-004", Severity.ERROR, "bref (reference semi-chord) is required for imported structures.",
                field_name="bref", hint="bref 0.5   # semi-chord in metres")
    if s == "openradioss_eig":
        _radioss_eig_rules(pc, vals, rep)
    if s == "calculix":
        _calculix_rules(pc, vals, rep)
    if s == "gmsh":
        _gmsh_rules(pc, vals, rep)
    if s in ("plate", "geometry", "stl", "gmsh") and vals.get("material") == "custom":
        for k in ("young", "poisson", "material_density"):
            if k not in vals:
                rep.add("CFG-REQ-005", Severity.ERROR, f"material custom needs {k}.", field_name=k)
    if s in ("geometry", "stl") and vals.get("fe_type") == "shell" and "thickness" not in vals:
        rep.add("CFG-DEFAULT-003", Severity.WARNING, "Shell FE of an STL without 'thickness': 2 mm is used.",
                field_name="thickness", promote_in_strict=True)


def _gmsh_rules(pc, vals, rep):
    """gmsh: module installed, file type, clamp_group / group_thickness names present in the file."""
    try:
        import gmsh  # noqa: F401
    except ImportError:
        rep.add("CFG-DEP-001", Severity.ERROR, "structure gmsh needs the gmsh Python module.",
                field_name="structure", location=pc.location("structure"), hint="pip install gmsh")
    exts = (".geo", ".msh", ".step", ".stp", ".iges", ".igs", ".brep", ".brp")
    src = None
    if "structure_file" in vals:
        src = Path(str(vals["structure_file"]))
        src = src if src.is_absolute() else pc.base_dir / src
        if src.suffix.lower() not in exts:
            rep.add("CFG-VALUE-007", Severity.ERROR,
                    f"structure gmsh reads {', '.join(exts)}; got {src.name}.",
                    field_name="structure_file", location=pc.location("structure_file"),
                    hint="for STL use structure stl")
    gt = vals.get("group_thickness")
    if gt is not None:
        gt = gt if isinstance(gt, list) else [gt]
        ok = len(gt) % 2 == 0 and all(isinstance(t, (int, float)) and not isinstance(t, bool) and t > 0
                                      for t in gt[1::2])
        if not ok:
            rep.add("CFG-VALUE-007", Severity.ERROR, "group_thickness needs pairs: name, thickness, ...",
                    field_name="group_thickness", location=pc.location("group_thickness"),
                    hint="group_thickness Skin, 0.002, Spar, 0.004")
    names = []
    for k in ("clamp_group", "group_thickness"):
        if k in vals:
            v = vals[k] if isinstance(vals[k], list) else [vals[k]]
            names += [(k, str(n)) for n in (v[0::2] if k == "group_thickness" else v)]
    if names and src is not None and src.is_file() and src.suffix.lower() in (".geo", ".msh"):
        with open(src, "rb") as fh:
            text = fh.read(4 << 20).decode("latin-1")
        if src.suffix.lower() == ".msh" and "$PhysicalNames" not in text:
            found = set()
        else:
            found = set(re.findall(r'"([^"\n]+)"', text))
        if src.suffix.lower() == ".geo" and re.search(r"(?m)^\s*Include\b", text):
            found = None                                       # names may come from included files
        for k, n in names:
            if found is not None and n not in found:
                rep.add("CFG-PATH-004", Severity.ERROR, f"{k}: physical group \"{n}\" is not defined in {src.name}.",
                        field_name=k, location=pc.location(k),
                        hint=f'add  Physical Curve("{n}") = {{...}};  (or Surface) to the .geo')
    if "clamp_group" not in vals and vals.get("clamp") not in (None, "none", False):
        rep.add("CFG-DEFAULT-001", Severity.NOTICE,
                "gmsh: no clamp_group - the nodes at the clamp plane (default y = 0) are clamped.",
                field_name="clamp_group", hint="clamp_group Root   # a physical group of the .geo")


def _calculix_rules(pc, vals, rep):
    """calculix: the .frd is ASCII with result blocks; frequency_file exists; ccx has no unit system."""
    def path(k):
        p = Path(str(vals[k]))
        return p if p.is_absolute() else pc.base_dir / p

    if "frequency_file" in vals and not path("frequency_file").exists():
        rep.add("CFG-PATH-001", Severity.ERROR, f"frequency_file does not exist: {path('frequency_file')}",
                field_name="frequency_file", location=pc.location("frequency_file"))
    if "structure_file" in vals and path("structure_file").is_file():
        src = path("structure_file")
        with open(src, "rb") as fh:
            head = fh.read(1 << 20)
        loc = pc.location("structure_file")
        if b"\x00" in head[:4096] or re.search(rb"(?m)^ {4}2C.{18}[ \d]{12}.{37}[23]\s*$", head):
            rep.add("CFG-PATH-003", Severity.ERROR, f"{src.name} is a binary .frd; FLAAP reads ASCII .frd.",
                    field_name="structure_file", location=loc,
                    hint="write the .frd in ASCII or convert it with cgx")
        elif len(head) < (1 << 20) and not re.search(rb"(?m)^ +100C", head):
            rep.add("CFG-PATH-002", Severity.ERROR, f"{src.name}: no result blocks (100C ... DISP) found.",
                    field_name="structure_file", location=loc,
                    hint="add  *NODE FILE  /  U  to the *FREQUENCY step")
    if float(vals.get("length_scale", 1.0)) == 1.0 and float(vals.get("mass_scale", 1.0)) == 1.0:
        rep.add("CFG-DEFAULT-001", Severity.NOTICE,
                "calculix: ccx has no unit system - the .frd is taken as SI (kg, m, s) with unit "
                "generalized mass; set length_scale / mass_scale for mm-tonne-s models.",
                field_name="length_scale", hint="length_scale 0.001  and  mass_scale 1000   # mm, tonne, s")


def _radioss_eig_rules(pc, vals, rep):
    """openradioss_eig: optional files exist, mode VTK files found, unit system known."""
    def path(k):
        p = Path(str(vals[k]))
        return p if p.is_absolute() else pc.base_dir / p

    for k in ("mode_files", "bulk_file", "frequency_file"):
        if k in vals and not path(k).exists():
            rep.add("CFG-PATH-001", Severity.ERROR, f"{k} does not exist: {path(k)}",
                    field_name=k, location=pc.location(k))
    k = "mode_files" if "mode_files" in vals else "structure_file"
    src = path(k) if k in vals else None
    if src is not None and src.exists():
        if src.is_dir():
            found = any(src.glob("*.vtk")) or any(src.glob("*.VTK"))
        elif src.suffix.lower() == ".vtk":
            found = True
        else:                                   # structure_file is the listing: VTK files beside it
            found = k == "structure_file" and any(src.parent.glob("*.vtk"))
        if not found:
            rep.add("CFG-PATH-002", Severity.ERROR, f"no Radioss mode VTK files (<run>A001.vtk ...) for {src}",
                    field_name=k, location=pc.location(k),
                    hint="anim_to_vtk <run>A001 > <run>A001.vtk ...; mode_files path/to/vtk_folder")
    if "radioss_units" in vals:
        u = vals["radioss_units"]
        loc = pc.location("radioss_units")
        if not (isinstance(u, list) and len(u) == 3):
            rep.add("CFG-VALUE-005", Severity.ERROR, "radioss_units needs three values: mass, length, time.",
                    field_name="radioss_units", location=loc, hint="radioss_units kg, mm, ms")
        else:
            from ..io.openradioss import UNIT_LENGTH, UNIT_MASS, UNIT_TIME, _unit
            for tab, x, what in zip((UNIT_MASS, UNIT_LENGTH, UNIT_TIME), u, ("mass", "length", "time")):
                try:
                    _unit(tab, str(x), what)
                except ValueError as e:
                    rep.add("CFG-VALUE-001", Severity.ERROR, f"radioss_units: {e}.", field_name="radioss_units",
                            location=loc)
            if re.search(r"\bMg\b", pc.raw_values.get("radioss_units", "")):
                rep.add("CFG-VALUE-001", Severity.WARNING,
                        "radioss_units 'Mg' is read as milligram (config values are lower-cased).",
                        field_name="radioss_units", location=loc, hint="use 't' or 'tonne' for megagram",
                        promote_in_strict=True)
        if "bulk_file" in vals:
            rep.add("CFG-CONFLICT-001", Severity.WARNING,
                    "radioss_units overrides the /BEGIN units of the starter deck (bulk_file).",
                    field_name="radioss_units", location=loc, hint="omit radioss_units to use the deck units",
                    promote_in_strict=True)
    elif "bulk_file" not in vals:
        rep.add("CFG-DEFAULT-001", Severity.NOTICE,
                "openradioss_eig without bulk_file or radioss_units: mode files taken as SI (kg, m, s) "
                "with unit generalized mass (Inorm = 0).", field_name="radioss_units",
                hint="bulk_file path/to/run_0000.rad  or  radioss_units kg, mm, ms")


def _transform_rules(pc, vals, rep):
    if not any(k in vals for k in ("struct_mirror", "struct_rotate", "struct_translate")):
        return
    from ..structures.transform import transform_from_values
    try:
        T = transform_from_values(vals)
    except ValueError as e:
        key = next((k for k in ("struct_rotate", "struct_translate", "struct_mirror") if k in vals),
                   "struct_rotate")
        rep.add("CFG-VALUE-007", Severity.ERROR, str(e), field_name=key, location=pc.location(key))
        return
    if vals.get("structure", "goland") == "demo_aircraft":
        rep.add("CFG-STRUCT-002", Severity.ERROR,
                "struct_mirror / struct_rotate / struct_translate do not apply to demo_aircraft.",
                field_name="structure", location=pc.location("structure"))
    custom = any(k.startswith("surface") and k[7:].isdigit() for k in vals)
    if not custom and vals.get("structure", "goland") in ("goland", "goland_ctrl", "agard445", "plate"):
        import numpy as np
        if not np.allclose(T.matrix[:, 0], [1.0, 0.0, 0.0]):
            rep.add("CFG-STRUCT-003", Severity.ERROR,
                    f"'{T.describe()}' turns the chord direction: the preset aero surface cannot follow.",
                    field_name="struct_rotate", location=pc.location("struct_rotate"),
                    hint="give the aero surface explicitly (surface1 ...)")
    if custom:
        rep.add("CFG-STRUCT-004", Severity.NOTICE,
                "The structure is re-oriented; the surfaceN keys are taken as written (aero frame: x downstream, "
                "y right, z up).", field_name="struct_rotate")


def _symmetry_rules(pc, vals, rep):
    """Half-model symmetry needs the planform at y >= 0."""
    sym = vals.get("symmetry", "symmetric")
    if sym not in ("symmetric", "antisymmetric"):
        return
    for k, val in vals.items():
        if INDEXED_PATTERNS["surface"][0].match(k) and isinstance(val, list) and len(val) == 11:
            try:
                ys = (float(val[2]), float(val[6]))
            except (TypeError, ValueError):
                continue
            if min(ys) < -1e-9:
                rep.add("CFG-SYM-001", Severity.WARNING,
                        f"symmetry {sym} images the surface about y = 0, but {k} reaches y = {min(ys):g}: "
                        "the image would overlap the model.",
                        field_name=k, location=pc.location(k),
                        hint="model the right half (y >= 0) or use symmetry none for a full model",
                        promote_in_strict=True)


def _indexed_rules(pc, vals, rep):
    surfaces = []
    for key, val in vals.items():
        for kind, (rx, fmt) in INDEXED_PATTERNS.items():
            if rx.match(key):
                n = len(fmt.split(","))
                lst = val if isinstance(val, list) else [val]
                if kind == "surface":
                    surfaces.append(str(lst[0]))
                if len(lst) != n:
                    rep.add("CFG-VALUE-005", Severity.ERROR, f"{key} needs {n} values: {fmt}.",
                            field_name=key, location=pc.location(key))
    for key, val in vals.items():
        if INDEXED_PATTERNS["control"][0].match(key) and isinstance(val, list) and len(val) == 5 and surfaces:
            if str(val[1]) not in surfaces:
                rep.add("CFG-VALUE-006", Severity.ERROR,
                        f"{key} refers to surface {val[1]!r}, defined surfaces: {', '.join(surfaces)}.",
                        field_name=key, location=pc.location(key))


def _flight_rules(pc, vals, rep):
    if "density" in vals and "altitude" in vals:
        rep.add("CFG-CONFLICT-001", Severity.WARNING, "Both density and altitude given; density is used.",
                field_name="altitude", location=pc.location("altitude"), promote_in_strict=True)


def _tsd_rules(pc, vals, rep):
    if vals.get("aero_method") == "tsd" or vals.get("model") == "tsd":
        m = vals.get("mach", 0.0)
        if not 0.5 <= m < 1.0:
            rep.add("CFG-RANGE-001", Severity.WARNING,
                    f"TSD is meant for transonic flow; mach {m} is outside 0.5 <= M < 1.",
                    field_name="mach", location=pc.location("mach"), promote_in_strict=True)
        nsurf = sum(1 for k in vals if k.startswith("surface") and k[7:].isdigit())
        if nsurf > 1:
            rep.add("CFG-TSD-001", Severity.ERROR, "The TSD solver handles a single lifting surface.",
                    field_name="surface2", location=pc.location("surface2"))
        from .builder import tsd_camber, tsd_thickness
        for key, fn in (("tsd_thickness", tsd_thickness), ("tsd_camber", tsd_camber)):
            if key in vals:
                try:
                    fn(vals[key])
                except ValueError as e:
                    rep.add("CFG-VALUE-001", Severity.ERROR, str(e), field_name=key, location=pc.location(key))
        for key in ("tsd_alpha", "alpha"):
            a = vals.get(key)
            if isinstance(a, (int, float)) and abs(a) > 5.0:
                rep.add("CFG-RANGE-006", Severity.WARNING,
                        f"{key} {a} deg: small-disturbance theory loses accuracy beyond a few degrees "
                        "(leading-edge suction peak, separation not modelled).",
                        field_name=key, location=pc.location(key), promote_in_strict=True)
        if "alphas" in vals and vals.get("model") == "tsd" and vals.get("mode", "steady") != "steady":
            rep.add("CFG-TSD-002", Severity.WARNING, "alphas (incidence sweep) is only used by mode steady; "
                    "the GAFs use tsd_alpha.", field_name="alphas", location=pc.location("alphas"),
                    promote_in_strict=True)
    am = vals.get("aero_method", "dlm")
    m = vals.get("mach", 0.0)
    if am in ("machbox", "piston"):
        if m <= 1.0:
            rep.add("CFG-RANGE-003", Severity.ERROR, f"aero_method {am} is supersonic: mach must be > 1.",
                    field_name="mach", location=pc.location("mach"), hint="Use aero_method dlm for M < 1.")
        elif am == "piston" and m < 1.6:
            rep.add("CFG-RANGE-004", Severity.WARNING,
                    f"Piston theory is a high-Mach approximation; at mach {m} use aero_method machbox.",
                    field_name="mach", location=pc.location("mach"), promote_in_strict=True)
    elif m >= 1.0:
        rep.add("CFG-RANGE-002", Severity.ERROR, "DLM/TSD/UVLM are subsonic/transonic only (mach < 1).",
                field_name="mach", location=pc.location("mach"),
                hint="Use aero_method machbox (or piston) for supersonic flow.")
    if (am == "uvlm" or vals.get("model") == "uvlm") and m > 0.3:
        rep.add("CFG-RANGE-005", Severity.WARNING,
                f"The UVLM is incompressible; mach {m} is ignored by the unsteady solver "
                "(the steady VLM applies Prandtl-Glauert).",
                field_name="mach", location=pc.location("mach"), promote_in_strict=True)


def _uvlm_rules(pc, vals, rep):
    if vals.get("mode") == "harmonic" and not (vals.get("pitch_amplitude_deg") or vals.get("plunge_amplitude")):
        rep.add("CFG-UVLM-001", Severity.ERROR, "mode harmonic needs pitch_amplitude_deg and/or plunge_amplitude.",
                field_name="pitch_amplitude_deg", location=pc.location("mode"))
    if vals.get("mode") == "gust" and vals.get("gust_shape", "1-cos") == "1-cos" and "gust_length" not in vals:
        rep.add("CFG-UVLM-002", Severity.ERROR, "A 1-cos gust needs gust_length (gradient H [m]).",
                field_name="gust_length", location=pc.location("gust_shape"))
    if vals.get("free_wake") and vals.get("mode") in ("steady", "flutter"):
        rep.add("CFG-UVLM-003", Severity.WARNING, f"free_wake is only used by the time-marching modes, "
                f"not by mode {vals.get('mode')}.", field_name="free_wake", location=pc.location("free_wake"),
                promote_in_strict=True)


def _choice(key, choices):
    def v(pc, vals, rep):
        if key in vals and str(vals[key]) not in choices:
            rep.add("CFG-VALUE-001", Severity.ERROR,
                    f"{key} {vals[key]!r} is not one of: {', '.join(choices)}.",
                    field_name=key, location=pc.location(key))
    return v


def _needs_demo(pc, vals, rep):
    m, md = vals.get("model"), vals.get("mode")
    if (m == "maneuver" or (m == "trim" and md == "static")) and vals.get("structure") != "demo_aircraft":
        rep.add("CFG-ROUTE-003", Severity.ERROR, f"{m}/{md} needs a free-flying model: structure demo_aircraft.",
                field_name="structure", location=pc.location("structure"),
                hint="Use model trim / mode divergence for clamped wings.")


COMMON_VALIDATORS = (
    _needs_demo,
    _choice("symmetry", ("none", "symmetric", "antisymmetric")),
    _choice("struct_mirror", ("none", "x", "y", "z")),
    _transform_rules,
    _choice("spline", ("auto", "ips", "tps", "rbf", "beam", "legacy")),
    _choice("aero_method", ("dlm", "tsd", "uvlm", "machbox", "piston")),
    _choice("gust_shape", ("sharp", "1-cos")),
    _choice("material", ("aluminum", "steel", "titanium", "cfrp_quasi", "custom")),
    _choice("fe_type", ("shell", "solid")),
    _choice("element", ("tet10", "tet4")),
    _choice("mesh_order", ("1", "2")),
    _choice("mesh_dim", ("2", "3")),
    _choice("nonlinearity", ("freeplay", "bilinear", "preloaded")),
    _positive(("nmodes", "kmax", "nk", "bref", "density", "velocity", "nlag", "thickness", "young",
               "material_density", "mesh_size", "length_scale", "mass_scale", "t_end", "dt", "machbox_nx",
               "wake_chords", "reduced_frequency", "ncycles", "gust_length", "wake_length", "chord_steps")),
    _range3(("velocities", "gust_gradients", "alphas")),
    _structure_rules, _indexed_rules, _symmetry_rules, _flight_rules, _tsd_rules,
)

ROUTES: dict[str, RouteSpec] = {}


def _route(name, description, modes, groups, defaults=None, required=None, validators=(), base=None,
           structures=STRUCTURES):
    base = BASE_DEFAULTS if base is None else base
    ROUTES[name] = RouteSpec(name, description, modes, groups, {**base, **(defaults or {})},
                             required or {}, (_choice("structure", structures),) + COMMON_VALIDATORS
                             + tuple(validators))


_route("modal", "Structural modal model (presets, plate FE, Nastran/ANSYS/OpenRadioss import, STL/STEP -> FE).",
       {"analysis": "frequencies, generalized masses and mode shapes"}, STRUCT_AERO)
_route("aero", "Unsteady aerodynamics: DLM (or TSD) generalized aerodynamic forces and pressures.",
       {"gaf": "GAF table Q(k), gust column and dCp per mode",
        "steady": "steady (k = 0) pressures per mode / unit incidence"}, STRUCT_AERO)
_route("flutter", "Flutter: PK, K (V-g) and p (state-space) methods with automatic identification.",
       {"pk": "PK method", "k": "K method (V-g)", "p": "p method on a Roger RFA state space",
        "all": "run the methods listed in 'methods'"}, STRUCT_AERO + FLIGHT,
       {"velocities": [10.0, 300.0, 146], "methods": ["pk", "k", "p"]})
_route("gust", "Gust loads on the Roger-RFA state space (CS-25 1-cos, von Karman turbulence).",
       {"discrete": "1-cos gust envelope over gradient lengths", "continuous": "von Karman PSD, A-bar, N0"},
       STRUCT_AERO + FLIGHT + (GUST_KEYS,), {"fg": 1.0, "turbulence_length": 762.0},
       {"discrete": ("velocity",), "continuous": ("velocity",)})
_route("trim", "Static aeroelastic trim, divergence and control effectiveness (demo_aircraft or wing).",
       {"static": "1-g (n_z) trim of the free-flying demo aircraft", "divergence": "divergence dynamic pressure"},
       STRUCT_AERO + FLIGHT + (TRIM_KEYS,), {"nz": 1.0, "control": "elevator"},
       {"static": ("velocity",)})
_route("maneuver", "CS-25.331(c) pitch manoeuvres in the time domain (demo_aircraft).",
       {"checked": "checked pitch manoeuvre scaled to nz_limit", "abrupt": "abrupt elevator step"},
       STRUCT_AERO + FLIGHT + (TRIM_KEYS,), {"nz_limit": 2.5, "control": "elevator", "structure": "demo_aircraft"},
       {"checked": ("velocity",), "abrupt": ("velocity",)})
_route("tsd", "Transonic small-disturbance solver (nonlinear steady flow with thickness, camber and incidence; "
       "harmonic GAFs about it).",
       {"steady": "nonlinear steady mean flow: upper / lower Cp, CL, shock positions (incidence sweep with alphas)",
        "gaf": "time-linearised TSD GAF table",
        "flutter": "TSD GAFs + PK flutter"}, STRUCT_AERO + FLIGHT + (TSD_KEYS,),
       {"tsd_thickness": "naca65a004", "mach": 0.9, "aero_method": "tsd", "velocities": [10.0, 300.0, 146]})
_route("uvlm", "Vortex-lattice aerodynamics: steady VLM and nonlinear time-domain UVLM (lift, drag, moment).",
       {"steady": "VLM incidence sweep: CL, CDi (Trefftz plane), Cm, span load, aerodynamic centre",
        "impulsive": "UVLM impulsive start at 'alpha': CL(t), CD(t), Cm(t) (Wagner-type build-up)",
        "harmonic": "UVLM harmonic pitch / plunge: CL, CD (thrust), Cm histories and first harmonic",
        "gust": "UVLM response to a sharp-edged or 1-cos gust (Kuessner-type)",
        "flutter": "linearised UVLM: GAFs + PK and discrete-time state-space eigenvalue sweep"},
       STRUCT_AERO + FLIGHT + (UVLM_KEYS,),
       {"alpha": 0.0, "alphas": [-4.0, 10.0, 8], "pitch_amplitude_deg": 0.0, "plunge_amplitude": 0.0,
        "phase_deg": 90.0, "ncycles": 4.0, "gust_shape": "1-cos", "free_wake": False, "wake_chords": 30.0,
        "aero_method": "uvlm", "velocities": [60.0, 260.0, 41], "chord_steps": 1},
       {"impulsive": ("velocity", "alpha"), "harmonic": ("velocity", "reduced_frequency"),
        "gust": ("velocity", "gust_velocity")},
       validators=(_uvlm_rules,))
NL_BASE = {"plots": True, "verbose": True, "structure": "section", "nmodes": 4}
_route("freeplay", "Control-surface freeplay limit-cycle oscillation (piecewise-linear hinge, exact switching).",
       {"describing_function": "LCO amplitude/frequency vs speed (harmonic balance) + time-marching check",
        "time": "time marching at 'velocity' (time history, phase plane, LCO detection)",
        "bifurcation": "time-marching continuation up and down the 'velocities' sweep (hysteresis)"},
       (NL_PLANT_KEYS, NONLINEAR_KEYS), {"nonlinearity": "freeplay", "damping": 0.01},
       {"time": ("velocity",)}, base=NL_BASE, structures=NL_PLANTS)
_route("ase", "Aeroservoelasticity: actuator + sensors + static output feedback on the RFA state space.",
       {"flutter_suppression": "open- vs closed-loop p-method flutter (V-g, root locus)",
        "margins": "loop Bode/Nyquist, gain/phase margins at 'velocity' and versus speed",
        "gla": "1-cos gust load alleviation (goland_flap: root bending, flap demand)"},
       (NL_PLANT_KEYS, ASE_KEYS, {"gust_gradients": GUST_KEYS["gust_gradients"], "t_end": NONLINEAR_KEYS["t_end"],
                                  "dt": NONLINEAR_KEYS["dt"]}),
       {"actuator_damping": 0.7, "damping": 0.02}, base=NL_BASE, structures=NL_PLANTS)


# ----------------------------------------------------------------- resolve
def resolve(pc: ParsedConfig, report: DiagnosticReport | None = None) -> tuple[dict, DiagnosticReport]:
    """Validate a parsed config against its route; returns (values with defaults, report)."""
    rep = report or DiagnosticReport()
    rep.extend(pc.diagnostics.diagnostics)
    vals = dict(pc.values)
    model = vals.get("model")
    if model is None:
        rep.add("CFG-REQ-001", Severity.ERROR, "Missing required key 'model'.", field_name="model",
                hint="model flutter   # one of " + ", ".join(ROUTES))
        return vals, rep
    if not isinstance(model, str) or model not in ROUTES:
        rep.add("CFG-ROUTE-001", Severity.ERROR, f"Unknown model {model!r}.", field_name="model",
                location=pc.location("model"), hint="Known: " + ", ".join(ROUTES))
        return vals, rep
    spec = ROUTES[model]
    mode = vals.get("mode")
    if mode is None:
        mode = next(iter(spec.modes))
        vals["mode"] = mode
        rep.add("CFG-DEFAULT-001", Severity.NOTICE, f"mode not given; using default {mode!r}.", field_name="mode")
    if mode not in spec.modes:
        rep.add("CFG-ROUTE-002", Severity.ERROR, f"model {model} has no mode {mode!r}.", field_name="mode",
                location=pc.location("mode"), hint="Modes: " + ", ".join(spec.modes))
        return vals, rep
    known = spec.known_keys
    for k in pc.values:
        if k in known or any(rx.match(k) for rx, _ in INDEXED_PATTERNS.values()):
            continue
        import difflib
        sug = difflib.get_close_matches(k, list(known), n=1)
        rep.add("CFG-KEY-001", Severity.WARNING, f"Unknown key {k!r} for model {model}; it is ignored.",
                field_name=k, location=pc.location(k), hint=f"Did you mean {sug[0]!r}?" if sug else None,
                promote_in_strict=True)
    for k in spec.required_by_mode.get(mode, ()):
        if k not in vals:
            rep.add("CFG-REQ-001", Severity.ERROR, f"Missing required key {k!r} for {model}/{mode}.",
                    field_name=k, hint=f"{k} <value>   # {known.get(k, ('', ''))[0]}")
    for k, d in spec.defaults.items():
        if k not in vals and k in known:
            vals[k] = d
            rep.add("CFG-DEFAULT-001", Severity.NOTICE, f"{k} not given; default {d!r} is used.", field_name=k)
    # type checks
    for k, v in pc.values.items():
        typ = known.get(k, ("", "any"))[1]
        bad = ((typ == "int" and not (isinstance(v, int) and not isinstance(v, bool)))
               or (typ == "float" and not (isinstance(v, (int, float)) and not isinstance(v, bool)))
               or (typ == "bool" and not isinstance(v, bool))
               or (typ == "vec" and not (isinstance(v, list) and len(v) in (2, 3))))
        if bad:
            rep.add("CFG-TYPE-001", Severity.ERROR, f"{k} should be of type {typ} (got {v!r}).",
                    field_name=k, location=pc.location(k))
    for vfun in spec.validators:
        vfun(pc, vals, rep)
    if "output_dir" not in vals:
        stem = pc.source_path.stem if pc.source_path else "run"
        vals["output_dir"] = str(Path("output") / stem)
        rep.add("CFG-DEFAULT-001", Severity.NOTICE, f"output_dir not given; using {vals['output_dir']}.",
                field_name="output_dir")
    return vals, rep


def describe(model: str) -> str:
    spec = ROUTES[model]
    lines = [f"model {model}: {spec.description}", "", "modes:"]
    lines += [f"  {m:<22s}{d}" for m, d in spec.modes.items()]
    req = {m: r for m, r in spec.required_by_mode.items() if r}
    if req:
        lines += ["", "required keys:"] + [f"  mode {m}: {', '.join(r)}" for m, r in req.items()]
    lines += ["", "keys:"]
    for k, (d, t) in spec.known_keys.items():
        dv = spec.defaults.get(k)
        lines.append(f"  {k:<26s}{t:<6s}{d}" + (f"  [default {dv}]" if dv is not None else ""))
    lines += ["", "indexed keys:"] + [f"  {n}N  {fmt}" for n, (_, fmt) in INDEXED_PATTERNS.items()]
    return "\n".join(lines)
