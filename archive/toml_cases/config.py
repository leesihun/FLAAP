"""TOML case files: build the model, run the requested analyses, write results.

Example (``cases/goland.toml``)::

    title = "Goland wing"

    [model]
    type = "goland"        # goland | agard445 | legacy | nastran
    nmodes = 4
    damping = 0.0          # modal damping ratio zeta

    [aero]
    mach = 0.0
    nchord = 6
    nspan = 12
    kmax = 2.0
    nk = 20

    [flutter]
    method = ["pk", "k", "p"]
    density = 1.225        # or: altitude = 0.0
    velocities = [60.0, 260.0, 101]   # start, stop, count

    [gust]
    velocity = 100.0
    altitude = 0.0
    gradients = [9.0, 107.0, 8]
    nlag = 4

    [output]
    directory = "results/goland"
    plots = true

``model.type = "nastran"`` reads ``model.f06`` (+ optional ``model.bulk``)
and needs one or more ``[[aero.surface]]`` tables (``name, root_le,
root_chord, tip_le, tip_chord, nchord, nspan``) plus ``aero.bref``.
``model.type = "legacy"`` runs a Fortran-FLAAP input directory
(``model.path``).
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

try:                                    # Python >= 3.11
    import tomllib as _toml
except ModuleNotFoundError:             # pragma: no cover
    try:
        import tomli as _toml
    except ModuleNotFoundError:
        _toml = None


def load(path: str) -> dict:
    if _toml is None:
        raise RuntimeError("TOML support needs Python >= 3.11 or the 'tomli' package")
    with open(path, "rb") as f:
        cfg = _toml.load(f)
    cfg.setdefault("_dir", os.path.dirname(os.path.abspath(path)))
    return cfg


def _range(v, default):
    if v is None:
        return np.asarray(default, float)
    if isinstance(v, (list, tuple)) and len(v) == 3 and isinstance(v[2], int):
        return np.linspace(float(v[0]), float(v[1]), int(v[2]))
    return np.asarray(v, float)


def _density(sec: dict) -> float:
    from .atmosphere import isa
    if "density" in sec:
        return float(sec["density"])
    return isa(float(sec.get("altitude", 0.0))).density


@dataclass
class CaseModel:
    modal: object
    aero: object                 # AeroModel
    title: str = ""
    extra: dict = field(default_factory=dict)


def build(cfg: dict) -> CaseModel:
    """Structural + aerodynamic model of a case dictionary."""
    from .aero.gaf import AeroModel, AeroModes
    from .geometry import Surface, build_mesh
    from .structures.spline import beam_spline, surface_spline
    m = cfg.get("model", {})
    a = cfg.get("aero", {})
    typ = m.get("type", "goland").lower()
    nmodes = int(m.get("nmodes", 4))
    zeta = float(m.get("damping", 0.0))
    mach = float(a.get("mach", 0.0))
    base = cfg.get("_dir", ".")
    if typ == "goland":
        from .models.goland import CHORD, goland_modal, goland_surface
        modal = goland_modal(nmodes, damping=zeta)
        mesh = build_mesh([goland_surface(int(a.get("nchord", 6)), int(a.get("nspan", 12)))], "symmetric")
        sp = beam_spline(mesh, modal.nodes)
        bref = CHORD / 2
    elif typ in ("agard445", "agard"):
        from .models.agard445 import B_REF, EXPERIMENT_FREQ_HZ, agard445_modal, agard445_surface
        modal = agard445_modal(nmodes)
        if m.get("experimental_frequencies", False):
            modal.omega[:4] = 2 * np.pi * EXPERIMENT_FREQ_HZ[:min(4, nmodes)]
        modal.damping[:] = zeta
        mesh = build_mesh([agard445_surface(int(a.get("nchord", 8)), int(a.get("nspan", 10)))], "symmetric")
        sp = surface_spline(mesh, modal.nodes, method="ips")
        bref = B_REF
    elif typ == "legacy":
        from .io.legacy import legacy_aero, read_case
        case = read_case(os.path.join(base, m.get("path", ".")))
        modal = case.modal(float(m.get("length_scale", 1.0)), m.get("nmodes"))
        modal.damping[:] = zeta
        am = legacy_aero(case, modal, a.get("spline", "legacy"), a.get("mach"))
        return CaseModel(modal, am, cfg.get("title", "legacy"), {"case": case})
    elif typ == "nastran":
        from .io.nastran import read_f06_modal
        bulk = m.get("bulk")
        modal = read_f06_modal(os.path.join(base, m["f06"]),
                               None if bulk is None else os.path.join(base, bulk),
                               nmodes=m.get("nmodes"), length_scale=float(m.get("length_scale", 1.0)),
                               mass_scale=float(m.get("mass_scale", 1.0)),
                               transform=m.get("transform"), origin=m.get("origin", (0, 0, 0)))
        modal.damping[:] = zeta
        surfs = [Surface(s["name"], tuple(s["root_le"]), float(s["root_chord"]), tuple(s["tip_le"]),
                         float(s["tip_chord"]), nchord=int(s.get("nchord", 6)), nspan=int(s.get("nspan", 10)))
                 for s in a["surface"]]
        mesh = build_mesh(surfs, a.get("symmetry", "symmetric"))
        sp = surface_spline(mesh, modal.nodes, method=a.get("spline", "ips"))
        bref = float(a["bref"])
    else:
        raise ValueError(f"unknown model type {typ!r}")
    modes = AeroModes.from_spline(sp, modal.dof_vectors(), modal.labels or
                                  [f"mode{i + 1}" for i in range(modal.nmodes)])
    return CaseModel(modal, AeroModel(mesh, modes, mach, bref), cfg.get("title", typ))


def gaf_table(cfg: dict, model: CaseModel, verbose=False):
    a = cfg.get("aero", {})
    if "k" in a:
        kv = np.asarray(a["k"], float)
    else:
        kmax, nk = float(a.get("kmax", 2.0)), int(a.get("nk", 20))
        kv = np.r_[0.0, np.geomspace(min(0.01, kmax / 100), kmax, nk - 1)]
    return model.aero.table(kv, gust=True, keep_cp=False, verbose=verbose)


def run(cfg: dict, outdir: str | None = None, verbose: bool = True) -> dict:
    """Run all analyses requested by ``cfg``; returns a result dictionary."""
    from .flutter import identify, k_method, pk_method
    out = outdir or os.path.join(cfg.get("_dir", "."), cfg.get("output", {}).get("directory", "results"))
    os.makedirs(out, exist_ok=True)
    plots = cfg.get("output", {}).get("plots", True)
    log = print if verbose else (lambda *a, **k: None)
    model = build(cfg)
    mm = model.modal
    log(f"== {model.title}")
    log(mm.table())
    tab = gaf_table(cfg, model, verbose=False)
    tab.save(os.path.join(out, "gaf_table.npz"))
    tab.write_text(os.path.join(out, "gaf_table.txt"))
    log(f"GAF table: {tab.k.size} reduced frequencies, M = {tab.mach}, {model.aero.mesh.n} boxes")
    results = {"model": model, "table": tab}
    nm = mm.nmodes
    tab_s = tab.subset(range(nm))
    # ------------------------------------------------------------ flutter
    if "flutter" in cfg:
        f = cfg["flutter"]
        rho = _density(f)
        V = _range(f.get("velocities"), np.linspace(10, 300, 146))
        methods = f.get("method", "pk")
        methods = [methods] if isinstance(methods, str) else list(methods)
        results["flutter"] = {}
        for meth in methods:
            if meth == "pk":
                r = pk_method(mm.M, mm.K, tab_s, rho, V, C=mm.C())
            elif meth == "k":
                r = k_method(mm.M, mm.K, tab_s, rho)
            elif meth == "p":
                from .aero.rfa import fit_roger
                from .timedomain import modal_system
                rfa = fit_roger(tab_s, int(f.get("nlag", 4)))
                r = modal_system(mm, rfa).flutter(V, rho)
            else:
                raise ValueError(f"unknown flutter method {meth!r}")
            identify(r)
            results["flutter"][meth] = r
            log(r.summary())
            d = os.path.join(out, f"flutter_{meth}")
            os.makedirs(d, exist_ok=True)
            r.write_legacy(d)
            if plots:
                _save_flutter_plot(r, os.path.join(out, f"flutter_{meth}.png"))
    # ------------------------------------------------------------ gust
    if "gust" in cfg:
        from . import gust as G
        from .aero.rfa import fit_roger
        from .timedomain import modal_system
        g = cfg["gust"]
        V = float(g["velocity"])
        alt = float(g.get("altitude", 0.0))
        rho = _density(g)
        rfa = fit_roger(tab_s, int(g.get("nlag", 4)))
        sysm = modal_system(mm, rfa)
        Hs = _range(g.get("gradients"), G.gust_gradients())
        fg = float(g.get("fg", 1.0))
        env = G.gust_envelope(sysm, V, rho, Hs, lambda H: G.design_gust_velocity(H, alt, fg))
        results["gust"] = env
        lines = ["label, H_critical [m], max, min"]
        for j, lab in enumerate(env["labels"]):
            lines.append(f"{lab}, {env['critical_H'][lab]:.2f}, {env['max'][:, j].max():.6e}, "
                         f"{env['min'][:, j].min():.6e}")
        with open(os.path.join(out, "gust_envelope.csv"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        log(f"gust envelope (V = {V} m/s, rho = {rho:.4f}): " + "; ".join(lines[1:4]))
        if plots:
            ss = sysm.state_space(V, rho)
            Hc = env["critical_H"][env["labels"][0]]
            gr = G.discrete_gust(ss, V, Hc, G.design_gust_velocity(Hc, alt, fg))
            _save_gust_plot(gr, os.path.join(out, "gust_time.png"))
    with open(os.path.join(out, "summary.txt"), "w") as fh:
        fh.write(f"{model.title}\n{mm.table()}\n")
        for meth, r in results.get("flutter", {}).items():
            fh.write(r.summary() + "\n")
    log(f"results written to {out}")
    return results


def _save_flutter_plot(r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ax = r.plot()
    ax[0].figure.savefig(path, dpi=110)
    plt.close(ax[0].figure)


def _save_gust_plot(gr, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .plotting import plot_gust
    ax = plot_gust(gr)
    ax[0].figure.savefig(path, dpi=110)
    plt.close(ax[0].figure)
