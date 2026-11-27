import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


@pytest.fixture(scope="session")
def root_dir():
    return ROOT


@pytest.fixture(scope="session")
def legacy_dir():
    """Original Fortran-FLAAP files (decks, FLAAP.exe results, 20.f06)."""
    return os.path.join(ROOT, "legacy_fortran")


@pytest.fixture(scope="session")
def goland():
    """Goland wing: 4-mode beam model, 6x12 DLM table (M = 0)."""
    from flaap.aero.gaf import AeroModel, AeroModes
    from flaap.geometry import build_mesh
    from flaap.models.goland import CHORD, goland_modal, goland_surface
    from flaap.structures.spline import beam_spline
    mm = goland_modal(4, damping=0.0)
    mesh = build_mesh([goland_surface(6, 12)], "symmetric")
    modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors())
    am = AeroModel(mesh, modes, 0.0, CHORD / 2)
    tab = am.table(np.r_[0.0, np.geomspace(0.01, 2.0, 20)], gust=True, keep_cp=True)
    return {"modal": mm, "mesh": mesh, "aero": am, "table": tab}


@pytest.fixture(scope="session")
def demo():
    """Demo flexible aircraft (rigid plunge/pitch + 4 wing modes + elevator)."""
    from flaap.aero.gaf import AeroModel
    from flaap.models.demo_aircraft import demo_aircraft
    ac, mesh, modes = demo_aircraft(4)
    tab = AeroModel(mesh, modes, 0.0, 0.9144).table(np.r_[0, np.geomspace(0.01, 2.0, 18)])
    return {"ac": ac, "mesh": mesh, "modes": modes, "table": tab}
