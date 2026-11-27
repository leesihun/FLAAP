"""Simple flexible "aircraft" for trim / manoeuvre / gust demonstrations.

Half model (symmetric): the Goland wing (elastic beam, clamped to the
fuselage) + a rigid horizontal tail with elevator + a lumped fuselage
mass.  Rigid plunge/pitch modes about the reference point are combined
with the fixed-root wing modes (:func:`flaap.maneuver.flexible_aircraft`).
All values are illustrative, not a real aircraft.
"""
from __future__ import annotations

from ..geometry import ControlSurface, Surface, build_mesh
from ..maneuver import flexible_aircraft
from ..structures.spline import beam_spline
from ..aero.gaf import AeroModes, rigid_aero_modes, control_aero_mode
from .goland import CHORD, SPAN, goland_modal

X_REF = 0.65                 # reference point / c.g. [m]
TAIL_LE = 6.0
TAIL_CHORD = 1.0
TAIL_SPAN = 2.0
BODY = {"m": 400.0, "x": 0.60, "Iyy": 2500.0}     # half fuselage


def demo_surfaces(nchord=6, nspan=12, tail_nchord=4, tail_nspan=4):
    wing = Surface("wing", (0, 0, 0), CHORD, (0, SPAN, 0), CHORD, nchord=nchord, nspan=nspan)
    tail = Surface("htail", (TAIL_LE, 0, 0), TAIL_CHORD, (TAIL_LE + 0.3, TAIL_SPAN, 0), 0.7 * TAIL_CHORD,
                   nchord=tail_nchord, nspan=tail_nspan,
                   control_surfaces=[ControlSurface("elevator", 0.7, (0.0, 1.0))])
    return [wing, tail]


def demo_aircraft(nmodes=4, damping=0.02, **mesh_kw):
    """Returns ``(aircraft, mesh, aero_modes)`` with variables
    ``[plunge, pitch, wing modes..., elevator]``."""
    elastic = goland_modal(nmodes, damping=damping)
    ac = flexible_aircraft(elastic, [X_REF, 0, 0], ("plunge", "pitch"), body_mass=BODY)
    mesh = build_mesh(demo_surfaces(**mesh_kw), "symmetric")
    wing = mesh.surface_id == 0
    sp = beam_spline(mesh, elastic.nodes, box_mask=wing)
    hc, hx, hf = sp.modal(elastic.dof_vectors())
    el = AeroModes(hc, hx, hf, ac.names[2:])
    modes = rigid_aero_modes(mesh, X_REF, ("plunge", "pitch")) + el + control_aero_mode(mesh, "elevator")
    return ac, mesh, modes
