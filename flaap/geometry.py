"""Lifting-surface geometry and aerodynamic box (panel) meshes.

Conventions (used throughout FLAAP-Python)
-----------------------------------------
* aerodynamic axes: ``x`` downstream, ``y`` to the right (span), ``z`` up
* SI units
* every lifting surface is a CAERO1-like trapezoid split into
  ``nchord x nspan`` boxes; each box carries

  - a doublet (bound vortex) line at 1/4 box chord,
  - a collocation (downwash) point at 3/4 box chord, mid-span,
  - a pressure / force point at the middle of the doublet line,
  - a normal ``n = (0, -sin(gamma), cos(gamma))`` (gamma = dihedral).

* the symmetry plane is ``y = 0`` (only the right half is modelled when
  ``symmetry`` is ``'symmetric'`` or ``'antisymmetric'``).

The Fortran routine ``aero_paneling_VLM.f90`` used a span-along-``X`` /
chord-along-``Y`` legacy frame with the trailing edge at ``Y = 0``.  Use
:func:`surface_from_legacy` to convert ``LEE_input.txt`` geometry; it
produces exactly the same box corners (``X1..X4``/``Y1..Y4``), 1/4-chord
pressure points (``CP``) and 3/4-chord vortex points (``lvpx``) as the
original, but expressed in the aerodynamic frame.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


def _spacing(n: int, kind: str | Sequence[float] = "linear") -> np.ndarray:
    """Return ``n+1`` division fractions in [0, 1]."""
    if not isinstance(kind, str):
        s = np.asarray(kind, dtype=float)
        if s.size != n + 1 or abs(s[0]) > 1e-12 or abs(s[-1] - 1) > 1e-12:
            raise ValueError("explicit spacing must have n+1 values from 0 to 1")
        return s
    i = np.arange(n + 1) / n
    if kind == "linear":
        return i
    if kind == "cosine":
        return 0.5 * (1 - np.cos(np.pi * i))
    if kind == "sine":          # clustered at the outboard / trailing end
        return np.sin(0.5 * np.pi * i)
    raise ValueError(f"unknown spacing '{kind}'")


@dataclass
class ControlSurface:
    """Trailing-edge control surface on a lifting surface.

    ``hinge`` is the hinge-line chord fraction (0 = LE, 1 = TE), ``eta`` the
    spanwise fraction range covered.  Boxes whose centre lies aft of the
    hinge line inside the span range belong to the surface.
    """

    name: str
    hinge: float = 0.75
    eta: tuple[float, float] = (0.0, 1.0)
    antisymmetric: bool = False   # informative: ailerons on a full model


@dataclass
class Surface:
    """Trapezoidal lifting surface (Nastran CAERO1 analogue)."""

    name: str
    root_le: Sequence[float]
    root_chord: float
    tip_le: Sequence[float]
    tip_chord: float
    nchord: int = 8
    nspan: int = 10
    chord_spacing: str | Sequence[float] = "linear"
    span_spacing: str | Sequence[float] = "linear"
    control_surfaces: list[ControlSurface] = field(default_factory=list)

    @property
    def span(self) -> float:
        r = np.asarray(self.root_le, float)
        t = np.asarray(self.tip_le, float)
        return float(np.hypot(t[1] - r[1], t[2] - r[2]))

    @property
    def area(self) -> float:
        return 0.5 * (self.root_chord + self.tip_chord) * self.span

    def corner_grid(self) -> np.ndarray:
        """(nspan+1, nchord+1, 3) array of lattice corner points."""
        r = np.asarray(self.root_le, float)
        t = np.asarray(self.tip_le, float)
        eta = _spacing(self.nspan, self.span_spacing)
        xi = _spacing(self.nchord, self.chord_spacing)
        le = r[None, :] + eta[:, None] * (t - r)[None, :]
        ch = self.root_chord + eta * (self.tip_chord - self.root_chord)
        pts = np.repeat(le[:, None, :], self.nchord + 1, axis=1).copy()
        pts[:, :, 0] += ch[:, None] * xi[None, :]
        return pts


@dataclass
class AeroMesh:
    """Flat arrays describing ``N`` aerodynamic boxes.

    ``corners[:, 0..3]`` = inboard-LE, outboard-LE, outboard-TE, inboard-TE.
    """

    corners: np.ndarray              # (N,4,3)
    surface_id: np.ndarray           # (N,)
    strip_id: np.ndarray             # (N,)  global strip index
    chord_index: np.ndarray          # (N,)  0 = leading box
    surfaces: list[Surface]
    symmetry: str = "none"           # 'none' | 'symmetric' | 'antisymmetric'
    control: dict[str, np.ndarray] = field(default_factory=dict)  # name -> (N,) hinge x or nan

    # ---------------------------------------------------------------- derived
    def __post_init__(self):
        c = self.corners
        self.n = c.shape[0]
        self.A = c[:, 0] + 0.25 * (c[:, 3] - c[:, 0])         # doublet line start
        self.B = c[:, 1] + 0.25 * (c[:, 2] - c[:, 1])         # doublet line end
        p34_in = c[:, 0] + 0.75 * (c[:, 3] - c[:, 0])
        p34_out = c[:, 1] + 0.75 * (c[:, 2] - c[:, 1])
        self.colloc = 0.5 * (p34_in + p34_out)                 # 3/4 chord
        self.force_point = 0.5 * (self.A + self.B)             # 1/4 chord
        self.center = c.mean(axis=1)
        dy = self.B[:, 1] - self.A[:, 1]
        dz = self.B[:, 2] - self.A[:, 2]
        self.gamma = np.arctan2(dz, dy)                        # dihedral
        self.e = 0.5 * np.hypot(dy, dz)                        # semi-width
        self.tan_sweep = (self.B[:, 0] - self.A[:, 0]) / (2 * self.e)
        self.chord = 0.5 * ((c[:, 3, 0] - c[:, 0, 0]) + (c[:, 2, 0] - c[:, 1, 0]))
        self.normal = np.stack([np.zeros(self.n), -np.sin(self.gamma),
                                np.cos(self.gamma)], axis=1)
        # area (planform projected on the box plane)
        d1 = c[:, 2] - c[:, 0]
        d2 = c[:, 3] - c[:, 1]
        self.area = 0.5 * np.linalg.norm(np.cross(d1, d2), axis=1)

    # ---------------------------------------------------------------- helpers
    @property
    def total_area(self) -> float:
        return float(self.area.sum())

    def strips(self):
        """Yield (strip index, box indices sorted LE->TE)."""
        for s in np.unique(self.strip_id):
            idx = np.where(self.strip_id == s)[0]
            yield s, idx[np.argsort(self.chord_index[idx])]

    def image(self) -> "AeroMesh":
        """Mirror image about ``y = 0`` with lift-positive orientation."""
        c = self.corners.copy()
        c[:, :, 1] *= -1
        c = c[:, [1, 0, 3, 2]]     # keep 'A' inboard->outboard ordering w.r.t. lift
        m = AeroMesh(c, self.surface_id, self.strip_id, self.chord_index,
                     self.surfaces, "none", dict(self.control))
        return m

    def interior_hinge_mask(self, name: str) -> np.ndarray:
        return np.isfinite(self.control[name])

    def summary(self) -> str:
        return (f"AeroMesh: {self.n} boxes, {len(np.unique(self.strip_id))} strips, "
                f"area={self.total_area:.4g} m^2, symmetry={self.symmetry}")


def build_mesh(surfaces: Sequence[Surface], symmetry: str = "none") -> AeroMesh:
    corners, sid, strip, cidx = [], [], [], []
    controls: dict[str, list[np.ndarray]] = {}
    strip0 = 0
    for k, s in enumerate(surfaces):
        g = s.corner_grid()
        ns, nc = s.nspan, s.nchord
        eta_edges = _spacing(ns, s.span_spacing)
        xi_edges = _spacing(nc, s.chord_spacing)
        hinge_x = {cs.name: np.full(ns * nc, np.nan) for cs in s.control_surfaces}
        b = 0
        for j in range(ns):
            for i in range(nc):
                corners.append(np.stack([g[j, i], g[j + 1, i], g[j + 1, i + 1], g[j, i + 1]]))
                sid.append(k)
                strip.append(strip0 + j)
                cidx.append(i)
                eta_c = 0.5 * (eta_edges[j] + eta_edges[j + 1])
                xi_c = 0.5 * (xi_edges[i] + xi_edges[i + 1])
                for cs in s.control_surfaces:
                    if cs.eta[0] - 1e-9 <= eta_c <= cs.eta[1] + 1e-9 and xi_c > cs.hinge:
                        # hinge x location at box mid-span (local chord)
                        le = 0.5 * (g[j, 0, 0] + g[j + 1, 0, 0])
                        ch = 0.5 * ((g[j, -1, 0] - g[j, 0, 0]) + (g[j + 1, -1, 0] - g[j + 1, 0, 0]))
                        hinge_x[cs.name][b] = le + cs.hinge * ch
                b += 1
        for name, hx in hinge_x.items():
            controls.setdefault(name, [])
            controls[name].append(hx)
        # surfaces without that control: pad
        strip0 += ns
    n_per = [s.nspan * s.nchord for s in surfaces]
    control = {}
    for name in controls:
        arr = []
        it = iter(controls[name])
        for s, n in zip(surfaces, n_per):
            if any(cs.name == name for cs in s.control_surfaces):
                arr.append(next(it))
            else:
                arr.append(np.full(n, np.nan))
        control[name] = np.concatenate(arr)
    return AeroMesh(np.array(corners), np.array(sid), np.array(strip), np.array(cidx),
                    list(surfaces), symmetry, control)


# --------------------------------------------------------------------------
# Legacy (Fortran) geometry
# --------------------------------------------------------------------------
@dataclass
class LegacyGeometry:
    """Planform as given in ``LEE_input.txt`` (span along X, chord along Y,
    trailing edge at Y = 0)."""

    span: float
    chord: float
    ibte: tuple[float, float]
    obte: tuple[float, float]
    oble: tuple[float, float]
    ible: tuple[float, float]

    def to_aero(self, x_leg, y_leg):
        """Legacy (X span, Y chord) -> aero (x, y)."""
        return self.ible[1] - np.asarray(y_leg), np.asarray(x_leg)

    def from_aero(self, x, y):
        return np.asarray(y), self.ible[1] - np.asarray(x)


def surface_from_legacy(geo: LegacyGeometry, n_span: int, n_chord: int,
                        name: str = "wing") -> Surface:
    """Build the aero ``Surface`` equivalent to ``aero_paneling_VLM``."""
    x_ible, y_ible = geo.ible
    x_oble, y_oble = geo.oble
    x_ibte, y_ibte = geo.ibte
    x_obte, y_obte = geo.obte
    root_le = (0.0, x_ible, 0.0)
    tip_le = (y_ible - y_oble, x_oble, 0.0)
    return Surface(name, root_le, abs(y_ible - y_ibte), tip_le, abs(y_oble - y_obte),
                   nchord=n_chord, nspan=n_span)


def legacy_panel_arrays(mesh: AeroMesh, geo: LegacyGeometry) -> dict[str, np.ndarray]:
    """Return the Fortran ``panel`` module arrays (X1..X4, Y1..Y4, CP, lvp).

    Legacy ordering: corner 1 = inboard TE side? The Fortran code used
    ``x1``/``x4`` along the chord measured from the trailing edge, with the
    1/4 point from ``x1`` being the pressure point (``CPX``) i.e. the 3/4-
    chord collocation point in the aero frame.
    """
    c = mesh.corners
    X = {}
    # legacy y (chord, from TE) of inboard/outboard TE and LE corners
    _, y1 = geo.from_aero(c[:, 3, 0], c[:, 3, 1])  # inboard TE
    _, y2 = geo.from_aero(c[:, 2, 0], c[:, 2, 1])  # outboard TE
    _, y3 = geo.from_aero(c[:, 1, 0], c[:, 1, 1])  # outboard LE
    _, y4 = geo.from_aero(c[:, 0, 0], c[:, 0, 1])  # inboard LE
    X["x1"], X["x2"], X["x3"], X["x4"] = y1, y2, y3, y4    # (chordwise, legacy 'x')
    X["y1"] = c[:, 3, 1]
    X["y2"] = c[:, 2, 1]
    X["y3"] = c[:, 1, 1]
    X["y4"] = c[:, 0, 1]
    X["cpx"] = geo.from_aero(mesh.colloc[:, 0], mesh.colloc[:, 1])[1]
    X["cpy"] = mesh.colloc[:, 1]
    X["lvpx"] = geo.from_aero(mesh.A[:, 0], mesh.A[:, 1])[1]
    X["lvpy"] = mesh.A[:, 1]
    X["rvpy"] = mesh.B[:, 1]
    X["area"] = mesh.area
    return X
