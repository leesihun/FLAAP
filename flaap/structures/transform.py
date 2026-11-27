"""Rigid re-orientation of a structural (modal) model: mirror, rotate, translate.

Config keys (applied in this order: mirror -> rotate -> translate)::

    struct_mirror     y                 # reflect about the plane normal to x | y | z (y = left <-> right)
    struct_rotate     z, 90             # axis, degrees (right-hand rule); pairs may be chained: z, 90, x, 180
    struct_translate  0, 1.5, 0         # shift [m] (after length_scale)

Why it exists: an imported FE model (Nastran, ANSYS, STL ...) is rarely modelled in the aerodynamic
frame (x downstream, y right, z up).  The transform brings it there without touching the file.

A modal model is ``nodes`` (Nn,3) + ``shapes`` (Nn,6,nm) = (ux, uy, uz, rx, ry, rz).  Under an orthogonal
map ``Q`` the node moves ``x' = Q x``; translations rotate like vectors, ``u' = Q u``; rotations are
pseudo-vectors, ``r' = det(Q) Q r`` (a reflection flips the tangential components).  Generalised masses
and frequencies are invariant.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field

import numpy as np

AXES = {"x": 0, "y": 1, "z": 2}


@dataclass
class Transform:
    mirror: str | None = None                       # None | 'x' | 'y' | 'z'
    rotations: list[tuple[str, float]] = field(default_factory=list)   # (axis, degrees)
    translate: tuple[float, float, float] = (0.0, 0.0, 0.0)

    @property
    def identity(self) -> bool:
        return self.mirror is None and not self.rotations and not any(self.translate)

    @property
    def matrix(self) -> np.ndarray:
        """Orthogonal 3x3 ``Q`` (mirror first, then the rotations in order)."""
        Q = np.eye(3)
        if self.mirror:
            Q[AXES[self.mirror], AXES[self.mirror]] = -1.0
        for ax, deg in self.rotations:
            Q = rotation_matrix(ax, deg) @ Q
        return Q

    @property
    def det(self) -> float:
        return float(np.sign(np.linalg.det(self.matrix)))

    def describe(self) -> str:
        parts = []
        if self.mirror:
            parts.append(f"mirror {self.mirror}")
        parts += [f"rotate {a} {d:g} deg" for a, d in self.rotations]
        if any(self.translate):
            parts.append("translate " + ", ".join(f"{t:g}" for t in self.translate))
        return "; ".join(parts) or "identity"


def rotation_matrix(axis: str, deg: float) -> np.ndarray:
    a = np.radians(float(deg))
    c, s = np.cos(a), np.sin(a)
    if abs(c) < 1e-12:                                  # exact for multiples of 90 deg
        c = 0.0
    if abs(s) < 1e-12:
        s = 0.0
    i = AXES[axis]
    j, k = (i + 1) % 3, (i + 2) % 3
    R = np.eye(3)
    R[j, j] = c
    R[j, k] = -s
    R[k, j] = s
    R[k, k] = c
    return R


def transform_from_values(v: dict) -> Transform:
    """Build a :class:`Transform` from the flat config values (raises ``ValueError`` when malformed)."""
    mir = v.get("struct_mirror")
    mir = None if mir in (None, False, "none", "", "off") else str(mir).lower()
    if mir is not None and mir not in AXES:
        raise ValueError(f"struct_mirror {mir!r}: use none, x, y or z")
    rot = v.get("struct_rotate")
    rots: list[tuple[str, float]] = []
    if rot not in (None, False, "none"):
        lst = rot if isinstance(rot, list) else [rot]
        if len(lst) % 2:
            raise ValueError("struct_rotate takes axis, degrees pairs (e.g. 'z, 90')")
        for a, d in zip(lst[0::2], lst[1::2]):
            a = str(a).lower()
            if a not in AXES:
                raise ValueError(f"struct_rotate axis {a!r}: use x, y or z")
            try:
                rots.append((a, float(d)))
            except (TypeError, ValueError):
                raise ValueError(f"struct_rotate angle {d!r} is not a number") from None
    tr = v.get("struct_translate")
    t = (0.0, 0.0, 0.0)
    if tr not in (None, False, "none"):
        lst = tr if isinstance(tr, list) else [tr]
        if len(lst) != 3:
            raise ValueError("struct_translate takes 3 values: dx, dy, dz")
        t = tuple(float(x) for x in lst)
    return Transform(mir, rots, t)


# ---------------------------------------------------------------------------------- nodes / modes
def transform_points(P, T: Transform) -> np.ndarray:
    P = np.asarray(P, float)
    return P @ T.matrix.T + np.asarray(T.translate, float)


def transform_modal(mm, T: Transform):
    """Return a transformed copy of a ``ModalModel`` (nodes, shapes); FE connectivity is not carried over.

    ``mm.fe`` is dropped because the stiffness matrices still refer to the original axes; routes that
    recover element loads from the FE model (``modal_load_matrices``) then stop with a clear message.
    """
    if T.identity:
        return mm
    if not hasattr(mm, "shapes") or not hasattr(mm, "nodes"):
        raise ValueError("structure_transform: this structure has no nodal mode shapes to re-orient")
    Q = T.matrix
    out = copy.copy(mm)
    out.nodes = transform_points(mm.nodes, T)
    sh = np.asarray(mm.shapes, float)
    new = np.empty_like(sh)
    new[:, 0:3, :] = np.einsum("ij,njm->nim", Q, sh[:, 0:3, :])
    new[:, 3:6, :] = T.det * np.einsum("ij,njm->nim", Q, sh[:, 3:6, :])
    out.shapes = new
    out.fe = None
    return out


# ---------------------------------------------------------------------------------- aero surfaces
def transform_surface(s, T: Transform):
    """Transform a trapezoidal :class:`~flaap.geometry.Surface` (preset planform).

    A ``Surface`` has its chords along +x, so the map must keep +x (the flow direction) and keep the span
    in the y-z plane; mirrors about y / z and a rotation about x qualify.  A mirror about y re-orders
    root/tip so that the span runs in +y again (control-surface ``eta`` ranges are flipped with it).
    Anything else raises ``ValueError`` -> give explicit ``surfaceN`` keys instead.
    """
    from ..geometry import ControlSurface, Surface
    Q = T.matrix
    if not np.allclose(Q[:, 0], [1.0, 0.0, 0.0]):
        raise ValueError(f"'{T.describe()}' turns the chord direction (x): a preset aero surface cannot follow; "
                         "give the aero surface explicitly (surface1 ...) or convert the preset in Studio")
    r = transform_points(s.root_le, T)
    t = transform_points(s.tip_le, T)
    rc, tc = s.root_chord, s.tip_chord
    ctrls = [ControlSurface(c.name, c.hinge, tuple(c.eta), c.antisymmetric) for c in s.control_surfaces]
    swapped = bool(t[1] < r[1] - 1e-12)
    if swapped:                                          # span now runs towards -y: swap root and tip
        r, t, rc, tc = t, r, tc, rc
        ctrls = [ControlSurface(c.name, c.hinge, (1.0 - c.eta[1], 1.0 - c.eta[0]), c.antisymmetric) for c in ctrls]
    span_sp = s.span_spacing
    if swapped and not (isinstance(span_sp, str) and span_sp in ("linear", "cosine")):
        from ..geometry import _spacing
        span_sp = 1.0 - _spacing(s.nspan, span_sp)[::-1]   # keep the clustering on the same physical end
    return Surface(s.name, tuple(float(x) for x in r), rc, tuple(float(x) for x in t), tc,
                   s.nchord, s.nspan, s.chord_spacing, span_sp, ctrls)
