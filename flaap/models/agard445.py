"""AGARD 445.6 wing, "weakened model 3" (Yates 1987, NASA TM-100492).

Planform (half model, right wing, root at the symmetry plane ``y = 0``,
root leading edge at the origin)::

    root chord 0.558677 m      tip chord 0.368580 m     taper 0.66
    semi-span  0.762 m         quarter-chord sweep 45 deg
    -> mid-chord sweep 43.15 deg, LE sweep 46.73 deg, TE sweep 34.6 deg

Note: the LE sweep that follows from the quarter-chord sweep, chords and
span above is 46.73 deg (tip LE at x = 0.8095 m); the value 46.31 deg that
is quoted in some papers is not consistent with these dimensions.  The
geometry here is defined by the 45 deg quarter-chord line.

Airfoil NACA 65A004 (symmetric, 4 % thick); the plate thickness at a node is
``chord(y) * t/c(x/c)`` from the tabulated 65A004 ordinates (pchip
interpolation).  The 65A004 thickness is exactly zero at the leading edge;
a floor ``t_min`` (default 0.1 mm, ~0.03 % of the tip chord; it lowers the
frequencies by < 0.1 % through the added LE mass) is applied.  Zero nodal thickness would not
make the element singular by itself (the Gauss points are interior), but
the floor keeps the LE nodal rotations well conditioned.

Equivalent orthotropic material of the NASA Nastran models of the weakened
model 3 (e.g. Kolonay 1996).  Yates (1987) itself lists E = 3.2455 GPa,
G = 0.4119 GPa, nu = 0.31 for his FE model; the set below gives a better
frequency match on the MITC4 plate::

    E1 = 3.1511 GPa (along the mid-chord / sweep line), E2 = 0.4162 GPa,
    G12 = 0.4392 GPa, nu12 = 0.31, rho = 381.98 kg/m^3

Structural model: MITC4 Reissner-Mindlin plate (:mod:`flaap.structures.plate`)
in the x-y plane with all DOFs of the root nodes clamped.  Without any
tuning the model gives (16x16 mesh / 32x32 mesh)::

    9.60 / 9.61, 39.97 / 39.75, 51.06 / 50.85, 96.84 / 95.74 Hz

i.e. within 1 % (bending), +4.9 % (1st torsion), +0.7 % (2nd bending) and
-1.7 % (2nd torsion) of the measured frequencies.  The structural panel
mass is 1.70 kg whereas the mass ratios of the flutter data imply 1.867 kg
(rho * mu * frustum volume; Yates gives a measured 1.8628 kg), so the plate
model is ~9 % lighter than the test article.

Natural frequencies, Yates (1987) NASA TM-100492 (1st bending, 1st torsion,
2nd bending, 2nd torsion):

* measured (ground vibration test): 9.60, 38.10, 50.70, 98.50 Hz
  (:data:`EXPERIMENT_FREQ_HZ`)
* Yates' FE model: 9.5992, 38.1650, 48.3482, 91.5448 Hz
  (:data:`YATES_FE_FREQ_HZ`; the "38.17 / 48.35 / 91.54 Hz" often quoted as
  experimental are these FE values).  Yates recommends g = 0.02 structural
  damping.

Flutter data, see :data:`FLUTTER_DATA` and its consistency notes.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator

from ..geometry import Surface
from ..structures.modal import FEModel, ModalModel
from ..structures.plate import OrthotropicMaterial, build_plate_model, plate_mesh_trapezoid

__all__ = ["ROOT_CHORD", "TIP_CHORD", "SEMISPAN", "SWEEP_QC_DEG", "SWEEP_LE_DEG",
           "SWEEP_MID_DEG", "TIP_LE_X", "B_REF", "OMEGA_ALPHA", "EXPERIMENT_FREQ_HZ",
           "YATES_FE_FREQ_HZ", "PUBLISHED_FE_FREQ_HZ", "FLUTTER_DATA", "MATERIAL", "naca65a004_thickness",
           "agard445_material", "agard445_thickness", "agard445_surface", "agard445_fe",
           "agard445_modal", "frustum_volume", "flutter_speed_index", "check_flutter_data"]

# ------------------------------------------------------------------ geometry
ROOT_CHORD = 0.558677
TIP_CHORD = 0.368580
SEMISPAN = 0.762
SWEEP_QC_DEG = 45.0
TIP_LE_X = SEMISPAN * np.tan(np.deg2rad(SWEEP_QC_DEG)) + 0.25 * (ROOT_CHORD - TIP_CHORD)
SWEEP_LE_DEG = float(np.rad2deg(np.arctan(TIP_LE_X / SEMISPAN)))
SWEEP_MID_DEG = float(np.rad2deg(np.arctan(
    (TIP_LE_X + 0.5 * TIP_CHORD - 0.5 * ROOT_CHORD) / SEMISPAN)))
B_REF = 0.5 * ROOT_CHORD          # root semichord used in the flutter speed index (0.2794 m)

EXPERIMENT_FREQ_HZ = np.array([9.60, 38.10, 50.70, 98.50])     # measured (Yates 1987)
YATES_FE_FREQ_HZ = np.array([9.5992, 38.1650, 48.3482, 91.5448])  # Yates' FE model
PUBLISHED_FE_FREQ_HZ = YATES_FE_FREQ_HZ                           # backwards-compatible alias
OMEGA_ALPHA = 239.3   # omega_alpha used by Yates to reduce FSI and omega_f/omega_alpha [rad/s]

# NACA 65A004 half-thickness ordinates, % chord (Abbott & von Doenhoff; the
# "A" modification gives straight sides aft of ~80 % chord, finite TE)
_XC = np.array([0.0, 0.5, 0.75, 1.25, 2.5, 5.0, 7.5, 10.0, 15.0, 20.0, 25.0, 30.0,
                35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0,
                90.0, 95.0, 100.0]) / 100.0
_YT = np.array([0.0, 0.311, 0.378, 0.481, 0.656, 0.877, 1.062, 1.216, 1.463, 1.649,
                1.790, 1.894, 1.962, 1.997, 1.996, 1.952, 1.867, 1.742, 1.584, 1.400,
                1.193, 0.966, 0.728, 0.490, 0.249, 0.009]) / 100.0
_YT_INTERP = PchipInterpolator(_XC, _YT)


def naca65a004_thickness(xc) -> np.ndarray:
    """Thickness-to-chord ratio t/c (full thickness) at chord fraction ``xc``."""
    xc = np.clip(np.asarray(xc, float), 0.0, 1.0)
    return 2.0 * _YT_INTERP(xc)


# ------------------------------------------------------------------ material
def agard445_material(**kw) -> OrthotropicMaterial:
    """Equivalent orthotropic material with axis 1 along the mid-chord line."""
    p = dict(E1=3.1511e9, E2=0.4162e9, G12=0.4392e9, nu12=0.31, rho=381.98,
             angle=OrthotropicMaterial.along_sweep(SWEEP_MID_DEG))
    p.update(kw)
    return OrthotropicMaterial(**p)


MATERIAL = agard445_material()


# ------------------------------------------------------------------ models
def _le_chord(y):
    eta = np.asarray(y, float) / SEMISPAN
    return eta * TIP_LE_X, ROOT_CHORD + eta * (TIP_CHORD - ROOT_CHORD)


def agard445_thickness(nodes: np.ndarray, t_min: float = 1e-4) -> np.ndarray:
    """Plate thickness at points (x, y): chord(y) * t/c(x/c), floored at t_min."""
    nodes = np.asarray(nodes, float)
    xle, c = _le_chord(nodes[:, 1])
    t = c * naca65a004_thickness((nodes[:, 0] - xle) / c)
    return np.maximum(t, t_min)


def agard445_surface(nchord: int = 8, nspan: int = 10, **kw) -> Surface:
    """Aerodynamic lifting surface (half model, root at y = 0; use
    ``build_mesh([...], symmetry='symmetric')``)."""
    return Surface("agard445", (0.0, 0.0, 0.0), ROOT_CHORD, (TIP_LE_X, SEMISPAN, 0.0),
                   TIP_CHORD, nchord=nchord, nspan=nspan, **kw)


def agard445_fe(n_chord: int = 16, n_span: int = 16, t_min: float = 1e-4,
                material: OrthotropicMaterial | None = None, chord_spacing="linear",
                membrane: bool = False) -> FEModel:
    """MITC4 plate model of the weakened model 3, clamped at the root."""
    nodes, quads = plate_mesh_trapezoid((0.0, 0.0), ROOT_CHORD, (TIP_LE_X, SEMISPAN),
                                        TIP_CHORD, n_chord, n_span,
                                        chord_spacing=chord_spacing)
    t = agard445_thickness(nodes, t_min)
    return build_plate_model(nodes, quads, t, material or MATERIAL, clamp_y=0.0,
                             membrane=membrane, name="AGARD 445.6 weakened 3")


def agard445_modal(nmodes: int = 4, **kw) -> ModalModel:
    """Normal modes of :func:`agard445_fe` (``kw`` passed to it)."""
    mm = agard445_fe(**kw).modes(nmodes, name="AGARD 445.6 weakened 3")
    labels = ["1st bending", "1st torsion", "2nd bending", "2nd torsion"]
    mm.labels = (labels + [f"mode {i + 1}" for i in range(4, nmodes)])[:nmodes]
    return mm


# ------------------------------------------------------------------ flutter data
def frustum_volume() -> float:
    """Volume of the conical frustum (root/tip chords as diameters, height =
    semispan) used to define the mass ratio mu = m / (rho V)."""
    R, r = 0.5 * ROOT_CHORD, 0.5 * TIP_CHORD
    return np.pi * SEMISPAN / 3.0 * (R * R + R * r + r * r)


def flutter_speed_index(V, mu, b=B_REF, omega_alpha=OMEGA_ALPHA):
    """FSI = V / (b omega_alpha sqrt(mu))."""
    return np.asarray(V) / (b * omega_alpha * np.sqrt(mu))


#: Experimental flutter points of the weakened model 3 in air, Yates (1987)
#: NASA TM-100492 Table II.  keys: Mach; values: ``rho`` [kg/m^3], ``V``
#: flutter speed [m/s], ``f`` flutter frequency [Hz]
#: (= ``freq_ratio`` x omega_alpha / 2 pi), ``FSI``, ``freq_ratio``
#: (omega_f / omega_alpha), ``mu`` (mass ratio).
#:
#: Consistency (see :func:`check_flutter_data`):
#:
#: * ``rho * mu`` = 14.315 kg/m^3 for all points (= panel mass / frustum
#:   volume, i.e. a panel mass of 1.867 kg).
#: * FSI recomputed from ``V, mu`` with b = 0.2794 m and omega_alpha =
#:   239.3 rad/s agrees with the tabulated FSI within 0.1 % at every Mach.
FLUTTER_DATA: dict[float, dict[str, float]] = {
    0.499: dict(rho=0.42776, V=172.46, FSI=0.4459, freq_ratio=0.5353, mu=33.465),
    0.678: dict(rho=0.20821, V=231.37, FSI=0.4174, freq_ratio=0.4722, mu=68.753),
    0.901: dict(rho=0.09947, V=296.69, FSI=0.3700, freq_ratio=0.4225, mu=143.920),
    0.954: dict(rho=0.06339, V=307.36, FSI=0.3059, freq_ratio=0.3807, mu=225.820),
    0.960: dict(rho=0.06339, V=309.01, FSI=0.3076, freq_ratio=0.3648, mu=225.820),
    1.072: dict(rho=0.05515, V=344.73, FSI=0.3201, freq_ratio=0.3623, mu=259.590),
    1.141: dict(rho=0.07834, V=364.33, FSI=0.4031, freq_ratio=0.4593, mu=182.740),
}
for _M, _d in FLUTTER_DATA.items():
    _d["f"] = round(_d["freq_ratio"] * OMEGA_ALPHA / (2 * np.pi), 3)
    _d["V_from_FSI"] = round(float(_d["FSI"] * B_REF * OMEGA_ALPHA * np.sqrt(_d["mu"])), 1)
del _M, _d


def check_flutter_data() -> str:
    """Table of the internal consistency checks of :data:`FLUTTER_DATA`."""
    vol = frustum_volume()
    lines = [f"frustum volume {vol:.5f} m^3, b_s = {B_REF:.4f} m, "
             f"omega_a = {OMEGA_ALPHA:.2f} rad/s",
             "   M      rho*mu   m_panel   FSI_tab  FSI_calc   err%   V/M"]
    for M, d in FLUTTER_DATA.items():
        fsi = float(flutter_speed_index(d["V"], d["mu"]))
        lines.append(f"{M:6.3f} {d['rho'] * d['mu']:9.4f} {d['rho'] * d['mu'] * vol:8.4f} "
                     f"{d['FSI']:8.4f} {fsi:8.4f} {100 * (fsi / d['FSI'] - 1):6.2f} "
                     f"{d['V'] / M:6.1f}")
    return "\n".join(lines)
