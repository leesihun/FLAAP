"""Goland wing (Goland 1945; Goland & Luke 1948) validation model.

Uniform cantilever wing, half-span 20 ft, chord 6 ft, elastic axis at
33 % chord, c.g. at 43 % chord::

    m    = 35.71 kg/m          I_ea = 8.64 kg m  (about the elastic axis)
    EI   = 9.77e6 N m^2        GJ   = 0.987e6 N m^2

Reference results (incompressible):

* natural frequencies (exact coupled bending-torsion beam, Goland 1945;
  tabulated in Savino et al., Energies 14 (2021) Table 2):
  48.13, 95.76, 243.79, 347.65 rad/s
* flutter, strip theory / Theodorsen at rho = 1.225 kg/m^3:
  V_f = 137.2 m/s, w_f = 70.7 rad/s (Goland 1945; Patil, Hodges & Cesnik
  2001 report 135.6 m/s / 70.2 rad/s)
* flutter with 3-D lifting-surface aerodynamics is higher than strip theory
  because of tip relief.  The published 3-D results are quoted at
  rho = 1.02 kg/m^3 (``DLM_REFERENCES_RHO102``): Murua et al. 2012 (SHARP
  UVLM) 165 m/s, Wang et al. 2006 163.8 m/s (UVLM) and 174.3 m/s (ZAERO),
  Savino et al. 2021 (DUST) 168.2 m/s (vortex lattice) and 174.2 m/s
  (surface panels).  FLAAP's DLM is compared with these at the same density.

The "modified Goland wing" (with tip store, Beran et al. 2004) used by the
legacy FLAAP input deck is a different configuration; see
``flaap.models.legacy``.
"""
from __future__ import annotations

import numpy as np

from ..geometry import Surface
from ..structures.beam import BeamProperties, build_beam_model
from ..structures.modal import FEModel, ModalModel

SPAN = 6.096
CHORD = 1.8288
X_EA = 0.33 * CHORD
X_CG = 0.43 * CHORD
# 8.64 kg m is the pitch inertia about the elastic axis (this reproduces the
# classical strip-theory flutter result 137 m/s / 70.7 rad/s)
PROPS = BeamProperties(EI=9.77e6, GJ=0.987e6, m=35.71,
                       I_cg=8.64 - 35.71 * (X_CG - X_EA) ** 2, x_ea=X_EA, x_cg=X_CG)

REFERENCE = {
    "omega_rad_s": [48.13, 95.76, 243.79, 347.65],
    "strip_flutter_speed": 137.2,
    "strip_flutter_omega": 70.7,
    "density": 1.225,
    "dlm_density": 1.02,
    "dlm_flutter_speed_range": (163.8, 174.3),
}

# published 3-D flutter results for the Goland wing at rho = 1.02 kg/m^3
# (source, method, V_f [m/s], f_f [Hz] or None)
DLM_REFERENCES_RHO102 = [
    ("Murua et al. 2012", "UVLM (SHARP)", 165.0, 10.98),
    ("Wang et al. 2006", "UVLM", 163.8, None),
    ("Wang et al. 2006", "ZAERO (DLM)", 174.3, None),
    ("Savino et al. 2021", "DUST vortex lattice", 168.2, 10.84),
    ("Savino et al. 2021", "DUST surface panels", 174.2, 11.06),
]


def goland_surface(nchord: int = 8, nspan: int = 16, spacing: str = "linear") -> Surface:
    return Surface("goland", (0.0, 0.0, 0.0), CHORD, (0.0, SPAN, 0.0), CHORD,
                   nchord=nchord, nspan=nspan, chord_spacing=spacing)


def goland_fe(nelem: int = 20, tip_store: dict | None = None) -> FEModel:
    y = np.linspace(0.0, SPAN, nelem + 1)
    pm = [tip_store] if tip_store else None
    return build_beam_model(y, PROPS, point_masses=pm, name="Goland wing")


def goland_modal(nmodes: int = 4, nelem: int = 20, damping: float = 0.0) -> ModalModel:
    return goland_fe(nelem).modes(nmodes, damping=damping)
