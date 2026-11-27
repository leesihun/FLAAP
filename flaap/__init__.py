"""FLAAP - Flutter and LCO Aeroelastic Analysis Program (Python edition).

Python port and extension of the original Fortran FLAAP code
(Intel Visual Fortran + MKL, 2017-2018).  The package provides

* lifting-surface paneling (port of ``aero_paneling_VLM.f90``)
* steady VLM and unsteady Doublet-Lattice aerodynamics (Albano-Rodden,
  quartic kernel approximation) - replaces the missing ``DLM.f90``
* generalized aerodynamic forces (port of ``calculation.f90``)
* Roger rational function approximation (``generalizedforce.f90``)
* flutter solutions: K (V-g), PK and p (state-space) methods with
  automatic flutter / divergence identification
* time-domain aeroservoelastic state-space models, discrete and
  continuous gust loads, time-domain maneuver loads
* static aeroelastic trim, divergence and control effectiveness
* a transonic small-disturbance (TSD) aerodynamic solver
* structural models: legacy modal files, MSC Nastran .f06, ANSYS and
  OpenRadioss readers, beam / plate / shell / solid FE, STL / STEP
  geometry -> FE modal analysis
* control-surface freeplay LCO (``nonlinear``) and aeroservoelasticity
  (``ase``)
* an AI-CAE4ALL-style config launcher (``launcher``)
* validation cases (Goland wing, AGARD 445.6, Conyers plate wing,
  Conner typical section, Leissa plate)
"""

__version__ = "2.0.0"

from .atmosphere import isa, Atmosphere

__all__ = ["Atmosphere", "isa", "__version__"]
