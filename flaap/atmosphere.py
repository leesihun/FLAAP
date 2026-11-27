"""International Standard Atmosphere (ISA 1976, 0-86 km geometric approx.).

Replaces the Fortran ``air_properties`` / ``altitute`` routines (module
``ATMOS``) which read ``atmosphere_input.txt`` and produced the density
``RHO`` and speed of sound.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

R_AIR = 287.05287       # J/(kg K)
GAMMA = 1.4
G0 = 9.80665            # m/s^2
R_EARTH = 6356766.0     # m (for geopotential altitude)

# (base geopotential altitude [m], base temperature [K], lapse rate [K/m])
_LAYERS = [
    (0.0, 288.15, -0.0065),
    (11000.0, 216.65, 0.0),
    (20000.0, 216.65, 0.001),
    (32000.0, 228.65, 0.0028),
    (47000.0, 270.65, 0.0),
    (51000.0, 270.65, -0.0028),
    (71000.0, 214.65, -0.002),
    (84852.0, 186.946, 0.0),
]
_P0 = 101325.0


def _base_pressures():
    p = [_P0]
    for i in range(1, len(_LAYERS)):
        h0, t0, lapse = _LAYERS[i - 1]
        h1 = _LAYERS[i][0]
        if lapse == 0.0:
            p.append(p[-1] * math.exp(-G0 * (h1 - h0) / (R_AIR * t0)))
        else:
            t1 = t0 + lapse * (h1 - h0)
            p.append(p[-1] * (t1 / t0) ** (-G0 / (lapse * R_AIR)))
    return p


_PB = _base_pressures()


@dataclass(frozen=True)
class Atmosphere:
    altitude: float      # geometric altitude [m]
    temperature: float   # K
    pressure: float      # Pa
    density: float       # kg/m^3
    speed_of_sound: float  # m/s
    viscosity: float     # Pa s (Sutherland)

    def velocity(self, mach: float) -> float:
        return mach * self.speed_of_sound

    def dynamic_pressure(self, velocity: float) -> float:
        return 0.5 * self.density * velocity ** 2


def isa(altitude: float = 0.0, delta_temp: float = 0.0) -> Atmosphere:
    """Standard atmosphere at geometric altitude ``altitude`` [m]."""
    h = R_EARTH * altitude / (R_EARTH + altitude)  # geopotential
    if h < -5000 or h > 86000:
        raise ValueError("altitude outside ISA model range (-5..86 km)")
    i = 0
    for j in range(len(_LAYERS)):
        if h >= _LAYERS[j][0]:
            i = j
    h0, t0, lapse = _LAYERS[i]
    pb = _PB[i]
    t = t0 + lapse * (h - h0)
    if lapse == 0.0:
        p = pb * math.exp(-G0 * (h - h0) / (R_AIR * t0))
    else:
        p = pb * (t / t0) ** (-G0 / (lapse * R_AIR))
    t_real = t + delta_temp
    rho = p / (R_AIR * t_real)
    a = math.sqrt(GAMMA * R_AIR * t_real)
    mu = 1.458e-6 * t_real ** 1.5 / (t_real + 110.4)
    return Atmosphere(altitude, t_real, p, rho, a, mu)


def altitude_from_density(rho: float) -> float:
    """Inverse of ISA density (bisection), used for matched-point flutter."""
    lo, hi = -4999.0, 85999.0
    if not (isa(hi).density <= rho <= isa(lo).density):
        raise ValueError("density outside ISA range")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if isa(mid).density > rho:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
