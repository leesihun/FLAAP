"""Gust load analysis: discrete 1-cos gusts and continuous turbulence.

Gust inputs are the vertical gust velocity ``w_g`` (positive *up*) at the
gust reference point ``x_0`` of the aerodynamic model (by default the most
forward point of the mesh); the penetration delay of every box is part of
the gust GAF column ``Q_g``.

Discrete gust (CS-25.341(a) / FAR 25.341(a))::

    w_g(s) = U_ds / 2 (1 - cos(pi s / H)),   0 <= s <= 2H,  s = V t
    U_ds   = U_ref F_g (H / 106.68 m)^(1/6)      (H = 9 .. 106.68 m)

``U_ref`` is 17.07 m/s EAS at sea level, linearly reduced to 13.41 m/s at
4572 m and 6.36 m/s at 18288 m (halved at V_D).

Continuous turbulence (CS-25.341(b)) with the von Karman spectrum
(L = 762 m)::

    Phi(W) = sigma^2 L/pi (1 + 8/3 (1.339 L W)^2) / (1 + (1.339 L W)^2)^(11/6)

``A_bar = sigma_y / sigma_w`` and ``N_0`` (zero crossings) are computed from
the frequency response of the time-domain model or directly from the
tabulated GAFs (:func:`frequency_response_gust`).  The design load is
``y_1g +/- U_sigma A_bar`` with ``U_sigma`` = 27.43 m/s TAS (sea level)
reducing to 24.08 m/s TAS at 7315 m.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import CubicSpline

from .atmosphere import isa

H_MAX = 106.68          # 350 ft
H_MIN = 9.144           # 30 ft
L_VK = 762.0            # 2500 ft


# ----------------------------------------------------------------------
# certification gust definitions
# ----------------------------------------------------------------------
def cs25_uref(altitude: float, at_vd: bool = False) -> float:
    """Reference gust velocity U_ref [m/s EAS] (CS-25.341(a)(5))."""
    h = float(altitude)
    if h <= 4572.0:
        u = 17.07 + (13.41 - 17.07) * h / 4572.0
    else:
        u = 13.41 + (6.36 - 13.41) * min(h - 4572.0, 18288.0 - 4572.0) / (18288.0 - 4572.0)
    return 0.5 * u if at_vd else u


def cs25_fg(zmo: float, mlw: float, mtow: float, mzfw: float) -> float:
    """Flight profile alleviation factor F_g at sea level (CS-25.341(a)(6))."""
    R1, R2 = mlw / mtow, mzfw / mtow
    Fgz = 1.0 - zmo / 76200.0
    Fgm = np.sqrt(R2 * np.tan(np.pi * R1 / 4.0))
    return 0.5 * (Fgz + Fgm)


def fg_at_altitude(fg_sl: float, altitude: float, zmo: float) -> float:
    """F_g increases linearly from the sea level value to 1.0 at Z_mo."""
    if zmo <= 0:
        return 1.0
    return float(fg_sl + (1.0 - fg_sl) * min(max(altitude / zmo, 0.0), 1.0))


def design_gust_velocity(H: float, altitude: float, fg: float = 1.0, at_vd: bool = False,
                         tas: bool = True) -> float:
    """U_ds for gust gradient ``H`` [m] (half the gust length)."""
    u = cs25_uref(altitude, at_vd) * fg * (H / H_MAX) ** (1.0 / 6.0)
    if tas:
        u *= np.sqrt(isa(0.0).density / isa(altitude).density)
    return u


def cs25_usigma(altitude: float, at_vd: bool = False) -> float:
    """Continuous turbulence reference intensity U_sigma_ref [m/s TAS]."""
    h = float(altitude)
    u = 27.43 + (24.08 - 27.43) * min(h, 7315.0) / 7315.0
    return 0.5 * u if at_vd else u


def gust_gradients(n: int = 8, hmin: float = H_MIN, hmax: float = H_MAX) -> np.ndarray:
    """Gust gradients spread geometrically over the certification range."""
    return np.geomspace(hmin, hmax, n)


def one_minus_cosine(t, V: float, H: float, Uds: float, t0: float = 0.0):
    """1-cos gust velocity and its time derivative at the gust reference."""
    t = np.asarray(t, float)
    s = V * (t - t0)
    on = (s >= 0) & (s <= 2 * H)
    wg = np.where(on, 0.5 * Uds * (1 - np.cos(np.pi * s / H)), 0.0)
    wgd = np.where(on, 0.5 * Uds * np.pi * V / H * np.sin(np.pi * s / H), 0.0)
    return wg, wgd


# ----------------------------------------------------------------------
# turbulence spectra
# ----------------------------------------------------------------------
def von_karman_psd(omega, V: float, sigma: float = 1.0, L: float = L_VK):
    """One-sided vertical gust PSD per rad/s (integrates to sigma^2)."""
    Om = np.asarray(omega, float) / V
    a = (1.339 * L * Om) ** 2
    return sigma ** 2 * L / np.pi * (1 + 8.0 / 3.0 * a) / (1 + a) ** (11.0 / 6.0) / V


def dryden_psd(omega, V: float, sigma: float = 1.0, L: float = L_VK):
    Om = np.asarray(omega, float) / V
    a = (L * Om) ** 2
    return sigma ** 2 * L / np.pi * (1 + 3 * a) / (1 + a) ** 2 / V


# ----------------------------------------------------------------------
# time-domain discrete gust
# ----------------------------------------------------------------------
@dataclass
class GustResult:
    t: np.ndarray
    y: np.ndarray            # (nt, nout)
    labels: list[str]
    wg: np.ndarray
    H: float
    Uds: float

    def peak(self, name: str):
        i = self.labels.index(name)
        return float(self.y[:, i].max()), float(self.y[:, i].min())


def discrete_gust(ss, V: float, H: float, Uds: float, t_end: float | None = None,
                  dt: float | None = None, sign: float = 1.0) -> GustResult:
    """Simulate a 1-cos gust on a state space from
    :meth:`flaap.timedomain.AeroelasticSystem.state_space` (with gust)."""
    if "wg" not in ss.inputs:
        raise ValueError("state space has no gust input")
    wmax = np.max(np.abs(ss.eig().imag)) if ss.n else 1.0
    dt = dt or min(2 * H / V / 50.0, 2 * np.pi / max(wmax, 1e-6) / 20.0)
    t_end = t_end or (2 * H / V + 2.0)
    t = np.arange(0.0, t_end + dt / 2, dt)
    wg, wgd = one_minus_cosine(t, V, H, sign * Uds)
    u = np.zeros((t.size, len(ss.inputs)))
    u[:, ss.inputs.index("wg")] = wg
    u[:, ss.inputs.index("wg_dot")] = wgd
    y, _ = ss.simulate(t, u)
    return GustResult(t, y, list(ss.outputs), wg, H, sign * Uds)


def gust_envelope(system, V: float, rho: float, gradients, uds_func, outputs=None,
                  t_end: float | None = None) -> dict:
    """Run a family of 1-cos gusts and return peak/min of each output.

    ``uds_func(H) -> U_ds`` (TAS).  Returns ``{'H': ..., 'max': (nH, nout),
    'min': ..., 'labels': ..., 'critical_H': {label: H}}``.
    """
    ss = system.state_space(V, rho, gust=True)
    labels = ss.outputs if outputs is None else list(outputs)
    idx = ss.index(labels)
    mx, mn = [], []
    for H in gradients:
        r = discrete_gust(ss, V, H, uds_func(H), t_end=t_end)
        mx.append(r.y[:, idx].max(axis=0))
        mn.append(r.y[:, idx].min(axis=0))
    mx, mn = np.array(mx), np.array(mn)
    amp = np.maximum(np.abs(mx), np.abs(mn))
    crit = {lab: float(np.asarray(gradients)[np.argmax(amp[:, j])]) for j, lab in enumerate(labels)}
    return {"H": np.asarray(gradients), "max": mx, "min": mn, "labels": labels, "critical_H": crit}


# ----------------------------------------------------------------------
# frequency domain (directly on tabulated GAFs, no RFA)
# ----------------------------------------------------------------------
def frequency_response_gust(M, K, C, table, V: float, rho: float, omega, load_rows=None,
                            inertia=None, stiffness=None):
    """Modal response per unit gust velocity ``w_g`` (m/s) at ``omega``.

    Returns ``eta (nw, nm)`` and, if ``load_rows=(rows (nk,nr,nm), rows_g (nk,nr))``
    is given, the station loads ``(nw, nr)`` (aero + ``inertia @ eta''`` +
    ``stiffness @ eta``).  Reduced frequencies beyond the table are clamped
    (keep the excitation spectrum inside the tabulated range).
    """
    omega = np.atleast_1d(np.asarray(omega, float))
    q = 0.5 * rho * V * V
    b = table.bref
    nm = M.shape[0]
    C = np.zeros_like(M) if C is None else C
    kmax = table.k.max()
    if load_rows is not None:
        rows, rows_g = load_rows
        sR = CubicSpline(table.k, rows, axis=0)
        sRg = CubicSpline(table.k, rows_g, axis=0)
    eta = np.zeros((omega.size, nm), complex)
    loads = None if load_rows is None else np.zeros((omega.size, rows.shape[1]), complex)
    for i, w in enumerate(omega):
        k = min(w * b / V, kmax)
        Q, Qg = table.interp(k)
        Z = -w * w * M + 1j * w * C + K - q * Q[:nm, :nm]
        e = np.linalg.solve(Z, q * Qg[:nm] / V)
        eta[i] = e
        if loads is not None:
            L = q * (sR(k)[:, :nm] @ e + sRg(k) / V)
            if inertia is not None:
                L = L + inertia @ (-w * w * e)
            if stiffness is not None:
                L = L + stiffness @ e
            loads[i] = L
    return eta, loads


def fft_response(t, wg, H_of_omega):
    """Time response to the gust history ``wg(t)`` given a transfer function
    ``H_of_omega(omega) -> (nw, nout)``; zero-padded FFT (linear convolution)."""
    t = np.asarray(t, float)
    dt = t[1] - t[0]
    n = 2 * t.size
    W = np.fft.rfft(wg, n)
    om = 2 * np.pi * np.fft.rfftfreq(n, dt)
    Hm = H_of_omega(om)
    Y = Hm * W[:, None]
    y = np.fft.irfft(Y, n, axis=0)[:t.size]
    return y


# ----------------------------------------------------------------------
# continuous turbulence
# ----------------------------------------------------------------------
def turbulence_from_ss(ss, V: float, outputs=None, omega=None, L: float = L_VK,
                       spectrum: str = "von_karman") -> dict:
    """A_bar and N_0 of outputs of a state space (inputs ``wg``, ``wg_dot``)."""
    omega = np.linspace(1e-3, 2.5 * np.max(np.abs(ss.eig().imag)), 4000) if omega is None else omega
    H = ss.freq_response(omega, inputs=["wg", "wg_dot"], outputs=outputs)
    Ht = H[:, :, 0] + 1j * omega[:, None] * H[:, :, 1]
    labels = ss.outputs if outputs is None else list(outputs)
    return turbulence_stats(omega, Ht, V, labels, L, spectrum)


def turbulence_stats(omega, Hw, V: float, labels, L: float = L_VK,
                     spectrum: str = "von_karman") -> dict:
    """``Hw`` (nw, nout) transfer functions output / w_g."""
    psd = von_karman_psd if spectrum == "von_karman" else dryden_psd
    Phi = psd(omega, V, 1.0, L)
    Sy = np.abs(Hw) ** 2 * Phi[:, None]
    m0 = trapezoid(Sy, omega, axis=0)
    m2 = trapezoid(Sy * omega[:, None] ** 2, omega, axis=0)
    Abar = np.sqrt(m0)
    N0 = np.sqrt(m2 / np.maximum(m0, 1e-300)) / (2 * np.pi)
    return {"labels": list(labels), "Abar": Abar, "N0": N0, "omega": omega, "psd": Sy}
