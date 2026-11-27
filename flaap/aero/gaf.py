"""Generalized aerodynamic forces (port of ``calculation.f90``).

For a set of "aero modes" (structural modes splined to the boxes, rigid
body modes, control-surface rotations) the generalized aerodynamic force
matrix at reduced frequency ``k = w b / U`` is::

    W(k)   = dh/dx + i (k / b) h            (normalwash at 3/4-chord, w/U)
    dCp(k) = D(k)^-1 W(k)                   (box pressure jumps)
    Q(k)   = H_f^T diag(A) dCp(k)           (generalized forces / q)

so that the aerodynamic generalized force on mode ``i`` is ``q Q_ij eta_j``.
Gust: for a vertical gust ``w_g = W_g exp(i w (t - (x - x_0)/U))`` the
equivalent normalwash column is ``-cos(gamma) exp(-i k (x - x_0) / b)``
per unit ``W_g / U`` (exactly the Fortran ``gwash`` with the gust
reference at the trailing edge replaced by an arbitrary ``x_0``).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.interpolate import CubicSpline

from ..geometry import AeroMesh
from .dlm import DLM


@dataclass
class AeroModes:
    """Mode shapes evaluated on the aero mesh."""

    h_c: np.ndarray       # (N,nm) normal displacement at collocation points
    hx_c: np.ndarray      # (N,nm) slope at collocation points
    h_f: np.ndarray       # (N,nm) normal displacement at force points
    names: list[str] = field(default_factory=list)

    @property
    def nm(self) -> int:
        return self.h_c.shape[1]

    def __add__(self, other: "AeroModes") -> "AeroModes":
        return AeroModes(np.hstack([self.h_c, other.h_c]), np.hstack([self.hx_c, other.hx_c]),
                         np.hstack([self.h_f, other.h_f]), self.names + other.names)

    def select(self, idx) -> "AeroModes":
        idx = list(np.atleast_1d(idx))
        return AeroModes(self.h_c[:, idx], self.hx_c[:, idx], self.h_f[:, idx],
                         [self.names[i] for i in idx] if self.names else [])

    @classmethod
    def from_spline(cls, spline, shapes, names=None):
        hc, hx, hf = spline.modal(shapes)
        nm = shapes.shape[-1]
        return cls(hc, hx, hf, list(names) if names else [f"mode{i + 1}" for i in range(nm)])


def rigid_aero_modes(mesh: AeroMesh, x_ref: float, kinds=("plunge", "pitch"), y_ref=0.0,
                     bref: float = 1.0) -> AeroModes:
    """Rigid / trim-parameter aero modes.

    * ``plunge``  h = 1 (up)
    * ``pitch``   h = -(x - x_ref) (nose up; at k = 0 identical to ``alpha``)
    * ``roll``    h = (y - y_ref) (right wing up)
    * ``alpha``   angle of attack without displacement, w/U = -1
    * ``pitch_rate``  non-dimensional q_hat = Q b_ref / U: w/U = -(x - x_ref) q_hat / b_ref
    * ``roll_rate``   non-dimensional p_hat = P b_ref / U: w/U = (y - y_ref) p_hat / b_ref

    Rate and ``alpha`` modes have no displacement (``h = 0``), so their
    normalwash is independent of k (quasi-steady parameters for trim).
    """
    nz = mesh.normal[:, 2]
    ny = mesh.normal[:, 1]
    xc, yc, zc = mesh.colloc[:, 0], mesh.colloc[:, 1], mesh.colloc[:, 2]
    xf, yf, zf = mesh.force_point[:, 0], mesh.force_point[:, 1], mesh.force_point[:, 2]
    zero = np.zeros(mesh.n)
    cols = {"c": [], "x": [], "f": []}
    for k in kinds:
        if k == "plunge":
            hc, hx, hf = nz.copy(), zero, nz.copy()
        elif k == "pitch":
            hc, hx, hf = -(xc - x_ref) * nz, -nz, -(xf - x_ref) * nz
        elif k == "roll":
            hc = (yc - y_ref) * nz - zc * ny
            hx = zero
            hf = (yf - y_ref) * nz - zf * ny
        elif k == "alpha":
            hc, hx, hf = zero, -nz, zero
        elif k == "pitch_rate":
            hc, hx, hf = zero, -(xc - x_ref) * nz / bref, zero
        elif k == "roll_rate":
            hc, hx, hf = zero, ((yc - y_ref) * nz - zc * ny) / bref, zero
        else:
            raise ValueError(k)
        cols["c"].append(hc)
        cols["x"].append(hx)
        cols["f"].append(hf)
    return AeroModes(np.array(cols["c"]).T, np.array(cols["x"]).T, np.array(cols["f"]).T, list(kinds))


def control_aero_mode(mesh: AeroMesh, name: str) -> AeroModes:
    """Rigid rotation of a control surface (TE down positive, 1 rad)."""
    xh = mesh.control[name]
    on = np.isfinite(xh)
    hc = np.where(on, -(mesh.colloc[:, 0] - np.nan_to_num(xh)), 0.0)
    hx = np.where(on, -1.0, 0.0)
    hf = np.where(on, -(mesh.force_point[:, 0] - np.nan_to_num(xh)), 0.0)
    # boxes straddling the hinge would need partial treatment; boxes are
    # assigned wholly by their centre (Nastran AESURF convention).
    return AeroModes(hc[:, None], hx[:, None], hf[:, None], [name])


@dataclass
class GAFTable:
    """Tabulated GAFs ``Q(k)``, gust columns ``Qg(k)`` and pressures."""

    k: np.ndarray            # (nk,)
    Q: np.ndarray            # (nk,nm,nm)
    Qg: np.ndarray | None    # (nk,nm)
    mach: float
    bref: float
    names: list[str]
    cp: np.ndarray | None = None    # (nk,N,nm)  box dCp per unit mode amplitude
    cpg: np.ndarray | None = None   # (nk,N)

    def interp(self, k: float):
        """Cubic interpolation of Q (and Qg) at reduced frequency ``k``."""
        k = float(k)
        if self.k.size == 1:
            return self.Q[0], (None if self.Qg is None else self.Qg[0])
        if not hasattr(self, "_spl"):
            self._spl = CubicSpline(self.k, self.Q, axis=0)
            self._splg = None if self.Qg is None else CubicSpline(self.k, self.Qg, axis=0)
        kk = min(max(k, self.k[0]), self.k[-1]) if k < self.k[0] else k
        Q = self._spl(kk)
        Qg = None if self._splg is None else self._splg(kk)
        return Q, Qg

    def subset(self, idx) -> "GAFTable":
        idx = list(np.atleast_1d(idx))
        return GAFTable(self.k, self.Q[:, idx][:, :, idx],
                        None if self.Qg is None else self.Qg[:, idx], self.mach, self.bref,
                        [self.names[i] for i in idx],
                        None if self.cp is None else self.cp[:, :, idx], self.cpg)

    def save(self, path: str):
        np.savez(path, k=self.k, Q=self.Q, Qg=self.Qg if self.Qg is not None else np.zeros(0),
                 mach=self.mach, bref=self.bref, names=np.array(self.names))

    @classmethod
    def load(cls, path: str) -> "GAFTable":
        d = np.load(path, allow_pickle=False)
        Qg = d["Qg"] if d["Qg"].size else None
        return cls(d["k"], d["Q"], Qg, float(d["mach"]), float(d["bref"]), list(d["names"]))

    def write_text(self, path: str):
        """Human-readable dump similar to the Fortran Q output."""
        with open(path, "w") as f:
            f.write(f"# FLAAP GAF table  M={self.mach}  bref={self.bref}  nm={len(self.names)}\n")
            f.write("# modes: " + " ".join(self.names) + "\n")
            for i, k in enumerate(self.k):
                f.write(f"k = {k:.6f}\n")
                for r in range(self.Q.shape[1]):
                    f.write(" ".join(f"{v.real:14.6e}{v.imag:+14.6e}j" for v in self.Q[i, r]) + "\n")
                if self.Qg is not None:
                    f.write("gust: " + " ".join(f"{v.real:14.6e}{v.imag:+14.6e}j" for v in self.Qg[i]) + "\n")


class AeroModel:
    """DLM aerodynamics attached to a set of aero modes."""

    def __init__(self, mesh: AeroMesh, modes: AeroModes, mach: float, bref: float,
                 x_gust0: float | None = None, dlm: DLM | None = None,
                 aic_provider=None):
        self.mesh = mesh
        self.modes = modes
        self.mach = mach
        self.bref = bref
        self.x_gust0 = float(mesh.corners[:, :, 0].min()) if x_gust0 is None else x_gust0
        self.dlm = dlm if (dlm is not None or aic_provider is not None) else DLM(mesh, mach)
        # optional alternative aerodynamic operator: function k -> D(k) (e.g. TSD-corrected)
        self.aic_provider = aic_provider

    def aic(self, k: float) -> np.ndarray:
        if self.aic_provider is not None:
            return self.aic_provider(k)
        return self.dlm.aic(k, self.bref)

    def normalwash(self, k: float) -> np.ndarray:
        return self.modes.hx_c + 1j * (k / self.bref) * self.modes.h_c

    def gust_normalwash(self, k: float) -> np.ndarray:
        x = self.mesh.colloc[:, 0]
        return -np.cos(self.mesh.gamma) * np.exp(-1j * (k / self.bref) * (x - self.x_gust0))

    def solve(self, k: float, gust: bool = True):
        D = self.aic(k)
        W = self.normalwash(k)
        rhs = np.hstack([W, self.gust_normalwash(k)[:, None]]) if gust else W
        cp = np.linalg.solve(D, rhs)
        HA = self.modes.h_f * self.mesh.area[:, None]
        Qall = HA.T @ cp
        # half models (symmetry): generalized forces refer to the modelled half only
        if gust:
            return Qall[:, :-1], Qall[:, -1], cp[:, :-1], cp[:, -1]
        return Qall, None, cp, None

    def table(self, kvals, gust: bool = True, keep_cp: bool = True, verbose=False) -> GAFTable:
        kvals = np.atleast_1d(np.asarray(kvals, float))
        nm = self.modes.nm
        Q = np.zeros((kvals.size, nm, nm), complex)
        Qg = np.zeros((kvals.size, nm), complex) if gust else None
        cp = np.zeros((kvals.size, self.mesh.n, nm), complex) if keep_cp else None
        cpg = np.zeros((kvals.size, self.mesh.n), complex) if (keep_cp and gust) else None
        for i, k in enumerate(kvals):
            q, qg, c, cg = self.solve(k, gust)
            Q[i] = q
            if gust:
                Qg[i] = qg
            if keep_cp:
                cp[i] = c
                if gust:
                    cpg[i] = cg
            if verbose:
                print(f"  GAF k={k:.4f} done")
        return GAFTable(kvals, Q, Qg, self.mach, self.bref, list(self.modes.names), cp, cpg)


def default_k_list(kmax: float = 2.0, n: int = 16) -> np.ndarray:
    """Reduced frequencies clustered at low k (good for RFA and PK)."""
    t = np.linspace(0, 1, n)
    return kmax * t ** 1.6
