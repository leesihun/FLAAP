"""Readers / runner for the input files of the original Fortran FLAAP.

Files (all list-directed, one value per line unless noted):

``LEE_input.txt``
    1 gust velocity, 2 gust acceleration, 3 wing "x" length (span),
    4 wing "y" length (chord), 5-12 planform corners ``(X, Y)`` of the
    inboard TE, outboard TE, outboard LE, inboard LE, 13 thickness,
    14 Mach number, 15-17 reduced frequency min / step / max (GAF table),
    18 number of Roger lag terms, 19 number of panels along X (span),
    20 number of panels along Y (chord), 21-23 reduced frequency
    min / increment / max (flutter search).

    **Legacy axes**: ``X`` is spanwise and ``Y`` chordwise with ``Y = 0``
    at the trailing edge (the LE has the larger ``Y``; the flow goes in
    ``-Y``).  FLAAP axes are ``x = Y_LE,root - Y`` (downstream),
    ``y = X``.

``M_K.txt``
    one line per mode: natural frequency [Hz], generalized mass,
    generalized stiffness.

``MODESHAPE_input.txt``
    line 1: number of rows ``3 Np``, line 2: number of modes; then ``3 Np``
    rows with one value per mode: for every grid point ``(T3, R1, R2)``
    (only the ``T3`` rows are used by the original program; rotations may
    be dummies).  Grid point order (``calculation.f90``): span stations
    ``i = 1..nX`` (the root row ``i = 0`` is clamped and omitted) outer,
    chord stations ``j = 0..nY`` inner - i.e. the aerodynamic lattice
    corners.

``atmosphere_input.txt``
    1 Mach number (informational), 2 altitude [m].

:func:`run_legacy` reproduces the program flow (lattice, corner-node
interpolation of the mode shapes, DLM GAFs, PK flutter) and writes
``Vg_save.dat``, ``Vf_save.dat`` and ``pk_results_mode*.dat``.
"""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np

from ..atmosphere import isa
from ..geometry import Surface, build_mesh
from ..structures.modal import ModalModel
from ..structures.spline import legacy_corner_matrices, surface_spline


def _numbers(path: str) -> list[float]:
    with open(path, encoding="utf-8", errors="replace") as f:
        out = []
        for tok in f.read().split():
            try:
                out.append(float(tok.replace("D", "E").replace("d", "e")))
            except ValueError:
                break
        return out


@dataclass
class LegacyInput:
    gust_velocity: float
    gust_accel: float
    length_x: float            # span
    length_y: float            # chord
    corners: np.ndarray        # (4,2) legacy (X, Y): in-TE, out-TE, out-LE, in-LE
    thickness: float
    mach: float
    k_table: tuple             # (kmin, kstep, kmax)
    nlag: int
    n_span: int
    n_chord: int
    k_search: tuple            # (kmin, kinc, kmax)

    # -------------------------------------------------------- geometry
    @property
    def x_le_root(self) -> float:
        return float(self.corners[3, 1])

    def to_flaap(self, X, Y):
        """Legacy (X span, Y chord-from-TE) -> FLAAP (x, y, 0)."""
        X = np.asarray(X, float)
        Y = np.asarray(Y, float)
        return np.stack([self.x_le_root - Y, X, np.zeros_like(X)], axis=-1)

    def surface(self, name: str = "wing") -> Surface:
        ite, ote, ole, ile = self.corners
        root_le = self.to_flaap(ile[0], ile[1])
        tip_le = self.to_flaap(ole[0], ole[1])
        return Surface(name, tuple(root_le), float(ile[1] - ite[1]), tuple(tip_le),
                       float(ole[1] - ote[1]), nchord=self.n_chord, nspan=self.n_span)

    def grid(self) -> np.ndarray:
        """(nX * (nY+1), 3) legacy structural grid in FLAAP axes, in file order."""
        ite, ote, ole, ile = self.corners
        pts = []
        for i in range(1, self.n_span + 1):
            s = i / self.n_span
            X = ite[0] + s * (ote[0] - ite[0])
            yte = ite[1] + s * (ote[1] - ite[1])
            yle = ile[1] + s * (ole[1] - ile[1])
            for j in range(self.n_chord + 1):
                pts.append(self.to_flaap(X, yte + (yle - yte) * j / self.n_chord))
        return np.array(pts)

    def k_list(self) -> np.ndarray:
        k0, dk, k1 = self.k_table
        if dk <= 0:
            return np.array([k0, k1])
        return np.round(np.arange(k0, k1 + 0.5 * dk, dk), 12)


def read_lee_input(path: str) -> LegacyInput:
    v = _numbers(path)
    if len(v) < 22:
        raise ValueError(f"{path}: expected >= 22 values, got {len(v)}")
    corners = np.array(v[4:12]).reshape(4, 2)
    ks = tuple(v[20:23]) if len(v) >= 23 else (v[20], v[21], v[16])
    return LegacyInput(v[0], v[1], v[2], v[3], corners, v[12], v[13], (v[14], v[15], v[16]),
                       int(v[17]), int(v[18]), int(v[19]), ks)


def read_m_k(path: str):
    """Returns ``(freq_hz, gen_mass, gen_stiffness)`` arrays."""
    a = np.atleast_2d(np.loadtxt(path))
    return a[:, 0], a[:, 1], a[:, 2]


def read_modeshape(path: str):
    """Returns ``S (3 Np, nm)`` (raw rows) of ``MODESHAPE_input.txt``."""
    with open(path, encoding="utf-8", errors="replace") as f:
        tok = f.read().split()
    nrow, nm = int(float(tok[0])), int(float(tok[1]))
    data = np.array([float(t) for t in tok[2:2 + nrow * nm]])
    if data.size != nrow * nm:
        raise ValueError(f"{path}: expected {nrow * nm} values, got {data.size}")
    return data.reshape(nrow, nm)


def read_atmosphere(path: str):
    v = _numbers(path)
    return (v[0] if v else 0.0), (v[1] if len(v) > 1 else 0.0)


@dataclass
class LegacyCase:
    inp: LegacyInput
    freq_hz: np.ndarray
    gen_mass: np.ndarray
    gen_stiff: np.ndarray
    S: np.ndarray
    mach_atm: float
    altitude: float
    directory: str = "."

    @property
    def nmodes(self) -> int:
        return self.freq_hz.size

    def modal(self, length_scale: float = 1.0, nmodes: int | None = None) -> ModalModel:
        """Modal model on the legacy grid (``uz`` only; rotations are not
        used by the original program).  ``omega`` from the generalized
        stiffness / mass (consistent with ``M_K.txt``)."""
        nm = self.nmodes if nmodes is None else nmodes
        nodes = self.inp.grid()
        npt = nodes.shape[0]
        z = self.S[0::3, :nm][:npt] * length_scale
        if z.shape[0] != npt:
            raise ValueError(f"MODESHAPE has {z.shape[0]} points, grid has {npt}")
        sh = np.zeros((npt, 6, nm))
        sh[:, 2, :] = z
        omega = np.sqrt(self.gen_stiff[:nm] / self.gen_mass[:nm])
        return ModalModel(nodes, sh, omega, self.gen_mass[:nm], name="legacy",
                          labels=[f"mode{i + 1}" for i in range(nm)])


def read_case(directory: str = ".") -> LegacyCase:
    j = lambda f: os.path.join(directory, f)  # noqa: E731
    inp = read_lee_input(j("LEE_input.txt"))
    f, m, k = read_m_k(j("M_K.txt"))
    S = read_modeshape(j("MODESHAPE_input.txt"))
    nm = min(f.size, S.shape[1])
    mach, alt = read_atmosphere(j("atmosphere_input.txt")) if os.path.exists(j("atmosphere_input.txt")) else (0.0, 0.0)
    return LegacyCase(inp, f[:nm], m[:nm], k[:nm], S[:, :nm], mach, alt, directory)


def legacy_aero(case: LegacyCase, modal: ModalModel, spline: str = "legacy", mach: float | None = None):
    """Lattice + aerodynamic modes of a legacy case.

    ``spline='legacy'``: the corner interpolation of ``calculation.f90``;
    ``'ips'``: infinite-plate spline through the grid points (+ clamped root).
    """
    from ..aero.gaf import AeroModel, AeroModes
    inp = case.inp
    mesh = build_mesh([inp.surface()], "symmetric")
    if spline == "legacy":
        sp = legacy_corner_matrices(mesh, modal.nodes)
        vecs = modal.dof_vectors()
    else:
        # add the clamped root row (zero displacement) to the spline points
        ite, ile = inp.corners[0], inp.corners[3]
        root = np.array([inp.to_flaap(ite[0], ite[1] + (ile[1] - ite[1]) * j / inp.n_chord)
                         for j in range(inp.n_chord + 1)])
        nodes = np.vstack([root, modal.nodes])
        sp = surface_spline(mesh, nodes, method=spline)
        vecs = np.vstack([np.zeros((6 * root.shape[0], modal.nmodes)), modal.dof_vectors()])
    modes = AeroModes.from_spline(sp, vecs, modal.labels)
    bref = 0.5 * float(inp.corners[3, 1] - inp.corners[0, 1])
    return AeroModel(mesh, modes, inp.mach if mach is None else mach, bref)


def run_legacy(directory: str = ".", velocities=None, method: str = "pk", output: str | None = None,
               spline: str = "legacy", length_scale: float = 1.0, kvals=None, density: float | None = None,
               nmodes: int | None = None, verbose: bool = False):
    """Run a Fortran-FLAAP case directory with the Python solver.

    Returns ``(result, case, table)``; writes the legacy result files to
    ``output`` (if given).
    """
    from ..flutter import identify, k_method, pk_method
    case = read_case(directory)
    modal = case.modal(length_scale, nmodes)
    am = legacy_aero(case, modal, spline)
    kv = case.inp.k_list() if kvals is None else np.asarray(kvals, float)
    table = am.table(kv, gust=True, keep_cp=False, verbose=verbose)
    rho = isa(case.altitude).density if density is None else density
    M = np.diag(modal.gen_mass)
    K = np.diag(modal.gen_mass * modal.omega ** 2)
    V = np.linspace(5.0, 300.0, 296) if velocities is None else np.asarray(velocities, float)
    if method == "pk":
        res = pk_method(M, K, table, rho, V)
    else:
        res = k_method(M, K, table, rho)
    identify(res)
    if output:
        os.makedirs(output, exist_ok=True)
        res.write_legacy(output)
    return res, case, table
