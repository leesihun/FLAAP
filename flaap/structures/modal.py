"""Structural models in modal form.

Two containers are used by every analysis module:

``FEModel``
    finite-element model with 6 DOF per node (ux, uy, uz, rx, ry, rz),
    sparse/dense ``K`` and ``M`` and a list of constrained DOFs.  It can
    compute normal modes and supplies the mass matrix for inertial loads
    (summation-of-forces) and root reactions (mode-displacement loads).

``ModalModel``
    the reduced model used by flutter / gust / trim: natural frequencies,
    generalized masses, modal damping and mode shapes at the structural
    nodes (translations and rotations).  Legacy FLAAP modal input files
    (``M_K.txt`` + ``MODESHAPE_input.txt``) and Nastran ``.f06`` eigenvector
    output are read into this container (see :mod:`flaap.io`).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass
class FEModel:
    nodes: np.ndarray                    # (Nn,3)
    K: object                            # (6Nn,6Nn) dense or sparse
    M: object
    fixed: np.ndarray = field(default_factory=lambda: np.zeros(0, int))  # constrained DOFs
    name: str = "fe"
    elements: dict = field(default_factory=dict)  # optional info (for plotting)

    @property
    def ndof(self) -> int:
        return 6 * self.nodes.shape[0]

    def free_dofs(self) -> np.ndarray:
        mask = np.ones(self.ndof, bool)
        mask[self.fixed] = False
        # also drop DOFs with no stiffness and no mass (e.g. unused drilling)
        Kd = self.K.diagonal() if sp.issparse(self.K) else np.diag(self.K)
        Md = self.M.diagonal() if sp.issparse(self.M) else np.diag(self.M)
        mask &= (np.abs(Kd) > 0) | (np.abs(Md) > 0)
        return np.where(mask)[0]

    def modes(self, n: int, name: str | None = None, damping: float = 0.0) -> "ModalModel":
        free = self.free_dofs()
        K = self.K[free][:, free] if not sp.issparse(self.K) else self.K.tocsr()[free][:, free]
        M = self.M[free][:, free] if not sp.issparse(self.M) else self.M.tocsr()[free][:, free]
        if sp.issparse(K) and K.shape[0] > 400:
            w2, v = spla.eigsh(K.tocsc(), k=n, M=M.tocsc(), sigma=0.0, which="LM")
        else:
            Kd = K.toarray() if sp.issparse(K) else K
            Md = M.toarray() if sp.issparse(M) else M
            w2, v = sla.eigh(Kd, Md, subset_by_index=[0, n - 1])
        order = np.argsort(w2)
        w2, v = w2[order], v[:, order]
        full = np.zeros((self.ndof, n))
        full[free] = v
        Mfull = self.M
        gm = np.einsum("in,in->n", full, (Mfull @ full))
        # mass-normalise
        full /= np.sqrt(gm)
        omega = np.sqrt(np.maximum(w2, 0.0))
        return ModalModel.from_dof_vectors(self.nodes, full, omega, np.ones(n),
                                           damping=np.full(n, damping), fe=self,
                                           name=name or self.name)


@dataclass
class ModalModel:
    nodes: np.ndarray            # (Nn,3) structural grid (aero frame)
    shapes: np.ndarray           # (Nn,6,nm)  ux,uy,uz,rx,ry,rz per mode
    omega: np.ndarray            # (nm,) natural circular frequencies [rad/s]
    gen_mass: np.ndarray         # (nm,)
    damping: np.ndarray | None = None   # (nm,) structural damping ratio zeta
    fe: FEModel | None = None
    name: str = "modal"
    labels: list[str] | None = None

    def __post_init__(self):
        self.omega = np.asarray(self.omega, float)
        self.gen_mass = np.asarray(self.gen_mass, float)
        if self.damping is None:
            self.damping = np.zeros_like(self.omega)
        self.damping = np.asarray(self.damping, float)

    # ------------------------------------------------------------------
    @classmethod
    def from_dof_vectors(cls, nodes, vecs, omega, gen_mass, **kw):
        nn = nodes.shape[0]
        sh = vecs.reshape(nn, 6, -1)
        return cls(np.asarray(nodes, float), sh, omega, gen_mass, **kw)

    @property
    def nmodes(self) -> int:
        return self.omega.size

    @property
    def freq_hz(self) -> np.ndarray:
        return self.omega / (2 * np.pi)

    @property
    def M(self) -> np.ndarray:
        return np.diag(self.gen_mass)

    @property
    def K(self) -> np.ndarray:
        return np.diag(self.gen_mass * self.omega ** 2)

    def C(self) -> np.ndarray:
        """Viscous modal damping matrix 2 zeta w m."""
        return np.diag(2 * self.damping * self.omega * self.gen_mass)

    def dof_vectors(self) -> np.ndarray:
        return self.shapes.reshape(-1, self.nmodes)

    def select(self, idx) -> "ModalModel":
        idx = np.atleast_1d(idx)
        return ModalModel(self.nodes, self.shapes[:, :, idx], self.omega[idx],
                          self.gen_mass[idx], self.damping[idx], self.fe, self.name,
                          [self.labels[i] for i in idx] if self.labels else None)

    def uz(self) -> np.ndarray:
        return self.shapes[:, 2, :]

    def table(self) -> str:
        lines = [" mode   f [Hz]     w [rad/s]   gen.mass     zeta"]
        for i in range(self.nmodes):
            lines.append(f"{i + 1:5d} {self.freq_hz[i]:10.4f} {self.omega[i]:11.4f} "
                         f"{self.gen_mass[i]:11.4e} {self.damping[i]:8.4f}")
        return "\n".join(lines)


def rigid_body_shapes(nodes: np.ndarray, ref: np.ndarray, kinds=("plunge", "pitch")) -> np.ndarray:
    """Rigid-body DOF vectors (Nn,6,nr) about point ``ref``.

    ``plunge`` = +z translation, ``pitch`` = nose-up rotation about +y
    (h = -(x - x_ref) theta), ``roll`` = rotation about +x (right wing up)
    """
    nodes = np.asarray(nodes, float)
    d = nodes - np.asarray(ref, float)[None, :]
    out = []
    for k in kinds:
        s = np.zeros((nodes.shape[0], 6))
        if k == "plunge":
            s[:, 2] = 1.0
        elif k == "pitch":
            # nose-up rotation theta about -y?  use h = -(x-xref)*theta, ry = -theta? keep
            # FE convention ry positive = rotation about +y => uz = -x*ry; nose-up = +ry
            s[:, 2] = -d[:, 0]
            s[:, 0] = d[:, 2]
            s[:, 4] = 1.0
        elif k == "roll":
            s[:, 2] = d[:, 1]
            s[:, 1] = -d[:, 2]
            s[:, 3] = 1.0
        elif k == "surge":
            s[:, 0] = 1.0
        else:
            raise ValueError(k)
        out.append(s)
    return np.stack(out, axis=2)
