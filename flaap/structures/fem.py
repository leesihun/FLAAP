"""Generic 3-D finite-element assembly and modal solution.

A :class:`FEMesh` holds node coordinates and a list of
:class:`ElementBlock` (element type + connectivity + material +
thickness).  :func:`assemble` builds sparse ``K`` and ``M`` in the FLAAP
6-DOF layout (``ux uy uz rx ry rz`` per node), :func:`build_fe_model`
adds boundary conditions and returns a :class:`~flaap.structures.modal.FEModel`,
and :func:`solve_modes` computes normal modes with
``scipy.sparse.linalg.eigsh`` in shift-invert mode and returns a
:class:`~flaap.structures.modal.ModalModel` ready for the splines and the
flutter solvers.

Element types (see :mod:`flaap.structures.shell` and
:mod:`flaap.structures.solid` for the formulations):

=========  =====================================  ==========================
kind       element                                 nodal DOFs used
=========  =====================================  ==========================
``quad4``  flat shell: MITC4 + QM6 + drilling     6
``tri3``   flat shell: DKT + Allman + drilling    6
``tet4``   linear tetrahedron                      3 (translations)
``tet10``  quadratic tetrahedron                   3
``hex8``   brick with incompatible modes           3
=========  =====================================  ==========================

DOFs without stiffness and mass (rotations of nodes attached only to
solids, unused nodes) are removed before the eigen-solution.  The rotations
of such nodes in the returned mode shapes are estimated from the
translation field (``rotations="lsq"``, default: least-squares local
displacement gradient, ``theta = curl(u) / 2``) or set to zero
(``rotations="zero"``).  The flutter splines of FLAAP use the
translations (IPS / TPS / RBF on ``uz``); the rotations matter only for
beam splines and are exact for rigid-body motions and linear fields.

Shift-invert: ``eigsh(K, k, M, sigma)`` with a small negative ``sigma``
(default ``-1`` (rad/s)^2) so that free-free structures (singular ``K``)
are handled; rigid-body modes come out with ``omega ~ 0`` and can be
dropped with ``skip_rigid=True``.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .modal import FEModel, ModalModel
from .plate import OrthotropicMaterial
from .shell import shell_quad_matrices, shell_tri_matrices
from .solid import (SOLID_FACES, hex8_matrices, rotations_from_translations, tet4_matrices,
                    tet10_matrices)

__all__ = ["Material", "MATERIALS", "get_material", "ElementBlock", "FEMesh", "assemble",
           "build_fe_model", "dofs_of_nodes", "select_nodes", "solve_modes", "modal_analysis",
           "surface_triangles", "SHELL_KINDS", "SOLID_KINDS", "NODES_PER_ELEMENT"]

SHELL_KINDS = ("quad4", "tri3")
SOLID_KINDS = ("tet4", "tet10", "hex8")
NODES_PER_ELEMENT = {"quad4": 4, "tri3": 3, "tet4": 4, "tet10": 10, "hex8": 8}


# ----------------------------------------------------------------- materials
@dataclass
class Material:
    """Isotropic linear-elastic material (SI: Pa, -, kg/m^3)."""

    E: float
    nu: float
    rho: float
    name: str = "material"

    @property
    def G(self) -> float:
        return self.E / (2 * (1 + self.nu))

    def ortho(self, **kw) -> OrthotropicMaterial:
        """Plane-stress shell material (:class:`OrthotropicMaterial`)."""
        return OrthotropicMaterial.isotropic(self.E, self.nu, self.rho, **kw)


#: Preset materials (typical handbook values, SI).  ``cfrp_quasi_iso`` is an
#: in-plane quasi-isotropic carbon/epoxy laminate treated as isotropic.
MATERIALS: dict[str, Material] = {
    "aluminium_7075": Material(71.7e9, 0.33, 2810.0, "Al 7075-T6"),
    "aluminium_6061": Material(68.9e9, 0.33, 2700.0, "Al 6061-T6"),
    "aluminium_2024": Material(73.1e9, 0.33, 2780.0, "Al 2024-T3"),
    "steel": Material(200e9, 0.29, 7850.0, "structural steel"),
    "stainless_steel": Material(193e9, 0.29, 8000.0, "stainless steel 304"),
    "titanium_6al4v": Material(113.8e9, 0.342, 4430.0, "Ti-6Al-4V"),
    "cfrp_quasi_iso": Material(50e9, 0.31, 1600.0, "CFRP quasi-isotropic"),
    "gfrp_quasi_iso": Material(20e9, 0.30, 1900.0, "GFRP quasi-isotropic"),
}
_ALIASES = {"aluminum_7075": "aluminium_7075", "al7075": "aluminium_7075", "aluminium": "aluminium_7075",
            "aluminum": "aluminium_7075", "al6061": "aluminium_6061", "aluminum_6061": "aluminium_6061",
            "al2024": "aluminium_2024", "aluminum_2024": "aluminium_2024", "titanium": "titanium_6al4v",
            "ti6al4v": "titanium_6al4v", "cfrp": "cfrp_quasi_iso", "gfrp": "gfrp_quasi_iso"}


def get_material(spec) -> Material | OrthotropicMaterial:
    """Material from a preset name, a dict ``{"E":, "nu":, "rho":}``, a
    :class:`Material` or an :class:`OrthotropicMaterial` (returned as is)."""
    if isinstance(spec, (Material, OrthotropicMaterial)):
        return spec
    if isinstance(spec, str):
        key = spec.strip().lower().replace(" ", "_").replace("-", "_")
        key = _ALIASES.get(key, key)
        if key not in MATERIALS:
            raise KeyError(f"unknown material '{spec}'; presets: {sorted(MATERIALS)}")
        return MATERIALS[key]
    if isinstance(spec, dict):
        d = {k.lower(): v for k, v in spec.items()}
        if "e1" in d:
            return OrthotropicMaterial(**{k: v for k, v in spec.items()})
        return Material(float(d["e"]), float(d["nu"]), float(d["rho"]), str(d.get("name", "material")))
    raise TypeError(f"cannot interpret material {spec!r}")


def _shell_mat(m) -> OrthotropicMaterial:
    return m if isinstance(m, OrthotropicMaterial) else m.ortho()


def _solid_mat(m) -> tuple[float, float, float]:
    if isinstance(m, OrthotropicMaterial):
        if not (np.isclose(m.E1, m.E2) and np.isclose(m.G12, m.E1 / (2 * (1 + m.nu12)), rtol=1e-3)):
            warnings.warn("orthotropic material on a solid element: isotropic E1, nu12 used")
        return m.E1, m.nu12, m.rho
    return m.E, m.nu, m.rho


# ----------------------------------------------------------------- mesh
@dataclass
class ElementBlock:
    """A group of elements of one type and material.

    ``thickness`` (shells only): scalar, per element (ne,) or nodal values
    per element (ne, nnode).  ``ref_dir``: global direction defining the
    shell material axis 1 (projected, then rotated by ``material.angle``).
    """

    kind: str
    conn: np.ndarray
    material: object
    thickness: object = None
    name: str = ""
    ref_dir: tuple = (1.0, 0.0, 0.0)
    lumped: bool | None = None

    def __post_init__(self):
        self.kind = self.kind.lower()
        if self.kind not in NODES_PER_ELEMENT:
            raise ValueError(f"unknown element kind '{self.kind}'")
        self.conn = np.atleast_2d(np.asarray(self.conn, int))
        if self.conn.shape[1] != NODES_PER_ELEMENT[self.kind]:
            raise ValueError(f"{self.kind}: expected {NODES_PER_ELEMENT[self.kind]} nodes per element")
        self.material = get_material(self.material)
        if self.kind in SHELL_KINDS and self.thickness is None:
            raise ValueError(f"shell block '{self.name}' needs a thickness")

    @property
    def n(self) -> int:
        return self.conn.shape[0]


@dataclass
class FEMesh:
    """Nodes (Nn,3) + element blocks.  ``node_ids`` optional external ids."""

    nodes: np.ndarray
    blocks: list = field(default_factory=list)
    node_ids: np.ndarray | None = None
    name: str = "mesh"

    def __post_init__(self):
        self.nodes = np.asarray(self.nodes, float)

    @property
    def nn(self) -> int:
        return self.nodes.shape[0]

    def add(self, kind, conn, material, thickness=None, name="", **kw) -> "FEMesh":
        self.blocks.append(ElementBlock(kind, conn, material, thickness, name, **kw))
        return self

    def conns(self, kinds=None) -> list[np.ndarray]:
        return [b.conn for b in self.blocks if kinds is None or b.kind in kinds]

    def surface_triangles(self) -> np.ndarray:
        return surface_triangles(self)

    def mass(self) -> float:
        K, M = assemble(self, stiffness=False)
        u = np.zeros(6 * self.nn)
        u[0::6] = 1.0
        return float(u @ (M @ u))

    def index_of(self, ids) -> np.ndarray:
        """Row indices of external node ids."""
        if self.node_ids is None:
            return np.asarray(ids, int)
        lut = {int(g): i for i, g in enumerate(self.node_ids)}
        return np.array([lut[int(g)] for g in np.atleast_1d(ids)], int)


def _coo(ndof, dofs, Ke):
    nd = dofs.shape[1]
    rows = np.repeat(dofs, nd, axis=1).ravel()
    cols = np.tile(dofs, (1, nd)).ravel()
    return rows, cols, Ke.ravel()


def assemble(mesh: FEMesh, lumped: bool = False, drill: float = 1e-3, stiffness: bool = True,
             chunk: int = 20000):
    """Sparse global ``K, M`` (6 Nn x 6 Nn, CSR) of all element blocks."""
    nn = mesh.nn
    ndof = 6 * nn
    rk, ck, vk, rm, cm, vm = [], [], [], [], [], []
    for b in mesh.blocks:
        lump = lumped if b.lumped is None else b.lumped
        for s in range(0, b.n, chunk):
            conn = b.conn[s:s + chunk]
            X = mesh.nodes[conn]
            if b.kind in SHELL_KINDS:
                t = np.asarray(b.thickness, float)
                if t.ndim >= 1 and t.shape[0] == b.n:
                    t = t[s:s + chunk]
                f = shell_quad_matrices if b.kind == "quad4" else shell_tri_matrices
                Ke, Me = f(X, t, _shell_mat(b.material), lumped=lump, drill=drill, ref_dir=b.ref_dir)
                comp = np.arange(6)
            else:
                E, nu, rho = _solid_mat(b.material)
                f = {"tet4": tet4_matrices, "tet10": tet10_matrices, "hex8": hex8_matrices}[b.kind]
                Ke, Me = f(X, E, nu, rho, lumped=lump)
                comp = np.arange(3)
            dofs = (6 * conn[:, :, None] + comp[None, None, :]).reshape(conn.shape[0], -1)
            if stiffness:
                r, c, v = _coo(ndof, dofs, Ke)
                rk.append(r), ck.append(c), vk.append(v)
            r, c, v = _coo(ndof, dofs, Me)
            rm.append(r), cm.append(c), vm.append(v)

    def build(r, c, v):
        if not r:
            return sp.csr_matrix((ndof, ndof))
        A = sp.coo_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))),
                          shape=(ndof, ndof)).tocsr()
        return (0.5 * (A + A.T)).tocsr()

    K = build(rk, ck, vk) if stiffness else None
    return K, build(rm, cm, vm)


# ----------------------------------------------------------------- boundary conditions
def dofs_of_nodes(nodes_idx, components="123456") -> np.ndarray:
    """Global DOF indices of node rows ``nodes_idx`` for the Nastran-style
    component string (``"123"`` = translations, ``"456"`` = rotations)."""
    comp = np.array(sorted({int(c) - 1 for c in str(components)}), int)
    idx = np.asarray(nodes_idx, int).ravel()
    return (6 * idx[:, None] + comp[None, :]).ravel()


def select_nodes(nodes: np.ndarray, *, below: tuple | None = None, above: tuple | None = None,
                 plane: tuple | None = None, box: tuple | None = None, ids=None, node_ids=None,
                 tol: float | None = None) -> np.ndarray:
    """Row indices of nodes satisfying all given criteria.

    * ``below=(axis, value)``: coordinate ``<= value + tol`` (axis 0/1/2 or
      "x"/"y"/"z"); ``above`` likewise ``>= value - tol``
    * ``plane=(point, normal)``: distance to the plane ``<= tol``
    * ``box=(lo, hi)``: inside the axis-aligned box (inclusive, + tol)
    * ``ids``: external node ids (requires ``node_ids``) or row indices
    ``tol`` defaults to 1e-6 of the model size.
    """
    nodes = np.asarray(nodes, float)
    size = max(float(np.ptp(nodes, axis=0).max()), 1e-12)
    tol = 1e-6 * size if tol is None else tol
    mask = np.ones(nodes.shape[0], bool)
    ax = lambda a: {"x": 0, "y": 1, "z": 2}.get(a, a) if isinstance(a, str) else int(a)  # noqa: E731
    if below is not None:
        mask &= nodes[:, ax(below[0])] <= below[1] + tol
    if above is not None:
        mask &= nodes[:, ax(above[0])] >= above[1] - tol
    if plane is not None:
        p, n = np.asarray(plane[0], float), np.asarray(plane[1], float)
        n = n / np.linalg.norm(n)
        mask &= np.abs((nodes - p) @ n) <= tol
    if box is not None:
        lo, hi = np.asarray(box[0], float), np.asarray(box[1], float)
        mask &= np.all((nodes >= lo - tol) & (nodes <= hi + tol), axis=1)
    if ids is not None:
        sel = np.zeros_like(mask)
        if node_ids is not None:
            lut = {int(g): i for i, g in enumerate(node_ids)}
            sel[[lut[int(g)] for g in np.atleast_1d(ids) if int(g) in lut]] = True
        else:
            sel[np.asarray(ids, int)] = True
        mask &= sel
    return np.where(mask)[0]


def build_fe_model(mesh: FEMesh, fixed=None, lumped: bool = False, drill: float = 1e-3,
                   name: str | None = None) -> FEModel:
    """Assemble ``mesh`` into an :class:`FEModel`.

    ``fixed``: constrained global DOF indices (see :func:`dofs_of_nodes`).
    DOFs without stiffness and mass (solid-node rotations) are dropped
    automatically by :meth:`FEModel.free_dofs` / :func:`solve_modes`.
    ``FEModel.elements`` holds ``{"mesh": mesh, kind: conn, ...}``.
    """
    K, M = assemble(mesh, lumped=lumped, drill=drill)
    fixed = np.unique(np.asarray([] if fixed is None else fixed, int))
    el = {"mesh": mesh}
    for b in mesh.blocks:
        el.setdefault(b.kind, [])
        el[b.kind].append(b.conn)
    for k in list(el):
        if k != "mesh":
            el[k] = np.concatenate(el[k], 0)
    return FEModel(mesh.nodes, K, M, fixed, name=name or mesh.name, elements=el)


# ----------------------------------------------------------------- eigen-solution
def solve_modes(fe: FEModel, nmodes: int, sigma: float = -1.0, rotations: str = "lsq",
                skip_rigid: bool = False, rigid_hz: float = 1e-2, damping: float = 0.0,
                name: str | None = None, mesh: FEMesh | None = None) -> ModalModel:
    """Normal modes of ``fe`` (mass-normalized, ``gen_mass = 1``).

    ``sigma`` [(rad/s)^2] shift for the shift-invert Lanczos solution
    (``eigsh``); a negative value keeps ``K - sigma M`` positive definite
    for free-free models.  ``skip_rigid``: drop modes below ``rigid_hz``
    (``nmodes`` elastic modes are then returned).  ``rotations``:
    ``"lsq"`` / ``"zero"`` for nodes whose rotation DOFs are not part of
    the model (solid nodes).  The returned :class:`ModalModel` has
    ``fe`` set (``fe.elements["mesh"]`` gives the mesh for plotting).
    """
    mesh = mesh if mesh is not None else fe.elements.get("mesh") if fe.elements else None
    ndof = fe.ndof
    K = fe.K.tocsr() if sp.issparse(fe.K) else sp.csr_matrix(fe.K)
    M = fe.M.tocsr() if sp.issparse(fe.M) else sp.csr_matrix(fe.M)
    act = np.ones(ndof, bool)
    act[fe.fixed] = False
    act &= (np.abs(K.diagonal()) > 0) | (np.abs(M.diagonal()) > 0)
    free = np.where(act)[0]
    Kf = K[free][:, free]
    Mf = M[free][:, free]
    nreq = nmodes + (6 if skip_rigid else 0)
    nreq = min(nreq, free.size - 1)
    if free.size <= 600:
        w2, v = sla.eigh(Kf.toarray(), Mf.toarray(), subset_by_index=[0, nreq - 1])
    else:
        w2, v = spla.eigsh(Kf.tocsc(), k=nreq, M=Mf.tocsc(), sigma=sigma, which="LM")
    order = np.argsort(w2)
    w2, v = w2[order], v[:, order]
    if skip_rigid:
        keep = np.sqrt(np.abs(w2)) / (2 * np.pi) >= rigid_hz
        w2, v = w2[keep][:nmodes], v[:, keep][:, :nmodes]
    nm = w2.size
    full = np.zeros((ndof, nm))
    full[free] = v
    gm = np.einsum("in,in->n", full, M @ full)
    full /= np.sqrt(gm)
    # sign convention: largest translation component positive.  Symmetric
    # structures give ties (e.g. both tip corners of a torsion mode, equal to
    # ~1e-12), so the first DOF within 1e-6 of the maximum is used: the sign
    # must not depend on round-off of the eigensolver.
    tr = full.reshape(-1, 6, nm)[:, :3, :].reshape(-1, nm)
    atr = np.abs(tr)
    pick = np.argmax(atr >= (1.0 - 1e-6) * atr.max(axis=0), axis=0)
    sgn = np.sign(tr[pick, np.arange(nm)])
    full *= np.where(sgn == 0, 1.0, sgn)
    shapes = full.reshape(-1, 6, nm)
    # rotations of nodes without rotational DOFs
    rot_act = act.reshape(-1, 6)[:, 3:].any(axis=1) | np.isin(
        np.arange(fe.nodes.shape[0]), fe.fixed // 6)   # clamped nodes keep 0
    missing = ~rot_act
    if rotations == "lsq" and missing.any() and mesh is not None:
        conns = mesh.conns()
        rot = rotations_from_translations(fe.nodes, conns, shapes[:, :3, :], nodes_mask=missing)
        shapes[missing, 3:, :] = rot[missing]
    elif rotations not in ("lsq", "zero"):
        raise ValueError("rotations must be 'lsq' or 'zero'")
    omega = np.sqrt(np.maximum(w2, 0.0))
    mm = ModalModel(np.asarray(fe.nodes, float), shapes, omega, np.ones(nm),
                    damping=np.full(nm, damping), fe=fe, name=name or fe.name,
                    labels=[f"mode{i + 1}" for i in range(nm)])
    return mm


def modal_analysis(mesh: FEMesh, nmodes: int, fixed=None, lumped: bool = False,
                   drill: float = 1e-3, **kw) -> ModalModel:
    """``build_fe_model`` + ``solve_modes`` in one call."""
    fe = build_fe_model(mesh, fixed=fixed, lumped=lumped, drill=drill)
    return solve_modes(fe, nmodes, mesh=mesh, **kw)


# ----------------------------------------------------------------- visualisation
def surface_triangles(mesh: FEMesh) -> np.ndarray:
    """Triangulated outer surface (nt,3) of node row indices, for plotting.

    Shell quads are split into two triangles; for solid blocks the faces
    that belong to exactly one element (the boundary) are returned with
    outward orientation (tet10 boundary faces are split into 4 linear
    triangles through the mid-edge nodes).
    """
    tris = []
    faces = []
    for b in mesh.blocks:
        c = b.conn
        if b.kind == "tri3":
            tris.append(c)
        elif b.kind == "quad4":
            tris.append(c[:, [0, 1, 2]])
            tris.append(c[:, [0, 2, 3]])
        else:
            for f in SOLID_FACES[b.kind]:
                fc = c[:, list(f)]
                if b.kind == "tet10":
                    # mid-edge node of each face edge
                    mids = []
                    from .solid import TET10_EDGES
                    for a, bb in ((f[0], f[1]), (f[1], f[2]), (f[2], f[0])):
                        k = TET10_EDGES.index(tuple(sorted((a, bb))))
                        mids.append(c[:, 4 + k])
                    fc = np.column_stack([fc] + mids)
                faces.append(fc)
    if faces:
        # group by size
        by = {}
        for f in faces:
            by.setdefault(f.shape[1], []).append(f)
        for k, lst in by.items():
            F = np.concatenate(lst, 0)
            key = np.sort(F[:, :3] if k == 6 else F, axis=1)
            _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
            B = F[cnt[inv.ravel()] == 1]
            if k == 3:
                tris.append(B)
            elif k == 4:
                tris.append(B[:, [0, 1, 2]])
                tris.append(B[:, [0, 2, 3]])
            else:      # tet10 face: corners a b c, mids ab bc ca
                a, b2, c2, ab, bc, ca = B.T
                tris += [np.stack(t, 1) for t in ((a, ab, ca), (ab, b2, bc), (ca, bc, c2), (ab, bc, ca))]
    if not tris:
        return np.zeros((0, 3), int)
    return np.concatenate(tris, 0)
