"""Normal modes straight from CAD / STL geometry.

Pipeline
--------
1. **Geometry input**

   * STL / OBJ / PLY / OFF (any format :mod:`trimesh` reads): a triangulated
     surface.  Duplicate vertices are merged; for the solid route the
     surface must be closed (watertight).
   * STEP / IGES / BREP: read with the gmsh OpenCASCADE kernel.  STEP and
     IGES carry units; by default gmsh is asked to convert to metres
     (``Geometry.OCCTargetUnit = "M"``).  BREP has no units (``length_scale``).

2. **Meshing**

   * ``kind="solid"``: gmsh Delaunay volume mesh of linear
     tetrahedra (``Mesh.MeshSizeMax = size``).  STL input is either used
     as-is (the STL triangles become the boundary mesh; ``remesh=False``)
     or reparametrised and remeshed at ``size`` (``remesh=True``,
     ``classifySurfaces`` + ``createGeometry``).  ``element="tet10"`` adds
     straight-sided mid-edge nodes (:func:`flaap.structures.solid.tet4_to_tet10`),
     which removes the shear locking of tet4 (see :mod:`flaap.structures.solid`).
   * ``kind="shell"``: the surface triangles themselves (``remesh=False``) or
     a gmsh 2-D remesh at ``size`` become DKT/Allman shell triangles
     (:mod:`flaap.structures.shell`) of thickness ``thickness`` - a thin-walled
     skin model of the closed surface (or a mid-surface model of an open one).

3. **Boundary conditions** via :func:`flaap.structures.fem.select_nodes`:
   ``clamp=dict(below=("y", 0.0))``, ``dict(plane=(p, n))``,
   ``dict(box=(lo, hi))`` or ``dict(ids=[...])`` (row indices); several dicts
   in a list are OR-ed.  ``components="123456"`` (default: all six; the
   rotations of solid nodes carry no stiffness and are ignored).

4. **Modal solution** (:func:`flaap.structures.fem.solve_modes`):
   shift-invert Lanczos, mass-normalised modes, rotations of solid nodes
   from a least-squares fit of the displacement gradient (``rotations="lsq"``)
   so that the shapes can be splined like any other FLAAP modal model.

The result :class:`GeometryModal` bundles the :class:`ModalModel`, the
:class:`FEMesh`, the clamped nodes and the boundary triangulation
(:func:`surface_triangulation`) for plotting (:func:`plot_mode`,
:func:`plot_modes`).  :func:`planform_nodes` picks a planform-uniform node
subset for :func:`flaap.structures.spline.surface_spline` (the dense IPS
system would otherwise be built on thousands of upper/lower-surface nodes).
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, field

import numpy as np

from .fem import FEMesh, build_fe_model, dofs_of_nodes, get_material, select_nodes, solve_modes
from .modal import ModalModel
from .solid import fix_tet10_order, tet4_to_tet10, tet_volumes

__all__ = ["GeometryModal", "load_surface", "mesh_geometry", "geometry_modal", "clamp_nodes",
           "surface_triangulation", "planform_nodes", "plot_mode", "plot_modes",
           "CAD_EXTENSIONS", "MESH_EXTENSIONS", "naca_wing_surface", "write_stl"]

CAD_EXTENSIONS = (".step", ".stp", ".iges", ".igs", ".brep", ".brp")
MESH_EXTENSIONS = (".stl", ".obj", ".ply", ".off", ".glb", ".gltf")


@dataclass
class GeometryModal:
    """Result of :func:`geometry_modal`."""

    modal: ModalModel
    mesh: FEMesh
    clamped: np.ndarray                      # row indices of clamped nodes
    info: dict = field(default_factory=dict)

    @property
    def surface(self) -> np.ndarray:
        """Boundary triangulation (nt,3) of node rows (for plotting)."""
        return surface_triangulation(self.mesh)

    def freq_hz(self) -> np.ndarray:
        return self.modal.freq_hz


# ----------------------------------------------------------------- gmsh helpers
class _Gmsh:
    """Context manager: initialise gmsh (quiet) unless already initialised."""

    def __init__(self, name="flaap", verbose=False):
        self.name, self.verbose = name, verbose

    def __enter__(self):
        import gmsh
        self.gmsh = gmsh
        self.own = not gmsh.isInitialized()
        if self.own:
            gmsh.initialize(interruptible=False)
        gmsh.option.setNumber("General.Terminal", 1 if self.verbose else 0)
        gmsh.option.setNumber("General.Verbosity", 5 if self.verbose else 1)
        gmsh.model.add(self.name)
        return gmsh

    def __exit__(self, *exc):
        try:
            self.gmsh.model.remove()
        finally:
            if self.own:
                self.gmsh.finalize()
        return False


def _ext(path) -> str:
    return os.path.splitext(str(path))[1].lower()


def _gmsh_nodes(gmsh):
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    xyz = np.asarray(xyz, float).reshape(-1, 3)
    tags = np.asarray(tags, np.int64)
    lut = np.full(int(tags.max()) + 1, -1, np.int64)
    lut[tags] = np.arange(tags.size)
    return xyz, lut


def _gmsh_elements(gmsh, etype, nper, lut):
    _, ntags = gmsh.model.mesh.getElementsByType(etype)
    if len(ntags) == 0:
        return np.zeros((0, nper), int)
    return lut[np.asarray(ntags, np.int64).reshape(-1, nper)].astype(int)


def _compress(nodes, conn):
    used, inv = np.unique(conn.ravel(), return_inverse=True)
    return nodes[used], inv.reshape(conn.shape)


def _set_sizes(gmsh, size, size_min=None):
    if size is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMax", float(size))
        gmsh.option.setNumber("Mesh.MeshSizeMin", float(size_min if size_min else 0.1 * size))
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Optimize", 1)


def _load_into_gmsh(gmsh, path, length_scale, occ_unit, remesh, angle_deg, verts=None, faces=None):
    """Bring geometry into the current gmsh model.  Returns True if CAD (OCC)."""
    ext = _ext(path) if path is not None else ".stl"
    if ext in CAD_EXTENSIONS:
        if occ_unit and ext not in (".brep", ".brp"):
            gmsh.option.setString("Geometry.OCCTargetUnit", occ_unit)
        gmsh.model.occ.importShapes(str(path))
        if length_scale != 1.0:
            gmsh.model.occ.dilate(gmsh.model.occ.getEntities(), 0, 0, 0,
                                  length_scale, length_scale, length_scale)
        gmsh.model.occ.synchronize()
        return True
    # triangulated surface -> discrete gmsh surface (via a temporary STL)
    if verts is None:
        verts, faces = load_surface(path, length_scale)
    import tempfile
    fd, tmp = tempfile.mkstemp(suffix=".stl")
    os.close(fd)
    try:
        _write_stl(tmp, verts, faces)
        gmsh.merge(tmp)
    finally:
        os.remove(tmp)
    gmsh.model.mesh.removeDuplicateNodes()
    if remesh:
        gmsh.model.mesh.classifySurfaces(np.deg2rad(angle_deg), True, True, np.pi)
        gmsh.model.mesh.createGeometry()
    gmsh.model.geo.synchronize()
    return False


def _write_stl(path, verts, faces):
    v = np.asarray(verts, float)[np.asarray(faces, int)]
    n = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-300)
    rec = np.zeros(len(v), dtype=[("n", "<f4", 3), ("v", "<f4", 9), ("a", "<u2")])
    rec["n"] = n
    rec["v"] = v.reshape(-1, 9)
    with open(path, "wb") as fh:
        fh.write(b"flaap binary stl".ljust(80, b" "))
        fh.write(np.uint32(len(v)).tobytes())
        fh.write(rec.tobytes())


def write_stl(path, verts, faces):
    """Write a binary STL (float32) of a triangulated surface."""
    _write_stl(path, verts, faces)


def naca_wing_surface(semispan: float, root_chord: float, tip_chord: float, sweep_le: float = 0.0,
                      tc: float = 0.12, n_chord: int = 24, n_span: int = 12,
                      root_le=(0.0, 0.0, 0.0)) -> tuple[np.ndarray, np.ndarray]:
    """Closed, watertight triangulation of a straight-tapered NACA 00xx wing.

    x chordwise, y spanwise (root at ``y = root_le[1]``), z thickness.
    ``n_chord`` points per surface (cosine spacing), ``n_span`` span stations,
    closed trailing edge (a4 = -0.1036), flat root/tip caps (centroid fans).
    ``sweep_le``: leading-edge sweep angle [deg].  Returns (vertices, faces)
    with outward normals.
    """
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n_chord + 1)))          # 0..1
    yt = 5 * tc * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3
                   - 0.1036 * x ** 4)
    yt[-1] = 0.0
    # closed loop: TE -> upper -> LE -> lower -> (TE)
    xu, zu = x[::-1], yt[::-1]
    xl, zl = x[1:-1], -yt[1:-1]
    loop = np.vstack([np.c_[xu, zu], np.c_[xl, zl]])           # (np,2), np = 2 n_chord
    npnt = loop.shape[0]
    eta = np.linspace(0, 1, n_span + 1)
    r0 = np.asarray(root_le, float)
    tan = np.tan(np.deg2rad(sweep_le))
    V = []
    for e in eta:
        c = root_chord + e * (tip_chord - root_chord)
        y = e * semispan
        V.append(np.c_[r0[0] + y * tan + c * loop[:, 0], np.full(npnt, r0[1] + y),
                       r0[2] + c * loop[:, 1]])
    V = np.vstack(V)
    F = []
    for j in range(n_span):
        a, b = j * npnt, (j + 1) * npnt
        for i in range(npnt):
            i1 = (i + 1) % npnt
            F.append([a + i, a + i1, b + i1])
            F.append([a + i, b + i1, b + i])
    for ring, cen_y, flip in ((0, 0, True), (n_span, n_span, False)):
        base = ring * npnt
        cen = V[base:base + npnt].mean(0)
        V = np.vstack([V, cen])
        k = V.shape[0] - 1
        for i in range(npnt):
            i1 = (i + 1) % npnt
            F.append([k, base + i, base + i1] if not flip else [k, base + i1, base + i])
    F = np.asarray(F, int)
    # orient outward (signed volume > 0)
    X = V[F]
    vol = np.einsum("ij,ij->i", X[:, 0], np.cross(X[:, 1], X[:, 2])).sum() / 6.0
    if vol < 0:
        F = F[:, ::-1]
    return V, F


# ----------------------------------------------------------------- public: input
def load_surface(path, length_scale: float = 1.0, size: float | None = None,
                 occ_unit: str | None = "M") -> tuple[np.ndarray, np.ndarray]:
    """Triangulated surface of a geometry file -> ``(vertices (n,3), faces (m,3))``.

    STL/OBJ/PLY/...: read with trimesh (vertices merged).  STEP/IGES/BREP:
    2-D gmsh mesh of all surfaces at ``size`` (default: 1/20 of the
    bounding-box diagonal).  Coordinates are multiplied by ``length_scale``.
    """
    ext = _ext(path)
    if ext in CAD_EXTENSIONS:
        with _Gmsh("surface") as gmsh:
            _load_into_gmsh(gmsh, path, length_scale, occ_unit, False, 40.0)
            lo = np.array(gmsh.model.getBoundingBox(-1, -1))
            diag = float(np.linalg.norm(lo[3:] - lo[:3]))
            _set_sizes(gmsh, size or diag / 20.0)
            gmsh.model.mesh.generate(2)
            xyz, lut = _gmsh_nodes(gmsh)
            tris = _gmsh_elements(gmsh, 2, 3, lut)
        return _compress(xyz, tris)
    import trimesh
    tm = trimesh.load(str(path), force="mesh", process=True)
    tm.merge_vertices()
    tm.remove_unreferenced_vertices()
    return np.asarray(tm.vertices, float) * length_scale, np.asarray(tm.faces, int)


def mesh_geometry(source, kind: str = "solid", element: str = "tet10", size: float | None = None,
                  thickness: float | None = None, material="aluminium_7075",
                  length_scale: float = 1.0, remesh: bool | None = None, angle_deg: float = 40.0,
                  occ_unit: str | None = "M", name: str | None = None,
                  verbose: bool = False) -> FEMesh:
    """Finite-element mesh of a geometry.

    ``source``: file path (STL/OBJ/PLY/STEP/IGES/BREP), a ``trimesh.Trimesh``
    or a ``(vertices, faces)`` tuple.  ``kind``: ``"solid"`` (tet4/tet10,
    ``element``) or ``"shell"`` (tri3 DKT/Allman shells, ``thickness``
    required).  ``size``: target element size [m] (default 1/15 of the
    smallest bounding-box edge for solids, 1/30 of the diagonal for shells).
    ``remesh`` (triangulated input only): remesh the surface at ``size``
    (default: True for solids, False for shells).
    """
    kind = kind.lower()
    if kind not in ("solid", "shell"):
        raise ValueError("kind must be 'solid' or 'shell'")
    path, verts, faces = None, None, None
    if isinstance(source, (str, os.PathLike)):
        path = source
        if _ext(path) not in CAD_EXTENSIONS:
            verts, faces = load_surface(path, length_scale)
    elif isinstance(source, tuple):
        verts, faces = np.asarray(source[0], float) * length_scale, np.asarray(source[1], int)
    else:   # trimesh-like
        verts, faces = np.asarray(source.vertices, float) * length_scale, np.asarray(source.faces, int)
    is_cad = path is not None and _ext(path) in CAD_EXTENSIONS
    name = name or (os.path.splitext(os.path.basename(str(path)))[0] if path else "geometry")
    mat = get_material(material)
    if remesh is None:
        remesh = kind == "solid"

    if kind == "shell":
        if thickness is None:
            raise ValueError("kind='shell' needs a thickness")
        if not is_cad and not remesh:
            nodes, tris = verts, faces
        else:
            with _Gmsh(name, int(verbose) > 1) as gmsh:
                _load_into_gmsh(gmsh, path, 1.0 if verts is not None else length_scale, occ_unit,
                                True, angle_deg, verts, faces)
                bb = np.array(gmsh.model.getBoundingBox(-1, -1))
                _set_sizes(gmsh, size or float(np.linalg.norm(bb[3:] - bb[:3])) / 30.0)
                gmsh.model.mesh.generate(2)
                xyz, lut = _gmsh_nodes(gmsh)
                tris = _gmsh_elements(gmsh, 2, 3, lut)
            nodes, tris = _compress(xyz, tris)
        # drop degenerate triangles
        X = nodes[tris]
        area = 0.5 * np.linalg.norm(np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]), axis=1)
        tris = tris[area > 1e-12 * max(area.max(), 1e-300)]
        nodes, tris = _compress(nodes, tris)
        return FEMesh(nodes, name=name).add("tri3", tris, mat, thickness, name="skin")

    element = element.lower()
    if element not in ("tet4", "tet10"):
        raise ValueError("solid element must be 'tet4' or 'tet10'")
    with _Gmsh(name, int(verbose) > 1) as gmsh:
        cad = _load_into_gmsh(gmsh, path, 1.0 if verts is not None else length_scale, occ_unit,
                              bool(remesh), angle_deg, verts, faces)
        bb = np.array(gmsh.model.getBoundingBox(-1, -1))
        ext = bb[3:] - bb[:3]
        h = size or float(np.min(ext[ext > 0])) / 15.0
        if not cad:
            surfs = [s[1] for s in gmsh.model.getEntities(2)]
            loop = gmsh.model.geo.addSurfaceLoop(surfs)
            gmsh.model.geo.addVolume([loop])
            gmsh.model.geo.synchronize()
        elif not gmsh.model.getEntities(3):
            raise ValueError(f"{path}: no solid volume in the CAD file (use kind='shell')")
        _set_sizes(gmsh, h)
        gmsh.model.mesh.generate(3)
        xyz, lut = _gmsh_nodes(gmsh)
        tets = _gmsh_elements(gmsh, 4, 4, lut)
    if tets.size == 0:
        raise RuntimeError("gmsh produced no tetrahedra (is the surface closed?)")
    nodes, tets = _compress(xyz, tets)
    neg = tet_volumes(nodes[tets]) < 0
    tets[neg, 1], tets[neg, 2] = tets[neg, 2].copy(), tets[neg, 1].copy()
    if element == "tet10":
        nodes, tets = tet4_to_tet10(nodes, tets)
        tets = fix_tet10_order(nodes, tets)
    return FEMesh(nodes, name=name).add(element, tets, mat, name="solid")


# ----------------------------------------------------------------- BCs, analysis
def clamp_nodes(nodes: np.ndarray, clamp) -> np.ndarray:
    """Row indices of the nodes selected by ``clamp`` (dict or list of dicts,
    keys as :func:`select_nodes`; a list is OR-ed; an int array = row indices)."""
    if clamp is None:
        return np.zeros(0, int)
    if isinstance(clamp, dict):
        clamp = [clamp]
    if isinstance(clamp, np.ndarray) or (isinstance(clamp, (list, tuple)) and clamp
                                         and not isinstance(clamp[0], dict)):
        return np.unique(np.asarray(clamp, int))
    sel = [select_nodes(nodes, **c) for c in clamp]
    return np.unique(np.concatenate(sel)) if sel else np.zeros(0, int)


def geometry_modal(source, nmodes: int = 10, kind: str = "solid", element: str = "tet10",
                   size: float | None = None, thickness: float | None = None,
                   material="aluminium_7075", clamp=None, components: str = "123456",
                   length_scale: float = 1.0, remesh: bool | None = None, lumped: bool = False,
                   rotations: str = "lsq", skip_rigid: bool = False, damping: float = 0.0,
                   mesh: FEMesh | None = None, verbose: bool = False, **mesh_kw) -> GeometryModal:
    """Geometry file -> mesh -> clamped normal modes.

    Parameters are those of :func:`mesh_geometry` plus ``nmodes``,
    ``clamp`` (see :func:`clamp_nodes`, e.g. ``dict(below=("y", 0.0))``),
    ``components`` (clamped DOFs, "123456" = all), ``lumped`` mass,
    ``rotations`` ("lsq"/"zero" for solid nodes), ``skip_rigid`` (free-free:
    return ``nmodes`` elastic modes) and modal ``damping``.  A prebuilt
    ``mesh`` skips the meshing step.

    Returns :class:`GeometryModal` (``.modal`` ModalModel with mass-normalised
    shapes (Nn,6,nm), ``gen_mass = 1``; ``.mesh``; ``.clamped``; ``.info``
    with node/element/DOF counts and timings).
    """
    t0 = time.perf_counter()
    if mesh is None:
        mesh = mesh_geometry(source, kind, element, size, thickness, material, length_scale,
                             remesh, verbose=verbose, **mesh_kw)
    t1 = time.perf_counter()
    cl = clamp_nodes(mesh.nodes, clamp)
    if clamp is not None and cl.size == 0:
        raise ValueError(f"clamp {clamp!r} selected no nodes")
    fe = build_fe_model(mesh, fixed=dofs_of_nodes(cl, components), lumped=lumped)
    t2 = time.perf_counter()
    modal = solve_modes(fe, nmodes, rotations=rotations, skip_rigid=skip_rigid,
                        damping=damping, mesh=mesh, name=mesh.name)
    t3 = time.perf_counter()
    nel = sum(b.n for b in mesh.blocks)
    ux = np.zeros(fe.M.shape[0])
    ux[0::6] = 1.0
    mass = float(ux @ (fe.M @ ux))
    info = dict(nodes=mesh.nn, elements=nel, kinds=[b.kind for b in mesh.blocks],
                dofs=int(6 * mesh.nn), clamped_nodes=int(cl.size), mass=mass,
                t_mesh=t1 - t0, t_assemble=t2 - t1, t_solve=t3 - t2, t_total=t3 - t0)
    if verbose:
        print(f"[geometry_modal] {mesh.name}: {mesh.nn} nodes, {nel} {info['kinds']} elements, "
              f"mass {info['mass']:.4g} kg, mesh {info['t_mesh']:.2f}s assemble "
              f"{info['t_assemble']:.2f}s solve {info['t_solve']:.2f}s")
        print("  f [Hz] =", np.array2string(modal.freq_hz, precision=3))
    return GeometryModal(modal, mesh, cl, info)


# ----------------------------------------------------------------- plotting
def surface_triangulation(mesh_or_result) -> np.ndarray:
    """Boundary triangles (nt,3) of node rows of a mesh / GeometryModal /
    ModalModel with ``fe.elements['mesh']`` (solid outer faces, shell elements)."""
    m = mesh_or_result
    if isinstance(m, GeometryModal):
        m = m.mesh
    elif isinstance(m, ModalModel):
        m = m.fe.elements["mesh"]
    return m.surface_triangles()


def planform_nodes(nodes: np.ndarray, n1: int = 20, n2: int = 12, axes=(1, 0),
                   exclude=None) -> np.ndarray:
    """A planform-uniform node subset (row indices) for surface splining.

    The (axes[0], axes[1]) bounding rectangle (default span y x chord x) is
    divided into ``n1 x n2`` cells and, per cell, the node closest to the
    cell centre in the plane is kept (so upper/lower-surface twins are not
    both used).  ``exclude``: row indices never selected (e.g. clamped).
    Returns a boolean mask (Nn,).
    """
    nodes = np.asarray(nodes, float)
    a, b = axes
    P = nodes[:, [a, b]]
    lo, hi = P.min(0), P.max(0)
    span = np.maximum(hi - lo, 1e-12)
    ij = np.minimum(((P - lo) / span * [n1, n2]).astype(int), [n1 - 1, n2 - 1])
    cell = ij[:, 0] * n2 + ij[:, 1]
    centre = lo + (ij + 0.5) / [n1, n2] * span
    d = np.linalg.norm((P - centre) / span, axis=1)
    if exclude is not None:
        d[np.asarray(exclude, int)] = np.inf
    order = np.lexsort((d, cell))
    first = np.ones(order.size, bool)
    first[1:] = cell[order][1:] != cell[order][:-1]
    pick = order[first]
    pick = pick[np.isfinite(d[pick])]
    mask = np.zeros(nodes.shape[0], bool)
    mask[pick] = True
    return mask


def _equal_axes(ax, X):
    c = 0.5 * (X.max(0) + X.min(0))
    r = 0.5 * np.ptp(X, axis=0).max()
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)


def plot_mode(modal: ModalModel, mode: int = 0, mesh: FEMesh | None = None, tris=None,
              scale: float | None = None, ax=None, undeformed: bool = True, cmap: str = "viridis",
              view=(25, -60), title: str | None = None):
    """3-D surface plot of mode ``mode`` (0-based) coloured by |u|.

    ``scale``: displacement scale (default: max |u| -> 10 % of model size).
    Returns the matplotlib axes.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    if mesh is None:
        mesh = modal.fe.elements["mesh"]
    tris = surface_triangulation(mesh) if tris is None else tris
    X = np.asarray(modal.nodes, float)
    u = modal.shapes[:, :3, mode]
    mag = np.linalg.norm(u, axis=1)
    size = np.ptp(X, axis=0).max()
    if scale is None:
        scale = 0.1 * size / max(mag.max(), 1e-300)
    Y = X + scale * u
    if ax is None:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection="3d")
    if undeformed:
        ax.add_collection3d(Poly3DCollection(X[tris], facecolor=(1, 1, 1, 0), edgecolor=(0.6, 0.6, 0.6, 0.3),
                                             linewidth=0.15))
    fc = plt.get_cmap(cmap)(mag[tris].mean(1) / max(mag.max(), 1e-300))
    pc = Poly3DCollection(Y[tris], facecolor=fc, edgecolor="k", linewidth=0.05)
    ax.add_collection3d(pc)
    _equal_axes(ax, np.vstack([X, Y]))
    ax.view_init(*view)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    lab = modal.labels[mode] if modal.labels else f"mode {mode + 1}"
    ax.set_title(title or f"{lab}: {modal.freq_hz[mode]:.3f} Hz", fontsize=9)
    return ax


def plot_modes(modal: ModalModel, modes=None, mesh: FEMesh | None = None, path: str | None = None,
               ncol: int = 3, suptitle: str | None = None, **kw):
    """Grid of :func:`plot_mode` panels; saves to ``path`` if given.  Returns the figure."""
    import matplotlib.pyplot as plt
    modes = list(range(min(6, modal.omega.size))) if modes is None else list(modes)
    mesh = modal.fe.elements["mesh"] if mesh is None else mesh
    tris = surface_triangulation(mesh)
    nrow = int(np.ceil(len(modes) / ncol))
    fig = plt.figure(figsize=(4.2 * min(ncol, len(modes)), 3.6 * nrow))
    for i, m in enumerate(modes):
        ax = fig.add_subplot(nrow, min(ncol, len(modes)), i + 1, projection="3d")
        plot_mode(modal, m, mesh=mesh, tris=tris, ax=ax, **kw)
    if suptitle:
        fig.suptitle(suptitle)
    fig.tight_layout()
    if path is not None:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        fig.savefig(path, dpi=130)
    return fig
