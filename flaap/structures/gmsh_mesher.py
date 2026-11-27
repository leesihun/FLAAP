"""gmsh-based mesher: ``.geo`` scripts, ``.msh`` meshes and CAD files -> FE mesh -> modes.

:mod:`flaap.structures.geometry_modal` uses gmsh internally to mesh STL / STEP
geometry with one automatic size. This module gives direct access to gmsh as
the mesher of a FLAAP FE model:

* ``.geo`` (gmsh script): the geometry, sizes, transfinite lines,
  ``Recombine`` and ``Physical`` groups defined in the script are used as
  written. If the script does not mesh by itself (no ``Mesh 2;`` / ``Mesh 3;``),
  the model is meshed in the highest dimension that has entities.
* ``.msh`` (any version gmsh reads: MSH 2.2 / 4.1, ASCII or binary): an
  existing mesh, used as is.
* ``.step`` / ``.stp`` / ``.iges`` / ``.igs`` / ``.brep``: OpenCASCADE import
  (STEP/IGES converted to metres), then meshed.

Mesher options (they apply when gmsh generates the mesh, i.e. not to a
``.msh``): ``size`` (``Mesh.MeshSizeMax``), ``order`` (1 or 2),
``recombine`` (quadrilateral surface elements, ``Mesh.RecombineAll``),
``dim`` (2 = shell, 3 = solid), ``algorithm`` (``Mesh.Algorithm``) and
``options`` (any further ``{"Mesh.Xxx": value}``). The options are applied
after the file is loaded and override the settings the script made.

Element mapping (FLAAP FE, :mod:`flaap.structures.fem`):

=====================  ==================  ==============================
gmsh type              FLAAP element       note
=====================  ==================  ==============================
2 tri3 / 3 quad4       ``tri3`` / ``quad4`` DKT/Allman, MITC4 shells (``thickness``)
4 tet4 / 11 tet10      ``tet4`` / ``tet10`` mid-edge order fixed (gmsh 8/9 swap)
5 hex8                 ``hex8``
9 tri6, 16/10 quad8/9  corner nodes         warning: no quadratic shells
17 hex20               corner nodes (hex8)  warning
6 prism, 7 pyramid     not supported        ValueError
=====================  ==================  ==============================

Only elements of the highest dimension present are used (shells for a 2-D
mesh, solids for a 3-D mesh). When the model has physical groups of that
dimension, only their elements are used, as gmsh saves them. Each group
becomes its own element block, so it can carry its own shell thickness.
Physical groups of any dimension can clamp nodes (``clamp_group``): for
example, a ``Physical Curve("Root")`` on the wing root edge.
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field

import numpy as np

from .fem import FEMesh, get_material
from .geometry_modal import CAD_EXTENSIONS, GeometryModal, _Gmsh, geometry_modal
from .solid import fix_tet10_order, tet_volumes

__all__ = ["GMSH_EXTENSIONS", "GmshMesh", "gmsh_mesh", "gmsh_femesh", "gmsh_modal", "group_nodes"]

GMSH_EXTENSIONS = (".geo", ".msh") + CAD_EXTENSIONS

# gmsh element type -> (FLAAP kind, number of nodes, corner nodes kept)
_TYPES = {2: ("tri3", 3, 3), 3: ("quad4", 4, 4), 4: ("tet4", 4, 4), 5: ("hex8", 8, 8), 11: ("tet10", 10, 10),
          9: ("tri3", 6, 3), 16: ("quad4", 8, 4), 10: ("quad4", 9, 4), 17: ("hex8", 20, 8)}
_UNSUPPORTED = {6: "6-node prism", 7: "5-node pyramid", 13: "18-node prism", 18: "15-node prism",
                14: "14-node pyramid", 19: "13-node pyramid", 12: "27-node hexahedron",
                29: "20-node tetrahedron"}


@dataclass
class GmshMesh:
    """Mesh read from gmsh.

    ``nodes`` (n,3) [m]; ``node_tags`` gmsh node tags; ``blocks`` list of
    ``(kind, conn, group)`` with ``conn`` row indices into ``nodes``;
    ``groups`` ``{name: row indices}`` of every physical group (all
    dimensions); ``dim`` the element dimension used; ``info`` counts.
    """

    nodes: np.ndarray
    node_tags: np.ndarray
    blocks: list
    groups: dict = field(default_factory=dict)
    dim: int = 2
    name: str = "gmsh"
    info: dict = field(default_factory=dict)

    @property
    def nelem(self) -> int:
        return int(sum(len(c) for _, c, _ in self.blocks))

    def kinds(self) -> list[str]:
        return sorted({k for k, _, _ in self.blocks})


def _ext(path) -> str:
    return os.path.splitext(str(path))[1].lower()


def _phys_name(gmsh, dim, tag):
    name = gmsh.model.getPhysicalName(dim, tag)
    return name if name else f"{['point', 'curve', 'surface', 'volume'][dim]}{tag}"


def _has_mesh(gmsh) -> bool:
    tags, _, _ = gmsh.model.mesh.getNodes()
    if len(tags) == 0:
        return False
    for d in (3, 2):
        types, _, _ = gmsh.model.mesh.getElements(d)
        if len(types):
            return True
    return False


def gmsh_mesh(path, *, size: float | None = None, size_min: float | None = None, order: int | None = None,
              recombine: bool | None = None, dim: int | None = None, algorithm: int | None = None,
              options: dict | None = None, length_scale: float = 1.0, occ_unit: str | None = "M",
              save: str | None = None, verbose: bool = False) -> GmshMesh:
    """Load / mesh ``path`` (``.geo``, ``.msh``, STEP/IGES/BREP) with gmsh.

    ``length_scale`` multiplies the node coordinates of a ``.geo`` / ``.msh``
    / BREP (STEP and IGES are converted to metres by OpenCASCADE). ``save``:
    also write the mesh to this ``.msh`` path (MSH 4.1 unless the extension
    says otherwise).
    """
    path = str(path)
    ext = _ext(path)
    if ext not in GMSH_EXTENSIONS:
        raise ValueError(f"{path}: gmsh input must be one of {', '.join(GMSH_EXTENSIONS)}")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    name = os.path.splitext(os.path.basename(path))[0]
    with _Gmsh(name, verbose) as gmsh:
        if ext in CAD_EXTENSIONS:
            if occ_unit and ext not in (".brep", ".brp"):
                gmsh.option.setString("Geometry.OCCTargetUnit", occ_unit)
            gmsh.model.occ.importShapes(path)
            gmsh.model.occ.synchronize()
            scale_after = 1.0 if ext not in (".brep", ".brp") else float(length_scale)
        else:
            gmsh.open(path)                   # .geo runs the script; .msh loads the mesh
            scale_after = float(length_scale)
        generated = False
        if not _has_mesh(gmsh):
            if ext == ".msh":
                raise ValueError(f"{path}: the .msh file has no 2-D or 3-D elements")
            if size is not None:
                gmsh.option.setNumber("Mesh.MeshSizeMax", float(size))
                gmsh.option.setNumber("Mesh.MeshSizeMin", float(size_min if size_min else 0.1 * size))
            if recombine:
                gmsh.option.setNumber("Mesh.RecombineAll", 1)
                if algorithm is None:
                    gmsh.option.setNumber("Mesh.Algorithm", 8)       # frontal-Delaunay for quads
            if algorithm is not None:
                gmsh.option.setNumber("Mesh.Algorithm", int(algorithm))
            for k, val in (options or {}).items():
                if isinstance(val, str):
                    gmsh.option.setString(k, val)
                else:
                    gmsh.option.setNumber(k, float(val))
            top = max([d for d in (3, 2, 1) if gmsh.model.getEntities(d)], default=0)
            d = int(dim) if dim else top
            if d < 2 or top < 2:
                raise ValueError(f"{path}: no surface or volume to mesh (FLAAP FE needs shells or solids)")
            gmsh.model.mesh.generate(d)
            if order and int(order) > 1:
                gmsh.model.mesh.setOrder(int(order))
            generated = True
        elif order and int(order) > 1:
            gmsh.model.mesh.setOrder(int(order))
        if save:
            gmsh.write(str(save))
        out = _extract(gmsh, path, dim, scale_after)
    out.name = name
    out.info["generated"] = generated
    if verbose:
        print(f"[gmsh] {name}: {len(out.nodes)} nodes, {out.nelem} elements {out.kinds()}, "
              f"groups {sorted(out.groups)}")
    return out


def _extract(gmsh, path, dim, scale) -> GmshMesh:
    tags, xyz, _ = gmsh.model.mesh.getNodes()
    tags = np.asarray(tags, np.int64)
    xyz = np.asarray(xyz, float).reshape(-1, 3) * scale
    lut = np.full(int(tags.max()) + 1, -1, np.int64)
    lut[tags] = np.arange(tags.size)

    have = [d for d in (3, 2) if len(gmsh.model.mesh.getElements(d)[0])]
    if not have:
        raise ValueError(f"{path}: no 2-D or 3-D elements in the mesh")
    d = int(dim) if dim and int(dim) in have else have[0]

    # physical groups (all dimensions) -> node rows
    groups, ent_group = {}, {}
    for gd, gt in gmsh.model.getPhysicalGroups():
        nm = _phys_name(gmsh, gd, gt)
        ntags, _ = gmsh.model.mesh.getNodesForPhysicalGroup(gd, gt)
        rows = lut[np.asarray(ntags, np.int64)] if len(ntags) else np.zeros(0, np.int64)
        groups[nm] = np.unique(np.concatenate([groups.get(nm, np.zeros(0, np.int64)), rows])).astype(int)
        if gd == d:
            for et in gmsh.model.getEntitiesForPhysicalGroup(gd, gt):
                ent_group.setdefault(int(et), nm)

    entities = [t for _, t in gmsh.model.getEntities(d)]
    if ent_group:                               # gmsh semantics: only physical elements
        entities = [t for t in entities if t in ent_group]
    acc: dict[tuple[str, str], list] = {}
    reduced = set()
    for et in entities:
        types, _, ntags = gmsh.model.mesh.getElements(d, et)
        grp = ent_group.get(et, "all")
        for ty, nt in zip(types, ntags):
            ty = int(ty)
            if ty in _UNSUPPORTED:
                raise ValueError(f"{path}: {_UNSUPPORTED[ty]} elements are not supported by the FLAAP FE "
                                 "(use tetrahedra / hexahedra, or do not recombine the volume mesh)")
            if ty not in _TYPES:
                name = gmsh.model.mesh.getElementProperties(ty)[0]
                raise ValueError(f"{path}: gmsh element type {ty} ({name}) is not supported")
            kind, npe, keep = _TYPES[ty]
            conn = lut[np.asarray(nt, np.int64).reshape(-1, npe)][:, :keep]
            if keep < npe:
                reduced.add(gmsh.model.mesh.getElementProperties(ty)[0])
            acc.setdefault((kind, grp), []).append(conn)
    if reduced:
        warnings.warn(f"{os.path.basename(str(path))}: {', '.join(sorted(reduced))} elements reduced to their "
                      "corner nodes (FLAAP has no quadratic shells / hex20)", stacklevel=3)

    # keep only the nodes the elements use
    blocks = [(k, np.vstack(c), g) for (k, g), c in acc.items()]
    used = np.unique(np.concatenate([c.ravel() for _, c, _ in blocks]))
    new = np.full(len(xyz), -1, np.int64)
    new[used] = np.arange(used.size)
    nodes = xyz[used]
    out_blocks = []
    for kind, conn, grp in blocks:
        conn = new[conn].astype(int)
        if kind == "tet10":
            conn = fix_tet10_order(nodes, conn)
        elif kind == "tet4":
            neg = tet_volumes(nodes[conn]) < 0
            conn[neg, 1], conn[neg, 2] = conn[neg, 2].copy(), conn[neg, 1].copy()
        out_blocks.append((kind, conn, grp))
    groups = {k: np.unique(new[r][new[r] >= 0]).astype(int) for k, r in groups.items()}
    info = {"nodes": int(len(nodes)), "elements": int(sum(len(c) for _, c, _ in out_blocks)),
            "dim": d, "kinds": sorted({k for k, _, _ in out_blocks}),
            "groups": {k: int(len(v)) for k, v in groups.items()}}
    return GmshMesh(nodes, tags[used], out_blocks, groups, d, info=info)


def group_nodes(gm: GmshMesh, names) -> np.ndarray:
    """Row indices of the nodes of physical group(s) ``names`` (str or list)."""
    if isinstance(names, str):
        names = [names]
    missing = [n for n in names if n not in gm.groups]
    if missing:
        raise ValueError(f"physical group(s) {missing} not found; the mesh has {sorted(gm.groups) or 'none'}")
    return np.unique(np.concatenate([gm.groups[n] for n in names])).astype(int)


def gmsh_femesh(gm: GmshMesh, material="aluminium_7075", thickness: float | None = None,
                group_thickness: dict | None = None, group_material: dict | None = None) -> FEMesh:
    """:class:`FEMesh` of a :class:`GmshMesh`: one block per (element kind, physical group).

    ``thickness``: shell thickness [m] (``group_thickness`` overrides it per group);
    ``group_material``: per-group material (preset name, dict or Material).
    """
    group_thickness = group_thickness or {}
    group_material = group_material or {}
    for g in list(group_thickness) + list(group_material):
        if g not in {grp for _, _, grp in gm.blocks}:
            raise ValueError(f"group '{g}' has no {gm.dim}-D elements "
                             f"(element groups: {sorted({grp for _, _, grp in gm.blocks})})")
    fe = FEMesh(gm.nodes, node_ids=gm.node_tags, name=gm.name)
    for kind, conn, grp in gm.blocks:
        mat = get_material(group_material.get(grp, material))
        t = group_thickness.get(grp, thickness)
        if kind in ("tri3", "quad4"):
            if t is None:
                raise ValueError(f"shell elements in group '{grp}' need a thickness")
            if kind == "tri3":                                 # drop degenerate triangles
                X = gm.nodes[conn]
                a = 0.5 * np.linalg.norm(np.cross(X[:, 1] - X[:, 0], X[:, 2] - X[:, 0]), axis=1)
                conn = conn[a > 1e-12 * max(a.max(), 1e-300)]
            fe.add(kind, conn, mat, float(t), name=grp)
        else:
            fe.add(kind, conn, mat, name=grp)
    return fe


def gmsh_modal(path, nmodes: int = 10, *, material="aluminium_7075", thickness: float | None = None,
               group_thickness: dict | None = None, group_material: dict | None = None,
               clamp=None, clamp_group=None, components: str = "123456", damping: float = 0.0,
               length_scale: float = 1.0, lumped: bool = False, skip_rigid: bool = False,
               verbose: bool = False, **mesh_kw) -> GeometryModal:
    """gmsh file -> FE mesh -> clamped normal modes.

    ``clamp_group``: physical group name(s) whose nodes are clamped (OR-ed with
    ``clamp``, see :func:`flaap.structures.geometry_modal.clamp_nodes`).
    ``mesh_kw`` go to :func:`gmsh_mesh` (``size``, ``order``, ``recombine``,
    ``dim``, ``algorithm``, ``options``, ``save``). Returns
    :class:`~flaap.structures.geometry_modal.GeometryModal`; ``info`` also
    holds the gmsh counts and physical groups.
    """
    from .geometry_modal import clamp_nodes
    gm = gmsh_mesh(path, length_scale=length_scale, verbose=verbose, **mesh_kw)
    fe = gmsh_femesh(gm, material, thickness, group_thickness, group_material)
    rows = []
    if clamp_group:
        rows.append(group_nodes(gm, clamp_group))
    if clamp is not None:
        rows.append(clamp_nodes(gm.nodes, clamp))
    cl = np.unique(np.concatenate(rows)).astype(int) if rows else None
    res = geometry_modal(None, nmodes, mesh=fe, clamp=cl, components=components, damping=damping,
                         lumped=lumped, skip_rigid=skip_rigid, verbose=verbose)
    res.info["gmsh"] = gm.info
    res.info["gmsh_mesh"] = gm
    return res
