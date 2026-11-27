"""OpenRadioss starter-deck reader (``*_0000.rad``) -> FE model -> normal modes.

OpenRadioss is an explicit-dynamics code and writes no modal results, so
FLAAP reads the *model* (nodes, elements, properties, materials, boundary
conditions) from the starter input and computes the normal modes with
its own elements (:mod:`flaap.structures.fem`: MITC4/DKT flat shells,
tet4/tet10/hex8 solids).

Supported keywords
------------------
=====================================  ==========================================
``/BEGIN``                             run name, version, input / work units
``/NODE``                              ``node_ID  X  Y  Z`` (I10, 3F20)
``/SHELL/part_ID``                     4-node shells (n3 == n4 -> triangle);
                                       optional element thickness (cols 81-100)
``/SH3N/part_ID``                      3-node shells
``/BRICK/part_ID``                     8-node bricks (degenerate bricks with
                                       repeated nodes are rejected)
``/TETRA4/part_ID``                    4-node tetrahedra
``/TETRA10/part_ID``                   10-node tetrahedra (corners, then mid
                                       nodes 12 23 31 14 24 34; continuation
                                       lines with blank first field)
``/PART/part_ID``                      title, ``prop_ID  mat_ID  [subset_ID]``
``/PROP/SHELL`` (``/PROP/TYPE1``)      thickness ``Thick`` (3rd data line)
``/PROP/SOLID`` (``TYPE14``),          accepted (no data needed)
``/PROP/TETRA4``, ``/PROP/TETRA10``
``/MAT/LAW1`` (``/MAT/ELAST``)         rho, E, nu
``/MAT/LAW2`` (``PLAS_JOHNS``,         rho, E, nu (elastic part only)
``PLAS_ZERIL``), ``/MAT/LAW36``
(``PLAS_TAB``), ``LAW27``, ``LAW44``
``/BCS/bcs_ID``                        ``Tra rot`` digits, ``skew_ID``, ``grnod_ID``
``/GRNOD/NODE``, ``/GRNOD/PART``,      node groups
``/GRNOD/GRNOD``, ``/GRNOD/BOX``
``/BOX/RECTA``                         axis-aligned box (two corner points)
``/SKEW/FIX``                          fixed skew (origin, X', Y'); BCs in a
                                       skew frame are applied exactly by a
                                       nodal rotation of ``K`` and ``M``
=====================================  ==========================================

Format rules: ``#`` lines are comments, keyword lines start with ``/``,
the first line after most option keywords is a title.  Data are in fixed
fields of 10 characters (integers) and 20 characters (reals); a line
that cannot be read with fixed fields falls back to whitespace splitting.
Unknown keywords are skipped (listed in ``RadiossDeck.skipped``).

Units: the deck values are in the ``/BEGIN`` *input* units (3rd line:
mass, length, time, e.g. ``kg mm ms`` or ``Mg mm s``) and are converted to
SI on reading.  Without ``/BEGIN`` (or units) SI is assumed unless
``units=(mass, length, time)`` is given.

Elastic properties of plasticity laws are used as is (modal analysis is
linear); composite shell properties are not supported.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field

import numpy as np

from ..structures.fem import ElementBlock, FEMesh, Material, build_fe_model, solve_modes
from ..structures.modal import ModalModel
from ..structures.solid import fix_tet10_order, tet_volumes

__all__ = ["RadiossDeck", "read_radioss", "radioss_mesh", "radioss_modal", "UNIT_MASS",
           "UNIT_LENGTH", "UNIT_TIME"]

UNIT_MASS = {"kg": 1.0, "g": 1e-3, "mg": 1e-6, "Mg": 1e3, "t": 1e3, "ton": 1e3, "tonne": 1e3,
             "lb": 0.45359237, "lbm": 0.45359237, "slug": 14.5939029, "ug": 1e-9}
UNIT_LENGTH = {"m": 1.0, "cm": 1e-2, "mm": 1e-3, "um": 1e-6, "mum": 1e-6, "km": 1e3, "dm": 0.1,
               "in": 0.0254, "inch": 0.0254, "ft": 0.3048, "ft.": 0.3048}
UNIT_TIME = {"s": 1.0, "ms": 1e-3, "us": 1e-6, "mus": 1e-6, "ns": 1e-9, "min": 60.0, "h": 3600.0}

_TITLED = ("PART", "PROP", "MAT", "BCS", "GRNOD", "SKEW", "BOX")


def _unit(table, s, what):
    s = s.strip()
    if not s:
        return 1.0
    if s in table:
        return table[s]
    low = {k.lower(): v for k, v in table.items() if k != "Mg"}
    if s.lower() in low and s != "MG":
        return low[s.lower()]
    if s.upper() == "MG":
        return 1e3 if s in ("MG", "Mg") else 1e-6
    raise ValueError(f"unknown {what} unit '{s}'")


def _fnum(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    s = s.replace("D", "E").replace("d", "e")
    try:
        return float(s)
    except ValueError:
        m = re.match(r"^([+-]?\d*\.?\d+|[+-]?\d+\.)([+-]\d+)$", s)   # 1.5-3 -> 1.5e-3
        if m:
            return float(m.group(1) + "e" + m.group(2))
        raise


def _fields(line: str, widths) -> list[str]:
    out, p = [], 0
    for w in widths:
        out.append(line[p:p + w])
        p += w
    return out


def _ints(line: str, n: int | None = None, width: int = 10) -> list[int]:
    """Fixed 10-char integer fields (blank -> 0); whitespace fallback."""
    line = line.rstrip()
    try:
        k = (len(line) + width - 1) // width if n is None else n
        vals = [int(f) if f.strip() else 0 for f in _fields(line, [width] * k)]
    except ValueError:
        vals = [int(float(t)) for t in line.replace(",", " ").split()]
    return vals


def _reals(line: str, n: int, first_int: int = 0) -> list[float]:
    """``first_int`` 10-char fields then 20-char real fields; fallback split."""
    line = line.rstrip()
    try:
        w = [10] * first_int + [20] * n
        f = _fields(line, w)
        return [_fnum(x) for x in f]
    except ValueError:
        t = line.replace(",", " ").split()
        return [_fnum(x) for x in t] + [0.0] * max(0, first_int + n - len(t))


@dataclass
class RadiossDeck:
    """Parsed starter deck (all values in SI)."""

    path: str = ""
    title: str = ""
    units: tuple = ("kg", "m", "s")
    scale: dict = field(default_factory=lambda: {"mass": 1.0, "length": 1.0, "time": 1.0})
    node_ids: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    nodes: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    elements: dict = field(default_factory=dict)    # (kind, part) -> (elem ids, node ids (ne, k), thick)
    parts: dict = field(default_factory=dict)       # part -> {"title", "prop", "mat"}
    props: dict = field(default_factory=dict)       # prop -> {"type", "thickness", "title"}
    mats: dict = field(default_factory=dict)        # mat -> Material (SI)
    bcs: list = field(default_factory=list)         # {"id", "tra", "rot", "skew", "grnod"}
    grnod: dict = field(default_factory=dict)       # id -> node ids (array)
    skews: dict = field(default_factory=dict)       # id -> 3x3 rows X', Y', Z'
    boxes: dict = field(default_factory=dict)
    skipped: list = field(default_factory=list)

    def summary(self) -> str:
        ne = {k: v[0].size for k, v in self.elements.items()}
        return (f"{self.title}: {self.nodes.shape[0]} nodes, elements {ne}, parts {sorted(self.parts)}, "
                f"materials {sorted(self.mats)}, {len(self.bcs)} BCS, units {self.units}")


def _blocks(lines):
    """Split into (keyword line, data lines) with comments removed."""
    cur, data = None, []
    for ln in lines:
        ln = ln.rstrip("\n\r")
        if ln.startswith("#") or ln.startswith("$"):
            continue
        if ln.startswith("/"):
            if cur is not None:
                yield cur, data
            cur, data = ln.strip(), []
            if cur.upper().startswith("/END") and not cur.upper().startswith("/END_"):
                cur = None
                break
        elif cur is not None:
            data.append(ln.expandtabs(8))
    if cur is not None:
        yield cur, data


def read_radioss(path: str, units=None) -> RadiossDeck:
    """Parse an OpenRadioss starter deck into a :class:`RadiossDeck` (SI)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    deck = RadiossDeck(path=path)
    blocks = list(_blocks(lines))
    # units first
    mu, lu, tu = units if units is not None else ("kg", "m", "s")
    for kw, data in blocks:
        if kw.upper().startswith("/BEGIN"):
            d = [x for x in data if x.strip()]
            deck.title = d[0].strip() if d else ""
            if units is None and len(d) >= 3:
                f = _fields(d[2].ljust(60), [20, 20, 20])
                mu, lu, tu = (f[0].strip() or "kg", f[1].strip() or "m", f[2].strip() or "s")
    deck.units = (mu, lu, tu)
    sm, sl, st = _unit(UNIT_MASS, mu, "mass"), _unit(UNIT_LENGTH, lu, "length"), _unit(UNIT_TIME, tu, "time")
    deck.scale = {"mass": sm, "length": sl, "time": st}
    s_E = sm / (sl * st * st)
    s_rho = sm / sl ** 3

    nid, xyz = [], []
    els: dict = {}
    grnod_raw = {}
    for kw, data in blocks:
        parts = [p for p in kw.split("/") if p]
        key = parts[0].upper()
        sub = parts[1].upper() if len(parts) > 1 else ""
        data_nb = [d for d in data if d.strip()]
        if key == "BEGIN":
            continue
        if key == "NODE":
            for d in data_nb:
                v = _reals(d, 3, first_int=1)
                nid.append(int(v[0]))
                xyz.append(v[1:4])
            continue
        if key in ("SHELL", "SH3N", "BRICK", "TETRA4", "TETRA10"):
            part = int(parts[1]) if len(parts) > 1 else 0
            kind = {"SHELL": "shell", "SH3N": "tri3", "BRICK": "hex8", "TETRA4": "tet4",
                    "TETRA10": "tet10"}[key]
            need = {"shell": 4, "tri3": 3, "hex8": 8, "tet4": 4, "tet10": 10}[kind]
            i = 0
            while i < len(data_nb):
                d = data_nb[i]
                v = _ints(d[:10 * (need + 1)] if kind != "shell" else d[:50], None)
                eid, nodes_ = v[0], [x for x in v[1:] if x != 0][:need]
                while len(nodes_) < need and i + 1 < len(data_nb) and not data_nb[i + 1][:10].strip():
                    i += 1
                    nodes_ += [x for x in _ints(data_nb[i])[1:] if x != 0]
                thick = 0.0
                if kind == "shell" and len(d.rstrip()) > 80:
                    try:
                        thick = _fnum(d[80:100]) * sl
                    except ValueError:
                        thick = 0.0
                if kind == "shell":
                    k2 = "quad4" if len(nodes_) == 4 and nodes_[2] != nodes_[3] else "tri3"
                    nodes_ = nodes_[:4] if k2 == "quad4" else nodes_[:3]
                else:
                    k2 = kind
                    if len(nodes_) < need:
                        raise ValueError(f"{kw}: element {eid} has {len(nodes_)} nodes, {need} expected")
                if k2 == "hex8" and len(set(nodes_)) < 8:
                    raise NotImplementedError(f"degenerate /BRICK {eid} (use /TETRA4 or /PENTA6)")
                els.setdefault((k2, part), []).append((eid, nodes_, thick))
                i += 1
            continue
        if key in _TITLED:
            title = data[0].strip() if data else ""
            body = [d for d in data[1:] if d.strip()]
        if key == "PART":
            v = _ints(body[0], 3)
            deck.parts[int(parts[1])] = {"title": title, "prop": v[0], "mat": v[1]}
        elif key == "PROP":
            pid = int(parts[2]) if len(parts) > 2 else int(parts[1])   # /PROP/TYPE/id[/unit]
            ptype = sub
            entry = {"type": ptype, "title": title, "thickness": None}
            if ptype in ("SHELL", "TYPE1"):
                if len(body) >= 3:
                    v = _reals(body[2], 2, first_int=2)
                    entry["thickness"] = v[2] * sl
            elif ptype in ("SH_SANDW", "TYPE11", "SH_COMP", "TYPE9", "PCOMPP", "TYPE51"):
                warnings.warn(f"/PROP/{ptype}/{pid}: composite shells not supported "
                              "(element thickness required)")
            deck.props[pid] = entry
        elif key == "MAT":
            mid = int(parts[2]) if len(parts) > 2 else 0
            law = sub
            elastic = ("LAW1", "ELAST", "LAW2", "PLAS_JOHNS", "PLAS_ZERIL", "LAW36", "PLAS_TAB",
                       "LAW27", "PLAS_BRIT", "LAW44", "COWPER", "LAW22", "DAMA", "LAW23")
            if law not in elastic:
                warnings.warn(f"/MAT/{law}/{mid}: unsupported law, trying rho / E / nu layout")
            rho = _reals(body[0], 1)[0] * s_rho
            v = _reals(body[1], 2)
            deck.mats[mid] = Material(v[0] * s_E, v[1], rho, title or f"mat{mid}")
        elif key == "BCS":
            bid = int(parts[1])
            ln = body[0]
            code = ln[:10]
            grp = re.findall(r"\d+", code)
            if len(grp) >= 2:
                tra, rot = grp[0].rjust(3, "0")[-3:], grp[1].rjust(3, "0")[-3:]
            elif grp and len(grp[0]) == 6:
                tra, rot = grp[0][:3], grp[0][3:]
            elif grp:
                tra, rot = grp[0].rjust(3, "0")[-3:], "000"
            else:
                tra, rot = "000", "000"
            rest = _ints(ln[10:].ljust(20), 2)
            deck.bcs.append({"id": bid, "title": title, "tra": tra, "rot": rot,
                             "skew": rest[0], "grnod": rest[1]})
        elif key == "GRNOD":
            gid = int(parts[2])
            vals = []
            for d in body:
                vals += [x for x in _ints(d) if x != 0]
            grnod_raw[gid] = (sub, vals)
        elif key == "SKEW":
            sid = int(parts[2])
            if sub not in ("FIX",):
                warnings.warn(f"/SKEW/{sub}/{sid}: only /SKEW/FIX supported; ignored")
                continue
            vals = []
            for d in body:
                vals += _reals(d, 3)
            X1 = np.array(vals[3:6], float)
            Y1 = np.array(vals[6:9], float)
            ex = X1 / np.linalg.norm(X1)
            ez = np.cross(ex, Y1)
            ez /= np.linalg.norm(ez)
            ey = np.cross(ez, ex)
            deck.skews[sid] = np.stack([ex, ey, ez])
        elif key == "BOX":
            bxid = int(parts[2])
            vals = []
            for d in body:
                vals += _reals(d, 3, first_int=0) if len(d.rstrip()) > 10 else []
            p = np.array(vals[-6:], float) * sl
            deck.boxes[bxid] = (np.minimum(p[:3], p[3:]), np.maximum(p[:3], p[3:]))
        elif key in ("NODE",):
            pass
        else:
            deck.skipped.append(kw)

    if not nid:
        raise ValueError(f"{path}: no /NODE data")
    deck.node_ids = np.array(nid, int)
    deck.nodes = np.array(xyz, float) * sl
    for k, lst in els.items():
        ids = np.array([e[0] for e in lst], int)
        conn = np.array([e[1] for e in lst], int)
        th = np.array([e[2] for e in lst], float)
        deck.elements[k] = (ids, conn, th)
    # node groups
    part_nodes = lambda p: np.unique(np.concatenate(  # noqa: E731
        [v[1].ravel() for k, v in deck.elements.items() if k[1] == p] or [np.zeros(0, int)]))

    def resolve(gid, depth=0):
        if gid in deck.grnod:
            return deck.grnod[gid]
        if gid not in grnod_raw or depth > 20:
            raise KeyError(f"/GRNOD {gid} not defined")
        typ, vals = grnod_raw[gid]
        if typ == "NODE":
            out = np.array(vals, int)
        elif typ == "PART":
            out = np.unique(np.concatenate([part_nodes(p) for p in vals]))
        elif typ == "GRNOD":
            out = np.unique(np.concatenate([resolve(g, depth + 1) for g in vals]))
        elif typ == "BOX":
            out = []
            for b in vals:
                lo, hi = deck.boxes[b]
                tol = 1e-9 * max(1.0, float(np.abs(deck.nodes).max()))
                m = np.all((deck.nodes >= lo - tol) & (deck.nodes <= hi + tol), axis=1)
                out.append(deck.node_ids[m])
            out = np.unique(np.concatenate(out))
        else:
            raise NotImplementedError(f"/GRNOD/{typ}")
        deck.grnod[gid] = np.asarray(out, int)
        return deck.grnod[gid]

    for g in grnod_raw:
        try:
            resolve(g)
        except (KeyError, NotImplementedError) as e:
            warnings.warn(str(e))
    return deck


def radioss_mesh(deck: RadiossDeck, tet10_fix: bool = True) -> FEMesh:
    """FE mesh (:class:`FEMesh`) of the deck: one block per (element type, part)."""
    lut = {int(g): i for i, g in enumerate(deck.node_ids)}
    mesh = FEMesh(deck.nodes.copy(), node_ids=deck.node_ids.copy(), name=deck.title or "radioss")
    for (kind, part), (ids, conn_ids, th) in sorted(deck.elements.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        conn = np.vectorize(lut.__getitem__)(conn_ids) if conn_ids.size else conn_ids
        conn = np.asarray(conn, int).reshape(conn_ids.shape)
        p = deck.parts.get(part)
        if p is None:
            raise KeyError(f"elements of undefined /PART {part}")
        if p["mat"] not in deck.mats:
            raise KeyError(f"/PART {part}: material {p['mat']} not defined or unsupported")
        mat = deck.mats[p["mat"]]
        thick = None
        if kind in ("quad4", "tri3"):
            tp = deck.props.get(p["prop"], {}).get("thickness")
            thick = np.where(th > 0, th, tp if tp is not None else 0.0)
            if np.any(thick <= 0):
                raise ValueError(f"/PART {part}: no shell thickness (prop {p['prop']})")
        elif kind == "tet10" and tet10_fix:
            conn = fix_tet10_order(mesh.nodes, conn)
        elif kind == "tet4":
            neg = tet_volumes(mesh.nodes[conn]) < 0
            conn[neg, 1], conn[neg, 2] = conn[neg, 2].copy(), conn[neg, 1].copy()
        elif kind == "hex8":
            X = mesh.nodes[conn]
            v = np.einsum("ei,ei->e", np.cross(X[:, 1] - X[:, 0], X[:, 3] - X[:, 0]), X[:, 4] - X[:, 0])
            neg = v < 0
            if np.any(neg):
                conn[neg] = conn[neg][:, [0, 3, 2, 1, 4, 7, 6, 5]]
        mesh.blocks.append(ElementBlock(kind, conn, mat, thick, name=p["title"] or f"part{part}"))
    return mesh


def _bc_dofs(deck: RadiossDeck, mesh: FEMesh):
    lut = {int(g): i for i, g in enumerate(deck.node_ids)}
    fixed, skew_nodes = [], {}
    for bc in deck.bcs:
        if bc["grnod"] not in deck.grnod:
            raise KeyError(f"/BCS/{bc['id']}: /GRNOD {bc['grnod']} not defined")
        rows = np.array([lut[int(g)] for g in deck.grnod[bc["grnod"]] if int(g) in lut], int)
        comp = [i for i, c in enumerate(bc["tra"] + bc["rot"]) if c == "1"]
        if bc["skew"] and len(comp) not in (0, 6):
            if bc["skew"] not in deck.skews:
                raise KeyError(f"/BCS/{bc['id']}: /SKEW {bc['skew']} not defined")
            for r in rows:
                skew_nodes[r] = deck.skews[bc["skew"]]
        fixed.append((6 * rows[:, None] + np.array(comp, int)[None, :]).ravel())
    fixed = np.unique(np.concatenate(fixed)) if fixed else np.zeros(0, int)
    return fixed, skew_nodes


def radioss_modal(path: str, nmodes: int = 10, units=None, lumped: bool = False,
                  rotations: str = "lsq", skip_rigid: bool = False, sigma: float = -1.0,
                  damping: float = 0.0, deck: RadiossDeck | None = None) -> ModalModel:
    """Normal modes of an OpenRadioss starter deck.

    Reads ``path`` (:func:`read_radioss`), builds the FE model
    (:func:`radioss_mesh`), applies the ``/BCS`` constraints (DOFs in a
    ``/SKEW`` frame by nodal transformation) and solves
    (:func:`flaap.structures.fem.solve_modes`).  Returns a mass-normalized
    :class:`ModalModel` in SI units with ``node_ids`` attribute (Radioss
    node ids) and ``fe.elements["mesh"]`` (the :class:`FEMesh`).
    """
    import scipy.sparse as sp

    deck = deck or read_radioss(path, units=units)
    mesh = radioss_mesh(deck)
    fixed, skew_nodes = _bc_dofs(deck, mesh)
    fe = build_fe_model(mesh, fixed=fixed, lumped=lumped)
    T = None
    if skew_nodes:
        n = 6 * mesh.nn
        rows, cols, vals = list(range(n)), list(range(n)), [1.0] * n
        diag = np.ones(n, bool)
        for r, R in skew_nodes.items():
            for blk in (0, 3):
                base = 6 * r + blk
                diag[base:base + 3] = False
                for a in range(3):
                    for b in range(3):
                        rows.append(base + a)
                        cols.append(base + b)
                        vals.append(R[b, a])       # u_global = R^T u_skew
        keep = diag[np.array(rows)] | (np.arange(len(rows)) >= n)
        T = sp.coo_matrix((np.array(vals)[keep], (np.array(rows)[keep], np.array(cols)[keep])),
                          shape=(n, n)).tocsr()
        fe.K = (T.T @ fe.K @ T).tocsr()
        fe.M = (T.T @ fe.M @ T).tocsr()
    mm = solve_modes(fe, nmodes, sigma=sigma, rotations=rotations, skip_rigid=skip_rigid,
                     damping=damping, name=deck.title or "radioss", mesh=mesh)
    if T is not None:
        U = T @ mm.shapes.reshape(-1, mm.nmodes)
        mm.shapes = U.reshape(mm.shapes.shape)
        if rotations == "lsq":      # re-estimate solid-node rotations in the global frame
            from ..structures.solid import rotations_from_translations
            shell = np.zeros(mesh.nn, bool)
            for b in mesh.blocks:
                if b.kind in ("quad4", "tri3"):
                    shell[b.conn.ravel()] = True
            miss = ~shell
            if miss.any():
                rot = rotations_from_translations(mesh.nodes, mesh.conns(), mm.shapes[:, :3, :], miss)
                mm.shapes[miss, 3:, :] = rot[miss]
    mm.node_ids = deck.node_ids.copy()
    return mm
