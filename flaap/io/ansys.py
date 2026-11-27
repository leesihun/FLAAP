"""ANSYS (MAPDL / Workbench Mechanical) modal results -> :class:`ModalModel`.

The structural side of a flutter analysis is often done in ANSYS.  This
module reads the *text* output that a modal run produces, so that no ANSYS
licence or binary-format library is needed, plus (optionally) the binary
``.rst`` file through ``ansys.mapdl.reader`` when that package happens to be
installed.

Supported layouts
-----------------
MAPDL listing (batch ``.out`` file, ``/OUTPUT`` file or a copy of the
Output Window); several of the blocks below may be in one file or in
separate files.  MAPDL page headers (form feed ``1`` in column 1, the
``***** ANSYS - ENGINEERING ANALYSIS SYSTEM  RELEASE ... *****`` banner, the
``VERSION=... CP=`` line and the title line) are removed before parsing,
so tables split over several pages are merged.

* frequencies, :func:`read_frequencies`

  - ``SET,LIST`` -> ``*****  INDEX OF DATA SETS ON RESULTS FILE  *****`` with
    columns ``SET  TIME/FREQ  LOAD STEP  SUBSTEP  CUMULATIVE`` (preferred);
  - the modal solver summary ``*** FREQUENCIES FROM BLOCK LANCZOS
    ITERATION ***`` / ``MODE    FREQUENCY (HERTZ)`` table, and any other
    ``MODE  FREQUENCY ...`` table such as the participation-factor or
    modal-mass summaries (first two columns, first occurrence per mode);
  - a Workbench "Tabular Data" export (tab, comma, semicolon or blank
    separated) with a header containing ``Frequency`` (e.g. ``Mode``,
    ``Frequency [Hz]``).

* node coordinates

  - ``NLIST`` (``LIST ALL SELECTED NODES.   DSYS= 0`` / ``NODE  X  Y  Z  THXY
    THYZ THZX``), :func:`read_nlist`;
  - a ``CDWRITE`` file (``.cdb``): ``NBLOCK`` with the Fortran format line
    (``(3i9,6e21.13e3)``, ``(3i8,6e20.13)``, ``(1i9,3e20.9e3)`` ...) parsed
    in fixed columns driven by that format line; omitted trailing fields
    are zero; the block ends at ``N,R5.3,LOC,-1`` or ``-1``.  Single
    ``N,id,x,y,z`` commands are also read.  :func:`read_cdb_nodes`;
    :func:`read_cdb_elements` reads ``EBLOCK`` (solid format) connectivity.

* mode shapes, :func:`read_prnsol`

  - ``PRNSOL,U`` (``NODE UX UY UZ USUM``), ``PRNSOL,ROT`` (``NODE ROTX ROTY
    ROTZ RSUM``) and ``PRNSOL,DOF`` blocks: ``PRINT U    NODAL SOLUTION PER
    NODE`` / ``***** POST1 NODAL DEGREE OF FREEDOM LISTING *****``;
  - the data set is taken from the ``LOAD STEP=  1  SUBSTEP=  n`` header of
    the block (also the ``SET COMMAND GOT LOAD STEP= ... SUBSTEP= ...``
    echo of ``SET,1,n``) and the frequency from ``FREQ=``;
  - column wrapping: when the requested components do not fit the page
    width MAPDL prints the node list again with the remaining columns
    (new ``NODE ...`` header line).  All column groups are merged by node
    id, so wrapped ``PRNSOL,DOF`` listings and separate U / ROT listings of
    the same set are combined;
  - numbers may run together (``-0.12345E-002-0.23456E-002``), both 2- and
    3-digit exponents and ``D`` exponents are accepted;
  - the ``MAXIMUM ABSOLUTE VALUES`` trailer ends a block.

Workbench Mechanical text export of a result object (right click ->
Export -> Text file), one file per mode, :func:`read_workbench_modal`:
``Node Number  X Location (m)  Y Location (m)  Z Location (m)
Directional Deformation (m)``.  Units in the column headers (m, mm, cm, um,
in, ft) are converted to metres.  A file may also carry several labelled
deformation columns (``Directional Deformation X (m)``, ``UX (mm)`` ...).
``Total Deformation`` is accepted with a warning: it is unsigned.

Binary results, :func:`read_rst_modal`: via ``ansys.mapdl.reader`` (only if
importable; it is never installed by FLAAP).

Normalisation and units
-----------------------
ANSYS modal shapes are mass-normalised by default (``MODOPT,...,,,ON`` /
Workbench default), so ``gen_mass = 1`` in the *solver* unit system.  Pass
``gen_mass`` (scalar or per mode) if the shapes are normalised otherwise
(e.g. ``MXPAND`` with ``NRMKEY = OFF``, max component = 1).  For MAPDL
listings ``length_scale`` / ``mass_scale`` convert solver units to SI as in
:func:`flaap.io.nastran.read_f06_modal`: coordinates and translations are
multiplied by ``length_scale`` (rotations unchanged, they are the same
physical motion), the generalised mass by ``mass_scale * length_scale**2``
(e.g. mm-tonne-s: ``length_scale=1e-3, mass_scale=1e3``).  Workbench
exports carry their length unit in the header and are converted
automatically (the mass normalisation is then that of the SI solution).
``transform`` (3x3) maps ANSYS global axes to the FLAAP aero frame (x
downstream, y right, z up) and ``origin`` is subtracted afterwards.

Synthetic sample writers
------------------------
:func:`write_mapdl_listing`, :func:`write_cdb` and :func:`write_workbench_modal`
write a :class:`ModalModel` in the above layouts.  They exist to produce
realistic test / example files (``examples/ansys/``) because ANSYS itself
cannot be run here; they are not meant to drive ANSYS.
"""
from __future__ import annotations

import datetime as _dt
import os
import re
import warnings

import numpy as np

from ..structures.modal import ModalModel

__all__ = ["read_frequencies", "read_nlist", "read_cdb_nodes", "read_cdb_elements",
           "read_prnsol", "read_mapdl_modal", "read_workbench_modal", "read_rst_modal",
           "read_ansys_modal", "write_mapdl_listing", "write_cdb", "write_workbench_modal"]

# --------------------------------------------------------------------- helpers
_FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
_FLOAT_RE = re.compile(_FLOAT)
# a number glued to the previous one: ...E-002-0.123  -> handled by finditer
_ROW_RE = re.compile(r"^\s*(\d+)\s*((?:\s*" + _FLOAT + r")+)\s*$")

_U_COMP = {"UX": 0, "UY": 1, "UZ": 2, "ROTX": 3, "ROTY": 4, "ROTZ": 5}
_KNOWN_COLS = set(_U_COMP) | {"USUM", "RSUM", "TEMP", "PRES", "VOLT", "MAG", "WARP"}

_LEN_UNITS = {"m": 1.0, "mm": 1e-3, "cm": 1e-2, "um": 1e-6, "µm": 1e-6, "μm": 1e-6,
              "in": 0.0254, "ft": 0.3048, "km": 1e3}


def _num(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def _floats(s: str) -> list[float]:
    """All numbers of a string, also when they run together."""
    return [_num(m.group(0)) for m in _FLOAT_RE.finditer(s)]


def _data_row(line: str):
    """``(id, [values])`` if the line is ``int  num num ...`` else ``None``."""
    m = _ROW_RE.match(line)
    if not m:
        return None
    return int(m.group(1)), _floats(m.group(2))


def _read_text(path: str) -> str:
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    if raw[:3] == b"\xef\xbb\xbf":
        return raw[3:].decode("utf-8", errors="replace")
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


_KEYLINE = re.compile(r"NODE|PRINT|\*\*\*|LOAD STEP|LIST|INDEX|MODE|SET|FREQ|THE FOLLOWING")


def _lines(path: str) -> list[str]:
    """Listing lines with MAPDL page headers removed."""
    raw = _read_text(path).splitlines()
    out: list[str] = []
    i, n = 0, len(raw)
    while i < n:
        ln = raw[i].rstrip("\r\n").replace("\f", "")
        if re.match(r"^1 *$|^1 +\*{5}", ln):
            # page eject: banner, product line, VERSION= line, blank, title
            j = i + 1
            while j < n and j <= i + 7 and "VERSION=" not in raw[j]:
                j += 1
            if j < n and "VERSION=" in raw[j]:
                k = j + 1
                while k < n and k <= j + 2 and not raw[k].strip():
                    k += 1
                if k < n and k <= j + 2 and not _KEYLINE.search(raw[k]) and _data_row(raw[k]) is None:
                    k += 1          # title line
                i = k
                continue
            i += 1
            continue
        if "ENGINEERING ANALYSIS SYSTEM" in ln or "VERSION=" in ln:
            i += 1
            continue
        out.append(ln)
        i += 1
    return out


def _paths(listing) -> list[str]:
    if isinstance(listing, (str, os.PathLike)):
        return [os.fspath(listing)]
    return [os.fspath(p) for p in listing]


# ----------------------------------------------------------------- frequencies
_SETLIST_ROW = re.compile(r"^\s*(\d+)\s+(" + _FLOAT + r")\s+(\d+)\s+(\d+)\s+(\d+)\s*$")
_FREQ_ROW = re.compile(r"^\s*(\d+)\s+(" + _FLOAT + r")(\s+.*)?$")


def _setlist(lines) -> list[tuple[int, float, int, int]]:
    """Rows ``(set, freq, load step, substep)`` of SET,LIST tables."""
    rows, active = [], False
    for ln in lines:
        if "INDEX OF DATA SETS" in ln:
            active = True
            continue
        if not active:
            continue
        m = _SETLIST_ROW.match(ln)
        if m:
            rows.append((int(m.group(1)), _num(m.group(2)), int(m.group(3)), int(m.group(4))))
        elif rows and ln.strip() and not re.match(r"^\s*SET\s+TIME", ln):
            active = False
    seen, out = set(), []
    for r in rows:
        if r[0] not in seen:
            seen.add(r[0])
            out.append(r)
    return sorted(out)


def _solver_freqs(lines) -> dict[int, float]:
    """``MODE  FREQUENCY ...`` tables of the modal solution output."""
    res: dict[int, float] = {}
    active, got = False, False
    for ln in lines:
        if re.search(r"\bMODE +FREQUENCY", ln):        # MAPDL prints upper case
            active, got = True, False
            continue
        if not active:
            continue
        if not ln.strip():
            continue
        m = _FREQ_ROW.match(ln)
        if m and _data_row(ln.split("!")[0]) is not None:
            res.setdefault(int(m.group(1)), _num(m.group(2)))
            got = True
        elif got or "*****" in ln:
            active = False
    return res


def _split_cells(line: str, delim: str | None) -> list[str]:
    if delim is None:
        return line.split()
    return [c.strip().strip('"') for c in line.split(delim)]


def _delimiter(header: str) -> str | None:
    for d in ("\t", ";", ","):
        if d in header:
            return d
    return None


def _table_freqs(lines) -> np.ndarray | None:
    """Workbench 'Tabular Data' style table with a ``Frequency`` column."""
    for h, hdr in enumerate(lines):
        if re.search(r"frequency", hdr, re.I) and not re.search(r"node number", hdr, re.I):
            delim = _delimiter(hdr)
            cells = _split_cells(hdr, delim)
            if delim is None:
                # join unit tokens to their name: 'Frequency [Hz]'
                merged = []
                for c in cells:
                    if merged and re.match(r"^[\[(].*[\])]$", c):
                        merged[-1] += " " + c
                    else:
                        merged.append(c)
                cells = merged
            col = next(i for i, c in enumerate(cells) if re.search(r"frequency", c, re.I))
            mcol = next((i for i, c in enumerate(cells) if re.fullmatch(r"\s*mode\s*", c, re.I)), None)
            rows = []
            for ln in lines[h + 1:]:
                if not ln.strip():
                    continue
                c = _split_cells(ln, delim)
                # tabular-data copies often have an extra leading row index
                shift = len(c) - len(cells) if len(c) > len(cells) else 0
                try:
                    f = _num(c[col + shift])
                    mo = _num(c[mcol + shift]) if mcol is not None else len(rows) + 1
                except (ValueError, IndexError):
                    if rows:
                        break
                    continue
                rows.append((mo, f))
            if rows:
                rows.sort(key=lambda r: r[0])
                return np.array([r[1] for r in rows])
    return None


def read_frequencies(path: str) -> np.ndarray:
    """Natural frequencies [Hz] ordered by data set / mode number.

    Tried in order: ``SET,LIST`` index table, ``MODE  FREQUENCY`` solver
    tables, Workbench frequency table.  Raises ``ValueError`` if none.
    """
    lines = _lines(path)
    sl = _setlist(lines)
    if sl:
        return np.array([r[1] for r in sl])
    sf = _solver_freqs(lines)
    if sf:
        return np.array([sf[k] for k in sorted(sf)])
    tf = _table_freqs(lines)
    if tf is not None:
        return tf
    raise ValueError(f"{path}: no ANSYS frequency table found")


# ----------------------------------------------------------------- nodes
def read_nlist(path: str) -> dict[int, np.ndarray]:
    """``{node: (x, y, z)}`` from ``NLIST`` output (rotation angles ignored)."""
    lines = _lines(path)
    nodes: dict[int, np.ndarray] = {}
    active = False
    for ln in lines:
        if "LIST ALL SELECTED NODES" in ln or re.match(r"^\s*NODE\s+X\s+Y\s+Z\b", ln):
            active = True
            continue
        if not active:
            continue
        if re.search(r"INDEX OF DATA SETS|NODAL SOLUTION|DEGREE OF FREEDOM LISTING|"
                     r"LIST ALL SELECTED ELEM|MODE\s+FREQUENCY|MAXIMUM ABSOLUTE", ln):
            active = False
            continue
        r = _data_row(ln)
        if r is not None and len(r[1]) >= 3:
            nodes[r[0]] = np.array(r[1][:3], float)
        elif r is not None and len(r[1]) in (1, 2):
            nodes[r[0]] = np.array(list(r[1]) + [0.0] * (3 - len(r[1])), float)
    if not nodes:
        raise ValueError(f"{path}: no NLIST node table found")
    return nodes


def _fortran_fields(fmt: str) -> list[tuple[str, int]]:
    """``(3i9,6e21.13e3)`` -> [('i',9)]*3 + [('e',21)]*6."""
    f = fmt.strip().strip("()").lower()
    out = []
    for item in f.split(","):
        m = re.match(r"^(\d*)([iefgd])(\d+)", item.strip())
        if not m:
            continue
        rep = int(m.group(1) or 1)
        out += [(m.group(2), int(m.group(3)))] * rep
    return out


def _fixed(line: str, fields) -> list[str]:
    out, pos = [], 0
    for _, w in fields:
        out.append(line[pos:pos + w])
        pos += w
    return out


def read_cdb_nodes(path: str) -> dict[int, np.ndarray]:
    """``{node: (x, y, z)}`` from a ``CDWRITE`` ``.cdb`` (NBLOCK and ``N,`` commands)."""
    lines = _read_text(path).splitlines()
    nodes: dict[int, np.ndarray] = {}
    i = 0
    while i < len(lines):
        ln = lines[i]
        up = ln.strip().upper()
        if up.startswith("NBLOCK"):
            fields = _fortran_fields(lines[i + 1])
            nint = sum(1 for k, _ in fields if k == "i")
            i += 2
            while i < len(lines):
                row = lines[i]
                u = row.strip().upper()
                if u.startswith("N,R5.3,LOC") or u.startswith("N,UNBL,LOC") or u == "-1" \
                        or not row.strip():
                    break
                cells = _fixed(row, fields)
                nid = int(cells[0])
                xyz = [0.0, 0.0, 0.0]
                for k, c in enumerate(cells[nint:nint + 3]):
                    if c.strip():
                        xyz[k] = _num(c)
                nodes[nid] = np.array(xyz)
                i += 1
        elif re.match(r"^N\s*,\s*\d+", ln.strip(), re.I):
            c = [s.strip() for s in ln.strip().split(",")]
            xyz = [_num(v) if v else 0.0 for v in (c[2:5] + ["", "", ""])[:3]]
            nodes[int(c[1])] = np.array(xyz)
        i += 1
    if not nodes:
        raise ValueError(f"{path}: no NBLOCK / N nodes found")
    return nodes


def read_cdb_elements(path: str) -> dict[int, tuple[int, np.ndarray]]:
    """``{element: (type number, node ids)}`` from the solid-format ``EBLOCK``
    of a ``.cdb`` (11 attribute fields, up to 8 nodes per line, continuation
    line for more).  Repeated trailing nodes (degenerate shapes) are kept."""
    lines = _read_text(path).splitlines()
    out: dict[int, tuple[int, np.ndarray]] = {}
    i = 0
    while i < len(lines):
        up = lines[i].strip().upper()
        if up.startswith("EBLOCK"):
            solid = "SOLID" in up
            fields = _fortran_fields(lines[i + 1])
            w = fields[0][1] if fields else 9
            i += 2
            while i < len(lines):
                row = lines[i]
                if row.strip() == "-1" or not row.strip():
                    break
                v = [int(row[k:k + w]) for k in range(0, len(row.rstrip()), w)]
                if solid:
                    nn, eid, etype = v[8], v[10], v[1]
                    conn = v[11:]
                    while len(conn) < nn and i + 1 < len(lines):
                        i += 1
                        r2 = lines[i]
                        conn += [int(r2[k:k + w]) for k in range(0, len(r2.rstrip()), w)]
                    out[eid] = (etype, np.array(conn[:nn], int))
                else:                    # non-solid (contact) format: id, sec, real, mat, ...
                    out[v[0]] = (v[1], np.array(v[5:], int))
                i += 1
        i += 1
    return out


# ----------------------------------------------------------------- PRNSOL
def read_prnsol(path: str) -> dict[int, dict]:
    """Nodal DOF listings (``PRNSOL,U / ROT / DOF``) per data set.

    Returns ``{set: {"ids": (n,), "U": (n,3), "ROT": (n,3) (if listed),
    "freq": float (if printed), "load_step": int, "substep": int}}``.  The
    set number comes from a ``SET,LIST`` table in the same file if present,
    else the substep (all load steps = 1, the modal case), else the order
    of appearance.
    """
    lines = _lines(path)
    blocks: dict[tuple[int, int], dict] = {}
    order: list[tuple[int, int]] = []
    ls, sub, freq = 1, None, None
    active, cols = False, None
    nanon = 0
    for ln in lines:
        m = re.search(r"LOAD STEP=\s*(\d+)\s+SUBSTEP=\s*(\d+)", ln)
        if m:
            ls, sub = int(m.group(1)), int(m.group(2))
            freq = None
        m = re.search(r"\bFREQ=\s*(" + _FLOAT + ")", ln)
        if m:
            freq = _num(m.group(1))
        if re.search(r"NODAL SOLUTION PER NODE|NODAL DEGREE OF FREEDOM LISTING", ln):
            if not active:
                cols = None
            active = True
            continue
        if re.search(r"MAXIMUM ABSOLUTE VALUES|LIST ALL SELECTED|INDEX OF DATA SETS", ln):
            active, cols = False, None
            continue
        if not active:
            continue
        hm = re.match(r"^\s*NODE((?:\s+[A-Z][A-Z0-9]*)+)\s*$", ln)
        if hm:
            names = hm.group(1).split()
            if all(n in _KNOWN_COLS for n in names):
                cols = names
            continue
        if cols is None:
            continue
        r = _data_row(ln)
        if r is None or len(r[1]) != len(cols):
            continue
        if sub is None:
            nanon += 1
            sub = nanon
        key = (ls, sub)
        if key not in blocks:
            blocks[key] = {"vals": {}, "has_u": False, "has_rot": False, "freq": freq}
            order.append(key)
        b = blocks[key]
        if freq is not None:
            b["freq"] = freq
        v = b["vals"].setdefault(r[0], np.zeros(6))
        for name, val in zip(cols, r[1]):
            c = _U_COMP.get(name)
            if c is not None:
                v[c] = val
                b["has_u" if c < 3 else "has_rot"] = True
    if not blocks:
        raise ValueError(f"{path}: no PRNSOL nodal DOF listing found")
    table = {(r[2], r[3]): r[0] for r in _setlist(lines)}
    if all(k in table for k in order):
        keymap = {k: table[k] for k in order}
    elif all(k[0] == 1 for k in order):
        keymap = {k: k[1] for k in order}
    else:
        keymap = {k: i + 1 for i, k in enumerate(order)}
    out = {}
    for k in order:
        b = blocks[k]
        ids = np.array(sorted(b["vals"]))
        a = np.array([b["vals"][i] for i in ids])
        d = {"ids": ids, "load_step": k[0], "substep": k[1]}
        if b["has_u"]:
            d["U"] = a[:, :3]
        if b["has_rot"]:
            d["ROT"] = a[:, 3:]
        if b["freq"] is not None:
            d["freq"] = b["freq"]
        out[keymap[k]] = d
    return out


# ----------------------------------------------------------------- assemble
def _finish(ids, coords, sh, freq_hz, gen_mass, damping, length_scale, mass_scale,
            transform, origin, labels, name) -> ModalModel:
    nm = sh.shape[2]
    T = np.eye(3) if transform is None else np.asarray(transform, float)
    nodes = (coords * length_scale) @ T.T - np.asarray(origin, float)
    out = np.empty_like(sh)
    out[:, :3] = np.einsum("ij,njm->nim", T, sh[:, :3]) * length_scale
    out[:, 3:] = np.einsum("ij,njm->nim", T, sh[:, 3:])
    gm = np.broadcast_to(np.asarray(gen_mass, float), (nm,)) * mass_scale * length_scale ** 2
    zeta = np.broadcast_to(np.asarray(damping, float), (nm,)).copy()
    return ModalModel(nodes, out, 2 * np.pi * np.asarray(freq_hz, float), gm.copy(),
                      damping=zeta, name=name, labels=labels)


def _node_source(nodes, paths) -> dict[int, np.ndarray]:
    if isinstance(nodes, dict):
        return {int(k): np.asarray(v, float) for k, v in nodes.items()}
    if nodes is not None:
        p = os.fspath(nodes)
        txt = _read_text(p)
        if p.lower().endswith(".cdb") or re.search(r"^\s*NBLOCK", txt, re.M | re.I):
            return read_cdb_nodes(p)
        return read_nlist(p)
    for p in paths:
        try:
            return read_nlist(p)
        except ValueError:
            pass
    raise ValueError("no node coordinates: pass nodes= (CDB, NLIST file or dict) "
                     "or include an NLIST listing")


def read_mapdl_modal(listing, nodes=None, nmodes: int | None = None, gen_mass=1.0,
                     skip_below_hz: float = 1e-3, length_scale: float = 1.0,
                     mass_scale: float = 1.0, transform: np.ndarray | None = None,
                     origin=(0.0, 0.0, 0.0), damping=0.0, node_ids=None,
                     frequencies=None) -> ModalModel:
    """Modal model from MAPDL text output.

    ``listing``: one path or a list of paths holding the frequency table
    (``SET,LIST`` or solver summary), the ``PRNSOL`` blocks and optionally
    the ``NLIST``.  ``nodes``: ``.cdb`` / NLIST path or ``{id: xyz}`` (default:
    NLIST found in ``listing``).  ``node_ids`` restricts / orders the nodes
    (default: nodes listed in the first kept set that have coordinates).
    ``frequencies`` [Hz] overrides the tables.  Sets below ``skip_below_hz``
    (rigid-body modes) are skipped.  ``gen_mass``: generalised mass in
    solver units (1 for mass-normalised shapes), scalar or per kept mode.
    Rotations are zero if no ``PRNSOL,ROT``/``DOF`` was printed.
    """
    paths = _paths(listing)
    sets: dict[int, dict] = {}
    freqs = None if frequencies is None else np.asarray(frequencies, float)
    for p in paths:
        try:
            for k, d in read_prnsol(p).items():
                if k in sets:
                    old = sets[k]
                    for f in ("U", "ROT"):
                        if f in d and f not in old:
                            # align by ids
                            pos = {g: r for r, g in enumerate(d["ids"])}
                            old[f] = np.array([d[f][pos[g]] if g in pos else np.zeros(3)
                                               for g in old["ids"]])
                else:
                    sets[k] = d
        except ValueError:
            pass
        if freqs is None:
            try:
                freqs = read_frequencies(p)
            except ValueError:
                pass
    if not sets:
        raise ValueError("no PRNSOL mode shapes found in " + ", ".join(paths))
    coords = _node_source(nodes, paths)

    def fr(k):
        if freqs is not None and k - 1 < freqs.size:
            return float(freqs[k - 1])
        if "freq" in sets[k]:
            return float(sets[k]["freq"])
        raise ValueError(f"no frequency for set {k}")

    keep = [k for k in sorted(sets) if fr(k) >= skip_below_hz]
    if nmodes is not None:
        keep = keep[:nmodes]
    if not keep:
        raise ValueError("no elastic modes found")
    ids = np.asarray(node_ids if node_ids is not None
                     else [g for g in sets[keep[0]]["ids"] if int(g) in coords], int)
    xyz = np.array([coords[int(g)] for g in ids])
    sh = np.zeros((ids.size, 6, len(keep)))
    for c, k in enumerate(keep):
        d = sets[k]
        pos = {int(g): r for r, g in enumerate(d["ids"])}
        rows = np.array([pos.get(int(g), -1) for g in ids])
        ok = rows >= 0
        if "U" in d:
            sh[ok, :3, c] = d["U"][rows[ok]]
        if "ROT" in d:
            sh[ok, 3:, c] = d["ROT"][rows[ok]]
    m = _finish(ids, xyz, sh, [fr(k) for k in keep], gen_mass, damping, length_scale,
                mass_scale, transform, origin, [f"mode{k}" for k in keep], "ansys")
    m.node_ids = ids
    return m


# ----------------------------------------------------------------- Workbench
def _unit_of(header: str) -> float:
    m = re.search(r"[\[(]\s*([^\])]+?)\s*[\])]\s*$", header.strip())
    if not m:
        return 1.0
    u = m.group(1).strip()
    if u not in _LEN_UNITS:
        raise ValueError(f"unknown length unit '{u}' in column '{header}'")
    return _LEN_UNITS[u]


def _axis_of(header: str) -> int | None:
    h = re.sub(r"[\[(][^\])]*[\])]", "", header)       # drop the unit
    m = re.search(r"\bU([XYZ])\b", h, re.I) or re.search(r"\b([XYZ])\b(?:\s*Axis)?", h)
    return "XYZ".index(m.group(1).upper()) if m else None


def _read_wb_table(path: str):
    """Parse a Workbench result export -> (ids, xyz [m], {col header: values [m]})."""
    lines = [ln for ln in _read_text(path).splitlines() if ln.strip()]
    h = next((i for i, ln in enumerate(lines) if re.search(r"node\s*number", ln, re.I)), None)
    if h is None:
        raise ValueError(f"{path}: no 'Node Number' header")
    delim = _delimiter(lines[h])
    if delim is None:
        raise ValueError(f"{path}: cannot determine the column separator")
    hdr = _split_cells(lines[h], delim)
    data = []
    for ln in lines[h + 1:]:
        c = _split_cells(ln, delim)
        try:
            data.append([_num(v) for v in c[:len(hdr)]])
        except ValueError:
            break
    a = np.array(data, float)
    ids = a[:, 0].astype(int)
    loc = {}
    defs = {}
    for j, name in enumerate(hdr[1:], 1):
        if re.search(r"location", name, re.I):
            ax = re.search(r"\b([XYZ])\s*Location", name, re.I)
            loc["XYZ".index(ax.group(1).upper())] = a[:, j] * _unit_of(name)
        elif re.search(r"deformation|displacement|\bU[XYZ]\b", name, re.I):
            defs[name] = a[:, j] * _unit_of(name)
    xyz = np.stack([loc.get(k, np.zeros(ids.size)) for k in range(3)], 1) if loc else None
    return ids, xyz, defs


def read_workbench_modal(mode_files, freq_file: str | None = None, frequencies=None,
                         direction="Z", gen_mass=1.0, nmodes: int | None = None,
                         length_scale: float = 1.0, mass_scale: float = 1.0,
                         transform: np.ndarray | None = None, origin=(0.0, 0.0, 0.0),
                         damping=0.0) -> ModalModel:
    """Modal model from Workbench Mechanical text exports, one entry per mode.

    ``mode_files[i]`` is a path (columns ``Node Number``, ``X/Y/Z Location``
    and one or more deformation columns) or a dict ``{"X": path, "Y": path,
    "Z": path}`` of directional-deformation exports.  A single unlabelled
    ``Directional Deformation`` column is assigned to ``direction`` (a
    letter, or one letter per mode).  ``Total Deformation`` is unsigned:
    accepted with a warning and assigned to ``direction``.  Frequencies [Hz]
    from ``frequencies`` or ``freq_file`` (Workbench frequency table).
    Lengths are converted to metres from the header units; ``length_scale``
    / ``mass_scale`` are applied on top as in :func:`read_mapdl_modal`.
    Rotations are not exported by Workbench and are set to zero.
    """
    if frequencies is not None:
        freqs = np.asarray(frequencies, float)
    elif freq_file is not None:
        freqs = read_frequencies(freq_file)
    else:
        raise ValueError("give frequencies= or freq_file=")
    files = list(mode_files)
    if nmodes is not None:
        files = files[:nmodes]
    dirs = [direction] * len(files) if isinstance(direction, str) else list(direction)
    ids0, xyz0 = None, None
    cols = []
    for m, entry in enumerate(files):
        u = {}
        items = entry.items() if isinstance(entry, dict) else [(dirs[m], entry)]
        for ax_given, p in items:
            ids, xyz, defs = _read_wb_table(os.fspath(p))
            if ids0 is None:
                ids0, xyz0 = ids, xyz
            pos = {g: r for r, g in enumerate(ids)}
            rows = np.array([pos.get(g, -1) for g in ids0])
            for name, v in defs.items():
                if re.search(r"total", name, re.I):
                    warnings.warn(f"{p}: 'Total Deformation' is unsigned - the mode shape "
                                  "sign (and nodal lines) are lost; export Directional "
                                  "Deformation instead", stacklevel=2)
                ax = _axis_of(name)
                if ax is None or re.search(r"total", name, re.I):
                    ax = "XYZ".index(str(ax_given).upper())
                vv = np.where(rows >= 0, v[np.maximum(rows, 0)], 0.0)
                u[ax] = vv
        cols.append(u)
    if xyz0 is None:
        raise ValueError("no X/Y/Z Location columns in the Workbench export")
    nm = len(files)
    sh = np.zeros((ids0.size, 6, nm))
    for m, u in enumerate(cols):
        for ax, v in u.items():
            sh[:, ax, m] = v
    if freqs.size < nm:
        raise ValueError(f"{freqs.size} frequencies for {nm} mode files")
    mm = _finish(ids0, xyz0, sh, freqs[:nm], gen_mass, damping, length_scale, mass_scale,
                 transform, origin, [f"mode{i + 1}" for i in range(nm)], "ansys-workbench")
    mm.node_ids = ids0
    return mm


def _workbench_dir(directory: str):
    files = sorted(os.listdir(directory))
    modes: dict[int, dict | str] = {}
    freq_file = None
    for f in files:
        p = os.path.join(directory, f)
        if not os.path.isfile(p) or not f.lower().endswith((".txt", ".csv", ".tsv", ".xls")):
            continue
        head = _read_text(p)[:2000]
        if re.search(r"node\s*number", head, re.I):
            m = re.search(r"(\d+)", f)
            if not m:
                continue
            k = int(m.group(1))
            a = re.search(r"[_\-\s]([XYZ])\.\w+$", f, re.I)
            if a:
                d = modes.setdefault(k, {})
                if isinstance(d, dict):
                    d[a.group(1).upper()] = p
            else:
                modes[k] = p
        elif re.search(r"frequency", head, re.I):
            freq_file = p
    return [modes[k] for k in sorted(modes)], freq_file


# ----------------------------------------------------------------- binary rst
def read_rst_modal(path: str, nmodes: int | None = None, gen_mass=1.0,
                   skip_below_hz: float = 1e-3, length_scale: float = 1.0,
                   mass_scale: float = 1.0, transform: np.ndarray | None = None,
                   origin=(0.0, 0.0, 0.0), damping=0.0) -> ModalModel:
    """Modal model from a binary ``.rst`` via ``ansys.mapdl.reader``.

    Raises ``ImportError`` if ``ansys-mapdl-reader`` is not installed (it is
    an optional dependency and never installed by FLAAP).
    """
    try:
        from ansys.mapdl import reader as pyreader  # type: ignore
    except ImportError as exc:   # pragma: no cover - depends on the environment
        raise ImportError("reading .rst files needs the optional package "
                          "'ansys-mapdl-reader' (pip install ansys-mapdl-reader); "
                          "alternatively export SET,LIST / NLIST / PRNSOL text "
                          "listings and use read_mapdl_modal") from exc
    rst = pyreader.read_binary(os.fspath(path))   # pragma: no cover
    freqs = np.asarray(rst.time_values, float)    # pragma: no cover
    keep = [i for i, f in enumerate(freqs) if f >= skip_below_hz]   # pragma: no cover
    if nmodes is not None:                         # pragma: no cover
        keep = keep[:nmodes]
    nnum = np.asarray(rst.mesh.nnum, int)          # pragma: no cover
    xyz = np.asarray(rst.mesh.nodes, float)[:, :3]  # pragma: no cover
    sh = np.zeros((nnum.size, 6, len(keep)))       # pragma: no cover
    for c, i in enumerate(keep):                   # pragma: no cover
        ids, disp = rst.nodal_displacement(i)
        names = list(rst.result_dof(i))
        pos = {g: r for r, g in enumerate(np.asarray(ids, int))}
        rows = np.array([pos.get(g, -1) for g in nnum])
        ok = rows >= 0
        for j, nmn in enumerate(names):
            comp = _U_COMP.get(str(nmn).upper())
            if comp is not None:
                sh[ok, comp, c] = np.asarray(disp)[rows[ok], j]
    mm = _finish(nnum, xyz, sh, freqs[keep], gen_mass, damping, length_scale,  # pragma: no cover
                 mass_scale, transform, origin, [f"mode{i + 1}" for i in keep], "ansys-rst")
    mm.node_ids = nnum                             # pragma: no cover
    return mm                                      # pragma: no cover


# ----------------------------------------------------------------- dispatcher
def read_ansys_modal(path, **kw) -> ModalModel:
    """Dispatch on the input: ``.rst``/``.rstp`` -> :func:`read_rst_modal`;
    a directory -> Workbench exports (``Node Number`` files + frequency
    table) or the MAPDL listings it contains; a list of paths or a text
    listing -> :func:`read_mapdl_modal` (pass ``nodes=`` for a ``.cdb``).
    """
    if not isinstance(path, (str, os.PathLike)):
        return read_mapdl_modal(path, **kw)
    p = os.fspath(path)
    ext = os.path.splitext(p)[1].lower()
    if os.path.isdir(p):
        modes, ff = _workbench_dir(p)
        if modes:
            if ff is not None and "frequencies" not in kw:
                kw.setdefault("freq_file", ff)
            return read_workbench_modal(modes, **kw)
        outs = [os.path.join(p, f) for f in sorted(os.listdir(p))
                if f.lower().endswith((".out", ".lis", ".txt", ".log"))]
        if not outs:
            raise ValueError(f"{p}: no ANSYS exports found")
        if "nodes" not in kw:
            cdb = [os.path.join(p, f) for f in sorted(os.listdir(p)) if f.lower().endswith(".cdb")]
            if cdb:
                kw["nodes"] = cdb[0]
        return read_mapdl_modal(outs, **kw)
    if ext in (".rst", ".rstp", ".rth"):
        return read_rst_modal(p, **kw)
    if ext == ".cdb":
        raise ValueError("a .cdb holds no results: pass the listing and nodes='model.cdb'")
    head = _read_text(p)[:4000]
    if re.search(r"node\s*number", head, re.I):
        return read_workbench_modal([p], **kw)
    return read_mapdl_modal(p, **kw)


# ======================================================================
# synthetic sample writers (for tests / examples)
# ======================================================================
_STAMP = _dt.datetime(2026, 11, 26, 10, 15, 32)


def _efmt(v: float, digits: int = 5, expd: int = 3) -> str:
    """MAPDL-style ``0.12345E-002`` (zero printed as ``0.0000``)."""
    if v == 0.0 or not np.isfinite(v):
        return "0.0000"
    a = abs(v)
    e = int(np.floor(np.log10(a))) + 1
    m = round(a / 10.0 ** e, digits)
    if m >= 1.0:
        m /= 10.0
        e += 1
    s = f"{m:.{digits}f}E{e:+0{expd + 1}d}"
    return ("-" if v < 0 else "") + s


class _Pager:
    """Accumulates MAPDL output with page headers every ``lpp`` lines."""

    def __init__(self, title: str, lpp: int = 56, release: str = "2023 R1", version: str = "23.1"):
        self.title, self.lpp = title, lpp
        self.rel, self.ver = release, version
        self.buf: list[str] = []
        self.on_page = 0
        self.cp = 1.0

    def page(self):
        self.cp += 0.047
        stamp = _STAMP.strftime("%H:%M:%S  %b %d, %Y").upper()
        self.buf += ["1",
                     f" ***** Ansys - ENGINEERING ANALYSIS SYSTEM  RELEASE {self.rel:<12s}"
                     f"{self.ver:>6s}     *****",
                     " Ansys Mechanical Enterprise                       ",
                     f" 00000000  VERSION=WINDOWS x64   {stamp} CP={self.cp:12.3f}",
                     "",
                     f" {self.title:<78s}",
                     ""]
        self.on_page = 7

    def add(self, *lines: str, header: list[str] | None = None):
        for ln in lines:
            if self.on_page >= self.lpp:
                self.page()
                if header:
                    self.buf += header
                    self.on_page += len(header)
            self.buf.append(ln)
            self.on_page += 1

    def need(self, n: int):
        """Start a new page unless ``n`` more lines fit (keeps headers with data)."""
        if self.on_page + n > self.lpp:
            self.page()

    def text(self) -> str:
        return "\n".join(self.buf) + "\n"


def _native(modal: ModalModel, unit_length: float, unit_mass: float):
    """Shapes / coordinates in the solver unit system, mass-normalised there."""
    gm_native = modal.gen_mass / (unit_mass * unit_length ** 2)
    a = 1.0 / np.sqrt(gm_native)
    sh = modal.shapes.copy() * a[None, None, :]
    sh[:, :3] /= unit_length
    return modal.nodes / unit_length, sh


def write_mapdl_listing(modal: ModalModel, path: str, node_ids=None,
                        title: str = "FLAAP synthetic modal model", rotations: bool = True,
                        setlist: bool = True, nlist: bool = True, solver_summary: bool = True,
                        dof: bool = False, wrap: int = 4, digits: int = 5,
                        unit_length: float = 1.0, unit_mass: float = 1.0,
                        lines_per_page: int = 56) -> None:
    """Write ``modal`` as an MAPDL batch-output listing (synthetic sample).

    Contents: the Block Lanczos frequency summary (``solver_summary``), the
    ``SET,LIST`` index (``setlist``), ``NLIST`` (``nlist``) and per mode a
    ``SET,1,n`` echo followed by ``PRNSOL,U`` and ``PRNSOL,ROT``
    (``rotations``), or by a ``PRNSOL,DOF`` listing wrapped after ``wrap``
    columns (``dof=True``).  Values use the MAPDL ``0.12345E-002`` format
    with ``digits`` significant digits.  Output is in the unit system
    ``unit_length`` [m] / ``unit_mass`` [kg] with shapes mass-normalised
    there (``modal`` is assumed SI).
    """
    nn = modal.nodes.shape[0]
    ids = np.arange(1, nn + 1) if node_ids is None else np.asarray(node_ids, int)
    xyz, sh = _native(modal, unit_length, unit_mass)
    f = modal.freq_hz
    P = _Pager(title, lines_per_page)
    P.page()
    if solver_summary:
        P.add(" *****  MAPDL SOLVE    COMMAND  *****", "",
              " *** NOTE ***                            CP =       0.844   TIME= 10:15:31",
              " The modal analysis is using the Block Lanczos mode extraction method.", "",
              " *** FREQUENCIES FROM BLOCK LANCZOS ITERATION ***", "",
              "  MODE    FREQUENCY (HERTZ)", "")
        for i, fi in enumerate(f):
            P.add(f"  {i + 1:4d}     {fi:<16.11f}")
        P.add("", " *** MAPDL BINARY FILE STATISTICS", "  BUFFER SIZE USED= 16384", "")
    P.add(" FINISH SOLUTION PROCESSING", "", "", " *****  ROUTINE COMPLETED  *****  CP =         1.203",
          "", " /POST1", "", " *****  MAPDL RESULTS INTERPRETATION (POST1)  *****", "")
    if setlist:
        P.need(8)
        P.add(" SET,LIST", "", " *****  INDEX OF DATA SETS ON RESULTS FILE  *****", "",
              "   SET   TIME/FREQ    LOAD STEP   SUBSTEP  CUMULATIVE")
        for i, fi in enumerate(f):
            P.add(f"  {i + 1:4d}{fi:12.4f}  {1:10d}{i + 1:10d}{i + 1:10d}")
        P.add("")
    if nlist:
        hdr = ["", " LIST ALL SELECTED NODES.   DSYS=      0", "",
               "    NODE        X                   Y                   Z             "
               "THXY     THYZ     THZX"]
        P.need(len(hdr) + 4)
        P.add(" NLIST", *hdr)
        for g, x in zip(ids, xyz):
            P.add(f"{g:9d}{x[0]:20.13f}{x[1]:20.13f}{x[2]:20.13f}"
                  f"{0.0:9.2f}{0.0:9.2f}{0.0:9.2f}", header=hdr)
        P.add("")

    def block(item: str, names: list[str], comps: list[np.ndarray], m: int):
        head = [f" PRINT {item:<4s} NODAL SOLUTION PER NODE", "",
                "  ***** POST1 NODAL DEGREE OF FREEDOM LISTING *****", "",
                f"  LOAD STEP={1:6d}  SUBSTEP={m + 1:6d}",
                f"   FREQ={f[m]:11.3f}      LOAD CASE=   0", "",
                "  THE FOLLOWING DEGREE OF FREEDOM RESULTS ARE IN THE GLOBAL COORDINATE SYSTEM",
                ""]
        colh = "    NODE" + "".join(f"{n:>9s}    " for n in names)
        P.need(len(head) + 4)
        P.add(*head, colh, header=head + [colh])
        A = np.stack(comps, 1)
        for g, row in zip(ids, A):
            P.add(f"{g:9d}" + "".join(f"{_efmt(v, digits):>13s}" for v in row),
                  header=head + [colh])
        return A

    def maxabs(A):
        imax = np.abs(A).argmax(0)
        P.add("", " MAXIMUM ABSOLUTE VALUES",
              " NODE " + "".join(f"{ids[i]:13d}" for i in imax)[1:],
              " VALUE" + "".join(f"{_efmt(A[i, j], digits):>13s}" for j, i in enumerate(imax)), "")

    for m in range(modal.nmodes):
        P.add("", f" SET,1,{m + 1}", "",
              f" USE LOAD STEP     1  SUBSTEP     {m + 1}  FOR LOAD CASE 0 ", "",
              f" SET COMMAND GOT LOAD STEP=     1  SUBSTEP=     {m + 1}  CUMULATIVE ITERATION="
              f"{m + 1:6d}",
              f"   FREQ={f[m]:11.3f}      LOAD CASE=   0", f"   TITLE= {title}", "")
        U = sh[:, :3, m]
        R = sh[:, 3:, m]
        groups = []
        if dof:
            names = ["UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"]
            vals = [U[:, 0], U[:, 1], U[:, 2], R[:, 0], R[:, 1], R[:, 2]]
            P.add(" PRNSOL,DOF", "")
            for s in range(0, 6, wrap):
                block("DOF", names[s:s + wrap], vals[s:s + wrap], m)
                P.add("")
            groups.append(np.stack(vals, 1))
        else:
            P.add(" PRNSOL,U", "")
            groups.append(block("U", ["UX", "UY", "UZ", "USUM"],
                                [U[:, 0], U[:, 1], U[:, 2], np.linalg.norm(U, axis=1)], m))
            if rotations:
                maxabs(groups.pop())
                P.add(" PRNSOL,ROT", "")
                groups.append(block("ROT", ["ROTX", "ROTY", "ROTZ", "RSUM"],
                                    [R[:, 0], R[:, 1], R[:, 2], np.linalg.norm(R, axis=1)], m))
        for A in groups:
            maxabs(A)
    P.add(" FINISH", "", " EXIT THE MAPDL POST1 DATABASE PROCESSOR", "",
          " ***** ROUTINE COMPLETED *****  CP =         1.531", "")
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(P.text())


def _cdb_e(v: float, width: int, prec: int, expd: int) -> str:
    s = np.format_float_scientific(float(v), precision=prec, exp_digits=expd, unique=False)
    return s.upper().rjust(width)


def write_cdb(model, path: str, fmt: str = "(3i9,6e21.13e3)", node_ids=None,
              elements=None, title: str = "FLAAP synthetic model",
              unit_length: float = 1.0) -> None:
    """Write nodes (and optional 4-node shell elements) as a ``CDWRITE``-style
    ``.cdb`` (synthetic sample).  ``model``: :class:`ModalModel`, ``FEModel``
    or a (n,3) array.  ``fmt``: NBLOCK format, e.g. ``(3i9,6e21.13e3)``
    (release >= 2020) or ``(3i8,6e20.13)`` (older).  ``elements`` (ne,4)
    zero-based connectivity written as SHELL181 in an EBLOCK (default:
    ``model.fe.elements['quad4']`` if available)."""
    nodes = np.asarray(getattr(model, "nodes", model), float) / unit_length
    nn = nodes.shape[0]
    ids = np.arange(1, nn + 1) if node_ids is None else np.asarray(node_ids, int)
    if elements is None:
        fe = getattr(model, "fe", None) or (model if hasattr(model, "elements") else None)
        if fe is not None and getattr(fe, "elements", None):
            elements = fe.elements.get("quad4")
    fields = _fortran_fields(fmt)
    iw = fields[0][1]
    ef = re.search(r"e(\d+)\.(\d+)(?:e(\d+))?", fmt.lower())
    ew, ep, ee = int(ef.group(1)), int(ef.group(2)), int(ef.group(3) or 2)
    L = [f"/COM,ANSYS RELEASE 2023 R1  BUILD 23.1      UP20221122       "
         f"{_STAMP.strftime('%H:%M:%S')}",
         "/PREP7", "/NOPR", f"/TITLE,{title}",
         "*IF,_CDRDOFF,EQ,1,THEN     !if solid model was read in",
         "_CDRDOFF=             !reset flag, numoffs already performed",
         "*ELSE              !offset database for the following FE model",
         f"NUMOFF,NODE,{int(ids.max()):10d}"]
    ne = 0 if elements is None else len(elements)
    if ne:
        L += [f"NUMOFF,ELEM,{ne:10d}", f"NUMOFF,MAT ,{1:10d}", f"NUMOFF,TYPE,{1:10d}"]
    L += ["*ENDIF", "*SET,_BUTTON ,  0.000000000000", "DOF,DELETE"]
    if ne:
        L += [f"ET,{1:8d},181", "KEYOPT,       1, 1,       0", "KEYOPT,       1, 3,       2"]
    L += [f"NBLOCK,6,SOLID,{int(ids.max()):10d},{nn:10d}", fmt]
    for g, x in zip(ids, nodes):
        nz = np.nonzero(x)[0]
        last = nz.max() + 1 if nz.size else 0         # CDWRITE omits trailing zeros
        L.append(f"{g:{iw}d}{0:{iw}d}{0:{iw}d}" +
                 "".join(_cdb_e(v, ew, ep, ee) for v in x[:last]))
    L.append(f"N,R5.3,LOC,{-1:{iw}d},")
    if ne:
        L += [f"EBLOCK,19,SOLID,{ne:10d},{ne:10d}", "(19i9)"]
        for e, q in enumerate(np.asarray(elements, int)):
            att = [1, 1, 1, 1, 0, 0, 0, 0, 4, 0, e + 1]
            L.append("".join(f"{v:9d}" for v in att + list(ids[q])))
        L.append(f"{-1:9d}")
    L += ["/GO", "FINISH"]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


def write_workbench_modal(modal: ModalModel, directory: str, directions=("Z",),
                          node_ids=None, unit: str = "m", total: bool = False,
                          sig: int = 6) -> list[str]:
    """Write Workbench Mechanical style text exports (synthetic sample):
    ``mode_<i>_<D>.txt`` per mode and direction (or ``mode_<i>_total.txt``
    for ``total=True``) and ``frequencies.txt`` (Tabular Data of the modal
    analysis).  Returns the list of mode files."""
    os.makedirs(directory, exist_ok=True)
    s = _LEN_UNITS[unit]
    nn = modal.nodes.shape[0]
    ids = np.arange(1, nn + 1) if node_ids is None else np.asarray(node_ids, int)

    def e(v):
        return np.format_float_scientific(float(v), precision=sig - 1, exp_digits=3,
                                          unique=False).lower()

    def num(v):
        return "0." if v == 0 else e(v)

    files = []
    xyz = modal.nodes / s
    for m in range(modal.nmodes):
        items = [("total", np.linalg.norm(modal.shapes[:, :3, m], axis=1))] if total else \
            [(d, modal.shapes[:, "XYZ".index(d), m]) for d in directions]
        for d, v in items:
            name = "Total Deformation" if total else "Directional Deformation"
            fn = os.path.join(directory, f"mode_{m + 1}_{d}.txt")
            with open(fn, "w", encoding="utf-8", newline="\r\n") as fh:
                fh.write(f"Node Number\tX Location ({unit})\tY Location ({unit})\t"
                         f"Z Location ({unit})\t{name} ({unit})\n")
                for g, x, val in zip(ids, xyz, v / s):
                    fh.write(f"{g}\t{num(x[0])}\t{num(x[1])}\t{num(x[2])}\t{num(val)}\n")
            files.append(fn)
    with open(os.path.join(directory, "frequencies.txt"), "w", encoding="utf-8",
              newline="\r\n") as fh:
        fh.write("\tMode\tFrequency [Hz]\n")
        for m, fz in enumerate(modal.freq_hz):
            fh.write(f"{m + 1}\t{m + 1}.\t{fz:.4f}\n")
    return files
