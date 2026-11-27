"""CalculiX ``*FREQUENCY`` results (``.frd`` + optional ``.dat``) -> :class:`ModalModel`.

CalculiX (``ccx``) writes nodal results of a ``*FREQUENCY`` step to the
``.frd`` file when the step holds ``*NODE FILE`` with ``U``.  This module
reads that file (ASCII; short or long node format), keeps the ``DISP`` blocks
of the eigenmodes and builds a modal model like the other importers.

What the ``.frd`` file holds (cgx manual "Result Format" [1], ccx ``frd.c``
and ``arpack.c`` [2])
-----------------------------------------------------------------------------
* ``1C`` model header, ``1U`` user lines, ``1P`` parameter lines.
* ``2C`` node block ``(1X,'   2','C',18X,I12,37X,I1)``: the last field is
  the format (0 short, 1 long, 2/3 binary; ``ccx`` writes 1 for ASCII).  Node
  records ``(1X,'-1',I5|I10,3E12.5)``, ``-3`` ends the block.
* ``3C`` element block, same header.  ``-1`` record: element, type, group,
  material; ``-2`` records: connectivity (15I5 short, 10I10 long).  Types:
  1 he8, 2 pe6, 3 te4, 4 he20, 5 pe15, 6 te10, 7 tr3, 8 tr6, 9 qu4, 10 qu8,
  11 be2, 12 be3.  Shell and beam elements appear *expanded* to solids.
* ``100C`` result block ``(1X,' 100','C',6A1,E12.5,I12,20A1,I2,I5,10A1,I2)``:
  set name, VALUE, number of nodes, text, ICTYPE (2 = frequency), step
  number, analysis name, format.  For a ``*FREQUENCY`` step VALUE is the
  eigenfrequency in cycles per time unit (``sqrt(lambda)/2 pi``).  Then
  ``-4  DISP  4  1`` and four ``-5`` component lines (D1 D2 D3 ALL; ALL is
  computed) and the data records ``(1X,I2,I5|I10,6E12.5)``; ``-3`` ends.
* ``9999`` ends the file.
* Eigenvectors are normalised to *unit generalised mass* (``arpack.c``
  rescales them so that phi^T M phi = 1).  Cyclic-symmetry runs also write
  ``DISPI`` (imaginary part); it is ignored.
* The ``.dat`` file holds the table ``E I G E N V A L U E   O U T P U T``
  (mode, eigenvalue, rad/time, cycles/time, imaginary); it gives the
  frequencies with more digits than the ``.frd`` header (E12.5).

Units are the solver's consistent units: ``length_scale`` converts lengths
to m, ``mass_scale`` masses to kg (unit generalised mass in solver units
becomes ``mass_scale * length_scale**2`` in SI).

Sources: [1] G. Dhondt, K. Wittig, CalculiX GraphiX manual, "Result Format"
(http://www.dhondt.de/cgx_2.21.pdf, section 11).  [2] CalculiX CrunchiX
source ``frd.c``, ``arpack.c`` (http://www.dhondt.de, GPL).  The sample
writers below produce files in that layout; they are not ccx output.
"""
from __future__ import annotations

import math
import os
import warnings

import numpy as np

from ..structures.modal import ModalModel
from .ansys import _finish

FRD_TYPES = {1: "he8", 2: "pe6", 3: "te4", 4: "he20", 5: "pe15", 6: "te10", 7: "tr3", 8: "tr6",
             9: "qu4", 10: "qu8", 11: "be2", 12: "be3"}
_NNODE = {1: 8, 2: 6, 3: 4, 4: 20, 5: 15, 6: 10, 7: 3, 8: 6, 9: 4, 10: 8, 11: 2, 12: 3}


def _f(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    try:
        return float(s)
    except ValueError:                      # Fortran "1.0-100" without E
        for i in range(len(s) - 1, 0, -1):
            if s[i] in "+-" and s[i - 1].isdigit():
                return float(s[:i] + "E" + s[i:])
        raise


def _floats12(s: str, n: int | None = None) -> list[float]:
    out = []
    for i in range(0, len(s.rstrip()), 12):
        out.append(_f(s[i:i + 12]))
    return out if n is None else out[:n]


def _id_width(line: str, nval: int | None, fmt: int | None) -> int:
    if fmt == 0:
        return 5
    if fmt == 1:
        return 10
    body = len(line.rstrip()) - 3
    for w in (10, 5):
        if (body - w) > 0 and (body - w) % 12 == 0 and (nval is None or (body - w) // 12 == nval):
            return w
    return 10


def _open_text(path: str) -> list[str]:
    with open(path, "rb") as fh:
        raw = fh.read()
    try:
        txt = raw.decode("ascii")
    except UnicodeDecodeError:
        txt = raw.decode("latin-1")
    return txt.splitlines()


def read_frd(path: str) -> dict:
    """Parse an ASCII ``.frd`` file.

    Returns ``{"nodes": {id: xyz}, "elements": {id: (type, conn)},
    "blocks": [{"name", "value", "ictype", "numstp", "analysis", "text",
    "components", "data": {id: array}}], "header": [lines]}``.
    """
    lines = _open_text(path)
    nodes: dict[int, np.ndarray] = {}
    elements: dict[int, tuple[int, np.ndarray]] = {}
    blocks: list[dict] = []
    header: list[str] = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        key = line[:6].strip()
        if key.startswith("9999"):
            break
        if key in ("1C", "1U") or key.startswith("1P") or key.startswith("1U"):
            header.append(line[5:].rstrip() if not key.startswith("1C") else line.rstrip())
            i += 1
            continue
        if key.startswith("2C"):
            fmt = int(line[73:75] or 1) if len(line) > 73 and line[73:75].strip() else None
            if fmt in (2, 3):
                raise ValueError(f"{path}: binary .frd (format {fmt}) is not supported; "
                                 "run ccx without *NODE FILE binary output (OUTPUT=3D ASCII)")
            i += 1
            while i < n and not lines[i].startswith(" -3"):
                rec = lines[i]
                if rec.startswith(" -1"):
                    w = _id_width(rec, 3, fmt)
                    nid = int(rec[3:3 + w])
                    nodes[nid] = np.array(_floats12(rec[3 + w:], 3), float)
                i += 1
            i += 1
            continue
        if key.startswith("3C"):
            fmt = int(line[73:75]) if len(line) > 74 and line[73:75].strip() else None
            if fmt in (2, 3):
                raise ValueError(f"{path}: binary .frd element block is not supported")
            i += 1
            cur, typ, conn = None, 0, []
            while i < n and not lines[i].startswith(" -3"):
                rec = lines[i]
                if rec.startswith(" -1"):
                    if cur is not None:
                        elements[cur] = (typ, np.array(conn, int))
                    parts = rec[3:].split()
                    cur, typ, conn = int(parts[0]), int(parts[1]), []
                elif rec.startswith(" -2"):
                    body = rec[3:]
                    w = 5 if fmt == 0 else 10
                    if fmt is None:
                        conn.extend(int(x) for x in body.split())
                    else:
                        conn.extend(int(body[j:j + w]) for j in range(0, len(body.rstrip()), w))
                i += 1
            if cur is not None:
                elements[cur] = (typ, np.array(conn, int))
            i += 1
            continue
        if key.startswith("100C"):
            s = line.ljust(80)
            blk = {"setname": s[6:12].strip(), "value": _f(s[12:24]), "numnod": int(s[24:36] or 0),
                   "text": s[36:56].strip(), "ictype": int(s[56:58] or 0) if s[56:58].strip() else 0,
                   "numstp": int(s[58:63]) if s[58:63].strip() else 0,
                   "analysis": s[63:73].strip(), "format": int(s[73:75]) if s[73:75].strip() else None,
                   "name": "", "components": [], "data": {}}
            if blk["format"] in (2, 3):
                raise ValueError(f"{path}: binary result block is not supported")
            i += 1
            last = None
            while i < n and not lines[i].startswith(" -3"):
                rec = lines[i]
                tag = rec[:3]
                if tag == " -4":
                    parts = rec[3:].split()
                    blk["name"] = parts[0] if parts else ""
                elif tag == " -5":
                    parts = rec[3:].split()
                    if parts:
                        blk["components"].append(parts[0])
                elif tag == " -1":
                    ncomp = sum(1 for c in blk["components"] if c.upper() != "ALL") or None
                    w = _id_width(rec, None, blk["format"])
                    nid = int(rec[3:3 + w])
                    vals = _floats12(rec[3 + w:])
                    blk["data"][nid] = vals[:ncomp] if ncomp and len(vals) > ncomp else vals
                    last = nid
                elif tag == " -2" and last is not None:
                    w = _id_width(rec, None, blk["format"])
                    blk["data"][last] = blk["data"][last] + _floats12(rec[3 + w:])
                i += 1
            blk["data"] = {k: np.asarray(v, float) for k, v in blk["data"].items()}
            blocks.append(blk)
            i += 1
            continue
        i += 1
    return {"nodes": nodes, "elements": elements, "blocks": blocks, "header": header}


def read_calculix_dat_frequencies(path: str) -> np.ndarray:
    """Frequencies [Hz] of the last ``E I G E N V A L U E   O U T P U T`` table of a ``.dat``."""
    lines = _open_text(path)
    tables, cur = [], None
    for line in lines:
        up = line.upper()
        if "E I G E N V A L U E" in up and "O U T P U T" in up:
            cur = []
            tables.append(cur)
            continue
        if cur is None:
            continue
        parts = line.split()
        if not parts:
            if cur:
                cur = None               # blank line after the rows ends the table
            continue
        if parts[0].isdigit() and len(parts) >= 2:
            try:
                lam = _f(parts[1])
            except ValueError:
                continue
            cur.append(math.sqrt(max(lam, 0.0)) / (2 * math.pi))
        elif cur:
            cur = None
    tables = [t for t in tables if t]
    if not tables:
        raise ValueError(f"{path}: no 'E I G E N V A L U E   O U T P U T' table")
    return np.asarray(tables[-1], float)


def _mode_blocks(frd: dict, step: int | None) -> list[dict]:
    disp = [b for b in frd["blocks"] if b["name"].upper() == "DISP"]
    if not disp:
        raise ValueError("no DISP result blocks in the .frd (add *NODE FILE / U to the *FREQUENCY step)")
    freq = [b for b in disp if b["ictype"] == 2] or disp
    groups, cur = [], []
    for b in freq:                       # a new *FREQUENCY step restarts at a lower frequency
        if cur and b["value"] < cur[-1]["value"] * (1 - 1e-6) - 1e-12:
            groups.append(cur)
            cur = []
        cur.append(b)
    if cur:
        groups.append(cur)
    return groups[-1 if step is None else step]


def read_calculix_modal(path: str, nmodes: int | None = None, frequency_file: str | None = None,
                        skip_below_hz: float = 1e-3, gen_mass=1.0, length_scale: float = 1.0,
                        mass_scale: float = 1.0, damping=0.0, transform=None, origin=(0.0, 0.0, 0.0),
                        step: int | None = None, frequencies=None, frd: dict | None = None) -> ModalModel:
    """Modal model from a CalculiX ``.frd`` (+ optional ``.dat`` frequencies).

    ``step``: index of the eigenmode group (default: the last ``*FREQUENCY``
    step).  ``frequency_file``: a ``.dat`` whose eigenvalue table overrides
    the E12.5 header frequencies (matched by mode order); default: the
    ``.dat`` next to the ``.frd`` if it exists.  Modes below
    ``skip_below_hz`` (rigid-body modes) are skipped.  Only nodes with
    results are kept; rotations are zero (``ccx`` writes translations).
    ``frd``: the already parsed :func:`read_frd` dictionary of ``path``.
    """
    frd = read_frd(path) if frd is None else frd
    blocks = _mode_blocks(frd, step)
    f_frd = np.array([b["value"] for b in blocks], float)
    if frequencies is not None:
        f_all = np.asarray(frequencies, float)
    else:
        if frequency_file is None:
            cand = os.path.splitext(path)[0] + ".dat"
            frequency_file = cand if os.path.exists(cand) else None
        f_all = f_frd
        if frequency_file:
            fd = read_calculix_dat_frequencies(frequency_file)
            n = len(f_frd)
            for lo in (0, len(fd) - n):            # same modes first, else the last n of the table
                if len(fd) >= n and lo >= 0 and np.allclose(fd[lo:lo + n], f_frd, rtol=2e-4, atol=1e-6):
                    f_all = fd[lo:lo + n]
                    break
            else:
                warnings.warn(f"{frequency_file}: eigenvalue table does not match the .frd frequencies; "
                              "the .frd values are used", stacklevel=2)
    keep = [j for j in range(len(blocks)) if f_all[j] >= skip_below_hz]
    if nmodes is not None:
        keep = keep[:int(nmodes)]
    if not keep:
        raise ValueError(f"{path}: no elastic modes above {skip_below_hz} Hz")
    ids = sorted(set(blocks[keep[0]]["data"]) & set(frd["nodes"]))
    if not ids:
        raise ValueError(f"{path}: the DISP blocks hold no node of the node block")
    coords = np.array([frd["nodes"][k] for k in ids])
    sh = np.zeros((len(ids), 6, len(keep)))
    for m, j in enumerate(keep):
        d = blocks[j]["data"]
        for r, k in enumerate(ids):
            v = d.get(k)
            if v is not None:
                nv = min(len(v), 6)
                sh[r, :nv, m] = v[:nv]
    gm = np.asarray(gen_mass, float)
    if gm.ndim and gm.size == len(blocks):
        gm = gm[keep]
    labels = [f"ccx mode {j + 1}" for j in keep]
    mm = _finish(np.asarray(ids), coords, sh, f_all[keep], gm, damping, length_scale, mass_scale,
                 transform, origin, labels, os.path.basename(path))
    mm.node_ids = np.asarray(ids)
    return mm


# ----------------------------------------------------------------- visualisation
_SOLID_FACES = {
    "he": ([(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)], []),
    "te": ([], [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]),
    "pe": ([(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)], [(0, 1, 2), (3, 4, 5)]),
}


def frd_surface(frd: dict, node_ids) -> dict:
    """Display mesh ``{"tris", "edges"}`` (0-based into ``node_ids``) of the
    element block: boundary faces of solids, corners of 2-D elements, beams
    as edges.  Elements referring to nodes outside ``node_ids`` are dropped."""
    index = {int(k): i for i, k in enumerate(node_ids)}
    faces: dict[tuple, tuple] = {}
    tris, edges = [], []
    for typ, conn in frd["elements"].values():
        name = FRD_TYPES.get(typ, "")
        c = [index.get(int(x)) for x in conn[:_NNODE.get(typ, len(conn))]]
        if name[:2] in _SOLID_FACES:
            quads, trs = _SOLID_FACES[name[:2]]
            for f in quads + trs:
                if all(c[q] is not None for q in f):
                    key = tuple(sorted(c[q] for q in f))
                    if key in faces:
                        faces[key] = None
                    else:
                        faces[key] = tuple(c[q] for q in f)
        elif name in ("tr3", "tr6") and all(x is not None for x in c[:3]):
            tris.append(c[:3])
        elif name in ("qu4", "qu8") and all(x is not None for x in c[:4]):
            tris += [[c[0], c[1], c[2]], [c[0], c[2], c[3]]]
        elif name in ("be2", "be3") and c[0] is not None and c[1] is not None:
            edges.append([c[0], c[1]])
    for f in faces.values():
        if f is None:
            continue
        if len(f) == 3:
            tris.append(list(f))
        else:
            tris += [[f[0], f[1], f[2]], [f[0], f[2], f[3]]]
    return {"tris": np.asarray(tris, int).reshape(-1, 3), "edges": np.asarray(edges, int).reshape(-1, 2)}


# ----------------------------------------------------------------- sample writers
def _e(v: float) -> str:
    s = f"{v:12.5E}"
    return s if len(s) == 12 else f"{v:.4E}".rjust(12)


def write_frd(modal: ModalModel, path: str, elements=None, long: bool = True, node_ids=None,
              unit_length: float = 1.0, unit_mass: float = 1.0, static_step: bool = False,
              rigid_modes: int = 0) -> None:
    """Write ``modal`` as a ccx-style ASCII ``.frd`` (node, element, DISP blocks).

    ``elements``: ``[(type_code, [node ids...]), ...]`` (frd type codes);
    ``unit_length`` / ``unit_mass``: SI value of the solver units (shapes are
    rescaled so that unit generalised mass holds in solver units).
    ``static_step`` adds a preceding static DISP block, ``rigid_modes`` zero-
    frequency modes - both must be filtered by the reader.
    """
    ids = np.arange(1, modal.nodes.shape[0] + 1) if node_ids is None else np.asarray(node_ids, int)
    w = 10 if long else 5
    fmt = 1 if long else 0
    xyz = modal.nodes / unit_length
    gm = modal.gen_mass / (unit_mass * unit_length ** 2)
    sh = modal.shapes[:, :3] / unit_length / np.sqrt(gm)[None, None, :]
    fhz = modal.omega / (2 * np.pi)
    out = ["    1C" + os.path.splitext(os.path.basename(path))[0][:6].ljust(6),
           "    1UDATE              26.november.2026",
           "    1UTIME              12:00:00",
           "    1UHOST",
           "    1UPGM               CalculiX",
           "    1UVERSION           Version 2.21",
           "    1UCOMPILETIME       sample",
           "    1UDIR",
           "    1UDBN",
           "    1UMAT    1          ALUMINIUM"]
    out.append(f"    2C{'':18s}{len(ids):12d}{'':37s}{fmt:1d}")
    for k, p in zip(ids, xyz):
        out.append(f" -1{k:{w}d}" + "".join(_e(x) for x in p))
    out.append(" -3")
    if elements:
        out.append(f"    3C{'':18s}{len(elements):12d}{'':37s}{fmt:1d}")
        for e, (typ, conn) in enumerate(elements, 1):
            out.append(f" -1{e:{w}d}{typ:5d}{0:5d}{1:5d}")
            per = 10 if long else 15
            conn = list(conn)
            for j in range(0, len(conn), per):
                out.append(" -2" + "".join(f"{int(c):{w}d}" for c in conn[j:j + per]))
        out.append(" -3")

    def block(value, numstp, ictype, analysis, data):
        out.append(f"  100CL{100 + numstp:5d}{_e(value)}{len(ids):12d}{'':20s}{ictype:2d}{numstp:5d}"
                   f"{analysis:10s}{fmt:2d}")
        out.append(" -4  DISP        4    1")
        for c in ("D1", "D2", "D3"):
            out.append(f" -5  {c:8s}{1:4d}{2:5d}{c[1]:>5s}{0:5d}")
        out.append(" -5  ALL         1    2    0    0    1ALL")
        for k, v in zip(ids, data):
            out.append(f" -1{k:{w}d}" + "".join(_e(x) for x in v))
        out.append(" -3")

    nstp = 0
    out.append("    1PSTEP" + f"{1:26d}{1:12d}{1:12d}")
    if static_step:
        nstp += 1
        block(1.0, nstp, 0, "", 1e-3 * np.ones((len(ids), 3)))
    rng = np.random.default_rng(1)
    for r in range(rigid_modes):
        nstp += 1
        out.append(f"    1PMODE{nstp:26d}")
        v = np.zeros((len(ids), 3))
        v[:, r % 3] = 1.0 + 0 * rng.random()
        block(0.0, nstp, 2, "MODAL", v)
    for m in range(sh.shape[2]):
        nstp += 1
        out.append(f"    1PMODE{m + 1 + rigid_modes:26d}")
        block(fhz[m], nstp, 2, "MODAL", sh[:, :, m])
    out.append("9999")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")


def write_calculix_dat(modal: ModalModel, path: str, rigid_modes: int = 0) -> None:
    """``.dat`` with the ccx eigenvalue table (lambda, rad/time, cycles/time)."""
    w = np.r_[np.zeros(rigid_modes), modal.omega]
    rows = ["", "     E I G E N V A L U E   O U T P U T", "",
            " MODE NO    EIGENVALUE                       FREQUENCY   ",
            "                                     REAL PART            IMAGINARY PART",
            "                             (RAD/TIME)      (CYCLES/TIME     (RAD/TIME)", ""]
    for j, om in enumerate(w, 1):
        rows.append(f"{j:7d}   {om * om:14.7E}   {om:14.7E}   {om / (2 * np.pi):14.7E}   {0.0:14.7E}")
    rows.append("")
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(rows) + "\n")


def write_calculix_inp(nodes: np.ndarray, quads, path: str, *, thickness: float, E: float, nu: float,
                       rho: float, clamped_nodes, nmodes: int = 6) -> None:
    """Shell (S4) plate deck for ``ccx`` with a ``*FREQUENCY`` step writing ``U`` to the ``.frd``.

    Written for users who want to run the solver themselves; FLAAP does not
    call ``ccx``.
    """
    L = ["** FLAAP sample: S4 plate, *FREQUENCY -> plate.frd / plate.dat", "*NODE, NSET=NALL"]
    L += [f"{i + 1}, {x:.8g}, {y:.8g}, {z:.8g}" for i, (x, y, z) in enumerate(nodes)]
    L.append("*ELEMENT, TYPE=S4, ELSET=EALL")
    L += [f"{e + 1}, " + ", ".join(str(int(c) + 1) for c in q) for e, q in enumerate(quads)]
    L.append("*NSET, NSET=ROOT")
    cl = [int(c) + 1 for c in clamped_nodes]
    L += [", ".join(str(c) for c in cl[j:j + 8]) for j in range(0, len(cl), 8)]
    L += ["*MATERIAL, NAME=ALU", "*ELASTIC", f"{E:.6g}, {nu:.6g}", "*DENSITY", f"{rho:.6g}",
          "*SHELL SECTION, ELSET=EALL, MATERIAL=ALU", f"{thickness:.6g}",
          "*BOUNDARY", "ROOT, 1, 6",
          "*STEP", "*FREQUENCY", f"{nmodes}", "*NODE FILE", "U", "*END STEP"]
    with open(path, "w", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
