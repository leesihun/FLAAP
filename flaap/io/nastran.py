"""MSC/NX Nastran readers: bulk-data GRIDs and SOL 103 ``.f06`` results.

Replaces ``read_f06.f90`` of the original program (which extracted the
T3 component of the eigenvectors from a fixed file layout) with a general
parser:

* :func:`read_grids` - ``GRID`` entries in small-field, large-field
  (``GRID*``) or free-field (comma) format, from a bulk data file or from
  the sorted bulk data echo of an ``.f06`` (``   34-        GRID ...``).
  Only basic coordinates (CP = 0 or blank) are supported.
* :func:`read_eigenvalues` - the ``R E A L   E I G E N V A L U E S`` table.
* :func:`read_eigenvectors` - ``R E A L   E I G E N V E C T O R   N O .``
  blocks (all six components per grid, continuation pages merged).
* :func:`read_f06_modal` - a :class:`~flaap.structures.modal.ModalModel`
  from the above (rigid-body modes skipped by frequency threshold,
  optional unit conversion and coordinate transform to the aero frame).
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from ..structures.modal import ModalModel

_ECHO = re.compile(r"^\s*\d+-\s{8}")


def _f(s: str) -> float:
    s = s.strip()
    if not s:
        return 0.0
    s = s.replace("D", "E").replace("d", "e")
    # Nastran shorthand exponents: 1.5-3 -> 1.5e-3
    m = re.match(r"^([+-]?\d*\.?\d*)([+-]\d+)$", s)
    if m and "e" not in s.lower():
        return float(m.group(1) + "e" + m.group(2))
    return float(s)


def _bulk_lines(path: str):
    """Yield bulk data lines (echo prefix removed, comments dropped)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            m = _ECHO.match(line)
            if m:
                line = line[m.end():]
            if line.startswith("$"):
                continue
            yield line


def read_grids(path: str) -> dict[int, np.ndarray]:
    """``{grid id: (x, y, z)}`` of all GRID entries."""
    lines = list(_bulk_lines(path))
    grids = {}
    i = 0
    while i < len(lines):
        ln = lines[i]
        key = ln[:8].replace(" ", "").upper()
        if "," in ln[:16] and ln.split(",")[0].strip().upper().startswith("GRID"):
            f = [t.strip() for t in ln.split(",")]
            if f[0].upper().rstrip("*") == "GRID":
                cp = f[2] if len(f) > 2 else ""
                if cp not in ("", "0"):
                    raise NotImplementedError("GRID with CP != 0")
                grids[int(f[1])] = np.array([_f(f[3]), _f(f[4]), _f(f[5]) if len(f) > 5 else 0.0])
        elif key == "GRID*":
            nxt = lines[i + 1] if i + 1 < len(lines) else ""
            gid = int(ln[8:24])
            cp = ln[24:40].strip()
            if cp not in ("", "0"):
                raise NotImplementedError("GRID with CP != 0")
            x, y = _f(ln[40:56]), _f(ln[56:72])
            z = _f(nxt[8:24]) if nxt.lstrip().startswith("*") else 0.0
            grids[gid] = np.array([x, y, z])
            i += 1
        elif key == "GRID":
            cp = ln[16:24].strip()
            if cp not in ("", "0"):
                raise NotImplementedError("GRID with CP != 0")
            grids[int(ln[8:16])] = np.array([_f(ln[24:32]), _f(ln[32:40]), _f(ln[40:48])])
        i += 1
    return grids


@dataclass
class Eigenvalues:
    mode: np.ndarray
    eigenvalue: np.ndarray
    omega: np.ndarray
    freq_hz: np.ndarray
    gen_mass: np.ndarray
    gen_stiff: np.ndarray


def read_eigenvalues(path: str) -> Eigenvalues:
    rows = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    for i, ln in enumerate(lines):
        if "R E A L   E I G E N V A L U E S" in ln:
            j = i + 1
            while j < len(lines) and not re.match(r"^\s*\d+\s+\d+\s+[-+\d.E]+", lines[j]):
                j += 1
            while j < len(lines):
                t = lines[j].split()
                try:
                    if len(t) < 7 or not t[0].isdigit():
                        raise ValueError
                    rows.append([float(x) for x in t[:7]])
                    j += 1
                except ValueError:
                    break
    if not rows:
        raise ValueError(f"{path}: no real eigenvalue table")
    a = np.array(rows)
    _, idx = np.unique(a[:, 0], return_index=True)
    a = a[np.sort(idx)]
    return Eigenvalues(a[:, 0].astype(int), a[:, 2], a[:, 3], a[:, 4], a[:, 5], a[:, 6])


def read_eigenvectors(path: str) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """``{mode: (grid ids (n,), values (n, 6))}`` (scalar points get T1 only)."""
    out: dict[int, dict[int, np.ndarray]] = {}
    mode = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            m = re.search(r"R E A L   E I G E N V E C T O R   N O \.\s+(\d+)", ln)
            if m:
                mode = int(m.group(1))
                out.setdefault(mode, {})
                continue
            if mode is None:
                continue
            if ln.startswith("1") and ("PAGE" in ln or "MSC" in ln or "NASTRAN" in ln.upper()):
                continue
            t = ln.split()
            if len(t) >= 3 and t[0].isdigit() and t[1] in ("G", "S"):
                vals = np.zeros(6)
                v = [float(x) for x in t[2:8]]
                vals[:len(v)] = v
                out[mode][int(t[0])] = vals
            elif t and not t[0].isdigit() and "POINT" not in ln and "CYCLES" not in ln \
                    and "EIGENVALUE" not in ln and "SUBCASE" not in ln and out[mode]:
                # any other text after the block ends the current eigenvector
                if not ln.strip().startswith(("0", "1")):
                    mode = None
    res = {}
    for k, d in out.items():
        ids = np.array(sorted(d))
        res[k] = (ids, np.array([d[i] for i in ids]))
    return res


def read_f06_modal(f06: str, bulk: str | None = None, nmodes: int | None = None,
                   skip_below_hz: float = 1e-3, length_scale: float = 1.0, mass_scale: float = 1.0,
                   transform: np.ndarray | None = None, origin=(0.0, 0.0, 0.0),
                   grid_ids=None) -> ModalModel:
    """Modal model from a SOL 103 ``.f06`` (+ GRIDs from ``bulk`` or the f06 echo).

    ``length_scale`` / ``mass_scale`` convert the Nastran units to SI (e.g.
    mm-tonne: ``length_scale=1e-3``, ``mass_scale=1e3``); mode shapes are
    rescaled to stay consistent with the generalized masses.
    ``transform`` (3x3) maps Nastran basic axes to the aero frame (applied
    to coordinates, translations and rotations); ``origin`` is subtracted
    from the coordinates after scaling.
    """
    ev = read_eigenvalues(f06)
    vecs = read_eigenvectors(f06)
    grids = read_grids(bulk or f06)
    keep = [i for i, f in enumerate(ev.freq_hz) if f >= skip_below_hz and int(ev.mode[i]) in vecs]
    if nmodes is not None:
        keep = keep[:nmodes]
    if not keep:
        raise ValueError("no elastic modes with eigenvectors found")
    ids0 = vecs[int(ev.mode[keep[0]])][0]
    ids = np.asarray(grid_ids if grid_ids is not None else [g for g in ids0 if g in grids])
    T = np.eye(3) if transform is None else np.asarray(transform, float)
    nodes = (np.array([grids[g] for g in ids]) * length_scale) @ T.T - np.asarray(origin, float)
    nm = len(keep)
    sh = np.zeros((ids.size, 6, nm))
    for c, i in enumerate(keep):
        gid, v = vecs[int(ev.mode[i])]
        pos = {g: r for r, g in enumerate(gid)}
        rows = np.array([v[pos[g]] for g in ids])
        sh[:, :3, c] = rows[:, :3] @ T.T * length_scale
        sh[:, 3:, c] = rows[:, 3:] @ T.T
    gm = ev.gen_mass[keep] * mass_scale * length_scale ** 2
    return ModalModel(nodes, sh, ev.omega[keep], gm, name="nastran",
                      labels=[f"mode{int(ev.mode[i])}" for i in keep])
