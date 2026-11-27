"""OpenRadioss / Radioss ``/EIG`` eigen results -> :class:`ModalModel`.

:mod:`flaap.io.openradioss` reads a *starter deck* and solves the modes with
FLAAP's own FE.  This module instead reads normal modes computed by the
Radioss engine itself (``/EIG`` in the starter + ``/IMPL/LINEAR`` in the
engine input) and post-converted with ``anim_to_vtk``.

What Radioss writes for ``/EIG`` (from public documentation and source)
-----------------------------------------------------------------------
* Input: ``/EIG/eig_ID/unit_ID`` (title; ``grnd_ID grnd_bc Trarot Ifile``;
  ``Nmod Inorm Cutfreq Freqmin Nbloc``; ``Incv Niter Ipri Tol``; optional
  external modes file name).  The engine needs ``/IMPL/LINEAR``.  ``Inorm =
  0`` (default) normalises the eigenvectors to *unit generalised mass*,
  ``Inorm = 1`` to a unit largest displacement.  Modes are computed with
  ARPACK; frequencies below ``Freqmin`` (default 0.001 Hz) are dropped.
  [1]
* Mode shapes: one animation file per mode, ``<run>A001`` = mode 1, ... ([1],
  HyperView note).  The animation header carries the mode's frequency in
  the "time" slot, labelled ``Frequency=`` / ``Eigenmode`` (``Static mode``
  and 0 for the static modes of a flexible-body ``/EIG``) -- ``genani.F``
  [2].  Nodal vectors requested by ``/ANIM/VECT/DISP`` and
  ``/ANIM/VECT/DROT`` are named ``Displacement`` and ``Rotational DOF``
  (``genani.F``, ``anim_dcod_key_0.F`` [2]).
* ``anim_to_vtk <run>A001 > <run>A001.vtk`` converts an animation file to
  legacy ASCII VTK [3]: ``# vtk DataFile Version 3.0`` / ``vtk output`` /
  ``ASCII`` / ``DATASET UNSTRUCTURED_GRID``, ``FIELD FieldData 2`` with
  ``TIME 1 1 double`` (the header "time", i.e. the frequency of an
  eigenmode file) and ``CYCLE 1 1 int``, ``POINTS n float`` (model
  coordinates), ``CELLS`` / ``CELL_TYPES`` (0-based point indices),
  ``POINT_DATA n`` with ``SCALARS NODE_ID int 1`` (Radioss node ids), the
  nodal scalars and ``VECTORS <name> float`` (blanks of the Radioss name
  replaced by ``_``: ``Displacement``, ``Rotational_DOF``), then
  ``CELL_DATA`` (``ELEMENT_ID``, ``PART_ID``, ``EROSION_STATUS`` ...).
  Numbers are written by C++ ``operator<<`` (6 significant digits,
  single-precision data).

Important: in the public OpenRadioss code base the engine eigen solver is
*not* included -- ``engine/stub/eig.F``, ``eig1.F``, ``eigp.F``,
``eigcond.F`` and ``imp_k_eig.F`` are stubs that stop the run, and the call
of ``IMP_EIGSOL`` in ``engine/source/engine/resol.F`` is compiled only with
the ``DNC`` macro, which the public build does not define [2].  The starter
accepts and echoes ``/EIG`` (``hm_read_eig.F``), but an open-source engine
produces no modes.  Eigen results therefore come from the commercial
(Altair) Radioss engine; their animation files are converted as above.  For
the same reason the layout of the engine-listing (``*_0001.out``) eigenvalue
table is not available in public source or documentation: the frequencies
are taken from the VTK ``TIME`` field, and :func:`read_radioss_frequencies`
parses a listing / frequency table heuristically (header with ``FREQ`` or
``EIGEN``, rows ``mode  value ...``; or a plain two-column table).

Sources
  [1] Altair Radioss Reference Guide, ``/EIG`` (Starter):
      https://help.altair.com/hwsolvers/rad/topics/solvers/rad/eig_starter_r.htm
  [2] OpenRadioss source, preserved as OpenCourant (AGPL-3.0, full history):
      https://github.com/OpenCourant/OpenCourant -- ``engine/stub/eig.F``,
      ``engine/source/engine/resol.F``, ``engine/source/output/anim/generate/
      genani.F`` and ``ani_txt.F``, ``engine/source/output/anim/reader/
      anim_dcod_key_0.F``, ``starter/source/general_controls/computation/
      hm_read_eig.F``.  (github.com/OpenRadioss is no longer public.)
  [3] anim_to_vtk converter: https://github.com/OpenCourant/Tools --
      ``output_converters/anim_to_vtk/src/anim_to_vtk.cpp``.

Units and normalisation
-----------------------
VTK coordinates, displacements and frequencies are in the solver unit
system (``/BEGIN`` *work* units; a frequency is in cycles per time unit,
e.g. 1/ms = kHz).  They are converted to SI with the starter deck units
(``deck=``) or ``units=(mass, length, time)`` (default SI).  Translations
are multiplied by the length scale, rotations are unchanged (same physical
motion), the generalised mass of unit-mass-normalised shapes becomes
``mass_scale * length_scale**2``.

``gen_mass``:

* ``"auto"`` (default): with a starter deck, ``phi^T M phi`` with FLAAP's
  consistent mass matrix of the deck (:func:`flaap.io.openradioss.radioss_mesh`
  + :func:`flaap.structures.fem.build_fe_model`); this is valid for both
  ``Inorm`` choices.  Without a deck the shapes are assumed to be
  ``Inorm = 0`` (unit generalised mass in solver units).
* ``"unit"``: ``Inorm = 0`` assumed even when a deck is given.
* ``"deck"``: always from the deck mass matrix (error without deck).
* a number / per-mode array in solver units.

Rotations: from ``Rotational_DOF`` when written (``/ANIM/VECT/DROT``).
Otherwise zero, except that with a deck the rotations of solid-element
nodes are estimated from the translations (local least squares, as
:func:`radioss_modal`); shell nodes then keep zero rotations (a warning
recommends ``/ANIM/VECT/DROT``).

Synthetic writer
----------------
:func:`write_radioss_eigen_vtk` writes a :class:`ModalModel` as
``<run>A###.vtk`` files in the ``anim_to_vtk`` layout above.  It exists to
make test / example files (``examples/openradioss/eig/``) because no
Radioss eigen results are available here; it does not imitate the engine
listing.
"""
from __future__ import annotations

import glob
import os
import re
import warnings

import numpy as np

from ..structures.modal import ModalModel
from .openradioss import (UNIT_LENGTH, UNIT_MASS, UNIT_TIME, RadiossDeck, _bc_dofs, _unit,
                          radioss_mesh, read_radioss)

__all__ = ["read_radioss_vtk", "read_radioss_frequencies", "find_radioss_mode_files",
           "read_radioss_eigen", "write_radioss_eigen_vtk"]

_NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?|[-+]?(?:nan|inf)"
_NUM_RE = re.compile(_NUM, re.I)
_MODE_RE = re.compile(r"A(\d{3,})(?:\.vtk)?$", re.I)


# ===================================================================== VTK
def read_radioss_vtk(path: str) -> dict:
    """Parse a legacy ASCII VTK file written by ``anim_to_vtk``.

    Returns ``{"time", "cycle", "points" (n,3), "node_ids" (n,), "point_scalars"
    {name: (n,)}, "vectors" {name: (n,3)}, "cells" [list of 0-based index
    lists], "cell_types" (nc,), "cell_scalars" {name: (nc,)}, "title"}``.
    Generic legacy-VTK reader for UNSTRUCTURED_GRID with FIELD, SCALARS,
    VECTORS, NORMALS, TENSORS and LOOKUP_TABLE sections.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    if len(lines) < 4 or not lines[0].lstrip().lower().startswith("# vtk datafile"):
        raise ValueError(f"{path}: not a legacy VTK file")
    if lines[2].strip().upper() != "ASCII":
        raise ValueError(f"{path}: only ASCII legacy VTK is supported (got {lines[2].strip()})")
    tok = " ".join(lines[3:]).split()
    out = {"title": lines[1].strip(), "time": None, "cycle": None, "points": None, "node_ids": None,
           "point_scalars": {}, "vectors": {}, "cells": [], "cell_types": np.zeros(0, int),
           "cell_scalars": {}}
    i, n = 0, len(tok)
    assoc, npt, ncell = None, 0, 0

    def take(k, dtype=float):
        nonlocal i
        if i + k > n:
            raise ValueError(f"{path}: unexpected end of file")
        a = np.array(tok[i:i + k], dtype=float)
        i += k
        return a.astype(dtype) if dtype is not float else a

    while i < n:
        kw = tok[i].upper()
        i += 1
        if kw == "DATASET":
            if tok[i].upper() != "UNSTRUCTURED_GRID":
                warnings.warn(f"{path}: DATASET {tok[i]} (UNSTRUCTURED_GRID expected)")
            i += 1
        elif kw == "FIELD":
            nf = int(tok[i + 1])
            i += 2
            for _ in range(nf):
                name, nc, nt = tok[i], int(tok[i + 1]), int(tok[i + 2])
                i += 4
                vals = take(nc * nt)
                key = name.upper()
                if assoc is None and key == "TIME" and vals.size:
                    out["time"] = float(vals[0])
                elif assoc is None and key == "CYCLE" and vals.size:
                    out["cycle"] = int(vals[0])
                elif assoc == "point" and nc == 1:
                    out["point_scalars"][name] = vals
                elif assoc == "point" and nc == 3:
                    out["vectors"][name] = vals.reshape(-1, 3)
                elif assoc == "cell" and nc == 1:
                    out["cell_scalars"][name] = vals
        elif kw == "POINTS":
            npt = int(tok[i])
            i += 2
            out["points"] = take(3 * npt).reshape(npt, 3)
        elif kw == "CELLS":
            nc, size = int(tok[i]), int(tok[i + 1])
            i += 2
            flat = take(size, int)
            cells, p = [], 0
            while p < flat.size:
                k = int(flat[p])
                cells.append(flat[p + 1:p + 1 + k].tolist())
                p += k + 1
            if len(cells) != nc:
                raise ValueError(f"{path}: CELLS count mismatch")
            out["cells"] = cells
        elif kw == "CELL_TYPES":
            k = int(tok[i])
            i += 1
            out["cell_types"] = take(k, int)
        elif kw == "POINT_DATA":
            assoc, npt = "point", int(tok[i])
            i += 1
        elif kw == "CELL_DATA":
            assoc, ncell = "cell", int(tok[i])
            i += 1
        elif kw == "SCALARS":
            name = tok[i]
            i += 2
            ncomp = 1
            if i < n and tok[i].isdigit():
                ncomp = int(tok[i])
                i += 1
            if i < n and tok[i].upper() == "LOOKUP_TABLE":
                i += 2
            cnt = (npt if assoc == "point" else ncell) * ncomp
            vals = take(cnt)
            if assoc == "point":
                if name.upper() == "NODE_ID":
                    out["node_ids"] = vals.astype(np.int64)
                else:
                    out["point_scalars"][name] = vals
            else:
                out["cell_scalars"][name] = vals
        elif kw in ("VECTORS", "NORMALS"):
            name = tok[i]
            i += 2
            cnt = npt if assoc == "point" else ncell
            vals = take(3 * cnt).reshape(cnt, 3)
            if assoc == "point":
                out["vectors"][name] = vals
        elif kw == "TENSORS":
            i += 2
            take(9 * (npt if assoc == "point" else ncell))
        elif kw == "LOOKUP_TABLE":            # a colour table: name size + 4*size values
            k = int(tok[i + 1])
            i += 2
            take(4 * k)
        elif kw == "METADATA":                # VTK >= 5: skip to the blank-terminated block
            while i < n and tok[i].upper() not in ("SCALARS", "VECTORS", "POINT_DATA", "CELL_DATA",
                                                   "FIELD", "TENSORS", "NORMALS"):
                i += 1
        else:
            raise ValueError(f"{path}: unexpected token {tok[i - 1]!r}")
    if out["points"] is None:
        raise ValueError(f"{path}: no POINTS section")
    if out["node_ids"] is None:
        warnings.warn(f"{path}: no NODE_ID scalars; point index + 1 used as node id")
        out["node_ids"] = np.arange(1, out["points"].shape[0] + 1)
    return out


def _vector(d: dict, *names):
    low = {k.lower(): v for k, v in d["vectors"].items()}
    for nm in names:
        if nm.lower() in low:
            return low[nm.lower()]
    return None


def find_radioss_mode_files(path) -> list[str]:
    """Mode VTK files from a directory, one ``<run>A###.vtk`` file (its
    siblings with the same run prefix are added) or a list of paths.
    Sorted by the animation number (``A001`` = mode 1)."""
    if isinstance(path, (list, tuple)):
        files = [os.fspath(p) for p in path]
    else:
        p = os.fspath(path)
        if os.path.isdir(p):
            files = glob.glob(os.path.join(p, "*.vtk")) + glob.glob(os.path.join(p, "*.VTK"))
            num = [f for f in files if _MODE_RE.search(os.path.basename(f))]
            files = num or files
        elif p.lower().endswith(".vtk"):
            m = _MODE_RE.search(os.path.basename(p))
            if m:
                stem = os.path.basename(p)[:m.start()]
                files = [f for f in glob.glob(os.path.join(os.path.dirname(p) or ".", glob.escape(stem)
                                                           + "A*.vtk"))
                         if _MODE_RE.fullmatch(os.path.basename(f)[len(stem):])]
            else:
                files = [p]
        else:
            raise ValueError(f"{p}: not a VTK file or directory")
    files = sorted(set(files), key=lambda f: (_mode_no(f) is None, _mode_no(f) or 0, f))
    if not files:
        raise ValueError(f"no VTK mode files found in {path}")
    return files


def _mode_no(f: str) -> int | None:
    m = _MODE_RE.search(os.path.basename(f))
    return int(m.group(1)) if m else None


# ===================================================================== listing
def read_radioss_frequencies(path: str) -> dict[int, float]:
    """``{mode: frequency}`` from an engine listing or a frequency table.

    Heuristic (the Radioss eigen listing layout is not public): the first
    table whose header line contains ``FREQ`` -- the value column is the
    number nearest to the header word -- with rows ``mode value ...``; else a
    header containing ``EIGEN`` (eigenvalues ``omega^2`` -> ``sqrt / 2 pi``);
    else rows of exactly two numbers ``mode frequency``.  Lines starting
    with ``#`` are comments.  Values are in the solver time unit.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [ln.rstrip("\n\r").expandtabs(8) for ln in fh]
    for word, conv in (("FREQ", None), ("EIGEN", "eig")):
        res = _header_table(lines, word, conv)
        if res:
            return res
    res = {}
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        t = s.replace(",", " ").replace(";", " ").split()
        if len(t) == 2 and t[0].isdigit() and _NUM_RE.fullmatch(t[1]):
            res.setdefault(int(t[0]), float(t[1].replace("D", "E").replace("d", "e")))
    if not res:
        raise ValueError(f"{path}: no frequency table found")
    return res


def _header_table(lines, word, conv):
    for h, head in enumerate(lines):
        up = head.upper()
        if word not in up:
            continue
        m = re.search(r"\S*" + word + r"\S*", up)
        center = (m.start() + m.end()) / 2
        res, started = {}, False
        for ln in lines[h + 1:]:
            s = ln.strip()
            if not s or set(s) <= set("-=*_ "):
                if started and not s:
                    break
                continue
            nums = list(_NUM_RE.finditer(ln))
            t0 = s.split()[0]
            if not (t0.isdigit() and len(nums) >= 2):
                if started:
                    break
                continue
            started = True
            best = min(nums[1:], key=lambda q: abs((q.start() + q.end()) / 2 - center))
            v = float(best.group().replace("D", "E").replace("d", "e"))
            if conv == "eig":
                v = np.sqrt(max(v, 0.0)) / (2 * np.pi)
            res.setdefault(int(t0), v)
        if res:
            return res
    return None


# ===================================================================== reader
def _unit_scales(deck, units):
    if deck is not None:
        return deck.scale["mass"], deck.scale["length"], deck.scale["time"], deck.units
    mu, lu, tu = units if units is not None else ("kg", "m", "s")
    return (_unit(UNIT_MASS, mu, "mass"), _unit(UNIT_LENGTH, lu, "length"), _unit(UNIT_TIME, tu, "time"),
            (mu, lu, tu))


def _resolve_sources(source, vtk_files, listing, deck):
    """(vtk files, listing path or None, deck path / object or None)."""
    if isinstance(source, (list, tuple)):
        vtk_files = vtk_files or source
        source = None
    if source is not None:
        p = os.fspath(source)
        if os.path.isdir(p):
            vtk_files = vtk_files or p
            if listing is None:
                outs = sorted(glob.glob(os.path.join(p, "*_0001.out")))
                listing = outs[0] if outs else None
            if deck is None:
                rad = sorted(glob.glob(os.path.join(p, "*_0000.rad")))
                deck = rad[0] if rad else None
        elif p.lower().endswith(".vtk"):
            vtk_files = vtk_files or p
        else:
            listing = listing or p
            if vtk_files is None:
                d = os.path.dirname(p) or "."
                run = re.sub(r"_\d{4}\.out$", "", os.path.basename(p), flags=re.I)
                cand = [f for f in glob.glob(os.path.join(d, glob.escape(run) + "A*.vtk"))
                        if _MODE_RE.fullmatch(os.path.basename(f)[len(run):])]
                vtk_files = cand or d
    if vtk_files is None:
        raise ValueError("no mode shapes: pass the VTK mode files (vtk_files= or a directory)")
    return find_radioss_mode_files(vtk_files), listing, deck


def read_radioss_eigen(source=None, vtk_files=None, deck=None, units=None, *, listing=None,
                       nmodes: int | None = None, gen_mass="auto", frequencies=None,
                       skip_below_hz: float = 1e-3, damping=0.0, rotations: str = "auto",
                       name: str = "openradioss_eig") -> ModalModel:
    """Normal modes computed by the Radioss engine (``/EIG``) -> :class:`ModalModel` (SI).

    ``source``: a directory holding the ``<run>A###.vtk`` mode files (an
    ``*_0001.out`` listing and an ``*_0000.rad`` starter deck in it are used
    when present), one mode file (siblings ``<run>A###.vtk`` are added), a
    list of VTK files, or the engine listing (VTK files then found next to
    it, or given by ``vtk_files``).  ``deck``: starter deck path or
    :class:`RadiossDeck` -- gives SI node coordinates, the unit system,
    generalised masses and the viewer mesh.  ``units``: ``(mass, length,
    time)`` of the run (default: the deck ``/BEGIN`` units, else SI; with a
    deck path it overrides ``/BEGIN`` as in :func:`read_radioss`).  ``listing``:
    optional listing / frequency table (solver units) overriding the VTK
    ``TIME`` frequencies; ``frequencies`` [Hz, SI] overrides both.  Modes
    below ``skip_below_hz`` (rigid-body / static modes, ``TIME = 0``) are
    skipped; ``nmodes`` keeps the first ones.  See the module docstring for
    ``gen_mass`` and ``rotations`` (``"auto"`` | ``"zero"``).

    The result has ``node_ids`` (Radioss ids), ``labels`` (``mode<k>``,
    ``k`` = animation number), ``source_files`` and, with a deck, ``fe``
    (FLAAP FE model of the deck, ``fe.elements["mesh"]``).
    """
    files, listing, deck = _resolve_sources(source, vtk_files, listing, deck)
    if isinstance(deck, (str, os.PathLike)):
        deck = read_radioss(os.fspath(deck), units=units)
    if deck is not None and not isinstance(deck, RadiossDeck):
        raise TypeError("deck must be a path or a RadiossDeck")
    sm, sl, st, unit_names = _unit_scales(deck, units)

    data = [read_radioss_vtk(f) for f in files]
    lst = None
    if listing is not None:
        try:
            lst = read_radioss_frequencies(os.fspath(listing))
        except ValueError as e:
            warnings.warn(f"listing ignored: {e}")
    fr_override = None if frequencies is None else np.atleast_1d(np.asarray(frequencies, float))

    modes = []                               # (label, f_hz, data)
    for c, (f, d) in enumerate(zip(files, data)):
        k = _mode_no(f) or c + 1
        if fr_override is not None:
            if c >= fr_override.size:
                break
            fhz = float(fr_override[c])
        elif lst is not None and k in lst:
            fhz = lst[k] / st
        elif d["time"] is not None:
            fhz = d["time"] / st
        else:
            raise ValueError(f"{f}: no frequency (no TIME field, no listing entry)")
        if fhz < skip_below_hz:
            continue
        modes.append((f"mode{k}", fhz, d))
    if nmodes is not None:
        modes = modes[:int(nmodes)]
    if not modes:
        raise ValueError("no elastic modes found in " + ", ".join(files))
    nm = len(modes)

    # node set and coordinates
    d0 = modes[0][2]
    if deck is not None:
        ids = deck.node_ids.copy()
        nodes = deck.nodes.copy()
    else:
        ids = np.asarray(d0["node_ids"], np.int64)
        nodes = np.asarray(d0["points"], float) * sl
    lut = {int(g): r for r, g in enumerate(ids)}
    sh = np.zeros((ids.size, 6, nm))
    have_rot = True
    for c, (_, _, d) in enumerate(modes):
        U = _vector(d, "Displacement", "DISPLACEMENT", "Displacements")
        if U is None:
            raise ValueError("mode file without 'Displacement' vectors (request /ANIM/VECT/DISP)")
        R = _vector(d, "Rotational_DOF", "Rotational DOF", "DROT")
        rows = np.array([lut.get(int(g), -1) for g in d["node_ids"]])
        ok = rows >= 0
        if deck is not None and c == 0:
            miss = ids.size - np.unique(rows[ok]).size
            if miss:
                warnings.warn(f"{miss} deck nodes have no mode data (zero shapes)")
        sh[rows[ok], :3, c] = U[ok] * sl
        if R is not None:
            sh[rows[ok], 3:, c] = R[ok]
        else:
            have_rot = False

    fe, mesh = None, None
    if deck is not None:
        try:
            from ..structures.fem import build_fe_model
            mesh = radioss_mesh(deck)
            fixed, _ = _bc_dofs(deck, mesh)
            fe = build_fe_model(mesh, fixed=fixed)
        except (KeyError, ValueError, NotImplementedError) as e:
            warnings.warn(f"starter deck FE model not built ({e}); no deck mass matrix")
    if not have_rot:
        if rotations == "auto" and mesh is not None:
            from ..structures.solid import rotations_from_translations
            shell = np.zeros(mesh.nn, bool)
            for b in mesh.blocks:
                if b.kind in ("quad4", "tri3"):
                    shell[b.conn.ravel()] = True
            solid = ~shell
            if fe is not None:                        # constrained nodes keep zero (as solve_modes)
                solid[np.unique(fe.fixed // 6)] = False
            if solid.any():
                rot = rotations_from_translations(mesh.nodes, mesh.conns(), sh[:, :3, :], solid)
                sh[solid, 3:, :] = rot[solid]
            if shell.any():
                warnings.warn("no 'Rotational_DOF' in the mode files: shell-node rotations set to zero "
                              "(add /ANIM/VECT/DROT to the engine input)")
        elif rotations not in ("auto", "zero"):
            raise ValueError("rotations must be 'auto' or 'zero'")

    # generalised masses
    gm_mode = gen_mass if isinstance(gen_mass, str) else None
    if gm_mode is None:
        gm = np.broadcast_to(np.asarray(gen_mass, float), (nm,)) * sm * sl ** 2
    elif gm_mode in ("auto", "deck"):
        if fe is not None:
            V = sh.reshape(-1, nm)
            gm = np.einsum("in,in->n", V, fe.M @ V)
        elif gm_mode == "deck":
            raise ValueError("gen_mass='deck' needs a starter deck (deck=)")
        else:
            gm = np.full(nm, sm * sl ** 2)
    elif gm_mode == "unit":
        gm = np.full(nm, sm * sl ** 2)
    else:
        raise ValueError("gen_mass must be 'auto', 'unit', 'deck' or numbers")

    zeta = np.broadcast_to(np.asarray(damping, float), (nm,)).copy()
    mm = ModalModel(nodes, sh, 2 * np.pi * np.array([m[1] for m in modes]), np.asarray(gm, float).copy(),
                    damping=zeta, fe=fe, name=(deck.title if deck is not None and deck.title else name),
                    labels=[m[0] for m in modes])
    mm.node_ids = ids
    mm.source_files = [f for f, d in zip(files, data) if any(d is m[2] for m in modes)]
    mm.units = unit_names
    return mm


# ===================================================================== writer
def _g(v) -> str:
    """C++ ``operator<<`` default formatting of a float (``%g``, 6 digits)."""
    return f"{float(np.float32(v)):g}"


def write_radioss_eigen_vtk(modal: ModalModel, directory: str, run: str = "run", deck=None,
                            units=None, node_ids=None, rotations: bool = True,
                            frequencies_solver=None) -> list[str]:
    """Write ``modal`` (SI) as ``<run>A001.vtk`` ... in the ``anim_to_vtk`` layout.

    Synthetic sample writer (tests / examples).  Shapes are rescaled to unit
    generalised mass in the solver unit system (``Inorm = 0``) and written
    in solver units (deck ``/BEGIN`` units or ``units``); the ``TIME`` field
    holds the frequency in cycles per solver time unit.  Cells (shells as
    4-node facets, solids as 8-node bricks / 4-node tetrahedra) come from
    the deck when given.  Returns the written paths.
    """
    if isinstance(deck, (str, os.PathLike)):
        deck = read_radioss(os.fspath(deck), units=units)
    sm, sl, st, _ = _unit_scales(deck, units)
    os.makedirs(directory, exist_ok=True)
    if node_ids is None:
        node_ids = getattr(modal, "node_ids", None)
    ids = np.asarray(node_ids if node_ids is not None else np.arange(1, modal.nodes.shape[0] + 1))
    X = modal.nodes / sl
    cells, ctypes, cids, cparts = [], [], [], []
    if deck is not None:
        lut = {int(g): r for r, g in enumerate(ids)}
        for (kind, part), (eids, conn, _) in sorted(deck.elements.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            for e, c in zip(eids, conn):
                c = [lut[int(g)] for g in c]
                if kind in ("quad4", "tri3"):
                    cells.append(c + [c[-1]] if kind == "tri3" else c)
                    ctypes.append(5 if kind == "tri3" else 9)
                elif kind in ("tet4", "tet10"):
                    cells.append(sorted(c[:4]))
                    ctypes.append(10)
                else:
                    cells.append(c[:8])
                    ctypes.append(12)
                cids.append(int(e))
                cparts.append(int(part))
    paths = []
    gm_solver = modal.gen_mass / (sm * sl ** 2)
    for m in range(modal.nmodes):
        scale = 1.0 / np.sqrt(gm_solver[m])
        U = modal.shapes[:, :3, m] * scale / sl
        R = modal.shapes[:, 3:, m] * scale
        f = (modal.freq_hz[m] * st if frequencies_solver is None else frequencies_solver[m])
        L = ["# vtk DataFile Version 3.0", "vtk output", "ASCII", "DATASET UNSTRUCTURED_GRID",
             "FIELD FieldData 2", "TIME 1 1 double", _g(f), "CYCLE 1 1 int", "0",
             f"POINTS {X.shape[0]} float"]
        L += [f"{_g(a)} {_g(b)} {_g(c)}" for a, b, c in X]
        L.append("")
        if cells:
            size = sum(len(c) + 1 for c in cells)
            L.append(f"CELLS {len(cells)} {size}")
            for c in cells:
                L.append(f"{len(c)} " + ("  " if len(c) == 8 else " ").join(str(x) for x in c))
            L.append("")
            L.append(f"CELL_TYPES {len(cells)}")
            L += [str(t) for t in ctypes]
        L.append("")
        L += [f"POINT_DATA {X.shape[0]}", "SCALARS NODE_ID int 1", "LOOKUP_TABLE default"]
        L += [str(int(g)) for g in ids]
        L.append("")
        L.append("VECTORS Displacement float")
        L += [f"{_g(a)} {_g(b)} {_g(c)}" for a, b, c in U]
        L.append("")
        if rotations:
            L.append("VECTORS Rotational_DOF float")
            L += [f"{_g(a)} {_g(b)} {_g(c)}" for a, b, c in R]
            L.append("")
        if cells:
            L += [f"CELL_DATA {len(cells)}", "SCALARS ELEMENT_ID int 1", "LOOKUP_TABLE default"]
            L += [str(e) for e in cids]
            L.append("")
            L += ["SCALARS PART_ID int 1", "LOOKUP_TABLE default"] + [str(p) for p in cparts]
            L.append("")
            L += ["SCALARS EROSION_STATUS int 1", "LOOKUP_TABLE default"] + ["0"] * len(cells)
            L.append("")
        p = os.path.join(directory, f"{run}A{m + 1:03d}.vtk")
        with open(p, "w", encoding="ascii", newline="\n") as fh:
            fh.write("\n".join(L) + "\n")
        paths.append(p)
    return paths
