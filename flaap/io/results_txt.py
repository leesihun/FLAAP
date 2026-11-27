"""FLAAP result text files - the human-readable output contract.

Every file the launcher writes to an output directory (``results.txt``,
``modes.txt``, ``flutter_pk.txt``, ``timehist_*.txt`` ...) uses the same
grammar as the input configs (AI-CAE4ALL style), plus table blocks::

    % comment line
    kind          flutter              # inline comment
    method        pk
    names         mode 1, mode 2, mode 3
    meta.V        100.0                # nested keys use dotted paths

    table points
         kind  branch  velocity  freq_hz         k  density
      flutter       1  153.2855 11.24785 0.4216087    1.225
    end

Rules

* ``key value`` per line; the first whitespace splits key and value.
* ``%`` starts a comment line, ``#`` (after a blank) an inline comment.
* A value with commas is a list; ``[x]`` is a one-element list, ``[]`` empty.
* ``true`` / ``false`` are booleans, ``none`` is a missing value; numbers are
  ``int`` / ``float``.  A string that would otherwise be read as a number,
  boolean or list is written in double quotes.
* ``table <name> [records|matrix]`` ... ``end`` holds a table: the first line
  is the column header, then one row per line (blank separated, strings with
  blanks quoted, ``nan`` = missing).  ``records`` tables are lists of dicts
  and ``matrix`` tables lists of rows inside the key tree (e.g. route summaries).

:func:`write_result` turns the contract dictionaries of the launcher into a
text file and :func:`read_result` returns the same dictionary again (numbers
in data tables keep 7 significant digits), so tools and FLAAP Studio read the
text files directly.  ``python -m flaap.io.results_txt FILE.txt`` prints a
file as JSON.
"""
from __future__ import annotations

import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

__all__ = ["TxtDoc", "parse_txt", "read_txt", "read_result", "write_result", "result_text", "KINDS",
           "FORMAT"]

FORMAT = "flaap-txt 1"
KINDS = ("results", "modes", "aero_mesh", "cp", "flutter", "timehist", "curves", "tsd_steady")

_INT = re.compile(r"^[+-]?\d+$")
_FLOAT = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_SPECIAL = {"true", "false", "none", "nan", "inf", "+inf", "-inf", "[]", "{}", "end", "table"}


# ----------------------------------------------------------------- quoting
def _quote(s: str) -> str:
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\r", "") + '"'


def _needs_quote(s: str, cell: bool) -> bool:
    if s == "" or s != s.strip() or _INT.match(s) or _FLOAT.match(s) or s.lower() in _SPECIAL:
        return True
    if s[0] in "[{%\"" or any(c in s for c in ',#"\n\r'):
        return True
    return cell and any(c.isspace() for c in s)


def _seg(s: str) -> str:
    """One segment of a dotted key / a table name / a column name."""
    s = str(s)
    if s == "" or any(c.isspace() or c in '.,#"%' for c in s) or s.lower() in _SPECIAL:
        return _quote(s)
    return s


def _fmt(v, cell: bool = False, sig: int | None = None) -> str:
    if v is None:
        return "nan" if cell else "none"
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        v = float(v)
        if not math.isfinite(v):
            return "nan" if cell else "none"
        return repr(v) if sig is None else f"{v:.{sig}g}"
    s = str(v)
    return _quote(s) if _needs_quote(s, cell) else s


def _fmt_value(v) -> str:
    """Right-hand side of a ``key value`` line."""
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            return "[]"
        items = ", ".join(_fmt(x) for x in v)
        return f"[{items}]" if len(v) == 1 else items
    if isinstance(v, dict) and not v:
        return "{}"
    return _fmt(v)


# ----------------------------------------------------------------- tokenizer
def _tokens(line: str, sep=None):
    """Split ``line`` at blanks (``sep=None``) or at ``sep`` outside double quotes.

    Returns (token, was_quoted) pairs; stops at an inline ``#`` comment."""
    out, buf, quoted, inq, i, n = [], [], False, False, 0, len(line)
    has = False
    while i < n:
        c = line[i]
        if inq:
            if c == "\\" and i + 1 < n:
                nx = line[i + 1]
                buf.append({"n": "\n", "t": "\t"}.get(nx, nx))
                i += 2
                continue
            if c == '"':
                inq = False
            else:
                buf.append(c)
        elif c == '"':
            inq, quoted, has = True, True, True
        elif c == "#" and (i == 0 or line[i - 1].isspace() or (sep and line[i - 1] == sep)):
            break
        elif (sep is None and c.isspace()) or (sep is not None and c == sep):
            if sep is not None or has:
                tok = "".join(buf)
                out.append((tok if quoted else tok.strip(), quoted))
            buf, quoted, has = [], False, False
        else:
            buf.append(c)
            has = has or not c.isspace()
        i += 1
    if has or (sep is not None and (out or buf)):
        tok = "".join(buf)
        out.append((tok if quoted else tok.strip(), quoted))
    return out


def _scalar(tok: str, quoted: bool, cell: bool = False):
    if quoted:
        return tok
    if _INT.match(tok):
        return int(tok)
    if _FLOAT.match(tok):
        return float(tok)
    lo = tok.lower()
    if lo in ("true", "false"):
        return lo == "true"
    if lo in ("none", "nan", "inf", "+inf", "-inf", "null"):
        return None
    if lo == "{}" and not cell:
        return {}
    return tok


def _value(rest: str):
    s = rest.strip()
    if s == "[]":
        return []
    bracket = s.startswith("[") and s.endswith("]")
    body = s[1:-1] if bracket else s
    parts = _tokens(body, sep=",")
    if bracket or len(parts) > 1:
        return [_scalar(t, q) for t, q in parts if t != "" or q]
    if not parts:
        return ""
    tok, q = parts[0]
    if q:
        # a quoted string, possibly followed by nothing else
        return tok
    return _scalar(tok, False)


def _split_key(line: str):
    """``key rest`` with quoted key segments; returns (segments, rest)."""
    segs, buf, inq, i, n = [], [], False, 0, len(line)
    while i < n:
        c = line[i]
        if inq:
            if c == "\\" and i + 1 < n:
                buf.append({"n": "\n"}.get(line[i + 1], line[i + 1]))
                i += 2
                continue
            if c == '"':
                inq = False
            else:
                buf.append(c)
        elif c == '"':
            inq = True
        elif c == ".":
            segs.append("".join(buf))
            buf = []
        elif c.isspace():
            break
        else:
            buf.append(c)
        i += 1
    segs.append("".join(buf))
    return segs, line[i:].strip()


# ----------------------------------------------------------------- reader
def parse_txt(text: str) -> dict:
    """Parse a result text file -> ``{"values": nested dict, "tables": {name: table}}``.

    A table is ``{"columns": [...], "rows": [[...], ...], "type": "", "comment": str}``;
    ``records`` / ``matrix`` tables are also placed into ``values`` at their dotted path."""
    values: dict = {}
    tables: dict = {}
    cur = None
    last_comment = ""
    for raw in text.lstrip("\ufeff").splitlines():
        s = raw.strip()
        if not s:
            continue
        if s.startswith("%"):
            if cur is None:
                last_comment = s.lstrip("% ").strip()
            continue
        if cur is not None:
            if s == "end":
                cur = None
                continue
            toks = _tokens(s)
            if not toks:
                continue
            if cur["columns"] is None:
                cur["columns"] = [t for t, _ in toks]
            else:
                cur["rows"].append([_scalar(t, q, True) for t, q in toks])
            continue
        if s.startswith("table ") or s.startswith("table\t"):
            segs, rest = _split_key(s[6:].strip())
            name = ".".join(segs) if len(segs) > 1 else segs[0]
            cur = {"columns": None, "rows": [], "type": rest.strip().lower(), "comment": last_comment,
                   "path": segs}
            tables[name] = cur
            last_comment = ""
            continue
        segs, rest = _split_key(s)
        if not rest.split("#", 1)[0].strip() and not rest.startswith('"'):
            val = ""
        else:
            val = _value(rest)
        _put(values, segs, val)
    for t in tables.values():
        t["columns"] = t["columns"] or []
        if t["type"] == "records":
            _put(values, t["path"], [dict(zip(t["columns"], r)) for r in t["rows"]])
        elif t["type"] == "matrix":
            _put(values, t["path"], [list(r) for r in t["rows"]])
    return {"values": values, "tables": tables}


def _put(d: dict, segs, val):
    for k in segs[:-1]:
        nxt = d.get(k)
        if not isinstance(nxt, dict):
            nxt = d[k] = {}
        d = nxt
    d[segs[-1]] = val


def read_txt(path) -> dict:
    return parse_txt(Path(path).read_text(encoding="utf-8"))


def _col(tab, name, default=None):
    try:
        i = tab["columns"].index(name)
    except ValueError:
        return default
    return [r[i] if i < len(r) else None for r in tab["rows"]]


def _cols(tab, start=0):
    """columns[start:] as a dict name -> list."""
    return {c: [r[i] if i < len(r) else None for r in tab["rows"]] for i, c in enumerate(tab["columns"]) if i >= start}


def _numbered(tables, prefix):
    """tables named ``prefix.1``, ``prefix.2`` ... in numeric order."""
    out = []
    for name, t in tables.items():
        if name.startswith(prefix + "."):
            tail = name[len(prefix) + 1:]
            if tail.isdigit():
                out.append((int(tail), t))
    return [t for _, t in sorted(out, key=lambda x: x[0])]


# ----------------------------------------------------------------- writer
class TxtDoc:
    """Build a result text file line by line."""

    def __init__(self, kind: str, title: str | None = None, file: str | None = None):
        self.lines: list[str] = []
        bar = "% " + "=" * 74
        self.lines += [bar, "% FLAAP result file" + (f"  {file}" if file else ""),
                       f"% {title}" if title else "% ", f"% written {time.strftime('%Y-%m-%d %H:%M:%S')}  "
                       f"(format {FORMAT}: key value lines + 'table ... end' blocks)", bar]
        self._kv: list[tuple[str, str]] = []
        self.kv("kind", kind)
        if title:
            self.kv("title", title)

    # key value ---------------------------------------------------------
    def kv(self, key: str, value):
        self._kv.append((_seg(key), _fmt_value(value)))

    def _kvpath(self, segs, value):
        self._kv.append((".".join(_seg(s) for s in segs), _fmt_value(value)))

    def _emit_kv(self):
        if not self._kv:
            return
        w = max(14, max(len(k) for k, _ in self._kv) + 2)
        self.lines += [f"{k:<{w}s}{v}" for k, v in self._kv]
        self._kv = []

    def comment(self, text: str = ""):
        self._emit_kv()
        self.lines += [f"% {t}" if t else "%" for t in str(text).splitlines() or [""]]

    def section(self, text: str):
        self._emit_kv()
        self.lines += ["", "% " + "-" * 74] + [f"% {t}" for t in text.splitlines()] + ["% " + "-" * 74]

    def blank(self):
        self._emit_kv()
        self.lines.append("")

    # generic tree --------------------------------------------------------
    def tree(self, segs, obj, sig_tables=10):
        """Write a nested dict / list structure under the dotted path ``segs``."""
        if isinstance(segs, str):
            segs = [segs]
        if isinstance(obj, np.ndarray):
            obj = obj.tolist()
        if isinstance(obj, dict):
            if not obj:
                self._kvpath(segs, {})
            for k, v in obj.items():
                self.tree(list(segs) + [str(k)], v, sig_tables)
            return
        if isinstance(obj, (list, tuple)) and obj:
            if all(isinstance(x, dict) for x in obj) and all(
                    not isinstance(v, (dict, list, tuple, np.ndarray)) for x in obj for v in x.values()):
                cols: list[str] = []
                for x in obj:
                    cols += [k for k in x if k not in cols]
                self._table_path(segs, "records", {c: [x.get(c) for x in obj] for c in cols}, sig=sig_tables)
                return
            if all(isinstance(x, (list, tuple, np.ndarray)) for x in obj):
                rows = [list(np.asarray(x).tolist()) if isinstance(x, np.ndarray) else list(x) for x in obj]
                if all(not isinstance(v, (list, dict, tuple)) for r in rows for v in r):
                    ncol = max(len(r) for r in rows)
                    cols = {f"c{j + 1}": [r[j] if j < len(r) else None for r in rows] for j in range(ncol)}
                    self._table_path(segs, "matrix", cols, sig=sig_tables)
                    return
            if any(isinstance(x, (dict, list, tuple)) for x in obj):
                self._kvpath(segs, json.dumps(obj, default=str))
                return
        self._kvpath(segs, obj)

    # tables ------------------------------------------------------------
    def table(self, name, columns: dict, note: str | None = None, sig: int = 7, kind: str = ""):
        segs = name if isinstance(name, (list, tuple)) else [name]
        self._table_path(segs, kind, columns, note=note, sig=sig)

    def _table_path(self, segs, kind, columns: dict, note=None, sig=7):
        self._emit_kv()
        self.lines.append("")
        if note:
            self.lines += [f"% {t}" for t in note.splitlines()]
        name = ".".join(_seg(s) for s in segs)
        self.lines.append(f"table {name}" + (f" {kind}" if kind else ""))
        names = [_seg(c) if any(ch.isspace() or ch in '"#%' for ch in str(c)) or str(c) == "" or
                 str(c).lower() in _SPECIAL else str(c) for c in columns]
        cells = []
        for v in columns.values():
            a = v if isinstance(v, np.ndarray) else list(v)
            if isinstance(a, np.ndarray) and a.dtype.kind == "f":
                cells.append(["nan" if not math.isfinite(x) else f"{x:.{sig}g}" for x in a.tolist()])
            elif isinstance(a, np.ndarray) and a.dtype.kind in "iu":
                cells.append([str(x) for x in a.tolist()])
            else:
                cells.append([_fmt(x, True, sig) for x in (a.tolist() if isinstance(a, np.ndarray) else a)])
        nrow = max((len(c) for c in cells), default=0)
        widths = [max([len(h)] + [len(x) for x in c]) for h, c in zip(names, cells)]
        self.lines.append("  " + "  ".join(f"{h:>{w}s}" for h, w in zip(names, widths)))
        for i in range(nrow):
            self.lines.append("  " + "  ".join(f"{(c[i] if i < len(c) else 'nan'):>{w}s}"
                                                for c, w in zip(cells, widths)))
        self.lines.append("end")

    def text(self) -> str:
        self._emit_kv()
        return "\n".join(self.lines) + "\n"

    def write(self, path):
        Path(path).write_text(self.text(), encoding="utf-8")


# ----------------------------------------------------------------- per-kind writers
def _np(x):
    """Nested lists with None -> float array with NaN."""
    def conv(v):
        if isinstance(v, list):
            return [conv(t) for t in v]
        return np.nan if v is None else v
    return np.asarray(conv(x), dtype=float)


def _w_results(doc: TxtDoc, d: dict):
    for k in ("model", "mode", "status", "error", "elapsed_s", "created", "config", "output_dir"):
        if k in d:
            doc.kv(k, d[k])
    extra = [k for k in d if k not in ("title", "model", "mode", "status", "error", "elapsed_s", "created",
                                       "config", "output_dir", "summary", "artifacts")]
    for k in extra:
        doc.tree([k], d[k])
    doc.section("summary: key results of the run (nested keys as summary.<name>.<field>)")
    doc.tree(["summary"], d.get("summary", {}))
    arts = d.get("artifacts") or []
    doc.section("artifacts: files written by the run")
    if arts:
        doc.tree(["artifacts"], arts)
    else:
        doc.kv("artifacts", [])


def _r_results(p) -> dict:
    v = dict(p["values"])
    v.pop("kind", None)
    v.setdefault("summary", {})
    v.setdefault("artifacts", [])
    if v.get("artifacts") == "":
        v["artifacts"] = []
    return v


def _w_modes(doc: TxtDoc, d: dict):
    labels = list(d.get("labels", []))
    f = _np(d.get("freq_hz", []))
    gm = _np(d.get("gen_mass", [np.nan] * len(f)))
    nodes = _np(d.get("nodes", np.zeros((0, 3))))
    shapes = _np(d.get("shapes", np.zeros((0, 0, 3))))
    doc.kv("name", d.get("name"))
    doc.kv("structure", d.get("structure"))
    doc.kv("n_modes", len(labels))
    doc.kv("n_nodes", int(nodes.shape[0]))
    doc.kv("node_numbering", "1-based (connectivity tables use the same numbers)")
    doc.kv("shape_scaling", "each mode normalised to max |u| = 1")
    doc.table("modes", {"mode": np.arange(1, len(labels) + 1), "label": labels, "freq_hz": f,
                        "omega_rad_s": 2 * np.pi * f, "gen_mass": gm},
              note="natural frequencies and generalized masses", sig=8)
    if nodes.size:
        doc.table("nodes", {"node": np.arange(1, nodes.shape[0] + 1), "x": nodes[:, 0], "y": nodes[:, 1],
                            "z": nodes[:, 2]}, note="structural node coordinates [m]")
    for key, cols in (("edges", ("n1", "n2")), ("quads", ("n1", "n2", "n3", "n4")), ("tris", ("n1", "n2", "n3"))):
        if key in d and len(d[key]):
            c = np.asarray(d[key], dtype=int) + 1
            doc.table(key, {n: c[:, j] for j, n in enumerate(cols)}, note=f"{key} for plotting (node numbers)")
    for i in range(shapes.shape[0]):
        lab = labels[i] if i < len(labels) else f"mode {i + 1}"
        doc.table(["shape", str(i + 1)], {"node": np.arange(1, shapes.shape[1] + 1), "ux": shapes[i, :, 0],
                                          "uy": shapes[i, :, 1], "uz": shapes[i, :, 2]},
                  note=f"mode {i + 1}: {lab}, f = {f[i]:.6g} Hz - nodal translations", sig=5)
    if "aero_h" in d:
        h = _np(d["aero_h"])
        for i in range(h.shape[0]):
            doc.table(["aero_h", str(i + 1)], {"box": np.arange(1, h.shape[1] + 1), **{
                f"h{j + 1}": h[i, :, j] for j in range(h.shape[2])}},
                note=f"mode {i + 1}: normal displacement at the 4 corners of each aero box", sig=5)


def _r_modes(p) -> dict:
    v, t = p["values"], p["tables"]
    m = t.get("modes", {"columns": [], "rows": []})
    d = {"name": v.get("name"), "labels": [str(x) for x in _col(m, "label", [])],
         "freq_hz": _col(m, "freq_hz", []), "gen_mass": _col(m, "gen_mass", [])}
    if "nodes" in t:
        nt = t["nodes"]
        d["nodes"] = [list(r[1:4]) for r in nt["rows"]]
    d["shapes"] = [[list(r[1:4]) for r in s["rows"]] for s in _numbered(t, "shape")]
    d["structure"] = v.get("structure")
    for key in ("edges", "quads", "tris"):
        if key in t:
            d[key] = [[int(x) - 1 for x in r] for r in t[key]["rows"]]
    ah = _numbered(t, "aero_h")
    if ah:
        d["aero_h"] = [[list(r[1:]) for r in s["rows"]] for s in ah]
    return d


def _w_aero_mesh(doc: TxtDoc, d: dict):
    c = _np(d["corners"])
    col = _np(d.get("colloc", np.zeros((c.shape[0], 3))))
    for k in ("surfaces", "symmetry", "bref", "mach"):
        doc.kv(k, d.get(k))
    ctrl = d.get("controls", {}) or {}
    doc.kv("controls", list(ctrl))
    doc.kv("n_boxes", int(c.shape[0]))
    cols = {"box": np.arange(1, c.shape[0] + 1), "surface": d.get("surface_id", [0] * c.shape[0]),
            "area": _np(d.get("area", [np.nan] * c.shape[0])), "xc": col[:, 0], "yc": col[:, 1], "zc": col[:, 2]}
    for j in range(c.shape[1]):
        for a, ax in enumerate("xyz"):
            cols[f"{ax}{j + 1}"] = c[:, j, a]
    for n, mask in ctrl.items():
        cols[f"ctrl:{n}"] = [bool(x) for x in mask]
    doc.table("boxes", cols, note="aerodynamic boxes: surface index (0-based), area [m^2], collocation point "
              "(xc, yc, zc) and corners 1-4 [m];\nctrl:<name> = box belongs to that control surface")


def _r_aero_mesh(p) -> dict:
    v, t = p["values"], p["tables"]
    b = t["boxes"]
    cols = b["columns"]
    corner_cols = [c for c in cols if re.fullmatch(r"[xyz]\d+", c)]
    ncorner = max(int(c[1:]) for c in corner_cols) if corner_cols else 0
    idx = {c: i for i, c in enumerate(cols)}
    corners = [[[r[idx[f"{ax}{j + 1}"]] for ax in "xyz"] for j in range(ncorner)] for r in b["rows"]]
    controls = {c[5:]: [bool(r[idx[c]]) for r in b["rows"]] for c in cols if c.startswith("ctrl:")}
    surfaces = v.get("surfaces")
    if isinstance(surfaces, str):
        surfaces = [surfaces]
    return {"corners": corners, "surface_id": [int(x) for x in _col(b, "surface", [])],
            "surfaces": surfaces or [], "symmetry": v.get("symmetry"),
            "colloc": [[r[idx["xc"]], r[idx["yc"]], r[idx["zc"]]] for r in b["rows"]],
            "area": _col(b, "area", []), "controls": controls, "bref": v.get("bref"), "mach": v.get("mach")}


def _as_list(x):
    if x is None or x == "":
        return []
    return x if isinstance(x, list) else [x]


def _w_cp(doc: TxtDoc, d: dict):
    k = _np(d["k"]).ravel()
    modes = list(d["modes"])
    re_, im_ = _np(d["re"]), _np(d["im"])
    for key in ("mach", "bref"):
        doc.kv(key, d.get(key))
    doc.kv("modes", modes)
    doc.kv("k", k.tolist())
    doc.kv("n_boxes", int(re_.shape[2]) if re_.ndim == 3 else 0)
    for ik in range(re_.shape[0]):
        cols = {"box": np.arange(1, re_.shape[2] + 1)}
        for im, m in enumerate(modes):
            cols[f"re:{m}"] = re_[ik, im]
            cols[f"im:{m}"] = im_[ik, im]
        doc.table(["cp", str(ik + 1)], cols, note=f"reduced frequency k = {k[ik]:.6g}: lifting pressure dCp "
                  "per unit generalized coordinate (real / imaginary part)", sig=5)


def _r_cp(p) -> dict:
    v, t = p["values"], p["tables"]
    modes = [str(m) for m in _as_list(v.get("modes"))]
    re_, im_ = [], []
    for tab in _numbered(t, "cp"):
        re_.append([_col(tab, f"re:{m}", []) for m in modes])
        im_.append([_col(tab, f"im:{m}", []) for m in modes])
    return {"k": _as_list(v.get("k")), "modes": modes, "re": re_, "im": im_, "mach": v.get("mach"),
            "bref": v.get("bref")}


def _w_flutter(doc: TxtDoc, d: dict):
    names = list(d.get("names", []))
    V, F, G, K = (_np(d[x]) for x in ("V", "freq", "g", "k"))
    if G.ndim == 1:
        G, F, K = G[:, None], F[:, None], K[:, None]
    nb = G.shape[1]
    per = V.ndim == 2 and nb > 1 and not np.allclose(np.nan_to_num(V - V[:, :1]), 0)
    for key in ("method", "method_title", "density", "bref", "g_struct"):
        doc.kv(key, d.get(key))
    doc.kv("names", names)
    doc.kv("n_velocities", int(G.shape[0]))
    doc.kv("v_per_branch", bool(per))
    pts = d.get("points") or []
    doc.section("identified flutter (g = 0 crossing) and divergence points;\n"
                "branch = 0-based branch index (its V-g table is branch.<branch+1>)")
    cols = ["kind", "branch", "velocity", "freq_hz", "k", "density"]
    for x in pts:
        cols += [c for c in x if c not in cols]
    doc.table("points", {c: [x.get(c) for x in pts] for c in cols}, kind="records", sig=8)
    doc.section("V-g / V-f per branch: velocity V [m/s], damping g, frequency [Hz], reduced frequency k")
    for b in range(nb):
        Vb = V[:, b] if V.ndim == 2 else V
        doc.table(["branch", str(b + 1)], {"V": Vb, "g": G[:, b], "freq_hz": F[:, b], "k": K[:, b]},
                  note=f"branch {b + 1}: {names[b] if b < len(names) else ''}")


def _r_flutter(p) -> dict:
    v, t = p["values"], p["tables"]
    br = _numbered(t, "branch")
    cols = [{c: _col(tb, c, []) for c in ("V", "g", "freq_hz", "k")} for tb in br]
    n = len(cols[0]["V"]) if cols else 0

    def mat(key):
        return [[c[key][i] for c in cols] for i in range(n)]
    V = mat("V") if (v.get("v_per_branch") or not cols) else [[cols[0]["V"][i]] * len(cols) for i in range(n)]
    return {"method": v.get("method"), "method_title": v.get("method_title"), "V": V, "freq": mat("freq_hz"),
            "g": mat("g"), "k": mat("k"), "names": [str(x) for x in _as_list(v.get("names"))],
            "density": v.get("density"), "bref": v.get("bref"), "g_struct": v.get("g_struct"),
            "points": v.get("points") or []}


def _w_timehist(doc: TxtDoc, d: dict):
    doc.tree(["meta"], d.get("meta") or {})
    doc.tree(["units"], d.get("units") or {})
    series = d.get("series", {})
    doc.kv("series", list(series))
    doc.table("timehist", {"t": _np(d["t"]), **{k: _np(s) for k, s in series.items()}},
              note="time [s] and one column per output signal", sig=7)


def _r_timehist(p) -> dict:
    v, t = p["values"], p["tables"]
    tab = t["timehist"]
    cols = _cols(tab, 1)
    return {"t": _col(tab, "t", []), "series": cols, "meta": v.get("meta") or {}, "units": v.get("units") or {}}


def _w_curves(doc: TxtDoc, d: dict):
    for k in ("xlabel", "ylabel", "logy"):
        doc.kv(k, d.get(k))
    doc.tree(["meta"], d.get("meta") or {})
    curves = d.get("curves", {})
    doc.kv("curves", list(curves))
    doc.table("curves", {"x": _np(d["x"]), **{k: _np(s) for k, s in curves.items()}},
              note=f"x = {d.get('xlabel', '')}; one column per curve ({d.get('ylabel', '')})", sig=7)


def _r_curves(p) -> dict:
    v, t = p["values"], p["tables"]
    tab = t["curves"]
    return {"x": _col(tab, "x", []), "curves": _cols(tab, 1), "xlabel": v.get("xlabel", ""),
            "ylabel": v.get("ylabel", ""), "logy": bool(v.get("logy")), "meta": v.get("meta") or {}}


def _w_tsd(doc: TxtDoc, d: dict):
    x, y, cp = _np(d["x"]), _np(d["y"]).ravel(), _np(d["cp"])
    for k in ("mach", "cp_star", "converged"):
        doc.kv(k, d.get(k))
    doc.kv("residuals", _np(d.get("residuals", [])).tolist())
    doc.kv("ny", int(cp.shape[0]))
    doc.kv("nx", int(cp.shape[1]))
    lifting = d.get("cp_lower") is not None
    if lifting:
        for k in ("alpha_deg", "CL"):
            doc.kv(k, d.get(k))
    J, I = np.meshgrid(np.arange(cp.shape[0]), np.arange(cp.shape[1]), indexing="ij")
    cols = {"j": J.ravel() + 1, "i": I.ravel() + 1, "y": np.repeat(y, cp.shape[1]), "x": x.ravel(),
            "cp": cp.ravel()}
    note = "upper-surface steady Cp on the wing grid (j = span station, i = chordwise point)"
    if lifting:
        cols["cp_lower"] = _np(d["cp_lower"]).ravel()
        note = "steady Cp on the wing grid, upper (cp) and lower (cp_lower) surface (j = span station, " \
               "i = chordwise point)"
    doc.table("grid", cols, note=note, sig=6)
    if lifting:
        doc.table("strips", {"j": np.arange(cp.shape[0]) + 1, "y": y, "cl": _np(d["cl"]),
                             "xs_upper": _np(d["shock_upper"]), "dcp_upper": _np(d["jump_upper"]),
                             "xs_lower": _np(d["shock_lower"]), "dcp_lower": _np(d["jump_lower"])},
                  note="section lift and shock position x/c (NaN = no shock) with the Cp rise across it", sig=6)


def _r_tsd(p) -> dict:
    v, t = p["values"], p["tables"]
    ny, nx = int(v["ny"]), int(v["nx"])
    g = t["grid"]
    xs, cps, ys = _col(g, "x"), _col(g, "cp"), _col(g, "y")
    d = {"x": [xs[j * nx:(j + 1) * nx] for j in range(ny)], "y": [ys[j * nx] for j in range(ny)],
         "cp": [cps[j * nx:(j + 1) * nx] for j in range(ny)], "cp_star": v.get("cp_star"),
         "mach": v.get("mach"), "converged": v.get("converged"), "residuals": _as_list(v.get("residuals"))}
    if "cp_lower" in g["columns"]:
        cl = _col(g, "cp_lower")
        d["cp_lower"] = [cl[j * nx:(j + 1) * nx] for j in range(ny)]
        d["alpha_deg"], d["CL"] = v.get("alpha_deg"), v.get("CL")
        s = t["strips"]
        for k_out, k_in in (("cl", "cl"), ("shock_upper", "xs_upper"), ("jump_upper", "dcp_upper"),
                            ("shock_lower", "xs_lower"), ("jump_lower", "dcp_lower")):
            d[k_out] = _col(s, k_in)
    return d


_WRITERS = {"results": _w_results, "modes": _w_modes, "aero_mesh": _w_aero_mesh, "cp": _w_cp,
            "flutter": _w_flutter, "timehist": _w_timehist, "curves": _w_curves, "tsd_steady": _w_tsd}
_READERS = {"results": _r_results, "modes": _r_modes, "aero_mesh": _r_aero_mesh, "cp": _r_cp,
            "flutter": _r_flutter, "timehist": _r_timehist, "curves": _r_curves, "tsd_steady": _r_tsd}


def result_text(kind: str, data: dict, title: str | None = None, file: str | None = None) -> str:
    doc = TxtDoc(kind, title if title is not None else data.get("title"), file)
    if kind in _WRITERS:
        _WRITERS[kind](doc, data)
    else:
        for k, v in data.items():
            if k != "title":
                doc.tree([k], v)
    return doc.text()


def write_result(path, kind: str, data: dict, title: str | None = None) -> Path:
    """Write a contract dictionary of ``kind`` as a result text file."""
    p = Path(path)
    p.write_text(result_text(kind, data, title, p.name), encoding="utf-8")
    return p


def read_result(path, kind: str | None = None) -> dict:
    """Read a result text file back into its contract dictionary."""
    p = parse_txt(Path(path).read_text(encoding="utf-8"))
    kind = kind or p["values"].get("kind")
    if kind in _READERS:
        d = _READERS[kind](p)
        if kind != "results" and "title" in p["values"]:
            d.setdefault("title", p["values"]["title"])
        return d
    v = dict(p["values"])
    v.pop("kind", None)
    return v


def _main(argv):
    if not argv:
        print(__doc__)
        return 2
    for a in argv:
        print(json.dumps(read_result(a), indent=1, default=str, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
