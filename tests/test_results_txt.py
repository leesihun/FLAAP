"""Round trip of the human-readable result text format (flaap.io.results_txt)."""
import math

import numpy as np
import pytest

from flaap.io.results_txt import TxtDoc, parse_txt, read_result, result_text, write_result


def _close(a, b, rel=1e-6):
    if isinstance(a, dict):
        assert set(a) == set(b), (set(a), set(b))
        for k in a:
            _close(a[k], b[k], rel)
    elif isinstance(a, (list, tuple, np.ndarray)):
        a, b = list(a), list(b)
        assert len(a) == len(b)
        for x, y in zip(a, b):
            _close(x, y, rel)
    elif isinstance(a, float) and isinstance(b, (int, float)) and b is not None:
        if math.isnan(a):
            assert b is None or math.isnan(b)
        else:
            assert b == pytest.approx(a, rel=rel, abs=1e-12)
    elif a is None or (isinstance(a, float) and math.isnan(a)):
        assert b is None or (isinstance(b, float) and math.isnan(b))
    else:
        assert a == b, (a, b)


def test_grammar_scalars_lists_strings():
    doc = TxtDoc("results", "case, with comma")
    doc.kv("a", 1.5)
    doc.kv("b", [1, 2, 3])
    doc.kv("one", [7.0])
    doc.kv("empty", [])
    doc.kv("flag", True)
    doc.kv("missing", None)
    doc.kv("numstr", "42")
    doc.kv("boolstr", "true")
    doc.kv("hash", "a # b")
    doc.kv("names", ["mode 1", "mode 2"])
    doc.tree(["summary"], {"flutter_pk": {"velocity": 153.28712345678, "mode": 2}})
    v = parse_txt(doc.text())["values"]
    assert v["title"] == "case, with comma"
    assert v["a"] == 1.5 and v["b"] == [1, 2, 3] and v["one"] == [7.0] and v["empty"] == []
    assert v["flag"] is True and v["missing"] is None
    assert v["numstr"] == "42" and v["boolstr"] == "true" and v["hash"] == "a # b"
    assert v["names"] == ["mode 1", "mode 2"]
    assert v["summary"]["flutter_pk"]["velocity"] == 153.28712345678   # exact repr round trip
    assert v["summary"]["flutter_pk"]["mode"] == 2


def test_results_records_and_nan():
    res = {"title": "t", "model": "flutter", "mode": "pk", "status": "ok", "error": None, "elapsed_s": 1.25,
           "summary": {"frequencies_hz": [1.0, 2.5], "rows": [{"V": 10.0, "kind": "lco", "A": None},
                                                            {"V": 12.0, "kind": "decay", "A": 0.5}]},
           "artifacts": [{"file": "modes.txt", "kind": "modes", "title": "Mode shapes and frequencies"}]}
    back = read_result_text(res, "results")
    assert back["summary"]["rows"][0]["A"] is None and back["summary"]["rows"][1]["kind"] == "decay"
    assert back["artifacts"][0]["title"] == "Mode shapes and frequencies"
    assert back["status"] == "ok" and back["error"] is None


def read_result_text(data, kind, tmp=None):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / f"{kind}.txt"
        write_result(p, kind, data)
        return read_result(p)


def test_modes_roundtrip():
    rng = np.random.default_rng(0)
    d = {"name": "plate", "structure": "plate", "labels": ["mode 1", "mode 2"], "freq_hz": [1.5, 7.25],
         "gen_mass": [1.0, 1.0], "nodes": rng.random((5, 3)).tolist(), "shapes": rng.random((2, 5, 3)).tolist(),
         "quads": [[0, 1, 2, 3]], "aero_h": rng.random((2, 3, 4)).tolist()}
    back = read_result_text(d, "modes")
    _close(d, back, rel=1e-4)


def test_flutter_roundtrip_shared_and_per_branch_V():
    V = np.linspace(10, 100, 7)
    g = np.column_stack([np.linspace(-0.1, 0.05, 7), np.full(7, -0.02)])
    d = {"method": "pk", "method_title": "PK method", "V": np.column_stack([V, V]).tolist(), "g": g.tolist(),
         "freq": np.ones((7, 2)).tolist(), "k": np.full((7, 2), 0.3).tolist(), "names": ["mode 1", "mode 2"],
         "density": 1.225, "bref": 0.9, "g_struct": 0.0,
         "points": [{"kind": "flutter", "branch": 0, "velocity": 80.0, "freq_hz": 1.0, "k": 0.3, "density": 1.225}]}
    back = read_result_text(d, "flutter")
    _close(d, back)
    d2 = dict(d, V=np.column_stack([V, V * 1.1]).tolist(), method="k")
    _close(d2, read_result_text(d2, "flutter"))


def test_timehist_curves_cp_tsd_aero_mesh():
    t = np.linspace(0, 1, 11)
    th = {"t": t.tolist(), "series": {"root.Mx": np.sin(t).tolist(), "eta:mode 1": np.cos(t).tolist()},
          "meta": {"H": 30.0, "label": "critical gust"}, "units": {"root.Mx": "N m"}}
    _close(th, read_result_text(th, "timehist"))
    cv = {"x": t.tolist(), "curves": {"A [deg]:up": [None] + t[1:].tolist()}, "xlabel": "V", "ylabel": "A",
          "logy": False, "meta": {}}
    _close(cv, read_result_text(cv, "curves"))
    cp = {"k": [0.1, 0.5], "modes": ["mode 1", "gust"], "re": np.ones((2, 2, 4)).tolist(),
          "im": np.zeros((2, 2, 4)).tolist(), "mach": 0.5, "bref": 1.0}
    _close(cp, read_result_text(cp, "cp"))
    ts = {"x": [[0.0, 0.5, 1.0], [0.0, 0.4, 0.8]], "y": [0.0, 1.0], "cp": [[-0.1, -0.5, 0.1], [-0.2, -0.4, 0.0]],
          "cp_star": -0.08, "mach": 0.954, "converged": True, "residuals": [1.0, 0.1, 0.01]}
    _close(ts, read_result_text(ts, "tsd_steady"))
    am = {"corners": np.arange(24, dtype=float).reshape(2, 4, 3).tolist(), "surface_id": [0, 0],
          "surfaces": ["wing"], "symmetry": True, "colloc": [[0.0, 0.1, 0.0], [0.2, 0.3, 0.0]], "area": [0.1, 0.2],
          "controls": {"aileron": [False, True]}, "bref": 1.0, "mach": 0.0}
    _close(am, read_result_text(am, "aero_mesh"))


def test_text_is_human_readable():
    txt = result_text("flutter", {"method": "pk", "V": [[10.0], [20.0]], "g": [[-0.1], [0.1]], "freq": [[1.0], [1.0]],
                                  "k": [[0.2], [0.1]], "names": ["mode 1"], "points": []}, "demo")
    assert "{" not in txt and "table branch.1" in txt and "end" in txt
    assert txt.splitlines()[5].startswith("kind")
