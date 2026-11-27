"""Launcher tests: flat-config grammar, preflight diagnostics, shipped configs, a real run."""
from pathlib import Path

import pytest

from flaap.launcher import parse_text, resolve, run_config
from flaap.launcher.diagnostics import Severity
from flaap.launcher.main import main as cli_main
from flaap.launcher.parser import write_config, parse_config
from flaap.io.results_txt import read_result

ROOT = Path(__file__).resolve().parents[1]


def _codes(rep, sev=None):
    return {d.code for d in rep.diagnostics if sev is None or d.severity is sev}


def _check(text):
    pc = parse_text(text)
    vals, rep = resolve(pc)
    return pc, vals, rep


def test_grammar_values():
    pc = parse_text("% comment\nmodel flutter   # inline\nmode pk\nvelocities 60, 260, 101\n"
                    "plots false\ntitle My Wing Case\nmach 0.3\n")
    v = pc.values
    assert v["model"] == "flutter" and v["mode"] == "pk"
    assert v["velocities"] == [60, 260, 101]
    assert v["plots"] is False
    assert v["title"] == "My Wing Case"          # text keys keep case
    assert v["mach"] == pytest.approx(0.3)


def test_duplicate_and_bom():
    pc = parse_text("﻿model modal\nmodel flutter\n")
    codes = _codes(pc.diagnostics)
    assert "CFG-FILE-005" in codes and "CFG-DUP-001" in codes


def test_unknown_key_suggestion():
    _, _, rep = _check("model flutter\nmode pk\nvelocites 60, 260, 101\n")
    d = [d for d in rep.diagnostics if d.code == "CFG-KEY-001"]
    assert d and d[0].severity is Severity.WARNING
    assert rep.has_errors(strict=True) and not rep.has_errors(strict=False)
    assert "velocities" in (d[0].hint or "") + d[0].message


def test_route_and_mode_errors():
    _, _, rep = _check("model flutterr\n")
    assert "CFG-ROUTE-001" in _codes(rep, Severity.ERROR)
    _, _, rep = _check("model flutter\nmode q\n")
    assert "CFG-ROUTE-002" in _codes(rep, Severity.ERROR)


def test_required_and_defaults():
    _, vals, rep = _check("model gust\nmode discrete\n")
    assert rep.has_errors()                      # velocity required
    _, vals, rep = _check("model flutter\n")
    assert not rep.has_errors()
    assert vals["mode"] == "pk" and vals["structure"] == "goland"
    assert vals["velocities"] == [10.0, 300.0, 146]


def test_value_checks():
    _, _, rep = _check("model flutter\nmode pk\nsymmetry sideways\n")
    assert rep.has_errors()
    _, _, rep = _check("model flutter\nmode pk\nnmodes -2\n")
    assert rep.has_errors()
    _, _, rep = _check("model tsd\nmode steady\nstructure agard445\nmach 1.2\n")
    assert rep.has_errors()


def test_file_structure_needs_file(tmp_path):
    _, _, rep = _check("model flutter\nmode pk\nstructure nastran\n")
    assert rep.has_errors()


def test_round_trip(tmp_path):
    _, vals, rep = _check("model flutter\nmode pk\nvelocities 60, 260, 101\ntitle Round Trip\n")
    p = tmp_path / "config_rt.txt"
    write_config(vals, p)
    pc2 = parse_config(p)
    vals2, rep2 = resolve(pc2)
    assert not rep2.has_errors()
    assert vals2["velocities"] == vals["velocities"] and vals2["title"] == "Round Trip"


def test_shipped_configs_pass_preflight():
    files = sorted((ROOT / "configs").rglob("config_*.txt"))
    assert len(files) >= 15
    for f in files:
        _, rep = resolve(parse_config(f))
        assert not rep.has_errors(), f"{f}: {[d.message for d in rep.diagnostics if d.severity is Severity.ERROR]}"


def test_cli_check_and_list(capsys):
    assert cli_main(["--list-models", "--no-color"]) == 0
    assert "flutter" in capsys.readouterr().out
    cfg = ROOT / "configs/flutter/ex1/config_goland_pk.txt"
    assert cli_main(["--config", str(cfg), "--check", "--no-color"]) == 0


def test_run_modal_and_flutter(tmp_path):
    res = run_config(ROOT / "configs/modal/ex1/config_goland_modal.txt", output_dir=tmp_path / "m", verbose=False)
    assert res["status"] == "ok"
    modes = read_result(tmp_path / "m/modes.txt")
    assert len(modes["freq_hz"]) == 4 and modes["freq_hz"][0] == pytest.approx(7.66, rel=0.01)
    assert len(modes["aero_h"]) == 4

    res = run_config(ROOT / "configs/flutter/ex1/config_goland_pk.txt", output_dir=tmp_path / "f", verbose=False)
    assert res["status"] == "ok"
    crit = res["summary"]["flutter_pk"]
    assert crit["kind"] == "flutter" and 140 < crit["velocity"] < 170
    fl = read_result(tmp_path / "f/flutter_pk.txt")
    assert fl["method"] == "pk" and len(fl["V"]) == 101
    assert (tmp_path / "f/results.txt").exists() and (tmp_path / "f/cp.txt").exists()


def test_assignment_syntax():
    pc = parse_text("model = flutter\nmode: pk\nvelocities=60,260,101\ntitle a = b\n")
    v = pc.values
    assert v["model"] == "flutter" and v["mode"] == "pk"
    assert v["velocities"] == [60, 260, 101]
    assert v["title"] == "a = b"                 # value text after the first word is kept verbatim
    assert "CFG-SYNTAX-003" in _codes(pc.diagnostics)
    _, _, rep = _check("model = flutter\nmode = pk\n")
    assert not rep.has_errors()


def test_freeplay_ase_preflight():
    _, vals, rep = _check("model freeplay\nmode time\n")
    assert rep.has_errors()                      # velocity required for a time history
    _, vals, rep = _check("model freeplay\n")
    assert not rep.has_errors()
    assert vals["mode"] == "describing_function" and vals["structure"] == "section"
    _, _, rep = _check("model freeplay\nnonlinearity cubic\n")
    assert rep.has_errors()
    _, _, rep = _check("model ase\nstructure goland\n")
    assert rep.has_errors()                      # only section / goland_flap plants
    _, vals, rep = _check("model ase\nmode margins\nstructure goland_flap\n")
    assert not rep.has_errors() and vals["actuator_damping"] == 0.7


def test_run_freeplay_and_ase(tmp_path):
    res = run_config(ROOT / "configs/freeplay/ex2/config_section_time.txt", output_dir=tmp_path / "fp", verbose=False)
    assert res["status"] == "ok"
    lco = res["summary"]["lco"]
    assert lco["kind"] == "lco" and lco["freq_hz"] == pytest.approx(2.54, rel=0.02)
    assert lco["amplitude_deg"] == pytest.approx(1.209, rel=0.03)
    th = read_result(tmp_path / "fp/timehist_lco.txt")
    assert len(th["t"]) > 100 and th["series"]

    res = run_config(ROOT / "configs/ase/ex1/config_section_afs.txt", output_dir=tmp_path / "ase", verbose=False)
    assert res["status"] == "ok"
    sm = res["summary"]
    assert sm["V_f_open"] == pytest.approx(9.245, rel=2e-3)
    assert sm["V_f_closed"] == pytest.approx(17.40, rel=5e-3)
    assert (tmp_path / "ase/flutter_open.txt").exists() and (tmp_path / "ase/flutter_closed.txt").exists()
