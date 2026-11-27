"""CalculiX ``.frd`` / ``.dat`` modal import (flaap.io.calculix) and the launcher ``calculix`` structure."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from flaap.io import calculix as C  # noqa: E402
from flaap.launcher.diagnostics import Severity  # noqa: E402
from flaap.structures.modal import ModalModel  # noqa: E402

CDIR = os.path.join(ROOT, "examples", "calculix")


@pytest.fixture(scope="module")
def ref():
    r = np.load(os.path.join(CDIR, "reference_modes.npz"))
    return ModalModel(r["nodes"], r["shapes"], 2 * np.pi * r["freq_hz"], np.ones(r["freq_hz"].size))


def _mac(a, b):
    a, b = a.ravel(), b.ravel()
    return (a @ b) ** 2 / ((a @ a) * (b @ b))


def _same(mm, ref, rtol_f=1e-6, atol_sh=1e-5):
    np.testing.assert_allclose(mm.freq_hz, ref.freq_hz[:mm.nmodes], rtol=rtol_f)
    np.testing.assert_allclose(mm.nodes, ref.nodes, atol=1e-6)
    for m in range(mm.nmodes):
        a = mm.shapes[:, :3, m] / np.sqrt(mm.gen_mass[m])
        b = ref.shapes[:, :3, m] / np.sqrt(ref.gen_mass[m])
        assert np.abs(a - b).max() < atol_sh * np.abs(b).max()
        assert _mac(a, b) > 0.999999
    assert np.all(mm.shapes[:, 3:] == 0)                               # ccx writes translations only


def test_frd_parse_sample():
    d = C.read_frd(os.path.join(CDIR, "plate.frd"))
    assert len(d["nodes"]) == 91 and len(d["elements"]) == 72
    assert {t for t, _ in d["elements"].values()} == {9}
    b = d["blocks"][0]
    assert b["name"] == "DISP" and b["ictype"] == 2 and b["analysis"] == "MODAL" and b["numstp"] == 1
    assert b["components"] == ["D1", "D2", "D3", "ALL"] and b["format"] == 1
    assert b["value"] == pytest.approx(3.37144) and len(b["data"]) == 91 and b["data"][1].size == 3
    vis = C.frd_surface(d, sorted(d["nodes"]))
    assert vis["tris"].shape == (144, 3) and vis["edges"].shape == (0, 2)


def test_frd_long_with_dat(ref):
    mm = C.read_calculix_modal(os.path.join(CDIR, "plate.frd"))          # plate.dat found beside it
    assert mm.nmodes == 6 and mm.labels[0] == "ccx mode 1"
    np.testing.assert_allclose(mm.gen_mass, 1.0)
    _same(mm, ref, rtol_f=1e-7)
    hdr = C.read_calculix_modal(os.path.join(CDIR, "plate.frd"), frequency_file=None, frequencies=None,
                                nmodes=3)
    assert hdr.nmodes == 3
    f = C.read_calculix_dat_frequencies(os.path.join(CDIR, "plate.dat"))
    np.testing.assert_allclose(f, ref.freq_hz, rtol=1e-7)


def test_frd_short_mm_static_and_rigid(ref):
    """Short node format, mm-tonne-s, a static step and 3 rigid modes before the elastic ones."""
    mm = C.read_calculix_modal(os.path.join(CDIR, "plate_mm_short.frd"), length_scale=1e-3, mass_scale=1e3)
    assert mm.nmodes == 6 and mm.labels[0] == "ccx mode 4"
    np.testing.assert_allclose(mm.gen_mass, 1e3 * 1e-6)                # Ms Ls^2
    _same(mm, ref, rtol_f=5e-5)                                         # E12.5 header frequencies
    raw = C.read_calculix_modal(os.path.join(CDIR, "plate_mm_short.frd"), skip_below_hz=-1.0)
    assert raw.nmodes == 9 and np.all(raw.freq_hz[:3] == 0)
    assert np.ptp(raw.nodes[:, 1]) == pytest.approx(1000.0)


def test_frd_dat_mismatch_and_two_steps(tmp_path, ref):
    p = tmp_path / "two.frd"
    C.write_frd(ref, str(p), long=False)
    txt = p.read_text().splitlines()
    body = [ln for ln in txt if ln.strip() != "9999"]
    i0 = next(i for i, ln in enumerate(body) if ln.startswith("    1PSTEP"))
    p.write_text("\n".join(body + body[i0:] + ["9999"]) + "\n")       # same modes again: 2nd step
    d = C.read_frd(str(p))
    assert len([b for b in d["blocks"] if b["name"] == "DISP"]) == 12
    a = C.read_calculix_modal(str(p))                                   # last step by default
    b = C.read_calculix_modal(str(p), step=0)
    assert a.nmodes == b.nmodes == 6
    dat = tmp_path / "other.dat"
    C.write_calculix_dat(ModalModel(ref.nodes, ref.shapes, 2 * ref.omega, ref.gen_mass), str(dat))
    with pytest.warns(UserWarning, match="does not match"):
        c = C.read_calculix_modal(str(p), frequency_file=str(dat))
    np.testing.assert_allclose(c.freq_hz, ref.freq_hz, rtol=5e-5)


def test_frd_solid_surface_and_binary(tmp_path):
    from flaap.structures.solid import box_mesh
    nodes, conn = box_mesh(1.0, 0.2, 0.1, 4, 2, 1, "hex8")
    nm = 2
    sh = np.zeros((len(nodes), 6, nm))
    sh[:, 2, 0] = nodes[:, 0] ** 2
    sh[:, 1, 1] = nodes[:, 0] ** 2
    mm = ModalModel(nodes, sh, np.array([10.0, 20.0]), np.ones(nm))
    p = tmp_path / "beam.frd"
    C.write_frd(mm, str(p), elements=[(1, [int(c) + 1 for c in e]) for e in conn])
    d = C.read_frd(str(p))
    vis = C.frd_surface(d, sorted(d["nodes"]))
    # 4x2x1 hex box: boundary quads = 2*(4*2) + 2*(4*1) + 2*(2*1) = 28 -> 56 triangles
    assert vis["tris"].shape == (56, 3)
    b = tmp_path / "bin.frd"
    b.write_text(f"    1Cbin\n    2C{'':18s}{8:12d}{'':37s}3\n")
    with pytest.raises(ValueError, match="binary"):
        C.read_frd(str(b))


def test_fortran_exponent_without_e():
    assert C._f("1.23456-100") == pytest.approx(1.23456e-100)
    assert C._f("-1.5+003") == pytest.approx(-1500.0)


def test_calculix_launcher_preflight_and_run(tmp_path):
    from flaap.launcher import run_config
    from flaap.launcher.parser import parse_text
    from flaap.launcher.specs import resolve

    def check(*lines):
        txt = "\n".join(("model modal", "mode analysis", "structure calculix") + lines) + "\n"
        _, rep = resolve(parse_text(txt, tmp_path / "config_x.txt"))
        return rep.diagnostics

    def codes(diags, sev=None):
        return {d.code for d in diags if sev is None or d.severity is sev}

    frd = os.path.join(CDIR, "plate.frd").replace("\\", "/")
    ok = check(f"structure_file {frd}")
    assert not codes(ok, Severity.ERROR) and any("no unit system" in d.message for d in ok)
    assert not any("no unit system" in d.message for d in check(f"structure_file {frd}", "length_scale 0.001"))
    assert "CFG-REQ-002" in codes(check(), Severity.ERROR)
    assert codes(check(f"structure_file {frd}", "frequency_file nope.dat"), Severity.ERROR) == {"CFG-PATH-001"}
    (tmp_path / "bin.frd").write_bytes(b"    1Cbin\n" + f"    2C{'':18s}{8:12d}{'':37s}3\n".encode() + b"\x00\x01")
    assert "CFG-PATH-003" in codes(check("structure_file bin.frd"), Severity.ERROR)
    (tmp_path / "geo.frd").write_text(f"    1Cgeo\n    2C{'':18s}{0:12d}{'':37s}1\n -3\n9999\n")
    assert "CFG-PATH-002" in codes(check("structure_file geo.frd"), Severity.ERROR)
    # flutter route needs surfaces + bref
    txt = "model flutter\nmode pk\nstructure calculix\nstructure_file " + frd + "\n"
    _, rep = resolve(parse_text(txt, tmp_path / "config_f.txt"))
    assert {"CFG-REQ-003", "CFG-REQ-004"} <= codes(rep.diagnostics, Severity.ERROR)
    res = run_config(os.path.join(ROOT, "configs", "modal", "ex8", "config_calculix_modal.txt"),
                     output_dir=tmp_path / "out", verbose=False, overrides={"plots": False})
    assert res["status"] == "ok" and int(res["summary"]["nmodes"]) == 6
