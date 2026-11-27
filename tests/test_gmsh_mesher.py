"""gmsh mesher (flaap.structures.gmsh_mesher) and the launcher ``gmsh`` structure."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

pytest.importorskip("gmsh")

from flaap.launcher.diagnostics import Severity  # noqa: E402
from flaap.structures import gmsh_mesher as G  # noqa: E402

GDIR = os.path.join(ROOT, "examples", "gmsh")
AL = dict(E=71.7e9, nu=0.33, rho=2810.0)
PLATE_HZ = [3.37144, 14.3225, 21.3217, 47.3130, 61.7587, 93.4401]     # examples/calculix reference


def test_geo_plate_quads_groups_and_modes():
    gm = G.gmsh_mesh(os.path.join(GDIR, "plate_wing.geo"))
    assert gm.dim == 2 and gm.kinds() == ["quad4"] and gm.nelem == 72 and len(gm.nodes) == 91
    root = G.group_nodes(gm, "Root")
    assert len(root) == 7 and np.allclose(gm.nodes[root, 1], 0.0)
    res = G.gmsh_modal(os.path.join(GDIR, "plate_wing.geo"), 6, material=AL, thickness=0.004, clamp_group="Root")
    np.testing.assert_allclose(res.modal.freq_hz, PLATE_HZ, rtol=2e-5)
    with pytest.raises(ValueError, match="not found"):
        G.group_nodes(gm, "Tip")


def test_mesher_options_tri_quad_order_and_msh_roundtrip(tmp_path):
    geo = tmp_path / "sq.geo"
    geo.write_text('SetFactory("OpenCASCADE");\nRectangle(1) = {0, 0, 0, 0.4, 1.0};\n'
                   'Physical Curve("Root") = {1};\nPhysical Surface("A") = {1};\n')
    tri = G.gmsh_mesh(str(geo), size=0.1)
    assert tri.kinds() == ["tri3"]
    quad = G.gmsh_mesh(str(geo), size=0.1, recombine=True)
    assert "quad4" in quad.kinds()
    msh = tmp_path / "sq.msh"
    with pytest.warns(UserWarning, match="corner nodes"):
        o2 = G.gmsh_mesh(str(geo), size=0.1, order=2, save=str(msh))
    assert o2.kinds() == ["tri3"] and msh.is_file()
    with pytest.warns(UserWarning, match="corner nodes"):
        back = G.gmsh_mesh(str(msh))                   # existing .msh is read as is (no remesh)
        scaled = G.gmsh_mesh(str(msh), length_scale=1e-3)
    assert back.info["generated"] is False and "Root" in back.groups
    assert np.ptp(scaled.nodes[:, 1]) == pytest.approx(1e-3)
    fe = G.gmsh_femesh(tri, AL, 0.002, group_thickness={"A": 0.003})
    assert fe.blocks[0].thickness == 0.003 and fe.blocks[0].name == "A"
    with pytest.raises(ValueError, match="no 2-D elements"):
        G.gmsh_femesh(tri, AL, 0.002, group_thickness={"Root": 0.003})


def test_solid_tet10_bar_against_beam_theory():
    res = G.gmsh_modal(os.path.join(GDIR, "box_beam.geo"), 2, material=AL, clamp_group="Root", order=2)
    gm = G.gmsh_mesh(os.path.join(GDIR, "box_beam.geo"), order=2)
    assert gm.kinds() == ["tet10"] and np.allclose(gm.nodes[gm.groups["Root"], 1], 0.0)
    E, rho, b, h = AL["E"], AL["rho"], 0.10, 0.02
    f_eb = 1.875104 ** 2 / (2 * np.pi) * np.sqrt(E * b * h ** 3 / 12 / (rho * b * h))
    assert res.modal.freq_hz[0] == pytest.approx(f_eb, rel=0.02)         # plate effect of b/h = 5: +0.8 %


def test_unsupported_and_bad_input(tmp_path):
    with pytest.raises(ValueError, match="gmsh input"):
        G.gmsh_mesh(os.path.join(ROOT, "examples", "geometry", "wing.stl"))
    geo = tmp_path / "line.geo"
    geo.write_text("Point(1) = {0,0,0};\nPoint(2) = {1,0,0};\nLine(1) = {1,2};\n")
    with pytest.raises(ValueError, match="no surface or volume"):
        G.gmsh_mesh(str(geo))
    wedge = tmp_path / "wedge.geo"
    wedge.write_text('SetFactory("OpenCASCADE");\nRectangle(1) = {0,0,0,1,1};\n'
                     'Extrude {0,0,0.1} { Surface{1}; Layers{2}; Recombine; }\n')
    with pytest.raises(ValueError, match="prism"):
        G.gmsh_mesh(str(wedge), size=0.5)


def test_gmsh_launcher_preflight_and_run(tmp_path):
    from flaap.launcher import run_config
    from flaap.launcher.parser import parse_text
    from flaap.launcher.specs import resolve

    geo = os.path.join(GDIR, "plate_wing.geo").replace("\\", "/")

    def diags(*lines):
        txt = "\n".join(("model modal", "mode analysis", "structure gmsh") + lines) + "\n"
        vals, rep = resolve(parse_text(txt, tmp_path / "config_g.txt"))
        return vals, rep.diagnostics

    def errors(d):
        return {x.code for x in d if x.severity is Severity.ERROR}

    vals, d = diags(f"structure_file {geo}", "clamp_group Root", "group_thickness Skin, 0.003")
    assert not errors(d) and vals["clamp_group"] == "Root" and vals["group_thickness"] == ["Skin", 0.003]
    assert errors(diags(f"structure_file {geo}", "clamp_group Tip")[1]) == {"CFG-PATH-004"}
    assert errors(diags(f"structure_file {geo}", "group_thickness Skin")[1]) == {"CFG-VALUE-007"}
    assert "CFG-VALUE-007" in errors(diags("structure_file ../../examples/geometry/wing.stl")[1])
    assert "CFG-VALUE-001" in errors(diags(f"structure_file {geo}", "mesh_order 3")[1])
    assert any("no clamp_group" in x.message for x in diags(f"structure_file {geo}")[1])
    res = run_config(os.path.join(ROOT, "configs", "modal", "ex9", "config_gmsh_plate.txt"),
                     output_dir=tmp_path / "out", verbose=False,
                     overrides={"plots": False, "mesh_save": str(tmp_path / "p.msh")})
    assert res["status"] == "ok" and (tmp_path / "p.msh").is_file()
    f = [float(x) for x in str(res["summary"]["frequencies_hz"]).replace("[", "").replace("]", "")
         .replace(",", " ").split()] if not isinstance(res["summary"]["frequencies_hz"], list) \
        else res["summary"]["frequencies_hz"]
    np.testing.assert_allclose(f, PLATE_HZ, rtol=2e-5)
