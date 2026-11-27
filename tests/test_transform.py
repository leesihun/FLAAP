"""Structure re-orientation (mirror / rotate / translate) and symmetric / antisymmetric aero images."""
import numpy as np
import pytest

from flaap.geometry import Surface, build_mesh
from flaap.launcher import parse_text, resolve
from flaap.launcher.builder import build_case
from flaap.structures.transform import (Transform, rotation_matrix, transform_from_values, transform_modal,
                                        transform_surface)

BASE = """model modal
structure goland
nmodes 3
nchord 4
nspan 6
"""


def _case(extra=""):
    pc = parse_text(BASE + extra)
    vals, rep = resolve(pc)
    assert not rep.has_errors(), [d.message for d in rep.diagnostics]
    return build_case(vals, pc.base_dir, log=lambda *_: None), vals


def test_rotation_matrix_exact():
    R = rotation_matrix("z", 90)
    assert np.array_equal(R @ [1, 0, 0], [0, 1, 0])
    assert np.allclose(rotation_matrix("x", 90) @ [0, 1, 0], [0, 0, 1])


def test_parse_transform_errors():
    with pytest.raises(ValueError):
        transform_from_values({"struct_mirror": "q"})
    with pytest.raises(ValueError):
        transform_from_values({"struct_rotate": ["z"]})
    with pytest.raises(ValueError):
        transform_from_values({"struct_translate": [1, 2]})
    T = transform_from_values({"struct_mirror": "y", "struct_rotate": ["z", 90, "x", 180], "struct_translate": [0, 1, 0]})
    assert T.mirror == "y" and len(T.rotations) == 2 and T.translate == (0, 1, 0)


def test_mirror_modal_pseudovector():
    (case, _), (case_m, _) = _case(), _case("struct_mirror y\n")
    a, b = case.modal, case_m.modal
    assert np.allclose(b.nodes[:, 1], -a.nodes[:, 1])
    # translation uy flips, ux/uz do not; rotation ry unchanged, rx/rz flip
    assert np.allclose(b.shapes[:, 0], a.shapes[:, 0]) and np.allclose(b.shapes[:, 1], -a.shapes[:, 1])
    assert np.allclose(b.shapes[:, 3], -a.shapes[:, 3]) and np.allclose(b.shapes[:, 4], a.shapes[:, 4])
    assert np.allclose(b.shapes[:, 5], -a.shapes[:, 5])
    assert np.allclose(b.omega, a.omega) and np.allclose(b.gen_mass, a.gen_mass)


def test_mirror_preset_surface_swaps_root_tip_and_controls():
    from flaap.geometry import ControlSurface
    s = Surface("w", (0, 0, 0), 1.0, (0.3, 2.0, 0), 0.5, 4, 6)
    s.control_surfaces.append(ControlSurface("ail", 0.8, (0.6, 1.0)))
    m = transform_surface(s, Transform(mirror="y"))
    assert m.root_le[1] == pytest.approx(-2.0) and m.tip_le[1] == pytest.approx(0.0)
    assert m.root_chord == 0.5 and m.tip_chord == 1.0
    assert m.control_surfaces[0].eta == pytest.approx((0.0, 0.4))
    with pytest.raises(ValueError):
        transform_surface(s, Transform(rotations=[("z", 90.0)]))


def test_mirrored_wing_gives_identical_gaf():
    """A mirrored (left) half wing with the symmetric image is the same aircraft: identical Q(k)."""
    from flaap.aero.gaf import AeroModel
    out = []
    for extra in ("", "struct_mirror y\n"):
        case, _ = _case(extra)
        am = AeroModel(case.mesh, case.modes, 0.0, case.bref)
        out.append(am.table([0.0, 0.3, 1.0], gust=False, keep_cp=False).Q)
    assert np.allclose(out[0], out[1], rtol=1e-8, atol=1e-10)


def test_rotated_plate_equals_reference():
    """A plate meshed along x (span along x) rotated by +90 deg about z has the reference modes."""
    ref = parse_text("model modal\nstructure plate\nplate_root_le 0,0,0\nplate_root_chord 0.4\n"
                     "plate_tip_le 0.1,1.0,0\nplate_tip_chord 0.3\nnmodes 3\nplate_nchord 4\nplate_nspan 6\n")
    # same plate but modelled with the span along -x: rotate by +90 deg about z maps (x, y) -> (-y, x);
    # build it by hand from the reference modal model instead (the FE builder is axis-aligned)
    vals, _ = resolve(ref)
    case = build_case(vals, ref.base_dir, log=lambda *_: None)
    mm = case.modal
    T = Transform(rotations=[("z", 90.0)])
    rot = transform_modal(mm, T)
    back = transform_modal(rot, Transform(rotations=[("z", -90.0)]))
    assert np.allclose(back.nodes, mm.nodes, atol=1e-12) and np.allclose(back.shapes, mm.shapes, atol=1e-12)
    assert np.allclose(rot.nodes[:, 0], -mm.nodes[:, 1]) and np.allclose(rot.nodes[:, 1], mm.nodes[:, 0])
    assert np.allclose(rot.shapes[:, 2], mm.shapes[:, 2])        # uz unaffected by a rotation about z


def test_transform_rules_in_preflight():
    pc = parse_text(BASE + "struct_rotate z, 90\n")
    _, rep = resolve(pc)
    assert any(d.code == "CFG-STRUCT-003" for d in rep.diagnostics)      # preset chord cannot turn
    pc = parse_text(BASE + "struct_mirror q\n")
    _, rep = resolve(pc)
    assert rep.has_errors()
    pc = parse_text(BASE + "struct_mirror y\nstruct_translate 0, 0, 1\n")
    _, rep = resolve(pc)
    assert not rep.has_errors()


def test_half_model_symmetry_warning():
    pc = parse_text("model aero\nmode gaf\nsurface1 wing, 0,-3,0, 1, 0,0,0, 1, 4, 6\nsymmetry symmetric\n")
    _, rep = resolve(pc)
    assert any(d.code == "CFG-SYM-001" for d in rep.diagnostics)


@pytest.mark.parametrize("sym,wash,moment", [("symmetric", lambda y: np.ones_like(y), False),
                                              ("antisymmetric", lambda y: y, True)])
def test_half_model_image_equals_full_model(sym, wash, moment):
    """Half model + image == full-span model: heave (symmetric) and roll (antisymmetric) loads."""
    from flaap.aero.dlm import DLM
    half = build_mesh([Surface("w", (0, 0, 0), 1.0, (0.2, 3.0, 0), 0.6, 3, 6)], sym)
    left = Surface("wl", (0.2, -3.0, 0), 0.6, (0, 0, 0), 1.0, 3, 6)
    right = Surface("wr", (0, 0, 0), 1.0, (0.2, 3.0, 0), 0.6, 3, 6)
    full = build_mesh([left, right], "none")

    def load(mesh, k):
        w = wash(mesh.colloc[:, 1]).astype(complex)
        cp = np.linalg.solve(DLM(mesh, 0.3).aic(k, 1.0), w)
        wt = mesh.area * (mesh.colloc[:, 1] if moment else 1.0)
        return np.sum(cp * wt)
    for k in (0.0, 0.4, 1.5):
        assert 2 * load(half, k) == pytest.approx(load(full, k), rel=3e-2)
