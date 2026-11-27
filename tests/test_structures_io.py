"""Structural I/O and 3-D FE tests.

* ANSYS modal reader round trip (synthetic samples in examples/ansys/)
* OpenRadioss starter decks (examples/openradioss/) vs analytic solutions
* Radioss /EIG mode files (anim_to_vtk VTK, synthetic samples in examples/openradioss/eig/)
* MITC4 / DKT shells vs Leissa, invariance, rigid-body modes
* tet10 / hex8 solids vs beam theory, rigid-body modes, LSQ rotations
* STL / STEP geometry -> modes (gmsh meshing)
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from flaap.io import ansys as A  # noqa: E402
from flaap.io import openradioss_results as RE  # noqa: E402
from flaap.io.openradioss import radioss_modal, read_radioss  # noqa: E402
from flaap.launcher.diagnostics import Severity  # noqa: E402
from flaap.structures.fem import (FEMesh, dofs_of_nodes, modal_analysis,  # noqa: E402
                                  select_nodes)
from flaap.structures.modal import ModalModel  # noqa: E402
from flaap.structures.plate import plate_mesh_trapezoid  # noqa: E402
from flaap.structures.reference import (LEISSA_CFFF, eb_cantilever_hz, plate_lambda,  # noqa: E402
                                        timoshenko_cantilever_hz)
from flaap.structures.solid import box_mesh, rotations_from_translations  # noqa: E402

ADIR = os.path.join(ROOT, "examples", "ansys")
RDIR = os.path.join(ROOT, "examples", "openradioss")
E_AL, NU_AL, RHO_AL = 70e9, 0.3, 2700.0
E_ST, NU_ST, RHO_ST = 210e9, 0.3, 7850.0


# ====================================================================== ANSYS
@pytest.fixture(scope="module")
def ansys_ref():
    r = np.load(os.path.join(ADIR, "reference_modes.npz"))
    return {k: r[k] for k in r.files}


def _mac(a, b):
    a, b = a.ravel(), b.ravel()
    return (a @ b) ** 2 / ((a @ a) * (b @ b))


def _check_shapes(mm, ref, comps=slice(0, 6), tol=2e-3):
    sh = ref["shapes"]
    for m in range(mm.nmodes):
        a, b = mm.shapes[:, comps, m], sh[:, comps, m]
        assert _mac(a, b) > 0.9999
        assert np.abs(a - b).max() < tol * np.abs(b).max()


def test_ansys_mapdl_listing_roundtrip(ansys_ref):
    out = os.path.join(ADIR, "plate_modal.out")
    np.testing.assert_allclose(A.read_frequencies(out), ansys_ref["freq_hz"], rtol=1e-4)
    mm = A.read_mapdl_modal(out, nodes=os.path.join(ADIR, "plate.cdb"))
    assert mm.nmodes == 6 and mm.name == "ansys" and mm.labels[0] == "mode1"
    np.testing.assert_allclose(mm.freq_hz, ansys_ref["freq_hz"], rtol=1e-4)
    np.testing.assert_allclose(mm.gen_mass, 1.0)
    np.testing.assert_allclose(mm.nodes, ansys_ref["nodes"], atol=1e-12)
    _check_shapes(mm, ansys_ref)
    mm2 = A.read_mapdl_modal(out)                         # NLIST inside the listing
    np.testing.assert_allclose(mm2.nodes, mm.nodes, atol=1e-12)
    nl, cdb = A.read_nlist(out), A.read_cdb_nodes(os.path.join(ADIR, "plate.cdb"))
    assert sorted(nl) == sorted(cdb)
    assert max(np.abs(nl[k] - cdb[k]).max() for k in nl) < 1e-12
    el = A.read_cdb_elements(os.path.join(ADIR, "plate.cdb"))
    assert len(el) == 72 and all(v[0] == 1 and v[1].size == 4 for v in el.values())


def test_ansys_prnsol_dof_wrapped(ansys_ref):
    p = os.path.join(ADIR, "plate_prnsol_dof.out")
    sets = A.read_prnsol(p)
    assert sorted(sets) == [1, 2, 3, 4, 5, 6] and "ROT" in sets[1] and "U" in sets[1]
    np.testing.assert_allclose(A.read_frequencies(p), ansys_ref["freq_hz"], rtol=1e-6)
    mm = A.read_mapdl_modal(p, nodes=os.path.join(ADIR, "plate.cdb"), nmodes=4)
    assert mm.nmodes == 4
    _check_shapes(mm, dict(ansys_ref, shapes=ansys_ref["shapes"][:, :, :4]))


def test_ansys_mm_units_and_gen_mass(ansys_ref):
    mm = A.read_mapdl_modal(os.path.join(ADIR, "plate_modal_mm.out"),
                            length_scale=1e-3, mass_scale=1e3)
    np.testing.assert_allclose(mm.nodes, ansys_ref["nodes"], atol=1e-9)
    np.testing.assert_allclose(mm.gen_mass, 1e-3)
    mm.shapes = mm.shapes / np.sqrt(mm.gen_mass)[None, None, :]
    _check_shapes(mm, ansys_ref)
    g = np.arange(1.0, 7.0)
    mm3 = A.read_mapdl_modal(os.path.join(ADIR, "plate_modal.out"), gen_mass=g, damping=0.01)
    np.testing.assert_allclose(mm3.gen_mass, g)
    np.testing.assert_allclose(mm3.damping, 0.01)
    assert A.read_mapdl_modal(os.path.join(ADIR, "plate_modal.out"), skip_below_hz=10.0).nmodes == 5


def test_ansys_workbench_roundtrip(ansys_ref, tmp_path):
    wb = os.path.join(ADIR, "workbench")
    files = [os.path.join(wb, f"mode_{i}_Z.txt") for i in range(1, 7)]
    mm = A.read_workbench_modal(files, freq_file=os.path.join(wb, "frequencies.txt"))
    np.testing.assert_allclose(mm.freq_hz, ansys_ref["freq_hz"], rtol=1e-4)
    np.testing.assert_allclose(mm.nodes, ansys_ref["nodes"], atol=1e-6)
    _check_shapes(mm, ansys_ref, comps=slice(2, 3), tol=1e-5)
    np.testing.assert_allclose(A.read_ansys_modal(wb).shapes, mm.shapes)
    ref = ModalModel(ansys_ref["nodes"], ansys_ref["shapes"][:, :, :2],
                     2 * np.pi * ansys_ref["freq_hz"][:2], np.ones(2))
    f_mm = A.write_workbench_modal(ref, str(tmp_path / "mm"), unit="mm")
    mm3 = A.read_workbench_modal(f_mm, frequencies=ref.freq_hz)
    np.testing.assert_allclose(mm3.nodes, ref.nodes, atol=1e-6)
    np.testing.assert_allclose(mm3.shapes[:, 2], ref.shapes[:, 2], atol=1e-5 * np.abs(ref.shapes).max())
    f_tot = A.write_workbench_modal(ref, str(tmp_path / "tot"), total=True)
    with pytest.warns(UserWarning, match="unsigned"):
        A.read_workbench_modal(f_tot, frequencies=ref.freq_hz)


def test_ansys_frequency_variants(ansys_ref, tmp_path):
    ref = ModalModel(ansys_ref["nodes"], ansys_ref["shapes"],
                     2 * np.pi * ansys_ref["freq_hz"], np.ones(6))
    p = tmp_path / "solver_only.out"
    A.write_mapdl_listing(ref, str(p), setlist=False, rotations=False)
    np.testing.assert_allclose(A.read_frequencies(str(p)), ref.freq_hz, rtol=1e-9)
    assert np.all(A.read_mapdl_modal(str(p)).shapes[:, 3:] == 0.0)
    csv = tmp_path / "freq.csv"
    csv.write_text('"Mode","Frequency [Hz]"\n1,3.5\n2,14.25\n', encoding="utf-8")
    np.testing.assert_allclose(A.read_frequencies(str(csv)), [3.5, 14.25])
    mm2 = A.read_ansys_modal([str(p)], nodes=os.path.join(ADIR, "plate.cdb"))
    np.testing.assert_allclose(mm2.freq_hz, ref.freq_hz, rtol=1e-9)
    c = tmp_path / "old.cdb"
    A.write_cdb(ref, str(c), fmt="(3i8,6e20.13)")
    old = A.read_cdb_nodes(str(c))
    assert max(np.abs(old[i + 1] - ref.nodes[i]).max() for i in range(len(old))) < 1e-12


def test_ansys_rst_missing_package_message():
    import importlib.util
    try:
        installed = importlib.util.find_spec("ansys.mapdl.reader") is not None
    except ModuleNotFoundError:
        installed = False
    if installed:
        pytest.skip("ansys-mapdl-reader installed (no sample .rst to test)")
    with pytest.raises(ImportError, match="ansys-mapdl-reader"):
        A.read_rst_modal("none.rst")


# ====================================================================== OpenRadioss
def _plate_lam(mm, n=5):
    return plate_lambda(mm.omega[:n], 1.0, 0.01, E_AL, NU_AL, RHO_AL)


@pytest.mark.parametrize("deck", ["plate_cantilever_0000.rad", "plate_tria_0000.rad"])
def test_radioss_plate_vs_leissa(deck):
    mm = radioss_modal(os.path.join(RDIR, deck), nmodes=5)
    lam = _plate_lam(mm)
    np.testing.assert_allclose(lam, LEISSA_CFFF, rtol=0.012)
    assert mm.node_ids.size == 289 and mm.fe.elements["mesh"].mass() == pytest.approx(27.0, rel=1e-9)


def test_radioss_deck_parsing_units():
    d = read_radioss(os.path.join(RDIR, "plate_cantilever_0000.rad"))
    assert d.units == ("kg", "mm", "ms")
    assert np.ptp(d.nodes[:, 0]) == pytest.approx(1.0)               # mm -> m
    m = d.mats[1]
    assert m.E == pytest.approx(E_AL) and m.rho == pytest.approx(RHO_AL) and m.nu == 0.3
    t = read_radioss(os.path.join(RDIR, "plate_tria_0000.rad"))
    assert t.units == ("Mg", "mm", "s") and t.mats[1].E == pytest.approx(E_AL)
    assert len(t.grnod[1]) == 17                                     # /GRNOD/BOX root nodes
    b = read_radioss(os.path.join(RDIR, "beam_tetra10_0000.rad"))
    (kind, _), (ids, conn, _) = next(iter(b.elements.items()))
    assert kind == "tet10" and conn.shape == (288, 10)


@pytest.mark.parametrize("deck", ["beam_brick_0000.rad", "beam_tetra10_0000.rad"])
def test_radioss_beam_vs_theory(deck):
    mm = radioss_modal(os.path.join(RDIR, deck), nmodes=4)
    L, b, h = 1.0, 0.06, 0.04
    fz = timoshenko_cantilever_hz(L, E_ST, NU_ST, b * h ** 3 / 12, RHO_ST, b * h, 2)
    fy = timoshenko_cantilever_hz(L, E_ST, NU_ST, h * b ** 3 / 12, RHO_ST, b * h, 2)
    ref = np.sort(np.r_[fz, fy])
    np.testing.assert_allclose(mm.freq_hz, ref, rtol=0.02)


def test_radioss_skew_bc(tmp_path):
    """Tip prop in a /SKEW frame == global prop of the unrotated beam."""
    sys.path.insert(0, RDIR)
    import make_decks as MD
    Lx, b, h = 1.0, 0.06, 0.04
    nodes, hexes = box_mesh(Lx, b, h, 12, 2, 2, "hex8")
    root = np.where(np.abs(nodes[:, 0]) < 1e-9)[0] + 1
    tip = np.where(np.abs(nodes[:, 0] - Lx) < 1e-9)[0] + 1
    c, s = np.cos(np.deg2rad(30.0)), np.sin(np.deg2rad(30.0))
    R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    def deck(path, X, skew):
        L = MD.header("skewbeam", ("kg", "m", "s")) + MD.nodes_block(X, 1.0)
        L += ["/BRICK/1"] + [f"{e + 1:10d}" + "".join(f"{n + 1:10d}" for n in q)
                             for e, q in enumerate(hexes)]
        L += MD.part(1, "beam", 1, 1) + MD.prop_solid(1, "s") + MD.law1(1, "steel", RHO_ST, E_ST, NU_ST)
        L += MD.grnod_node(1, "root", root) + MD.grnod_node(2, "tip", tip)
        L += MD.bcs(1, "root", 1, "111 000")
        if skew:
            L += ["/SKEW/FIX/7", "beam frame", MD._r(0) * 3, MD._r(c) + MD._r(s) + MD._r(0),
                  MD._r(-s) + MD._r(c) + MD._r(0)]
        L += MD.bcs(2, "tip prop", 2, "010 000", skew=7 if skew else 0)
        MD.write(str(path), L)

    deck(tmp_path / "a_0000.rad", nodes, False)
    deck(tmp_path / "b_0000.rad", nodes @ R.T, True)
    fa = radioss_modal(str(tmp_path / "a_0000.rad"), 4)
    fb = radioss_modal(str(tmp_path / "b_0000.rad"), 4)
    np.testing.assert_allclose(fb.freq_hz, fa.freq_hz, rtol=1e-8)
    # mode shapes are the rotated ones (translations)
    for m in range(4):
        assert _mac(fb.shapes[:, :3, m], fa.shapes[:, :3, m] @ R.T) > 0.9999


# ====================================================================== Radioss /EIG results
EDIR = os.path.join(RDIR, "eig")
PLATE = os.path.join(RDIR, "plate_cantilever_0000.rad")


@pytest.fixture(scope="module")
def plate_modes():
    return radioss_modal(PLATE, nmodes=6)


def _same_modes(mm, ref, rtol_f=1e-5, mac=0.99999):
    np.testing.assert_allclose(mm.freq_hz, ref.freq_hz[:mm.nmodes], rtol=rtol_f)
    for m in range(mm.nmodes):
        assert _mac(mm.shapes[:, :, m], ref.shapes[:, :, m]) > mac
        a = mm.shapes[:, :, m] / np.sqrt(mm.gen_mass[m])               # unit-mass shapes, same sign
        b = ref.shapes[:, :, m] / np.sqrt(ref.gen_mass[m])
        assert np.abs(a - b).max() < 1e-4 * np.abs(b).max()


def test_radioss_eig_sample_set(plate_modes):
    """examples/openradioss/eig (synthetic anim_to_vtk files, kg-mm-ms) + starter deck."""
    mm = RE.read_radioss_eigen(EDIR, deck=PLATE)
    assert mm.nmodes == 6 and mm.labels == [f"mode{k}" for k in range(1, 7)]
    assert mm.node_ids.size == 289 and mm.fe.elements["mesh"].mass() == pytest.approx(27.0)
    assert np.ptp(mm.nodes[:, 0]) == pytest.approx(1.0)                # mm -> m
    # unit generalised mass in kg-mm-ms == 1e-6 kg m^2; deck mass matrix agrees
    np.testing.assert_allclose(mm.gen_mass, 1e-6, rtol=1e-5)
    _same_modes(mm, plate_modes)
    # same result through one mode file / the explicit file list, without deck (units given)
    one = RE.read_radioss_eigen(os.path.join(EDIR, "plate_cantileverA003.vtk"), units=("kg", "mm", "ms"))
    assert one.nmodes == 6 and one.fe is None
    np.testing.assert_allclose(one.nodes, mm.nodes, atol=1e-12)
    np.testing.assert_allclose(one.shapes, mm.shapes, atol=1e-12)
    assert np.all(one.gen_mass == pytest.approx(1e-6))
    # no unit information: solver units taken as SI (frequency = TIME = cycles per ms)
    si = RE.read_radioss_eigen(sorted(RE.find_radioss_mode_files(EDIR))[:3], skip_below_hz=0.0)
    np.testing.assert_allclose(si.freq_hz, plate_modes.freq_hz[:3] * 1e-3, rtol=1e-5)
    d = RE.read_radioss_vtk(os.path.join(EDIR, "plate_cantileverA001.vtk"))
    assert len(d["cells"]) == 256 and set(d["cell_types"]) == {9}
    assert set(d["vectors"]) == {"Displacement", "Rotational_DOF"} and d["cell_scalars"]["PART_ID"][0] == 1


def test_radioss_eig_tonne_mm_s_roundtrip(tmp_path):
    deck = os.path.join(RDIR, "plate_tria_0000.rad")                 # /BEGIN Mg mm s
    ref = radioss_modal(deck, nmodes=4)
    RE.write_radioss_eigen_vtk(ref, str(tmp_path), "tria", deck=deck)
    head = (tmp_path / "triaA001.vtk").read_text().splitlines()
    assert head[5] == "TIME 1 1 double" and float(head[6]) == pytest.approx(ref.freq_hz[0], rel=1e-5)
    a = RE.read_radioss_eigen(str(tmp_path), deck=deck)
    b = RE.read_radioss_eigen(str(tmp_path), units=("Mg", "mm", "s"))
    for mm in (a, b):
        np.testing.assert_allclose(mm.gen_mass, 1e3 * 1e-6, rtol=1e-4)  # Ms Ls^2
        _same_modes(mm, ref)
    with pytest.raises(ValueError, match="needs a starter deck"):
        RE.read_radioss_eigen(str(tmp_path), gen_mass="deck")


def test_radioss_eig_listing_rigid_and_rotations(tmp_path, plate_modes):
    ref = plate_modes
    fs = ref.freq_hz[:4] * 1e-3                                       # 1/ms
    RE.write_radioss_eigen_vtk(ref.__class__(ref.nodes, ref.shapes[:, :, :4], ref.omega[:4], ref.gen_mass[:4]),
                              str(tmp_path), "plate", deck=PLATE, node_ids=ref.node_ids,
                              frequencies_solver=[0.0, *fs[1:]], rotations=False)
    # rigid/static mode (TIME = 0) skipped; shell rotations missing -> zero + warning
    with pytest.warns(UserWarning, match="Rotational_DOF"):
        mm = RE.read_radioss_eigen(str(tmp_path), deck=PLATE)
    assert mm.labels == ["mode2", "mode3", "mode4"] and np.all(mm.shapes[:, 3:, :] == 0)
    np.testing.assert_allclose(mm.freq_hz, ref.freq_hz[1:4], rtol=1e-5)
    # listings: FREQ column (picked by header position), EIGENVALUE only, plain two-column table
    w = 2 * np.pi * fs
    lst = tmp_path / "plate_0001.out"
    lst.write_text("  EIGEN ANALYSIS\n\n   MODE     EIGENVALUE      FREQUENCY      GEN. MASS\n"
                   + "".join(f"{i + 1:7d} {w[i] ** 2:14.6E} {fs[i]:14.6E} {1.0:14.6E}\n" for i in range(4))
                   + "\n   NORMAL TERMINATION\n")
    got = RE.read_radioss_frequencies(str(lst))
    np.testing.assert_allclose([got[k] for k in (1, 2, 3, 4)], fs, rtol=1e-6)
    eig = tmp_path / "eig.txt"
    eig.write_text("MODE  EIGENVALUE\n" + "".join(f"{i + 1} {w[i] ** 2:.8e}\n" for i in range(4)))
    np.testing.assert_allclose(list(RE.read_radioss_frequencies(str(eig)).values()), fs, rtol=1e-6)
    tab = tmp_path / "freq.txt"
    tab.write_text("# mode  f\n" + "".join(f"{i + 1}, {2 * fs[i]:.8e}\n" for i in range(4)))
    assert RE.read_radioss_frequencies(str(tab))[2] == pytest.approx(2 * fs[1])
    # listing as source: VTK files found beside it (run name), listing frequencies override TIME
    with pytest.warns(UserWarning):
        a = RE.read_radioss_eigen(str(lst), deck=PLATE)
        b = RE.read_radioss_eigen(str(tmp_path), deck=PLATE, listing=str(tab), nmodes=2)
    assert a.nmodes == 4 and a.freq_hz[0] == pytest.approx(ref.freq_hz[0], rel=1e-6)
    np.testing.assert_allclose(b.freq_hz, 2 * ref.freq_hz[:2], rtol=1e-6)


def test_radioss_eig_solid_rotations(tmp_path):
    deck = os.path.join(RDIR, "beam_brick_0000.rad")
    ref = radioss_modal(deck, nmodes=3)                               # lsq rotations
    RE.write_radioss_eigen_vtk(ref, str(tmp_path), "beam", deck=deck, rotations=False)
    mm = RE.read_radioss_eigen(os.path.join(str(tmp_path), "beamA002.vtk"), deck=deck)
    assert mm.nmodes == 3
    for m in range(3):
        assert _mac(mm.shapes[:, 3:, m], ref.shapes[:, 3:, m]) > 0.9999   # rotations re-estimated


def test_radioss_eig_launcher_preflight(tmp_path):
    from flaap.launcher.parser import parse_text
    from flaap.launcher.specs import resolve

    def check(*lines):
        txt = "\n".join(("model modal", "mode analysis", "structure openradioss_eig") + lines) + "\n"
        _, rep = resolve(parse_text(txt, tmp_path / "config_x.txt"))
        return [d for d in rep.diagnostics if d.code != "CFG-DEFAULT-001" or "openradioss" in d.message]

    def codes(diags, sev=None):
        return {d.code for d in diags if sev is None or d.severity is sev}

    e, deck = EDIR.replace("\\", "/"), PLATE.replace("\\", "/")
    assert not codes(check(f"structure_file {e}", f"bulk_file {deck}"), Severity.ERROR)
    (tmp_path / "empty").mkdir()
    bad = check("structure_file empty", "bulk_file missing_0000.rad")
    assert codes(bad, Severity.ERROR) == {"CFG-PATH-001", "CFG-PATH-002"}
    assert [d.field for d in bad if d.code == "CFG-PATH-001"] == ["bulk_file"]
    (tmp_path / "run_0001.out").write_text("x\n")
    assert "CFG-PATH-002" in codes(check("structure_file run_0001.out"))
    assert "CFG-PATH-002" not in codes(check("structure_file run_0001.out", f"mode_files {e}"))
    mg = check(f"structure_file {e}", "radioss_units Mg, mm, s")
    assert any(d.severity is Severity.WARNING and "milligram" in d.message for d in mg)
    assert "CFG-VALUE-001" in codes(check(f"structure_file {e}", "radioss_units kg, furlong, s"), Severity.ERROR)
    assert any("without bulk_file" in d.message for d in check(f"structure_file {e}"))


# ====================================================================== shells
def _square_plate(n, kind):
    nodes, quads = plate_mesh_trapezoid((0, 0), 1.0, (0, 1.0), 1.0, n, n)
    nodes = np.c_[nodes, np.zeros(len(nodes))] if nodes.shape[1] == 2 else nodes
    conn = quads if kind == "quad4" else np.vstack([quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]])
    return FEMesh(nodes).add(kind, conn, dict(E=E_AL, nu=NU_AL, rho=RHO_AL), 0.01)


@pytest.mark.parametrize("kind", ["quad4", "tri3"])
def test_shell_plate_vs_leissa(kind):
    mesh = _square_plate(12, kind)
    root = select_nodes(mesh.nodes, below=("y", 0.0))
    mm = modal_analysis(mesh, 5, fixed=dofs_of_nodes(root))
    np.testing.assert_allclose(_plate_lam(mm), LEISSA_CFFF, rtol=0.02)


def test_shell_rotation_invariance_and_rigid_modes():
    mesh = _square_plate(6, "tri3")
    root = select_nodes(mesh.nodes, below=("y", 0.0))
    f0 = modal_analysis(mesh, 4, fixed=dofs_of_nodes(root)).freq_hz
    ang = np.deg2rad([25.0, -40.0, 70.0])
    from scipy.spatial.transform import Rotation
    R = Rotation.from_euler("xyz", ang).as_matrix()
    for kind in ("tri3", "quad4"):
        m = _square_plate(6, kind)
        f_ref = modal_analysis(m, 4, fixed=dofs_of_nodes(select_nodes(m.nodes, below=("y", 0.0)))).freq_hz
        m.nodes = m.nodes @ R.T + [3.0, -1.0, 2.0]
        f_rot = modal_analysis(m, 4, fixed=dofs_of_nodes(root)).freq_hz
        np.testing.assert_allclose(f_rot, f_ref, rtol=1e-6)
    assert f0[0] > 0
    for kind in ("tri3", "quad4"):           # free-free: 6 rigid-body modes
        f = modal_analysis(_square_plate(6, kind), 8).freq_hz
        assert np.all(f[:6] < 1e-2) and f[6] > 1.0


# ====================================================================== solids
def _beam(kind, nx, ny=2, nz=2, L=1.0, b=0.06, h=0.04):
    nodes, conn = box_mesh(L, b, h, nx, ny, nz, kind)
    return FEMesh(nodes).add(kind, conn, dict(E=E_ST, nu=NU_ST, rho=RHO_ST))


@pytest.mark.parametrize("kind,nx", [("hex8", 20), ("tet10", 10)])
def test_solid_beam_vs_theory(kind, nx):
    L, b, h = 1.0, 0.06, 0.04
    mesh = _beam(kind, nx)
    assert mesh.mass() == pytest.approx(RHO_ST * L * b * h, rel=1e-10)
    root = select_nodes(mesh.nodes, below=("x", 0.0))
    mm = modal_analysis(mesh, 3, fixed=dofs_of_nodes(root, "123"))
    eb = np.sort(np.r_[eb_cantilever_hz(L, E_ST, b * h ** 3 / 12, RHO_ST, b * h, 2),
                       eb_cantilever_hz(L, E_ST, h * b ** 3 / 12, RHO_ST, b * h, 1)])
    ti = np.sort(np.r_[timoshenko_cantilever_hz(L, E_ST, NU_ST, b * h ** 3 / 12, RHO_ST, b * h, 2),
                       timoshenko_cantilever_hz(L, E_ST, NU_ST, h * b ** 3 / 12, RHO_ST, b * h, 1)])
    np.testing.assert_allclose(mm.freq_hz, ti, rtol=0.02)
    np.testing.assert_allclose(mm.freq_hz, eb, rtol=0.025)
    # rotations at the tip = bending slope (LSQ estimate) - finite and nonzero
    assert np.all(np.isfinite(mm.shapes)) and np.abs(mm.shapes[:, 3:, 0]).max() > 0


@pytest.mark.parametrize("kind", ["tet10", "hex8"])
def test_solid_free_free_rigid_modes(kind):
    f = modal_analysis(_beam(kind, 6, 2, 2, 0.5), 8).freq_hz
    assert np.all(f[:6] < 1e-2) and f[6] > 10.0


def test_rotations_from_translations_rigid():
    nodes, conn = box_mesh(1.0, 0.3, 0.2, 4, 2, 2, "tet10")
    th = np.array([0.01, -0.02, 0.03])
    u = np.cross(th, nodes - [0.2, 0.1, 0.0])                 # small rigid rotation
    rot = rotations_from_translations(nodes, [conn], u[:, :, None])
    np.testing.assert_allclose(rot[:, :, 0], np.broadcast_to(th, rot[:, :, 0].shape), atol=1e-12)


# ====================================================================== geometry
def test_stl_box_beam_tet10_vs_theory(tmp_path):
    pytest.importorskip("gmsh")
    trimesh = pytest.importorskip("trimesh")
    from flaap.structures.geometry_modal import geometry_modal
    L, b, h = 1.0, 0.06, 0.04
    tm = trimesh.creation.box(extents=(L, b, h))
    tm.apply_translation((L / 2, 0, 0))
    p = str(tmp_path / "beam.stl")
    tm.export(p)
    r = geometry_modal(p, 3, kind="solid", element="tet10", size=0.03,
                       material=dict(E=E_ST, nu=NU_ST, rho=RHO_ST), clamp=dict(below=("x", 0.0)))
    ti = np.sort(np.r_[timoshenko_cantilever_hz(L, E_ST, NU_ST, b * h ** 3 / 12, RHO_ST, b * h, 2),
                       timoshenko_cantilever_hz(L, E_ST, NU_ST, h * b ** 3 / 12, RHO_ST, b * h, 1)])
    np.testing.assert_allclose(r.modal.freq_hz, ti, rtol=0.02)
    assert r.info["mass"] == pytest.approx(RHO_ST * L * b * h, rel=1e-6)
    assert r.surface.shape[1] == 3 and r.clamped.size > 0


def test_step_box_beam_and_stl_plate_shell(tmp_path):
    gmsh = pytest.importorskip("gmsh")
    from flaap.structures.geometry_modal import geometry_modal, write_stl
    L, b, h = 1.0, 0.06, 0.04
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("b")
    gmsh.model.occ.addBox(0, -b / 2 * 1e3, -h / 2 * 1e3, L * 1e3, b * 1e3, h * 1e3)   # mm
    gmsh.model.occ.synchronize()
    ps = str(tmp_path / "beam.step")
    gmsh.write(ps)
    gmsh.finalize()
    r = geometry_modal(ps, 2, kind="solid", size=0.03, material="steel",
                       clamp=dict(below=("x", 0.0)))
    assert np.ptp(r.mesh.nodes[:, 0]) == pytest.approx(L)              # STEP mm -> m
    fz = timoshenko_cantilever_hz(L, 200e9, 0.29, b * h ** 3 / 12, 7850.0, b * h, 1)[0]
    assert r.modal.freq_hz[0] == pytest.approx(fz, rel=0.02)
    # open flat STL (2 triangles) -> gmsh remesh -> DKT shells vs Leissa
    pp = str(tmp_path / "plate.stl")
    write_stl(pp, np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], float),
              np.array([[0, 1, 2], [0, 2, 3]]))
    r = geometry_modal(pp, 5, kind="shell", thickness=0.01, size=0.08, remesh=True,
                       material=dict(E=E_AL, nu=NU_AL, rho=RHO_AL), clamp=dict(below=("y", 0.0)))
    np.testing.assert_allclose(_plate_lam(r.modal), LEISSA_CFFF, rtol=0.02)


def test_stl_wing_modes(tmp_path):
    pytest.importorskip("gmsh")
    trimesh = pytest.importorskip("trimesh")
    from flaap.structures.geometry_modal import (geometry_modal, naca_wing_surface,
                                                 planform_nodes, surface_triangulation, write_stl)
    V, F = naca_wing_surface(1.2, 0.4, 0.25, 10.0, 0.12, 16, 8)
    tm = trimesh.Trimesh(V, F, process=False)
    assert tm.is_watertight and tm.is_winding_consistent and tm.volume > 0
    p = str(tmp_path / "wing.stl")
    write_stl(p, V, F)
    solid = geometry_modal(p, 4, kind="solid", element="tet10", size=0.05, remesh=False,
                           material="aluminium_7075", clamp=dict(below=("y", 0.0)))
    shell = geometry_modal(p, 4, kind="shell", thickness=0.002, material="aluminium_7075",
                           clamp=dict(below=("y", 0.0)))
    assert solid.info["mass"] == pytest.approx(tm.volume * 2810.0, rel=1e-6)
    for r in (solid, shell):
        f = r.modal.freq_hz
        assert np.all(np.diff(f) >= 0) and 1.0 < f[0] < 100.0
        # first mode: flapwise bending (tip uz dominant)
        tip = r.mesh.nodes[:, 1] > 1.15
        u = np.abs(r.modal.shapes[tip, :3, 0]).mean(0)
        assert u[2] > 3 * max(u[0], u[1])
        np.testing.assert_allclose(r.modal.M, np.eye(4), atol=1e-8)
        tri = surface_triangulation(r)
        assert tri.max() < r.mesh.nn
    m = planform_nodes(solid.modal.nodes, 10, 4)
    assert 30 <= m.sum() <= 40
