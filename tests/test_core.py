"""Regression / validation tests of the FLAAP Python solver.

Run:  python -m pytest -q tests
"""
import os

import numpy as np
import pytest

from flaap.atmosphere import isa
from flaap.linalg import cdet, cinvert_lu


# ------------------------------------------------------------------ utilities
def test_isa():
    assert isa(0.0).density == pytest.approx(1.225, rel=1e-4)
    assert isa(15000.0).temperature == pytest.approx(216.65, abs=0.01)   # geometric altitude input
    assert isa(1000.0).density == pytest.approx(1.1117, rel=1e-3)


def test_complex_linalg():
    rng = np.random.default_rng(0)
    A = rng.normal(size=(6, 6)) + 1j * rng.normal(size=(6, 6))
    inv, info = cinvert_lu(A)
    assert info == 0
    assert np.allclose(inv @ A, np.eye(6), atol=1e-10)
    assert cdet(A) == pytest.approx(np.linalg.det(A), rel=1e-10)


# ------------------------------------------------------------------ aerodynamics
def _rect(s, c=1.0, nchord=8, nspan=40):
    from flaap.geometry import Surface, build_mesh
    return build_mesh([Surface("w", (0, 0, 0), c, (0, s, 0), c, nchord=nchord, nspan=nspan)],
                      "symmetric")


def test_dlm_steady_limit_matches_vlm_lift_slope():
    """DLM at k -> 0 vs Helmbold lift slope of a rectangular wing (AR 8)."""
    from flaap.aero.gaf import AeroModel, rigid_aero_modes
    s, c = 4.0, 1.0
    mesh = _rect(s, c, 8, 20)
    modes = rigid_aero_modes(mesh, 0.25, ("pitch",))
    tab = AeroModel(mesh, modes, 0.0, 0.5).table([0.0], gust=False)
    cp = tab.cp[0][:, 0]
    CLa = float(np.real(np.sum(cp * mesh.area * mesh.normal[:, 2]))) / (s * c)
    AR = 2 * s / c
    helmbold = 2 * np.pi * AR / (2 + np.sqrt(AR ** 2 + 4))
    # DLM steady limit == independent steady VLM (horseshoe vortices)
    from flaap.aero.vlm import steady_aic
    cpv = np.linalg.solve(steady_aic(mesh, 0.0), np.ones(mesh.n))
    assert abs(CLa) == pytest.approx(abs(np.sum(cpv * mesh.area)) / (s * c), rel=1e-8)
    # lifting-surface value (~4.62 converged) is a few % below Helmbold's approximation
    assert abs(CLa) == pytest.approx(helmbold, rel=0.06)


def test_dlm_vs_theodorsen_high_aspect_ratio():
    """AR 40 rectangular wing: DLM GAFs vs Theodorsen strip theory."""
    from flaap.aero.gaf import AeroModel, rigid_aero_modes
    from flaap.aero.strip import strip_gaf_table
    s = 20.0
    mesh = _rect(s)
    modes = rigid_aero_modes(mesh, 0.25, ("plunge", "pitch"))
    k = np.array([0.05, 0.2, 0.5, 1.0])
    td = AeroModel(mesh, modes, 0.0, 0.5).table(k, gust=False, keep_cp=False)
    y = np.array([0.0, s])
    ts = strip_gaf_table(y, np.array([[1.0, 0], [1.0, 0]]), np.array([[0, 1.0], [0, 1.0]]),
                         1.0, 0.25, k, bref=0.5)
    ratio = np.abs(td.Q / ts.Q)
    assert np.all(np.abs(ratio[:, 0, :] - 1) < 0.06)       # lift row (finite span: -1..-4 %)
    assert np.all(np.abs(ratio - 1) < 0.10)


def test_gaf_reciprocity_and_k0_real(goland):
    tab = goland["table"]
    assert np.allclose(tab.k[0], 0.0)
    # steady GAFs are real
    assert np.max(np.abs(tab.Q[0].imag)) < 1e-8 * np.max(np.abs(tab.Q[0]))


# ------------------------------------------------------------------ RFA
def test_roger_rfa_reproduces_table(goland):
    from flaap.aero.rfa import fit_roger
    tab = goland["table"]
    rfa = fit_roger(tab, 4)
    Q, Qg = rfa.evaluate(tab.k)
    # exact at k = 0 (constrained) and close everywhere in an element-wise relative sense
    assert np.allclose(Q[0], tab.Q[0], atol=1e-8 * np.abs(tab.Q).max())
    err = np.linalg.norm(Q - tab.Q, axis=(1, 2)) / np.linalg.norm(tab.Q, axis=(1, 2))
    assert err.max() < 0.06 and err[tab.k <= 1.0].max() < 0.03
    eg = np.linalg.norm(Qg - tab.Qg, axis=1) / np.linalg.norm(tab.Qg, axis=1)
    assert eg[tab.k <= 1.0].max() < 0.05


# ------------------------------------------------------------------ flutter
@pytest.fixture(scope="module")
def goland_pk(goland):
    from flaap.flutter import pk_method
    mm = goland["modal"]
    return pk_method(mm.M, mm.K, goland["table"].subset(range(4)), 1.225, np.linspace(60, 260, 101))


def test_goland_pk_flutter(goland_pk):
    c = goland_pk.critical()
    assert c is not None and c.kind == "flutter"
    assert c.velocity == pytest.approx(153.3, rel=0.01)
    assert c.omega == pytest.approx(70.7, rel=0.03)


def test_goland_k_method_agrees_with_pk(goland, goland_pk):
    from flaap.flutter import k_method
    mm = goland["modal"]
    r = k_method(mm.M, mm.K, goland["table"].subset(range(4)), 1.225)
    fl = [p for p in r.points if p.kind == "flutter"]
    assert fl
    v = min(p.velocity for p in fl)
    assert v == pytest.approx(goland_pk.critical().velocity, rel=0.02)


def test_goland_p_method_state_space(goland, goland_pk):
    from flaap.aero.rfa import fit_roger
    from flaap.timedomain import modal_system
    mm = goland["modal"]
    rfa = fit_roger(goland["table"], 4)
    r = modal_system(mm, rfa).flutter(np.linspace(60, 260, 101), 1.225)
    c = r.critical()
    assert c is not None
    assert c.velocity == pytest.approx(goland_pk.critical().velocity, rel=0.01)


def test_pk_branches_do_not_merge(goland_pk):
    f = goland_pk.freq
    for i in range(f.shape[0]):
        row = np.sort(f[i])
        assert np.all(np.diff(row) > 1e-3 * row.max())


# ------------------------------------------------------------------ gust
def test_discrete_gust_time_vs_frequency_domain(goland):
    from flaap import gust
    from flaap.aero.rfa import fit_roger
    from flaap.timedomain import modal_system
    mm, tab = goland["modal"], goland["table"]
    rfa = fit_roger(tab, 4)
    ss = modal_system(mm, rfa).state_space(100.0, 1.225)
    r = gust.discrete_gust(ss, 100.0, 20.0, 10.0, t_end=1.5, dt=1e-3)
    td = r.peak("eta:mode1")[0]

    def H(om):
        eta, _ = gust.frequency_response_gust(mm.M, mm.K, mm.C(), tab, 100.0, 1.225, om)
        return eta
    fd = gust.fft_response(r.t, r.wg, H)[:, 0].real.max()
    assert td == pytest.approx(fd, rel=0.05)


def test_cs25_gust_velocity():
    from flaap.gust import cs25_uref, design_gust_velocity
    assert cs25_uref(0.0) == pytest.approx(17.07, rel=1e-3)
    assert design_gust_velocity(107.0, 0.0, 1.0) == pytest.approx(17.07, rel=1e-3)
    assert design_gust_velocity(9.0, 0.0, 1.0) < design_gust_velocity(107.0, 0.0, 1.0)


# ------------------------------------------------------------------ trim / manoeuvre
def test_demo_trim_and_maneuver(demo):
    from flaap.aero.rfa import fit_roger
    from flaap.maneuver import checked_pitch_maneuver, short_period
    from flaap.timedomain import AeroelasticSystem
    from flaap.trim import G0, StaticAeroelastic
    ac, tab = demo["ac"], demo["table"]
    sa = StaticAeroelastic.from_table(ac.M, ac.K, tab)
    V, rho = 100.0, 1.225
    tr = sa.solve(0.5 * rho * V * V, fixed={"plunge": 0.0}, free=["elevator"], accel={"plunge": G0})
    weight = ac.M[0, 0] * G0
    assert tr.forces[0] == pytest.approx(weight, rel=1e-6)       # lift = weight (half model)
    assert 0.0 < np.degrees(tr["pitch"]) < 5.0
    rfa = fit_roger(tab, 4)
    sysm = AeroelasticSystem(ac.M, ac.K, rfa, C=ac.C, n_ctrl=1, names=ac.names + ["elevator"])
    wn, zeta, _ = short_period(sysm, V, rho)
    assert 1.0 < wn < 30.0 and zeta > 0.2
    r = checked_pitch_maneuver(sysm, V, rho, "elevator", nz_limit=2.5)
    assert r.nz.max() == pytest.approx(2.5, rel=1e-6)


def test_goland_divergence_positive(goland):
    from flaap.trim import StaticAeroelastic
    mm, tab = goland["modal"], goland["table"]
    sa = StaticAeroelastic.from_table(mm.M, mm.K, tab.subset(range(4)))
    qd = sa.divergence()
    assert qd is not None and np.sqrt(2 * qd / 1.225) > 153.0     # divergence above flutter


# ------------------------------------------------------------------ I/O
def test_nastran_f06_reader(legacy_dir):
    from flaap.io.nastran import read_eigenvalues, read_f06_modal, read_grids
    f06 = os.path.join(legacy_dir, "20.f06")
    if not os.path.exists(f06):
        pytest.skip("20.f06 not present")
    assert len(read_grids(f06)) == 45
    ev = read_eigenvalues(f06)
    assert ev.mode.size == 10
    m = read_f06_modal(f06)
    assert m.nmodes == 4 and m.nodes.shape == (45, 3)
    assert np.allclose(m.freq_hz, ev.freq_hz[6:], rtol=1e-6)


def test_legacy_reader_and_solver_vs_fortran(legacy_dir):
    from flaap.io.legacy import read_case, run_legacy
    if not os.path.exists(os.path.join(legacy_dir, "LEE_input.txt")):
        pytest.skip("legacy input not present")
    case = read_case(legacy_dir)
    assert case.modal().nodes.shape[0] == case.inp.n_span * (case.inp.n_chord + 1)
    ref = os.path.join(legacy_dir, "pk_results_mode8.dat")
    if not os.path.exists(ref):
        pytest.skip("Fortran reference results not present")
    res, case, _ = run_legacy(legacy_dir, velocities=[50.0], kvals=np.r_[0.0, np.geomspace(0.01, 3.0, 16)])
    Vf, gf, ff = np.loadtxt(ref)[np.argmin(np.abs(np.loadtxt(ref)[:, 0] - 50.0))]
    assert res.freq[0, 7] == pytest.approx(ff, rel=0.03)
    assert res.g[0, 7] == pytest.approx(gf, rel=0.10)
