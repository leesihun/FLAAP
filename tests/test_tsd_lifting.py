"""Lifting (angle of attack + camber) nonlinear TSD mean flow."""
import numpy as np
import pytest

from flaap.aero.tsd import TSDSolver, naca4_thickness, rigid_mode_funcs, tsd_gaf_table
from flaap.aero.tsd_lifting import (LiftingTSDSolver, naca4_camber, naca4_section,
                                    shock_positions)
from flaap.geometry import Surface

COARSE = dict(nchord=16, nspan=8, nz=18, n_up=10, n_down=10, n_out=8)


@pytest.fixture(scope="module")
def wing():
    return Surface("w", (0, 0, 0), 1.0, (0, 4.0, 0), 1.0, nchord=10, nspan=10)


def test_naca4_camber_and_section():
    zc = naca4_camber(0.02, 0.4)
    xi = np.linspace(0, 1, 201)
    assert zc(xi).max() == pytest.approx(0.02, rel=1e-4)
    assert xi[np.argmax(zc(xi))] == pytest.approx(0.4, abs=0.01)
    h = 1e-6
    x = np.array([0.1, 0.3, 0.5, 0.9])
    np.testing.assert_allclose(zc.slope(x), (zc(x + h) - zc(x - h)) / (2 * h), atol=1e-7)
    t, c = naca4_section("NACA2412")
    assert c.m == 0.02 and c.p == 0.4 and t(0.3) > 0
    assert naca4_section("0012")[1] is None
    with pytest.raises(ValueError):
        naca4_section("2x12")


def test_shock_detector():
    xi = np.linspace(0, 1, 41)[None, :]
    cp = np.where(xi < 0.6, -0.9, -0.1)
    xs, dj = shock_positions(xi, cp, -0.5)
    assert xs[0] == pytest.approx(0.6125, abs=0.03) and dj[0] == pytest.approx(0.8)
    xs, dj = shock_positions(xi, np.full_like(xi, -0.2), -0.5)
    assert np.isnan(xs[0]) and dj[0] == 0


def test_zero_incidence_reduces_to_symmetric_solver(wing):
    """alpha = 0, no camber: U = L, the steady flow and the GAFs equal TSDSolver's."""
    th = naca4_thickness(0.12)
    s0 = TSDSolver(wing, 0.8, th, **COARSE)
    s1 = LiftingTSDSolver(wing, 0.8, th, **COARSE)
    r0, r1 = s0.steady(), s1.steady()
    assert r1.converged and r1.n_supersonic == r0.n_supersonic > 0
    np.testing.assert_allclose(r1.phi, r0.phi, atol=1e-12)
    np.testing.assert_allclose(r1.phi_lower, r1.phi, atol=1e-12)
    assert abs(r1.CL) < 1e-12
    funcs, names = rigid_mode_funcs(0.25)
    ta = tsd_gaf_table(wing, 0.8, [0.0, 0.3], 0.5, funcs, names, thickness=th, solver=s0)
    tb = tsd_gaf_table(wing, 0.8, [0.0, 0.3], 0.5, funcs, names, thickness=th, solver=s1)
    np.testing.assert_allclose(tb.Q, ta.Q, atol=1e-10 * np.abs(ta.Q).max())


def test_small_incidence_matches_linear_lift(wing):
    """Flat plate, M = 0.5: CL/alpha equals the linear (k = 0) solution, and is
    close to the Helmbold estimate for the aspect ratio 8 wing."""
    lin = TSDSolver(wing, 0.5, None, **COARSE)
    a = np.radians(1.0)
    phi = lin.solve_normalwash(0.0, np.full((lin.nspan * lin.nchord, 1), -a))
    jump = 2 * phi.reshape(*lin.shape, -1)[lin.iw1, :lin.nspan, 0, 0].real
    cl_lin = 2 * jump * lin.gmap / lin.b[:lin.nspan]
    r = LiftingTSDSolver(wing, 0.5, None, alpha_deg=1.0, **COARSE).steady()
    assert r.converged
    np.testing.assert_allclose(r.cl, cl_lin, rtol=2e-3)
    AR, beta = 8.0, np.sqrt(1 - 0.25)
    cla = 2 * np.pi * AR / (2 + np.sqrt(AR ** 2 * beta ** 2 + 4))
    assert r.CL / a == pytest.approx(cla, rel=0.08)


def test_incidence_sign_swaps_surfaces_and_moves_shock(wing):
    th = naca4_thickness(0.12)
    rs = {a: LiftingTSDSolver(wing, 0.8, th, alpha_deg=a, **COARSE).steady() for a in (-1.0, 0.0, 1.0)}
    assert all(r.converged for r in rs.values())
    rp, rm, r0 = rs[1.0], rs[-1.0], rs[0.0]
    np.testing.assert_allclose(rp.cp, rm.cp_lower, atol=1e-10)
    np.testing.assert_allclose(rp.cp_lower, rm.cp, atol=1e-10)
    assert rp.CL == pytest.approx(-rm.CL, rel=1e-10) and rp.CL > 0
    # nonlinear: transonic lift slope above the incompressible-like linear value
    assert rp.n_supersonic > r0.n_supersonic > rp.n_supersonic_lower
    # upper shock aft of the zero-incidence one and stronger; lower one weaker
    assert rp.shock_upper[0] >= r0.shock_upper[0]
    assert rp.jump_upper[0] > r0.jump_upper[0] > rp.jump_lower[0]


def test_pseudo_transient_and_warm_start_reach_newton_solution(wing):
    """The PTC fallback (forced with maxiter=0) and a warm-started sweep step
    converge to the same flow as plain Newton."""
    th = naca4_thickness(0.12)
    sol = LiftingTSDSolver(wing, 0.8, th, alpha_deg=1.0, **COARSE)
    rn = sol.steady()
    rp = sol.steady(maxiter=0)
    assert rn.converged and rp.converged
    np.testing.assert_allclose(rp.phi_lower, rn.phi_lower, atol=1e-8)
    np.testing.assert_allclose(rp.cl, rn.cl, rtol=1e-7)
    sol.set_alpha(0.0)
    r0 = sol.steady()
    sol.set_alpha(1.0)
    rw = sol.steady(warm_start=True)
    assert rw.converged and rw.iterations < rn.iterations
    np.testing.assert_allclose(rw.cl, rn.cl, rtol=1e-7)
    assert abs(r0.CL) < 1e-10


def test_camber_gives_lift_at_zero_incidence(wing):
    t, c = naca4_section("2412")
    r = LiftingTSDSolver(wing, 0.6, t, camber=c, **COARSE).steady()
    assert r.converged and r.CL > 0.1
    # thin-airfoil zero-lift angle of NACA 2412 is about -2.1 deg
    r2 = LiftingTSDSolver(wing, 0.6, t, camber=c, alpha_deg=-2.1, **COARSE).steady()
    assert abs(r2.CL) < 0.1 * r.CL


def test_lifting_gafs_finite_and_damped(wing):
    th = naca4_thickness(0.12)
    funcs, names = rigid_mode_funcs(0.25)
    sol = LiftingTSDSolver(wing, 0.8, th, alpha_deg=1.5, **COARSE)
    tab = tsd_gaf_table(wing, 0.8, [0.0, 0.2], 0.5, funcs, names, thickness=th, solver=sol)
    assert np.all(np.isfinite(tab.Q)) and np.abs(tab.Q).max() > 0
    assert np.allclose(tab.Q[0].imag, 0, atol=1e-10 * np.abs(tab.Q[0]).max())
    assert tab.Q[1, 0, 0].imag < 0
    # differs from the zero-incidence mean flow
    t0 = tsd_gaf_table(wing, 0.8, [0.0, 0.2], 0.5, funcs, names, thickness=th, **COARSE)
    assert np.abs(tab.Q - t0.Q).max() > 1e-3 * np.abs(t0.Q).max()


def test_launcher_alpha_sweep_and_checks(tmp_path):
    from flaap.io.results_txt import read_result
    from flaap.launcher import parse_text, resolve, run_config
    cfg = tmp_path / "config_sweep.txt"
    cfg.write_text("model tsd\nmode steady\nsurface1 wing, 0, 0, 0, 1, 0, 4, 0, 1, 8, 10\nsymmetry symmetric\n"
                   "mach 0.8\ntsd_thickness naca2412\ntsd_camber naca2412\nalphas -1, 1, 3\n"
                   "tsd_nx 12\ntsd_ny 6\nplots false\n", encoding="utf-8")
    res = run_config(cfg, output_dir=tmp_path / "o", verbose=False)
    assert res["status"] == "ok"
    s = res["summary"]
    assert s["alpha_deg"] == [-1.0, 0.0, 1.0] and s["converged"]
    assert s["CL"][0] < s["CL"][1] < s["CL"][2] and s["CL"][1] > 0          # cambered: lift at alpha 0
    d = read_result(tmp_path / "o/tsd_steady.txt")
    assert d["alpha_deg"] == 1.0 and len(d["cp_lower"]) == 6 and len(d["cl"]) == 6
    assert (tmp_path / "o/tsd_steady_a01.txt").exists() and (tmp_path / "o/curves_tsd_alpha_sweep.txt").exists()
    # preflight: bad camber name is an error, large incidence a warning
    _, rep = resolve(parse_text("model tsd\nmode steady\nstructure agard445\nmach 0.9\ntsd_camber naca2x\n"))
    assert rep.has_errors()
    _, rep = resolve(parse_text("model tsd\nmode steady\nstructure agard445\nmach 0.9\ntsd_alpha 8\n"))
    assert not rep.has_errors() and any("CFG-RANGE-006" in d.code for d in rep.diagnostics)
