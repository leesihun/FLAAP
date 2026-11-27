"""Tests of flaap.nonlinear (freeplay / LCO) and flaap.ase (aeroservoelasticity).

All cases use the Vasconcellos/Conner pitch-plunge-flap typical section
(exact Theodorsen RFA), so the module runs in a few seconds.
"""
import numpy as np
import pytest

from flaap.ase import (ASEModel, Actuator, Sensor, connect, feedback, gain, lead_compensator,
                       lowpass, notch, pid, section_sensor, series, sumblk, tf)
from flaap.models.typical_section import section_system, vasconcellos2016
from flaap.nonlinear import (Bilinear, Cubic, Freeplay, HingeNonlinearity, NonlinearAeroelastic,
                             PiecewiseLinear, bifurcation_sweep, describing_function_lco,
                             freeplay_keq, preloaded_freeplay)

D0 = np.radians(1.0)


@pytest.fixture(scope="module")
def section():
    sec = vasconcellos2016()
    hinge_sys, _ = section_system(sec)
    ctrl_sys, _ = section_system(sec, flap="control")
    model = NonlinearAeroelastic(hinge_sys, "beta", Freeplay(sec.K_beta, D0))
    df = describing_function_lco(model, sec.rho, np.linspace(1.0, 10.0, 37))
    return {"sec": sec, "hinge": hinge_sys, "ctrl": ctrl_sys, "model": model, "df": df}


# ----------------------------------------------------------------------
# nonlinear: describing functions
# ----------------------------------------------------------------------
@pytest.mark.parametrize("nl", [Freeplay(2.0, 0.01), Bilinear(2.0, 0.01, 0.5),
                                Freeplay(2.0, 0.01, offset=0.004), Cubic(2.0, 5e3),
                                PiecewiseLinear([-0.02, 0.0, 0.01], [3.0, 1.0, 0.2, 2.0])])
def test_describing_function_analytic_vs_quadrature(nl):
    A = np.array([0.0, 0.003, 0.01, 0.0137, 0.05, 0.3])
    for B in (0.0, 0.004, -0.012):
        N1, N0 = nl.describing_function(A, B)
        q1, q0 = HingeNonlinearity.describing_function(nl, A, B, n=20000)
        assert np.allclose(N1[1:], q1[1:], rtol=1e-6, atol=1e-6 * nl.K)
        assert np.allclose(N0, q0, rtol=1e-6, atol=1e-7 * nl.K)


def test_freeplay_keq_closed_form():
    fp = Freeplay(3.0, 0.02)
    r = np.array([0.5, 1.0, 1.5, 3.0, 50.0])
    assert np.allclose(fp.keq(r * 0.02) / 3.0, freeplay_keq(r), atol=1e-12)
    assert freeplay_keq(50.0) > 0.97 and freeplay_keq(1.0) == 0.0
    pf = preloaded_freeplay(3.0, 0.02, 0.05)
    assert not pf.is_odd and fp.is_odd
    assert pf.moment(0.0) == pytest.approx(-3.0 * 0.03)


# ----------------------------------------------------------------------
# nonlinear: linear limits, time marching, LCO
# ----------------------------------------------------------------------
def test_flutter_limits(section):
    df, sec = section["df"], section["sec"]
    lin = section["hinge"].flutter(np.linspace(1, 15, 80), sec.rho).critical()
    assert df.V_fK == pytest.approx(lin.velocity, rel=5e-3)          # K = K_beta
    assert df.V_fK == pytest.approx(9.2, rel=0.02)                    # published 9.2 m/s
    assert 2.5 < df.V_f0 < 3.5                                        # zero hinge stiffness
    st = np.array([df.V[i] for b in df.branches() for i in b if df.stable[i]])
    assert st.min() < df.V_f0 and st.max() < df.V_fK + 0.3           # LCO band (with subcritical fold)


def test_marching_methods_agree(section):
    m, sec = section["model"], section["sec"]
    x0 = {"delta": 3 * D0, "alpha": 0.01}
    ra = m.simulate(6.0, sec.rho, 4.0, x0=x0, dt=0.005, method="exact")
    rb = m.simulate(6.0, sec.rho, 4.0, x0=x0, dt=0.005, method="events")
    rc = m.simulate(6.0, sec.rho, 4.0, x0=x0, dt=0.005, method="smooth")
    scale = np.max(np.abs(ra.delta))
    assert ra.nswitch > 5
    assert np.max(np.abs(ra.delta - rb.delta)) < 1e-5 * scale
    assert np.max(np.abs(ra.delta - rc.delta)) < 0.05 * scale


@pytest.mark.parametrize("V", [3.5, 5.0, 7.0, 8.5])
def test_time_marching_vs_describing_function(section, V):
    m, sec, df = section["model"], section["sec"], section["df"]
    r = m.simulate(V, sec.rho, 60.0, x0={"delta": 3 * D0}, dt=0.005)
    info = r.lco()
    assert info.kind == "lco"
    A_df = df.amplitude(V)
    w_df = df.frequency(V)
    assert A_df.size == 1
    assert info.amplitude == pytest.approx(A_df[0], rel=0.10)
    assert info.omega == pytest.approx(w_df[0], rel=0.03)


def test_amplitude_scales_with_freeplay(section):
    sec, sysm = section["sec"], section["hinge"]
    m1 = section["model"]
    m2 = NonlinearAeroelastic(sysm, "beta", Freeplay(sec.K_beta, 2 * D0))
    r1 = m1.simulate(5.0, sec.rho, 40.0, x0={"delta": 3 * D0}, dt=0.005)
    r2 = m2.simulate(5.0, sec.rho, 40.0, x0={"delta": 6 * D0}, dt=0.005)
    assert r2.lco().amplitude / r1.lco().amplitude == pytest.approx(2.0, rel=1e-3)
    assert r2.lco("alpha").amplitude / r1.lco("alpha").amplitude == pytest.approx(2.0, rel=1e-3)


def test_no_lco_below_fold(section):
    m, sec, df = section["model"], section["sec"], section["df"]
    st = [df.V[i] for b in df.branches() for i in b if df.stable[i]]
    V = 0.8 * min(st)                                 # below the fold (and V_f0)
    assert V < df.V_f0
    for a0 in (1.5, 4.0, 10.0):
        r = m.simulate(V, sec.rho, 40.0, x0={"delta": a0 * D0, "alpha": 0.02}, dt=0.005)
        assert r.lco().kind == "decay"


def test_bifurcation_sweep(section):
    m, sec = section["model"], section["sec"]
    out = bifurcation_sweep(m, [1.5, 4.0, 7.0], sec.rho, 40.0, x0={"delta": 3 * D0}, dt=0.005)
    assert list(out["kind"]) == ["decay", "lco", "lco"]
    assert out["amp"][2] > out["amp"][1] > D0


# ----------------------------------------------------------------------
# ASE
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def ase(section):
    sec = section["sec"]
    act = Actuator.second_order("beta", 2 * np.pi * 15.0, 0.7)
    sens = [Sensor("pitch", [0.0, 1.0], "disp"), Sensor("pitch_rate", [0.0, 1.0], "rate"),
            Sensor("plunge", [1.0, 0.0], "disp"), Sensor("plunge_rate", [1.0, 0.0], "rate"),
            section_sensor(sec, 0.5 * sec.b, "acc", "acc_te")]
    return ASEModel(section["ctrl"], [act], sens, sec.rho)


LAW = ["pitch", "pitch_rate", "plunge", "plunge_rate"]


def test_tf_blocks():
    s = 1j * np.array([0.3, 2.0, 17.0])
    T = tf([2.0, 3.0, 1.0], [1.0, 0.5, 4.0])
    H = T.freq_response(s.imag)[:, 0, 0]
    assert np.allclose(H, np.polyval([2, 3, 1], s) / np.polyval([1, 0.5, 4], s))
    P = pid(2.0, 0.5, 0.1, tau_d=0.01)
    Hp = P.freq_response(s.imag)[:, 0, 0]
    assert np.allclose(Hp, 2.0 + 0.5 / s + 0.1 * s / (0.01 * s + 1))
    L = lead_compensator(40.0, 10.0, 2.0)
    HL = L.freq_response([10.0])[0, 0, 0]
    assert abs(HL) == pytest.approx(2.0) and np.degrees(np.angle(HL)) == pytest.approx(40.0)
    N = notch(10.0)
    assert abs(N.freq_response([10.0])[0, 0, 0]) == pytest.approx(0.04, rel=1e-6)
    # series / sum / feedback with an algebraic loop: y = G/(1 - G k) for u -> y
    G = series(lowpass(5.0, 1, input="u", output="v"), gain(3.0, "v", "y"))
    cl = connect([G, sumblk("u", ["r", "fb"]), gain(0.5, "y", "fb")], ["r"], ["y"])
    Hcl = cl.freq_response([0.0, 4.0])[:, 0, 0]
    Hg = 3 * 5 / (1j * np.array([0, 4.0]) + 5)
    assert np.allclose(Hcl, Hg / (1 - 0.5 * Hg))
    fb = feedback(gain(2.0, "a", "b"), gain(0.25, "b", "a"), inputs=[], outputs=["b"])
    assert fb.D.size == 0


@pytest.mark.parametrize("order", [2, 3])
def test_actuator_derivatives(order):
    act = (Actuator.second_order("d", 60.0, 0.6) if order == 2
           else Actuator.third_order("d", 60.0, 0.6, 0.01))
    w = np.array([0.5, 7.0, 45.0, 300.0])
    H = act.state_space().freq_response(w)
    assert np.allclose(H[:, 0, 0], act.response(w))
    assert np.allclose(H[:, 1, 0], 1j * w * H[:, 0, 0])
    assert np.allclose(H[:, 2, 0], (1j * w) ** 2 * H[:, 0, 0])
    # time domain: delta'' from the actuator equals the derivative of delta'
    t = np.linspace(0, 0.3, 3001)
    u = np.sin(40 * t)[:, None]
    y, _ = act.state_space().simulate(t, u)
    assert np.allclose(np.gradient(y[:, 1], t)[5:-5], y[5:-5, 2], atol=2e-3 * np.abs(y[:, 2]).max())


def test_zero_gain_reproduces_open_loop(section, ase):
    sec = section["sec"]
    V = np.linspace(2.0, 14.0, 61)
    ref = section["ctrl"].flutter(V, sec.rho).critical().velocity
    for K in (None, gain(np.zeros((1, 4)), LAW, "beta_cmd")):
        r = ase.flutter(V, K)
        assert r.critical().velocity == pytest.approx(ref, rel=1e-6)
    assert ase.flutter_speed((2.0, 14.0), gain(np.zeros((1, 4)), LAW, "beta_cmd")) == pytest.approx(ref, rel=2e-3)
    # the structural part of the closed loop with K = 0 is the open-loop plant
    lam_ol = np.linalg.eigvals(section["ctrl"].state_space(8.0, sec.rho, gust=False).A)
    lam_cl = np.linalg.eigvals(ase.closed_loop(8.0, gain(np.zeros((1, 4)), LAW, "beta_cmd"), gust=False).A)
    assert all(np.min(np.abs(lam_cl - l)) < 1e-8 * max(1, abs(l)) for l in lam_ol)


def test_flutter_suppression_and_margins(section, ase):
    K = gain(np.array([[1.6, 0.16, -18.0, -2.3]]), LAW, "beta_cmd")
    vf0 = ase.flutter_speed((2.0, 30.0), None)
    vf1 = ase.flutter_speed((2.0, 30.0), K)
    assert vf0 == pytest.approx(9.24, rel=0.01)
    assert vf1 > 1.6 * vf0
    for V in (9.0, 12.0, 15.0):
        m = ase.margins(V, K)
        assert m.closed_loop_stable and m.nyquist_stable
        # gain margins: scaling the loop by the margin factor puts a root on the axis
        for f in (m.gm_lower, m.gm_upper):
            if 0 < f < np.inf:
                lo = ase.closed_loop(V, K, 0.97 * f, gust=False, extra_outputs=[]).A
                hi = ase.closed_loop(V, K, 1.03 * f, gust=False, extra_outputs=[]).A
                s = np.max(np.linalg.eigvals(lo).real) * np.max(np.linalg.eigvals(hi).real)
                assert s < 0
    m = ase.margins(12.0, K)
    assert m.P == 2                                      # open loop is beyond flutter
    m_bad = ase.margins(19.0, K)
    assert not m_bad.closed_loop_stable and not m_bad.nyquist_stable


def test_acceleration_feedback_algebraic_loop(section, ase):
    """An accelerometer sees delta'' (feed-through from the actuator command
    for a 2nd-order actuator), so a static gain creates an algebraic loop."""
    K = gain(0.01, "acc_te", "beta_cmd")
    P = ase.plant(6.0, gust=False, extra_outputs=[])
    ia, ic = P.outputs.index("acc_te"), P.inputs.index("beta_cmd")
    assert abs(P.D[ia, ic]) > 0
    cl = ase.closed_loop(6.0, K, gust=False, extra_outputs=[])
    # manual closure: u = k y = k (C x + D u)  ->  u = k C x / (1 - k D)
    k = 0.01
    A_ref = P.A + P.B[:, [ic]] @ (k * P.C[[ia]] / (1 - k * P.D[ia, ic]))
    assert np.allclose(np.sort_complex(np.linalg.eigvals(cl.A)), np.sort_complex(np.linalg.eigvals(A_ref)))
