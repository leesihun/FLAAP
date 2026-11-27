"""Tests for the MITC4 plate element and the AGARD 445.6 model."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flaap.models import agard445 as A  # noqa: E402
from flaap.structures.modal import rigid_body_shapes  # noqa: E402
from flaap.structures.plate import (OrthotropicMaterial, build_plate_model,  # noqa: E402
                                    plate_mesh_trapezoid)

E, NU, RHO = 70e9, 0.3, 2700.0
LEISSA = np.array([3.492, 8.525, 21.43])   # CFFF square plate, nu = 0.3


def _square(n, h, a=1.0, **kw):
    nodes, quads = plate_mesh_trapezoid((0, 0), a, (0, a), a, n, n)
    return build_plate_model(nodes, quads, h, OrthotropicMaterial.isotropic(E, NU, RHO), **kw)


def _lam(fe, h, a=1.0, n=3):
    D = E * h ** 3 / (12 * (1 - NU ** 2))
    return fe.modes(n).omega * a * a * np.sqrt(RHO * h / D)


def test_leissa_cantilever_square():
    lam = _lam(_square(12, 0.01), 0.01)
    np.testing.assert_allclose(lam, LEISSA, rtol=0.015)


def test_thin_plate_no_shear_locking():
    # coarse mesh, very thin plate: MITC4 must not lock
    lam_thin = _lam(_square(6, 1e-3), 1e-3)
    lam_mod = _lam(_square(6, 1e-2), 1e-2)
    np.testing.assert_allclose(lam_thin, lam_mod, rtol=5e-3)
    np.testing.assert_allclose(lam_thin[:2], LEISSA[:2], rtol=0.01)


def test_mesh_convergence():
    err = [abs(_lam(_square(n, 0.01), 0.01)[2] / 21.29 - 1) for n in (4, 8, 16)]
    assert err[0] > err[1] > err[2] and err[2] < 0.01


def test_rigid_body_modes_sign_convention():
    """Free plate: plunge, nose-up pitch (uz = -dx*ry) and roll (uz = dy*rx)
    are strain free -> rotation DOFs follow rx = dw/dy, ry = -dw/dx."""
    fe = _square(5, 0.02, a=2.0, clamp_y=None, membrane=True)
    rb = rigid_body_shapes(fe.nodes, np.array([0.3, 0.7, 0.0]),
                           ("plunge", "pitch", "roll", "surge"))
    U = rb.reshape(-1, rb.shape[-1])
    U[fe.fixed] = 0.0          # rz is constrained (and unused by rigid motions)
    kmax = abs(fe.K).max()
    assert np.abs(fe.K @ U).max() < 1e-9 * kmax
    # wrong-sign pitch must have strain energy
    bad = U[:, 1].copy()
    bad[4::6] *= -1
    assert np.abs(fe.K @ bad).max() > 1e-3 * kmax


def test_kirchhoff_rotations_match_slopes():
    nc = ns = 16
    a = 1.0
    nodes, quads = plate_mesh_trapezoid((0, 0), a, (0, a), a, nc, ns)
    fe = build_plate_model(nodes, quads, 0.005, OrthotropicMaterial.isotropic(E, NU, RHO))
    sh = fe.modes(2).shapes                      # (Nn,6,nm)
    for m in range(2):
        w = sh[:, 2, m].reshape(ns + 1, nc + 1)  # [span j, chord i]
        rx = sh[:, 3, m].reshape(ns + 1, nc + 1)
        ry = sh[:, 4, m].reshape(ns + 1, nc + 1)
        dwdy = np.gradient(w, a / ns, axis=0)
        dwdx = np.gradient(w, a / nc, axis=1)
        s = np.abs(sh[:, 3:5, m]).max()
        assert np.abs(rx[2:-2, 2:-2] - dwdy[2:-2, 2:-2]).max() < 0.05 * s
        assert np.abs(ry[2:-2, 2:-2] + dwdx[2:-2, 2:-2]).max() < 0.05 * s


def test_orthotropic_rotation_and_matrices():
    m0 = OrthotropicMaterial(10e9, 1e9, 0.5e9, 0.3, 500.0, angle=0.0)
    m90 = OrthotropicMaterial(10e9, 1e9, 0.5e9, 0.3, 500.0, angle=90.0)
    Q0, Q90 = m0.Q(), m90.Q()
    assert np.isclose(Q0[0, 0], Q90[1, 1]) and np.isclose(Q0[1, 1], Q90[0, 0])
    assert abs(Q90[0, 2]) < 1e-6 * Q0[0, 0]
    iso = OrthotropicMaterial.isotropic(E, NU, RHO)
    iso_r = OrthotropicMaterial.isotropic(E, NU, RHO, angle=33.0)
    np.testing.assert_allclose(iso_r.Q(), iso.Q(), atol=1e-6 * E)
    # tapered plate: symmetric matrices and exact translational mass
    nodes, quads = plate_mesh_trapezoid((0, 0), 1.0, (0.3, 2.0), 0.5, 4, 6)
    t = 0.01 + 0.005 * nodes[:, 1]
    fe = build_plate_model(nodes, quads, t, m0, clamp_y=None)
    assert abs(fe.K - fe.K.T).max() < 1e-12 * abs(fe.K).max()
    u = np.zeros(fe.ndof)
    u[2::6] = 1.0
    # exact: rho * int t dA, with t linear in y and chord c(y) = 1 - 0.25 y
    exact = 500.0 * (0.01 * 2 + 0.0025 * 2 ** 2 / 2 - 0.00125 * 2 ** 3 / 3)
    assert np.isclose(u @ (fe.M @ u), exact, rtol=1e-10)


# ------------------------------------------------------------------ AGARD 445.6
@pytest.fixture(scope="module")
def agard_modes():
    return A.agard445_modal(4)


def test_agard_geometry():
    assert np.isclose(A.SWEEP_LE_DEG, 46.73, atol=0.01)
    assert np.isclose(A.SWEEP_MID_DEG, 43.16, atol=0.01)
    s = A.agard445_surface(6, 8)
    assert np.isclose(s.span, 0.762) and np.isclose(s.area, 0.35329, rtol=1e-4)
    # quarter-chord line swept 45 deg
    assert np.isclose(s.tip_le[0] + 0.25 * s.tip_chord - 0.25 * s.root_chord, s.span)
    assert np.isclose(A.naca65a004_thickness(0.4), 0.04, rtol=0.01)


def test_agard_frequencies(agard_modes):
    f = agard_modes.freq_hz
    np.testing.assert_allclose(f, A.EXPERIMENT_FREQ_HZ, rtol=0.06)   # measured GVT
    np.testing.assert_allclose(f, A.YATES_FE_FREQ_HZ, rtol=0.07)
    assert abs(f[0] / 9.60 - 1) < 0.02


def test_agard_mode_types(agard_modes):
    nodes = agard_modes.nodes
    tip = np.where(np.isclose(nodes[:, 1], A.SEMISPAN))[0]
    le, te = tip[np.argmin(nodes[tip, 0])], tip[np.argmax(nodes[tip, 0])]
    w = agard_modes.uz()
    assert w[le, 0] * w[te, 0] > 0          # 1st bending
    assert w[le, 1] * w[te, 1] < 0          # 1st torsion
    np.testing.assert_allclose(agard_modes.gen_mass, 1.0)


def test_agard_flutter_data_consistency():
    vol = A.frustum_volume()
    for M, d in A.FLUTTER_DATA.items():
        assert np.isclose(d["rho"] * d["mu"] * vol, 1.867, rtol=1e-3)
        fsi = A.flutter_speed_index(d["V"], d["mu"])
        assert abs(fsi / d["FSI"] - 1) < 0.002          # Yates Table II is self-consistent
        assert abs(d["V_from_FSI"] / d["V"] - 1) < 0.002
