"""Steady VLM / nonlinear UVLM / linearised UVLM against classical results.

References: lifting-line span efficiency of an elliptic wing (e = 1),
thin-aerofoil 2 pi, Wagner (Jones approximation), Theodorsen, Garrick's
plunge thrust.  Coarse meshes keep the suite fast; the converged numbers
are in ``validation/uvlm_validation.py``.
"""
from pathlib import Path

import numpy as np
import pytest
from scipy.special import hankel2

from flaap.aero.gaf import rigid_aero_modes
from flaap.aero.uvlm import UVLM, Kinematics, LinearUVLM, harmonic, steady_vlm
from flaap.geometry import AeroMesh, Surface, build_mesh

ROOT = Path(__file__).resolve().parents[1]
B = 0.5                                     # semi-chord of the 2-D (A = 1000) wing


def _strip2d(nc):
    return build_mesh([Surface("w", (0, 0, 0), 2 * B, (0, 500, 0), 2 * B, nchord=nc, nspan=1)], "symmetric")


def _elliptic(A, nc, ns, c0=1.0):
    bh = A * np.pi * c0 / 8
    y = bh * np.sin(np.linspace(0, np.pi / 2, ns + 1))
    c = c0 * np.sqrt(np.clip(1 - (y / bh) ** 2, 0, None))
    c[-1] = max(c[-1], 1e-4 * c0)
    xle = 0.25 * c0 - 0.25 * c
    xi = np.linspace(0, 1, nc + 1)
    corners, strip, ci = [], [], []
    for j in range(ns):
        for i in range(nc):
            def p(jj, ii):
                return [xle[jj] + xi[ii] * c[jj], y[jj], 0.0]
            corners.append([p(j, i), p(j + 1, i), p(j + 1, i + 1), p(j, i + 1)])
            strip.append(j)
            ci.append(i)
    n = len(corners)
    return AeroMesh(np.array(corners), np.zeros(n, int), np.array(strip), np.array(ci),
                    [Surface("w", (0, 0, 0), c0, (0, bh, 0), 0.0, nc, ns)], "symmetric")


def _theodorsen_C(k):
    return hankel2(1, k) / (hankel2(1, k) + 1j * hankel2(0, k))


def _theodorsen_cl(k, a, alpha_hat, h_hat, U=1.0):
    """Complex CL amplitude (h positive up, pitch about x = b a nose up)."""
    w = k * U / B
    L = np.pi * B ** 2 * (w ** 2 * h_hat + U * 1j * w * alpha_hat + B * a * w ** 2 * alpha_hat) \
        + 2 * np.pi * U * B * _theodorsen_C(k) * (-1j * w * h_hat + U * alpha_hat + B * (0.5 - a) * 1j * w * alpha_hat)
    return L / (0.5 * U ** 2 * 2 * B)


def _fit(t, y, w):
    sel = t >= t[-1] - 2 * np.pi / w + 1e-9
    A = np.column_stack([np.ones(sel.sum()), np.cos(w * t[sel]), np.sin(w * t[sel])])
    c = np.linalg.lstsq(A, y[sel], rcond=None)[0]
    return c[0], c[1] - 1j * c[2]


# ---------------------------------------------------------------- steady VLM
def test_vlm_elliptic_wing():
    r = steady_vlm(_elliptic(8.0, 6, 16), [2.0, 6.0])
    assert r.aspect_ratio == pytest.approx(8.0, rel=0.01)       # cosine-spaced chords
    assert r.e == pytest.approx(1.0, abs=0.01)                    # elliptic loading
    assert r.CL_alpha == pytest.approx(4.906, rel=0.03)           # Helmbold
    assert r.CDi[1] == pytest.approx(r.CL[1] ** 2 / (np.pi * 8.0), rel=0.02)


def test_vlm_2d_limit_and_rectangle():
    r = steady_vlm(_strip2d(4), [3.0])
    assert r.CL_alpha == pytest.approx(2 * np.pi, rel=2e-3)
    m = build_mesh([Surface("w", (0, 0, 0), 1.0, (0, 2.0, 0), 1.0, nchord=6, nspan=20)], "symmetric")
    r = steady_vlm(m, [4.0], x_ref=0.25)
    assert 3.5 < r.CL_alpha < 3.8 and 0.97 < r.e < 1.0           # A = 4 rectangle
    assert abs(r.x_ac - 0.25) < 0.05                              # a.c. near c/4


# ---------------------------------------------------------------- nonlinear UVLM
def test_uvlm_wagner_indicial_lift():
    nc, U, alpha = 8, 1.0, np.radians(2.0)
    r = UVLM(_strip2d(nc)).simulate(U, 6.0, Kinematics(alpha0=alpha), rho=1.0)
    phi = r.CL / (2 * np.pi * alpha)
    s = r.s
    jones = 1 - 0.165 * np.exp(-0.0455 * s) - 0.335 * np.exp(-0.3 * s)
    sel = s >= 2.0
    assert np.max(np.abs(phi[sel] - jones[sel])) < 0.03
    assert phi[-1] == pytest.approx(jones[-1], abs=0.015)


def test_uvlm_harmonic_plunge_lift_and_garrick_thrust():
    nc, U, k, h0 = 8, 1.0, 0.5, 0.05 * B
    w = k * U / B
    r = UVLM(_strip2d(nc)).simulate(U, 3 * 2 * np.pi / w, Kinematics(x_ref=0.5 * B, plunge=harmonic(h0, w)), rho=1.0)
    _, cl = _fit(r.t, r.CL, w)
    ref = _theodorsen_cl(k, -0.5, 0.0, -1j * h0)                  # h = h0 sin(wt) = Re(-i h0 e^{iwt})
    assert abs(cl) == pytest.approx(abs(ref), rel=0.09)           # 1st order in nc: +7 % at nc = 8
    assert np.degrees(np.angle(cl / ref)) == pytest.approx(0.0, abs=4.0)
    cd0, _ = _fit(r.t, r.CD, w)
    C = _theodorsen_C(k)
    ct = np.pi * k ** 2 * (h0 / B) ** 2 * abs(C) ** 2              # Garrick
    assert -cd0 == pytest.approx(ct, rel=0.05)


# ---------------------------------------------------------------- linearised UVLM
def test_linear_uvlm_gaf_vs_theodorsen():
    m = _strip2d(8)
    md = rigid_aero_modes(m, 0.25, ("plunge", "pitch"), bref=B)
    kv = np.array([0.1, 0.5])
    tab = LinearUVLM(m, md, B, wake_chords=60).gaf_table(kv, gust=False)
    S = m.total_area
    for i, k in enumerate(kv):
        for j, (ah, hh) in enumerate(((0.0, 1.0), (1.0, 0.0))):  # plunge, pitch columns
            ref = _theodorsen_cl(k, -0.5, ah, hh)
            assert abs(tab.Q[i, 0, j] / S - ref) / abs(ref) < 0.07


# ---------------------------------------------------------------- launcher
def test_uvlm_preflight_rules():
    from flaap.launcher import parse_text, resolve
    from flaap.launcher.diagnostics import Severity

    def codes(text, sev=Severity.ERROR):
        _, rep = resolve(parse_text(text))
        return {d.code for d in rep.diagnostics if d.severity is sev}

    assert not codes("model uvlm\nmode steady\n")
    assert "CFG-UVLM-001" in codes("model uvlm\nmode harmonic\nvelocity 10\nreduced_frequency 0.2\n")
    assert "CFG-UVLM-002" in codes("model uvlm\nmode gust\nvelocity 10\ngust_velocity 1\n")
    assert not codes("model uvlm\nmode gust\nvelocity 10\ngust_velocity 1\ngust_shape sharp\n")
    assert codes("model uvlm\nmode impulsive\n")                  # velocity / alpha required
    assert "CFG-RANGE-003" in codes("model flutter\nmach 0.8\naero_method machbox\n")
    assert "CFG-RANGE-004" in codes("model flutter\nmach 1.2\naero_method piston\n", Severity.WARNING)
    assert not codes("model flutter\nmach 1.2\naero_method machbox\n")
    assert "CFG-RANGE-005" in codes("model flutter\nmach 0.5\naero_method uvlm\n", Severity.WARNING)


def test_run_uvlm_steady_and_pitch(tmp_path):
    from flaap.io.results_txt import read_result
    from flaap.launcher import run_config
    res = run_config(ROOT / "configs/uvlm/ex1/config_rect_wing_steady.txt", output_dir=tmp_path / "s",
                     verbose=False, overrides={"plots": False})
    assert res["status"] == "ok"
    sm = res["summary"]
    assert sm["CL_alpha_per_rad"] == pytest.approx(4.644, rel=0.01) and 0.95 < sm["span_efficiency"] < 0.98
    st = read_result(tmp_path / "s/vlm_steady.txt")
    assert len(st["CL"]) == 8

    res = run_config(ROOT / "configs/uvlm/ex4/config_pitch_theodorsen.txt", output_dir=tmp_path / "p",
                     verbose=False, overrides={"plots": False, "verbose": False})
    assert res["status"] == "ok"
    ref = _theodorsen_cl(0.5, -0.5, -1j * np.radians(1.0), 0.0)   # pitch = th0 sin(wt)
    assert res["summary"]["CL_1st_harmonic_amplitude"] == pytest.approx(abs(ref), rel=0.06)
    th = read_result(tmp_path / "p/timehist_uvlm.txt")
    assert {"CL", "CD", "Cm"} <= set(th["series"])
