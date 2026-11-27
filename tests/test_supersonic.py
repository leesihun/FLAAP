"""Supersonic aerodynamics (piston theory, Mach box) against linear theory."""
import numpy as np
import pytest
from scipy.special import ellipe

from flaap.aero.gaf import AeroModes, rigid_aero_modes
from flaap.aero.supersonic import (MachBox, garrick_rubinow_2d, piston_theory_table,
                                   supersonic_gaf_table)
from flaap.geometry import Surface, build_mesh


def _alpha_field(p):
    return np.zeros((len(p), 1)), -np.ones((len(p), 1))


def _cl_alpha(surfaces, mach, nx, symmetry="symmetric"):
    mb = MachBox(surfaces, mach, nx=nx, symmetry=symmetry)
    _, d = mb.table(_alpha_field, 1.0, [0.0], gust=False)
    a = mb.area[mb.wing]
    return float(np.real(d[0, :, 0] @ a) / a.sum())


def test_piston_ackeret():
    """Steady Van Dyke piston theory = Ackeret 2-D: CL_alpha = 4 / beta."""
    M = 2.0
    mesh = build_mesh([Surface("w", (0, 0, 0), 1.0, (0, 1.0, 0), 1.0, nchord=8, nspan=2)])
    modes = rigid_aero_modes(mesh, 0.25, ("plunge", "alpha"))
    tab = piston_theory_table(mesh, modes, M, 0.5, [0.0], gust=False)
    assert tab.Q[0, 0, 1].real / mesh.total_area == pytest.approx(4 / np.sqrt(M * M - 1), rel=1e-12)


@pytest.mark.parametrize("aspect", [2.0, 40.0])
def test_machbox_rectangle_lift_slope(aspect):
    """Rectangular wing, supersonic tips: CL_alpha = 4/beta (1 - 1/(2 beta A))."""
    M = 2.0
    b = np.sqrt(M * M - 1)
    cla = _cl_alpha([Surface("w", (0, 0, 0), 1.0, (0, aspect / 2, 0), 1.0)], M, 40)
    assert cla == pytest.approx(4 / b * (1 - 1 / (2 * b * aspect)), rel=4e-3)


def test_machbox_delta_subsonic_le():
    """Delta wing with subsonic leading edges (m = beta tan(eps) = 0.69):
    CL_alpha = 2 pi tan(eps) / E(sqrt(1 - m^2))."""
    M, tane = 2.0, 0.4
    m = np.sqrt(M * M - 1) * tane
    cla = _cl_alpha([Surface("d", (0, 0, 0), 1.0, (1.0, tane, 0), 0.0)], M, 40)
    assert cla == pytest.approx(2 * np.pi * tane / ellipe(1 - m * m), rel=6e-3)


def test_machbox_symmetric_half_equals_full_span():
    M = 1.5
    half = _cl_alpha([Surface("w", (0, 0, 0), 1.0, (0.3, 1.0, 0), 0.6)], M, 30, "symmetric")
    full = _cl_alpha([Surface("w", (0.3, -1.0, 0), 0.6, (0, 0, 0), 1.0),
                      Surface("w", (0, 0, 0), 1.0, (0.3, 1.0, 0), 0.6)], M, 30, "none")
    assert half == pytest.approx(full, rel=2e-3)


@pytest.mark.parametrize("mach, k", [(1.2, 0.5), (2.0, 1.0)])
def test_machbox_2d_unsteady_vs_garrick_rubinow(mach, k):
    """Root strip of a very wide rectangular wing vs the exact 2-D integral."""
    b, c, xr = 0.5, 1.0, 0.5
    mb = MachBox([Surface("w", (0, 0, 0), c, (0, 50 * c, 0), c)], mach, nx=40, symmetry="symmetric")
    kq = k / b
    pc, pf = mb.wing_points("center"), mb.wing_points("segment")
    W = np.stack([1j * kq * np.ones(len(pc)), -1 - 1j * kq * (pc[:, 0] - xr)], -1)
    _, d = mb.solve(k, b, W)
    sel = mb.Y[mb.wing] < mb.dy
    h = np.stack([np.ones(len(pf)), -(pf[:, 0] - xr)], -1)[sel]
    Q = np.einsum("xi,xj,x->ij", h, d[sel], mb.seg_len[mb.wing][sel])
    ref = garrick_rubinow_2d(mach, k, b, c, xr)
    assert np.abs(Q - ref).max() < 3e-3 * np.abs(ref).max()


# ------------------------------------------------------------------ panel flutter
_BL = [4.730040745, 7.853204624, 10.99560784, 14.13716549, 17.27875966, 20.42035225]


def _panel_shapes(x, kind, n):
    h, hx, K = [], [], []
    for i in range(n):
        if kind == "ss":
            a = (i + 1) * np.pi
            h.append(np.sin(a * x))
            hx.append(a * np.cos(a * x))
            K.append(a ** 4 / 2)
        else:
            b = _BL[i]
            s = (np.cosh(b) - np.cos(b)) / (np.sinh(b) - np.sin(b))
            h.append(np.cosh(b * x) - np.cos(b * x) - s * (np.sinh(b * x) - np.sin(b * x)))
            hx.append(b * (np.sinh(b * x) + np.sin(b * x) - s * (np.cosh(b * x) - np.cos(b * x))))
            K.append(b ** 4)
    return np.array(h).T, np.array(hx).T, np.diag(K), np.eye(n) * (0.5 if kind == "ss" else 1.0)


def panel_lambda_cr(kind, n=6, nc=60, mach=2.0):
    """Coalescence boundary lambda = 2 q a^3 / (beta D), one side exposed."""
    mesh = build_mesh([Surface("p", (0, 0, 0), 1.0, (0, 1.0, 0), 1.0, nchord=nc, nspan=1)])
    hc, hxc, K, Mm = _panel_shapes(mesh.colloc[:, 0], kind, n)
    hf = _panel_shapes(mesh.force_point[:, 0], kind, n)[0]
    tab = piston_theory_table(mesh, AeroModes(hc, hxc, hf), mach, 1.0, [0.0], gust=False)
    A = 0.5 * tab.Q[0].real
    beta = np.sqrt(mach * mach - 1)

    def unstable(lam):
        ev = np.linalg.eigvals(np.linalg.solve(Mm, K - 0.5 * lam * beta * A))
        return np.abs(ev.imag).max() > 1e-8 * np.abs(ev).max()

    lo, hi = 10.0, 2000.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    return 0.5 * (lo + hi)


@pytest.mark.parametrize("kind, ref", [("ss", 343.3), ("clamped", 636.6)])
def test_panel_flutter_piston(kind, ref):
    """2-D panel flutter, Hedgepeth (1957) / Dowell (1975): lambda_cr = 343.3
    (simply supported), 636.6 (clamped)."""
    assert panel_lambda_cr(kind) == pytest.approx(ref, rel=2e-3)


def test_supersonic_table_feeds_pk():
    """Dispatcher returns DLM-convention tables usable by pk_method."""
    from flaap.flutter import pk_method
    mesh = build_mesh([Surface("w", (0, 0, 0), 1.0, (0.2, 2.0, 0), 0.8, nchord=6, nspan=6)],
                      "symmetric")
    modes = rigid_aero_modes(mesh, 0.4, ("plunge", "pitch"))
    kv = np.r_[1e-4, 0.1, 0.3, 0.6]
    tp = supersonic_gaf_table(mesh, modes, 2.0, 0.5, kv, method="piston")
    tm = supersonic_gaf_table(mesh, modes, 2.0, 0.5, kv, method="machbox", nx=30)
    for t in (tp, tm):
        assert t.Q.shape == (4, 2, 2) and t.Qg.shape == (4, 2)
        assert t.cp.shape == (4, mesh.n, 2) and t.cpg.shape == (4, mesh.n)
    # steady pitch lift: Mach box within the piston (2-D) value, tip losses
    assert 0.8 < tm.Q[0, 0, 1].real / tp.Q[0, 0, 1].real < 1.0
    res = pk_method(np.diag([10.0, 0.5]), np.diag([10.0 * 20 ** 2, 0.5 * 50 ** 2]), tm, 0.5,
                    np.linspace(10, 200, 20))
    assert np.all(np.isfinite(res.freq)) and np.all(np.isfinite(res.g))
    with pytest.raises(ValueError):
        supersonic_gaf_table(mesh, modes, 0.8, 0.5, kv, method="machbox")
