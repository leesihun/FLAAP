"""Tests for the frequency-domain transonic small-disturbance solver
(:mod:`flaap.aero.tsd`)."""
import numpy as np
import pytest
import scipy.sparse as sp

from flaap.aero.gaf import AeroModel, GAFTable, rigid_aero_modes
from flaap.aero.tsd import (TSDSolver, _block_solve, cp_critical, modes_from_structure,
                            naca4_thickness, naca65a004, rigid_mode_funcs, tsd_gaf_table)
from flaap.geometry import Surface, build_mesh

COARSE = dict(nchord=16, nspan=8, nz=18, n_up=10, n_down=10, n_out=8)


def test_block_solver_matches_sparse_direct():
    rng = np.random.default_rng(0)
    P, nb = 7, 6
    N = P * nb
    K = sp.lil_matrix((N, N), dtype=complex)
    for c in range(nb):
        for d in (-2, -1, 0, 1):
            if 0 <= c + d < nb:
                blk = rng.normal(size=(P, P)) + 1j * rng.normal(size=(P, P))
                if d == 0:
                    blk += 10 * np.eye(P)
                K[c * P:(c + 1) * P, (c + d) * P:(c + d + 1) * P] = blk
    K = K.tocsr()
    R = rng.normal(size=(N, 3))
    x = _block_solve(K, R, P, nb)
    np.testing.assert_allclose(K @ x, R, atol=1e-10)
    # check-pointed variant (memory budget exceeded)
    np.testing.assert_allclose(_block_solve(K, R, P, nb, max_bytes=0), x, atol=1e-12)
    x1 = _block_solve(K, R[:, 0], P, nb)
    assert x1.shape == (N,)
    np.testing.assert_allclose(x1, x[:, 0], atol=1e-12)


def test_flat_plate_matches_dlm_subsonic():
    """(a) zero thickness, M = 0.5, rectangular AR 4: rigid plunge/pitch GAFs
    within 10 % of the doublet-lattice method at k = 0 and 0.3."""
    c, s, M = 1.0, 2.0, 0.5
    surf = Surface("w", (0, 0, 0), c, (0, s, 0), c, nchord=10, nspan=16)
    bref, xref = 0.5 * c, 0.25 * c
    mesh = build_mesh([surf], "symmetric")
    kvals = [0.0, 0.3]
    dlm = AeroModel(mesh, rigid_aero_modes(mesh, xref), M, bref).table(kvals, gust=False)
    funcs, names = rigid_mode_funcs(xref)
    tab = tsd_gaf_table(surf, M, kvals, bref, funcs, names, **COARSE)
    assert isinstance(tab, GAFTable)
    assert tab.Q.shape == (2, 2, 2) and tab.cp.shape == (2, 16 * 8, 2)
    assert tab.names == names
    for i in range(len(kvals)):
        Qd, Qt = dlm.Q[i], tab.Q[i]
        scale = np.abs(Qd).max()
        # every term within 10 % of the matrix scale ...
        assert np.abs(Qt - Qd).max() < 0.10 * scale
        # ... and the lift terms individually within 10 %
        for a, b in [(0, 1)] + ([(0, 0)] if kvals[i] > 0 else []):
            assert abs(Qt[a, b] - Qd[a, b]) < 0.10 * abs(Qd[a, b]), (kvals[i], a, b)
    # steady pitch: centre of pressure close to the DLM one (edge correction)
    xcp_dlm = xref - dlm.Q[0, 1, 1].real / dlm.Q[0, 0, 1].real
    xcp_tsd = xref - tab.Q[0, 1, 1].real / tab.Q[0, 0, 1].real
    assert abs(xcp_tsd - xcp_dlm) < 0.01 * c
    # plunge aerodynamic damping
    assert tab.Q[1, 0, 0].imag < 0


def test_steady_transonic_shock():
    """(b) NACA 0012 on a high-aspect-ratio wing at M = 0.8: supersonic pocket
    terminated by a shock."""
    M = 0.8
    surf = Surface("w", (0, 0, 0), 1.0, (0, 4.0, 0), 1.0, nchord=10, nspan=10)
    sol = TSDSolver(surf, M, naca4_thickness(0.12), nchord=32, nspan=6, nz=18,
                    n_up=10, n_down=10, n_out=8)
    st = sol.steady()
    assert st.converged and st.residuals[-1] < 1e-8
    assert st.n_supersonic > 20
    cp = st.cp[0]                                    # root strip
    assert cp.min() < st.cp_star < 0
    assert np.isclose(st.cp_star, cp_critical(M, tsd=True))
    i_min = int(np.argmin(cp))
    # shock: pressure rise across <= 3 nodes downstream of the minimum, from
    # supersonic to subsonic
    jump = cp[i_min + 3] - cp[i_min]
    assert jump > 0.15
    assert cp[i_min + 3] > st.cp_star
    assert 0.15 < st.x[0, i_min] < 0.8
    # symmetric-flow sanity: stagnation-like compression at the nose
    assert cp[0] > 0


def test_agard_transonic_gafs_finite():
    """(c) AGARD 445.6 (NACA 65A004 mean flow) at M = 0.96 with structural modes."""
    from flaap.models.agard445 import B_REF, agard445_modal, agard445_surface
    mm = agard445_modal(2, n_chord=8, n_span=8)
    funcs, names = modes_from_structure(mm)
    assert names == ["1st bending", "1st torsion"]
    surf = agard445_surface()
    tab = tsd_gaf_table(surf, 0.96, [0.0, 0.2], B_REF, funcs, names, thickness=naca65a004(),
                        **COARSE)
    sol = tab.tsd_solver
    mf = sol.mean_flow
    assert mf.converged and mf.n_supersonic > 0
    assert np.all(np.isfinite(tab.Q)) and np.all(np.isfinite(tab.cp))
    assert np.abs(tab.Q).max() > 0
    assert np.allclose(tab.Q[0].imag, 0, atol=1e-10 * np.abs(tab.Q[0]).max())
    # bending is aerodynamically damped
    assert tab.Q[1, 0, 0].imag < 0
    # spline modes reproduce the structural normal displacement at the nodes
    h, _ = funcs[0](mm.nodes[:, 0], mm.nodes[:, 1])
    np.testing.assert_allclose(h, mm.shapes[:, 2, 0], atol=1e-6 * np.abs(mm.shapes[:, 2, 0]).max())


def test_rigid_mode_conventions():
    funcs, names = rigid_mode_funcs(0.3, ("plunge", "pitch"))
    x = np.array([0.0, 0.3, 1.0])
    h, hx = funcs[1](x, 0 * x)
    np.testing.assert_allclose(h, -(x - 0.3))
    np.testing.assert_allclose(hx, -1.0)
    h, hx = funcs[0](x, 0 * x)
    np.testing.assert_allclose(h, 1.0)
    np.testing.assert_allclose(hx, 0.0)


def test_input_validation():
    surf = Surface("w", (0, 0, 0), 1.0, (0, 2.0, 0), 1.0, nchord=4, nspan=4)
    with pytest.raises(ValueError):
        TSDSolver(surf, 1.05)
    with pytest.raises(NotImplementedError):
        TSDSolver(surf, 0.5, symmetry="antisymmetric")
