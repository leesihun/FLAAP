"""Pitch-plunge-flap typical section with Theodorsen aerodynamics.

Three-degree-of-freedom airfoil (plunge ``h``, pitch ``alpha``, trailing-
edge flap ``beta``) in incompressible flow.  This is the model used by
Conner, Tang, Dowell & Virgin (1997, J. Fluids Struct. 11, 89-109) for
control-surface freeplay and its later replications.  Example: the
Vasconcellos, Abdelkefi, Hajj, Almeida & Marques (2016) rig (J. Vib.
Control 22(6), 1633-1644), whose parameters are reproduced in
:func:`vasconcellos2016`.

FLAAP conventions
-----------------
* Generalised coordinates ``xi = [h, alpha, beta]``:

  - ``h``: plunge at the elastic axis, positive **up**, in m;
  - ``alpha``: pitch, positive nose up, in rad;
  - ``beta``: flap rotation about the hinge, positive trailing edge
    down, in rad.

* Geometry in semichords ``b`` measured from mid-chord. ``a b`` is the
  elastic axis and ``c b`` the hinge line.
* Harmonic motion ``e^{i w t}``, reduced frequency ``k = w b / U``.
* Generalised aerodynamic forces ``F = q Q(k) xi`` with
  ``q = rho U^2 / 2``.  ``Q`` is per unit span, multiplied by ``span``.
  This is the same convention as :class:`flaap.aero.gaf.GAFTable`, so
  the RFA and state-space machinery of FLAAP applies unchanged.

Aerodynamics (Theodorsen 1935)
------------------------------
The equations below are written with the plunge ``w`` positive **down**
(Theodorsen's convention). ``F = [-L, M_alpha, M_beta]`` is the
generalised force in ``[w, alpha, beta]``::

    F = M_nc xi'' + B_nc xi' + K_nc xi + 2 pi rho U b C(k) R f(t)

    M_nc = -rho b^2 [[pi,      -pi a b,               -T1 b               ],
                     [-pi a b,  pi b^2 (1/8 + a^2),  -(T7 + (c-a) T1) b^2],
                     [-T1 b,   -(T7 + (c-a) T1) b^2, -T3 b^2 / pi         ]]
    B_nc = -rho b^2 U [[0, pi,                         -T4                                ],
                       [0, pi (1/2 - a) b,             (T1 - T8 - (c-a) T4 + T11/2) b      ],
                       [0, (-2 T9 - T1 + T4 (a - 1/2)) b, -T4 T11 b / (2 pi)              ]]
    K_nc = -rho b^2 U^2 [[0, 0, 0], [0, 0, T4 + T10], [0, 0, (T5 - T4 T10) / pi]]
    R    = [-1, b (a + 1/2), -b T12 / (2 pi)]
    f(t) = U alpha + w' + b (1/2 - a) alpha' + U T10 beta / pi + b T11 beta' / (2 pi)

These are Conner et al. (1997) eqs. (2)-(4) and Vasconcellos et al.
(2016) eqs. (5)-(11). For FLAAP's plunge-up coordinate the matrices
are transformed with ``T = diag(-1, 1, 1)``: ``Q_up = T Q_w T``.

Gust column (vertical gust ``w_g``, positive up, per unit ``w_g / U``,
gust reference at the leading edge): ``4 pi b S(k) e^{-ik} R``, with
``S`` the Sears function referenced to mid-chord.  The hinge-moment row
of the gust column uses the circulatory distribution ``R``.  This is an
approximation: the exact gust pressure distribution differs slightly
for the flap.

Structure, coordinates ``[h (up), alpha, beta]``::

    M = [[m_T,   -S_a,               -S_b],
         [-S_a,   I_a,  I_b + b (c - a) S_b],
         [-S_b,   I_b + b (c - a) S_b, I_b]]
    K = diag(K_h, K_a, K_b)

    S_a = m_W b x_a,   S_b = m_W b x_b,   I_a = m_W b^2 r_a^2,
    I_b = m_W b^2 r_b^2, K_h = m_W w_h^2, K_a = I_a w_a^2, K_b = I_b w_b^2.

Here ``m_W`` is the wing (pitching) mass and ``m_T`` the total plunging
mass (wing plus the moving suspension parts).
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from ..aero.gaf import GAFTable
from ..aero.strip import sears, theodorsen


# ----------------------------------------------------------------------
# Theodorsen T-functions
# ----------------------------------------------------------------------
def theodorsen_T(c: float, a: float) -> dict:
    """Theodorsen's geometric functions ``T1 .. T14`` for hinge ``c`` and
    elastic axis ``a``, both in semichords from mid-chord (NACA Rep. 496)."""
    s = np.sqrt(1.0 - c * c)
    ac = np.arccos(c)
    T = {}
    T[1] = -s * (2 + c * c) / 3 + c * ac
    T[2] = c * (1 - c * c) - s * (1 + c * c) * ac + c * ac * ac
    T[3] = -(1 / 8 + c * c) * ac * ac + c * s * ac * (7 + 2 * c * c) / 4 - (1 - c * c) * (5 * c * c + 4) / 8
    T[4] = -ac + c * s
    T[5] = -(1 - c * c) - ac * ac + 2 * c * s * ac
    T[6] = T[2]
    T[7] = -(1 / 8 + c * c) * ac + c * s * (7 + 2 * c * c) / 8
    T[8] = -s * (2 * c * c + 1) / 3 + c * ac
    T[9] = 0.5 * (s ** 3 / 3 + a * T[4])
    T[10] = s + ac
    T[11] = ac * (1 - 2 * c) + s * (2 - c)
    T[12] = s * (2 + c) - ac * (2 * c + 1)
    T[13] = 0.5 * (-T[7] - (c - a) * T[1])
    T[14] = 1 / 16 + 0.5 * a * c
    return T


def theodorsen_flap_matrices(b: float, a: float, c: float):
    """Non-dimensional Theodorsen matrices in Theodorsen's plunge-down
    coordinates ``[w, alpha, beta]``.

    Returns ``(Mh, Bh, Kh, R, S1, S2)`` with ``M_nc = -rho b^2 Mh``,
    ``B_nc = -rho b^2 U Bh``, ``K_nc = -rho b^2 U^2 Kh`` and
    ``f = U S1 . xi + S2 . xi'``.
    """
    T = theodorsen_T(c, a)
    pi = np.pi
    Mh = np.array([
        [pi, -pi * a * b, -T[1] * b],
        [-pi * a * b, pi * b * b * (1 / 8 + a * a), -(T[7] + (c - a) * T[1]) * b * b],
        [-T[1] * b, 2 * T[13] * b * b, -T[3] * b * b / pi]])
    Bh = np.array([
        [0.0, pi, -T[4]],
        [0.0, pi * (0.5 - a) * b, (T[1] - T[8] - (c - a) * T[4] + 0.5 * T[11]) * b],
        [0.0, (-2 * T[9] - T[1] + T[4] * (a - 0.5)) * b, -T[4] * T[11] * b / (2 * pi)]])
    Kh = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, T[4] + T[10]], [0.0, 0.0, (T[5] - T[4] * T[10]) / pi]])
    R = np.array([-1.0, b * (a + 0.5), -b * T[12] / (2 * pi)])
    S1 = np.array([0.0, 1.0, T[10] / pi])
    S2 = np.array([1.0, b * (0.5 - a), b * T[11] / (2 * pi)])
    return Mh, Bh, Kh, R, S1, S2


def theodorsen_flap_gaf(b: float, a: float, c: float, kvals, span: float = 1.0,
                        gust: bool = True, names=("h", "alpha", "beta")) -> GAFTable:
    """GAF table ``Q(k)`` (nk, 3, 3) of the pitch-plunge-flap section.

    Coordinates are ``[h (up), alpha, beta]`` and ``bref = b``.  Forces
    per unit ``q`` are multiplied by ``span``.  The flap enters exactly
    like a control-surface aero mode (TE down, rad), so a column subset
    of this table can drive :class:`flaap.timedomain.AeroelasticSystem`
    with the flap as a prescribed control input.
    """
    k = np.atleast_1d(np.asarray(kvals, float))
    Mh, Bh, Kh, R, S1, S2 = theodorsen_flap_matrices(b, a, c)
    T = np.diag([-1.0, 1.0, 1.0])
    Q = np.zeros((k.size, 3, 3), complex)
    Qg = np.zeros((k.size, 3), complex) if gust else None
    C = theodorsen(k)
    S = sears(k)
    for i, kk in enumerate(k):
        Qw = (2 * kk * kk * Mh - 2j * kk * b * Bh - 2 * b * b * Kh
              + 4 * np.pi * b * C[i] * np.outer(R, S1 + 1j * (kk / b) * S2))
        Q[i] = span * (T @ Qw @ T)
        if gust:
            # an upward gust acts like a downward plunge velocity: f += w_g
            Qg[i] = span * (T @ (4 * np.pi * b * S[i] * np.exp(-1j * kk) * R))
    return GAFTable(k, Q, Qg, 0.0, b, list(names))


# ----------------------------------------------------------------------
# structural model
# ----------------------------------------------------------------------
@dataclass
class TypicalSection:
    """Pitch-plunge-flap typical section (per unit span unless ``span`` is set).

    Frequencies are the *uncoupled* ones as defined by Conner et al.
    (1997) and Vasconcellos et al. (2016): ``K_h = m_W w_h^2``,
    ``K_a = I_a w_a^2`` and ``K_b = I_b w_b^2``.
    """

    b: float                 # semichord [m]
    a: float                 # elastic axis [semichords from mid-chord]
    c: float                 # hinge line [semichords from mid-chord]
    m_W: float               # wing (pitching) mass [kg] (per span)
    m_T: float               # total plunging mass [kg]
    x_alpha: float           # (c.g. - e.a.) / b, positive aft
    x_beta: float            # flap static moment about hinge / (m_W b)
    r_alpha2: float          # I_a / (m_W b^2)
    r_beta2: float           # I_b / (m_W b^2)
    omega_h: float           # [rad/s]
    omega_alpha: float
    omega_beta: float
    span: float = 1.0
    rho: float = 1.225
    zeta: tuple = (0.0, 0.0, 0.0)    # modal-free viscous damping ratios on the diagonal

    # ------------------------------------------------------------------
    @property
    def S_alpha(self):
        return self.m_W * self.b * self.x_alpha

    @property
    def S_beta(self):
        return self.m_W * self.b * self.x_beta

    @property
    def I_alpha(self):
        return self.m_W * self.b ** 2 * self.r_alpha2

    @property
    def I_beta(self):
        return self.m_W * self.b ** 2 * self.r_beta2

    @property
    def K_h(self):
        return self.m_W * self.omega_h ** 2

    @property
    def K_alpha(self):
        return self.I_alpha * self.omega_alpha ** 2

    @property
    def K_beta(self):
        return self.I_beta * self.omega_beta ** 2

    @property
    def mass_ratio(self) -> float:
        """``mu = m_W / (pi rho b^2)`` (per unit span)."""
        return self.m_W / (np.pi * self.rho * self.b ** 2)

    def with_(self, **kw) -> "TypicalSection":
        return replace(self, **kw)

    def matrices(self, K_beta: float | None = None):
        """``(M, C, K)`` (3x3) in ``[h, alpha, beta]``, multiplied by ``span``.

        ``K_beta`` overrides the nominal hinge stiffness (0 for the
        freeplay / zero-stiffness model).  Damping is ``2 zeta_i sqrt(K_ii
        M_ii)`` on the diagonal (if the stiffness is zero, the nominal
        hinge stiffness is used for the flap damping).
        """
        b, a, c = self.b, self.a, self.c
        Sa, Sb, Ia, Ib = self.S_alpha, self.S_beta, self.I_alpha, self.I_beta
        Iab = Ib + b * (c - a) * Sb
        M = np.array([[self.m_T, -Sa, -Sb],
                      [-Sa, Ia, Iab],
                      [-Sb, Iab, Ib]]) * self.span
        kb = self.K_beta if K_beta is None else K_beta
        K = np.diag([self.K_h, self.K_alpha, kb]) * self.span
        knom = np.array([self.K_h, self.K_alpha, self.K_beta]) * self.span
        C = np.diag(2 * np.asarray(self.zeta, float) * np.sqrt(knom * np.diag(M)))
        return M, C, K

    def gaf_table(self, kvals=None, gust: bool = True) -> GAFTable:
        if kvals is None:
            kvals = np.r_[0.0, np.geomspace(0.01, 4.0, 30)]
        return theodorsen_flap_gaf(self.b, self.a, self.c, kvals, self.span, gust)


CONNER1997_DAMPING = (0.0113, 0.01626, 0.0115)     # (zeta_h, zeta_alpha, zeta_beta)


def vasconcellos2016(**kw) -> TypicalSection:
    """Typical section of Vasconcellos et al. (2016), Table 1.

    This is a 3-DOF pitch-plunge-flap airfoil with control-surface
    freeplay, a replication of the Conner/Tang/Dowell (1997) experiment.
    ``rho = 1.0844 kg/m^3``. The reported linear flutter speed
    (Theodorsen/Wagner, no freeplay) is ``U_f = 9.2 m/s``, and the
    experimental freeplay is +-3.57 deg.  The paper's aerodynamic
    matrices treat ``m_W`` as mass per unit span (``mu = m_W / (pi rho
    b^2) = 28.2``), so ``span = 1``.

    The structural damping ratios default to the values measured by
    Conner et al. (1997), :data:`CONNER1997_DAMPING`.  With them the
    linear flutter speed (Jones-Wagner RFA) is 9.23 m/s (9.14 m/s
    undamped; 9.19 m/s with the exact Theodorsen function, undamped).
    """
    p = dict(b=0.125, a=-0.5, c=0.5, m_W=1.51, m_T=3.67, x_alpha=0.66, x_beta=0.0028,
             r_alpha2=0.533, r_beta2=0.0055, omega_h=30.72, omega_alpha=17.16, omega_beta=62.58,
             span=1.0, rho=1.0844, zeta=CONNER1997_DAMPING)
    p.update(kw)
    return TypicalSection(**p)


#: Literature values for the Conner / Vasconcellos section.  Only the linear
#: flutter speed is compared quantitatively; the published LCO amplitudes
#: belong to the experimental rig (freeplay +-3.57 deg, measured damping
#: matrix and a fitted 5th-order pitch-stiffness polynomial, experimental
#: U_f 15.2 m/s, numerical 15.1 m/s, LCO onset ~0.9 U_f), which FLAAP does
#: not model, so the freeplay LCO results are cross-checked internally
#: (time marching vs describing function).
VASCONCELLOS2016_REFERENCE = {
    "U_f_linear": 9.2,          # m/s, numerical, linear (no freeplay)
    "freeplay_deg": 3.57,       # experimental rig only (not used by FLAAP's validation)
}


JONES = ((0.165, 0.0455), (0.335, 0.3))     # C(p) ~ 1 - sum a p / (p + beta), p = i k


def fit_theodorsen_lags(nlag: int = 4, kmax: float = 3.0, nk: int = 200):
    """Least-squares lag approximation ``C(p) ~ 1 - sum_l a_l p / (p + beta_l)``
    of the Theodorsen function on ``0 < k <= kmax``, with the poles
    optimised (log-parametrised, so they stay positive).  Returns
    ``((a_1, beta_1), ...)``.  Four lags give a maximum error of about
    1e-4 on ``|C|``; Jones' two-lag fit (:data:`JONES`) gives about 1e-2.
    """
    from scipy.optimize import least_squares
    k = np.geomspace(1e-3, kmax, nk)
    Ck = theodorsen(k)
    p = 1j * k

    def coeffs(lb):
        b = np.exp(lb)
        B = np.array([p / (p + bb) for bb in b]).T
        Ar = np.vstack([B.real, B.imag])
        dr = np.r_[(1 - Ck).real, (1 - Ck).imag]
        a, *_ = np.linalg.lstsq(Ar, dr, rcond=None)
        return a, b, B

    def res(lb):
        a, b, B = coeffs(lb)
        r = 1 - B @ a - Ck
        return np.r_[r.real, r.imag]

    lb0 = np.log(np.geomspace(0.04, 0.8 * kmax, nlag))
    sol = least_squares(res, lb0)
    a, b, _ = coeffs(sol.x)
    order = np.argsort(b)
    return tuple((float(a[i]), float(b[i])) for i in order)


def theodorsen_rfa(sec_or_b, a: float | None = None, c: float | None = None, span: float = 1.0,
                   lags=JONES, gust: bool = True, kvals=None):
    """Exact Roger-form realisation of the Theodorsen flap GAFs.

    The non-circulatory terms are already quadratic in ``p = i k``
    (``A0, A1, A2``).  The circulatory term is
    ``4 pi b C(p) R (S1 + p S2 / b)``. Replacing ``C(p)`` by
    ``1 - sum_l a_l p / (p + beta_l)`` and using
    ``p^2 / (p + beta) = p - beta p / (p + beta)`` gives exactly::

        A0 += K R S1
        A1 += K R S2 / b (1 - sum a_l)
        A_lag,l = -K a_l R (S1 - beta_l S2 / b),      K = 4 pi b span

    so the time-domain model is the classical Theodorsen/Wagner model of
    Conner et al. (1997) with ``lags = JONES``.  ``lags`` is either a
    sequence of ``(a_l, beta_l)`` pairs or an integer ``n`` (fitted with
    :func:`fit_theodorsen_lags`).  The gust column is fitted with the
    same lag roots (:func:`flaap.aero.rfa.fit_rows`).
    """
    from ..aero.rfa import RogerRFA, fit_rows
    if isinstance(sec_or_b, TypicalSection):
        sec = sec_or_b
        b, a, c, span = sec.b, sec.a, sec.c, sec.span
    else:
        b = float(sec_or_b)
    if isinstance(lags, (int, np.integer)):
        lags = fit_theodorsen_lags(int(lags))
    al = np.array([l[0] for l in lags], float)
    bl = np.array([l[1] for l in lags], float)
    Mh, Bh, Kh, R, S1, S2 = theodorsen_flap_matrices(b, a, c)
    T = np.diag([-1.0, 1.0, 1.0])
    Kc = 4 * np.pi * b * span
    A = np.zeros((3 + bl.size, 3, 3))
    # Q_w(p) = -2 p^2 Mh - 2 p b Bh - 2 b^2 Kh + circulatory   (k^2 = -p^2, i k = p)
    A[0] = -2 * b * b * Kh * span + Kc * np.outer(R, S1)
    A[1] = -2 * b * Bh * span + Kc * np.outer(R, S2) / b * (1 - al.sum())
    A[2] = -2 * Mh * span
    for l in range(bl.size):
        A[3 + l] = -Kc * al[l] * np.outer(R, S1 - bl[l] * S2 / b)
    A = np.einsum("ij,tjk,kl->til", T, A, T)
    Ag = None
    if gust:
        k = np.r_[0.0, np.geomspace(0.01, 4.0, 40)] if kvals is None else np.asarray(kvals, float)
        tab = theodorsen_flap_gaf(b, a, c, k, span, gust=True)
        Ag = fit_rows(k, tab.Qg, bl, constrain_k0=True, with_acc=False)
    return RogerRFA(A, bl, Ag, b, 0.0, ["h", "alpha", "beta"], kmax=np.inf)


def section_system(sec: TypicalSection, flap: str = "structural", lags=JONES,
                   K_beta: float | None = None, gust: bool = True):
    """Linear :class:`~flaap.timedomain.AeroelasticSystem` of the section.

    ``flap='structural'``: three structural DOFs ``[h, alpha, beta]``
    with hinge stiffness ``K_beta`` (default nominal), used for freeplay
    and flutter.  ``flap='control'``: two structural DOFs ``[h, alpha]``
    and the flap as a prescribed control input ``beta``, with inertial
    coupling ``M_ctrl``, used for ASE.  ``lags``: Theodorsen lag
    approximation (see :func:`theodorsen_rfa`).
    Returns ``(system, rfa)``.
    """
    from ..timedomain import AeroelasticSystem
    rfa = theodorsen_rfa(sec, lags=lags, gust=gust)
    M, C, K = sec.matrices(K_beta)
    if flap == "structural":
        return AeroelasticSystem(M, K, rfa, C=C), rfa
    if flap == "control":
        return (AeroelasticSystem(M[:2, :2], K[:2, :2], rfa, C=C[:2, :2], n_ctrl=1,
                                  M_ctrl=M[:2, 2:3]), rfa)
    raise ValueError(flap)
