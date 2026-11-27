"""Goland wing with a trailing-edge flap (3-D freeplay / ASE example).

The Goland wing (:mod:`flaap.models.goland`) is given an outboard
trailing-edge flap (hinge at 75 % chord, 50-83 % span by default).  The
flap rotation ``delta`` (TE down, rad) is modelled with the rigid
control-surface aero mode :func:`flaap.aero.gaf.control_aero_mode`.  Its
generalised hinge moment is the corresponding row of the GAF matrix.
Two linear systems are built from the same GAF table / RFA:

* ``hinge``: the flap is an elastic DOF ``[q_1..q_n, delta]`` with
  hinge stiffness ``K_h``, hinge inertia ``I_h`` and inertial coupling
  ``M_{i,delta}`` to the wing modes.  This is the model for freeplay
  (:mod:`flaap.nonlinear`).
* ``ctrl``: the flap is a prescribed control input (``n_ctrl = 1``,
  ``M_ctrl = M_{i,delta}``) driven by an actuator.  This is the model
  for aeroservoelasticity (:mod:`flaap.ase`).

Flap inertia, per unit span and about the hinge line ``x_h``: mass
``m_f``, static moment ``S_h = m_f (x_cg,f - x_h)`` (positive aft) and
inertia ``I_h``.  A wing mode ``i`` with beam deflection ``w_i(y)`` (up)
and rotation ``ry_i(y)`` (nose up) moves a flap point by
``z = w_i - (x - x_ea) ry_i - (x - x_h) delta``.  Hence::

    M_{i,delta} = int_flap [ -w_i S_h + ry_i (I_h + (x_h - x_ea) S_h) ] dy
    M_{delta,delta} = I_h l_f

The flap mass is taken to be part of Goland's distributed mass, so the
wing modes are unchanged.  The default hinge stiffness gives an
uncoupled flap frequency ``omega_h = 2 x`` the first torsion frequency.
All numbers are illustrative, not a published configuration.
"""
from __future__ import annotations

import numpy as np

from ..aero.gaf import AeroModel, AeroModes, control_aero_mode
from ..aero.rfa import fit_roger
from ..geometry import ControlSurface, Surface, build_mesh
from ..structures.spline import beam_spline
from ..timedomain import AeroelasticSystem
from .goland import CHORD, SPAN, X_EA, goland_modal

FLAP_DEFAULTS = dict(hinge=0.75, eta=(0.5, 0.8333), m_f=1.8, x_cg_frac=0.4)


def flap_inertia(modal, hinge: float = 0.75, eta=(0.5, 0.8333), m_f: float = 1.8,
                 x_cg_frac: float = 0.4, nq: int = 200):
    """Flap inertia terms: returns ``(M_qd (n,), I_h_total, S_h, I_h)``.

    ``m_f`` is the flap mass per span (kg/m), with its c.g. at
    ``x_cg_frac`` of the flap chord aft of the hinge.  ``I_h`` is about
    the hinge, for a uniform plate: ``m_f c_f^2 / 3``.
    """
    x_h = hinge * CHORD
    c_f = (1 - hinge) * CHORD
    S_h = m_f * x_cg_frac * c_f
    I_h = m_f * c_f ** 2 / 3.0
    y = np.linspace(eta[0] * SPAN, eta[1] * SPAN, nq)
    ny = modal.nodes[:, 1]
    o = np.argsort(ny)
    sh = modal.shapes                     # (nn, 6, nm)
    Mqd = np.zeros(modal.nmodes)
    for i in range(modal.nmodes):
        w = np.interp(y, ny[o], sh[o, 2, i])
        ry = np.interp(y, ny[o], sh[o, 4, i])
        Mqd[i] = (getattr(np, 'trapezoid', None) or np.trapz)(-w * S_h + ry * (I_h + (x_h - X_EA) * S_h), y)
    l_f = (eta[1] - eta[0]) * SPAN
    return Mqd, I_h * l_f, S_h, I_h


def goland_flap(nmodes: int = 4, nchord: int = 8, nspan: int = 12, hinge: float = 0.75,
                eta=(0.5, 0.8333), m_f: float = 1.8, x_cg_frac: float = 0.4,
                omega_flap: float | None = None, kvals=None, nlag: int = 4, gust: bool = True,
                zeta: float = 0.0, zeta_flap: float | None = None, keep_cp: bool = True) -> dict:
    """Build the Goland-flap model.

    Returns a dict with:

    * ``hinge``: :class:`AeroelasticSystem` with the flap as an elastic DOF
      (last coordinate ``'flap'``);
    * ``ctrl``: :class:`AeroelasticSystem` with the flap as a control input;
    * ``table``, ``rfa``, ``modal``, ``mesh``, ``aero``;
    * ``K_h``, ``I_h``, ``M_qd``;
    * ``rho = 1.225``.

    ``omega_flap``: uncoupled flap frequency ``sqrt(K_h / I_h)`` (default
    2 x the first torsion frequency, about 194 rad/s).  ``zeta``: modal
    damping of the wing modes; ``zeta_flap`` (default ``zeta``) is the
    damping of the hinge.
    """
    mm = goland_modal(nmodes, damping=zeta)
    surf = Surface("goland", (0.0, 0.0, 0.0), CHORD, (0.0, SPAN, 0.0), CHORD, nchord=nchord,
                   nspan=nspan, control_surfaces=[ControlSurface("flap", hinge, tuple(eta))])
    mesh = build_mesh([surf], "symmetric")
    modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors()) + control_aero_mode(mesh, "flap")
    am = AeroModel(mesh, modes, 0.0, CHORD / 2)
    if kvals is None:
        kvals = np.r_[0.0, np.geomspace(0.01, 2.0, 16)]
    tab = am.table(kvals, gust=gust, keep_cp=keep_cp)
    rfa = fit_roger(tab, nlag)
    Mqd, Ih, S_h, I_h = flap_inertia(mm, hinge, eta, m_f, x_cg_frac)
    n = mm.nmodes
    Mw = np.diag(mm.gen_mass)
    Kw = np.diag(mm.gen_mass * mm.omega ** 2)
    Cw = np.diag(2 * mm.damping * mm.gen_mass * mm.omega)
    # first torsion: the mode with the largest tip rotation relative to bending
    if omega_flap is None:
        tors = [i for i in range(n) if np.max(np.abs(mm.shapes[:, 4, i])) > np.max(np.abs(mm.shapes[:, 2, i])) / CHORD]
        omega_flap = 2.0 * (mm.omega[tors[0]] if tors else mm.omega[1])
    K_h = Ih * omega_flap ** 2
    zf = zeta if zeta_flap is None else zeta_flap
    M5 = np.zeros((n + 1, n + 1))
    M5[:n, :n] = Mw
    M5[:n, n] = M5[n, :n] = Mqd
    M5[n, n] = Ih
    K5 = np.zeros((n + 1, n + 1))
    K5[:n, :n] = Kw
    K5[n, n] = K_h
    C5 = np.zeros((n + 1, n + 1))
    C5[:n, :n] = Cw
    C5[n, n] = 2 * zf * np.sqrt(K_h * Ih)
    names = [f"mode{i + 1}" for i in range(n)] + ["flap"]
    hinge_sys = AeroelasticSystem(M5, K5, rfa, C=C5, names=names)
    ctrl_sys = AeroelasticSystem(Mw, Kw, rfa, C=Cw, n_ctrl=1, M_ctrl=Mqd[:, None], names=names)
    return dict(hinge=hinge_sys, ctrl=ctrl_sys, table=tab, rfa=rfa, modal=mm, mesh=mesh, aero=am,
                K_h=K_h, I_h=Ih, M_qd=Mqd, omega_flap=omega_flap, rho=1.225, S_h=S_h, I_h_span=I_h)
