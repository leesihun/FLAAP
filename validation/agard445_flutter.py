"""Validation: AGARD 445.6 weakened model 3 flutter boundary (DLM + PK).

Structure: MITC4 orthotropic plate (:mod:`flaap.models.agard445`), first
4 modes.  Aerodynamics: doublet lattice at each test Mach number
(subsonic linear theory - transonic/supersonic points are run for
reference only; the DLM cannot capture the transonic dip and is not valid
for M > 1).

Mass ratio: the plate model is lighter than the test article
(1.70 vs 1.867 kg), so the PK solution is run at the *equivalent* density
``rho_eq = rho_test * m_model / m_test`` which reproduces the test mass ratio
``mu``; the flutter speed index ``FSI = V / (b_s w_alpha sqrt(mu))`` is then
directly comparable with the experiment.

Run:  python validation/agard445_flutter.py [--fine] [--plot out.png]
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from flaap.aero.gaf import AeroModel, AeroModes, default_k_list  # noqa: E402
from flaap.flutter import identify, pk_method  # noqa: E402
from flaap.geometry import build_mesh  # noqa: E402
from flaap.models.agard445 import (B_REF, EXPERIMENT_FREQ_HZ, FLUTTER_DATA, OMEGA_ALPHA,  # noqa: E402
                                   agard445_fe, agard445_surface, frustum_volume)
from flaap.structures.spline import surface_spline  # noqa: E402


def run(fine: bool = False, machs=None, verbose=True, exp_freq: bool = False):
    """``exp_freq``: replace the FE natural frequencies by the measured ones
    (FE mode shapes kept) - the usual practice in published comparisons."""
    fe = agard445_fe(n_chord=16, n_span=16)
    modal = fe.modes(4, name="AGARD 445.6")
    if exp_freq:
        modal.omega = 2 * np.pi * EXPERIMENT_FREQ_HZ.copy()
    # panel mass from a rigid plunge field (valid for consistent mass matrices)
    u = np.zeros(fe.ndof)
    u[2::6] = 1.0
    m_model = float(u @ (fe.M @ u))
    m_test = 14.315 * frustum_volume()
    nc, ns = (12, 16) if fine else (8, 10)
    mesh = build_mesh([agard445_surface(nc, ns)], "symmetric")
    sp = surface_spline(mesh, modal.nodes, method="ips")
    modes = AeroModes.from_spline(sp, modal.dof_vectors(), [f"m{i + 1}" for i in range(4)])
    kv = default_k_list(1.2, 14)
    kv[0] = 1e-4
    rows = []
    if verbose:
        print("model frequencies [Hz]:", np.round(modal.freq_hz, 2), " exp:", EXPERIMENT_FREQ_HZ)
        print(f"panel mass model {m_model:.4f} kg, test {m_test:.4f} kg")
        print("   M    FSI_exp  FSI_dlm   err%   wf/wa_exp  wf/wa_dlm   V_dlm[m/s]  (rho_eq)")
    for M, d in FLUTTER_DATA.items():
        if machs is not None and M not in machs:
            continue
        t0 = time.perf_counter()
        am = AeroModel(mesh, modes, min(M, 0.99), B_REF)
        table = am.table(kv, gust=False, keep_cp=False)
        rho_eq = d["rho"] * m_model / m_test
        V = np.linspace(50.0, 700.0, 261)
        res = pk_method(modal.M, modal.K, table, rho_eq, V)
        identify(res)
        fl = [p for p in res.points if p.kind == "flutter"]
        c = min(fl, key=lambda p: p.velocity) if fl else None
        if c is None:
            rows.append((M, d["FSI"], np.nan, d["freq_ratio"], np.nan, np.nan))
            continue
        fsi = c.velocity / (B_REF * OMEGA_ALPHA * np.sqrt(d["mu"]))
        fr = 2 * np.pi * c.freq_hz / OMEGA_ALPHA
        rows.append((M, d["FSI"], fsi, d["freq_ratio"], fr, c.velocity))
        if verbose:
            tag = "  (M>1: DLM not valid)" if M > 1 else ("  (transonic)" if M > 0.85 else "")
            print(f"{M:6.3f}  {d['FSI']:7.4f}  {fsi:7.4f}  {100 * (fsi / d['FSI'] - 1):6.1f}"
                  f"   {d['freq_ratio']:8.4f}  {fr:9.4f}   {c.velocity:9.1f}  ({rho_eq:.4f})"
                  f"{tag}   [{time.perf_counter() - t0:.1f} s]")
    return np.array(rows, float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", action="store_true")
    ap.add_argument("--plot", default=None)
    ap.add_argument("--exp-freq", action="store_true", help="use the measured natural frequencies")
    a = ap.parse_args()
    r = run(a.fine, exp_freq=a.exp_freq)
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(r[:, 0], r[:, 1], "ko-", label="experiment (Yates 1987)")
        ax[0].plot(r[:, 0], r[:, 2], "bs--", label="FLAAP DLM + PK")
        ax[0].set_xlabel("Mach"); ax[0].set_ylabel("flutter speed index"); ax[0].legend()
        ax[1].plot(r[:, 0], r[:, 3], "ko-", r[:, 0], r[:, 4], "bs--")
        ax[1].set_xlabel("Mach"); ax[1].set_ylabel(r"$\omega_f/\omega_\alpha$")
        for x in ax:
            x.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=120)
        print("saved", a.plot)


if __name__ == "__main__":
    main()
