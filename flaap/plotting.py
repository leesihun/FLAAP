"""Matplotlib helpers (lattice, mode shapes, pressures, gust histories)."""
from __future__ import annotations

import numpy as np


def _ax(ax, projection=None):
    import matplotlib.pyplot as plt
    if ax is not None:
        return ax
    fig = plt.figure(figsize=(7, 5))
    return fig.add_subplot(111, projection=projection)


def plot_mesh(mesh, ax=None, values=None, cmap="coolwarm", label=None):
    """Plan view of the box lattice, optionally coloured by ``values`` (N,)."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    ax = _ax(ax)
    polys = mesh.corners[:, :, :2]
    pc = PolyCollection(polys, edgecolors="k", linewidths=0.3,
                        facecolors="none" if values is None else None, cmap=cmap)
    if values is not None:
        pc.set_array(np.real(values))
        plt.colorbar(pc, ax=ax, label=label)
    ax.add_collection(pc)
    ax.plot(mesh.colloc[:, 0], mesh.colloc[:, 1], "r.", ms=1.5)
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    return ax


def plot_pressure(mesh, cp, ax=None, part="real"):
    """Box pressure jump ``dCp`` (complex per-mode column or real)."""
    v = getattr(np, part)(cp) if np.iscomplexobj(cp) else cp
    return plot_mesh(mesh, ax, v, label=f"dCp ({part})")


def plot_mode(modal, i: int, ax=None, scale: float | None = None):
    """3-D scatter of the vertical displacement of mode ``i``."""
    ax = _ax(ax, "3d")
    uz = modal.shapes[:, 2, i]
    s = scale if scale is not None else 0.2 * np.ptp(modal.nodes[:, :2]) / max(np.abs(uz).max(), 1e-30)
    p = modal.nodes
    ax.scatter(p[:, 0], p[:, 1], p[:, 2] + s * uz, c=uz, cmap="viridis", s=8)
    lab = modal.labels[i] if modal.labels else f"mode {i + 1}"
    ax.set_title(f"{lab}: {modal.freq_hz[i]:.3f} Hz")
    return ax


def plot_gust(result, labels=None, ax=None):
    """Gust velocity and selected outputs of a :class:`flaap.gust.GustResult`."""
    import matplotlib.pyplot as plt
    labels = labels or [l for l in result.labels if l.startswith("eta:")][:4]
    if ax is None:
        fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    ax[0].plot(result.t, result.wg, "k")
    ax[0].set_ylabel("w_g [m/s]")
    ax[0].set_title(f"1-cos gust H = {result.H:.1f} m, U_ds = {result.Uds:.2f} m/s")
    for lab in labels:
        ax[1].plot(result.t, result.y[:, result.labels.index(lab)], label=lab)
    ax[1].legend(fontsize=7)
    ax[1].set_xlabel("t [s]")
    for a in ax:
        a.grid(alpha=0.3)
    return ax


def plot_fsi(mach, fsi_exp, fsi_calc, ax=None, labels=("experiment", "FLAAP")):
    ax = _ax(ax)
    ax.plot(mach, fsi_exp, "ko-", label=labels[0])
    ax.plot(mach, fsi_calc, "bs--", label=labels[1])
    ax.set_xlabel("Mach")
    ax.set_ylabel("flutter speed index")
    ax.grid(alpha=0.3)
    ax.legend()
    return ax
