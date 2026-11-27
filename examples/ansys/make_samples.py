"""Write synthetic ANSYS modal output files from a known FLAAP FE model.

ANSYS cannot be run here, so the sample files in this folder are produced
from a cantilever aluminium plate (0.5 m chord x 1.0 m span x 4 mm, clamped
at y = 0) solved with the MITC4 plate element of
:mod:`flaap.structures.plate`.  The modes are written in MAPDL / Workbench
layouts by the synthetic writers of :mod:`flaap.io.ansys`; reading them back
must reproduce ``reference_modes.npz`` (see tests/test_structures_io.py).

Files:
  plate_modal.out       MAPDL batch output: Block Lanczos frequency summary,
                        SET,LIST, NLIST, SET,1,n + PRNSOL,U + PRNSOL,ROT per mode
  plate.cdb             CDWRITE file, NBLOCK (3i9,6e21.13e3) + SHELL181 EBLOCK
  plate_prnsol_dof.out  PRNSOL,DOF listing wrapped after 4 columns (no SET,LIST,
                        frequencies from the solver summary; nodes from plate.cdb)
  plate_modal_mm.out    same as plate_modal.out in the mm-tonne-s unit system
  workbench/            Workbench Mechanical exports mode_<i>_Z.txt + frequencies.txt
  reference_modes.npz   nodes, shapes (Nn,6,nm), freq_hz of the source model

Run:  python examples/ansys/make_samples.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from flaap.io.ansys import write_cdb, write_mapdl_listing, write_workbench_modal  # noqa: E402
from flaap.structures.plate import (OrthotropicMaterial, build_plate_model,  # noqa: E402
                                    plate_mesh_trapezoid)

NMODES = 6


def reference_model():
    nodes, quads = plate_mesh_trapezoid((0.0, 0.0), 0.5, (0.0, 1.0), 0.5, 6, 12)
    mat = OrthotropicMaterial.isotropic(71.7e9, 0.33, 2810.0)      # Al 7075-T6
    fe = build_plate_model(nodes, quads, 0.004, mat, clamp_y=0.0, name="cantilever plate")
    return fe.modes(NMODES, name="cantilever plate")


def main(out=HERE):
    mm = reference_model()
    title = "FLAAP cantilever plate 0.5x1.0x0.004 Al7075 - modal"
    write_mapdl_listing(mm, os.path.join(out, "plate_modal.out"), title=title)
    write_cdb(mm, os.path.join(out, "plate.cdb"), title=title)
    write_mapdl_listing(mm, os.path.join(out, "plate_prnsol_dof.out"), title=title,
                        setlist=False, nlist=False, dof=True, wrap=4)
    write_mapdl_listing(mm, os.path.join(out, "plate_modal_mm.out"), title=title + " (mm,t,s)",
                        unit_length=1e-3, unit_mass=1e3)
    write_workbench_modal(mm, os.path.join(out, "workbench"), directions=("Z",))
    np.savez(os.path.join(out, "reference_modes.npz"), nodes=mm.nodes, shapes=mm.shapes,
             freq_hz=mm.freq_hz)
    print(mm.table())


if __name__ == "__main__":
    main()
