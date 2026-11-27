"""Write synthetic CalculiX modal result files from a known FLAAP FE model.

CalculiX (ccx) is not run here.  The samples come from the same cantilever
aluminium plate as examples/ansys (0.5 m chord x 1.0 m span x 4 mm, clamped
at y = 0, MITC4 plate element of :mod:`flaap.structures.plate`), written in
the ccx layouts by the sample writers of :mod:`flaap.io.calculix`.  Reading
them back must reproduce ``reference_modes.npz`` (tests/test_structures_io.py).

Files:
  plate.frd           ASCII .frd, long node format: node block, qu4 element
                      block, one DISP block per mode (ICTYPE 2, VALUE = Hz)
  plate.dat           eigenvalue table (more digits than the .frd header)
  plate_mm_short.frd  short node format, mm-tonne-s units, a static DISP step
                      first and 3 zero-frequency modes (must be skipped)
  plate.inp           S4 shell deck with *FREQUENCY + *NODE FILE U, for
                      running ccx yourself (ccx plate -> plate.frd / plate.dat)
  reference_modes.npz nodes, shapes (Nn,6,nm), freq_hz of the source model

Run:  python examples/calculix/make_samples.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))

from flaap.io.calculix import write_calculix_dat, write_calculix_inp, write_frd  # noqa: E402
from flaap.structures.plate import (OrthotropicMaterial, build_plate_model,  # noqa: E402
                                    plate_mesh_trapezoid)

NMODES = 6


def reference_model():
    nodes, quads = plate_mesh_trapezoid((0.0, 0.0), 0.5, (0.0, 1.0), 0.5, 6, 12)
    mat = OrthotropicMaterial.isotropic(71.7e9, 0.33, 2810.0)      # Al 7075-T6
    fe = build_plate_model(nodes, quads, 0.004, mat, clamp_y=0.0, name="cantilever plate")
    return fe.modes(NMODES, name="cantilever plate"), np.asarray(nodes), np.asarray(quads)


def main(out=HERE):
    mm, nodes, quads = reference_model()
    elements = [(9, [int(c) + 1 for c in q]) for q in quads]        # frd type 9 = qu4
    write_frd(mm, os.path.join(out, "plate.frd"), elements=elements)
    write_calculix_dat(mm, os.path.join(out, "plate.dat"))
    write_frd(mm, os.path.join(out, "plate_mm_short.frd"), elements=elements, long=False,
              unit_length=1e-3, unit_mass=1e3, static_step=True, rigid_modes=3)
    xyz = nodes if nodes.shape[1] == 3 else np.c_[nodes, np.zeros(len(nodes))]
    root = np.flatnonzero(np.abs(xyz[:, 1]) < 1e-9)
    write_calculix_inp(xyz, quads, os.path.join(out, "plate.inp"), thickness=0.004, E=71.7e9, nu=0.33,
                       rho=2810.0, clamped_nodes=root, nmodes=NMODES)
    np.savez(os.path.join(out, "reference_modes.npz"), nodes=mm.nodes, shapes=mm.shapes,
             freq_hz=mm.freq_hz)
    print(mm.table())


if __name__ == "__main__":
    main()
