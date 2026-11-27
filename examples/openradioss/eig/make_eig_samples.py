"""Write synthetic Radioss ``/EIG`` mode files (``anim_to_vtk`` VTK layout).

The public OpenRadioss engine has no eigen solver (``/EIG`` is a stub
there) and no Radioss eigen results are available here, so these sample
files are *synthetic*: the cantilever plate deck ``../plate_cantilever_0000.rad``
(kg-mm-ms) is solved with FLAAP's own FE (:func:`flaap.io.openradioss.radioss_modal`)
and the modes are written by :func:`flaap.io.openradioss_results.write_radioss_eigen_vtk`
in the layout ``anim_to_vtk`` produces from ``<run>A001 ...`` animation
files of a ``/EIG`` run (unit generalised mass, ``Inorm = 0``, solver
units, frequency [1/ms] in the ``TIME`` field, ``Displacement`` and
``Rotational_DOF`` vectors).

Files:
  plate_cantileverA001.vtk ... A006.vtk   modes 1-6

Run:  python examples/openradioss/eig/make_eig_samples.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", ".."))

from flaap.io.openradioss import radioss_modal  # noqa: E402
from flaap.io.openradioss_results import write_radioss_eigen_vtk  # noqa: E402

DECK = os.path.join(HERE, "..", "plate_cantilever_0000.rad")
NMODES = 6


def main(out=HERE):
    mm = radioss_modal(DECK, nmodes=NMODES)
    paths = write_radioss_eigen_vtk(mm, out, "plate_cantilever", deck=DECK)
    print(mm.table())
    return paths


if __name__ == "__main__":
    main()
