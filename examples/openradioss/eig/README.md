# Radioss `/EIG` mode files: synthetic sample set

**These files are synthetic and written in the OpenRadioss format. They are
not real solver output.**

The public OpenRadioss engine cannot run `/EIG`. Its eigen solver
(`engine/stub/eig.F` and related files) is a stub, so only the commercial
Altair Radioss engine produces eigen results. To exercise the reader, the
cantilever plate deck `../plate_cantilever_0000.rad` (1 m x 1 m x 10 mm Al,
kg-mm-ms units) was solved with FLAAP's FE (`radioss_modal`). The modes were
then written by `flaap.io.openradioss_results.write_radioss_eigen_vtk` in the
layout that `anim_to_vtk` produces from the `<run>A001`, `<run>A002`, ...
animation files of a `/EIG` run.

The VTK files contain:

- `TIME`: the mode frequency, in 1/ms (solver time unit).
- `POINTS`: node coordinates in mm.
- `NODE_ID`: Radioss node ids.
- `Displacement` and `Rotational_DOF`: the mode shapes, normalized to unit
  generalized mass in solver units (`Inorm = 0`).
- `CELLS`: the shell quads, with `ELEMENT_ID`, `PART_ID` and `EROSION_STATUS`
  cell data.

| file | mode | f [Hz] |
|------|------|--------|
| plate_cantileverA001.vtk | 1 | 8.5172 |
| plate_cantileverA002.vtk | 2 | 20.876 |
| plate_cantileverA003.vtk | 3 | 52.640 |
| plate_cantileverA004.vtk | 4 | 67.070 |
| plate_cantileverA005.vtk | 5 | 76.436 |
| plate_cantileverA006.vtk | 6 | 133.95 |

No engine listing (`*_0001.out`) is included. The layout of its eigenvalue
table is not documented publicly, so none was fabricated.
`read_radioss_eigen` takes the frequencies from the VTK `TIME` field. A
listing or a plain `mode frequency` table can be passed as `listing=`
(`frequency_file` in a config).

To use the samples:

```python
from flaap.io import read_radioss_eigen
mm = read_radioss_eigen("examples/openradioss/eig", deck="examples/openradioss/plate_cantilever_0000.rad")
```

The launcher equivalent is `configs/modal/ex7/config_openradioss_eig.txt`.

To regenerate the files, run `python examples/openradioss/eig/make_eig_samples.py`.

Format sources:

- the Altair Radioss Reference Guide, `/EIG` (Starter)
- OpenCourant/OpenCourant (the preserved OpenRadioss source): `genani.F`
  and `engine/stub/eig.F`
- OpenCourant/Tools: `output_converters/anim_to_vtk/src/anim_to_vtk.cpp`

The `flaap/io/openradioss_results.py` docstring has the details.
