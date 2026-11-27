# gmsh mesher examples

With `structure gmsh`, gmsh becomes the mesher of the FLAAP FE model
(`flaap.structures.gmsh_mesher`). It accepts three kinds of input:

- **`.geo` script.** Sizes, `Transfinite`, `Recombine`, `Physical` groups and
  OpenCASCADE booleans are used as written. If the script does not mesh by
  itself, it is meshed in its highest dimension.
- **`.msh` file.** Any MSH version gmsh reads (2.2 or 4.1, ASCII or binary).
  The mesh is used as is.
- **CAD file.** STEP / IGES / BREP, imported through OpenCASCADE and meshed.

| file | model |
|------|-------|
| `plate_wing.geo` | 0.5 × 1.0 m aluminium plate, transfinite 6 × 12, `Recombine` → 72 quad4 (MITC4) shells, `Physical Curve("Root")` |
| `box_beam.geo` | OCC box 0.10 × 1.0 × 0.02 m, unstructured tetrahedra (`mesh_order 2` → tet10), `Physical Surface("Root")` |

## Config keys

```
structure        gmsh
structure_file   plate_wing.geo        # .geo | .msh | .step/.stp | .iges/.igs | .brep
thickness        0.004                 # shells; group_thickness Skin, 0.004, Spar, 0.006
clamp_group      Root                  # physical group(s) clamped (else: clamp y, 0.0)
mesh_size        0.02                  # Mesh.MeshSizeMax (when gmsh meshes)
mesh_order       2                     # 1: tri3/quad4/tet4/hex8, 2: tet10
recombine        true                  # quad4 surface elements
mesh_dim         2                     # 2 shells | 3 solids (default: highest entity dimension)
mesh_algorithm   8                     # Mesh.Algorithm
mesh_save        out/plate.msh         # also write the FE mesh
```

How a mesh becomes the FE model:

- **Element dimension.** Only elements of the highest dimension are used.
- **Physical groups.** When the model has physical groups of that dimension,
  only their elements are used (gmsh's own save rule). Each group becomes an
  element block, so `group_thickness` can give each group its own thickness.
- **Element types.**
  - tri6 / quad8 / quad9 / hex20 are reduced to their corner nodes, with a
    warning.
  - Prisms and pyramids are rejected, since the FLAAP FE has no such
    elements.

Preflight checks the following:

- the gmsh module is importable (`CFG-DEP-001`);
- the file type is accepted (`CFG-VALUE-007`);
- `clamp_group` / `group_thickness` names appear in the `.geo` / `.msh`
  (`CFG-PATH-004`).

## Command line

```
python -m flaap mesh examples/gmsh/plate_wing.geo -o plate.msh
python -m flaap mesh part.step --size 0.01 --order 2
```

These commands print the element blocks, the bounding box and the physical groups.

## Configs and validation

- `configs/modal/ex9/config_gmsh_plate.txt`: f = 3.3714, 14.3225, 21.3217,
  47.313, 61.759, 93.440 Hz. These equal the plate FE reference of
  `examples/calculix` (Δf 3·10⁻⁵ %).
- `configs/modal/ex10/config_gmsh_solid.txt`: tet10 bar, f1 = 16.45 Hz
  against 16.32 Hz from Euler-Bernoulli (+0.8 %, plate effect of b/h = 5).
- `configs/flutter/ex12/config_gmsh_pk.txt`: DLM + PK, V_f = 70.585 m/s. This
  is the same as the CalculiX `.frd` model of `flutter/ex11`.
- `validation/structures_validation.py` → `gmsh`
