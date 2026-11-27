# CalculiX `*FREQUENCY` results: sample set

**These `.frd` / `.dat` files are synthetic. They are written in the CalculiX
format, but they are not real `ccx` output.** No `ccx` binary is bundled. The
cantilever plate (0.5 m × 1.0 m × 4 mm Al, 6 × 12 MITC4 shells, the same plate
as `../ansys`) was solved with FLAAP's FE. The modes were then written by
`flaap.io.calculix.write_frd` / `write_calculix_dat` in the layout that `ccx`
uses. Run `python make_samples.py` to regenerate them.

| file | contents |
|------|----------|
| `plate.frd` | long node format (I10), kg-m-s, 6 modes, `qu4` elements |
| `plate.dat` | `E I G E N V A L U E   O U T P U T` table for `plate.frd` |
| `plate_mm_short.frd` | short node format (I5), mm-tonne-s, a static step first and 3 rigid-body modes before the 6 elastic ones |
| `plate.inp` | the equivalent `ccx` input deck (`S4` shells, `*FREQUENCY`, `*NODE FILE U`) |
| `reference_modes.npz` | the FLAAP FE modes, used as the reference by the tests |

Frequencies [Hz]: 3.3714, 14.3225, 21.3217, 47.3130, 61.7587, 93.4401.

## Using a real CalculiX run

```
ccx plate            # reads plate.inp, writes plate.frd + plate.dat
```

The step needs `*NODE FILE` with `U`, so that the `DISP` blocks are written.
The `.frd` must be ASCII. This is the `ccx` default. A binary `.frd` (cgx
`OUTPUT` binary, format flag 2/3) is rejected in preflight with
CFG-PATH-003.

```
structure        calculix
structure_file   plate.frd
frequency_file   plate.dat      # optional: the .dat next to the .frd is used by default
length_scale     0.001          # model in mm-tonne-s
mass_scale       1000
```

What the reader does with the file:

- **Mode selection.** The reader takes the `DISP` blocks of the last
  `*FREQUENCY` step. A new step starts where the frequency drops. Modes below
  `skip_below_hz` (1 mHz, so the rigid-body modes) are skipped.
- **Frequencies.** The `.dat` eigenvalue table (`f = √λ / 2π`) replaces the
  `E12.5` header frequencies when it matches them. It can match at the start
  of the table or at its tail. When it does not match, the reader warns and
  keeps the `.frd` values.
- **Units.** `ccx` has no unit system. Modes are mass-normalised in the model's
  units, so the generalized mass in SI is `mass_scale · length_scale²`.
- **Shapes.** `ccx` writes translations only, so the rotational DOFs are zero.
- **Display.** The display surface is built from the `.frd` elements:
  - solid boundary faces;
  - shells split into triangles;
  - beams as edges.

Configs: `configs/modal/ex8/config_calculix_modal.txt` (modal) and
`configs/flutter/ex11/config_calculix_pk.txt` (DLM + PK, mm-t-s short file).
Validation: `validation/structures_validation.py` → `calculix`.
