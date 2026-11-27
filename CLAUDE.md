# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

FLAAP (Flutter and LCO Aeroelastic Analysis Program) is a Python port and extension of an original Fortran program (kept unchanged in `legacy_fortran/`). The Python package is `flaap/`. Python >= 3.9; core deps are numpy, scipy and matplotlib. Optional deps: `gmsh`/`trimesh` (STL/STEP geometry), `reportlab`/`pillow` (white paper), `pytest`. The repo is not a git repository.

## Commands

```bash
python -m pytest -q tests                        # regression/validation tests
python -m pytest -q tests/test_tsd.py            # one file
python -m pytest -q tests/test_core.py -k name   # one test

python FLAAP_main.py --config configs/flutter/ex1/config_goland_pk.txt   # run a case (= python -m flaap run CFG)
python FLAAP_main.py --config <cfg> --check      # preflight only (also --dry-run, --explain-config, --strict)
python FLAAP_main.py --list-models               # routes
python FLAAP_main.py --describe freeplay         # keys/defaults/modes of one route

python tools/run_all_configs.py [--filter ase] [--skip-slow]   # run every configs/**/config_*.txt -> output/run_all_summary.txt
python tools/gen_config_doc.py                   # regenerate docs/CONFIGURATION.md (do not hand-edit the per-route tables)
python validation/report.py [chapter ...]        # validation/REPORT.md; full run ~40 min; chapters: structures dlm strip tsd_flutter tsd_lifting supersonic uvlm freeplay ase gust trim gaps
python docs/whitepaper/build_whitepaper.py       # docs/FLAAP_Whitepaper.pdf (Korean); reads output/ from run_all_configs and validation/
python studio/start_studio.py                    # FLAAP Studio GUI at http://127.0.0.1:8765
python -m flaap legacy legacy_fortran            # run the original Fortran deck (replaces FLAAP.exe)
```

There is no linter config. `output/` is generated and git-ignored.

## Architecture

**Config-driven launcher** (`flaap/launcher/`) is the main entry point; the same flow serves `FLAAP_main.py`, `python -m flaap run`, `launcher.run_config()` and Studio:

1. `parser.py` reads a flat `key value` file (`key = value` also works; `%` comment lines, `#` inline comments). Indexed keys (`surfaceN`, `controlN`, `stationN`) are handled by the patterns in `specs.py`.
2. `specs.py` is the **single source of truth** for keys, defaults and validation: key catalogues, one `RouteSpec` per `model` (registered via `_route(...)` into `ROUTES`), and validator functions that emit `Diagnostic`s (`diagnostics.py`). `docs/CONFIGURATION.md` is generated from it.
3. `builder.py` turns resolved values into a `Case` (structure/modal model + aero mesh/model): `build_structure` dispatches on the `structure` key (goland, agard445, plate, nastran, ansys, openradioss, openradioss_eig, calculix, stl/geometry, gmsh, legacy, ...).
4. `runner.py` executes the route (`ROUTE_FUNCS` maps model name -> `route_*`) and writes the output contract through `Output`. `ext_routes.py` (freeplay, ase) and `vlm_routes.py` (uvlm) hold the later routes.

Routes (`model`): modal, aero, flutter, gust, trim, maneuver, tsd, freeplay, ase, uvlm. Each has `mode`s documented via `--describe`.

**Adding a route or key** usually means touching `specs.py` (keys, defaults, validators), the route function in `runner.py`/`ext_routes.py`/`vlm_routes.py`, `ROUTE_FUNCS`, a shipped example under `configs/<route>/exN/config_*.txt`, then re-running `tools/gen_config_doc.py`. Studio builds its model form from `/api/models`, which is derived from the same specs.

**Output contract**: every result file is plain text in the config grammar (`key value` lines plus `table ... end` blocks): `results.txt`, `modes.txt`, `cp.txt`, `flutter_<method>.txt`, `timehist_*.txt`, `curves_*.txt`, `tsd_steady.txt`, `vlm_steady.txt`, `figures/*.png`. Read back with `flaap.io.results_txt.read_result`. Studio and the white paper consume these, so changing the format affects both (see `docs/OUTPUT_FORMAT.md`).

**Numerical core** (`flaap/`):
- `geometry` (lifting surfaces, boxes, symmetry) -> `aero/` (`dlm`, `vlm`, `uvlm`, `supersonic`, `tsd`, `tsd_lifting`, `strip`) -> `aero/gaf.py` (`AeroModel`, `AeroModes`: normalwash, GAF tables Q(k), gust column) -> `aero/rfa.py` (Roger/minimum-state RFA).
- `structures/` (beam, MITC4 plate, shell, solid FE, modal model, splines, gmsh mesher, STL/STEP -> FE modal) and `io/` (Nastran `.f06`, ANSYS, OpenRadioss, CalculiX, legacy Fortran decks, results txt) produce a modal model (M, K, mode shapes at nodes). Splines map structural nodes to aero boxes.
- Analyses on top: `flutter` (K, PK, p methods), `timedomain` (state space), `gust`, `trim`, `maneuver`, `loads`, `nonlinear` (freeplay LCO), `ase`.
- The DLM, main program and PK solver sources were missing from the Fortran snapshot (only `.obj`), so those were reimplemented from the program flow and literature; the README package table maps modules to Fortran origins.

**Studio** (`studio/`): stdlib HTTP server + vanilla JS ES modules, no external packages, binds to 127.0.0.1. It never modifies `flaap/` or `configs/`; its runs go to `studio_runs/<name>/` and `output/studio_<name>/`.

`archive/toml_cases/` is the superseded TOML interface (not imported). `legacy_fortran/` is read-only reference material.

## Conventions

- Axes: x downstream, y right, z up; SI units. Harmonic time dependence e^{iωt}; reduced frequency k = ωb/U. ΔCp positive for lift; rotation r_y positive nose-up.
- Config paths in `structure_file` etc. are resolved relative to the config file's directory (`pc.base_dir`).
- Docs and README are bilingual in places (Korean summaries); the white paper is Korean.

## Validation reporting

Validation sections (in `validation/REPORT.md`, the README validation summary, the white paper, and any new validation script output) must be titled by the **capability validated** (e.g. "linear DLM flutter", "TSD transonic flutter"), not by the test case (Goland, AGARD 445.6, ...), which appears as a sub-item. Reference priority: 1) experiment, 2) other independent analyses/codes, 3) analytic formulas only if nothing better exists. Code-internal checks (limits, solver-vs-solver, I/O round trips) go in a separate "consistency" part and are not presented as validation. State explicitly when a reference level is missing.
