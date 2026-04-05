# FLAAP

FLAAP (Flutter and LCO Aeroelastic Analysis Program) is a Windows-based Fortran project for flutter analysis using structural data from MSC Nastran and aerodynamic modeling based on a DLM / 2D potential formulation. The codebase combines structural modal inputs, aerodynamic panel generation, generalized aerodynamic force calculation, and export of flutter result tables such as V-g and V-f curves.

This repository appears to be a research/development snapshot rather than a fully packaged release. It includes source files, a prebuilt executable, example input/output data, and many compiler artifacts.

## What Is In This Repository?

- `FLAAP.exe`: prebuilt Windows executable.
- `*.f90`: core Fortran routines for aerodynamic paneling, generalized aerodynamic forces, utilities, and legacy/experimental modules.
- `20.f06`, `20.dat`: example MSC Nastran files included in the repo.
- `LEE_input.txt`: primary geometry / aerodynamic analysis input.
- `MODESHAPE_input.txt`: modal input data used by the solver.
- `M_K.txt`: structural matrix data.
- `atmosphere_input.txt`: Mach / altitude input.
- `Vg_save.dat`, `Vf_save.dat`, `pk_results_mode*.dat`: sample flutter output files.
- `Conyer's wing`, `Conyer's wing data file`, `Modified Goland wing data file`: additional example cases and data sets.

## Workflow Overview

1. Prepare structural data from MSC Nastran.
   Supply the modal and structural data needed by FLAAP. The repository includes `20.f06` and `20.dat` as example Nastran outputs.
2. Define wing geometry and aerodynamic settings.
   `LEE_input.txt` stores geometry, paneling, reduced-frequency range, Mach number, and Roger approximation settings.
3. Provide atmosphere data.
   `atmosphere_input.txt` provides Mach number and altitude.
4. Run the flutter analysis.
   The solver builds aerodynamic paneling, computes generalized aerodynamic forces, and assembles the flutter equations.
5. Review results.
   `saveV_g.f90` writes `Vg_save.dat` and `Vf_save.dat`; the repository also contains example `pk_results_mode*.dat` files.

## Key Source Files

- `aero_paneling_VLM.f90`: generates aerodynamic panel geometry for the wing planform.
- `calculation.f90`: computes generalized aerodynamic forces and intermediate aerodynamic quantities.
- `saveV_g.f90`: writes velocity-damping and velocity-frequency data files.
- `read_f06.f90`: older helper routine for extracting mode-shape data from Nastran output. The source comments say it is not used in the current program.
- `STATESPACE.f90`: placeholder / future work for state-space or LCO-related analysis; comments say it is not used in the current program.
- `vg_method_old.f90`: older V-g / flutter-development routine.

## Expected Input Files

### `LEE_input.txt`

This file describes the wing and aerodynamic discretization. The companion `LEE_input_readme.txt` indicates that it contains:

1. Gust velocity
2. Gust acceleration
3. Wing `x` length
4. Wing `y` length
5. Planform corner coordinates
6. Thickness
7. Mach number
8. Reduced-frequency range settings
9. Number of Roger lag states
10. Number of aerodynamic panels in `x`
11. Number of aerodynamic panels in `y`

### `MODESHAPE_input.txt`

Mode-shape input used by the solver. Based on the comments in `read_f06.f90`, this data is tied to Nastran modal output, but the current workflow appears to generate or preprocess it outside the visible source tree.

### `M_K.txt`

Structural matrix data used by the flutter solution.

### `atmosphere_input.txt`

`atmosphere_input_readme.txt` shows that this file contains:

1. Mach number
2. Altitude

## Running The Included Executable

The repository includes `FLAAP.exe`, so the simplest workflow is:

1. Keep the required input files in the same working directory as the executable.
2. Replace the sample inputs with your case data while preserving the expected file names.
3. Run `FLAAP.exe`.
4. Inspect `Vg_save.dat`, `Vf_save.dat`, and any `pk_results_mode*.dat` outputs.

Because this codebase uses fixed filenames in the visible routines and the main program source is not included here, it is safest to preserve the original filenames unless you also update the source.

If the executable does not start on a clean machine, you may need the matching Intel Fortran / MKL runtime libraries from the original build environment.

## Build Notes

This snapshot was historically built on Windows with Intel Visual Fortran and Intel MKL. `BuildLog.htm` references:

- Intel Visual Fortran Compiler XE 15.0.1
- x64 Release build
- Intel MKL libraries
- Visual Studio integration on Windows

The original Visual Studio solution / project files are not included in this repository, and the build log references `FLAAP_main.f90`, which is also not present in the current snapshot. Because of that, rebuilding from source may require recovering missing project files or missing source units from the original development environment.

## Current Limitations

- The repository includes many compiled artifacts (`.obj`, `.mod`, manifest files`) alongside source files.
- Some source files are explicitly marked as legacy or not used.
- A full clean rebuild is not documented in this snapshot.
- No repository-wide license file is currently included.

## Example Cases

The repository includes multiple example case folders:

- `Conyer's wing`
- `Conyer's wing data file`
- `Modified Goland wing data file`

These folders appear to provide alternative input sets for testing or comparison studies.

## Reference

Project repo: [leesihun/FLAAP](https://github.com/leesihun/FLAAP)
