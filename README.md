# FLAAP — Flutter and LCO Aeroelastic Analysis Program (Python edition)

FLAAP began as a Fortran program (Intel Visual Fortran + MKL, 2017–2018) for DLM-based flutter analysis. This repository now also contains a **complete Python port**, `flaap/`, with substantial extensions:

- DLM unsteady aerodynamics and a transonic small-disturbance (TSD) solver.
- **VLM / UVLM**: steady vortex lattice (CL, induced drag, span load), nonlinear time-domain UVLM (impulsive start, harmonic pitch / plunge, sharp-edged and 1-cos gusts, free wake) giving lift / drag / moment histories, and a linearised UVLM for GAFs and discrete-time state-space flutter.
- **Supersonic aerodynamics (M > 1)**: piston theory and a Mach-box method (Ackeret / Evvard limits) for supersonic flutter.
- K, PK and p flutter methods with automatic flutter and divergence identification.
- RFA state-space models with discrete / continuous gust and CS-25 manoeuvre loads.
- Static aeroelastic trim, divergence and control effectiveness.
- **Control-surface freeplay LCO**: exact piecewise-linear time marching, describing function, bifurcation sweeps.
- **Aeroservoelasticity (ASE)**: actuator / sensor models, active flutter suppression, gust load alleviation, Nyquist / Bode margins.
- **Structures**: beam, MITC4 plate, flat shell (quad4 / tri3) and solid (tet4 / tet10) FE.
- **Model import**: legacy Fortran decks, MSC Nastran `.f06`, **ANSYS** modal results (MAPDL listings, PRNSOL, Workbench text, `.cdb`, `.rst`), **OpenRadioss / Radioss** starter decks and `/EIG` modal results (`anim_to_vtk` VTK files), **CalculiX** `*FREQUENCY` results (ASCII `.frd` + `.dat`), **STL / STEP geometry + material → FE modal analysis**, and a **gmsh mesher** (`.geo` scripts, `.msh` meshes, STEP → shell / solid FE with physical-group clamps and thicknesses).
- **AI-CAE4ALL-style launcher**: plain `key value` config files with preflight diagnostics. Results are human-readable `key value` text files in the same grammar (`results.txt`, `modes.txt`, `flutter_*.txt`, …).
- **FLAAP Studio**: a local browser GUI to draw planforms CAD-style, run cases and view modes, ΔCp, V-g / V-f, time histories and curves.
- **White paper** ([docs/FLAAP_Whitepaper.pdf](docs/FLAAP_Whitepaper.pdf)): theory, capabilities and validation of every feature.

The original Fortran sources, executable and inputs are kept unchanged. The previous README is kept as [README_fortran_original.md](README_fortran_original.md).

> 한국어 요약: 기존 Fortran FLAAP을 Python으로 완전히 포팅하고 다음을 추가했습니다.
> - DLM/TSD 공력, K/PK/p 플러터 자동 판별, RFA 시간영역 돌풍·기동 하중, 정적 트림/발산
> - VLM(정상)·UVLM(시간영역 양력/항력/모멘트: 급출발, 조화 피치/플런지, 돌풍, 자유후류)·선형 UVLM 플러터, 초음속 공력(피스톤 이론, Mach box)
> - 조종면 freeplay LCO (시간적분 + 기술함수 + 분기 선도), ASE (능동 플러터 억제, 돌풍 하중 완화, 안정여유)
> - 쉘/고체 FE, Nastran·ANSYS·OpenRadioss(starter 덱 + /EIG 모달 결과)·CalculiX(.frd) 연동, STL/STEP 형상 + 물성 → 모달 해석 → 플러터, gmsh 메셔(.geo/.msh/STEP → 쉘·고체 FE, Physical 그룹 구속)
> - AI-CAE4ALL 형식 평문 설정 파일 + 사전검증 (`FLAAP_main.py --config ...`), 입력과 같은 key-value 문법의 사람이 읽기 쉬운 TXT 결과 파일
> - FLAAP Studio GUI (CAD형 평면도 편집, 모드 애니메이션, ΔCp, V-g/V-f, 시간응답)
> - 이론·검증 백서 PDF (`docs/FLAAP_Whitepaper.pdf`)

## Quick start

```bash
pip install numpy scipy matplotlib         # core
pip install gmsh trimesh                   # optional: STL / STEP geometry -> FE modal
pip install reportlab pillow pytest        # optional: white paper PDF, tests

# AI-CAE4ALL-style launcher (configs/<route>/exN/config_*.txt)
python FLAAP_main.py --config configs/flutter/ex1/config_goland_pk.txt
python FLAAP_main.py --config <cfg> --check          # preflight only (also --dry-run, --explain-config, --strict)
python FLAAP_main.py --list-models                   # routes: modal aero flutter gust trim maneuver tsd freeplay ase uvlm
python FLAAP_main.py --describe freeplay             # keys, defaults and modes of one route
python tools/run_all_configs.py [--filter ase]       # run every shipped config -> output/run_all_summary.txt

python studio/start_studio.py                        # FLAAP Studio GUI at http://127.0.0.1:8765
python docs/whitepaper/build_whitepaper.py           # rebuild docs/FLAAP_Whitepaper.pdf from output/ + validation/

python -m pytest -q tests                            # regression / validation tests
python validation/report.py                         # full validation report -> validation/REPORT.md
python -m flaap legacy legacy_fortran               # original Fortran deck (replaces FLAAP.exe) -> output/legacy_legacy_fortran
python -m flaap run configs/flutter/ex1/config_goland_pk.txt   # = FLAAP_main.py --config ...
python -m flaap mesh examples/gmsh/plate_wing.geo -o plate.msh  # gmsh mesher: element blocks + physical groups
```

A config file is one `key value` per line (`key = value` also works); `%` starts a comment line and `#` an inline comment:

```text
% STL wing: geometry + material -> tet10 modal -> DLM + PK flutter
model        flutter
mode         pk
structure         stl
structure_file    ../../../examples/geometry/wing.stl
fe_type           solid
element           tet10
mesh_size         0.03
material          custom
young             0.4e9
poisson           0.3
material_density  120
clamp             y, 0.0
mach              0.3
bref              0.1625
velocities        10, 300, 146
surface1          wing, 0, 0, 0, 0.40, 0.2116, 1.2, 0, 0.25, 8, 12
```

Every key, default and diagnostic code is documented in [docs/CONFIGURATION.md](docs/CONFIGURATION.md) (generated by `tools/gen_config_doc.py`); the output files (`results.txt`, `modes.txt`, `cp.txt`, `flutter_*.txt`, `timehist_*.txt`, `curves_*.txt`, read back with `flaap.io.results_txt.read_result`) in [docs/OUTPUT_FORMAT.md](docs/OUTPUT_FORMAT.md). Studio usage is in [studio/README.md](studio/README.md).

### Shipped example configs

| route | examples |
|---|---|
| `modal` | Goland beam, cantilever plate, AGARD 445.6, ANSYS listing, OpenRadioss deck, STL wing (tet10), Radioss /EIG results (ex7), CalculiX `.frd` + `.dat` (ex8), gmsh `.geo` plate quad4 (ex9) and OCC bar tet10 (ex10) |
| `aero` | Goland GAF table, steady plate with aileron |
| `flutter` | Goland PK / all methods, AGARD 445.6, plate, Nastran `.f06`, legacy deck, STL wing, AGARD 445.6 at M 1.141 (Mach box, ex8), plate at M 2 (piston theory vs Mach box, ex9), AGARD 445.6 at M 0.499 with measured frequencies (ex10), CalculiX `.frd` plate in mm-t-s (ex11), gmsh `.geo` plate (ex12) |
| `gust` | Goland discrete 1-cos envelope, continuous von Kármán turbulence |
| `trim` / `maneuver` | demo aircraft 1-g trim, Goland divergence, checked / abrupt pitch manoeuvres |
| `tsd` | AGARD 445.6 steady M 0.954, TSD + PK flutter, NACA 0012 wing incidence sweep at M 0.8 (lifting mean flow: CL and shock positions vs alpha, ex3), AGARD 445.6 flutter about a lifting mean flow at alpha 1 deg (ex4) |
| `freeplay` | typical-section DF + time marching, time history at 6 m/s, Goland wing + flap DF |
| `ase` | typical-section flutter suppression, loop margins, Goland + flap gust load alleviation |
| `uvlm` | steady VLM rectangle (ex1), impulsive start (ex2), plunge thrust (ex3), pitch vs Theodorsen (ex4), 1-cos gust (ex5), Goland linear-UVLM flutter (ex6) |

Python API example (Goland wing):

```python
import numpy as np
from flaap.models.goland import goland_modal, goland_surface, CHORD
from flaap.geometry import build_mesh
from flaap.structures.spline import beam_spline
from flaap.aero.gaf import AeroModel, AeroModes
from flaap.flutter import pk_method
from flaap.aero.rfa import fit_roger
from flaap.timedomain import modal_system
from flaap import gust

mm   = goland_modal(4, damping=0.02)
mesh = build_mesh([goland_surface(6, 12)], "symmetric")
am   = AeroModel(mesh, AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors()), 0.0, CHORD / 2)
tab  = am.table(np.r_[0, np.geomspace(0.01, 2, 19)])          # GAFs Q(k) + gust column
res  = pk_method(mm.M, mm.K, tab, 1.225, np.linspace(60, 260, 101), C=mm.C())
print(res.summary())                                            # flutter ~153 m/s, 70.7 rad/s

sysm = modal_system(mm, fit_roger(tab, 4))                      # Roger RFA state space
ss   = sysm.state_space(100.0, 1.225)                           # inputs: w_g, w_g'
r    = gust.discrete_gust(ss, 100.0, H=20.0, Uds=10.0)          # 1-cos gust time history
```

## Package layout

| module | content | Fortran origin |
|---|---|---|
| `flaap.geometry` | trapezoidal lifting surfaces, control surfaces, box lattice, symmetry | `aero_paneling_VLM.f90` |
| `flaap.aero.vlm` | steady vortex lattice (horseshoes, Prandtl–Glauert) | `aero_paneling_VLM.f90` |
| `flaap.aero.dlm` | doublet-lattice method (Albano–Rodden, quartic kernel, Laschka) | `DLM.obj` (source missing; reimplemented) |
| `flaap.aero.gaf` | aero modes, normalwash, GAF tables `Q(k)`, gust column, pressures | `calculation.f90` |
| `flaap.aero.rfa` | Roger RFA (auto-weighted, constrained), minimum-state RFA | `generalizedforce.obj`, `STATESPACE.f90` |
| `flaap.aero.strip` | Theodorsen / Sears strip theory | – (new) |
| `flaap.aero.uvlm` | steady VLM (Trefftz drag, span load), nonlinear UVLM (time-domain CL / CD / Cm, gusts, free wake), linearised UVLM (GAFs, discrete-time state-space flutter) | – (new) |
| `flaap.aero.supersonic` | piston theory and Mach-box supersonic GAFs | – (new) |
| `flaap.aero.tsd` | 3-D transonic small-disturbance solver (steady nonlinear + time-linearised harmonic) | – (new) |
| `flaap.aero.tsd_lifting` | lifting TSD: nonlinear mean flow with incidence and NACA 4-digit camber (upper / lower fields, wake jump), CL and shock positions vs alpha, GAFs about the lifting mean flow | – (new) |
| `flaap.flutter` | K (V-g), PK and p methods, mode tracking, flutter / divergence / hump identification, matched point, legacy output files | `flutter_k_method`, `flutter_P-K_method`, `VG`, `saveV_g.f90`, `vg_method_old.f90` |
| `flaap.timedomain` | aeroelastic state space (`A, B, C, D`) with gust and control inputs and load outputs | `generate_A_matrices`, `STATESPACE.f90` |
| `flaap.gust` | CS-25 1-cos gusts, gust envelope, von Kármán / Dryden turbulence (`A_bar`, `N0`), frequency-domain check | – (new) |
| `flaap.trim` | static aeroelastic trim, divergence, control effectiveness, reversal | – (new) |
| `flaap.maneuver` | rigid + elastic aircraft model, short period, CS-25.331(c) checked and abrupt pitch manoeuvres | – (new) |
| `flaap.loads` | monitor stations, force-summation and mode-displacement loads | – (new) |
| `flaap.structures` | modal model, beam FE, MITC4 orthotropic plate FE, flat shell quad4 / tri3, solid tet4 / tet10 / hex8, generic FE assembly, STL / STEP → gmsh → FE modal, gmsh mesher for `.geo` / `.msh` / CAD (`gmsh_mesher`), IPS / TPS / RBF / beam / legacy splines | `TPS_lib.f90`, `Thin_plate_spline.f90`, `calculation.f90` |
| `flaap.nonlinear` | freeplay / bilinear / preloaded hinge springs, exact-switching time marching, LCO detection, describing function, bifurcation sweeps | – (new) |
| `flaap.ase` | actuators, sensors, static / dynamic controllers, closed-loop flutter, Nyquist / Bode margins, gust load alleviation | – (new) |
| `flaap.io.ansys`, `flaap.io.openradioss`, `flaap.io.openradioss_results`, `flaap.io.calculix` | ANSYS modal results (MAPDL listing, PRNSOL, Workbench, `.cdb`, `.rst`), OpenRadioss starter decks → FE / modal model, Radioss `/EIG` mode shapes + frequencies from `anim_to_vtk` files, CalculiX `*FREQUENCY` `.frd` (short / long ASCII) + `.dat` eigenvalue table | – (new) |
| `flaap.launcher` | AI-CAE4ALL-style config parser, preflight diagnostics, route specs, runners, output contract (`FLAAP_main.py`) | – (new) |
| `flaap.io.legacy` | `LEE_input.txt`, `M_K.txt`, `MODESHAPE_input.txt`, `atmosphere_input.txt` and a full legacy run | `input*.obj`, `FLAAP_main.obj` |
| `flaap.io.nastran` | GRID (small, large, free field), `.f06` eigenvalues and eigenvectors → modal model | `read_f06.f90` |
| `flaap.linalg`, `flaap.atmosphere` | complex LU inverse and determinant, eigen sort, MAC; ISA atmosphere | `c_inverse`, `C_det1`, `c_eigen`, `c_matrix_sort_imag`, `transpose_vector`, `air_properties`, `altitute` |
| `flaap.models` | Goland wing, AGARD 445.6 weakened model 3, demo flexible aircraft | – |
| `flaap.launcher`, `flaap.cli` | key-value `.txt` configs (parser, preflight, routes, output contract) and the `python -m flaap` command line | – |
| `studio/` | FLAAP Studio local web GUI (stdlib HTTP server + vanilla JS) | – (new) |
| `docs/whitepaper/` | white paper generator (reportlab) | – (new) |

Conventions:
- Axes: x downstream, y right, z up. SI units.
- Harmonic time dependence e^{iωt}; reduced frequency k = ωb/U.
- The pressure jump ΔCp is positive for lift.
- Rotation r_y is positive nose-up.

## Input files

Every input is a flat `key value` text file (`configs/<route>/exN/*.txt`, see [docs/CONFIGURATION.md](docs/CONFIGURATION.md)). `python -m flaap run CONFIG.txt` is the same as `python FLAAP_main.py --config CONFIG.txt`. The model files that a config points to (`structure_file`, `.f06`, `.cdb`, `.vtk`, `.stl`, ...) are read in their native formats.

## Validation (summary — full tables and figures in [validation/REPORT.md](validation/REPORT.md))

Each item is named after the **capability** validated; the test case (Goland wing, AGARD 445.6, ...) is the evidence underneath it. References are used in priority order **experiment → other independent analyses/codes → analytic solutions**; code-internal checks (limits, solver-vs-solver, file round trips) are listed as *consistency* and are not counted as validation. Missing levels are stated.

| capability | experiment | other analyses | analytic | headline |
|---|---|---|---|---|
| Structural modes (plate/shell/beam/solid FE) | AGARD 445.6 measured frequencies | Leissa plate, Yates FE | Goland beam, Timoshenko | plate +0.0 / +5.0 / +0.8 / −1.6 % vs test; shells ≤ 0.6 % vs Leissa; beams ≤ 0.15 % |
| Linear DLM flutter (subsonic) | AGARD 445.6, M ≤ 0.9 | Goland 3-D (5 published results), Fortran FLAAP.exe | none | FSI +1.2 … +3.3 % vs test; Goland within published range |
| Strip-theory (Theodorsen) flutter | none | Goland 1945, Patil 2001 | – | V_f −0.2 % |
| TSD transonic flutter | AGARD 445.6, M 0.5–0.96 | not yet compared | n/a | transonic dip: +0.2 / −3.5 % (DLM +10.6 / +8.6 %) |
| Lifting TSD mean flow (shock position) | **gap** | AGARD-AR-211 Euler (NACA 0012) | Helmbold | upper shock 0.035 c ahead of Euler at equal cl |
| Supersonic (Mach box, piston) | AGARD 445.6, M 1.07 / 1.14 | Dowell panel flutter | Ackeret, Garrick–Rubinow | +8.5 / +13.9 % vs test; ≤ 0.05 % vs Dowell |
| VLM / UVLM aero and UVLM flutter | none | Murua 2012 UVLM (Goland) | Wagner, Küssner, Theodorsen, Garrick | V_f 167.0 vs 165.0 m/s (+1.2 %) |
| Freeplay LCO | **gap** (LCO amplitudes) | Vasconcellos 2016 linear V_f | none | V_f 9.23 vs 9.2 m/s |
| ASE, gust, trim / manoeuvre | none | none | none | consistency only |

**Linear DLM flutter.** *Experiment* — AGARD 445.6 weakened 3 (Yates 1987, NASA TM-100492; measured frequencies 9.60 / 38.10 / 50.70 / 98.50 Hz with FE mode shapes):

| Mach | FSI exp | DLM | err DLM |
|---|---|---|---|
| 0.499 | 0.4459 | 0.4512 | +1.2 % |
| 0.678 | 0.4174 | 0.4313 | +3.3 % |
| 0.901 | 0.3700 | 0.3749 | +1.3 % |
| 0.954 | 0.3059 | 0.3383 | +10.6 % (transonic, see TSD) |
| 0.960 | 0.3076 | 0.3340 | +8.6 % (transonic, see TSD) |

*Other analyses* — Goland wing at ρ = 1.02: DLM 6×12 / 10×20 164.6 / 166.9 m/s vs Murua 2012 165.0, Wang 2006 163.8 / 174.3, Savino 2021 168.2 / 174.2 m/s; legacy deck vs Fortran FLAAP.exe: modes 8–10 agree to a few % in damping and 1–3 % in frequency. *Consistency* — PK / K / p (Roger RFA) 153.3 / 153.3 / 153.1 m/s on the same GAF table.

**TSD transonic flutter.** *Experiment* — same AGARD 445.6 points: TSD (NACA 65A004 nonlinear mean flow + time-linearised unsteady solution) +1.1 / −2.6 / **+0.2** / −3.5 % at M 0.499 / 0.901 / 0.954 / 0.960, recovering the transonic dip that DLM misses. *Consistency* — zero-thickness TSD reproduces DLM. The often-quoted "experimental" 38.17 / 48.35 / 91.54 Hz are Yates' FE values.

**Lifting TSD mean flow** (`validation/tsd_lifting_validation.py`). *Experiment* — none compared yet (transonic NACA 0012 tunnel data carry viscous/wall effects). *Other analysis* — NACA 0012, M 0.8: at section cl 0.35 the TSD upper shock (x/c 0.585) is 0.035 c ahead of the AGARD-AR-211 Euler value 0.62. *Analytic* — CLα vs Helmbold −4 %. *Consistency* — alpha = 0 equals the symmetric solver to 1e-15; small alpha equals linear TSD. Flutter about the lifting mean flow (AGARD 445.6, M 0.954, mean alpha 0 / 1 / 2°) lowers the FSI; only alpha = 0 has a test point.

**Other capabilities** (details per chapter in REPORT.md):
- **Structural modes**: Leissa CFFF plate (shells ≤ 0.6 % at 16×16), Goland beam ±0.15 %, tet10 / hex8 within 0.7–1.4 % of Timoshenko (tet4 locks, +20–28 %); I/O consistency: OpenRadioss, ANSYS (Δf 0.001 %, MAC 1), CalculiX `.frd` (Δf ≤ 0.0003 %, same PK flutter speed), gmsh `.geo`/`.msh`, STL / STEP → tet10.
- **Supersonic**: AGARD 445.6 Mach box +8.5 / +13.9 % at M 1.072 / 1.141 (linear potential theory; piston theory not applicable at these Mach numbers); panel flutter λ within 0.05 % of Dowell; Ackeret exact, Garrick–Rubinow within 0.03 %.
- **VLM / UVLM**: Goland flutter Richardson 167.0 m/s vs Murua's UVLM 165.0 (+1.2 %); Wagner within 1.7 %, Küssner within 1.8 % for s ≥ 4, harmonic pitch → Theodorsen (Richardson ratio 1.0006 / 1.0016), Garrick thrust within 1.6 %; consistency: 3-D impulsive start → steady VLM, prescribed vs free wake.
- **Freeplay LCO**: linear V_f 9.232 vs 9.2 m/s (Vasconcellos 2016); LCO amplitudes not compared with the Conner (1997) rig (different stiffness model). Consistency: describing function vs time marching within 3.8 % in amplitude and 1.5 % in frequency, subcritical hysteresis.
- **ASE / gust / trim** (consistency only): K = 0 closed loop equals the open-loop V_f; Nyquist agrees with closed-loop eigenvalues; time-domain RFA gust equals the frequency-domain response within about 0.1 % at the peaks; 1-g trim lift equals the weight, and the checked manoeuvre reaches n_z = 2.5.
- **Unsteady aero functions** (`tests/test_aero_functions.py`): Theodorsen / Sears within 2·10⁻⁴ of tabulated values.

## Proposed next features

1. **Higher-order supersonic aerodynamics:** ZONA51 / supersonic doublet-point lifting surface (the Mach box and piston theory are in), thickness effects at low supersonic Mach.
2. **g-method / damping-consistent flutter** and matched-point flutter-boundary sweeps over Mach and altitude.
3. **Extending TSD:** viscous (boundary-layer) coupling, Ames modified TSD coefficients, twist, time-marching TSD for shock motion and LCO (nonzero mean incidence and camber are in: `flaap.aero.tsd_lifting`).
4. **More nonlinearities:** cubic / hysteretic hinge stiffness, friction, multiple nonlinear hinges, harmonic balance with higher harmonics.
5. **Dynamic control-law synthesis** (H∞ / LQG) and MIMO singular-value margins on the existing ASE state space.
6. **Correction factors:** W2GJ / Giesing downwash weighting and a body / interference model.
7. **Sensitivities and optimization:** analytic flutter-speed derivatives, mass-balance and stiffness tailoring.
8. **More interfaces:** Nastran OP2 / punch, ZAERO / SOL 145-146 GAF exchange, CFD (Euler / RANS) GAFs.
9. **Uncertainty quantification** of the flutter boundary (Monte Carlo, polynomial chaos).
10. **Test correlation:** GVT import, MAC / frequency matching and model updating against measured modes.
11. **More literature validation:** Vasconcellos (2016) experimental freeplay rig (measured damping matrix, 5th-order pitch stiffness), BACT / Benchmark Active Controls for ASE.
12. **Radioss binary animation reader** (A### files without `anim_to_vtk`) and a Radioss `/EIG` round trip on a real run.
13. **Studio batch studies:** parameter sweeps (Mach, altitude, mass, stiffness) with a per-run TXT/PDF report.

## Legacy Fortran program

The original files are kept unchanged in `legacy_fortran/`:

- `src/*.f90`: sources;
- `build/*.obj`, `*.mod`, `BuildLog.htm`: build products;
- `FLAAP.exe`, the Conyers input deck, `Vg_save.dat` / `Vf_save.dat` / `pk_results_mode*.dat` and the Nastran sample `20.f06`, in the folder itself;
- `cases/`: the Conyers and modified Goland deck folders;
- `old_strange_Versions/`.

Only their location changed. The source of the main program, the DLM and the PK solver was never included in this snapshot (only `.obj` files exist). Those parts were reimplemented from the program flow, the remaining sources and the literature. `python -m flaap legacy <dir>` reads the same input files and writes `Vg_save.dat`, `Vf_save.dat` and `pk_results_mode*.dat` in the original format.
