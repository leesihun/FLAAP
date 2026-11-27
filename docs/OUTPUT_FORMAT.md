# FLAAP output format

Every launcher run writes one folder (`output_dir`, default `output/<config stem>`).
Results are **plain-text files in the same `key value` grammar as the input configs**
(AI-CAE4ALL style). You can read them in any editor, diff them and grep them. Coordinates
are in metres; x points downstream, y to the right (span) and z up.

FLAAP Studio, the tests and the white paper all read these files with one parser,
`flaap.io.results_txt.read_result(path)`. It returns the same dictionary as the former JSON
contract. To get JSON as well, add `output_json true` to the config. To convert one file
on the command line, run `python -m flaap.io.results_txt output/<run>/flutter_pk.txt`. It prints JSON.

## Text format (`flaap-txt 1`)

```text
% ==========================================================================
% FLAAP result file  flutter_pk.txt
% Goland wing - PK flutter
% written 2026-11-26 14:02:11  (format flaap-txt 1: key value lines + 'table ... end' blocks)
% ==========================================================================
kind          flutter
title         "Goland wing - PK flutter"
method        pk
density       1.225
names         "mode 1", "mode 2", "mode 3", "mode 4"

% branch = 0-based branch index (its V-g table is branch.<branch+1>)
table points records
      kind  branch  velocity  freq_hz        k  density
   flutter       1  153.2871  11.2702  0.4253    1.225
end

% branch 2: mode 2
table branch.2
        V          g   freq_hz         k
       60  -0.01523     7.712    0.2563
      ...
end
```

### Rules

* **Lines and comments.** Each line holds one `key value` pair. `%` starts a comment line and `#` starts an inline comment.
* **Nesting.** A dotted key stands for a nested value. `summary.flutter_pk.velocity 153.29` reads back as
  `{"summary": {"flutter_pk": {"velocity": 153.29}}}`.
* **Values.**
  * Comma-separated values form a list. `[x]` is a one-element list and `[]` is an empty list.
  * `true` and `false` are booleans.
  * `none` marks a missing value.
  * A string is put in `"..."` quotes in two cases: when it would otherwise read as a number, boolean or list, and when it contains `,`, `#` or `"`.
* **Tables.** A `table <name> [records|matrix]` … `end` block holds an array.
  * The first line names the columns. The rows below are separated by whitespace.
  * In a table, `nan` marks a missing value, and a string with spaces is quoted.
  * A `records` table reads back as a list of dicts. A `matrix` table reads back as a list of rows. A plain table holds column data, which the reader for that file kind interprets.
* **Numbered tables.** `shape.N`, `aero_h.N`, `cp.N` and `branch.N` give one table per mode, reduced frequency or branch. N starts at 1.
* **Precision.** Scalars are written at full precision, so they round-trip exactly. Data tables use 5–8 significant digits.
* **Node numbers in `modes.txt`.** Nodes are numbered from 1. The connectivity tables use the same numbers, and the reader converts them back to 0-based indices.

## Always written

| file | content |
|---|---|
| `config_resolved.txt` | the full resolved config (user keys and defaults), which runs as is |
| `results.txt` | the run record (see below) |
| `log.txt` | solver log |
| `figures/*.png` | matplotlib figures (unless `plots false`) |

### `results.txt` (kind `results`)

```text
kind          results
title         "Goland wing - PK flutter"
model         flutter
mode          pk
status        ok                 % or error; then `error` holds the traceback
elapsed_s     4.6
created       "2026-11-26 11:23:18"
config        C:/.../config_goland_pk.txt
summary.frequencies_hz      7.67, 15.24, 38.85, 55.36
summary.flutter_pk.kind     flutter
summary.flutter_pk.velocity 153.29

table artifacts records
           file   kind  title
      modes.txt  modes  "Mode shapes and frequencies"
  ...
end
```

`artifacts[].kind` is one of `config, modes, aero_mesh, cp, gaf, flutter, timehist,
curves, tsd_steady, table, text, figure`. Viewers choose how to display a file from this kind.

## Structural and aerodynamic model

### `modes.txt` (kind `modes`)

| item | content |
|---|---|
| `name`, `structure`, `n_modes`, `n_nodes` | model name, structure source and sizes |
| `table modes` | `mode label freq_hz omega_rad_s gen_mass` |
| `table nodes` | `node x y z`: structural node coordinates |
| `table edges` / `quads` / `tris` | connectivity in 1-based node numbers, whichever exists |
| `table shape.N` | `node ux uy uz`: translations of mode N, normalised to max \|u\| = 1 |
| `table aero_h.N` | `box h1 h2 h3 h4`: normal displacement of the 4 corners of each aero box, from the spline |

### `aero_mesh.txt` (kind `aero_mesh`)

The key-value lines are `surfaces`, `symmetry`, `bref`, `mach`, `controls` and `n_boxes`.
Then `table boxes` has the columns `box surface area xc yc zc x1 y1 z1 … x4 y4 z4 ctrl:<name>`:

* The corners are ordered root-LE, tip-LE, tip-TE, root-TE.
* `surface` is the 0-based surface index.
* `(xc, yc, zc)` is the 3/4-chord collocation point.
* `ctrl:<name>` is `true` when the box is part of that control surface.

## Unsteady aerodynamics

### `cp.txt` (kind `cp`)

The key-value lines are `mach`, `bref`, `modes`, `k` (a list) and `n_boxes`. `modes` holds nm + 1 labels; the last one is `gust` when a gust column is present.
Then there is one `table cp.N` for each reduced frequency k_N, with the columns `box re:<mode> im:<mode> …`.
Each value is the lifting-pressure jump ΔCp per unit generalized coordinate. For the gust column it is per unit gust angle.

### `gaf_table.npz` / `gaf_table.txt` (kind `gaf`)

These files hold the complex generalized aerodynamic forces, normalised by the dynamic pressure (F = q Q(k) η):
`k` (nk), `Q` (nk, nm, nm) and the gust column `Qg` (nk, nm). They also store `mach`, `bref` and `names`.
The table has the same layout for every `aero_method` (`dlm`, `tsd`, `uvlm`, `machbox`, `piston`).

### `vlm_steady.txt` (kind `vlm_steady`)

Steady VLM polar (`model uvlm`, `mode steady`): `alpha_deg`, `CL`, `CDi` (Trefftz plane), `Cm` (about `x_ref`,
nose up, divided by `c_ref`), the scalars `CL_alpha` [1/rad], `span_efficiency`, `aspect_ratio`, `s_ref`,
`c_ref`, `x_ref`, `x_ac`, `mach`, and the span loads `y_strip`, `cl_c` (nalpha x nstrip, local cl c / c_ref)
and `cl_local` (local section cl).

### `uvlm_wake.txt` (kind `uvlm_wake`)

Final wake of a time-domain UVLM run: `shape` (rows, columns, 3 per surface) and `surface<N>`, the wake node
coordinates (row-major, x y z per node). `free_wake` tells whether the wake was force-free.

### `tsd_steady.txt` (kind `tsd_steady`)

The key-value lines are `mach`, `cp_star` (the sonic Cp), `converged`, `residuals` (a list), `ny` and `nx`.
Then `table grid` has the columns `j i y x cp`. Each row is the upper-surface Cp at chordwise point i of span station j.

## Flutter

### `flutter_<method>.txt` (kind `flutter`), method = `pk | k | p`

| item | content |
|---|---|
| `method`, `method_title` | short key and full title |
| `density`, `bref`, `g_struct`, `names`, `n_velocities` | flight and reference data, branch names |
| `v_per_branch` | `true` when each branch has its own velocity grid (K method) |
| `table points records` | crossings: `kind (flutter/divergence) branch velocity freq_hz k density` |
| `table branch.N` | `V g freq_hz k` for branch N: velocity [m/s], damping g (= 2ζ), frequency [Hz], reduced frequency |

The folder `flutter_<method>/` holds the legacy Fortran-style text files (`Vg_save.dat`, `Vf_save.dat`, `pk_results_mode<i>.dat`).
In `results.txt`, `summary.flutter_<method>.*` is the critical (lowest) point.

## Time domain and curves

### `timehist_<name>.txt` (kind `timehist`)

The key-value lines are `meta.*`, `units.*` and `series`. `meta` holds case data such as the gust length H, the velocity or the control channel.
Then `table timehist` has the columns `t <series 1> <series 2> …`.
This file comes from the gust route (critical 1-cos response), the maneuver route, the freeplay route (LCO time history) and the ASE route (closed loop).

### `curves_<name>.txt` (kind `curves`)

The key-value lines are `xlabel`, `ylabel`, `logy`, `meta.*` and `curves`. Then `table curves` has the columns `x <curve 1> …`.
These files hold, for example, the gust envelope vs. H, turbulence PSDs, the describing-function LCO amplitude vs.
velocity, and ASE loop-gain and margin curves.

## Route summaries (`results.txt` → `summary`)

| route | summary keys |
|---|---|
| modal | `frequencies_hz`, `nmodes`, `n_boxes` |
| aero | `nk`, `mach`, `bref`, `Q0_diag` |
| flutter | `density`, `flutter_<method>` |
| gust discrete | `critical_H`, `peak` per output |
| gust continuous | `A_bar`, `N0` per output |
| trim static | `alpha_deg`, `<control>_deg`, `loads` |
| trim divergence | `q_divergence`, `V_divergence` |
| maneuver | `short_period_wn`, `short_period_zeta`, `peak_control_deg`, `peak_nz`, `peak_loads` |
| tsd steady | `converged`, `iterations`, `n_supersonic`, `cp_min`, `cp_star` |
| freeplay describing_function | `V_f_linear` (K = K_h), `V_f0` (K = 0), `V_fold` (lowest stable LCO speed), `time_marching` rows `{V, kind, A_tm_deg (peak), A1_tm_deg (1st harmonic), f_tm_hz, A_df_deg, f_df_hz}` |
| freeplay time | `lco` `{V, kind (decay/lco/growing), amplitude_deg, amplitude_h1_deg, mean_deg, freq_hz, growth, switches}` |
| freeplay bifurcation | `bifurcation` `{first_lco_up, kinds_up, kinds_down}` (amplitudes in `curves_bifurcation.txt`) |
| ase flutter_suppression | `flutter_open`, `flutter_closed`, `V_f_open`, `V_f_closed`, `increase` (fraction) |
| ase margins | `margins` `{V, P, N_ccw, nyquist_stable, closed_loop_stable, gm_lower_db, gm_upper_db, pm_deg, w_pm}`, `margins_vs_V` |
| ase gla | `gla` `{V, Uds, rows: [{H, Mx_open, Mx_closed, Mx_reduction, Fz_open, Fz_closed, flap_peak_deg}]}`, `V_f_open`, `V_f_closed` |
| uvlm (all) | `s_ref`, `c_ref`, `x_ref`, `n_panels` |
| uvlm steady | `CL_alpha_per_rad`, `CL_alpha_per_deg`, `span_efficiency`, `aspect_ratio`, `x_ac`, `CL_max_alpha`, `CDi_max_alpha` |
| uvlm impulsive / harmonic / gust | `velocity`, `density`, `dt`, `nsteps`, `CL_final`, `CD_final`, `Cm_final` |
| uvlm impulsive | `CL_steady_vlm`, `CDi_steady_vlm`, `CL_final_over_steady` |
| uvlm harmonic | `CL_mean`, `CL_1st_harmonic_amplitude`, `CL_1st_harmonic_phase_deg` (vs sin ωt), `CD_mean_last_cycle`, `CT_mean_last_cycle` (= -CD, thrust), `Cm_1st_harmonic_amplitude` |
| uvlm gust | `dCL_gust_peak`, `t_peak`, `dL_gust_peak`, `dMy_gust_peak`, `CL_at_peak` (increments over a gust-free run when `alpha` ≠ 0) |
| uvlm flutter | `density`, `flutter_pk` (PK on linearised-UVLM GAFs), `flutter_uvlm_state_space` `{velocity, freq_hz, mode}` |

All freeplay / ase summaries also carry `plant`, `rho` and `defaults` (plant-specific
defaults that were applied because the key was not in the config).

### Files written by the freeplay, ase and uvlm routes

| route / mode | files |
|---|---|
| freeplay describing_function | `curves_lco_df.txt` (groups `A/delta0:...`, `f [Hz]:...` per branch, `nan` off-branch), `curves_lco_tm_vs_df.txt`, `timehist_lco.txt`, `figures/lco_bifurcation.png` |
| freeplay time | `timehist_lco.txt` (hinge, pitch, plunge or modal coordinates), `figures/phase_plane.png` |
| freeplay bifurcation | `curves_bifurcation.txt` (groups `A [deg]:sweep up/down`, `f [Hz]:sweep up/down`) |
| ase flutter_suppression | `flutter_open.txt`, `flutter_closed.txt` (+ legacy folders), `figures/ase_vg_rootlocus.png` |
| ase margins | `curves_bode.txt` (groups `|L| [dB]:loop`, `phase [deg]:loop`; `meta` = the margins), `figures/nyquist_bode.png` |
| ase gla | `timehist_gla.txt` (open / closed root Mx, Fz and flap angle for each H) |
| uvlm steady | `vlm_steady.txt`, `curves_vlm_lift.txt` (CL, Cm vs α), `curves_vlm_polar.txt` (CDi vs CL + elliptic CL²/πA), `curves_vlm_spanload.txt`, `aero_mesh.txt` |
| uvlm impulsive / harmonic / gust | `timehist_uvlm.txt` (CL, CD, Cm, L [N], D [N], My [N m], plus `pitch_deg`, `plunge`, `gust_at_LE` and, for gusts at α ≠ 0, `dCL_gust`, `dL_gust`, `dMy_gust`), `curves_uvlm_spanload.txt`, `uvlm_wake.txt`, `figures/uvlm_wake.png` |
| uvlm flutter | `modes.txt`, `gaf_table.*`, `flutter_pk.txt`, `curves_uvlm_ss_damping.txt`, `curves_uvlm_ss_frequency.txt` |

Curve labels of the form `group:name` are plotted together in Studio, one axis per group.
