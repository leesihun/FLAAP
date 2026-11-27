# FLAAP Studio

로컬 브라우저 GUI for FLAAP. It lets you build or edit a configuration, validate it, run it, and look at every `output/*` result. It runs fully offline: standard-library Python plus vanilla JS ES modules, with no CDN and no external packages.

## Start / 실행

```
python studio/start_studio.py                 # http://127.0.0.1:8765 를 브라우저로 엽니다
python studio/start_studio.py --port 9000 --no-browser
```

The server binds to `127.0.0.1` only. Press Ctrl+C to stop it. It never modifies `flaap/` or `configs/`.

- Studio runs write their config to `studio_runs/<name>/config_<name>.txt`.
- Their results go to `output/studio_<name>/`.

## Model tab

- **Config library**: open any file under `configs/`, or start a new model with **New**.
- **Planform editor** (centre). Edit surfaces, controls and stations:
  - Use **Select** (V) to pick and drag the corner handles, the hinge/eta handles of controls, or whole objects.
  - Use **Surface** (S) to add a surface with four clicks.
  - Use **Station** (M) to add a station.
  - **F** fits the view and **Del** deletes the selection.
  - The view shows DLM boxes, the mirror image, the plate FE grid and the clamp line.
  - Preset structures (goland, agard445, …) appear as grey "ghost" surfaces. You can convert them to editable `surfaceN` keys.
- **Properties** (right): numeric properties of the selection, plus the model/mode route form generated from `/api/models`.
- **Bottom panel**:
  - *Config text* stays in two-way sync with the editor.
  - *Diagnostics* shows the preflight report. Click a row to jump to the line.
  - *Run log* streams the log of a run.
  - **Validate** / **Run** / **Cancel**. **Open results** jumps to the results of the run.

## Results tab

Pick a run on the left. A view appears only when the run has the files it needs.

| View | Source files |
|---|---|
| Summary | `results.txt`: stat tiles, flattened summary, artifact list |
| Modes 3D | `modes.txt` (+ `aero_mesh.txt`, `aero_h`): animated mode shapes, orbit/zoom/pan, amplitude, mirror |
| ΔCp | `cp.txt` + `aero_mesh.txt`: Re/Im/\|ΔCp\| per k and per mode, diverging colour bar |
| Flutter | every `flutter_*.txt` (globbed, not hard-coded): method selector, linked V-g / V-f plots, critical points of all methods |
| Time | `timehist_*.txt`: small multiples grouped by name prefix, series toggles |
| Curves | `curves_*.txt`: log x/y toggles, grouped small multiples |
| TSD | `tsd_steady.txt`: chordwise Cp per station, planform Cp map, convergence |
| Figures | `figures/*.png` gallery (click to enlarge) |
| Files | `*.txt` / `*.csv` text viewer |

The page URL tracks the open view. You can link a view directly, e.g. `http://127.0.0.1:8765/#results/flutter_goland_all/flutter`.

## API (JSON)

**GET**
- `/api/health`
- `/api/models`
- `/api/presets`
- `/api/configs`
- `/api/config?path=configs/...`
- `/api/job?id=&since=`
- `/api/jobs`
- `/api/outputs`
- `/api/file?path=output/...`
  - Only files under `output/` are served, and path traversal is rejected.

**POST**
- `/api/validate` `{text, source}`
- `/api/run` `{text, name, source}`
- `/api/cancel` `{id}`

## Files

```
studio/start_studio.py      server (stdlib only)
studio/web/index.html       layout
studio/web/css/studio.css   AI-CAE4ALL palette
studio/web/js/              main, api, util, config (parse/generate), editor (planform),
                            model_tab, results_tab, plot (SVG charts), view3d (canvas 3-D), colormap
```
