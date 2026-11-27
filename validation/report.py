"""Run the validation cases and write ``validation/REPORT.md`` (+ figures).

The report is organised by the *capability* that is validated (e.g. "linear
DLM flutter", "TSD transonic flutter"), not by test case.  Inside every
chapter the references are given in priority order:

1. experiment,
2. other independent analyses / codes,
3. analytic solutions (only where nothing better exists),

followed by self-consistency checks (code-vs-code limits, solver-vs-solver,
I/O round trips), which are reported separately and are not counted as
validation.  Missing reference levels are stated explicitly.

Run:  python validation/report.py                      (all chapters, ~40 min)
      python validation/report.py uvlm supersonic      (only those, spliced into REPORT.md)
      python validation/report.py --rerun uvlm         (also re-run the script behind a log-based chapter)
Chapters: structures dlm strip tsd_flutter tsd_lifting supersonic uvlm freeplay ase gust trim gaps
"""
from __future__ import annotations

import os
import re
import runpy
import subprocess
import sys
import time
from functools import lru_cache

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from flaap import gust as G  # noqa: E402
from flaap.aero.gaf import AeroModel, AeroModes  # noqa: E402
from flaap.aero.rfa import fit_roger  # noqa: E402
from flaap.flutter import k_method, pk_method  # noqa: E402
from flaap.geometry import build_mesh  # noqa: E402
from flaap.structures.spline import beam_spline  # noqa: E402
from flaap.timedomain import AeroelasticSystem, modal_system  # noqa: E402

FIG = os.path.join(HERE, "figures")
RERUN = False                      # set by --rerun: re-run the scripts behind log-based chapters

EXP, OTH, ANA, CON = "1순위 실험", "2순위 타 해석", "3순위 해석해", "자기일관성 (검증 근거 아님)"


# ======================================================================= helpers
def _pct(a, b):
    return 100.0 * (a / b - 1.0)


def _head(md, num, title, target, refs, verdict, summary):
    """Chapter heading: what is validated, references by priority, verdict.

    ``refs`` = [(level, text)], ``summary`` = (exp, other, analytic, result) for
    the overview table (stored as an HTML comment so that spliced chapters keep it)."""
    md.append(f"## {num}. {title}\n")
    md.append(f"**검증 대상:** {target}\n")
    md.append("**비교 기준 (우선순위 순):**\n")
    for lvl, txt in refs:
        md.append(f"- **{lvl}** — {txt}")
    md.append(f"\n**판정:** {verdict}\n")
    md.append("<!-- summary | " + " | ".join([title] + [s.replace("|", "/") for s in summary]) + " -->\n")


def _sub(md, num, level, title):
    md.append(f"### {num} [{level}] {title}\n")


def _rows_table(md, rows, ref_label="ref", extras=("min_mac", "nmodes", "nodes", "f_hz", "f_ref_hz")):
    """Rows of the structures script (dicts) -> markdown table."""
    md.append(f"| case | quantity | FLAAP | {ref_label} | err % | notes |\n|---|---|---|---|---|---|")
    for r in rows:
        ref = f"{r['reference']:.5g}" if "reference" in r else ""
        err = f"{r['err_pct']:+.3f}" if "err_pct" in r else ""
        notes = ", ".join(f"{k} {r[k]:.4g}" if isinstance(r[k], float) else f"{k} {r[k]}"
                          for k in extras if k in r)
        md.append(f"| {r['case']} | {r['quantity']} | {r['computed']:.5g} | {ref} | {err} | {notes} |")
    md.append("")


def _log(script, txt, rerun=None):
    """Text output of a validation script, cached in validation/<txt>."""
    path = os.path.join(HERE, txt)
    if (RERUN if rerun is None else rerun) or not os.path.exists(path):
        env = dict(os.environ, PYTHONPATH=ROOT, MPLBACKEND="Agg")
        p = subprocess.run([sys.executable, os.path.join(HERE, script)], cwd=ROOT, env=env,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if p.returncode:
            raise RuntimeError(f"{script} failed:\n{p.stderr[-3000:]}")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(p.stdout)
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def _chunks(text):
    """Split a log into its numbered parts ('1-3.', '4.', ...) -> {key: text}."""
    out, key = {"pre": []}, "pre"
    for line in text.splitlines():
        m = re.match(r"^(\d+(?:-\d+)?)\. ", line)
        if m:
            key = m.group(1)
            out[key] = []
        out[key].append(line)
    return {k: "\n".join(v).rstrip() for k, v in out.items()}


def _code(md, text):
    md.append("```\n" + text.rstrip() + "\n```\n")


def _find(pattern, text, flags=re.M):
    return [tuple(float(g) for g in m.groups()) for m in re.finditer(pattern, text, flags)]


# ============================================================ shared Goland data
@lru_cache(None)
def _goland():
    from flaap.models.goland import CHORD, REFERENCE, goland_modal, goland_surface
    mm = goland_modal(4)
    kv = np.r_[0, np.geomspace(0.01, 2, 19)]

    def dlm_table(nc, ns):
        mesh = build_mesh([goland_surface(nc, ns)], "symmetric")
        modes = AeroModes.from_spline(beam_spline(mesh, mm.nodes), mm.dof_vectors())
        return AeroModel(mesh, modes, 0.0, CHORD / 2).table(kv)
    return dict(mm=mm, tab=dlm_table(6, 12), dlm_table=dlm_table, CHORD=CHORD, REFERENCE=REFERENCE)


@lru_cache(None)
def _structures_run():
    S = runpy.run_path(os.path.join(HERE, "structures_validation.py"))
    return S["run"](verbose=False, plot=True, flutter=True), S


# ===================================================================== chapters
def structures(md):
    from flaap.models.agard445 import EXPERIMENT_FREQ_HZ, agard445_fe
    out, S = _structures_run()
    g = _goland()
    mm, ref = g["mm"], g["REFERENCE"]
    ffe = agard445_fe(n_chord=16, n_span=16).modes(4).freq_hz
    fexp = np.asarray(EXPERIMENT_FREQ_HZ, float)
    yates = np.array([9.60, 38.17, 48.35, 91.54])
    e_fl = np.round(_pct(ffe, fexp), 1) + 0.0          # +0.0 avoids printing "-0.0"
    e_ya = _pct(yates, fexp)
    lei = [r for r in out["shell_leissa"] if r["case"] in ("quad4 16x16", "tri3 16x16")]
    e_lei = max(abs(r["err_pct"]) for r in lei)
    e_gol = np.abs(_pct(mm.omega, np.asarray(ref["omega_rad_s"]))).max()
    sol, sol_case = {}, {}                              # finest (last-listed) mesh of each element type
    for r in out["solid_beam"]:
        el = r["case"].split()[0]
        if sol_case.get(el) != r["case"]:
            sol_case[el], sol[el] = r["case"], []
        sol[el].append(abs(r["err_pct"]))
    _head(md, 1, "구조 고유진동 해석 검증 (판·쉘·보·고체 FE)",
          "`flaap.structures`의 MITC4 판, DKT/Allman 쉘, Euler-Bernoulli 보, hex8/tet10/tet4 고체 요소로 "
          "구한 고유진동수·모드형상 (플러터 해석의 구조 입력).",
          [(EXP, "AGARD 445.6 weakened model 3 측정 고유진동수 (Yates 1987, NASA TM-100492)."),
           (OTH, "Yates (1987) FE 모델 진동수; Leissa (1973) CFFF 판 고유치 (Ritz 해석)."),
           (ANA, "Goland (1945) 보 정확해; Timoshenko 보 이론 (고체 요소 보)."),
           (CON, "강체 모드, OpenRadioss/ANSYS/CalculiX/gmsh/STL/STEP 입력 경로 왕복.")],
          f"실험 대비 MITC4 판 진동수 오차 {e_fl[0]:+.1f} / {e_fl[1]:+.1f} / {e_fl[2]:+.1f} / {e_fl[3]:+.1f} % "
          f"(1차 굽힘/1차 비틀림/2차 굽힘/2차 비틀림; 같은 실험 대비 Yates FE 모델 {e_ya[1]:+.1f} / {e_ya[2]:+.1f} / "
          f"{e_ya[3]:+.1f} %). Leissa 대비 16x16 쉘 최대 {e_lei:.2f} %, Goland 정확해 대비 보 최대 {e_gol:.2f} %, "
          f"Timoshenko 대비 {sol_case['hex8']} 최대 {max(sol['hex8']):.1f} %, {sol_case['tet10']} 최대 "
          f"{max(sol['tet10']):.1f} %, {sol_case['tet4']} {max(sol['tet4']):.0f} % (선형 사면체의 잠김: tet10 사용 권장).",
          ("AGARD 445.6 진동수", "Leissa, Yates FE", "Goland 보, Timoshenko",
           f"실험 대비 {e_fl.min():+.1f}~{e_fl.max():+.1f} %; 해석해 대비 <= "
           f"{max(e_gol, max(sol['hex8']), max(sol['tet10'])):.1f} % (tet4 제외)"))

    _sub(md, "1.1", EXP, "AGARD 445.6 판 날개 고유진동수 (MITC4 판 16x16)")
    md.append("재료: Kolonay (1996) 등가 직교이방성 물성 (진동수 맞춤 보정 없음).\n")
    md.append("| 진동수 [Hz] | 1차 굽힘 | 1차 비틀림 | 2차 굽힘 | 2차 비틀림 |\n|---|---|---|---|---|")
    md.append("| **실험 (Yates 1987)** | " + " | ".join(f"{f:.2f}" for f in fexp) + " |")
    md.append("| FLAAP MITC4 16x16 | " + " | ".join(f"{f:.2f}" for f in ffe) + " |")
    md.append("| FLAAP 오차 | " + " | ".join(f"{e:+.1f} %" for e in e_fl) + " |")
    md.append("| (타 해석) Yates FE 모델 | " + " | ".join(f"{f:.2f}" for f in yates) + " |")
    md.append("| Yates FE 오차 | " + " | ".join(f"{e:+.1f} %" for e in e_ya) + " |\n")

    _sub(md, "1.2", OTH, "직사각형 외팔 판 (CFFF, a/h = 100) 고유치: Leissa (1973)")
    md.append("lambda = omega a^2 sqrt(rho h / D).\n")
    md.append(f"| 요소 / 격자 | " + " | ".join(f"lambda{i + 1}" for i in range(5)) + " |\n|---|---|---|---|---|---|")
    leissa = [r["reference"] for r in out["shell_leissa"][:5]]
    md.append("| **Leissa (1973)** | " + " | ".join(f"{v:.3f}" for v in leissa) + " |")
    cases = []
    for r in out["shell_leissa"]:
        if r["case"] not in cases:
            cases.append(r["case"])
    for c in cases:
        rr = [r for r in out["shell_leissa"] if r["case"] == c]
        md.append(f"| {c} | " + " | ".join(f"{r['computed']:.3f} ({r['err_pct']:+.2f} %)" for r in rr) + " |")
    md.append("\n![shell convergence](figures/structures_shell_convergence.png)\n")

    _sub(md, "1.3", ANA, "보: Goland (1945) 정확해, Timoshenko 보 이론 (고체 요소)")
    md.append("| Goland 날개 모드 | FLAAP 보 FE [rad/s] | 정확해 (Goland 1945) [rad/s] | 오차 |\n|---|---|---|---|")
    for i, (w, we) in enumerate(zip(mm.omega, ref["omega_rad_s"])):
        md.append(f"| {i + 1} | {w:.2f} | {we:.2f} | {_pct(w, we):+.2f} % |")
    md.append("\n정확해 수치: Savino et al., *Energies* 14 (2021) Table 2.\n")
    md.append("강철 외팔 보 (고체 요소) vs Timoshenko 보 이론:\n")
    _rows_table(md, out["solid_beam"], "Timoshenko", extras=("euler_bernoulli", "dofs"))
    md.append("![beam modes](figures/structures_beam_modes.png)\n")

    _sub(md, "1.4", CON, "강체 모드와 입력 경로 (외부 FE 결과·형상 파일)")
    md.append("기준값은 위 1.2/1.3과 같은 Leissa/Timoshenko 값 또는 원본 FE 결과 자체 (읽기 왕복).\n")
    for key, title in (("rigid_modes", "자유-자유 강체 모드 (0 Hz 이어야 함)"),
                       ("openradioss", "OpenRadioss 입력 덱 -> FLAAP 고유해석"),
                       ("ansys", "ANSYS MAPDL/Workbench 결과 읽기 (원본 대비)"),
                       ("calculix", "CalculiX .frd/.dat 읽기 (원본 대비) 및 플러터"),
                       ("gmsh", "gmsh 메셔 (.geo/.msh)"),
                       ("geometry", "STL/STEP 형상 -> 자동 격자 -> 고유해석"),
                       ("wing_flutter", "STL 날개 형상 -> tet10 모드 -> PK 플러터 (시연, 기준값 없음)")):
        if key in out:
            md.append(f"**{title}**\n")
            _rows_table(md, out[key], extras=("first_elastic", "min_mac", "nmodes", "nodes", "f_hz",
                                              "f_ref_hz", "freq_hz"))
    md.append("![plate modes](figures/structures_plate_modes.png)\n")


def dlm(md):
    from flaap.io.legacy import run_legacy
    g = _goland()
    mm, tab, ref = g["mm"], g["tab"], g["REFERENCE"]
    from flaap.models.goland import DLM_REFERENCES_RHO102
    mod = runpy.run_path(os.path.join(HERE, "agard445_flutter.py"))
    rm = mod["run"](exp_freq=True, verbose=False)
    rf = mod["run"](exp_freq=False, verbose=False)
    sub = rm[:, 0] < 0.93
    trans = (rm[:, 0] > 0.93) & (rm[:, 0] < 1.0)
    em = _pct(rm[:, 2], rm[:, 1])
    rho1 = ref["dlm_density"]
    V1 = np.linspace(80, 300, 221)
    c3 = {}
    for nc, ns in ((6, 12), (10, 20)):
        t = tab if (nc, ns) == (6, 12) else g["dlm_table"](nc, ns)
        c3[(nc, ns)] = pk_method(mm.M, mm.K, t, rho1, V1).critical()
    refv = [r[2] for r in DLM_REFERENCES_RHO102]
    e3 = _pct(c3[(10, 20)].velocity, np.mean(refv))
    _head(md, 2, "선형 DLM 플러터 검증 (아음속, DLM + PK)",
          "이중격자법(DLM) 비정상 공력 + 스플라인 + PK 플러터 해석의 플러터 속도·진동수 "
          "(`aero_method dlm`, `mode pk`).",
          [(EXP, "AGARD 445.6 weakened model 3 플러터 경계 (Yates 1987), M 0.499-0.960."),
           (OTH, "Goland 날개 3-D 결과 (Murua 2012 UVLM, Wang 2006 UVLM/ZAERO, Savino 2021 DUST); "
                 "원본 Fortran FLAAP.exe (같은 입력 덱)."),
           (ANA, "없음 — 3-D 비정상 양력면의 해석해는 없음 (2-D 극한은 3장 스트립 이론)."),
           (CON, "같은 GAF 표에서 PK / K / p(Roger RFA) 해법 비교, DLM 격자 수렴.")],
          f"아음속 선형 영역 (M <= 0.9) 실험 대비 FSI 오차 {em[sub].min():+.1f}~{em[sub].max():+.1f} % "
          f"(측정 진동수 사용). 천음속 M 0.954/0.960에서는 {em[trans][0]:+.1f} / {em[trans][1]:+.1f} %로 "
          f"천음속 딥을 놓침 — 선형 이론의 한계이며 4장 TSD가 보완. Goland 3-D: 문헌 범위 "
          f"{min(refv):.1f}-{max(refv):.1f} m/s 안 ({c3[(10, 20)].velocity:.1f} m/s, 평균 대비 {e3:+.1f} %).",
          ("AGARD 445.6 (M<=0.9)", "Goland 3-D 문헌 5건, Fortran", "없음",
           f"실험 대비 {em[sub].min():+.1f}~{em[sub].max():+.1f} %; 타 해석 평균 대비 {e3:+.1f} %"))

    _sub(md, "2.1", EXP, "AGARD 445.6 플러터 경계 (Yates 1987)")
    md.append("FSI, wf/wa는 b = 0.2794 m, omega_alpha = 239.3 rad/s 기준. 시험 질량비를 맞춤 "
              "(rho_eq = rho * m_model / m_test). 구조: 1.1의 MITC4 판 모드형상.\n")
    md.append("| Mach | FSI 실험 | FLAAP (측정 진동수) | 오차 % | FLAAP (FE 진동수) | 오차 % | wf/wa 실험 | wf/wa FLAAP |\n"
              "|---|---|---|---|---|---|---|---|")
    for a, b in zip(rm, rf):
        if a[0] > 1:
            continue
        md.append(f"| {a[0]:.3f} | {a[1]:.4f} | {a[2]:.4f} | {_pct(a[2], a[1]):+.1f} | {b[2]:.4f} | "
                  f"{_pct(b[2], b[1]):+.1f} | {a[3]:.4f} | {a[4]:.4f} |")
    md.append("\nM > 1 점 (1.072, 1.141)은 DLM 적용 범위 밖 — 6장 초음속 Mach box에서 검증.\n")
    fig, ax = plt.subplots(figsize=(6, 4))
    k = rm[:, 0] < 1
    ax.plot(rm[k, 0], rm[k, 1], "ko-", label="experiment (Yates 1987)")
    ax.plot(rm[k, 0], rm[k, 2], "bs--", label="FLAAP DLM + PK (measured freq.)")
    ax.plot(rf[k, 0], rf[k, 2], "c^:", label="FLAAP DLM + PK (FE freq.)")
    ax.set_xlabel("Mach"); ax.set_ylabel("FSI"); ax.grid(alpha=0.3); ax.legend()
    fig.savefig(os.path.join(FIG, "agard_fsi.png"), dpi=110)
    plt.close(fig)
    md.append("![AGARD FSI](figures/agard_fsi.png)\n")

    _sub(md, "2.2", OTH, f"Goland 날개 3-D 플러터 (rho = {rho1} kg/m^3, 문헌 양력면/패널 해석)")
    md.append("| 출처 | 방법 | V_f [m/s] | f_f [Hz] |\n|---|---|---|---|")
    for (nc, ns), c in c3.items():
        md.append(f"| **FLAAP** | DLM {nc}x{ns}, PK | {c.velocity:.1f} | {c.omega / (2 * np.pi):.2f} |")
    for src, meth, vf, ff in DLM_REFERENCES_RHO102:
        md.append(f"| {src} | {meth} | {vf:.1f} | {'' if ff is None else f'{ff:.2f}'} |")
    md.append(f"\n문헌 범위 {min(refv):.1f}-{max(refv):.1f} m/s (평균 {np.mean(refv):.1f}); FLAAP "
              f"{c3[(6, 12)].velocity:.1f} / {c3[(10, 20)].velocity:.1f} m/s.\n")

    legacy_dir = os.path.join(ROOT, "legacy_fortran")
    if os.path.exists(os.path.join(legacy_dir, "pk_results_mode7.dat")):
        res, case, _ = run_legacy(legacy_dir, velocities=np.linspace(5, 300, 296),
                                  kvals=np.r_[0.0, np.geomspace(0.01, 3.0, 20)])
        _sub(md, "2.3", OTH, "원본 Fortran FLAAP.exe (Conyers 판 날개, 같은 입력 덱)")
        md.append("Python 레거시 스플라인 + DLM + PK vs Fortran `pk_results_mode*.dat` (진동수 Hz).\n")
        md.append("| mode | V [m/s] | g Fortran | g FLAAP-py | f Fortran | f FLAAP-py |\n|---|---|---|---|---|---|")
        fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
        for m in (7, 8, 9, 10):
            rf_ = np.loadtxt(os.path.join(legacy_dir, f"pk_results_mode{m}.dat"))
            b = m - 1
            for V in (50.0, 150.0, 250.0):
                i = np.argmin(np.abs(rf_[:, 0] - V))
                j = np.argmin(np.abs(res.V[:, b] - V))
                md.append(f"| {m} | {V:.0f} | {rf_[i, 1]:+.4f} | {res.g[j, b]:+.4f} | {rf_[i, 2]:.2f} | "
                          f"{res.freq[j, b]:.2f} |")
            ax[0].plot(rf_[:, 0], rf_[:, 1], "k-", lw=0.8)
            ax[0].plot(res.V[:, b], res.g[:, b], "--", label=f"mode {m}")
            ax[1].plot(rf_[:, 0], rf_[:, 2], "k-", lw=0.8)
            ax[1].plot(res.V[:, b], res.freq[:, b], "--")
        ax[0].legend(); ax[0].set_ylabel("g"); ax[1].set_ylabel("f [Hz]"); ax[1].set_xlabel("V [m/s]")
        ax[0].set_title("solid: Fortran FLAAP.exe, dashed: Python")
        fig.savefig(os.path.join(FIG, "legacy_vs_fortran.png"), dpi=110)
        plt.close(fig)
        md.append("\n![legacy](figures/legacy_vs_fortran.png)\n")
        md.append("모드 8-10은 감쇠 수 %, 진동수 1-3 % 이내 일치. 약 120 m/s 이상에서 모드 7은 모드 8과 "
                  "합체 영역에 들어가 (Fortran 파일 자체도 250 m/s에서 g = -9.6) 두 코드의 모드 추적 "
                  "라벨이 달라짐. Fortran 출력에는 모드 7-10만 있어, Python이 예측한 저속 모드 2 불안정 "
                  f"({res.critical()}) 은 교차 확인할 수 없음.\n")

    _sub(md, "2.4", CON, "플러터 해법 간 비교와 격자 수렴 (Goland, rho = 1.225)")
    V = np.linspace(60, 260, 101)
    rho0 = ref["density"]
    md.append("| 해법 (같은 DLM 6x12 GAF 표) | V_f [m/s] | w_f [rad/s] |\n|---|---|---|")
    for name, r in (("PK", pk_method(mm.M, mm.K, tab, rho0, V)), ("K (V-g)", k_method(mm.M, mm.K, tab, rho0)),
                    ("p (Roger RFA, 4 lags)", modal_system(mm, fit_roger(tab, 4)).flutter(V, rho0))):
        c = r.critical()
        md.append(f"| {name} | {c.velocity:.1f} | {c.omega:.1f} |")
        if name == "PK":
            ax = r.plot(title="Goland wing - DLM + PK")
            ax[0].figure.savefig(os.path.join(FIG, "goland_vg.png"), dpi=110)
            plt.close("all")
    md.append("\n![Goland V-g](figures/goland_vg.png)\n")


def strip(md):
    from flaap.aero.strip import strip_gaf_table
    g = _goland()
    mm, ref, CHORD = g["mm"], g["REFERENCE"], g["CHORD"]
    rho0 = ref["density"]
    ny = mm.nodes[:, 1]
    o = np.argsort(ny)
    st = strip_gaf_table(ny[o], mm.shapes[o, 2, :], mm.shapes[o, 4, :], CHORD, 0.33,
                         np.r_[0, np.geomspace(0.005, 2, 40)], bref=CHORD / 2)
    cs = pk_method(mm.M, mm.K, st, rho0, np.linspace(40, 300, 521)).critical()
    ev, ew = _pct(cs.velocity, ref["strip_flutter_speed"]), _pct(cs.omega, ref["strip_flutter_omega"])
    _head(md, 3, "스트립 이론 (Theodorsen) 플러터 검증",
          "Theodorsen 2-D 비정상 공력을 날개 길이 방향 스트립에 적용한 플러터 해석 (`aero_method strip`).",
          [(EXP, "없음 — Goland 날개는 해석용 기준 모델로 시험 데이터 없음."),
           (OTH, "Goland (1945) 스트립 이론 플러터 해, Patil, Hodges & Cesnik (2001)."),
           (ANA, "Theodorsen 함수 자체가 공력 모델 (별도 비교 불필요)."),
           (CON, "없음.")],
          f"Goland (1945) 대비 V_f {ev:+.1f} %, w_f {ew:+.1f} %.",
          ("없음", "Goland 1945, Patil 2001", "-", f"V_f {ev:+.1f} %, w_f {ew:+.1f} %"))
    _sub(md, "3.1", OTH, f"Goland 날개, rho = {rho0} kg/m^3")
    md.append("| 출처 | V_f [m/s] | w_f [rad/s] |\n|---|---|---|")
    md.append(f"| **FLAAP 스트립 이론 + PK** | {cs.velocity:.1f} | {cs.omega:.1f} |")
    md.append(f"| Goland (1945) | {ref['strip_flutter_speed']} | {ref['strip_flutter_omega']} |")
    md.append("| Patil, Hodges & Cesnik (2001) | 135.6 | 70.2 |\n")
    md.append("같은 모델에서 DLM (3-D 날개끝 효과 포함)은 2.4의 153 m/s로 스트립 이론보다 높음.\n")


def tsd_flutter(md):
    mod = runpy.run_path(os.path.join(HERE, "agard445_tsd_flutter.py"))
    r = mod["run"]((0.499, 0.901, 0.954, 0.960), flat=True, verbose=False)
    ed, et = _pct(r[:, 2], r[:, 1]), _pct(r[:, 4], r[:, 1])
    ef = np.abs(_pct(r[:, 3], r[:, 2])).max()
    tr = r[:, 0] > 0.93
    _head(md, 4, "TSD 천음속 플러터 검증 (비선형 TSD 평균류 + 시간선형화 비정상 해 + PK)",
          "`flaap.aero.tsd`: 두께에 의한 천음속 비선형 정상 평균류(충격파 포함) 위의 시간선형화 "
          "비정상 공력 (`aero_method tsd`) 으로 구한 플러터 경계, 특히 천음속 딥.",
          [(EXP, "AGARD 445.6 weakened model 3 (Yates 1987), M 0.499-0.960, NACA 65A004."),
           (OTH, "미비교 — 문헌의 CAP-TSD/Euler 결과와의 정량 비교는 아직 없음."),
           (ANA, "해당 없음 (천음속 비선형 해석해 없음)."),
           (CON, "두께 0 TSD ('TSD flat') = 선형 DLM 이어야 함.")],
          f"천음속 M 0.954/0.960 실험 대비 TSD {et[tr][0]:+.1f} / {et[tr][1]:+.1f} % (같은 조건 DLM "
          f"{ed[tr][0]:+.1f} / {ed[tr][1]:+.1f} %) — 천음속 딥 재현. 아음속 M 0.499/0.901은 "
          f"{et[0]:+.1f} / {et[1]:+.1f} %. TSD flat vs DLM 최대 {ef:.1f} %.",
          ("AGARD 445.6 (M 0.5-0.96)", "미비교", "-",
           f"천음속 {et[tr][0]:+.1f} / {et[tr][1]:+.1f} % (DLM {ed[tr][0]:+.1f} / {ed[tr][1]:+.1f} %)"))
    _sub(md, "4.1", EXP, "AGARD 445.6 천음속 플러터 경계 (Yates 1987; 측정 진동수, FE 모드형상)")
    md.append("TSD 격자 51x20x22, NACA 65A004 두께.\n")
    md.append("| Mach | FSI 실험 | TSD | 오차 % | DLM | 오차 % | wf/wa 실험 | wf/wa TSD |\n"
              "|---|---|---|---|---|---|---|---|")
    for row, a, b in zip(r, et, ed):
        md.append(f"| {row[0]:.3f} | {row[1]:.4f} | {row[4]:.4f} | {a:+.1f} | {row[2]:.4f} | {b:+.1f} | "
                  f"{row[5]:.4f} | {row[7]:.4f} |")
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(r[:, 0], r[:, 1], "ko-", label="experiment (Yates 1987)")
    ax.plot(r[:, 0], r[:, 2], "bs--", label="DLM + PK")
    ax.plot(r[:, 0], r[:, 4], "r^-", label="TSD + PK")
    ax.set_xlabel("Mach"); ax.set_ylabel("FSI"); ax.grid(alpha=0.3); ax.legend()
    fig.savefig(os.path.join(FIG, "agard_tsd_fsi.png"), dpi=110)
    plt.close(fig)
    md.append("\n두께에 의한 천음속 평균류가 M 0.95 부근 플러터 속도를 낮춰, 선형 DLM이 놓치는 천음속 딥을 재현.\n")
    md.append("![AGARD TSD](figures/agard_tsd_fsi.png)\n")
    _sub(md, "4.2", CON, "두께 0 TSD vs DLM")
    md.append("| Mach | DLM | TSD flat | 차이 % |\n|---|---|---|---|")
    for row in r:
        md.append(f"| {row[0]:.3f} | {row[2]:.4f} | {row[3]:.4f} | {_pct(row[3], row[2]):+.1f} |")
    md.append("")


def tsd_lifting(md):
    TL = runpy.run_path(os.path.join(HERE, "tsd_lifting_validation.py"))
    c = TL["consistency"](verbose=False)
    rows, res = TL["naca0012_sweep"](verbose=False)
    ag = TL["agard_alpha"](verbose=False)
    TL["make_figure"](rows, res, ag, os.path.join(FIG, "tsd_lifting.png"))
    AR = TL["AR211"]
    cl = rows[:, 2]
    xs_eq = float(np.interp(AR["CL"], cl, rows[:, 3]))
    a_eq = float(np.interp(AR["CL"], cl, rows[:, 0]))
    conv = all(r.converged for r in res)
    e_cla = _pct(c["CLa"], c["CLa_helmbold"])
    _head(md, 5, "양력 TSD 평균류 검증 (받음각·캠버, 충격파 위치)",
          "`flaap.aero.tsd_lifting`: 받음각/캠버가 있는 비선형 TSD 정상 해 — 윗면/아랫면 충격파 위치와 "
          "세기, 양력, 그리고 그 평균류 위의 플러터 (`tsd_alpha`, `tsd_camber`).",
          [(EXP, "충격파 위치: 미비교 (실험 공백). NACA 0012 천음속 2-D 풍동 자료(예: Harris 1981, NASA "
                 "TM-81927)는 점성·벽 보정이 섞여 있어 아직 비교하지 않음. 받음각이 있는 AGARD 445.6 "
                 "플러터 실험은 없음 (alpha = 0 점만 Yates 1987, 4장)."),
           (OTH, "AGARD-AR-211 (1985) 비점성 Euler 합의값, NACA 0012, M 0.8, alpha 1.25 deg."),
           (ANA, "아음속 선형 극한 양력 기울기: Helmbold 식."),
           (CON, "alpha = 0이면 대칭 TSD와 동일, +-alpha 면 교환 대칭, 선형 TSD 극한.")],
          f"같은 단면 cl = {AR['CL']}에서 윗면 충격파 x/c = {xs_eq:.3f} (Euler {AR['xs_upper']}; "
          f"TSD 충격파가 {abs(xs_eq - AR['xs_upper']):.3f} c {'앞' if xs_eq < AR['xs_upper'] else '뒤'}). "
          f"Euler의 약한 아랫면 충격파는 TSD에서 alpha > 0.5 deg이면 사라짐. CL_alpha Helmbold 대비 "
          f"{e_cla:+.1f} %. 모든 받음각 수렴: {'예' if conv else '아니오'}.",
          ("미비교 (공백)", "AGARD-AR-211 Euler", "Helmbold",
           f"충격파 x/c {xs_eq:.3f} vs Euler {AR['xs_upper']}"))

    _sub(md, "5.1", EXP, "실험 비교 없음")
    md.append("충격파 위치에 대한 실험 비교는 아직 없음 (위 비교 기준 참조). 받음각이 있는 AGARD 445.6의 "
              "플러터 결과(5.4)도 실험이 없는 경향 예측임.\n")

    _sub(md, "5.2", OTH, "NACA 0012, M = 0.8: AGARD-AR-211 Euler 대비")
    md.append("종횡비 16 직사각형 날개의 뿌리 단면 (준 2-D), 시위 방향 40점.\n")
    md.append("| alpha [deg] | CL | 뿌리 cl | 윗면 충격파 x/c | dCp 윗면 | 아랫면 충격파 x/c | 수렴 |\n"
              "|---|---|---|---|---|---|---|")
    for row, r in zip(rows, res):
        lo = "-" if np.isnan(row[5]) else f"{row[5]:.3f}"
        md.append(f"| {row[0]:g} | {row[1]:.3f} | {row[2]:.3f} | {row[3]:.3f} | {row[4]:.2f} | {lo} | "
                  f"{'y' if r.converged else 'N'} |")
    md.append(f"| **Euler (AR-211), alpha {AR['alpha']}** | | **~{AR['CL']}** | **~{AR['xs_upper']}** | | "
              f"**~{AR['xs_lower']} (약함)** | |")
    md.append(f"\n같은 단면 cl = {AR['CL']} (alpha ~ {a_eq:.2f} deg; 유한 날개 뿌리 단면은 내리흐름 때문에 "
              f"같은 cl에 더 큰 alpha가 필요)에서 TSD 윗면 충격파 x/c = {xs_eq:.3f}. 약한 아랫면 충격파는 "
              "TSD에서 alpha ~ 0.5 deg 이상에서 사라짐 (고전 TSD 계수, Ames 항 없음; 앞전 근처 얇은 날개 "
              "경계조건).\n")

    _sub(md, "5.3", ANA, "선형 극한 양력 기울기 (M 0.5, 평판, A = 8)")
    md.append(f"| | CL_alpha [1/rad] |\n|---|---|\n| FLAAP 양력 TSD | {c['CLa']:.3f} |\n"
              f"| Helmbold | {c['CLa_helmbold']:.3f} |\n| 오차 | {e_cla:+.1f} % |\n")

    _sub(md, "5.4", CON, "정확 일치 검사와 받음각 평균류 위 플러터 경향")
    md.append("| 검사 | 결과 |\n|---|---|")
    md.append(f"| alpha = 0, 캠버 없음 vs 대칭 `TSDSolver` | max abs dphi {c['alpha0_phi']:.1e}, GAF 상대차 "
              f"{c['alpha0_gaf']:.1e} |")
    md.append(f"| alpha = +-1 deg: Cp 윗면(+) vs 아랫면(-) | {c['sign_cp']:.1e} (CL 반대칭 {c['sign_CL']:.1e}) |")
    md.append(f"| M 0.5 평판, alpha 1 deg vs 선형 (k = 0) TSD | 단면 cl 상대차 {c['linear_cl']:.1e} |\n")
    md.append("AGARD 445.6, M = 0.954, 평균 받음각에 따른 플러터 (alpha = 0만 실험 있음):\n")
    md.append("| 평균 alpha [deg] | CL | FSI | wf/wa |\n|---|---|---|---|")
    for row in ag:
        e = " (실험 0.3059)" if row[0] == 0 else ""
        md.append(f"| {row[0]:g} | {row[1]:.3f} | {row[3]:.4f}{e} | {row[4]:.3f} |")
    md.append(f"\n평균 받음각이 FSI를 낮춤 ({_pct(ag[-1, 3], ag[0, 3]):+.1f} % at {ag[-1, 0]:g} deg): 윗면 "
              "충격파가 강해져 충격파 운동 하중이 불안정 쪽으로 기여.\n")
    md.append("![lifting TSD](figures/tsd_lifting.png)\n")


def supersonic(md):
    txt = _log("supersonic_validation.py", "supersonic_results.txt")
    ch = _chunks(txt)
    try:
        mb = _find(r"^\s+(1\.\d+) Mach box nx = 80\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)%", ch["6"])[:2]
        e_exp = f"{mb[0][3]:+.1f} / {mb[1][3]:+.1f} %"
        e_pan = max(abs(v[0]) for v in _find(r"\s([+-]\d+\.\d+) %", ch["5"]))
        e_st = max(abs(v[0]) for v in _find(r"\s([+-]\d+\.\d+) %", ch["1-3"]))
        e_gr = max(v[3] for v in _find(r"^\s+(\d\.\d+)\s+(\d\.\d+)\s+[+-][\d.]+\s+[+-][\d.]+j\s+([\d.]+)%"
                                       r"\s+([\d.]+)%", ch["4"]))
    except (KeyError, IndexError, ValueError):
        e_exp, e_pan, e_st, e_gr = "(로그 참조)", float("nan"), float("nan"), float("nan")
    _head(md, 6, "초음속 공력·플러터 검증 (Mach box, piston theory)",
          "`flaap.aero.supersonic` (`aero_method machbox | piston`): 초음속 정상/비정상 압력과 플러터.",
          [(EXP, "AGARD 445.6 weakened model 3, M 1.072 / 1.141 (Yates 1987)."),
           (OTH, "패널 플러터 임계 동압: Hedgepeth (1957) / Dowell (1975) 수렴 Galerkin 해."),
           (ANA, "정상 양력 기울기 (Ackeret, 초음속 날개끝, Stewart 삼각익), 2-D 진동 평판 "
                 "(Garrick-Rubinow 적분), 2-모드 Galerkin 패널 식."),
           (CON, "없음.")],
          f"실험 대비 Mach box (nx 80, 측정 진동수) FSI {e_exp} — M 1.07-1.14의 낮은 초음속에서 "
          f"두께·점성을 뺀 선형 포텐셜 이론의 한계; piston theory는 이 Mach에서 적용 불가 (대부분 플러터 "
          f"없음, 1차 4/M은 M 1.141에서 FSI 3.6배). "
          f"패널 플러터 타 해석 대비 최대 {e_pan:.2f} %, 해석해 대비 정상 양력 최대 {e_st:.2f} %, "
          f"Garrick-Rubinow 대비 최대 {e_gr:.3f} %.",
          ("AGARD 445.6 (M 1.07/1.14)", "Dowell 패널 플러터", "Ackeret, Garrick-Rubinow 등",
           f"실험 대비 {e_exp}; 해석해 < {max(e_st, e_gr):.1f} %"))
    _sub(md, "6.1", EXP, "AGARD 445.6 초음속 플러터 경계 (Yates 1987)")
    _code(md, ch.get("6", ""))
    _sub(md, "6.2", OTH, "패널 플러터 임계 동압 (Hedgepeth / Dowell)")
    _code(md, ch.get("5", ""))
    _sub(md, "6.3", ANA, "정상 양력 기울기와 2-D 진동 평판")
    _code(md, ch.get("1-3", "") + "\n\n" + ch.get("4", ""))
    md.append("![supersonic](figures/supersonic_validation.png)\n")


def uvlm(md):
    txt = _log("uvlm_validation.py", "uvlm_results.txt")
    ch = _chunks(txt)
    try:
        vf = _find(r"Richardson nc -> inf: V_f\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+) %", ch["7"])[0]
        wag = _find(r"^\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)%\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)%",
                    ch["2"])
        e_w = max(abs(v[3]) for v in wag)
        e_k = max(abs(v[6]) for v in wag if v[0] >= 4)
        rich = [v[0] for v in _find(r"Richardson .*ratio ([\d.]+)", ch["3"])]
        e_th = 100 * max(abs(v - 1) for v in rich)
        e_ga = max(abs(v[0]) for v in _find(r"\s([+-]\d+\.\d+)%\s*$", ch["4"]))
        s_vf = f"{vf[0]:.1f} m/s vs Murua UVLM {vf[1]:.1f} m/s ({vf[2]:+.1f} %)"
    except (KeyError, IndexError, ValueError):
        s_vf, e_w, e_k, e_th, e_ga = "(로그 참조)", float("nan"), float("nan"), float("nan"), float("nan")
    _head(md, 7, "VLM / UVLM 공력 및 UVLM 플러터 검증",
          "`flaap.aero.uvlm`: 정상 VLM, 비선형 시간영역 UVLM (규정/자유 후류), 주파수영역 선형화 UVLM "
          "GAF와 그것을 쓴 플러터 (`model uvlm`, `aero_method uvlm`).",
          [(EXP, "없음 — 비교한 실험 자료 없음."),
           (OTH, "Goland 날개 플러터: Murua et al. (2012) UVLM, Wang (1999); 직사각형 날개 양력면 해석 값."),
           (ANA, "Helmbold 양력 기울기, 타원익 스팬 효율, Wagner/Küssner 함수, Theodorsen 조화 양력, "
                 "Garrick 추력."),
           (CON, "3-D 충격 출발의 정상 극한, 규정 vs 자유 후류, 주파수영역 vs 시간영역 상태공간.")],
          f"UVLM 플러터 (Richardson 외삽) {s_vf}. 해석해 대비: Wagner 최대 {e_w:.2f} %, Küssner (s >= 4) "
          f"{e_k:.2f} %, Theodorsen 외삽 {e_th:.2f} %, Garrick 추력 최대 {e_ga:.2f} %.",
          ("없음", "Murua 2012 (Goland)", "Helmbold, Wagner, Theodorsen 등", s_vf))
    _sub(md, "7.1", OTH, "선형화 UVLM + PK 플러터: Goland 날개 (Murua 2012)")
    _code(md, ch.get("7", ""))
    _sub(md, "7.2", ANA, "정상 VLM: 양력 기울기, 유도항력")
    _code(md, ch.get("1", ""))
    _sub(md, "7.3", ANA, "비선형 UVLM 지시 응답: Wagner / Küssner")
    _code(md, ch.get("2", ""))
    _sub(md, "7.4", ANA, "조화 피치·플런지: Theodorsen 양력, Garrick 추력")
    _code(md, ch.get("3", "") + "\n\n" + ch.get("4", ""))
    _sub(md, "7.5", ANA, "선형화 UVLM GAF vs Theodorsen")
    _code(md, ch.get("5", ""))
    _sub(md, "7.6", CON, "3-D 충격 출발: 정상 VLM 극한, 규정 vs 자유 후류")
    _code(md, ch.get("6", ""))
    md.append("![uvlm](figures/uvlm_validation.png)\n")


def _split_dashes(txt):
    parts = re.split(r"^-{20,}\s*$", txt, flags=re.M)
    return [p.strip("\n") for p in parts]


def freeplay(md):
    txt = _log("freeplay_lco.py", "freeplay_results.txt")
    parts = _split_dashes(txt)
    try:
        vfl = _find(r"linear flutter speed .*?:\s+([\d.]+) m/s\s+\(published ([\d.]+)", txt)[0]
        df = _find(r"^\s+[\d.]+\s+lco\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([+-][\d.]+)%\s+[\d.]+\s+[\d.]+\s+([+-][\d.]+)%",
                   parts[0])
        eA, ew = max(abs(v[0]) for v in df), max(abs(v[1]) for v in df)
        s = f"{vfl[0]:.2f} m/s vs {vfl[1]:.1f} m/s"
    except (IndexError, ValueError):
        vfl, eA, ew, s = None, float("nan"), float("nan"), "(로그 참조)"
    _head(md, 8, "조종면 프리플레이 비선형 LCO 검증",
          "`flaap.nonlinear`: 조종면 프리플레이의 한계 사이클 진동(LCO) — 기술함수(DF)법과 구간선형 "
          "시간 적분, 분기도·이력 현상.",
          [(EXP, "LCO 진폭 미비교 — Conner et al. (1997) 시험 장치(프리플레이 +-3.57 deg, 측정 감쇠 행렬, "
                 "5차 피치 강성 다항식)는 FLAAP 모델과 달라 비교하지 않음."),
           (OTH, "Vasconcellos et al. (2016) 선형 플러터 속도 (프리플레이 없음, 9.2 m/s)."),
           (ANA, "없음."),
           (CON, "DF vs 시간 적분 (진폭·진동수), 진폭의 delta0 비례성, 상·하향 속도 스윕 이력, "
                 "Goland + 플랩 3-D 모델의 DF vs 시간 적분.")],
          f"선형 플러터 속도 {s}. 자기일관성: DF vs 시간 적분 LCO 진폭 최대 {eA:.1f} %, 진동수 {ew:.1f} %.",
          ("미비교 (공백)", "Vasconcellos 2016 V_f", "없음", f"V_f {s}"))
    _sub(md, "8.1", OTH, "Conner/Vasconcellos 3자유도 단면: 선형 플러터 속도")
    if vfl:
        md.append(f"| | V_f [m/s] |\n|---|---|\n| FLAAP (Jones-Wagner RFA, Conner 감쇠) | {vfl[0]:.3f} |\n"
                  f"| Vasconcellos et al. (2016) | {vfl[1]:.1f} |\n| 오차 | {_pct(vfl[0], vfl[1]):+.1f} % |\n")
    _sub(md, "8.2", CON, "기술함수 vs 시간 적분, 분기도 (단면 / Goland + 플랩)")
    _code(md, parts[0].split("\n", 3)[-1] if len(parts) > 1 else txt)
    if len(parts) > 1:
        _code(md, parts[1])
    md.append("![bifurcation](figures/freeplay_bifurcation.png)\n")


def ase(md):
    txt = _log("ase_flutter_suppression.py", "ase_results.txt")
    parts = _split_dashes(txt)
    try:
        v_ol = _find(r"open-loop p-method .*?V_f = ([\d.]+)", txt)[0][0]
        v_k0 = _find(r"K = 0:\s+V_f = ([\d.]+)", txt)[0][0]
        laws = _find(r"law [AB]:.*?V_f = ([\d.]+) m/s \(([+-][\d.]+) %\)", txt)
        mrows = re.findall(r"^\s+[\d.]+\s+\d+\s+\d+\s+(True|False)\s+(True|False)", txt, re.M)
        agree = all(a == b for a, b in mrows)
        s = (f"K = 0 폐루프 V_f {v_k0:.4f} = 개루프 p-method {v_ol:.4f} m/s; Nyquist 판정과 폐루프 고유값 "
             f"판정 {'모든' if agree else '일부'} 속도에서 일치 ({len(mrows)}점). 제어 법칙 A/B로 V_f "
             f"{laws[0][1]:+.0f} / {laws[1][1]:+.0f} %.")
    except (IndexError, ValueError):
        s = "(로그 참조)"
    _head(md, 9, "공력서보탄성 (ASE) 플러터 억제·돌풍하중 완화 — 자기일관성 검사",
          "`flaap.ase`: 작동기/센서/제어 법칙을 포함한 폐루프 플러터, 이득·위상 여유, Nyquist, GLA.",
          [(EXP, "없음."),
           (OTH, "없음 — 같은 제어 법칙의 문헌 결과 없음."),
           (ANA, "없음."),
           (CON, "이득 0 = 개루프 플러터, Nyquist vs 고유값, 이득 여유 x0.97/x1.03에서 근 교차 확인.")],
          s, ("없음", "없음", "없음", "자기일관성만"))
    _sub(md, "9.1", CON, "Conner/Vasconcellos 단면 플러터 억제")
    _code(md, parts[0].split("\n", 3)[-1] if len(parts) > 1 else txt)
    if len(parts) > 1:
        _sub(md, "9.2", CON, "Goland + 플랩 돌풍하중 완화 (GLA)")
        _code(md, parts[1])
    md.append("![ASE](figures/ase_vg_rootlocus.png)\n")


def gust(md):
    g = _goland()
    from flaap.models.goland import goland_modal
    tab = g["tab"]
    mm2 = goland_modal(4, damping=0.02)
    ss = modal_system(mm2, fit_roger(tab, 4)).state_space(100.0, 1.225)
    r = G.discrete_gust(ss, 100.0, 20.0, 10.0, t_end=1.5, dt=1e-3)

    def H(om):
        return G.frequency_response_gust(mm2.M, mm2.K, mm2.C(), tab, 100.0, 1.225, om)[0]
    yf = G.fft_response(r.t, r.wg, H).real
    peaks = []
    for i in range(2):
        mx, mn = r.peak(f"eta:mode{i + 1}")
        big = mn if abs(mn) > abs(mx) else mx
        ref = yf[:, i].min() if abs(mn) > abs(mx) else yf[:, i].max()
        peaks.append(_pct(big, ref))
    _head(md, 10, "돌풍 응답 검증 (1-cos 이산 돌풍) — 자기일관성 검사",
          "`flaap.gust`: RFA 상태공간 시간영역 돌풍 응답.",
          [(EXP, "없음."), (OTH, "없음."), (ANA, "없음."),
           (CON, "같은 GAF 표의 주파수영역 해 (표 보간 GAF + FFT) 와 비교.")],
          f"주요 모드 1/2 최대 응답 TD vs FD 차이 {peaks[0]:+.2f} / {peaks[1]:+.2f} %.",
          ("없음", "없음", "없음", "TD vs FD 일치"))
    _sub(md, "10.1", CON, "Goland 날개: 시간영역 (RFA) vs 주파수영역 (GAF + FFT)")
    md.append("1-cos 돌풍 H = 20 m, U = 10 m/s, V = 100 m/s, zeta = 2 %\n")
    md.append("| output | TD max | FD max | TD min | FD min |\n|---|---|---|---|---|")
    for i in range(4):
        lab = f"eta:mode{i + 1}"
        mx, mn = r.peak(lab)
        md.append(f"| {lab} | {mx:.4e} | {yf[:, i].max():.4e} | {mn:.4e} | {yf[:, i].min():.4e} |")
    fig, ax = plt.subplots(figsize=(7, 3.5))
    for i in range(2):
        ax.plot(r.t, r.y[:, r.labels.index(f"eta:mode{i + 1}")], label=f"TD mode {i + 1}")
        ax.plot(r.t, yf[:, i], "--", label=f"FD mode {i + 1}")
    ax.legend(); ax.grid(alpha=0.3); ax.set_xlabel("t [s]")
    fig.savefig(os.path.join(FIG, "goland_gust.png"), dpi=110)
    plt.close(fig)
    md.append("\n![gust](figures/goland_gust.png)\n")


def trim(md):
    from flaap.loads import MonitorStation, modal_load_matrices, station_output_group
    from flaap.maneuver import checked_pitch_maneuver, short_period
    from flaap.models.demo_aircraft import demo_aircraft
    from flaap.trim import G0, StaticAeroelastic, load_rows_from_table
    ac, mesh, modes = demo_aircraft(4)
    tab = AeroModel(mesh, modes, 0.0, 0.9144).table(np.r_[0, np.geomspace(0.01, 2.0, 18)])
    st = [MonitorStation("root", [0, 0.05, 0], [0, 1, 0], ("Fz", "Mx", "My"))]
    rows = load_rows_from_table(tab, mesh, st)
    inert = {lab: r for lab, r in zip(st[0].labels, modal_load_matrices(ac, st[0])[0])}
    sa = StaticAeroelastic.from_table(ac.M, ac.K, tab, load_rows=rows, inertia_rows=inert)
    V, rho = 100.0, 1.225
    tr = sa.solve(0.5 * rho * V * V, fixed={"plunge": 0.0}, free=["elevator"], accel={"plunge": G0})
    W = ac.M[0, 0] * G0
    rfa = fit_roger(tab, 4)
    sysm = AeroelasticSystem(ac.M, ac.K, rfa, C=ac.C, n_ctrl=1, names=ac.names + ["elevator"])
    station_output_group(sysm, tab, mesh, ac, st, "both")
    wn, z, _ = short_period(sysm, V, rho)
    r = checked_pitch_maneuver(sysm, V, rho, "elevator", nz_limit=2.5)
    _head(md, 11, "정적 공탄성 트림·기동하중 — 자기일관성 검사",
          "`flaap.trim`, `flaap.maneuver`, `flaap.loads`: 유연 기체 1-g 트림, CS-25.331(c)(2) 점검 피치 "
          "기동, 모니터 단면 하중.",
          [(EXP, "없음."), (OTH, "없음 — 시연 기체 (기준 해 없음)."), (ANA, "없음."),
           (CON, "트림 양력 = 반 기체 무게, 기동 중 n_z 한계 도달.")],
          f"트림 양력 {tr.forces[0]:.1f} N = 반 무게 {W:.1f} N (상대차 {abs(tr.forces[0] / W - 1):.1e}); 기동 "
          f"n_z 최대 {r.nz.max():.2f} (한계 2.5).",
          ("없음", "없음", "없음", "자기일관성만"))
    _sub(md, "11.1", CON, "시연 유연 기체")
    md.append(f"1-g 트림 V = {V} m/s: alpha = {np.degrees(tr['pitch']):.3f} deg, 승강타 = "
              f"{np.degrees(tr['elevator']):.3f} deg, 양력 = {tr.forces[0]:.1f} N (반 무게 {W:.1f} N), "
              "뿌리 하중 " + ", ".join(f"{k} = {v:.1f}" for k, v in tr.loads.items()) + "\n")
    md.append(f"단주기: wn = {wn:.2f} rad/s, zeta = {z:.3f}. n_z = 2.5 점검 기동: 승강타 최대 "
              f"{np.degrees(np.abs(r.delta).max()):.2f} deg, 뿌리 Fz 증분 (FSM) {r.peak('root.Fz.fsm')[0]:.0f} N, "
              f"Mx {r.peak('root.Mx.fsm')[0]:.0f} N m\n")
    fl = sysm.flutter(np.linspace(40, 250, 106), rho)
    md.append("자유 비행 플러터 (p-method): " + (str(fl.critical()) if fl.critical() else "없음") + "\n")
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(7, 5))
    ax[0].plot(r.t, np.degrees(r.delta)); ax[0].set_ylabel("elevator [deg]")
    ax[1].plot(r.t, r.nz); ax[1].set_ylabel("n_z"); ax[1].set_xlabel("t [s]")
    for a in ax:
        a.grid(alpha=0.3)
    fig.savefig(os.path.join(FIG, "maneuver.png"), dpi=110)
    plt.close(fig)
    md.append("![maneuver](figures/maneuver.png)\n")


def gaps(md):
    md.append("## 12. 실험 비교가 없는 항목 (향후 과제)\n")
    md.append("| 기능 | 현재 최상위 기준 | 필요한 실험/자료 |\n|---|---|---|")
    md.append("| 양력 TSD 충격파 위치 (5장) | 타 해석 (AGARD-AR-211 Euler) | NACA 0012 등 천음속 2-D/3-D 풍동 "
              "압력 분포 (벽·점성 보정 고려) |")
    md.append("| 받음각이 있는 천음속 플러터 (5.4) | 없음 (경향) | 평균 받음각을 바꾼 천음속 플러터 시험 |")
    md.append("| TSD 천음속 플러터 (4장) | 실험 있음; 타 해석 미비교 | 문헌 CAP-TSD / Euler AGARD 445.6 결과 |")
    md.append("| 스트립 이론·UVLM 플러터 (3, 7장) | 타 해석 | 3-D 날개 플러터 풍동 시험 (예: AGARD 445.6 "
              "저 Mach로 UVLM 비교) |")
    md.append("| 프리플레이 LCO 진폭 (8장) | 자기일관성 | Conner (1997) 시험 장치 그대로의 모델 (+-3.57 deg, "
              "비선형 강성) |")
    md.append("| ASE, 돌풍, 트림·기동 (9-11장) | 자기일관성 | 시험 또는 독립 코드 (예: NASTRAN SOL 144/146) "
              "결과 |\n")


SECTIONS = {"structures": structures, "dlm": dlm, "strip": strip, "tsd_flutter": tsd_flutter,
            "tsd_lifting": tsd_lifting, "supersonic": supersonic, "uvlm": uvlm, "freeplay": freeplay,
            "ase": ase, "gust": gust, "trim": trim, "gaps": gaps}


# ====================================================================== assembly
def _split_sections(text):
    """REPORT.md text -> (preamble lines, [(heading, block text)])."""
    pre, secs = [], []
    for line in text.splitlines():
        if line.startswith("## "):
            secs.append([line, []])
        elif secs:
            secs[-1][1].append(line)
        else:
            pre.append(line)
    return pre, [(h, "\n".join([h] + b)) for h, b in secs]


def _section_key(heading):
    num = heading[3:].split(" ", 1)[0].rstrip(".")
    digits = "".join(ch for ch in num if ch.isdigit())
    return (int(digits or 99), num)


def _overview(blocks):
    """Chapter 0: principles + summary table rebuilt from the chapters' summary comments."""
    md = ["## 0. 검증 원칙과 요약\n",
          "- 각 장의 제목은 **무엇을 검증했는지 (해석 기능)** 이다. Goland 날개, AGARD 445.6 등 사례는 그 "
          "아래 항목이다.",
          "- 비교 기준 우선순위: **1순위 실험 -> 2순위 독립된 다른 해석/코드 -> 3순위 해석해**. 상위 기준이 "
          "없으면 그 사실과 이유를 적는다.",
          "- 코드 내부 비교 (극한 일치, 해법 간 일치, 입력 경로 왕복) 는 **자기일관성** 으로 따로 표시하며 "
          "검증 근거로 세지 않는다.\n",
          "| 장 | 검증 대상 | 1순위 실험 | 2순위 타 해석 | 3순위 해석해 | 대표 결과 |",
          "|---|---|---|---|---|---|"]
    for h in sorted(blocks, key=_section_key):
        m = re.search(r"<!-- summary \| (.*?) -->", blocks[h])
        if m:
            cells = [c.strip() for c in m.group(1).split(" | ")]
            md.append(f"| {_section_key(h)[1]} | " + " | ".join(cells) + " |")
    return "\n".join(md) + "\n"


def main(argv=None):
    global RERUN
    names = list(sys.argv[1:] if argv is None else argv)
    if "--rerun" in names:
        RERUN = True
        names.remove("--rerun")
    unknown = [n for n in names if n not in SECTIONS]
    if unknown:
        raise SystemExit(f"unknown section(s) {unknown}; choose from {list(SECTIONS)}")
    os.makedirs(FIG, exist_ok=True)
    report = os.path.join(HERE, "REPORT.md")
    stamp = f"Generated {time.strftime('%Y-%m-%d %H:%M')} by `validation/report.py`."
    blocks = {}
    if names and os.path.exists(report):           # splice the regenerated chapters into the report
        with open(report, encoding="utf-8") as fh:
            _, secs = _split_sections(fh.read())
        blocks = {h: b for h, b in secs if _section_key(h)[0] != 0}
    for name in names or list(SECTIONS):
        t0 = time.perf_counter()
        md = []
        SECTIONS[name](md)
        b = "\n".join(md)
        num = _section_key(b.split("\n", 1)[0])[0]
        blocks = {h: v for h, v in blocks.items() if _section_key(h)[0] != num}
        blocks[b.split("\n", 1)[0]] = b
        print(f"{name}: {time.perf_counter() - t0:.1f} s", flush=True)
    text = "\n".join(["# FLAAP (Python) 검증 보고서\n", stamp + "\n", _overview(blocks)]
                     + [blocks[h].rstrip() + "\n" for h in sorted(blocks, key=_section_key)])
    with open(report, "w", encoding="utf-8") as fh:
        fh.write(text)
    print("wrote", report)


if __name__ == "__main__":
    main()
