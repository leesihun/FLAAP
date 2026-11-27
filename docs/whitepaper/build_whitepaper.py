"""Build the FLAAP technical whitepaper (Korean) as a PDF.

    python tools/run_all_configs.py          # (re)generate output/* used for figures
    python validation/report.py              # (optional) regenerate validation/REPORT.md
    python docs/whitepaper/build_whitepaper.py

Output: docs/FLAAP_Whitepaper.pdf.  Figures are regenerated from the run outputs
into docs/whitepaper/figures, equations are rendered with matplotlib mathtext.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
EQD = HERE / "eq"
OUT = ROOT / "output"
from flaap.io.results_txt import read_result  # noqa: E402
VAL = ROOT / "validation"
PDF = ROOT / "docs" / "FLAAP_Whitepaper.pdf"
FIG.mkdir(exist_ok=True)
EQD.mkdir(exist_ok=True)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402,F401
from matplotlib import mathtext  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402

from reportlab.lib import colors  # noqa: E402
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY  # noqa: E402
from reportlab.lib.pagesizes import A4  # noqa: E402
from reportlab.lib.styles import ParagraphStyle  # noqa: E402
from reportlab.lib.units import mm  # noqa: E402
from reportlab.pdfbase import pdfmetrics  # noqa: E402
from reportlab.pdfbase.ttfonts import TTFont  # noqa: E402
from reportlab.platypus import (BaseDocTemplate, CondPageBreak, Frame, Image, KeepTogether,  # noqa: E402
                                NextPageTemplate, PageBreak, PageTemplate, Paragraph, Preformatted,
                                Spacer, Table, TableStyle)
from reportlab.platypus.tableofcontents import TableOfContents  # noqa: E402

# --------------------------------------------------------------------------- style
INK = "#17211e"
GREEN = "#126b57"
AMBER = "#b76b2a"
RED = "#b84d48"
BLUE = "#3f6f91"
GREY = "#6b7570"
PAPER = "#f3f1eb"
PANEL = "#fbfaf6"
ACID = "#ccf48d"
SERIES = [GREEN, AMBER, BLUE, RED, "#7a5c99", "#4b8b3b", "#a0526d", "#3a3a3a"]

plt.rcParams.update({
    "font.family": ["Malgun Gothic", "DejaVu Sans"], "axes.unicode_minus": False, "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.edgecolor": "#55605b", "axes.linewidth": 0.7, "axes.grid": True, "grid.color": "#d9d6cc",
    "grid.linewidth": 0.6, "axes.prop_cycle": matplotlib.cycler(color=SERIES), "legend.fontsize": 8,
    "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 200, "mathtext.fontset": "cm",
    "axes.spines.top": False, "axes.spines.right": False,
})

FONT_DIR = Path(r"C:\Windows\Fonts")
pdfmetrics.registerFont(TTFont("Malgun", str(FONT_DIR / "malgun.ttf")))
pdfmetrics.registerFont(TTFont("MalgunBd", str(FONT_DIR / "malgunbd.ttf")))
pdfmetrics.registerFont(TTFont("Consolas", str(FONT_DIR / "consola.ttf")))
from reportlab.pdfbase.pdfmetrics import registerFontFamily  # noqa: E402

registerFontFamily("Malgun", normal="Malgun", bold="MalgunBd", italic="Malgun", boldItalic="MalgunBd")

PAGE_W, PAGE_H = A4
MARGIN = 20 * mm
TEXT_W = PAGE_W - 2 * MARGIN

S = {
    "body": ParagraphStyle("body", fontName="Malgun", fontSize=9.6, leading=15.2, alignment=TA_JUSTIFY,
                           textColor=colors.HexColor(INK), spaceAfter=5),
    "small": ParagraphStyle("small", fontName="Malgun", fontSize=8.2, leading=11.5,
                            textColor=colors.HexColor(GREY)),
    "caption": ParagraphStyle("caption", fontName="Malgun", fontSize=8.2, leading=11.5, alignment=TA_CENTER,
                              textColor=colors.HexColor(GREY), spaceBefore=3, spaceAfter=10),
    "h1": ParagraphStyle("h1", fontName="MalgunBd", fontSize=17, leading=23, textColor=colors.HexColor(GREEN),
                         spaceBefore=6, spaceAfter=10),
    "h2": ParagraphStyle("h2", fontName="MalgunBd", fontSize=12.5, leading=17, textColor=colors.HexColor(INK),
                         spaceBefore=12, spaceAfter=6),
    "h3": ParagraphStyle("h3", fontName="MalgunBd", fontSize=10.5, leading=14, textColor=colors.HexColor(GREEN),
                         spaceBefore=8, spaceAfter=4),
    "bullet": ParagraphStyle("bullet", fontName="Malgun", fontSize=9.6, leading=14.6, leftIndent=12,
                             bulletIndent=2, textColor=colors.HexColor(INK), spaceAfter=2),
    "code": ParagraphStyle("code", fontName="Consolas", fontSize=7.9, leading=10.4, textColor=colors.HexColor(INK),
                           backColor=colors.HexColor(PANEL), borderColor=colors.HexColor("#d9d6cc"),
                           borderWidth=0.6, borderPadding=6, leftIndent=6, rightIndent=6, spaceBefore=4,
                           spaceAfter=10),
    "cell": ParagraphStyle("cell", fontName="Malgun", fontSize=8, leading=10.6, textColor=colors.HexColor(INK)),
    "cellb": ParagraphStyle("cellb", fontName="MalgunBd", fontSize=8, leading=10.6, textColor=colors.white),
    "toc1": ParagraphStyle("toc1", fontName="MalgunBd", fontSize=10.5, leading=17, leftIndent=0),
    "toc2": ParagraphStyle("toc2", fontName="Malgun", fontSize=9.2, leading=14, leftIndent=16),
    "box": ParagraphStyle("box", fontName="Malgun", fontSize=9.2, leading=14, textColor=colors.HexColor(INK),
                          backColor=colors.HexColor("#e9f3ec"), borderColor=colors.HexColor(GREEN),
                          borderWidth=0.8, borderPadding=8, spaceBefore=6, spaceAfter=12),
}

# --------------------------------------------------------------------------- helpers
_eq_no = [0]
_fig_no = [0]
_tab_no = [0]


def P(text, style="body"):
    return Paragraph(text, S[style])


def bullets(items):
    return [Paragraph(t, S["bullet"], bulletText="•") for t in items]


def code(text):
    return Preformatted(text.strip("\n"), S["code"])


def eq(tex, size=11.5, number=True):
    """Display equation rendered by mathtext, centred, with an equation number."""
    _eq_no[0] += 1
    name = EQD / f"eq{_eq_no[0]:03d}.png"
    mathtext.math_to_image(f"${tex}$", str(name), dpi=300, prop=FontProperties(size=size), format="png")
    from PIL import Image as PILImage
    w, h = PILImage.open(name).size
    img = Image(str(name), width=w * 72 / 300, height=h * 72 / 300)
    num = Paragraph(f"({_eq_no[0]})" if number else "", S["small"])
    t = Table([["", img, num]], colWidths=[18 * mm, TEXT_W - 36 * mm, 18 * mm])
    t.setStyle(TableStyle([("ALIGN", (1, 0), (1, 0), "CENTER"), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                           ("ALIGN", (2, 0), (2, 0), "RIGHT"), ("TOPPADDING", (0, 0), (-1, -1), 3),
                           ("BOTTOMPADDING", (0, 0), (-1, -1), 3)]))
    return t


def figure(path, caption, width=None, max_h=105 * mm):
    path = Path(path)
    if not path.exists():
        return [P(f"<i>(그림 없음: {path.name})</i>", "small")]
    _fig_no[0] += 1
    from PIL import Image as PILImage
    w, h = PILImage.open(path).size
    width = width or TEXT_W
    height = width * h / w
    if height > max_h:
        height = max_h
        width = height * w / h
    return [KeepTogether([Image(str(path), width=width, height=height),
                          Paragraph(f"<b>그림 {_fig_no[0]}.</b> {caption}", S["caption"])])]


def table(rows, col_w=None, caption=None, header=True, align_right_from=None, font=8, keep=True):
    _tab_no[0] += 1 if caption else 0
    st_c = ParagraphStyle("c", parent=S["cell"], fontSize=font, leading=font * 1.32)
    st_h = ParagraphStyle("h", parent=S["cellb"], fontSize=font, leading=font * 1.32)
    data = []
    for i, r in enumerate(rows):
        data.append([Paragraph(str(c), st_h if (header and i == 0) else st_c) for c in r])
    ncol = len(rows[0])
    col_w = col_w or [TEXT_W / ncol] * ncol
    if sum(col_w) > TEXT_W * 1.001:
        col_w = [c * TEXT_W / sum(col_w) for c in col_w]
    t = Table(data, colWidths=col_w, repeatRows=1 if header else 0, hAlign="CENTER")
    style = [("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cfcabd")),
             ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
             ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
             ("ROWBACKGROUNDS", (0, 1 if header else 0), (-1, -1),
              [colors.white, colors.HexColor(PANEL)])]
    if header:
        style.append(("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#122820")))
    t.setStyle(TableStyle(style))
    out = []
    if caption:
        out.append(Paragraph(f"<b>표 {_tab_no[0]}.</b> {caption}", S["caption"]))
    out.append(t)
    out.append(Spacer(1, 8))
    return [KeepTogether(out)] if (keep and len(rows) < 24) else out


def load(rel):
    """Read a FLAAP result text file (output/<run>/<name>.txt) into its contract dict."""
    p = OUT / rel
    if p.suffix == ".json":
        p = p.with_suffix(".txt")
    return read_result(p) if p.exists() else None


def mdinline(x):
    """Markdown inline text -> reportlab paragraph markup."""
    import re
    x = x.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    x = re.sub(r"`([^`]+)`", r"<font face='Consolas'>\1</font>", x)
    x = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", x)
    x = re.sub(r"(?<![\w*])\*([^*\n]+)\*(?![\w*])", r"<i>\1</i>", x)
    return x


def md_chapter(md_text, n):
    """Text of REPORT.md chapter ``## n. ...`` (without the next chapter)."""
    import re
    m = re.search(rf"^## {n}\. .*$", md_text, flags=re.M)
    if not m:
        return ""
    nxt = re.search(r"^## ", md_text[m.end():], flags=re.M)
    return md_text[m.start(): m.end() + (nxt.start() if nxt else len(md_text))]


def val_basis(md_text, n):
    """Capability / reference-level table of validation chapter n, parsed from REPORT.md."""
    import re
    ch = md_chapter(md_text, n)
    if not ch:
        return []
    rows = [["구분", "내용"]]
    m = re.search(r"\*\*검증 대상:\*\*\s*(.+)", ch)
    if m:
        rows.append(["<b>검증 대상</b>", mdinline(m.group(1))])
    for lab, pat in (("1순위 실험", "1순위 실험"), ("2순위 타 해석", "2순위 타 해석"), ("3순위 해석해", "3순위 해석해"),
                     ("자기일관성 (검증 근거 아님)", "자기일관성")):
        m = re.search(rf"^- \*\*{pat}[^*]*\*\* — (.+)$", ch, flags=re.M)
        if m:
            rows.append([f"<b>{lab}</b>", mdinline(m.group(1))])
    m = re.search(r"\*\*판정:\*\*\s*(.+)", ch)
    if m:
        rows.append(["<b>판정</b>", mdinline(m.group(1))])
    return table(rows, [34 * mm, 136 * mm], font=7.8) if len(rows) > 1 else []


def val_summary_rows(md_text):
    """Chapter summary table of REPORT.md (<!-- summary | ... --> comments), numbered as in this chapter."""
    import re
    rows = [["절", "검증 대상 (기능)", "1순위 실험", "2순위 타 해석", "3순위 해석해", "대표 결과"]]
    for i, ln in enumerate(re.findall(r"^<!-- summary \| (.*) -->$", md_text, flags=re.M), 1):
        f = [c.strip() for c in ln.split(" | ")]
        if len(f) == 5:
            rows.append([f"4.{i}"] + [mdinline(c) for c in f])
    return rows if len(rows) > 1 else None


def md_table(md_text, after_heading, nth=0):
    """Return the nth markdown table (list of rows) after a heading containing `after_heading`."""
    i = md_text.find(after_heading)
    if i < 0:
        return None
    lines = md_text[i:].splitlines()
    tables, cur = [], []
    for ln in lines[1:]:
        if ln.startswith("#") and tables:
            break
        if ln.strip().startswith("|"):
            if set(ln.replace("|", "").strip()) <= set("-: "):
                continue
            cur.append([mdinline(c.strip()) for c in ln.strip().strip("|").split("|")])
        elif cur:
            tables.append(cur)
            cur = []
    if cur:
        tables.append(cur)
    return tables[nth] if len(tables) > nth else None


def fmt(x, nd=1):
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "-"


# --------------------------------------------------------------------------- figures
def fig_architecture():
    fig, ax = plt.subplots(figsize=(9.2, 4.3))
    ax.set_axis_off()
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 48)

    def box(x, y, w, h, title, sub, fc, tc="white"):
        ax.add_patch(matplotlib.patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                                       fc=fc, ec="none"))
        ax.text(x + w / 2, y + h * 0.64, title, ha="center", va="center", color=tc, fontsize=8.6, weight="bold")
        ax.text(x + w / 2, y + h * 0.28, sub, ha="center", va="center", color=tc, fontsize=6.6)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color="#55605b", lw=1.0))

    box(1, 34, 21, 11, "Input", "flat config / Studio GUI", "#122820")
    box(1, 19, 21, 11, "Structure", "Goland, AGARD, plate FE\nNastran / ANSYS / Radioss / ccx\nSTL -> shell/solid FE", GREEN)
    box(1, 3, 21, 11, "Aero geometry", "surfaces, controls,\nsymmetry, stations", GREEN)
    box(27, 19, 18, 11, "Spline", "IPS / TPS / RBF\nbeam / legacy", BLUE)
    box(27, 3, 18, 11, "Aerodynamics", "VLM + DLM (subsonic)\nTSD (transonic)", BLUE)
    box(50, 11, 18, 11, "GAF  Q(k)", "modes, gust,\ncontrols, dCp", AMBER)
    box(50, 30, 18, 11, "Roger RFA", "state space\nx' = A x + B u", AMBER)
    box(73, 37, 25, 9, "Flutter", "PK / K / p, identification", "#2f3b36")
    box(73, 26.5, 25, 9, "Gust & manoeuvre", "1-cos, von Karman, CS-25.331", "#2f3b36")
    box(73, 16, 25, 9, "Static aeroelastic", "trim, divergence, reversal", "#2f3b36")
    box(73, 5.5, 25, 9, "Freeplay LCO / ASE", "describing fn, servo loops", "#2f3b36")
    box(73, -4, 25, 8, "Output contract", "JSON + PNG -> Studio, PDF", ACID, tc=INK)
    arrow(11.5, 34, 11.5, 30.5)
    arrow(22.5, 24.5, 26.5, 24.5)
    arrow(22.5, 8.5, 26.5, 8.5)
    arrow(36, 19, 36, 14.5)
    arrow(45.5, 8.5, 49.5, 14)
    arrow(45.5, 24.5, 49.5, 19)
    arrow(59, 22.5, 59, 29.5)
    for yy in (41.5, 31, 20.5, 10):
        arrow(68.5, 35.5 if yy > 25 else 16.5, 72.5, yy)
    ax.set_ylim(-5, 48)
    fig.tight_layout()
    p = FIG / "architecture.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def _deformed_mesh_axes(ax, corners, h, scale, title, cmap="RdBu_r"):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    c = corners.copy()
    c[:, :, 2] += h * scale
    vmax = np.abs(h).max() or 1
    col = plt.get_cmap(cmap)(0.5 + 0.5 * h.mean(axis=1) / vmax)
    ax.add_collection3d(Poly3DCollection(c, facecolors=col, edgecolor="#2a2a2a", linewidths=0.12))
    pts = c.reshape(-1, 3)
    lo, hi = pts.min(0), pts.max(0)
    ext = np.maximum(hi - lo, 0.22 * (hi - lo).max())
    ctr = 0.5 * (lo + hi)
    for setter, i in ((ax.set_xlim, 0), (ax.set_ylim, 1), (ax.set_zlim, 2)):
        setter(ctr[i] - ext[i] / 2, ctr[i] + ext[i] / 2)
    ax.set_box_aspect(tuple(ext))
    ax.set_axis_off()
    ax.view_init(26, -128)
    ax.set_title(title, fontsize=8.5)


def fig_modes(run, name, nmax=4):
    m = load(f"{run}/modes.txt")
    a = load(f"{run}/aero_mesh.txt")
    if m is None or a is None or "aero_h" not in m:
        return None
    corners = np.array(a["corners"])
    H = np.array(m["aero_h"])
    n = min(nmax, H.shape[0])
    span = np.ptp(corners.reshape(-1, 3), axis=0).max()
    fig = plt.figure(figsize=(9.0, 2.6 if n <= 4 else 4.8))
    for i in range(n):
        ax = fig.add_subplot(1, n, i + 1, projection="3d")
        _deformed_mesh_axes(ax, corners, H[i], 0.16 * span, f"{m['labels'][i]}\n{m['freq_hz'][i]:.2f} Hz")
    fig.subplots_adjust(left=0, right=1, bottom=0, top=0.9, wspace=0)
    p = FIG / f"modes_{name}.png"
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_cp(run, name, k_target=0.1, modes=(0, 1)):
    cp = load(f"{run}/cp.txt")
    a = load(f"{run}/aero_mesh.txt")
    if cp is None or a is None:
        return None
    from matplotlib.collections import PolyCollection
    k = np.array(cp["k"])
    ik = int(np.argmin(np.abs(k - k_target)))
    corners = np.array(a["corners"])
    poly = corners[:, :, [1, 0]]                      # span to the right, chord downwards
    fig, axs = plt.subplots(1, 2 * len(modes), figsize=(9.2, 3.0))
    for j, im in enumerate(modes):
        for part, (key, lab) in enumerate((("re", "Re"), ("im", "Im"))):
            ax = axs[2 * j + part]
            v = np.array(cp[key])[ik, im]
            vmax = np.abs(v).max() or 1
            pc = PolyCollection(poly, array=v, cmap="RdBu_r", edgecolor="#555", linewidths=0.1)
            pc.set_clim(-vmax, vmax)
            ax.add_collection(pc)
            ax.autoscale_view()
            ax.invert_yaxis()
            ax.set_aspect("equal")
            ax.grid(False)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_title(f"{cp['modes'][im]}  {lab} ΔCp\nk = {k[ik]:.3f}", fontsize=8)
            fig.colorbar(pc, ax=ax, fraction=0.05, pad=0.02, shrink=0.8).ax.tick_params(labelsize=6)
    fig.tight_layout()
    p = FIG / f"cp_{name}.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_flutter_compare(run, name, vmax=None):
    fl = {m: load(f"{run}/flutter_{m}.txt") for m in ("pk", "k", "p")}
    fl = {m: d for m, d in fl.items() if d}
    if not fl:
        return None
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.6, 5.6), sharex=True)
    styles = {"pk": dict(lw=1.6, ls="-"), "k": dict(lw=0, marker="o", ms=2.2), "p": dict(lw=1.0, ls="--")}
    for m, d in fl.items():
        V = np.array(d["V"], dtype=float)
        g = np.array(d["g"], dtype=float)
        f = np.array(d["freq"], dtype=float)
        for b in range(g.shape[1]):
            c = SERIES[b % len(SERIES)]
            lab = f"{d['names'][b]}" if m == "pk" else None
            Vb = V[:, b] if V.ndim == 2 else V
            a1.plot(Vb, g[:, b], color=c, label=lab, **styles[m])
            a2.plot(Vb, f[:, b], color=c, **styles[m])
        for pnt in d["points"]:
            if pnt["kind"] == "flutter" and pnt["velocity"]:
                a1.plot(pnt["velocity"], 0, marker="*", ms=11, color=RED, mec="white", zorder=5)
    ref = fl.get("pk") or next(iter(fl.values()))
    vr = np.array(ref["V"], dtype=float).ravel()
    fr = np.array(ref["freq"], dtype=float)
    a2.set_xlim(vr.min(), vmax or vr.max())
    a2.set_ylim(0, 1.25 * np.nanmax(fr))
    a1.axhline(0, color="#333", lw=0.8)
    a1.set_ylim(-0.5, 0.25)
    a1.set_ylabel("damping g")
    a2.set_ylabel("frequency [Hz]")
    a2.set_xlabel("V [m/s]")
    h1, = a1.plot([], [], color="#444", lw=1.6, label="PK (line)")
    h2, = a1.plot([], [], color="#444", lw=0, marker="o", ms=3, label="K (dots)")
    h3, = a1.plot([], [], color="#444", lw=1, ls="--", label="p (dashed)")
    a1.legend(ncol=4, loc="lower left")
    fig.tight_layout()
    p = FIG / f"vg_{name}.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_tsd_steady(run="tsd_agard_steady"):
    d = load(f"{run}/tsd_steady.txt")
    if d is None:
        return None
    x = np.array(d["x"])
    y = np.array(d["y"])
    cp = np.array(d["cp"])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 3.4), gridspec_kw={"width_ratios": [1.15, 1]})
    Y = np.repeat(y[:, None], x.shape[1], axis=1)
    cs = a1.tricontourf(Y.ravel(), x.ravel(), cp.ravel(), levels=24, cmap="RdBu")
    a1.tricontour(Y.ravel(), x.ravel(), cp.ravel(), levels=[d["cp_star"]], colors=[INK], linewidths=1.0)
    a1.invert_yaxis()
    a1.set_aspect("equal")
    a1.set_xlabel("y [m]")
    a1.set_ylabel("x [m]")
    a1.set_title(f"upper-surface Cp, M = {d['mach']}  (black: Cp*)")
    a1.grid(False)
    fig.colorbar(cs, ax=a1, shrink=0.85)
    for j in np.linspace(0, len(y) - 1, 4).astype(int):
        xc = (x[j] - x[j].min()) / max(np.ptp(x[j]), 1e-12)
        a2.plot(xc, cp[j], label=f"y = {y[j]:.2f} m")
    a2.axhline(d["cp_star"], color=RED, ls="--", lw=0.9, label="Cp* (sonic)")
    a2.invert_yaxis()
    a2.set_xlabel("x / c")
    a2.set_ylabel("Cp")
    a2.legend()
    a2.set_title("chordwise Cp (NACA 65A004, alpha = 0)")
    fig.tight_layout()
    p = FIG / "tsd_steady.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_timehist(run, fname, keys, name, title=None, ncol=1):
    d = load(f"{run}/{fname}")
    if d is None:
        return None
    t = np.array(d["t"], dtype=float)
    keys = [k for k in keys if k in d["series"]]
    n = len(keys)
    fig, axs = plt.subplots(n, 1, figsize=(7.6, 1.35 * n + 0.6), sharex=True, squeeze=False)
    for i, k in enumerate(keys):
        a = axs[i, 0]
        a.plot(t, np.array(d["series"][k], dtype=float), color=SERIES[i % len(SERIES)], lw=1.2)
        a.set_ylabel(k, fontsize=7.5)
    axs[-1, 0].set_xlabel("t [s]")
    if title:
        axs[0, 0].set_title(title)
    fig.tight_layout()
    p = FIG / f"th_{name}.png"
    fig.savefig(p)
    plt.close(fig)
    return p


def fig_curves(run, fname, keys, name, title=None, normalise=False):
    d = load(f"{run}/{fname}")
    if d is None:
        return None
    x = np.array(d["x"], dtype=float)
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    for i, k in enumerate(keys or list(d["curves"])):
        if k not in d["curves"]:
            continue
        y = np.array(d["curves"][k], dtype=float)
        if normalise:
            y = y / (np.nanmax(np.abs(y)) or 1)
        ax.plot(x, y, label=k, color=SERIES[i % len(SERIES)], lw=1.3, marker="o" if x.size < 30 else None, ms=3)
    if d.get("logy"):
        ax.set_yscale("log")
        ax.set_xscale("log")
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(matplotlib.ticker.LogFormatter(labelOnlyBase=True))
    ax.set_xlabel(d.get("xlabel", ""))
    ax.set_ylabel(d.get("ylabel", "") + (" (normalised)" if normalise else ""))
    if title:
        ax.set_title(title)
    ax.legend(ncol=2)
    fig.tight_layout()
    p = FIG / f"cv_{name}.png"
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- page templates
class Doc(BaseDocTemplate):
    def __init__(self, path):
        super().__init__(str(path), pagesize=A4, leftMargin=MARGIN, rightMargin=MARGIN, topMargin=22 * mm,
                         bottomMargin=18 * mm, title="FLAAP 기술 백서", author="FLAAP project",
                         subject="Aeroelastic analysis: capabilities, theory, validation")
        frame = Frame(MARGIN, 18 * mm, TEXT_W, PAGE_H - 40 * mm, id="f")
        cover = Frame(0, 0, PAGE_W, PAGE_H, id="c", leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate("cover", [cover], onPage=self._cover),
                               PageTemplate("normal", [frame], onPage=self._normal)])

    def _cover(self, c, doc):
        c.saveState()
        c.setFillColor(colors.HexColor("#122820"))
        c.rect(0, 0, PAGE_W, PAGE_H, stroke=0, fill=1)
        c.setFillColor(colors.HexColor(ACID))
        c.rect(MARGIN, PAGE_H - 62 * mm, 26 * mm, 2.2 * mm, stroke=0, fill=1)
        c.setFillColor(colors.white)
        c.setFont("MalgunBd", 40)
        c.drawString(MARGIN, PAGE_H - 85 * mm, "FLAAP")
        c.setFont("MalgunBd", 20)
        c.drawString(MARGIN, PAGE_H - 100 * mm, "공탄성 해석 프로그램 기술 백서")
        c.setFont("Malgun", 11.5)
        c.setFillColor(colors.HexColor("#c9d6cf"))
        c.drawString(MARGIN, PAGE_H - 112 * mm, "Flutter, Loads & Aeroelastic Analysis Program (Python edition)")
        c.drawString(MARGIN, PAGE_H - 119 * mm, "기능 · 이론 · 검증 결과")
        img = FIG / "cover_modes.png"
        if img.exists():
            c.drawImage(str(img), MARGIN - 4 * mm, 52 * mm, width=PAGE_W - 2 * MARGIN + 8 * mm,
                        height=(PAGE_W - 2 * MARGIN + 8 * mm) * 0.42, mask="auto")
        c.setFont("Malgun", 9)
        c.setFillColor(colors.HexColor("#9fb3aa"))
        c.drawString(MARGIN, 30 * mm, "DLM / VLM · TSD · PK / K / p flutter · Roger RFA · gust · trim · manoeuvre")
        c.drawString(MARGIN, 24 * mm, "control-surface freeplay · aeroservoelasticity · Nastran / ANSYS / OpenRadioss / CalculiX / STL")
        c.drawString(MARGIN, 16 * mm, time.strftime("%Y년 %m월"))
        c.restoreState()

    def _normal(self, c, doc):
        c.saveState()
        c.setStrokeColor(colors.HexColor("#cfcabd"))
        c.setLineWidth(0.5)
        c.line(MARGIN, PAGE_H - 15 * mm, PAGE_W - MARGIN, PAGE_H - 15 * mm)
        c.setFont("MalgunBd", 8)
        c.setFillColor(colors.HexColor(GREEN))
        c.drawString(MARGIN, PAGE_H - 13 * mm, "FLAAP")
        c.setFont("Malgun", 8)
        c.setFillColor(colors.HexColor(GREY))
        c.drawString(MARGIN + 13 * mm, PAGE_H - 13 * mm, "공탄성 해석 프로그램 기술 백서")
        c.drawRightString(PAGE_W - MARGIN, 10 * mm, f"{doc.page}")
        c.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph):
            st = flowable.style.name
            if st in ("h1", "h2"):
                lvl = 0 if st == "h1" else 1
                txt = flowable.getPlainText()
                key = f"h{id(flowable)}"
                self.canv.bookmarkPage(key)
                self.canv.addOutlineEntry(txt, key, level=lvl, closed=lvl > 0)
                self.notify("TOCEntry", (lvl, txt, self.page, key))


def H1(t):
    return [CondPageBreak(60 * mm), Paragraph(t, S["h1"])]


def H2(t):
    return [CondPageBreak(70 * mm), Paragraph(t, S["h2"])]


def H3(t):
    return [CondPageBreak(75 * mm), Paragraph(t, S["h3"])]


# --------------------------------------------------------------------------- content
def runs_summary():
    p = OUT / "run_all_summary.txt"
    if not p.exists():
        return {}
    out = {}
    for r in read_result(p).get("runs") or []:
        r = dict(r)
        rt = Path(r["output_dir"]) / "results.txt" if r.get("output_dir") else None
        r["summary"] = read_result(rt, "results").get("summary", {}) if rt and rt.exists() else {}
        out[r["config"]] = r
    return out


def sec_summary(R, md):
    s = []
    s += H1("요약")
    s.append(P(
        "FLAAP(Flutter, Loads &amp; Aeroelastic Analysis Program)는 원래 Fortran으로 작성된 아음속 공탄성 해석 "
        "코드를 Python으로 완전히 재구현하고, 천음속·초음속 공력, VLM/UVLM 시간영역 공력, 시간영역 하중, 정적 공탄성, 조종면 비선형성, "
        "공력-서보-탄성(ASE), 상용 FE 결과 연동 및 형상 기반 모달 해석까지 확장한 통합 공탄성 해석 "
        "환경이다. 이 문서는 (1) FLAAP이 현재 수행하는 기능, (2) 각 기능의 이론적 근거, (3) 풍동 시험, 독립 해석, 해석해에 대한 "
        "검증 결과의 순서로 구성되며 검증 결과로 끝난다."))
    s.append(P("<b>FLAAP으로 할 수 있는 일</b>", "h3"))
    s += table([
        ["활용 목적", "사용 기능 (model / mode)", "얻는 결과"],
        ["플러터·발산 속도 예측<br/>(아음속)", "flutter / pk, k, p (DLM)", "플러터 속도·진동수, V-g / V-f 곡선, 발산 속도"],
        ["천음속 플러터 경계와 dip", "tsd / steady, gaf, flutter", "충격파 위치·CL이 포함된 정상 해, 천음속 플러터 속도"],
        ["초음속 플러터", "flutter (Mach box, 피스톤 이론)", "M &gt; 1 플러터 속도, 패널 플러터"],
        ["저속 하중과 와류 후류 효과", "uvlm / steady, impulsive, harmonic, gust, flutter", "CL, 유도항력, 스팬 하중, 비선형 시간응답, UVLM 플러터"],
        ["돌풍 하중 (CS-25.341)", "gust / discrete, continuous", "1-cos 돌풍 포락선, 임계 돌풍 길이, 난류 응답 Ā, N<sub>0</sub>"],
        ["기동 하중 (CS-25.331(c))", "maneuver / checked, abrupt", "하중배수, 조종 입력, 단면 하중 시간응답"],
        ["1-g 트림, 조종 효율·역전", "trim / static, divergence", "트림 자세와 조종면 각, 발산 동압, 역전 속도"],
        ["조종면 freeplay 허용 한도", "freeplay / describing_function, time, bifurcation", "LCO 진폭-속도 곡선, 이력(hysteresis), 시간 응답"],
        ["플러터 억제·돌풍 완화 제어 설계", "ase / flutter_suppression, margins, gla", "폐루프 플러터 속도, Nyquist / Bode 안정여유, 하중 저감"],
        ["실제 형상 모델 해석", "modal + structure = nastran, ansys, openradioss, calculix, stl, gmsh",
         "외부 FE 결과 또는 형상 + 물성에서 구조 모달 모델 구성, 이후 모든 해석에 연결"],
        ["원본 Fortran 결과 재현", "python -m flaap legacy", "원본 입력덱 실행, 원본과 같은 형식의 출력"],
    ], [42 * mm, 62 * mm, 66 * mm], caption="활용 목적별 FLAAP 기능 (이론은 3장, 입력·출력 키는 2.4, 2.6절)", font=7.8)
    s.append(P("<b>현재 다루지 않는 것.</b> 플러터 이후의 LCO는 조종면 힌지 비선형(freeplay, bilinear, preloaded)에 한정된다. "
               "공력 비선형(충격파 진동, 실속)이나 날개 구조의 기하 비선형에 의한 LCO, 점성 효과는 아직 포함하지 않는다(2.9절)."))
    s.append(CondPageBreak(110 * mm))
    s.append(P("<b>검증 요약</b> (검증 대상 기능별, 비교 기준 우선순위: 실험 &gt; 타 해석 &gt; 해석해)", "h3"))
    rows = val_summary_rows(md)
    if rows:
        s += table(rows, [9 * mm, 42 * mm, 24 * mm, 27 * mm, 24 * mm, 44 * mm], font=7.0,
                   caption="검증 결과 요약 (상세는 4장)")
    return s


def sec_intro():
    s = []
    s += H1("1. 개요")
    s.append(P(
        "공탄성(aeroelasticity)은 구조의 탄성 변형, 관성력, 공기력이 상호작용하는 현상이다. 항공기 설계에서는 "
        "플러터(동적 불안정), 발산(정적 불안정), 돌풍과 기동에 의한 동하중, 조종면 효율 저하, 그리고 제어계가 "
        "개입하는 공력-서보-탄성(ASE) 문제를 인증 단계(CS/FAR 25.629, 25.341, 25.331)에서 모두 다뤄야 한다."))
    s.append(P(
        "원본 FLAAP은 DLM 기반 일반화 공기력 계산과 PK 플러터 해석을 수행하는 Fortran 프로그램이었으며, 일부 "
        "루틴은 목적 파일만 남아 있었다. Python 판은 다음을 목표로 재설계되었다."))
    s += bullets([
        "원본 입력덱(LEE_input.txt, M_K.txt, MODESHAPE_input.txt)을 그대로 읽어 동일한 결과를 재현할 것",
        "공력(DLM/TSD/UVLM/초음속) - 구조(FE/외부 모달) - 해석(플러터/하중/트림/ASE)을 독립 모듈로 분리하고, 동일한 "
        "일반화 공기력 표 Q(k)를 모든 해석이 공유할 것",
        "평문 설정 파일 하나로 해석을 정의하고 사전검증(preflight) 후 실행하며, 결과는 입력과 같은 key-value 문법의 TXT 계약으로 "
        "출력하여 GUI와 보고서가 같은 데이터를 쓰도록 할 것",
        "모든 기능을 공개 벤치마크로 검증하고, 그 과정을 재실행 가능한 스크립트로 남길 것",
    ])
    s.append(P("<b>문서 구성.</b> 2장은 현재 기능(해석 경로, 공력·구조 모델, 입출력, GUI, 예제, 적용 범위), 3장은 각 기능의 이론, "
               "4장은 검증 결과이다. 4장은 기능별로 비교 기준(실험, 타 해석, 해석해)과 자기일관성 검사를 구분하여 보인다."))
    s.append(P("좌표계와 단위: x 하류, y 우측(스팬), z 위, SI 단위. 조화 운동은 e<super>iωt</super>, 무차원 "
               "진동수 k = ωb/V (b: 기준 반코드). 조종면 변위는 뒷전 아래(TE down)가 양이다."))
    return s


def sec_func_overview():
    s = []
    s += H1("2. 기능")
    s.append(P(
        "이 장은 FLAAP이 현재 수행할 수 있는 해석을 정리한다. 모든 해석은 하나의 평문 설정 파일의 <font face='Consolas'>model</font>"
        "(해석 경로)과 <font face='Consolas'>mode</font>(경로 내 연산)로 선택하며, 구조 모델, 공력 모델, 해석 모듈은 "
        "서로 독립적으로 조합된다. 이론적 배경은 3장, 정확도 검증은 4장에 있다."))
    s += H2("2.1 해석 경로와 연산")
    s += table([
        ["model", "mode", "수행 내용 / 주요 출력"],
        ["modal", "analysis", "구조 모달 모델 구성: 진동수, 일반화 질량, 모드형상(프리셋, 판·쉘·고체 FE, Nastran / ANSYS / OpenRadioss / "
                              "CalculiX 입력, STL / STEP / gmsh 형상 + 물성)"],
        ["aero", "gaf, steady", "비정상 GAF 표 Q(k)(모드별, 돌풍 열, 조종면)와 ΔCp 분포, k = 0 정상 압력"],
        ["flutter", "pk, k, p, all", "PK / K(V-g) / p 방법 플러터. 플러터·발산 자동 식별, V-g / V-f 곡선, 원본 형식 출력"],
        ["gust", "discrete, continuous", "CS-25.341 1-cos 돌풍 포락선과 임계 돌풍 길이, von Kármán 난류 PSD, Ā, N<sub>0</sub>"],
        ["trim", "static, divergence", "1-g 정적 공탄성 트림, 발산 동압, 조종 효율 / 역전, 단면 하중"],
        ["maneuver", "checked, abrupt", "CS-25.331(c) 체크 / 급격 피치 기동: 하중배수, 조종 입력, 단면 하중 시간응답"],
        ["tsd", "steady, gaf, flutter", "3-D 천음속 TSD: 두께·캠버·받음각 비선형 평균유동(충격파 위치, CL, 받음각 스윕), 시간선형화 GAF, "
                                        "TSD 기반 PK 플러터(받음각 평균유동 포함)"],
        ["uvlm", "steady, impulsive, harmonic, gust, flutter", "정상 VLM(CL, 유도항력, 스팬하중, 공력중심), 비선형 UVLM 시간영역 "
                                                              "(급출발, 조화 피치·플런지, 돌풍, 자유후류), 선형 UVLM GAF + PK / 이산시간 플러터"],
        ["freeplay", "describing_function, time, bifurcation", "조종면 freeplay LCO: 기술함수 진폭-속도 곡선, 시간적분, 분기·이력 스윕"],
        ["ase", "flutter_suppression, margins, gla", "작동기·센서·정적 출력 피드백 폐루프: 플러터 억제, Nyquist / Bode 안정여유, "
                                                    "돌풍 하중 완화"],
    ], [20 * mm, 38 * mm, 112 * mm], caption="해석 경로(model)와 연산(mode)", font=7.8)
    s += H3("공력 모델 (aero_method)")
    s += table([
        ["aero_method", "적용 범위", "특징"],
        ["dlm (기본)", "아음속 선형 M &lt; 1", "Albano-Rodden 이중격자법, 조종면·돌풍 열 포함"],
        ["tsd", "천음속 (M ≈ 0.9-1.0)", "비선형 평균유동 위 시간선형화, 천음속 dip 재현, 받음각·캠버 평균유동 선택"],
        ["uvlm", "비압축성 M = 0", "선형화 UVLM 상태공간(GAF + 이산시간 플러터)"],
        ["machbox / piston", "초음속 M &gt; 1", "Mach box(Evvard 끝단 효과 포함) / 피스톤 이론(1·3차, Van Dyke)"],
        ["(Python API) strip", "2-D 극한", "Theodorsen / Sears 띠이론: 해석해 비교와 단면 모델"],
    ], [34 * mm, 42 * mm, 94 * mm], caption="공력 모델 선택", font=7.8)
    s += H3("구조 입력 (structure)")
    s += table([
        ["structure", "내용"],
        ["goland, goland_ctrl, agard445, demo_aircraft", "내장 모델: Goland 보(조종면 포함 선택), AGARD 445.6 약화 모델 3, 날개 + 꼬리 데모 항공기"],
        ["plate", "MITC4 직교이방성 판 사다리꼴 날개(요소 수, 두께, 재료 지정)"],
        ["nastran, ansys", "SOL 103 .f06 고유벡터 + GRID / ANSYS 모달 listing, PRNSOL, Workbench 텍스트, .cdb, .rst"],
        ["openradioss, openradioss_eig", "입력덱(/NODE, /SHELL, /SH3N, /BRICK, /TETRA)에서 FE 재구성 후 고유치 해석 / Radioss /EIG 결과(VTK)"],
        ["calculix", "*FREQUENCY .frd(short / long ASCII) + .dat 고유치 표"],
        ["stl, geometry, gmsh", "STL / STEP 형상 + 물성 → 쉘 또는 tet10 / gmsh .geo, .msh, CAD → quad4·tri3·tet10(Physical 그룹 구속, 두께)"],
        ["legacy", "원본 Fortran 입력덱(LEE_input.txt, M_K.txt, MODESHAPE_input.txt, atmosphere_input.txt)"],
    ], [58 * mm, 112 * mm], caption="구조 모델 입력 경로", font=7.8)
    s.append(P(
        "구조-공력 연성은 스플라인(ips, tps, rbf, beam, legacy, auto)으로 선택하고, 상태공간 해석에는 Roger RFA(지연항 수 nlag)를 쓴다. "
        "파이썬 API(<font face='Consolas'>flaap.*</font>)로 같은 기능을 직접 호출할 수 있으며, "
        "<font face='Consolas'>python -m flaap</font> 명령줄(run, legacy, mesh, f06, validate)과 "
        "FLAAP Studio GUI도 같은 실행 경로를 쓴다."))
    return s

def _esc(x):
    return str(x).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def sec_input_keys():
    """2.4 input keys, generated from flaap.launcher.specs (single source of truth)."""
    from flaap.launcher.specs import BASE_DEFAULTS, COMMON_KEYS, INDEXED_PATTERNS, ROUTES
    from flaap.launcher.parser import format_value
    s = []
    s += H2("2.4 입력 키 설명")
    s.append(P(
        "입력 키는 launcher의 단일 정의(<font face='Consolas'>flaap/launcher/specs.py</font>)에서 자동으로 가져온 것이며 "
        "전체 목록은 docs/CONFIGURATION.md와 <font face='Consolas'>--describe &lt;model&gt;</font>으로도 볼 수 있다. "
        "아래 표는 각 키를 처음 정의하는 경로 순서(modal, aero, flutter, ...)로 한 번씩만 싣고, "
        "'사용 경로' 열에 그 키가 유효한 경로를 적었다. 키 이름과 값은 소문자로 읽히며 경로 키(structure_file 등)만 대소문자를 유지한다. "
        "목록 값은 쉼표로 구분한다(예: <font face='Consolas'>velocities 60, 260, 101</font>은 시작, 끝, 개수). "
        "키 설명 문구는 specs.py의 원문(영문)을 그대로 쓴다."))
    s += H3("모든 경로 공통 키")
    kor = {
        "model": "해석 경로 선택: modal, aero, flutter, gust, trim, maneuver, tsd, uvlm, freeplay, ase (필수)",
        "mode": "경로 안의 연산(생략하면 그 경로의 첫 번째 mode)",
        "title": "케이스 제목(대소문자 유지, results.txt에 기록)",
        "output_dir": "출력 폴더(기본 output/&lt;설정 파일 이름&gt;)",
        "plots": "PNG 그림 생성 여부(기본 true)",
        "verbose": "진행 상황 출력 여부(기본 true)",
        "output_json": "모든 결과 파일을 JSON으로도 저장(기본은 텍스트만)",
    }
    rows = [["키", "유형", "설명"]]
    for k, (desc, typ) in COMMON_KEYS.items():
        rows.append([f"<font face='Consolas'>{k}</font>", typ, kor.get(k, _esc(desc))])
    s += table(rows, [30 * mm, 14 * mm, 126 * mm], font=7.6, caption="모든 경로 공통 키")
    s += H3("인덱스 키: 양력면, 조종면, 하중 단면")
    idx_kor = {
        "surface": "양력면 정의. 키 이름 surface1, surface2, ... 값은 이름, 뿌리 앞전 좌표(rlex, rley, rlez), 뿌리 코드, "
                   "끝 앞전 좌표(tlex, tley, tlez), 끝 코드, 코드 방향 분할 수(nchord), 스팬 방향 분할 수(nspan)",
        "control": "조종면 정의. 이름, 소속 양력면 이름, 힌지선 위치(코드 비 0-1), 스팬 시작과 끝 비율(eta0, eta1)",
        "station": "하중 회수 단면. 이름, 단면 위치(px, py, pz), 단면 법선 축(ax, ay, az)",
    }
    rows = [["키 패턴", "값 (순서대로)", "설명"]]
    for k, (_, fields) in INDEXED_PATTERNS.items():
        rows.append([f"<font face='Consolas'>{k}N</font>", _esc(fields), idx_kor.get(k, "")])
    s += table(rows, [24 * mm, 56 * mm, 90 * mm], font=7.4, caption="인덱스 키(surfaceN, controlN, stationN)")

    order = ["modal", "aero", "flutter", "gust", "trim", "maneuver", "tsd", "uvlm"]
    keys_of, defs_of = {}, {}
    for n in order:
        sp = ROUTES[n]
        ks = {}
        for g in sp.key_groups:
            ks.update(g)
        for ck in COMMON_KEYS:
            ks.pop(ck, None)
        keys_of[n] = ks
        d = dict(BASE_DEFAULTS)
        d.update(sp.defaults)
        defs_of[n] = d
    seen = set()
    titles = {
        "modal": "구조 모델과 공통 해석 키 (modal 경로에서 정의, 대부분의 경로가 공유)",
        "aero": "공력 모델과 GAF 키 (aero 경로에서 처음 정의)",
        "flutter": "플러터 해석 키",
        "gust": "돌풍 해석 키",
        "trim": "트림, 발산 키",
        "maneuver": "기동 해석 키",
        "tsd": "천음속 TSD 키",
        "uvlm": "VLM / UVLM 키",
    }
    for n in order:
        new = [k for k in keys_of[n] if k not in seen]
        if not new:
            continue
        s += H3(titles[n])
        rows = [["키", "유형", "기본값", "사용 경로", "설명"]]
        for k in new:
            desc, typ = keys_of[n][k]
            users = [r for r in order if k in keys_of[r]]
            use = "전체(modal-uvlm)" if len(users) == len(order) else ", ".join(users)
            dv = [defs_of[u].get(k) for u in users]
            dv = [format_value(v) for v in dv if v is not None]
            dtxt = "" if not dv else (dv[0] if len(set(dv)) == 1 else "경로별 상이")
            rows.append([f"<font face='Consolas'>{k}</font>", typ, _esc(dtxt), use, _esc(desc)])
        seen.update(new)
        s += table(rows, [27 * mm, 12 * mm, 20 * mm, 24 * mm, 87 * mm], font=6.9,
                   caption=titles[n] + f" ({len(new)}개)", keep=False)
    for n, ttl in (("freeplay", "조종면 freeplay LCO 키 (freeplay 경로 전용 플랜트)"),
                   ("ase", "공력-서보-탄성 키 (ase 경로 전용 플랜트)")):
        sp = ROUTES[n]
        d = dict(BASE_DEFAULTS)
        d.update(sp.defaults)
        ks = {}
        for g in sp.key_groups:
            ks.update(g)
        for ck in COMMON_KEYS:
            ks.pop(ck, None)
        rows = [["키", "유형", "기본값", "설명"]]
        for k, (desc, typ) in ks.items():
            dv = d.get(k)
            rows.append([f"<font face='Consolas'>{k}</font>", typ, _esc(format_value(dv)) if dv is not None else "", _esc(desc)])
        s += H3(ttl)
        s += table(rows, [30 * mm, 12 * mm, 22 * mm, 106 * mm], font=7.0, caption=ttl)
    rows = [["model", "mode", "필수 키"]]
    for n, sp in ROUTES.items():
        for m, r in sp.required_by_mode.items():
            if r:
                rows.append([n, m, _esc(", ".join(r))])
    if len(rows) > 1:
        s += H3("mode별 필수 키")
        s += table(rows, [26 * mm, 40 * mm, 104 * mm], font=7.6, caption="경로, mode별로 반드시 지정해야 하는 키")
    return s


def sec_output_keys():
    """2.6 output files and their keys."""
    s = []
    s += H2("2.6 출력 키 설명")
    s.append(P(
        "모든 출력 파일은 2.5절의 문법으로 쓰이며, 아래 표의 키는 <font face='Consolas'>flaap.io.results_txt.read_result</font>로 "
        "읽는다. 각 파일의 첫 줄 <font face='Consolas'>kind</font> 값이 파일 종류이고, "
        "뷰어(Studio, 본 백서)는 results.txt의 artifacts 표에 적힌 kind로 표시 방법을 정한다."))
    s += H3("출력 폴더에 항상 생기는 파일")
    s += table([
        ["파일", "kind", "내용"],
        ["results.txt", "results", "실행 기록: 상태, 요약값(summary.*), 산출물 목록(artifacts 표)"],
        ["config_resolved.txt", "config", "기본값을 채운 전체 설정(그대로 다시 실행 가능)"],
        ["log.txt", "-", "솔버 진행 로그"],
        ["figures/*.png", "figure", "matplotlib 그림(plots false이면 생략)"],
    ], [38 * mm, 22 * mm, 110 * mm], font=7.8, caption="항상 생성되는 파일")
    s += H3("results.txt의 키")
    s += table([
        ["키", "설명"],
        ["<font face='Consolas'>kind</font>", "항상 results"],
        ["<font face='Consolas'>title, model, mode</font>", "케이스 제목과 실행한 경로, mode"],
        ["<font face='Consolas'>status</font>", "ok 또는 error. error이면 error 키에 트레이스백"],
        ["<font face='Consolas'>elapsed_s, created, config</font>", "계산 시간 [s], 생성 시각, 설정 파일 경로"],
        ["<font face='Consolas'>summary.*</font>", "경로별 요약값(아래 표). 중첩 값은 점으로 이은 키(summary.flutter_pk.velocity)"],
        ["<font face='Consolas'>table artifacts</font>", "산출물 목록: file, kind, title. kind는 config, modes, aero_mesh, cp, gaf, flutter, "
                                                        "timehist, curves, tsd_steady, table, text, figure"],
    ], [58 * mm, 112 * mm], font=7.8, caption="results.txt 키")
    s += H3("해석 결과 파일과 주요 키")
    s += table([
        ["파일 (kind)", "키 / 표", "설명"],
        ["modes.txt (modes)", "name, structure, n_modes, n_nodes", "모델 이름, 구조 입력 종류와 크기"],
        ["", "table modes", "mode, label, freq_hz, omega_rad_s, gen_mass: 진동수와 일반화 질량"],
        ["", "table nodes / edges / quads / tris", "구조 절점 좌표 x y z와 연결 정보(1부터 시작하는 절점 번호)"],
        ["", "table shape.N", "N번째 모드의 절점 변위 ux uy uz(최대값 1로 정규화)"],
        ["", "table aero_h.N", "N번째 모드의 공력 박스 네 모서리의 법선 변위 h1-h4(스플라인 결과)"],
        ["aero_mesh.txt (aero_mesh)", "surfaces, symmetry, bref, mach, controls, n_boxes", "공력 격자 요약"],
        ["", "table boxes", "box, surface, area, 배치점 xc yc zc, 네 모서리 x1..z4, ctrl:&lt;이름&gt;(조종면 소속 여부)"],
        ["cp.txt (cp)", "mach, bref, modes, k, n_boxes", "압력 분포 조건. modes 마지막 값이 gust이면 돌풍 열 포함"],
        ["", "table cp.N", "N번째 환산진동수의 박스별 ΔCp(re:&lt;모드&gt;, im:&lt;모드&gt;). 일반화 좌표 단위당 압력차"],
        ["gaf_table.txt (gaf)", "k, Q, Qg, mach, bref, names", "복소 GAF Q(k)(F = qQη)와 돌풍 열. aero_method와 무관한 같은 배치"],
        ["flutter_&lt;pk|k|p&gt;.txt (flutter)", "method, density, bref, g_struct, names, n_velocities", "해석 방법, 밀도, 기준 길이, 구조 감쇠, 가지 이름"],
        ["", "table points", "플러터, 발산 교차점: kind, branch, velocity, freq_hz, k, density"],
        ["", "table branch.N", "N번째 가지의 V, g(= 2ζ), freq_hz, k"],
        ["timehist_&lt;이름&gt;.txt (timehist)", "meta.*, units.*, series", "사례 정보(돌풍 길이 H, 속도 등), 단위, 계열 이름"],
        ["", "table timehist", "t와 각 계열의 시간 응답"],
        ["curves_&lt;이름&gt;.txt (curves)", "xlabel, ylabel, logy, meta.*, curves", "포락선, PSD, LCO 진폭-속도, 안정여유 곡선 등의 축 정보"],
        ["", "table curves", "x와 각 곡선의 값"],
        ["tsd_steady.txt (tsd_steady)", "mach, cp_star, converged, residuals, ny, nx", "TSD 정상 해의 조건, 음속 Cp*, 수렴 이력"],
        ["", "table grid", "j, i, y, x, cp: 윗면 Cp 분포"],
        ["vlm_steady.txt (vlm_steady)", "alpha_deg, CL, CDi, Cm, CL_alpha, span_efficiency, x_ac", "정상 VLM 극곡선과 양력 기울기, 스팬 효율, 공기역학 중심"],
        ["", "y_strip, cl_c, cl_local", "스팬 하중(스트립 y, 국소 cl c / c_ref, 국소 cl)"],
        ["uvlm_wake.txt (uvlm_wake)", "shape, surfaceN, free_wake", "UVLM 최종 후류 절점 좌표와 자유후류 여부"],
    ], [38 * mm, 54 * mm, 78 * mm], font=7.0, caption="결과 파일별 키", keep=False)
    s += H3("results.txt의 summary 키 (경로, mode별)")
    s += table([
        ["경로 / mode", "summary 키"],
        ["modal", "frequencies_hz, nmodes, n_boxes"],
        ["aero", "nk, mach, bref, Q0_diag"],
        ["flutter", "density, flutter_&lt;방법&gt;(velocity, frequency_hz, kind 등 임계점)"],
        ["gust discrete", "critical_H(임계 돌풍 길이), 출력별 peak"],
        ["gust continuous", "출력별 A_bar, N0"],
        ["trim static", "alpha_deg, &lt;조종면&gt;_deg, loads"],
        ["trim divergence", "q_divergence, V_divergence"],
        ["maneuver", "short_period_wn, short_period_zeta, peak_control_deg, peak_nz, peak_loads"],
        ["tsd steady", "converged, iterations, n_supersonic, cp_min, cp_star"],
        ["freeplay describing_function", "V_f_linear, V_f0, V_fold, time_marching(속도별 LCO 판정과 진폭 비교)"],
        ["freeplay time", "lco(V, kind: decay / lco / growing, amplitude_deg, mean_deg, freq_hz, growth, switches)"],
        ["freeplay bifurcation", "bifurcation(first_lco_up, kinds_up, kinds_down). 진폭은 curves_bifurcation.txt"],
        ["ase flutter_suppression", "flutter_open, flutter_closed, V_f_open, V_f_closed, increase"],
        ["ase margins", "margins(V, gm_lower_db, gm_upper_db, pm_deg, nyquist_stable 등), margins_vs_V"],
        ["ase gla", "gla(V, Uds, 행별 H, Mx_open, Mx_closed, Mx_reduction, flap_peak_deg), V_f_open, V_f_closed"],
        ["uvlm (모든 mode)", "s_ref, c_ref, x_ref, n_panels"],
        ["uvlm steady", "CL_alpha_per_rad, CL_alpha_per_deg, span_efficiency, aspect_ratio, x_ac, CL_max_alpha, CDi_max_alpha"],
        ["uvlm impulsive, harmonic, gust", "velocity, density, dt, nsteps, CL_final, CD_final, Cm_final"],
        ["uvlm impulsive", "CL_steady_vlm, CDi_steady_vlm, CL_final_over_steady"],
        ["uvlm harmonic", "CL_mean, CL_1st_harmonic_amplitude, CL_1st_harmonic_phase_deg, CD_mean_last_cycle, CT_mean_last_cycle"],
        ["uvlm gust", "dCL_gust_peak, t_peak, dL_gust_peak, dMy_gust_peak, CL_at_peak"],
        ["uvlm flutter", "density, flutter_pk, flutter_uvlm_state_space(velocity, freq_hz, mode)"],
    ], [50 * mm, 120 * mm], font=7.2, caption="summary 키 (자세한 정의는 docs/OUTPUT_FORMAT.md)", keep=False)
    s += H3("freeplay, ase 경로가 쓰는 파일")
    s += table([
        ["경로 / mode", "파일"],
        ["freeplay describing_function", "curves_lco_df.txt, curves_lco_tm_vs_df.txt, timehist_lco.txt, figures/lco_bifurcation.png"],
        ["freeplay time", "timehist_lco.txt(힌지, 피치, 플런지 응답), figures/phase_plane.png"],
        ["freeplay bifurcation", "curves_bifurcation.txt(속도 상승/하강 스윕의 진폭과 진동수)"],
        ["ase flutter_suppression", "flutter_open.txt, flutter_closed.txt, figures/ase_vg_rootlocus.png"],
        ["ase margins", "curves_bode.txt(|L|, 위상), figures/nyquist_bode.png"],
        ["ase gla", "timehist_gla.txt(개루프, 폐루프 뿌리 Mx, Fz, 조종면 각)"],
    ], [50 * mm, 120 * mm], font=7.6, caption="freeplay, ase 출력 파일")
    return s



def sec_arch(R):
    s = []
    s += H2("2.2 모듈 구성과 해석 흐름")
    s += figure(fig_architecture(), "FLAAP 해석 흐름. 모든 해석은 동일한 모달 모델과 GAF 표 Q(k)를 공유한다.")
    s += table([
        ["패키지 / 모듈", "역할"],
        ["flaap.structures", "beam, plate(MITC4), shell(quad4/tri3), solid(tet4/tet10), modal, spline"],
        ["flaap.io", "legacy 입력덱, Nastran .f06, ANSYS 모달, OpenRadioss 입력덱, Radioss /EIG 결과(anim_to_vtk), CalculiX .frd/.dat 판독, TXT 결과 입출력"],
        ["flaap.geometry / aero", "양력면 격자, VLM, DLM, GAF, Roger/최소상태 RFA, 띠이론, TSD, UVLM, 초음속(피스톤/Mach box)"],
        ["flaap.flutter", "K / PK / p 방법과 플러터·발산 식별, legacy 출력"],
        ["flaap.timedomain, gust, maneuver, loads", "상태공간, 돌풍, 기동, 단면하중(FSM/MDM)"],
        ["flaap.trim", "정적 트림, 발산, 조종 효율/역전"],
        ["flaap.nonlinear, flaap.ase", "조종면 freeplay LCO, 공력-서보-탄성"],
        ["flaap.launcher", "평문 설정 파서, 사전검증 진단, 라우트 실행, 결과 계약"],
        ["studio/", "FLAAP Studio: 브라우저 기반 CAD형 모델러와 결과 뷰어"],
    ], [52 * mm, 118 * mm], caption="주요 모듈")
    s += H2("2.3 입력: 평문 설정 파일")
    s.append(P(
        "AI-CAE4ALL 형식을 따라 한 줄에 <font face='Consolas'>key value</font> 하나를 쓴다. "
        "<font face='Consolas'>model</font>은 해석 경로(modal, aero, flutter, gust, trim, maneuver, tsd, freeplay, "
        "ase), <font face='Consolas'>mode</font>는 경로 내 연산을 고른다. % 는 주석 줄, # 은 행내 주석이며 쉼표로 "
        "목록을 쓴다. 양력면·조종면·하중 단면은 <font face='Consolas'>surfaceN / controlN / stationN</font> 인덱스 "
        "키로 정의한다."))
    s.append(code("""
% Goland wing flutter, DLM + PK (sea level)
model        flutter
mode         pk
structure    goland              # or plate / nastran / ansys / openradioss / stl / legacy
nmodes       4
mach         0.0
kmax         2.0
nk           20
density      1.225
velocities   60, 260, 101        # start, stop, count [m/s]
surface1     wing, 0,0,0, 0.5, 0.35,1.0,0, 0.3, 10, 14     # (custom planform)
control1     aileron, wing, 0.75, 0.55, 0.95              # hinge 75 % c, 55-95 % span
"""))
    s.append(P(
        "실행 전 사전검증이 파일 인코딩(BOM), 중복 키, 미지 키(철자 제안 포함), 형식, 값 범위, 파일 존재, "
        "경로별 필수 키, 상호 충돌(예: density와 altitude 동시 지정), TSD 적용 범위 등을 검사하고, "
        "CFG-xxx-nnn 코드와 수정 힌트를 출력한다. 오류가 있으면 해석을 시작하지 않는다(종료코드 1)."))
    s.append(code("""
python FLAAP_main.py --config configs/flutter/ex1/config_goland_pk.txt
python FLAAP_main.py --config <cfg> --check | --dry-run | --explain-config [--strict]
python FLAAP_main.py --list-models | --describe flutter | --audit-configs configs
python studio/start_studio.py          # GUI (http://127.0.0.1:8765)
"""))
    s += sec_input_keys()
    s += H2("2.5 출력 계약")
    s.append(P(
        "각 실행은 하나의 출력 폴더를 만든다. 결과는 입력 설정과 같은 문법(한 줄에 key value, % 주석)의 "
        "사람이 읽을 수 있는 텍스트 파일이다. results.txt(상태, 요약, 산출물 목록), config_resolved.txt(기본값 "
        "포함 전체 설정), modes.txt(진동수·모드형상·공력 격자 변위), aero_mesh.txt, cp.txt(k별·모드별 ΔCp), "
        "flutter_&lt;pk|k|p&gt;.txt(V-g/V-f), timehist_*.txt(시간응답), curves_*.txt(포락선·PSD 등), "
        "tsd_steady.txt 및 figures/*.png이다. 중첩 값은 점(.)으로 이은 키(summary.flutter_pk.velocity)로, "
        "배열은 'table 이름 ... end' 블록(첫 줄 열 이름, 이후 공백 구분 행)으로 쓴다. 누락 값은 none / nan이다. "
        "GUI, 테스트, 본 백서는 모두 같은 파서(flaap.io.results_txt.read_result)로 이 파일만 읽어 그림을 만든다. "
        "기존 JSON이 필요하면 설정에 output_json true를 추가한다."))
    s.append(code("""
kind          results
model         flutter
status        ok
summary.flutter_pk.velocity     153.29
summary.flutter_pk.frequency_hz 11.27

table branch.1                  % flutter_pk.txt: V-g / V-f of branch 1
        V          g   freq_hz         k
       60  -0.01523     7.712   0.2563
      ...
end
"""))
    s += sec_output_keys()
    shots = sorted((ROOT / "studio" / "screenshots").glob("*.png")) if (ROOT / "studio" / "screenshots").exists() else []
    if shots:
        s += H2("2.7 FLAAP Studio (GUI)")
        s.append(P(
            "FLAAP Studio는 표준 라이브러리 HTTP 서버와 의존성 없는 JavaScript로 구현된 로컬 GUI이다. 평면도 위에서 "
            "양력면과 조종면을 CAD처럼 그리고(격자 스냅, 핸들 편집, 수치 속성창), 설정 파일과 양방향으로 변환하며, "
            "사전검증과 실행 로그를 확인하고, 결과를 3-D 모드 애니메이션, ΔCp 분포, V-g/V-f, 시간응답 그래프로 본다."))
        cap = {"01": "모델러: 평면도 위 양력면/조종면 그리기와 속성 편집", "02": "플러터 결과 V-g / V-f",
               "03": "Goland 날개 3-D 모드형상 애니메이션", "04": "진동 압력분포 ΔCp", "05": "돌풍 시간응답",
               "06": "돌풍 포락선/PSD 곡선", "07": "TSD 정상해 Cp 분포", "13": "Nastran .f06 모드형상",
               "14": "기동 시간응답"}
        for sp in shots:
            if sp.stem[:2] in cap:
                s += figure(sp, f"FLAAP Studio - {cap[sp.stem[:2]]}", max_h=85 * mm)
    return s


def sec_theory_structure():
    s = []
    s += H1("3. 이론")
    s += H2("3.1 구조 동역학과 모달 축약")
    s.append(P("FE 이산화된 구조의 비감쇠 자유진동은 일반화 고유치 문제로 귀결되며, 하위 n개 모드로 변위를 "
               "근사하면 모달 운동방정식을 얻는다."))
    s.append(eq(r"\mathbf{K}\,\phi_i=\omega_i^2\,\mathbf{M}\,\phi_i,\qquad \mathbf{u}(t)\approx\mathbf{\Phi}\,\eta(t)"))
    s.append(eq(r"\mathbf{M}_{hh}\ddot\eta+\mathbf{C}_{hh}\dot\eta+\mathbf{K}_{hh}\eta=q_\infty\,\mathbf{Q}(k)\,\eta+\mathbf{f}_{ext},"
                r"\qquad \mathbf{M}_{hh}=\mathbf{\Phi}^T\mathbf{M}\mathbf{\Phi},\ \ \mathbf{K}_{hh}=\mathbf{\Phi}^T\mathbf{K}\mathbf{\Phi}"))
    s.append(P("모달 감쇠는 C<sub>hh</sub> = diag(2ζ<sub>i</sub> ω<sub>i</sub> m<sub>i</sub>)로 둔다. FLAAP은 다음 요소로 "
               "구조 모델을 직접 만들거나 외부 모달 결과를 읽는다."))
    s += table([
        ["요소 / 소스", "정식화", "용도"],
        ["보 (Goland)", "Hermite 3차 굽힘 + 선형 St-Venant 비틀림, 질량중심 편심 관성 결합", "막대형 날개, 해석해 검증"],
        ["MITC4 판", "Reissner-Mindlin, 혼합보간 전단(Bathe-Dvorkin)으로 잠김 제거, 직교이방성", "AGARD 445.6, 사다리꼴 판"],
        ["평판 쉘 quad4 / tri3", "MITC4 굽힘 + 막 + drilling, 요소 국부좌표 → 전역 6자유도 변환", "STL 중립면 쉘 모델"],
        ["고체 tet4 / tet10", "등매개 사면체, 일관 질량", "덩어리형 부품, STEP 고체"],
        ["Nastran / ANSYS / CalculiX", ".f06 SOL 103 고유벡터, ANSYS 모달 출력, CalculiX .frd 판독", "상용·오픈소스 FE 모달 연동"],
        ["OpenRadioss", "입력덱(/NODE, /SHELL, /BRICK, /TETRA, 물성)으로 FE 재구성 후 고유치 해석", "충돌해석 모델 재사용"],
        ["Radioss /EIG 결과", "anim_to_vtk 변환 VTK(Displacement, Rotational_DOF)의 모드형상 + 진동수(TIME 또는 엔진 리스팅), 선택적으로 입력덱 질량행렬로 일반화 질량 계산", "상용 Radioss 고유치 해석 결과 직접 사용"],
        ["STL / STEP", "gmsh 격자 → 쉘(두께 지정) 또는 tet10 → 고유치 해석", "형상 + 물성 → 모달"],
    ], [34 * mm, 90 * mm, 46 * mm], caption="구조 모델 소스")
    s += H2("3.2 구조-공력 연성 (스플라인)")
    s.append(P("구조 절점 변위를 공력 격자의 힘점(1/4 코드)과 콜로케이션점(3/4 코드)으로 보간한다. 기본은 "
               "무한판 스플라인(IPS, Harder &amp; Desmarais 1972)으로, 무한 판의 처짐 해를 기저로 쓴다."))
    s.append(eq(r"w(x,y)=a_0+a_1x+a_2y+\sum_{i=1}^{N}F_i\,r_i^2\ln r_i^2,\qquad r_i^2=(x-x_i)^2+(y-y_i)^2"))
    s.append(P("기울기 ∂w/∂x는 해석적으로 미분해 수직세(normalwash)에 쓴다. 이 밖에 TPS/다중이차 RBF(원본 "
               "Thin_plate_spline.f90 대체), 막대 모델용 보 스플라인(굽힘+비틀림 강체 단면 가정), 원본의 "
               "역거리 4-corner 보간(legacy)을 제공한다. 공력 격자 변위는 결과 파일(modes.txt의 aero_h 표)로도 "
               "출력되어 GUI 3-D 모드 표시에 쓰인다."))
    return s


def sec_theory_aero():
    s = []
    s += H2("3.3 이중격자법 (DLM)")
    s.append(P("양력면을 사다리꼴 상자로 나누고, 각 상자의 1/4 코드선에 압력 이중항(doublet) 선을, 3/4 코드 "
               "중앙에 콜로케이션점을 둔다(Albano &amp; Rodden 1969). 콜로케이션점의 수직세와 상자 압력차 ΔCp는 "
               "공력 영향계수 행렬 D로 연결된다."))
    s.append(eq(r"\frac{w_r}{V}=\sum_s D_{rs}(k,M)\,\Delta C_{p,s},\qquad \mathbf{D}=\mathbf{D}_0+\mathbf{D}_1(k)"))
    s.append(P("정상 부분 D<sub>0</sub>는 말굽와류 VLM으로(압축성은 Prandtl-Glauert 신장), 진동 증분 D<sub>1</sub>은 "
               "Landahl(1967) 커널을 Laschka 11항 지수근사와 Rodden-Taylor-McIntosh(1998)의 스팬방향 4차 근사로 적분한다."))
    s.append(eq(r"D_{1,rs}=\frac{\Delta x_s}{8\pi}\int_{-e}^{e}\left[\frac{P_1}{r_1^2}+\frac{P_2}{r_1^4}\right]d\eta,"
                r"\qquad P_j=\left(K_j\,e^{-i\omega x_0/V}-K_{j0}\right)T_j"))
    s.append(P("모드 j의 수직세와 일반화 공기력(GAF)은 다음과 같다. h는 공력 격자에서의 모드 수직변위이다."))
    s.append(eq(r"\frac{w_j}{V}=\frac{\partial h_j}{\partial x}+\frac{ik}{b}h_j,\qquad "
                r"\Delta\mathbf{C}_p=\mathbf{D}^{-1}\mathbf{W},\qquad \mathbf{Q}(k)=\mathbf{H}_f^{T}\,\mathrm{diag}(A)\,\Delta\mathbf{C}_p"))
    s.append(P("수직 돌풍 w<sub>g</sub> e<super>iω(t-(x-x<sub>0</sub>)/V)</super>에 대해서는 침투 지연을 포함한 수직세 열 "
               "-cos γ e<super>-ik(x-x<sub>0</sub>)/b</super>로 돌풍 GAF 열 Q<sub>g</sub>를, 조종면은 힌지선 회전에 따른 "
               "강체 회전 모드를 추가해 Q<sub>δ</sub>를 얻는다. 대칭/반대칭 조건은 y=0 거울상 상자로 처리한다."))
    s.append(P("<b>띠이론 (Theodorsen).</b> 3-D 효과를 무시한 2-D 극한으로, 각 스팬 위치의 단면에 비정상 판 이론을 적용한다. "
               "플런지 h와 피치 α(탄성축 기준)에 대한 단위 스팬당 양력과 모멘트는 Theodorsen 함수 C(k)로 표현되며, "
               "FLAAP은 이를 같은 모달 GAF 표 형식으로 만들어 PK 해석에 넣는다. Sears 함수는 돌풍 열에 쓴다."))
    s.append(eq(r"C(k)=\frac{H_1^{(2)}(k)}{H_1^{(2)}(k)+iH_0^{(2)}(k)}"))
    s.append(eq(r"L'=\pi\rho b^2\left(\ddot h+V\dot\alpha-ba\ddot\alpha\right)+2\pi\rho Vb\,C(k)\left[\dot h+V\alpha+b\left(\frac{1}{2}-a\right)\dot\alpha\right]"))
    s += H2("3.4 천음속 미소교란 (TSD)")
    s.append(P("DLM은 선형 이론이므로 M ≈ 0.9-1.0에서 충격파와 두께 효과가 만드는 '천음속 dip'을 예측하지 못한다. "
               "FLAAP의 TSD 솔버는 CAP-TSD/XTRAN3S 계열의 3-D 주파수영역 시간선형화 정식화를 따른다. 정상 "
               "평균유동은 보존형 비선형 TSD 방정식을 Murman 형 형식 전환 차분(완전 보존 충격점 연산자)과 "
               "뉴턴 반복으로 푼다."))
    s.append(eq(r"\frac{\partial}{\partial x}\left[(1-M^2)\phi_x-\frac{\gamma+1}{2}M^2\phi_x^2\right]+\phi_{yy}+\phi_{zz}=0,"
                r"\qquad C_p=-2\phi_x"))
    s.append(P("비정상 섭동 φ = φ<sub>0</sub> + φ<sub>1</sub> e<super>i k̄ t</super>는 평균유동에 대해 선형화된다. "
               "계수 A가 국소 Mach 수 정보를 담아, 초음속 영역과 충격파 위치가 진동 압력에 반영된다."))
    s.append(eq(r"\frac{\partial}{\partial x}\left[A\,\phi_{1,x}\right]+\phi_{1,yy}+\phi_{1,zz}+M^2\bar k^2\phi_1-2i\bar kM^2\phi_{1,x}=0,"
                r"\quad A=1-M^2-(\gamma+1)M^2\phi_{0,x}"))
    s.append(P("φ<sub>0</sub> = 0이면 정확한 선형 압축성 포텐셜 방정식이 되므로 두께 0, 아음속에서 DLM과 일치해야 하며 "
               "이를 단위시험으로 확인한다. 날개면에는 Neumann 조건 φ<sub>1,z</sub> = w/V, 후류에는 순환의 대류 조건, "
               "원방에는 복소 좌표 신장 PML(Berenger 1994)을 둔다. 결과 ΔCp를 같은 스플라인으로 적분해 "
               "TSD 기반 Q(k)를 만들면 이후의 모든 해석(플러터, 돌풍)이 그대로 동작한다."))
    s.append(P("<b>양력 TSD 평균유동 (받음각, 캠버).</b> 대칭 두께 문제는 z → -z 대칭이라 윗면 해만 풀면 된다. 받음각 α와 캠버 "
               "z<sub>c</sub>(x, y)가 있으면 이 대칭이 깨지므로 윗면 포텐셜 U(x, y, z)와 거울상 아랫면 포텐셜 L(x, y, ζ) = φ(x, y, -ζ)를 "
               "같은 반공간 격자에 두고, 각자 자기 계수 A와 형식 전환 스위치를 가진 TSD 방정식을 푼다. 두 장은 z = 0 평면에서만 결합된다. "
               "날개면에서는 얇은 날개 경계조건(윗면/아랫면 기울기)을, 날개 밖 상류·측면에서는 연속(U = L)과 법선 유속 연속을, "
               "후류에서는 법선 유속 연속과 뒷전 순환 U - L의 대류를 둔다."))
    s.append(eq(r"U_z=\frac{t'}{2}+z_c'-\alpha,\quad L_\zeta=\frac{t'}{2}-z_c'+\alpha\ \ (\mathrm{wing}),\qquad "
                r"(U-L)(x)=(U-L)_{TE}\,e^{-i\bar k(x-x_{TE})}\ \ (\mathrm{wake})"))
    s.append(P("압력은 Cp<sub>u</sub> = -2U<sub>x</sub>, Cp<sub>l</sub> = -2L<sub>x</sub>이고 압력차는 ΔCp = 2(U - L)<sub>x</sub>"
               "(비정상에서는 2[(U<sub>1</sub> - L<sub>1</sub>)<sub>x</sub> + i k̄ (U<sub>1</sub> - L<sub>1</sub>)])이다. "
               "(U - L)/2가 대칭 솔버의 포텐셜 역할을 하므로 GAF, ΔCp, 뒷전 포텐셜 계산을 그대로 재사용한다. "
               "α = 0이고 캠버가 없으면 대칭 솔버와 기계 정밀도(10<super>-15</super>)로 일치해야 하며, 받음각에 따라 윗면 충격파는 뒤로 "
               "이동하고 약한 아랫면 충격파는 사라진다. 이 평균유동에 대해 선형화하면 받음각이 있는 천음속 플러터를 구할 수 있다."))
    return s


def sec_theory_flutter():
    s = []
    s += H2("3.5 플러터 해석: K, PK, p 방법")
    s.append(P("<b>K 방법(V-g).</b> 조화 운동을 가정하고 인공 구조감쇠 g를 도입하여 k마다 복소 고유치를 푼다. "
               "g가 0을 지나 양이 되는 점이 플러터이다."))
    s.append(eq(r"\left[\mathbf{K}(1+ig)-\omega^2\left(\mathbf{M}+\frac{\rho b^2}{2k^2}\mathbf{Q}(k)\right)\right]\eta=0"))
    s.append(P("<b>PK 방법.</b> 속도마다 각 가지(branch)에 대해 k를 반복 갱신하여 공기력의 실수부를 강성으로, "
               "허수부를 감쇠로 넣은 고유치 문제를 맞춘다(British/Nastran PK). 감쇠 g = 2γ가 물리적 의미를 가진다."))
    s.append(eq(r"\left[\mathbf{M}p^2+\left(\mathbf{C}-\frac{\rho bV}{2k}\mathbf{Q}^I\right)p+\left(\mathbf{K}-q\,\mathbf{Q}^R\right)\right]\eta=0,"
                r"\qquad p=\omega(\gamma+i)"))
    s.append(P("<b>p 방법.</b> 아래의 Roger RFA로 만든 상태공간 행렬 A(V)의 고유치를 직접 구한다. 반복이 없고 임의 감쇠에서도 "
               "정확하며, 시간영역 모델과 같은 행렬을 쓴다. 세 방법 모두 모드 추적(고유벡터 상관 + 진동수 연속성)과 "
               "g의 부호 변화 보간으로 플러터(ω &gt; 0)와 발산(ω → 0)을 자동 식별하고, 원본 형식의 "
               "Vg_save.dat / Vf_save.dat도 출력한다."))
    s += H2("3.6 유리함수 근사(RFA)와 상태공간 모델")
    s.append(P("표로 주어진 Q(ik)를 Laplace 영역으로 확장하기 위해 Roger(1977) 근사를 최소자승으로 맞춘다. "
               "지연근 β<sub>l</sub>은 원본 기본값 1.7 k<sub>max</sub>(l/(n<sub>lag</sub>+1))<super>2</super>이다. "
               "k = 0 정확 일치 구속, Karpel 최소상태 근사도 제공한다."))
    s.append(eq(r"\mathbf{Q}(\bar p)\approx\mathbf{A}_0+\mathbf{A}_1\bar p+\mathbf{A}_2\bar p^2+\sum_{l=1}^{n_{lag}}\mathbf{A}_{l+2}\frac{\bar p}{\bar p+\beta_l},"
                r"\qquad \bar p=\frac{sb}{V}"))
    s.append(P("지연 상태 y<sub>l</sub> = p̄/(p̄+β<sub>l</sub>) ξ 를 도입하면 1차 상태방정식이 된다. 입력은 돌풍 "
               "w<sub>g</sub>, ẇ<sub>g</sub>와 조종면 δ, δ̇, δ̈이며, 단면 하중과 상자 압력도 같은 상태의 선형 출력이다."))
    s.append(eq(r"\dot{\mathbf{x}}=\mathbf{A}(V)\mathbf{x}+\mathbf{B}\mathbf{u},\quad \mathbf{x}=[\eta,\dot\eta,y_1..y_L,z_1..z_L],"
                r"\quad \dot y_l=\dot\xi-\frac{V}{b}\beta_l y_l"))
    return s


def sec_theory_loads():
    s = []
    s += H2("3.7 돌풍 하중")
    s.append(P("<b>이산 돌풍 (CS-25.341(a)).</b> 1-cos 형상의 수직 돌풍을 기준점에서 정의하고, 상자별 침투 지연은 Q<sub>g</sub>에 "
               "포함된다. 돌풍 길이 H = 9-107 m를 훑어 각 출력의 최대값 포락선과 임계 H를 구한다."))
    s.append(eq(r"w_g(s)=\frac{U_{ds}}{2}\left(1-\cos\frac{\pi s}{H}\right),\quad 0\leq s\leq 2H,\qquad "
                r"U_{ds}=U_{ref}\,F_g\left(\frac{H}{106.68}\right)^{1/6}"))
    s.append(P("<b>연속 난류 (CS-25.341(b)).</b> von Kármán 스펙트럼(L = 762 m)에 대해 주파수응답 H<sub>y</sub>(ω)로 "
               "응답 PSD를 구하고, 설계 하중 y<sub>1g</sub> ± U<sub>σ</sub> Ā와 영교차 수 N<sub>0</sub>를 계산한다."))
    s.append(eq(r"\Phi_w(\Omega)=\sigma^2\frac{L}{\pi}\frac{1+\frac{8}{3}(1.339L\Omega)^2}{\left[1+(1.339L\Omega)^2\right]^{11/6}},"
                r"\qquad \bar A=\frac{\sigma_y}{\sigma_w},\quad N_0=\frac{1}{2\pi}\frac{\sigma_{\dot y}}{\sigma_y}"))
    s += H2("3.8 정적 공탄성: 트림과 발산")
    s.append(P("k = 0의 GAF로 탄성·강체 모드와 트림 변수(받음각, 피치 각속도, 조종면 등) ξ를 연결한다. 자유 변수 수를 "
               "방정식 수와 같게 지정하여 푸는 Nastran SOL 144형 모달 정식화이다."))
    s.append(eq(r"\mathbf{K}\eta-q\,\mathbf{Q}_0\,\xi=-\mathbf{M}\,\mathbf{a},\qquad a_{plunge}=n_zg"))
    s.append(P("발산 동압은 det(K - q Q<sub>0,ee</sub>) = 0의 최소 양의 근이며, 조종 효율은 탄성/강체 조종미계수 비, "
               "역전 동압은 효율이 0이 되는 q로 구한다."))
    s += H2("3.9 기동 하중 (CS-25.331(c))")
    s.append(P("강체(플런지, 피치) 모드와 탄성 모드를 전체 FE 질량행렬로 결합하고 RFA 상태공간에 조종면을 처방 입력으로 "
               "넣는다. 받음각은 피치 모드 수직세 α = θ - ż/V로 생기며 단주기 운동이 자연히 포함된다. (c)(1) 급격 입력과 "
               "(c)(2) 체크 기동 δ(t) = δ<sub>1</sub> sin ωt (0 ≤ t ≤ 3π/2ω)을 구현하고, δ<sub>1</sub>은 한계 하중배수에 "
               "정확히 도달하도록 선형 스케일한다. 단면 하중은 힘 합산법(FSM: 공기력 + 관성력)과 모드 변위법(MDM)으로 복원한다."))
    s.append(eq(r"\mathbf{F}_{FSM}=q\,\Delta C_p A\,\mathbf{n}-\mathbf{M}_{fe}\mathbf{\Phi}\ddot\eta,\qquad "
                r"\mathbf{F}_{MDM}=\mathbf{K}_{fe}\mathbf{\Phi}\eta"))
    return s


def sec_theory_nonlinear_ase():
    s = []
    s += H2("3.10 조종면 freeplay와 한계주기진동(LCO)")
    s.append(P("조종면 힌지의 유격(freeplay) δ<sub>fp</sub>는 힌지 모멘트-회전 관계를 구간선형으로 만든다. 선형 플러터 "
               "속도 이하에서도 유한 진폭의 LCO가 생길 수 있어 CS-25.629와 AC 25.629-1B가 별도 평가를 요구한다."))
    s.append(eq(r"M_h=K_h\,f(\delta),\qquad f(\delta)=\delta-\delta_{fp}\ (\delta>\delta_{fp}),\quad 0\ "
                r"(|\delta|\leq\delta_{fp}),\quad \delta+\delta_{fp}\ (\delta<-\delta_{fp})"))
    s.append(P("<b>시간적분.</b> 조종면 회전을 자유도로 갖는 RFA 상태공간에 비선형 힌지 복원력을 피드백으로 넣고 "
               "적분하여 LCO 진폭과 진동수를 얻는다. <b>기술함수(describing function).</b> 진폭 A의 조화 운동에 대한 "
               "1차 고조파 등가 강성은 다음과 같으며, 등가 힌지 강성으로 선형 플러터 해석을 반복하면 진폭-속도 LCO "
               "곡선을 빠르게 얻는다. 두 방법의 일치가 검증 지표이다."))
    s.append(eq(r"\frac{K_{eq}(A)}{K_h}=1-\frac{2}{\pi}\left[\arcsin\frac{\delta_{fp}}{A}+\frac{\delta_{fp}}{A}\sqrt{1-\frac{\delta_{fp}^2}{A^2}}\right],\qquad A>\delta_{fp}"))
    s += H2("3.11 공력-서보-탄성 (ASE)")
    s.append(P("조종면 작동기를 2차 시스템으로, 센서(가속도/각속도)를 공탄성 상태의 선형 출력으로 두고 제어기 K(s)를 "
               "결합해 폐루프 상태공간을 만든다. 폐루프 고유치의 속도 추적으로 능동 플러터 억제(AFS) 효과를, "
               "돌풍 입력 응답으로 돌풍 하중 완화(GLA)를, 루프 전달함수의 Nyquist/Bode로 이득·위상 여유를 평가한다."))
    s.append(eq(r"\frac{\delta(s)}{\delta_c(s)}=\frac{\omega_a^2}{s^2+2\zeta_a\omega_as+\omega_a^2},\qquad "
                r"y=\mathbf{C}_s\mathbf{x},\qquad \delta_c=-K(s)\,y"))
    s.append(eq(r"\dot{\mathbf{x}}_{cl}=\mathbf{A}_{cl}(V,K)\,\mathbf{x}_{cl}+\mathbf{B}_{g}w_g,\qquad L(s)=K(s)\,G_{y\delta}(s)"))
    return s


def sec_theory_vlm_supersonic():
    s = []
    s += H2("3.12 와류격자법 (VLM)과 비정상 와류격자법 (UVLM)")
    s.append(P("<b>정상 VLM.</b> 양력면을 패널로 나누고 1/4 코드에 말굽와류, 3/4 코드에 콜로케이션점을 둔다. "
               "비관통 조건으로 순환 Γ를 구하고, 쿠타-주코프스키 정리로 양력을, Trefftz 평면의 후류 수직세로 "
               "유도항력을 계산한다(압축성은 Prandtl-Glauert 신장)."))
    s.append(eq(r"\sum_s A_{rs}\Gamma_s=-\mathbf{V}_\infty\cdot\mathbf{n}_r,\qquad "
                r"L=\rho V_\infty\sum_s\Gamma_s\Delta y_s,\qquad D_i=-\frac{\rho}{2}\int_{-b/2}^{b/2}\Gamma(y)\,w_{T}(y)\,dy"))
    s.append(P("<b>비선형 UVLM.</b> Katz &amp; Plotkin(2001)의 와류고리 정식화를 따른다. 패널마다 와류고리 세기 Γ를 두고, "
               "매 시간 스텝 Δt = Δx/U마다 뒷전 고리가 후류로 방출된다(Kelvin 정리). 날개 운동(피치, 플런지)과 돌풍 "
               "수직세 w<sub>g</sub>(x - Ut)는 콜로케이션 조건에 들어가고, 후류는 고정(평판) 또는 자유(국소 유도속도로 "
               "대류, roll-up)로 둘 수 있다. 하중은 Simpson, Palacios &amp; Murua(2013)처럼 각 속박 와류 선분에 국소 속도"
               "(자유류 + 물체 운동 + 돌풍 + 유도속도)로 쿠타-주코프스키 정리를 적용하고 비정상 항 ρ(∂Γ/∂t)A n을 더해 구한다. "
               "비압축성(M = 0) 정식화이다."))
    s.append(eq(r"\sum_s A_{rs}\Gamma_s^{n}=-\left(\mathbf{V}_\infty+\mathbf{v}_{w}^{n}-\mathbf{v}_{b}^{n}+w_g\mathbf{e}_z\right)\cdot\mathbf{n}_r,"
                r"\qquad \mathbf{F}=\sum_{\mathrm{seg}}\rho\,\Gamma_{\mathrm{seg}}\,\mathbf{V}_{\mathrm{loc}}\times\mathbf{l}_{\mathrm{seg}}+\sum_s\rho\,\frac{\partial\Gamma_s}{\partial t}A_s\mathbf{n}_s"))
    s.append(P("양력은 자유류 수직 성분, 항력은 자유류 방향 성분이다(x 하류, z 위). 항력에는 유도항력과 앞전 흡입력의 "
               "비정상 효과가 포함되므로 조화 플런지에서 평균 추력(C<sub>T</sub> = -C<sub>D</sub>)이 나온다(Garrick 1936). "
               "결과는 CL, CD, CY, Cm, L, D, M<sub>y</sub>의 시간 이력과 스트립 양력 분포로 저장된다."))
    s.append(P("<b>선형 UVLM.</b> 평판 형상과 길이 UΔt의 고정 평판 후류 주위로 선형화하면(Murua, Palacios &amp; Graham 2012) "
               "이산시간 상태공간 모델이 된다. 상태는 날개·후류 순환이고 입력은 모드 변위·속도와 돌풍이다. "
               "z = e<super>iωΔt</super>에서의 전달함수로 GAF Q(k)를 만들어 PK 해석에 쓰거나, Newmark 구조 적분과 결합한 "
               "단일 이산시간 계 E x<super>n+1</super> = F x<super>n</super>의 고유치가 단위원을 벗어나는 속도로 플러터를 직접 찾는다."))
    s.append(eq(r"\boldsymbol{\Gamma}^{n+1}=\mathbf{A}_a\boldsymbol{\Gamma}^{n}+\mathbf{B}_a\mathbf{u}^{n+1},\qquad "
                r"\mathbf{Q}(k)=\mathbf{C}_a\left(z\mathbf{I}-\mathbf{A}_a\right)^{-1}\mathbf{B}_a+\mathbf{D}_a,\quad z=e^{2ik\Delta x/c}"))
    s += H2("3.13 초음속 공력 (M > 1)")
    s.append(P("DLM은 아음속 전용이므로 M &gt; 1에서는 두 가지 방법을 쓴다. <b>피스톤 이론</b>(Lighthill 1953, "
               "Ashley &amp; Zartarian 1956)은 국소 압력이 국소 수직세에만 의존한다는 고 Mach·고 진동수 근사이다. "
               "1차(4/M), Van Dyke 형(Ackeret 4/β 강성과 (M²-2)/(M²-1) 감쇠 계수), 평균 두께 기울기를 넣은 3차 형을 제공한다."))
    s.append(eq(r"\Delta C_p=-\frac{4}{\beta}\left[\frac{\partial h}{\partial x}+\frac{M^2-2}{M^2-1}\,\frac{ik}{b}h\right],"
                r"\qquad \beta=\sqrt{M^2-1}"))
    s.append(P("<b>Mach box 법</b>(Pines, Dugundji &amp; Neuringer 1955; Zartarian &amp; Hsu 1955)은 선형 초음속 포텐셜 "
               "방정식의 기본해(초음속 소스)를 대각선이 Mach 선과 나란한 상자(Δx = βΔy)로 적분한다. 날개 옆 diaphragm 상자는 "
               "φ = 0(압력차 없음) 조건으로 미지 수직세를 풀고(Evvard 1950), 뒷전 뒤 후류 상자는 Δp = 0으로 둔다. "
               "상자의 영향은 하류 행에만 미치므로 행 단위의 명시적 전진으로 풀리며, 아음속 앞전과 측면 끝단 효과가 포함된다."))
    s.append(eq(r"\phi(x,y)=-\frac{1}{\pi}\iint_{\mathrm{Mach\ cone}}\frac{w(\xi,\eta)\,e^{-i\lambda(x-\xi)}\cos(\kappa R)}{R}\,d\xi\,d\eta,"
                r"\quad R=\sqrt{(x-\xi)^2-\beta^2(y-\eta)^2},\ \lambda=\frac{\omega M^2}{U\beta^2},\ \kappa=\frac{\omega M}{U\beta^2}"))
    s.append(P("두 방법 모두 DLM과 같은 AeroModes, 스플라인, GAF 표 형식을 출력하므로 aero_method = piston 또는 machbox로 "
               "바꾸기만 하면 PK, K, p, 돌풍, 상태공간 해석이 그대로 동작한다. 사전검증은 M ≤ 1에서 machbox/piston을 "
               "오류로, M &lt; 1.6에서 피스톤 이론을 경고로, M &gt; 0.3의 UVLM을 경고로 표시한다."))
    return s


def val_intro(R, md):
    s = []
    s += H1("4. 검증")
    s.append(P("모든 결과는 저장소의 설정 파일과 스크립트로 재생성된다(configs/*, validation/report.py, "
               "tools/run_all_configs.py). 표의 수치는 이 백서를 만들 때의 검증 보고서(validation/REPORT.md)와 실행 결과에서 "
               "자동으로 읽어 온다."))
    s.append(P("<b>검증 원칙.</b> 각 절의 제목은 <b>어떤 해석 기능을 검증했는지</b>이며, Goland 날개나 AGARD 445.6 같은 사례는 "
               "그 아래 항목이다. 비교 기준은 우선순위 순으로 <b>1순위 실험, 2순위 독립된 다른 해석·코드, 3순위 해석해</b>를 쓰고, "
               "상위 기준이 없으면 그 사실과 이유를 적는다. 코드 내부 비교(극한 일치, 해법 간 일치, 입력 경로 왕복)는 "
               "<b>자기일관성</b>으로 따로 표시하며 검증 근거로 세지 않는다."))
    rows = val_summary_rows(md)
    if rows:
        s += table(rows, [9 * mm, 42 * mm, 24 * mm, 27 * mm, 24 * mm, 44 * mm], font=7.0,
                   caption="검증 대상별 비교 기준과 대표 결과 (4.1-4.11)")
    return s


def val_dlm_head(R, md):
    s = []
    s += H2("4.2 선형 DLM 플러터 검증 (아음속, DLM + PK)")
    s += val_basis(md, 2)
    return s


def val_goland(R, md):
    s = []
    s += H3("[2순위 타 해석] Goland 날개 3-D 플러터 (ρ = 1.02 kg/m³)")
    s.append(P("Goland(1945) 균일 외팔보 날개(반스팬 6.096 m, 코드 1.829 m, 탄성축 33 % 코드, 무게중심 43 % 코드)는 "
               "굽힘-비틀림 플러터의 대표적 벤치마크이다. 20개 보 요소 모델의 고유진동수는 Goland(1945) 해석해"
               "(Savino et al. 2021 Table 2에 수록)와 ±0.15 % 이내로 일치한다. 플러터는 두 가지로 비교한다. "
               "(1) Goland 원 논문과 같은 띠이론(Theodorsen)을 같은 보 모드에 적용해 ρ = 1.225 kg/m³에서 137.2 m/s / "
               "70.7 rad/s와 비교하고, (2) 3-D 양력면 결과는 문헌들이 사용한 ρ = 1.02 kg/m³에서 각 문헌값"
               "(Murua 2012, Wang 2006, Savino 2021)과 개별 비교한다."))
    gt = md_table(md, "### 2.2 ", 0)
    if gt:
        s += table(gt, [40 * mm, 50 * mm, 35 * mm, 30 * mm],
                   caption="Goland 날개 3-D 공력 플러터, ρ = 1.02 kg/m³ (문헌값과 개별 비교)")
    fm = fig_modes("modal_goland", "goland")
    if fm:
        s += figure(fm, "Goland 날개 진동 모드 (보 모델 → 보 스플라인 → DLM 격자 6×12)", max_h=55 * mm)
    fv = fig_flutter_compare("flutter_goland_all", "goland")
    if fv:
        s += figure(fv, "Goland 날개 V-g / V-f. 실선 PK, 점 K, 파선 p 방법. ★: 플러터점. 세 방법이 0.2 % 이내로 일치하며 "
                        "DLM의 3-D 팁 효과로 띠이론(137 m/s)보다 높다 (해면 밀도).", max_h=110 * mm)
    fc = fig_cp("flutter_goland_pk", "goland", 0.3, (1, 2))
    if fc:
        s += figure(fc, "Goland 날개 진동 압력분포 ΔCp (k ≈ 0.3, 2·3차 모드 단위 일반화 좌표당, 실수/허수부)", max_h=60 * mm)
    return s


def val_agard(R, md):
    s = []
    s += H3("[1순위 실험] AGARD 445.6 플러터 경계 (Yates 1987)")
    s.append(P("AGARD 445.6 약화 모델 3은 45° 후퇴 NACA 65A004 날개의 풍동 플러터 시험(Yates 1987, NASA TM-100492)이다. "
               "기준값은 TM-100492 Table II의 공기 중 시험점이며, FSI와 진동수비는 Yates와 같이 b<sub>s</sub> = 0.2794 m, "
               "ω<sub>α</sub> = 239.3 rad/s로 무차원화한다(이 값으로 표의 V, μ, FSI가 0.1 % 이내로 서로 일치함을 확인). "
               "흔히 '실험값'으로 인용되는 38.17 / 48.35 / 91.54 Hz는 Yates의 FE 모델 값이고, 지상진동시험 측정값은 "
               "9.60 / 38.10 / 50.70 / 98.50 Hz이다. 시험의 질량비를 맞추기 위해 ρ<sub>eq</sub> = ρ m<sub>model</sub>/m<sub>test</sub>를 "
               "쓰고 플러터 속도 지수 FSI = V<sub>f</sub>/(b<sub>s</sub> ω<sub>α</sub> √μ)로 비교한다."))
    fm = fig_modes("modal_agard", "agard")
    if fm:
        s += figure(fm, "AGARD 445.6 모드 (MITC4 직교이방성 판, IPS 스플라인)", max_h=55 * mm)
    at = md_table(md, "### 2.1 ", 0)
    if at:
        s += table(at, None, caption="AGARD 445.6 FSI, 측정 진동수 + FE 모드형상 (DLM + PK)")
    s += figure(VAL / "figures" / "agard_fsi.png", "AGARD 445.6 FSI와 진동수비. M ≤ 0.9에서 1-3 % 일치, M > 0.95 천음속 영역에서 "
                                                    "선형 DLM은 dip을 놓친다.", max_h=80 * mm)
    fc = fig_cp("flutter_agard_M0901", "agard", 0.1, (0, 1))
    if fc:
        s += figure(fc, "AGARD 445.6, M = 0.901, k ≈ 0.1의 모드 1·2 ΔCp", max_h=60 * mm)
    return s


def val_tsd(R, md):
    s = []
    s += H2("4.4 TSD 천음속 플러터 검증 (비선형 TSD 평균류 + 시간선형화 비정상 해 + PK)")
    s += val_basis(md, 4)
    s += H3("[1순위 실험] AGARD 445.6 천음속 플러터 경계 (Yates 1987)")
    tt = md_table(md, "### 4.1 ", 0)
    if tt:
        s += table(tt, None, caption="DLM vs TSD 플러터 (측정 진동수)", font=7.2)
    s += figure(VAL / "figures" / "agard_tsd_fsi.png", "AGARD 445.6 FSI: 실험, DLM, TSD. TSD는 M = 0.954에서 실험과 +0.2 %, "
                                                        "M 0.901-0.960에서 -3.5~+1.1 %로 천음속 dip을 재현한다.", max_h=75 * mm)
    t2 = md_table(md, "### 4.2 ", 0)
    if t2:
        s += H3("[자기일관성] 두께 0 TSD vs DLM")
        s += table(t2, None, caption="두께 0 TSD와 DLM의 FSI 비교 (검증 근거 아님)", font=7.6)
    ft = fig_tsd_steady()
    if ft:
        s += figure(ft, "TSD 정상해, M = 0.954, NACA 65A004 두께 유동. 검은 선은 음속선(Cp*)으로 날개 상면에 초음속 영역이 형성된다.",
                    max_h=70 * mm)
    return s


def val_legacy(R, md):
    s = []
    s += H3("[2순위 타 해석] 원본 Fortran FLAAP.exe (Conyers 판 날개, 같은 입력 덱)")
    s.append(P("동일 입력덱(Conyers 판 날개, 4-corner legacy 스플라인 + DLM + PK)에서 Python 판과 Fortran FLAAP.exe의 "
               "pk_results_mode*.dat를 비교했다. 저속 구간 모드 8-10은 감쇠 수 %, 진동수 1-3 %로 일치한다. 120 m/s 이상에서 "
               "모드 7-8이 결합하는 구간은 두 프로그램의 가지 추적 규칙 차이로 레이블이 달라진다."))
    lt = md_table(md, "### 2.3 ", 0)
    if lt:
        s += table(lt, None, caption="원본 Fortran 대비 PK 결과", font=7.4)
    s += figure(VAL / "figures" / "legacy_vs_fortran.png", "원본 Fortran FLAAP(점) vs Python(선) V-g / V-f", max_h=85 * mm)
    return s


def val_gust(R, md):
    s = []
    s += H2("4.10 돌풍 응답 검증 (1-cos 이산 돌풍, 연속 난류) — 자기일관성 검사")
    s += val_basis(md, 10)
    s.append(P("같은 1-cos 돌풍(H = 20 m, U = 10 m/s, V = 100 m/s)에 대해 RFA 상태공간 시간적분과, 표로 된 GAF를 "
               "직접 쓰는 주파수영역 해(FFT)를 비교하면 최대/최소 응답이 1 % 이내로 일치한다. 이는 RFA 근사와 시간적분이 "
               "정확함을 보인다."))
    gt = md_table(md, "### 10.1 ", 0)
    if gt:
        s += table(gt, None, caption="돌풍 응답: 시간영역 vs 주파수영역", font=7.6)
    s += figure(VAL / "figures" / "goland_gust.png", "Goland 날개 1-cos 돌풍 응답: 시간영역(RFA)과 주파수영역(GAF + FFT)", max_h=80 * mm)
    fe = fig_curves("gust_goland_discrete", "curves_gust_envelope.txt", ["root.Fz", "root.Mx", "root.My", "eta:mode 1"],
                    "gust_env", "CS-25.341(a) 1-cos 돌풍 포락선 (정규화)", normalise=True)
    if fe:
        s += figure(fe, "돌풍 길이 H에 대한 최대 응답 포락선 (V = 100 m/s, 해면)", max_h=65 * mm)
    fp = fig_curves("gust_goland_continuous", "curves_turbulence_psd.txt", ["eta:mode 1", "eta:mode 2", "etadd:mode 1"],
                    "psd", "von Kármán 난류 응답 PSD")
    if fp:
        s += figure(fp, "연속 난류 응답 PSD (L = 762 m)", max_h=65 * mm)
    return s


def val_trim(R, md):
    s = []
    s += H2("4.11 정적 공탄성 트림·기동하중 — 자기일관성 검사")
    s += val_basis(md, 11)
    s.append(P("날개 + 꼬리날개 데모 항공기(자유-자유 탄성 모드 + 강체 플런지/피치)로 1-g 트림과 CS-25.331(c)(2) 체크 "
               "기동을 계산했다. 트림 양력이 반중량과 정확히 일치하고, 체크 기동은 지정 한계 하중배수 2.5에 도달하도록 "
               "조종 입력이 스케일된다."))
    tr = R.get("trim/ex1/config_demo_trim.txt", {}).get("summary", {})
    mn = R.get("maneuver/ex1/config_demo_checked.txt", {}).get("summary", {})
    dv = R.get("trim/ex2/config_goland_divergence.txt", {}).get("summary", {})
    rows = [["항목", "값"]]
    if tr:
        rows += [["1-g 트림 받음각", f"{fmt(tr.get('alpha_deg'), 3)} deg"],
                 ["1-g 트림 승강타", f"{fmt(tr.get('elevator_deg'), 3)} deg"],
                 ["날개뿌리 굽힘모멘트 Mx (1 g)", f"{fmt((tr.get('loads') or {}).get('root.Mx'), 0)} N·m"]]
    if mn:
        rows += [["단주기 ω<sub>n</sub> / ζ", f"{fmt(mn.get('short_period_wn'), 2)} rad/s / {fmt(mn.get('short_period_zeta'), 3)}"],
                 ["체크 기동 최대 승강타", f"{fmt(mn.get('peak_control_deg'), 2)} deg"],
                 ["체크 기동 증분 Mx (FSM)", f"{fmt((mn.get('peak_loads') or {}).get('root.Mx.fsm'), 0)} N·m"]]
    if dv:
        rows += [["Goland 날개 발산 속도 (해면)", f"{fmt(dv.get('V_divergence'))} m/s"]]
    if len(rows) > 1:
        s += table(rows, [80 * mm, 60 * mm], caption="정적 공탄성 및 기동 결과")
    fth = fig_timehist("maneuver_demo_checked", "timehist_maneuver.txt",
                       ["elevator [deg]", "n_z", "root.Mx.fsm", "root.My.fsm"], "maneuver",
                       "CS-25.331(c)(2) checked manoeuvre, V = 100 m/s")
    if fth:
        s += figure(fth, "체크 피치 기동 시간응답: 승강타, 하중배수, 날개뿌리 굽힘/비틀림 하중", max_h=95 * mm)
    return s


def val_freeplay(R, md):
    s = []
    s += H2("4.8 조종면 프리플레이 비선형 LCO 검증")
    s += val_basis(md, 8)
    s.append(P("Conner et al.(1997)의 3자유도(h, α, β) 전형 단면(Vasconcellos 2016 매개변수, δ<sub>fp</sub> = 1°)에서 "
               "시간적분(TM)과 기술함수(DF) 해를 비교했다. 선형 플러터 속도는 9.232 m/s(문헌 9.2 m/s)이며, freeplay는 "
               "K<sub>eq</sub> → 0인 3.05 m/s부터 LCO를 만들고 DF 안정 가지는 1.88 m/s의 아임계 fold까지 이어진다. "
               "TM 피크 진폭은 DF보다 3-4 % 크고(고조파 성분), 진동수는 1.5 % 이내로 일치한다. 진폭은 δ<sub>fp</sub>에 "
               "정확히 비례한다(δ<sub>fp</sub> 2° / 1° 비 2.0000). 문헌과 정량 비교한 것은 선형 플러터 속도뿐이며, LCO 진폭·진동수는 "
               "TM-DF 상호검증이다. Vasconcellos(2016)의 실험 LCO(δ<sub>fp</sub> ±3.57°, 측정 감쇠행렬, 5차 피치 강성 다항식, "
               "실험 U<sub>f</sub> 15.2 m/s)는 실험 장치 고유의 비선형 강성을 포함하므로 향후 과제로 남긴다."))
    fs = R.get("freeplay/ex1/config_section_df.txt", {}).get("summary", {})
    rows = [["V [m/s]", "TM A [deg]", "TM A<sub>1</sub> [deg]", "DF A [deg]", "ΔA", "TM f [Hz]", "DF f [Hz]", "Δf"]]
    for r in fs.get("time_marching", []):
        if r.get("A_df_deg"):
            rows.append([fmt(r["V"], 1), fmt(r["A_tm_deg"], 3), fmt(r.get("A1_tm_deg"), 3), fmt(r["A_df_deg"], 3),
                         f"{100 * (r['A_tm_deg'] / r['A_df_deg'] - 1):+.1f} %", fmt(r["f_tm_hz"], 3),
                         fmt(r["f_df_hz"], 3), f"{100 * (r['f_tm_hz'] / r['f_df_hz'] - 1):+.1f} %"])
    if len(rows) > 1:
        s += table(rows, None, font=7.4,
                   caption=f"전형 단면 freeplay LCO, config freeplay/ex1: TM vs DF (선형 V<sub>f</sub> "
                           f"{fmt(fs.get('V_f_linear'), 3)}, K = 0 V<sub>f</sub> {fmt(fs.get('V_f0'), 3)}, "
                           f"fold {fmt(fs.get('V_fold'), 2)} m/s). A: 피크, A<sub>1</sub>: 1차 고조파")
    s += figure(VAL / "figures" / "freeplay_keq.png",
                "기술함수 등가 강성 K<sub>eq</sub>(A)/K<sub>h</sub>와 이에 따른 선형 플러터 속도", max_h=60 * mm)
    s += figure(VAL / "figures" / "freeplay_bifurcation.png",
                "분기 선도: DF 가지(선)와 TM 상승/하강 연속 스윕(점). 2-3 m/s에서 상승 스윕은 감쇠, 하강 스윕은 "
                "LCO를 유지하는 이력(hysteresis)을 재현한다.", max_h=75 * mm)
    s += figure(VAL / "figures" / "freeplay_time_history.png",
                "LCO 시간응답: 조종면 회전이 유격 경계(±δ<sub>fp</sub>)를 넘나들며 일정 진폭으로 수렴", max_h=70 * mm)
    s += figure(VAL / "figures" / "freeplay_phase_plane.png", "LCO 위상평면", max_h=60 * mm)
    fg = R.get("freeplay/ex3/config_goland_flap_df.txt", {}).get("summary", {})
    s.append(P("<b>Goland 날개 + 플랩.</b> 50-83 % 스팬, 75 % 코드 힌지 플랩(K<sub>h</sub> = 9340 N·m/rad)을 붙인 Goland "
               f"날개(DLM 모달 공력)에서 δ<sub>fp</sub> = 0.5°일 때 선형 V<sub>f</sub> = {fmt(fg.get('V_f_linear'))} m/s, "
               f"K = 0 V<sub>f</sub> = {fmt(fg.get('V_f0'))} m/s이다. 145 m/s에서 TM/DF 진폭비 1.066/1.062, "
               "진동수 76.48/76.49 rad/s로 0.4 % 이내 일치한다. 150 m/s 이상에서는 고조파 기여로 TM 진폭이 더 크다."))
    s += figure(VAL / "figures" / "goland_freeplay_lco.png",
                "Goland + 플랩 freeplay LCO 진폭/진동수: DF(선) vs TM(점)", max_h=70 * mm)
    return s


def val_ase(R, md):
    s = []
    s += H2("4.9 공력서보탄성 (ASE) 플러터 억제·돌풍하중 완화 — 자기일관성 검사")
    s += val_basis(md, 9)
    s.append(P("같은 전형 단면에 2차 작동기(15 Hz, ζ 0.7)와 h, α, dh/dt, dα/dt 센서를 두고 정적 출력 피드백을 걸었다. "
               "K = 0 폐루프는 개루프 p 방법 V<sub>f</sub> = 9.245 m/s를 정확히 재현하여 결합이 올바름을 보인다. "
               "피치+피치율 법칙(A)은 12.97 m/s(+40 %), 전 출력 피드백(B)은 17.40 m/s(+88 %)로 플러터 속도를 높인다."))
    fa = R.get("ase/ex1/config_section_afs.txt", {}).get("summary", {})
    fmg = R.get("ase/ex2/config_section_margins.txt", {}).get("summary", {})
    rows = [["항목", "값"],
            ["개루프 V<sub>f</sub>", f"{fmt(fa.get('V_f_open'), 2)} m/s"],
            ["폐루프 V<sub>f</sub> (법칙 B)", f"{fmt(fa.get('V_f_closed'), 2)} m/s "
                                         f"(+{fmt(100 * (fa.get('increase') or 0), 1)} %)"]]
    mg = fmg.get("margins") or {}
    if mg:
        rows += [[f"V = {fmt(mg.get('V'), 1)} m/s 개루프 불안정 극점 P / Nyquist N<sub>ccw</sub>",
                  f"{mg.get('P')} / {mg.get('N_ccw')}"],
                 ["하한 이득여유 GM<sub>low</sub>", f"{fmt(mg.get('gm_lower_db'), 2)} dB"],
                 ["위상여유 PM", f"{fmt(mg.get('pm_deg'), 1)} deg @ {fmt(mg.get('w_pm'), 1)} rad/s"]]
    s += table(rows, [95 * mm, 60 * mm], caption="능동 플러터 억제와 안정여유 (전형 단면, 법칙 B, config ase/ex1-2)")
    s.append(P("개루프가 불안정한 12-17 m/s에서도 Nyquist 판정(N<sub>ccw</sub> = P)과 폐루프 고유치가 일치한다. "
               "17 m/s에서 위상여유가 0.9°로 줄고 18 m/s에서 판정이 불안정으로 바뀌어 V<sub>f</sub> = 17.40 m/s와 정합한다."))
    s += figure(VAL / "figures" / "ase_vg_rootlocus.png", "개루프/폐루프 V-g와 속도 근궤적", max_h=75 * mm)
    s += figure(VAL / "figures" / "ase_bode_nyquist.png", "루프 전달함수 L(iω)의 Bode/Nyquist 선도 (법칙 B)", max_h=75 * mm)
    gl = R.get("ase/ex3/config_goland_gla.txt", {}).get("summary", {})
    s.append(P("<b>돌풍 하중 완화(GLA).</b> Goland + 플랩 날개에서 날개뿌리 굽힘 변형률 게이지 신호를 플랩에 되먹임"
               "(k = -2×10<super>-6</super> rad/N·m, V = 100 m/s, U<sub>ds</sub> = 5 m/s)하면 1-cos 돌풍의 뿌리 "
               f"굽힘모멘트가 12-25 % 줄고, 플러터 속도도 {fmt(gl.get('V_f_open'))} → {fmt(gl.get('V_f_closed'))} m/s로 "
               "높아진다."))
    rows = [["H [m]", "Mx 개루프 [N·m]", "Mx GLA [N·m]", "Mx 감소", "Fz 감소", "플랩 최대 [deg]"]]
    for r in (gl.get("gla") or {}).get("rows", []):
        fz = 1 - r["Fz_closed"] / r["Fz_open"] if r.get("Fz_open") else None
        rows.append([fmt(r.get("H"), 1), fmt(r.get("Mx_open"), 0), fmt(r.get("Mx_closed"), 0),
                     f"{fmt(100 * r['Mx_reduction'], 1)} %", f"{fmt(100 * fz, 1) if fz is not None else '-'} %",
                     fmt(r.get("flap_peak_deg"), 2)])
    if len(rows) > 1:
        s += table(rows, None, caption="Goland + 플랩 GLA 결과 (config ase/ex3)", font=7.6)
    s += figure(VAL / "figures" / "ase_gla_goland.png", "GLA 개루프/폐루프 돌풍 응답", max_h=75 * mm)
    return s


def val_struct(R, md):
    s = []
    s += H2("4.1 구조 고유진동 해석 검증 (판·쉘·보·고체 FE와 외부 모델 입력)")
    s += val_basis(md, 1)
    s += H3("[1순위 실험] AGARD 445.6 판 날개 고유진동수 (MITC4 판 16×16)")
    t11 = md_table(md, "### 1.1 ", 0)
    if t11:
        s += table(t11, None, caption="AGARD 445.6 고유진동수 [Hz]: 측정(Yates 1987), FLAAP, Yates FE 모델", font=7.8)
    s += H3("[2순위 타 해석] 직사각형 외팔 판 (CFFF, a/h = 100): Leissa (1973)")
    t12 = md_table(md, "### 1.2 ", 0)
    if t12:
        s += table(t12, None, caption="무차원 고유치 λ = ωa²√(ρh/D): 쉘 요소 vs Leissa", font=7.2)
    s += H3("[3순위 해석해] Goland 보 정확해와 Timoshenko 보 이론 (고체 요소)")
    t13 = md_table(md, "### 1.3 ", 0)
    if t13:
        s += table(t13, None, caption="Goland 날개 보 FE 고유진동수 vs 정확해(Goland 1945)", font=7.8)
    t13b = md_table(md, "### 1.3 ", 1)
    if t13b:
        s += table([r[:5] for r in t13b], None, caption="외팔 강철 보: 고체 요소 vs Timoshenko 보", font=7.4)
    s += H3("[자기일관성] 입력 경로 왕복과 외부 FE 결과")
    s.append(P("새 구조 요소와 입력 경로를 해석해와 교차 검증으로 확인했다. 외팔 정사각 판의 무차원 진동수 "
               "λ = ωa²√(ρh/D)는 Leissa(1969) 값과, 외팔 직사각 보는 Timoshenko 보 이론과 비교한다 "
               "(validation/structures_validation.py)."))
    s += table([
        ["검증 항목", "기준", "결과"],
        ["쉘 quad4 / tri3 16×16, 외팔 판 λ<sub>1-5</sub>", "Leissa (1969)", "오차 0.7 % 이내"],
        ["고체 hex8 20×2×2 / tet10 16×2×2, 외팔 보 f<sub>1-3</sub>", "Timoshenko 보 (수렴 1-D FE)", "+0.7~1.4 % / +0.7~0.9 %"],
        ["tet4 40×4×4 (참고)", "Timoshenko 보", "+19~28 % (전단 잠김, tet10 권장)"],
        ["자유-자유 강체 모드 6개", "0 Hz", "최대 0.003 Hz"],
        ["OpenRadioss 입력덱 (/SHELL, /SH3N, /BRICK, /TETRA10)", "Leissa / Timoshenko", "판 0.6 %, 보 1.5 % 이내"],
        ["ANSYS 모달 (MAPDL listing, PRNSOL, mm-t 단위, Workbench)", "기준 FE 모달", "Δf 0.001 %, MAC = 1.000"],
        ["Radioss /EIG VTK (합성 샘플 A001-A006)", "입력덱 FE 모달", "진동수 일치(Δf &lt; 0.001 %)"],
        ["CalculiX .frd (long SI + .dat / short mm-t-s, 정적 스텝·강체 모드 포함)", "기준 FE 모달",
         "Δf ≤ 0.0003 %, MAC = 1.000"],
        ["CalculiX .frd → DLM + PK 플러터 (flutter/ex11)", "같은 판의 FLAAP FE", "V<sub>f</sub> 70.585 m/s 일치"],
        ["gmsh .geo 판 (transfinite 6×12, Recombine → quad4)", "판 FE 기준 모드", "Δf 3·10<super>-5</super> %, .msh 왕복 일치"],
        ["gmsh .geo OCC 막대 (tet10, Physical Root 구속)", "Euler-Bernoulli 1차 굽힘", "+0.8 % (b/h = 5 판 효과)"],
        ["gmsh .geo 판 → DLM + PK (flutter/ex12)", "CalculiX .frd 모델 (ex11)", "V<sub>f</sub> 70.585 m/s 일치"],
        ["STL 상자 → tet10, h = 0.02 m", "Timoshenko 보", "+0.26~0.34 %"],
        ["STEP 상자 (mm 단위) → tet10", "Timoshenko 보", "+0.27~0.35 %"],
        ["STL 판 → gmsh tri3 쉘", "Leissa (1969)", "오차 0.74 % 이내"],
    ], [78 * mm, 36 * mm, 56 * mm], caption="구조 요소와 모델 입력 검증")
    s += figure(VAL / "figures" / "structures_shell_convergence.png", "외팔 판 쉘 요소 수렴 (Leissa 대비 오차)", max_h=70 * mm)
    s += figure(VAL / "figures" / "structures_plate_modes.png", "외팔 판 모드형상 (quad4 16×16)", max_h=65 * mm)
    s += figure(VAL / "figures" / "structures_beam_modes.png", "외팔 보 고체 요소 모드형상", max_h=65 * mm)
    ma = R.get("modal/ex4/config_ansys_modal.txt", {}).get("summary", {})
    mr = R.get("modal/ex5/config_openradioss_modal.txt", {}).get("summary", {})
    ms = R.get("modal/ex6/config_stl_modal.txt", {}).get("summary", {})
    me = R.get("modal/ex7/config_openradioss_eig.txt", {}).get("summary", {})
    mc = R.get("modal/ex8/config_calculix_modal.txt", {}).get("summary", {})
    mg = R.get("modal/ex9/config_gmsh_plate.txt", {}).get("summary", {})
    mgs = R.get("modal/ex10/config_gmsh_solid.txt", {}).get("summary", {})
    fl = R.get("flutter/ex7/config_stl_wing_pk.txt", {}).get("summary", {})

    def flist(sm, n=6):
        return ", ".join(fmt(f, 2) for f in (sm.get("frequencies_hz") or [])[:n]) or "-"
    rows = [["예제 (configs/)", "입력", "진동수 [Hz]"],
            ["modal/ex4 ANSYS", "plate_modal.out + plate.cdb", flist(ma)],
            ["modal/ex5 OpenRadioss", "plate_cantilever_0000.rad", flist(mr)],
            ["modal/ex6 STL", "wing.stl + E 0.4 GPa, ν 0.3, ρ 120 kg/m³, tet10", flist(ms)],
            ["modal/ex7 Radioss /EIG", "plate_cantileverA001-A006.vtk + 입력덱(질량)", flist(me)],
            ["modal/ex8 CalculiX", "plate.frd + plate.dat (ccx *FREQUENCY)", flist(mc)],
            ["modal/ex9 gmsh 판", "plate_wing.geo (quad4, clamp_group Root)", flist(mg)],
            ["modal/ex10 gmsh 고체", "box_beam.geo (OCC, tet10)", flist(mgs)]]
    s += table(rows, [38 * mm, 62 * mm, 70 * mm], caption="외부 모델 입력 예제의 모달 결과", font=7.4)
    pk = fl.get("flutter_pk") or {}
    s.append(P("<b>형상 → 모달 → 플러터.</b> NACA 단면 STL 날개(반스팬 1.2 m, 코드 0.40-0.21 m)에 물성만 주면 gmsh로 "
               "tet10 격자(7672 절점)를 만들어 고유치 해석을 하고, 평면 위 절점 201개로 IPS 스플라인하여 DLM + PK 해석을 "
               f"이어서 수행한다. M = 0.3, 해면에서 플러터 속도 {fmt(pk.get('velocity'))} m/s, "
               f"{fmt(pk.get('freq_hz'), 2)} Hz(굽힘-비틀림 결합)이며 전체 계산은 약 11 s이다(config flutter/ex7)."))
    s += figure(VAL / "figures" / "stl_wing_modes.png", "STL 날개 tet10 모드형상", max_h=70 * mm)
    s += figure(VAL / "figures" / "stl_wing_vg.png", "STL 날개 V-g / V-f (DLM + PK, M = 0.3)", max_h=80 * mm)
    return s


def val_uvlm(R, md):
    s = []
    s += H2("4.7 VLM / UVLM 공력 및 UVLM 플러터 검증")
    s += val_basis(md, 7)
    s.append(P("정상 VLM, 비선형 UVLM, 선형 UVLM을 고전 해석해와 문헌값으로 검증했다(validation/uvlm_validation.py, "
               "tests/test_uvlm.py). 2차원 결과는 가로세로비 1000의 띠 날개로 얻었다."))
    s += table([
        ["검증 항목", "기준", "결과"],
        ["타원 날개 A 4 / 8, 12×48: 스팬 효율 e", "양력선 e = 1", "1.0001 / 1.0009"],
        ["타원 날개 A 8: CL<sub>α</sub>", "Helmbold 4.906", "4.812-4.832 (-1.5~-1.9 %)"],
        ["A 1000 띠 날개: CL<sub>α</sub>", "얇은 익형 2π", "6.2748 (-0.13 %)"],
        ["직사각 A 4: CL<sub>α</sub>, e, x<sub>ac</sub>/c", "양력면 이론 3.6-3.7, 0.98-0.99", "3.628, 0.991, 0.232"],
        ["Wagner 지시 양력 (nc 20)", "Wagner 함수", "s ≥ 2에서 0.7 %, s ≥ 10에서 0.02 % 이내"],
        ["Küssner 돌풍 양력 (nc 12)", "Küssner 함수", "s ≥ 4에서 1.8 %, s = 20에서 0.06 %"],
        ["조화 피치 |CL| (k 0.5 / 1), nc → ∞", "Theodorsen", "비 1.0006 / 1.0016 (Richardson)"],
        ["조화 플런지 평균 추력 (k 0.1-1)", "Garrick (1936)", "-0.05~-1.6 %"],
        ["선형 UVLM GAF (nc 32, k ≤ 1)", "Theodorsen", "최대 1.8 % (nc에 1차 수렴)"],
        ["3-D 직사각 A 4 급출발 (α 5°)", "정상 VLM CL", "99.6 % 도달, 자유후류 차이 &lt; 0.1 %"],
    ], [70 * mm, 46 * mm, 54 * mm], caption="VLM / UVLM 검증 요약")
    s.append(P("<b>Goland 날개 플러터(ρ = 1.02 kg/m³).</b> 같은 4개 보 모드로 선형 UVLM + PK를 격자 세분화하면 "
               "8×12 / 16×12 / 32×12에서 172.7 / 170.0 / 168.5 m/s이며 Richardson 외삽값 167.0 m/s는 Murua et al.(2012)의 "
               "UVLM 결과 165.0 m/s와 1.2 % 차이이다. 같은 8×12 격자의 이산시간 상태공간 고유치 해석은 173.2 m/s, "
               "10.67 Hz로 PK(172.7 m/s, 10.63 Hz)와 일치하여 주파수영역과 시간영역 정식화의 일관성을 보인다. "
               "같은 모드의 DLM 8×12 + PK는 165.1 m/s이다."))
    s += figure(VAL / "figures" / "uvlm_validation.png",
                "VLM / UVLM 검증: 타원 날개 스팬 하중, Wagner/Küssner, 피치 수렴, Garrick 추력, GAF 오차, Goland 감쇠",
                max_h=120 * mm)
    rows = [["예제 (configs/uvlm)", "내용", "주요 결과"]]
    for cfg, name, f in [
        ("uvlm/ex1/config_rect_wing_steady.txt", "정상 VLM, A 8 직사각",
         lambda m: f"CLα {fmt(m.get('CL_alpha_per_rad'), 3)}/rad, e {fmt(m.get('span_efficiency'), 4)}"),
        ("uvlm/ex2/config_rect_wing_impulsive.txt", "급출발, A 4, α 5°",
         lambda m: f"CL {fmt(m.get('CL_final'), 4)} = 정상 VLM의 {fmt(100 * (m.get('CL_final_over_steady') or 0), 1)} %"),
        ("uvlm/ex3/config_plunge_thrust.txt", "2-D 플런지 k 0.5",
         lambda m: f"C_T {m.get('CT_mean_last_cycle', 0):.3e} (Garrick 7.466e-04)"),
        ("uvlm/ex4/config_pitch_theodorsen.txt", "2-D 피치 k 0.5, 1°",
         lambda m: f"|CL| {fmt(m.get('CL_1st_harmonic_amplitude'), 4)} (Theodorsen 0.0800)"),
        ("uvlm/ex5/config_gust_1cos.txt", "1-cos 돌풍, V 20 m/s",
         lambda m: f"ΔCL peak {fmt(m.get('dCL_gust_peak'), 4)} @ {fmt(m.get('t_peak'), 3)} s"),
        ("uvlm/ex6/config_goland_uvlm_flutter.txt", "Goland 선형 UVLM",
         lambda m: f"PK {fmt((m.get('flutter_pk') or {}).get('velocity'), 1)} m/s, 상태공간 "
                   f"{fmt((m.get('flutter_uvlm_state_space') or {}).get('velocity'), 1)} m/s"),
    ]:
        sm = R.get(cfg, {}).get("summary", {})
        rows.append([cfg.split("/")[1], name, f(sm) if sm else "-"])
    s += table(rows, [24 * mm, 46 * mm, 100 * mm], caption="UVLM 예제의 시간영역 결과", font=7.4)
    p = fig_timehist("uvlm_gust", "timehist_uvlm.txt", ["CL", "CD", "Cm"], "uvlm_gust_hist")
    if p:
        s += figure(p, "1-cos 돌풍 UVLM 양력/항력/모멘트 계수 시간 이력 (configs/uvlm/ex5)", max_h=70 * mm)
    return s


def val_ss(R, md):
    s = []
    s += H2("4.6 초음속 공력·플러터 검증 (Mach box, piston theory)")
    s += val_basis(md, 6)
    s.append(P("피스톤 이론과 Mach box를 정상 해석해, 2-D 진동 평판 엄밀해, 패널 플러터 문헌값, AGARD 445.6 초음속 시험과 "
               "비교했다(validation/supersonic_validation.py)."))
    s += table([
        ["검증 항목", "기준", "결과"],
        ["피스톤 이론 (Van Dyke), M 2 정상 양력", "Ackeret 4/β", "0.00 %"],
        ["Mach box 직사각 A 40 / A 2, nx 80", "Ackeret / Evvard 끝단 보정", "+0.01 % / +0.15 %"],
        ["Mach box 삼각 날개 (아음속 앞전, m 0.69)", "원뿔 유동 해", "+0.17 %"],
        ["2-D 진동 평판 GAF, M 1.2 / 2, k ≤ 1", "Garrick-Rubinow", "최대 0.024 % (nx 80)"],
        ["패널 플러터 λ<sub>cr</sub>: 2-D 단순지지 / 고정", "Hedgepeth 343.3 / Dowell 636.6", "343.30 / 636.51"],
        ["패널 플러터 λ<sub>cr</sub>: 정사각 단순지지", "Dowell 512.6", "512.36"],
        ["AGARD 445.6 M 1.072 / 1.141 FSI (측정 진동수)", "Yates 0.3201 / 0.4031", "0.3474 / 0.4593 (+8.5 / +13.9 %)"],
    ], [70 * mm, 50 * mm, 50 * mm], caption="초음속 공력 검증 요약")
    s.append(P("AGARD 445.6의 초음속 점에서 피스톤 이론은 적용 범위 밖(M 1.07-1.14)이어서 플러터를 찾지 못하거나 크게 "
               "벗어난다. 선형 포텐셜 이론인 Mach box는 시험보다 8-14 % 높으며, 이는 낮은 초음속에서 두께와 점성 효과가 "
               "선형 이론에 없기 때문이다(DLM의 M 0.954 오차 +10.6 %와 같은 성격). 예제 flutter/ex9의 M 2 평판은 "
               "피스톤 이론 588.3 m/s, Mach box 538.3 m/s로 고 Mach 근사 차이가 약 9 %이다."))
    s += figure(VAL / "figures" / "supersonic_validation.png", "초음속 공력 검증", max_h=110 * mm)
    return s


def val_strip(R, md):
    s = []
    s += H2("4.3 스트립 이론 (Theodorsen) 플러터 검증")
    s += val_basis(md, 3)
    s += H3("[2순위 타 해석] Goland 날개, ρ = 1.225 kg/m³")
    t = md_table(md, "### 3.1 ", 0)
    if t:
        s += table(t, [70 * mm, 40 * mm, 40 * mm], caption="Goland 날개 플러터, 띠이론 + PK: Goland (1945), Patil et al. (2001)")
    s.append(P("같은 모델에서 3-D 날개끝 효과를 포함한 DLM은 153 m/s로 띠이론(137 m/s)보다 높다(4.2절 자기일관성 비교). "
               "띠이론 검증은 Goland 1945가 쓴 공력 모델과 같은 가정에서 구조 모달 모델과 PK 해석 체인이 올바름을 확인하는 데 의미가 있다."))
    return s


def val_lifting(R, md):
    s = []
    s += H2("4.5 양력 TSD 평균류 검증 (받음각·캠버, 충격파 위치)")
    s += val_basis(md, 5)
    s += H3("[2순위 타 해석] NACA 0012, M = 0.8: AGARD-AR-211 Euler 합의값 대비")
    t = md_table(md, "### 5.2 ", 0)
    if t:
        s += table(t, None, caption="NACA 0012 직사각 날개(A = 16) 뿌리 단면의 받음각 스윕 vs Euler", font=7.4)
    s += H3("[3순위 해석해] 선형 극한 양력 기울기 (M 0.5, 평판, A = 8)")
    t = md_table(md, "### 5.3 ", 0)
    if t:
        s += table(t, [60 * mm, 50 * mm], caption="CL<sub>α</sub>: 양력 TSD vs Helmbold", font=7.8)
    s += H3("[자기일관성] 정확 일치 검사와 받음각 평균유동 위 플러터 경향")
    t = md_table(md, "### 5.4 ", 0)
    if t:
        s += table(t, [70 * mm, 90 * mm], caption="양력 TSD 일치 검사 (검증 근거 아님)", font=7.6)
    t = md_table(md, "### 5.4 ", 1)
    if t:
        s += table(t, None, caption="AGARD 445.6, M = 0.954: 평균 받음각에 따른 플러터 (α = 0만 실험 있음)", font=7.8)
    s += figure(VAL / "figures" / "tsd_lifting.png",
                "양력 TSD: 받음각 스윕의 Cp, 충격파 위치, CL과 받음각 평균유동 위 플러터", max_h=105 * mm)
    return s


def val_gaps(R, md):
    s = []
    s += H2("4.12 실험 비교가 없는 항목과 향후 과제")
    s.append(P("검증 근거의 강도는 기능마다 다르다. 아래 항목은 현재 1순위(실험) 비교가 없거나 자기일관성에만 의존하므로, "
               "필요한 실험 또는 독립 코드 자료를 함께 적는다. 이 표가 현재 검증의 한계이며 다음 검증의 우선순위이다."))
    t = md_table(md, "## 12. ", 0)
    if t:
        s += table(t, [46 * mm, 44 * mm, 80 * mm], caption="실험 비교 공백과 필요한 자료", font=7.6)
    s += [Spacer(1, 6)]
    s.append(P("<b>종합.</b> 구조 FE(AGARD 445.6 판 진동수), 선형 DLM 플러터(AGARD 445.6 M ≤ 0.9, +1.2~+3.3 %), TSD 천음속 플러터"
               "(M 0.954에서 +0.2 %, 천음속 dip 재현), 초음속 Mach box(M 1.07 / 1.14에서 +8.5 / +13.9 %)는 시험 측정값(진동 시험, 풍동 플러터 시험)과 직접 비교했다. "
               "UVLM, 스트립 이론, 패널 플러터, VLM은 문헌의 독립된 해석 또는 해석해와 일치한다. 프리플레이 LCO 진폭, ASE, 돌풍, "
               "트림·기동은 자기일관성 수준이며, 양력 TSD의 충격파 위치는 Euler 해와만 비교했다."))
    return s


def sec_examples(R):
    s = []
    s += H2("2.8 예제 설정과 실행 결과")
    s.append(P("저장소의 configs/ 폴더에는 모든 경로의 예제가 들어 있으며, tools/run_all_configs.py가 이를 일괄 실행한다. "
               "아래 표는 이 백서 작성 시의 실행 결과이다."))
    rows = [["설정 파일", "상태", "시간 [s]", "주요 결과"]]
    for cfg, r in sorted(R.items()):
        sm = r.get("summary") or {}
        key = ""
        for k in ("flutter_pk", "flutter_p", "flutter_k"):
            if sm.get(k):
                key = f"{k[8:]}: {sm[k]['kind']} V = {fmt(sm[k]['velocity'])} m/s, f = {fmt(sm[k]['freq_hz'], 2)} Hz"
                break
        if not key:
            if "V_divergence" in sm:
                key = f"V_div = {fmt(sm['V_divergence'])} m/s"
            elif "alpha_deg" in sm:
                key = f"α = {fmt(sm['alpha_deg'], 3)}°, δ_e = {fmt(sm.get('elevator_deg'), 3)}°"
            elif "peak_nz" in sm:
                key = f"n_z,max = {fmt(sm['peak_nz'], 2)}, δ_max = {fmt(sm.get('peak_control_deg'), 2)}°"
            elif "cp_min" in sm:
                key = f"Cp_min = {fmt(sm['cp_min'], 3)}, 초음속 격자 {sm.get('n_supersonic')}"
            elif "critical_H" in sm:
                key = f"root.Mx 임계 H = {fmt(sm['critical_H'].get('root.Mx'), 0)} m"
            elif "A_bar" in sm:
                key = f"A-bar(root.Mx) = {fmt(sm['A_bar'].get('root.Mx'), 0)}"
            elif "V_f_closed" in sm:
                key = f"V_f 개루프 {fmt(sm['V_f_open'], 2)} → 폐루프 {fmt(sm['V_f_closed'], 2)} m/s"
                red = [100 * r["Mx_reduction"] for r in (sm.get("gla") or {}).get("rows", [])]
                if red:
                    key += f", GLA Mx -{fmt(min(red))}~{fmt(max(red))} %"
            elif "margins" in sm:
                mg = sm["margins"]
                key = f"V = {fmt(mg['V'])}: GM_low {fmt(mg.get('gm_lower_db'), 2)} dB, PM {fmt(mg.get('pm_deg'))}°"
            elif "V_f_linear" in sm:
                key = f"V_f {fmt(sm['V_f_linear'], 2)}, LCO 시작 {fmt(sm.get('V_f0'), 2)} m/s"
            elif "lco" in sm:
                key = f"V = {fmt(sm['lco']['V'])}: {sm['lco']['kind']}, A = {fmt(sm['lco'].get('amplitude_deg'), 3)}°"
            elif "CL_alpha_per_rad" in sm:
                key = f"CLα = {fmt(sm['CL_alpha_per_rad'], 3)}/rad, e = {fmt(sm.get('span_efficiency'), 4)}"
            elif "dCL_gust_peak" in sm:
                key = f"ΔCL_gust,max = {fmt(sm['dCL_gust_peak'], 4)} @ t = {fmt(sm.get('t_peak'), 3)} s"
            elif "CT_mean_last_cycle" in sm:
                key = (f"|CL|₁ = {fmt(sm.get('CL_1st_harmonic_amplitude'), 4)}, "
                       f"C_T = {sm['CT_mean_last_cycle']:.2e}")
            elif "CL_final_over_steady" in sm:
                key = f"CL = {fmt(sm['CL_final'], 4)} ({fmt(100 * sm['CL_final_over_steady'], 1)} % 정상 VLM)"
            elif "frequencies_hz" in sm:
                key = "f = " + ", ".join(fmt(f, 2) for f in sm["frequencies_hz"][:4]) + " Hz"
        s_ok = "정상" if r["status"] == "ok" else r["status"]
        rows.append([cfg, s_ok, fmt(r.get("elapsed_s")), key])
    s += table(rows, [62 * mm, 14 * mm, 14 * mm, 80 * mm], caption="예제 설정 실행 결과", font=7.2)
    return s


def sec_limits():
    s = []
    s += H2("2.9 적용 범위와 한계")
    s += bullets([
        "DLM은 아음속 선형 이론이다. M &gt; 0.9에서는 TSD 경로를, M &gt; 1에서는 Mach box(또는 M ≳ 1.6에서 피스톤 이론)를 써야 한다. "
        "초음속 방법은 선형 이론이라 낮은 초음속에서 두께·점성 효과를 포함하지 않는다.",
        "UVLM/VLM은 비압축성, 얇은 양력면, 부착 유동 가정이다. 박리, 실속, 두께 효과는 포함하지 않는다.",
        "TSD 솔버는 단일 양력면과 고전 TSD 계수를 쓴다(Ames 항 없음). 받음각·캠버가 있는 평균유동(양력 TSD)에서 윗면 충격파는 "
        "Euler 해보다 약 0.035 c 앞에 놓이고 약한 아랫면 충격파는 α &gt; 0.5°에서 사라진다. 점성(충격파-경계층 간섭)은 "
        "포함하지 않으며 충격파 위치의 실험 비교는 아직 없다.",
        "Roger RFA의 정확도는 k 범위와 지연항 수에 좌우된다. 사전검증은 kmax와 nlag를 경고 대상으로 다룬다.",
        "외부 모달 결과의 단위는 사용자가 length_scale, mass_scale로 맞춰야 한다(예: in → m 0.0254).",
        "STL 형상의 쉘 모델은 중립면 추출과 균일 두께를 가정한다. 복잡한 조립체는 tet10 고체 또는 상용 FE를 권장한다.",
        "freeplay LCO는 단일 조종면 힌지의 구간선형 강성을 다룬다. 마찰, 다중 비선형성은 시간적분 경로로만 가능하다.",
    ])
    return s


def sec_refs():
    s = []
    s += H2("3.14 참고문헌")
    refs = [
        "Albano, E., Rodden, W.P., A doublet-lattice method for calculating lift distributions on oscillating surfaces in subsonic flows, AIAA J. 7(2), 1969.",
        "Rodden, W.P., Taylor, P.F., McIntosh, S.C., Further refinement of the subsonic doublet-lattice method, J. Aircraft 35(5), 1998.",
        "Harder, R.L., Desmarais, R.N., Interpolation using surface splines, J. Aircraft 9(2), 1972.",
        "Roger, K.L., Airplane math modeling methods for active control design, AGARD CP-228, 1977.",
        "Karpel, M., Design for active flutter suppression and gust alleviation using state-space aeroelastic modeling, J. Aircraft 19(3), 1982.",
        "Hassig, H.J., An approximate true damping solution of the flutter equation by determinant iteration, J. Aircraft 8(11), 1971.",
        "Goland, M., The flutter of a uniform cantilever wing, J. Appl. Mech. 12(4), 1945.",
        "Yates, E.C., AGARD standard aeroelastic configurations for dynamic response I - Wing 445.6, AGARD R-765, 1985; "
        "NASA TM-100492, 1987.",
        "Patil, M.J., Hodges, D.H., Cesnik, C.E.S., Limit-cycle oscillations in high-aspect-ratio wings, J. Fluids Struct. 15, 2001.",
        "Wang, Z., Chen, P.C., Liu, D.D., Mook, D.T., Patil, M.J., Time domain nonlinear aeroelastic analysis for HALE wings, "
        "AIAA 2006-1640.",
        "Murua, J., Palacios, R., Graham, J.M.R., Applications of the unsteady vortex-lattice method in aircraft "
        "aeroelasticity and flight dynamics, Prog. Aerosp. Sci. 55, 2012.",
        "Savino, A. et al., Goland wing flutter with the mid-fidelity solver DUST (Table 2: exact natural frequencies, "
        "3-D flutter results), Energies 14, 2021.",
        "Bisplinghoff, R.L., Ashley, H., Halfman, R.L., Aeroelasticity, Addison-Wesley, 1955.",
        "Sears, W.R., Some aspects of non-stationary airfoil theory and its practical application, J. Aeronaut. Sci. 8(3), 1941.",
        "Batina, J.T., Efficient algorithm for solution of the unsteady transonic small-disturbance equation, J. Aircraft 25(7), 1988.",
        "Murman, E.M., Cole, J.D., Calculation of plane steady transonic flows, AIAA J. 9(1), 1971.",
        "Bathe, K.J., Dvorkin, E.N., A four-node plate bending element based on Mindlin/Reissner plate theory and a mixed interpolation, IJNME 21, 1985.",
        "Tang, D., Dowell, E.H., Flutter and limit-cycle oscillations of a wing with control surface freeplay, J. Aircraft / AIAA papers, 1990s-2010s.",
        "Conner, M.D., Tang, D.M., Dowell, E.H., Virgin, L.N., Nonlinear behavior of a typical airfoil section with control surface freeplay, J. Fluids Struct. 11, 1997.",
        "Katz, J., Plotkin, A., Low-Speed Aerodynamics, 2nd ed., Cambridge University Press, 2001.",
        "Simpson, R.J.S., Palacios, R., Murua, J., Induced-drag calculations in the unsteady vortex lattice method, AIAA J. 51(7), 2013.",
        "Wagner, H., Über die Entstehung des dynamischen Auftriebes von Tragflügeln, ZAMM 5(1), 1925.",
        "Küssner, H.G., Zusammenfassender Bericht über den instationären Auftrieb von Flügeln, Luftfahrtforschung 13, 1936.",
        "Garrick, I.E., Propulsion of a flapping and oscillating airfoil, NACA Report 567, 1936.",
        "Ashley, H., Zartarian, G., Piston theory - a new aerodynamic tool for the aeroelastician, J. Aeronaut. Sci. 23(12), 1956.",
        "Pines, S., Dugundji, J., Neuringer, J., Aerodynamic flutter derivatives for a flexible wing with supersonic and subsonic edges, J. Aeronaut. Sci. 22(10), 1955.",
        "Evvard, J.C., Use of source distributions for evaluating theoretical aerodynamics of thin finite wings at supersonic speeds, NACA Report 951, 1950.",
        "Garrick, I.E., Rubinow, S.I., Flutter and oscillating air-force calculations for an airfoil in a two-dimensional supersonic flow, NACA Report 846, 1946.",
        "Dowell, E.H., Aeroelasticity of Plates and Shells, Noordhoff, 1975.",
        "EASA CS-25 Amendment 27: 25.331, 25.341, 25.629; FAA AC 25.629-1B.",
    ]
    for i, r in enumerate(refs, 1):
        s.append(Paragraph(f"[{i}] {r}", S["small"]))
    return s


def cover_figure():
    m = load("modal_agard/modes.txt")
    a = load("modal_agard/aero_mesh.txt")
    if m is None or a is None:
        return
    corners = np.array(a["corners"])
    H = np.array(m["aero_h"])
    span = np.ptp(corners.reshape(-1, 3), axis=0).max()
    fig = plt.figure(figsize=(10, 4.2), facecolor="#122820")
    for i in range(min(3, H.shape[0])):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d", facecolor="#122820")
        _deformed_mesh_axes(ax, corners, H[i], 0.2 * span, "", cmap="YlGn")
        ax.text2D(0.5, 0.03, f"{m['labels'][i]}  {m['freq_hz'][i]:.1f} Hz", transform=ax.transAxes, ha="center",
                  color="#c9d6cf", fontsize=9)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1, wspace=0)
    fig.savefig(FIG / "cover_modes.png", facecolor=fig.get_facecolor())
    plt.close(fig)


def build():
    R = runs_summary()
    md = (VAL / "REPORT.md").read_text(encoding="utf-8") if (VAL / "REPORT.md").exists() else ""
    cover_figure()
    toc = TableOfContents()
    toc.levelStyles = [S["toc1"], S["toc2"]]
    story = [NextPageTemplate("normal"), PageBreak(), Paragraph("목차", S["h1"]), toc, PageBreak()]
    story += sec_summary(R, md)
    story += [PageBreak()]
    story += sec_intro()
    story += sec_func_overview()
    story += sec_arch(R)
    story += sec_examples(R)
    story += sec_limits()
    story += [PageBreak()]
    story += sec_theory_structure()
    story += sec_theory_aero()
    story += sec_theory_flutter()
    story += sec_theory_loads()
    story += sec_theory_nonlinear_ase()
    story += sec_theory_vlm_supersonic()
    story += sec_refs()
    story += [PageBreak()]
    for fn in (val_intro, val_struct, val_dlm_head, val_agard, val_goland, val_legacy, val_strip, val_tsd,
               val_lifting, val_ss, val_uvlm, val_freeplay, val_ase, val_gust, val_trim, val_gaps):
        story += fn(R, md)
    doc = Doc(PDF)
    doc.multiBuild(story)
    print("wrote", PDF, f"({PDF.stat().st_size / 1e6:.1f} MB, {_fig_no[0]} figures, {_eq_no[0]} equations)")


if __name__ == "__main__":
    build()
