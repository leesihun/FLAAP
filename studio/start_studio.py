"""FLAAP Studio - local browser GUI for the FLAAP launcher.

    python studio/start_studio.py [--port 8765] [--no-browser]

Python standard library only (http.server + threading + subprocess + json).
Serves ``studio/web`` and a small JSON API; binds 127.0.0.1 only.

API
    GET  /api/models                 route specs (modes, key groups, defaults, choices)
    GET  /api/presets                preset planforms (goland, agard445, demo_aircraft, ...)
    GET  /api/configs                configs/**/*.txt and studio_runs/*/config_*.txt
    GET  /api/config?path=           raw config text
    POST /api/validate {text, source?}           preflight diagnostics
    POST /api/run      {text, name, source?}     save + launch FLAAP_main.py, returns job id
    GET  /api/job?id=&since=         job status + log lines from index ``since``
    GET  /api/jobs                   all jobs of this session
    POST /api/cancel   {id}          terminate a running job
    GET  /api/outputs                output/* run dirs with results.txt
    GET  /api/file?path=             raw PNG / TXT / CSV / JSON below output/ only
    GET  /api/result?path=           a FLAAP result text file (results.txt, modes.txt,
                                     flutter_*.txt ...) parsed by flaap.io.results_txt

FLAAP writes every result as a plain-text file (``key value`` lines and
``table ... end`` blocks); the viewers read those files through /api/result.
"""
from __future__ import annotations

import argparse
import json
import mimetypes
import os
import re
import subprocess
import sys
import threading
import time
import traceback
import uuid
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

STUDIO = Path(__file__).resolve().parent
ROOT = STUDIO.parent
WEB = STUDIO / "web"
OUTPUT = ROOT / "output"
CONFIGS = ROOT / "configs"
RUNS = ROOT / "studio_runs"
MAIN = ROOT / "FLAAP_main.py"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from flaap.launcher import specs as S  # noqa: E402
from flaap.io.results_txt import read_result  # noqa: E402
from flaap.launcher.parser import PATH_KEYS, parse_text  # noqa: E402

MIME = {".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8",
        ".mjs": "text/javascript; charset=utf-8", ".css": "text/css; charset=utf-8",
        ".json": "application/json; charset=utf-8", ".png": "image/png", ".svg": "image/svg+xml",
        ".txt": "text/plain; charset=utf-8", ".csv": "text/plain; charset=utf-8", ".ico": "image/x-icon",
        ".woff2": "font/woff2"}
FILE_EXT = {".json", ".png", ".txt", ".csv"}
MAX_BODY = 4 * 1024 * 1024

CHOICES = {
    "model": list(S.ROUTES),
    "structure": list(S.STRUCTURES),
    "symmetry": ["none", "symmetric", "antisymmetric"],
    "spline": ["auto", "ips", "tps", "rbf", "beam", "legacy"],
    "aero_method": ["dlm", "tsd"],
    "material": ["aluminum", "steel", "titanium", "cfrp_quasi", "custom"],
    "fe_type": ["shell", "solid"],
    "sensor": ["tip_accel", "root_bending", "modal"],
}


def _group_name(g) -> str:
    named = [("Common", S.COMMON_KEYS), ("Structure", S.STRUCTURE_KEYS), ("Aerodynamics", S.AERO_KEYS),
             ("Flight", S.FLIGHT_KEYS), ("Gust", S.GUST_KEYS), ("Trim / manoeuvre", S.TRIM_KEYS),
             ("TSD", S.TSD_KEYS), ("Nonlinear / freeplay", S.NONLINEAR_KEYS), ("ASE", S.ASE_KEYS)]
    for n, d in named:
        if g is d:
            return n
    return "Other"


def _jsonable(v):
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    return str(v)


def models_payload() -> dict:
    routes = {}
    for name, spec in S.ROUTES.items():
        groups = []
        for g in (S.COMMON_KEYS,) + tuple(spec.key_groups):
            keys = [{"key": k, "desc": d, "type": t, "default": _jsonable(spec.defaults.get(k)),
                     "choices": CHOICES.get(k)} for k, (d, t) in g.items()]
            groups.append({"name": _group_name(g), "keys": keys})
        routes[name] = {"description": spec.description, "modes": dict(spec.modes), "groups": groups,
                        "defaults": _jsonable(spec.defaults),
                        "required_by_mode": {m: list(r) for m, r in spec.required_by_mode.items()}}
    return {"routes": routes, "order": list(S.ROUTES),
            "indexed": {k: fmt for k, (_, fmt) in S.INDEXED_PATTERNS.items()},
            "path_keys": sorted(PATH_KEYS), "choices": CHOICES}


_PRESETS = None


def presets_payload() -> dict:
    """Preset planforms of the built-in structures (needs numpy; computed once)."""
    global _PRESETS
    if _PRESETS is not None:
        return _PRESETS

    def surf(s):
        return {"name": s.name, "rle": [float(x) for x in s.root_le], "rc": float(s.root_chord),
                "tle": [float(x) for x in s.tip_le], "tc": float(s.tip_chord), "nchord": int(s.nchord),
                "nspan": int(s.nspan),
                "controls": [{"name": c.name, "hinge": float(c.hinge), "eta0": float(c.eta[0]),
                              "eta1": float(c.eta[1])} for c in s.control_surfaces]}
    out = {}
    try:
        from flaap.models.goland import goland_surface
        out["goland"] = [surf(goland_surface(6, 12))]
        g = surf(goland_surface(6, 12))
        g["controls"] = [{"name": "aileron", "hinge": 0.8, "eta0": 0.6, "eta1": 1.0}]
        out["goland_ctrl"] = [g]
    except Exception:  # noqa: BLE001
        traceback.print_exc()
    try:
        from flaap.models.agard445 import agard445_surface
        out["agard445"] = [surf(agard445_surface(8, 10))]
    except Exception:  # noqa: BLE001
        traceback.print_exc()
    try:
        from flaap.models.demo_aircraft import demo_surfaces
        out["demo_aircraft"] = [surf(s) for s in demo_surfaces()]
    except Exception:  # noqa: BLE001
        traceback.print_exc()
    _PRESETS = {"presets": out}
    return _PRESETS


# ----------------------------------------------------------------- config helpers
def _rel(p: Path) -> str:
    return p.resolve().relative_to(ROOT).as_posix()


def _inside(p: Path, base: Path) -> bool:
    try:
        p = p.resolve()
        base = base.resolve()
        return os.path.commonpath([str(p), str(base)]) == str(base)
    except (ValueError, OSError):
        return False


def list_configs() -> list[dict]:
    items = []
    for base, group, pat in ((CONFIGS, "configs", "*.txt"), (RUNS, "studio_runs", "config_*.txt")):
        if not base.is_dir():
            continue
        for p in sorted(base.rglob(pat)):
            if p.is_file():
                items.append({"path": _rel(p), "group": group, "name": p.name,
                              "folder": p.parent.relative_to(base).as_posix(),
                              "mtime": p.stat().st_mtime})
    return items


def safe_config_path(rel: str | None) -> Path | None:
    if not rel:
        return None
    p = (ROOT / rel)
    if p.suffix.lower() != ".txt":
        return None
    if _inside(p, CONFIGS) or _inside(p, RUNS):
        return p.resolve()
    return None


def absolutize_paths(text: str, base_dir: Path) -> str:
    """Rewrite relative path-key values (except output_dir) to absolute paths so a
    config loaded from configs/** still finds its files when saved elsewhere."""
    out = []
    for line in text.splitlines():
        s = line.strip()
        if s and not s.startswith("%"):
            body, _, comment = line.partition("#")
            parts = body.split(None, 1)
            if len(parts) == 2 and parts[0].lower() in PATH_KEYS and parts[0].lower() != "output_dir":
                val = parts[1].strip()
                if val and "," not in val and not Path(val).is_absolute():
                    newp = (base_dir / val).resolve()
                    line = f"{parts[0]:<16s}{newp.as_posix()}" + (f"   #{comment}" if comment else "")
        out.append(line)
    return "\n".join(out) + "\n"


def validate_text(text: str, source: str | None) -> dict:
    src = safe_config_path(source)
    fake = (src.parent if src else RUNS) / "config_studio_check.txt"
    pc = parse_text(text, fake)
    try:
        vals, rep = S.resolve(pc)
    except Exception as exc:  # resolver crash on malformed input -> report it as a diagnostic
        line = None
        for i, ln in enumerate(text.splitlines(), 1):
            toks = ln.split("%")[0].split("#")[0].split()
            if len(toks) > 1 and (toks[1] in ("=", ":") or toks[0].endswith((":", "="))):
                line = i
                break
        diag = {"code": "STUDIO-PARSE", "severity": "error",
                "message": f"could not resolve config ({type(exc).__name__}: {exc})",
                "hint": "FLAAP configs use whitespace-separated 'key value' lines (no '=' or ':')."}
        if line:
            diag["location"] = {"path": fake.name, "line": line, "column": 1}
        return {"summary": {"errors": 1, "warnings": 0, "notices": 0}, "diagnostics": [diag],
                "resolved": {}, "defaulted": [], "ok": False}
    d = rep.to_dict()
    for diag in d["diagnostics"]:
        if "location" in diag:
            diag["location"]["path"] = Path(diag["location"]["path"]).name
    defaults = sorted(k for k in vals if k not in pc.values)
    d["resolved"] = {k: _jsonable(vals[k]) for k in sorted(vals) if k != "output_dir"}
    d["defaulted"] = defaults
    d["ok"] = not rep.has_errors()
    return d


# ----------------------------------------------------------------- jobs
class Job:
    def __init__(self, name: str, cfg: Path, outdir: Path):
        self.id = uuid.uuid4().hex[:12]
        self.name = name
        self.cfg = cfg
        self.outdir = outdir
        self.lines: list[str] = []
        self.status = "starting"
        self.returncode = None
        self.started = time.time()
        self.finished = None
        self.proc: subprocess.Popen | None = None
        self.lock = threading.Lock()

    def start(self):
        cmd = [sys.executable, "-u", str(MAIN), "--config", str(self.cfg), "--output-dir", str(self.outdir),
               "--no-color"]
        env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8", NO_COLOR="1",
                   MPLBACKEND="Agg")
        self._add(f"$ python FLAAP_main.py --config {_rel(self.cfg)} --output-dir {_rel(self.outdir)}")
        flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        self.proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     stdin=subprocess.DEVNULL, env=env, creationflags=flags)
        self.status = "running"
        threading.Thread(target=self._pump, daemon=True).start()

    def _add(self, line: str):
        with self.lock:
            self.lines.append(line)
            if len(self.lines) > 20000:
                del self.lines[:5000]

    def _pump(self):
        assert self.proc and self.proc.stdout
        for raw in iter(self.proc.stdout.readline, b""):
            self._add(raw.decode("utf-8", "replace").rstrip("\r\n"))
        rc = self.proc.wait()
        self.returncode = rc
        self.finished = time.time()
        if self.status != "cancelled":
            self.status = {0: "done", 1: "preflight_failed", 3: "run_failed"}.get(rc, "failed")
        self._add(f"[studio] process exited with code {rc} ({self.status})")

    def to_dict(self, since: int = 0) -> dict:
        with self.lock:
            lines = self.lines[since:]
            n = len(self.lines)
        has_res = (self.outdir / "results.txt").is_file() and self.finished is not None
        return {"id": self.id, "name": self.name, "status": self.status, "returncode": self.returncode,
                "started": self.started, "finished": self.finished,
                "elapsed": round((self.finished or time.time()) - self.started, 1),
                "config": _rel(self.cfg), "output_dir": _rel(self.outdir) if _inside(self.outdir, ROOT) else str(self.outdir),
                "has_results": has_res, "lines": lines, "next": n}


JOBS: dict[str, Job] = {}


def start_job(text: str, name: str, source: str | None) -> Job:
    name = re.sub(r"[^A-Za-z0-9_-]+", "_", (name or "").strip()).strip("_")[:60] or time.strftime("run_%H%M%S")
    src = safe_config_path(source)
    if src is not None:
        text = absolutize_paths(text, src.parent)
    rd = RUNS / name
    rd.mkdir(parents=True, exist_ok=True)
    cfg = rd / f"config_{name}.txt"
    cfg.write_text(text if text.endswith("\n") else text + "\n", encoding="utf-8")
    job = Job(name, cfg, OUTPUT / f"studio_{name}")
    JOBS[job.id] = job
    job.start()
    return job


# ----------------------------------------------------------------- outputs
def list_outputs() -> list[dict]:
    runs = []
    if not OUTPUT.is_dir():
        return runs
    for d in sorted(OUTPUT.iterdir()):
        rj = d / "results.txt"
        if not (d.is_dir() and rj.is_file()):
            continue
        try:
            r = read_result(rj, "results")
        except Exception as exc:  # noqa: BLE001
            r = {"status": "unreadable", "error": str(exc)}
        files = sorted(p.name for p in d.iterdir() if p.is_file() and p.suffix.lower() in FILE_EXT)
        figs = sorted("figures/" + p.name for p in (d / "figures").glob("*.png")) if (d / "figures").is_dir() else []
        runs.append({"name": d.name, "path": _rel(d), "title": r.get("title"), "model": r.get("model"),
                     "mode": r.get("mode"), "status": r.get("status"), "error": r.get("error"),
                     "created": r.get("created"), "elapsed_s": r.get("elapsed_s"),
                     "artifacts": r.get("artifacts", []), "files": files, "figures": figs,
                     "mtime": rj.stat().st_mtime})
    return runs


def safe_output_file(rel: str | None) -> Path | None:
    if not rel or "\x00" in rel:
        return None
    p = (ROOT / rel)
    if p.suffix.lower() not in FILE_EXT:
        return None
    if not _inside(p, OUTPUT):
        return None
    p = p.resolve()
    return p if p.is_file() else None


# ----------------------------------------------------------------- HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = "FLAAPStudio/1.0"

    def log_message(self, fmt, *args):  # quieter console
        if getattr(self.server, "verbose", False):
            sys.stderr.write("[studio] " + (fmt % args) + "\n")

    # ---- responses
    def _send(self, code: int, body: bytes, ctype: str, extra: dict | None = None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)

    def _json(self, data, code: int = 200):
        self._send(code, json.dumps(data, allow_nan=False, default=str).encode("utf-8"), MIME[".json"])

    def _err(self, code: int, msg: str):
        self._json({"error": msg}, code)

    def _body(self) -> dict:
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_BODY:
            raise ValueError("request body too large")
        raw = self.rfile.read(n) if n else b"{}"
        d = json.loads(raw.decode("utf-8") or "{}")
        if not isinstance(d, dict):
            raise ValueError("JSON object expected")
        return d

    def _host_ok(self) -> bool:
        host = (self.headers.get("Host") or "").split(":")[0]
        return host in ("127.0.0.1", "localhost", "")

    # ---- routing
    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        if not self._host_ok():
            return self._err(403, "forbidden host")
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path.startswith("/api/"):
                return self._api_get(u.path, q)
            return self._static(u.path)
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            return self._err(500, f"{type(exc).__name__}: {exc}")

    def do_POST(self):
        if not self._host_ok():
            return self._err(403, "forbidden host")
        u = urlparse(self.path)
        try:
            body = self._body()
        except Exception as exc:  # noqa: BLE001
            return self._err(400, f"bad request: {exc}")
        try:
            if u.path == "/api/validate":
                return self._json(validate_text(str(body.get("text", "")), body.get("source")))
            if u.path == "/api/run":
                text = str(body.get("text", ""))
                if not text.strip():
                    return self._err(400, "empty config")
                job = start_job(text, str(body.get("name", "")), body.get("source"))
                return self._json({"id": job.id, "name": job.name, "config": _rel(job.cfg),
                                   "output_dir": _rel(job.outdir)})
            if u.path == "/api/cancel":
                job = JOBS.get(str(body.get("id")))
                if not job:
                    return self._err(404, "no such job")
                if job.proc and job.proc.poll() is None:
                    job.status = "cancelled"
                    job.proc.terminate()
                return self._json({"id": job.id, "status": job.status})
            return self._err(404, "unknown endpoint")
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            return self._err(500, f"{type(exc).__name__}: {exc}")

    def _api_get(self, path: str, q: dict):
        if path == "/api/models":
            return self._json(models_payload())
        if path == "/api/presets":
            return self._json(presets_payload())
        if path == "/api/configs":
            return self._json({"configs": list_configs()})
        if path == "/api/config":
            p = safe_config_path(q.get("path"))
            if p is None or not p.is_file():
                return self._err(404, "config not found (only configs/** and studio_runs/** .txt)")
            return self._json({"path": _rel(p), "text": p.read_text(encoding="utf-8", errors="replace")})
        if path == "/api/job":
            job = JOBS.get(q.get("id", ""))
            if not job:
                return self._err(404, "no such job")
            try:
                since = max(0, int(q.get("since", 0)))
            except ValueError:
                since = 0
            return self._json(job.to_dict(since))
        if path == "/api/jobs":
            return self._json({"jobs": [dict(j.to_dict(len(j.lines)), lines=[]) for j in JOBS.values()]})
        if path == "/api/outputs":
            return self._json({"outputs": list_outputs()})
        if path == "/api/file":
            p = safe_output_file(q.get("path"))
            if p is None:
                return self._err(404, "file not found or not allowed (output/** .txt/.png/.csv/.json only)")
            return self._send(200, p.read_bytes(), MIME.get(p.suffix.lower(), "application/octet-stream"))
        if path == "/api/result":
            p = safe_output_file(q.get("path"))
            if p is None or p.suffix.lower() != ".txt":
                return self._err(404, "result file not found or not allowed (output/**/*.txt only)")
            return self._json(read_result(p))
        if path == "/api/health":
            return self._json({"ok": True, "root": str(ROOT), "python": sys.version.split()[0]})
        return self._err(404, "unknown endpoint")

    def _static(self, path: str):
        rel = path.lstrip("/") or "index.html"
        p = WEB / rel
        if not _inside(p, WEB) or not p.is_file():
            if "." not in Path(rel).name:
                p = WEB / "index.html"
            else:
                return self._err(404, "not found")
        ctype = MIME.get(p.suffix.lower()) or mimetypes.guess_type(p.name)[0] or "application/octet-stream"
        self._send(200, p.read_bytes(), ctype)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="FLAAP Studio local GUI server")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--no-browser", action="store_true")
    ap.add_argument("--verbose", action="store_true", help="log every HTTP request")
    a = ap.parse_args(argv)
    srv = ThreadingHTTPServer(("127.0.0.1", a.port), Handler)
    srv.daemon_threads = True
    srv.verbose = a.verbose
    url = f"http://127.0.0.1:{a.port}/"
    print(f"FLAAP Studio on {url}  (root {ROOT})  - Ctrl+C to stop", flush=True)
    if not a.no_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for j in JOBS.values():
            if j.proc and j.proc.poll() is None:
                j.proc.terminate()
        srv.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
