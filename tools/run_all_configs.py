"""Run every shipped config (configs/**/config_*.txt) and tabulate the outcome.

    python tools/run_all_configs.py [--filter flutter] [--skip-slow]

Writes output/run_all_summary.txt (a flaap-txt table; docs/whitepaper/build_whitepaper.py reads it
and then each run's results.txt).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from flaap.launcher import run_config  # noqa: E402
from flaap.io.results_txt import TxtDoc, read_result  # noqa: E402

SLOW = ("tsd/ex2", "uvlm/ex2", "uvlm/ex5", "uvlm/ex6", "flutter/ex8")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--filter", default="", help="substring of the config path")
    ap.add_argument("--skip-slow", action="store_true")
    a = ap.parse_args(argv)
    rows = []
    files = sorted((ROOT / "configs").rglob("config_*.txt"))
    for f in files:
        rel = f.relative_to(ROOT / "configs").as_posix()
        if a.filter and a.filter not in rel:
            continue
        if a.skip_slow and any(rel.startswith(s) for s in SLOW):
            continue
        t0 = time.time()
        try:
            res = run_config(f, verbose=False)
            status, err = res["status"], res.get("error")
            summ = res.get("summary", {})
            outdir = res.get("output_dir")
        except Exception as e:  # preflight error
            status, err, summ, outdir = "preflight-error", str(e)[:500], {}, None
        dt = time.time() - t0
        rows.append({"config": rel, "status": status, "elapsed_s": round(dt, 2), "error": err,
                     "summary": summ, "output_dir": outdir})
        print(f"{status:>16s} {dt:7.1f}s  {rel}" + (f"\n{'':>26s}{str(err).splitlines()[-1][:150]}" if err else ""),
              flush=True)
    out = ROOT / "output" / "run_all_summary.txt"
    out.parent.mkdir(exist_ok=True)
    merged = {}
    if out.exists() and a.filter:
        old = read_result(out).get("runs") or []
        merged = {r["config"]: r for r in old}
    merged.update({r["config"]: r for r in rows})
    runs = list(merged.values())
    doc = TxtDoc("run_all_summary", "every shipped config run by tools/run_all_configs.py")
    doc.kv("created", time.strftime("%Y-%m-%d %H:%M"))
    doc.kv("configs", len(runs))
    doc.kv("failed", sum(r["status"] != "ok" for r in runs))
    one_line = lambda e: str(e).strip().splitlines()[-1][:200] if e else None  # noqa: E731
    cols = {k: [r.get(k) for r in runs] for k in ("config", "status", "elapsed_s", "output_dir")}
    cols["error"] = [one_line(r.get("error")) for r in runs]
    doc.table("runs", cols, note="summary of each run: <output_dir>/results.txt", kind="records")
    doc.write(out)
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"{len(rows)} configs, {len(bad)} failed -> {out}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
