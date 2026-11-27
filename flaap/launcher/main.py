"""FLAAP launcher (AI-CAE4ALL style).

    python FLAAP_main.py --config configs/flutter/ex1/config_goland_pk.txt
    python FLAAP_main.py --config <cfg> --check            # preflight only
    python FLAAP_main.py --config <cfg> --dry-run          # preflight + resolved config, no solve
    python FLAAP_main.py --config <cfg> --explain-config   # every key, value and source
    python FLAAP_main.py --list-models
    python FLAAP_main.py --describe flutter
    python FLAAP_main.py --audit-configs configs           # preflight every config under a folder
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from .diagnostics import Severity, render_report
from .parser import format_value, parse_config
from .specs import ROUTES, describe, resolve


def _color(no_color: bool) -> bool:
    return (not no_color) and sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def preflight(path, strict=False):
    pc = parse_config(path)
    vals, rep = resolve(pc)
    return pc, vals, rep


def explain(pc, vals) -> str:
    lines = [f"{'key':<26s}{'value':<44s}source"]
    for k in sorted(vals):
        if k in pc.values:
            src = f"config line {pc.locations[k].line}"
        else:
            src = "default"
        lines.append(f"{k:<26s}{format_value(vals[k])[:42]:<44s}{src}")
    return "\n".join(lines)


def audit(folder: Path, strict: bool, color: bool) -> int:
    files = sorted(folder.rglob("config_*.txt"))
    bad = 0
    for f in files:
        _, _, rep = preflight(f, strict)
        e = rep.count(Severity.ERROR, strict=strict)
        w = rep.count(Severity.WARNING, strict=strict)
        flag = "FAIL" if e else "ok  "
        bad += bool(e)
        print(f"  {flag} {f.relative_to(folder)}  ({e} errors, {w} warnings)")
        if e:
            for d in rep.sorted(strict=strict):
                if d.effective_severity(strict) is Severity.ERROR:
                    print(f"        [{d.code}] {d.message}")
    print(f"{len(files)} configs, {bad} failing")
    return 1 if bad else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="FLAAP_main.py", description="FLAAP aeroelastic launcher (flat key-value configs)")
    ap.add_argument("--config", "-c", help="flat text config file")
    ap.add_argument("--check", action="store_true", help="preflight only")
    ap.add_argument("--dry-run", action="store_true", help="preflight and print the resolved config")
    ap.add_argument("--explain-config", action="store_true", help="show each key, value and source")
    ap.add_argument("--show-defaults", action="store_true", help="list defaulted keys in the report")
    ap.add_argument("--strict", action="store_true", help="promote selected warnings to errors")
    ap.add_argument("--json-report", help="write the preflight report as JSON")
    ap.add_argument("--list-models", action="store_true")
    ap.add_argument("--describe", metavar="MODEL")
    ap.add_argument("--audit-configs", metavar="DIR")
    ap.add_argument("--output-dir", help="override output_dir")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--no-color", action="store_true")
    a = ap.parse_args(argv)
    color = _color(a.no_color)
    if a.list_models:
        for n, s in ROUTES.items():
            print(f"{n:<10s} {s.description}")
            print(f"{'':<10s} modes: {', '.join(s.modes)}")
        return 0
    if a.describe:
        if a.describe not in ROUTES:
            print(f"unknown model {a.describe!r}; known: {', '.join(ROUTES)}")
            return 2
        print(describe(a.describe))
        return 0
    if a.audit_configs:
        return audit(Path(a.audit_configs), a.strict, color)
    if not a.config:
        ap.print_help()
        return 2
    pc, vals, rep = preflight(a.config, a.strict)
    if a.output_dir:
        vals["output_dir"] = a.output_dir
    print(render_report(rep, strict=a.strict, show_notices=a.show_defaults or a.explain_config, color=color))
    if a.json_report:
        rep.write_json(a.json_report, strict=a.strict)
    if a.explain_config and not rep.has_errors(strict=a.strict):
        print()
        print(explain(pc, vals))
    if rep.has_errors(strict=a.strict):
        return 1
    if a.check or a.explain_config:
        return 0
    if a.dry_run:
        print()
        for k in sorted(vals):
            print(f"{k:<26s}{format_value(vals[k])}")
        return 0
    from .runner import run
    res = run(vals, pc.base_dir, verbose=not a.quiet and vals.get("verbose", True), config_path=pc.source_path)
    from ..io.results_txt import TxtDoc
    doc = TxtDoc("results")
    doc.kv("status", res["status"])
    doc.kv("elapsed_s", res["elapsed_s"])
    doc.tree(["summary"], res["summary"])
    txt = doc.text().split("\n", 6)[-1]               # drop the file banner and the kind line
    print(txt[:4000] + (" ...\n" if len(txt) > 4000 else ""), end="")
    print(f"% full results: {Path(res['output_dir']) / 'results.txt'}")
    return 0 if res["status"] == "ok" else 3


if __name__ == "__main__":
    sys.exit(main())
