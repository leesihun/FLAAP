"""Preflight diagnostics (AI-CAE4ALL style): code, severity, field, location, hint."""
from __future__ import annotations

import json
from dataclasses import dataclass, field as dc_field
from enum import Enum
from pathlib import Path
from typing import Any, Iterable


class Severity(str, Enum):
    ERROR = "error"
    WARNING = "warning"
    NOTICE = "notice"


@dataclass(frozen=True)
class SourceLocation:
    path: Path
    line: int
    column: int = 1

    def to_dict(self) -> dict[str, Any]:
        return {"path": str(self.path), "line": self.line, "column": self.column}


@dataclass(frozen=True)
class Diagnostic:
    code: str
    severity: Severity
    message: str
    field: str | None = None
    location: SourceLocation | None = None
    hint: str | None = None
    promote_in_strict: bool = False

    def effective_severity(self, strict: bool) -> Severity:
        if strict and self.promote_in_strict and self.severity is Severity.WARNING:
            return Severity.ERROR
        return self.severity

    def to_dict(self, strict: bool = False) -> dict[str, Any]:
        d: dict[str, Any] = {"code": self.code, "severity": self.effective_severity(strict).value,
                             "message": self.message}
        if self.field is not None:
            d["field"] = self.field
        if self.location is not None:
            d["location"] = self.location.to_dict()
        if self.hint is not None:
            d["hint"] = self.hint
        return d


@dataclass
class DiagnosticReport:
    diagnostics: list[Diagnostic] = dc_field(default_factory=list)

    def add(self, code: str, severity: Severity, message: str, *, field_name: str | None = None,
            location: SourceLocation | None = None, hint: str | None = None,
            promote_in_strict: bool = False) -> None:
        self.diagnostics.append(Diagnostic(code, severity, message, field_name, location, hint,
                                           promote_in_strict))

    def extend(self, diags: Iterable[Diagnostic]) -> None:
        self.diagnostics.extend(diags)

    def count(self, severity: Severity, *, strict: bool = False) -> int:
        return sum(d.effective_severity(strict) is severity for d in self.diagnostics)

    def has_errors(self, *, strict: bool = False) -> bool:
        return self.count(Severity.ERROR, strict=strict) > 0

    def sorted(self, *, strict: bool = False) -> list[Diagnostic]:
        rank = {Severity.ERROR: 0, Severity.WARNING: 1, Severity.NOTICE: 2}
        return sorted(self.diagnostics, key=lambda d: (rank[d.effective_severity(strict)],
                                                       d.location.line if d.location else 0, d.code))

    def to_dict(self, *, strict: bool = False) -> dict[str, Any]:
        return {"summary": {"errors": self.count(Severity.ERROR, strict=strict),
                            "warnings": self.count(Severity.WARNING, strict=strict),
                            "notices": self.count(Severity.NOTICE, strict=strict)},
                "diagnostics": [d.to_dict(strict) for d in self.sorted(strict=strict)]}

    def write_json(self, path, *, strict: bool = False) -> None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_dict(strict=strict), indent=2), encoding="utf-8")


_RESET = "\x1b[0m"
_COLOR = {Severity.ERROR: "\x1b[31m", Severity.WARNING: "\x1b[33m", Severity.NOTICE: "\x1b[2m"}


def render_report(report: DiagnosticReport, *, strict=False, show_notices=False, color=False) -> str:
    def paint(t, a):
        return f"{a}{t}{_RESET}" if color else t
    e = report.count(Severity.ERROR, strict=strict)
    w = report.count(Severity.WARNING, strict=strict)
    n = report.count(Severity.NOTICE, strict=strict)
    verdict = paint("FAILED", "\x1b[31m") if e else paint("PASSED", "\x1b[32m")
    lines = [f"Preflight: {verdict} ({e} errors, {w} warnings, {n} notices)"]
    diags = report.sorted(strict=strict)
    for sev, title in ((Severity.ERROR, "ERRORS"), (Severity.WARNING, "WARNINGS"), (Severity.NOTICE, "NOTICES")):
        sel = [d for d in diags if d.effective_severity(strict) is sev]
        if not sel:
            continue
        if sev is Severity.NOTICE and not show_notices:
            lines.append(f"{title}: {len(sel)} (use --show-defaults or --explain-config to expand)")
            continue
        lines += ["", paint(f"{title} ({len(sel)})", _COLOR[sev])]
        for d in sel:
            where = f" ({d.location.path.name}:{d.location.line})" if d.location else ""
            fld = f" [{d.field}]" if d.field else ""
            lines.append(f"  {paint('[' + d.code + ']', _COLOR[sev])}{fld}{where} {d.message}")
            if d.hint:
                lines.append(f"    Hint: {d.hint}")
    return "\n".join(lines)
