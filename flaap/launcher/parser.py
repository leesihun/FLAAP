"""Flat ``key value`` config parser (AI-CAE4ALL grammar).

* one ``key value`` pair per line, the first whitespace splits key and value;
* blank lines and lines starting with ``%`` are ignored, ``#`` starts an
  inline comment;
* comma- (or space-) separated values become a list, a single value stays a
  scalar;
* ``true`` / ``false`` become booleans, numbers (also ``1e-4``) become
  ``int`` / ``float``;
* keys and string values are lowercased except for *path keys*
  (:data:`PATH_KEYS`) and *text keys* (:data:`TEXT_KEYS`), which keep case;
* duplicate keys and a UTF-8 BOM are errors.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .diagnostics import DiagnosticReport, Severity, SourceLocation

_KEY_RE = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_ASSIGN_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)\s*[=:]\s*(.+)")
_INT_RE = re.compile(r"^[+-]?\d+$")
_FLOAT_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")

PATH_KEYS = frozenset({"output_dir", "structure_file", "bulk_file", "frequency_file", "geometry_file",
                       "legacy_dir", "mode_files", "gaf_file", "log_file"})
TEXT_KEYS = frozenset({"title"})
NAME_KEYS = frozenset({"clamp_group", "group_thickness"})   # case-sensitive names; numbers still converted


@dataclass(frozen=True)
class ConfigEntry:
    key: str
    raw_value: str
    value: Any
    location: SourceLocation


@dataclass
class ParsedConfig:
    source_path: Path | None
    values: dict[str, Any] = field(default_factory=dict)
    raw_values: dict[str, str] = field(default_factory=dict)
    locations: dict[str, SourceLocation] = field(default_factory=dict)
    entries: list[ConfigEntry] = field(default_factory=list)
    diagnostics: DiagnosticReport = field(default_factory=DiagnosticReport)

    def location(self, key: str) -> SourceLocation | None:
        return self.locations.get(key.lower())

    @property
    def base_dir(self) -> Path:
        return self.source_path.parent if self.source_path else Path.cwd()


def _scalar(s: str, keep_case: bool):
    if _INT_RE.match(s):
        return int(s)
    if _FLOAT_RE.match(s):
        return float(s)
    if s.lower() in ("true", "false"):
        return s.lower() == "true"
    return s if keep_case else s.lower()


def parse_names(value: str) -> Any:
    """Like :func:`parse_value` but names keep their case (gmsh physical group names)."""
    value = value.strip()
    parts = [p.strip() for p in value.split(",")] if "," in value else value.split()
    if len(parts) == 1:
        return _scalar(parts[0], True)
    return [_scalar(p, True) for p in parts if p]


def parse_value(value: str, keep_case: bool = False) -> Any:
    value = value.strip()
    if "," in value:
        parts = [p.strip() for p in value.split(",")]
    elif " " in value or "\t" in value:
        parts = value.split()
    else:
        return _scalar(value, keep_case)
    if keep_case and not all(_FLOAT_RE.match(p) for p in parts):
        return value if "," not in value else [p for p in parts]
    return [_scalar(p, keep_case) for p in parts]


def format_value(v: Any) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (list, tuple)):
        return ", ".join(format_value(x) for x in v)
    if isinstance(v, float):
        return f"{v:.10g}"
    return str(v)


def parse_text(text: str, source: Path | None = None) -> ParsedConfig:
    parsed = ParsedConfig(source)
    src = source or Path("<string>")
    if text.startswith("﻿"):
        parsed.diagnostics.add("CFG-FILE-005", Severity.ERROR, "Config starts with a UTF-8 byte-order mark.",
                               location=SourceLocation(src, 1),
                               hint="Save as UTF-8 without BOM (PowerShell 5.1 'Out-File -Encoding utf8' adds one).")
        text = text.lstrip("﻿")
    for ln, line in enumerate(text.splitlines(), start=1):
        s = line.strip()
        if not s or s.startswith("%"):
            continue
        body = s.split("#", 1)[0].strip()
        if not body:
            continue
        loc = SourceLocation(src, ln)
        m = _ASSIGN_RE.fullmatch(body)
        if m:                                   # tolerate "key = value" / "key: value"
            parts = [m.group(1), m.group(2)]
            parsed.diagnostics.add("CFG-SYNTAX-003", Severity.NOTICE,
                                   f"'{m.group(1)} {body[len(m.group(1)):].lstrip()[0]} ...' accepted; "
                                   "the FLAAP syntax is 'key value'.", location=loc)
        else:
            parts = body.split(None, 1)
        if len(parts) != 2:
            parsed.diagnostics.add("CFG-SYNTAX-001", Severity.ERROR, "Expected a key followed by a value.",
                                   location=loc, hint="Use: key value")
            continue
        raw_key, raw_value = parts
        if not _KEY_RE.fullmatch(raw_key):
            parsed.diagnostics.add("CFG-SYNTAX-002", Severity.ERROR, f"Invalid config key {raw_key!r}.",
                                   location=loc)
            continue
        key = raw_key.lower()
        if key in parsed.values:
            first = parsed.locations[key]
            parsed.diagnostics.add("CFG-DUP-001", Severity.ERROR,
                                   f"Duplicate key {key!r} (first defined on line {first.line}).",
                                   field_name=key, location=loc,
                                   hint="Each key may appear once; number repeated items (surface1, surface2 ...).")
            continue
        keep = key in PATH_KEYS or key in TEXT_KEYS
        val = parse_names(raw_value) if key in NAME_KEYS else parse_value(raw_value, keep)
        parsed.values[key] = val
        parsed.raw_values[key] = raw_value
        parsed.locations[key] = loc
        parsed.entries.append(ConfigEntry(key, raw_value, val, loc))
    return parsed


def parse_config(path) -> ParsedConfig:
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        pc = ParsedConfig(p)
        pc.diagnostics.add("CFG-FILE-001", Severity.ERROR, f"Config file does not exist: {p}",
                           hint="Pass an existing flat text config with --config.")
        return pc
    try:
        text = p.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        pc = ParsedConfig(p)
        pc.diagnostics.add("CFG-FILE-003", Severity.ERROR, f"Config is not valid UTF-8: {exc}")
        return pc
    return parse_text(text, p)


def write_config(values: dict, path, header: str | None = None, order=None) -> None:
    """Write a value dictionary back as a flat config file."""
    keys = list(order or values.keys())
    width = max((len(k) for k in keys), default=8) + 2
    lines = []
    if header:
        lines += [f"% {h}" for h in header.splitlines()]
    for k in keys:
        if k in values and not k.startswith("_"):
            lines.append(f"{k:<{width}s}{format_value(values[k])}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
