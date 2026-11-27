"""AI-CAE4ALL-style launcher: flat ``key value`` configs, preflight diagnostics,
route specs and a fixed output contract (see ``docs/CONFIGURATION.md`` and
``docs/OUTPUT_FORMAT.md``)."""
from .diagnostics import Diagnostic, DiagnosticReport, Severity, render_report
from .parser import parse_config, parse_text, write_config
from .specs import ROUTES, resolve

__all__ = ["Diagnostic", "DiagnosticReport", "Severity", "render_report", "parse_config", "parse_text",
           "write_config", "ROUTES", "resolve", "run_config"]


def run_config(path, output_dir=None, verbose=True, overrides=None):
    """Preflight + run a config file from Python; returns the results.txt content as a dict.

    ``overrides`` (dict) replaces resolved values after the preflight, e.g. ``{"plots": False}``.
    """
    from .runner import run
    pc = parse_config(path)
    vals, rep = resolve(pc)
    if rep.has_errors():
        raise ValueError(render_report(rep))
    if output_dir:
        vals["output_dir"] = str(output_dir)
    vals.update(overrides or {})
    return run(vals, pc.base_dir, verbose=verbose, config_path=pc.source_path)
