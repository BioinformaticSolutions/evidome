"""``run``: every step in order, the way Paralome runs a pipeline.

    crosswalk  ->  enrich  ->  regulators  ->  report  ->  figures

Locally, one after another on this machine. With ``--backend slurm``, the
crosswalk runs here (about a minute), then three jobs are chained with
``afterok`` dependencies: the enrichment array (one task per comparison), its
finalize job, and a ``run-finish`` job that runs regulators and report. The
command returns with the job IDs; the evidence bundle and workbook appear
when the last job ends. Every step is the same code the individual commands
run, with the same options, so a chained run and a step-by-step run write the
same files.
"""

from __future__ import annotations

import shlex
import sys
from dataclasses import replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence


def _figures(crosswalk_dir: Path, figure_settings, threads, log) -> List[str]:
    """The last step; never fails the run for want of matplotlib (the evidence is complete without it)."""
    if figure_settings is None:
        return []
    from .figures import make_figures

    result = make_figures(crosswalk_dir, figure_settings, threads=threads, log=log)
    return [f"wrote {result['figures']} figures in {result['dir']}"] if result.get("figures") else \
           [f"figures: {result.get('skipped') or 'nothing to draw'}"]


def run_local(crosswalk_dir: Path, reference_dir: Path, *, release: int, settings, sources, regulator_options: Dict,
              threads: Optional[int], figure_settings=None, log=print) -> List[str]:
    from .cli import enrich_lines
    from .enrich import enrich
    from .regulators import regulator_lines, run_regulators
    from .report import build_report

    lines = enrich_lines(enrich(crosswalk_dir, reference_dir, release=release, settings=settings, sources=sources,
                                threads=threads, log=log))
    lines += regulator_lines(run_regulators(crosswalk_dir, reference_dir, release=release, alpha=settings.alpha,
                                            threads=threads, log=log, **regulator_options))
    paths = build_report(crosswalk_dir, log=log)
    lines += [f"wrote {paths[k]}" for k in ("evidence", "archive", "workbook")]
    return lines + _figures(crosswalk_dir, figure_settings, threads, log)


def submit_chain(crosswalk_dir: Path, reference_dir: Path, *, options, release: int, settings, sources,
                 regulator_options: Dict, figure_settings=None, submit=None, log=print) -> Dict[str, str]:
    """Enrichment array + finalize (enrich.plan_slurm), then run-finish after them."""
    from .enrich import plan_slurm
    from .identity import NAME
    from .slurm import script_text, submit as sbatch_submit

    submit = submit or sbatch_submit
    ids = plan_slurm(crosswalk_dir, reference_dir, options=options, release=release, settings=settings,
                     sources=sources, submit=submit, extra={"regulators": dict(regulator_options),
                            "figures": figure_settings.as_dict() if figure_settings else None}, log=log)
    work = Path(ids["plan"]).parent
    n = len(_comparisons_in(ids["plan"]))
    finish_options = replace(options, cpus_per_task=max(1, min(options.cpus_per_task, n)), mem="16G")
    text = script_text(f"{NAME}-finish", finish_options, work / "logs" / "finish.%j.log",
                       f"{shlex.quote(sys.executable)} -m {NAME}.cli run-finish --plan {shlex.quote(str(ids['plan']))}")
    (work / "finish.sbatch").write_text(text)
    ids["finish"] = submit(work / "finish.sbatch", dependency=f"afterok:{ids['finalize']}")
    log(f"run-finish job {ids['finish']} (regulators, report, figures) after finalize {ids['finalize']}")
    return ids


def _comparisons_in(plan_path) -> List[str]:
    import json

    return [c["name"] for c in json.loads(Path(plan_path).read_text())["comparisons"]]


def run_finish(plan_path: Path, *, threads: Optional[int] = None, log=print) -> List[str]:
    """(run by SLURM) regulators, report and figures for a submitted plan."""
    from .enrich import _read_plan
    from .regulators import regulator_lines, run_regulators
    from .report import build_report

    plan, settings, _ = _read_plan(plan_path)
    crosswalk_dir, reference_dir = Path(plan["crosswalk_dir"]), Path(plan["reference_dir"])
    lines = regulator_lines(run_regulators(crosswalk_dir, reference_dir, release=plan["release"], alpha=settings.alpha,
                                           threads=threads, log=log, **(plan.get("regulators") or {})))
    paths = build_report(crosswalk_dir, log=log)
    lines += [f"wrote {paths[k]}" for k in ("evidence", "archive", "workbook")]
    from .figures import FigureSettings

    figures = plan.get("figures", {})
    return lines + _figures(crosswalk_dir, None if figures is None else FigureSettings.from_dict(figures), threads, log)
