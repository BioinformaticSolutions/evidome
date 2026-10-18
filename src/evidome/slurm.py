"""Submitting work to SLURM.

Modelled on Paralome's ``pipeline/slurm.py`` -- same shape of options, same
lesson about quoting -- but written here, not imported: the two programs share
ideas, never code. Only ``sbatch`` is called, with ``--parsable`` so
the job ID is the whole answer.

Lesson carried over from Paralome: every value that can contain a space (log
paths above all) is shell-quoted in the #SBATCH lines, because SLURM splits
those lines on whitespace like a shell does -- an unquoted space there loses
the job before it writes a single log line.
"""

from __future__ import annotations

import shlex
import shutil
import subprocess
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional

from .errors import ProjectError


@dataclass(frozen=True)
class SlurmOptions:
    cpus_per_task: int = 40
    mem: str = "4G"   # the largest comparison measured peaked at 0.73 GB
    partition: Optional[str] = None
    account: Optional[str] = None
    qos: Optional[str] = None
    nodelist: Optional[str] = None
    time: Optional[str] = None

    def small(self) -> "SlurmOptions":
        """The same placement for a light follow-up job (finalize): 1 CPU, 8G."""
        return replace(self, cpus_per_task=1, mem="8G")


def script_text(job_name: str, options: SlurmOptions, log: Path, command: str, *, array: Optional[str] = None) -> str:
    q = shlex.quote
    lines = ["#!/usr/bin/env bash",
             f"#SBATCH --job-name={q(job_name)}",
             f"#SBATCH --cpus-per-task={options.cpus_per_task}",
             f"#SBATCH --mem={q(options.mem)}",
             f"#SBATCH --output={q(str(log))}",
             f"#SBATCH --error={q(str(log))}"]
    for flag, value in (("partition", options.partition), ("account", options.account), ("qos", options.qos),
                        ("nodelist", options.nodelist), ("time", options.time)):
        if value:
            lines.append(f"#SBATCH --{flag}={q(value)}")
    if array:
        lines.append(f"#SBATCH --array={array}")
    lines += ["", "set -euo pipefail", command, ""]
    return "\n".join(lines)


def submit(script: Path, *, dependency: Optional[str] = None, runner=subprocess.run) -> str:
    """``sbatch --parsable [--dependency=...] script`` -> job ID. EVD-0023 if sbatch is missing or refuses."""
    if runner is subprocess.run and shutil.which("sbatch") is None:
        raise ProjectError(23, "sbatch is not on PATH on this machine")
    args = ["sbatch", "--parsable"] + ([f"--dependency={dependency}"] if dependency else []) + [str(script)]
    result = runner(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise ProjectError(23, f"{' '.join(args)}: {(result.stderr or result.stdout).strip()}")
    return result.stdout.strip().split(";", 1)[0]
