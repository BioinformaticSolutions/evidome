"""Ask Paralome, the sister program, where its shared things are.

Paralome is asked first. If Paralome is not there, or ``setup site`` has not
been run yet, that is reported as it is, never replaced by a guess.

Paralome resolves its genome registry through four steps (site.json, an
environment variable, the user's config, a built-in default). Re-deriving that
here would be a second implementation of one contract, free to disagree with
Paralome the first time either changed. So this module runs Paralome's own
``setup where --json`` (Paralome >= 0.4.dev45) and uses the answer -- for the
communal reference root and the genome registry alike.

Paralome is run as an external program, the way Paralome runs samtools;
it is never imported. Candidates, first answer wins: the PARALOME environment
variable (prefixed with this program's name) if set, otherwise ``paralome``
then ``paralome-dev`` on PATH.

Every outcome is a plain sentence the CLI prints, never a guess:
- Paralome is not installed;
- Paralome is installed but did not answer (older than 0.4.dev45, or it failed);
- Paralome answered, but ``paralome setup site`` has not been run on this
  machine -- nothing it reports is then used;
- Paralome answered for a set-up machine -- its shared root and registry are used.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, List, Optional

from .identity import env_var

#: Program names looked for on PATH, in order.
CANDIDATES = ("paralome", "paralome-dev")
WHERE_ARGS = ("setup", "where", "--json")
WHERE_SCHEMA_PREFIX = "paralome.where."
TIMEOUT_SECONDS = 60


@dataclass
class SisterReport:
    installed: bool = False
    answered: bool = False
    configured: bool = False
    command: Optional[str] = None
    version: Optional[str] = None
    shared_root: Optional[Path] = None
    genome_registry: Optional[Path] = None
    genome_registry_source: Optional[str] = None
    site_config: Optional[str] = None
    notes: List[str] = field(default_factory=list)

    @property
    def usable(self) -> bool:
        """True only for a set-up machine: the one case its paths are used."""
        return self.answered and self.configured

    def summary(self) -> str:
        if not self.installed:
            return "Paralome is not installed (looked for " + ", ".join(CANDIDATES) + " on PATH)"
        if not self.answered:
            return f"Paralome found ({self.command}) but did not answer `setup where`: " + "; ".join(self.notes)
        if not self.configured:
            return (f"Paralome {self.version} ({self.command}) is installed, but `paralome setup site` has not "
                    f"been run on this machine (no {self.site_config}); its locations are not used")
        return f"Paralome {self.version} ({self.command}): shared root {self.shared_root}"


def _candidates(explicit: Optional[str]) -> List[str]:
    chosen = explicit or os.environ.get(env_var("PARALOME"), "").strip() or None
    if chosen:
        return [chosen]
    found = []
    for name in CANDIDATES:
        path = shutil.which(name)
        if path and path not in found:
            found.append(path)
    return found


def _ask(command: str) -> dict:
    result = subprocess.run([command, *WHERE_ARGS], capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
    if result.returncode != 0:
        tail = (result.stderr or result.stdout).strip().splitlines()[-1:] or ["no output"]
        raise RuntimeError(f"exit {result.returncode}: {tail[0]}")
    payload = json.loads(result.stdout)
    if not isinstance(payload, dict) or not str(payload.get("schema", "")).startswith(WHERE_SCHEMA_PREFIX):
        raise RuntimeError("answer is not a `setup where` report")
    return payload


def ask_paralome(explicit: Optional[str] = None, runner: Callable[[str], dict] = _ask) -> SisterReport:
    """Run Paralome's ``setup where --json`` and turn the answer into a report."""
    report = SisterReport()
    commands = _candidates(explicit)
    if not commands:
        return report
    report.installed = True
    for command in commands:
        report.command = command
        try:
            payload = runner(command)
        except FileNotFoundError:
            report.notes.append(f"{command}: not found")
            continue
        except subprocess.TimeoutExpired:
            report.notes.append(f"{command}: no answer within {TIMEOUT_SECONDS}s")
            continue
        except (RuntimeError, ValueError, OSError) as exc:
            report.notes.append(f"{command}: {exc} (Paralome older than 0.4.dev45?)")
            continue
        report.answered = True
        report.version = payload.get("paralome_version")
        report.configured = bool(payload.get("site_configured"))
        report.site_config = payload.get("site_config")
        if payload.get("shared_root"):
            report.shared_root = Path(payload["shared_root"])
        if payload.get("genome_registry"):
            report.genome_registry = Path(payload["genome_registry"])
        report.genome_registry_source = payload.get("genome_registry_source")
        return report
    report.installed = any(Path(c).exists() or shutil.which(c) for c in commands)
    return report


_CACHED: Optional[SisterReport] = None


def sister(refresh: bool = False) -> SisterReport:
    """The answer for this process -- Paralome is asked at most once per run."""
    global _CACHED
    if _CACHED is None or refresh:
        _CACHED = ask_paralome()
    return _CACHED
