""""Is this stage already done?" -- answered here, and only here.

This module exists because of Paralome's most expensive recurring bug. There,
"is this stage done" was answered in four independent places (the stage's own
checkpoint check, the fan-out wrapper, the index payload filter, and the
orchestrator/supervisor) that never consulted each other. A fix at one layer
left the other three wrong, so fixes shipped one layer at a time, each
"verified", each insufficient.

The rule here: every caller that needs to know whether a stage can be skipped
-- the stage itself, any runner, any status command, any future supervisor --
calls ``reuse_decision``. Nobody reads ``complete: true`` out of a checkpoint
and draws their own conclusion.

``reuse_decision`` checks three things, in this order, and returns the FIRST
reason to recompute as a sentence (never a bare bool -- "recomputing because X"
is actionable; "recomputing" is not):

1. version   -- a checkpoint written by any other program version is stale.
                Safe because the version only bumps on a behaviour change.
2. parameters -- every parameter that affects the output must match.
3. outputs   -- every recorded output file must still exist. Deleting results
                and re-running is the routine case, not the edge case.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

from .jsonio import VERSION_KEY, read_json, write_json
from .version import __version__

COMPLETE_KEY = "complete"
PARAMS_KEY = "parameters"
OUTPUTS_KEY = "outputs"


@dataclass(frozen=True)
class ReuseDecision:
    reuse: bool
    reason: str
    checkpoint: Dict[str, Any]

    def __bool__(self) -> bool:  # `if decision:` reads as "can reuse"
        return self.reuse


def reuse_decision(
    checkpoint_path: Path,
    parameters: Mapping[str, Any],
    *,
    version: Optional[str] = None,
) -> ReuseDecision:
    """Whether the checkpoint at ``checkpoint_path`` may stand in for a re-run."""
    path = Path(checkpoint_path)
    if not path.exists():
        return ReuseDecision(False, f"no checkpoint at {path}", {})
    previous = read_json(path)
    if not isinstance(previous, dict):      # valid JSON, but not a checkpoint
        return ReuseDecision(False, f"{path} is not a checkpoint object", {})
    if not previous.get(COMPLETE_KEY):
        return ReuseDecision(False, "checkpoint is not marked complete", previous)

    running = version if version is not None else __version__
    recorded = previous.get(VERSION_KEY)
    if recorded != running:
        return ReuseDecision(False, f"written by version {recorded}, running {running}", previous)

    recorded_params = previous.get(PARAMS_KEY, {})
    for key in sorted(set(parameters) | set(recorded_params)):
        if _normalise(parameters.get(key)) != _normalise(recorded_params.get(key)):
            return ReuseDecision(
                False,
                f"parameter {key!r} changed: {recorded_params.get(key)!r} -> {parameters.get(key)!r}",
                previous,
            )

    for name, output in sorted(previous.get(OUTPUTS_KEY, {}).items()):
        if output and not Path(output).exists():
            return ReuseDecision(False, f"recorded output {name!r} no longer exists: {output}", previous)

    return ReuseDecision(True, "version, parameters and outputs all match", previous)


def write_checkpoint(
    checkpoint_path: Path,
    parameters: Mapping[str, Any],
    outputs: Mapping[str, Optional[str]],
    **extra: Any,
) -> Path:
    """Record a completed stage. ``outputs`` maps output kind -> path."""
    payload: Dict[str, Any] = dict(extra)
    payload[COMPLETE_KEY] = True
    payload[PARAMS_KEY] = {k: _normalise(v) for k, v in parameters.items()}
    payload[OUTPUTS_KEY] = {k: (str(v) if v is not None else None) for k, v in outputs.items()}
    return write_json(checkpoint_path, payload)


def _normalise(value: Any) -> Any:
    """Make recorded and requested parameters comparable after a JSON round trip."""
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, tuple):
        return [_normalise(v) for v in value]
    if isinstance(value, list):
        return [_normalise(v) for v in value]
    if isinstance(value, Mapping):
        return {str(k): _normalise(v) for k, v in value.items()}
    return value
