"""The one JSON writer and reader. Every JSON artifact is version-stamped.

Paralome lesson, twice over:

1. JSON writing existed in several near-identical copies before it was pulled
   into one module. One implementation here from day one.
2. Paralome made version stamping OPT-IN (``stamp=True``) to avoid disturbing
   existing results, and so for months "which code wrote this?" could not be
   answered from the artifact itself. Here there are no legacy results, so
   stamping is MANDATORY: ``write_json`` refuses to write an unstamped payload
   and there is no flag to turn that off.

The stamp is two keys, ``generated_by`` (the display name) and
``program_version``, so which code wrote an output file can be read from the
file itself.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Union

from .errors import ProjectError
from .identity import DISPLAY_NAME
from .version import __version__

VERSION_KEY = "program_version"
GENERATOR_KEY = "generated_by"

PathLike = Union[str, Path]


def stamp(payload: Mapping[str, Any], *, version: Optional[str] = None) -> Dict[str, Any]:
    """``payload`` plus who wrote it. Returns a copy; never mutates."""
    stamped = dict(payload)
    stamped[GENERATOR_KEY] = DISPLAY_NAME
    stamped[VERSION_KEY] = version if version is not None else __version__
    return stamped


def dumps(payload: Mapping[str, Any]) -> str:
    """Deterministic, newline-terminated JSON (sorted keys, 2-space indent)."""
    return json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"


def _read_umask() -> int:
    try:
        with open("/proc/self/status", encoding="ascii") as handle:
            for line in handle:
                if line.startswith("Umask:"):
                    return int(line.split()[1], 8)
    except (OSError, ValueError):
        pass
    mask = os.umask(0o022)      # read once at import, before any thread exists
    os.umask(mask)
    return mask


#: mkstemp creates files readable by their owner only; JSON is written with the permissions an ordinary file
#: would get, so a shared reference directory's manifests are readable by everyone using it.
_FILE_MODE = 0o666 & ~_read_umask()


def write_json(path: PathLike, payload: Mapping[str, Any]) -> Path:
    """Stamp, then write atomically. The ONLY way JSON reaches disk.

    Atomic: written to a sibling temp file and ``os.replace``d, so a reader --
    another stage, the supervisor, a person with ``cat`` -- never sees half a
    checkpoint.
    """
    target = Path(path)
    # Always re-stamped: a payload read back, edited and written again names the version that wrote it last, so
    # reference.json names the release that last indexed into it, not the one that created it.
    stamped = stamp(payload)
    if not stamped.get(VERSION_KEY):
        raise ProjectError(2, f"refusing to write {target} with an empty {VERSION_KEY}")
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(dumps(stamped))
        os.chmod(tmp, _FILE_MODE)
        os.replace(tmp, target)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return target


def read_json(path: PathLike) -> Dict[str, Any]:
    """Read a JSON artifact. A corrupt file is a numbered error, not a traceback."""
    target = Path(path)
    try:
        return json.loads(target.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ProjectError(1, f"{target}: {exc}") from exc
    except FileNotFoundError:
        raise
    except (OSError, UnicodeDecodeError) as exc:     # a clean error, not a raw traceback
        raise ProjectError(1, f"{target}: {getattr(exc, 'strerror', None) or exc}") from exc
