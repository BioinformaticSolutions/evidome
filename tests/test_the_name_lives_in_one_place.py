"""The project name must not leak into source files outside identity.py.

The program's name is defined in one place, identity.py. Inside src/ the
program learns its own name from identity.py and imports itself relatively; a
literal copy of the name anywhere else in src/ would be a second definition
that could drift from the first.

pyproject.toml and the package directory carry the same name and are exempt;
so are tests (which import the package by name) and docs.
"""

from pathlib import Path

from evidome import identity

SRC = Path(identity.__file__).resolve().parent


def test_name_appears_only_in_identity():
    offenders = []
    for path in SRC.rglob("*.py"):
        if path.name == "identity.py":
            continue
        text = path.read_text(encoding="utf-8")
        for token in (identity.NAME, identity.DISPLAY_NAME, identity.ENV_PREFIX):
            if token in text:
                offenders.append(f"{path.name}: {token}")
    assert not offenders, f"project name hard-coded outside identity.py: {offenders}"


def test_src_imports_itself_relatively():
    offenders = [
        path.name
        for path in SRC.rglob("*.py")
        if path.name != "identity.py"  # its docstring describes the rule
        and (
            f"import {identity.NAME}" in path.read_text(encoding="utf-8")
            or f"from {identity.NAME}" in path.read_text(encoding="utf-8")
        )
    ]
    assert not offenders, f"absolute self-imports (use relative imports): {offenders}"


def test_package_directory_matches_identity():
    assert SRC.name == identity.NAME
