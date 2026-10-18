"""The version has one source, and every registered error has a number and a fix."""

from pathlib import Path

import pytest

from evidome import cli, errors
from evidome.version import __version__

REPO = Path(__file__).resolve().parent.parent


def test_version_file_mirrors_version_py():
    assert (REPO / "VERSION").read_text().strip() == __version__


def test_every_error_has_a_title_and_a_fix():
    for code, spec in errors.REGISTRY.items():
        assert spec.code == code
        assert spec.title.strip() and spec.fix.strip(), errors.format_code(code)


def test_unregistered_code_cannot_be_raised():
    with pytest.raises(ValueError):
        errors.ProjectError(9999, "nope")


def test_error_message_carries_its_code():
    err = errors.ProjectError(1, "x")
    assert str(err).startswith(errors.format_code(1))


def test_version_report_says_where_the_code_came_from():
    report = "\n".join(cli.version_report())
    assert __version__ in report and "code imported from" in report


def test_version_report_carries_the_licence_notice():
    """The same AGPL section 13 notice Paralome shows."""
    report = "\n".join(cli.version_report())
    assert "GNU Affero General Public License, version 3" in report and "see NOTICE" in report


def test_an_error_raised_in_a_worker_process_arrives_as_itself():
    """Indexing runs in worker processes; a numbered error must survive the trip back."""
    import pickle

    err = pickle.loads(pickle.dumps(errors.ProjectError(11, "gene_orthologs.gz is truncated")))
    assert err.code == 11 and "truncated" in str(err) and err.fix
