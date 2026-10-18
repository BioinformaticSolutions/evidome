"""Paralome is asked where its things are, and every outcome is reported.

Paralome is the first place asked. If Paralome is not installed, or its
`setup site` has not been run yet, the program reports which of the two it is.

These run real processes: small fake `paralome` / `paralome-dev` programs on
PATH that answer `setup where --json` the way Paralome 0.4.dev45 does, fail the
way an older Paralome does, or are absent.
"""

import json
import os
import stat

import pytest

from evidome import sister as sister_module


def _fake(bin_dir, name, body):
    path = bin_dir / name
    if body == "__ANSWER__":
        answer = bin_dir / f"{name}.json"
        answer.write_text(_ANSWERS.pop() + "\n")
        body = f'while IFS= read -r line; do echo "$line"; done < "{answer}"'
    path.write_text("#!/bin/sh\n" + body + "\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


_ANSWERS = []


def _answer(**overrides):
    """A shell body printing a `setup where` report using builtins only: PATH in
    these tests holds nothing but the fakes, so even `cat` is absent -- which is
    also what keeps a real Paralome on the test machine out of them."""
    payload = {"schema": "paralome.where.v1", "paralome_version": "0.4.dev45",
               "site_configured": True, "site_config": "/usr/share/paralome/site.json",
               "shared_root": "/storage/paralome", "genome_registry": "/storage/paralome/registry",
               "genome_registry_source": "this machine's site.json"}
    payload.update(overrides)
    _ANSWERS.append(json.dumps(payload))
    return "__ANSWER__"


@pytest.fixture
def path_with(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    monkeypatch.setenv("PATH", str(bin_dir))
    monkeypatch.delenv("EVIDOME_PARALOME", raising=False)
    return bin_dir


def test_not_installed_is_said_plainly(path_with):
    report = sister_module.ask_paralome()
    assert not report.installed and not report.usable
    assert "not installed" in report.summary()


def test_a_set_up_machine_supplies_root_and_registry(path_with):
    _fake(path_with, "paralome", _answer())
    report = sister_module.ask_paralome()
    assert report.usable and str(report.shared_root) == "/storage/paralome"
    assert str(report.genome_registry) == "/storage/paralome/registry"
    assert report.version == "0.4.dev45"


def test_setup_site_not_run_is_reported_and_nothing_is_used(path_with):
    _fake(path_with, "paralome", _answer(site_configured=False, shared_root=None,
                                         genome_registry="/usr/share/paralome/genomes.d"))
    report = sister_module.ask_paralome()
    assert report.answered and not report.configured and not report.usable
    assert "setup site` has not been run" in report.summary()


def test_an_older_paralome_that_cannot_answer_is_reported(path_with):
    _fake(path_with, "paralome", 'echo "No such command \'where\'." >&2; exit 2')
    report = sister_module.ask_paralome()
    assert report.installed and not report.answered
    assert "did not answer" in report.summary() and "0.4.dev45" in report.summary()


def test_the_development_build_answers_when_production_cannot(path_with):
    """A machine can have both: an old `paralome` and a current `paralome-dev`."""
    _fake(path_with, "paralome", "exit 2")
    _fake(path_with, "paralome-dev", _answer(paralome_version="0.4.dev45"))
    report = sister_module.ask_paralome()
    assert report.usable and report.command.endswith("paralome-dev")


def test_an_explicit_choice_wins(path_with, tmp_path, monkeypatch):
    _fake(path_with, "paralome", _answer(shared_root="/from/path"))
    chosen = _fake(tmp_path, "my-paralome", _answer(shared_root="/from/choice"))
    monkeypatch.setenv("EVIDOME_PARALOME", str(chosen))
    assert str(sister_module.ask_paralome().shared_root) == "/from/choice"


def test_a_reply_that_is_not_a_where_report_is_rejected(path_with):
    _fake(path_with, "paralome", 'echo \'{"hello": 1}\'')
    report = sister_module.ask_paralome()
    assert not report.answered and "not a `setup where` report" in report.summary()
