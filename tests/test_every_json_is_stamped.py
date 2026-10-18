"""Every JSON artifact records which program and version wrote it.

Paralome made this opt-in and spent months unable to answer "which code wrote
this file?" from the file. Here it is mandatory; these tests keep it that way.
"""

import json

import pytest

from evidome import jsonio
from evidome.errors import ProjectError
from evidome.identity import DISPLAY_NAME
from evidome.version import __version__


def test_write_json_stamps_version_and_generator(tmp_path):
    path = jsonio.write_json(tmp_path / "a.json", {"x": 1})
    data = json.loads(path.read_text())
    assert data[jsonio.VERSION_KEY] == __version__
    assert data[jsonio.GENERATOR_KEY] == DISPLAY_NAME


def test_write_json_does_not_mutate_the_payload(tmp_path):
    payload = {"x": 1}
    jsonio.write_json(tmp_path / "a.json", payload)
    assert payload == {"x": 1}


def test_output_is_deterministic_and_newline_terminated(tmp_path):
    a = jsonio.write_json(tmp_path / "a.json", {"b": 2, "a": 1}).read_text()
    b = jsonio.write_json(tmp_path / "b.json", {"a": 1, "b": 2}).read_text()
    assert a == b and a.endswith("\n")


def test_empty_version_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(jsonio, "__version__", "")
    with pytest.raises(ProjectError) as info:
        jsonio.write_json(tmp_path / "a.json", {"x": 1})
    assert info.value.code == 2


def test_a_rewritten_file_names_the_version_that_rewrote_it(tmp_path):
    """Read back, edited and written again, a file carries the stamp of the version that last wrote it."""
    path = tmp_path / "reference.json"
    jsonio.write_json(path, {jsonio.VERSION_KEY: "0.0.1", "species": {}})
    assert jsonio.read_json(path)[jsonio.VERSION_KEY] == jsonio.__version__


def test_no_temp_file_left_behind(tmp_path):
    jsonio.write_json(tmp_path / "a.json", {"x": 1})
    assert [p.name for p in tmp_path.iterdir()] == ["a.json"]


def test_corrupt_json_is_a_numbered_error(tmp_path):
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(ProjectError) as info:
        jsonio.read_json(bad)
    assert info.value.code == 1
