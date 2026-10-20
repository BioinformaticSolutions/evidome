"""1.1: EVD-0015 says whether the account is outside the reference's group, or
the folder is not group-writable, instead of always asking for the modes to be
changed -- they can be right while the account is simply not in the group.

Groups and writability are injected: the suite may run as root.
"""

from __future__ import annotations

import os

from evidome.reference import unwritable_reason


def test_outside_the_group_names_the_group_not_the_modes(tmp_path):
    folder = tmp_path / "evidome"
    folder.mkdir(mode=0o2775)
    gid = folder.stat().st_gid
    text = unwritable_reason(folder / "registered_genomes" / "Mouse_GRCm39.e112",
                             groups=[gid + 1], can_write=False)
    assert "is not in group" in text and "usermod -aG" in text
    assert "modes are not the problem" in text
    assert str(folder) in text                       # the nearest folder that exists


def test_in_the_group_but_not_group_writable_names_the_chmod(tmp_path):
    folder = tmp_path / "evidome"
    folder.mkdir()
    folder.chmod(0o755)
    text = unwritable_reason(folder / "x", groups=[folder.stat().st_gid], can_write=False)
    assert "not group-writable" in text and "0755" in text and "chmod -R g+rwX" in text


def test_neither_leaves_the_general_text(tmp_path):
    folder = tmp_path / "evidome"
    folder.mkdir()
    folder.chmod(0o2775)
    assert unwritable_reason(folder, groups=[folder.stat().st_gid], can_write=False) == ""
    assert unwritable_reason(folder, groups=[], can_write=True) == ""


def test_the_run_stops_with_the_group_named(tmp_path, monkeypatch):
    from evidome import reference
    from evidome.errors import ProjectError

    blocked = tmp_path / "shared"
    blocked.mkdir(mode=0o2775)
    real = reference.unwritable_reason
    monkeypatch.setattr(reference, "unwritable_reason",
                        lambda path: real(path, groups=[os.stat(blocked).st_gid + 1], can_write=False))
    monkeypatch.setattr(reference.Path, "mkdir",
                        lambda *a, **k: (_ for _ in ()).throw(PermissionError(13, "Permission denied")))
    try:
        reference.ensure_writable(blocked / "evidome")
    except ProjectError as raised:
        assert "EVD-0015" in str(raised) and "is not in group" in str(raised)
    else:
        raise AssertionError("ensure_writable did not stop")
