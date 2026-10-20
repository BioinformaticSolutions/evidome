"""Ensembl 112's homology folders publish no MD5SUM or CHECKSUMS.

A mouse Mouse_GRCm39.e112 run stopped with EVD-0011 "no MD5SUM or CHECKSUMS
entry" for Compara.112.protein_default.homologies.tsv.gz, and every stratum
downloaded the 140 MB file again. A gzip file there is accepted when it
decompresses whole with a valid CRC-32 and length.
"""
import gzip
import hashlib

import pytest

from evidome.errors import ProjectError
from evidome.reference import verify


def _gz(path, text=b"gene\thomolog\n" * 1000):
    path.write_bytes(gzip.compress(text))
    return path


def test_a_gzip_in_a_folder_without_checksums_is_checked_whole_and_accepted(tmp_path):
    path = _gz(tmp_path / "Compara.112.protein_default.homologies.tsv.gz")
    result = verify(path, tmp_path.as_uri())
    size = len(b"gene\thomolog\n") * 1000
    assert result["method"] == "gzip-crc" and result["value"].endswith(f" {size}")


def test_a_truncated_gzip_there_is_refused(tmp_path):
    path = _gz(tmp_path / "x.tsv.gz")
    path.write_bytes(path.read_bytes()[:-6])
    with pytest.raises(ProjectError) as raised:
        verify(path, tmp_path.as_uri())
    assert raised.value.code == 11


def test_a_corrupted_gzip_there_is_refused(tmp_path):
    path = _gz(tmp_path / "x.tsv.gz")
    data = bytearray(path.read_bytes())
    data[-8] ^= 0xFF          # the trailer's CRC-32
    path.write_bytes(bytes(data))
    with pytest.raises(ProjectError):
        verify(path, tmp_path.as_uri())


def test_an_uncompressed_file_there_is_still_refused(tmp_path):
    path = tmp_path / "x.tsv"
    path.write_text("gene\thomolog\n")
    with pytest.raises(ProjectError):
        verify(path, tmp_path.as_uri())


def test_a_folder_that_publishes_checksums_must_list_the_file(tmp_path):
    path = _gz(tmp_path / "x.tsv.gz")
    (tmp_path / "MD5SUM").write_text("0" * 32 + "  other.tsv.gz\n")
    with pytest.raises(ProjectError):
        verify(path, tmp_path.as_uri())


def test_a_listed_file_still_uses_the_published_checksum(tmp_path):
    path = _gz(tmp_path / "x.tsv.gz")
    (tmp_path / "MD5SUM").write_text(f"{hashlib.md5(path.read_bytes()).hexdigest()}  x.tsv.gz\n")
    assert verify(path, tmp_path.as_uri())["method"] == "md5"
