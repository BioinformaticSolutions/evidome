"""Tests never see the machine they run on.

The program asks Paralome, the sister program, where its shared root and genome
registry are. On a machine where Paralome is installed and set up, a test would
otherwise silently pick up that machine's real genomes. Every test starts on a
machine with no Paralome installed unless it builds its own answer.
"""

import pytest


@pytest.fixture(autouse=True)
def no_paralome(monkeypatch):
    from evidome import sister

    report = sister.SisterReport()          # installed=False: "not installed"
    monkeypatch.setattr(sister, "_CACHED", report)
    monkeypatch.setattr(sister, "sister", lambda refresh=False: report)
    return report


@pytest.fixture(autouse=True)
def no_network_to_ncbi(monkeypatch, tmp_path):
    """NCBI's real server is never contacted from a test; the crosswalk mirror points this at its own copy."""
    from evidome import reference

    monkeypatch.setattr(reference, "DEFAULT_NCBI_BASE_URL", (tmp_path / "no-ncbi-mirror").as_uri())


@pytest.fixture(autouse=True)
def no_network_to_ensembl_rest(monkeypatch):
    """Ensembl's REST service is never contacted from a test (the crosswalk asks it which species an
    unknown input is); a test that needs an answer replaces reference.species_of_ids."""
    from evidome import reference

    def refuse(*args, **kwargs):
        from evidome.errors import ProjectError

        raise ProjectError(10, "Ensembl REST is not contacted from tests")

    monkeypatch.setattr(reference, "species_of_ids", refuse)
