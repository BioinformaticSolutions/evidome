"""1.1: a standalone install knows nothing of Paralome's shared genomes or indexes.

Without Paralome the user says where the indexes live, or should be built; given nothing, the run stops and
names --reference-dir and --gtf. --gtf is the run's annotation, and a standalone --reference-dir is remembered
for the account's later runs.
"""

from __future__ import annotations

import pytest

from evidome import reference
from evidome.cli import main
from evidome.crosswalk import build_crosswalk, ensure_genome_species
from evidome.errors import ProjectError
from evidome.reference import fetch, resolve_reference_dir

from test_registered_genomes_identify_the_input import IDS, _fetch, _gtf, mirror  # noqa: F401 (fixture)
from test_reference_and_crosswalk import M


def test_nothing_given_stops_and_names_what_to_give():
    with pytest.raises(ProjectError) as info:
        resolve_reference_dir(None)
    assert info.value.code == 13
    assert "--reference-dir" in str(info.value) and "--gtf" in str(info.value)


def test_a_standalone_reference_dir_is_remembered(tmp_path):
    first, _ = resolve_reference_dir(tmp_path / "my_indexes")
    again, source = resolve_reference_dir(None)
    assert again == (tmp_path / "my_indexes").resolve() and "last --reference-dir" in source
    other, _ = resolve_reference_dir(tmp_path / "elsewhere")
    assert resolve_reference_dir(None)[0] == (tmp_path / "elsewhere").resolve()


def test_the_run_says_so_on_the_command_line(tmp_path, capsys):
    table = tmp_path / "t.csv"
    table.write_text("Gene,P-value,Q-value,Fold_change\nENSMMUG00000000001,0.01,0.1,2\n")
    assert main([str(table), "--out", str(tmp_path / "o")]) == 2
    err = capsys.readouterr().err
    assert "EVD-0013" in err and "--reference-dir" in err and "--gtf" in err


def test_gtf_on_the_run_is_the_inputs_genome_with_no_registry(mirror, tmp_path):
    base, ref = mirror
    _fetch(ref, base)
    gtf = _gtf(tmp_path / "anywhere/Macaca_mulatta.Mmul_10.113.gtf", [(g, "", "protein_coding") for g in IDS])
    genome = ensure_genome_species(IDS, ref, human_rel=116, gtf=gtf, fetch_species=False, log=lambda *_: None)
    assert genome.name == "Macaca_mulatta.Mmul_10.113" and genome.source == "--gtf"
    call, _ = build_crosswalk(IDS, [None] * len(IDS), ref, use_registry=False, gtf=gtf)
    assert call.genome.name == "Macaca_mulatta.Mmul_10.113" and call.orthology_release == 113
    assert (ref / "registered_genomes" / "Macaca_mulatta.Mmul_10.113" / "genes.tsv").exists()


def test_a_gtf_that_does_not_hold_the_input_is_refused(mirror, tmp_path):
    base, ref = mirror
    _fetch(ref, base)
    gtf = _gtf(tmp_path / "Macaca_mulatta.Mmul_10.113.gtf", [("ENSMMUG00000099999", "", "protein_coding")])
    with pytest.raises(ProjectError) as info:
        ensure_genome_species(IDS, ref, human_rel=116, gtf=gtf, fetch_species=False, log=lambda *_: None)
    assert info.value.code == 14 and "--gtf" in str(info.value)


def _stub_reference_steps(monkeypatch, calls):
    """Each reference step records that it ran and leaves the file the run checks for."""
    from evidome import genesets, pathway_signs, reference, regulators
    from evidome.reference import release_dir

    def touch(path):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n")

    monkeypatch.setattr(reference, "fetch", lambda ref, species, release=None, **_: (
        calls.append(("fetch", release)), touch(release_dir(ref, release) / "human_genes.tsv")))
    monkeypatch.setattr(genesets, "fetch_genesets", lambda ref, release=None, **_: (
        calls.append(("genesets", release)), touch(release_dir(ref, release) / "genesets/gene_sets.tsv")))
    monkeypatch.setattr(regulators, "fetch_regulators", lambda ref, release=None, **_: (
        calls.append(("regulators", release)), touch(release_dir(ref, release) / "regulators/regulons.tsv")))
    monkeypatch.setattr(pathway_signs, "fetch_pathway_signs", lambda ref, release=None, **_: (
        calls.append(("pathway-signs", release)), touch(release_dir(ref, release) / "genesets/pathway_signs.json")))
    from evidome import releases

    monkeypatch.setattr(releases, "newest_release", lambda base_url=None: 116)


def _args(**kw):
    from types import SimpleNamespace

    return SimpleNamespace(**{"release": None, "no_fetch_species": False, "threads": 1, **kw})


def test_the_run_builds_what_an_empty_reference_lacks_in_order(tmp_path, monkeypatch):
    from evidome.cli import _ensure_reference

    calls = []
    _stub_reference_steps(monkeypatch, calls)
    _ensure_reference(tmp_path / "ref", _args(), lambda *_: None)
    assert calls == [("fetch", 116), ("genesets", 116), ("regulators", 116), ("pathway-signs", 116)]
    calls.clear()
    _ensure_reference(tmp_path / "ref", _args(), lambda *_: None)      # a later run: nothing to build
    assert calls == []


def test_offline_the_run_reports_what_is_missing_and_fetches_nothing(tmp_path, monkeypatch):
    from evidome.cli import _ensure_reference

    calls = []
    _stub_reference_steps(monkeypatch, calls)
    with pytest.raises(ProjectError) as info:
        _ensure_reference(tmp_path / "ref", _args(no_fetch_species=True), lambda *_: None)
    assert info.value.code == 19 and calls == []
