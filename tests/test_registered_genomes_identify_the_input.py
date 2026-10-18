"""The machine's genome registry identifies the input before Ensembl does.

Fixture: a registry laid out as Paralome writes it (<shared_root>/registry/<user>.json,
{"genomes": {name: {"gtf": path, ...}}}), plain uncompressed GTFs as they sit in
Paralome genome folders, and the Ensembl mirror from the crosswalk tests.
"""

import json

import pytest

from evidome import genomes, reference
from evidome.crosswalk import STATUS_NOT_IN_ENSEMBL, build_crosswalk, summarize
from evidome.reference import fetch

from test_reference_and_crosswalk import M, mirror  # noqa: F401  (fixture)


def _gtf(path, genes):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!genome-build x\n" + "\n".join(
        f'1\tensembl\tgene\t1\t9\t.\t+\t.\tgene_id "{g}"; gene_name "{s}"; gene_biotype "{b}";' for g, s, b in genes) + "\n")
    return path


@pytest.fixture
def machine(tmp_path, monkeypatch):
    root = tmp_path / "communal"
    # Paralome's own answer to `setup where` on a set-up machine.
    from evidome import sister

    report = sister.SisterReport(installed=True, answered=True, configured=True, command="paralome",
                                 version="0.4.dev45", shared_root=root, genome_registry=root / "registry",
                                 genome_registry_source="this machine's site.json")
    monkeypatch.setattr(sister, "sister", lambda refresh=False: report)
    # The genome the macaque data was quantified against: everything Ensembl
    # knows, plus two genes Ensembl 116's annotation does not have.
    macaque = [(f"ENSMMUG0000000000{i}", f"G{i}", "protein_coding") for i in (1, 2, 3, 4, 6, 7)] + [
        ("ENSMMUG00000000008", "", "lncRNA"),
        ("ENSMMUG00000000020", "NEW20", "protein_coding"),
        ("ENSMMUG00000000021", "", "processed_pseudogene")]
    mgtf = _gtf(root / "genomes/Macaca_mulatta_10.e113/Macaca_mulatta.Mmul_10.113.gtf", macaque)
    sgtf = _gtf(root / "genomes/Ovis_aries_Ramb_v2/Ovis_aries.ARS-UI_Ramb_v2.0.gtf",
                [("ENSOARG00000000001", "S1", "protein_coding")])
    (root / "registry").mkdir(parents=True)
    (root / "registry/owner.json").write_text(json.dumps({"registered_by": "owner", "genomes": {
        "Macaca_mulatta_10.e113": {"gtf": str(mgtf), "registered_by": "owner"},
        "Ovis_aries_Ramb_v2": {"gtf": str(sgtf), "registered_by": "owner"}}}))
    # A copy published by someone else must not outrank the owner's entry.
    (root / "registry/alice.json").write_text(json.dumps({"registered_by": "alice", "genomes": {
        "Macaca_mulatta_10.e113": {"gtf": "/nowhere.gtf", "registered_by": "owner"}}}))
    return root, mgtf, sgtf


IDS = [f"ENSMMUG0000000000{i}" for i in (1, 2, 3, 4, 6, 7, 8)] + ["ENSMMUG00000000020", "ENSMMUG00000000021"]


def test_registry_is_the_one_paralome_reports(machine):
    root, mgtf, _ = machine
    assert genomes.registry_dir() == root / "registry"
    assert genomes.registered_genomes()["Macaca_mulatta_10.e113"]["gtf"] == str(mgtf)  # owner wins


def test_registered_genome_names_the_input_and_ensembl_gap_is_reported(machine, mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    call, rows = build_crosswalk(IDS, [None] * len(IDS), ref)
    assert call.genome.name == "Macaca_mulatta_10.e113" and call.genome.fraction == 1.0
    assert call.species == M and call.matched == 7                  # 7/9 in Ensembl: fine, >= 50%
    by = {r["input_gene_id"]: r for r in rows}
    assert by["ENSMMUG00000000020"]["mapping_status"] == STATUS_NOT_IN_ENSEMBL
    assert by["ENSMMUG00000000021"]["input_biotype"] == "processed_pseudogene"   # biotype from the registered genome
    summary = summarize(call, rows)
    assert summary["genome"]["name"] == "Macaca_mulatta_10.e113"
    assert summary["input_genes_by_status"][STATUS_NOT_IN_ENSEMBL] == 2
    # each NCBI/symbol note is printed once: the two note lists are separate objects
    assert call.notes is not call.genome.notes
    assert not any("NCBI" in n or "symbol" in n for n in call.genome.notes)


def test_only_genomes_sharing_the_prefix_are_indexed_and_the_index_is_reused(machine, mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    build_crosswalk(IDS, [None] * len(IDS), ref)
    cache = ref / "registered_genomes"
    assert (cache / "Macaca_mulatta_10.e113/genes.tsv").exists()
    assert not (cache / "Ovis_aries_Ramb_v2").exists()                # ENSOARG: never parsed
    first = (cache / "Macaca_mulatta_10.e113/checkpoint.json").read_text()
    build_crosswalk(IDS, [None] * len(IDS), ref)
    assert (cache / "Macaca_mulatta_10.e113/checkpoint.json").read_text() == first  # reused, not rewritten


def test_a_changed_gtf_is_reindexed(machine, mirror):
    root, mgtf, _ = machine
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    build_crosswalk(IDS, [None] * len(IDS), ref)
    mgtf.write_text(mgtf.read_text() + '1\te\tgene\t1\t9\t.\t+\t.\tgene_id "ENSMMUG00000000022"; gene_biotype "lncRNA";\n')
    build_crosswalk(IDS, [None] * len(IDS), ref)
    table = genomes.read_gene_table(ref / "registered_genomes/Macaca_mulatta_10.e113/genes.tsv")
    assert "ENSMMUG00000000022" in table


def test_no_registry_falls_back_to_ensembl_with_the_strict_threshold(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    call, _ = build_crosswalk(IDS[:7], [None] * 7, ref)
    assert call.genome is None and call.species == M and any("no registered genomes" in n for n in call.notes)
    assert any("Paralome is not installed" in n for n in call.notes)   # said, not guessed


def test_withdrawn_genome_is_used_but_flagged(machine, mirror):
    root, _, _ = machine
    base, ref = mirror
    data = json.loads((root / "registry/owner.json").read_text())
    data["genomes"]["Macaca_mulatta_10.e113"]["withdrawn"] = {"severity": "fatal", "reason": "wrong strand"}
    (root / "registry/owner.json").write_text(json.dumps(data))
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    call, _ = build_crosswalk(IDS, [None] * len(IDS), ref)
    assert call.genome.withdrawn and any("WITHDRAWN" in n for n in call.genome.notes)


def test_standalone_install_uses_its_own_registry(mirror, tmp_path):
    """No Paralome on the machine -- a genome registered in the program's
    own registry identifies the input exactly as the sister registry would."""
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    gtf = _gtf(tmp_path / "anywhere/macaque.gtf",
               [(g, "", "protein_coding") for g in IDS])
    own = genomes.own_registry_dir(ref)
    own.mkdir(parents=True)
    (own / "lab.json").write_text(json.dumps({"genomes": {"macaque_lab": {"gtf": str(gtf)}}}))
    call, rows = build_crosswalk(IDS, [None] * len(IDS), ref)
    assert call.genome.name == "macaque_lab" and call.genome.source == genomes.SOURCE_OWN
    assert summarize(call, rows)["genome"]["source"] == genomes.SOURCE_OWN


def test_sister_registry_takes_precedence_on_a_shared_name(machine, mirror):
    root, mgtf, _ = machine
    base, ref = mirror
    own = genomes.own_registry_dir(ref)
    own.mkdir(parents=True)
    (own / "x.json").write_text(json.dumps({"genomes": {"Macaca_mulatta_10.e113": {"gtf": "/elsewhere.gtf"}}}))
    merged = genomes.all_registered_genomes(ref)
    assert merged["Macaca_mulatta_10.e113"]["gtf"] == str(mgtf)
    assert merged["Macaca_mulatta_10.e113"]["_source"] == genomes.SOURCE_SISTER
