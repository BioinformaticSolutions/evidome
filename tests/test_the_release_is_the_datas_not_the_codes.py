"""1.1: the Ensembl release is decided by the data, never fixed in the code.

Sheep counted against Rambouillet v2.0 (Ensembl 113) matched no species: Evidome read everything from Ensembl
116, and Ensembl 115 had replaced the Rambouillet v2.0 gene set. A registered genome now names its own release
(from its GTF's name), its orthologs come from that release, and the human genes they reach are carried by
stable ID to the human-side release the gene sets and regulators are read at.
"""

from __future__ import annotations

import csv

import pytest

from evidome import reference
from evidome.crosswalk import build_crosswalk, ensure_genome_species
from evidome.errors import ProjectError
from evidome.reference import fetch, load_manifest
from evidome.releases import UNFETCHED, human_release, newest_release, release_of_gtf, species_of_gtf

from test_registered_genomes_identify_the_input import IDS, _fetch, _gtf, machine, mirror  # noqa: F401 (fixtures)
from test_reference_and_crosswalk import M


@pytest.mark.parametrize("name,release,species", [
    ("Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.113.gtf", 113, "ovis_aries_rambouillet"),
    ("Macaca_mulatta.Mmul_10.113.gtf", 113, "macaca_mulatta"),
    ("Mus_musculus.GRCm39.112.gtf.gz", 112, "mus_musculus"),
    ("Homo_sapiens.GRCh38.81.gtf", 81, "homo_sapiens"),
    ("Macaca_mulatta.Mmul_10.116.chr.gtf.gz", 116, "macaca_mulatta"),
])
def test_a_registered_annotation_names_its_release_and_species(name, release, species):
    assert release_of_gtf(name) == release
    assert species_of_gtf(name) == species


def test_an_annotation_without_a_release_in_its_name_is_refused_not_guessed():
    with pytest.raises(ProjectError) as info:
        release_of_gtf("/genomes/Ovis_aries.ARS-UI_Ramb_v2.0.gtf", "Ovis_aries_Ramb_v2")
    assert info.value.code == 36 and "Ovis_aries_Ramb_v2" in str(info.value)


def test_no_release_is_fixed_in_the_code():
    from pathlib import Path

    import evidome

    src = Path(evidome.__file__).parent
    for path in src.glob("*.py"):
        text = path.read_text()
        assert "DEFAULT_RELEASE" not in text, path.name
        assert "release: int = 116" not in text and "default=116" not in text, path.name


def test_the_human_release_is_the_references_own_and_reading_never_downloads(mirror):
    base, ref = mirror
    assert human_release(ref) == UNFETCHED                     # empty: a reader is told, nothing is fetched
    assert newest_release(base) == 116                          # the newest the mirror publishes
    fetch(ref, [M], release=113, base_url=base, log=lambda *_: None)
    assert human_release(ref) == 113                            # what is on disk
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    assert human_release(ref) == 116                            # the newest on disk
    assert human_release(ref, "113") == 113                     # --release, when given


def test_orthologs_come_from_the_registered_genomes_release_and_reach_the_human_release(machine, mirror):
    base, ref = mirror
    _fetch(ref, base)
    call, rows = build_crosswalk(IDS, [None] * len(IDS), ref)
    assert call.genome.name == "Macaca_mulatta_10.e113"
    assert call.orthology_release == 113 and call.human_release == 116
    assert any("Ensembl 113 Compara" in n and "carried to Ensembl 116" in n for n in call.notes)
    assert any(r["human_gene_id"] for r in rows)


def test_a_human_gene_the_human_release_no_longer_has_is_left_out_and_counted(machine, mirror):
    base, ref = mirror
    _fetch(ref, base)
    _, before = build_crosswalk(IDS, [None] * len(IDS), ref)
    retired = next(r["human_gene_id"] for r in before if r["human_gene_id"])
    path = ref / "ensembl-116" / "human_genes.tsv"
    with path.open(newline="") as handle:
        kept = [row for row in csv.DictReader(handle, delimiter="\t") if row["human_gene_id"] != retired]
        fields = list(kept[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)
    call, after = build_crosswalk(IDS, [None] * len(IDS), ref)
    assert retired not in {r["human_gene_id"] for r in after}
    note = next(n for n in call.notes if "carried to Ensembl 116" in n)
    assert "0 pair(s)" not in note


def test_the_run_fetches_the_registered_genomes_own_release(machine, mirror, monkeypatch):
    base, ref = mirror
    monkeypatch.setattr(reference, "DEFAULT_BASE_URL", base)
    fetch(ref, [M], base_url=base, log=lambda *_: None)          # the human side only, as on a fresh machine
    assert not (ref / "ensembl-113").exists()
    genome = ensure_genome_species(IDS, ref, human_rel=116, log=lambda *_: None)
    assert genome.name == "Macaca_mulatta_10.e113"
    assert M in load_manifest(ref, 113)["species"]


# --- Ensembl's new platform ---------------------------------------------------

SHEEP = "ovis_aries_gca016772045v2"
SHEEP_IDS = ["ENSOARG00285000001", "ENSOARG00285000002", "ENSOARG00285000003"]
PLATFORM_HEADER = ("ref_species\tref_assembly\tquery_species\tquery_assembly\tref_gene_stable_id\tref_gene_name\t"
                   "query_gene_stable_id\tquery_gene_name\thomology_type\tquery_perc_id\tquery_perc_cov")


def _platform(tmp_path, own_rows_have_human=True):
    """The platform's FTP tree, as listed on 2026-10-22, with a sheep dataset and a human one."""
    import gzip
    import hashlib

    base = tmp_path / "ensemblorganisms"

    def homology(path, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(gzip.compress(("\n".join([PLATFORM_HEADER] + rows) + "\n").encode()))
        (path.parent / "md5sum.txt").write_text(f"{hashlib.md5(path.read_bytes()).hexdigest()}  {path.name}\n")

    sheep_rows = [f"bos_taurus\tARS-UCD1.2\t{SHEEP}\tARS-UI_Ramb_v3.0\tENSBTAG1\t\t{SHEEP_IDS[2]}\t\thomolog_bbh\t91.0\t92.0"]
    human_rows = [f"homo_sapiens\tGRCh38\t{SHEEP}\tARS-UI_Ramb_v3.0\tENSG00000000001\tTSPAN6\t{SHEEP_IDS[0]}\tTSPAN6\t"
                  "homolog_rbbh\t94.0\t100.0",
                  f"homo_sapiens\tGRCh38\t{SHEEP}\tARS-UI_Ramb_v3.0\tENSG00000000002\tPARA1\t{SHEEP_IDS[1]}\t\t"
                  "homolog_bbh\t58.0\t100.0"]
    homology(base / "GCA/016/772/045/2/ensembl/2024_12/homology/2026_04_09/homology.tsv.gz",
             sheep_rows + (human_rows if own_rows_have_human else []))
    flipped = [f"{SHEEP}\tARS-UI_Ramb_v3.0\thomo_sapiens\tGRCh38\t{SHEEP_IDS[0]}\t\tENSG00000000001\tTSPAN6\t"
               "homolog_rbbh\t94.0\t100.0"]
    homology(base / "GCA/000/001/405/29/ensembl/2026_04/homology/2026_07_13/homology.tsv.gz", flipped)
    gtf = _gtf(tmp_path / "genomes/Sheep_Ramb3/genes.gtf",
               [(g, f"S{i}", "protein_coding") for i, g in enumerate(SHEEP_IDS)])
    gtf.write_text("#!genome-build-accession GCA_016772045.2\n#!genebuild-last-updated 2024-12\n" + gtf.read_text())
    return base.as_uri(), gtf


def test_a_platform_annotation_names_its_dataset():
    assert release_of_gtf("/x/ensemblorganisms/GCA/016/772/045/2/ensembl/2024_12/geneset/genes.gtf.gz") == \
        "GCA_016772045.2_ensembl_2024_12"


def test_a_platform_annotation_is_read_from_its_headers(tmp_path):
    _, gtf = _platform(tmp_path)
    assert release_of_gtf(gtf) == "GCA_016772045.2_ensembl_2024_12"


@pytest.mark.parametrize("own_rows_have_human", [True, False])
def test_a_platform_genome_is_mapped_through_the_platforms_homologies(mirror, tmp_path, monkeypatch,
                                                                      own_rows_have_human):
    import json

    from evidome import genomes, sister

    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)                      # the human side, release 116
    platform, gtf = _platform(tmp_path, own_rows_have_human)
    monkeypatch.setattr(reference, "DEFAULT_PLATFORM_URL", platform)
    registry = tmp_path / "registry"
    registry.mkdir()
    (registry / "owner.json").write_text(json.dumps({"registered_by": "owner", "genomes": {
        "Sheep_Ramb3": {"gtf": str(gtf), "registered_by": "owner"}}}))
    monkeypatch.setattr(genomes, "registry_dir", lambda: registry)
    monkeypatch.setattr(sister, "sister", lambda refresh=False: sister.SisterReport(installed=False))
    genome = ensure_genome_species(SHEEP_IDS, ref, human_rel=116, log=lambda *_: None)
    assert genome.name == "Sheep_Ramb3"
    call, rows = build_crosswalk(SHEEP_IDS, [None] * 3, ref)
    assert call.species == SHEEP
    assert call.orthology_release == "GCA_016772045.2_ensembl_2024_12" and call.human_release == 116
    mapped = {r["input_gene_id"]: r["human_gene_id"] for r in rows if r["human_gene_id"]}
    assert mapped[SHEEP_IDS[0]] == "ENSG00000000001"
    if own_rows_have_human:
        assert mapped[SHEEP_IDS[1]] == "ENSG00000000002"
        assert next(r for r in rows if r["input_gene_id"] == SHEEP_IDS[0])["homology_type"] == "homolog_rbbh"


def test_a_release_that_does_not_gzip_its_homologies_is_read_as_published(tmp_path):
    """Release 113 publishes Compara.113.*_default.homologies.tsv uncompressed, with CHECKSUMS only; 116 gzips
    them. The folder's listing names the file (2026-10-22: 113 was asked for the .gz and answered 404)."""
    import gzip

    from test_reference_and_crosswalk import _mirror_for

    base, ref = _mirror_for(tmp_path, releases=(113,))
    from pathlib import Path
    from urllib.parse import urlparse

    root = Path(urlparse(base).path) / "release-113/tsv/ensembl-compara/homologies"
    for folder in root.iterdir():
        for gz in folder.glob("*.tsv.gz"):
            plain = gz.with_suffix("")
            plain.write_bytes(gzip.decompress(gz.read_bytes()))
            gz.unlink()
        (folder / "MD5SUM").unlink()
        from evidome.reference import bsd_sum

        (folder / "CHECKSUMS").write_text("".join(f"{a} {b} {p.name}\n" for p in folder.glob("*.tsv")
                                                  for a, b in [bsd_sum(p)]))
    manifest = fetch(ref, [M], release=113, base_url=base, log=lambda *_: None)
    import json

    data = json.loads(manifest.read_text())
    assert data["species"][M]["files"][f"{M}/protein"]["url"].endswith("Compara.113.protein_default.homologies.tsv")
    assert data["species"][M]["counts"]["orthologs_to_human"] > 0
