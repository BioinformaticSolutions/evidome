"""Gene-set reference data: four sources, one index of human Ensembl genes."""

import csv
import gzip
import json

import pytest

from evidome.errors import ProjectError
from evidome.genesets import fetch_genesets, load_gene_sets, parse_obo
from evidome.reference import fetch

from test_reference_and_crosswalk import M, _mirror_for

# NCBI human GeneID -> Ensembl, beside the mirror's orthology rows. 99 has no gene in Ensembl 116.
HUMAN_G2E = "".join(f"9606\t{g}\t{e}\tNM\t-\t-\t-\n" for g, e in (
    ("2", "ENSG00000000002"), ("4", "ENSG00000000004"), ("5", "ENSG00000000005.3"),
    ("10", "ENSG00000000010"), ("99", "ENSG00000099999")))
GENE2GO = ("#tax_id\tGeneID\tGO_ID\tEvidence\tQualifier\tGO_term\tPubMed\tCategory\n"
           "9606\t1\tGO:0000003\tIDA\tinvolved_in\tleaf\t-\tProcess\n"
           "9606\t2\tGO:0000002\tIEA\tinvolved_in\tmiddle\t-\tProcess\n"
           "9606\t4\tGO:0000003\tIDA\tNOT involved_in\tleaf\t-\tProcess\n"
           "9606\t5\tGO:0000009\tIDA\tinvolved_in\tleaf (by alt id)\t-\tProcess\n"
           "9606\t10\tGO:0000004\tIDA\tinvolved_in\tobsolete\t-\tProcess\n"
           "9606\t10\tGO:0000005\tIDA\tenables\tbinding\t-\tFunction\n"
           "10090\t12345\tGO:0000003\tIDA\tinvolved_in\tmouse\t-\tProcess\n")
OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: process root
namespace: biological_process

[Term]
id: GO:0000002
name: middle process
namespace: biological_process
is_a: GO:0000001 ! process root

[Term]
id: GO:0000003
name: leaf process
namespace: biological_process
alt_id: GO:0000009
relationship: part_of GO:0000002 ! middle process

[Term]
id: GO:0000004
name: old process
namespace: biological_process
is_obsolete: true

[Term]
id: GO:0000005
name: binding
namespace: molecular_function

[Typedef]
id: part_of
name: part of
"""
REACTOME_MAPPING = "".join("\t".join(r) + "\n" for r in (
    ("ENSG00000000001", "R-HSA-2", "https://x", "Child pathway", "TAS", "Homo sapiens"),
    ("ENSG00000000001", "R-HSA-1", "https://x", "Top pathway", "TAS", "Homo sapiens"),
    ("ENSG00000000002", "R-HSA-1", "https://x", "Top pathway", "IEA", "Homo sapiens"),
    ("ENSMMUG00000000001", "R-MMU-2", "https://x", "Child pathway", "IEA", "Macaca mulatta"),
    ("ENSG00000099999", "R-HSA-3", "https://x", "Orphan", "TAS", "Homo sapiens"),
))
REACTOME_PATHWAYS = "R-HSA-1\tTop pathway\tHomo sapiens\nR-HSA-2\tChild pathway\tHomo sapiens\n" \
                    "R-MMU-2\tChild pathway\tMacaca mulatta\n"
REACTOME_RELATIONS = "R-HSA-1\tR-HSA-2\nR-MMU-1\tR-MMU-2\n"
WP_GMT = ("Apoptosis%WikiPathways_20260910%WP254%Homo sapiens\thttps://www.wikipathways.org/WP254\t1\t9\t99\n"
          "Nothing mapped%WikiPathways_20260910%WP999%Homo sapiens\thttps://x\t99\n")
HALLMARK_GMT = "HALLMARK_APOPTOSIS\thttps://www.gsea-msigdb.org/x\t1\t2\t4\n"


@pytest.fixture
def sources(tmp_path):
    base, ref = _mirror_for(tmp_path)
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    ncbi = tmp_path / "ncbi"
    with gzip.open(ncbi / "gene2ensembl.gz", "rt") as h:
        g2e = h.read()
    # the ortholog fetch already copied gene2ensembl; the gene-set step must REUSE that copy
    raw_copy = ref / "ensembl-116/raw/ncbi/gene2ensembl.gz"
    with gzip.open(raw_copy, "wt") as h:
        h.write(g2e + HUMAN_G2E)
    with gzip.open(ncbi / "gene2go.gz", "wt") as h:
        h.write(GENE2GO)
    reactome = tmp_path / "reactome"
    reactome.mkdir()
    (reactome / "Ensembl2Reactome_All_Levels.txt").write_text(REACTOME_MAPPING)
    (reactome / "ReactomePathways.txt").write_text(REACTOME_PATHWAYS)
    (reactome / "ReactomePathwaysRelation.txt").write_text(REACTOME_RELATIONS)
    (tmp_path / "go-basic.obo").write_text(OBO)
    wp = tmp_path / "wp"
    wp.mkdir()
    (wp / "wikipathways-20260810-gmt-Homo_sapiens.gmt").write_text("Old%WikiPathways_20260810%WP1%Homo sapiens\tx\t1\n")
    (wp / "wikipathways-20260910-gmt-Homo_sapiens.gmt").write_text(WP_GMT)
    (wp / "wikipathways-20260910-gmt-Mus_musculus.gmt").write_text("")
    msig = tmp_path / "msigdb" / "2026.1.Hs"
    msig.mkdir(parents=True)
    (msig / "h.all.v2026.1.Hs.entrez.gmt").write_text(HALLMARK_GMT)
    urls = {"reactome": reactome.as_uri(), "go_obo": (tmp_path / "go-basic.obo").as_uri(),
            "wikipathways_dir": wp.as_uri(), "msigdb": (tmp_path / "msigdb").as_uri(), "ncbi": ncbi.as_uri()}
    return ref, urls


def _index(ref, urls, **kw):
    fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None, **kw)
    return load_gene_sets(ref)


def test_every_source_lands_as_human_ensembl_sets(sources):
    info, members = _index(*sources)
    assert {i["source"] for i in info.values()} == {"reactome", "go_bp", "go_mf", "wikipathways", "hallmark"}
    assert members["MSIGDB:HALLMARK_APOPTOSIS"] == ["ENSG00000000001", "ENSG00000000002", "ENSG00000000004"]
    assert members["WP:WP254"] == ["ENSG00000000001", "ENSG00000000009"]       # the newest release; 99 dropped
    assert info["WP:WP254"]["name"] == "Apoptosis" and info["WP:WP254"]["source_ids"] == "3"
    assert "WP:WP999" not in info and "WP:WP1" not in info


def test_go_is_propagated_through_is_a_and_part_of_without_not_or_obsolete(sources):
    _, members = _index(*sources)
    assert members["GO:0000003"] == ["ENSG00000000001", "ENSG00000000005"]      # alt id folded; NOT excluded
    assert members["GO:0000002"] == ["ENSG00000000001", "ENSG00000000002", "ENSG00000000005"]  # part_of
    assert members["GO:0000001"] == members["GO:0000002"]                       # is_a
    assert "GO:0000004" not in members and members["GO:0000005"] == ["ENSG00000000010"]


def test_reactome_is_human_only_with_top_level_category(sources):
    info, members = _index(*sources)
    assert members["REACTOME:R-HSA-2"] == ["ENSG00000000001"]
    assert info["REACTOME:R-HSA-2"]["category"] == "Top pathway"
    assert info["REACTOME:R-HSA-1"]["category"] == ""                            # a top level itself
    assert "REACTOME:R-HSA-3" not in info and not any(k.startswith("REACTOME:R-MMU") for k in info)


def test_the_manifest_records_licences_checksums_and_drops(sources):
    ref, urls = sources
    path = fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None)
    manifest = json.loads(path.read_text())
    assert manifest["licences"]["reactome/mapping"]["licence"] == "CC0 1.0"
    assert manifest["licences"]["hallmark/gmt"]["licence"] == "CC BY 4.0"
    assert manifest["files"]["go/obo"]["checksum"]["value"]
    assert manifest["counts"]["wikipathways"]["sets_with_no_mapped_gene"] == 1
    assert manifest["program_version"]


def test_downloads_are_reused_and_the_ortholog_copy_of_gene2ensembl_is_shared(sources):
    ref, urls = sources
    fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None)
    lines = []
    fetch_genesets(ref, urls=urls, threads=1, log=lines.append)
    assert all(l.split()[0] == "reused" for l in lines if l.startswith("  ") and "md5" in l)
    assert not (ref / "ensembl-116/genesets/raw/ncbi").exists()


def test_parallel_and_sequential_indexes_are_identical(sources, tmp_path):
    ref, urls = sources
    fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None)
    one = (ref / "ensembl-116/genesets/members.tsv").read_bytes()
    fetch_genesets(ref, urls=urls, threads=4, log=lambda *_: None)
    assert (ref / "ensembl-116/genesets/members.tsv").read_bytes() == one


def test_a_bad_gmt_is_prv_0018(sources, tmp_path):
    ref, urls = sources
    wp = tmp_path / "wp" / "wikipathways-20260910-gmt-Homo_sapiens.gmt"
    wp.write_text("not a gmt\n")
    with pytest.raises(ProjectError) as info:
        fetch_genesets(ref, urls=urls, threads=1, refresh=True, log=lambda *_: None)
    assert info.value.code == 18


def test_no_index_is_prv_0019(tmp_path):
    with pytest.raises(ProjectError) as info:
        load_gene_sets(tmp_path)
    assert info.value.code == 19


def test_obo_cycles_do_not_hang(tmp_path):
    obo = tmp_path / "c.obo"
    obo.write_text("[Term]\nid: GO:1\nname: a\nnamespace: biological_process\nis_a: GO:2\n\n"
                   "[Term]\nid: GO:2\nname: b\nnamespace: biological_process\nis_a: GO:1\n")
    terms, parents, _ = parse_obo(obo)
    assert parents["GO:1"] == {"GO:2"} and len(terms) == 2
