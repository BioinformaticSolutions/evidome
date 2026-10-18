"""Fetch -> verify -> index -> identify species -> crosswalk, against a local mirror
laid out exactly like the real Ensembl e116 FTP site:
one directory per species, protein + ncrna homology files, MD5SUM beside them,
tab-separated with the verified header.
"""

import csv
import gzip
import hashlib
import json

import pytest

from evidome.cli import main
from evidome.crosswalk import build_crosswalk
from evidome.errors import ProjectError
from evidome.reference import bsd_sum, fetch

HEADER = ("gene_stable_id protein_stable_id species identity homology_type homology_gene_stable_id "
          "homology_protein_stable_id homology_species homology_identity dn ds goc_score wga_coverage "
          "is_high_confidence homology_id").split()

M = "macaca_mulatta"
H = "homo_sapiens"
# Shaped like the real e116 files: the macaque
# files hold (almost) no human rows -- the macaque-human pairs sit in the
# HUMAN files, usually with human first. Ensembl's README says row order is
# arbitrary, so one pair here is stored macaque-first, in the macaque file,
# AND repeated (same homology_id) in the human file.
# (gene, species, identity, type, partner, partner_species, partner_identity, goc, wga, high_conf, homology_id)
MACAQUE_PROTEIN = [
    ("ENSMMUG00000000006", M, "85.0", "ortholog_one2one", "ENSPPRG00000000006", "panthera_pardus", "84.0", "NULL", "NULL", "0", "10"),
    ("ENSMMUG00000000001", M, "99.0", "ortholog_one2one", "ENSPPRG00000000001", "panthera_pardus", "98.0", "NULL", "NULL", "0", "11"),
]
MACAQUE_NCRNA = [
    ("ENSMMUG00000000007", M, "92.0", "ortholog_one2one", "ENSG00000000007", H, "91.0", "NULL", "90.0", "1", "70"),
]
HUMAN_PROTEIN = [
    # one-to-one, human first: macaque identity is the SECOND identity column
    ("ENSG00000000001", H, "96.0", "ortholog_one2one", "ENSMMUG00000000001", M, "95.1", "100", "98.5", "1", "1"),
    # two macaque paralogs -> one human gene (many-to-one collision)
    ("ENSG00000000002", H, "88.0", "ortholog_one2many", "ENSMMUG00000000002", M, "90.0", "75", "90.0", "1", "2"),
    ("ENSG00000000002", H, "60.0", "ortholog_one2many", "ENSMMUG00000000003", M, "70.0", "NULL", "NULL", "0", "3"),
    # one macaque gene -> two human genes
    ("ENSG00000000004", H, "80.0", "ortholog_one2many", "ENSMMUG00000000004", M, "80.0", "50", "60.0", "1", "4"),
    ("ENSG00000000005", H, "77.0", "ortholog_one2many", "ENSMMUG00000000004", M, "78.0", "50", "60.0", "0", "5"),
    # a paralog to human is not an ortholog; human vs other species is noise
    ("ENSG00000000009", H, "41.0", "other_paralog", "ENSMMUG00000000006", M, "40.0", "NULL", "NULL", "0", "6"),
    ("ENSG00000000001", H, "90.0", "ortholog_one2one", "ENSPPRG00000000001", "panthera_pardus", "89.0", "NULL", "NULL", "1", "12"),
]
HUMAN_NCRNA = [
    # the same pair as MACAQUE_NCRNA, same homology_id, stored human-first
    ("ENSG00000000007", H, "91.0", "ortholog_one2one", "ENSMMUG00000000007", M, "92.0", "NULL", "90.0", "1", "70"),
]
MACAQUE_ANNOTATION = [(f"ENSMMUG0000000000{i}", "protein_coding") for i in (1, 2, 3, 4, 6, 7)] + [
    ("ENSMMUG00000000008", "lncRNA"),   # in the annotation, in no homology file at all
]
GTF_GENES = [("ENSG00000000001", "TSPAN6", "protein_coding"), ("ENSG00000000002", "PARA1", "protein_coding"),
             ("ENSG00000000004", "DUPA", "protein_coding"), ("ENSG00000000005", "DUPB", "protein_coding"),
             ("ENSG00000000007", "LINC0007", "lncRNA"), ("ENSG00000000009", "PAR9", "protein_coding"),
             ("ENSG00000000010", "H4C1", "protein_coding"),
             ("ENSG00000000011", "TWIN", "protein_coding"), ("ENSG00000000012", "TWIN", "protein_coding")]

# NCBI Gene, shaped like gene/DATA: many rows per gene (one per
# transcript), versioned Ensembl IDs on some rows, other species mixed in.
# Macaque gene 6 has only a PARALOG in Compara; NCBI calls it an ortholog of
# human gene 9. Macaque gene 1 is an NCBI ortholog too, but Compara already
# maps it, so NCBI must not add to it.
GENE2ENSEMBL = ("#tax_id\tGeneID\tEnsembl_gene_identifier\tRNA_nucleotide_accession.version\t"
                "Ensembl_rna_identifier\tprotein_accession.version\tEnsembl_protein_identifier\n"
                "9544\t700006\tENSMMUG00000000006\tXM_1.1\tENSMMUT1\tXP_1.1\tENSMMUP1\n"
                "9544\t700006\tENSMMUG00000000006\tXM_2.1\tENSMMUT2\tXP_2.1\tENSMMUP2\n"
                "9544\t700001\tENSMMUG00000000001\tXM_3.1\t-\t-\t-\n"
                "9606\t9\tENSG00000000009.12\tNM_9.1\tENST9\tNP_9.1\tENSP9\n"
                "9606\t1\tENSG00000000001\tNM_1.1\tENST1\tNP_1.1\tENSP1\n"
                "10090\t12345\tENSMUSG00000000009\tNM_x\t-\t-\t-\n")
GENE_ORTHOLOGS = ("#tax_id\tGeneID\trelationship\tOther_tax_id\tOther_GeneID\n"
                  "9606\t9\tOrtholog\t9544\t700006\n"
                  "9606\t1\tOrtholog\t9544\t700001\n"
                  "9606\t9\tOrtholog\t10090\t12345\n")


def _gz(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as h:
        h.write(text)
    return path


def _homology_text(rows):
    lines = ["\t".join(HEADER)]
    for g, sp, ident, htype, hg, hsp, hident, goc, wga, hc, hid in rows:
        lines.append("\t".join([g, g.replace("G", "T", 1), sp, ident, htype, hg, hg.replace("G", "P", 1), hsp,
                                hident, "NULL", "NULL", goc, wga, hc, hid]))
    return "\n".join(lines) + "\n"


def _species_dir(base, species, protein, ncrna):
    hom = base / "release-116/tsv/ensembl-compara/homologies" / species
    p = _gz(hom / "Compara.116.protein_default.homologies.tsv.gz", _homology_text(protein))
    n = _gz(hom / "Compara.116.ncrna_default.homologies.tsv.gz", _homology_text(ncrna))
    (hom / "MD5SUM").write_text("".join(f"{hashlib.md5(f.read_bytes()).hexdigest()}  {f.name}\n" for f in (p, n)))
    return hom


@pytest.fixture
def mirror(tmp_path):
    return _mirror_for(tmp_path)


def _mirror_for(tmp_path):
    """The mirror and an empty reference directory under ``tmp_path`` (shared with other test files)."""
    from evidome import reference

    ncbi = tmp_path / "ncbi"
    _gz(ncbi / "gene2ensembl.gz", GENE2ENSEMBL)
    _gz(ncbi / "gene_orthologs.gz", GENE_ORTHOLOGS)
    reference.DEFAULT_NCBI_BASE_URL = ncbi.as_uri()   # restored by the conftest guard after each test
    base = tmp_path / "pub"
    (base / "release-116").mkdir(parents=True)
    (base / "release-116" / "species_EnsemblVertebrates.txt").write_text(SPECIES_TABLE)
    _species_dir(base, M, MACAQUE_PROTEIN, MACAQUE_NCRNA)
    _species_dir(base, H, HUMAN_PROTEIN, HUMAN_NCRNA)
    gtf_dir = base / "release-116/gtf/homo_sapiens"
    gtf_lines = ["#!genome-build GRCh38.p14"] + [
        f'1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "{g}"; gene_version "1"; gene_name "{s}"; gene_biotype "{b}";'
        for g, s, b in GTF_GENES] + ['1\tensembl\ttranscript\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001"; transcript_id "ENST1";']
    gtf = _gz(gtf_dir / "Homo_sapiens.GRCh38.116.gtf.gz", "\n".join(gtf_lines) + "\n")
    s, b = bsd_sum(gtf)  # this directory only has CHECKSUMS: exercises the fallback
    (gtf_dir / "CHECKSUMS").write_text(f"{s} {b} {gtf.name}\n")
    # The macaque's full annotation, with the decoy variants Ensembl ships for
    # some species (e.g. sheep). It knows genes the homology files never mention.
    mgtf_dir = base / "release-116/gtf" / M
    mgtf_lines = [f'1\tensembl\tgene\t1\t100\t.\t+\t.\tgene_id "{g}"; gene_biotype "{b}";'
                  for g, b in MACAQUE_ANNOTATION]
    mgtf = _gz(mgtf_dir / "Macaca_mulatta.Mmul_10.116.gtf.gz", "\n".join(mgtf_lines) + "\n")
    for decoy in ("Macaca_mulatta.Mmul_10.116.chr.gtf.gz", "Macaca_mulatta.Mmul_10.116.abinitio.gtf.gz"):
        _gz(mgtf_dir / decoy, "")
    s, b = bsd_sum(mgtf)
    (mgtf_dir / "CHECKSUMS").write_text(f"{s} {b} {mgtf.name}\n")
    return base.as_uri(), tmp_path / "reference"


def test_fetch_indexes_only_orthologs_to_human(mirror):
    base, ref = mirror
    manifest = json.loads(fetch(ref, [M], base_url=base, log=lambda *_: None).read_text())
    counts = manifest["species"][M]["counts"]
    assert counts["orthologs_to_human"] == 6
    assert counts["duplicate_rows_skipped"] == 1            # homology_id 70, stored twice
    assert counts["rows_human_first"] == 5 and counts["rows_species_first"] == 1
    assert manifest["species"][M]["files"][f"{H}/protein"]["checksum"]["method"] == "md5"
    assert manifest["species"][M]["indexed_by"] == manifest["program_version"]
    assert counts["annotated_genes"] == 7 and counts["genes_with_any_homology"] == 6 and counts["genes"] == 7
    assert manifest["species"][M]["files"][f"{M}/gtf"]["url"].endswith("Macaca_mulatta.Mmul_10.116.gtf.gz")
    assert manifest["human"]["checksum"]["method"] == "bsd-sum"
    assert manifest["program_version"]
    genes = (ref / "ensembl-116/species" / M / "genes.txt").read_text().split()
    assert "ENSMMUG00000000006" in genes  # known gene, just no human ortholog


def test_bad_checksum_is_refused(mirror):
    base, ref = mirror
    md5 = ref.parent / "pub/release-116/tsv/ensembl-compara/homologies" / M / "MD5SUM"
    md5.write_text(md5.read_text().replace(md5.read_text()[:4], "0000", 1))
    with pytest.raises(ProjectError) as info:
        fetch(ref, [M], base_url=base, log=lambda *_: None)
    assert info.value.code == 11


def test_crosswalk_keeps_every_gene_and_exposes_ambiguity(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    ids = [f"ENSMMUG0000000000{i}" for i in (1, 2, 3, 4, 6, 7, 8, 9)]
    call, rows = build_crosswalk(ids, [None] * len(ids), ref, min_fraction=0.8)
    assert call.species == M and call.matched == 7
    by = {}
    for r in rows:
        by.setdefault(r["input_gene_id"], []).append(r)
    assert set(by) == set(ids)                                   # nobody dropped
    assert by["ENSMMUG00000000001"][0]["human_symbol"] == "TSPAN6"
    assert by["ENSMMUG00000000002"][0]["input_genes_sharing_human_gene"] == 2
    assert len(by["ENSMMUG00000000004"]) == 2 and by["ENSMMUG00000000004"][0]["human_targets_of_input_gene"] == 2
    # Compara holds only a paralog for gene 6 (never an ortholog); NCBI reaches it
    assert by["ENSMMUG00000000006"][0]["homology_source"] == "ncbi"
    assert by["ENSMMUG00000000006"][0]["human_symbol"] == "PAR9"
    assert by["ENSMMUG00000000007"][0]["homology_source"] == "ncrna"
    # annotated but without any homology: a known gene with no human ortholog, labelled by biotype
    assert by["ENSMMUG00000000008"][0]["mapping_status"] == "no human ortholog in reference"
    assert by["ENSMMUG00000000008"][0]["input_biotype"] == "lncRNA"
    assert by["ENSMMUG00000000009"][0]["mapping_status"] == "gene ID not in reference"
    assert by["ENSMMUG00000000003"][0]["goc_score"] == ""         # NULL -> empty


def test_species_below_threshold_is_prv_0014(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    with pytest.raises(ProjectError) as info:
        build_crosswalk(["ENSMMUG00000000001", "X1", "X2", "X3"], [None] * 4, ref)
    assert info.value.code == 14


def test_human_input_maps_to_itself(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    call, rows = build_crosswalk(["ENSG00000000001", "ENSG00000000007"], ["TSPAN6", None], ref)
    assert call.species == "homo_sapiens" and rows[1]["human_symbol"] == "LINC0007"


def test_cli_end_to_end(mirror, tmp_path, capsys):
    base, ref = mirror
    assert main(["reference", "fetch", "--species", M, "--reference-dir", str(ref), "--base-url", base]) == 0
    table = tmp_path / "t.csv"
    table.write_text("Gene,Common_name,A_ANOVA_p,A_ANOVA_q,A_F\n"
                     "ENSMMUG00000000001,TSPAN6,0.01,0.1,5\nENSMMUG00000000002,,0.5,0.9,1\n")
    out = tmp_path / "out"
    assert main(["crosswalk", str(table), "--out", str(out), "--reference-dir", str(ref)]) == 0
    rows = list(csv.DictReader((out / "ortholog_crosswalk.csv").open()))
    assert [r["human_symbol"] for r in rows] == ["TSPAN6", "PARA1"]
    summary = json.loads((out / "ortholog_crosswalk.json").read_text())
    assert summary["species"] == M and summary["program_version"]
    assert "macaca_mulatta" in capsys.readouterr().out


def test_missing_reference_is_prv_0013(tmp_path, capsys):
    table = tmp_path / "t.csv"
    table.write_text("Gene,A_ANOVA_p,A_ANOVA_q,A_F\nENSMMUG00000000001,0.01,0.1,5\n")
    assert main(["crosswalk", str(table), "--out", str(tmp_path / "o"), "--reference-dir", str(tmp_path / "none")]) == 2
    assert "EVD-0013" in capsys.readouterr().err


def test_parallel_and_sequential_fetch_give_identical_indexes(mirror, tmp_path):
    base, _ = mirror
    one, many = tmp_path / "seq", tmp_path / "par"
    fetch(one, [M], base_url=base, threads=1, log=lambda *_: None)
    fetch(many, [M], base_url=base, threads=8, log=lambda *_: None)
    for rel in ("ensembl-116/human_genes.tsv", f"ensembl-116/species/{M}/orthologs_to_human.tsv",
                f"ensembl-116/species/{M}/genes.txt"):
        assert (one / rel).read_bytes() == (many / rel).read_bytes()


def test_rerun_reuses_verified_downloads(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    manifest = ref / "ensembl-116/reference.json"
    data = json.loads(manifest.read_text()); del data["species"][M]; manifest.write_text(json.dumps(data))
    lines = []
    fetch(ref, [M], base_url=base, log=lines.append)
    assert any(f"reused {M}/Compara.116.protein_default" in l for l in lines)   # species shown


def test_unwritable_reference_dir_is_prv_0015_not_a_traceback(tmp_path, capsys):
    blocked = tmp_path / "blocked"
    blocked.write_text("a file, so nothing can be created beneath it")
    assert main(["reference", "fetch", "--species", M, "--reference-dir", str(blocked / "ref")]) == 2
    err = capsys.readouterr().err
    assert "EVD-0015" in err and "Traceback" not in err


def _paralome(monkeypatch, **fields):
    from evidome import sister

    report = sister.SisterReport(**fields)
    monkeypatch.setattr(sister, "sister", lambda refresh=False: report)
    return report


def test_reference_dir_follows_paralomes_shared_root(tmp_path, monkeypatch):
    """Reference data lives in the same communal place as Paralome's genomes,
    and the location is asked of Paralome, not re-derived."""
    from evidome import reference
    from evidome.identity import NAME, env_var

    _paralome(monkeypatch, installed=True, answered=True, configured=True, command="/usr/bin/paralome",
              version="0.4.dev45", shared_root=tmp_path / "communal")
    monkeypatch.setenv(env_var("REFERENCE"), str(tmp_path / "ignored"))
    path, source = reference.resolve_reference_dir()
    assert path == tmp_path / "communal" / NAME and "shared root" in source   # machine wins over env
    assert reference.resolve_reference_dir(tmp_path / "x")[0] == tmp_path / "x"  # explicit flag wins


def test_unset_up_paralome_is_reported_not_used(tmp_path, monkeypatch):
    """Installed but `setup site` never run: say so, do not use its fallbacks."""
    from evidome import reference
    from evidome.identity import env_var

    report = _paralome(monkeypatch, installed=True, answered=True, configured=False, command="paralome",
                       version="0.4.dev45", site_config="/usr/share/paralome/site.json",
                       shared_root=None, genome_registry=tmp_path / "default_registry")
    monkeypatch.delenv(env_var("REFERENCE"), raising=False)
    with pytest.raises(ProjectError) as info:
        reference.resolve_reference_dir()
    assert info.value.code == 13 and "setup site` has not been run" in str(info.value)
    assert not report.usable


def test_reference_dir_falls_back_to_env_then_errors(tmp_path, monkeypatch):
    from evidome import reference
    from evidome.identity import env_var

    monkeypatch.setenv(env_var("REFERENCE"), str(tmp_path / "env"))
    path, source = reference.resolve_reference_dir()
    assert path == tmp_path / "env" and "not installed" in source
    monkeypatch.delenv(env_var("REFERENCE"))
    with pytest.raises(ProjectError) as info:
        reference.resolve_reference_dir()
    assert info.value.code == 13 and "not installed" in str(info.value)


def test_orientation_swaps_identities(mirror):
    """Human-first rows: the macaque identity is the second column; it must land in identity_pct."""
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    _, rows = build_crosswalk(["ENSMMUG00000000001"], [None], ref, min_fraction=0.5)
    assert (rows[0]["identity_pct"], rows[0]["human_identity_pct"]) == ("95.1", "96.0")


def test_an_index_from_an_older_version_is_rebuilt_from_the_same_downloads(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    manifest = ref / "ensembl-116/reference.json"
    data = json.loads(manifest.read_text()); del data["species"][M]["indexed_by"]; manifest.write_text(json.dumps(data))
    lines = []
    fetch(ref, [M], base_url=base, log=lines.append)
    assert any("re-indexing" in l for l in lines)
    assert all("fetched" not in l for l in lines)                 # everything reused
    assert json.loads(manifest.read_text())["species"][M]["counts"]["orthologs_to_human"] == 6


SPECIES_TABLE = ("#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\t"
                 "microarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id\n"
                 "Sheep\tovis_aries\tEnsemblVertebrates\t9940\tARS-UI_Ramb_v3.0\tGCA_016772045.2\tENS01\t\t\t\tY\t\t\tovis_aries_core_116_3\t1\n"
                 "Sheep - Texel\tovis_aries_texel\tEnsemblVertebrates\t9940\tOar_v3.1\tGCA_000298735.1\tx\tY\tN\tN\tN\tN\tY\tovis_aries_texel_core_116_31\t1\n"
                 "Macaque\tmacaca_mulatta\tEnsemblVertebrates\t9544\tMmul_10\tGCA_003339765.3\tx\tY\tY\tY\tY\tY\tY\tmacaca_mulatta_core_116_10\t1\n")


def test_species_list_marks_orthology_and_filters_by_assembly(mirror, capsys):
    """A user is not expected to guess Ensembl's folder names: the list says which species have orthology."""
    base, _ = mirror
    from urllib.parse import urlparse
    from pathlib import Path as P
    (P(urlparse(base).path) / "release-116" / "species_EnsemblVertebrates.txt").write_text(SPECIES_TABLE)
    _species_dir(P(urlparse(base).path), "ovis_aries", [], [])
    from evidome.reference import list_species
    rows = {r["species"]: r for r in list_species(base_url=base)}
    assert rows["ovis_aries"]["orthology"] == "yes" and rows["ovis_aries_texel"]["orthology"] == "no"
    assert [r["species"] for r in list_species(base_url=base, text="ramb")] == ["ovis_aries"]
    assert main(["reference", "species", "sheep", "--base-url", base]) == 0
    out = capsys.readouterr().out
    assert "ARS-UI_Ramb_v3.0" in out and "ovis_aries_texel" not in out and "1 more genome" in out


def test_a_species_without_a_single_main_gtf_is_prv_0016(mirror):
    base, ref = mirror
    from urllib.parse import urlparse
    from pathlib import Path as P
    (P(urlparse(base).path) / "release-116/gtf" / M / "Macaca_mulatta.Mmul_10b.116.gtf.gz").write_bytes(b"")
    with pytest.raises(ProjectError) as info:
        fetch(ref, [M], base_url=base, log=lambda *_: None)
    assert info.value.code == 16


# ---------------------------------------------------------------------------
# the second and third sources
# ---------------------------------------------------------------------------

def test_ncbi_orthologs_are_indexed_in_either_orientation_for_this_taxon_only(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    rows = list(csv.DictReader((ref / f"ensembl-116/species/{M}/ncbi_orthologs_to_human.tsv").open(), delimiter="\t"))
    assert [(r["gene_id"], r["human_gene_id"]) for r in rows] == [
        ("ENSMMUG00000000001", "ENSG00000000001"), ("ENSMMUG00000000006", "ENSG00000000009")]  # version stripped
    manifest = json.loads((ref / "ensembl-116/reference.json").read_text())["species"][M]
    assert manifest["taxonomy_id"] == "9544" and manifest["counts"]["ncbi_orthologs_to_human"] == 2
    assert manifest["files"]["ncbi/gene2ensembl.gz"]["checksum"]["value"]


def test_compara_wins_where_it_has_an_answer(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    _, rows = build_crosswalk(["ENSMMUG00000000001"], [None], ref, min_fraction=0.5)
    assert len(rows) == 1 and rows[0]["homology_source"] == "protein"


def test_symbol_match_only_when_unambiguous_and_nothing_else_reached_the_gene(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    ids = ["ENSMMUG00000000001", "ENSMMUG00000000007", "ENSMMUG00000000008", "ENSMMUG00000000003"]
    # gene 1 is named after a DIFFERENT human gene: Compara's answer stands; 8 matches H4C1 (any case);
    # a gene named after two human genes is not matched; gene 7 has Compara, its symbol is ignored
    call, rows = build_crosswalk(ids, ["PAR9", "twin", "h4c1", None], ref, min_fraction=0.5)
    by = {r["input_gene_id"]: r for r in rows}
    assert by["ENSMMUG00000000001"]["human_symbol"] == "TSPAN6"
    assert by["ENSMMUG00000000008"]["human_symbol"] == "H4C1" and by["ENSMMUG00000000008"]["homology_source"] == "symbol"
    assert by["ENSMMUG00000000007"]["homology_source"] == "ncrna"
    _, rows = build_crosswalk(["ENSMMUG00000000008"], ["TWIN"], ref, min_fraction=0.5)
    assert rows[0]["mapping_status"] == "no human ortholog in reference"
    call, _ = build_crosswalk(["ENSMMUG00000000008"], ["TWIN"], ref, min_fraction=0.5)
    assert any("shared by several human genes" in n for n in call.notes)


def test_sources_can_be_switched_off(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    _, rows = build_crosswalk(["ENSMMUG00000000006", "ENSMMUG00000000008"], [None, "H4C1"], ref,
                              min_fraction=0.5, use_ncbi=False, use_symbols=False)
    assert all(r["mapping_status"] == "no human ortholog in reference" for r in rows)


def test_an_index_without_ncbi_still_crosswalks_and_says_so(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    (ref / f"ensembl-116/species/{M}/ncbi_orthologs_to_human.tsv").rename(ref / "moved.tsv")
    call, rows = build_crosswalk(["ENSMMUG00000000006"], [None], ref, min_fraction=0.5)
    assert rows[0]["mapping_status"] == "no human ortholog in reference"
    assert any("run `reference fetch" in n for n in call.notes)


def test_ncbi_downloads_are_kept_unless_refreshed(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    lines = []
    fetch(ref, [M], base_url=base, refresh_ncbi=True, log=lines.append)
    assert any("fetched ncbi/gene2ensembl.gz" in l for l in lines)
    data = json.loads((ref / "ensembl-116/reference.json").read_text()); del data["species"][M]["indexed_by"]
    (ref / "ensembl-116/reference.json").write_text(json.dumps(data))
    lines = []
    fetch(ref, [M], base_url=base, log=lines.append)
    assert any("reused ncbi/gene2ensembl.gz" in l for l in lines)


def test_a_truncated_ncbi_file_is_prv_0011(mirror):
    base, ref = mirror
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    raw = ref / "ensembl-116/raw/ncbi/gene_orthologs.gz"
    raw.write_bytes(raw.read_bytes()[:-8])
    data = json.loads((ref / "ensembl-116/reference.json").read_text()); del data["species"][M]["indexed_by"]
    (ref / "ensembl-116/reference.json").write_text(json.dumps(data))
    with pytest.raises(ProjectError) as info:
        fetch(ref, [M], base_url=base, threads=1, log=lambda *_: None)
    assert info.value.code == 11 and "--refresh-ncbi" in str(info.value)


def test_every_download_names_the_program(mirror, monkeypatch):
    """reactome.org refuses Python's default client name (403) and accepts a named one."""
    import urllib.request

    from evidome.identity import DISPLAY_NAME
    from evidome.version import __version__

    seen = []
    real = urllib.request.urlopen

    def spy(request, *a, **kw):
        seen.append(request.get_header("User-agent") if isinstance(request, urllib.request.Request) else None)
        return real(request, *a, **kw)

    monkeypatch.setattr(urllib.request, "urlopen", spy)
    base, ref = mirror
    fetch(ref, [M], base_url=base, threads=1, log=lambda *_: None)
    assert seen and all(ua == f"{DISPLAY_NAME}/{__version__} (bioinformatics reference-data download)" for ua in seen)
