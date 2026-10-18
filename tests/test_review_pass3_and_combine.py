"""Input checks, reference-download errors, and many per-comparison tables woven into one (``combine``)."""

import argparse
import csv
import gzip
import io
import os
import zipfile

import pytest

from evidome import cli, combine as combine_mod, jsonio, parallel, reference
from evidome.cli import main
from evidome.collisions import resolve
from evidome.errors import ProjectError
from evidome.utility import score
from evidome.wide_stats import read_wide_statistics

from test_evidence_bundle_describes_itself import COLLECTRI, SIGNALLING, TABLE  # noqa: F401
from test_gene_sets_are_indexed_to_human_ensembl import sources  # noqa: F401  (fixture)
from test_one_command_runs_everything import OPTS, world  # noqa: F401  (fixture)

BODY = ("Gene,X_Raw_FPKM_Control_mean,X_Raw_FPKM_Treatment_mean,X_Fold_change,X_P-value,X_Q-value,X_T-statistic\n"
        "ENSMMUG1,1,2,2,0.01,0.05,{t}\n")


# --- combine ----------------------------------------------------------------------------------------------------

SC_HEADER = "Gene,Common_name,Control_mean,Treatment_mean,Fold_change,P-value,Q-value,S1__0,S2__0\n"


def _sc_folder(tmp_path):
    folder = tmp_path / "sc"
    folder.mkdir()
    (folder / "Gq_vs_WT_0_naive_cell_t_rna.csv").write_text(
        SC_HEADER + "ENSMUSG01,Aa,1.0,2.0,1.0,0.001,0.01,1,2\n.,Septin2,1,1,0,0.9,0.95,1,1\n"
                    "ENSMUSG02.3,Bb,2.0,1.0,-1.0,0.02,0.04,2,1\n")
    (folder / "Gq_vs_WT_1_naive_cell_t_rna.csv").write_text(
        SC_HEADER + "ENSMUSG02,Bb,1.5,1.0,-0.5,0.2,0.3,1,1\nENSMUSG03,,1.0,1.0,0.0,0.9,0.9,1,1\n")
    return folder


def test_a_folder_of_comparisons_becomes_one_wide_table(tmp_path):
    folder = _sc_folder(tmp_path)
    out = combine_mod.combine(combine_mod.find_tables(folder), tmp_path / "wide.csv",
                              drop=["_naive_cell_t_rna"], log=lambda *_: None)
    table = read_wide_statistics(out)
    assert sorted(c.name for c in table.comparisons) == ["Gq_vs_WT_0", "Gq_vs_WT_1"]
    assert table.gene_ids == ["ENSMUSG01", "ENSMUSG02", "ENSMUSG03"]       # '.' left out, versions joined
    rows = list(csv.DictReader(line for line in out.read_text().splitlines() if not line.startswith("#")))
    assert rows[0]["Gq_vs_WT_0_P-value"] == "0.001" and rows[0]["Gq_vs_WT_1_P-value"] == ""   # copied as text
    assert rows[1]["Gene"] == "ENSMUSG02.3" and rows[2]["Common_name"] == ""
    assert "S1__0" not in "".join(rows[0])


def test_the_whole_run_takes_a_folder(world, tmp_path):
    """Split the wide fixture into one file per comparison: the folder gives the same results as the table."""
    ref, table = world
    lines = TABLE.strip().splitlines()
    body = [line for line in lines if not line.startswith("#")]
    header, rows = body[0].split(","), [r.split(",") for r in body[1:]]
    folder = tmp_path / "per"
    folder.mkdir()
    for name in ("A", "F1"):
        keep = [i for i, c in enumerate(header) if c.startswith(f"{name}_")]
        with (folder / f"{name}.csv").open("w", newline="") as h:
            w = csv.writer(h, lineterminator="\n")
            w.writerow(header[:2] + [header[i][len(name) + 1:] for i in keep] + ["sample_1"])
            w.writerows(r[:2] + [r[i] for i in keep] + ["7"] for r in rows)
    whole, split = tmp_path / "whole", tmp_path / "split"
    assert main([str(table), "--out", str(whole), "--reference-dir", str(ref)] + OPTS) == 0
    assert main([str(folder), "--out", str(split), "--reference-dir", str(ref)] + OPTS) == 0
    assert (split / combine_mod.COMBINED_NAME).exists()
    for name in ("enrichment/A/ora.csv", "enrichment/A/gsea.csv", "regulators/A/regulators.csv"):
        assert (whole / name).read_bytes() == (split / name).read_bytes(), name


def test_names_that_collide_or_no_files_are_numbered_errors(tmp_path):
    folder = _sc_folder(tmp_path)
    with pytest.raises(ProjectError) as info:
        combine_mod.combine(combine_mod.find_tables(folder), tmp_path / "w.csv", drop=["_0_naive_cell_t_rna",
                                                                                         "_1_naive_cell_t_rna"])
    assert info.value.code == 31
    with pytest.raises(ProjectError) as info:
        combine_mod.find_tables(folder, "*.tsv")
    assert info.value.code == 31


def test_combine_and_inspect_from_the_command_line(tmp_path, capsys):
    folder = _sc_folder(tmp_path)
    assert main(["combine", str(folder), "--out", str(tmp_path / "w.csv"), "--drop", "_naive_cell_t_rna"]) == 0
    assert main(["inspect", str(folder), "--drop", "_naive_cell_t_rna"]) == 0
    assert "Gq_vs_WT_1" in capsys.readouterr().out


# --- input tables ------------------------------------------------------------------------------------------------

def test_an_infinite_statistic_is_refused_but_an_infinite_fold_change_is_not(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text(BODY.format(t="Inf"))
    with pytest.raises(ProjectError) as info:
        read_wide_statistics(path)
    assert info.value.code == 7 and "infinite" in str(info.value)
    path.write_text(BODY.format(t="3").replace(",2,0.01", ",inf,0.01"))
    assert read_wide_statistics(path).comparisons


def test_a_missing_table_is_a_numbered_error(tmp_path):
    for bad in (tmp_path / "nope.csv", tmp_path):
        with pytest.raises(ProjectError) as info:
            read_wide_statistics(bad)
        assert info.value.code == 30


def test_a_blank_line_before_the_header_is_skipped(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text("# generated by hand\n\n" + BODY.format(t="3"))
    assert read_wide_statistics(path).gene_ids == ["ENSMMUG1"]


# --- options ------------------------------------------------------------------------------------------------------

def test_a_bad_gsea_option_stops_the_run_before_the_crosswalk_writes(world, tmp_path):
    ref, table = world
    out = tmp_path / "study"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref), "--permutations", "0"]) == 2
    assert not (out / "ortholog_crosswalk.csv").exists()


@pytest.mark.parametrize("extra", [["--alpha", "2"], ["--alpha", "0"], ["--max-targets", "0"],
                                   ["--min-targets", "0"]])
def test_alpha_and_target_limits_are_checked(world, tmp_path, extra):
    ref, table = world
    assert main([str(table), "--out", str(tmp_path / "s"), "--reference-dir", str(ref)] + OPTS + extra) == 2


def test_enrich_and_regulators_take_the_crosswalk_alpha(tmp_path):
    (tmp_path / "comparison_summary.json").write_text('{"alpha": 0.01}')
    args = argparse.Namespace(alpha=None, crosswalk_dir=tmp_path)
    assert cli._alpha(args) == 0.01


# --- shared files and downloads --------------------------------------------------------------------------------------

def test_json_is_readable_by_others_like_any_file(tmp_path):
    old = os.umask(0o022)
    try:
        path = jsonio.write_json(tmp_path / "x.json", {"a": 1})
    finally:
        os.umask(old)
    assert path.stat().st_mode & 0o044 == 0o044 or jsonio._FILE_MODE & 0o044 == 0


class _Short(io.BytesIO):
    headers = {"Content-Length": "100"}


def test_a_download_cut_short_is_refused_and_not_kept(tmp_path, monkeypatch):
    monkeypatch.setattr(reference, "_open_url", lambda url: _Short(b"x" * 40))
    with pytest.raises(ProjectError) as info:
        reference._download("http://example.org/f", tmp_path / "f")
    assert info.value.code == 10 and not (tmp_path / "f").exists() and not (tmp_path / "f.part").exists()


def test_a_corrupt_gzip_is_a_numbered_error(tmp_path):
    path = tmp_path / "gene2ensembl.gz"
    data = bytearray(gzip.compress(("#tax_id\tGeneID\tEnsembl_gene_identifier\n" +
                                    "".join(f"9606\t{i}\tENSG{i:011d}\n" for i in range(5000))).encode()))
    data[200:208] = b"\xff" * 8
    path.write_bytes(bytes(data))
    with pytest.raises(ProjectError) as info:
        list(reference._ncbi_rows(path, ("GeneID",)))
    assert info.value.code == 11


def test_threads_default_to_the_cpus_this_process_may_use(monkeypatch):
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {0}, raising=False)
    assert parallel.resolve_threads(None) == 1


# --- statistics and descriptions --------------------------------------------------------------------------------------

def test_a_regulator_without_signed_targets_is_judged_by_name():
    out = score({"collectri": [(0.01, "TF1", frozenset()), (0.02, "TF2", frozenset())],
                 "signalling_via_tf": [(0.01, "TF1", frozenset())]}, ["collectri", "signalling_via_tf"], 0.5)
    assert out["collectri"]["unique_findings"] == 1 and out["signalling_via_tf"]["unique_findings"] == 0
    assert out["collectri"]["independent_findings"] == 2


def test_a_gene_switched_off_is_the_largest_ratio_change():
    b = resolve("H", "c", "contrast", {"A": {"P-value": 0.01, "Fold_change": 0.0},
                                       "B": {"P-value": 0.01, "Fold_change": 1.1}}, fc_ratio=True)
    assert b.picks["Fold_change"] == "A"


def test_an_inhibited_regulator_lists_the_targets_that_went_against_it(monkeypatch):
    from evidome import report

    reg = {"comparison": "C", "network": "collectri", "regulator": "TF4", "regulator_kind": "tf",
           "targets_in_universe": "20", "activity_q": "0.001", "activity_t": "-9", "activity_p": "0.0001",
           "consistent_targets": "", "inconsistent_targets": "G1;G2", "overlap_q": "0.5", "overlap_p": "0.4"}
    monkeypatch.setattr(report, "_enrichment", lambda src, name: ([], iter(())))
    monkeypatch.setattr(report, "_regulator_rows", lambda src: ([], iter([reg])))
    src = argparse.Namespace(comparisons=[{"comparison": "C"}])
    rows = report.findings_rows(src)
    assert rows[0]["direction"] == "down" and rows[0]["genes"] == "G1;G2"


def test_signs_multiply_down_a_chain_of_controls(tmp_path):
    from evidome.pathway_signs import pathway_patterns, read_biopax
    from test_pathways_get_activity_patterns_from_reactome import _owl

    extra = "".join(
        f'<bp:Control rdf:ID="{cid}"><bp:controlled rdf:resource="#N1"/><bp:controller rdf:resource="#{pid}"/>'
        f'<bp:controlType>{kind}</bp:controlType></bp:Control>'
        f'<bp:Protein rdf:ID="{pid}"><bp:entityReference rdf:resource="#PR_{acc}"/></bp:Protein>'
        f'<bp:ProteinReference rdf:ID="PR_{acc}"><bp:xref rdf:resource="#X_{acc}"/></bp:ProteinReference>'
        f'<bp:UnificationXref rdf:ID="X_{acc}"><bp:db>UniProt</bp:db><bp:id>{acc}</bp:id></bp:UnificationXref>'
        for cid, pid, kind, acc in (("N3", "ACT", "ACTIVATION", "P00011"), ("N4", "INH", "INHIBITION", "P00012")))
    path = tmp_path / "biopax.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("Homo_sapiens.owl", _owl().replace("</rdf:RDF>", extra + "</rdf:RDF>"))
    patterns, _ = pathway_patterns(read_biopax(path))
    # N1 inhibits R1: an activator of N1 brakes the pathway, an inhibitor of N1 releases the brake
    assert patterns["R-HSA-100"]["P00011"] == -1 and patterns["R-HSA-100"]["P00012"] == 1


def test_reruns_write_the_same_zip_and_workbook(world, tmp_path):
    ref, table = world
    out = tmp_path / "study"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref)] + OPTS) == 0
    first = [(out / n).read_bytes() for n in ("study_combined_results.zip", "study_workbook.xlsx")]
    assert main(["report", str(out)]) == 0
    assert [(out / n).read_bytes() for n in ("study_combined_results.zip", "study_workbook.xlsx")] == first


# --- pass 4 ---------------------------------------------------------------------------------------------------------

from evidome.enrich import enrichment_score, fold_change_convention, null_es  # noqa: E402
from test_reference_and_crosswalk import M, _mirror_for  # noqa: E402


def test_a_new_species_is_found_and_fetched_by_the_crosswalk(tmp_path, monkeypatch):
    base, ref = _mirror_for(tmp_path)
    monkeypatch.setattr(reference, "DEFAULT_BASE_URL", base)
    asked = []
    monkeypatch.setattr(reference, "species_of_ids", lambda ids, **kw: asked.append(list(ids)) or {M: 2})
    table = tmp_path / "t.csv"
    table.write_text("Gene,Common_name,A_ANOVA_p,A_ANOVA_q,A_F\n"
                     "ENSMMUG00000000001,TSPAN6,0.01,0.1,5\nENSMMUG00000000002,,0.5,0.9,1\n")
    assert main(["crosswalk", str(table), "--out", str(tmp_path / "o"), "--reference-dir", str(ref)]) == 0
    assert asked and M in jsonio.read_json(ref / "ensembl-116" / "reference.json")["species"]
    asked.clear()          # the second run finds it downloaded and asks nobody
    assert main(["crosswalk", str(table), "--out", str(tmp_path / "o2"), "--reference-dir", str(ref)]) == 0
    assert not asked


def test_without_network_the_missing_species_error_stands(tmp_path, capsys):
    table = tmp_path / "t.csv"
    table.write_text("Gene,A_ANOVA_p,A_ANOVA_q,A_F\nENSMMUG00000000001,0.01,0.1,5\n")
    assert main(["crosswalk", str(table), "--out", str(tmp_path / "o"), "--reference-dir", str(tmp_path / "r")]) == 2
    err = capsys.readouterr().err
    assert "EVD-0013" in err and "could not ask Ensembl" in err


def test_species_with_mixed_ids_are_not_guessed(tmp_path, monkeypatch):
    from evidome.crosswalk import ensure_species

    monkeypatch.setattr(reference, "species_of_ids", lambda ids, **kw: {"mus_musculus": 5, "rattus_norvegicus": 5})
    with pytest.raises(ProjectError) as info:
        ensure_species([f"ENSMUSG{i:011d}" for i in range(10)], tmp_path / "r", log=lambda *_: None)
    assert info.value.code == 14 and "several species" in str(info.value)


def test_a_log2_fold_change_keeps_its_direction_and_a_one_is_up():
    assert fold_change_convention([-1.2, 0.3, 1.0, 0.0]) == "log"
    from evidome.enrich import ComparisonTable  # noqa: F401
    from evidome import enrich

    path_rows = [{"human_gene_id": f"G{i}", "human_symbol": "", "tested": "yes", "P-value": "0.01",
                  "Fold_change": v} for i, v in enumerate(["1.0", "-0.4", "0.3", "0.0"])]
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "c.csv")
        with open(p, "w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=list(path_rows[0]))
            w.writeheader()
            w.writerows(path_rows)
        t = enrich.read_comparison(p, "c", "contrast")
    assert [t.direction[f"G{i}"] for i in range(4)] == [1, -1, 1, 0]


def test_an_exact_gsea_tie_goes_to_the_positive_side():
    w = [1.0] * 39
    assert enrichment_score(w, [5, 20, 33])[0] > 0 and null_es(w, [5, 20, 33]) > 0


def test_combine_refuses_to_overwrite_an_input_and_skips_its_own_output(tmp_path):
    folder = _sc_folder(tmp_path)
    victim = folder / "Gq_vs_WT_0_naive_cell_t_rna.csv"
    with pytest.raises(ProjectError) as info:
        combine_mod.find_tables(folder, out=victim)
    assert info.value.code == 31
    combine_mod.combine(combine_mod.find_tables(folder), folder / "wide.csv", log=lambda *_: None)
    (folder / "._hidden.csv").write_text("junk")
    (folder / "ortholog_crosswalk.csv").write_text("junk")
    assert len(combine_mod.find_tables(folder)) == 2           # wide.csv, the hidden file, run outputs skipped
    with pytest.raises(ProjectError) as info:
        combine_mod.combine(combine_mod.find_tables(folder), tmp_path, log=lambda *_: None)
    assert info.value.code == 31


def test_combine_honours_id_column_headers_and_placeholder_symbols(tmp_path):
    folder = tmp_path / "f"
    folder.mkdir()
    (folder / "A.csv").write_text("# generated by paralome 0.5\n#   X: pairwise, control=WT, treatment=Gq\n"
                                  "ens,Common_name,Control_mean,Treatment_mean,Fold_change,P-value,Q-value\n"
                                  "ENSMUSG01,ENSMUSG01,1,2,1,0.01,0.05\n")
    (folder / "B.csv").write_text("ens,Common_name,Control_mean,Treatment_mean,Fold_change,P-value,Q-value\n"
                                  "ENSMUSG01,Real,1,2,1,0.01,0.05\n")
    out = combine_mod.combine(combine_mod.find_tables(folder), tmp_path / "w.csv", id_column="ens",
                              log=lambda *_: None)
    table = read_wide_statistics(out)
    a = {c.name: c for c in table.comparisons}["A"]
    assert (a.control, a.treatment) == ("WT", "Gq") and table.symbols == ["Real"]
    assert "paralome 0.5" in table.generator


def test_a_repeated_column_in_a_per_comparison_file_is_refused(tmp_path):
    folder = tmp_path / "f"
    folder.mkdir()
    (folder / "A.csv").write_text("Gene,Control_mean,Treatment_mean,Fold_change,P-value,Q-value,P-value\n"
                                  "ENSMUSG01,1,2,1,0.01,0.05,0.9\n")
    with pytest.raises(ProjectError) as info:
        combine_mod.combine(combine_mod.find_tables(folder), tmp_path / "w.csv", log=lambda *_: None)
    assert info.value.code == 8


def test_inspect_on_a_folder_leaves_nothing_behind(tmp_path, capsys):
    import glob
    import tempfile

    before = set(glob.glob(os.path.join(tempfile.gettempdir(), "evidome_*")))
    assert main(["inspect", str(_sc_folder(tmp_path))]) == 0
    assert set(glob.glob(os.path.join(tempfile.gettempdir(), "evidome_*"))) == before
    assert "(combined)" in capsys.readouterr().out


def test_rows_without_an_id_are_left_out_of_a_wide_table_too(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text(BODY.format(t="3") + ".,1,2,2,0.01,0.05,3\n.,1,2,2,0.01,0.05,3\n")
    table = read_wide_statistics(path)
    assert table.gene_ids == ["ENSMMUG1"] and any("2 row(s)" in n for n in table.notes)


def test_utility_options_are_range_checked(tmp_path):
    assert main(["utility", str(tmp_path), "--q", "5"]) == 2


def test_a_damaged_biopax_download_names_refresh(tmp_path):
    from evidome.pathway_signs import read_biopax

    bad = tmp_path / "biopax.zip"
    bad.write_text("<html>not a zip</html>")
    with pytest.raises(ProjectError) as info:
        read_biopax(bad)
    assert info.value.code == 18 and "--refresh" in str(info.value)


def test_a_few_known_ids_do_not_download_a_species(tmp_path, monkeypatch):
    from evidome.crosswalk import ensure_species

    monkeypatch.setattr(reference, "species_of_ids", lambda ids, **kw: {"macaca_mulatta": 1})
    fetched = []
    monkeypatch.setattr(reference, "fetch", lambda *a, **k: fetched.append(a))
    with pytest.raises(ProjectError) as info:
        ensure_species([f"X{i}" for i in range(60)], tmp_path / "r", log=lambda *_: None)
    assert info.value.code == 14 and not fetched


def test_the_rest_sample_reaches_the_end_of_the_table(monkeypatch):
    import json as _json
    sent = {}

    class _Resp(io.BytesIO):
        def __init__(self):
            super().__init__(b"{}")

    def fake_open(request, timeout=None):
        sent["ids"] = _json.loads(request.data)["ids"]
        return _Resp()

    monkeypatch.undo()          # the conftest guard replaces species_of_ids; this test needs the real one
    monkeypatch.setattr(reference.urllib.request, "urlopen", fake_open)
    reference.species_of_ids([f"G{i}" for i in range(149)])
    assert len(sent["ids"]) == 50 and sent["ids"][-1] in {f"G{i}" for i in range(140, 149)}


def test_two_genome_names_that_clean_up_alike_keep_separate_caches():
    from evidome.genomes import _safe

    assert _safe("sheep v2") != _safe("sheep_v2") and _safe("sheep_v2") == "sheep_v2"
