"""Defects that were found and fixed, each pinned so it cannot come back."""

import csv
import gzip

import pytest

from evidome import enrich as enrich_mod
from evidome.enrich import Settings, enrich, enrichment_score, null_es
from evidome.errors import ProjectError
from evidome.reference import index_gtf
from evidome.wide_stats import DIRECTION_ROLE, read_wide_statistics
from evidome.xlsx import _number
from test_enrichment_is_exact_and_reproducible import _world


# --- F-test directions in a three-factor design --------------------------------------------------------------

FTESTS = ["diet_WD_vs_CD", "drug_MET_vs_VEH", "sex_M_vs_F", "diet_x_drug", "diet_x_sex", "diet_x_drug_x_sex"]
CONTRASTS = {"WD_vs_CD_main_effect": 1.0, "MET_vs_VEH_main_effect": 2.0, "M_vs_F_main_effect": 3.0,
             "WD_vs_CD_within_VEH_contrast": 0.5, "WD_vs_CD_within_MET_contrast": 2.5,
             "M_vs_F_within_CD_contrast": -1.0, "M_vs_F_within_WD_contrast": 4.0}


def _three_factor(tmp_path):
    cols = ["Gene"]
    for f in FTESTS:
        cols += [f"{f}_F", f"{f}_ANOVA_p", f"{f}_ANOVA_q"]
    for c in CONTRASTS:
        cols += [f"{c}_Raw_FPKM_Control_mean", f"{c}_Raw_FPKM_Treatment_mean", f"{c}_Fold_change", f"{c}_Effect",
                 f"{c}_P-value", f"{c}_Q-value"]
    row = ["ENSMMUG0"] + ["5", "0.01", "0.05"] * len(FTESTS)
    for e in CONTRASTS.values():
        row += ["1", "2", "2", str(e), "0.01", "0.05"]
    path = tmp_path / "stats.csv"
    path.write_text(",".join(cols) + "\n" + ",".join(row) + "\n")
    return {c.name: c for c in read_wide_statistics(path).comparisons}, read_wide_statistics(path)


def test_each_interaction_takes_the_factor_pair_it_names(tmp_path):
    by, table = _three_factor(tmp_path)
    by = {c.name: c for c in table.comparisons}
    assert by["diet_x_drug"].direction_from == ("WD_vs_CD_within_MET_contrast (Effect) minus "
                                                "WD_vs_CD_within_VEH_contrast (Effect)")
    assert table.column(by["diet_x_drug"], DIRECTION_ROLE) == [2.0]
    assert by["diet_x_sex"].direction_from == ("M_vs_F_within_WD_contrast (Effect) minus "
                                               "M_vs_F_within_CD_contrast (Effect)")
    assert table.column(by["diet_x_sex"], DIRECTION_ROLE) == [5.0]


def test_a_three_way_interaction_gets_no_two_way_direction(tmp_path):
    _, table = _three_factor(tmp_path)
    three = {c.name: c for c in table.comparisons}["diet_x_drug_x_sex"]
    assert three.direction_from is None and DIRECTION_ROLE not in three.columns


# --- input tables ----------------------------------------------------------------------------------------------

BODY = ("Gene,X_Raw_FPKM_Control_mean,X_Raw_FPKM_Treatment_mean,X_Fold_change,X_P-value,X_Q-value\n"
        "ENSMMUG1,1,2,2,0.01,0.05\n")


def test_a_byte_order_mark_does_not_hide_the_id_column(tmp_path):
    path = tmp_path / "bom.csv"
    path.write_bytes(b"\xef\xbb\xbf" + BODY.encode())
    table = read_wide_statistics(path)
    assert table.id_column == "Gene" and table.gene_ids == ["ENSMMUG1"]


def test_a_file_that_is_not_utf8_is_a_numbered_error(tmp_path):
    path = tmp_path / "latin1.csv"
    path.write_bytes(BODY.encode() + b"ENSMMUG2,1,2,2,0.01,0.05\xe9\n")
    with pytest.raises(ProjectError) as info:
        read_wide_statistics(path)
    assert info.value.code == 7 and "UTF-8" in str(info.value)


# --- GENCODE annotations -----------------------------------------------------------------------------------------

def test_a_gencode_gtf_is_indexed_without_versions_and_with_its_biotypes(tmp_path):
    gtf = tmp_path / "g.gtf.gz"
    with gzip.open(gtf, "wt") as h:
        h.write('chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSMUSG00000000001.4"; gene_type "protein_coding"; '
                'gene_name "Gnai3";\n')
    out = tmp_path / "genes.tsv"
    assert index_gtf(gtf, out, ("gene_id", "symbol", "biotype")) == {"ENSMUSG00000000001"}
    rows = list(csv.reader(out.open(), delimiter="\t"))
    assert rows[1] == ["ENSMUSG00000000001", "Gnai3", "protein_coding"]


# --- GSEA ------------------------------------------------------------------------------------------------------

def test_observed_and_null_scores_break_ties_the_same_way():
    for w, hits in (([1.0] * 4, [0, 3]), ([3, 2, 1, 1, 2, 3], [0, 5])):
        assert enrichment_score(w, hits)[0] == null_es([abs(x) for x in w], hits) > 0


def test_every_adaptive_round_draws_new_permutations(tmp_path, monkeypatch):
    """Rounds 2 and 3 get different block seeds, even when round 2 draws fewer than 10,000."""
    tags = []
    real = enrich_mod._run_null_jobs

    def spy(jobs, *args, **kwargs):
        tags.extend((name, size, tag) for name, size, _, tag in jobs)
        return real(jobs, *args, **kwargs)

    monkeypatch.setattr(enrich_mod, "_run_null_jobs", spy)
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(permutations=100, max_permutations=100000, min_exceedances=10_000_000),
           threads=1, log=lambda *_: None)
    assert len(tags) == len(set(tags))


# --- the workbook ------------------------------------------------------------------------------------------------

def test_an_underscored_value_is_written_as_text_not_a_number():
    assert _number("1_2") is None and _number("12") == "12"


# --- pass 2 -------------------------------------------------------------------------------------------------------

import argparse  # noqa: E402

from evidome import cli  # noqa: E402
from evidome.checkpoint import reuse_decision  # noqa: E402
from evidome.collisions import resolve  # noqa: E402
from evidome.enrich import ComparisonTable, pathway_activation  # noqa: E402
from evidome.reference import fetch  # noqa: E402
from test_reference_and_crosswalk import M, _mirror_for  # noqa: E402


def test_pathway_activation_counts_only_the_sets_own_members():
    t = ComparisonTable("C", "contrast", ["A", "B", "X"], {"A": "A", "B": "B", "X": "X"},
                        {"A": 0.01, "B": 0.01, "X": 0.01}, {"A": 1, "B": 1, "X": 1}, {})
    t.direction_from = "Fold_change"
    out = pathway_activation(t, {"A": 1, "B": 1, "X": 1}, {"A", "B", "X"}, members={"A", "B"})
    assert out["pattern_genes"] == 2 and "X" not in str(out.get("pattern_consistent_genes", ""))


def test_a_misspelt_gene_set_source_is_refused():
    with pytest.raises(SystemExit):
        cli.build_parser().parse_args(["enrich", "x", "--sources", "Hallmark"])


def test_the_whole_run_names_a_missing_index_before_running(tmp_path):
    with pytest.raises(ProjectError) as info:
        cli._require_indexes(tmp_path, 116)
    assert info.value.code == 19


def test_options_before_the_table_go_to_the_whole_run():
    assert cli._route(["--threads", "2", "t.csv", "--out", "x"])[0] == cli.PIPELINE
    assert cli._route(["--version"]) == ["--version"]


def test_a_ratio_fold_change_is_judged_on_a_log_scale():
    b = resolve("H", "c", "contrast", {"A": {"P-value": 0.01, "Fold_change": 0.25},
                                       "B": {"P-value": 0.02, "Fold_change": 1.1}}, fc_ratio=True)
    assert b.winner == "A" and b.picks["Fold_change"] == "A" and not b.tie


def test_a_contrast_and_an_f_test_may_share_a_name(tmp_path):
    cols = ["Gene", "MET_vs_VEH_Raw_FPKM_Control_mean", "MET_vs_VEH_Raw_FPKM_Treatment_mean", "MET_vs_VEH_Fold_change",
            "MET_vs_VEH_P-value", "MET_vs_VEH_Q-value", "MET_vs_VEH_F", "MET_vs_VEH_ANOVA_p", "MET_vs_VEH_ANOVA_q"]
    path = tmp_path / "s.csv"
    path.write_text(",".join(cols) + "\nENSMMUG1,1,2,2,0.01,0.05,5,0.01,0.05\n")
    names = sorted((c.name, c.kind) for c in read_wide_statistics(path).comparisons)
    assert names == [("MET_vs_VEH", "contrast"), ("MET_vs_VEH_ANOVA", "ftest")]


def test_a_repeated_column_is_refused(tmp_path):
    path = tmp_path / "s.csv"
    path.write_text("Gene,A_Raw_FPKM_Control_mean,A_Raw_FPKM_Treatment_mean,A_Fold_change,A_P-value,A_Q-value,A_P-value\n"
                    "ENSMMUG1,1,2,2,0.01,0.05,0.9\n")
    with pytest.raises(ProjectError) as info:
        read_wide_statistics(path)
    assert info.value.code == 8


def test_a_checkpoint_that_is_not_an_object_is_not_reused(tmp_path):
    path = tmp_path / "checkpoint.json"
    path.write_text("[]")
    assert not reuse_decision(path, {})


def test_out_of_range_numbers_are_a_numbered_error():
    args = argparse.Namespace(alpha=0.05, min_size=10, max_size=500, permutations=0, seed=1,
                              max_permutations=0, min_exceedances=10)
    with pytest.raises(ProjectError) as info:
        cli._settings(args)
    assert info.value.code == 29


def test_database_value_states_its_own_cluster_defaults():
    args = cli.build_parser().parse_args(["database-value", "x"])
    assert (args.cpus_per_task, args.mem) == (16, "32G") and not hasattr(args, "figure_top")
    assert not hasattr(cli.build_parser().parse_args(["reference", "pathway-signs"]), "threads")


@pytest.mark.parametrize("lost", ["human_genes.tsv", f"species/{M}/ncbi_orthologs_to_human.tsv",
                                  f"species/{M}/species_genes.tsv", f"species/{M}/genes.txt"])
def test_a_deleted_index_file_is_rebuilt(tmp_path, lost):
    base, ref = _mirror_for(tmp_path)
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    (ref / "ensembl-116" / lost).unlink()
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    assert (ref / "ensembl-116" / lost).exists()
