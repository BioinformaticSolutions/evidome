"""Per-comparison human-gene tables: crosswalk joined to statistics,
many-to-one collisions decided by the majority vote per comparison, nothing dropped."""

import csv
import json

from evidome.human_tables import build_human_tables, write_human_tables
from evidome.wide_stats import read_wide_statistics

TABLE = (
    "Gene,Common_name,"
    "A_Raw_FPKM_Control_mean,A_Raw_FPKM_Treatment_mean,A_Fold_change,A_Effect,A_T-statistic,A_P-value,A_Q-value,"
    "B_F,B_ANOVA_p,B_ANOVA_q\n"
    # M2 beats M3 in A on every metric; M3 beats M2 in B on both.
    "M1,one,1,2,2.0,1.0,3.0,0.01,0.05,9,0.001,0.01\n"
    "M2,two,1,4,4.0,2.0,5.0,0.001,0.01,1,0.6,0.9\n"
    "M3,three,1,2,1.5,0.5,1.0,0.3,0.5,20,0.0001,0.001\n"
    "M4,four,1,3,3.0,1.5,4.0,0.02,0.1,5,0.03,0.2\n"
    "M6,six,1,5,5.0,2.5,6.0,0.004,0.02,7,0.01,0.1\n"
    "M7,seven,.,.,.,.,.,.,.,.,.,.\n"
    "M8,eight,.,.,.,.,.,.,.,.,.,.\n"
)


def _cw(gid, sym, biotype, hid="", hsym="", status="mapped", targets=1, sharing=1, htype="ortholog_one2one",
        source="protein"):
    return {"input_gene_id": gid, "input_symbol": sym, "input_biotype": biotype, "species": "x",
            "human_gene_id": hid, "human_symbol": hsym, "human_biotype": "protein_coding" if hid else "",
            "homology_type": htype if hid else "", "is_high_confidence": "1" if hid else "",
            "identity_pct": "90" if hid else "", "human_identity_pct": "89" if hid else "",
            "goc_score": "", "wga_coverage": "", "homology_source": source if hid else "",
            "human_targets_of_input_gene": targets if hid else 0,
            "input_genes_sharing_human_gene": sharing if hid else 0, "mapping_status": status}


CROSSWALK = [
    _cw("M1", "one", "protein_coding", "H1", "HONE"),
    _cw("M2", "two", "protein_coding", "H2", "HTWO", sharing=2, htype="ortholog_many2many"),
    _cw("M3", "three", "protein_coding", "H2", "HTWO", sharing=2, htype="ortholog_many2many"),
    _cw("M4", "four", "protein_coding", "H4", "HFOUR", targets=2, htype="ortholog_one2many"),
    _cw("M4", "four", "protein_coding", "H5", "HFIVE", targets=2, htype="ortholog_one2many"),
    _cw("M6", "six", "protein_coding", status="no human ortholog in reference"),
    _cw("M7", "seven", "lncRNA", status="no human ortholog in reference"),
    _cw("M8", "eight", "protein_coding", "H8", "HEIGHT"),
]


def _run(tmp_path, threads=1):
    path = tmp_path / "t.csv"
    path.write_text(TABLE)
    table = read_wide_statistics(path)
    return table, build_human_tables(table, CROSSWALK, threads=threads)


def _by_human(result):
    return {r["human_gene_id"]: r for r in result.human_rows}


def test_every_human_gene_reached_gets_one_row_per_comparison(tmp_path):
    _, results = _run(tmp_path)
    assert [r.name for r in results] == ["A", "B"]
    for result in results:
        assert sorted(_by_human(result)) == ["H1", "H2", "H4", "H5", "H8"]


def test_the_collision_winner_is_decided_per_comparison(tmp_path):
    _, (a, b) = _run(tmp_path)
    assert _by_human(a)["H2"]["source_gene_id"] == "M2" and _by_human(a)["H2"]["decided_by"] == "vote"
    assert _by_human(b)["H2"]["source_gene_id"] == "M3"
    assert _by_human(a)["H2"]["P-value"] == "0.001" and _by_human(b)["H2"]["ANOVA_p"] == "0.0001"


def test_the_full_ballot_is_written_one_row_per_candidate(tmp_path):
    _, (a, _) = _run(tmp_path)
    rows = {r["candidate_gene_id"]: r for r in a.ballot_rows}
    assert set(rows) == {"M2", "M3"}
    assert rows["M2"]["votes"] == 4 and rows["M2"]["winner"] == "yes"
    assert rows["M2"]["metrics_won"] == "P-value;T-statistic;Fold_change;Effect"
    assert rows["M3"]["votes"] == 0 and rows["M3"]["winner"] == "no"


def test_one_to_many_feeds_every_target_and_says_so(tmp_path):
    _, (a, _) = _run(tmp_path)
    h = _by_human(a)
    assert h["H4"]["source_gene_id"] == h["H5"]["source_gene_id"] == "M4"
    assert h["H4"]["source_human_targets"] == 2 and h["H4"]["homology_type"] == "ortholog_one2many"


def test_untested_human_gene_is_kept_and_labelled(tmp_path):
    _, (a, _) = _run(tmp_path)
    h8 = _by_human(a)["H8"]
    assert h8["tested"] == "no" and h8["source_gene_id"] == "" and h8["P-value"] == ""
    assert h8["decided_by"] == "none tested"


def test_summary_counts_where_the_significant_genes_went(tmp_path):
    _, (a, b) = _run(tmp_path)
    s = a.summary
    assert s["input_genes_tested"] == 5 and s["input_genes_significant"] == 4   # M1 M2 M4 M6
    assert s["significant_mapped"] == 3 and s["significant_unmapped"] == 1
    assert s["significant_unmapped_by_biotype"] == {"protein_coding": 1}
    assert s["significant_unmapped_by_status"] == {"no human ortholog in reference": 1}
    assert s["significant_mapped_but_outvoted"] == 0
    assert s["human_genes_significant"] == 4                                     # H1 H2 H4 H5
    # In B, M2 (p 0.6) is not significant; M3 wins H2. Nobody significant is outvoted.
    assert b.summary["input_genes_significant"] == 4 and b.summary["significant_mapped_but_outvoted"] == 0


def test_a_significant_loser_is_counted_as_outvoted(tmp_path):
    path = tmp_path / "t.csv"
    path.write_text(TABLE.replace("M3,three,1,2,1.5,0.5,1.0,0.3,", "M3,three,1,2,1.5,0.5,1.0,0.04,"))
    table = read_wide_statistics(path)
    a = build_human_tables(table, CROSSWALK, threads=1)[0]
    assert a.summary["significant_mapped_but_outvoted"] == 1


def test_parallel_and_sequential_give_identical_tables(tmp_path):
    _, one = _run(tmp_path, threads=1)
    _, many = _run(tmp_path, threads=4)
    assert [(r.human_rows, r.ballot_rows, r.summary) for r in one] == \
           [(r.human_rows, r.ballot_rows, r.summary) for r in many]


def test_files_are_written_per_comparison_with_their_own_statistic_columns(tmp_path):
    table, results = _run(tmp_path)
    write_human_tables(results, table, tmp_path / "out")
    a = list(csv.DictReader((tmp_path / "out/comparisons/A.human_genes.csv").open()))
    b_header = (tmp_path / "out/comparisons/B.human_genes.csv").open().readline().strip().split(",")
    assert len(a) == 5 and "T-statistic" in a[0] and "Q-value" in a[0]
    assert b_header[-3:] == ["F", "ANOVA_p", "ANOVA_q"]
    assert (tmp_path / "out/collisions/A.ballots.csv").exists()


def test_crosswalk_command_writes_the_tables_and_summary(tmp_path, capsys):
    """End to end through the CLI, against the mirrored reference used by the crosswalk tests."""
    from test_reference_and_crosswalk import M, _mirror_for

    from evidome.cli import main

    base, ref = _mirror_for(tmp_path)
    assert main(["reference", "fetch", "--species", M, "--reference-dir", str(ref), "--base-url", base]) == 0
    table = tmp_path / "t.csv"
    table.write_text("Gene,Common_name,A_ANOVA_p,A_ANOVA_q,A_F\n"
                     "ENSMMUG00000000001,TSPAN6,0.01,0.1,5\n"
                     "ENSMMUG00000000002,,0.5,0.9,1\n"
                     "ENSMMUG00000000003,,0.02,0.2,4\n"
                     "ENSMMUG00000000008,,0.001,0.01,9\n")
    out = tmp_path / "out"
    assert main(["crosswalk", str(table), "--out", str(out), "--reference-dir", str(ref), "--threads", "1"]) == 0
    rows = {r["human_symbol"]: r for r in csv.DictReader((out / "comparisons/A.human_genes.csv").open())}
    assert rows["PARA1"]["source_gene_id"] == "ENSMMUG00000000003"      # p 0.02 and F 4 beat p 0.5 and F 1
    summary = json.loads((out / "comparison_summary.json").read_text())
    a = summary["comparisons"][0]
    assert a["significant_unmapped_by_biotype"] == {"lncRNA": 1} and summary["program_version"]
    assert a["significant_mapped_by_source"] == {"protein": 2}
    assert rows["PARA1"]["homology_source"] == "protein"
    assert "significant but unmapped -- lncRNA 1" in capsys.readouterr().out


def test_a_compara_candidate_takes_precedence_over_ncbi_and_symbol(tmp_path):
    """M3 reached H2 only through NCBI. In B it has the better statistics, but M2 is Compara's."""
    path = tmp_path / "t.csv"
    path.write_text(TABLE)
    table = read_wide_statistics(path)
    crosswalk = [dict(r) for r in CROSSWALK]
    for r in crosswalk:
        if r["input_gene_id"] == "M3":
            r.update(homology_source="ncbi", homology_type="ortholog (NCBI)")
    a, b = build_human_tables(table, crosswalk, threads=1)
    h2 = _by_human(b)["H2"]
    assert h2["source_gene_id"] == "M2" and h2["decided_by"] == "Compara precedence"
    assert h2["n_tested_candidates"] == 2
    ballots = {r["candidate_gene_id"]: r for r in b.ballot_rows}
    assert ballots["M3"]["excluded"] == "Compara candidate present" and ballots["M3"]["votes"] == ""
    assert ballots["M2"]["excluded"] == "" and ballots["M2"]["winner"] == "yes"
    # M3 is significant in B (p 0.0001) and was set aside: counted as excluded, not as outvoted
    assert b.summary["significant_excluded_by_compara_precedence"] == 1
    assert b.summary["significant_mapped_but_outvoted"] == 0


def test_weaker_sources_still_compete_among_themselves(tmp_path):
    path = tmp_path / "t.csv"
    path.write_text(TABLE)
    table = read_wide_statistics(path)
    crosswalk = [dict(r) for r in CROSSWALK]
    for r in crosswalk:
        if r["input_gene_id"] in ("M2", "M3"):
            r["homology_source"] = "ncbi" if r["input_gene_id"] == "M2" else "symbol"
    _, b = build_human_tables(table, crosswalk, threads=1)
    assert _by_human(b)["H2"]["source_gene_id"] == "M3" and _by_human(b)["H2"]["decided_by"] == "vote"
    assert all(r["excluded"] == "" for r in b.ballot_rows)
