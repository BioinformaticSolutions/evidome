"""Upstream-regulator evidence: signed regulons from OmniPath, overlap p, activation z, activity t."""

import csv
import json

import pytest

from evidome.enrich import read_comparison
from evidome.errors import ProjectError
from evidome.regulators import (analyse_comparison, binom_two_sided_p, build_regulons, fetch_regulators,
                                  load_regulons, omnipath_urls)

from test_enrichment_is_exact_and_reproducible import _gene, _world

HEADER = ("source\ttarget\tsource_genesymbol\ttarget_genesymbol\tis_directed\tis_stimulation\tis_inhibition\t"
          "consensus_direction\tconsensus_stimulation\tconsensus_inhibition\tsources\tcuration_effort\n")


def _row(src, tgt, s_sym, t_sym, stim, inhib):
    return {"source": src, "target": tgt, "source_genesymbol": s_sym, "target_genesymbol": t_sym,
            "consensus_stimulation": str(stim), "consensus_inhibition": str(inhib), "sources": "CollecTRI"}


SYMBOLS = {"A": "EA", "B": "EB", "C": "EC", "TF1": "ETF1", "TF2": "ETF2", "KIN": "EKIN"}


def test_signs_complexes_conflicts_and_unmapped_targets():
    collectri = [_row("P1", "x", "TF1", "A", True, False), _row("P1", "x", "TF1", "B", False, True),
                 _row("P1", "x", "TF1", "C", True, True),                              # both: unsigned
                 _row("P1", "x", "TF1", "A", False, True),                             # conflicts with +: 0
                 _row("COMPLEX:P1_P2", "x", "TF1_TF2", "B", True, False),
                 _row("P1", "x", "TF1", "NOPE", True, False), _row("P1", "x", "TF1", "DUP", True, False)]
    rows, counts = build_regulons(collectri, [], SYMBOLS, {"DUP"})
    got = {(r[1], r[4]): r[6] for r in rows}
    assert got[("TF1", "EA")] == 0 and got[("TF1", "EB")] == -1 and got[("TF1", "EC")] == 0
    assert [r[2] for r in rows if r[1] == "TF1_TF2"] == ["complex"] and rows[0][3] in ("", "ETF1")
    assert counts["collectri_targets_unmapped"] == 1 and counts["collectri_targets_ambiguous"] == 1


def test_signalling_regulators_multiply_signs_through_their_tfs_and_cancel():
    collectri = [_row("P1", "x", "TF1", "A", True, False), _row("P1", "x", "TF1", "B", True, False),
                 _row("P2", "x", "TF2", "A", True, False), _row("P2", "x", "TF2", "C", False, True)]
    omnipath = [_row("K", "P1", "KIN", "TF1", True, False),          # KIN activates TF1
                _row("K", "P2", "KIN", "TF2", False, True),          # KIN inhibits TF2
                _row("P1", "P1", "TF1", "TF1", True, False)]         # self-loop: ignored
    rows, counts = build_regulons(collectri, omnipath, SYMBOLS, set())
    via = {r[4]: (r[6], r[7]) for r in rows if r[0] == "signalling_via_tf"}
    # A: +1*+1 + (-1)*(+1) = 0 -> dropped; B: +1; C: (-1)*(-1) = +1
    assert set(via) == {"EB", "EC"} and via["EB"][0] == 1 and via["EC"][0] == 1
    assert via["EB"][1] == "via TF1;TF2" and counts["signalling_regulators"] == 1


def _regs(**defs):
    return {("collectri", name): {"kind": "tf", "gene": gene, "targets": targets}
            for name, (gene, targets) in defs.items()}


def test_planted_activation_and_repression_are_read_the_right_way(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    up, down = [_gene(i) for i in range(20)], [_gene(i) for i in range(365, 380)]
    regs = _regs(ACT=(_gene(0), {**{g: 1 for g in up}, **{g: -1 for g in down}}),
                 REP=("", {g: -1 for g in up}),
                 NOISE=("", {_gene(i): (1 if i % 2 else -1) for i in range(100, 160)}),
                 TINY=("", {_gene(0): 1}))
    rows = {r["regulator"]: r for r in analyse_comparison(t, regs, alpha=0.05, min_targets=5, max_targets=1000)}
    act, rep, noise = rows["ACT"], rows["REP"], rows["NOISE"]
    assert "TINY" not in rows
    assert act["n_consistent"] == 35 and act["n_inconsistent"] == 0 and act["activation_z"] > 5
    assert act["activity_t"] > 5 and act["activity_q"] < 1e-6 and act["overlap_q"] < 1e-6
    assert rep["activation_z"] < -4 and rep["activity_t"] < -4 and rep["inconsistent_targets"].startswith("S0;")
    assert noise["activity_p"] > 0.01
    assert act["regulator_in_data"] == "yes" and act["regulator_direction"] == "up"
    assert rep["regulator_in_data"] == "no"


def test_ftests_get_overlap_only(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "F1.human_genes.csv", "F1", "ftest")
    regs = _regs(ACT=("", {_gene(i): 1 for i in range(20)}))
    (row,) = analyse_comparison(t, regs, alpha=0.05, min_targets=5, max_targets=1000)
    assert row["overlap_q"] < 1e-6 and "activity_t" not in row and "activation_z" not in row


def test_indirect_hubs_above_max_targets_are_not_tested(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    regs = {("signalling_via_tf", "HUB"): {"kind": "signalling", "gene": "", "targets": {_gene(i): 1 for i in range(300)}},
            ("collectri", "BIG"): {"kind": "tf", "gene": "", "targets": {_gene(i): 1 for i in range(300)}}}
    rows = analyse_comparison(t, regs, alpha=0.05, min_targets=5, max_targets=100)
    assert [r["regulator"] for r in rows] == ["BIG"]


def test_binomial_is_symmetric_and_exact():
    assert binom_two_sided_p(5, 10) == pytest.approx(1.0)
    assert binom_two_sided_p(0, 10) == pytest.approx(2 / 1024)
    assert binom_two_sided_p(10, 10) == binom_two_sided_p(0, 10)


def test_urls_ask_for_commercial_records_only():
    assert all("license=commercial" in u and "organisms=9606" in u for u in omnipath_urls().values())


def test_no_index_is_prv_0025(tmp_path):
    with pytest.raises(ProjectError) as info:
        load_regulons(tmp_path)
    assert info.value.code == 25


def test_a_changed_omnipath_header_is_prv_0018(tmp_path):
    from test_reference_and_crosswalk import M, _mirror_for
    from evidome.reference import fetch

    base, ref = _mirror_for(tmp_path)
    fetch(ref, [M], base_url=base, log=lambda *_: None)
    bad = tmp_path / "bad.tsv"
    bad.write_text("from\tto\nA\tB\n")
    with pytest.raises(ProjectError) as info:
        fetch_regulators(ref, urls={"collectri": bad.as_uri(), "omnipath": bad.as_uri()}, threads=1,
                         log=lambda *_: None)
    assert info.value.code == 18


# ---------------------------------------------------------------------------
# bias-corrected z, via / driving TFs
# ---------------------------------------------------------------------------

def test_bias_corrected_z_matches_kraemer_by_hand(tmp_path):
    import math

    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    sig = [g for g in t.genes if t.p[g] < 0.05]
    up = sum(1 for g in sig if t.direction[g] > 0)
    down = sum(1 for g in sig if t.direction[g] < 0)
    targets = {**{_gene(i): 1 for i in range(20)}, **{_gene(i): 1 for i in range(370, 375)},
               **{_gene(i): -1 for i in range(200, 210)}}
    (row,) = analyse_comparison(t, _regs(R=("", targets)), alpha=0.05, min_targets=5, max_targets=1000)
    mu = (up - down) / (up + down) * (25 - 10) / 35
    n, c, i = row["n_consistent"] + row["n_inconsistent"], row["n_consistent"], row["n_inconsistent"]
    assert row["bias_term"] == pytest.approx(mu, abs=1e-6)
    assert row["activation_z_bias_corrected"] == pytest.approx((c - i - n * mu) / math.sqrt(n * (1 - mu * mu)),
                                                              abs=1e-6)
    assert row["activation_z"] == pytest.approx((c - i) / math.sqrt(n), abs=1e-6)


def test_signalling_regulators_name_the_tfs_driving_them(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    regs = {("collectri", "TFA"): {"kind": "tf", "gene": "", "via": [], "targets": {_gene(i): 1 for i in range(12)}},
            ("collectri", "TFB"): {"kind": "tf", "gene": "", "via": [], "targets": {_gene(i): 1 for i in range(30, 34)}},
            ("collectri", "TFC"): {"kind": "tf", "gene": "", "via": [], "targets": {_gene(i): 1 for i in range(100, 105)}},
            ("signalling_via_tf", "KIN"): {"kind": "signalling", "gene": "", "via": ["TFA", "TFB", "TFC"],
                                           "targets": {_gene(i): 1 for i in list(range(12)) + list(range(30, 34))}}}
    rows = {r["regulator"]: r for r in analyse_comparison(t, regs, alpha=0.05, min_targets=3, max_targets=1000)}
    assert rows["KIN"]["via_tfs"] == "TFA;TFB;TFC"
    assert rows["KIN"]["driving_tfs"].startswith("TFA(12);TFB(4)")
    assert "via_tfs" not in rows["TFA"]


def test_a_signalling_only_run_still_knows_the_tfs(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    regs = {("collectri", "TFA"): {"kind": "tf", "gene": "", "via": [], "targets": {_gene(i): 1 for i in range(12)}},
            ("signalling_via_tf", "KIN"): {"kind": "signalling", "gene": "", "via": ["TFA"],
                                           "targets": {_gene(i): 1 for i in range(12)}}}
    rows = analyse_comparison(t, regs, alpha=0.05, min_targets=3, max_targets=1000, networks={"signalling_via_tf"})
    assert [r["regulator"] for r in rows] == ["KIN"] and rows[0]["driving_tfs"] == "TFA(12)"
