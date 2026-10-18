"""A burden test's evidence: no direction-of-effect statistics, and an ORA null weighted by detectability.

Under a burden test's null a gene reaches p < 0.05 with a probability set by how many qualifying alleles the
reference shows for it -- orders of magnitude apart between genes -- so a plain hypergeometric ORA finds gene
sets rich in "detectable" genes enriched by construction. The producer gives that probability per gene
(Null_selection_probability) and ORA uses Wallenius' noncentral hypergeometric with it, as goseq does. Proved
here: the Wallenius tail against the exact draw-by-draw urn, and the ORA's calibration under a simulated
weighted null (Monte Carlo), where the plain hypergeometric is not calibrated.
"""

import csv
import json
import math
import random

import pytest

from evidome import figures
from evidome.cli import main
from evidome.enrich import (NULL_HYPERGEOMETRIC, NULL_UNADJUSTED, NULL_WALLENIUS, ComparisonTable, Settings,
                            detectability_odds, hypergeom_sf, ora_null, overlap_p, wallenius_pmf, wallenius_sf)
from evidome.errors import ProjectError
from evidome.figures import FigureSettings, draw, panels_for
from evidome.regulators import REGULATOR_COLUMNS
from evidome.wide_stats import NULL_SELECTION_ALPHA, NULL_SELECTION_ROLE, read_wide_statistics

from test_a_burden_table_is_read_as_a_burden_test import COMPARISONS, HEADER, ROWS, _burden_table
from test_evidence_bundle_describes_itself import COLLECTRI, SIGNALLING  # noqa: F401
from test_gene_sets_are_indexed_to_human_ensembl import sources  # noqa: F401  (fixture)


# ---------------------------------------------------------------------------
# Wallenius, proved against the urn it describes
# ---------------------------------------------------------------------------

def _urn(m1, m2, n, omega):
    """Exact Wallenius distribution: n balls drawn one at a time, each weighted (m1 balls omega, m2 balls 1)."""
    dist = {0: 1.0}
    for drawn in range(n):
        nxt = {}
        for red, prob in dist.items():
            left_red, left_white = m1 - red, m2 - (drawn - red)
            total = omega * left_red + left_white
            if left_red:
                nxt[red + 1] = nxt.get(red + 1, 0.0) + prob * omega * left_red / total
            if left_white:
                nxt[red] = nxt.get(red, 0.0) + prob * left_white / total
        dist = nxt
    return dist


@pytest.mark.parametrize("m1, m2, n, omega", [
    (5, 7, 6, 3.0), (10, 50, 8, 0.2), (30, 300, 40, 5.0), (12, 1000, 25, 0.05), (3, 3, 3, 10.0),
    (20, 5000, 10, 100.0), (200, 18000, 400, 1.7), (15, 18000, 150, 0.3), (100, 100, 150, 0.01),
])
def test_the_wallenius_tail_matches_the_exact_urn(m1, m2, n, omega):
    exact = _urn(m1, m2, n, omega)
    for k in range(0, min(n, m1) + 2):
        tail = sum(v for x, v in exact.items() if x >= k)
        got = wallenius_sf(k, m1 + m2, m1, n, omega)
        if tail > 1e-250:
            assert got == pytest.approx(tail, rel=1e-8, abs=0.0), (k, got, tail)
        else:
            assert got < 1e-240
    for x, v in exact.items():
        if v > 1e-250:
            assert wallenius_pmf(x, n, m1, m2, omega) == pytest.approx(v, rel=1e-8)


def test_odds_of_one_is_the_hypergeometric_exactly():
    for k, N, K, n in ((3, 18000, 100, 300), (0, 50, 10, 5), (7, 60, 12, 20), (12, 40, 12, 30)):
        assert wallenius_sf(k, N, K, n, 1.0) == hypergeom_sf(k, N, K, n)
        # and the integral route (odds a hair from 1) agrees with it
        assert wallenius_sf(k, N, K, n, 1.0 + 1e-12) == pytest.approx(hypergeom_sf(k, N, K, n), rel=1e-6)


def test_the_odds_are_the_mean_selection_probability_inside_over_outside():
    weight = {"a": 0.2, "b": 0.4, "c": 0.01, "d": 0.03}
    assert detectability_odds(["a", "b"], list(weight), weight) == pytest.approx(0.3 / 0.02)
    assert detectability_odds(["a", "b"], list(weight), weight, total=0.64) == pytest.approx(0.3 / 0.02)
    assert detectability_odds(list(weight), list(weight), weight) == 1.0          # nothing outside


# ---------------------------------------------------------------------------
# Calibration under the weighted null (Monte Carlo)
# ---------------------------------------------------------------------------

def _universe(rng, n_genes=1500):
    """Selection probabilities spread ~450x (as on a real burden reference), log-uniform in [1e-4, 0.045]."""
    genes = [f"ENSG{i:011d}" for i in range(n_genes)]
    pi = {g: math.exp(rng.uniform(math.log(1e-4), math.log(0.045))) for g in genes}
    return genes, pi


def _rejection_rates(seed, biased, sims=1500):
    rng = random.Random(seed)
    genes, pi = _universe(rng)
    ranked = sorted(genes, key=lambda g: -pi[g])
    # a set of 80: the most detectable genes (biased), or a random draw (unbiased)
    members = ranked[:300:4][:80] if biased else rng.sample(genes, 80)
    table = ComparisonTable("b", "burden", genes, {g: g for g in genes}, {}, {}, {}, selection=pi)
    plain = weighted = 0
    for _ in range(sims):
        chosen = {g for g in genes if rng.random() < pi[g]}
        k = sum(1 for g in members if g in chosen)
        n = len(chosen)
        plain += overlap_p(table, NULL_HYPERGEOMETRIC, k, len(genes), members, n) < 0.05
        weighted += overlap_p(table, NULL_WALLENIUS, k, len(genes), members, n) < 0.05
    return plain / sims, weighted / sims


def test_under_the_weighted_null_the_plain_ora_is_biased_and_the_wallenius_ora_is_not():
    """Genes selected independently with their own probabilities (the burden null); a set of the most
    detectable genes. The hypergeometric calls it enriched in most simulations; Wallenius at about the nominal
    rate (it is an approximation: a band, not an equality)."""
    plain, weighted = _rejection_rates(2026, biased=True)
    assert plain > 0.25, plain                    # nominal 0.05; measured 0.35
    assert 0.005 <= weighted <= 0.08, weighted    # measured 0.028 (a discrete test sits a little under 0.05)


def test_an_unbiased_set_is_calibrated_either_way():
    plain, weighted = _rejection_rates(7, biased=False)
    assert plain <= 0.08 and weighted <= 0.08, (plain, weighted)


# ---------------------------------------------------------------------------
# The input contract
# ---------------------------------------------------------------------------

ROLE = NULL_SELECTION_ROLE
SEL_COLUMNS = ["gene_id", "gene_name"] + [f"{c}_{r}" for c in COMPARISONS
                                          for r in ("Log2_fold_change", "P-value", "Q-value", ROLE)]
SELECTION = (0.03, 0.001, 0.02, 0.004, 0.0002, 0.01, 0.04)


def _rows_with_selection(selection=SELECTION):
    out = []
    for row, w in zip(ROWS, selection):
        gid, sym, l1, p1, q1, l2, p2, q2 = row
        w2 = "." if p2 == "." else w
        out.append((gid, sym, l1, p1, q1, w, l2, p2, q2, w2))
    return out


def test_a_burden_test_may_carry_its_null_selection_probabilities(tmp_path):
    table = read_wide_statistics(_burden_table(tmp_path / "b.tsv", columns=SEL_COLUMNS, rows=_rows_with_selection()))
    for comp in table.comparisons:
        assert comp.kind == "burden" and ROLE in comp.columns
    assert table.column(table.comparisons[0], ROLE)[1] == 0.001
    assert NULL_SELECTION_ALPHA == Settings().alpha == 0.05


def test_the_role_is_refused_on_anything_but_a_burden_test(tmp_path):
    header = "".join(f"# {c}: control=gnomAD, treatment=probands\n" for c in COMPARISONS)
    with pytest.raises(ProjectError) as err:
        read_wide_statistics(_burden_table(tmp_path / "c.tsv", header=header, columns=SEL_COLUMNS,
                                           rows=_rows_with_selection()))
    assert err.value.code == 8 and "burden test only" in str(err.value)


@pytest.mark.parametrize("bad, words", [(".", "no dominant_ptv_Null_selection_probability"),
                                        (1.5, "is not a probability")])
def test_a_tested_gene_without_a_valid_probability_is_refused(tmp_path, bad, words):
    selection = (bad,) + SELECTION[1:]
    with pytest.raises(ProjectError) as err:
        read_wide_statistics(_burden_table(tmp_path / "d.tsv", columns=SEL_COLUMNS,
                                           rows=_rows_with_selection(selection)))
    assert err.value.code == 7 and words in str(err.value)


def test_which_null_a_burden_comparison_gets():
    genes = ["g1", "g2"]
    burden = ComparisonTable("b", "burden", genes, {}, {}, {}, {}, selection={"g1": 0.1, "g2": 0.2})
    assert ora_null(burden, 0.05) == NULL_WALLENIUS
    assert ora_null(burden, 0.01) == NULL_UNADJUSTED          # probabilities are for p < 0.05
    assert ora_null(ComparisonTable("b", "burden", genes, {}, {}, {}, {}), 0.05) == NULL_UNADJUSTED
    assert ora_null(ComparisonTable("c", "contrast", genes, {}, {}, {}, {}), 0.05) == NULL_HYPERGEOMETRIC


# ---------------------------------------------------------------------------
# End to end: the evidence for a burden table
# ---------------------------------------------------------------------------

OPTS = ["--threads", "1", "--min-species-match", "0.5", "--min-size", "1", "--permutations", "100",
        "--max-permutations", "1000", "--min-targets", "2", "--figure-formats", "svg"]


@pytest.fixture
def weighted_run(sources, tmp_path):
    from evidome.genesets import fetch_genesets
    from evidome.regulators import fetch_regulators

    ref, urls = sources
    fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None)
    (tmp_path / "c.tsv").write_text(COLLECTRI)
    (tmp_path / "s.tsv").write_text(SIGNALLING)
    fetch_regulators(ref, urls={"collectri": (tmp_path / "c.tsv").as_uri(), "omnipath": (tmp_path / "s.tsv").as_uri()},
                     threads=1, log=lambda *_: None)
    out = tmp_path / "burden_evidence"
    table = _burden_table(tmp_path / "gene_burden.tsv", columns=SEL_COLUMNS, rows=_rows_with_selection())
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref)] + OPTS) == 0
    return out


def test_the_burden_ora_is_weighted_and_says_so(weighted_run):
    out = weighted_run
    with (out / "comparisons/dominant_ptv.human_genes.csv").open() as h:
        rows = {r["human_gene_id"]: r for r in csv.DictReader(h)}
    assert rows["ENSG00000000002"][ROLE] == "0.001"
    with (out / "enrichment/dominant_ptv/ora.csv").open() as h:
        ora = list(csv.DictReader(h))
    assert ora and {r["null_model"] for r in ora} == {NULL_WALLENIUS}
    enrichment = json.loads((out / "enrichment/enrichment.json").read_text())
    entry = {c["comparison"]: c for c in enrichment["comparisons"]}["dominant_ptv"]
    assert entry["ora_null"] == NULL_WALLENIUS and "ora_null_note" not in entry
    assert entry["direction_statistics"].startswith("not reported: a burden test measures damage counts")
    assert entry["gsea"] is None and "confounded" in entry["gsea_skipped"]
    assert "Wallenius" in enrichment["settings"]["ora"]
    # no activation z on a burden comparison's ORA
    assert all(r["pathway_activation_z"] == "" and r["pattern_consistent"] == "" for r in ora)


def test_burden_regulators_report_overlap_only(weighted_run):
    with (weighted_run / "regulators/dominant_ptv/regulators.csv").open() as h:
        rows = list(csv.DictReader(h))
    assert rows
    for r in rows:
        assert r["overlap_p"] != "" and r["overlap_null"] == NULL_WALLENIUS
        for blank in ("activation_z", "z_p", "bias_term", "activation_z_bias_corrected", "activity_t",
                      "activity_p", "activity_q", "regulator_direction", "n_consistent", "n_inconsistent"):
            assert r[blank] == "", (blank, r[blank])
    manifest = json.loads((weighted_run / "regulators/regulators.json").read_text())
    entry = {c["comparison"]: c for c in manifest["comparisons"]}["dominant_ptv"]
    assert "damage counts, not expression direction" in entry["direction_statistics"]


def test_the_evidence_bundle_describes_the_burden_choices(weighted_run):
    evidence = weighted_run / "evidence"
    with (evidence / "comparisons.csv").open() as h:
        comps = list(csv.DictReader(h))
    assert {c["ora_null"] for c in comps} == {NULL_WALLENIUS}
    assert all("damage counts" in c["direction_statistics"] for c in comps)
    manifest = json.loads((evidence / "manifest.json").read_text())
    described = {c["name"] for f in manifest["files"] for c in f["columns"]}
    assert {"null_model", "overlap_null", ROLE, "ora_null", "direction_statistics"} <= described
    assert any("detectability" in line for line in manifest["reading_guide"])


def test_a_burden_figure_is_one_panel_ranked_by_p_without_activation_words(weighted_run):
    if not figures.available():
        pytest.skip("matplotlib not installed")
    svgs = sorted((weighted_run / "figures").glob("*.svg"))
    assert svgs
    for svg in svgs:
        text = svg.read_text()
        assert "burden test: probands vs gnomAD" in text
        for word in ("Predicted activated", "Predicted inhibited", "Activation z", "Increased", "Decreased"):
            assert word not in text, (svg.name, word)
        with svg.with_suffix(".csv").open() as h:
            data = list(csv.DictReader(h))
        assert {r["side"] for r in data} == {"ranked by p"}
        ps = [float(r["p"]) for r in data]
        assert ps == sorted(ps)


def test_a_burden_panel_holds_every_side_in_one_ranking(tmp_path):
    d = tmp_path / "study"
    path = d / "regulators/dominant_ptv/regulators.csv"
    path.parent.mkdir(parents=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, REGULATOR_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows([{"network": "collectri", "regulator": f"R{i}", "overlap_p": 10 ** -i, "activation_z": z}
                     for i, z in ((3, 2.5), (5, -3.0), (1, ""))])
    [panel] = panels_for(d, "dominant_ptv", 25, kind="burden")
    assert panel.single and panel.left == [] and panel.left_title == "" and panel.right_title == ""
    assert [i.label for i in panel.right] == ["R5", "R3", "R1"]
    assert all(i.score is None for i in panel.right)
    if figures.available():
        written = draw(panel, tmp_path / "figs", FigureSettings(formats=("svg",)), version="t")
        assert {p.suffix for p in written} == {".svg", ".csv"}
