"""Enrichment: exact ORA, GSEA as Subramanian 2005, reproducible whatever --threads is.

The implementation was validated once against scipy (hypergeometric, BH) and
GSEApy 1.3.1 (ES identical to 1e-15, NES within permutation noise); these
tests pin the arithmetic without either as a dependency.
"""

import csv
import json
import math
import random

import pytest

from evidome.cli import main
from evidome.enrich import (Settings, bh, enrich, enrichment_score, fold_change_convention, hypergeom_sf,
                              read_comparison)
from evidome.errors import ProjectError


def _brute_sf(k, N, K, n):
    total = math.comb(N, n)
    return sum(math.comb(K, x) * math.comb(N - K, n - x) for x in range(k, min(K, n) + 1)) / total


def test_hypergeometric_tail_is_exact_on_both_sides_of_the_mode():
    rng = random.Random(0)
    for _ in range(400):
        N = rng.randint(5, 300)
        K, n = rng.randint(0, N), rng.randint(0, N)
        k = rng.randint(0, min(K, n) + 1)
        assert hypergeom_sf(k, N, K, n) == pytest.approx(_brute_sf(k, N, K, n), rel=1e-9, abs=1e-300)


def test_bh_matches_the_textbook_example():
    assert [round(q, 4) for q in bh([0.01, 0.04, 0.03, 0.005])] == [0.02, 0.04, 0.04, 0.02]


def test_enrichment_score_by_hand():
    # 4 genes, weights 4 3 2 1, set = positions 0 and 3: +4/5, -1/2, -1/2 (to -0.2), +1/5 (to 0)
    es, peak = enrichment_score([4, 3, 2, 1], [0, 3])
    assert es == pytest.approx(0.8) and peak == 0
    es, peak = enrichment_score([4, 3, 2, 1], [2, 3])      # bottom of the list: -1 before the first hit
    assert es == pytest.approx(-1.0) and peak == -1


@pytest.mark.parametrize("values,expected", [([2.0, -3.0, 1.0], "signed-times"), ([2.0, 0.5, 1.0], "ratio")])
def test_fold_change_convention_is_read_from_the_column(values, expected):
    assert fold_change_convention(values) == expected


def test_a_log_fold_change_column_is_read_as_log():
    """Negative values with values between -1 and 1 are a log fold change (Paralome single-cell),
    read rather than refused as EVD-0020."""
    assert fold_change_convention([-2.0, 0.5, 0.0, 1.0]) == "log"


# ---------------------------------------------------------------------------
# end to end: a crosswalk output directory and a gene-set index, built by hand
# ---------------------------------------------------------------------------

N_GENES = 400
MAXP = 20000            # keeps the adaptive rounds quick; the default cap is 100,000


def _gene(i):
    return f"ENSG{i:011d}"


def _world(tmp_path, with_t=True, fc_convention="signed-times"):
    rng = random.Random(5)
    cw = tmp_path / "cw"
    (cw / "comparisons").mkdir(parents=True)
    # genes 0..39 strongly up, 360..399 strongly down, the rest noise; the last 10 untested
    rows = []
    for i in range(N_GENES):
        if i >= 390:
            rows.append({"human_gene_id": _gene(i), "human_symbol": f"S{i}", "tested": "no"})
            continue
        t = 6 - i * 0.01 if i < 40 else (-6 + (399 - i) * 0.01 if i >= 350 else rng.gauss(0, 1))
        p = min(1.0, 2 * math.exp(-abs(t) * 1.5))
        fc = 2 ** (t / 3)
        if fc_convention == "signed-times" and fc < 1:
            fc = -1 / fc
        row = {"human_gene_id": _gene(i), "human_symbol": f"S{i}", "tested": "yes", "P-value": repr(p),
               "Fold_change": repr(fc), "F": "", "ANOVA_p": repr(p)}
        if with_t:
            row["T-statistic"] = repr(t)
        rows.append(row)
    cols = ["human_gene_id", "human_symbol", "tested", "P-value", "Fold_change", "F", "ANOVA_p"] + \
           (["T-statistic"] if with_t else [])
    for name in ("C", "F1"):
        with (cw / "comparisons" / f"{name}.human_genes.csv").open("w", newline="") as h:
            w = csv.DictWriter(h, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
    (cw / "comparison_summary.json").write_text(json.dumps({"comparisons": [
        {"comparison": "C", "kind": "contrast"}, {"comparison": "F1", "kind": "ftest"}]}))

    ref = tmp_path / "ref"
    gs = ref / "ensembl-116" / "genesets"
    gs.mkdir(parents=True)
    sets = {
        "TOP": list(range(0, 30)),
        "BOTTOM": list(range(360, 390)),
        "NOISE": list(range(100, 160)),
        "TINY": [0, 1, 2],                                  # below --min-size
        "UNTESTED_HEAVY": list(range(385, 400)),            # only 5 tested members: below min size
    }
    with (gs / "gene_sets.tsv").open("w") as h:
        h.write("set_id\tsource\tname\tcategory\tn_genes\tsource_ids\n")
        for sid, genes in sets.items():
            src = "hallmark" if sid != "NOISE" else "reactome"
            h.write(f"X:{sid}\t{src}\t{sid.title()}\t\t{len(genes)}\t{len(genes)}\n")
    with (gs / "members.tsv").open("w") as h:
        h.write("set_id\thuman_gene_id\n")
        for sid, genes in sets.items():
            for g in genes:
                h.write(f"X:{sid}\t{_gene(g)}\n")
    (gs / "genesets.json").write_text(json.dumps({"indexed_by": "test", "files": {}, "licences": {}}))
    return cw, ref


def _rows(path):
    return {(r.get("direction", ""), r["set_id"]): r for r in csv.DictReader(path.open())}


def test_ora_finds_the_planted_sets_in_the_right_direction(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(max_permutations=MAXP, permutations=200), threads=1, log=lambda *_: None)
    ora = _rows(cw / "enrichment" / "C" / "ora.csv")
    assert float(ora[("up", "X:TOP")]["q"]) < 1e-10 and float(ora[("down", "X:TOP")]["p"]) == 1.0
    assert float(ora[("down", "X:BOTTOM")]["q"]) < 1e-10
    assert float(ora[("all", "X:NOISE")]["p"]) > 0.05
    assert ora[("up", "X:TOP")]["universe"] == "390"             # untested genes are not in the universe
    assert ("up", "X:TINY") not in ora and ("up", "X:UNTESTED_HEAVY") not in ora
    assert ora[("up", "X:TOP")]["overlap_genes"].startswith("S0;S1;")


def test_gsea_ranks_by_t_and_finds_both_ends(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(max_permutations=MAXP, permutations=500), threads=1, log=lambda *_: None)
    g = _rows(cw / "enrichment" / "C" / "gsea.csv")
    top, bottom, noise = g[("", "X:TOP")], g[("", "X:BOTTOM")], g[("", "X:NOISE")]
    assert float(top["nes"]) > 2 and float(top["p"]) < 0.01      # at the floor: 1 / (same-sign nulls + 1)
    assert float(bottom["nes"]) < -2 and float(bottom["es"]) < 0
    assert float(noise["p"]) > 0.05
    assert top["rank_metric"] == "T-statistic" and int(top["leading_edge_size"]) == 30


def test_without_t_the_rank_is_signed_log_p_and_direction_comes_from_fold_change(tmp_path):
    cw, ref = _world(tmp_path, with_t=False)
    enrich(cw, ref, settings=Settings(max_permutations=MAXP, permutations=200), threads=1, log=lambda *_: None)
    manifest = json.loads((cw / "enrichment" / "enrichment.json").read_text())
    c = manifest["comparisons"][0]
    assert c["direction_from"] == "Fold_change (signed-times)" and c["rank_metric"].startswith("sign(")
    assert float(_rows(cw / "enrichment" / "C" / "ora.csv")[("down", "X:BOTTOM")]["q"]) < 1e-10


def test_a_ratio_fold_change_gives_the_same_directions(tmp_path):
    cw, ref = _world(tmp_path, with_t=False, fc_convention="ratio")
    enrich(cw, ref, settings=Settings(max_permutations=MAXP, permutations=200), threads=1, log=lambda *_: None)
    manifest = json.loads((cw / "enrichment" / "enrichment.json").read_text())
    assert manifest["comparisons"][0]["direction_from"] == "Fold_change (ratio)"
    assert float(_rows(cw / "enrichment" / "C" / "ora.csv")[("down", "X:BOTTOM")]["q"]) < 1e-10


def test_ftests_get_ora_all_only_and_no_gsea(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(max_permutations=MAXP, permutations=100), threads=1, log=lambda *_: None)
    assert {d for d, _ in _rows(cw / "enrichment" / "F1" / "ora.csv")} == {"all"}
    assert not (cw / "enrichment" / "F1" / "gsea.csv").exists()
    manifest = json.loads((cw / "enrichment" / "enrichment.json").read_text())
    assert manifest["comparisons"][1]["gsea_skipped"].startswith("F-test")
    assert manifest["comparisons"][0]["sets_outside_size_range"] == {"hallmark": 2}


def test_results_do_not_depend_on_threads(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, out=tmp_path / "one", settings=Settings(max_permutations=MAXP, permutations=200), threads=1, log=lambda *_: None)
    enrich(cw, ref, out=tmp_path / "many", settings=Settings(max_permutations=MAXP, permutations=200), threads=6, log=lambda *_: None)
    for name in ("C/ora.csv", "C/gsea.csv", "F1/ora.csv"):
        assert (tmp_path / "one" / name).read_bytes() == (tmp_path / "many" / name).read_bytes()


def test_a_different_seed_changes_only_the_permutation_results(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, out=tmp_path / "a", settings=Settings(max_permutations=MAXP, permutations=200, seed=1), threads=1, log=lambda *_: None)
    enrich(cw, ref, out=tmp_path / "b", settings=Settings(max_permutations=MAXP, permutations=200, seed=2), threads=1, log=lambda *_: None)
    assert (tmp_path / "a/C/ora.csv").read_bytes() == (tmp_path / "b/C/ora.csv").read_bytes()
    a, b = _rows(tmp_path / "a/C/gsea.csv"), _rows(tmp_path / "b/C/gsea.csv")
    assert a[("", "X:NOISE")]["es"] == b[("", "X:NOISE")]["es"]


def test_cli_prints_the_counts_and_refuses_a_non_crosswalk_dir(tmp_path, capsys):
    cw, ref = _world(tmp_path)
    assert main(["enrich", str(cw), "--reference-dir", str(ref), "--permutations", "100", "--threads", "1"]) == 0
    out = capsys.readouterr().out
    assert "ORA   hallmark/all" in out and "GSEA  skipped: F-test" in out
    assert main(["enrich", str(tmp_path), "--reference-dir", str(ref)]) == 2
    assert "EVD-0021" in capsys.readouterr().err


def test_read_comparison_keeps_only_tested_genes(tmp_path):
    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    assert len(t.genes) == 390 and t.direction[_gene(0)] == 1 and t.direction[_gene(380)] == -1


# ---------------------------------------------------------------------------
# adaptive permutations
# ---------------------------------------------------------------------------

def test_only_sizes_with_extreme_sets_get_more_permutations(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(permutations=200, max_permutations=MAXP), threads=1, log=lambda *_: None)
    g = _rows(cw / "enrichment" / "C" / "gsea.csv")
    top, noise = g[("", "X:TOP")], g[("", "X:NOISE")]
    assert int(top["permutations"]) == MAXP                    # 200 -> 2,000 -> 20,000
    assert float(top["p"]) < 1 / 5000                          # far below the 200-permutation floor
    assert int(noise["permutations"]) == 200                   # never extended: plenty of nulls beat it


def test_capping_at_the_first_round_switches_adaptivity_off(tmp_path):
    cw, ref = _world(tmp_path)
    enrich(cw, ref, settings=Settings(permutations=200, max_permutations=200), threads=1, log=lambda *_: None)
    assert {r["permutations"] for r in _rows(cw / "enrichment" / "C" / "gsea.csv").values()} == {"200"}


def test_adaptive_results_do_not_depend_on_threads(tmp_path):
    cw, ref = _world(tmp_path)
    for name, threads in (("one", 1), ("many", 5)):
        enrich(cw, ref, out=tmp_path / name, settings=Settings(permutations=300, max_permutations=30000),
               threads=threads, log=lambda *_: None)
    assert (tmp_path / "one/C/gsea.csv").read_bytes() == (tmp_path / "many/C/gsea.csv").read_bytes()


def test_a_null_job_yields_exactly_its_permutations(tmp_path):
    from evidome.enrich import Null, _null_task, gsea_inputs

    cw, _ = _world(tmp_path)
    t = read_comparison(cw / "comparisons" / "C.human_genes.csv", "C", "contrast")
    weights = gsea_inputs(t)[1]
    (size, pos, neg), = _null_task(("C", weights, [(60, 100, "1:C:60")]))
    null = Null()
    null.add(pos, neg, 100)
    assert null.permutations == 100 and len(null.pos) + len(null.neg) == 100


# ---------------------------------------------------------------------------
# the faster null kernel
# ---------------------------------------------------------------------------

def test_the_null_kernel_matches_the_reference_enrichment_score():
    from evidome.enrich import null_es, random_positions

    rng = random.Random(4)
    weights = sorted((rng.gauss(0, 2) for _ in range(3000)), reverse=True)
    abs_w = [abs(w) for w in weights]
    for size in (1, 2, 15, 150, 2999):
        for _ in range(50):
            hits = random_positions(rng, len(weights), size)
            assert null_es(abs_w, hits) == pytest.approx(enrichment_score(weights, hits)[0], abs=1e-12)


def test_random_positions_are_distinct_sorted_and_uniform():
    from collections import Counter

    from evidome.enrich import random_positions

    rng = random.Random(9)
    counts = Counter()
    for _ in range(20000):
        pos = random_positions(rng, 50, 10)
        assert pos == sorted(set(pos)) and len(pos) == 10 and 0 <= pos[0] and pos[-1] < 50
        counts.update(pos)
    assert 3700 < min(counts.values()) and max(counts.values()) < 4300          # expected 4000 each
    assert random_positions(random.Random("a"), 1000, 30) == random_positions(random.Random("a"), 1000, 30)
    assert random_positions(rng, 12, 12) == list(range(12))                     # k == n still terminates
