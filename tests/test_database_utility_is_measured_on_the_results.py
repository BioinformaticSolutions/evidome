"""Each database's utility: what it finds that the others do not."""

import csv
import json

from evidome.cli import main
from evidome.enrich import ORA_COLUMNS
from evidome.regulators import REGULATOR_COLUMNS
from evidome.utility import independent, jaccard, source_utility


def _write(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def _ora(source, set_id, genes, q):
    return {"direction": "all", "source": source, "set_id": set_id, "name": set_id, "q": q,
            "overlap_genes": ";".join(genes)}


def _run(tmp_path):
    d = tmp_path / "run"
    d.mkdir()
    (d / "comparison_summary.json").write_text(json.dumps({"comparisons": [{"comparison": "A"}]}))
    abc, abd, xyz = list("ABCDEFGH"), list("ABCDEFGI"), list("XYZW")
    _write(d / "enrichment/A/ora.csv", ORA_COLUMNS, [
        _ora("reactome", "R1", abc, 0.001),
        _ora("go_bp", "GO1", abd, 0.001),            # the same finding as R1, said by GO
        _ora("go_bp", "GO2", abc, 0.002),            # ... and again by GO: one independent GO finding
        _ora("hallmark", "H1", xyz, 0.01),           # new
        _ora("wikipathways", "W1", abc, 0.5),        # not significant
        {**_ora("reactome", "R1", abc, 0.0001), "direction": "up"},   # only "all" counts
    ])
    _write(d / "regulators/A/regulators.csv", REGULATOR_COLUMNS, [
        {"network": "collectri", "regulator": "TP53", "overlap_q": 0.001, "consistent_targets": "A;B;C",
         "inconsistent_targets": "D"},
        {"network": "signalling_via_tf", "regulator": "ATM", "overlap_q": 0.001, "consistent_targets": "A;B;C",
         "inconsistent_targets": ""},
        {"network": "signalling_via_tf", "regulator": "EGFR", "overlap_q": 0.2, "consistent_targets": "Q"},
    ])
    return d


def test_jaccard_and_independence():
    a, b, c = frozenset("ABCD"), frozenset("ABCE"), frozenset("XY")
    assert jaccard(a, b) == 0.6 and jaccard(a, c) == 0
    assert independent([(0.01, "1", a), (0.02, "2", b), (0.03, "3", c)], 0.5) == 2


def test_utility_scores_what_each_database_adds(tmp_path):
    rows = {(r["family"], r["source"], r["comparison"]): r for r in source_utility(_run(tmp_path))}
    go = rows[("pathways", "go_bp", "ALL")]
    assert (go["findings"], go["independent_findings"], go["unique_findings"]) == (2, 1, 0)
    assert go["utility"] == 0.0 and go["unique_genes"] == 1 and go["marginal_gain"] == round(1 / 13, 4)  # gene I
    hall = rows[("pathways", "hallmark", "ALL")]
    assert (hall["findings"], hall["unique_findings"], hall["utility"]) == (1, 1, 1.0)
    assert hall["marginal_gain"] == round(4 / 13, 4)
    assert rows[("pathways", "wikipathways", "ALL")]["findings"] == 0          # tested, found nothing
    assert ("pathways", "go_mf", "ALL") not in rows                            # not tested: no row
    sig = rows[("regulators", "signalling_via_tf", "ALL")]
    assert sig["findings"] == 1 and sig["unique_findings"] == 0 and sig["unique_genes"] == 0


def test_the_command_prints_and_writes(tmp_path, capsys):
    d = _run(tmp_path)
    assert main(["utility", str(d)]) == 0
    text = capsys.readouterr().out
    assert "database utility" in text and "hallmark" in text
    written = list(csv.DictReader((d / "source_utility.csv").open()))
    assert {r["comparison"] for r in written} == {"ALL", "A"}
