"""Goodness of fit and each database's Shapley share."""

import csv
import json
import random

import pytest

pytest.importorskip("sklearn")

from evidome import database_value as dvm   # noqa: E402
from evidome.cli import main                 # noqa: E402


def _world(tmp_path, n=600, seed=3):
    rng = random.Random(seed)
    genes = [f"ENSG{i:011d}" for i in range(n)]
    sets = {"REACTOME:R1": ("reactome", genes[:40]), "HALLMARK_H1": ("hallmark", genes[40:80])}
    for k in range(12):                                   # GO: random sets, no signal
        sets[f"GO:{k}"] = ("go_bp", rng.sample(genes, 30))
    sets["REACTOME:R1_PARENT"] = ("reactome", genes[:40] + rng.sample(genes[80:], 5))  # the child, restated
    targets = genes[100:130]                               # a regulator: 15 up, 15 down targets
    y = []
    for i, g in enumerate(genes):
        v = rng.gauss(0, 1)
        v += 2.0 if i < 40 else 0
        v -= 1.5 if 40 <= i < 80 else 0
        v += 1.5 if 100 <= i < 115 else (-1.5 if 115 <= i < 130 else 0)
        y.append(v)
    run = tmp_path / "run"
    (run / "comparisons").mkdir(parents=True)
    (run / "comparison_summary.json").write_text(json.dumps({"comparisons": [{"comparison": "A", "kind": "contrast"}]}))
    with (run / "comparisons/A.human_genes.csv").open("w", newline="") as h:
        w = csv.writer(h)
        w.writerow(["human_gene_id", "human_symbol", "tested", "P-value", "T-statistic"])
        for g, v in zip(genes, y):
            w.writerow([g, g[-4:], "yes", 0.5, v])
    ref = tmp_path / "ref"
    gs = ref / "ensembl-116/genesets"
    gs.mkdir(parents=True)
    with (gs / "gene_sets.tsv").open("w") as h:
        h.write("set_id\tsource\tname\tcategory\tn_genes\tsource_ids\n")
        for sid, (src, mem) in sets.items():
            h.write(f"{sid}\t{src}\t{sid}\t\t{len(mem)}\t{len(mem)}\n")
    with (gs / "members.tsv").open("w") as h:
        h.write("set_id\thuman_gene_id\n")
        for sid, (_, mem) in sets.items():
            h.writelines(f"{sid}\t{g}\n" for g in mem)
    rd = ref / "ensembl-116/regulators"
    rd.mkdir(parents=True)
    with (rd / "regulons.tsv").open("w") as h:
        h.write("network\tregulator\tregulator_kind\tregulator_human_gene_id\ttarget_human_gene_id\ttarget_symbol\t"
                "sign\tevidence\n")
        for i, g in enumerate(targets):
            h.write(f"collectri\tTF1\ttf\t\t{g}\t\t{1 if i < 15 else -1}\tcurated\n")
        for g in rng.sample(genes, 20):
            h.write(f"collectri\tTF2\ttf\t\t{g}\t\t1\tcurated\n")
    return run, ref


SMALL = dvm.ValueSettings(folds=3, screen=20, bootstrap=6, permutations=4, alphas=8)


def test_shapley_is_exact_and_adds_up():
    v = {(): 0.0, ("a",): 0.2, ("b",): 0.1, ("a", "b"): 0.25}
    phi = dvm.shapley(v, ["a", "b"])
    assert phi == pytest.approx({"a": 0.175, "b": 0.075}) and sum(phi.values()) == pytest.approx(0.25)


def test_the_databases_that_carry_signal_get_the_credit(tmp_path):
    run, ref = _world(tmp_path)
    lines = dvm.database_value(run, ref, settings=SMALL, threads=2, log=lambda *_: None)
    out = run / "database_value"
    share = {r["database"]: float(r["shapley_r2"]) for r in csv.DictReader((out / "database_value.csv").open())}
    fit = next(csv.DictReader((out / "goodness_of_fit.csv").open()))
    total = float(fit["heldout_r2_all"])
    assert total > 0.3 and sum(share.values()) == pytest.approx(total, abs=1e-4)
    assert min(share["reactome"], share["hallmark"], share["collectri"]) > 5 * max(share["go_bp"], 0.01)
    assert float(fit["null_95th_r2"]) < 0.05 < total and float(fit["null_p"]) == pytest.approx(0.2)
    curve = list(csv.DictReader((out / "curves.csv").open()))
    assert [int(r["databases"]) for r in curve] == [1, 2, 3, 4]
    assert float(curve[-1]["mean_r2"]) == pytest.approx(total, abs=1e-6)
    stable = {r["feature"]: float(r["selection_frequency"]) for r in csv.DictReader((out / "A/stable_features.csv").open())}
    assert stable.get("REACTOME:R1", 0) + stable.get("REACTOME:R1_PARENT", 0) >= 1.0 and stable["collectri:TF1"] == 1.0
    assert any("goodness of fit" in line for line in lines)


def test_the_command_runs_locally_and_draws_the_curve(tmp_path, capsys):
    run, ref = _world(tmp_path)
    assert main(["database-value", str(run), "--reference-dir", str(ref), "--folds", "3", "--screen", "20",
                 "--bootstrap", "4", "--permutations", "3", "--threads", "2", "--figure-formats", "svg"]) == 0
    assert "goodness of fit" in capsys.readouterr().out
    assert (run / "figures/A_database_value.svg").exists()


def test_without_scikit_learn_it_says_so(tmp_path, monkeypatch, capsys):
    run, ref = _world(tmp_path)
    monkeypatch.setattr(dvm, "available", lambda: False)
    assert main(["database-value", str(run), "--reference-dir", str(ref)]) == 2
    assert "EVD-0028" in capsys.readouterr().err
