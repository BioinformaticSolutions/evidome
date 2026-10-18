"""'<name> TABLE --out DIR' runs every step; there is no 'run' command to type."""

import json

import pytest

from evidome.cli import COMMANDS, PIPELINE, _route, build_parser, main

from test_evidence_bundle_describes_itself import COLLECTRI, SIGNALLING, TABLE  # noqa: F401
from test_gene_sets_are_indexed_to_human_ensembl import sources  # noqa: F401  (fixture)


def test_a_table_is_routed_to_the_pipeline_and_commands_are_not():
    assert _route(["stats.csv", "--out", "x"])[0] == PIPELINE
    assert _route(["enrich", "x"]) == ["enrich", "x"]
    assert _route(["--version"]) == ["--version"] and _route([]) == []
    help_text = build_parser().format_help()
    entries = [line.split()[0] for line in help_text.splitlines()
               if line.startswith("    ") and not line.startswith("     ")]
    assert PIPELINE not in help_text and "run" not in entries and "TABLE --out DIR" in help_text
    assert "run-finish" not in help_text and "enrich-task" not in help_text


@pytest.fixture
def world(sources, tmp_path):
    from evidome.genesets import fetch_genesets
    from evidome.regulators import fetch_regulators

    ref, urls = sources
    fetch_genesets(ref, urls=urls, threads=1, log=lambda *_: None)
    (tmp_path / "c.tsv").write_text(COLLECTRI)
    (tmp_path / "s.tsv").write_text(SIGNALLING)
    fetch_regulators(ref, urls={"collectri": (tmp_path / "c.tsv").as_uri(), "omnipath": (tmp_path / "s.tsv").as_uri()},
                     threads=1, log=lambda *_: None)
    table = tmp_path / "stats.csv"
    table.write_text(TABLE)
    return ref, table


OPTS = ["--threads", "1", "--min-species-match", "0.5", "--min-size", "1", "--permutations", "100",
        "--max-permutations", "1000", "--min-targets", "2", "--figure-formats", "svg"]


def test_one_command_writes_everything_locally(world, tmp_path, capsys):
    ref, table = world
    out = tmp_path / "study"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref)] + OPTS) == 0
    for path in ("ortholog_crosswalk.csv", "enrichment/enrichment.json", "regulators/regulators.json",
                 "evidence/manifest.json", "study_combined_results.zip", "study_workbook.xlsx"):
        assert (out / path).exists(), path
    text = capsys.readouterr().out
    assert "gene sets with q < 0.05" in text and "regulators with q < 0.05" in text
    assert "figures:" in text or "figures in" in text     # this tiny fixture has nothing two-sided to draw


def test_the_same_files_as_step_by_step(world, tmp_path):
    ref, table = world
    one, steps = tmp_path / "one", tmp_path / "steps"
    assert main([str(table), "--out", str(one), "--reference-dir", str(ref)] + OPTS) == 0
    common = ["--reference-dir", str(ref), "--threads", "1"]
    assert main(["crosswalk", str(table), "--out", str(steps), "--min-species-match", "0.5"] + common) == 0
    assert main(["enrich", str(steps), "--min-size", "1", "--permutations", "100", "--max-permutations", "1000"]
                + common) == 0
    assert main(["regulators", str(steps), "--min-targets", "2"] + common) == 0
    for name in ("enrichment/A/gsea.csv", "enrichment/F1/ora.csv", "regulators/A/regulators.csv"):
        assert (one / name).read_bytes() == (steps / name).read_bytes(), name


def test_slurm_chains_finish_after_finalize(world, tmp_path, monkeypatch, capsys):
    from evidome import slurm

    calls = []

    def fake(script, dependency=None, runner=None):
        calls.append((script.name, dependency))
        return str(500 + len(calls))

    monkeypatch.setattr(slurm, "submit", fake)
    ref, table = world
    out = tmp_path / "study"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref), "--backend", "slurm",
                 "--cpus-per-task", "8"] + OPTS) == 0
    assert calls == [("task.sbatch", None), ("finalize.sbatch", "afterok:501"), ("finish.sbatch", "afterok:502")]
    printed = [line for line in capsys.readouterr().out.splitlines() if line.startswith("{")]
    assert json.loads(printed[-1]) == {"submitted_job_ids": ["501", "502", "503"]}   # for Paralome's supervisor
    finish = (out / "enrichment/_slurm/finish.sbatch").read_text()
    assert "run-finish --plan" in finish and "#SBATCH --cpus-per-task=2" in finish       # 2 comparisons
    plan = json.loads((out / "enrichment/_slurm/plan.json").read_text())
    assert plan["regulators"]["min_targets"] == 2
    assert plan["figures"]["formats"] == ["svg"] and plan["figures"]["colors"] == "#0044CC,#CC0000"

    # what the three jobs will run, run here in order
    from evidome.enrich import run_slurm_finalize, run_slurm_task
    from evidome.pipeline import run_finish

    plan_path = out / "enrichment/_slurm/plan.json"
    for i in range(2):
        run_slurm_task(plan_path, i, threads=1, log=lambda *_: None)
    run_slurm_finalize(plan_path)
    lines = run_finish(plan_path, threads=1, log=lambda *_: None)
    assert (out / "study_combined_results.zip").exists() and any("workbook" in l for l in lines)
    assert any(l.startswith(("figures", "wrote")) and "figure" in l for l in lines)


def test_every_listed_command_is_known_to_the_router():
    parser = build_parser()
    listed = {a.dest for a in parser._subparsers._group_actions[0]._choices_actions}  # noqa: SLF001
    assert listed <= set(COMMANDS)
