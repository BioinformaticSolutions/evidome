"""The SLURM backend: one array task per comparison, then finalize -- same files as a local run."""

import json
import subprocess

import pytest

from evidome.enrich import Settings, enrich, plan_slurm, run_slurm_finalize, run_slurm_task
from evidome.errors import ProjectError
from evidome.parallel import resolve_threads
from evidome.slurm import SlurmOptions, script_text, submit

from test_enrichment_is_exact_and_reproducible import MAXP, _world

SETTINGS = Settings(permutations=200, max_permutations=MAXP)


class FakeSbatch:
    def __init__(self):
        self.calls = []

    def __call__(self, script, dependency=None):
        self.calls.append((script, dependency))
        return str(1000 + len(self.calls))


def _submitted(tmp_path, out_name="slurm_out"):
    cw, ref = _world(tmp_path)
    fake = FakeSbatch()
    ids = plan_slurm(cw, ref, options=SlurmOptions(cpus_per_task=8), out=tmp_path / out_name, settings=SETTINGS,
                     submit=fake, log=lambda *_: None)
    return cw, ref, fake, ids


def test_an_array_of_one_task_per_comparison_then_finalize_after_it(tmp_path):
    _, _, fake, ids = _submitted(tmp_path)
    (task, dep1), (final, dep2) = fake.calls
    assert dep1 is None and dep2 == "afterok:1001" and ids == {**ids, "array": "1001", "finalize": "1002"}
    text = task.read_text()
    assert "#SBATCH --array=0-1" in text and "#SBATCH --cpus-per-task=8" in text
    assert 'enrich-task --plan' in text and '"$SLURM_ARRAY_TASK_ID"' in text
    assert "#SBATCH --cpus-per-task=1" in final.read_text() and "enrich-finalize" in final.read_text()


def test_slurm_tasks_write_exactly_what_a_local_run_writes(tmp_path):
    cw, ref, _, ids = _submitted(tmp_path)
    for index in (1, 0):                                          # any order, as SLURM may run them
        run_slurm_task(ids["plan"], index, threads=3, log=lambda *_: None)
    manifest = run_slurm_finalize(ids["plan"])
    enrich(cw, ref, out=tmp_path / "local", settings=SETTINGS, threads=1, log=lambda *_: None)
    for name in ("C/ora.csv", "C/gsea.csv", "F1/ora.csv"):
        assert (tmp_path / "slurm_out" / name).read_bytes() == (tmp_path / "local" / name).read_bytes()
    a = json.loads(manifest.read_text())["comparisons"]
    b = json.loads((tmp_path / "local" / "enrichment.json").read_text())["comparisons"]
    assert a == b


def test_finalize_names_the_comparison_that_did_not_finish(tmp_path):
    _, _, _, ids = _submitted(tmp_path)
    run_slurm_task(ids["plan"], 0, threads=1, log=lambda *_: None)
    with pytest.raises(ProjectError) as info:
        run_slurm_finalize(ids["plan"])
    assert info.value.code == 22 and "F1" in str(info.value)


def test_a_node_running_another_version_refuses_the_task(tmp_path):
    """A compute node running an older install than the head node is reported, loudly."""
    _, _, _, ids = _submitted(tmp_path)
    plan = json.loads(open(ids["plan"]).read())
    plan["planned_by"] = "0.0.old"
    open(ids["plan"], "w").write(json.dumps(plan))
    with pytest.raises(ProjectError) as info:
        run_slurm_task(ids["plan"], 0, log=lambda *_: None)
    assert info.value.code == 22 and "0.0.old" in str(info.value)


def test_paths_with_spaces_are_quoted_in_sbatch_lines(tmp_path):
    text = script_text("job", SlurmOptions(partition="big mem"), tmp_path / "my logs" / "x.log", "echo hi")
    assert "#SBATCH --output='" in text and "#SBATCH --partition='big mem'" in text


def test_submit_reads_the_parsable_job_id_and_reports_refusals():
    ok = lambda args, **kw: subprocess.CompletedProcess(args, 0, "4242;cluster\n", "")
    assert submit("s.sbatch", runner=ok) == "4242"
    seen = []
    submit("s.sbatch", dependency="afterok:1", runner=lambda args, **kw: seen.append(args) or ok(args))
    assert seen[0][:3] == ["sbatch", "--parsable", "--dependency=afterok:1"]
    bad = lambda args, **kw: subprocess.CompletedProcess(args, 1, "", "invalid partition")
    with pytest.raises(ProjectError) as info:
        submit("s.sbatch", runner=bad)
    assert info.value.code == 23 and "invalid partition" in str(info.value)


def test_threads_follow_the_slurm_allocation(monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "40")
    assert resolve_threads(0) == 40 and resolve_threads(6) == 6
    monkeypatch.delenv("SLURM_CPUS_PER_TASK")
    assert resolve_threads(0) >= 1


def test_cli_without_sbatch_is_prv_0023_not_a_traceback(tmp_path, monkeypatch, capsys):
    from evidome.cli import main

    cw, ref = _world(tmp_path)
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    assert main(["enrich", str(cw), "--reference-dir", str(ref), "--backend", "slurm"]) == 2
    err = capsys.readouterr().err
    assert "EVD-0023" in err and "Traceback" not in err
