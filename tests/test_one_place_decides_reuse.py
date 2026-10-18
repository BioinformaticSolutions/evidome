"""Checkpoint reuse is decided by version, then parameters, then outputs.

Each case below is a way Paralome handed back a stale result:
  - a new version with unchanged parameters (dev42 -> dev43 -> dev44),
  - a changed parameter with the same version (Patch 100),
  - a deleted output beside a checkpoint still saying complete (0.0.134).
"""

from evidome import checkpoint


def _complete(tmp_path, params=None, version=None):
    out = tmp_path / "result.tsv"
    out.write_text("x\n")
    cp = tmp_path / "stage.json"
    checkpoint.write_checkpoint(cp, params or {"alpha": 0.05}, {"table": str(out)})
    return cp, out


def test_matching_checkpoint_is_reused(tmp_path):
    cp, _ = _complete(tmp_path)
    decision = checkpoint.reuse_decision(cp, {"alpha": 0.05})
    assert decision.reuse, decision.reason


def test_missing_checkpoint_recomputes(tmp_path):
    decision = checkpoint.reuse_decision(tmp_path / "none.json", {})
    assert not decision and "no checkpoint" in decision.reason


def test_other_version_recomputes_even_with_identical_parameters(tmp_path):
    cp, _ = _complete(tmp_path)
    decision = checkpoint.reuse_decision(cp, {"alpha": 0.05}, version="9.9.9")
    assert not decision and "version" in decision.reason


def test_changed_parameter_recomputes_and_names_it(tmp_path):
    cp, _ = _complete(tmp_path)
    decision = checkpoint.reuse_decision(cp, {"alpha": 0.1})
    assert not decision and "'alpha'" in decision.reason


def test_added_parameter_recomputes(tmp_path):
    cp, _ = _complete(tmp_path)
    decision = checkpoint.reuse_decision(cp, {"alpha": 0.05, "species": "macaque"})
    assert not decision and "'species'" in decision.reason


def test_deleted_output_recomputes(tmp_path):
    cp, out = _complete(tmp_path)
    out.unlink()
    decision = checkpoint.reuse_decision(cp, {"alpha": 0.05})
    assert not decision and "no longer exists" in decision.reason


def test_tuple_and_list_parameters_compare_equal(tmp_path):
    cp, _ = _complete(tmp_path, params={"species": ("human", "macaque")})
    assert checkpoint.reuse_decision(cp, {"species": ["human", "macaque"]})
