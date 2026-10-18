"""Many-to-one collisions: every metric votes, majority wins, ties are deterministic."""

from evidome.collisions import resolve


def c(p, t=None, fc=None, eff=None):
    return {"P-value": p, "T-statistic": t, "Fold_change": fc, "Effect": eff}


def test_clear_majority_wins():
    b = resolve("ENSG1", "X", "contrast", {
        "A": c(0.001, -4.0, -1.8, -0.8),
        "B": c(0.2, -1.0, -3.0, -1.5),   # bigger fold change, weaker evidence
    })
    # p and t pick A; FC and Effect pick B -> 2-2 tie, broken by smallest p
    assert b.tie and b.decided_by == "tie-break" and b.winner == "A"


def test_three_of_four_metrics():
    b = resolve("ENSG1", "X", "contrast", {
        "A": c(0.001, 4.0, 1.8, 0.9),
        "B": c(0.2, 1.0, 3.0, 0.5),
    })
    assert not b.tie and b.winner == "A" and b.votes == {"A": 3, "B": 1}
    assert b.picks["Fold_change"] == "B"


def test_absent_metric_does_not_vote():
    b = resolve("ENSG1", "X", "contrast", {"A": c(0.01, fc=1.5), "B": c(0.02, fc=2.5)})
    assert set(b.picks) == {"P-value", "Fold_change"} and b.tie and b.winner == "A"


def test_untested_candidate_cannot_win():
    b = resolve("ENSG1", "X", "contrast", {"A": c(None, fc=9.9), "B": c(0.5, 0.1, 1.01, 0.01)})
    assert b.winner == "B" and b.decided_by == "only tested candidate" and b.candidates == ["A", "B"]


def test_nobody_tested():
    b = resolve("ENSG1", "X", "contrast", {"A": c(None), "B": c(None)})
    assert b.winner is None and b.decided_by == "none tested"


def test_ftest_uses_p_and_f():
    b = resolve("ENSG1", "X", "ftest", {
        "A": {"ANOVA_p": 0.3, "F": 1.1}, "B": {"ANOVA_p": 0.01, "F": 7.0}})
    assert b.winner == "B" and b.votes == {"A": 0, "B": 2}


def test_metric_that_cannot_separate_abstains():
    b = resolve("ENSG1", "X", "contrast", {"A": c(0.05, 2.0, 1.5), "B": c(0.05, 3.0, 1.5)})
    assert b.picks == {"T-statistic": "B"} and b.winner == "B"


def test_three_way_is_deterministic_regardless_of_input_order():
    cands = {"C": c(0.01, 3.0, 2.0), "A": c(0.01, 3.0, 2.0), "B": c(0.02, 1.0, 1.2)}
    assert resolve("G", "X", "contrast", cands).winner == "A"
    assert resolve("G", "X", "contrast", dict(reversed(list(cands.items())))).winner == "A"
