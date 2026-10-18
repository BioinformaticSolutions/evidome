"""Many-to-one ortholog collisions, resolved by a majority vote of every metric.

When several input genes map to the same human gene (e.g. two macaque paralogs
of one human gene), a gene-level analysis needs one value per human gene. Every
metric is calculated and each casts a vote; the majority wins.

How the vote works, per human gene, PER COMPARISON (the winner may differ
between comparisons -- one paralog can carry the diet effect and another the
drug effect):

- Only candidates TESTED in that comparison may vote or win.
- Each available metric votes for the candidate it ranks best:
    directional contrast: smallest P-value, largest |T-statistic|,
                          largest |Fold_change| (or |Log2_fold_change|, the
                          one the block has), largest |Effect|
                          (a burden test votes as a contrast)
    F-test:               smallest ANOVA p, largest F
  A metric absent from the file (e.g. no T-statistic on a plain t-test slice)
  simply does not vote. Q-value does not vote: it is a monotone transform of
  P-value across all genes, so it would only duplicate P-value's vote.
  |Fold_change| means the size of the change on a log scale: for a ratio
  column (1 = no change) it is |log(ratio)|, so 0.25 is as large a change as
  4, and 0 (switched off) is the largest of all.
- Most votes wins. Ties are broken by smallest P-value, then largest
  |Fold_change| (contrasts), then gene ID -- so the result is deterministic.

What the agent should know, and the output records: the metrics are not
independent. P-value and |t| move together, as do |Fold_change| and |Effect|,
so in practice the vote is usually "significance axis vs. magnitude axis",
and a 2-2 split is broken by P-value. Every collision's full ballot -- each
metric's pick, the vote count, whether it was a tie -- is kept, so a reader can
see when the choice was close.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

#: metric role -> (direction, uses absolute value). "min" = smaller is better.
CONTRAST_METRICS: Tuple[Tuple[str, str, bool], ...] = (
    ("P-value", "min", False),
    ("T-statistic", "max", True),
    ("Fold_change", "max", True),
    ("Log2_fold_change", "max", True),     # a block has this or Fold_change, never both
    ("Effect", "max", True),
)
FTEST_METRICS: Tuple[Tuple[str, str, bool], ...] = (
    ("ANOVA_p", "min", False),
    ("F", "max", False),
)


@dataclass
class Ballot:
    human_gene: str
    comparison: str
    candidates: List[str]  # every input gene mapping here, tested or not
    tested: List[str]
    winner: Optional[str]
    votes: Dict[str, int] = field(default_factory=dict)  # candidate -> votes
    picks: Dict[str, str] = field(default_factory=dict)  # metric -> candidate it chose
    tie: bool = False
    decided_by: str = "vote"  # "vote", "tie-break", "only tested candidate", "none tested"


def resolve(
    human_gene: str,
    comparison: str,
    kind: str,
    candidates: Mapping[str, Mapping[str, Optional[float]]],
    fc_ratio: bool = False,
) -> Ballot:
    """Pick one input gene for ``human_gene`` in one comparison.

    ``candidates`` maps input gene ID -> {metric role -> value or None}, where
    the roles are the comparison's column roles (``P-value``, ``Fold_change``,
    ... or ``ANOVA_p``, ``F``). A candidate is tested when its p is not None.

    ``fc_ratio``: the comparison's Fold_change is a plain ratio (0.5 = halved) rather than Paralome's signed
    "times" (-2 = halved). Its magnitude is then |log(FC)|, so a four-fold drop (0.25) outranks a 1.1-fold rise;
    by |FC| it would lose to it.
    """
    def magnitude(role: str, v: float) -> Optional[float]:
        if role == "Fold_change" and fc_ratio:
            # 0 is a gene switched off: the largest possible change
            return abs(math.log(v)) if v > 0 else (math.inf if v == 0 else None)
        return abs(v)
    metrics = FTEST_METRICS if kind == "ftest" else CONTRAST_METRICS
    p_role = metrics[0][0]
    names = sorted(candidates)
    tested = [g for g in names if candidates[g].get(p_role) is not None]
    ballot = Ballot(human_gene, comparison, names, tested, None)

    if not tested:
        ballot.decided_by = "none tested"
        return ballot
    if len(tested) == 1:
        ballot.winner = tested[0]
        ballot.decided_by = "only tested candidate"
        return ballot

    votes = {g: 0 for g in tested}
    for role, direction, absolute in metrics:
        scored = [(g, candidates[g].get(role)) for g in tested]
        scored = [(g, magnitude(role, v) if absolute else v) for g, v in scored if v is not None]
        scored = [(g, v) for g, v in scored if v is not None]
        if not scored:
            continue  # metric absent for this comparison: no vote
        best = min(v for _, v in scored) if direction == "min" else max(v for _, v in scored)
        leaders = [g for g, v in scored if v == best]
        if len(leaders) != 1:
            continue  # the metric cannot tell them apart: it abstains
        ballot.picks[role] = leaders[0]
        votes[leaders[0]] += 1
    ballot.votes = votes

    top = max(votes.values())
    leaders = [g for g in tested if votes[g] == top]
    if len(leaders) == 1:
        ballot.winner = leaders[0]
        return ballot

    ballot.tie = True
    ballot.decided_by = "tie-break"

    def tie_key(g: str):
        values = candidates[g]
        p = values.get(p_role)
        role = "Fold_change" if values.get("Fold_change") is not None else "Log2_fold_change"
        fc = values.get(role)
        size = magnitude(role, fc) if fc is not None else None
        return (p, -size if size is not None else 0.0, g)

    ballot.winner = sorted(leaders, key=tie_key)[0]
    return ballot
