"""Gene-set enrichment on the per-comparison human-gene tables.

Two methods, implemented here in the standard library (validated against
GSEApy; no runtime dependency):

Over-representation (ORA), every comparison
    Universe: the human genes TESTED in that comparison. Selected: those with
    raw p < alpha (the same threshold as the crosswalk summary). One-sided
    hypergeometric test per gene set; contrasts are tested three ways (up,
    down, all), F-tests only "all" unless their contrasts gave them a
    direction (wide_stats.derive_ftest_directions), then all three.
    BH q-values within (comparison, direction, source).

Preranked GSEA, every comparison with a direction (contrasts; F-tests directed
by their contrasts, ranked by sign(Direction_estimate) x -log10(ANOVA p))
    Genes ranked by the T-statistic where the table has one, else by
    sign x -log10(p). Enrichment score as Subramanian et al. 2005 (weight 1).
    Significance from gene-set permutation: for a set of size s the null is
    the ES of random size-s sets from the same ranked list, so ONE null per
    distinct size serves every set of that size (as fgsea's simple mode
    does); nulls are computed in parallel. NES = ES / mean |null ES| of the
    same sign; p = (#null at least as extreme + 1) / (#null of that sign + 1);
    BH q within (comparison, source). F-tests that their contrasts give no direction: skipped.

    Adaptive permutations. With a fixed 1,000 permutations the smallest
    reachable p is ~0.002, and BH across 13,784 GO sets then needs ~550 sets
    at that floor before any reaches q < 0.05 -- so GSEA could find nothing in
    GO or WikiPathways for that reason alone.
    So: every size starts with ``permutations``; any size where some set is
    beaten by fewer than ``min_exceedances`` nulls gets ten times as many,
    repeatedly, up to ``max_permutations``. Extra permutations come in fixed
    blocks of 10,000, each with its own seed, so the result does not depend on
    --threads. Each row records the permutations its p rests on.

Direction of a gene: sign of the T-statistic if present, else of Effect, else
of Fold_change read in the convention the column itself shows -- signed
"times" (Paralome bulk: +2 up, -2 down, 1 unchanged), a plain ratio (0.5
down), or a log fold change (Paralome single-cell: log2, 0 unchanged; a
column with negative values and values between -1 and 1); the log form is
the most common of all.
A Log2_fold_change column is declared log2 and is read as log2,
never judged from its values: a burden table whose genes all have a little
more than expected (every value in (0, 1)) would otherwise read as ratios
below 1, i.e. down. A burden test (kind 'burden': observed against expected
carriers) is read like a contrast; up = more carriers than expected.

Burden tests (kind 'burden')
    A gene-burden comparison is not an expression change, and two things
    follow. (1) No direction-of-effect statistic is reported for it: the
    pathway activation z (and the regulators' activation z, bias-corrected z,
    activity t and the regulator's own direction) are left blank, because a
    burden test measures damage counts, not expression direction; the up /
    down / all ORA split is kept (up = more carriers than expected). (2) The
    ORA null is weighted by DETECTABILITY. Under the burden test's null, a
    gene's chance of reaching p < alpha depends on how many qualifying alleles
    gnomAD shows for it -- it differs between genes by orders of magnitude --
    so gene sets rich in such genes are "enriched" by construction. The
    producer gives, per gene, that probability (Null_selection_probability,
    wide_stats.py) and the ORA p of a set is then the upper tail of Wallenius'
    noncentral hypergeometric distribution, with odds = the mean probability
    of the set's genes over the mean of the other genes of the universe (the
    approximation goseq uses for gene-length bias; Young et al. 2010, Genome
    Biology 11:R14). Each ORA row says which null its p comes from
    (null_model); a burden comparison without the column -- or run at an alpha
    the probabilities were not computed for -- keeps the hypergeometric and is
    marked "unadjusted for detectability". GSEA is not run for a burden test:
    its ranking (sign x -log10 p) is confounded with the same detectability,
    and no sound weighted version is implemented.

Sets are tested when their size WITHIN THE UNIVERSE is in [min_size,
max_size] (for GSEA: within the ranked list, which under a sign x -log10(p)
metric leaves out tested genes with no direction; a T-statistic ranks every
gene that has one, T = 0 included); ORA's sets outside the range are counted
per source in the manifest (sets_outside_size_range), GSEA's are not. Every tested
set is reported -- nothing is filtered by q.

Same inputs, same seed, same version: byte-identical output, whatever
--threads is (seeds are derived from the seed, the comparison and the size).
"""

from __future__ import annotations

import csv
import math
import random
from array import array
from itertools import accumulate, repeat
from operator import mul, sub
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .wide_stats import NULL_SELECTION_ALPHA, NULL_SELECTION_ROLE, is_directional, p_role_of

DEFAULT_MIN_SIZE = 10
DEFAULT_MAX_SIZE = 500
DEFAULT_PERMUTATIONS = 1000
DEFAULT_SEED = 1
DEFAULT_MAX_PERMUTATIONS = 100000
DEFAULT_MIN_EXCEEDANCES = 10

ORA_COLUMNS = ("comparison", "direction", "source", "set_id", "name", "category", "set_size", "universe",
               "selected", "overlap", "expected", "fold_enrichment", "p", "q", "overlap_genes",
               "pattern_genes", "pattern_consistent", "pattern_inconsistent", "pathway_activation_z",
               "pattern_consistent_genes", "pattern_inconsistent_genes", "null_model")

#: The ORA null each row's p comes from (column null_model).
NULL_HYPERGEOMETRIC = "hypergeometric"
NULL_WALLENIUS = "Wallenius noncentral hypergeometric, weighted by detectability"
NULL_UNADJUSTED = "hypergeometric (unadjusted for detectability)"
#: Why a burden comparison has no direction-of-effect statistic, and no GSEA.
BURDEN_NO_DIRECTION = "not reported: a burden test measures damage counts, not expression direction"
BURDEN_NO_GSEA = ("burden test: GSEA ranks genes by sign x -log10 p, which is confounded with each gene's "
                  "detectability (genes with more common qualifying alleles reach small p more easily); not run")
GSEA_COLUMNS = ("comparison", "source", "set_id", "name", "category", "set_size", "es", "nes", "p", "q",
                "permutations", "leading_edge_size", "leading_edge_genes", "rank_metric")


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _log_choose(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k), X ~ Hypergeometric(population N, successes K, draws n). Exact, summed from the smaller tail."""
    lo, hi = max(0, n - (N - K)), min(K, n)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0

    def pmf(x: int) -> float:
        return math.exp(_log_choose(K, x) + _log_choose(N - K, n - x) - _log_choose(N, n))

    mode = int((n + 1) * (K + 1) / (N + 2))
    if k > mode:  # upper tail is the small one: sum k..hi with the term recurrence
        term = pmf(k)
        total, x = term, k
        while x < hi:
            term *= (K - x) * (n - x) / ((x + 1) * (N - K - n + x + 1))
            x += 1
            total += term
            if term < total * 1e-17:
                break
        return min(1.0, total)
    # lower tail (x < k) is the small one: 1 - sum lo..k-1, summed downward from k-1
    term = pmf(k - 1)
    total, x = term, k - 1
    while x > lo:
        term *= x * (N - K - n + x) / ((K - x + 1) * (n - x + 1))
        x -= 1
        total += term
        if term < total * 1e-17:
            break
    return max(0.0, 1.0 - total)


def _gauss_legendre(order: int) -> List[Tuple[float, float]]:
    """(node, weight) pairs of Gauss-Legendre quadrature on [-1, 1], by Newton's method on P_order."""
    nodes = []
    for i in range(1, order + 1):
        t = math.cos(math.pi * (i - 0.25) / (order + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, t
            for k in range(2, order + 1):
                p0, p1 = p1, ((2 * k - 1) * t * p1 - (k - 1) * p0) / k
            dp = order * (t * p1 - p0) / (t * t - 1.0)
            step = p1 / dp
            t -= step
            if abs(step) < 1e-16:
                break
        nodes.append((t, 2.0 / ((1.0 - t * t) * dp * dp)))
    return nodes


GAUSS_NODES = _gauss_legendre(24)


def _wallenius_log_integral(x: int, n: int, m1: int, m2: int, omega: float) -> float:
    """log of Fog's integral for the Wallenius pmf, d * int_0^1 (1 - t^(omega/d))^x (1 - t^(1/d))^(n-x) dt.

    With t = exp(d v) the integrand is exp(g(v)), g(v) = d v + x log(1 - e^(omega v)) + (n-x) log(1 - e^v),
    v < 0, which is concave (a sum of concave terms) and so has one peak. The peak is found by bisection on
    g' (decreasing), and exp(g - g_max) is integrated outwards from it by Gauss-Legendre quadrature in steps
    of its own width, until a step adds nothing at double precision.
    """
    d = omega * (m1 - x) + (m2 - (n - x))
    if d <= 0:
        return 0.0 if x == m1 and n - x == m2 else -math.inf

    def g(v: float) -> float:
        total = d * v
        if x:
            a = -math.expm1(omega * v)
            if a <= 0:
                return -math.inf
            total += x * math.log(a)
        if n - x:
            b = -math.expm1(v)
            if b <= 0:
                return -math.inf
            total += (n - x) * math.log(b)
        return total

    def slope(v: float) -> float:
        out = d
        if x:
            e = math.exp(omega * v)
            out -= x * omega * e / -math.expm1(omega * v)
        if n - x:
            e = math.exp(v)
            out -= (n - x) * e / -math.expm1(v)
        return out

    if n == 0:
        return 0.0                      # d * int_0^1 t^0 ... over v: d * int exp(d v) dv = 1
    lo, hi = -1.0, -1e-300
    while slope(lo) < 0:
        lo *= 2.0
        if lo < -1e6:
            break
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if slope(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-15 * max(1.0, abs(lo)):
            break
    peak = 0.5 * (lo + hi)
    top = g(peak)
    # width from the curvature at the peak (Laplace); never wider than the distance to 0 allows
    curve = 0.0
    if x:
        e = math.exp(omega * peak)
        curve += x * omega * omega * e / (math.expm1(omega * peak) ** 2)
    if n - x:
        e = math.exp(peak)
        curve += (n - x) * e / (math.expm1(peak) ** 2)
    width = 1.0 / math.sqrt(curve) if curve > 0 else 1.0 / d
    width = min(width, max(abs(peak), 1e-300))

    def f(v: float) -> float:
        value = g(v)
        return 0.0 if value == -math.inf else math.exp(value - top)

    def piece(a: float, b: float) -> float:
        """Gauss-Legendre over [a, b] (GAUSS_NODES points)."""
        half, mid = 0.5 * (b - a), 0.5 * (b + a)
        return half * math.fsum(w * f(mid + half * t) for t, w in GAUSS_NODES)

    area = 0.0
    # left of the peak: v -> -infinity (g falls at least linearly)
    a = peak
    for _ in range(100000):
        b = a - width
        part = piece(b, a)
        area += part
        a = b
        if part <= area * 1e-17:
            break
        width *= 1.5 if part < area * 1e-3 else 1.0
    # right of the peak: up to v = 0, in steps that shrink with the distance left
    a = peak
    step = width
    for _ in range(100000):
        room = -a
        if room <= 0:
            break
        b = a + min(step, room) if room > step * 1.5 else 0.0
        part = piece(a, b)
        area += part
        a = b
        if b >= 0.0 or part <= area * 1e-17:
            break
    return math.log(d) + top + math.log(area) if area > 0 else -math.inf


def wallenius_pmf(x: int, n: int, m1: int, m2: int, omega: float) -> float:
    """P(X = x): x of the n balls drawn one by one without replacement are from the m1 (weight omega each)
    rather than the m2 (weight 1 each) -- Wallenius' noncentral hypergeometric distribution."""
    lo, hi = max(0, n - m2), min(n, m1)
    if x < lo or x > hi:
        return 0.0
    if omega == 1.0:
        return math.exp(_log_choose(m1, x) + _log_choose(m2, n - x) - _log_choose(m1 + m2, n))
    log_value = _log_choose(m1, x) + _log_choose(m2, n - x) + _wallenius_log_integral(x, n, m1, m2, omega)
    return math.exp(log_value) if log_value > -745 else 0.0


#: The odds of a weighted ORA are kept within these bounds (a set, or the rest of the universe, whose genes all
#: have probability 0 would give 0 or infinite odds).
ODDS_BOUNDS = (1e-12, 1e12)


def wallenius_sf(k: int, N: int, K: int, n: int, omega: float) -> float:
    """P(X >= k), X ~ Wallenius(K weighted omega among N, n drawn). omega = 1 is hypergeom_sf, exactly."""
    if omega == 1.0:
        return hypergeom_sf(k, N, K, n)
    omega = min(max(omega, ODDS_BOUNDS[0]), ODDS_BOUNDS[1])
    m1, m2 = K, N - K
    lo, hi = max(0, n - m2), min(n, m1)
    if k <= lo:
        return 1.0
    if k > hi:
        return 0.0
    here, below = wallenius_pmf(k, n, m1, m2, omega), wallenius_pmf(k - 1, n, m1, m2, omega)
    if here <= below:
        # at or past the mode: the upper tail is the small one, its terms falling from k
        total, x, term = here, k, here
        while x < hi:
            x += 1
            term = wallenius_pmf(x, n, m1, m2, omega)
            total += term
            if term <= total * 1e-17:
                break
        return min(1.0, total)
    # before the mode: 1 - the lower tail, its terms falling from k-1 downwards
    total, x, term = below, k - 1, below
    while x > lo:
        x -= 1
        term = wallenius_pmf(x, n, m1, m2, omega)
        total += term
        if term <= total * 1e-17:
            break
    return min(1.0, max(0.0, 1.0 - total))


def detectability_odds(members: Sequence[str], universe: Sequence[str], weight: Dict[str, float],
                       total: Optional[float] = None) -> float:
    """Mean null selection probability of the set's genes over that of the rest of the universe (goseq's odds).

    ``members`` are within ``universe``; ``total`` (the universe's sum of weights) may be passed to save the sum.
    """
    inside = [weight[g] for g in members]
    rest = len(universe) - len(inside)
    if not inside or rest <= 0:
        return 1.0
    total = math.fsum(weight[g] for g in universe) if total is None else total
    mean_in, mean_out = math.fsum(inside) / len(inside), max(total - math.fsum(inside), 0.0) / rest
    if mean_out <= 0:
        return ODDS_BOUNDS[1] if mean_in > 0 else 1.0
    return min(max(mean_in / mean_out, ODDS_BOUNDS[0]), ODDS_BOUNDS[1])


def bh(pvalues: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg q-values, in input order."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * n / rank)
        q[i] = running
    return q


def enrichment_score(weights: Sequence[float], hits: Sequence[int]) -> Tuple[float, int]:
    """Weighted KS running sum (Subramanian 2005, p=1). ``hits`` = sorted rank positions; returns (ES, peak index).

    Computed with exactly null_es's expressions, so an observed set and an identical null draw score the same to
    the last bit; different arithmetic could differ by one ulp on a tie -- enough to give them opposite signs. A
    positive and a negative extreme of equal size: the positive wins, in both."""
    N, s = len(weights), len(hits)
    ws = [abs(weights[h]) for h in hits]
    nr = sum(ws)
    if s == 0 or s == N or nr == 0:
        return 0.0, -1
    inv, miss = 1.0 / nr, 1.0 / (N - s)
    after = list(map(sub, map(mul, accumulate(ws), repeat(inv)), map(mul, map(sub, hits, range(s)), repeat(miss))))
    before = list(map(sub, after, map(mul, ws, repeat(inv))))
    top, bottom = max(after), min(before)
    if top >= -bottom - TIE_TOLERANCE:
        return top, after.index(top)
    return bottom, before.index(bottom) - 1


#: Extremes of equal size are equal to within rounding: 7/36 summed two ways can differ in the last bit, which
#: would give an exact tie the negative sign. Running sums are at most 1 in size.
TIE_TOLERANCE = 1e-12


#: Permutations added beyond the first round come in blocks of this size, each
#: with its own seed -- a fixed partition, so results never depend on --threads.
PERMUTATION_BLOCK = 10000


def null_es(abs_weights: Sequence[float], hits: Sequence[int]) -> float:
    """``enrichment_score(...)[0]`` for the permutation null, with the per-gene work in C.

    Same arithmetic (checked to 1e-16 against the loop): the running sum just
    after hit i is cum_i / NR - (h_i - i) / (N - s); just before it, that minus
    w_i / NR. The loop costs ~1.2x this.
    """
    N, s = len(abs_weights), len(hits)
    ws = [abs_weights[h] for h in hits]
    nr = sum(ws)
    if s == 0 or s == N or nr == 0:
        return 0.0
    inv, miss = 1.0 / nr, 1.0 / (N - s)
    after = list(map(sub, map(mul, accumulate(ws), repeat(inv)), map(mul, map(sub, hits, range(s)), repeat(miss))))
    top = max(after)
    bottom = min(map(sub, after, map(mul, ws, repeat(inv))))
    return top if top >= -bottom - TIE_TOLERANCE else bottom


def random_positions(rng: random.Random, n: int, k: int) -> List[int]:
    """``sorted(rng.sample(range(n), k))`` in distribution, with the random bits drawn in one call.

    Drawing random sets is most of GSEA's time. ``random.sample`` makes one Python-level call per position; this
    draws 32-bit words in bulk, reduces them mod n (bias < n / 2**32, ~4e-6
    for 18,444 genes), keeps first occurrences, and draws more on the rare
    shortfall. About 2x faster. Deterministic for a given generator state.
    """
    m = k + k // 8 + 8
    while True:
        words = array("I")
        words.frombytes(rng.randbytes(4 * m))
        unique = list(dict.fromkeys([w % n for w in words]))
        if len(unique) >= k:
            return sorted(unique[:k])
        m *= 2


def _null_task(task: Tuple[str, List[float], List[Tuple[int, int, str]]]) -> List[Tuple[int, List[float], List[float]]]:
    """Process-pool worker: for each (size, permutations, seed tag) job, the null ES split by sign."""
    comparison, weights, jobs = task
    abs_weights = [abs(w) for w in weights]
    n = len(weights)
    out = []
    for size, permutations, tag in jobs:
        rng = random.Random(tag)
        pos, neg = array("d"), array("d")
        for _ in range(permutations):
            es = null_es(abs_weights, random_positions(rng, n, size))
            (pos if es >= 0 else neg).append(es)
        out.append((size, pos, neg))   # arrays: 8 bytes per value in memory and across processes
    return out


class Null:
    """One set size's null distribution, sorted, for counting exceedances by bisection."""

    def __init__(self) -> None:
        self.pos = array("d")
        self.neg = array("d")
        self.permutations = 0

    def add(self, pos, neg, permutations: int) -> None:
        self.pos = array("d", sorted(self.pos + array("d", pos)))
        self.neg = array("d", sorted(self.neg + array("d", neg)))
        self.permutations += permutations

    def exceedances(self, es: float) -> Tuple[int, int]:
        """(nulls at least as extreme as ``es``, nulls of the same sign)."""
        import bisect

        if es >= 0:
            return len(self.pos) - bisect.bisect_left(self.pos, es), len(self.pos)
        return bisect.bisect_right(self.neg, es), len(self.neg)

    def scale(self, es: float) -> float:
        same = self.pos if es >= 0 else self.neg
        return sum(abs(v) for v in same) / len(same) if same else 0.0


def _first_tag(seed: int, comparison: str, size: int) -> str:
    return f"{seed}:{comparison}:{size}"


def _block_tag(seed: int, comparison: str, size: int, block: int) -> str:
    return f"{seed}:{comparison}:{size}:block{block}"


# ---------------------------------------------------------------------------
# reading the crosswalk output
# ---------------------------------------------------------------------------

@dataclass
class ComparisonTable:
    name: str
    kind: str
    genes: List[str]                       # tested human genes, sorted
    symbols: Dict[str, str]
    p: Dict[str, float]
    direction: Dict[str, int]              # +1 / -1 / 0 (contrasts; F-tests directed by their contrasts)
    metric: Dict[str, float]               # ranking metric (the same)
    direction_from: str = ""
    metric_name: str = ""
    notes: List[str] = field(default_factory=list)
    #: A burden test's per-gene null selection probability (wide_stats.NULL_SELECTION_ROLE), if given.
    selection: Dict[str, float] = field(default_factory=dict)
    selection_total: Optional[float] = None


def _num(text: str) -> Optional[float]:
    text = (text or "").strip()
    if not text:
        return None
    value = float(text)
    return None if math.isnan(value) else value


def fold_change_convention(values: Sequence[float]) -> str:
    """Judged from the column itself: 'ratio' (0.5 = halved; never negative, some values between 0 and 1),
    'log' (-1 = halved in log2; negative values and values between -1 and 1), else 'signed-times' (Paralome
    bulk: -2 = halved, no value between -1 and 1). One rule for the whole program (human_tables uses it too)."""
    negative = any(v < 0 for v in values)
    fraction = any(0 < abs(v) < 1 for v in values)
    if fraction and not negative:
        return "ratio"
    return "log" if fraction else "signed-times"


def read_comparison(path: Path, name: str, kind: str) -> ComparisonTable:
    p_col = p_role_of(kind)
    with Path(path).open(encoding="utf-8", newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r.get("tested") == "yes"]
    rows.sort(key=lambda r: r["human_gene_id"])
    table = ComparisonTable(name, kind, [r["human_gene_id"] for r in rows],
                            {r["human_gene_id"]: r["human_symbol"] or r["human_gene_id"] for r in rows},
                            {r["human_gene_id"]: float(r[p_col]) for r in rows}, {}, {})
    header = rows[0].keys() if rows else ()
    if NULL_SELECTION_ROLE in header:
        table.selection = {r["human_gene_id"]: float(r[NULL_SELECTION_ROLE]) for r in rows
                           if _num(r[NULL_SELECTION_ROLE]) is not None}
    if not is_directional(kind):
        # an F-test has a direction when its contrasts gave it one (wide_stats.derive_ftest_directions)
        if "Direction_estimate" in header and any(_num(r["Direction_estimate"]) is not None for r in rows):
            sign = {r["human_gene_id"]: _num(r["Direction_estimate"]) for r in rows}
            table.direction_from = "Direction_estimate"
            for g in table.genes:
                v = sign.get(g)
                table.direction[g] = 0 if v is None or v == 0 else (1 if v > 0 else -1)
                if table.direction[g]:
                    table.metric[g] = table.direction[g] * -math.log10(max(table.p[g], 1e-300))
            table.metric_name = "sign(Direction_estimate) x -log10(ANOVA p)"
        return table
    if "T-statistic" in header and any(_num(r["T-statistic"]) is not None for r in rows):
        source, sign = "T-statistic", {r["human_gene_id"]: _num(r["T-statistic"]) for r in rows}
    elif "Effect" in header and any(_num(r["Effect"]) is not None for r in rows):
        source, sign = "Effect", {r["human_gene_id"]: _num(r["Effect"]) for r in rows}
    elif "Log2_fold_change" in header:
        # declared log2: 0 is unchanged, whatever the values look like -- never re-guessed
        source, sign = "Log2_fold_change", {r["human_gene_id"]: _num(r["Log2_fold_change"]) for r in rows}
    else:
        fc = {r["human_gene_id"]: _num(r["Fold_change"]) for r in rows}
        convention = fold_change_convention([v for v in fc.values() if v is not None])
        pivot = 1.0 if convention == "ratio" else 0.0
        source = f"Fold_change ({convention})"
        # 1 is "unchanged" for a ratio and for signed times, but a two-fold rise in log2
        sign = {g: (None if v is None or (v == 1.0 and convention != "log")
                    else v - pivot if convention == "ratio" else v)
                for g, v in fc.items()}
    table.direction_from = source
    for g in table.genes:
        v = sign.get(g)
        table.direction[g] = 0 if v is None or v == 0 else (1 if v > 0 else -1)
    if source == "T-statistic":
        table.metric_name = "T-statistic"
        table.metric = {g: sign[g] for g in table.genes if sign.get(g) is not None}
    else:
        table.metric_name = f"sign({source}) x -log10(p)"
        for g in table.genes:
            if table.direction[g]:
                table.metric[g] = table.direction[g] * -math.log10(max(table.p[g], 1e-300))
    return table


# ---------------------------------------------------------------------------
# the analysis
# ---------------------------------------------------------------------------

@dataclass
class Settings:
    alpha: float = 0.05
    min_size: int = DEFAULT_MIN_SIZE
    max_size: int = DEFAULT_MAX_SIZE
    permutations: int = DEFAULT_PERMUTATIONS
    seed: int = DEFAULT_SEED
    max_permutations: int = DEFAULT_MAX_PERMUTATIONS
    min_exceedances: int = DEFAULT_MIN_EXCEEDANCES


def _in_universe(members: Dict[str, List[str]], universe: set) -> Dict[str, List[str]]:
    return {sid: [g for g in genes if g in universe] for sid, genes in members.items()}


def ora_null(table: ComparisonTable, alpha: float) -> str:
    """Which null a comparison's ORA (and regulator overlap) p-values come from."""
    if table.kind != "burden":
        return NULL_HYPERGEOMETRIC
    usable = (table.selection and all(g in table.selection for g in table.genes)
              and math.isclose(alpha, NULL_SELECTION_ALPHA))
    return NULL_WALLENIUS if usable else NULL_UNADJUSTED


def ora_null_note(table: ComparisonTable, alpha: float) -> str:
    """Why a burden comparison's ORA is unadjusted, in words ('' when it is not)."""
    if ora_null(table, alpha) != NULL_UNADJUSTED:
        return ""
    if not table.selection:
        return f"the input gives no {NULL_SELECTION_ROLE} for this burden comparison"
    if not math.isclose(alpha, NULL_SELECTION_ALPHA):
        return f"{NULL_SELECTION_ROLE} is P(p < {NULL_SELECTION_ALPHA}); this run selects at alpha {alpha}"
    return f"some tested genes have no {NULL_SELECTION_ROLE}"


def overlap_p(table: ComparisonTable, null: str, k: int, N: int, members: Sequence[str], n: int) -> float:
    """One-sided ORA p of k selected among a set's members (all within the universe): hypergeometric, or
    Wallenius by detectability."""
    if not n:
        return 1.0
    if null == NULL_WALLENIUS:
        if k <= max(0, n - (N - len(members))):
            return 1.0
        if table.selection_total is None:
            table.selection_total = math.fsum(table.selection[g] for g in table.genes)
        odds = detectability_odds(members, table.genes, table.selection, table.selection_total)
        return wallenius_sf(k, N, len(members), n, odds)
    return hypergeom_sf(k, N, len(members), n)


def run_ora(table: ComparisonTable, info: Dict[str, Dict[str, object]], members: Dict[str, List[str]],
            settings: Settings) -> Tuple[List[Dict[str, object]], Dict[str, int]]:
    universe = set(table.genes)
    N = len(universe)
    within = _in_universe(members, universe)
    testable = {sid: g for sid, g in within.items() if settings.min_size <= len(g) <= settings.max_size}
    skipped: Dict[str, int] = defaultdict(int)
    for sid, g in within.items():
        if sid not in testable:
            skipped[str(info[sid]["source"])] += 1
    selections = {"all": {g for g in table.genes if table.p[g] < settings.alpha}}
    if table.direction_from:          # a contrast, or an F-test given a direction by its contrasts
        selections["up"] = {g for g in selections["all"] if table.direction[g] > 0}
        selections["down"] = {g for g in selections["all"] if table.direction[g] < 0}
    null = ora_null(table, settings.alpha)
    rows: List[Dict[str, object]] = []
    for direction in ("up", "down", "all"):
        if direction not in selections:
            continue
        chosen = selections[direction]
        n = len(chosen)
        by_source: Dict[str, List[Dict[str, object]]] = defaultdict(list)
        for sid in sorted(testable):
            genes = testable[sid]
            K = len(genes)
            hit = sorted(g for g in genes if g in chosen)
            k = len(hit)
            expected = K * n / N if N else 0.0
            by_source[str(info[sid]["source"])].append({
                "comparison": table.name, "direction": direction, "source": info[sid]["source"], "set_id": sid,
                "name": info[sid]["name"], "category": info[sid]["category"], "set_size": K, "universe": N,
                "selected": n, "overlap": k, "expected": round(expected, 4),
                "fold_enrichment": round(k / expected, 4) if expected else "",
                "p": overlap_p(table, null, k, N, genes, n),
                "overlap_genes": ";".join(sorted(table.symbols[g] for g in hit)),
                "null_model": null,
            })
        for source in sorted(by_source):
            group = by_source[source]
            for row, q in zip(group, bh([float(r["p"]) for r in group])):
                row["q"] = q
                # no activation z for a burden test (BURDEN_NO_DIRECTION)
                if direction == "all" and table.kind != "burden":
                    row.update(pathway_activation(table, info[row["set_id"]].get("pattern"), chosen,
                                                  members=testable[row["set_id"]]))
            rows.extend(group)
    return rows, dict(skipped)


def pathway_activation(table: ComparisonTable, pattern: Optional[Dict[str, int]], significant,
                       members=None) -> Dict[str, object]:
    """Pathway activation z over the set's significant genes that have an activity pattern.

    consistent: pattern sign x the gene's direction = +1. z = (consistent - inconsistent) / sqrt(N), given when
    N >= PATTERN_MIN; blank columns for a set with no pattern or a comparison with no direction."""
    from .pathway_signs import PATTERN_MIN

    if not pattern:
        return {}
    # Only the set's own tested members. The pattern comes from BioPAX proteins mapped through UniProt, the
    # members from Ensembl2Reactome; where the two mappings differ, genes outside the set would otherwise be
    # counted and could carry the whole z.
    tested = [g for g in pattern if g in table.p and (members is None or g in members)]
    out: Dict[str, object] = {"pattern_genes": len(tested)}
    if not table.direction_from:
        return out
    consistent = sorted(table.symbols[g] for g in tested if g in significant and pattern[g] * table.direction[g] > 0)
    inconsistent = sorted(table.symbols[g] for g in tested if g in significant and pattern[g] * table.direction[g] < 0)
    n = len(consistent) + len(inconsistent)
    out.update({"pattern_consistent": len(consistent), "pattern_inconsistent": len(inconsistent),
                "pathway_activation_z": round((len(consistent) - len(inconsistent)) / math.sqrt(n), 6)
                if n >= PATTERN_MIN else "",
                "pattern_consistent_genes": ";".join(consistent),
                "pattern_inconsistent_genes": ";".join(inconsistent)})
    return out


def gsea_inputs(table: ComparisonTable) -> Tuple[List[str], List[float]]:
    ranked = sorted(table.metric, key=lambda g: (-table.metric[g], g))
    return ranked, [table.metric[g] for g in ranked]


def observed_scores(table: ComparisonTable, members, settings: Settings) -> Dict[str, Tuple[int, float, int, List[int]]]:
    """{set_id: (size, ES, peak, hit positions)} for every set in the size range."""
    ranked, weights = gsea_inputs(table)
    position = {g: i for i, g in enumerate(ranked)}
    scores = {}
    for sid in sorted(members):
        hits = sorted(position[g] for g in members[sid] if g in position)
        if settings.min_size <= len(hits) <= settings.max_size:
            es, peak = enrichment_score(weights, hits)
            scores[sid] = (len(hits), es, peak, hits)
    return scores


def sizes_needing_more(scores, nulls: Dict[int, Null], settings: Settings) -> List[int]:
    """Sizes where some set is beaten by fewer than ``min_exceedances`` nulls and the cap is not reached."""
    need = set()
    for size, es, _, _ in scores.values():
        null = nulls[size]
        if size in need or null.permutations >= settings.max_permutations:
            continue
        if null.exceedances(es)[0] < settings.min_exceedances:
            need.add(size)
    return sorted(need)


def run_gsea(table: ComparisonTable, info, members, settings: Settings,
             nulls: Dict[int, Null], scores=None) -> List[Dict[str, object]]:
    ranked, weights = gsea_inputs(table)
    scores = scores if scores is not None else observed_scores(table, members, settings)
    rows_by_source: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for sid in sorted(scores):
        s, es, peak, hits = scores[sid]
        null = nulls[s]
        extreme, same = null.exceedances(es)
        scale = null.scale(es)
        nes = es / scale if scale else 0.0
        p = (extreme + 1) / (same + 1)
        if es >= 0:
            edge = [ranked[h] for h in hits[: peak + 1]]
        else:
            edge = [ranked[h] for h in hits[peak + 1:]]
        rows_by_source[str(info[sid]["source"])].append({
            "comparison": table.name, "source": info[sid]["source"], "set_id": sid, "name": info[sid]["name"],
            "category": info[sid]["category"], "set_size": s, "es": round(es, 6), "nes": round(nes, 6), "p": p,
            "permutations": null.permutations, "leading_edge_size": len(edge),
            "leading_edge_genes": ";".join(table.symbols[g] for g in edge), "rank_metric": table.metric_name,
        })
    rows: List[Dict[str, object]] = []
    for source in sorted(rows_by_source):
        group = rows_by_source[source]
        for row, q in zip(group, bh([float(r["p"]) for r in group])):
            row["q"] = q
        rows.extend(group)
    return rows


def _write(path: Path, columns: Sequence[str], rows: Sequence[Dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (repr(v) if isinstance(v, float) else v) for k, v in row.items()})
    tmp.replace(path)
    return path


def _comparisons(crosswalk_dir: Path) -> List[Tuple[str, str]]:
    """[(name, kind)] in the crosswalk's order; EVD-0021 if this is not a crosswalk output directory."""
    from .jsonio import read_json

    summary_path = Path(crosswalk_dir) / "comparison_summary.json"
    if not summary_path.exists():
        raise ProjectError(21, f"{summary_path} not found")
    comps = [(c["comparison"], c["kind"]) for c in read_json(summary_path)["comparisons"]]
    for name, _ in comps:
        path = Path(crosswalk_dir) / "comparisons" / f"{name}.human_genes.csv"
        if not path.exists():
            raise ProjectError(21, f"{path} not found")
    return comps


def enrich_comparison(name: str, kind: str, crosswalk_dir: Path, info, members, out: Path, settings: Settings,
                      budget: int, log=print) -> Dict[str, object]:
    """Everything for ONE comparison -- ORA, GSEA with its adaptive rounds -- written under ``out/<name>/``.

    The unit of work for both backends: locally the comparisons run
    one after another, each using every CPU; under SLURM each is one array
    task on whichever node SLURM gives it. Same seeds either way, so the files
    are identical.
    """
    from .jsonio import write_json

    t = read_comparison(Path(crosswalk_dir) / "comparisons" / f"{name}.human_genes.csv", name, kind)
    ora_rows, skipped = run_ora(t, info, members, settings)
    ora_path = _write(out / t.name / "ora.csv", ORA_COLUMNS, ora_rows)
    entry: Dict[str, object] = {
        "comparison": t.name, "kind": t.kind, "universe": len(t.genes),
        "selected": sum(1 for g in t.genes if t.p[g] < settings.alpha),
        "sets_outside_size_range": skipped, "ora": str(ora_path.relative_to(out)),
        "ora_sets_q_below_0.05": _count_q(ora_rows), "ora_null": ora_null(t, settings.alpha)}
    if ora_null_note(t, settings.alpha):
        entry["ora_null_note"] = ora_null_note(t, settings.alpha)
    if t.kind == "burden":
        entry["direction_statistics"] = BURDEN_NO_DIRECTION
    if t.metric and t.kind != "burden":
        weights = {t.name: gsea_inputs(t)[1]}
        scores = observed_scores(t, members, settings)
        nulls: Dict[str, Dict[int, Null]] = {t.name: defaultdict(Null)}
        first = [(t.name, size, settings.permutations, _first_tag(settings.seed, t.name, size))
                 for size in sorted({v[0] for v in scores.values()})]
        _run_null_jobs(first, weights, nulls, budget, log, f"{t.name}: GSEA round 1")
        round_no = 1
        # One block counter per set size across ALL rounds. Recomputing it from the count drawn so far would give
        # round 3 the same block numbers -- and so the same seeds and the same null draws -- as round 2 whenever
        # round 2 drew fewer than PERMUTATION_BLOCK.
        next_block: Dict[int, int] = defaultdict(int)
        while True:
            jobs = []
            for size in sizes_needing_more(scores, nulls[t.name], settings):
                have = nulls[t.name][size].permutations
                target = min(settings.max_permutations, have * 10)
                block = next_block[size]
                while have < target:
                    n = min(PERMUTATION_BLOCK, target - have)
                    jobs.append((t.name, size, n, _block_tag(settings.seed, t.name, size, block)))
                    have += n
                    block += 1
                next_block[size] = block
            if not jobs:
                break
            round_no += 1
            _run_null_jobs(jobs, weights, nulls, budget, log, f"{t.name}: GSEA round {round_no}")
        gsea_rows = run_gsea(t, info, members, settings, nulls[t.name], scores)
        gsea_path = _write(out / t.name / "gsea.csv", GSEA_COLUMNS, gsea_rows)
        entry.update({"gsea": str(gsea_path.relative_to(out)), "direction_from": t.direction_from,
                      "rank_metric": t.metric_name, "ranked_genes": len(t.metric), "gsea_rounds": round_no,
                      "gsea_permutations": sum(n.permutations for n in nulls[t.name].values()),
                      "gsea_sets_q_below_0.05": _count_q(gsea_rows)})
    else:
        # an earlier run's gsea.csv would otherwise be read as this run's
        (out / t.name / "gsea.csv").unlink(missing_ok=True)
        entry["gsea"] = None
        entry["gsea_skipped"] = ("F-test: no direction to rank by" if t.kind == "ftest" else
                                 BURDEN_NO_GSEA if t.kind == "burden" else "no ranking metric")
    write_json(out / t.name / "result.json", entry)
    return entry


def finalize(crosswalk_dir: Path, reference_dir: Path, out: Path, release: int, settings: Settings,
             comparisons: Sequence[Tuple[str, str]], sources_used: Sequence[str]) -> Path:
    """enrichment.json from every comparison's result.json; EVD-0022 naming any comparison not finished."""
    from .genesets import genesets_dir
    from .jsonio import read_json, write_json

    missing = [name for name, _ in comparisons if not (out / name / "result.json").exists()]
    if missing:
        raise ProjectError(22, f"no result for {', '.join(missing)} under {out}")
    report = []
    for name, _ in comparisons:
        entry = read_json(out / name / "result.json")
        report.append({k: v for k, v in entry.items() if k not in ("program_version", "generated_by")})
    gs_manifest = read_json(genesets_dir(reference_dir, release) / "genesets.json")
    signs_path = genesets_dir(reference_dir, release) / "pathway_signs.json"
    signs = read_json(signs_path) if signs_path.exists() else {}
    return write_json(out / "enrichment.json", {
        # resolved, so a SLURM run and a local run record the same path
        "crosswalk": str(Path(crosswalk_dir).resolve()), "reference_dir": str(reference_dir),
        "ensembl_release": release,
        "gene_sets": {"indexed_by": gs_manifest.get("indexed_by"), "files": gs_manifest.get("files"),
                      "licences": gs_manifest.get("licences"), "sources_used": sorted(sources_used)},
        # which sign rule the activation z used, so an index built under another rule shows it
        "pathway_signs": ({"indexed_by": signs.get("indexed_by"), "rule": signs.get("rule")} if signs else None),
        "settings": {"alpha": settings.alpha, "min_size": settings.min_size, "max_size": settings.max_size,
                     "permutations": settings.permutations, "seed": settings.seed,
                     "max_permutations": settings.max_permutations,
                     "min_exceedances": settings.min_exceedances,
                     "ora": "one-sided hypergeometric; universe = human genes tested in the comparison; "
                            "BH within comparison x direction x source; a burden test with "
                            f"{NULL_SELECTION_ROLE}: Wallenius noncentral hypergeometric upper tail, odds = mean "
                            "selection probability of the set's genes / of the rest (goseq's approximation); "
                            "without it: hypergeometric, unadjusted for detectability (column null_model)",
                     "burden": f"no direction-of-effect statistics ({BURDEN_NO_DIRECTION}); no GSEA",
                     "gsea": "preranked, weight 1, gene-set permutation with one null per set size; "
                             "adaptive: a size gets 10x more permutations (blocks of 10,000, fixed seeds) "
                             "while any set of that size is beaten by fewer than min_exceedances nulls, "
                             "up to max_permutations; BH within comparison x source"},
        "comparisons": report,
    })


def _gene_sets(reference_dir: Path, release: int, sources, log):
    from .genesets import load_gene_sets

    info, members = load_gene_sets(reference_dir, release, sources)
    used = sorted({str(i["source"]) for i in info.values()})
    log(f"gene sets: {len(info)} from {', '.join(used)}")
    from .genesets import genesets_dir
    from .pathway_signs import stale_signs_note

    note = stale_signs_note(genesets_dir(reference_dir, release))
    if note:
        log(note)
    return info, members, used


def enrich(crosswalk_dir: Path, reference_dir: Path, *, out: Optional[Path] = None, release: int = 116,
           settings: Optional[Settings] = None, sources: Optional[Sequence[str]] = None,
           threads: Optional[int] = None, log=print) -> Path:
    """ORA and preranked GSEA for every comparison, on this machine. Returns the manifest path."""
    from .parallel import resolve_threads

    settings = settings or Settings()
    budget = resolve_threads(threads)
    crosswalk_dir = Path(crosswalk_dir)
    comparisons = _comparisons(crosswalk_dir)
    out = Path(out) if out else crosswalk_dir / "enrichment"
    info, members, used = _gene_sets(reference_dir, release, sources, log)
    for name, kind in comparisons:
        enrich_comparison(name, kind, crosswalk_dir, info, members, out, settings, budget, log)
    return finalize(crosswalk_dir, reference_dir, out, release, settings, comparisons, used)


# ---------------------------------------------------------------------------
# the SLURM backend
# ---------------------------------------------------------------------------

PLAN_NAME = "plan.json"


def plan_slurm(crosswalk_dir: Path, reference_dir: Path, *, options, out: Optional[Path] = None,
               release: int = 116, settings: Optional[Settings] = None, sources: Optional[Sequence[str]] = None,
               submit=None, extra: Optional[Dict[str, object]] = None, log=print) -> Dict[str, str]:
    """Write the plan and two sbatch scripts, submit an array (one task per comparison) and a finalize job.

    Returns {"array": job id, "finalize": job id, "plan": path}. Nothing is
    computed here: the login node only checks the inputs exist.
    """
    import dataclasses
    import shlex
    import sys

    from .genesets import genesets_dir
    from .identity import NAME
    from .jsonio import write_json
    from .slurm import script_text, submit as sbatch_submit
    from .version import __version__

    submit = submit or sbatch_submit
    settings = settings or Settings()
    crosswalk_dir = Path(crosswalk_dir).resolve()
    reference_dir = Path(reference_dir).resolve()
    comparisons = _comparisons(crosswalk_dir)
    if not (genesets_dir(reference_dir, release) / "gene_sets.tsv").exists():
        raise ProjectError(19, f"no gene-set index at {genesets_dir(reference_dir, release)}")
    out = (Path(out) if out else crosswalk_dir / "enrichment").resolve()
    work = out / "_slurm"
    (work / "logs").mkdir(parents=True, exist_ok=True)
    plan_path = write_json(work / PLAN_NAME, {
        "planned_by": __version__, "crosswalk_dir": str(crosswalk_dir), "reference_dir": str(reference_dir),
        "out": str(out), "release": release, "settings": dataclasses.asdict(settings),
        "sources": list(sources) if sources else None,
        "comparisons": [{"name": n, "kind": k} for n, k in comparisons],
        **(extra or {}),
    })
    command = f"{shlex.quote(sys.executable)} -m {NAME}.cli"
    n = len(comparisons)
    task = script_text(f"{NAME}-enrich", options, work / "logs" / "task.%A_%a.log",
                       f"{command} enrich-task --plan {shlex.quote(str(plan_path))} --index \"$SLURM_ARRAY_TASK_ID\"",
                       array=f"0-{n - 1}")
    final = script_text(f"{NAME}-enrich-finalize", options.small(), work / "logs" / "finalize.%j.log",
                        f"{command} enrich-finalize --plan {shlex.quote(str(plan_path))}")
    (work / "task.sbatch").write_text(task)
    (work / "finalize.sbatch").write_text(final)
    array_id = submit(work / "task.sbatch")
    final_id = submit(work / "finalize.sbatch", dependency=f"afterok:{array_id}")
    log(f"submitted {n} comparison task(s) as array job {array_id} ({options.cpus_per_task} CPUs each), "
        f"finalize job {final_id} after it")
    return {"array": array_id, "finalize": final_id, "plan": str(plan_path), "logs": str(work / "logs")}


def _read_plan(plan_path: Path):
    from .jsonio import read_json
    from .version import __version__

    plan = read_json(Path(plan_path))
    if plan["planned_by"] != __version__:
        raise ProjectError(22, f"this node runs {__version__} but the plan was written by {plan['planned_by']}: "
                               "install the same version on every node")
    settings = Settings(**plan["settings"])
    comparisons = [(c["name"], c["kind"]) for c in plan["comparisons"]]
    return plan, settings, comparisons


def run_slurm_task(plan_path: Path, index: int, *, threads: Optional[int] = None, log=print) -> Dict[str, object]:
    """One array task: the ``index``-th comparison of the plan, with the CPUs SLURM allocated."""
    from .parallel import resolve_threads

    plan, settings, comparisons = _read_plan(plan_path)
    if not 0 <= index < len(comparisons):
        raise ProjectError(22, f"task index {index} outside 0..{len(comparisons) - 1}")
    name, kind = comparisons[index]
    info, members, _ = _gene_sets(Path(plan["reference_dir"]), plan["release"], plan["sources"], log)
    return enrich_comparison(name, kind, Path(plan["crosswalk_dir"]), info, members, Path(plan["out"]), settings,
                             resolve_threads(threads), log)


def run_slurm_finalize(plan_path: Path, log=print) -> Path:
    from .genesets import load_gene_sets

    plan, settings, comparisons = _read_plan(plan_path)
    info, _ = load_gene_sets(Path(plan["reference_dir"]), plan["release"], plan["sources"])
    return finalize(Path(plan["crosswalk_dir"]), Path(plan["reference_dir"]), Path(plan["out"]), plan["release"],
                    settings, comparisons, sorted({str(i["source"]) for i in info.values()}))


def _run_null_jobs(jobs, weights, nulls, budget: int, log, label: str) -> None:
    """Run (comparison, size, permutations, tag) jobs in processes, balanced by cost, merged in job order."""
    from .parallel import map_ordered

    if not jobs:
        return
    tasks_n = max(1, min(budget, len(jobs)))
    bins: List[List[Tuple[str, int, int, str]]] = [[] for _ in range(tasks_n)]
    load = [0] * tasks_n
    for job in sorted(jobs, key=lambda j: (-j[1] * j[2], j[0], j[3])):
        i = load.index(min(load))
        bins[i].append(job)
        load[i] += job[1] * job[2]
    tasks = []
    for b in bins:
        by_comp: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        for comp, size, n, tag in b:
            by_comp[comp].append((size, n, tag))
        tasks.extend((comp, weights[comp], items) for comp, items in sorted(by_comp.items()))
    log(f"{label}: {len(jobs)} null job(s), {sum(j[2] for j in jobs)} permutations, "
        f"{len(tasks)} task(s) on up to {min(budget, len(tasks))} process(es)")
    for (comp, _, _), results in zip(tasks, map_ordered(_null_task, tasks, threads=budget, kind="process")):
        for size, pos, neg in results:
            nulls[comp][size].add(pos, neg, len(pos) + len(neg))


def _count_q(rows: Sequence[Dict[str, object]]) -> Dict[str, int]:
    counts: Dict[str, int] = defaultdict(int)
    for r in rows:
        if float(r["q"]) < 0.05:
            key = f"{r['source']}/{r['direction']}" if "direction" in r else str(r["source"])
            counts[key] += 1
    return dict(sorted(counts.items()))
