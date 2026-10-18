"""Upstream-regulator evidence.

Reference data (``reference regulators``), from the OmniPath web service with
``license=commercial`` -- OmniPath keeps only records supported by resources
whose licences allow for-profit use -- human only:

    collectri            CollecTRI: signed transcription factor -> target gene
                         regulons (Mueller-Dott et al. 2023, NAR 51:10934)
    omnipath             OmniPath's literature-curated, directed, signed
                         protein -> protein signalling network

Header verified against the live service (2026-09-27): source, target,
source_genesymbol, target_genesymbol, is_directed, is_stimulation,
is_inhibition, consensus_direction, consensus_stimulation,
consensus_inhibition, then the requested fields.

Two regulator networks are indexed into ``ensembl-<release>/regulators/``:

    collectri            each transcription factor (or TF complex) and its
                         direct target genes, sign = consensus stimulation
                         (+1) or inhibition (-1); neither or both = 0 (the
                         edge counts for overlap, not for direction)
    signalling_via_tf    each signalling protein X with a signed edge to a
                         CollecTRI transcription factor T: X's targets are T's
                         targets, sign = sign(X->T) x sign(T->gene), summed
                         over every T; a gene whose paths cancel is dropped.
                         INDIRECT evidence, labelled as such everywhere.

Targets are mapped to human Ensembl genes of the release by exact, unambiguous
symbol (the service's gene symbols); unmapped and ambiguous ones are counted.

Analysis (``regulators <crosswalk dir>``), per comparison and regulator,
among the human genes tested in that comparison (the universe):

    overlap_p      exact hypergeometric: are the regulator's targets
                   over-represented among significant genes (raw p < alpha)?
                   Every comparison.
    activation_z   contrasts: over significant targets with a signed edge,
                   consistent = sign(edge) x direction(gene) = +1;
                   z = (consistent - inconsistent) / sqrt(their sum).
                   z_p: exact two-sided binomial
                   p of that split.
    activation_z_bias_corrected
                   contrasts: the correction for a dataset
                   whose significant genes lean one way (Kraemer et al. 2014,
                   Bioinformatics 30:523): mu = mu_data x mu_regulator, where
                   mu_data = (up - down) / (up + down) over the comparison's
                   significant genes and mu_regulator = (activating -
                   inhibiting) / (their sum) over the regulator's signed
                   targets among tested genes; z_b = sum(x_i - mu) /
                   sqrt(N (1 - mu^2)). bias_term reports mu.
    activity_t     contrasts: univariate linear model over ALL tested genes,
                   rank metric (the GSEA metric) ~ edge sign (0 for
                   non-targets); t of the slope and its two-sided p (CollecTRI
                   authors' recommended method, decoupler "ulm"). No
                   significance threshold involved.

BH q-values within (comparison, network) for overlap_p and activity_p.

Burden tests (kind 'burden'): overlap_p only. activation_z, z_p,
bias_term, activation_z_bias_corrected, activity_t/p/q, the consistent /
inconsistent counts and regulator_direction are left blank -- a burden test
measures damage counts, not expression direction (regulators.json says so per
comparison). overlap_p uses the same null as the comparison's ORA (enrich.py):
Wallenius, weighted by each gene's null selection probability, when the input
gives one; otherwise hypergeometric, unadjusted for detectability. Column
overlap_null names it on every row.

Signalling regulators also carry via_tfs -- every transcription
factor they reach -- and driving_tfs -- those of them with significant targets
in the comparison, as TF(count), most first. Many signalling proteins act
through the same few TFs, so their results are largely the same evidence;
driving_tfs lets a reader collapse them onto the TFs behind them.
Regulators are tested with >= min_targets targets in the universe and, for
the indirect network, <= max_targets (hub kinases reach most of the genome).
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .errors import ProjectError

OMNIPATH_BASE = "https://omnipathdb.org"
NETWORKS = ("collectri", "signalling_via_tf")
REGULON_COLUMNS = ("network", "regulator", "regulator_kind", "regulator_human_gene_id", "target_human_gene_id",
                   "target_symbol", "sign", "evidence")
REGULATOR_COLUMNS = ("comparison", "network", "regulator", "regulator_kind", "targets_in_universe",
                     "targets_signed", "significant_targets", "overlap_p", "overlap_q", "n_consistent",
                     "n_inconsistent", "activation_z", "z_p", "bias_term", "activation_z_bias_corrected",
                   "activity_t", "activity_p", "activity_q", "via_tfs", "driving_tfs",
                     "regulator_in_data", "regulator_direction", "regulator_p", "rank_metric",
                     "consistent_targets", "inconsistent_targets", "overlap_null")
DEFAULT_MIN_TARGETS = 5
DEFAULT_MAX_TARGETS = 1000
_NEEDED = ("source", "target", "source_genesymbol", "target_genesymbol", "consensus_stimulation",
           "consensus_inhibition")


def regulators_dir(reference_dir: Path, release: int) -> Path:
    from .reference import release_dir

    return release_dir(reference_dir, release) / "regulators"


def omnipath_urls(base: str = OMNIPATH_BASE) -> Dict[str, str]:
    common = "organisms=9606&genesymbols=yes&license=commercial&format=tsv"
    return {
        "collectri": f"{base}/interactions?datasets=collectri&{common}&fields=sources,curation_effort",
        "omnipath": f"{base}/interactions?datasets=omnipath&{common}&directed=yes&fields=sources,type",
    }


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def _betacf(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta function (Lentz), as Numerical Recipes 6.4."""
    tiny, qab, qap, qam = 1e-300, a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h


def betainc(a: float, b: float, x: float) -> float:
    """Regularised incomplete beta I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_two_sided_p(t: float, df: float) -> float:
    if df <= 0 or math.isnan(t):
        return float("nan")
    return betainc(df / 2.0, 0.5, df / (df + t * t))


def binom_two_sided_p(k: int, n: int) -> float:
    """Exact two-sided binomial p for k successes of n at 1/2: sum of outcomes no likelier than k."""
    if n == 0:
        return 1.0
    logp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2) for i in range(n + 1)]
    cut = logp[k] + 1e-9
    return min(1.0, sum(math.exp(v) for v in logp if v <= cut))


class Ulm:
    """Univariate linear model y ~ x over a fixed universe, with sparse x (non-targets are 0)."""

    def __init__(self, y: Dict[str, float]):
        self.y = y
        self.n = len(y)
        self.sy = sum(y.values())
        self.syy = sum(v * v for v in y.values())

    def t(self, weights: Dict[str, int]) -> Tuple[Optional[float], Optional[float]]:
        w = {g: s for g, s in weights.items() if g in self.y and s}
        n = self.n
        if n < 3 or not w:
            return None, None
        sx = float(sum(w.values()))
        sxx = float(sum(s * s for s in w.values()))
        sxy = sum(s * self.y[g] for g, s in w.items())
        vx = sxx - sx * sx / n
        if vx <= 0:
            return None, None
        cxy = sxy - sx * self.sy / n
        vy = self.syy - self.sy * self.sy / n
        slope = cxy / vx
        sse = max(vy - slope * cxy, 0.0)
        if sse == 0.0:
            return None, None
        se = math.sqrt(sse / (n - 2) / vx)
        t = slope / se
        return t, t_two_sided_p(t, n - 2)


# ---------------------------------------------------------------------------
# reference data
# ---------------------------------------------------------------------------

def _sign(row: Dict[str, str]) -> int:
    stim = row["consensus_stimulation"].strip().lower() in ("true", "1")
    inhib = row["consensus_inhibition"].strip().lower() in ("true", "1")
    return 1 if stim and not inhib else (-1 if inhib and not stim else 0)


def _read_omnipath(path: Path) -> List[Dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in _NEEDED if c not in (reader.fieldnames or [])]
        if missing:
            raise ProjectError(18, f"{Path(path).name}: no column {', '.join(missing)}; header is {reader.fieldnames}")
        rows = list(reader)
    if not rows:
        raise ProjectError(18, f"{Path(path).name}: no interactions")
    return rows


def _symbol_index(human_genes: Path) -> Tuple[Dict[str, str], Set[str]]:
    by_symbol: Dict[str, List[str]] = defaultdict(list)
    with Path(human_genes).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["human_symbol"]:
                by_symbol[row["human_symbol"].upper()].append(row["human_gene_id"])
    unique = {s: ids[0] for s, ids in by_symbol.items() if len(ids) == 1}
    ambiguous = {s for s, ids in by_symbol.items() if len(ids) > 1}
    return unique, ambiguous


def build_regulons(collectri_rows, omnipath_rows, symbols: Dict[str, str], ambiguous: Set[str]):
    """-> (rows for regulons.tsv, counts). Pure function of the two downloads and the symbol index."""
    counts: Dict[str, int] = defaultdict(int)
    direct: Dict[str, Dict[str, int]] = defaultdict(dict)       # TF name -> {target gene: sign}
    evidence: Dict[Tuple[str, str], str] = {}
    kind: Dict[str, str] = {}
    tf_by_uniprot: Dict[str, str] = {}
    for r in collectri_rows:
        tf, target = r["source_genesymbol"], r["target_genesymbol"].upper()
        gene = symbols.get(target)
        counts["collectri_edges"] += 1
        if gene is None:
            counts["collectri_targets_ambiguous" if target in ambiguous else "collectri_targets_unmapped"] += 1
            continue
        s = _sign(r)
        previous = direct[tf].get(gene)
        direct[tf][gene] = s if previous is None or previous == s else 0     # conflicting records: unsigned
        evidence[(tf, gene)] = r.get("sources", "")
        kind[tf] = "complex" if r["source"].startswith("COMPLEX:") else "tf"
        if not r["source"].startswith("COMPLEX:"):
            tf_by_uniprot[r["source"]] = tf
    rows = []
    for tf in sorted(direct):
        reg_gene = symbols.get(tf.upper(), "") if kind[tf] == "tf" else ""
        for gene, s in sorted(direct[tf].items()):
            rows.append(("collectri", tf, kind[tf], reg_gene, gene, "", s, evidence[(tf, gene)]))
    counts["collectri_regulators"] = len(direct)

    # indirect: X -(s1)-> TF -(s2)-> gene; sum of s1*s2 over TFs, sign kept, zero dropped
    via: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    tfs_used: Dict[str, Set[str]] = defaultdict(set)
    for r in omnipath_rows:
        tf = tf_by_uniprot.get(r["target"])
        s1 = _sign(r)
        x = r["source_genesymbol"]
        if tf is None or s1 == 0 or r["source"].startswith("COMPLEX:") or x == tf or not x:
            continue
        tfs_used[x].add(tf)
        for gene, s2 in direct[tf].items():
            if s2:
                via[x][gene] += s1 * s2
    for x in sorted(via):
        reg_gene = symbols.get(x.upper(), "")
        path = ";".join(sorted(tfs_used[x]))
        for gene, total in sorted(via[x].items()):
            if total:
                rows.append(("signalling_via_tf", x, "signalling", reg_gene, gene, "", 1 if total > 0 else -1,
                             f"via {path}"))
    # regulators with at least one target row (ones whose paths all cancel are not counted)
    counts["signalling_regulators"] = sum(1 for x in via if any(via[x].values()))
    counts["regulon_rows"] = len(rows)
    return rows, dict(counts)


def fetch_regulators(reference_dir: Path, *, release: int = 116, refresh: bool = False,
                     urls: Optional[Dict[str, str]] = None, threads: Optional[int] = None, log=print) -> Path:
    from datetime import datetime, timezone

    from .jsonio import write_json
    from .parallel import MAX_DOWNLOADS, map_ordered, resolve_threads
    from .reference import _download_ncbi, ensure_writable, human_genes_path, load_manifest, release_dir
    from .version import __version__

    load_manifest(reference_dir, release)
    human = human_genes_path(reference_dir, release)       # before downloading
    out = ensure_writable(regulators_dir(reference_dir, release))
    wanted = urls or omnipath_urls()
    tasks = [(wanted[k], str(out / "raw" / f"omnipath_{k}.tsv"), refresh) for k in ("collectri", "omnipath")]
    got = dict(zip(("collectri", "omnipath"),
                   map_ordered(_download_ncbi, tasks, threads=min(resolve_threads(threads), MAX_DOWNLOADS),
                               kind="thread")))
    for k, d in got.items():
        log(f"  {'reused' if d['reused'] else 'fetched'} OmniPath {k} ({d['bytes'] / 1e6:.1f} MB, md5 recorded)")
    symbols, ambiguous = _symbol_index(human)
    rows, counts = build_regulons(_read_omnipath(Path(str(got["collectri"]["path"]))),
                                  _read_omnipath(Path(str(got["omnipath"]["path"]))), symbols, ambiguous)
    with (out / "regulons.tsv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(REGULON_COLUMNS)
        w.writerows(rows)
    path = write_json(out / "regulators.json", {
        "ensembl_release": release, "indexed_by": __version__,
        "indexed": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": {k: {f: d[f] for f in ("url", "checksum", "bytes", "downloaded")} for k, d in got.items()},
        "licence": "OmniPath web service with license=commercial: only records supported by resources whose "
                   "licences allow for-profit use (OmniPath's own per-resource licence tracking).",
        "attribution": "OmniPath (Tuerei et al.; omnipathdb.org); CollecTRI (Mueller-Dott et al. 2023, "
                       "Nucleic Acids Research 51:10934).",
        "counts": counts,
        "rules": {"sign": "consensus_stimulation only = +1, consensus_inhibition only = -1, else 0",
                  "signalling_via_tf": "X->TF (signed, OmniPath) then TF->gene (CollecTRI): sum of sign products; "
                                       "zero sums dropped",
                  "targets": "exact unambiguous human symbol in the Ensembl release; others counted"},
    })
    log(f"  collectri          {counts.get('collectri_regulators', 0):>6} regulators  "
        f"({counts.get('collectri_targets_unmapped', 0)} target records unmapped, "
        f"{counts.get('collectri_targets_ambiguous', 0)} ambiguous)")
    log(f"  signalling_via_tf  {counts.get('signalling_regulators', 0):>6} regulators")
    return path


def load_regulons(reference_dir: Path, release: int = 116, networks: Optional[Iterable[str]] = None):
    """{(network, regulator): {"kind", "gene", "targets": {gene: sign}}}; EVD-0025 if not indexed."""
    path = regulators_dir(reference_dir, release) / "regulons.tsv"
    if not path.exists():
        raise ProjectError(25, f"no regulator index at {path.parent}")
    wanted = set(networks) if networks else set(NETWORKS)
    regulons: Dict[Tuple[str, str], Dict[str, object]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["network"] not in wanted:
                continue
            key = (row["network"], row["regulator"])
            entry = regulons.setdefault(key, {"kind": row["regulator_kind"], "gene": row["regulator_human_gene_id"],
                                               "targets": {}, "via": []})
            if row["network"] == "signalling_via_tf" and not entry["via"] and row["evidence"].startswith("via "):
                entry["via"] = row["evidence"][4:].split(";")
            entry["targets"][row["target_human_gene_id"]] = int(row["sign"])  # type: ignore[index]
    return regulons


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def analyse_comparison(table, regulons, *, alpha: float, min_targets: int, max_targets: int,
                       networks: Optional[Set[str]] = None) -> List[Dict[str, object]]:
    """Every regulator for one comparison (an ``enrich.ComparisonTable``)."""
    from .enrich import bh, ora_null, overlap_p

    universe = set(table.genes)
    N = len(universe)
    significant = {g for g in table.genes if table.p[g] < alpha}
    n_sig = len(significant)
    # has a direction: a contrast, or an F-test directed by its contrasts. Not a burden test: its
    # direction-of-effect statistics are not reported (enrich.BURDEN_NO_DIRECTION).
    contrast = bool(table.metric) and table.kind != "burden"
    null = ora_null(table, alpha)
    # Over ALL tested genes, as the module docstring says: a gene with no direction (fold change exactly 1) is a
    # tested gene with no signed signal, so it enters as 0 -- as in database-value -- rather than leaving the
    # regression.
    ulm = Ulm({g: table.metric.get(g, 0.0) for g in table.genes}) if contrast else None
    up = sum(1 for g in significant if table.direction.get(g, 0) > 0)
    down = sum(1 for g in significant if table.direction.get(g, 0) < 0)
    mu_data = (up - down) / (up + down) if up + down else 0.0
    tf_significant: Dict[str, int] = {}
    for (network, regulator), entry in regulons.items():
        if network == "collectri":
            tf_significant[regulator] = sum(1 for g in entry["targets"] if g in significant)  # type: ignore[union-attr]
    rows_by_network: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for (network, regulator), entry in sorted(regulons.items()):
        if networks is not None and network not in networks:
            continue
        targets = {g: s for g, s in entry["targets"].items() if g in universe}  # type: ignore[union-attr]
        K = len(targets)
        if K < min_targets or (network == "signalling_via_tf" and K > max_targets):
            continue
        sig_targets = [g for g in targets if g in significant]
        row: Dict[str, object] = {
            "comparison": table.name, "network": network, "regulator": regulator,
            "regulator_kind": entry["kind"], "targets_in_universe": K,
            "targets_signed": sum(1 for s in targets.values() if s),
            "significant_targets": len(sig_targets),
            "overlap_p": overlap_p(table, null, len(sig_targets), N, list(targets), n_sig),
            "rank_metric": table.metric_name if contrast else "",
            "overlap_null": null,
        }
        via = entry.get("via") or []
        if network == "signalling_via_tf":
            driving = sorted(((tf_significant.get(tf, 0), tf) for tf in via if tf_significant.get(tf, 0)),
                             key=lambda x: (-x[0], x[1]))
            row.update({"via_tfs": ";".join(via), "driving_tfs": ";".join(f"{tf}({n})" for n, tf in driving)})
        gene = entry["gene"]
        if gene and gene in universe:
            row.update({"regulator_in_data": "yes", "regulator_p": table.p[gene],
                        "regulator_direction": {1: "up", -1: "down", 0: ""}[table.direction.get(gene, 0)]
                        if contrast else ""})
        else:
            row["regulator_in_data"] = "no"
        if contrast:
            consistent, inconsistent = [], []
            for g in sig_targets:
                s, d = targets[g], table.direction.get(g, 0)
                if s and d:
                    (consistent if s * d > 0 else inconsistent).append(table.symbols[g])
            n = len(consistent) + len(inconsistent)
            t, tp = ulm.t(targets)  # type: ignore[union-attr]
            act = sum(1 for s in targets.values() if s > 0)
            inh = sum(1 for s in targets.values() if s < 0)
            mu = mu_data * ((act - inh) / (act + inh) if act + inh else 0.0)
            z_b = ""
            if n and abs(mu) < 1:
                z_b = round((len(consistent) - len(inconsistent) - n * mu) / math.sqrt(n * (1 - mu * mu)), 6)
            row.update({
                "n_consistent": len(consistent), "n_inconsistent": len(inconsistent),
                "activation_z": round((len(consistent) - len(inconsistent)) / math.sqrt(n), 6) if n else "",
                "z_p": binom_two_sided_p(len(consistent), n) if n else "",
                "bias_term": round(mu, 6), "activation_z_bias_corrected": z_b,
                "activity_t": round(t, 6) if t is not None else "", "activity_p": tp if tp is not None else "",
                "consistent_targets": ";".join(sorted(consistent)),
                "inconsistent_targets": ";".join(sorted(inconsistent)),
            })
        rows_by_network[network].append(row)
    rows: List[Dict[str, object]] = []
    for network in NETWORKS:
        group = rows_by_network.get(network, [])
        for row, q in zip(group, bh([float(r["overlap_p"]) for r in group])):
            row["overlap_q"] = q
        with_p = [r for r in group if r.get("activity_p", "") != ""]
        for row, q in zip(with_p, bh([float(r["activity_p"]) for r in with_p])):
            row["activity_q"] = q
        rows.extend(group)
    return rows


def run_regulators(crosswalk_dir: Path, reference_dir: Path, *, out: Optional[Path] = None, release: int = 116,
                   alpha: float = 0.05, min_targets: int = DEFAULT_MIN_TARGETS,
                   max_targets: int = DEFAULT_MAX_TARGETS, networks: Optional[Sequence[str]] = None,
                   threads: Optional[int] = None, log=print) -> Path:
    from .enrich import _comparisons, _write, read_comparison
    from .jsonio import read_json, write_json
    from .parallel import map_ordered, resolve_threads

    crosswalk_dir = Path(crosswalk_dir)
    comparisons = _comparisons(crosswalk_dir)
    out = Path(out) if out else crosswalk_dir / "regulators"
    regulons = load_regulons(reference_dir, release)       # all: signalling rows need the TF regulons
    wanted = set(networks) if networks else set(NETWORKS)
    log(f"regulators: {sum(1 for k in regulons if k[0] in wanted)} across {', '.join(sorted(wanted))}")
    tasks = [(crosswalk_dir / "comparisons" / f"{name}.human_genes.csv", name, kind, regulons, alpha, min_targets,
              max_targets, wanted) for name, kind in comparisons]
    results = map_ordered(_task, tasks, threads=resolve_threads(threads), kind="process")
    report = []
    for (name, kind), rows in zip(comparisons, results):
        path = _write(out / name / "regulators.csv", REGULATOR_COLUMNS, rows)
        counts: Dict[str, int] = defaultdict(int)
        for r in rows:
            if float(r["overlap_q"]) < 0.05:
                counts[f"{r['network']}/overlap"] += 1
            if r.get("activity_q", "") != "" and float(r["activity_q"]) < 0.05:
                counts[f"{r['network']}/activity"] += 1
        entry = {"comparison": name, "kind": kind, "regulators_tested": len(rows),
                 "file": str(path.relative_to(out)), "q_below_0.05": dict(sorted(counts.items()))}
        if kind == "burden":
            from .enrich import BURDEN_NO_DIRECTION

            entry["direction_statistics"] = BURDEN_NO_DIRECTION
            entry["overlap_null"] = rows[0]["overlap_null"] if rows else ""
        report.append(entry)
    ref_manifest = read_json(regulators_dir(reference_dir, release) / "regulators.json")
    return write_json(out / "regulators.json", {
        "crosswalk": str(Path(crosswalk_dir).resolve()), "reference_dir": str(reference_dir),
        "ensembl_release": release,
        "regulator_data": {k: ref_manifest.get(k) for k in ("indexed_by", "files", "licence", "attribution",
                                                              "counts", "rules")},
        "settings": {"alpha": alpha, "min_targets": min_targets, "max_targets_signalling_via_tf": max_targets,
                     "networks": sorted(wanted),
                     "overlap": "one-sided hypergeometric of significant genes among targets; universe = tested "
                                "human genes; BH within comparison x network; a burden test with null selection "
                                "probabilities: Wallenius, weighted by detectability, as its ORA (overlap_null)",
                     "burden": "overlap only: no activation z, bias-corrected z, activity t or regulator direction "
                               "(a burden test measures damage counts, not expression direction)",
                     "activation_z": "(consistent - inconsistent) / sqrt(n) over significant signed targets; "
                                     "z_p exact two-sided binomial",
                     "activation_z_bias_corrected": "Kraemer et al. 2014: mu = mu_data x mu_regulator; "
                                                    "sum(x - mu) / sqrt(n (1 - mu^2))",
                     "activity": "univariate linear model, rank metric ~ edge sign over all tested genes; t and "
                                 "two-sided p on n-2 df; BH within comparison x network"},
        "comparisons": report,
    })


def _task(task):
    from .enrich import read_comparison

    path, name, kind, regulons, alpha, min_targets, max_targets, networks = task
    table = read_comparison(path, name, kind)
    return analyse_comparison(table, regulons, alpha=alpha, min_targets=min_targets, max_targets=max_targets,
                              networks=networks)


def regulator_lines(manifest_path: Path) -> List[str]:
    from .jsonio import read_json

    manifest = read_json(manifest_path)
    lines = ["regulators with q < 0.05 (all tested regulators are in the files; nothing is filtered)"]
    for c in manifest["comparisons"]:
        found = ", ".join(f"{k} {v}" for k, v in c["q_below_0.05"].items()) or "none"
        lines.append(f"  {c['comparison']:<34} tested {c['regulators_tested']:>5}   {found}")
    lines.append(f"wrote {Path(manifest_path).parent}")
    return lines
