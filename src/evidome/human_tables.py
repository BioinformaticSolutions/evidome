"""Per-comparison human-gene tables: the crosswalk joined to the statistics.

Enrichment and regulator analyses run on HUMAN gene sets, so every comparison
in the input needs one value per human gene. This module produces that, and
records how each value was chosen:

- Every human gene the crosswalk reached gets one row per comparison -- tested
  or not (classify, don't exclude).
- A human gene reached by one input gene takes that gene's statistics.
- A human gene reached by several input genes (many-to-one) takes the winner
  of ``collisions.resolve`` -- the majority vote of every metric, among the
  candidates tested in THAT comparison. The winner may differ
  between comparisons. Every such ballot is written out, one row per
  candidate.
- Compara precedence: when any tested candidate is a Compara pair,
  only Compara candidates vote; NCBI and symbol candidates stay in the ballot,
  marked excluded. Without it, a weaker-source
  gene can win about half of the mixed collisions -- mostly retrocopies of
  ribosomal and hnRNP genes standing in for the gene Compara had chosen.
- An input gene with several human targets (one-to-many) contributes its
  statistics to each of them; ``source_human_targets`` says so.
- Human input maps to itself.

The per-comparison summary answers the question an agent asks first: of the
genes that changed, how many reached a human gene -- and of those that did not,
what are they (by biotype and mapping status)? It also counts significant input
genes that mapped but were outvoted in every collision they took part in: their
signal is not represented in the human table, and a reader should know.

"Significant" here is raw p < ``alpha`` (default 0.05), the convention of the
input's producer; it is a counting threshold for the summary only. Every row
carries its own p and q, and nothing is filtered by it.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .collisions import CONTRAST_METRICS, FTEST_METRICS, resolve
from .wide_stats import CONTRAST_ROLES, FTEST_ROLES, WideStatistics, is_directional, p_role_of

DEFAULT_ALPHA = 0.05

HUMAN_TABLE_IDENTITY = (
    "human_gene_id", "human_symbol", "human_biotype", "tested", "significant",
    "source_gene_id", "source_symbol", "source_biotype",
    "homology_source", "homology_type", "is_high_confidence", "identity_pct", "human_identity_pct",
    "source_human_targets", "n_candidates", "n_tested_candidates", "candidates",
    "decided_by", "tie",
)
BALLOT_COLUMNS = (
    "human_gene_id", "human_symbol", "candidate_gene_id", "candidate_symbol", "homology_source", "tested",
    "votes", "metrics_won", "winner", "tie", "decided_by", "excluded",
)
COMPARA_SOURCES = frozenset({"protein", "ncrna"})
EXCLUDED_BY_COMPARA = "Compara candidate present"


def stat_roles(kind: str, available: Sequence[str]) -> List[str]:
    """The comparison's statistic columns, in the input contract's order."""
    order = CONTRAST_ROLES if is_directional(kind) else FTEST_ROLES
    return [r for r in order if r in available]


def _fmt(value: Optional[float]) -> str:
    return "" if value is None else repr(float(value))


@dataclass
class Job:
    """Everything one comparison needs; plain data, so it pickles to a worker process."""
    name: str
    kind: str
    control: Optional[str]
    treatment: Optional[str]
    roles: List[str]
    values: Dict[str, List[Optional[float]]]  # role -> per-input-gene values
    gene_ids: List[str]
    genes: Dict[str, Dict[str, object]]  # input gene -> symbol, biotype, status, targets
    pairs: Dict[Tuple[str, str], Dict[str, str]]  # (input gene, human gene) -> ortholog evidence
    humans: Dict[str, Tuple[str, str, List[str]]]  # human gene -> (symbol, biotype, sorted candidates)
    alpha: float
    direction_from: Optional[str] = None


@dataclass
class Result:
    name: str
    human_rows: List[Dict[str, object]]
    ballot_rows: List[Dict[str, object]]
    summary: Dict[str, object]


def prepare_jobs(table: WideStatistics, crosswalk_rows: Sequence[Dict[str, object]], *,
                 alpha: float = DEFAULT_ALPHA) -> List[Job]:
    genes: Dict[str, Dict[str, object]] = {}
    pairs: Dict[Tuple[str, str], Dict[str, str]] = {}
    candidates: Dict[str, set] = defaultdict(set)
    human_info: Dict[str, Tuple[str, str]] = {}
    for r in crosswalk_rows:
        gid = str(r["input_gene_id"])
        genes.setdefault(gid, {
            "symbol": str(r["input_symbol"] or ""), "biotype": str(r["input_biotype"] or ""),
            "status": str(r["mapping_status"]), "targets": int(r["human_targets_of_input_gene"] or 0),
        })
        hid = str(r["human_gene_id"] or "")
        if not hid:
            continue
        candidates[hid].add(gid)
        human_info[hid] = (str(r["human_symbol"] or ""), str(r["human_biotype"] or ""))
        pairs[(gid, hid)] = {k: str(r[k] or "") for k in
                             ("homology_source", "homology_type", "is_high_confidence", "identity_pct",
                              "human_identity_pct")}
    humans = {h: (human_info[h][0], human_info[h][1], sorted(candidates[h])) for h in sorted(candidates)}
    jobs = []
    for comp in table.comparisons:
        roles = stat_roles(comp.kind, comp.columns)
        jobs.append(Job(comp.name, comp.kind, comp.control, comp.treatment, roles,
                        {role: table.column(comp, role) for role in roles},
                        list(table.gene_ids), genes, pairs, humans, alpha, comp.direction_from))
    return jobs


def run_job(job: Job) -> Result:
    """One comparison's human table, ballots and summary. Module-level: runs in a worker process."""
    p_role = p_role_of(job.kind)
    metric_roles = [m[0] for m in (CONTRAST_METRICS if is_directional(job.kind) else FTEST_METRICS)]
    from .enrich import fold_change_convention

    # The one reading of the Fold_change convention (enrich.fold_change_convention): a ratio is voted on
    # |log ratio|; signed times and log fold changes on |value|. A declared Log2_fold_change is log2 and is never
    # judged from its values: all of them between 0 and 1 would read as a ratio.
    fcs = [v for v in job.values.get("Fold_change", ()) if v is not None]
    fc_ratio = is_directional(job.kind) and "Fold_change" in job.values and fold_change_convention(fcs) == "ratio"
    index = {g: i for i, g in enumerate(job.gene_ids)}

    def stats(gid: str) -> Dict[str, Optional[float]]:
        i = index.get(gid)
        return {role: (None if i is None else job.values[role][i]) for role in job.roles}

    human_rows: List[Dict[str, object]] = []
    ballot_rows: List[Dict[str, object]] = []
    winners: set = set()
    decided: Dict[str, int] = defaultdict(int)
    human_tested = human_significant = 0

    excluded_significant: set = set()
    for hid, (hsym, hbio, cands) in job.humans.items():
        per = {g: stats(g) for g in cands}
        source = {g: job.pairs.get((g, hid), {}).get("homology_source", "") for g in cands}
        tested = [g for g in cands if per[g].get(p_role) is not None]
        excluded: set = set()
        if any(source[g] in COMPARA_SOURCES for g in tested):
            excluded = {g for g in cands if source[g] not in COMPARA_SOURCES and source[g]}
        voters = [g for g in cands if g not in excluded]
        ballot = resolve(hid, job.name, job.kind, {g: {r: per[g].get(r) for r in metric_roles} for g in voters},
                         fc_ratio=fc_ratio)
        if excluded & set(tested) and ballot.decided_by == "only tested candidate":
            ballot.decided_by = "Compara precedence"
        for g in excluded:
            p = per[g].get(p_role)
            if p is not None and p < job.alpha:
                excluded_significant.add(g)
        row: Dict[str, object] = {c: "" for c in HUMAN_TABLE_IDENTITY}
        row.update({"human_gene_id": hid, "human_symbol": hsym, "human_biotype": hbio,
                    "n_candidates": len(cands), "n_tested_candidates": len(tested),
                    "candidates": ";".join(cands), "decided_by": ballot.decided_by,
                    "tie": "yes" if ballot.tie else "no"})
        if len(cands) == 1 and ballot.decided_by == "only tested candidate":
            row["decided_by"] = "single candidate"
        winner = ballot.winner
        if winner is None:
            row.update({"tested": "no", "significant": "no"})
            row.update({role: "" for role in job.roles})
        else:
            winners.add(winner)
            decided[str(row["decided_by"])] += 1
            human_tested += 1
            p = per[winner].get(p_role)
            sig = p is not None and p < job.alpha
            human_significant += sig
            info = job.genes.get(winner, {})
            evidence = job.pairs.get((winner, hid), {})
            row.update({"tested": "yes", "significant": "yes" if sig else "no",
                        "source_gene_id": winner, "source_symbol": info.get("symbol", ""),
                        "source_biotype": info.get("biotype", ""),
                        "source_human_targets": info.get("targets", ""), **evidence})
            row.update({role: _fmt(per[winner].get(role)) for role in job.roles})
        human_rows.append(row)

        if len(cands) > 1:
            won: Dict[str, List[str]] = defaultdict(list)
            for metric, pick in ballot.picks.items():
                won[pick].append(metric)
            for g in cands:
                b: Dict[str, object] = {
                    "human_gene_id": hid, "human_symbol": hsym, "candidate_gene_id": g,
                    "candidate_symbol": job.genes.get(g, {}).get("symbol", ""),
                    "homology_source": source[g],
                    "tested": "yes" if g in tested else "no",
                    "votes": ballot.votes.get(g, 0) if g in ballot.tested else "",
                    "metrics_won": ";".join(won.get(g, [])),
                    "winner": "yes" if g == winner else "no",
                    "tie": "yes" if ballot.tie else "no", "decided_by": ballot.decided_by,
                    "excluded": EXCLUDED_BY_COMPARA if g in excluded else "",
                }
                b.update({role: _fmt(per[g].get(role)) for role in job.roles})
                ballot_rows.append(b)

    # The summary, from the input side.
    p_values = job.values[p_role]
    tested_ids = [g for g, p in zip(job.gene_ids, p_values) if p is not None]
    sig_ids = [g for g, p in zip(job.gene_ids, p_values) if p is not None and p < job.alpha]
    mapped = {g for (g, _h) in job.pairs}

    def unmapped_breakdown(ids: Sequence[str]) -> Tuple[Dict[str, int], Dict[str, int]]:
        by_status: Dict[str, int] = defaultdict(int)
        by_biotype: Dict[str, int] = defaultdict(int)
        for g in ids:
            info = job.genes.get(g, {})
            by_status[str(info.get("status", "")) or "(not in crosswalk)"] += 1
            by_biotype[str(info.get("biotype", "")) or "(not annotated)"] += 1
        order = lambda d: dict(sorted(d.items(), key=lambda kv: (-kv[1], kv[0])))  # noqa: E731
        return order(by_status), order(by_biotype)

    sig_unmapped = [g for g in sig_ids if g not in mapped]
    sources: Dict[str, str] = {}
    for (g, _h), evidence in sorted(job.pairs.items()):
        sources.setdefault(g, evidence.get("homology_source", "") or "(identity)")
    sig_by_source: Dict[str, int] = defaultdict(int)
    for g in sig_ids:
        if g in mapped:
            sig_by_source[sources[g]] += 1
    sig_status, sig_biotype = unmapped_breakdown(sig_unmapped)
    collisions_tested = sum(1 for r in human_rows if int(r["n_tested_candidates"]) > 1)
    summary = {
        "comparison": job.name, "kind": job.kind, "control": job.control, "treatment": job.treatment,
        "direction_from": job.direction_from, "alpha": job.alpha, "significance": f"raw {p_role} < {job.alpha}",
        "input_genes_tested": len(tested_ids),
        "input_genes_tested_and_mapped": sum(g in mapped for g in tested_ids),
        "input_genes_significant": len(sig_ids),
        "significant_mapped": len(sig_ids) - len(sig_unmapped),
        "significant_mapped_by_source": dict(sorted(sig_by_source.items())),
        "significant_mapped_but_outvoted": sum(g in mapped and g not in winners and g not in excluded_significant
                                               for g in sig_ids),
        "significant_excluded_by_compara_precedence": sum(g not in winners for g in excluded_significant),
        "significant_unmapped": len(sig_unmapped),
        "significant_unmapped_by_status": sig_status,
        "significant_unmapped_by_biotype": sig_biotype,
        "human_genes_reached": len(job.humans),
        "human_genes_tested": human_tested,
        "human_genes_significant": human_significant,
        "human_genes_with_tested_collisions": collisions_tested,
        "human_genes_by_decision": dict(sorted(decided.items())),
    }
    return Result(job.name, human_rows, ballot_rows, summary)


def build_human_tables(table: WideStatistics, crosswalk_rows: Sequence[Dict[str, object]], *,
                       alpha: float = DEFAULT_ALPHA, threads: Optional[int] = None) -> List[Result]:
    """Every comparison, in parallel (one process per comparison), results in input order."""
    from .parallel import map_ordered, resolve_threads

    jobs = prepare_jobs(table, crosswalk_rows, alpha=alpha)
    return map_ordered(run_job, jobs, threads=resolve_threads(threads), kind="process")


def _write_csv(path: Path, columns: Sequence[str], rows: Sequence[Dict[str, object]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)
    return path


def write_human_tables(results: Sequence[Result], table: WideStatistics, out: Path) -> List[Path]:
    """``comparisons/<name>.human_genes.csv`` and ``collisions/<name>.ballots.csv``."""
    kinds = {c.name: c for c in table.comparisons}
    written = []
    for result in results:
        comp = kinds[result.name]
        roles = stat_roles(comp.kind, comp.columns)
        written.append(_write_csv(Path(out) / "comparisons" / f"{result.name}.human_genes.csv",
                                  HUMAN_TABLE_IDENTITY + tuple(roles), result.human_rows))
        written.append(_write_csv(Path(out) / "collisions" / f"{result.name}.ballots.csv",
                                  BALLOT_COLUMNS + tuple(roles), result.ballot_rows))
    return written


def summary_lines(results: Sequence[Result]) -> List[str]:
    """The per-comparison table printed by ``crosswalk``."""
    if not results:
        return []
    alpha = results[0].summary["alpha"]
    lines = [f"per comparison (significant = raw p < {alpha})",
             f"  {'comparison':<34} {'tested':>7} {'signif':>7} {'->human':>7} {'via NCBI':>8} "
             f"{'via sym':>7} {'outvoted':>8} {'excluded':>8} {'unmapped':>8} {'human signif':>12} {'collisions':>10}"]
    for r in results:
        s = r.summary
        via = s["significant_mapped_by_source"]
        lines.append(f"  {s['comparison']:<34} {s['input_genes_tested']:>7} {s['input_genes_significant']:>7} "
                     f"{s['significant_mapped']:>7} {via.get('ncbi', 0):>8} {via.get('symbol', 0):>7} "
                     f"{s['significant_mapped_but_outvoted']:>8} "
                     f"{s['significant_excluded_by_compara_precedence']:>8} "
                     f"{s['significant_unmapped']:>8} {s['human_genes_significant']:>12} "
                     f"{s['human_genes_with_tested_collisions']:>10}")
    for r in results:
        s = r.summary
        if s["significant_unmapped"]:
            parts = ", ".join(f"{k} {v}" for k, v in list(s["significant_unmapped_by_biotype"].items())[:5])
            lines.append(f"  {s['comparison']}: significant but unmapped -- {parts}")
    return lines
