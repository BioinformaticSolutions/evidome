"""How much each database adds that the others do not.

Too many competing or redundant databases add noise. This is a utility score
for each database: how much new information each one adds.

Measured on a finished run, from the results themselves -- what matters is
what a database adds to THIS study's evidence, not how its library overlaps
in the abstract. Two families, each database compared with the others of its
family:

    pathways     reactome, go_bp, go_mf, go_cc, wikipathways, hallmark (ORA,
                 all significant genes); a finding = a gene set with q < --q;
                 its evidence = the significant genes it contains
    regulators   collectri, signalling_via_tf; a finding = a regulator with
                 overlap q < --q; its evidence = its significant signed
                 targets (consistent + inconsistent). F-tests without a
                 direction have no such list and are not counted here.

Per database, per comparison and pooled over comparisons:

    findings              significant findings
    independent_findings  after merging its own near-duplicates: findings
                          whose evidence overlaps at Jaccard >= --jaccard are
                          one (greedy, lowest q first) -- a database that
                          says the same thing ten ways counts once
    unique_findings       findings with no near-duplicate (Jaccard >=
                          --jaccard) among the OTHER databases' findings
    unique_fraction       unique_findings / findings
    evidence_genes        genes in any of its findings
    unique_genes          of those, genes in no other database's findings
    marginal_gain         unique_genes / all evidence genes of the family:
                          the share of the family's evidence that would be
                          lost if this database were dropped
    utility               unique_fraction x independence (independent /
                          findings): 1 = every finding is new and distinct,
                          0 = everything it says is said elsewhere or by
                          itself already

Nothing here changes an analysis; it describes one. Written to
``<run>/source_utility.csv`` by the ``utility`` command.
"""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from .errors import ProjectError

PATHWAY_SOURCES = ("reactome", "go_bp", "go_mf", "go_cc", "wikipathways", "hallmark")
REGULATOR_NETWORKS = ("collectri", "signalling_via_tf")
UTILITY_COLUMNS = ("family", "source", "comparison", "findings", "independent_findings", "unique_findings",
                   "unique_fraction", "evidence_genes", "unique_genes", "marginal_gain", "utility")

Finding = Tuple[float, str, FrozenSet[str]]          # (q, id, evidence genes)


def jaccard(a: FrozenSet[str], b: FrozenSet[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def _num(text: object) -> Optional[float]:
    try:
        return float(str(text))
    except ValueError:
        return None


def _genes(*fields: str) -> FrozenSet[str]:
    return frozenset(g for f in fields for g in (f or "").split(";") if g)


def pathway_findings(ora: Path, q_max: float) -> Dict[str, List[Finding]]:
    out: Dict[str, List[Finding]] = defaultdict(list)
    with Path(ora).open(newline="", encoding="utf-8") as handle:
        for r in csv.DictReader(handle):
            q = _num(r.get("q"))
            if r.get("direction") == "all" and q is not None and q < q_max:
                genes = _genes(r.get("overlap_genes", ""))
                if genes:
                    out[r["source"]].append((q, r["set_id"], genes))
    return out


def regulator_findings(path: Path, q_max: float) -> Dict[str, List[Finding]]:
    out: Dict[str, List[Finding]] = defaultdict(list)
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for r in csv.DictReader(handle):
            q = _num(r.get("overlap_q"))
            if q is not None and q < q_max:
                # A finding is defined by its q: one whose significant targets are all unsigned has no evidence
                # genes, and is then new by definition -- it is still counted.
                genes = _genes(r.get("consistent_targets", ""), r.get("inconsistent_targets", ""))
                out[r["network"]].append((q, r["regulator"], genes))
    return out


def _directed(path: Path) -> bool:
    """Whether a regulators.csv belongs to a comparison with a direction: only then are activity_t and the
    signed (consistent / inconsistent) target lists ever filled."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return any(r.get(c) not in (None, "") for r in csv.DictReader(handle)
                   for c in ("activity_t", "consistent_targets", "inconsistent_targets"))


def _identity(finding: Finding) -> FrozenSet[str]:
    """What a finding is compared by: its evidence genes, or -- when it has none (a regulator whose significant
    targets are all unsigned) -- its own name, so the same regulator in two networks is one finding and two
    different ones are two (an empty set would match nothing, and every such finding would count as new)."""
    _, name, genes = finding
    return genes or frozenset({f"\0{name}"})


def independent(findings: Sequence[Finding], threshold: float) -> int:
    """Greedy clustering, lowest q first: a finding within ``threshold`` Jaccard of a kept one joins it."""
    kept: List[FrozenSet[str]] = []
    for f in sorted(findings, key=lambda f: (f[0], f[1])):
        genes = _identity(f)
        if not any(jaccard(genes, k) >= threshold for k in kept):
            kept.append(genes)
    return len(kept)


def score(by_source: Dict[str, List[Finding]], sources: Sequence[str], threshold: float) -> Dict[str, Dict[str, object]]:
    """One comparison (or a pool): the columns above for every source in ``sources``."""
    evidence = {s: frozenset().union(*(g for _, _, g in by_source.get(s, []))) for s in sources}
    everything = frozenset().union(*evidence.values()) if evidence else frozenset()
    result: Dict[str, Dict[str, object]] = {}
    for s in sources:
        mine = by_source.get(s, [])
        others = [_identity(f) for o in sources if o != s for f in by_source.get(o, [])]
        unique = sum(1 for f in mine if not any(jaccard(_identity(f), h) >= threshold for h in others))
        indep = independent(mine, threshold)
        elsewhere = frozenset().union(*(evidence[o] for o in sources if o != s))
        unique_genes = len(evidence[s] - elsewhere)
        n = len(mine)
        result[s] = {
            "findings": n, "independent_findings": indep, "unique_findings": unique,
            "unique_fraction": round(unique / n, 4) if n else "",
            "evidence_genes": len(evidence[s]), "unique_genes": unique_genes,
            "marginal_gain": round(unique_genes / len(everything), 4) if everything else "",
            "utility": round((unique / n) * (indep / n), 4) if n else "",
        }
    return result


def _pool(rows: List[Dict[str, object]]) -> Dict[str, object]:
    """Pooled over comparisons: counts summed, fractions recomputed from the sums."""
    f = sum(int(r["findings"]) for r in rows)
    i = sum(int(r["independent_findings"]) for r in rows)
    u = sum(int(r["unique_findings"]) for r in rows)
    e = sum(int(r["evidence_genes"]) for r in rows)
    ug = sum(int(r["unique_genes"]) for r in rows)
    fam = sum(int(r["_family_genes"]) for r in rows)
    return {"findings": f, "independent_findings": i, "unique_findings": u,
            "unique_fraction": round(u / f, 4) if f else "", "evidence_genes": e, "unique_genes": ug,
            "marginal_gain": round(ug / fam, 4) if fam else "", "utility": round((u / f) * (i / f), 4) if f else ""}


def source_utility(run_dir: Path, *, q_max: float = 0.05, threshold: float = 0.5) -> List[Dict[str, object]]:
    run_dir = Path(run_dir)
    summary = run_dir / "comparison_summary.json"
    if not summary.exists():
        raise ProjectError(21, f"{summary} not found")
    comparisons = [c["comparison"] for c in json.loads(summary.read_text())["comparisons"]]
    rows: List[Dict[str, object]] = []
    per_family: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for comp in comparisons:
        for family, path, reader, sources in (
                ("pathways", run_dir / "enrichment" / comp / "ora.csv", pathway_findings, PATHWAY_SOURCES),
                ("regulators", run_dir / "regulators" / comp / "regulators.csv", regulator_findings,
                 REGULATOR_NETWORKS)):
            if not path.exists():
                continue
            if family == "regulators" and not _directed(path):
                continue      # F-tests without a direction are not counted here, as the module docstring says
            found = reader(path, q_max)
            listed = [s for s in sources if s in found or _source_tested(path, s, family)]
            scored = score(found, listed, threshold)
            family_genes = len(frozenset().union(*(g for s in listed for _, _, g in found.get(s, []))))
            for s in listed:
                row = {"family": family, "source": s, "comparison": comp, **scored[s], "_family_genes": family_genes}
                rows.append(row)
                per_family[(family, s)].append(row)
    pooled = [{"family": fam, "source": s, "comparison": "ALL", **_pool(rs)} for (fam, s), rs in per_family.items()]
    order = {s: i for i, s in enumerate(PATHWAY_SOURCES + REGULATOR_NETWORKS)}
    pooled.sort(key=lambda r: (r["family"], order.get(str(r["source"]), 99)))
    for r in rows:
        r.pop("_family_genes", None)
    return pooled + rows


def _source_tested(path: Path, source: str, family: str) -> bool:
    """A database with no significant finding still gets a row (utility 0 findings) if it was tested."""
    key = "source" if family == "pathways" else "network"
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return any(r.get(key) == source for r in csv.DictReader(handle))


def write_utility(run_dir: Path, rows: List[Dict[str, object]], out: Optional[Path] = None) -> Path:
    out = Path(out) if out else Path(run_dir) / "source_utility.csv"
    tmp = out.with_name(out.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=UTILITY_COLUMNS, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    tmp.replace(out)
    return out


def utility_lines(rows: List[Dict[str, object]], q_max: float, threshold: float) -> List[str]:
    lines = [f"database utility, pooled over comparisons (finding: q < {q_max}; near-duplicate: gene Jaccard >= "
             f"{threshold})",
             f"  {'family':<11}{'database':<19}{'findings':>9}{'independent':>12}{'unique':>8}{'unique%':>9}"
             f"{'marginal gain':>15}{'utility':>9}"]
    for r in rows:
        if r["comparison"] != "ALL":
            continue
        pct = f"{100 * r['unique_fraction']:.0f}%" if r["unique_fraction"] != "" else "-"
        gain = f"{100 * r['marginal_gain']:.1f}%" if r["marginal_gain"] != "" else "-"
        util = f"{r['utility']:.2f}" if r["utility"] != "" else "-"
        lines.append(f"  {r['family']:<11}{r['source']:<19}{r['findings']:>9}{r['independent_findings']:>12}"
                     f"{r['unique_findings']:>8}{pct:>9}{gain:>15}{util:>9}")
    return lines
