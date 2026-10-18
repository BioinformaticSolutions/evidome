"""Gene-set reference data: download, map to human Ensembl genes, index.

Four sources, all free for academic and commercial use (checked 2026-09-27):

    source        licence     file(s)                                         IDs
    reactome      CC0         Ensembl2Reactome_All_Levels.txt,                 Ensembl
                              ReactomePathways.txt, ReactomePathwaysRelation.txt
    go_bp/mf/cc   CC BY 4.0   NCBI gene2go.gz + go-basic.obo                   NCBI GeneID
    wikipathways  CC0         wikipathways-<date>-gmt-Homo_sapiens.gmt         NCBI GeneID
    hallmark      CC BY 4.0   h.all.v<release>.entrez.gmt (MSigDB)             NCBI GeneID

NCBI GeneIDs are translated with NCBI's own gene2ensembl (already used for
orthology). Every member must be a human gene in the Ensembl release
the crosswalk uses; members that are not are dropped and COUNTED per set
(``source_ids`` vs ``n_genes``), never silently.

GO: direct annotations from gene2go, "NOT" qualifiers excluded, then
propagated to every ancestor through is_a and part_of (the standard true-path
rule); obsolete terms dropped, alternative IDs folded into their primary term.
Reactome: the "All_Levels" file already assigns genes to every ancestor
pathway; each pathway's top-level categories come from the hierarchy file.

None of these sources publishes checksums. As for NCBI: the local copy is kept
and used until ``--refresh``; its md5, size and download time are recorded.

On disk, beside the ortholog data::

    ensembl-116/genesets/
      raw/...                   downloads
      gene_sets.tsv             set_id, source, name, category, n_genes, source_ids
      members.tsv               set_id, human_gene_id  (one row per membership)
      genesets.json             sources, URLs, md5s, licences, counts (stamped)
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

from .errors import ProjectError

REACTOME_BASE = "https://reactome.org/download/current"
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"
WIKIPATHWAYS_DIR = "https://data.wikipathways.org/current/gmt"
MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
DEFAULT_MSIGDB_RELEASE = "2026.1.Hs"

SOURCES = ("reactome", "go_bp", "go_mf", "go_cc", "wikipathways", "hallmark")
LICENCES = {
    "reactome": ("CC0 1.0", "Reactome (https://reactome.org); attribution encouraged, not required."),
    "go": ("CC BY 4.0", "Gene Ontology Consortium (http://geneontology.org), annotations via NCBI Gene gene2go."),
    "wikipathways": ("CC0 1.0", "WikiPathways (https://www.wikipathways.org)."),
    "hallmark": ("CC BY 4.0", "MSigDB Hallmark gene sets, Liberzon et al. 2015, Cell Systems 1:417; "
                              "Broad Institute / UC San Diego (https://www.gsea-msigdb.org)."),
}
GO_NAMESPACES = {"biological_process": "go_bp", "molecular_function": "go_mf", "cellular_component": "go_cc"}
GENE_SET_COLUMNS = ("set_id", "source", "name", "category", "n_genes", "source_ids")
_WP_FILE = re.compile(r"^wikipathways-(\d{8})-gmt-Homo_sapiens\.gmt$")


def genesets_dir(reference_dir: Path, release: int) -> Path:
    from .reference import release_dir

    return release_dir(reference_dir, release) / "genesets"


# ---------------------------------------------------------------------------
# parsers -- each returns {set_id: (source, name, category, set of source IDs)}
# ---------------------------------------------------------------------------

Sets = Dict[str, Tuple[str, str, str, Set[str]]]


def _text(path: Path):
    from .reference import _open_tsv

    return _open_tsv(path)


def parse_gmt(path: Path, source: str, prefix: str) -> Sets:
    """GMT: name <tab> description <tab> member IDs. WikiPathways names are 'Name%WikiPathways_date%WPid%species'."""
    sets: Sets = {}
    with _text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0]:
                continue
            raw_name = fields[0]
            if "%" in raw_name:
                parts = raw_name.split("%")
                name, set_id = parts[0], (parts[2] if len(parts) > 2 else parts[0])
            else:
                name = set_id = raw_name
            members = {m.strip() for m in fields[2:] if m.strip()}
            sets[f"{prefix}:{set_id}"] = (source, name, "", members)
    if not sets:
        raise ProjectError(18, f"{path.name}: no gene sets read (expected GMT: name, description, members)")
    return sets


def parse_reactome(mapping: Path, pathways: Path, relations: Path) -> Sets:
    """Ensembl2Reactome_All_Levels: ensembl_id, pathway_id, url, pathway_name, evidence, species (no header)."""
    names: Dict[str, str] = {}
    with _text(pathways) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "Homo sapiens":
                names[f[0]] = f[1]
    parents: Dict[str, Set[str]] = defaultdict(set)
    with _text(relations) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0].startswith("R-HSA-"):
                parents[f[1]].add(f[0])

    def roots(pid: str) -> Set[str]:
        seen, stack, found = set(), [pid], set()
        while stack:
            p = stack.pop()
            if p in seen:
                continue
            seen.add(p)
            if not parents.get(p):
                found.add(p)
            stack.extend(parents.get(p, ()))
        return found

    members: Dict[str, Set[str]] = defaultdict(set)
    rows = 0
    with _text(mapping) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rows += 1
            if f[5] == "Homo sapiens" and f[0].startswith("ENSG"):
                members[f[1]].add(f[0].split(".", 1)[0])
                names.setdefault(f[1], f[3])
    if not rows:
        raise ProjectError(18, f"{mapping.name}: no rows with the expected 6 tab-separated columns")
    sets: Sets = {}
    for pid, genes in members.items():
        category = ";".join(sorted(names.get(r, r) for r in roots(pid) if r != pid))
        sets[f"REACTOME:{pid}"] = ("reactome", names.get(pid, pid), category, genes)
    return sets


def parse_obo(path: Path) -> Tuple[Dict[str, Tuple[str, str]], Dict[str, Set[str]], Dict[str, str]]:
    """go-basic.obo -> ({term: (name, namespace)}, {term: parents via is_a/part_of}, {alt_id: term})."""
    terms: Dict[str, Tuple[str, str]] = {}
    parents: Dict[str, Set[str]] = defaultdict(set)
    alt: Dict[str, str] = {}
    current: Dict[str, object] = {}

    def flush():
        if current.get("id") and not current.get("obsolete"):
            tid = str(current["id"])
            terms[tid] = (str(current.get("name", "")), str(current.get("namespace", "")))
            parents[tid].update(current.get("parents", set()))  # type: ignore[arg-type]
            for a in current.get("alt", []):  # type: ignore[union-attr]
                alt[a] = tid

    in_term = False
    with _text(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("["):
                if in_term:
                    flush()
                in_term = line == "[Term]"
                current = {"parents": set(), "alt": []}
                continue
            if not in_term or ": " not in line:
                continue
            key, value = line.split(": ", 1)
            value = value.split(" ! ", 1)[0].strip()
            if key == "id":
                current["id"] = value
            elif key == "name":
                current["name"] = value
            elif key == "namespace":
                current["namespace"] = value
            elif key == "is_obsolete" and value == "true":
                current["obsolete"] = True
            elif key == "alt_id":
                current["alt"].append(value)  # type: ignore[union-attr]
            elif key == "is_a":
                current["parents"].add(value.split()[0])  # type: ignore[union-attr]
            elif key == "relationship" and value.startswith("part_of "):
                current["parents"].add(value.split()[1])  # type: ignore[union-attr]
    if in_term:
        flush()
    if not terms:
        raise ProjectError(18, f"{path.name}: no [Term] stanzas read")
    return terms, parents, alt


def parse_go(gene2go: Path, obo: Path) -> Sets:
    """Human GO annotations (gene2go, taxon 9606), NOT excluded, propagated to ancestors."""
    from .reference import HUMAN_TAXON, _ncbi_rows

    terms, parents, alt = parse_obo(obo)
    direct: Dict[str, Set[str]] = defaultdict(set)
    for f, col in _ncbi_rows(gene2go, ("tax_id", "GeneID", "GO_ID", "Qualifier"),
                             keep=lambda line: line.startswith(HUMAN_TAXON + "\t")):
        if f[col["tax_id"]] != HUMAN_TAXON or f[col["Qualifier"]].startswith("NOT"):
            continue
        term = alt.get(f[col["GO_ID"]], f[col["GO_ID"]])
        if term in terms:
            direct[term].add(f[col["GeneID"]])
    ancestors: Dict[str, Set[str]] = {}

    def up(term: str) -> Set[str]:
        if term not in ancestors:
            ancestors[term] = set()  # guard against cycles
            found = set()
            for p in parents.get(term, ()):
                if p in terms:
                    found.add(p)
                    found |= up(p)
            ancestors[term] = found
        return ancestors[term]

    members: Dict[str, Set[str]] = defaultdict(set)
    for term, genes in direct.items():
        for t in {term} | up(term):
            members[t] |= genes
    sets: Sets = {}
    for term, genes in members.items():
        name, namespace = terms[term]
        source = GO_NAMESPACES.get(namespace)
        if source:
            sets[term] = (source, name, namespace, genes)
    return sets


# ---------------------------------------------------------------------------
# download, map, index
# ---------------------------------------------------------------------------

def _human_entrez_to_ensembl(gene2ensembl: Path) -> Dict[str, Set[str]]:
    from .reference import HUMAN_TAXON, _ncbi_rows

    table: Dict[str, Set[str]] = defaultdict(set)
    for f, col in _ncbi_rows(gene2ensembl, ("tax_id", "GeneID", "Ensembl_gene_identifier"),
                             keep=lambda line: line.startswith(HUMAN_TAXON + "\t")):
        ens = f[col["Ensembl_gene_identifier"]]
        if f[col["tax_id"]] == HUMAN_TAXON and ens and ens != "-":
            table[f[col["GeneID"]]].add(ens.split(".", 1)[0])
    return table


def _parse_task(task: Tuple[str, Dict[str, str]]) -> Sets:
    """Process-pool worker: parse one source's raw files."""
    what, paths = task
    if what == "reactome":
        return parse_reactome(Path(paths["mapping"]), Path(paths["pathways"]), Path(paths["relations"]))
    if what == "go":
        return parse_go(Path(paths["gene2go"]), Path(paths["obo"]))
    if what == "wikipathways":
        return parse_gmt(Path(paths["gmt"]), "wikipathways", "WP")
    if what == "hallmark":
        return parse_gmt(Path(paths["gmt"]), "hallmark", "MSIGDB")
    raise ValueError(what)


def wikipathways_file(directory_url: str, local_dir: Optional[Path] = None) -> str:
    """The newest monthly WikiPathways GMT: the newest already in ``local_dir`` when there is one (so a
    re-index without --refresh keeps its month and runs offline), else the newest in the listing."""
    from .reference import _list_directory

    if local_dir is not None and Path(local_dir).is_dir():
        local = sorted(p.name for p in Path(local_dir).iterdir() if _WP_FILE.match(p.name))
        if local:
            return local[-1]
    names = [n for n in _list_directory(directory_url) if _WP_FILE.match(n)]
    if not names:
        raise ProjectError(10, f"{directory_url}: no wikipathways-<date>-gmt-Homo_sapiens.gmt in the listing")
    return sorted(names)[-1]


def fetch_genesets(reference_dir: Path, *, release: int = 116, msigdb_release: str = DEFAULT_MSIGDB_RELEASE,
                   refresh: bool = False, threads: Optional[int] = None, urls: Optional[Dict[str, str]] = None,
                   log=print) -> Path:
    """Download, map to this release's human Ensembl genes, and index every source. Returns genesets.json."""
    from .jsonio import write_json
    from .parallel import MAX_DOWNLOADS, map_ordered, resolve_threads
    from .reference import (_download_ncbi, ensure_writable, load_manifest, ncbi_base_url, release_dir,
                            human_genes_path)
    from .version import __version__

    budget = resolve_threads(threads)
    load_manifest(reference_dir, release)  # EVD-0013 unless `reference fetch` has run
    human_path = human_genes_path(reference_dir, release)
    human = set()
    with human_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            human.add(row["human_gene_id"])
    out = ensure_writable(genesets_dir(reference_dir, release))
    raw = out / "raw"
    u = dict(urls or {})
    wp_dir = u.get("wikipathways_dir", WIKIPATHWAYS_DIR)
    ncbi = u.get("ncbi", ncbi_base_url())
    wanted = {
        "reactome/mapping": f"{u.get('reactome', REACTOME_BASE)}/Ensembl2Reactome_All_Levels.txt",
        "reactome/pathways": f"{u.get('reactome', REACTOME_BASE)}/ReactomePathways.txt",
        "reactome/relations": f"{u.get('reactome', REACTOME_BASE)}/ReactomePathwaysRelation.txt",
        "go/obo": u.get("go_obo", GO_OBO_URL),
        "go/gene2go": f"{ncbi}/gene2go.gz",
        "ncbi/gene2ensembl": f"{ncbi}/gene2ensembl.gz",
        "wikipathways/gmt": f"{wp_dir}/{wikipathways_file(wp_dir, None if refresh else raw / 'wikipathways')}",
        "hallmark/gmt": f"{u.get('msigdb', MSIGDB_BASE)}/{msigdb_release}/h.all.v{msigdb_release}.entrez.gmt",
    }
    # gene2ensembl is shared with the ortholog data: reuse that copy.
    ortholog_ncbi = release_dir(reference_dir, release) / "raw" / "ncbi" / "gene2ensembl.gz"
    tasks = []
    for key, url in wanted.items():
        dest = ortholog_ncbi if key == "ncbi/gene2ensembl" else raw / key.split("/")[0] / url.rsplit("/", 1)[1]
        # The shared gene2ensembl is never replaced from here: reference.json records it as the source
        # of every species' NCBI index, and only 'reference fetch --refresh-ncbi' may change it.
        tasks.append((url, str(dest), refresh and key != "ncbi/gene2ensembl"))
    log(f"gene sets: {len(tasks)} file(s), {min(budget, MAX_DOWNLOADS)} at a time")
    got = dict(zip(wanted, map_ordered(_download_ncbi, tasks, threads=min(budget, MAX_DOWNLOADS), kind="thread")))
    for key, d in got.items():
        log(f"  {'reused' if d['reused'] else 'fetched'} {key} ({Path(str(d['path'])).name}, md5 recorded)")

    entrez = _human_entrez_to_ensembl(Path(str(got["ncbi/gene2ensembl"]["path"])))
    parse_tasks = [
        ("reactome", {"mapping": got["reactome/mapping"]["path"], "pathways": got["reactome/pathways"]["path"],
                      "relations": got["reactome/relations"]["path"]}),
        ("go", {"gene2go": got["go/gene2go"]["path"], "obo": got["go/obo"]["path"]}),
        ("wikipathways", {"gmt": got["wikipathways/gmt"]["path"]}),
        ("hallmark", {"gmt": got["hallmark/gmt"]["path"]}),
    ]
    log(f"gene sets: parsing {len(parse_tasks)} source(s) with up to {min(budget, len(parse_tasks))} process(es)")
    parsed = map_ordered(_parse_task, parse_tasks, threads=budget, kind="process")

    counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    rows: List[Tuple[str, str, str, str, int, int]] = []
    membership: List[Tuple[str, str]] = []
    for (what, _), sets in zip(parse_tasks, parsed):
        for set_id, (source, name, category, ids) in sorted(sets.items()):
            if what == "reactome":
                genes = {g for g in ids if g in human}
            else:
                genes = {e for g in ids for e in entrez.get(g, ()) if e in human}
            counts[source]["sets"] += 1
            counts[source]["source_ids"] += len(ids)
            counts[source]["memberships"] += len(genes)
            if not genes:
                counts[source]["sets_with_no_mapped_gene"] += 1
                continue
            rows.append((set_id, source, name, category, len(genes), len(ids)))
            membership.extend((set_id, g) for g in sorted(genes))

    with (out / "gene_sets.tsv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(GENE_SET_COLUMNS)
        w.writerows(rows)
    with (out / "members.tsv").open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(("set_id", "human_gene_id"))
        w.writerows(membership)

    licence_key = {"reactome": "reactome", "go": "go", "wikipathways": "wikipathways", "hallmark": "hallmark",
                   "ncbi": "go"}
    manifest = {
        "ensembl_release": release,
        "indexed_by": __version__,
        "indexed": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "msigdb_release": msigdb_release,
        "files": {k: {f: d[f] for f in ("url", "checksum", "bytes", "downloaded")} for k, d in got.items()},
        "licences": {k: {"licence": LICENCES[licence_key[k.split("/")[0]]][0],
                         "attribution": LICENCES[licence_key[k.split("/")[0]]][1]} for k in got},
        "counts": {s: dict(counts[s]) for s in SOURCES if s in counts},
        "rules": {
            "go": "gene2go, human, 'NOT' qualifiers excluded, propagated through is_a and part_of; "
                  "obsolete terms dropped, alt_ids folded",
            "members": "mapped to human Ensembl genes of this release; unmapped members counted, not kept",
        },
    }
    path = write_json(out / "genesets.json", manifest)
    for source in SOURCES:
        c = counts.get(source)
        if c:
            log(f"  {source:<13} {c['sets']:>6} sets  {c['memberships']:>8} memberships  "
                f"({c.get('sets_with_no_mapped_gene', 0)} sets with no gene in Ensembl {release})")
    return path


def load_gene_sets(reference_dir: Path, release: int = 116,
                   sources: Optional[Iterable[str]] = None) -> Tuple[Dict[str, Dict[str, object]], Dict[str, List[str]]]:
    """({set_id: info}, {set_id: sorted human gene IDs}) from the index; EVD-0019 if absent."""
    out = genesets_dir(reference_dir, release)
    if not (out / "gene_sets.tsv").exists():
        raise ProjectError(19, f"no gene-set index at {out}")
    wanted = set(sources) if sources else None
    info: Dict[str, Dict[str, object]] = {}
    with (out / "gene_sets.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if wanted is None or row["source"] in wanted:
                info[row["set_id"]] = row
    members: Dict[str, List[str]] = defaultdict(list)
    with (out / "members.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["set_id"] in info:
                members[row["set_id"]].append(row["human_gene_id"])
    # Reactome activity patterns, where `reference pathway-signs` has built them
    from .pathway_signs import load_patterns

    for sid, pattern in load_patterns(out, info).items():
        info[sid]["pattern"] = pattern
    return info, dict(members)
