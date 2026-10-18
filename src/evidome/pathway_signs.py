"""Reactome activity patterns: which way each gene is expected to move when a pathway is active.

A pathway's activation is scored from an activity pattern -- for each
molecule of the pathway, whether it is expected up (+) or down (-) when the
pathway is active -- and the dataset's direction for each significant
molecule:

    z = (consistent - inconsistent) / sqrt(N)

over the N significant dataset molecules that have a pattern, and no z where
too few do. The patterns come from Reactome's own reaction logic, which records for
every reaction what goes in, what comes out, what catalyses it, and what
regulates it positively or NEGATIVELY. From it, for each Reactome pathway P
(and everything under it):

    +1   a gene whose product is an input, output, catalyst or positive
         regulator of a reaction in P -- it takes part in P running
    -1   a gene whose product is a NEGATIVE regulator (BioPAX controlType
         INHIBITION*) of a reaction in P -- it brakes P: the "repressing a
         repressor" case, which a count of up and down genes gets wrong
     -   a gene with both roles in P: no pattern (counted, not guessed)

A control on a control takes the product of the signs down the chain:
an activator of an inhibition brakes P (-1), and an inhibitor of
an inhibition releases the brake (+1).

Proteins are followed through complexes (component) and sets
(memberPhysicalEntity / memberEntityReference) to UniProt accessions, and
those to human Ensembl genes with Ensembl's own UniProt cross-references for
the release (an accession on more than one gene is dropped as ambiguous).
Small molecules have no gene and are ignored.

Sources (both free for any use): Reactome's BioPAX Level 3 export
(``biopax.zip``, member ``Homo_sapiens.owl``; CC0) and Ensembl's
``Homo_sapiens.GRCh38.<release>.uniprot.tsv.gz``. BioPAX is a formal standard
(Level 3); the parser reads it by that standard, streaming, and refuses a file
in which it finds no pathways or no negative regulation (EVD-0018), reporting
what it did find.

On disk, beside the gene sets::

    ensembl-<release>/genesets/
      raw/reactome/biopax.zip, raw/ensembl/<uniprot tsv>
      pathway_signs.tsv     set_id, human_gene_id, sign (+1/-1)
      pathway_signs.json    sources, counts (stamped)

Used by ``enrich`` for Reactome sets when present (``genesets.load_gene_sets``
attaches each set's pattern).
"""

from __future__ import annotations

import csv
import json
import gzip
import re
import zipfile
import zlib
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple
from xml.etree.ElementTree import iterparse

from .errors import ProjectError

REACTOME_BIOPAX = "https://reactome.org/download/current/biopax.zip"
SIGNS_COLUMNS = ("set_id", "human_gene_id", "sign")
RDF = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}"
#: Minimum significant genes with a pattern for a pathway z.
PATTERN_MIN = 4

_INTERACTION_LINKS = {"left", "right", "participant", "product", "template"}
_CONTROL_TAGS = {"Catalysis", "Control", "Modulation", "TemplateReactionRegulation"}
_ENTITY_LINKS = {"entityReference", "component", "memberPhysicalEntity", "memberEntityReference", "xref"}
_KEEP_LINKS = _INTERACTION_LINKS | _ENTITY_LINKS | {"pathwayComponent", "controlled", "controller"}
_KEEP_TEXT = {"controlType", "db", "id"}
_ISOFORM = re.compile(r"-\d+$")


def uniprot_tsv_url(base: str, release: int) -> str:
    return f"{base}/release-{release}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{release}.uniprot.tsv.gz"


def _local(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _ref(value: Optional[str]) -> Optional[str]:
    """'#Protein1', 'http://...#Protein1' and 'Protein1' are one object: its fragment."""
    if not value:
        return None
    return value.rsplit("#", 1)[-1]


# ---------------------------------------------------------------------------
# BioPAX, streamed
# ---------------------------------------------------------------------------

class Model:
    """The few BioPAX objects and properties the patterns need, by object ID."""

    def __init__(self) -> None:
        self.kind: Dict[str, str] = {}
        self.links: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
        self.text: Dict[str, Dict[str, str]] = defaultdict(dict)

    @classmethod
    def read(cls, handle) -> "Model":
        model = cls()
        depth = 0
        root = None
        for event, elem in iterparse(handle, events=("start", "end")):
            if event == "start":
                if root is None:
                    root = elem
                depth += 1
                continue
            depth -= 1
            if depth != 1:              # only the objects directly under rdf:RDF; their children are read with them
                continue
            oid = _ref(elem.get(f"{RDF}ID") or elem.get(f"{RDF}about"))
            if oid:
                model.kind[oid] = _local(elem.tag)
                for child in elem:
                    name = _local(child.tag)
                    target = _ref(child.get(f"{RDF}resource"))
                    if target and name in _KEEP_LINKS:
                        model.links[oid][name].append(target)
                    elif name in _KEEP_TEXT and child.text:
                        model.text[oid][name] = child.text.strip()
            if root is not None:
                root.clear()            # the object is read; drop it and everything before it
        return model


def _uniprots(model: Model, oid: str, memo: Dict[str, Set[str]], seen: Optional[Set[str]] = None) -> Set[str]:
    """UniProt accessions reached from an entity or reference through components, members and references."""
    if oid in memo:
        return memo[oid]
    seen = seen if seen is not None else set()
    if oid in seen:
        return set()
    seen.add(oid)
    found: Set[str] = set()
    links = model.links.get(oid, {})
    for x in links.get("xref", ()):
        t = model.text.get(x, {})
        if model.kind.get(x) == "UnificationXref" and t.get("db", "").lower().startswith("uniprot") and t.get("id"):
            found.add(_ISOFORM.sub("", t["id"].strip()))
    for name in ("entityReference", "component", "memberPhysicalEntity", "memberEntityReference"):
        for target in links.get(name, ()):
            found |= _uniprots(model, target, memo, seen)
    memo[oid] = found
    return found


def _stable_id(model: Model, oid: str) -> Optional[str]:
    for x in model.links.get(oid, {}).get("xref", ()):
        t = model.text.get(x, {})
        if t.get("db", "").lower() == "reactome" and t.get("id", "").startswith("R-HSA-"):
            return t["id"].split(".", 1)[0]
    return None


def pathway_patterns(model: Model) -> Tuple[Dict[str, Dict[str, int]], Dict[str, int]]:
    """{R-HSA pathway: {uniprot: +1/-1}} and counts. Conflicting roles within a pathway are dropped."""
    controls_on: Dict[str, List[str]] = defaultdict(list)
    counts: Dict[str, int] = defaultdict(int)
    for oid, kind in model.kind.items():
        if kind in _CONTROL_TAGS:
            ctype = model.text.get(oid, {}).get("controlType", "ACTIVATION").upper()
            counts["controls_inhibition" if ctype.startswith("INHIBITION") else "controls_activation"] += 1
            for target in model.links.get(oid, {}).get("controlled", ()):
                controls_on[target].append(oid)
    memo: Dict[str, Set[str]] = {}
    closure: Dict[str, Set[str]] = {}

    def contents(pid: str, stack: Set[str]) -> Set[str]:
        """Every pathway and interaction under ``pid``, itself included."""
        if pid in closure:
            return closure[pid]
        out = {pid}
        for c in model.links.get(pid, {}).get("pathwayComponent", ()):
            if c in stack:
                continue
            out |= contents(c, stack | {c}) if model.kind.get(c) == "Pathway" else {c}
        closure[pid] = out
        return out

    patterns: Dict[str, Dict[str, int]] = {}
    for pid, kind in model.kind.items():
        if kind != "Pathway":
            continue
        stable = _stable_id(model, pid)
        if not stable:
            counts["pathways_without_stable_id"] += 1
            continue
        counts["pathways"] += 1
        plus: Set[str] = set()
        minus: Set[str] = set()
        for item in contents(pid, {pid}):
            k = model.kind.get(item)
            if k in _CONTROL_TAGS:
                continue            # reached through controls_on, with its target
            if k != "Pathway":
                for name in _INTERACTION_LINKS:
                    for e in model.links.get(item, {}).get(name, ()):
                        plus |= _uniprots(model, e, memo)
            # controls on the step, and controls on those controls (e.g. a negative regulator of a catalysis).
            # The sign multiplies down the chain: an inhibitor of an inhibition releases a brake and is +1, an
            # activator of an inhibition is -1.
            pending = [(c, 1) for c in controls_on.get(item, ())]
            done: Set[Tuple[str, int]] = set()
            while pending:
                ctl, outer = pending.pop()
                if (ctl, outer) in done:
                    continue
                done.add((ctl, outer))
                ctype = model.text.get(ctl, {}).get("controlType", "ACTIVATION").upper()
                sign = outer * (-1 if ctype.startswith("INHIBITION") else 1)
                pending.extend((c, sign) for c in controls_on.get(ctl, ()))
                genes: Set[str] = set()
                for e in model.links.get(ctl, {}).get("controller", ()):
                    genes |= _uniprots(model, e, memo)
                (minus if sign < 0 else plus).update(genes)
        both = plus & minus
        counts["conflicting_roles_dropped"] += len(both)
        pattern = {u: 1 for u in plus - both}
        pattern.update({u: -1 for u in minus - both})
        if pattern:
            patterns[stable] = pattern
    return patterns, dict(counts)


def _damaged(path: Path, exc: BaseException) -> ProjectError:
    return ProjectError(18, f"{Path(path).name} is damaged or not what was expected ({exc}); run 'reference "
                            "pathway-signs --refresh' to download it again")


def read_biopax(path: Path) -> Model:
    """The human model from Reactome's biopax.zip (member Homo_sapiens.owl), or a plain .owl file. A damaged
    download is EVD-0018 naming --refresh, not a raw zipfile or XML error on every rerun."""
    from xml.etree.ElementTree import ParseError

    try:
        return _read_biopax(Path(path))
    except (zipfile.BadZipFile, ParseError, EOFError, OSError) as exc:
        raise _damaged(path, exc) from exc


def _read_biopax(path: Path) -> Model:
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = [n for n in zf.namelist() if n.lower().rsplit("/", 1)[-1] == "homo_sapiens.owl"]
            if not names:
                raise ProjectError(18, f"{path.name}: no Homo_sapiens.owl among {len(zf.namelist())} members "
                                       f"(first: {', '.join(zf.namelist()[:5])})")
            with zf.open(names[0]) as handle:
                return Model.read(handle)
    with path.open("rb") as handle:
        return Model.read(handle)


def uniprot_to_ensembl(path: Path) -> Tuple[Dict[str, str], int]:
    """Ensembl's UniProt cross-references -> {accession: gene}; accessions on two or more genes are dropped."""
    try:
        return _uniprot_to_ensembl(path)
    except (EOFError, OSError, UnicodeDecodeError, zlib.error) as exc:
        raise _damaged(path, exc) from exc


def _uniprot_to_ensembl(path: Path) -> Tuple[Dict[str, str], int]:
    genes: Dict[str, Set[str]] = defaultdict(set)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for needed in ("gene_stable_id", "xref"):
            if needed not in col:
                raise ProjectError(18, f"{Path(path).name}: no {needed!r} column (header: {', '.join(header[:9])})")
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(col["gene_stable_id"], col["xref"]):
                genes[_ISOFORM.sub("", f[col["xref"]].strip())].add(f[col["gene_stable_id"]].split(".", 1)[0])
    unique = {acc: next(iter(g)) for acc, g in genes.items() if len(g) == 1}
    return unique, len(genes) - len(unique)


def build_signs(biopax: Path, uniprot_tsv: Path, human: Set[str], reactome_sets: Set[str]) -> Tuple[List[Tuple[str, str, int]], Dict[str, int]]:
    """(rows, counts). Only pathways that are gene sets in the index (REACTOME:<id>) and genes of this release."""
    model = read_biopax(biopax)
    patterns, found = pathway_patterns(model)
    counts: Dict[str, int] = defaultdict(int, found)
    if not counts.get("pathways") or not counts.get("controls_inhibition"):
        raise ProjectError(18, f"{Path(biopax).name}: not the BioPAX layout expected -- pathways with a Reactome "
                               f"ID {counts.get('pathways', 0)}, negative regulations "
                               f"{counts.get('controls_inhibition', 0)}, objects read {len(model.kind)}")
    to_gene, ambiguous = uniprot_to_ensembl(uniprot_tsv)
    counts["uniprot_ambiguous_dropped"] = ambiguous
    rows: List[Tuple[str, str, int]] = []
    unmapped: Set[str] = set()
    for stable, pattern in sorted(patterns.items()):
        set_id = f"REACTOME:{stable}"
        if set_id not in reactome_sets:
            counts["pathways_not_in_gene_sets"] += 1
            continue
        by_gene: Dict[str, Set[int]] = defaultdict(set)
        for acc, sign in pattern.items():
            gene = to_gene.get(acc)
            if gene is None or gene not in human:
                unmapped.add(acc)
                continue
            by_gene[gene].add(sign)
        kept = 0
        for gene, signs in sorted(by_gene.items()):
            if len(signs) == 1:
                rows.append((set_id, gene, signs.pop()))
                kept += 1
            else:
                counts["conflicting_roles_dropped"] += 1       # two proteins of one gene, opposite roles
        if kept:
            counts["pathways_with_pattern"] += 1
    counts["uniprot_unmapped"] = len(unmapped)
    counts["pattern_rows"] = len(rows)
    counts["negative_pattern_rows"] = sum(1 for r in rows if r[2] < 0)
    return rows, dict(counts)


def fetch_pathway_signs(reference_dir: Path, *, release: int = 116, refresh: bool = False,
                        urls: Optional[Dict[str, str]] = None, log=print) -> Path:
    """Download Reactome BioPAX and Ensembl's UniProt table, build pathway_signs.tsv. Needs `reference genesets`."""
    from .genesets import genesets_dir
    from .jsonio import write_json
    from .reference import DEFAULT_BASE_URL, _download_ncbi, ensure_writable, human_genes_path, release_dir
    from .version import __version__

    out = ensure_writable(genesets_dir(reference_dir, release))
    if not (out / "gene_sets.tsv").exists():
        raise ProjectError(19, f"no gene-set index at {out}")
    human: Set[str] = set()
    with human_genes_path(reference_dir, release).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            human.add(row["human_gene_id"])
    with (out / "gene_sets.tsv").open(encoding="utf-8", newline="") as handle:
        reactome_sets = {r["set_id"] for r in csv.DictReader(handle, delimiter="\t") if r["source"] == "reactome"}
    u = dict(urls or {})
    wanted = {"reactome/biopax": u.get("biopax", REACTOME_BIOPAX),
              "ensembl/uniprot": u.get("uniprot", uniprot_tsv_url(u.get("ensembl", DEFAULT_BASE_URL), release))}
    got = {}
    for key, url in wanted.items():
        got[key] = _download_ncbi((url, str(out / "raw" / key.split("/")[0] / url.rsplit("/", 1)[1]), refresh))
        log(f"  {'reused' if got[key]['reused'] else 'fetched'} {key} ({Path(str(got[key]['path'])).name})")
    log("pathway patterns: reading Reactome BioPAX (streamed)")
    rows, counts = build_signs(Path(str(got["reactome/biopax"]["path"])), Path(str(got["ensembl/uniprot"]["path"])),
                               human, reactome_sets)
    tmp = out / "pathway_signs.tsv.tmp"
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(SIGNS_COLUMNS)
        w.writerows(rows)
    tmp.replace(out / "pathway_signs.tsv")
    manifest = {
        "ensembl_release": release, "indexed_by": __version__,
        "indexed": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": {k: {f: d[f] for f in ("url", "checksum", "bytes", "downloaded")} for k, d in got.items()},
        "licences": {"reactome/biopax": "CC0 1.0 (Reactome)", "ensembl/uniprot": "Ensembl: no restrictions"},
        "rule": "+1 input/output/catalyst/positive regulator of a reaction in the pathway; -1 negative regulator "
                "(controlType INHIBITION*); a control on a control takes the product of the signs down the chain "
                "(an inhibitor of an inhibition is +1); both roles: dropped",
        "counts": counts,
    }
    path = write_json(out / "pathway_signs.json", manifest)
    log(f"  {counts.get('pathways_with_pattern', 0)} Reactome pathways with a pattern; "
        f"{counts.get('pattern_rows', 0)} gene signs ({counts.get('negative_pattern_rows', 0)} negative); "
        f"{counts.get('conflicting_roles_dropped', 0)} conflicting roles dropped")
    return path


#: The first version whose sign rule multiplies down chains of controls.
SIGN_RULE_SINCE = 39


def stale_signs_note(directory: Path) -> Optional[str]:
    """A warning when pathway_signs.tsv was built under the older one-level sign rule, by a pre-release
    version older than SIGN_RULE_SINCE."""
    path = Path(directory) / "pathway_signs.json"
    if not path.exists():
        return None
    try:
        indexed_by = str(json.loads(path.read_text(encoding="utf-8")).get("indexed_by", ""))
    except (OSError, ValueError, AttributeError):
        return None
    m = re.match(r"^0\.1\.dev(\d+)$", indexed_by)
    if m and int(m.group(1)) < SIGN_RULE_SINCE:
        return (f"note: the Reactome activity patterns were built by {indexed_by}, under the older sign rule; "
                "run 'reference pathway-signs' again (it rebuilds them from the files already downloaded) before "
                "trusting the activation z")
    return None


def load_patterns(directory: Path, set_ids: Iterable[str]) -> Dict[str, Dict[str, int]]:
    """{set_id: {human gene: sign}} for the given sets, from pathway_signs.tsv; {} if it was never built."""
    path = Path(directory) / "pathway_signs.tsv"
    if not path.exists():
        return {}
    wanted = set(set_ids)
    out: Dict[str, Dict[str, int]] = defaultdict(dict)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["set_id"] in wanted:
                out[row["set_id"]][row["human_gene_id"]] = int(row["sign"])
    return dict(out)
