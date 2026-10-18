"""The genomes the sister program (Paralome) has registered on this machine.

Paralome keeps a machine-wide genome registry beside its genomes, under the
same communal root this program already uses, visible to every node. Each registered genome records
the exact GTF its data was aligned and quantified against. That annotation
DEFINES an input table's gene IDs, so it is the right thing to identify the
input against:

- the match against the right genome is (near) total, so the report can name
  the genome, not just the species;
- it holds the genome the user actually used (e.g. sheep Rambouillet v2.0),
  where Ensembl 116 has moved on (Rambouillet v3.0);
- with Ensembl 116 still supplying orthology, the difference between the two
  annotations becomes a measured, reported number instead of an unexplained
  loss.

Where the registry is comes from Paralome itself (``setup where``);
the registry files are then read as plain JSON -- Paralome is never imported.
Registry files look like
``{"schema": ..., "registered_by": user, "genomes": {name: {"gtf": path, ...}}}``;
when two files publish the same name, the entry in its registrant's own file
wins, as in Paralome.

Each genome's GTF is indexed once into the reference directory
(``registered_genomes/<name>/genes.tsv``) and reused until the GTF itself
changes (path, size or modification time) or the program version changes --
decided by ``checkpoint.reuse_decision``, the one authority on reuse.
Only genomes whose gene IDs share the input's prefix are indexed at all, found
by reading each GTF's first gene line.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .checkpoint import reuse_decision, write_checkpoint
from .reference import _open_tsv, index_gtf

REGISTRY_KEY = "site_registry"
REGISTRY_DIRNAME = "registry"
WITHDRAWN_KEY = "withdrawn"

_LEADING_LETTERS = re.compile(r"^([A-Za-z]+)")


@dataclass
class GenomeCall:
    name: str
    gtf: str
    matched: int
    total: int
    ranking: List[Tuple[str, int]]
    withdrawn: Optional[Mapping] = None
    notes: List[str] = field(default_factory=list)
    source: str = ""

    @property
    def fraction(self) -> float:
        return self.matched / self.total if self.total else 0.0


def registry_dir() -> Optional[Path]:
    """Paralome's genome registry -- as Paralome reports it -- or None.

    Only for a machine where ``paralome setup site`` has been run: if Paralome
    is not there, or setup site has not been run yet, that is reported, not
    replaced by a guess.
    """
    from .sister import sister

    paralome = sister()
    return paralome.genome_registry if paralome.usable else None


def registered_genomes(directory: Optional[Path] = None) -> Dict[str, Dict]:
    """Every genome registered on this machine, merged across the drop-in files.

    A missing, unreadable or malformed file is skipped. A name published by two
    files resolves to the entry in its registrant's own file.
    """
    directory = directory if directory is not None else registry_dir()
    if directory is None or not Path(directory).is_dir():
        return {}
    merged: Dict[str, Dict] = {}
    owned: Dict[str, bool] = {}
    for path in sorted(Path(directory).glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        genomes = payload.get("genomes") if isinstance(payload, dict) else None
        if not isinstance(genomes, dict):
            continue
        for name, entry in genomes.items():
            if not isinstance(entry, dict):
                continue
            owner = str(entry.get("registered_by") or "")
            owns = not owner or owner == path.stem
            if name not in merged or (owns and not owned.get(name, False)):
                merged[name] = entry
                owned[name] = owns
    return merged


#: The program's OWN registry, inside its reference directory. Same
#: file format as the sister program's. It is what makes a standalone install --
#: no Paralome on the machine -- a first-class case rather than a fallback:
#: register a genome by writing ``<reference_dir>/genome_registry/<user>.json``.
OWN_REGISTRY_DIRNAME = "genome_registry"
SOURCE_SISTER = "Paralome genome registry"
SOURCE_OWN = "own genome registry"


def own_registry_dir(reference_dir: Path) -> Path:
    return Path(reference_dir) / OWN_REGISTRY_DIRNAME


def all_registered_genomes(reference_dir: Path, sister_dir: Optional[Path] = None) -> Dict[str, Dict]:
    """Every genome from every source, each entry tagged with ``_source``.

    Sources, in precedence order (a name defined by both keeps the first):
    the machine's sister registry (the communal record every node sees), then
    this program's own registry. New sources -- e.g. a curated catalogue for a
    standalone deployment -- slot in here and nowhere else.
    """
    merged: Dict[str, Dict] = {}
    for source, directory in ((SOURCE_SISTER, sister_dir if sister_dir is not None else registry_dir()),
                              (SOURCE_OWN, own_registry_dir(reference_dir))):
        for name, entry in registered_genomes(directory).items():
            if name not in merged:
                merged[name] = dict(entry, _source=source)
    return merged


def leading_prefix(gene_id: str) -> str:
    match = _LEADING_LETTERS.match(gene_id)
    return match.group(1) if match else ""


def gtf_first_gene_id(gtf: Path) -> Optional[str]:
    """The first gene_id in a GTF, read without parsing the whole file."""
    try:
        with _open_tsv(Path(gtf)) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if len(fields) >= 9 and fields[2] == "gene":
                    match = re.search(r'gene_id "([^"]+)"', fields[8])
                    return match.group(1) if match else None
    except OSError:
        return None
    return None


def _safe(name: str) -> str:
    """A folder name for a genome's cache. A name that had to be changed gets a short hash of the original, so
    'sheep v2' and 'sheep_v2' do not share one cache; plain names keep their existing folder."""
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", name)
    if safe == name:
        return safe
    import hashlib

    return f"{safe}-{hashlib.sha1(name.encode('utf-8')).hexdigest()[:8]}"


def _index_one(task: Tuple[str, str, str]) -> Tuple[str, Optional[str], str]:
    """Process-pool worker: (name, gtf, cache_dir) -> (name, genes.tsv or None, note)."""
    name, gtf, cache = task
    gtf_path = Path(gtf)
    try:
        stat = gtf_path.stat()
    except OSError as exc:
        return name, None, f"{name}: GTF not readable ({exc.strerror or exc})"
    cache_dir = Path(cache)
    checkpoint = cache_dir / "checkpoint.json"
    genes = cache_dir / "genes.tsv"
    params = {"gtf": str(gtf_path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if reuse_decision(checkpoint, params):
        return name, str(genes), "reused"
    cache_dir.mkdir(parents=True, exist_ok=True)
    count = len(index_gtf(gtf_path, genes, ("gene_id", "symbol", "biotype")))
    write_checkpoint(checkpoint, params, {"genes": str(genes)}, genome=name, genes=count)
    return name, str(genes), f"indexed {count} genes"


def read_gene_table(path: Path) -> Dict[str, Tuple[str, str]]:
    import csv

    table: Dict[str, Tuple[str, str]] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            table[row["gene_id"]] = (row["symbol"], row["biotype"])
    return table


def identify_registered_genome(
    gene_ids: Sequence[str],
    reference_dir: Path,
    *,
    min_fraction: float = 0.9,
    threads: Optional[int] = None,
    registry: Optional[Path] = None,
) -> Tuple[Optional[GenomeCall], Dict[str, Tuple[str, str]], List[str]]:
    """The registered genome whose annotation holds the input IDs, its gene table, and notes.

    Returns (None, {}, notes) when the machine has no registry, no registered
    genome shares the input's ID prefix, or none reaches ``min_fraction`` -- the
    caller then falls back to Ensembl's own annotation.
    """
    from .parallel import map_ordered, resolve_threads

    notes: List[str] = []
    from .sister import sister

    if registry is None and not sister().usable:
        notes.append(sister().summary())
    genomes = all_registered_genomes(reference_dir, registry)
    if not genomes:
        return None, {}, notes + ["no registered genomes found (no Paralome registry, nothing in the own registry)"]
    ids = set(gene_ids)
    wanted_prefix = Counter(leading_prefix(g) for g in ids).most_common(1)[0][0] if ids else ""

    tasks = []
    for name, entry in sorted(genomes.items()):
        gtf = entry.get("gtf")
        if not isinstance(gtf, str) or not gtf:
            continue
        first = gtf_first_gene_id(Path(gtf))
        if first is None:
            notes.append(f"{name}: GTF {gtf} not readable")
            continue
        if leading_prefix(first) != wanted_prefix:
            continue
        tasks.append((name, gtf, str(Path(reference_dir) / "registered_genomes" / _safe(name))))
    if not tasks:
        return None, {}, notes + [f"no registered genome uses {wanted_prefix}* gene IDs"]

    results = map_ordered(_index_one, tasks, threads=resolve_threads(threads), kind="process")
    ranking = []
    tables: Dict[str, Dict[str, Tuple[str, str]]] = {}
    for name, genes, note in results:
        if genes is None:
            notes.append(note)
            continue
        tables[name] = read_gene_table(Path(genes))
        ranking.append((name, len(ids & set(tables[name]))))
    if not ranking:
        return None, {}, notes
    ranking.sort(key=lambda item: (-item[1], item[0]))
    best, matched = ranking[0]
    call = GenomeCall(best, genomes[best]["gtf"], matched, len(ids), ranking,
                      withdrawn=genomes[best].get(WITHDRAWN_KEY), notes=notes,
                      source=genomes[best].get("_source", ""))
    if call.fraction < min_fraction:
        notes.append(f"best registered genome {best} holds {matched}/{len(ids)} IDs "
                     f"({call.fraction:.1%}), below {min_fraction:.0%}")
        return None, {}, notes
    if call.withdrawn:
        call.notes.append(f"WARNING: registered genome {best} has been WITHDRAWN: {call.withdrawn}")
    return call, tables[best], call.notes
