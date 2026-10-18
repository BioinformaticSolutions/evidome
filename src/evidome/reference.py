"""Reference data: download, verify and index Ensembl 116 Compara.

Runs on a machine with network access; the data is downloaded where it is used.

Layout observed on the real FTP site, not assumed::

    {base}/release-116/tsv/ensembl-compara/homologies/{species}/
        Compara.116.protein_default.homologies.tsv.gz
        Compara.116.ncrna_default.homologies.tsv.gz
        MD5SUM   CHECKSUMS   README.gene_trees.tsv_dumps.txt

One directory per species (223 in e116), no assembly subdirectory, tab-separated,
header:

    gene_stable_id protein_stable_id species identity homology_type
    homology_gene_stable_id homology_protein_stable_id homology_species
    homology_identity dn ds goc_score wga_coverage is_high_confidence homology_id

EACH PAIR IS STORED ONCE, AND ITS ORIENTATION IS ARBITRARY. Reading only a
species' own files indexes macaque with 0 orthologs to human: the macaque
file has 3.3 million rows against 146 partner species and not one of them is
human. The macaque-human pairs live in the HUMAN directory's files. Ensembl's
own README says why the code must not assume either side: "the order of the
first and second homologous genes within each row is arbitrary". So a species'
orthologs to human are collected from its own files AND the human files, each
pair is oriented (input species on one side, human on the other, identities
swapped to match), and duplicates are dropped by Ensembl's unique homology_id.
The human files (~110 MB + 33 MB in e116) are downloaded once and shared by
every species.

Human gene names come from the same release's human GTF, so every name in a
report belongs to the same Ensembl release as the orthology that produced it.

On disk, under the reference directory::

    ensembl-116/
      raw/{species}/...                    downloads, kept for re-indexing
      species/{species}/orthologs_to_human.tsv   ortholog rows only, to human
      species/{species}/genes.txt          every gene ID of that species seen
                                           in its homology files (membership)
      human_genes.tsv                      human gene ID, symbol, biotype
      reference.json                       what was fetched, from where, the
                                           checksums, row counts (stamped)

Every download is checksum-verified against the directory's MD5SUM (or, if
there is none, CHECKSUMS with the BSD ``sum`` algorithm). A file that cannot be
verified is refused, not used.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import http.client
import io
import os
import shutil
import urllib.request
import zlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .jsonio import read_json, write_json
from .wide_stats import strip_version

DEFAULT_BASE_URL = "https://ftp.ensembl.org/pub"
#: NCBI Gene's bulk files. Read through ``ncbi_base_url()`` at call
#: time, so tests can point it at a local mirror.
DEFAULT_NCBI_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA"
#: Ensembl's REST service: which species a gene ID belongs to, for any species Ensembl has.
DEFAULT_REST_URL = "https://rest.ensembl.org"
NCBI_FILES = ("gene2ensembl.gz", "gene_orthologs.gz")
HUMAN_TAXON = "9606"
NCBI_ORTHOLOG_COLUMNS = ("gene_id", "human_gene_id", "ncbi_gene_id", "human_ncbi_gene_id")
DEFAULT_RELEASE = 116
HOMOLOGY_KINDS = ("protein", "ncrna")
HUMAN = "homo_sapiens"

ORTHOLOG_COLUMNS = (
    "gene_id", "human_gene_id", "homology_type", "is_high_confidence",
    "identity", "human_identity", "goc_score", "wga_coverage", "source",
)


# ---------------------------------------------------------------------------
# where the reference data lives
# ---------------------------------------------------------------------------

def resolve_reference_dir(explicit: Optional[Path] = None) -> Tuple[Path, str]:
    """Where the reference data lives, and which rule decided it.

    1. ``--reference-dir`` on the command line -- an explicit, one-command
       choice (testing, a one-off).
    2. The communal root Paralome reports for a set-up machine (asked, never
       re-derived -- see sister.py): ``<shared_root>/<name>/``, a sibling of
       Paralome's ``genomes/`` and ``registry/``. This program is
       Paralome's sister, so its reference data lives in the same communal
       place as the genomes. When the machine has decided, the environment
       variable is ignored -- the rule Paralome applies ("site.json wins
       outright").
    3. The REFERENCE environment variable (prefixed with the program name).
    """
    from .identity import NAME, env_var
    from .sister import sister

    if explicit is not None:
        return Path(explicit), "--reference-dir"
    paralome = sister()
    if paralome.usable and paralome.shared_root is not None:
        return paralome.shared_root / NAME, f"Paralome's shared root ({paralome.command})"
    value = os.environ.get(env_var("REFERENCE"), "").strip()
    if value:
        return Path(value), f"{env_var('REFERENCE')} ({paralome.summary()})"
    raise ProjectError(
        13,
        f"no --reference-dir, no ${env_var('REFERENCE')}, and no shared root from Paralome: {paralome.summary()}",
    )


def release_dir(reference_dir: Path, release: int = DEFAULT_RELEASE) -> Path:
    return Path(reference_dir) / f"ensembl-{release}"


def human_genes_path(reference_dir: Path, release: int = DEFAULT_RELEASE) -> Path:
    """The human gene index, or EVD-0013 if it is missing (never a raw
    FileNotFoundError)."""
    path = release_dir(reference_dir, release) / "human_genes.tsv"
    if not path.exists():
        raise ProjectError(13, f"no human gene index at {path}")
    return path


def homology_url(base: str, release: int, species: str, kind: str) -> str:
    return f"{base}/release-{release}/tsv/ensembl-compara/homologies/{species}/Compara.{release}.{kind}_default.homologies.tsv.gz"


def human_gtf_url(base: str, release: int) -> Tuple[str, str]:
    name = f"Homo_sapiens.GRCh38.{release}.gtf.gz"
    return f"{base}/release-{release}/gtf/{HUMAN}/{name}", name


# ---------------------------------------------------------------------------
# download + verify
# ---------------------------------------------------------------------------

def _open_url(url: str):
    """The one way this program opens a URL: it names itself.

    reactome.org answers 403 to Python's default
    client name ("Python-urllib/3.9") and 200 to a named client -- checked
    with curl both ways. Naming the program and its version is also what
    public data servers ask of scripted downloads. No personal contact
    details are sent.
    """
    from .identity import DISPLAY_NAME
    from .version import __version__

    request = urllib.request.Request(url, headers={
        "User-Agent": f"{DISPLAY_NAME}/{__version__} (bioinformatics reference-data download)"})
    return urllib.request.urlopen(request)


def _download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    try:
        with _open_url(url) as response, tmp.open("wb") as handle:
            announced = response.headers.get("Content-Length") if getattr(response, "headers", None) else None
            shutil.copyfileobj(response, handle, length=1 << 20)
        # A connection closed early ends the read without an error; without this check the short file would be
        # kept and reused on every later run. Only Ensembl files have published checksums to catch it.
        if announced and announced.isdigit() and tmp.stat().st_size != int(announced):
            raise ProjectError(10, f"{url}: received {tmp.stat().st_size} of {announced} bytes; the connection "
                                   "closed early")
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise ProjectError(10, f"{url}: {exc}") from exc
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    tmp.replace(dest)
    return dest


def _fetch_text(url: str) -> Optional[str]:
    try:
        with _open_url(url) as response:
            return response.read().decode("utf-8", "replace")
    except OSError:
        return None


def bsd_sum(path: Path) -> Tuple[int, int]:
    """``sum -r``: BSD 16-bit rotating checksum and 1 KiB block count."""
    checksum = 0
    size = 0
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            size += len(chunk)
            for byte in chunk:
                checksum = ((checksum >> 1) + ((checksum & 1) << 15) + byte) & 0xFFFF
    return checksum, (size + 1023) // 1024


def md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_sums(text: str) -> Dict[str, List[str]]:
    table: Dict[str, List[str]] = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            table[parts[-1].lstrip("*")] = parts[:-1]
    return table


def verify(path: Path, directory_url: str) -> Dict[str, str]:
    """Verify ``path`` against the MD5SUM (preferred) or CHECKSUMS in its directory."""
    name = Path(path).name
    md5_text = _fetch_text(f"{directory_url}/MD5SUM")
    if md5_text:
        expected = _parse_sums(md5_text).get(name)
        if expected:
            actual = md5_file(path)
            if actual != expected[0]:
                raise ProjectError(11, f"{name}: md5 {actual}, MD5SUM says {expected[0]}")
            return {"method": "md5", "value": actual}
    sums_text = _fetch_text(f"{directory_url}/CHECKSUMS")
    if sums_text:
        expected = _parse_sums(sums_text).get(name)
        if expected and len(expected) >= 2:
            checksum, blocks = bsd_sum(path)
            if (checksum, blocks) != (int(expected[0]), int(expected[1])):
                raise ProjectError(11, f"{name}: sum {checksum} {blocks}, CHECKSUMS says {' '.join(expected)}")
            return {"method": "bsd-sum", "value": f"{checksum} {blocks}"}
    raise ProjectError(11, f"{name}: no MD5SUM or CHECKSUMS entry in {directory_url}")


# ---------------------------------------------------------------------------
# indexing
# ---------------------------------------------------------------------------

def _open_tsv(path: Path):
    """Text handle on a gzipped or plain file (Paralome's genome GTFs are uncompressed)."""
    path = Path(path)
    with path.open("rb") as probe:
        gzipped = probe.read(2) == b"\x1f\x8b"
    if gzipped:
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return path.open(encoding="utf-8", newline="")


def index_homologies(raw_files: Iterable[Tuple[str, Path]], species: str, out_dir: Path) -> Dict[str, int]:
    """Orthologs between ``species`` and human, from any files, in either orientation.

    ``raw_files`` is (kind, path) for the species' own homology files and the
    human ones. A row counts when one side is ``species`` and the other is human;
    it is written with the ``species`` gene first. Also lists every ``species``
    gene seen on either side of any row (for membership).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    genes = set()
    seen_ids = set()
    kept = 0
    duplicates = 0
    by_type: Dict[str, int] = {}
    oriented = {"species_first": 0, "human_first": 0}
    ortholog_path = out_dir / "orthologs_to_human.tsv"
    with ortholog_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(ORTHOLOG_COLUMNS)
        for kind, path in raw_files:
            with _open_tsv(path) as handle:
                reader = csv.reader(handle, delimiter="\t")
                header = next(reader, [])
                required = ("gene_stable_id", "species", "identity", "homology_type", "homology_gene_stable_id",
                            "homology_species", "homology_identity", "is_high_confidence")
                missing = [c for c in required if c not in header]
                if missing:
                    raise ProjectError(12, f"{Path(path).name}: missing columns {', '.join(missing)}")
                col = {name: i for i, name in enumerate(header)}
                g, s, ht, hg, hs = (col["gene_stable_id"], col["species"], col["homology_type"],
                                    col["homology_gene_stable_id"], col["homology_species"])
                ident, hident, hc = col["identity"], col["homology_identity"], col["is_high_confidence"]
                goc, wga, hid = col.get("goc_score"), col.get("wga_coverage"), col.get("homology_id")
                # Millions of rows: index positions, not a dict per row.
                for row in reader:
                    first, second = row[s], row[hs]
                    if first == species:
                        genes.add(row[g])
                    if second == species:
                        genes.add(row[hg])
                    if first == species and second == HUMAN:
                        mine, theirs, my_ident, their_ident = row[g], row[hg], row[ident], row[hident]
                        orientation = "species_first"
                    elif first == HUMAN and second == species:
                        mine, theirs, my_ident, their_ident = row[hg], row[g], row[hident], row[ident]
                        orientation = "human_first"
                    else:
                        continue
                    if not row[ht].startswith("ortholog"):
                        continue
                    key = row[hid] if hid is not None else (mine, theirs, row[ht])
                    if key in seen_ids:
                        duplicates += 1
                        continue
                    seen_ids.add(key)
                    writer.writerow([
                        mine, theirs, row[ht], row[hc], my_ident, their_ident,
                        row[goc] if goc is not None else "", row[wga] if wga is not None else "", kind,
                    ])
                    kept += 1
                    oriented[orientation] += 1
                    by_type[row[ht]] = by_type.get(row[ht], 0) + 1
    (out_dir / "genes.txt").write_text("\n".join(sorted(genes)) + "\n", encoding="utf-8")
    return {"genes": len(genes), "orthologs_to_human": kept, "duplicate_rows_skipped": duplicates,
            **{f"rows_{k}": v for k, v in oriented.items()},
            **{f"type:{k}": v for k, v in sorted(by_type.items())}}


def index_gtf(gtf: Path, out: Path, columns: Tuple[str, str, str]) -> set:
    """Gene ID, symbol and biotype for every gene line of a GTF. Returns the gene IDs."""
    ids = set()
    with _open_tsv(gtf) as handle, out.open("w", newline="", encoding="utf-8") as dest:
        writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = {}
            for part in fields[8].split(";"):
                part = part.strip()
                if " " in part:
                    key, value = part.split(" ", 1)
                    attrs[key] = value.strip('"')
            # GENCODE versions its gene IDs (ENSG...1.2) and calls the biotype gene_type; input IDs are read with
            # their versions stripped, so a versioned index would match nothing.
            gene_id = strip_version(attrs.get("gene_id", ""))[0]
            biotype = attrs.get("gene_biotype") or attrs.get("gene_type", "")
            writer.writerow((gene_id, attrs.get("gene_name", ""), biotype))
            ids.add(gene_id)
    return ids


def index_human_gtf(gtf: Path, out: Path) -> int:
    """Human gene ID -> symbol and biotype, from the release's own GTF."""
    return len(index_gtf(gtf, out, ("human_gene_id", "human_symbol", "human_biotype")))


def species_gtf_url(base: str, release: int, species: str) -> str:
    """The species' main gene annotation for the release, found by listing its GTF directory.

    Needed because a species' homology files only mention genes that HAVE a
    homolog. For macaque, the homology-derived gene list holds 29,563 genes
    while a real expression table has 35,432; the other 5,869 are genes Compara
    builds no tree for (lncRNAs, pseudogenes, ...). Membership has to be judged against
    the full annotation. The directory also holds ``.chr.``, ``.abinitio.`` and
    patch variants (as for sheep); the main file is the one without.
    """
    import re

    directory = f"{base}/release-{release}/gtf/{species}"
    pattern = re.compile(rf"^[A-Z][A-Za-z0-9_]*\.[^/]+\.{release}\.gtf\.gz$")
    candidates = [name for name in _list_directory(directory)
                  if pattern.match(name) and not any(tag in name for tag in (".chr.", ".abinitio.", ".chr_patch"))]
    if len(candidates) != 1:
        raise ProjectError(16, f"{directory}: expected one main GTF, found {candidates or 'none'}")
    return f"{directory}/{candidates[0]}"


def _download_verified(task: Tuple[str, str]) -> Dict[str, object]:
    """Download ``url`` to ``dest`` unless a verified copy is already there."""
    url, dest = task
    path = Path(dest)
    directory = url.rsplit("/", 1)[0]
    if path.exists():
        try:
            return {"url": url, "path": str(path), "checksum": verify(path, directory), "reused": True}
        except ProjectError:
            pass  # stale or partial: download again
    _download(url, path)
    return {"url": url, "path": str(path), "checksum": verify(path, directory), "reused": False}


def _index_task(task: Tuple[str, str, List[Tuple[str, str]], str]) -> Tuple[str, Dict[str, int]]:
    """Process-pool worker: index the human GTF, or one species' homologies, annotation and NCBI orthologs."""
    what, species, raw_files, out = task
    if what == "human":
        return what, {"genes": index_human_gtf(Path(raw_files[0][1]), Path(out))}
    out_dir = Path(out)
    homology_files = [(k, Path(f)) for k, f in raw_files if k in HOMOLOGY_KINDS]
    counts = index_homologies(homology_files, species, out_dir)
    annotated = set()
    for kind, f in raw_files:
        if kind == "gtf":
            annotated = index_gtf(Path(f), out_dir / "species_genes.tsv", ("gene_id", "symbol", "biotype"))
    genes_path = out_dir / "genes.txt"
    homology_genes = {line for line in genes_path.read_text(encoding="utf-8").splitlines() if line}
    union = annotated | homology_genes
    genes_path.write_text("\n".join(sorted(union)) + "\n", encoding="utf-8")
    counts["annotated_genes"] = len(annotated)
    counts["genes_with_any_homology"] = len(homology_genes)
    counts["genes"] = len(union)
    ncbi = dict(raw_files)
    if "gene2ensembl" in ncbi:
        counts.update(index_ncbi_orthologs(Path(ncbi["gene2ensembl"]), Path(ncbi["gene_orthologs"]),
                                           ncbi["taxon"], out_dir / "ncbi_orthologs_to_human.tsv"))
    return species, counts


# ---------------------------------------------------------------------------
# NCBI Gene orthologs: the second source
# ---------------------------------------------------------------------------

#: Recorded for files whose provider publishes no checksum (NCBI, Reactome, GO, WikiPathways, MSigDB).
UNPUBLISHED_CHECKSUM = "md5 recorded at download (the source publishes none)"


def ncbi_base_url() -> str:
    return DEFAULT_NCBI_BASE_URL


def _download_ncbi(task: Tuple[str, str, bool]) -> Dict[str, object]:
    """Download an NCBI file unless a copy is already there (or ``refresh``).

    NCBI publishes no checksum for these files and rewrites them daily, so the
    local copy is the reference until the user asks for a newer one: its md5,
    size and download time are recorded for provenance. Integrity is checked
    when it is decompressed for indexing (a truncated gzip fails there).
    """
    url, dest, refresh = task
    path = Path(dest)
    reused = path.exists() and not refresh
    if not reused:
        _download(url, path)
    stamp = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")
    return {"url": url, "path": str(path), "reused": reused, "bytes": path.stat().st_size, "downloaded": stamp,
            "checksum": {"method": UNPUBLISHED_CHECKSUM, "value": md5_file(path)}}


def _ncbi_rows(path: Path, needed: Tuple[str, ...], keep=None):
    """Rows of an NCBI gene/DATA table as lists, with the column index; ``keep(line)`` pre-filters cheaply."""
    try:
        with _open_tsv(path) as handle:
            header = handle.readline().lstrip("#").rstrip("\n").split("\t")
            col = {name: i for i, name in enumerate(header)}
            missing = [n for n in needed if n not in col]
            if missing:
                raise ProjectError(17, f"{path.name}: no column {', '.join(missing)}; header is {header}")
            width = len(header)
            for line in handle:
                if keep is None or keep(line):
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < width:      # garbage that still decompressed: not an IndexError
                        raise EOFError(f"a row has {len(fields)} of {width} fields")
                    yield fields, col
    except (EOFError, OSError, zlib.error, UnicodeDecodeError) as exc:
        raise ProjectError(11, f"{path.name} is truncated or corrupt ({exc}); run `reference fetch` again "
                               "with --refresh-ncbi") from exc


def index_ncbi_orthologs(gene2ensembl: Path, gene_orthologs: Path, taxon: str, out: Path) -> Dict[str, int]:
    """Ensembl gene pairs (species, human) that NCBI Gene calls orthologs.

    NCBI's orthologs are between NCBI GeneIDs; ``gene2ensembl`` translates both
    sides to Ensembl gene IDs (one GeneID can match several Ensembl genes and
    vice versa -- every combination is kept, as Compara's one-to-many pairs
    are). Orientation is not assumed, as for Compara.
    """
    ensembl: Dict[str, Dict[str, set]] = {taxon: {}, HUMAN_TAXON: {}}
    prefixes = (taxon + "\t", HUMAN_TAXON + "\t")
    for f, col in _ncbi_rows(gene2ensembl, ("tax_id", "GeneID", "Ensembl_gene_identifier"),
                             keep=lambda line: line.startswith(prefixes)):
        tax, gene, ens = f[col["tax_id"]], f[col["GeneID"]], f[col["Ensembl_gene_identifier"]]
        if tax in ensembl and ens and ens != "-":
            ensembl[tax].setdefault(gene, set()).add(ens.split(".", 1)[0])
    pairs = set()
    wanted = {taxon, HUMAN_TAXON}
    for f, col in _ncbi_rows(gene_orthologs, ("tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"),
                             keep=lambda line: line.startswith(prefixes)):
        if f[col["relationship"]] != "Ortholog":
            continue
        a, b = (f[col["tax_id"]], f[col["GeneID"]]), (f[col["Other_tax_id"]], f[col["Other_GeneID"]])
        if {a[0], b[0]} != wanted:
            continue
        sp_gene, human_gene = (a[1], b[1]) if a[0] == taxon else (b[1], a[1])
        pairs.add((sp_gene, human_gene))
    rows = set()
    for sp_gene, human_gene in pairs:
        for sp_ens in ensembl[taxon].get(sp_gene, ()):
            for human_ens in ensembl[HUMAN_TAXON].get(human_gene, ()):
                rows.add((sp_ens, human_ens, sp_gene, human_gene))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(NCBI_ORTHOLOG_COLUMNS)
        writer.writerows(sorted(rows))
    return {"ncbi_ortholog_pairs_to_human": len(pairs), "ncbi_genes_with_ensembl_id": len(ensembl[taxon]),
            "ncbi_orthologs_to_human": len(rows)}


def taxonomy_ids(*, release: int = DEFAULT_RELEASE, base_url: str = DEFAULT_BASE_URL) -> Dict[str, str]:
    """Ensembl species name -> NCBI taxonomy ID, from the release's own species table."""
    body = _fetch_text(f"{base_url}/release-{release}/{SPECIES_TABLE}")
    if body is None:
        raise ProjectError(10, f"{base_url}/release-{release}/{SPECIES_TABLE}: not available")
    lines = [l for l in body.splitlines() if l.strip()]
    header = lines[0].lstrip("#").split("\t")
    col = {name: i for i, name in enumerate(header)}
    table = {}
    for line in lines[1:]:
        f = line.split("\t")
        if max(col["species"], col["taxonomy_id"]) < len(f):
            table[f[col["species"]].strip()] = f[col["taxonomy_id"]].strip()
    return table


def ensure_writable(directory: Path) -> Path:
    """Create ``directory`` if needed and prove it can be written, before any work starts.

    An uncreatable reference directory would otherwise surface as a page of
    pathlib traceback from inside a download thread. Here it is one numbered
    sentence, raised before anything is downloaded.
    """
    directory = Path(directory)
    try:
        directory.mkdir(parents=True, exist_ok=True)
        probe = directory / ".write_test"
        probe.write_text("ok")
        probe.unlink()
    except OSError as exc:
        raise ProjectError(15, f"{directory}: {exc.strerror or exc}") from exc
    return directory


def fetch(
    reference_dir: Path,
    species: Iterable[str],
    *,
    release: int = DEFAULT_RELEASE,
    base_url: str = DEFAULT_BASE_URL,
    threads: Optional[int] = None,
    refresh_ncbi: bool = False,
    ncbi_url: Optional[str] = None,
    log=print,
) -> Path:
    """Download, verify and index what ``species`` needs. Returns the manifest path.

    Downloads run in threads (at most parallel.MAX_DOWNLOADS at once); indexing,
    which is decompression-bound, runs one process per species. The
    manifest is written once, by this process, after all workers finish -- no
    worker ever writes it. Re-running resumes: a raw file that is already there
    and still verifies is reused, not downloaded again.
    """
    from .parallel import MAX_DOWNLOADS, map_ordered, resolve_threads

    budget = resolve_threads(threads)
    root = release_dir(reference_dir, release)
    ensure_writable(root)
    manifest_path = root / "reference.json"
    manifest = read_json(manifest_path) if manifest_path.exists() else {"release": release, "species": {}}
    manifest.setdefault("species", {})

    from .version import __version__

    wanted = []
    for sp in species:
        sp = sp.strip().lower()
        if sp in wanted:
            continue                    # listed twice
        if sp == HUMAN:
            log(f"{sp}: nothing to fetch")
            continue
        entry = manifest["species"].get(sp)
        # The index files must still be there as well, so a deleted index is rebuilt, not reported "already
        # indexed".
        index_present = all((root / "species" / sp / f).exists()
                            for f in ("genes.txt", "orthologs_to_human.tsv", "species_genes.tsv",
                                      "ncbi_orthologs_to_human.tsv"))
        if entry and entry.get("indexed_by") == __version__ and not refresh_ncbi and index_present:
            log(f"{sp}: already indexed by {__version__}")
            continue
        if entry:
            log(f"{sp}: indexed by {entry.get('indexed_by', 'an earlier version')}, "
                f"running {__version__} -- re-indexing (downloads are reused)")
        wanted.append(sp)

    # The human index too is rebuilt when it is missing from the manifest, from disk, or was built by another
    # version, so a deleted human_genes.tsv is rebuilt rather than left for the crosswalk to crash on.
    human_entry = manifest.get("human") or {}
    need_human = (not human_entry or not (root / "human_genes.tsv").exists()
                  or human_entry.get("indexed_by") != __version__)
    downloads: List[Tuple[str, str]] = []
    if need_human:
        url, name = human_gtf_url(base_url, release)
        downloads.append((url, str(root / "raw" / HUMAN / name)))
    if wanted:
        # The species-human pairs may sit in either species' files (see module
        # docstring), so the human homology files are always part of the set.
        for kind in HOMOLOGY_KINDS:
            url = homology_url(base_url, release, HUMAN, kind)
            downloads.append((url, str(root / "raw" / HUMAN / Path(url).name)))
    gtf_urls: Dict[str, str] = {}
    for sp in wanted:
        for kind in HOMOLOGY_KINDS:
            url = homology_url(base_url, release, sp, kind)
            downloads.append((url, str(root / "raw" / sp / Path(url).name)))
        gtf_urls[sp] = species_gtf_url(base_url, release, sp)
        downloads.append((gtf_urls[sp], str(root / "raw" / sp / Path(gtf_urls[sp]).name)))
    if not downloads:
        return manifest_path

    taxa: Dict[str, str] = {}
    ncbi_tasks: List[Tuple[str, str, bool]] = []
    if wanted:
        known = taxonomy_ids(release=release, base_url=base_url)
        for sp in wanted:
            if sp not in known:
                raise ProjectError(10, f"{sp}: not in Ensembl {release}'s species table; no taxonomy ID for NCBI")
            taxa[sp] = known[sp]
        ncbi_tasks = [(f"{ncbi_url or ncbi_base_url()}/{name}", str(root / "raw" / "ncbi" / name), refresh_ncbi)
                      for name in NCBI_FILES]

    log(f"downloading {len(downloads) + len(ncbi_tasks)} file(s), "
        f"{min(budget, MAX_DOWNLOADS, len(downloads))} at a time")
    fetched = {d["url"]: d for d in map_ordered(_download_verified, downloads,
                                                threads=min(budget, MAX_DOWNLOADS), kind="thread")}
    ncbi = {Path(d["path"]).name: d for d in map_ordered(_download_ncbi, ncbi_tasks,
                                                          threads=min(budget, MAX_DOWNLOADS), kind="thread")}
    for d in list(fetched.values()) + list(ncbi.values()):
        path = Path(str(d["path"]))
        check = "no published checksum; md5 recorded" if d["checksum"]["method"] == UNPUBLISHED_CHECKSUM \
            else f"{d['checksum']['method']} ok"
        log(f"  {'reused' if d['reused'] else 'fetched'} {path.parent.name}/{path.name} ({check})")

    tasks = []
    if need_human:
        url, name = human_gtf_url(base_url, release)
        tasks.append(("human", HUMAN, [("gtf", fetched[url]["path"])], str(root / "human_genes.tsv")))
    for sp in wanted:
        raw = [(kind, fetched[homology_url(base_url, release, name, kind)]["path"])
               for name in (sp, HUMAN) for kind in HOMOLOGY_KINDS]
        raw.append(("gtf", fetched[gtf_urls[sp]]["path"]))
        raw += [("gene2ensembl", ncbi["gene2ensembl.gz"]["path"]),
                ("gene_orthologs", ncbi["gene_orthologs.gz"]["path"]), ("taxon", taxa[sp])]
        tasks.append(("species", sp, raw, str(root / "species" / sp)))
    # The index files are rewritten in place, so their manifest entries go first: a run killed or failing
    # part-way must not leave truncated files under an entry saying this version had indexed them, for the
    # next fetch to keep as "already indexed".
    dropped = [sp for sp in wanted if manifest["species"].pop(sp, None) is not None]
    if need_human and manifest.pop("human", None) is not None:
        dropped.append(HUMAN)
    if dropped and manifest_path.exists():
        write_json(manifest_path, manifest)
    log(f"indexing {len(tasks)} item(s) with up to {min(budget, len(tasks))} process(es)")
    results = map_ordered(_index_task, tasks, threads=budget, kind="process")

    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for (what, sp, raw, _), (_, counts) in zip(tasks, results):
        if what == "human":
            d = fetched[human_gtf_url(base_url, release)[0]]
            manifest["human"] = {"url": d["url"], "checksum": d["checksum"], "genes": counts["genes"], "fetched": now,
                                 "indexed_by": __version__}
        else:
            files = {f"{name}/{kind}": {"url": fetched[homology_url(base_url, release, name, kind)]["url"],
                                        "checksum": fetched[homology_url(base_url, release, name, kind)]["checksum"]}
                     for name in (sp, HUMAN) for kind in HOMOLOGY_KINDS}
            files[f"{sp}/gtf"] = {"url": gtf_urls[sp], "checksum": fetched[gtf_urls[sp]]["checksum"]}
            for name, d in ncbi.items():
                files[f"ncbi/{name}"] = {k: d[k] for k in ("url", "checksum", "bytes", "downloaded")}
            manifest["species"][sp] = {"files": files, "counts": counts, "fetched": now,
                                       "taxonomy_id": taxa[sp], "indexed_by": __version__}
            log(f"{sp}: {counts['annotated_genes']} annotated genes ({counts['genes_with_any_homology']} with any "
                f"homology), {counts['orthologs_to_human']} ortholog rows to human")
            log(f"{sp}: NCBI (taxon {taxa[sp]}): {counts['ncbi_genes_with_ensembl_id']} NCBI genes with an Ensembl "
                f"ID, {counts['ncbi_orthologs_to_human']} Ensembl-to-human ortholog rows")
    write_json(manifest_path, manifest)
    return manifest_path


# ---------------------------------------------------------------------------
# which species exist, under what name, on which assembly
# ---------------------------------------------------------------------------

#: Ensembl's own per-release species table. Verified for
#: e116: tab-separated, header
#: "#name species division taxonomy_id assembly assembly_accession genebuild ..."
SPECIES_TABLE = "species_EnsemblVertebrates.txt"


def _list_directory(url: str) -> List[str]:
    """Entry names in an FTP-over-HTTP directory listing (or a local mirror)."""
    import re

    if url.startswith("file:"):
        from urllib.parse import urlparse
        from urllib.request import url2pathname

        local = Path(url2pathname(urlparse(url).path))
        return sorted(p.name for p in local.iterdir()) if local.is_dir() else []
    text = _fetch_text(url + "/")
    if text is None:
        raise ProjectError(10, f"{url}/: could not list the directory")
    return sorted({m.rstrip("/") for m in re.findall(r'href="([^"?/][^"?]*)"', text) if not m.startswith("..")})


def species_of_ids(gene_ids: Sequence[str], *, rest_url: Optional[str] = None, sample: int = 50) -> Dict[str, int]:
    """{Ensembl species name: IDs found} for up to ``sample`` IDs spread through the list, asked of Ensembl's
    REST service in one request (POST /lookup/id). Unknown or retired IDs are not counted."""
    import json

    from .identity import DISPLAY_NAME
    from .version import __version__

    ids = list(dict.fromkeys(gene_ids))
    if not ids:
        return {}
    k = min(sample, len(ids))
    picked = [ids[i * len(ids) // k] for i in range(k)]      # spread through the whole table
    url = f"{(rest_url or DEFAULT_REST_URL).rstrip('/')}/lookup/id"
    request = urllib.request.Request(
        url, data=json.dumps({"ids": picked}).encode("utf-8"), method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": f"{DISPLAY_NAME}/{__version__} (bioinformatics reference-data download)"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            answer = json.loads(response.read().decode("utf-8"))
    except (OSError, ValueError, http.client.HTTPException) as exc:
        raise ProjectError(10, f"{url}: {exc}") from exc
    counts: Dict[str, int] = {}
    for record in (answer or {}).values():
        if isinstance(record, dict) and record.get("species"):
            counts[str(record["species"])] = counts.get(str(record["species"]), 0) + 1
    return counts


def species_with_orthology(*, release: int = DEFAULT_RELEASE, base_url: Optional[str] = None) -> List[str]:
    return _list_directory(f"{base_url or DEFAULT_BASE_URL}/release-{release}/tsv/ensembl-compara/homologies")


def list_species(*, release: int = DEFAULT_RELEASE, base_url: str = DEFAULT_BASE_URL,
                 text: Optional[str] = None) -> List[Dict[str, str]]:
    """Every species in the release's table, marked with whether Compara has orthology for it.

    Only a species WITH orthology can be fetched: Ensembl lists many genomes
    (for sheep in e116: thirteen, mostly breed assemblies) but computes gene-tree
    homologies for one reference per species. ``text`` filters, case-insensitive,
    on common name, species name, assembly or accession -- so a user can search
    by the assembly they know ("Ramb") rather than guess a folder name.
    """
    body = _fetch_text(f"{base_url}/release-{release}/{SPECIES_TABLE}")
    if body is None:
        raise ProjectError(10, f"{base_url}/release-{release}/{SPECIES_TABLE}: not available")
    with_orthology = set(_list_directory(f"{base_url}/release-{release}/tsv/ensembl-compara/homologies"))
    lines = [l for l in body.splitlines() if l.strip()]
    header = lines[0].lstrip("#").split("\t")
    col = {name: i for i, name in enumerate(header)}
    needle = (text or "").lower()
    rows = []
    for line in lines[1:]:
        f = line.split("\t")
        get = lambda key: f[col[key]].strip() if key in col and col[key] < len(f) else ""
        row = {
            "species": get("species"), "name": get("name"), "assembly": get("assembly"),
            "accession": get("assembly_accession"),
            "orthology": "yes" if get("species") in with_orthology else "no",
        }
        if needle and not any(needle in v.lower() for v in row.values()):
            continue
        rows.append(row)
    rows.sort(key=lambda r: (r["orthology"] != "yes", r["name"].lower()))
    return rows


def load_manifest(reference_dir: Path, release: int = DEFAULT_RELEASE) -> Dict:
    path = release_dir(reference_dir, release) / "reference.json"
    if not path.exists():
        raise ProjectError(13, f"no reference data at {path.parent}")
    return read_json(path)
