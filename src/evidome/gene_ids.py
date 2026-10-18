"""Input gene IDs that are not Ensembl IDs -- gene symbols or NCBI GeneIDs -- translated to the species'
Ensembl genes before the crosswalk.

Everything after the crosswalk works on Ensembl gene IDs, so a table keyed by another kind of ID is translated
first, gene by gene, and the rest of the run is unchanged. The kind is recognised from the ID column itself:

    NCBI GeneID      90% or more of the IDs are whole numbers (NCBI Gene's GeneID, "Entrez ID")
    Ensembl gene ID  at least half look like Ensembl stable gene IDs (ENS...G, FBgn, WBGene), or, with
                     --species, at least half are in that species' downloaded gene list (yeast's YAL001C)
    gene symbol      anything else

Symbols carry no species (macaque and human symbols are nearly the same), so a symbol table needs --species.
An NCBI GeneID belongs to exactly one taxon, so the species is found from the IDs, as it is for Ensembl IDs:
from the downloaded species' NCBI tables first, else from NCBI's gene2ensembl (downloaded once, shared).

Translation, per input ID, never guessed:

    symbol      the species' Ensembl annotation (its GTF's gene_name) -- an exact match first; failing that, a
                match ignoring case that is the only one. Current symbols only: an outdated symbol is not in the
                annotation and is left untranslated.
    NCBI ID     NCBI's gene2ensembl rows for the species' taxon.

An ID that names no gene, or more than one, is left untranslated, as are two input IDs that would become the
same Ensembl gene. An untranslated gene stays in every output under its given ID, with the reason as its
mapping status in the crosswalk; the crosswalk records, for every gene, the ID as given and how it was
translated.
"""

from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .errors import ProjectError

KIND_ENSEMBL = "Ensembl gene ID"
KIND_NCBI = "NCBI GeneID"
KIND_SYMBOL = "gene symbol"

#: How a gene's input ID became its Ensembl gene ID (crosswalk column id_translation).
HOW_ENSEMBL = "Ensembl gene ID (as given)"
HOW_SYMBOL_EXACT = "symbol, exact"
HOW_SYMBOL_CASE = "symbol, ignoring case"
HOW_NCBI = "NCBI GeneID (gene2ensembl)"

_ENSEMBL_LIKE = re.compile(r"^(ENS[A-Z]*G\d{6,}|FBgn\d+|WBGene\d+)$")
#: Share of the IDs that must be whole numbers for the table to be read as NCBI GeneIDs.
NCBI_SHARE = 0.9


@dataclass
class Translation:
    kind: str
    #: Set when the IDs were translated: the species given (symbols) or found from the IDs (NCBI GeneIDs).
    species: Optional[str] = None
    #: Per gene, keyed by the ID the rest of the run uses (the Ensembl ID, or the given ID if untranslated).
    given: Dict[str, str] = field(default_factory=dict)
    how: Dict[str, str] = field(default_factory=dict)
    #: Given ID -> why it was not translated (its crosswalk mapping status).
    untranslated: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    @property
    def translated(self) -> bool:
        return self.kind != KIND_ENSEMBL

    def summary(self) -> Dict[str, object]:
        reasons = Counter(self.untranslated.values())
        return {"kind": self.kind, "species": self.species,
                "translated": sum(1 for g in self.given if g not in self.untranslated) if self.translated else None,
                "untranslated": len(self.untranslated),
                "untranslated_by_reason": dict(sorted(reasons.items(), key=lambda kv: -kv[1])),
                "by_method": dict(Counter(h for g, h in self.how.items() if g not in self.untranslated))}


def id_kind(gene_ids: Sequence[str]) -> str:
    """What the IDs look like: KIND_NCBI, KIND_ENSEMBL or KIND_SYMBOL (membership can still turn a 'symbol'
    table into Ensembl IDs, see ``translate_ids``)."""
    ids = [g for g in gene_ids if g]
    if not ids:
        return KIND_ENSEMBL
    if sum(g.isdigit() for g in ids) >= NCBI_SHARE * len(ids):
        return KIND_NCBI
    if 2 * sum(bool(_ENSEMBL_LIKE.match(g)) for g in ids) >= len(ids):
        return KIND_ENSEMBL
    return KIND_SYMBOL


def normalise_species(name: str) -> str:
    """'Macaca mulatta' or 'macaca_mulatta' -> 'macaca_mulatta' (Ensembl's production name)."""
    return "_".join(name.strip().lower().split())


# ---------------------------------------------------------------------------
# the species' own gene lists
# ---------------------------------------------------------------------------

def _species_annotation(reference_dir: Path, release: int, species: str) -> List[Tuple[str, str]]:
    """[(Ensembl gene ID, symbol)] of the species, from its downloaded annotation."""
    from .reference import HUMAN, human_genes_path, release_dir

    if species == HUMAN:
        path, id_col, sym_col = human_genes_path(reference_dir, release), "human_gene_id", "human_symbol"
    else:
        path, id_col, sym_col = release_dir(reference_dir, release) / "species" / species / "species_genes.tsv", \
            "gene_id", "symbol"
    if not path.exists():
        raise ProjectError(13, f"{species}: no gene annotation at {path}; run 'reference fetch --species {species}'")
    with path.open(encoding="utf-8", newline="") as handle:
        return [(row[id_col], row.get(sym_col) or "") for row in csv.DictReader(handle, delimiter="\t")]


def _downloaded(reference_dir: Path, release: int, species: str) -> bool:
    from .reference import HUMAN, human_genes_path, load_manifest, release_dir

    if species == HUMAN:
        return human_genes_path(reference_dir, release).exists()
    if not (release_dir(reference_dir, release) / "reference.json").exists():
        return False
    return species in load_manifest(reference_dir, release).get("species", {}) and \
        (release_dir(reference_dir, release) / "species" / species / "species_genes.tsv").exists()


def _ensure_downloaded(reference_dir: Path, release: int, species: str, *, fetch_species: bool,
                       threads: Optional[int], log) -> None:
    if _downloaded(reference_dir, release, species):
        return
    if not fetch_species:
        raise ProjectError(13, f"{species} is not downloaded under {reference_dir} and --no-fetch-species is set; "
                               f"run 'reference fetch --species {species}'")
    from . import reference

    log(f"downloading the reference data for {species} (once; later runs reuse it)")
    reference.fetch(reference_dir, [species], release=release, base_url=reference.DEFAULT_BASE_URL,
                    threads=threads, log=log)


# ---------------------------------------------------------------------------
# NCBI GeneIDs
# ---------------------------------------------------------------------------

def _gene2ensembl(reference_dir: Path, release: int, *, fetch_species: bool, log) -> Path:
    """NCBI's gene2ensembl, the copy 'reference fetch' keeps (downloaded here if there is none)."""
    from . import reference

    path = reference.release_dir(reference_dir, release) / "raw" / "ncbi" / "gene2ensembl.gz"
    if path.exists():
        return path
    if not fetch_species:
        raise ProjectError(13, f"no copy of NCBI's gene2ensembl at {path} and --no-fetch-species is set; run "
                               "'reference fetch' for any species to download it")
    log("downloading NCBI's gene2ensembl (once; later runs reuse it)")
    reference._download_ncbi((f"{reference.ncbi_base_url()}/gene2ensembl.gz", str(path), False))
    return path


def _ncbi_to_ensembl(path: Path, taxa: Optional[Set[str]] = None,
                     genes: Optional[Set[str]] = None) -> Dict[Tuple[str, str], Set[str]]:
    """{(taxon, GeneID): {Ensembl gene IDs}} from gene2ensembl, for the given taxa and/or GeneIDs."""
    from .reference import _ncbi_rows

    prefixes = tuple(t + "\t" for t in taxa) if taxa else ()

    def keep(line: str) -> bool:
        # The file is tens of millions of lines: reject on the taxon and the GeneID before splitting a row
        if prefixes and not line.startswith(prefixes):
            return False
        return genes is None or line.split("\t", 2)[1] in genes

    table: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    for f, col in _ncbi_rows(path, ("tax_id", "GeneID", "Ensembl_gene_identifier"), keep=keep):
        gene, ens = f[col["GeneID"]], f[col["Ensembl_gene_identifier"]]
        if genes is not None and gene not in genes:
            continue
        if ens and ens != "-":
            table[(f[col["tax_id"]], gene)].add(ens.split(".", 1)[0])
    return table


def _taxon_of(reference_dir: Path, release: int, species: str) -> str:
    from . import reference

    if species == reference.HUMAN:
        return reference.HUMAN_TAXON
    if (reference.release_dir(reference_dir, release) / "reference.json").exists():
        entry = reference.load_manifest(reference_dir, release).get("species", {}).get(species)
        if entry and entry.get("taxonomy_id"):
            return str(entry["taxonomy_id"])
    taxa = reference.taxonomy_ids(release=release, base_url=reference.DEFAULT_BASE_URL)
    if species not in taxa:
        raise ProjectError(10, f"{species}: not in Ensembl {release}'s species table")
    return taxa[species]


def _species_of_taxon(reference_dir: Path, release: int, taxon: str) -> str:
    from . import reference

    if taxon == reference.HUMAN_TAXON:
        return reference.HUMAN
    if (reference.release_dir(reference_dir, release) / "reference.json").exists():
        for sp, entry in reference.load_manifest(reference_dir, release).get("species", {}).items():
            if str(entry.get("taxonomy_id")) == taxon:
                return sp
    for sp, tax in reference.taxonomy_ids(release=release, base_url=reference.DEFAULT_BASE_URL).items():
        if tax == taxon:
            return sp
    raise ProjectError(14, f"the NCBI GeneIDs belong to NCBI taxon {taxon}, which is not a species in Ensembl "
                           f"{release}")


def _ncbi_species(ids: Sequence[str], g2e: Path) -> Tuple[str, int, int]:
    """(taxon, IDs in it, IDs found in any taxon): the taxon most of the found IDs belong to."""
    found = _ncbi_to_ensembl(g2e, genes=set(ids))
    by_taxon = Counter(taxon for taxon, _ in found)
    if not by_taxon:
        raise ProjectError(14, f"none of the {len(set(ids))} NCBI GeneIDs (e.g. {ids[0]}) has an Ensembl gene in "
                               "NCBI's gene2ensembl; check the ID column")
    total = sum(by_taxon.values())
    taxon, n = by_taxon.most_common(1)[0]
    if 2 * n <= total:
        raise ProjectError(14, "the NCBI GeneIDs belong to several taxa ("
                               + ", ".join(f"{t} {c}" for t, c in by_taxon.most_common(5)) + "); split the table")
    return taxon, n, total


# ---------------------------------------------------------------------------
# the translation
# ---------------------------------------------------------------------------

def translate_ids(gene_ids: Sequence[str], raw_ids: Sequence[str], reference_dir: Path, *, release: int = 116,
                  species: Optional[str] = None, fetch_species: bool = True, threads: Optional[int] = None,
                  log=print) -> Tuple[Translation, List[str]]:
    """(the translation, the IDs the rest of the run uses -- one per input gene, in order).

    ``gene_ids`` are the table's IDs (Ensembl version suffixes stripped), ``raw_ids`` the same as written."""
    species = normalise_species(species) if species else None
    kind = id_kind(gene_ids)
    if kind == KIND_SYMBOL and species and _downloaded(reference_dir, release, species):
        known = {gid for gid, _ in _species_annotation(reference_dir, release, species)}
        if 2 * sum(g in known for g in gene_ids) >= len(gene_ids):
            kind = KIND_ENSEMBL         # the species' own stable IDs, in a form the pattern does not know
    if kind == KIND_ENSEMBL:
        return Translation(kind, None, dict(zip(gene_ids, raw_ids)), {g: HOW_ENSEMBL for g in gene_ids}), \
            list(gene_ids)

    if kind == KIND_SYMBOL:
        if not species:
            example = next((g for g in gene_ids if g), "")
            raise ProjectError(32, f"the gene IDs look like gene symbols (e.g. {example!r}), which do not say "
                                   "which species they belong to")
        _ensure_downloaded(reference_dir, release, species, fetch_species=fetch_species, threads=threads, log=log)
        candidates = _by_symbol(_species_annotation(reference_dir, release, species))
        mapped = {g: candidates(g) for g in gene_ids}
        what = "symbol"
    else:
        g2e = _gene2ensembl(reference_dir, release, fetch_species=fetch_species, log=log)
        if species:
            taxon = _taxon_of(reference_dir, release, species)
        else:
            taxon, n, total = _ncbi_species(gene_ids, g2e)
            species = _species_of_taxon(reference_dir, release, taxon)
            log(f"NCBI GeneIDs: {n} of the {total} found in gene2ensembl belong to {species} (taxon {taxon})")
        _ensure_downloaded(reference_dir, release, species, fetch_species=fetch_species, threads=threads, log=log)
        table = _ncbi_to_ensembl(g2e, taxa={taxon}, genes=set(gene_ids))
        mapped = {g: (sorted(table.get((taxon, g), ())), HOW_NCBI) for g in gene_ids}
        what = "NCBI GeneID"

    translation = Translation(kind, species)
    claims: Dict[str, List[str]] = defaultdict(list)
    for g in gene_ids:
        targets, _ = mapped[g]
        if len(targets) == 1:
            claims[targets[0]].append(g)
    ids: List[str] = []
    for g, raw in zip(gene_ids, raw_ids):
        targets, how = mapped[g]
        # The status is one of a few fixed sentences, so the crosswalk summary counts genes by reason; the
        # particulars (which genes) go in id_translation.
        if not targets:
            reason = (f"symbol not in the {species} annotation" if kind == KIND_SYMBOL
                      else f"NCBI GeneID with no {species} Ensembl gene")
            detail = ""
        elif len(targets) > 1:
            reason = f"{what} names several {species} Ensembl genes"
            detail = f"names {len(targets)}: {';'.join(targets[:10])}"
        elif len(claims[targets[0]]) > 1:
            reason = f"same {species} Ensembl gene as another input ID"
            others = [o for o in claims[targets[0]] if o != g]
            detail = f"{targets[0]}, also named by {';'.join(others[:10])}"
        else:
            ids.append(targets[0])
            translation.given[targets[0]] = raw
            translation.how[targets[0]] = how
            continue
        ids.append(g)
        translation.given[g] = raw
        translation.how[g] = detail
        translation.untranslated[g] = f"not translated: {reason}"
    done = len(gene_ids) - len(translation.untranslated)
    translation.notes.append(f"{kind}s translated to {species} Ensembl genes: {done} of {len(gene_ids)}")
    return translation, ids


def _by_symbol(annotation: Iterable[Tuple[str, str]]):
    """symbol -> ([Ensembl IDs], how): exact first, else ignoring case."""
    exact: Dict[str, List[str]] = defaultdict(list)
    folded: Dict[str, List[str]] = defaultdict(list)
    for gid, sym in annotation:
        if sym:
            exact[sym].append(gid)
            folded[sym.casefold()].append(gid)

    def candidates(symbol: str) -> Tuple[List[str], str]:
        if symbol in exact:
            return sorted(set(exact[symbol])), HOW_SYMBOL_EXACT
        return sorted(set(folded.get(symbol.casefold(), ()))), HOW_SYMBOL_CASE

    return candidates
