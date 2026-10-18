"""Species identification by membership and the input -> human crosswalk.

Orthology comes from up to three sources, strictly in order; a later
source is consulted only for a gene the earlier ones left without any human
gene, and ``homology_source`` names the source of every pair:

1. Ensembl Compara gene trees (``protein`` / ``ncrna``) -- typed, scored;
2. NCBI Gene orthologs (``ncbi``), translated to Ensembl IDs by NCBI's own
   gene2ensembl -- reaches many unnamed Ensembl gene models Compara left out;
3. an exact symbol match (``symbol``), case-insensitive, only when the symbol
   belongs to exactly one human gene -- recovers named multi-copy genes
   (histones, ...) that gene trees decline to call. The symbol is evidence
   here, never an identifier: the input is still keyed by its gene IDs.

The crosswalk has one row per (input gene, human ortholog) pair, and one row
for every input gene with no human ortholog -- nothing is dropped (classify,
don't exclude). Each row carries the evidence behind the mapping so
an agent can weigh it: ortholog type, Compara's high-confidence flag, percent
identity both ways, gene-order and whole-genome-alignment scores, and two
counts that expose ambiguity -- how many human genes this input gene maps to,
and how many input genes share this human gene (the many-to-one collisions
resolved later, per comparison).
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .gene_ids import Translation
from .genomes import GenomeCall, identify_registered_genome
from .reference import HUMAN, load_manifest, release_dir

CROSSWALK_COLUMNS = (
    "input_gene_id", "input_symbol", "input_biotype", "species",
    "human_gene_id", "human_symbol", "human_biotype",
    "homology_type", "is_high_confidence", "identity_pct", "human_identity_pct",
    "goc_score", "wga_coverage", "homology_source",
    "human_targets_of_input_gene", "input_genes_sharing_human_gene",
    "mapping_status", "input_id_given", "id_translation",
)

STATUS_MAPPED = "mapped"
SOURCE_NCBI = "ncbi"
SOURCE_SYMBOL = "symbol"
SOURCE_ORDER = ("protein", "ncrna", SOURCE_NCBI, SOURCE_SYMBOL)
STATUS_HUMAN = "human input (identity)"
STATUS_NO_ORTHOLOG = "no human ortholog in reference"
STATUS_NOT_IN_REFERENCE = "gene ID not in reference"
#: In the registered genome's annotation (the one the data was quantified
#: against) but absent from Ensembl's annotation of the species: the two
#: annotations differ. Not the same thing as "no human ortholog".
STATUS_NOT_IN_ENSEMBL = "in registered genome, not in Ensembl reference"


@dataclass
class SpeciesCall:
    species: str
    matched: int
    total: int
    ranking: List[Tuple[str, int]]
    genome: Optional[GenomeCall] = None
    notes: List[str] = field(default_factory=list)

    @property
    def fraction(self) -> float:
        return self.matched / self.total if self.total else 0.0


def _read_ids(path: Path) -> set:
    return {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}


def _human_genes(root: Path) -> Dict[str, Tuple[str, str]]:
    table = {}
    path = root / "human_genes.tsv"
    if not path.exists():
        raise ProjectError(13, f"no human gene index at {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            table[row["human_gene_id"]] = (row["human_symbol"], row["human_biotype"])
    return table


def identify_species(gene_ids: Sequence[str], reference_dir: Path, *, release: int = 116,
                     min_fraction: float = 0.9) -> SpeciesCall:
    """The indexed species whose gene list contains the most input IDs."""
    root = release_dir(reference_dir, release)
    manifest = load_manifest(reference_dir, release)
    ids = set(gene_ids)
    ranking = []
    if (root / "human_genes.tsv").exists():
        ranking.append((HUMAN, len(ids & set(_human_genes(root)))))
    for sp in manifest.get("species", {}):
        genes = root / "species" / sp / "genes.txt"
        if not genes.exists():          # one species' lost index must not fail every input
            raise ProjectError(13, f"{sp} is listed in {root / 'reference.json'} but {genes} is missing; run "
                                   f"'reference fetch --species {sp}' to rebuild it")
        ranking.append((sp, len(ids & _read_ids(genes))))
    ranking.sort(key=lambda item: (-item[1], item[0]))
    if not ranking:
        raise ProjectError(13, f"no species indexed under {root}")
    best, matched = ranking[0]
    call = SpeciesCall(best, matched, len(ids), ranking)
    if call.fraction < min_fraction:
        raise ProjectError(
            14,
            f"best match is {best} with {matched}/{len(ids)} IDs ({call.fraction:.1%}), below "
            f"{min_fraction:.0%}; indexed: {', '.join(f'{s} {n}' for s, n in ranking[:5])}",
        )
    return call


def ensure_species(gene_ids: Sequence[str], reference_dir: Path, *, release: int = 116, min_fraction: float = 0.5,
                   threads: Optional[int] = None, log=print) -> Optional[str]:
    """Download the input's species when no downloaded species holds its IDs, so a new species is
    found, converted and tested in one run.

    Ensembl's REST service names the species of a sample of the IDs (any species Ensembl has, not only the
    prefixes this program knows); a clear majority is fetched exactly as 'reference fetch --species' would,
    once, into the shared reference directory. Returns the species fetched, or None when nothing was needed.
    Needs network access; --no-fetch-species turns it off."""
    from . import reference
    from .reference import fetch, species_of_ids, species_with_orthology

    try:
        identify_species(gene_ids, reference_dir, release=release, min_fraction=min_fraction)
        return None
    except ProjectError as err:
        if err.code not in (13, 14) or "is missing" in str(err):
            raise
        original = err
    log("no downloaded species holds these gene IDs: asking Ensembl which species they belong to")
    try:
        found = species_of_ids(gene_ids, rest_url=reference.DEFAULT_REST_URL)
    except ProjectError as err:
        # No network here: the original answer stands, with the reason nothing was fetched
        raise ProjectError(original.code, f"{original.message}; "
                                          f"could not ask Ensembl which species this is ({err.message})") from None
    asked = min(len(set(gene_ids)), 50)
    if not found:
        raise ProjectError(14, f"Ensembl does not know any of {asked} sampled gene IDs (e.g. {gene_ids[0]}); "
                               "the input may not use Ensembl IDs")
    species, n = max(sorted(found.items()), key=lambda kv: kv[1])
    if n * 2 <= asked and len(found) == 1:
        raise ProjectError(14, f"Ensembl knows only {n} of {asked} sampled gene IDs ({species}); too few to "
                               "download a species for -- the input may not use Ensembl IDs")
    if n * 2 <= asked:
        raise ProjectError(14, "the sampled gene IDs belong to several species ("
                               + ", ".join(f"{s} {c}" for s, c in sorted(found.items(), key=lambda kv: -kv[1]))
                               + "); split the table by species")
    log(f"Ensembl: {n} of {asked} sampled IDs belong to {species}")
    manifest_path = release_dir(reference_dir, release) / "reference.json"
    if manifest_path.exists() and (species in load_manifest(reference_dir, release).get("species", {})
                                   or (species == HUMAN and (release_dir(reference_dir, release) /
                                                             "human_genes.tsv").exists())):
        # Already downloaded: the IDs match it poorly (another release or assembly), and fetching again would
        # change nothing
        raise ProjectError(original.code, f"{original.message}; Ensembl confirms the IDs are {species}, which is "
                                          "already downloaded -- the IDs may come from another release or "
                                          "assembly")
    if species != HUMAN:
        with_orthology = species_with_orthology(release=release, base_url=reference.DEFAULT_BASE_URL)
        if with_orthology and species not in with_orthology:
            raise ProjectError(14, f"{species} has no Compara orthology in Ensembl {release}; run 'reference "
                                   f"species' to find the assembly that has it, and give the table that genome's IDs")
    log(f"downloading the reference data for {species} (once; later runs reuse it)")
    fetch(reference_dir, [species], release=release, base_url=reference.DEFAULT_BASE_URL, threads=threads, log=log)
    return species


def build_crosswalk(gene_ids: Sequence[str], symbols: Sequence[Optional[str]], reference_dir: Path, *,
                    release: int = 116, min_fraction: float = 0.9, min_ensembl_fraction: float = 0.5,
                    use_registry: bool = True, registry: Optional[Path] = None,
                    threads: Optional[int] = None, use_ncbi: bool = True,
                    use_symbols: bool = True,
                    translation: Optional[Translation] = None) -> Tuple[SpeciesCall, List[Dict[str, object]]]:
    """Identify the genome (registered first) and species, then map every gene to human.

    With a registered genome holding >= ``min_fraction`` of the input IDs, that
    genome is the input's identity and Ensembl only has to supply orthology, so
    the Ensembl species needs only ``min_ensembl_fraction`` -- the shortfall is
    reported per gene as STATUS_NOT_IN_ENSEMBL, not treated as failure. Without
    one, Ensembl's annotation must itself hold ``min_fraction``.

    ``translation`` (gene_ids.translate_ids): when the input's IDs were symbols or NCBI GeneIDs, ``gene_ids`` are
    already the species' Ensembl genes (or the given ID, for a gene left untranslated), the species is the one
    the translation was made in, and the registry is not consulted. Every row records the ID as given and how
    it was translated.
    """
    if translation is not None and translation.translated:
        root = release_dir(reference_dir, release)
        known = set(_human_genes(root)) if translation.species == HUMAN else \
            _read_ids(root / "species" / str(translation.species) / "genes.txt")
        matched = sum(g in known for g in gene_ids)
        call = SpeciesCall(str(translation.species), matched, len(set(gene_ids)), [(str(translation.species), matched)])
        call.notes = [f"species from the {translation.kind}s' translation; registry not consulted"] + translation.notes
        genome, genome_genes = None, {}
    else:
        genome, genome_genes, notes = (None, {}, ["registry not consulted"])
        if use_registry:
            genome, genome_genes, notes = identify_registered_genome(
                gene_ids, reference_dir, min_fraction=min_fraction, threads=threads, registry=registry)
        call = identify_species(gene_ids, reference_dir, release=release,
                                min_fraction=min_ensembl_fraction if genome else min_fraction)
        call.genome = genome
        call.notes = list(notes)  # a copy: later notes must not also land in genome.notes
    root = release_dir(reference_dir, release)
    human = _human_genes(root)
    rows: List[Dict[str, object]] = []
    untranslated = translation.untranslated if translation is not None else {}

    def provenance(gid: str) -> Dict[str, str]:
        if translation is None:
            return {"input_id_given": gid, "id_translation": ""}
        return {"input_id_given": translation.given.get(gid, gid), "id_translation": translation.how.get(gid, "")}

    if call.species == HUMAN:
        for gid, sym in zip(gene_ids, symbols):
            known = gid in human
            hs, hb = human.get(gid, ("", ""))
            # The registered genome is what the data was quantified against: its biotype wins, and a
            # gene it holds that this Ensembl release lacks is "not in Ensembl", not "not in the reference".
            in_genome = gid in genome_genes
            ib = genome_genes[gid][1] if in_genome and genome_genes[gid][1] else hb
            rows.append({
                "input_gene_id": gid, "input_symbol": sym or "", "input_biotype": ib, "species": HUMAN,
                "human_gene_id": gid if known else "", "human_symbol": hs, "human_biotype": hb,
                "homology_type": "", "is_high_confidence": "", "identity_pct": "", "human_identity_pct": "",
                "goc_score": "", "wga_coverage": "", "homology_source": "",
                "human_targets_of_input_gene": 1 if known else 0,
                "input_genes_sharing_human_gene": 1 if known else 0,
                "mapping_status": (STATUS_HUMAN if known else untranslated.get(gid) or
                                   (STATUS_NOT_IN_ENSEMBL if in_genome else STATUS_NOT_IN_REFERENCE)),
                **provenance(gid),
            })
        return call, rows

    species_dir = root / "species" / call.species
    wanted = set(gene_ids)
    known_genes = _read_ids(species_dir / "genes.txt")
    biotype: Dict[str, str] = {}
    species_symbols: Dict[str, str] = {}
    annotation = species_dir / "species_genes.tsv"
    if annotation.exists():
        with annotation.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                biotype[row["gene_id"]] = row["biotype"]
                species_symbols[row["gene_id"]] = row.get("symbol", "")
    # The registered genome is the annotation the data was quantified against:
    # its biotypes win where both know a gene.
    for gid, (_, kind) in genome_genes.items():
        if kind:
            biotype[gid] = kind
    by_input: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    if not (species_dir / "orthologs_to_human.tsv").exists():
        raise ProjectError(13, f"{species_dir / 'orthologs_to_human.tsv'} is missing; run 'reference fetch "
                               f"--species {species_dir.name}' to rebuild it")
    with (species_dir / "orthologs_to_human.tsv").open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["gene_id"] in wanted:
                by_input[row["gene_id"]].append(row)

    # Source 2: NCBI Gene orthologs, only for genes Compara left without a human gene.
    ncbi_path = species_dir / "ncbi_orthologs_to_human.tsv"
    if not use_ncbi:
        call.notes.append("NCBI orthologs not consulted (--no-ncbi)")
    elif not ncbi_path.exists():
        call.notes.append(f"NCBI orthologs not indexed for {call.species}: run `reference fetch --species "
                          f"{call.species}` again with this version to add them")
    else:
        with ncbi_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gid = row["gene_id"]
                # gene2ensembl is many-to-many: one Ensembl gene under two NCBI GeneIDs gives the same Ensembl pair
                # twice in the index; one pair is one hit.
                if any(h["human_gene_id"] == row["human_gene_id"] for h in by_input.get(gid, ())):
                    continue
                if gid in wanted and row["human_gene_id"] in human and not _compara(by_input.get(gid)):
                    by_input[gid].append({
                        "gene_id": gid, "human_gene_id": row["human_gene_id"], "homology_type": "ortholog (NCBI)",
                        "is_high_confidence": "", "identity": "", "human_identity": "", "goc_score": "",
                        "wga_coverage": "", "source": SOURCE_NCBI,
                    })

    # Source 3: an exact, unambiguous symbol, only for genes still without a human gene.
    if use_symbols:
        by_symbol: Dict[str, List[str]] = defaultdict(list)
        for hid, (hsym, _) in human.items():
            if hsym:
                by_symbol[hsym.upper()].append(hid)
        ambiguous = 0
        for gid, sym in zip(gene_ids, symbols):
            if by_input.get(gid):
                continue
            label = sym or genome_genes.get(gid, ("", ""))[0] or species_symbols.get(gid, "")
            matches = by_symbol.get(label.upper(), []) if label else []
            if len(matches) == 1:
                by_input[gid].append({
                    "gene_id": gid, "human_gene_id": matches[0], "homology_type": "symbol match",
                    "is_high_confidence": "", "identity": "", "human_identity": "", "goc_score": "",
                    "wga_coverage": "", "source": SOURCE_SYMBOL,
                })
            elif len(matches) > 1:
                ambiguous += 1
        if ambiguous:
            call.notes.append(f"{ambiguous} input gene(s) have a symbol shared by several human genes; "
                              "not matched by symbol")
    else:
        call.notes.append("symbol matching not used (--no-symbol-match)")

    sharing: Dict[str, int] = defaultdict(int)
    for gid, hits in by_input.items():
        for human_id in {h["human_gene_id"] for h in hits}:
            sharing[human_id] += 1

    for gid, sym in zip(gene_ids, symbols):
        hits = sorted(by_input.get(gid, ()), key=lambda h: h["human_gene_id"])
        if not hits:
            if gid in untranslated:
                status = untranslated[gid]
            elif gid in known_genes:
                status = STATUS_NO_ORTHOLOG
            elif gid in genome_genes:
                status = STATUS_NOT_IN_ENSEMBL
            else:
                status = STATUS_NOT_IN_REFERENCE
            rows.append({c: "" for c in CROSSWALK_COLUMNS} | {
                "input_gene_id": gid, "input_symbol": sym or "", "input_biotype": biotype.get(gid, ""),
                "species": call.species,
                "human_targets_of_input_gene": 0, "input_genes_sharing_human_gene": 0, "mapping_status": status,
                **provenance(gid),
            })
            continue
        targets = len({h["human_gene_id"] for h in hits})
        for h in hits:
            hs, hb = human.get(h["human_gene_id"], ("", ""))
            rows.append({
                "input_gene_id": gid, "input_symbol": sym or "", "input_biotype": biotype.get(gid, ""),
                "species": call.species,
                "human_gene_id": h["human_gene_id"], "human_symbol": hs, "human_biotype": hb,
                "homology_type": h["homology_type"], "is_high_confidence": h["is_high_confidence"],
                "identity_pct": h["identity"], "human_identity_pct": h["human_identity"],
                "goc_score": _null(h["goc_score"]), "wga_coverage": _null(h["wga_coverage"]),
                "homology_source": h["source"],
                "human_targets_of_input_gene": targets,
                "input_genes_sharing_human_gene": sharing[h["human_gene_id"]],
                "mapping_status": STATUS_MAPPED,
                **provenance(gid),
            })
    return call, rows


def _compara(hits) -> bool:
    return any(h["source"] not in (SOURCE_NCBI, SOURCE_SYMBOL) for h in hits or ())


def _null(value: str) -> str:
    return "" if value in ("NULL", "\\N") else value


def write_crosswalk(rows: List[Dict[str, object]], path: Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CROSSWALK_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def summarize(call: SpeciesCall, rows: List[Dict[str, object]]) -> Dict[str, object]:
    inputs = {r["input_gene_id"] for r in rows}
    by_status: Dict[str, set] = defaultdict(set)
    types: Dict[str, int] = defaultdict(int)
    by_source: Dict[str, set] = defaultdict(set)
    for r in rows:
        if r["homology_source"]:
            by_source[str(r["homology_source"])].add(r["input_gene_id"])
    for r in rows:
        by_status[str(r["mapping_status"])].add(r["input_gene_id"])
        if r["homology_type"]:
            types[str(r["homology_type"])] += 1
    humans = {r["human_gene_id"] for r in rows if r["human_gene_id"]}
    unmapped_by_biotype: Dict[str, int] = defaultdict(int)
    for r in rows:
        if r["mapping_status"] != STATUS_MAPPED and r["mapping_status"] != STATUS_HUMAN:
            unmapped_by_biotype[str(r["input_biotype"]) or "(not annotated)"] += 1
    collisions = {r["human_gene_id"] for r in rows if int(r["input_genes_sharing_human_gene"] or 0) > 1}
    genome = call.genome
    return {
        "genome": None if genome is None else {
            "name": genome.name, "gtf": genome.gtf, "source": genome.source,
            "matched": genome.matched, "total": genome.total, "fraction": round(genome.fraction, 4),
            "withdrawn": genome.withdrawn, "ranking": genome.ranking[:5],
        },
        "identification_notes": call.notes,
        "species": call.species,
        "species_match": {"matched": call.matched, "total": call.total, "fraction": round(call.fraction, 4),
                          "ranking": call.ranking[:5]},
        "input_genes": len(inputs),
        "input_genes_by_status": {k: len(v) for k, v in sorted(by_status.items())},
        "mapped_input_genes_by_source": {k: len(by_source[k]) for k in SOURCE_ORDER if k in by_source},
        "human_genes_reached": len(humans),
        "human_genes_with_many_to_one_collisions": len(collisions),
        "pairs_by_homology_type": dict(sorted(types.items())),
        "unmapped_input_genes_by_biotype": dict(sorted(unmapped_by_biotype.items(), key=lambda kv: -kv[1])),
    }
