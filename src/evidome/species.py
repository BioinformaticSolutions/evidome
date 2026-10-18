"""Which species (and later, which assembly) a set of gene IDs belongs to.

Two stages, deliberately:

1. PREFIX -- a fast first guess from the Ensembl stable-ID prefix (this module). Only a
   guess, for three reasons documented in Ensembl's own prefix table (e116):
   one prefix can cover many assemblies (ENSOAR: several sheep breeds; ENSMUS:
   mouse strains; ENSCAF: dog breeds), human has no species code (ENSG...),
   and non-vertebrate Ensembl genomes do not use ENS prefixes at all (octopus
   genes in Ensembl Metazoa look like ``Ocbimv22027483m.g``).
2. MEMBERSHIP -- the decision, made elsewhere: the input IDs are looked up in the
   registered genomes (``genomes.identify_registered_genome``) and in the
   downloaded Ensembl gene lists (``crosswalk.identify_species``), and the genome that
   holds (nearly) all of them wins, with the match rate reported.

This module's answer is labelled as a guess wherever it is shown (``inspect``); the
crosswalk replaces it with the membership result.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

#: Ensembl gene-ID prefix (species code + "G") -> species. Source: Ensembl
#: release 116 stable-ID prefix table. Deliberately short: a starter set for
#: the first-pass guess, not a registry. The membership stage replaces it.
GENE_PREFIXES: Dict[str, str] = {
    "ENSG": "Homo sapiens",
    "ENSMUSG": "Mus musculus",
    "ENSRNOG": "Rattus norvegicus",
    "ENSCAFG": "Canis lupus familiaris",
    "ENSOARG": "Ovis aries",
    "ENSMMUG": "Macaca mulatta",
    "ENSMFAG": "Macaca fascicularis",
    "ENSMNEG": "Macaca nemestrina",
    "ENSBTAG": "Bos taurus",
    "ENSSSCG": "Sus scrofa",
}

_ENSEMBL_GENE = re.compile(r"^(ENS[A-Z]*G)\d{11}$")


@dataclass
class PrefixGuess:
    species: Optional[str]
    prefix: Optional[str]
    matched: int
    total: int
    prefix_counts: List[Tuple[str, int]]

    @property
    def fraction(self) -> float:
        return self.matched / self.total if self.total else 0.0


def gene_prefix(gene_id: str) -> Optional[str]:
    match = _ENSEMBL_GENE.match(gene_id)
    return match.group(1) if match else None


def guess_species_by_prefix(gene_ids: Iterable[str]) -> PrefixGuess:
    counts: Counter = Counter()
    total = 0
    for gid in gene_ids:
        total += 1
        counts[gene_prefix(gid) or "(not an ENS gene ID)"] += 1
    ranked = counts.most_common()
    top_prefix, top_n = ranked[0] if ranked else (None, 0)
    species = GENE_PREFIXES.get(top_prefix) if top_prefix else None
    return PrefixGuess(species, top_prefix if species else None, top_n if species else 0, total, ranked)
