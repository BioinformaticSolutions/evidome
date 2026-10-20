"""Which Ensembl release the reference data come from: decided by the data, never fixed in the code.

1.1. This program read every reference file from Ensembl 116. Data counted against an annotation Ensembl
has since replaced -- sheep Rambouillet v2.0, Ensembl 113, whose genes Ensembl 115 replaced with the
Rambouillet v3.0 gene set -- then matched no species at all. A release fixed in the code builds that failure
in for every annotation that is replaced later, so there is no release constant. Two releases are decided:

- The human release: the newest release in the reference directory that holds the gene sets; without one,
  the newest release Ensembl publishes (a step that downloads) or the newest on disk (a step that only reads).
  Gene sets and regulators are read at this release.
- The orthology release: the release of the annotation the input was counted against, read from the
  registered genome's GTF; the human release when the input matches no registered genome. The species'
  orthologs to human come from Compara of this release, and the human genes they reach are carried to the
  human release by stable ID (crosswalk.carry_to_human_release).

A release is an Ensembl release number (classic Ensembl, up to 116), or a dataset of Ensembl's new platform,
which numbers nothing: <assembly accession>_<provider>_<YYYY_MM>, as its FTP site keys it
(ensemblorganisms/<GCA>/NNN/NNN/NNN/<version>/<provider>/<YYYY_MM>/), with homologies by best reciprocal BLAST
hit in <dataset>/homology/<YYYY_MM_DD>/homology.tsv.gz.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

from .errors import ProjectError

Release = Union[int, str]

_RELEASE_DIR = re.compile(r"^ensembl-(\d+)$")
_FTP_RELEASE = re.compile(r"^release-(\d+)$")
#: Ensembl's classic GTF names: Species.Assembly.<release>.gtf, optionally .chr/.abinitio and .gz
_GTF_NAME = re.compile(r"^(?P<species>[A-Z][a-z]+(?:_[a-z0-9]+)+)\.(?P<assembly>.+)\.(?P<release>\d{2,3})"
                       r"(?:\.(?:chr|abinitio|chr_patch_hapl_scaff))?\.gtf(?:\.gz)?$")


def _number(release: Release) -> Release:
    text = str(release).strip()
    return int(text) if text.isdigit() else text


def releases_in(reference_dir: Path, required: str = "") -> List[int]:
    """Releases with a folder in the reference directory (and holding ``required``, when given), oldest first."""
    found = []
    root = Path(reference_dir)
    if root.is_dir():
        for entry in root.iterdir():
            match = _RELEASE_DIR.match(entry.name)
            if match and entry.is_dir() and (not required or (entry / required).exists()):
                found.append(int(match.group(1)))
    return sorted(found)


def newest_release(base_url: Optional[str] = None) -> int:
    """The newest release Ensembl publishes, read from its FTP site's listing (or a mirror's)."""
    from . import reference

    base = base_url or reference.DEFAULT_BASE_URL
    numbers = [int(m.group(1)) for name in reference._list_directory(base) for m in [_FTP_RELEASE.match(name)] if m]
    if not numbers:
        raise ProjectError(10, f"{base}/: no release-N folders listed, so the newest Ensembl release is not known; "
                               "give --release")
    return max(numbers)


#: What a reader is given when the reference directory holds no release yet: its own "not found" error then
#: names the folder, and nothing is downloaded to answer a question about what is on disk.
UNFETCHED = "unfetched"


def human_release(reference_dir: Optional[Path], explicit: Optional[Release] = None,
                  base_url: Optional[str] = None, *, offline: bool = True) -> Release:
    """The release gene sets, regulators and the human gene index are read at (see the module docstring).

    ``offline``: a reader (the default) never goes to the network; with nothing on disk it gets UNFETCHED. A
    step that downloads passes offline=False and gets the newest release Ensembl publishes."""
    if explicit is not None:
        return _number(explicit)
    # The human side is where the gene sets are. A release on disk only because a registered genome's
    # orthologs came from it (an older one) does not become the human side.
    present = releases_in(Path(reference_dir), "genesets/gene_sets.tsv") if reference_dir is not None else []
    if present:
        return present[-1]
    if not offline:
        return newest_release(base_url)
    for required in ("human_genes.tsv", ""):
        present = releases_in(Path(reference_dir), required) if reference_dir is not None else []
        if present:
            return present[-1]
    return UNFETCHED


#: Ensembl's new platform: .../ensemblorganisms/<GCA|GCF>/NNN/NNN/NNN/V/<provider>/<YYYY_MM>/geneset/genes.gtf.gz
_PLATFORM_GTF = re.compile(r"/(?P<provider>[A-Za-z][^/]*)/(?P<date>\d{4}_\d{2})/geneset/"
                           r"genes(?:-including_alt)?\.gtf(?:\.[bg]z)?$")
_ACCESSION = re.compile(r"^(?P<kind>GC[AF])_(?P<digits>\d{9})\.(?P<version>\d+)$")
#: The human reference assembly's accession, whatever its version.
HUMAN_ACCESSION = "GCA_000001405"
#: A dataset of Ensembl's new platform as a release identity: <accession>_<provider>_<YYYY_MM>.
_PLATFORM_RELEASE = re.compile(r"^(?P<accession>GC[AF]_\d{9}\.\d+)_(?P<provider>[A-Za-z][^_]*(?:_[A-Za-z][^_]*)*)_"
                               r"(?P<date>\d{4}_\d{2})$")


def _gtf_headers(gtf: Union[str, Path]) -> dict:
    """The GTF's ``#!key value`` header lines (Ensembl writes the assembly accession and annotation date)."""
    import gzip

    path = Path(str(gtf))
    headers = {}
    try:
        opener = gzip.open if path.read_bytes()[:2] == b"\x1f\x8b" else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                if line.startswith("#!"):
                    key, _, value = line[2:].strip().partition(" ")
                    headers[key] = value.strip()
    except OSError:
        pass
    return headers


def platform_release(gtf: Union[str, Path]) -> Optional[str]:
    """The new-platform dataset a GTF belongs to, as a release identity, or None when it is not one.

    Read from the path when the GTF keeps the FTP layout (.../<provider>/<YYYY_MM>/geneset/genes.gtf.gz), else
    from its headers: ``genome-build-accession`` and ``genebuild-last-updated`` (provider: ensembl)."""
    headers = _gtf_headers(gtf)
    accession = headers.get("genome-build-accession", "")
    match = _PLATFORM_GTF.search(str(gtf))
    if match:
        parts = Path(str(gtf)).parts
        try:
            i = parts.index("ensemblorganisms")
            kind, a, b, c, version = parts[i + 1:i + 6]
            accession = f"{kind}_{a}{b}{c}.{version}"
        except (ValueError, IndexError):
            pass
        if _ACCESSION.match(accession or ""):
            return f"{accession}_{match.group('provider')}_{match.group('date')}"
        return None
    if Path(str(gtf)).name.startswith("genes") and _ACCESSION.match(accession or ""):
        date = headers.get("genebuild-last-updated", "")[:7].replace("-", "_")
        if re.match(r"^\d{4}_\d{2}$", date):
            return f"{accession}_ensembl_{date}"
    return None


def is_platform(release: Release) -> bool:
    return bool(_PLATFORM_RELEASE.match(str(release)))


def platform_parts(release: Release) -> Tuple[str, str, str]:
    """(accession, provider, YYYY_MM) of a new-platform release identity."""
    match = _PLATFORM_RELEASE.match(str(release))
    if not match:
        raise ValueError(f"not a new-platform release: {release}")
    return match.group("accession"), match.group("provider"), match.group("date")


def platform_dataset_url(release: Release, base: str) -> str:
    accession, provider, date = platform_parts(release)
    m = _ACCESSION.match(accession)
    d = m.group("digits")
    return f"{base}/{m.group('kind')}/{d[0:3]}/{d[3:6]}/{d[6:9]}/{m.group('version')}/{provider}/{date}"


def release_of_gtf(gtf: Union[str, Path], genome: str = "") -> Release:
    """The Ensembl release a registered genome's annotation comes from, read from the GTF's own name."""
    match = _GTF_NAME.match(Path(str(gtf)).name)
    if match:
        return int(match.group("release"))
    platform = platform_release(gtf)
    if platform:
        return platform
    raise ProjectError(36, f"registered genome {genome or Path(str(gtf)).name}: its GTF {gtf} does not carry an "
                           "Ensembl release in its file name (Species.Assembly.<release>.gtf)")


def species_of_gtf(gtf: Union[str, Path]) -> Optional[str]:
    """Ensembl's species folder name for a classic GTF: Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.113.gtf ->
    ovis_aries_rambouillet. None when the name is not Ensembl's."""
    match = _GTF_NAME.match(Path(str(gtf)).name)
    return match.group("species").lower() if match else None
