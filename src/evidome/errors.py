"""Numbered errors that say how to fix themselves.

Paralome lesson: an error a user can quote ("PLM-0304") is an error that can be
looked up, searched for, and pinned by a test. An error that is
only prose gets paraphrased, and the paraphrase is not searchable.

Every error the program raises on purpose is registered here, once, with the
sentence that tells the user what to do about it. Codes are never reused and
never renumbered; a retired code stays in the registry, marked retired.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

from .identity import ERROR_PREFIX, env_var


@dataclass(frozen=True)
class ErrorSpec:
    code: int
    title: str
    fix: str
    retired: bool = False


#: code -> spec. Add new codes at the end; never reuse a number.
REGISTRY: Dict[int, ErrorSpec] = {}


def register(code: int, title: str, fix: str, *, retired: bool = False) -> ErrorSpec:
    if code in REGISTRY:
        raise ValueError(f"error code {code} is already registered")
    spec = ErrorSpec(code, title, fix, retired)
    REGISTRY[code] = spec
    return spec


def format_code(code: int) -> str:
    return f"{ERROR_PREFIX}-{code:04d}"


def _rebuild(code: int, message: str, context: Dict[str, Any]) -> "ProjectError":
    return ProjectError(code, message, **context)


class ProjectError(Exception):
    """An expected, explained failure. Carries its code and its fix."""

    def __init__(self, code: int, message: str, **context: Any) -> None:
        spec = REGISTRY.get(code)
        if spec is None:
            raise ValueError(f"unregistered error code {code}; add it to errors.REGISTRY")
        self.code = code
        self.spec = spec
        self.context = context
        self.message = message
        super().__init__(f"{format_code(code)} {spec.title}: {message}")

    def __reduce__(self):
        # Raised inside a worker process (parallel.map_ordered), the error is
        # pickled back to the parent; without this it arrives as a TypeError
        # traceback instead of its numbered sentence.
        return (_rebuild, (self.code, self.message, self.context))

    @property
    def fix(self) -> str:
        return self.spec.fix


def lookup(code: int) -> Optional[ErrorSpec]:
    return REGISTRY.get(code)


# ---------------------------------------------------------------------------
# The registry. 0001-0099: infrastructure (checkpoints, JSON, versions).
# ---------------------------------------------------------------------------

register(
    1,
    "Unreadable checkpoint or manifest",
    "The JSON file exists but cannot be read or is not valid JSON. If the message says permission "
    "denied, ask whoever wrote it to make it readable (chmod a+r). Otherwise delete that one file "
    "and re-run the same command; the stage will recompute.",
)
register(
    2,
    "Unstamped output",
    "A JSON artifact was about to be written without the program version that "
    "produced it. This is a programming error, not a user error -- report it.",
)
register(
    3,
    "Input lacks an Ensembl gene ID column",
    "Every input table must carry a gene ID column: Ensembl gene IDs (ENSG..., ENSMUSG..., ENSMMUG...), NCBI "
    "GeneIDs, or gene symbols (with --species). Name it with --id-column if it is not one of the recognised "
    "names.",
)

# 4-9: the input table (wide_stats.py).
register(
    4,
    "Duplicate gene ID",
    "Each gene must appear once. Version suffixes are stripped before this check, "
    "so ENSG00000139618.16 and ENSG00000139618.17 count as the same gene.",
)
register(
    5,
    "No comparisons found",
    "The table must contain at least one comparison block: {name}_P-value (with its "
    "sibling columns) or {name}_ANOVA_p. The input format is described in docs/INPUT_FORMAT.md, shipped "
    "with the package.",
)
register(
    6,
    "Incomplete comparison block",
    "A comparison is missing a required column. Directional comparisons (contrasts and burden "
    "tests) need P-value, Q-value and exactly one of Fold_change or Log2_fold_change; group means "
    "are optional, but a pair (Raw_FPKM_Control_mean + Raw_FPKM_Treatment_mean, or Control_mean + "
    "Treatment_mean) is given whole or not at all. F-tests need F, ANOVA_p and ANOVA_q.",
)
register(
    7,
    "Malformed row or value",
    "Every row must have as many fields as the header, and every statistic must be a "
    "finite number or a missing marker ('.' or empty); only a fold change may be infinite. "
    "The file must be UTF-8 text.",
)
register(
    8,
    "Unrecognised column",
    "Every column must be the gene ID, the symbol, or part of a comparison block. "
    "Remove the column or rename it to fit the grammar; unexplained columns are not "
    "passed through, because an agent would have to guess what they mean. A column name "
    "may appear only once: rename or remove the repeated copy.",
)
register(
    9,
    "Header and columns disagree",
    "The '#' comment header describes comparisons that do not match the columns. "
    "The file was probably edited by hand; regenerate it from its source.",
)

# 10 onwards: reference data, analysis stages, the cluster and figures, in the order they were added.
register(
    10,
    "Download failed",
    "Check the network from this machine and the URL. Reference data is fetched on "
    "a machine with internet access; re-running the same fetch command resumes where it stopped.",
)
register(
    11,
    "Checksum mismatch or unverifiable download",
    "The downloaded file does not match the directory's MD5SUM/CHECKSUMS, or neither "
    "lists it. Delete the file under raw/ and re-run the fetch. Unverified files are never used.",
)
register(
    12,
    "Unexpected homology file layout",
    "The Compara file lacks columns this reader relies on. The Ensembl release or file "
    "format differs from what was verified (e116); read the file's header before changing code.",
)
register(
    13,
    "No reference data",
    "If the message says where it looked and found nothing, run 'reference fetch "
    "--species <ensembl_species_name>' first. If it says there is no reference directory "
    "at all: on a Paralome machine an administrator runs 'paralome setup site' (the "
    "reference data then lives beside Paralome's genomes); otherwise set the "
    f"{env_var('REFERENCE')} environment variable or pass --reference-dir.",
)
register(
    14,
    "Species not identified",
    "Too few input gene IDs were found in any indexed species. The crosswalk asks Ensembl which "
    "species the IDs belong to and fetches it (unless --no-fetch-species); if that did not settle it, "
    "fetch the right species with 'reference fetch --species' (the prefix guess from 'inspect' is a "
    "hint), or lower --min-species-match if the input deliberately mixes species.",
)
register(
    15,
    "Reference directory cannot be created or written",
    f"Point --reference-dir (or the {env_var('REFERENCE')} environment variable) at a directory you "
    "can write to -- ideally on a filesystem every node sees, since later stages read "
    "it from compute nodes. Its parent must exist or be creatable by you.",
)
register(
    16,
    "Species annotation not found",
    "The species' GTF directory in this Ensembl release did not contain exactly one main "
    "annotation file. Check the species name with 'reference species', and read the "
    "directory listing before changing the file-name rule.",
)
register(
    17,
    "Unexpected NCBI Gene file layout",
    "gene2ensembl or gene_orthologs lacks a column this reader relies on. NCBI has changed "
    "the file; read its header line on the machine that fetched it before changing code.",
)
register(
    18,
    "Unexpected gene-set file layout",
    "A gene-set source file could not be read in the format verified for it (GMT, Reactome's "
    "6-column mapping, or GO's OBO). The provider has changed the file; read its first lines on "
    "the machine that fetched it before changing code, or re-fetch with --refresh. A damaged download "
    "(a truncated or HTML file saved under the expected name) is fixed by --refresh alone.",
)
register(
    19,
    "No gene-set index",
    "Run 'reference genesets' on a machine with network access first; it needs 'reference "
    "fetch' to have run for the same Ensembl release.",
)
register(
    20,
    "Fold_change convention cannot be read",
    "Retired: a column with negative values and values between -1 and 1 is read as a "
    "log fold change. Kept so the number is never reused.",
)
register(
    21,
    "Not a crosswalk output directory",
    "Point the command (enrich, regulators, report, figures, utility, database-value) at the --out "
    "directory of a 'crosswalk' run: it needs comparison_summary.json, "
    "ortholog_crosswalk.csv and comparisons/<name>.human_genes.csv.",
)
register(
    22,
    "Enrichment work incomplete or inconsistent",
    "Under SLURM: read the task logs in <out>/_slurm/logs for the comparison named. If the "
    "message names two versions, a node runs a different installed version -- install the "
    "same wheel everywhere (check with --version on each node) and submit again.",
)
register(
    23,
    "SLURM submission failed",
    "Run on a machine where sbatch works (the cluster's login node), check the partition, "
    "account and QoS options against 'sinfo' and 'sacctmgr', or use --backend local.",
)
register(
    24,
    "Evidence column without a description",
    "Every column in the evidence bundle must be described for the agent reading it. Add the "
    "column to report.COLUMN_TEXT (a program change, not a user error).",
)
register(
    25,
    "No regulator index",
    "Run 'reference regulators' on a machine with network access first; it needs 'reference "
    "fetch' to have run for the same Ensembl release.",
)
register(
    26,
    "Figure option not understood",
    "--colors takes exactly two colours, low then high, separated by a comma: names "
    "matplotlib knows (blue, firebrick) or #RRGGBB values ('#0044CC,#CC0000'). "
    "--figure-formats takes a comma list of pdf, svg, png, tiff; --figure-top at least 1; "
    "--dpi at least 72; --exclude a comma list of nc, zero, unsigned.",
)
register(
    27,
    "Figures need matplotlib",
    "Install it for the Python this program runs under (python3 -m pip install matplotlib), "
    "then run 'figures' on the output folder. The evidence bundle does not depend on it.",
)
register(
    28,
    "The database-value model needs numpy, scipy and scikit-learn",
    "Install them for the Python this program runs under (python3 -m pip install scikit-learn, which brings "
    "numpy and scipy), on every node that runs the jobs. The rest of the program does not depend on them.",
)
register(
    29,
    "Numeric option out of range",
    "--alpha greater than 0 and less than 1; --min-species-match above 0 and --min-ensembl-match at least 0, "
    "both at most 1; utility: --q and --jaccard above 0 and at most 1. GSEA: --permutations at least 1, --max-permutations at least "
    "--permutations, --min-size at least 1 and no more than --max-size, --min-exceedances at least 1. "
    "Regulators: --min-targets at least 1 and no more than --max-targets. database-value: --folds at least 2, "
    "--screen at least 1, --bootstrap and --permutations at least 0 (0 skips them).",
)
register(
    30,
    "Input table cannot be read",
    "Check the path: it must name an existing, readable table, or a folder of per-comparison tables where "
    "the command takes one (inspect, crosswalk, the whole run, combine).",
)
register(
    31,
    "Per-comparison tables cannot be combined",
    "Each file must be one comparison's gene table with a gene ID column (Gene, gene_id, Ensembl_ID, "
    "ensembl_gene_id, or R's unnamed row names) and either P-value (or ANOVA_p) with its sibling columns, or "
    "the columns of a DESeq2, edgeR, limma or Seurat FindMarkers results table (docs/INPUT_FORMAT.md). "
    "Comparison names come from the "
    "file names: use --pattern to choose the files and --drop to remove shared text from the names so each "
    "is unique and made of letters, digits and _ . + -.",
)
register(
    32,
    "Gene symbols need a species",
    "Gene symbols do not say which species they belong to (macaque and human symbols are nearly the same). "
    "Name the species with --species, as Ensembl names it, e.g. --species macaca_mulatta, --species "
    "mus_musculus or --species homo_sapiens. Ensembl gene IDs and NCBI GeneIDs need no --species.",
)
