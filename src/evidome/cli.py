"""Command-line entry point.

Deliberately almost empty: the commands come out of the objectives discussion,
not ahead of it. What is settled already is only how the program identifies
itself -- ``--version`` reports the version of the code that is ACTUALLY
IMPORTED and the file it was imported from, next to the installed package
metadata, and says so loudly if they disagree. (Paralome's ``--version`` read
metadata only, which is the wrong question when you need to know what code a
compute node will run.)
"""

from __future__ import annotations

import argparse
import json
import sys
from importlib.metadata import PackageNotFoundError, version as installed_version
from pathlib import Path
from typing import List, Optional

from . import version as _version_module
from .identity import DISPLAY_NAME, NAME, env_var
from .version import __version__


def version_report() -> List[str]:
    lines = [f"{DISPLAY_NAME} {__version__}"]
    lines.append(f"  code imported from {Path(_version_module.__file__).parent}")
    lines.append(f"  python {sys.executable}")
    try:
        metadata = installed_version(NAME)
    except PackageNotFoundError:
        metadata = None
    if metadata is None:
        lines.append("  installed metadata: none (running from a source tree)")
    elif metadata != __version__:
        lines.append(
            f"  WARNING: installed metadata says {metadata}, but the imported code is "
            f"{__version__}. A stale copy is shadowing the install."
        )
    # AGPL section 13 asks a program with an interactive interface to show its
    # Appropriate Legal Notices; the command people run to ask what they are
    # running is where they belong. Same wording as Paralome, which is
    # licensed the same way.
    lines += [
        "",
        "Copyright (C) 2026 Kenneth L. Jones",
        "Licensed under the GNU Affero General Public License, version 3.",
        "This program comes with ABSOLUTELY NO WARRANTY.",
        "Full terms: <https://www.gnu.org/licenses/agpl-3.0.html>",
        "A separate commercial licence is available; see NOTICE.",
    ]
    return lines


def inspect_report(path: Path, id_column: Optional[str] = None, label: Optional[str] = None) -> List[str]:
    """Plain-text summary of an input table: what the program understood from it."""
    from .gene_ids import KIND_ENSEMBL, KIND_NCBI, id_kind
    from .species import guess_species_by_prefix
    from .wide_stats import read_wide_statistics

    table = read_wide_statistics(path, id_column=id_column)
    kind = id_kind(table.gene_ids)
    guess = guess_species_by_prefix(table.gene_ids)
    with_symbol = sum(s is not None for s in table.symbols)
    lines = [label or f"{path}"]
    lines.append(f"  produced by      {table.generator or '(not stated)'}")
    lines.append(f"  genes            {len(table.gene_ids)}  (ID column {table.id_column!r}; "
                 f"{table.versions_stripped} version suffixes stripped)")
    lines.append(f"  with a symbol    {with_symbol}  ({len(table.gene_ids) - with_symbol} without)")
    translated = "  -- translated to the species' Ensembl genes at the crosswalk; "
    lines.append(f"  gene IDs         {kind}s" + ("" if kind == KIND_ENSEMBL else translated + (
        "the species is found from them" if kind == KIND_NCBI else "needs --species")))
    if kind != KIND_ENSEMBL:
        pass                # a prefix says nothing about a symbol's or an NCBI GeneID's species
    elif guess.species:
        lines.append(f"  species (guess)  {guess.species}  -- prefix {guess.prefix} on "
                     f"{guess.matched}/{guess.total} IDs; not yet confirmed against reference gene lists")
    else:
        lines.append(f"  species (guess)  unknown prefix; most common: {guess.prefix_counts[:3]}")
    lines.append(f"  comparisons      {len(table.comparisons)}")
    for comp in table.comparisons:
        p = table.column(comp, comp.p_role)
        sig = sum(v is not None and v < 0.05 for v in p)
        if comp.kind in ("contrast", "burden"):
            extra = [r for r in ("Effect", "T-statistic", "Log2_fold_change", "Null_selection_probability")
                     if r in comp.columns]
            what = ("burden test: " if comp.kind == "burden" else "") + \
                f"{comp.treatment or '?'} vs {comp.control or '?'}" + (f"  [+{', '.join(extra)}]" if extra else "")
        elif comp.same_test_as:
            what = f"F-test = {comp.same_test_as} (skipped)"
        elif comp.direction_from:
            what = "F-test, direction from its contrasts"
        else:
            what = "F-test, no direction"
        lines.append(f"    {comp.name:<34} {what:<40} tested {comp.tested:>6}  p<0.05 {sig:>5}")
    for note in table.notes:
        lines.append(f"  note: {note}")
    return lines


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=NAME, usage=f"{NAME} TABLE --out DIR [options]   |   {NAME} COMMAND ...",
        description=f"{DISPLAY_NAME}: functional evidence reports for an interpreting agent. "
                    f"'{NAME} TABLE --out DIR' runs everything in order -- crosswalk, enrich, regulators, report "
                    f"(add --backend slurm to spread it over the cluster; '{NAME} TABLE --help' lists the options). "
                    "The commands below run one step each.")
    parser.add_argument("--version", action="store_true", help="report the version of the code actually running, and where it lives")
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    inspect = sub.add_parser("inspect", help="read an input statistics table and report what was understood from it")
    inspect.add_argument("table", type=Path, help="wide gene-statistics CSV (e.g. Paralome's all_genes_clean_statistics.csv), "
                                                  "one per-comparison table (also a DESeq2, edgeR, limma or Seurat "
                                                  "FindMarkers results table), or a folder of them; see "
                                                  "docs/INPUT_FORMAT.md in the package")
    inspect.add_argument("--id-column", help="name of the Ensembl gene ID column, if not one of the recognised names")
    _folder_options(inspect)

    comb = sub.add_parser("combine", help="weave a folder of per-comparison gene tables (one comparison per file, "
                                          "e.g. Paralome's single-cell results, or DESeq2, edgeR, limma or Seurat "
                                          "FindMarkers results tables) into one wide table")
    comb.add_argument("folder", type=Path, help="folder holding the per-comparison tables")
    comb.add_argument("--out", type=Path, required=True, help="the wide CSV to write")
    comb.add_argument("--id-column", help="name of the files' Ensembl gene ID column, if not one of the "
                                          "recognised names")
    _folder_options(comb)

    reference = sub.add_parser("reference", help="download and index reference data (run where there is internet access)")
    ref_sub = reference.add_subparsers(dest="reference_command")
    fetch = ref_sub.add_parser("fetch", help="download, verify and index Ensembl Compara orthology for species")
    fetch.add_argument("--species", action="append", required=True,
                       help="Ensembl species directory name, e.g. macaca_mulatta (repeatable)")
    _reference_options(fetch)
    fetch.add_argument("--base-url", default=None, help="mirror of https://ftp.ensembl.org/pub (testing)")
    fetch.add_argument("--refresh-ncbi", action="store_true",
                       help="download NCBI's gene2ensembl and gene_orthologs again (NCBI updates them daily; "
                            "by default the copy already downloaded is kept, for reproducibility)")
    fetch.add_argument("--ncbi-base-url", default=None,
                       help="mirror of https://ftp.ncbi.nlm.nih.gov/gene/DATA (testing)")

    genesets = ref_sub.add_parser("genesets", help="download and index the gene sets enrichment uses "
                                  "(Reactome, GO, WikiPathways, MSigDB Hallmark)")
    genesets.add_argument("--refresh", action="store_true",
                          help="download every source again (by default the copies already downloaded are kept)")
    genesets.add_argument("--msigdb-release", default=None, help="MSigDB human release (default 2026.1.Hs)")
    _reference_options(genesets)

    regs = ref_sub.add_parser("regulators", help="download and index upstream-regulator networks from OmniPath "
                              "(CollecTRI TF regulons, signed signalling), commercial-use records only")
    regs.add_argument("--refresh", action="store_true", help="download again (by default the copies are kept)")
    _reference_options(regs)

    signs = ref_sub.add_parser("pathway-signs", help="download Reactome's BioPAX and build activity patterns (which "
                               "genes are expected up or down when a pathway is active) for pathway activation z-scores")
    signs.add_argument("--refresh", action="store_true", help="download again (by default the copies are kept)")
    _reference_options(signs, threads=False)

    species = ref_sub.add_parser("species", help="list the species names --species accepts, with their assemblies")
    species.add_argument("text", nargs="?", help="filter: common name, species name, assembly or accession (e.g. sheep, Ramb)")
    species.add_argument("--all", action="store_true", help="also list genomes without orthology (they cannot be fetched)")
    species.add_argument("--release", type=int, default=116)
    species.add_argument("--base-url", default=None, help="mirror of https://ftp.ensembl.org/pub (testing)")

    crosswalk = sub.add_parser("crosswalk", help="identify the species of an input table and map every gene to human")
    crosswalk.add_argument("table", type=Path, help="wide gene-statistics CSV, one per-comparison table (also a DESeq2, "
                                                    "edgeR, limma or Seurat FindMarkers results table), or a folder of "
                                                    "them; see docs/INPUT_FORMAT.md in the package")
    crosswalk.add_argument("--out", type=Path, required=True, help="output directory")
    crosswalk.add_argument("--alpha", type=float, default=0.05,
                           help="raw p threshold used only to COUNT significant genes in the per-comparison "
                                "summary; nothing is filtered by it (default 0.05)")
    _crosswalk_options(crosswalk)
    _reference_options(crosswalk)

    enrich = sub.add_parser("enrich", help="gene-set enrichment (ORA and preranked GSEA) for every comparison "
                                           "of a crosswalk run")
    enrich.add_argument("crosswalk_dir", type=Path, help="the --out directory of a crosswalk run")
    enrich.add_argument("--alpha", type=float, default=None,
                        help="raw p threshold selecting genes for over-representation (default: the crosswalk "
                             "run's --alpha, so the counts and the selection agree)")
    _enrich_options(enrich)
    _backend_options(enrich)
    _reference_options(enrich)

    regulators = sub.add_parser("regulators", help="upstream-regulator evidence for every comparison of a crosswalk "
                                                   "run (overlap p, activation z, activity t)")
    regulators.add_argument("crosswalk_dir", type=Path, help="the --out directory of a crosswalk run")
    regulators.add_argument("--alpha", type=float, default=None,
                            help="raw p threshold for 'significant' targets (overlap, activation z; default: the "
                                 "crosswalk run's --alpha)")
    _regulator_options(regulators)
    _reference_options(regulators)

    # Not a listed command (there is no 'run' command): '<name> TABLE --out DIR' is routed here by main().
    run = sub.add_parser(PIPELINE, prog=f"{NAME} TABLE",
                         description="Everything, in order: crosswalk, enrich, regulators, report -- locally, or with "
                                     "--backend slurm as a chain of SLURM jobs.")
    run.add_argument("table", type=Path, help="wide gene-statistics CSV (e.g. Paralome's all_genes_clean_statistics.csv), "
                                              "one per-comparison table (also a DESeq2, edgeR, limma or Seurat "
                                              "FindMarkers results table), or a folder of them (see 'combine'); "
                                              "see docs/INPUT_FORMAT.md in the package")
    run.add_argument("--out", type=Path, required=True, help="output directory (the crosswalk folder)")
    run.add_argument("--alpha", type=float, default=0.05,
                     help="raw p threshold for 'significant' genes in every step (default 0.05); nothing is filtered")
    _crosswalk_options(run)
    _enrich_options(run)
    _regulator_options(run)
    _figure_options(run, pipeline=True)
    _backend_options(run)
    _reference_options(run)

    report = sub.add_parser("report", help="write the evidence bundle (CSV + manifest + zip) and a tabbed workbook "
                                           "for a crosswalk folder")
    report.add_argument("crosswalk_dir", type=Path, help="the --out directory of a crosswalk run (after enrich)")
    report.add_argument("--out", type=Path, default=None, help="evidence folder (default <crosswalk_dir>/evidence)")
    report.add_argument("--workbook", type=Path, default=None,
                        help="workbook path (default <crosswalk_dir>/<name>_workbook.xlsx)")

    figures = sub.add_parser("figures", help="draw the top-25 mirrored bar charts (regulators, pathways) for a "
                                             "crosswalk folder again, e.g. in other colours")
    figures.add_argument("crosswalk_dir", type=Path, help="the --out directory of a finished run")
    figures.add_argument("--out", type=Path, default=None, help="figure folder (default <crosswalk_dir>/figures)")
    _figure_options(figures, pipeline=False)
    figures.add_argument("--threads", type=int, default=0, help="CPU budget (default: all CPUs)")

    util = sub.add_parser("utility", help="score how much each database (gene-set source, regulator network) adds "
                                          "that the others do not, from a finished run's own results")
    util.add_argument("crosswalk_dir", type=Path, help="the --out directory of a finished run")
    util.add_argument("--q", type=float, default=0.05, help="a finding is a result with q below this (default 0.05)")
    util.add_argument("--jaccard", type=float, default=0.5,
                      help="two findings are near-duplicates when their significant genes overlap at this "
                           "Jaccard index or more (default 0.5)")
    util.add_argument("--out", type=Path, default=None, help="CSV path (default <crosswalk_dir>/source_utility.csv)")

    dv = sub.add_parser("database-value", help="how much of each comparison's gene-to-gene variation the databases "
                                               "explain together (goodness of fit) and what each one contributes, "
                                               "over every order of entry (Shapley); needs scikit-learn")
    dv.add_argument("crosswalk_dir", type=Path, help="the --out directory of a finished run")
    dv.add_argument("--out", type=Path, default=None, help="output folder (default <crosswalk_dir>/database_value)")
    dv.add_argument("--folds", type=int, default=5, help="cross-validation folds over genes (default 5)")
    dv.add_argument("--screen", type=int, default=300,
                    help="candidate columns kept per database inside each fold (default 300)")
    dv.add_argument("--bootstrap", type=int, default=50, help="bootstrap refits for feature stability (default 50)")
    dv.add_argument("--permutations", type=int, default=20,
                    help="shuffled-response refits for the goodness-of-fit null (default 20)")
    dv.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    dv.add_argument("--sources", nargs="+", default=None, choices=_gene_set_sources(),
                    help="gene-set sources (default all)")
    dv.add_argument("--networks", nargs="+", default=None, choices=("collectri", "signalling_via_tf"),
                    help="regulator networks (default both)")
    _figure_options(dv, pipeline=False, curves_only=True)
    _backend_options(dv, cpus=16, mem="32G")
    _reference_options(dv)
    dvt = sub.add_parser("database-value-task")    # run by SLURM, not listed
    dvt.add_argument("--plan", type=Path, required=True)
    dvt.add_argument("--index", type=int, required=True)
    dvt.add_argument("--threads", type=int, default=0)
    dvf = sub.add_parser("database-value-finalize")
    dvf.add_argument("--plan", type=Path, required=True)

    task = sub.add_parser("enrich-task")    # run by SLURM, not listed
    task.add_argument("--plan", type=Path, required=True)
    task.add_argument("--index", type=int, required=True)
    task.add_argument("--threads", type=int, default=0)
    final = sub.add_parser("enrich-finalize")
    final.add_argument("--plan", type=Path, required=True)
    finish = sub.add_parser("run-finish")   # run by SLURM, not listed
    finish.add_argument("--plan", type=Path, required=True)
    finish.add_argument("--threads", type=int, default=0)
    return parser


def _folder_options(parser: argparse.ArgumentParser) -> None:
    from .combine import DEFAULT_PATTERN

    parser.add_argument("--pattern", default=DEFAULT_PATTERN,
                        help="when TABLE is a folder: which files are the per-comparison tables, one glob or "
                             f"several separated by commas (default {DEFAULT_PATTERN})")
    parser.add_argument("--drop", action="append", default=[],
                        help="when TABLE is a folder: text to remove from the file names to make the comparison "
                             "names, repeatable (e.g. --drop _naive_cell_t_rna)")


def _crosswalk_options(parser: argparse.ArgumentParser) -> None:
    _folder_options(parser)
    parser.add_argument("--id-column")
    parser.add_argument("--species",
                        help="the input's species, as Ensembl names it (e.g. macaca_mulatta, mus_musculus, "
                             "homo_sapiens); needed when the gene IDs are gene symbols, which do not say which "
                             "species they belong to. Ensembl IDs and NCBI GeneIDs name their species themselves")
    parser.add_argument("--min-species-match", type=float, default=0.9,
                        help="fraction of input IDs that must be found in one genome (default 0.9)")
    parser.add_argument("--min-ensembl-match", type=float, default=0.5,
                        help="when a registered genome identified the input, the fraction that must also be in "
                             "Ensembl's annotation for orthology to be meaningful (default 0.5); below the lower of "
                             "this and --min-species-match, the input's species is looked up and fetched (see "
                             "--no-fetch-species)")
    parser.add_argument("--no-ncbi", action="store_true",
                        help="do not use NCBI Gene orthologs for genes Compara leaves without a human gene")
    parser.add_argument("--no-symbol-match", action="store_true",
                        help="do not match remaining genes to a human gene by an exact, unambiguous symbol")
    parser.add_argument("--keep-duplicate-ftests", action="store_true",
                        help="also analyse a main-effect F-test whose main-effect contrast is in the table (by "
                             "default it is skipped: with two levels per factor they are the same test)")
    parser.add_argument("--no-fetch-species", action="store_true",
                        help="do not download a species the reference directory lacks (by default Ensembl is asked "
                             "which species the IDs belong to and it is fetched once, as 'reference fetch' would)")
    parser.add_argument("--no-registry", action="store_true",
                        help="identify against Ensembl's annotation only, ignoring the machine's genome registry")


def _enrich_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-size", type=int, default=10, help="smallest gene set tested, counted within the "
                                                               "comparison's tested genes (default 10)")
    parser.add_argument("--max-size", type=int, default=500, help="largest gene set tested (default 500)")
    parser.add_argument("--permutations", type=int, default=1000, help="GSEA permutations per set size (default 1000)")
    parser.add_argument("--seed", type=int, default=1, help="random seed for GSEA permutations (default 1)")
    parser.add_argument("--max-permutations", type=int, default=100000,
                        help="adaptive GSEA: a set size whose sets are beaten by too few nulls gets 10x more "
                             "permutations, up to this many (default 100000; equal to --permutations turns it off)")
    parser.add_argument("--min-exceedances", type=int, default=10,
                        help="adaptive GSEA: extend while a set is beaten by fewer nulls than this (default 10)")
    # choices: a misspelt source would otherwise select no gene sets at all and the run would still succeed.
    parser.add_argument("--sources", nargs="+", default=None, choices=_gene_set_sources(),
                        help="gene-set sources to use (default all): reactome go_bp go_mf go_cc wikipathways hallmark")


def _regulator_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--min-targets", type=int, default=5,
                        help="fewest targets among tested genes for a regulator to be tested (default 5)")
    parser.add_argument("--max-targets", type=int, default=1000,
                        help="most targets for an indirect (signalling_via_tf) regulator (default 1000)")
    parser.add_argument("--networks", nargs="+", default=None, choices=("collectri", "signalling_via_tf"),
                        help="regulator networks to use (default both)")


def _gene_set_sources():
    from .genesets import SOURCES

    return SOURCES


def _figure_options(parser: argparse.ArgumentParser, *, pipeline: bool, curves_only: bool = False) -> None:
    """curves_only: database-value draws its own curves and reads only colours, formats and dpi; the bar-figure
    options are not offered there."""
    from .figures import DEFAULT_COLORS, DEFAULT_DPI, DEFAULT_FORMATS, DEFAULT_TOP

    parser.add_argument("--colors", default=DEFAULT_COLORS,
                        help="figures: two colours, low (inhibited/decreased) then high (activated/increased), "
                             f"comma-separated: names or #RRGGBB (default '{DEFAULT_COLORS}', blue and red)")
    if not curves_only:
        parser.add_argument("--figure-top", type=int, default=DEFAULT_TOP,
                            help=f"figures: rows on each side (default {DEFAULT_TOP})")
    parser.add_argument("--figure-formats", default=",".join(DEFAULT_FORMATS),
                        help=f"figures: any of pdf,svg,png,tiff (default {','.join(DEFAULT_FORMATS)})")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                        help=f"figures: resolution of png/tiff (default {DEFAULT_DPI})")
    if curves_only:
        return
    parser.add_argument("--pathway-score", default="activation", choices=("activation", "direction"),
                        help="figures: what colours a pathway -- activation (the activation z from Reactome's "
                             "activity patterns; pathways with none are NC, grey) or direction ((up - down) / sqrt(n) "
                             "for every set). Default activation; runs without patterns fall back to direction")
    parser.add_argument("--exclude", default="",
                        help="figures: leave out results with no sign -- nc (not calculable, grey), zero (a score "
                             "of 0, white) or unsigned (both); comma-separated. The next results by p take their "
                             "places. Default: none left out")
    if pipeline:
        parser.add_argument("--no-figures", action="store_true", help="do not draw figures")


def _figure_settings(args):
    """FigureSettings from the options, checked now -- before a long run, not at its end. None: --no-figures."""
    from .errors import ProjectError
    from .figures import FigureSettings, parse_colors, parse_exclude, parse_formats

    if getattr(args, "no_figures", False):
        return None
    from .figures import DEFAULT_TOP

    parse_colors(args.colors)
    top = getattr(args, "figure_top", DEFAULT_TOP)
    if top < 1 or args.dpi < 72:
        named = "--figure-top must be at least 1 and --dpi" if hasattr(args, "figure_top") else "--dpi must be"
        raise ProjectError(26, f"{named} at least 72")
    return FigureSettings(colors=args.colors, top=top, formats=parse_formats(args.figure_formats),
                          dpi=args.dpi, exclude=parse_exclude(getattr(args, "exclude", "")),
                          pathway_score=getattr(args, "pathway_score", "activation"))


def _backend_options(parser: argparse.ArgumentParser, *, cpus: int = 40, mem: str = "4G") -> None:
    """Defaults per command: database-value's are 16 CPUs and 32G. They are set here, at parsing, so its --help
    shows them and an explicit --cpus-per-task 40 or --mem 4G is never silently replaced."""
    parser.add_argument("--backend", choices=("local", "slurm"), default="local",
                        help="local: run here with --threads; slurm: one array task per comparison, then a "
                             "finalize job (default local)")
    parser.add_argument("--cpus-per-task", type=int, default=cpus,
                        help=f"slurm: CPUs per comparison task (default {cpus})")
    parser.add_argument("--mem", default=mem, help=f"slurm: memory per comparison task (default {mem})")
    parser.add_argument("--partition", default=None, help="slurm: partition")
    parser.add_argument("--account", default=None, help="slurm: account")
    parser.add_argument("--qos", default=None, help="slurm: quality of service")
    parser.add_argument("--nodelist", default=None, help="slurm: restrict to these nodes")
    parser.add_argument("--time", default=None, help="slurm: time limit per task, e.g. 2:00:00")


def _alpha(args) -> float:
    """--alpha, checked; enrich and regulators inherit the crosswalk run's when it is not given, so the
    report's stated alpha is the one ORA selected genes with."""
    from .errors import ProjectError

    if args.alpha is None:
        summary = Path(args.crosswalk_dir) / "comparison_summary.json"
        try:
            args.alpha = float(json.loads(summary.read_text(encoding="utf-8")).get("alpha", 0.05))
        except (OSError, ValueError, TypeError, AttributeError):
            args.alpha = 0.05
    if not 0 < args.alpha < 1:
        raise ProjectError(29, f"--alpha {args.alpha}")
    return args.alpha


def _check_targets(args) -> None:
    from .errors import ProjectError

    if args.min_targets < 1 or args.max_targets < args.min_targets:
        raise ProjectError(29, f"--min-targets {args.min_targets}, --max-targets {args.max_targets}")


def _settings(args):
    from .enrich import Settings
    from .errors import ProjectError

    _alpha(args)
    # --permutations 0 would give every GSEA row p = 1 and exit 0.
    if (args.permutations < 1 or args.max_permutations < args.permutations or args.min_size < 1
            or args.min_size > args.max_size or args.min_exceedances < 1):
        raise ProjectError(29, f"--permutations {args.permutations}, --max-permutations {args.max_permutations}, "
                               f"--min-size {args.min_size}, --max-size {args.max_size}, "
                               f"--min-exceedances {args.min_exceedances}")

    return Settings(args.alpha, args.min_size, args.max_size, args.permutations, args.seed,
                    args.max_permutations, args.min_exceedances)


def _slurm_options(args):
    from .slurm import SlurmOptions

    return SlurmOptions(args.cpus_per_task, args.mem, args.partition, args.account, args.qos, args.nodelist,
                        args.time)


def _reference_options(parser: argparse.ArgumentParser, *, threads: bool = True) -> None:
    parser.add_argument("--reference-dir", type=Path, default=None,
                        help="reference data directory (default: <machine shared root>/"
                             f"{NAME}, else ${env_var('REFERENCE')})")
    parser.add_argument("--release", type=int, default=116, help="Ensembl release (default 116)")
    if not threads:     # pathway-signs is one sequential parse; it takes no --threads
        return
    parser.add_argument("--threads", type=int, default=0,
                        help="TOTAL CPU budget for this command (default: all CPUs; 1 = sequential)")


def _reference_dir(args) -> Path:
    from .reference import resolve_reference_dir

    path, source = resolve_reference_dir(args.reference_dir)
    # Absolute: it is recorded in enrichment.json and regulators.json, where a relative path would make a later
    # `report` from another directory find no gene-set index, and local and SLURM runs would record it
    # differently, since the SLURM plan always resolves it.
    path = Path(path).resolve()
    print(f"reference data: {path}  (from {source})", file=sys.stderr)
    return path


def _require_indexes(reference_dir: Path, release: int) -> None:
    """The one-command run needs both indexes; say so before anything runs or is submitted.

    Looked for only in the run-finish step -- on a compute node, after enrichment -- a missing regulator index
    would let a SLURM run submit cleanly and fail later, and Paralome, which records EVD-0025 as 'not set up
    here', would see an ordinary failure instead."""
    from .errors import ProjectError
    from .genesets import genesets_dir
    from .regulators import regulators_dir

    if not (genesets_dir(reference_dir, release) / "gene_sets.tsv").exists():
        raise ProjectError(19, f"no gene-set index at {genesets_dir(reference_dir, release)}")
    if not (regulators_dir(reference_dir, release) / "regulons.tsv").exists():
        raise ProjectError(25, f"no regulator index at {regulators_dir(reference_dir, release)}")


def run_crosswalk(args, reference_dir: Optional[Path] = None) -> List[str]:
    from .crosswalk import build_crosswalk, summarize, write_crosswalk
    from .human_tables import build_human_tables, summary_lines, write_human_tables
    from .jsonio import write_json
    from .wide_stats import read_wide_statistics

    from .wide_stats import drop_duplicate_ftests

    from .errors import ProjectError

    if not 0 < args.alpha < 1 or not 0 < args.min_species_match <= 1 or not 0 <= args.min_ensembl_match <= 1:
        raise ProjectError(29, f"--alpha {args.alpha}, --min-species-match {args.min_species_match}, "
                               f"--min-ensembl-match {args.min_ensembl_match}")
    from .combine import table_or_folder

    table_path, folder, id_column = table_or_folder(args.table, args.out, pattern=args.pattern, drop=args.drop,
                                                    id_column=args.id_column,
                                                    log=lambda line: print(line, file=sys.stderr))
    table = read_wide_statistics(table_path, id_column=id_column)
    skipped = [] if args.keep_duplicate_ftests else drop_duplicate_ftests(table)
    reference_dir = reference_dir or _reference_dir(args)
    from .gene_ids import KIND_SYMBOL, normalise_species, translate_ids

    # Symbols and NCBI GeneIDs become the species' Ensembl genes here; Ensembl IDs pass through unchanged
    translation, gene_ids = translate_ids(table.gene_ids, table.input_gene_ids, reference_dir, release=args.release,
                                          species=args.species, fetch_species=not args.no_fetch_species,
                                          threads=args.threads, log=lambda line: print(line, file=sys.stderr))
    symbols = list(table.symbols)
    if translation.kind == KIND_SYMBOL:      # the given symbol is the gene's symbol, where the table had no other
        symbols = [s or g for s, g in zip(symbols, table.gene_ids)]
    table.gene_ids, table.symbols = gene_ids, symbols
    if not translation.translated and not args.no_fetch_species:
        from .crosswalk import ensure_species

        # the loosest threshold the crosswalk itself would accept: a table it can map never triggers a download
        ensure_species(table.gene_ids, reference_dir, release=args.release,
                       min_fraction=min(args.min_species_match, args.min_ensembl_match),
                       threads=args.threads, log=lambda line: print(line, file=sys.stderr))
    call, rows = build_crosswalk(table.gene_ids, table.symbols, reference_dir,
                                 release=args.release, min_fraction=args.min_species_match,
                                 min_ensembl_fraction=args.min_ensembl_match,
                                 use_registry=not args.no_registry, threads=args.threads,
                                 use_ncbi=not args.no_ncbi, use_symbols=not args.no_symbol_match,
                                 translation=translation)
    if args.species and normalise_species(args.species) != call.species:
        raise ProjectError(14, f"--species {args.species}, but the gene IDs are {call.species}'s "
                               f"({call.matched}/{call.total})")
    csv_path = write_crosswalk(rows, args.out / "ortholog_crosswalk.csv")
    summary = summarize(call, rows)
    summary["input_ids"] = translation.summary()
    summary.update({"input_table": str(Path(table_path).resolve()), "input_produced_by": table.generator,
                    **({"input_folder": str(Path(folder).resolve())} if folder else {}),
                    "reference": f"Ensembl {args.release} Compara"})
    json_path = write_json(args.out / "ortholog_crosswalk.json", summary)
    results = build_human_tables(table, rows, alpha=args.alpha, threads=args.threads)
    table_paths = write_human_tables(results, table, args.out)
    comparisons_path = write_json(args.out / "comparison_summary.json", {
        "input_table": str(Path(table_path).resolve()), **({"input_folder": str(Path(folder).resolve())} if folder else {}),
        "species": call.species,
        "genome": None if call.genome is None else call.genome.name,
        "alpha": args.alpha, "comparisons": [r.summary for r in results],
        "skipped_comparisons": [{"comparison": f, "kind": "ftest", "same_test_as": c} for f, c in skipped],
    })
    lines = []
    if call.genome is not None:
        g = call.genome
        lines.append(f"genome           {g.name}  ({g.matched}/{g.total} IDs, {g.fraction:.1%}; {g.source})")
        lines.append(f"  annotation     {g.gtf}")
    else:
        for note in call.notes:
            if not ("NCBI" in note or "symbol" in note):
                lines.append(f"genome           not from the registry: {note}")
    for note in (call.genome.notes if call.genome else []):
        lines.append(f"  note: {note}")
    for note in call.notes:
        if "NCBI" in note or "symbol" in note:
            lines.append(f"  note: {note}")
    if translation.translated:
        ids = translation.summary()
        lines.append(f"input IDs        {translation.kind}s, translated to {call.species} Ensembl genes: "
                     f"{ids['translated']} of {call.total}")
        for reason, n in ids["untranslated_by_reason"].items():
            lines.append(f"  {reason:<60} {n:>6}")
    lines.append(f"Ensembl species  {call.species}  ({call.matched}/{call.total} IDs in Ensembl {args.release}, "
                 f"{call.fraction:.1%})")
    for status, n in summary["input_genes_by_status"].items():
        lines.append(f"  {status:<34} {n:>6}")
    for source, n in summary["mapped_input_genes_by_source"].items():
        lines.append(f"    via {source:<30} {n:>6}")
    lines.append(f"human genes reached                {summary['human_genes_reached']:>6}")
    lines.append(f"  with many-to-one collisions      {summary['human_genes_with_many_to_one_collisions']:>6}")
    for kind, n in summary["pairs_by_homology_type"].items():
        lines.append(f"  pairs {kind:<28} {n:>6}")
    if summary["unmapped_input_genes_by_biotype"]:
        lines.append("unmapped input genes, by biotype")
        for kind, n in list(summary["unmapped_input_genes_by_biotype"].items())[:8]:
            lines.append(f"  {kind:<34} {n:>6}")
    lines += summary_lines(results)
    lines.append(f"wrote {csv_path}")
    lines.append(f"wrote {json_path}")
    lines.append(f"wrote {len(table_paths) // 2} human-gene tables in {args.out / 'comparisons'}, "
                 f"ballots in {args.out / 'collisions'}")
    for f, c in skipped:
        lines.append(f"skipped F-test {f}: the same test as {c}, which is analysed in its place")
    lines.append(f"wrote {comparisons_path}")
    return lines


def submitted_json(ids, keys) -> str:
    """One machine-readable line naming every submitted job.

    A caller that waits on the jobs -- Paralome's pipeline supervisor, which
    runs this program as its last stage -- reads ``{"submitted_job_ids": [...]}``
    from the output, the same form its own stages print. The human lines above
    stay for people.
    """
    import json

    return json.dumps({"submitted_job_ids": [str(ids[k]) for k in keys]})


def enrich_lines(manifest_path: Path) -> List[str]:
    """What ``enrich`` prints: per comparison, how many sets pass q < 0.05, by source."""
    from .jsonio import read_json

    manifest = read_json(manifest_path)
    lines = [f"gene sets with q < 0.05 (all tested sets are in the files; nothing is filtered)"]
    for c in manifest["comparisons"]:
        lines.append(f"  {c['comparison']}  (universe {c['universe']}, selected {c['selected']})")
        ora = c["ora_sets_q_below_0.05"]
        lines.append("    ORA   " + (", ".join(f"{k} {v}" for k, v in ora.items()) or "none"))
        if c.get("gsea"):
            gsea = c["gsea_sets_q_below_0.05"]
            lines.append("    GSEA  " + (", ".join(f"{k} {v}" for k, v in gsea.items()) or "none")
                         + f"   [ranked by {c['rank_metric']}]")
        else:
            lines.append(f"    GSEA  skipped: {c['gsea_skipped']}")
    lines.append(f"wrote {manifest_path.parent}")
    return lines


#: Internal name of the whole-pipeline parser; users never type it (``<name> TABLE --out DIR``).
PIPELINE = "__pipeline__"
COMMANDS = ("inspect", "combine", "reference", "crosswalk", "enrich", "regulators", "report", "figures", "utility",
            "database-value", "database-value-task", "database-value-finalize", "enrich-task",
            "enrich-finalize", "run-finish")


def _route(argv: List[str]) -> List[str]:
    """``<name> TABLE ...`` -> the pipeline parser; anything else is a command, an option, or help."""
    if argv and not argv[0].startswith("-") and argv[0] not in COMMANDS:
        return [PIPELINE] + argv
    # Options before TABLE (`<name> --threads 2 TABLE --out DIR`) belong to the pipeline too: the top level has
    # no options but --version and --help. Otherwise they would reach the command parser, whose error would
    # name '2' as an invalid command and list the internal ones.
    if argv and argv[0].startswith("-") and argv[0] not in ("-h", "--help", "--version"):
        return [PIPELINE] + argv
    return argv


def main(argv: Optional[List[str]] = None) -> int:
    from .errors import ProjectError

    argv = _route(list(sys.argv[1:] if argv is None else argv))
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.version:
        print("\n".join(version_report()))
        return 0
    try:
        if args.command == "inspect":
            from .combine import table_or_folder

            import tempfile

            with tempfile.TemporaryDirectory(prefix=f"{NAME}_inspect_") as scratch:
                table, folder, id_column = table_or_folder(args.table, Path(scratch), pattern=args.pattern,
                                                      drop=args.drop, id_column=args.id_column,
                                                      log=lambda line: print(line, file=sys.stderr))
                print("\n".join(inspect_report(table, id_column,
                                                label=f"{folder} (combined)" if folder else None)))
            return 0
        if args.command == "combine":
            from .combine import combine, find_tables

            combine(find_tables(args.folder, args.pattern, out=args.out), args.out, drop=args.drop,
                    id_column=args.id_column)
            print(f"next: {NAME} inspect {args.out}")
            return 0
        if args.command == "reference" and args.reference_command == "fetch":
            from .reference import DEFAULT_BASE_URL, fetch

            manifest = fetch(_reference_dir(args), args.species, release=args.release,
                             base_url=args.base_url or DEFAULT_BASE_URL, threads=args.threads,
                             refresh_ncbi=args.refresh_ncbi, ncbi_url=args.ncbi_base_url)
            print(f"reference manifest: {manifest}")
            return 0
        if args.command == "reference" and args.reference_command == "genesets":
            from .genesets import DEFAULT_MSIGDB_RELEASE, fetch_genesets

            path = fetch_genesets(_reference_dir(args), release=args.release,
                                  msigdb_release=args.msigdb_release or DEFAULT_MSIGDB_RELEASE,
                                  refresh=args.refresh, threads=args.threads)
            print(f"gene-set manifest: {path}")
            return 0
        if args.command == "reference" and args.reference_command == "species":
            from .reference import DEFAULT_BASE_URL, list_species

            rows = list_species(release=args.release, base_url=args.base_url or DEFAULT_BASE_URL, text=args.text)
            shown = [r for r in rows if args.all or r["orthology"] == "yes"]
            hidden = len(rows) - len(shown)
            print(f"{'--species':<34} {'common name':<30} {'assembly':<22} {'accession':<17} orthology")
            for r in shown:
                print(f"{r['species']:<34} {r['name'][:30]:<30} {r['assembly']:<22} {r['accession']:<17} {r['orthology']}")
            if not shown:
                print("(no species with orthology matched)")
            if hidden:
                print(f"{hidden} more genome(s) matched but have no Compara orthology and cannot be fetched; --all shows them")
            print(f"Ensembl release {args.release}. Pick the row whose assembly matches the annotation your data "
                  f"was aligned to; crosswalk checks the choice against your gene IDs.")
            return 0
        if args.command == "enrich":
            from .enrich import Settings, enrich, plan_slurm

            settings = _settings(args)
            if args.backend == "slurm":
                ids = plan_slurm(args.crosswalk_dir, _reference_dir(args), options=_slurm_options(args),
                                 release=args.release, settings=settings, sources=args.sources)
                print(f"array job {ids['array']}, finalize job {ids['finalize']}")
                print(submitted_json(ids, ("array", "finalize")))
                print(f"logs: {ids['logs']}")
                print(f"the summary is printed at the end of the finalize log")
                return 0
            path = enrich(args.crosswalk_dir, _reference_dir(args), release=args.release,
                          settings=settings, sources=args.sources, threads=args.threads,
                          log=lambda line: print(line, file=sys.stderr))
            print("\n".join(enrich_lines(path)))
            return 0
        if args.command == "reference" and args.reference_command == "pathway-signs":
            from .pathway_signs import fetch_pathway_signs

            path = fetch_pathway_signs(_reference_dir(args), release=args.release, refresh=args.refresh)
            print(f"pathway-sign manifest: {path}")
            return 0
        if args.command == "reference" and args.reference_command == "regulators":
            from .regulators import fetch_regulators

            path = fetch_regulators(_reference_dir(args), release=args.release, refresh=args.refresh,
                                    threads=args.threads)
            print(f"regulator manifest: {path}")
            return 0
        if args.command == "regulators":
            from .regulators import regulator_lines, run_regulators

            _check_targets(args)
            path = run_regulators(args.crosswalk_dir, _reference_dir(args), release=args.release,
                                  alpha=_alpha(args), min_targets=args.min_targets, max_targets=args.max_targets,
                                  networks=args.networks, threads=args.threads,
                                  log=lambda line: print(line, file=sys.stderr))
            print("\n".join(regulator_lines(path)))
            return 0
        if args.command == "report":
            from .report import build_report

            paths = build_report(args.crosswalk_dir, out=args.out, workbook=args.workbook)
            for key in ("evidence", "archive", "workbook"):
                print(f"wrote {paths[key]}")
            return 0
        if args.command in ("database-value", "database-value-task"):
            from .database_value import single_thread_maths

            single_thread_maths()          # before numpy loads
        if args.command == "database-value":
            from . import database_value as dvm

            # --folds 1 and --screen 0 would otherwise end in a raw scikit-learn ValueError.
            if args.folds < 2 or args.screen < 1 or args.bootstrap < 0 or args.permutations < 0:
                raise ProjectError(29, f"--folds {args.folds}, --screen {args.screen}, --bootstrap {args.bootstrap}, "
                                       f"--permutations {args.permutations}")
            settings = dvm.ValueSettings(folds=args.folds, screen=args.screen, bootstrap=args.bootstrap,
                                         permutations=args.permutations, seed=args.seed,
                                         **dvm.run_limits(args.crosswalk_dir))
            figure_settings = _figure_settings(args)
            log = lambda line: print(line, file=sys.stderr)  # noqa: E731
            if args.backend == "slurm":
                options = _slurm_options(args)
                ids = dvm.plan_slurm(args.crosswalk_dir, _reference_dir(args), options=options, release=args.release,
                                     out=args.out, settings=settings, sources=args.sources, networks=args.networks,
                                     figure_settings=figure_settings, log=log)
                print(f"array job {ids['array']}, finalize job {ids['finalize']}")
                print(submitted_json(ids, ("array", "finalize")))
                print(f"logs: {ids['logs']}")
                return 0
            print("\n".join(dvm.database_value(args.crosswalk_dir, _reference_dir(args), release=args.release,
                                               out=args.out, settings=settings, sources=args.sources,
                                               networks=args.networks, threads=args.threads,
                                               figure_settings=figure_settings, log=log)))
            return 0
        if args.command == "database-value-task":
            from .database_value import run_task

            entry = run_task(args.plan, args.index, threads=args.threads)
            print(f"{entry['fit']['comparison']}: held-out R2 {entry['fit']['heldout_r2_all']}")
            return 0
        if args.command == "database-value-finalize":
            from .database_value import run_finalize

            print("\n".join(run_finalize(args.plan)))
            return 0
        if args.command == "utility":
            from .utility import source_utility, utility_lines, write_utility

            if not (0 < args.q <= 1 and 0 < args.jaccard <= 1):
                raise ProjectError(29, f"--q {args.q}, --jaccard {args.jaccard}")
            rows = source_utility(args.crosswalk_dir, q_max=args.q, threshold=args.jaccard)
            print("\n".join(utility_lines(rows, args.q, args.jaccard)))
            print(f"wrote {write_utility(args.crosswalk_dir, rows, args.out)}  (per comparison and pooled)")
            return 0
        if args.command == "figures":
            from .figures import make_figures

            result = make_figures(args.crosswalk_dir, _figure_settings(args), out=args.out, threads=args.threads,
                                  required=True)
            print(f"wrote {result['figures']} figures in {result['dir']}")
            return 0
        if args.command == PIPELINE:
            from .pipeline import run_local, submit_chain

            # every option is checked before anything runs, so a bad GSEA option is not found only after the
            # crosswalk has already overwritten the folder
            figure_settings = _figure_settings(args)
            settings = _settings(args)
            _check_targets(args)
            reference = _reference_dir(args)
            _require_indexes(reference, args.release)
            log = lambda line: print(line, file=sys.stderr)  # noqa: E731
            print("\n".join(run_crosswalk(args, reference)))
            regulator_options = {"min_targets": args.min_targets, "max_targets": args.max_targets,
                                 "networks": args.networks}
            if args.backend == "slurm":
                ids = submit_chain(args.out, reference, options=_slurm_options(args), release=args.release,
                                   settings=settings, sources=args.sources,
                                   regulator_options=regulator_options, figure_settings=figure_settings, log=log)
                print(f"array job {ids['array']}, finalize job {ids['finalize']}, run-finish job {ids['finish']}")
                print(submitted_json(ids, ("array", "finalize", "finish")))
                print(f"logs: {ids['logs']}")
                print(f"when job {ids['finish']} ends: {args.out}/{args.out.resolve().name}_combined_results.zip and the workbook are written; "
                      f"its log ends with the summaries")
                return 0
            print("\n".join(run_local(args.out, reference, release=args.release, settings=settings,
                                       sources=args.sources, regulator_options=regulator_options,
                                       figure_settings=figure_settings, threads=args.threads, log=log)))
            return 0
        if args.command == "run-finish":
            from .pipeline import run_finish

            print("\n".join(run_finish(args.plan, threads=args.threads)))
            return 0
        if args.command == "enrich-task":
            from .enrich import run_slurm_task

            entry = run_slurm_task(args.plan, args.index, threads=args.threads)
            print(f"{entry['comparison']}: done")
            return 0
        if args.command == "enrich-finalize":
            from .enrich import run_slurm_finalize

            print("\n".join(enrich_lines(run_slurm_finalize(args.plan))))
            return 0
        if args.command == "crosswalk":
            print("\n".join(run_crosswalk(args)))
            return 0
    except ProjectError as err:
        print(str(err), file=sys.stderr)
        print(f"  fix: {err.fix}", file=sys.stderr)
        return 2
    parser.print_help()
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
