"""Many tables, one per comparison, woven into the one wide table this program reads.

Some producers (Paralome single-cell runs among them) write one full gene table per
comparison and no wide table; the same files may go to other tools, and this program must see
exactly the same statistics. Each input is one comparison::

    Gene,Common_name,Control_mean,Treatment_mean,Fold_change,P-value,Q-value,<one column per sample>

and becomes a block of the wide table, named after its file::

    Gene,Common_name,<name>_Control_mean,<name>_Treatment_mean,<name>_Fold_change,<name>_P-value,<name>_Q-value

Rules:

- Statistics are copied as text, never recomputed or reformatted.
- Every column whose name is a comparison role is kept (the means, Fold_change,
  P-value, Q-value, Effect, T-statistic; or F, ANOVA_p, ANOVA_q for an
  F-test). Anything else -- the per-sample columns -- is left out, because the
  wide table refuses columns it cannot place.
- A gene absent from one table is blank in that comparison's block.
- Rows with no gene ID ('.', 'NA' or empty) are left out and counted; the
  crosswalk maps by stable ID only.
- The comparison name is the file name without its extension, with each
  ``--drop`` text removed. It must be letters, digits and ``_ . + -``, and
  unique.
- A file's own '#' header is carried over: who wrote it, and its control and
  treatment when it names one contrast.
- ``--id-column`` names the ID column of the files; Paralome writes the gene ID
  as the symbol where there is none, and that is not kept as a symbol.
- Hidden files, the tables a run writes into its --out folder, and earlier
  combined tables are never inputs; an --out that is one of the inputs is
  refused.

Results tables written by the usual differential-expression tools are read as they come, one comparison
per file, by their own column names (``TOOL_FORMATS``): DESeq2 (``log2FoldChange``, ``pvalue``,
``padj``), edgeR (``logFC``, ``PValue``, ``FDR``), limma (``logFC``, ``P.Value``, ``adj.P.Val``, ``t``)
and Seurat's FindMarkers (``avg_log2FC``, ``p_val``, ``p_val_adj``). Only the named columns are
kept; the tool's other columns are left out like the per-sample ones. The gene ID column may be R's
row names: a blank first column name (``write.csv``), or a header line one field shorter than the
rows (``write.table``).

A single per-comparison file -- a tool's results table, or a table with unprefixed ``P-value`` /
``ANOVA_p`` columns -- given in place of the wide table is combined on its own the same way.

``<name> combine FOLDER --out wide.csv`` writes the table; a folder given in
place of the table to ``crosswalk`` or the one-command run is combined the
same way, into ``<out>/input_combined.csv``, and ``inspect`` combines it in a
temporary folder that is removed afterwards.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .wide_stats import (CONTRAST_ROLES, FTEST_REQUIRED, ID_COLUMN_NAMES, MISSING_TOKENS, SYMBOL_COLUMN_NAMES,
                         _parse_header, delimiter_of, is_directional, strip_version)

#: Results tables of the usual differential-expression tools, recognised by their column names: (tool,
#: {its column: the role it fills}). A file matches a tool when it has every listed column and none of this
#: program's own role columns. The fold changes are all log2. Other columns of the tool (baseMean, lfcSE,
#: stat, logCPM, LR, F, AveExpr, B, pct.1, pct.2, ...) are not carried.
TOOL_FORMATS: Tuple[Tuple[str, Dict[str, str]], ...] = (
    ("DESeq2", {"log2FoldChange": "Log2_fold_change", "pvalue": "P-value", "padj": "Q-value"}),
    ("edgeR", {"logFC": "Log2_fold_change", "PValue": "P-value", "FDR": "Q-value"}),
    ("limma", {"logFC": "Log2_fold_change", "P.Value": "P-value", "adj.P.Val": "Q-value", "t": "T-statistic"}),
    ("Seurat FindMarkers", {"avg_log2FC": "Log2_fold_change", "p_val": "P-value", "p_val_adj": "Q-value"}),
)
#: R's row names: the ID column of a table written by write.csv has no name.
ROW_NAMES_COLUMN = ""

COMBINED_NAME = "input_combined.csv"
#: Which files of a folder are the tables: comma-separated globs (--pattern).
DEFAULT_PATTERN = "*.csv,*.tsv,*.txt"
#: Tables a run writes at the top of its --out folder: never inputs, should --out be the input folder itself.
RUN_OUTPUTS = frozenset({COMBINED_NAME, "ortholog_crosswalk.csv", "source_utility.csv"})
_NAME = re.compile(r"^[A-Za-z0-9_.+-]+$")


def comparison_name(path: Path, drop: Sequence[str] = ()) -> str:
    name = Path(path).stem
    for text in drop:
        name = name.replace(text, "")
    return name.strip("_")


def _written_by_us(path: Path) -> bool:
    """A table this program wrote (a previous combine): its first line names it."""
    from .identity import DISPLAY_NAME

    try:
        with path.open(encoding="utf-8-sig", errors="replace") as handle:
            first = handle.readline()
    except OSError:
        return False
    return first.startswith("# generated by") and f"{DISPLAY_NAME} " in first and " combine" in first


def find_tables(folder: Path, pattern: str = DEFAULT_PATTERN, out: Optional[Path] = None) -> List[Path]:
    """The per-comparison tables in ``folder``: matching ``pattern`` (one glob, or several separated by commas),
    not hidden (macOS writes ._name copies on shared storage), not a table this program wrote. ``out``, the file
    about to be written, must not be one of them: that would overwrite an input."""
    folder = Path(folder)
    if not folder.is_dir():
        raise ProjectError(30, f"cannot read {folder}: no such folder")
    globs = [g.strip() for g in pattern.split(",") if g.strip()]
    files = sorted({p for g in globs for p in folder.glob(g) if p.is_file() and not p.name.startswith(".")
                    and p.name not in RUN_OUTPUTS and not _written_by_us(p)})
    if out is not None:
        target = Path(out).resolve()
        if any(p.resolve() == target for p in files):
            raise ProjectError(31, f"--out {out} is one of the input tables; writing it would destroy that "
                                   "comparison -- write the combined table somewhere else")
    if not files:
        raise ProjectError(31, f"no file matching {pattern!r} in {folder}")
    return files


def tool_of(header: Sequence[str]) -> Optional[Tuple[str, Dict[str, str]]]:
    """(tool, {its column: role}) when ``header`` is one of TOOL_FORMATS' results tables, else None. A file that
    has this program's own role columns is read by those, whatever else it holds."""
    if "P-value" in header or "ANOVA_p" in header:
        return None
    for tool, columns in TOOL_FORMATS:
        if all(c in header for c in columns):
            return tool, columns
    return None


def _header_and_rows(lines: Sequence[str]) -> Tuple[List[str], List[List[str]]]:
    """The column names and data rows of a table's non-comment lines. A header one field shorter than every
    row is R's write.table with row names: the missing first name is the row-name (gene ID) column."""
    reader = csv.reader(lines, delimiter=delimiter_of(lines[0]) if lines else ",")
    header = [c.strip() for c in next(reader, [])]
    rows = list(reader)
    if header and rows and all(len(r) == len(header) + 1 for r in rows):
        header = [ROW_NAMES_COLUMN] + header
    return header, rows


def is_single_comparison(path: Path) -> bool:
    """Whether the file ``path`` holds one comparison's table rather than the wide table: a tool's results
    table, or unprefixed P-value / ANOVA_p columns."""
    try:
        with Path(path).open(newline="", encoding="utf-8-sig") as handle:
            lines = [line for line in handle.read().splitlines() if line.strip() and not line.startswith("#")]
    except (OSError, UnicodeDecodeError):
        return False     # the wide-table reader reports what is wrong with it
    header, _ = _header_and_rows(lines[:2])
    return tool_of(header) is not None or "P-value" in header or "ANOVA_p" in header


class _Table:
    def __init__(self, path: Path, name: str):
        self.path, self.name = path, name
        self.roles: List[str] = []
        self.rows: Dict[str, List[str]] = {}
        self.symbols: Dict[str, str] = {}
        self.ids: Dict[str, str] = {}
        self.no_id = 0
        self.generator: Optional[str] = None
        self.arms: Optional[Tuple[str, str]] = None
        self.test: Optional[str] = None     # 'burden' when the file's header says test=burden
        self.tool: Optional[str] = None     # the tool whose results table this is (TOOL_FORMATS)


def _read(path: Path, name: str, id_column: Optional[str] = None) -> _Table:
    table = _Table(path, name)
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            raw = handle.read().splitlines()
    except UnicodeDecodeError:
        raise ProjectError(7, f"{path} is not UTF-8 text; save it as UTF-8") from None
    except OSError as exc:
        raise ProjectError(30, f"cannot read {path}: {exc.strerror or exc}") from None
    comments = [line for line in raw if line.startswith("#")]
    lines = [line for line in raw if line.strip() and not line.startswith("#")]
    # What the file's own header says survives the combining: who wrote it, and which arm is control.
    table.generator, declared = _parse_header(comments)
    arms = [(v, kind) for (_, kind), v in declared.items() if is_directional(kind) and v[0]]
    table.arms = arms[0][0] if len(arms) == 1 else None
    table.test = arms[0][1] if len(arms) == 1 and arms[0][1] != "contrast" else None
    header, data = _header_and_rows(lines)
    repeated = sorted({c for c in header if header.count(c) > 1})
    if repeated:
        raise ProjectError(8, f"{path.name}: column name(s) appear more than once: {', '.join(repeated)}")
    if id_column is not None:
        if id_column not in header:
            raise ProjectError(3, f"{path.name}: requested ID column {id_column!r} is not in the file")
        id_col = id_column
    else:
        id_col = next((c for c in ID_COLUMN_NAMES if c in header), None)
        if id_col is None and header and header[0] == ROW_NAMES_COLUMN:
            id_col = ROW_NAMES_COLUMN
    if id_col is None:
        raise ProjectError(31, f"{path.name}: none of {', '.join(ID_COLUMN_NAMES)} is a column and the first column "
                               "is not R's unnamed row names (name the ID column with --id-column)")
    symbol_col = next((c for c in SYMBOL_COLUMN_NAMES if c in header and c != id_col), None)
    tool = tool_of(header)
    if tool is not None:
        table.tool, source = tool[0], {role: column for column, role in tool[1].items()}
        table.roles = [r for r in CONTRAST_ROLES if r in source]
    else:
        # An F-test's Direction_estimate is derived by the reader, never taken from a file
        roles = FTEST_REQUIRED if "ANOVA_p" in header and "P-value" not in header else CONTRAST_ROLES
        table.roles = [r for r in roles if r in header]
        source = {r: r for r in table.roles}
    if not any(r in table.roles for r in ("P-value", "ANOVA_p")):
        known = "; ".join(f"{name}: {', '.join(columns)}" for name, columns in TOOL_FORMATS)
        raise ProjectError(31, f"{path.name}: no P-value or ANOVA_p column, and not a results table of a known "
                               f"tool ({known}); is this a per-comparison gene table?")
    at = {r: header.index(source[r]) for r in table.roles}
    at.update({c: header.index(c) for c in [id_col] + ([symbol_col] if symbol_col else [])})
    for row_no, row in enumerate(data, start=1):
        if len(row) != len(header):
            raise ProjectError(7, f"{path.name}, data row {row_no}: {len(row)} fields, header has {len(header)}")
        raw_id = row[at[id_col]].strip()
        if raw_id in MISSING_TOKENS:
            table.no_id += 1
            continue
        stem, _ = strip_version(raw_id)
        if stem in table.rows:
            raise ProjectError(4, f"{path.name}, data row {row_no}: gene {stem} appears more than once")
        table.ids[stem] = raw_id
        table.rows[stem] = [row[at[r]].strip() for r in table.roles]
        symbol = row[at[symbol_col]].strip() if symbol_col else ""
        # Paralome writes the ID where a gene has no symbol: that is not a symbol, and must not hide another
        # file's real one
        table.symbols[stem] = "" if symbol in MISSING_TOKENS or symbol in (raw_id, stem) else symbol
    return table


def combine(files: Sequence[Path], out: Path, *, drop: Sequence[str] = (), id_column: Optional[str] = None,
            log=print) -> Path:
    """Write the wide table for ``files`` (one comparison each) to ``out``. Returns ``out``."""
    from .identity import DISPLAY_NAME
    from .version import __version__

    out = Path(out)
    if out.is_dir():
        raise ProjectError(31, f"--out {out} is a folder; give the file to write, e.g. {out / 'wide.csv'}")
    files = [Path(f) for f in files]
    names = [comparison_name(f, drop) for f in files]
    for name, path in zip(names, files):
        if not _NAME.match(name):
            raise ProjectError(31, f"{path.name}: comparison name {name!r} is empty or has spaces or punctuation; "
                                   "rename the file or change --drop")
    repeated = sorted({n for n in names if names.count(n) > 1})
    if repeated:
        raise ProjectError(31, f"more than one file gives the comparison name {', '.join(repeated)}; change --drop")

    tables = []
    order: List[str] = []
    symbols: Dict[str, str] = {}
    ids: Dict[str, str] = {}
    for name, path in zip(names, files):
        t = _read(path, name, id_column)
        tables.append(t)
        for gene in t.rows:
            if gene not in ids:
                order.append(gene)
                ids[gene] = t.ids[gene]
            if not symbols.get(gene) and t.symbols.get(gene):
                symbols[gene] = t.symbols[gene]
        note = f"; {t.no_id} row(s) without a gene ID left out" if t.no_id else ""
        tool = f", {t.tool} results" if t.tool else ""
        log(f"  {name:<40} {len(t.rows):>7} genes  ({path.name}{tool}{note})")

    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    folders = sorted({str(f.resolve().parent) for f in files})
    sources = sorted({t.generator for t in tables if t.generator})
    made_by = f"{'; '.join(sources)}; " if sources else ""
    try:
        with tmp.open("w", newline="", encoding="utf-8") as handle:
            handle.write(f"# generated by {made_by}{DISPLAY_NAME} {__version__} combine, from {len(files)} "
                         f"per-comparison table(s) in {', '.join(folders)}\n")
            for t in tables:
                if t.arms and "ANOVA_p" not in t.roles:
                    test = f", test={t.test}" if t.test else ""
                    handle.write(f"#   {t.name}: per-comparison table {t.path.name}, control={t.arms[0]}, "
                                 f"treatment={t.arms[1]}{test}\n")
                elif t.tool:
                    handle.write(f"#   {t.name}: per-comparison table {t.path.name}, {t.tool} results\n")
                else:
                    handle.write(f"#   {t.name}: per-comparison table {t.path.name}\n")
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["Gene", "Common_name"] + [f"{t.name}_{r}" for t in tables for r in t.roles])
            for gene in order:
                row: List[str] = [ids[gene], symbols.get(gene, "")]
                for t in tables:
                    row += t.rows.get(gene) or [""] * len(t.roles)
                writer.writerow(row)
        tmp.replace(out)
    except OSError as exc:
        tmp.unlink(missing_ok=True)
        raise ProjectError(30, f"cannot write {out}: {exc.strerror or exc}") from None
    log(f"combined {len(files)} table(s): {len(order)} genes -> {out}")
    return out


def table_or_folder(path: Path, out_dir: Path, *, pattern: str = DEFAULT_PATTERN, drop: Sequence[str] = (),
                    id_column: Optional[str] = None, log=print) -> Tuple[Path, Optional[Path], Optional[str]]:
    """(the wide table to read, the folder or single-comparison file it was combined from or None, the ID column
    to read it by). A wide table is returned as it is; a folder, or a file holding one comparison
    (``is_single_comparison``), is combined into ``out_dir/input_combined.csv``, whose ID column is always
    'Gene' -- --id-column applies to the per-comparison files."""
    path = Path(path)
    target = Path(out_dir) / COMBINED_NAME
    if not path.is_dir():
        if not is_single_comparison(path):
            return path, None, id_column
        log(f"{path} holds one comparison: reading it as a per-comparison table")
        return combine([path], target, drop=drop, id_column=id_column, log=log), path, None
    log(f"{path} is a folder: combining its per-comparison tables")
    files = find_tables(path, pattern, out=target)
    return combine(files, target, drop=drop, id_column=id_column, log=log), path, None
