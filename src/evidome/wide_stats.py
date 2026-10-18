"""Reader for the wide gene-statistics table -- the program's input contract.

One row per gene, one block of columns per comparison, optional ``#`` comment
lines at the top. Comma- or tab-separated: the column header line decides (a tab in
it means tab-separated). Paralome's ``all_genes_clean_statistics.csv`` is the first
producer of this format, but nothing here depends on Paralome: any
table that follows the column grammar below is a valid input.

Column grammar
--------------
Gene identity:
    an Ensembl gene ID column -- mandatory. Recognised by name (see
    ``ID_COLUMN_NAMES``) or named explicitly by the caller.
    an optional symbol column (see ``SYMBOL_COLUMN_NAMES``). A symbol equal to
    the gene ID means "no symbol" (Paralome falls back to the ID).

Directional comparison ("contrast") -- every column prefixed ``{name}_``:
    P-value, Q-value                         (required)
    exactly one of                           (required):
        Fold_change                          (in the producer's convention --
                                              signed times, a ratio or a log
                                              fold change; read from the values)
        Log2_fold_change                     (declared log2: 0 unchanged, +1
                                              doubled, -1 halved; never
                                              re-guessed from the values)
    one pair of group means                  (optional; a half
                                              pair is refused):
        Raw_FPKM_Control_mean, Raw_FPKM_Treatment_mean   (bulk RNA, FPKM)
        Control_mean, Treatment_mean                     (any other unit --
                                              pseudobulk counts, cell-level
                                              expression, accessibility; the
                                              producer's header says which)
    Effect, T-statistic                      (optional; present for model
                                              contrasts, absent for t-tests)

Burden test -- a contrast whose header line ends ``test=burden``::

    # dominant_ptv: control=gnomAD, treatment=probands, test=burden

is read as kind ``burden``: a gene-burden test of observed carriers (treatment,
e.g. probands) against the number expected from a reference (control, e.g.
gnomAD). Its columns are a contrast's, usually Log2_fold_change (log2 observed
/ expected), P-value and Q-value with no means. It has a direction like a
contrast: up = more carriers than expected. ``test=`` with any other value is
refused, never ignored.

A burden test -- and only a burden test -- may also carry::

    {name}_Null_selection_probability        (optional)

the probability, under the burden test's null, that the gene reaches
p < ``NULL_SELECTION_ALPHA`` (0.05, the ORA threshold's default): a gene with
many common qualifying alleles can reach it, a gene whose alleles are all
extremely rare almost never can. Over-representation analysis of a burden
comparison uses it as the weights of its null (enrich.py: Wallenius'
noncentral hypergeometric, as goseq does). When present, every gene with a
P-value in that comparison must have it, within [0, 1]. On any other kind of
comparison it is refused.

Different statistical tests write different subsets of the optional columns;
nothing downstream depends on which test produced a block. The means are
carried to the outputs as given and never enter a statistic.

Direction-less test ("F-test") -- every column prefixed ``{name}_``:
    F, ANOVA_p, ANOVA_q                      (all required)

Missing values are ``.`` or empty. A gene whose P-value is missing in a
comparison was not tested in it (Paralome writes all-zero genes that way) --
recorded as untested, never silently treated as p = 1.

Every column must be claimed by the identity columns or by exactly one
comparison block. An unclaimed column is reported, not ignored: a column the
reader does not understand is a column whose meaning an agent would have to
guess.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError

MISSING_TOKENS = frozenset({"", ".", "NA", "NaN", "nan"})

ID_COLUMN_NAMES = ("Gene", "gene_id", "Ensembl_ID", "ensembl_gene_id")
SYMBOL_COLUMN_NAMES = ("Common_name", "gene_name", "Symbol", "symbol")

CONTRAST_REQUIRED = ("P-value", "Q-value")
#: Exactly one per contrast: Fold_change in the producer's convention (read from its values), or
#: Log2_fold_change, declared log2 and never re-guessed.
FOLD_CHANGE_ROLES = ("Fold_change", "Log2_fold_change")
#: Group means: optional (a burden test has none), but a pair is complete or absent; which pair says
#: what unit they are in.
MEAN_PAIRS = (("Raw_FPKM_Control_mean", "Raw_FPKM_Treatment_mean"), ("Control_mean", "Treatment_mean"))
#: A burden test's per-gene probability of reaching p < NULL_SELECTION_ALPHA under its null (the weights
#: of enrich's detectability-adjusted ORA). Only a burden comparison may carry it.
NULL_SELECTION_ROLE = "Null_selection_probability"
NULL_SELECTION_ALPHA = 0.05
CONTRAST_OPTIONAL = ("Effect", "T-statistic", NULL_SELECTION_ROLE)
#: Every role a contrast block can hold, in the order outputs list them (bulk RNA's order unchanged).
CONTRAST_ROLES = (MEAN_PAIRS[0] + MEAN_PAIRS[1] + FOLD_CHANGE_ROLES + ("P-value", "Q-value") + CONTRAST_OPTIONAL)
#: Kinds with a direction of their own, read the same way downstream: a two-group contrast, and a burden test
#: (observed against expected carriers). An F-test ("ftest") is the only other kind.
DIRECTIONAL_KINDS = ("contrast", "burden")
#: The values a header's ``test=`` token may take, and the kind each makes.
HEADER_TESTS = {"burden": "burden"}
FTEST_REQUIRED = ("F", "ANOVA_p", "ANOVA_q")
#: Not in any input: the signed estimate that gives an F-test a direction, taken from its own contrasts
#: (``derive_ftest_directions``). Listed after the F-test's own columns.
DIRECTION_ROLE = "Direction_estimate"
FTEST_ROLES = FTEST_REQUIRED + (DIRECTION_ROLE,)

#: Ensembl stable IDs may carry a version suffix (ENSG00000139618.17). The
#: version is stripped for matching and the fact is recorded.
_VERSION_SUFFIX = re.compile(r"^(?P<stem>.+?)\.(?P<version>\d+)$")

# "#   WD_vs_CD_main_effect: control=CD, treatment=WD"
# "#   WD_vs_CD: pairwise, control=CD, treatment=WD"     (a leading description is allowed)
_HEADER_CONTRAST = re.compile(r"^\s*(?P<name>\S+):\s*(?:[^=]*?,\s*)?control=(?P<control>[^,]+),\s*treatment=(?P<treatment>.+?)\s*$")
# "#   dominant_ptv: control=gnomAD, treatment=probands, test=burden"   (the test token, last)
_HEADER_TEST = re.compile(r",\s*test=(?P<test>[^,\s]*)\s*$")
# "#   WD_main_effect: ANOVA model contrast (pooled-error estimate, ...)" -- described, arms not named
_HEADER_DESCRIBED = re.compile(r"^\s*(?P<name>\S+):\s*(?P<text>\S.*?)\s*$")
_HEADER_FTEST = re.compile(r"^\s*(?P<name>\S+):\s*ANOVA F-test")
_HEADER_GENERATOR = re.compile(r"^\s*generated by (?P<generator>.+?)\s*$")


@dataclass
class Comparison:
    name: str
    kind: str  # "contrast", "burden" (a contrast-like burden test) or "ftest"
    columns: Dict[str, str]  # role (e.g. "P-value") -> column name in the file
    control: Optional[str] = None
    treatment: Optional[str] = None
    tested: int = 0
    direction_from: Optional[str] = None   # F-tests: how their direction was taken from their contrasts
    same_test_as: Optional[str] = None     # a main-effect F-test: the main-effect contrast that is the same test
    source_name: Optional[str] = None      # the name in the file, when an F-test shares it with a contrast

    @property
    def p_role(self) -> str:
        return p_role_of(self.kind)

    @property
    def fold_change_role(self) -> Optional[str]:
        """'Fold_change' or 'Log2_fold_change' (whichever the block has); None for an F-test."""
        return next((r for r in FOLD_CHANGE_ROLES if r in self.columns), None)


def is_directional(kind: str) -> bool:
    """A contrast or a burden test: a comparison with a direction of its own."""
    return kind in DIRECTIONAL_KINDS


def p_role_of(kind: str) -> str:
    return "P-value" if is_directional(kind) else "ANOVA_p"


def delimiter_of(header_line: str) -> str:
    """Tab when the column header line holds a tab, else comma."""
    return "\t" if "\t" in header_line else ","


@dataclass
class WideStatistics:
    path: Path
    header_lines: List[str]
    generator: Optional[str]
    id_column: str
    symbol_column: Optional[str]
    gene_ids: List[str]  # version-stripped
    input_gene_ids: List[str]  # exactly as in the file
    symbols: List[Optional[str]]
    comparisons: List[Comparison]
    values: Dict[str, List[Optional[float]]]  # column -> per-gene values
    versions_stripped: int = 0
    notes: List[str] = field(default_factory=list)

    def column(self, comparison: Comparison, role: str) -> List[Optional[float]]:
        return self.values[comparison.columns[role]]


def strip_version(gene_id: str) -> Tuple[str, bool]:
    match = _VERSION_SUFFIX.match(gene_id)
    if match and match.group("stem").upper().startswith("ENS"):
        return match.group("stem"), True
    return gene_id, False


def _to_float(raw: str, column: str, line: int) -> Optional[float]:
    text = raw.strip()
    if text in MISSING_TOKENS:
        return None
    try:
        value = float(text)
    except ValueError as exc:
        raise ProjectError(7, f"line {line}, column {column!r}: {raw!r} is not a number") from exc
    if math.isnan(value):
        return None
    # An infinite T-statistic (R writes Inf for a zero-variance gene) would become an infinite GSEA weight and
    # regulator response and turn every score it touched into nan. Only a fold change may be infinite: a ratio
    # of a gene switched on from zero is a real, signed answer.
    if math.isinf(value) and not column.endswith(("_Fold_change", "_Log2_fold_change")):
        raise ProjectError(7, f"line {line}, column {column!r}: {raw!r} is infinite; only a fold change may be "
                              "infinite -- leave the cell blank or give a finite value")
    return value


def _parse_header(lines: Sequence[str]) -> Tuple[Optional[str], Dict[str, Tuple[str, Optional[str], Optional[str]]]]:
    generator = None
    # Keyed by (name, kind): a contrast and an F-test may share a name (Paralome names a pairwise
    # test <treatment>_vs_<control> and a main-effect F-test whose levels restate the factor the same way).
    declared: Dict[Tuple[str, str], Tuple[Optional[str], Optional[str]]] = {}
    for line in lines:
        body = line.lstrip("#")
        m = _HEADER_GENERATOR.match(body)
        if m and generator is None:
            generator = m.group("generator")
            continue
        # ", test=<kind>" after the arms names the kind of test (only 'burden' so far). It is checked
        # here, not ignored: an unknown test, or one on a line without arms, is refused.
        kind = "contrast"
        t = _HEADER_TEST.search(body)
        if t:
            if t.group("test") not in HEADER_TESTS:
                raise ProjectError(9, f"header line {line.strip()!r}: test={t.group('test')!r} is not a known test "
                                      f"(known: {', '.join(sorted(HEADER_TESTS))})")
            kind = HEADER_TESTS[t.group("test")]
            body = body[: t.start()]
        m = _HEADER_CONTRAST.match(body)
        if m:
            declared[(m.group("name"), kind)] = (m.group("control").strip(), m.group("treatment").strip())
            continue
        if t:
            raise ProjectError(9, f"header line {line.strip()!r}: test={t.group('test')} needs the comparison's "
                                  "control= and treatment= before it")
        m = _HEADER_FTEST.match(body)
        if m:
            declared[(m.group("name"), "ftest")] = (None, None)
            continue
        m = _HEADER_DESCRIBED.match(body)
        if m and not any(name == m.group("name") for name, _ in declared):
            declared[(m.group("name"), "described")] = (None, None)
    return generator, declared


def _find_blocks(columns: Sequence[str]) -> Tuple[List[Comparison], List[str]]:
    """Group columns into comparison blocks by suffix. Returns (blocks, unclaimed)."""
    colset = set(columns)
    blocks: List[Comparison] = []
    claimed = set()
    for col in columns:
        if col.endswith("_P-value"):
            name = col[: -len("_P-value")]
            roles = {r: f"{name}_{r}" for r in CONTRAST_ROLES if f"{name}_{r}" in colset}
            missing = [f"{name}_{r}" for r in CONTRAST_REQUIRED if r not in roles]
            fold = [r for r in FOLD_CHANGE_ROLES if r in roles]
            if not fold:
                missing.append(" or ".join(f"{name}_{r}" for r in FOLD_CHANGE_ROLES))
            elif len(fold) > 1:     # two readings of one number: which is meant would be a guess
                raise ProjectError(6, f"comparison {name!r} has both {' and '.join(f'{name}_{r}' for r in fold)}; "
                                      "give exactly one")
            # The means are optional (a burden test has none), but half a pair is still refused
            for pair in MEAN_PAIRS:
                if (pair[0] in roles) != (pair[1] in roles):
                    missing.append(f"{name}_{pair[1] if pair[0] in roles else pair[0]}")
            if missing:
                raise ProjectError(6, f"comparison {name!r} lacks {', '.join(missing)}")
            blocks.append(Comparison(name, "contrast", roles))
            claimed.update(roles.values())
        elif col.endswith("_ANOVA_p"):
            name = col[: -len("_ANOVA_p")]
            roles = {r: f"{name}_{r}" for r in FTEST_REQUIRED if f"{name}_{r}" in colset}
            missing = [r for r in FTEST_REQUIRED if r not in roles]
            if missing:
                raise ProjectError(6, f"comparison {name!r} lacks {', '.join(f'{name}_{r}' for r in missing)}")
            blocks.append(Comparison(name, "ftest", roles))
            claimed.update(roles.values())
    unclaimed = [c for c in columns if c not in claimed]
    return blocks, unclaimed


def read_wide_statistics(path: Path, *, id_column: Optional[str] = None) -> WideStatistics:
    path = Path(path)
    # utf-8-sig: Excel's "CSV UTF-8" starts the file with a byte-order mark, which otherwise becomes part of the
    # first column's name ('\ufeffGene') and the ID column would be reported missing.
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            raw_lines = handle.read().splitlines()
    except UnicodeDecodeError as exc:
        line = exc.object[: exc.start].count(b"\n") + 1
        raise ProjectError(7, f"{path} is not UTF-8 text (line {line}, byte 0x{exc.object[exc.start]:02x}); "
                              "save it as UTF-8") from None
    except OSError as exc:      # missing, a directory, unreadable
        raise ProjectError(30, f"cannot read {path}: {exc.strerror or exc}") from None
    header_lines = []
    skipped = 0
    for line in raw_lines:
        if line.startswith("#"):
            header_lines.append(line)
        elif not line.strip():      # a blank line before the header is skipped, not read as an empty header row
            skipped += 1
        else:
            break
    body = raw_lines[len(header_lines) + skipped:]
    reader = csv.reader(body, delimiter=delimiter_of(body[0]) if body else ",")
    try:
        columns = next(reader)
    except StopIteration:
        raise ProjectError(5, f"{path}: no column header after the comment lines") from None
    columns = [c.strip() for c in columns]

    if id_column is None:
        id_column = next((c for c in ID_COLUMN_NAMES if c in columns), None)
        if id_column is None:
            raise ProjectError(3, f"{path}: none of {', '.join(ID_COLUMN_NAMES)} is a column; name it explicitly")
    elif id_column not in columns:
        raise ProjectError(3, f"{path}: requested ID column {id_column!r} is not in the file")
    symbol_column = next((c for c in SYMBOL_COLUMN_NAMES if c in columns and c != id_column), None)

    duplicated = sorted({c for c in columns if columns.count(c) > 1})
    if duplicated:     # otherwise the last copy would silently win and the comparison be read twice
        raise ProjectError(8, f"{path}: column name(s) appear more than once: {', '.join(duplicated)}")
    data_columns = [c for c in columns if c not in (id_column, symbol_column)]
    comparisons, unclaimed = _find_blocks(data_columns)
    if not comparisons:
        raise ProjectError(5, f"{path}: no comparison blocks recognised (no *_P-value or *_ANOVA_p columns)")
    if unclaimed:
        raise ProjectError(8, f"{path}: columns not part of any comparison block: {', '.join(unclaimed)}")

    generator, declared = _parse_header(header_lines)
    # A contrast-shaped block its header calls a burden test is one; the kind is the header's word.
    for comp in comparisons:
        if comp.kind == "contrast" and (comp.name, "burden") in declared:
            comp.kind = "burden"
    notes: List[str] = []
    for comp in comparisons:
        if (comp.name, comp.kind) in declared:
            comp.control, comp.treatment = declared[(comp.name, comp.kind)]
        elif (comp.name, "described") in declared:
            continue
        elif any(name == comp.name for name, _ in declared):
            other = next(kind for name, kind in declared if name == comp.name)
            if any(c.name == comp.name and c.kind == other for c in comparisons):
                # the header describes the other test of a shared name
                notes.append(f"comparison {comp.name!r} ({comp.kind}) is not described in the comment header")
                continue
            raise ProjectError(9, f"comparison {comp.name!r}: header says {other}, columns say {comp.kind}")
        elif header_lines:
            notes.append(f"comparison {comp.name!r} is not described in the comment header")
    for comp in comparisons:
        if NULL_SELECTION_ROLE in comp.columns and comp.kind != "burden":
            raise ProjectError(8, f"{path}: column {comp.columns[NULL_SELECTION_ROLE]!r} belongs to a burden test "
                                  f"only, and comparison {comp.name!r} is a {comp.kind} (its header line needs "
                                  "', test=burden')")
    present = {(c.name, c.kind) for c in comparisons}
    for name, kind in declared:
        if kind != "described" and (name, kind) not in present:
            raise ProjectError(9, f"header describes comparison {name!r} but the file has no columns for it")
    # A name shared by a contrast and an F-test: the F-test is renamed <name>_ANOVA everywhere downstream (file
    # names, ballots, the evidence), keeping its name in the file for matching it to its contrasts.
    contrast_names = {c.name for c in comparisons if is_directional(c.kind)}
    for comp in comparisons:
        if comp.kind == "ftest" and comp.name in contrast_names:
            comp.source_name, comp.name = comp.name, f"{comp.name}_ANOVA"
            notes.append(f"F-test {comp.source_name!r} shares its name with a contrast: reported as {comp.name!r}")

    index = {c: i for i, c in enumerate(columns)}
    input_ids: List[str] = []
    gene_ids: List[str] = []
    symbols: List[Optional[str]] = []
    values: Dict[str, List[Optional[float]]] = {c: [] for c in data_columns}
    stripped = 0
    first_data_line = len(header_lines) + skipped + 2
    no_id = 0
    for offset, row in enumerate(reader):
        line = first_data_line + offset
        if not row or all(not cell.strip() for cell in row):
            continue
        if len(row) != len(columns):
            raise ProjectError(7, f"line {line}: {len(row)} fields, header has {len(columns)}")
        raw_id = row[index[id_column]].strip()
        if raw_id in MISSING_TOKENS:      # '.', 'NA', empty: no ID to map by
            no_id += 1
            continue
        stem, was_versioned = strip_version(raw_id)
        stripped += was_versioned
        input_ids.append(raw_id)
        gene_ids.append(stem)
        symbol = row[index[symbol_column]].strip() if symbol_column else ""
        symbols.append(symbol if symbol and symbol not in (raw_id, stem) and symbol not in MISSING_TOKENS else None)
        for col in data_columns:
            values[col].append(_to_float(row[index[col]], col, line))

    if no_id:
        notes.append(f"{no_id} row(s) with no gene ID ('.', 'NA' or empty) left out")
    if not gene_ids:
        raise ProjectError(3, f"{path}: no row has a gene ID")
    seen: Dict[str, int] = {}
    for gid in gene_ids:
        seen[gid] = seen.get(gid, 0) + 1
    dups = [g for g, n in seen.items() if n > 1]
    if dups:
        raise ProjectError(4, f"{len(dups)} gene ID(s) appear more than once, e.g. {dups[0]}")

    table = WideStatistics(path, header_lines, generator, id_column, symbol_column, gene_ids, input_ids,
                           symbols, comparisons, values, stripped, notes)
    for comp in comparisons:
        comp.tested = sum(v is not None for v in table.column(comp, comp.p_role))
        if NULL_SELECTION_ROLE in comp.columns:
            _check_selection(table, comp)
    derive_ftest_directions(table)
    return table


def _check_selection(table: WideStatistics, comp: Comparison) -> None:
    """A burden test's null selection probabilities: one per tested gene, each within [0, 1]."""
    column = comp.columns[NULL_SELECTION_ROLE]
    for gene, p, w in zip(table.input_gene_ids, table.column(comp, comp.p_role), table.values[column]):
        if p is None:
            continue
        if w is None:
            raise ProjectError(7, f"{table.path}: gene {gene} has a {comp.name}_P-value but no {column}; the "
                                  "detectability-weighted ORA needs one for every tested gene")
        if not 0.0 <= w <= 1.0:
            raise ProjectError(7, f"{table.path}: gene {gene}, column {column!r}: {w!r} is not a probability")


def derive_ftest_directions(table: WideStatistics) -> None:
    """Give each F-test the direction its own contrasts show.

    Every factor of a factorial ANOVA has two levels, so each F-test has one degree of freedom and one direction:

    main effect  ``<factor>_<alt>_vs_<ref>`` (or ``<alt>_vs_<ref>``): the direction of its main-effect contrast,
                 ``<alt>_vs_<ref>_main_effect`` -- the same comparison, estimated from the same model. Its Effect
                 (log2) where given, else its T-statistic.
    interaction  ``<A>_x_<B>``: the difference of differences, Effect(<A contrast> within <B alt>) minus
                 Effect(<A contrast> within <B ref>) -- the interaction coefficient itself; positive when A's effect
                 is larger (more up) at B's alternative level. Needs both simple-effect contrasts and their Effect.

    The value is added as a derived column, role ``Direction_estimate``, and ``direction_from`` says what it came
    from. An F-test whose contrasts are not in the table keeps no direction, as before. Names follow the grammar
    the input's producer writes (Paralome); nothing is inferred from a name that does not match exactly.
    """
    contrasts = {c.name: c for c in table.comparisons if c.kind == "contrast"}

    def signed(comp: Comparison) -> Tuple[Optional[str], Optional[List[Optional[float]]]]:
        for role in ("Effect", "T-statistic"):
            if role in comp.columns:
                return role, table.column(comp, role)
        return None, None

    main_stems = {name[: -len("_main_effect")]: c for name, c in contrasts.items() if name.endswith("_main_effect")}
    # A main-effect F-test is named <factor>_<stem> (e.g. diet_WD_vs_CD): that is the only place the factor's own
    # name meets its levels, and an interaction F-test is named by factors (diet_x_metformin).
    factor_stems: Dict[str, str] = {}
    for comp in table.comparisons:
        fname = comp.source_name or comp.name
        if comp.kind == "ftest" and "_x_" not in fname:
            for st in sorted(main_stems, key=len, reverse=True):
                if fname.endswith("_" + st):
                    factor_stems[fname[: -len(st) - 1]] = st
                    break
    for comp in table.comparisons:
        if comp.kind != "ftest":
            continue
        fname = comp.source_name or comp.name
        stem = next((st for st in sorted(main_stems, key=len, reverse=True)
                     if fname == st or fname.endswith("_" + st)), None)
        values: Optional[List[Optional[float]]] = None
        if stem is not None:
            comp.same_test_as = main_stems[stem].name
            role, values = signed(main_stems[stem])
            if values is not None:
                comp.direction_from = f"{main_stems[stem].name} ({role})"
        elif "_x_" in fname:
            values, how = _interaction_direction(table, contrasts, main_stems, fname, factor_stems)
            if values is not None:
                comp.direction_from = how
        if values is None:
            continue
        column = f"{comp.name}__{DIRECTION_ROLE}"
        table.values[column] = list(values)
        comp.columns[DIRECTION_ROLE] = column
        table.notes.append(f"F-test {comp.name!r}: direction taken from {comp.direction_from}")


def drop_duplicate_ftests(table: WideStatistics) -> List[Tuple[str, str]]:
    """Remove each main-effect F-test whose main-effect contrast is in the table.

    With two levels per factor they are one test (F = t^2 in a balanced design; the database-utility table showed
    identical results), so analysing both reports every finding twice. The contrast is kept -- it carries the
    direction, fold change and effect size -- and the F-test is listed as skipped, with the contrast standing in
    for it. Interaction F-tests have no contrast of their own and are kept. Returns [(F-test, contrast)]."""
    dropped = [(c.name, c.same_test_as) for c in table.comparisons if c.kind == "ftest" and c.same_test_as]
    table.comparisons = [c for c in table.comparisons if not (c.kind == "ftest" and c.same_test_as)]
    for name, contrast in dropped:
        table.notes.append(f"F-test {name!r} not analysed: it is the same test as {contrast!r}, analysed in its place")
    return dropped


def _interaction_direction(table: WideStatistics, contrasts: Dict[str, Comparison],
                           main_stems: Dict[str, Comparison], name: str = "",
                           factor_stems: Optional[Dict[str, str]] = None) -> Tuple[Optional[List[Optional[float]]], str]:
    """Direction of the interaction F-test ``name`` (<A>_x_<B>): Effect(within B alt) - Effect(within B ref).

    The pair of factors is the one the F-test names (names are matched exactly), so in a three-factor design one
    interaction never takes its sign from another. A pair is accepted only if it is the one the
    name identifies -- through the main-effect F-tests' <factor>_<stem> names -- or, when those are absent, if it is
    the only candidate pair in the table. Three-way and unmatched names get no direction."""
    parts = name.split("_x_") if name else []
    wanted: Optional[set] = None
    if len(parts) == 2 and factor_stems and all(p in factor_stems for p in parts):
        wanted = {factor_stems[parts[0]], factor_stems[parts[1]]}
    elif len(parts) != 2:
        return None, ""
    candidates = []
    for stem_a in sorted(main_stems):
        within = {}
        for cname, c in contrasts.items():
            if cname.startswith(stem_a + "_within_") and cname.endswith("_contrast") and "Effect" in c.columns:
                within[cname[len(stem_a) + len("_within_"): -len("_contrast")]] = c
        if len(within) != 2:
            continue
        for stem_b in main_stems:
            if stem_b == stem_a:
                continue
            for alt, ref in ((a, b) for a in within for b in within if a != b):
                if stem_b == f"{alt}_vs_{ref}":
                    candidates.append((stem_a, stem_b, within[alt], within[ref]))
    if wanted is not None:
        candidates = [c for c in candidates if {c[0], c[1]} == wanted]
    elif len({frozenset((c[0], c[1])) for c in candidates}) != 1:
        return None, ""
    if not candidates:
        return None, ""
    _, _, hi_c, lo_c = candidates[0]
    hi, lo = table.column(hi_c, "Effect"), table.column(lo_c, "Effect")
    diff = [None if x is None or y is None else x - y for x, y in zip(hi, lo)]
    return diff, f"{hi_c.name} (Effect) minus {lo_c.name} (Effect)"
