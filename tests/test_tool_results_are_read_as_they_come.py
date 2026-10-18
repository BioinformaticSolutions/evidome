"""Results tables of DESeq2, edgeR, limma and Seurat's FindMarkers are read as they come, one comparison each."""

import csv

import pytest

from evidome import combine as combine_mod
from evidome.cli import main
from evidome.errors import ProjectError
from evidome.wide_stats import read_wide_statistics

from test_evidence_bundle_describes_itself import COLLECTRI, SIGNALLING  # noqa: F401
from test_gene_sets_are_indexed_to_human_ensembl import sources  # noqa: F401  (fixture)
from test_one_command_runs_everything import OPTS, world  # noqa: F401  (fixture)

# (gene, log2 fold change, p, adjusted p), the genes of the shared fixture
GENES = [("ENSMMUG00000000001", "2", "0.001", "0.01"), ("ENSMMUG00000000002", "1.585", "0.01", "0.05"),
         ("ENSMMUG00000000003", "0.263", "0.3", "0.5"), ("ENSMMUG00000000004", "-1", "0.02", "0.1"),
         ("ENSMMUG00000000006", "1", "0.04", "0.2"), ("ENSMMUG00000000007", "0", "0.9", "0.9"),
         ("ENSMMUG00000000008", "-1.585", "0.005", "0.03")]

# Each tool's header as R's write.csv writes it (row names first, unnamed), and how a row is filled
TOOLS = {
    "DESeq2": (["", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"],
               lambda g, fc, p, q: [g, "100", fc, "0.2", "3.1", p, q]),
    "edgeR": (["", "logFC", "logCPM", "F", "PValue", "FDR"],
              lambda g, fc, p, q: [g, fc, "5", "9.5", p, q]),
    "limma": (["", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"],
              lambda g, fc, p, q: [g, fc, "7", "4.2", p, q, "1.1"]),
    "Seurat FindMarkers": (["", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"],
                           lambda g, fc, p, q: [g, p, fc, "0.5", "0.2", q]),
}


def _write(path, header, rows, delimiter=","):
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter=delimiter, quoting=csv.QUOTE_NONNUMERIC if delimiter == "," else csv.QUOTE_MINIMAL,
                   lineterminator="\n").writerow(header)
        csv.writer(handle, delimiter=delimiter, lineterminator="\n").writerows(rows)


def _tool_file(path, tool):
    header, fill = TOOLS[tool]
    _write(path, header, [fill(*g) for g in GENES])
    return path


def _rows(path):
    return list(csv.DictReader(line for line in path.read_text().splitlines() if not line.startswith("#")))


@pytest.mark.parametrize("tool", sorted(TOOLS))
def test_each_tools_table_is_read_by_its_own_column_names(tmp_path, tool):
    folder = tmp_path / "results"
    folder.mkdir()
    _tool_file(folder / "WD_vs_CD.csv", tool)
    out = combine_mod.combine(combine_mod.find_tables(folder), tmp_path / "wide.csv", log=lambda *_: None)
    table = read_wide_statistics(out)
    (comparison,) = table.comparisons
    assert comparison.name == "WD_vs_CD" and comparison.kind == "contrast"
    assert table.gene_ids == [g for g, *_ in GENES]
    rows = _rows(out)
    assert [(r["WD_vs_CD_Log2_fold_change"], r["WD_vs_CD_P-value"], r["WD_vs_CD_Q-value"]) for r in rows] == \
        [tuple(g[1:]) for g in GENES]                                                 # copied as text
    carried = {c[len("WD_vs_CD_"):] for c in rows[0] if c.startswith("WD_vs_CD_")}
    assert carried == ({"Log2_fold_change", "P-value", "Q-value", "T-statistic"} if tool == "limma"
                       else {"Log2_fold_change", "P-value", "Q-value"})               # the tool's other columns left out
    assert f"WD_vs_CD.csv, {tool} results" in out.read_text().splitlines()[1]


def test_row_names_from_write_table_and_tabs_are_read(tmp_path):
    """R's write.table(sep="\\t") leaves the row-name column out of the header line."""
    path = tmp_path / "deseq2.tsv"
    header, fill = TOOLS["DESeq2"]
    _write(path, header[1:], [fill(*g) for g in GENES], delimiter="\t")
    out = combine_mod.combine([path], tmp_path / "wide.csv", log=lambda *_: None)
    assert read_wide_statistics(out).gene_ids == [g for g, *_ in GENES]


def test_a_named_id_column_is_used_before_the_row_names(tmp_path):
    path = tmp_path / "A.csv"
    _write(path, ["gene_id", "log2FoldChange", "pvalue", "padj"], [[g, fc, p, q] for g, fc, p, q in GENES])
    out = combine_mod.combine([path], tmp_path / "wide.csv", log=lambda *_: None)
    assert read_wide_statistics(out).gene_ids == [g for g, *_ in GENES]


def test_our_own_columns_win_over_a_tools(tmp_path):
    """A per-comparison table with P-value is read by this program's names even if it also has logFC."""
    path = tmp_path / "A.csv"
    _write(path, ["Gene", "Log2_fold_change", "P-value", "Q-value", "logFC", "PValue", "FDR"],
           [[g, fc, p, q, "9", "9", "9"] for g, fc, p, q in GENES])
    rows = _rows(combine_mod.combine([path], tmp_path / "wide.csv", log=lambda *_: None))
    assert rows[0]["A_P-value"] == "0.001" and "A_PValue" not in rows[0]


def test_an_unknown_table_names_the_formats_it_would_have_read(tmp_path):
    path = tmp_path / "A.csv"
    _write(path, ["", "avg_logFC", "p_val", "p_val_adj"], [[g, fc, p, q] for g, fc, p, q in GENES])
    with pytest.raises(ProjectError) as info:
        combine_mod.combine([path], tmp_path / "wide.csv", log=lambda *_: None)
    assert info.value.code == 31 and "DESeq2: log2FoldChange" in str(info.value) and "avg_log2FC" in str(info.value)


def test_one_tool_file_given_as_the_table_is_one_comparison(tmp_path, capsys):
    path = _tool_file(tmp_path / "WD_vs_CD.csv", "edgeR")
    assert combine_mod.is_single_comparison(path)
    assert main(["inspect", str(path)]) == 0
    assert "WD_vs_CD" in capsys.readouterr().out


def test_a_wide_table_is_not_taken_for_one_comparison(world):
    _, table = world
    assert not combine_mod.is_single_comparison(table)


def test_the_whole_run_reads_a_deseq2_file_as_it_reads_the_same_numbers_in_a_wide_table(world, tmp_path):
    ref, _ = world
    deseq2 = _tool_file(tmp_path / "A.csv", "DESeq2")
    wide = tmp_path / "wide.csv"
    _write(wide, ["Gene", "A_Log2_fold_change", "A_P-value", "A_Q-value"], [list(g) for g in GENES])
    from_wide, from_tool = tmp_path / "from_wide", tmp_path / "from_tool"
    assert main([str(wide), "--out", str(from_wide), "--reference-dir", str(ref)] + OPTS) == 0
    assert main([str(deseq2), "--out", str(from_tool), "--reference-dir", str(ref)] + OPTS) == 0
    assert (from_tool / combine_mod.COMBINED_NAME).exists()
    for name in ("enrichment/A/ora.csv", "enrichment/A/gsea.csv", "regulators/A/regulators.csv"):
        assert (from_wide / name).read_bytes() == (from_tool / name).read_bytes(), name


# --- the input-format document ------------------------------------------------------------------------------------

import re  # noqa: E402
from pathlib import Path  # noqa: E402

import evidome  # noqa: E402

DOC = Path(evidome.__file__).resolve().parent / "docs" / "INPUT_FORMAT.md"
_EXAMPLE = re.compile(r"<!-- example: (?P<name>\S+) shape=(?P<shape>wide|single) -->\n```\n(?P<body>.*?)```", re.S)


def test_the_input_format_document_ships_with_the_package():
    assert DOC.is_file() and "# Evidome input format" in DOC.read_text()


def test_every_example_in_the_document_is_read_as_it_says(tmp_path):
    examples = list(_EXAMPLE.finditer(DOC.read_text()))
    assert len(examples) >= 4
    for number, m in enumerate(examples):
        folder = tmp_path / str(number)
        folder.mkdir()
        path = folder / f"{m['name']}.csv"
        path.write_text(m["body"])
        assert combine_mod.is_single_comparison(path) == (m["shape"] == "single"), m["name"]
        table_path, _, id_column = combine_mod.table_or_folder(path, folder / "out", log=lambda *_: None)
        table = read_wide_statistics(table_path, id_column=id_column)
        assert table.comparisons and len(table.gene_ids) >= 3, m["name"]


def test_the_document_lists_every_tool_format_the_code_reads():
    text = DOC.read_text()
    for tool, columns in combine_mod.TOOL_FORMATS:
        row = next((line for line in text.splitlines() if line.startswith(f"| {tool.split()[0]} ")), "")
        assert row, tool
        for column in columns:
            assert f"`{column}`" in row, (tool, column)
