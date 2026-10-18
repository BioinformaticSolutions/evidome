"""Gene symbols (with --species) and NCBI GeneIDs are translated to the species' Ensembl genes before the
crosswalk; an ID that names no gene, or several, stays in the outputs untranslated, with the reason."""

import csv

import pytest

from evidome import gene_ids as G
from evidome.cli import main
from evidome.reference import release_dir

from test_evidence_bundle_describes_itself import COLLECTRI, SIGNALLING, TABLE  # noqa: F401
from test_gene_sets_are_indexed_to_human_ensembl import sources  # noqa: F401  (fixture)
from test_one_command_runs_everything import OPTS, world  # noqa: F401  (fixture)

M = "macaca_mulatta"
# The macaque annotation's symbols: Dupa in another case, TWIN on two genes
SYMBOLS = {"ENSMMUG00000000001": "TSPAN6", "ENSMMUG00000000002": "PARA1", "ENSMMUG00000000003": "PARA1L",
           "ENSMMUG00000000004": "Dupa", "ENSMMUG00000000006": "TWIN", "ENSMMUG00000000007": "TWIN",
           "ENSMMUG00000000008": "H4C1"}


def _symbols_in_annotation(ref):
    path = release_dir(ref, 116) / "species" / M / "species_genes.tsv"
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, ["gene_id", "symbol", "biotype"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows({**r, "symbol": SYMBOLS.get(r["gene_id"], r["symbol"])} for r in rows)


def _rekeyed(tmp_path, name, key):
    """The shared fixture's table with its Ensembl IDs replaced by key(id) and no symbol column."""
    lines = TABLE.strip().splitlines()
    out = []
    for line in lines:
        if line.startswith("#"):
            out.append(line)
            continue
        f = line.split(",")
        out.append(",".join(["Gene"] + f[2:]) if f[0] == "Gene" else ",".join([key(f[0])] + f[2:]))
    path = tmp_path / name
    path.write_text("\n".join(out) + "\n")
    return path


def _crosswalk(path):
    with path.open() as handle:
        return {r["input_id_given"]: r for r in csv.DictReader(handle)}


def test_the_kind_of_id_is_read_from_the_column():
    assert G.id_kind(["ENSG00000139618", "ENSMUSG00000000001"]) == G.KIND_ENSEMBL
    assert G.id_kind(["FBgn0000008", "WBGene00000001"]) == G.KIND_ENSEMBL
    assert G.id_kind(["7157", "672", "675"]) == G.KIND_NCBI
    assert G.id_kind(["TP53", "BRCA1", "Gapdh"]) == G.KIND_SYMBOL
    assert G.normalise_species("Macaca mulatta") == M


def test_symbols_without_a_species_are_refused_with_the_fix(world, tmp_path, capsys):
    ref, _ = world
    table = _rekeyed(tmp_path, "sym.csv", lambda g: "SYM" + g[-2:])
    assert main([str(table), "--out", str(tmp_path / "o"), "--reference-dir", str(ref)] + OPTS) != 0
    err = capsys.readouterr().err
    assert "EVD-0032" in err and "--species" in err


def test_symbols_with_their_species_give_the_same_results_as_the_ensembl_ids(world, tmp_path):
    ref, ensembl_table = world
    _symbols_in_annotation(ref)
    unique = dict(SYMBOLS, ENSMMUG00000000006="G6", ENSMMUG00000000007="G7")
    path = release_dir(ref, 116) / "species" / M / "species_genes.tsv"
    text = path.read_text().replace("\tTWIN\t", "\tG6\t", 1).replace("\tTWIN\t", "\tG7\t", 1)
    path.write_text(text)
    by_symbol = _rekeyed(tmp_path, "sym.csv", lambda g: unique[g].upper())      # DUPA: Dupa ignoring case
    by_id = _rekeyed(tmp_path, "ids.csv", lambda g: g)
    a, b = tmp_path / "a", tmp_path / "b"
    assert main([str(by_id), "--out", str(a), "--reference-dir", str(ref)] + OPTS) == 0
    assert main([str(by_symbol), "--out", str(b), "--reference-dir", str(ref), "--species", "Macaca mulatta"]
                + OPTS) == 0
    for name in ("enrichment/A/ora.csv", "enrichment/A/gsea.csv", "regulators/A/regulators.csv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    rows = _crosswalk(b / "ortholog_crosswalk.csv")
    assert rows["TSPAN6"]["input_gene_id"] == "ENSMMUG00000000001"
    assert rows["TSPAN6"]["id_translation"] == G.HOW_SYMBOL_EXACT
    assert rows["DUPA"]["input_gene_id"] == "ENSMMUG00000000004"
    assert rows["DUPA"]["id_translation"] == G.HOW_SYMBOL_CASE
    assert rows["TSPAN6"]["input_symbol"] == "TSPAN6"


def test_a_symbol_naming_no_gene_or_several_is_kept_untranslated_with_the_reason(world, tmp_path):
    ref, _ = world
    _symbols_in_annotation(ref)
    table = _rekeyed(tmp_path, "sym.csv", lambda g: {"ENSMMUG00000000007": "NOSUCH"}.get(g, SYMBOLS[g]))
    out = tmp_path / "o"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref), "--species", M] + OPTS) == 0
    rows = _crosswalk(out / "ortholog_crosswalk.csv")
    assert rows["TWIN"]["input_gene_id"] == "TWIN"
    assert rows["TWIN"]["mapping_status"] == f"not translated: symbol names several {M} Ensembl genes"
    assert rows["TWIN"]["id_translation"].startswith("names 2: ENSMMUG00000000006;ENSMMUG00000000007")
    assert rows["NOSUCH"]["mapping_status"] == f"not translated: symbol not in the {M} annotation"
    assert rows["PARA1"]["mapping_status"] == "mapped"
    human = list(csv.DictReader((out / "comparisons" / "A.human_genes.csv").open()))
    assert human                                         # the run goes on with the genes that were translated


def test_ncbi_gene_ids_find_their_species_and_are_translated(world, tmp_path):
    ref, _ = world
    ncbi = {"ENSMMUG00000000001": "700001", "ENSMMUG00000000006": "700006"}
    table = _rekeyed(tmp_path, "ncbi.csv", lambda g: ncbi.get(g, "8" + g[-2:]))   # the rest: unknown GeneIDs
    out = tmp_path / "o"
    assert main([str(table), "--out", str(out), "--reference-dir", str(ref)] + OPTS) == 0
    rows = _crosswalk(out / "ortholog_crosswalk.csv")
    assert rows["700001"]["input_gene_id"] == "ENSMMUG00000000001" and rows["700001"]["species"] == M
    assert rows["700001"]["id_translation"] == G.HOW_NCBI and rows["700001"]["mapping_status"] == "mapped"
    assert rows["802"]["mapping_status"] == f"not translated: NCBI GeneID with no {M} Ensembl gene"


def test_human_ncbi_gene_ids_are_translated_too(world, tmp_path):
    ref, _ = world
    table = tmp_path / "human.csv"
    table.write_text("Gene,A_Log2_fold_change,A_P-value,A_Q-value\n2,1,0.01,0.05\n4,-1,0.02,0.1\n10,2,0.001,0.01\n")
    out = tmp_path / "o"
    assert main(["crosswalk", str(table), "--out", str(out), "--reference-dir", str(ref), "--threads", "1",
                 "--min-species-match", "0.5"]) == 0
    rows = _crosswalk(out / "ortholog_crosswalk.csv")
    assert rows["2"]["input_gene_id"] == "ENSG00000000002" and rows["2"]["species"] == "homo_sapiens"
    assert rows["10"]["mapping_status"] == "human input (identity)"


def test_a_species_that_contradicts_the_ensembl_ids_is_refused(world, tmp_path, capsys):
    ref, table = world
    assert main([str(table), "--out", str(tmp_path / "o"), "--reference-dir", str(ref), "--species",
                 "homo_sapiens"] + OPTS) != 0
    assert "EVD-0014" in capsys.readouterr().err


def test_inspect_says_what_kind_of_ids_a_table_has(tmp_path, capsys):
    path = tmp_path / "t.csv"
    path.write_text("Gene,A_Log2_fold_change,A_P-value,A_Q-value\nTP53,1,0.01,0.05\nBRCA1,-1,0.02,0.1\n")
    assert main(["inspect", str(path)]) == 0
    out = capsys.readouterr().out
    assert "gene symbols" in out and "needs --species" in out


def test_a_folder_reads_csv_tsv_and_txt_by_default(tmp_path):
    from evidome import combine as C

    folder = tmp_path / "f"
    folder.mkdir()
    for name, sep in (("A.csv", ","), ("B.tsv", "\t"), ("C.txt", "\t")):
        (folder / name).write_text(sep.join(["Gene", "Log2_fold_change", "P-value", "Q-value"]) + "\n"
                                   + sep.join(["ENSG00000000003", "1", "0.01", "0.05"]) + "\n")
    (folder / "D.json").write_text("{}")
    assert [p.name for p in C.find_tables(folder)] == ["A.csv", "B.tsv", "C.txt"]
    assert [p.name for p in C.find_tables(folder, "*.tsv")] == ["B.tsv"]
