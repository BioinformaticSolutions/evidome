# Evidome

Evidome turns a gene-level differential expression result into a structured,
self-describing **functional evidence report**: gene-set enrichment (Reactome,
Gene Ontology, WikiPathways, MSigDB Hallmark; ORA and preranked GSEA) and
upstream-regulator evidence (CollecTRI and OmniPath), written for an AI agent
to read and interpret. The program itself makes no interpretive claims.

It works across species. It identifies the species from the input's Ensembl
gene IDs, maps every gene to its human ortholog (Ensembl Compara), records how
each mapping was made, and resolves many-to-one orthologs by majority vote
across the statistics.

Evidome is standalone software. Paralome runs it as an external tool, but
nothing about Evidome depends on Paralome.

## Install

```
python3 -m pip install .            # from a copy of this repository
```

Requires Python 3.9 or later. pip installs the two required packages,
matplotlib and scikit-learn (see `NOTICE`).

## Reference data

Fetched once per machine, which needs internet access for this step:

```
evidome reference genesets
evidome reference regulators
evidome reference pathway-signs
evidome reference fetch --species <name>      # per species
```

The first three are shared by every species. A run on a species that has not
been fetched fetches it by itself.

## Use

```
evidome TABLE --out DIR                   # every step, on this machine
evidome TABLE --out DIR --backend slurm   # every step, as SLURM jobs
evidome --help                            # the individual steps
```

## Input

TABLE is one of:

- a wide statistics table: one row per gene, an Ensembl gene ID column, and a
  block of columns per comparison (`<name>_P-value`, `<name>_Q-value`, and
  `<name>_Log2_fold_change` or `<name>_Fold_change`);
- one comparison's table, or a folder of them (one comparison per file, named
  after the file), which is combined into a wide table first;
- results tables from DESeq2, edgeR, limma or Seurat's FindMarkers, written
  from R as they come (`write.csv(...)`), one comparison per file, read by
  their own column names.

Gene IDs may be Ensembl gene IDs or NCBI GeneIDs, whose species is found from
the IDs, or gene symbols, which need the species named: `--species
macaca_mulatta`. Symbols and NCBI GeneIDs are translated to the species'
Ensembl genes; an ID that names no gene, or several, is kept and reported,
never guessed. A folder is read for its `.csv`, `.tsv` and `.txt` files.
`evidome inspect TABLE` shows what Evidome reads from a table without running
anything. The full format -- every column, the comment header, F-tests,
gene-burden tests, what each column is used for, what is refused and why --
is in [`src/evidome/docs/INPUT_FORMAT.md`](https://github.com/BioinformaticSolutions/evidome/blob/main/src/evidome/docs/INPUT_FORMAT.md),
which is also installed with the package.

## Output

The evidence zip (`<out>_combined_results.zip`) holds every table as CSV, a
manifest that describes each file and column, and
`Evidome_Viewer_v<version>.html`, a page for browsing the results. Beside the
zip, `<out>_workbook.xlsx` shows the same tables as one workbook. Open the
page in a web browser (unzipped beside the files, or on its own) and choose
the zip or the unzipped folder. It lists the findings, ORA, GSEA, regulators
and genes of each comparison, with a pop-up per gene set, regulator and gene
that exports as TSV. The page reads the files on your computer and uploads
nothing.

## Develop

```
python3 -m pytest            # from the repository root
PYTHONPATH=src python3 -m evidome.cli --version
```

## Licence

Copyright (C) 2026 Kenneth L. Jones. Dual licensed: the GNU Affero General
Public License, version 3 (see `LICENSE`), or a separate commercial licence.
See `NOTICE`, which also lists the reference data sources and their licences.
