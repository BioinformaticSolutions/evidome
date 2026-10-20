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
python3 -m pip install evidome
```

Requires Python 3.9 or later. pip installs the two required packages,
matplotlib and scikit-learn (see `NOTICE`).

## Reference data

Name the directory the reference data and indexes live in, or should be built
in, with `--reference-dir DIR`. It is remembered for your later runs, so it is
given once. (With Paralome set up on the machine, they live beside Paralome's
genomes and nothing needs giving.)

Fetched once, which needs internet access for this step:

```
evidome reference fetch --species <name> --reference-dir DIR   # first; then once per species
evidome reference genesets
evidome reference regulators
evidome reference pathway-signs
```

`fetch` comes first: the gene sets and regulators are indexed to the human
genes it downloads, and the pathway signs are built from the gene sets. Gene
sets, regulators and pathway signs are shared by every species. A run on a
species that has not been fetched fetches it by itself.

On a shared machine the reference directory belongs to a group: it is created
group-writable, and its group is passed on to everything made in it. A run can
write there -- a species not yet fetched, or the first index of a registered
genome -- so every account that runs Evidome must be in that group.

To identify your data by the exact annotation they were counted against
(rather than by Ensembl's own), give that annotation with the run:

```
evidome TABLE --out DIR --gtf <annotation.gtf>
```

Its orthologs then come from the Ensembl release the annotation belongs to,
read from the name Ensembl gave the GTF (`Species.Assembly.<release>.gtf`), or
from the dataset of Ensembl's new platform (`genes.gtf.gz`). Or register it
once, for every run:

```
evidome reference genome add <name> --gtf <annotation.gtf>
evidome reference genome list
evidome reference genome remove <name>    # only the account that added it
```

On a machine with Paralome, the genomes registered in Paralome are used as
they are.

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

The evidence zip (`<out>_evidome_output.zip`) holds every table as CSV and a
manifest that describes each file and column. No Excel workbook is
written. To browse the results, open the zip (or the unzipped folder) in the Evidence Browser of
[OmicView](https://github.com/BioinformaticSolutions/omicview), hosted or
downloaded. It lists the findings, ORA, GSEA, regulators and genes of each
comparison, with a pop-up per gene set, regulator and gene that exports as TSV.
It reads the files on your computer and uploads nothing.

## Develop

```
python3 -m pytest            # from the repository root
PYTHONPATH=src python3 -m evidome.cli --version
```

## Licence

Copyright (C) 2026 Kenneth L. Jones. Dual licensed: the GNU Affero General
Public License, version 3 (see `LICENSE`), or a separate commercial licence.
See `NOTICE`, which also lists the reference data sources and their licences.
