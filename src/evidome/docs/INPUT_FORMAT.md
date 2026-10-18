# Evidome input format

Evidome reads gene-level differential results: one row per gene, with a
p-value and a direction of change for each comparison. This document describes
every input it accepts, what each column is used for, and what is refused and
why.

Before a run, check a table with:

```
evidome inspect TABLE
```

It reads the table exactly as a run would and prints what it understood: the
program that wrote it, the number of genes and the ID column, how many genes
have a symbol, a first guess at the species, and each comparison with its kind,
its arms, how many genes it tested and how many have p < 0.05, followed by any
notes. It writes nothing; a folder is combined in a temporary folder that is
removed afterwards.

## Contents

1. [The three shapes of input](#1-the-three-shapes-of-input)
2. [Text, delimiters and missing values](#2-text-delimiters-and-missing-values)
3. [Gene identifiers](#3-gene-identifiers)
4. [The wide table](#4-the-wide-table)
5. [The comment header](#5-the-comment-header)
6. [Per-comparison tables](#6-per-comparison-tables)
7. [Results tables from DESeq2, edgeR, limma and Seurat](#7-results-tables-from-deseq2-edger-limma-and-seurat)
8. [F-tests and their direction](#8-f-tests-and-their-direction)
9. [Gene-burden tests](#9-gene-burden-tests)
10. [What each column is used for](#10-what-each-column-is-used-for)
11. [What is refused](#11-what-is-refused)
12. [Error codes for input](#12-error-codes-for-input)

## 1. The three shapes of input

`TABLE`, given to `evidome TABLE --out DIR`, `evidome crosswalk` or
`evidome inspect`, may be:

| Shape | What it is | How Evidome reads it |
|---|---|---|
| **Wide table** | One file, one row per gene, a block of columns per comparison, each column prefixed with the comparison's name (section 4). Paralome writes this. | As it is. |
| **Per-comparison table** | One comparison per file, the same column names without the prefix (section 6). A single file, or a folder of them. | Combined into one wide table first, `DIR/input_combined.csv`; the comparison is named after the file. |
| **A tool's results table** | DESeq2, edgeR, limma or Seurat FindMarkers results, written from R as they come (section 7). A single file, or a folder of them. | Read by the tool's own column names and combined the same way. |

A file is taken for a single comparison when its columns are not prefixed
(it has a `P-value` or `ANOVA_p` column), or when it is one of the tools'
results tables. Otherwise it is read as a wide table.

`evidome combine FOLDER --out FILE` writes the combined wide table on its own,
for inspection or reuse.

## 2. Text, delimiters and missing values

- **Encoding:** UTF-8. A byte-order mark (Excel's "CSV UTF-8") is accepted. Any
  other encoding is refused (EVD-0007). Excel workbooks (`.xlsx`) are not read;
  save the sheet as CSV.
- **Delimiter:** tab when the column header line contains a tab, otherwise
  comma. Space-separated files are not read.
- **Quoting:** standard CSV quoting (`"a,b"`) is understood, as R's
  `write.csv` writes it.
- **Comment lines:** lines starting with `#` before the column header are the
  comment header (section 5). Blank lines before the header are skipped.
- **Missing values:** an empty cell, `.`, `NA`, `NaN` or `nan`. A gene whose
  p-value is missing in a comparison was **not tested** in it: it is left out of
  that comparison's universe, never treated as p = 1.
- **Numbers:** any form Python reads as a number (`0.001`, `1e-300`, `-2.5`).
  `Inf` and `-Inf` are accepted only in a fold-change column (a gene switched
  on from zero has a real, signed answer). An infinite p-value, q-value,
  T-statistic or effect is refused (EVD-0007): leave the cell empty instead.
- **Rows:** every row must have as many fields as the header (EVD-0007). A row
  with no gene ID (empty, `.` or `NA`) is left out and counted in the notes.

## 3. Gene identifiers

Every table has one gene ID column, holding one of three kinds of ID:

| Kind | Example | Species |
|---|---|---|
| Ensembl stable gene ID | `ENSMMUG00000000001`, `ENSG00000139618.17` | found from the IDs |
| NCBI GeneID ("Entrez ID") | `7157` | found from the IDs |
| Gene symbol | `TP53`, `Gapdh` | **must be named with `--species`** |

- **Column name:** one of `Gene`, `gene_id`, `Ensembl_ID` or
  `ensembl_gene_id`, whatever kind of ID it holds. Any other name is given
  with `--id-column NAME`. In a per-comparison or tool table, R's unnamed
  row-name column is also taken as the ID column (section 7).
- **Which kind:** read from the column itself. If 90% or more of the IDs are
  whole numbers, they are NCBI GeneIDs. If at least half look like Ensembl
  stable gene IDs (`ENS...G...`, `FBgn...`, `WBGene...`), they are Ensembl
  IDs. Otherwise they are symbols; with `--species`, IDs that are mostly that
  species' own Ensembl gene IDs in another form (yeast's `YAL001C`) are read as
  Ensembl IDs. `evidome inspect` says which kind it found.
- **One kind per table.** An ID of another kind in the column is simply one
  that is not found.
- **Each gene appears once** (EVD-0004). Ensembl version suffixes are stripped
  first: `ENSG00000139618.16` and `ENSG00000139618.17` are the same gene.

### Ensembl IDs

Used as they are. The species is the one whose downloaded gene list holds at
least 90% of the IDs (`--min-species-match`). A species not yet downloaded is
looked up with Ensembl and fetched the first time it is seen
(`--no-fetch-species` turns that off). Any species Ensembl annotates can be
used, human included.

### Symbols and NCBI GeneIDs

These are translated, gene by gene, to the species' Ensembl genes before
anything else happens; the rest of the run is the same as for Ensembl IDs.

- **Gene symbols need `--species`**, because a symbol does not say which
  species it belongs to: macaque and human symbols are nearly all the same.
  Give the species as Ensembl names it: `--species macaca_mulatta`,
  `--species mus_musculus`, `--species homo_sapiens` (`"Macaca mulatta"` is
  read the same way). Without it the run stops with EVD-0032.
- **A symbol** is looked up in the species' Ensembl annotation (the `gene_name`
  of its GTF, for the Ensembl release in use): an exact match first; failing
  that, a match that ignores case, if it is the only one. Only current
  symbols are matched; an outdated symbol is not in the annotation.
- **An NCBI GeneID** is looked up in NCBI's `gene2ensembl`, which is
  downloaded once and kept with the reference data. Each GeneID belongs to
  one taxon, so the species is found from the IDs: the taxon that most of them
  belong to. `--species` may be given to name it instead.
- **Never guessed.** An ID that names no gene of the species, or more than one,
  is left untranslated, and so are two input IDs that would become the same
  Ensembl gene. An untranslated gene is not dropped: it stays in every output
  under the ID it was given, and the crosswalk gives the reason as its
  `mapping_status` (`not translated: symbol not in the macaca_mulatta
  annotation`, `not translated: symbol names several macaca_mulatta Ensembl
  genes`, ...). Its symbol can still reach a human gene through the symbol
  match below.
- **What was done is recorded.** Every crosswalk row has `input_id_given`, the
  ID exactly as the table gave it, and `id_translation`: `symbol, exact`,
  `symbol, ignoring case`, `NCBI GeneID (gene2ensembl)`, `Ensembl gene ID (as
  given)`, or, for an untranslated ID, the genes it named.
- `--species` given with Ensembl IDs or NCBI GeneIDs must agree with the
  species the IDs belong to (EVD-0014).

### Symbols beside the IDs

A separate symbol column is optional: one of `Common_name`, `gene_name`,
`Symbol` or `symbol`. A symbol equal to the gene ID means "no symbol". When
the IDs are symbols, they are the genes' symbols too. Symbols are carried to
the outputs and, for genes no ortholog source maps, used for an exact match to
a human gene that is the only one with that symbol (turned off with
`--no-symbol-match`).

### Orthology to human

In order: Ensembl Compara gene trees, then NCBI Gene orthologs (`--no-ncbi`
turns these off), then the symbol match above. Every mapping, and every gene
left without one, is listed in `ortholog_crosswalk.csv` with its evidence.

<!-- example: symbols shape=wide -->
```
Gene,WD_vs_CD_Log2_fold_change,WD_vs_CD_P-value,WD_vs_CD_Q-value
TSPAN6,2.0,0.001,0.01
PARA1,-1.2,0.02,0.10
H4C1,0.1,0.80,0.90
```

Run with `--species macaca_mulatta` (or whichever species the data is from).

## 4. The wide table

One row per gene. Columns: the gene ID, an optional symbol, and one **block**
per comparison. Every column of a block is named `<comparison>_<role>`.

<!-- example: minimal-wide shape=wide -->
```
Gene,Common_name,WD_vs_CD_Log2_fold_change,WD_vs_CD_P-value,WD_vs_CD_Q-value
ENSMMUG00000000001,TSPAN6,2.0,0.001,0.01
ENSMMUG00000000002,,-1.2,0.02,0.10
ENSMMUG00000000003,,0.1,0.80,0.90
```

### Comparison names

A comparison's name is everything before the role suffix: in
`diet_WD_vs_CD_P-value` the comparison is `diet_WD_vs_CD`. Names may contain
underscores. A block is found by its `_P-value` column (a contrast) or its
`_ANOVA_p` column (an F-test); its other columns are looked up by the same
name.

### A contrast (two groups, with a direction)

| Role | Required | Meaning |
|---|---|---|
| `P-value` | yes | Raw p-value of the test. |
| `Q-value` | yes | Adjusted p-value (any method; carried, see section 10). |
| `Log2_fold_change` | exactly one of these two | log2 fold change: 0 unchanged, +1 doubled, -1 halved. Always read as log2. |
| `Fold_change` | | A fold change in the producer's own convention, judged from the column's values (below). |
| `Control_mean`, `Treatment_mean` | optional, as a pair | Group means, any unit. |
| `Raw_FPKM_Control_mean`, `Raw_FPKM_Treatment_mean` | optional, as a pair | Group means in FPKM (Paralome bulk RNA). |
| `T-statistic` | optional | The test's t (or moderated t). |
| `Effect` | optional | The model's effect estimate, log2 (Paralome's model contrasts). |

**Positive means up in the treatment group.** A fold change, T-statistic or
effect above zero says the gene is higher in the treatment: the arm named by
`treatment=` in the comment header (section 5), or the numerator of the
comparison.

**How a `Fold_change` column is read.** The whole column is judged once:

| The column holds | Read as | Unchanged | Example |
|---|---|---|---|
| no negative values, and some between 0 and 1 | a ratio | 1 | 0.5 = halved, 2 = doubled |
| negative values, and some between -1 and 1 | a log fold change | 0 | -1 = halved (log2) |
| neither of the above | signed "times" (Paralome bulk) | 1 or -1 | -2 = halved, +2 = doubled |

A `Log2_fold_change` column is never judged this way. Use it whenever the
values are log2: a column of small positive log2 values (every gene a little
up) would otherwise read as ratios below 1, that is, down.

Giving both `Fold_change` and `Log2_fold_change` in one block is refused
(EVD-0006), as is half of a pair of means.

### An F-test (several groups, no direction of its own)

| Role | Required | Meaning |
|---|---|---|
| `F` | yes | The F statistic. |
| `ANOVA_p` | yes | Raw p-value. |
| `ANOVA_q` | yes | Adjusted p-value. |

An F-test gets a direction only from its own contrasts in the same table
(section 8).

<!-- example: full-wide shape=wide -->
```
# generated by my_pipeline 2.1
#   WD_vs_CD: pairwise, control=CD, treatment=WD
#   diet: ANOVA F-test
Gene,Common_name,WD_vs_CD_Control_mean,WD_vs_CD_Treatment_mean,WD_vs_CD_Fold_change,WD_vs_CD_T-statistic,WD_vs_CD_P-value,WD_vs_CD_Q-value,diet_F,diet_ANOVA_p,diet_ANOVA_q
ENSMUSG00000000001,Gnai3,10.2,20.5,2.01,5.1,0.001,0.01,26.0,0.001,0.01
ENSMUSG00000000003,Pbsn,3.1,1.4,-2.21,-3.9,0.004,0.03,15.2,0.004,0.03
ENSMUSG00000000028,Cdc45,5.0,5.2,1.04,0.3,0.77,0.88,0.1,0.77,0.88
ENSMUSG00000000031,H19,0.0,4.0,Inf,6.2,0.0005,0.008,38.4,0.0005,0.008
```

### Every column must be claimed

Each column must be the gene ID, the symbol, or part of a block. Any other
column (sample values, annotations, extra statistics) is refused, by name
(EVD-0008): an unexplained column is one whose meaning an AI agent reading the
evidence would have to guess. Remove it, or use a per-comparison table, where
other columns are left out by design (section 6). A column name may appear
only once.

## 5. The comment header

Optional `#` lines before the column header describe the table. Each line is
one of:

| Line | Meaning |
|---|---|
| `# generated by <program and version>` | Who wrote the table; recorded in the evidence. |
| `#   <name>: control=<A>, treatment=<B>` | Names the two arms of contrast `<name>`. |
| `#   <name>: <description>, control=<A>, treatment=<B>` | The same, with a description first. |
| `#   <name>: control=<A>, treatment=<B>, test=burden` | Contrast-shaped block that is a gene-burden test (section 9). |
| `#   <name>: ANOVA F-test ...` | Block `<name>` is an F-test. |
| `#   <name>: <any other text>` | A free description of `<name>`. |
| any other `#` line | Ignored. |

The arms are recorded in the evidence, so an agent knows which way "up"
points. A header line that names a comparison the table has no columns for,
or calls a block a different kind than its columns show, is refused
(EVD-0009). `test=` with any value other than `burden` is refused, never
ignored. A comparison the header does not mention is read from its columns,
with a note.

## 6. Per-comparison tables

One comparison per file. The columns are the roles of section 4 **without**
the name prefix, and the comparison is named after the file:

<!-- example: WD_vs_CD shape=single -->
```
# WD_vs_CD: control=CD, treatment=WD
Gene,Common_name,Control_mean,Treatment_mean,Log2_fold_change,P-value,Q-value,S1,S2,S3,S4
ENSMMUG00000000001,TSPAN6,1.0,4.0,2.0,0.001,0.01,1,1,4,4
ENSMMUG00000000002,,3.0,1.3,-1.2,0.02,0.10,3,3,1,1
ENSMMUG00000000003,,2.0,2.1,0.1,0.80,0.90,2,2,2,2
```

- **Name:** the file name without its extension (`WD_vs_CD.csv` is
  `WD_vs_CD`). Letters, digits and `_ . + -` only, and unique across the
  folder. `--drop TEXT` (repeatable) removes text shared by every file name,
  e.g. `--drop _naive_cell_t_rna`.
- **Which files:** in a folder, every `*.csv`, `*.tsv` and `*.txt` file by
  default. `--pattern` chooses others: one glob or several separated by commas
  (`--pattern '*.tsv'`, `--pattern 'WD_*.csv,CD_*.csv'`). Every chosen file
  must be a table: a `README.txt` in the folder is read as one, and refused.
  Hidden files and tables Evidome wrote itself are skipped.
- **Other columns** (per-sample values, annotations) are left out: only the
  role columns, the ID and the symbol are kept.
- **Statistics are copied as text,** never recomputed or reformatted.
- **The file's own `#` header** is kept: who wrote it, and its control and
  treatment when it names one contrast.
- A gene missing from one file is untested in that comparison.

## 7. Results tables from DESeq2, edgeR, limma and Seurat

A results table written from R as it comes is read by its own column names,
one comparison per file. Only the columns below are used; the tool's other
columns are left out.

| Tool | Fold change (log2) | P-value | Adjusted p (Q-value) | Also used |
|---|---|---|---|---|
| DESeq2 `results()` | `log2FoldChange` | `pvalue` | `padj` | |
| edgeR `topTags()` | `logFC` | `PValue` | `FDR` | |
| limma `topTable()` | `logFC` | `P.Value` | `adj.P.Val` | `t` as the T-statistic |
| Seurat `FindMarkers()` | `avg_log2FC` | `p_val` | `p_val_adj` | |

Left out: DESeq2's `baseMean`, `lfcSE`, `stat`; edgeR's `logCPM`, `LR`, `F`;
limma's `AveExpr`, `B`; Seurat's `pct.1`, `pct.2`; and any annotation column
other than a recognised symbol column.

**Gene IDs** may be Ensembl IDs, NCBI GeneIDs or gene symbols (section 3;
symbols need `--species`). R writes them as row names:
`write.csv` writes the row-name column with an empty name, and
`write.table` leaves it out of the header line altogether. Both are read as
the ID column. A named column (`gene_id`, or one given with `--id-column`) is
used instead when there is one.

**Writing the files from R,** one comparison per file, named for the
comparison:

```r
# DESeq2
write.csv(as.data.frame(results(dds, contrast = c("diet", "WD", "CD"))), "WD_vs_CD.csv")
# edgeR
write.csv(topTags(glmQLFTest(fit, coef = "dietWD"), n = Inf)$table, "WD_vs_CD.csv")
# limma
write.csv(topTable(fit2, coef = "WD_vs_CD", number = Inf), "WD_vs_CD.csv")
# Seurat (features are usually symbols: run with --species)
write.csv(FindMarkers(obj, ident.1 = "WD", ident.2 = "CD",
                      logfc.threshold = 0, min.pct = 0), "WD_vs_CD.csv")
```

<!-- example: WD_vs_CD shape=single -->
```
"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"ENSMMUG00000000001",120.5,2.0,0.3,6.7,0.001,0.01
"ENSMMUG00000000002",88.1,-1.2,0.4,-3.0,0.02,0.1
"ENSMMUG00000000003",45.0,0.1,0.5,0.2,0.8,0.9
```

Things to know:

- **Which way is up:** the sign is the tool's. Positive is higher in DESeq2's
  contrast numerator (`"WD"` above), in edgeR's and limma's coefficient, and in
  Seurat's `ident.1`. To record the arms, put a header line at the top of the
  file before the column names: `# WD_vs_CD: control=CD, treatment=WD`.
- **Every tested gene, not only the significant ones.** Over-representation
  uses the genes in the file as the universe of tested genes, so write the full
  table: `n = Inf` (edgeR), `number = Inf` (limma), and for Seurat
  `logfc.threshold = 0, min.pct = 0`, which otherwise drop genes before
  testing.
- **DESeq2's NA values:** a gene with `pvalue` NA (all counts zero, or an
  outlier) is untested. A gene with a `pvalue` but `padj` NA (independent
  filtering) is tested; its Q-value is missing.
- **Seurat's `p_val_adj`** is Bonferroni-adjusted; it is carried as the
  Q-value as given.
- **One coefficient per file.** An edgeR or limma table holding several
  coefficients (`logFC.A`, `logFC.B`, ...) is not one of these formats: write
  each coefficient to its own file.
- **Not read:** Seurat's `FindAllMarkers()` output (many clusters in one
  file; write one `FindMarkers()` per comparison instead), and Seurat before
  version 4, whose `avg_logFC` is a natural-log fold change.
- A table that has Evidome's own `P-value` column is read by Evidome's names,
  whatever tool columns it also holds.

## 8. F-tests and their direction

An F-test has no direction, so over-representation tests it once ("all")
and GSEA skips it, unless its own contrasts in the same table give it one.
Two name patterns are recognised (as Paralome writes them):

- **Main effect** `<factor>_<alt>_vs_<ref>` (or `<alt>_vs_<ref>`): takes the
  direction of the contrast `<alt>_vs_<ref>_main_effect`, from its `Effect`
  where given, else its `T-statistic`. With two levels this F-test and that
  contrast are the same test, so by default the F-test is not analysed and the
  contrast stands in for it (`--keep-duplicate-ftests` analyses both).
- **Interaction** `<A>_x_<B>`: Effect of the A contrast within B's alternative
  level minus Effect within B's reference level (contrasts named
  `<A stem>_within_<level>_contrast`, each with an `Effect` column).

An F-test whose name matches neither pattern, or whose contrasts are not in
the table, is analysed without a direction.

## 9. Gene-burden tests

A contrast-shaped block whose header line ends `, test=burden` is a gene-burden
test: observed carriers of qualifying variants (treatment, e.g. probands)
against the number expected from a reference (control, e.g. gnomAD).

```
#   dominant_ptv: control=gnomAD, treatment=probands, test=burden
```

Its columns are a contrast's, usually `Log2_fold_change` (log2 observed /
expected), `P-value` and `Q-value`, with no means. Up means more carriers than
expected. A burden test, and only a burden test, may carry
`<name>_Null_selection_probability`: for each gene, the chance under the
test's null of reaching p < 0.05. When present, and the run uses the default
`--alpha 0.05`, over-representation weights its null by it; at another alpha
the unweighted null is used and the results say so. It must be given for every
tested gene, within [0, 1].
No GSEA and no activation z-scores are computed for a burden test.

## 10. What each column is used for

| Column | Used for |
|---|---|
| Gene ID | Species identification, translation to Ensembl genes (symbols, NCBI GeneIDs), orthology to human, the key of every output row. |
| Symbol | Carried; a last-resort orthology match (section 3). |
| `P-value` / `ANOVA_p` | Over-representation: a gene is selected when raw p < `--alpha` (default 0.05). GSEA ranking (sign x -log10 p) when there is no T-statistic. A vote when several input genes map to one human gene. |
| `Q-value` / `ANOVA_q` | Carried to the outputs. Not used in any statistic and not a vote: it is a monotone transform of the p-value. |
| `T-statistic` | The gene's direction (first choice) and the GSEA ranking metric. A collision vote. |
| `Effect` | The direction when there is no T-statistic; an F-test's direction (section 8). A collision vote. |
| `Log2_fold_change` / `Fold_change` | The direction when there is neither T-statistic nor Effect. A collision vote, on its magnitude. |
| Means | Carried to the outputs. Never enter a statistic. |
| `F` | Carried. A collision vote for F-tests. |
| `Null_selection_probability` | The weights of a burden test's over-representation null (section 9). |

When several input genes map to the same human gene (paralogs), each column
above marked "vote" casts a vote for the input gene it ranks best, and the
majority decides which input gene stands for the human gene in that
comparison. Every ballot is written to `collision_ballots.csv`.

## 11. What is refused

| Refused | Why |
|---|---|
| Gene symbols without `--species` | A symbol does not say which species it belongs to. |
| A table mixing species | Each table is translated and mapped in one species; split it. |
| A gene listed twice | Which row is the gene would be a guess. |
| A column no block claims (wide table) | An agent would have to guess its meaning. |
| `Fold_change` and `Log2_fold_change` in one block | Two readings of one number. |
| Half a pair of means | The pair is given whole or not at all. |
| An infinite statistic other than a fold change | It would turn every score it touches into nan. |
| A header that disagrees with the columns | The table was probably edited by hand. |
| `Null_selection_probability` outside a burden test | It has a meaning only for a burden test's null. |
| A non-UTF-8 file, `.xlsx`, a space-separated file | Not read as text in a known form. |

## 12. Error codes for input

| Code | Meaning |
|---|---|
| EVD-0003 | No gene ID column (or the one named with `--id-column` is absent). |
| EVD-0004 | A gene ID appears more than once. |
| EVD-0005 | No comparison block found (no `*_P-value` or `*_ANOVA_p` column). |
| EVD-0006 | A block lacks a required column, or has both fold-change columns. |
| EVD-0007 | A row with the wrong number of fields, a value that is not a number, an infinite statistic, or a file that is not UTF-8. |
| EVD-0008 | A column no block claims, or a column name used twice. |
| EVD-0009 | The comment header disagrees with the columns. |
| EVD-0014 | The species could not be identified from the IDs, or `--species` contradicts them. |
| EVD-0030 | The table or folder cannot be read. |
| EVD-0031 | Per-comparison tables cannot be combined: no ID column, not one comparison's table, or names that are empty, repeated or not letters, digits and `_ . + -`. |
| EVD-0032 | The gene IDs are symbols and no `--species` was given. |

Every error prints its code, the file and, where it applies, the line or
column it is about, and a sentence saying how to fix it.
