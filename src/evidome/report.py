"""The evidence bundle for an interpreting agent, and a workbook for people.

``report <crosswalk dir>`` gathers what ``crosswalk`` and ``enrich`` wrote into
one self-describing folder, ``<crosswalk dir>/evidence/``:

    manifest.json               provenance, file index with every column
                                described, glossary, reading guide
    comparisons.csv             one row per comparison: design and counts
    findings.csv                every enrichment and regulator result with q < 0.05 (index)
    enrichment_ora.csv          every over-representation result
    enrichment_gsea.csv         every GSEA result
    regulators.csv              every upstream regulator tested
    human_gene_statistics.csv   one row per (comparison, human gene)
    ortholog_crosswalk.csv      every input gene -> human, with evidence
    collision_ballots.csv       every many-to-one vote, one row per candidate
    gene_sets.csv               id, source, name, category, size of each set used

plus ``<crosswalk dir name>_combined_results.zip`` (the same folder, compressed, unzipping to a
folder of that name, with ``<Name>_Viewer_v<version>.html`` beside the files: a page for people that
opens the zip it came in and browses it, named and stamped with the program version that wrote it)
and a tabbed workbook
``<crosswalk dir name>_workbook.xlsx``. The CSVs are the record; the workbook
is a view of them, and says on its README sheet what it leaves out.

The program states what each number is and how it was made -- never what it
means biologically. The reading guide is written in that spirit:
definitions, directions and caveats, no conclusions.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .xlsx import FIXED_ZIP_TIME

CONTRAST_STATS = ("Raw_FPKM_Control_mean", "Raw_FPKM_Treatment_mean", "Control_mean", "Treatment_mean",
                  "Fold_change", "Log2_fold_change", "Effect", "T-statistic", "P-value", "Q-value",
                  "Null_selection_probability")
#: Statistic columns kept in the long tables only when some comparison has them (alternatives and optional roles).
OPTIONAL_STATS = ("Control_mean", "Treatment_mean", "Raw_FPKM_Control_mean", "Raw_FPKM_Treatment_mean",
                  "Direction_estimate", "Fold_change", "Log2_fold_change", "Null_selection_probability")
FTEST_STATS = ("F", "ANOVA_p", "ANOVA_q", "Direction_estimate")

# ---------------------------------------------------------------------------
# one description per column, used by the manifest AND the workbook
# ---------------------------------------------------------------------------

COLUMN_TEXT: Dict[str, str] = {
    # identity
    "comparison": "Name of the comparison, as in the input table's column prefix.",
    "kind": "'contrast' (two groups, has a direction), 'burden' (a gene-burden test: carriers observed in the "
            "treatment group, e.g. probands, against the number expected from the control, e.g. gnomAD allele "
            "frequencies; has a direction, up = more than expected) or 'ftest' (ANOVA F-test; it has a direction "
            "only when its own contrasts give it one -- see direction_from).",
    "control": "Control group of a contrast, from the input header; for a burden test, the reference the "
               "expected count comes from (e.g. gnomAD).",
    "status": "'analysed', or why a comparison in the input was not: a main-effect F-test is the same test as its "
              "main-effect contrast (two levels per factor), which is analysed in its place.",
    "direction_from": "F-tests only: the contrast(s) the F-test's direction was taken from (its main-effect "
                      "contrast; for an interaction, the difference of its two simple-effect contrasts). Blank: "
                      "a contrast (its own direction), or an F-test whose contrasts are not in the input.",
    "Direction_estimate": "F-tests only: the signed estimate giving the F-test a direction, from its contrasts "
                          "(see comparisons.csv direction_from): main effect -- that contrast's Effect (log2) or "
                          "T-statistic; interaction -- Effect within one level minus Effect within the other. "
                          "Derived by this program; not in the input.",
    "treatment": "Treatment group of a contrast, from the input header; for a burden test, the cohort whose "
                 "carriers were counted (e.g. probands).",
    "human_gene_id": "Human Ensembl gene ID (release in manifest provenance).",
    "human_symbol": "Human gene symbol from the same Ensembl release.",
    "human_biotype": "Human gene biotype (protein_coding, lncRNA, ...).",
    "input_gene_id": "The input species' Ensembl gene ID, version suffix removed (translated from a symbol or "
                     "NCBI GeneID where the input gave one; the given ID where it could not be translated).",
    "input_symbol": "The input table's symbol for the gene, if it had one.",
    "input_biotype": "Biotype of the input gene in the annotation it was quantified against.",
    "species": "Ensembl species of the input genes.",
    "source_gene_id": "Input gene whose statistics represent this human gene in this comparison.",
    "source_symbol": "Symbol of that input gene.",
    "source_biotype": "Biotype of that input gene.",
    "candidate_gene_id": "An input gene competing for this human gene.",
    "candidate_symbol": "Its symbol.",
    # orthology evidence
    "homology_source": "Where the input->human pair came from: protein / ncrna (Ensembl Compara gene trees, "
                       "strongest), ncbi (NCBI Gene orthologs), symbol (exact unambiguous symbol, weakest).",
    "homology_type": "ortholog_one2one / one2many / many2many (Compara), 'ortholog (NCBI)', or 'symbol match'.",
    "is_high_confidence": "Compara's high-confidence flag (1/0); empty for NCBI and symbol pairs.",
    "identity_pct": "Percent of the input protein identical to the human one (Compara).",
    "human_identity_pct": "Percent of the human protein identical to the input one (Compara).",
    "goc_score": "Compara gene-order conservation score (0-100), empty if not computed.",
    "wga_coverage": "Compara whole-genome-alignment coverage (0-100), empty if not computed.",
    "human_targets_of_input_gene": "How many human genes this input gene maps to (>1 = one-to-many).",
    "input_genes_sharing_human_gene": "How many input genes map to this human gene (>1 = collision).",
    "mapping_status": "mapped / no human ortholog in reference / in registered genome, not in Ensembl "
                      "reference / gene ID not in reference / human input (identity) / not translated: <why> "
                      "(an input symbol or NCBI GeneID that named no gene of the species, several, or the same "
                      "gene as another input ID).",
    "input_id_given": "The gene's ID exactly as the input table gave it: an Ensembl ID (with any version), an "
                      "NCBI GeneID or a gene symbol.",
    "id_translation": "How the given ID became input_gene_id: 'Ensembl gene ID (as given)', 'symbol, exact', "
                      "'symbol, ignoring case' (the only gene of the species with that symbol in another case) or "
                      "'NCBI GeneID (gene2ensembl)'; for an untranslated ID, the genes it named, if any.",
    "source_human_targets": "Human genes the source gene maps to; >1 means its statistics also represent others.",
    "n_candidates": "Input genes mapping to this human gene.",
    "n_tested_candidates": "Of those, how many were tested in this comparison.",
    "candidates": "All candidate input gene IDs, ';'-separated.",
    "decided_by": "How the source gene was chosen: single candidate, vote, tie-break, only tested candidate, "
                  "Compara precedence, none tested.",
    "tie": "yes if the vote tied and a tie-break decided.",
    "tested": "yes if the gene had a p-value in this comparison.",
    "significant": "yes if raw p < alpha (manifest settings) -- a counting threshold, nothing is filtered by it.",
    "votes": "Metrics this candidate won (majority vote of p, |t|, |fold change|, |effect|; or p and F).",
    "metrics_won": "Which metrics it won, ';'-separated.",
    "winner": "yes for the candidate that represents the human gene in this comparison.",
    "excluded": "Why a candidate could not vote (a Compara candidate was present), else empty.",
    # statistics
    "Raw_FPKM_Control_mean": "Mean FPKM in the control group (input table).",
    "Raw_FPKM_Treatment_mean": "Mean FPKM in the treatment group (input table).",
    "Control_mean": "Mean in the control group, in the unit the input table's header states (not FPKM; "
                    "e.g. pseudobulk counts or accessibility). Blank for inputs that give FPKM means.",
    "Treatment_mean": "Mean in the treatment group, in the unit the input table's header states (not FPKM).",
    "Fold_change": "Fold change as the input reports it: Paralome bulk signed 'times' (+2 up, -2 down, 1 "
                   "unchanged), Paralome single-cell log2 (+1 up, -1 down, 0 unchanged), or a plain ratio (0.5 down).",
    "Log2_fold_change": "log2 fold change, declared as log2 by the input (0 unchanged, +1 doubled, -1 halved); "
                        "for a burden test, log2(observed / expected carriers). Read as given, never re-judged.",
    "Effect": "Model effect size (log2 scale) from the input, where present.",
    "T-statistic": "t statistic of the contrast (sign = direction), where present.",
    "P-value": "Raw p-value of the contrast.",
    "Q-value": "FDR q-value of the contrast, as computed by the input's producer.",
    "Null_selection_probability": "Burden tests only, from the input: the probability, under the burden test's "
                                  "null, that the gene reaches p < 0.05 -- its detectability. The weights of the "
                                  "burden ORA's null (null_model).",
    "F": "ANOVA F statistic.",
    "ANOVA_p": "ANOVA p-value.",
    "ANOVA_q": "ANOVA FDR q-value.",
    # enrichment
    "method": "ORA (over-representation of significant genes), GSEA (ranked-list enrichment), "
              "REGULATOR_ACTIVITY (activity_t) or REGULATOR_OVERLAP (overlap_p).",
    "direction": "ORA: which significant genes were tested -- up, down, or all. Findings: for GSEA, the sign "
                 "of NES (up = towards the top of the ranking = higher in treatment); for REGULATOR_ACTIVITY, "
                 "the sign of activity_t; REGULATOR_OVERLAP rows are 'all'.",
    "source": "Gene-set collection (reactome, go_bp, go_mf, go_cc, wikipathways, hallmark); in findings, the "
              "regulator network for REGULATOR_* rows.",
    "set_id": "Gene-set identifier (REACTOME:R-HSA-..., GO:..., WP:WP..., MSIGDB:HALLMARK_...); in findings, "
              "<network>:<regulator> for REGULATOR_* rows.",
    "name": "Gene-set name; in findings, the regulator for REGULATOR_* rows.",
    "category": "Reactome: top-level pathway(s); GO: namespace; else empty. In findings, the regulator kind for "
                "REGULATOR_* rows.",
    "set_size": "Genes of the set among the genes tested (ORA universe / GSEA ranked list).",
    "universe": "Human genes tested in the comparison: the ORA background.",
    "selected": "Significant genes in this direction: the ORA draw.",
    "overlap": "Selected genes in the set.",
    "expected": "Overlap expected by chance: set_size x selected / universe.",
    "fold_enrichment": "overlap / expected.",
    "p": "ORA: exact one-sided hypergeometric p. GSEA: permutation p, (exceeding nulls + 1) / (nulls + 1).",
    "q": "Benjamini-Hochberg q within (comparison, direction, source) for ORA, (comparison, source) for GSEA.",
    "overlap_genes": "Human symbols of the overlapping genes, ';'-separated.",
    "pattern_genes": "ORA, Reactome only: the set's tested genes with an activity pattern (expected up or down "
                     "when the pathway is active, from Reactome's reaction logic).",
    "pattern_consistent": "Significant pattern genes whose direction matches the pattern.",
    "pattern_inconsistent": "Significant pattern genes whose direction opposes the pattern.",
    "pathway_activation_z": "Pathway activation z: (consistent - inconsistent) / sqrt(N) over significant "
                            "pattern genes; blank when N < 4, no pattern, or no direction. Positive: the pathway "
                            "looks active; a negative regulator going down counts toward activation.",
    "pattern_consistent_genes": "Symbols of the consistent pattern genes.",
    "null_model": "The null the ORA p comes from: 'hypergeometric'; for a burden test with "
                  "Null_selection_probability, 'Wallenius noncentral hypergeometric, weighted by detectability' "
                  "(odds = mean selection probability of the set's genes / of the rest of the universe, as goseq); "
                  "for a burden test without it, 'hypergeometric (unadjusted for detectability)'.",
    "overlap_null": "The null overlap_p comes from, as null_model for ORA.",
    "ora_null": "The comparison's ORA null (see null_model).",
    "direction_statistics": "Burden tests: why no direction-of-effect statistic (activation z, activity t, "
                            "pathway activation z, regulator direction) is reported for them.",
    "pattern_inconsistent_genes": "Symbols of the inconsistent pattern genes.",
    "es": "GSEA enrichment score (weighted running sum maximum, Subramanian 2005).",
    "nes": "Normalised enrichment score: ES / mean |null ES| of the same sign. Positive = up, negative = down.",
    "permutations": "Permutations the GSEA p rests on (adaptive: more where the set was extreme).",
    "leading_edge_size": "Genes in the GSEA leading edge.",
    "leading_edge_genes": "Leading-edge human symbols, ';'-separated, in rank order.",
    "rank_metric": "What GSEA ranked genes by (T-statistic, or sign x -log10 p).",
    "effect": "Findings: fold_enrichment (ORA), NES (GSEA), activity_t (REGULATOR_ACTIVITY), activation_z "
              "(REGULATOR_OVERLAP, contrasts).",
    "genes": "Findings: overlap genes (ORA), leading-edge genes (GSEA), targets that moved the way the activity "
             "says (REGULATOR_ACTIVITY: consistent targets when up, inconsistent when down), "
             "significant signed targets (REGULATOR_OVERLAP; blank for a comparison without a direction, whose "
             "targets carry no sign).",
    "rank": "Findings: rank by q, then p, within (comparison, method, direction, source).",
    "n_genes": "Genes in the set after mapping to human Ensembl genes of this release.",
    "source_ids": "Members the source listed before mapping (difference = members not in this release).",
    # comparisons table
    "alpha": "Raw-p threshold used for counting significant genes in this table (crosswalk --alpha; enrich and "
             "regulators use it too unless given their own -- see enrichment_settings and regulator_settings).",
    "input_genes_tested": "Input genes with a p-value in this comparison.",
    "input_genes_significant": "Tested input genes with raw p < alpha.",
    "significant_mapped": "Significant input genes that reached a human gene (any source).",
    "significant_via_protein": "... of which through Compara protein trees.",
    "significant_via_ncrna": "... through Compara ncRNA trees.",
    "significant_via_ncbi": "... through NCBI Gene orthologs only.",
    "significant_via_symbol": "... through an exact symbol only.",
    "significant_mapped_but_outvoted": "Significant genes that lost every collision they were in.",
    "significant_excluded_by_compara_precedence": "Significant genes set aside because a Compara candidate "
                                                  "competed for the same human gene.",
    "significant_unmapped": "Significant input genes that reached no human gene.",
    "significant_unmapped_by_biotype": "Their biotypes, 'biotype=count;...'.",
    "significant_unmapped_by_status": "Their mapping statuses, 'status=count;...'.",
    "human_genes_reached": "Human genes reached by any input gene.",
    "human_genes_tested": "Human genes whose representative gene was tested.",
    "human_genes_significant": "Human genes whose representative gene has raw p < alpha.",
    "human_genes_with_tested_collisions": "Human genes with >= 2 tested candidates.",
    "ora_universe": "ORA background size (= human genes tested).",
    "ora_selected": "Significant human genes used by ORA (all directions).",
    "gsea_rank_metric": "GSEA ranking metric, empty for F-tests without a direction.",
    "gsea_permutations": "Total GSEA permutations run for this comparison.",
    "ora_sets_q_below_0_05": "ORA results with q < 0.05, 'source/direction=count;...'.",
    "gsea_sets_q_below_0_05": "GSEA results with q < 0.05, 'source=count;...'.",
    "regulators_q_below_0_05": "Regulators with q < 0.05, 'network/overlap|activity=count;...'.",
    # regulators
    "network": "collectri (transcription factor -> target genes, direct) or signalling_via_tf (signalling "
               "protein -> transcription factor -> its targets, INDIRECT).",
    "regulator": "Regulator gene symbol, or a TF complex written as its members (e.g. JUN_FOS).",
    "regulator_kind": "tf, complex (TF complex) or signalling (upstream protein, indirect network).",
    "targets_in_universe": "The regulator's target genes among the genes tested in the comparison.",
    "targets_signed": "Of those, targets with a known direction of effect (activation or repression).",
    "significant_targets": "Targets with raw p < alpha in the comparison.",
    "overlap_p": "Exact one-sided hypergeometric p: targets over-represented among significant genes.",
    "overlap_q": "BH q of overlap_p within (comparison, network).",
    "n_consistent": "Significant signed targets changing as the regulator's ACTIVATION would predict.",
    "n_inconsistent": "Significant signed targets changing the opposite way.",
    "activation_z": "(n_consistent - n_inconsistent) / sqrt(their sum): positive fits activation of the regulator, "
                    "negative its inhibition (without bias correction).",
    "z_p": "Exact two-sided binomial p of the consistent / inconsistent split.",
    "bias_term": "mu = mu_data x mu_regulator: how much the comparison's up/down imbalance and the regulator's "
                 "activating/inhibiting imbalance alone would push the z (Kraemer et al. 2014).",
    "activation_z_bias_corrected": "Activation z corrected for bias_term: sum(x - mu) / sqrt(n (1 - mu^2)) "
                                   "(the bias-corrected z).",
    "via_tfs": "signalling_via_tf only: every transcription factor this regulator reaches.",
    "driving_tfs": "signalling_via_tf only: those TFs with significant targets here, as TF(count), most first. "
                   "Regulators sharing driving TFs are largely the same evidence.",
    "activity_t": "t of the slope in a linear model of the rank metric on the regulator's target signs over all "
                  "tested genes: positive = targets move as activation predicts.",
    "activity_p": "Two-sided p of activity_t.",
    "activity_q": "BH q of activity_p within (comparison, network).",
    "regulator_in_data": "yes if the regulator's own gene was tested in the comparison.",
    "regulator_direction": "Direction of the regulator's own gene (up / down), if tested.",
    "regulator_p": "Raw p of the regulator's own gene, if tested.",
    "consistent_targets": "Symbols of the consistent significant targets, ';'-separated.",
    "inconsistent_targets": "Symbols of the inconsistent significant targets, ';'-separated.",
}

FILES: Dict[str, str] = {
    "comparisons.csv": "One row per comparison: design, how many genes were tested and significant, where the "
                       "significant genes went, and how many gene sets reached q < 0.05.",
    "findings.csv": "Index: every ORA, GSEA and regulator result with q < 0.05, ranked within (comparison, method, "
                    "direction, source). A convenience view; the full results are in enrichment_ora.csv, "
                    "enrichment_gsea.csv and regulators.csv.",
    "enrichment_ora.csv": "Every over-representation result: each gene set in the size range, each direction.",
    "enrichment_gsea.csv": "Every GSEA result (comparisons with a direction): each gene set in the size range.",
    "regulators.csv": "Every upstream regulator tested per comparison: overlap of its targets with significant "
                      "genes, activation z and activity t (comparisons with a direction: contrasts, or F-tests "
                      "directed by their contrasts), its targets that agree and disagree.",
    "human_gene_statistics.csv": "One row per (comparison, human gene reached): the representative input gene, "
                                 "how it was chosen, and its statistics. Untested genes are included.",
    "ortholog_crosswalk.csv": "Every input gene and its human ortholog(s) with the evidence; unmapped genes "
                              "included with their status.",
    "collision_ballots.csv": "Every many-to-one collision per comparison, one row per candidate.",
    "gene_sets.csv": "Every gene set of the sources used: identifier, name, category, size.",
}

GLOSSARY: Dict[str, str] = {
    "ORA": "Over-representation analysis: are significant genes more frequent in a gene set than chance allows, "
           "given the genes tested? Exact hypergeometric test.",
    "GSEA": "Gene-set enrichment analysis (preranked): do a set's genes sit towards the top or bottom of the "
            "whole ranked gene list, significant or not?",
    "NES": "Normalised enrichment score; comparable across sets of different sizes. Sign gives direction.",
    "leading edge": "The genes of a set that drive its GSEA score (before the running-sum peak).",
    "q-value": "Benjamini-Hochberg false discovery rate estimate within the stated family of tests.",
    "universe": "The background for ORA: human genes tested in that comparison, not the whole genome.",
    "collision": "Several input genes mapping to one human gene; one is chosen per comparison by vote.",
    "Compara": "Ensembl's gene-tree orthology (the release is in provenance: ensembl_release), the primary "
               "ortholog source.",
    "contrast": "A two-group comparison with a direction (treatment vs control).",
    "burden test": "A gene-burden comparison (kind 'burden'): how many of the treatment cohort (e.g. probands) "
                   "carry qualifying variants in a gene, against the number expected from the control reference "
                   "(e.g. gnomAD allele frequencies). Its ORA splits genes like a contrast (up = more carriers "
                   "than expected, down = fewer; the informative side is up), but no direction-of-effect statistic "
                   "is reported for it -- it measures damage counts, not expression direction -- and its ORA null "
                   "is weighted by each gene's detectability when the input gives Null_selection_probability.",
    "F-test": "An ANOVA test. It has a direction only when its own contrasts in the table give it one "
              "(direction_from); then it is analysed like a contrast, otherwise only ORA of all significant "
              "genes applies.",
}

READING_GUIDE: List[str] = [
    "The files are evidence, not conclusions: nothing is filtered by significance except findings.csv, "
    "which is an index of the q < 0.05 rows of the full tables.",
    "Direction: 'up' means higher in the treatment group than the control group of that comparison "
    "(GSEA: positive NES). An F-test has a direction only when its own contrasts give it one "
    "(comparisons.csv direction_from); otherwise it has none. For a burden test (kind 'burden'), 'up' means "
    "more carriers in the cohort than the reference predicts (log2 observed/expected > 0); the up results are the "
    "ones that answer a burden question, down results (fewer carriers than expected) are reported as computed. "
    "A burden test has no activation z, activity t, pathway activation z or regulator direction (blank, "
    "direction_statistics says why) and no GSEA; its figures are one panel ranked by p.",
    "Burden tests and detectability: under the null, a gene whose qualifying alleles are common in the reference "
    "reaches p < 0.05 far more easily than one whose alleles are all very rare, so plain ORA favours gene sets "
    "rich in such genes. With Null_selection_probability in the input the ORA p is the Wallenius noncentral "
    "hypergeometric tail weighted by it (null_model says so on every row); rows marked 'unadjusted for "
    "detectability' are not, and should be read with that bias in mind.",
    "ORA and GSEA answer different questions: ORA uses only genes past the raw-p threshold; GSEA uses the whole "
    "ranking. A set can appear in one and not the other.",
    "Gene sets overlap heavily, especially GO terms (a gene annotated to a term counts for every broader parent "
    "term). Many significant sets can reflect the same few genes; compare overlap_genes / leading_edge_genes.",
    "q-values are computed within one source (and direction) at a time; they are not comparable across sources "
    "of very different size (50 Hallmark sets vs ~14,000 GO BP terms).",
    "GSEA p-values are permutation estimates; 'permutations' gives their resolution (smallest possible p is "
    "about 2 / permutations, because only nulls of the same sign count).",
    "Human-gene statistics come from the input species' genes through orthology. homology_source ranks the "
    "evidence: protein/ncrna (Compara) > ncbi > symbol. NCBI and symbol pairs carry no identity scores.",
    "When several input genes map to one human gene, one represents it per comparison (collision_ballots.csv "
    "shows the vote); others' statistics are not in human_gene_statistics.csv for that gene.",
    "Input genes with no human ortholog (many lncRNA and unnamed gene models) cannot contribute to any gene "
    "set; comparisons.csv counts them per comparison.",
    "Fold_change follows the input's convention (Paralome bulk: signed 'times', -2 = halved; Paralome "
    "single-cell: log2, -1 = halved; or a plain ratio, 0.5 = halved). Log2_fold_change is log2 as declared "
    "(a comparison has one or the other).",
    "Regulators: a regulator's own expression need not change for its targets to; activation_z and activity_t "
    "describe whether its KNOWN targets moved as its activation would predict. signalling_via_tf regulators are "
    "inferred through a transcription factor and are weaker, indirect evidence.",
    "Many signalling_via_tf regulators act through the same transcription factors: compare driving_tfs "
    "before counting them as independent findings, and read them against the collectri rows of those TFs.",
    "activation_z assumes significant genes are as likely to go up as down; where they are not, compare "
    "activation_z_bias_corrected (bias_term shows how large the correction is).",
    "Regulator networks come from OmniPath filtered to records whose sources allow commercial use; a regulator "
    "absent from the files may simply be absent from that knowledge.",
]


# ---------------------------------------------------------------------------
# reading what crosswalk and enrich wrote
# ---------------------------------------------------------------------------

def _rows(path: Path) -> Iterator[Dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def _flat(d: Optional[Dict[str, object]]) -> str:
    return ";".join(f"{k}={v}" for k, v in (d or {}).items())


class Sources:
    """Locations of everything the report reads, with clear errors for what is missing."""

    def __init__(self, crosswalk_dir: Path):
        # Resolved: the bundle is named after this folder, and `report .` would otherwise name it '_combined_results'.
        self.dir = Path(crosswalk_dir).resolve()
        self.summary_path = self.dir / "comparison_summary.json"
        self.crosswalk_csv = self.dir / "ortholog_crosswalk.csv"
        self.crosswalk_json = self.dir / "ortholog_crosswalk.json"
        self.enrichment_dir = self.dir / "enrichment"
        self.enrichment_json = self.enrichment_dir / "enrichment.json"
        for needed in (self.summary_path, self.crosswalk_csv):
            if not needed.exists():
                raise ProjectError(21, f"{needed} not found")
        self.summary = json.loads(self.summary_path.read_text())
        self.crosswalk = json.loads(self.crosswalk_json.read_text()) if self.crosswalk_json.exists() else {}
        self.enrichment = json.loads(self.enrichment_json.read_text()) if self.enrichment_json.exists() else None
        self.regulators_dir = self.dir / "regulators"
        path = self.regulators_dir / "regulators.json"
        self.regulators = json.loads(path.read_text()) if path.exists() else None
        self.comparisons: List[Dict[str, object]] = self.summary["comparisons"]


def comparisons_rows(src: Sources) -> Tuple[List[str], List[Dict[str, object]]]:
    enrich_by = {c["comparison"]: c for c in (src.enrichment or {}).get("comparisons", [])}
    columns = ["comparison", "kind", "status", "control", "treatment", "direction_from", "alpha", "input_genes_tested",
               "input_genes_significant", "significant_mapped", "significant_via_protein", "significant_via_ncrna",
               "significant_via_ncbi", "significant_via_symbol", "significant_mapped_but_outvoted",
               "significant_excluded_by_compara_precedence", "significant_unmapped",
               "significant_unmapped_by_biotype", "significant_unmapped_by_status", "human_genes_reached",
               "human_genes_tested", "human_genes_significant", "human_genes_with_tested_collisions",
               "ora_universe", "ora_selected", "gsea_rank_metric", "gsea_permutations",
               "ora_sets_q_below_0_05", "gsea_sets_q_below_0_05", "regulators_q_below_0_05", "ora_null",
               "direction_statistics"]
    regs_by = {c["comparison"]: c for c in (src.regulators or {}).get("comparisons", [])}
    rows = []
    for c in src.comparisons:
        e = enrich_by.get(c["comparison"], {})
        rg = regs_by.get(c["comparison"], {})
        via = c.get("significant_mapped_by_source", {}) or {}
        rows.append({
            "comparison": c["comparison"], "kind": c.get("kind"), "status": "analysed",
            "control": c.get("control") or "",
            "treatment": c.get("treatment") or "", "direction_from": c.get("direction_from") or "",
            "alpha": c.get("alpha"),
            "input_genes_tested": c.get("input_genes_tested"), "input_genes_significant": c.get("input_genes_significant"),
            "significant_mapped": c.get("significant_mapped"),
            "significant_via_protein": via.get("protein", 0), "significant_via_ncrna": via.get("ncrna", 0),
            "significant_via_ncbi": via.get("ncbi", 0), "significant_via_symbol": via.get("symbol", 0),
            "significant_mapped_but_outvoted": c.get("significant_mapped_but_outvoted"),
            "significant_excluded_by_compara_precedence": c.get("significant_excluded_by_compara_precedence", ""),
            "significant_unmapped": c.get("significant_unmapped"),
            "significant_unmapped_by_biotype": _flat(c.get("significant_unmapped_by_biotype")),
            "significant_unmapped_by_status": _flat(c.get("significant_unmapped_by_status")),
            "human_genes_reached": c.get("human_genes_reached"), "human_genes_tested": c.get("human_genes_tested"),
            "human_genes_significant": c.get("human_genes_significant"),
            "human_genes_with_tested_collisions": c.get("human_genes_with_tested_collisions"),
            "ora_universe": e.get("universe", ""), "ora_selected": e.get("selected", ""),
            "gsea_rank_metric": e.get("rank_metric", ""), "gsea_permutations": e.get("gsea_permutations", ""),
            "ora_sets_q_below_0_05": _flat(e.get("ora_sets_q_below_0.05")),
            "gsea_sets_q_below_0_05": _flat(e.get("gsea_sets_q_below_0.05")),
            "regulators_q_below_0_05": _flat(rg.get("q_below_0.05")),
            "ora_null": e.get("ora_null", ""),
            "direction_statistics": e.get("direction_statistics") or rg.get("direction_statistics") or "",
        })
    for s in src.summary.get("skipped_comparisons", []):
        rows.append({"comparison": s["comparison"], "kind": s.get("kind"),
                     "status": f"not analysed: the same test as {s['same_test_as']}"})
    return columns, rows


def _stat_columns(kinds: Iterable[str]) -> List[str]:
    kinds = set(kinds)
    return ([c for c in CONTRAST_STATS] if kinds & {"contrast", "burden"} else []) + \
           ([c for c in FTEST_STATS] if "ftest" in kinds else [])


def _long(src: Sources, folder: str, suffix: str, leading: Sequence[str]) -> Tuple[List[str], Iterator[Dict[str, object]]]:
    """Concatenate one per-comparison CSV family, with a comparison column and the union of statistic columns."""
    paths = [src.dir / folder / f"{c['comparison']}{suffix}" for c in src.comparisons]
    present = set()
    for path in paths:
        if path.exists():
            with path.open(newline="", encoding="utf-8") as handle:
                present.update(next(csv.reader(handle), []))
    # the union of statistic columns, less the alternatives no comparison used (e.g. FPKM vs other-unit means)
    stats = [c for c in _stat_columns(c.get("kind", "contrast") for c in src.comparisons)
             if c in present or c not in OPTIONAL_STATS]
    columns = ["comparison"] + list(leading) + stats

    def rows():
        for c in src.comparisons:
            path = src.dir / folder / f"{c['comparison']}{suffix}"
            if path.exists():
                for r in _rows(path):
                    yield {"comparison": c["comparison"], **r}
    return columns, rows()


def _enrichment(src: Sources, name: str) -> Tuple[List[str], Iterator[Dict[str, str]]]:
    from .enrich import GSEA_COLUMNS, ORA_COLUMNS

    columns = list(ORA_COLUMNS if name == "ora" else GSEA_COLUMNS)

    def rows():
        if src.enrichment is None:
            return
        # only the files enrichment.json lists, so an earlier run's leftover gsea.csv is never picked up
        listed = {e.get("comparison"): e for e in src.enrichment.get("comparisons", [])}
        for c in src.comparisons:
            entry = listed.get(c["comparison"])
            if entry is not None and not entry.get(name):
                continue
            path = src.enrichment_dir / c["comparison"] / f"{name}.csv"
            if path.exists():
                yield from _rows(path)
    return columns, rows()


FINDINGS_COLUMNS = ["comparison", "method", "direction", "source", "rank", "set_id", "name", "category",
                    "set_size", "effect", "p", "q", "genes"]


def findings_rows(src: Sources, q_max: float = 0.05) -> List[Dict[str, object]]:
    found: Dict[Tuple[str, str, str, str], List[Dict[str, object]]] = {}
    for method in ("ora", "gsea"):
        for r in _enrichment(src, method)[1]:
            if float(r["q"]) >= q_max:
                continue
            if method == "ora":
                direction, effect, genes = r["direction"], r["fold_enrichment"], r["overlap_genes"]
            else:
                direction = "up" if float(r["nes"]) >= 0 else "down"
                effect, genes = r["nes"], r["leading_edge_genes"]
            key = (r["comparison"], method.upper(), direction, r["source"])
            found.setdefault(key, []).append({
                "comparison": r["comparison"], "method": method.upper(), "direction": direction,
                "source": r["source"], "set_id": r["set_id"], "name": r["name"], "category": r["category"],
                "set_size": r["set_size"], "effect": effect, "p": r["p"], "q": r["q"], "genes": genes})
    for r in _regulator_rows(src)[1]:
        base = {"comparison": r["comparison"], "source": r["network"], "set_id": f"{r['network']}:{r['regulator']}",
                "name": r["regulator"], "category": r["regulator_kind"], "set_size": r["targets_in_universe"]}
        if r.get("activity_q") not in (None, "") and float(r["activity_q"]) < q_max:
            direction = "up" if float(r["activity_t"]) >= 0 else "down"
            # the targets that moved the way the activity says (for an inhibited regulator, the targets moving
            # against its sign)
            genes = r["consistent_targets"] if direction == "up" else r["inconsistent_targets"]
            found.setdefault((r["comparison"], "REGULATOR_ACTIVITY", direction, r["network"]), []).append(
                {**base, "method": "REGULATOR_ACTIVITY", "direction": direction, "effect": r["activity_t"],
                 "p": r["activity_p"], "q": r["activity_q"], "genes": genes})
        if float(r["overlap_q"]) < q_max:
            genes = ";".join(x for x in (r.get("consistent_targets", ""), r.get("inconsistent_targets", "")) if x)
            found.setdefault((r["comparison"], "REGULATOR_OVERLAP", "all", r["network"]), []).append(
                {**base, "method": "REGULATOR_OVERLAP", "direction": "all", "effect": r.get("activation_z", ""),
                 "p": r["overlap_p"], "q": r["overlap_q"], "genes": genes})
    order = {c["comparison"]: i for i, c in enumerate(src.comparisons)}
    rows = []
    for key in sorted(found, key=lambda k: (order.get(k[0], 0), k[1], k[2], k[3])):
        group = sorted(found[key], key=lambda r: (float(r["q"]), float(r["p"]), r["set_id"]))
        for i, r in enumerate(group, 1):
            r["rank"] = i
            rows.append(r)
    return rows


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------

def _write_csv(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> int:
    n = 0
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for r in rows:
            writer.writerow(r)
            n += 1
    tmp.replace(path)
    return n


#: The evidence browser, a single web page shipped in the package and copied into every
#: evidence zip. Opened in a web browser, it is pointed at the zip it came in and reads it in place, on the
#: user's computer: tables, filters, and a pop-up per gene set, regulator and gene, each exporting as TSV.
#: In the zip it is named for the program and the version that wrote it: <Name>_Viewer_v<version>.html.
VIEWER = "evidence_browser.html"


def viewer_name() -> str:
    from .identity import DISPLAY_NAME
    from .version import __version__

    return f"{DISPLAY_NAME}_Viewer_v{__version__}.html"
VIEWER_TEXT = ("A web page for people, not evidence: open it in a web browser and choose this zip; it reads the "
               "CSV files on the user's computer (nothing is uploaded) and lists them with filters, and a pop-up "
               "per gene set, regulator and gene. It adds nothing the CSV files do not hold.")


def viewer_page() -> bytes:
    """The page's bytes, from the installed package, stamped with the program and version writing it."""
    from importlib.resources import files

    from .identity import DISPLAY_NAME
    from .version import __version__

    page = (files(__package__) / "web" / VIEWER).read_bytes()
    blank = b'<meta name="written-by" content="">'
    if page.count(blank) != 1:
        raise ProjectError(24, f"the packaged {VIEWER} has no written-by stamp to fill")
    return page.replace(blank, f'<meta name="written-by" content="{DISPLAY_NAME} {__version__}">'.encode())


def build_report(crosswalk_dir: Path, *, out: Optional[Path] = None, workbook: Optional[Path] = None,
                 log=print) -> Dict[str, Path]:
    from .crosswalk import CROSSWALK_COLUMNS
    from .human_tables import BALLOT_COLUMNS, HUMAN_TABLE_IDENTITY
    from .identity import DISPLAY_NAME, NAME
    from .jsonio import write_json
    from .version import __version__

    src = Sources(crosswalk_dir)
    out = Path(out) if out else src.dir / "evidence"
    out.mkdir(parents=True, exist_ok=True)
    counts: Dict[str, int] = {}
    columns_of: Dict[str, List[str]] = {}

    def emit(name: str, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
        columns_of[name] = list(columns)
        counts[name] = _write_csv(out / name, columns, rows)
        log(f"  {name:<28} {counts[name]:>9} rows")

    emit("comparisons.csv", *comparisons_rows(src))
    emit("findings.csv", FINDINGS_COLUMNS, findings_rows(src))
    emit("enrichment_ora.csv", *_enrichment(src, "ora"))
    emit("enrichment_gsea.csv", *_enrichment(src, "gsea"))
    emit("regulators.csv", *_regulator_rows(src))
    emit("human_gene_statistics.csv", *_long(src, "comparisons", ".human_genes.csv", HUMAN_TABLE_IDENTITY))
    emit("ortholog_crosswalk.csv", CROSSWALK_COLUMNS, _rows(src.crosswalk_csv))
    emit("collision_ballots.csv", *_long(src, "collisions", ".ballots.csv", BALLOT_COLUMNS))
    gene_set_columns = ["set_id", "source", "name", "category", "n_genes", "source_ids"]
    emit("gene_sets.csv", gene_set_columns, _gene_set_rows(src))

    undescribed = sorted({c for cols in columns_of.values() for c in cols} - set(COLUMN_TEXT))
    if undescribed:  # a column an agent would have to guess: refuse, as the input reader does
        raise ProjectError(24, f"columns without a description: {', '.join(undescribed)}")

    enrichment = src.enrichment or {}
    manifest = {
        "schema": f"{NAME}.evidence.v1",
        "about": f"{DISPLAY_NAME} evidence bundle: functional evidence for the comparisons of one input table, "
                 "for an interpreting agent. The program makes no interpretive claims; see reading_guide.",
        "provenance": {
            "input_table": src.summary.get("input_table") or src.crosswalk.get("input_table"),
            "input_folder": src.summary.get("input_folder") or src.crosswalk.get("input_folder"),
            "input_produced_by": src.crosswalk.get("input_produced_by"),
            "genome": src.crosswalk.get("genome"),
            "species": src.crosswalk.get("species") or src.summary.get("species"),
            "orthology_reference": src.crosswalk.get("reference"),
            "input_genes": src.crosswalk.get("input_genes"),
            "input_genes_by_status": src.crosswalk.get("input_genes_by_status"),
            "mapped_input_genes_by_source": src.crosswalk.get("mapped_input_genes_by_source"),
            "crosswalk_written_by": src.crosswalk.get("program_version"),
            "enrichment_written_by": enrichment.get("program_version"),
            "ensembl_release": enrichment.get("ensembl_release"),
            "gene_sets": enrichment.get("gene_sets"),
            "enrichment_settings": enrichment.get("settings"),
            "significance_threshold_alpha": src.summary.get("alpha"),
            "regulators_written_by": (src.regulators or {}).get("program_version"),
            "regulator_data": (src.regulators or {}).get("regulator_data"),
            "regulator_settings": (src.regulators or {}).get("settings"),
        },
        "files": [{"path": name, "rows": counts[name], "description": FILES[name],
                   "columns": [{"name": c, "description": COLUMN_TEXT[c]} for c in columns_of[name]]}
                  for name in FILES],
        "glossary": GLOSSARY,
        "reading_guide": READING_GUIDE,
    }
    # The page in the zip, named here so an agent reading the manifest knows what it is.
    manifest["viewer"] = {"path": viewer_name(), "description": VIEWER_TEXT}
    notes = []
    if src.enrichment is None:
        notes.append("no enrichment in this crosswalk folder: run 'enrich' first for ORA/GSEA tables")
    if src.regulators is None:
        notes.append("no regulator analysis in this crosswalk folder: run 'regulators' first for regulators.csv")
    if src.enrichment is not None and counts.get("gene_sets.csv", 0) == 0:
        notes.append(f"gene_sets.csv is empty: the gene-set index recorded by enrichment "
                     f"({src.enrichment.get('reference_dir')}) was not found from here")
    if notes:
        manifest["notes"] = notes
    manifest_path = write_json(out / "manifest.json", manifest)

    # Named after the crosswalk folder, and so is the folder inside it: two studies' bundles
    # can sit, and be unzipped, side by side.
    stem = f"{src.dir.name}_combined_results"
    archive = src.dir / f"{stem}.zip"
    tmp = archive.with_name(archive.name + ".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for name in ["manifest.json"] + list(FILES):
            # fixed timestamps: the same results make the same bytes
            info = zipfile.ZipInfo(f"{stem}/{name}", date_time=FIXED_ZIP_TIME)
            info.compress_type, info.external_attr = zipfile.ZIP_DEFLATED, 0o644 << 16
            with (out / name).open("rb") as src_file, zf.open(info, "w", force_zip64=True) as dest:
                shutil.copyfileobj(src_file, dest, 1 << 20)
        info = zipfile.ZipInfo(f"{stem}/{viewer_name()}", date_time=FIXED_ZIP_TIME)
        info.compress_type, info.external_attr = zipfile.ZIP_DEFLATED, 0o644 << 16
        zf.writestr(info, viewer_page())
    tmp.replace(archive)
    log(f"  {archive.name:<28} {archive.stat().st_size / 1e6:>8.1f} MB")

    book = Path(workbook) if workbook else src.dir / f"{src.dir.name}_workbook.xlsx"
    _workbook(book, out, manifest, columns_of)
    log(f"  {book.name:<28} {book.stat().st_size / 1e6:>8.1f} MB")
    return {"evidence": out, "manifest": manifest_path, "archive": archive, "workbook": book}


def _regulator_rows(src: Sources) -> Tuple[List[str], Iterator[Dict[str, str]]]:
    from .regulators import REGULATOR_COLUMNS

    def rows():
        if src.regulators is None:
            return
        for c in src.comparisons:
            path = src.regulators_dir / c["comparison"] / "regulators.csv"
            if path.exists():
                yield from _rows(path)
    return list(REGULATOR_COLUMNS), rows()


def _gene_set_rows(src: Sources) -> Iterator[Dict[str, str]]:
    if src.enrichment is None:
        return
    used = set((src.enrichment.get("gene_sets") or {}).get("sources_used") or [])
    ref, release = src.enrichment.get("reference_dir"), src.enrichment.get("ensembl_release")
    if not ref:
        return
    from .genesets import genesets_dir

    path = genesets_dir(Path(ref), int(release)) / "gene_sets.tsv"
    if not path.exists():
        return
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source"] in used:
                yield row


# ---------------------------------------------------------------------------
# the workbook: a view for people
# ---------------------------------------------------------------------------

WORKBOOK_P_MAX = 0.05
TEXT_COLUMNS = ("via_tfs", "driving_tfs", "regulator", "network", "consistent_targets", "inconsistent_targets", "comparison", "set_id", "name", "category", "human_gene_id", "human_symbol", "input_gene_id",
                "input_symbol", "source_gene_id", "source_symbol", "candidate_gene_id", "candidate_symbol",
                "candidates", "overlap_genes", "leading_edge_genes", "genes", "control", "treatment")


def _workbook(path: Path, evidence: Path, manifest: Dict[str, object], columns_of: Dict[str, List[str]]) -> Path:
    from .xlsx import Sheet, write_workbook

    def sheet(title: str, name: str, keep=None) -> Sheet:
        columns = columns_of[name]

        def rows():
            for r in _rows(evidence / name):
                if keep is None or keep(r):
                    yield [r.get(c, "") for c in columns]
        return Sheet(title, columns, rows(), text_columns=TEXT_COLUMNS)

    def p_ok(r):
        try:
            return float(r["p"]) < WORKBOOK_P_MAX
        except (KeyError, ValueError):
            return False

    def reg_ok(r):
        values = [r.get("overlap_p", ""), r.get("activity_p", "")]
        return any(v not in ("", None) and float(v) < WORKBOOK_P_MAX for v in values)

    readme = [[f"{manifest['about']}"], [""],
              ["This workbook is a view of the evidence folder beside it (the CSV files are the record)."],
              [f"The ORA, GSEA and Regulators sheets show only results with p < {WORKBOOK_P_MAX} (Regulators: "
               "overlap or activity); every result is in evidence/enrichment_ora.csv, enrichment_gsea.csv and "
               "regulators.csv."], [""],
              ["Provenance"]]
    for key, value in (manifest["provenance"] or {}).items():  # type: ignore[union-attr]
        if isinstance(value, (dict, list)):
            value = json.dumps(value, sort_keys=True)[:2000]
        readme.append([f"{key}: {value}"])
    readme += [[""], ["Reading guide"]] + [[f"- {line}"] for line in manifest["reading_guide"]]  # type: ignore
    readme += [[""], ["Glossary"]] + [[f"{k}: {v}"] for k, v in manifest["glossary"].items()]  # type: ignore
    readme += [[""], ["Columns"]]
    for f in manifest["files"]:  # type: ignore[union-attr]
        readme.append([f"{f['path']} -- {f['description']}"])
        readme += [[f"    {c['name']}: {c['description']}"] for c in f["columns"]]
    # Excel holds 1,048,575 data rows per sheet. The README is written first, so a sheet that may be cut is named
    # here from the row counts already known, so the cut is counted and reported.
    from .xlsx import MAX_ROWS

    # Only sheets that hold every row of their CSV can be judged from its count; the p < 0.05 sheets are judged
    # after writing (the stderr report below), since their CSVs' counts overstate them.
    unfiltered = {"comparisons.csv", "findings.csv", "human_gene_statistics.csv", "ortholog_crosswalk.csv",
                  "collision_ballots.csv"}
    too_long = [(f["path"], f["rows"]) for f in manifest["files"]  # type: ignore
                if f["path"] in unfiltered and int(f["rows"]) > MAX_ROWS]
    if too_long:
        readme[2:2] = [[f"NOTE: {name} has {rows:,} rows; its sheet holds at most the first {MAX_ROWS:,}. "
                        f"The full table is evidence/{name}."] for name, rows in too_long]
    sheets = [
        Sheet("README", ["README"], readme, header=False, widths=[160], wrap=True),
        sheet("Comparisons", "comparisons.csv"),
        sheet("Findings q<0.05", "findings.csv"),
        sheet("ORA p<0.05", "enrichment_ora.csv", p_ok),
        sheet("GSEA p<0.05", "enrichment_gsea.csv", p_ok),
        sheet("Regulators p<0.05", "regulators.csv", reg_ok),
        sheet("Human genes", "human_gene_statistics.csv"),
        sheet("Crosswalk", "ortholog_crosswalk.csv"),
        sheet("Collisions", "collision_ballots.csv"),
    ]
    written = write_workbook(path, sheets)
    for s in sheets:
        if s.cut:
            print(f"  workbook sheet {s.name!r}: {s.cut:,} rows beyond Excel's limit left out "
                  f"(all of them are in the evidence CSV)", file=sys.stderr)
    return written
