"""Pathway activation from Reactome's reaction logic."""

import csv
import gzip
import zipfile

import pytest

from evidome import figures
from evidome.enrich import ORA_COLUMNS, ComparisonTable, pathway_activation
from evidome.errors import ProjectError
from evidome.pathway_signs import build_signs, pathway_patterns, read_biopax, uniprot_to_ensembl

BP = "http://www.biopax.org/release/biopax-level3.owl#"


def _owl(about: bool = False) -> str:
    """A pathway in BioPAX Level 3 as Reactome writes it: a reaction with inputs, a catalysis, and a negative
    regulation whose controller is a complex; a sub-pathway with a reaction; a set (memberPhysicalEntity)."""
    base = "http://www.reactome.org/biopax/Homo_sapiens#" if about else ""
    ident = (lambda i: f'rdf:about="{base}{i}"') if about else (lambda i: f'rdf:ID="{i}"')
    ref = lambda i: f'rdf:resource="{base if about else "#"}{i}"'   # noqa: E731
    prot = lambda pid, acc: (f'<bp:Protein {ident(pid)}><bp:entityReference {ref("PR_" + acc)}/></bp:Protein>'   # noqa: E731
                             f'<bp:ProteinReference {ident("PR_" + acc)}><bp:xref {ref("X_" + acc)}/></bp:ProteinReference>'
                             f'<bp:UnificationXref {ident("X_" + acc)}><bp:db>UniProt</bp:db><bp:id>{acc}</bp:id>'
                             f'</bp:UnificationXref>')
    return f"""<?xml version="1.0"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:bp="{BP}">
<bp:Pathway {ident("P1")}><bp:displayName>Signalling</bp:displayName>
  <bp:pathwayComponent {ref("R1")}/><bp:pathwayComponent {ref("P2")}/><bp:xref {ref("XP1")}/></bp:Pathway>
<bp:UnificationXref {ident("XP1")}><bp:db>Reactome</bp:db><bp:id>R-HSA-100</bp:id></bp:UnificationXref>
<bp:Pathway {ident("P2")}><bp:pathwayComponent {ref("R2")}/><bp:xref {ref("XP2")}/></bp:Pathway>
<bp:UnificationXref {ident("XP2")}><bp:db>Reactome</bp:db><bp:id>R-HSA-200.3</bp:id></bp:UnificationXref>
<bp:BiochemicalReaction {ident("R1")}><bp:left {ref("A")}/><bp:right {ref("B")}/></bp:BiochemicalReaction>
<bp:BiochemicalReaction {ident("R2")}><bp:left {ref("SET")}/></bp:BiochemicalReaction>
<bp:Catalysis {ident("C1")}><bp:controlled {ref("R1")}/><bp:controller {ref("K")}/>
  <bp:controlType>ACTIVATION</bp:controlType></bp:Catalysis>
<bp:Control {ident("N1")}><bp:controlled {ref("R1")}/><bp:controller {ref("CPLX")}/>
  <bp:controlType>INHIBITION</bp:controlType></bp:Control>
<bp:Control {ident("N2")}><bp:controlled {ref("C1")}/><bp:controller {ref("PH")}/>
  <bp:controlType>INHIBITION-COMPETITIVE</bp:controlType></bp:Control>
<bp:Complex {ident("CPLX")}><bp:component {ref("I1")}/><bp:component {ref("I2")}/></bp:Complex>
<bp:Protein {ident("SET")}><bp:memberPhysicalEntity {ref("M1")}/><bp:memberPhysicalEntity {ref("M2")}/></bp:Protein>
<bp:SmallMolecule {ident("ATP")}/>
{prot("A", "P00001")}{prot("B", "P00002-2")}{prot("K", "P00003")}{prot("I1", "P00004")}{prot("I2", "P00005")}
{prot("PH", "P00006")}{prot("M1", "P00007")}{prot("M2", "P00008")}
</rdf:RDF>"""


@pytest.mark.parametrize("about", [False, True])
def test_signs_follow_reactions_complexes_sets_and_negative_regulation(tmp_path, about):
    path = tmp_path / "biopax.zip"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("Homo_sapiens.owl", _owl(about))
        zf.writestr("Mus_musculus.owl", "<x/>")
    patterns, counts = pathway_patterns(read_biopax(path))
    assert patterns["R-HSA-100"] == {"P00001": 1, "P00002": 1, "P00003": 1, "P00004": -1, "P00005": -1,
                                     "P00006": -1, "P00007": 1, "P00008": 1}   # sub-pathway included
    assert patterns["R-HSA-200"] == {"P00007": 1, "P00008": 1}               # version suffix dropped
    assert counts["controls_inhibition"] == 2 and counts["pathways"] == 2


def test_uniprot_to_genes_and_the_index(tmp_path):
    zpath = tmp_path / "biopax.zip"
    with zipfile.ZipFile(zpath, "w") as zf:
        zf.writestr("Homo_sapiens.owl", _owl())
    tsv = tmp_path / "Homo_sapiens.GRCh38.116.uniprot.tsv.gz"
    lines = ["gene_stable_id\ttranscript_stable_id\tprotein_stable_id\txref\tdb_name"]
    for i, acc in enumerate(["P00001", "P00002", "P00003", "P00004", "P00005", "P00006", "P00007"], 1):
        lines.append(f"ENSG{i:011d}\tENST1\tENSP1\t{acc}\tUniprot/SWISSPROT")
    lines.append("ENSG00000000099\tENST9\tENSP9\tP00007\tUniprot/SPTREMBL")    # P00007 on two genes: ambiguous
    with gzip.open(tsv, "wt") as h:
        h.write("\n".join(lines) + "\n")
    to_gene, ambiguous = uniprot_to_ensembl(tsv)
    assert ambiguous == 1 and "P00007" not in to_gene
    human = {f"ENSG{i:011d}" for i in range(1, 8)}
    rows, counts = build_signs(zpath, tsv, human, {"REACTOME:R-HSA-100"})
    signs = {g: s for _, g, s in rows}
    assert signs == {"ENSG00000000001": 1, "ENSG00000000002": 1, "ENSG00000000003": 1, "ENSG00000000004": -1,
                     "ENSG00000000005": -1, "ENSG00000000006": -1}
    assert counts["pathways_not_in_gene_sets"] == 1 and counts["negative_pattern_rows"] == 3


def test_a_file_without_the_expected_layout_is_refused(tmp_path):
    path = tmp_path / "x.owl"
    path.write_text('<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"/>')
    with pytest.raises(ProjectError) as err:
        build_signs(path, path, set(), set())
    assert err.value.code == 18 and "negative regulations 0" in str(err.value)


def _table(direction):
    genes = sorted(direction)
    return ComparisonTable("A", "contrast", genes, {g: g.lower() for g in genes}, {g: 0.01 for g in genes},
                           dict(direction), {g: float(d) for g, d in direction.items()}, direction_from="T-statistic")


def test_activation_z_counts_a_brake_going_down_as_activation():
    pattern = {"G1": 1, "G2": 1, "G3": 1, "BRAKE": -1, "OFF": 1}
    t = _table({"G1": 1, "G2": 1, "G3": -1, "BRAKE": -1, "OFF": 1})
    out = pathway_activation(t, pattern, {"G1", "G2", "G3", "BRAKE"})       # OFF: not significant
    assert (out["pattern_consistent"], out["pattern_inconsistent"]) == (3, 1)
    assert out["pathway_activation_z"] == 1.0 and out["pattern_consistent_genes"] == "brake;g1;g2"
    few = pathway_activation(t, pattern, {"G1", "G2", "BRAKE"})
    assert few["pathway_activation_z"] == "" and few["pattern_genes"] == 5        # N = 3 < 4: no z
    assert pathway_activation(t, None, {"G1"}) == {}


def _ora_rows():
    rows = []
    for sid, source, z in (("REACTOME:R1", "reactome", "2.5"), ("REACTOME:R2", "reactome", ""),
                           ("W1", "wikipathways", None), ("REACTOME:R3", "reactome", "-3.0")):
        for d, ov in (("up", 5), ("down", 1), ("all", 6)):
            row = {"direction": d, "source": source, "set_id": sid, "name": sid, "overlap": ov, "p": 1e-6 if d == "all" else 0.5}
            if d == "all" and z is not None:
                row.update({"pattern_genes": 8, "pathway_activation_z": z, "pattern_consistent": 5,
                            "pattern_inconsistent": 1})
            rows.append(row)
    return rows


def test_the_pathway_figure_uses_the_activation_z_and_greys_what_has_no_pattern(tmp_path):
    path = tmp_path / "ora.csv"
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, ORA_COLUMNS, extrasaction="ignore")
        w.writeheader()
        w.writerows(_ora_rows())
    panel = figures.pathway_panel(path, "A", 25)
    scores = {i.id: i.score for i in panel.left + panel.right}
    assert scores == {"REACTOME:R1": 2.5, "REACTOME:R2": None, "W1": None, "REACTOME:R3": -3.0}
    assert panel.right_title == "Predicted activated" and "Activation" in panel.score_name
    old = figures.pathway_panel(path, "A", 25, score="direction")
    assert {i.id: round(i.score, 3) for i in old.right}["W1"] == round(4 / 6 ** 0.5, 3)


def test_reference_pathway_signs_builds_the_index_and_load_gene_sets_attaches_it(tmp_path):
    from evidome.genesets import load_gene_sets
    from evidome.pathway_signs import fetch_pathway_signs

    ref = tmp_path / "ref"
    rel = ref / "ensembl-116"
    (rel / "genesets").mkdir(parents=True)
    (rel / "human_genes.tsv").write_text("human_gene_id\n" + "".join(f"ENSG{i:011d}\n" for i in range(1, 8)))
    (rel / "genesets/gene_sets.tsv").write_text("set_id\tsource\tname\tcategory\tn_genes\tsource_ids\n"
                                               "REACTOME:R-HSA-100\treactome\tSignalling\t\t3\t3\n")
    (rel / "genesets/members.tsv").write_text("set_id\thuman_gene_id\n" + "".join(
        f"REACTOME:R-HSA-100\tENSG{i:011d}\n" for i in (1, 4, 6)))
    src = tmp_path / "mirror"
    src.mkdir()
    with zipfile.ZipFile(src / "biopax.zip", "w") as zf:
        zf.writestr("Homo_sapiens.owl", _owl())
    with gzip.open(src / "uniprot.tsv.gz", "wt") as h:
        h.write("gene_stable_id\txref\n" + "".join(f"ENSG{i:011d}\tP{i:05d}\n" for i in range(1, 8)))
    manifest = fetch_pathway_signs(ref, urls={"biopax": (src / "biopax.zip").as_uri(),
                                              "uniprot": (src / "uniprot.tsv.gz").as_uri()}, log=lambda *_: None)
    assert manifest.exists()
    info, _ = load_gene_sets(ref)
    assert info["REACTOME:R-HSA-100"]["pattern"]["ENSG00000000004"] == -1
