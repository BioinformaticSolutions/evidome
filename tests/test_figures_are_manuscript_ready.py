"""The mirrored top-25 figure: what it selects, what it draws, and how it is written."""

import csv
import json
import math

import pytest

from evidome import figures
from evidome.cli import main
from evidome.enrich import GSEA_COLUMNS, ORA_COLUMNS
from evidome.errors import ProjectError
from evidome.figures import FigureSettings, Item, make_figures, parse_colors, parse_formats, split_top
from evidome.regulators import REGULATOR_COLUMNS


def _write(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


@pytest.fixture
def run(tmp_path):
    d = tmp_path / "study"
    d.mkdir()
    (d / "comparison_summary.json").write_text(json.dumps({"comparisons": [{"comparison": "A"}, {"comparison": "F1"}]}))
    regs = []
    for i in range(30):
        regs.append({"network": "collectri", "regulator": f"UP{i}", "overlap_p": 10 ** -(40 - i), "activation_z": 2.5 + i / 10})
        regs.append({"network": "collectri", "regulator": f"DN{i}", "overlap_p": 10 ** -(30 - i / 2), "activation_z": -3})
    regs += [{"network": "collectri", "regulator": "WEAKZ", "overlap_p": 1e-50, "activation_z": 1.9},
             {"network": "collectri", "regulator": "NOTSIG", "overlap_p": 0.2, "activation_z": 5},
             {"network": "collectri", "regulator": "ZERO_P", "overlap_p": 0.0, "activation_z": 4},
             {"network": "signalling_via_tf", "regulator": "TNF", "overlap_p": 1e-4, "activation_z": 2.2},
             {"network": "signalling_via_tf", "regulator": "NA", "overlap_p": 1e-4, "activation_z": ""}]
    _write(d / "regulators/A/regulators.csv", REGULATOR_COLUMNS, regs)
    ora = []
    for sid, up, down, p in (("S1", 9, 0, 1e-8), ("S2", 1, 8, 1e-6), ("S3", 3, 3, 1e-9), ("S4", 0, 0, 1)):
        for direction, overlap in (("up", up), ("down", down), ("all", up + down)):
            ora.append({"direction": direction, "source": "reactome", "set_id": sid, "name": f"Pathway {sid}",
                        "overlap": overlap, "p": p if direction == "all" else 0.5})
    _write(d / "enrichment/A/ora.csv", ORA_COLUMNS, ora)
    _write(d / "enrichment/A/gsea.csv", GSEA_COLUMNS,
           [{"source": "hallmark", "set_id": "H1", "name": "HALLMARK_ONE", "nes": 1.7, "p": 0.001},
            {"source": "hallmark", "set_id": "H2", "name": "HALLMARK_TWO", "nes": -1.4, "p": 0.01}])
    # F1 is an F-test: ORA has only direction "all" -- no pathway figure for it
    _write(d / "enrichment/F1/ora.csv", ORA_COLUMNS,
           [{"direction": "all", "source": "reactome", "set_id": "S1", "name": "x", "overlap": 5, "p": 1e-5}])
    return d


def test_selection_is_the_lowest_p_on_each_side_with_no_cutoff():
    items = [Item("a", "a", 1e-5, 2.0), Item("b", "b", 1e-9, -2.5), Item("c", "c", 0.3, 0.4),
             Item("d", "d", 0.9, 6.0), Item("e", "e", 1e-7, 3.0), Item("z", "z", 1e-12, 0.0),
             Item("n", "n", 0.8, -0.1)]
    neg, pos = split_top(items, top=3)
    assert [i.label for i in neg] == ["z", "b", "n"]     # z: score 0, no sign, lowest p: on both sides for now
    assert [i.label for i in pos] == ["z", "e", "a"]     # a is not significant and is still drawn


def test_no_direction_places_by_p_and_is_drawn_grey(run):
    items = [Item("u1", "u1", 1e-9, None), Item("u2", "u2", 1e-8, None), Item("r", "r", 1e-7, 3.0),
             Item("l", "l", 1e-6, -3.0), Item("far", "far", 0.5, None)]
    neg, pos = split_top(items, top=2)
    assert [i.label for i in neg] == ["u1", "u2"] and [i.label for i in pos] == ["u1", "u2"]
    make_figures(run, FigureSettings(formats=("svg",)), threads=1, log=lambda *_: None)
    svg = (run / "figures/A_pathways.svg").read_text()
    assert "not calculable" in svg and ">NC<" in svg and "#b4b4b4" in svg.lower()


def test_zero_is_white_and_not_calculable_is_grey(tmp_path):
    from evidome.figures import FigureSettings, Panel, draw

    panel = Panel("t", "pathways", "t", "x", "z", "Decreased", "Increased",
                  [Item("zero", "zero", 1e-5, 0.0), Item("nc", "nc", 1e-4, None)],
                  [Item("up", "up", 1e-5, 2.0), Item("up2", "up2", 1e-4, 1.0)])
    draw(panel, tmp_path, FigureSettings(formats=("svg",)), version="t")
    svg = (tmp_path / "t.svg").read_text().lower()
    assert "fill: #ffffff" in svg and "#b4b4b4" in svg and "0.00" in svg
    assert list(figures.csv.DictReader((tmp_path / "t.csv").open()))[1]["score"] == ""


def test_one_regulator_figure_across_networks(run):
    [reg] = [p for p in figures.panels_for(run, "A", 25) if p.kind == "regulators"]
    assert reg.stem == "A_regulators"
    assert len(reg.right) == 25 and len(reg.left) == 25
    assert reg.right[0].label == "ZERO_P" and reg.right[0].neg_log10_p == figures.LOG_P_CAP
    labels = [i.label for i in reg.left + reg.right]
    assert len(labels) == len(set(labels))


def test_a_regulator_in_both_networks_appears_once_with_its_lower_p(tmp_path):
    path = tmp_path / "r.csv"
    _write(path, REGULATOR_COLUMNS, [
        {"network": "collectri", "regulator": "STAT3", "overlap_p": 1e-4, "activation_z": 3},
        {"network": "signalling_via_tf", "regulator": "STAT3", "overlap_p": 1e-9, "activation_z": 2.5},
        {"network": "collectri", "regulator": "E2F1", "overlap_p": 1e-6, "activation_z": -3}])
    panel = figures.regulator_panel(path, "A", 25)
    assert [(i.label, i.p) for i in panel.right] == [("STAT3", 1e-9)]
    assert panel.right[0].detail.startswith("signalling_via_tf")


def test_one_pathway_figure_across_sources_each_row_naming_its_source(run):
    [path] = [p for p in figures.panels_for(run, "A", 25) if p.kind == "pathways"]
    assert path.stem == "A_pathways"
    assert [i.id for i in path.right] == ["S3", "S1", "S4"]          # S3 (0) and S4 (NC) on both sides
    assert [i.id for i in path.left] == ["S3", "S2", "S4"]
    assert path.right[1].label == "Pathway S1" and path.right[1].tag == "Reactome"
    scores = {i.id: i.score for i in path.left + path.right}
    assert scores["S1"] == 3.0 and round(scores["S2"], 4) == round(-7 / 3, 4)
    assert scores["S3"] == 0 and scores["S4"] is None          # 3 up / 3 down; no significant gene: both grey
    assert not [p for p in figures.panels_for(run, "F1", 25) if p.kind == "pathways"]


def test_two_figures_per_comparison_at_most(run):
    assert sorted(p.stem for p in figures.panels_for(run, "A", 25)) == ["A_pathways", "A_regulators"]


def test_every_format_is_written_at_manuscript_quality(run):
    result = make_figures(run, FigureSettings(formats=("pdf", "svg", "png", "tiff"), dpi=300), threads=1,
                          log=lambda *_: None)
    stem = run / "figures" / "A_regulators"
    for ext in ("pdf", "svg", "png", "tiff", "csv"):
        assert stem.with_suffix(f".{ext}").stat().st_size > 0, ext
    pdf = stem.with_suffix(".pdf").read_bytes()
    assert b"/FontFile2" in pdf                          # TrueType embedded: editable, no substitution
    svg = stem.with_suffix(".svg").read_text()
    assert "<text" in svg and "ESR1" not in svg and "UP0" in svg      # text kept as text
    from PIL import Image
    with Image.open(stem.with_suffix(".png")) as im:
        assert round(im.info["dpi"][0]) == 300 and abs(im.size[0] - 7.0 * 300) <= 2
    rows = list(csv.DictReader(stem.with_suffix(".csv").open()))
    assert rows[0]["label"] == "DN0" and rows[0]["side"] == "negative" and len(rows) == 50
    assert result["figures"] == len(list((run / "figures").rglob("*.pdf")))


def test_colours_are_the_pair_given(run):
    make_figures(run, FigureSettings(colors="#00FF00,#FF00FF", formats=("svg",)), threads=1, log=lambda *_: None)
    svg = (run / "figures/A_regulators.svg").read_text().lower()
    assert "#00ff00" in svg and "#ff00ff" in svg          # bar edges are the pure end colours
    default = FigureSettings()
    assert parse_colors(default.colors) == ("#0044CC", "#CC0000")


@pytest.mark.parametrize("text", ["red", "red,green,blue", "red,notacolour", ""])
def test_a_colour_list_that_is_not_two_colours_is_prv_0026(text):
    with pytest.raises(ProjectError) as err:
        parse_colors(text)
    assert err.value.code == 26


def test_formats_are_checked():
    assert parse_formats("PDF, png") == ("pdf", "png")
    with pytest.raises(ProjectError):
        parse_formats("pdf,jpg")


def test_bad_colours_stop_the_pipeline_before_it_starts(tmp_path, capsys):
    assert main([str(tmp_path / "missing.csv"), "--out", str(tmp_path / "o"), "--colors", "blue"]) == 2
    assert "EVD-0026" in capsys.readouterr().err
    assert not (tmp_path / "o").exists()


def test_without_matplotlib_the_run_skips_and_the_command_says_why(run, monkeypatch, capsys):
    monkeypatch.setattr(figures, "available", lambda: False)
    notes = []
    assert make_figures(run, log=notes.append)["files"] == [] and "skipped" in notes[0]
    assert main(["figures", str(run)]) == 2
    assert "EVD-0027" in capsys.readouterr().err


def test_the_figures_command_redraws(run, capsys):
    assert main(["figures", str(run), "--colors", "navy,darkorange", "--figure-formats", "svg", "--threads", "1"]) == 0
    assert "figures in" in capsys.readouterr().out
    assert (run / "figures/A_pathways.svg").exists()


def test_both_sides_always_have_the_same_number_of_bars(run):
    for panel in figures.panels_for(run, "A", 25):
        assert len(panel.left) == len(panel.right) > 0, panel.stem
    [few] = [p for p in figures.panels_for(run, "A", 3) if p.kind == "regulators"]
    assert [i.label for i in few.right] == ["ZERO_P", "WEAKZ", "UP0"] and len(few.left) == 3   # |z| < 2 is drawn too


def test_a_long_pathway_name_is_cut_but_its_source_never_is():
    text = figures._wrap("Respiratory electron transport, ATP synthesis by chemiosmotic coupling, and heat "
                         "production by uncoupling proteins.", "Reactome")
    lines = text.split("\n")
    assert len(lines) == 2 and lines[-1].endswith("… [Reactome]") and all(len(l) <= 40 for l in lines)
    assert figures._wrap("Adipogenesis", "WikiPathways") == "Adipogenesis [WikiPathways]"


def test_go_terms_are_not_drawn(tmp_path):
    path = tmp_path / "ora.csv"
    rows = []
    for source, sid in (("go_bp", "GO:1"), ("go_cc", "GO:2"), ("reactome", "R1"), ("hallmark", "HALLMARK_HYPOXIA")):
        for direction, overlap in (("up", 9), ("down", 0), ("all", 9)):
            rows.append({"direction": direction, "source": source, "set_id": sid, "name": sid, "overlap": overlap,
                         "p": 1e-9 if direction == "all" else 0.5})
    _write(path, ORA_COLUMNS, rows)
    panel = figures.pathway_panel(path, "A", 25)
    assert sorted((i.label, i.tag) for i in panel.right) == [("HYPOXIA", "Hallmark"), ("R1", "Reactome")]


def test_axes_end_on_round_numbers_and_every_figure_of_a_kind_has_one_shape():
    assert [figures.nice_ceiling(v) for v in (0.7, 1.0, 13, 236, 351)] == pytest.approx([0.7, 1.0, 14, 250, 400])
    assert len({geo["panel"] for geo in figures.LAYOUT.values()}) == 1          # one axis scale for both kinds
    assert figures.LAYOUT["pathways"]["label"] > figures.LAYOUT["regulators"]["label"]


def test_every_figure_has_the_same_height_however_many_bars(run):
    from PIL import Image

    make_figures(run, FigureSettings(formats=("png",), dpi=100), threads=1, log=lambda *_: None)
    sizes = {p.stem: Image.open(p).size for p in (run / "figures").glob("*.png")}
    assert len({h for _, h in sizes.values()}) == 1                             # pathways has 1 bar a side here
    assert sizes["A_pathways"][0] > sizes["A_regulators"][0]



def test_exclude_takes_unsigned_results_off_and_the_next_by_p_fill_in(run):
    from evidome.figures import parse_exclude

    assert parse_exclude("") == () and parse_exclude("unsigned") == ("nc", "zero")
    assert parse_exclude("NC, zero") == ("nc", "zero") and parse_exclude("zero") == ("zero",)
    with pytest.raises(ProjectError):
        parse_exclude("grey")
    ids = lambda panel: ([i.id for i in panel.left], [i.id for i in panel.right])   # noqa: E731
    [both] = [p for p in figures.panels_for(run, "A", 25) if p.kind == "pathways"]
    [no_nc] = [p for p in figures.panels_for(run, "A", 25, ("nc",)) if p.kind == "pathways"]
    [neither] = [p for p in figures.panels_for(run, "A", 25, ("nc", "zero")) if p.kind == "pathways"]
    assert ids(both) == (["S3", "S2", "S4"], ["S3", "S1", "S4"])
    assert ids(no_nc) == (["S3", "S2"], ["S3", "S1"])
    assert ids(neither) == (["S2"], ["S1"])


def test_exclude_travels_to_the_slurm_finish_job_and_the_figures_command(run, capsys):
    assert FigureSettings.from_dict(FigureSettings(exclude=("nc",)).as_dict()).exclude == ("nc",)
    assert main(["figures", str(run), "--exclude", "unsigned", "--figure-formats", "csv"]) == 2   # csv is not a format
    assert main(["figures", str(run), "--exclude", "unsigned", "--figure-formats", "svg", "--threads", "1"]) == 0
    rows = list(csv.DictReader((run / "figures/A_pathways.csv").open()))
    assert {r["id"] for r in rows} == {"S1", "S2"}


def test_an_f_test_draws_no_all_grey_regulator_figure(run):
    # F-tests have no direction, so regulators.csv has no activation z for any regulator
    _write(run / "regulators/F1/regulators.csv", REGULATOR_COLUMNS,
           [{"network": "collectri", "regulator": f"R{i}", "overlap_p": 10 ** -i, "activation_z": ""}
            for i in range(1, 40)])
    assert figures.panels_for(run, "F1", 25) == []
    assert sorted(p.stem for p in figures.panels_for(run, "A", 25)) == ["A_pathways", "A_regulators"]
