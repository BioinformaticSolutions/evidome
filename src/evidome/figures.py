"""Figures: the mirrored top-25 bar chart, for regulators and for pathways.

The top 25 activated and inhibited regulators and pathways, drawn to manuscript
quality in design and image resolution. Blue/red by default; ``--colors A,B``
sets the pair (low, high).

One layout for every figure, read from the tables the program already wrote
(nothing is recomputed here):

    [ bars grow left  | names ][ names | bars grow right ]
      negative z, lowest p       positive z, lowest p
                  colour bar: the signed score, centred on white

Bar length is -log10 p, labelled at the bar end; the signed score is printed
inside the bar near the centre, and sets the colour. Rows are the 25 lowest p
on each side -- no significance or score cutoff, so the figure shows the true
distribution, significant or not; a thin dashed line marks p = 0.05 for the
reader. A side is decided by the score's sign. A
score of exactly 0 is a real answer, drawn white (the colour scale's centre);
a result whose score is not calculable (NC) is
drawn grey, its score shown as NC, and still places wherever its p
earns it; having no side of its own, it is drawn on both sides. ``--exclude nc``,
``zero`` or ``unsigned`` takes them off the figure; the next results by p take
their places. A comparison where no result has a calculable score
(an F-test whose contrasts give it no direction) gets no figure: it would be all grey.
Always as many bars on one side as on the other: on the rare input with
fewer than 25 on a side, the shorter side sets the count. A burden test
is the exception: it has no direction-of-effect score (a burden
test measures damage counts, not expression direction -- no activation z is
computed for it), so its figures are ONE panel: the ``--figure-top`` results
ranked by p, bars in the high colour, no "activated / inhibited" side titles,
no score and no colour bar; the title names the test.

What each family plots -- chosen so the score is the one the table already
reports, never a new one presented as the program's finding:

Two figures per comparison, the top hits across every database -- one for
pathways and one for regulators:

regulators  (CollecTRI and signalling_via_tf together; a regulator tested in
            both appears once, with its lower-p row -- the CSV names the network)
    p = overlap_p (hypergeometric: are the regulator's targets over-
    represented among significant genes);
    score = activation_z (the activation z-score). A regulator with no
    signed target among the significant genes has no z: NC, grey.
pathways    (Reactome, WikiPathways and Hallmark together, by over-
            representation; GO is left out as too generic to stand in front of
            pathways. Each row ends with its source, e.g.
            "Interferon alpha/beta signaling [Reactome]"; Hallmark names lose
            their HALLMARK_ prefix and underscores)
    p = ORA p for all significant genes (direction "all");
    score, --pathway-score activation (the default) =
    pathway_activation_z: Reactome's own activity pattern says which
    genes should rise and which fall when the pathway is active, and the z
    counts how many of the set's significant genes agree. Sets without a
    pattern (WikiPathways, Hallmark, Reactome sets with too few signed genes)
    are NC, grey; titles read "Predicted inhibited" / "Predicted activated".
    score, --pathway-score direction = (up - down) / sqrt(overlap), counting
    the set's significant genes by direction -- the z-score of GOplot (Walter
    et al. 2015, Bioinformatics 31:2912): which way the genes moved, not that
    the pathway is "activated"; titles read "Decreased" / "Increased". A set
    with no significant gene has no score (NC, grey); one with as many up as
    down scores 0 (white).
ORA p values are comparable across sources (each is an exact hypergeometric
test on the same gene universe), so one ranking across them is sound. GSEA
results stay in the tables; they are not drawn.

Output, in ``<out>/figures/``, per comparison (<stem> = <comparison>_regulators
or <comparison>_pathways):
    <stem>.pdf   vector, fonts embedded as TrueType (editable text)
    <stem>.svg   vector, text kept as text
    <stem>.png   raster at --dpi (default 600)
    <stem>.csv   the rows plotted, with p, -log10 p and the score -- the
                 figure's source data, as journals ask for
    (tiff on request: LZW-compressed raster at --dpi)

Rendering uses matplotlib, a required dependency. Should it be missing (an
install made without dependencies), the step is skipped with a note and the
evidence is unaffected.
Proportions are fixed (``LAYOUT``), not fitted to the labels: every figure has
the same height (``--figure-top`` rows of one height, filled from the top) and
the same bar-panel width; pathways have a wider name column, so they are wider
(regulators 7.0 in, pathways 10.1 in) and line up with regulators side by side
at one height. Names wrap (pathways, two lines) or shorten to fit their column. The axis ends on a round number just past the longest bar and its
label, shared by both sides. The score sits at the bar's base in black, on a
lighter shade of the bar's colour. Text 7 pt in Arial or
its metric twins (Helvetica, Liberation Sans), then DejaVu Sans.
"""

from __future__ import annotations

import csv
import math
import re
import textwrap
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .errors import ProjectError
from .identity import DISPLAY_NAME, NAME

DEFAULT_COLORS = "#0044CC,#CC0000"     # a pure blue low, a pure red high
DEFAULT_FORMATS = ("pdf", "svg", "png")
KNOWN_FORMATS = ("pdf", "svg", "png", "tiff")
DEFAULT_DPI = 600
DEFAULT_TOP = 25
REFERENCE_P = 0.05          # drawn as a dashed line; selects nothing
LOG_P_CAP = 300.0                     # p = 0 (below double precision) is drawn at 300 and labelled ">300"

#: Fixed geometry, in inches, so figures keep one shape whatever their labels or values: every
#: figure is the same HEIGHT -- one row height, and always ``--figure-top`` rows (25), filled from the top -- and
#: the same bar-panel width, so the axes of a regulators and a pathways figure are drawn at one scale. Pathways
#: take a wider name column, so that figure is wider; set side by side at one height, the two line up.
ROW_IN = 0.22
LAYOUT = {"regulators": {"label": 0.85, "panel": 2.50},
          "pathways": {"label": 2.40, "panel": 2.50}}
MARGIN_IN = 0.10
CENTRE_GAP_IN = 0.10
LABEL_PAD_IN = 0.07          # between a name and its panel
P_LABEL_ROOM_IN = 0.30       # past the longest bar, for its -log10 p label
SHADE_MIX = 0.45             # the score box: 45% of the bar's colour, the rest white
NO_DIRECTION_FILL = "#B4B4B4"  # a score that is not calculable (NC); a score of 0 is white, on the scale
NO_DIRECTION_EDGE = "#7A7A7A"

FONT = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]
FIGURE_COLUMNS = ("side", "rank", "label", "source", "id", "p", "neg_log10_p", "score", "score_name", "detail")


@dataclass
class FigureSettings:
    colors: str = DEFAULT_COLORS
    top: int = DEFAULT_TOP
    formats: Tuple[str, ...] = DEFAULT_FORMATS
    dpi: int = DEFAULT_DPI
    exclude: Tuple[str, ...] = ()        # any of EXCLUDABLE
    pathway_score: str = "activation"    # PATHWAY_SCORES

    def as_dict(self) -> Dict[str, object]:
        return {"colors": self.colors, "top": self.top, "formats": list(self.formats), "dpi": self.dpi,
                "exclude": list(self.exclude), "pathway_score": self.pathway_score}

    @classmethod
    def from_dict(cls, d: Optional[Dict[str, object]]) -> "FigureSettings":
        d = d or {}
        return cls(colors=str(d.get("colors", DEFAULT_COLORS)), top=int(d.get("top", DEFAULT_TOP)),
                   formats=tuple(d.get("formats", DEFAULT_FORMATS)), dpi=int(d.get("dpi", DEFAULT_DPI)),
                   exclude=tuple(d.get("exclude", ())), pathway_score=str(d.get("pathway_score", "activation")))



#: With these and save_reproducibly, the same figure gives the same bytes on every run.
REPRODUCIBLE_RC = {"svg.hashsalt": NAME}


def save_reproducibly(fig, path: Path, fmt: str, dpi: int) -> None:
    """Save without a creation date, so a rerun writes identical files (shared by every figure)."""
    from .version import __version__ as version

    meta = {"pdf": {"Creator": f"{DISPLAY_NAME} {version}", "CreationDate": None, "ModDate": None},
            "svg": {"Creator": f"{DISPLAY_NAME} {version}", "Date": None},
            "png": {"Software": f"{DISPLAY_NAME} {version}"}}
    kwargs = {"dpi": dpi, "facecolor": "white"}
    if meta.get(fmt):
        kwargs["metadata"] = meta[fmt]
    if fmt == "tiff":
        kwargs["pil_kwargs"] = {"compression": "tiff_lzw"}
    fig.savefig(path, format=fmt, **kwargs)

def available() -> bool:
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return False
    return True


def parse_colors(text: str) -> Tuple[str, str]:
    """'navy,firebrick' -> ('navy', 'firebrick'): exactly two colours matplotlib understands. EVD-0026 otherwise."""
    parts = [p.strip() for p in str(text).split(",") if p.strip()]
    if len(parts) != 2:
        raise ProjectError(26, f"--colors needs exactly two colours, low then high (e.g. 'blue,red'); got {text!r}")
    if available():
        from matplotlib.colors import to_rgb
        for p in parts:
            try:
                to_rgb(p)
            except ValueError:
                raise ProjectError(26, f"--colors: {p!r} is not a colour name or #RRGGBB value") from None
    return parts[0], parts[1]


#: What ``--exclude`` can take off a figure, for a reader who does not want them: results with no sign. ``nc`` -- not calculable (grey); ``zero`` -- a score of
#: exactly 0 (white); ``unsigned`` -- both. Their places go to the next results by p.
EXCLUDABLE = ("nc", "zero", "unsigned")


def parse_exclude(text: Optional[str]) -> Tuple[str, ...]:
    parts = tuple(p.strip().lower() for p in str(text or "").split(",") if p.strip())
    bad = [p for p in parts if p not in EXCLUDABLE]
    if bad:
        raise ProjectError(26, f"--exclude: {', '.join(bad)} (choose from {', '.join(EXCLUDABLE)})")
    if "unsigned" in parts:
        return ("nc", "zero")
    return tuple(sorted(set(parts)))


def _drop(items: List["Item"], exclude: Sequence[str]) -> List["Item"]:
    return [i for i in items if not (("nc" in exclude and i.score is None) or
                                     ("zero" in exclude and i.score is not None and i.score == 0))]


def parse_formats(text: str) -> Tuple[str, ...]:
    parts = tuple(p.strip().lower().lstrip(".") for p in str(text).split(",") if p.strip())
    bad = [p for p in parts if p not in KNOWN_FORMATS]
    if bad or not parts:
        raise ProjectError(26, f"--figure-formats: {', '.join(bad) or 'nothing given'} "
                               f"(choose from {', '.join(KNOWN_FORMATS)})")
    return parts


# ---------------------------------------------------------------------------
# what to plot
# ---------------------------------------------------------------------------

@dataclass
class Item:
    label: str
    id: str
    p: float
    score: Optional[float]    # None: not calculable (NC); 0 is an answer
    detail: str = ""
    tag: str = ""             # the gene-set source, shown after the name and never cut off

    @property
    def directed(self) -> bool:
        return self.score is not None and self.score != 0

    @property
    def neg_log10_p(self) -> float:
        return LOG_P_CAP if self.p <= 0 else min(LOG_P_CAP, -math.log10(self.p))


@dataclass
class Panel:
    """One figure's content: the two sides, and the words that go on it."""
    stem: str
    kind: str            # regulators | pathways
    title: str
    x_label: str
    score_name: str
    left_title: str
    right_title: str
    left: List[Item] = field(default_factory=list)
    right: List[Item] = field(default_factory=list)
    #: One panel ranked by p, no sides and no score (a burden test); the items are in ``right``.
    single: bool = False


def _num(text: object) -> Optional[float]:
    try:
        v = float(str(text))
    except ValueError:
        return None
    return None if math.isnan(v) else v


def _rows(path: Path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def split_top(items: Sequence[Item], top: int) -> Tuple[List[Item], List[Item]]:
    """(negative side, positive side), ``top`` each, no cutoff.

    Every result is taken in order of p (ties: larger |score|, then name). A result with a signed score goes to its
    own side while that side has room. One with no sign -- a score of exactly 0 (drawn white) or a score that is not
    calculable (drawn grey, NC) -- goes on BOTH sides, wherever each has room, since it
    belongs to neither."""
    order = lambda i: (i.p, -abs(i.score or 0.0), i.label)  # noqa: E731
    neg: List[Item] = []
    pos: List[Item] = []
    for item in sorted(items, key=order):
        if len(neg) >= top and len(pos) >= top:
            break
        if item.directed:
            sides = [neg if item.score < 0 else pos]
        else:
            sides = [neg, pos]
        for side in sides:
            if len(side) < top:
                side.append(item)
    return neg, pos


#: How a gene-set source is named on a figure row.
#: Sources the pathways figure draws from. GO terms are left out as too generic to put in front of
#: more informative pathways; they stay in the tables.
FIGURE_SOURCES = ("reactome", "wikipathways", "hallmark")
SOURCE_LABEL = {"reactome": "Reactome", "go_bp": "GO BP", "go_mf": "GO MF", "go_cc": "GO CC",
                "wikipathways": "WikiPathways", "hallmark": "Hallmark"}


def _single(stem: str, kind: str, title: str, x_label: str, items: List[Item], top: int) -> Panel:
    """One panel ranked by p (ties: name), no sides, no score: a burden test's figure."""
    ranked = sorted(items, key=lambda i: (i.p, i.label))[:top]
    for item in ranked:
        item.score = None
    return Panel(stem, kind, title, x_label, "", "", "", [], ranked, single=True)


def regulator_panel(path: Path, comparison: str, top: int, exclude: Sequence[str] = (),
                    single: bool = False) -> Panel:
    """Every network in one figure. A regulator tested in both networks keeps its lowest-p row once.

    ``single`` (a burden test): one panel by overlap p, no activation z."""
    if single:
        best_p: Dict[str, Item] = {}
        for r in _rows(path):
            p = _num(r.get("overlap_p"))
            if p is None:
                continue
            item = Item(r["regulator"], r["regulator"], p, None,
                        f"{r['network']}; {r.get('significant_targets', '')} significant targets")
            if r["regulator"] not in best_p or p < best_p[r["regulator"]].p:
                best_p[r["regulator"]] = item
        return _single(f"{comparison}_regulators", "regulators", f"{comparison} — upstream regulators",
                       "−log$_{10}$(p-value of overlap)", list(best_p.values()), top)
    best: Dict[str, Item] = {}
    for r in _rows(path):
        p, z = _num(r.get("overlap_p")), _num(r.get("activation_z"))
        if p is None:
            continue
        item = Item(r["regulator"], r["regulator"], p, z,
                    f"{r['network']}; {r.get('n_consistent', '')} consistent / "
                    f"{r.get('n_inconsistent', '')} inconsistent targets")
        held = best.get(item.label)
        if held is None or (item.p, -abs(item.score or 0.0)) < (held.p, -abs(held.score or 0.0)):
            best[item.label] = item
    neg, pos = split_top(_drop(list(best.values()), exclude), top)
    return Panel(f"{comparison}_regulators", "regulators", f"{comparison} — upstream regulators",
                 "−log$_{10}$(p-value of overlap)", "Activation z-score",
                 "Predicted inhibited", "Predicted activated", neg, pos)


def _display_name(name: str, source: str) -> str:
    """MSigDB's HALLMARK_TNFA_SIGNALING_VIA_NFKB is drawn as 'TNFA SIGNALING VIA NFKB' (the row says [Hallmark]);
    the CSV keeps the set ID. Every other name is drawn as its source spells it."""
    if source == "hallmark" and name.upper().startswith("HALLMARK_"):
        return name[len("HALLMARK_"):].replace("_", " ")
    return name


#: What colours a pathway bar. ``activation``: the pathway activation z from Reactome's
#: activity patterns; a set with no pattern, or too few pattern genes, is NC (grey). ``direction``: (up - down) / sqrt(n), every set. When the run has no patterns at
#: all (``reference pathway-signs`` not built), ``activation`` falls back to ``direction`` and the legend says so.
PATHWAY_SCORES = ("activation", "direction")


def pathway_panel(path: Path, comparison: str, top: int, exclude: Sequence[str] = (),
                  score: str = "activation", single: bool = False) -> Panel:
    """Every gene-set source in one figure, by over-representation; each row names its source.

    ``single`` (a burden test): one panel by the ORA p of all significant genes, no score."""
    if single:
        items = []
        for r in _rows(path):
            if r.get("direction") != "all" or r["source"] not in FIGURE_SOURCES:
                continue
            p = _num(r.get("p"))
            if p is None:
                continue
            items.append(Item(_display_name(r["name"], r["source"]), r["set_id"], p, None,
                              f"{r['source']}; {r.get('overlap') or 0} of {r.get('set_size') or 0}; "
                              f"{r.get('null_model') or ''}", SOURCE_LABEL.get(r["source"], r["source"])))
        return _single(f"{comparison}_pathways", "pathways", f"{comparison} — pathways",
                       "−log$_{10}$(p-value of overlap)", items, top)
    rows: Dict[Tuple[str, str], Dict[str, Dict[str, str]]] = {}
    patterned = False
    for r in _rows(path):
        rows.setdefault((r["source"], r["set_id"]), {})[r["direction"]] = r
        patterned = patterned or bool(r.get("pattern_genes"))
    activation = score == "activation" and patterned
    items: List[Item] = []
    for (source, set_id), d in rows.items():
        if source not in FIGURE_SOURCES:
            continue
        if not {"all", "up", "down"} <= set(d):
            continue            # F-test comparisons have no direction: no pathway figure
        n = int(d["all"]["overlap"] or 0)
        up, down = int(d["up"]["overlap"] or 0), int(d["down"]["overlap"] or 0)
        p = _num(d["all"]["p"])
        if p is None:
            continue
        if activation:
            value = _num(d["all"].get("pathway_activation_z"))
            detail = (f"{source}; pattern: {d['all'].get('pattern_consistent') or 0} consistent / "
                      f"{d['all'].get('pattern_inconsistent') or 0} inconsistent")
        else:
            value = (up - down) / math.sqrt(n) if n else None
            detail = f"{source}; {up} up / {down} down of {n}"
        items.append(Item(_display_name(d["all"]["name"], source), set_id, p, value, detail,
                          SOURCE_LABEL.get(source, source)))
    neg, pos = split_top(_drop(items, exclude), top)
    if activation:
        return Panel(f"{comparison}_pathways", "pathways", f"{comparison} — pathways",
                     "−log$_{10}$(p-value of overlap)", "Activation z-score  (Reactome activity pattern)",
                     "Predicted inhibited", "Predicted activated", neg, pos)
    return Panel(f"{comparison}_pathways", "pathways", f"{comparison} — pathways",
                 "−log$_{10}$(p-value of overlap)", "z-score  (up − down) / √n",
                 "Decreased", "Increased", neg, pos)


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

def _wrap(label: str, tag: str = "", width: int = 40, lines: int = 2) -> str:
    """At most ``lines`` lines of ``width``; a long name is cut with an ellipsis, the [tag] always kept whole."""
    suffix = f" [{tag}]" if tag else ""
    parts = textwrap.wrap(label + suffix, width) or [label]
    if len(parts) <= lines:
        return "\n".join(parts)
    parts = textwrap.wrap(label, width)[:lines]
    room = width - len(suffix) - 1
    parts[-1] = parts[-1][:max(room, 8)].rstrip() + "…" + suffix
    return "\n".join(parts)


def _fit(label: str, tag: str, room_in: float, inches, lines: int) -> str:
    """The widest wrap of ``label`` (+ ``[tag]``) whose every line fits ``room_in`` inches, in at most ``lines``."""
    for width in range(48, 11, -1):
        text = _wrap(label, tag, width, lines)
        if all(inches(line) <= room_in for line in text.split("\n")):
            return text
    return _wrap(label, tag, 12, lines)


def nice_ceiling(value: float) -> float:
    """The smallest round number (1, 1.2, 1.4, 1.5, ... 9 x 10^k) at or above ``value``: where the axis ends."""
    if value <= 0:
        return 1.0
    k = 10 ** math.floor(math.log10(value))
    return next(m * k for m in (1, 1.2, 1.4, 1.5, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10)
                if m * k >= value - 1e-12)


def _fmt_p(v: float) -> str:
    return ">300" if v >= LOG_P_CAP else f"{v:.1f}"


def _fmt_score(v: Optional[float]) -> str:
    if v is None:
        return "NC"
    return f"{v:.2f}"


def draw(panel: Panel, out_dir: Path, settings: FigureSettings, *, version: str) -> List[Path]:
    """Render one panel in every requested format, plus its source-data CSV. Returns the files written."""
    if panel.single:
        return _draw_single(panel, out_dir, settings)
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
    from matplotlib.figure import Figure
    from matplotlib.cm import ScalarMappable

    low, high = parse_colors(settings.colors)
    cmap = LinearSegmentedColormap.from_list(NAME, [to_rgb(low), (1, 1, 1), to_rgb(high)])
    scores = [abs(i.score) for i in panel.left + panel.right if i.score is not None]
    zmax = max(2.0, math.ceil(max(scores))) if scores else 2.0      # the scale ends just past the largest score
    norm = Normalize(-zmax, zmax)
    xmax = max([i.neg_log10_p for i in panel.left + panel.right] + [1.0])

    pathway = panel.kind != "regulators"
    n = max(settings.top, len(panel.left), len(panel.right), 1)       # empty rows keep the height fixed
    geo = LAYOUT[panel.kind]
    label_w, panel_w, row_in = geo["label"], geo["panel"], ROW_IN
    width = MARGIN_IN * 2 + CENTRE_GAP_IN + 2 * (label_w + panel_w)
    rc = {**REPRODUCIBLE_RC, "font.family": "sans-serif", "font.sans-serif": FONT, "font.size": 7, "axes.linewidth": 0.5,
          "xtick.major.width": 0.5, "xtick.major.size": 2.5, "xtick.labelsize": 6.5, "pdf.fonttype": 42,
          "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": NAME, "mathtext.default": "regular"}
    written: List[Path] = []
    import logging
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    with matplotlib.rc_context(rc):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.ticker import MaxNLocator

        probe = Figure(figsize=(width, 2))
        renderer = FigureCanvasAgg(probe).get_renderer()

        def inches(text: str) -> float:
            t = probe.text(0, 0, text, fontsize=7)
            w = t.get_window_extent(renderer).width / probe.dpi
            t.remove()
            return w

        room = label_w - LABEL_PAD_IN
        labels_l = [_fit(i.label, i.tag, room, inches, 2 if pathway else 1) for i in panel.left]
        labels_r = [_fit(i.label, i.tag, room, inches, 2 if pathway else 1) for i in panel.right]

        top_in, bottom_in = 0.55, 0.95
        height = top_in + n * row_in + bottom_in
        fig = Figure(figsize=(width, height))
        fig.patch.set_facecolor("white")
        h = n * row_in / height
        y0 = bottom_in / height
        ax_l = fig.add_axes([MARGIN_IN / width, y0, panel_w / width, h])
        ax_r = fig.add_axes([(MARGIN_IN + panel_w + 2 * label_w + CENTRE_GAP_IN) / width, y0, panel_w / width, h])
        # both panels share one axis, ending at a round number with room for the p labels past the longest bar
        xlim = nice_ceiling(xmax * (1 + P_LABEL_ROOM_IN / panel_w * 1.1))
        per_inch = xlim / panel_w                   # data units per inch along a panel
        base = 0.03 * per_inch                      # the score box starts this far from the centre line

        for ax, items, labels, side in ((ax_l, panel.left, labels_l, "left"), (ax_r, panel.right, labels_r, "right")):
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_xlim((xlim, 0) if side == "left" else (0, xlim))
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
            for spine in ("top", "left", "right"):
                ax.spines[spine].set_visible(False)
            ax.tick_params(axis="y", length=0, pad=3)
            for y in range(n + 1):
                ax.axhline(y - 0.5, color="#E3E3E3", linewidth=0.4, zorder=0)
            ys = list(range(len(items)))
            ha = "right" if side == "left" else "left"
            for y, item in zip(ys, items):
                if item.score is None:              # not calculable: grey, "NC"
                    color, edge = to_rgb(NO_DIRECTION_FILL), to_rgb(NO_DIRECTION_EDGE)
                else:                               # 0 is an answer: white, the centre of the scale
                    color = cmap(norm(item.score))
                    edge = (to_rgb(low) if item.score < 0 else to_rgb(high)) if item.score else \
                        to_rgb(NO_DIRECTION_EDGE)
                ax.barh(y, item.neg_log10_p, height=0.72, color=color, edgecolor=edge, linewidth=0.4, zorder=2)
                # the score: black, on a lighter shade of the bar's own colour, at the bar's base;
                # the p label goes past the bar end, or past the score box if the bar is shorter
                score_text = _fmt_score(item.score)
                shade = tuple(SHADE_MIX * c + (1 - SHADE_MIX) for c in color[:3])
                ax.text(base, y, score_text, va="center", ha=ha, fontsize=6, color="black", zorder=3,
                        bbox={"boxstyle": "square,pad=0.18", "facecolor": shade, "edgecolor": "none"})
                box_end = base + (inches(score_text) * 6 / 7 + 0.05) * per_inch
                ax.text(max(item.neg_log10_p, box_end) + 0.03 * per_inch, y, _fmt_p(item.neg_log10_p),
                        va="center", ha=ha, fontsize=6, color="#333333", zorder=3)
            ax.set_yticks(ys)
            ax.set_yticklabels(labels, fontsize=7, linespacing=0.95)
            if side == "left":
                ax.yaxis.tick_right()
                for t in ax.get_yticklabels():
                    t.set_horizontalalignment("left")
                ax.yaxis.set_tick_params(pad=3)
            ax.set_xlabel(panel.x_label, fontsize=7, labelpad=2)
            ax.set_title(panel.left_title if side == "left" else panel.right_title, fontsize=7.5,
                         fontweight="bold", pad=4, loc="right" if side == "left" else "left",
                         color=to_rgb(low) if side == "left" else to_rgb(high))
            ref = -math.log10(REFERENCE_P)
            if ref < xlim:
                ax.axvline(ref, color="#808080", linewidth=0.5, linestyle=(0, (2.5, 2)), zorder=2.5)
        fig.text(0.5, 1 - 0.14 / height, panel.title, ha="center", va="top", fontsize=8)

        cax = fig.add_axes([(width / 2 - 0.9) / width, 0.30 / height, 1.8 / width, 0.09 / height])
        bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
        bar.outline.set_linewidth(0.4)
        step = 2 if zmax > 4 else 1
        bar.set_ticks([t for t in range(-int(zmax), int(zmax) + 1) if t % step == 0])
        bar.ax.tick_params(labelsize=6, width=0.4, length=2, pad=1.5)
        bar.set_label(panel.score_name, fontsize=6.5, labelpad=2)
        if any(i.score is None for i in panel.left + panel.right):
            key = fig.add_axes([(width / 2 + 1.05) / width, 0.30 / height, 0.16 / width, 0.09 / height])
            key.set_facecolor(NO_DIRECTION_FILL)
            for spine in key.spines.values():
                spine.set_linewidth(0.4)
                spine.set_edgecolor(NO_DIRECTION_EDGE)
            key.set_xticks([])
            key.set_yticks([])
            key.text(1.25, 0.5, "NC  not calculable", transform=key.transAxes, va="center", ha="left", fontsize=6)

        out_dir.mkdir(parents=True, exist_ok=True)
        for fmt in settings.formats:
            path = out_dir / f"{panel.stem}.{fmt}"
            save_reproducibly(fig, path, fmt, settings.dpi)
            written.append(path)

    data = out_dir / f"{panel.stem}.csv"
    with data.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=FIGURE_COLUMNS, lineterminator="\n")
        w.writeheader()
        for side, items in (("negative", panel.left), ("positive", panel.right)):
            for rank, i in enumerate(items, 1):
                w.writerow({"side": side, "rank": rank, "label": i.label, "source": i.tag, "id": i.id, "p": repr(i.p),
                            "neg_log10_p": round(i.neg_log10_p, 4), "score": "" if i.score is None else round(i.score, 4),
                            "score_name": panel.score_name, "detail": i.detail})
    written.append(data)
    return written


def _draw_single(panel: Panel, out_dir: Path, settings: FigureSettings) -> List[Path]:
    """A burden test's figure: one panel, rows ranked by p from the top, bar length -log10 p in the
    high colour, the p label past the bar; no sides, no score, no colour bar. Same row height, fonts, p = 0.05
    line and file set as the two-sided figures."""
    import logging

    import matplotlib
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.colors import to_rgb
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    _low, high = parse_colors(settings.colors)
    items = panel.right
    n = max(settings.top, len(items), 1)
    geo = LAYOUT[panel.kind]
    label_w, panel_w = geo["label"], geo["panel"]
    width = MARGIN_IN * 2 + label_w + panel_w
    rc = {**REPRODUCIBLE_RC, "font.family": "sans-serif", "font.sans-serif": FONT, "font.size": 7, "axes.linewidth": 0.5,
          "xtick.major.width": 0.5, "xtick.major.size": 2.5, "xtick.labelsize": 6.5, "pdf.fonttype": 42,
          "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": NAME, "mathtext.default": "regular"}
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    written: List[Path] = []
    with matplotlib.rc_context(rc):
        probe = Figure(figsize=(width, 2))
        renderer = FigureCanvasAgg(probe).get_renderer()

        def inches(text: str) -> float:
            t = probe.text(0, 0, text, fontsize=7)
            w = t.get_window_extent(renderer).width / probe.dpi
            t.remove()
            return w

        labels = [_fit(i.label, i.tag, label_w - LABEL_PAD_IN, inches, 2 if panel.kind != "regulators" else 1)
                  for i in items]
        top_in, bottom_in = 0.45, 0.55
        height = top_in + n * ROW_IN + bottom_in
        fig = Figure(figsize=(width, height))
        fig.patch.set_facecolor("white")
        ax = fig.add_axes([(MARGIN_IN + label_w) / width, bottom_in / height, panel_w / width, n * ROW_IN / height])
        xmax = max([i.neg_log10_p for i in items] + [1.0])
        xlim = nice_ceiling(xmax * (1 + P_LABEL_ROOM_IN / panel_w * 1.1))
        per_inch = xlim / panel_w
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xlim(0, xlim)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]))
        for spine in ("top", "left", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis="y", length=0, pad=3)
        for y in range(n + 1):
            ax.axhline(y - 0.5, color="#E3E3E3", linewidth=0.4, zorder=0)
        color = to_rgb(high)
        for y, item in enumerate(items):
            ax.barh(y, item.neg_log10_p, height=0.72, color=color, edgecolor=color, linewidth=0.4, zorder=2)
            ax.text(item.neg_log10_p + 0.03 * per_inch, y, _fmt_p(item.neg_log10_p), va="center", ha="left",
                    fontsize=6, color="#333333", zorder=3)
        ax.set_yticks(list(range(len(items))))
        ax.set_yticklabels(labels, fontsize=7, linespacing=0.95)
        ax.set_xlabel(panel.x_label, fontsize=7, labelpad=2)
        ref = -math.log10(REFERENCE_P)
        if ref < xlim:
            ax.axvline(ref, color="#808080", linewidth=0.5, linestyle=(0, (2.5, 2)), zorder=2.5)
        fig.text(0.5, 1 - 0.12 / height, panel.title, ha="center", va="top", fontsize=8)
        out_dir.mkdir(parents=True, exist_ok=True)
        for fmt in settings.formats:
            path = out_dir / f"{panel.stem}.{fmt}"
            save_reproducibly(fig, path, fmt, settings.dpi)
            written.append(path)
    data = out_dir / f"{panel.stem}.csv"
    with data.open("w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=FIGURE_COLUMNS, lineterminator="\n")
        w.writeheader()
        for rank, i in enumerate(items, 1):
            w.writerow({"side": "ranked by p", "rank": rank, "label": i.label, "source": i.tag, "id": i.id,
                        "p": repr(i.p), "neg_log10_p": round(i.neg_log10_p, 4), "score": "", "score_name": "",
                        "detail": i.detail})
    written.append(data)
    return written


_SAFE = re.compile(r"[^A-Za-z0-9._+-]+")


def panels_for(crosswalk_dir: Path, comparison: str, top: int, exclude: Sequence[str] = (),
               pathway_score: str = "activation", kind: Optional[str] = None) -> List[Panel]:
    """At most two figures per comparison: regulators, and pathways. ``kind`` is the
    comparison's kind (comparison_summary.json); a burden test draws the side(s) it has."""
    d = Path(crosswalk_dir)
    panels: List[Panel] = []
    single = kind == "burden"           # one panel ranked by p, no direction score
    reg = d / "regulators" / comparison / "regulators.csv"
    if reg.exists():
        panels.append(regulator_panel(reg, comparison, top, exclude, single=single))
    ora = d / "enrichment" / comparison / "ora.csv"
    if ora.exists():
        panels.append(pathway_panel(ora, comparison, top, exclude, pathway_score, single=single))
    if single:
        panels = [p for p in panels if p.right]
        for p in panels:
            p.stem = _SAFE.sub("_", p.stem)
        return panels
    # The same number of bars on each side, so the figure reads symmetrically:
    # each side keeps its k lowest-p rows, k = the shorter side. A one-sided result draws no figure (its rows stay
    # in the tables). A burden test was returned above: one panel ranked by p.
    for p in panels:
        k = min(len(p.left), len(p.right))
        p.left, p.right = p.left[:k], p.right[:k]
    panels = [p for p in panels if p.left or p.right]
    # A comparison where no result has a calculable score -- an F-test, which has no direction, so no activation
    # z for any regulator -- would be an all-grey figure that says nothing a table does not: none is drawn.
    panels = [p for p in panels if any(i.score is not None for i in p.left + p.right)]
    for p in panels:
        p.stem = _SAFE.sub("_", p.stem)
    return panels


def _draw_job(job) -> List[str]:
    panel, out_dir, settings, version = job
    return [str(p) for p in draw(panel, Path(out_dir), settings, version=version)]


def make_figures(crosswalk_dir: Path, settings: Optional[FigureSettings] = None, *, out: Optional[Path] = None,
                 threads: Optional[int] = None, required: bool = False, log=print) -> Dict[str, object]:
    """Every figure for every comparison in a crosswalk folder. Skips (with a note) when matplotlib is absent,
    unless ``required`` (the ``figures`` command), where that is EVD-0027."""
    import json

    from .parallel import map_ordered, resolve_threads
    from .version import __version__

    settings = settings or FigureSettings()
    parse_colors(settings.colors)
    crosswalk_dir = Path(crosswalk_dir)
    out = Path(out) if out else crosswalk_dir / "figures"
    if not available():
        if required:
            raise ProjectError(27, "matplotlib is not installed for this Python")
        log("  figures skipped: matplotlib is not installed for this Python (pip install matplotlib)")
        return {"skipped": "matplotlib not installed", "files": []}
    summary = crosswalk_dir / "comparison_summary.json"
    if not summary.exists():
        raise ProjectError(21, f"{summary} not found")
    described = json.loads(summary.read_text())["comparisons"]
    jobs = []
    for c in described:
        for p in panels_for(crosswalk_dir, c["comparison"], settings.top, settings.exclude, settings.pathway_score,
                            kind=c.get("kind")):
            if c.get("direction_from"):          # an F-test: say where its direction came from
                p.title += f"  (direction: {c['direction_from']})"
            if c.get("kind") == "burden":        # say what was compared
                p.title += f"  (burden test: {c.get('treatment') or 'cases'} vs {c.get('control') or 'reference'})"
            jobs.append((p, str(out), settings, __version__))
    results = map_ordered(_draw_job, jobs, threads=resolve_threads(threads)) if jobs else []
    files = [f for r in results for f in r]
    log(f"  figures/                     {len(jobs):>9} figures ({', '.join(settings.formats)} + csv), "
        f"colours {settings.colors}")
    return {"figures": len(jobs), "files": files, "dir": str(out)}
