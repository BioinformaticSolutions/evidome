"""A minimal .xlsx writer, standard library only.

Enough for a review workbook: several sheets, a bold frozen header row with an
autofilter, numbers stored as numbers, text as inline strings. Rows are
streamed into the zip, so a large sheet is never held as XML in memory.

Excel's limits are enforced, not discovered by the reader: at most 1,048,575
data rows per sheet (the caller is told how many were cut), 32,767 characters
per cell (longer text is cut and marked), sheet names of 31 characters without
[]:*?/\\ .
"""

from __future__ import annotations

import math
import re
import zipfile
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape

MAX_ROWS = 1_048_575
MAX_CELL = 32_767
_ILLEGAL = re.compile("[\x00-\x08\x0b\x0c\x0e-\x1f]")
_SHEET_BAD = re.compile(r"[\[\]:*?/\\]")

_CONTENT_TYPES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>
{sheets}
</Types>"""
_ROOT_RELS = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>"""
_STYLES = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>
<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>
<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="3"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/><xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1"><alignment wrapText="1" vertical="top"/></xf></cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>"""


def column_letter(index: int) -> str:
    """0 -> A, 25 -> Z, 26 -> AA."""
    letters = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        letters = chr(65 + rem) + letters
    return letters


def sheet_name(name: str, used: set) -> str:
    base = _SHEET_BAD.sub("_", name)[:31] or "Sheet"
    candidate, n = base, 2
    while candidate.lower() in used:
        suffix = f"~{n}"
        candidate, n = base[: 31 - len(suffix)] + suffix, n + 1
    used.add(candidate.lower())
    return candidate


def _number(value: object) -> Optional[str]:
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return repr(value) if math.isfinite(value) else None
    text = str(value).strip()
    # "_" too: Python's float() reads "1_2" as 12, Excel does not, and such a cell makes Excel offer a repair.
    if not text or text[0].isalpha() or text.startswith(("0x", "+")) or "_" in text:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    if not math.isfinite(number) or (len(text) > 1 and text[0] == "0" and text[1].isdigit()):
        return None  # keep leading-zero codes as text
    return text


def _cell(ref: str, value: object, style: int, text_columns: bool) -> str:
    if value is None or value == "":
        return ""
    if not text_columns:
        number = _number(value)
        if number is not None:
            style_attr = f' s="{style}"' if style else ""
            return f'<c r="{ref}"{style_attr}><v>{number}</v></c>'
    text = _ILLEGAL.sub("", str(value))
    if len(text) > MAX_CELL:
        text = text[: MAX_CELL - 30] + " [cut: Excel cell limit]"
    style_attr = f' s="{style}"' if style else ""
    return f'<c r="{ref}" t="inlineStr"{style_attr}><is><t xml:space="preserve">{escape(text)}</t></is></c>'


class Sheet:
    def __init__(self, name: str, columns: Sequence[str], rows: Iterable[Sequence[object]],
                 text_columns: Sequence[str] = (), widths: Optional[Sequence[int]] = None, wrap: bool = False,
                 header: bool = True):
        self.name, self.columns, self.rows = name, list(columns), rows
        self.text_columns = set(text_columns)
        self.widths, self.wrap, self.header = widths, wrap, header
        self.written = 0
        self.cut = 0


#: Every member of a zip written here carries this time, so the same content gives the same bytes.
FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


def _member(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, date_time=FIXED_ZIP_TIME)
    info.compress_type, info.external_attr = zipfile.ZIP_DEFLATED, 0o644 << 16
    return info


def write_workbook(path: Path, sheets: List[Sheet]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    used: set = set()
    names = [sheet_name(s.name, used) for s in sheets]
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(_member("[Content_Types].xml"), _CONTENT_TYPES.format(sheets="\n".join(
            f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" '
            f'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for i in range(len(sheets)))))
        zf.writestr(_member("_rels/.rels"), _ROOT_RELS)
        zf.writestr(_member("xl/styles.xml"), _STYLES)
        zf.writestr(_member("xl/workbook.xml"),
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
                    + "".join(f'<sheet name="{escape(n, {chr(34): "&quot;"})}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
                              for i, n in enumerate(names))
                    + "</sheets></workbook>")
        zf.writestr(_member("xl/_rels/workbook.xml.rels"),
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                    + "".join(f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/'
                              f'officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                              for i in range(len(sheets)))
                    + f'<Relationship Id="rId{len(sheets) + 1}" Type="http://schemas.openxmlformats.org/'
                      f'officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>')
        for i, sheet in enumerate(sheets):
            _write_sheet(zf, f"xl/worksheets/sheet{i + 1}.xml", sheet)
    tmp.replace(path)
    return path


def _write_sheet(zf: zipfile.ZipFile, member: str, sheet: Sheet) -> None:
    ncol = max(1, len(sheet.columns))
    last = column_letter(ncol - 1)
    text_idx = {i for i, c in enumerate(sheet.columns) if c in sheet.text_columns}
    with zf.open(_member(member), "w", force_zip64=True) as raw:
        def put(text: str) -> None:
            raw.write(text.encode("utf-8"))

        put('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">')
        if sheet.header:
            put('<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" '
                'state="frozen"/></sheetView></sheetViews>')
        widths = sheet.widths or [min(60, max(10, len(c) + 2)) for c in sheet.columns]
        put("<cols>" + "".join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>'
                               for i, w in enumerate(widths)) + "</cols><sheetData>")
        row_no = 0
        if sheet.header:
            row_no = 1
            put('<row r="1">' + "".join(_cell(f"{column_letter(i)}1", c, 1, True)
                                        for i, c in enumerate(sheet.columns)) + "</row>")
        style = 2 if sheet.wrap else 0
        for values in sheet.rows:
            if sheet.written >= MAX_ROWS:
                sheet.cut += 1
                continue
            row_no += 1
            sheet.written += 1
            put(f'<row r="{row_no}">' + "".join(
                _cell(f"{column_letter(i)}{row_no}", v, style, i in text_idx) for i, v in enumerate(values))
                + "</row>")
        put("</sheetData>")
        if sheet.header and row_no > 1:
            put(f'<autoFilter ref="A1:{last}{row_no}"/>')
        put("</worksheet>")
