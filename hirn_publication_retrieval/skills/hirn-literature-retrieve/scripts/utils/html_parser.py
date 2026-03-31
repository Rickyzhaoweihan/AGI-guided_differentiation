"""HTML parser for HIRN publication entries.

Parses the AJAX HTML response from the HIRN WordPress site
(hirnetwork.org/published-works) and extracts structured publication
metadata from each ``<div class="row">`` block.

Uses only Python standard-library modules (re, html.parser).
"""
from __future__ import annotations

import re
from html.parser import HTMLParser
from typing import List


# ---------------------------------------------------------------------------
# Regex helpers
# ---------------------------------------------------------------------------
_PMID_URL_RE = re.compile(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)")
_DOI_RE = re.compile(r"doi:\s*(10\.\S+?)[\s.<]")
_PMID_TEXT_RE = re.compile(r"PMID:\s*(\d+)")
_JOURNAL_RE = re.compile(r"\.\s*([A-Z][^.]+?)\.\s*\d{4}")


# ---------------------------------------------------------------------------
# Lightweight HTML-to-text helper using stdlib HTMLParser
# ---------------------------------------------------------------------------
class _TagStripper(HTMLParser):
    """Strip HTML tags and return plain text."""

    def __init__(self) -> None:
        super().__init__()
        self._pieces: list[str] = []

    def handle_data(self, data: str) -> None:
        self._pieces.append(data)

    def get_text(self) -> str:
        return "".join(self._pieces)


def _strip_tags(html: str) -> str:
    s = _TagStripper()
    s.feed(html)
    return s.get_text()


# ---------------------------------------------------------------------------
# Row-level field extractor using HTMLParser
# ---------------------------------------------------------------------------
class _RowParser(HTMLParser):
    """Parse a single ``<div class="row">`` block and extract fields."""

    def __init__(self) -> None:
        super().__init__()
        self.date: str = ""
        self.consortia_raw: str = ""
        self.title: str = ""
        self.pmid: str = ""
        self.citation_html_pieces: list[str] = []

        # Internal state
        self._col_sm2_count = 0
        self._in_col_sm2 = False
        self._in_pub_title = False
        self._in_title_link = False
        self._title_link_href: str = ""
        self._in_citation_p = False
        self._in_col_sm4 = False
        self._col_sm4_depth = 0
        self._pub_title_seen = False

    # -- tag handlers -------------------------------------------------------

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_dict = dict(attrs)
        cls = attr_dict.get("class", "") or ""

        if tag == "div" and "col-sm-2" in cls:
            self._col_sm2_count += 1
            self._in_col_sm2 = True

        elif tag == "div" and "col-sm-4" in cls:
            self._in_col_sm4 = True
            self._col_sm4_depth = 1

        elif self._in_col_sm4 and tag == "div":
            self._col_sm4_depth += 1

        elif tag == "p" and "pub-title" in cls:
            self._in_pub_title = True

        elif tag == "a" and self._in_pub_title:
            self._in_title_link = True
            href = attr_dict.get("href", "")
            self._title_link_href = href or ""
            m = _PMID_URL_RE.search(self._title_link_href)
            if m:
                self.pmid = m.group(1)

        elif tag == "p" and self._in_col_sm4 and self._pub_title_seen and not self._in_pub_title:
            self._in_citation_p = True

        elif tag == "br" and self._in_citation_p:
            # Preserve line breaks as newlines for later parsing
            self.citation_html_pieces.append("\n")

    def handle_endtag(self, tag: str) -> None:
        if tag == "div" and self._in_col_sm2:
            self._in_col_sm2 = False

        elif tag == "div" and self._in_col_sm4:
            self._col_sm4_depth -= 1
            if self._col_sm4_depth <= 0:
                self._in_col_sm4 = False

        elif tag == "p" and self._in_pub_title:
            self._in_pub_title = False
            self._pub_title_seen = True

        elif tag == "a" and self._in_title_link:
            self._in_title_link = False

        elif tag == "p" and self._in_citation_p:
            self._in_citation_p = False

    def handle_data(self, data: str) -> None:
        if self._in_col_sm2:
            if self._col_sm2_count == 1:
                self.date += data.strip()
            elif self._col_sm2_count == 2:
                self.consortia_raw += data.strip()

        elif self._in_title_link:
            self.title += data

        elif self._in_citation_p:
            self.citation_html_pieces.append(data)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_hirn_publications(html: str) -> list[dict]:
    """Parse HIRN AJAX HTML and return a list of publication dicts.

    Each dict has keys: ``pmid``, ``doi``, ``title``, ``authors``,
    ``journal``, ``date``, ``consortia``.

    Parameters
    ----------
    html:
        Raw HTML string returned by the HIRN AJAX endpoint.

    Returns
    -------
    list[dict]
        Parsed publications. Returns an empty list when *html* is empty
        or contains no ``<div class="row">`` blocks.
    """
    if not html or not html.strip():
        return []

    # Split the HTML into individual row blocks.
    # We split on '<div class="row">' and discard the first (empty) piece.
    row_chunks = re.split(r'<div\s+class="row">', html)

    publications: list[dict] = []

    for chunk in row_chunks:
        # Skip empty / preamble chunks that don't contain publication data
        if "col-sm-4" not in chunk:
            continue

        # Re-add the opening tag so the parser sees valid structure
        row_html = '<div class="row">' + chunk

        parser = _RowParser()
        parser.feed(row_html)

        # -- title -----------------------------------------------------------
        title = " ".join(parser.title.split()).strip().rstrip(".")

        # -- consortia -------------------------------------------------------
        consortia = [c.strip() for c in parser.consortia_raw.split(",") if c.strip()]

        # -- citation text ---------------------------------------------------
        citation_text = "".join(parser.citation_html_pieces)
        citation_text = " ".join(citation_text.split())  # normalise whitespace

        # -- doi -------------------------------------------------------------
        doi = ""
        m = _DOI_RE.search(citation_text)
        if m:
            doi = m.group(1).rstrip(".")

        # -- pmid (fallback to citation text if URL didn't work) -------------
        pmid = parser.pmid
        if not pmid:
            m = _PMID_TEXT_RE.search(citation_text)
            if m:
                pmid = m.group(1)

        # -- authors ---------------------------------------------------------
        # Authors are everything before the first ". <Capital>" which starts
        # the journal name.
        authors = ""
        author_match = re.match(r"(.+?)\.\s+(?=[A-Z])", citation_text)
        if author_match:
            authors = author_match.group(1).strip()

        # -- journal ---------------------------------------------------------
        journal = ""
        j_match = _JOURNAL_RE.search(citation_text)
        if j_match:
            journal = j_match.group(1).strip()

        publications.append(
            {
                "pmid": pmid,
                "doi": doi,
                "title": title,
                "authors": authors,
                "journal": journal,
                "date": parser.date,
                "consortia": consortia,
            }
        )

    return publications
