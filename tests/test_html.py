"""Tests for the HTML report.

Markdown cannot draw the position track, and raw pipes are hard to read in an
editor. This renderer exists for the file a reader actually opens, so the tests
cover what that implies: valid structure, no external assets, both colour
schemes, and the same numbers as the other two renderers.
"""

from __future__ import annotations

import re
from html.parser import HTMLParser

import pytest

from abscore.html import render_html, render_html_document
from abscore.report import analyse, column_note, render_markdown, render_text
from tests.fixtures import ADALIMUMAB_VH, TRASTUZUMAB_VH, TRASTUZUMAB_VL

CACHE = "tests/data/iedb"
VOID = {"meta", "br", "hr", "img", "input", "link"}


def run(sequence, **kwargs):
    return analyse(sequence, cache_dir=CACHE, offline=True, **kwargs)


class _Balance(HTMLParser):
    """Enough of a parser to catch a tag left open by a formatting slip."""

    def __init__(self):
        super().__init__()
        self.stack: list[str] = []
        self.mismatched: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag not in VOID:
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if self.stack and self.stack[-1] == tag:
            self.stack.pop()
        else:
            self.mismatched.append(tag)


@pytest.fixture(scope="module")
def scored():
    return run(TRASTUZUMAB_VH, with_epitopes=True)


@pytest.fixture(scope="module")
def page(scored):
    return render_html(scored)


class TestThePageIsWellFormedAndSelfContained:
    def test_tags_are_balanced(self, page):
        check = _Balance()
        check.feed(page)
        assert not check.mismatched, check.mismatched
        assert not check.stack, check.stack

    def test_nothing_is_fetched_from_elsewhere(self, page):
        """It has to open from a file, offline, with no assets beside it."""
        assert "src=" not in page
        assert "href=" not in page
        assert "http://" not in page and "https://" not in page

    def test_both_colour_schemes_are_declared(self, page):
        """The media query covers the OS setting; the data-theme scope covers a
        theme toggle, and the toggle has to win either way."""
        assert "prefers-color-scheme: dark" in page
        assert ':root[data-theme="dark"]' in page

    def test_the_position_track_is_embedded(self, page):
        """The thing markdown structurally cannot do."""
        assert "ab-track" in page
        assert "ab-cell" in page

    def test_every_section_is_present(self, page):
        assert re.findall(r"<h2>([^<]+)</h2>", page) == [
            "Chain",
            "Regions",
            "Humanness",
            "MHC class II",
            "Caveats",
        ]


class TestSequenceDataIsEscaped:
    def test_markup_in_a_name_cannot_inject(self):
        """A FASTA header is attacker-controlled text from a file."""
        page = render_html(run(TRASTUZUMAB_VH, with_epitopes=False), title="<script>x</script>")
        assert "<script>x</script>" not in page
        assert "&lt;script&gt;" in page


class TestTheRenderersAgree:
    """Three renderings of one report. The numbers must not drift apart."""

    @pytest.mark.parametrize("sequence", [TRASTUZUMAB_VH, TRASTUZUMAB_VL, ADALIMUMAB_VH])
    def test_the_headline_numbers_match_across_all_three(self, sequence):
        report = run(sequence, with_epitopes=False)
        renderings = (
            render_html(report),
            render_markdown(report),
            render_text(report),
        )
        wanted = [
            report.germline.gene,
            f"{report.germline.identity:.1%}",
            f"{report.humanness.score:.1%}",
        ]
        for value in wanted:
            for rendering in renderings:
                assert value in rendering, (value, rendering[:80])

    def test_the_caveat_count_matches_across_all_three(self):
        report = run(TRASTUZUMAB_VL, with_epitopes=False)
        count = len(report.warnings())
        assert f"Caveats ({count})" in render_html(report)
        assert f"CAVEATS ({count})" in render_text(report)

    def test_the_core_counts_match_the_text_rendering(self, scored):
        cores = scored.epitopes.core_counts_by_definition()
        page = render_html(scored)
        text = render_text(scored)
        for definition in ("germline", "repertoire", "proteome"):
            unmatched = str(cores[f"{definition}_unmatched"])
            assert unmatched in page
            assert unmatched in text


class TestABatchBecomesOnePage:
    def test_each_chain_gets_its_own_heading(self):
        reports = [
            ("humanized", run(TRASTUZUMAB_VH, with_epitopes=False)),
            ("human", run(ADALIMUMAB_VH, with_epitopes=False)),
        ]
        page = render_html_document(reports)
        assert re.findall(r"<h1>([^<]+)</h1>", page) == ["humanized", "human"]

    def test_a_single_chain_batch_renders_as_one_report(self):
        page = render_html_document([("query", run(TRASTUZUMAB_VH, with_epitopes=False))])
        assert page.count("<h2>Chain</h2>") == 1

    def test_the_stylesheets_do_not_repeat_per_chain(self):
        """The track ships its own CSS, which a batch would otherwise repeat
        once per chain -- ~1.5 kB each, in a file someone emails. The document
        emits it once and renders the tracks bare."""
        one = render_html_document([("a", run(TRASTUZUMAB_VH, with_epitopes=False))])
        two = render_html_document(
            [
                ("a", run(TRASTUZUMAB_VH, with_epitopes=False)),
                ("b", run(ADALIMUMAB_VH, with_epitopes=False)),
            ]
        )
        for page in (one, two):
            assert page.count("<style>") == 1
            assert page.count(".fields .k") == 1
        # Same count of track-style rules whether there is one chain or two.
        assert one.count(".ab-track {") == two.count(".ab-track {")
        # ...and the tracks themselves still render.
        assert two.count('<div class="ab-track">') == 2


class TestAnUnavailableEpitopeSection:
    def test_it_says_so_and_keeps_the_rest(self):
        page = render_html(run(TRASTUZUMAB_VH, with_epitopes=False))
        assert "Not available" in page
        assert "<h2>Humanness</h2>" in page
        assert "computed offline" in page


class TestEveryColumnCarriesItsDefinition:
    """A reader should not have to ask what a column means.

    The definitions live in `report.COLUMN_NOTES` and reach this page as
    `title` attributes. The risk being guarded against is silent drift: a
    column renamed here, or a new one added, stops matching its entry and the
    tooltip quietly disappears with nothing failing.
    """

    def test_every_header_in_the_page_has_a_title(self):
        page = render_html(run(TRASTUZUMAB_VH))
        # One regex capturing the tag and its label together, so the two can
        # never be matched against each other out of step.
        # The lookahead is load-bearing: `<th[^>]*>` also matches `<thead>`.
        headers = re.findall(r"(<th(?=[\s>])[^>]*>)(.*?)</th>", page, flags=re.DOTALL)
        assert headers, "no table headers found; the selector is wrong"
        untitled = [label for tag, label in headers if "title=" not in tag]
        assert not untitled, f"columns with no definition on hover: {untitled}"

    def test_the_definitions_are_the_shared_ones(self):
        page = render_html(run(TRASTUZUMAB_VH))
        assert 'title="exact 9-mer is present in that reference set"' in page
        assert "lookup result, not a risk score" in page

    def test_a_reference_suffixed_label_still_finds_its_note(self):
        """The app's batch column is "Unmatched (germline)", not "Unmatched"."""
        assert column_note("Unmatched (germline)") == column_note("Not found")
        assert column_note("Matched") == column_note("Found")

    def test_an_unknown_label_renders_a_plain_header(self):
        assert column_note("Something Nobody Defined") is None
