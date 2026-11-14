"""Tests for the validation figure.

A figure cannot be unit-tested for being readable, so these test the things that
silently ruin one and that a glance at a thumbnail will not catch: geometry
escaping the canvas, a group quietly dropped, a median drawn in the wrong place,
and jitter that moves between runs.

The geometry test exists because of a real false alarm. A preview renderer
cropped the right-hand column and it looked like the figure overflowed; the
coordinates said otherwise. Measuring beats looking at a thumbnail.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest

from abscore.figure import (
    DOT_RADIUS,
    PANEL_HEIGHT,
    ROW_GAP,
    SERIES_COLOURS,
    Y_CEILING,
    Y_FLOOR,
    render_figure_document,
    render_figure_svg,
)
from abscore.validation import CHAIN_TYPES, FIGURE_ORIGINS, METRICS, ScoredChain

SVG = "{http://www.w3.org/2000/svg}"


def chain(record_id, origin, chain_type, v=0.9, fr=0.9, oasis=0.9):
    return ScoredChain(
        record_id=record_id,
        name=f"name-{record_id}",
        origin=origin,
        chain_type=chain_type,
        v_identity=v,
        framework_identity=fr,
        oasis=oasis,
        numbered_slots=120,
    )


@pytest.fixture
def rows():
    """At least one chain in every drawn (origin, chain type) cell, plus one
    chain the figure must not draw."""
    out = []
    for c_index, chain_type in enumerate(CHAIN_TYPES):
        for o_index, origin in enumerate(FIGURE_ORIGINS):
            for repeat in range(3):
                value = 0.95 - 0.05 * o_index - 0.01 * repeat
                out.append(
                    chain(
                        f"{chain_type}{o_index}{repeat}",
                        origin,
                        chain_type,
                        v=value,
                        fr=value,
                        oasis=value - 0.1 * c_index,
                    )
                )
    out.append(chain("vet1", "Caninized", "H", v=0.75))
    out.append(chain("mixed1", "Human/Humanized", "K", v=0.99))
    return out


@pytest.fixture
def tree(rows):
    return ET.fromstring(render_figure_svg(rows))


class TestItIsValidSvg:
    def test_it_parses(self, tree):
        assert tree.tag == f"{SVG}svg"

    def test_the_viewbox_matches_the_declared_size(self, tree):
        _, _, width, height = tree.get("viewBox").split()
        assert tree.get("width") == width
        assert tree.get("height") == height

    def test_it_is_labelled_for_a_screen_reader(self, tree):
        assert tree.get("role") == "img"
        assert tree.find(f"{SVG}title") is not None
        assert tree.find(f"{SVG}desc") is not None

    def test_nothing_is_loaded_from_outside(self, rows):
        svg = render_figure_svg(rows)
        for forbidden in ("http://", "https://", "<image", "xlink:href", "@import"):
            assert forbidden not in svg.replace("http://www.w3.org/2000/svg", "")


class TestWrapping:
    def test_it_breaks_at_words(self):
        from abscore.figure import _wrap

        assert _wrap("one two three four", 9) == ["one two", "three", "four"]

    def test_a_word_longer_than_the_budget_is_kept_whole(self):
        from abscore.figure import _wrap

        assert _wrap("supercalifragilistic x", 5) == ["supercalifragilistic", "x"]

    def test_empty_text_gives_no_lines(self):
        from abscore.figure import _wrap

        assert _wrap("   ", 10) == []


class TestNothingEscapesTheCanvas:
    """The check the eye cannot do."""

    def test_every_drawn_coordinate_is_inside_the_viewbox(self, tree):
        _, _, width, height = (float(v) for v in tree.get("viewBox").split())
        offenders = []
        for element in tree.iter():
            tag = element.tag.removeprefix(SVG)
            if tag == "circle":
                xs = [float(element.get("cx")) - DOT_RADIUS, float(element.get("cx")) + DOT_RADIUS]
                ys = [float(element.get("cy")) - DOT_RADIUS, float(element.get("cy")) + DOT_RADIUS]
            elif tag == "line":
                xs = [float(element.get("x1")), float(element.get("x2"))]
                ys = [float(element.get("y1")), float(element.get("y2"))]
            elif tag == "rect":
                left = float(element.get("x"))
                xs = [left, left + float(element.get("width"))]
                ys = [
                    float(element.get("y")),
                    float(element.get("y")) + float(element.get("height")),
                ]
            else:
                continue
            if min(xs) < 0 or max(xs) > width or min(ys) < 0 or max(ys) > height:
                offenders.append((tag, xs, ys))
        assert not offenders, offenders[:5]

    def test_no_note_line_runs_off_the_canvas(self, tree):
        """Anchors being inside is not enough -- a long string starts inside
        and ends outside, and SVG crops it without a word. The first version of
        this figure lost the end of both footnotes that way.

        The budget is the same characters-per-pixel estimate the wrapper uses,
        so this catches a note that grows past what the wrapper was told.
        """
        from abscore.figure import NOTE_CHAR_WIDTH

        _, _, width, _ = (float(v) for v in tree.get("viewBox").split())
        for element in tree.iter(f"{SVG}text"):
            if element.get("class") != "note":
                continue
            end = float(element.get("x")) + len(element.text) * NOTE_CHAR_WIDTH
            assert end <= width, (element.text, end, width)

    def test_the_canvas_grows_with_the_notes(self):
        """The notes are sized to the canvas, not trimmed to fit it."""
        from abscore.figure import render_figure_svg as render

        tree = ET.fromstring(render([chain("a", "Human", "H")]))
        _, _, _, height = (float(v) for v in tree.get("viewBox").split())
        lowest = max(float(t.get("y")) for t in tree.iter(f"{SVG}text") if t.get("class") == "note")
        assert lowest < height

    def test_text_anchors_sit_inside_the_canvas(self, tree):
        _, _, width, height = (float(v) for v in tree.get("viewBox").split())
        for element in tree.iter(f"{SVG}text"):
            assert 0 <= float(element.get("x")) <= width
            assert 0 <= float(element.get("y")) <= height


class TestEveryChainIsDrawnOnce:
    def test_one_point_per_chain_per_metric(self, rows, tree):
        drawn = len(list(tree.iter(f"{SVG}circle")))
        plotted = [r for r in rows if r.origin in FIGURE_ORIGINS]
        assert drawn == len(plotted) * len(METRICS)

    def test_labels_outside_the_figure_origins_are_left_out(self, rows, tree):
        """Scored but not plotted. The row counts must agree with the points,
        so a veterinary or mixed-origin chain cannot inflate an n."""
        titles = " ".join(t.text or "" for t in tree.iter(f"{SVG}title"))
        assert "Caninized" not in titles
        assert "Human/Humanized" not in titles

    def test_each_point_names_its_antibody(self, tree):
        circle = next(iter(tree.iter(f"{SVG}circle")))
        assert circle.find(f"{SVG}title").text.startswith("name-")

    def test_an_empty_group_keeps_its_place_on_the_axis(self):
        """The panel has no chimeric or murine lambda chains. Those cells stay
        on the axis, named and empty, because an absent group and a group with
        no separation look identical once the category is dropped."""
        svg = render_figure_svg([chain("a", "Human", "L", v=0.9)])
        for origin in FIGURE_ORIGINS:
            assert f">{origin}</text>" in svg

    def test_group_sizes_are_left_to_the_table(self):
        """Counts were drawn on the figure and then removed: under every panel
        they crowded the points, and drawn once they read as column totals
        when they described one row. The separation table carries them."""
        tree = ET.fromstring(render_figure_svg([chain("a", "Human", "L", v=0.9)]))
        drawn = [element.text or "" for element in tree.iter(f"{SVG}text")]
        assert not [text for text in drawn if "n=" in text], drawn

    def test_a_group_size_still_reaches_the_hover_text(self):
        """Removed from the drawing, kept on the median's tooltip: it costs no
        space there and answers "how many is this rule summarising?" at the
        moment the question arises."""
        tree = ET.fromstring(render_figure_svg([chain("a", "Human", "L", v=0.9)]))
        titles = [
            line.find(f"{SVG}title").text
            for line in tree.iter(f"{SVG}line")
            if line.find(f"{SVG}title") is not None
        ]
        assert any("n=1" in title for title in titles), titles


class TestTheMedianIsTheMedian:
    def test_the_rule_sits_at_the_median_value(self):
        values = [0.4, 0.6, 0.9]  # median 0.6
        rows = [chain(f"h{i}", "Human", "H", v=v) for i, v in enumerate(values)]
        tree = ET.fromstring(render_figure_svg(rows))
        medians = [
            line
            for line in tree.iter(f"{SVG}line")
            if line.get("class") == "median" and line.find(f"{SVG}title") is not None
        ]
        assert len(medians) == len(METRICS)

        y = float(medians[0].get("y1"))
        circles = sorted(
            (c for c in tree.iter(f"{SVG}circle") if float(c.get("cy")) > 0),
            key=lambda c: float(c.get("cy")),
        )
        middle = [c for c in circles if abs(float(c.get("cy")) - y) < 0.6]
        assert middle, "the median rule does not line up with the median point"

    def test_an_even_group_averages_the_two_middle_values(self):
        rows = [chain(f"h{i}", "Human", "H", v=v) for i, v in enumerate([0.4, 0.6, 0.8, 1.0])]
        tree = ET.fromstring(render_figure_svg(rows))
        median_line = next(
            line for line in tree.iter(f"{SVG}line") if line.get("class") == "median"
        )
        title = median_line.find(f"{SVG}title").text
        assert "0.700" in title


class TestTheScale:
    def test_the_floor_is_below_every_real_panel_score(self):
        """0.245 is the lowest score any committed panel chain produces. If a
        future change pushes something below the floor it would be drawn
        outside its panel, which this keeps honest."""
        assert Y_FLOOR < 0.245

    def test_a_point_at_the_ceiling_sits_on_the_panel_top(self):
        tree = ET.fromstring(render_figure_svg([chain("a", "Human", "H", v=Y_CEILING)]))
        circle = next(iter(tree.iter(f"{SVG}circle")))
        top = min(float(line.get("y1")) for line in tree.iter(f"{SVG}line"))
        assert abs(float(circle.get("cy")) - top) < 0.5

    def test_the_scale_is_shared_by_every_panel(self):
        """One value must land at one height everywhere, or reading across a
        row compares nothing."""
        same = [
            chain("a", "Human", "H", v=0.8, fr=0.8, oasis=0.8),
            chain("b", "Human", "K", v=0.8, fr=0.8, oasis=0.8),
        ]
        tree = ET.fromstring(render_figure_svg(same))
        heights = {
            round(float(c.get("cy")) % (PANEL_HEIGHT + ROW_GAP), 1)
            for c in tree.iter(f"{SVG}circle")
        }
        assert len(heights) == 1


class TestJitterIsReproducible:
    def test_two_renderings_are_byte_identical(self, rows):
        assert render_figure_svg(rows) == render_figure_svg(rows)

    def test_a_chain_keeps_its_offset_when_its_neighbours_change(self, rows):
        first = render_figure_svg(rows)
        trimmed = render_figure_svg(rows[:-5])
        marker = "<title>name-H00"
        assert marker in first and marker in trimmed

        # The circle before each title carries cx; compare the first occurrence.
        def cx_of(svg):
            index = svg.index(marker)
            return svg[:index].rsplit('cx="', 1)[1].split('"')[0]

        assert cx_of(first) == cx_of(trimmed)


class TestColour:
    def test_three_series_only(self):
        """Small multiples use the all-pairs pairlist, where only the first
        three palette slots clear the separation floors in both modes."""
        assert len(SERIES_COLOURS) == len(METRICS) == 3

    def test_both_colour_schemes_are_declared(self, rows):
        svg = render_figure_svg(rows)
        assert "prefers-color-scheme: dark" in svg
        assert '[data-theme="dark"]' in svg
        for light, dark in SERIES_COLOURS:
            assert light in svg and dark in svg

    def test_identity_never_rests_on_colour_alone(self, rows):
        """Each column carries its metric name as a direct label."""
        svg = render_figure_svg(rows)
        for _, short, _ in METRICS:
            assert f">{short}</text>" in svg


class TestTheDocument:
    def test_it_is_self_contained(self, rows):
        page = render_figure_document(rows, table_html="<table></table>")
        assert page.startswith("<!doctype html>")
        for forbidden in ("http://", "https://", "<script"):
            assert forbidden not in page.replace("http://www.w3.org/2000/svg", "")

    def test_the_table_is_included_beside_the_figure(self, rows):
        page = render_figure_document(rows, table_html="<table id='sep'></table>")
        assert "id='sep'" in page
        assert page.index("<svg") < page.index("id='sep'")
