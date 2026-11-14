"""Tests for the panel comparison.

These use hand-built rows rather than the real panel, so they run offline and in
milliseconds. The real-panel assertions -- that human chains actually separate
from murine ones, and that OASis does not beat plain germline identity -- live in
`tests/test_panel.py`, which is marked slow.

What is checked here is the arithmetic and the reporting, because a validation
section is the part of a project a reader is least able to check for themselves.
An AUC computed on the wrong subset, or a tie printed as a win, would be
invisible.
"""

from __future__ import annotations

import pytest

from abscore.validation import (
    METRICS,
    SeparationRow,
    auc,
    highest_scoring,
    overlapping_chains,
    render_separation_html,
    render_separation_text,
    separation_table,
)
from tests.test_figure import chain


class TestAuc:
    def test_perfect_separation_is_one(self):
        assert auc([0.9, 0.8], [0.3, 0.2]) == 1.0

    def test_reversed_separation_is_zero(self):
        assert auc([0.2], [0.9]) == 0.0

    def test_no_information_is_a_half(self):
        assert auc([0.5, 0.5], [0.5, 0.5]) == 0.5

    def test_a_tie_counts_as_half_a_win(self):
        """Not a win. Germline identity is quantised -- a V gene of ~98
        positions gives identities on a coarse grid, so exact ties between a
        human and a murine chain are common rather than exotic. Counting them
        as wins would inflate every number in the table."""
        assert auc([0.5], [0.5]) == 0.5
        assert auc([0.9, 0.5], [0.5]) == 0.75

    def test_an_empty_side_is_an_error_not_a_default(self):
        """Returning 0.5, or 0.0, for an unmeasurable comparison would put a
        number in the table where there is no evidence."""
        with pytest.raises(ValueError):
            auc([], [0.5])
        with pytest.raises(ValueError):
            auc([0.5], [])


class TestSeparationTable:
    @pytest.fixture
    def rows(self):
        out = []
        for i in range(6):
            out.append(chain(f"h{i}", "Human", "H", v=0.95, fr=0.95, oasis=0.90))
            out.append(chain(f"m{i}", "Mouse", "H", v=0.70, fr=0.70, oasis=0.60))
            out.append(chain(f"hk{i}", "Human", "K", v=0.95, fr=0.95, oasis=0.90))
            out.append(chain(f"mk{i}", "Chimeric", "K", v=0.70, fr=0.70, oasis=0.60))
        return out

    def test_there_is_a_pooled_row_and_one_per_chain_type(self, rows):
        table = separation_table(rows)
        assert [row.group for row in table] == ["all", "H", "K", "L"]

    def test_each_row_counts_only_its_own_chain_type(self, rows):
        by_group = {row.group: row for row in separation_table(rows)}
        assert (by_group["H"].positives, by_group["H"].negatives) == (6, 6)
        assert (by_group["all"].positives, by_group["all"].negatives) == (12, 12)

    def test_chimeric_chains_are_negatives(self, rows):
        """A chimeric antibody carries a murine V domain, which is what both
        metrics measure. Putting it with the positives would be scoring the
        constant region, which neither metric sees."""
        by_group = {row.group: row for row in separation_table(rows)}
        assert by_group["K"].negatives == 6

    def test_a_group_too_small_to_measure_reports_counts_and_no_auc(self, rows):
        """The real panel has no chimeric or murine lambda chains at all. That
        row must print with its zero rather than disappear."""
        by_group = {row.group: row for row in separation_table(rows)}
        assert by_group["L"].aucs == {}
        assert by_group["L"].positives == 0

    def test_the_minimum_group_size_is_enforced_on_both_sides(self):
        rows = [chain(f"h{i}", "Human", "H", v=0.9) for i in range(20)]
        rows += [chain("m0", "Mouse", "H", v=0.5)]
        by_group = {row.group: row for row in separation_table(rows)}
        assert by_group["H"].aucs == {}, "one negative is not a measurement"

    def test_every_metric_is_scored_when_the_groups_are_big_enough(self, rows):
        by_group = {row.group: row for row in separation_table(rows)}
        assert set(by_group["H"].aucs) == {attribute for attribute, _, _ in METRICS}


class TestBestMetric:
    def test_it_names_the_highest(self):
        row = SeparationRow("all", "all", 5, 5, {"v_identity": 0.9, "oasis": 0.8})
        assert row.best_metric == "v_identity"

    def test_an_unmeasured_row_has_no_best(self):
        assert SeparationRow("L", "lambda", 0, 0, {}).best_metric is None


class TestRendering:
    @pytest.fixture
    def table(self):
        return [
            SeparationRow("all", "all chains", 400, 152, {a: 0.9 for a, _, _ in METRICS}),
            SeparationRow("L", "lambda only", 49, 0, {}),
        ]

    def test_the_text_table_marks_the_best_metric(self):
        table = [SeparationRow("all", "all chains", 5, 5, {"v_identity": 0.99, "oasis": 0.98})]
        text = render_separation_text(table)
        assert "0.99000*" in text
        assert "0.98000 " in text

    def test_five_decimals_survive(self):
        """The real kappa row separates by less than 1e-4. Rounding to three or
        four decimals would print a tie as a tie when it is a win, or the
        reverse."""
        table = [SeparationRow("K", "kappa only", 5, 5, {"v_identity": 0.994911})]
        assert "0.99491" in render_separation_text(table)

    def test_an_unmeasured_row_says_so_in_both_renderers(self, table):
        assert "too few to measure" in render_separation_text(table)
        assert "too few to measure" in render_separation_html(table)

    def test_the_text_table_fits_the_terminal(self, table):
        assert all(len(line) <= 88 for line in render_separation_text(table).splitlines())

    def test_the_html_table_carries_a_caption_explaining_auc(self, table):
        html = render_separation_html(table)
        assert "<caption>" in html
        assert "0.5 is no information" in html

    def test_the_html_row_count_matches_the_table(self, table):
        assert render_separation_html(table).count("<tr>") == len(table) + 1


class TestOverlap:
    def test_it_finds_positives_below_the_best_negative(self):
        rows = [
            chain("good", "Human", "H", v=0.95),
            chain("low", "Human", "H", v=0.70),
            chain("mouse", "Mouse", "H", v=0.80),
        ]
        assert [c.record_id for c in overlapping_chains(rows)] == ["low"]

    def test_a_tie_with_the_best_negative_counts_as_overlap(self):
        rows = [chain("tie", "Human", "H", v=0.80), chain("mouse", "Mouse", "H", v=0.80)]
        assert [c.record_id for c in overlapping_chains(rows)] == ["tie"]

    def test_chimeric_chains_do_not_set_the_ceiling_by_default(self):
        """One chimeric kappa chain in the panel reaches 0.915 V identity, near
        the human median. Including it would sweep a quarter of the human
        chains into this list and say nothing about the metric."""
        rows = [
            chain("human", "Human", "K", v=0.90),
            chain("chimeric", "Chimeric", "K", v=0.95),
            chain("mouse", "Mouse", "K", v=0.80),
        ]
        assert overlapping_chains(rows) == []

    def test_it_is_ordered_worst_first(self):
        rows = [
            chain("a", "Human", "H", v=0.78),
            chain("b", "Human", "H", v=0.70),
            chain("mouse", "Mouse", "H", v=0.80),
        ]
        assert [c.record_id for c in overlapping_chains(rows)] == ["b", "a"]

    def test_no_negatives_means_no_claim(self):
        rows = [chain("a", "Human", "L", v=0.70)]
        assert overlapping_chains(rows) == []


class TestHighestScoring:
    def test_it_names_the_chain(self):
        rows = [chain("a", "Mouse", "H", v=0.7), chain("b", "Mouse", "H", v=0.81)]
        assert highest_scoring(rows, origins=frozenset({"Mouse"})).record_id == "b"

    def test_an_absent_group_gives_nothing_rather_than_a_guess(self):
        rows = [chain("a", "Human", "H", v=0.9)]
        assert highest_scoring(rows, origins=frozenset({"Mouse"})) is None
