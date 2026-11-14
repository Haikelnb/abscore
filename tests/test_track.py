"""Tests for the IMGT position track.

The track exists to answer one question a table cannot: "CDR2 is 56-65, so why
are there three letters?" It draws every slot, so a short loop, a structural
gap, a truncated end and a real deletion are visibly different things.

That distinction is the whole value and it is also the easiest thing to get
wrong. The first version of this renderer drew a short CDR with the same mark
as a framework hole, which would have flagged every antibody ever made as
defective.
"""

from __future__ import annotations

import pytest

from abscore.number import ANCHOR_POSITIONS, STRUCTURAL_FRAMEWORK_GAPS, number_sequence
from abscore.track import _rows, render_track_html
from tests.fixtures import TRASTUZUMAB_VH, TRASTUZUMAB_VL


def rows(sequence: str):
    return _rows(number_sequence(sequence)[0])


def kinds(sequence: str) -> set[str]:
    return {kind for _, _, _, cells in rows(sequence) for kind, _ in cells}


class TestWhatEachSlotIs:
    def test_a_clean_chain_has_no_unexpected_gap(self):
        """Loops and structural gaps are expected; `gap` means a finding."""
        for sequence in (TRASTUZUMAB_VH, TRASTUZUMAB_VL):
            assert "gap" not in kinds(sequence)

    def test_a_short_cdr_loop_is_marked_as_a_loop_not_a_gap(self):
        """A kappa CDR2 fills 3 of 10 slots. That is normal, and the seven
        empty slots must not read as damage."""
        cdr2 = next(r for r in rows(TRASTUZUMAB_VL) if r[0] == "CDR2")
        _, low, high, cells = cdr2
        assert (low, high) == (56, 65)
        assert [kind for kind, _ in cells].count("residue") == 3
        assert [kind for kind, _ in cells].count("loop") == 7
        assert "gap" not in {kind for kind, _ in cells}

    def test_the_loop_fills_outward_from_both_ends(self):
        """IMGT leaves the MIDDLE empty rather than shifting residues along --
        which is what keeps position 23 meaning the same thing everywhere."""
        _, low, _, cells = next(r for r in rows(TRASTUZUMAB_VL) if r[0] == "CDR2")
        filled = [low + i for i, (kind, _) in enumerate(cells) if kind == "residue"]
        assert filled == [56, 57, 65]

    def test_structural_framework_slots_are_marked_structural(self):
        """Kappa is missing IMGT 73, 81 and 82 in nearly every chain, so FR3
        shows 36 of 39 on a chain the pipeline reports as clean."""
        domain = number_sequence(TRASTUZUMAB_VL)[0]
        expected = STRUCTURAL_FRAMEWORK_GAPS[domain.chain_type]
        found = {
            low + i
            for _, low, _, cells in _rows(domain)
            for i, (kind, _) in enumerate(cells)
            if kind == "structural"
        }
        assert found == expected

    def test_slots_beyond_the_numbered_span_are_outside_not_deleted(self):
        """The bug this project already shipped once: a read that stopped early
        reported every missing FR4 slot as a deletion."""
        _, low, _, cells = next(r for r in rows(TRASTUZUMAB_VL) if r[0] == "FR4")
        outside = [low + i for i, (kind, _) in enumerate(cells) if kind == "outside"]
        assert outside == [128]
        assert "gap" not in {kind for kind, _ in cells}

    def test_a_real_framework_deletion_is_marked_as_a_gap(self):
        """And it lands where the pipeline's own warning says it does."""
        domain = number_sequence(TRASTUZUMAB_VH)[0]
        mutant = number_sequence(domain.seq[:75] + domain.seq[76:])[0]
        flagged = {
            low + i
            for _, low, _, cells in _rows(mutant)
            for i, (kind, _) in enumerate(cells)
            if kind == "gap"
        }
        assert flagged, "expected an unexpected-gap cell"
        assert any(str(p) in w for p in flagged for w in mutant.warnings if "deletion" in w)


class TestTheRenderedFragment:
    @pytest.fixture
    def html(self):
        return render_track_html(number_sequence(TRASTUZUMAB_VL)[0])

    def test_every_region_is_labelled_so_colour_is_never_the_only_cue(self, html):
        for region in ("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"):
            assert f">{region}<" in html

    def test_conserved_anchors_are_marked(self, html):
        assert "ab-anchor" in html
        for anchor in ANCHOR_POSITIONS:
            assert f"IMGT {anchor}" in html

    def test_dark_mode_is_declared_under_both_scopes(self, html):
        """The media query covers the OS setting; the data-theme scope covers a
        theme toggle, and the toggle has to win either way."""
        assert "prefers-color-scheme: dark" in html
        assert ':root[data-theme="dark"]' in html

    def test_the_legend_distinguishes_the_one_mark_that_is_a_finding(self, html):
        assert "unexpected" in html
        assert "no deletion is reported" in html

    def test_residues_are_escaped_rather_than_interpolated_raw(self):
        """Sequence data reaches HTML, so it goes through escaping."""
        html = render_track_html(number_sequence(TRASTUZUMAB_VH)[0])
        assert "<script" not in html.lower()
        assert "title=" in html
