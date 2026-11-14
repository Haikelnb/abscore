"""Tests for the IEDB client's response handling.

The theme is that **a wrong number must never look like a good result.** A
dropped prediction row surfaces as a smaller count, which
reads exactly like a more human antibody. So the parser refuses a damaged
response rather than returning what it could salvage, and the cache -- written
only after a successful parse -- never stores one.

No network: these build responses by hand and read the committed fixtures.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from abscore.iedb import (
    EmptyResponseError,
    IEDBError,
    MalformedResponseError,
    _parse_tsv,
)

HEADER = "allele\tseq_num\tstart\tend\tpeptide\tcore_peptide\tscore\trank\n"
GOOD = "HLA-DRB1*01:01\t1\t1\t15\tAAAAAAAAAAAAAAA\tAAAAAAAAA\t0.5\t2.0"
SECOND = "HLA-DRB1*01:01\t1\t2\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\t1.0"

COMMITTED = sorted(Path("tests/data/iedb").glob("*.tsv"))


class TestAWellFormedResponse:
    def test_every_row_becomes_a_binding(self):
        bindings = _parse_tsv(HEADER + GOOD + "\n" + SECOND + "\n")
        assert len(bindings) == 2
        assert bindings[0].allele == "HLA-DRB1*01:01"
        assert bindings[0].peptide == "AAAAAAAAAAAAAAA"
        assert bindings[0].rank == 2.0

    def test_blank_lines_are_not_damage(self):
        bindings = _parse_tsv(HEADER + GOOD + "\n\n" + SECOND + "\n")
        assert len(bindings) == 2

    @pytest.mark.parametrize("path", COMMITTED, ids=lambda p: p.name[:12])
    def test_every_committed_fixture_still_parses(self, path):
        """These are real IEDB responses. If strictness rejects one of them,
        the rule is wrong -- not the fixture."""
        assert _parse_tsv(path.read_text())


class TestADamagedResponseIsRefused:
    """Each of these used to return one binding instead of two, silently."""

    @pytest.mark.parametrize(
        ("label", "row"),
        [
            ("missing rank", "HLA-DRB1*01:01\t1\t2\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\t"),
            ("unparseable rank", "HLA-DRB1*01:01\t1\t2\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\tabc"),
            ("missing allele", "\t1\t2\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\t1.0"),
            ("missing peptide", "HLA-DRB1*01:01\t1\t2\t16\t\tBBBBBBBBB\t0.4\t1.0"),
            ("truncated row", "HLA-DRB1*01:01\t1\t2"),
            ("unparseable start", "HLA-DRB1*01:01\t1\tx\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\t1.0"),
        ],
    )
    def test_one_bad_row_fails_the_whole_response(self, label, row):
        with pytest.raises(MalformedResponseError):
            _parse_tsv(HEADER + GOOD + "\n" + row + "\n")

    def test_the_error_names_the_line_and_what_was_missing(self):
        bad = "HLA-DRB1*01:01\t1\t2\t16\tBBBBBBBBBBBBBBB\tBBBBBBBBB\t0.4\t"
        with pytest.raises(MalformedResponseError, match=r"line 3 is missing rank"):
            _parse_tsv(HEADER + GOOD + "\n" + bad + "\n")

    def test_it_is_not_retried(self):
        """EmptyResponseError is transient and retried; this is not. The
        distinction is what stops a changed endpoint from being hammered."""
        assert issubclass(MalformedResponseError, IEDBError)
        assert not issubclass(MalformedResponseError, EmptyResponseError)


class TestAnEmptyResponse:
    def test_no_rows_is_the_transient_failure(self):
        with pytest.raises(EmptyResponseError):
            _parse_tsv(HEADER)

    def test_an_empty_body_is_the_transient_failure(self):
        with pytest.raises(EmptyResponseError):
            _parse_tsv("")

    def test_rows_that_are_all_unusable_is_not_an_empty_response(self):
        """Content that cannot be read points to a server change."""
        with pytest.raises(MalformedResponseError):
            _parse_tsv(HEADER + "\t\t\t\t\t\t\tx\n")


class TestColumnNameDrift:
    """IEDB has renamed columns between versions, so each field is looked up by
    any of its known names."""

    def test_percentile_rank_is_accepted_for_rank(self):
        header = "allele\tstart\tend\tsequence\tcore\tscore\tpercentile_rank\n"
        row = "HLA-DRB1*01:01\t1\t15\tAAAAAAAAAAAAAAA\tAAAAAAAAA\t0.5\t2.0"
        bindings = _parse_tsv(header + row + "\n")
        assert bindings[0].rank == 2.0
        assert bindings[0].peptide == "AAAAAAAAAAAAAAA"

    def test_a_response_with_no_rank_column_at_all_is_refused(self):
        header = "allele\tstart\tend\tpeptide\tcore_peptide\tscore\n"
        row = "HLA-DRB1*01:01\t1\t15\tAAAAAAAAAAAAAAA\tAAAAAAAAA\t0.5"
        with pytest.raises(MalformedResponseError, match="rank"):
            _parse_tsv(header + row + "\n")
