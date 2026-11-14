"""Tests for repertoire-based (OASis-style) humanness.

No network: `promb` bundles its peptide table with the package.
"""

import pytest

from abscore.humanness import (
    DATABASES,
    OAS_THRESHOLD_MIN_DONOR_FRACTION,
    OAS_THRESHOLD_NAME,
    HumannessError,
    get_database,
    score_domain,
    score_sequence,
)
from abscore.number import CDR_REGIONS, IMGT_REGIONS, number_single_chain
from tests.fixtures import (
    ABAGOVOMAB_VH,
    ABELACIMAB_VL_LAMBDA,
    ADALIMUMAB_VH,
    ELIPOVIMAB_VH,
    PEMBROLIZUMAB_VH,
    RITUXIMAB_VH,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
)

OAS_PARTIAL_VH = (
    "SETLSLTCTVSGGSMSSYYWAWIRQPAGKRLEWIGRIYSSGAPTTIPPSRVESPISIDMSKSQFSLNL"
    "NSVTAADTAVYYCARDPAWFDPWGRGILVTVSS"
)


class TestDatabase:
    def test_oasis_uses_9mers(self):
        assert get_database("human-oas").peptide_length == 9

    def test_database_is_cached(self):
        assert get_database("human-oas") is get_database("human-oas")

    def test_unknown_database_raises(self):
        with pytest.raises(HumannessError, match="unknown database"):
            get_database("not-a-database")

    def test_proteome_databases_are_available_for_the_self_filter(self):
        """epitopes.py needs a broader "self" than Ig germline alone."""
        assert "human-swissprot" in DATABASES
        assert "human-reference" in DATABASES


class TestAgreementWithPromb:
    """Our window-chopping must reproduce promb's own score exactly.

    We chop windows ourselves rather than calling promb's scorer, so that each
    peptide can carry its IMGT positions. This pins that the extra bookkeeping
    does not change the number.
    """

    @pytest.mark.parametrize(
        "sequence",
        [
            ADALIMUMAB_VH,
            TRASTUZUMAB_VH,
            TRASTUZUMAB_VL,
            RITUXIMAB_VH,
            ABAGOVOMAB_VH,
            ELIPOVIMAB_VH,
            ABELACIMAB_VL_LAMBDA,
        ],
    )
    def test_score_matches_promb_exactly(self, sequence):
        domain = number_single_chain(sequence)
        assert score_domain(domain).score == pytest.approx(
            get_database().compute_peptide_content(domain.seq), abs=1e-12
        )


class TestScoreOrdering:
    def test_human_above_humanized_above_chimeric_above_murine(self):
        human = score_sequence(ADALIMUMAB_VH).score
        humanized = score_sequence(TRASTUZUMAB_VH).score
        chimeric = score_sequence(RITUXIMAB_VH).score
        murine = score_sequence(ABAGOVOMAB_VH).score
        assert human > humanized > chimeric > murine

    def test_score_is_a_fraction(self):
        result = score_sequence(TRASTUZUMAB_VH)
        assert 0.0 <= result.score <= 1.0
        assert result.score == result.n_found / result.n_peptides


class TestWindowing:
    def test_peptide_count_follows_the_window(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        result = score_domain(domain)
        assert result.n_peptides == len(domain.seq) - result.peptide_length + 1

    def test_every_peptide_is_the_right_length(self):
        result = score_sequence(TRASTUZUMAB_VH)
        assert all(len(p.sequence) == result.peptide_length for p in result.peptides)

    def test_every_peptide_carries_one_position_per_residue(self):
        result = score_sequence(TRASTUZUMAB_VH)
        assert all(len(p.positions) == len(p.sequence) for p in result.peptides)

    def test_windows_are_chopped_in_imgt_order_not_input_order(self):
        """Matters for CDR3 insertions, where the two orders differ."""
        domain = number_single_chain(ELIPOVIMAB_VH)
        result = score_domain(domain)
        rebuilt = result.peptides[0].sequence + "".join(p.sequence[-1] for p in result.peptides[1:])
        assert rebuilt == domain.seq


class TestRegionLocalisation:
    def test_non_human_peptides_concentrate_in_cdrs(self):
        """The prediction, asserted rather than eyeballed."""
        counts = score_sequence(TRASTUZUMAB_VH).counts_by_region()
        cdr_rate = sum(counts[r][0] for r in CDR_REGIONS) / sum(counts[r][1] for r in CDR_REGIONS)
        fr = ("FR1", "FR2", "FR3", "FR4")
        fr_rate = sum(counts[r][0] for r in fr) / sum(counts[r][1] for r in fr)
        assert cdr_rate > fr_rate

    def test_region_counts_are_not_a_partition(self):
        """A peptide spanning two regions is counted in both, by design."""
        result = score_sequence(TRASTUZUMAB_VH)
        total = sum(t for _, t in result.counts_by_region().values())
        assert total > result.n_peptides

    def test_boundary_spanning_peptides_straddle_a_framework_cdr_boundary(self):
        result = score_sequence(TRASTUZUMAB_VH)
        assert result.boundary_spanning_peptides
        for peptide in result.boundary_spanning_peptides:
            assert peptide.touches_cdr
            assert any(r not in CDR_REGIONS for r in peptide.regions)
            assert not peptide.found

    def test_a_peptide_inside_one_region_does_not_span_a_boundary(self):
        result = score_sequence(TRASTUZUMAB_VH)
        interior = [p for p in result.peptides if len(p.regions) == 1]
        assert interior
        assert not any(p.spans_boundary for p in interior)


class TestPerResidueFlags:
    def test_flag_count_is_bounded_by_the_window_length(self):
        result = score_sequence(TRASTUZUMAB_VH)
        assert all(1 <= n <= result.peptide_length for n in result.flags_by_position().values())

    def test_only_flagged_positions_appear(self):
        result = score_sequence(ADALIMUMAB_VH)
        covered = {p for pep in result.non_human_peptides for p in pep.positions}
        assert set(result.flags_by_position()) == covered

    def test_worst_positions_are_ranked_by_flag_count(self):
        worst = score_sequence(TRASTUZUMAB_VH).worst_positions(8)
        assert [c for _, c, _ in worst] == sorted((c for _, c, _ in worst), reverse=True)

    def test_a_fully_human_region_has_no_flags(self):
        result = score_sequence(ELIPOVIMAB_VH)
        assert result.counts_by_region()["FR1"][0] == 0
        flagged_fr1 = [p for p in result.flags_by_position() if 1 <= p[0] <= 26]
        assert flagged_fr1 == []


class TestKnownLimitations:
    """Behaviour that is a real weakness of the method, pinned so the README
    can describe it accurately rather than from memory."""

    def test_long_cdr3_makes_a_human_antibody_indistinguishable_from_a_chimeric(self):
        """Elipovimab is labelled Human yet scores level with a chimeric antibody.

        CDR3 is junctional, so it is near-unique to one individual and its
        9-mers are mostly absent from other donors' repertoires. A long CDR3
        therefore reads as non-human. Germline identity cannot see this at all,
        because a V gene does not encode CDR3.

        The two scores sit within a fraction of a point of each other, so the
        assertion is indistinguishability rather than an ordering -- which is
        itself the failure: the score cannot tell a fully human antibody from one
        with entirely murine variable domains.
        """
        elipovimab = score_sequence(ELIPOVIMAB_VH)
        chimeric = score_sequence(RITUXIMAB_VH)
        assert abs(elipovimab.score - chimeric.score) < 0.02

        # And it is CDR3 that does it.
        bad, total = elipovimab.counts_by_region()["CDR3"]
        assert bad / total > 0.8

    def test_excluding_cdr3_lifts_the_long_cdr3_chain_most(self):
        def without_cdr3(sequence):
            result = score_sequence(sequence)
            kept = [p for p in result.peptides if "CDR3" not in p.regions]
            return sum(1 for p in kept if p.found) / len(kept) - result.score

        assert without_cdr3(ELIPOVIMAB_VH) > without_cdr3(TRASTUZUMAB_VH)

    def test_homopolymer_scores_as_human(self):
        """promb's own warning: "more human doesn't always mean good"."""
        db = get_database()
        assert db.contains("GGGGGGGGG")


class TestProvenanceAndWarnings:
    def test_provenance_names_the_threshold_and_its_limit(self):
        provenance = score_sequence(TRASTUZUMAB_VH).provenance()
        assert provenance["threshold_name"] == OAS_THRESHOLD_NAME == "relaxed"
        assert provenance["threshold_min_donor_fraction"] == OAS_THRESHOLD_MIN_DONOR_FRACTION
        # The published ADA correlation used the medium threshold instead.
        assert "medium" in provenance["note"]

    def test_provenance_records_versions_and_source(self):
        provenance = score_sequence(TRASTUZUMAB_VH).provenance()
        assert provenance["promb_version"]
        assert provenance["peptide_length"] == 9
        assert "2019" in provenance["oas_snapshot"]
        assert "promb asserts nothing" in provenance["oas_snapshot"]
        assert "10perc_subjects" in provenance["reference_file"]

    def test_complete_chain_has_no_warnings(self):
        assert score_sequence(TRASTUZUMAB_VH).warnings() == ()

    def test_partial_domain_warns_that_the_score_is_not_comparable(self):
        result = score_sequence(OAS_PARTIAL_VH)
        assert any("not directly comparable" in w for w in result.warnings())

    def test_pembrolizumab_scores_as_human_despite_numbering_as_mouse(self):
        """Separates two different questions that are easy to conflate.

        ANARCI's species call is a numbering nearest-neighbour; the humanness
        score is a repertoire lookup. Pembrolizumab numbers best against the
        mouse HMM while scoring in the humanized range.
        """
        domain = number_single_chain(PEMBROLIZUMAB_VH)
        assert domain.species == "mouse"
        assert score_domain(domain).score > score_sequence(ABAGOVOMAB_VH).score


class TestFlaggedRuns:
    """Contiguous flagged stretches, reported as observation not diagnosis.

    An earlier version of this module called every boundary-spanning non-human
    peptide a "junction peptide", implying a graft seam. On a real OAS record
    with a compensated indel that labelled a 14-peptide frameshift artefact as 14
    junctions. Position does not establish cause, so the API reports runs and
    leaves the interpretation to the reader.
    """

    def test_a_clean_region_produces_no_run(self):
        runs = score_sequence(ELIPOVIMAB_VH).flagged_runs()
        assert all(not (start[0] <= 26 and end[0] <= 26) for start, end, _, _ in runs)

    def test_runs_are_contiguous_and_ordered(self):
        runs = score_sequence(TRASTUZUMAB_VH).flagged_runs()
        assert runs
        for start, end, length, _ in runs:
            assert start[0] <= end[0]
            assert length >= 1
        starts = [start[0] for start, _, _, _ in runs]
        assert starts == sorted(starts)

    def test_runs_cover_every_flagged_position_exactly_once(self):
        result = score_sequence(TRASTUZUMAB_VH)
        total = sum(length for _, _, length, _ in result.flagged_runs())
        assert total == len(result.flags_by_position())

    def test_each_run_names_the_regions_it_spans(self):
        for _, _, _, regions in score_sequence(TRASTUZUMAB_VH).flagged_runs():
            assert regions
            assert all(r in IMGT_REGIONS for r in regions)

    def test_a_frameshifted_record_gives_one_long_run_not_many_junctions(self):
        """The case that motivated dropping the junction label."""
        result = score_sequence(OAS_PARTIAL_VH)
        longest = max(result.flagged_runs(), key=lambda run: run[2])
        start, end, length, regions = longest
        assert length >= 15
        assert "CDR2" in regions and "FR3" in regions
        # Many boundary-spanning peptides, yet they form one stretch.
        assert len(result.boundary_spanning_peptides) > 10

    def test_spans_boundary_is_geometry_only(self):
        """Documented as a fact about position, with nothing reading cause from it."""
        import inspect

        from abscore.humanness import Peptide

        doc = inspect.getdoc(Peptide.spans_boundary.fget)
        assert "A geometric fact" in doc


class TestALowComplexityInsertRaisesTheScore:
    """A limitation of the published definition, pinned rather than patched.

    OASis is the fraction of a chain's 9-mer windows found in the reference
    peptide set. A run of one residue contributes many identical windows, so a
    long glycine insert adds windows that all match the same reference peptide
    and the fraction climbs. Inserting 100 glycines into a heavy-chain CDR3
    takes OASis from roughly 70% to roughly 81% while germline V identity does
    not move, the insert being in CDR3.

    The score is deliberately left alone. It is the published OASis definition,
    and quietly adding a complexity penalty would make the number incomparable
    to the literature it is taken from. What this test does is keep the
    behaviour documented and honest: it is recorded in
    `provenance.KNOWN_SENSITIVITIES["low_complexity_inflates_humanness"]`, and
    if the arithmetic ever changes, this fails and that entry has to be
    rewritten.

    It says nothing about whether such a domain would fold, bind, or be
    tolerated. A high OASis on a synthetic repeat is a statement about peptide
    bookkeeping, not about an antibody.
    """

    @staticmethod
    def _with_glycines(count):
        """Insert into CDR3, which germline V identity does not cover."""
        cdr3_start = TRASTUZUMAB_VH.index("SRWGGDGFYAMDY")
        return TRASTUZUMAB_VH[: cdr3_start + 2] + "G" * count + TRASTUZUMAB_VH[cdr3_start + 2 :]

    def test_the_repeated_window_is_in_the_reference_set(self):
        """The mechanism, not just the symptom: every added window is this one
        peptide, and it is genuinely present in the human reference."""
        from abscore.humanness import get_database

        assert get_database().contains("GGGGGGGGG")

    def test_the_score_climbs_with_the_length_of_the_repeat(self):
        baseline = score_domain(number_single_chain(TRASTUZUMAB_VH)).score
        inflated = score_domain(number_single_chain(self._with_glycines(100))).score
        assert inflated > baseline + 0.05, (baseline, inflated)

    def test_germline_identity_is_untouched_by_it(self):
        """Which is why the two scores disagree on such a sequence, and why
        reporting both is worth the column."""
        from abscore.germline import assign_closest_germline

        baseline = assign_closest_germline(number_single_chain(TRASTUZUMAB_VH))
        inflated = assign_closest_germline(number_single_chain(self._with_glycines(100)))
        assert inflated.identity == baseline.identity

    def test_the_limitation_is_recorded_in_provenance(self):
        from abscore.provenance import KNOWN_SENSITIVITIES

        entry = KNOWN_SENSITIVITIES["low_complexity_inflates_humanness"]
        assert "does not establish" in entry["residual"]
