"""Tests for germline V assignment and identity.

No network: ANARCI's germline set is bundled with the package.
"""

import pytest

from abscore.germline import (
    SCOPES,
    GermlineError,
    assign_closest_germline,
    compare_to_germline,
    germline_data_sha256,
    germline_positions,
    germline_provenance,
    v_germlines,
)
from abscore.number import CDR_REGIONS, number_single_chain
from tests.fixtures import (
    ABAGOVOMAB_VH,
    ABELACIMAB_VL_LAMBDA,
    ADALIMUMAB_VH,
    RITUXIMAB_VH,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
)

OAS_PARTIAL_VH = (
    "SETLSLTCTVSGGSMSSYYWAWIRQPAGKRLEWIGRIYSSGAPTTIPPSRVESPISIDMSKSQFSLNL"
    "NSVTAADTAVYYCARDPAWFDPWGRGILVTVSS"
)


def assign(sequence):
    return assign_closest_germline(number_single_chain(sequence))


class TestGermlineReference:
    def test_every_allele_is_a_128_slot_imgt_string(self):
        """What makes identity a slot comparison rather than an alignment."""
        for chain in ("H", "K", "L"):
            for sequence in v_germlines(chain).values():
                assert len(sequence) == 128

    def test_human_sets_are_non_empty_for_all_three_chain_types(self):
        assert len(v_germlines("H")) > 100
        assert len(v_germlines("K")) > 50
        assert len(v_germlines("L")) > 50

    def test_v_alleles_stop_before_fr4(self):
        """A V gene does not encode FR4 -- that is the J gene.

        Human V alleles reach IMGT 106-108; J alleles start at 115. So positions
        107-114 are encoded by neither, being junctional.
        """
        for chain in ("H", "K", "L"):
            for sequence in v_germlines(chain).values():
                last = max(germline_positions(sequence))
                assert last <= 108, f"V allele reaches IMGT {last}"

    def test_germline_positions_drops_gaps(self):
        sequence = v_germlines("H")["IGHV4-4*07"]
        positions = germline_positions(sequence)
        assert "-" not in positions.values()
        assert positions[104] == "C"  # 2nd-CYS

    def test_unknown_chain_type_raises(self):
        with pytest.raises(GermlineError, match="no germline V set"):
            v_germlines("Z")


class TestScopes:
    def test_framework_scope_excludes_fr4_because_fr4_is_j_encoded(self):
        assert "FR4" not in SCOPES["framework"]
        assert "FR4" not in SCOPES["V"]

    def test_cdr_scope_excludes_cdr3_because_cdr3_is_junctional(self):
        assert "CDR3" not in SCOPES["cdr"]
        assert "CDR3" not in SCOPES["V"]
        assert set(SCOPES["cdr"]) == {"CDR1", "CDR2"}

    def test_v_scope_is_framework_plus_cdr1_and_cdr2(self):
        assert set(SCOPES["V"]) == set(SCOPES["framework"]) | set(SCOPES["cdr"])


class TestProvenance:
    def test_germline_data_hash_is_stable_within_a_run(self):
        assert germline_data_sha256() == germline_data_sha256()

    def test_provenance_pins_the_germline_data_and_admits_the_gap(self):
        provenance = germline_provenance()
        assert provenance["germline_data_matches_pin"] is True
        assert provenance["anarci_version"]
        # ANARCI records no IMGT release; the report must say so rather than guess.
        assert provenance["imgt_release"] == "not recorded by ANARCI"


class TestIdentityOrdering:
    """Identity must track the engineering history of each antibody."""

    def test_human_scores_above_humanized_above_murine(self):
        human = assign(ADALIMUMAB_VH).identity
        humanized = assign(TRASTUZUMAB_VH).identity
        murine = assign(ABAGOVOMAB_VH).identity
        assert human > humanized > murine

    def test_humanized_framework_is_human_while_its_cdrs_are_not(self):
        """The signature of CDR grafting: human scaffold, murine loops."""
        result = assign(TRASTUZUMAB_VH)
        assert result.scopes["framework"].identity > 0.85
        assert result.scopes["cdr"].identity < 0.60

    def test_framework_scope_separates_humanized_from_chimeric_better_than_v(self):
        """Why framework-only identity is worth reporting separately.

        A chimeric antibody's V-domain is entirely murine, so its framework is
        murine too; a humanized one has a human framework. Averaging the CDRs in
        shrinks that gap.
        """
        humanized = assign(TRASTUZUMAB_VH)
        chimeric = assign(RITUXIMAB_VH)
        v_gap = humanized.scopes["V"].identity - chimeric.scopes["V"].identity
        fr_gap = humanized.scopes["framework"].identity - chimeric.scopes["framework"].identity
        assert fr_gap > v_gap

    def test_chimeric_framework_is_not_human(self):
        assert assign(RITUXIMAB_VH).scopes["framework"].identity < 0.80


class TestGeneAssignment:
    @pytest.mark.parametrize(
        ("sequence", "prefix"),
        [
            (TRASTUZUMAB_VH, "IGHV"),
            (TRASTUZUMAB_VL, "IGKV"),
            (ABELACIMAB_VL_LAMBDA, "IGLV"),
        ],
    )
    def test_gene_matches_the_chain_type(self, sequence, prefix):
        assert assign(sequence).gene.startswith(prefix)

    def test_allele_level_ties_alone_do_not_make_a_call_ambiguous(self):
        """Alleles of one gene differ by a residue or two; the gene is still right."""
        result = assign(ADALIMUMAB_VH)
        assert [c.gene for c in result.contenders] == [result.gene_without_allele]
        assert not result.is_ambiguous
        assert result.ambiguous_genes == ()

    def test_an_exact_tie_between_different_genes_is_reported_as_ambiguous(self):
        """Rituximab's heavy chain fits IGHV1-46 and IGHV1-3 equally well.

        Nothing in the sequence separates them, which is what the caveat is for.
        """
        result = assign(RITUXIMAB_VH)
        assert result.is_ambiguous
        assert [c.gene for c in result.ambiguous_genes] == ["IGHV1-3"]
        assert result.ambiguous_genes[0].positions_behind == pytest.approx(0.0)
        assert any("ambiguous" in w for w in result.warnings())

    def test_a_gene_one_position_behind_is_not_ambiguous(self):
        """Trastuzumab: IGHV3-53 trails IGHV3-66 by one position.

        79 matches beats 78, so the call is not in doubt and must not be
        reported as if it were. An earlier threshold of 2 positions flagged
        this, which conflated "cannot be separated" with "separated by one
        residue" and crowded out the caveats that matter. The margin is still
        available, as `nearest_other_gene`.
        """
        result = assign(TRASTUZUMAB_VH)
        assert not result.is_ambiguous
        assert result.ambiguous_genes == ()
        assert result.nearest_other_gene is not None
        assert result.nearest_other_gene.gene == "IGHV3-53"
        assert result.nearest_other_gene.positions_behind == pytest.approx(1.0)

    def test_the_nearest_other_gene_is_found_even_when_an_allele_masks_it(self):
        """The bug this replaced: looking only at the next-best allele.

        Trastuzumab's runner-up allele is another IGHV3-66 at a dead tie, so a
        check on the next-best allele alone never sees IGHV3-53 at all. Alleles
        are collapsed to genes first, which is why the margin is visible.
        """
        result = assign(TRASTUZUMAB_VH)
        assert result.runner_up_is_same_gene  # the masking condition
        assert result.nearest_other_gene.gene == "IGHV3-53"

    def test_a_one_position_gap_is_not_lost_to_floating_point(self):
        """`positions_behind` is a float from an identity difference, and a gene
        exactly one position behind came out as 0.9999999999999964 -- which
        slipped under a threshold of 1.0, so "ties only" still admitted
        one-position gaps. Abagovomab was the case that showed it.
        """
        result = assign(ABAGOVOMAB_VH)
        assert not result.is_ambiguous
        assert result.nearest_other_gene.positions_behind == pytest.approx(1.0)

    def test_contenders_always_include_the_winning_gene(self):
        for sequence in (TRASTUZUMAB_VH, ADALIMUMAB_VH, ABELACIMAB_VL_LAMBDA):
            result = assign(sequence)
            assert result.gene_without_allele in {c.gene for c in result.contenders}

    def test_contenders_are_ordered_best_first(self):
        contenders = assign(ABAGOVOMAB_VH).contenders
        assert [c.identity for c in contenders] == sorted(
            (c.identity for c in contenders), reverse=True
        )

    def test_ambiguity_threshold_is_in_positions_not_percent(self):
        """A looser threshold admits more genes; a strict one admits only the winner."""
        from abscore.germline import assign_closest_germline as assign_with

        domain = number_single_chain(TRASTUZUMAB_VH)
        strict = assign_with(domain, ambiguity_positions=0.5)
        loose = assign_with(domain, ambiguity_positions=6.0)
        assert not strict.is_ambiguous
        assert len(loose.contenders) > len(strict.contenders)

    def test_every_contender_reports_how_far_behind_it_is(self):
        result = assign(ABAGOVOMAB_VH)
        assert all(0.0 <= c.positions_behind < 2.0 for c in result.contenders)
        winner = next(c for c in result.contenders if c.gene == result.gene_without_allele)
        assert winner.positions_behind == pytest.approx(0.0, abs=1e-9)

    def test_all_alleles_of_the_chain_type_are_considered(self):
        result = assign(TRASTUZUMAB_VH)
        assert result.n_alleles_considered == len(v_germlines("H"))

    def test_naming_a_gene_directly_gives_the_same_scores(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        best = assign_closest_germline(domain)
        direct = compare_to_germline(domain, best.gene)
        assert direct.scopes["V"] == best.scopes["V"]

    def test_unknown_gene_raises(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        with pytest.raises(GermlineError, match="not in the human H V set"):
            compare_to_germline(domain, "IGHV999-9*01")


class TestDenominators:
    def test_complete_domain_collapses_both_denominators(self):
        scope = assign(TRASTUZUMAB_VH).scopes["V"]
        assert scope.is_full_coverage
        assert scope.coverage == 1.0
        assert scope.identity == scope.identity_over_germline

    def test_partial_domain_separates_the_two_denominators(self):
        """The case the choice of denominator was made for."""
        scope = assign(OAS_PARTIAL_VH).scopes["V"]
        assert not scope.is_full_coverage
        assert scope.coverage < 1.0
        # Identity over positions present is HIGHER than over the germline span,
        # because absent positions count as mismatches in the latter.
        assert scope.identity > scope.identity_over_germline

    def test_partial_domain_carries_a_coverage_warning(self):
        result = assign(OAS_PARTIAL_VH)
        assert any("coverage" in w for w in result.warnings())

    def test_describe_states_the_denominator(self):
        """A bare percentage is not reportable; the denominator travels with it."""
        text = assign(TRASTUZUMAB_VH).scopes["V"].describe()
        assert "/" in text and "positions" in text

    def test_describe_mentions_coverage_only_when_incomplete(self):
        assert "coverage" not in assign(TRASTUZUMAB_VH).scopes["V"].describe()
        assert "coverage" in assign(OAS_PARTIAL_VH).scopes["V"].describe()


class TestMismatches:
    def test_mismatch_count_is_consistent_with_identity(self):
        result = assign(TRASTUZUMAB_VH)
        scope = result.scopes["V"]
        present = [m for m in result.mismatches if not m.is_absence]
        assert len(present) == scope.compared - scope.matches

    def test_mismatches_split_into_framework_and_cdr(self):
        result = assign(TRASTUZUMAB_VH)
        assert len(result.mismatches_in_framework()) + len(result.mismatches_in_cdrs()) == len(
            result.mismatches
        )

    def test_humanized_mismatches_concentrate_in_cdrs_relative_to_length(self):
        """8 CDR positions out of ~15 vs 8 framework out of ~80."""
        result = assign(TRASTUZUMAB_VH)
        cdr_rate = len(result.mismatches_in_cdrs()) / result.scopes["cdr"].compared
        fr_rate = len(result.mismatches_in_framework()) / result.scopes["framework"].compared
        assert cdr_rate > fr_rate

    def test_absent_positions_are_recorded_as_absences_not_substitutions(self):
        result = assign(OAS_PARTIAL_VH)
        absences = [m for m in result.mismatches if m.is_absence]
        assert absences
        assert all(m.query_residue is None for m in absences)

    def test_every_mismatch_carries_a_region(self):
        result = assign(TRASTUZUMAB_VH)
        valid = set(CDR_REGIONS) | {"FR1", "FR2", "FR3"}
        assert all(m.region in valid for m in result.mismatches)
