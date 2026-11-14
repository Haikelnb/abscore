"""Tests for MHC class II predicted-binder burden against a human reference.

No network: every IEDB response these need is committed under `tests/data/iedb`
and `offline=True` makes a cache miss an error rather than a silent fetch. See
that directory's README for provenance.
"""

import pytest

from abscore.epitopes import (
    BINDER_RANK,
    FAR,
    JUNCTIONAL_WINDOW,
    SELF_DEFINITIONS,
    UNMATCHED_THRESHOLD,
    Epitope,
    EpitopeBurden,
    epitope_burden_for_sequence,
    germline_peptide_source,
)
from abscore.iedb import DEFAULT_ALLELES, DEFAULT_METHOD, Binding, OfflineError
from abscore.number import CDR_REGIONS
from tests.fixtures import ABAGOVOMAB_VH, ADALIMUMAB_VH, RITUXIMAB_VH, TRASTUZUMAB_VH

CACHE = "tests/data/iedb"


def burden(sequence, **kwargs):
    return epitope_burden_for_sequence(sequence, cache_dir=CACHE, offline=True, **kwargs)


@pytest.fixture(scope="module")
def humanized():
    return burden(TRASTUZUMAB_VH)


@pytest.fixture(scope="module")
def human():
    return burden(ADALIMUMAB_VH)


@pytest.fixture(scope="module")
def murine():
    return burden(ABAGOVOMAB_VH)


class TestOfflineReplay:
    def test_cached_requests_need_no_network(self, humanized):
        assert humanized.n_predictions > 0

    def test_an_uncached_sequence_raises_rather_than_fetching(self):
        """Keeps the suite honest: a new fixture cannot quietly hit the network."""
        with pytest.raises(OfflineError, match="no cached IEDB response"):
            burden(RITUXIMAB_VH)

    def test_replay_is_deterministic(self):
        assert burden(TRASTUZUMAB_VH).counts_by_definition() == (
            burden(TRASTUZUMAB_VH).counts_by_definition()
        )


class TestGermlineSelfPeptides:
    def test_the_self_peptide_set_is_built_from_v_and_j_segments(self):
        """FR4 is J-encoded, so a framework peptide there is no less self."""
        peptides = germline_peptide_source(9)
        assert len(peptides) > 5000
        # WGQGTLVTV sits in FR4, which only a J segment can supply.
        assert "WGQGTLVTV" in peptides

    def test_peptides_are_the_requested_length(self):
        assert all(len(p) == 9 for p in germline_peptide_source(9))

    def test_a_conserved_framework_3_motif_is_self(self):
        """The motif that makes naive counting wrong, asserted directly."""
        assert "YLQMNSLRA" in germline_peptide_source(9)


class TestSelfFiltering:
    def test_filtering_removes_binders(self, humanized):
        counts = humanized.counts_by_definition()
        assert counts["raw"] > counts["germline_unmatched"] > 0

    def test_the_strongest_raw_binder_is_a_germline_framework_peptide(self, humanized):
        """The worked example behind this module.

        The top predicted binder in trastuzumab across the eight-allele panel is
        a conserved framework 3 peptide that human IGHV3 germline genes encode.
        Counted naively it is the molecule's top risk; it is self.
        """
        top = humanized.binders[0]
        assert top.core == "YLQMNSLRA"
        assert top.regions == ("FR3",)
        assert top.matches_reference("germline")
        assert top not in humanized.unmatched("germline")

    def test_every_self_definition_is_evaluated_for_every_binder(self, humanized):
        for epitope in humanized.binders:
            assert set(epitope.distances) == set(SELF_DEFINITIONS)
            assert all(d in (0, 1, FAR) for d in epitope.distances.values())

    def test_repertoire_is_the_most_permissive_definition(self, humanized):
        """OAS contains somatically mutated peptides, so it excuses more as self."""
        counts = humanized.counts_by_definition()
        assert counts["repertoire_unmatched"] < counts["germline_unmatched"]

    def test_outside_the_junctional_window_proteome_matches_germline(self, humanized):
        """Ig framework peptides are in the proteome because they are expressed.

        So for an antibody the proteome filter catches nothing the germline
        filter misses -- except in the junctional window, where the germline
        definition abstains and the proteome definition still applies.
        """
        germline = set(humanized.unique_cores("germline"))
        proteome = set(humanized.unique_cores("proteome"))
        junctional = {e.core for e in humanized.unclassifiable("germline")}
        assert germline == proteome - junctional


class TestSeparationByOrigin:
    """How the counts behave across three fixed sequences of known origin.

    These are regression checks on four committed responses, and nothing more.
    Species origin and MHC II presentation are different endpoints: a fully
    human sequence can carry many predicted hits, so agreement with origin here
    is not evidence that either count is the correct immunogenicity ordering.
    """

    def test_the_raw_hit_count_does_not_track_origin(self, humanized, murine):
        """The murine antibody yields FEWER raw predicted hits than the
        humanized one, so the raw count alone orders them the other way."""
        assert murine.count() < humanized.count()

    def test_the_unmatched_count_separates_them_by_origin(self, human, humanized, murine):
        assert murine.count("germline") > humanized.count("germline") > human.count("germline")

    def test_no_hit_in_the_murine_antibody_matches_human_germline(self, murine):
        """Its V-domain is entirely murine, so nothing matches human germline.

        Every hit is therefore either unmatched or, in the junctional window,
        unclassifiable -- never matched.
        """
        counts = murine.counts_by_definition()
        assert counts["germline_unmatched"] + counts["germline_unclassifiable"] == counts["raw"]
        assert not [e for e in murine.binders if e.matches_reference("germline")]


class TestLocalisation:
    def test_unmatched_burden_concentrates_in_cdrs(self, humanized):
        """The prediction from the tolerance argument, as a number.

        CDR-touching share of binders rises as the self filter tightens.
        """

        def cdr_share(pool):
            return sum(1 for e in pool if e.touches_cdr) / len(pool)

        raw = cdr_share(humanized.binders)
        germline = cdr_share(humanized.unmatched("germline"))
        repertoire = cdr_share(humanized.unmatched("repertoire"))
        assert raw < germline <= repertoire

    def test_filtering_removes_far_more_from_framework_than_from_cdrs(self, humanized):
        raw = humanized.counts_by_region(None)
        kept = humanized.counts_by_region("germline")
        assert raw["FR3"] - kept["FR3"] > 20
        # Framework loses far more than CDR2 does. CDR2 is not untouched --
        # a couple of somatically mutated cores sit one residue from germline
        # there too -- but the drop is small next to framework's.
        assert (raw["FR3"] - kept["FR3"]) > 5 * (raw["CDR2"] - kept["CDR2"])

    def test_every_epitope_is_located_in_imgt_coordinates(self, humanized):
        for epitope in humanized.binders:
            assert len(epitope.positions) == len(epitope.peptide)
            assert epitope.regions

    def test_the_core_can_be_located_inside_the_peptide(self, humanized):
        located = [e for e in humanized.binders if e.core_positions]
        assert located
        for epitope in located:
            assert len(epitope.core_positions) == len(epitope.core)


class TestCounting:
    def test_unique_cores_collapse_one_peptide_binding_many_alleles(self, humanized):
        """Two different questions: how many risky peptides, how many presentations."""
        assert len(humanized.unique_cores(None)) < humanized.count()

    def test_per_allele_counts_cover_the_whole_panel(self, humanized):
        counts = humanized.counts_by_allele("germline")
        assert set(counts) == set(DEFAULT_ALLELES)
        assert sum(counts.values()) == humanized.count("germline")

    def test_a_caller_chosen_rank_cutoff_subsets_the_binders(self, humanized):
        """No built-in "strong binder" tier: no rank cutoff could be sourced."""
        strict = humanized.binders_below_rank(2.0)
        assert set(strict) <= set(humanized.binders)
        assert all(e.rank <= 2.0 for e in strict)

    def test_all_binders_clear_the_rank_cutoff(self, humanized):
        assert all(e.rank <= BINDER_RANK for e in humanized.binders)

    def test_binders_are_a_subset_of_predictions(self, humanized):
        assert humanized.count() < humanized.n_predictions


class TestWeighting:
    def test_weighting_requires_caller_supplied_frequencies(self, humanized):
        """No frequency table is shipped; frequencies are population-specific."""
        frequencies = dict.fromkeys(DEFAULT_ALLELES, 0.1)
        weighted = humanized.weighted_burden(frequencies, "germline")
        assert weighted == pytest.approx(humanized.count("germline") * 0.1)

    def test_frequencies_for_alleles_outside_the_panel_are_rejected(self, humanized):
        with pytest.raises(KeyError, match="not in the panel"):
            humanized.weighted_burden({"HLA-DRB1*04:05": 0.2})

    def test_provenance_states_that_weighting_is_not_applied_by_default(self, humanized):
        note = humanized.provenance()["frequency_weighting"]
        assert "not applied by default" in note


class TestProvenance:
    def test_the_predictor_and_service_are_recorded(self, humanized):
        record = humanized.provenance()
        assert record["method"] == DEFAULT_METHOD
        assert "iedb" in str(record["url"]).lower()
        assert record["peptide_length"] == 15

    def test_the_allele_panel_and_its_source_are_recorded(self, humanized):
        record = humanized.provenance()
        assert len(record["alleles"]) == 8
        assert "Prihoda" in record["allele_panel_source"]
        # The widely cited nine-allele set is NOT what the paper used.
        assert "04:05" in record["allele_panel_source"]

    def test_the_licence_position_is_stated(self, humanized):
        note = humanized.provenance()["licence_note"]
        assert "never" in note and "downloaded" in note

    def test_the_missing_model_version_is_disclosed(self, humanized):
        assert "not reported" in humanized.provenance()["software_version"]

    def test_rank_rather_than_affinity_is_justified(self, humanized):
        assert "percentile rank" in humanized.provenance()["threshold_basis"]

    def test_each_self_definition_is_documented(self, humanized):
        definitions = humanized.provenance()["self_definitions"]
        assert set(definitions) == set(SELF_DEFINITIONS)
        assert all(definitions.values())

    def test_the_recorded_threshold_is_the_one_the_counts_used(self, humanized):
        """Provenance alone is not enough: an earlier version reported the
        stored threshold while every count still used the module default, so a
        report could advertise a threshold of 3 and return the answer for 2.
        This asserts membership; the text is checked elsewhere."""
        from dataclasses import replace

        for threshold in (1, 2, 3):
            altered = replace(humanized, mismatch_threshold=threshold)
            expected = sum(1 for e in altered.epitopes if e.is_unmatched("germline", threshold))
            assert altered.provenance()["unmatched_threshold"] == threshold
            assert altered.count("germline") == expected
            assert len(altered.unmatched("germline")) == expected
            assert altered.counts_by_definition()["germline_unmatched"] == expected
            assert sum(altered.counts_by_allele("germline").values()) == expected

    def test_a_stricter_threshold_can_only_shrink_the_count(self, humanized):
        from dataclasses import replace

        counts = [replace(humanized, mismatch_threshold=t).count("germline") for t in (1, 2, 3)]
        assert counts == sorted(counts, reverse=True), counts

    def test_an_explicit_threshold_still_overrides_for_a_sensitivity_analysis(self, humanized):
        """`threshold=1` is exact matching, and must stay available."""
        assert len(humanized.unmatched("germline", 1)) >= len(humanized.unmatched("germline"))

    def test_overridden_settings_are_recorded_rather_than_the_defaults(self, humanized):
        """Provenance that reports the default while the run used an override is
        worse than no provenance: it is a false record of how a number was made.
        Constructed by `replace` rather than by a second prediction, because the
        plumbing under test is what the burden reports, given its own fields."""
        from dataclasses import replace

        altered = replace(
            humanized,
            alleles=("HLA-DRB1*15:01",),
            peptide_length=13,
            mismatch_threshold=3,
        )
        record = altered.provenance()
        assert record["alleles"] == ["HLA-DRB1*15:01"]
        assert record["peptide_length"] == 13
        assert record["unmatched_threshold"] == 3
        assert record["alleles_are_default_panel"] is False
        # The default panel is still described, and still flagged as the default.
        assert "Prihoda" in record["allele_panel_source"]
        assert humanized.provenance()["alleles_are_default_panel"] is True


class TestRegionsAreNotAPartition:
    def test_a_peptide_spanning_two_regions_counts_in_both(self, humanized):
        counts = humanized.counts_by_region(None)
        assert sum(counts.values()) > humanized.count()

    def test_cdr_regions_are_recognised(self, humanized):
        touching = [e for e in humanized.binders if e.touches_cdr]
        assert touching
        assert all(any(r in CDR_REGIONS for r in e.regions) for e in touching)


class TestJunctionalWindow:
    """The germline definition of self has a boundary, and it must admit it.

    CDR3 spans the V-D-J junction. A germline D gene has no single protein
    sequence -- it is used in any of six reading frames and trimmed by a variable
    number of nucleotides at each end -- so no germline protein reference exists
    for the junctional window. Neither ANARCI nor IMGT's amino-acid reference
    ships D protein sequences, which is why.

    Calling those peptides "not germline-encoded, therefore unmatched" would
    overstate the burden exactly where the metric is most sensitive. On
    adalimumab it overstated it by 17 of 35.
    """

    def test_the_window_sits_between_v_and_j_coverage(self):
        """Derived from germline coverage at import time."""
        from anarci.germlines import all_germlines

        for chain, (low, high) in JUNCTIONAL_WINDOW.items():
            v_end = max(
                max(i + 1 for i, c in enumerate(s) if c != "-")
                for s in all_germlines["V"][chain]["human"].values()
            )
            j_start = min(
                min(i + 1 for i, c in enumerate(s) if c != "-")
                for s in all_germlines["J"][chain]["human"].values()
            )
            assert low == v_end + 1
            assert high == j_start - 1

    def test_junctional_peptides_are_not_counted_as_germline_unmatched(self, humanized):
        for epitope in humanized.unclassifiable("germline"):
            assert epitope.is_junctional
            assert not epitope.is_unmatched("germline")

    def test_junctional_peptides_overlap_cdr3(self, human):
        assert human.unclassifiable("germline")
        for epitope in human.unclassifiable("germline"):
            assert "CDR3" in epitope.regions

    def test_the_repertoire_definition_still_applies_in_the_window(self, human):
        """It is built from real rearranged sequences, so it can see CDR3."""
        junctional = human.unclassifiable("germline")
        assert any(
            e.matches_reference("repertoire") or e.is_unmatched("repertoire") for e in junctional
        )

    def test_abstaining_lowers_the_count_rather_than_raising_it(self, human):
        """The fix is conservative: it can only remove from the burden."""
        counts = human.counts_by_definition()
        assert counts["germline_unclassifiable"] > 0
        assert counts["germline_unmatched"] < counts["raw"]

    def test_the_ordering_by_origin_survives_the_abstention(self, human, humanized, murine):
        assert murine.count("germline") > humanized.count("germline") > human.count("germline")

    def test_the_window_is_recorded_in_provenance(self, humanized):
        record = humanized.provenance()
        assert record["junctional_window_imgt"]["H"] == [108, 114]
        assert "IC50" in record["binder_cutoff_source"]


class TestDistanceAnnotation:
    """Distance to the reference is recorded, and is not used as a discount.

    The headline count is everything with no exact match. How far the nearest
    reference peptide sits is kept alongside, because it is informative -- but
    one mismatch is not treated as safer than two, since whether a substitution
    matters depends on which position it hits, which this does not model.

    The chain below is why the annotation exists at all: a human lambda chain at
    95.4% germline identity, carrying four somatic mutations, produces 28
    unmatched hits whose cores are *every one* exactly one residue from a
    germline peptide. That is worth seeing. It is not grounds for reporting zero.
    """

    LAMBDA_FOUR_MUTATIONS = (
        "SYVLTQPPSVSVAPGQTARITCGGDNIGDKSVHWYQQKTGQAPVLVVYDDSDRPSGIPERFSGSNSGNTATL"
        "TITRVEAGDEADYYCQVWDITSDHVVFGGGTQLTVL"
    )

    def test_distances_are_zero_one_or_far(self, humanized):
        for epitope in humanized.binders:
            assert epitope.distance_to_reference("germline") in (0, 1, FAR)

    def test_an_exact_germline_peptide_matches_at_distance_zero(self, humanized):
        top = humanized.binders[0]
        assert top.core == "YLQMNSLRA"
        assert top.distance_to_reference("germline") == 0
        assert top.matches_reference("germline")
        assert not top.is_one_substitution("germline")

    def test_one_substitution_hits_are_counted_as_unmatched(self, humanized):
        """The change this class exists to pin: a near-reference core is
        unmatched, and stays counted."""
        near = set(humanized.one_substitution("germline"))
        assert near
        assert near <= set(humanized.unmatched("germline"))

    def test_the_subsets_sum_back_to_the_unmatched_count(self, humanized):
        counts = humanized.counts_by_definition()
        for definition in SELF_DEFINITIONS:
            assert (
                counts[f"{definition}_one_substitution"] + counts[f"{definition}_distant"]
                == counts[f"{definition}_unmatched"]
            )

    def test_the_three_categories_partition_the_predicted_hits(self, humanized):
        counts = humanized.counts_by_definition()
        for definition in SELF_DEFINITIONS:
            assert (
                counts[f"{definition}_matched"]
                + counts[f"{definition}_unmatched"]
                + counts[f"{definition}_unclassifiable"]
                == counts["raw"]
            )

    def test_the_default_threshold_counts_every_non_exact_match(self, humanized):
        assert humanized.mismatch_threshold == UNMATCHED_THRESHOLD == 1

    def test_the_exploratory_threshold_is_available_but_is_not_the_default(self, humanized):
        exploratory = humanized.distant("germline")
        headline = humanized.unmatched("germline")
        assert set(exploratory) < set(headline), "distant must be a strict subset here"
        assert humanized.count("germline") == len(headline)

    def test_a_distance_beyond_one_is_never_displayed_as_an_exact_two(self, humanized):
        labels = {e.distance_label("germline") for e in humanized.binders}
        assert labels <= {"0", "1", "2+", "n/a"}
        assert "2" not in labels


class TestOnRealChains:
    """Pinned against sequences where the distance annotation is the whole story."""

    @pytest.fixture(scope="class")
    @classmethod
    def lambda_chain(cls):
        return burden(TestDistanceAnnotation.LAMBDA_FOUR_MUTATIONS)

    def test_a_near_germline_human_chain_reports_its_hits_rather_than_hiding_them(
        self, lambda_chain
    ):
        """28 unmatched, all of them one substitution from germline.

        An earlier version reported 0 here, by excluding one-substitution
        variants. That was a hypothesis about tolerance doing the work of a
        measurement. The hits are now counted and their distance is shown.
        """
        counts = lambda_chain.counts_by_definition()
        assert counts["germline_unmatched"] == 28
        assert counts["germline_one_substitution"] == 28
        assert counts["germline_distant"] == 0

    def test_every_core_of_that_chain_is_one_residue_from_germline(self, lambda_chain):
        from abscore.epitopes import germline_peptide_source

        source = germline_peptide_source(9)
        cores = {e.core for e in lambda_chain.one_substitution("germline")}
        assert cores
        for core in cores:
            assert core not in source
            best = min(sum(a != b for a, b in zip(core, p, strict=True)) for p in source)
            assert best == 1

    def test_the_murine_chain_is_genuinely_distant_not_a_near_variant(self, murine):
        """Its framework sits well clear of anything human: every
        unmatched hit is two or more residues from any germline peptide."""
        counts = murine.counts_by_definition()
        assert counts["germline_one_substitution"] == 0
        assert counts["germline_unmatched"] == counts["germline_distant"]

    def test_unmatched_counts_separate_these_chains_by_origin(
        self, human, humanized, murine, lambda_chain
    ):
        """An observation about this panel, NOT a claim that origin order is the
        correct immunogenicity order. Species origin and MHC II presentation are
        different endpoints, and a human sequence can carry many predicted hits.
        This is a regression check on four fixed sequences.
        """
        assert murine.count("germline") > humanized.count("germline")
        assert humanized.count("germline") > lambda_chain.count("germline")
        assert lambda_chain.count("germline") > human.count("germline")

    def test_provenance_records_the_threshold_and_why(self, lambda_chain):
        record = lambda_chain.provenance()
        assert record["unmatched_threshold"] == 1
        assert "95.4%" in record["threshold_rationale"]
        assert "_distant" in record["threshold_rationale"]
        assert "matched + unmatched + unclassifiable" in record["count_reconciliation"]


class TestSyntheticDistanceCases:
    """Built by hand, so the distances are exactly what the test says.

    The real chains exercise the same logic, but their distances come out of
    the germline lookup; here they are set directly, which is what makes this a
    test of the counting rule rather than of the reference data.
    """

    @staticmethod
    def _epitope(core: str, distance: int) -> Epitope:
        binding = Binding(
            allele="HLA-DRB1*01:01",
            peptide=core + "AAAAAA",
            core=core,
            start=1,
            end=15,
            score=0.5,
            rank=1.0,
        )
        return Epitope(
            binding=binding,
            positions=tuple((i, " ") for i in range(1, 16)),
            self_distances=tuple((d, distance) for d in SELF_DEFINITIONS),
            chain_type="H",
        )

    @pytest.fixture
    def mixed(self):
        """One exact match, one at distance 1, one at 2+."""
        return EpitopeBurden(
            epitopes=(
                self._epitope("CCCCCCCCC", 0),
                self._epitope("AAAAAAAAA", 1),
                self._epitope("BBBBBBBBB", FAR),
            ),
            alleles=("HLA-DRB1*01:01",),
            method="netmhciipan_el",
            peptide_length=15,
            binder_rank=10.0,
            chain_type="H",
            n_predictions=3,
        )

    def test_both_the_one_mismatch_and_the_two_plus_core_are_unmatched(self, mixed):
        """The headline behaviour of this module, stated at its simplest."""
        cores = {e.core for e in mixed.unmatched("germline")}
        assert cores == {"AAAAAAAAA", "BBBBBBBBB"}

    def test_the_exact_match_is_not_unmatched(self, mixed):
        assert [e.core for e in mixed.matched("germline")] == ["CCCCCCCCC"]

    def test_the_exploratory_count_drops_the_one_mismatch_core(self, mixed):
        assert [e.core for e in mixed.distant("germline")] == ["BBBBBBBBB"]

    def test_the_categories_partition_and_the_subsets_sum_back(self, mixed):
        counts = mixed.counts_by_definition()
        assert counts["germline_matched"] == 1
        assert counts["germline_unmatched"] == 2
        assert counts["germline_unclassifiable"] == 0
        assert counts["germline_matched"] + counts["germline_unmatched"] == counts["raw"]
        assert (
            counts["germline_one_substitution"] + counts["germline_distant"]
            == counts["germline_unmatched"]
        )

    def test_a_distance_beyond_one_is_labelled_2plus_not_2(self, mixed):
        labels = [e.distance_label("germline") for e in mixed.epitopes]
        assert labels == ["0", "1", "2+"]
