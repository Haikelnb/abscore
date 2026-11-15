"""Tests for IMGT numbering and region assignment.

No network: ANARCI runs locally against bundled HMMs.
"""

import pytest

from abscore.number import (
    ANCHOR_POSITIONS,
    CDR_REGIONS,
    FRAMEWORK_REGIONS,
    IMGT_REGIONS,
    NUCLEOTIDE_MIN_LENGTH,
    STRUCTURAL_FRAMEWORK_GAPS,
    HmmerNotFoundError,
    NumberingError,
    imgt_sort_key,
    number_sequence,
    number_single_chain,
    region_of_position,
)
from tests.fixtures import (
    ABELACIMAB_VL_LAMBDA,
    ELIPOVIMAB_VH,
    ELIPOVIMAB_VH_CDR3,
    RITUXIMAB_VH,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
    VH_WITH_TAIL,
)


class TestRegionBoundaries:
    def test_regions_tile_1_to_128_without_overlap(self):
        covered = []
        for lo, hi in IMGT_REGIONS.values():
            covered.extend(range(lo, hi + 1))
        assert sorted(covered) == list(range(1, 129))

    def test_region_names_partition_into_fr_and_cdr(self):
        assert set(FRAMEWORK_REGIONS) | set(CDR_REGIONS) == set(IMGT_REGIONS)

    @pytest.mark.parametrize(
        ("position", "expected"),
        [
            (1, "FR1"),
            (26, "FR1"),
            (27, "CDR1"),
            (38, "CDR1"),
            (39, "FR2"),
            (55, "FR2"),
            (56, "CDR2"),
            (65, "CDR2"),
            (66, "FR3"),
            (104, "FR3"),
            (105, "CDR3"),
            (117, "CDR3"),
            (118, "FR4"),
            (128, "FR4"),
        ],
    )
    def test_boundary_positions(self, position, expected):
        assert region_of_position(position) == expected

    def test_positions_outside_v_domain_have_no_region(self):
        assert region_of_position(0) is None
        assert region_of_position(129) is None

    def test_accepts_position_with_insertion_code(self):
        assert region_of_position((111, "A")) == "CDR3"


class TestHeavyChain:
    def test_calls_heavy_chain(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert domain.chain_type == "H"
        assert domain.is_heavy and not domain.is_light

    def test_every_input_residue_is_numbered(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert domain.seq == TRASTUZUMAB_VH
        assert len(domain.numbering) == len(TRASTUZUMAB_VH)

    def test_regions_reconstruct_the_domain_in_order(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        rebuilt = "".join(domain.regions[r] for r in IMGT_REGIONS)
        assert rebuilt == domain.seq

    def test_framework_and_cdrs_are_disjoint_and_complete(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert len(domain.framework_seq()) + len(domain.cdr_seq()) == len(domain.seq)

    def test_assigns_a_human_germline_v_gene(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert domain.v_gene is not None
        assert domain.v_gene.gene.startswith("IGHV")
        assert 0.0 <= domain.v_gene.identity <= 1.0

    def test_cdr3_starts_after_the_conserved_cysteine(self):
        # IMGT position 104 is the conserved second cysteine of the V-domain;
        # CDR3 begins at 105. This is a load-bearing invariant for CDR masking.
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert domain.numbering[(104, " ")] == "C"


class TestLightChain:
    def test_calls_light_chain(self):
        domain = number_single_chain(TRASTUZUMAB_VL)
        assert domain.chain_type in ("K", "L")
        assert domain.is_light and not domain.is_heavy

    def test_assigns_a_light_germline_v_gene(self):
        domain = number_single_chain(TRASTUZUMAB_VL)
        assert domain.v_gene is not None
        assert domain.v_gene.gene.startswith(("IGKV", "IGLV"))


class TestVDomainTrimming:
    def test_constant_region_is_excluded_from_numbering(self):
        """A full-length chain must be trimmed to its V-domain before scoring."""
        domain = number_single_chain(VH_WITH_TAIL)
        assert len(domain.seq) < len(VH_WITH_TAIL)
        assert domain.seq == VH_WITH_TAIL[domain.query_start : domain.query_end]

    def test_query_slice_maps_back_onto_the_input(self):
        domain = number_single_chain(VH_WITH_TAIL)
        assert VH_WITH_TAIL[domain.query_start : domain.query_end] == domain.seq


class TestInputValidation:
    def test_empty_sequence_is_rejected(self):
        with pytest.raises(NumberingError, match="empty"):
            number_single_chain("")

    def test_whitespace_and_case_are_normalised(self):
        spaced = f"  {TRASTUZUMAB_VH[:60].lower()}\n{TRASTUZUMAB_VH[60:].lower()}  "
        assert number_single_chain(spaced).seq == TRASTUZUMAB_VH

    @pytest.mark.parametrize("bad", ["-", "*", "@", "1"])
    def test_non_amino_acid_characters_are_rejected(self, bad):
        with pytest.raises(NumberingError, match="not amino acids"):
            number_single_chain(TRASTUZUMAB_VH[:50] + bad + TRASTUZUMAB_VH[50:])

    @pytest.mark.parametrize("code", ["X", "B", "Z"])
    def test_a_single_ambiguity_code_is_tolerated_and_warned(self, code):
        """Refusing a whole chain over one uncertain residue is disproportionate.

        Rejecting any sequence containing an X made 3,469 of 17,786 real OAS
        repertoire records (19.5%) unscorable.
        """
        sequence = TRASTUZUMAB_VH[:50] + code + TRASTUZUMAB_VH[51:]
        domain = number_single_chain(sequence)
        assert domain.chain_type == "H"
        assert any("ambiguity code" in w for w in domain.warnings)

    def test_too_many_ambiguity_codes_inside_the_domain_is_refused(self):
        mangled = TRASTUZUMAB_VH[:40] + "X" * 12 + TRASTUZUMAB_VH[52:]
        with pytest.raises(NumberingError, match="above the"):
            number_single_chain(mangled)

    def test_ambiguity_outside_the_extracted_domain_does_not_refuse_it(self):
        """The limit measures the sequence that gets scored.

        Ten X residues at the N-terminus push the domain boundary past them,
        so nine fall outside it and the extracted domain carries one. Refusing
        that would be refusing a clean domain over residues no score touches.
        The truncation it does cause is reported on its own.
        """
        mangled = "X" * 10 + TRASTUZUMAB_VH[10:]
        domain = number_single_chain(mangled)
        assert domain.seq.count("X") == 1
        assert domain.residues_outside_domain == 9

    def test_the_ambiguity_limit_is_a_fraction_of_the_sequence(self):
        from abscore.number import MAX_AMBIGUOUS_FRACTION

        assert 0.0 < MAX_AMBIGUOUS_FRACTION < 0.1

    def test_non_antibody_sequence_is_rejected(self):
        # Human serum albumin N-terminus: a real protein from elsewhere.
        albumin = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS"
        with pytest.raises(NumberingError, match="no antibody V-domain"):
            number_single_chain(albumin)


class TestMultipleDomains:
    def test_scfv_yields_two_domains(self):
        scfv = TRASTUZUMAB_VH + "GGGGSGGGGSGGGGS" + TRASTUZUMAB_VL
        domains = number_sequence(scfv)
        assert len(domains) == 2
        assert {d.chain_type for d in domains} == {"H", "K"}

    def test_single_chain_helper_refuses_an_scfv(self):
        """Guards against silently scoring only the first domain of a construct."""
        scfv = TRASTUZUMAB_VH + "GGGGSGGGGSGGGGS" + TRASTUZUMAB_VL
        with pytest.raises(NumberingError, match="2 V-domains"):
            number_single_chain(scfv)

    def test_domains_do_not_overlap_in_the_query(self):
        scfv = TRASTUZUMAB_VH + "GGGGSGGGGSGGGGS" + TRASTUZUMAB_VL
        first, second = sorted(number_sequence(scfv), key=lambda d: d.query_start)
        assert first.query_end <= second.query_start


class TestCompleteness:
    """A partial V-domain warns; it does not refuse."""

    def test_complete_chains_have_no_warnings(self):
        for seq in (TRASTUZUMAB_VH, TRASTUZUMAB_VL):
            domain = number_single_chain(seq)
            assert domain.warnings == ()
            assert domain.is_complete

    def test_all_conserved_anchors_present_in_complete_chains(self):
        heavy = number_single_chain(TRASTUZUMAB_VH)
        light = number_single_chain(TRASTUZUMAB_VL)
        assert heavy.anchors[23] == "C"
        assert heavy.anchors[41] == "W"
        assert heavy.anchors[104] == "C"
        assert heavy.anchors[118] == "W"  # J-TRP in heavy chains
        assert light.anchors[118] == "F"  # J-PHE in light chains
        assert set(heavy.anchors) == set(ANCHOR_POSITIONS)

    def test_c_terminal_truncation_warns_and_still_scores(self):
        domain = number_single_chain(TRASTUZUMAB_VH[:100])
        assert not domain.is_complete
        assert domain.numbering  # still numbered, with a warning
        assert any("FR4" in w for w in domain.warnings)
        assert any("anchor" in w and "118" in w for w in domain.warnings)

    def test_n_terminal_truncation_warns_even_though_all_regions_populated(self):
        """The case that region-emptiness alone cannot catch."""
        domain = number_single_chain(TRASTUZUMAB_VH[20:])
        assert all(domain.regions[r] for r in IMGT_REGIONS)  # nothing empty
        assert not domain.is_complete
        assert any("starts at IMGT" in w for w in domain.warnings)

    def test_heavy_chains_end_at_128_and_light_at_127(self):
        assert number_single_chain(TRASTUZUMAB_VH).imgt_span == (1, 128)
        assert number_single_chain(TRASTUZUMAB_VL).imgt_span == (1, 127)

    def test_gapped_positions_are_not_treated_as_missing(self):
        """IMGT leaves slots empty wherever a loop is short; that is not partial."""
        domain = number_single_chain(TRASTUZUMAB_VL)
        present = {p[0] for p in domain.numbering}
        assert len(set(range(1, 128)) - present) > 10  # many gaps, by design
        assert domain.is_complete


class TestInsertionOrdering:
    """IMGT numbers CDR3 outwards from the top of the loop.

    Insertions at 111 and below run ascending AFTER the base position; insertions
    at 112 and above run descending BEFORE it. A naive (number, letter) sort gets
    this backwards and transposes residues in long CDR3 loops without erroring.
    """

    def test_sort_key_matches_imgt_order_at_the_111_112_boundary(self):
        expected = [
            (111, " "),
            (111, "A"),
            (112, "A"),
            (112, " "),
        ]
        assert sorted(expected, key=imgt_sort_key) == expected

    def test_sort_key_handles_multi_letter_insertions_both_sides(self):
        expected = [
            (111, " "),
            (111, "A"),
            (111, "B"),
            (111, "C"),
            (112, "C"),
            (112, "B"),
            (112, "A"),
            (112, " "),
        ]
        assert sorted(expected, key=imgt_sort_key) == expected

    def test_naive_sort_would_be_wrong(self):
        """Pins the bug this function exists to prevent."""
        positions = [(112, " "), (112, "A")]
        assert sorted(positions) == [(112, " "), (112, "A")]  # naive
        assert sorted(positions, key=imgt_sort_key) == [(112, "A"), (112, " ")]

    def test_long_cdr3_is_reconstructed_in_the_right_order(self):
        domain = number_single_chain(ELIPOVIMAB_VH)
        assert domain.regions["CDR3"] == ELIPOVIMAB_VH_CDR3
        assert len(domain.regions["CDR3"]) == 26

    def test_long_cdr3_chain_still_reconstructs_the_input(self):
        domain = number_single_chain(ELIPOVIMAB_VH)
        assert domain.seq == ELIPOVIMAB_VH
        assert domain.is_complete

    def test_insertion_codes_are_present_and_span_both_sides(self):
        domain = number_single_chain(ELIPOVIMAB_VH)
        inserted = sorted({p[0] for p in domain.numbering if p[1] != " "})
        assert inserted == [111, 112]


class TestConservedAnchors:
    def test_all_five_imgt_anchor_positions_are_checked(self):
        # IMGT lists five conserved residues: 1st-CYS 23, CONSERVED-TRP 41,
        # hydrophobic 89, 2nd-CYS 104, J-PHE/J-TRP 118.
        assert set(ANCHOR_POSITIONS) == {23, 41, 89, 104, 118}

    def test_position_89_is_hydrophobic_in_every_fixture(self):
        for seq in (TRASTUZUMAB_VH, TRASTUZUMAB_VL, RITUXIMAB_VH, ELIPOVIMAB_VH):
            domain = number_single_chain(seq)
            assert domain.numbering[(89, " ")] in ANCHOR_POSITIONS[89]

    def test_mutated_anchor_warns_about_misalignment_not_truncation(self):
        """An anchor present but wrong is a different failure from one absent."""
        # Find the RAW index of IMGT position 23 -- they are not the same number,
        # which is the entire reason this module exists.
        reference = number_single_chain(TRASTUZUMAB_VH)
        ordered = sorted(reference.numbering, key=imgt_sort_key)
        raw_index = ordered.index((23, " "))
        assert TRASTUZUMAB_VH[raw_index] == "C"

        mutated = TRASTUZUMAB_VH[:raw_index] + "S" + TRASTUZUMAB_VH[raw_index + 1 :]
        domain = number_single_chain(mutated)
        assert any("misaligned" in w for w in domain.warnings)
        assert not any("missing conserved" in w for w in domain.warnings)


class TestNucleotideInput:
    """DNA passes the residue check (A/C/G/T are amino acids) so it needs its own.

    The sequence here is a real nucleotide VH fragment supplied during testing. It
    previously failed with "no antibody V-domain found", which pointed at the
    wrong problem.
    """

    DNA_VH = (
        "CTTCGGAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATGAGTAGTTATTACTGGGCCTGGATCCGG"
        "CAGCCCGCCGGGAAGCGACTGGAGTGGATTGGTCGTATCTATTCCAGTGGGGCACCAACTACAATCCCTCCCTC"
        "AAGAGTCGAGTCACCCATCTCAATAGACATGTCCAAGAGCCAATTCTCCCTGAACCTGAACTCTGTGACCGCCG"
        "CGGACACGGCCGTTTATTACTGTGCGAGAGATCCCGCCTGGTTCGACCCCTGGGGCCGAGGAATCCTGGTCACC"
        "GTCTCCTCAGGTAAATCGTGGCCGTGGGTCGCGTCTCTCAAGGGCACACAGGGGGGATAGGG"
    )

    def test_nucleotide_sequence_is_diagnosed_as_such(self):
        with pytest.raises(NumberingError, match="nucleotide sequence"):
            number_single_chain(self.DNA_VH)

    def test_error_does_not_claim_the_sequence_is_not_an_antibody(self):
        """It is an antibody -- in frame 2. The message must not say otherwise."""
        with pytest.raises(NumberingError) as excinfo:
            number_single_chain(self.DNA_VH)
        assert "no antibody V-domain" not in str(excinfo.value)

    def test_the_translated_protein_does_number(self):
        """Frame 2 of DNA_VH, confirming the input was a real antibody chain."""
        protein = (
            "SETLSLTCTVSGGSMSSYYWAWIRQPAGKRLEWIGRIYSSGAPTTIPPSRVESPISIDMSKSQFSLNL"
            "NSVTAADTAVYYCARDPAWFDPWGRGILVTVSS"
        )
        domain = number_single_chain(protein)
        assert domain.chain_type == "H"
        assert domain.v_gene is not None
        assert domain.v_gene.gene.startswith("IGHV4")
        # N-terminally truncated: starts mid-FR1, so it must warn.
        assert domain.imgt_span[0] > 1
        assert any("starts at IMGT" in w for w in domain.warnings)

    def test_a_short_acgt_run_is_not_mistaken_for_dna(self):
        """Guards the length floor: a short peptide must not trip the check."""
        short = "ACGT" * 5  # 20 chars, below the floor
        assert len(short) < NUCLEOTIDE_MIN_LENGTH
        with pytest.raises(NumberingError, match="no antibody V-domain"):
            number_single_chain(short)

    def test_a_real_protein_chain_is_not_mistaken_for_dna(self):
        assert number_single_chain(TRASTUZUMAB_VH).chain_type == "H"


class TestFrameworkDeletions:
    """A gap inside framework is a structural deletion; a gap in a CDR is not.

    CDR1-IMGT is 5-12 residues, so a short loop legitimately leaves 31-34 empty.
    Framework is fixed-length, so a missing framework position means something
    real -- except at the slots IMGT leaves empty by design per chain type.
    """

    OAS_FRAMESHIFTED_VH = (
        "SETLSLTCTVSGGSMSSYYWAWIRQPAGKRLEWIGRIYSSGAPTTIPPSRVESPISIDMSKSQFSLNL"
        "NSVTAADTAVYYCARDPAWFDPWGRGILVTVSS"
    )

    def test_framework_deletion_is_flagged(self):
        domain = number_single_chain(self.OAS_FRAMESHIFTED_VH)
        assert any("deletion inside framework" in w and "74" in w for w in domain.warnings)

    def test_position_73_is_a_structural_gap_not_a_deletion(self):
        """73 is absent in 100% of the 1142-domain panel, so it must not fire.

        OAS's own ANARCI_status lists 73 under "Deletions", which overstates it.
        """
        domain = number_single_chain(self.OAS_FRAMESHIFTED_VH)
        assert 73 not in {p[0] for p in domain.numbering}
        deletion_warnings = [w for w in domain.warnings if "deletion inside framework" in w]
        assert len(deletion_warnings) == 1
        assert "73" not in deletion_warnings[0]

    @pytest.mark.parametrize(
        "sequence",
        [
            TRASTUZUMAB_VH,
            TRASTUZUMAB_VL,
            RITUXIMAB_VH,
            ELIPOVIMAB_VH,
            ABELACIMAB_VL_LAMBDA,
        ],
    )
    def test_complete_chains_raise_no_deletion_warning(self, sequence):
        domain = number_single_chain(sequence)
        assert not any("deletion inside framework" in w for w in domain.warnings)
        assert domain.is_complete

    def test_cdr_gaps_are_never_reported_as_deletions(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        absent = set(range(1, 129)) - {p[0] for p in domain.numbering}
        in_cdrs = {p for p in absent if region_of_position(p) in CDR_REGIONS}
        assert in_cdrs  # there are some, by design
        assert domain.warnings == ()

    def test_structural_gaps_differ_by_chain_type(self):
        """81/82 are empty in every kappa and lambda chain, never in heavy."""
        assert 81 in STRUCTURAL_FRAMEWORK_GAPS["K"]
        assert 81 not in STRUCTURAL_FRAMEWORK_GAPS["H"]
        heavy = number_single_chain(TRASTUZUMAB_VH)
        kappa = number_single_chain(TRASTUZUMAB_VL)
        assert 81 in {p[0] for p in heavy.numbering}
        assert 81 not in {p[0] for p in kappa.numbering}


class TestLambdaChain:
    """Closes the gap left by the heavy/kappa-only fixtures."""

    def test_lambda_is_called_as_l(self):
        assert number_single_chain(ABELACIMAB_VL_LAMBDA).chain_type == "L"

    def test_lambda_ends_at_127_like_kappa(self):
        assert number_single_chain(ABELACIMAB_VL_LAMBDA).imgt_span == (1, 127)

    def test_lambda_carries_j_phe_not_j_trp(self):
        assert number_single_chain(ABELACIMAB_VL_LAMBDA).anchors[118] == "F"

    def test_lambda_v_gene_is_iglv(self):
        domain = number_single_chain(ABELACIMAB_VL_LAMBDA)
        assert domain.v_gene is not None
        assert domain.v_gene.gene.startswith("IGLV")


class TestTruncationIsNotADeletion:
    """Regression for a bug found on OAS repertoire data, not on therapeutics.

    The framework-deletion scan was bounded at the first observed position but
    not the last, so every missing FR4 position in a read that stopped early was
    reported as a "deletion inside framework" -- on top of the truncation warning
    that already described it. 19% of 14,317 OAS sequences were affected.
    Marketed therapeutics are complete, so the therapeutic panel could not
    surface this.
    """

    def test_c_terminal_truncation_reports_truncation_not_deletions(self):
        domain = number_single_chain(TRASTUZUMAB_VH[:100])
        assert not domain.is_complete
        assert any("ends at IMGT" in w for w in domain.warnings)
        assert not any("deletion inside framework" in w for w in domain.warnings)

    def test_n_terminal_truncation_reports_truncation_not_deletions(self):
        domain = number_single_chain(TRASTUZUMAB_VH[20:])
        assert any("starts at IMGT" in w for w in domain.warnings)
        assert not any("deletion inside framework" in w for w in domain.warnings)

    def test_truncation_at_both_ends_reports_neither_as_a_deletion(self):
        domain = number_single_chain(TRASTUZUMAB_VH[20:100])
        assert not any("deletion inside framework" in w for w in domain.warnings)

    def test_an_interior_gap_is_still_reported(self):
        """The signal must survive the fix: position 74 is interior to the span."""
        partial = (
            "SETLSLTCTVSGGSMSSYYWAWIRQPAGKRLEWIGRIYSSGAPTTIPPSRVESPISIDMSKSQFSLNL"
            "NSVTAADTAVYYCARDPAWFDPWGRGILVTVSS"
        )
        domain = number_single_chain(partial)
        assert any("deletion inside framework" in w and "74" in w for w in domain.warnings)


class TestMissingSystemDependency:
    """A missing `hmmscan` must say what to install, not leak errno 2.

    ANARCI shells out to HMMER, a C program installed outside Python, so a
    correct `uv sync` still leaves it missing. The bare
    `FileNotFoundError: 'hmmscan'` that surfaced before named no remedy, and it
    is the first thing anyone running this on a clean machine will hit.
    """

    def test_a_missing_hmmscan_names_the_dependency_and_how_to_get_it(self, monkeypatch):
        def no_hmmscan(*args, **kwargs):
            raise FileNotFoundError(2, "No such file or directory", "hmmscan")

        monkeypatch.setattr("abscore.number.run_anarci", no_hmmscan)
        with pytest.raises(HmmerNotFoundError) as caught:
            number_sequence(TRASTUZUMAB_VH)
        message = str(caught.value)
        assert "hmmscan" in message
        assert "brew install hmmer" in message
        assert "system dependency" in message

    def test_it_is_not_a_numbering_error(self):
        """The sequence is not what is wrong, so it must not be reported as a
        bad-input error -- the CLI exits 2 for those."""
        assert not issubclass(HmmerNotFoundError, NumberingError)

    def test_an_unrelated_missing_file_is_not_swallowed(self, monkeypatch):
        def other_failure(*args, **kwargs):
            raise FileNotFoundError(2, "No such file or directory", "some_other_thing")

        monkeypatch.setattr("abscore.number.run_anarci", other_failure)
        with pytest.raises(FileNotFoundError):
            number_sequence(TRASTUZUMAB_VH)


class TestDomainOrderDoesNotChangeTheResult:
    """A TCR beside an antibody must not depend on which came first.

    ANARCI numbers T-cell receptor variable domains too, and those are out of
    scope here, so they are skipped. The metadata list it returns is parallel
    to the RAW domain list; indexing it by the count of ACCEPTED domains means
    that after one skip, every later domain is read with its predecessor's
    metadata. The antibody was then seen as a TCR and skipped as well, and the
    whole input was rejected as "not an antibody" -- but only when the TCR came
    first. The same two domains in the other order scored normally.

    A result that depends on the order of the input is the shape of this bug,
    so that is what these assert.
    """

    # Synthetic TCR alpha variable domain, built from ANARCI's own bundled
    # human TRAV/TRAJ reference. Not a real receptor and not claimed to fold.
    TCR_ALPHA = (
        "AQSVTQLDSHVSVSEGTPVLLRCNYSSSVPPYLFWYVQYPNQGLQLLLKYTSAATLVKGINGFEAEFKKSE"
        "TSFHLTKPSAHMSDAAEYFCAVSEDSSYKLIFGSGTRLLVRP"
    )

    @pytest.fixture(scope="class")
    @classmethod
    def both_orders(cls):
        from abscore.number import number_sequence

        return (
            number_sequence(cls.TCR_ALPHA + TRASTUZUMAB_VH),
            number_sequence(TRASTUZUMAB_VH + cls.TCR_ALPHA),
        )

    def test_the_antibody_is_found_in_both_orders(self, both_orders):
        for domains in both_orders:
            assert [d.chain_type for d in domains] == ["H"]

    def test_the_extracted_domain_is_identical_in_both_orders(self, both_orders):
        first, second = both_orders
        assert first[0].seq == second[0].seq == TRASTUZUMAB_VH

    def test_the_metadata_follows_the_antibody_not_its_neighbour(self, both_orders):
        """The real failure: the H domain was read with the TCR's metadata.

        Asserting that the construct merely succeeds would pass even if the
        germline call came from the wrong domain, so this checks the call.
        """
        bare = number_single_chain(TRASTUZUMAB_VH)
        for domains in both_orders:
            assert domains[0].v_gene is not None
            assert domains[0].v_gene.gene == bare.v_gene.gene
            assert domains[0].species == bare.species

    def test_the_excluded_domain_is_disclosed(self, both_orders):
        """Skipping it silently would leave the caller unaware it was there."""
        for domains in both_orders:
            assert any("non-antibody variable domain" in w for w in domains[0].warnings)

    def test_a_lone_tcr_is_still_refused(self):
        """Disclosure is not acceptance: with no antibody there is nothing to
        score, and that must stay an error rather than an empty report."""
        with pytest.raises(NumberingError, match="not an antibody"):
            number_single_chain(self.TCR_ALPHA)


class TestWhatWasScoredIsStated:
    def test_a_bare_domain_says_nothing_extra(self):
        domain = number_single_chain(TRASTUZUMAB_VH)
        assert domain.residues_outside_domain == 0
        assert domain.input_length == len(TRASTUZUMAB_VH)

    def test_a_flanked_domain_reports_what_fell_outside(self):
        domain = number_single_chain(TRASTUZUMAB_VH + "GGGGSHHHHHH")
        assert domain.residues_outside_domain == 11

    def test_the_lambda_boundary_shift_is_pinned(self):
        """A tag pulls one residue into the lambda domain and moves the score.

        Both numbers are correct over the sequence each was computed on. What
        would be wrong is the report giving no way to tell the two runs apart,
        so this pins the shift and the disclosure together. If the HMM or its
        version changes this, the figures in `provenance.KNOWN_SENSITIVITIES`
        need rewriting -- which is the point of pinning it.
        """
        from abscore.humanness import score_domain

        bare = number_single_chain(ABELACIMAB_VL_LAMBDA)
        tagged = number_single_chain(ABELACIMAB_VL_LAMBDA + "GGGGSHHHHHH")
        assert len(tagged.seq) == len(bare.seq) + 1
        assert score_domain(tagged).score != score_domain(bare).score
        assert tagged.residues_outside_domain == 10
