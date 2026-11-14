"""Tests for the end-to-end report.

Offline throughout: epitope predictions come from the committed cache under
tests/data/iedb, and every other score is computed locally.
"""

import pytest

from abscore.number import NumberingError
from abscore.report import (
    analyse,
    escape_markdown,
    render_markdown,
    render_summary,
    render_text,
)
from tests.fixtures import (
    ABAGOVOMAB_VH,
    ADALIMUMAB_VH,
    ALBUMIN_FRAGMENT,
    RITUXIMAB_VH,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
)

CACHE = "tests/data/iedb"


def run(sequence, **kwargs):
    return analyse(sequence, cache_dir=CACHE, offline=True, **kwargs)


@pytest.fixture(scope="module")
def humanized():
    return run(TRASTUZUMAB_VH)


class TestEndToEnd:
    def test_one_call_produces_every_score(self, humanized):
        assert humanized.germline.identity is not None
        assert humanized.humanness.score is not None
        assert humanized.epitopes is not None

    def test_light_chains_work_too(self):
        report = run(TRASTUZUMAB_VL, with_epitopes=False)
        assert report.chain_type == "K"
        assert "kappa" in report.chain_name

    def test_a_non_antibody_is_refused(self):
        with pytest.raises(NumberingError):
            run(ALBUMIN_FRAGMENT)

    def test_an_scfv_is_refused_with_a_useful_message(self):
        scfv = TRASTUZUMAB_VH + "GGGGS" * 3 + TRASTUZUMAB_VL
        with pytest.raises(NumberingError, match="one chain at a time"):
            run(scfv)


class TestGracefulDegradation:
    """A free public service must not be a hard dependency."""

    def test_an_unavailable_predictor_does_not_lose_the_rest(self):
        report = run(RITUXIMAB_VH)  # no cached epitopes for this fixture
        assert report.epitopes is None
        assert report.epitope_error
        # Everything computed offline still present.
        assert report.germline.identity is not None
        assert report.humanness.score is not None

    def test_the_markdown_says_which_service_was_unreachable(self):
        text = render_markdown(run(RITUXIMAB_VH))
        assert "Not available" in text
        assert "IEDB" in text
        assert "computed offline and is unaffected" in text

    def test_epitopes_can_be_skipped_outright(self):
        report = run(TRASTUZUMAB_VH, with_epitopes=False)
        assert report.epitopes is None
        assert report.epitope_error is None

    def test_the_summary_marks_epitopes_unavailable(self):
        assert "binders=unavailable" in render_summary(run(RITUXIMAB_VH))


class TestNoBareNumbers:
    """Every score must carry its basis."""

    def test_germline_identity_states_its_denominator_and_gene(self, humanized):
        text = render_markdown(humanized)
        scope = humanized.germline.scopes["V"]
        assert f"{scope.matches}/{scope.compared} positions" in text
        assert humanized.germline.gene in text

    def test_oasis_states_its_threshold_and_peptide_count(self, humanized):
        text = render_markdown(humanized)
        assert f"{humanized.humanness.n_found}/{humanized.humanness.n_peptides}" in text
        assert "OAS donors" in text

    def test_the_epitope_section_states_method_alleles_and_cutoff(self, humanized):
        text = render_markdown(humanized)
        assert humanized.epitopes.method in text
        assert "HLA-DRB1" in text
        assert "percentile rank" in text

    def test_framework_scope_is_explained_rather_than_assumed(self, humanized):
        assert "FR4 is J-encoded" in render_markdown(humanized)

    def test_cdr3_exclusion_from_germline_is_explained(self, humanized):
        assert "junctional in origin" in render_markdown(humanized)


class TestHonestFraming:
    def test_the_numbering_hmm_species_is_kept_out_of_the_report(self, humanized):
        """It fed no score and misled on every reading.

        ANARCI picks whichever of its eight germline profiles fits best, purely
        to assign coordinates. The winner is a weak signal even on clean input
        -- pembrolizumab, a real humanized therapeutic, scores mouse 184.1
        against human 180.5 -- and on a corrupted sequence it is noise: an OAS
        record with a compensated indel came back alpaca. It needed a caption
        disclaiming it wherever it appeared, so it is reported nowhere.

        It stays in run.json, where traceability is the point.
        """
        for text in (render_markdown(humanized), render_text(humanized)):
            assert "numbering nearest neighbour" not in text
            assert "Numbering HMM" not in text

    def test_a_non_human_profile_call_never_reaches_the_report(self):
        """The sharp case. This OAS record carries a compensated indel, so 34
        residues translate out of frame and ANARCI's human profile drops from
        bitscore 141 to 80 -- below its own threshold, after which the alpaca
        profile wins. A reader must never see "alpaca" on a human B-cell read.
        """
        report = run(
            "GGSLRLSLCSLWISPSVAMLCTGSARLQARGWSGWQLYHMMGSNKYYADSREGRFTISRDNSKNTLY"
            "LQMNSLRAEDTAVYYCASLTDRHIVLIRGDWYFESLGAVAPWVTVSS",
            with_epitopes=False,
        )
        assert report.domain.species == "alpaca", "fixture no longer exercises the case"
        for text in (render_markdown(report), render_text(report)):
            assert "alpaca" not in text

    def test_all_three_self_definitions_are_shown(self, humanized):
        text = render_markdown(humanized)
        for definition in ("germline", "repertoire", "proteome"):
            assert f"*{definition}*" in text

    def test_the_junctional_abstention_is_explained(self, humanized):
        text = render_markdown(humanized)
        assert "unclassifiable, *germline*" in text
        assert "no single protein sequence" in text

    def test_flagged_runs_are_presented_without_a_cause(self, humanized):
        """A grafted CDR, heavy somatic mutation and a sequencing artefact all
        produce a run, so the report states the observation and stops."""
        text = render_markdown(humanized)
        assert "an observation" in text
        assert "carries no cause" in text

    def test_the_limits_section_names_what_sequence_cannot_show(self, humanized):
        text = render_markdown(humanized)
        for phrase in ("aggregation", "dose", "route of administration", "HLA type"):
            assert phrase in text
        assert "where to spend an ex vivo" in text

    def test_an_ambiguous_gene_call_is_surfaced(self):
        """Rituximab's heavy chain is an exact tie between IGHV1-46 and IGHV1-3."""
        report = run(RITUXIMAB_VH, with_epitopes=False)
        assert report.germline.is_ambiguous
        assert "cannot be separated from" in render_markdown(report)


class TestWarningsTravel:
    def test_warnings_are_collected_from_every_module(self):
        report = run(RITUXIMAB_VH, with_epitopes=False)
        assert report.warnings()
        assert all(isinstance(w, str) for w in report.warnings())

    def test_a_truncated_chain_reports_one_finding_not_four(self):
        """Partial domain, reduced coverage, fewer peptides and fewer windows
        are one fact stated four times. The consequences are folded into the
        cause so an independent caveat is not buried under repetitions."""
        report = run(TRASTUZUMAB_VH[17:], with_epitopes=False)
        partial = [w for w in report.warnings() if w.startswith("partial V-domain")]
        assert len(partial) == 1
        assert "Consequences:" in partial[0]
        assert "coverage" in partial[0]
        assert not [w for w in report.warnings() if w.startswith("germline identity computed at")]
        assert not [w for w in report.warnings() if w.startswith("humanness computed on a partial")]

    def test_a_complete_chain_keeps_its_caveats_unfolded(self):
        """Folding applies only when there is a partial-domain cause to fold
        into -- a coverage warning on a complete chain stands on its own."""
        report = run(RITUXIMAB_VH, with_epitopes=False)
        assert not any("Consequences:" in w for w in report.warnings())

    def test_warnings_are_deduplicated(self):
        report = run(ABAGOVOMAB_VH)
        assert len(report.warnings()) == len(set(report.warnings()))

    def test_every_warning_reaches_the_markdown(self, humanized):
        """Escaped, because allele names contain markdown-active asterisks --
        see TestAlleleNamesSurviveRendering."""
        text = render_markdown(humanized)
        for message in humanized.warnings():
            assert escape_markdown(message) in text

    def test_a_clean_chain_says_so_rather_than_showing_nothing(self):
        report = run(ADALIMUMAB_VH)
        assert report.warnings() == ()
        assert "No sequence-level warnings" in render_markdown(report)


class TestProvenance:
    def test_provenance_covers_every_section_that_produced_a_number(self, humanized):
        record = humanized.provenance()
        assert set(record) == {"germline_identity", "humanness", "epitopes"}

    def test_an_unavailable_service_is_recorded_as_such(self):
        record = run(RITUXIMAB_VH).provenance()
        assert record["epitopes"]["status"] == "unavailable"
        assert record["epitopes"]["reason"]


class TestSummaryLine:
    def test_summary_is_one_line_of_key_values(self, humanized):
        line = render_summary(humanized)
        assert "\n" not in line
        for key in ("chain=", "germline=", "V=", "FR=", "OASis="):
            assert key in line

    def test_summary_reports_the_epitope_units_separately(self, humanized):
        line = render_summary(humanized)
        assert "hits=" in line and "cores=" in line and "unmatched_germline=" in line


class TestAlleleNamesSurviveRendering:
    """Every IMGT allele name contains an asterisk, which markdown eats.

    `IGHV3-30-3*03 ... IGHV3-30*17` puts a matched pair of asterisks in one
    string, so a renderer italicises the text between them and deletes both
    allele numbers. A caveat that names the wrong gene is worse than no caveat,
    so the escaping is tested rather than assumed.
    """

    def test_asterisks_in_a_message_are_escaped(self):
        message = "ambiguous: IGHV3-30-3*03 cannot be separated from IGHV3-30*17"
        assert escape_markdown(message).count("\\*") == 2

    def test_escaping_is_idempotent_in_effect_not_silently_lossy(self):
        """A backslash must itself be escaped, or round-tripping mangles it."""
        assert escape_markdown("a\\b") == "a\\\\b"

    def test_the_rendered_caveat_keeps_every_allele_number(self):
        """End to end on a real tie: no allele number is lost on the way out.

        Rituximab's heavy chain ties IGHV1-46 against IGHV1-3, so its caveat
        genuinely names two alleles and genuinely contains a matched pair of
        asterisks -- the exact shape that used to be swallowed. Checked on the
        caveat line itself, because bold markers and code spans elsewhere in the
        report are legitimate markdown.
        """
        report = run(RITUXIMAB_VH, with_epitopes=False)
        assert report.germline.is_ambiguous
        alleles = (report.germline.gene, report.germline.ambiguous_genes[0].best_allele)

        text = render_markdown(report)

        # The caveat line protects allele names by escaping.
        caveat = next(line for line in text.splitlines() if "call is ambiguous" in line)
        for allele in alleles:
            assert allele.replace("*", "\\*") in caveat, f"{allele} corrupted in: {caveat}"
        assert "*" not in caveat.replace("\\*", ""), caveat

        # The chain section protects them with a code span instead, which is
        # the idiomatic markdown for an identifier.
        chain = next(line for line in text.splitlines() if "cannot be separated from:" in line)
        assert f"`{alleles[1]}`" in chain, chain


class TestTheRegionTableExplainsItself:
    """ "CDR2 | 56-65 | SAS" invites the question this class answers.

    The IMGT column is a fixed slot range; the sequence is what is present. A
    kappa CDR2 fills 3 of its 10 slots, which looks like missing data until the
    table says how many slots were filled.
    """

    def test_a_kappa_cdr2_reports_three_of_ten_slots(self):
        text = render_markdown(run(TRASTUZUMAB_VL, with_epitopes=False))
        row = next(line for line in text.splitlines() if line.startswith("| CDR2 "))
        assert "56-65" in row
        assert "3/10" in row

    def test_a_complete_framework_region_reports_full(self):
        text = render_markdown(run(TRASTUZUMAB_VL, with_epitopes=False))
        row = next(line for line in text.splitlines() if line.startswith("| FR2 "))
        assert "17/17" in row

    def test_structural_gaps_are_named_so_they_are_not_read_as_deletions(self):
        """FR3 shows 36/39 on a clean kappa chain, because IMGT 73, 81 and 82
        are empty in essentially every kappa. Without that said, the table
        implies three deletions on a chain the pipeline reports as clean."""
        report = run(TRASTUZUMAB_VL, with_epitopes=False)
        text = render_markdown(report)
        row = next(line for line in text.splitlines() if line.startswith("| FR3 "))
        assert "36/39" in row
        assert not [w for w in report.warnings() if "deletion" in w]
        for position in ("73", "81", "82"):
            assert position in text
        assert "not deletions" in text or "no deletion is reported" in text

    def test_the_kabat_boundary_difference_is_stated(self):
        """The commonest way to misread a short CDR-L2."""
        text = render_markdown(run(TRASTUZUMAB_VL, with_epitopes=False))
        assert "Kabat" in text


class TestTerminalRendering:
    """Markdown tables are unreadable in a terminal, so the report renders twice.

    Not a degraded copy: the text form carries the same numbers in the same
    order, with the region track as ASCII rather than a pipe table.
    """

    @pytest.mark.parametrize("with_epitopes", [False, True])
    def test_no_markdown_syntax_leaks_into_the_text_form(self, with_epitopes):
        """Literal asterisks and pipes are exactly what makes markdown
        unreadable on a terminal.

        Parametrised over the epitope section because the first version of this
        test only covered a report without one, and a bold marker duly shipped
        inside the section it skipped.
        """
        text = render_text(run(TRASTUZUMAB_VH, with_epitopes=with_epitopes))
        assert "**" not in text
        assert "|" not in text
        assert "---" not in text

    def test_the_sections_are_in_the_same_order_as_the_markdown(self):
        report = run(TRASTUZUMAB_VL, with_epitopes=False)
        text = render_text(report)
        order = [s for s in ("CHAIN", "REGIONS", "HUMANNESS", "MHC CLASS II", "CAVEATS")]
        positions = [text.index(s) for s in order]
        assert positions == sorted(positions)

    def test_the_headline_numbers_match_the_markdown(self):
        report = run(TRASTUZUMAB_VL, with_epitopes=False)
        text = render_text(report)
        assert report.germline.gene in text
        assert f"{report.germline.identity:.1%}" in text
        assert f"{report.humanness.score:.1%}" in text

    def test_the_region_track_marks_a_short_loop_as_normal(self):
        """A kappa CDR2 shows its three residues at the ends with dots between,
        and no '!' -- which is reserved for a real framework gap."""
        text = render_text(run(TRASTUZUMAB_VL, with_epitopes=False))
        row = next(line for line in text.splitlines() if line.strip().startswith("CDR2"))
        assert "SA" in row and row.strip().endswith(("3/10",))
        assert "!" not in row

    def test_a_real_framework_deletion_shows_the_finding_mark(self):
        from abscore.number import number_sequence

        domain = number_sequence(TRASTUZUMAB_VH)[0]
        text = render_text(run(domain.seq[:75] + domain.seq[76:], with_epitopes=False))
        fr3 = next(line for line in text.splitlines() if line.strip().startswith("FR3"))
        assert "!" in fr3

    def test_an_unavailable_epitope_section_says_the_rest_is_unaffected(self):
        text = render_text(run(TRASTUZUMAB_VL, with_epitopes=False))
        assert "Not available" in text
        assert "computed offline" in text

    def test_a_chain_with_no_caveats_says_so_rather_than_showing_nothing(self):
        text = render_text(run(ADALIMUMAB_VH, with_epitopes=False))
        assert "complete, unambiguous" in text

    def test_lines_stay_within_terminal_width(self):
        """A wrapped line is the thing this renderer exists to avoid."""
        for sequence in (TRASTUZUMAB_VL, TRASTUZUMAB_VH, ADALIMUMAB_VH):
            text = render_text(run(sequence, with_epitopes=False))
            too_long = [line for line in text.splitlines() if len(line) > 100]
            assert not too_long, too_long


class TestTheReferenceTableReadsInPlainWords:
    """ "Unmatched" is precise and opaque. The terminal table says "found" and
    "not found", names each reference set in plain words, and counts cores
    rather than hits.
    """

    def _block(self, report) -> str:
        text = render_text(report, self_definition="germline")
        return text[text.index("MHC CLASS II") : text.index("CAVEATS")]

    @pytest.fixture
    def block(self):
        return self._block(run(TRASTUZUMAB_VH, with_epitopes=True))

    def test_the_columns_are_plain_english(self, block):
        assert "found" in block and "not found" in block
        assert "unmatched" not in block.lower()
        assert "non-self" not in block.lower()

    def test_each_reference_set_is_named_rather_than_keyed(self, block):
        for label in ("germline Ig genes", "OAS repertoire", "human proteome"):
            assert label in block

    def test_the_rows_count_cores_and_say_so(self, block):
        """A hit is one peptide-allele pair, so hits overstate how many
        distinct peptides there are -- YLQMNSLRA alone is 27 of 86 hits."""
        assert "distinct 9-mer cores" in block
        assert "These rows count cores" in block

    def test_the_row_reconciles_to_the_core_total(self):
        report = run(TRASTUZUMAB_VH, with_epitopes=True)
        cores = report.epitopes.core_counts_by_definition()
        for definition in ("germline", "repertoire", "proteome"):
            total = (
                cores[f"{definition}_matched"]
                + cores[f"{definition}_unmatched"]
                + cores[f"{definition}_unclassifiable"]
            )
            assert total == cores["distinct_cores"], definition

    def test_the_subsets_sum_to_not_found(self):
        report = run(TRASTUZUMAB_VH, with_epitopes=True)
        cores = report.epitopes.core_counts_by_definition()
        for definition in ("germline", "repertoire", "proteome"):
            assert (
                cores[f"{definition}_one_substitution"] + cores[f"{definition}_distant"]
                == cores[f"{definition}_unmatched"]
            ), definition

    def test_cores_are_never_more_numerous_than_hits(self):
        report = run(TRASTUZUMAB_VH, with_epitopes=True)
        hits = report.epitopes.counts_by_definition()
        cores = report.epitopes.core_counts_by_definition()
        for key in ("germline_matched", "germline_unmatched", "germline_distant"):
            assert cores[key] <= hits[key], key

    def test_not_found_is_presented_as_a_lookup_result(self, block):
        """The reference tells you whether the 9-mer is in it. Everything past
        that -- tolerance, presentation, a response -- is unaddressed."""
        assert "lookup result" in block
        assert "stays on it" in block

    def test_the_selected_reference_is_marked(self, block):
        assert "germline Ig genes *" in block


class TestSeveralPartialSignalsBecomeOneFinding:
    """`completeness_warnings` runs four independent checks for a partial
    domain, so a badly truncated chain trips all four. They describe one fact.

    An earlier version attached the consequences paragraph to each of them,
    printing it four times on exactly the chains that needed careful reading.
    """

    #: Heavily truncated and misaligned: trips every partial-domain signal and
    #: carries a cysteine at the conserved IMGT 41 tryptophan.
    BADLY_TRUNCATED = (
        "GGSLRLSLCSLWISPSVAMLCTGSARLQARGWSGWQLYHMMGSNKYYADSREGRFTISRDNSKNTLYLQMNSLRAEDT"
        "AVYYCASLTDRHIVLIRGDWYFESLGAVAPWVTVSS"
    )

    @pytest.fixture
    def report(self):
        return run(self.BADLY_TRUNCATED, with_epitopes=False)

    def test_all_four_signals_fire_on_this_chain(self, report):
        """If they stop firing, this fixture no longer tests the merge."""
        partial = [w for w in report.domain.warnings if w.startswith("partial V-domain")]
        assert len(partial) >= 4, partial

    def test_they_arrive_as_one_caveat(self, report):
        partial = [w for w in report.warnings() if w.startswith("partial V-domain")]
        assert len(partial) == 1

    def test_every_signal_survives_the_merge(self, report):
        """Merging must not drop evidence: all four specifics stay readable."""
        merged = next(w for w in report.warnings() if w.startswith("partial V-domain"))
        for fragment in (
            "anchor position(s) 23, 118",
            "FR1, FR4",
            "starts at IMGT 27",
            "ends at IMGT 117",
        ):
            assert fragment in merged, fragment

    def test_the_consequences_appear_exactly_once(self, report):
        assert sum(w.count("Consequences:") for w in report.warnings()) == 1

    def test_the_independent_caveats_stay_separate(self, report):
        """A misaligned anchor and an ambiguous gene call are their own
        findings and must not be swallowed by the merge."""
        warnings = report.warnings()
        assert any("41=C" in w for w in warnings)
        assert any("ambiguous" in w for w in warnings)
        assert len(warnings) == 3

    def test_the_rendered_count_matches_the_findings(self, report):
        text = render_text(report)
        assert "CAVEATS (3)" in text


class TestAConsequenceNeedsItsCause:
    """A complete domain carrying one odd residue must not be called partial.

    The bug this pins: `humanness` and `epitopes` appended "computed on a
    partial domain" whenever numbering said *anything at all*. An unexpected
    residue at a conserved anchor says a residue is WRONG, not that sequence is
    MISSING, so a full-length chain with one substitution reported two caveats
    claiming its scores were computed over a reduced chain. Both were false:
    every window was scored.

    It mattered because the caveat block is what a reader uses to decide how
    far to trust the numbers, and the one caveat worth reading -- an unexpected
    residue at the V-J junction anchor -- was buried between two that were not.

    The sequence below is a real OAS lambda record. It numbers to a complete
    IMGT 1-127 span and carries C at position 118, where J-PHE/J-TRP belongs.
    """

    OAS_LAMBDA_ODD_ANCHOR = (
        "QFLLTQPPSTSGTPGQTVTISCSGRTSNIGKNTVNWYQQLPGTAPKLLIYNDNQRPSGVPARFSGSESGTSASLAIS"
        "GLQSDDEAEYYCLAWDDSLNGWACGGGTKLTVL"
    )

    @pytest.fixture(scope="class")
    def report(self):
        return run(self.OAS_LAMBDA_ODD_ANCHOR, with_epitopes=False)

    def test_the_chain_really_is_complete(self, report):
        """If this stops holding, the fixture no longer tests anything."""
        assert report.domain.imgt_span == (1, 127)
        assert not [w for w in report.domain.warnings if w.startswith("partial V-domain")], (
            report.domain.warnings
        )

    def test_the_odd_anchor_is_still_reported(self, report):
        """The finding that deserved the attention the noise was taking."""
        assert any("118=C" in w for w in report.warnings())

    def test_humanness_does_not_claim_a_partial_domain(self, report):
        assert not [w for w in report.humanness.warnings() if "partial domain" in w]

    def test_the_epitope_count_does_not_claim_a_partial_domain(self):
        report = run(self.OAS_LAMBDA_ODD_ANCHOR)
        if report.epitopes is None:
            pytest.skip("no cached IEDB response for this chain")
        assert not [w for w in report.epitopes.warnings() if "partial domain" in w]

    def test_the_wrong_anchor_is_the_only_caveat(self, report):
        assert len(report.warnings()) == 1, report.warnings()

    def test_every_window_was_scored(self, report):
        """The claim the caveats made was that fewer windows were used. It was
        arithmetically false, and this says so in numbers."""
        residues = len(report.domain.numbering)
        assert report.humanness.n_peptides == residues - report.humanness.peptide_length + 1


class TestAnOrphanedConsequenceIsDropped:
    """The backstop, tested directly on the merge.

    Fixing the two modules clears today's case. This pins the layer underneath:
    if a future module appends a bare "computed on a partial domain" with no
    partial-domain warning behind it, the report drops it rather than printing
    a claim about the input that nothing supports.
    """

    class _Appends:
        """Stands in for a module that appends one extra warning.

        The score objects are frozen, so the swap happens on the report, which
        is the seam the merge actually reads through."""

        def __init__(self, wrapped, extra):
            self._wrapped = wrapped
            self._extra = extra

        def warnings(self):
            return (*self._wrapped.warnings(), self._extra)

        def __getattr__(self, name):
            return getattr(self._wrapped, name)

    def _warnings(self, report, extra):
        original = report.humanness
        report.humanness = self._Appends(original, extra)
        try:
            return report.warnings()
        finally:
            report.humanness = original

    @pytest.fixture(scope="class")
    def complete(self):
        return run(ADALIMUMAB_VH, with_epitopes=False)

    def test_a_restatement_with_no_cause_is_dropped(self, complete):
        orphan = "humanness computed on a partial domain: the score is over fewer peptides."
        assert orphan not in self._warnings(complete, orphan)

    def test_a_measured_consequence_is_kept(self, complete):
        """`germline identity computed at X% coverage` carries a figure taken
        from the alignment, so it stands on its own and is not a restatement."""
        measured = "germline identity computed at 62.0% coverage of the V gene."
        assert measured in self._warnings(complete, measured)

    def test_an_unrelated_caveat_is_untouched(self, complete):
        other = "something else entirely worth saying."
        assert other in self._warnings(complete, other)

    def test_classification_does_not_depend_on_message_order(self, complete):
        """An earlier version only recognised a consequence once a partial
        warning had already been seen, so the result depended on which module
        spoke first."""
        orphan = "epitope burden computed on a partial domain: fewer windows."
        assert orphan not in self._warnings(complete, orphan)
