"""Tests for run provenance.

The point of these is not that the fields exist, but that the record stays
honest: that it names what was used, admits what is unknown, and quantifies the
limitations rather than reassuring about them. A provenance file that quietly
drops an inconvenient field is worse than none.
"""

import json

import pytest

from abscore.provenance import (
    KNOWN_SENSITIVITIES,
    VERIFIED_AGAINST,
    file_digest,
    hmmer_version,
    panel_digests,
    run_metadata,
    settings,
    tool_versions,
    write_run_json,
)


class TestToolVersions:
    def test_every_tool_that_produces_a_number_is_recorded(self):
        tools = tool_versions()
        for name in ("abscore", "anarci", "promb", "biopython", "python"):
            assert tools[name], name

    def test_hmmer_is_recorded_even_though_it_is_not_a_python_package(self):
        """ANARCI shells out to it, so the numbering depends on its version."""
        assert hmmer_version() is not None
        assert tool_versions()["hmmer"] == hmmer_version()

    def test_hmmer_version_looks_like_a_version(self):
        assert hmmer_version()[0].isdigit()


class TestGermlineProvenance:
    def test_germline_data_is_pinned_by_hash(self):
        germline = run_metadata()["germline"]
        assert len(germline["germline_data_sha256"]) == 64
        assert germline["germline_data_matches_pin"] is True

    def test_the_missing_imgt_release_is_stated_not_omitted(self):
        """ANARCI ships no IMGT release number. The record must say so.

        Leaving the field out would let a reader assume it was checked.
        """
        assert run_metadata()["germline"]["imgt_release"] == "not recorded by ANARCI"


class TestSettings:
    def test_thresholds_are_written_out_rather_than_implied(self):
        recorded = settings()
        assert recorded["humanness_threshold_min_donor_fraction"] == 0.10
        assert recorded["humanness_threshold_name"] == "relaxed"
        assert 0 < recorded["max_ambiguous_fraction"] < 0.1

    def test_region_boundaries_are_recorded_so_a_cdr_claim_is_interpretable(self):
        """ "4 epitopes in CDRs" is meaningless without the scheme and bounds."""
        recorded = settings()
        assert recorded["numbering_scheme"] == "imgt"
        assert recorded["imgt_regions"]["CDR3"] == [105, 117]
        assert recorded["imgt_regions"]["FR1"] == [1, 26]

    def test_the_identity_denominator_is_stated(self):
        text = settings()["germline_identity_denominator"]
        assert "coverage" in text and "both query and germline" in text

    def test_germline_scopes_record_that_fr4_and_cdr3_are_excluded(self):
        scopes = settings()["germline_scopes"]
        assert "FR4" not in scopes["V"]
        assert "CDR3" not in scopes["V"]

    def test_structural_gap_rule_is_recorded_per_chain_type(self):
        gaps = settings()["structural_framework_gaps"]
        assert gaps["H"] == [10, 73]
        assert gaps["K"] == [73, 81, 82]
        assert gaps["L"] == [10, 73, 81, 82]


class TestVerifiedAgainst:
    def test_each_load_bearing_convention_names_a_source_and_a_date(self):
        for key, entry in VERIFIED_AGAINST.items():
            assert entry["source"], key
            assert entry["retrieved"], key
            assert entry["note"], key

    def test_imgt_boundaries_cite_imgt_itself(self):
        assert "imgt.org" in VERIFIED_AGAINST["imgt_region_boundaries"]["source"]

    def test_the_fr4_boundary_discrepancy_is_disclosed(self):
        """IMGT's own pages disagree (118-128 vs 118-129). Say which was used."""
        note = VERIFIED_AGAINST["imgt_region_boundaries"]["note"]
        assert "118-129" in note and "128 is used here" in note

    def test_the_structural_gap_rule_admits_it_is_not_shown_to_generalise(self):
        note = VERIFIED_AGAINST["structural_framework_gaps"]["note"]
        assert "not shown to generalise" in note

    def test_the_independent_oas_check_is_recorded_with_its_numbers(self):
        note = VERIFIED_AGAINST["numbering_vs_oas"]["note"]
        assert "99.908" in note
        assert "15,506" in note


class TestKnownSensitivities:
    def test_each_sensitivity_is_quantified_not_merely_asserted(self):
        for key, entry in KNOWN_SENSITIVITIES.items():
            assert entry["measured"], key
            assert entry["residual"], key

    def test_the_germline_set_sensitivity_carries_its_measured_numbers(self):
        entry = KNOWN_SENSITIVITIES["germline_set_choice"]
        assert "7.67%" in entry["measured"]
        assert "0.26%" in entry["residual"]
        assert "ambiguous" in entry["mitigated"]

    def test_the_threshold_mismatch_with_published_work_is_disclosed(self):
        entry = KNOWN_SENSITIVITIES["oasis_threshold"]
        assert "relaxed" in entry["measured"]
        assert "medium" in entry["residual"]
        assert "not comparable" in entry["residual"]

    def test_the_reference_set_quality_problem_is_disclosed(self):
        entry = KNOWN_SENSITIVITIES["oas_reference_quality"]
        assert "artefact" in entry["residual"]
        assert "Not quantified" in entry["residual"]


class TestPanelDigests:
    def test_the_reference_panel_is_hashed_so_a_run_can_be_tied_to_it(self):
        digests = panel_digests()
        assert digests["imgt_mab_db.fa"]
        assert all(len(d) == 64 for d in digests.values() if d)

    def test_a_missing_file_gives_none_rather_than_raising(self, tmp_path):
        assert file_digest(tmp_path / "nope.txt") is None


class TestWriteRunJson:
    def test_writes_valid_json_and_returns_the_path(self, tmp_path):
        path = write_run_json(tmp_path / "out")
        assert path.name == "run.json"
        json.loads(path.read_text())

    def test_creates_the_directory(self, tmp_path):
        path = write_run_json(tmp_path / "a" / "b")
        assert path.exists()

    def test_caller_inputs_are_merged_at_the_top_level(self, tmp_path):
        path = write_run_json(tmp_path, inputs={"sequences": 3, "source": "cli"})
        record = json.loads(path.read_text())
        assert record["inputs"] == {"sequences": 3, "source": "cli"}

    def test_the_record_is_serialisable_without_coercion(self):
        """Guards against a frozenset or Path leaking in and breaking json.dumps."""
        json.dumps(run_metadata())

    @pytest.mark.parametrize(
        "key",
        [
            "run_utc",
            "environment",
            "tools",
            "germline",
            "settings",
            "reference_panel",
            "verified_against",
            "known_sensitivities",
        ],
    )
    def test_required_sections_are_present(self, key):
        assert key in run_metadata()

    def test_timestamp_is_utc(self):
        assert run_metadata()["run_utc"].endswith("+00:00")


class TestEveryClaimStatesWhetherItCanBeReproduced:
    """The rule this file enforces: no number without its status.

    Several figures recorded here were measured during development against
    inputs the repository does not ship -- OAS publishes no licence for its
    sequence data, so it was read locally and is not redistributed. Those
    numbers are true and dated; a reader cannot rerun them from a checkout.
    Presenting one as a repository result is the single failure mode provenance
    exists to prevent, so each entry says which it is.
    """

    def test_every_verified_convention_declares_its_status(self):
        for name, entry in VERIFIED_AGAINST.items():
            assert "reproducible_from_repo" in entry, f"{name} does not say"
            assert entry["reproducible_from_repo"].strip()

    def test_every_sensitivity_declares_its_status(self):
        for name, entry in KNOWN_SENSITIVITIES.items():
            assert "reproducible_from_repo" in entry, f"{name} does not say"
            assert entry["reproducible_from_repo"].strip()

    def test_the_oas_cross_check_is_marked_as_not_reproducible_here(self):
        """It is the strongest number in the file and the least reproducible."""
        status = VERIFIED_AGAINST["numbering_vs_oas"]["reproducible_from_repo"]
        assert status.lower().startswith("no")
        assert "licence" in status or "license" in status

    def test_the_status_reaches_run_json(self):
        record = run_metadata()
        for entry in record["verified_against"].values():
            assert "reproducible_from_repo" in entry
