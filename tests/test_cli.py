"""Tests for the command line interface.

Offline throughout: every invocation passes --offline with the committed cache,
so nothing here touches the network.
"""

import json

import pytest

from abscore.cli import build_parser, main

CACHE = "tests/data/iedb"
TRASTUZUMAB_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)
ABAGOVOMAB_VH = (
    "QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKG"
    "KATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS"
)
ALBUMIN = "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKS"


def run(*argv):
    return main(["score", "--cache-dir", CACHE, "--offline", "--quiet", *argv])


class TestExitCodes:
    def test_a_scored_chain_exits_zero(self):
        assert run("--seq", TRASTUZUMAB_VH) == 0

    def test_a_non_antibody_exits_two(self, capsys):
        """Bad input is the caller's problem and must be distinguishable."""
        assert run("--seq", ALBUMIN, "--no-epitopes") == 2
        assert "no antibody V-domain" in capsys.readouterr().err

    def test_a_missing_system_dependency_exits_three_with_instructions(self, monkeypatch):
        """The likeliest failure on a machine that has never run this.

        HMMER is a system dependency; `uv sync` cannot install it. ANARCI
        shells out to `hmmscan`, so without it nothing numbers. This used to
        reach the terminal as an unhandled traceback with the useful message
        buried in the last line -- the worst presentation for the one error a
        first-time user is most likely to hit, and the one with a one-line fix.

        Exit 3, not 2: nothing is wrong with what the caller passed.
        """
        from abscore import cli
        from abscore.number import HmmerNotFoundError

        def missing(*args, **kwargs):
            raise HmmerNotFoundError()

        monkeypatch.setattr(cli, "cmd_score", missing)
        assert run("--seq", TRASTUZUMAB_VH) == 3

    def test_the_instructions_name_every_platform(self, capsys, monkeypatch):
        from abscore import cli
        from abscore.number import HmmerNotFoundError

        monkeypatch.setattr(
            cli, "cmd_score", lambda *a, **k: (_ for _ in ()).throw(HmmerNotFoundError())
        )
        run("--seq", TRASTUZUMAB_VH)
        err = capsys.readouterr().err
        for platform in ("brew install hmmer", "apt-get install hmmer", "bioconda"):
            assert platform in err

    def test_an_unavailable_service_does_not_fail_the_run(self):
        """IEDB being down is not the caller's fault; the report still stands."""
        assert run("--seq", ABAGOVOMAB_VH.replace("QVKLQ", "EVQLV")) == 0


class TestOutputStreams:
    def test_results_go_to_stdout_and_progress_to_stderr(self, capsys):
        """Checked on the DEFAULT format, since that is the path a reader hits:
        the report on stdout, the running log on stderr."""
        main(["score", "--seq", TRASTUZUMAB_VH, "--cache-dir", CACHE, "--offline"])
        captured = capsys.readouterr()
        assert "CHAIN" in captured.out
        assert "abscore:" in captured.err
        # So `abscore score ... > file` writes a clean file.
        assert "abscore:" not in captured.out

    def test_quiet_suppresses_progress_only(self, capsys):
        run("--seq", TRASTUZUMAB_VH)
        captured = capsys.readouterr()
        assert "CHAIN" in captured.out
        assert "abscore:" not in captured.err


class TestFormats:
    def test_summary_is_one_line(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "summary")
        assert len(capsys.readouterr().out.strip().splitlines()) == 1

    def test_markdown_renders_the_full_report(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "markdown")
        out = capsys.readouterr().out
        assert "## Humanness" in out

    def test_json_is_valid_and_carries_provenance(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "json")
        payload = json.loads(capsys.readouterr().out)
        assert payload["germline"]["gene"].startswith("IGHV")
        assert payload["humanness"]["oasis_identity"] > 0
        assert payload["provenance"]["humanness"]["promb_version"]

    def test_json_keeps_the_basis_beside_every_score(self, capsys):
        """A number without its denominator is not reportable."""
        run("--seq", TRASTUZUMAB_VH, "--format", "json")
        payload = json.loads(capsys.readouterr().out)
        scope = payload["germline"]["scopes"]["V"]
        assert {"matches", "compared", "germline_positions", "coverage"} <= set(scope)
        humanness = payload["humanness"]
        assert {"peptides_found", "peptides_total", "threshold"} <= set(humanness)


class TestInputModes:
    def test_fasta_file_gives_one_report_per_record(self, tmp_path, capsys):
        path = tmp_path / "two.fa"
        path.write_text(f">a\n{TRASTUZUMAB_VH}\n>b\n{ABAGOVOMAB_VH}\n")
        assert run("--fasta", str(path), "--format", "summary") == 0
        lines = capsys.readouterr().out.strip().splitlines()
        assert len(lines) == 2
        assert lines[0].startswith("a\t")
        assert lines[1].startswith("b\t")

    def test_a_bare_sequence_file_without_headers_works(self, tmp_path, capsys):
        path = tmp_path / "one.txt"
        path.write_text(TRASTUZUMAB_VH)
        assert run("--fasta", str(path), "--format", "summary") == 0
        assert "chain=H" in capsys.readouterr().out

    def test_seq_and_fasta_together_are_rejected(self, tmp_path):
        path = tmp_path / "x.fa"
        path.write_text(f">a\n{TRASTUZUMAB_VH}\n")
        with pytest.raises(SystemExit):
            main(["score", "--seq", TRASTUZUMAB_VH, "--fasta", str(path)])

    def test_one_bad_record_does_not_lose_the_others(self, tmp_path, capsys):
        path = tmp_path / "mixed.fa"
        path.write_text(f">good\n{TRASTUZUMAB_VH}\n>bad\n{ALBUMIN}\n")
        assert run("--fasta", str(path), "--format", "summary") == 2  # reports the failure
        assert "good\t" in capsys.readouterr().out  # but still scores the good one


class TestOutDirectory:
    def test_markdown_and_run_json_are_written(self, tmp_path):
        path = tmp_path / "out"
        run("--seq", TRASTUZUMAB_VH, "--format", "markdown", "--out", str(path))
        assert (path / "query.md").exists()
        assert (path / "run.json").exists()
        json.loads((path / "run.json").read_text())

    def test_run_json_records_the_inputs_and_options(self, tmp_path):
        path = tmp_path / "out"
        run("--seq", TRASTUZUMAB_VH, "--format", "json", "--out", str(path))
        record = json.loads((path / "run.json").read_text())
        assert record["inputs"]["sequences"] == 1
        assert record["inputs"]["offline"] is True
        assert record["inputs"]["self_definition"] == "germline"

    def test_fasta_names_become_filenames(self, tmp_path):
        source = tmp_path / "two.fa"
        source.write_text(f">tras VH\n{TRASTUZUMAB_VH}\n>abago/VH\n{ABAGOVOMAB_VH}\n")
        out = tmp_path / "out"
        run("--fasta", str(source), "--format", "markdown", "--out", str(out))
        assert (out / "tras_VH.md").exists()
        assert (out / "abago_VH.md").exists()  # unsafe characters replaced


class TestOptions:
    def test_no_epitopes_skips_the_section(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--no-epitopes", "--format", "summary")
        assert "binders=unavailable" in capsys.readouterr().out

    def test_self_definition_is_selectable(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "json", "--self-definition", "repertoire")
        payload = json.loads(capsys.readouterr().out)
        assert payload["epitopes"]["reference_reported"] == "repertoire"

    def test_a_single_allele_can_be_requested(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "json", "--alleles", "HLA-DRB1*01:01")
        payload = json.loads(capsys.readouterr().out)
        assert list(payload["epitopes"]["by_allele"]) == ["HLA-DRB1*01:01"]

    def test_the_predictor_choices_match_what_iedb_accepts(self):
        from abscore.iedb import METHODS

        parser = build_parser()
        action = next(
            a
            for a in parser._subparsers._actions[-1].choices["score"]._actions
            if a.dest == "method"
        )
        assert tuple(action.choices) == METHODS


class TestProvenanceCommand:
    def test_it_prints_valid_json_with_the_sensitivities(self, capsys):
        assert main(["provenance"]) == 0
        record = json.loads(capsys.readouterr().out)
        assert record["tools"]["anarci"]
        assert "known_sensitivities" in record
        assert "reference_distance" in record["known_sensitivities"]


class TestAppParsingMatchesCli:
    """The app parses pasted input with the same rules as the CLI.

    A batch pasted into the app and the same batch piped to `abscore score`
    must split into the same records, or the two interfaces disagree about what
    the user gave them.
    """

    @staticmethod
    def _app_parse(text):
        source = __import__("pathlib").Path("app.py").read_text()
        namespace: dict = {}
        body = source[source.index("def parse_input") : source.index("EXAMPLES = {")]
        exec(body, {}, namespace)  # noqa: S102 - importing app.py needs streamlit
        return namespace["parse_input"](text)

    @staticmethod
    def _cli_parse(text, tmp_path):
        import argparse

        from abscore.cli import _read_sequences

        path = tmp_path / "in.fa"
        path.write_text(text)
        args = argparse.Namespace(seq=None, fasta=str(path))
        return _read_sequences(args)

    @pytest.mark.parametrize(
        "text",
        [
            TRASTUZUMAB_VH,
            f">a\n{TRASTUZUMAB_VH}\n>b\n{ABAGOVOMAB_VH}\n",
            f">a (human)\n{TRASTUZUMAB_VH[:60]}\n{TRASTUZUMAB_VH[60:]}\n",
            f">\n{TRASTUZUMAB_VH}\n",
        ],
    )
    def test_both_interfaces_split_input_identically(self, text, tmp_path):
        assert self._app_parse(text) == self._cli_parse(text, tmp_path)

    def test_a_bare_sequence_loses_its_whitespace(self):
        assert self._app_parse(" EVQLV\nESGGG ") == [("query", "EVQLVESGGG")]

    def test_empty_input_yields_no_records(self):
        assert self._app_parse("   \n  ") == []


class TestAppPercentScaling:
    """The batch table's identity columns are scaled to percent in the data.

    Streamlit's `NumberColumn(format=...)` is printf, so a 0-1 fraction under
    "%.1f%%" renders 0.85 as "0.9%", which is a wrong number. The
    scaling therefore happens in the row, and this test is what keeps it there.
    """

    @staticmethod
    def _app_percent():
        source = __import__("pathlib").Path("app.py").read_text()
        namespace: dict = {}
        body = source[source.index("def _percent") : source.index("st.title")]
        exec(body, {}, namespace)  # noqa: S102 - importing app.py needs streamlit
        return namespace["_percent"]

    def test_a_fraction_becomes_a_number_of_percent(self):
        assert self._app_percent()(0.85) == pytest.approx(85.0)

    def test_an_absent_score_stays_absent(self):
        """None must survive as None. Zero would read as a real measurement."""
        assert self._app_percent()(None) is None

    def test_the_batch_table_scales_every_column_it_formats_as_a_percentage(self):
        """Guards the pairing itself -- a future column formatted as a
        percentage but handed a raw fraction is the same bug again."""
        source = __import__("pathlib").Path("app.py").read_text()
        table = source[source.index("if len(results) > 1:") : source.index("Download all reports")]
        code = [line for line in table.splitlines() if not line.lstrip().startswith("#")]
        formatted = {line.split('"')[1] for line in code if "NumberColumn" in line and "%%" in line}
        scaled = {line.split('"')[1] for line in code if "_percent(" in line and '"' in line}
        assert formatted, "expected percentage-formatted columns in the batch table"
        assert formatted == scaled


class TestTextFormatAndSave:
    """The pairing the terminal view is built for: readable on screen, complete
    in the file. A text run and its saved markdown must not disagree."""

    def test_text_format_prints_the_report_without_markdown_syntax(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "text", "--no-epitopes")
        out = capsys.readouterr().out
        assert "CHAIN" in out and "HUMANNESS" in out
        assert "|" not in out and "**" not in out

    def test_save_writes_the_markdown_with_its_tables(self, capsys, tmp_path):
        path = tmp_path / "report.md"
        run("--seq", TRASTUZUMAB_VH, "--format", "text", "--no-epitopes", "--save", str(path))
        assert path.exists()
        written = path.read_text()
        assert written.startswith("#")
        assert "| Region | IMGT slots | Filled | Sequence |" in written
        # The terminal got the text form; the markdown went to the file.
        assert "| Region |" not in capsys.readouterr().out

    def test_both_views_carry_the_same_headline_numbers(self, capsys, tmp_path):
        path = tmp_path / "report.md"
        run("--seq", TRASTUZUMAB_VH, "--format", "text", "--no-epitopes", "--save", str(path))
        terminal = capsys.readouterr().out
        for value in ("IGHV3-66*04", "83.2%", "70.5%"):
            assert value in terminal, value
            assert value in path.read_text(), value

    def test_save_works_with_a_batch(self, capsys, tmp_path):
        path = tmp_path / "all.md"
        fasta = tmp_path / "two.fa"
        fasta.write_text(f">a\n{TRASTUZUMAB_VH}\n>b\n{TRASTUZUMAB_VH}\n")
        run("--fasta", str(fasta), "--format", "text", "--no-epitopes", "--save", str(path))
        written = path.read_text()
        assert written.count("# ") >= 2 or "---" in written
        assert "a" in capsys.readouterr().out


class TestTheNumberingHmmStaysMachineReadable:
    """Dropped from the human-readable output, kept in the JSON.

    The species ANARCI picked feeds no score -- every one hardcodes human -- so
    it misleads in a report. It is still what produced the coordinates, which is
    exactly what run.json is for.
    """

    def test_the_json_still_records_which_profile_numbered_the_chain(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "json", "--no-epitopes")
        payload = json.loads(capsys.readouterr().out)
        assert payload["numbering_hmm_species"]

    def test_the_text_report_does_not_mention_it(self, capsys):
        run("--seq", TRASTUZUMAB_VH, "--format", "text", "--no-epitopes")
        assert "Numbering HMM" not in capsys.readouterr().out
