"""Tests for the Streamlit front end, by executing it against a fake Streamlit.

Why this exists: `app.py` is a script rather than a set of functions, so importing it
proves nothing and starting the server proves only that the server starts. The
page's branches -- a batch table, a failed record, an unavailable service, a
chain with no caveats -- are reached by *running the script*, and a bad f-string
or dict key in any of them is a crash the user sees and the test suite does not.

The stub below records what the page asked to render and answers its widgets
from a scripted dict. It deliberately enforces one piece of real Streamlit
behaviour: `icon=` is validated, because Streamlit rejects anything that is not
a single emoji or a `:material/name:`, and a stub that accepts every icon hides
exactly the crash it should catch. A dingbat check mark ("✓") shipped this
way once -- reachable only from a chain with zero caveats, which is why it
survived every other test.

No network: epitope scenarios patch the prediction call to read the committed
responses under `tests/data/iedb`.
"""

from __future__ import annotations

import contextlib
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest
from streamlit.string_util import validate_icon_or_emoji

import abscore.epitopes as epitopes_module
from abscore.iedb import IEDBError, predict_mhcii
from tests.fixtures import (
    ABAGOVOMAB_VH,
    ADALIMUMAB_VH,
    ALBUMIN_FRAGMENT,
    RITUXIMAB_VH,
    TRASTUZUMAB_VH,
    TRASTUZUMAB_VL,
    VH_WITH_TAIL,
)

APP = Path(__file__).resolve().parent.parent / "app.py"
CACHE = "tests/data/iedb"

# Streamlit's icon validator imports `streamlit.emojis` and
# `streamlit.material_icon_names` lazily, and the stub replaces
# sys.modules["streamlit"], which would break those imports. Warm them here.
for _warm in ("⚠️", ":material/check_circle:"):
    validate_icon_or_emoji(_warm)

# Widget labels, as the app spells them.
SEQ = "Sequence, or FASTA for a batch"
UPLOAD = "Or upload a FASTA file"
# Widget *keys* rather than labels: these controls live under the epitope heading and
# are read from session_state before that section is drawn.
EPITOPES = "with_epitopes"
SELF = "self_definition"
ALLELES = "alleles"

OUTPUT_CALLS = frozenset(
    {
        "title",
        "caption",
        "subheader",
        "header",
        "write",
        "markdown",
        "info",
        "warning",
        "error",
        "success",
        "code",
        "json",
        "dataframe",
        "divider",
        "html",
        "download_button",
        "set_page_config",
        "metric",
        "table",
        "text",
    }
)


class Stop(Exception):
    """Raised for st.stop(), which ends the script rather than the process."""


class FakeStreamlit:
    """Records rendered output; answers widgets from `answers`."""

    def __init__(self, answers: dict):
        self.answers = answers
        self.calls: list[tuple[str, object]] = []
        self.sidebar = self
        self.column_config = _ColumnConfig()
        # The epitope controls are keyed widgets rendered below the results but
        # read above them, via session_state. The same answers dict serves both:
        # entries keyed by widget label answer the widget, entries keyed by
        # widget key answer the session_state read.
        self.session_state = dict(answers)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def __getattr__(self, name):
        if name in OUTPUT_CALLS:
            return self._record(name)
        raise AttributeError(f"app.py used an unstubbed Streamlit API: st.{name}")

    def _record(self, kind):
        def fn(*args, **kwargs):
            if kwargs.get("icon") is not None:
                validate_icon_or_emoji(kwargs["icon"])
            self.calls.append((kind, args[0] if args else None))

        return fn

    # --- widgets ---
    def _answer(self, label, default):
        return self.answers.get(label, default)

    def checkbox(self, label, value=False, **kw):
        return self._answer(label, value)

    def selectbox(self, label, options, index=0, **kw):
        options = list(options)
        return self._answer(label, options[index] if options else None)

    def multiselect(self, label, options, default=None, **kw):
        return self._answer(label, list(default if default is not None else []))

    def text_area(self, label, value="", **kw):
        return self._answer(label, value)

    def file_uploader(self, label, **kw):
        return self._answer(label, None)

    # --- layout ---
    def columns(self, spec, **kw):
        return [self] * (spec if isinstance(spec, int) else len(spec))

    @contextmanager
    def container(self, **kw):
        yield self

    @contextmanager
    def expander(self, label, **kw):
        self.calls.append(("expander", label))
        yield self

    def progress(self, value, text=None):
        return _Bar()

    def stop(self):
        raise Stop

    def cache_data(self, fn=None, **kw):
        return (lambda f: f) if fn is None else fn

    # --- what a test asserts against ---
    @property
    def page(self) -> str:
        return "\n".join(str(arg) for _, arg in self.calls if isinstance(arg, str))

    def of_kind(self, kind: str) -> list:
        return [arg for k, arg in self.calls if k == kind]

    @property
    def tables(self) -> list:
        return [arg for arg in self.of_kind("dataframe") if isinstance(arg, list)]


class _Bar:
    def progress(self, *a, **k): ...

    def empty(self): ...


class _ColumnConfig:
    def NumberColumn(self, *a, **k):  # noqa: N802 - mirrors Streamlit's API
        return {"format": k.get("format")}


def render(answers: dict, monkeypatch) -> FakeStreamlit:
    """Execute app.py with these widget answers, and return what it rendered."""
    stub = FakeStreamlit(answers)
    monkeypatch.setitem(sys.modules, "streamlit", stub)
    with contextlib.suppress(Stop):
        exec(compile(APP.read_text(), str(APP), "exec"), {"__name__": "__main__"})  # noqa: S102
    return stub


@pytest.fixture
def cached_epitopes(monkeypatch):
    """Serve epitope predictions from the committed responses, never the network."""

    def from_cache(sequence, **kwargs):
        kwargs["cache_dir"] = CACHE
        kwargs["offline"] = True
        return predict_mhcii(sequence, **kwargs)

    monkeypatch.setattr(epitopes_module, "predict_mhcii", from_cache)


class TestEveryBranchRenders:
    """Each of these reaches a different part of the page."""

    def test_a_single_chain(self, monkeypatch):
        stub = render({SEQ: TRASTUZUMAB_VH}, monkeypatch)
        assert "Humanness" in stub.page
        assert "Chain" in stub.page

    def test_a_chain_with_no_caveats_reaches_the_success_branch(self, monkeypatch):
        """The branch that a dingbat icon used to crash. Adalimumab is the only
        fixture with zero caveats, so it is the only one that gets here."""
        stub = render({SEQ: ADALIMUMAB_VH}, monkeypatch)
        assert stub.of_kind("success"), "expected the no-caveats branch to render"
        assert not stub.of_kind("warning")

    def test_a_chain_with_caveats_reaches_the_warning_branch(self, monkeypatch):
        """Rituximab's heavy chain ties two germline genes, so it has a caveat.

        Abagovomab was used here until the ambiguity rule was tightened to
        exact ties; it now has none, which is correct -- its runner-up gene is a
        position behind, so nothing ties.
        """
        stub = render({SEQ: RITUXIMAB_VH}, monkeypatch)
        assert stub.of_kind("warning")
        assert not stub.of_kind("success")

    def test_a_batch(self, monkeypatch):
        text = f">human\n{ADALIMUMAB_VH}\n>murine\n{ABAGOVOMAB_VH}\n"
        stub = render({SEQ: text}, monkeypatch)
        assert "2 chains" in stub.page
        assert [row["Name"] for row in stub.tables[0]] == ["human", "murine"]

    def test_epitopes(self, monkeypatch, cached_epitopes):
        stub = render({SEQ: TRASTUZUMAB_VH, EPITOPES: True}, monkeypatch)
        assert "MHC class II epitope burden" in stub.page

    @pytest.mark.parametrize("definition", ["germline", "repertoire", "proteome"])
    def test_each_self_definition(self, definition, monkeypatch, cached_epitopes):
        stub = render({SEQ: TRASTUZUMAB_VH, EPITOPES: True, SELF: definition}, monkeypatch)
        assert f"Unmatched ({definition})" in stub.page

    def test_an_unavailable_iedb_degrades_one_panel_only(self, monkeypatch):
        """The page's central robustness claim, and it is testable."""

        def dead(*a, **k):
            raise IEDBError("simulated outage")

        monkeypatch.setattr(epitopes_module, "predict_mhcii", dead)
        stub = render({SEQ: TRASTUZUMAB_VH, EPITOPES: True}, monkeypatch)
        assert any("IEDB was not reachable" in str(e) for e in stub.of_kind("error"))
        # Everything computed offline still rendered.
        assert "Humanness" in stub.page
        assert "Regions" in stub.page


class TestInputsTheUserWillActuallyPaste:
    def test_nothing_pasted_prompts_rather_than_fails(self, monkeypatch):
        stub = render({SEQ: "   \n  "}, monkeypatch)
        assert "Paste a sequence" in stub.page

    def test_a_non_antibody_is_rejected_with_a_reason(self, monkeypatch):
        stub = render({SEQ: ALBUMIN_FRAGMENT}, monkeypatch)
        assert any("no antibody V-domain" in str(e) for e in stub.of_kind("error"))

    def test_one_bad_record_does_not_lose_the_good_ones(self, monkeypatch):
        text = f">good\n{TRASTUZUMAB_VH}\n>bad\n{ALBUMIN_FRAGMENT}\n"
        stub = render({SEQ: text}, monkeypatch)
        assert stub.of_kind("error")
        assert "Humanness" in stub.page

    def test_two_domains_in_one_record_say_so(self, monkeypatch):
        scfv = TRASTUZUMAB_VH + "GGGGSGGGGSGGGGS" + TRASTUZUMAB_VL
        stub = render({SEQ: scfv}, monkeypatch)
        assert any("2 V-domains" in str(e) for e in stub.of_kind("error"))

    def test_a_constant_region_tail_is_scored_on_the_v_domain_alone(self, monkeypatch):
        """Pasting a full-length chain is fine: only the V-domain is numbered.

        This asserted a caveat until the ambiguity rule was tightened, and it
        passed for the wrong reason -- the caveat came from the V-domain's gene
        call. The trailing residues raise no warning at all, which this states
        plainly.
        """
        stub = render({SEQ: VH_WITH_TAIL}, monkeypatch)
        assert "Humanness" in stub.page
        assert stub.of_kind("success"), "a clean V-domain, tail or no tail"

    def test_lowercase_and_whitespace_are_accepted(self, monkeypatch):
        messy = TRASTUZUMAB_VH[:40].lower() + "\n " + TRASTUZUMAB_VH[40:]
        stub = render({SEQ: messy}, monkeypatch)
        assert "Humanness" in stub.page

    def test_an_uploaded_file_is_used(self, monkeypatch):
        class Upload:
            def getvalue(self):
                return f">uploaded\n{ADALIMUMAB_VH}\n".encode()

        stub = render({UPLOAD: Upload()}, monkeypatch)
        assert "Humanness" in stub.page

    def test_deselecting_every_allele_falls_back_to_the_panel(self, monkeypatch):
        stub = render({SEQ: TRASTUZUMAB_VH, ALLELES: []}, monkeypatch)
        assert "Humanness" in stub.page


class TestTheBatchTableReportsPercentages:
    """Streamlit's NumberColumn format is printf, so a 0-1 fraction under
    "%.1f%%" renders 0.85 as "0.9%", which is a wrong number."""

    def test_identities_are_scaled_to_percent(self, monkeypatch):
        text = f">human\n{ADALIMUMAB_VH}\n>murine\n{ABAGOVOMAB_VH}\n"
        rows = render({SEQ: text}, monkeypatch).tables[0]
        for row in rows:
            for column in ("V identity", "FR identity", "OASis"):
                assert 1.0 < row[column] <= 100.0, f"{column} looks like a fraction"

    def test_the_more_human_chain_scores_higher(self, monkeypatch):
        """Not a validation of the scores -- a check that the table is not
        transposing or mislabelling its rows."""
        text = f">human\n{ADALIMUMAB_VH}\n>murine\n{ABAGOVOMAB_VH}\n"
        human, murine = render({SEQ: text}, monkeypatch).tables[0]
        assert human["V identity"] > murine["V identity"]
        assert human["OASis"] > murine["OASis"]
