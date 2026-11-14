"""Streamlit front end for abscore.

Run with:  uv run streamlit run app.py

Two decisions worth stating, because they are what keep this usable outside a
demo:

The scores that matter run **offline** -- germline identity, repertoire
humanness and germline assignment need no network. Only the epitope section calls
IEDB, and it is requested explicitly instead of on every keystroke, so a slow
or unavailable public service degrades one panel instead of the page.

Nothing is shown without its basis. Every figure is accompanied by its
denominator or reference set, and every caveat the pipeline raised is displayed
above the scores instead of below them -- a partial domain or an
ambiguous germline call changes how the numbers should be read, so it is not
a footnote.
"""

from __future__ import annotations

import streamlit as st

from abscore.epitopes import SELF_DEFINITIONS
from abscore.iedb import DEFAULT_ALLELES, DEFAULT_METHOD, METHODS, IEDBError
from abscore.number import NumberingError
from abscore.report import analyse, escape_markdown, extraction_note, render_markdown


def parse_input(text: str) -> list[tuple[str, str]]:
    """Split pasted text into (name, sequence) pairs.

    Accepts a bare sequence, or multi-FASTA. The same parsing the CLI uses, so
    a batch pasted here and a batch piped to `abscore score` behave identically.
    """
    text = text.strip()
    if not text:
        return []
    if ">" not in text:
        return [("query", "".join(text.split()))]

    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].strip() or f"seq{len(records) + 1}", []
        else:
            parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


EXAMPLES = {
    "paste your own": "",
    "Trastuzumab VH (humanized)": (
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
        "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    ),
    "Adalimumab VH (human)": (
        "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEG"
        "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS"
    ),
    "Abagovomab VH (murine)": (
        "QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKG"
        "KATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS"
    ),
    "Trastuzumab VL (kappa)": (
        "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSR"
        "SGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
    ),
    "Batch of four (FASTA)": (
        ">adalimumab_VH (human)\n"
        "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEG"
        "RFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS\n"
        ">trastuzumab_VH (humanized)\n"
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
        "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS\n"
        ">rituximab_VH (chimeric)\n"
        "QVQLQQPGAELVKPGASVKMSCKASGYTFTSYNMHWVKQTPGRGLEWIGAIYPGNGDTSYNQKFKG"
        "KATLTADKSSSTAYMQLSSLTSEDSAVYYCARSTYYGGDWYFNVWGAGTTVTVSA\n"
        ">abagovomab_VH (murine)\n"
        "QVKLQESGAELARPGASVKLSCKASGYTFTNYWMQWVKQRPGQGLDWIGAIYPGDGNTRYTHKFKG"
        "KATLTADKSSSTAYMQLSSLASEDSGVYYCARGEGNYAWFAYWGQGTTVTVSS"
    ),
}

st.set_page_config(page_title="abscore", layout="wide")


@st.cache_data(show_spinner=False)
def _run(sequence: str, with_epitopes: bool, method: str, alleles: tuple[str, ...]):
    """Cached so changing a display option does not re-run the pipeline."""
    return analyse(
        sequence,
        with_epitopes=with_epitopes,
        method=method,
        alleles=alleles,
    )


def _pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def _percent(value: float | None) -> float | None:
    """A 0-1 fraction as a number of percent, for a table column."""
    return None if value is None else value * 100


st.title("abscore")
st.caption(
    "Immunogenicity triage for an antibody heavy or light chain: germline "
    "identity, repertoire humanness, and MHC class II predicted binders "
    "unmatched in a human reference."
)

# --- epitope settings ---------------------------------------------------------
# Read from session_state here, rendered under the epitope heading far below.
# A Streamlit script runs top to bottom, so these values are needed before the
# section they belong to is drawn. Keyed widgets put the answer in session_state
# and changing one reruns the script, so the read above and the widget below
# stay in step; on the first run the keys are absent and these defaults apply.
with_epitopes = st.session_state.get("with_epitopes", False)
self_definition = st.session_state.get("self_definition", SELF_DEFINITIONS[0])
method = st.session_state.get("method", DEFAULT_METHOD)
alleles = st.session_state.get("alleles", list(DEFAULT_ALLELES))

choice = st.selectbox("Example", list(EXAMPLES))
uploaded = st.file_uploader(
    "Or upload a FASTA file", type=("fa", "fasta", "txt"), label_visibility="collapsed"
)
default = EXAMPLES[choice]
if uploaded is not None:
    default = uploaded.getvalue().decode("utf-8", errors="replace")

text = st.text_area(
    "Sequence, or FASTA for a batch",
    value=default,
    height=160,
    placeholder=(
        "Paste one heavy or light chain variable domain, or several in FASTA "
        "format. Whitespace and lowercase are fine."
    ),
)

records = parse_input(text)
if not records:
    st.info("Paste a sequence, pick an example, or upload a FASTA file to begin.")
    st.stop()

if with_epitopes and len(records) > 1:
    st.caption(
        f"{len(records)} chains x {len(alleles) or len(DEFAULT_ALLELES)} alleles. "
        "Uncached predictions take roughly 20 seconds per chain."
    )

# --- score every record, keeping failures visible ---------------------------
results: list[tuple[str, object]] = []
failures: list[tuple[str, str]] = []
bar = st.progress(0.0, text="Scoring…") if len(records) > 1 else None
for index, (name, sequence) in enumerate(records, start=1):
    if bar:
        bar.progress((index - 1) / len(records), text=f"Scoring {name} ({index}/{len(records)})")
    try:
        results.append(
            (name, _run(sequence, with_epitopes, method, tuple(alleles) or DEFAULT_ALLELES))
        )
    except NumberingError as exc:
        failures.append((name, str(exc)))
    except IEDBError as exc:  # pragma: no cover - analyse() captures these
        failures.append((name, str(exc)))
if bar:
    bar.empty()

if failures:
    with st.container(border=True):
        st.subheader(f"{len(failures)} sequence(s) could not be scored")
        for name, reason in failures:
            st.error(f"**{name}**: {reason}", icon=":material/close:")

if not results:
    st.stop()

# --- batch: a sortable table, then a chosen chain in full -------------------
if len(results) > 1:
    st.subheader(f"{len(results)} chains")
    st.dataframe(
        [
            {
                "Name": name,
                "Chain": r.chain_type,
                "Germline V": r.germline.gene,
                # Scaled to percent in the row, because Streamlit's
                # NumberColumn format is printf, so "%.1f%%" on a 0-1 fraction
                # renders 0.85 as "0.9%". Sorting is unaffected either way.
                "V identity": _percent(r.germline.identity),
                "FR identity": _percent(r.germline.framework_identity),
                "OASis": _percent(r.humanness.score),
                "Hits": (r.epitopes.count() if r.epitopes else None),
                "Cores": (len(r.epitopes.unique_cores(None)) if r.epitopes else None),
                f"Unmatched ({self_definition})": (
                    r.epitopes.count(self_definition) if r.epitopes else None
                ),
                "Caveats": len(r.warnings()),
            }
            for name, r in results
        ],
        hide_index=True,
        use_container_width=True,
        column_config={
            "V identity": st.column_config.NumberColumn(format="%.1f%%"),
            "FR identity": st.column_config.NumberColumn(format="%.1f%%"),
            "OASis": st.column_config.NumberColumn(format="%.1f%%"),
        },
    )
    st.caption(
        "Click a column header to sort. Identities are percentages; see a "
        "single chain below for the denominator behind each one."
    )

    st.download_button(
        "Download all reports (markdown)",
        data="\n\n---\n\n".join(
            render_markdown(r, title=name, self_definition=self_definition) for name, r in results
        ),
        file_name="abscore_reports.md",
        mime="text/markdown",
    )

    st.divider()
    chosen = st.selectbox(
        "Show one chain in full",
        [name for name, _ in results],
        help="Every number below carries its basis; the table above is a summary.",
    )
    report = dict(results)[chosen]
    st.subheader(chosen)
else:
    report = results[0][1]

domain = report.domain
germline = report.germline

# --- caveats first: they change how the numbers should be read ---------------
warnings = report.warnings()
if warnings:
    with st.container(border=True):
        st.subheader(f"{len(warnings)} caveat(s) on these numbers")
        for message in warnings:
            # Allele names contain '*', which markdown reads as emphasis.
            st.warning(escape_markdown(message), icon="⚠️")
else:
    st.success(
        "Complete, unambiguous V-domain. No sequence-level caveats.",
        icon=":material/check_circle:",
    )

# --- the chain ---------------------------------------------------------------
left, right = st.columns([2, 3])
with left:
    st.subheader("Chain")
    st.write(f"**Type** {report.chain_name} ({domain.chain_type})")
    span = f"{domain.imgt_span[0]}-{domain.imgt_span[1]}"
    st.write(f"**Numbered** {len(domain.numbering)} residues, IMGT {span}")
    extracted = extraction_note(domain)
    if extracted:
        st.write(f"**Scored on** {extracted}")
    st.write(f"**Closest germline V** `{germline.gene}`")
    if germline.ambiguous_genes:
        st.write(
            "**Indistinguishable from** "
            + ", ".join(f"`{c.best_allele}`" for c in germline.ambiguous_genes)
        )
with right:
    st.subheader("Regions")
    from abscore.number import IMGT_REGIONS, STRUCTURAL_FRAMEWORK_GAPS
    from abscore.track import render_track_html

    st.html(render_track_html(domain))

    # The table view stays available: a picture is not an accessible substitute
    # for the values, and a reader may want to copy a region's sequence.
    with st.expander("Regions as a table"):
        st.dataframe(
            [
                {
                    "Region": name,
                    "IMGT slots": f"{low}-{high}",
                    "Filled": f"{len(domain.regions.get(name, ''))}/{high - low + 1}",
                    "Sequence": domain.regions.get(name, ""),
                }
                for name, (low, high) in IMGT_REGIONS.items()
            ],
            hide_index=True,
            use_container_width=True,
        )
        gaps = ", ".join(
            str(p) for p in sorted(STRUCTURAL_FRAMEWORK_GAPS.get(domain.chain_type, ()))
        )
        st.caption(
            f"Framework slots {gaps} (and the last slot, 128) are empty in nearly "
            "every chain of this type, so no deletion is reported. These are IMGT "
            "boundaries: Kabat puts four more residues in CDR-L2 that appear here "
            "at the start of FR3."
        )

st.divider()

# --- humanness ---------------------------------------------------------------
st.subheader("Humanness")
scope_v = germline.scopes["V"]
scope_fr = germline.scopes["framework"]
scope_cdr = germline.scopes["cdr"]

cols = st.columns(4)
cols[0].metric("Germline identity, V", _pct(scope_v.identity))
cols[0].caption(f"{scope_v.matches}/{scope_v.compared} positions vs {germline.gene}")
cols[1].metric("Germline identity, framework", _pct(scope_fr.identity))
cols[1].caption("FR1+FR2+FR3 only; FR4 is J-encoded")
cols[2].metric("Germline identity, CDR1+2", _pct(scope_cdr.identity))
cols[2].caption("CDR3 excluded; junctional in origin")
cols[3].metric("OASis identity", _pct(report.humanness.score))
cols[3].caption(
    f"{report.humanness.n_found}/{report.humanness.n_peptides} 9-mers in "
    f"≥{report.humanness.threshold_min_donor_fraction:.0%} of OAS donors"
)

if not scope_v.is_full_coverage and scope_v.coverage is not None:
    st.warning(
        f"Coverage {scope_v.coverage:.1%}: {scope_v.compared} of {scope_v.covered} "
        "germline positions observed. Identity is not comparable with "
        "full-coverage values.",
        icon="⚠️",
    )

runs = report.humanness.flagged_runs()
if runs:
    with st.expander(f"{len(runs)} stretch(es) not in the OAS reference"):
        st.dataframe(
            [
                {
                    "IMGT span": f"{s[0]}{s[1].strip()}-{e[0]}{e[1].strip()}",
                    "Positions flagged": n,
                    "Regions": ", ".join(r),
                }
                for s, e, n, r in sorted(runs, key=lambda r: -r[2])
            ],
            hide_index=True,
            use_container_width=True,
        )
        st.caption("An observation. The cause stays open.")

st.divider()

# --- epitopes ----------------------------------------------------------------
st.subheader("MHC class II epitope burden")

with st.container(border=True):
    st.checkbox(
        "Predict MHC class II epitopes",
        key="with_epitopes",
        help=(
            "Calls the IEDB Analysis Resource: about 20 seconds for 8 alleles, "
            "then cached. Every other number on this page is computed offline."
        ),
    )
    settings = st.columns(2)
    with settings[0]:
        st.selectbox(
            "Reference set",
            SELF_DEFINITIONS,
            key="self_definition",
            help=(
                "Which human reference a predicted core is looked up in. "
                "germline: human germline Ig V and J alleles. "
                "repertoire: human antibody repertoires (OAS). "
                "proteome: anywhere in the human proteome."
            ),
        )
    with settings[1]:
        st.selectbox(
            "Predictor",
            METHODS,
            key="method",
            help=(
                "Run on IEDB's servers; nothing is downloaded. Only the default "
                "is pre-cached, so changing it costs a live call."
            ),
        )
    st.multiselect(
        "HLA-DRB1 alleles",
        DEFAULT_ALLELES,
        default=list(DEFAULT_ALLELES),
        key="alleles",
        help=(
            "MHC II grooves differ between people, so a peptide presented to one "
            "patient may not be presented to another. These eight are the panel "
            "from Prihoda et al., MAbs 2022. The count sums across all eight, so "
            "it describes a population rather than the molecule alone."
        ),
    )

if report.epitopes is None:
    if report.epitope_error:
        st.error(
            "IEDB was not reachable. Every other number on this page is "
            "computed offline and is unaffected.",
            icon="🔌",
        )
        st.caption(report.epitope_error)
    else:
        st.info("Not requested. Enable **Predict MHC class II epitopes** in the sidebar.")
else:
    burden = report.epitopes
    counts = burden.counts_by_definition()
    st.caption(
        f"{len(burden.alleles)} HLA-DRB1 alleles · `{burden.method}` · "
        f"{burden.peptide_length}-mers · percentile rank ≤ {burden.binder_rank:g}"
    )
    cols = st.columns(4)
    cols[0].metric("Predicted hits", counts["raw"])
    cols[0].caption("one per peptide-allele pair")
    cols[1].metric("Distinct cores", counts["distinct_cores"])
    cols[1].caption("one 9-mer, however many alleles present it")
    cols[2].metric(f"Unmatched ({self_definition})", counts[f"{self_definition}_unmatched"])
    cols[2].caption(f"no exact match in the {self_definition} reference")
    cols[3].metric("Matched", counts[f"{self_definition}_matched"])
    cols[3].caption("in the reference; lower priority, still on the list")

    st.caption(
        f"Of the {counts[f'{self_definition}_unmatched']} unmatched: "
        f"{counts[f'{self_definition}_one_substitution']} are one substitution "
        f"from the reference, {counts[f'{self_definition}_distant']} are 2+ away. "
        "Distance is an annotation; the effect of a substitution turns on "
        "which position changed."
    )

    unclassifiable = counts[f"{self_definition}_unclassifiable"]
    if unclassifiable:
        st.caption(
            f"{unclassifiable} junctional peptide(s) unclassifiable by "
            f"{self_definition}; abstained on."
        )

    unmatched = burden.unmatched(self_definition)
    if unmatched:
        st.dataframe(
            [
                {
                    "Rank": e.rank,
                    "Allele": e.allele,
                    "Peptide": e.peptide,
                    "Core": e.core,
                    "Region": ", ".join(e.regions),
                    "Distance": e.distance_label(self_definition),
                }
                for e in unmatched[:20]
            ],
            hide_index=True,
            use_container_width=True,
        )
        st.caption(
            (
                "Lower rank = more likely presented (an eluted-ligand score; "
                "`netmhciipan_ba` gives affinity)."
                if burden.method.endswith("_el")
                else "Lower rank = tighter predicted binding."
            )
            + " Top 20 shown."
        )

st.divider()

with st.expander("Full markdown report"):
    st.code(
        render_markdown(report, self_definition=self_definition),
        language="markdown",
    )

with st.expander("Provenance: what produced every number"):
    st.json(report.provenance())

st.caption(
    "Computed from sequence alone. Use it as triage; the full report states "
    "what clinical immunogenicity also depends on."
)
