"""The report as a self-contained HTML page.

Markdown was the only file format here for a while, and it has one structural
limit: it cannot draw the position track. A pipe table can say CDR2 fills 3 of
10 slots; only the track shows *which* 10 and where the gaps fall, which is the
clearest thing in the report. Raw markdown pipes are also hard to read unless
something renders them, and a reader opening `report.md` in an editor gets the
pipes.

So this writes one file that opens in any browser with no assets and no
network: styles inline, the track inline, light and dark both. Markdown stays
for the places that want it, such as a README or a PR comment.

Two decisions:

**Numbers lead, reasoning folds away.** Each explanatory passage sits in a
`<details>` block. The caveats and the counts are visible on open; the three
paragraphs explaining how to read a distance are one click away. In the
markdown report that prose sits between every table, which buries the numbers
it is meant to support.

**No templating engine.** The page is small and the escaping is explicit, so a
dependency would buy nothing. Every value that reaches the page goes through
`escape` -- sequence data included, since it arrives from a file.
"""

from __future__ import annotations

from html import escape

from .epitopes import REFERENCE_LABELS, SELF_DEFINITIONS, SelfDefinition
from .number import IMGT_REGIONS, STRUCTURAL_FRAMEWORK_GAPS
from .report import ChainReport, column_note, extraction_note
from .track import render_track_html, track_styles


def _th(label: str, *, num: bool = False) -> str:
    """One table header, carrying its definition as a hover tooltip.

    The definitions come from `report.COLUMN_NOTES` so this page and the app
    explain a column the same way. A label with no entry there renders as a
    plain header rather than failing: a missing tooltip is a smaller problem
    than a crash, and an empty `title` would show as an empty box on hover.
    """
    note = column_note(label)
    classes = " class='num'" if num else ""
    title = f' title="{escape(note)}"' if note else ""
    return f"<th{classes}{title}>{escape(label)}</th>"


_STYLE = """
:root {
  --bg: #fcfcfb;
  --panel: #ffffff;
  --ink: #0b0b0b;
  --muted: #52514e;
  --line: #e2e1dc;
  --accent: #2a78d6;
  --warn-bg: #fdf6e3;
  --warn-line: #eda100;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg: #141413;
    --panel: #1a1a19;
    --ink: #f0efec;
    --muted: #c3c2b7;
    --line: #35352f;
    --accent: #3987e5;
    --warn-bg: #2a2416;
    --warn-line: #c98500;
  }
}
:root[data-theme="dark"] {
  --bg: #141413; --panel: #1a1a19; --ink: #f0efec; --muted: #c3c2b7;
  --line: #35352f; --accent: #3987e5; --warn-bg: #2a2416; --warn-line: #c98500;
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 2rem 1rem 4rem; background: var(--bg); color: var(--ink);
  font: 15px/1.6 ui-sans-serif, system-ui, -apple-system, "Segoe UI", sans-serif;
}
main { max-width: 56rem; margin: 0 auto; }
h1 { font-size: 1.5rem; margin: 0 0 0.25rem; }
h2 {
  font-size: 1.05rem; margin: 2.5rem 0 0.75rem; padding-bottom: 0.35rem;
  border-bottom: 1px solid var(--line); letter-spacing: 0.02em;
}
.sub { color: var(--muted); margin: 0 0 2rem; font-size: 0.9rem; }
table { border-collapse: collapse; width: 100%; margin: 0.5rem 0 0.75rem; font-size: 0.9rem; }
caption { text-align: left; color: var(--muted); font-size: 0.85rem; padding-bottom: 0.4rem; }
th, td { text-align: left; padding: 0.45rem 0.7rem; border-bottom: 1px solid var(--line); }
th {
  font-weight: 600; color: var(--muted); font-size: 0.78rem;
  text-transform: uppercase; letter-spacing: 0.04em; white-space: nowrap;
}
td.num, th.num { text-align: right; font-variant-numeric: tabular-nums; }
/* A column header carrying a definition. The dotted rule is the affordance:
   a native title tooltip is invisible until hovered, so without it a reader
   never learns the explanation is there. cursor:help says the same thing a
   second way, for anyone who does not read the underline as a signal. */
th[title] {
  cursor: help; text-decoration: underline dotted var(--muted);
  text-underline-offset: 3px;
}
tbody tr:last-child td { border-bottom: none; }
.mono, code { font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: 0.86em; }
code {
  background: var(--panel); border: 1px solid var(--line);
  border-radius: 3px; padding: 0 0.25em;
}
.fields { list-style: none; padding: 0; margin: 0.5rem 0; }
.fields li { display: flex; gap: 1rem; padding: 0.3rem 0; border-bottom: 1px solid var(--line); }
.fields li:last-child { border-bottom: none; }
.fields .k { color: var(--muted); min-width: 13rem; }
.lead { font-size: 1.05rem; margin: 0.25rem 0 0.75rem; }
.lead strong { font-variant-numeric: tabular-nums; }
details {
  margin: 0.5rem 0 1rem; padding: 0.6rem 0.9rem; background: var(--panel);
  border: 1px solid var(--line); border-radius: 6px;
}
details > summary {
  cursor: pointer; color: var(--accent); font-size: 0.88rem; list-style: none;
}
details > summary::-webkit-details-marker { display: none; }
details > summary::before { content: "▸ "; }
details[open] > summary::before { content: "▾ "; }
details p { margin: 0.7rem 0 0; color: var(--muted); font-size: 0.9rem; }
.caveat {
  background: var(--warn-bg); border-left: 3px solid var(--warn-line);
  padding: 0.6rem 0.9rem; margin: 0.5rem 0; border-radius: 0 4px 4px 0;
  font-size: 0.9rem;
}
.ok { color: var(--muted); font-size: 0.9rem; }
.unavailable { color: var(--muted); }
footer { margin-top: 3rem; color: var(--muted); font-size: 0.82rem; }
"""


def _percent(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def _details(summary: str, *paragraphs: str) -> str:
    body = "".join(f"<p>{p}</p>" for p in paragraphs)
    return f"<details><summary>{escape(summary)}</summary>{body}</details>"


def _chain(report: ChainReport) -> str:
    domain, germline = report.domain, report.germline
    rows = [
        ("Type", f"{escape(report.chain_name)} ({escape(domain.chain_type)})"),
        (
            "Numbered",
            f"{len(domain.numbering)} residues, IMGT {domain.imgt_span[0]}-{domain.imgt_span[1]}",
        ),
        ("Closest germline V", f"<code>{escape(germline.gene)}</code>"),
    ]
    extracted = extraction_note(domain)
    if extracted:
        rows.insert(2, ("Scored on", escape(extracted)))
    if germline.ambiguous_genes:
        others = ", ".join(
            f"<code>{escape(c.best_allele)}</code>" for c in germline.ambiguous_genes
        )
        rows.append(("Indistinguishable from", others))
    elif germline.nearest_other_gene is not None:
        other = germline.nearest_other_gene
        rows.append(
            (
                "Next closest gene",
                f"<code>{escape(other.best_allele)}</code> "
                f"({other.positions_behind:.0f} position(s) behind)",
            )
        )
    fields = "".join(
        f'<li><span class="k">{escape(k)}</span><span>{v}</span></li>' for k, v in rows
    )
    return f'<h2>Chain</h2><ul class="fields">{fields}</ul>'


def _regions(report: ChainReport) -> str:
    domain = report.domain
    body = []
    for name, (low, high) in IMGT_REGIONS.items():
        residues = domain.regions.get(name, "")
        body.append(
            f"<tr><td>{escape(name)}</td>"
            f'<td class="num">{low}-{high}</td>'
            f'<td class="num">{len(residues)}/{high - low + 1}</td>'
            f'<td class="mono">{escape(residues)}</td></tr>'
        )
    gaps = ", ".join(str(p) for p in sorted(STRUCTURAL_FRAMEWORK_GAPS.get(domain.chain_type, ())))
    return (
        "<h2>Regions</h2>"
        + render_track_html(domain, include_styles=False)
        + "<table><thead><tr>"
        + _th("Region")
        + _th("IMGT slots", num=True)
        + _th("Filled", num=True)
        + _th("Sequence")
        + "</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
        + _details(
            "How to read the slot counts",
            "IMGT reserves a fixed number of slots per region and fills them outward from "
            "both ends, so a shorter loop leaves the <em>middle</em> slots empty while "
            "everything else stays put. That fixed frame is what makes position 23 mean the "
            "same thing in every antibody.",
            f"Framework slots {escape(gaps)} are empty in essentially every "
            f"{escape(domain.chain_type)} chain, along with the final slot 128 in most "
            "chains, so no deletion is reported for them. An <em>unexpected</em> framework "
            "gap is reported as a caveat.",
            "These are the IMGT boundaries. Kabat assigns four more residues to CDR-L2, "
            "which appear here at the start of FR3.",
        )
    )


def _humanness(report: ChainReport) -> str:
    germline, humanness = report.germline, report.humanness
    v, fr, cdr = (germline.scopes[s] for s in ("V", "framework", "cdr"))
    rows = [
        ("Germline identity, V", _percent(v.identity), f"{v.matches}/{v.compared} positions"),
        ("&nbsp;&nbsp;framework", _percent(fr.identity), "FR1+FR2+FR3; the J gene encodes FR4"),
        ("&nbsp;&nbsp;CDR1+CDR2", _percent(cdr.identity), "CDR3 excluded: junctional in origin"),
        (
            "OASis identity",
            _percent(humanness.score),
            f"{humanness.n_found}/{humanness.n_peptides} 9-mers in "
            f"&ge;{humanness.threshold_min_donor_fraction:.0%} of OAS donors",
        ),
    ]
    body = "".join(
        f'<tr><td>{k}</td><td class="num">{escape(value)}</td>'
        f'<td class="unavailable">{basis}</td></tr>'
        for k, value, basis in rows
    )
    return (
        "<h2>Humanness</h2><table><thead><tr>"
        + _th("Measure")
        + _th("Value", num=True)
        + _th("Basis")
        + f"</tr></thead><tbody>{body}</tbody></table>"
        + _details(
            "Why germline identity leads",
            "It separates therapeutic antibodies by origin at least as well as the "
            "repertoire score on the committed panel, and it is far simpler to defend. "
            "OASis is reported beside it because it covers CDR3, which a V gene does not "
            "encode and germline identity therefore cannot reach.",
        )
    )


def _epitopes(report: ChainReport, self_definition: SelfDefinition) -> str:
    burden = report.epitopes
    if burden is None:
        reason = report.epitope_error or "not requested"
        return (
            "<h2>MHC class II</h2>"
            f'<p class="unavailable">Not available: {escape(str(reason))}. '
            "Every other number on this page is computed offline and is unaffected.</p>"
        )

    hits = burden.counts_by_definition()
    cores = burden.core_counts_by_definition()
    rows = []
    for definition in SELF_DEFINITIONS:
        label = REFERENCE_LABELS[definition]
        mark = " *" if definition == self_definition else ""
        absent = cores[f"{definition}_unclassifiable"]
        rows.append(
            f"<tr><td>{escape(label + mark)}</td>"
            f'<td class="num">{cores[f"{definition}_matched"]}</td>'
            f'<td class="num"><strong>{cores[f"{definition}_unmatched"]}</strong></td>'
            f'<td class="num">{cores[f"{definition}_one_substitution"]}</td>'
            f'<td class="num">{cores[f"{definition}_distant"]}</td>'
            f'<td class="num">{absent if absent else "-"}</td></tr>'
        )
    return (
        "<h2>MHC class II</h2>"
        f'<p class="lead"><strong>{cores["distinct_cores"]}</strong> distinct 9-mer cores '
        f"predicted to bind, across <strong>{hits['raw']}</strong> peptide-allele "
        "pairs.</p>"
        f'<p class="unavailable">{len(burden.alleles)} HLA-DRB1 alleles &middot; '
        f"<code>{escape(burden.method)}</code> &middot; {burden.peptide_length}-mers "
        f"&middot; percentile rank &le; {burden.binder_rank:g}</p>"
        "<table><caption>Is each core found in a human reference set?</caption>"
        "<thead><tr>"
        + _th("Reference set")
        + _th("Found", num=True)
        + _th("Not found", num=True)
        + _th("1 away", num=True)
        + _th("2+ away", num=True)
        + _th("No entry", num=True)
        + "</tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table>"
        f'<p class="unavailable">found + not found + no entry = {cores["distinct_cores"]} on '
        f"every row. The * marks the reference the rest of this page follows.</p>"
        + _details(
            "What these columns mean",
            "<strong>Not found</strong> means that exact 9-mer does not occur in the "
            "reference set. That is a lookup result. A core the set contains drops down the "
            "priority list and stays on it; a core the set lacks has still to be shown to "
            "provoke a response.",
            "<strong>1 away</strong> and <strong>2+ away</strong> split the not-found cores "
            "by how far the nearest reference peptide sits, and sum to it. The 2+ column "
            "leaves out one-substitution variants and stays exploratory: the effect of a "
            "substitution turns on which position changed, and the number of changes "
            "settles nothing on its own.",
            "<strong>No entry</strong> means germline holds no reference for a junctional "
            "CDR3 core, so the report abstains there. It records an abstention.",
            f"These rows count cores. One core can be presented by several alleles and sit "
            f"in several overlapping windows, so the {hits['raw']} hits overstate how many "
            "distinct peptides are in play.",
        )
    )


def _caveats(report: ChainReport) -> str:
    warnings = report.warnings()
    if not warnings:
        return '<h2>Caveats</h2><p class="ok">None: a complete, unambiguous V-domain.</p>'
    body = "".join(f'<div class="caveat">{escape(w)}</div>' for w in warnings)
    return f"<h2>Caveats ({len(warnings)})</h2>{body}"


def render_html(
    report: ChainReport,
    *,
    title: str | None = None,
    self_definition: SelfDefinition = "germline",
) -> str:
    """Render one report as a self-contained HTML page."""
    heading = title or "Antibody chain immunogenicity triage"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(heading)}</title>
<style>{_STYLE}\n{track_styles()}</style>
</head><body><main>
<h1>{escape(heading)}</h1>
<p class="sub">Computed from sequence alone. Use it as triage ahead of an assay.</p>
{_chain(report)}
{_regions(report)}
{_humanness(report)}
{_epitopes(report, self_definition)}
{_caveats(report)}
<footer>
Every number here is computed from sequence alone. Clinical anti-drug antibody risk is
driven at least as much by aggregation, dose, route of administration, formulation and
the patient's own HLA type and immune status. A predicted MHC II binder is several steps
from a clinical response: the peptide must be generated by proteolytic processing, be
presented, and meet a T cell that exists in that patient and escaped tolerance. Treat
this as triage that decides where to spend an ex vivo T cell assay or MHC-associated
peptide proteomics run.
</footer>
</main></body></html>"""


def render_html_document(
    reports: list[tuple[str, ChainReport]],
    *,
    self_definition: SelfDefinition = "germline",
) -> str:
    """One page for a batch, each chain under its own heading."""
    if len(reports) == 1:
        name, report = reports[0]
        return render_html(
            report,
            title=None if name == "query" else name,
            self_definition=self_definition,
        )
    sections = "\n".join(
        f"<h1>{escape(name)}</h1>"
        + _chain(report)
        + _regions(report)
        + _humanness(report)
        + _epitopes(report, self_definition)
        + _caveats(report)
        for name, report in reports
    )
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Antibody chain immunogenicity triage &mdash; {len(reports)} chains</title>
<style>{_STYLE}\n{track_styles()}</style>
</head><body><main>
<p class="sub">Computed from sequence alone. Use it as triage ahead of an assay.</p>
{sections}
</main></body></html>"""
