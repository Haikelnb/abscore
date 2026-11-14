"""The validation figure: every panel chain's score, drawn as a point.

Why a strip plot of individual points rather than a box plot or a bar of means.
The interesting part of this comparison is not where the groups sit, which three
numbers could state; it is **how far they overlap**. A box plot hides the four
human chains that score below the most human-like murine chain, and a bar of
group means hides everything. Drawing roughly a thousand points means the
overlap, the outliers and the group sizes are all visible at once, and an
unusual antibody can be found by looking rather than by querying.

Design decisions, each deliberate:

**Small multiples: metric across, chain type down.** Nine panels sharing one y
axis. Reading down a column compares chain types within one metric; reading
across a row answers the question this figure exists for -- whether OASis
separates better than plain germline identity. Both comparisons are then
adjacent rather than remembered.

**One shared y scale, 0.2 to 1.0.** All three metrics are fractions on the same
footing, so giving each panel its own range would make a 0.03 spread look like a
0.3 one. The floor is 0.2 rather than 0 because no chain scores below 0.245, and
rather than 0.5 because clipping a real point to gain axis resolution is the one
thing a figure about overlap must not do.

**Colour carries the metric, not the group.** Three hues, taken in fixed order
from the first three slots of the validated reference palette -- the three that
clear the all-pairs separation floors in both light and dark, which is the
pairlist small multiples require. Colour is redundant with the column title
directly above each panel, so identity never rests on colour alone; a reader
with any form of colour vision deficiency loses nothing. Origin is carried by
position on the category axis, and the groups are ordered by expected humanness
so the trend reads left to right.

**Medians, drawn as a rule rather than a marker.** A wide flat tick across each
group reads as a summary of that group instead of competing with its points.

**Points are jittered deterministically.** The offset comes from a hash of the
record id, so the same chain lands in the same place on every run and two
renderings of the same panel are byte-identical. A random jitter would make the
figure unreproducible for no gain.
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from html import escape

from .validation import (
    CHAIN_LABELS,
    CHAIN_TYPES,
    FIGURE_ORIGINS,
    METRICS,
    NEGATIVE_ORIGINS,
    ScoredChain,
)

#: Slots 1-3 of the reference categorical palette, light and dark steps. These
#: three validate on the all-pairs pairlist in both modes, which is the one
#: small multiples need; a fourth slot would not.
SERIES_COLOURS: tuple[tuple[str, str], ...] = (
    ("#2a78d6", "#3987e5"),  # blue
    ("#eb6834", "#d95926"),  # orange
    ("#1baf7a", "#199e70"),  # aqua
)

Y_FLOOR = 0.2
Y_CEILING = 1.0
Y_TICKS = (0.2, 0.4, 0.6, 0.8, 1.0)

PANEL_WIDTH = 232
PANEL_HEIGHT = 150
COLUMN_GAP = 14
ROW_GAP = 32
MARGIN_LEFT = 64
MARGIN_TOP = 54
MARGIN_RIGHT = 44
#: Space under the last panel for its category labels, before the notes.
CATEGORY_BLOCK = 26
NOTE_LINE_HEIGHT = 14
NOTE_BASELINE = 12

#: Average advance width of the note font, in px per character. Measured on the
#: shipped stack at 10.5px and rounded up, because a note that overflows the
#: canvas is cropped silently by whatever renders the file -- which has already
#: happened once here.
NOTE_CHAR_WIDTH = 6.1

DOT_RADIUS = 2.4
JITTER_SPAN = 0.62  # fraction of a slot width the jitter is allowed to use


def _jitter(record_id: str, chain_type: str) -> float:
    """A stable offset in [-0.5, 0.5], scaled by JITTER_SPAN at draw time."""
    digest = hashlib.sha256(f"{record_id}/{chain_type}".encode()).digest()
    return int.from_bytes(digest[:4], "big") / 0xFFFFFFFF - 0.5


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def _figure_styles() -> str:
    series = "\n".join(
        f"  --series-{index + 1}: {light};" for index, (light, _) in enumerate(SERIES_COLOURS)
    )
    series_dark = " ".join(
        f"--series-{index + 1}: {dark};" for index, (_, dark) in enumerate(SERIES_COLOURS)
    )
    return f""".ab-fig {{
  --surface: #fcfcfb;
  --ink: #0b0b0b;
  --muted: #52514e;
  --grid: #e3e2dd;
{series}
  font-family: ui-sans-serif, system-ui, -apple-system, sans-serif;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .ab-fig {{
    --surface: #1a1a19; --ink: #f0efec; --muted: #c3c2b7; --grid: #333330;
    {series_dark}
  }}
}}
:root[data-theme="dark"] .ab-fig {{
  --surface: #1a1a19; --ink: #f0efec; --muted: #c3c2b7; --grid: #333330;
  {series_dark}
}}
.ab-fig .frame {{ fill: var(--surface); }}
.ab-fig .grid {{ stroke: var(--grid); stroke-width: 1; }}
.ab-fig .axis {{ stroke: var(--grid); stroke-width: 1; }}
.ab-fig .dot {{ fill-opacity: 0.55; stroke: var(--surface); stroke-width: 0.6; }}
.ab-fig .median {{ stroke: var(--ink); stroke-width: 2; stroke-linecap: round; }}
.ab-fig .col-title {{ fill: var(--ink); font-size: 12.5px; font-weight: 600; }}
.ab-fig .col-sub {{ fill: var(--muted); font-size: 10.5px; }}
.ab-fig .row-title {{ fill: var(--ink); font-size: 12px; font-weight: 600; }}
.ab-fig .tick {{ fill: var(--muted); font-size: 10px; font-variant-numeric: tabular-nums; }}
.ab-fig .cat {{ fill: var(--ink); font-size: 10.5px; }}
.ab-fig .note {{ fill: var(--muted); font-size: 10.5px; }}
.ab-fig .ceiling {{ stroke: var(--grid); stroke-width: 1; stroke-dasharray: 3 3; }}"""


def _panel(
    rows: Sequence[ScoredChain],
    attribute: str,
    series: int,
    x0: float,
    y0: float,
    *,
    show_categories: bool,
) -> list[str]:
    """One panel: four jittered groups, their medians, and the shared grid."""
    parts: list[str] = []
    colour = f"var(--series-{series + 1})"

    def y_of(value: float) -> float:
        fraction = (value - Y_FLOOR) / (Y_CEILING - Y_FLOOR)
        return y0 + PANEL_HEIGHT - fraction * PANEL_HEIGHT

    for tick in Y_TICKS:
        parts.append(
            f'<line class="grid" x1="{x0:.1f}" x2="{x0 + PANEL_WIDTH:.1f}" '
            f'y1="{y_of(tick):.1f}" y2="{y_of(tick):.1f}"/>'
        )

    slot = PANEL_WIDTH / len(FIGURE_ORIGINS)
    for index, origin in enumerate(FIGURE_ORIGINS):
        centre = x0 + slot * (index + 0.5)
        group = [r for r in rows if r.origin == origin]

        for chain in group:
            value = getattr(chain, attribute)
            offset = _jitter(chain.record_id, chain.chain_type) * slot * JITTER_SPAN
            parts.append(
                f'<circle class="dot" fill="{colour}" cx="{centre + offset:.1f}" '
                f'cy="{y_of(value):.1f}" r="{DOT_RADIUS}">'
                f"<title>{escape(chain.name)} ({escape(chain.record_id)}), "
                f"{escape(origin)} {CHAIN_LABELS[chain.chain_type]}: {value:.3f}</title>"
                f"</circle>"
            )

        if group:
            median = _median([getattr(r, attribute) for r in group])
            half = slot * 0.32
            parts.append(
                f'<line class="median" x1="{centre - half:.1f}" x2="{centre + half:.1f}" '
                f'y1="{y_of(median):.1f}" y2="{y_of(median):.1f}">'
                f"<title>{escape(origin)} median {median:.3f} (n={len(group)})</title>"
                f"</line>"
            )

        # Group sizes are deliberately absent. Printed under every panel they
        # were thirty-six numbers competing with the points for attention, and
        # printed once at the foot of the figure they read as column totals
        # when they described only the bottom row. The separation table beside
        # the figure carries the counts instead, where a reader can line them
        # up against the comparison they belong to.
        if show_categories:
            parts.append(
                f'<text class="cat" x="{centre:.1f}" y="{y0 + PANEL_HEIGHT + 16:.1f}" '
                f'text-anchor="middle">{escape(origin)}</text>'
            )

    parts.append(
        f'<line class="axis" x1="{x0:.1f}" x2="{x0 + PANEL_WIDTH:.1f}" '
        f'y1="{y0 + PANEL_HEIGHT:.1f}" y2="{y0 + PANEL_HEIGHT:.1f}"/>'
    )
    return parts


def _wrap(text: str, limit: int) -> list[str]:
    """Greedy word wrap at a character budget.

    SVG has no text layout, so a long string runs off the canvas and the
    renderer crops it without complaint. Measuring text properly needs the font
    metrics; a conservative characters-per-line budget is enough here and keeps
    the figure dependency-free.
    """
    lines: list[str] = []
    current = ""
    for word in text.split():
        candidate = f"{current} {word}".strip()
        if current and len(candidate) > limit:
            lines.append(current)
            current = word
        else:
            current = candidate
    if current:
        lines.append(current)
    return lines


def render_figure_svg(rows: Sequence[ScoredChain], *, include_styles: bool = True) -> str:
    """The nine-panel figure as a standalone SVG element.

    Styles are scoped to `.ab-fig` and declared under three scopes -- local
    root, prefers-color-scheme, and an explicit data-theme block -- so the
    figure follows a page's theme toggle as well as the OS setting, and so
    several figures in one document can share one stylesheet.
    """
    width = MARGIN_LEFT + len(METRICS) * PANEL_WIDTH + (len(METRICS) - 1) * COLUMN_GAP
    width += MARGIN_RIGHT
    negatives = sorted(NEGATIVE_ORIGINS)
    note_limit = int((width - 16) / NOTE_CHAR_WIDTH)
    note_lines = [
        line
        for note in (
            (
                f"Each point is one chain; a flat rule is the group median. Origins run "
                f"left to right by expected humanness. One shared y scale "
                f"({Y_FLOOR:g}-{Y_CEILING:g}); no point is clipped, and points below the "
                f"{' / '.join(negatives).lower()} range are kept in place."
            ),
        )
        for line in _wrap(note, note_limit)
    ]

    # The canvas grows with the notes rather than the notes being trimmed to a
    # fixed canvas.
    height = MARGIN_TOP + len(CHAIN_TYPES) * PANEL_HEIGHT + (len(CHAIN_TYPES) - 1) * ROW_GAP
    height += CATEGORY_BLOCK + len(note_lines) * NOTE_LINE_HEIGHT + NOTE_BASELINE

    parts: list[str] = [f'<rect class="frame" x="0" y="0" width="{width}" height="{height}"/>']

    for column, (_, short, long) in enumerate(METRICS):
        x0 = MARGIN_LEFT + column * (PANEL_WIDTH + COLUMN_GAP)
        parts.append(
            f'<rect x="{x0:.1f}" y="20" width="11" height="11" rx="2" '
            f'fill="var(--series-{column + 1})"/>'
        )
        parts.append(f'<text class="col-title" x="{x0 + 17:.1f}" y="29.5">{escape(short)}</text>')
        parts.append(f'<text class="col-sub" x="{x0:.1f}" y="45">{escape(long)}</text>')

    for row, chain in enumerate(CHAIN_TYPES):
        y0 = MARGIN_TOP + row * (PANEL_HEIGHT + ROW_GAP)
        # Only the origins the figure draws, so the row count matches the
        # points. The panel also holds mixed labels ("Human/Humanized") and
        # veterinary ones, which are scored but not plotted.
        subset = [r for r in rows if r.chain_type == chain and r.origin in FIGURE_ORIGINS]
        last_row = row == len(CHAIN_TYPES) - 1

        parts.append(
            f'<text class="row-title" x="8" y="{y0 + 12:.1f}">{escape(CHAIN_LABELS[chain])}</text>'
        )
        for tick in Y_TICKS:
            fraction = (tick - Y_FLOOR) / (Y_CEILING - Y_FLOOR)
            y = y0 + PANEL_HEIGHT - fraction * PANEL_HEIGHT
            parts.append(
                f'<text class="tick" x="{MARGIN_LEFT - 7}" y="{y + 3.5:.1f}" '
                f'text-anchor="end">{tick:.1f}</text>'
            )

        for column, (attribute, _, _) in enumerate(METRICS):
            x0 = MARGIN_LEFT + column * (PANEL_WIDTH + COLUMN_GAP)
            parts.extend(_panel(subset, attribute, column, x0, y0, show_categories=last_row))

    first_note = height - NOTE_BASELINE - (len(note_lines) - 1) * NOTE_LINE_HEIGHT
    for index, line in enumerate(note_lines):
        y = first_note + index * NOTE_LINE_HEIGHT
        parts.append(f'<text class="note" x="8" y="{y:.0f}">{escape(line)}</text>')

    styles = f"<style>{_figure_styles()}</style>" if include_styles else ""
    body = "\n".join(parts)
    return (
        f'<svg class="ab-fig" xmlns="http://www.w3.org/2000/svg" role="img" '
        f'viewBox="0 0 {width} {height}" width="{width}" height="{height}" '
        f'aria-labelledby="ab-fig-title ab-fig-desc">'
        f'<title id="ab-fig-title">Humanness scores by engineering origin and chain type'
        f"</title>"
        f'<desc id="ab-fig-desc">Nine panels: three humanness metrics across, three chain '
        f"types down. Within each panel, one jittered point per antibody chain, grouped by "
        f"curated engineering origin, with the group median drawn as a flat rule. The "
        f"accompanying separation table gives the same comparison numerically.</desc>"
        f"{styles}{body}</svg>"
    )


def render_figure_document(rows: Sequence[ScoredChain], *, table_html: str = "") -> str:
    """A standalone HTML page holding the figure and, optionally, the table.

    Self-contained: no external stylesheet, font or script, so it opens offline
    and can be attached to an email unchanged.
    """
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>abscore validation</title>
<style>
:root {{ --page: #fcfcfb; --ink: #0b0b0b; --muted: #52514e; --edge: #d9d8d3; }}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --page: #141413; --ink: #f0efec; --muted: #c3c2b7; --edge: #3a3a37;
  }}
}}
:root[data-theme="dark"] {{
  --page: #141413; --ink: #f0efec; --muted: #c3c2b7; --edge: #3a3a37;
}}
body {{
  background: var(--page); color: var(--ink); margin: 0;
  font: 15px/1.6 ui-sans-serif, system-ui, -apple-system, sans-serif;
}}
main {{ max-width: 920px; margin: 0 auto; padding: 2rem 16px 4rem; }}
h1 {{ font-size: 1.4rem; margin: 0 0 0.25rem; }}
p.lede {{ color: var(--muted); margin: 0 0 1.5rem; }}
figure {{ margin: 0 0 2rem; overflow-x: auto; }}
table {{ border-collapse: collapse; font-variant-numeric: tabular-nums; font-size: 13.5px; }}
th, td {{ border-bottom: 1px solid var(--edge); padding: 0.4rem 0.75rem; text-align: right; }}
th:first-child, td:first-child {{ text-align: left; }}
thead th {{ color: var(--muted); font-weight: 600; }}
caption {{ caption-side: bottom; color: var(--muted); font-size: 12px;
           text-align: left; padding-top: 0.6rem; }}
</style></head>
<body><main>
<h1>Humanness scores against curated engineering origin</h1>
<p class="lede">Every chain of the committed therapeutic panel that numbers and
receives a germline assignment. Nothing in the scoring is fitted to these
labels.</p>
<figure>{render_figure_svg(rows)}</figure>
{table_html}
</main></body></html>
"""
