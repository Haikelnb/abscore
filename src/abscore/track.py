"""A position track: the sequence laid out against its IMGT slots.

Why this exists. A table row reading `CDR2 | 56-65 | SAS` invites exactly one
question -- ten slots, three letters, where did the rest go? The answer is that
IMGT reserves a fixed number of slots per region and fills them outward from
both ends, so a short loop leaves a hole in the *middle* rather than shifting
everything along. That is the whole point of the scheme: position 23 is the
first cysteine in every antibody regardless of how long its loops are. A count
of filled slots states that; a picture of the slots shows it, and a reader
stops wondering.

Design notes, since each is a decision rather than a default:

**Two colours.** Framework and CDR, taken from the first two
categorical slots of the validated reference palette (blue and orange; both
modes pass the lightness band, chroma floor, CVD separation, normal-vision floor
and 3:1 contrast against their surface). Seven regions would need seven hues,
which is a rainbow for data that only has one meaningful split. Region identity
is carried by the direct label on every row, never by colour alone.

**One row per region, in sequence order.** A single 128-cell track is wider than
any screen and wraps arbitrarily; wrapping at region boundaries is the natural
break and puts each region's slot range beside it.

**Anchors are marked.** IMGT 23, 41, 89, 104 and 118 are the conserved positions
the numbering is pinned to, so they get a ring. Seeing the first cysteine land
on 23 is the fastest confirmation that a chain numbered correctly.

**Gaps are drawn.** An empty slot is a thin outline. Structural
gaps -- the slots empty in nearly every chain of a type -- are shown differently
from an unexpected gap, because conflating them is the bug this project has hit
more than once: an earlier version reported every missing FR4 slot of a
truncated read as a deletion.
"""

from __future__ import annotations

from html import escape

from .number import (
    ANCHOR_POSITIONS,
    CDR_REGIONS,
    IMGT_REGIONS,
    STRUCTURAL_FRAMEWORK_GAPS,
    Domain,
    imgt_sort_key,
)

#: Slots 1 and 2 of the reference categorical palette, light and dark steps.
#: Validated as a pair in both modes; see the module docstring.
FRAMEWORK_COLOUR = ("#2a78d6", "#3987e5")
CDR_COLOUR = ("#eb6834", "#d95926")

CELL = 19  # px per slot; a residue letter is legible at 11px inside this.
GAP = 2  # px between cells, the surface showing through


def _rows(domain: Domain) -> list[tuple[str, int, int, list[tuple[str, str | None]]]]:
    """Per region: name, slot bounds, and one (kind, residue) pair per slot.

    `kind` is what the cell is, so the renderer never re-derives it:

    - "residue"    a residue sits in this slot
    - "loop"       an empty slot inside a CDR. Normal and expected: CDR length
                   varies, IMGT allocates a fixed maximum, and a shorter loop
                   leaves the middle slots empty.
    - "structural" an empty framework slot that is empty in nearly every chain
                   of this type (IMGT 73 in all, 81/82 in kappa and lambda, 10
                   in heavy and lambda).
    - "outside"    beyond the numbered span -- truncation at an end.
    - "gap"        an empty framework slot that is none of the above. Framework
                   is fixed-length, so this is a structural deletion and the
                   pipeline raises a caveat for it.

    Distinguishing "loop" from "gap" is the point. Drawing a short CDR the same
    way as a framework hole would mark every antibody ever made as defective,
    and conflating a truncated end with a deletion is a bug this project has
    already shipped once.
    """
    structural = STRUCTURAL_FRAMEWORK_GAPS.get(domain.chain_type, frozenset())
    present: dict[int, str] = {}
    for position, residue in sorted(domain.numbering.items(), key=lambda kv: imgt_sort_key(kv[0])):
        # Insertion codes share a slot with their base position; the base
        # residue is what the track shows, and the count carries the rest.
        present.setdefault(position[0], residue)

    low_seen, high_seen = domain.imgt_span

    out = []
    for name, (low, high) in IMGT_REGIONS.items():
        in_cdr = name in CDR_REGIONS
        cells: list[tuple[str, str | None]] = []
        for slot in range(low, high + 1):
            residue = present.get(slot)
            if residue is not None:
                cells.append(("residue", residue))
            elif slot < low_seen or slot > high_seen:
                cells.append(("outside", None))
            elif in_cdr:
                cells.append(("loop", None))
            elif slot in structural:
                cells.append(("structural", None))
            else:
                cells.append(("gap", None))
        out.append((name, low, high, cells))
    return out


def _insertion_note(domain: Domain) -> str:
    """CDR3 insertions share a slot, so say so rather than hiding them."""
    bases: dict[int, int] = {}
    for position in domain.numbering:
        if position[1].strip():
            bases[position[0]] = bases.get(position[0], 0) + 1
    if not bases:
        return ""
    described = ", ".join(f"{slot} (+{count})" for slot, count in sorted(bases.items()))
    return (
        f" Insertion codes share a slot with their base position, so {described} "
        "holds more residues than the one cell shown."
    )


def track_styles() -> str:
    """The track's CSS, so a document holding several tracks can emit it once.

    A batch of fifty chains repeating this block fifty times is fifty copies of
    the same rules. Browsers cope, but it is ~1.5 kB per chain of pure
    duplication in a file someone emails.
    """
    fr_light, fr_dark = FRAMEWORK_COLOUR
    cdr_light, cdr_dark = CDR_COLOUR
    return f""".ab-track {{
  --fr: {fr_light};
  --cdr: {cdr_dark};
  --ink: #0b0b0b;
  --muted: #52514e;
  --edge: #c9c8c2;
  font-family: ui-monospace, "SF Mono", Menlo, monospace;
  margin: 0.25rem 0 0.5rem;
  overflow-x: auto;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .ab-track {{
    --fr: {fr_dark}; --cdr: {cdr_dark};
    --ink: #f0efec; --muted: #c3c2b7; --edge: #4a4a46;
  }}
}}
:root[data-theme="dark"] .ab-track {{
  --fr: {fr_dark}; --cdr: {cdr_dark};
  --ink: #f0efec; --muted: #c3c2b7; --edge: #4a4a46;
}}
.ab-row {{ display: flex; align-items: center; gap: 0.5rem; margin-bottom: {GAP}px; }}
.ab-row-fr {{ --role: var(--fr); }}
.ab-row-cdr {{ --role: var(--cdr); }}
.ab-label {{
  width: 5.5rem; flex: none; color: var(--ink); font-size: 12px;
  display: flex; align-items: center; gap: 0.4rem;
}}
.ab-swatch {{
  width: 10px; height: 10px; flex: none; border-radius: 2px; background: var(--role);
}}
.ab-span {{ width: 4.5rem; flex: none; color: var(--muted); font-size: 11px; text-align: right; }}
.ab-cells {{ display: flex; gap: {GAP}px; flex: none; }}
.ab-count {{ color: var(--muted); font-size: 11px; padding-left: 0.5rem; }}
.ab-cell {{
  width: {CELL}px; height: {CELL}px; flex: none; border-radius: 3px;
  display: flex; align-items: center; justify-content: center;
  font-size: 11px; line-height: 1;
}}
.ab-residue {{ background: var(--role); color: #fff; }}
.ab-gap, .ab-loop, .ab-structural, .ab-outside {{ background: none; }}
.ab-loop {{ border: 1px solid color-mix(in oklab, var(--role) 45%, transparent); }}
.ab-structural {{ border: 1px dotted var(--edge); }}
.ab-gap {{ border: 2px dashed var(--role); }}
.ab-outside {{ border: 1px solid transparent; }}
.ab-anchor {{ outline: 2px solid var(--ink); outline-offset: 1px; }}
.ab-legend {{
  color: var(--muted); font-size: 11px; margin-top: 0.4rem;
  font-family: ui-sans-serif, system-ui, sans-serif;
}}"""


def render_track_html(domain: Domain, *, include_styles: bool = True) -> str:
    """The track as a self-contained HTML fragment with its own scoped styles.

    Colours are declared as custom properties under three scopes -- the local
    root, a prefers-color-scheme block and an explicit data-theme block -- so
    the page's theme toggle wins over the OS setting either way.
    """
    rows = _rows(domain)
    anchors = set(ANCHOR_POSITIONS)

    parts: list[str] = []
    for name, low, high, cells in rows:
        is_cdr = name in CDR_REGIONS
        filled = sum(1 for kind, _ in cells if kind == "residue")
        role = "cdr" if is_cdr else "fr"

        cell_html = []
        for offset, (kind, residue) in enumerate(cells):
            slot = low + offset
            classes = [f"ab-cell ab-{kind}"]
            if slot in anchors:
                classes.append("ab-anchor")
            label = escape(residue) if residue else ""
            title = f"IMGT {slot}"
            if residue:
                title += f": {residue}"
            elif kind == "loop":
                title += ": unused slot, this loop is shorter than its allocation"
            elif kind == "structural":
                title += ": empty in nearly every chain of this type"
            elif kind == "outside":
                title += ": outside the numbered span"
            else:
                title += ": unexpected gap"
            if slot in anchors:
                title += f", conserved anchor ({'/'.join(ANCHOR_POSITIONS[slot])})"
            cell_html.append(
                f'<span class="{" ".join(classes)}" title="{escape(title)}">{label}</span>'
            )

        parts.append(
            f'<div class="ab-row ab-row-{role}">'
            f'<span class="ab-label"><span class="ab-swatch"></span>{escape(name)}</span>'
            f'<span class="ab-span">{low}-{high}</span>'
            f'<span class="ab-cells">{"".join(cell_html)}</span>'
            f'<span class="ab-count">{filled}/{len(cells)}</span>'
            f"</div>"
        )

    note = _insertion_note(domain)
    styles = f"<style>{track_styles()}</style>" if include_styles else ""

    return f"""{styles}
<div class="ab-track">
{"".join(parts)}
<div class="ab-legend">
Filled cell = a residue at that slot. Faint outline in a CDR = an unused slot
because that loop is shorter than its allocation, which is normal. Dotted
outline = a framework slot empty in nearly every {escape(domain.chain_type)}
chain, so no deletion is reported. <b>Bold dashed outline = an unexpected
framework gap</b>, the only one of these that is a finding. Blank = outside the
numbered span. Ringed = a conserved anchor (IMGT
{", ".join(str(a) for a in sorted(anchors))}).{escape(note)}
</div>
</div>
""".strip()
