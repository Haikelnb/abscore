"""Scoring the labelled panel, and comparing the two humanness metrics on it.

What this answers
-----------------
Two questions a reader is right to ask of any humanness score:

1. Does it separate antibodies whose engineering origin is known? Human
   therapeutics should score above murine and chimeric ones.
2. Does the more elaborate method earn its complexity? OASis 9-mer donor
   prevalence is a much bigger apparatus than "percent identity to the closest
   germline V gene", so the comparison is worth making explicitly.

The answer to (2) on this panel is that it does not: plain germline V identity
separates human from murine and chimeric chains at least as well as OASis
everywhere the panel can test it. That is a negative result, and it is reported
as one. `tests/test_panel.py` asserts it, so a future change that reverses it
fails the suite instead of quietly invalidating the write-up.

Separation is measured as **AUC** -- the probability that a randomly chosen
positive chain outscores a randomly chosen negative one, ties counting a half.
It needs no threshold, which matters because no threshold here is calibrated.

What this does not establish
----------------------------
Every number comes from one panel, and that panel also supplied the structural
framework-gap whitelist, so the two are not independent. Origin labels are human
curation of engineering history rather than anything derived from sequence, and
nothing in the scoring is fitted to them -- which is what makes the separation
meaningful at all. But it is corroboration on a single source. See
`panels/imgt_mab_db/SOURCE.md`.

Nor is origin a clinical endpoint. A human-origin antibody can still be
immunogenic, and the panel's own exceptions are visible in the figure: a few
human-labelled chains score below the most human-like murine chain, mostly
antibodies that reached their function through heavy somatic hypermutation.
Those points are drawn, not dropped.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from html import escape

from .germline import assign_closest_germline
from .humanness import score_domain
from .number import NumberingError, number_sequence
from .panel import PanelEntry, load_panel

#: The metrics compared, as (attribute, short label, long label).
METRICS: tuple[tuple[str, str, str], ...] = (
    ("v_identity", "V identity", "identity to closest germline V"),
    ("framework_identity", "FR identity", "identity across framework only"),
    ("oasis", "OASis", "9-mers seen in human repertoires"),
)

#: Origin labels on the figure's category axis, ordered by expected humanness.
#: Veterinary labels are excluded here: n=18 across both, too few to plot as
#: distributions. They are checked in `tests/test_panel.py` instead.
FIGURE_ORIGINS: tuple[str, ...] = ("Human", "Humanized", "Chimeric", "Mouse")

#: Chain types, in the order the figure stacks them.
CHAIN_TYPES: tuple[str, ...] = ("H", "K", "L")

CHAIN_LABELS = {"H": "heavy", "K": "kappa", "L": "lambda"}

#: The separation tasks reported. Chimeric chains carry a murine V domain, so
#: they belong with the negatives for a V-identity task.
POSITIVE_ORIGINS = frozenset({"Human"})
NEGATIVE_ORIGINS = frozenset({"Mouse", "Chimeric"})


@dataclass(frozen=True)
class ScoredChain:
    """One panel chain with its labels and all three metrics."""

    record_id: str
    name: str
    origin: str
    chain_type: str
    v_identity: float
    framework_identity: float
    oasis: float
    numbered_slots: int


def score_panel(
    entries: Iterable[PanelEntry] | None = None,
    *,
    progress: Callable[[str], None] | None = None,
) -> list[ScoredChain]:
    """Number and score every panel chain that can be numbered and assigned.

    Chains are skipped only for mechanical reasons, each of which is a property
    of the sequence and not of its score: numbering fails, the record holds
    more than one V domain, or no germline assignment is possible. Nothing is
    skipped for scoring badly.
    """
    rows: list[ScoredChain] = []
    source = list(load_panel() if entries is None else entries)
    for index, entry in enumerate(source, start=1):
        if progress and index % 100 == 0:
            progress(f"scored {index}/{len(source)} panel chains")
        try:
            domains = number_sequence(entry.sequence)
        except NumberingError:
            continue
        if len(domains) != 1:
            continue
        domain = domains[0]
        assignment = assign_closest_germline(domain)
        if assignment.identity is None or assignment.framework_identity is None:
            continue
        rows.append(
            ScoredChain(
                record_id=entry.record_id,
                name=entry.name,
                origin=entry.origin,
                chain_type=domain.chain_type,
                v_identity=assignment.identity,
                framework_identity=assignment.framework_identity,
                oasis=score_domain(domain).score,
                numbered_slots=len(domain.numbering),
            )
        )
    return rows


def auc(positives: Sequence[float], negatives: Sequence[float]) -> float:
    """Probability a random positive outscores a random negative; ties count 0.5.

    Computed directly rather than through a rank formula: the panel is ~1000
    chains, the quadratic cost is irrelevant, and the definition is then visible
    in the code.
    """
    if not positives or not negatives:
        raise ValueError("AUC needs at least one positive and one negative")
    wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in positives for n in negatives)
    return wins / (len(positives) * len(negatives))


@dataclass(frozen=True)
class SeparationRow:
    """One (chain type) x (all metrics) comparison.

    `aucs` is keyed by metric attribute. `positives`/`negatives` are group sizes,
    carried because an AUC computed on 13 negatives deserves to be read next to
    that 13.
    """

    group: str
    label: str
    positives: int
    negatives: int
    aucs: dict[str, float]

    @property
    def best_metric(self) -> str | None:
        if not self.aucs:
            return None
        return max(self.aucs, key=lambda attribute: self.aucs[attribute])


def separation_table(
    rows: Sequence[ScoredChain],
    *,
    positive: frozenset[str] = POSITIVE_ORIGINS,
    negative: frozenset[str] = NEGATIVE_ORIGINS,
    min_group: int = 5,
) -> list[SeparationRow]:
    """Separation per chain type, plus a pooled row.

    A chain type with fewer than `min_group` chains on either side is reported
    with its counts and no AUCs, because an AUC on three negatives is a number
    that invites more confidence than it holds. The panel contains no chimeric
    or murine lambda chains at all, so lambda lands in exactly that case --
    which is itself worth printing rather than omitting.
    """
    groups: list[tuple[str, str, list[ScoredChain]]] = [
        ("all", "all chains", list(rows)),
        *(
            (chain, f"{CHAIN_LABELS[chain]} only", [r for r in rows if r.chain_type == chain])
            for chain in CHAIN_TYPES
        ),
    ]

    table: list[SeparationRow] = []
    for group, label, subset in groups:
        positives = [r for r in subset if r.origin in positive]
        negatives = [r for r in subset if r.origin in negative]
        aucs: dict[str, float] = {}
        if len(positives) >= min_group and len(negatives) >= min_group:
            for attribute, _, _ in METRICS:
                aucs[attribute] = auc(
                    [getattr(r, attribute) for r in positives],
                    [getattr(r, attribute) for r in negatives],
                )
        table.append(
            SeparationRow(
                group=group,
                label=label,
                positives=len(positives),
                negatives=len(negatives),
                aucs=aucs,
            )
        )
    return table


def render_separation_text(table: Sequence[SeparationRow], *, width: int = 88) -> str:
    """The comparison table, laid out for a terminal.

    Five decimal places, which looks excessive until you see the kappa row: the
    margin between the best and second-best metric there is smaller than four
    decimals can show, and rounding it away would turn a tie into a win.
    """
    header = f"{'chain type':<14}{'human':>7}{'murine/':>10}" + "".join(
        f"{short:>13}" for _, short, _ in METRICS
    )
    second = f"{'':<14}{'(n)':>7}{'chimeric':>10}" + "".join(f"{'(AUC)':>13}" for _ in METRICS)

    lines = ["-" * width, header, second, "-" * width]
    for row in table:
        cells = f"{row.label:<14}{row.positives:>7}{row.negatives:>10}"
        if row.aucs:
            best = row.best_metric
            for attribute, _, _ in METRICS:
                # A metric absent from this row prints blank. The table is a
                # fixed grid of columns and a row that measured fewer of them
                # must still line up under the right headings.
                if attribute not in row.aucs:
                    cells += f"{'-':>13}"
                    continue
                value = f"{row.aucs[attribute]:.5f}"
                cells += f"{value + ('*' if attribute == best else ' '):>13}"
        else:
            cells += f"{'too few to measure':>39}"
        lines.append(cells)
    lines.append("-" * width)
    lines.append("* best separator in that row. AUC 1.0 = perfect, 0.5 = no information.")
    return "\n".join(lines)


def render_separation_html(table: Sequence[SeparationRow]) -> str:
    """The same table as HTML, for the figure's page.

    The figure's light-mode aqua sits below 3:1 against the surface, so a table
    view of the same comparison is not optional decoration -- it is what makes
    the figure's third series readable for everyone.
    """
    head = "".join(f"<th>{escape(short)} AUC</th>" for _, short, _ in METRICS)
    body: list[str] = []
    for row in table:
        cells = f"<td>{escape(row.label)}</td><td>{row.positives}</td><td>{row.negatives}</td>"
        if row.aucs:
            best = row.best_metric
            for attribute, _, _ in METRICS:
                if attribute not in row.aucs:
                    cells += "<td>-</td>"
                    continue
                mark = " &#9733;" if attribute == best else ""
                cells += f"<td>{row.aucs[attribute]:.5f}{mark}</td>"
        else:
            cells += f'<td colspan="{len(METRICS)}">too few to measure</td>'
        body.append(f"<tr>{cells}</tr>")

    return (
        "<table><caption>Separation of human-origin chains from murine and chimeric "
        "ones, as AUC: the probability that a randomly chosen human chain outscores a "
        "randomly chosen murine or chimeric one. 1.0 is perfect, 0.5 is no information. "
        "&#9733; marks the best separator in each row. Five decimals because the kappa "
        "margin is smaller than four would show.</caption>"
        f"<thead><tr><th>chain type</th><th>human (n)</th>"
        f"<th>murine/chimeric (n)</th>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table>"
    )


def overlapping_chains(
    rows: Sequence[ScoredChain],
    *,
    attribute: str = "v_identity",
    positive: frozenset[str] = POSITIVE_ORIGINS,
    negative: frozenset[str] = frozenset({"Mouse"}),
) -> list[ScoredChain]:
    """Positives scoring at or below the highest-scoring negative.

    These are the chains an AUC of 0.99 does not show you. They are kept in the
    figure and listed in the report, because the condition that produces them --
    a genuinely human antibody carrying heavy somatic hypermutation -- is a
    limitation of germline-identity scoring rather than noise in the panel.

    The default negative set is murine chains alone, narrower than the AUC
    task's. A chimeric antibody pairs a murine V domain with human constant
    regions, but the panel holds one chimeric kappa chain at 0.915 V identity --
    close to the human median. Using the full negative set would put the ceiling
    there and sweep a quarter of the human chains into this list, which says
    something about that one chimeric chain and nothing about the metric.
    Murine chains top out at 0.817, and that is a boundary worth measuring
    against. `chimeric_ceiling()` reports the other case on its own terms.
    """
    negatives = [getattr(r, attribute) for r in rows if r.origin in negative]
    if not negatives:
        return []
    ceiling = max(negatives)
    return sorted(
        (r for r in rows if r.origin in positive and getattr(r, attribute) <= ceiling),
        key=lambda r: getattr(r, attribute),
    )


def highest_scoring(
    rows: Sequence[ScoredChain],
    *,
    origins: frozenset[str],
    attribute: str = "v_identity",
) -> ScoredChain | None:
    """The chain in `origins` with the highest score on `attribute`.

    Used to name the single most human-like negative, because "the pooled AUC is
    0.994 rather than 1.0" is a claim about specific antibodies and should be
    traceable to them.
    """
    candidates = [r for r in rows if r.origin in origins]
    if not candidates:
        return None
    return max(candidates, key=lambda r: getattr(r, attribute))
