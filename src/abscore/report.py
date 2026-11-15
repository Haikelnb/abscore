"""End to end: one sequence in, one report out.

This is the only module that composes the others, and the only place that
decides what a reader sees first. Three principles it follows:

**No bare numbers.** Every score carries its denominator, its reference set, and
any warning attached to the domain it came from. A percentage with no stated
basis is not reportable.

**The simple baseline is the headline.** Germline identity matches or beats the
repertoire score at separating therapeutic antibodies by origin on the committed
panel (AUC 0.996 vs 0.989 for human against murine and chimeric), so it leads.
The repertoire score is reported beside it because it sees CDR3, which germline
identity structurally cannot.

The epitope section is optional and degrades. It needs a live call to IEDB, and a
free public service should not be a hard dependency in a tool someone else runs:
if it is unavailable the rest of the report is produced with that section marked
unavailable, naming the service.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Literal

from .epitopes import REFERENCE_LABELS, SELF_DEFINITIONS, EpitopeBurden, epitope_burden
from .germline import GermlineAssignment, assign_closest_germline
from .humanness import HumannessScore, score_domain
from .iedb import IEDBError
from .number import PARTIAL_DOMAIN_PREFIX, Domain, number_sequence

SelfDefinition = Literal["germline", "repertoire", "proteome"]

#: One short definition per table column, shared by every renderer that can show
#: it on hover: the HTML report puts it in a `title` attribute, the Streamlit app
#: passes it to `column_config(help=...)`.
#:
#: It lives here rather than in each renderer because the same column must not be
#: explained two different ways in two places. A reader who hovers "Not found" in
#: the browser and hovers it again in the app has to get the same sentence, or the
#: definitions are worse than none.
#:
#: Each entry is one phrase. It states what the number counts, not why it matters:
#: the longer reasoning stays in the prose under each table, which is where a
#: reader who wants it will look. A tooltip that has to be read twice has failed.
COLUMN_NOTES: dict[str, str] = {
    # Regions
    "Region": "IMGT framework or CDR region of this V-domain",
    "IMGT slots": "the fixed position range IMGT reserves for this region",
    "Filled": "residues present out of slots reserved; a shorter loop leaves the middle empty",
    "Sequence": "the residues occupying those slots, in order",
    # Humanness
    "Measure": "which humanness number this row reports",
    "Value": "percent identity, or percent of 9-mers, depending on the measure",
    "Basis": "what it was computed over, and the denominator behind it",
    # Epitope reference-set table
    "Reference set": "the human peptide set each predicted core is looked up in",
    "Found": "exact 9-mer is present in that reference set",
    "Not found": "exact 9-mer is absent from it; a lookup result, not a risk score",
    "1 away": "of the not-found, those whose nearest reference peptide differs at one position",
    "2+ away": "of the not-found, those differing at two or more; exploratory, never a headline",
    "No entry": "the reference cannot classify this core either way, so the report abstains",
    # Batch table
    "Name": "the FASTA record this row was scored from",
    "Chain": "H heavy, K kappa, L lambda",
    "Germline V": "closest human germline V gene",
    "V identity": "percent identity to that gene across the whole V region",
    "FR identity": "the same, framework only (FR1+FR2+FR3); FR4 is J-encoded",
    "OASis": "percent of the chain's 9-mers seen in human repertoires",
    "Hits": "predicted binders, one per peptide-allele pair",
    "Cores": "distinct 9-mer binding cores, however many alleles present each",
    "Caveats": "findings this chain raised; see its full report",
    # Flagged runs, and the regions a predicted peptide spans.
    #
    # Plural on purpose. A 15-mer peptide is long enough to cross two region
    # boundaries, and usually does: IMGT reserves 12 slots for CDR1 and 10 for
    # CDR2, but a short loop fills far fewer, so no window can sit inside the
    # loop alone and every one of them reaches into the flanking frameworks.
    # On trastuzumab's heavy chain 18 of 77 hits span three regions. Three
    # regions in a row is the expected shape of a short CDR, not an anomaly,
    # and this note has to say so or a reader reads it as one.
    "IMGT span": "the stretch of positions, by IMGT number",
    "Positions flagged": "positions inside the span covered by at least one unseen 9-mer",
    "Regions": "every region this stretch covers; a short CDR forces its flanking frameworks in",
    # Strongest hits
    "Rank": "percentile rank from the predictor; lower means more likely presented",
    "Allele": "the HLA-DRB1 allele this prediction is for",
    "Peptide": "the submitted 15-mer window",
    "Core": "the 9 residues predicted to sit in the groove",
    "Distance": "residues to the nearest reference peptide: 0, 1, or 2+",
}


def column_note(label: str) -> str | None:
    """The short definition for a column label, if one is defined.

    Falls back to a `(reference)`-suffixed label, so "Unmatched (germline)"
    finds the "Unmatched" entry rather than going unexplained.
    """
    if label in COLUMN_NOTES:
        return COLUMN_NOTES[label]
    base = label.split(" (")[0].strip()
    return COLUMN_NOTES.get(base)


COLUMN_NOTES["Unmatched"] = COLUMN_NOTES["Not found"]
COLUMN_NOTES["Matched"] = COLUMN_NOTES["Found"]


@dataclass
class ChainReport:
    """Everything computed for one antibody chain."""

    sequence: str
    domain: Domain
    germline: GermlineAssignment
    humanness: HumannessScore
    epitopes: EpitopeBurden | None = None
    epitope_error: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def chain_type(self) -> str:
        return self.domain.chain_type

    @property
    def chain_name(self) -> str:
        return {"H": "heavy", "K": "light (kappa)", "L": "light (lambda)"}.get(
            self.chain_type, self.chain_type
        )

    def warnings(self) -> tuple[str, ...]:
        """Every caveat attached to any number in this report, as findings.

        Two kinds of repetition are collapsed here, because a reader counting
        caveats should be counting findings.

        **Several signals, one cause.** `completeness_warnings` runs four
        independent checks for a partial domain -- missing anchors, empty
        regions, where the numbering starts, where it ends -- because no single
        one is sufficient. A badly truncated chain trips all four, which is four
        warnings about one fact. They are merged into one caveat that lists
        every signal that fired.

        **Several consequences, one cause.** Reduced germline coverage, a
        humanness score over fewer peptides and an epitope count over fewer
        windows all follow from the domain being partial, so they are attached
        to that caveat instead of standing beside it.

        An earlier version attached the consequences to *each* partial-domain
        warning, which printed the same paragraph four times on exactly the
        chains that needed reading most carefully.

        **A consequence without its cause is dropped.** The two restatements
        below say only "this number was computed over a reduced chain". If
        nothing reported a reduced chain, that is not a hedge, it is a false
        claim about the input, and it appears in the one block a reader uses to
        decide how much to trust the numbers. This is a backstop: the modules
        that raise them now test the same prefix, so an orphan means a new
        consequence was added without a cause to hang it on.

        Each module still raises its own warning, because each is consumed on
        its own elsewhere. The merging happens here, where they meet.
        """
        #: Pure restatements of "the domain is partial". Meaningless alone.
        restatements = (
            "humanness computed on a partial domain",
            "epitope burden computed on a partial domain",
        )
        #: A measurement that happens to follow from partiality. It carries a
        #: coverage figure computed from the alignment, so it stands on its own
        #: and is kept even when no partial-domain cause was reported.
        measured = ("germline identity computed at",)

        partial: list[str] = []
        others: list[str] = []
        derived: list[str] = []
        for message in (
            *self.domain.warnings,
            *self.germline.warnings(),
            *self.humanness.warnings(),
            *(self.epitopes.warnings() if self.epitopes else ()),
        ):
            if message in partial or message in others or message in derived:
                continue
            if message.startswith(PARTIAL_DOMAIN_PREFIX):
                partial.append(message)
            elif message.startswith(restatements + measured):
                # Classified without consulting `partial`, so the result does
                # not depend on which module spoke first.
                derived.append(message)
            else:
                others.append(message)

        if not partial:
            kept = [m for m in derived if not m.startswith(restatements)]
            return tuple(others + kept)

        # One caveat for the partial domain, listing every signal that fired.
        signals = [m[len(PARTIAL_DOMAIN_PREFIX) :].rstrip(".") for m in partial]
        merged = PARTIAL_DOMAIN_PREFIX + "; ".join(signals) + "."
        if derived:
            merged += " Consequences: " + " ".join(derived)
        return tuple([merged, *others])

    def provenance(self) -> dict[str, object]:
        from .germline import germline_provenance

        record: dict[str, object] = {
            "germline_identity": germline_provenance(),
            "humanness": self.humanness.provenance(),
        }
        if self.epitopes:
            record["epitopes"] = self.epitopes.provenance()
        elif self.epitope_error:
            record["epitopes"] = {"status": "unavailable", "reason": self.epitope_error}
        return record


def analyse(
    sequence: str,
    *,
    with_epitopes: bool = True,
    self_definition: SelfDefinition = "germline",
    cache_dir: str | None = None,
    offline: bool = False,
    progress: Callable[[str], None] | None = None,
    **epitope_kwargs,
) -> ChainReport:
    """Run the whole pipeline on one chain.

    Raises NumberingError if the sequence is not a single antibody V-domain.
    Epitope prediction failures are captured and reported, so the rest of the
    report survives an unavailable service.
    """

    def step(message: str) -> None:
        if progress:
            progress(message)

    step("numbering")
    domains = number_sequence(sequence)
    if len(domains) != 1:
        kinds = ", ".join(d.chain_type for d in domains)
        from .number import NumberingError

        raise NumberingError(
            f"found {len(domains)} V-domains ({kinds}). Analyse one chain at a time."
        )
    domain = domains[0]

    step(f"numbered: {domain.chain_type} chain, IMGT {domain.imgt_span[0]}-{domain.imgt_span[1]}")

    step("germline assignment")
    germline = assign_closest_germline(domain)
    step(f"germline: {germline.gene}, V {germline.identity:.1%}")

    step("humanness")
    humanness = score_domain(domain)
    step(f"humanness: OASis {humanness.score:.1%}")

    burden: EpitopeBurden | None = None
    error: str | None = None
    if with_epitopes:
        step("epitopes")
        kwargs = dict(epitope_kwargs)
        if cache_dir is not None:
            kwargs["cache_dir"] = cache_dir
        if progress is not None:
            kwargs["progress"] = progress
        try:
            burden = epitope_burden(domain, offline=offline, **kwargs)
            # Cores, not hits, for the unmatched figure. `count('germline')`
            # returns peptide-allele pairs, and printing that beside a hit
            # count invites reading both as the same unit. The table this line
            # precedes counts cores, so this one does too.
            cores = burden.core_counts_by_definition()
            step(
                f"epitopes: {burden.count()} hits, {cores['distinct_cores']} cores, "
                f"{cores['germline_unmatched']} not in germline"
            )
        except IEDBError as exc:
            error = str(exc)
            step("epitopes: unavailable")

    return ChainReport(
        sequence=domain.seq,
        domain=domain,
        germline=germline,
        humanness=humanness,
        epitopes=burden,
        epitope_error=error,
    )


def escape_markdown(text: str) -> str:
    """Escape the characters that make an IMGT allele name render as italics.

    Every IMGT allele name contains an asterisk -- `IGHV3-30-3*03` -- so a
    message naming two of them puts a matched pair of asterisks in the string,
    and any markdown renderer reads the text between them as emphasis. The
    observed result: "IGHV3-30-3*03 cannot be separated from IGHV3-30*17"
    displayed as "IGHV3-30-303 cannot be separated from IGHV3-30" in italics,
    with both allele numbers silently deleted. A caveat about which gene was
    called is worse than useless if it corrupts the gene names.

    Warning strings stay plain text, because they also go to stderr and to
    JSON, where backslashes would be noise. Escaping therefore happens at each
    display boundary that renders markdown -- here and in the Streamlit app.
    """
    for character in ("\\", "*", "_", "`"):
        text = text.replace(character, "\\" + character)
    return text


def _percent(value: float | None) -> str:
    return "n/a" if value is None else f"{value:.1%}"


def extraction_note(domain) -> str | None:
    """How much of the input was scored, when that is not all of it.

    Every score is computed over the extracted V-domain. For a bare variable
    domain the input and the scored sequence are the same and this returns
    None. For a full-length chain, a tagged construct or anything with flanks
    they differ, and a report that does not say so gives no way to tell which
    residues produced the numbers.

    It also exposes a boundary effect that is otherwise silent: the numbering
    HMM decides where the domain ends, and that decision can move by a residue
    when the flanking sequence changes, which shifts the score slightly. The
    extracted length is what makes two such runs distinguishable.
    """
    outside = domain.residues_outside_domain
    if outside <= 0:
        return None
    return (
        f"{len(domain.seq)} of {domain.input_length} input residues "
        f"({outside} outside the domain were not scored)"
    )


def render_markdown(
    report: ChainReport,
    *,
    title: str | None = None,
    self_definition: SelfDefinition = "germline",
) -> str:
    """Render a report as markdown.

    Ordered so that stopping after any section still leaves something honest:
    what was found, then the scores with their bases, then where the problems
    are, then every caveat, then provenance.
    """
    domain = report.domain
    germline = report.germline
    lines: list[str] = []

    lines.append(f"# {title or 'Antibody chain immunogenicity triage'}")
    lines.append("")

    # --- what this is -------------------------------------------------------
    lines.append("## Chain")
    lines.append("")
    lines.append(f"- **Type:** {report.chain_name} chain ({domain.chain_type})")
    span = f"{domain.imgt_span[0]}-{domain.imgt_span[1]}"
    lines.append(f"- **Residues numbered:** {len(domain.numbering)}, IMGT {span}")
    extracted = extraction_note(domain)
    if extracted:
        lines.append(f"- **Scored on:** {extracted}")
    lines.append(f"- **Closest human germline V:** `{germline.gene}`")
    if germline.ambiguous_genes:
        others = ", ".join(f"`{c.best_allele}`" for c in germline.ambiguous_genes)
        lines.append(f"  - cannot be separated from: {others}")
    lines.append("")

    lines.append("| Region | IMGT slots | Filled | Sequence |")
    lines.append("|---|---|---|---|")
    from .number import IMGT_REGIONS, STRUCTURAL_FRAMEWORK_GAPS

    for name, (low, high) in IMGT_REGIONS.items():
        residues = domain.regions.get(name, "")
        lines.append(f"| {name} | {low}-{high} | {len(residues)}/{high - low + 1} | `{residues}` |")
    lines.append("")
    lines.append(
        "**How to read Filled.** IMGT reserves a fixed number of slots per "
        "region and fills them outward from both ends, so a shorter loop leaves "
        "the *middle* slots empty while everything else stays put -- that "
        "fixed frame is what makes position 23 mean the same thing in every "
        f"antibody. This chain's CDR2 fills "
        f"{len(domain.regions.get('CDR2', ''))} of "
        f"{IMGT_REGIONS['CDR2'][1] - IMGT_REGIONS['CDR2'][0] + 1} slots; a light "
        "chain CDR2 is typically 3 residues, at positions 56, 57 and 65."
    )
    lines.append("")
    lines.append(
        f"Framework regions also have slots that are empty in essentially every "
        f"chain of a given type: IMGT 73 in all of them, 81 and 82 in kappa and "
        f"lambda, 10 in heavy and lambda, and the final slot 128 in most chains. "
        f"This chain is missing "
        f"{', '.join(str(p) for p in sorted(STRUCTURAL_FRAMEWORK_GAPS.get(domain.chain_type, ())))}"
        " for that reason, which is why no deletion is reported. An *unexpected* "
        "framework gap is reported as a caveat instead."
    )
    lines.append("")
    lines.append(
        "One boundary to be aware of: these are the IMGT definitions. Kabat "
        "assigns four more residues to CDR-L2 which appear here at the start of "
        "FR3, so a CDR2 that looks short against Kabat is not missing anything."
    )
    lines.append("")

    # --- scores -------------------------------------------------------------
    lines.append("## Humanness")
    lines.append("")
    lines.append("| Measure | Value | Basis |")
    lines.append("|---|---|---|")
    v = germline.scopes["V"]
    fr = germline.scopes["framework"]
    cdr = germline.scopes["cdr"]
    lines.append(
        f"| **Germline identity, V** | {_percent(v.identity)} | "
        f"{v.matches}/{v.compared} positions vs {germline.gene} |"
    )
    lines.append(
        f"| **Germline identity, framework** | {_percent(fr.identity)} | "
        f"FR1+FR2+FR3 only; FR4 is J-encoded |"
    )
    lines.append(
        f"| Germline identity, CDR1+CDR2 | {_percent(cdr.identity)} | "
        f"CDR3 excluded: it is junctional in origin |"
    )
    lines.append(
        f"| **OASis identity** | {_percent(report.humanness.score)} | "
        f"{report.humanness.n_found}/{report.humanness.n_peptides} 9-mers found in "
        f"≥{report.humanness.threshold_min_donor_fraction:.0%} of OAS donors |"
    )
    if not v.is_full_coverage and v.coverage is not None:
        lines.append(
            f"| Coverage | {_percent(v.coverage)} | {v.compared} of {v.covered} "
            "germline positions observed; identity not comparable with full-coverage values |"
        )
    lines.append("")
    lines.append(
        "Germline identity leads because it separates therapeutic antibodies by "
        "origin at least as well as the repertoire score on the committed panel. "
        "OASis is reported beside it because it covers CDR3, which a V gene does "
        "not encode and germline identity therefore cannot see."
    )
    lines.append("")

    # --- where it stops looking human ---------------------------------------
    runs = report.humanness.flagged_runs()
    if runs:
        lines.append("### Stretches not found in the OAS reference")
        lines.append("")
        lines.append("| IMGT span | Positions flagged | Regions |")
        lines.append("|---|---|---|")
        for start, end, length, regions in sorted(runs, key=lambda r: -r[2]):
            lines.append(
                f"| {start[0]}{start[1].strip()}-{end[0]}{end[1].strip()} | "
                f"{length} | {', '.join(regions)} |"
            )
        lines.append("")
        lines.append(
            f"A run is where the sequence stops matching the reference set -- "
            f"9-mers carried by at least "
            f"{report.humanness.threshold_min_donor_fraction:.0%} of OAS donors. "
            "Absent from that set is not the same as non-human: a peptide one "
            "donor carries and most do not is absent from it and still human. "
            "A flagged position is covered by at least one unmatched 9-mer, "
            "which leaves open whether every 9-mer covering it is unmatched. A "
            "run is an observation and carries no cause: "
            "a grafted CDR, heavy somatic mutation and a sequencing artefact "
            "all produce one."
        )
        lines.append("")

    # --- epitopes -----------------------------------------------------------
    lines.append("## MHC class II epitope burden")
    lines.append("")
    if report.epitopes is None:
        lines.append(
            "**Not available.** This section needs a live prediction from the "
            "IEDB Analysis Resource, which was not reachable."
        )
        if report.epitope_error:
            lines.append("")
            lines.append(f"> {report.epitope_error}")
        lines.append("")
        lines.append("Every other number in this report is computed offline and is unaffected.")
    else:
        burden = report.epitopes
        counts = burden.counts_by_definition()
        lines.append(
            f"Predicted against {len(burden.alleles)} HLA-DRB1 alleles with "
            f"`{burden.method}`, {burden.peptide_length}-mers, percentile rank "
            f"≤ {burden.binder_rank:g}."
        )
        if burden.method.endswith("_el"):
            lines.append("")
            lines.append(
                "The rank is over an eluted-ligand score, which covers "
                "presentation rather than binding affinity alone; "
                "`netmhciipan_ba` is the affinity output."
            )
        lines.append("")
        lines.append("| Count | Predicted hits | Distinct 9-mer cores |")
        lines.append("|---|---|---|")
        lines.append(f"| All predicted binders | {counts['raw']} | {counts['distinct_cores']} |")
        for definition in SELF_DEFINITIONS:
            lines.append(
                f"| **Unmatched, *{definition}*** | {counts[f'{definition}_unmatched']} | "
                f"{len(burden.unique_cores(definition))} |"
            )
            lines.append(
                f"| &nbsp;&nbsp;matched in reference, *{definition}* | "
                f"{counts[f'{definition}_matched']} | - |"
            )
            lines.append(
                f"| &nbsp;&nbsp;unmatched, nearest 1 substitution away, *{definition}* | "
                f"{counts[f'{definition}_one_substitution']} | - |"
            )
            lines.append(
                f"| &nbsp;&nbsp;unmatched, nearest 2+ away *(exploratory)*, *{definition}* | "
                f"{counts[f'{definition}_distant']} | - |"
            )
            if counts[f"{definition}_unclassifiable"]:
                lines.append(
                    f"| &nbsp;&nbsp;unclassifiable, *{definition}* | "
                    f"{counts[f'{definition}_unclassifiable']} | - |"
                )
        lines.append("")
        lines.append(
            "**The units differ.** A predicted hit is one peptide-allele pair; "
            "a distinct core is one 9-mer however many alleles present it. "
            "Neither is a count of experimentally confirmed epitopes."
        )
        lines.append("")
        lines.append(
            "**What unmatched means.** The exact 9-mer does not occur in the "
            "reference set. Tolerance deletes or suppresses most T cells "
            "reactive to self peptides, so a core the reference contains drops "
            "down the priority list. It stays on the list: central and "
            "peripheral tolerance are both incomplete, self-reactive T cells "
            "persist in healthy people, and a fully human sequence still "
            "carries many predicted hits."
        )
        lines.append("")
        lines.append(
            "**How to read the distance.** The two sub-rows split the "
            "unmatched count by how far the nearest reference peptide sits, and "
            "they sum back to it. The 2+ row leaves out one-substitution "
            "variants. It appears here because that assumption is common in the "
            "literature, and it is marked *exploratory* because the effect of a "
            "substitution turns on which position changed and whether a T-cell "
            "receptor contacts it. The number of changes alone settles nothing. "
            "Distances beyond one go uncomputed, so 2+ always means at least "
            "two."
        )
        lines.append("")
        lines.append(
            "**Totals reconcile.** Per reference set, matched + unmatched + "
            "unclassifiable equals the predicted hits."
        )
        lines.append("")
        lines.append(
            "Peptides in the junctional window are *unclassifiable* by the "
            "germline definition, which has no reference for them: a germline D "
            "gene has no single protein sequence. Unclassifiable records an "
            "abstention. A low germline count therefore leaves that stretch of "
            "CDR3 open, and CDR3 is where a non-germline sequence is most "
            "likely to sit. The repertoire and proteome definitions do cover "
            "these peptides."
        )
        lines.append("")
        one_sub = counts["germline_one_substitution"]
        unmatched_total = counts["germline_unmatched"]
        if one_sub and unmatched_total:
            share = one_sub / unmatched_total
            if share >= 0.5:
                lines.append(
                    f"> **{share:.0%} of the unmatched germline count ({one_sub} of "
                    f"{unmatched_total}) is one substitution from a germline "
                    "peptide.** Consistent with somatic mutation of an otherwise "
                    "germline sequence, though a mismatch count alone cannot "
                    "establish that origin, and proximity alone does not make "
                    "these cores harmless."
                )
                lines.append("")

        unmatched = burden.unmatched(self_definition)
        if unmatched:
            lines.append(f"### Strongest unmatched hits (*{self_definition}* reference)")
            lines.append("")
            lines.append("| Rank | Allele | Peptide | Core | Regions | Distance |")
            lines.append("|---|---|---|---|---|---|")
            for epitope in unmatched[:10]:
                lines.append(
                    f"| {epitope.rank:g} | {epitope.allele} | `{epitope.peptide}` | "
                    f"`{epitope.core}` | {', '.join(epitope.regions)} | "
                    f"{epitope.distance_label(self_definition)} |"
                )
            lines.append("")
            by_region = burden.counts_by_region(self_definition)
            lines.append(
                "Unmatched hits by region: "
                + ", ".join(f"{region} {count}" for region, count in sorted(by_region.items()))
            )
            lines.append("")

    # --- caveats ------------------------------------------------------------
    warnings = report.warnings()
    lines.append("## Caveats on these numbers")
    lines.append("")
    if warnings:
        for message in warnings:
            lines.append(f"- {escape_markdown(message)}")
    else:
        lines.append("- No sequence-level warnings: a complete, unambiguous V-domain.")
    lines.append("")
    lines.append(
        "**What none of this establishes.** Every number here is computed from "
        "sequence alone. Clinical anti-drug antibody risk is driven at least as "
        "much by aggregation, dose, route of administration, formulation and the "
        "patient's own HLA type and immune status. A predicted MHC II binder is "
        "several steps from a clinical response: the peptide must actually be "
        "generated by proteolytic processing, actually be presented, and meet a "
        "T cell that exists in that patient and escaped tolerance. Treat this as "
        "triage that decides where to spend an ex vivo T cell assay or "
        "MHC-associated peptide proteomics run."
    )
    lines.append("")

    return "\n".join(lines)


def render_summary(report: ChainReport) -> str:
    """One line per score, for a CLI or a log."""
    germline = report.germline
    parts = [
        f"chain={report.chain_type}",
        f"germline={germline.gene}",
        f"V={_percent(germline.identity)}",
        f"FR={_percent(germline.framework_identity)}",
        f"OASis={_percent(report.humanness.score)}",
    ]
    if report.epitopes:
        # `hits` counts peptide-allele pairs; the other two count distinct
        # 9-mer cores. Mixing the units here once printed `cores=18
        # unmatched_germline=18` for a chain with 3 unmatched cores, which
        # reads as "every core is unmatched". The key name is unchanged so
        # anything parsing this line keeps working.
        counts = report.epitopes.counts_by_definition()
        cores = report.epitopes.core_counts_by_definition()
        parts.append(f"hits={counts['raw']}")
        parts.append(f"cores={cores['distinct_cores']}")
        parts.append(f"unmatched_germline={cores['germline_unmatched']}")
    else:
        parts.append("binders=unavailable")
    if report.warnings():
        parts.append(f"warnings={len(report.warnings())}")
    return " ".join(parts)


# --- terminal rendering ------------------------------------------------------
#
# Markdown pipe tables are unreadable in a terminal: the pipes do not align, the
# bold markers are literal asterisks, and a wide row wraps into nonsense. So the
# full report is rendered twice -- as markdown for a file, where the tables
# render, and as aligned plain text for the terminal. Same numbers, same order,
# neither one a degraded copy of the other.

TEXT_WIDTH = 88
_GLYPHS = {"residue": None, "loop": "·", "structural": "·", "outside": " ", "gap": "!"}


def _rule(title: str, width: int = TEXT_WIDTH) -> list[str]:
    return ["", title.upper(), "─" * min(width, max(len(title), 20))]


def _field(label: str, value: str, note: str = "", pad: int = 22) -> str:
    line = f"  {label:<{pad}} {value}"
    return f"{line}   {note}" if note else line


def _wrap(text: str, indent: str = "  ", width: int = TEXT_WIDTH) -> list[str]:
    import textwrap

    return textwrap.wrap(text, width=width, initial_indent=indent, subsequent_indent=indent)


def render_text(
    report: ChainReport,
    *,
    title: str | None = None,
    self_definition: SelfDefinition = "germline",
) -> str:
    """Render a report as aligned plain text for a terminal.

    The region view is the ASCII form of the app's position track: one row per
    region, one character per IMGT slot, so a short CDR loop and a framework
    deletion look different here too. `!` is the only mark that is a finding.
    """
    from .track import _rows

    domain = report.domain
    germline = report.germline
    out: list[str] = []

    if title:
        out += [title, "=" * min(len(title), TEXT_WIDTH)]

    out += _rule("Chain")
    out.append(_field("Type", f"{report.chain_name} ({domain.chain_type})"))
    span = f"{domain.imgt_span[0]}-{domain.imgt_span[1]}"
    out.append(_field("Numbered", f"{len(domain.numbering)} residues, IMGT {span}"))
    extracted = extraction_note(domain)
    if extracted:
        out.append(_field("Scored on", extracted))
    out.append(_field("Closest germline V", germline.gene))
    if germline.ambiguous_genes:
        others = ", ".join(c.best_allele for c in germline.ambiguous_genes)
        out.append(_field("Indistinguishable from", others))
    elif germline.nearest_other_gene is not None:
        other = germline.nearest_other_gene
        out.append(
            _field(
                "Next closest gene",
                f"{other.best_allele} ({other.positions_behind:.0f} position(s) behind)",
            )
        )

    out += _rule("Regions")
    for name, low, high, cells in _rows(domain):
        track = "".join(r if r else _GLYPHS[k] for k, r in cells)
        filled = sum(1 for k, _ in cells if k == "residue")
        out.append(f"  {name:<5} {f'{low}-{high}':>8}  {track:<40} {filled}/{len(cells)}")
    out += _wrap(
        "One character per IMGT slot. A shorter loop leaves the middle slots "
        "empty, shown as ·, which is normal; '!' is an unexpected framework "
        "gap and the only mark here that is a finding."
    )

    out += _rule("Humanness")
    v, fr, cdr = (germline.scopes[s] for s in ("V", "framework", "cdr"))
    out.append(
        _field("Germline identity, V", _percent(v.identity), f"{v.matches}/{v.compared} positions")
    )
    out.append(_field("  framework", _percent(fr.identity), "FR1+FR2+FR3; FR4 is J-encoded"))
    out.append(_field("  CDR1+CDR2", _percent(cdr.identity), "CDR3 is junctional, excluded"))
    humanness = report.humanness
    out.append(
        _field(
            "OASis identity",
            _percent(humanness.score),
            f"{humanness.n_found}/{humanness.n_peptides} 9-mers in "
            f"≥{humanness.threshold_min_donor_fraction:.0%} of donors",
        )
    )

    out += _rule("MHC class II")
    if report.epitopes is None:
        reason = report.epitope_error or "not requested"
        out += _wrap(f"Not available: {reason}")
        out += _wrap("Every other number above is computed offline and is unaffected.")
    else:
        burden = report.epitopes
        hits = burden.counts_by_definition()
        cores = burden.core_counts_by_definition()
        out += _wrap(
            f"{cores['distinct_cores']} distinct 9-mer cores predicted to bind, "
            f"across {hits['raw']} peptide-allele pairs."
        )
        out += _wrap(
            f"{len(burden.alleles)} HLA-DRB1 alleles \u00b7 {burden.method} \u00b7 "
            f"{burden.peptide_length}-mers \u00b7 percentile rank \u2264 "
            f"{burden.binder_rank:g}"
        )
        out.append("")
        out += _wrap("Is each core found in a human reference set?")
        out.append("")
        out.append(
            f"  {'reference set':<20}{'found':>6}{'not found':>12}"
            f"{'1 away':>10}{'2+ away':>10}{'no entry':>11}"
        )
        for definition in SELF_DEFINITIONS:
            label = REFERENCE_LABELS[definition]
            if definition == self_definition:
                label += " *"
            absent = cores[f"{definition}_unclassifiable"]
            out.append(
                f"  {label:<20}"
                f"{cores[f'{definition}_matched']:>6}"
                f"{cores[f'{definition}_unmatched']:>12}"
                f"{cores[f'{definition}_one_substitution']:>10}"
                f"{cores[f'{definition}_distant']:>10}"
                f"{(absent if absent else '-'):>11}"
            )
        out.append("")
        out += _wrap(
            "'Not found' means that exact 9-mer does not occur in the reference "
            "set. That is a lookup result. A core the set contains drops down "
            "the priority list and stays on it; a core the set lacks has still "
            "to be shown to provoke a response."
        )
        out += _wrap(
            "'1 away' and '2+ away' split the not-found cores by how far the "
            "nearest reference peptide sits, and sum to it. '2+ away' leaves "
            "out one-substitution variants and stays exploratory: the effect of "
            "a substitution turns on which position changed. 'no entry' means "
            "germline holds no reference for a junctional CDR3 core, so the "
            "report abstains there."
        )
        out += _wrap(f"found + not found + no entry = {cores['distinct_cores']} on every row.")
        out += _wrap(
            f"These rows count cores. One core can be presented by several "
            f"alleles and sit in several overlapping windows, so the "
            f"{hits['raw']} hits overstate how many distinct peptides are in "
            f"play."
        )
        if self_definition:
            out += _wrap("* the reference the rest of this report follows.")

    warnings = report.warnings()
    out += _rule(f"Caveats ({len(warnings)})" if warnings else "Caveats")
    if not warnings:
        out.append("  None: a complete, unambiguous V-domain.")
    for message in warnings:
        wrapped = _wrap(message, indent="    ")
        wrapped[0] = "  - " + wrapped[0].lstrip()
        out += wrapped

    out.append("")
    return "\n".join(out)
