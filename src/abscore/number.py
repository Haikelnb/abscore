"""IMGT numbering, region assignment and germline V/J calls for antibody chains.

Everything downstream of this module works in IMGT coordinates. A raw sequence index
fails as a coordinate for antibodies: CDR loops vary in length, so index 50 in
one chain is not the structurally equivalent position in another. IMGT numbering
assigns fixed position slots and inserts gaps at defined positions, which makes
positions comparable across chains, loci and species.

Numbering is done by ANARCI (profile HMMs over IMGT germline alignments), which also
reports the chain type, the best-matching species, and the closest germline V and J
genes with their identities.

References
----------
- IMGT unique numbering: Lefranc M-P et al., Dev Comp Immunol 2003;27:55-77.
- ANARCI: Dunbar J, Deane CM. Bioinformatics 2016;32:298-300.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Literal

from anarci import run_anarci

# IMGT unique numbering, V-DOMAIN. Inclusive position ranges.
# Lefranc M-P et al. 2003; see also https://www.imgt.org/IMGTScientificChart/
IMGT_REGIONS: dict[str, tuple[int, int]] = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

FRAMEWORK_REGIONS = ("FR1", "FR2", "FR3", "FR4")
CDR_REGIONS = ("CDR1", "CDR2", "CDR3")

# Conserved anchor positions of the IMGT V-DOMAIN. These four residues are
# near-invariant across loci and species, which makes them a biological test for
# whether a complete domain is present -- unlike counting residues, which cannot
# work: IMGT inserts gaps wherever a loop is short, so a complete V-domain
# legitimately leaves 8-22 of the 128 slots empty.
#
#   23  1st-CYS            cysteine
#   41  CONSERVED-TRP      tryptophan
#   89  hydrophobic        (IMGT lists this as "hydrophobic amino acid" rather
#                          residue, so it is checked against a set)
#  104  2nd-CYS            cysteine, immediately before CDR3
#  118  J-PHE or J-TRP     tryptophan in heavy, phenylalanine in light
#
# Verified against https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
# (retrieved 2026-11-12) and against 9 marketed therapeutic chains, see
# tests/fixtures.py.
_HYDROPHOBIC = ("A", "V", "L", "I", "P", "F", "M", "W", "C", "G", "Y")
ANCHOR_POSITIONS: dict[int, tuple[str, ...]] = {
    23: ("C",),
    41: ("W",),
    89: _HYDROPHOBIC,
    104: ("C",),
    118: ("W", "F"),
}

# Last IMGT position expected for a complete domain, by chain type. Heavy chains
# run to 128; kappa and lambda light chains end at 127.
#
# NOTE: IMGT's own pages are inconsistent about where FR4-IMGT ends -- the
# numbering overview gives 118-128, while one amino-acid table gives 118-129.
# 128 is used here because it matches the protein display header and, more
# importantly, because it is what ANARCI emits for every complete heavy chain in
# the test panel.
_LAST_POSITION = {"H": 128, "K": 127, "L": 127}

#: Marks a warning that says *sequence is missing* -- the chain is truncated,
#: a region is empty, an anchor position has no residue at all.
#:
#: Downstream modules key off this prefix to decide whether their own numbers
#: were computed over a reduced chain, so it is load-bearing rather than
#: cosmetic. A warning about a residue being WRONG (an unexpected anchor, an
#: ambiguity code, a framework deletion inside the observed span) must not
#: carry it: nothing is missing in those cases and no count is reduced.
PARTIAL_DOMAIN_PREFIX = "partial V-domain: "

# Thresholds for the "this looks like DNA" check in _clean(). A real V-domain
# is ~110-130 residues, so the length floor is well below any genuine input, and
# a genuine peptide of that length made almost entirely of Ala/Cys/Gly/Thr does
# not occur in antibody variable domains.
NUCLEOTIDE_THRESHOLD = 0.9
NUCLEOTIDE_MIN_LENGTH = 60

#: The 20 standard amino acids.
STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")

#: Ambiguity codes that repertoire sequencing produces. X is any residue; B and Z
#: are Asx and Glx. Tolerated in small amounts, never counted as a match.
AMBIGUITY_CODES = frozenset("XBZ")

#: How much of a sequence may be ambiguity codes before it is refused.
MAX_AMBIGUOUS_FRACTION = 0.02

# Framework positions that IMGT leaves empty by design, per chain type.
#
# CDR positions are freely gapped -- CDR1-IMGT is 5-12 residues, so a short loop
# legitimately leaves 31-34 empty. Framework is fixed-length, so a missing
# framework position is normally a real structural deletion worth flagging. But
# not always: a few framework slots are *always* empty for a given chain type,
# and flagging those would fire on every chain.
#
# Derived empirically from the 1142 complete V-domains in the IMGT/mAb-DB
# therapeutic panel (Merck/BioPhi-2021-publication, MIT); a position is listed
# here if it is absent in >=10% of complete domains of that chain type. In
# practice every position listed was absent in >=95% of them:
#
#   H (n=574): 10 (100%), 73 (100%)
#   K (n=500): 73 (100%), 81 (100%), 82 (100%)
#   L (n=68):  10 (100%), 73 (100%), 81 (96%), 82 (96%)
#
# Position 73 is absent in 100% of all 1142 chains regardless of type. Note that
# OAS's ANARCI_status lists 73 under "Deletions", which overstates it -- it is a
# structural gap. See STRUCTURAL_FRAMEWORK_GAPS.
STRUCTURAL_FRAMEWORK_GAPS: dict[str, frozenset[int]] = {
    "H": frozenset({10, 73}),
    "K": frozenset({73, 81, 82}),
    "L": frozenset({10, 73, 81, 82}),
}

# ANARCI chain_type codes we accept. H = heavy, K = kappa light, L = lambda light.
ChainType = Literal["H", "K", "L"]
_ACCEPTED_CHAINS = ("H", "K", "L")

# IMGT position, with an insertion code. ' ' means no insertion.
Position = tuple[int, str]


class NumberingError(ValueError):
    """Raised when a sequence cannot be numbered as an antibody V-domain."""


class HmmerNotFoundError(RuntimeError):
    """Raised when HMMER's `hmmscan` is not on PATH.

    A separate error because it is not the caller's sequence that is wrong, it
    is the install. ANARCI shells out to `hmmscan`, which is a C program
    installed outside Python, so a correct `uv sync` still leaves it missing and
    the bare `FileNotFoundError: 'hmmscan'` that surfaces says nothing about
    what to do.
    """

    def __init__(self) -> None:
        super().__init__(
            "HMMER's `hmmscan` was not found on PATH. ANARCI needs it to assign "
            "antibody numbering, and it is a system dependency that the Python "
            "lockfile cannot provide.\n"
            "  macOS:  brew install hmmer\n"
            "  Debian: apt-get install hmmer\n"
            "  conda:  conda install -c bioconda hmmer\n"
            "Then check it is visible: hmmscan -h"
        )


@dataclass(frozen=True)
class GermlineCall:
    """Closest germline gene for one segment, as called by ANARCI."""

    species: str
    gene: str
    identity: float

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"{self.gene} ({self.species}, {self.identity:.1%})"


@dataclass(frozen=True)
class Domain:
    """One numbered antibody V-domain.

    Attributes
    ----------
    chain_type:
        "H", "K" or "L" as called by ANARCI.
    species:
        Best-matching species for the numbering HMM. This is *not* a humanness
        score, but it is a useful sanity check: a chain expected to be human that
        numbers best against a mouse HMM should be looked at before anything else
        is computed from it.
    numbering:
        IMGT position -> residue, gaps excluded. Positions carry an insertion code.
    seq:
        The V-domain residues only, in order, with gaps removed. This is the
        substring of the input that was numbered -- not the full input.
    query_start, query_end:
        Zero-based, end-exclusive slice of the *input* sequence that this domain
        occupies. Lets a caller map results back onto a full-length chain.
    v_gene, j_gene:
        Closest germline calls, or None if ANARCI did not assign them.
    """

    chain_type: ChainType
    species: str
    numbering: dict[Position, str]
    seq: str
    query_start: int
    query_end: int
    v_gene: GermlineCall | None = None
    j_gene: GermlineCall | None = None
    hmm_bitscore: float | None = None
    hmm_evalue: float | None = None
    regions: dict[str, str] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    #: Length of the cleaned input this domain was extracted from. Every score
    #: is computed over `seq`, so when this is larger the report has to say so:
    #: otherwise a full-length chain, a tagged construct and a bare V-domain
    #: all print the same way while being scored on different sequences.
    input_length: int = 0

    @property
    def residues_outside_domain(self) -> int:
        """Input residues that were not scored. Zero for a bare V-domain."""
        return max(0, self.input_length - len(self.seq))

    @property
    def is_complete(self) -> bool:
        """True if no completeness warning was raised for this domain.

        A partial domain is still numbered and still scored -- callers are warned,
        not refused -- but any score computed from it is over fewer positions than
        the score it will be compared against, so it must carry the warning
        through to the report.
        """
        return not self.warnings

    @property
    def imgt_span(self) -> tuple[int, int]:
        """Lowest and highest IMGT position present."""
        numbers = [p[0] for p in self.numbering]
        return (min(numbers), max(numbers))

    @property
    def anchors(self) -> dict[int, str | None]:
        """The conserved anchor residues found at each anchor position."""
        return {pos: self.numbering.get((pos, " ")) for pos in ANCHOR_POSITIONS}

    @property
    def is_heavy(self) -> bool:
        return self.chain_type == "H"

    @property
    def is_light(self) -> bool:
        return self.chain_type in ("K", "L")

    def region_of(self, position: Position) -> str | None:
        """Return the IMGT region name containing `position`, or None."""
        return region_of_position(position)

    def framework_seq(self) -> str:
        """Concatenated framework residues, in IMGT order."""
        return "".join(self.regions.get(r, "") for r in FRAMEWORK_REGIONS)

    def cdr_seq(self) -> str:
        """Concatenated CDR residues, in IMGT order."""
        return "".join(self.regions.get(r, "") for r in CDR_REGIONS)

    def positions_in(self, region: str) -> list[Position]:
        """IMGT positions present in this domain that fall in `region`."""
        if region not in IMGT_REGIONS:
            raise KeyError(f"unknown IMGT region: {region!r}")
        lo, hi = IMGT_REGIONS[region]
        return sorted(p for p in self.numbering if lo <= p[0] <= hi)


def region_of_position(position: Position | int) -> str | None:
    """Map an IMGT position to its region name.

    Positions outside 1-128 (which ANARCI does not emit for a V-domain) return None.
    """
    number = position[0] if isinstance(position, tuple) else position
    for name, (lo, hi) in IMGT_REGIONS.items():
        if lo <= number <= hi:
            return name
    return None


def imgt_sort_key(position: Position) -> tuple[int, int]:
    """Order IMGT positions in true sequence order, insertion codes included.

    Sorting naively by (number, letter) is WRONG. IMGT numbers CDR3 outwards from
    the top of the loop, so insertions are added in the order 111, 111A, 111B ...
    and then 112B, 112A, 112 -- that is, insertions at 112 and above run
    *descending* and sit *before* their base position. ANARCI emits them in the
    correct order, e.g.:

        105 106 107 108 109 110 111 111A 112A 112 113 114 115 116 117

    A naive sort puts 112 before 112A and silently transposes residues in exactly
    the chains that matter most -- those with long CDR3 loops.

    Reference: IMGT unique numbering, CDR3-IMGT gap/insertion order
    (https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html).
    """
    number, code = position
    if code == " " or code == "":
        return (number, 0)
    rank = sum((ord(c) - ord("A") + 1) for c in code)
    # Positions up to 111 take insertions after the base; 112 and above take them
    # before, in reverse order.
    return (number, rank) if number <= 111 else (number, -rank)


def _clean(sequence: str) -> str:
    """Strip whitespace and uppercase. Rejects anything outside the residue set."""
    seq = "".join(sequence.split()).upper()
    if not seq:
        raise NumberingError("empty sequence")
    bad = sorted(set(seq) - STANDARD_RESIDUES - AMBIGUITY_CODES)
    if bad:
        raise NumberingError(
            f"sequence contains characters that are not amino acids: {''.join(bad)}. "
            "Gaps, stop codons and other symbols are not supported."
        )

    # Ambiguity codes are tolerated up to a limit rather than refused outright.
    # Repertoire sequencing routinely produces them: rejecting any sequence
    # containing an X made 3,469 of 17,786 real OAS records (19.5%) unscorable,
    # which is disproportionate for one uncertain residue in ~120.
    #
    # That limit is NOT applied here. It belongs to the extracted domain, which
    # does not exist yet at this point, so it is enforced in
    # `_reject_ambiguous_domains` after numbering. Applied to the raw input it
    # measured the wrong sequence: see that function.

    # A nucleotide sequence passes the residue check above, because A, C, G and T
    # are all valid amino acid letters (Ala, Cys, Gly, Thr). Without this check
    # it fails later as "no antibody V-domain found", which points at the wrong
    # problem. Translation is deliberately not offered: a nucleotide sequence has
    # three reading frames and picking one silently would mean returning a
    # confident score from a guess.
    nucleotide_fraction = sum(seq.count(base) for base in "ACGTUN") / len(seq)
    if len(seq) >= NUCLEOTIDE_MIN_LENGTH and nucleotide_fraction >= NUCLEOTIDE_THRESHOLD:
        raise NumberingError(
            f"looks like a nucleotide sequence ({nucleotide_fraction:.0%} ACGTU/N over "
            f"{len(seq)} characters). This tool scores amino acid sequences; translate "
            "it in the correct reading frame first."
        )
    return seq


def completeness_warnings(
    numbering: dict[Position, str],
    regions: dict[str, str],
    chain_type: str,
) -> tuple[str, ...]:
    """Check whether a numbered domain looks like a complete V-domain.

    Returns a tuple of human-readable warnings; empty means complete. Nothing
    here refuses a sequence -- a partial domain is a legitimate input (a reviewer
    pasting a fragment, a sequencing read that did not reach the end) and gets
    scored. The caller's job is to carry these warnings into the report so that a
    number computed over part of a domain is never mistaken for one computed over
    all of it.

    Four independent signals, because no one of them is sufficient:

    1. Conserved anchors (23, 41, 89, 104, 118). Biological, and the strongest test.
    2. Empty regions. Catches C-terminal truncation.
    3. Span. Catches truncation at either end. Signal 2 misses N-terminal
       truncation entirely: dropping the first 20 residues of a chain still
       leaves every region populated, because FR1 is 26 positions long.
    4. Framework deletions -- gaps *interior* to the observed span. Framework is
       fixed-length, so an interior gap is a structural deletion rather than a
       short loop, except at the slots IMGT leaves empty per chain type.

    Signals 3 and 4 must not overlap, and getting that wrong was a real bug:
    scanning for framework deletions past the last observed residue reported
    every missing FR4 position as a deletion, so a read that simply stopped early
    came back with a dozen spurious "deletion" warnings on top of the truncation
    warning that already described it. Marketed therapeutics are complete, so the
    therapeutic panel could not surface this; it took 14,317 OAS repertoire
    sequences, where C-terminal truncation is routine, for it to appear.
    """
    messages: list[str] = []

    ambiguous = sorted(
        (position for position, residue in numbering.items() if residue in AMBIGUITY_CODES),
        key=imgt_sort_key,
    )
    if ambiguous:
        messages.append(
            "sequence contains "
            + str(len(ambiguous))
            + " ambiguity code(s) at IMGT position(s) "
            + ", ".join(f"{p[0]}{p[1].strip()}" for p in ambiguous)
            + ". These never count as a germline match and never occur in a human "
            "peptide, so both scores are depressed slightly."
        )

    missing_anchors = []
    wrong_anchors = []
    for position, expected in ANCHOR_POSITIONS.items():
        residue = numbering.get((position, " "))
        if residue is None:
            missing_anchors.append(position)
        elif residue not in expected:
            wrong_anchors.append(f"{position}={residue} (expected {'/'.join(expected)})")

    if missing_anchors:
        messages.append(
            PARTIAL_DOMAIN_PREFIX
            + "missing conserved IMGT anchor position(s) "
            + ", ".join(str(p) for p in missing_anchors)
            + ". Scores are computed over the positions present only."
        )
    if wrong_anchors:
        messages.append(
            "unexpected residue at conserved IMGT anchor position(s): "
            + "; ".join(wrong_anchors)
            + ". The numbering may be misaligned."
        )

    empty = [name for name in IMGT_REGIONS if not regions.get(name)]
    if empty:
        messages.append(PARTIAL_DOMAIN_PREFIX + "no residues numbered in " + ", ".join(empty) + ".")

    # Deletions inside framework. Framework is fixed-length, so a gap there is a
    # structural deletion rather than a short loop -- except at the positions
    # IMGT leaves empty by design for this chain type.
    #
    # The scan is bounded by the FIRST and LAST observed positions. Anything
    # outside that span is truncation, which the span checks below report on
    # their own; counting it here as well both mislabels it and double-reports it.
    present = {p[0] for p in numbering}
    structural = STRUCTURAL_FRAMEWORK_GAPS.get(chain_type, frozenset())
    first_observed = min(present, default=1)
    last_observed = max(present, default=0)
    framework_deletions = sorted(
        position
        for region in FRAMEWORK_REGIONS
        for position in range(IMGT_REGIONS[region][0], IMGT_REGIONS[region][1] + 1)
        if first_observed <= position <= last_observed
        and position not in present
        and position not in structural
    )
    if framework_deletions:
        messages.append(
            "deletion inside framework at IMGT position(s) "
            + ", ".join(str(p) for p in framework_deletions)
            + ". Framework is fixed-length, so this reads as a structural deletion "
            "rather than a short loop. Per-position comparisons will be offset "
            "against germline here."
        )

    numbers = [p[0] for p in numbering]
    if numbers:
        first, last = min(numbers), max(numbers)
        if first > 1:
            messages.append(
                f"{PARTIAL_DOMAIN_PREFIX}numbering starts at IMGT {first} instead of 1 "
                f"({first - 1} N-terminal position(s) absent)."
            )
        expected_last = _LAST_POSITION.get(chain_type)
        if expected_last is not None and last < expected_last:
            messages.append(
                f"{PARTIAL_DOMAIN_PREFIX}numbering ends at IMGT {last}, expected "
                f"{expected_last} for chain type {chain_type}."
            )

    return tuple(messages)


def _germline_call(details: dict, segment: str) -> GermlineCall | None:
    """Pull one germline segment call out of an ANARCI details dict."""
    germlines = details.get("germlines") or {}
    entry = germlines.get(segment)
    if not entry:
        return None
    # ANARCI shape: [(species, gene), identity]
    try:
        (species, gene), identity = entry[0], entry[1]
    except (TypeError, ValueError, IndexError):
        return None
    return GermlineCall(species=str(species), gene=str(gene), identity=float(identity))


def number_sequence(
    sequence: str,
    *,
    scheme: str = "imgt",
    assign_germline: bool = True,
    allowed_species: list[str] | None = None,
) -> list[Domain]:
    """Number every antibody V-domain found in `sequence`.

    A full-length heavy chain contains one V-domain followed by constant domains;
    only the V-domain is numbered and returned, so pasting a full-length chain is
    fine. A single-chain construct (e.g. an scFv) yields more than one domain.

    Raises
    ------
    NumberingError
        If no antibody V-domain is found, or the sequence is not clean protein.
    """
    seq = _clean(sequence)

    kwargs: dict = {"scheme": scheme, "assign_germline": assign_germline}
    if allowed_species is not None:
        kwargs["allowed_species"] = allowed_species

    try:
        _, numbered, details, _ = run_anarci([("query", seq)], **kwargs)
    except FileNotFoundError as exc:
        # ANARCI invokes `hmmscan` by name; a missing binary arrives here as a
        # bare FileNotFoundError with no guidance.
        if "hmmscan" in str(exc):
            raise HmmerNotFoundError() from exc
        raise

    domains_numbered = numbered[0]
    domains_details = details[0]
    if not domains_numbered:
        raise NumberingError(
            "no antibody V-domain found. Expected an IgG heavy (VH) or light "
            "(VL kappa/lambda) variable domain sequence."
        )

    domains: list[Domain] = []
    excluded: list[str] = []
    # The metadata list is parallel to the RAW domain list, so it is indexed by
    # the raw position. An earlier version indexed it by `len(domains)` -- the
    # count of domains ACCEPTED so far -- which stays put whenever a domain is
    # skipped. One skipped TCR therefore shifted every later domain onto the
    # wrong metadata, and a TCR followed by an antibody was rejected as "not an
    # antibody" while the same two domains in the other order scored fine. A
    # result that depends on the order of the input is the symptom to watch for.
    for index, (numbering, start, end) in enumerate(domains_numbered):
        detail = domains_details[index] if index < len(domains_details) else {}
        chain_type = detail.get("chain_type")
        if chain_type not in _ACCEPTED_CHAINS:
            # ANARCI also numbers TCR chains (A/B/G/D). Those are out of scope,
            # but a domain dropped in silence is a domain the caller does not
            # know was there, so each one is recorded and disclosed below.
            excluded.append(str(chain_type))
            continue

        occupied = {pos: res for pos, res in numbering if res != "-"}
        regions: dict[str, str] = {}
        for name, (lo, hi) in IMGT_REGIONS.items():
            residues = [
                occupied[p] for p in sorted(occupied, key=imgt_sort_key) if lo <= p[0] <= hi
            ]
            regions[name] = "".join(residues)

        domains.append(
            Domain(
                chain_type=chain_type,
                species=str(detail.get("species", "unknown")),
                numbering=occupied,
                seq="".join(occupied[p] for p in sorted(occupied, key=imgt_sort_key)),
                query_start=int(start),
                query_end=int(end) + 1,
                v_gene=_germline_call(detail, "v_gene"),
                j_gene=_germline_call(detail, "j_gene"),
                hmm_bitscore=detail.get("bitscore"),
                hmm_evalue=detail.get("evalue"),
                regions=regions,
                warnings=completeness_warnings(occupied, regions, chain_type),
            )
        )

    if not domains:
        raise NumberingError(
            "a variable domain was found but it is not an antibody heavy or light "
            "chain (ANARCI called a non-IG chain type)."
        )

    _reject_ambiguous_domains(domains)
    domains = [replace(d, input_length=len(seq)) for d in domains]

    if excluded:
        # Policy: score the antibody domain and say what was left out. The
        # alternative -- refusing any construct that also holds a non-antibody
        # domain -- would reject a legitimate input over a neighbour the tool
        # has nothing to say about. Either policy is defensible; a policy that
        # depends on which domain came first is not.
        note = (
            f"input also contains {len(excluded)} non-antibody variable domain(s) "
            f"(ANARCI chain type {', '.join(sorted(set(excluded)))}). They are "
            "outside this tool's scope and were excluded from every score."
        )
        domains = [replace(d, warnings=(*d.warnings, note)) for d in domains]

    return domains


def _reject_ambiguous_domains(domains: list[Domain]) -> None:
    """Apply the ambiguity limit to each EXTRACTED domain.

    The limit belongs to the sequence that gets scored, not to whatever was
    pasted in around it. Checked on the whole input it was trivially bypassed:
    a VH with 12 of its 120 residues replaced by X was refused at 10%, and the
    identical domain behind a 650-residue unrelated prefix passed at 1.5% --
    scoring exactly the 120 residues that had just been called too uncertain.
    The denominator has to be the thing the number is computed over.

    The converse case changes too, and for the same reason: a clean domain
    inside a long ambiguous construct is now scored, because uncertainty
    outside the extracted domain does not make the extracted domain uncertain.
    """
    for domain in domains:
        n_ambiguous = sum(domain.seq.count(code) for code in AMBIGUITY_CODES)
        if not n_ambiguous:
            continue
        fraction = n_ambiguous / len(domain.seq)
        if fraction > MAX_AMBIGUOUS_FRACTION:
            raise NumberingError(
                f"the extracted {domain.chain_type} domain is {fraction:.1%} ambiguity "
                f"codes ({n_ambiguous} of {len(domain.seq)} residues), above the "
                f"{MAX_AMBIGUOUS_FRACTION:.1%} limit. Per-position comparison would "
                "not be meaningful. The limit applies to the scored domain, so "
                "padding the input does not change this."
            )


def number_single_chain(sequence: str, **kwargs) -> Domain:
    """Number a sequence expected to contain exactly one antibody V-domain.

    Raises NumberingError if more than one domain is present, so that a caller
    cannot silently score only the first domain of an scFv.
    """
    domains = number_sequence(sequence, **kwargs)
    if len(domains) > 1:
        kinds = ", ".join(d.chain_type for d in domains)
        raise NumberingError(
            f"found {len(domains)} V-domains ({kinds}). Pass one chain at a time, "
            "or use number_sequence() to handle each domain explicitly."
        )
    return domains[0]
