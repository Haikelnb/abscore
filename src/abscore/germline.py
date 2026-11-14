"""Closest human germline V gene, and per-position identity to it.

This is the plain baseline that the repertoire-based humanness score (see
`humanness.py`) has to beat to justify itself. If a therapeutic panel separates
by origin almost as well on germline identity as on a 9-mer repertoire score,
that is the honest finding and it belongs in the report.

What is compared, and why not more
----------------------------------
A V gene does not encode the whole variable domain, so germline identity cannot
be computed over the whole domain. Measured against ANARCI's bundled IMGT
germline set:

- human V alleles carry residues up to IMGT 106-108 -- the `AR`/`AA` that opens
  CDR3 is V-encoded;
- human J alleles start at IMGT 115, so the J gene encodes FR4 (118-128);
- IMGT 107-114 is encoded by neither: it is junctional, built from the D segment
  and from nucleotides added at random during V(D)J recombination.

So "framework identity to the V gene" means FR1 + FR2 + FR3. Including FR4 would
be comparing a V gene against a region the V gene does not encode. Rather than
hard-coding a cut-off, the germline allele's own gap pattern is taken as the
authority for which positions it covers.

Denominators
------------
Two are reported, because for a partial domain they differ and only one answers
the question the caller asked:

- `identity`           matches / positions present in BOTH query and germline.
                       The primary number: "of what I can see, how germline is
                       it?" This is also what OAS's `v_identity` and ANARCI's own
                       germline identity report.
- `identity_over_germline`
                       matches / all germline positions in scope, counting
                       absent query positions as mismatches. Needed when every
                       member of a panel must share a denominator.
- `coverage`           positions compared / germline positions in scope.

For a complete domain all of this collapses: coverage is 1.0 and the two
identities are equal. They diverge only for partial domains, which is exactly
when a single number would mislead -- a truncated but perfectly human chain
scores low on `identity_over_germline` because of its length, rather than because it is
non-human. Report `identity` with `coverage` beside it, and never compare
identities across different coverage without saying so.

References
----------
- IMGT unique numbering: Lefranc M-P et al., Dev Comp Immunol 2003;27:55-77.
- ANARCI (germline set and numbering): Dunbar J, Deane CM.
  Bioinformatics 2016;32:298-300.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import cache

from anarci.germlines import all_germlines

from .number import AMBIGUITY_CODES, CDR_REGIONS, IMGT_REGIONS, Domain, Position

# Scopes that germline identity is reported over. Each is a set of IMGT region
# names; the actual positions compared are these regions intersected with the
# positions the germline allele itself covers.
#
# FR4 is deliberately absent from "framework": it is encoded by the J gene.
# CDR3 is deliberately absent from "cdr": its origin is junctional.
SCOPES: dict[str, tuple[str, ...]] = {
    "V": ("FR1", "CDR1", "FR2", "CDR2", "FR3"),
    "framework": ("FR1", "FR2", "FR3"),
    "cdr": ("CDR1", "CDR2"),
}

# ANARCI ships no IMGT release number for its germline set, so the data is pinned
# by content hash instead. A changed hash means the germline reference moved and
# any stored identity is no longer comparable.
_EXPECTED_GERMLINE_SHA256 = "91ad0e1177cc52c429d38dc6821422682de22ea24303d4da76f7a3ff48df88e0"


class GermlineError(ValueError):
    """Raised when no germline allele is available for a chain."""


def germline_data_sha256() -> str:
    """Content hash of ANARCI's bundled germline set, for provenance."""
    blob = json.dumps(all_germlines, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def germline_provenance() -> dict[str, str | bool]:
    """What produced the germline numbers, for `run.json`."""
    from importlib.metadata import version

    digest = germline_data_sha256()
    return {
        "germline_source": "ANARCI bundled IMGT germline set (anarci.germlines.all_germlines)",
        "anarci_version": version("anarci"),
        "germline_data_sha256": digest,
        "germline_data_matches_pin": digest == _EXPECTED_GERMLINE_SHA256,
        "imgt_release": "not recorded by ANARCI",
    }


@cache
def v_germlines(chain_type: str, species: str = "human") -> dict[str, str]:
    """Germline V alleles for a chain type, as IMGT-gapped 128-character strings.

    Every allele is already in IMGT numbering, so position i of the string is
    IMGT position i+1. That is what makes identity a slot-by-slot comparison
    rather than an alignment with gap penalties to choose.
    """
    try:
        alleles = all_germlines["V"][chain_type][species]
    except KeyError as exc:
        raise GermlineError(
            f"no germline V set for chain_type={chain_type!r} species={species!r}"
        ) from exc
    return dict(alleles)


def germline_positions(allele_sequence: str) -> dict[int, str]:
    """IMGT position -> residue for the positions this allele actually covers."""
    return {i + 1: residue for i, residue in enumerate(allele_sequence) if residue != "-"}


def _scope_positions(scope: str) -> set[int]:
    try:
        regions = SCOPES[scope]
    except KeyError as exc:
        raise KeyError(f"unknown scope {scope!r}; expected one of {sorted(SCOPES)}") from exc
    positions: set[int] = set()
    for region in regions:
        low, high = IMGT_REGIONS[region]
        positions.update(range(low, high + 1))
    return positions


@dataclass(frozen=True)
class ScopedIdentity:
    """Identity to one germline allele over one scope.

    `compared` counts positions present in both query and germline; `covered`
    counts the germline positions in scope. They differ only for a partial
    domain.
    """

    scope: str
    matches: int
    compared: int
    covered: int

    @property
    def identity(self) -> float | None:
        """matches / positions compared. The primary number. None if nothing overlapped."""
        if self.compared == 0:
            return None
        return self.matches / self.compared

    @property
    def identity_over_germline(self) -> float | None:
        """matches / germline positions in scope; absences count as mismatches."""
        if self.covered == 0:
            return None
        return self.matches / self.covered

    @property
    def coverage(self) -> float | None:
        """Fraction of the germline's positions in scope that the query covers."""
        if self.covered == 0:
            return None
        return self.compared / self.covered

    @property
    def is_full_coverage(self) -> bool:
        return self.compared == self.covered and self.covered > 0

    def describe(self) -> str:
        """One line, with the denominator stated. Never print a bare percentage."""
        if self.identity is None:
            return f"{self.scope}: no overlap with germline"
        text = f"{self.scope}: {self.identity:.1%} ({self.matches}/{self.compared} positions"
        if not self.is_full_coverage:
            text += f", coverage {self.coverage:.1%} of {self.covered} germline positions"
        return text + ")"


@dataclass(frozen=True)
class Mismatch:
    """One position where the query differs from, or is absent against, germline."""

    position: Position
    region: str
    germline_residue: str
    query_residue: str | None  # None means the position is absent from the query

    @property
    def is_absence(self) -> bool:
        return self.query_residue is None


@dataclass(frozen=True)
class GeneContender:
    """A gene whose best allele comes close enough to the winner to matter."""

    gene: str
    best_allele: str
    identity: float
    positions_behind: float

    def __str__(self) -> str:  # pragma: no cover - display only
        return f"{self.best_allele} ({self.identity:.1%}, {self.positions_behind:.0f} behind)"


@dataclass(frozen=True)
class GermlineAssignment:
    """The closest germline V allele for a domain, with identities and mismatches."""

    gene: str
    species: str
    chain_type: str
    scopes: dict[str, ScopedIdentity]
    mismatches: tuple[Mismatch, ...]
    runner_up: str | None = None
    runner_up_identity: float | None = None
    n_alleles_considered: int = 0
    contenders: tuple[GeneContender, ...] = ()
    #: Best-fitting allele of the nearest gene that is NOT the called one, kept
    #: regardless of the ambiguity threshold. `contenders` holds only ties, so
    #: without this the margin to the next gene would be invisible as soon as
    #: the call is unambiguous -- which is exactly when a reader wants it.
    nearest_other_gene: GeneContender | None = None

    @property
    def identity(self) -> float | None:
        """Identity over the V scope -- the headline number."""
        return self.scopes["V"].identity

    @property
    def framework_identity(self) -> float | None:
        """Identity over FR1+FR2+FR3.

        Arguably the fairer measure of "is the scaffold human?", because CDRs are
        *supposed* to differ from germline -- that is what affinity maturation
        does -- so penalising CDR divergence partly penalises the antibody for
        working.
        """
        return self.scopes["framework"].identity

    @property
    def margin(self) -> float | None:
        """Identity lead over the runner-up allele.

        A small margin means the gene call is not decided by the data. Alleles of
        one gene differ by a residue or two, so a near-tie between alleles of the
        same gene is expected and harmless; a near-tie between different genes is
        not, and callers should surface it.
        """
        if self.identity is None or self.runner_up_identity is None:
            return None
        return self.identity - self.runner_up_identity

    @property
    def gene_without_allele(self) -> str:
        return self.gene.split("*")[0]

    @property
    def runner_up_is_same_gene(self) -> bool:
        if self.runner_up is None:
            return False
        return self.runner_up.split("*")[0] == self.gene_without_allele

    @property
    def margin_positions(self) -> float | None:
        """The runner-up lead expressed in positions rather than percent."""
        if self.margin is None:
            return None
        return self.margin * self.scopes["V"].compared

    @property
    def ambiguous_genes(self) -> tuple[GeneContender, ...]:
        """Other *genes* close enough that the assignment could flip.

        Alleles of one gene differ by a residue or two, so a near-tie between
        alleles of the same gene is expected and harmless -- the gene is still
        called correctly. A near-tie between different genes is not: the
        assignment, and every per-position mismatch derived from it, could change
        on one residue.

        Note that looking only at the single next-best allele is not enough. The
        immediate runner-up is often another allele of the winning gene, which
        hides a different gene behind it, so alleles are collapsed to genes
        before the threshold is applied.

        At the default threshold this reports exact ties only. Trastuzumab shows
        why that is the right line: IGHV3-53 trails it by one position, and the
        two genes differ from each other at exactly one position (IMGT 81, D in
        IGHV3-66 and H in IGHV3-53), where trastuzumab has D. The call rests on
        one residue, but it rests on it correctly, so the gene is reported
        without an ambiguity caveat and the margin is available separately.
        """
        return tuple(c for c in self.contenders if c.gene != self.gene_without_allele)

    @property
    def is_ambiguous(self) -> bool:
        """True when another gene is within the ambiguity threshold."""
        return bool(self.ambiguous_genes)

    def warnings(self) -> tuple[str, ...]:
        """Caveats that must travel with these numbers into a report."""
        messages: list[str] = []
        if self.ambiguous_genes:
            others = ", ".join(
                f"{c.best_allele} ({c.positions_behind:.0f} position(s) behind)"
                for c in self.ambiguous_genes
            )
            messages.append(
                f"germline gene call is ambiguous: {self.gene} cannot be separated from "
                f"{len(self.ambiguous_genes)} other gene(s) -- {others}. Per-position "
                "mismatches could change with a different assignment."
            )
        v = self.scopes["V"]
        if not v.is_full_coverage and v.coverage is not None:
            messages.append(
                f"germline identity computed at {v.coverage:.1%} coverage "
                f"({v.compared} of {v.covered} germline positions). Not comparable "
                "with identities computed at full coverage."
            )
        return tuple(messages)

    def mismatches_in_cdrs(self) -> tuple[Mismatch, ...]:
        return tuple(m for m in self.mismatches if m.region in CDR_REGIONS)

    def mismatches_in_framework(self) -> tuple[Mismatch, ...]:
        return tuple(m for m in self.mismatches if m.region not in CDR_REGIONS)


def compare_to_germline(
    domain: Domain,
    gene: str,
    *,
    species: str = "human",
) -> GermlineAssignment:
    """Score `domain` against one named germline V allele."""
    alleles = v_germlines(domain.chain_type, species)
    if gene not in alleles:
        raise GermlineError(
            f"germline {gene!r} not in the {species} {domain.chain_type} V set "
            f"({len(alleles)} alleles available)"
        )
    germline = germline_positions(alleles[gene])

    # Insertion codes cannot align to a germline slot: the germline set has no
    # insertions, so a query insertion has no counterpart to compare against.
    # Those positions are excluded from every scope rather than counted as
    # mismatches, since the absence sits in the reference.
    query = {p[0]: residue for p, residue in domain.numbering.items() if p[1] == " "}

    scopes: dict[str, ScopedIdentity] = {}
    for scope in SCOPES:
        in_scope = _scope_positions(scope)
        covered = sorted(p for p in germline if p in in_scope)
        compared = [p for p in covered if p in query]
        # An ambiguity code never matches, even against itself: X means "we do
        # not know" rather than "the same unknown residue".
        matches = sum(
            1 for p in compared if query[p] == germline[p] and query[p] not in AMBIGUITY_CODES
        )
        scopes[scope] = ScopedIdentity(
            scope=scope, matches=matches, compared=len(compared), covered=len(covered)
        )

    v_scope = _scope_positions("V")
    mismatches = tuple(
        Mismatch(
            position=(position, " "),
            region=_region_name(position),
            germline_residue=germline[position],
            query_residue=query.get(position),
        )
        for position in sorted(germline)
        if position in v_scope and query.get(position) != germline[position]
    )

    return GermlineAssignment(
        gene=gene,
        species=species,
        chain_type=domain.chain_type,
        scopes=scopes,
        mismatches=mismatches,
        n_alleles_considered=1,
    )


def assign_closest_germline(
    domain: Domain,
    *,
    species: str = "human",
    scope: str = "V",
    ambiguity_positions: float = 1.0,
) -> GermlineAssignment:
    """Find the germline V allele closest to `domain` and score against it.

    Selection is on identity over `scope`, with the number of positions compared
    as the tie-break, so that an allele which covers more of the query wins a tie
    against one that happens to match over fewer positions.

    `ambiguity_positions` sets how close another gene must come before the call
    is reported as ambiguous, in positions rather than percent -- a percentage
    threshold would mean different numbers of residues on domains of different
    length, and for a partial domain it would mean fewer.

    The default of 1.0 means **ties only**: a gene is a contender only if it
    matches at exactly as many positions as the winner. A gene that is behind by
    even one position is not reported as ambiguous, because it is not: 79
    matches beats 78, and the better fit really is the better fit.

    An earlier default of 2.0 also flagged a gene trailing by a single position,
    which conflated two different facts. "Cannot be separated" should mean
    literally cannot. A thin winning margin is worth knowing too, but it is a
    number -- see `margin_positions` -- and reporting it as a caveat overstated
    the uncertainty and crowded out the caveats that mattered.
    """
    alleles = v_germlines(domain.chain_type, species)
    if not alleles:
        raise GermlineError(f"no {species} {domain.chain_type} V alleles available")

    scored: list[tuple[float, int, str]] = []
    for gene in alleles:
        result = compare_to_germline(domain, gene, species=species)
        identity = result.scopes[scope].identity
        if identity is None:
            continue
        scored.append((identity, result.scopes[scope].compared, gene))

    if not scored:
        raise GermlineError(
            f"no germline allele overlapped the domain in scope {scope!r}; the "
            "domain may not cover any V-encoded position"
        )

    scored.sort(reverse=True)
    best_identity, _, best_gene = scored[0]
    assignment = compare_to_germline(domain, best_gene, species=species)

    runner_up = runner_up_identity = None
    for identity, _, gene in scored[1:]:
        if gene != best_gene:
            runner_up, runner_up_identity = gene, identity
            break

    # Collapse alleles to genes, keeping each gene's best allele, then keep the
    # genes within the threshold. One allele of the winning gene is always kept,
    # so `contenders` always includes the winner itself.
    compared = assignment.scopes[scope].compared
    best_per_gene: dict[str, tuple[float, str]] = {}
    for identity, _, gene in scored:
        name = gene.split("*")[0]
        if name not in best_per_gene or identity > best_per_gene[name][0]:
            best_per_gene[name] = (identity, gene)

    by_gene = sorted(
        (
            GeneContender(
                gene=name,
                best_allele=allele,
                identity=identity,
                positions_behind=(best_identity - identity) * compared,
            )
            for name, (identity, allele) in best_per_gene.items()
        ),
        key=lambda c: (-c.identity, c.gene),
    )
    called_gene = assignment.gene.split("*")[0]
    # Compared with a tolerance, because `positions_behind` is a float derived
    # from an identity difference: a gene exactly one position behind came out
    # as 0.9999999999999964 and slipped under a threshold of 1.0, so "ties only"
    # silently still admitted one-position gaps.
    contenders = tuple(c for c in by_gene if c.positions_behind < ambiguity_positions - 1e-6)
    nearest_other = next((c for c in by_gene if c.gene != called_gene), None)

    return GermlineAssignment(
        gene=assignment.gene,
        species=assignment.species,
        chain_type=assignment.chain_type,
        scopes=assignment.scopes,
        mismatches=assignment.mismatches,
        runner_up=runner_up,
        runner_up_identity=runner_up_identity,
        n_alleles_considered=len(scored),
        contenders=contenders,
        nearest_other_gene=nearest_other,
    )


def _region_name(position: int) -> str:
    for name, (low, high) in IMGT_REGIONS.items():
        if low <= position <= high:
            return name
    return "?"
