"""Repertoire-based humanness: how much of a chain do real humans actually make?

This is the OASis method (Prihoda et al., MAbs 2022), computed through `promb`.
Where `germline.py` asks "how close is this to a human germline *gene*", this asks
"has any human *actually produced* these peptides". The two differ in a way that
matters, and the difference is the point of reporting both:

- A germline gene does not cover CDR3, so germline identity is structurally blind
  to the most non-human part of the molecule.
- Somatic hypermutation moves every real human antibody away from its germline,
  so germline identity penalises a well-matured, perfectly human antibody.
- Only the peptide view sees boundary-spanning peptides -- the 9-mers straddling
  a CDR/framework seam. Those can be sequences no human has produced while every
  residue on either side is human. Averaged identity dilutes them to nothing;
  peptide lookup finds them.

Method
------
Slide a 9-mer window along the V-domain. For each peptide, ask whether it occurs
in human antibody repertoires. The score is the fraction that do -- OASis
identity.

The reference is a precomputed table rather than a live query: BioPhi derived it once
from 118,713,869 OAS sequences across 231 donors in 26 studies (OAS snapshot
November 2019), recording for each 9-mer the fraction of *donors* whose
repertoire contains it.

Donors are the unit, deliberately: one donor's clonal expansion can
contribute thousands of near-identical sequences, so a sequence-weighted count
would let a single person's immune history dominate. Counting donors asks the
question actually meant -- how many separate human immune systems have made this
peptide.

Threshold
---------
BioPhi defines four prevalence thresholds: loose (>=1% of donors), relaxed
(>=10%), medium (>=50%) and strict (>=90%). **`promb` ships only the relaxed
set**, so that is what this module computes. The paper's headline correlation
with clinical ADA (R = -0.53) used *medium*, which needs BioPhi's full 22 GB
database. Any comparison against published OASis numbers has to name the
threshold, so the threshold travels in the provenance dict and into the report.

References
----------
- Prihoda D et al. BioPhi. MAbs 2022;14(1):2020203. doi:10.1080/19420862.2021.2020203
- Kovaltsuk A et al. Observed Antibody Space. J Immunol 2018;201:2502-2509.
- promb (MIT): https://github.com/MSDLLCpapers/promb
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from functools import cache

from .number import (
    AMBIGUITY_CODES,
    CDR_REGIONS,
    PARTIAL_DOMAIN_PREFIX,
    Domain,
    Position,
    imgt_sort_key,
    region_of_position,
)

#: Reference databases `promb` ships. "human-oas" is antibody repertoire data and
#: is the right one for humanness. The other two are general human protein
#: sequence, and are here because `epitopes.py` needs a broader definition of
#: "self" than Ig germline alone -- a peptide identical to a stretch of human
#: albumin is also self.
DATABASES = {
    "human-oas": "human antibody repertoires (OAS), peptides in >=10% of donors",
    "human-swissprot": "human SwissProt proteome",
    "human-reference": "human reference proteome UP000005640",
}

#: The donor-prevalence threshold behind the bundled "human-oas" table.
#: BioPhi's name for it, and the numeric cutoff.
OAS_THRESHOLD_NAME = "relaxed"
OAS_THRESHOLD_MIN_DONOR_FRACTION = 0.10


class HumannessError(ValueError):
    """Raised when a chain cannot be scored."""


#: Databases promb ships as whole sequences rather than as a peptide list. These
#: need an explicit peptide length at load time, because they are chopped on the
#: fly; "human-oas" is already a list of 9-mers and infers its own.
FASTA_DATABASES = ("human-swissprot", "human-reference")


@cache
def get_database(name: str = "human-oas", peptide_length: int | None = None):
    """Load a `promb` peptide database, once per (name, peptide_length).

    `promb` prints a banner and a load message on import and on load; both are
    swallowed so that library use stays quiet and CLI output stays parseable.
    """
    if name not in DATABASES:
        raise HumannessError(f"unknown database {name!r}; expected one of {sorted(DATABASES)}")
    if name in FASTA_DATABASES and peptide_length is None:
        raise HumannessError(
            f"{name!r} is a whole-sequence database, so peptide_length is required "
            "(peptides are chopped on load). Only 'human-oas' infers its own length."
        )
    from promb import init_db

    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        if peptide_length is None:
            return init_db(name, verbose=False)
        return init_db(name, peptide_length=peptide_length, verbose=False)


@dataclass(frozen=True)
class Peptide:
    """One sliding-window peptide, located in IMGT coordinates.

    `positions` is what makes the result reportable: it lets a non-human peptide
    be placed in a region, so the report can say *where* a chain stops looking
    human rather than only how much of it does.
    """

    sequence: str
    start_index: int
    positions: tuple[Position, ...]
    found: bool

    @property
    def regions(self) -> tuple[str, ...]:
        """Distinct IMGT regions this peptide spans, in order."""
        seen: list[str] = []
        for position in self.positions:
            region = region_of_position(position)
            if region is not None and region not in seen:
                seen.append(region)
        return tuple(seen)

    @property
    def spans_boundary(self) -> bool:
        """True if the peptide straddles a framework/CDR boundary.

        A geometric fact. It is tempting to read a
        boundary-spanning non-human peptide as a graft junction -- the seam where
        a murine CDR meets a human framework -- but position does not establish
        cause. A frameshifted window, a heavily mutated loop and a genuine graft
        seam all produce boundary-spanning peptides, and on a real OAS record this
        property labelled a 14-peptide frameshift artefact. Nothing downstream
        interprets it, and the report states which positions are flagged rather
        than why.
        """
        regions = self.regions
        if len(regions) < 2:
            return False
        return any(r in CDR_REGIONS for r in regions) and any(r not in CDR_REGIONS for r in regions)

    @property
    def touches_cdr(self) -> bool:
        return any(r in CDR_REGIONS for r in self.regions)


@dataclass(frozen=True)
class HumannessScore:
    """OASis-style humanness for one chain, with every peptide retained."""

    peptides: tuple[Peptide, ...]
    database: str
    peptide_length: int
    threshold_name: str
    threshold_min_donor_fraction: float
    chain_type: str
    domain_warnings: tuple[str, ...] = ()

    @property
    def n_peptides(self) -> int:
        return len(self.peptides)

    @property
    def n_found(self) -> int:
        return sum(1 for p in self.peptides if p.found)

    @property
    def score(self) -> float | None:
        """Fraction of peptides present in the reference. OASis identity."""
        if not self.peptides:
            return None
        return self.n_found / self.n_peptides

    @property
    def non_human_peptides(self) -> tuple[Peptide, ...]:
        return tuple(p for p in self.peptides if not p.found)

    @property
    def boundary_spanning_peptides(self) -> tuple[Peptide, ...]:
        """Non-human peptides that straddle a framework/CDR boundary."""
        return tuple(p for p in self.non_human_peptides if p.spans_boundary)

    def flagged_runs(self) -> tuple[tuple[Position, Position, int, tuple[str, ...]], ...]:
        """Contiguous stretches covered only by peptides absent from the reference.

        Reported as (first position, last position, length, regions spanned).

        Three limits on how far this can be read, all of them real:

        *Absent from the reference is not "non-human".* The reference holds
        peptides meeting the relaxed threshold -- present in at least 10% of OAS
        donors -- so a peptide one person carries and most do not is absent from
        it while being perfectly human.

        *A flagged position is covered by at least one unmatched window, not
        only by unmatched windows.* Each residue sits in up to nine 9-mers, and
        the flag marks the position as participating in an unmatched one.

        *The residues in a run are not shown to be the cause.* A run is an
        observation about where the sequence stops matching the reference, with
        no claim about why. A long run and a single flagged position are
        different findings and should not be summarised as one count.

        Adjacency is judged in SEQUENCE order, never by IMGT position arithmetic.
        IMGT numbering is not dense -- a complete domain legitimately leaves 8-22
        of its 128 slots empty, and a chain may also carry a real deletion. Two
        residues that sit side by side in the protein can therefore be numbered
        72 and 75. Testing `position - previous <= 1` splits one stretch into
        several: on a frameshifted OAS record it reported 62-72 and 75-85 as
        separate findings when they are one 22-residue stretch, because IMGT 73
        (a structural gap) and 74 (a deletion) are absent.
        """
        flagged = self.flags_by_position()
        if not flagged:
            return ()

        # Index every numbered position by its place in the sequence, so that
        # "next residue along" means the next residue in the chain.
        order = {
            position: index
            for index, position in enumerate(
                sorted(
                    {p for peptide in self.peptides for p in peptide.positions}, key=imgt_sort_key
                )
            )
        }
        ordered = sorted(flagged, key=lambda p: order[p])

        runs: list[tuple[Position, Position, int, tuple[str, ...]]] = []
        group = [ordered[0]]
        for position in ordered[1:]:
            if order[position] - order[group[-1]] == 1:
                group.append(position)
                continue
            runs.append(self._describe_run(group))
            group = [position]
        runs.append(self._describe_run(group))
        return tuple(runs)

    @staticmethod
    def _describe_run(
        group: list[Position],
    ) -> tuple[Position, Position, int, tuple[str, ...]]:
        regions: list[str] = []
        for position in group:
            region = region_of_position(position)
            if region is not None and region not in regions:
                regions.append(region)
        return (group[0], group[-1], len(group), tuple(regions))

    def counts_by_region(self) -> dict[str, tuple[int, int]]:
        """region -> (non-human peptides touching it, total peptides touching it).

        A peptide spanning two regions is counted in both, so the totals sum to
        more than `n_peptides`. That is intended: the question is "how much of
        the evidence against this region is there". They are not a partition.
        """
        counts: dict[str, list[int]] = {}
        for peptide in self.peptides:
            for region in peptide.regions:
                entry = counts.setdefault(region, [0, 0])
                entry[1] += 1
                if not peptide.found:
                    entry[0] += 1
        return {region: (bad, total) for region, (bad, total) in counts.items()}

    def flags_by_position(self) -> dict[Position, int]:
        """IMGT position -> how many non-human peptides cover it.

        The per-residue view. A residue in the middle of a non-human stretch is
        covered by up to `peptide_length` non-human windows; one at the edge by
        fewer. The count is therefore a crude measure of how deep inside
        non-human sequence a position sits.
        """
        flags: dict[Position, int] = {}
        for peptide in self.non_human_peptides:
            for position in peptide.positions:
                flags[position] = flags.get(position, 0) + 1
        return flags

    def worst_positions(self, limit: int = 10) -> list[tuple[Position, int, str | None]]:
        """The positions covered by the most non-human peptides."""
        flags = self.flags_by_position()
        ranked = sorted(flags.items(), key=lambda item: (-item[1], imgt_sort_key(item[0])))
        return [
            (position, count, region_of_position(position)) for position, count in ranked[:limit]
        ]

    def provenance(self) -> dict[str, object]:
        """What produced this number, for `run.json`."""
        from importlib.metadata import version

        return {
            "method": "OASis (9-mer repertoire prevalence)",
            "method_reference": "Prihoda et al., MAbs 2022, doi:10.1080/19420862.2021.2020203",
            "implementation": "promb",
            "promb_version": version("promb"),
            "database": self.database,
            "database_description": DATABASES.get(self.database, "?"),
            "peptide_length": self.peptide_length,
            "threshold_name": self.threshold_name,
            "threshold_min_donor_fraction": self.threshold_min_donor_fraction,
            "reference_file": "OASis_9mers_v1_10perc_subjects.txt.gz (bundled with promb)",
            "human_proteome_file": "uniprotkb_proteome_UP000005640_2025_04_21.fasta.gz",
            "oas_snapshot": (
                "not stated by promb. The bundled table is OASis 9-mers v1, which "
                "BioPhi built from an OAS snapshot accessed November 2019 "
                "(118,713,869 sequences, 231 subjects, 26 studies) -- inferred from "
                "the file name matching BioPhi's v1 release; promb asserts nothing."
            ),
            "note": (
                "promb ships only the relaxed (>=10% of donors) peptide set. "
                "Published OASis correlations with clinical ADA used the medium "
                "(>=50%) threshold, which requires BioPhi's full 22 GB database."
            ),
        }

    def warnings(self) -> tuple[str, ...]:
        """Caveats that must travel with this number."""
        messages = list(self.domain_warnings)
        if self.n_peptides == 0:
            messages.append(
                f"chain is shorter than the {self.peptide_length}-mer window; no score."
            )
        if any(m.startswith(PARTIAL_DOMAIN_PREFIX) for m in self.domain_warnings):
            messages.append(
                "humanness computed on a partial domain: the score is over fewer "
                "peptides than a complete chain and is not directly comparable."
            )
        return tuple(messages)


def _windows(domain: Domain, length: int) -> list[Peptide]:
    """Slide a window over the domain, carrying IMGT positions along.

    The residues are taken in IMGT order rather than raw input order so that a
    chain with CDR3 insertions is chopped in true sequence order.
    """
    positions = sorted(domain.numbering, key=imgt_sort_key)
    residues = [domain.numbering[p] for p in positions]
    peptides: list[Peptide] = []
    for start in range(len(residues) - length + 1):
        peptides.append(
            Peptide(
                sequence="".join(residues[start : start + length]),
                start_index=start,
                positions=tuple(positions[start : start + length]),
                found=False,
            )
        )
    return peptides


def score_domain(
    domain: Domain,
    *,
    database: str = "human-oas",
) -> HumannessScore:
    """Score one numbered V-domain against a peptide reference."""
    db = get_database(database)
    length = db.peptide_length

    peptides = tuple(
        Peptide(
            sequence=window.sequence,
            start_index=window.start_index,
            positions=window.positions,
            # A window containing an ambiguity code cannot match a reference
            # peptide. Short-circuiting keeps it out of the lookup rather than
            # relying on the reference never containing an X.
            found=(
                False if set(window.sequence) & AMBIGUITY_CODES else db.contains(window.sequence)
            ),
        )
        for window in _windows(domain, length)
    )

    is_oas = database == "human-oas"
    return HumannessScore(
        peptides=peptides,
        database=database,
        peptide_length=length,
        threshold_name=OAS_THRESHOLD_NAME if is_oas else "exact match",
        threshold_min_donor_fraction=OAS_THRESHOLD_MIN_DONOR_FRACTION if is_oas else 0.0,
        chain_type=domain.chain_type,
        domain_warnings=domain.warnings,
    )


def score_sequence(sequence: str, *, database: str = "human-oas") -> HumannessScore:
    """Number `sequence` and score its single V-domain."""
    from .number import number_single_chain

    return score_domain(number_single_chain(sequence), database=database)
