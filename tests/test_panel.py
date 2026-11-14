"""Panel-wide tests, including the held-out check on the framework-gap rule.

These are slower than the unit tests (they number ~1148 sequences), so they are
marked `panel` and can be deselected with `-m "not panel"`.

The point of the split test
---------------------------
`STRUCTURAL_FRAMEWORK_GAPS` in `abscore.number` was derived by observing which
framework positions are empty across this panel. Asserting that the same panel
then produces no false warnings is circular: the rule was fitted to that data, so
the assertion cannot fail.

So the rule is re-derived here from a random half of the panel and checked
against the held-out half. That is a real test. What it shows is that the rule is
stable under resampling *within* this panel -- not that it generalises to
sequences from another source. See `panels/imgt_mab_db/SOURCE.md`.
"""

import random
from collections import Counter, defaultdict

import pytest

from abscore.germline import assign_closest_germline
from abscore.humanness import score_domain
from abscore.number import (
    FRAMEWORK_REGIONS,
    IMGT_REGIONS,
    STRUCTURAL_FRAMEWORK_GAPS,
    NumberingError,
    number_sequence,
)
from abscore.panel import load_ada, load_panel

pytestmark = pytest.mark.panel

LAST_POSITION = {"H": 128, "K": 127, "L": 127}
FRAMEWORK_POSITIONS = [
    position
    for region in FRAMEWORK_REGIONS
    for position in range(IMGT_REGIONS[region][0], IMGT_REGIONS[region][1] + 1)
]


@pytest.fixture(scope="module")
def domains():
    """Every panel chain that numbers, as (entry, domain) pairs."""
    result = []
    for entry in load_panel():
        try:
            parsed = number_sequence(entry.sequence)
        except NumberingError:
            continue
        if len(parsed) == 1:
            result.append((entry, parsed[0]))
    return result


def _scored_chains(domains):
    """Reuse the already-numbered `domains` fixture rather than calling
    `score_panel`, which would number all 1148 sequences a second time."""
    from abscore.validation import ScoredChain

    rows = []
    for entry, domain in domains:
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


def _derive_gap_rule(subset, min_fraction=0.10):
    """Which framework positions are empty in >= min_fraction of chains, per type."""
    totals = Counter()
    absent = defaultdict(Counter)
    for _, domain in subset:
        totals[domain.chain_type] += 1
        present = {p[0] for p in domain.numbering}
        for position in FRAMEWORK_POSITIONS:
            if position <= LAST_POSITION[domain.chain_type] and position not in present:
                absent[domain.chain_type][position] += 1
    return {
        chain: frozenset(
            position
            for position, count in absent[chain].items()
            if count / totals[chain] >= min_fraction
        )
        for chain in totals
    }


def _unexpected_framework_gaps(domain, rule):
    """Framework positions empty in this domain that the rule does not excuse."""
    present = {p[0] for p in domain.numbering}
    first = min(present)
    allowed = rule.get(domain.chain_type, frozenset())
    return [
        position
        for position in FRAMEWORK_POSITIONS
        if first <= position <= LAST_POSITION[domain.chain_type]
        and position not in present
        and position not in allowed
    ]


class TestPanelIntegrity:
    def test_panel_size(self):
        assert len(load_panel()) == 1148

    def test_almost_every_sequence_numbers(self, domains):
        assert len(domains) >= 1146

    def test_chain_order_in_the_fasta_is_not_consistent(self, domains):
        """Pins why chain type must come from numbering; the suffix lies."""
        by_suffix = defaultdict(Counter)
        for entry, domain in domains:
            by_suffix[entry.record_id.split("_")[-1]][domain.chain_type] += 1
        assert by_suffix["seq1"]["H"] > 0
        assert by_suffix["seq1"]["K"] + by_suffix["seq1"]["L"] > 0

    def test_ada_panel_loads(self):
        """217, matching the BioPhi paper exactly.

        The paper states "only 217/218 therapeutics were included" (catumaxomab's
        sequence was unavailable). Note that `wc -l` reports 228 lines for this
        file: some Reference fields contain embedded newlines, so the line count
        is not the record count.
        """
        ada = load_ada()
        assert len(ada) == 217
        assert all(0.0 <= v <= 100.0 for v in ada.values())


class TestFrameworkGapRuleHeldOut:
    """The rule must hold on data it was not derived from."""

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_rule_derived_on_half_the_panel_holds_on_the_other_half(self, domains, seed):
        indices = list(range(len(domains)))
        random.Random(seed).shuffle(indices)
        midpoint = len(indices) // 2
        train = [domains[i] for i in indices[:midpoint]]
        held_out = [domains[i] for i in indices[midpoint:]]

        rule = _derive_gap_rule(train)
        offenders = [
            (entry.record_id, _unexpected_framework_gaps(domain, rule))
            for entry, domain in held_out
            if _unexpected_framework_gaps(domain, rule)
        ]
        # Some chains do carry real framework deletions, so a small number of
        # flags on held-out data is correct. What would be wrong is the rule
        # differing between halves, tested below.
        assert len(offenders) / len(held_out) < 0.02, offenders[:5]

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_the_derived_rule_is_the_one_shipped_in_the_code(self, domains, seed):
        indices = list(range(len(domains)))
        random.Random(seed).shuffle(indices)
        rule = _derive_gap_rule([domains[i] for i in indices[: len(indices) // 2]])
        for chain in ("H", "K", "L"):
            assert rule[chain] == STRUCTURAL_FRAMEWORK_GAPS[chain], chain

    def test_the_rule_is_identical_on_both_halves(self, domains):
        """Stability under resampling: the rule is a property of the chain type,
        not of which antibodies happened to be sampled."""
        indices = list(range(len(domains)))
        random.Random(7).shuffle(indices)
        midpoint = len(indices) // 2
        first = _derive_gap_rule([domains[i] for i in indices[:midpoint]])
        second = _derive_gap_rule([domains[i] for i in indices[midpoint:]])
        assert first == second


class TestOriginSeparation:
    """Validation: do the scores track the known engineering history?

    Not circular -- nothing is fitted to these labels, and the labels are human
    curation of engineering history rather than anything derived from sequence.
    But it is a single panel, so it corroborates and nothing more.
    """

    @pytest.fixture(scope="class")
    @classmethod
    def scored(cls, domains):
        rows = []
        for entry, domain in domains:
            if not entry.is_single_origin and not entry.is_veterinary:
                continue
            assignment = assign_closest_germline(domain)
            if assignment.identity is None:
                continue
            rows.append(
                {
                    "origin": entry.origin,
                    "v": assignment.identity,
                    "fr": assignment.framework_identity,
                    "oasis": score_domain(domain).score,
                }
            )
        return rows

    @staticmethod
    def _auc(positives, negatives):
        wins = sum(1.0 if p > n else (0.5 if p == n else 0.0) for p in positives for n in negatives)
        return wins / (len(positives) * len(negatives))

    def _split(self, scored, metric, positive, negative):
        return (
            [r[metric] for r in scored if r["origin"] in positive],
            [r[metric] for r in scored if r["origin"] in negative],
        )

    @pytest.mark.parametrize("metric", ["v", "fr", "oasis"])
    def test_every_metric_separates_human_from_murine(self, scored, metric):
        positives, negatives = self._split(scored, metric, {"Human"}, {"Mouse", "Chimeric"})
        assert self._auc(positives, negatives) > 0.95

    def test_plain_germline_identity_is_not_beaten_by_oasis(self, scored):
        """The honest-baseline finding, as a test rather than a claim.

        If a future change makes OASis the better separator, this test fails and
        the README's central claim has to be rewritten -- which is the point.
        """
        germline = self._auc(*self._split(scored, "v", {"Human"}, {"Mouse", "Chimeric"}))
        oasis = self._auc(*self._split(scored, "oasis", {"Human"}, {"Mouse", "Chimeric"}))
        assert germline >= oasis

    @pytest.mark.parametrize("chain_type", ["H", "K"])
    def test_the_baseline_finding_holds_within_each_chain_type(self, domains, chain_type):
        """The pooled comparison above is dominated by heavy chains, so it can
        hold while being false for one chain type. On kappa chains the margin
        is about 4e-4 -- a real but tiny win for germline identity -- and a
        pooled test would not notice it reversing.

        Lambda is not parametrised because the panel holds no chimeric or
        murine lambda chains at all, so there is nothing to separate. That
        absence is reported rather than papered over; see
        `abscore.validation.separation_table`.
        """
        from abscore.validation import separation_table

        rows = _scored_chains(domains)
        by_group = {row.group: row for row in separation_table(rows)}
        row = by_group[chain_type]
        assert row.positives >= 50 and row.negatives >= 50, (row.positives, row.negatives)
        assert row.aucs["v_identity"] >= row.aucs["oasis"], row.aucs

    def test_the_panel_cannot_say_anything_about_lambda_chains(self, domains):
        """A limit worth pinning: every lambda chain in the panel is human or
        humanized, so no lambda separation is measurable. If a future panel
        adds murine lambda chains this fails, and the figure's empty cells and
        the README's caveat both need revisiting."""
        from abscore.validation import separation_table

        by_group = {row.group: row for row in separation_table(_scored_chains(domains))}
        assert by_group["L"].negatives == 0
        assert by_group["L"].aucs == {}

    def test_framework_identity_separates_humanized_from_chimeric_best(self, scored):
        """A humanized antibody has a human framework; a chimeric one does not."""
        task = ({"Humanized"}, {"Mouse", "Chimeric"})
        framework = self._auc(*self._split(scored, "fr", *task))
        oasis = self._auc(*self._split(scored, "oasis", *task))
        assert framework > oasis

    def test_veterinary_antibodies_score_below_human(self, scored):
        """Veterinary classes sit below human, and that is all this shows.

        An earlier version of this test asserted that caninized antibodies are
        the *least* human class of all, with a 0.10 slack. That claim is false
        on this panel -- caninized median V identity is 0.751 against 0.701 for
        murine -- and the slack was wide enough to hide it. Which is the point:
        a passing assertion with a tolerance bolted on is not evidence.

        So this asserts only the defensible direction, and the classes are
        small: n=12 caninized, n=6 felinized. Nothing here supports a claim
        about how veterinary and murine antibodies rank against each other.
        """
        import statistics

        human = [r["v"] for r in scored if r["origin"] == "Human"]
        veterinary = [r["v"] for r in scored if r["origin"] in {"Caninized", "Felinized"}]
        assert len(veterinary) >= 10, "too few to say anything at all"
        assert statistics.median(veterinary) < statistics.median(human)
