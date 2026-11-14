"""Run provenance: what produced every number, written to `run.json`.

The brief this project answers asks for every number to be traceable to its
source, with tool and version. That is only true if it is written down, so each
run emits a `run.json` next to its output recording the environment, the version
of every tool and reference dataset used, and the settings in force.

Two things here are deliberately awkward to read, and both are honest:

- `germline.imgt_release` is `"not recorded by ANARCI"`. ANARCI ships no IMGT
  release number with its germline set, so the set is pinned by content hash
  instead. A hash identifies the data exactly even without a version label, but
  it cannot tell you whether that data is current.
- `known_sensitivities` records measured consequences of the choices above,
  rather than asserting they do not matter.

Every entry below carries `reproducible_from_repo`. Some of these measurements
were made during development against inputs this repository does not ship --
the OAS data units behind the numbering comparison are tens of megabytes and
nothing in this checkout recomputes from them, so they were read locally and
are not redistributed. Those numbers are true and dated, and a reader cannot rerun
them from this checkout. Saying which is which is the point: an unreproducible
number presented as a repository result would be the one thing this file exists
to prevent.
"""

from __future__ import annotations

import json
import platform
import re
import subprocess
import sys
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .germline import SCOPES, germline_provenance
from .humanness import DATABASES, OAS_THRESHOLD_MIN_DONOR_FRACTION, OAS_THRESHOLD_NAME
from .number import (
    ANCHOR_POSITIONS,
    IMGT_REGIONS,
    MAX_AMBIGUOUS_FRACTION,
    STRUCTURAL_FRAMEWORK_GAPS,
)

#: Where each load-bearing convention was verified, and when. A number whose
#: definition came from recall rather than a source does not belong in a report.
VERIFIED_AGAINST = {
    "imgt_region_boundaries": {
        "source": "https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html",
        "retrieved": "2026-11-12",
        "note": "FR1 1-26, CDR1 27-38, FR2 39-55, CDR2 56-65, FR3 66-104, "
        "CDR3 105-117, FR4 118-128. IMGT's own pages give FR4 as 118-128 in the "
        "numbering overview and 118-129 in one amino-acid table; 128 is used here.",
        "reproducible_from_repo": "N/A -- a cited definition rather than a measurement. "
        "The boundaries as implemented are asserted in tests/test_number.py.",
    },
    "imgt_conserved_anchors": {
        "source": "https://imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html",
        "retrieved": "2026-11-12",
        "note": "1st-CYS 23, CONSERVED-TRP 41, hydrophobic 89, 2nd-CYS 104, J-PHE or J-TRP 118.",
        "reproducible_from_repo": "N/A -- a cited definition. The anchors as "
        "implemented are asserted in tests/test_number.py.",
    },
    "oasis_thresholds": {
        "source": "Prihoda et al., MAbs 2022, doi:10.1080/19420862.2021.2020203",
        "retrieved": "2026-11-12",
        "note": "loose >=1%, relaxed >=10%, medium >=50%, strict >=90% of OAS "
        "subjects. promb bundles the relaxed set only.",
        "reproducible_from_repo": "N/A -- a published definition. Which "
        "threshold this tool actually uses is recorded in settings and asserted "
        "in tests/test_humanness.py.",
    },
    "structural_framework_gaps": {
        "source": "derived from panels/imgt_mab_db (1146 complete V-domains)",
        "retrieved": "2026-11-12",
        "note": "Derived on a random half of the panel and verified on the "
        "held-out half; see tests/test_panel.py. Stable under resampling within "
        "this panel; not shown to generalise to other populations.",
        "reproducible_from_repo": "Yes. The panel is committed under "
        "panels/imgt_mab_db and tests/test_panel.py rederives the rule and "
        "checks it on the held-out half. Note the halves may split the two "
        "chains of one antibody, so they are not fully independent.",
    },
    "numbering_vs_oas": {
        "source": "OAS data units SRR8283601_1_Heavy_Bulk (Briney et al., Nature "
        "2019) and 1279049_1_Light_Bulk (Jaffe et al., 2022)",
        "retrieved": "2026-11-12",
        "note": "15,506 sequences numbered independently by OAS's own ANARCI run: "
        "99.908% position-level agreement (1,656,469/1,658,002), 96.098% of "
        "sequences identical end to end, zero region-boundary mismatches.",
        "reproducible_from_repo": "NO. OAS data is CC BY-4.0, so these two data "
        "units could be redistributed with attribution, but they are tens of "
        "megabytes and nothing here recomputes from them, so they were read "
        "locally and are not shipped. The figures are a dated development measurement, "
        "not a result this checkout recomputes. Committed instead: the numbering "
        "itself, its region-boundary tests, and the insertion-ordering cases the "
        "check turned up.",
    },
}

#: Measured consequences of known limitations. Each entry gives a number.
KNOWN_SENSITIVITIES = {
    "germline_set_choice": {
        "measured": "7.67% of 1148 panel chains (88) get a different germline "
        "GENE call from ANARCI's set than from the IMGT snapshot vendored in "
        "Merck/BioPhi-2021-publication (249 vs 403 human IGHV alleles; neither "
        "is a superset of the other).",
        "mitigated": "80 of those 88 are already reported as an ambiguous gene "
        "call; 5 more are duplicated-locus paralogs (e.g. IGHV3-23 vs "
        "IGHV3-23D) with the same sequence under a different name.",
        "residual": "3 chains (0.26%) change gene without being flagged, all "
        "cases where one set holds an allele the other lacks.",
        "reproducible_from_repo": "Partly. The panel is committed; the IMGT "
        "snapshot compared against is not vendored here, so the 7.67% figure is "
        "a dated development measurement. What the suite does check is that an "
        "ambiguous call is reported whenever genes are indistinguishable.",
    },
    "domain_extraction_boundary": {
        "measured": "Scores are computed over the V-domain the numbering HMM "
        "extracts, not over the submitted sequence. Where that boundary falls "
        "can move by a residue when the flanking sequence changes: appending "
        "GGGGSHHHHHH to one lambda fixture pulls the first added residue into "
        "IMGT 128, lengthening the extracted domain from 110 to 111 residues "
        "and moving OASis from 83.3% to 82.5%. The same suffixes on a heavy "
        "fixture changed neither the extracted domain nor the scores.",
        "mitigated": "Every report states the extracted length against the "
        "input length, and the JSON carries input_residues, "
        "extracted_residues and query_span, so a score is attributable to a "
        "specific stretch of the input.",
        "residual": "Scores from differently-flanked copies of one domain are "
        "not guaranteed identical. Both numbers are correct over the sequence "
        "each was computed on; the sequences differ.",
        "reproducible_from_repo": "Yes -- tests/test_number.py checks that the "
        "extracted length is reported whenever it is shorter than the input, "
        "and the lambda boundary shift is pinned as a regression test.",
    },
    "low_complexity_inflates_humanness": {
        "measured": "OASis counts 9-mer windows, so a long repeat of one "
        "residue contributes many identical windows. Inserting 100 glycines "
        "into a heavy-chain CDR3 raised OASis from 70.5% to 81.1% (120 to 220 "
        "residues), because GGGGGGGGG occurs in the reference peptide set and "
        "every repeated window matches it. Germline V identity was unchanged, "
        "the insertion being in CDR3.",
        "mitigated": "Nothing in the scoring is changed to suppress this. The "
        "score is the published OASis definition and altering it silently "
        "would make the number incomparable to the literature it comes from.",
        "residual": "A high OASis does not establish a realistic loop length, "
        "a plausible structure, or low immunogenicity. Synthetic or "
        "low-complexity domains are outside the range this tool was checked "
        "on, which is therapeutic antibodies of ordinary length.",
        "reproducible_from_repo": "Yes -- the arithmetic is asserted directly "
        "in tests/test_humanness.py from the committed reference set.",
    },
    "oasis_threshold": {
        "measured": "Scores use the relaxed threshold (>=10% of OAS subjects), "
        "because that is the peptide set promb bundles.",
        "residual": "The published OASis correlation with clinical ADA "
        "(R = -0.53) used the medium threshold (>=50%), which needs BioPhi's "
        "full 22 GB database. Scores here are not comparable to published "
        "medium-threshold numbers.",
        "reproducible_from_repo": "Yes -- the threshold in force is recorded in "
        "settings and asserted by the suite. The published comparison figure is "
        "cited from the paper; it is not recomputed here.",
    },
    "oas_reference_quality": {
        "measured": "OAS contains records with compensated indels that leave a "
        "window of the amino acid sequence out of frame relative to germline, "
        "while OAS's own flags report them as productive and not frameshifted.",
        "residual": "Some peptides counted as human in the reference set are "
        "sequencing or annotation artefacts. Not quantified.",
        "reproducible_from_repo": "No -- an observation from reading OAS "
        "records locally, which are not redistributed here. Stated because it "
        "bears on what the humanness score means. It is an observation, not a result "
        "this repository produces.",
    },
    "reference_distance": {
        "measured": "Distance to the nearest reference peptide varies sharply "
        "between chains and is not interchangeable with the count. A human "
        "lambda chain at 95.4% germline identity, carrying four somatic "
        "mutations, gives 28 unmatched hits whose distinct 9-mer cores are "
        "every one exactly one residue from a germline peptide. A murine "
        "control gives 60 unmatched hits, none of them within one residue.",
        "mitigated": "The headline count is every core with no exact match "
        "(threshold 1), so neither chain is silently discounted, and the "
        "distance is reported as an annotation (0, 1, 2+). The count excluding "
        "one-substitution variants is reported separately as *_distant.",
        "residual": "Whether a one-substitution core is actually lower risk is "
        "NOT established here: it depends on which position changed and whether "
        "a T-cell receptor contacts it, neither of which is modelled. The "
        "*_distant count assumes it does and is therefore exploratory only. An "
        "earlier version made that assumption the default headline, which "
        "reported 0 for the lambda chain above.",
        "reproducible_from_repo": "Yes -- both chains replay from the committed "
        "responses under tests/data/iedb and the counts are asserted there.",
    },
    "junctional_window": {
        "measured": "The germline definition of self cannot classify the "
        "junctional window (IMGT 108-114 heavy, 109-115 light): a germline D "
        "gene has no single protein sequence, being used in any of six reading "
        "frames and trimmed by a variable number of nucleotides at each end.",
        "mitigated": "Peptides whose binding core overlaps that window are "
        "reported as unclassifiable by germline rather than as unmatched. An "
        "earlier version counted them as unmatched, which overstated the burden "
        "on adalimumab by 17 of 35. Unclassifiable records an abstention and gets "
        "its own category, so that matched "
        "+ unmatched + unclassifiable reconciles to the predicted hits.",
        "residual": "The repertoire definition does apply there, being built "
        "from real rearranged sequences, so CDR3 risk is assessed under that "
        "definition rather than left unassessed. A low germline unmatched count "
        "still leaves the junctional part of CDR3 open.",
        "reproducible_from_repo": "Yes -- the window, the abstention and the "
        "adalimumab figure all replay from committed responses.",
    },
    "epitope_panel_validation": {
        "measured": "Eight requests per chain means roughly 9,200 requests to "
        "validate epitope burden across the 1,146-chain panel. Multi-FASTA "
        "batching reduces that to about 920: batches of 10 full-length chains "
        "returned complete results in 2-7 s (measured 2026-11-12).",
        "mitigated": "Empty responses are retried with backoff, and batched "
        "cache warming skips failed batches rather than aborting a run.",
        "residual": "Epitope burden is demonstrated on individual antibodies and "
        "has NOT been validated for origin separation across the full panel, "
        "unlike the humanness and germline scores. Any epitope count is reported "
        "without a reference distribution, so a single count has nothing to be "
        "compared against.",
        "reproducible_from_repo": "The committed responses under tests/data/iedb "
        "replay four chains offline, which demonstrates integration and the "
        "counting rules -- not predictive validity, and not panel-scale "
        "behaviour. The request-timing figures are dated live measurements.",
    },
    "iedb_service_reliability": {
        "measured": "Measured 2026-11-12 against the live endpoint with a 120 aa "
        "chain: 83% first-attempt success (5/6 one minute apart, the failure an "
        "HTTP 500 in 0.3 s); no rate limiting found at this scale (48/48 across "
        "gaps of 10, 5, 3, 2, 1 and 0 seconds, back-to-back included, with no "
        "rate-limit headers); successful responses were byte-identical every "
        "time. The failure mode is a transient outage window in which requests "
        "fail regardless of shape or spacing, so failures are correlated in time.",
        "mitigated": "Up to four retries with a wide exponential backoff, since "
        "correlated failures defeat a tight one. An empty 200 is retried like a "
        "transport error, because that is how a failure often presents. Every "
        "successful response is cached, so a result need never be re-fetched.",
        "residual": "Scoring an arbitrary sequence depends on a free public "
        "service being up at that moment. Committed results replay from cache "
        "and do not depend on it. The service reports no model version, so the "
        "retrieval date is the only handle on which model produced a number.",
        "reproducible_from_repo": "No -- these are dated measurements against "
        "the live endpoint and would have to be remeasured. Offline replay from "
        "the committed responses is reproducible and is what the suite uses.",
    },
    "ambiguity_codes": {
        "measured": f"Sequences are accepted up to {MAX_AMBIGUOUS_FRACTION:.0%} "
        "ambiguity codes (X/B/Z). Rejecting any sequence containing one made "
        "3,469 of 17,786 real OAS records (19.5%) unscorable.",
        "residual": "An ambiguity code never counts as a germline match and "
        "never occurs in a reference peptide, so both scores are depressed "
        "slightly on sequences containing one. Reported as a warning.",
        "reproducible_from_repo": "Partly. The tolerance and its consequences "
        "are asserted by the suite; the 19.5%-of-OAS figure that motivated the "
        "tolerance is a dated development measurement on data not shipped here.",
    },
}


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def hmmer_version() -> str | None:
    """HMMER's version, as reported by the binary ANARCI shells out to.

    ANARCI's numbering depends on it, so it belongs in the record even though it
    ships outside Python packaging.
    """
    try:
        result = subprocess.run(
            ["hmmscan", "-h"], capture_output=True, text=True, timeout=20, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return None
    match = re.search(r"HMMER\s+([0-9][^;\s]*)", result.stdout)
    return match.group(1) if match else None


def file_digest(path: Path) -> str | None:
    """SHA-256 of a file, or None if it is not there."""
    import hashlib

    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def panel_digests() -> dict[str, str | None]:
    """Hashes of the committed reference panel, so a run can be tied to it."""
    from .panel import ADA_PATH, FASTA_PATH, LABELS_PATH

    return {path.name: file_digest(path) for path in (FASTA_PATH, LABELS_PATH, ADA_PATH)}


def tool_versions() -> dict[str, str | None]:
    return {
        "abscore": _package_version("abscore"),
        "anarci": _package_version("anarci"),
        "promb": _package_version("promb"),
        "biopython": _package_version("biopython"),
        "hmmer": hmmer_version(),
        "python": sys.version.split()[0],
    }


def settings() -> dict[str, object]:
    """The thresholds and conventions in force, written out rather than implied."""
    return {
        "numbering_scheme": "imgt",
        "imgt_regions": {name: list(bounds) for name, bounds in IMGT_REGIONS.items()},
        "conserved_anchors": {
            str(position): list(residues) for position, residues in ANCHOR_POSITIONS.items()
        },
        "structural_framework_gaps": {
            chain: sorted(positions) for chain, positions in STRUCTURAL_FRAMEWORK_GAPS.items()
        },
        "max_ambiguous_fraction": MAX_AMBIGUOUS_FRACTION,
        "germline_scopes": {name: list(regions) for name, regions in SCOPES.items()},
        "germline_identity_denominator": (
            "positions present in both query and germline; coverage reported "
            "separately, and identity over the full germline span reported as a "
            "secondary field"
        ),
        "humanness_databases": DATABASES,
        "humanness_threshold_name": OAS_THRESHOLD_NAME,
        "humanness_threshold_min_donor_fraction": OAS_THRESHOLD_MIN_DONOR_FRACTION,
    }


def run_metadata(**extra: object) -> dict[str, object]:
    """Everything needed to reproduce and audit a run.

    `extra` is merged at the top level, for a caller to record its own inputs
    (which sequences, which options).
    """
    record: dict[str, object] = {
        "run_utc": datetime.now(UTC).isoformat(timespec="seconds"),
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python_implementation": platform.python_implementation(),
        },
        "tools": tool_versions(),
        "germline": germline_provenance(),
        "settings": settings(),
        "reference_panel": panel_digests(),
        "verified_against": VERIFIED_AGAINST,
        "known_sensitivities": KNOWN_SENSITIVITIES,
    }
    record.update(extra)
    return record


def write_run_json(directory: str | Path, **extra: object) -> Path:
    """Write `run.json` into `directory`, creating it if needed."""
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    path = target / "run.json"
    path.write_text(json.dumps(run_metadata(**extra), indent=2, sort_keys=False) + "\n")
    return path
