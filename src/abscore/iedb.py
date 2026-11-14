"""Cached client for the IEDB MHC class II binding prediction service.

Why the service and not the software
------------------------------------
IEDB's strongest class II predictor is NetMHCIIpan, which is also what the BioPhi
paper used (netMHCIIpan 3.1). Its *downloadable* form carries the CBS/DTU
NetMHC academic licence: free only for publicly funded academic institutions, no
redistribution, and -- the clause that matters for a public code sample --
benchmark results may not be published without written consent.

None of that applies here, because nothing is downloaded. The sequence is sent to
IEDB's own server, which exists to answer exactly this question, and the answer
comes back. No restricted binary enters this repository. The method name and
service URL are recorded with every result so a reader knows which predictor
produced which number.

Service reliability, as measured
--------------------------------
Measured against the live endpoint on 2026-11-12, with a 120 aa heavy chain:

- **Availability.** Six single requests one minute apart: five returned 106 rows
  in 1.5-1.9 s, one failed with HTTP 500 in 0.3 s. 83% first-attempt success.
- **Rate limiting: none found at this scale.** Eight requests (one chain's
  allele panel) at gaps of 10, 5, 3, 2, 1 and 0 seconds: 48/48 succeeded,
  including back-to-back. No Retry-After or rate-limit headers. Request spacing
  is not the variable.
- **Batching works.** Multi-FASTA requests of 2, 3, 5 and 10 full-length chains
  all returned complete results, every `seq_num` present, in 2-7 s.

The failure mode is a transient outage window. Throttling and payload limits are
both ruled out: during one such window earlier the same day, single requests,
batched requests and an eight-allele request all failed regardless of shape or
spacing, and all of them worked once it passed. Failures are therefore
correlated in time, which is why retries use a wide backoff.

An earlier version of this module recorded "no consistent batch size limit" and
treated batching as undependable. That was wrong: those probes all ran inside
one outage window. The measurement is recorded here rather than the inference.

Note also that a failure arrives as HTTP 200 with no prediction rows as often as
it arrives as HTTP 500, so an empty body is retried like a transport error.

Reproducibility
---------------
Every response is cached on disk, keyed by a hash of the request. A repeated
request is served from the cache without touching the network, and `offline=True`
refuses to make a request at all, so a committed result can be shown to replay
from cache rather than silently re-fetching. This matters because IEDB updates
its models: without a frozen cache, a result committed today could stop
reproducing after a server-side change, with nothing to show what moved.

References
----------
- IEDB Analysis Resource: https://www.iedb.org
- Reynisson B et al. NetMHCpan-4.1 and NetMHCIIpan-4.0. Nucleic Acids Res 2020.
- Vita R et al. The Immune Epitope Database (IEDB): 2018 update.
  Nucleic Acids Res 2019;47:D339-D343.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

#: The legacy tools endpoint. Chosen over the next-generation API because it is
#: synchronous -- one POST returns the full TSV -- where the NG API requires
#: submitting a pipeline and polling for a result. For a few hundred peptides the
#: legacy endpoint answers in seconds and needs no polling loop.
IEDB_MHCII_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhcii/"

#: Methods the endpoint accepts. The default matches what the BioPhi paper used,
#: so results are comparable to the published analysis.
#:
#: `netmhciipan_el` and `netmhciipan_ba` are the current recommended predictors
#: (IEDB labels them "recommended ... 2023.09"). `consensus`, `nn_align` and
#: `smm_align` are also DTU-derived. `comblib` and `tepitope` are matrix methods
#: of independent origin, kept available for anyone who wants a predictor with no
#: DTU lineage at all, at the cost of accuracy.
METHODS = (
    "netmhciipan_el",
    "netmhciipan_ba",
    "consensus",
    "nn_align",
    "smm_align",
    "comblib",
    "tepitope",
)
DEFAULT_METHOD = "netmhciipan_el"

#: The eight HLA-DRB1 alleles used in the BioPhi paper's MHC II analysis, which
#: in turn followed Amimeur et al. Using the same panel keeps our epitope counts
#: comparable to theirs.
#:
#: NOTE: these are the paper's eight. A nine-allele set including DRB1*04:05,
#: *08:02 and *13:02 is also widely cited; it is not what the paper used.
DEFAULT_ALLELES = (
    "HLA-DRB1*01:01",
    "HLA-DRB1*03:01",
    "HLA-DRB1*04:01",
    "HLA-DRB1*07:01",
    "HLA-DRB1*08:01",
    "HLA-DRB1*11:01",
    "HLA-DRB1*13:01",
    "HLA-DRB1*15:01",
)

#: Class II peptides are 13-25 residues with a 9-residue core in the binding
#: groove. 15-mers are the convention: long enough to contain any core, short
#: enough that the predictor can locate it.
DEFAULT_PEPTIDE_LENGTH = 15

DEFAULT_CACHE_DIR = Path("cache") / "iedb"
DEFAULT_TIMEOUT = 300
USER_AGENT = "abscore/0.1 (antibody immunogenicity triage; python-urllib)"


class IEDBError(RuntimeError):
    """Raised when IEDB cannot be reached or returns something unusable."""


class OfflineError(IEDBError):
    """Raised when a request is not cached and the network is disallowed."""


class EmptyResponseError(IEDBError):
    """Raised when IEDB answers with no prediction rows.

    This is the service's characteristic failure mode and it is transient: the
    same request can return an empty body, then the full result twenty seconds
    later, with HTTP 200 or 500 either way. It is therefore retried, unlike a
    response that has content but cannot be parsed, which indicates a real
    server-side change and is raised immediately.
    """


class MalformedResponseError(IEDBError):
    """Raised when a response has rows but one of them is unusable.

    Deliberately fatal. A dropped row surfaces as a *smaller number*, which
    reads exactly like a better antibody, so nothing announces the loss. The
    cache is written only after a successful parse, so strictness here also
    stops a truncated or changed response from being stored and replayed
    forever as though it were complete.

    Not retried. An empty body is a hiccup; content that fails to parse means
    the endpoint's contract changed, and repeating the request will not fix it.
    """


@dataclass(frozen=True)
class Binding:
    """One predicted peptide-allele binding.

    `rank` is the percentile rank: where this peptide's score falls against a
    background set of random peptides for this allele. Thresholding on rank
    rather than on predicted affinity is deliberate -- each allele's score
    distribution is different, so a single nanomolar cutoff means different things
    for different alleles, while a rank cutoff means the same thing for all of
    them.
    """

    allele: str
    peptide: str
    core: str
    start: int
    end: int
    score: float
    rank: float

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def _request_key(sequence: str, allele: str, method: str, length: int) -> str:
    payload = json.dumps(
        {
            "url": IEDB_MHCII_URL,
            "sequence": sequence,
            "allele": allele,
            "method": method,
            "length": length,
        },
        sort_keys=True,
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def _parse_tsv(text: str) -> tuple[Binding, ...]:
    """Parse IEDB's TSV response.

    Column names have changed between IEDB versions, so each field is looked up
    by any of its known names rather than by position.
    """
    if not text.strip() or not text.strip().splitlines()[1:]:
        raise EmptyResponseError(
            "IEDB returned no prediction rows. This is its usual transient failure and is retried."
        )

    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise EmptyResponseError("IEDB returned no columns")

    def pick(row: dict, *names: str) -> str | None:
        for name in names:
            if name in row and row[name] not in (None, ""):
                return row[name]
        return None

    bindings: list[Binding] = []
    for number, row in enumerate(reader, start=2):  # row 1 is the header
        if all(value in (None, "") for value in row.values()):
            continue  # a blank line; damage is handled below

        peptide = pick(row, "peptide", "sequence")
        allele = pick(row, "allele")
        rank = pick(row, "rank", "percentile_rank", "adjusted_rank")
        missing = [
            name
            for name, value in (("peptide", peptide), ("allele", allele), ("rank", rank))
            if value is None
        ]
        if missing:
            raise MalformedResponseError(
                f"IEDB response line {number} is missing {', '.join(missing)}. "
                f"Columns seen: {', '.join(reader.fieldnames or [])}. "
                "Refusing to return a partial count."
            )
        try:
            bindings.append(
                Binding(
                    allele=allele,
                    peptide=peptide,
                    core=pick(row, "core_peptide", "core", "core_sequence") or "",
                    start=int(pick(row, "start", "seq_start") or 0),
                    end=int(pick(row, "end", "seq_end") or 0),
                    score=float(pick(row, "score", "ic50", "affinity") or "nan"),
                    rank=float(rank),
                )
            )
        except ValueError as exc:
            raise MalformedResponseError(
                f"IEDB response line {number} has an unparseable value: {exc}. "
                "Refusing to return a partial count."
            ) from exc

    if not bindings:
        head = text.strip().splitlines()[:3]
        raise IEDBError(
            "IEDB returned no usable predictions. First lines of the response: " + " | ".join(head)
        )
    return tuple(bindings)


def predict_mhcii(
    sequence: str,
    *,
    alleles: tuple[str, ...] = DEFAULT_ALLELES,
    method: str = DEFAULT_METHOD,
    peptide_length: int = DEFAULT_PEPTIDE_LENGTH,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    offline: bool = False,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = 4,
    pause: float = 0.5,
    backoff: float = 3.0,
    progress: Callable[[str], None] | None = None,
) -> tuple[Binding, ...]:
    """Predict MHC class II binding for every window of `sequence`.

    Returns one `Binding` per peptide-allele pair.

    One request is made per allele, instead of one listing every allele. The
    endpoint accepts a comma-separated allele list and did serve an eight-allele
    request during development, but later returned an empty body for the same
    request shape while a single-allele request still worked. Per-allele requests
    are also better cached: adding a ninth allele re-uses the eight responses
    already on disk instead of invalidating them.

    Responses are cached by request hash; `offline=True` raises rather than
    making a request. `pause` is the delay inserted between successive requests,
    to avoid hammering a free public service.

    `progress` is called with a short status line per allele, saying whether the
    response came from cache or the network. The library stays silent unless a
    caller asks: this exists because eight sequential requests take about twenty
    seconds, and a tool that prints nothing for twenty seconds looks broken.
    """
    if method not in METHODS:
        raise IEDBError(f"unknown method {method!r}; expected one of {list(METHODS)}")
    if not alleles:
        raise IEDBError("no alleles requested")

    sequence = "".join(sequence.split()).upper()
    alleles = tuple(alleles)

    results: list[Binding] = []
    for index, allele in enumerate(alleles):
        if progress:
            cached = (
                Path(cache_dir) / f"{_request_key(sequence, allele, method, peptide_length)}.tsv"
            ).exists()
            source = "cached" if cached else "fetching"
            progress(f"allele {index + 1}/{len(alleles)} {allele} [{source}]")
        results.extend(
            _predict_one(
                sequence,
                allele,
                method=method,
                peptide_length=peptide_length,
                cache_dir=cache_dir,
                offline=offline,
                timeout=timeout,
                retries=retries,
                pause=pause if index else 0.0,
                backoff=backoff,
            )
        )
    return tuple(results)


def _predict_one(
    sequence: str,
    allele: str,
    *,
    method: str,
    peptide_length: int,
    cache_dir: str | Path,
    offline: bool,
    timeout: int,
    retries: int,
    pause: float,
    backoff: float,
) -> tuple[Binding, ...]:
    """One sequence against one allele, cached by request hash."""
    cache_path = Path(cache_dir) / f"{_request_key(sequence, allele, method, peptide_length)}.tsv"
    if cache_path.exists():
        return _parse_tsv(cache_path.read_text())

    if offline:
        raise OfflineError(
            f"no cached IEDB response for method={method}, allele={allele}, "
            f"{len(sequence)} residues, and offline=True. "
            f"Expected cache file: {cache_path}"
        )

    # The endpoint does not percent-decode its form fields: an allele sent as
    # "HLA-DRB1%2A01%3A01" yields an empty body with HTTP 500, while the same
    # allele sent literally as "HLA-DRB1*01:01" works. So "*" and ":" are left
    # unescaped. Everything else is still encoded normally.
    body = urllib.parse.urlencode(
        {
            "method": method,
            "sequence_text": sequence,
            "allele": allele,
            "length": str(peptide_length),
        },
        safe="*:",
        quote_via=urllib.parse.quote,
    ).encode()

    if pause:
        time.sleep(pause)

    last_error: Exception | None = None
    for attempt in range(retries + 1):
        request = urllib.request.Request(
            IEDB_MHCII_URL,
            data=body,
            headers={
                "User-Agent": USER_AGENT,
                "Content-Type": "application/x-www-form-urlencoded",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            bindings = _parse_tsv(text)
            break
        except (
            urllib.error.URLError,
            TimeoutError,
            OSError,
            EmptyResponseError,
        ) as exc:
            # Transport failures and empty responses are both transient here.
            # A parse failure on a response that HAS rows is not retried: it is
            # raised by _parse_tsv as IEDBError, which this does not catch.
            last_error = exc
            if attempt < retries:
                time.sleep(backoff * (2**attempt))
    else:
        raise IEDBError(
            f"IEDB request failed for allele {allele} after {retries + 1} attempts: {last_error}"
        )

    # Written only after a successful parse, so a bad response is never cached.
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(text)
    return bindings


def warm_cache(
    sequences: dict[str, str],
    *,
    alleles: tuple[str, ...] = DEFAULT_ALLELES,
    method: str = DEFAULT_METHOD,
    peptide_length: int = DEFAULT_PEPTIDE_LENGTH,
    batch_size: int = 10,
    cache_dir: str | Path = DEFAULT_CACHE_DIR,
    timeout: int = DEFAULT_TIMEOUT,
    retries: int = 4,
    pause: float = 1.0,
    backoff: float = 3.0,
) -> dict[str, int]:
    """Populate the cache for many sequences using multi-FASTA requests.

    The endpoint accepts several sequences in one request as FASTA and tags each
    row with a 1-based `seq_num`, so one request can cover a batch. That matters
    for scoring a panel: eight requests per chain means ~9,200 requests for a
    1,146-chain panel, well past a reasonable load on a free public service.
    At `batch_size=10` the same work takes ~920.

    The response is split by `seq_num` and written to the SAME per-sequence cache
    files that `predict_mhcii` reads, so this only ever warms the cache --
    nothing else changes, scoring one sequence behaves identically, and adding a
    sequence later re-uses everything already fetched.

    Returns a count of cache files written per allele. Batches that fail after
    retries are skipped and reported in the result as a shortfall rather than
    raising, so one bad batch does not lose the rest of a long run.

    Batches of 10 full-length chains were served complete in 2-7 s when the
    service was healthy (see the module docstring). Failed batches are skipped
    rather than raised, and anything this does not fetch is simply fetched one
    sequence at a time by `predict_mhcii` later, so a bad batch costs time and
    nothing else.
    """
    items = [(name, "".join(seq.split()).upper()) for name, seq in sequences.items()]
    written: dict[str, int] = dict.fromkeys(alleles, 0)

    for allele in alleles:
        pending = [
            (name, seq)
            for name, seq in items
            if not (
                Path(cache_dir) / f"{_request_key(seq, allele, method, peptide_length)}.tsv"
            ).exists()
        ]
        for start in range(0, len(pending), batch_size):
            batch = pending[start : start + batch_size]
            try:
                per_sequence = _predict_batch(
                    batch,
                    allele,
                    method=method,
                    peptide_length=peptide_length,
                    timeout=timeout,
                    retries=retries,
                    backoff=backoff,
                )
            except IEDBError:
                continue
            for index, (_, seq) in enumerate(batch, start=1):
                text = per_sequence.get(index)
                if not text:
                    continue
                path = Path(cache_dir) / f"{_request_key(seq, allele, method, peptide_length)}.tsv"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
                written[allele] += 1
            time.sleep(pause)
    return written


def _predict_batch(
    batch: list[tuple[str, str]],
    allele: str,
    *,
    method: str,
    peptide_length: int,
    timeout: int,
    retries: int,
    backoff: float,
) -> dict[int, str]:
    """One multi-FASTA request; returns seq_num -> a single-sequence TSV.

    Each returned TSV is rebuilt with the original header and that sequence's
    rows only, so it is indistinguishable from a single-sequence response and
    can be cached under the per-sequence key.
    """
    fasta = "\n".join(f">s{i}\n{seq}" for i, (_, seq) in enumerate(batch, start=1))
    # ">" and newline are left unescaped for the same reason as "*" and ":":
    # this endpoint does not percent-decode its form fields.
    body = urllib.parse.urlencode(
        {
            "method": method,
            "sequence_text": fasta,
            "allele": allele,
            "length": str(peptide_length),
        },
        safe="*:\n>",
        quote_via=urllib.parse.quote,
    ).encode()

    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            request = urllib.request.Request(
                IEDB_MHCII_URL,
                data=body,
                headers={
                    "User-Agent": USER_AGENT,
                    "Content-Type": "application/x-www-form-urlencoded",
                },
            )
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            lines = [line for line in text.strip().splitlines() if line.strip()]
            if len(lines) < 2:
                raise EmptyResponseError("IEDB returned no prediction rows for this batch")
            header, rows = lines[0], lines[1:]
            columns = header.split("\t")
            if "seq_num" not in columns:
                raise IEDBError("batched response has no seq_num column")
            index_of = columns.index("seq_num")
            grouped: dict[int, list[str]] = {}
            for row in rows:
                fields = row.split("\t")
                if len(fields) <= index_of:
                    continue
                try:
                    number = int(fields[index_of])
                except ValueError:
                    continue
                grouped.setdefault(number, []).append(row)
            if not grouped:
                raise EmptyResponseError("batched response grouped to nothing")
            return {n: "\n".join([header, *rs]) + "\n" for n, rs in grouped.items()}
        except (
            urllib.error.URLError,
            TimeoutError,
            OSError,
            EmptyResponseError,
        ) as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(backoff * (2**attempt))
    raise IEDBError(f"batched IEDB request failed for {allele}: {last_error}")


def service_provenance(
    method: str = DEFAULT_METHOD,
    *,
    peptide_length: int = DEFAULT_PEPTIDE_LENGTH,
    alleles: Sequence[str] = DEFAULT_ALLELES,
) -> dict[str, object]:
    """What produced the predictions, for `run.json`.

    Every setting is a parameter, because provenance that reports the default
    while the run used an override is worse than no provenance at all. Callers
    pass what the run actually used; `allele_panel_source` explains the default
    panel and stays correct as a note on it either way.
    """
    return {
        "service": "IEDB Analysis Resource, MHC class II binding",
        "url": IEDB_MHCII_URL,
        "method": method,
        "peptide_length": peptide_length,
        "alleles": list(alleles),
        "alleles_are_default_panel": tuple(alleles) == DEFAULT_ALLELES,
        "allele_panel_source": (
            "The eight HLA-DRB1 alleles used in Prihoda et al., MAbs 2022, "
            "following Amimeur et al. A nine-allele set including DRB1*04:05, "
            "*08:02 and *13:02 is also widely cited but is not what that paper used."
        ),
        "method_basis": (
            "netmhciipan_el is the default by choice; IEDB's own default was "
            "not simply accepted. NetMHCIIpan 4.x has two heads: BA, trained on measured "
            "binding affinities, and EL, trained on peptides recovered from "
            "MHC II by mass spectrometry. EL scores presentation, which "
            "subsumes processing and transport, and a peptide that binds well "
            "but is never processed out of the protein is not an epitope. "
            "Consequence: a rank from this method is not an affinity, and "
            "nothing here should be read in nM. netmhciipan_ba is available."
        ),
        "threshold_basis": (
            "percentile rank in place of predicted affinity: each allele's score "
            "distribution differs, so one nanomolar cutoff is not comparable "
            "across alleles while one rank cutoff is"
        ),
        "licence_note": (
            "NetMHCIIpan is called through IEDB's hosted service, never "
            "downloaded. The CBS/DTU NetMHC academic licence governs the "
            "downloadable software -- including a restriction on publishing "
            "benchmark results -- and no such software is vendored here."
        ),
        "software_version": (
            "not reported by the IEDB endpoint; responses are cached so a run "
            "can be replayed against the model that produced it"
        ),
    }
