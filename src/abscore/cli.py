"""Command line interface.

Design notes, since both are deliberate:

Progress goes to **stderr**, results to **stdout**. So `abscore score --seq ... >
report.md` writes a clean file while the status lines still appear on the
terminal. Eight sequential IEDB requests take about twenty seconds, and a tool
that prints nothing for twenty seconds looks broken.

Exit codes distinguish *bad input* from *a service being down*. A sequence that
is not an antibody is the caller's problem (exit 2); IEDB being unavailable is
not, and does not fail the run at all -- the report is produced with that section
marked unavailable, because every other number is computed offline.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .epitopes import SELF_DEFINITIONS
from .iedb import DEFAULT_ALLELES, DEFAULT_METHOD, METHODS
from .number import HmmerNotFoundError, NumberingError
from .provenance import run_metadata, write_run_json
from .report import ChainReport, analyse, render_markdown, render_summary, render_text

EXIT_OK = 0
EXIT_BAD_INPUT = 2
#: A dependency the Python lockfile cannot install is missing. Separated from
#: EXIT_BAD_INPUT because nothing is wrong with what the caller passed, and
#: from an unhandled crash because the fix is one install command.
EXIT_MISSING_DEPENDENCY = 3


def _log(message: str) -> None:
    print(f"abscore: {message}", file=sys.stderr, flush=True)


def _read_sequences(args: argparse.Namespace) -> list[tuple[str, str]]:
    """Collect (name, sequence) pairs from --seq, --fasta or stdin."""
    if args.seq:
        return [("query", args.seq)]

    if args.fasta:
        text = Path(args.fasta).read_text()
    elif not sys.stdin.isatty():
        text = sys.stdin.read()
    else:
        raise SystemExit("no input: pass --seq SEQUENCE, --fasta FILE, or pipe FASTA on stdin")

    if ">" not in text:
        return [("query", text)]

    sequences: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                sequences.append((name, "".join(parts)))
            name, parts = line[1:].strip() or f"seq{len(sequences) + 1}", []
        else:
            parts.append(line)
    if name is not None:
        sequences.append((name, "".join(parts)))
    return sequences


def _analyse_one(name: str, sequence: str, args: argparse.Namespace) -> ChainReport | None:
    try:
        return analyse(
            sequence,
            with_epitopes=not args.no_epitopes,
            cache_dir=args.cache_dir,
            offline=args.offline,
            progress=None if args.quiet else _log,
            alleles=tuple(args.alleles) if args.alleles else DEFAULT_ALLELES,
            method=args.method,
        )
    except NumberingError as exc:
        print(f"abscore: {name}: {exc}", file=sys.stderr)
        return None


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def cmd_score(args: argparse.Namespace) -> int:
    sequences = _read_sequences(args)
    if not sequences:
        print("abscore: no sequences found in input", file=sys.stderr)
        return EXIT_BAD_INPUT

    reports: list[tuple[str, ChainReport]] = []
    failures = 0
    for name, sequence in sequences:
        if len(sequences) > 1 and not args.quiet:
            _log(f"--- {name}")
        report = _analyse_one(name, sequence, args)
        if report is None:
            failures += 1
            continue
        reports.append((name, report))

    if not reports:
        return EXIT_BAD_INPUT

    out_dir = Path(args.out) if args.out else None

    # Names are shown whenever the INPUT held more than one record, regardless of how
    # more than one succeeded: a FASTA where one record fails must not silently
    # drop the names from the rest, or the caller cannot match rows to records.
    named = len(sequences) > 1
    if args.format == "summary":
        for name, report in reports:
            prefix = f"{name}\t" if named else ""
            print(f"{prefix}{render_summary(report)}")
    elif args.format == "text":
        for index, (name, report) in enumerate(reports):
            if index:
                print()
            print(
                render_text(
                    report,
                    title=None if name == "query" else name,
                    self_definition=args.self_definition,
                )
            )
    elif args.format == "markdown":
        for name, report in reports:
            text = render_markdown(
                report,
                title=None if name == "query" else name,
                self_definition=args.self_definition,
            )
            if out_dir:
                path = out_dir / f"{_slug(name)}.md"
                _write(path, text)
                _log(f"wrote {path}")
            else:
                print(text)
                if named:
                    print("\n---\n")
    elif args.format == "json":
        payload = [
            {"name": name, **_as_dict(report, args.self_definition)} for name, report in reports
        ]
        text = json.dumps(payload if named else payload[0], indent=2)
        if out_dir:
            path = out_dir / "scores.json"
            _write(path, text)
            _log(f"wrote {path}")
        else:
            print(text)

    # --save writes the full report alongside whatever went to the terminal, so
    # a text run on screen and the file never disagree. The format follows the
    # extension: .html renders tables and the position track in a browser with
    # no assets, which markdown structurally cannot do; .md stays for a README
    # or a PR comment.
    if args.save:
        path = Path(args.save)
        if path.suffix.lower() in (".html", ".htm"):
            from .html import render_html_document

            text = render_html_document(reports, self_definition=args.self_definition)
        else:
            text = "\n\n---\n\n".join(
                render_markdown(
                    report,
                    title=None if name == "query" else name,
                    self_definition=args.self_definition,
                )
                for name, report in reports
            )
        _write(path, text)
        _log(f"wrote {path}")

    if out_dir:
        path = write_run_json(
            out_dir,
            inputs={
                "sequences": len(reports),
                "failed": failures,
                "self_definition": args.self_definition,
                "epitopes_requested": not args.no_epitopes,
                "offline": args.offline,
            },
        )
        _log(f"wrote {path}")

    return EXIT_OK if not failures else EXIT_BAD_INPUT


def _slug(name: str) -> str:
    keep = [c if c.isalnum() or c in "-_" else "_" for c in name]
    return "".join(keep).strip("_") or "query"


def _as_dict(report: ChainReport, self_definition: str) -> dict:
    """Flat, machine-readable scores. Every number keeps its basis alongside."""
    germline = report.germline
    scopes = {
        name: {
            "identity": scope.identity,
            "identity_over_germline": scope.identity_over_germline,
            "coverage": scope.coverage,
            "matches": scope.matches,
            "compared": scope.compared,
            "germline_positions": scope.covered,
        }
        for name, scope in germline.scopes.items()
    }
    payload: dict = {
        "chain_type": report.chain_type,
        "imgt_span": list(report.domain.imgt_span),
        # Every score below is computed over `extracted_residues`, not over the
        # input. These three make that attributable when they differ.
        "input_residues": report.domain.input_length,
        "extracted_residues": len(report.domain.seq),
        "query_span": [report.domain.query_start, report.domain.query_end],
        "regions": dict(report.domain.regions),
        "numbering_hmm_species": report.domain.species,
        "germline": {
            "gene": germline.gene,
            "ambiguous": germline.is_ambiguous,
            "indistinguishable_from": [c.best_allele for c in germline.ambiguous_genes],
            "scopes": scopes,
        },
        "humanness": {
            "oasis_identity": report.humanness.score,
            "peptides_found": report.humanness.n_found,
            "peptides_total": report.humanness.n_peptides,
            "threshold": report.humanness.threshold_name,
            "min_donor_fraction": report.humanness.threshold_min_donor_fraction,
        },
        "warnings": list(report.warnings()),
        "provenance": report.provenance(),
    }
    if report.epitopes:
        burden = report.epitopes
        payload["epitopes"] = {
            "counts": burden.counts_by_definition(),
            "reference_reported": self_definition,
            "unmatched_threshold": burden.mismatch_threshold,
            "predictor": burden.method,
            "alleles": list(burden.alleles),
            "peptide_length": burden.peptide_length,
            "binder_percentile_rank_cutoff": burden.binder_rank,
            "by_allele": burden.counts_by_allele(self_definition),
            "by_region": burden.counts_by_region(self_definition),
            "distance_annotations": [
                {
                    "allele": e.allele,
                    "peptide": e.peptide,
                    "core": e.core,
                    "rank": e.rank,
                    "regions": list(e.regions),
                    "distance": e.distance_label(self_definition),
                }
                for e in burden.unmatched(self_definition)
            ],
        }
    else:
        payload["epitopes"] = {"status": "unavailable", "reason": report.epitope_error}
    return payload


def cmd_validate(args: argparse.Namespace) -> int:
    """Score the labelled panel, print the comparison, optionally draw it."""
    from .validation import (
        highest_scoring,
        overlapping_chains,
        render_separation_html,
        render_separation_text,
        score_panel,
        separation_table,
    )

    if not args.quiet:
        _log("scoring the committed panel; this numbers every chain and takes a minute")
    rows = score_panel(progress=None if args.quiet else _log)
    table = separation_table(rows)

    print(f"{len(rows)} panel chains scored\n")
    print(render_separation_text(table))

    print("\nWhere the groups overlap\n" + "-" * 24)
    murine = highest_scoring(rows, origins=frozenset({"Mouse"}))
    if murine is not None:
        print(
            f"Most human-like murine chain: {murine.name} ({murine.record_id}, "
            f"{murine.chain_type}) at V={murine.v_identity:.3f}."
        )
    overlap = overlapping_chains(rows)
    print(
        f"{len(overlap)} human-origin chains score at or below it. They are kept in "
        "the figure and in this count:"
    )
    for chain in overlap:
        print(
            f"  {chain.name:<24}{chain.record_id:<14}"
            f"{chain.chain_type}  V={chain.v_identity:.3f}  "
            f"FR={chain.framework_identity:.3f}  OASis={chain.oasis:.3f}"
        )
    chimeric = highest_scoring(rows, origins=frozenset({"Chimeric"}))
    if chimeric is not None:
        print(
            f"\nSeparately, the most human-like chimeric chain is {chimeric.name} "
            f"({chimeric.record_id}, {chimeric.chain_type}) at V="
            f"{chimeric.v_identity:.3f}, close to the human median. One chain, and it "
            "is most of the distance between the pooled AUC above and 1.0."
        )

    if args.figure:
        from .figure import render_figure_document, render_figure_svg

        path = Path(args.figure)
        if path.suffix.lower() == ".svg":
            path.write_text(render_figure_svg(rows))
        else:
            path.write_text(render_figure_document(rows, table_html=render_separation_html(table)))
        _log(f"wrote {path} ({path.stat().st_size:,} bytes)")

    return EXIT_OK


def cmd_provenance(args: argparse.Namespace) -> int:
    print(json.dumps(run_metadata(), indent=2))
    return EXIT_OK


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="abscore",
        description=(
            "Immunogenicity triage for an antibody heavy or light chain: "
            "germline identity, repertoire humanness, MHC class II predicted "
            "binders unmatched in a human reference."
        ),
        epilog=(
            "Scores are computed from sequence alone. Use them as triage ahead of "
            "an assay. Every report carries its "
            "own limits section."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    score = sub.add_parser("score", help="score one or more chains")
    source = score.add_argument_group("input (one of)")
    source.add_argument("--seq", help="an amino acid sequence")
    source.add_argument("--fasta", help="a FASTA file; one report per record")
    score.add_argument(
        "--format",
        choices=("summary", "text", "markdown", "json"),
        default="text",
        help=(
            "text is the full report laid out for a terminal (default); "
            "summary is one line per chain, for scripting; markdown is the same "
            "report with tables, for a file; json is machine-readable"
        ),
    )
    score.add_argument(
        "--save",
        metavar="PATH",
        help=(
            "also write the full report here. The format follows the extension: "
            ".html renders the tables and the position track in a browser, .md "
            "gives markdown. Readable on screen, complete in the file."
        ),
    )
    score.add_argument(
        "--out",
        help="directory to write reports and run.json into, instead of stdout",
    )
    score.add_argument(
        "--self-definition",
        choices=SELF_DEFINITIONS,
        default="germline",
        help="which human reference the epitope cores are compared against (default: germline)",
    )
    score.add_argument(
        "--method",
        choices=METHODS,
        default=DEFAULT_METHOD,
        help=f"IEDB class II predictor (default: {DEFAULT_METHOD})",
    )
    score.add_argument(
        "--alleles",
        nargs="+",
        metavar="HLA",
        help=f"HLA-DRB1 alleles (default: the {len(DEFAULT_ALLELES)}-allele panel)",
    )
    score.add_argument(
        "--no-epitopes",
        action="store_true",
        help="skip MHC II prediction entirely; everything else is offline",
    )
    score.add_argument(
        "--offline",
        action="store_true",
        help="never contact IEDB; fail the epitope section unless it is cached",
    )
    score.add_argument("--cache-dir", help="where IEDB responses are cached")
    score.add_argument("--quiet", action="store_true", help="no progress on stderr")
    score.set_defaults(func=cmd_score)

    validate = sub.add_parser(
        "validate",
        help="score the labelled panel and compare the humanness metrics on it",
        description=(
            "Scores every chain of the committed therapeutic panel, then reports how "
            "well each metric separates human-origin chains from murine and chimeric "
            "ones -- split by chain type, because the pooled answer is driven by heavy "
            "chains. Offline: no IEDB call, no network."
        ),
    )
    validate.add_argument(
        "--figure",
        metavar="PATH",
        help=(
            "write the figure here. A .html path gives a self-contained page with the "
            "table beside it; a .svg path gives the figure alone, for embedding."
        ),
    )
    validate.add_argument("--quiet", action="store_true", help="no progress on stderr")
    validate.set_defaults(func=cmd_validate)

    prov = sub.add_parser(
        "provenance",
        help="print tool versions, settings, verified sources and known sensitivities",
    )
    prov.set_defaults(func=cmd_provenance)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "command", None) == "score" and args.seq and args.fasta:
        parser.error("pass --seq or --fasta, not both")

    try:
        return args.func(args)
    except HmmerNotFoundError as exc:
        # The likeliest failure on a machine that has never run this before,
        # and the one a traceback serves worst: the message already says what
        # to install, so it is printed rather than buried under a stack.
        print(f"abscore: {exc}", file=sys.stderr)
        return EXIT_MISSING_DEPENDENCY


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
