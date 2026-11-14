# abscore

Immunogenicity triage for a single antibody heavy or light chain, from sequence
alone. Give it a V-domain; it returns how close that domain sits to human
germline, how much of it has been seen in human antibody repertoires, and which
of its predicted MHC class II binding cores have no match in a human reference.

Every number is computed from public data, and every number that cannot be
reproduced from this checkout says so.

```
$ abscore score --seq EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKG...

CHAIN
────────────────────
  Type                   heavy (H)
  Numbered               120 residues, IMGT 1-128
  Closest germline V     IGHV3-66*04
  Next closest gene      IGHV3-53*04 (1 position(s) behind)

HUMANNESS
────────────────────
  Germline identity, V   83.2%   79/95 positions
    framework            90.0%   FR1+FR2+FR3; FR4 is J-encoded
    CDR1+CDR2            46.7%   CDR3 is junctional, excluded
  OASis identity         70.5%   79/112 9-mers in ≥10% of donors

MHC CLASS II
────────────────────
  17 distinct 9-mer cores predicted to bind, across 86 peptide-allele pairs.
  8 HLA-DRB1 alleles · netmhciipan_el · 15-mers · percentile rank ≤ 10

  Is each core found in a human reference set?

  reference set        found   not found    1 away   2+ away   no entry
  germline Ig genes *      6          10         3         7          1
  OAS repertoire          11           6         2         4          -
  human proteome           6          11         3         8          -
```

---

## Install

The one thing a Python lockfile cannot install is **HMMER**. ANARCI assigns
antibody numbering by shelling out to `hmmscan`, a compiled C program. Install
it first, or nothing will number.

### macOS and Linux

```bash
brew install hmmer           # macOS
sudo apt-get install hmmer   # Debian/Ubuntu
# conda install -c bioconda hmmer

git clone <this repository> && cd abscore
uv sync --frozen --all-extras
uv run abscore score --seq EVQLVESGGG...
```

`--frozen` installs exactly the pinned versions in `uv.lock` and refuses to
change it. ([uv](https://docs.astral.sh/uv/) replaces pip and venv; install it
with `brew install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh`.)

If `hmmscan` is missing, the tool says so and exits 3 with the install command
for your platform, rather than failing somewhere deeper.

### Windows

HMMER has no native Windows build, so on Windows this runs inside Linux: either
**WSL2** (`wsl --install`, then follow the Linux steps) or **Docker**.
This is a constraint of ANARCI, which every Python antibody-numbering project
inherits.

### Docker: any platform, nothing to install

```bash
docker build -t abscore .
docker run --rm abscore score --seq EVQLVESGGG...
cat chains.fasta | docker run --rm -i abscore score
docker run --rm -v "$PWD/out:/out" abscore score --seq EVQ... --save /out/report.html
```

The image carries HMMER, Python 3.11 and the locked dependencies (637 MB). Two
other targets exist: `--target app` adds the Streamlit interface
(`docker run -p 8501:8501 abscore-app`), and `--target test` runs the whole
suite inside the container and fails the build if anything is red.

The CLI image scores a known chain during the build, so an image whose
numbering is broken cannot be produced.

---

## Use

```bash
abscore score --seq EVQLVESGGG...              # one chain
abscore score --fasta chains.fa                # one report per record
cat chains.fa | abscore score                  # or on stdin

abscore score --seq ... --save report.html     # a browser-readable page
abscore score --seq ... --format json          # machine-readable
abscore score --seq ... --format summary       # one line per chain

abscore validate                               # the panel comparison below
abscore provenance                             # versions, settings, limits
uv run streamlit run app.py                    # the interactive version
```

The terminal report is the default. `--save` picks its format from the
extension: `.html` renders the tables and a colour-coded IMGT position track,
`.md` gives markdown.

Progress goes to stderr and results to stdout, so `abscore score ... > out.txt`
writes a clean file while the status lines stay on screen.

Exit codes: `0` scored, `2` bad input, `3` a missing system dependency. IEDB
being unavailable is none of those. The report is still produced, with the
epitope section marked unavailable, because every other number is computed
offline.

---

## What it computes

### Numbering

IMGT numbering via [ANARCI](https://github.com/oxpig/ANARCI). IMGT assigns 128
fixed slots per V-domain and fills them outward from both ends, so a short loop
leaves a gap in the *middle* rather than shifting everything along. That is what
makes position 23 the first cysteine in every antibody regardless of loop
length, and why the region table prints `8/12` for a CDR1 rather than just the
residues.

Five conserved anchors are checked (IMGT 23 Cys, 41 Trp, 89 hydrophobic, 104
Cys, 118 Phe/Trp). An unexpected residue at one of them is reported as a
finding; it means the numbering may be misaligned, or the sequence may carry a
real substitution there, and the report does not guess which.

### Humanness, two ways

**Germline identity.** Percent identity to the closest human germline V gene,
reported three ways: over the whole V region, over framework only (FR1+FR2+FR3;
FR4 is J-encoded), and over CDR1+CDR2. CDR3 is junctional and excluded, since
there is no germline sequence to compare it against. The denominator is
positions present in both query and germline, with coverage reported separately
so a partial domain cannot inflate the figure.

**OASis identity.** The fraction of the chain's 9-mer windows that appear in
human antibody repertoires, as published in
[BioPhi](https://doi.org/10.1080/19420862.2021.2020203) (Prihoda et al., *mAbs*
2022) and computed through `promb`.

Both are reported because they disagree in informative ways, and because on this
project's own validation the simpler one is not worse (below).

### MHC class II epitope burden

15-mer windows are submitted to **IEDB's hosted NetMHCIIpan service**, against a
default panel of 8 HLA-DRB1 alleles covering common supertypes, at a percentile
rank cutoff of 10. Predicted binders are reduced to their distinct 9-mer binding
cores. One core can be presented by several alleles and sit in several
overlapping windows, so the raw hit count overstates how many distinct peptides
are in play.

Each core is then looked up in three human reference sets: germline Ig genes,
the OAS repertoire, and the human proteome. The output reports, per set, how many
cores are **found**, **not found**, and how far the not-found ones sit from the
nearest reference peptide.

A core the reference contains drops down the priority list. A core it lacks has
*still to be shown* to provoke a response. It is a lookup result rather than a
risk score. Cores one substitution away are reported separately and those two or more
away are marked exploratory, because the effect of a substitution depends on
which position changed.

The predictions come from a hosted service, so they are cached locally and the
exact method, allele panel, peptide length and cutoff are recorded in every run.

---

## Validation

`abscore validate` scores every chain of the committed therapeutic panel (1,148
chains from IMGT/mAb-DB) and asks how well each metric separates antibodies of
known human origin from murine and chimeric ones. Nothing in the scoring is
fitted to those origin labels, which are human curation of engineering history.

Separation is reported as **AUC**: the probability that a randomly chosen human
chain outscores a randomly chosen murine or chimeric one. It needs no threshold,
which matters because no threshold here is calibrated.

```
chain type      human   murine/   V identity  FR identity        OASis
                  (n)  chimeric        (AUC)        (AUC)        (AUC)
----------------------------------------------------------------------
all chains        400       152     0.99415*     0.99074      0.98641
heavy only        201        76     0.99764*     0.99241      0.99103
kappa only        150        76     0.99491*     0.99254      0.99452
lambda only        49         0         too few to measure
```

### The result worth reporting is a negative one

**Plain germline identity separates at least as well as OASis, everywhere this
panel can measure it.** OASis is by far the bigger apparatus, a 9-mer prevalence
database built from millions of repertoire sequences, and on this task it does
not beat percent identity to the nearest germline gene.

That is asserted as a test, not as a sentence:
`test_plain_germline_identity_is_not_beaten_by_oasis` fails the build if it ever
stops being true, and a per-chain-type version fails if it stops being true for
heavy or kappa chains separately. The pooled row is heavy-chain-driven; on kappa
chains the margin is **0.00039**, about four pair-swaps out of 11,400, which is
why the table prints five decimals.

OASis still earns its place: it is the only one of the two that says *which*
peptides are unusual, which is what the epitope section needs.

### Where the groups overlap

An AUC of 0.994 hides the tails, so the figure draws every chain as a point and
the command names the exceptions:

```
Most human-like murine chain: tenatumomab (8832_seq1, K) at V=0.817.
4 human-origin chains score at or below it:
  elipovimab    10988_seq2  L  V=0.709  FR=0.771  OASis=0.392
  elipovimab    10988_seq1  H  V=0.779  FR=0.825  OASis=0.524
  navivumab     9976_seq1   H  V=0.781  FR=0.800  OASis=0.695
  tosatoxumab   9757_seq1   H  V=0.812  FR=0.863  OASis=0.743
```

All four are anti-infective antibodies. A plausible reading is that antibodies
reaching their function through heavy somatic hypermutation drift far from their
own germline while remaining entirely human. That would be a limitation of
germline-identity scoring rather than a labelling error. **That reading is not
verified here**: the discovery route of these four was not checked against the
literature, and it is offered as a hypothesis, not a finding.

Separately, one chimeric chain (vibecotamab, V=0.915) sits near the human median
and accounts for most of the distance between the pooled AUC and 1.0.

`abscore validate --figure validation.html` writes the figure and this table as
a self-contained page.

### What the validation does not establish

- **One panel.** It is also the panel the structural framework-gap rule was
  derived from, so the two are not independent. That rule is re-derived from a
  random half and checked on the held-out half, which shows it is stable under
  resampling *within* this panel, and no more than that.
- **Lambda chains are untested.** The panel contains 49 human and 10 humanized
  lambda chains and **zero** chimeric or murine ones, so no lambda separation is
  measurable. A test asserts this, so if a future panel adds them, it fails.
- **Origin is not a clinical endpoint.** A human-origin antibody can still be
  immunogenic.
- **The MHC class II section has no validation at all.** Predicted binders are
  not compared against measured epitopes anywhere in this project.

---

## What this is not

It is a **triage** tool: a way to rank candidates and decide what to look at
first, ahead of an assay.

It does not predict immunogenicity. No score here has been shown to correlate
with clinical anti-drug antibody rates. That is not because the data is
unavailable, but because establishing it is a validation exercise this project
did not attempt. A high humanness score does not mean a molecule is safe, and a
core with no germline match has not been shown to provoke anything. Those statements
appear in the tool's own output for the same reason they appear here.

---

## Limitations

`abscore provenance` prints these in full, with what was measured, what is
mitigated and what remains. The ones that most affect how a number should be
read:

| | |
|---|---|
| **OASis threshold** | `promb` ships only the relaxed peptide set (≥10% of donors). The published OASis correlation with clinical ADA used the medium set (≥50%), which needs BioPhi's full 22 GB database. **Scores here are not comparable to published medium-threshold numbers.** |
| **Germline set choice** | 7.67% of panel chains get a different germline *gene* call from ANARCI's set than from a separate IMGT snapshot. Most are already flagged as ambiguous; ~0.26% change gene without being flagged. |
| **Extraction boundary** | Scores are computed over the extracted V-domain, not the submitted sequence. Where the HMM ends the domain can move by a residue when flanking sequence changes: appending a His-tag to one lambda chain moved OASis from 83.3% to 82.5%. Every report states the extracted length against the input length. |
| **Low-complexity inflation** | OASis counts 9-mer windows, so a long single-residue repeat contributes many identical matching windows. Inserting 100 glycines into a CDR3 raised OASis from 70.5% to 81.1%. The published formula is deliberately left unmodified; the behaviour is pinned by a test instead. |
| **OAS record quality** | OAS contains records with compensated indels that leave part of the translation out of frame while OAS's own flags report them as productive. |
| **IEDB availability** | Predictions come from a hosted service. Responses are cached, a malformed response is refused rather than partially parsed, and an unavailable service marks the section rather than failing the run. |

---

## Data and licences

All inputs are public and used within their terms. Nothing proprietary is
included.

| Source | Use | Terms |
|---|---|---|
| **IMGT** germline sequences (via ANARCI) | numbering, germline assignment | CC BY-4.0 |
| **IMGT/mAb-DB** therapeutic panel, via [`Merck/BioPhi-2021-publication`](https://github.com/Merck/BioPhi-2021-publication) | validation panel (1,148 chains) | MIT, redistributed with `LICENSE.biophi` and SHA-256 checksums in `panels/imgt_mab_db/SOURCE.md` |
| **OAS** repertoire sequences | numbering spot-checks | OAS publishes no licence or terms, so sequences were read locally and **are not redistributed here** |
| **`promb`** peptide database | OASis humanness | MIT |
| **IEDB** NetMHCIIpan | class II prediction | called as a hosted service, so the CBS/DTU-licensed software is not redistributed |

---

## Development

```bash
uv sync --frozen --all-extras
uv run pytest                      # 558 tests, offline
uv run pytest -m "not panel"       # skip the 1,148-chain panel tests
uv run ruff check . && uv run ruff format --check .
```

The suite touches no network: every IEDB response it needs is committed under
`tests/data/iedb` and every invocation passes `--offline`.

CI runs the same suite on clean amd64 Linux, plus the container build, on every
push. That is what checks these install instructions actually work. A passing
suite on the author's machine proves only that the author's machine works.

### Layout

```
src/abscore/
  number.py       IMGT numbering, anchors, completeness checks
  germline.py     closest germline V gene, identity by scope
  humanness.py    OASis 9-mer scoring
  epitopes.py     binding cores, reference-set lookup, distances
  panel.py        loading the committed therapeutic panel
  iedb.py         the hosted predictor: requests, caching, strict parsing
  report.py       one report, three renderers (text, markdown, JSON)
  html.py         self-contained HTML report
  track.py        the colour-coded IMGT position track
  figure.py       the validation figure
  validation.py   panel scoring and the metric comparison
  provenance.py   versions, settings, verified sources, known limits
  cli.py          score · validate · provenance
app.py            Streamlit interface
panels/           committed therapeutic panel (+ SOURCE.md, licence, checksums)
tests/            558 tests, including committed IEDB fixtures
```

---

## Scope

Built as a take-home exercise, to roughly a day of work, using public data only.
The design choices that follow from that (one panel, no clinical validation, a
hosted predictor rather than a local one) are stated where they matter rather
than left for a reader to discover.
