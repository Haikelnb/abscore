# abscore

Immunogenicity triage for a single antibody heavy or light chain, from sequence
alone. Give it a V-domain; it returns how close that domain sits to human
germline, how much of it has been seen in human antibody repertoires, and which
of its predicted MHC class II binding cores have no match in a human reference.

It ranks candidates ahead of an assay. It does not predict immunogenicity: no
score here has been shown to correlate with clinical anti-drug antibody rates,
and establishing that is a validation exercise this project did not attempt.

Every number is computed from public data, and every number that cannot be
reproduced from this checkout says so.

```
$ abscore score --offline --cache-dir tests/data/iedb \
    --fasta examples/adalimumab.fasta

CHAIN
────────────────────
  Type                   heavy (H)
  Numbered               121 residues, IMGT 1-128
  Closest germline V     IGHV3-9*01
  Next closest gene      IGHV3-20*04 (7 position(s) behind)

HUMANNESS
────────────────────
  Germline identity, V   93.8%   90/96 positions
    framework            95.0%   FR1+FR2+FR3; FR4 is J-encoded
    CDR1+CDR2            87.5%   CDR3 is junctional, excluded
  OASis identity         85.8%   97/113 9-mers in ≥10% of donors

MHC CLASS II
────────────────────
  18 distinct 9-mer cores predicted to bind, across 77 peptide-allele pairs.
  8 HLA-DRB1 alleles · netmhciipan_el · 15-mers · percentile rank ≤ 10

  Is each core found in a human reference set?

  reference set        found   not found    1 away   2+ away   no entry
  germline Ig genes *     11           3         1         2          4
  OAS repertoire          15           3         1         2          -
  human proteome          11           7         1         6          -
```

That is the first of the two records, adalimumab's heavy chain; the kappa chain
follows it, with its MHC class II section marked unavailable: only the heavy
chain's predictions are in the committed cache, and `--offline` never calls out.
The command otherwise reproduces this from a clean checkout with no network. The
full report also prints a per-region IMGT track and a caveats section.

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

uv run abscore score --fasta examples/adalimumab.fasta --no-epitopes --format summary
```

That last command is fully offline and should print:

```
adalimumab_heavy	chain=H germline=IGHV3-9*01 V=93.8% FR=95.0% OASis=85.8% binders=unavailable
adalimumab_kappa	chain=K germline=IGKV1-27*03 V=97.7% FR=98.7% OASis=97.0% binders=unavailable
```

`--frozen` installs exactly the pinned versions in `uv.lock` and refuses to
change it. ([uv](https://docs.astral.sh/uv/) replaces pip and venv; install it
with `brew install uv` or `curl -LsSf https://astral.sh/uv/install.sh | sh`.)

If `hmmscan` is missing, the tool says so and exits 3 with the install command
for your platform, rather than failing somewhere deeper.

### Windows, and Docker

HMMER has no native Windows build, so on Windows this runs inside Linux: either
**WSL2** (`wsl --install`, then follow the Linux steps) or Docker. This is a
constraint of ANARCI, which every Python antibody-numbering project inherits.

```bash
docker build -t abscore .
docker run --rm abscore score --fasta examples/adalimumab.fasta --no-epitopes
```

The image carries HMMER, Python 3.11 and the locked dependencies (637 MB), and
scores a known chain during the build, so an image whose numbering is broken
cannot be produced. `--target app` adds the Streamlit interface and
`--target test` runs the suite inside the container.

---

## Use

```bash
abscore score --seq EVQ...VSS                  # one chain
abscore score --fasta chains.fa                # one report per record
cat chains.fa | abscore score                  # or on stdin

abscore score --fasta chains.fa --format summary   # one line per record
abscore score --fasta chains.fa --format json      # machine-readable
abscore score --fasta chains.fa --save report.html # all records, one page
abscore score --fasta chains.fa --out results/     # one file per record

abscore validate                               # the panel comparison below
abscore provenance                             # versions, settings, limits
uv run streamlit run app.py                    # the interactive version
```

`--save` picks its format from the extension: `.html` renders the tables and a
colour-coded IMGT position track, `.md` gives markdown.

Each record is submitted to the MHC class II predictor separately, so a large
FASTA is network-bound. For batch work use `--no-epitopes`, or a warm
`--cache-dir`.

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
leaves a gap in the *middle* rather than shifting everything along; that is why
the region table prints `8/12` for a CDR1 rather than just the residues. Five
conserved anchors are checked (IMGT 23 Cys, 41 Trp, 89 hydrophobic, 104 Cys,
118 Phe/Trp). An unexpected residue at one is reported as a finding, since it
means either the numbering is misaligned or the sequence carries a real
substitution, and the report does not guess which.

### Humanness, two ways

**Germline identity.** Percent identity to the closest human germline V gene,
reported over the whole V region, over framework only (FR1+FR2+FR3; FR4 is
J-encoded), and over CDR1+CDR2. CDR3 is junctional and excluded, since there is
no germline sequence to compare it against. The denominator is positions present
in both query and germline, with coverage reported separately so a partial
domain cannot inflate the figure.

**OASis identity.** The fraction of the chain's 9-mer windows that appear in
human antibody repertoires, as published in
[BioPhi](https://doi.org/10.1080/19420862.2021.2020203) (Prihoda et al., *mAbs*
2022) and computed through `promb`.

Both are reported because they answer related but different questions: a V-gene
identity does not describe the junctional region, and a peptide prevalence does
not say which gene a chain came from.

### MHC class II epitope burden

15-mer windows are submitted to **IEDB's hosted NetMHCIIpan service**
(`netmhciipan_el`, eluted-ligand likelihood rather than binding affinity),
against a default panel of 8 HLA-DRB1 alleles covering common supertypes, at a
percentile rank cutoff of 10. Predicted binders are reduced to their distinct
9-mer binding cores. One core can be presented by several alleles and sit in
several overlapping windows, so the raw hit count overstates how many distinct
peptides are in play.

Each core is then looked up in three human reference sets: germline Ig genes,
the OAS repertoire, and the human proteome. The output reports, per set, how
many cores are **found**, **not found**, and how far the not-found ones sit from
the nearest reference peptide. This is an exact sequence annotation with unknown
relevance to patient-level tolerance. A core the reference contains drops down
the priority list; a core it lacks has *still to be shown* to provoke a
response. Distance is not a risk scale, so cores two or more substitutions away
are marked exploratory.

Three reference sets are exposed rather than one because the choice of reference
changes the answer, and that sensitivity is the point. Predictions come from a
hosted service, so they are cached locally and the exact method, allele panel,
peptide length and cutoff are recorded in every run.

---

## Validation

`abscore validate` scores every chain of the committed therapeutic panel (1,148
chains from IMGT/mAb-DB) and asks how well each metric separates antibodies of
known human origin from murine and chimeric ones.

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

Both scores discriminate origin well on this panel. Germline V identity had the
highest observed AUC; OASis also provides local peptide and junctional-region
information that a V-gene identity cannot. The differences between columns are
observed point estimates on one panel, not established superiority of one
measure over another: on kappa chains the margin is **0.00039**, about four
pair-swaps out of 11,400, which is why the table prints five decimals.

An AUC of 0.994 hides the tails, so the figure draws every chain as a point:

![Humanness scores by engineering origin and chain type](docs/validation.svg)

Three metrics across, three chain types down, with each group's median as a flat
rule. The overlap is small and nameable: one murine chain (tenatumomab,
V=0.817) outscores four human-origin chains, the lowest being elipovimab's
lambda chain at V=0.709, and one chimeric chain (vibecotamab, V=0.915) sits near
the human median and is most of the distance between the pooled AUC and 1.0.
`abscore validate` prints those by name; `--figure validation.html` writes the
figure and this table as a self-contained page, `--figure validation.svg`
regenerates the image above.

### What the validation does not establish

- **Origin is not a clinical endpoint.** Separating human from murine chains is
  not predicting anti-drug antibody responses. A human-origin antibody can still
  be immunogenic.
- **One panel.** It is also the panel the structural framework-gap rule was
  derived from, so the two are not independent. That rule is re-derived from a
  random half and checked on the held-out half, which shows it is stable under
  resampling *within* this panel, and no more than that.
- **No lambda comparison.** Origin discrimination could not be evaluated for
  lambda: the panel holds 49 human and 10 humanized lambda chains and **zero**
  chimeric or murine ones. Lambda implementation behaviour is tested; the
  comparison is not measurable. A test asserts this, so if a future panel adds
  negatives, it fails.
- **MHC class II is tested, not validated.** The integration is tested against
  stored IEDB responses, but this project does not independently check its
  predictions or reference annotations against experimental binding,
  presentation or T-cell assays.

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
| **IMGT/mAb-DB** therapeutic panel, via [`Merck/BioPhi-2021-publication`](https://github.com/Merck/BioPhi-2021-publication) | validation panel (1,148 chains), `examples/` | MIT, redistributed with `LICENSE.biophi` and SHA-256 checksums in `panels/imgt_mab_db/SOURCE.md` |
| **[OAS](https://opig.stats.ox.ac.uk/webapps/oas/)** repertoire sequences | numbering spot-checks | CC BY-4.0. Redistributable with attribution, but the two data units used are tens of megabytes and nothing here recomputes from them, so they were read locally and **are not shipped** |
| **`promb`** peptide database | OASis humanness | MIT |
| **IEDB** NetMHCIIpan | class II prediction | called as a hosted service, so the CBS/DTU-licensed software is not redistributed |
