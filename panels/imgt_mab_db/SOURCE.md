# Therapeutic antibody panel — source and provenance

Redistributed here under the MIT License (see `LICENSE.biophi`). Retrieved **2026-11-12**.

## Where these files came from

[`Merck/BioPhi-2021-publication`](https://github.com/Merck/BioPhi-2021-publication),
the analysis repository for:

> Prihoda D, Maamary J, Waight A, Juan V, Fayadat-Dilman L, Svozil D, Bitton DA.
> **BioPhi: A platform for antibody design, humanization, and humanness evaluation
> based on natural antibody repertoires and deep learning.**
> *mAbs* 2022;14(1):2020203. doi:[10.1080/19420862.2021.2020203](https://doi.org/10.1080/19420862.2021.2020203)

That repository derived them in turn from **IMGT/mAb-DB** (<https://www.imgt.org>),
and credits **Thera-SAbDab** (OPIG) in its acknowledgements.

| File | Upstream path | SHA-256 |
|---|---|---|
| `imgt_mab_db.fa` | `data/tasks/humanness/imgt_mab_db.fa` | `0ea81cee903f925c53ffc894ee479ef94e5b0d8e937b987a161ce5178a38ab75` |
| `imgt_mab_db.tsv` | `data/tasks/humanness/imgt_mab_db.tsv` | `d31ae16b6a56bbfb3ca72f1c8c33ec47b68abd1c1eb0d88a3bf02e8c7249a997` |
| `immunogenicity.tsv` | `data/tasks/immunogenicity/immunogenicity.tsv` | `417c8b70790e4b1c2a8a053435466210c248972307abbf52d562964809aaf8a7` |
| `LICENSE.biophi` | `LICENSE` | — |

## Contents

### `imgt_mab_db.tsv` — 564 therapeutics

`INN_Number`, `Antibody` (INN name), `Species` (origin label).

| `Species` label | n antibodies |
|---|---|
| Humanized | 229 |
| Human | 198 |
| Chimeric | 63 |
| Humanized/Chimeric | 41 |
| Mouse | 13 |
| Human/Humanized | 7 |
| Caninized | 6 |
| Human/Chimeric | 4 |
| Felinized | 3 |

The compound labels (`Humanized/Chimeric`, `Human/Humanized`, …) describe pairs
whose two chains have different origins. They are kept rather than collapsed,
because forcing them into one class would invent a label the source does not make.
The veterinary labels (`Caninized`, `Felinized`) are likewise kept: they are a
useful negative control, since a caninized antibody should score as *less* human
than a murine one, and does.

### `imgt_mab_db.fa` — 1148 sequences

FASTA keyed `<INN_Number>_seq1` / `_seq2`. **Chain order is not fixed** — for some
entries `seq1` is the heavy chain, for others the light. Chain type must be
determined by numbering, never assumed from the suffix. 1146 of the 1148 number
as an antibody V-domain.

### `immunogenicity.tsv` — 217 antibodies with clinical ADA

`Name`, `Immunogenicity` (reported ADA %), `Species`, `Reference` (a literature
URL per antibody; some contain embedded newlines, so `wc -l` reports 228 lines for
217 records). Curated in the Hu-mAb study (Marks et al., *Bioinformatics*
2021) and used in the BioPhi paper's ADA correlation, where OASis medium reached
R = −0.53 and plain germline content R = −0.50.

## What this panel is used for here

1. **Deriving the structural framework-gap whitelist** in `abscore.number` —
   which framework IMGT positions are empty by design per chain type. Derived on
   a random half and checked against the held-out half; see
   `tests/test_panel.py`.
2. **Validation**: whether germline identity and OASis separate the origin
   classes, and which separates better.
3. **Test fixtures** in `tests/fixtures.py`.

## Limitation to state when citing results from it

Items 1–3 all draw on this single source, so a test fixture, the rule it checks,
and the panel it is validated against are not independent of one another. The
split in `tests/test_panel.py` shows the gap rule is stable under resampling
*within* this panel; it does not show it generalises to sequences from elsewhere.
The one genuinely independent check in this project is an OAS record numbered by
OAS's own ANARCI run, reproduced position-for-position
(`tests/test_number.py::TestFrameworkDeletions`).
