# Cached IEDB responses, for offline tests

Raw TSV responses from the IEDB Analysis Resource MHC class II binding endpoint,
retrieved **2026-11-12** and committed so the test suite runs with no network.

    POST https://tools-cluster-interface.iedb.org/tools_api/mhcii/
    method=netmhciipan_el   length=15   allele=<one HLA-DRB1 allele>

One file per (sequence, allele, method, length), named by the SHA-256 of the
request. Four chains × eight HLA-DRB1 alleles = 32 files:

| Fixture | INN | Origin label |
|---|---|---|
| `TRASTUZUMAB_VH` | 7637 | Humanized |
| `ADALIMUMAB_VH` | 7860 | Human |
| `ABAGOVOMAB_VH` | 8659 | Mouse |

A fourth cached chain is a lambda light chain pulled from OAS during testing
(four somatic mutations from IGLV3-21*02), kept because it is the case that
forced the near-self tier: under exact matching it reported 28 non-self binders,
all of them one residue from a germline peptide.

Sequences are from `panels/imgt_mab_db/` (MIT); see `tests/fixtures.py`.

## Why these are committed rather than fetched

The predictions come from a live public service that updates its models. Without
a frozen copy, a result committed today could quietly stop reproducing after a
server-side change, with nothing to show what moved. With it, `offline=True`
proves a number was replayed rather than re-fetched.

The flip side is that these files record one version of one model at one date.
IEDB does not report a software version in the response, so the retrieval date is
the only handle on which model produced them.

## On the predictor and its licence

These were produced by **NetMHCIIpan 4.1 EL**, run on IEDB's servers. Nothing was
downloaded.

The restrictive CBS/DTU NetMHC academic licence, which limits use to publicly
funded academic institutions and restricts publishing benchmark results, governs
the *downloadable* software. No such software is in this repository. What is
committed here is the output of a public web service for three published
therapeutic antibody sequences.

## Regenerating

Delete a file and re-run with `offline=False`; it will be re-fetched. To rebuild
all of them:

```python
from abscore.epitopes import epitope_burden_for_sequence
import tests.fixtures as F

for name in ("TRASTUZUMAB_VH", "ADALIMUMAB_VH", "ABAGOVOMAB_VH"):
    epitope_burden_for_sequence(getattr(F, name), cache_dir="tests/data/iedb")
```

Expect the numbers in `tests/test_epitopes.py` to need updating if IEDB's model
has changed since 2026-11-12. That the tests would fail is the point.
