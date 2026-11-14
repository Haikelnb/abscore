"""Loading the committed therapeutic antibody panel.

The panel is the only labelled data in the project, and it is used for three
things that are not independent of each other: deriving the structural
framework-gap whitelist, validating the scores, and supplying test fixtures. See
`panels/imgt_mab_db/SOURCE.md` for provenance and for what that non-independence
does and does not permit.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

PANEL_DIR = Path(__file__).resolve().parents[2] / "panels" / "imgt_mab_db"
FASTA_PATH = PANEL_DIR / "imgt_mab_db.fa"
LABELS_PATH = PANEL_DIR / "imgt_mab_db.tsv"
ADA_PATH = PANEL_DIR / "immunogenicity.tsv"

#: Origin labels for antibodies whose two chains share one origin.
SINGLE_ORIGIN_LABELS = ("Human", "Humanized", "Chimeric", "Mouse")

#: Veterinary labels. A caninized antibody should score as less human than a
#: murine one, which makes these a useful negative control rather than noise.
VETERINARY_LABELS = ("Caninized", "Felinized")


class PanelError(FileNotFoundError):
    """Raised when the committed panel is missing."""


@dataclass(frozen=True)
class PanelEntry:
    """One chain of one therapeutic antibody.

    `record_id` is the FASTA key, e.g. `7637_seq1`. The `seq1`/`seq2` suffix does
    NOT indicate chain type -- the source files are not consistent about which
    chain comes first, so chain type must come from numbering.
    """

    record_id: str
    inn_number: str
    name: str
    origin: str
    sequence: str

    @property
    def is_single_origin(self) -> bool:
        return self.origin in SINGLE_ORIGIN_LABELS

    @property
    def is_veterinary(self) -> bool:
        return self.origin in VETERINARY_LABELS


def _require(path: Path) -> Path:
    if not path.exists():
        raise PanelError(
            f"panel file missing: {path}. See panels/imgt_mab_db/SOURCE.md for how to restore it."
        )
    return path


@lru_cache(maxsize=1)
def load_labels() -> dict[str, tuple[str, str]]:
    """INN number -> (name, origin label)."""
    with _require(LABELS_PATH).open() as handle:
        return {
            row["INN_Number"]: (row["Antibody"], row["Species"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


@lru_cache(maxsize=1)
def load_sequences() -> dict[str, str]:
    """FASTA record id -> amino acid sequence."""
    sequences: dict[str, str] = {}
    key: str | None = None
    parts: list[str] = []
    with _require(FASTA_PATH).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if key is not None:
                    sequences[key] = "".join(parts)
                key, parts = line[1:], []
            else:
                parts.append(line)
    if key is not None:
        sequences[key] = "".join(parts)
    return sequences


@lru_cache(maxsize=1)
def load_panel() -> tuple[PanelEntry, ...]:
    """Every panel chain, sorted by record id so iteration order is stable."""
    labels = load_labels()
    entries = []
    for record_id, sequence in sorted(load_sequences().items()):
        inn = record_id.split("_")[0]
        name, origin = labels.get(inn, ("?", "?"))
        entries.append(
            PanelEntry(
                record_id=record_id,
                inn_number=inn,
                name=name,
                origin=origin,
                sequence=sequence,
            )
        )
    return tuple(entries)


@lru_cache(maxsize=1)
def load_ada() -> dict[str, float]:
    """Antibody name -> reported clinical ADA percentage.

    227 antibodies, curated in the Hu-mAb study and used for the BioPhi ADA
    correlation. ADA rates are not comparable across assays, trials, dosing
    routes or patient populations, so these are a weak outcome measure and any
    correlation against them has to be reported as such.
    """
    with _require(ADA_PATH).open() as handle:
        return {
            row["Name"]: float(row["Immunogenicity"])
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("Immunogenicity")
        }
