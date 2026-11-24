"""Occurrence TSV must not repeat core DwC column names from metadata."""

from __future__ import annotations

import csv
from pathlib import Path

from bin.build_darwin_core import (
    OCCURRENCE_FIELDS,
    build_dna_derived_data_table,
    build_occurrence_table,
    occurrence_id_for,
    occurrence_metadata_fields,
)


def test_occurrence_metadata_fields_excludes_core_occurrence_columns() -> None:
    metadata = {
        "555111": {
            "sampleID": "555111",
            "eventID": "555111",
            "basisOfRecord": "MaterialSample",
            "locality": "Nice",
            "lib_layout": "single",
            "decimalLatitude": "43.7",
            "assay_name": "MiFish",
        }
    }
    fields = occurrence_metadata_fields(metadata)
    assert "eventID" not in fields
    assert "basisOfRecord" not in fields
    assert "sampleID" not in fields
    assert "lib_layout" not in fields  # DNA-derived extension field
    assert "assay_name" not in fields  # used only in occurrenceID
    assert "locality" in fields
    assert "decimalLatitude" in fields
    assert not (set(fields) & set(OCCURRENCE_FIELDS))


def test_occurrence_id_includes_assay_name() -> None:
    assert occurrence_id_for("OBIS_EDNA_1", "abc", "MiFish") == "OBIS_EDNA_1_MiFish_abc"
    assert occurrence_id_for("OBIS_EDNA_1", "abc", "  ") == "OBIS_EDNA_1_abc"


def test_occurrence_and_dna_share_assay_qualified_id(tmp_path: Path) -> None:
    sintax = [{"ASV_ID": "abc", "Genus": "Gadus", "Species": ""}]
    dada2 = [{"ASV_ID": "abc", "sequence": "ACGT", "OBIS_EDNA_1": "12"}]
    metadata = {
        "OBIS_EDNA_1": {
            "assay_name": "GC",
            "locality": "Reef",
            "target_gene": "ITS2",
        }
    }
    fields = occurrence_metadata_fields(metadata)
    occurrence_path = tmp_path / "occurrence.tsv"
    index = build_occurrence_table(
        sintax,
        [],
        dada2,
        metadata,
        fields,
        {},
        occurrence_path,
        clean_prefix=False,
    )
    dna_path = tmp_path / "dnaderiveddata.tsv"
    build_dna_derived_data_table(
        dada2,
        index,
        metadata,
        ("target_gene",),
        dna_path,
    )

    with occurrence_path.open(newline="") as handle:
        occurrence_rows = list(csv.DictReader(handle, delimiter="\t"))
    with dna_path.open(newline="") as handle:
        dna_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert occurrence_rows[0]["occurrenceID"] == "OBIS_EDNA_1_GC_abc"
    assert "assay_name" not in occurrence_rows[0]
    assert occurrence_rows[0]["locality"] == "Reef"
    assert dna_rows[0]["occurrenceID"] == "OBIS_EDNA_1_GC_abc"
    assert dna_rows[0]["DNA_sequence"] == "ACGT"
    assert dna_rows[0]["target_gene"] == "ITS2"
