"""Occurrence TSV must not repeat core DwC column names from metadata."""

from __future__ import annotations

from bin.build_darwin_core import OCCURRENCE_FIELDS, occurrence_metadata_fields


def test_occurrence_metadata_fields_excludes_core_occurrence_columns() -> None:
    metadata = {
        "555111": {
            "sampleID": "555111",
            "eventID": "555111",
            "basisOfRecord": "MaterialSample",
            "locality": "Nice",
            "lib_layout": "single",
            "decimalLatitude": "43.7",
        }
    }
    fields = occurrence_metadata_fields(metadata)
    assert "eventID" not in fields
    assert "basisOfRecord" not in fields
    assert "sampleID" not in fields
    assert "lib_layout" not in fields  # DNA-derived extension field
    assert "locality" in fields
    assert "decimalLatitude" in fields
    assert not (set(fields) & set(OCCURRENCE_FIELDS))
