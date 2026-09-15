"""Regression tests for invalid MIDORI-style species labels."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from bin.build_darwin_core import identification_remarks, merge_asv_taxonomy, species_label
from bin.worms_match import clean_missing_rank_values, enrich_rows_with_worms, finest_taxon

# From MIDORI header:
# >HQ270472;tax=d:Eukaryota,p:NA,c:Dinophyceae,o:Gymnodiniales,f:Gymnodiniaceae,g:Gymnodinium,s:Gymnodinium_sp._NVA/RUS/2008
GYMNODINIUM_ROW = {
    "Kingdom": "Eukaryota",
    "Phylum": "NA",
    "Class": "Dinophyceae",
    "Order": "Gymnodiniales",
    "Family": "Gymnodiniaceae",
    "Genus": "Gymnodinium",
    "Species": "Gymnodinium_sp._NVA/RUS/2008",
}

XFAIL_REASON = "species placeholder filtering not implemented"


@pytest.mark.xfail(strict=True, reason=XFAIL_REASON)
def test_finest_taxon_skips_gymnodinium_sp_placeholder() -> None:
    name, rank = finest_taxon(GYMNODINIUM_ROW)
    assert (name, rank) == ("Gymnodinium", "Genus")


@pytest.mark.xfail(strict=True, reason=XFAIL_REASON)
def test_species_label_rejects_gymnodinium_sp_placeholder() -> None:
    assert species_label({"Species": "Gymnodinium_sp._NVA/RUS/2008"}) == ""


@pytest.mark.xfail(strict=True, reason=XFAIL_REASON)
def test_merge_asv_taxonomy_omits_gymnodinium_sp_placeholder() -> None:
    merged = merge_asv_taxonomy(GYMNODINIUM_ROW, GYMNODINIUM_ROW)
    assert merged["species"] == ""


def test_finest_taxon_falls_back_to_genus_when_species_na() -> None:
    name, rank = finest_taxon({"Genus": "Gymnodinium", "Species": "NA"})
    assert (name, rank) == ("Gymnodinium", "Genus")


def test_clean_missing_rank_values_strips_na_from_rank_columns() -> None:
    row = {
        "Kingdom": "Eukaryota",
        "Phylum": "NA",
        "Class": "n/a",
        "Order": "unassigned",
        "Family": "Gymnodiniaceae",
        "Genus": "Gymnodinium",
        "Species": "NA",
        "confidence": "1.00",
    }
    cleaned = clean_missing_rank_values([row])[0]
    assert cleaned["Phylum"] == ""
    assert cleaned["Class"] == ""
    assert cleaned["Order"] == ""
    assert cleaned["Species"] == ""
    assert cleaned["Kingdom"] == "Eukaryota"
    assert cleaned["Family"] == "Gymnodiniaceae"
    assert cleaned["Genus"] == "Gymnodinium"
    assert cleaned["confidence"] == "1.00"


def test_clean_missing_rank_values_handles_lowercase_occurrence_ranks() -> None:
    row = {
        "kingdom": "Eukaryota",
        "phylum": "NA",
        "genus": "Gymnodinium",
        "species": "unassigned",
        "occurrenceID": "s1_asv1",
    }
    cleaned = clean_missing_rank_values([row])[0]
    assert cleaned["phylum"] == ""
    assert cleaned["species"] == ""
    assert cleaned["occurrenceID"] == "s1_asv1"


def test_merge_asv_taxonomy_strips_na_rank_values() -> None:
    sintax_row = {"Kingdom": "Eukaryota", "Phylum": "NA", "Genus": "Gadus", "Species": ""}
    vsearch_row = {"Genus": "Gadus", "Species": "morhua"}
    merged = merge_asv_taxonomy(sintax_row, vsearch_row)
    assert merged["phylum"] == ""
    assert merged["species"] == "morhua"


def test_merge_asv_taxonomy_uses_vsearch_species_when_genera_match() -> None:
    sintax_row = {"Genus": "Gadus", "Species": ""}
    vsearch_row = {"Genus": "Gadus", "Species": "morhua"}
    merged = merge_asv_taxonomy(sintax_row, vsearch_row)
    assert merged["species"] == "morhua"
    assert merged["identificationRemarks"] == ""


def test_merge_asv_taxonomy_identification_remarks_includes_hits_for_species() -> None:
    sintax_row = {"Genus": "Gadus", "Species": "Gadus_morhua"}
    vsearch_row = {"Genus": "Gadus", "Species": "morhua"}
    merged = merge_asv_taxonomy(
        sintax_row,
        vsearch_row,
        raw_sintax="g:Gadus(1.00),s:Gadus_morhua(0.90)",
        hits=["OP288117", "KJ763881"],
    )
    assert merged["species"] == "morhua"
    assert merged["identificationRemarks"] == (
        "SINTAX species assignment was: Gadus_morhua; "
        "SINTAX confidence: g:Gadus(1.00),s:Gadus_morhua(0.90); "
        "VSEARCH hits: OP288117,KJ763881"
    )


def test_merge_asv_taxonomy_omits_hits_when_species_not_assigned() -> None:
    sintax_row = {"Genus": "Eukaryota", "Species": ""}
    vsearch_row = {"Genus": "Gymnodinium", "Species": "Gymnodinium_sp._NVA/RUS/2008"}
    merged = merge_asv_taxonomy(
        sintax_row, vsearch_row, hits=["OP288117", "KJ763881"]
    )
    assert merged["species"] == ""
    assert merged["identificationRemarks"] == ""


def test_identification_remarks_includes_hits() -> None:
    assert identification_remarks("", hits=["OP288117", "KJ763881"]) == (
        "VSEARCH hits: OP288117,KJ763881"
    )


def test_accession_from_blast6_target_and_load_hits(tmp_path: Path) -> None:
    from bin.build_darwin_core import (
        accession_from_blast6_target,
        load_vsearch_hits_by_asv,
    )

    assert (
        accession_from_blast6_target(
            "OP288117;tax=d:Eukaryota,s:Strombidium_biarmatum"
        )
        == "OP288117"
    )
    assert accession_from_blast6_target("*") == ""

    hits_path = tmp_path / "ASV_tax_vsearch_lca.user.txt"
    hits_path.write_text(
        "asv1\tOP288117;tax=d:Eukaryota\t100.0\t10\t0\t0\t1\t10\t1\t10\t-1\t0\n"
        "asv1\tKJ763881;tax=d:Eukaryota\t100.0\t10\t0\t0\t1\t10\t1\t10\t-1\t0\n"
        "asv1\tOP288117;tax=d:Eukaryota\t100.0\t10\t0\t0\t1\t10\t1\t10\t-1\t0\n"
        "asv2\t*\t0.0\t0\t0\t0\t0\t0\t0\t0\t-1\t0\n"
    )
    assert load_vsearch_hits_by_asv(hits_path) == {
        "asv1": ["OP288117", "KJ763881"],
    }


def test_merge_asv_taxonomy_identification_remarks_includes_sintax_species_and_raw() -> None:
    sintax_row = {"Genus": "Homo", "Species": "Homo sapiens_9606"}
    raw_sintax = (
        "k:Eukaryota_2759(1.00),p:Chordata_7711(1.00),c:Mammalia_40674(1.00),"
        "o:Primates_9443(1.00),f:Hominidae_9604(1.00),g:Homo_9605(1.00),"
        "s:Homo sapiens_9606(1.00)"
    )
    merged = merge_asv_taxonomy(sintax_row, None, raw_sintax=raw_sintax)
    assert merged["identificationRemarks"] == (
        "SINTAX species assignment was: Homo sapiens_9606; "
        "SINTAX confidence: k:Eukaryota_2759(1.00),p:Chordata_7711(1.00),c:Mammalia_40674(1.00),"
        "o:Primates_9443(1.00),f:Hominidae_9604(1.00),g:Homo_9605(1.00),"
        "s:Homo sapiens_9606(1.00)"
    )


def test_identification_remarks_includes_confidence_without_species() -> None:
    raw_sintax = "k:Eukaryota_2759(1.00),p:Chordata_7711(1.00)"
    assert identification_remarks("", raw_sintax) == (
        "SINTAX confidence: k:Eukaryota_2759(1.00),p:Chordata_7711(1.00)"
    )


def test_identification_remarks_omits_raw_when_missing() -> None:
    assert identification_remarks("Homo sapiens_9606") == (
        "SINTAX species assignment was: Homo sapiens_9606"
    )


def test_merge_asv_taxonomy_drops_vsearch_species_when_genera_differ() -> None:
    sintax_row = {"Genus": "Eukaryota", "Species": ""}
    vsearch_row = {"Genus": "Gymnodinium", "Species": "Gymnodinium_sp._NVA/RUS/2008"}
    merged = merge_asv_taxonomy(sintax_row, vsearch_row)
    assert merged["species"] == ""
    assert merged["identificationRemarks"] == ""


def test_merge_asv_taxonomy_drops_vsearch_species_when_sintax_genus_missing() -> None:
    sintax_row = {"Kingdom": "Eukaryota", "Genus": "", "Species": ""}
    vsearch_row = {"Genus": "Gymnodinium", "Species": "Gymnodinium_sp._NVA/RUS/2008"}
    merged = merge_asv_taxonomy(sintax_row, vsearch_row)
    assert merged["species"] == ""
    assert merged["identificationRemarks"] == ""


@pytest.mark.xfail(strict=True, reason=XFAIL_REASON)
def test_enrich_rows_with_worms_skips_gymnodinium_sp_placeholder() -> None:
    def fake_match(taxa, worms_db):
        return {
            (name, rank): {
                "matched_aphiaid": 109452 if rank == "Genus" else None,
                "matched_name": name.split("_")[0] if rank == "Genus" else name,
                "canonical": name.split("_")[0] if rank == "Genus" else name,
            }
            for name, rank in taxa
        }

    with patch("bin.worms_match.match_distinct_taxa", side_effect=fake_match):
        enriched = enrich_rows_with_worms([GYMNODINIUM_ROW], worms_db=Path("/dev/null"))

    assert enriched[0]["taxonRank"] == "genus"
    assert enriched[0]["scientificName"] == "Gymnodinium"
    assert enriched[0]["scientificNameID"] == "urn:lsid:marinespecies.org:taxname:109452"
