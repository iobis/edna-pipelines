"""Regression tests for invalid MIDORI-style species labels."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from bin.build_darwin_core import merge_asv_taxonomy, species_label
from bin.worms_match import enrich_rows_with_worms, finest_taxon

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
