"""Regression tests for invalid MIDORI-style species labels."""

from __future__ import annotations

import pytest

from bin.build_darwin_core import merge_asv_taxonomy, species_label
from bin.worms_match import finest_taxon

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
