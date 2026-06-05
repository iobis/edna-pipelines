#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path


INPUT_RANK_COLUMNS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
OUTPUT_RANK_COLUMNS = tuple(rank.lower() for rank in INPUT_RANK_COLUMNS)
OCCURRENCE_FIELDS = (
    "occurrenceID",
    "eventID",
    "basisOfRecord",
    "taxonID",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "identificationRemarks",
    "organismQuantity",
    "organismQuantityType",
    "sample_id"
)
DNA_DERIVED_DATA_FIELDS = (
    "occurrenceID",
    "DNA_sequence",
)
DNA_DERIVED_DATA_METADATA_FIELDS = (
    "lib_layout",
    "seq_meth",
    "samp_name",
    "samp_size",
    "target_gene",
    "ref_db",
    "tax_class",
    "pcr_primer_forward",
    "pcr_primer_reverse",
    "pcr_primer_reference",
    "env_medium",
)
METADATA_SAMPLE_ID_COLUMNS = ("sampleID", "sample_id", "sample")
SAMPLE_ID_PREFIX = "s_"


def clean_sample_id(sample_id: str, *, clean_prefix: bool) -> str:
    if clean_prefix and sample_id.startswith(SAMPLE_ID_PREFIX):
        return sample_id[len(SAMPLE_ID_PREFIX) :]
    return sample_id


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    by_sample: dict[str, dict[str, str]] = {}
    for row in load_tsv(path):
        sample_id = ""
        for column in METADATA_SAMPLE_ID_COLUMNS:
            sample_id = (row.get(column) or "").strip()
            if sample_id:
                break
        if sample_id:
            by_sample[sample_id] = row
    return by_sample


def metadata_dna_derived_fields(
    metadata_by_sample: dict[str, dict[str, str]],
) -> tuple[str, ...]:
    available: set[str] = set()
    for row in metadata_by_sample.values():
        available.update(row.keys())
    return tuple(
        field for field in DNA_DERIVED_DATA_METADATA_FIELDS if field in available
    )


def find_dada2_table(ampliseq_root: Path) -> Path:
    path = ampliseq_root / "dada2" / "DADA2_table.tsv"
    if not path.is_file():
        raise FileNotFoundError(f"DADA2 table not found: {path}")
    return path


def find_pipeline_params_json(ampliseq_root: Path) -> Path:
    info_dir = ampliseq_root / "pipeline_info"
    if not info_dir.is_dir():
        raise FileNotFoundError(f"pipeline_info/ missing under {ampliseq_root}")

    matches = sorted(
        info_dir.glob("params_*.json"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    if not matches:
        raise FileNotFoundError(f"no params_*.json found in {info_dir}")
    return matches[0]


def find_worms_matched(output_dir: Path, method: str) -> Path:
    path = output_dir / "worms" / method / f"worms_matched.{method}.tsv"
    if not path.is_file():
        raise FileNotFoundError(
            f"WoRMS-matched TSV not found: {path} (run WORMS_MATCH first)"
        )
    return path


def write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def species_label(row: dict[str, str]) -> str:
    return (row.get("Species") or "").strip()


def identification_remarks(sintax_species: str, vsearch_species: str) -> str:
    remarks: list[str] = []
    if sintax_species and not vsearch_species:
        remarks.append(
            "SINTAX species assignment removed (not replaced by VSEARCH): "
            f"{sintax_species}"
        )
    elif vsearch_species:
        if sintax_species:
            remarks.append(
                f"VSEARCH species assignment: {vsearch_species}; "
                f"SINTAX species was: {sintax_species}"
            )
        else:
            remarks.append(
                f"VSEARCH species assignment: {vsearch_species}; "
                "SINTAX had no species assignment"
            )
    return "; ".join(remarks)


def merge_asv_taxonomy(
    sintax_row: dict[str, str], vsearch_row: dict[str, str] | None
) -> dict[str, str]:
    """SINTAX is the taxonomic basis, species level comes from VSEARCH only."""

    vsearch_row = vsearch_row or {}
    sintax_species = species_label(sintax_row)
    vsearch_species = species_label(vsearch_row)

    merged: dict[str, str] = {
        "kingdom": (sintax_row.get("Kingdom") or "").strip(),
        "phylum": (sintax_row.get("Phylum") or "").strip(),
        "class": (sintax_row.get("Class") or "").strip(),
        "order": (sintax_row.get("Order") or "").strip(),
        "family": (sintax_row.get("Family") or "").strip(),
        "genus": (sintax_row.get("Genus") or "").strip(),
        "species": vsearch_species,
    }

    merged["identificationRemarks"] = identification_remarks(sintax_species, vsearch_species)
    return merged


def index_rows_by_asv(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    by_id: dict[str, dict[str, str]] = {}
    for row in rows:
        asv_id = (row.get("ASV_ID") or "").strip()
        if asv_id:
            by_id[asv_id] = row
    return by_id


def build_occurrence_table(
    sintax_rows: list[dict[str, str]],
    vsearch_rows: list[dict[str, str]],
    dada2_rows: list[dict[str, str]],
    metadata_by_sample: dict[str, dict[str, str]],
    pipeline_params: dict,
    output_path: Path,
    *,
    clean_prefix: bool,
) -> list[tuple[str, str, str]]:
    """
    Build Darwin Core Occurrence table. Occurrence IDs are constructed as {sample_id}_{asv_id}.
    Non-DNA-derived metadata columns are added to the output.
    """

    _ = pipeline_params

    vsearch_by_asv = index_rows_by_asv(vsearch_rows)
    sintax_by_asv = index_rows_by_asv(sintax_rows)

    # Collect all metadata field names (excluding the various sample ID keys)
    metadata_fields: set[str] = set()
    for meta in metadata_by_sample.values():
        metadata_fields.update(meta.keys())
    metadata_fields.difference_update(set(METADATA_SAMPLE_ID_COLUMNS) | {"", *DNA_DERIVED_DATA_METADATA_FIELDS})

    occurrence_rows: list[dict[str, str]] = []
    occurrence_index: list[tuple[str, str, str]] = []

    # Iterate over DADA2 rows to get sample x ASV abundance
    for dada2_row in dada2_rows:
        asv_id = (dada2_row.get("ASV_ID") or "").strip()
        if not asv_id:
            continue

        sintax_row = sintax_by_asv.get(asv_id)
        if not sintax_row:
            # Skip ASVs without SINTAX taxonomy
            continue
        merged = merge_asv_taxonomy(sintax_row, vsearch_by_asv.get(asv_id))

        for sample_id, raw_value in dada2_row.items():
            if sample_id in {"ASV_ID", "sequence"}:
                continue
            abundance = (raw_value or "").strip()
            if not abundance or abundance == "0":
                continue

            output_sample_id = clean_sample_id(sample_id, clean_prefix=clean_prefix)
            occurrence_id = f"{output_sample_id}_{asv_id}"
            meta_row = metadata_by_sample.get(output_sample_id, metadata_by_sample.get(sample_id, {}))

            row: dict[str, str] = {
                "occurrenceID": occurrence_id,
                "eventID": "",
                "basisOfRecord": "",
                "taxonID": "",
                "organismQuantity": abundance,
                "organismQuantityType": "",
                "sample_id": output_sample_id,
                **{rank: merged.get(rank, "") for rank in OUTPUT_RANK_COLUMNS},
                "identificationRemarks": merged.get("identificationRemarks", ""),
            }

            # Attach all metadata fields for this sample
            for field in metadata_fields:
                row[field] = (meta_row.get(field) or "").strip()

            occurrence_rows.append(row)
            occurrence_index.append((occurrence_id, asv_id, output_sample_id))

    # Extend the fixed Darwin Core fields with all metadata columns
    occurrence_fieldnames = list(OCCURRENCE_FIELDS) + sorted(metadata_fields)
    write_tsv(output_path, tuple(occurrence_fieldnames), occurrence_rows)
    return occurrence_index


def build_dna_derived_data_table(
    dada2_rows: list[dict[str, str]],
    occurrence_index: list[tuple[str, str, str]],
    metadata_by_sample: dict[str, dict[str, str]],
    dna_metadata_fields: tuple[str, ...],
    output_path: Path,
) -> None:
    """
    Build DNA derived data table aligned 1:1 with Occurrence records.

    Each row uses the same occurrenceID as Occurrence, the ASV sequence from the
    DADA2 table, and DNA-derived metadata fields from the matching sample row.
    """

    dada2_by_asv = index_rows_by_asv(dada2_rows)
    dna_rows: list[dict[str, str]] = []

    for occurrence_id, asv_id, sample_id in occurrence_index:
        if not asv_id:
            continue
        sequence = (dada2_by_asv.get(asv_id, {}).get("sequence") or "").strip()
        meta_row = metadata_by_sample.get(sample_id, {})
        row: dict[str, str] = {
            "occurrenceID": occurrence_id,
            "DNA_sequence": sequence,
        }
        for field in dna_metadata_fields:
            row[field] = (meta_row.get(field) or "").strip()
        dna_rows.append(row)

    fieldnames = DNA_DERIVED_DATA_FIELDS + dna_metadata_fields
    write_tsv(output_path, fieldnames, dna_rows)


def build_darwin_core(
    ampliseq_results: Path,
    output_dir: Path,
    metadata_path: Path,
    *,
    clean_prefix: bool = False,
) -> Path:
    """Read ampliseq + worms under output_dir, write Darwin Core tables to publishing/."""

    if not ampliseq_results.is_dir():
        raise FileNotFoundError(f"ampliseq results not found: {ampliseq_results}")
    if not output_dir.is_dir():
        raise FileNotFoundError(f"output directory not found: {output_dir}")

    with find_pipeline_params_json(ampliseq_results).open() as handle:
        pipeline_params = json.load(handle)

    sintax_rows = load_tsv(find_worms_matched(output_dir, "sintax"))
    vsearch_rows = load_tsv(find_worms_matched(output_dir, "vsearch"))
    dada2_rows = load_tsv(find_dada2_table(ampliseq_results))
    metadata_by_sample = load_metadata(metadata_path)
    dna_metadata_fields = metadata_dna_derived_fields(metadata_by_sample)

    dwc_dir = output_dir / "publishing"
    dwc_dir.mkdir(parents=True, exist_ok=True)
    occurrence_index = build_occurrence_table(
        sintax_rows,
        vsearch_rows,
        dada2_rows,
        metadata_by_sample,
        pipeline_params,
        dwc_dir / "occurrence.tsv",
        clean_prefix=clean_prefix,
    )
    build_dna_derived_data_table(
        dada2_rows,
        occurrence_index,
        metadata_by_sample,
        dna_metadata_fields,
        dwc_dir / "dnaderiveddata.tsv",
    )
    return dwc_dir


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ampliseq-results", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="pipeline output directory",
    )
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--clean-prefix",
        action="store_true",
        help="Remove s_ prefix from sample IDs in Darwin Core output",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.metadata.is_file():
        print(f"metadata not found: {args.metadata}", file=sys.stderr)
        return 1
    try:
        build_darwin_core(
            args.ampliseq_results,
            args.output,
            args.metadata,
            clean_prefix=args.clean_prefix,
        )
    except (FileNotFoundError, json.JSONDecodeError, OSError) as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
