#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import os
import sys
from pathlib import Path
from aphiasync.sqlite import match as worms_sqlite_match


# Most-specific rank first (used by finest_taxon).
FINEST_RANK_COLUMNS = (
    "Species",
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum",
    "Kingdom",
)
WORMS_OUTPUT_COLUMNS = ("scientificName", "scientificNameID", "taxonRank")
WORMS_LSID_PREFIX = "urn:lsid:marinespecies.org:taxname:"
METHOD_PREFIX = {
    "sintax": "ASV_tax_sintax",
    "vsearch": "ASV_tax_vsearch_lca",
}
METHOD_SUBDIR = {
    "sintax": "sintax",
    "vsearch": "vsearch_lca",
}
PLACEHOLDER_NAMES = frozenset({"na", "unassigned", ""})


def load_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, [dict(row) for row in reader]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def find_taxonomy_tsv(input_path: Path, method: str) -> Path:
    """
    Resolve ampliseq taxonomy TSV from a file path, method folder, or ampliseq outdir.
    Picks the single ASV_tax_*.tsv matching the method prefix.
    """

    if input_path.is_file():
        return input_path

    if not input_path.is_dir():
        raise FileNotFoundError(f"input not found: {input_path}")

    prefix = METHOD_PREFIX[method]
    subdir = METHOD_SUBDIR[method]
    search_dirs: list[Path] = [input_path]
    method_dir = input_path / subdir
    if method_dir.is_dir() and method_dir not in search_dirs:
        search_dirs.insert(0, method_dir)

    for directory in search_dirs:
        matches = sorted(
            path
            for path in directory.iterdir()
            if path.is_file()
            and path.name.startswith(prefix)
            and path.suffix == ".tsv"
            and ".raw." not in path.name
        )
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            names = ", ".join(path.name for path in matches)
            raise FileNotFoundError(f"multiple {prefix}*.tsv files in {directory}: {names}")

    raise FileNotFoundError(f"no {prefix}*.tsv (non-.raw) found under {input_path}")


def aphiaid_to_lsid(aphiaid: int | str | None) -> str:
    if aphiaid is None or aphiaid == "":
        return ""
    return f"{WORMS_LSID_PREFIX}{int(aphiaid)}"


def worms_match_columns(worms: dict[str, str | int | None]) -> tuple[str, str]:
    """Return (scientificName, scientificNameID) for output from match results."""

    matched_aphiaid = worms.get("matched_aphiaid")
    if matched_aphiaid is not None:
        matched_name = worms.get("matched_name") or worms.get("canonical") or ""
        return matched_name, aphiaid_to_lsid(matched_aphiaid)
    return worms.get("canonical") or "", ""


def finest_taxon(row: dict[str, str], ranks: tuple[str, ...] = FINEST_RANK_COLUMNS) -> tuple[str, str]:
    """Return (scientific_name, taxon_rank) for the most specific non-empty rank."""

    for rank in ranks:
        # Accept both Darwin Core lowercase ranks and upstream title-case ranks.
        name = (row.get(rank) or row.get(rank.lower()) or row.get(rank.upper()) or "").strip()
        if name and name.lower() not in PLACEHOLDER_NAMES:
            return name, rank
    return "", ""


def match_distinct_taxa(
    taxa: set[tuple[str, str]],
    worms_db: Path,
) -> dict[tuple[str, str], dict[str, str]]:
    """
    Match each distinct scientific name once via aphiasync SQLite lookup.
    Rank is kept on each cache key for downstream use; matching uses names only.
    """

    if not worms_db.is_file():
        raise FileNotFoundError(f"WoRMS SQLite database not found: {worms_db}")

    distinct_names = sorted({name for name, _rank in taxa})
    if not distinct_names:
        return {}

    os.environ["WORMS_DB_PATH"] = str(worms_db)
    name_results = worms_sqlite_match(distinct_names)

    cache: dict[tuple[str, str], dict[str, str]] = {}
    for scientific_name, taxon_rank in taxa:
        hit = name_results.get(scientific_name, {})
        matched_aphiaid = hit.get("aphiaid")
        canonical = hit.get("canonical")
        record = hit.get("record") or {}
        matched_name = record.get("valid_name") or record.get("scientificname") or ""
        cache[(scientific_name, taxon_rank)] = {
            "matched_aphiaid": matched_aphiaid,
            "matched_name": matched_name,
            "canonical": canonical or "",
        }
    return cache


def enrich_rows_with_worms(
    rows: list[dict[str, str]],
    worms_db: Path,
    ranks: tuple[str, ...] = FINEST_RANK_COLUMNS,
) -> list[dict[str, str]]:
    """Add or overwrite scientificName and scientificNameID from finest rank + WoRMS."""

    distinct_taxa: set[tuple[str, str]] = set()
    for row in rows:
        scientific_name, taxon_rank = finest_taxon(row, ranks=ranks)
        if scientific_name:
            distinct_taxa.add((scientific_name, taxon_rank))

    worms_cache = match_distinct_taxa(distinct_taxa, worms_db=worms_db)

    enriched: list[dict[str, str]] = []
    for row in rows:
        out_row = dict(row)
        scientific_name, taxon_rank = finest_taxon(row, ranks=ranks)
        if not scientific_name:
            out_row["scientificName"] = ""
            out_row["scientificNameID"] = ""
            out_row["taxonRank"] = ""
        else:
            scientific_name_out, scientific_name_id = worms_match_columns(
                worms_cache[(scientific_name, taxon_rank)]
            )
            out_row["scientificName"] = scientific_name_out
            out_row["scientificNameID"] = scientific_name_id
            out_row["taxonRank"] = taxon_rank.lower()
        enriched.append(out_row)
    return enriched


def output_fieldnames(input_fieldnames: list[str]) -> list[str]:
    fieldnames = list(input_fieldnames)
    for column in WORMS_OUTPUT_COLUMNS:
        if column not in fieldnames:
            fieldnames.append(column)
    return fieldnames


def process_tsv_table(
    input_path: Path,
    output_path: Path,
    worms_db: Path,
    ranks: tuple[str, ...] = FINEST_RANK_COLUMNS,
) -> None:
    """Read a TSV with taxonomic rank columns; write it with WoRMS match columns added."""

    input_fieldnames, rows = load_tsv(input_path)
    enriched = enrich_rows_with_worms(rows, worms_db=worms_db, ranks=ranks)
    write_tsv(output_path, output_fieldnames(input_fieldnames), enriched)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--method",
        choices=("sintax", "vsearch"),
        help="ampliseq taxonomy method (required when --input is a directory)",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="TSV with rank columns, or ampliseq outdir / method folder for --method",
    )
    parser.add_argument("--output", type=Path, required=True, help="WoRMS-enriched TSV")
    parser.add_argument(
        "--worms-db",
        type=Path,
        default=os.environ.get("WORMS_DB_PATH"),
        help="WoRMS SQLite database",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.worms_db:
        print("Missing WoRMS database: pass --worms-db", file=sys.stderr)
        return 1
    args.worms_db = Path(args.worms_db)
    input_path = Path(args.input)
    try:
        if input_path.is_file():
            table_path = input_path
        else:
            if not args.method:
                print("--method is required when --input is a directory", file=sys.stderr)
                return 1
            table_path = find_taxonomy_tsv(input_path, args.method)
        process_tsv_table(table_path, args.output, worms_db=args.worms_db)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
