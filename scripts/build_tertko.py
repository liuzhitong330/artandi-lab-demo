#!/usr/bin/env python3
"""Reproduce the processed-data subset and descriptive checks used by the demo.

No rows are deduplicated, imputed, winsorized or omitted. All descriptive
summaries refer to deposited molecule-level measurements, not biological
replicates. A read-length restriction is a stress test, not validated QC.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
from itertools import combinations
import json
import math
from pathlib import Path
import statistics
from urllib.request import Request, urlopen
from zipfile import ZipFile

import openpyxl


SOURCE_URL = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-49007-4/MediaObjects/41467_2024_49007_MOESM10_ESM.zip"
SOURCE_ARCHIVE_SHA256 = "9526cdddcbbf793ba0e17fbc9007c82d2f50073e31b0e946e1de4e97c3e26c9b"
WORKBOOK_SHA256 = "26115393784129d05c68d92138c4a417df1b5cee9a308976065e37cdbcff4d24"
SOURCE_MEMBER = "Source Data File.xlsx"
SHEET = "hesc_tertko"
DOI = "10.1038/s41467-024-49007-4"
COLUMNS = ["row", "day", "telomere_bp", "read_bp", "chromosome", "arm", "direction"]
EXPECTED_HEADERS = ["chromosome", "telomere_length", "read_length", "arm", "direction", "day", "Genotype"]
EXPECTED_COUNTS = {66: 377, 78: 350, 98: 968, 105: 915}
THRESHOLDS = [1000, 1500, 2000, 2500, 3000, 4000]
MIN_READ_BP = [0, 1000, 4000, 8000]


def sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def quantile(values: list[float], p: float) -> float | None:
    """Linear interpolation at (n-1)*p (R type 7 / NumPy default)."""
    if not values:
        return None
    values = sorted(values)
    pos = (len(values) - 1) * p
    lower = math.floor(pos)
    upper = math.ceil(pos)
    return values[lower] + (values[upper] - values[lower]) * (pos - lower)


def describe(rows: list[list], day: int, threshold: int = 2000, min_read: int = 0) -> dict:
    selected = [r for r in rows if r[1] == day and r[3] >= min_read]
    values = [r[2] for r in selected]
    n = len(values)
    short_n = sum(v < threshold for v in values)
    return {
        "day": day,
        "n": n,
        "mean": statistics.mean(values) if n else None,
        "median": statistics.median(values) if n else None,
        "q1": quantile(values, 0.25),
        "q3": quantile(values, 0.75),
        f"short_lt_{threshold}_n": short_n,
        f"short_lt_{threshold}_fraction": short_n / n if n else None,
    }


def contrast(rows: list[list], earlier: int, later: int, threshold: int, min_read: int) -> dict:
    a = describe(rows, earlier, threshold, min_read)
    b = describe(rows, later, threshold, min_read)
    key = f"short_lt_{threshold}_fraction"
    result = {
        "earlier_day": earlier,
        "later_day": later,
        "threshold_bp": threshold,
        "minimum_read_bp": min_read,
        "earlier_n": a["n"],
        "later_n": b["n"],
        "earlier_retained_fraction": a["n"] / sum(r[1] == earlier for r in rows),
        "later_retained_fraction": b["n"] / sum(r[1] == later for r in rows),
        "mean_change_bp": b["mean"] - a["mean"] if a["n"] and b["n"] else None,
        "median_change_bp": b["median"] - a["median"] if a["n"] and b["n"] else None,
        "short_fraction_change_pp": 100 * (b[key] - a[key]) if a["n"] and b["n"] else None,
    }
    result["mean_median_opposite_signs"] = (
        result["mean_change_bp"] * result["median_change_bp"] < 0
        if a["n"] and b["n"] else None
    )
    return result


def get_source_workbook(local_workbook: Path | None) -> bytes:
    if local_workbook:
        content = local_workbook.read_bytes()
    else:
        request = Request(SOURCE_URL, headers={"User-Agent": "artandi-demo-reproducibility/1.0"})
        with urlopen(request, timeout=60) as response:
            archive = response.read(50_000_001)
        if len(archive) > 50_000_000:
            raise ValueError("Archive exceeds the 50 MB safety bound")
        if sha256(archive) != SOURCE_ARCHIVE_SHA256:
            raise ValueError("Source archive SHA256 mismatch; stop rather than silently updating data")
        # Read only the exact expected member into memory. Never extract a
        # remote archive to disk or trust archive-supplied destination paths.
        with ZipFile(io.BytesIO(archive)) as bundle:
            candidates = [n for n in bundle.namelist() if n == SOURCE_MEMBER]
            if len(candidates) != 1:
                raise ValueError("Expected workbook member is missing or ambiguous")
            info = bundle.getinfo(SOURCE_MEMBER)
            if info.file_size > 20_000_000:
                raise ValueError("Workbook exceeds the 20 MB safety bound")
            content = bundle.read(SOURCE_MEMBER)
    if sha256(content) != WORKBOOK_SHA256:
        raise ValueError("Workbook SHA256 mismatch; stop rather than silently updating data")
    return content


def parse_workbook(content: bytes) -> tuple[list[list], dict]:
    workbook = openpyxl.load_workbook(io.BytesIO(content), read_only=True, data_only=True)
    sheet = workbook[SHEET]
    iterator = sheet.iter_rows(values_only=True)
    title = next(iterator)
    headers = list(next(iterator))
    if headers != EXPECTED_HEADERS:
        raise ValueError(f"Unexpected schema: {headers!r}")
    rows = []
    originals = []
    invalid = []
    telomere_exceeds_read = []
    for physical_row, values in enumerate(iterator, start=3):
        if len(values) != 7:
            raise ValueError(f"Unexpected column count at worksheet row {physical_row}")
        chromosome, telomere, read_length, arm, direction, day, genotype = values
        originals.append(tuple(values))
        if genotype != "TERT KO":
            raise ValueError(f"Unexpected genotype at worksheet row {physical_row}: {genotype!r}")
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0
               for v in (telomere, read_length)):
            invalid.append(physical_row)
        if telomere > read_length:
            telomere_exceeds_read.append(physical_row)
        if day not in EXPECTED_COUNTS:
            raise ValueError(f"Unexpected passage day at worksheet row {physical_row}: {day!r}")
        rows.append([physical_row, day, telomere, read_length, chromosome, arm, direction])
    workbook.close()
    if invalid:
        raise ValueError(f"Invalid lengths at rows {invalid}; no automatic row dropping is permitted")
    if len(rows) != 2610 or dict(Counter(r[1] for r in rows)) != EXPECTED_COUNTS:
        raise ValueError("Source row count differs from the pinned dataset")
    original_counts = Counter(originals)
    repeated_groups = [[rows[i][0] for i, row in enumerate(originals) if row == item]
                       for item, count in original_counts.items() if count > 1]
    audit = {
        "worksheet_title": title[0],
        "source_headers": headers,
        "worksheet_first_data_row": 3,
        "worksheet_last_data_row": rows[-1][0],
        "total_rows": len(rows),
        "rows_by_day": dict(sorted(Counter(r[1] for r in rows).items())),
        "nonpositive_or_nonfinite_length_rows": invalid,
        "telomere_exceeds_read_length_rows": telomere_exceeds_read,
        "duplicate_visible_rows": len(originals) - len(original_counts),
        "repeated_visible_row_groups": repeated_groups,
        "duplicate_policy": "Preserved. Identical seven-field values cannot establish duplicate molecules because no source molecule ID is provided.",
        "minimum_read_bp": min(r[3] for r in rows),
        "maximum_read_bp": max(r[3] for r in rows),
        "minimum_telomere_bp": min(r[2] for r in rows),
        "maximum_telomere_bp": max(r[2] for r in rows),
        "source_day_discrepancy": "The workbook labels the final sample day 105, as does the article narrative. The Figure 2 caption says 108. Values are retained as deposited, without silently harmonizing labels.",
    }
    if telomere_exceeds_read:
        # Preserve a failed check rather than silently correcting or excluding
        # the affected records. This pinned workbook currently has no failures.
        print(f"WARNING: telomere length exceeds read length in rows {telomere_exceeds_read}")
    return rows, audit


def build(local_workbook: Path | None, output: Path) -> tuple[dict, dict]:
    content = get_source_workbook(local_workbook)
    rows, audit = parse_workbook(content)
    days = sorted(EXPECTED_COUNTS)
    source = {
        "url": SOURCE_URL,
        "doi": DOI,
        "sheet": SHEET,
        "sha256": WORKBOOK_SHA256,
        "archive_sha256": SOURCE_ARCHIVE_SHA256,
        "archive_member": SOURCE_MEMBER,
        "license": "CC BY 4.0",
        "license_url": "https://creativecommons.org/licenses/by/4.0/",
        "attribution": "Sanchez et al. (2024), Nature Communications, doi:10.1038/s41467-024-49007-4. Processed source data, transformed into a row-preserving subset and descriptive summaries by Cathy Liu with AI assistance.",
        "notes": [
            "Processed telomere measurements from the published source workbook, not raw nanopore signal.",
            "Row identifiers are physical Excel worksheet rows, not original molecule IDs.",
            "All 2610 rows are retained, including 25 repeated visible records; no duplicate-molecule inference is possible.",
            "Four deposited passage distributions; independent biological replicate and clone IDs are not supplied in this sheet.",
            "Final sample uses deposited day 105; the Figure 2 caption says 108.",
            "A short-tail threshold is exploratory, not a validated biological or clinical cutoff. Read-length restrictions are sensitivity tests, not recommended QC.",
            "Quartiles use linear interpolation at (n-1)*p. Short counts use strict <; read retention uses >=.",
        ],
    }
    payload = {"source": source, "columns": COLUMNS, "rows": rows, "days": days,
               "summary": [describe(rows, day) for day in days],
               "duplicate_visible_rows": audit["duplicate_visible_rows"]}
    checks = {"source": source, "audit": audit, "thresholds_bp": THRESHOLDS,
              "minimum_read_bp": MIN_READ_BP,
              "contrasts": [contrast(rows, a, b, threshold, min_read)
                            for a, b in combinations(days, 2)
                            for threshold in THRESHOLDS for min_read in MIN_READ_BP]}
    output.mkdir(parents=True, exist_ok=True)
    (output / "tertko.json").write_text(json.dumps(payload, separators=(",", ":"), allow_nan=False) + "\n", encoding="utf-8")
    (output / "analysis_checks.json").write_text(json.dumps(checks, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (output / "tertko.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)
    return payload, checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, help="Use a local copy of the pinned workbook instead of downloading it")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    args = parser.parse_args()
    payload, checks = build(args.workbook, args.output)
    print(json.dumps({"summary": payload["summary"], "audit": checks["audit"],
                      "contrast_cases": len(checks["contrasts"])}, indent=2))


if __name__ == "__main__":
    main()
