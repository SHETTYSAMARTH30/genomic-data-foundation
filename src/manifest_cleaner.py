from __future__ import annotations

import csv
import logging
from datetime import datetime
from dateutil import parser as dateutil_parser

logger = logging.getLogger(__name__)

_PLATFORM_MAP = {
    "novaseq6000": "NovaSeq 6000",
    "novaseqx": "NovaSeq X",
}


def normalize_date(raw: str) -> str | None:
    """Return ISO date string "YYYY-MM-DD" or None for empty/unparseable input."""
    raw = " ".join(raw.strip().split())
    if not raw:
        return None
    # Try unambiguous ISO formats first to avoid dayfirst misinterpretation.
    for fmt in ("%Y-%m-%d", "%Y/%m/%d"):
        try:
            return datetime.strptime(raw, fmt).date().isoformat()
        except ValueError:
            continue
    # Fall back to dateutil with dayfirst=True for DD/MM/YYYY and "Mon DD YYYY".
    try:
        return dateutil_parser.parse(raw, dayfirst=True).date().isoformat()
    except Exception:
        logger.warning("Could not parse date: %r", raw)
        return None


def normalize_purity(raw: str) -> float | None:
    """Return tumor purity normalised to [0.0, 1.0] or None."""
    raw = raw.strip()
    if not raw:
        return None
    try:
        if raw.endswith("%"):
            return float(raw[:-1]) / 100.0
        val = float(raw)
        # Values > 1.0 are assumed to be percentages (e.g. 63 → 0.63).
        return val / 100.0 if val > 1.0 else val
    except ValueError:
        logger.warning("Could not parse tumor_purity: %r", raw)
        return None


def normalize_sex(raw: str) -> str | None:
    """Return "M", "F", or None."""
    v = raw.strip().lower()
    if not v:
        return None
    if v in ("f", "female"):
        return "F"
    if v in ("m", "male"):
        return "M"
    logger.warning("Unknown sex value: %r", raw)
    return None


def normalize_platform(raw: str) -> str:
    """Return canonical platform name: "NovaSeq 6000" or "NovaSeq X"."""
    key = raw.strip().lower().replace(" ", "").replace("-", "")
    return _PLATFORM_MAP.get(key, raw.strip())


def normalize_qc_status(raw: str) -> str | None:
    """Return uppercased QC status or None for blank input."""
    v = raw.strip().upper()
    return v if v else None


def _completeness(row: dict) -> int:
    """Count non-None, non-empty fields — used to select the most complete duplicate."""
    return sum(1 for v in row.values() if v is not None and v != "")


def clean_manifest(csv_path: str) -> list[dict]:
    """
    Read a sample manifest CSV, normalise every field, deduplicate on
    (sample_id, batch_id) keeping the most complete row, and return a list
    of normalised dicts.
    """
    raw_rows: list[dict] = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for raw in csv.DictReader(f):
            raw_rows.append({
                "sample_id":           raw.get("sample_id", "").strip(),
                "batch_id":            raw.get("batch_id", "").strip(),
                "patient_id":          raw.get("patient_id", "").strip() or None,
                "collection_date":     normalize_date(raw.get("collection_date", "")),
                "tissue":              raw.get("tissue", "").strip().lower() or None,
                "diagnosis":           raw.get("diagnosis", "").strip() or None,
                "disease_group":       raw.get("disease_group", "").strip() or None,
                "tumor_purity":        normalize_purity(raw.get("tumor_purity", "")),
                "sex_reported":        normalize_sex(raw.get("sex_reported", "")),
                "sex_inferred":        normalize_sex(raw.get("sex_inferred", "")),
                "sequencing_platform": normalize_platform(raw.get("sequencing_platform", "")),
                "qc_status":           normalize_qc_status(raw.get("qc_status", "")),
                "notes":               raw.get("notes", "").strip() or None,
                # library_prep is present in batch 2 manifests; absent columns
                # are returned as None by dict.get with a default of None.
                "library_prep":        (raw.get("library_prep") or "").strip() or None,
            })

    seen: dict[tuple, dict] = {}
    for row in raw_rows:
        key = (row["sample_id"], row["batch_id"])
        if key in seen:
            logger.warning(
                "Duplicate manifest row for %s/%s — keeping most complete",
                row["sample_id"], row["batch_id"],
            )
            if _completeness(row) > _completeness(seen[key]):
                seen[key] = row
        else:
            seen[key] = row

    return list(seen.values())
