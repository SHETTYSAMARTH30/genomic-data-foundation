from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

_CSQ_FIELDS = ["gene", "consequence", "impact", "hgvsc", "hgvsp", "exon", "mane_select"]


def parse_info(info_str: str) -> dict:
    """Parse a VCF INFO field string into a dict.

    KEY=VALUE pairs are stored as str values.
    FLAG fields (no '=') are stored as True.
    """
    result = {}
    for part in info_str.split(";"):
        if "=" in part:
            key, _, val = part.partition("=")
            result[key.strip()] = val
        elif part.strip():
            result[part.strip()] = True
    return result


def parse_csq(csq_str: str) -> dict:
    """Parse a pipe-delimited CSQ annotation string.

    Batch 1 has 6 fields; batch 2 has 7 (the 7th is mane_select).
    Empty fields become None.
    """
    parts = csq_str.split("|")
    result = {}
    for i, field in enumerate(_CSQ_FIELDS):
        val = parts[i] if i < len(parts) else ""
        result[field] = val if val else None
    return result


def normalize_vaf(vaf: float) -> float:
    """Normalize VAF to the 0.0–1.0 range.

    Batch 2 VCFs encode VAF as a percentage (e.g. 56.14).  Any value > 1.0
    is divided by 100 to produce a fraction.
    """
    return vaf / 100.0 if vaf > 1.0 else vaf


def _parse_ad(ad_str: str) -> tuple[int | None, int | None]:
    """Parse AD field 'ref,alt' into (ref_count, alt_count)."""
    try:
        ref, alt = ad_str.split(",", 1)
        return int(ref), int(alt)
    except (ValueError, AttributeError):
        return None, None


def _safe_int(val: str | None) -> int | None:
    """Safely parse an integer from a string, handling None and '.' values."""
    if val is None or val == ".":
        return None
    try:
        return int(val)
    except ValueError:
        return None


def _float_from_info(info: dict, key: str) -> float | None:
    """Extract and parse a float value from the INFO dict.

    Returns None if the key is not present, is a flag (True), or cannot be parsed.
    """
    v = info.get(key)
    if v is None or v is True:
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def parse_vcf(vcf_path: str, sample_id: str, batch_id: str) -> list[dict]:
    """Parse a VCF file and return a list of variant dicts.

    Malformed data rows (not exactly 10 tab-separated columns) are skipped
    with a WARNING log message.

    Schema keys returned per variant:
        sample_id, batch_id, chrom, pos, vcf_id, ref, alt, qual, filter,
        caller, max_pop_af, hotspot, ccf, gnomad_af_popmax,
        gene, consequence, impact, hgvsc, hgvsp, exon, mane_select,
        gt, ad_ref, ad_alt, dp, vaf
    """
    variants: list[dict] = []
    path = Path(vcf_path)

    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.rstrip("\n")
            if line.startswith("#") or not line:
                continue

            parts = line.split("\t")
            if len(parts) != 10:
                logger.warning(
                    "%s line %d: expected 10 columns, got %d — skipping",
                    path.name,
                    lineno,
                    len(parts),
                )
                continue

            chrom, pos, vcf_id, ref, alt, qual, filter_, info_str, fmt_str, sample_str = parts

            info = parse_info(info_str)
            csq = parse_csq(info.get("CSQ", ""))

            fmt_keys = fmt_str.split(":")
            fmt_vals = sample_str.split(":")

            if len(fmt_keys) != len(fmt_vals):
                logger.warning(
                    "%s line %d: FORMAT has %d fields but sample column has %d values — skipping",
                    path.name, lineno, len(fmt_keys), len(fmt_vals),
                )
                continue

            fmt = dict(zip(fmt_keys, fmt_vals))

            ad_ref, ad_alt = _parse_ad(fmt.get("AD", ""))

            try:
                vaf = normalize_vaf(float(fmt["VAF"]))
            except (KeyError, ValueError):
                vaf = None

            try:
                qual_val: float | None = float(qual) if qual != "." else None
            except ValueError:
                qual_val = None

            variants.append({
                "sample_id":        sample_id,
                "batch_id":         batch_id,
                "chrom":            chrom,
                "pos":              int(pos),
                "vcf_id":           vcf_id if vcf_id != "." else None,
                "ref":              ref,
                "alt":              alt,
                "qual":             qual_val,
                "filter":           filter_,
                "caller":           info.get("CALLER"),
                "max_pop_af":       _float_from_info(info, "MAX_POP_AF"),
                "hotspot":          bool(info.get("HOTSPOT", False)),
                "ccf":              _float_from_info(info, "CCF"),
                "gnomad_af_popmax": _float_from_info(info, "GNOMAD_AF_POPMAX"),
                **csq,
                "gt":     fmt.get("GT"),
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                "dp":     _safe_int(fmt.get("DP")),
                "vaf":    vaf,
            })

    return variants
