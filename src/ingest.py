# src/ingest.py
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from schema import init_db, upsert_samples, upsert_variants
from manifest_cleaner import clean_manifest
from vcf_parser import parse_vcf

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

DEFAULT_DB = "curated/genomics.duckdb"


def ingest_batch(batch_path: str, db_path: str = DEFAULT_DB) -> None:
    batch_dir = Path(batch_path)
    if not batch_dir.is_dir():
        logger.error("Batch directory not found: %s", batch_path)
        sys.exit(1)

    manifest_path = batch_dir / "sample_manifest.csv"
    if not manifest_path.exists():
        logger.error("No manifest at %s", manifest_path)
        sys.exit(1)

    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = init_db(db_path)

    logger.info("Cleaning manifest: %s", manifest_path)
    sample_rows = clean_manifest(str(manifest_path))
    upsert_samples(conn, sample_rows)
    logger.info("  %d sample rows loaded", len(sample_rows))

    manifest_ids = {r["sample_id"] for r in sample_rows}
    vcf_files = sorted(batch_dir.glob("*.somatic.vcf"))
    logger.info("Found %d VCF file(s)", len(vcf_files))

    total = 0
    skipped = 0
    for vcf_path in vcf_files:
        sample_id = vcf_path.name.split(".")[0]
        batch_id = batch_dir.name

        if sample_id not in manifest_ids:
            logger.warning("  %s: no manifest entry — skipping", vcf_path.name)
            skipped += 1
            continue

        try:
            variants = parse_vcf(str(vcf_path), sample_id, batch_id)
            upsert_variants(conn, variants)
            total += len(variants)
            logger.info("  %s: %d variants", vcf_path.name, len(variants))
        except Exception as exc:
            logger.error("  %s: FAILED — %s", vcf_path.name, exc)
            skipped += 1

    vcf_sample_ids = {vcf_path.stem.split(".")[0] for vcf_path in vcf_files}
    for sid in manifest_ids - vcf_sample_ids:
        logger.warning("sample %s has manifest entry but no VCF file", sid)

    conn.close()
    logger.info("Done. %d variants loaded, %d file(s) skipped.", total, skipped)


def main() -> None:
    parser = argparse.ArgumentParser(description="Ingest a genomic batch into DuckDB")
    parser.add_argument("--batch", required=True,
                        help="Batch directory (e.g. data/batch_2026_01)")
    parser.add_argument("--db", default=DEFAULT_DB,
                        help=f"DuckDB output path (default: {DEFAULT_DB})")
    args = parser.parse_args()
    ingest_batch(args.batch, args.db)


if __name__ == "__main__":
    main()
