# tests/test_integration.py
from __future__ import annotations

import pytest
import duckdb
from pathlib import Path
from ingest import ingest_batch


@pytest.fixture
def db_path(tmp_path):
    return str(tmp_path / "test.duckdb")


def test_batch1_loads_samples(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count = conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    conn.close()
    # Batch 1 has 16 manifest rows; S-0011 is a duplicate → 15 unique samples
    assert count == 15


def test_batch1_loads_variants(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()
    # 13 VCF files loaded (14 present minus S-0008 which has no VCF)
    # Each file has ~50+ variants; total well above 0
    assert count > 100


def test_batch1_idempotent(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count_after_first = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count_after_second = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()
    assert count_after_first == count_after_second


def test_batch2_loads_additively(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count_b1 = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()

    ingest_batch("data/batch_2026_02", db_path)
    conn = duckdb.connect(db_path)
    count_b2 = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    conn.close()

    assert count_b2 > count_b1


def test_batch1_queries_unchanged_after_batch2(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    count_before = conn.execute(
        "SELECT COUNT(*) FROM variants WHERE batch_id = 'batch_2026_01'"
    ).fetchone()[0]
    conn.close()

    ingest_batch("data/batch_2026_02", db_path)
    conn = duckdb.connect(db_path)
    count_after = conn.execute(
        "SELECT COUNT(*) FROM variants WHERE batch_id = 'batch_2026_01'"
    ).fetchone()[0]
    conn.close()

    assert count_before == count_after


def test_all_vaf_in_range(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    ingest_batch("data/batch_2026_02", db_path)
    conn = duckdb.connect(db_path)
    bad = conn.execute(
        "SELECT COUNT(*) FROM variants WHERE vaf IS NOT NULL AND (vaf < 0 OR vaf > 1)"
    ).fetchone()[0]
    conn.close()
    assert bad == 0, f"{bad} variants have VAF outside [0, 1]"


def test_all_purity_in_range(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    conn = duckdb.connect(db_path)
    bad = conn.execute(
        "SELECT COUNT(*) FROM samples WHERE tumor_purity IS NOT NULL AND (tumor_purity < 0 OR tumor_purity > 1)"
    ).fetchone()[0]
    conn.close()
    assert bad == 0


def test_batch2_samples_have_library_prep(db_path):
    ingest_batch("data/batch_2026_02", db_path)
    conn = duckdb.connect(db_path)
    missing = conn.execute(
        "SELECT COUNT(*) FROM samples WHERE batch_id = 'batch_2026_02' AND library_prep IS NULL"
    ).fetchone()[0]
    conn.close()
    assert missing == 0


def test_kras_query_runs(db_path):
    ingest_batch("data/batch_2026_01", db_path)
    ingest_batch("data/batch_2026_02", db_path)
    conn = duckdb.connect(db_path)
    rows = conn.execute("""
        SELECT s.sample_id, s.tissue, v.hgvsp, v.vaf
        FROM variants v
        JOIN samples s USING (sample_id)
        WHERE v.gene = 'KRAS' AND v.filter = 'PASS'
        ORDER BY v.vaf DESC
    """).fetchall()
    conn.close()
    assert len(rows) > 0
