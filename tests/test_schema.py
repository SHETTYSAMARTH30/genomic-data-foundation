import pytest
import duckdb
from schema import init_db, upsert_samples, upsert_variants


def _sample_row(**overrides):
    base = {
        "sample_id": "S-0001", "batch_id": "batch_test",
        "patient_id": "P-0001", "collection_date": "2026-01-05",
        "tissue": "lung", "diagnosis": "Lung adenocarcinoma",
        "disease_group": "Thoracic", "tumor_purity": 0.32,
        "sex_reported": "F", "sex_inferred": "F",
        "sequencing_platform": "NovaSeq X", "qc_status": "PASS",
        "notes": None, "library_prep": None,
    }
    return {**base, **overrides}


def _variant_row(**overrides):
    base = {
        "sample_id": "S-0001", "batch_id": "batch_test",
        "chrom": "chr1", "pos": 1000, "vcf_id": None,
        "ref": "A", "alt": "G", "qual": 100.0,
        "filter": "PASS", "caller": "mutect2",
        "max_pop_af": 0.0, "hotspot": True,
        "ccf": None, "gnomad_af_popmax": None,
        "gene": "GENE1", "consequence": "missense_variant",
        "impact": "MODERATE", "hgvsc": "c.100A>G", "hgvsp": "p.Arg34Gly",
        "exon": "1/5", "mane_select": None,
        "gt": "0/1", "ad_ref": 50, "ad_alt": 30, "dp": 80, "vaf": 0.375,
    }
    return {**base, **overrides}


def test_init_db_creates_tables(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    tables = {t[0] for t in conn.execute("SHOW TABLES").fetchall()}
    assert "samples" in tables
    assert "variants" in tables
    conn.close()


def test_upsert_samples_inserts_row(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    count = upsert_samples(conn, [_sample_row()])
    assert count == 1
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 1
    conn.close()


def test_upsert_variants_inserts_row(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    count = upsert_variants(conn, [_variant_row()])
    assert count == 1
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_upsert_samples_idempotent(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    upsert_samples(conn, [_sample_row()])
    upsert_samples(conn, [_sample_row()])  # same PK — must be ignored
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 1
    conn.close()


def test_upsert_variants_idempotent(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    upsert_variants(conn, [_variant_row()])
    upsert_variants(conn, [_variant_row()])
    assert conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 1
    conn.close()


def test_upsert_samples_different_batches(tmp_path):
    conn = init_db(str(tmp_path / "test.duckdb"))
    upsert_samples(conn, [_sample_row(batch_id="batch_a")])
    upsert_samples(conn, [_sample_row(batch_id="batch_b")])
    assert conn.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == 2
    conn.close()
