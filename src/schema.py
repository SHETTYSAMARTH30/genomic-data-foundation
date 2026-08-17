from datetime import datetime, timezone
import duckdb

SAMPLES_DDL = """
CREATE TABLE IF NOT EXISTS samples (
    sample_id           VARCHAR NOT NULL,
    batch_id            VARCHAR NOT NULL,
    patient_id          VARCHAR,
    collection_date     DATE,
    tissue              VARCHAR,
    diagnosis           VARCHAR,
    disease_group       VARCHAR,
    tumor_purity        FLOAT,
    sex_reported        VARCHAR,
    sex_inferred        VARCHAR,
    sequencing_platform VARCHAR,
    qc_status           VARCHAR,
    notes               VARCHAR,
    library_prep        VARCHAR,
    ingested_at         TIMESTAMP,
    PRIMARY KEY (sample_id, batch_id)
)
"""

VARIANTS_DDL = """
CREATE TABLE IF NOT EXISTS variants (
    sample_id        VARCHAR NOT NULL,
    batch_id         VARCHAR NOT NULL,
    chrom            VARCHAR NOT NULL,
    pos              BIGINT  NOT NULL,
    vcf_id           VARCHAR,
    ref              VARCHAR NOT NULL,
    alt              VARCHAR NOT NULL,
    qual             FLOAT,
    filter           VARCHAR,
    caller           VARCHAR,
    max_pop_af       FLOAT,
    hotspot          BOOLEAN,
    ccf              FLOAT,
    gnomad_af_popmax FLOAT,
    gene             VARCHAR,
    consequence      VARCHAR,
    impact           VARCHAR,
    hgvsc            VARCHAR,
    hgvsp            VARCHAR,
    exon             VARCHAR,
    mane_select      VARCHAR,
    gt               VARCHAR,
    ad_ref           INTEGER,
    ad_alt           INTEGER,
    dp               INTEGER,
    vaf              FLOAT,
    ingested_at      TIMESTAMP,
    PRIMARY KEY (sample_id, batch_id, chrom, pos, ref, alt)
)
"""

_SAMPLE_COLS = [
    "sample_id", "batch_id", "patient_id", "collection_date", "tissue",
    "diagnosis", "disease_group", "tumor_purity", "sex_reported", "sex_inferred",
    "sequencing_platform", "qc_status", "notes", "library_prep", "ingested_at",
]

_VARIANT_COLS = [
    "sample_id", "batch_id", "chrom", "pos", "vcf_id", "ref", "alt", "qual",
    "filter", "caller", "max_pop_af", "hotspot", "ccf", "gnomad_af_popmax",
    "gene", "consequence", "impact", "hgvsc", "hgvsp", "exon", "mane_select",
    "gt", "ad_ref", "ad_alt", "dp", "vaf", "ingested_at",
]


def init_db(db_path: str) -> duckdb.DuckDBPyConnection:
    conn = duckdb.connect(db_path)
    conn.execute(SAMPLES_DDL)
    conn.execute(VARIANTS_DDL)
    return conn


def upsert_samples(conn: duckdb.DuckDBPyConnection, rows: list[dict]) -> int:
    if not rows:
        return 0
    now = datetime.now(timezone.utc)
    placeholders = ", ".join("?" * len(_SAMPLE_COLS))
    sql = f"""
        INSERT INTO samples ({', '.join(_SAMPLE_COLS)})
        VALUES ({placeholders})
        ON CONFLICT (sample_id, batch_id) DO NOTHING
    """
    for row in rows:
        values = [row.get(c) for c in _SAMPLE_COLS[:-1]] + [now]
        conn.execute(sql, values)
    return len(rows)


def upsert_variants(conn: duckdb.DuckDBPyConnection, rows: list[dict]) -> int:
    if not rows:
        return 0
    now = datetime.now(timezone.utc)
    placeholders = ", ".join("?" * len(_VARIANT_COLS))
    sql = f"""
        INSERT INTO variants ({', '.join(_VARIANT_COLS)})
        VALUES ({placeholders})
        ON CONFLICT (sample_id, batch_id, chrom, pos, ref, alt) DO NOTHING
    """
    for row in rows:
        values = [row.get(c) for c in _VARIANT_COLS[:-1]] + [now]
        conn.execute(sql, values)
    return len(rows)
