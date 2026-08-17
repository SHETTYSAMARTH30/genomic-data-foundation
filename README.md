# Genomic Data Foundation — First Slice

## Quick Start

```bash
pip install -r requirements.txt

python src/ingest.py --batch data/batch_2026_01
python src/ingest.py --batch data/batch_2026_02
```

Both commands are safe to re-run; the pipeline is idempotent.

## Running Tests

```bash
pytest tests/ -v
```

## Running the Example Queries

```bash
python -c "
import duckdb
conn = duckdb.connect('curated/genomics.duckdb')
print(conn.sql(open('queries/01_variants_in_gene_by_tissue.sql').read()).df().to_string())
"
```

Or open an interactive session:

```bash
python -c "import duckdb; conn = duckdb.connect('curated/genomics.duckdb')"
```

---

## Schema

**`samples`** — grain: one row per `(sample_id, batch_id)`

Holds clinical and administrative metadata for each sample. All values from the manifest are normalized before loading (see Cleaning section below).

**`variants`** — grain: one row per `(sample_id, chrom, pos, ref, alt)`

One row per mutation per sample. Consequence annotations (gene, HGVS notation, impact tier) are inlined from the VCF CSQ field. `vaf` is always stored in the 0.0–1.0 fraction range regardless of how the source file expressed it.

---

## Key Decisions and Trade-offs

**DuckDB over SQLite.** Columnar storage makes cohort-level `GROUP BY` queries fast. DuckDB's single-file format is just as simple to share as SQLite. For 1000× scale, swap the DuckDB file for a cloud warehouse (BigQuery, Redshift, Snowflake) — the SQL is portable.

**Custom VCF parser over cyvcf2/pysam.** The files are standard VCFv4.2 with no exotic features. A 60-line custom reader has zero compiled dependencies and is auditable by a reviewer in five minutes. cyvcf2 would matter at scale (millions of variants) but is overkill here.

**Two tables (moderately normalized).** A single flat table would duplicate sample metadata for every variant row. Full normalization (separate callers, consequences, filters tables) adds joins without research value. Two tables with `sample_id` as the join key is the right balance.

**Additive schema migration.** New columns from batch_2026_02 (`library_prep`, `ccf`, `gnomad_af_popmax`, `mane_select`) are included in the initial DDL as nullable columns. Batch_2026_01 rows carry NULL for these fields. No DDL changes are required when batch_2026_02 arrives, and no batch_2026_01 queries break.

**Idempotency via `ON CONFLICT DO NOTHING`.** Running the pipeline twice produces identical results. Primary keys are `(sample_id, batch_id)` for samples and `(sample_id, chrom, pos, ref, alt)` for variants.

---

## Data Cleaning

### Manifest

| Issue | Fix |
|---|---|
| 6 different date formats | Cascade through ISO formats then `dateutil.parser.parse(dayfirst=True)` |
| Mixed tissue case/whitespace | `strip().lower()` |
| Purity as `0.32` or `63%` | Strip `%`, divide by 100 if >1.0 |
| Sex as `f`, `female`, `Female`, `m`, etc. | Map to `F` / `M` |
| Platform spelling variants | Map to canonical `NovaSeq 6000` / `NovaSeq X` |
| QC status mixed case | `.upper()` |
| S-0011 appears twice | Keep most complete row; log WARNING |
| S-0008 in manifest, no VCF | Load manifest row; log that no variants were loaded |

### VCFs

| Issue | Fix |
|---|---|
| Batch_2026_02 VAF in 0–100 range | If `vaf > 1.0`, divide by 100 |
| `HOTSPOT` is a FLAG (no `=value`) | Detect by absence of `=`; map to `True` |
| Batch_2026_02 CSQ has 7 fields vs 6 | Parse by position; 7th field → `mane_select`, else `None` |
| S-0020 line 101 truncated (8 columns) | Skip row with WARNING; continue loading rest of file |
| `.` for missing `ID` or `QUAL` | Map to `NULL` |

---

## What I Deliberately Did Not Do

- **No Parquet layer** — DuckDB on a single file is sufficient for this dataset and a clean fit for the submission format.
- **No workflow orchestration** (Airflow, Prefect) — a CLI script is simpler for a reviewer to understand and run.
- **No Dockerfile** — dependencies are minimal enough that `pip install -r requirements.txt` is all that's needed.
- **No multi-consequence VCF handling** — CSQ has one entry per variant in these files; promoting to a separate table would add joins with no benefit here.
- **No ingestion_log table** — errors are logged to stderr; in production I would write to a table and emit CloudWatch metrics.

---

## At 1000× Scale

- Swap the DuckDB file for a cloud warehouse. The SQL queries are portable.
- Replace the manifest CSV cleaner with a Lambda triggered by S3 `PutObject` events on the landing bucket.
- Replace the VCF parser loop with a Spark job (EMR/Dataproc) that reads from S3 and writes Parquet to the curated bucket.
- Use Flyway or Alembic for versioned schema migrations instead of the current DDL-on-connect approach.
- Partition the variants table by `batch_id` and `chrom` for query performance.
- Route errors to a dead-letter queue and a monitoring dashboard (CloudWatch, Grafana).

---

## Production S3 Architecture (brief)

```
[Sequencer pipeline]
      ↓ s3:PutObject
[s3://genomics-landing/batch_YYYY_MM/]     ← pipeline service account only
      ↓ S3 event → Lambda → Step Function
  Step 1: validate files (quarantine malformed to /rejected/)
  Step 2: clean + transform (Fargate task)
  Step 3: load into curated store
  Step 4: data quality checks
      ↓
[s3://genomics-curated/]                   ← researcher role: s3:GetObject only
```

Each bucket has a separate IAM resource policy. Researchers never see the landing bucket.

---

## Governance

**De-identification:** `patient_id` is pseudonymous (P-0001), not de-identified. Re-identification is possible given diagnosis + tissue + collection_date. In production: date-shift collection_date by a per-patient offset; scrub the `notes` column (free text, highest PHI risk); document the pseudonymization approach in a data use agreement.

**Access control:** Landing zone is service-account only. Curated zone is researcher read-only (`SELECT`, no `INSERT/UPDATE/DELETE`). If the data use agreement requires it, create a researcher-facing VIEW that omits `patient_id`.

**AI assistant documentation:** A `data_dictionary.md` (see `docs/`) explains each column's meaning, units, valid values, and what NULL means in nullable columns. This is the contract that makes an NLP layer trustworthy — without it, an assistant can't know that `ccf IS NULL` means "batch_2026_01 didn't report this field" rather than "this sample has no cancer cells."

---

## AI Tool Usage

This project was built with Claude Code. Claude was used for:
- Design brainstorming and schema decisions
- Drafting all module code
- Writing tests (all verified by running them against real data)

I verified and understood every line before committing. In the live session I can walk through and extend any part of the code.
