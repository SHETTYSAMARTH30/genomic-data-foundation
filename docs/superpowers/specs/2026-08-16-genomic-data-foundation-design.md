# Design Spec: Genomic Data Foundation — First Slice

**Date:** 2026-08-16
**Status:** Approved
**Scope:** Batch ingestion pipeline, curated schema, example queries, schema evolution, governance

---

## 1. Context & Goals

The sequencing pipeline drops per-sample VCF files and a manifest CSV into a local landing area (`data/`). Researchers need a clean, queryable store. A natural-language assistant will eventually sit on top of it — so the foundation must be trustworthy enough for a model to reason over without asking us to interpret ambiguities.

This slice covers:
- Ingesting two batches (batch_2026_01, batch_2026_02) with messy manifests and evolving VCF schemas
- A curated DuckDB store that is idempotent and additive across batches
- 2–3 research-style queries
- Brief governance notes

---

## 2. Technology Choices

| Concern | Choice | Reason |
|---|---|---|
| Language | Python | Standard for genomics; fits `requirements.txt` |
| Curated store | DuckDB | Columnar, OLAP-optimized, single-file, native Python API |
| VCF parser | Custom line-by-line reader | Zero compiled dependencies; files have no exotic VCF features |
| Schema style | Moderately normalized (2 tables) | Clean joins for SQL users; simple enough for NLP to navigate |

---

## 3. Directory Layout

```
candidate_bundle/
├── data/                          ← IMMUTABLE landing zone (never modified)
│   ├── batch_2026_01/
│   └── batch_2026_02/
├── curated/                       ← Pipeline output (gitignored)
│   └── genomics.duckdb
├── src/
│   ├── ingest.py                  ← CLI entrypoint
│   ├── vcf_parser.py              ← Pure VCF parser (no I/O side effects)
│   ├── manifest_cleaner.py        ← CSV normalization
│   └── schema.py                  ← DDL + upsert helpers + migration
├── queries/
│   ├── 01_variants_in_gene_by_tissue.sql
│   ├── 02_variant_count_per_gene.sql
│   └── 03_high_impact_with_sample_metadata.sql
├── tests/
│   ├── test_vcf_parser.py
│   └── test_manifest_cleaner.py
├── requirements.txt
└── README.md
```

**Run command:**
```bash
pip install -r requirements.txt
python src/ingest.py --batch data/batch_2026_01
python src/ingest.py --batch data/batch_2026_02   # additive; batch-1 queries unaffected
```

---

## 4. Schema Design

### 4.1 `samples` table

**Grain:** One row per `(sample_id, batch_id)`.

Why this grain: a sample is uniquely identified by its ID within a batch. The same patient could theoretically appear across batches (though not in this dataset), so `patient_id` is not the PK.

| Column | Type | Notes |
|---|---|---|
| `sample_id` | VARCHAR | e.g. S-0001 |
| `batch_id` | VARCHAR | e.g. batch_2026_01 |
| `patient_id` | VARCHAR | Pseudonymous identifier |
| `collection_date` | DATE | Normalized from 6 messy formats |
| `tissue` | VARCHAR | Lowercased + stripped (lung, colon, breast…) |
| `diagnosis` | VARCHAR | Free text |
| `disease_group` | VARCHAR | e.g. Thoracic, Gastrointestinal |
| `tumor_purity` | FLOAT | Normalized to 0.0–1.0 (63% → 0.63) |
| `sex_reported` | VARCHAR | Normalized to M/F |
| `sex_inferred` | VARCHAR | Normalized to M/F |
| `sequencing_platform` | VARCHAR | Normalized (novaseq 6000 → NovaSeq 6000) |
| `qc_status` | VARCHAR | Normalized to PASS/FAIL |
| `notes` | VARCHAR | Raw, may contain PHI — see governance |
| `library_prep` | VARCHAR | **Nullable** — batch_2026_02 only |
| `ingested_at` | TIMESTAMP | Pipeline run timestamp |

**Primary key:** `(sample_id, batch_id)`

**Duplicate handling:** S-0011 appears twice in batch_2026_01 with conflicting `tumor_purity`. Resolution: keep the row with a decimal-format purity (not percentage string); log a WARNING with both rows.

**Missing VCF:** S-0008 is in the batch_2026_01 manifest but has no corresponding VCF. The manifest row is loaded (sample metadata is valid). No variant rows are created for it. Logged as a WARNING.

---

### 4.2 `variants` table

**Grain:** One row per `(sample_id, chrom, pos, ref, alt)`.

Why this grain: a variant is uniquely identified by its genomic locus + alleles within a sample. Multiple consequences per variant are not present in these files (CSQ has one entry per variant row); if that changes, CSQ can be promoted to its own table.

| Column | Type | Notes |
|---|---|---|
| `sample_id` | VARCHAR | FK → samples |
| `batch_id` | VARCHAR | Denormalized for filtering without join |
| `chrom` | VARCHAR | e.g. chr1 |
| `pos` | BIGINT | |
| `vcf_id` | VARCHAR | rsID, nullable |
| `ref` | VARCHAR | |
| `alt` | VARCHAR | |
| `qual` | FLOAT | Nullable (`.` in VCF → NULL) |
| `filter` | VARCHAR | PASS, weak_evidence, low_depth, germline_risk |
| `caller` | VARCHAR | e.g. mutect2, strelka2, mutect2_strelka2 |
| `max_pop_af` | FLOAT | |
| `hotspot` | BOOLEAN | Flag field — present → TRUE, absent → FALSE |
| `ccf` | FLOAT | **Nullable** — Cancer cell fraction; batch_2026_02 only |
| `gnomad_af_popmax` | FLOAT | **Nullable** — batch_2026_02 only |
| `gene` | VARCHAR | From CSQ field 1 |
| `consequence` | VARCHAR | From CSQ field 2 (e.g. missense_variant) |
| `impact` | VARCHAR | HIGH / MODERATE / LOW / MODIFIER |
| `hgvsc` | VARCHAR | Coding sequence change notation |
| `hgvsp` | VARCHAR | Protein change notation |
| `exon` | VARCHAR | Exon number/total |
| `mane_select` | VARCHAR | **Nullable** — 7th CSQ field; batch_2026_02 only |
| `gt` | VARCHAR | Genotype (0/1, 1/1) |
| `ad_ref` | INTEGER | Ref allele read depth |
| `ad_alt` | INTEGER | Alt allele read depth |
| `dp` | INTEGER | Total read depth |
| `vaf` | FLOAT | **Normalized to 0.0–1.0** (batch_2026_02 reports 0–100; divide by 100) |
| `ingested_at` | TIMESTAMP | |

**Primary key:** `(sample_id, chrom, pos, ref, alt)`

---

## 5. Pipeline Logic

### 5.1 Manifest cleaning (`manifest_cleaner.py`)

1. **Dates** — cascade through 6 known formats with `dateutil.parser`; NULL + WARNING if no format matches
2. **Tissue** — `strip().lower()`
3. **Tumor purity** — strip `%`; if value > 1.0, divide by 100; blank → NULL
4. **Sex** — `{f, F, female, Female} → F`, `{m, M, male, Male} → M`
5. **Platform** — lowercase + deduplicate spelling variants into canonical names
6. **QC status** — `.upper()`; blank → NULL
7. **Deduplicate** — group by `(sample_id, batch_id)`; keep most complete row; log conflicts

### 5.2 VCF parsing (`vcf_parser.py`)

1. Parse `##INFO` and `##FORMAT` header lines → build field registry (makes the parser schema-aware per file)
2. For each data row: split by tab; parse INFO as `key=value` pairs; handle FLAG fields (HOTSPOT has no `=value`)
3. Parse CSQ: `|`-delimited; map by position; length-check to handle 6 vs 7 fields gracefully
4. Parse FORMAT/sample column
5. **VAF normalization:** if `vaf > 1.0`, divide by 100

### 5.3 Idempotency (`schema.py`)

- `INSERT OR IGNORE` on primary keys — re-running skips already-loaded rows without error or duplication
- **Schema migration for batch_2026_02:** run `ALTER TABLE variants ADD COLUMN IF NOT EXISTS ccf FLOAT` (and similarly for other new columns) before loading. DuckDB supports this natively. Existing batch_2026_01 rows keep NULL in new columns — batch_2026_01 queries are unaffected.

### 5.4 Error handling

Each file is wrapped in a try/except. Failures are logged to `stderr` with the filename and error. The pipeline continues with remaining files. In production this would write to an `ingestion_log` table and emit a metric.

---

## 6. Example Queries

### 6.1 Which samples carry a KRAS variant and what tissue are they from?
```sql
SELECT s.sample_id, s.tissue, s.diagnosis, s.tumor_purity,
       v.consequence, v.hgvsp, v.vaf
FROM variants v
JOIN samples s USING (sample_id)
WHERE v.gene = 'KRAS'
  AND v.filter = 'PASS'
ORDER BY v.vaf DESC;
```

### 6.2 Count PASS variants per gene across the cohort
```sql
SELECT gene,
       COUNT(DISTINCT sample_id) AS n_samples,
       COUNT(*) AS n_variants,
       ROUND(AVG(vaf), 4) AS mean_vaf
FROM variants
WHERE filter = 'PASS'
GROUP BY gene
ORDER BY n_samples DESC;
```

### 6.3 All HIGH-impact variants with sample metadata, sorted by VAF
```sql
SELECT s.sample_id, s.tissue, s.disease_group,
       v.gene, v.consequence, v.hgvsp, v.vaf,
       v.hotspot, v.caller
FROM variants v
JOIN samples s USING (sample_id)
WHERE v.impact = 'HIGH'
  AND v.filter = 'PASS'
  AND s.qc_status = 'PASS'
ORDER BY v.vaf DESC;
```

---

## 7. Handling Batch_2026_02 — Schema Evolution

Changes from batch_2026_01:
1. **New manifest column** `library_prep` — nullable in `samples`; batch_2026_01 rows have NULL
2. **New INFO fields** `CCF`, `GNOMAD_AF_POPMAX` — nullable columns in `variants`
3. **Extended CSQ** — 7th field `MANE_SELECT`; parser handles by position with length check
4. **VAF as percentage** (0–100 instead of 0–1) — normalized in parser
5. **One malformed file** — caught by try/except; logged; other files continue

Approach: additive schema migration via `ALTER TABLE … ADD COLUMN IF NOT EXISTS`. No existing queries touch the new columns, so they remain valid.

**The malformed file:** `S-0020.somatic.vcf` line 101 is truncated mid-row — the CSQ annotation and the entire FORMAT/sample genotype columns are missing, leaving only 8 of the expected 10 tab-delimited columns. The file header is valid; the issue is a truncated data row (likely a write failure or disk issue during sequencer output).

**Detecting and handling the malformed file in production:**
- Row-level validation during parsing: check each data row has exactly 10 tab-delimited columns; log a per-row WARNING and skip that row (don't abort the whole file)
- File-level validation before parsing: confirm the file is readable, has a `##fileformat=VCFv4` header, and a `#CHROM` column header
- On file-level failure: move to `rejected/` prefix; write to `ingestion_log`; fire an alert (SNS/PagerDuty)
- For partial files (valid header, bad rows): load good rows, skip bad rows, surface a count of skipped rows
- Do not block the rest of the batch

---

## 8. Governance

### De-identification
- `patient_id` (P-0001) is pseudonymous, not de-identified — a re-identification attack via diagnosis + tissue + date is possible
- `notes` field is highest risk (free text; could contain names, procedure details)
- `collection_date` should be date-shifted by a per-patient offset before researcher access
- Real-world: IRB-approved de-identification protocol; document what was done and when

### Access Control
- **Landing zone:** pipeline service account only — no researcher access
- **Curated zone:** researcher role gets `SELECT` only; no `INSERT/UPDATE/DELETE`
- **Column-level:** create a researcher-facing VIEW that omits `patient_id` if required by the data use agreement
- **S3 equivalent:** separate buckets with bucket policies; researchers assume a role with `s3:GetObject` on curated only

### AI Assistant Documentation
A `data_dictionary.md` (committed alongside the schema) that explains:
- Each column's meaning, units, and valid values
- What NULL means in each nullable column (e.g., `ccf IS NULL` means the sample is from batch_2026_01 which did not report CCF)
- Known data quality issues (some samples have NULL tumor_purity; batch_2026_01 VAF was in 0–1, already normalized)
- The pseudonymization approach and what researchers should not infer from `patient_id`

This is the contract between the data and the NLP layer.

---

## 9. What Was Deliberately Not Done

- No Parquet layer (overkill for synthetic data at this scale)
- No Airflow/workflow orchestration (single CLI script is sufficient and cleaner for review)
- No multi-consequence VCF handling (not present in these files)
- No full normalization of `caller` into a separate table (only 3 distinct values)
- No Dockerfile (noted as nice-to-have, not required)

---

## 10. Scale Considerations (1000x)

- Swap DuckDB file for DuckDB on S3 or migrate to BigQuery/Redshift/Snowflake
- VCF parsing becomes a Spark or dbt job; manifest cleaning becomes a Lambda triggered by S3 `PutObject`
- Landing zone → curated zone as separate S3 buckets with IAM boundary policies
- Schema migrations via a proper migration tool (Flyway, Alembic) with versioned scripts
- `ingestion_log` table becomes a proper observability layer (CloudWatch, DataDog)
- Partitioning variants table by `batch_id` and `chrom` for query performance

---

## 11. Production S3 Workflow (brief)

```
[Sequencer pipeline]
      ↓ s3:PutObject
[s3://raw-genomics/landing/batch_YYYY_MM/]   ← service account only
      ↓ S3 event → Lambda (validates + triggers Step Function)
[Step Function]
  Step 1: validate manifest + VCFs (quarantine malformed to /rejected/)
  Step 2: clean + transform (ECS Fargate task running this pipeline)
  Step 3: load into curated store
  Step 4: run data quality checks; fail the pipeline if critical checks fail
      ↓
[s3://curated-genomics/genomics.parquet/]    ← researcher role: GetObject only
      ↓
[Athena / researcher SQL client]
```

Each zone's bucket has a separate IAM resource policy. The pipeline role has write access to curated; researchers never see the landing bucket.
