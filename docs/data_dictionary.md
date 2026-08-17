# Data Dictionary — Genomic Data Foundation

This document is the contract between the curated database and any downstream
consumer — whether a SQL user, a BI tool, or an NLP assistant. Every column is
defined here: its meaning, its units, its valid values, and exactly what NULL
means. An AI assistant must not reason over this database without reading this
file first.

---

## Table: `samples`

**Grain:** one row per `(sample_id, batch_id)`.
A sample is one sequenced tumour specimen. The same patient (`patient_id`) could
theoretically appear across batches, so `batch_id` is part of the primary key.

| Column | Type | Nullable | Description |
|---|---|---|---|
| `sample_id` | VARCHAR | No | Lab-assigned specimen identifier. Format: `S-NNNN`. Unique within a batch; may repeat across batches for re-sequenced specimens. Part of the primary key. |
| `batch_id` | VARCHAR | No | Ingestion batch identifier. Format: `batch_YYYY_MM`. Reflects when the sequencing run was processed, not when the sample was collected. Part of the primary key. |
| `patient_id` | VARCHAR | Yes | Pseudonymous patient identifier. Format: `P-NNNN`. **Not de-identified** — see Governance below. NULL means the manifest row had no patient ID. |
| `collection_date` | DATE | Yes | Date the tumour specimen was collected. Normalised from six different input formats (ISO, slash-separated, natural language). NULL means the date was blank or unparseable in the manifest. |
| `tissue` | VARCHAR | Yes | Tumour tissue of origin, lowercased and stripped. Examples: `lung`, `colon`, `breast`, `pancreas`, `stomach`. NULL means the manifest field was blank. |
| `diagnosis` | VARCHAR | Yes | Free-text clinical diagnosis as recorded in the manifest. Examples: `Lung adenocarcinoma`, `Invasive ductal carcinoma`. Not normalised — do not group by this column without accounting for spelling variation. |
| `disease_group` | VARCHAR | Yes | Coarser disease category. Examples: `Thoracic`, `Gastrointestinal`, `Breast`. Useful for cohort-level grouping without relying on free-text `diagnosis`. NULL means the manifest field was blank. |
| `tumor_purity` | FLOAT | Yes | Estimated fraction of tumour cells in the specimen. **Always stored as a fraction in [0.0, 1.0]**, regardless of how the source manifest expressed it (`63%` → `0.63`). NULL means the manifest field was blank or unparseable — it does NOT mean zero purity. |
| `sex_reported` | VARCHAR | Yes | Biological sex as reported in the manifest. Normalised: `F` (female) or `M` (male). NULL means the field was blank. |
| `sex_inferred` | VARCHAR | Yes | Biological sex inferred from sequencing data. Same normalisation as `sex_reported`. NULL means not computed or not reported. |
| `sequencing_platform` | VARCHAR | Yes | Sequencing instrument used. Normalised to canonical names: `NovaSeq 6000`, `NovaSeq X`. NULL means the manifest field was blank. |
| `qc_status` | VARCHAR | Yes | Quality-control outcome. Normalised to uppercase: `PASS` or `FAIL`. NULL means the field was blank — treat as unknown, not as PASS. |
| `notes` | VARCHAR | Yes | Free-text notes from the manifest. **Highest PHI risk** — may contain names, procedure details, or other identifying information. Do not expose in researcher-facing views without scrubbing. NULL means the field was blank. |
| `library_prep` | VARCHAR | Yes | DNA library preparation protocol. **NULL means the sample is from `batch_2026_01`**, which did not report this field. A non-NULL value indicates `batch_2026_02` or later. Do not interpret NULL as "unknown protocol". |
| `ingested_at` | TIMESTAMP | No | UTC timestamp when this row was written to the database by the pipeline. Not a clinical date — for audit and lineage purposes only. |

---

## Table: `variants`

**Grain:** one row per `(sample_id, batch_id, chrom, pos, ref, alt)`.
A variant is one somatic mutation called in one sample. The locus
`(chrom, pos, ref, alt)` uniquely identifies the mutation within a sample.

### Identity and location

| Column | Type | Nullable | Description |
|---|---|---|---|
| `sample_id` | VARCHAR | No | Links to `samples.sample_id`. |
| `batch_id` | VARCHAR | No | Denormalised from the ingestion batch. Included in the primary key so that a re-sequenced sample in a later batch does not overwrite earlier results. Also useful for filtering without joining to `samples`. |
| `chrom` | VARCHAR | No | Chromosome. Format: `chr1`, `chr2`, …, `chrX`, `chrY`. Taken verbatim from the VCF CHROM column. |
| `pos` | BIGINT | No | 1-based genomic position of the variant on the chromosome. Taken verbatim from the VCF POS column. |
| `vcf_id` | VARCHAR | Yes | dbSNP or other public variant identifier (rsID). NULL means the VCF had `.` in the ID column, indicating no public ID was assigned — it does NOT mean the variant is novel. |
| `ref` | VARCHAR | No | Reference allele at this position (the "normal" DNA sequence). |
| `alt` | VARCHAR | No | Alternate (mutant) allele observed in the tumour. |

### Call quality

| Column | Type | Nullable | Description |
|---|---|---|---|
| `qual` | FLOAT | Yes | Phred-scaled quality score for the variant call. Higher is more confident. NULL means the VCF had `.` in the QUAL column — the caller did not assign a score. |
| `filter` | VARCHAR | Yes | Filter outcome from the variant caller. Common values: `PASS` (high-confidence call), `weak_evidence`, `low_depth`, `germline_risk`. Only `PASS` variants should be used in clinical or research analyses. NULL means no filter was applied. |
| `caller` | VARCHAR | Yes | Variant-calling software that produced this call. Values observed: `mutect2`, `strelka2`, `mutect2_strelka2` (called by both). NULL means the caller was not recorded. |

### Population and clinical annotations

| Column | Type | Nullable | Description |
|---|---|---|---|
| `max_pop_af` | FLOAT | Yes | Maximum allele frequency of this variant across population databases. Range: [0.0, 1.0]. High values (e.g. > 0.01) suggest the variant may be a common germline polymorphism rather than a somatic mutation. |
| `hotspot` | BOOLEAN | Yes | Whether this position is a known cancer mutation hotspot. `TRUE` means the position appears recurrently across cancer patients in external databases. `FALSE` means it was not flagged. Derived from the `HOTSPOT` FLAG field in the VCF INFO column — a FLAG with no value means True; absence means False. |
| `ccf` | FLOAT | Yes | Cancer Cell Fraction — the estimated fraction of cancer cells in the tumour that carry this specific mutation. Range: [0.0, 1.0]. **NULL means the sample is from `batch_2026_01`**, which did not report CCF. A non-NULL value indicates `batch_2026_02` or later. Do not treat NULL as zero CCF. |
| `gnomad_af_popmax` | FLOAT | Yes | gnomAD population-maximum allele frequency — the highest allele frequency of this variant across any gnomAD population group. Range: [0.0, 1.0]. **NULL means the sample is from `batch_2026_01`**, which did not report this field. |

### Consequence annotation (from VCF CSQ field)

These columns are parsed from the `CSQ` INFO field, which encodes the predicted
functional effect of the variant on the gene's transcript.

| Column | Type | Nullable | Description |
|---|---|---|---|
| `gene` | VARCHAR | Yes | HGNC gene symbol. Examples: `KRAS`, `TP53`, `BRCA1`. NULL means no gene annotation was present in the CSQ field. |
| `consequence` | VARCHAR | Yes | Sequence Ontology consequence term. Examples: `missense_variant`, `frameshift_variant`, `synonymous_variant`, `stop_gained`, `splice_donor_variant`, `intron_variant`. Describes what the mutation does to the transcript. |
| `impact` | VARCHAR | Yes | High-level impact tier assigned by the annotation tool. Values: `HIGH` (likely loss of function), `MODERATE` (potentially alters protein), `LOW` (likely benign), `MODIFIER` (non-coding or unclear). |
| `hgvsc` | VARCHAR | Yes | HGVS coding-sequence notation. Format: `c.<position><change>`. Example: `c.35G>T`. Describes the change at the DNA level within the coding sequence. NULL means not applicable or not annotated. |
| `hgvsp` | VARCHAR | Yes | HGVS protein notation. Format: `p.<Ref><position><Alt>`. Example: `p.Gly12Val`. Describes the resulting amino acid change. NULL means the variant does not affect the protein (e.g. synonymous, intronic) or was not annotated. |
| `exon` | VARCHAR | Yes | Exon number and total exon count. Format: `<exon>/<total>`. Example: `2/3`. NULL means the variant is not in an exon (intronic, UTR, etc.). |
| `mane_select` | VARCHAR | Yes | MANE Select transcript identifier — the single representative transcript for this gene agreed upon by NCBI and Ensembl. Example: `ENST00000256078.10`. **NULL means the sample is from `batch_2026_01`**, which did not report this field. |

### Genotype and read support

| Column | Type | Nullable | Description |
|---|---|---|---|
| `gt` | VARCHAR | Yes | Genotype call. `0/1` = heterozygous (one copy mutant), `1/1` = homozygous (both copies mutant). Taken verbatim from the GT field in the VCF sample column. |
| `ad_ref` | INTEGER | Yes | Number of sequencing reads supporting the reference (normal) allele at this position. NULL means the AD field was absent or unparseable. |
| `ad_alt` | INTEGER | Yes | Number of sequencing reads supporting the alternate (mutant) allele at this position. NULL means the AD field was absent or unparseable. |
| `dp` | INTEGER | Yes | Total sequencing read depth at this position (`ad_ref + ad_alt` approximately). Higher depth means a more reliable call. NULL means the DP field was absent or had the value `.`. |
| `vaf` | FLOAT | Yes | Variant Allele Fraction — the proportion of reads at this position carrying the mutant allele. **Always stored as a fraction in [0.0, 1.0]**, regardless of how the source VCF expressed it. `batch_2026_02` reported VAF as a percentage (0–100); the pipeline divides by 100 before storing. NULL means VAF was not computable. |

### Audit

| Column | Type | Nullable | Description |
|---|---|---|---|
| `ingested_at` | TIMESTAMP | No | UTC timestamp when this row was written to the database. For audit and lineage purposes only. |

---

## Known Data Quality Issues

| Issue | Affected rows | Detail |
|---|---|---|
| `tumor_purity` NULL | S-0003 | Purity was blank in the batch_2026_01 manifest. Treat as unknown, not zero. |
| S-0008 has no variants | S-0008 | The sample appears in the batch_2026_01 manifest but no VCF file was delivered. The sample row exists; no variant rows exist. This is a known data delivery gap, not a pipeline error. |
| S-0020 line 101 skipped | S-0020 | One row in S-0020.somatic.vcf was truncated (8 columns instead of 10). That row was skipped with a WARNING; the remaining 63 variants from the file were loaded normally. |
| `ccf`, `gnomad_af_popmax`, `mane_select`, `library_prep` NULL for batch_2026_01 | All batch_2026_01 rows | These fields were not reported by the batch_2026_01 pipeline. NULL here means "not available for this batch", not "missing data". |

---

## Governance Notes

**Pseudonymisation:** `patient_id` (format: `P-NNNN`) is a pseudonym, not a
de-identified identifier. A re-identification attack is possible if an adversary
knows a patient's cancer type, tissue, and approximate collection date. Do not
share the raw database outside the research team.

**PHI risk:** The `notes` column is free text from the manifest and may
accidentally contain names, dates of birth, procedure codes, or other
identifying details. Any researcher-facing view should omit or scrub this column.

**Date shifting:** `collection_date` has not been shifted. In a production
deployment, dates should be shifted by a per-patient random offset before
researcher access to prevent record linkage against hospital admission databases.

**Access model:** The pipeline service account writes to the database. Researchers
should have read-only access (`SELECT` only — no `INSERT`, `UPDATE`, or `DELETE`).
If the data use agreement requires it, create a view that omits `patient_id` and
`notes`.

**For NLP assistants:** Before generating SQL from a natural-language question,
verify which columns are relevant, their units, and the NULL semantics above.
In particular: `NULL` in `ccf`, `gnomad_af_popmax`, `mane_select`, and
`library_prep` always means "batch_2026_01 did not report this field" — never
filter these columns with `IS NOT NULL` unless you intend to exclude all
batch_2026_01 data.
