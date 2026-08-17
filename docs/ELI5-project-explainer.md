# What Are We Building? (Plain Language Explainer)

> This document explains the project in everyday language — no coding background needed.

---

## The Big Picture

Imagine a hospital lab that sequences patients' tumor DNA and produces a file for each patient after the analysis. These files contain a list of the genetic "typos" found in that patient's tumor (called **variants** or **mutations**). The lab also keeps a spreadsheet (the **manifest**) tracking who each patient is, what cancer they have, and details about their sample.

The problem: **these files pile up in a folder** and nobody has an easy way to search across all of them at once. A researcher can't easily ask "which of our lung cancer patients has a mutation in the KRAS gene?" — they'd have to open each file individually.

**Our job is to build the infrastructure that makes these questions easy to answer.**

---

## What Does the Data Look Like?

**The manifest (a spreadsheet):** One row per patient sample. Contains things like:
- Patient ID (anonymized, e.g. "P-0001")
- What kind of cancer they have ("Lung adenocarcinoma")
- What tissue the sample came from ("Lung")
- How pure the tumor sample was (tumor purity — a higher number means fewer healthy cells mixed in)
- Whether the sample passed quality checks

**The VCF files (one per sample):** "VCF" stands for Variant Call Format — it's a standard text file format used in genetics. Think of it as a list of every genetic "typo" found in that patient's tumor, along with details like:
- Where on the genome the typo is (chromosome and position)
- What the normal DNA letter is, and what the mutant letter is
- Which gene the typo falls in (e.g. KRAS, TP53, BRCA1)
- How confident the lab is that this is a real mutation
- What fraction of the tumor cells carry this mutation (variant allele fraction / VAF)

We have **two batches** of data:
- **Batch 1 (January 2026):** 14 samples
- **Batch 2 (February 2026):** 5 more samples (with some extra information not in Batch 1)

---

## What Problems Do We Need to Solve?

### Problem 1: The Data Is Messy

The manifest spreadsheet was filled in by humans (or inconsistent software), so:
- Dates are written in 6 different formats ("Jan 5 2026", "2026-01-07", "17/01/2026"…)
- "Lung" is sometimes written as "LUNG", "lung", " Lung" (with a space)
- Tumor purity is sometimes "0.32" (a fraction) and sometimes "63%" (a percentage)
- One patient ("S-0011") appears twice with slightly different information
- One sample ("S-0008") is listed in the spreadsheet but its data file is missing

We need to **clean all of this up automatically** — not by hand — so the database always has consistent, trustworthy data.

### Problem 2: The Two Batches Aren't Identical

Batch 2 has some new information that Batch 1 doesn't:
- A new column in the manifest: `library_prep` (how the DNA was prepared for sequencing)
- New fields in the VCF files: cancer cell fraction (CCF) and an additional population frequency database
- The variant fraction (VAF) numbers in Batch 2 are written as percentages (e.g. "56.14") instead of decimals ("0.56") — the same number, just expressed differently

We need to **handle these differences automatically** so that adding Batch 2 doesn't break any of the queries that were already working for Batch 1.

### Problem 3: No Single Place to Query

Right now the data is split across 19 separate files. We need a **single database** where a researcher can run one query and get answers that span the whole cohort.

---

## What Are We Building, Exactly?

### 1. A Pipeline (the "factory")

A Python script that reads the raw files and loads them into a clean database. Think of it as a factory line:

```
Raw messy files → [Clean & normalize] → [Validate] → [Load] → Clean database
```

Key properties:
- **Safe to run twice:** Running the pipeline again won't create duplicate rows or corrupt the data
- **Treats raw files as sacred:** The original files are never touched or modified
- **Handles errors gracefully:** If one file is broken, the rest still load; the broken one is logged

### 2. A Database (the "library")

We're using **DuckDB** — a lightweight but powerful analytical database that lives in a single file on your laptop. No server to set up. Think of it like a very smart spreadsheet that can answer complex questions instantly.

The database has two main tables:

**Samples table** — one row per patient sample:
- Who the patient is, what cancer they have, the tissue type, quality metrics, and batch info
- All messy values from the manifest are cleaned and standardized here

**Variants table** — one row per mutation per sample:
- The genomic location, the gene affected, the type of mutation
- How confident we are it's real, how common it is in the general population
- Whether it's a known cancer "hotspot" (a mutation that appears in many cancer patients)
- All the numeric values normalized to consistent units

### 3. Example Research Queries

Once the database is loaded, researchers can ask questions like:

- **"Which samples have a mutation in the KRAS gene, and what tissue are they from?"**
  *(Useful for finding patients who might respond to KRAS-targeting drugs)*

- **"What are the most commonly mutated genes across our entire cohort?"**
  *(Gives a bird's-eye view of the mutation landscape)*

- **"Show me all high-severity mutations in patients whose samples passed quality control"**
  *(Filters to the most clinically relevant findings)*

---

## Who Benefits and How?

| Who | Before | After |
|---|---|---|
| **Researcher** | Opens 19 files manually, copies data into Excel | Runs one SQL query, gets results in seconds |
| **Clinician** | Can't easily see which patients share a mutation | Queries the cohort in real time |
| **AI assistant** | No structured data to read from | Reads a clean, well-documented database |
| **Data engineer** | New batches break old queries | Batches load additively; old queries keep working |

---

## The Privacy Angle

The patient IDs in these files are already pseudonymized (e.g. "P-0001" instead of a real name). But that's not the same as fully de-identified — if you know someone's cancer type, tissue, and collection date, you might be able to figure out who they are.

In a real-world deployment, we'd also:
- Slightly shift the collection dates (so you can't match them to hospital records)
- Carefully review the "notes" column (free text that might accidentally contain a name or procedure detail)
- Give researchers read-only access to the clean database — they can query it but can't delete or modify anything

---

## The Future: A Natural-Language Assistant

Eventually, a researcher should be able to ask a question in plain English — like "how many lung cancer patients in our cohort carry a TP53 mutation?" — and get an answer without writing SQL.

For that to work reliably, the underlying database needs to be clean, well-documented, and consistent. Every column needs a clear definition. Every NULL value needs an explanation. Every unit needs to be standardized.

That's why this project matters: **we're not just loading data, we're building a foundation that a future AI can trust**.

---

## What We're Not Doing (and Why)

| Not doing | Why |
|---|---|
| Modifying the raw data files | They're the source of truth; we always want to be able to re-run from scratch |
| Setting up a cloud server | The exercise is designed to run locally; no cost, no account needed |
| Handling every possible VCF format | These files are standard VCFv4.2; a custom reader is cleaner than pulling in a complex library |
| Building a user interface | Out of scope; researchers interact via SQL queries today |

---

## How to Run It (for the technically curious)

```bash
# Install dependencies
pip install -r requirements.txt

# Load Batch 1
python src/ingest.py --batch data/batch_2026_01

# Load Batch 2 (safe to run after Batch 1; nothing breaks)
python src/ingest.py --batch data/batch_2026_02

# Open the database and run a query
python -c "import duckdb; print(duckdb.connect('curated/genomics.duckdb').sql('SELECT gene, COUNT(*) FROM variants GROUP BY gene ORDER BY 2 DESC LIMIT 10').df())"
```

That's it. One file, one command per batch, results you can trust.
