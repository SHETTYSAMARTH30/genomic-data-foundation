# Take-Home Exercise: Genomic Data Foundation (first slice)

**Time.** You have 3 days to turn this around, but please spend only about **two to three hours** of
focused work on it — we've deliberately kept the build small. We care about judgment, clarity, and
robustness far more than completeness; better to do less, well-reasoned and documented, than to
gold-plate. If you run short, describe the rest in your README rather than grinding; we'll go deeper
together in a live session.

**No paid infrastructure required.** Use a free, local engine — DuckDB, SQLite, or Parquet-on-disk
are all perfect. Please don't stand up any paid cloud services; we don't expect you to spend money
or set up a cloud account for this.

**AI tools: encouraged.** We build with AI-assisted tooling every day, so use whatever you normally
would (Claude Code, Copilot, Cursor). In your README, just note briefly how you used it and where you
didn't trust it. In a later live session we'll walk through your code and extend it, so build things
you understand.

## Context

Our sequencing pipeline produces per-sample result files plus a manifest and drops them into a
landing area (think an S3 prefix; here it's just a local folder). Researchers need to query this data
reliably, and eventually a natural-language assistant will answer questions on top of it — and that
assistant is only as reliable as the foundation underneath it. This exercise is the first slice of
that foundation.

You've been given two batches under `data/`:

- `batch_2026_01/` — one somatic VCF per sample (`S-0001.somatic.vcf`, …) plus a somewhat messy
  `sample_manifest.csv`.
- `batch_2026_02/` — a later drop that is **not identical in shape** to the first.

The VCFs are uncompressed VCFv4.2 (in production they arrive bgzipped and indexed). Any parsing
approach is fine — a library, or your own reader; there's nothing exotic in these files.

The data is small and synthetic; treat it as representative of something that will grow by orders of
magnitude.

## What we'd like you to build

1. **Ingest into a governed layout.** Treat the provided `data/` folder as an immutable landing zone
   (raw, exactly as the pipeline dropped it — don't modify it in place). Load, validate, and clean
   the batch-1 files, and land the results in a separate curated zone: a queryable store of your
   choice. Make the process re-runnable; running it twice should not corrupt or duplicate data.

2. **Model.** Design the schema researchers will query. In your README, state the **grain** of each
   table and why. Keep in mind the schema has to serve a SQL user today and a natural-language
   assistant later.

3. **Query.** Write 2–3 example queries that answer research-style questions. Include at least one
   that joins variants to sample metadata, and one that aggregates across the cohort (for example,
   "which samples carry a variant in gene X and what tissue are they from," or "count variants per
   gene"). The specific questions are yours to choose.

4. **Handle a changed batch.** Now ingest `batch_2026_02/`, which adds new INFO fields and a new
   manifest column. Handle them **additively** — without breaking your batch-1 queries and without
   hand-editing the data. The batch also includes one malformed file; you don't need to build
   handling for it, just describe in your README how you'd detect and deal with it in a real
   pipeline.

**Governance note (brief).** In a few sentences: how would you handle de-identification and access
control before researchers touch this, and what would you document so an AI assistant could answer
questions on this data reliably?

**README.** Your key decisions and tradeoffs, what you deliberately did *not* do, what you'd change
at 1000x scale, and a short note on how you'd run the landing-to-curated step as an S3-to-S3 workflow
in production, with each zone's access fenced appropriately.

## What we're looking for

Clear schema and grain decisions; a pipeline that's reproducible and doesn't fall over on messy or
changing input; documentation good enough that someone (or an agent) could trust the data without
asking you; and sound governance instincts. Clean, readable, tested code matters more than volume of
features. We're at least as interested in *how you reason*, in the README, as in what you finish.

## How to submit, and how we'll run it

Send a link to a Git repo (a private repo shared with us is fine) or a zip. Include a README that
lets us run the whole thing end to end with a single documented command from a clean checkout. Keep
dependencies minimal (a `requirements.txt` is fine; a Dockerfile is a nice-to-have, not required).
Because you're using a local engine, there should be nothing for us to provision — we clone, follow
your README, and run.

If anything is ambiguous, make a reasonable assumption, note it, and move on.
