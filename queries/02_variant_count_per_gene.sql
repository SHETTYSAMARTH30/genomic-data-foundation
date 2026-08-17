-- Count PASS variants per gene across the entire cohort.
-- Shows which genes are most frequently mutated and in how many samples.
SELECT
    gene,
    COUNT(DISTINCT sample_id) AS n_samples,
    COUNT(*)                  AS n_variants,
    ROUND(AVG(vaf), 4)        AS mean_vaf
FROM variants
WHERE filter = 'PASS'
GROUP BY gene
ORDER BY n_samples DESC, n_variants DESC;
