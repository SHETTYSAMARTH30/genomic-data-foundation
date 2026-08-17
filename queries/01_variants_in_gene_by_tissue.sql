-- Which samples carry a KRAS variant and what tissue are they from?
-- Change 'KRAS' to any gene of interest.
SELECT
    s.sample_id,
    s.tissue,
    s.diagnosis,
    ROUND(s.tumor_purity, 2)  AS tumor_purity,
    v.consequence,
    v.hgvsp,
    ROUND(v.vaf, 4)           AS vaf
FROM variants v
JOIN samples s USING (sample_id)
WHERE v.gene = 'KRAS'
  AND v.filter = 'PASS'
ORDER BY v.vaf DESC;
