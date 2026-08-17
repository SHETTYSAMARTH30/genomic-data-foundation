-- All HIGH-impact variants in QC-passing samples, sorted by VAF.
-- HIGH-impact includes stop_gained, frameshift_variant, splice_acceptor_variant.
SELECT
    s.sample_id,
    s.tissue,
    s.disease_group,
    v.gene,
    v.consequence,
    v.hgvsp,
    ROUND(v.vaf, 4) AS vaf,
    v.hotspot,
    v.caller
FROM variants v
JOIN samples s USING (sample_id)
WHERE v.impact = 'HIGH'
  AND v.filter = 'PASS'
  AND s.qc_status = 'PASS'
ORDER BY v.vaf DESC;
