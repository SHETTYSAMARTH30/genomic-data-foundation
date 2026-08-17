from __future__ import annotations

import pytest
from vcf_parser import parse_info, parse_csq, normalize_vaf, parse_vcf


# --- parse_info ---

def test_parse_info_basic_fields():
    result = parse_info("CALLER=mutect2;MAX_POP_AF=0.001")
    assert result["CALLER"] == "mutect2"
    assert result["MAX_POP_AF"] == "0.001"

def test_parse_info_flag_field_present():
    result = parse_info("CALLER=mutect2;MAX_POP_AF=0.0;HOTSPOT;CSQ=GENE|mis|MOD|c.1A>G|p.X1Y|1/5")
    assert result["HOTSPOT"] is True

def test_parse_info_flag_field_absent():
    result = parse_info("CALLER=strelka2;MAX_POP_AF=0.001")
    assert "HOTSPOT" not in result

def test_parse_info_csq_preserved():
    csq_val = "GENE1|missense_variant|MODERATE|c.100A>G|p.Arg34Gly|1/5"
    result = parse_info(f"CALLER=mutect2;CSQ={csq_val}")
    assert result["CSQ"] == csq_val


# --- parse_csq ---

def test_parse_csq_six_fields():
    csq = "ARID1A|missense_variant|MODERATE|c.234A>G|p.Gln79Met|6/20"
    r = parse_csq(csq)
    assert r["gene"] == "ARID1A"
    assert r["consequence"] == "missense_variant"
    assert r["impact"] == "MODERATE"
    assert r["hgvsc"] == "c.234A>G"
    assert r["hgvsp"] == "p.Gln79Met"
    assert r["exon"] == "6/20"
    assert r["mane_select"] is None

def test_parse_csq_seven_fields():
    csq = "ARID1A|inframe_insertion|MODERATE|c.2865_2866insGGT|p.Asp956Gly|4/20|ENST00081104912.7"
    r = parse_csq(csq)
    assert r["mane_select"] == "ENST00081104912.7"

def test_parse_csq_empty_hgvsp_and_exon():
    # UTR variant: HGVSp and exon fields are empty in batch 1
    csq = "ARID1A|5_prime_UTR_variant|MODIFIER|c.-199G>T||"
    r = parse_csq(csq)
    assert r["gene"] == "ARID1A"
    assert r["hgvsp"] is None
    assert r["exon"] is None


# --- normalize_vaf ---

def test_normalize_vaf_fraction_unchanged():
    assert normalize_vaf(0.3390) == pytest.approx(0.3390)

def test_normalize_vaf_percentage_divided():
    assert normalize_vaf(56.14) == pytest.approx(0.5614)

def test_normalize_vaf_exactly_one():
    # 1.0 is not > 1.0, so it stays (100% VAF, rare but valid)
    assert normalize_vaf(1.0) == pytest.approx(1.0)

def test_normalize_vaf_small_percentage():
    assert normalize_vaf(9.66) == pytest.approx(0.0966)


# --- parse_vcf (integration against real files) ---

def test_parse_vcf_batch1_returns_variants():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    assert len(variants) > 0

def test_parse_vcf_batch1_sample_id_set():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    assert all(v["sample_id"] == "S-0001" for v in variants)

def test_parse_vcf_batch1_vaf_in_range():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    for v in variants:
        assert v["vaf"] is None or 0.0 <= v["vaf"] <= 1.0, f"VAF out of range: {v['vaf']}"

def test_parse_vcf_hotspot_is_boolean():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    assert all(isinstance(v["hotspot"], bool) for v in variants)

def test_parse_vcf_batch1_ccf_is_none():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    assert all(v["ccf"] is None for v in variants)

def test_parse_vcf_batch2_vaf_normalized():
    # Batch 2 reports VAF as percentage — must be normalized to 0–1
    variants = parse_vcf("data/batch_2026_02/S-0016.somatic.vcf", "S-0016", "batch_2026_02")
    for v in variants:
        assert v["vaf"] is None or 0.0 <= v["vaf"] <= 1.0, f"VAF out of range: {v['vaf']}"

def test_parse_vcf_batch2_has_ccf():
    variants = parse_vcf("data/batch_2026_02/S-0016.somatic.vcf", "S-0016", "batch_2026_02")
    assert any(v["ccf"] is not None for v in variants)

def test_parse_vcf_skips_malformed_row():
    # S-0020.somatic.vcf line 101 has only 8 columns — should be skipped
    variants = parse_vcf("data/batch_2026_02/S-0020.somatic.vcf", "S-0020", "batch_2026_02")
    assert len(variants) > 0
    # The truncated row is chr17:43122442 — must not appear
    bad = [v for v in variants if v["chrom"] == "chr17" and v["pos"] == 43122442]
    assert bad == []

def test_parse_vcf_dot_id_becomes_none():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    # Most variants have "." as ID — should be None
    dot_ids = [v for v in variants if v["vcf_id"] == "."]
    assert dot_ids == []

def test_parse_vcf_dot_qual_becomes_none():
    variants = parse_vcf("data/batch_2026_01/S-0001.somatic.vcf", "S-0001", "batch_2026_01")
    dot_quals = [v for v in variants if v["qual"] == "."]
    assert dot_quals == []
