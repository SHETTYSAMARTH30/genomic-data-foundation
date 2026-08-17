import pytest
from manifest_cleaner import (
    normalize_date, normalize_purity, normalize_sex,
    normalize_platform, normalize_qc_status, clean_manifest,
)


# --- normalize_date ---

def test_normalize_date_iso():
    assert normalize_date("2026-01-07") == "2026-01-07"

def test_normalize_date_slash_iso():
    assert normalize_date("2026/01/09") == "2026-01-09"

def test_normalize_date_day_first():
    assert normalize_date("17/01/2026") == "2026-01-17"

def test_normalize_date_day_first_zero_padded():
    assert normalize_date("04/01/2026") == "2026-01-04"

def test_normalize_date_month_name_padded():
    assert normalize_date("Jan 25 2026") == "2026-01-25"

def test_normalize_date_month_name_single_digit():
    assert normalize_date("Jan 5 2026") == "2026-01-05"

def test_normalize_date_feb():
    assert normalize_date("Feb 14 2026") == "2026-02-14"

def test_normalize_date_empty():
    assert normalize_date("") is None

def test_normalize_date_whitespace_only():
    assert normalize_date("   ") is None


# --- normalize_purity ---

def test_normalize_purity_decimal():
    assert normalize_purity("0.32") == pytest.approx(0.32)

def test_normalize_purity_decimal_high():
    assert normalize_purity("0.772") == pytest.approx(0.772)

def test_normalize_purity_percent():
    assert normalize_purity("63%") == pytest.approx(0.63)

def test_normalize_purity_percent_int():
    assert normalize_purity("56%") == pytest.approx(0.56)

def test_normalize_purity_empty():
    assert normalize_purity("") is None

def test_normalize_purity_all_in_range():
    for raw in ["0.32", "0.772", "63%", "57%", "0.70", "36%"]:
        result = normalize_purity(raw)
        assert result is not None and 0.0 <= result <= 1.0, f"Out of range for {raw!r}: {result}"


# --- normalize_sex ---

def test_normalize_sex_female_variants():
    for v in ["f", "F", "female", "Female"]:
        assert normalize_sex(v) == "F", f"Expected F for {v!r}"

def test_normalize_sex_male_variants():
    for v in ["m", "M", "male", "Male"]:
        assert normalize_sex(v) == "M", f"Expected M for {v!r}"

def test_normalize_sex_empty():
    assert normalize_sex("") is None

def test_normalize_sex_whitespace():
    assert normalize_sex("  ") is None


# --- normalize_platform ---

def test_normalize_platform_already_canonical_6000():
    assert normalize_platform("NovaSeq 6000") == "NovaSeq 6000"

def test_normalize_platform_lowercase_6000():
    assert normalize_platform("novaseq 6000") == "NovaSeq 6000"

def test_normalize_platform_no_space_6000():
    assert normalize_platform("NovaSeq6000") == "NovaSeq 6000"

def test_normalize_platform_novaseq_x():
    assert normalize_platform("NovaSeq X") == "NovaSeq X"

def test_normalize_platform_novaseq_x_lowercase():
    assert normalize_platform("novaseq x") == "NovaSeq X"


# --- normalize_qc_status ---

def test_normalize_qc_status_pass_lowercase():
    assert normalize_qc_status("pass") == "PASS"

def test_normalize_qc_status_pass_upper():
    assert normalize_qc_status("PASS") == "PASS"

def test_normalize_qc_status_empty():
    assert normalize_qc_status("") is None


# --- clean_manifest integration ---

def test_clean_manifest_deduplicates_s0011(tmp_path):
    csv_content = (
        "sample_id,patient_id,batch_id,collection_date,tissue,diagnosis,"
        "disease_group,tumor_purity,sex_reported,sex_inferred,sequencing_platform,qc_status,notes\n"
        "S-0011,P-0011,batch_2026_01,Jan 25 2026,lung,Lung adenocarcinoma,"
        "Thoracic,57%,f,F,NovaSeq 6000,PASS,\n"
        "S-0011,P-0011,batch_2026_01,Jan 25 2026,lung,Lung adenocarcinoma,"
        "Thoracic,0.64,f,F,NovaSeq 6000,PASS,\n"
    )
    p = tmp_path / "manifest.csv"
    p.write_text(csv_content)
    rows = clean_manifest(str(p))
    s0011 = [r for r in rows if r["sample_id"] == "S-0011"]
    assert len(s0011) == 1
    # Design spec: keep the decimal-format purity row (0.64), not the percentage row (57% -> 0.57)
    assert s0011[0]["tumor_purity"] == pytest.approx(0.64)


def test_clean_manifest_batch1_all_purities_in_range():
    rows = clean_manifest("data/batch_2026_01/sample_manifest.csv")
    purities = [r["tumor_purity"] for r in rows if r["tumor_purity"] is not None]
    assert all(0.0 <= p <= 1.0 for p in purities)


def test_clean_manifest_batch1_all_dates_parseable():
    rows = clean_manifest("data/batch_2026_01/sample_manifest.csv")
    # S-0012 has empty date; all others should parse
    dated = [r for r in rows if r["collection_date"] is not None]
    assert len(dated) >= 12


def test_clean_manifest_batch2_has_library_prep():
    rows = clean_manifest("data/batch_2026_02/sample_manifest.csv")
    assert all(r.get("library_prep") is not None for r in rows)


def test_clean_manifest_tissue_lowercased():
    rows = clean_manifest("data/batch_2026_01/sample_manifest.csv")
    for r in rows:
        if r["tissue"]:
            assert r["tissue"] == r["tissue"].lower(), f"tissue not lowercase: {r['tissue']!r}"
            assert r["tissue"] == r["tissue"].strip(), f"tissue has whitespace: {r['tissue']!r}"
