import sys
from pathlib import Path
import pandas as pd

# Ensure project root and workflow/scripts are importable
ROOT_DIR = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = ROOT_DIR / "workflow" / "scripts"
for p in (str(ROOT_DIR), str(SCRIPTS_DIR)):
    if p not in sys.path:
        sys.path.insert(0, p)

import fetch_macro
import macro_dashboard


def test_clean_accession_preserves_bvbrc_numeric_ids():
    """Ensures standard INSDC accessions strip versions while BV-BRC numeric IDs stay intact."""
    assert fetch_macro.clean_accession("NC_045512.2") == "NC_045512"
    assert fetch_macro.clean_accession("ENA|OR287049|OR287049.1") == "OR287049"
    assert fetch_macro.clean_accession("OZ200095.1") == "OZ200095"
    assert fetch_macro.clean_accession("2697049.140730") == "2697049.140730"
    assert fetch_macro.clean_accession("2697049.140730.con.0001") == "2697049.140730.con.0001"

    assert fetch_macro.is_insdc_accession("OZ200095") is True
    assert fetch_macro.is_insdc_accession("NC_045512") is True
    assert fetch_macro.is_insdc_accession("2697049.140730") is False


def test_sanitize_date_rejects_epoch_zero_and_corrupt_years():
    """Rejects Unix epoch-0 sentinel dates (1969-12-31 / 1970-01-01) while keeping valid historical dates."""
    assert fetch_macro.sanitize_date("1969-12-31T18:00:00Z", is_db_timestamp=False) == ""
    assert fetch_macro.sanitize_date("1969-12-31", is_db_timestamp=False) == ""
    assert fetch_macro.sanitize_date("1970-01-01", is_db_timestamp=True) == ""
    assert fetch_macro.sanitize_date("1900-01-01", is_db_timestamp=False) == ""
    assert fetch_macro.sanitize_date("2099-05-01", is_db_timestamp=False) == ""
    assert fetch_macro.sanitize_date("1975-06-01", is_db_timestamp=True) == ""

    # Valid biological collection dates (including historical pathogens) must survive
    assert fetch_macro.sanitize_date("1918-10-15", is_db_timestamp=False) == "1918-10-15"
    assert fetch_macro.sanitize_date("2025-04-15T00:00:00Z", is_db_timestamp=False) == "2025-04-15"


def test_parse_spatiotemporal_normalization_and_iso3():
    """Verifies alias mapping, ISO-3 codes, plural /humans/ title parsing, and sentinel date fallback."""
    assert fetch_macro.parse_spatiotemporal_fields(
        "SARS-CoV-2/human/USA/CA-123/2026",
        raw_country="United States of America: California",
        raw_date="2026-08-15",
        fallback_date="2026-09-01",
    ) == ("USA", "2026-08-15", "USA")

    assert fetch_macro.parse_spatiotemporal_fields(
        "Severe acute respiratory syndrome coronavirus 2 SARS-CoV-2/humans/USA/WA-PHL-023800/2022",
        raw_country="",
        raw_date="2022-06-30",
        fallback_date="",
    ) == ("USA", "2022-06-30", "USA")

    assert fetch_macro.parse_spatiotemporal_fields(
        "monopartite assembly",
        raw_country="United Kingdom: England",
        raw_date="1969-12-31",
        fallback_date="2026-10-12",
    ) == ("United Kingdom", "2026-10-12", "GBR")

    assert fetch_macro.parse_spatiotemporal_fields(
        "monopartite assembly",
        raw_country="Turkiye: Istanbul",
        raw_date="2025-03-10",
        fallback_date="",
    ) == ("Turkey", "2025-03-10", "TUR")

    # Diamond Princess / Kanagawa DDBJ isolate maps to Japan; unknown 3-letter codes stay Unknown
    assert fetch_macro.parse_spatiotemporal_fields(
        "Severe acute respiratory syndrome coronavirus 2 SARS-CoV-2/Hu/DP/Kng/19-027 RNA, complete genome",
        raw_country="",
        raw_date="2020-02-10",
        fallback_date="",
    ) == ("Japan", "2020-02-10", "JPN")

    assert fetch_macro.parse_spatiotemporal_fields(
        "Severe acute respiratory syndrome coronavirus 2 SARS-CoV-2/Hu/LAB/XYZ/2020, complete genome",
        raw_country="",
        raw_date="2020-02-10",
        fallback_date="",
    ) == ("Unknown", "2020-02-10", "")


def test_enrich_missing_via_ncbi_mocked(monkeypatch):
    """Mocks NCBI esummary JSON fallback to confirm Unknown INSDC rows resolve while numeric IDs are skipped."""
    captured_ids = []

    class DummyResp:
        status_code = 200

        def json(self):
            return {
                "result": {
                    "uids": ["1001"],
                    "1001": {
                        "caption": "OZ200095",
                        "title": "Severe acute respiratory syndrome coronavirus 2 genome assembly",
                        "subtype": "geo_loc_name|collection_date",
                        "subname": "Denmark|2024-08-26",
                    },
                }
            }

    def fake_post(url, data=None, timeout=None):
        captured_ids.append(data["id"])
        return DummyResp()

    monkeypatch.setattr(fetch_macro.requests, "post", fake_post)

    rows = [
        {"Accession": "OZ200095", "Country": "Unknown", "ISO3": "", "Collection_Date": "Unknown"},
        {"Accession": "2697049.140730", "Country": "Unknown", "ISO3": "", "Collection_Date": "2021-02-13"},
        {"Accession": "MW565837", "Country": "USA", "ISO3": "USA", "Collection_Date": "2021-01-20"},
    ]
    enriched = fetch_macro.enrich_missing_via_ncbi(rows, email="test@example.com")

    assert captured_ids == ["OZ200095"]
    assert enriched[0]["Country"] == "Denmark"
    assert enriched[0]["ISO3"] == "DNK"
    assert enriched[0]["Collection_Date"] == "2024-08-26"
    assert enriched[1]["Country"] == "Unknown"
    assert enriched[2]["Country"] == "USA"


def test_parse_row_selection_ranges_and_edge_cases():
    """Tests inclusive ranges, reversed ranges, deduplication, 'all', clamping, and invalid tokens."""
    assert macro_dashboard.parse_row_selection("0-5, 10", 100) == ([0, 1, 2, 3, 4, 5, 10], [])
    assert macro_dashboard.parse_row_selection("0-2, 1-4, 2", 10) == ([0, 1, 2, 3, 4], [])
    assert macro_dashboard.parse_row_selection("5-2", 10) == ([2, 3, 4, 5], [])
    assert macro_dashboard.parse_row_selection("all", 4) == ([0, 1, 2, 3], [])
    assert macro_dashboard.parse_row_selection("", 10) == ([], [])

    idx_clamped, warns = macro_dashboard.parse_row_selection("0-2, 8-20, bad_token, 50", 10)
    assert idx_clamped == [0, 1, 2, 8, 9]
    assert len(warns) == 3


def test_build_choropleth_figure_scales_and_empty_handling():
    """Verifies Log10/Linear Choropleth generation and None return when no valid ISO-3 rows exist."""
    df = pd.DataFrame([
        {"Country": "USA", "ISO3": "USA"},
        {"Country": "USA", "ISO3": "USA"},
        {"Country": "Denmark", "ISO3": "DNK"},
        {"Country": "Unknown", "ISO3": ""},
    ])
    fig_log = macro_dashboard.build_choropleth_figure(df, use_log_scale=True)
    fig_lin = macro_dashboard.build_choropleth_figure(df, use_log_scale=False)
    assert fig_log is not None and len(fig_log.data) == 1
    assert fig_lin is not None and len(fig_lin.data) == 1

    unknown_only = pd.DataFrame([{"Country": "Unknown", "ISO3": ""}])
    assert macro_dashboard.build_choropleth_figure(unknown_only, use_log_scale=True) is None


def test_flush_macro_results_isolation(tmp_path, monkeypatch):
    """Confirms Macro flush deletes only macro CSV + filter keys, preserving sidebar inputs and Deep outputs."""
    macro_csv = tmp_path / "macro" / "macro_metadata.csv"
    deep_csv = tmp_path / "summary" / "sample_report.csv"
    macro_csv.parent.mkdir(parents=True)
    deep_csv.parent.mkdir(parents=True)
    macro_csv.write_text("Accession\nOZ200095\n")
    deep_csv.write_text("Accession\nNC_045512\n")

    fake_state = {
        "macro_filter_countries": ["Denmark"],
        "macro_row_selection": "0-5",
        "macro_last_logs": "done",
        "macro_db": "BV-BRC",
    }
    monkeypatch.setattr(macro_dashboard, "MACRO_RESULTS_PATH", str(macro_csv))
    monkeypatch.setattr(macro_dashboard.st, "session_state", fake_state)

    assert macro_dashboard.flush_macro_results() is True
    assert not macro_csv.exists()
    assert deep_csv.exists()
    assert "macro_filter_countries" not in fake_state
    assert "macro_row_selection" not in fake_state
    assert "macro_last_logs" not in fake_state
    assert fake_state["macro_db"] == "BV-BRC"
    assert macro_dashboard.flush_macro_results() is False