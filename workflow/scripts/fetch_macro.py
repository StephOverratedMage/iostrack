import os
import io
import re
import sys
import time
import yaml
import random
import calendar
import datetime
import requests
import pandas as pd
from Bio import Entrez

NCBI_ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
ENA_PORTAL_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_PORTAL_COUNT = "https://www.ebi.ac.uk/ena/portal/api/count"
ENA_FIELDS = "accession,sample_accession,description,base_count,first_public,last_updated,country,collection_date"
BVBRC_GENOME_URL = "https://www.bv-brc.org/api/genome/"
BVBRC_HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/x-www-form-urlencoded",
}
BVBRC_FIELDS = (
    "select(genome_id,genbank_accessions,genome_name,genome_length,"
    "isolation_country,geographic_location,collection_date,collection_year,date_inserted,date_modified)"
)

MISSING_TOKENS = {
    "", "nan", "none", "null", "unknown", "missing",
    "not collected", "not applicable", "not provided", "na", "n/a",
}
SENTINEL_DATES = {"1969-12-31", "1970-01-01", "1900-01-01", "0000-00-00", "0001-01-01"}

COUNTRY_ALIASES = {
    "USA": "USA", "UNITED STATES": "USA", "UNITED STATES OF AMERICA": "USA",
    "US": "USA", "U.S.A.": "USA", "U.S.": "USA",
    "CHN": "China", "PRC": "China", "CHINA": "China", "PEOPLE'S REPUBLIC OF CHINA": "China",
    "CHL": "Chile", "CHILE": "Chile",
    "GBR": "United Kingdom", "UK": "United Kingdom", "UNITED KINGDOM": "United Kingdom",
    "GREAT BRITAIN": "United Kingdom", "ENGLAND": "United Kingdom",
    "SCOTLAND": "United Kingdom", "WALES": "United Kingdom", "NORTHERN IRELAND": "United Kingdom",
    "DEU": "Germany", "GERMANY": "Germany", "DEUTSCHLAND": "Germany",
    "FRA": "France", "FRANCE": "France", "ITA": "Italy", "ITALY": "Italy",
    "ESP": "Spain", "SPAIN": "Spain", "CAN": "Canada", "CANADA": "Canada",
    "AUS": "Australia", "AUSTRALIA": "Australia", "JPN": "Japan", "JAPAN": "Japan",
    "KOR": "South Korea", "SOUTH KOREA": "South Korea", "REPUBLIC OF KOREA": "South Korea",
    "KOREA, SOUTH": "South Korea", "KOREA": "South Korea",
    "BRA": "Brazil", "BRAZIL": "Brazil", "BRASIL": "Brazil",
    "IND": "India", "INDIA": "India", "ZAF": "South Africa", "SOUTH AFRICA": "South Africa",
    "MEX": "Mexico", "MEXICO": "Mexico",
    "RUS": "Russia", "RUSSIA": "Russia", "RUSSIAN FEDERATION": "Russia",
    "ARG": "Argentina", "COL": "Colombia", "PER": "Peru",
    "NLD": "Netherlands", "THE NETHERLANDS": "Netherlands", "HOLLAND": "Netherlands",
    "CHE": "Switzerland", "SWE": "Sweden", "BEL": "Belgium", "AUT": "Austria",
    "POL": "Poland", "DNK": "Denmark", "NOR": "Norway", "FIN": "Finland",
    "PRT": "Portugal", "GRC": "Greece", "TUR": "Turkey", "TURKIYE": "Turkey",
    "ISR": "Israel", "EGY": "Egypt", "NGA": "Nigeria", "KEN": "Kenya",
    "THA": "Thailand", "VNM": "Vietnam", "VIET NAM": "Vietnam",
    "SGP": "Singapore", "MYS": "Malaysia", "PHL": "Philippines",
    "IDN": "Indonesia", "NZL": "New Zealand", "PAK": "Pakistan", "BGD": "Bangladesh",
    "SAU": "Saudi Arabia", "ARE": "UAE", "UNITED ARAB EMIRATES": "UAE",
    "CZE": "Czechia", "CZECH REPUBLIC": "Czechia", "CZECHIA": "Czechia",
    "TWN": "Taiwan", "TAIWAN": "Taiwan", "TAIWAN, PROVINCE OF CHINA": "Taiwan",
    "HKG": "Hong Kong", "HONG KONG": "Hong Kong",
    "IRN": "Iran", "IRAN": "Iran", "IRAN, ISLAMIC REPUBLIC OF": "Iran",
    "IRL": "Ireland", "ROU": "Romania", "HUN": "Hungary", "BGR": "Bulgaria",
    "HRV": "Croatia", "SRB": "Serbia", "SVK": "Slovakia", "SVN": "Slovenia",
    "UKR": "Ukraine", "LTU": "Lithuania", "LVA": "Latvia", "EST": "Estonia",
    "ISL": "Iceland", "LUX": "Luxembourg", "CYP": "Cyprus", "MLT": "Malta",
    "MAR": "Morocco", "TUN": "Tunisia", "DZA": "Algeria", "ETH": "Ethiopia",
    "GHA": "Ghana", "UGA": "Uganda", "SEN": "Senegal", "CMR": "Cameroon",
    "COD": "Democratic Republic of the Congo", "DEMOCRATIC REPUBLIC OF THE CONGO": "Democratic Republic of the Congo",
    "DRC": "Democratic Republic of the Congo",
    "CIV": "Cote d'Ivoire", "COTE D'IVOIRE": "Cote d'Ivoire", "IVORY COAST": "Cote d'Ivoire",
    "ECU": "Ecuador", "URY": "Uruguay", "PRY": "Paraguay", "BOL": "Bolivia",
    "VEN": "Venezuela", "CRI": "Costa Rica", "PAN": "Panama", "GTM": "Guatemala",
    "DOM": "Dominican Republic", "PRI": "Puerto Rico", "JAM": "Jamaica",
    "LKA": "Sri Lanka", "NPL": "Nepal", "KHM": "Cambodia", "LAO": "Laos",
    "MMR": "Myanmar", "KAZ": "Kazakhstan", "GEO": "Georgia", "JOR": "Jordan",
    "LBN": "Lebanon", "QAT": "Qatar", "KWT": "Kuwait", "OMN": "Oman", "BHR": "Bahrain",
}

REGION_TO_COUNTRY = {
    "WASHINGTON": "USA", "CALIFORNIA": "USA", "NEW YORK": "USA", "TEXAS": "USA",
    "FLORIDA": "USA", "MASSACHUSETTS": "USA", "SEATTLE": "USA", "HAWAII": "USA",
    "WISCONSIN": "USA", "MINNESOTA": "USA", "PENNSYLVANIA": "USA", "MARYLAND": "USA",
    "SOUTH CAROLINA": "USA", "UTAH": "USA", "ILLINOIS": "USA", "MICHIGAN": "USA",
    "WUHAN": "China", "BEIJING": "China", "SHANGHAI": "China", "GUANGDONG": "China",
    "LONDON": "United Kingdom", "ONTARIO": "Canada", "BRITISH COLUMBIA": "Canada",
    "QUEBEC": "Canada", "VICTORIA": "Australia", "NEW SOUTH WALES": "Australia",
    "SYDNEY": "Australia", "TOKYO": "Japan", "OSAKA": "Japan", "SAO PAULO": "Brazil",
}

COUNTRY_TO_ISO3 = {
    "Afghanistan": "AFG", "Albania": "ALB", "Algeria": "DZA", "Andorra": "AND",
    "Angola": "AGO", "Argentina": "ARG", "Armenia": "ARM", "Australia": "AUS",
    "Austria": "AUT", "Azerbaijan": "AZE", "Bahamas": "BHS", "Bahrain": "BHR",
    "Bangladesh": "BGD", "Barbados": "BRB", "Belarus": "BLR", "Belgium": "BEL",
    "Belize": "BLZ", "Benin": "BEN", "Bhutan": "BTN", "Bolivia": "BOL",
    "Bosnia And Herzegovina": "BIH", "Botswana": "BWA", "Brazil": "BRA",
    "Brunei": "BRN", "Bulgaria": "BGR", "Burkina Faso": "BFA", "Burundi": "BDI",
    "Cambodia": "KHM", "Cameroon": "CMR", "Canada": "CAN", "Cape Verde": "CPV",
    "Central African Republic": "CAF", "Chad": "TCD", "Chile": "CHL", "China": "CHN",
    "Colombia": "COL", "Comoros": "COM", "Congo": "COG",
    "Democratic Republic of the Congo": "COD", "Costa Rica": "CRI", "Croatia": "HRV",
    "Cuba": "CUB", "Cyprus": "CYP", "Czechia": "CZE", "Denmark": "DNK",
    "Djibouti": "DJI", "Dominican Republic": "DOM", "Ecuador": "ECU", "Egypt": "EGY",
    "El Salvador": "SLV", "Equatorial Guinea": "GNQ", "Eritrea": "ERI", "Estonia": "EST",
    "Eswatini": "SWZ", "Ethiopia": "ETH", "Fiji": "FJI", "Finland": "FIN",
    "France": "FRA", "Gabon": "GAB", "Gambia": "GMB", "Georgia": "GEO",
    "Germany": "DEU", "Ghana": "GHA", "Greece": "GRC", "Guatemala": "GTM",
    "Guinea": "GIN", "Guinea-Bissau": "GNB", "Guyana": "GUY", "Haiti": "HTI",
    "Honduras": "HND", "Hong Kong": "HKG", "Hungary": "HUN", "Iceland": "ISL",
    "India": "IND", "Indonesia": "IDN", "Iran": "IRN", "Iraq": "IRQ",
    "Ireland": "IRL", "Israel": "ISR", "Italy": "ITA", "Cote d'Ivoire": "CIV",
    "Jamaica": "JAM", "Japan": "JPN", "Jordan": "JOR", "Kazakhstan": "KAZ",
    "Kenya": "KEN", "Kuwait": "KWT", "Kyrgyzstan": "KGZ", "Laos": "LAO",
    "Latvia": "LVA", "Lebanon": "LBN", "Lesotho": "LSO", "Liberia": "LBR",
    "Libya": "LBY", "Lithuania": "LTU", "Luxembourg": "LUX", "Madagascar": "MDG",
    "Malawi": "MWI", "Malaysia": "MYS", "Maldives": "MDV", "Mali": "MLI",
    "Malta": "MLT", "Mauritania": "MRT", "Mauritius": "MUS", "Mexico": "MEX",
    "Moldova": "MDA", "Monaco": "MCO", "Mongolia": "MNG", "Montenegro": "MNE",
    "Morocco": "MAR", "Mozambique": "MOZ", "Myanmar": "MMR", "Namibia": "NAM",
    "Nepal": "NPL", "Netherlands": "NLD", "New Zealand": "NZL", "Nicaragua": "NIC",
    "Niger": "NER", "Nigeria": "NGA", "North Korea": "PRK", "North Macedonia": "MKD",
    "Norway": "NOR", "Oman": "OMN", "Pakistan": "PAK", "Palestine": "PSE",
    "Panama": "PAN", "Papua New Guinea": "PNG", "Paraguay": "PRY", "Peru": "PER",
    "Philippines": "PHL", "Poland": "POL", "Portugal": "PRT", "Puerto Rico": "PRI",
    "Qatar": "QAT", "Romania": "ROU", "Russia": "RUS", "Rwanda": "RWA",
    "Saudi Arabia": "SAU", "Senegal": "SEN", "Serbia": "SRB", "Sierra Leone": "SLE",
    "Singapore": "SGP", "Slovakia": "SVK", "Slovenia": "SVN", "Somalia": "SOM",
    "South Africa": "ZAF", "South Korea": "KOR", "South Sudan": "SSD", "Spain": "ESP",
    "Sri Lanka": "LKA", "Sudan": "SDN", "Suriname": "SUR", "Sweden": "SWE",
    "Switzerland": "CHE", "Syria": "SYR", "Taiwan": "TWN", "Tajikistan": "TJK",
    "Tanzania": "TZA", "Thailand": "THA", "Timor-Leste": "TLS", "Togo": "TGO",
    "Trinidad And Tobago": "TTO", "Tunisia": "TUN", "Turkey": "TUR",
    "Turkmenistan": "TKM", "UAE": "ARE", "Uganda": "UGA", "Ukraine": "UKR",
    "United Kingdom": "GBR", "USA": "USA", "Uruguay": "URY", "Uzbekistan": "UZB",
    "Venezuela": "VEN", "Vietnam": "VNM", "Yemen": "YEM", "Zambia": "ZMB",
    "Zimbabwe": "ZWE",
}
COUNTRY_TO_ISO3_UPPER = {k.upper(): v for k, v in COUNTRY_TO_ISO3.items()}


def emit_progress(pct, message):
    """Emits a structured, unbuffered progress line consumed by macro_dashboard.py."""
    clamped = max(0, min(100, int(round(pct))))
    clean_msg = str(message).replace("\n", " ").strip()
    print(f"PROGRESS|{clamped}|{clean_msg}", flush=True)


def clean_accession(raw_id):
    """Extracts clean base accession while preserving BV-BRC numeric genome IDs."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    if match:
        return match.group(1)
    if re.match(r"^\d+\.\d+", cleaned):
        return cleaned
    return cleaned.split(".")[0]


def is_insdc_accession(acc):
    """Returns True if acc is a standard INSDC/GenBank/ENA accession (not a numeric BV-BRC genome_id)."""
    return bool(re.fullmatch(r"[A-Z]{1,6}_?\d{4,}", str(acc or "").strip()))


def clean_meta_str(val):
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    text = str(val).strip()
    return "" if text.lower() in MISSING_TOKENS else text


def sanitize_date(val, is_db_timestamp=False):
    cleaned = clean_meta_str(val)
    if not cleaned:
        return ""
    date_part = cleaned.split("T")[0] if "T" in cleaned else cleaned
    norm = date_part.replace("/", "-")
    if norm in SENTINEL_DATES:
        return ""
    m_yr = re.match(r"^(\d{4})", norm)
    if m_yr:
        yr = int(m_yr.group(1))
        min_yr = 1980 if is_db_timestamp else 1900
        if yr < min_yr or yr > datetime.date.today().year + 1:
            return ""
    return date_part


def normalize_country(raw_val):
    cleaned = clean_meta_str(raw_val)
    if not cleaned:
        return "Unknown"
    base = cleaned.split(":")[0].strip(" ,;()[]")
    if not base or base.lower() in MISSING_TOKENS:
        return "Unknown"
    upper = base.upper()
    if upper in COUNTRY_ALIASES:
        return COUNTRY_ALIASES[upper]
    if upper in REGION_TO_COUNTRY:
        return REGION_TO_COUNTRY[upper]
    return base if upper in ("USA", "UK", "UAE") else base.title()


def country_to_iso3(country_name):
    if not country_name or country_name == "Unknown":
        return ""
    return COUNTRY_TO_ISO3_UPPER.get(str(country_name).strip().upper(), "")


def parse_spatiotemporal_fields(title, raw_country="", raw_date="", fallback_date=""):
    country = normalize_country(raw_country)

    if country == "Unknown" and "/DP/KNG/" in str(title).upper():
        country = "Japan"

    if country == "Unknown" and "/" in str(title):
        slash_tokens = [tok.strip(" ,;()[]") for tok in str(title).split("/")[1:]]
        skip_terms = {
            "HUMAN", "HUMANS", "HU", "HOMO SAPIENS", "SWINE", "BAT", "AVIAN", "CAT", "DOG",
            "MINK", "DEER", "RODENT", "MONKEY", "CAMEL", "MAMMAL", "ENV",
            "ENVIRONMENT", "WASTEWATER", "ISOLATE", "STRAIN", "COMPLETE",
            "GENOME", "SARS-COV-2", "HCOV-19",
        }
        for tok in slash_tokens:
            cand = tok.upper()
            if not cand or cand in skip_terms or cand.isdigit():
                continue
            if cand in COUNTRY_ALIASES or cand in REGION_TO_COUNTRY or cand in COUNTRY_TO_ISO3_UPPER:
                country = normalize_country(cand)
                break

    if country == "Unknown":
        title_upper = str(title).upper()
        for code, mapped in COUNTRY_ALIASES.items():
            if len(code) > 3 and code in title_upper:
                country = mapped
                break
        if country == "Unknown":
            for reg, mapped in REGION_TO_COUNTRY.items():
                if reg in title_upper:
                    country = mapped
                    break

    col_date = sanitize_date(raw_date, is_db_timestamp=False)
    if not col_date:
        m_ymd = re.search(r"/(20[0-2][0-9])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]*/", str(title))
        if m_ymd:
            col_date = f"{m_ymd.group(1)}-{m_ymd.group(2)}-{m_ymd.group(3)}"
        else:
            m_yymmdd = re.search(r"-(2[0-6])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]{2,}/20\1", str(title))
            if m_yymmdd:
                col_date = f"20{m_yymmdd.group(1)}-{m_yymmdd.group(2)}-{m_yymmdd.group(3)}"

    if not col_date:
        col_date = sanitize_date(fallback_date, is_db_timestamp=True)

    return country, col_date or "Unknown", country_to_iso3(country)


def resolve_taxid(organism_input):
    """Resolves an organism name or numeric TaxID to an NCBI Taxonomy ID string."""
    text = str(organism_input or "").strip()
    if not text:
        return None
    if text.isdigit():
        return text
    m_tx = re.search(r"(?:tax_tree|taxon_id|taxon_lineage_ids)\s*[\(,]\s*(\d+)", text)
    if m_tx:
        return m_tx.group(1)
    clean_org = re.sub(r"\[Organism[^\]]*\]", "", text, flags=re.IGNORECASE).strip(' "\'')
    try:
        with Entrez.esearch(db="taxonomy", term=f'"{clean_org}"[Scientific Name]', retmax=3) as h:
            res = Entrez.read(h)
        if res.get("IdList"):
            return str(res["IdList"][0])
        with Entrez.esearch(db="taxonomy", term=clean_org, retmax=3) as h:
            res = Entrez.read(h)
        if res.get("IdList"):
            return str(res["IdList"][0])
    except Exception as e:
        print(f"⚠️ Warning: Taxonomy lookup failed for '{clean_org}' ({e})", flush=True)
    return None


def interleave_lists(buckets):
    out, seen = [], set()
    max_len = max((len(b) for b in buckets), default=0)
    for idx in range(max_len):
        for b in buckets:
            if idx < len(b) and b[idx] not in seen:
                seen.add(b[idx])
                out.append(b[idx])
    return out


def enrich_missing_via_ncbi(rows, email, max_enrich=1000):
    """
    Batched NCBI esummary JSON fallback for BV-BRC/ENA rows that have
    Country=='Unknown' or Collection_Date=='Unknown' but carry a standard INSDC accession.
    """
    need_accs = []
    seen = set()
    for r in rows:
        if r.get("Country") == "Unknown" or r.get("Collection_Date") == "Unknown" or r.get("_date_is_fallback"):
            acc = str(r.get("Accession", "")).strip()
            if is_insdc_accession(acc) and acc not in seen:
                seen.add(acc)
                need_accs.append(acc)
                if len(need_accs) >= max_enrich:
                    break

    if not need_accs:
        return rows

    batch_size = 200
    total_batches = max(1, (len(need_accs) + batch_size - 1) // batch_size)
    emit_progress(80, f"Enriching {len(need_accs):,} incomplete record(s) via NCBI fallback (0/{total_batches} batches)...")

    ncbi_lookup = {}
    for b_idx, i in enumerate(range(0, len(need_accs), batch_size), start=1):
        chunk = need_accs[i : i + batch_size]
        payload = {
            "db": "nucleotide",
            "id": ",".join(chunk),
            "retmode": "json",
            "email": email,
        }
        try:
            resp = requests.post(NCBI_ESUMMARY_URL, data=payload, timeout=(10, 30))
            if resp.status_code == 200:
                result_dict = resp.json().get("result", {})
                for uid in result_dict.get("uids", []):
                    s = result_dict.get(str(uid), {})
                    if not isinstance(s, dict):
                        continue
                    acc = clean_accession(s.get("caption") or s.get("accessionversion") or "")
                    if not acc:
                        continue
                    sub_types = str(s.get("subtype", "")).split("|")
                    sub_names = str(s.get("subname", "")).split("|")
                    sub_dict = {
                        k.strip().lower(): v.strip()
                        for k, v in zip(sub_types, sub_names)
                        if k.strip() and v.strip()
                    }
                    raw_country = (
                        sub_dict.get("geo_loc_name")
                        or sub_dict.get("geo-loc-name")
                        or sub_dict.get("country")
                        or ""
                    )
                    raw_date = sub_dict.get("collection_date") or sub_dict.get("collection-date") or ""
                    c, d, iso = parse_spatiotemporal_fields(
                        title=s.get("title", ""),
                        raw_country=raw_country,
                        raw_date=raw_date,
                        fallback_date="",
                    )
                    ncbi_lookup[acc] = (c, d, iso)
        except Exception as e:
            print(f"⚠️ Warning: NCBI fallback batch skipped ({e})", flush=True)

        pct = 80 + int(15 * (b_idx / total_batches))
        emit_progress(pct, f"NCBI enrichment fallback: batch {b_idx}/{total_batches} complete...")
        if i + batch_size < len(need_accs):
            time.sleep(0.35)

    resolved_geo = 0
    resolved_date = 0
    for r in rows:
        acc = str(r.get("Accession", "")).strip()
        if acc not in ncbi_lookup:
            continue
        c, d, iso = ncbi_lookup[acc]
        if r.get("Country") == "Unknown" and c != "Unknown":
            r["Country"] = c
            r["ISO3"] = iso
            resolved_geo += 1
        if r.get("Collection_Date") == "Unknown" and d != "Unknown":
            r["Collection_Date"] = d
            resolved_date += 1

    print(f"✅ NCBI fallback resolved {resolved_geo} country tag(s) and {resolved_date} date(s).", flush=True)
    return rows


# ---------------------- NCBI MACRO FETCHER ----------------------

def ncbi_esearch(term, retmax=500, retstart=0, sort="Date Released", mindate=None, maxdate=None, count_only=False):
    kwargs = {"db": "nucleotide", "term": term}
    if count_only:
        kwargs["rettype"] = "count"
    else:
        kwargs["retmax"] = int(retmax)
        kwargs["retstart"] = int(retstart)
        if sort:
            kwargs["sort"] = sort
    if mindate and maxdate:
        kwargs["mindate"] = mindate
        kwargs["maxdate"] = maxdate
        kwargs["datetype"] = "pdat"
    try:
        with Entrez.esearch(**kwargs) as handle:
            rec = Entrez.read(handle)
        return int(rec.get("Count", 0)), list(rec.get("IdList", []))
    except Exception as e:
        print(f"⚠️ Warning: NCBI esearch sub-query skipped ({e})", flush=True)
        return 0, []


def find_ncbi_min_year(base_query, current_year, floor_year=1980):
    total, _ = ncbi_esearch(base_query, count_only=True)
    if total == 0:
        return current_year, 0
    lo, hi = floor_year, current_year
    while lo < hi:
        mid = (lo + hi) // 2
        cnt, _ = ncbi_esearch(
            base_query, mindate=f"{floor_year}/01/01", maxdate=f"{mid}/12/31", count_only=True
        )
        if cnt > 0:
            hi = mid
        else:
            lo = mid + 1
    return lo, total


def fetch_macro_ncbi(organism_input, target_quota, min_len, max_len, strategy, random_seed, email):
    text = str(organism_input).strip()
    if text.isdigit():
        ncbi_query = f"txid{text}[Organism:exp]"
    elif "[" in text:
        ncbi_query = text
    else:
        ncbi_query = f'"{text}"[Organism]'

    base_query = f"({ncbi_query}) AND {min_len}:{max_len}[SLEN]"
    today = datetime.date.today()
    current_year, current_month = today.year, today.month
    uids = []

    if strategy == "Newest First (Deterministic)":
        emit_progress(10, f"Querying NCBI for newest {target_quota:,} accessions...")
        _, uids = ncbi_esearch(base_query, retmax=target_quota, sort="Date Released")
        emit_progress(40, f"Found {len(uids):,} NCBI UID(s). Starting metadata download...")

    elif strategy == "Longest Sequences First":
        emit_progress(10, f"Querying NCBI for longest {target_quota:,} accessions...")
        _, uids = ncbi_esearch(base_query, retmax=target_quota, sort="Sequence Length")
        emit_progress(40, f"Found {len(uids):,} NCBI UID(s). Starting metadata download...")

    elif strategy == "Oldest First (Deterministic)":
        emit_progress(8, "Probing earliest publication year in NCBI...")
        min_year, total = find_ncbi_min_year(base_query, current_year)
        emit_progress(15, f"Scanning NCBI chronologically from {min_year} ({total:,} total hits)...")
        seen = set()
        for yr in range(min_year, current_year + 1):
            yr_cnt, _ = ncbi_esearch(base_query, mindate=f"{yr}/01/01", maxdate=f"{yr}/12/31", count_only=True)
            if yr_cnt == 0:
                continue
            if yr_cnt <= 2000:
                _, yr_ids = ncbi_esearch(
                    base_query, retmax=2000, sort="Date Released",
                    mindate=f"{yr}/01/01", maxdate=f"{yr}/12/31"
                )
                for u in reversed(yr_ids):
                    if u not in seen:
                        seen.add(u)
                        uids.append(u)
                        if len(uids) >= target_quota:
                            break
            else:
                max_m = current_month if yr == current_year else 12
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    m_start, m_end = f"{yr}/{m:02d}/01", f"{yr}/{m:02d}/{last_day:02d}"
                    m_cnt, _ = ncbi_esearch(base_query, mindate=m_start, maxdate=m_end, count_only=True)
                    if m_cnt == 0:
                        continue
                    need = min(500, target_quota - len(uids) + 50)
                    retstart = min(max(0, m_cnt - need), 9500)
                    _, m_ids = ncbi_esearch(
                        base_query, retmax=need, retstart=retstart, sort="Date Released",
                        mindate=m_start, maxdate=m_end
                    )
                    for u in reversed(m_ids):
                        if u not in seen:
                            seen.add(u)
                            uids.append(u)
                            if len(uids) >= target_quota:
                                break
                    pct = 15 + int(25 * min(1.0, len(uids) / max(1, target_quota)))
                    emit_progress(pct, f"Collected {len(uids):,}/{target_quota:,} oldest UIDs (through {yr}-{m:02d})...")
                    if len(uids) >= target_quota:
                        break
            pct = 15 + int(25 * min(1.0, len(uids) / max(1, target_quota)))
            emit_progress(pct, f"Collected {len(uids):,}/{target_quota:,} oldest UIDs (through {yr})...")
            if len(uids) >= target_quota:
                break

    else:
        rng = random.Random(random_seed)
        emit_progress(8, "Probing NCBI archive date span for stratified sampling...")
        min_year, total = find_ncbi_min_year(base_query, current_year)
        emit_progress(14, f"Stratified sampling across {min_year}-{current_year} ({total:,} total hits)...")
        if total <= min(target_quota * 2, 9500):
            _, all_ids = ncbi_esearch(base_query, retmax=min(total, 9500), sort="Date Released")
            all_ids.sort()
            rng.shuffle(all_ids)
            uids = all_ids[:target_quota]
        else:
            ym_windows = [
                (yr, m)
                for yr in range(min_year, current_year + 1)
                for m in range(1, (current_month if yr == current_year else 12) + 1)
            ]
            rng.shuffle(ym_windows)
            n_win = min(len(ym_windows), max(10, min(20, target_quota // 120 + 8)))
            per_win = max(50, int((target_quota * 1.4) // n_win) + 1)
            buckets = []
            for yr, m in ym_windows:
                if len(buckets) >= n_win:
                    break
                last_day = calendar.monthrange(yr, m)[1]
                m_start, m_end = f"{yr}/{m:02d}/01", f"{yr}/{m:02d}/{last_day:02d}"
                m_cnt, _ = ncbi_esearch(base_query, mindate=m_start, maxdate=m_end, count_only=True)
                if m_cnt == 0:
                    continue
                max_start = max(0, min(m_cnt - per_win, 9000))
                retstart = rng.randint(0, max_start) if max_start > 0 else 0
                _, m_ids = ncbi_esearch(
                    base_query, retmax=per_win, retstart=retstart,
                    sort="Date Released", mindate=m_start, maxdate=m_end
                )
                if m_ids:
                    m_ids.sort()
                    rng.shuffle(m_ids)
                    buckets.append(m_ids)
                    pct = 14 + int(26 * (len(buckets) / max(1, n_win)))
                    emit_progress(pct, f"Sampled NCBI time window {len(buckets)}/{n_win} ({yr}-{m:02d})...")
            uids = interleave_lists(buckets)[:target_quota]

    uids = uids[:target_quota]
    if not uids:
        return []

    batch_size = 250
    total_batches = max(1, (len(uids) + batch_size - 1) // batch_size)
    emit_progress(42, f"Downloading NCBI JSON summaries for {len(uids):,} records ({total_batches} batches)...")

    rows = []
    seen_accs = set()
    for b_idx, i in enumerate(range(0, len(uids), batch_size), start=1):
        chunk = uids[i : i + batch_size]
        payload = {
            "db": "nucleotide",
            "id": ",".join(chunk),
            "retmode": "json",
            "email": email,
        }
        raw_json = None
        for _ in range(2):
            try:
                resp = requests.post(NCBI_ESUMMARY_URL, data=payload, timeout=(10, 35))
                if resp.status_code == 200:
                    raw_json = resp.json()
                    break
                time.sleep(0.8)
            except Exception:
                time.sleep(0.8)
        if raw_json:
            result_dict = raw_json.get("result", {})
            for uid in result_dict.get("uids", []):
                s = result_dict.get(str(uid), {})
                if not isinstance(s, dict):
                    continue
                raw_acc = s.get("caption") or s.get("accessionversion") or ""
                clean_acc = clean_accession(raw_acc)
                if not clean_acc or clean_acc in seen_accs:
                    continue
                seen_accs.add(clean_acc)

                title = s.get("title", "")
                update_date = sanitize_date(s.get("updatedate", ""), is_db_timestamp=True)
                create_date = sanitize_date(s.get("createdate", ""), is_db_timestamp=True) or update_date

                sub_types = str(s.get("subtype", "")).split("|")
                sub_names = str(s.get("subname", "")).split("|")
                sub_dict = {
                    k.strip().lower(): v.strip()
                    for k, v in zip(sub_types, sub_names)
                    if k.strip() and v.strip()
                }
                raw_country = sub_dict.get("geo_loc_name") or sub_dict.get("geo-loc-name") or sub_dict.get("country") or ""
                raw_date = sub_dict.get("collection_date") or sub_dict.get("collection-date") or ""
                country, col_date, iso3 = parse_spatiotemporal_fields(title, raw_country, raw_date, create_date)

                rows.append({
                    "Accession": clean_acc,
                    "Title": title,
                    "Length_bp": int(s.get("slen") or 0),
                    "UpdateDate": update_date,
                    "Collection_Date": col_date,
                    "Country": country,
                    "ISO3": iso3,
                    "Data_Source": "NCBI",
                })

        pct = 42 + int(53 * (b_idx / total_batches))
        emit_progress(pct, f"Downloaded NCBI metadata batch {b_idx}/{total_batches} ({len(rows):,} records)...")
        time.sleep(0.35)

    return rows


# ---------------------- ENA MACRO FETCHER ----------------------

def ena_count(query, timeout=(8, 15), retries=2, return_none_on_error=False):
    for attempt in range(retries):
        try:
            resp = requests.get(
                ENA_PORTAL_COUNT,
                params={"result": "sequence", "query": query},
                timeout=timeout,
            )
            if resp.status_code == 200:
                digits = re.search(r"(\d+)", resp.text)
                return int(digits.group(1)) if digits else 0
        except Exception:
            pass
        if attempt + 1 < retries:
            time.sleep(0.6)
    return None if return_none_on_error else 0


def ena_search_df(query, limit=1000, timeout=(8, 25)):
    """Queries ENA Portal /search WITHOUT offset (ENA Portal rejects offset with HTTP 400)."""
    params = {
        "result": "sequence",
        "query": query,
        "limit": int(limit),
        "format": "tsv",
        "fields": ENA_FIELDS,
    }
    try:
        resp = requests.get(ENA_PORTAL_SEARCH, params=params, timeout=timeout)
        if resp.status_code != 200:
            print(f"⚠️ Warning: ENA HTTP {resp.status_code} for query slice: {resp.text[:120].strip()}", flush=True)
            return pd.DataFrame()
        txt = resp.text.strip()
        if not txt:
            return pd.DataFrame()
        df = pd.read_csv(io.StringIO(txt), sep="\t")
        if "accession" in df.columns and not df.empty:
            df["accession"] = df["accession"].astype(str).apply(clean_accession)
            return df.drop_duplicates(subset=["accession"])
    except Exception as e:
        print(f"⚠️ Warning: ENA sub-query skipped ({e})", flush=True)
    return pd.DataFrame()


def find_ena_min_year(base_query, current_year, floor_year=1980):
    """Fast min-year probe (checks 2019 first so post-2019 taxa take <=3 calls)."""
    total = ena_count(base_query)
    if total == 0:
        return current_year, 0

    if current_year > 2019 and ena_count(f"{base_query} AND first_public<=2019-12-31") == 0:
        lo, hi = 2020, current_year
    else:
        lo, hi = floor_year, min(2019, current_year)

    while lo < hi:
        mid = (lo + hi) // 2
        if ena_count(f"{base_query} AND first_public<={mid}-12-31") > 0:
            hi = mid
        else:
            lo = mid + 1
    return lo, total


def interleave_dfs(df_list):
    if not df_list:
        return pd.DataFrame()
    rows, seen = [], set()
    max_rows = max(len(d) for d in df_list)
    for idx in range(max_rows):
        for d in df_list:
            if idx < len(d):
                row = d.iloc[idx]
                acc = row["accession"]
                if acc not in seen:
                    seen.add(acc)
                    rows.append(row)
    return pd.DataFrame(rows).reset_index(drop=True) if rows else pd.DataFrame()


def enrich_ena_via_biosample(rows):
    """
    Batched ENA Portal (result=sample) lookup using linked BioSample accessions (e.g., SAMEA6786471).
    Retrieves true sample-origin country and collection_date when omitted on the sequence record.
    """
    need_samples = sorted({
        r["_sample_acc"]
        for r in rows
        if r.get("_sample_acc") and (
            r.get("Country") == "Unknown"
            or r.get("Collection_Date") == "Unknown"
            or r.get("_date_is_fallback")
        )
    })
    if not need_samples:
        return rows

    batch_size = 50
    total_batches = max(1, (len(need_samples) + batch_size - 1) // batch_size)
    emit_progress(77, f"Checking {len(need_samples):,} linked ENA BioSample(s) for sample origin...")

    sample_lookup = {}
    for b_idx, i in enumerate(range(0, len(need_samples), batch_size), start=1):
        chunk = need_samples[i : i + batch_size]
        q = "(" + " OR ".join(f'sample_accession="{s}"' for s in chunk) + ")"
        payload = {
            "result": "sample",
            "query": q,
            "fields": "sample_accession,country,collection_date",
            "format": "tsv",
            "limit": len(chunk) + 10,
        }
        try:
            resp = requests.post(ENA_PORTAL_SEARCH, data=payload, timeout=(8, 20))
            if resp.status_code == 200 and resp.text.strip():
                sdf = pd.read_csv(io.StringIO(resp.text.strip()), sep="\t")
                if "sample_accession" in sdf.columns:
                    for _, sr in sdf.iterrows():
                        sacc = clean_meta_str(sr.get("sample_accession", ""))
                        if sacc:
                            sample_lookup[sacc] = (
                                sr.get("country", ""),
                                sr.get("collection_date", ""),
                            )
        except Exception as e:
            print(f"⚠️ Warning: ENA BioSample batch {b_idx}/{total_batches} skipped ({e})", flush=True)

    resolved_geo = 0
    resolved_date = 0
    for r in rows:
        sacc = r.get("_sample_acc", "")
        if sacc not in sample_lookup:
            continue
        raw_c, raw_d = sample_lookup[sacc]
        if r.get("Country") == "Unknown":
            c = normalize_country(raw_c)
            if c != "Unknown":
                r["Country"] = c
                r["ISO3"] = country_to_iso3(c)
                resolved_geo += 1
        if r.get("Collection_Date") == "Unknown" or r.get("_date_is_fallback"):
            d = sanitize_date(raw_d, is_db_timestamp=False)
            if d:
                r["Collection_Date"] = d
                r["_date_is_fallback"] = False
                resolved_date += 1

    print(f"✅ ENA BioSample lookup resolved {resolved_geo} country tag(s) and {resolved_date} date(s).", flush=True)
    return rows

def fetch_macro_ena(organism_input, target_quota, min_len, max_len, strategy, random_seed):
    text = str(organism_input).strip()
    if any(tok in text for tok in ("tax_tree(", "tax_eq(", "=", " AND ")):
        ena_query = text
    else:
        emit_progress(6, f"Resolving NCBI TaxID for '{text}'...")
        taxid = resolve_taxid(text)
        if not taxid:
            print(f"❌ Could not resolve NCBI TaxID for ENA query '{text}'.", file=sys.stderr, flush=True)
            return []
        ena_query = f"tax_tree({taxid})"

    base_query = f"{ena_query} AND base_count>={min_len} AND base_count<={max_len}"
    today = datetime.date.today()
    current_year, current_month = today.year, today.month
    meta_df = pd.DataFrame()

    if strategy == "Newest First (Deterministic)":
        windows = (7, 30, 90, 365, None)
        frames = []
        for idx, days_back in enumerate(windows, start=1):
            label = f"last {days_back} days" if days_back is not None else "all dates"
            emit_progress(10 + int(60 * (idx - 1) / len(windows)), f"Querying ENA rolling recency window ({label})...")
            q = (
                f"{base_query} AND last_updated>={(today - datetime.timedelta(days=days_back)).isoformat()}"
                if days_back is not None
                else base_query
            )
            part = ena_search_df(q, limit=target_quota)
            if not part.empty:
                frames.append(part)
                merged = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"])
                if len(merged) >= target_quota:
                    meta_df = merged
                    break
        if meta_df.empty and frames:
            meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"])
        if not meta_df.empty and "last_updated" in meta_df.columns:
            meta_df = meta_df.sort_values(by=["last_updated", "accession"], ascending=[False, True]).head(target_quota)

    elif strategy == "Oldest First (Deterministic)":
        emit_progress(10, "Probing earliest public year in ENA...")
        min_year, total = find_ena_min_year(base_query, current_year)
        emit_progress(18, f"Scanning ENA chronologically from {min_year} ({total:,} total hits)...")
        frames = []
        collected = 0
        t_start = time.monotonic()
        for yr in range(min_year, current_year + 1):
            if time.monotonic() - t_start > 35:
                break
            yr_q = f"{base_query} AND first_public>={yr}-01-01 AND first_public<={yr}-12-31"
            yr_cnt = ena_count(yr_q)
            if yr_cnt == 0:
                continue
            if yr_cnt <= target_quota:
                part = ena_search_df(yr_q, limit=target_quota)
                if not part.empty:
                    frames.append(part)
            else:
                max_m = current_month if yr == current_year else 12
                for m in range(1, max_m + 1):
                    if time.monotonic() - t_start > 35:
                        break
                    last_day = calendar.monthrange(yr, m)[1]
                    m_q = f"{base_query} AND first_public>={yr}-{m:02d}-01 AND first_public<={yr}-{m:02d}-{last_day:02d}"
                    part = ena_search_df(m_q, limit=target_quota)
                    if not part.empty:
                        frames.append(part)
                    if frames:
                        collected = len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]))
                        pct = 18 + int(55 * min(1.0, collected / max(1, target_quota)))
                        emit_progress(pct, f"Collected {collected:,}/{target_quota:,} ENA records (through {yr}-{m:02d})...")
                        if collected >= target_quota:
                            break
            if frames:
                collected = len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]))
                pct = 18 + int(55 * min(1.0, collected / max(1, target_quota)))
                emit_progress(pct, f"Collected {collected:,}/{target_quota:,} ENA records (through {yr})...")
            if collected >= target_quota:
                break
        if frames:
            meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"])
            sort_col = "first_public" if "first_public" in meta_df.columns else "last_updated"
            meta_df = meta_df.sort_values(by=[sort_col, "accession"], ascending=[True, True]).head(target_quota)

    elif strategy == "Longest Sequences First":
        pool_cap = max(target_quota * 5, 1500)
        total_avail = ena_count(base_query, return_none_on_error=True)

        if total_avail is not None and total_avail <= pool_cap:
            emit_progress(55, f"Fetching all {total_avail:,} matching ENA records for exact length ranking...")
            meta_df = ena_search_df(base_query, limit=max(total_avail + 50, target_quota * 2))
        else:
            lo, hi = min_len, max_len
            best_low, best_cnt = min_len, total_avail
            for step in range(1, 23):
                if lo > hi:
                    break
                mid = (lo + hi) // 2
                pct = 10 + min(52, int(52 * (step / 14)))
                emit_progress(pct, f"Converging ENA length threshold (step {step}: >= {mid:,} bp)...")
                cnt = ena_count(
                    f"{ena_query} AND base_count>={mid} AND base_count<={max_len}",
                    return_none_on_error=True,
                )
                if cnt is None:
                    print(f"⚠️ Warning: ENA count timed out at >= {mid:,} bp; using >= {best_low:,} bp.", flush=True)
                    break
                if cnt >= target_quota:
                    best_low = mid
                    best_cnt = cnt
                    if cnt <= pool_cap:
                        break
                    lo = mid + 1
                else:
                    hi = mid - 1

            if best_cnt is not None and best_cnt <= pool_cap:
                emit_progress(65, f"Fetching complete upper tail ({best_cnt:,} ENA records >= {best_low:,} bp)...")
                meta_df = ena_search_df(
                    f"{ena_query} AND base_count>={best_low} AND base_count<={max_len}",
                    limit=max(best_cnt + 50, target_quota * 2),
                )
            else:
                emit_progress(64, f"Fetching strictly longer ENA tail (> {best_low:,} bp) + tied records at {best_low:,} bp...")
                frames = []
                if best_low < max_len:
                    strict_df = ena_search_df(
                        f"{ena_query} AND base_count>={best_low + 1} AND base_count<={max_len}",
                        limit=max(target_quota * 2, 1000),
                    )
                    if not strict_df.empty:
                        frames.append(strict_df)
                tie_df = ena_search_df(
                    f"{ena_query} AND base_count>={best_low} AND base_count<={best_low}",
                    limit=pool_cap,
                )
                if not tie_df.empty:
                    frames.append(tie_df)
                meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]) if frames else pd.DataFrame()

        if not meta_df.empty:
            meta_df["_bc"] = pd.to_numeric(meta_df["base_count"], errors="coerce").fillna(0)
            meta_df = meta_df.sort_values(by=["_bc", "accession"], ascending=[False, True]).drop(columns=["_bc"]).head(target_quota)

    else:
        rng = random.Random(random_seed)
        emit_progress(10, "Probing ENA archive date span...")
        min_year, total = find_ena_min_year(base_query, current_year)
        emit_progress(16, f"Stratified ENA sampling across {min_year}-{current_year} ({total:,} total hits)...")

        if total > 0 and total <= min(target_quota * 3, 4000):
            emit_progress(35, f"Downloading ENA metadata pool ({total:,} records)...")
            meta_df = ena_search_df(base_query, limit=min(total, 4000))
            if not meta_df.empty:
                meta_df = meta_df.sort_values(by=["accession"]).sample(frac=1.0, random_state=random_seed).head(target_quota)
        else:
            use_3day = (total == 0 or total > 50000)
            date_windows = []
            for yr in range(min_year, current_year + 1):
                max_m = current_month if yr == current_year else 12
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    if use_3day:
                        for d_start in range(1, last_day + 1, 3):
                            d_end = min(last_day, d_start + 2)
                            date_windows.append((f"{yr}-{m:02d}-{d_start:02d}", f"{yr}-{m:02d}-{d_end:02d}"))
                    else:
                        date_windows.append((f"{yr}-{m:02d}-01", f"{yr}-{m:02d}-{last_day:02d}"))

            rng.shuffle(date_windows)
            n_target_buckets = min(10, max(6, target_quota // 200 + 5))
            per_win = max(30, int((target_quota * 1.6) // n_target_buckets) + 1)
            buckets = []
            t_start = time.monotonic()
            attempts = 0
            max_attempts = n_target_buckets + 6

            for w_start, w_end in date_windows:
                if len(buckets) >= n_target_buckets or attempts >= max_attempts or (time.monotonic() - t_start) > 35:
                    break
                attempts += 1
                pct = 18 + int(52 * (attempts / max(1, max_attempts)))
                emit_progress(pct, f"Querying ENA date window {attempts}/{max_attempts} ({w_start} to {w_end})...")
                w_q = f"{base_query} AND first_public>={w_start} AND first_public<={w_end}"
                part = ena_search_df(w_q, limit=per_win, timeout=(8, 15))
                if not part.empty:
                    part = part.sort_values(by=["accession"]).sample(frac=1.0, random_state=rng.randint(0, 10**6))
                    buckets.append(part)

            meta_df = interleave_dfs(buckets).head(target_quota)

            if len(meta_df) < target_quota:
                emit_progress(72, "Topping up ENA cohort to target count...")
                fb = ena_search_df(base_query, limit=target_quota)
                if not fb.empty:
                    meta_df = pd.concat([meta_df, fb], ignore_index=True).drop_duplicates(subset=["accession"]).head(target_quota)

    if meta_df.empty:
        return []

    emit_progress(75, f"Normalizing spatiotemporal fields for {len(meta_df):,} ENA records...")
    rows = []
    for _, r in meta_df.iterrows():
        acc = clean_accession(r.get("accession", ""))
        title = clean_meta_str(r.get("description", "")) or acc
        upd = sanitize_date(r.get("last_updated", ""), is_db_timestamp=True) or sanitize_date(r.get("first_public", ""), is_db_timestamp=True)
        fb = sanitize_date(r.get("first_public", ""), is_db_timestamp=True) or upd
        raw_date = sanitize_date(r.get("collection_date", ""), is_db_timestamp=False)
        country, col_date, iso3 = parse_spatiotemporal_fields(title, r.get("country", ""), raw_date, fb)
        rows.append({
            "Accession": acc,
            "Title": title,
            "Length_bp": int(pd.to_numeric(r.get("base_count", 0), errors="coerce") or 0),
            "UpdateDate": upd,
            "Collection_Date": col_date,
            "Country": country,
            "ISO3": iso3,
            "Data_Source": "ENA",
            "_sample_acc": clean_meta_str(r.get("sample_accession", "")),
            "_date_is_fallback": (not bool(raw_date)),
        })

    rows = enrich_ena_via_biosample(rows)
    return rows


# ---------------------- BV-BRC MACRO FETCHER ----------------------

def parse_content_range_total(header_val):
    if not header_val:
        return 0
    m = re.search(r"/(\d+)", str(header_val))
    return int(m.group(1)) if m else 0


def bvbrc_genome_query(rql_query, timeout=(10, 45)):
    try:
        resp = requests.post(BVBRC_GENOME_URL, data=rql_query, headers=BVBRC_HEADERS, timeout=timeout)
        total = parse_content_range_total(resp.headers.get("Content-Range", ""))
        if resp.status_code == 200:
            data = resp.json()
            return 200, (data if isinstance(data, list) else []), total
        print(
            f"⚠️ Warning: BV-BRC HTTP {resp.status_code}: {resp.text[:180].strip()}",
            flush=True,
        )
        return resp.status_code, [], total
    except Exception as e:
        print(f"⚠️ Warning: BV-BRC genome query failed ({e})", flush=True)
        return 0, [], 0


def fetch_macro_bvbrc(organism_input, target_quota, min_len, max_len, strategy, random_seed):
    text = str(organism_input).strip()
    if any(tok in text for tok in ("eq(", "in(", "keyword(")):
        tax_filter = text
        alt_tax_filter = None
    else:
        emit_progress(6, f"Resolving NCBI TaxID for '{text}'...")
        taxid = resolve_taxid(text)
        if not taxid:
            print(f"❌ Could not resolve NCBI TaxID for BV-BRC query '{text}'.", file=sys.stderr, flush=True)
            return []
        tax_filter = f"eq(taxon_id,{taxid})"
        alt_tax_filter = f"eq(taxon_lineage_ids,{taxid})"

    emit_progress(12, "Probing BV-BRC genome archive size...")
    base_filter = f"{tax_filter}&ge(genome_length,{min_len})&le(genome_length,{max_len})"
    _, _, total_count = bvbrc_genome_query(f"{base_filter}&select(genome_id)&limit(1)", timeout=(8, 25))
    if total_count == 0 and alt_tax_filter:
        base_filter = f"{alt_tax_filter}&ge(genome_length,{min_len})&le(genome_length,{max_len})"
        _, _, total_count = bvbrc_genome_query(f"{base_filter}&select(genome_id)&limit(1)", timeout=(8, 25))

    if total_count == 0:
        return []

    items = []
    if strategy == "Newest First (Deterministic)":
        emit_progress(25, f"Fetching newest {target_quota:,} BV-BRC genome records ({total_count:,} available)...")
        _, items, _ = bvbrc_genome_query(
            f"{base_filter}&{BVBRC_FIELDS}&sort(-date_inserted,+genome_id)&limit({target_quota})"
        )
    elif strategy == "Oldest First (Deterministic)":
        emit_progress(25, f"Fetching oldest {target_quota:,} BV-BRC genome records ({total_count:,} available)...")
        status, items, _ = bvbrc_genome_query(
            f"{base_filter}&gt(date_inserted,1980-01-01T00%3A00%3A00.000Z)&{BVBRC_FIELDS}&sort(+date_inserted,+genome_id)&limit({target_quota})"
        )
        if status != 200 or not items:
            emit_progress(35, "Retrying BV-BRC oldest query without server-side date filter...")
            _, items, _ = bvbrc_genome_query(
                f"{base_filter}&{BVBRC_FIELDS}&sort(+date_inserted,+genome_id)&limit({max(target_quota * 2, 200)})"
            )
            valid_items = [
                it for it in items
                if not str(it.get("date_inserted", "")).startswith(("1969", "1970", "1900", "0000"))
            ]
            items = (valid_items or items)[:target_quota]
    elif strategy == "Longest Sequences First":
        emit_progress(25, f"Fetching longest {target_quota:,} BV-BRC genome records ({total_count:,} available)...")
        _, items, _ = bvbrc_genome_query(
            f"{base_filter}&{BVBRC_FIELDS}&sort(-genome_length,+genome_id)&limit({target_quota})"
        )
    else:
        rng = random.Random(random_seed)
        if total_count <= min(target_quota * 2, 5000):
            emit_progress(30, f"Downloading BV-BRC metadata pool ({total_count:,} records)...")
            _, items, _ = bvbrc_genome_query(f"{base_filter}&{BVBRC_FIELDS}&sort(+genome_id)&limit({min(total_count, 5000)})")
            items.sort(key=lambda x: str(x.get("genome_id", "")))
            rng.shuffle(items)
            items = items[:target_quota]
        else:
            n_slices = 10
            slice_size = max(50, int((target_quota * 1.3) // n_slices) + 1)
            safe_span = max(0, min(total_count - slice_size, 150000))
            sort_modes = ["+genome_id", "-genome_id", "+date_inserted", "-date_inserted"]
            buckets = []
            for idx in range(n_slices):
                sort_expr = sort_modes[idx % len(sort_modes)]
                off = rng.randint(0, safe_span) if safe_span > 0 else 0
                pct = 18 + int(54 * (idx / n_slices))
                emit_progress(pct, f"Fetching BV-BRC archive slice {idx + 1}/{n_slices} ({sort_expr}, offset {off:,})...")
                _, sl, _ = bvbrc_genome_query(f"{base_filter}&{BVBRC_FIELDS}&sort({sort_expr})&limit({slice_size},{off})")
                if sl:
                    sl.sort(key=lambda x: str(x.get("genome_id", "")))
                    rng.shuffle(sl)
                    buckets.append(sl)
            seen_gids = set()
            max_len_b = max((len(b) for b in buckets), default=0)
            for idx in range(max_len_b):
                for b in buckets:
                    if idx < len(b):
                        gid = str(b[idx].get("genome_id", ""))
                        if gid and gid not in seen_gids:
                            seen_gids.add(gid)
                            items.append(b[idx])
            items = items[:target_quota]

    emit_progress(75, f"Normalizing spatiotemporal fields for {len(items):,} BV-BRC records...")
    rows = []
    seen_accs = set()
    for it in items:
        gb_raw = clean_meta_str(it.get("genbank_accessions", ""))
        first_gb = re.split(r"[,\s;]+", gb_raw)[0] if gb_raw else ""
        raw_id = first_gb or str(it.get("genome_id", ""))
        acc = clean_accession(raw_id)
        if not acc or acc in seen_accs:
            continue
        seen_accs.add(acc)

        title = clean_meta_str(it.get("genome_name", "")) or acc
        update_date = (
            sanitize_date(it.get("date_modified", ""), is_db_timestamp=True)
            or sanitize_date(it.get("date_inserted", ""), is_db_timestamp=True)
        )
        raw_date_field = clean_meta_str(it.get("collection_date", ""))
        if raw_date_field:
            raw_col_date = sanitize_date(raw_date_field, is_db_timestamp=False)
        else:
            raw_col_date = sanitize_date(it.get("collection_year", ""), is_db_timestamp=False)

        fallback_date = sanitize_date(it.get("date_inserted", ""), is_db_timestamp=True) or update_date
        raw_country = clean_meta_str(it.get("isolation_country", "")) or clean_meta_str(it.get("geographic_location", ""))
        country, col_date, iso3 = parse_spatiotemporal_fields(title, raw_country, raw_col_date, fallback_date)

        rows.append({
            "Accession": acc,
            "Title": title,
            "Length_bp": int(it.get("genome_length") or 0),
            "UpdateDate": update_date,
            "Collection_Date": col_date,
            "Country": country,
            "ISO3": iso3,
            "Data_Source": "BV-BRC",
        })
    return rows


def main():
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else "config/macro_config.yaml"
    with open(cfg_path, "r") as f:
        config = yaml.safe_load(f) or {}

    os.makedirs("results/macro", exist_ok=True)
    email = config.get("ncbi_email") or "test@example.com"
    Entrez.email = email

    data_source = str(config.get("data_source", "NCBI")).upper()
    organism_input = str(config.get("organism_query", "Severe acute respiratory syndrome coronavirus 2")).strip()
    target_quota = int(config.get("max_records", 1000))
    min_len = int(config.get("min_length", 10000))
    max_len = int(config.get("max_length", 500000))
    strategy = config.get("selection_strategy", "Newest First (Deterministic)")
    random_seed = int(config.get("random_seed", 42))

    emit_progress(2, f"Starting Macro Surveillance ({data_source} | N={target_quota:,} | {strategy})...")

    if data_source == "ENA":
        rows = fetch_macro_ena(organism_input, target_quota, min_len, max_len, strategy, random_seed)
    elif data_source == "BV-BRC":
        rows = fetch_macro_bvbrc(organism_input, target_quota, min_len, max_len, strategy, random_seed)
    else:
        rows = fetch_macro_ncbi(organism_input, target_quota, min_len, max_len, strategy, random_seed, email)

    if rows and data_source in ("BV-BRC", "ENA"):
        rows = enrich_missing_via_ncbi(rows, email=email)

    cols = ["Accession", "Title", "Length_bp", "UpdateDate", "Collection_Date", "Country", "ISO3", "Data_Source"]
    out_df = pd.DataFrame(rows, columns=cols) if rows else pd.DataFrame(columns=cols)
    if out_df.empty:
        print("❌ No metadata records returned for the selected Macro query.", file=sys.stderr, flush=True)
        sys.exit(1)

    emit_progress(97, f"Writing {len(out_df):,} records to results/macro/macro_metadata.csv...")
    out_path = "results/macro/macro_metadata.csv"
    out_df.to_csv(out_path, index=False)
    emit_progress(100, f"Completed! Saved {len(out_df):,} metadata records.")


if __name__ == "__main__":
    main()