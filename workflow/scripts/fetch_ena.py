import os
import io
import re
import sys
import yaml
import random
import calendar
import datetime
import requests
import pandas as pd
from Bio import Entrez, SeqIO


def clean_accession(raw_id):
    """Extracts clean base accession from NCBI or ENA."""
    cleaned = str(raw_id).strip().strip("|")
    if "|" in cleaned:
        cleaned = cleaned.split("|")[-1]
    match = re.search(r"([A-Z]{1,6}_?[0-9]{4,})", cleaned)
    return match.group(1) if match else cleaned.split(".")[0]


ISO_COUNTRY_MAP = {
    "USA": "USA", "UNITED STATES": "USA", "US": "USA", "U.S.A.": "USA",
    "CHN": "China", "PRC": "China", "CHINA": "China",
    "CHL": "Chile", "CHILE": "Chile",
    "GBR": "United Kingdom", "UK": "United Kingdom", "ENGLAND": "United Kingdom",
    "SCOTLAND": "United Kingdom", "WALES": "United Kingdom", "NORTHERN IRELAND": "United Kingdom",
    "DEU": "Germany", "GERMANY": "Germany",
    "FRA": "France", "FRANCE": "France",
    "ITA": "Italy", "ITALY": "Italy",
    "ESP": "Spain", "SPAIN": "Spain",
    "CAN": "Canada", "CANADA": "Canada",
    "AUS": "Australia", "AUSTRALIA": "Australia",
    "JPN": "Japan", "JAPAN": "Japan",
    "KOR": "South Korea", "SOUTH KOREA": "South Korea",
    "BRA": "Brazil", "BRAZIL": "Brazil",
    "IND": "India", "INDIA": "India",
    "ZAF": "South Africa", "SOUTH AFRICA": "South Africa",
    "MEX": "Mexico", "MEXICO": "Mexico",
    "RUS": "Russia", "RUSSIA": "Russia",
    "ARG": "Argentina", "COL": "Colombia", "PER": "Peru",
    "NLD": "Netherlands", "CHE": "Switzerland", "SWE": "Sweden",
    "BEL": "Belgium", "AUT": "Austria", "POL": "Poland", "DNK": "Denmark",
    "NOR": "Norway", "FIN": "Finland", "PRT": "Portugal", "GRC": "Greece",
    "TUR": "Turkey", "ISR": "Israel", "EGY": "Egypt", "NGA": "Nigeria",
    "KEN": "Kenya", "THA": "Thailand", "VNM": "Vietnam", "SGP": "Singapore",
    "MYS": "Malaysia", "PHL": "Philippines", "IDN": "Indonesia", "NZL": "New Zealand",
    "PAK": "Pakistan", "BGD": "Bangladesh", "SAU": "Saudi Arabia", "ARE": "UAE",
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
    "KANAGAWA": "Japan", "YOKOHAMA": "Japan",
}

MISSING_TOKENS = {
    "", "nan", "none", "null", "unknown", "missing",
    "not collected", "not applicable", "not provided", "na", "n/a",
}

ENA_PORTAL_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_PORTAL_COUNT = "https://www.ebi.ac.uk/ena/portal/api/count"
ENA_FIELDS = "accession,description,base_count,first_public,last_updated,country,collection_date"


def clean_meta_str(val):
    """Returns empty string for NaN/null/placeholder metadata values."""
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    text = str(val).strip()
    return "" if text.lower() in MISSING_TOKENS else text


def normalize_country(raw_val):
    cleaned = clean_meta_str(raw_val)
    if not cleaned:
        return "Unknown"
    base = cleaned.split(":")[0].strip(" ,;()[]")
    if not base or base.lower() in MISSING_TOKENS:
        return "Unknown"
    upper = base.upper()
    if upper in ISO_COUNTRY_MAP:
        return ISO_COUNTRY_MAP[upper]
    if upper in REGION_TO_COUNTRY:
        return REGION_TO_COUNTRY[upper]
    return base if upper in ("USA", "UK", "UAE") else base.title()


def parse_spatiotemporal_fields(title, raw_country="", raw_date="", fallback_date=""):
    """Normalizes structured country/date fields with title-parsing and date fallbacks."""
    country = normalize_country(raw_country)

    if country == "Unknown" and "/" in str(title):
        slash_tokens = [tok.strip(" ,;()[]") for tok in str(title).split("/")[1:]]
        skip_terms = {
            "HUMAN", "HOMO SAPIENS", "SWINE", "BAT", "AVIAN", "CAT", "DOG",
            "MINK", "DEER", "RODENT", "MONKEY", "CAMEL", "MAMMAL", "ENV",
            "ENVIRONMENT", "WASTEWATER", "ISOLATE", "STRAIN", "COMPLETE",
            "GENOME", "SARS-COV-2", "HCOV-19",
        }
        for tok in slash_tokens:
            cand = tok.upper()
            if not cand or cand in skip_terms or cand.isdigit():
                continue
            if cand in ISO_COUNTRY_MAP or cand in REGION_TO_COUNTRY:
                country = normalize_country(cand)
                break

    if country == "Unknown":
        title_upper = str(title).upper()
        for code, mapped in ISO_COUNTRY_MAP.items():
            if len(code) > 3 and code in title_upper:
                country = mapped
                break
        if country == "Unknown":
            for reg, mapped in REGION_TO_COUNTRY.items():
                if reg in title_upper:
                    country = mapped
                    break

    collection_date = clean_meta_str(raw_date)
    if collection_date and "T" in collection_date:
        collection_date = collection_date.split("T")[0]

    if not collection_date:
        m_ymd = re.search(r"/(20[0-2][0-9])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]*/", str(title))
        if m_ymd:
            collection_date = f"{m_ymd.group(1)}-{m_ymd.group(2)}-{m_ymd.group(3)}"
        else:
            m_yymmdd = re.search(r"-(2[0-6])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]{2,}/20\1", str(title))
            if m_yymmdd:
                collection_date = f"20{m_yymmdd.group(1)}-{m_yymmdd.group(2)}-{m_yymmdd.group(3)}"

    if not collection_date:
        fb = clean_meta_str(fallback_date)
        collection_date = fb.split("T")[0] if "T" in fb else fb

    return country, collection_date or "Unknown"


def load_local_fasta_isolates(config, seen_accessions, min_len, max_n_percent):
    """Loads user-uploaded local FASTA sequences and applies QC filters."""
    local_path = config.get("local_fasta_path", "")
    local_records = []
    local_meta = []

    if not local_path or not os.path.exists(local_path):
        return local_records, local_meta

    print(f"📂 Reading uploaded local FASTA file: {local_path}")
    for rec in SeqIO.parse(local_path, "fasta"):
        raw_header = rec.description or rec.id
        sanitized = re.sub(r"[^A-Za-z0-9_.-]", "_", str(rec.id).strip())
        clean_id = clean_accession(sanitized)
        seq_str = str(rec.seq).upper()
        if not seq_str:
            continue

        n_percent = (seq_str.count("N") / len(seq_str)) * 100.0
        if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
            seen_accessions.add(clean_id)
            rec.id = clean_id
            rec.name = clean_id
            rec.description = ""
            local_records.append(rec)
            local_meta.append({
                "Accession": clean_id,
                "Title": f"[Local FASTA] {raw_header}",
                "Length_bp": len(seq_str),
                "UpdateDate": "Local",
                "Collection_Date": "Local",
                "Country": "Local Upload",
            })

    print(f"✅ Accepted {len(local_records)} local FASTA sequence(s) after QC.")
    return local_records, local_meta


def get_nucleotide_ids_for_assembly(assembly_acc):
    """Finds linked nucleotide UIDs for a GCA/GCF assembly without duplicating RefSeq+GenBank or grabbing WGS masters."""
    search_handle = Entrez.esearch(db="assembly", term=f"{assembly_acc}[Assembly Accession]")
    search_record = Entrez.read(search_handle)
    if not search_record["IdList"]:
        search_handle = Entrez.esearch(db="assembly", term=assembly_acc)
        search_record = Entrez.read(search_handle)
        if not search_record["IdList"]:
            return None

    asm_id = search_record["IdList"][0]
    link_handle = Entrez.elink(dbfrom="assembly", db="nucleotide", id=asm_id)
    link_record = Entrez.read(link_handle)

    try:
        link_sets = {
            ls["LinkName"]: [link["Id"] for link in ls["Link"]]
            for ls in link_record[0].get("LinkSetDb", [])
        }
        if assembly_acc.upper().startswith("GCF_"):
            preferred_order = ["assembly_nuccore_refseq", "assembly_nuccore_insdc"]
        else:
            preferred_order = ["assembly_nuccore_insdc", "assembly_nuccore_refseq"]

        for name in preferred_order:
            if name in link_sets and link_sets[name]:
                return link_sets[name]

        for name, ids in link_sets.items():
            if "wgsmaster" not in name and ids:
                return ids
    except (IndexError, KeyError):
        pass
    return None


def ena_count(query):
    """Queries ENA Portal /count endpoint and returns integer hit count."""
    try:
        resp = requests.get(
            ENA_PORTAL_COUNT,
            params={"result": "sequence", "query": query},
            timeout=25,
        )
        if resp.status_code == 200:
            digits = re.search(r"(\d+)", resp.text)
            return int(digits.group(1)) if digits else 0
    except Exception:
        pass
    return 0


def ena_search_df(query, limit=250, offset=0):
    """Queries ENA Portal /search endpoint and returns a cleaned DataFrame."""
    params = {
        "result": "sequence",
        "query": query,
        "limit": int(limit),
        "format": "tsv",
        "fields": ENA_FIELDS,
    }

    try:
        resp = requests.get(ENA_PORTAL_SEARCH, params=params, timeout=30)
        if resp.status_code != 200:
            return pd.DataFrame()
        txt = resp.text.strip()
        if not txt:
            return pd.DataFrame()
        df = pd.read_csv(io.StringIO(txt), sep="\t")
        if "accession" in df.columns and not df.empty:
            df["accession"] = df["accession"].astype(str).apply(clean_accession)
            return df.drop_duplicates(subset=["accession"])
    except Exception as e:
        print(f"⚠️ Warning: ENA sub-query skipped ({e})")
    return pd.DataFrame()


def find_ena_min_year(base_query, current_year, floor_year=1980):
    """Binary searches ENA first_public year to locate earliest public records (<=6 count calls)."""
    total = ena_count(base_query)
    if total == 0:
        return current_year, 0

    lo, hi = floor_year, current_year
    while lo < hi:
        mid = (lo + hi) // 2
        cnt = ena_count(f"({base_query}) AND first_public<={mid}-12-31")
        if cnt > 0:
            hi = mid
        else:
            lo = mid + 1
    return lo, total


def interleave_dataframes(df_list):
    """Interleaves candidate DataFrames row-by-row while deduplicating on accession."""
    if not df_list:
        return pd.DataFrame()
    rows = []
    seen = set()
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


def build_ena_candidate_df(ena_query, min_len, max_len, max_ref_len, target_quota, selection_strategy, random_seed, seen_accessions):
    """Builds a deterministic ENA candidate DataFrame for all 7 selection strategies."""
    base_query = f"({ena_query}) AND base_count>={min_len} AND base_count<={max_len}"
    today = datetime.date.today()
    current_year = today.year
    current_month = today.month

    is_two_stage = selection_strategy in (
        "Newest + Best Alignment (Two-Stage)",
        "Fewest SNPs (Two-Stage)",
    )
    db_fetch_limit = max(target_quota * 2, 60) if selection_strategy == "Fewest SNPs (Two-Stage)" else (
        target_quota * 2 if is_two_stage else target_quota
    )
    desired_pool = max(db_fetch_limit * 6, 300)

    if selection_strategy in ("Newest First (Deterministic)", "Newest + Best Alignment (Two-Stage)"):
        print(f"🌍 Querying ENA Portal API with rolling recency window ({selection_strategy})...")
        frames = []
        collected = 0
        for days_back in (7, 30, 90, 365, None):
            if days_back is not None:
                cutoff = (today - datetime.timedelta(days=days_back)).isoformat()
                q = f"({base_query}) AND last_updated>={cutoff}"
            else:
                q = base_query
            df_part = ena_search_df(q, limit=desired_pool)
            if not df_part.empty:
                df_part = df_part[~df_part["accession"].isin(seen_accessions)]
                frames.append(df_part)
                collected = len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]))
            if collected >= desired_pool:
                break

        meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]) if frames else pd.DataFrame()
        if not meta_df.empty and "last_updated" in meta_df.columns:
            meta_df = meta_df.sort_values(by=["last_updated", "accession"], ascending=[False, True]).reset_index(drop=True)
        return meta_df, db_fetch_limit

    if selection_strategy == "Oldest First (Deterministic)":
        min_year, total = find_ena_min_year(base_query, current_year)
        print(f"🕰️ Scanning ENA chronologically from first_public year {min_year} (total hits: {total:,})...")
        frames = []
        collected = 0
        for yr in range(min_year, current_year + 1):
            yr_q = f"({base_query}) AND first_public>={yr}-01-01 AND first_public<={yr}-12-31"
            yr_cnt = ena_count(yr_q)
            if yr_cnt == 0:
                continue
            if yr_cnt <= desired_pool:
                df_part = ena_search_df(yr_q, limit=desired_pool)
                if not df_part.empty:
                    df_part = df_part[~df_part["accession"].isin(seen_accessions)]
                    frames.append(df_part)
            else:
                max_m = current_month if yr == current_year else 12
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    m_q = f"({base_query}) AND first_public>={yr}-{m:02d}-01 AND first_public<={yr}-{m:02d}-{last_day:02d}"
                    df_part = ena_search_df(m_q, limit=desired_pool)
                    if not df_part.empty:
                        df_part = df_part[~df_part["accession"].isin(seen_accessions)]
                        frames.append(df_part)
                    if frames and len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"])) >= desired_pool:
                        break
            if frames:
                collected = len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]))
            if collected >= desired_pool:
                break

        meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]) if frames else pd.DataFrame()
        if meta_df.empty:
            meta_df = ena_search_df(base_query, limit=desired_pool)
        if not meta_df.empty:
            sort_col = "first_public" if "first_public" in meta_df.columns else "last_updated"
            meta_df = meta_df.sort_values(by=[sort_col, "accession"], ascending=[True, True]).reset_index(drop=True)
        return meta_df, db_fetch_limit

    if selection_strategy == "Longest Sequences First":
        print("📏 Locating longest ENA length threshold via server-side base_count binary search...")
        target_band = max(target_quota * 4, 150)
        lo, hi = min_len, max_len
        best_low = min_len
        for _ in range(8):
            if lo > hi:
                break
            mid = (lo + hi) // 2
            cnt = ena_count(f"({ena_query}) AND base_count>={mid} AND base_count<={max_len}")
            if cnt >= target_band:
                best_low = mid
                lo = mid + 1
            else:
                hi = mid - 1

        meta_df = ena_search_df(
            f"({ena_query}) AND base_count>={best_low} AND base_count<={max_len}",
            limit=desired_pool,
        )
        if meta_df.empty or len(meta_df) < target_quota:
            fallback_df = ena_search_df(base_query, limit=desired_pool)
            meta_df = pd.concat([meta_df, fallback_df], ignore_index=True).drop_duplicates(subset=["accession"])
        if not meta_df.empty:
            meta_df = meta_df[~meta_df["accession"].isin(seen_accessions)]
            meta_df["_bc"] = pd.to_numeric(meta_df["base_count"], errors="coerce").fillna(0)
            meta_df = meta_df.sort_values(by=["_bc", "accession"], ascending=[False, True]).drop(columns=["_bc"]).reset_index(drop=True)
        return meta_df, db_fetch_limit

    if selection_strategy == "Closest to Reference Length":
        print(f"📏 Querying progressive ENA base_count bands centered around {max_ref_len:,} bp...")
        frames = []
        for delta in (50, 250, 1000, None):
            if delta is None:
                lo_b, hi_b = min_len, max_len
            else:
                lo_b = max(min_len, max_ref_len - delta)
                hi_b = min(max_len, max_ref_len + delta)
            band_q = f"({ena_query}) AND base_count>={lo_b} AND base_count<={hi_b}"
            df_part = ena_search_df(band_q, limit=desired_pool)
            if not df_part.empty:
                df_part = df_part[~df_part["accession"].isin(seen_accessions)]
                frames.append(df_part)
            if frames and len(pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"])) >= desired_pool:
                break

        meta_df = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["accession"]) if frames else pd.DataFrame()
        if not meta_df.empty:
            meta_df["_len_diff"] = (
                pd.to_numeric(meta_df["base_count"], errors="coerce").fillna(0) - max_ref_len
            ).abs()
            meta_df = meta_df.sort_values(by=["_len_diff", "accession"], ascending=[True, True]).drop(columns=["_len_diff"]).reset_index(drop=True)
        return meta_df, db_fetch_limit

    if selection_strategy == "Fewest SNPs (Two-Stage)":
        min_year, total = find_ena_min_year(base_query, current_year)
        print(f"🎯 Stage 1 ({selection_strategy}): Sampling multi-era ENA pool from {min_year}-{current_year} (total: {total:,})...")
        all_years = list(range(min_year, current_year + 1))
        if len(all_years) <= 6:
            sample_years = all_years
        else:
            step = (len(all_years) - 1) / 5.0
            sample_years = sorted({all_years[int(round(i * step))] for i in range(6)})

        per_bucket = max(25, (db_fetch_limit * 3) // max(1, len(sample_years)))
        buckets = []
        for idx, yr in enumerate(sample_years):
            if idx == 0:
                max_m = current_month if yr == current_year else 12
                earliest_frames = []
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    m_q = f"({base_query}) AND first_public>={yr}-{m:02d}-01 AND first_public<={yr}-{m:02d}-{last_day:02d}"
                    m_cnt = ena_count(m_q)
                    if m_cnt == 0:
                        continue
                    df_m = ena_search_df(m_q, limit=per_bucket, offset=0)
                    if not df_m.empty:
                        earliest_frames.append(df_m[~df_m["accession"].isin(seen_accessions)])
                    if m_cnt > per_bucket * 2:
                        mid_off = min(m_cnt // 2, 2000)
                        df_mid = ena_search_df(m_q, limit=per_bucket, offset=mid_off)
                        if not df_mid.empty:
                            earliest_frames.append(df_mid[~df_mid["accession"].isin(seen_accessions)])
                    if earliest_frames and sum(len(x) for x in earliest_frames) >= per_bucket * 2:
                        break
                if earliest_frames:
                    buckets.append(pd.concat(earliest_frames, ignore_index=True).drop_duplicates(subset=["accession"]))
            else:
                yr_q = f"({base_query}) AND first_public>={yr}-01-01 AND first_public<={yr}-12-31"
                yr_cnt = ena_count(yr_q)
                if yr_cnt == 0:
                    continue
                off = min(yr_cnt // 2, 3000) if yr_cnt > per_bucket * 2 else 0
                df_y = ena_search_df(yr_q, limit=per_bucket, offset=off)
                if not df_y.empty:
                    buckets.append(df_y[~df_y["accession"].isin(seen_accessions)])

        meta_df = interleave_dataframes(buckets)
        if meta_df.empty:
            meta_df = ena_search_df(base_query, limit=desired_pool)
        return meta_df, db_fetch_limit

    # Random Sample (Reproducible Seed): stratified across first_public (year, month) windows + seeded offsets
    rng = random.Random(random_seed)
    min_year, total = find_ena_min_year(base_query, current_year)
    print(f"🎲 Sampling ENA across {min_year}-{current_year} with reproducible seed={random_seed} (total: {total:,})...")

    if total <= 1000:
        meta_df = ena_search_df(base_query, limit=1000)
        if not meta_df.empty:
            meta_df = meta_df[~meta_df["accession"].isin(seen_accessions)]
            meta_df = meta_df.sort_values(by=["accession"], ascending=[True])
            meta_df = meta_df.sample(frac=1.0, random_state=random_seed).reset_index(drop=True)
        return meta_df, db_fetch_limit

    ym_windows = [
        (yr, m)
        for yr in range(min_year, current_year + 1)
        for m in range(1, (current_month if yr == current_year else 12) + 1)
    ]
    rng.shuffle(ym_windows)

    target_windows = min(12, max(6, target_quota // 3))
    per_win = max(20, (target_quota * 4) // target_windows)
    buckets = []

    for yr, m in ym_windows:
        if len(buckets) >= target_windows:
            break
        last_day = calendar.monthrange(yr, m)[1]
        m_q = f"({base_query}) AND first_public>={yr}-{m:02d}-01 AND first_public<={yr}-{m:02d}-{last_day:02d}"
        m_cnt = ena_count(m_q)
        if m_cnt == 0:
            continue
        max_off = max(0, min(m_cnt - per_win, 10000))
        off = rng.randint(0, max_off) if max_off > 0 else 0
        df_m = ena_search_df(m_q, limit=per_win, offset=off)
        if not df_m.empty:
            df_m = df_m[~df_m["accession"].isin(seen_accessions)]
            if not df_m.empty:
                df_m = df_m.sort_values(by=["accession"], ascending=[True])
                df_m = df_m.sample(frac=1.0, random_state=rng.randint(0, 10**6)).reset_index(drop=True)
                buckets.append(df_m)

    meta_df = interleave_dataframes(buckets)
    if meta_df.empty:
        meta_df = ena_search_df(base_query, limit=desired_pool)
        if not meta_df.empty:
            meta_df = meta_df[~meta_df["accession"].isin(seen_accessions)]
            meta_df = meta_df.sort_values(by=["accession"], ascending=[True])
            meta_df = meta_df.sample(frac=1.0, random_state=random_seed).reset_index(drop=True)
    return meta_df, db_fetch_limit


def main():
    with open("config/config.yaml", "r") as f:
        config = yaml.safe_load(f)

    os.makedirs("results/raw", exist_ok=True)
    Entrez.email = config.get("ncbi_email") or "test@example.com"

    ref_acc = str(config["reference_accession"]).strip()
    base_ref = clean_accession(ref_acc)
    print(f"🧬 Fetching Reference Genome & GenBank Annotations: {ref_acc}")

    fetch_id = ref_acc
    if ref_acc.upper().startswith(("GCF_", "GCA_")):
        print(f"🔗 Assembly detected. Resolving {ref_acc} to nucleotide sequences...")
        nuc_uids = get_nucleotide_ids_for_assembly(ref_acc)
        if not nuc_uids:
            print(f"❌ Failed to find nucleotide sequences for assembly {ref_acc}", file=sys.stderr)
            sys.exit(1)
        fetch_id = ",".join(nuc_uids)

    try:
        with Entrez.efetch(db="nucleotide", id=fetch_id, rettype="gbwithparts", retmode="text") as handle:
            gb_text = handle.read()

        with open("results/raw/reference.gb", "w") as gb_out:
            gb_out.write(gb_text)

        ref_records = list(SeqIO.parse(io.StringIO(gb_text), "genbank"))
        if not ref_records:
            raise ValueError("No GenBank records found in reference handle.")

        original_ref_accs = []
        for rec in ref_records:
            cid = clean_accession(rec.id)
            original_ref_accs.append(cid)
            rec.id = cid
            rec.name = cid
            rec.description = ""

        SeqIO.write(ref_records, "results/raw/reference.fasta", "fasta")
        max_ref_len = max(len(r.seq) for r in ref_records)
        print(f"✅ Loaded {len(ref_records)} reference segment(s)/contig(s) (largest: {max_ref_len:,} bp).")

    except Exception as e:
        print(f"❌ Failed to fetch NCBI Reference: {e}", file=sys.stderr)
        sys.exit(1)

    # 2. Build ENA Query & Fetch Candidate Metadata Pool
    search_mode = config.get("search_mode", "By Organism & Keyword")
    selection_strategy = config.get("selection_strategy", "Newest First (Deterministic)")
    random_seed = int(config.get("random_seed", 42))
    ena_query = config.get("ena_query", "tax_tree(11250)")
    target_quota = int(config.get("max_sequences", 20))
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))

    seen_accessions = set(original_ref_accs + [base_ref])
    local_records, local_meta = load_local_fasta_isolates(
        config, seen_accessions, min_len, max_n_percent
    )

    db_records = []
    db_meta_rows = []
    meta_cols = ["Accession", "Title", "Length_bp", "UpdateDate", "Collection_Date", "Country"]

    if search_mode != "Local FASTA Only":
        if search_mode == "By Accession List":
            raw_accs = [x.strip() for x in re.split(r"[,\s]+", str(config.get("accession_list", ""))) if x.strip()]
            db_fetch_limit = max(target_quota, len(raw_accs))
            print(f"🌍 Querying ENA Portal API (By Accession List): {ena_query}")
            meta_df = ena_search_df(ena_query, limit=max(db_fetch_limit * 2, 200))
            if not meta_df.empty:
                meta_df = meta_df[~meta_df["accession"].isin(seen_accessions)]
        else:
            meta_df, db_fetch_limit = build_ena_candidate_df(
                ena_query=ena_query,
                min_len=min_len,
                max_len=max_len,
                max_ref_len=max_ref_len,
                target_quota=target_quota,
                selection_strategy=selection_strategy,
                random_seed=random_seed,
                seen_accessions=seen_accessions,
            )

        if meta_df.empty and not local_records:
            print(f"❌ ENA returned NO sequences for query: {ena_query}", file=sys.stderr)
            sys.exit(1)

        if not meta_df.empty:
            meta_lookup = {row["accession"]: row for _, row in meta_df.iterrows()}
            candidate_ids = meta_df["accession"].tolist()

            print(f"🧬 Actively hunting for up to {db_fetch_limit} pristine ENA sequence(s) from {len(candidate_ids)} candidate(s)...")
            chunk_size = 30

            for i in range(0, len(candidate_ids), chunk_size):
                chunk = candidate_ids[i : i + chunk_size]
                browser_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{','.join(chunk)}"
                try:
                    fasta_response = requests.get(browser_url, timeout=30)
                except Exception:
                    continue

                if fasta_response.status_code != 200:
                    continue

                records = list(SeqIO.parse(io.StringIO(fasta_response.text), "fasta"))
                rec_by_id = {clean_accession(r.id): r for r in records}
                for cid in chunk:
                    record = rec_by_id.get(cid)
                    if record is None:
                        continue

                    seq_str = str(record.seq).upper()
                    if len(seq_str) == 0:
                        continue

                    n_count = seq_str.count("N")
                    n_percent = (n_count / len(seq_str)) * 100.0

                    if cid not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                        seen_accessions.add(cid)
                        record.id = cid
                        record.name = cid
                        record.description = ""
                        db_records.append(record)

                        row = meta_lookup.get(cid, {})
                        title = clean_meta_str(row.get("description", "")) or cid
                        upd = clean_meta_str(row.get("last_updated", "")) or clean_meta_str(row.get("first_public", ""))
                        country, col_date = parse_spatiotemporal_fields(
                            title=title,
                            raw_country=row.get("country", ""),
                            raw_date=row.get("collection_date", ""),
                            fallback_date=clean_meta_str(row.get("first_public", "")) or upd,
                        )
                        db_meta_rows.append({
                            "Accession": cid,
                            "Title": title,
                            "Length_bp": int(row.get("base_count") or len(seq_str)),
                            "UpdateDate": upd,
                            "Collection_Date": col_date,
                            "Country": country,
                        })

                    if len(db_records) >= db_fetch_limit:
                        break
                if len(db_records) >= db_fetch_limit:
                    break

    valid_records = local_records + db_records
    if not valid_records:
        print("❌ No ENA or local sequences passed the QC filters.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")

    all_meta = local_meta + db_meta_rows
    final_meta_df = pd.DataFrame(all_meta, columns=meta_cols) if all_meta else pd.DataFrame(columns=meta_cols)
    final_meta_df.to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine sequence(s) and metadata.")


if __name__ == "__main__":
    main()