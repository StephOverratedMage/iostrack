import os
import io
import re
import sys
import json
import yaml
import random
import calendar
import datetime
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


def normalize_country(raw_val):
    if not raw_val:
        return "Unknown"
    base = str(raw_val).split(":")[0].strip(" ,;()[]")
    if not base:
        return "Unknown"
    upper = base.upper()
    if upper in ISO_COUNTRY_MAP:
        return ISO_COUNTRY_MAP[upper]
    if upper in REGION_TO_COUNTRY:
        return REGION_TO_COUNTRY[upper]
    return base if upper in ("USA", "UK", "UAE") else base.title()


def parse_spatiotemporal_from_text(title, sub_type_str="", sub_name_str="", fallback_date=""):
    """Extracts (country, collection_date) from SubType/SubName or FASTA/GenBank title."""
    sub_types = str(sub_type_str).split("|")
    sub_names = str(sub_name_str).split("|")
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

    collection_date = (
        sub_dict.get("collection_date")
        or sub_dict.get("collection-date")
        or ""
    )
    if not collection_date:
        m_ymd = re.search(r"/(20[0-2][0-9])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]*/", str(title))
        if m_ymd:
            collection_date = f"{m_ymd.group(1)}-{m_ymd.group(2)}-{m_ymd.group(3)}"
        else:
            m_yymmdd = re.search(r"-(2[0-6])(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])[0-9]{2,}/20\1", str(title))
            if m_yymmdd:
                collection_date = f"20{m_yymmdd.group(1)}-{m_yymmdd.group(2)}-{m_yymmdd.group(3)}"

    if not collection_date:
        collection_date = str(fallback_date or "")

    return country, collection_date


def load_local_fasta_isolates(config, seen_accessions, min_len, max_n_percent):
    """Loads user-uploaded local FASTA sequences and applies QC filters."""
    local_path = config.get("local_fasta_path", "")
    local_records = []
    local_meta = {}

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
            local_meta[clean_id] = {
                "Accession": clean_id,
                "Title": f"[Local FASTA] {raw_header}",
                "Length_bp": len(seq_str),
                "UpdateDate": "Local",
                "Collection_Date": "Local",
                "Country": "Local Upload",
            }

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


def ncbi_esearch(term, retmax=500, retstart=0, sort="Date Released", mindate=None, maxdate=None, datetype="pdat", count_only=False):
    """Wrapper around Entrez.esearch returning (total_count, id_list)."""
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
        kwargs["datetype"] = datetype
    try:
        with Entrez.esearch(**kwargs) as handle:
            rec = Entrez.read(handle)
        return int(rec.get("Count", 0)), list(rec.get("IdList", []))
    except Exception as e:
        print(f"⚠️ Warning: NCBI esearch sub-query skipped ({e})")
        return 0, []


def find_ncbi_min_year(base_query, current_year, floor_year=1980):
    """Binary searches [floor_year, current_year] for the earliest release year with hits (<=6 lightweight count calls)."""
    total, _ = ncbi_esearch(base_query, count_only=True)
    if total == 0:
        return current_year, 0

    lo, hi = floor_year, current_year
    while lo < hi:
        mid = (lo + hi) // 2
        cnt, _ = ncbi_esearch(
            base_query,
            mindate=f"{floor_year}/01/01",
            maxdate=f"{mid}/12/31",
            count_only=True,
        )
        if cnt > 0:
            hi = mid
        else:
            lo = mid + 1
    return lo, total


def interleave_buckets(buckets):
    """Deterministically interleaves lists of UIDs while preserving deduplication."""
    out = []
    seen = set()
    max_len = max((len(b) for b in buckets), default=0)
    for idx in range(max_len):
        for b in buckets:
            if idx < len(b) and b[idx] not in seen:
                seen.add(b[idx])
                out.append(b[idx])
    return out


def build_ncbi_candidate_pool(ncbi_query, exclude_str, min_len, max_len, max_ref_len, target_quota, selection_strategy, random_seed):
    """Constructs a deterministic candidate UID list and Stage 1 fetch limit for all 7 selection strategies."""
    base_query = f"({ncbi_query}) AND {min_len}:{max_len}[SLEN] {exclude_str}".strip()
    today = datetime.date.today()
    current_year = today.year
    current_month = today.month

    if selection_strategy == "Newest First (Deterministic)":
        print(f"🔍 Executing NCBI Query ({selection_strategy}): {base_query}")
        _, id_list = ncbi_esearch(base_query, retmax=max(target_quota * 5, 500), sort="Date Released")
        return id_list, target_quota

    if selection_strategy == "Newest + Best Alignment (Two-Stage)":
        db_fetch_limit = target_quota * 2
        print(f"🎯 Stage 1 ({selection_strategy}): Fetching newest 2x pool ({db_fetch_limit} target)...")
        _, id_list = ncbi_esearch(base_query, retmax=max(db_fetch_limit * 4, 500), sort="Date Released")
        return id_list, db_fetch_limit

    if selection_strategy == "Longest Sequences First":
        print(f"📏 Executing NCBI Server-Side Sort ({selection_strategy}): {base_query}")
        _, id_list = ncbi_esearch(base_query, retmax=max(target_quota * 5, 500), sort="Sequence Length")
        return id_list, target_quota

    if selection_strategy == "Closest to Reference Length":
        db_fetch_limit = max(target_quota * 3, 90)
        desired_uids = max(db_fetch_limit * 2, 200)
        print(f"📏 Querying progressive NCBI [SLEN] bands centered around {max_ref_len:,} bp...")
        id_list = []
        seen = set()
        for delta in (50, 250, 1000, None):
            if delta is None:
                lo, hi = min_len, max_len
            else:
                lo = max(min_len, max_ref_len - delta)
                hi = min(max_len, max_ref_len + delta)
            band_query = f"({ncbi_query}) AND {lo}:{hi}[SLEN] {exclude_str}".strip()
            _, band_ids = ncbi_esearch(band_query, retmax=desired_uids, sort="Date Released")
            for uid in band_ids:
                if uid not in seen:
                    seen.add(uid)
                    id_list.append(uid)
            if len(id_list) >= desired_uids:
                break
        return id_list, db_fetch_limit

    if selection_strategy == "Oldest First (Deterministic)":
        min_year, total = find_ncbi_min_year(base_query, current_year)
        print(f"🕰️ Scanning NCBI archive chronologically from earliest year {min_year} (total hits: {total:,})...")
        desired_uids = max(target_quota * 5, 250)
        id_list = []
        seen = set()

        for yr in range(min_year, current_year + 1):
            yr_cnt, _ = ncbi_esearch(base_query, mindate=f"{yr}/01/01", maxdate=f"{yr}/12/31", count_only=True)
            if yr_cnt == 0:
                continue
            if yr_cnt <= 500:
                _, yr_ids = ncbi_esearch(
                    base_query, retmax=500, sort="Date Released",
                    mindate=f"{yr}/01/01", maxdate=f"{yr}/12/31"
                )
                for uid in reversed(yr_ids):
                    if uid not in seen:
                        seen.add(uid)
                        id_list.append(uid)
            else:
                max_m = current_month if yr == current_year else 12
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    m_start = f"{yr}/{m:02d}/01"
                    m_end = f"{yr}/{m:02d}/{last_day:02d}"
                    m_cnt, _ = ncbi_esearch(base_query, mindate=m_start, maxdate=m_end, count_only=True)
                    if m_cnt == 0:
                        continue
                    retstart = min(max(0, m_cnt - 250), 9500)
                    _, m_ids = ncbi_esearch(
                        base_query, retmax=250, retstart=retstart, sort="Date Released",
                        mindate=m_start, maxdate=m_end
                    )
                    for uid in reversed(m_ids):
                        if uid not in seen:
                            seen.add(uid)
                            id_list.append(uid)
                    if len(id_list) >= desired_uids:
                        break
            if len(id_list) >= desired_uids:
                break

        if not id_list:
            _, id_list = ncbi_esearch(base_query, retmax=desired_uids, sort="Date Released")
            id_list.reverse()
        return id_list, target_quota

    if selection_strategy == "Fewest SNPs (Two-Stage)":
        db_fetch_limit = max(target_quota * 2, 60)
        min_year, total = find_ncbi_min_year(base_query, current_year)
        print(f"🎯 Stage 1 ({selection_strategy}): Sampling multi-era NCBI pool from {min_year}-{current_year} (total: {total:,})...")
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
                earliest_ids = []
                for m in range(1, max_m + 1):
                    last_day = calendar.monthrange(yr, m)[1]
                    m_start, m_end = f"{yr}/{m:02d}/01", f"{yr}/{m:02d}/{last_day:02d}"
                    m_cnt, _ = ncbi_esearch(base_query, mindate=m_start, maxdate=m_end, count_only=True)
                    if m_cnt == 0:
                        continue
                    retstart = min(max(0, m_cnt - per_bucket), 9500)
                    _, m_ids = ncbi_esearch(
                        base_query, retmax=per_bucket, retstart=retstart,
                        sort="Date Released", mindate=m_start, maxdate=m_end
                    )
                    earliest_ids.extend(reversed(m_ids))
                    if len(earliest_ids) >= per_bucket * 2:
                        break
                if earliest_ids:
                    buckets.append(earliest_ids)
            else:
                _, yr_ids = ncbi_esearch(
                    base_query, retmax=per_bucket, sort="Date Released",
                    mindate=f"{yr}/01/01", maxdate=f"{yr}/12/31"
                )
                if yr_ids:
                    buckets.append(yr_ids)

        id_list = interleave_buckets(buckets)
        if not id_list:
            _, id_list = ncbi_esearch(base_query, retmax=max(db_fetch_limit * 4, 300), sort="Date Released")
        return id_list, db_fetch_limit

    # Random Sample (Reproducible Seed): stratified across (year, month) windows + seeded retstart
    rng = random.Random(random_seed)
    min_year, total = find_ncbi_min_year(base_query, current_year)
    print(f"🎲 Sampling NCBI across {min_year}-{current_year} with reproducible seed={random_seed} (total: {total:,})...")

    if total <= 5000:
        _, id_list = ncbi_esearch(base_query, retmax=5000, sort="Date Released")
        id_list.sort()
        rng.shuffle(id_list)
        return id_list, target_quota

    ym_windows = [
        (yr, m)
        for yr in range(min_year, current_year + 1)
        for m in range(1, (current_month if yr == current_year else 12) + 1)
    ]
    rng.shuffle(ym_windows)

    target_windows = min(12, max(6, target_quota // 3))
    per_win = max(15, (target_quota * 4) // target_windows)
    buckets = []

    for yr, m in ym_windows:
        if len(buckets) >= target_windows:
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

    id_list = interleave_buckets(buckets)
    if not id_list:
        _, id_list = ncbi_esearch(base_query, retmax=5000, sort="Date Released")
        id_list.sort()
        rng.shuffle(id_list)
    return id_list, target_quota


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

    # 2. Get a pool of wild isolates (and/or local FASTA uploads)
    search_mode = config.get("search_mode", "By Organism & Keyword")
    selection_strategy = config.get("selection_strategy", "Newest First (Deterministic)")
    random_seed = int(config.get("random_seed", 42))
    ncbi_query = config.get("search_term", '"Respiratory syncytial virus"[Organism]')
    min_len = int(config.get("min_length", 10000))
    max_len = max(50000, int(max_ref_len * 1.5))
    max_n_percent = float(config.get("max_ambiguous_n_percent", 1.0))
    target_quota = int(config.get("max_sequences", 20))

    seen_accessions = set(original_ref_accs + [base_ref])
    local_records, local_meta_by_acc = load_local_fasta_isolates(
        config, seen_accessions, min_len, max_n_percent
    )
    db_records = []

    if search_mode != "Local FASTA Only":
        if search_mode == "By Accession List":
            raw_accs = [x.strip() for x in re.split(r"[,\s]+", str(config.get("accession_list", ""))) if x.strip()]
            id_list = []
            for acc in raw_accs:
                if acc.upper().startswith(("GCF_", "GCA_")):
                    asm_uids = get_nucleotide_ids_for_assembly(acc)
                    if asm_uids:
                        id_list.extend(asm_uids)
                else:
                    _, acc_ids = ncbi_esearch(f"{acc}[Accession]", retmax=5, sort=None)
                    id_list.extend(acc_ids)
            db_fetch_limit = max(target_quota, len(id_list))
            print(f"🔍 Resolved {len(id_list)} nucleotide record(s) from Accession List.")
        else:
            exclude_accs = original_ref_accs[:10]
            exclude_str = " ".join([f"NOT {acc}[Accession]" for acc in exclude_accs])
            id_list, db_fetch_limit = build_ncbi_candidate_pool(
                ncbi_query=ncbi_query,
                exclude_str=exclude_str,
                min_len=min_len,
                max_len=max_len,
                max_ref_len=max_ref_len,
                target_quota=target_quota,
                selection_strategy=selection_strategy,
                random_seed=random_seed,
            )

        print(f"📦 Prepared {len(id_list)} candidate UID(s) in NCBI.")
        if not id_list and not local_records:
            print("❌ NCBI returned no sequences for query.", file=sys.stderr)
            sys.exit(1)

        chunk_size = 30
        print("🧬 Actively hunting for pristine sequences...")

        for i in range(0, len(id_list), chunk_size):
            chunk_uids = id_list[i : i + chunk_size]
            try:
                with Entrez.efetch(db="nucleotide", id=",".join(chunk_uids), rettype="fasta", retmode="text") as fetch_handle:
                    records = list(SeqIO.parse(fetch_handle, "fasta"))
            except Exception:
                continue

            for rec in records:
                clean_id = clean_accession(rec.id)
                seq_str = str(rec.seq).upper()

                if len(seq_str) == 0:
                    continue

                n_count = seq_str.count("N")
                n_percent = (n_count / len(seq_str)) * 100.0

                if clean_id not in seen_accessions and len(seq_str) >= min_len and n_percent <= max_n_percent:
                    seen_accessions.add(clean_id)
                    rec.id = clean_id
                    rec.name = clean_id
                    rec.annotations["orig_title"] = rec.description
                    rec.annotations["n_percent"] = n_percent
                    rec.description = ""
                    db_records.append(rec)

                if len(db_records) >= db_fetch_limit:
                    break

            if len(db_records) >= db_fetch_limit:
                break

        if search_mode == "By Organism & Keyword" and selection_strategy == "Closest to Reference Length":
            db_records.sort(
                key=lambda r: (
                    abs(len(r.seq) - max_ref_len),
                    round(r.annotations.get("n_percent", 0.0), 4),
                    r.id,
                )
            )
            db_records = db_records[:target_quota]
            print(f"📏 Selected top {len(db_records)} isolate(s) closest to reference length ({max_ref_len:,} bp).")

    valid_records = local_records + db_records
    if not valid_records:
        print("❌ No sequences passed the QC filters during fetching.", file=sys.stderr)
        sys.exit(1)

    SeqIO.write(valid_records, "results/raw/isolates.fasta", "fasta")

    db_accs = [r.id for r in valid_records if r.id not in local_meta_by_acc]
    meta_by_acc = dict(local_meta_by_acc)

    if db_accs:
        print(f"📊 Fetching metadata and spatiotemporal fields for {len(db_accs)} NCBI sequence(s)...")
        try:
            with Entrez.esummary(db="nucleotide", id=",".join(db_accs), retmode="json") as summary_handle:
                raw_json = json.load(summary_handle)

            result_dict = raw_json.get("result", {})
            uids = result_dict.get("uids", [])

            for uid in uids:
                s = result_dict.get(str(uid), {})
                if not isinstance(s, dict):
                    continue

                raw_acc = s.get("caption") or s.get("accessionversion") or ""
                clean_acc = clean_accession(raw_acc)
                title = s.get("title", "")
                update_date = s.get("updatedate", "")
                create_date = s.get("createdate", "") or update_date

                country, collection_date = parse_spatiotemporal_from_text(
                    title=title,
                    sub_type_str=s.get("subtype", ""),
                    sub_name_str=s.get("subname", ""),
                    fallback_date=create_date,
                )

                meta_by_acc[clean_acc] = {
                    "Accession": clean_acc,
                    "Title": title,
                    "Length_bp": int(s.get("slen") or 0),
                    "UpdateDate": update_date,
                    "Collection_Date": collection_date,
                    "Country": country,
                }
        except Exception as e:
            print(f"⚠️ Warning: esummary JSON metadata extraction triggered ({e})")

    meta_data = []
    for r in valid_records:
        if r.id in meta_by_acc:
            row = meta_by_acc[r.id]
            if not row.get("Length_bp"):
                row["Length_bp"] = len(r.seq)
            meta_data.append(row)
        else:
            orig_title = r.annotations.get("orig_title", r.id)
            country, collection_date = parse_spatiotemporal_from_text(title=orig_title)
            meta_data.append({
                "Accession": r.id,
                "Title": orig_title,
                "Length_bp": len(r.seq),
                "UpdateDate": "",
                "Collection_Date": collection_date,
                "Country": country,
            })

    pd.DataFrame(meta_data).to_csv("results/raw/metadata.csv", index=False)
    print(f"✅ Success! Saved {len(valid_records)} pristine sequence(s) and metadata.")


if __name__ == "__main__":
    main()